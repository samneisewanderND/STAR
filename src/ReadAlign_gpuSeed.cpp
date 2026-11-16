#include "ReadAlign.h"
#include "readLoad.h"
#include "SequenceFuns.h"
#include "ErrorWarning.h"
#include <chrono>

void ReadAlign::gpuSeedPrepareBatch()
{//pre-load the next batch of reads from the input streams, run their seed search on the GPU, rewind the streams
    GpuSeedBatch &B=*gpuSeed;

    gpuSeedRangePush("gpuSeed prepare batch");
    auto t0=std::chrono::steady_clock::now();

    B.nReads=0;
    B.next=0;
    B.nRec=0;
    B.reads.clear();
    B.seq.clear();
    B.chains.clear();

    for (uint32 im=0; im<P.readNends; im++) {
        if (!readInStream[im]->good()) {//nothing to pre-load; keep the stream state untouched for oneRead()
            gpuSeedRangePop();
            return;
        };
    };

    streampos pos0[MAX_N_MATES];
    for (uint32 im=0; im<P.readNends; im++)
        pos0[im]=readInStream[im]->tellg();

    uint64 iReadAll1=0;
    uint32 readFilesIndex1=0;
    char readFilter1='N';

    while (B.nReads < B.maxReads) {
        //load the read exactly as ReadAlign::oneRead() does
        int readStatus[MAX_N_MATES]={-1,-1,-1};
        bool statusConsistent=true;
        for (uint32 im=0; im<P.readNends; im++) {
            readStatus[im] = readLoad(*(readInStream[im]), P, B.readLength[im], B.readLengthOriginal[im], B.readName[im], B.seq0[im], B.seqNum[im],
                                      B.qual[im], clipMates[im], iReadAll1, readFilesIndex1, readFilter1, B.readNameExtra[im]);
            if (readStatus[im] != readStatus[0]) {
                statusConsistent=false; //oneRead() will report the error when it reaches this read
                break;
            };
        };
        if (!statusConsistent || readStatus[0]==-1)
            break;

        GpuSeedReadInfo ri;
        ri.seqOff=B.seq.size();
        ri.chainStart=B.chains.size();
        ri.valid=true;

        char *R=B.seqNum[0];
        uint Lread1;
        if (P.readNmates==2) {//combine two mates together
            Lread1=B.readLength[0]+B.readLength[1]+1;
            if (Lread1>DEF_readSeqLengthMax) {//oneRead() will report the error when it reaches this read
                ri.valid=false;
            } else {
                R[B.readLength[0]]=MARK_FRAG_SPACER_BASE;
                complementSeqNumbers(B.seqNum[1], R+B.readLength[0]+1, B.readLength[1]);
                for (uint ii=0;ii<B.readLength[1]/2;ii++) {
                    swap(R[Lread1-ii-1],R[ii+B.readLength[0]+1]);
                };
            };
        } else {
            Lread1=B.readLength[0];
        };
        ri.Lread=Lread1;

        if (ri.valid) {
            B.seq.insert(B.seq.end(), R, R+Lread1);

            //enumerate seed-search chains in the order of ReadAlign::mapOneRead()
            uint Nsplit1 = Lread1>0 ? qualitySplit(R, Lread1, P.maxNsplit, P.seedSplitMin, B.splitR) : 0;
            if (Nsplit1>0) {
                uint seedSearchStartLmax=min(P.seedSearchStartLmax, (uint) (P.seedSearchStartLmaxOverLread*(Lread1-1)));
                for (uint ip=0; ip<Nsplit1; ip++) {
                    uint Nstart = P.seedSearchStartLmax>0 && seedSearchStartLmax<B.splitR[1][ip] ? B.splitR[1][ip]/seedSearchStartLmax+1 : 1;
                    uint Lstart = B.splitR[1][ip]/Nstart;
                    uint64 fwd0=B.chains.size(); //forward, istart=0 chain of this piece
                    for (uint iDir=0; iDir<2; iDir++) {
                        for (uint istart=0; istart<Nstart; istart++) {
                            GpuSeedChain c;
                            c.readOff=(uint32_t) ri.seqOff;
                            c.pieceStart=(uint32_t) B.splitR[0][ip];
                            c.pieceLen=(uint32_t) B.splitR[1][ip];
                            c.iFrag=(uint32_t) B.splitR[2][ip];
                            c.Lstart=(uint32_t) Lstart;
                            c.iDir=(uint32_t) iDir;
                            c.istart=(uint32_t) istart;
                            c.partner = (iDir==1 && istart==0) ? (uint32_t) fwd0 : GPU_SEED_NO_PARTNER;
                            B.chains.push_back(c);
                        };
                    };
                };
            };
        };

        ri.chainN=B.chains.size()-ri.chainStart;
        B.reads.push_back(ri);
        B.nReads++;
    };

    for (uint32 im=0; im<P.readNends; im++) {//rewind the streams to the start of the batch
        readInStream[im]->clear();
        readInStream[im]->seekg(pos0[im]);
    };

    auto t1=std::chrono::steady_clock::now();
    B.stats.timePrepare += std::chrono::duration<double>(t1-t0).count();
    gpuSeedRangePop();

    if (B.nReads==0)
        return;

    if (B.seq.size()>=(1LLU<<32)) {//offsets are 32-bit on the GPU
        for (auto &ri : B.reads) ri.valid=false;
        return;
    };

    if (B.dev==NULL)
        B.dev=gpuSeedDevBatchCreate();

    gpuSeedRangePush("gpuSeed GPU search");
    string err;
    uint64_t nRec1=0;
    bool runOK=gpuSeedDevRun(B.dev, B.seq.data(), B.seq.size(), B.chains.data(), B.chains.size(), (uint32_t) P.seedMapMin, &B.out, &B.rec, &nRec1, err);
    B.nRec=nRec1;
    if (!runOK) {
        exitWithError("EXITING because of FATAL ERROR in the GPU seed search: "+err+"\nSOLUTION: re-run with --gpuSeedSearch No\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_RUNTIME, P);
    };
    gpuSeedRangePop();

    B.stats.timeGpu += std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count();
    B.stats.batches++;
    B.stats.reads += B.nReads;
    B.stats.chains += B.chains.size();
    B.stats.records += B.nRec;
};

bool ReadAlign::gpuSeedReplay()
{//replay the GPU seed search results for the current read through storeAligns(), in the same order as the CPU loop in mapOneRead().
 //Returns false if the CPU seed search has to be used for this read; in this case nothing has been stored.
    if (gpuSeed==NULL || gpuSeedReadInd<0)
        return false;

    GpuSeedBatch &B=*gpuSeed;
    uint64 iR=(uint64) gpuSeedReadInd;
    gpuSeedReadInd=-1;
    if (iR>=B.nReads)
        return false;

    const GpuSeedReadInfo &ri=B.reads[iR];
    if (!ri.valid || ri.Lread!=Lread || (Lread>0 && memcmp(B.seq.data()+ri.seqOff, Read1[0], Lread)!=0)) {
        B.stats.readsCpu++;
        return false;
    };

    uint seedSearchStartLmax = Nsplit>0 ? min(P.seedSearchStartLmax, (uint) (P.seedSearchStartLmaxOverLread*(Lread-1))) : 0;

    //1st pass: check that the chains match this read's pieces and that all chains that will be replayed were completed on the GPU
    uint64 ic=ri.chainStart, icEnd=ri.chainStart+ri.chainN;
    for (uint ip=0; ip<Nsplit; ip++) {
        uint Nstart = P.seedSearchStartLmax>0 && seedSearchStartLmax<splitR[1][ip] ? splitR[1][ip]/seedSearchStartLmax+1 : 1;
        uint Lstart = splitR[1][ip]/Nstart;
        bool flagDirMap=true;
        for (uint iDir=0; iDir<2; iDir++) {
            for (uint istart=0; istart<Nstart; istart++, ic++) {
                if (ic>=icEnd) {
                    B.stats.readsCpu++;
                    return false;
                };
                const GpuSeedChain &c=B.chains[ic];
                const GpuSeedChainOut &o=B.out[ic];
                if (c.pieceStart!=splitR[0][ip] || c.pieceLen!=splitR[1][ip] || c.iFrag!=splitR[2][ip] || c.Lstart!=Lstart || c.iDir!=iDir || c.istart!=istart) {
                    B.stats.readsCpu++;
                    return false;
                };
                if (flagDirMap || istart>0) {
                    if (o.status!=GPU_SEED_ST_OK || o.recOff+o.recN>B.nRec) {
                        B.stats.readsCpu++;
                        return false;
                    };
                    if (iDir==0 && istart==0 && o.nSearch>0 && splitR[0][ip]+o.firstLbest==splitR[1][ip])
                        flagDirMap=false;
                };
            };
        };
    };
    if (ic!=icEnd) {
        B.stats.readsCpu++;
        return false;
    };

    //2nd pass: store the seeds
    ic=ri.chainStart;
    for (uint ip=0; ip<Nsplit; ip++) {
        uint Nstart = P.seedSearchStartLmax>0 && seedSearchStartLmax<splitR[1][ip] ? splitR[1][ip]/seedSearchStartLmax+1 : 1;
        bool flagDirMap=true;
        for (uint iDir=0; iDir<2; iDir++) {
            for (uint istart=0; istart<Nstart; istart++, ic++) {
                if (!flagDirMap && istart==0)
                    continue;
                const GpuSeedChainOut &o=B.out[ic];
                for (uint64 ir=o.recOff; ir<o.recOff+o.recN; ir++) {
                    const GpuSeedRec &r=B.rec[ir];
                    uint indStartEnd[2]={r.saStart, r.saEnd};
                    storeAligns(iDir, r.shift, indStartEnd[1]-indStartEnd[0]+1, r.L, indStartEnd, splitR[2][ip]);
                };
                if (iDir==0 && istart==0 && o.nSearch>0 && splitR[0][ip]+o.firstLbest==splitR[1][ip])
                    flagDirMap=false;
            };
        };
    };

    B.stats.readsGpu++;
    return true;
};
