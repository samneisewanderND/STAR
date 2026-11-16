// CUDA implementation of STAR's seed search (maximum mappable prefix search in the suffix array).
//
// Each GPU thread runs one "chain": the sequential series of maxMappableLength2strands() calls that
// ReadAlign::mapOneRead() performs for one (piece, direction, start) combination. The device functions
// below are line-by-line ports of the host functions in ReadAlign_maxMappableLength2strands.cpp and
// SuffixArrayFuns.cpp, so the results (and the order of the storeAligns calls replayed on the host)
// are identical to the CPU path.
//
// Testing without a GPU: compiling with -DGPU_SEED_HOST_EMULATION (make CUDA=1 NVCCFLAGS="-O3 -DGPU_SEED_HOST_EMULATION")
// runs the same chain code (seedOneChain) on the host instead of launching kernels. For validation only.

#include "GpuSeedDevice.h"

#include <cuda_runtime.h>
#include <nvtx3/nvToolsExt.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <sstream>

typedef unsigned long long u64;
typedef long long i64;

#define HD __host__ __device__
#ifdef __CUDA_ARCH__
    #define LDG(p) __ldg(p)
#else
    #define LDG(p) (*(p))
#endif

#define GPU_SEED_BLOCK 128
#define GPU_SEED_POOL_PER_CHAIN_INIT 4

/////////////////////////////////////////////////////////////////////////////// device genome
struct DevGenome {
    const signed char *G;
    const u64 *SA, *SAi;
    u64 saMask, saiMask;
    u64 nSA, nGenome;
    u64 GstrandMask, SAiMarkNmaskC, SAiMarkNmask, SAiMarkAbsentMaskC;
    uint32_t saWordLength, saiWordLength, GstrandBit, gSAindexNbases, gSAsparseD;
    u64 SAindexStart[GPU_SEED_MAX_SAINDEX_NBASES+1];
};

static DevGenome g_dev;
static bool g_devLoaded=false;
#ifndef GPU_SEED_HOST_EMULATION
static void *g_devGbase=NULL, *g_devSA=NULL, *g_devSAi=NULL;
#endif

/////////////////////////////////////////////////////////////////////////////// helpers
#define CUDA_TRY(call) do { cudaError_t e_=(call); if (e_!=cudaSuccess) { \
    std::ostringstream o_; o_ << #call << " failed: " << cudaGetErrorString(e_) << " (" << __FILE__ << ":" << __LINE__ << ")"; \
    err=o_.str(); return false; } } while(0)

HD __forceinline__ u64 packedGet(const u64 *__restrict__ w, uint32_t wl, u64 mask, u64 i)
{//same bits as PackedArray::operator[]: [i*wl, i*wl+wl)
    u64 b=i*wl;
    u64 q=b>>6;
    uint32_t sh=(uint32_t)(b&63);
    u64 v=LDG(w+q)>>sh;
    if (sh+wl>64)
        v |= LDG(w+q+1)<<(64-sh);
    return v & mask;
};

HD __forceinline__ u64 saGet(const DevGenome &g, u64 i)  {return packedGet(g.SA,  g.saWordLength,  g.saMask,  i);};
HD __forceinline__ u64 saiGet(const DevGenome &g, u64 i) {return packedGet(g.SAi, g.saiWordLength, g.saiMask, i);};

HD __forceinline__ signed char complBase(signed char x)
{//complementSeqNumbers
    return (x>=0 && x<=3) ? (signed char)(3-x) : x;
};

HD __forceinline__ u64 medianU64(u64 a, u64 b)
{
    return a/2 + b/2 + (a%2 + b%2)/2;
};

/////////////////////////////////////////////////////////////////////////////// ports of SuffixArrayFuns.cpp
//r is the numeric read sequence (Read1[0]); the complementary sequence (Read1[1]) is computed on the fly
HD u64 dCompareSeqToGenome(const DevGenome &g, const signed char *r, u64 S, u64 N, u64 L, u64 iSA, bool dirR, bool &compRes)
{
    u64 SAstr=saGet(g, iSA);
    bool dirG = (SAstr>>g.GstrandBit) == 0;
    SAstr &= g.GstrandMask;
    const signed char *G=g.G;

    if (dirR && dirG) {//forward on read, forward on genome
        const signed char *s = r + S + L;
        const signed char *gp = G + (i64)(SAstr + L);
        for (i64 ii=0; (u64) ii < N-L; ii++) {
            signed char sc=s[ii], gc=gp[ii];
            if (sc!=gc) {
                compRes = sc>gc;
                return ii+L;
            };
        };
        return N;
    } else if (dirR && !dirG) {
        const signed char *s = r + S + L;
        const signed char *gp = G + (i64)(g.nGenome-1-SAstr-L);
        for (i64 ii=0; (u64) ii < N-L; ii++) {
            signed char sc=complBase(s[ii]), gc=gp[-ii];
            if (sc!=gc) {
                compRes = !(sc>gc || gc>3);
                return ii+L;
            };
        };
        return N;
    } else if (!dirR && dirG) {
        const signed char *s = r + (i64)(S - L);
        const signed char *gp = G + (i64)(SAstr + L);
        for (i64 ii=0; (u64) ii < N-L; ii++) {
            signed char sc=complBase(s[-ii]), gc=gp[ii];
            if (sc!=gc) {
                compRes = sc>gc;
                return ii+L;
            };
        };
        return N;
    } else {//!dirR && !dirG
        const signed char *s = r + (i64)(S - L);
        const signed char *gp = G + (i64)(g.nGenome-1-SAstr-L);
        for (i64 ii=0; (u64) ii < N-L; ii++) {
            signed char sc=s[-ii], gc=gp[-ii];
            if (sc!=gc) {
                compRes = !(sc>gc || gc>3);
                return ii+L;
            };
        };
        return N;
    };
};

HD u64 dFindMultRange(const DevGenome &g, u64 i3, u64 L3, u64 i1, u64 L1, u64 i1a, u64 L1a, u64 i1b, u64 L1b,
                              const signed char *r, bool dirR, u64 S)
{
    bool compRes;
    if (L1<L3) {
        L1b=L1; i1b=i1; i1a=i3;
    } else {
        if (L1a<L1) {
            L1b=L1a; i1b=i1a; i1a=i1;
        };
    };
    while ( (i1b+1<i1a) | (i1b>i1a+1) ) {
        u64 i1c=medianU64(i1a,i1b);
        u64 L1c=dCompareSeqToGenome(g,r,S,L3,L1b,i1c,dirR,compRes);
        if (L1c==L3) {
            i1a=i1c;
        } else {
            i1b=i1c; L1b=L1c;
        };
    };
    return i1a;
};

HD u64 dMaxMappableLength(const DevGenome &g, const signed char *r, u64 S, u64 N, u64 i1, u64 i2, bool dirR, u64 &L, u64 *indStartEnd)
{
    bool compRes=false;
    u64 L1,L2,i3,L3,L1a,L1b,L2a,L2b,i1a,i1b,i2a,i2b;

    L1=dCompareSeqToGenome(g,r,S,N,L,i1,dirR,compRes);
    L2=dCompareSeqToGenome(g,r,S,N,L,i2,dirR,compRes);
    L= L1<L2 ? L1 : L2;

    L1a=L1;L1b=L1;i1a=i1;i1b=i1;
    L2a=L2;L2b=L2;i2a=i2;i2b=i2;

    i3=i1;L3=L1;
    while (i1+1<i2) {//main binary search loop
        i3=medianU64(i1,i2);
        L3=dCompareSeqToGenome(g,r,S,N,L,i3,dirR,compRes);

        if (L3==N) break;

        if (compRes) {
            if (L3>L1) {
               L1b=L1a; L1a=L1; i1b=i1a; i1a=i1;
            };
            i1=i3;L1=L3;
        } else {
            if (L3>L2) {
               L2b=L2a; L2a=L2; i2b=i2a; i2a=i2;
            };
            i2=i3;L2=L3;
        };
        L= L1<L2 ? L1 : L2;
    };

    if (L3<N) {
        if (L1>L2) {
            i3=i1;L3=L1;
        } else {
            i3=i2;L3=L2;
        };
    };

    i1=dFindMultRange(g,i3,L3,i1,L1,i1a,L1a,i1b,L1b,r,dirR,S);
    i2=dFindMultRange(g,i3,L3,i2,L2,i2a,L2a,i2b,L2b,r,dirR,S);

    L=L3;
    indStartEnd[0]=i1;
    indStartEnd[1]=i2;
    return i2-i1+1;
};

/////////////////////////////////////////////////////////////////////////////// port of ReadAlign::maxMappableLength2strands
//appends the storeAligns records to rec[n...]; returns GPU_SEED_ST_OK or an error status
HD uint32_t dMaxMappableLength2strands(const DevGenome &g, const signed char *r, u64 pieceStartIn, u64 pieceLengthIn, uint32_t iDir,
                                               GpuSeedRec *rec, uint32_t &n, u64 &maxLbest)
{
    u64 indStartEndAll[GPU_SEED_MAX_SPARSED][2], maxLall[GPU_SEED_MAX_SPARSED];
    maxLbest=0;
    bool dirR = iDir==0;

    u64 nDist = pieceLengthIn < (u64)g.gSAsparseD ? pieceLengthIn : (u64)g.gSAsparseD;
    for (u64 iDist=0; iDist<nDist; iDist++) {
        u64 pieceStart;
        u64 pieceLength=pieceLengthIn-iDist;

        u64 Lmax = (u64)g.gSAindexNbases < pieceLength ? (u64)g.gSAindexNbases : pieceLength;
        u64 ind1=0;
        if (dirR) {
            pieceStart=pieceStartIn+iDist;
            for (u64 ii=0;ii<Lmax;ii++) {
                ind1 <<=2LLU;
                ind1 += (u64) r[pieceStart+ii];
            };
        } else {
            pieceStart=pieceStartIn-iDist;
            for (u64 ii=0;ii<Lmax;ii++) {
                ind1 <<=2LLU;
                ind1 += ( 3-((u64) r[pieceStart-ii]) );
            };
        };

        //find SA boundaries
        u64 Lind=Lmax;
        u64 iSA1=0, iSA2;
        while (Lind>0) {
            iSA1=saiGet(g, g.SAindexStart[Lind-1]+ind1);
            if ((iSA1 & g.SAiMarkAbsentMaskC) == 0) {
                break;
            } else {
                --Lind;
                ind1 = ind1 >> 2;
            };
        };
        if (Lind==0)
            return GPU_SEED_ST_ERROR; //undefined on the CPU as well; let the CPU handle this read

        bool iSA2good = true;
        if (g.SAindexStart[Lind-1]+ind1+1 < g.SAindexStart[Lind]) {
            iSA2 = saiGet(g, g.SAindexStart[Lind-1]+ind1+1);
            if ( (iSA2 & g.SAiMarkAbsentMaskC) == 0) {
                iSA2 = (iSA2 & g.SAiMarkNmask) - 1;
            } else {
                iSA2 = g.nSA-1;
                iSA2good = false;
            };
        } else {
            iSA2=g.nSA-1;
            iSA2good = false;
        };

        u64 maxL, indStartEnd[2];
        bool iSA1noN = (iSA1 & g.SAiMarkNmaskC)==0;
        if (Lind < g.gSAindexNbases && iSA1noN && iSA2good) {//no need for SA search
            indStartEnd[0]=iSA1;
            indStartEnd[1]=iSA2;
            maxL=Lind;
        } else if (iSA1==iSA2 && iSA1noN && iSA2good) {//unique align already, just find maxL
            indStartEnd[0]=indStartEnd[1]=iSA1;
            bool comparRes;
            maxL=dCompareSeqToGenome(g, r, pieceStart, pieceLength, Lind, iSA1, dirR, comparRes);
        } else {//SA search
            if (iSA2good && iSA1noN) {
                maxL = Lind;
            } else {
                maxL=0;
            };
            dMaxMappableLength(g, r, pieceStart, pieceLength, iSA1 & g.SAiMarkNmask, iSA2, dirR, maxL, indStartEnd);
        };

        if (maxL+iDist > maxLbest)
            maxLbest=maxL+iDist;
        indStartEndAll[iDist][0]=indStartEnd[0];
        indStartEndAll[iDist][1]=indStartEnd[1];
        maxLall[iDist]=maxL;
    };

    for (u64 iDist=0; iDist<nDist; iDist++) {//same order as storeAligns calls on the CPU
        if ( (maxLall[iDist]+iDist) == maxLbest) {
            if (n>=GPU_SEED_CHAIN_REC_MAX)
                return GPU_SEED_ST_OVERFLOW;
            GpuSeedRec &a=rec[n++];
            a.saStart = indStartEndAll[iDist][0];
            a.saEnd   = indStartEndAll[iDist][1];
            a.shift   = (uint32_t) (dirR ? pieceStartIn+iDist : pieceStartIn-iDist);
            a.L       = (uint32_t) maxLall[iDist];
            a.Lbest   = (uint32_t) maxLbest;
            a.endOfSearch = 0;
        };
    };
    rec[n-1].endOfSearch=1;
    return GPU_SEED_ST_OK;
};

/////////////////////////////////////////////////////////////////////////////// kernel
//phase 1: all chains except reverse istart=0 chains; phase 2: reverse istart=0 chains, which are only
//searched if the forward istart=0 chain of the same piece did not map the full piece (flagDirMap in mapOneRead)
HD void seedOneChain(const DevGenome &g, const signed char *__restrict__ seq, const GpuSeedChain *__restrict__ chains,
                     GpuSeedChainOut *__restrict__ out, GpuSeedRec *__restrict__ pool, u64 *poolCounter, u64 poolCap,
                     u64 ic, uint32_t seedMapMin, int phase)
{
    const GpuSeedChain d=chains[ic];
    bool isRev0 = d.partner!=GPU_SEED_NO_PARTNER;
    if ( (phase==1) == isRev0 )
        return;

    GpuSeedChainOut o;
    o.recOff=0; o.recN=0; o.nSearch=0; o.firstLbest=0; o.status=GPU_SEED_ST_OK;

    if (phase==2) {
        const GpuSeedChainOut po=out[d.partner];
        if (po.status!=GPU_SEED_ST_OK || (po.nSearch>0 && (u64)d.pieceStart+po.firstLbest==(u64)d.pieceLen)) {
            o.status=GPU_SEED_ST_SKIPPED;
            out[ic]=o;
            return;
        };
    };

    const signed char *r = seq + d.readOff;
    GpuSeedRec rec[GPU_SEED_CHAIN_REC_MAX];
    uint32_t n=0;

    u64 pieceStart=d.pieceStart, pieceLen=d.pieceLen, Lstart=d.Lstart, istart=d.istart;
    u64 Lmapped=0;
    while ( istart*Lstart + Lmapped + seedMapMin < pieceLen ) {
        u64 Shift = d.iDir==0 ? ( pieceStart + istart*Lstart + Lmapped ) : ( pieceStart + pieceLen - istart*Lstart-1-Lmapped);
        u64 seedLength = pieceLen - Lmapped - istart*Lstart;
        u64 L;
        uint32_t st=dMaxMappableLength2strands(g, r, Shift, seedLength, d.iDir, rec, n, L);
        if (st!=GPU_SEED_ST_OK) {
            o.status=st;
            break;
        };
        if (o.nSearch==0)
            o.firstLbest=(uint32_t)L;
        o.nSearch++;
        Lmapped+=L;
    };

    if (o.status==GPU_SEED_ST_OK && n>0) {
        #ifdef __CUDA_ARCH__
            u64 off=atomicAdd(poolCounter, (u64)n);
        #else
            u64 off=*poolCounter;
            *poolCounter+=n;
        #endif
        if (off+n>poolCap) {
            o.status=GPU_SEED_ST_POOL_FULL;
        } else {
            for (uint32_t k=0; k<n; k++)
                pool[off+k]=rec[k];
            o.recOff=off;
            o.recN=n;
        };
    };
    out[ic]=o;
};

__global__ void __launch_bounds__(GPU_SEED_BLOCK)
kSeedChains(const DevGenome g, const signed char *__restrict__ seq, const GpuSeedChain *__restrict__ chains,
            GpuSeedChainOut *__restrict__ out, GpuSeedRec *__restrict__ pool, u64 *poolCounter, u64 poolCap,
            u64 nChains, uint32_t seedMapMin, int phase)
{
    u64 ic = (u64)blockIdx.x*blockDim.x + threadIdx.x;
    if (ic<nChains)
        seedOneChain(g, seq, chains, out, pool, poolCounter, poolCap, ic, seedMapMin, phase);
};

#ifdef GPU_SEED_HOST_EMULATION
/////////////////////////////////////////////////////////////////////////////// host emulation (validation only)
#include <vector>
static std::vector<u64> g_emuSA, g_emuSAi;

bool gpuSeedCompiled() {return true;};

void gpuSeedGenomeFree() {g_devLoaded=false;};

bool gpuSeedGenomeUpload(const GpuSeedGenomeHost &gh, std::string &info, std::string &)
{
    g_emuSA.assign((gh.SAlengthByte+7)/8+2, 0);
    memcpy(g_emuSA.data(), gh.SA, gh.SAlengthByte);
    g_emuSAi.assign((gh.SAilengthByte+7)/8+2, 0);
    memcpy(g_emuSAi.data(), gh.SAi, gh.SAilengthByte);

    g_dev.G = (const signed char*) gh.G;
    g_dev.SA = g_emuSA.data();
    g_dev.SAi = g_emuSAi.data();
    g_dev.saWordLength = gh.SAwordLength;
    g_dev.saiWordLength = gh.SAiwordLength;
    g_dev.saMask  = (~0LLU) >> (64-gh.SAwordLength);
    g_dev.saiMask = (~0LLU) >> (64-gh.SAiwordLength);
    g_dev.nSA = gh.nSA;
    g_dev.nGenome = gh.nGenome;
    g_dev.GstrandBit = gh.GstrandBit;
    g_dev.GstrandMask = gh.GstrandMask;
    g_dev.SAiMarkNmaskC = gh.SAiMarkNmaskC;
    g_dev.SAiMarkNmask = gh.SAiMarkNmask;
    g_dev.SAiMarkAbsentMaskC = gh.SAiMarkAbsentMaskC;
    g_dev.gSAindexNbases = gh.gSAindexNbases;
    g_dev.gSAsparseD = gh.gSAsparseD;
    for (uint32_t ii=0; ii<=GPU_SEED_MAX_SAINDEX_NBASES; ii++)
        g_dev.SAindexStart[ii] = ii<=gh.gSAindexNbases ? gh.genomeSAindexStart[ii] : 0;
    g_devLoaded=true;
    info="HOST EMULATION of the GPU seed search (GPU_SEED_HOST_EMULATION), validation only";
    return true;
};

struct GpuSeedDevBatch {
    std::vector<GpuSeedChainOut> out;
    std::vector<GpuSeedRec> rec;
    u64 poolPerChain;
};

GpuSeedDevBatch* gpuSeedDevBatchCreate()
{
    GpuSeedDevBatch *b=new GpuSeedDevBatch;
    b->poolPerChain=GPU_SEED_POOL_PER_CHAIN_INIT;
    return b;
};

void gpuSeedDevBatchDestroy(GpuSeedDevBatch *b) {delete b;};

bool gpuSeedDevRun(GpuSeedDevBatch *b, const char *seq, uint64_t, const GpuSeedChain *chains, uint64_t nChains,
                   uint32_t seedMapMin, const GpuSeedChainOut **out, const GpuSeedRec **rec, uint64_t *nRec, std::string &err)
{
    if (!g_devLoaded) {
        err="genome index was not uploaded";
        return false;
    };
    u64 poolCap=nChains*b->poolPerChain;
    b->out.resize(nChains);
    b->rec.resize(poolCap+1);
    u64 counter=0;
    for (int phase=1; phase<=2; phase++)
        for (u64 ic=0; ic<nChains; ic++)
            seedOneChain(g_dev, (const signed char*) seq, chains, b->out.data(), b->rec.data(), &counter, poolCap, ic, seedMapMin, phase);
    if (counter>poolCap && b->poolPerChain<GPU_SEED_CHAIN_REC_MAX)
        b->poolPerChain*=2;
    *out=b->out.data();
    *rec=b->rec.data();
    *nRec=std::min(counter, poolCap);
    return true;
};

void gpuSeedRangePush(const char *name) {nvtxRangePushA(name);};
void gpuSeedRangePop() {nvtxRangePop();};

#else //real GPU implementation

/////////////////////////////////////////////////////////////////////////////// host: genome
bool gpuSeedCompiled() {return true;};

static bool copyToDevice(void *dst, const void *src, u64 nBytes, std::string &err)
{
    const u64 block=1LLU<<30;
    for (u64 ii=0; ii<nBytes; ii+=block) {
        u64 n1=std::min(block, nBytes-ii);
        CUDA_TRY(cudaMemcpy((char*)dst+ii, (const char*)src+ii, n1, cudaMemcpyHostToDevice));
    };
    return true;
};

static bool uploadPacked(void *&dPtr, const char *hPtr, u64 lengthByte, std::string &err)
{//allocate with zeroed padding so that the 2nd 64-bit word of the last element can always be read
    u64 nAlloc=(lengthByte+7)/8*8+16;
    CUDA_TRY(cudaMalloc(&dPtr, nAlloc));
    CUDA_TRY(cudaMemset(dPtr, 0, nAlloc));
    return copyToDevice(dPtr, hPtr, lengthByte, err);
};

void gpuSeedGenomeFree()
{
    if (g_devGbase!=NULL) cudaFree(g_devGbase);
    if (g_devSA!=NULL)    cudaFree(g_devSA);
    if (g_devSAi!=NULL)   cudaFree(g_devSAi);
    g_devGbase=g_devSA=g_devSAi=NULL;
    g_devLoaded=false;
};

bool gpuSeedGenomeUpload(const GpuSeedGenomeHost &gh, std::string &info, std::string &err)
{
    gpuSeedGenomeFree();

    if (gh.gSAsparseD<1 || gh.gSAsparseD>GPU_SEED_MAX_SPARSED) {
        err="genomeSAsparseD="+std::to_string(gh.gSAsparseD)+" is not supported (max "+std::to_string(GPU_SEED_MAX_SPARSED)+")";
        return false;
    };
    if (gh.gSAindexNbases<1 || gh.gSAindexNbases>GPU_SEED_MAX_SAINDEX_NBASES) {
        err="genomeSAindexNbases="+std::to_string(gh.gSAindexNbases)+" is not supported (max "+std::to_string(GPU_SEED_MAX_SAINDEX_NBASES)+")";
        return false;
    };
    if (gh.SAwordLength>57 || gh.SAiwordLength>57) {
        err="packed suffix array word length is too large";
        return false;
    };

    auto t0=std::chrono::steady_clock::now();

    int dev=0;
    CUDA_TRY(cudaGetDevice(&dev));
    cudaDeviceProp prop;
    CUDA_TRY(cudaGetDeviceProperties(&prop, dev));
    size_t memFree=0, memTotal=0;
    CUDA_TRY(cudaMemGetInfo(&memFree, &memTotal));

    u64 nG=gh.padL+gh.nGenome+gh.padR;
    u64 need=nG + gh.SAlengthByte + gh.SAilengthByte;
    if (need + (512LLU<<20) > memFree) {//leave some room for per-thread batch buffers
        std::ostringstream o;
        o << "not enough GPU memory on device " << dev << " (" << prop.name << "): genome index needs " << need/1000000 << " MB, free "
          << memFree/1000000 << " MB. SOLUTION: use a GPU with more memory, or a sparser suffix array (genomeGenerate --genomeSAsparseD 2)";
        err=o.str();
        return false;
    };

    nvtxRangePushA("gpuSeed genome upload");
    CUDA_TRY(cudaMalloc(&g_devGbase, nG));
    if (!copyToDevice(g_devGbase, gh.G-gh.padL, nG, err)) return false;
    if (!uploadPacked(g_devSA,  gh.SA,  gh.SAlengthByte,  err)) return false;
    if (!uploadPacked(g_devSAi, gh.SAi, gh.SAilengthByte, err)) return false;
    nvtxRangePop();

    g_dev.G = (const signed char*) g_devGbase + gh.padL;
    g_dev.SA = (const u64*) g_devSA;
    g_dev.SAi = (const u64*) g_devSAi;
    g_dev.saWordLength = gh.SAwordLength;
    g_dev.saiWordLength = gh.SAiwordLength;
    g_dev.saMask  = (~0LLU) >> (64-gh.SAwordLength);
    g_dev.saiMask = (~0LLU) >> (64-gh.SAiwordLength);
    g_dev.nSA = gh.nSA;
    g_dev.nGenome = gh.nGenome;
    g_dev.GstrandBit = gh.GstrandBit;
    g_dev.GstrandMask = gh.GstrandMask;
    g_dev.SAiMarkNmaskC = gh.SAiMarkNmaskC;
    g_dev.SAiMarkNmask = gh.SAiMarkNmask;
    g_dev.SAiMarkAbsentMaskC = gh.SAiMarkAbsentMaskC;
    g_dev.gSAindexNbases = gh.gSAindexNbases;
    g_dev.gSAsparseD = gh.gSAsparseD;
    for (uint32_t ii=0; ii<=GPU_SEED_MAX_SAINDEX_NBASES; ii++)
        g_dev.SAindexStart[ii] = ii<=gh.gSAindexNbases ? gh.genomeSAindexStart[ii] : 0;

    g_devLoaded=true;

    double dt=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::ostringstream o;
    o << "GPU device " << dev << ": " << prop.name << ", " << memTotal/1000000 << " MB total; uploaded genome index: "
      << need/1000000 << " MB (G=" << nG/1000000 << " MB, SA=" << gh.SAlengthByte/1000000 << " MB, SAi=" << gh.SAilengthByte/1000000
      << " MB) in " << dt << " s";
    info=o.str();
    return true;
};

/////////////////////////////////////////////////////////////////////////////// host: per-thread batches
struct GpuSeedDevBatch {
    cudaStream_t stream;
    bool streamCreated;

    //pinned host buffers
    char *hSeq; u64 hSeqCap;
    GpuSeedChain *hChain; u64 hChainCap;
    GpuSeedChainOut *hOut; u64 hOutCap;
    GpuSeedRec *hRec; u64 hRecCap;
    u64 *hCounter;

    //device buffers
    char *dSeq; u64 dSeqCap;
    GpuSeedChain *dChain; u64 dChainCap;
    GpuSeedChainOut *dOut; u64 dOutCap;
    GpuSeedRec *dRec; u64 dRecCap;
    u64 *dCounter;

    u64 poolPerChain;
};

template <class T> static bool growHost(T *&p, u64 &cap, u64 need, std::string &err)
{
    if (need<=cap) return true;
    if (p!=NULL) cudaFreeHost(p);
    p=NULL; cap=0;
    u64 n=std::max(need, need/4*5);
    CUDA_TRY(cudaMallocHost((void**)&p, n*sizeof(T)));
    cap=n;
    return true;
};

template <class T> static bool growDev(T *&p, u64 &cap, u64 need, std::string &err)
{
    if (need<=cap) return true;
    if (p!=NULL) cudaFree(p);
    p=NULL; cap=0;
    u64 n=std::max(need, need/4*5);
    CUDA_TRY(cudaMalloc((void**)&p, n*sizeof(T)));
    cap=n;
    return true;
};

GpuSeedDevBatch* gpuSeedDevBatchCreate()
{
    GpuSeedDevBatch *b=new GpuSeedDevBatch;
    memset(b, 0, sizeof(*b));
    b->poolPerChain=GPU_SEED_POOL_PER_CHAIN_INIT;
    return b;
};

void gpuSeedDevBatchDestroy(GpuSeedDevBatch *b)
{
    if (b==NULL) return;
    if (b->streamCreated) cudaStreamDestroy(b->stream);
    if (b->hSeq) cudaFreeHost(b->hSeq);
    if (b->hChain) cudaFreeHost(b->hChain);
    if (b->hOut) cudaFreeHost(b->hOut);
    if (b->hRec) cudaFreeHost(b->hRec);
    if (b->hCounter) cudaFreeHost(b->hCounter);
    if (b->dSeq) cudaFree(b->dSeq);
    if (b->dChain) cudaFree(b->dChain);
    if (b->dOut) cudaFree(b->dOut);
    if (b->dRec) cudaFree(b->dRec);
    if (b->dCounter) cudaFree(b->dCounter);
    delete b;
};

bool gpuSeedDevRun(GpuSeedDevBatch *b, const char *seq, uint64_t nSeq, const GpuSeedChain *chains, uint64_t nChains,
                   uint32_t seedMapMin, const GpuSeedChainOut **out, const GpuSeedRec **rec, uint64_t *nRec, std::string &err)
{
    if (!g_devLoaded) {
        err="genome index was not uploaded to the GPU";
        return false;
    };

    if (!b->streamCreated) {
        CUDA_TRY(cudaStreamCreateWithFlags(&b->stream, cudaStreamNonBlocking));
        CUDA_TRY(cudaMallocHost((void**)&b->hCounter, sizeof(u64)));
        CUDA_TRY(cudaMalloc((void**)&b->dCounter, sizeof(u64)));
        b->streamCreated=true;
    };

    *nRec=0;
    if (nChains==0) {
        *out=b->hOut;
        *rec=b->hRec;
        return true;
    };

    u64 poolCap=nChains*b->poolPerChain;

    if (!growHost(b->hSeq, b->hSeqCap, nSeq+1, err)) return false;
    if (!growHost(b->hChain, b->hChainCap, nChains, err)) return false;
    if (!growHost(b->hOut, b->hOutCap, nChains, err)) return false;
    if (!growHost(b->hRec, b->hRecCap, poolCap, err)) return false;
    if (!growDev(b->dSeq, b->dSeqCap, nSeq+1, err)) return false;
    if (!growDev(b->dChain, b->dChainCap, nChains, err)) return false;
    if (!growDev(b->dOut, b->dOutCap, nChains, err)) return false;
    if (!growDev(b->dRec, b->dRecCap, poolCap, err)) return false;

    memcpy(b->hSeq, seq, nSeq);
    memcpy(b->hChain, chains, nChains*sizeof(GpuSeedChain));

    cudaStream_t st=b->stream;
    CUDA_TRY(cudaMemcpyAsync(b->dSeq, b->hSeq, nSeq, cudaMemcpyHostToDevice, st));
    CUDA_TRY(cudaMemcpyAsync(b->dChain, b->hChain, nChains*sizeof(GpuSeedChain), cudaMemcpyHostToDevice, st));
    CUDA_TRY(cudaMemsetAsync(b->dCounter, 0, sizeof(u64), st));

    u64 nBlocks=(nChains+GPU_SEED_BLOCK-1)/GPU_SEED_BLOCK;
    for (int phase=1; phase<=2; phase++) {
        kSeedChains<<<(unsigned int)nBlocks, GPU_SEED_BLOCK, 0, st>>>(g_dev, (const signed char*)b->dSeq, b->dChain, b->dOut,
                                                                     b->dRec, b->dCounter, poolCap, nChains, seedMapMin, phase);
        CUDA_TRY(cudaGetLastError());
    };

    CUDA_TRY(cudaMemcpyAsync(b->hOut, b->dOut, nChains*sizeof(GpuSeedChainOut), cudaMemcpyDeviceToHost, st));
    CUDA_TRY(cudaMemcpyAsync(b->hCounter, b->dCounter, sizeof(u64), cudaMemcpyDeviceToHost, st));
    CUDA_TRY(cudaStreamSynchronize(st));

    u64 nUsed=std::min(*b->hCounter, poolCap);
    if (nUsed>0) {
        CUDA_TRY(cudaMemcpyAsync(b->hRec, b->dRec, nUsed*sizeof(GpuSeedRec), cudaMemcpyDeviceToHost, st));
        CUDA_TRY(cudaStreamSynchronize(st));
    };

    if (*b->hCounter>poolCap && b->poolPerChain<GPU_SEED_CHAIN_REC_MAX)
        b->poolPerChain*=2; //pool was too small: chains that did not fit fall back to CPU; enlarge for the next batches

    *out=b->hOut;
    *rec=b->hRec;
    *nRec=nUsed;
    return true;
};

/////////////////////////////////////////////////////////////////////////////// profiling ranges
void gpuSeedRangePush(const char *name) {nvtxRangePushA(name);};
void gpuSeedRangePop() {nvtxRangePop();};

#endif //GPU_SEED_HOST_EMULATION
