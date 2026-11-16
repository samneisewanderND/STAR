#include "GpuSeed.h"
#include "Genome.h"
#include "ErrorWarning.h"
#include "TimeFunctions.h"
#include <pthread.h>

static bool g_gpuSeedActive=false;

static struct {
    bool valid=false;
    const char *G=NULL, *SA=NULL, *SAi=NULL;
    uint64 nGenome=0, nSA=0, SAlengthByte=0;
} g_uploaded; //what is currently on the GPU

static pthread_mutex_t g_statsMutex=PTHREAD_MUTEX_INITIALIZER;
static GpuSeedStats g_statsTotal;

void GpuSeedStats::add(const GpuSeedStats &s)
{
    batches+=s.batches; reads+=s.reads; readsGpu+=s.readsGpu; readsCpu+=s.readsCpu;
    chains+=s.chains; records+=s.records;
    timePrepare+=s.timePrepare; timeGpu+=s.timeGpu;
};

GpuSeedBatch::GpuSeedBatch(Parameters &Pin) : maxReads(Pin.gpuSeed.batchReads), nReads(0), next(0),
                                              out(NULL), rec(NULL), nRec(0), dev(NULL), P(Pin)
{
    readName.resize(P.readNends);
    seq0.resize(P.readNends);
    seqNum.resize(P.readNends);
    qual.resize(P.readNends);
    readNameExtra.resize(P.readNends);
    for (uint32 im=0; im<P.readNends; im++) {
        readName[im] = new char[DEF_readNameLengthMax];
        seq0[im]     = new char[DEF_readSeqLengthMax+1];
        seqNum[im]   = new char[DEF_readSeqLengthMax+1];
        qual[im]     = new char[DEF_readSeqLengthMax+1];
    };
    for (int ii=0; ii<3; ii++)
        splitR[ii] = new uint[P.maxNsplit];
};

GpuSeedBatch::~GpuSeedBatch()
{
    gpuSeedDevBatchDestroy(dev);
    for (uint32 im=0; im<readName.size(); im++) {
        delete[] readName[im];
        delete[] seq0[im];
        delete[] seqNum[im];
        delete[] qual[im];
    };
    for (int ii=0; ii<3; ii++)
        delete[] splitR[ii];
};

static void statsWrite(ostream &o, const GpuSeedStats &s)
{
    o << "batches=" << s.batches << "; reads pre-loaded=" << s.reads << "; reads seeded on GPU=" << s.readsGpu
      << "; reads seeded on CPU (fallback)=" << s.readsCpu << "; chains=" << s.chains << "; seed records=" << s.records
      << "; pre-load time=" << s.timePrepare << " s; GPU wait time=" << s.timeGpu << " s";
};

void GpuSeedBatch::logStats(ostream &logOut, int iThread)
{
    logOut << "GPU seed search, thread #" << iThread << ": ";
    statsWrite(logOut, stats);
    logOut << '\n';

    pthread_mutex_lock(&g_statsMutex);
    g_statsTotal.add(stats);
    pthread_mutex_unlock(&g_statsMutex);
    stats=GpuSeedStats();
};

void gpuSeedLogSummary(Parameters &P)
{
    if (!P.gpuSeed.yes)
        return;
    pthread_mutex_lock(&g_statsMutex);
    ostringstream o;
    o << "GPU seed search summary (all threads, all mapping passes): ";
    statsWrite(o, g_statsTotal);
    pthread_mutex_unlock(&g_statsMutex);
    P.inOut->logMain << o.str() << endl;
    *P.inOut->logStdOut << timeMonthDayTime() << " ..... " << o.str() << endl;
};

bool gpuSeedActive()
{
    return g_gpuSeedActive;
};

void gpuSeedRelease()
{
    gpuSeedGenomeFree();
    g_uploaded.valid=false;
    g_gpuSeedActive=false;
};

void gpuSeedGenomeSync(Genome &mapGen, Parameters &P)
{
    if (!P.gpuSeed.yes) {
        g_gpuSeedActive=false;
        return;
    };

    string notSupported;
    if (P.pGe.gType==101) {
        notSupported="--genomeType SuperTranscriptome";
    } else if (P.seedSearchLmax>0) {
        notSupported="--seedSearchLmax > 0";
    } else if (mapGen.pGe.gSAsparseD > GPU_SEED_MAX_SPARSED) {
        notSupported="--genomeSAsparseD > "+to_string(GPU_SEED_MAX_SPARSED);
    } else if (mapGen.pGe.gSAindexNbases > GPU_SEED_MAX_SAINDEX_NBASES) {
        notSupported="--genomeSAindexNbases > "+to_string(GPU_SEED_MAX_SAINDEX_NBASES);
    };
    if (!notSupported.empty()) {
        warningMessage("--gpuSeedSearch Yes is not supported with "+notSupported+". The seed search will run on the CPU.",
                       std::cerr, P.inOut->logMain, P);
        g_gpuSeedActive=false;
        return;
    };

    if (g_uploaded.valid && g_uploaded.G==mapGen.G && g_uploaded.SA==mapGen.SA.charArray && g_uploaded.SAi==mapGen.SAi.charArray
        && g_uploaded.nGenome==mapGen.nGenome && g_uploaded.nSA==mapGen.nSA && g_uploaded.SAlengthByte==mapGen.SA.lengthByte) {
        g_gpuSeedActive=true; //same genome is already on the GPU
        return;
    };

    *P.inOut->logStdOut << timeMonthDayTime() << " ..... loading genome index into GPU memory\n" << flush;
    P.inOut->logMain    << timeMonthDayTime() << " ..... loading genome index into GPU memory\n" << flush;

    GpuSeedGenomeHost gh;
    //the genome sequence is padded by genomeLoad() with GENOME_spacingChar on both sides; the seed search may read into the padding
    const uint64 gPad=200;
    gh.G = mapGen.G;
    gh.nGenome = mapGen.nGenome;
    gh.padL = (mapGen.G1!=NULL && mapGen.G1<=mapGen.G && (uint64)(mapGen.G-mapGen.G1)<=gPad) ? (uint64)(mapGen.G-mapGen.G1) : 0;
    gh.padR = gPad;
    gh.SA = mapGen.SA.charArray;
    gh.SAlengthByte = mapGen.SA.lengthByte;
    gh.SAwordLength = (uint32) mapGen.SA.wordLength;
    gh.nSA = mapGen.nSA;
    gh.SAi = mapGen.SAi.charArray;
    gh.SAilengthByte = mapGen.SAi.lengthByte;
    gh.SAiwordLength = (uint32) mapGen.SAi.wordLength;
    gh.GstrandBit = mapGen.GstrandBit;
    gh.GstrandMask = mapGen.GstrandMask;
    gh.SAiMarkNmaskC = mapGen.SAiMarkNmaskC;
    gh.SAiMarkNmask = mapGen.SAiMarkNmask;
    gh.SAiMarkAbsentMaskC = mapGen.SAiMarkAbsentMaskC;
    gh.gSAindexNbases = (uint32) mapGen.pGe.gSAindexNbases;
    gh.gSAsparseD = (uint32) mapGen.pGe.gSAsparseD;
    gh.genomeSAindexStart = (const uint64_t*) mapGen.genomeSAindexStart;

    g_uploaded.valid=false;
    string info, err;
    if (!gpuSeedGenomeUpload(gh, info, err)) {
        exitWithError("EXITING because of FATAL ERROR: could not load the genome index into GPU memory for --gpuSeedSearch Yes: "+err+"\n",
                      std::cerr, P.inOut->logMain, EXIT_CODE_MEMORY_ALLOCATION, P);
    };

    g_uploaded.valid=true;
    g_uploaded.G=mapGen.G;
    g_uploaded.SA=mapGen.SA.charArray;
    g_uploaded.SAi=mapGen.SAi.charArray;
    g_uploaded.nGenome=mapGen.nGenome;
    g_uploaded.nSA=mapGen.nSA;
    g_uploaded.SAlengthByte=mapGen.SA.lengthByte;
    g_gpuSeedActive=true;

    P.inOut->logMain << timeMonthDayTime() << " ..... " << info << '\n' << flush;
    *P.inOut->logStdOut << timeMonthDayTime() << " ..... finished loading genome index into GPU memory\n" << flush;
};
