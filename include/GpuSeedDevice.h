#ifndef H_GpuSeedDevice
#define H_GpuSeedDevice

// Interface between the host-side GPU seeding logic (C++) and the CUDA device module.
// This header is shared by nvcc and g++ translation units, so it must only use
// standard fixed-width types: do NOT include IncludeDefine.h here (it #defines uint etc.).

#include <stdint.h>
#include <string>

#define GPU_SEED_MAX_SPARSED        16 //max supported --genomeSAsparseD
#define GPU_SEED_MAX_SAINDEX_NBASES 24 //max supported --genomeSAindexNbases
#define GPU_SEED_CHAIN_REC_MAX      32 //max storeAligns records per chain; more -> read falls back to CPU
#define GPU_SEED_NO_PARTNER         0xFFFFFFFFu

//chain status values
#define GPU_SEED_ST_OK          1 //chain searched, records valid
#define GPU_SEED_ST_SKIPPED     2 //reverse-direction istart=0 chain skipped since the forward search mapped the full piece
#define GPU_SEED_ST_OVERFLOW    3 //too many records for one chain
#define GPU_SEED_ST_ERROR       4 //unexpected state (e.g. no SAi prefix of length 1)
#define GPU_SEED_ST_POOL_FULL   5 //record pool exhausted

//description of the genome index arrays in host memory, used to upload them to the GPU
struct GpuSeedGenomeHost {
    const char *G;              //genome sequence, G[-padL ... nGenome+padR-1] must be valid
    uint64_t nGenome, padL, padR;

    const char *SA;             //packed suffix array
    uint64_t SAlengthByte, nSA;
    uint32_t SAwordLength;

    const char *SAi;            //packed suffix array prefix index
    uint64_t SAilengthByte;
    uint32_t SAiwordLength;

    uint32_t GstrandBit;
    uint64_t GstrandMask;
    uint64_t SAiMarkNmaskC, SAiMarkNmask, SAiMarkAbsentMaskC;

    uint32_t gSAindexNbases, gSAsparseD;
    const uint64_t *genomeSAindexStart; //gSAindexNbases+1 values
};

//one seed-search chain: the sequential series of maximum mappable prefix searches
//performed in ReadAlign::mapOneRead for one (piece, direction, start) combination
struct GpuSeedChain {
    uint32_t readOff;     //offset of the read numeric sequence in the batch sequence buffer
    uint32_t pieceStart;  //splitR[0][ip]
    uint32_t pieceLen;    //splitR[1][ip]
    uint32_t iFrag;       //splitR[2][ip], only used on the host
    uint32_t Lstart;
    uint32_t iDir;
    uint32_t istart;
    uint32_t partner;     //for reverse istart=0 chains: index of the forward istart=0 chain of the same piece
};

struct GpuSeedChainOut {
    uint64_t recOff;      //offset of the first record in the record pool
    uint32_t recN;        //number of records
    uint32_t status;      //GPU_SEED_ST_*
    uint32_t nSearch;     //number of maxMappableLength2strands calls
    uint32_t firstLbest;  //maxLbest of the first search
};

//one ReadAlign::storeAligns call
struct GpuSeedRec {
    uint64_t saStart, saEnd; //indStartEnd; Nrep=saEnd-saStart+1
    uint32_t shift;          //read position passed to storeAligns
    uint32_t L;              //mapped length passed to storeAligns
    uint32_t Lbest;          //maxLbest of the search this record belongs to
    uint32_t endOfSearch;    //1 for the last record of a search
};

struct GpuSeedDevBatch; //opaque per-thread device context

//true if this STAR binary was compiled with CUDA support
bool gpuSeedCompiled();

//upload genome index to the GPU (frees the previous copy if any). Returns false and sets err on failure.
//info receives a human-readable description of the device/upload
bool gpuSeedGenomeUpload(const GpuSeedGenomeHost &gh, std::string &info, std::string &err);
void gpuSeedGenomeFree();

GpuSeedDevBatch* gpuSeedDevBatchCreate();
void gpuSeedDevBatchDestroy(GpuSeedDevBatch *b);

//run seed search for nChains chains on the GPU (blocking for the calling thread).
//On success, *out points to nChains chain results and *rec to *nRec records; both remain valid until the next call.
bool gpuSeedDevRun(GpuSeedDevBatch *b, const char *seq, uint64_t nSeq, const GpuSeedChain *chains, uint64_t nChains,
                   uint32_t seedMapMin, const GpuSeedChainOut **out, const GpuSeedRec **rec, uint64_t *nRec, std::string &err);

//NVTX ranges for profiling with Nsight Systems (no-ops without CUDA)
void gpuSeedRangePush(const char *name);
void gpuSeedRangePop();

#endif
