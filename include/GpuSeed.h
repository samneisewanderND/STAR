#ifndef H_GpuSeed
#define H_GpuSeed

// Host-side part of the GPU seed search (--gpuSeedSearch Yes).
//
// Flow, per mapping thread (see ReadAlignChunk::mapChunk and ReadAlign_gpuSeed.cpp):
//   1. ReadAlign::gpuSeedPrepareBatch() pre-loads up to --gpuSeedBatchReads reads from the input stream,
//      reproduces the read preparation of ReadAlign::oneRead()/mapOneRead() (mates merging, qualitySplit),
//      enumerates the seed-search chains, runs them on the GPU, and rewinds the input stream.
//   2. The normal per-read processing continues; ReadAlign::mapOneRead() replays the GPU results through
//      storeAligns() in the original order (ReadAlign::gpuSeedReplay()), or falls back to the CPU search.

#include "IncludeDefine.h"
#include "Parameters.h"
#include "GpuSeedDevice.h"

class Genome;

struct GpuSeedReadInfo {
    uint64 seqOff;      //offset of Read1[0] in GpuSeedBatch::seq
    uint64 chainStart;  //first chain of this read in GpuSeedBatch::chains
    uint64 chainN;      //number of chains
    uint64 Lread;
    bool valid;         //false: this read is mapped with the CPU seed search
};

struct GpuSeedStats {
    uint64 batches=0, reads=0, readsGpu=0, readsCpu=0, chains=0, records=0;
    double timePrepare=0, timeGpu=0;
    void add(const GpuSeedStats &s);
};

class GpuSeedBatch {
public:
    GpuSeedBatch(Parameters &Pin);
    ~GpuSeedBatch();

    uint64 maxReads;    //--gpuSeedBatchReads
    uint64 nReads;      //reads in the current batch
    uint64 next;        //next batch read to be processed by ReadAlign::oneRead()

    vector<GpuSeedReadInfo> reads;
    vector<char> seq;               //Read1[0] of all reads
    vector<GpuSeedChain> chains;
    const GpuSeedChainOut *out;     //GPU results, one per chain
    const GpuSeedRec *rec;          //record pool
    uint64 nRec;

    //scratch arrays for pre-loading reads (same sizes as in ReadAlign)
    vector<char*> readName, seq0, seqNum, qual;
    vector<string> readNameExtra;
    uint readLength[MAX_N_MATES], readLengthOriginal[MAX_N_MATES];
    uint *splitR[3];

    GpuSeedStats stats;
    GpuSeedDevBatch *dev;

    void reset() {nReads=0; next=0;};
    void logStats(ostream &logOut, int iThread); //also adds stats to the global summary

private:
    Parameters &P;
};

//uploads the genome index to the GPU if it is not there yet or changed (e.g. after on-the-fly junction insertion).
//Must be called by the main thread before the mapping threads are spawned.
void gpuSeedGenomeSync(Genome &mapGen, Parameters &P);
//true if the GPU seed search is enabled and supported for the current genome/parameters
bool gpuSeedActive();
//writes the summary of all threads to Log.out and stdout
void gpuSeedLogSummary(Parameters &P);
//releases the GPU copy of the genome
void gpuSeedRelease();

#endif
