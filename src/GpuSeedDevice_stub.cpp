// Stub of the CUDA seed-search module, linked when STAR is built without CUDA (make CUDA=0, the default).
#include "GpuSeedDevice.h"

static const char *noCudaMsg="this STAR binary was compiled without CUDA support. SOLUTION: re-compile with \"make CUDA=1\", or use --gpuSeedSearch No";

bool gpuSeedCompiled() {return false;};

bool gpuSeedGenomeUpload(const GpuSeedGenomeHost &, std::string &, std::string &err)
{
    err=noCudaMsg;
    return false;
};

void gpuSeedGenomeFree() {};

GpuSeedDevBatch* gpuSeedDevBatchCreate() {return NULL;};

void gpuSeedDevBatchDestroy(GpuSeedDevBatch *) {};

bool gpuSeedDevRun(GpuSeedDevBatch *, const char *, uint64_t, const GpuSeedChain *, uint64_t,
                   uint32_t, const GpuSeedChainOut **, const GpuSeedRec **, uint64_t *, std::string &err)
{
    err=noCudaMsg;
    return false;
};

void gpuSeedRangePush(const char *) {};
void gpuSeedRangePop() {};
