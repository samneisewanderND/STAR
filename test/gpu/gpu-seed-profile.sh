#!/bin/bash
#SBATCH -A cis261613-gpu
#SBATCH -p gpu
#SBATCH --job-name=STAR-gpu-seed
#SBATCH --output=STAR-gpu-seed-%j.out
#SBATCH --error=STAR-gpu-seed-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32   # GPU nodes are shared: 128 cores / 4 GPUs = 32 cores per GPU
#SBATCH --mem=80G
#SBATCH --time=00:30:00

# Aligns the human subsample to the GRCh38 index with CPU seeding and with GPU seeding,
# and appends the runtimes to a CSV file (one row per run).
# Build first with CUDA:
#   module load modtree/gpu gcc/11.2.0 cuda/12.8.0 && make -j16 STAR CUDA=1 CXX=g++
# Optional: MODES="Yes" sbatch gpu-seed-profile.sh

set -euo pipefail
module load modtree/gpu gcc/11.2.0 cuda/12.8.0

PROOT=/home/x-sneisewander/hpds-rnaseq/STAR
GENOME=/anvil/scratch/x-sneisewander/GRCh38/index
READS1=/anvil/scratch/x-sneisewander/genomes/human/subsample_1.fastq
READS2=/anvil/scratch/x-sneisewander/genomes/human/subsample_2.fastq
THREADS=32
MODES=${MODES:-"No Yes"}   # values of --gpuSeedSearch to run
OUT=/anvil/scratch/x-sneisewander/star-gpu-seed/t${THREADS}_${SLURM_JOB_ID}
RESULTS=${PROOT}/test/gpu/gpu-seed-results.csv

mkdir -p "${OUT}"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -n 1)

if [[ ! -s "${RESULTS}" ]]; then
    echo "job_id,gpu_seed_search,threads,gpu,wall_s,mapping_speed_Mreads_per_hour,reads_seeded_gpu,reads_seeded_cpu_fallback,gpu_wait_s" > "${RESULTS}"
fi

for MODE in ${MODES}; do
    echo "=== --gpuSeedSearch ${MODE}, ${THREADS} threads"
    PREFIX=${OUT}/gpu${MODE}_

    SECONDS=0
    ${PROOT}/bin/STAR \
        --runThreadN ${THREADS} \
        --genomeDir ${GENOME} \
        --readFilesIn ${READS1} ${READS2} \
        --outSAMtype BAM Unsorted \
        --gpuSeedSearch ${MODE} \
        --outFileNamePrefix ${PREFIX}
    WALL=${SECONDS}

    # STAR's mapping speed excludes the genome loading; for GPU runs it includes the index upload to the GPU
    SPEED=$(grep "Mapping speed" ${PREFIX}Log.final.out | cut -f2)

    # GPU seeding statistics (summed over all threads), empty for CPU runs
    SUMMARY=$(grep "GPU seed search summary" ${PREFIX}Log.out || true)
    GPU_READS=$(echo "${SUMMARY}" | sed -n 's/.*reads seeded on GPU=\([0-9]*\).*/\1/p')
    CPU_READS=$(echo "${SUMMARY}" | sed -n 's/.*reads seeded on CPU (fallback)=\([0-9]*\).*/\1/p')
    GPU_WAIT=$(echo "${SUMMARY}" | sed -n 's/.*GPU wait time=\([0-9.e+-]*\) s.*/\1/p')

    echo "${SLURM_JOB_ID},${MODE},${THREADS},${GPU_NAME},${WALL},${SPEED},${GPU_READS},${CPU_READS},${GPU_WAIT}" >> "${RESULTS}"
    echo "wall time: ${WALL} s, mapping speed: ${SPEED} M reads/hour"

    rm -f ${PREFIX}Aligned.out.bam   # keep the logs, drop the large BAM
done

echo "results appended to ${RESULTS}"
