#!/bin/bash
#SBATCH -A cis261613
#SBATCH -p shared
#SBATCH --job-name=STAR-threads
#SBATCH --output=STAR-threads-%j.out
#SBATCH --error=STAR-threads-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=01:30:00

set -euo pipefail

PROOT=/home/x-sneisewander/hpds-rnaseq/STAR
READS=${SCRATCH}/genomes/human/subsample      # trimmed FASTQs (see below)
OUT=${SCRATCH}/star-threads/t${THREADS}_${SLURM_JOB_ID}/
RESULTS=${PROOT}/test/cpu/thread-results.csv

mkdir -p "${OUT}"
echo "CPUs allocated: ${THREADS}"

SECONDS=0
${PROOT}/bin/STAR \
    --runThreadN ${THREADS} \
    --genomeDir ${SCRATCH}/GRCh38/index \
    --readFilesIn ${READS}_1.fastq ${READS}_2.fastq \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${OUT}
WALL=${SECONDS}

# STAR's own numbers exclude index loading, which doesn't scale with threads
SPEED=$(grep "Mapping speed" ${OUT}Log.final.out | cut -f2)
echo "${THREADS},${WALL},${SPEED},${SLURM_JOB_ID}" >> "${RESULTS}"

rm -f ${OUT}Aligned.out.bam   # keep Log.final.out, drop the big BAM
