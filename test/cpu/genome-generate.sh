#!/bin/bash
#SBATCH -A cis261613
#SBATCH -p shared
#SBATCH --job-name=STAR-sequence
#SBATCH --output=STAR-generate-genome.out
#SBATCH --error=STAR-generate-genome.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=03:00:00

PROOT=/home/x-sneisewander/hpds-rnaseq/STAR
GENOME=human/GCF_000001405.40_GRCh38.p14_genomic

echo "CPUs allocated: " ${SLURM_CPUS_PER_TASK}
mkdir -p ${SCRATCH}/GRCh38/

${PROOT}/bin/STAR \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM 60000000000 \
	--runThreadN 16 \
	--genomeDir ${SCRATCH}/GRCh38/index \
	--genomeFastaFiles ${PROOT}/test/genomes/${GENOME}.fna \
	--sjdbGTFfile ${PROOT}/test/genomes/${GENOME}.gtf \
	--sjdbOverhang 149 \
	--genomeSAindexNbases 14 \
	--outFileNamePrefix ${SCRATCH}/GRCh38/
