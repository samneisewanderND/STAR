#!/bin/bash
# Map the same reads with CPU seeding (--gpuSeedSearch No) and GPU seeding (--gpuSeedSearch Yes)
# and check that the alignments are identical.
#
# usage: compare-gpu-seed.sh STAR_BINARY GENOME_DIR OUT_DIR READS_1 [READS_2] [-- extra STAR options]
# e.g.:  compare-gpu-seed.sh bin/STAR $SCRATCH/GRCh38/index $SCRATCH/cmp r_1.fq r_2.fq -- --runThreadN 16

set -euo pipefail

STAR=$1; GENOME=$2; OUT=$3; shift 3
READS=()
while [[ $# -gt 0 && "$1" != "--" ]]; do READS+=("$1"); shift; done
[[ $# -gt 0 ]] && shift
EXTRA=("$@")

mkdir -p "${OUT}"

for mode in No Yes; do
    echo "=== --gpuSeedSearch ${mode}"
    SECONDS=0
    "${STAR}" --genomeDir "${GENOME}" --readFilesIn "${READS[@]}" \
        --outSAMtype SAM --outSAMattributes NH HI AS nM NM MD jM jI \
        --outFileNamePrefix "${OUT}/${mode}_" --gpuSeedSearch ${mode} "${EXTRA[@]}"
    echo "wall time: ${SECONDS} s"
    grep -E "Mapping speed|Uniquely mapped reads %|% of reads mapped to multiple loci" "${OUT}/${mode}_Log.final.out"
    grep "GPU seed search summary" "${OUT}/${mode}_Log.out" || true
done

# the read order in the SAM depends on thread scheduling: compare sorted alignment lines (header excluded)
for mode in No Yes; do
    grep -v '^@' "${OUT}/${mode}_Aligned.out.sam" | LC_ALL=C sort -S 2G > "${OUT}/${mode}_sorted.sam"
done

status=0
if cmp -s "${OUT}/No_sorted.sam" "${OUT}/Yes_sorted.sam"; then
    echo "PASS: alignments are identical ($(wc -l < "${OUT}/No_sorted.sam") SAM lines)"
else
    echo "FAIL: alignments differ"
    diff "${OUT}/No_sorted.sam" "${OUT}/Yes_sorted.sam" | head -20
    status=1
fi
if cmp -s "${OUT}/No_SJ.out.tab" "${OUT}/Yes_SJ.out.tab"; then
    echo "PASS: SJ.out.tab identical"
else
    echo "FAIL: SJ.out.tab differs"
    status=1
fi
exit ${status}
