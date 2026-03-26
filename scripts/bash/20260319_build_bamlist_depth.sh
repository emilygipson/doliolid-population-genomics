#!/bin/bash
#SBATCH --job-name=bamlist_depth_120
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --dependency=afterok:43720322
#SBATCH --output=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.out
#SBATCH --error=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
# =============================================================================
# 20260319_build_bamlist_depth.sh
#
# PURPOSE:
#   1. Build the 120-sample BAM list (99 existing + 21 new)
#   2. Verify all BAMs exist and are indexed
#   3. Run ANGSD depth scan on all 120 samples
#   4. Compute per-sample read counts from flagstat for downsampling decision
#
# DEPENDS ON: 20260319_map_markdup_21.sh (job 43719593)
# =============================================================================
set -euo pipefail

module purge
module load SAMtools/1.21-GCC-13.3.0
module load angsd/0.940-GCC-12.3.0

WORKDIR=/scratch/eeg37520/nuc_120_mitomatch_redo
EXISTING_BAMLIST=/scratch/eeg37520/doliolid_popgen/dedupped_redo/angsd/onePerGene_markdup_top100.bamlist.txt
EXISTING_BAMDIR=/scratch/eeg37520/ReferenceVersionAssessments/onePerGene/picard_markdup
NEW_BAMDIR=${WORKDIR}/markdup
BAMLISTDIR=${WORKDIR}/bamlists
ANGSDDIR=${WORKDIR}/angsd
REF=/scratch/eeg37520/ReferenceVersionAssessments/onePerGene/index/final_ref.cds.onePerGene.fasta

mkdir -p ${BAMLISTDIR} ${ANGSDDIR}

echo "=== Step 1: Build 120-sample BAM list ==="
echo "Start: $(date)"

# Take the existing 100-sample list, remove DD_6_10
grep -v "DD_6_10" ${EXISTING_BAMLIST} > ${BAMLISTDIR}/nuc120_all.bamlist.txt

# Add the 6 DL5 samples (new BAMs from our mapping job)
for s in DL5_02 DL5_04 DL5_09 DL5_10 DL5_23 DL5_24; do
    echo "${NEW_BAMDIR}/${s}.opg.markdup.bam" >> ${BAMLISTDIR}/nuc120_all.bamlist.txt
done

# Add the 15 EG samples (new BAMs from our mapping job)
for s in EG1 EG2 EG3 EG5 EG6 EG7 EG8 EG9 EG10 EG11 EG12 EG13 EG14 EG15 EG16; do
    echo "${NEW_BAMDIR}/${s}.opg.markdup.bam" >> ${BAMLISTDIR}/nuc120_all.bamlist.txt
done

NBAMS=$(wc -l < ${BAMLISTDIR}/nuc120_all.bamlist.txt)
echo "BAM list has ${NBAMS} entries"

if [[ ${NBAMS} -ne 120 ]]; then
    echo "ERROR: Expected 120 BAMs, got ${NBAMS}"
    exit 1
fi

echo "=== Step 2: Verify all BAMs exist and are indexed ==="

MISSING=0
while read BAM; do
    if [[ ! -f "${BAM}" ]]; then
        echo "MISSING BAM: ${BAM}"
        MISSING=$((MISSING + 1))
    elif [[ ! -f "${BAM}.bai" ]]; then
        echo "MISSING INDEX: ${BAM}.bai — indexing now"
        samtools index ${BAM}
    fi
done < ${BAMLISTDIR}/nuc120_all.bamlist.txt

if [[ ${MISSING} -gt 0 ]]; then
    echo "ERROR: ${MISSING} BAM files missing. Cannot proceed."
    exit 1
fi

echo "All 120 BAMs verified."

echo "=== Step 3: Per-sample flagstat summary ==="

# Collect mapped read counts for downsampling decision
echo -e "sample\ttotal_reads\tmapped_reads\tmapped_pct\tproperly_paired" > ${BAMLISTDIR}/nuc120_flagstat_summary.tsv

while read BAM; do
    SAMPLE=$(basename ${BAM} .opg.markdup.bam)
    STATS=$(samtools flagstat ${BAM})
    TOTAL=$(echo "${STATS}" | head -1 | awk '{print $1}')
    MAPPED=$(echo "${STATS}" | grep "mapped (" | head -1 | awk '{print $1}')
    MAPPED_PCT=$(echo "${STATS}" | grep "mapped (" | head -1 | sed 's/.*(\([0-9.]*\)%.*/\1/')
    PAIRED=$(echo "${STATS}" | grep "properly paired" | awk '{print $1}')
    echo -e "${SAMPLE}\t${TOTAL}\t${MAPPED}\t${MAPPED_PCT}\t${PAIRED}" >> ${BAMLISTDIR}/nuc120_flagstat_summary.tsv
done < ${BAMLISTDIR}/nuc120_all.bamlist.txt

echo "Flagstat summary written to ${BAMLISTDIR}/nuc120_flagstat_summary.tsv"

echo "=== Step 4: ANGSD depth scan (all 120 samples) ==="

angsd -b ${BAMLISTDIR}/nuc120_all.bamlist.txt \
    -ref ${REF} \
    -out ${ANGSDDIR}/nuc120_depth \
    -doCounts 1 \
    -doDepth 1 \
    -maxDepth 200 \
    -minMapQ 30 \
    -minQ 20 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -nThreads 16

echo "Depth scan complete."
echo "=== Step 5: Summary stats ==="

# Quick look at the depth distribution
echo "Per-sample depth file: ${ANGSDDIR}/nuc120_depth.depthSample"
echo "Global depth file: ${ANGSDDIR}/nuc120_depth.depthGlobal"

# Print the top and bottom of the flagstat summary sorted by mapped reads
echo ""
echo "=== Samples by mapped read count (lowest 10) ==="
sort -t$'\t' -k3 -n ${BAMLISTDIR}/nuc120_flagstat_summary.tsv | head -11

echo ""
echo "=== Samples by mapped read count (highest 10) ==="
sort -t$'\t' -k3 -rn ${BAMLISTDIR}/nuc120_flagstat_summary.tsv | head -10

echo ""
echo "Done: $(date)"
