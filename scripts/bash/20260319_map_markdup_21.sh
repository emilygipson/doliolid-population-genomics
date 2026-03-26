#!/bin/bash
#SBATCH --job-name=map_markdup_21
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-21
#SBATCH --output=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%A_%a.out
#SBATCH --error=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
# =============================================================================
# 20260319_map_markdup_21.sh
#
# PURPOSE:
#   Map 21 new samples (6 DL5 + 15 EG) to onePerGene reference with BWA-MEM,
#   then run Picard MarkDuplicates. Matches the pipeline used for the existing
#   100 samples: bwa mem → samtools sort → samtools fixmate → samtools sort →
#   Picard MarkDuplicates (REMOVE_DUPLICATES=true).
#
# INPUT:
#   Sample list: /scratch/eeg37520/nuc_120_mitomatch_redo/metadata/new21_samples.txt
#   Format: SAMPLE_ID<tab>R1_PATH<tab>R2_PATH
#
# OUTPUT:
#   Sorted BAMs: /scratch/eeg37520/nuc_120_mitomatch_redo/mapping/
#   Markdup BAMs: /scratch/eeg37520/nuc_120_mitomatch_redo/markdup/
# =============================================================================
set -euo pipefail

module purge
module load BWA/0.7.18-GCCcore-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load picard/3.3.0-Java-17

REF=/scratch/eeg37520/ReferenceVersionAssessments/onePerGene/index/final_ref.cds.onePerGene.fasta
SAMPLELIST=/scratch/eeg37520/nuc_120_mitomatch_redo/metadata/new21_samples.txt
MAPDIR=/scratch/eeg37520/nuc_120_mitomatch_redo/mapping
MDDIR=/scratch/eeg37520/nuc_120_mitomatch_redo/markdup

mkdir -p ${MAPDIR} ${MDDIR}

# Parse sample info for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLELIST} | cut -f1)
R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLELIST} | cut -f2)
R2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLELIST} | cut -f3)

echo "=== Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE} ==="
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Start: $(date)"

# Verify input files exist
if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
    echo "ERROR: Input files not found for ${SAMPLE}"
    exit 1
fi

# Step 1: BWA-MEM → coordinate-sorted BAM
echo "--- Step 1: BWA-MEM + sort ---"
bwa mem -t 4 ${REF} ${R1} ${R2} | \
    samtools sort -@ 4 -o ${MAPDIR}/${SAMPLE}.onePerGene.sorted.bam

# Step 2: Name sort for fixmate
echo "--- Step 2: Name sort ---"
samtools sort -n -@ 4 \
    -o ${MDDIR}/${SAMPLE}.namesorted.bam \
    ${MAPDIR}/${SAMPLE}.onePerGene.sorted.bam

# Step 3: Fixmate (add mate score tags)
echo "--- Step 3: Fixmate ---"
samtools fixmate -m \
    ${MDDIR}/${SAMPLE}.namesorted.bam \
    ${MDDIR}/${SAMPLE}.fixmate.bam

# Step 4: Position sort after fixmate
echo "--- Step 4: Position sort ---"
samtools sort -@ 4 \
    -o ${MDDIR}/${SAMPLE}.possorted.bam \
    ${MDDIR}/${SAMPLE}.fixmate.bam

# Step 5: Picard MarkDuplicates
echo "--- Step 5: Picard MarkDuplicates ---"
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=${MDDIR}/${SAMPLE}.possorted.bam \
    O=${MDDIR}/${SAMPLE}.opg.markdup.bam \
    M=${MDDIR}/${SAMPLE}.markdup.metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=LENIENT

# Step 6: Index final BAM
echo "--- Step 6: Index ---"
samtools index ${MDDIR}/${SAMPLE}.opg.markdup.bam

# Step 7: Flagstat for QC
echo "--- Step 7: Flagstat ---"
samtools flagstat ${MDDIR}/${SAMPLE}.opg.markdup.bam > ${MDDIR}/${SAMPLE}.flagstat.txt

# Clean up intermediates
rm -f ${MDDIR}/${SAMPLE}.namesorted.bam \
      ${MDDIR}/${SAMPLE}.fixmate.bam \
      ${MDDIR}/${SAMPLE}.possorted.bam

echo "Done: ${SAMPLE}"
echo "End: $(date)"
