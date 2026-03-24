#!/bin/bash
#SBATCH --job-name=beagle_118
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.out
#SBATCH --error=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
# =============================================================================
# 20260320_angsd_beagle_118.sh
#
# PURPOSE: Beagle genotype likelihoods for PCAngsd PCA + kinship
# FIX: setMaxDepth=300 (was 50), minInd=24 (was 50)
# Matches parameters from old working run (minInd=50 for 100 samples
# was ~50%; 24 for 118 = ~20%. Old run also used setMaxDepth=300.)
# =============================================================================
set -euo pipefail

module purge
module load angsd/0.940-GCC-12.3.0

WORKDIR=/scratch/eeg37520/nuc_120_mitomatch_redo
BAMLIST=${WORKDIR}/bamlists/nuc120_all.bamlist.txt
REF=/scratch/eeg37520/ReferenceVersionAssessments/onePerGene/index/final_ref.cds.onePerGene.fasta
OUTDIR=${WORKDIR}/angsd
OUTPRE=${OUTDIR}/nuc118_beagle_minInd24

mkdir -p ${OUTDIR}

NSAMP=$(wc -l < ${BAMLIST})
echo "=== ANGSD beagle GL: ${NSAMP} samples ==="
echo "BAM list: ${BAMLIST} (${NSAMP} samples)"
echo "minInd: 24 (20% of 118)"
echo "setMaxDepth: 300"
echo "Start: $(date)"

angsd -b ${BAMLIST} \
    -ref ${REF} \
    -out ${OUTPRE} \
    -GL 1 \
    -doGlf 2 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -doCounts 1 \
    -setMinDepth 2 \
    -setMaxDepth 300 \
    -minMapQ 30 \
    -minQ 20 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -baq 1 \
    -minInd 24 \
    -SNP_pval 1e-6 \
    -minMaf 0.05 \
    -nThreads 16

echo ""
echo "=== Beagle output ==="
echo "File: ${OUTPRE}.beagle.gz"
SNP_COUNT=$(zcat ${OUTPRE}.beagle.gz | tail -n +2 | wc -l)
echo "SNPs retained: ${SNP_COUNT}"

echo ""
echo "Done: $(date)"
