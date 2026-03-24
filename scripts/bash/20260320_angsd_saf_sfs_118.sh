#!/bin/bash
#SBATCH --job-name=saf_sfs_118
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
# 20260320_angsd_saf_sfs_118.sh
#
# PURPOSE: Global SAF + folded SFS + diversity stats for 118 samples
# FIX: setMaxDepth=300 (was 50, matching old 100-sample runs)
# =============================================================================
set -euo pipefail

module purge
module load angsd/0.940-GCC-12.3.0

WORKDIR=/scratch/eeg37520/nuc_120_mitomatch_redo
BAMLIST=${WORKDIR}/bamlists/nuc120_all.bamlist.txt
REF=/scratch/eeg37520/ReferenceVersionAssessments/onePerGene/index/final_ref.cds.onePerGene.fasta
OUTDIR=${WORKDIR}/angsd
OUTPRE=${OUTDIR}/nuc118

mkdir -p ${OUTDIR}

NSAMP=$(wc -l < ${BAMLIST})
MININD=$((NSAMP * 20 / 100))

echo "=== ANGSD SAF + SFS: ${NSAMP} samples ==="
echo "BAM list: ${BAMLIST} (${NSAMP} samples)"
echo "minInd: ${MININD} (20% of ${NSAMP})"
echo "setMaxDepth: 300"
echo "Start: $(date)"

# --- SAF estimation ---
angsd -b ${BAMLIST} \
    -ref ${REF} \
    -out ${OUTPRE}_saf \
    -GL 1 \
    -doSaf 1 \
    -anc ${REF} \
    -doCounts 1 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -setMinDepth 2 \
    -setMaxDepth 300 \
    -minMapQ 30 \
    -minQ 20 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -baq 1 \
    -minInd ${MININD} \
    -nThreads 16

echo "SAF complete: $(date)"

# --- Folded SFS ---
realSFS ${OUTPRE}_saf.saf.idx -fold 1 -P 16 > ${OUTPRE}_global.sfs
echo "SFS complete: $(date)"

# --- Diversity stats ---
realSFS saf2theta ${OUTPRE}_saf.saf.idx \
    -sfs ${OUTPRE}_global.sfs \
    -fold 1 \
    -outname ${OUTPRE}_thetas

thetaStat do_stat ${OUTPRE}_thetas.thetas.idx

echo ""
echo "=== Diversity results ==="
cat ${OUTPRE}_thetas.thetas.idx.pestPG

echo ""
echo "=== SFS ==="
cat ${OUTPRE}_global.sfs

echo ""
echo "Done: $(date)"
