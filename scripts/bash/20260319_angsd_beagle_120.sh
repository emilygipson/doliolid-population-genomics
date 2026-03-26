#!/bin/bash
#SBATCH --job-name=beagle_120
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.out
#SBATCH --error=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
# =============================================================================
# 20260319_angsd_beagle_120.sh
#
# PURPOSE:
#   Generate beagle-format genotype likelihoods for PCAngsd PCA and
#   downstream analyses (kinship, admixture). Uses SNP_pval filter to
#   retain only polymorphic sites. minInd=50 for stricter missingness
#   filtering on the PCA to reduce coverage artifacts.
#
# DEPENDS ON: 20260319_build_bamlist_depth.sh (job 43719648)
#
# NOTE: Same downsampling caveat as the SAF script.
# =============================================================================
set -euo pipefail

module purge
module load angsd/0.940-GCC-12.3.0

WORKDIR=/scratch/eeg37520/nuc_120_mitomatch_redo
BAMLIST=${WORKDIR}/bamlists/nuc120_all.bamlist.txt
REF=/scratch/eeg37520/ReferenceVersionAssessments/onePerGene/index/final_ref.cds.onePerGene.fasta
OUTDIR=${WORKDIR}/angsd
OUTPRE=${OUTDIR}/nuc120_beagle_minInd50

mkdir -p ${OUTDIR}

NBAMS=$(wc -l < ${BAMLIST})
if [[ ${NBAMS} -ne 118 ]]; then
    echo "ERROR: BAM list has ${NBAMS} entries, expected 118"
    exit 1
fi

echo "=== ANGSD beagle GL: 118 samples ==="
echo "BAM list: ${BAMLIST} (${NBAMS} samples)"
echo "minInd: 50 (stricter for PCA)"
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
    -setMaxDepth 50 \
    -minMapQ 30 \
    -minQ 20 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -baq 1 \
    -minInd 50 \
    -SNP_pval 1e-6 \
    -minMaf 0.05 \
    -nThreads 16

echo "=== Beagle output ==="
echo "File: ${OUTPRE}.beagle.gz"
NSITES=$(zcat ${OUTPRE}.beagle.gz | tail -n +2 | wc -l)
echo "SNPs retained: ${NSITES}"
echo ""
echo "Done: $(date)"
