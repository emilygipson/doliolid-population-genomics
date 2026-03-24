#!/bin/bash
#SBATCH --job-name=bloom_fst_118
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.out
#SBATCH --error=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
# =============================================================================
# 20260320_bloom_fst_118.sh
#
# PURPOSE: Bloom vs non-bloom Fst for 118-sample nuclear dataset
# FIX: setMaxDepth=300 (was 50), fixed BAM list construction to handle
#      .opg.markdup.ds.bam suffix from downsampled EG samples
#
# BLOOM (n=80): DD7, DD11, DD17, DD21, DM18, EG2,6,7,8,9,10,11,12,13,14,15,16
# NON-BLOOM (n=38): DD4, DD6, DD8, DD9, DD10, DD15, DD19, DL5, EG1,3,5
# =============================================================================
set -euo pipefail

module purge
module load angsd/0.940-GCC-12.3.0

WORKDIR=/scratch/eeg37520/nuc_120_mitomatch_redo
BAMLIST=${WORKDIR}/bamlists/nuc120_all.bamlist.txt
REF=/scratch/eeg37520/ReferenceVersionAssessments/onePerGene/index/final_ref.cds.onePerGene.fasta
FSTDIR=${WORKDIR}/fst
BAMLISTDIR=${WORKDIR}/bamlists

mkdir -p ${FSTDIR}

echo "=== Step 1: Build bloom and non-bloom BAM lists ==="
echo "Start: $(date)"

> ${BAMLISTDIR}/nuc118_bloom.bamlist.txt
> ${BAMLISTDIR}/nuc118_nonbloom.bamlist.txt

while read BAM; do
    # Strip path and handle both suffixes:
    #   .opg.markdup.ds.bam (downsampled)
    #   .opg.markdup.bam (standard)
    SAMPLE=$(basename "${BAM}")
    SAMPLE=${SAMPLE%.bam}
    SAMPLE=${SAMPLE%.ds}
    SAMPLE=${SAMPLE%.opg.markdup}

    BLOOM=0

    case ${SAMPLE} in
        DD_7_*|DD_11_*|DD_17_*|DD_21_*) BLOOM=1 ;;
        DD_4_*|DD_6_*|DD_8_*|DD_9_*|DD_10_*|DD_15_*|DD_19_*) BLOOM=0 ;;
        DM18_*) BLOOM=1 ;;
        DL5_*) BLOOM=0 ;;
        EG2|EG6|EG7|EG8|EG9|EG10|EG11|EG12|EG13|EG14|EG15|EG16) BLOOM=1 ;;
        EG1|EG3|EG5) BLOOM=0 ;;
        *) echo "WARNING: Unknown sample ${SAMPLE}, skipping"; continue ;;
    esac

    if [[ ${BLOOM} -eq 1 ]]; then
        echo "${BAM}" >> ${BAMLISTDIR}/nuc118_bloom.bamlist.txt
    else
        echo "${BAM}" >> ${BAMLISTDIR}/nuc118_nonbloom.bamlist.txt
    fi
done < ${BAMLIST}

NBLOOM=$(wc -l < ${BAMLISTDIR}/nuc118_bloom.bamlist.txt)
NNONBLOOM=$(wc -l < ${BAMLISTDIR}/nuc118_nonbloom.bamlist.txt)
echo "Bloom samples: ${NBLOOM}"
echo "Non-bloom samples: ${NNONBLOOM}"

if [[ ${NBLOOM} -ne 80 ]]; then
    echo "ERROR: Expected 80 bloom, got ${NBLOOM}. Aborting."
    exit 1
fi
if [[ ${NNONBLOOM} -ne 38 ]]; then
    echo "ERROR: Expected 38 non-bloom, got ${NNONBLOOM}. Aborting."
    exit 1
fi

MININD_BLOOM=$((NBLOOM * 20 / 100))
MININD_NONBLOOM=$((NNONBLOOM * 20 / 100))
echo "minInd bloom: ${MININD_BLOOM} (20% of ${NBLOOM})"
echo "minInd non-bloom: ${MININD_NONBLOOM} (20% of ${NNONBLOOM})"
echo "setMaxDepth: 300"

echo "=== Step 2: SAF estimation — bloom group ==="
echo "$(date)"

angsd -b ${BAMLISTDIR}/nuc118_bloom.bamlist.txt \
    -ref ${REF} \
    -out ${FSTDIR}/bloom \
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
    -minInd ${MININD_BLOOM} \
    -nThreads 16

echo "Bloom SAF complete: $(date)"

echo "=== Step 3: SAF estimation — non-bloom group ==="
echo "$(date)"

angsd -b ${BAMLISTDIR}/nuc118_nonbloom.bamlist.txt \
    -ref ${REF} \
    -out ${FSTDIR}/nonbloom \
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
    -minInd ${MININD_NONBLOOM} \
    -nThreads 16

echo "Non-bloom SAF complete: $(date)"

echo "=== Step 4: Per-group 1D folded SFS ==="

realSFS ${FSTDIR}/bloom.saf.idx -fold 1 -P 16 > ${FSTDIR}/bloom.sfs
realSFS ${FSTDIR}/nonbloom.saf.idx -fold 1 -P 16 > ${FSTDIR}/nonbloom.sfs

echo "1D SFS complete: $(date)"

echo "=== Step 5: 2D SFS (bloom x non-bloom) ==="

realSFS ${FSTDIR}/bloom.saf.idx ${FSTDIR}/nonbloom.saf.idx \
    -fold 1 -P 16 > ${FSTDIR}/bloom_nonbloom_2d.sfs

echo "2D SFS complete: $(date)"

echo "=== Step 6: Fst estimation ==="

realSFS fst index ${FSTDIR}/bloom.saf.idx ${FSTDIR}/nonbloom.saf.idx \
    -sfs ${FSTDIR}/bloom_nonbloom_2d.sfs \
    -fold 1 \
    -fstout ${FSTDIR}/bloom_nonbloom_fst

realSFS fst stats ${FSTDIR}/bloom_nonbloom_fst.fst.idx

echo ""
echo "=== Per-group diversity stats ==="

echo "--- Bloom group ---"
realSFS saf2theta ${FSTDIR}/bloom.saf.idx \
    -sfs ${FSTDIR}/bloom.sfs \
    -fold 1 \
    -outname ${FSTDIR}/bloom_thetas
thetaStat do_stat ${FSTDIR}/bloom_thetas.thetas.idx
cat ${FSTDIR}/bloom_thetas.thetas.idx.pestPG

echo ""
echo "--- Non-bloom group ---"
realSFS saf2theta ${FSTDIR}/nonbloom.saf.idx \
    -sfs ${FSTDIR}/nonbloom.sfs \
    -fold 1 \
    -outname ${FSTDIR}/nonbloom_thetas
thetaStat do_stat ${FSTDIR}/nonbloom_thetas.thetas.idx
cat ${FSTDIR}/nonbloom_thetas.thetas.idx.pestPG

echo ""
echo "Done: $(date)"
