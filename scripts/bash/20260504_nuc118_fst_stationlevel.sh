#!/bin/bash
#SBATCH --job-name=nuc118_fst_stationlevel
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
#SBATCH --output=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/nuc118_fst_stationlevel_%j.out
#SBATCH --error=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/nuc118_fst_stationlevel_%j.err

set -euo pipefail

module load angsd/0.940-GCC-12.3.0

WORKDIR=/scratch/eeg37520/nuc_120_mitomatch_redo
BAMLIST_DIR=${WORKDIR}/bamlists
OUTDIR=${WORKDIR}/fst_stationlevel_20260504
mkdir -p ${OUTDIR}

BLOOM_BAMS=${BAMLIST_DIR}/nuc118_bloom_stationlevel_20260504.bamlist.txt
NONBLOOM_BAMS=${BAMLIST_DIR}/nuc118_nonbloom_stationlevel_20260504.bamlist.txt

REF=$(awk '{for(i=1;i<NF;i++) if($i=="-ref") {print $(i+1); exit}}' ${WORKDIR}/fst/bloom.arg)
echo "Reference (from original bloom.arg): ${REF}"
if [ ! -f "${REF}" ]; then
    echo "ERROR: reference file not found at ${REF}"
    exit 1
fi

N_BLOOM=$(wc -l < ${BLOOM_BAMS})
N_NONBLOOM=$(wc -l < ${NONBLOOM_BAMS})
MININD_BLOOM=$(( N_BLOOM / 5 ))
MININD_NONBLOOM=$(( N_NONBLOOM / 5 ))

echo "=== Station-level partition ==="
echo "Bloom n = ${N_BLOOM}, minInd = ${MININD_BLOOM}"
echo "Non-bloom n = ${N_NONBLOOM}, minInd = ${MININD_NONBLOOM}"
echo

echo "=== Step 1: bloom SAF ==="
angsd \
    -bam ${BLOOM_BAMS} \
    -ref ${REF} \
    -anc ${REF} \
    -out ${OUTDIR}/bloom_stationlevel \
    -GL 1 \
    -doSaf 1 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -doCounts 1 \
    -setMinDepth 2 \
    -setMaxDepth 50 \
    -minMapQ 30 \
    -minQ 20 \
    -baq 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -minInd ${MININD_BLOOM} \
    -nThreads 8

echo
echo "=== Step 2: non-bloom SAF ==="
angsd \
    -bam ${NONBLOOM_BAMS} \
    -ref ${REF} \
    -anc ${REF} \
    -out ${OUTDIR}/nonbloom_stationlevel \
    -GL 1 \
    -doSaf 1 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -doCounts 1 \
    -setMinDepth 2 \
    -setMaxDepth 50 \
    -minMapQ 30 \
    -minQ 20 \
    -baq 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -minInd ${MININD_NONBLOOM} \
    -nThreads 8

echo
echo "=== Step 3: per-group folded 1D SFS ==="
realSFS ${OUTDIR}/bloom_stationlevel.saf.idx -fold 1 -P 8 > ${OUTDIR}/bloom_stationlevel.sfs
realSFS ${OUTDIR}/nonbloom_stationlevel.saf.idx -fold 1 -P 8 > ${OUTDIR}/nonbloom_stationlevel.sfs

echo
echo "=== Step 4: joint 2D folded SFS ==="
realSFS \
    ${OUTDIR}/bloom_stationlevel.saf.idx \
    ${OUTDIR}/nonbloom_stationlevel.saf.idx \
    -fold 1 -P 8 \
    > ${OUTDIR}/bloom_nonbloom_stationlevel_2d.sfs

echo
echo "=== Step 5a: Fst index, whichFst 0 (Reynolds 1983) ==="
realSFS fst index \
    ${OUTDIR}/bloom_stationlevel.saf.idx \
    ${OUTDIR}/nonbloom_stationlevel.saf.idx \
    -sfs ${OUTDIR}/bloom_nonbloom_stationlevel_2d.sfs \
    -fold 1 \
    -fstout ${OUTDIR}/bloom_nonbloom_stationlevel_fst_whichFst0 \
    -whichFst 0

echo
echo "=== Step 5b: Fst index, whichFst 1 (Bhatia, recommended for small samples) ==="
realSFS fst index \
    ${OUTDIR}/bloom_stationlevel.saf.idx \
    ${OUTDIR}/nonbloom_stationlevel.saf.idx \
    -sfs ${OUTDIR}/bloom_nonbloom_stationlevel_2d.sfs \
    -fold 1 \
    -fstout ${OUTDIR}/bloom_nonbloom_stationlevel_fst_whichFst1 \
    -whichFst 1

echo
echo "=== Step 6a: Fst stats, whichFst 0 ==="
realSFS fst stats ${OUTDIR}/bloom_nonbloom_stationlevel_fst_whichFst0.fst.idx \
    > ${OUTDIR}/bloom_nonbloom_stationlevel_fst_whichFst0_stats.txt 2>&1

echo "=== Step 6b: Fst stats, whichFst 1 ==="
realSFS fst stats ${OUTDIR}/bloom_nonbloom_stationlevel_fst_whichFst1.fst.idx \
    > ${OUTDIR}/bloom_nonbloom_stationlevel_fst_whichFst1_stats.txt 2>&1

echo
echo "=== HEADLINE RESULTS ==="
echo "Bloom n = ${N_BLOOM}, Non-bloom n = ${N_NONBLOOM}"
echo
echo "--- whichFst 0 (Reynolds 1983) ---"
cat ${OUTDIR}/bloom_nonbloom_stationlevel_fst_whichFst0_stats.txt
echo
echo "--- whichFst 1 (Bhatia) ---"
cat ${OUTDIR}/bloom_nonbloom_stationlevel_fst_whichFst1_stats.txt
echo
echo "All outputs in ${OUTDIR}/"
ls -la ${OUTDIR}/
