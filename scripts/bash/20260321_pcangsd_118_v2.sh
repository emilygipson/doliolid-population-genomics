#!/bin/bash
#SBATCH --job-name=pcangsd_118
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.out
#SBATCH --error=/scratch/eeg37520/nuc_120_mitomatch_redo/logs/%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
# =============================================================================
# 20260321_pcangsd_118_v2.sh
#
# PURPOSE:
#   Rerun PCAngsd with increased iteration limit (500) to ensure convergence.
#   v1 hit the 100-iteration default and did not converge (RMSE=0.000288).
# =============================================================================
set -euo pipefail

module purge
module load PCAngsd/1.2-gfbf-2023a

WORKDIR=/scratch/eeg37520/nuc_120_mitomatch_redo
BEAGLE=${WORKDIR}/angsd/nuc118_beagle_minInd24.beagle.gz
OUTDIR=${WORKDIR}/pcangsd
OUTPRE=${OUTDIR}/nuc118_pcangsd_minInd24_v2

mkdir -p ${OUTDIR}

echo "=== PCAngsd v2: 118 samples, --iter 500 ==="
echo "Input: ${BEAGLE}"
echo "Start: $(date)"

pcangsd -b ${BEAGLE} \
    -o ${OUTPRE} \
    -t 8 \
    --iter 500

echo "=== Output ==="
echo "Covariance matrix: ${OUTPRE}.cov"
echo ""
echo "Done: $(date)"
