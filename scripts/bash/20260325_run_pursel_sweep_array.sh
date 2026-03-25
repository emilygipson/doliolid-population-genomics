#!/bin/bash
#SBATCH --job-name=pursel_sweep
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --array=1-11
#SBATCH --output=/scratch/eeg37520/doliolid_slim/logs/20260325_pursel_sweep_%a.out
#SBATCH --error=/scratch/eeg37520/doliolid_slim/logs/20260325_pursel_sweep_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

# =============================================================================
# PURIFYING SELECTION BLOOM SWEEP — 11 DFE parameter combinations
# =============================================================================
# All at combo 4 life cycle: K=5000, OOZ_SURVIVAL=0.50, NURSE_MORT=0.05
# MU values must be calibrated before submission.
#
# Parameter file format (space-delimited, no header):
#   PURSEL_ID  DFE_TYPE  SEL_COEFF  GAMMA_SHAPE  MU
#
# Emily Gipson, UGA mfflab
# Created: 2026-03-25
# =============================================================================

cd /scratch/eeg37520/doliolid_slim

source /apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate slim_env

mkdir -p /scratch/eeg37520/doliolid_slim/logs
mkdir -p /scratch/eeg37520/doliolid_slim/sweeps

# Read parameters for this array task
PARAMFILE=/scratch/eeg37520/doliolid_slim/pursel_calibration/20260325_pursel_sweep_params.txt
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMFILE)

PURSEL_ID=$(echo $LINE | awk '{print $1}')
DFE_TYPE=$(echo $LINE | awk '{print $2}')
SEL_COEFF=$(echo $LINE | awk '{print $3}')
GAMMA_SHAPE=$(echo $LINE | awk '{print $4}')
MU=$(echo $LINE | awk '{print $5}')

echo "========================================"
echo "  Pursel Bloom Sweep: ID ${PURSEL_ID}"
echo "  DFE_TYPE=${DFE_TYPE}"
echo "  SEL_COEFF=${SEL_COEFF}"
echo "  GAMMA_SHAPE=${GAMMA_SHAPE}"
echo "  MU=${MU}"
echo "  Array task: ${SLURM_ARRAY_TASK_ID}"
echo "  Job ID: ${SLURM_JOB_ID}"
echo "  Start: $(date)"
echo "========================================"

# Build the SLiM command
SLIM_CMD="slim -d PURSEL_ID=${PURSEL_ID} -d \"DFE_TYPE='${DFE_TYPE}'\" -d SEL_COEFF=${SEL_COEFF} -d MU=${MU}"

# Add GAMMA_SHAPE only for gamma DFE
if [ "$DFE_TYPE" = "gamma" ]; then
    SLIM_CMD="${SLIM_CMD} -d GAMMA_SHAPE=${GAMMA_SHAPE}"
fi

SLIM_CMD="${SLIM_CMD} /scratch/eeg37520/doliolid_slim/20260325_pursel_bloom_sweep.slim"

echo "Running: ${SLIM_CMD}"
eval ${SLIM_CMD}

echo "========================================"
echo "  Pursel sweep ${PURSEL_ID} complete: $(date)"
echo "========================================"
