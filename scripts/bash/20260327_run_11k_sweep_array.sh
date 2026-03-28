#!/bin/bash
#SBATCH --job-name=sweep_11k
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --array=1-18
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
#SBATCH --output=/scratch/eeg37520/doliolid_slim/logs/20260327_11k_sweep_%a.out
#SBATCH --error=/scratch/eeg37520/doliolid_slim/logs/20260327_11k_sweep_%a.err

# =============================================================================
# 11K NEUTRAL BLOOM SWEEP — 18 parameter combinations, 7 sample sizes
# =============================================================================
# DO NOT SUBMIT until 11k calibrations complete and sweep param file is built.
#
# Emily Gipson, UGA mfflab
# Created: 2026-03-27
# =============================================================================

cd /scratch/eeg37520/doliolid_slim

source /apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate slim_env

# Parameter file: combo_id K_NURSES ooz_survival nurse_mortality CALIBRATED_mu
PARAMFILE=/scratch/eeg37520/doliolid_slim/sweeps/20260327_11k_sweep_params.txt

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMFILE)
COMBO_ID=$(echo $LINE | awk '{print $1}')
K_NURSES=$(echo $LINE | awk '{print $2}')
OOZ_SURV=$(echo $LINE | awk '{print $3}')
NURSE_MORT=$(echo $LINE | awk '{print $4}')
MU=$(echo $LINE | awk '{print $5}')

echo "=== 11K Sweep: combo ${COMBO_ID} ==="
echo "  K_NURSES=${K_NURSES} OOZ_SURV=${OOZ_SURV} NURSE_MORT=${NURSE_MORT} MU=${MU}"

slim \
    -d K_NURSES=${K_NURSES} \
    -d MU=${MU} \
    -d OOZ_SURVIVAL=${OOZ_SURV} \
    -d NURSE_MORT=${NURSE_MORT} \
    -d COMBO_ID=${COMBO_ID} \
    -d N_REPS=500 \
    /scratch/eeg37520/doliolid_slim/20260327_bloom_sweep_11k.slim

echo "=== 11K Combo ${COMBO_ID} sweep complete ==="
