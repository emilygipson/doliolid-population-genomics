#!/bin/bash
#SBATCH --job-name=cal_11k
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --array=1-18
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
#SBATCH --output=/scratch/eeg37520/doliolid_slim/logs/20260327_11k_cal_%a.out
#SBATCH --error=/scratch/eeg37520/doliolid_slim/logs/20260327_11k_cal_%a.err

# =============================================================================
# NEUTRAL 11K ADAPTIVE CALIBRATION — 18 parameter combinations
# =============================================================================

cd /scratch/eeg37520/doliolid_slim

source /apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate slim_env

PARAMFILE=/scratch/eeg37520/doliolid_slim/20260327_11k_calibration_params.txt

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMFILE)
COMBO_ID=$(echo $LINE | awk '{print $1}')
K_NURSES=$(echo $LINE | awk '{print $2}')
OOZ_SURV=$(echo $LINE | awk '{print $3}')
NURSE_MORT=$(echo $LINE | awk '{print $4}')
MU=$(echo $LINE | awk '{print $5}')

echo "=== 11K Calibration: combo ${COMBO_ID} ==="
echo "  K=${K_NURSES} OOZ=${OOZ_SURV} MORT=${NURSE_MORT} MU=${MU}"

slim \
    -d K_NURSES=${K_NURSES} \
    -d MU=${MU} \
    -d OOZ_SURVIVAL=${OOZ_SURV} \
    -d NURSE_MORTALITY=${NURSE_MORT} \
    -d COMBO_ID=${COMBO_ID} \
    -d "LOGDIR='/scratch/eeg37520/doliolid_slim/calibration_logs'" \
    /scratch/eeg37520/doliolid_slim/20260327_neutral_calibration_11k.slim

echo "=== Combo ${COMBO_ID} calibration complete ==="
