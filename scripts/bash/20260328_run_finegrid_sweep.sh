#!/bin/bash
#SBATCH --job-name=finegrid
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=168:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
#SBATCH --output=/scratch/eeg37520/doliolid_slim/logs/finegrid_%A_%a.out
#SBATCH --error=/scratch/eeg37520/doliolid_slim/logs/finegrid_%A_%a.err
#SBATCH --array=1-18

source /apps/eb/Miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate slim_env

cd /scratch/eeg37520/doliolid_slim

PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sweeps/20260327_11k_sweep_params.txt)
COMBO_ID=$(echo $PARAMS | awk '{print $1}')
K_NURSES=$(echo $PARAMS | awk '{print $2}')
OOZ_SURVIVAL=$(echo $PARAMS | awk '{print $3}')
NURSE_MORT=$(echo $PARAMS | awk '{print $4}')
MU=$(echo $PARAMS | awk '{print $5}')

echo "Task ${SLURM_ARRAY_TASK_ID}: combo=${COMBO_ID} K=${K_NURSES} ooz=${OOZ_SURVIVAL} mort=${NURSE_MORT} mu=${MU}"

slim \
  -d COMBO_ID=${COMBO_ID} \
  -d K_NURSES=${K_NURSES} \
  -d OOZ_SURVIVAL=${OOZ_SURVIVAL} \
  -d NURSE_MORT=${NURSE_MORT} \
  -d MU=${MU} \
  20260328_bloom_sweep_11k_finegrid.slim
