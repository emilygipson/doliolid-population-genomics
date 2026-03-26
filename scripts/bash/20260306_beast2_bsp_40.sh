#!/bin/bash
#SBATCH --job-name=beast2_bsp40
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim/BEAST/logs/beast2_bsp40_%j.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim/BEAST/logs/beast2_bsp40_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

module purge
module load Beast/2.7.7-GCC-12.3.0

cd /scratch/eeg37520/doliolid_popgen/mitobim/BEAST

echo "=========================================="
echo "BEAST2 BSP — 40 sequence subsample"
echo "Model: HKY+G4, strict clock, rate=3.0e-8"
echo "Chain: 50M generations"
echo "Start: $(date)"
echo "=========================================="

beast -threads 8 dgeg_bsp_40.xml

echo "=========================================="
echo "[$(date)] Done"
echo "=========================================="
