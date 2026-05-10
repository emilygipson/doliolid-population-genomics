#!/bin/bash
#SBATCH --job-name=phi_test_119
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim/BEAST/may2026_mito119_BSP_redo/logs/phi_test_119_%j.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim/BEAST/may2026_mito119_BSP_redo/logs/phi_test_119_%j.err

set -euo pipefail

module load PhiPack/2016.06.14-GCC-13.2.0

WORKDIR=/scratch/eeg37520/doliolid_popgen/mitobim/BEAST/may2026_mito119_BSP_redo
ALN=${WORKDIR}/all_119_mitobim_filtered.fasta
OUTDIR=${WORKDIR}/phi_test_20260510

mkdir -p ${OUTDIR}
mkdir -p ${WORKDIR}/logs
cd ${OUTDIR}

# -f input fasta
# -p 1000 permutations for PHI p-value
# -o also report NSS and Max-Chi^2 (bonus diagnostics)
# -v verbose
Phi -f ${ALN} -p 1000 -o -v > phi_results_119_20260510.txt 2>&1

echo "PHI test done. Results in ${OUTDIR}/phi_results_119_20260510.txt"
