#!/bin/bash
#SBATCH --job-name=bsp40_rates_v2
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-2
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim/BEAST/rate_sensitivity/bsp40_rate_v2_%A_%a.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim/BEAST/rate_sensitivity/bsp40_rate_v2_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

module purge
module load Beast/2.7.7-GCC-12.3.0

cd /scratch/eeg37520/doliolid_popgen/mitobim/BEAST/rate_sensitivity

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
    RATE="1.5E-8"
    LABEL="slow"
elif [ $SLURM_ARRAY_TASK_ID -eq 2 ]; then
    RATE="6.0E-8"
    LABEL="fast"
fi

# Create XML with modified rate, unique filenames, and matched log frequencies
sed -e "s/3.0E-8/${RATE}/" \
    -e 's/fileName="subsample_40_aligned.log"/fileName="bsp40_'"${LABEL}"'.log"/' \
    -e 's/fileName="subsample_40_aligned"/fileName=""/' \
    -e 's/fileName="$(filebase)-$(tree).trees"/fileName="bsp40_'"${LABEL}"'.trees"/' \
    -e 's/logEvery="1000" mode="tree"/logEvery="5000" mode="tree"/' \
    ../dgeg_bsp_40.xml > dgeg_bsp_40_${LABEL}.xml

echo "=========================================="
echo "BEAST2 BSP 40-seq, rate = ${RATE} (${LABEL})"
echo "Start: $(date)"
echo "=========================================="

beast -overwrite -threads 8 dgeg_bsp_40_${LABEL}.xml

echo "=========================================="
echo "[$(date)] Done"
echo "=========================================="
