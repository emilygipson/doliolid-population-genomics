#!/bin/bash
#SBATCH --job-name=getorg_all
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-109
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim/logs/getorg_%A_%a.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim/logs/getorg_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

module purge
module load GetOrganelle/1.7.7.1-foss-2023a

READDIR=/scratch/eeg37520/Oct2025_133_proc/clean_dedup_files
OUTBASE=/scratch/eeg37520/doliolid_popgen/mitobim/annotation_qc/getorganelle_all
SEED=/scratch/eeg37520/doliolid_popgen/mitobim/annotation_qc/cox1_seed.fasta
SAMPLE_LIST=/scratch/eeg37520/doliolid_popgen/mitobim/annotation_qc/getorg_sample_list.txt

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)
R1=${READDIR}/${SAMPLE}_unclassified__1.fq
R2=${READDIR}/${SAMPLE}_unclassified__2.fq
OUTDIR=${OUTBASE}/${SAMPLE}

echo "=========================================="
echo "GetOrganelle: ${SAMPLE}"
echo "Task: ${SLURM_ARRAY_TASK_ID} of 109"
echo "Start: $(date)"
echo "=========================================="

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: Read files not found for ${SAMPLE}"
    exit 1
fi

get_organelle_from_reads.py \
    -1 $R1 \
    -2 $R2 \
    -o $OUTDIR \
    -t 8 \
    -s $SEED \
    -F animal_mt \
    -R 20 \
    --max-reads 5E7

echo "=========================================="
echo "[$(date)] Done: ${SAMPLE}"
echo "=========================================="
