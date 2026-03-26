#!/bin/bash
#SBATCH --job-name=mitobim_conspec
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-105
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/mitobim_%A_%a.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/mitobim_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

module purge
module load MITObim/1.9.1

# Create logs directory
mkdir -p /scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs

cd /scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific

# Get sample list (105 samples, excluding DL5_25)
SAMPLE_LIST=/scratch/eeg37520/doliolid_popgen/mitobim/sample_names_105.txt
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "ERROR: Sample list not found at $SAMPLE_LIST"
    exit 1
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)
if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "=========================================="
echo "MITObim assembly: $SAMPLE"
echo "Reference: DD_21_05 (conspecific)"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Start: $(date)"
echo "=========================================="

# Conspecific reference
REF_FASTA=/scratch/eeg37520/doliolid_popgen/mitobim/getorg_complete/DD_21_05_mt.fasta

# Input reads
READ_DIR=/scratch/eeg37520/doliolid_popgen/mitobim/mt_reads/fastq
R1=${READ_DIR}/${SAMPLE}_R1.fastq.gz
R2=${READ_DIR}/${SAMPLE}_R2.fastq.gz
SINGLETON=${READ_DIR}/${SAMPLE}_singleton.fastq.gz

# Check input files exist
if [ ! -f "$R1" ] || [ ! -f "$R2" ] || [ ! -f "$SINGLETON" ]; then
    echo "ERROR: Missing input files for $SAMPLE"
    echo "R1: $R1"
    echo "R2: $R2" 
    echo "Singleton: $SINGLETON"
    exit 1
fi

# Create sample directory
mkdir -p ${SAMPLE}
cd ${SAMPLE}

# Run MITObim
MITObim.pl -sample $SAMPLE \
          -ref DD_21_05_ref \
          -readpool $R1 $R2 $SINGLETON \
          -quick $REF_FASTA \
          --clean \
          --end 50

echo "=========================================="
echo "[$SAMPLE] $(date) Done"
echo "=========================================="
