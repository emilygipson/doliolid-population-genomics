#!/bin/bash
#SBATCH --job-name=conspec_divstats
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/align_diversity_%j.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/align_diversity_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu

# Load modules — UPDATE MAFFT VERSION FROM module spider OUTPUT
module load MAFFT/7.520-GCC-12.3.0-with-extensions
module load R/4.3.2-foss-2023a

WORKDIR=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific
ASSEMBLYDIR=$WORKDIR/assemblies
OUTDIR=$WORKDIR/analysis
mkdir -p $OUTDIR

echo "=========================================="
echo "Conspecific assembly alignment & diversity"
echo "Start: $(date)"
echo "=========================================="

# Step 1: Collect all 105 noIUPAC fastas into one multi-fasta
echo "[$(date)] Collecting assemblies..."
> $OUTDIR/all_105_conspecific.fasta
for f in $(find $ASSEMBLYDIR -name "*noIUPAC.fasta" | sort); do
  sample=$(echo $f | sed 's|.*/assemblies/||' | cut -d'/' -f1)
  echo ">$sample" >> $OUTDIR/all_105_conspecific.fasta
  grep -v "^>" $f | tr -d '\n' >> $OUTDIR/all_105_conspecific.fasta
  echo "" >> $OUTDIR/all_105_conspecific.fasta
done
echo "[$(date)] Collected $(grep -c '^>' $OUTDIR/all_105_conspecific.fasta) sequences"

# Step 2: Align with MAFFT
echo "[$(date)] Running MAFFT alignment..."
mafft --auto --thread 8 $OUTDIR/all_105_conspecific.fasta > $OUTDIR/all_105_conspecific_aligned.fasta
echo "[$(date)] Alignment complete"

# Step 3: Diversity stats in R
echo "[$(date)] Computing diversity statistics..."
cat > $OUTDIR/diversity_comparison.R << 'RSCRIPT'
library(ape)
library(pegas)

aln <- read.FASTA("/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/analysis/all_105_conspecific_aligned.fasta")
cat("Sequences loaded:", length(aln), "\n")

# Alignment dimensions
aln_matrix <- as.matrix(as.DNAbin(aln))
cat("Alignment length:", ncol(aln_matrix), "\n")

# Nucleotide diversity
pi_val <- nuc.div(aln)
cat("Nucleotide diversity (pi):", pi_val, "\n")

# Haplotypes
h <- haplotype(aln)
nhaps <- nrow(h)
cat("Number of haplotypes:", nhaps, "\n")

# Haplotype diversity
hdiv <- hap.div(aln)
cat("Haplotype diversity (Hd):", hdiv, "\n")

# Tajima's D
seg <- length(seg.sites(aln))
cat("Segregating sites:", seg, "\n")
taj <- tajima.test(aln)
cat("Tajima's D:", taj$D, "\n")
cat("Tajima's D p-value:", taj$Pval.normal, "\n")

# Summary comparison table
cat("\n========================================\n")
cat("CONSPECIFIC vs ORIGINAL COMPARISON\n")
cat("========================================\n")
cat(sprintf("%-30s %-15s %-15s\n", "Metric", "Original", "Conspecific"))
cat(sprintf("%-30s %-15s %-15s\n", "-----", "--------", "-----------"))
cat(sprintf("%-30s %-15s %-15.5f\n", "Nucleotide diversity (pi)", "0.00397", pi_val))
cat(sprintf("%-30s %-15s %-15d\n", "Haplotypes", "89", nhaps))
cat(sprintf("%-30s %-15s %-15.3f\n", "Haplotype diversity (Hd)", "0.998", hdiv))
cat(sprintf("%-30s %-15s %-15.3f\n", "Tajima's D", "-2.300", taj$D))
cat(sprintf("%-30s %-15s %-15d\n", "Segregating sites", "117", seg))
cat("========================================\n")

sink(file = "/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/analysis/diversity_results.txt")
cat("Conspecific MITObim Diversity Results\n")
cat("Date:", format(Sys.time()), "\n\n")
cat(sprintf("Sequences: %d\n", length(aln)))
cat(sprintf("Alignment length: %d\n", ncol(aln_matrix)))
cat(sprintf("Segregating sites: %d\n", seg))
cat(sprintf("Nucleotide diversity (pi): %.5f\n", pi_val))
cat(sprintf("Haplotypes: %d\n", nhaps))
cat(sprintf("Haplotype diversity (Hd): %.3f\n", hdiv))
cat(sprintf("Tajima's D: %.3f\n", taj$D))
cat(sprintf("Tajima's D p-value: %.4f\n", taj$Pval.normal))
sink()

cat("\nResults saved to diversity_results.txt\n")
RSCRIPT

Rscript $OUTDIR/diversity_comparison.R

echo "=========================================="
echo "[$(date)] Done"
echo "=========================================="
