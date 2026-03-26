## 20260319_120sample_diversity_stats.R
## Diversity statistics on the expanded 120-sample conspecific MITObim alignment
## Requires: all_120_mitobim_aligned.fasta (from job 43694746)
## Run locally in RStudio after downloading alignment from cluster

library(ape)
library(pegas)

## ---- Read alignment ----
aln <- read.FASTA("all_120_mitobim_aligned.fasta")
cat("Samples read:", length(aln), "\n")
cat("Sample names:\n")
print(labels(aln))

## ---- Check alignment dimensions ----
mat <- as.character(as.matrix(aln))
cat("Raw alignment: ", nrow(mat), "samples x", ncol(mat), "columns\n")

## ---- Strict filtering: keep only columns with A/C/G/T in ALL samples ----
valid <- c("a", "c", "g", "t")
keep <- apply(mat, 2, function(col) all(col %in% valid))
cat("Columns retained:", sum(keep), "of", ncol(mat), 
    "(", round(100 * sum(keep) / ncol(mat), 1), "%)\n")

mat_filt <- mat[, keep]
cat("Filtered alignment:", nrow(mat_filt), "x", ncol(mat_filt), "\n")

## ---- Write filtered alignment ----
aln_filt <- as.DNAbin(mat_filt)
write.FASTA(aln_filt, "all_120_mitobim_filtered.fasta")
cat("Filtered alignment written to all_120_mitobim_filtered.fasta\n")

## ---- Basic diversity stats ----
## Exclude DM18_27 for main stats (divergent lineage)
exclude <- "DM18_27"
keep_samples <- !labels(aln_filt) %in% exclude
aln_119 <- aln_filt[keep_samples]
cat("\nStats computed on", length(aln_119), "samples (excluding", exclude, ")\n")

## Segregating sites
seg <- seg.sites(aln_119)
cat("Segregating sites:", length(seg), "\n")

## Haplotypes
h <- haplotype(aln_119)
nhaps <- nrow(h)
cat("Unique haplotypes:", nhaps, "/", length(aln_119), "\n")

## Haplotype diversity
hd <- hap.div(aln_119)
cat("Haplotype diversity (Hd):", round(hd, 4), "\n")

## Nucleotide diversity
pi_val <- nuc.div(aln_119)
cat("Nucleotide diversity (pi):", round(pi_val, 5), "\n")

## Tajima's D
tajd <- tajima.test(aln_119)
cat("Tajima's D:", round(tajd$D, 3), " p =", round(tajd$Pval.normal, 4), "\n")

## Fu's Fs
## pegas doesn't have Fu's Fs directly — compute from haplotype frequencies
## Using the R.utils approach
cat("\n--- Fu's Fs requires manual computation or external tool ---\n")
cat("Use DnaSP or compute from haplotype frequencies\n")

## Pairwise differences
d <- dist.dna(aln_119, model = "N", pairwise.deletion = TRUE)
cat("\nPairwise nucleotide differences:\n")
cat("  Mean:", round(mean(d), 2), "\n")
cat("  Min:", min(d), "\n")
cat("  Max:", max(d), "\n")
cat("  Median:", median(d), "\n")

## Singletons
mat_119 <- as.character(as.matrix(aln_119))
seg_cols <- seg.sites(aln_119)
n_singletons <- 0
for (s in seg_cols) {
    col <- mat_119[, s]
    tab <- table(col)
    if (any(tab == 1)) n_singletons <- n_singletons + 1
}
cat("\nSingletons:", n_singletons, "of", length(seg_cols), 
    "(", round(100 * n_singletons / length(seg_cols), 1), "%)\n")

## ---- Summary table ----
cat("\n\n===== TAB-DELIMITED SUMMARY (paste into Excel) =====\n")
cat("Metric\tValue\n")
cat("Samples\t", length(aln_119), "\n", sep = "")
cat("Alignment length (filtered)\t", ncol(mat_filt), "\n", sep = "")
cat("Segregating sites\t", length(seg), "\n", sep = "")
cat("Singletons\t", n_singletons, " (", round(100 * n_singletons / length(seg), 1), "%)\n", sep = "")
cat("Unique haplotypes\t", nhaps, "/", length(aln_119), "\n", sep = "")
cat("Hd\t", round(hd, 4), "\n", sep = "")
cat("Pi\t", round(pi_val, 5), "\n", sep = "")
cat("Tajima's D\t", round(tajd$D, 3), " (p = ", round(tajd$Pval.normal, 4), ")\n", sep = "")
cat("Mean pairwise differences\t", round(mean(d), 2), "\n", sep = "")
cat("Min pairwise differences\t", min(d), "\n", sep = "")

## ---- Also compute stats on full 120 including DM18_27 ----
cat("\n\n===== FULL 120 (including DM18_27) =====\n")
seg_full <- seg.sites(aln_filt)
cat("Segregating sites:", length(seg_full), "\n")
h_full <- haplotype(aln_filt)
cat("Unique haplotypes:", nrow(h_full), "/", length(aln_filt), "\n")
cat("Hd:", round(hap.div(aln_filt), 4), "\n")
cat("Pi:", round(nuc.div(aln_filt), 5), "\n")

## ---- Nurse-only Tajima's D ----
## You'll need to define nurse sample names
## Update this list with the full set of nurses in the 120
nurse_samples <- c(
    ## Original 105 nurses
    "DD_11_04", "DD_11_22", "DD_8_06", "DD_19_25",
    "DD_21_04", "DD_21_05", "DD_21_06", "DD_21_08",
    "DL5_02", "DL5_04", "DL5_05", "DL5_06", "DL5_07", 
    "DL5_08", "DL5_09", "DL5_10", "DL5_23", "DL5_24",
    "DM18_22",
    ## EG nurses (excluding DM18_27 which is already excluded)
    "DD_8_2", "DD_11_3", "DD_11_6", "DD_11_23",
    "DM18_21", "DM18_23", "DM18_25", "DM18_26"
)

nurse_in_data <- nurse_samples[nurse_samples %in% labels(aln_119)]
cat("\nNurse samples found:", length(nurse_in_data), "\n")
if (length(nurse_in_data) >= 4) {
    aln_nurse <- aln_119[labels(aln_119) %in% nurse_in_data]
    tajd_nurse <- tajima.test(aln_nurse)
    cat("Nurse Tajima's D:", round(tajd_nurse$D, 3), 
        " p =", round(tajd_nurse$Pval.normal, 4), "\n")
}
