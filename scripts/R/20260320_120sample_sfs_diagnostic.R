## 20260319_120sample_sfs_diagnostic.R
## Site frequency spectrum diagnostic analysis on expanded 120-sample dataset
## Progressive singleton/doubleton removal, genomic distribution, SFS shape
## Requires: all_120_mitobim_filtered.fasta
## Run locally in RStudio

library(ape)
library(pegas)

## ---- Read filtered alignment ----
aln <- read.FASTA("all_120_mitobim_filtered.fasta")

## Remove DM18_27 and DD_9_29
exclude <- c("DM18_27", "DD_9_29")
keep <- !labels(aln) %in% exclude
aln <- aln[keep]
n <- length(aln)
cat("Samples:", n, "\n")

mat <- as.character(as.matrix(aln))
cat("Filtered sites:", ncol(mat), "\n")

## ---- Compute folded SFS ----
seg_idx <- seg.sites(aln)
cat("Segregating sites:", length(seg_idx), "\n")

## Count minor allele frequency at each segregating site
seg_mat <- mat[, seg_idx]
sfs <- numeric(floor(n / 2))

for (j in 1:ncol(seg_mat)) {
    col <- seg_mat[, j]
    tab <- table(col)
    ## Minor allele count = minimum count among alleles present
    min_count <- min(tab)
    ## Fold: if min_count > n/2, use n - min_count
    folded_count <- min(min_count, n - min_count)
    if (folded_count >= 1 && folded_count <= length(sfs)) {
        sfs[folded_count] <- sfs[folded_count] + 1
    }
}

cat("\nFolded SFS (first 10 classes):\n")
cat("Freq\tCount\tPct\n")
total_snps <- sum(sfs)
for (i in 1:min(10, length(sfs))) {
    cat(i, "\t", sfs[i], "\t", round(100 * sfs[i] / total_snps, 1), "%\n", sep = "")
}

## ---- Expected SFS under constant-size neutral model ----
## E[xi] = theta / i for unfolded; folded is sum of complementary classes
expected_unfolded <- 1 / (1:(n - 1))
expected_unfolded <- expected_unfolded / sum(expected_unfolded)  ## normalize

## Fold the expected SFS
expected_folded <- numeric(floor(n / 2))
for (i in 1:length(expected_folded)) {
    if (i < n - i) {
        expected_folded[i] <- expected_unfolded[i] + expected_unfolded[n - i]
    } else if (i == n - i) {
        expected_folded[i] <- expected_unfolded[i]
    }
}
expected_folded <- expected_folded / sum(expected_folded)  ## normalize

## Compare observed vs expected
cat("\nObserved vs Expected (neutral constant-size):\n")
cat("Freq\tObserved\tObs_pct\tExpected_pct\tRatio\n")
for (i in 1:min(5, length(sfs))) {
    obs_pct <- 100 * sfs[i] / total_snps
    exp_pct <- 100 * expected_folded[i]
    ratio <- obs_pct / exp_pct
    cat(i, "\t", sfs[i], "\t", round(obs_pct, 1), "%\t", 
        round(exp_pct, 1), "%\t", round(ratio, 2), "x\n", sep = "")
}

## ---- Progressive removal analysis ----
## Remove singletons, then doubletons, etc., and recompute Tajima's D each time
cat("\n--- Progressive removal of rare variant classes ---\n")
cat("Removed\tRemaining_SNPs\tTajima_D\tp\n")

## Identify frequency class of each segregating site
site_freq_class <- integer(ncol(seg_mat))
for (j in 1:ncol(seg_mat)) {
    col <- seg_mat[, j]
    tab <- table(col)
    min_count <- min(tab)
    site_freq_class[j] <- min(min_count, n - min_count)
}

## Progressive removal
for (remove_up_to in 0:4) {
    if (remove_up_to == 0) {
        keep_sites <- rep(TRUE, ncol(mat))
    } else {
        ## Keep only segregating sites with freq class > remove_up_to
        ## plus all invariant sites
        is_seg <- rep(FALSE, ncol(mat))
        is_seg[seg_idx] <- TRUE
        
        ## For segregating sites, check if their freq class exceeds threshold
        seg_keep <- site_freq_class > remove_up_to
        
        keep_sites <- rep(TRUE, ncol(mat))
        seg_counter <- 0
        for (col_i in 1:ncol(mat)) {
            if (is_seg[col_i]) {
                seg_counter <- seg_counter + 1
                if (!seg_keep[seg_counter]) {
                    keep_sites[col_i] <- FALSE
                }
            }
        }
    }
    
    mat_sub <- mat[, keep_sites]
    aln_sub <- as.DNAbin(mat_sub)
    
    remaining <- length(seg.sites(aln_sub))
    
    if (remaining > 0) {
        tajd <- tajima.test(aln_sub)
        label <- if (remove_up_to == 0) "None" else paste0("freq<=", remove_up_to)
        cat(label, "\t", remaining, "\t", round(tajd$D, 3), "\t", 
            round(tajd$Pval.normal, 4), "\n", sep = "")
    }
}

## ---- Singleton genomic distribution ----
cat("\n--- Singleton genomic distribution ---\n")

singleton_positions <- seg_idx[site_freq_class == 1]
n_singletons <- length(singleton_positions)
cat("Total singletons:", n_singletons, "\n")

## 500-bp windows
window_size <- 500
max_pos <- ncol(mat)
n_windows <- ceiling(max_pos / window_size)
window_counts <- integer(n_windows)

for (pos in singleton_positions) {
    w <- ceiling(pos / window_size)
    if (w <= n_windows) window_counts[w] <- window_counts[w] + 1
}

cat("Mean singletons per 500bp window:", round(mean(window_counts), 1), "\n")
cat("SD:", round(sd(window_counts), 1), "\n")
cat("Max:", max(window_counts), "\n")
threshold <- mean(window_counts) + 3 * sd(window_counts)
cat("Hotspot threshold (mean + 3SD):", round(threshold, 1), "\n")
cat("Windows exceeding threshold:", sum(window_counts > threshold), "\n")

## ---- Per-sample singleton counts ----
cat("\n--- Per-sample singleton counts ---\n")
sample_singletons <- integer(n)
names(sample_singletons) <- rownames(mat)

for (j in which(site_freq_class == 1)) {
    col <- seg_mat[, j]
    tab <- table(col)
    minor_allele <- names(tab)[which.min(tab)]
    carriers <- which(seg_mat[, j] == minor_allele)
    for (c in carriers) {
        sample_singletons[c] <- sample_singletons[c] + 1
    }
}

cat("Mean singletons per sample:", round(mean(sample_singletons), 1), "\n")
cat("Max:", max(sample_singletons), "(", names(which.max(sample_singletons)), ")\n")
cat("Min:", min(sample_singletons), "\n")
cat("Samples with 0 singletons:", sum(sample_singletons == 0), "\n")

## ---- Figures ----
## SFS barplot
png("sfs_120_observed_vs_expected.png", width = 8, height = 5, units = "in", res = 300)
max_show <- min(15, length(sfs))
barplot_data <- rbind(
    sfs[1:max_show] / total_snps * 100,
    expected_folded[1:max_show] * 100
)
bp <- barplot(barplot_data, beside = TRUE, 
              names.arg = 1:max_show,
              col = c("steelblue", "gray70"),
              main = paste0("Folded SFS (n=", n, " samples)"),
              xlab = "Minor allele count",
              ylab = "Percentage of SNPs",
              legend.text = c("Observed", "Expected (neutral)"),
              args.legend = list(x = "topright"))
dev.off()
cat("\nFigure saved: sfs_120_observed_vs_expected.png\n")

## Singleton distribution across genome
png("singleton_distribution_120.png", width = 10, height = 4, units = "in", res = 300)
plot(1:n_windows, window_counts, type = "h", col = "steelblue",
     xlab = "Genomic window (500 bp)", ylab = "Singleton count",
     main = paste0("Singleton Distribution Across Mitogenome (n=", n, ")"))
abline(h = threshold, col = "red", lty = 2)
text(n_windows * 0.8, threshold + 1, paste("Threshold:", round(threshold, 1)), col = "red")
dev.off()
cat("Figure saved: singleton_distribution_120.png\n")

## ---- Tab-delimited SFS table ----
cat("\n\n===== TAB-DELIMITED SFS TABLE =====\n")
cat("Minor_allele_count\tObserved_count\tObserved_pct\tExpected_pct\tExcess_ratio\n")
for (i in 1:min(10, length(sfs))) {
    obs_pct <- 100 * sfs[i] / total_snps
    exp_pct <- 100 * expected_folded[i]
    ratio <- obs_pct / exp_pct
    cat(i, "\t", sfs[i], "\t", round(obs_pct, 1), "\t", round(exp_pct, 1), "\t", 
        round(ratio, 2), "\n", sep = "")
}
