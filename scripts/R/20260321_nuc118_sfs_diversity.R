## 20260321_nuc118_sfs_diversity.R
## SFS plots and diversity stats for 118-sample nuclear dataset
## Compares global, bloom, and non-bloom SFS and diversity
## Updates mito-nuclear pi ratio
## Run locally in RStudio
##
## INPUT FILES (all in working directory):
##   nuc118_global.sfs                       -- folded SFS (global, 118 samples)
##   bloom.sfs                               -- folded SFS (bloom group)
##   nonbloom.sfs                            -- folded SFS (non-bloom group)
##   nuc118_thetas.thetas.idx.pestPG         -- per-contig diversity (global)
##   bloom_thetas.thetas.idx.pestPG          -- per-contig diversity (bloom)
##   nonbloom_thetas.thetas.idx.pestPG       -- per-contig diversity (non-bloom)

## ============================================================
## PART 1: DIVERSITY STATS FROM THETAS FILES
## ============================================================

## ---- Read per-contig diversity and compute genome-wide summaries ----
read_pestPG <- function(file) {
    df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
    return(df)
}

global_thetas <- read_pestPG("nuc118_thetas.thetas.idx.pestPG")
bloom_thetas <- read_pestPG("bloom_thetas.thetas.idx.pestPG")
nonbloom_thetas <- read_pestPG("nonbloom_thetas.thetas.idx.pestPG")

## Genome-wide diversity: sum tW and tP, divide by total nSites
summarize_diversity <- function(df, label) {
    total_tW <- sum(df$tW, na.rm = TRUE)
    total_tP <- sum(df$tP, na.rm = TRUE)
    total_sites <- sum(df$nSites, na.rm = TRUE)
    n_contigs <- nrow(df)
    pi <- total_tP / total_sites
    thetaW <- total_tW / total_sites
    ## Genome-wide Tajima's D is not simply the mean of per-contig D values.
    ## It must be computed from the genome-wide tW and tP.
    ## D = (tP - tW) / sqrt(Var), but we don't have the variance term from pestPG.
    ## Report the per-site values and note that D should be computed from the SFS.
    cat(sprintf("\n--- %s ---\n", label))
    cat(sprintf("  Contigs: %d\n", n_contigs))
    cat(sprintf("  Total sites: %s\n", format(total_sites, big.mark = ",")))
    cat(sprintf("  Theta_W (total): %.2f\n", total_tW))
    cat(sprintf("  Theta_Pi (total): %.2f\n", total_tP))
    cat(sprintf("  Pi (per-site): %.6f\n", pi))
    cat(sprintf("  Theta_W (per-site): %.6f\n", thetaW))
    return(data.frame(label = label, n_contigs = n_contigs,
                      total_sites = total_sites,
                      tW = total_tW, tP = total_tP,
                      pi = pi, thetaW = thetaW))
}

cat("===== NUCLEAR DIVERSITY STATS (118 samples) =====\n")
sum_global <- summarize_diversity(global_thetas, "Global (118)")
sum_bloom <- summarize_diversity(bloom_thetas, "Bloom (80)")
sum_nonbloom <- summarize_diversity(nonbloom_thetas, "Non-bloom (38)")

## ---- Compare to previous 100-sample results ----
cat("\n--- Comparison to Previous 100-sample Results ---\n")
cat(sprintf("  Previous pi: 0.0041\n"))
cat(sprintf("  Current pi:  %.4f\n", sum_global$pi))
cat(sprintf("  Change: %+.4f (%.1f%%)\n",
            sum_global$pi - 0.0041,
            (sum_global$pi - 0.0041) / 0.0041 * 100))
cat(sprintf("  Previous sites: 5,608,457\n"))
cat(sprintf("  Current sites:  %s\n", format(sum_global$total_sites, big.mark = ",")))

## ---- Mito-nuclear pi ratio ----
## Mito pi from 119-sample dataset (excluding DM18_27): 0.00641
mito_pi <- 0.00641
nuc_pi <- sum_global$pi
ratio <- mito_pi / nuc_pi
cat(sprintf("\n--- Mito-Nuclear Pi Ratio ---\n"))
cat(sprintf("  Mito pi: %.5f (119 samples, 11,016 bp)\n", mito_pi))
cat(sprintf("  Nuclear pi: %.6f (118 samples, %s sites)\n", nuc_pi,
            format(sum_global$total_sites, big.mark = ",")))
cat(sprintf("  Ratio (mito/nuclear): %.2f\n", ratio))
cat(sprintf("  Expected for invertebrates: ~5-10x (highly variable)\n"))

## ---- Bloom vs non-bloom comparison ----
cat(sprintf("\n--- Bloom vs Non-bloom Nuclear Diversity ---\n"))
cat(sprintf("  Bloom pi:     %.6f\n", sum_bloom$pi))
cat(sprintf("  Non-bloom pi: %.6f\n", sum_nonbloom$pi))
cat(sprintf("  Ratio (bloom/non-bloom): %.3f\n", sum_bloom$pi / sum_nonbloom$pi))
cat(sprintf("  Fst (weighted): 0.002\n"))
cat(sprintf("  Fst (unweighted): 0.007\n"))

## ---- Write diversity summary table ----
div_summary <- rbind(sum_global, sum_bloom, sum_nonbloom)
div_summary$pi_rounded <- round(div_summary$pi, 6)
write.table(div_summary, "20260321_nuc118_diversity_summary.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nSaved: 20260321_nuc118_diversity_summary.tsv\n")

## ============================================================
## PART 2: SFS PLOTS
## ============================================================

## ---- Read SFS files ----
## ANGSD folded SFS: space-delimited, n+1 entries for n diploid samples (2n+1 entries)
## First entry is the count of monomorphic sites
read_sfs <- function(file) {
    sfs <- scan(file, what = numeric(), quiet = TRUE)
    return(sfs)
}

sfs_global <- read_sfs("nuc118_global.sfs")
sfs_bloom <- read_sfs("bloom.sfs")
sfs_nonbloom <- read_sfs("nonbloom.sfs")

cat(sprintf("\nGlobal SFS entries: %d (expected for 118 diploid: %d)\n",
            length(sfs_global), 2 * 118 + 1))
cat(sprintf("Bloom SFS entries: %d\n", length(sfs_bloom)))
cat(sprintf("Non-bloom SFS entries: %d\n", length(sfs_nonbloom)))

## ---- Plot 1: Global folded SFS (skip monomorphic bin) ----
## For a folded SFS of n diploid individuals, entries go from 0 to n
## We skip bin 0 (monomorphic) and plot bins 1 to n
n_global <- (length(sfs_global) - 1) / 2  ## number of diploid individuals
sfs_folded_global <- sfs_global[2:(n_global + 1)]

png("20260321_sfs_global.png", width = 8, height = 5, units = "in", res = 300)
par(mar = c(5, 5, 3, 2))
barplot(sfs_folded_global,
        names.arg = 1:length(sfs_folded_global),
        xlab = "Minor allele count",
        ylab = "Number of sites",
        main = sprintf("Nuclear Folded SFS (118 samples, %s total sites)",
                       format(round(sum(sfs_global)), big.mark = ",")),
        col = "steelblue", border = NA,
        las = 1)
dev.off()
cat("Saved: 20260321_sfs_global.png\n")

## ---- Plot 2: Zoom on low-frequency end (bins 1-20) ----
png("20260321_sfs_global_zoom.png", width = 8, height = 5, units = "in", res = 300)
par(mar = c(5, 5, 3, 2))
n_show <- min(20, length(sfs_folded_global))
barplot(sfs_folded_global[1:n_show],
        names.arg = 1:n_show,
        xlab = "Minor allele count",
        ylab = "Number of sites",
        main = "Nuclear Folded SFS - Low Frequency End",
        col = "steelblue", border = NA,
        las = 1)
dev.off()
cat("Saved: 20260321_sfs_global_zoom.png\n")

## ---- Plot 3: Bloom vs non-bloom SFS comparison (normalized) ----
n_bloom <- (length(sfs_bloom) - 1) / 2
n_nonbloom <- (length(sfs_nonbloom) - 1) / 2
sfs_folded_bloom <- sfs_bloom[2:(n_bloom + 1)]
sfs_folded_nonbloom <- sfs_nonbloom[2:(n_nonbloom + 1)]

## Normalize to proportions for comparison (different sample sizes)
prop_bloom <- sfs_folded_bloom / sum(sfs_folded_bloom)
prop_nonbloom <- sfs_folded_nonbloom / sum(sfs_folded_nonbloom)

## Show first 20 bins of each
n_compare <- 20

png("20260321_sfs_bloom_comparison.png", width = 8, height = 5, units = "in", res = 300)
par(mar = c(5, 5, 3, 2))
plot(1:min(n_compare, length(prop_bloom)),
     prop_bloom[1:min(n_compare, length(prop_bloom))],
     type = "b", pch = 19, col = "#E41A1C", lwd = 2,
     xlab = "Minor allele count",
     ylab = "Proportion of segregating sites",
     main = "SFS Shape Comparison: Bloom vs Non-bloom",
     ylim = c(0, max(c(prop_bloom[1:min(n_compare, length(prop_bloom))],
                       prop_nonbloom[1:min(n_compare, length(prop_nonbloom))])) * 1.1))
lines(1:min(n_compare, length(prop_nonbloom)),
      prop_nonbloom[1:min(n_compare, length(prop_nonbloom))],
      type = "b", pch = 17, col = "#377EB8", lwd = 2)
legend("topright",
       legend = c(sprintf("Bloom (n=%d)", n_bloom),
                  sprintf("Non-bloom (n=%d)", n_nonbloom)),
       col = c("#E41A1C", "#377EB8"),
       pch = c(19, 17), lwd = 2, bty = "n")
dev.off()
cat("Saved: 20260321_sfs_bloom_comparison.png\n")

## ---- Summary table for paper ----
cat("\n\n===== TAB-DELIMITED SUMMARY FOR PAPER =====\n")
cat("Metric\tGlobal\tBloom\tNon-bloom\n")
cat(sprintf("N samples\t118\t80\t38\n"))
cat(sprintf("Sites analyzed\t%s\t%s\t%s\n",
            format(sum_global$total_sites, big.mark = ","),
            format(sum_bloom$total_sites, big.mark = ","),
            format(sum_nonbloom$total_sites, big.mark = ",")))
cat(sprintf("Pi\t%.6f\t%.6f\t%.6f\n", sum_global$pi, sum_bloom$pi, sum_nonbloom$pi))
cat(sprintf("Theta_W (per-site)\t%.6f\t%.6f\t%.6f\n",
            sum_global$thetaW, sum_bloom$thetaW, sum_nonbloom$thetaW))
cat(sprintf("Fst (weighted)\t0.002\t-\t-\n"))
cat(sprintf("Mito pi\t0.00641\t-\t-\n"))
cat(sprintf("Mito/nuclear pi ratio\t%.2f\t-\t-\n", ratio))
