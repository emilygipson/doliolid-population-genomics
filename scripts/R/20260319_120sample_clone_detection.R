## 20260319_120sample_clone_detection.R
## Clone detection analysis on expanded 120-sample dataset
## Special focus on enriched within-tow groups from EG samples
## Requires: all_120_mitobim_filtered.fasta
## Run locally in RStudio

library(ape)
library(pegas)
library(poppr)

## ---- Read filtered alignment ----
aln <- read.FASTA("all_120_mitobim_filtered.fasta")
cat("Samples:", length(aln), "\n")

## ---- Remove DM18_27 and EG4 ----
exclude <- c("DM18_27", "DD_9_29")
keep <- !labels(aln) %in% exclude
aln <- aln[keep]
n <- length(aln)
cat("After exclusions:", n, "samples\n")

## ---- Pairwise distance matrix ----
d <- dist.dna(aln, model = "N", pairwise.deletion = TRUE)
d_mat <- as.matrix(d)
diag(d_mat) <- NA

cat("\nPairwise distance summary:\n")
cat("  Min:", min(d_mat, na.rm = TRUE), "\n")
cat("  Mean:", round(mean(d_mat, na.rm = TRUE), 2), "\n")
cat("  Median:", median(d_mat, na.rm = TRUE), "\n")
cat("  Max:", max(d_mat, na.rm = TRUE), "\n")

## ---- Find closest pairs ----
min_d <- min(d_mat, na.rm = TRUE)
closest <- which(d_mat == min_d, arr.ind = TRUE)
cat("\nClosest pairs (", min_d, "differences):\n")
seen <- c()
for (i in 1:nrow(closest)) {
    pair <- sort(c(closest[i, 1], closest[i, 2]))
    key <- paste(pair, collapse = "-")
    if (!key %in% seen) {
        seen <- c(seen, key)
        cat("  ", rownames(d_mat)[pair[1]], " vs ", rownames(d_mat)[pair[2]], "\n")
    }
}

## ---- Pairwise distance histogram ----
png("clone_threshold_histogram_120.png", width = 8, height = 5, units = "in", res = 300)
hist(as.numeric(d), breaks = 50, main = "Pairwise Nucleotide Differences (119 samples)",
     xlab = "Nucleotide differences", ylab = "Frequency", col = "steelblue", border = "white")
abline(v = min_d, col = "red", lwd = 2, lty = 2)
text(min_d + 5, par("usr")[4] * 0.9, paste("Min =", min_d), col = "red", adj = 0)
dev.off()
cat("\nHistogram saved: clone_threshold_histogram_120.png\n")

## ---- poppr clone analysis ----
## Convert to genind for poppr
## poppr needs loci, so we use each variable site as a "locus"
mat <- as.character(as.matrix(aln))
seg_idx <- seg.sites(aln)
mat_seg <- mat[, seg_idx]

## Convert to genind-compatible format
## Each segregating site = one locus with alleles a,c,g,t
## Build a data frame of alleles
geno_df <- as.data.frame(mat_seg, stringsAsFactors = FALSE)
colnames(geno_df) <- paste0("L", 1:ncol(geno_df))

## Create genclone object via df2genind
gi <- df2genind(geno_df, ploidy = 1, ind.names = labels(aln))
gc <- as.genclone(gi)

## MLG counts at various thresholds
cat("\nMLG counts at distance thresholds:\n")
cat("Threshold\tMLGs\n")
for (thresh in c(0, 1, 2, 3, 5, 7, 10, 15, 20)) {
    mlg_count <- mlg(gc, quiet = TRUE)
    ## Use filter to collapse at threshold
    if (thresh == 0) {
        cat(thresh, "\t", mlg_count, "\n")
    } else {
        dd <- diss.dist(gc)
        mlg.filter(gc) <- thresh
        filt_count <- mlg(gc, quiet = TRUE)
        mlg.filter(gc) <- 0  ## reset
        cat(thresh, "\t", filt_count, "\n")
    }
}

## ---- Within-tow analyses ----
## DD11 cast 3 (now 12 samples with EG6, EG7, EG8)
dd11_cast3 <- c(
    "DD_11_04", "DD_11_05", "DD_11_06", "DD_11_07", "DD_11_08",
    "DD_11_11", "DD_11_22", "DD_11_24", "DD_11_25",
    "DD_11_3", "DD_11_6", "DD_11_23"
)
dd11_in <- dd11_cast3[dd11_cast3 %in% labels(aln)]
cat("\n--- DD11 cast 3 ---\n")
cat("Samples found:", length(dd11_in), "\n")
if (length(dd11_in) >= 2) {
    d_dd11 <- as.matrix(dist.dna(aln[dd11_in], model = "N", pairwise.deletion = TRUE))
    diag(d_dd11) <- NA
    cat("Min pairwise diff:", min(d_dd11, na.rm = TRUE), "\n")
    cat("Mean pairwise diff:", round(mean(d_dd11, na.rm = TRUE), 2), "\n")
    cat("All unique:", length(dd11_in) == length(unique(apply(as.character(as.matrix(aln[dd11_in])), 1, paste, collapse = ""))), "\n")
}

## DM18 station 4 / 45m (now includes EG9, EG10)
## Original: DM18_1 through DM18_20 minus some, plus EG9=DM18_2, EG10=DM18_10
## We know original had 17 gono from stn4. EG9 and EG10 are also stn4.
## Since EG samples use their real IDs (DM18_2, DM18_10), they should already
## be in the alignment under those names.
## The original 17 from stn4/45m need to be defined from metadata.
## For now, identify all DM18 samples and let Emily subset by station.
dm18_samples <- labels(aln)[grepl("^DM18_", labels(aln))]
cat("\n--- All DM18 samples in dataset ---\n")
cat("Count:", length(dm18_samples), "\n")
cat("Names:", paste(dm18_samples, collapse = ", "), "\n")

## DM18 station 5 / 25m (DM18_22 + EG11-14 = DM18_21, DM18_23, DM18_25, DM18_26)
## DM18_27 already excluded
stn5_samples <- c("DM18_22", "DM18_21", "DM18_23", "DM18_25", "DM18_26")
stn5_in <- stn5_samples[stn5_samples %in% labels(aln)]
cat("\n--- DM18 station 5 / 25m ---\n")
cat("Samples found:", length(stn5_in), "\n")
if (length(stn5_in) >= 2) {
    d_stn5 <- as.matrix(dist.dna(aln[stn5_in], model = "N", pairwise.deletion = TRUE))
    diag(d_stn5) <- NA
    cat("Min pairwise diff:", min(d_stn5, na.rm = TRUE), "\n")
    cat("Mean pairwise diff:", round(mean(d_stn5, na.rm = TRUE), 2), "\n")
    n_unique <- length(unique(apply(as.character(as.matrix(aln[stn5_in])), 1, paste, collapse = "")))
    cat("Unique haplotypes:", n_unique, "/", length(stn5_in), "\n")
}

## ---- Summary ----
cat("\n\n===== CLONE DETECTION SUMMARY =====\n")
cat("Total samples:", n, "\n")
cat("Naive MLGs (unique haplotypes):", length(unique(apply(mat, 1, paste, collapse = ""))), "/", n, "\n")
cat("Min pairwise distance:", min(d_mat, na.rm = TRUE), "\n")
cat("No clonal pairs if min pairwise distance >> 0\n")
