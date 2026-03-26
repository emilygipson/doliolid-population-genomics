## 20260319_120sample_pergene_diversity.R
## Per-gene diversity and neutrality tests on expanded 120-sample dataset
## Uses unfiltered alignment + gene coordinates from DD_21_05 GFF
## Requires: all_120_mitobim_aligned.fasta (unfiltered)
## Run locally in RStudio

library(ape)
library(pegas)

## ---- Read UNFILTERED alignment ----
aln <- read.FASTA("all_120_mitobim_aligned.fasta")
mat <- as.character(as.matrix(aln))
cat("Raw alignment:", nrow(mat), "x", ncol(mat), "\n")

## ---- Remove DM18_27 and DD_9_29 ----
exclude <- c("DM18_27", "DD_9_29")
keep <- !rownames(mat) %in% exclude
mat <- mat[keep, ]
cat("After exclusions:", nrow(mat), "samples\n")

## ---- Gene coordinates (from DD_21_05 GFF, 1-based) ----
## These are coordinates on the UNROTATED DD_21_05 reference (15,329 bp)
## The MITObim alignment uses DD_21_05 as reference, so the alignment
## coordinates should roughly correspond (with gap offsets)
genes <- data.frame(
    gene = c("nad4", "nad6", "cox3", "cox2", "cob", "nad4l", "nad3",
             "nad5", "atp6", "cox1", "atp8", "nad1", "nad2"),
    start = c(230, 3106, 4451, 5369, 6018, 7342, 7587,
              8455, 10218, 10930, 12593, 12898, 13931),
    end = c(1558, 3591, 5245, 6034, 7097, 7602, 7934,
            10140, 10847, 12480, 12739, 13806, 14938),
    stringsAsFactors = FALSE
)
genes$length <- genes$end - genes$start + 1

## ---- Map reference coordinates to alignment columns ----
## Find the DD_21_05 row (should be in the alignment as it was the reference)
ref_name <- "DD_21_05"
if (!ref_name %in% rownames(mat)) {
    cat("WARNING: DD_21_05 not found in alignment. Available names:\n")
    print(head(rownames(mat), 20))
    cat("Trying to find a close match...\n")
    ref_idx <- grep("DD_21_05", rownames(mat))
    if (length(ref_idx) > 0) {
        ref_name <- rownames(mat)[ref_idx[1]]
        cat("Using:", ref_name, "\n")
    } else {
        stop("Cannot find reference sample in alignment")
    }
}

ref_seq <- mat[ref_name, ]

## Build mapping: reference position (ungapped) -> alignment column
ref_pos <- 0
ref_to_aln <- integer(0)
for (col in 1:length(ref_seq)) {
    if (ref_seq[col] %in% c("a", "c", "g", "t")) {
        ref_pos <- ref_pos + 1
        ref_to_aln[ref_pos] <- col
    }
}
cat("Reference ungapped length:", ref_pos, "\n")

## ---- Per-gene analysis ----
results <- data.frame(
    Gene = character(), Sites = integer(), Clean_sites = integer(),
    Pct_kept = numeric(), Haps = integer(), Hd = numeric(),
    Pi = numeric(), TajD = numeric(), TajD_p = numeric(),
    stringsAsFactors = FALSE
)

for (i in 1:nrow(genes)) {
    g <- genes$gene[i]
    gstart <- genes$start[i]
    gend <- genes$end[i]
    
    ## Map to alignment columns
    if (gstart > length(ref_to_aln) || gend > length(ref_to_aln)) {
        cat("WARNING:", g, "coordinates exceed reference length. Skipping.\n")
        next
    }
    
    aln_start <- ref_to_aln[gstart]
    aln_end <- ref_to_aln[gend]
    
    if (is.na(aln_start) || is.na(aln_end)) {
        cat("WARNING:", g, "could not map coordinates. Skipping.\n")
        next
    }
    
    ## Extract gene region
    gene_mat <- mat[, aln_start:aln_end]
    raw_cols <- ncol(gene_mat)
    
    ## Filter: keep only columns where all samples have ACGT
    valid <- c("a", "c", "g", "t")
    keep_cols <- apply(gene_mat, 2, function(col) all(col %in% valid))
    gene_mat_clean <- gene_mat[, keep_cols, drop = FALSE]
    clean_cols <- ncol(gene_mat_clean)
    pct_kept <- round(100 * clean_cols / raw_cols, 1)
    
    if (clean_cols < 10) {
        cat(g, ": Only", clean_cols, "clean sites. Skipping.\n")
        next
    }
    
    ## Convert to DNAbin
    gene_dna <- as.DNAbin(gene_mat_clean)
    
    ## Stats
    nhaps <- nrow(haplotype(gene_dna))
    hd <- hap.div(gene_dna)
    pi_val <- nuc.div(gene_dna)
    tajd <- tajima.test(gene_dna)
    
    results <- rbind(results, data.frame(
        Gene = g, Sites = clean_cols, Raw_sites = raw_cols,
        Pct_kept = pct_kept, Haps = nhaps, 
        Hd = round(hd, 3), Pi = round(pi_val, 5),
        TajD = round(tajd$D, 3), TajD_p = round(tajd$Pval.normal, 3),
        stringsAsFactors = FALSE
    ))
    
    cat(sprintf("%-6s %4d/%4d (%5.1f%%) %3d haps  Hd=%.3f  Pi=%.5f  D=%.3f (p=%.3f)\n",
                g, clean_cols, raw_cols, pct_kept, nhaps, hd, pi_val, tajd$D, tajd$Pval.normal))
}

## ---- Summary ----
cat("\nTotal coding sites (clean):", sum(results$Sites), "\n")
cat("All genes negative Tajima's D:", all(results$TajD < 0), "\n")
cat("Genes with p < 0.05:", sum(results$TajD_p < 0.05), "/", nrow(results), "\n")

## ---- Tab-delimited output for Excel ----
cat("\n\n===== TAB-DELIMITED TABLE =====\n")
cat("Gene\tClean_sites\tRaw_sites\tPct_kept\tHaplotypes\tHd\tPi\tTajima_D\tp\n")
for (i in 1:nrow(results)) {
    r <- results[i, ]
    cat(r$Gene, "\t", r$Sites, "\t", r$Raw_sites, "\t", r$Pct_kept, "\t",
        r$Haps, "\t", r$Hd, "\t", r$Pi, "\t", r$TajD, "\t", r$TajD_p, "\n", sep = "")
}
