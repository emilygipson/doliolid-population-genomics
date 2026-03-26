# =============================================================================
# Per-gene dN/dS analysis — D. gegenbauri mitogenomes
# 59 complete GetOrganelle assemblies, 13 protein-coding genes
# Emily Garner, UGA, 2026-03-06
# =============================================================================

library(seqinr)
library(ape)

genes <- c("cox1", "cox2", "cox3", "cob", "nad1", "nad2", "nad3",
           "nad4", "nad4L", "nad5", "nad6", "atp6", "atp8")

results <- data.frame(
  gene = character(),
  n_seqs = integer(),
  length_bp = integer(),
  length_codons = integer(),
  pi = numeric(),
  dN = numeric(),
  dS = numeric(),
  omega = numeric(),
  stringsAsFactors = FALSE
)

for (gene in genes) {
  cat("Processing", gene, "...\n")
  
  fname <- paste0(gene, "_aligned.fasta")
  if (!file.exists(fname)) {
    cat("  WARNING: file not found, skipping\n")
    next
  }
  
  aln <- read.alignment(fname, format = "fasta")
  dna <- read.dna(fname, format = "fasta")
  
  n_seqs <- aln$nb
  len_bp <- nchar(aln$seq[[1]])
  len_codons <- len_bp %/% 3
  pi_val <- nuc.div(dna)
  
  kaks_result <- kaks(aln)
  
  dn_vals <- as.numeric(kaks_result$ka)
  ds_vals <- as.numeric(kaks_result$ks)
  
  dn_vals <- dn_vals[!is.na(dn_vals) & is.finite(dn_vals)]
  ds_vals <- ds_vals[!is.na(ds_vals) & is.finite(ds_vals)]
  
  mean_dn <- mean(dn_vals)
  mean_ds <- mean(ds_vals)
  omega <- ifelse(mean_ds > 0, mean_dn / mean_ds, NA)
  
  results <- rbind(results, data.frame(
    gene = gene, n_seqs = n_seqs, length_bp = len_bp,
    length_codons = len_codons, pi = pi_val,
    dN = mean_dn, dS = mean_ds, omega = omega,
    stringsAsFactors = FALSE
  ))
  
  cat("  ", n_seqs, "seqs,", len_codons, "codons, dN =", round(mean_dn, 5),
      ", dS =", round(mean_ds, 5), ", omega =", round(omega, 4), "\n")
}

cat("\n===== Per-gene diversity and selection =====\n\n")
print(results, digits = 4)

# =============================================================================
# Figure
# =============================================================================

results$gene <- factor(results$gene, levels = results$gene[order(results$omega)])

pdf("per_gene_dNdS.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5, 6, 3, 1))

barplot(results$omega[order(results$omega)],
        names.arg = results$gene[order(results$omega)],
        horiz = TRUE, las = 1,
        xlab = expression(omega ~ "(dN/dS)"),
        main = "Selection on mitochondrial genes",
        col = "#2e6da4",
        border = NA, xlim = c(0, max(results$omega, na.rm = TRUE) * 1.1))
abline(v = 1, lty = 2, col = "gray40")

barplot(results$pi[order(results$omega)],
        names.arg = results$gene[order(results$omega)],
        horiz = TRUE, las = 1,
        xlab = expression(pi ~ "(nucleotide diversity)"),
        main = "Nucleotide diversity per gene",
        col = "#2e6da4",
        border = NA)
dev.off()

# =============================================================================
# Summary statistics for text
# =============================================================================

cat("\nSummary for paper:\n")
cat("  omega range:", round(min(results$omega, na.rm=TRUE), 3), "-",
    round(max(results$omega, na.rm=TRUE), 3), "\n")
cat("  Most constrained:", results$gene[which.min(results$omega)],
    "(omega =", round(min(results$omega, na.rm=TRUE), 3), ")\n")
cat("  Least constrained:", results$gene[which.max(results$omega)],
    "(omega =", round(max(results$omega, na.rm=TRUE), 3), ")\n")
cat("  Mean omega across all genes:", round(mean(results$omega, na.rm=TRUE), 3), "\n")
cat("  pi range:", round(min(results$pi), 4), "-", round(max(results$pi), 4), "\n")
