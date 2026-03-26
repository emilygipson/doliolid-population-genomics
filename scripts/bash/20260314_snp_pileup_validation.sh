#!/bin/bash
#SBATCH --job-name=snp_pileup_val
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eeg37520@uga.edu
#SBATCH --output=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/snp_pileup_val_20260314.out
#SBATCH --error=/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/logs/snp_pileup_val_20260314.err

module load R/4.4.1-gfbf-2023b

Rscript - <<'EOF'
.libPaths(c("/scratch/eeg37520/Rlibs", .libPaths()))
library(ape)

# --- Paths ---
aln_unfiltered_path <- "/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/analysis/all_105_conspecific_aligned.fasta"
aln_filtered_path   <- "/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/analysis/all_105_conspecific_filtered.fasta"
pileup_dir          <- "/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/remap_pileup/pileup"
outdir              <- "/scratch/eeg37520/doliolid_popgen/mitobim_redo_conspecific/snp_validation"
dir.create(outdir, showWarnings = FALSE)

# --- Load alignments ---
cat("Loading alignments...\n")
aln_uf  <- as.character(as.matrix(read.FASTA(aln_unfiltered_path)))
aln_f   <- as.character(as.matrix(read.FASTA(aln_filtered_path)))

# Remove DM18_27
aln_uf <- aln_uf[rownames(aln_uf) != "DM18_27", ]
aln_f  <- aln_f[rownames(aln_f)   != "DM18_27", ]

n         <- nrow(aln_f)
nsites_f  <- ncol(aln_f)
nsites_uf <- ncol(aln_uf)
samples   <- rownames(aln_f)

cat("Samples:", n, "\n")
cat("Filtered sites:", nsites_f, "\n")
cat("Unfiltered sites:", nsites_uf, "\n")

# --- Map filtered columns to unfiltered columns ---
cat("Mapping filtered columns to unfiltered columns...\n")

valid_bases <- c("a", "c", "g", "t")

filtered_to_unfiltered <- integer(nsites_f)
fi <- 1L
for (ui in seq_len(nsites_uf)) {
  if (fi > nsites_f) break
  col <- tolower(aln_uf[, ui])
  if (all(col %in% valid_bases)) {
    filtered_to_unfiltered[fi] <- ui
    fi <- fi + 1L
  }
}

cat("Mapped", sum(filtered_to_unfiltered > 0), "filtered columns to unfiltered columns\n")

# --- Find all segregating sites and variant carriers ---
cat("Identifying segregating sites and variant carriers...\n")

snp_records <- list()

for (fi in seq_len(nsites_f)) {
  col <- tolower(aln_f[, fi])
  if (!all(col %in% valid_bases)) next
  tbl <- table(col)
  if (length(tbl) < 2) next

  majority <- names(sort(tbl, decreasing = TRUE))[1]
  for (allele in names(tbl)) {
    if (allele == majority) next
    carriers <- samples[col == allele]
    for (carrier in carriers) {
      snp_records[[length(snp_records) + 1]] <- list(
        filtered_col   = fi,
        unfiltered_col = filtered_to_unfiltered[fi],
        sample         = carrier,
        ref_allele     = majority,
        var_allele     = allele,
        allele_count   = tbl[[allele]]
      )
    }
  }
}

cat("Total variant carrier instances:", length(snp_records), "\n")

# --- Pre-compute cumulative non-gap counts ---
cat("Computing per-sample reference positions...\n")

cum_nongap <- matrix(0L, nrow = n, ncol = nsites_uf)
rownames(cum_nongap) <- samples

for (si in seq_len(n)) {
  seq_row <- tolower(aln_uf[si, ])
  is_base <- seq_row %in% valid_bases
  cum_nongap[si, ] <- cumsum(is_base)
}

cat("Reference position mapping complete\n")

# --- Load all pileup allele_freq files ---
cat("Loading pileup data...\n")

pileup_data <- list()
for (s in samples) {
  pf <- file.path(pileup_dir, paste0(s, ".allele_freq.txt"))
  if (file.exists(pf)) {
    pd <- read.table(pf, header = FALSE, sep = "\t",
                     col.names = c("contig", "pos", "ref_base", "depth",
                                   "n_match", "cons_freq"),
                     stringsAsFactors = FALSE)
    pileup_data[[s]] <- pd
  } else {
    cat("WARNING: pileup file not found for", s, "\n")
  }
}

cat("Loaded pileup data for", length(pileup_data), "samples\n")

# --- Validate each SNP ---
cat("Validating SNPs against pileup data...\n")

results <- data.frame(
  filtered_col     = integer(),
  unfiltered_col   = integer(),
  sample           = character(),
  ref_allele       = character(),
  var_allele       = character(),
  allele_count     = integer(),
  ref_pos          = integer(),
  pileup_depth     = integer(),
  pileup_cons_freq = numeric(),
  var_freq         = numeric(),
  status           = character(),
  stringsAsFactors = FALSE
)

for (rec in snp_records) {
  s  <- rec$sample
  ui <- rec$unfiltered_col
  fi <- rec$filtered_col

  if (ui == 0) next
  if (!s %in% names(pileup_data)) next

  ref_pos <- cum_nongap[s, ui]
  if (ref_pos == 0) next

  pd  <- pileup_data[[s]]
  row <- pd[pd$pos == ref_pos, ]

  if (nrow(row) == 0) {
    status           <- "MISSING_FROM_PILEUP"
    pileup_depth     <- NA
    pileup_cons_freq <- NA
    var_freq         <- NA
  } else {
    pileup_depth     <- row$depth[1]
    pileup_cons_freq <- row$cons_freq[1]
    var_freq         <- 1 - pileup_cons_freq
    if (pileup_depth < 5) {
      status <- "LOW_DEPTH"
    } else if (var_freq >= 0.10) {
      status <- "VALIDATED"
    } else if (var_freq >= 0.05) {
      status <- "WEAK_SUPPORT"
    } else {
      status <- "NOT_SUPPORTED"
    }
  }

  results <- rbind(results, data.frame(
    filtered_col     = fi,
    unfiltered_col   = ui,
    sample           = s,
    ref_allele       = rec$ref_allele,
    var_allele       = rec$var_allele,
    allele_count     = as.integer(rec$allele_count),
    ref_pos          = ref_pos,
    pileup_depth     = pileup_depth,
    pileup_cons_freq = pileup_cons_freq,
    var_freq         = var_freq,
    status           = status,
    stringsAsFactors = FALSE
  ))
}

# --- Summary ---
cat("\n=== SNP Pileup Validation Summary ===\n")
cat("Total variant instances evaluated:", nrow(results), "\n")
print(table(results$status))

cat("\nSingleton-specific results:\n")
singletons_only <- results[results$allele_count == 1, ]
cat("Singleton instances:", nrow(singletons_only), "\n")
print(table(singletons_only$status))

cat("\nNot-supported singletons (var_freq < 5%):\n")
bad <- singletons_only[!is.na(singletons_only$status) &
                         singletons_only$status == "NOT_SUPPORTED", ]
if (nrow(bad) > 0) {
  print(bad[, c("sample", "filtered_col", "ref_pos", "pileup_depth",
                "pileup_cons_freq", "var_freq")])
} else {
  cat("  None\n")
}

# --- Write output ---
write.table(results,
            file.path(outdir, "snp_pileup_validation_20260314.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(singletons_only,
            file.path(outdir, "singleton_validation_20260314.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nOutput written to:", outdir, "\n")
cat("Done.\n")
EOF
