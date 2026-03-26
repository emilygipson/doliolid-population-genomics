## 20260319_120sample_bloom_comparison.R
## Bloom vs non-bloom population structure analysis on expanded 120-sample dataset
## Requires: all_120_mitobim_filtered.fasta (from diversity stats script)
## Run locally in RStudio

library(ape)
library(pegas)

## ---- Read filtered alignment ----
aln <- read.FASTA("all_120_mitobim_filtered.fasta")
cat("Samples:", length(aln), "\n")

## ---- Remove DM18_27 (divergent lineage) and EG4 if present ----
exclude <- c("DM18_27", "DD_9_29")
keep <- !labels(aln) %in% exclude
aln <- aln[keep]
cat("After exclusions:", length(aln), "samples\n")

## ---- Define bloom status ----
## Bloom cruises: DD7, DD11, DD17, DD21, DM18
## Non-bloom cruises: DD4, DD6, DD8, DD9, DD10, DD15, DD18, DD19, DL5
bloom_cruises <- c("DD7", "DD11", "DD17", "DD21", "DM18")

get_cruise <- function(name) {
    if (grepl("^DM18_", name)) return("DM18")
    if (grepl("^DL5_", name)) return("DL5")
    if (grepl("^DD_", name)) {
        parts <- strsplit(name, "_")[[1]]
        return(paste0("DD", parts[2]))
    }
    return(NA)
}

sample_names <- labels(aln)
cruises <- sapply(sample_names, get_cruise)
bloom_status <- ifelse(cruises %in% bloom_cruises, "bloom", "non-bloom")

cat("\nBloom:", sum(bloom_status == "bloom"), "\n")
cat("Non-bloom:", sum(bloom_status == "non-bloom"), "\n")

## Check for any unassigned
if (any(is.na(bloom_status))) {
    cat("WARNING: Unassigned samples:\n")
    print(sample_names[is.na(bloom_status)])
}

## ---- Pi by group ----
aln_bloom <- aln[bloom_status == "bloom"]
aln_nonbloom <- aln[bloom_status == "non-bloom"]

pi_bloom <- nuc.div(aln_bloom)
pi_nonbloom <- nuc.div(aln_nonbloom)
pi_ratio <- pi_bloom / pi_nonbloom

cat("\nPi bloom:", round(pi_bloom, 5), "\n")
cat("Pi non-bloom:", round(pi_nonbloom, 5), "\n")
cat("Pi ratio (bloom/non-bloom):", round(pi_ratio, 3), "\n")

## ---- Hudson's Snn ----
## Need combined distance matrix and group vector
d <- dist.dna(aln, model = "N", pairwise.deletion = TRUE)
pop_factor <- factor(bloom_status)

## Snn test
snn <- Snn(d, pop_factor)
cat("\nSnn:", round(snn[1], 3), "\n")

## Permutation test for Snn
n_perm <- 10000
snn_obs <- snn[1]
snn_perm <- numeric(n_perm)
for (i in 1:n_perm) {
    perm_pop <- sample(pop_factor)
    snn_perm[i] <- Snn(d, perm_pop)[1]
}
p_snn <- mean(snn_perm >= snn_obs)
cat("Snn p-value (", n_perm, "permutations):", p_snn, "\n")

## ---- AMOVA ----
## Bloom vs non-bloom
amova_bloom <- amova(d ~ pop_factor, nperm = 10000)
cat("\n--- AMOVA: Bloom vs Non-bloom ---\n")
print(amova_bloom)
pct_among_bloom <- amova_bloom$varcomp[1, 1] / sum(amova_bloom$varcomp[, 1]) * 100
cat("% variation among groups:", round(pct_among_bloom, 2), "%\n")

## ---- AMOVA by year ----
get_year <- function(name) {
    cruise <- get_cruise(name)
    if (cruise == "DM18") return("2022")
    if (cruise == "DL5") return("2024")
    if (grepl("^DD", cruise)) {
        dd_num <- as.numeric(gsub("DD", "", cruise))
        if (dd_num <= 15) return("2016")
        if (dd_num <= 21) return("2017")
    }
    return(NA)
}

years <- factor(sapply(sample_names, get_year))
cat("\nSamples by year:\n")
print(table(years))

amova_year <- amova(d ~ years, nperm = 10000)
cat("\n--- AMOVA: By Year ---\n")
print(amova_year)
pct_among_year <- amova_year$varcomp[1, 1] / sum(amova_year$varcomp[, 1]) * 100
cat("% variation among years:", round(pct_among_year, 2), "%\n")

## ---- Summary table ----
cat("\n\n===== TAB-DELIMITED SUMMARY =====\n")
cat("Metric\tValue\n")
cat("Bloom n\t", sum(bloom_status == "bloom"), "\n", sep = "")
cat("Non-bloom n\t", sum(bloom_status == "non-bloom"), "\n", sep = "")
cat("Snn\t", round(snn_obs, 3), " (p = ", p_snn, ")\n", sep = "")
cat("Pi bloom\t", round(pi_bloom, 5), "\n", sep = "")
cat("Pi non-bloom\t", round(pi_nonbloom, 5), "\n", sep = "")
cat("Pi ratio (bloom/non-bloom)\t", round(pi_ratio, 3), "\n", sep = "")
cat("AMOVA bloom % among\t", round(pct_among_bloom, 2), "%\n", sep = "")
cat("AMOVA year % among\t", round(pct_among_year, 2), "%\n", sep = "")
