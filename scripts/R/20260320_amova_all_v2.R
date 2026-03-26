## 20260320_amova_all_v2.R
## AMOVA: bloom vs non-bloom, by year, DD11 vs DL5
## Requires: all_120_mitobim_filtered.fasta
## Run locally in RStudio

library(ape)
library(pegas)

## ---- Read filtered alignment ----
aln <- read.FASTA("all_120_mitobim_filtered.fasta")
mat <- as.character(as.matrix(aln))

## Remove DM18_27
exclude <- "DM18_27"
keep_samp <- !labels(aln) %in% exclude
aln_119 <- as.DNAbin(mat[keep_samp, ])
cat("Working dataset:", nrow(as.matrix(aln_119)), "samples\n")

## ---- Assign cruise, bloom status, year ----
eg_cruise <- c(
    EG1 = "DD4", EG2 = "DD7", EG3 = "DD8", EG5 = "DD10",
    EG6 = "DD11", EG7 = "DD11", EG8 = "DD11",
    EG9 = "DM18", EG10 = "DM18", EG11 = "DM18", EG12 = "DM18",
    EG13 = "DM18", EG14 = "DM18", EG15 = "DM18", EG16 = "DM18"
)

bloom_cruises <- c("DD7", "DD11", "DD17", "DD21", "DM18")

get_cruise <- function(name) {
    if (name %in% names(eg_cruise)) return(eg_cruise[name])
    if (grepl("^DM18_", name)) return("DM18")
    if (grepl("^DL5_", name)) return("DL5")
    if (grepl("^DD_", name)) {
        parts <- strsplit(name, "_")[[1]]
        return(paste0("DD", parts[2]))
    }
    return(NA)
}

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

sample_names <- labels(aln_119)
cruises <- sapply(sample_names, get_cruise)
bloom_status <- factor(ifelse(cruises %in% bloom_cruises, "bloom", "non-bloom"))
years <- factor(sapply(sample_names, get_year))

cat("\nBloom status:\n")
print(table(bloom_status))
cat("\nYear:\n")
print(table(years))

## ---- Distance matrix ----
d <- dist.dna(aln_119, model = "N", pairwise.deletion = TRUE)

## ---- AMOVA 1: Bloom vs Non-bloom ----
cat("\n========== AMOVA: Bloom vs Non-bloom ==========\n")
amova_bloom <- pegas::amova(d ~ bloom_status, nperm = 10000)
print(amova_bloom)
pct_bloom <- amova_bloom$varcomp[1, 1] / sum(amova_bloom$varcomp[, 1]) * 100
p_bloom <- amova_bloom$varcomp[1, 2]
cat("% among groups:", round(pct_bloom, 2), "%\n")
cat("p =", p_bloom, "\n")

## ---- AMOVA 2: By Year ----
cat("\n========== AMOVA: By Year ==========\n")
amova_year <- pegas::amova(d ~ years, nperm = 10000)
print(amova_year)
pct_year <- amova_year$varcomp[1, 1] / sum(amova_year$varcomp[, 1]) * 100
p_year <- amova_year$varcomp[1, 2]
cat("% among years:", round(pct_year, 2), "%\n")
cat("p =", p_year, "\n")

## ---- AMOVA 3: DD11 vs DL5 ----
cat("\n========== AMOVA: DD11 vs DL5 ==========\n")
dd11_dl5 <- cruises %in% c("DD11", "DL5")
aln_dd11_dl5 <- aln_119[dd11_dl5, ]
d_dd11_dl5 <- dist.dna(aln_dd11_dl5, model = "N", pairwise.deletion = TRUE)
cruise_factor <- factor(cruises[dd11_dl5])

cat("DD11:", sum(cruise_factor == "DD11"), "\n")
cat("DL5:", sum(cruise_factor == "DL5"), "\n")

amova_dd11dl5 <- pegas::amova(d_dd11_dl5 ~ cruise_factor, nperm = 10000)
print(amova_dd11dl5)
pct_dd11dl5 <- amova_dd11dl5$varcomp[1, 1] / sum(amova_dd11dl5$varcomp[, 1]) * 100
p_dd11dl5 <- amova_dd11dl5$varcomp[1, 2]
cat("% among groups:", round(pct_dd11dl5, 2), "%\n")
cat("p =", p_dd11dl5, "\n")

## ---- Summary table ----
cat("\n\n===== TAB-DELIMITED SUMMARY =====\n")
cat("Comparison\tGroups\t% among\tp\n")
cat("Bloom vs Non-bloom\t", paste0("bloom(", sum(bloom_status == "bloom"), ") vs non-bloom(", sum(bloom_status == "non-bloom"), ")"),
    "\t", round(pct_bloom, 2), "%\t", p_bloom, "\n", sep = "")
cat("By Year\t", paste(names(table(years)), collapse = "/"),
    "\t", round(pct_year, 2), "%\t", p_year, "\n", sep = "")
cat("DD11 vs DL5\t", paste0("DD11(", sum(cruise_factor == "DD11"), ") vs DL5(", sum(cruise_factor == "DL5"), ")"),
    "\t", round(pct_dd11dl5, 2), "%\t", p_dd11dl5, "\n", sep = "")
