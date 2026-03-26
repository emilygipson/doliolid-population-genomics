## 20260319_pairwise_distance_figures.R
## Pairwise distance histograms for lab meeting and paper
## Uses pre-computed distances from the 120-sample alignment
## Run locally in RStudio

## ---- Pre-computed pairwise distances ----
## These are nucleotide differences from the 11,016-site filtered alignment
## 119 samples (DM18_27 excluded)

## Full dataset (119 samples, 7021 pairs)
full_dists <- read.delim("pairwise_dists_119.tsv")$distance

## DD11 cast 3 (12 samples from one net haul, 66 pairs)
dd11_dists <- read.delim("dd11_cast3_dists.tsv")$distance

## DM18 station 5 / 25m (5 nurses, 10 pairs)
stn5_dists <- read.delim("dm18_stn5_dists.tsv")$distance

## ---- Figure 1: Full dataset pairwise histogram ----
## For paper and slide 7
png("pairwise_hist_full_119.png", width = 8, height = 5, units = "in", res = 300)
par(mar = c(5, 4, 3, 1))
hist(full_dists, breaks = seq(0, max(full_dists) + 5, by = 2), 
     col = "steelblue", border = "white",
     main = "", xlab = "Pairwise nucleotide differences",
     ylab = "Number of pairs", las = 1)
mtext("119 D. gegenbauri mitogenomes | 11,016 filtered sites", 
      side = 3, line = 0.5, cex = 0.9)
abline(v = mean(full_dists), col = "red", lwd = 2, lty = 2)
text(mean(full_dists) + 3, par("usr")[4] * 0.9, 
     paste0("Mean = ", round(mean(full_dists), 1)), 
     col = "red", adj = 0, cex = 0.85)
abline(v = min(full_dists), col = "darkorange", lwd = 2, lty = 3)
text(min(full_dists) + 2, par("usr")[4] * 0.8, 
     paste0("Min = ", min(full_dists)), 
     col = "darkorange", adj = 0, cex = 0.85)
dev.off()
cat("Saved: pairwise_hist_full_119.png\n")

## ---- Figure 2: Within-tow comparison panel ----
## For paper and slide 6
png("within_tow_pairwise_panel.png", width = 10, height = 4, units = "in", res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

## Panel A: DD11 cast 3
hist(dd11_dists, breaks = seq(0, max(dd11_dists) + 5, by = 3),
     col = "#E69F00", border = "white",
     main = "DD11 cast 3 (12 samples, 1 net)",
     xlab = "Pairwise differences", ylab = "Pairs", las = 1, cex.main = 0.95)
abline(v = min(dd11_dists), col = "red", lwd = 2, lty = 2)
text(min(dd11_dists) + 2, par("usr")[4] * 0.85,
     paste0("Min = ", min(dd11_dists)), col = "red", adj = 0, cex = 0.8)
mtext("12/12 unique", side = 3, line = -1.5, cex = 0.75, font = 2)

## Panel B: DM18 station 5
hist(stn5_dists, breaks = seq(0, max(stn5_dists) + 5, by = 5),
     col = "#56B4E9", border = "white",
     main = "DM18 stn 5 (5 nurses, 1 net)",
     xlab = "Pairwise differences", ylab = "Pairs", las = 1, cex.main = 0.95)
abline(v = min(stn5_dists), col = "red", lwd = 2, lty = 2)
text(min(stn5_dists) + 2, par("usr")[4] * 0.85,
     paste0("Min = ", min(stn5_dists)), col = "red", adj = 0, cex = 0.8)
mtext("5/5 unique", side = 3, line = -1.5, cex = 0.75, font = 2)

## Panel C: Full dataset for context
hist(full_dists, breaks = seq(0, max(full_dists) + 5, by = 2),
     col = "gray60", border = "white",
     main = "Full dataset (119 samples)",
     xlab = "Pairwise differences", ylab = "Pairs", las = 1, cex.main = 0.95)
abline(v = 0, col = "red", lwd = 2)
text(5, par("usr")[4] * 0.85, "Zero = clones", col = "red", adj = 0, cex = 0.8)
mtext("119/119 unique", side = 3, line = -1.5, cex = 0.75, font = 2)

dev.off()
cat("Saved: within_tow_pairwise_panel.png\n")

## ---- Figure 3: Simplified version for slide 6 ----
## Just DD11 cast 3 with annotation
png("dd11_cast3_pairwise.png", width = 7, height = 5, units = "in", res = 300)
par(mar = c(5, 4, 3, 1))
hist(dd11_dists, breaks = seq(0, max(dd11_dists) + 5, by = 3),
     col = "#E69F00", border = "white",
     main = "", xlab = "Pairwise nucleotide differences",
     ylab = "Number of pairs", las = 1, xlim = c(0, 100))
mtext("DD11 cast 3: 12 individuals from one net haul during active bloom",
      side = 3, line = 1, cex = 0.9)
mtext("3 nurses + 2 phorozooids + 5 gonozooids + 2 unknown (original) + 3 nurses (EG)",
      side = 3, line = 0, cex = 0.75, col = "gray40")
abline(v = 0, col = "red", lwd = 3)
text(3, par("usr")[4] * 0.9, "If clones existed,\npairs would cluster here",
     col = "red", adj = 0, cex = 0.85)
arrows(20, par("usr")[4] * 0.7, min(dd11_dists), par("usr")[4] * 0.3,
       col = "darkorange", lwd = 2)
text(21, par("usr")[4] * 0.75, paste0("Closest pair = ", min(dd11_dists), " differences"),
     col = "darkorange", adj = 0, cex = 0.85)
dev.off()
cat("Saved: dd11_cast3_pairwise.png\n")

## ---- Summary stats ----
cat("\n=== Summary ===\n")
cat("Full 119: mean =", round(mean(full_dists), 1), 
    ", min =", min(full_dists), ", max =", max(full_dists), "\n")
cat("DD11 cast 3: mean =", round(mean(dd11_dists), 1),
    ", min =", min(dd11_dists), ", max =", max(dd11_dists), "\n")
cat("DM18 stn5: mean =", round(mean(stn5_dists), 1),
    ", min =", min(stn5_dists), ", max =", max(stn5_dists), "\n")
