#!/usr/bin/env Rscript
# 20260312_bsp_plot.R
# Clean BSP plot from Tracer TSV export

library(ggplot2)

# Read Tracer BSP export (tab-delimited)
# Tracer exports columns: Time, Mean, Median, Upper, Lower (or similar)
bsp <- read.delim("dgeg_bsp_105.tsv", header = TRUE)

# Check column names and adjust if needed
cat("Columns found:", paste(colnames(bsp), collapse = ", "), "\n")
cat("Rows:", nrow(bsp), "\n")
head(bsp)

# Tracer typically exports: Time, Mean, Median, Upper, Lower
# Rename if needed — adjust these to match your actual column names
colnames(bsp) <- tolower(colnames(bsp))

# Try to identify columns
time_col <- grep("time", colnames(bsp), value = TRUE)[1]
mean_col <- grep("mean", colnames(bsp), value = TRUE)[1]
median_col <- grep("median", colnames(bsp), value = TRUE)[1]
upper_col <- grep("upper", colnames(bsp), value = TRUE)[1]
lower_col <- grep("lower", colnames(bsp), value = TRUE)[1]

cat("\nUsing columns:\n")
cat("  Time:", time_col, "\n")
cat("  Mean:", mean_col, "\n")
cat("  Median:", median_col, "\n")
cat("  Upper:", upper_col, "\n")
cat("  Lower:", lower_col, "\n")

df <- data.frame(
  time = bsp[[time_col]],
  mean = bsp[[mean_col]],
  median = bsp[[median_col]],
  upper = bsp[[upper_col]],
  lower = bsp[[lower_col]]
)

# Plot
p <- ggplot(df, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = median), color = "steelblue", linewidth = 1.2) +
  scale_y_log10(labels = scales::comma) +
  annotation_logticks(sides = "l") +
  labs(
    x = "Time (generations before present)",
    y = expression(N[e]),
    title = expression(paste("Bayesian Skyline Plot — 105 ", italic("D. gegenbauri"), " mitogenomes")),
    subtitle = expression(paste("Mutation rate: 3.0 × ", 10^{-8}, " per site per generation"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave("bsp_105_plot.pdf", p, width = 10, height = 6)
ggsave("bsp_105_plot.png", p, width = 10, height = 6, dpi = 300)

cat("\nPlots saved: bsp_105_plot.pdf, bsp_105_plot.png\n")

# Print some key values
cat("\nKey BSP values:\n")
cat("Most recent Ne (median):", format(df$median[nrow(df)], big.mark = ","), "\n")
cat("Most ancient Ne (median):", format(df$median[1], big.mark = ","), "\n")
cat("Expansion ratio:", round(df$median[nrow(df)] / df$median[1], 1), "x\n")
