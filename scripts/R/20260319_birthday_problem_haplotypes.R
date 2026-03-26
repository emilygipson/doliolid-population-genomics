# 20260319_birthday_problem_haplotypes.R
# Birthday problem framing: probability of sampling a duplicate haplotype
# as a function of sample size, for different numbers of circulating haplotypes

library(ggplot2)

# Range of sample sizes
n_samples <- 1:500

# Different assumptions about number of circulating haplotypes (k)
# Conservative (k = 10,000) to realistic (k = 200,000)
k_values <- c(10000, 50000, 200000)
k_labels <- c("k = 10,000\n(very conservative)",
               "k = 50,000\n(conservative)",
               "k = 200,000\n(realistic for Ne ~ 500,000)")

# Birthday problem: P(at least one match) ≈ 1 - exp(-n^2 / (2k))
results <- data.frame()
for (i in seq_along(k_values)) {
  k <- k_values[i]
  p_match <- 1 - exp(-n_samples^2 / (2 * k))
  results <- rbind(results, data.frame(
    n = n_samples,
    p = p_match,
    k = k_labels[i],
    stringsAsFactors = FALSE
  ))
}

results$k <- factor(results$k, levels = k_labels)

# Find the 50% threshold for each k
n_at_50 <- sapply(k_values, function(k) {
  ceiling(sqrt(2 * k * log(2)))
})

threshold_df <- data.frame(
  k = factor(k_labels, levels = k_labels),
  n50 = n_at_50,
  label = paste0("n = ", format(n_at_50, big.mark = ","))
)

p <- ggplot(results, aes(x = n, y = p, color = k)) +
  geom_line(linewidth = 1.2) +
  # Our sample size
  geom_vline(xintercept = 120, linetype = "dashed", color = "grey40",
             linewidth = 0.6) +
  annotate("text", x = 130, y = 0.95,
           label = "Our sample\n(n = 120)", size = 3.5, hjust = 0,
           color = "grey40") +
  # 50% line
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey60",
             linewidth = 0.4) +
  # 50% threshold labels
  geom_point(data = threshold_df, aes(x = n50, y = 0.5, color = k),
             size = 3, shape = 16) +
  geom_text(data = threshold_df, aes(x = n50, y = 0.55, label = label, color = k),
            size = 3, hjust = 0, nudge_x = 5, show.legend = FALSE) +
  scale_color_manual(values = c("#E64B35", "#4DBBD5", "#00A087"),
                     name = "Circulating haplotypes") +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = c(0, 100, 120, 200, 300, 400, 500)) +
  labs(x = "Number of individuals sampled",
       y = "Probability of at least one shared haplotype pair",
       title = "How many samples before you expect a duplicate?",
       subtitle = expression(paste("Ne" %~% "400,000 - 600,000 from nuclear data; ",
                                    pi, " = 0.0066 across 11,016 sites"))) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    panel.grid.minor = element_blank()
  )

ggsave("20260319_birthday_problem_haplotypes.png", p, width = 9, height = 6,
       dpi = 300, bg = "white")
