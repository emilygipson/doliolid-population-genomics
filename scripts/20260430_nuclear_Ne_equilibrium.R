# Nuclear Ne equilibrium estimates from ANGSD genome-wide diversity
# Inputs: genome-wide tW and tP from nuc118_thetas.thetas.idx.pestPG
# Author: Emily Gipson
# Date: 2026-04-30

.libPaths(c("/scratch/eeg37520/Rlibs", .libPaths()))

# Inputs from section 4 of 20260430_nuc118_recon report
tW_per_site <- 0.004750
tP_per_site <- 0.004462
nSites      <- 6528505
TajimaD     <- -0.5551

# Mutation rate range (per site per generation)
# Note: lower mu produces higher Ne. Labels describe the Ne result, not mu.
mu_values <- c(low_Ne = 3.0e-9, mid_Ne = 2.5e-9, high_Ne = 2.0e-9)

# Compute Ne under Wright-Fisher equilibrium
ne_table <- data.frame(
  estimator       = rep(c("Ne = tW / 4mu", "Ne = pi / 4mu"), each = length(mu_values)),
  theta_per_site  = rep(c(tW_per_site, tP_per_site), each = length(mu_values)),
  mu_label        = rep(names(mu_values), times = 2),
  mu              = rep(mu_values, times = 2),
  Ne              = NA_real_
)
ne_table$Ne <- ne_table$theta_per_site / (4 * ne_table$mu)

# Round Ne to nearest 1000 for reporting
ne_table$Ne_rounded <- round(ne_table$Ne, -3)

# Add metadata rows (input values used)
meta <- data.frame(
  estimator       = c("nSites_analyzed", "Tajima_D_weighted_approx", "tW_per_site_input", "pi_per_site_input"),
  theta_per_site  = c(NA, NA, tW_per_site, tP_per_site),
  mu_label        = NA,
  mu              = NA,
  Ne              = c(nSites, TajimaD, NA, NA),
  Ne_rounded      = c(nSites, TajimaD, NA, NA)
)

out <- rbind(meta, ne_table)

# Write CSV
out_file <- "20260430_nuclear_Ne_equilibrium.csv"
write.csv(out, out_file, row.names = FALSE)

# Print to console for inspection
cat("Inputs:\n")
cat(sprintf("  tW per site:  %.6f\n", tW_per_site))
cat(sprintf("  pi per site:  %.6f\n", tP_per_site))
cat(sprintf("  nSites:       %d\n", nSites))
cat(sprintf("  Tajima's D:   %.4f (nSites-weighted approx)\n\n", TajimaD))

cat("Ne estimates (rounded to nearest 1000):\n\n")
print_table <- ne_table[, c("estimator", "mu_label", "mu", "Ne_rounded")]
print(print_table, row.names = FALSE)

cat(sprintf("\nWritten to: %s\n", out_file))
