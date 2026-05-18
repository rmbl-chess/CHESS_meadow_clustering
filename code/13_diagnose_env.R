# 13_diagnose_env.R — diagnose whether environmental signal already separates
# the Architecture B spec clusters, or whether the clusters are env-mixed
# (in which case env would help if added as clustering features).
#
# Approach:
#   1. Profile each spec cluster's snow-free DOY distribution: mean, sd,
#      IQR, range.
#   2. Visualize as a boxplot ordered by mean DOY.
#   3. Compute an env-coherence score: between-cluster variance over total
#      variance (eta²; ANOVA-style). High eta² = clusters already separate
#      on env. Low eta² = clusters span the env gradient indiscriminately,
#      and adding env as a clustering feature is likely to help.
#   4. Also check sub-cluster level: does sub-clustering by composition
#      reveal env structure that the spec clusters lacked?
#
# Inputs:
#   data/derived/final_clusters_B.rds
#   data/derived/environment.rds
# Outputs:
#   output/figures/env_per_cluster.png   boxplot of snow_free_doy by label
#   stdout: per-cluster summary + eta² scores

suppressPackageStartupMessages({
  library(tidyverse)
})

fc  <- readRDS("data/derived/final_clusters_B.rds")
env <- readRDS("data/derived/environment.rds")

asg <- fc$assignments |>
  dplyr::inner_join(env, by = c("site_number", "Year"))

cat(sprintf("Joined %d sites with both cluster + env\n", nrow(asg)))

eta_sq <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- factor(g[ok])
  grand <- mean(x)
  ss_total   <- sum((x - grand)^2)
  ss_between <- sum(table(g) * (tapply(x, g, mean) - grand)^2)
  ss_between / ss_total
}

# --- Spec cluster level -----------------------------------------------------
spec_profile <- asg |>
  group_by(spec_cluster) |>
  summarise(
    n       = dplyr::n(),
    mean    = mean(snow_free_doy, na.rm = TRUE),
    sd      = sd(snow_free_doy, na.rm = TRUE),
    median  = median(snow_free_doy, na.rm = TRUE),
    q25     = quantile(snow_free_doy, .25, na.rm = TRUE),
    q75     = quantile(snow_free_doy, .75, na.rm = TRUE),
    min     = min(snow_free_doy, na.rm = TRUE),
    max     = max(snow_free_doy, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(fc$spec_summary |>
              dplyr::select(spec_cluster, indicator_species),
            by = "spec_cluster") |>
  arrange(mean)

cat("\n=== Snow-free DOY per spectral cluster (ordered by mean) ===\n")
print(spec_profile, n = Inf, width = Inf)

eta_spec  <- eta_sq(asg$snow_free_doy, asg$spec_cluster)
eta_final <- eta_sq(asg$snow_free_doy, asg$final_label)
overall_var <- var(asg$snow_free_doy, na.rm = TRUE)
overall_sd  <- sd(asg$snow_free_doy, na.rm = TRUE)
cat(sprintf("\nOverall snow_free_doy: sd=%.1f days, range=%.0f-%.0f\n",
            overall_sd, min(asg$snow_free_doy, na.rm=TRUE),
            max(asg$snow_free_doy, na.rm=TRUE)))
cat(sprintf("Eta² (variance explained by spec_cluster):   %.3f\n", eta_spec))
cat(sprintf("Eta² (variance explained by final_label):    %.3f\n", eta_final))
cat("\nInterpretation:\n")
cat("  >0.5  clusters already separate strongly on env (env adds little)\n")
cat("  0.2-0.5  partial env separation (env may sharpen some clusters)\n")
cat("  <0.2  clusters span env indiscriminately (env should help)\n")

# --- Plot -------------------------------------------------------------------
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

spec_order <- spec_profile$spec_cluster
asg_plot <- asg |>
  mutate(spec_cluster = factor(spec_cluster, levels = spec_order),
         spec_label = sprintf("%s — %s",
                              spec_cluster,
                              stringr::str_replace_all(
                                fc$spec_summary$indicator_species[match(
                                  as.character(spec_cluster),
                                  fc$spec_summary$spec_cluster)],
                                "_", " "))) |>
  mutate(spec_label = factor(spec_label,
                             levels = unique(spec_label[order(spec_cluster)])))

p <- ggplot(asg_plot, aes(x = spec_label, y = snow_free_doy)) +
  geom_violin(aes(fill = spec_label), alpha = 0.4, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.4, size = 0.7) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = "Snow-free DOY (1993–2022 mean)",
       title = "Snow-free date climatology per spec cluster",
       subtitle = sprintf("Eta² (spec) = %.3f   Eta² (final) = %.3f   Overall SD = %.1f days",
                          eta_spec, eta_final, overall_sd)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/figures/env_per_cluster.png", p,
       width = 11, height = 6, dpi = 150)
cat("\nWrote output/figures/env_per_cluster.png\n")
