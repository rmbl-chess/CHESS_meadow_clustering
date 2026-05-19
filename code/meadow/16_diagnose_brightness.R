# 15_diagnose_brightness.R — show how the k_spec=12 clusters distribute along
# PC1 (brightness/albedo, 95% of spectral variance) and PC2 (greenness,
# 3.34%). Parallels 13_diagnose_env.R but for the dominant spectral axes.
#
# A clean separation on PC1 means brightness contributes to discrimination;
# a clean separation on PC2 means greenness/chlorophyll absorption does;
# overlapping clusters on PC1 mean other features (PCs 3+, env) are doing
# the work.
#
# Inputs:  data/derived/final_clusters_B.rds, .../spectral_features.rds
# Outputs: output/figures/pc1_pc2_per_cluster.png  + per-cluster summary

suppressPackageStartupMessages({
  library(tidyverse)
})

fc        <- readRDS("data/derived/final_clusters_B.rds")
spec_feat <- readRDS("data/derived/spectral_features.rds")$features

asg <- fc$assignments |>
  dplyr::inner_join(spec_feat |> dplyr::select(site_number, Year,
                                               spec_PC01, spec_PC02),
                    by = c("site_number", "Year"))

eta_sq <- function(x, g) {
  ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- factor(g[ok])
  grand <- mean(x)
  sum(table(g) * (tapply(x, g, mean) - grand)^2) / sum((x - grand)^2)
}

eta_pc1 <- eta_sq(asg$spec_PC01, asg$spec_cluster)
eta_pc2 <- eta_sq(asg$spec_PC02, asg$spec_cluster)

per_cl <- asg |>
  group_by(spec_cluster) |>
  summarise(
    n        = dplyr::n(),
    pc1_mean = mean(spec_PC01),
    pc1_sd   = sd(spec_PC01),
    pc2_mean = mean(spec_PC02),
    pc2_sd   = sd(spec_PC02),
    .groups  = "drop"
  ) |>
  left_join(fc$spec_summary |> dplyr::select(spec_cluster, indicator_species),
            by = "spec_cluster") |>
  arrange(pc1_mean)

cat(sprintf("Eta²(PC1) = %.3f   Eta²(PC2) = %.3f\n",
            eta_pc1, eta_pc2))
cat(sprintf("(reference) Eta²(snow_free_doy) at this k_spec = ~0.70\n\n"))
print(per_cl, n = Inf, width = Inf)

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

cluster_order <- per_cl$spec_cluster
asg_plot <- asg |>
  mutate(
    spec_cluster = factor(spec_cluster, levels = cluster_order),
    label = sprintf("%s — %s",
                    spec_cluster,
                    stringr::str_replace_all(
                      fc$spec_summary$indicator_species[match(
                        as.character(spec_cluster),
                        fc$spec_summary$spec_cluster)],
                      "_", " "))
  ) |>
  mutate(label = factor(label, levels = unique(label[order(spec_cluster)])))

p1 <- ggplot(asg_plot, aes(x = label, y = spec_PC01)) +
  geom_violin(aes(fill = label), alpha = 0.4, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.4, size = 0.7) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  labs(x = NULL, y = "PC1 (brightness, 95% of spectral variance)",
       title = "PC1 distribution per spec cluster",
       subtitle = sprintf("Eta²(PC1) = %.3f  vs  Eta²(snow_free_doy) = 0.70",
                          eta_pc1)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p2 <- ggplot(asg_plot, aes(x = label, y = spec_PC02)) +
  geom_violin(aes(fill = label), alpha = 0.4, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.4, size = 0.7) +
  scale_fill_brewer(palette = "Set3", guide = "none") +
  labs(x = NULL, y = "PC2 (greenness/red-edge, 3.34%)",
       title = "PC2 distribution per spec cluster",
       subtitle = sprintf("Eta²(PC2) = %.3f", eta_pc2)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("output/figures/pc1_per_cluster.png", p1,
       width = 12, height = 5, dpi = 150)
ggsave("output/figures/pc2_per_cluster.png", p2,
       width = 12, height = 5, dpi = 150)
cat("\nWrote output/figures/pc1_per_cluster.png and pc2_per_cluster.png\n")
