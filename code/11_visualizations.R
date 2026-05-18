# 11_visualizations.R — produce figures summarizing Architecture B clusters.
#
# Outputs (output/figures/, gitignored — move to docs/figures/ to commit):
#   spectra_per_spec_cluster.png   8 lines, mean brightness-normalized
#                                   reflectance per spec cluster
#   spectra_per_final_label.png    17 lines, faceted by spec cluster, so each
#                                   panel shows the sub-cluster spectra (and
#                                   visually confirms they overlap)
#   confusion_matrix.png           Row-normalized RF CV confusion at the 17-
#                                   label level, ordered by spec cluster
#   dendrogram.png                 Spectral Ward (PCs 1-12) with k=8 cut
#   composition_profile.png        Mean Hellinger per final label, top genera
#
# Inputs: data/derived/{final_clusters_B,spectral_clusters,veg_spectra,
#                       spectral_features,composition_genus}.rds

suppressPackageStartupMessages({
  library(tidyverse)
})

fc <- readRDS("data/derived/final_clusters_B.rds")
sc <- readRDS("data/derived/spectral_clusters.rds")
vs <- readRDS("data/derived/veg_spectra.rds")
cg <- readRDS("data/derived/composition_genus.rds")

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

joined <- vs$joined
wl     <- vs$wavelengths

# Indicator colour palette — keep consistent across all figures by spec cluster.
spec_levels <- sort(unique(fc$assignments$spec_cluster))
spec_palette <- setNames(
  scales::hue_pal()(length(spec_levels)),
  spec_levels
)

# Add interpretive labels for the legend.
spec_label_lookup <- fc$spec_summary |>
  dplyr::transmute(
    spec_cluster,
    legend = sprintf("%s — %s (n=%d)", spec_cluster, indicator_genus, n_sites)
  ) |>
  tibble::deframe()

# ============================================================================
# (1) Mean reflectance per spec cluster
# ============================================================================
rfl_cols <- grep("^rfl_band_", names(joined), value = TRUE)
band_num_lookup <- tibble(
  band       = rfl_cols,
  band_num   = as.integer(stringr::str_extract(rfl_cols, "\\d+$"))
) |>
  left_join(wl, by = c("band_num" = "band_number"))

joined_lab <- joined |>
  inner_join(fc$assignments |>
               dplyr::select(site_number, Year, spec_cluster, final_label),
             by = c("site_number", "Year"))

long_spec <- joined_lab |>
  dplyr::select(spec_cluster, final_label, dplyr::all_of(rfl_cols)) |>
  pivot_longer(dplyr::all_of(rfl_cols), names_to = "band", values_to = "rfl") |>
  left_join(band_num_lookup |> dplyr::select(band, center_wavelength_nm),
            by = "band")

mean_per_spec <- long_spec |>
  group_by(spec_cluster, center_wavelength_nm) |>
  summarise(mean_rfl = mean(rfl, na.rm = TRUE), .groups = "drop") |>
  mutate(legend = spec_label_lookup[spec_cluster])

water_bands <- tibble(
  xmin = c(1340, 1800, 2400),
  xmax = c(1450, 1950, 2510),
  label = c("H2O", "H2O", "H2O")
)

p1 <- ggplot(mean_per_spec,
             aes(x = center_wavelength_nm, y = mean_rfl, colour = spec_cluster,
                 group = spec_cluster)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = spec_palette,
                      labels = spec_label_lookup,
                      name = "Spec cluster — indicator (n)") +
  labs(x = "Wavelength (nm)",
       y = "Brightness-normalized reflectance",
       title = "Per-spec-cluster mean spectrum (Architecture B, k=8)",
       subtitle = "Grey bands = water absorption regions (masked in feature engineering)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")
ggsave("output/figures/spectra_per_spec_cluster.png", p1,
       width = 10, height = 6, dpi = 150)

# ============================================================================
# (2) Sub-cluster spectra faceted by spec cluster
# ============================================================================
mean_per_final <- long_spec |>
  group_by(spec_cluster, final_label, center_wavelength_nm) |>
  summarise(mean_rfl = mean(rfl, na.rm = TRUE), .groups = "drop")

p2 <- ggplot(mean_per_final,
             aes(x = center_wavelength_nm, y = mean_rfl, colour = final_label,
                 group = final_label)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ spec_cluster, ncol = 4) +
  labs(x = "Wavelength (nm)",
       y = "Brightness-normalized reflectance",
       title = "Sub-cluster spectra within each spec cluster",
       subtitle = "Sub-clusters share their parent's spectral signature by construction; "
                  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))
ggsave("output/figures/spectra_per_final_label.png", p2,
       width = 12, height = 7, dpi = 150)

# ============================================================================
# (3) Confusion matrix (row-normalized) at the 17-label level
# ============================================================================
cm <- fc$final_eval$confusion
cm_long <- as_tibble(as.table(cm))
names(cm_long) <- c("truth", "pred", "n")
cm_long <- cm_long |>
  group_by(truth) |>
  mutate(prop = if (sum(n) > 0) n / sum(n) else 0) |>
  ungroup()

# Order labels by spec_cluster + final_label for visual blocking.
final_order <- fc$final_summary |>
  arrange(spec_cluster, final_label) |>
  pull(final_label)
cm_long <- cm_long |>
  mutate(truth = factor(truth, levels = final_order),
         pred  = factor(pred,  levels = final_order))

p3 <- ggplot(cm_long, aes(x = pred, y = truth, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = ifelse(prop >= 0.05, sprintf("%.2f", prop), "")),
            size = 2.5, colour = "white") +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1)) +
  scale_y_discrete(limits = rev) +
  labs(x = "Predicted label", y = "True label",
       title = "RF CV confusion matrix (row-normalized)",
       subtitle = "Within-spec-cluster sub-classes are predicted as the parent's modal sub.",
       fill = "P(pred | truth)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
ggsave("output/figures/confusion_matrix.png", p3,
       width = 9, height = 8, dpi = 150)

# ============================================================================
# (4) Spectral Ward dendrogram with k=8 cut
# ============================================================================
png("output/figures/dendrogram.png", width = 1400, height = 700, res = 110)
op <- par(mar = c(2, 4, 3, 1))
hc <- sc$variant_A$hclust
plot(hc, labels = FALSE,
     main = "Spectral Ward hierarchical (PCs 1-12), cut at k=8",
     sub = "", xlab = "", ylab = "Ward height")
rect.hclust(hc, k = 8, border = spec_palette[spec_levels])
par(op)
dev.off()

# ============================================================================
# (5) Composition profile (top genera per final label)
# ============================================================================
nonsp_cover <- c("Other_Forb_cover", "Other_Graminoid_cover", "NPV_cover",
                 "Bare_cover", "Other_Moss_Lichen_cover",
                 "Other_Deciduous_Shrub_cover")
hell_named <- cg$hellinger |> dplyr::select(-dplyr::any_of(nonsp_cover))
hell_long_named <- hell_named |>
  pivot_longer(-c(site_number, Year), names_to = "feature", values_to = "h") |>
  mutate(feature = stringr::str_replace(feature, "_cover$", "")) |>
  inner_join(fc$assignments |>
               dplyr::select(site_number, Year, spec_cluster, final_label),
             by = c("site_number", "Year"))

top_genera_per_label <- hell_long_named |>
  group_by(spec_cluster, final_label, feature) |>
  summarise(mean_h = mean(h), .groups = "drop") |>
  group_by(spec_cluster, final_label) |>
  slice_max(mean_h, n = 6, with_ties = FALSE) |>
  ungroup()

p5 <- ggplot(top_genera_per_label,
             aes(x = reorder(feature, mean_h), y = mean_h, fill = spec_cluster)) +
  geom_col() +
  scale_fill_manual(values = spec_palette, guide = "none") +
  facet_wrap(~ final_label, scales = "free_y", ncol = 4) +
  coord_flip() +
  labs(x = NULL, y = "Mean Hellinger value (named genera only)",
       title = "Top 6 named genera per final label") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.y = element_text(size = 7))
ggsave("output/figures/composition_profile.png", p5,
       width = 13, height = 10, dpi = 150)

cat("Wrote:\n")
for (f in c("spectra_per_spec_cluster.png", "spectra_per_final_label.png",
            "confusion_matrix.png", "dendrogram.png", "composition_profile.png")) {
  p <- file.path("output/figures", f)
  cat(sprintf("  %s  (%s)\n", p,
              if (file.exists(p)) sprintf("%.1f KB", file.size(p)/1024) else "missing!"))
}
