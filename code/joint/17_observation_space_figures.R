# 17_observation_space_figures.R — situate field-sampling rounds (2018 / 2025 /
# 2026) against the random basin pixels we use to guide additional sampling, in
# both spectral (PCA) and species-composition space.
#
# Three figures:
#   obs_spectral_space.pdf     PC ordination panels: random basin pixels (grey)
#                              + training sites colored by round, shaped by
#                              meadow/shrub. "Where do field obs sit vs the
#                              basin cloud, and is 2026 covering new territory?"
#   obs_spectral_leverage.pdf  PC2-PC3 with random pixels colored by leverage
#                              (the sampling-priority surface) + training overlay
#                              by round. "Did the rounds target high-leverage
#                              spectral gaps, and where are the gaps now?"
#   obs_species_space.pdf      PCoA on meadow species-Hellinger composition,
#                              training colored by round, with each class
#                              centroid sized by basin prevalence
#                              (predicted_n_pixels). "Which communities each
#                              round sampled vs where the basin area is."
#
# Inputs (all on the current 2025-fit PCA basis):
#   data/derived/joint_training_set.rds        (spec_PC01..20, Year, class_type)
#   data/derived/inference_predictions.csv     (random pixel PCs, current basis,
#                                               written by 03 — NOT the stale 6k file)
#   data/derived/sampling_priority.gpkg        (leverage per pixel)
#   data/derived/composition_species.rds       (meadow Hellinger)
#   data/derived/final_clusters_B.rds          (meadow class membership)
#   data/derived/punch_list.csv                (predicted_n_pixels per class)
# Outputs: docs/figures/obs_{spectral_space,spectral_leverage,species_space}.pdf
#
# Base ggplot only (no patchwork/ggnewscale/ggrepel) to avoid new renv deps.

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})
dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)

round_cols <- c("2018" = "#1b9e77", "2025" = "#7570b3", "2026" = "#d95f02")
round_lab  <- function(v) factor(v, levels = c(2018L, 2025L, 2026L),
                                 labels = names(round_cols))

# --- Training sites in spectral space --------------------------------------
jt <- readRDS("data/derived/joint_training_set.rds")$training
pc_cols <- sprintf("spec_PC%02d", 1:20)
train <- jt |>
  dplyr::transmute(site_number, Year, round = round_lab(Year),
                   class_type = stringr::str_to_title(class_type),
                   dplyr::across(dplyr::all_of(pc_cols)))

# --- Random basin pixels (PCs + leverage) ----------------------------------
# inference_predictions.csv carries the CURRENT-basis PCs (03 refreshes them
# from the mosaics) — do NOT use the cached 6k file, whose PCs are stale-basis.
pix <- readr::read_csv("data/derived/inference_predictions.csv",
                       show_col_types = FALSE) |>
  dplyr::mutate(key = paste(round(x_utm, 1), round(y_utm, 1), sep = "_"))
sp <- sf::st_read("data/derived/sampling_priority.gpkg", quiet = TRUE) |>
  sf::st_drop_geometry() |>
  dplyr::mutate(key = paste(round(x_utm, 1), round(y_utm, 1), sep = "_")) |>
  dplyr::select(key, leverage)
rand <- pix |> dplyr::inner_join(sp, by = "key") |>
  dplyr::select(dplyr::all_of(pc_cols), leverage)
cat(sprintf("Random pixels joined to leverage: %d\n", nrow(rand)))

# Reshape PC pairs into a long, facet-able frame.
pairs <- tibble::tribble(
  ~panel,            ~xpc,        ~ypc,
  "PC1 vs PC2",      "spec_PC01", "spec_PC02",
  "PC2 vs PC3",      "spec_PC02", "spec_PC03",
  "PC2 vs PC4",      "spec_PC02", "spec_PC04",
  "PC3 vs PC4",      "spec_PC03", "spec_PC04")
long_pairs <- function(df, extra = character(0)) {
  purrr::pmap_dfr(pairs, function(panel, xpc, ypc) {
    df |> dplyr::transmute(panel = panel, x = .data[[xpc]], y = .data[[ypc]],
                           dplyr::across(dplyr::all_of(extra)))
  }) |> dplyr::mutate(panel = factor(panel, levels = pairs$panel))
}
rand_long  <- long_pairs(rand)
train_long <- long_pairs(train, c("round", "class_type"))

# ============================================================================
# FIGURE 1 — spectral PC ordination, rounds over the basin cloud
# ============================================================================
fig1 <- ggplot() +
  geom_point(data = rand_long, aes(x, y), colour = "grey80",
             size = 0.35, alpha = 0.45) +
  geom_point(data = train_long, aes(x, y, colour = round, shape = class_type),
             size = 1.5, alpha = 0.9) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_colour_manual(values = round_cols, name = "Field round") +
  scale_shape_manual(values = c(Meadow = 16, Shrub = 17), name = "Type") +
  labs(x = NULL, y = NULL,
       title = "Field-sampling rounds in spectral (PCA) space vs the random basin pixels",
       subtitle = paste("Grey = 5,064 basin pixels guiding sampling. Coloured = 955 training sites by",
                        "round (circle meadow, triangle shrub). PC1 = brightness.")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave("docs/figures/obs_spectral_space.pdf", fig1, width = 11, height = 9,
       device = cairo_pdf)
cat("Wrote docs/figures/obs_spectral_space.pdf\n")

# ============================================================================
# FIGURE 2 — leverage surface (PC2-PC3) with rounds overlaid
# ============================================================================
# Random pixels use `colour` = leverage; training uses filled shapes (21 circle
# / 24 triangle) with `fill` = round -> two independent scales, no ggnewscale.
fig2 <- ggplot() +
  geom_point(data = rand, aes(spec_PC02, spec_PC03, colour = leverage),
             size = 0.8, alpha = 0.75) +
  scale_colour_viridis_c(option = "magma", direction = -1, trans = "sqrt",
                         name = "Leverage\n(priority)") +
  geom_point(data = train, aes(spec_PC02, spec_PC03, fill = round,
                               shape = class_type),
             size = 2, colour = "grey15", stroke = 0.3, alpha = 0.95) +
  scale_fill_manual(values = round_cols, name = "Field round") +
  scale_shape_manual(values = c(Meadow = 21, Shrub = 24), name = "Type") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = "spec_PC02", y = "spec_PC03",
       title = "Sampling-priority (leverage) surface with field rounds overlaid",
       subtitle = paste("Random basin pixels shaded by leverage (darker = higher priority).",
                        "Did the rounds reach the high-leverage gaps?")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave("docs/figures/obs_spectral_leverage.pdf", fig2, width = 9.5, height = 7.5,
       device = cairo_pdf)
cat("Wrote docs/figures/obs_spectral_leverage.pdf\n")

# ============================================================================
# FIGURE 3 — species composition PCoA (meadow), rounds + basin-prevalence
# ============================================================================
comp <- readRDS("data/derived/composition_species.rds")$hellinger
sp_cols <- setdiff(names(comp), c("site_number", "Year"))
pco <- stats::cmdscale(stats::dist(as.matrix(comp[, sp_cols])), k = 2, eig = TRUE)
ve  <- round(100 * pco$eig[1:2] / sum(pco$eig[pco$eig > 0]), 1)
comp_xy <- comp |>
  dplyr::mutate(Dim1 = pco$points[, 1], Dim2 = pco$points[, 2],
                round = round_lab(Year))

asg <- readRDS("data/derived/final_clusters_B.rds")$assignments |>
  dplyr::select(site_number, Year, final_label)
punch <- readr::read_csv("data/derived/punch_list.csv", show_col_types = FALSE) |>
  dplyr::select(final_label, predicted_n_pixels)
cent <- comp_xy |>
  dplyr::inner_join(asg, by = c("site_number", "Year")) |>
  dplyr::group_by(final_label) |>
  dplyr::summarise(Dim1 = mean(Dim1), Dim2 = mean(Dim2), .groups = "drop") |>
  dplyr::left_join(punch, by = "final_label") |>
  dplyr::mutate(predicted_n_pixels = tidyr::replace_na(predicted_n_pixels, 0))

fig3 <- ggplot() +
  geom_point(data = comp_xy, aes(Dim1, Dim2, colour = round),
             size = 1.6, alpha = 0.8) +
  stat_ellipse(data = comp_xy, aes(Dim1, Dim2, colour = round),
               type = "norm", linewidth = 0.7) +
  geom_point(data = cent, aes(Dim1, Dim2, size = predicted_n_pixels),
             shape = 21, fill = NA, colour = "grey25", stroke = 0.6) +
  # Label only the highest-basin-area classes (the sampling targets) to declutter.
  geom_text(data = dplyr::slice_max(cent, predicted_n_pixels, n = 12),
            aes(Dim1, Dim2, label = final_label),
            size = 2.7, fontface = "bold", colour = "grey15",
            check_overlap = TRUE, vjust = -0.8) +
  scale_colour_manual(values = round_cols, name = "Field round") +
  scale_size_area(max_size = 14, name = "Basin pixels\npredicted (area)") +
  labs(x = sprintf("PCoA 1 (%.1f%%)", ve[1]),
       y = sprintf("PCoA 2 (%.1f%%)", ve[2]),
       title = "Meadow species-composition space (Hellinger PCoA) by field round",
       subtitle = paste("Each point a plot; ellipses = per-round 95% coverage. Open circles = class",
                        "centroids sized by basin area predicted to that class (the sampling target).")) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave("docs/figures/obs_species_space.pdf", fig3, width = 10.5, height = 8.5,
       device = cairo_pdf)
cat("Wrote docs/figures/obs_species_space.pdf\n")
