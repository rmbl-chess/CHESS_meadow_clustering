# 12_diagnostics.R — meadow-phase diagnostic suite.
#
# Consolidates seven earlier per-purpose scripts (12 visualizations,
# 14 env, 15 K-sweep, 16 brightness, 18 year-effect, 21 inference QC,
# 22 size distribution) into one file so the meadow phase has fewer
# moving parts. Each section is independent; common data loads are
# shared at the top.
#
# Outputs (output/figures/, gitignored — move to docs/figures/ to commit):
#   spectra_per_spec_cluster.png    mean reflectance per spec cluster
#   spectra_per_final_label.png     sub-cluster spectra by parent
#   confusion_matrix.png            row-normalized RF CV confusion
#   dendrogram.png                  Ward dendrogram with k=8 cut
#   composition_profile.png         top species per final label
#   env_per_cluster.png             snow-free DOY violins by cluster
#   pc1_per_cluster.png             PC1 (brightness) by spec cluster
#   pc2_per_cluster.png             PC2 (greenness) by spec cluster
#   size_distribution.png           cluster sizes ordered by DOY
#   year_effect_pcs.png             2018-2025 systematic shifts
#   year_effect_pc_loadings.png     wavelength loadings for top PCs
# Plus per-section stdout summaries and a kspec_sweep.rds artifact.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  if (!requireNamespace("diptest", quietly = TRUE)) {
    renv::install("diptest")
  }
  library(diptest)
})

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# --- Common data loads ----------------------------------------------------
fc        <- readRDS("data/derived/final_clusters_B.rds")
sc        <- readRDS("data/derived/spectral_clusters.rds")
vs        <- readRDS("data/derived/veg_spectra.rds")
cs        <- readRDS("data/derived/composition_species.rds")
env       <- readRDS("data/derived/environment.rds")
spec_meta <- readRDS("data/derived/spectral_features.rds")
spec_feat <- spec_meta$features

variant_used <- "variant_G"     # active variant matches 10_cluster_spectra.R
spec_variant <- sc[[variant_used]]

joined <- vs$joined
wl     <- vs$wavelengths

# Consistent colour palette across all spec-cluster figures.
spec_levels  <- sort(unique(fc$assignments$spec_cluster))
spec_palette <- setNames(scales::hue_pal()(length(spec_levels)), spec_levels)
spec_label_lookup <- fc$spec_summary |>
  dplyr::transmute(
    spec_cluster,
    legend = sprintf("%s — %s (n=%d)", spec_cluster,
                     stringr::str_replace_all(indicator_species, "_", " "),
                     n_sites)
  ) |>
  tibble::deframe()

water_bands <- tibble::tibble(
  xmin = c(1340, 1800, 2400),
  xmax = c(1450, 1950, 2510)
)

eta_sq <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- factor(g[ok])
  grand <- mean(x)
  sum(table(g) * (tapply(x, g, mean) - grand)^2) / sum((x - grand)^2)
}

# ============================================================================
# SECTION A. Cluster output visualizations
#   (mean spectra, sub-cluster spectra, confusion matrix, dendrogram,
#    composition profile)
# ============================================================================

rfl_cols <- grep("^rfl_band_", names(joined), value = TRUE)
band_num_lookup <- tibble::tibble(
  band     = rfl_cols,
  band_num = as.integer(stringr::str_extract(rfl_cols, "\\d+$"))
) |>
  dplyr::left_join(wl, by = c("band_num" = "band_number"))

joined_lab <- joined |>
  dplyr::inner_join(
    fc$assignments |> dplyr::select(site_number, Year,
                                     spec_cluster, final_label),
    by = c("site_number", "Year")
  )

long_spec <- joined_lab |>
  dplyr::select(spec_cluster, final_label, dplyr::all_of(rfl_cols)) |>
  tidyr::pivot_longer(dplyr::all_of(rfl_cols),
                      names_to = "band", values_to = "rfl") |>
  dplyr::left_join(band_num_lookup |>
                     dplyr::select(band, center_wavelength_nm),
                   by = "band")

mean_per_spec <- long_spec |>
  dplyr::group_by(spec_cluster, center_wavelength_nm) |>
  dplyr::summarise(mean_rfl = mean(rfl, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(legend = spec_label_lookup[spec_cluster])

p_spec <- ggplot(mean_per_spec,
                 aes(x = center_wavelength_nm, y = mean_rfl,
                     colour = spec_cluster, group = spec_cluster)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = spec_palette, labels = spec_label_lookup,
                      name = "Spec cluster — indicator (n)") +
  labs(x = "Wavelength (nm)", y = "Brightness-normalized reflectance",
       title = sprintf("Per-spec-cluster mean spectrum (%s, k=8)",
                       variant_used),
       subtitle = "Grey bands = water absorption (masked in features)") +
  theme_minimal(base_size = 11)
ggsave("output/figures/spectra_per_spec_cluster.png", p_spec,
       width = 10, height = 6, dpi = 150)

mean_per_final <- long_spec |>
  dplyr::group_by(spec_cluster, final_label, center_wavelength_nm) |>
  dplyr::summarise(mean_rfl = mean(rfl, na.rm = TRUE), .groups = "drop")

p_subspec <- ggplot(mean_per_final,
                    aes(x = center_wavelength_nm, y = mean_rfl,
                        colour = final_label, group = final_label)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ spec_cluster, ncol = 4) +
  labs(x = "Wavelength (nm)", y = "Brightness-normalized reflectance",
       title = "Sub-cluster spectra within each spec cluster") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
ggsave("output/figures/spectra_per_final_label.png", p_subspec,
       width = 12, height = 7, dpi = 150)

# Confusion matrix
cm <- fc$final_eval$confusion
cm_long <- tibble::as_tibble(as.table(cm))
names(cm_long) <- c("truth", "pred", "n")
cm_long <- cm_long |>
  dplyr::group_by(truth) |>
  dplyr::mutate(prop = if (sum(n) > 0) n / sum(n) else 0) |>
  dplyr::ungroup()
final_order <- fc$final_summary |>
  dplyr::arrange(spec_cluster, final_label) |>
  dplyr::pull(final_label)
cm_long <- cm_long |>
  dplyr::mutate(truth = factor(truth, levels = final_order),
                pred  = factor(pred,  levels = final_order))

p_cm <- ggplot(cm_long, aes(x = pred, y = truth, fill = prop)) +
  geom_tile() +
  geom_text(aes(label = ifelse(prop >= 0.05, sprintf("%.2f", prop), "")),
            size = 2.5, colour = "white") +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1)) +
  scale_y_discrete(limits = rev) +
  labs(x = "Predicted label", y = "True label",
       title = "RF CV confusion matrix (row-normalized)",
       fill = "P(pred | truth)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
ggsave("output/figures/confusion_matrix.png", p_cm,
       width = 9, height = 8, dpi = 150)

# Ward dendrogram
png("output/figures/dendrogram.png", width = 1400, height = 700, res = 110)
op <- par(mar = c(2, 4, 3, 1))
hc <- spec_variant$hclust
plot(hc, labels = FALSE,
     main = sprintf("Spectral Ward hierarchical (%s), cut at k=8",
                    paste(spec_variant$pc_cols[c(1, length(spec_variant$pc_cols))],
                          collapse = "..")),
     sub = "", xlab = "", ylab = "Ward height")
rect.hclust(hc, k = 8, border = spec_palette[spec_levels])
par(op)
dev.off()

# Composition profile (top genera per final label)
nonsp_cover <- c("Other_Forb_cover", "Other_Graminoid_cover", "NPV_cover",
                 "Bare_cover", "Other_Moss_Lichen_cover",
                 "Other_Deciduous_Shrub_cover")
hell_named <- cs$hellinger |> dplyr::select(-dplyr::any_of(nonsp_cover))
hell_long_named <- hell_named |>
  tidyr::pivot_longer(-c(site_number, Year),
                      names_to = "feature", values_to = "h") |>
  dplyr::mutate(feature = stringr::str_replace(feature, "_cover$", ""),
                feature = stringr::str_replace_all(feature, "_", " ")) |>
  dplyr::inner_join(
    fc$assignments |> dplyr::select(site_number, Year,
                                     spec_cluster, final_label),
    by = c("site_number", "Year")
  )

top_species_per_label <- hell_long_named |>
  dplyr::group_by(spec_cluster, final_label, feature) |>
  dplyr::summarise(mean_h = mean(h), .groups = "drop") |>
  dplyr::group_by(spec_cluster, final_label) |>
  dplyr::slice_max(mean_h, n = 6, with_ties = FALSE) |>
  dplyr::ungroup()

p_comp <- ggplot(top_species_per_label,
                 aes(x = reorder(feature, mean_h), y = mean_h,
                     fill = spec_cluster)) +
  geom_col() +
  scale_fill_manual(values = spec_palette, guide = "none") +
  facet_wrap(~ final_label, scales = "free_y", ncol = 4) +
  coord_flip() +
  labs(x = NULL, y = "Mean Hellinger value (named species only)",
       title = "Top 6 named species per final label") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold", size = 9),
        axis.text.y = element_text(size = 7))
ggsave("output/figures/composition_profile.png", p_comp,
       width = 14, height = 10, dpi = 150)

cat("[A] Wrote spectra_per_spec_cluster, spectra_per_final_label,",
    " confusion_matrix, dendrogram, composition_profile\n")

# ============================================================================
# SECTION B. Environment (snow-free DOY) per spec cluster
# ============================================================================
asg_env <- fc$assignments |>
  dplyr::inner_join(env, by = c("site_number", "Year"))
cat(sprintf("\n[B] %d sites with both cluster + env\n", nrow(asg_env)))

spec_profile <- asg_env |>
  dplyr::group_by(spec_cluster) |>
  dplyr::summarise(
    n      = dplyr::n(),
    mean   = mean(snow_free_doy, na.rm = TRUE),
    sd     = sd(snow_free_doy, na.rm = TRUE),
    median = median(snow_free_doy, na.rm = TRUE),
    q25    = quantile(snow_free_doy, .25, na.rm = TRUE),
    q75    = quantile(snow_free_doy, .75, na.rm = TRUE),
    min    = min(snow_free_doy, na.rm = TRUE),
    max    = max(snow_free_doy, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::left_join(fc$spec_summary |>
                     dplyr::select(spec_cluster, indicator_species),
                   by = "spec_cluster") |>
  dplyr::arrange(mean)

eta_spec  <- eta_sq(asg_env$snow_free_doy, asg_env$spec_cluster)
eta_final <- eta_sq(asg_env$snow_free_doy, asg_env$final_label)
overall_sd <- sd(asg_env$snow_free_doy, na.rm = TRUE)

cat("\n=== Snow-free DOY per spectral cluster (ordered by mean) ===\n")
print(spec_profile, n = Inf, width = Inf)
cat(sprintf("\nOverall snow_free_doy: sd=%.1f days\n", overall_sd))
cat(sprintf("Eta² (spec_cluster):   %.3f\n", eta_spec))
cat(sprintf("Eta² (final_label):    %.3f\n", eta_final))

spec_order <- spec_profile$spec_cluster
asg_env_plot <- asg_env |>
  dplyr::mutate(
    spec_cluster = factor(spec_cluster, levels = spec_order),
    spec_label = sprintf("%s — %s",
                         spec_cluster,
                         stringr::str_replace_all(
                           fc$spec_summary$indicator_species[match(
                             as.character(spec_cluster),
                             fc$spec_summary$spec_cluster)],
                           "_", " "))
  ) |>
  dplyr::mutate(spec_label = factor(spec_label,
                                    levels = unique(spec_label[order(spec_cluster)])))

p_env <- ggplot(asg_env_plot, aes(x = spec_label, y = snow_free_doy)) +
  geom_violin(aes(fill = spec_label), alpha = 0.4, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.4, size = 0.7) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = "Snow-free DOY (1993–2022 mean)",
       title = "Snow-free date climatology per spec cluster",
       subtitle = sprintf("Eta²(spec) = %.3f   Eta²(final) = %.3f   SD = %.1f d",
                          eta_spec, eta_final, overall_sd)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("output/figures/env_per_cluster.png", p_env,
       width = 11, height = 6, dpi = 150)
cat("[B] Wrote env_per_cluster.png\n")

# ============================================================================
# SECTION C. PC1 / PC2 (brightness / greenness) per spec cluster
# ============================================================================
asg_pcs <- fc$assignments |>
  dplyr::inner_join(
    spec_feat |> dplyr::select(site_number, Year, spec_PC01, spec_PC02),
    by = c("site_number", "Year")
  )
eta_pc1 <- eta_sq(asg_pcs$spec_PC01, asg_pcs$spec_cluster)
eta_pc2 <- eta_sq(asg_pcs$spec_PC02, asg_pcs$spec_cluster)

per_cl_pcs <- asg_pcs |>
  dplyr::group_by(spec_cluster) |>
  dplyr::summarise(
    n        = dplyr::n(),
    pc1_mean = mean(spec_PC01), pc1_sd = sd(spec_PC01),
    pc2_mean = mean(spec_PC02), pc2_sd = sd(spec_PC02),
    .groups  = "drop"
  ) |>
  dplyr::left_join(fc$spec_summary |>
                     dplyr::select(spec_cluster, indicator_species),
                   by = "spec_cluster") |>
  dplyr::arrange(pc1_mean)

cat(sprintf("\n[C] Eta²(PC1) = %.3f   Eta²(PC2) = %.3f   (cf. Eta²(DOY) = %.3f)\n",
            eta_pc1, eta_pc2, eta_spec))
print(per_cl_pcs, n = Inf, width = Inf)

cluster_order_pcs <- per_cl_pcs$spec_cluster
asg_pcs_plot <- asg_pcs |>
  dplyr::mutate(
    spec_cluster = factor(spec_cluster, levels = cluster_order_pcs),
    label = sprintf("%s — %s",
                    spec_cluster,
                    stringr::str_replace_all(
                      fc$spec_summary$indicator_species[match(
                        as.character(spec_cluster),
                        fc$spec_summary$spec_cluster)],
                      "_", " "))
  ) |>
  dplyr::mutate(label = factor(label,
                               levels = unique(label[order(spec_cluster)])))

make_pc_plot <- function(pc_col, eta, title, subtitle) {
  ggplot(asg_pcs_plot, aes(x = label, y = .data[[pc_col]])) +
    geom_violin(aes(fill = label), alpha = 0.4, scale = "width",
                draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_jitter(width = 0.18, height = 0, alpha = 0.4, size = 0.7) +
    scale_fill_brewer(palette = "Set3", guide = "none") +
    labs(x = NULL, y = title, subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
ggsave("output/figures/pc1_per_cluster.png",
       make_pc_plot("spec_PC01", eta_pc1,
                    "PC1 (brightness, ~95% spectral variance)",
                    sprintf("Eta²(PC1) = %.3f", eta_pc1)),
       width = 12, height = 5, dpi = 150)
ggsave("output/figures/pc2_per_cluster.png",
       make_pc_plot("spec_PC02", eta_pc2,
                    "PC2 (greenness / red edge, ~3% spectral variance)",
                    sprintf("Eta²(PC2) = %.3f", eta_pc2)),
       width = 12, height = 5, dpi = 150)
cat("[C] Wrote pc1_per_cluster.png + pc2_per_cluster.png\n")

# ============================================================================
# SECTION D. 2018 inference-quality report
# ============================================================================
clusters_with_inf <- asg_env |>
  dplyr::filter(source == "inferred_2018") |>
  dplyr::distinct(final_label) |> dplyr::pull(final_label)

per_cluster_inf <- purrr::map_dfr(clusters_with_inf, function(cl) {
  rows <- asg_env |> dplyr::filter(final_label == cl)
  a <- rows |> dplyr::filter(source == "clustered")
  i <- rows |> dplyr::filter(source == "inferred_2018")
  tibble::tibble(
    label = cl, n_anch = nrow(a), n_inf = nrow(i),
    anch_doy = round(mean(a$snow_free_doy), 1),
    inf_doy  = round(mean(i$snow_free_doy), 1),
    doy_shift   = round(mean(i$snow_free_doy) - mean(a$snow_free_doy), 1),
    med_doy_diff = round(stats::median(abs(i$inference_doy_diff_days)), 1),
    p90_doy_diff = round(stats::quantile(abs(i$inference_doy_diff_days), 0.9), 1),
    med_hell     = round(stats::median(i$inference_hell_distance), 2),
    p90_hell     = round(stats::quantile(i$inference_hell_distance, 0.9), 2)
  )
}) |> dplyr::arrange(dplyr::desc(p90_doy_diff))

cat("\n[D] === Per-cluster inference quality (sorted by p90 DOY mismatch) ===\n")
print(per_cluster_inf, n = Inf, width = Inf)

conf_tally <- asg_env |>
  dplyr::filter(source == "inferred_2018") |>
  dplyr::count(inference_confidence) |>
  dplyr::mutate(pct = round(100 * n / sum(n), 1))
cat("\n[D] === Inference confidence breakdown ===\n")
print(conf_tally)

cat("\n[D] === Top 10 worst-fit inferred sites ===\n")
worst <- asg_env |>
  dplyr::filter(source == "inferred_2018") |>
  dplyr::arrange(dplyr::desc(inference_distance)) |>
  dplyr::slice_head(n = 10) |>
  dplyr::select(site_number, final_label, snow_free_doy,
                inference_doy_diff_days, inference_hell_distance,
                inference_distance, inference_gap, inference_confidence)
print(worst, n = Inf, width = Inf)

# ============================================================================
# SECTION E. Cluster size distribution by tier + source
# ============================================================================
desc <- readr::read_csv("data/small_reference/label_community_names.csv",
                        show_col_types = FALSE)

sizes <- fc$assignments |>
  dplyr::count(final_label, source) |>
  tidyr::pivot_wider(names_from = source, values_from = n, values_fill = 0L) |>
  # Ensure both columns exist even if a source is empty after the year-
  # correction (e.g., no `inferred_2018` rows when all 2018 sites cluster
  # directly).
  (\(df) {
    if (!"clustered"     %in% names(df)) df$clustered     <- 0L
    if (!"inferred_2018" %in% names(df)) df$inferred_2018 <- 0L
    df
  })() |>
  dplyr::rename(n_anchor = clustered, n_inferred = inferred_2018) |>
  dplyr::mutate(n_total = n_anchor + n_inferred)

size_joined <- desc |>
  dplyr::select(final_label, recall, tier, snow_free_doy_mean) |>
  dplyr::left_join(sizes, by = "final_label") |>
  dplyr::mutate(tier = factor(tier, levels = c("strong", "marginal", "weak")))

plot_size <- size_joined |>
  dplyr::arrange(snow_free_doy_mean) |>
  dplyr::mutate(final_label = factor(final_label, levels = final_label)) |>
  tidyr::pivot_longer(c(n_anchor, n_inferred),
                      names_to = "source", values_to = "n") |>
  dplyr::mutate(source = dplyr::recode(source,
                                       n_anchor = "Clustered (2025)",
                                       n_inferred = "Inferred (2018)"))

p_size <- ggplot(plot_size,
                 aes(x = final_label, y = n, fill = source)) +
  geom_col(width = 0.75) +
  geom_text(
    data = size_joined |>
      dplyr::arrange(snow_free_doy_mean) |>
      dplyr::mutate(final_label = factor(final_label, levels = final_label)),
    aes(y = n_total + 1.5, x = final_label,
        label = sprintf("%.2f", recall), fill = NULL),
    size = 2.9, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("Clustered (2025)" = "#3a8fb7",
                               "Inferred (2018)"  = "#dca066"),
                    name = "Source") +
  labs(x = NULL, y = "Site count",
       title = "Cluster size distribution, ordered by snow-free DOY",
       subtitle = "Number above bar = CV recall (clustered-2025 only)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "top",
        panel.grid.major.x = element_blank())
ggsave("output/figures/size_distribution.png", p_size,
       width = 11, height = 6, dpi = 150)

cat("\n[E] === Size summary by tier ===\n")
print(size_joined |>
        dplyr::group_by(tier) |>
        dplyr::summarise(
          n_clusters     = dplyr::n(),
          anchor_total   = sum(n_anchor),
          inferred_total = sum(n_inferred),
          grand_total    = sum(n_total),
          min_size       = min(n_total),
          median_size    = stats::median(n_total),
          max_size       = max(n_total),
          .groups        = "drop"
        ) |>
        dplyr::arrange(tier))
cat("[E] Wrote size_distribution.png\n")

# ============================================================================
# SECTION F. K-sweep on variant_G (env coherence, dip stat, RF accuracy)
# ============================================================================
spec_cols <- grep("^(spec_PC|ndvi|ndwi|pri|red_edge|cai|ndli)",
                  names(spec_feat), value = TRUE)
# Match 10's classifier: drop year-shift features PC3 and PRI.
spec_cols <- setdiff(spec_cols, c("spec_PC03", "pri"))
joined_for_sweep <- spec_variant$assignments |>
  dplyr::inner_join(spec_feat, by = c("site_number", "Year")) |>
  dplyr::inner_join(env,       by = c("site_number", "Year"))
X_sweep <- as.matrix(joined_for_sweep[, c(spec_cols, "snow_free_doy")])

eval_rf_cv <- function(labels, X, n_folds = 5, seed = 42, n_trees = 500) {
  y <- factor(labels); n <- length(y)
  set.seed(seed)
  fold <- integer(n)
  for (lvl in levels(y)) {
    idx <- which(y == lvl)
    fold[idx] <- ((sample(seq_along(idx)) - 1L) %% n_folds) + 1L
  }
  preds <- factor(rep(NA_character_, n), levels = levels(y))
  for (f in seq_len(n_folds)) {
    tr <- which(fold != f); te <- which(fold == f)
    bad <- names(which(table(y[tr]) < 2))
    keep_tr <- tr[!(y[tr] %in% bad)]
    df_tr <- data.frame(X[keep_tr, , drop = FALSE], .y = y[keep_tr])
    fit <- ranger::ranger(.y ~ ., data = df_tr, num.trees = n_trees,
                          classification = TRUE, seed = seed, verbose = FALSE)
    preds[te] <- predict(fit, data.frame(X[te, , drop = FALSE]))$predictions
  }
  ok <- !is.na(preds); truth <- y[ok]; pred <- preds[ok]
  cm <- table(truth = truth, pred = pred)
  list(accuracy = mean(truth == pred),
       recall   = diag(cm) / rowSums(cm))
}

sweep_results <- purrr::map_dfr(sc$ks, function(k) {
  k_col  <- sprintf("k%02d", k)
  labels <- joined_for_sweep[[k_col]]
  doy    <- joined_for_sweep$snow_free_doy
  per_cl <- tibble::tibble(cluster = labels, doy = doy) |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      n   = dplyr::n(),
      sd  = sd(doy, na.rm = TRUE),
      iqr = IQR(doy, na.rm = TRUE),
      dip = if (sum(!is.na(doy)) >= 4)
              suppressWarnings(diptest::dip.test(doy))$statistic
            else NA_real_,
      .groups = "drop"
    )
  rf  <- eval_rf_cv(labels, X_sweep)
  rec <- as.numeric(rf$recall)
  tibble::tibble(
    k                = k,
    eta_sq_doy       = eta_sq(doy, labels),
    mean_within_iqr  = mean(per_cl$iqr, na.rm = TRUE),
    n_bimodal        = sum(per_cl$dip > 0.08, na.rm = TRUE),
    cv_accuracy      = rf$accuracy,
    n_strong         = sum(rec >= 0.80,                na.rm = TRUE),
    n_marginal       = sum(rec >= 0.50 & rec < 0.80,   na.rm = TRUE),
    n_weak           = sum(rec <  0.50,                na.rm = TRUE),
    n_tiny           = sum(per_cl$n  < 6),
    smallest_cluster = min(per_cl$n)
  )
})

cat("\n[F] === Variant G k_spec sweep ===\n")
cat("(higher eta_sq, lower iqr/dip; cv_accuracy tracks the chosen K)\n")
print(sweep_results |>
        dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))),
      n = Inf, width = Inf)
saveRDS(sweep_results, "data/derived/kspec_sweep.rds")
cat("[F] Wrote data/derived/kspec_sweep.rds\n")

# ============================================================================
# SECTION G. Year-effect on spectral PCs (2018 vs 2025 systematic shifts)
# ============================================================================
csh <- cs$hellinger
hell_cols <- setdiff(names(csh), c("site_number", "Year"))
H <- as.matrix(csh[, hell_cols])

idx_18 <- csh$Year == 2018L
idx_25 <- csh$Year == 2025L
H_18 <- H[idx_18, , drop = FALSE]; H_25 <- H[idx_25, , drop = FALSE]
keys_18 <- csh |> dplyr::filter(Year == 2018L) |>
  dplyr::select(site_number, Year)
keys_25 <- csh |> dplyr::filter(Year == 2025L) |>
  dplyr::select(site_number, Year)

sq_norms_18 <- rowSums(H_18^2); sq_norms_25 <- rowSums(H_25^2)
cross_d2 <- outer(sq_norms_18, sq_norms_25, "+") - 2 * H_18 %*% t(H_25)
cross_d2[cross_d2 < 0] <- 0
cross_d <- sqrt(cross_d2)
best_idx <- apply(cross_d, 1, which.min)
best_d   <- cross_d[cbind(seq_along(best_idx), best_idx)]

matches <- tibble::tibble(
  site_18 = keys_18$site_number,
  site_25 = keys_25$site_number[best_idx],
  hell_dist = best_d
)
threshold <- stats::median(matches$hell_dist)
matches_close <- matches |> dplyr::filter(hell_dist < threshold)
cat(sprintf("\n[G] All pairs: %d  Close pairs (Hell dist < %.3f, median): %d\n",
            nrow(matches), threshold, nrow(matches_close)))

spec_18 <- spec_feat |> dplyr::filter(Year == 2018L) |>
  dplyr::select(site_number, dplyr::starts_with("spec_PC"),
                ndvi, ndwi, pri, red_edge_slope, cai, ndli)
spec_25 <- spec_feat |> dplyr::filter(Year == 2025L) |>
  dplyr::select(site_number, dplyr::starts_with("spec_PC"),
                ndvi, ndwi, pri, red_edge_slope, cai, ndli)
pairs <- matches_close |>
  dplyr::inner_join(spec_18, by = c("site_18" = "site_number")) |>
  dplyr::inner_join(spec_25, by = c("site_25" = "site_number"),
                    suffix = c("_18", "_25"))

all_features <- c(
  grep("^spec_PC", names(spec_feat), value = TRUE),
  intersect(c("ndvi", "ndwi", "pri", "red_edge_slope", "cai", "ndli"),
            names(spec_feat))
)
feature_cols <- all_features[
  paste0(all_features, "_18") %in% names(pairs) &
  paste0(all_features, "_25") %in% names(pairs)
]

diffs <- purrr::map_dfc(feature_cols, function(f) {
  v18 <- pairs[[paste0(f, "_18")]]
  v25 <- pairs[[paste0(f, "_25")]]
  tibble::tibble(!!f := v25 - v18)
})
feat_sds <- vapply(feature_cols,
                   function(f) sd(spec_feat[[f]], na.rm = TRUE),
                   numeric(1))
per_feat <- purrr::map_dfr(feature_cols, function(f) {
  d <- diffs[[f]]
  if (all(is.na(d))) return(NULL)
  tt <- stats::t.test(d)
  tibble::tibble(
    feature       = f,
    mean_diff     = mean(d, na.rm = TRUE),
    sd_diff       = sd(d, na.rm = TRUE),
    t_stat        = unname(tt$statistic),
    p_value       = tt$p.value,
    feature_sd    = feat_sds[f],
    abs_effect_sd = abs(mean(d, na.rm = TRUE)) / feat_sds[f]
  )
}) |> dplyr::arrange(dplyr::desc(abs_effect_sd))

cat("\n[G] === Per-feature systematic shift (2025 - 2018), |effect/SD| sorted ===\n")
print(per_feat |>
        dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3))),
      n = Inf, width = Inf)

top_features <- per_feat |> dplyr::slice_head(n = 12) |> dplyr::pull(feature)
plot_df_year <- diffs |>
  dplyr::select(dplyr::all_of(top_features)) |>
  tidyr::pivot_longer(dplyr::everything(),
                      names_to = "feature", values_to = "diff") |>
  dplyr::mutate(feature = factor(feature, levels = top_features))

p_year <- ggplot(plot_df_year, aes(x = feature, y = diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_violin(fill = "steelblue", alpha = 0.4, scale = "width",
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.3, size = 0.7) +
  labs(x = NULL, y = "Value (2025) − Value (2018)",
       title = "Year effect on spectral features (matched pairs by composition)",
       subtitle = sprintf("%d matched pairs (Hellinger dist < median).",
                          nrow(pairs))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("output/figures/year_effect_pcs.png", p_year,
       width = 12, height = 6, dpi = 150)

# Wavelength loadings for top-shifting PCs
top_pcs <- per_feat |>
  dplyr::filter(grepl("^spec_PC", feature)) |>
  dplyr::slice_head(n = 6) |> dplyr::pull(feature)
pca <- spec_meta$pca; wl_nm <- spec_meta$keep_wl
loadings_long <- tibble::tibble(wavelength_nm = wl_nm)
for (pc_name in top_pcs) {
  pc_idx <- as.integer(stringr::str_extract(pc_name, "\\d+"))
  loadings_long[[pc_name]] <- pca$rotation[, pc_idx]
}
loadings_long <- loadings_long |>
  tidyr::pivot_longer(-wavelength_nm, names_to = "PC", values_to = "loading") |>
  dplyr::mutate(PC = factor(PC, levels = top_pcs))

p_load <- ggplot(loadings_long,
                 aes(x = wavelength_nm, y = loading, colour = PC)) +
  geom_rect(data = water_bands, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ PC, ncol = 2, scales = "free_y") +
  labs(x = "Wavelength (nm)", y = "PCA loading",
       title = "Wavelength loadings for top year-shifting PCs") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))
ggsave("output/figures/year_effect_pc_loadings.png", p_load,
       width = 12, height = 7, dpi = 150)

saveRDS(list(matches = matches, matches_close = matches_close,
             diffs = diffs, per_feature = per_feat,
             threshold = threshold),
        "data/derived/year_effect.rds")
cat("[G] Wrote year_effect_pcs.png, year_effect_pc_loadings.png, year_effect.rds\n")

cat("\n=== Diagnostic suite complete ===\n")
