# 10_subcluster_composition.R — Architecture B step 2: within each spectral
# cluster (variant F, k=12 = PCs 2-12 z-scaled + DOY z-scaled), sub-cluster sites by
# Hellinger composition AT SPECIES LEVEL. Coherent spec clusters stay as
# single training classes; heterogeneous ones split into 2-3 species-level
# sub-types by Ward on species Hellinger.
#
# Note: sub-clusters within a spectral cluster share the same spectral
# signature, so the downstream AOP classifier will not be able to fully
# distinguish them. The final RF eval here is diagnostic — it shows exactly
# which sub-classes are spectrally confused.
#
# Inputs:
#   data/derived/spectral_clusters.rds     (variant_C$assignments at k=8)
#   data/derived/composition_species.rds   (species-level Hellinger)
#   data/derived/spectral_features.rds     (RF input features)
# Outputs:
#   data/derived/final_clusters_B.rds      same structure as before, with
#                                          species-level indicators.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
})

sc           <- readRDS("data/derived/spectral_clusters.rds")
comp_species <- readRDS("data/derived/composition_species.rds")
spec_feat    <- readRDS("data/derived/spectral_features.rds")$features
env          <- readRDS("data/derived/environment.rds")

primary_k        <- 26
primary_variant  <- "variant_G"  # PCs 2,4-12 + DOY z-scaled (drop PC1 brightness + PC3 year-shift)
# do_subclustering = FALSE means spec_cluster is the training label and the
# sub_cluster column is kept only as ecological metadata. Sub-classes within a
# parent share the parent's spectral signature, so the AOP classifier can't
# tell them apart -- forcing them apart just creates classes with recall ~0.
# Setting this to TRUE restores the old behavior of splitting heterogeneous
# parents into sub-labels.
do_subclustering <- FALSE
k_col            <- sprintf("k%02d", primary_k)
spec_summary    <- sc[[primary_variant]]$characterizations[[k_col]]

asg <- sc[[primary_variant]]$assignments |>
  dplyr::select(site_number, Year, spec_cluster = dplyr::all_of(k_col)) |>
  mutate(spec_cluster = sprintf("S%02d", spec_cluster))

# --- Identify heterogeneous spec clusters ----------------------------------
dom_threshold <- 0.40
het_threshold <- 0.60   # relaxed slightly so cluster 7 (het 0.58) stays coherent

spec_summary <- spec_summary |>
  mutate(spec_cluster = sprintf("S%02d", cluster),
         coherent = dominance >= dom_threshold | heterogeneity <= het_threshold)

cat("\nCoherence classification (Variant F, species-level):\n")
print(spec_summary |> dplyr::select(spec_cluster, n_sites, dominance,
                                    heterogeneity, indicator_species,
                                    indicator_genus, coherent),
      n = Inf, width = Inf)

incoherent <- spec_summary$spec_cluster[!spec_summary$coherent]

# --- Sub-cluster the incoherent ones by species-level Hellinger Ward -------
hell <- comp_species$hellinger
nonsp_cover <- c("Other_Forb_cover", "Other_Graminoid_cover", "NPV_cover",
                 "Bare_cover", "Other_Moss_Lichen_cover",
                 "Other_Deciduous_Shrub_cover")
named_cols <- setdiff(names(hell), c("site_number", "Year", nonsp_cover))

asg$sub_cluster <- NA_character_
sub_details <- list()

for (sc_id in incoherent) {
  sites_in <- asg |> filter(spec_cluster == sc_id) |>
    dplyr::select(site_number, Year)
  hell_subset <- hell |>
    inner_join(sites_in, by = c("site_number", "Year")) |>
    dplyr::select(site_number, Year, dplyr::all_of(named_cols))
  if (nrow(hell_subset) < 6) next
  n_sub <- if (nrow(hell_subset) >= 60) 3L else 2L

  d  <- dist(as.matrix(hell_subset[, named_cols]), method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  cuts <- cutree(hc, k = n_sub)
  sub_tags <- letters[cuts]

  # Assign back to asg
  for (i in seq_along(sub_tags)) {
    mask <- asg$site_number == hell_subset$site_number[i] &
            asg$Year        == hell_subset$Year[i]
    asg$sub_cluster[mask] <- sub_tags[i]
  }
  sub_details[[sc_id]] <- list(hclust = hc, k = n_sub, n = nrow(hell_subset))
}

# Orphans: sites in heterogeneous spec clusters whose composition was too
# sparse for sub-clustering (only rare/zero named genera). Absorb them into
# the modal (most populous) sub-cluster of their parent — keeps them in the
# training set with a reasonable spec-cluster-consistent label rather than
# creating singleton classes.
modal_sub <- asg |>
  filter(!is.na(sub_cluster)) |>
  count(spec_cluster, sub_cluster) |>
  group_by(spec_cluster) |>
  slice_max(n, n = 1, with_ties = FALSE) |>
  ungroup() |>
  dplyr::select(spec_cluster, modal_sub = sub_cluster)

asg <- asg |>
  left_join(modal_sub, by = "spec_cluster") |>
  mutate(
    sub_cluster = if_else(
      is.na(sub_cluster) & spec_cluster %in% incoherent,
      modal_sub, sub_cluster
    )
  ) |>
  dplyr::select(-modal_sub) |>
  mutate(final_label = if (do_subclustering) {
    if_else(is.na(sub_cluster),
            spec_cluster,
            paste0(spec_cluster, ".", sub_cluster))
  } else {
    spec_cluster
  },
  source = "clustered")  # both 2018 + 2025 cluster directly after the
                          # year-effect radiometric correction (see
                          # code/joint/13_year_effect_analysis.R).

# --- Hellinger + DOY nearest-cluster fallback for 2018 sites without ------
# spectra ------------------------------------------------------------------
# Once 2018 AOP spectra have the per-band radiometric correction applied
# in code/meadow/04_join_spectra.R, 2018 sites cluster directly alongside
# 2025 (variant G in 10_cluster_spectra.R is no longer 2025-only). This
# block now ONLY handles the rare 2018 sites that lack spectra or env
# coverage and therefore drop out of feat_env_G. Each such site is
# assigned to its closest direct-clustered cluster in a JOINT space:
#   - full-species Hellinger (NOT the rare-trimmed comp_species — that
#     would zero out alpine sites whose indicators are rare species like
#     Sibbaldia procumbens or Mertensia lanceolata, letting them match
#     any cluster with a sparse trimmed centroid at distance ~0)
#   - snow-free DOY, z-scaled and weighted by `doy_weight` (default 1.0)
#     so a 1-SD (~21 day) DOY mismatch contributes the same as a 1-unit
#     Hellinger mismatch.
doy_weight <- 1.0   # raise to be stricter about DOY agreement
cover_combined <- readRDS("data/derived/cover_combined.rds")
nonsp_cover <- c("Other_Forb_cover", "Other_Graminoid_cover", "NPV_cover",
                 "Bare_cover", "Other_Moss_Lichen_cover",
                 "Other_Deciduous_Shrub_cover")

# Build full-species Hellinger for ALL sites (no rare-trim).
cov_cols     <- grep("_cover$", names(cover_combined), value = TRUE)
species_cols <- setdiff(cov_cols, nonsp_cover)
cov_mat <- as.matrix(cover_combined[, species_cols])
totals  <- rowSums(cov_mat)
hell_mat <- sqrt(cov_mat / pmax(totals, 1e-9))
keep <- totals > 0   # drop sites with no named-species cover

hell_full <- dplyr::bind_cols(
  cover_combined[keep, c("site_number", "Year")],
  tibble::as_tibble(hell_mat[keep, , drop = FALSE])
) |>
  dplyr::inner_join(env, by = c("site_number", "Year"))

doy_sd_global <- sd(hell_full$snow_free_doy)

hell_2025_full <- hell_full |>
  dplyr::filter(Year == 2025L) |>
  dplyr::inner_join(asg |> dplyr::select(site_number, Year, spec_cluster),
                    by = c("site_number", "Year"))
clusters_present <- sort(unique(hell_2025_full$spec_cluster))

centroid_mat <- vapply(clusters_present, function(cl) {
  rows <- hell_2025_full |> dplyr::filter(spec_cluster == cl)
  colMeans(as.matrix(rows[, species_cols]))
}, numeric(length(species_cols)))
centroid_mat <- t(centroid_mat); rownames(centroid_mat) <- clusters_present

centroid_doy <- hell_2025_full |>
  dplyr::group_by(spec_cluster) |>
  dplyr::summarise(mean_doy = mean(snow_free_doy), .groups = "drop") |>
  dplyr::arrange(match(spec_cluster, clusters_present))
centroid_doy_vec <- centroid_doy$mean_doy
names(centroid_doy_vec) <- centroid_doy$spec_cluster

# Skip 2018 sites that already got a direct cluster assignment.
already_clustered_2018 <- asg |>
  dplyr::filter(Year == 2018L) |>
  dplyr::pull(site_number) |> unique()
hell_2018_full <- hell_full |>
  dplyr::filter(Year == 2018L, !site_number %in% already_clustered_2018)
cat(sprintf("Hellinger fallback: %d 2018 sites without direct cluster (of %d total 2018 in cover)\n",
            nrow(hell_2018_full),
            sum(hell_full$Year == 2018L)))
if (nrow(hell_2018_full) > 0) {
  H_2018   <- as.matrix(hell_2018_full[, species_cols])
  doy_2018 <- hell_2018_full$snow_free_doy

  # Hellinger Euclidean squared (n_2018 x n_clusters)
  sq_2018   <- rowSums(H_2018^2)
  sq_cent   <- rowSums(centroid_mat^2)
  hell_d2   <- outer(sq_2018, sq_cent, "+") - 2 * H_2018 %*% t(centroid_mat)
  hell_d2[hell_d2 < 0] <- 0
  hell_d    <- sqrt(hell_d2)

  # DOY z-difference squared
  doy_diff   <- outer(doy_2018, centroid_doy_vec, "-")
  doy_z_diff <- doy_diff / doy_sd_global
  doy_d2     <- (doy_z_diff * doy_weight)^2

  combined_d <- sqrt(hell_d2 + doy_d2)
  rownames(combined_d) <- paste(hell_2018_full$site_number,
                                hell_2018_full$Year, sep = "_")

  best_idx     <- apply(combined_d, 1, which.min)
  dist_to_best <- combined_d[cbind(seq_along(best_idx), best_idx)]
  sorted_d     <- t(apply(combined_d, 1, sort))
  dist_gap     <- sorted_d[, 2] - sorted_d[, 1]

  asg_2018 <- tibble::tibble(
    site_number             = hell_2018_full$site_number,
    Year                    = hell_2018_full$Year,
    spec_cluster            = clusters_present[best_idx],
    sub_cluster             = NA_character_,
    final_label             = clusters_present[best_idx],
    source                  = "inferred_2018",
    inference_distance      = dist_to_best,
    inference_gap           = dist_gap,
    inference_hell_distance = hell_d[cbind(seq_along(best_idx), best_idx)],
    inference_doy_diff_days = doy_diff[cbind(seq_along(best_idx), best_idx)]
  ) |>
    dplyr::mutate(inference_confidence = dplyr::case_when(
      inference_distance < 0.90 ~ "high",
      inference_distance < 1.05 ~ "medium",
      TRUE                      ~ "low"
    ))
  cat(sprintf("Inferred 2018 labels for %d fallback sites\n",
              nrow(asg_2018)))
} else {
  asg_2018 <- tibble::tibble(
    site_number = integer(0), Year = integer(0),
    spec_cluster = character(0), sub_cluster = character(0),
    final_label = character(0), source = character(0),
    inference_distance = double(0), inference_gap = double(0),
    inference_hell_distance = double(0),
    inference_doy_diff_days = double(0),
    inference_confidence    = character(0)
  )
  cat("All 2018 sites in cover already have direct cluster assignments; ",
      "no Hellinger fallback needed.\n", sep = "")
}
if (nrow(asg_2018) > 0) {
cat(sprintf("  combined inference distance: median=%.3f, IQR %.3f-%.3f\n",
              stats::median(asg_2018$inference_distance),
              stats::quantile(asg_2018$inference_distance, 0.25),
              stats::quantile(asg_2018$inference_distance, 0.75)))
  cat(sprintf("  Hellinger-only distance:     median=%.3f, IQR %.3f-%.3f\n",
              stats::median(asg_2018$inference_hell_distance),
              stats::quantile(asg_2018$inference_hell_distance, 0.25),
              stats::quantile(asg_2018$inference_hell_distance, 0.75)))
  cat(sprintf("  DOY mismatch (days):         median=%.1f, IQR %.1f-%.1f\n",
              stats::median(abs(asg_2018$inference_doy_diff_days)),
              stats::quantile(abs(asg_2018$inference_doy_diff_days), 0.25),
              stats::quantile(abs(asg_2018$inference_doy_diff_days), 0.75)))
}  # close `if (nrow(asg_2018) > 0)`

# Add the inference columns to asg (NA for clustered rows) and stack.
asg <- asg |>
  dplyr::mutate(inference_distance = NA_real_, inference_gap = NA_real_) |>
  dplyr::bind_rows(asg_2018)

# --- Monotypic-stand species overrides --------------------------------------
# Some species form dense monotypic stands worth mapping as their own
# classes. A site with >= monotypic_threshold cover of a candidate species
# gets its final_label overridden to a monotypic class code (Mxx_...).
# Each site can only trigger one override (cover sums to 100 per site).
monotypic_threshold <- 70
monotypic_table <- tibble::tribble(
  ~species_col,                 ~label,
  "Veratrum_tenuipetalum_cover",  "Veratrum tenuipetalum",
  "Ligusticum_porteri_cover",     "Ligusticum porteri",
  "Caltha_leptosepala_cover",     "Caltha leptosepala",
  "Corydalis_caseana_cover",      "Corydalis caseana",
  "Osmorhiza_occidentalis_cover", "Osmorhiza occidentalis"
)

# Pull just the candidate cover columns from the joined veg+spectra table.
cover_for_override <- readRDS("data/derived/cover_combined.rds") |>
  dplyr::select(site_number, Year,
                dplyr::any_of(monotypic_table$species_col))

asg <- asg |>
  dplyr::left_join(cover_for_override, by = c("site_number", "Year")) |>
  dplyr::mutate(monotypic_species = NA_character_)

n_overridden_per_species <- integer(nrow(monotypic_table))
names(n_overridden_per_species) <- monotypic_table$label
for (i in seq_len(nrow(monotypic_table))) {
  col <- monotypic_table$species_col[i]
  lbl <- monotypic_table$label[i]
  if (!col %in% names(asg)) next
  mask <- !is.na(asg[[col]]) & asg[[col]] >= monotypic_threshold
  n_overridden_per_species[i] <- sum(mask)
  asg$final_label[mask]       <- lbl
  asg$monotypic_species[mask] <- sub("_cover$", "", col)
  asg$spec_cluster[mask]      <- lbl   # also override spec_cluster for consistency
}
cat(sprintf("\nMonotypic overrides applied (threshold = %d%% cover):\n",
            monotypic_threshold))
print(tibble::tibble(label = names(n_overridden_per_species),
                     n_sites = n_overridden_per_species))

# Drop the temporary cover columns (we keep monotypic_species as the flag).
asg <- asg |> dplyr::select(-dplyr::any_of(monotypic_table$species_col))

# --- RF eval on the CLUSTERED 2025 set only (NOT the inferred 2018) -------
# Including inferred labels in CV would be circular -- they were assigned
# by composition similarity, not learned from spectra.
spec_cols <- grep("^(spec_PC|ndvi|ndwi|pri|red_edge|cai|ndli)",
                  names(spec_feat), value = TRUE)
# Full feature set is fine now -- 2025-only clustering means the cross-year
# shift in PC3 and PRI is no longer a concern.
asg_clustered <- asg |> dplyr::filter(source == "clustered")
joined <- asg_clustered |>
  inner_join(spec_feat, by = c("site_number", "Year")) |>
  inner_join(env,       by = c("site_number", "Year"))
rf_feature_cols <- c(spec_cols, "snow_free_doy")
X <- as.matrix(joined[, rf_feature_cols])

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
  recall <- diag(cm) / rowSums(cm)
  list(accuracy = mean(truth == pred), recall = recall, confusion = cm)
}

final_eval <- eval_rf_cv(joined$final_label, X)

# --- Per-final-label characterization ---------------------------------------
# Characterization uses CLUSTERED (2025) sites only: the cluster is defined
# by 2025 composition + spectra, so its description should reflect that.
# Inferred 2018 assignments are observational; including them here would
# mix what defines the cluster with what was assigned to it.
hell_long <- comp_species$hellinger |>
  pivot_longer(-c(site_number, Year), names_to = "feature", values_to = "h") |>
  mutate(feature = stringr::str_replace(feature, "_cover$", "")) |>
  inner_join(asg |> dplyr::filter(source == "clustered") |>
                    dplyr::select(site_number, Year, final_label),
             by = c("site_number", "Year"))

nonsp_short <- stringr::str_replace(nonsp_cover, "_cover$", "")

per_label <- hell_long |>
  group_by(final_label, feature) |>
  summarise(mean_h = mean(h), .groups = "drop") |>
  group_by(final_label) |>
  arrange(desc(mean_h), .by_group = TRUE) |>
  summarise(
    indicator_species = {
      sp <- feature[!feature %in% nonsp_short]
      if (length(sp) == 0) "(no named)" else sp[1]
    },
    top_features = paste(head(feature, 5), collapse = ", "),
    .groups = "drop"
  )

sizes <- asg |> count(final_label, name = "n")
year_b <- asg |> count(final_label, Year) |>
  pivot_wider(names_from = Year, values_from = n, names_prefix = "n_",
              values_fill = 0L)

final_summary <- sizes |>
  left_join(year_b, by = "final_label") |>
  mutate(spec_cluster = stringr::str_extract(final_label, "^S\\d+"),
         recall = as.numeric(final_eval$recall[final_label])) |>
  left_join(per_label, by = "final_label") |>
  arrange(spec_cluster, final_label)

cat(sprintf("\n=== FINAL CLUSTERING (Architecture B, k_spec=%d) ===\n", primary_k))
cat(sprintf("Total final labels: %d   CV accuracy: %.3f (chance = %.3f)\n",
            nlevels(factor(asg$final_label)), final_eval$accuracy,
            1 / nlevels(factor(asg$final_label))))
print(final_summary, n = Inf, width = Inf)

cat("\n=== Spectral-cluster-level rollup (training labels if you use spec only) ===\n")
spec_roll <- asg |>
  group_by(spec_cluster) |>
  summarise(n = dplyr::n(),
            n_sub = dplyr::n_distinct(stats::na.omit(sub_cluster)),
            .groups = "drop") |>
  left_join(spec_summary |> dplyr::select(spec_cluster, indicator_species,
                                          indicator_genus,
                                          dominance, heterogeneity, coherent),
            by = "spec_cluster") |>
  arrange(desc(n))
print(spec_roll, n = Inf, width = Inf)

saveRDS(list(
  assignments    = asg,
  final_summary  = final_summary,
  spec_summary   = spec_summary,
  final_eval     = final_eval,
  sub_details    = sub_details,
  primary_k      = primary_k,
  dom_threshold  = dom_threshold,
  het_threshold  = het_threshold
), "data/derived/final_clusters_B.rds")
