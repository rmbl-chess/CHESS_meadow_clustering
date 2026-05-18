# 10_subcluster_composition.R — Architecture B step 2: within each spectral
# cluster (variant C, k=8 = drop PC1 brightness), sub-cluster sites by
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

primary_k       <- 8
primary_variant <- "variant_D"   # PCs 2-12 + z-scaled snow-free DOY
k_col           <- sprintf("k%02d", primary_k)
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

cat("\nCoherence classification (Variant C, species-level):\n")
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
  mutate(final_label = if_else(is.na(sub_cluster),
                               spec_cluster,
                               paste0(spec_cluster, ".", sub_cluster)))

# --- RF eval on final labels (5-fold CV) ------------------------------------
spec_cols <- grep("^(spec_PC|ndvi|ndwi|pri|red_edge|cai|ndli)",
                  names(spec_feat), value = TRUE)
# At deployment, the classifier will have spectra (per pixel) + env rasters
# (per pixel). Use both here so CV accuracy reflects the deployment setup.
joined <- asg |>
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
hell_long <- comp_species$hellinger |>
  pivot_longer(-c(site_number, Year), names_to = "feature", values_to = "h") |>
  mutate(feature = stringr::str_replace(feature, "_cover$", "")) |>
  inner_join(asg |> dplyr::select(site_number, Year, final_label),
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
