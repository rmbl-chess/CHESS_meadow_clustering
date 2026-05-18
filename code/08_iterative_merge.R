# 08_iterative_merge.R — iteratively merge clusters until every cluster has
# CV-RF recall in the spectral feature space >= threshold (0.4). Final output
# is a set of training-sample-ready clusters that are both compositionally
# coherent and spectrally mappable from AOP imagery.
#
# Algorithm (starting from genus k=15):
#   1. Fit 5-fold CV RF on the 26 spectral features (20 PCs + 6 indices).
#   2. If min(recall) >= threshold, stop.
#   3. Identify the weakest cluster.
#   4. Find its merge partner by COMPOSITION SIMILARITY: pick the cluster
#      whose Hellinger centroid is closest to the weakest cluster's centroid.
#      Using spectral confusion here is wrong — it makes spectrally-broad
#      clusters into sinks that absorb compositionally-unrelated weak
#      clusters. Spectral signal still drives the recall test (which is what
#      decides _whether_ to merge); composition decides _into what_.
#   5. Merge: smaller absorbed into larger.
#   6. Repeat. Safety floor at min_k = 4 and max_iter = 20.
#
# Inputs:
#   data/derived/clusters_genus.rds      (initial k=15 labels)
#   data/derived/composition_genus.rds   (Hellinger for dominant-genus naming)
#   data/derived/spectral_features.rds   (26-feature design matrix)
# Outputs:
#   data/derived/final_clusters.rds      list with:
#     final_assignments  tibble(site_number, Year, final_cluster, original_k15)
#     final_summary      tibble(final_cluster, n_sites, n_2018, n_2025,
#                               recall, indicator_genus, top_features)
#     final_eval         RF CV results on the final labels
#     spec_means         per-cluster mean spectral feature vector
#     merge_history      list of merge events for reproducibility

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
})

cl_genus   <- readRDS("data/derived/clusters_genus.rds")
comp_genus <- readRDS("data/derived/composition_genus.rds")
spec_feat  <- readRDS("data/derived/spectral_features.rds")$features

spec_cols <- grep("^(spec_PC|ndvi|ndwi|pri|red_edge|cai|ndli)",
                  names(spec_feat), value = TRUE)

joined <- cl_genus$assignments |>
  dplyr::select(site_number, Year, cluster = k15) |>
  dplyr::inner_join(spec_feat, by = c("site_number", "Year"))
joined$cluster <- factor(joined$cluster)

X <- as.matrix(joined[, spec_cols])

# Hellinger feature matrix aligned to `joined` row order. Drop the non-species
# columns (NPV, Bare, Other_Forb, etc.) from the distance computation: they
# appear in nearly every cluster's centroid and drown out the distinctive
# named-genus signal. Non-species cover still participates in the clustering
# upstream (06) and in the spectral test; just not in centroid distance.
hell_full <- comp_genus$hellinger
nonsp_set <- c("Other_Forb_cover", "Other_Graminoid_cover", "NPV_cover",
               "Bare_cover", "Other_Moss_Lichen_cover",
               "Other_Deciduous_Shrub_cover")
hell_feat_cols <- setdiff(names(hell_full),
                          c("site_number", "Year", nonsp_set))
H <- as.matrix(
  joined |>
    dplyr::select(site_number, Year) |>
    dplyr::left_join(hell_full, by = c("site_number", "Year")) |>
    dplyr::select(dplyr::all_of(hell_feat_cols))
)

# Composition centroid (mean Hellinger vector) per current cluster label.
cluster_centroid_dist <- function(labels) {
  lv <- levels(labels)
  centroids <- vapply(lv,
                      function(l) colMeans(H[labels == l, , drop = FALSE]),
                      numeric(ncol(H)))
  centroids <- t(centroids)
  rownames(centroids) <- lv
  as.matrix(dist(centroids, method = "euclidean"))
}

eval_rf_cv <- function(labels, X, n_folds = 5, seed = 42, n_trees = 500) {
  y <- factor(labels)
  n <- length(y)
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
                          classification = TRUE, seed = seed,
                          verbose = FALSE)
    preds[te] <- predict(fit, data.frame(X[te, , drop = FALSE]))$predictions
  }
  ok <- !is.na(preds)
  truth <- y[ok]; pred <- preds[ok]
  cm <- table(truth = truth, pred = pred)
  recall <- diag(cm) / rowSums(cm)
  cm_norm <- sweep(cm, 1, rowSums(cm), "/")
  sym <- (cm_norm + t(cm_norm)) / 2; diag(sym) <- NA_real_
  list(accuracy = mean(truth == pred),
       recall = recall, confusion = cm,
       confusion_norm = cm_norm, pairwise_symmetric = sym)
}

threshold <- 0.4
min_k     <- 4
max_iter  <- 20
labels    <- joined$cluster
merge_history <- list()
iter <- 0

repeat {
  iter <- iter + 1
  res <- eval_rf_cv(labels, X)
  weakest_lvl    <- names(which.min(res$recall))
  weakest_recall <- res$recall[weakest_lvl]
  current_k      <- nlevels(labels)

  cat(sprintf("\n--- iter %d  k=%d  acc=%.3f  weakest=%s (recall=%.3f) ---\n",
              iter, current_k, res$accuracy, weakest_lvl, weakest_recall))

  if (weakest_recall >= threshold) {
    cat("All clusters meet threshold. Stopping.\n"); break
  }
  if (current_k <= min_k) {
    cat("Hit min_k floor. Stopping.\n"); break
  }
  if (iter > max_iter) {
    cat("Hit max_iter. Stopping.\n"); break
  }

  # Merge partner by composition similarity (Hellinger centroid distance).
  cdm <- cluster_centroid_dist(labels)
  comp_dists <- cdm[weakest_lvl, ]
  comp_dists[weakest_lvl] <- Inf
  partner_lvl  <- names(which.min(comp_dists))
  partner_dist <- comp_dists[partner_lvl]
  partner_conf <- res$pairwise_symmetric[weakest_lvl, partner_lvl]

  n_w <- sum(labels == weakest_lvl); n_p <- sum(labels == partner_lvl)
  if (n_w >= n_p) { from <- partner_lvl; to <- weakest_lvl
  } else          { from <- weakest_lvl; to <- partner_lvl }

  cat(sprintf("Merge %s (n=%d) -> %s (n=%d); comp_dist=%.3f sym_conf=%.3f\n",
              from, sum(labels == from), to, sum(labels == to),
              partner_dist, partner_conf))

  merge_history[[iter]] <- tibble::tibble(
    iter = iter, k_before = current_k, acc_before = res$accuracy,
    weakest = weakest_lvl, weakest_recall = weakest_recall,
    partner = partner_lvl, comp_dist = partner_dist, sym_conf = partner_conf,
    merged_from = from, merged_to = to,
    n_from = n_w, n_to = n_p
  )

  labels[labels == from] <- to
  labels <- droplevels(labels)
}

# Relabel: M01..MNN ordered by descending size.
sizes <- sort(table(labels), decreasing = TRUE)
relabel_map  <- setNames(sprintf("M%02d", seq_along(sizes)), names(sizes))
final_labels <- factor(relabel_map[as.character(labels)],
                       levels = sprintf("M%02d", seq_along(sizes)))

final_res <- eval_rf_cv(final_labels, X)

final_assignments <- tibble::tibble(
  site_number   = joined$site_number,
  Year          = joined$Year,
  final_cluster = final_labels,
  original_k15  = joined$cluster
)

# Dominant features per final cluster — by mean Hellinger value.
hell <- comp_genus$hellinger
hell_long <- hell |>
  pivot_longer(-c(site_number, Year), names_to = "feature", values_to = "h") |>
  inner_join(final_assignments |> dplyr::select(site_number, Year, final_cluster),
             by = c("site_number", "Year")) |>
  mutate(feature = stringr::str_replace(feature, "_cover$", ""))

nonsp_set <- c("Other_Forb", "Other_Graminoid", "NPV", "Bare",
               "Other_Moss_Lichen", "Other_Deciduous_Shrub")

top_features_per_cluster <- hell_long |>
  group_by(final_cluster, feature) |>
  summarise(mean_h = mean(h), .groups = "drop") |>
  group_by(final_cluster) |>
  arrange(desc(mean_h), .by_group = TRUE) |>
  summarise(
    indicator_genus = {
      sp <- feature[!feature %in% nonsp_set]
      if (length(sp) == 0) "(open / generalist)" else sp[1]
    },
    top_features = paste(head(feature, 5), collapse = ", "),
    .groups = "drop"
  )

# Year breakdown per cluster.
year_breakdown <- final_assignments |>
  count(final_cluster, Year) |>
  pivot_wider(names_from = Year, values_from = n, names_prefix = "n_",
              values_fill = 0L)

final_summary <- final_assignments |>
  count(final_cluster, name = "n_sites") |>
  left_join(year_breakdown, by = "final_cluster") |>
  mutate(recall = as.numeric(final_res$recall[as.character(final_cluster)])) |>
  left_join(top_features_per_cluster, by = "final_cluster") |>
  arrange(desc(n_sites))

# Per-cluster mean of spectral features (useful for downstream plotting).
spec_means <- final_assignments |>
  inner_join(spec_feat, by = c("site_number", "Year")) |>
  group_by(final_cluster) |>
  summarise(across(all_of(spec_cols), mean), .groups = "drop")

cat("\n\n=== FINAL CLUSTERS ===\n")
cat(sprintf("k=%d   CV accuracy=%.3f (chance=%.3f)   threshold=%.2f\n",
            nlevels(final_labels), final_res$accuracy,
            1 / nlevels(final_labels), threshold))
print(final_summary, n = Inf, width = Inf)

saveRDS(list(
  final_assignments = final_assignments,
  final_summary     = final_summary,
  final_eval        = final_res,
  spec_means        = spec_means,
  merge_history     = dplyr::bind_rows(merge_history),
  relabel_map       = relabel_map,
  threshold         = threshold,
  min_k             = min_k
), "data/derived/final_clusters.rds")
