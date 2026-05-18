# 08_iterative_merge.R — preserve-strong, merge-weak iterative refinement of
# the genus k=15 clusters. Goal: training-sample-ready clusters that are both
# compositionally coherent and spectrally mappable from AOP imagery.
#
# Algorithm (starting from genus k=15):
#   1. Fit 5-fold CV RF on the 26 spectral features. Partition clusters into:
#        S = strong (recall >= strong_threshold), preserved as "fixed."
#        W = weak (recall <  strong_threshold), the merge candidates.
#   2. Compute per-cluster Hellinger centroids using NAMED-GENUS features
#      only (the 6 non-species categories are nearly universal and would
#      drown out distinctiveness in the distance metric).
#   3. While min(recall) < threshold and k > min_k:
#        a. Identify the weakest current W cluster.
#        b. Find merge partner by composition distance — restricted to other
#           W clusters. Only fall back to S if no W partners remain.
#        c. Merge: smaller absorbed into larger.
#        d. Re-fit RF. A merged W cluster may now meet threshold (effectively
#           "graduates" out of W); a strong cluster that absorbs a weak one
#           may slip below threshold (which the next iteration will handle).
#
# Inputs / outputs same as before; merge_history adds `partner_set` column
# (W or S) so the merge structure is auditable.

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

# Hellinger feature matrix, named-genus features only (see header).
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
  cm_norm <- sweep(cm, 1, rowSums(cm), "/")
  sym <- (cm_norm + t(cm_norm)) / 2; diag(sym) <- NA_real_
  list(accuracy = mean(truth == pred), recall = recall,
       confusion = cm, confusion_norm = cm_norm, pairwise_symmetric = sym)
}

cluster_centroid_dist <- function(labels) {
  lv <- levels(labels)
  centroids <- vapply(lv,
                      function(l) colMeans(H[labels == l, , drop = FALSE]),
                      numeric(ncol(H)))
  centroids <- t(centroids); rownames(centroids) <- lv
  as.matrix(dist(centroids, method = "euclidean"))
}

threshold        <- 0.4
strong_threshold <- 0.5
min_k            <- 4
max_iter         <- 30

# Initial RF: partition strong/weak.
init_res <- eval_rf_cv(joined$cluster, X)
strong_set <- names(init_res$recall)[init_res$recall >= strong_threshold]
weak_set   <- setdiff(levels(joined$cluster), strong_set)

cat(sprintf("Initial: strong (recall >= %.2f) = {%s} ; weak = {%s}\n",
            strong_threshold,
            paste(strong_set, collapse = ","),
            paste(weak_set, collapse = ",")))

labels <- joined$cluster
merge_history <- list()
iter <- 0

repeat {
  iter <- iter + 1
  res <- eval_rf_cv(labels, X)

  cat(sprintf("\n--- iter %d  k=%d  acc=%.3f ---\n",
              iter, nlevels(labels), res$accuracy))

  if (min(res$recall) >= threshold) {
    cat("All clusters meet threshold. Stopping.\n"); break
  }
  if (nlevels(labels) <= min_k) {
    cat("Hit min_k floor. Stopping.\n"); break
  }
  if (iter > max_iter) {
    cat("Hit max_iter. Stopping.\n"); break
  }

  # Dynamic classification: any cluster currently below threshold is a merge
  # candidate, regardless of original strong/weak status. This catches the
  # case where a strong cluster gets diluted by absorbing a weak one.
  weak_now    <- names(which(res$recall <  threshold))
  strong_now  <- names(which(res$recall >= threshold))
  weakest_lvl <- names(which.min(res$recall))
  weakest_recall <- res$recall[weakest_lvl]

  cdm <- cluster_centroid_dist(labels)
  cdm[weakest_lvl, weakest_lvl] <- Inf

  # Prefer to merge with another weak cluster (so two weaks may "graduate"
  # together). Only fall back to strong if no weak partner remains, so strong
  # clusters stay intact unless absolutely necessary.
  candidates_w <- setdiff(weak_now, weakest_lvl)
  candidates_s <- strong_now
  if (length(candidates_w) > 0) {
    pool <- candidates_w; pool_label <- "W"
  } else {
    pool <- candidates_s; pool_label <- "S"
  }
  if (length(pool) == 0) {
    cat(sprintf("Weakest %s has no merge partner. Stopping.\n", weakest_lvl)); break
  }

  # Keep names explicitly: matrix indexing with a length-1 column vector
  # returns an unnamed scalar, which silently breaks the rest of the merge.
  pool_dists  <- setNames(as.numeric(cdm[weakest_lvl, pool]), pool)
  partner_lvl <- names(which.min(pool_dists))
  partner_dist <- pool_dists[partner_lvl]
  partner_conf <- res$pairwise_symmetric[weakest_lvl, partner_lvl]

  n_w <- sum(labels == weakest_lvl); n_p <- sum(labels == partner_lvl)
  if (n_w >= n_p) { from <- partner_lvl; to <- weakest_lvl
  } else          { from <- weakest_lvl; to <- partner_lvl }

  cat(sprintf("Weakest=%s (recall=%.3f).  Merge %s (n=%d) -> %s (n=%d)  [%s; comp_dist=%.3f, sym_conf=%.3f]\n",
              weakest_lvl, weakest_recall, from, sum(labels == from),
              to, sum(labels == to), pool_label, partner_dist, partner_conf))

  merge_history[[iter]] <- tibble::tibble(
    iter = iter, k_before = nlevels(labels), acc_before = res$accuracy,
    weakest = weakest_lvl, weakest_recall = weakest_recall,
    partner = partner_lvl, partner_set = pool_label,
    comp_dist = partner_dist, sym_conf = partner_conf,
    merged_from = from, merged_to = to, n_from = n_w, n_to = n_p
  )

  # Apply merge. Whichever label survives keeps its strong/weak classification:
  # a weak cluster absorbed into a strong one inherits strong status (and the
  # strong cluster will be re-evaluated next iteration; if its recall drops
  # below the strong threshold, it just stays in `strong_set` but the global
  # threshold check will catch it if needed).
  labels[labels == from] <- to
  labels <- droplevels(labels)
}

# --- Relabel and final summary (same as before) -----------------------------
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

hell_long <- comp_genus$hellinger |>
  pivot_longer(-c(site_number, Year), names_to = "feature", values_to = "h") |>
  inner_join(final_assignments |> dplyr::select(site_number, Year, final_cluster),
             by = c("site_number", "Year")) |>
  mutate(feature = stringr::str_replace(feature, "_cover$", ""))

nonsp_short <- c("Other_Forb", "Other_Graminoid", "NPV", "Bare",
                 "Other_Moss_Lichen", "Other_Deciduous_Shrub")

top_features_per_cluster <- hell_long |>
  group_by(final_cluster, feature) |>
  summarise(mean_h = mean(h), .groups = "drop") |>
  group_by(final_cluster) |>
  arrange(desc(mean_h), .by_group = TRUE) |>
  summarise(
    indicator_genus = {
      sp <- feature[!feature %in% nonsp_short]
      if (length(sp) == 0) "(open / generalist)" else sp[1]
    },
    top_features = paste(head(feature, 5), collapse = ", "),
    .groups = "drop"
  )

year_breakdown <- final_assignments |>
  count(final_cluster, Year) |>
  pivot_wider(names_from = Year, values_from = n, names_prefix = "n_",
              values_fill = 0L)

final_summary <- final_assignments |>
  count(final_cluster, name = "n_sites") |>
  left_join(year_breakdown, by = "final_cluster") |>
  mutate(recall = as.numeric(final_res$recall[as.character(final_cluster)]),
         was_strong = sapply(as.character(original_k15_for_cluster <-
           {tapply(as.character(final_assignments$original_k15),
                   final_assignments$final_cluster, function(x) unique(x))}),
           function(orig_set) any(orig_set %in% strong_set))) |>
  left_join(top_features_per_cluster, by = "final_cluster") |>
  arrange(desc(n_sites))

spec_means <- final_assignments |>
  inner_join(spec_feat, by = c("site_number", "Year")) |>
  group_by(final_cluster) |>
  summarise(across(all_of(spec_cols), mean), .groups = "drop")

cat("\n\n=== FINAL CLUSTERS ===\n")
cat(sprintf("k=%d   CV accuracy=%.3f (chance=%.3f)   threshold=%.2f   strong_threshold=%.2f\n",
            nlevels(final_labels), final_res$accuracy,
            1 / nlevels(final_labels), threshold, strong_threshold))
print(final_summary, n = Inf, width = Inf)

saveRDS(list(
  final_assignments = final_assignments,
  final_summary     = final_summary,
  final_eval        = final_res,
  spec_means        = spec_means,
  merge_history     = dplyr::bind_rows(merge_history),
  relabel_map       = relabel_map,
  strong_set_orig   = strong_set,
  weak_set_orig     = weak_set,
  threshold         = threshold,
  strong_threshold  = strong_threshold,
  min_k             = min_k
), "data/derived/final_clusters.rds")
