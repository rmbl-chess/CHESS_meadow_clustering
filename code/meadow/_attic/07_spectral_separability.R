# 07_spectral_separability.R — for each composition-cluster solution, test
# whether the clusters are distinguishable in the spectral feature space.
#
# Method: stratified 5-fold cross-validated LDA on the 26 spectral features
# (20 PCs + 6 narrow-band indices). Report overall accuracy, per-cluster
# recall, and the symmetric pairwise confusion matrix (so we can identify
# merge candidates). LDA is the right first pass — brightness-normalized
# vegetation spectra are largely linear in the relevant feature space and
# LDA has no hyperparameters to tune. Drop in RF later if separability looks
# bad.
#
# Inputs:
#   data/derived/clusters_species.rds, .../clusters_genus.rds
#   data/derived/spectral_features.rds
# Outputs (per granularity x k):
#   data/derived/separability_<gran>_k<k>.rds  list with accuracy, recall,
#   confusion matrix, pairwise symmetric confusion, top merge candidates.

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
})
# MASS::select masks dplyr::select — use dplyr:: explicitly below.

cl_species <- readRDS("data/derived/clusters_species.rds")
cl_genus   <- readRDS("data/derived/clusters_genus.rds")
spec_feat  <- readRDS("data/derived/spectral_features.rds")$features

spec_cols <- grep("^(spec_PC|ndvi|ndwi|pri|red_edge|cai|ndli)",
                  names(spec_feat), value = TRUE)

eval_separability <- function(cl, k, n_folds = 5, seed = 42) {
  k_col <- sprintf("k%02d", k)
  joined <- cl$assignments |>
    dplyr::inner_join(spec_feat, by = c("site_number", "Year"))

  X <- as.matrix(joined[, spec_cols])
  y <- factor(joined[[k_col]])
  n <- nrow(X)

  set.seed(seed)
  # Stratified folds: shuffle within each class, then assign mod n_folds.
  fold <- integer(n)
  for (lvl in levels(y)) {
    idx <- which(y == lvl)
    fold[idx] <- ((sample(seq_along(idx)) - 1L) %% n_folds) + 1L
  }

  preds <- factor(rep(NA_character_, n), levels = levels(y))
  for (f in seq_len(n_folds)) {
    tr <- which(fold != f); te <- which(fold == f)
    # Skip classes with <2 train samples in this fold (LDA needs >=2).
    bad <- names(which(table(y[tr]) < 2))
    keep_tr <- tr[!(y[tr] %in% bad)]
    fit <- tryCatch(MASS::lda(X[keep_tr, ], y[keep_tr]), error = function(e) NULL)
    if (is.null(fit)) next
    preds[te] <- predict(fit, X[te, ])$class
  }

  ok    <- !is.na(preds)
  truth <- y[ok]; pred <- preds[ok]
  acc   <- mean(truth == pred)

  cm     <- table(truth = truth, pred = pred)
  recall <- diag(cm) / rowSums(cm)

  # Symmetric pairwise confusion: avg of P(pred j | true i) and P(pred i | true j).
  cm_norm <- sweep(cm, 1, rowSums(cm), "/")
  sym     <- (cm_norm + t(cm_norm)) / 2
  diag(sym) <- NA_real_

  # Tidy pair table, deduped (a < b by integer cluster id).
  pair_df <- as_tibble(as.table(sym), .name_repair = "minimal")
  names(pair_df) <- c("a", "b", "sym_conf")
  pair_df <- pair_df |>
    filter(!is.na(sym_conf)) |>
    mutate(a_i = as.integer(as.character(a)),
           b_i = as.integer(as.character(b))) |>
    filter(a_i < b_i) |>
    dplyr::select(cluster_a = a, cluster_b = b, sym_conf) |>
    arrange(desc(sym_conf))

  list(accuracy           = acc,
       n_samples          = sum(ok),
       n_classes          = nlevels(y),
       recall             = recall,
       cluster_sizes      = as.integer(table(y)),
       confusion          = cm,
       confusion_norm     = cm_norm,
       pairwise_symmetric = sym,
       merge_candidates   = pair_df)
}

print_summary <- function(res, label) {
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("Sites: %d   Classes: %d   CV accuracy: %.3f (chance = %.3f)\n",
              res$n_samples, res$n_classes, res$accuracy, 1 / res$n_classes))
  cat("\nPer-cluster recall (sorted):\n")
  print(round(sort(res$recall), 3))
  cat("\nTop merge candidates (top 8 pairs by symmetric confusion):\n")
  print(res$merge_candidates |> dplyr::slice_head(n = 8), n = 8)
}

for (gran in c("species", "genus")) {
  cl <- if (gran == "species") cl_species else cl_genus
  for (k in c(10, 15, 20)) {
    res <- eval_separability(cl, k)
    print_summary(res, sprintf("%s k=%d", gran, k))
    saveRDS(res, sprintf("data/derived/separability_%s_k%02d.rds", gran, k))
  }
}
