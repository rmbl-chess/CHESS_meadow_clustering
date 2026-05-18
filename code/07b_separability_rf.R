# 07b_separability_rf.R — RF version of 07. Same CV folds and reporting as
# the LDA pass, so per-cluster recall is directly comparable. If RF gives
# substantially higher recall on the LDA-weak clusters, the spectral feature
# space has non-linear structure LDA can't capture (and we'd want RF for
# downstream mapping). If recall is similar, the weak clusters are genuinely
# spectrally overlapping and merging is the right move.
#
# Inputs / outputs mirror 07 with `_rf` suffix on the saved files.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
})

cl_species <- readRDS("data/derived/clusters_species.rds")
cl_genus   <- readRDS("data/derived/clusters_genus.rds")
spec_feat  <- readRDS("data/derived/spectral_features.rds")$features

spec_cols <- grep("^(spec_PC|ndvi|ndwi|pri|red_edge|cai|ndli)",
                  names(spec_feat), value = TRUE)

eval_separability_rf <- function(cl, k, n_folds = 5, seed = 42, n_trees = 500) {
  k_col <- sprintf("k%02d", k)
  joined <- cl$assignments |>
    dplyr::inner_join(spec_feat, by = c("site_number", "Year"))

  X <- as.matrix(joined[, spec_cols])
  y <- factor(joined[[k_col]])
  n <- nrow(X)

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
    fit <- ranger::ranger(.y ~ ., data = df_tr,
                          num.trees = n_trees,
                          probability = FALSE,
                          classification = TRUE,
                          seed = seed)
    df_te <- data.frame(X[te, , drop = FALSE])
    preds[te] <- predict(fit, df_te)$predictions
  }

  ok    <- !is.na(preds)
  truth <- y[ok]; pred <- preds[ok]
  acc   <- mean(truth == pred)
  cm    <- table(truth = truth, pred = pred)
  recall <- diag(cm) / rowSums(cm)
  cm_norm <- sweep(cm, 1, rowSums(cm), "/")
  sym <- (cm_norm + t(cm_norm)) / 2
  diag(sym) <- NA_real_
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
       confusion          = cm,
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
    res <- eval_separability_rf(cl, k)
    print_summary(res, sprintf("%s k=%d (RF)", gran, k))
    saveRDS(res, sprintf("data/derived/separability_%s_k%02d_rf.rds", gran, k))
  }
}
