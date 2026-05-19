# 14_sweep_kspec.R — sweep k_spec for variant F and report per-k metrics that
# tell us where bimodality on snow-free DOY disappears:
#   - n_clusters
#   - Eta²(snow_free_doy)                   higher = more env-coherent
#   - max within-cluster Hartigan dip stat   lower = less bimodal
#   - mean within-cluster IQR of DOY (days)  lower = tighter env band
#   - 5-fold CV RF accuracy on PC + env      higher = more mappable
#
# Uses variant_F (PCs 2-12 z + DOY z), characterizations already computed in 09.
# Adds RF eval here so we can see accuracy at every k.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  if (!requireNamespace("diptest", quietly = TRUE)) {
    renv::install("diptest")
  }
  library(diptest)
})

sc        <- readRDS("data/derived/spectral_clusters.rds")
env       <- readRDS("data/derived/environment.rds")
spec_feat <- readRDS("data/derived/spectral_features.rds")$features

variant <- sc$variant_G
ks <- sc$ks

spec_cols <- grep("^(spec_PC|ndvi|ndwi|pri|red_edge|cai|ndli)",
                  names(spec_feat), value = TRUE)
# Match 10's classifier: drop year-shift features PC3 and PRI.
spec_cols <- setdiff(spec_cols, c("spec_PC03", "pri"))
joined <- variant$assignments |>
  dplyr::inner_join(spec_feat, by = c("site_number", "Year")) |>
  dplyr::inner_join(env,       by = c("site_number", "Year"))
X <- as.matrix(joined[, c(spec_cols, "snow_free_doy")])

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
  ok <- !is.na(preds)
  truth <- y[ok]; pred <- preds[ok]
  cm <- table(truth = truth, pred = pred)
  list(
    accuracy   = mean(truth == pred),
    recall     = diag(cm) / rowSums(cm)
  )
}

eta_sq <- function(x, g) {
  ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- factor(g[ok])
  grand <- mean(x)
  sum(table(g) * (tapply(x, g, mean) - grand)^2) / sum((x - grand)^2)
}

results <- purrr::map_dfr(ks, function(k) {
  k_col  <- sprintf("k%02d", k)
  labels <- joined[[k_col]]
  doy    <- joined$snow_free_doy

  per_cl <- tibble::tibble(cluster = labels, doy = doy) |>
    group_by(cluster) |>
    summarise(
      n    = dplyr::n(),
      sd   = sd(doy, na.rm = TRUE),
      iqr  = IQR(doy, na.rm = TRUE),
      dip  = if (sum(!is.na(doy)) >= 4)
               suppressWarnings(diptest::dip.test(doy))$statistic
             else NA_real_,
      .groups = "drop"
    )

  rf  <- eval_rf_cv(labels, X)
  rec <- as.numeric(rf$recall)

  tibble::tibble(
    k                  = k,
    eta_sq_doy         = eta_sq(doy, labels),
    mean_within_iqr    = mean(per_cl$iqr, na.rm = TRUE),
    n_bimodal          = sum(per_cl$dip > 0.08, na.rm = TRUE),
    cv_accuracy        = rf$accuracy,
    n_strong           = sum(rec >= 0.80,        na.rm = TRUE),
    n_marginal         = sum(rec >= 0.50 & rec < 0.80, na.rm = TRUE),
    n_weak             = sum(rec <  0.50,        na.rm = TRUE),
    n_tiny             = sum(per_cl$n  < 6),
    smallest_cluster   = min(per_cl$n)
  )
})

cat("\n=== Variant F k_spec sweep on snow-free DOY ===\n")
cat("(higher eta_sq, lower mean_within_*, lower max_dip, more clusters non-bimodal\n")
cat(" is the direction we want; bimodal cutoff is Hartigan dip > 0.08)\n\n")
print(results |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3))),
      n = Inf, width = Inf)

saveRDS(results, "data/derived/kspec_sweep.rds")
