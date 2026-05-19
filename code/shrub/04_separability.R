# 33_shrub_separability.R — class-imbalance + spectral-separability analysis
# for the canonical shrub training set. Mirrors the meadow separability
# work (05_preprocess + 07b_separability_rf):
#
#   1. Apply water-band mask, then PCA on the brightness-normalized spectra.
#   2. RF stratified-CV per-class recall (singleton classes have no test
#      fold; we still train on them).
#   3. Confusion matrix of out-of-fold predictions.
#   4. Centroid-based hierarchical clustering on per-class mean spectra
#      (Ward linkage in PC space) to visualize taxa similarity.
#
# Inputs:
#   data/derived/shrub_veg_spectra.rds
# Outputs:
#   data/derived/shrub_separability.rds   list(recall, confusion, dendro, pca)
#   output/shrub_class_counts.png
#   output/shrub_confusion.png
#   output/shrub_centroid_dendro.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(ggplot2)
})

vs <- readRDS("data/derived/shrub_veg_spectra.rds")
joined <- vs$joined
wl     <- vs$wavelengths

dir.create("output", showWarnings = FALSE)

# --- 1. Water-band mask + PCA --------------------------------------------
rfl_cols  <- grep("^rfl_band_", names(joined), value = TRUE)
band_nums <- as.integer(stringr::str_extract(rfl_cols, "\\d+$"))
band_wl   <- wl$center_wavelength_nm[match(band_nums, wl$band_number)]
water_mask <- (band_wl >= 1340 & band_wl <= 1450) |
              (band_wl >= 1800 & band_wl <= 1950) |
              (band_wl >  2400)
keep_cols <- rfl_cols[!water_mask]
keep_wl   <- band_wl[!water_mask]
cat(sprintf("Spectral: %d bands retained after water mask (of %d total)\n",
            length(keep_cols), length(rfl_cols)))

spec_mat <- as.matrix(joined[, keep_cols])
# Drop any rows with NA (rare; would block PCA).
good <- complete.cases(spec_mat)
if (sum(!good) > 0) {
  cat(sprintf("Dropping %d rows with NA reflectance\n", sum(!good)))
  joined   <- joined[good, ]
  spec_mat <- spec_mat[good, ]
}

pca <- prcomp(spec_mat, center = TRUE, scale. = FALSE)
n_pc <- 20
spec_pcs <- pca$x[, seq_len(n_pc), drop = FALSE]
colnames(spec_pcs) <- sprintf("PC%02d", seq_len(n_pc))
cat(sprintf("PCA: PC1-%d explain %.1f%% of variance\n",
            n_pc, 100 * sum(pca$sdev[seq_len(n_pc)]^2) / sum(pca$sdev^2)))

# --- 2. Class summary ----------------------------------------------------
classes <- joined |>
  dplyr::mutate(label = canonical_binomial) |>
  dplyr::count(label, name = "n") |>
  dplyr::arrange(dplyr::desc(n))
cat(sprintf("\nClasses: %d. Range: N = %d .. %d. Median: %d\n",
            nrow(classes), min(classes$n), max(classes$n),
            median(classes$n)))
print(as.data.frame(classes))

# Class-count plot
p_counts <- ggplot(classes,
                   aes(x = forcats::fct_reorder(label, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  labs(x = NULL, y = "N records",
       title = "Shrub training-class counts (canonical binomials)",
       subtitle = sprintf("Total = %d; %d classes",
                          nrow(joined), nrow(classes))) +
  theme_minimal(base_size = 10)
ggsave("output/shrub_class_counts.png", p_counts,
       width = 7, height = 8, dpi = 150)

# --- 3. RF stratified CV on PCs ------------------------------------------
y <- factor(joined$canonical_binomial)
X <- spec_pcs
n <- nrow(X)

# Drop classes with N < 2 from CV (need at least one in train AND test).
keep <- y %in% names(which(table(y) >= 2))
y_cv <- droplevels(y[keep])
X_cv <- X[keep, , drop = FALSE]
n_cv <- nrow(X_cv)
n_folds <- 5
cat(sprintf("\nRF CV: %d sites usable (>=2 per class); %d classes\n",
            n_cv, nlevels(y_cv)))

set.seed(42)
fold <- integer(n_cv)
for (lvl in levels(y_cv)) {
  idx <- which(y_cv == lvl)
  fold[idx] <- ((sample(seq_along(idx)) - 1L) %% n_folds) + 1L
}

preds <- factor(rep(NA_character_, n_cv), levels = levels(y_cv))
for (f in seq_len(n_folds)) {
  tr <- which(fold != f); te <- which(fold == f)
  # Skip if any class missing from training fold this round.
  good_levels <- names(which(table(y_cv[tr]) >= 1))
  drop_te <- te[!(as.character(y_cv[te]) %in% good_levels)]
  te <- setdiff(te, drop_te)
  if (length(tr) == 0 || length(te) == 0) next
  fit <- ranger::ranger(
    x = X_cv[tr, , drop = FALSE],
    y = y_cv[tr],
    num.trees = 500,
    classification = TRUE,
    seed = 42 + f
  )
  preds[te] <- predict(fit, X_cv[te, , drop = FALSE])$predictions
}

# Per-class recall
recall_df <- tibble::tibble(
  truth = y_cv,
  pred  = preds
) |>
  dplyr::filter(!is.na(pred)) |>
  dplyr::group_by(truth) |>
  dplyr::summarise(n_test = dplyr::n(),
                   recall = mean(pred == truth),
                   .groups = "drop") |>
  dplyr::rename(label = truth) |>
  dplyr::left_join(classes, by = "label") |>
  dplyr::mutate(label = as.character(label)) |>
  dplyr::arrange(dplyr::desc(recall), dplyr::desc(n))

cat("\n=== RF CV per-class recall ===\n")
print(as.data.frame(recall_df))

overall_acc <- mean(preds == y_cv, na.rm = TRUE)
cat(sprintf("\nOverall RF CV accuracy: %.2f%%\n", 100 * overall_acc))

# Confusion matrix (raw counts) — full matrix, only on the predicted rows
conf <- tibble::tibble(truth = y_cv, pred = preds) |>
  dplyr::filter(!is.na(pred)) |>
  dplyr::count(truth, pred, name = "n") |>
  tidyr::complete(truth, pred, fill = list(n = 0L))

# Order labels by N (descending) for the heatmap.
label_order <- classes$label
p_conf <- ggplot(conf,
                 aes(x = factor(pred, levels = label_order),
                     y = factor(truth, levels = rev(label_order)),
                     fill = n)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(n > 0, n, "")), size = 2.5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "Predicted", y = "Truth",
       title = "RF CV confusion matrix (canonical binomials)",
       subtitle = sprintf("Overall accuracy = %.1f%% (5-fold stratified CV)",
                          100 * overall_acc)) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  coord_fixed()
ggsave("output/shrub_confusion.png", p_conf,
       width = 11, height = 11, dpi = 150)

# --- 4. Centroid-based hclust on per-class mean PCs ---------------------
centroids <- as.data.frame(spec_pcs) |>
  dplyr::mutate(label = as.character(y)) |>
  dplyr::group_by(label) |>
  dplyr::summarise(dplyr::across(dplyr::everything(), mean),
                   .groups = "drop")
cent_mat <- as.matrix(centroids[, -1])
rownames(cent_mat) <- centroids$label
# Weight features by RF importance-ish proxy: use SD across classes per PC.
# For visualization keep it simple — Euclidean on standardized PC means.
cent_z <- scale(cent_mat)
d_cent <- dist(cent_z)
hc <- hclust(d_cent, method = "ward.D2")

png("output/shrub_centroid_dendro.png",
    width = 11 * 150, height = 7 * 150, res = 150)
par(mar = c(11, 4, 4, 1))
plot(hc, hang = -1, cex = 0.85,
     main = "Shrub canonical binomials — Ward dendrogram on PC1-20 centroids",
     xlab = "", sub = "")
dev.off()

# --- 5. Persist ---------------------------------------------------------
saveRDS(list(
  recall      = recall_df,
  confusion   = conf,
  overall_acc = overall_acc,
  centroids   = centroids,
  dendro      = hc,
  pca         = list(rotation = pca$rotation[, seq_len(n_pc)],
                     center   = pca$center,
                     keep_wl  = keep_wl,
                     keep_cols = keep_cols,
                     var_explained = cumsum(pca$sdev^2) / sum(pca$sdev^2))
), "data/derived/shrub_separability.rds")
cat("\nWrote data/derived/shrub_separability.rds\n")
cat("Wrote output/shrub_class_counts.png\n")
cat("Wrote output/shrub_confusion.png\n")
cat("Wrote output/shrub_centroid_dendro.png\n")
