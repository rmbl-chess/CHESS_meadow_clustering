# 34_shrub_label_set.R — define the final shrub label set by aggregating
# spectrally-inseparable taxa identified in 33_shrub_separability.R.
#
# Aggregation strategy (justified by the centroid dendrogram and the
# RF CV confusion matrix):
#
#   1. Salix species collapse to a single genus-level class "Salix sp." —
#      all 12 binomials fall into one dendrogram cluster at k=8, and
#      nearly every Salix → Salix-other off-diagonal in the confusion
#      matrix bleeds within the genus.
#   2. Ribes species collapse to "Ribes sp." — same pattern, smaller N.
#   3. Juniperus communis + J. scopulorum collapse to "Juniperus sp."
#      (J. scopulorum N=4 and confuses primarily with J. communis).
#   4. Artemisia cana (N=2) is dropped; we keep only A. tridentata —
#      the woody Artemisia is what the woody_taxa.csv reference targets.
#   5. Classes with N < 3 after the genus-collapses are dropped from the
#      training set (kept in the records for documentation but not the
#      `final_label` set).
#
# All rules live in `agg_rules` and `min_n_per_class` below — edit those
# to change the outcome. Re-running re-fits the RF CV with the new
# labels and reports per-class recall.
#
# Inputs:
#   data/derived/shrub_records_canonical.rds
#   data/derived/shrub_veg_spectra.rds
#   data/derived/shrub_separability.rds
# Outputs:
#   data/derived/shrub_label_crosswalk.csv  (canonical_binomial -> final_label)
#   data/derived/shrub_training_set.rds     (final training table)
#   data/derived/shrub_training_set.csv
#   output/shrub_final_class_counts.png
#   output/shrub_final_confusion.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(ggplot2)
})

records <- readRDS("data/derived/shrub_records_canonical.rds")
vs      <- readRDS("data/derived/shrub_veg_spectra.rds")
joined  <- vs$joined
wl      <- vs$wavelengths

genus_collapse <- c("Salix", "Ribes", "Juniperus")
drop_binomials <- c("Artemisia cana")
min_n_per_class <- 3L

# --- Build the crosswalk -------------------------------------------------
crosswalk <- joined |>
  dplyr::distinct(canonical_binomial, canonical_genus) |>
  dplyr::mutate(
    final_label = dplyr::case_when(
      canonical_binomial %in% drop_binomials       ~ NA_character_,
      canonical_genus    %in% genus_collapse       ~ paste(canonical_genus, "sp."),
      TRUE                                          ~ canonical_binomial
    )
  )

# Apply minimum-N filter after the collapse.
class_n <- joined |>
  dplyr::left_join(crosswalk, by = c("canonical_binomial", "canonical_genus")) |>
  dplyr::count(final_label, name = "n_total") |>
  dplyr::filter(!is.na(final_label))
keep_labels <- class_n |> dplyr::filter(n_total >= min_n_per_class) |>
  dplyr::pull(final_label)

crosswalk <- crosswalk |>
  dplyr::mutate(final_label = dplyr::if_else(final_label %in% keep_labels,
                                              final_label, NA_character_))

cat("=== Final shrub label crosswalk ===\n")
print(crosswalk |>
        dplyr::left_join(joined |> dplyr::count(canonical_binomial,
                                                 name = "n_records"),
                         by = "canonical_binomial") |>
        dplyr::arrange(final_label, dplyr::desc(n_records)) |>
        as.data.frame())

# --- Apply crosswalk to the training data --------------------------------
train <- joined |>
  dplyr::left_join(crosswalk |> dplyr::select(canonical_binomial, final_label),
                   by = "canonical_binomial") |>
  dplyr::filter(!is.na(final_label))

cat(sprintf("\nFinal training sites: %d (was %d)\n",
            nrow(train), nrow(joined)))

class_counts <- train |> dplyr::count(final_label, name = "n") |>
  dplyr::arrange(dplyr::desc(n))
cat(sprintf("Final classes: %d. Range: N = %d .. %d. Median: %d\n",
            nrow(class_counts), min(class_counts$n), max(class_counts$n),
            median(class_counts$n)))
print(as.data.frame(class_counts))

# --- Re-fit RF CV on the new label set -----------------------------------
sep_old <- readRDS("data/derived/shrub_separability.rds")
keep_cols <- sep_old$pca$keep_cols
pca_rot   <- sep_old$pca$rotation
pca_ctr   <- sep_old$pca$center

# Project the training spectra into the same PC space (using the PCA fit
# in 33, so PCs are comparable across runs).
spec_mat <- as.matrix(train[, keep_cols])
spec_mat_ctr <- sweep(spec_mat, 2, pca_ctr, FUN = "-")
X <- spec_mat_ctr %*% pca_rot
colnames(X) <- sprintf("PC%02d", seq_len(ncol(X)))
y <- factor(train$final_label)
n <- nrow(X)
n_folds <- 5

# Drop classes with N < 2 from CV.
keep <- y %in% names(which(table(y) >= 2))
y_cv <- droplevels(y[keep])
X_cv <- X[keep, , drop = FALSE]
cat(sprintf("\nRF CV on final labels: %d sites; %d classes\n",
            nrow(X_cv), nlevels(y_cv)))

set.seed(42)
fold <- integer(nrow(X_cv))
for (lvl in levels(y_cv)) {
  idx <- which(y_cv == lvl)
  fold[idx] <- ((sample(seq_along(idx)) - 1L) %% n_folds) + 1L
}

run_cv <- function(use_weights) {
  preds <- factor(rep(NA_character_, nrow(X_cv)), levels = levels(y_cv))
  for (f in seq_len(n_folds)) {
    tr <- which(fold != f); te <- which(fold == f)
    good_levels <- names(which(table(y_cv[tr]) >= 1))
    te <- te[as.character(y_cv[te]) %in% good_levels]
    if (length(tr) == 0 || length(te) == 0) next
    # Inverse-frequency class weights: balance rare classes against
    # the Salix mega-class.
    cw <- NULL
    if (use_weights) {
      tab <- table(y_cv[tr])
      cw  <- as.numeric(sum(tab) / (length(tab) * tab))
      names(cw) <- names(tab)
    }
    fit <- ranger::ranger(
      x = X_cv[tr, , drop = FALSE], y = y_cv[tr],
      num.trees = 500, classification = TRUE,
      class.weights = cw, seed = 42 + f
    )
    preds[te] <- predict(fit, X_cv[te, , drop = FALSE])$predictions
  }
  preds
}

preds_unw <- run_cv(use_weights = FALSE)
preds_w   <- run_cv(use_weights = TRUE)

summarise_recall <- function(preds, name) {
  tibble::tibble(truth = y_cv, pred = preds) |>
    dplyr::filter(!is.na(pred)) |>
    dplyr::group_by(truth) |>
    dplyr::summarise(n_test = dplyr::n(),
                     recall = mean(pred == truth), .groups = "drop") |>
    dplyr::rename(final_label = truth) |>
    dplyr::mutate(model = name)
}
recall_df <- dplyr::bind_rows(
  summarise_recall(preds_unw, "unweighted"),
  summarise_recall(preds_w,   "balanced")
) |>
  dplyr::left_join(class_counts, by = "final_label") |>
  dplyr::mutate(final_label = as.character(final_label)) |>
  tidyr::pivot_wider(id_cols = c(final_label, n),
                     names_from = model,
                     values_from = c(recall, n_test),
                     names_glue = "{model}_{.value}") |>
  dplyr::arrange(dplyr::desc(n))

overall_acc_unw <- mean(preds_unw == y_cv, na.rm = TRUE)
overall_acc_w   <- mean(preds_w   == y_cv, na.rm = TRUE)

cat("\n=== RF CV recall on final labels (unweighted vs balanced) ===\n")
print(as.data.frame(recall_df))
cat(sprintf(
  "\nOverall accuracy:  unweighted = %.1f%%  balanced = %.1f%%  (raw binomials = %.1f%%)\n",
  100 * overall_acc_unw, 100 * overall_acc_w, 100 * sep_old$overall_acc
))

# The two confusion matrices: keep the balanced one for the plot since
# that's the operational classifier downstream.
preds       <- preds_w
overall_acc <- overall_acc_w

# --- Plots ---------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
p_counts <- ggplot(class_counts,
                   aes(x = forcats::fct_reorder(final_label, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  labs(x = NULL, y = "N records",
       title = "Shrub final-label class counts",
       subtitle = sprintf("Total = %d; %d classes (min N = %d)",
                          nrow(train), nrow(class_counts), min_n_per_class)) +
  theme_minimal(base_size = 11)
ggsave("output/shrub_final_class_counts.png", p_counts,
       width = 7, height = 5, dpi = 150)

conf <- tibble::tibble(truth = y_cv, pred = preds) |>
  dplyr::filter(!is.na(pred)) |>
  dplyr::count(truth, pred, name = "n") |>
  tidyr::complete(truth, pred, fill = list(n = 0L))
label_order <- class_counts$final_label
p_conf <- ggplot(conf,
                 aes(x = factor(pred, levels = label_order),
                     y = factor(truth, levels = rev(label_order)),
                     fill = n)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(n > 0, n, "")), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(x = "Predicted", y = "Truth",
       title = "RF CV confusion — final shrub label set",
       subtitle = sprintf("Overall accuracy = %.1f%% (5-fold stratified CV)",
                          100 * overall_acc)) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  coord_fixed()
ggsave("output/shrub_final_confusion.png", p_conf,
       width = 9, height = 9, dpi = 150)

# --- Persist -------------------------------------------------------------
readr::write_csv(crosswalk, "data/derived/shrub_label_crosswalk.csv")
saveRDS(list(training      = train,
             crosswalk     = crosswalk,
             class_counts  = class_counts,
             recall        = recall_df,
             confusion     = conf,
             overall_acc   = overall_acc),
        "data/derived/shrub_training_set.rds")
readr::write_csv(
  train |> dplyr::select(site_number, Year, sampling_area,
                         canonical_binomial, canonical_genus, final_label,
                         vegetation_height_cm, n_pixels),
  "data/derived/shrub_training_set.csv"
)
cat("\nWrote data/derived/shrub_label_crosswalk.csv\n")
cat("Wrote data/derived/shrub_training_set.{rds,csv}\n")
cat("Wrote output/shrub_final_class_counts.png\n")
cat("Wrote output/shrub_final_confusion.png\n")
