# 39_landscape_distance.R — assess how far each inference pixel is from
# the closest labeled training centroid in feature space. Goal: identify
# the regions of the landscape where the classifier is being asked to
# extrapolate beyond what the training set actually covers, then surface
# (a) per-class "extrapolation rate" and (b) a spatial map of novelty.
#
# Method (see the design discussion in the README):
#   1. Standardize joint-training features so each PC / index / DOY / CHM
#      contributes comparable scale.
#   2. Pooled within-class covariance Sigma (stable across small classes).
#   3. Per inference pixel + per class centroid:
#        d^2 = (x - mu_k)^T Sigma^-1 (x - mu_k)
#      Mahalanobis distance is sqrt(d^2).
#   4. nearest_class = argmin_k d ; nearest_d = min_k d.
#   5. Calibrate OOD threshold on TRAINING pixels' own min_d (held out
#      against pooled centroids): pixels with min_d > q95 are flagged OOD.
#   6. Aggregations:
#      - class_novelty: per nearest_class, count + fraction OOD
#      - hex_novelty:   250 m hex grid, fraction OOD + modal class
#
# Inputs:
#   data/derived/joint_training_set.rds         (training, feature_cols)
#   data/derived/inference_predictions.csv      (inference pixels w/ CHM)
# Outputs:
#   data/derived/inference_pixel_distances.csv  (per pixel)
#   data/derived/novelty_by_class.csv
#   data/derived/novelty_by_hex.gpkg

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})

js       <- readRDS("data/derived/joint_training_set.rds")
training <- js$training
feature_cols <- js$feature_cols
infer    <- readr::read_csv("data/derived/inference_predictions.csv",
                            show_col_types = FALSE)

X_tr  <- as.matrix(training[, feature_cols])
X_inf <- as.matrix(infer[, feature_cols])
y_tr  <- factor(training$final_label)
cat(sprintf("Training: %d sites, %d features, %d classes\n",
            nrow(X_tr), ncol(X_tr), nlevels(y_tr)))
cat(sprintf("Inference: %d pixels\n", nrow(X_inf)))

# --- 1. Standardize features --------------------------------------------
mu_tr <- colMeans(X_tr)
sd_tr <- apply(X_tr, 2, stats::sd)
sd_tr[sd_tr == 0] <- 1
Xtr_z  <- sweep(sweep(X_tr,  2, mu_tr, "-"), 2, sd_tr, "/")
Xinf_z <- sweep(sweep(X_inf, 2, mu_tr, "-"), 2, sd_tr, "/")

# --- 2. Pooled within-class covariance ---------------------------------
# Centroid per class on z-scaled features
centroids <- tibble::as_tibble(Xtr_z) |>
  dplyr::mutate(label = as.character(y_tr)) |>
  dplyr::group_by(label) |>
  dplyr::summarise(dplyr::across(dplyr::everything(), mean),
                   .groups = "drop")
cent_mat <- as.matrix(centroids[, -1])
rownames(cent_mat) <- centroids$label

# Residuals (x - mu_class) for the pooled covariance
mu_per_site <- cent_mat[match(as.character(y_tr), centroids$label), ,
                        drop = FALSE]
residuals_tr <- Xtr_z - mu_per_site
# Use a small ridge (lambda * I) for numerical stability with high-dim PCs.
Sigma   <- (t(residuals_tr) %*% residuals_tr) / (nrow(Xtr_z) - nrow(cent_mat))
lambda  <- 1e-3 * mean(diag(Sigma))
Sigma_r <- Sigma + diag(lambda, ncol(Sigma))
Sigma_inv <- solve(Sigma_r)

# --- 3. Mahalanobis distance to every centroid --------------------------
mahala_to_centroids <- function(X_z, cent_mat, Sigma_inv) {
  d2 <- matrix(NA_real_, nrow = nrow(X_z), ncol = nrow(cent_mat))
  for (k in seq_len(nrow(cent_mat))) {
    diff <- sweep(X_z, 2, cent_mat[k, ], "-")
    d2[, k] <- rowSums((diff %*% Sigma_inv) * diff)
  }
  colnames(d2) <- rownames(cent_mat)
  sqrt(pmax(d2, 0))
}

cat("Computing distances ... ")
t0 <- Sys.time()
d_tr  <- mahala_to_centroids(Xtr_z,  cent_mat, Sigma_inv)
d_inf <- mahala_to_centroids(Xinf_z, cent_mat, Sigma_inv)
cat(sprintf("done (%.1fs)\n",
            as.numeric(Sys.time() - t0, units = "secs")))

# --- 4. Nearest + second nearest class per pixel -----------------------
top2 <- function(D) {
  # Returns nearest class, nearest_d, second class, second_d, margin
  nearest_idx <- apply(D, 1, which.min)
  nearest_d   <- D[cbind(seq_len(nrow(D)), nearest_idx)]
  # Mask out nearest to find second
  D2 <- D
  D2[cbind(seq_len(nrow(D)), nearest_idx)] <- Inf
  second_idx <- apply(D2, 1, which.min)
  second_d   <- D2[cbind(seq_len(nrow(D)), second_idx)]
  tibble::tibble(
    nearest_class = colnames(D)[nearest_idx],
    nearest_d     = nearest_d,
    second_class  = colnames(D)[second_idx],
    second_d      = second_d,
    margin        = second_d - nearest_d
  )
}
tr_dist  <- top2(d_tr)
inf_dist <- top2(d_inf)

# --- 5. OOD threshold from training-side own nearest-class distance ----
ood_threshold <- stats::quantile(tr_dist$nearest_d, 0.95, na.rm = TRUE)
cat(sprintf("OOD threshold (training 95th pct of nearest_d): %.2f\n",
            ood_threshold))
inf_out <- dplyr::bind_cols(
  infer |> dplyr::select(x_utm, y_utm, domain, snow_free_doy,
                         canopy_height_m, predicted_label),
  inf_dist
) |>
  dplyr::mutate(ood_flag = nearest_d > ood_threshold)

cat(sprintf("\nInference pixels flagged OOD: %d / %d (%.1f%%)\n",
            sum(inf_out$ood_flag, na.rm = TRUE), nrow(inf_out),
            100 * mean(inf_out$ood_flag, na.rm = TRUE)))

# --- 6a. Per-class novelty summary ---------------------------------------
class_novelty <- inf_out |>
  dplyr::count(nearest_class, name = "n_assigned") |>
  dplyr::left_join(
    inf_out |> dplyr::group_by(nearest_class) |>
      dplyr::summarise(n_ood       = sum(ood_flag, na.rm = TRUE),
                       median_d    = stats::median(nearest_d, na.rm = TRUE),
                       median_margin = stats::median(margin, na.rm = TRUE),
                       .groups     = "drop"),
    by = "nearest_class"
  ) |>
  dplyr::mutate(pct_ood = 100 * n_ood / n_assigned) |>
  dplyr::arrange(dplyr::desc(n_ood))
cat("\n=== Per-class novelty (most extrapolation = top) ===\n")
print(as.data.frame(class_novelty))

# --- 6b. Spatial: 250 m hex grid summary --------------------------------
pts <- sf::st_as_sf(inf_out, coords = c("x_utm", "y_utm"),
                    crs = 32613, remove = FALSE)
# Per-domain hex grids so domains stay separate (they don't overlap anyway)
make_hex <- function(pts_dom) {
  bbox <- sf::st_bbox(pts_dom)
  hex  <- sf::st_make_grid(sf::st_as_sfc(bbox),
                           cellsize = 250, square = FALSE) |>
    sf::st_sf() |>
    sf::st_set_crs(32613) |>
    dplyr::mutate(hex_id = dplyr::row_number())
  joined <- sf::st_join(pts_dom, hex, join = sf::st_within)
  summ <- joined |> sf::st_drop_geometry() |>
    dplyr::group_by(hex_id) |>
    dplyr::summarise(
      n_pixels       = dplyr::n(),
      n_ood          = sum(ood_flag, na.rm = TRUE),
      pct_ood        = 100 * mean(ood_flag, na.rm = TRUE),
      median_d       = stats::median(nearest_d, na.rm = TRUE),
      modal_class    = names(sort(table(nearest_class), decreasing = TRUE))[1],
      .groups        = "drop"
    )
  hex |> dplyr::inner_join(summ, by = "hex_id")
}
hex_by_domain <- pts |> dplyr::group_split(domain) |>
  purrr::map(make_hex) |>
  dplyr::bind_rows()
cat(sprintf("\nHex grid: %d cells with >=1 pixel\n", nrow(hex_by_domain)))

# --- 7. Persist ---------------------------------------------------------
readr::write_csv(inf_out, "data/derived/inference_pixel_distances.csv")
readr::write_csv(class_novelty, "data/derived/novelty_by_class.csv")
sf::st_write(hex_by_domain, "data/derived/novelty_by_hex.gpkg",
             delete_dsn = TRUE, quiet = TRUE)
cat("\nWrote data/derived/inference_pixel_distances.csv\n")
cat("Wrote data/derived/novelty_by_class.csv\n")
cat("Wrote data/derived/novelty_by_hex.gpkg\n")
