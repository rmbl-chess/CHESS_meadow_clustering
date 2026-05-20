# 09_inference.R — joint meadow + shrub classification map per domain.
#
# For each AOP domain (ALMO, CRBU, UPTA) this script:
#   1. Builds a 22-band 3 m raster stack:
#        bands  1-20  spec_PC01..PC20    (JPL PC mosaic, same basis as the
#                                          deployed aop_classifier_pca.csv)
#        band     21  snow_free_doy      (R4D061, resampled to 3 m)
#        band     22  canopy_height_m    (3 m max CHM from 01_canopy_height.R)
#   2. Refits the joint RF on the 858 training sites using ONLY those 22
#      features (drops the 6 narrow-band indices) so train and inference
#      share the feature space. Unweighted RF — matches 03_predict.
#   3. Predicts class + max class-probability per pixel via terra::predict
#      block-by-block.
#   4. Masks both outputs to NA wherever canopy_height_m > TREE_THRESHOLD_M
#      (4 m, per project convention — excludes tree pixels).
#   5. Writes one classified COG + one confidence COG per domain.
#
# Inputs:
#   data/derived/joint_training_set.rds
#   data/derived/aop_chm_3m/{DOMAIN}_chm_max_3m.tif       (from 01)
#   /Users/ian/Library/CloudStorage/.../JPL_delivered/PrincipalComponents/
#      {DOMAIN}_pc_mosaic_3m_v1.tif                       (20-band JPL mosaic)
#   R4D061 via rSDP                                       (snow-free DOY 27 m)
# Outputs:
#   data/derived/aop_classified/{DOMAIN}_class_3m_v1.tif
#   data/derived/aop_classified/{DOMAIN}_confidence_3m_v1.tif
#   data/derived/aop_classified/class_lookup.csv          (int code -> label)

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(ranger)
  library(rSDP)
})

terra::terraOptions(progress = 0)

TREE_THRESHOLD_M <- 4.0
N_PC <- 20

domains <- c("ALMO", "CRBU", "UPTA")
pc_dir  <- "/Users/ian/Library/CloudStorage/GoogleDrive-ibreckhe@gmail.com/My Drive/BreckheimerLab2025/Projects/CHESS/Data/AOP_mosaics/JPL_delivered/PrincipalComponents"
chm_dir <- "data/derived/aop_chm_3m"
out_dir <- "data/derived/aop_classified"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pc_paths  <- setNames(file.path(pc_dir,  sprintf("%s_pc_mosaic_3m_v1.tif", domains)), domains)
chm_paths <- setNames(file.path(chm_dir, sprintf("%s_chm_max_3m.tif",      domains)), domains)
stopifnot(all(file.exists(unlist(pc_paths))))
stopifnot(all(file.exists(unlist(chm_paths))))

# --- 1. Refit joint RF on 22 features -----------------------------------
js <- readRDS("data/derived/joint_training_set.rds")
features_22 <- c(sprintf("spec_PC%02d", seq_len(N_PC)),
                 "snow_free_doy", "canopy_height_m")
train <- js$training |>
  dplyr::filter(!is.na(snow_free_doy), !is.na(canopy_height_m))
X_tr <- as.matrix(train[, features_22])
y_tr <- factor(train$final_label)
cat(sprintf("RF training: %d sites x %d features, %d classes\n",
            nrow(X_tr), ncol(X_tr), nlevels(y_tr)))

# Probability=TRUE so we get class probabilities; argmax gives the class
# label, and rowMax gives the per-pixel confidence in one fit.
fit_joint <- ranger::ranger(
  x = X_tr, y = y_tr,
  num.trees   = 1000,
  classification = TRUE,
  probability = TRUE,
  seed = 42
)
cat(sprintf("OOB prediction error: %.1f%%\n",
            100 * fit_joint$prediction.error))

# Write the integer-code <-> label lookup once (uint8 class raster, levels
# 1..n start at 1 to match terra writeRaster conventions).
class_levels <- levels(y_tr)
class_lookup <- tibble::tibble(
  class_code  = seq_along(class_levels),
  final_label = class_levels
)
readr::write_csv(class_lookup, file.path(out_dir, "class_lookup.csv"))
cat(sprintf("Wrote %s (%d classes)\n",
            file.path(out_dir, "class_lookup.csv"), nrow(class_lookup)))

# --- 2. Pre-load and project R4D061 once --------------------------------
cat("\nLoading R4D061 (snow-free DOY) from rSDP ...\n")
r4d061 <- rSDP::sdp_get_raster("R4D061")
if (terra::nlyr(r4d061) > 1) r4d061 <- r4d061[[1]]

# --- 3. Per-domain inference --------------------------------------------
# Single combined predict: class label + max class probability in one
# RF traversal. Halves runtime vs two passes. num.threads = 0 means
# ranger uses all available cores per call.
n_cores <- max(1L, parallel::detectCores() - 1L)
predict_class_conf <- function(model, data, ...) {
  probs <- predict(model, data, num.threads = n_cores)$predictions
  if (is.null(dim(probs))) probs <- matrix(probs, nrow = 1)
  cbind(class = max.col(probs),
        conf  = apply(probs, 1, max))
}

run_domain <- function(dom) {
  cat(sprintf("\n=== %s ===\n", dom))

  pc_stack <- terra::rast(pc_paths[[dom]])
  if (terra::nlyr(pc_stack) < N_PC) {
    stop(sprintf("%s PC mosaic has %d bands, need >= %d",
                 dom, terra::nlyr(pc_stack), N_PC))
  }
  pc_stack <- pc_stack[[seq_len(N_PC)]]
  names(pc_stack) <- sprintf("spec_PC%02d", seq_len(N_PC))
  cat(sprintf("  PC mosaic: %d x %d px @ %.0fm, %d bands\n",
              terra::nrow(pc_stack), terra::ncol(pc_stack),
              terra::res(pc_stack)[1], terra::nlyr(pc_stack)))

  chm <- terra::rast(chm_paths[[dom]])
  chm <- terra::resample(chm, pc_stack, method = "near")
  names(chm) <- "canopy_height_m"

  cat("  Resampling R4D061 to 3 m grid ... ")
  t0 <- Sys.time()
  doy <- terra::project(r4d061, pc_stack, method = "near", threads = TRUE)
  names(doy) <- "snow_free_doy"
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))

  stack22 <- c(pc_stack, doy, chm)
  stopifnot(terra::nlyr(stack22) == length(features_22))
  names(stack22) <- features_22
  cat(sprintf("  22-band stack ready\n"))

  out_class <- file.path(out_dir, sprintf("%s_class_3m_v1.tif",       dom))
  out_conf  <- file.path(out_dir, sprintf("%s_confidence_3m_v1.tif", dom))
  tmp_2band <- tempfile(pattern = sprintf("%s_pred_", dom), fileext = ".tif")

  cat(sprintf("  Predicting class + confidence (one pass, %d threads) ... ",
              n_cores))
  t0 <- Sys.time()
  pred_2band <- terra::predict(
    stack22, fit_joint, fun = predict_class_conf,
    na.rm = TRUE, index = 1:2,
    filename = tmp_2band, overwrite = TRUE,
    wopt = list(datatype = "FLT4S",
                filetype = "GTiff",
                gdal = c("COMPRESS=DEFLATE", "LEVEL=6",
                         "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512",
                         "BIGTIFF=IF_SAFER"))
  )
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))

  # Tree mask AND write final COGs (class as uint8, confidence as float32).
  cat(sprintf("  Tree-masking (CHM_max > %.1f m -> NA) + writing COGs ... ",
              TREE_THRESHOLD_M))
  t0 <- Sys.time()
  pass <- chm <= TREE_THRESHOLD_M    # TRUE = keep, FALSE/NA = mask

  cls_r  <- terra::mask(pred_2band[[1]], pass, maskvalues = c(0, NA),
                        filename = out_class, overwrite = TRUE,
                        wopt = list(datatype = "INT1U", NAflag = 255,
                                    filetype = "COG",
                                    gdal = c("COMPRESS=DEFLATE", "LEVEL=6",
                                             "BLOCKSIZE=512",
                                             "OVERVIEW_RESAMPLING=NEAREST")))
  conf_r <- terra::mask(pred_2band[[2]], pass, maskvalues = c(0, NA),
                        filename = out_conf, overwrite = TRUE,
                        wopt = list(datatype = "FLT4S",
                                    filetype = "COG",
                                    gdal = c("COMPRESS=DEFLATE", "LEVEL=6",
                                             "PREDICTOR=YES", "BLOCKSIZE=512",
                                             "OVERVIEW_RESAMPLING=AVERAGE")))
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))
  unlink(tmp_2band)

  cat(sprintf("  Wrote %s (%.1f MB)\n", basename(out_class),
              file.size(out_class) / 1e6))
  cat(sprintf("  Wrote %s (%.1f MB)\n", basename(out_conf),
              file.size(out_conf)  / 1e6))
}

for (dom in domains) run_domain(dom)
cat("\nDone.\n")
