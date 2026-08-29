# 22_temporal_gating.R — margin/hysteresis gating of the 2018 CRBU map
# against the 2025 map, so only confident spectral shifts flip a pixel's
# class between years (approach chosen in 21_change_stability_diagnosis.R:
# bare/truly-stable ground flipped 49% under independent per-year argmax —
# half the apparent change was decision-boundary noise; a Δ-magnitude gate
# was rejected there as too blunt).
#
# Rule: a 2018 pixel KEEPS the 2025 map's class unless the 2018 probability
# of its own argmax class beats the 2018 probability of the 2025 class by a
# margin tau:   gated = if (own18 != c25 & P18(own18) - P18(c25) > tau)
#                        own18 else c25
# tau is CALIBRATED at the year-effect points: smallest tau on a 0.01 grid
# that drops the BARE (truly-stable) flip rate to <= 5%; the script reports
# how much vegetated change survives.
#
# Must run AFTER 09_inference.R so both years' maps come from the SAME RF
# (levels here match class_lookup.csv because the refit is deterministic:
# same training table, features and seed as 09).
#
# Inputs:
#   data/derived/joint_training_set.rds
#   data/derived/year_effect_points.csv                (12_year_effect_points)
#   data/derived/aop_pc_maps_mosaic/CRBU_pc_mosaic_{2018,2025}.tif
#   data/derived/aop_chm_3m/CRBU_chm_max_3m.tif
#   data/derived/aop_classified/CRBU_class_3m_v1.tif   (2025, from 09)
#   data/derived/aop_classified/CRBU_2018_class_3m_v1.tif  (mask pattern)
#   data/derived/aop_classified/class_lookup.csv
#   R4D061 via rSDP
# Outputs:
#   data/derived/aop_classified/CRBU_2018_class_3m_v1_gated.tif
#   data/derived/aop_classified/CRBU_2018_margin_3m_v1.tif
#   data/derived/temporal_gating_calibration.csv       (flip rate vs tau)

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(ranger)
  library(rSDP)
})
terra::terraOptions(progress = 0)

N_PC <- 20
features_22 <- c(sprintf("spec_PC%02d", seq_len(N_PC)),
                 "snow_free_doy", "canopy_height_m")
BARE_FLIP_TARGET <- 0.05
out_dir <- "data/derived/aop_classified"

# --- 1. Probability RF (deterministic twin of 09's classifier) -------------
js <- readRDS("data/derived/joint_training_set.rds")
train <- js$training |>
  dplyr::filter(!is.na(snow_free_doy), !is.na(canopy_height_m))
X_tr <- train[, features_22]
y_tr <- factor(train$final_label)
cat(sprintf("Probability RF: %d sites x %d features, %d classes\n",
            nrow(X_tr), length(features_22), nlevels(y_tr)))
fit_prob <- ranger::ranger(x = X_tr, y = y_tr, num.trees = 500,
                           probability = TRUE, seed = 42, verbose = FALSE)
lvls <- levels(y_tr)

lookup <- readr::read_csv(file.path(out_dir, "class_lookup.csv"),
                          show_col_types = FALSE)
stopifnot(identical(lookup$final_label, lvls))  # 09 and 22 must agree

# --- 2. Calibrate tau at the year-effect points ----------------------------
pts <- readr::read_csv("data/derived/year_effect_points.csv",
                       show_col_types = FALSE) |>
  dplyr::filter(domain == "CRBU")
v <- terra::vect(pts, geom = c("x_utm", "y_utm"), crs = "EPSG:32613")

pc25 <- terra::rast("data/derived/aop_pc_maps_mosaic/CRBU_pc_mosaic_2025.tif")[[1:N_PC]]
pc18 <- terra::rast("data/derived/aop_pc_maps_mosaic/CRBU_pc_mosaic_2018.tif")[[1:N_PC]]
chm  <- terra::rast("data/derived/aop_chm_3m/CRBU_chm_max_3m.tif")

feat_at_points <- function(pc) {
  m <- terra::extract(pc, v, ID = FALSE)
  names(m) <- sprintf("spec_PC%02d", seq_len(N_PC))
  m$snow_free_doy   <- pts$snow_free_doy
  m$canopy_height_m <- terra::extract(chm, v, ID = FALSE)[[1]]
  m
}
f25 <- feat_at_points(pc25)
f18 <- feat_at_points(pc18)
ok <- stats::complete.cases(f25) & stats::complete.cases(f18)
cat(sprintf("Calibration points with both-year features: %d of %d\n",
            sum(ok), nrow(pts)))
p25 <- predict(fit_prob, f25[ok, ])$predictions
p18 <- predict(fit_prob, f18[ok, ])$predictions
c25   <- max.col(p25)
own18 <- max.col(p18)
margin <- p18[cbind(seq_len(nrow(p18)), own18)] -
          p18[cbind(seq_len(nrow(p18)), c25)]
is_bare <- pts$point_type[ok] == "non_vegetated"
flips   <- own18 != c25

cal <- purrr::map_dfr(seq(0, 0.40, by = 0.01), function(tau) {
  gated_flip <- flips & margin > tau
  tibble::tibble(tau = tau,
                 bare_flip_rate = mean(gated_flip[is_bare]),
                 veg_flip_rate  = mean(gated_flip[!is_bare]),
                 raw_bare_flip  = mean(flips[is_bare]),
                 raw_veg_flip   = mean(flips[!is_bare]))
})
readr::write_csv(cal, "data/derived/temporal_gating_calibration.csv")
tau_star <- cal$tau[which(cal$bare_flip_rate <= BARE_FLIP_TARGET)[1]]
if (is.na(tau_star)) stop("No tau on the grid meets the bare flip target.")
sel <- cal[cal$tau == tau_star, ]
cat(sprintf(paste0(
  "Raw flip rates: bare %.1f%%, veg %.1f%%\n",
  "tau* = %.2f  ->  bare flip %.1f%% (target <= %.0f%%), veg flip %.1f%% ",
  "(%.0f%% of raw veg change retained)\n"),
  100 * sel$raw_bare_flip, 100 * sel$raw_veg_flip,
  tau_star, 100 * sel$bare_flip_rate, 100 * BARE_FLIP_TARGET,
  100 * sel$veg_flip_rate, 100 * sel$veg_flip_rate / sel$raw_veg_flip))

# --- 3. Apply the gate to the CRBU 2018 map --------------------------------
cat("\nBuilding the 2018 22-band stack + 2025-class band ...\n")
names(pc18) <- sprintf("spec_PC%02d", seq_len(N_PC))
chm3 <- terra::resample(chm, pc18, method = "near")
names(chm3) <- "canopy_height_m"
r4d061 <- rSDP::sdp_get_raster("R4D061")
if (terra::nlyr(r4d061) > 1) r4d061 <- r4d061[[1]]
doy <- terra::project(r4d061, pc18, method = "near", threads = TRUE)
names(doy) <- "snow_free_doy"

c25_r <- terra::rast(file.path(out_dir, "CRBU_class_3m_v1.tif"))
c25_r <- terra::resample(c25_r, pc18, method = "near")
# 0 = "no 2025 prior" sentinel so na.rm doesn't drop pixels the 2025 mask
# removed but 2018 kept; the gate then falls back to the 2018 argmax.
c25_r <- terra::subst(c25_r, NA, 0)
names(c25_r) <- "class25"

stack23 <- c(pc18, doy, chm3, c25_r)
n_cores <- max(1L, parallel::detectCores() - 1L)
gate_fun <- function(model, data, ...) {
  probs <- predict(model, data[, features_22], num.threads = n_cores)$predictions
  if (is.null(dim(probs))) probs <- matrix(probs, nrow = 1)
  own  <- max.col(probs)
  c25v <- as.integer(round(data$class25))
  has_prior <- c25v >= 1 & c25v <= ncol(probs)
  p_own   <- probs[cbind(seq_len(nrow(probs)), own)]
  p_prior <- rep(NA_real_, nrow(probs))
  p_prior[has_prior] <- probs[cbind(which(has_prior), c25v[has_prior])]
  margin <- p_own - p_prior
  gated <- ifelse(!has_prior, own,
                  ifelse(own != c25v & margin > tau_star, own, c25v))
  cbind(gated = gated, margin = ifelse(has_prior, margin, NA_real_))
}

out_tmp <- tempfile(pattern = "gated_", fileext = ".tif")
cat(sprintf("Predicting gated 2018 classes (tau = %.2f) ...\n", tau_star))
t0 <- Sys.time()
gated2 <- terra::predict(stack23, fit_prob, fun = gate_fun,
                         na.rm = TRUE, index = 1:2,
                         filename = out_tmp, overwrite = TRUE,
                         wopt = list(datatype = "FLT4S", filetype = "GTiff",
                                     gdal = c("COMPRESS=DEFLATE", "TILED=YES",
                                              "BIGTIFF=IF_SAFER")))
cat(sprintf("done (%.1f min)\n", as.numeric(Sys.time() - t0, units = "mins")))

# Mask to the deployed 2018 map's footprint (09's CHM + landcover mask).
mask18 <- terra::rast(file.path(out_dir, "CRBU_2018_class_3m_v1.tif"))
gated_class <- terra::mask(terra::as.int(gated2[[1]]), mask18)
out_gated  <- file.path(out_dir, "CRBU_2018_class_3m_v1_gated.tif")
out_margin <- file.path(out_dir, "CRBU_2018_margin_3m_v1.tif")
terra::writeRaster(gated_class, out_gated, overwrite = TRUE,
                   datatype = "INT1U",
                   gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
terra::writeRaster(terra::mask(gated2[[2]], mask18), out_margin,
                   overwrite = TRUE, datatype = "FLT4S",
                   gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

# --- 4. Before/after change accounting -------------------------------------
c25m <- terra::mask(terra::resample(
  terra::rast(file.path(out_dir, "CRBU_class_3m_v1.tif")), mask18,
  method = "near"), mask18)
chg_raw   <- mean(terra::values(mask18  != c25m), na.rm = TRUE)
chg_gated <- mean(terra::values(gated_class != c25m), na.rm = TRUE)
cat(sprintf("\n2018-vs-2025 changed-pixel fraction: raw %.1f%% -> gated %.1f%%\n",
            100 * chg_raw, 100 * chg_gated))
cat(sprintf("Wrote %s + %s\n", out_gated, out_margin))
