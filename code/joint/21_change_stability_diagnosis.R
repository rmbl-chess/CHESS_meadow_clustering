# 21_change_stability_diagnosis.R — why do meadow classes shift so much between
# the 2018 and 2025 CRBU inference maps, and how stable is the classifier at
# KNOWN-STABLE locations?
#
# Motivation: the 2018 vs 2025 CRBU classified maps show large, ecologically
# implausible class shifts. Per-pixel argmax classification is unstable at
# decision boundaries, so small year-to-year spectral perturbations (residual
# radiometric drift after the NDVI correction, phenology, sensor noise) flip
# borderline pixels and manufacture "change". We want only large, real shifts
# to flip pixels.
#
# Diagnostic: sample the 20 PCs at the year-effect point locations from the
# ACTUAL 2018 and 2025 CRBU PC mosaics (2018 already carries the radiometric
# correction). At a fixed location DOY (R4D061 climatology) and CHM are
# identical across years, so any class flip is driven PURELY by the spectral
# (PC) difference. Non-vegetated (bare) points are treated as the truly-stable
# noise floor (per user: assume non-veg stable, veg not).
#
# KEY FINDINGS (06_29 / 55-class classifier):
#   - Bare (stable) flip rate = 49% between 2018 and 2025 -> ~half of the
#     observed map "change" is pure noise. Meadow/shrub flips at 57%, only
#     ~8 pts above the noise floor => most change is spurious.
#   - Flip rate rises with spectral change Delta (31% -> 85% across quintiles),
#     so gating on change magnitude is viable in principle, BUT
#   - Bare is out-of-domain for the vegetation-fit PCA, so its Delta has a huge
#     tail (median 5, 95th 38). A Delta-gate calibrated on the bare 95th
#     over-suppresses (freezes 99.6% of vegetated pixels). => a naive
#     Delta-gate is too blunt; the noise floor from bare is unrepresentative of
#     vegetated noise.
#
# DECISION / NEXT STEP (not yet built): approach #1 = margin / hysteresis in
# probability space. Flip a pixel only when the new class beats the prior
# year's class by a probability margin tau. Lives in the classifier's own
# confidence space, so it sidesteps the out-of-domain Delta problem. To build:
#   1. Re-fit the deployed RF with probability = TRUE (or emit top-2 class probs
#      from 09).
#   2. At the year-effect points, compute flip rate vs classifier margin; find
#      the margin tau that drives the BARE flip rate to ~5%.
#   3. Check how much vegetated change survives above tau (real change kept).
#   4. Apply at inference: keep prior-year class unless P_new - P_prior > tau.
#
# Inputs:  data/derived/year_effect_points.csv
#          data/derived/aop_pc_maps_mosaic/CRBU_pc_mosaic_{2018,2025}.tif
#          data/derived/aop_chm_3m/CRBU_chm_max_3m.tif
#          data/derived/joint_training_set.rds
# Output:  console diagnostic (no files written).

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(ranger)
})

pc_cols <- sprintf("spec_PC%02d", 1:20)
feat    <- c(pc_cols, "snow_free_doy", "canopy_height_m")

pts <- read.csv("data/derived/year_effect_points.csv") |>
  dplyr::filter(doy_band != "late")                 # drop late-melt (snow risk)

xy  <- cbind(pts$x_utm, pts$y_utm)
m18 <- terra::extract(terra::rast("data/derived/aop_pc_maps_mosaic/CRBU_pc_mosaic_2018.tif"), xy)
m25 <- terra::extract(terra::rast("data/derived/aop_pc_maps_mosaic/CRBU_pc_mosaic_2025.tif"), xy)
chm <- terra::extract(terra::rast("data/derived/aop_chm_3m/CRBU_chm_max_3m.tif"), xy)[[1]]
names(m18) <- names(m25) <- pc_cols
ok  <- stats::complete.cases(m18) & stats::complete.cases(m25) & !is.na(chm)
pts <- pts[ok, ]; m18 <- m18[ok, ]; m25 <- m25[ok, ]; chm <- chm[ok]
cat(sprintf("Year-effect points with full mosaic+CHM coverage: %d\n", nrow(pts)))

# Deployed RF (22 features, unweighted; mirrors 09_inference.R).
jt  <- readRDS("data/derived/joint_training_set.rds")$training
fit <- ranger::ranger(x = as.matrix(jt[, feat]), y = factor(jt$final_label),
                      num.trees = 1000, classification = TRUE, seed = 42,
                      num.threads = 0)
mk  <- function(pc) { X <- as.data.frame(pc)
  X$snow_free_doy <- pts$snow_free_doy; X$canopy_height_m <- chm
  as.matrix(X[, feat]) }
c18  <- as.character(predict(fit, mk(m18))$predictions)
c25  <- as.character(predict(fit, mk(m25))$predictions)
flip <- c18 != c25

# z-scaled spectral change Delta (training PC SDs).
sd_tr <- apply(as.matrix(jt[, pc_cols]), 2, sd)
delta <- sqrt(rowSums(sweep(as.matrix(m25) - as.matrix(m18), 2, sd_tr, "/")^2))

diag <- tibble::tibble(cover_class = pts$cover_class, flip = flip, delta = delta)
cat("\n=== Flip rate 2018->2025 at stable locations (current 55-class RF) ===\n")
print(diag |> dplyr::group_by(cover_class) |>
        dplyr::summarise(n = dplyr::n(), flip_rate = round(mean(flip), 3),
                         delta_median = round(median(delta), 2), .groups = "drop"))

bare <- pts$cover_class == "bare"; veg <- pts$cover_class == "meadow_shrub"
tau  <- as.numeric(quantile(delta[bare], 0.95))
cat(sprintf("\nNoise floor (bare): flip %.1f%%, Delta median %.2f, 95th %.2f\n",
            100 * mean(flip[bare]), median(delta[bare]), tau))
cat(sprintf("Meadow/shrub: flip %.1f%% (excess over bare ~%.1f pts = candidate real change)\n",
            100 * mean(flip[veg]), 100 * (mean(flip[veg]) - mean(flip[bare]))))

cat("\n=== Flip rate by Delta quintile (meadow/shrub) — flips rise with Delta ===\n")
print(tibble::tibble(delta = delta[veg], flip = flip[veg]) |>
        dplyr::mutate(bin = cut(delta, quantile(delta, seq(0, 1, .2)), include.lowest = TRUE)) |>
        dplyr::group_by(bin) |>
        dplyr::summarise(n = dplyr::n(), flip_rate = round(mean(flip), 3), .groups = "drop"))

cat("\n=== Naive Delta-gate at bare-95th is too blunt (over-suppresses veg) ===\n")
for (nm in c("bare", "veg")) {
  s <- if (nm == "bare") bare else veg
  raw <- mean(flip[s]); gated <- mean(flip[s] & delta[s] > tau)
  cat(sprintf("  %-4s: flip %.1f%% -> %.1f%% under Delta>tau  (%.0f%% of flips removed)\n",
              nm, 100 * raw, 100 * gated, 100 * (raw - gated) / raw))
}
cat("\n=> Pursue approach #1 (margin/hysteresis) — see header for the plan.\n")
