# 07_resample_training_from_mosaic.R — re-source the joint classifier's
# spectral predictors (spec_PC01..20) by SAMPLING the deployed PC mosaics
# at each training crown, instead of R-projecting field-crown spectra.
#
# Why: training and inference must share one spectral feature space. The
# old design ran two projection paths (R for field spectra, Python for AOP
# pixels) that drifted apart on a basis refit -> meadow<->shrub swaps. By
# sampling training predictors from the SAME mosaics inference uses, the two
# are identical by construction, and CV recall becomes an honest estimate of
# deployed 3 m accuracy (crown-mean-of-1m-spectra vs an actual 3 m cell).
#
# Sampling unit: crown-overlap AREA-WEIGHTED mean of the 3 m cells the crown
# polygon covers (terra::extract(..., exact = TRUE)); for a crown smaller
# than one cell this collapses to the covering cell's value.
#
# Year/domain -> mosaic:
#   2025 crowns -> data/derived/aop_pc_maps_mosaic/{ALMO,UPTA}_pc_mosaic.tif,
#                  CRBU_pc_mosaic_2025.tif
#   2018 crowns -> CRBU_pc_mosaic_2018.tif  (NDVI-corrected; correction baked
#                  into the mosaic, so the classifier's 2018 inputs come
#                  through the identical path as 2018 deployment)
#
# Labels, snow_free_doy and canopy_height_m are reused from the existing
# joint_training_set.rds; only spec_PC01..20 are replaced. The R-side field
# PCA projection (05) and the shrub re-projection (02) drop out.
#
# Inputs:
#   data/derived/joint_training_set.rds        (labels, DOY, CHM, class_type)
#   data/derived/crowns_2025.gpkg, crowns_2018.gpkg
#   data/derived/aop_pc_maps_mosaic/*.tif
# Deployment note: the team ships the FIELD-SPECTRA classifier maps (they
# look more confident at 3 m); this script is the honest-CV report plus a
# ready-to-deploy mosaic-sampled set for if/when we switch the source of
# truth. It does NOT overwrite the canonical joint_training_set.rds.
#
# Outputs:
#   data/derived/cv_field_vs_mosaic_3m.csv      honest 3 m CV: per-class recall,
#                                               field-spectra vs mosaic-sampled
#   data/derived/joint_training_set_mosaic.rds  mosaic-sampled training set
#                                               (alternative; not deployed)

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(ranger)
})
terra::terraOptions(progress = 0)

N_PC <- 20
pc_cols <- sprintf("spec_PC%02d", seq_len(N_PC))
mos_dir <- "data/derived/aop_pc_maps_mosaic"

# --- 1. Existing training table (labels + DOY + CHM) ---------------------
js <- readRDS("data/derived/joint_training_set.rds")
tr <- js$training
cat(sprintf("Field-spectra training set: %d sites\n", nrow(tr)))

# --- 2. Crown geometry + domain per (site_number, Year) ------------------
c25 <- sf::st_read("data/derived/crowns_2025.gpkg", quiet = TRUE) |>
  dplyr::mutate(Year = 2025L) |>
  dplyr::select(site_number, Year, domain)
c18 <- sf::st_read("data/derived/crowns_2018.gpkg", quiet = TRUE) |>
  sf::st_transform(32613) |>
  dplyr::mutate(Year = 2018L, domain = "CRBU") |>
  dplyr::select(site_number, Year, domain)
crowns <- dplyr::bind_rows(c25, c18)

mosaic_path <- function(year, domain) {
  if (year == 2025L) {
    fn <- if (domain == "CRBU") "CRBU_pc_mosaic_2025.tif"
          else sprintf("%s_pc_mosaic.tif", domain)
  } else {
    fn <- "CRBU_pc_mosaic_2018.tif"     # 2018 is CRBU-only
  }
  file.path(mos_dir, fn)
}

# --- 3. Sample 20 PCs per crown (area-weighted mean over covered cells) ---
groups <- crowns |> sf::st_drop_geometry() |> dplyr::distinct(Year, domain)
sampled <- purrr::pmap_dfr(groups, function(Year, domain) {
  path <- mosaic_path(Year, domain)
  if (!file.exists(path)) {
    message(sprintf("  %d/%s: mosaic missing (%s) — skipping", Year, domain, path))
    return(NULL)
  }
  sub <- crowns |> dplyr::filter(Year == !!Year, domain == !!domain)
  if (nrow(sub) == 0) return(NULL)
  r  <- terra::rast(path)[[seq_len(N_PC)]]
  vv <- terra::vect(sub)
  # Crown-overlap area-weighted mean. For crowns smaller than a 3 m cell
  # the exact-coverage polygon summary returns NA (no cell center covered),
  # so fall back to the covering cell via the centroid point read.
  poly <- terra::extract(r, vv, fun = mean, na.rm = TRUE,
                         exact = TRUE, ID = FALSE)
  cent <- terra::extract(r, terra::centroids(vv), ID = FALSE)
  names(poly) <- pc_cols; names(cent) <- pc_cols
  v <- poly
  for (cc in pc_cols) v[[cc]] <- dplyr::coalesce(poly[[cc]], cent[[cc]])
  n_fb <- sum(is.na(poly[[1]]) & !is.na(cent[[1]]))
  cat(sprintf("  %d/%s: sampled %d crowns from %s (%d via centroid fallback)\n",
              Year, domain, nrow(sub), basename(path), n_fb))
  dplyr::bind_cols(
    sub |> sf::st_drop_geometry() |> dplyr::select(site_number, Year),
    v
  )
})

# A few sites have multiple crowns — average their per-crown PC vectors.
sampled_site <- sampled |>
  dplyr::group_by(site_number, Year) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(pc_cols),
                                 ~ mean(.x, na.rm = TRUE)),
                   .groups = "drop")

# --- 4. Swap spec_PC into the training table -----------------------------
tr_new <- tr |>
  dplyr::select(-dplyr::any_of(pc_cols)) |>
  dplyr::inner_join(sampled_site, by = c("site_number", "Year")) |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(pc_cols), ~ !is.na(.x)))

dropped <- nrow(tr) - nrow(tr_new)
cat(sprintf("Mosaic-sampled training set: %d sites (%d dropped to nodata/missing)\n",
            nrow(tr_new), dropped))

# --- 5. Honest CV: field-spectra vs mosaic-sampled, 22-feature (=09) ------
feat22 <- c(pc_cols, "snow_free_doy", "canopy_height_m")

cv_recall <- function(df) {
  df <- df |> dplyr::filter(!is.na(snow_free_doy), !is.na(canopy_height_m))
  X <- as.matrix(df[, feat22]); y <- factor(df$final_label)
  set.seed(42)
  fold <- integer(nrow(X))
  for (lvl in levels(y)) {
    idx <- which(y == lvl)
    fold[idx] <- ((sample(seq_along(idx)) - 1L) %% 5L) + 1L
  }
  preds <- factor(rep(NA_character_, nrow(X)), levels = levels(y))
  for (f in 1:5) {
    tr_i <- which(fold != f); te_i <- which(fold == f)
    if (!length(tr_i) || !length(te_i)) next
    fit <- ranger::ranger(x = X[tr_i, , drop = FALSE], y = y[tr_i],
                          num.trees = 800, classification = TRUE,
                          seed = 42 + f)
    preds[te_i] <- predict(fit, X[te_i, , drop = FALSE])$predictions
  }
  list(overall = mean(preds == y, na.rm = TRUE),
       per_class = tibble::tibble(final_label = y, pred = preds) |>
         dplyr::filter(!is.na(pred)) |>
         dplyr::group_by(final_label) |>
         dplyr::summarise(n = dplyr::n(),
                          recall = mean(pred == final_label), .groups = "drop"))
}

# Field-spectra baseline uses the ORIGINAL spec_PC from js$training.
cat("\nCV (22 features, unweighted, 5-fold) — field-spectra vs mosaic-sampled:\n")
cv_field <- cv_recall(tr)
cv_mos   <- cv_recall(tr_new)
cat(sprintf("  overall accuracy:  field-spectra = %.1f%%   mosaic-sampled = %.1f%%\n",
            100 * cv_field$overall, 100 * cv_mos$overall))

recall_cmp <- cv_field$per_class |>
  dplyr::rename(n_field = n, recall_field = recall) |>
  dplyr::full_join(cv_mos$per_class |>
                     dplyr::rename(n_mos = n, recall_mos = recall),
                   by = "final_label") |>
  dplyr::mutate(delta = recall_mos - recall_field) |>
  dplyr::arrange(delta)
cat("\nBiggest per-class recall drops (field -> mosaic):\n")
print(utils::head(as.data.frame(recall_cmp), 12), row.names = FALSE)

# --- 6. Persist ----------------------------------------------------------
# Non-destructive: the canonical joint_training_set.rds (field-spectra) is
# left untouched so it stays in sync with the deployed maps.
recall_cmp |>
  dplyr::mutate(overall_field  = round(cv_field$overall, 4),
                overall_mosaic = round(cv_mos$overall, 4)) |>
  readr::write_csv("data/derived/cv_field_vs_mosaic_3m.csv")
cat(sprintf("\nWrote data/derived/cv_field_vs_mosaic_3m.csv (field %.1f%% vs mosaic %.1f%%)\n",
            100 * cv_field$overall, 100 * cv_mos$overall))

js_new <- js
js_new$training         <- tr_new
js_new$spectra_source   <- "mosaic_sampled"
js_new$recall_22_mosaic <- cv_mos$per_class
js_new$overall_acc_22   <- list(field = cv_field$overall, mosaic = cv_mos$overall)
saveRDS(js_new, "data/derived/joint_training_set_mosaic.rds")
cat("Wrote data/derived/joint_training_set_mosaic.rds (alternative; not deployed)\n")
