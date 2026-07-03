# 20_wavelength_ablation.R — does excluding the ~2000 nm importance peak (and
# widening the water-band gaps) cost the classifier anything?
#
# Fixed-basis input ablation: the ~2000 nm feature sits just above the
# 1800-1950 nm water gap and may be an atmospheric-correction artifact. To test
# its value WITHOUT refitting the PCA (which would change the basis and cascade
# to the inference mosaics), we neutralise the excluded bands' contribution to
# each PC score:
#     modified_PC_j = PC_j - sum_{b in excluded} (refl[b]-center[b]) * v_j[b]
# i.e. we set the excluded reflectance to the basis mean, keeping rotation +
# center identical, then retrain the RF and compare OOB accuracy / balanced
# recall. This measures how much the CURRENT classifier leans on those bands
# (not what a from-scratch narrower-gap PCA would do -- that's a separate,
# basis-changing experiment).
#
# Inputs:  joint_training_set.rds, spectral_features.rds, veg_spectra.rds,
#          shrub_veg_spectra.rds
# Output:  data/derived/wavelength_ablation.csv (scenario table)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
})

jt <- readRDS("data/derived/joint_training_set.rds")$training
sf <- readRDS("data/derived/spectral_features.rds")
keep_wl <- sf$keep_wl; center <- sf$pca$center; R20 <- sf$pca$rotation[, 1:20]

# --- Reconstruct 348-band spectra aligned to keep_wl, per training site ------
# Meadow: veg_spectra on its own grid (direct water mask). Shrub:
# shrub_veg_spectra, nearest-band-matched to keep_wl (mirrors 02_training
# lines 51-60; the match can repeat a shrub band, so use matrix indexing).
mk_S <- function(joined, wldf, nearest) {
  rc <- grep("^rfl_band_", names(joined), value = TRUE)
  M  <- as.matrix(joined[, rc])                       # sites x 426
  bn <- as.integer(stringr::str_extract(rc, "\\d+$"))
  wl <- wldf$center_wavelength_nm[match(bn, wldf$band_number)]
  if (nearest) {
    idx <- vapply(keep_wl, function(w) which.min(abs(wl - w)), integer(1))
  } else {
    wm  <- (wl >= 1340 & wl <= 1450) | (wl >= 1800 & wl <= 1950) | (wl > 2400)
    idx <- which(!wm); stopifnot(length(idx) == length(keep_wl))
  }
  Sm <- M[, idx, drop = FALSE]; colnames(Sm) <- paste0("k", seq_along(keep_wl))
  cbind(joined[, c("site_number", "Year")], as.data.frame(Sm))
}
vs <- readRDS("data/derived/veg_spectra.rds")
sv <- readRDS("data/derived/shrub_veg_spectra.rds")
specW <- dplyr::bind_rows(mk_S(vs$joined, vs$wavelengths, FALSE),
                          mk_S(sv$joined, sv$wavelengths, TRUE))

kcols <- paste0("k", seq_along(keep_wl))
m <- jt |>
  dplyr::select(site_number, Year, final_label, snow_free_doy, canopy_height_m) |>
  dplyr::inner_join(specW, by = c("site_number", "Year"))
S  <- as.matrix(m[, kcols])
Sc <- sweep(S, 2, center, FUN = "-")
fullPC <- Sc %*% R20

# Validate: reconstruction must reproduce the deployed PCs.
stored <- as.matrix(jt[match(paste(m$site_number, m$Year),
                             paste(jt$site_number, jt$Year)),
                       sprintf("spec_PC%02d", 1:20)])
recon_cor <- min(vapply(1:20, function(j) stats::cor(fullPC[, j], stored[, j]), numeric(1)))
cat(sprintf("Reconstruction vs deployed PCs: min corr %.4f, max abs diff %.5f\n",
            recon_cor, max(abs(fullPC - stored))))
stopifnot(recon_cor > 0.999)

y   <- factor(m$final_label)
doy <- m$snow_free_doy; chm <- m$canopy_height_m

# --- Evaluate a feature set (OOB accuracy + balanced recall, avg over seeds) -
eval_feats <- function(modPC, seeds = 1:5) {
  X <- cbind(modPC, snow_free_doy = doy, canopy_height_m = chm)
  ok_rows <- stats::complete.cases(X)
  a <- c(); br <- c(); rec_acc <- NULL
  for (s in seeds) {
    fit <- ranger::ranger(x = X[ok_rows, ], y = droplevels(y[ok_rows]),
                          num.trees = 1000, classification = TRUE, seed = s,
                          num.threads = 0)
    pr <- fit$predictions; yy <- droplevels(y[ok_rows])
    ok <- !is.na(pr)
    a  <- c(a, mean(pr[ok] == yy[ok]))
    cm <- table(truth = yy[ok], pred = pr[ok])
    rec <- diag(cm) / rowSums(cm)
    br <- c(br, mean(rec, na.rm = TRUE))
    rec_acc <- if (is.null(rec_acc)) rec else rec_acc + rec
  }
  list(acc = mean(a), bal_recall = mean(br), recall = rec_acc / length(seeds))
}

# --- Exclusion scenarios (extra wavelength zones to neutralise) --------------
# Base water gaps: 1340-1450, 1800-1950, >2400. keep_wl are the retained bands.
scenarios <- list(
  "baseline (current gaps)"       = numeric(0),
  "exclude 2000 nm peak (1950-2050)" = list(c(1950, 2050)),
  "widen gaps +100 nm (50/side)"  = list(c(1290, 1340), c(1450, 1500),
                                         c(1750, 1800), c(1950, 2000), c(2350, 2400)),
  "widen gaps +200 nm (100/side)" = list(c(1240, 1340), c(1450, 1550),
                                         c(1700, 1800), c(1950, 2050), c(2300, 2400))
)
in_zones <- function(zones) {
  if (length(zones) == 0) return(integer(0))
  which(Reduce(`|`, lapply(zones, function(z) keep_wl >= z[1] & keep_wl <= z[2])))
}

base <- NULL; rows <- list()
for (nm in names(scenarios)) {
  E <- in_zones(scenarios[[nm]])
  modPC <- if (length(E)) fullPC - Sc[, E, drop = FALSE] %*% R20[E, , drop = FALSE] else fullPC
  res <- eval_feats(modPC)
  if (is.null(base)) base <- res
  rows[[nm]] <- tibble::tibble(
    scenario = nm, bands_excluded = length(E),
    oob_accuracy = round(res$acc, 4),
    bal_recall = round(res$bal_recall, 4),
    d_accuracy = round(res$acc - base$acc, 4),
    d_bal_recall = round(res$bal_recall - base$bal_recall, 4))
  cat(sprintf("%-34s bands -%3d | acc %.3f (%+.4f) | bal.recall %.3f (%+.4f)\n",
              nm, length(E), res$acc, res$acc - base$acc,
              res$bal_recall, res$bal_recall - base$bal_recall))
  if (nm == names(scenarios)[2]) res2000 <- res   # for per-class delta
}
tab <- dplyr::bind_rows(rows)
readr::write_csv(tab, "data/derived/wavelength_ablation.csv")

# Which classes lose the most recall when the 2000 nm peak is excluded?
cat("\nClasses most hurt by excluding the 2000 nm peak (recall drop):\n")
drop <- sort(base$recall - res2000$recall, decreasing = TRUE)
print(round(utils::head(drop[drop > 0], 8), 3))
cat("\nWrote data/derived/wavelength_ablation.csv\n")
