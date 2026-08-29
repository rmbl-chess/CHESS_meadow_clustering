# 23_willow_change.R — willow-shrubland cover change across CRBU, 2018 -> 2025.
#
# Change-detection test of the gated temporal product (22_temporal_gating.R):
# quantify willow (Salix other + Salix wolfii) cover in the gated 2018 map vs
# the 2025 map, overall and spatially. Expectation (user): a slight overall
# willow decrease consistent with long-term drying, with real spatial
# variation in trend direction. The RAW (ungated) 2018 map is summarized too,
# as the calibration contrast.
#
# Inputs:
#   data/derived/aop_classified/CRBU_class_3m_v1.tif           (2025)
#   data/derived/aop_classified/CRBU_2018_class_3m_v1.tif      (raw 2018)
#   data/derived/aop_classified/CRBU_2018_class_3m_v1_gated.tif
#   data/derived/aop_classified/class_lookup.csv
#   R4D061 via rSDP (snow-free DOY, for melt-band stratification)
# Outputs:
#   data/derived/willow_change_summary.csv     (overall + by melt band)
#   data/derived/willow_change_300m.tif        (delta willow fraction, 300 m)
#   docs/figures/willow_change_crbu.png/.pdf   (map + distribution)

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(rSDP)
})
terra::terraOptions(progress = 0)

WILLOW <- c("Salix other", "Salix wolfii")
AGG_M  <- 300           # aggregation cell for the trend map (m)
out_dir <- "data/derived/aop_classified"

lookup <- readr::read_csv(file.path(out_dir, "class_lookup.csv"),
                          show_col_types = FALSE)
willow_codes <- lookup$class_code[lookup$final_label %in% WILLOW]
stopifnot(length(willow_codes) == 2)

c25  <- terra::rast(file.path(out_dir, "CRBU_class_3m_v1.tif"))
c18r <- terra::rast(file.path(out_dir, "CRBU_2018_class_3m_v1.tif"))
c18g <- terra::rast(file.path(out_dir, "CRBU_2018_class_3m_v1_gated.tif"))
c18r <- terra::resample(c18r, c25, method = "near")
c18g <- terra::resample(c18g, c25, method = "near")

# Common analysis mask: classified in BOTH years (gated product's domain).
valid <- !is.na(c25) & !is.na(c18g)
w25  <- terra::mask(c25  %in% willow_codes, valid, maskvalues = FALSE)
w18g <- terra::mask(c18g %in% willow_codes, valid, maskvalues = FALSE)
w18r <- terra::mask(c18r %in% willow_codes, valid, maskvalues = FALSE)

px_ha <- prod(terra::res(c25)) / 1e4
n_valid <- terra::global(valid, "sum", na.rm = TRUE)[[1]]
area_of <- function(r) terra::global(r, "sum", na.rm = TRUE)[[1]] * px_ha

a25 <- area_of(w25); a18g <- area_of(w18g); a18r <- area_of(w18r)
cat(sprintf("CRBU common classified area: %.0f ha (%.1fM px)\n",
            n_valid * px_ha, n_valid / 1e6))
cat(sprintf("Willow area 2018 (raw):   %8.1f ha\n", a18r))
cat(sprintf("Willow area 2018 (gated): %8.1f ha\n", a18g))
cat(sprintf("Willow area 2025:         %8.1f ha\n", a25))
cat(sprintf("Gated change 2018->2025:  %+.1f ha (%+.2f%%)   [raw: %+.1f ha, %+.2f%%]\n",
            a25 - a18g, 100 * (a25 - a18g) / a18g,
            a25 - a18r, 100 * (a25 - a18r) / a18r))

# --- Transitions (gated) ---------------------------------------------------
loss_r <- w18g & !w25          # willow -> something else
gain_r <- !w18g & w25          # something else -> willow
n_loss <- area_of(loss_r); n_gain <- area_of(gain_r)
cat(sprintf("Gross willow loss: %.1f ha; gross gain: %.1f ha (net %+.1f ha)\n",
            n_loss, n_gain, n_gain - n_loss))
top_of <- function(mask_r, class_r, label) {
  # ifel: NA cells of the logical mask must NOT leak through (mask() with
  # maskvalues=FALSE keeps NA cells and inflates the freq badly).
  keep <- terra::ifel(mask_r, 1, NA)
  f <- terra::freq(terra::mask(class_r, keep)) |>
    tibble::as_tibble() |>
    dplyr::left_join(lookup, by = c(value = "class_code")) |>
    dplyr::arrange(dplyr::desc(count)) |> head(4)
  cat(sprintf("  %s: %s\n", label,
      paste(sprintf("%s %.0f ha", f$final_label, f$count * px_ha), collapse = ", ")))
}
top_of(loss_r, c25,  "willow LOSS becomes (2025 class)")
top_of(gain_r, c18g, "willow GAIN was    (2018 class)")

# --- 300 m trend surface ---------------------------------------------------
fact <- round(AGG_M / terra::res(c25)[1])
f18 <- terra::aggregate(w18g, fact, fun = "mean", na.rm = TRUE)
f25 <- terra::aggregate(w25,  fact, fun = "mean", na.rm = TRUE)
cover <- terra::aggregate(valid, fact, fun = "mean", na.rm = TRUE)
delta <- terra::mask(f25 - f18, cover >= 0.25, maskvalues = FALSE)
names(delta) <- "delta_willow_fraction"
terra::writeRaster(delta, "data/derived/willow_change_300m.tif",
                   overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))

# --- Stratify by snow-free DOY (melt band = elevation/moisture gradient) ---
r4 <- rSDP::sdp_get_raster("R4D061"); if (terra::nlyr(r4) > 1) r4 <- r4[[1]]
doy300 <- terra::project(r4, delta, method = "bilinear")
strat <- tibble::tibble(
  delta = terra::values(delta)[, 1],
  f18   = terra::values(terra::mask(f18, cover >= 0.25, maskvalues = FALSE))[, 1],
  doy   = terra::values(doy300)[, 1]) |>
  dplyr::filter(!is.na(delta), !is.na(doy))
bands <- strat |>
  dplyr::mutate(melt_band = cut(doy, c(-Inf, 120, 145, Inf),
                                labels = c("early (<120)", "mid (120-145)",
                                           "late (>145)"))) |>
  dplyr::group_by(melt_band) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   mean_willow_2018 = mean(f18),
                   mean_delta = mean(delta),
                   frac_cells_losing = mean(delta < -0.01),
                   frac_cells_gaining = mean(delta > 0.01),
                   .groups = "drop")
print(bands, width = Inf)

summary_tbl <- dplyr::bind_rows(
  tibble::tibble(stratum = "CRBU total", n_cells = nrow(strat),
                 willow_2018_gated_ha = a18g, willow_2025_ha = a25,
                 delta_ha = a25 - a18g, delta_pct = 100 * (a25 - a18g) / a18g,
                 gross_loss_ha = n_loss, gross_gain_ha = n_gain,
                 willow_2018_raw_ha = a18r),
)
readr::write_csv(summary_tbl, "data/derived/willow_change_summary.csv")
readr::write_csv(bands, "data/derived/willow_change_by_meltband.csv")

# --- Figure: trend map + cell distribution ---------------------------------
# Diverging BrBG (brown = willow loss / drying, teal = gain, neutral ~0).
dl <- as.data.frame(delta, xy = TRUE) |> dplyr::rename(d = 3)
lim <- max(abs(quantile(dl$d, c(0.01, 0.99))))
p_map <- ggplot(dl, aes(x, y, fill = d)) +
  geom_raster() +
  scale_fill_gradient2(low = "#8C510A", mid = "#F5F5F2", high = "#01665E",
                       limits = c(-lim, lim), oob = scales::squish,
                       labels = scales::percent,
                       name = expression(Delta ~ "willow fraction")) +
  coord_equal(expand = FALSE) +
  labs(title = "Willow shrubland change, CRBU 2018 → 2025",
       subtitle = sprintf("300 m cells; gated 2018 map (τ = 0.14). Net %+.0f ha (%+.1f%%); brown = loss, teal = gain",
                          a25 - a18g, 100 * (a25 - a18g) / a18g),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))
p_hist <- ggplot(dl, aes(d)) +
  geom_histogram(binwidth = 0.01, fill = "grey55", color = NA) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey20") +
  geom_vline(xintercept = mean(dl$d), linewidth = 0.4,
             color = "#8C510A", linetype = "dashed") +
  scale_x_continuous(labels = scales::percent, limits = c(-lim, lim)) +
  labs(x = expression(Delta ~ "willow fraction per 300 m cell"), y = "cells",
       subtitle = sprintf("cell-level mean %+.2f%% (dashed)", 100 * mean(dl$d))) +
  theme_minimal(base_size = 9)
# Compose with base grid (no gridExtra/patchwork in this renv).
draw_both <- function() {
  grid::grid.newpage()
  lay <- grid::grid.layout(2, 1, heights = grid::unit(c(4, 1.3), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(p_map,  vp = grid::viewport(layout.pos.row = 1))
  print(p_hist, vp = grid::viewport(layout.pos.row = 2))
  grid::popViewport()
}
png("docs/figures/willow_change_crbu.png", width = 8, height = 10,
    units = "in", res = 300, bg = "white")
draw_both(); dev.off()
pdf("docs/figures/willow_change_crbu.pdf", width = 8, height = 10)
draw_both(); dev.off()
cat("Wrote docs/figures/willow_change_crbu.{png,pdf}\n")
