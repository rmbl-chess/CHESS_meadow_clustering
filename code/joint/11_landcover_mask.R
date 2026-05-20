# 11_landcover_mask.R — drop predictions outside meadow + shrub landcover.
#
# The joint classifier predicts at every pixel the tree mask doesn't
# exclude, but the actual sampling domain (where the model was trained)
# was meadow + shrub. Outside those landcover categories we'd just be
# extrapolating onto rock, forest understory, pasture, water, etc.
# This script applies a hard landcover mask derived from the SDP R3D018
# 1 m landcover product (class table from the product XML):
#
#   1  evergreen trees and shrubs         drop
#   2  deciduous trees > 2 m              drop
#   3  meadow, grassland and subshrub     KEEP
#   4  persistent open water              drop
#   5  persistent snow and ice            drop
#   6  rock, bare soil, sparse veg        drop
#   7  building or structure              drop
#   8  paved or impervious                drop
#   9  irrigated pasture / cultivated     drop
#   10 deciduous shrubs <= 2 m            KEEP
#   11 evergreen forest understory        drop
#   12 deciduous forest understory        drop
#
# Per-domain workflow:
#   1. Crop R3D018 (1 m) to the domain's 3 m class raster extent.
#   2. Reclassify to binary: 1 if class in {3, 10}, 0 otherwise.
#   3. Aggregate by factor 3 with the modal statistic so each 3 m cell
#      is "1" iff a majority of its 9 underlying 1 m landcover pixels
#      were meadow or shrub.
#   4. Snap to the existing class raster's grid (nearest resample), then
#      mask `_class_3m_v1.tif` and `_confidence_3m_v1.tif` in place.
#
# Inputs:
#   data/raw/SDP/UG_landcover_1m_v4.tif
#   data/derived/aop_classified/{DOMAIN}_class_3m_v1.tif        (overwritten)
#   data/derived/aop_classified/{DOMAIN}_confidence_3m_v1.tif   (overwritten)
# Outputs:
#   data/derived/aop_chm_3m/{DOMAIN}_landcover_mask_3m.tif      (the mask)
#   masked versions of class + confidence rasters in place

suppressPackageStartupMessages({
  library(terra)
})
terra::terraOptions(progress = 0)

keep_classes <- c(3L, 10L)
domains <- c("ALMO", "CRBU", "UPTA")
in_dir  <- "data/derived/aop_classified"
mask_dir <- "data/derived/aop_chm_3m"
dir.create(mask_dir, showWarnings = FALSE, recursive = TRUE)

lc_1m <- terra::rast("data/raw/SDP/UG_landcover_1m_v4.tif")
cat(sprintf("R3D018 (1 m landcover): %d x %d px\n",
            terra::nrow(lc_1m), terra::ncol(lc_1m)))

for (dom in domains) {
  cat(sprintf("\n=== %s ===\n", dom))
  class_path <- file.path(in_dir, sprintf("%s_class_3m_v1.tif",       dom))
  conf_path  <- file.path(in_dir, sprintf("%s_confidence_3m_v1.tif", dom))
  if (!file.exists(class_path)) {
    cat(sprintf("  skipping — %s not found\n", basename(class_path)))
    next
  }
  cls_r  <- terra::rast(class_path)
  conf_r <- terra::rast(conf_path)

  # Crop R3D018 to a generous buffered extent (gives aggregate a clean
  # block boundary to work with), then reclassify to binary.
  buf <- terra::ext(cls_r)
  cat(sprintf("  cropping landcover to domain ... "))
  t0 <- Sys.time()
  lc_dom <- terra::crop(lc_1m, buf)
  cat(sprintf("done (%.1fs, %d x %d px)\n",
              as.numeric(Sys.time() - t0, units = "secs"),
              terra::nrow(lc_dom), terra::ncol(lc_dom)))

  # Binary reclass: keep_classes -> 1, everything else -> 0 (NA stays NA).
  cat("  reclassifying to keep / drop ... ")
  t0 <- Sys.time()
  rcl <- matrix(c(1, 1, 0,
                  2, 2, 0,
                  3, 3, 1,
                  4, 9, 0,
                  10, 10, 1,
                  11, 12, 0),
                ncol = 3, byrow = TRUE)
  mask_1m <- terra::classify(lc_dom, rcl,
                              include.lowest = TRUE, right = NA,
                              others = 0L,
                              datatype = "INT1U")
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))

  cat("  aggregating to 3 m via modal ... ")
  t0 <- Sys.time()
  mask_3m <- terra::aggregate(mask_1m, fact = 3, fun = "modal", na.rm = TRUE)
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))

  # Snap to the existing class raster's exact grid.
  mask_3m <- terra::resample(mask_3m, cls_r, method = "near")
  mask_path <- file.path(mask_dir, sprintf("%s_landcover_mask_3m.tif", dom))
  terra::writeRaster(
    mask_3m, mask_path, overwrite = TRUE,
    datatype = "INT1U", NAflag = 255,
    filetype = "COG",
    gdal = c("COMPRESS=DEFLATE", "LEVEL=6",
             "BLOCKSIZE=512", "OVERVIEW_RESAMPLING=NEAREST")
  )

  # Apply: keep where mask == 1, else NA. Overwrite the existing rasters.
  cat("  masking class + confidence ... ")
  t0 <- Sys.time()
  cls_masked <- terra::mask(cls_r, mask_3m, maskvalues = c(0, NA))
  terra::writeRaster(
    cls_masked, class_path, overwrite = TRUE,
    datatype = "INT1U", NAflag = 255,
    filetype = "COG",
    gdal = c("COMPRESS=DEFLATE", "LEVEL=6",
             "BLOCKSIZE=512", "OVERVIEW_RESAMPLING=NEAREST")
  )
  conf_masked <- terra::mask(conf_r, mask_3m, maskvalues = c(0, NA))
  terra::writeRaster(
    conf_masked, conf_path, overwrite = TRUE,
    datatype = "FLT4S",
    filetype = "COG",
    gdal = c("COMPRESS=DEFLATE", "LEVEL=6",
             "PREDICTOR=YES", "BLOCKSIZE=512",
             "OVERVIEW_RESAMPLING=AVERAGE")
  )
  cat(sprintf("done (%.1fs)\n",
              as.numeric(Sys.time() - t0, units = "secs")))

  # Quick stat: what fraction of valid pixels survived?
  before_valid <- terra::global(cls_r,        fun = "notNA")[1, 1]
  after_valid  <- terra::global(cls_masked,   fun = "notNA")[1, 1]
  cat(sprintf("  valid pixels: %.1fM -> %.1fM (%.0f%% retained)\n",
              before_valid / 1e6, after_valid / 1e6,
              100 * after_valid / before_valid))
}

cat("\nDone. Re-run code/joint/10_label_rasters.R to refresh the\n",
    "_labeled.tif + .qml outputs with the masked class data.\n", sep = "")
