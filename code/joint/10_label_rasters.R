# 10_label_rasters.R — attach a Raster Attribute Table (class code,
# canonical label, narrative description, class_type) to each domain's
# class raster from 09_inference.R, so QGIS shows the human-readable
# label automatically when the raster is loaded.
#
# Does NOT re-run prediction. Just opens each existing .tif, sets the
# categorical levels in memory, and writes a sister "_labeled.tif" with
# the RAT embedded. The original (numeric uint8) raster is left in place.
#
# Inputs:
#   data/derived/aop_classified/{DOMAIN}_class_3m_v1.tif
#   data/derived/aop_classified/class_lookup.csv
#   data/derived/class_summary_table.csv
# Outputs:
#   data/derived/aop_classified/{DOMAIN}_class_3m_v1_labeled.tif

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
})
terra::terraOptions(progress = 0)

lookup <- readr::read_csv("data/derived/aop_classified/class_lookup.csv",
                          show_col_types = FALSE)
csum   <- readr::read_csv("data/derived/class_summary_table.csv",
                          show_col_types = FALSE) |>
  dplyr::select(final_label, description, class_type)

# RAT schema: must have a `value` column (numeric, matches raster values).
# Other columns become attributes; QGIS picks one as the "active"
# category for symbology.
rat <- lookup |>
  dplyr::left_join(csum, by = "final_label") |>
  dplyr::transmute(value       = class_code,
                   final_label = final_label,
                   description = dplyr::coalesce(description, final_label),
                   class_type  = class_type)

cat(sprintf("RAT: %d classes\n", nrow(rat)))
print(head(rat, 5))

domains <- c("ALMO", "CRBU", "UPTA")
in_dir  <- "data/derived/aop_classified"

for (dom in domains) {
  in_path  <- file.path(in_dir, sprintf("%s_class_3m_v1.tif",         dom))
  out_path <- file.path(in_dir, sprintf("%s_class_3m_v1_labeled.tif", dom))
  if (!file.exists(in_path)) {
    cat(sprintf("Skipping %s — %s not found\n", dom, in_path))
    next
  }
  cat(sprintf("\n%s: attaching RAT and writing %s ...\n",
              dom, basename(out_path)))
  r <- terra::rast(in_path)
  levels(r) <- rat
  # Make `description` the active category so QGIS labels with the
  # human-readable narrative by default; `final_label` is still in the
  # attribute table for filter/joins.
  terra::activeCat(r) <- "description"

  t0 <- Sys.time()
  terra::writeRaster(
    r, out_path, overwrite = TRUE,
    datatype = "INT1U", NAflag = 255,
    filetype = "COG",
    gdal = c("COMPRESS=DEFLATE", "LEVEL=6",
             "BLOCKSIZE=512", "OVERVIEW_RESAMPLING=NEAREST")
  )
  cat(sprintf("  done (%.1fs, %.1f MB)\n",
              as.numeric(Sys.time() - t0, units = "secs"),
              file.size(out_path) / 1e6))
}

cat("\nDone.  Load *_labeled.tif in QGIS — the raster will display\n",
    "class names from the `description` column by default.\n", sep = "")
