# 10_label_rasters.R — attach a Raster Attribute Table (class code, short
# label, full description, class_type, ecological category) AND a per-class
# color table to each domain's class raster from 09_inference.R.
#
# Does NOT re-run prediction. Opens each existing .tif, sets categories +
# colortable in memory, writes a sister "_labeled.tif" COG. The original
# numeric uint8 raster is left in place.
#
# Color scheme: 3x3 grid on (moisture, elevation) with muted earth tones.
#
#   hue        = moisture     dry -> mesic -> wet
#                              warm brown -> olive -> slate blue
#   lightness  = elevation    montane -> subalpine -> alpine
#                              dark -> medium -> light
#   chroma     = constant     ~moderate, avoids vivid colors
#
# Disturbed (Taraxacum-dominated, S20) gets a single distinct earth-tone
# (muted mustard) so it stands out without breaking the saturation budget.
#
# Within a cell, classes share the same color — within-cell distinction
# comes from short_label rather than hue.
#
# Inputs:
#   data/derived/aop_classified/{DOMAIN}_class_3m_v1.tif
#   data/derived/aop_classified/class_lookup.csv
#   data/derived/class_summary_table.csv
#   data/small_reference/class_categories.csv
# Outputs:
#   data/derived/aop_classified/{DOMAIN}_class_3m_v1_labeled.tif

suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(grDevices)
})
terra::terraOptions(progress = 0)

lookup <- readr::read_csv("data/derived/aop_classified/class_lookup.csv",
                          show_col_types = FALSE)
csum   <- readr::read_csv("data/derived/class_summary_table.csv",
                          show_col_types = FALSE) |>
  dplyr::select(final_label, short_label, description, class_type)
cats   <- readr::read_csv("data/small_reference/class_categories.csv",
                          show_col_types = FALSE)

# --- 1. Build the RAT (one row per class code) ---------------------------
rat <- lookup |>
  dplyr::left_join(csum, by = "final_label") |>
  dplyr::left_join(cats, by = "final_label") |>
  dplyr::transmute(
    value       = class_code,
    short_label = dplyr::coalesce(short_label, final_label),
    description = dplyr::coalesce(description, final_label),
    final_label = final_label,
    class_type  = class_type,
    moisture    = moisture,
    elevation   = elevation
  )
stopifnot(all(!is.na(rat$moisture)),
          all(!is.na(rat$elevation)))
cat(sprintf("RAT: %d classes\n  moisture x elevation cell counts:\n",
            nrow(rat)))
print(rat |> dplyr::count(moisture, elevation, name = "n") |>
        tidyr::pivot_wider(names_from = elevation, values_from = n,
                           values_fill = 0L) |>
        as.data.frame())

# --- 2. 3x3 (moisture x elevation) muted earth-tone grid -----------------
# Equal-chroma, lightness varies along the elevation axis, hue along
# the moisture axis. Disturbed gets one off-grid color.
hex_grid <- matrix(
  c(
    # montane     subalpine    alpine
    "#8a6f3e",   "#b5915f",   "#d4b89e",   # dry   (warm brown -> tan)
    "#6b7e4e",   "#94a872",   "#b8c9a2",   # mesic (olive -> sage -> light green)
    "#4a6a85",   "#7095ae",   "#a0bbd0"    # wet   (slate -> steel -> light blue)
  ),
  nrow = 3, byrow = TRUE,
  dimnames = list(moisture  = c("dry", "mesic", "wet"),
                  elevation = c("montane", "subalpine", "alpine"))
)
disturbed_hex <- "#c5b83d"   # muted mustard

assign_color <- function(moisture, elevation) {
  out <- rep(NA_character_, length(moisture))
  is_dist <- moisture == "disturbed"
  out[is_dist] <- disturbed_hex
  if (any(!is_dist)) {
    out[!is_dist] <- hex_grid[cbind(moisture[!is_dist],
                                    elevation[!is_dist])]
  }
  out
}
rat$color_hex <- assign_color(rat$moisture, rat$elevation)
stopifnot(all(!is.na(rat$color_hex)))

# --- 3. Apply to each domain raster -------------------------------------
domains <- c("ALMO", "CRBU", "UPTA")
in_dir  <- "data/derived/aop_classified"

for (dom in domains) {
  in_path  <- file.path(in_dir, sprintf("%s_class_3m_v1.tif",         dom))
  out_path <- file.path(in_dir, sprintf("%s_class_3m_v1_labeled.tif", dom))
  if (!file.exists(in_path)) {
    cat(sprintf("Skipping %s — %s not found\n", dom, in_path))
    next
  }
  cat(sprintf("\n%s: attaching RAT + colortable -> %s ...\n",
              dom, basename(out_path)))
  r <- terra::rast(in_path)

  # RAT (drop color_hex; it's not part of the categories table — it goes
  # via coltab below).
  levels(r) <- rat |> dplyr::select(-color_hex)
  # Make `short_label` the default category so QGIS picks that up first.
  terra::activeCat(r) <- "short_label"

  # Colortable: requires a data.frame of value + rgba columns.
  rgb_mat <- t(grDevices::col2rgb(rat$color_hex, alpha = FALSE))
  ct <- data.frame(value = rat$value,
                   red   = rgb_mat[, "red"],
                   green = rgb_mat[, "green"],
                   blue  = rgb_mat[, "blue"],
                   alpha = 255L)
  terra::coltab(r) <- ct

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

# --- 4. Persist the categories + colors as a sidecar CSV ----------------
readr::write_csv(rat, "data/derived/aop_classified/class_lookup_labeled.csv")
cat(sprintf("\nWrote data/derived/aop_classified/class_lookup_labeled.csv\n"))

# --- 5. Write a QGIS QML sidecar per domain so QGIS auto-applies the
#       paletted renderer + labels regardless of how it reads the
#       embedded color table.
write_qml <- function(rat, qml_path) {
  entries <- paste(
    sprintf(
      '        <paletteEntry color="%s" value="%d" alpha="255" label="%s"/>',
      rat$color_hex, rat$value,
      # Encode any XML-unsafe chars in labels.
      gsub("\"", "&quot;",
        gsub("'", "&apos;",
          gsub("<", "&lt;",
            gsub(">", "&gt;",
              gsub("&", "&amp;", rat$short_label)))))
    ),
    collapse = "\n"
  )
  qml <- sprintf(
'<!DOCTYPE qgis PUBLIC \'http://mrcc.com/qgis.dtd\' \'SYSTEM\'>
<qgis version="3.x" styleCategories="AllStyleCategories">
  <pipe>
    <rasterrenderer band="1" type="paletted" nodataColor="" opacity="1" alphaBand="-1">
      <rasterTransparency/>
      <colorPalette>
%s
      </colorPalette>
    </rasterrenderer>
    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>
    <huesaturation colorizeOn="0" saturation="0"/>
    <rasterresampler/>
    <resamplingStage>resamplingFilter</resamplingStage>
  </pipe>
  <blendMode>0</blendMode>
</qgis>
', entries)
  writeLines(qml, qml_path)
}

for (dom in domains) {
  tif <- file.path(in_dir, sprintf("%s_class_3m_v1_labeled.tif", dom))
  if (!file.exists(tif)) next
  qml <- sub("\\.tif$", ".qml", tif)
  write_qml(rat, qml)
  cat(sprintf("Wrote %s\n", basename(qml)))
}

cat("\nLoad *_labeled.tif in QGIS — the sibling .qml auto-applies the\n",
    "paletted renderer with `short_label` labels per class.\n", sep = "")
