# 10_label_rasters.R — attach a Raster Attribute Table (class code, short
# label, full description, class_type, ecological category) AND a per-class
# color table to each domain's class raster from 09_inference.R.
#
# Does NOT re-run prediction. Opens each existing .tif, sets categories +
# colortable in memory, writes a sister "_labeled.tif" COG. The original
# numeric uint8 raster is left in place.
#
# Color scheme: each class gets a base hue from its ecological category,
# with within-category lightness variation so adjacent same-category
# classes are still visually distinguishable. Categories (10):
#
#   shrub_riparian      Salix complex, Alnus, Betula, Dasiphora, Cornus
#   shrub_mesic         Amelanchier, Sambucus, Symphoricarpos, Lonicera,
#                       Ribes, Prunus, Holodiscus
#   shrub_dry           Artemisia tridentata, Juniperus, Purshia, A. cana
#   wet_meadow          sedge / Caltha / Veratrum / Mertensia wetlands
#   tall_forb           Ligusticum / Veratrum / Corydalis tall forb stands
#   dwarf_subshrub      Vaccinium cespitosum dwarf-shrub meadows
#   sparse_rocky        Heterotheca / Deschampsia sparse rocky meadows
#   bunchgrass_meadow   Festuca / Oxytropis / Balsamorhiza grass-forb
#   dry_shrub_steppe    Chrysothamnus / Artemisia / Purshia low-elev steppe
#   disturbed           Taraxacum disturbed herbland
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
    category    = category
  )
stopifnot(all(!is.na(rat$category)))   # every class must be categorized
cat(sprintf("RAT: %d classes across %d categories\n",
            nrow(rat), dplyr::n_distinct(rat$category)))

# --- 2. Per-category base colors + within-category lightness variation ---
# Hand-picked hues with enough hue-spacing that ecologically distinct
# categories are unambiguous.
category_hex <- c(
  shrub_riparian    = "#2c7fb8",   # blue (water-edge)
  shrub_mesic       = "#6a3d9a",   # purple
  shrub_dry         = "#b15928",   # warm brown
  wet_meadow        = "#0570b0",   # deep blue
  tall_forb         = "#e31a1c",   # red-magenta (showy)
  dwarf_subshrub    = "#f4a582",   # salmon
  sparse_rocky      = "#969696",   # gray
  bunchgrass_meadow = "#7fbc41",   # yellow-green
  dry_shrub_steppe  = "#d4a017",   # gold-tan
  disturbed         = "#fee08b"    # pale yellow (anomaly)
)
stopifnot(all(unique(rat$category) %in% names(category_hex)))

# Within each category, fan the base color slightly through HCL lightness
# so adjacent same-category classes are still distinguishable in QGIS.
hue_fan <- function(base_hex, n) {
  if (n <= 1) return(base_hex)
  base_hcl <- as.numeric(grDevices::convertColor(
    t(grDevices::col2rgb(base_hex) / 255),
    from = "sRGB", to = "Lab"))      # placeholder; we use HSL via colorspace
  # Simpler: use grDevices::adjustcolor with varied alpha-blend-toward-white
  # for distinct but related shades.
  bases <- colorRampPalette(c(
    grDevices::adjustcolor(base_hex, red.f = 0.7, green.f = 0.7, blue.f = 0.7),
    base_hex,
    grDevices::adjustcolor(base_hex, red.f = 1.2, green.f = 1.2, blue.f = 1.2)
  ))(n)
  bases
}

rat <- rat |>
  dplyr::arrange(category, final_label) |>
  dplyr::group_by(category) |>
  dplyr::mutate(color_hex = hue_fan(category_hex[category[1]], dplyr::n())) |>
  dplyr::ungroup() |>
  dplyr::arrange(value)

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
