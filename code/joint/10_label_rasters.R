# 10_label_rasters.R — attach a Raster Attribute Table (class code, short
# label, full description, class_type, ecological category) AND a per-class
# color table to each domain's class raster from 09_inference.R.
#
# Does NOT re-run prediction. Opens each existing .tif, sets categories +
# colortable in memory, writes a sister "_labeled.tif" COG. The original
# numeric uint8 raster is left in place.
#
# Color scheme: 4 physiognomy/habitat families, each rendered as an HCL ramp
# so the HUE shows the broad structural category and the within-family
# lightness/hue ramp (ordered by elevation then snow-free DOY) distinguishes
# the communities inside it:
#
#   shrubland     blue - grey - purple
#   grassland     tan - yellow - brown      (graminoid-dominated)
#   forb_meadow   green - olive - brown     (forb-dominated)
#   wetland       pink - red
#
# Category is derived from class_type + moisture + dominant-indicator
# lifeform (shrub -> shrubland; wet -> wetland; graminoid top indicator ->
# grassland; else forb_meadow). An optional `category` column in
# class_categories.csv overrides the derived value per class.
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
# Provisional (top-1) NatureServe community per class from the crosswalk
# draft (20_ecosystem_crosswalk.R) — surfaced on the map so the draft
# classes can be reviewed spatially before final candidate selection.
ns_path  <- "data/derived/natureserve_candidates.csv"
ns_draft <- if (file.exists(ns_path)) {
  readr::read_csv(ns_path, show_col_types = FALSE) |>
    dplyr::filter(rank == 1) |>
    dplyr::select(final_label, ns_community = community, ns_grank = grank)
} else {
  tibble::tibble(final_label = character(), ns_community = character(),
                 ns_grank = character())
}

# --- 1. Build the RAT (one row per class code) ---------------------------
# Parent-inheritance fallback: a split sub-class (e.g. S01.a) with no explicit
# class_categories row inherits its PARENT's (S01) moisture/elevation, so the
# categories table need not be re-curated every time the splits change. The
# color *category* is still derived per sub-class from its own indicator
# (below), so e.g. S01.c Artemisia still resolves to shrubland. Explicit
# sub-class rows in class_categories.csv override the inherited values.
cats_resolved <- tibble::tibble(final_label = lookup$final_label,
                                parent = sub("\\..*$", "", lookup$final_label)) |>
  dplyr::left_join(cats, by = "final_label") |>
  dplyr::left_join(cats |> dplyr::transmute(parent = final_label,
                                            moisture_p = moisture,
                                            elevation_p = elevation),
                   by = "parent") |>
  dplyr::transmute(final_label,
                   moisture  = dplyr::coalesce(moisture, moisture_p),
                   elevation = dplyr::coalesce(elevation, elevation_p))

rat <- lookup |>
  dplyr::left_join(csum, by = "final_label") |>
  dplyr::left_join(cats_resolved, by = "final_label") |>
  dplyr::left_join(ns_draft, by = "final_label") |>
  dplyr::transmute(
    value        = class_code,
    short_label  = dplyr::coalesce(short_label, final_label),
    description  = dplyr::coalesce(description, final_label),
    final_label  = final_label,
    class_type   = class_type,
    moisture     = moisture,
    elevation    = elevation,
    ns_community = ns_community,
    ns_grank     = ns_grank,
    # on-map label: class short label + provisional NatureServe community
    map_label    = dplyr::if_else(
      is.na(ns_community),
      dplyr::coalesce(short_label, final_label),
      paste0(dplyr::coalesce(short_label, final_label), " — ", ns_community))
  )
stopifnot(all(!is.na(rat$moisture)),
          all(!is.na(rat$elevation)))
cat(sprintf("RAT: %d classes\n  moisture x elevation cell counts:\n",
            nrow(rat)))
print(rat |> dplyr::count(moisture, elevation, name = "n") |>
        tidyr::pivot_wider(names_from = elevation, values_from = n,
                           values_fill = 0L) |>
        as.data.frame())

# --- 2. Physiognomy/habitat family + within-family HCL ramp --------------
# Graminoid genera (grasses, sedges, rushes) -> grassland; everything else
# herbaceous -> forb_meadow. Used only when class_type is meadow and the
# class is not wet.
gram_genera <- c("Festuca", "Leucopoa", "Deschampsia", "Calamagrostis",
  "Bromus", "Bromopsis", "Poa", "Elymus", "Pascopyrum", "Achnatherum",
  "Hesperostipa", "Stipa", "Danthonia", "Phleum", "Trisetum", "Koeleria",
  "Muhlenbergia", "Agrostis", "Vahlodea", "Helictotrichon",
  "Carex", "Kobresia", "Juncus", "Eleocharis", "Luzula")
# Woody dominants among the meadow-clustered classes (e.g. Vaccinium heath,
# Purshia/Chrysothamnus shrub-steppe) are physiognomically shrubland even
# though their class_type is "meadow".
shrub_genera <- c("Vaccinium", "Artemisia", "Purshia", "Chrysothamnus",
  "Ericameria", "Amelanchier", "Salix", "Ribes", "Dasiphora", "Betula",
  "Alnus", "Lonicera", "Prunus", "Juniperus", "Sambucus", "Symphoricarpos",
  "Cornus", "Holodiscus", "Rhus", "Shepherdia", "Pentaphylloides")
genus1 <- function(x) sub(" .*$", "", x)

# dominant indicator + snow-free DOY (for lifeform call + within-family order)
aux <- readr::read_csv("data/small_reference/label_community_names.csv",
                       show_col_types = FALSE) |>
  dplyr::select(final_label, top_indicator, snow_free_doy_mean)
rat <- rat |> dplyr::left_join(aux, by = "final_label")

# optional per-class override from class_categories.csv $category
override <- if ("category" %in% names(cats))
  stats::setNames(cats$category, cats$final_label) else character(0)
rat <- rat |> dplyr::mutate(category = dplyr::case_when(
  final_label %in% names(override) &
    !is.na(override[final_label])          ~ unname(override[final_label]),
  class_type == "shrub"                    ~ "shrubland",
  genus1(top_indicator) %in% shrub_genera  ~ "shrubland",
  moisture == "wet"                        ~ "wetland",
  genus1(top_indicator) %in% gram_genera   ~ "grassland",
  TRUE                                     ~ "forb_meadow"))

# HCL ramps per family: (hue, chroma, lightness) start -> end. Wetland hue
# 350 -> 372 wraps through red (hcl() takes hue mod 360).
fam_hcl <- list(
  shrubland   = list(h = c(255, 300), c = c(30, 46), l = c(46, 78)),
  grassland   = list(h = c(82,  48),  c = c(46, 62), l = c(83, 52)),
  forb_meadow = list(h = c(140, 100), c = c(38, 54), l = c(68, 40)),
  wetland     = list(h = c(350, 372), c = c(45, 66), l = c(73, 50)))
fam_colors <- function(cat, n) {
  p <- fam_hcl[[cat]]; f <- if (n == 1) 0.5 else seq(0, 1, length.out = n)
  grDevices::hcl(h = p$h[1] + f * (p$h[2] - p$h[1]),
                 c = p$c[1] + f * (p$c[2] - p$c[1]),
                 l = p$l[1] + f * (p$l[2] - p$l[1]))
}
elev_ord <- c(montane = 1L, subalpine = 2L, alpine = 3L, disturbed = 2L)
rat <- rat |>
  dplyr::group_by(category) |>
  dplyr::arrange(dplyr::coalesce(elev_ord[elevation], 2L),
                 dplyr::coalesce(snow_free_doy_mean, 0), final_label,
                 .by_group = TRUE) |>
  dplyr::mutate(color_hex =
                  fam_colors(dplyr::first(category), dplyr::n())[dplyr::row_number()]) |>
  dplyr::ungroup() |>
  dplyr::arrange(value) |>
  dplyr::select(-top_indicator, -snow_free_doy_mean)   # keep `category` in RAT
stopifnot(all(!is.na(rat$color_hex)), all(!is.na(rat$category)))

cat("\nClasses per physiognomy family:\n")
print(rat |> dplyr::count(category, name = "n") |> as.data.frame())

# --- 3. Apply to each class raster present ------------------------------
# Discover from the class COGs that exist (ALMO/CRBU/UPTA + e.g. CRBU_2018),
# so new year/domain inference runs are labeled automatically.
in_dir  <- "data/derived/aop_classified"
domains <- sub("_class_3m_v1\\.tif$", "",
               list.files(in_dir, pattern = "_class_3m_v1\\.tif$"))

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
              gsub("&", "&amp;", rat$map_label)))))
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
