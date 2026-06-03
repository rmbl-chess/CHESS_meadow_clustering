# 20_ecosystem_crosswalk.R — propose a NatureServe IVC community crosswalk for
# each meadow class by scoring its IndVal diagnostic assemblage against the
# cached Colorado community catalog (natureserve_fetch.py). Emits a top-3
# DRAFT for human curation; the curated result lives in
# data/small_reference/class_ecosystem_crosswalk.csv and feeds 19/06/10.
#
# Matching is transparent (no black box): per class x community we combine
#   - species F1: recall of the class's IndVal-weighted diagnostic species
#     captured by the community's name-species, x precision (fraction of the
#     community's named species the class explains). The IVC community NAME is
#     the diagnostic assemblage (Group-level compositions are empty), so we
#     match against name_species.
#   - ecology agreement: class (moisture, elevation) vs cues parsed from the
#     community name + summary (fen/marsh/wet -> wet; grassland/dry -> dry;
#     alpine/subalpine/montane -> elevation band).
#   Candidates are restricted to communities occurring in Colorado.
#
# Species names are reconciled to NatureServe nomenclature first (Veratrum
# tenuipetalum -> V. californicum, etc.) via species_natureserve_crosswalk.csv,
# so no class searches a name NatureServe doesn't index.
#
# Inputs:
#   data/derived/natureserve_cache.json                     (from fetch tool)
#   data/derived/label_descriptions.csv                     (IndVal assemblages)
#   data/small_reference/species_natureserve_crosswalk.csv  (reconciliation)
#   data/small_reference/class_categories.csv               (moisture x elevation)
# Output:
#   data/derived/natureserve_candidates.csv   (top-3 per class, with evidence)

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

cache <- jsonlite::fromJSON("data/derived/natureserve_cache.json",
                            simplifyVector = FALSE)
xwalk <- readr::read_csv("data/small_reference/species_natureserve_crosswalk.csv",
                         show_col_types = FALSE)
desc  <- readr::read_csv("data/derived/label_descriptions.csv",
                         show_col_types = FALSE)
cats  <- readr::read_csv("data/small_reference/class_categories.csv",
                         show_col_types = FALSE)

# --- reconciliation map (our_name -> list(name, level)) ------------------
# natureserve_name may hold multiple ";"-separated aliases (e.g. Veratrum
# californicum; Veratrum viride — the CO community uses viride).
recon <- setNames(
  Map(function(n, l) list(names = trimws(strsplit(n, ";")[[1]]), level = l),
      xwalk$natureserve_name, xwalk$match_level),
  xwalk$our_name)
reconcile <- function(sp) {
  if (!is.null(recon[[sp]])) recon[[sp]] else list(names = sp, level = "species")
}
genus_of <- function(x) sub(" .*$", "", x)

# --- class diagnostic assemblages: species + IndVal weight ---------------
# Parse "Festuca thurberi (cov=7.8%, freq=36%, IV=1.6); ..." -> name + IV.
parse_indicators <- function(s) {
  if (is.na(s) || !nzchar(s)) return(tibble(sp = character(), iv = numeric()))
  parts <- str_split(s, ";\\s*")[[1]]
  nm <- str_trim(str_replace(parts, "\\s*\\(.*$", ""))
  iv <- as.numeric(str_match(parts, "IV=([0-9.]+)")[, 2])
  tibble(sp = nm, iv = ifelse(is.na(iv), 0, iv)) |>
    filter(nzchar(sp), str_detect(sp, "^[A-Z][a-z]+ "))
}

# --- community catalog -> tidy frame (Colorado only) ---------------------
eco <- purrr::imap_dfr(cache$ecosystems, function(e, id) {
  tibble(
    id          = id,
    name        = e$scientificName %||% NA_character_,
    level       = e$classificationLevel %||% NA_character_,
    grank       = e$roundedGRank %||% NA_character_,
    in_co       = isTRUE(e$in_colorado),
    url         = e$nsxUrl %||% NA_character_,
    species     = list(unlist(e$name_species)),
    text        = paste(e$scientificName %||% "", e$summary %||% "",
                        e$conceptSentence %||% "")
  )
}) |> filter(in_co)

# Drop tree-overstory communities: every meadow/shrub class here is treeless,
# so Forest/Woodland units are physiognomic false matches (they get picked up
# via a shared understory forb, e.g. an aspen forest with Veratrum beneath).
n_all <- nrow(eco)
eco <- eco |> filter(!str_detect(name, "\\b(Forest|Woodland)\\b"))
cat(sprintf("Catalog: %d Colorado communities (%d Forest/Woodland dropped)\n",
            nrow(eco), n_all - nrow(eco)))

# moisture / elevation cues from community text
cue_moisture <- function(t) {
  t <- tolower(t)
  dplyr::case_when(
    str_detect(t, "fen|marsh|wet meadow|swamp|riparian|wetland|seep") ~ "wet",
    str_detect(t, "dry|xeric|sagebrush|shrub-steppe|grassland")        ~ "dry",
    TRUE                                                               ~ "mesic")
}
cue_elev <- function(t) {
  t <- tolower(t)
  dplyr::case_when(
    str_detect(t, "alpine")                 ~ "alpine",
    str_detect(t, "subalpine|upper montane") ~ "subalpine",
    str_detect(t, "montane|foothill|lower")  ~ "montane",
    TRUE                                     ~ NA_character_)
}
eco <- eco |> mutate(moist_cue = vapply(text, cue_moisture, ""),
                     elev_cue  = vapply(text, cue_elev, ""))

# --- score each class against the catalog --------------------------------
score_class <- function(lbl) {
  ind <- parse_indicators(desc$indicators[desc$final_label == lbl][1])
  if (nrow(ind) == 0) return(NULL)
  # reconcile class species to NatureServe names + level
  rc <- lapply(ind$sp, reconcile)
  ind$ns_names <- lapply(rc, `[[`, "names")
  ind$ns_level <- vapply(rc, `[[`, "", "level")
  tot_iv <- sum(ind$iv)
  cmoist <- cats$moisture[cats$final_label == lbl][1]
  celev  <- cats$elevation[cats$final_label == lbl][1]

  scored <- eco |> rowwise() |> mutate(
    hit = list({
      cs <- unlist(species)
      mapply(function(nms, lvl)
        if (lvl == "genus") any(genus_of(cs) %in% nms) else any(nms %in% cs),
        ind$ns_names, ind$ns_level)
    }),
    matched     = list(ind$sp[unlist(hit)]),
    # a genus-only token (e.g. native Bromopsis -> genus Bromus) must NOT
    # carry a candidate by itself -- require >=1 species-level hit. So a
    # genuinely B. inermis-heavy class (Bromopsis inermis -> Bromus inermis,
    # species-level) still matches its ruderal grassland, but native
    # Bromopsis sp. does not grab it.
    has_species = any(unlist(hit) & ind$ns_level == "species"),
    recall    = { h <- unlist(hit)
                  w <- ifelse(ind$ns_level == "genus", 0.5, 1.0)  # genus down-weighted
                  if (tot_iv > 0) sum(ind$iv[h] * w[h]) / tot_iv else 0 },
    precision = { cs <- unlist(species)
                  if (length(cs)) length(unlist(matched)) / length(cs) else 0 },
    f1        = if (recall + precision > 0)
                  2 * recall * precision / (recall + precision) else 0,
    eco_bonus = 0.15 * (!is.na(cmoist) & moist_cue == cmoist) +
                0.15 * (!is.na(elev_cue) & !is.na(celev) & elev_cue == celev),
    score     = f1 + eco_bonus
  ) |> ungroup() |>
    filter(has_species, lengths(matched) > 0) |>
    arrange(desc(score)) |>
    slice_head(n = 3)
  if (nrow(scored) == 0) return(NULL)
  scored |> transmute(
    final_label = lbl, rank = row_number(),
    community = name, ecosystem_id = id, level, grank, url,
    matched_species = vapply(matched, paste, "", collapse = "; "),
    recall = round(recall, 2), precision = round(precision, 2),
    score = round(score, 2), moist_cue, elev_cue,
    class_moisture = cmoist, class_elevation = celev)
}

labels <- sort(unique(desc$final_label))
candidates <- purrr::map_dfr(labels, score_class)
no_match <- setdiff(labels, unique(candidates$final_label))

readr::write_csv(candidates, "data/derived/natureserve_candidates.csv")
cat(sprintf("\nWrote data/derived/natureserve_candidates.csv: %d classes with candidates, %d with none\n",
            dplyr::n_distinct(candidates$final_label), length(no_match)))
if (length(no_match)) cat("  no candidate:", paste(no_match, collapse = ", "), "\n")
cat("\nTop candidate per class (preview):\n")
print(candidates |> filter(rank == 1) |>
        select(final_label, community, grank, score, recall, matched_species) |>
        as.data.frame(), row.names = FALSE)
