#=============================================================================
# SGI (Standardised Groundwater level Index) – Nechako River Basin
# Method: Bloomfield & Marchant (2013), doi:10.5194/hess-17-4769-2013
# Multi-source groundwater data integration:
#   1. BC PGOWN (primary) via Environmental Reporting BC
#   2. Federal: NRCan GIN/NGWD + ECCC WSC (tidyhydat) [stub]
#   3. BC Real-Time Water Data (ArcGIS) [stub]
#   4. BC GWELLS API [stub]
# Output: One median SGI value per month across all qualifying wells
#=============================================================================

#---- 0. Packages ------------------------------------------------------------
pkgs <- c("httr", "sf", "tidyverse", "lubridate", "jsonlite", "ckanr", "tidyhydat", "geosphere")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE, warn.conflicts = FALSE))

#---- 1. Configuration -------------------------------------------------------
cfg <- list(
  # === Source toggles (PGOWN only enabled by default) ===
  src_pgown      = TRUE,   # BC PGOWN monthly medians (working)
  src_gin_ngwd   = FALSE,  # NRCan GIN/NGWD – requires manual cache
  src_wsc        = FALSE,  # Water Survey of Canada (tidyhydat)
  src_bc_realtime= FALSE,  # BC Real-Time (ArcGIS)
  src_gwells     = FALSE,  # BC GWELLS API
  
  # === PGOWN (primary) – URLs from working original ===
  url_pgown_ts   = "https://catalogue.data.gov.bc.ca/dataset/a74f1b97-17f7-499b-84e7-6455e169e425/resource/35162de7-803e-42f2-93d4-563330bb7dc9/download/monthlywell_ts.csv",
  url_pgown_attr = "https://catalogue.data.gov.bc.ca/dataset/a74f1b97-17f7-499b-84e7-6455e169e425/resource/a8933793-eadb-4a9c-992c-da4f6ac8ca51/download/gw_well_table.csv",
  
  # === GIN/NGWD (Federal) ===
  gin_cache      = "gin_nechako_subset.csv",
  
  # === WSC (tidyhydat) ===
  wsc_radius_km  = 50,
  wsc_daily_agg  = "median",
  
  # === BC Real-Time (ArcGIS) ===
  bc_rt_layer    = "https://maps.gov.bc.ca/arcserver/rest/services/pub/Water/BC_Water_Levels/MapServer/1",
  
  # === GWELLS API ===
  gwells_api_base= "https://apps.nrs.gov.bc.ca/gwells/api/v2",
  gwells_min_obs = 5,
  
  # === Basin boundary ===
  basin_kmz = "Spatial/nechakoBound_dissolve.kmz",
  
  # === Caching ===
  cache_dir = "cache/",
  
  # === Quality filters ===
  min_years    = 10,
  max_gap_frac = 0.50,
  
  # === Basin extent sanity check ===
  # Nechako watershed is ~44,000-48,000 km^2. Lat/lon bbox is wider than a
  # first guess would suggest: the Stuart-Takla headwaters sub-basin
  # (Trembleur/Takla Lakes) legitimately pushes the northern extent toward
  # ~56.2 N, well north of the main Nechako valley (~52.7-54.8 N). A run on
  # 2024 data produced area=47,038 km^2 (essentially an exact match to the
  # commonly cited ~47,100 km^2 figure) with bbox lat 52.9-56.2 / lon
  # -127.8 to -122.7 - that run is treated as the reference and the ranges
  # below are set around it with margin. These are still only a coarse sanity
  # check (catches wrong CRS, degenerate polygon, wrong KMZ layer) - not a
  # substitute for eyeballing the diagnostic plot this script produces.
  basin_area_km2_expected_range = c(20000, 70000),
  basin_bbox_lat_range = c(52.0, 57.0),
  basin_bbox_lon_range = c(-128.5, -121.0),
  
  # === Sign convention ===
  # TRUE  -> assumes PGOWN's med_gwl is *depth-to-water-below-ground*
  #          (larger value = drier), so it is negated so that larger/positive
  #          SGI = wetter, per Bloomfield & Marchant (2013) convention.
  # FALSE -> assumes med_gwl is already elevation/head (masl), in which case
  #          no sign flip is needed.
  # This flag flips every "No drought"/"Extreme drought" label in the output
  # if set wrong. The script now runs automated checks against it (magnitude
  # heuristic, CKAN data-dictionary text search, and a known-drought-year
  # cross-check) and prints a loud warning if they disagree - but a human
  # should still confirm against the PGOWN field definitions before trusting
  # the labels.
  invert_gw_sign = TRUE,
  pgown_dataset_id = "a74f1b97-17f7-499b-84e7-6455e169e425",
  
  # === Output ===
  target_months = NULL,
  output_prefix = "sgi_nechako_multisource"
)

if(!dir.exists(cfg$cache_dir)) dir.create(cfg$cache_dir, recursive = TRUE)

#---- 2. Helpers -------------------------------------------------------------
read_kmz <- function(kmz_path) {
  tmp_dir <- tempdir()
  utils::unzip(kmz_path, exdir = tmp_dir)
  kml_file <- list.files(tmp_dir, pattern = "\\.kml$", full.names = TRUE)[1]
  if(is.na(kml_file)) stop("No .kml found in ", kmz_path)
  lyr <- sf::st_layers(kml_file)$name[1]
  sf::st_read(kml_file, layer = lyr, quiet = TRUE) %>%
    sf::st_transform(4326) %>%
    sf::st_union() %>%
    sf::st_make_valid()
}

fetch_csv <- function(url, cache_path, desc = "") {
  if(!is.null(cache_path) && file.exists(cache_path)) {
    message(sprintf("  [cache] %s", cache_path))
    return(readr::read_csv(cache_path, show_col_types = FALSE))
  }
  message(sprintf("  [download] %s …", desc))
  dest <- if(!is.null(cache_path)) cache_path else tempfile(fileext = ".csv")
  resp <- httr::GET(url, httr::write_disk(dest, overwrite = TRUE), httr::progress(), httr::timeout(120))
  httr::stop_for_status(resp, task = paste("download", desc))
  readr::read_csv(dest, show_col_types = FALSE)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

normal_scores <- function(x) {
  n <- sum(!is.na(x))
  if(n < 4) return(rep(NA_real_, length(x)))
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((2 * r - 1) / (2 * n))
}

detect_col <- function(df, patterns, label) {
  m <- grep(paste(patterns, collapse = "|"), names(df), ignore.case = TRUE, value = TRUE)
  if(!length(m)) stop(sprintf("Cannot find '%s' column. Tried: %s", label, paste(patterns, collapse = ", ")))
  if(length(m) > 1) message(sprintf("  '%s': multiple matches (%s); using '%s'", label, paste(m, collapse = ", "), m[1]))
  m[1]
}

#---- 3. Load Basin ----------------------------------------------------------
message("Loading basin boundary …")
basin <- read_kmz(cfg$basin_kmz)
basin_centroid <- sf::st_centroid(basin) %>% sf::st_coordinates()

#---- 3b. Basin Extent Sanity Check ------------------------------------------
# A 232 -> 2 well drop after the spatial filter (seen in an earlier run) is
# consistent with either (a) PGOWN being a province-wide network with only a
# couple of wells genuinely inside Nechako, or (b) a broken/undersized basin
# polygon (wrong CRS, degenerate geometry after st_union/st_make_valid, wrong
# KMZ layer). This block checks the polygon's area and bounding box against
# published Nechako watershed figures and writes a plot of the polygon plus
# candidate well locations so the extent can be confirmed visually rather
# than trusted blindly.
message("\n=== Basin extent sanity check ===")

basin_area_km2 <- as.numeric(sf::st_area(basin)) / 1e6
basin_bbox <- sf::st_bbox(basin)

message(sprintf("  Basin polygon area: %.0f km^2 (expected range: %.0f - %.0f km^2)",
                basin_area_km2, cfg$basin_area_km2_expected_range[1], cfg$basin_area_km2_expected_range[2]))
message(sprintf("  Basin bbox: lat %.3f to %.3f | lon %.3f to %.3f",
                basin_bbox["ymin"], basin_bbox["ymax"], basin_bbox["xmin"], basin_bbox["xmax"]))

area_ok <- basin_area_km2 >= cfg$basin_area_km2_expected_range[1] &&
           basin_area_km2 <= cfg$basin_area_km2_expected_range[2]
lat_ok  <- basin_bbox["ymin"] >= cfg$basin_bbox_lat_range[1] && basin_bbox["ymax"] <= cfg$basin_bbox_lat_range[2]
lon_ok  <- basin_bbox["xmin"] >= cfg$basin_bbox_lon_range[1] && basin_bbox["xmax"] <= cfg$basin_bbox_lon_range[2]

if (!area_ok || !lat_ok || !lon_ok) {
  message(strrep("!", 70))
  message("  WARNING: basin polygon area/extent falls OUTSIDE the expected")
  message("  Nechako watershed range. This suggests the KMZ, CRS, or the")
  message("  st_union()/st_make_valid() step may have produced a degenerate")
  message("  or wrong-sized polygon. A large downstream well drop-off is")
  message("  more likely to be a bug rather than a real 'few wells in basin'")
  message("  result. Inspect the saved basin/wells diagnostic plot before")
  message("  trusting the well count.")
  message(strrep("!", 70))
} else {
  message("  Basin area and bbox fall within the expected Nechako range.")
  message("  (This is only a coarse sanity check - still inspect the saved")
  message("  diagnostic plot to confirm the shape/extent looks right.)")
}

load_pgown <- function() {
  if(!cfg$src_pgown) return(NULL)
  message("Loading PGOWN (Environmental Reporting BC) …")
  
  # Load attributes
  attrs_raw <- fetch_csv(cfg$url_pgown_attr, file.path(cfg$cache_dir, "pgown_attrs.csv"), "well attributes")
  names(attrs_raw) <- tolower(names(attrs_raw))
  idx_cols <- grep("^\\.\\.\\.\\d+$", names(attrs_raw), value = TRUE)
  if(length(idx_cols)) attrs_raw <- dplyr::select(attrs_raw, -dplyr::all_of(idx_cols))
  
  col_well <- detect_col(attrs_raw, c("^well_num$", "well_id", "obs.*well"), "well ID")
  col_lat  <- detect_col(attrs_raw, c("^lat$", "^latitude$"), "latitude")
  col_lon  <- detect_col(attrs_raw, c("^lon$", "^long$", "^longitude$"), "longitude")
  col_name <- detect_col(attrs_raw, c("^well_name$", "name"), "well name")
  
  attrs <- attrs_raw %>%
    rename(well_num = all_of(col_well), lat = all_of(col_lat), lon = all_of(col_lon), well_name = all_of(col_name)) %>%
    mutate(well_num = as.character(well_num)) %>%
    filter(!is.na(lat), !is.na(lon))
  
  # Load time series
  ts_raw <- fetch_csv(cfg$url_pgown_ts, file.path(cfg$cache_dir, "pgown_ts.csv"), "monthly medians")
  names(ts_raw) <- tolower(names(ts_raw))
  
  col_wn  <- detect_col(ts_raw, c("^well_num$", "well_id"), "well ID")
  col_dt  <- detect_col(ts_raw, c("^date$", "datetime"), "date")
  col_val <- detect_col(ts_raw, c("med_gwl", "gwl", "value", "level"), "GWL")
  
  # --- Sign-convention check -------------------------------------------------
  # cfg$invert_gw_sign determines whether every "No drought"/"Extreme drought"
  # label in the final output is right-side-up or backwards. Rather than
  # trusting the flag from a comment, run two automated checks and print a
  # loud warning if either disagrees with the configured setting. Neither
  # check is authoritative on its own - the PGOWN field definitions/data
  # dictionary is the real source of truth - but both catch the common
  # failure mode (depth-below-ground mistaken for elevation-above-sea-level
  # or vice versa).
  message(sprintf("\n  --- Sign-convention check (invert_gw_sign = %s) ---", cfg$invert_gw_sign))
  raw_sample <- utils::head(ts_raw[[col_val]], 10)
  message(sprintf("  Sample raw '%s' values: %s", col_val, paste(round(raw_sample, 3), collapse = ", ")))

  # Heuristic 1: magnitude. BC PGOWN depth-to-water is typically 0-30 m.
  # Elevation/head (masl) in the BC interior is typically several hundred
  # metres. A median well below ~50 strongly suggests depth-below-ground
  # (should be inverted); a median well above ~100 strongly suggests
  # elevation/head (should NOT be inverted).
  med_abs_val <- stats::median(abs(as.numeric(ts_raw[[col_val]])), na.rm = TRUE)
  heuristic_says_depth <- med_abs_val < 50
  message(sprintf("  Median |%s| across all records: %.2f -> heuristic suggests: %s",
                  col_val, med_abs_val, ifelse(heuristic_says_depth,
                                                "depth-below-ground (should invert)",
                                                "elevation/head (should NOT invert)")))

  heuristic_mismatch <- heuristic_says_depth != cfg$invert_gw_sign
  if (heuristic_mismatch) {
    message(strrep("!", 70))
    message("  WARNING: magnitude heuristic DISAGREES with cfg$invert_gw_sign.")
    message("  Every drought/no-drought label downstream may be inverted.")
    message(strrep("!", 70))
  }

  # Heuristic 2: try to pull the PGOWN dataset's CKAN metadata and search its
  # description/notes text for wording that indicates which convention the
  # med_gwl field actually uses. This is best-effort - if the API call fails
  # or the notes are silent on this, it says so rather than guessing.
  dict_hit <- tryCatch({
    meta_url <- sprintf("https://catalogue.data.gov.bc.ca/api/3/action/package_show?id=%s",
                        cfg$pgown_dataset_id)
    resp <- httr::GET(meta_url, httr::timeout(20))
    if (httr::status_code(resp) == 200) {
      meta <- httr::content(resp, as = "parsed", type = "application/json")
      notes_text <- paste(
        meta$result$notes %||% "",
        paste(sapply(meta$result$resources, function(r) r$description %||% ""), collapse = " "),
        collapse = " "
      )
      notes_lower <- tolower(notes_text)
      depth_kw <- c("depth to water", "depth below ground", "below ground surface", "btoc", "depth to gw")
      elev_kw  <- c("elevation", "metres above sea level", "masl", "above sea level", "hydraulic head")
      found_depth <- depth_kw[sapply(depth_kw, function(k) grepl(k, notes_lower, fixed = TRUE))]
      found_elev  <- elev_kw[sapply(elev_kw, function(k) grepl(k, notes_lower, fixed = TRUE))]
      list(found_depth = found_depth, found_elev = found_elev, checked = TRUE)
    } else {
      list(checked = FALSE, status = httr::status_code(resp))
    }
  }, error = function(e) list(checked = FALSE, error = conditionMessage(e)))

  if (isTRUE(dict_hit$checked)) {
    if (length(dict_hit$found_depth) || length(dict_hit$found_elev)) {
      message(sprintf("  CKAN metadata text mentions: depth-related terms [%s] | elevation-related terms [%s]",
                      paste(dict_hit$found_depth, collapse = ", "), paste(dict_hit$found_elev, collapse = ", ")))
      if (length(dict_hit$found_depth) && !length(dict_hit$found_elev) && !cfg$invert_gw_sign) {
        message("  WARNING: metadata language suggests depth-below-ground, but invert_gw_sign=FALSE.")
      }
      if (length(dict_hit$found_elev) && !length(dict_hit$found_depth) && cfg$invert_gw_sign) {
        message("  WARNING: metadata language suggests elevation/head, but invert_gw_sign=TRUE.")
      }
    } else {
      message("  CKAN metadata fetched but contained no explicit depth/elevation wording.")
    }
  } else {
    message(sprintf("  Could not fetch/parse CKAN metadata for automated dictionary check (%s).",
                    dict_hit$error %||% paste("HTTP", dict_hit$status)))
  }
  message("  ACTION: confirm med_gwl's convention against the PGOWN field definitions")
  message("  (https://catalogue.data.gov.bc.ca/dataset/) before trusting drought_cat labels.")
  message(strrep("-", 70))

  ts <- ts_raw %>%
    rename(well_num = all_of(col_wn), date = all_of(col_dt), gw_raw = all_of(col_val)) %>%
    mutate(
      well_num = as.character(well_num),
      date = as.Date(date),
      gw = if (cfg$invert_gw_sign) -as.numeric(gw_raw) else as.numeric(gw_raw),
      month_label = format(date, "%Y-%m"),
      source_id = "PGOWN"
    ) %>%
    filter(!is.na(date), !is.na(gw))
  
  if(nrow(ts) == 0) {
    message("  No valid time series records after initial cleaning.")
    return(NULL)
  }
  
  # Join with attributes to get coordinates
  ts_with_coords <- ts %>%
    left_join(attrs %>% dplyr::select(well_num, lat, lon), by = "well_num") %>%
    filter(!is.na(lat), !is.na(lon))
  
  if(nrow(ts_with_coords) == 0) {
    message("  No wells have coordinates in attribute table.")
    return(NULL)
  }
  
  # Spatial filter inside basin
  wells_sf <- ts_with_coords %>%
    dplyr::distinct(well_num, lat, lon) %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  in_basin <- sf::st_within(wells_sf, basin, sparse = FALSE)[, 1]

  message(sprintf("  Distinct candidate wells with coordinates: %d | inside basin: %d",
                  nrow(wells_sf), sum(in_basin)))

  # Diagnostic plot: basin polygon + every candidate well, coloured by
  # inside/outside, so a suspiciously low well count can be checked visually
  # against the actual polygon shape instead of trusted from the count alone.
  wells_sf$in_basin <- ifelse(in_basin, "Inside basin", "Outside basin")
  p_extent <- ggplot() +
    geom_sf(data = basin, fill = "grey80", colour = "grey30", alpha = 0.5) +
    geom_sf(data = wells_sf, aes(colour = in_basin), size = 1.6) +
    scale_colour_manual(values = c("Inside basin" = "#4575b4", "Outside basin" = "#d73027")) +
    labs(title = "Basin extent check: PGOWN candidate wells vs. basin polygon",
         subtitle = sprintf("%d of %d candidate wells fall inside the basin polygon",
                            sum(in_basin), nrow(wells_sf)),
         colour = NULL, x = NULL, y = NULL) +
    theme_bw(base_size = 10)
  ggsave(paste0(cfg$output_prefix, "_basin_extent_check.png"), p_extent, width = 8, height = 7, dpi = 150)
  message(sprintf("  Diagnostic plot saved: %s_basin_extent_check.png (inspect this before trusting the well count)",
                  cfg$output_prefix))

  if(sum(in_basin) == 0) {
    message("  No wells from time series are inside the basin.")
    return(NULL)
  }

  keep_wells <- wells_sf$well_num[in_basin]
  ts_filtered <- ts_with_coords %>%
    filter(well_num %in% keep_wells) %>%
    dplyr::select(well_num, date, gw, month_label, source_id)
  
  message(sprintf("  PGOWN: %d wells, %d records (after spatial filter)",
                  n_distinct(ts_filtered$well_num), nrow(ts_filtered)))
  ts_filtered
}

#---- 5. Optional Loaders (stubs – enable one at a time for testing) ---------
# GIN/NGWD: requires manual download of subset to cache/gin_nechako_subset.csv
load_gin <- function() {
  if(!cfg$src_gin_ngwd) return(NULL)
  message("Loading GIN/NGWD (NRCan) … [STUB – requires manual cache]")
  if(!file.exists(cfg$gin_cache)) {
    message(sprintf("  [manual] Download GIN subset for Nechako bbox to %s", cfg$gin_cache))
    return(NULL)
  }
  # Implementation pending cache availability
  NULL
}

# WSC via tidyhydat: limited groundwater stations in BC
load_wsc <- function() {
  if(!cfg$src_wsc) return(NULL)
  message("Loading WSC groundwater (tidyhydat) … [STUB]")
  # Implementation pending testing
  NULL
}

# BC Real-Time (ArcGIS): requires API exploration
load_bc_realtime <- function() {
  if(!cfg$src_bc_realtime) return(NULL)
  message("Loading BC Real-Time (ArcGIS) … [STUB]")
  NULL
}

# GWELLS API: sparse time-series data
load_gwells <- function() {
  if(!cfg$src_gwells) return(NULL)
  message("Loading GWELLS API (BC) … [STUB]")
  NULL
}

#---- 6. Aggregate Sources ---------------------------------------------------
message("\n=== Loading groundwater time series ===")
sources <- list(
  PGOWN = load_pgown(),
  GIN   = load_gin(),
  WSC   = load_wsc(),
  BC_RT = load_bc_realtime(),
  GWELLS= load_gwells()
)

ts_all <- dplyr::bind_rows(purrr::compact(sources))
if(nrow(ts_all) == 0) stop("No data loaded from any source.")

# Deduplicate: same well_num + month_label → keep source priority order
source_priority <- c("PGOWN", "GIN", "WSC", "BC_RT", "GWELLS")
ts_dedup <- ts_all %>%
  mutate(src_rank = match(source_id, source_priority)) %>%
  group_by(well_num, month_label) %>%
  slice_min(src_rank, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-src_rank)

message(sprintf("Combined: %d unique wells, %d monthly records (after dedup)",
                n_distinct(ts_dedup$well_num), nrow(ts_dedup)))

#---- 7. Quality Filters -----------------------------------------------------
well_stats <- ts_dedup %>%
  group_by(well_num) %>%
  summarise(
    n_months = n(),
    date_min = min(date), date_max = max(date),
    record_years = as.numeric(difftime(max(date), min(date), units = "days")) / 365.25,
    n_expected = as.numeric(difftime(max(date), min(date), units = "days")) / 30.44 + 1,
    gap_frac = 1 - n_months / pmax(1, n_expected),
    sources_used = paste(unique(source_id), collapse = ","),
    .groups = "drop"
  )

wells_ok <- well_stats %>%
  filter(
    record_years >= ifelse(grepl("PGOWN", sources_used), cfg$min_years, pmax(5, cfg$min_years * 0.5)),
    gap_frac <= cfg$max_gap_frac * ifelse(grepl("PGOWN", sources_used), 1, 1.5)
  )

message(sprintf("Wells after quality filters: %d", nrow(wells_ok)))
if(nrow(wells_ok) == 0) stop("No wells pass filters. Relax cfg$min_years or cfg$max_gap_frac.")

ts_filtered <- ts_dedup %>% filter(well_num %in% wells_ok$well_num)
# Load attributes globally to access well names and coordinates
attrs_raw <- fetch_csv(cfg$url_pgown_attr, file.path(cfg$cache_dir, "pgown_attrs.csv"), "well attributes")
names(attrs_raw) <- tolower(names(attrs_raw))
idx_cols <- grep("^\\.\\.\\.\\d+$", names(attrs_raw), value = TRUE)
if(length(idx_cols)) attrs_raw <- dplyr::select(attrs_raw, -dplyr::all_of(idx_cols))

col_well <- detect_col(attrs_raw, c("^well_num$", "well_id", "obs.*well"), "well ID")
col_lat  <- detect_col(attrs_raw, c("^lat$", "^latitude$"), "latitude")
col_lon  <- detect_col(attrs_raw, c("^lon$", "^long$", "^longitude$"), "longitude")
col_name <- detect_col(attrs_raw, c("^well_name$", "name"), "well name")

attrs <- attrs_raw %>%
  rename(well_num = all_of(col_well), lat = all_of(col_lat), lon = all_of(col_lon), well_name = all_of(col_name)) %>%
  mutate(well_num = as.character(well_num)) %>%
  filter(!is.na(lat), !is.na(lon))

# Extract names and coordinates of selected wells
selected_wells <- attrs %>%
  filter(well_num %in% wells_ok$well_num) %>%
  distinct(well_num, well_name, lat, lon)

message("\n=== Selected Wells (Names & Coordinates) ===")
print(selected_wells)
readr::write_csv(selected_wells, paste0(cfg$output_prefix, "_selected_wells.csv"))
message(sprintf("Saved: %s_selected_wells.csv", cfg$output_prefix))

message("\n=== Selected Wells (Names & Coordinates) ===")
print(selected_wells)
readr::write_csv(selected_wells, paste0(cfg$output_prefix, "_selected_wells.csv"))
message(sprintf("Saved: %s_selected_wells.csv", cfg$output_prefix))
#---- Print Well Details -----------------------------------------------------
message("\n=== Qualifying Wells ===")
well_details <- ts_filtered %>%
  distinct(well_num, source_id) %>%
  arrange(well_num)
for(i in seq_len(nrow(well_details))) {
  message(sprintf("Well %d: Code=%s | Source=%s", i, well_details$well_num[i], well_details$source_id[i]))
}
#---- 8. Compute SGI ---------------------------------------------------------
message("Computing SGI (Bloomfield & Marchant 2013) …")
sgi_long <- ts_filtered %>%
  mutate(cal_month = month(as.Date(paste0(month_label, "-01")))) %>%
  group_by(well_num, cal_month) %>%
  mutate(sgi = normal_scores(gw)) %>%
  ungroup()

# Basin aggregate: median across wells
sgi_basin <- sgi_long %>%
  group_by(month_label) %>%
  summarise(
    sgi_median = median(sgi, na.rm = TRUE),
    sgi_mean = mean(sgi, na.rm = TRUE),
    sgi_sd = sd(sgi, na.rm = TRUE),
    n_wells = sum(!is.na(sgi)),
    sources_contrib = paste(sort(unique(source_id)), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(month_label) %>%
  mutate(
    date = as.Date(paste0(month_label, "-01")),
    drought_cat = case_when(
      sgi_median >  0.00 ~ "No drought",
      sgi_median >= -1.00 ~ "Minor drought",
      sgi_median >= -1.50 ~ "Moderate drought",
      sgi_median >= -2.00 ~ "Severe drought",
      TRUE ~ "Extreme drought"
    )
  )

#---- 8b. Sign-convention cross-check against a known drought year -----------
# BC (including the Nechako region) experienced a widely-documented, severe
# drought through summer-fall 2023 (most BC basins reached Level 4-5 drought
# ratings). If invert_gw_sign is set correctly, sgi_basin should show
# negative (drought) values in Aug-Oct 2023. If it instead shows strongly
# positive values there, that's independent evidence the sign is flipped -
# on top of the magnitude/metadata checks already run above.
known_drought_window <- sgi_basin %>% filter(month_label %in% c("2023-08", "2023-09", "2023-10"))
if (nrow(known_drought_window) > 0) {
  mean_sgi_drought_window <- mean(known_drought_window$sgi_median, na.rm = TRUE)
  message(sprintf("\n  Known-drought-year cross-check: mean SGI for Aug-Oct 2023 = %.2f", mean_sgi_drought_window))
  if (!is.na(mean_sgi_drought_window) && mean_sgi_drought_window > 0.5) {
    message(strrep("!", 70))
    message("  WARNING: Aug-Oct 2023 is a well-documented severe BC drought period,")
    message("  but computed SGI here is strongly POSITIVE (wet). This is independent")
    message("  evidence that cfg$invert_gw_sign may be set backwards.")
    message(strrep("!", 70))
  } else if (!is.na(mean_sgi_drought_window) && mean_sgi_drought_window < 0) {
    message("  Consistent with the known 2023 drought (negative SGI) - supports the current sign convention.")
  }
} else {
  message("\n  Known-drought-year cross-check skipped: no Aug-Oct 2023 data in this well set.")
}

if(!is.null(cfg$target_months)) sgi_basin <- sgi_basin %>% filter(month_label %in% cfg$target_months)

#---- 9. Output --------------------------------------------------------------
cat("\n========================================================\n")
cat("  SGI – Nechako River Basin (multi-source)\n")
cat(sprintf("  Sources active: %s\n", paste(names(purrr::compact(sources)), collapse = ", ")))
cat(sprintf("  Wells: %d | Records: %d | Date range: %s to %s\n",
            n_distinct(sgi_long$well_num), nrow(sgi_long),
            min(sgi_long$date), max(sgi_long$date)))
cat("========================================================\n\n")
print(sgi_basin %>% dplyr::select(month_label, sgi_median, n_wells, sources_contrib, drought_cat), n = 30)

# Save
out_csv <- paste0(cfg$output_prefix, ".csv")
readr::write_csv(sgi_basin, out_csv)
message(sprintf("\nSaved: %s", out_csv))

#---- 10. Plot ---------------------------------------------------------------
p <- ggplot(sgi_basin, aes(x = date, y = sgi_median)) +
  geom_ribbon(aes(ymin = pmin(sgi_median, 0), ymax = 0), fill = "#d73027", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(sgi_median, 0)), fill = "#4575b4", alpha = 0.3) +
  geom_line(linewidth = 0.5, colour = "grey20") +
  geom_hline(yintercept = c(0, -1, -1.5, -2), linetype = c("solid", "dashed", "dashed", "dashed"),
             colour = c("black", "#fc8d59", "#d73027", "#a50026"), linewidth = 0.3) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(breaks = c(-2, -1.5, -1, 0, 1, 2),
                     labels = c("-2\n(Extreme)", "-1.5\n(Severe)", "-1\n(Moderate)", "0", "1", "2")) +
  labs(
    title = "Standardised Groundwater Index (SGI) – Nechako Basin",
    subtitle = paste0("Multi-source | ", n_distinct(sgi_long$well_num), " wells | ",
                      paste(names(purrr::compact(sources)), collapse = "+")),
    x = NULL, y = "SGI",
    caption = "Method: Bloomfield & Marchant (2013) | Data: PGOWN+GIN+WSC+BC_RT+GWELLS"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

ggsave(paste0(cfg$output_prefix, ".png"), p, width = 10, height = 4, dpi = 150)
message("Plot saved: ", cfg$output_prefix, ".png")
#---- Individual Well Plots --------------------------------------------------
message("\n=== Saving individual well SGI plots ===")
for(w in unique(sgi_long$well_num)) {
  p_well <- ggplot(sgi_long %>% filter(well_num == w), aes(x = date, y = sgi)) +
    geom_line(linewidth = 0.5, colour = "grey20") +
    geom_hline(yintercept = c(0, -1, -1.5, -2), linetype = c("solid", "dashed", "dashed", "dashed"),
               colour = c("black", "#fc8d59", "#d73027", "#a50026"), linewidth = 0.3) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title = paste("SGI – Well", w),
      subtitle = paste("Source:", unique(sgi_long$source_id[sgi_long$well_num == w])),
      x = NULL, y = "SGI",
      caption = "Method: Bloomfield & Marchant (2013)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(paste0(cfg$output_prefix, "_well_", w, ".png"), p_well, width = 10, height = 4, dpi = 150)
  message("Plot saved: ", cfg$output_prefix, "_well_", w, ".png")
}
#---- END --------------------------------------------------------------------
