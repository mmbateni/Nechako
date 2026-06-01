#=============================================================================
# SGI (Standardised Groundwater level Index) – Nechako River Basin
# Method: Bloomfield & Marchant (2013), doi:10.5194/hess-17-4769-2013
# Multi-source groundwater data integration:
#   1. BC PGOWN (primary) via Environmental Reporting BC
#   2. Federal: NRCan GIN/NGWD + ECCC WSC (tidyhydat) [stub]
#   3. BC Real-Time Water Data (ArcGIS) [stub]
#   4. BC GWELLS API [IMPLEMENTED]
# Spatial layers:
#   - Basin boundary (KMZ)
#   - BC Aquifer boundaries (WFS – WHSE_WATER_MANAGEMENT.GW_AQUIFERS_CLASSIFICATION_SVW)
# Output: Median SGI per month across all qualifying wells;
#         also stratified by aquifer unit
#=============================================================================

#---- 0. Packages ------------------------------------------------------------
pkgs <- c("httr", "sf", "tidyverse", "lubridate", "jsonlite",
          "ckanr", "tidyhydat", "geosphere")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)
invisible(lapply(pkgs, library, character.only = TRUE, warn.conflicts = FALSE))

#---- 1. Configuration -------------------------------------------------------
cfg <- list(
  # === Source toggles ===
  src_pgown       = TRUE,   # BC PGOWN monthly medians (working)
  src_gin_ngwd    = FALSE,  # NRCan GIN/NGWD – requires manual cache
  src_wsc         = FALSE,  # Water Survey of Canada (tidyhydat)
  src_bc_realtime = FALSE,  # BC Real-Time (ArcGIS)
  src_gwells      = TRUE,   # BC GWELLS API [NOW ENABLED]

  # === PGOWN (primary) ===
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
  # Endpoint returns paginated JSON; we filter spatially by basin bbox
  gwells_api_base  = "https://apps.nrs.gov.bc.ca/gwells/api/v2",
  gwells_min_obs   = 5,        # minimum water-level records per well
  gwells_page_size = 1000,     # results per page (API max)

  # === Aquifer boundaries (BC WFS) ===
  # WHSE_WATER_MANAGEMENT.GW_AQUIFERS_CLASSIFICATION_SVW served via OGC WFS
  # Filtered to basin bbox at fetch time; cache locally as GeoPackage
  aquifer_wfs_url  = paste0(
    "https://openmaps.gov.bc.ca/geo/pub/",
    "WHSE_WATER_MANAGEMENT.GW_AQUIFERS_CLASSIFICATION_SVW/wfs?",
    "service=WFS&version=2.0.0&request=GetFeature",
    "&typeName=pub:WHSE_WATER_MANAGEMENT.GW_AQUIFERS_CLASSIFICATION_SVW",
    "&outputFormat=application/json&srsName=EPSG:4326"
  ),
  aquifer_cache    = "cache/aquifers_nechako.gpkg",  # local GeoPackage cache

  # === Basin boundary ===
  basin_kmz = "Spatial/nechakoBound_dissolve.kmz",

  # === Caching ===
  cache_dir = "cache/",

  # === Quality filters ===
  min_years    = 10,
  max_gap_frac = 0.50,

  # === Output ===
  target_months = NULL,
  output_prefix = "sgi_nechako_multisource"
)

if (!dir.exists(cfg$cache_dir)) dir.create(cfg$cache_dir, recursive = TRUE)

#---- 2. Helpers -------------------------------------------------------------
read_kmz <- function(kmz_path) {
  tmp_dir  <- tempdir()
  utils::unzip(kmz_path, exdir = tmp_dir)
  kml_file <- list.files(tmp_dir, pattern = "\\.kml$", full.names = TRUE)[1]
  if (is.na(kml_file)) stop("No .kml found in ", kmz_path)
  lyr <- sf::st_layers(kml_file)$name[1]
  sf::st_read(kml_file, layer = lyr, quiet = TRUE) %>%
    sf::st_transform(4326) %>%
    sf::st_union() %>%
    sf::st_make_valid()
}

fetch_csv <- function(url, cache_path, desc = "") {
  if (!is.null(cache_path) && file.exists(cache_path)) {
    message(sprintf("  [cache] %s", cache_path))
    return(readr::read_csv(cache_path, show_col_types = FALSE))
  }
  message(sprintf("  [download] %s …", desc))
  dest <- if (!is.null(cache_path)) cache_path else tempfile(fileext = ".csv")
  resp <- httr::GET(url, httr::write_disk(dest, overwrite = TRUE),
                    httr::progress(), httr::timeout(120))
  httr::stop_for_status(resp, task = paste("download", desc))
  readr::read_csv(dest, show_col_types = FALSE)
}

normal_scores <- function(x) {
  n <- sum(!is.na(x))
  if (n < 4) return(rep(NA_real_, length(x)))
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((2 * r - 1) / (2 * n))
}

detect_col <- function(df, patterns, label) {
  m <- grep(paste(patterns, collapse = "|"), names(df),
            ignore.case = TRUE, value = TRUE)
  if (!length(m))
    stop(sprintf("Cannot find '%s' column. Tried: %s",
                 label, paste(patterns, collapse = ", ")))
  if (length(m) > 1)
    message(sprintf("  '%s': multiple matches (%s); using '%s'",
                    label, paste(m, collapse = ", "), m[1]))
  m[1]
}

#---- 3. Load Basin ----------------------------------------------------------
message("Loading basin boundary …")
basin          <- read_kmz(cfg$basin_kmz)
basin_bbox     <- sf::st_bbox(basin)          # xmin, ymin, xmax, ymax (WGS84)
basin_centroid <- sf::st_centroid(basin) %>% sf::st_coordinates()

#---- 4. Load Aquifer Boundaries (NEW) ---------------------------------------
load_aquifers <- function() {
  # Return cached GeoPackage if available
  if (file.exists(cfg$aquifer_cache)) {
    message(sprintf("  [cache] aquifer boundaries: %s", cfg$aquifer_cache))
    aq <- sf::st_read(cfg$aquifer_cache, quiet = TRUE)
    return(aq)
  }

  message("  [download] BC aquifer boundaries (WFS) …")

  # Build a CQL bbox filter so only aquifers overlapping the basin are fetched
  # CQL spatial filter: BBOX(GEOMETRY, minx, miny, maxx, maxy)
  bb  <- basin_bbox
  cql <- sprintf("BBOX(GEOMETRY,%f,%f,%f,%f)", bb["xmin"], bb["ymin"],
                 bb["xmax"], bb["ymax"])
  url <- paste0(cfg$aquifer_wfs_url, "&CQL_FILTER=", utils::URLencode(cql, reserved = TRUE))

  resp <- httr::GET(url, httr::timeout(120))
  httr::stop_for_status(resp, task = "fetch aquifer WFS")

  geojson_text <- httr::content(resp, as = "text", encoding = "UTF-8")
  aq_raw <- sf::st_read(geojson_text, quiet = TRUE) %>%
    sf::st_transform(4326) %>%
    sf::st_make_valid()

  if (nrow(aq_raw) == 0) {
    warning("  No aquifer polygons returned for basin bbox – aquifer grouping skipped.")
    return(NULL)
  }

  # Keep only polygons that actually intersect the basin (WFS bbox is a rough filter)
  in_basin_aq <- sf::st_intersects(aq_raw, basin, sparse = FALSE)[, 1]
  aq <- aq_raw[in_basin_aq, ]

  # Standardise key attribute names (WFS layer uses uppercase)
  names(aq) <- tolower(names(aq))
  # Typical fields: aquifer_id, aquifer_name, aquifer_subtype_code, vulnerability, type_description
  # Rename for convenience
  if ("aquifer_id" %in% names(aq))
    aq <- dplyr::rename(aq, aq_id = aquifer_id)
  if ("aquifer_name" %in% names(aq))
    aq <- dplyr::rename(aq, aq_name = aquifer_name)
  if ("aquifer_subtype_code" %in% names(aq))
    aq <- dplyr::rename(aq, aq_subtype = aquifer_subtype_code)
  if ("type_description" %in% names(aq))
    aq <- dplyr::rename(aq, aq_type_desc = type_description)

  message(sprintf("  Aquifers in basin: %d polygons", nrow(aq)))
  sf::st_write(aq, cfg$aquifer_cache, delete_dsn = TRUE, quiet = TRUE)
  aq
}

message("\nLoading aquifer boundaries …")
aquifers <- load_aquifers()   # sf or NULL

#---- 5. PGOWN Loader --------------------------------------------------------
load_pgown <- function() {
  if (!cfg$src_pgown) return(NULL)
  message("Loading PGOWN (Environmental Reporting BC) …")

  attrs_raw <- fetch_csv(cfg$url_pgown_attr,
                         file.path(cfg$cache_dir, "pgown_attrs.csv"),
                         "well attributes")
  names(attrs_raw) <- tolower(names(attrs_raw))
  idx_cols <- grep("^\\.\\.\\.\\d+$", names(attrs_raw), value = TRUE)
  if (length(idx_cols)) attrs_raw <- dplyr::select(attrs_raw, -dplyr::all_of(idx_cols))

  col_well <- detect_col(attrs_raw, c("^well_num$", "well_id", "obs.*well"), "well ID")
  col_lat  <- detect_col(attrs_raw, c("^lat$", "^latitude$"), "latitude")
  col_lon  <- detect_col(attrs_raw, c("^lon$", "^long$", "^longitude$"), "longitude")

  attrs <- attrs_raw %>%
    rename(well_num = all_of(col_well), lat = all_of(col_lat), lon = all_of(col_lon)) %>%
    mutate(well_num = as.character(well_num)) %>%
    filter(!is.na(lat), !is.na(lon))

  ts_raw <- fetch_csv(cfg$url_pgown_ts,
                      file.path(cfg$cache_dir, "pgown_ts.csv"),
                      "monthly medians")
  names(ts_raw) <- tolower(names(ts_raw))

  col_wn  <- detect_col(ts_raw, c("^well_num$", "well_id"), "well ID")
  col_dt  <- detect_col(ts_raw, c("^date$", "datetime"), "date")
  col_val <- detect_col(ts_raw, c("med_gwl", "gwl", "value", "level"), "GWL")

  ts <- ts_raw %>%
    rename(well_num = all_of(col_wn), date = all_of(col_dt), gw_raw = all_of(col_val)) %>%
    mutate(
      well_num    = as.character(well_num),
      date        = as.Date(date),
      gw          = -as.numeric(gw_raw),
      month_label = format(date, "%Y-%m"),
      source_id   = "PGOWN"
    ) %>%
    filter(!is.na(date), !is.na(gw))

  if (nrow(ts) == 0) { message("  No valid PGOWN records."); return(NULL) }

  ts_coords <- ts %>%
    left_join(attrs %>% dplyr::select(well_num, lat, lon), by = "well_num") %>%
    filter(!is.na(lat), !is.na(lon))

  if (nrow(ts_coords) == 0) { message("  No PGOWN wells with coordinates."); return(NULL) }

  wells_sf <- sf::st_as_sf(ts_coords, coords = c("lon", "lat"), crs = 4326)
  in_basin <- sf::st_within(wells_sf, basin, sparse = FALSE)[, 1]
  if (sum(in_basin) == 0) { message("  No PGOWN wells inside basin."); return(NULL) }

  out <- ts_coords[in_basin, ] %>%
    dplyr::select(well_num, date, gw, month_label, source_id)
  message(sprintf("  PGOWN: %d wells, %d records", n_distinct(out$well_num), nrow(out)))
  out
}

#---- 6. GWELLS Loader (IMPLEMENTED) -----------------------------------------
#
# Strategy:
#   (a) Fetch all wells inside the basin bounding box from the GWELLS wells
#       endpoint (paginated JSON).
#   (b) Spatially filter to the basin polygon.
#   (c) For each qualifying well, fetch its water-level records from the
#       /wells/{well_tag_number}/waterlevels/ endpoint.
#   (d) Aggregate to monthly medians and return in the common format.
#
# Notes:
#   • GWELLS water-level records are field measurements, not continuous
#     monitoring – coverage is typically sparse.  cfg$gwells_min_obs sets
#     the minimum number of observations required before a well is kept.
#   • PGOWN already provides continuous monthly medians for the observation
#     well network; GWELLS adds registered (non-PGOWN) wells that may have
#     periodic measurements.
#   • The function respects cfg$src_gwells = FALSE so it can be disabled.
#
load_gwells <- function() {
  if (!cfg$src_gwells) return(NULL)
  message("Loading GWELLS API (BC) …")

  bb <- basin_bbox   # named vector: xmin ymin xmax ymax

  #-- 6a. Fetch well list inside bbox (paginated) ----------------------------
  wells_cache <- file.path(cfg$cache_dir, "gwells_wells_nechako.rds")
  if (file.exists(wells_cache)) {
    message("  [cache] GWELLS well list")
    wells_df <- readRDS(wells_cache)
  } else {
    wells_list <- list()
    page       <- 1
    repeat {
      url <- sprintf(
        "%s/wells/?format=json&sw_lat=%f&sw_long=%f&ne_lat=%f&ne_long=%f&page=%d&page_size=%d",
        cfg$gwells_api_base,
        bb["ymin"], bb["xmin"], bb["ymax"], bb["xmax"],
        page, cfg$gwells_page_size
      )
      resp <- httr::GET(url, httr::timeout(60))
      if (httr::status_code(resp) != 200) {
        warning(sprintf("  GWELLS wells API returned HTTP %d on page %d – stopping.",
                        httr::status_code(resp), page))
        break
      }
      body <- jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"),
                                 flatten = TRUE)
      if (is.null(body$results) || length(body$results) == 0) break
      wells_list[[page]] <- as.data.frame(body$results)
      message(sprintf("    page %d: %d wells (total so far: %d / %s)",
                      page, nrow(wells_list[[page]]),
                      sum(sapply(wells_list, nrow)),
                      ifelse(is.null(body$count), "?", body$count)))
      if (is.null(body$`next`) || is.na(body$`next`)) break
      page <- page + 1
      Sys.sleep(0.3)   # be polite to the API
    }

    if (length(wells_list) == 0) {
      message("  No GWELLS wells found in basin bbox.")
      return(NULL)
    }
    wells_df <- dplyr::bind_rows(wells_list)
    saveRDS(wells_df, wells_cache)
  }

  # Harmonise coordinate column names (API returns latitude/longitude)
  names(wells_df) <- tolower(names(wells_df))
  if (!all(c("latitude", "longitude") %in% names(wells_df))) {
    warning("  GWELLS wells response missing latitude/longitude columns.")
    return(NULL)
  }
  wells_df <- wells_df %>%
    filter(!is.na(latitude), !is.na(longitude)) %>%
    mutate(latitude  = as.numeric(latitude),
           longitude = as.numeric(longitude))

  #-- 6b. Spatial filter to basin polygon ------------------------------------
  wells_sf   <- sf::st_as_sf(wells_df, coords = c("longitude", "latitude"), crs = 4326)
  in_basin   <- sf::st_within(wells_sf, basin, sparse = FALSE)[, 1]
  wells_basin <- wells_df[in_basin, ]

  message(sprintf("  GWELLS wells inside basin: %d", nrow(wells_basin)))
  if (nrow(wells_basin) == 0) return(NULL)

  # well_tag_number is the unique key for the water-levels endpoint
  if (!"well_tag_number" %in% names(wells_basin)) {
    warning("  GWELLS response has no 'well_tag_number' field.")
    return(NULL)
  }

  #-- 6c. Fetch water-level records per well ---------------------------------
  wl_cache <- file.path(cfg$cache_dir, "gwells_waterlevels_nechako.rds")
  if (file.exists(wl_cache)) {
    message("  [cache] GWELLS water-level records")
    wl_all <- readRDS(wl_cache)
  } else {
    wl_list <- list()
    n_wells <- nrow(wells_basin)
    for (i in seq_len(n_wells)) {
      wtn <- wells_basin$well_tag_number[i]
      url <- sprintf("%s/wells/%s/waterlevels/?format=json&page_size=10000",
                     cfg$gwells_api_base, wtn)
      resp <- httr::GET(url, httr::timeout(30))
      if (httr::status_code(resp) != 200) next
      body <- tryCatch(
        jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"),
                           flatten = TRUE),
        error = function(e) NULL
      )
      if (is.null(body) || is.null(body$results) || length(body$results) == 0) next
      df <- as.data.frame(body$results)
      df$well_tag_number <- wtn
      wl_list[[length(wl_list) + 1]] <- df
      if (i %% 20 == 0)
        message(sprintf("    fetched water levels: %d / %d wells", i, n_wells))
      Sys.sleep(0.2)
    }
    if (length(wl_list) == 0) {
      message("  No GWELLS water-level records found in basin.")
      return(NULL)
    }
    wl_all <- dplyr::bind_rows(wl_list)
    saveRDS(wl_all, wl_cache)
  }

  #-- 6d. Clean and aggregate to monthly medians ----------------------------
  names(wl_all) <- tolower(names(wl_all))

  # Identify date and depth columns
  date_col  <- grep("^date$|^surveyed_date$|^date_surveyed$|^measurement_date$",
                    names(wl_all), value = TRUE)[1]
  depth_col <- grep("^depth_below_ground_surface$|^static_level$|^water_depth$|^level$",
                    names(wl_all), value = TRUE)[1]

  if (is.na(date_col) || is.na(depth_col)) {
    warning(sprintf(
      "  GWELLS water-level response missing expected columns. Found: %s",
      paste(names(wl_all), collapse = ", ")))
    return(NULL)
  }

  wl_clean <- wl_all %>%
    rename(date_raw = all_of(date_col), depth_raw = all_of(depth_col)) %>%
    mutate(
      date       = as.Date(date_raw),
      gw_raw     = as.numeric(depth_raw),
      gw         = -gw_raw           # depth-to-water → higher = wetter, matches PGOWN sign
    ) %>%
    filter(!is.na(date), !is.na(gw)) %>%
    mutate(month_label = format(date, "%Y-%m"))

  # Drop wells with fewer than cfg$gwells_min_obs observations
  obs_per_well <- wl_clean %>%
    group_by(well_tag_number) %>%
    summarise(n_obs = n(), .groups = "drop") %>%
    filter(n_obs >= cfg$gwells_min_obs)

  wl_clean <- wl_clean %>% filter(well_tag_number %in% obs_per_well$well_tag_number)

  if (nrow(wl_clean) == 0) {
    message(sprintf("  No GWELLS wells pass the minimum %d observations filter.",
                    cfg$gwells_min_obs))
    return(NULL)
  }

  # Monthly median per well
  wl_monthly <- wl_clean %>%
    group_by(well_tag_number, month_label) %>%
    summarise(
      gw   = median(gw, na.rm = TRUE),
      date = min(date),
      .groups = "drop"
    ) %>%
    mutate(
      well_num  = as.character(well_tag_number),
      source_id = "GWELLS"
    ) %>%
    dplyr::select(well_num, date, gw, month_label, source_id)

  message(sprintf("  GWELLS: %d wells, %d monthly records (min %d obs filter applied)",
                  n_distinct(wl_monthly$well_num), nrow(wl_monthly),
                  cfg$gwells_min_obs))
  wl_monthly
}

#---- 7. Optional Loaders (stubs) --------------------------------------------
load_gin <- function() {
  if (!cfg$src_gin_ngwd) return(NULL)
  message("Loading GIN/NGWD (NRCan) … [STUB – requires manual cache]")
  if (!file.exists(cfg$gin_cache)) {
    message(sprintf("  [manual] Download GIN subset for Nechako bbox to %s", cfg$gin_cache))
    return(NULL)
  }
  NULL
}

load_wsc <- function() {
  if (!cfg$src_wsc) return(NULL)
  message("Loading WSC groundwater (tidyhydat) … [STUB]")
  NULL
}

load_bc_realtime <- function() {
  if (!cfg$src_bc_realtime) return(NULL)
  message("Loading BC Real-Time (ArcGIS) … [STUB]")
  NULL
}

#---- 8. Aggregate Sources ---------------------------------------------------
message("\n=== Loading groundwater time series ===")
sources <- list(
  PGOWN  = load_pgown(),
  GIN    = load_gin(),
  WSC    = load_wsc(),
  BC_RT  = load_bc_realtime(),
  GWELLS = load_gwells()
)

ts_all <- dplyr::bind_rows(purrr::compact(sources))
if (nrow(ts_all) == 0) stop("No data loaded from any source.")

# Deduplicate: same well_num + month → keep highest-priority source
source_priority <- c("PGOWN", "GIN", "WSC", "BC_RT", "GWELLS")
ts_dedup <- ts_all %>%
  mutate(src_rank = match(source_id, source_priority)) %>%
  group_by(well_num, month_label) %>%
  slice_min(src_rank, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-src_rank)

message(sprintf("Combined: %d unique wells, %d monthly records (after dedup)",
                n_distinct(ts_dedup$well_num), nrow(ts_dedup)))

#---- 9. Quality Filters -----------------------------------------------------
well_stats <- ts_dedup %>%
  group_by(well_num) %>%
  summarise(
    n_months     = n(),
    date_min     = min(date),
    date_max     = max(date),
    record_years = as.numeric(difftime(max(date), min(date), units = "days")) / 365.25,
    n_expected   = as.numeric(difftime(max(date), min(date), units = "days")) / 30.44 + 1,
    gap_frac     = 1 - n_months / pmax(1, n_expected),
    sources_used = paste(unique(source_id), collapse = ","),
    .groups = "drop"
  )

wells_ok <- well_stats %>%
  filter(
    record_years >= ifelse(grepl("PGOWN", sources_used),
                           cfg$min_years, pmax(5, cfg$min_years * 0.5)),
    gap_frac     <= cfg$max_gap_frac * ifelse(grepl("PGOWN", sources_used), 1, 1.5)
  )

message(sprintf("Wells after quality filters: %d", nrow(wells_ok)))
if (nrow(wells_ok) == 0)
  stop("No wells pass filters. Relax cfg$min_years or cfg$max_gap_frac.")

ts_filtered <- ts_dedup %>% filter(well_num %in% wells_ok$well_num)

#---- 10. Aquifer Spatial Join (NEW) -----------------------------------------
# Join each qualifying well to the aquifer polygon it falls within.
# Wells not inside any mapped aquifer get aq_id = NA ("Unmapped").
assign_aquifer <- function(ts_df, aq_sf) {
  if (is.null(aq_sf) || nrow(aq_sf) == 0) {
    message("  Aquifer layer unavailable – skipping aquifer assignment.")
    return(ts_df %>% mutate(aq_id = NA_integer_, aq_name = "Unmapped",
                            aq_subtype = NA_character_, aq_type_desc = NA_character_))
  }

  # Unique well coordinates from the time-series data
  # We need coordinates: re-join from the original source data
  # (ts_filtered only has well_num, date, gw, month_label, source_id)
  # Reconstruct from the deduped data's spatial context
  well_pts <- ts_df %>%
    distinct(well_num, source_id)

  # For PGOWN wells, retrieve coords from the cached attribute table
  pgown_attrs <- tryCatch(
    readr::read_csv(file.path(cfg$cache_dir, "pgown_attrs.csv"),
                    show_col_types = FALSE) %>%
      rename_with(tolower) %>%
      transmute(well_num = as.character(well_num),
                lat      = as.numeric(lat),
                lon      = as.numeric(lon)),
    error = function(e) tibble(well_num = character(), lat = numeric(), lon = numeric())
  )

  # For GWELLS wells, retrieve coords from cached well list
  gwells_locs <- tryCatch({
    wdf <- readRDS(file.path(cfg$cache_dir, "gwells_wells_nechako.rds"))
    wdf %>%
      mutate(well_num = as.character(well_tag_number),
             lat      = as.numeric(latitude),
             lon      = as.numeric(longitude)) %>%
      dplyr::select(well_num, lat, lon)
  }, error = function(e) tibble(well_num = character(), lat = numeric(), lon = numeric()))

  coords_all <- dplyr::bind_rows(pgown_attrs, gwells_locs) %>%
    filter(!is.na(lat), !is.na(lon)) %>%
    distinct(well_num, .keep_all = TRUE)

  well_coords <- well_pts %>%
    left_join(coords_all, by = "well_num") %>%
    filter(!is.na(lat), !is.na(lon))

  if (nrow(well_coords) == 0) {
    warning("  No well coordinates available for aquifer join.")
    return(ts_df %>% mutate(aq_id = NA_integer_, aq_name = "Unmapped",
                            aq_subtype = NA_character_, aq_type_desc = NA_character_))
  }

  wells_sf_aq <- sf::st_as_sf(well_coords, coords = c("lon", "lat"), crs = 4326)

  # Identify which fields to carry across from aquifer layer
  aq_cols <- intersect(c("aq_id", "aq_name", "aq_subtype", "aq_type_desc"), names(aq_sf))
  aq_lookup <- aq_sf %>% dplyr::select(dplyr::all_of(aq_cols))

  # Spatial join: each well gets the aquifer it falls within
  joined <- sf::st_join(wells_sf_aq, aq_lookup, join = sf::st_within, left = TRUE) %>%
    sf::st_drop_geometry()

  # A well can fall in >1 polygon (overlapping aquifer boundaries are rare but possible)
  # Keep first match per well
  joined_uniq <- joined %>%
    group_by(well_num) %>%
    slice(1) %>%
    ungroup()

  # Fill in missing aquifer fields
  for (col in c("aq_id", "aq_name", "aq_subtype", "aq_type_desc")) {
    if (!col %in% names(joined_uniq)) joined_uniq[[col]] <- NA
  }
  joined_uniq <- joined_uniq %>%
    mutate(aq_name = ifelse(is.na(aq_name), "Unmapped", aq_name))

  aq_assign <- joined_uniq %>% dplyr::select(well_num, aq_id, aq_name,
                                              aq_subtype, aq_type_desc)
  n_mapped <- sum(!is.na(aq_assign$aq_id))
  message(sprintf("  Wells mapped to an aquifer: %d / %d", n_mapped, nrow(aq_assign)))

  ts_df %>% left_join(aq_assign, by = "well_num")
}

message("\nAssigning wells to aquifer polygons …")
ts_filtered <- assign_aquifer(ts_filtered, aquifers)

#---- Print Well Details -----------------------------------------------------
message("\n=== Qualifying Wells ===")
well_details <- ts_filtered %>%
  distinct(well_num, source_id, aq_id, aq_name) %>%
  arrange(well_num)
for (i in seq_len(nrow(well_details))) {
  message(sprintf("Well %d: Code=%s | Source=%s | Aquifer=%s (ID=%s)",
                  i,
                  well_details$well_num[i],
                  well_details$source_id[i],
                  well_details$aq_name[i],
                  ifelse(is.na(well_details$aq_id[i]), "NA", well_details$aq_id[i])))
}

#---- 11. Compute SGI --------------------------------------------------------
message("Computing SGI (Bloomfield & Marchant 2013) …")
sgi_long <- ts_filtered %>%
  mutate(cal_month = month(as.Date(paste0(month_label, "-01")))) %>%
  group_by(well_num, cal_month) %>%
  mutate(sgi = normal_scores(gw)) %>%
  ungroup()

#-- 11a. Basin-wide aggregate ------------------------------------------------
sgi_basin <- sgi_long %>%
  group_by(month_label) %>%
  summarise(
    sgi_median     = median(sgi, na.rm = TRUE),
    sgi_mean       = mean(sgi, na.rm = TRUE),
    sgi_sd         = sd(sgi, na.rm = TRUE),
    n_wells        = sum(!is.na(sgi)),
    sources_contrib = paste(sort(unique(source_id)), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(month_label) %>%
  mutate(
    date       = as.Date(paste0(month_label, "-01")),
    drought_cat = case_when(
      sgi_median >  0.00 ~ "No drought",
      sgi_median >= -1.00 ~ "Minor drought",
      sgi_median >= -1.50 ~ "Moderate drought",
      sgi_median >= -2.00 ~ "Severe drought",
      TRUE                ~ "Extreme drought"
    )
  )

if (!is.null(cfg$target_months))
  sgi_basin <- sgi_basin %>% filter(month_label %in% cfg$target_months)

#-- 11b. Per-aquifer SGI aggregate (NEW) ------------------------------------
# Only computed where aquifer data is available
sgi_by_aquifer <- NULL
if (!all(is.na(sgi_long$aq_id))) {
  sgi_by_aquifer <- sgi_long %>%
    group_by(month_label, aq_id, aq_name) %>%
    summarise(
      sgi_median  = median(sgi, na.rm = TRUE),
      sgi_mean    = mean(sgi, na.rm = TRUE),
      sgi_sd      = sd(sgi, na.rm = TRUE),
      n_wells     = sum(!is.na(sgi)),
      .groups = "drop"
    ) %>%
    arrange(aq_id, month_label) %>%
    mutate(
      date       = as.Date(paste0(month_label, "-01")),
      drought_cat = case_when(
        sgi_median >  0.00 ~ "No drought",
        sgi_median >= -1.00 ~ "Minor drought",
        sgi_median >= -1.50 ~ "Moderate drought",
        sgi_median >= -2.00 ~ "Severe drought",
        TRUE                ~ "Extreme drought"
      )
    )

  if (!is.null(cfg$target_months))
    sgi_by_aquifer <- sgi_by_aquifer %>% filter(month_label %in% cfg$target_months)
}

#---- 12. Output -------------------------------------------------------------
cat("\n========================================================\n")
cat("  SGI – Nechako River Basin (multi-source)\n")
cat(sprintf("  Sources active: %s\n",
            paste(names(purrr::compact(sources)), collapse = ", ")))
cat(sprintf("  Wells: %d | Records: %d | Date range: %s to %s\n",
            n_distinct(sgi_long$well_num), nrow(sgi_long),
            min(sgi_long$date), max(sgi_long$date)))
if (!is.null(aquifers) && !is.null(sgi_by_aquifer)) {
  cat(sprintf("  Aquifers intersecting basin: %d | Aquifers with SGI data: %d\n",
              nrow(aquifers),
              n_distinct(sgi_by_aquifer$aq_id[!is.na(sgi_by_aquifer$aq_id)])))
}
cat("========================================================\n\n")

cat("--- Basin-wide SGI ---\n")
print(sgi_basin %>%
        dplyr::select(month_label, sgi_median, n_wells, sources_contrib, drought_cat),
      n = 30)

if (!is.null(sgi_by_aquifer)) {
  cat("\n--- Per-aquifer SGI (most recent 12 months) ---\n")
  recent_months <- tail(sort(unique(sgi_by_aquifer$month_label)), 12)
  print(sgi_by_aquifer %>%
          filter(month_label %in% recent_months, !is.na(aq_id)) %>%
          dplyr::select(month_label, aq_id, aq_name, sgi_median, n_wells, drought_cat) %>%
          arrange(month_label, aq_id),
        n = 60)
}

# Save CSVs
out_basin <- paste0(cfg$output_prefix, "_basin.csv")
readr::write_csv(sgi_basin, out_basin)
message(sprintf("\nSaved basin SGI: %s", out_basin))

if (!is.null(sgi_by_aquifer)) {
  out_aq <- paste0(cfg$output_prefix, "_by_aquifer.csv")
  readr::write_csv(sgi_by_aquifer, out_aq)
  message(sprintf("Saved per-aquifer SGI: %s", out_aq))
}

#---- 13. Plots --------------------------------------------------------------
#-- 13a. Basin-wide SGI ------------------------------------------------------
p_basin <- ggplot(sgi_basin, aes(x = date, y = sgi_median)) +
  geom_ribbon(aes(ymin = pmin(sgi_median, 0), ymax = 0),
              fill = "#d73027", alpha = 0.3) +
  geom_ribbon(aes(ymin = 0, ymax = pmax(sgi_median, 0)),
              fill = "#4575b4", alpha = 0.3) +
  geom_line(linewidth = 0.5, colour = "grey20") +
  geom_hline(yintercept = c(0, -1, -1.5, -2),
             linetype = c("solid", "dashed", "dashed", "dashed"),
             colour   = c("black", "#fc8d59", "#d73027", "#a50026"),
             linewidth = 0.3) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(
    breaks = c(-2, -1.5, -1, 0, 1, 2),
    labels = c("-2\n(Extreme)", "-1.5\n(Severe)", "-1\n(Moderate)", "0", "1", "2")
  ) +
  labs(
    title    = "Standardised Groundwater Index (SGI) – Nechako Basin",
    subtitle = paste0("Multi-source | ", n_distinct(sgi_long$well_num), " wells | ",
                      paste(names(purrr::compact(sources)), collapse = "+")),
    x = NULL, y = "SGI",
    caption  = "Method: Bloomfield & Marchant (2013) | Data: PGOWN+GWELLS"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(paste0(cfg$output_prefix, "_basin.png"), p_basin,
       width = 10, height = 4, dpi = 150)
message("Plot saved: ", cfg$output_prefix, "_basin.png")

#-- 13b. Per-aquifer SGI facet plot (NEW) ------------------------------------
if (!is.null(sgi_by_aquifer) && nrow(sgi_by_aquifer) > 0) {
  # Only plot aquifers that have a mapped ID and at least some data
  aq_plot_df <- sgi_by_aquifer %>%
    filter(!is.na(aq_id), !is.na(sgi_median)) %>%
    mutate(aq_label = paste0(aq_name, " (ID ", aq_id, ")"))

  if (nrow(aq_plot_df) > 0) {
    n_aq      <- n_distinct(aq_plot_df$aq_id)
    plot_h    <- max(3, 2 + n_aq * 1.6)

    p_aq <- ggplot(aq_plot_df, aes(x = date, y = sgi_median)) +
      geom_ribbon(aes(ymin = pmin(sgi_median, 0), ymax = 0),
                  fill = "#d73027", alpha = 0.35) +
      geom_ribbon(aes(ymin = 0, ymax = pmax(sgi_median, 0)),
                  fill = "#4575b4", alpha = 0.35) +
      geom_line(linewidth = 0.45, colour = "grey20") +
      geom_hline(yintercept = c(0, -1, -1.5, -2),
                 linetype   = c("solid", "dashed", "dashed", "dashed"),
                 colour     = c("black", "#fc8d59", "#d73027", "#a50026"),
                 linewidth  = 0.3) +
      facet_wrap(~aq_label, ncol = 1, scales = "free_x") +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      scale_y_continuous(breaks = c(-2, -1, 0, 1, 2)) +
      labs(
        title    = "SGI by Aquifer Unit – Nechako Basin",
        subtitle = paste0("BC Aquifer boundaries (GW_AQUIFERS_CLASSIFICATION_SVW) | ",
                          n_aq, " aquifer(s)"),
        x = NULL, y = "SGI",
        caption  = "Method: Bloomfield & Marchant (2013)"
      ) +
      theme_bw(base_size = 9) +
      theme(plot.title        = element_text(face = "bold"),
            panel.grid.minor  = element_blank(),
            strip.text        = element_text(face = "bold", size = 7))

    ggsave(paste0(cfg$output_prefix, "_by_aquifer.png"), p_aq,
           width = 10, height = plot_h, dpi = 150)
    message("Plot saved: ", cfg$output_prefix, "_by_aquifer.png")
  }
}

#-- 13c. Individual well SGI plots -------------------------------------------
message("\n=== Saving individual well SGI plots ===")
for (w in unique(sgi_long$well_num)) {
  wd         <- sgi_long %>% filter(well_num == w)
  aq_label   <- unique(wd$aq_name)[1]
  src_label  <- paste(unique(wd$source_id), collapse = "/")

  p_well <- ggplot(wd, aes(x = date, y = sgi)) +
    geom_line(linewidth = 0.5, colour = "grey20") +
    geom_hline(yintercept = c(0, -1, -1.5, -2),
               linetype   = c("solid", "dashed", "dashed", "dashed"),
               colour     = c("black", "#fc8d59", "#d73027", "#a50026"),
               linewidth  = 0.3) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    labs(
      title    = paste("SGI – Well", w),
      subtitle = paste0("Source: ", src_label,
                        " | Aquifer: ", ifelse(is.na(aq_label), "Unmapped", aq_label)),
      x = NULL, y = "SGI",
      caption  = "Method: Bloomfield & Marchant (2013)"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  ggsave(paste0(cfg$output_prefix, "_well_", w, ".png"), p_well,
         width = 10, height = 4, dpi = 150)
  message("Plot saved: ", cfg$output_prefix, "_well_", w, ".png")
}
#---- END --------------------------------------------------------------------
