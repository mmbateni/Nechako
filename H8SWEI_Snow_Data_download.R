# ==============================================================================
# Nechako_Snow_Data_download.R  (revised acquisition pipeline)
# ==============================================================================
# This script combines two complementary Nechako Basin workflows:
#
#   PART A - GROUND STATION DATA
#     Screens BC Automated Snow Weather Stations (ASWS) against the Nechako
#     watershed boundary and compiles historical station SWE records into
#     daily.csv (station_id, lon, lat, date, swe_mm) for H8SWEI_Snow.R.
#
#   PART B - REVISED DATA ACQUISITION
#     Acquires the datasets used by the revised H8 pipeline: CanSWE, Sentinel-1
#     / C-SNOW / RCM SAR, CMC monthly SWE, GlobSnow SWE, CanSISE ensemble, and
#     OpenLandMap snow climatology. ERA5-Land SWE and AMSR SWE are assumed to
#     already exist locally and are not re-downloaded here.
#     NOTE: CMC (NSIDC-0447) does not produce a daily SWE product at all -
#     CMC only estimates SWE at monthly resolution (Oct-Jun each year), so
#     this uses CMC's monthly SWE file, not a daily one. See the CMC
#     acquisition block further down for details.
# ==============================================================================
#
library(sf)
library(terra)
library(dplyr)
library(readr)
library(lubridate)
library(appeears)
library(httr)
library(jsonlite)
library(dotenv)
library(ncdf4)    # PART C - reading Crocus-ERA5 NetCDF attributes/metadata
library(ggplot2)  # PART C - Crocus-ERA5 vs ERA5-Land sanity-check plot
library(parallel) # PART B - parallel decompression/extraction of cached GlobSnow files (ships with base R, no install needed)
# NOTE: 'elevatr' is NOT loaded here on purpose - it's only needed for the
# small number of stations that fall outside the local DEM mosaic's coverage
# (see the "FALLBACK: fetch elevation/aspect ... from a public source" block
# in Part A). It's loaded lazily via requireNamespace() there so a normal run
# where every station is covered by the local DEM never needs it installed.
# install.packages("elevatr") if you want the fallback available.

# ==============================================================================
# WARNING HANDLING
# ==============================================================================
# This session has options(warn = 2) active (set globally, via .Rprofile, or
# by another package loaded earlier) - it promotes EVERY R warning into a
# fatal error, including routine, expected ones this script is written to
# tolerate: readr/vroom's column-type-guessing notices, and sf's "attribute
# variables are assumed to be spatially constant throughout all geometries"
# on every st_intersection() call, to name the two already hit. Both of
# those are informational, not data problems - but under warn=2 either one
# stops the whole run. sf/terra/readr all emit warnings like this routinely,
# so rather than hunt down and suppressWarnings() each one individually as it
# surfaces, this restores R's normal warning behaviour (report, don't stop)
# for the rest of this session.
# If you specifically want warnings to stay fatal (e.g. a stricter CI-style
# run where you want to catch every single one), comment this out - but
# expect to need suppressWarnings() around most sf/terra/readr calls below.
old_warn_opt <- getOption("warn")
if (!is.null(old_warn_opt) && old_warn_opt >= 2) {
  message(
    "[OK] Detected options(warn = ", old_warn_opt, ") - downgrading to ",
    "options(warn = 1) for this script's run so routine sf/terra/readr ",
    "warnings are reported, not treated as fatal errors. (Was set before ",
    "this script started - check .Rprofile or an earlier-loaded package if ",
    "you didn't set this deliberately.)"
  )
}
# UPDATE: previously this only forced warn=1 when warn=2 was already active -
# under R's own DEFAULT (warn=0), every warning() call in this script (CMC
# auth failures, C-SNOW "directory not found"/"no files found", GlobSnow 404s,
# etc.) is buffered and only printed as a batch at the very end of the run -
# or not at all if the script errors out, is interrupted, or is run inside an
# environment (RStudio "Source", Rscript, a scheduled task) that doesn't
# flush/display the end-of-run warning summary. That's exactly what made
# several real failures above look like silent no-ops in the console log: the
# warning() explaining WHY was sitting in R's buffer, not printed at the point
# of failure. Setting warn=1 unconditionally (regardless of what it was set to
# before) makes every warning() print immediately, right where it happens, in
# real time - matching how message()/cat() already behave elsewhere in this
# script. This is strictly more visibility, never less: it does not change
# whether anything is treated as fatal (that's controlled by warn=2 vs
# anything below it, not by 0 vs 1), so it does not turn any of the ordinary
# sf/terra/readr informational warnings into stops.
options(warn = 1)

# ==============================================================================
# WORKING DIRECTORY + CREDENTIALS - set FIRST, before any part runs
# ==============================================================================
# Previously this was set partway through Part B's CHANGELOG section (see
# the old "dotenv::load_dot_env(...)" / "setwd(...)" lines that used to sit
# just above the CREDENTIALS block below), which meant PART A RAN - AND
# WROTE daily.csv - USING WHATEVER DIRECTORY R HAPPENED TO ALREADY BE IN,
# not this project folder. If that ambient directory ever differed from
# where H8SWEI_Snow.R expects daily.csv (its own working directory), Part
# A's output could silently miss H8 entirely, even though both scripts
# "looked" consistent on paper.
#
# Moved here so every part of this script, including Part A, runs from the
# same folder. IMPORTANT: this path must match H8SWEI_Snow.R's own setwd()
# call EXACTLY - both scripts must point at the same folder for daily.csv,
# the DEM, and the basin KMZ to be found by both.
# ==============================================================================
# COMMON PATHS
# ==============================================================================
setwd("D:/Nechako_Drought/Nechako")
SNOW_DATA_DIR <- "snow_data"
DAILY_CSV_PATH <- file.path(SNOW_DATA_DIR, "daily.csv")
dir.create(SNOW_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dotenv::load_dot_env(".env")  # EARTHDATA_USER / EARTHDATA_PASS - used by Part B only
# Also optional in .env: CDSE_CLIENT_ID / CDSE_CLIENT_SECRET - an OAuth
# client created under a Copernicus Data Space Ecosystem / Sentinel Hub
# account (https://shapps.dataspace.copernicus.eu/dashboard/#/account/settings
# -> "OAuth clients" -> Create), used only as a fallback Sentinel-1 source in
# acquire_sentinel1_asf() when ASF's search API fails or returns nothing. Not
# required if you're fine with Sentinel-1 acquisition being skipped on an ASF
# outage. NOTE: this is NOT the same login as an ECMWF/Climate Data Store
# account (cds.climate.copernicus.eu) - that's a different service despite
# the similar name and won't work here.
# Create a scratch folder on D: (inside or near your project) and point terra at it
scratch_path <- "D:/Nechako_Drought/Nechako/terra_scratch"
dir.create(scratch_path, recursive = TRUE, showWarnings = FALSE)
terra::terraOptions(tempdir = normalizePath(scratch_path))
Sys.setenv(TMPDIR = "D:/Nechako_Drought/Nechako/terra_scratch")
# ==============================================================================
# PART TOGGLES - run/skip whole parts without hand-commenting code out
# ==============================================================================
# H8SWEI_Snow.R's only HARD dependency on this script is Part A (daily.csv,
# used for basin-wide bias correction). Nechako_Snow_Bridge.R produce optional
# cross-validation series and can take hours (AppEEARS processing + large
# NetCDF downloads/prep) - see the "run as a background batch job" notes in
# the Part B and Part C headers below. Defaults (TRUE) preserve the original
# "run everything" behaviour. Set RUN_SATELLITE_AND_CROCUS_PARTS <- FALSE
# for a fast run that only produces what H8SWEI_Snow.R actually needs.
#
# Both can also be supplied as environment variables (e.g. from a batch/
# shell wrapper script that chains this script and H8SWEI_Snow.R together)
# without editing this file - see run_nechako_pipeline.bat / .sh.
RUN_PART_A                     <- as.logical(Sys.getenv("RUN_PART_A", "TRUE"))
RUN_SATELLITE_AND_CROCUS_PARTS <- as.logical(Sys.getenv("RUN_SATELLITE_AND_CROCUS_PARTS", "TRUE"))
if (is.na(RUN_PART_A)) RUN_PART_A <- TRUE
if (is.na(RUN_SATELLITE_AND_CROCUS_PARTS)) RUN_SATELLITE_AND_CROCUS_PARTS <- TRUE
cat(sprintf("[OK] RUN_PART_A = %s, RUN_SATELLITE_AND_CROCUS_PARTS (B+C+D) = %s\n",
            RUN_PART_A, RUN_SATELLITE_AND_CROCUS_PARTS))

# ==============================================================================
# PART A - GROUND STATION DATA: SPATIAL/TOPO SCREENING + SWE TIME SERIES
# ==============================================================================
if (RUN_PART_A) {
  
  # ------------------------------------------------------------------------------
  # Helper: read a boundary from either .kml or .kmz
  # ------------------------------------------------------------------------------
  # GDAL's KML driver (used by sf::st_read) cannot open a .kmz file directly -
  # a .kmz is just a zip archive with a .kml (usually "doc.kml") inside it, and
  # passing the .kmz path straight to st_read() fails with an error like:
  #   "Cannot open ... ; The source could be corrupt or not supported."
  # This helper reads the .kml straight out of the zip via GDAL's /vsizip/
  # virtual filesystem (no manual unzip needed), and falls back to physically
  # extracting it to a temp folder if /vsizip/ isn't available in this GDAL build.
  read_kml_or_kmz <- function(path) {
    if (!grepl("\\.kmz$", path, ignore.case = TRUE)) {
      return(st_read(path, quiet = TRUE))
    }
    
    kml_inside <- utils::unzip(path, list = TRUE)$Name
    kml_inside <- kml_inside[grepl("\\.kml$", kml_inside, ignore.case = TRUE)]
    if (length(kml_inside) == 0) stop("No .kml file found inside ", path)
    kml_inside <- kml_inside[1]
    
    vsi_path <- file.path("/vsizip", path, kml_inside)
    result <- tryCatch(st_read(vsi_path, quiet = TRUE), error = function(e) NULL)
    if (!is.null(result)) return(result)
    
    # Fallback: physically unzip to a temp dir, then read the extracted .kml
    tmp_dir <- file.path(tempdir(), "kmz_extract")
    dir.create(tmp_dir, showWarnings = FALSE)
    utils::unzip(path, files = kml_inside, exdir = tmp_dir)
    st_read(file.path(tmp_dir, kml_inside), quiet = TRUE)
  }
  
  # ------------------------------------------------------------------------------
  # 1. Spatial & Topographical Screening
  # ------------------------------------------------------------------------------
  
  # Load the Nechako watershed boundary
  # NOTE: update this path to wherever you saved the boundary file. Your upload
  # was "Nechako_Basin.kmz" - adjust the filename below if yours differs.
  nechako_bound_path <- "Spatial/Nechako_Basin.kmz"
  nechako_bound <- read_kml_or_kmz(nechako_bound_path)
  
  # Download and load the BC ASWS Stations
  # NOTE: the "openmaps.gov.bc.ca/kml/geo/layers/..." endpoint used previously
  # is a Google-Earth PREVIEW file - it only contains two <GroundOverlay> raster
  # images (a rendered map + legend), not <Placemark> point geometries. That's
  # why st_read() on it parsed successfully but reported "0 features" - there
  # was nothing to parse.
  # The actual vector data comes from BC's WFS (Web Feature Service) instead.
  # This queries the DataBC WFS directly and returns real station point
  # geometries as GeoJSON, which st_read() can read like any other file/URL.
  stations_wfs_url <- paste0(
    "https://openmaps.gov.bc.ca/geo/pub/wfs?",
    "SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature",
    "&typeName=pub:WHSE_WATER_MANAGEMENT.SSL_SNOW_ASWS_STNS_SP",
    "&outputFormat=json&SRSNAME=EPSG:4326"
  )
  stations <- st_read(stations_wfs_url, quiet = TRUE)
  
  if (nrow(stations) == 0) {
    stop(
      "WFS query returned 0 station features. Check the layer name is still\n",
      "current at https://catalogue.data.gov.bc.ca/dataset/automated-snow-weather-station-locations\n",
      "or try the bcdata R package (install.packages('bcdata')) instead:\n",
      "  library(bcdata)\n",
      "  stations <- bcdc_query_geodata('automated-snow-weather-station-locations') %>% collect()"
    )
  }
  
  # ------------------------------------------------------------------------------
  # Identify which WFS attribute actually holds the station ID (Pillow_ID)
  # ------------------------------------------------------------------------------
  # An earlier version of this script assumed the station ID lived in a column
  # called "Name" - that assumption dates back to an older KML-based workflow,
  # where GDAL's KML driver auto-synthesizes generic "Name"/"description"
  # columns from <Placemark><name>. Now that stations come from the WFS as
  # GeoJSON, the real BCGW attribute schema is exposed instead, and there is no
  # guarantee a "Name" column exists or is populated - it silently existed but
  # was blank for every feature, so every station "matched" on empty string
  # instead of raising an error, and every join against the raw archive failed
  # with 0 rows.
  # Instead of guessing the column name, this finds it: load the raw ASWS
  # archive now (its Pillow_ID values, e.g. "_1A01P", are known-good), then
  # pick whichever character column in `stations` has the most overlap with
  # those IDs once both sides are normalized (leading "_" stripped, upper-
  # cased, whitespace-trimmed).
  normalize_stn_id <- function(x) toupper(trimws(sub("^_", "", as.character(x))))
  
  # IMPORTANT: this raw archive file must NOT be named "daily.csv" - this
  # script's own OUTPUT is written to "daily.csv" at the end, and if the raw
  # archive used the same name in the same working directory, the write step
  # would overwrite your only copy of the source data. Rename your downloaded
  # archive (or point this path at wherever you keep it) before running.
  RAW_ASWS_ARCHIVE_PATH <- "asws_province_daily_raw.csv"
  
  if (!file.exists(RAW_ASWS_ARCHIVE_PATH)) {
    stop(
      "Raw ASWS archive not found at '", RAW_ASWS_ARCHIVE_PATH, "'. Download ",
      "it from https://catalogue.data.gov.bc.ca/dataset/archive-automated-snow-weather-station-data ",
      "(historical) and/or https://catalogue.data.gov.bc.ca/dataset/current-season-automated-snow-weather-station-data ",
      "(current water year), save it under this name, and re-run."
    )
  }
  
  # NOTE: read_csv()/vroom guesses each column's type from only the first
  # 1000 rows by default. This archive is a province-wide, multi-year dump
  # that gets appended to every season, so a type that looked consistent in
  # the first 1000 rows (e.g. "value" reading as clean numeric, or "Date"
  # reading as one consistent format) can easily be contradicted further
  # down the file - e.g. a "T" (trace) or "M" (missing) quality-flag string
  # in `value`, or a stray disclaimer/footer row appended by the export tool.
  # That produces per-cell parsing failures which read_csv normally reports
  # as a warning - but if this session has options(warn = 2) set (globally,
  # via .Rprofile, or from another package), that warning gets promoted to
  # the fatal "(converted from warning) One or more parsing issues" error
  # that stops the whole script before it ever reaches Part A's later steps.
  #
  # Fix: (1) guess_max = Inf forces column-type guessing to scan the WHOLE
  # file instead of just the first 1000 rows, so a type is far less likely
  # to be misguessed in the first place; (2) explicit col_types keeps
  # Pillow_ID/Date/value as character on the way in, so no numeric/date
  # parser can fail (and throw) during the read itself - date/numeric
  # conversion happens afterwards, in mutate(), where individual bad values
  # become NA (dropped by the is.finite(swe_mm) filter below) instead of
  # halting the whole read.
  raw_archive_raw <- read_csv(
    RAW_ASWS_ARCHIVE_PATH,
    col_types = cols(.default = col_character()),
    guess_max = Inf
  )
  
  parse_problems <- problems(raw_archive_raw)
  if (nrow(parse_problems) > 0) {
    warning(
      nrow(parse_problems), " parsing issue(s) found in ", RAW_ASWS_ARCHIVE_PATH,
      " (shown below) - reading all columns as character avoided a hard stop, ",
      "but it's worth checking these rows aren't silently losing data:\n",
      paste(utils::capture.output(print(parse_problems, n = 20)), collapse = "\n")
    )
  }
  
  raw_archive <- raw_archive_raw %>%
    filter(variable == "SWE") %>%
    mutate(
      # normalize IDs: strip a leading underscore, uppercase, trim whitespace,
      # so "_1D06P", "1d06p", and " 1D06P " all match the WFS station-ID
      # column consistently, whichever form either source happens to use
      pillow_id_norm = normalize_stn_id(Pillow_ID),
      date = as.Date(Date),
      swe_mm = suppressWarnings(as.numeric(value))
    ) %>%
    filter(is.finite(swe_mm)) %>%
    select(pillow_id_norm, date, swe_mm)
  
  n_dropped_nonnumeric <- sum(
    !is.na(raw_archive_raw$value[raw_archive_raw$variable == "SWE"]) &
      is.na(suppressWarnings(as.numeric(raw_archive_raw$value[raw_archive_raw$variable == "SWE"])))
  )
  if (n_dropped_nonnumeric > 0) {
    message(
      n_dropped_nonnumeric, " SWE row(s) had a non-numeric `value` (e.g. a ",
      "quality flag like \"T\"/\"M\" rather than a number) and were dropped ",
      "by the is.finite(swe_mm) filter above."
    )
  }
  
  archive_ids <- unique(raw_archive$pillow_id_norm)
  
  geom_col <- attr(stations, "sf_column")
  candidate_cols <- names(stations)[
    vapply(stations, function(col) is.character(col) || is.factor(col), logical(1))
  ]
  candidate_cols <- setdiff(candidate_cols, geom_col)
  
  col_match_counts <- vapply(
    candidate_cols,
    function(col_name) length(intersect(normalize_stn_id(stations[[col_name]]), archive_ids)),
    integer(1)
  )
  
  message(
    "WFS station-ID column detection - matches against ", length(archive_ids),
    " known Pillow_IDs from the raw archive: ",
    paste(sprintf("%s (%d)", candidate_cols, col_match_counts), collapse = ", ")
  )
  
  if (length(col_match_counts) == 0 || max(col_match_counts) == 0) {
    stop(
      "Could not find a WFS attribute column whose values match any Pillow_ID ",
      "in ", RAW_ASWS_ARCHIVE_PATH, ". Columns available on `stations`: ",
      paste(names(stations), collapse = ", "),
      ". Inspect these against the raw archive's Pillow_ID values and set ",
      "the join column manually."
    )
  }
  
  station_id_col <- names(which.max(col_match_counts))
  message(
    "Using WFS column '", station_id_col, "' as the station ID (",
    max(col_match_counts), "/", length(archive_ids), " known Pillow_IDs matched)."
  )
  stations$station_id <- as.character(stations[[station_id_col]])
  stations$station_id_norm <- normalize_stn_id(stations$station_id)
  
  # Create a spatial buffer (100 km - 100% larger than the original 50 km) around
  # the Nechako boundary
  nechako_buffer <- st_buffer(nechako_bound, dist = 100000)
  
  # Filter for candidate stations located within the buffer
  # NOTE: st_intersection() against a MULTI-FEATURE buffer (nechako_bound, and
  # therefore nechako_buffer, can have more than one polygon part/row - see
  # the nlyr(basin_elev) > 1 handling further down, which already anticipates
  # this for the DEM mask) produces ONE ROW PER (station, buffer-part) pair
  # for any station overlapping more than one part - silently duplicating
  # that station in stations_near. This is exactly what produced
  # "4A12P, 4A12P, 4A13P, 4A13P, 4A12P, 4A13P" in one run: 2 real stations,
  # tripled each because the buffer had 3 parts. Dissolve the buffer into a
  # single geometry first so each station appears at most once, regardless of
  # how many parts the basin polygon has.
  # st_union() applies stricter s2 geometry-validity checks than the plain
  # st_intersection() this replaced - it can reject a degenerate/duplicate
  # vertex (e.g. "Edge 936 is degenerate (duplicate vertex)") in a
  # KML-derived polygon that st_buffer() itself passed through without
  # complaint. Repair with st_make_valid() first; if s2 still rejects it,
  # fall back to unioning with s2 turned off (GEOS/planar), which tolerates
  # this kind of near-degenerate geometry better - then restore whatever the
  # s2 setting was before this block regardless of which path was taken.
  nechako_buffer <- st_make_valid(nechako_buffer)
  nechako_buffer_dissolved <- tryCatch(
    st_union(nechako_buffer),
    error = function(e) {
      message(
        "st_union() failed under s2 (", conditionMessage(e), ") - retrying ",
        "with s2 turned off (GEOS/planar union) instead."
      )
      old_s2 <- sf_use_s2()
      sf_use_s2(FALSE)
      on.exit(sf_use_s2(old_s2), add = TRUE)
      st_union(st_make_valid(nechako_buffer))
    }
  )
  # st_intersection() always emits "attribute variables are assumed to be
  # spatially constant throughout all geometries" whenever the left-hand
  # sf object carries non-geometry columns (station name/ID/etc., as
  # `stations` does here) - it's sf telling you those attributes get
  # copied through unchanged onto the clipped geometry, which is exactly
  # what's wanted for point-in-polygon clipping like this. It's routine,
  # not a data problem, but this session has options(warn = 2) (or
  # equivalent) active - see the read_csv() fix above for the same issue -
  # so without suppressWarnings() here it stops the whole script.
  stations_near <- suppressWarnings(st_intersection(stations, nechako_buffer_dissolved))
  
  # Safety net: dedupe by normalized station ID in case duplicates come from
  # somewhere else entirely (e.g. the WFS itself returning duplicate rows for
  # a station) - dissolving the buffer above should already prevent the
  # geometry-driven case, but this makes the "one row per station" invariant
  # hold regardless of the cause.
  n_dupe_stations <- sum(duplicated(stations_near$station_id_norm))
  if (n_dupe_stations > 0) {
    warning(
      n_dupe_stations, " station(s) still appeared more than once in ",
      "stations_near even after dissolving the buffer - deduplicating by ",
      "station_id_norm (keeping the first match). This shouldn't normally ",
      "happen post-dissolve; if you see this warning, it's worth checking ",
      "`stations` directly for literal duplicate rows returned by the WFS."
    )
    stations_near <- stations_near[!duplicated(stations_near$station_id_norm), ]
  }
  
  # ------------------------------------------------------------------------------
  # Helper: sanity-check a DEM/elevation raster before it's trusted for screening
  # ------------------------------------------------------------------------------
  # Added after a run where the mosaicked CDEM DEM silently produced
  # elevation = 0 for every single candidate station, which then trivially
  # passed the +/-300 m window (basin_mean_elev was ALSO ~0 from the same
  # underlying raster) - the elevation/aspect screen reported "0 excluded,
  # 0 failed" for all 24 candidates, i.e. it wasn't discriminating anything.
  # The two live suspects are (a) 0 being used as a NoData sentinel in the
  # source tiles that terra isn't recognizing as NA, and/or (b) the
  # mosaic()/terrain()/mask() calls running right at the edge of available
  # disk space (each was estimated to need ~16GB against only 14-15GB free
  # in that run) and silently writing back truncated/placeholder output.
  # This check makes either failure mode a loud stop() instead of a silent
  # "successful" but meaningless screening result.
  check_dem_sane <- function(r, label, min_plausible = -50, max_plausible = 4000, fast = FALSE) {
    if (fast) {
      # Sampled estimate, not an exhaustive scan - used for the pre-aggregation
      # check, where aggregate() is about to read every cell anyway. Doing an
      # exact global() here first would mean reading the whole (potentially
      # large, native-resolution) raster from disk TWICE before it's ever
      # shrunk down. A few thousand regularly-spaced sample points is plenty
      # to catch "this is degenerate" without that extra full pass; the exact
      # check still runs (fast = FALSE, the default) right after aggregation,
      # on the much smaller raster that downstream code actually uses.
      samp <- tryCatch(
        terra::spatSample(r, size = 20000, method = "regular", na.rm = TRUE, values = TRUE),
        error = function(e) NULL
      )
      if (is.null(samp) || nrow(samp) == 0) {
        cat(sprintf("[DEM sanity check: %s] could not sample - skipping fast pre-check (exact check still runs after aggregation).\n", label))
        return(invisible(TRUE))
      }
      v <- samp[[1]]
      r_min <- min(v, na.rm = TRUE); r_max <- max(v, na.rm = TRUE); r_mean <- mean(v, na.rm = TRUE)
      n_cells <- length(v)
      n_na <- sum(is.na(v))
      n_zero <- sum(v == 0, na.rm = TRUE)
      cat(sprintf(
        "[DEM sanity check: %s, SAMPLED n=%d] min=%.1f max=%.1f mean=%.1f | %d NA, %d exactly 0\n",
        label, n_cells, r_min, r_max, r_mean, n_na, n_zero
      ))
    } else {
      stats_r <- as.numeric(global(r, c("min", "max", "mean"), na.rm = TRUE))
      r_min <- stats_r[1]; r_max <- stats_r[2]; r_mean <- stats_r[3]
      n_cells <- as.numeric(ncell(r))
      n_na    <- as.numeric(global(is.na(r), "sum"))[1]
      n_zero  <- tryCatch(as.numeric(global(r == 0, "sum", na.rm = TRUE))[1], error = function(e) NA_real_)
      cat(sprintf(
        "[DEM sanity check: %s] min=%.1f max=%.1f mean=%.1f | %.0f/%.0f cells NA, %s/%.0f cells exactly 0\n",
        label, r_min, r_max, r_mean, n_na, n_cells,
        if (is.na(n_zero)) "NA" else sprintf("%.0f", n_zero), n_cells
      ))
    }
    degenerate <- !is.finite(r_min) || !is.finite(r_max) ||
      (r_max - r_min) < 1 ||                     # essentially flat/constant everywhere
      (r_min < min_plausible || r_max > max_plausible)
    if (degenerate) {
      stop(
        "[DEM sanity check FAILED: ", label, "] The elevation raster is degenerate ",
        "(near-constant, all/mostly 0, or outside a plausible ", min_plausible, " to ",
        max_plausible, " m range for this region) - min=", round(r_min, 1),
        " max=", round(r_max, 1), " mean=", round(r_mean, 1), ". ",
        "This is exactly the failure signature that previously produced a silent ",
        "'0 excluded, 0 failed' elevation/aspect screen (every station trivially ",
        "passed because the DEM read back as ~0 everywhere). Likely causes:\n",
        "  1) The source tile(s) use 0 as a NoData sentinel that terra isn't ",
        "recognizing as NA - check NAflag(tile) on each source tile; set it ",
        "explicitly (e.g. NAflag(tile) <- 0) before mosaicking if confirmed.\n",
        "  2) mosaic()/terrain()/mask() ran short on scratch disk space (the ",
        "prior run warned ~16GB needed vs only 14-15GB free) and wrote back ",
        "truncated output. Free up disk space, point terra::terraOptions(tempdir=)",
        " at a drive with more room, and/or raise DEM_AGGREGATE_FACTOR below to ",
        "shrink the raster before terrain()/mask() ever touch it.\n",
        "Fix the underlying issue and re-run - this stop() is deliberate so daily.csv ",
        "is never built from a screen that isn't actually screening anything.",
        call. = FALSE
      )
    }
    invisible(TRUE)
  }
  
  # Load the Nechako DEM
  # ------------------------------------------------------------------------------
  # UPDATE: this now loads the ALREADY-MOSAICKED DEM directly from a previous
  # run instead of re-discovering BC_Elevation-*.tif tiles and re-mosaicking
  # them (and instead of falling back to the ESRI Grid or an elevatr
  # download). No downloading and no mosaicking happens in this run - the
  # tile-discovery / spatSample-per-tile / mosaic() / ESRI-Grid / elevatr
  # logic that used to live here has been removed entirely; it's the exact
  # same file this script itself writes out further down (DEM_MOSAIC_OUT_PATH,
  # "Spatial/nechako_dem_mosaic.tif"), so re-deriving it here would just be
  # redoing already-finished work.
  # If you DO need to rebuild the mosaic from tiles again (e.g. new/updated
  # BC_Elevation-*.tif tiles have been added to Spatial/), restore the
  # previous tile-discovery/mosaic block, or simply delete
  # DEM_MOSAIC_PATH below and re-run an earlier version of this script.
  DEM_MOSAIC_PATH <- "Spatial/nechako_dem_mosaic.tif"
  if (!file.exists(DEM_MOSAIC_PATH)) {
    stop(
      "[DEM] Expected an already-mosaicked DEM at '", DEM_MOSAIC_PATH, "' ",
      "from a previous run, but it was not found. This script has been set ",
      "to load that file directly rather than re-download/re-mosaic tiles. ",
      "Either place the mosaicked DEM at that path, or restore the original ",
      "tile-discovery/mosaic/elevatr-fallback logic if you need to build it ",
      "from scratch."
    )
  }
  message("[OK] Loading already-mosaicked DEM from '", DEM_MOSAIC_PATH, "' (no download, no mosaicking this run).")
  dem <- rast(DEM_MOSAIC_PATH)
  
  # NOTE: extract(dem, ...) names its output column after the DEM's band name -
  # it is only called "lyr.1" for certain unnamed in-memory rasters. A DEM read
  # from a real GeoTIFF (or produced by elevatr) usually keeps its own name
  # (e.g. the source filename, or "elevation"), so a hardcoded "$lyr.1" below
  # would silently return NULL instead of erroring - which is what caused:
  #   Error in `stopifnot()`: object 'elevation' not found
  # (stations_near$elevation <- NULL is a silent no-op: it never created the
  # column at all). Pin the DEM down to one band with a known name so the
  # extraction is no longer guessing.
  if (terra::nlyr(dem) > 1) {
    message("Note: DEM has ", terra::nlyr(dem), " layers - using layer 1.")
    dem <- dem[[1]]
  }
  names(dem) <- "elevation"
  
  # ---- CROP to the basin+buffer BEFORE any sanity check or aggregation ----
  # At this point `dem` can still be an entire province-wide mosaic (the GEE
  # export this pipeline reads from exports a whole-BC bounding box) - only a
  # small slice of it, around the Nechako basin, is ever used downstream.
  # Running the sanity check / aggregate() on the full province instead of
  # this slice caused three separate symptoms in one run:
  #   (a) check_dem_sane() FALSE-POSITIVE: its plausible range (-50 to 4000 m,
  #       set for the Nechako region specifically) correctly rejects a
  #       legitimate CDEM value elsewhere in BC (documented province-wide
  #       range is -226 to 5944 m - coastal depressions, Coast/Rocky Mountain
  #       peaks far from Nechako) as if it were degenerate/bad data.
  #   (b) the recurring "Estimated disk space needed ... 16GB, Available: 9GB"
  #       warning in this pipeline's log - from mosaicking/aggregating the
  #       ENTIRE province every run instead of just the area actually needed.
  #   (c) the high "exactly 0" fraction seen in the sanity-check output
  #       (e.g. 58.6% of cells) - most of a province-wide export is Pacific
  #       Ocean / non-BC land that clipToCollection() upstream in the GEE
  #       script masks to 0. That's expected and harmless out there; it just
  #       never needed to be loaded here in the first place.
  # Crop to the 100 km buffer's bounding box with a small margin so nothing
  # right at the buffer's edge is lost to extent-rounding.
  nechako_buffer_demcrs_precrop <- st_transform(nechako_buffer, st_crs(dem))
  crop_ext <- ext(vect(nechako_buffer_demcrs_precrop))
  CROP_MARGIN_DEG <- 0.05  # ~5 km margin at this latitude - just edge-rounding insurance
  crop_ext <- crop_ext + CROP_MARGIN_DEG
  dem_full_extent <- ext(dem)  # kept only to report the before/after below
  dem <- crop(dem, crop_ext)
  names(dem) <- "elevation"
  cat(sprintf(
    "[OK] Cropped DEM from the full mosaic extent (xmin=%.2f xmax=%.2f ymin=%.2f ymax=%.2f) down to the basin+buffer area (xmin=%.2f xmax=%.2f ymin=%.2f ymax=%.2f) before any further processing.\n",
    dem_full_extent$xmin, dem_full_extent$xmax, dem_full_extent$ymin, dem_full_extent$ymax,
    ext(dem)$xmin, ext(dem)$xmax, ext(dem)$ymin, ext(dem)$ymax
  ))
  
  # ---- Sanity-check the DEM BEFORE spending time on terrain()/mask()/extraction ----
  # fast = TRUE: sampled, not exhaustive - this is a quick pre-check; the
  # non-sampled check_dem_sane() call below (post-aggregation, or at native
  # resolution if aggregation is disabled) is the exact/exhaustive one.
  # Now runs on the CROPPED (basin+buffer-only) raster, so the -50/4000
  # plausibility window is actually checking the right region.
  check_dem_sane(dem, "post-mosaic/load, pre-processing", fast = TRUE)
  
  # ---- Optional: aggregate to a coarser resolution before terrain()/mask() ----
  # This step originally existed because mosaic()/terrain()/mask() were being
  # run on the FULL PROVINCE-WIDE mosaic (see the crop-before-aggregate note
  # above) - that run warned it needed ~16GB against only 14-15GB free, and
  # aggregating by DEM_AGGREGATE_FACTOR^2 fewer cells was the fix for that
  # memory pressure. Now that the crop to basin+buffer happens BEFORE this
  # point (a 100 km buffer around one basin, not all of BC), the raster
  # terrain()/mask() actually runs on is already small enough that aggregation
  # is no longer solving a real memory problem - it would just be discarding
  # native-resolution precision for no benefit, on a raster this size. Left
  # here as a toggle (set > 1) in case a future run ever needs a much larger
  # buffer/area and hits memory pressure again.
  DEM_AGGREGATE_FACTOR <- 1  # 1 = disabled (native resolution); the basin+buffer crop above already keeps this raster small enough that aggregation isn't needed
  if (DEM_AGGREGATE_FACTOR > 1) {
    message(sprintf(
      "Aggregating DEM by factor %d (%dx fewer cells) before terrain()/mask() to reduce disk/memory pressure...",
      DEM_AGGREGATE_FACTOR, DEM_AGGREGATE_FACTOR^2
    ))
    dem <- aggregate(dem, fact = DEM_AGGREGATE_FACTOR, fun = "mean", na.rm = TRUE)
    names(dem) <- "elevation"
    check_dem_sane(dem, "post-aggregation")
  } else {
    # Aggregation disabled - this IS the raster terrain()/extraction will use,
    # so run the exact (non-sampled) check here instead of relying only on
    # the fast pre-check above.
    check_dem_sane(dem, "post-mosaic/load, native resolution (aggregation disabled)")
  }
  
  # Reproject the basin boundary/buffer once to the DEM's CRS and reuse these
  # below, instead of leaving nechako_bound/nechako_buffer in their original
  # (KML/WGS84) CRS and relying on downstream calls to happen to still work.
  nechako_bound_demcrs <- st_transform(nechako_bound, st_crs(dem))
  nechako_buffer_demcrs <- st_transform(nechako_buffer, st_crs(dem))
  
  # DIAGNOSTIC: confirm the DEM actually covers the 100 km buffer, not just the
  # basin polygon - this is the condition that determines whether buffer-zone
  # stations even have a chance at getting real elevation/aspect values below.
  dem_ext <- ext(dem)
  buffer_ext <- ext(vect(nechako_buffer_demcrs))
  covers_buffer <- dem_ext$xmin <= buffer_ext$xmin && dem_ext$xmax >= buffer_ext$xmax &&
    dem_ext$ymin <= buffer_ext$ymin && dem_ext$ymax >= buffer_ext$ymax
  if (covers_buffer) {
    message("DEM extent fully covers the 100 km buffer - buffer-zone stations can be screened.")
  } else {
    warning(
      "DEM extent does NOT fully cover the 100 km buffer around the basin. ",
      "Stations outside the DEM's coverage will still get NA elevation/aspect ",
      "and be dropped by the filter() below unless listed in ",
      "KNOWN_NECHAKO_STATION_IDS. Add more BC_Elevation-*.tif tiles to ",
      "Spatial/ to close the gap."
    )
  }
  
  # ---- INTERNAL GAP CHECK: does every cell INSIDE the buffer actually have data? ----
  # covers_buffer above only compares outer BOUNDING BOXES - it happily passes
  # even if a tile is missing from the MIDDLE of the mosaic, because the two
  # outermost tiles alone can already satisfy "DEM extent >= buffer extent".
  # That's a real risk here: the number of GEE export shards depends on total
  # raster size and isn't fixed (this run produced 3 tiles where an earlier
  # run of the same export produced 6), so a partial download (some shards
  # missing/not yet synced from Drive) would silently pass covers_buffer while
  # still leaving a hole in the middle of the buffer. A station falling in
  # that hole gets NA elevation - not fatal (the n_na_elev/n_na_aspect
  # diagnostic below will still report it), but this check catches it earlier
  # and tells you it's a genuine internal gap rather than "a station just
  # happens to fall outside the DEM's edge".
  buffer_footprint <- rasterize(vect(nechako_buffer_demcrs), dem, field = 1, background = NA)
  n_footprint_cells <- as.numeric(global(!is.na(buffer_footprint), "sum"))[1]
  dem_na_in_buffer <- is.na(dem) & !is.na(buffer_footprint)
  n_na_in_buffer <- as.numeric(global(dem_na_in_buffer, "sum"))[1]
  pct_na_in_buffer <- if (n_footprint_cells > 0) 100 * n_na_in_buffer / n_footprint_cells else NA_real_
  
  cat(sprintf(
    "[DEM buffer coverage check] %.0f / %.0f cells inside the 100 km buffer have DEM data (%.2f%% missing).\n",
    n_footprint_cells - n_na_in_buffer, n_footprint_cells, pct_na_in_buffer
  ))
  
  BUFFER_GAP_TOLERANCE_PCT <- 0.5  # allow a small sliver (e.g. coastline/edge rounding), not a real hole
  if (is.na(pct_na_in_buffer) || pct_na_in_buffer > BUFFER_GAP_TOLERANCE_PCT) {
    # Report roughly WHERE the gap is (centroid of the missing cells), so it's
    # obvious whether this is an edge sliver or a genuinely missing tile.
    gap_cells <- dem_na_in_buffer
    gap_cells[gap_cells == 0] <- NA
    gap_centroid <- tryCatch({
      xy <- terra::xyFromCell(gap_cells, which(!is.na(terra::values(gap_cells))))
      if (nrow(xy) > 0) colMeans(xy) else c(NA, NA)
    }, error = function(e) c(NA, NA))
    
    stop(
      "[DEM buffer coverage check FAILED] ", round(pct_na_in_buffer, 2), "% of the 100 km buffer ",
      "area has NO DEM data - this looks like a missing/incomplete tile in the MIDDLE of the ",
      "mosaic, not just an edge effect (the covers_buffer bounding-box check above can miss this ",
      "entirely). Approximate center of the missing area (DEM CRS coords): ",
      sprintf("x=%.2f, y=%.2f", gap_centroid[1], gap_centroid[2]), ". ",
      "Check that ALL exported BC_Elevation-*.tif shards from the GEE export finished downloading/",
      "syncing into Spatial/ before running this script - the number of shards is not fixed run-to-",
      "run (it depends on the export's pixel dimensions), so don't assume 'I have N tiles, that's ",
      "normal.' If this really is just a rounding sliver at the buffer's edge, raise ",
      "BUFFER_GAP_TOLERANCE_PCT above instead of adding more tiles.",
      call. = FALSE
    )
  }
  
  # ---- Persist the DEM to a fixed path so H8SWEI_Snow.R can reuse it ----
  # H8's build_regions() already does `dem <- crop(dem, project(basin,
  # crs(dem)))` internally before building elevation x aspect regions - so
  # handing it this FULL 100 km-buffer mosaic (rather than a separate
  # basin-only version) is exactly right: Part A needs the buffer extent (to
  # screen neighbouring-basin stations too), H8 only ever needs the basin
  # subset of the same raster, and build_regions() does that cropping itself.
  # One shared file removes the "two DEMs, two paths, must match exactly"
  # failure mode that previously left H8's DEM_PATH pointing at a file that
  # neither script ever produced.
  # NOTE: even though `dem` was loaded directly from this same path above (no
  # download/re-mosaic this run), this write is NOT a no-op: `dem` has since
  # been cropped to the basin+buffer extent (and would be aggregated too, if
  # DEM_AGGREGATE_FACTOR > 1) - so this persists that cropped/validated
  # version back to disk for H8SWEI_Snow.R to reuse, same as before.
  #
  # UPDATE: writing straight to DEM_MOSAIC_OUT_PATH now fails with
  # "[writeRaster] source and target filename cannot be the same" - `dem` is
  # a SpatRaster still backed by (reading from) that exact file on disk, since
  # it was loaded directly from DEM_MOSAIC_PATH above rather than built fresh
  # in memory from tiles like before. terra refuses to let a raster overwrite
  # its own open source file mid-session. Fix: write the cropped/validated
  # result to a temp file first, then replace the original with it (file.rename
  # is a same-volume move/overwrite, not a re-read/re-write of dem) - so
  # DEM_MOSAIC_OUT_PATH ends up holding the new cropped version either way,
  # same as a direct overwrite would have, just without terra objecting to the
  # self-overwrite.
  DEM_MOSAIC_OUT_PATH <- "Spatial/nechako_dem_mosaic.tif"
  # NOTE: the temp filename must still END in ".tif" - writeRaster() guesses
  # the GDAL driver from the file extension, and a name like
  # "nechako_dem_mosaic.tif.tmp" (extension ".tmp") fails with
  # "[writeRaster] cannot guess file type from filename". Inserting the
  # "tmp" marker BEFORE the ".tif" (not after) keeps the real extension last.
  dem_tmp_out_path <- file.path(
    dirname(DEM_MOSAIC_OUT_PATH),
    sub("\\.tif$", "_tmp.tif", basename(DEM_MOSAIC_OUT_PATH))
  )
  writeRaster(dem, dem_tmp_out_path, overwrite = TRUE)
  # Read the just-written temp file back in and swap `dem` over to point at
  # it BEFORE removing/replacing the original - this drops dem's connection
  # to DEM_MOSAIC_OUT_PATH so the file.rename() below isn't moving a file
  # out from under a raster that's still reading from it (terrain()/mask()/
  # extract() calls further below all use `dem`, so it must stay valid).
  dem <- rast(dem_tmp_out_path)
  names(dem) <- "elevation"
  # Explicitly drop any lingering GDAL handle terra may still hold open on
  # DEM_MOSAIC_OUT_PATH (the object that used to point at it was just
  # reassigned above, but terra/GDAL file handles on Windows aren't always
  # released the instant an R object is overwritten) - without this,
  # file.rename() below can fail to replace a file Windows still considers
  # in use, even though no R variable references it anymore.
  invisible(gc())
  ok_swap <- file.rename(dem_tmp_out_path, DEM_MOSAIC_OUT_PATH)
  if (!ok_swap) {
    # file.rename() can fail across filesystems/drives (rare here, since both
    # paths are under Spatial/, but not impossible on some network drives) -
    # fall back to copy + remove instead of leaving the validated DEM stranded
    # under the temp filename.
    file.copy(dem_tmp_out_path, DEM_MOSAIC_OUT_PATH, overwrite = TRUE)
    file.remove(dem_tmp_out_path)
  }
  # Re-point `dem` at the final path now that the swap is done, so it's not
  # left referencing a filename (the temp one) that no longer exists.
  dem <- rast(DEM_MOSAIC_OUT_PATH)
  names(dem) <- "elevation"
  cat(sprintf(
    "[OK] Wrote validated DEM to '%s' (buffer extent, %s deg/m resolution) - set DEM_PATH to this SAME path in H8SWEI_Snow.R.\n",
    DEM_MOSAIC_OUT_PATH, paste(round(res(dem), 5), collapse = " x ")
  ))
  
  # Calculate slope aspect from the DEM (in degrees)
  aspect_raster <- terrain(dem, v = "aspect", unit = "degrees")
  
  # Reproject candidate stations to match the DEM coordinate reference system
  stations_near <- st_transform(stations_near, crs(dem))
  
  # Extract Elevation and Aspect for each candidate station
  stations_near$elevation <- extract(dem, stations_near)$elevation
  stations_near$aspect <- extract(aspect_raster, stations_near)$aspect
  
  # DIAGNOSTIC: any candidate station falling outside the DEM's extent (e.g. if
  # it's still the basin-only ESRI Grid, or a tile gap in the mosaicked CDEM
  # tiles) will get elevation = NA here. In a dplyr::filter(), "NA >= x"
  # evaluates to NA, and filter() treats NA rows as "drop silently" - so
  # without this check, such stations vanish with no explanation and it looks
  # like the whole basin failed the criteria.
  n_na_elev <- sum(is.na(stations_near$elevation))
  n_na_aspect <- sum(is.na(stations_near$aspect))
  message(sprintf(
    "DEM extraction: %d/%d candidate stations got NA elevation, %d/%d got NA aspect (outside the DEM's coverage).",
    n_na_elev, nrow(stations_near), n_na_aspect, nrow(stations_near)
  ))
  if (n_na_elev > 0 || n_na_aspect > 0) {
    coverage_note <- if (covers_buffer) {
      paste0(
        "The DEM does cover the full buffer, so these particular ",
        "stations are outside the buffer itself (or fall in a small nodata ",
        "gap between tiles) rather than missing tile coverage entirely - ",
        "worth checking them individually."
      )
    } else {
      paste0(
        "The DEM does not fully cover the buffer (see warning above) - ",
        "add more BC_Elevation-*.tif tiles to Spatial/ to close the gap if ",
        "you want these stations screened too."
      )
    }
    message(
      "These NA stations are outside the DEM's coverage and cannot be screened ",
      "by elevation/aspect at all. They will NOT survive the filter() below ",
      "(NA comparisons drop the row) unless they are also listed in ",
      "KNOWN_NECHAKO_STATION_IDS further down. ",
      coverage_note
    )
  }
  if (nrow(stations_near) > 0) {
    print(summary(stations_near$elevation))
  }
  
  # ---- Degenerate-value check on the EXTRACTED station values themselves ----
  # check_dem_sane() above should already stop() before we get here, but this
  # is a second, independent guard on what actually feeds the elevation/aspect
  # filter below - so an issue introduced anywhere between raster-build and
  # extraction (reprojection, station-buffer geometry, etc.) can't slip
  # through either. This is the check that would have directly caught the
  # previous bug: all 24 stations reading exactly elevation = 0.
  valid_elev <- stations_near$elevation[!is.na(stations_near$elevation)]
  if (length(valid_elev) > 0 && (max(valid_elev) - min(valid_elev)) < 1) {
    stop(
      "Extracted station elevations are degenerate (all ~", round(mean(valid_elev), 1),
      " m, range < 1 m across ", length(valid_elev), " station(s)) - this is the same ",
      "failure signature as the previous silent '0 excluded, 0 failed' screening bug. ",
      "The elevation/aspect filter below would trivially pass or fail EVERY station ",
      "instead of actually discriminating between them. Stopping instead of writing ",
      "a daily.csv built on a non-functional screen - inspect the DEM sanity-check ",
      "messages above (and try plot(dem)) before re-running.",
      call. = FALSE
    )
  }
  
  # Determine Nechako Basin's mean elevation
  basin_elev <- mask(dem, vect(nechako_bound_demcrs))
  
  # NOTE: global(..., "mean") returns one row per raster LAYER. If dem/basin_elev
  # ends up with more than one layer (e.g. because nechako_bound has multiple
  # polygon features, or the source DEM/mosaic wasn't fully flattened to a
  # single band), basin_mean_elev would silently become a vector instead of a
  # single number - which is what caused:
  #   Error in `stopifnot()`: ..1 must be of size 24 or 1, not size 100.
  # further down in the filter() step. Force it down to one layer/one number
  # here instead, with a heads-up if that happens.
  if (terra::nlyr(basin_elev) > 1) {
    message(
      "Note: DEM/basin mask has ", terra::nlyr(basin_elev), " layers - ",
      "using layer 1 for the basin mean elevation. If that's not the layer ",
      "you expect, check nechako_bound (nrow: ", nrow(nechako_bound), ") and dem (nlyr: ", terra::nlyr(dem), ")."
    )
    basin_elev <- basin_elev[[1]]
  }
  basin_mean_elev <- as.numeric(global(basin_elev, "mean", na.rm = TRUE)$mean)[1]
  if (!is.finite(basin_mean_elev) || abs(basin_mean_elev) < 1) {
    stop(
      "basin_mean_elev computed as ", basin_mean_elev, " - implausible for the ",
      "Nechako basin (expected several hundred meters). The mask()/mean step ",
      "produced a degenerate result even though the earlier DEM checks passed - ",
      "inspect `basin_elev` (e.g. plot(basin_elev)) before trusting the screen.",
      call. = FALSE
    )
  }
  cat(sprintf("[OK] Basin mean elevation: %.1f m (used for the +/-300m screening window below)\n", basin_mean_elev))
  
  # Define filtering parameters based on your criteria
  elev_window <- 300      # +/- 300m elevation window
  aspect_min <- 90        # Replace with your predominant Nechako slope minimum
  aspect_max <- 270       # Replace with your predominant Nechako slope maximum
  
  # ----------------------------------------------------------------------------
  # FALLBACK: fetch elevation/aspect for DEM-coverage-gap stations from a
  # public source (AWS Terrain Tiles via elevatr), instead of just excluding
  # them or blindly force-including them via KNOWN_NECHAKO_STATION_IDS.
  # ----------------------------------------------------------------------------
  # WHY THIS EXISTS: dem_coverage_excluded_stations.csv from a prior run
  # showed 2 stations (1A05P, 4A10P) dropped here NOT because they're outside
  # the basin - a KMZ-based check confirmed both sit 80-90 km from the basin
  # boundary, comfortably inside the 100 km buffer (17.8 km and 11.1 km
  # margin respectively) - but because the local DEM mosaic's own extent
  # clips a real sliver of the buffer (missing ~2.5-4.3 km past its edge).
  # Re-exporting/re-mosaicking more BC_Elevation-*.tif tiles is the "real"
  # fix, but that's a manual GIS step outside this script's control. This
  # fallback closes the gap automatically, using the exact same DEM source
  # (elevatr AWS Terrain Tiles) already proven working elsewhere in this
  # project (00Elevation_Basin.R's overview-map DEM step).
  #
  # SCOPE: only ever touches stations that are ALREADY NA after the primary
  # local DEM extraction above - it never overrides a real value the local
  # DEM produced, and it never runs at all if every station already has
  # elevation/aspect. Toggle with CSNOW_FALLBACK... no wait, use:
  #   Sys.setenv(DEM_FALLBACK_ELEVATR = "FALSE")  # to disable
  # before sourcing this script, if you'd rather these stations stay
  # excluded/NA and be handled manually via KNOWN_NECHAKO_STATION_IDS instead.
  USE_ELEVATR_FALLBACK <- isTRUE(as.logical(Sys.getenv("DEM_FALLBACK_ELEVATR", "TRUE")))
  
  fetch_elevation_aspect_elevatr <- function(lon, lat, station_id,
                                             patch_deg = 0.02,
                                             elevatr_zoom = 11) {
    # patch_deg: half-width (deg) of the small local raster pulled around each
    # point so aspect can be computed from real neighbouring cells - a single
    # point has no neighbourhood, so elevation alone isn't enough for aspect.
    # ~0.02 deg (~1-2 km at this latitude) is small enough to stay cheap per
    # station while still giving terrain() a real 3x3+ neighbourhood at z=11
    # (elevatr's z=11 AWS tiles are ~30m/px, well under patch_deg's ~2km span).
    pt_sf <- sf::st_as_sf(
      data.frame(station_id = station_id, x = lon, y = lat),
      coords = c("x", "y"), crs = 4326
    )
    patch_bbox_sf <- sf::st_as_sf(
      data.frame(x = c(lon - patch_deg, lon + patch_deg),
                 y = c(lat - patch_deg, lat + patch_deg)),
      coords = c("x", "y"), crs = 4326
    )
    r <- tryCatch(
      elevatr::get_elev_raster(locations = patch_bbox_sf, z = elevatr_zoom,
                               clip = "bbox", src = "aws",
                               neg_to_na = FALSE, verbose = FALSE),
      error = function(e) {
        message("    [FAIL] elevatr download failed for station ", station_id,
                ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(r)) return(list(elevation = NA_real_, aspect = NA_real_, source = "elevatr_failed"))
    
    r <- terra::rast(r)
    names(r) <- "elevation"
    elev_val <- tryCatch(
      terra::extract(r, terra::vect(pt_sf))$elevation,
      error = function(e) NA_real_
    )
    aspect_r <- tryCatch(
      terra::terrain(r, v = "aspect", unit = "degrees"),
      error = function(e) NULL
    )
    aspect_val <- if (!is.null(aspect_r)) {
      tryCatch(terra::extract(aspect_r, terra::vect(pt_sf))$aspect, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    
    if (length(elev_val) == 0) elev_val <- NA_real_
    if (length(aspect_val) == 0) aspect_val <- NA_real_
    list(elevation = elev_val[1], aspect = aspect_val[1], source = "elevatr_aws_z11_fallback")
  }
  
  n_na_before_fallback <- sum(is.na(stations_near$elevation) | is.na(stations_near$aspect))
  fallback_log <- data.frame(
    station_id = character(0), lon = numeric(0), lat = numeric(0),
    elevation = numeric(0), aspect = numeric(0), source = character(0),
    stringsAsFactors = FALSE
  )
  
  if (USE_ELEVATR_FALLBACK && n_na_before_fallback > 0) {
    if (!requireNamespace("elevatr", quietly = TRUE)) {
      message(
        "  [SKIP] DEM_FALLBACK_ELEVATR is TRUE but the 'elevatr' package isn't ",
        "installed - install.packages('elevatr') to enable the public-source ",
        "elevation/aspect fallback for the ", n_na_before_fallback, " station(s) ",
        "with NA elevation/aspect. Continuing without it; those stations will ",
        "still be excluded below unless listed in KNOWN_NECHAKO_STATION_IDS."
      )
    } else {
      message(sprintf(
        "  [OK] Attempting public-source (elevatr/AWS Terrain Tiles) elevation/aspect fallback for %d station(s) outside the local DEM's coverage...",
        n_na_before_fallback
      ))
      na_idx <- which(is.na(stations_near$elevation) | is.na(stations_near$aspect))
      # Fallback runs in WGS84 lon/lat regardless of the local DEM's CRS, since
      # elevatr/AWS tiles are always fetched in geographic coordinates.
      stations_near_wgs84 <- suppressWarnings(sf::st_transform(stations_near, 4326))
      na_coords <- sf::st_coordinates(stations_near_wgs84[na_idx, ])
      
      for (k in seq_along(na_idx)) {
        i <- na_idx[k]
        sid <- stations_near$station_id[i]
        lon_i <- na_coords[k, "X"]; lat_i <- na_coords[k, "Y"]
        res <- fetch_elevation_aspect_elevatr(lon_i, lat_i, sid)
        fallback_log <- rbind(fallback_log, data.frame(
          station_id = sid, lon = lon_i, lat = lat_i,
          elevation = res$elevation, aspect = res$aspect, source = res$source,
          stringsAsFactors = FALSE
        ))
        if (!is.na(res$elevation)) stations_near$elevation[i] <- res$elevation
        if (!is.na(res$aspect))    stations_near$aspect[i]    <- res$aspect
        if (!is.na(res$elevation) && !is.na(res$aspect)) {
          message(sprintf("    [OK] %s: elevation=%.1f m, aspect=%.1f deg (elevatr AWS fallback).", sid, res$elevation, res$aspect))
        } else {
          message(sprintf("    [SKIP] %s: elevatr fallback did not produce a usable value either (elevation=%s, aspect=%s) - still excluded below.",
                          sid, ifelse(is.na(res$elevation), "NA", "ok"), ifelse(is.na(res$aspect), "NA", "ok")))
        }
      }
      n_filled <- sum(!is.na(fallback_log$elevation) & !is.na(fallback_log$aspect))
      message(sprintf(
        "  [OK] Public-source fallback filled %d/%d previously-NA station(s). These are flagged elevation_source = 'elevatr_aws_z11_fallback' below so they remain distinguishable from local-DEM values.",
        n_filled, n_na_before_fallback
      ))
      if (nrow(fallback_log) > 0) {
        write_csv(fallback_log, file.path(SNOW_DATA_DIR, "dem_coverage_elevatr_fallback_log.csv"))
        message(sprintf("  [OK] Fallback attempt log (all attempts, including failures) written to '%s'.",
                        file.path(SNOW_DATA_DIR, "dem_coverage_elevatr_fallback_log.csv")))
      }
    }
  } else if (!USE_ELEVATR_FALLBACK && n_na_before_fallback > 0) {
    message(sprintf(
      "  [INFO] %d station(s) have NA elevation/aspect from the local DEM, but DEM_FALLBACK_ELEVATR=FALSE - skipping the public-source fallback (set it to TRUE, or unset it, to enable).",
      n_na_before_fallback
    ))
  }
  
  # Track which stations were filled by the fallback vs. the local DEM, so the
  # exclusion report/CSV below can tell the two apart instead of implying every
  # remaining exclusion is still a raw DEM gap.
  stations_near$elevation_source <- ifelse(
    stations_near$station_id %in% fallback_log$station_id[!is.na(fallback_log$elevation) & !is.na(fallback_log$aspect)],
    "elevatr_aws_z11_fallback", "local_dem_mosaic"
  )
  
  # Filter the stations
  # NOTE: explicit !is.na() guards below make the NA-drop behavior visible/
  # intentional instead of relying on filter()'s implicit "NA row => dropped"
  # rule, which is what silently zeroed out every buffer-zone station before
  # (they have elevation/aspect = NA because they're outside the DEM's extent,
  # not because they failed the numeric test).
  # NOTE: this now runs AFTER the elevatr fallback above, so any station the
  # fallback successfully filled is no longer NA here and will be evaluated
  # normally against elev_window/aspect_min/aspect_max - it's just no longer
  # guaranteed-excluded the way it was before the fallback existed.
  stations_outside_dem <- stations_near %>%
    filter(is.na(elevation) | is.na(aspect))
  
  # FIX (point 4): previously only a COUNT of DEM-coverage-excluded stations
  # was ever printed ("2 excluded for being outside the DEM's extent") - the
  # station IDs themselves were computed into `stations_outside_dem` but
  # never surfaced anywhere, so there was no way to know WHICH stations were
  # lost, or to tell "genuinely outside the basin" apart from "just missing
  # a DEM tile" without re-running with extra debugging. Print the IDs and
  # persist them to a CSV so they're inspectable and actionable - e.g. to
  # decide whether any belong in KNOWN_NECHAKO_STATION_IDS below. This does
  # NOT change which stations survive the filter (that behaviour is
  # unchanged) - it only makes the existing exclusion visible.
  # UPDATE: reflects post-fallback state - a station only appears here now if
  # BOTH the local DEM AND the elevatr fallback (when enabled) failed to
  # produce a value, or the fallback was disabled/unavailable.
  if (nrow(stations_outside_dem) > 0) {
    excluded_ids <- paste(stations_outside_dem$station_id, collapse = ", ")
    message(sprintf(
      "  [DEM coverage gap] station(s) STILL dropped from screening after the elevatr fallback (NA elevation/aspect - not evaluated, not necessarily outside the basin): %s",
      excluded_ids
    ))
    excluded_coords <- suppressWarnings(st_coordinates(st_transform(stations_outside_dem, 4326)))
    excluded_review <- data.frame(
      station_id = stations_outside_dem$station_id,
      lon = excluded_coords[, "X"],
      lat = excluded_coords[, "Y"],
      reason = "NA elevation/aspect - outside local DEM raster coverage AND elevatr public-source fallback failed/disabled this run",
      dem_covers_full_buffer = covers_buffer,
      elevatr_fallback_attempted = USE_ELEVATR_FALLBACK,
      stringsAsFactors = FALSE
    )
    write_csv(excluded_review, file.path(SNOW_DATA_DIR, "dem_coverage_excluded_stations.csv"))
    message(sprintf(
      "  [DEM coverage gap] Written to '%s' for review - add any of these to KNOWN_NECHAKO_STATION_IDS below ONLY if you've confirmed they belong in the basin; this list alone is not evidence of basin membership, just of a raster gap.",
      file.path(SNOW_DATA_DIR, "dem_coverage_excluded_stations.csv")
    ))
  }
  
  surviving_stations <- stations_near %>%
    filter(
      !is.na(elevation), !is.na(aspect),
      elevation >= (basin_mean_elev - elev_window),
      elevation <= (basin_mean_elev + elev_window),
      aspect >= aspect_min,
      aspect <= aspect_max
    )
  
  if (any(surviving_stations$elevation_source == "elevatr_aws_z11_fallback")) {
    fallback_survivors <- paste(
      surviving_stations$station_id[surviving_stations$elevation_source == "elevatr_aws_z11_fallback"],
      collapse = ", "
    )
    message(sprintf(
      "  [OK] %s passed elevation/aspect screening using the elevatr public-source fallback (local DEM had no coverage there). elevation_source column in daily.csv/surviving_stations marks these - worth a sanity spot-check against real terrain before fully trusting them alongside local-DEM-screened stations.",
      fallback_survivors
    ))
  }
  
  cat(sprintf(
    "Screening result: %d station(s) passed elevation/aspect criteria, %d excluded for being outside the DEM's extent (no elevation data), %d evaluated but failed the numeric criteria.\n",
    nrow(surviving_stations),
    nrow(stations_outside_dem),
    nrow(stations_near) - nrow(surviving_stations) - nrow(stations_outside_dem)
  ))
  print(paste("Found", nrow(surviving_stations), "stations matching spatial & topo criteria."))
  
  # ------------------------------------------------------------------------------
  # 1b. Force-include the known in-basin Nechako stations
  # ------------------------------------------------------------------------------
  # The elevation/aspect screen above is meant to vet ADDITIONAL neighbouring
  # candidates, not to re-vet the stations that are already known to be
  # physically inside the Nechako basin (they were the original
  # station_swe_obs.csv stations). If any of them happen to sit outside the
  # elev_window/aspect window used above (e.g. a valley-bottom pillow just
  # outside a ridge-driven aspect range), the filter would silently drop them
  # and you'd lose ground-truth stations instead of only gaining candidates.
  # Filled in from the reference daily.csv (5 distinct stations: 1B01P, 1B02P,
  # 1B08P, 4B15P, 4B16P) - same ID format as in the raw archive CSV below,
  # leading underscore included:
  KNOWN_NECHAKO_STATION_IDS <- c("_1B01P", "_1B02P", "_1B08P", "_4B15P", "_4B16P")
  if (length(KNOWN_NECHAKO_STATION_IDS) == 0) {
    warning(
      "KNOWN_NECHAKO_STATION_IDS is empty - if your known in-basin ",
      "stations don't happen to pass the elevation/aspect filter above, ",
      "they will be silently dropped from daily.csv. Fill this vector in ",
      "with their Pillow_IDs before relying on the output."
    )
  }
  
  # ------------------------------------------------------------------------------
  # 2. Format the Time-Series Data
  # ------------------------------------------------------------------------------
  # Rather than an unbuilt per-station API call, this joins the province-wide
  # BC ASWS archive (raw_archive, already loaded above - see "Identify which
  # WFS attribute actually holds the station ID") to the elevation/aspect-
  # screened station list.
  
  # Union: elevation/aspect-screened candidates + the known in-basin anchors,
  # de-duplicated by normalized station ID
  candidate_ids <- union(
    surviving_stations$station_id_norm,
    normalize_stn_id(KNOWN_NECHAKO_STATION_IDS)
  )
  stations_to_use <- stations_near %>% filter(station_id_norm %in% candidate_ids)
  
  if (nrow(stations_to_use) == 0) {
    stop(
      "No candidate stations (screened or known-in-basin) to process - check ",
      "elev_window/aspect_min/aspect_max above and KNOWN_NECHAKO_STATION_IDS."
    )
  }
  
  compiled_swe <- list()
  skipped <- character(0)
  
  for (i in seq_len(nrow(stations_to_use))) {
    
    stn_id <- stations_to_use$station_id[i]
    stn_id_norm <- stations_to_use$station_id_norm[i]
    
    station_swe <- raw_archive %>% filter(pillow_id_norm == stn_id_norm)
    
    if (nrow(station_swe) == 0) {
      skipped <- c(skipped, stn_id)
      next
    }
    
    # WGS84 lat/lon for the H8 mapping
    coords <- st_coordinates(st_transform(stations_to_use[i, ], 4326))
    
    formatted_data <- station_swe %>%
      transmute(
        station_id = stn_id,
        lon = coords[1, "X"],
        lat = coords[1, "Y"],
        date = date,
        swe_mm = swe_mm
      )
    
    compiled_swe[[stn_id]] <- formatted_data
  }
  
  if (length(skipped) > 0) {
    message(
      "No SWE records found in the raw archive for ", length(skipped),
      " station(s), skipped: ", paste(skipped, collapse = ", ")
    )
  }
  
  # Combine all formatted station data into a single dataframe
  final_daily_csv <- bind_rows(compiled_swe)
  
  n_distinct_stations <- length(unique(final_daily_csv$station_id))
  
  if (nrow(final_daily_csv) == 0) {
    known_ids_note <- if (length(KNOWN_NECHAKO_STATION_IDS) == 0) {
      paste0(
        "MOST LIKELY CAUSE: KNOWN_NECHAKO_STATION_IDS (further up) is empty. ",
        "Every candidate that reached this join came ONLY from the elevation/",
        "aspect screen (", length(unique(stations_to_use$station_id_norm)),
        " distinct station(s): ", paste(unique(stations_to_use$station_id), collapse = ", "),
        ") - none of your 3 already-known, already-verified in-basin stations ",
        "were force-included, because that vector was left empty. If those 3 ",
        "stations are the ones with real archived SWE data, this is exactly ",
        "what an empty daily.csv looks like. Fill in KNOWN_NECHAKO_STATION_IDS ",
        "with their Pillow_IDs and re-run before chasing anything else.\n\n"
      )
    } else {
      ""
    }
    stop(
      known_ids_note,
      "final_daily_csv has 0 rows after joining - daily.csv would be empty. ",
      "Other things worth checking: that Pillow_ID values in ", RAW_ASWS_ARCHIVE_PATH,
      " actually match the WFS '", station_id_col, "' column (see the ",
      "'WFS station-ID column detection' message above), and that the ",
      "'variable == \"SWE\"' rows exist for these specific stations (it's ",
      "possible for a station to exist in the WFS/archive without ever ",
      "reporting SWE specifically)."
    )
  }
  
  if (n_distinct_stations < 4) {
    warning(
      "Only ", n_distinct_stations, " distinct station(s) made it into daily.csv. ",
      "H8SWEI_Snow.R's LOSO_CV_MIN_STATIONS gate needs at least 4 distinct ",
      "stations before it will auto-select Tier 2/3 over Tier 1 - Tier 1 will ",
      "still be applied safely below this, but the tiered comparison won't run."
    )
  }
  
  # Save to the root directory for the pipeline to read
  write_csv(final_daily_csv, DAILY_CSV_PATH)
  print(sprintf(
    "Saved daily.csv: %d rows, %d distinct stations, columns station_id/lon/lat/date/swe_mm.",
    nrow(final_daily_csv), n_distinct_stations
  ))
  cat(sprintf("[OK] Part A complete. daily.csv is ready in %s for H8SWEI_Snow.R's bias correction.\n",
              normalizePath(".")))
  
} # end if (RUN_PART_A)

# ==============================================================================
# PART B - REVISED DATA ACQUISITION FOR THE H8 PIPELINE
# ==============================================================================
if (RUN_SATELLITE_AND_CROCUS_PARTS) {
  
  # ----------------------------------------------------------------------------
  # Shared helpers
  # ----------------------------------------------------------------------------
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  reauthenticate <- function() invisible(TRUE)
  
  ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
  }
  
  # Moved below the ensure_dir() definition above - it used to sit right at
  # the top of this block, before ensure_dir() was defined, which threw
  # "could not find function 'ensure_dir'" as soon as Part B started.
  ensure_dir("snow_data")
  
  # ----------------------------------------------------------------------------
  # C-SNOW (Lievens et al. Sentinel-1 snow depth) - local copy, and whether to
  # still also pull raw Sentinel-1 GRD via ASF/CDSE
  # ----------------------------------------------------------------------------
  # C-SNOW (https://ees.kuleuven.be/project/c-snow) is Lievens et al.'s own
  # PRECOMPUTED Sentinel-1 snow-depth product (empirical VH/VV change-detection
  # retrieval, Lievens et al. 2019/2022) - not raw Sentinel-1 backscatter. A
  # local copy already exists under CSNOW_LOCAL_DIR below (one NetCDF per day,
  # "snd_YYYYMMDD.nc", West_US_Canada tile), so there is no need to also
  # download and locally reprocess raw Sentinel-1 GRD scenes through
  # acquire_sentinel1_asf()/acquire_sentinel1_cdse() further down - that raw-GRD
  # path has no terrain-correction step of its own in this script, and would
  # add a MUCH larger local footprint (potentially many TB over the full
  # 2014-present window) than the ~1-2MB/day this precomputed product already
  # ships at. RUN_SENTINEL1_RAW_ASF therefore defaults to FALSE.
  # Set it back to TRUE (env var or here) only if C-SNOW turns out not to cover
  # the Nechako basin adequately - e.g. if the coverage/validation diagnostics
  # printed by acquire_csnow_local() below (fraction of valid basin pixels per
  # date, and its agreement with daily.csv ground-truth SWE/depth once you've
  # checked that yourself) suggest it's not usable here and you need to build
  # your own retrieval from raw SAR instead.
  CSNOW_LOCAL_DIR <- Sys.getenv(
    "CSNOW_LOCAL_DIR",
    unset = "D:/Nechako_Drought/Nechako/snow_data/csnow/West_US_Canada"
  )
  RUN_SENTINEL1_RAW_ASF <- as.logical(Sys.getenv("RUN_SENTINEL1_RAW_ASF", "FALSE"))
  if (is.na(RUN_SENTINEL1_RAW_ASF)) RUN_SENTINEL1_RAW_ASF <- FALSE
  cat(sprintf(
    "[OK] CSNOW_LOCAL_DIR = '%s', RUN_SENTINEL1_RAW_ASF = %s (raw Sentinel-1 GRD download via ASF/CDSE)\n",
    CSNOW_LOCAL_DIR, RUN_SENTINEL1_RAW_ASF
  ))
  
  normalize_id <- function(x) {
    toupper(trimws(gsub("^_", "", as.character(x))))
  }
  
  first_match <- function(x, patterns) {
    idx <- which(vapply(patterns, function(p) any(grepl(p, x, ignore.case = TRUE)), logical(1)))
    if (length(idx) == 0) return(NA_character_)
    hit <- patterns[idx[1]]
    x[which(grepl(hit, x, ignore.case = TRUE))[1]]
  }
  
  maybe_parse_date <- function(x) {
    if (inherits(x, "Date")) return(x)
    if (inherits(x, "POSIXt")) return(as.Date(x))
    x <- as.character(x)
    out <- suppressWarnings(lubridate::parse_date_time(
      x,
      orders = c("Y-m-d", "Y/m/d", "m/d/Y", "d/m/Y", "Ymd", "Y-m-d H:M:S", "Y/m/d H:M:S", "m/d/Y H:M:S"),
      tz = "UTC"
    ))
    as.Date(out)
  }
  
  write_earthdata_netrc <- function(user, pass) {
    if (!nzchar(user) || !nzchar(pass)) return(invisible(FALSE))
    netrc_path <- path.expand("~/.netrc")
    lines <- c(
      sprintf("machine urs.earthdata.nasa.gov login %s password %s", user, pass),
      sprintf("machine cmr.earthdata.nasa.gov login %s password %s", user, pass),
      sprintf("machine n5eil01u.ecs.nsidc.org login %s password %s", user, pass),
      # BUGFIX: daacdata.apps.nsidc.org was missing from this list entirely.
      # acquire_nsidc_https_dir() (used for CMC/NSIDC-0447 and CanSISE/
      # NSIDC-0668) authenticates against daacdata.apps.nsidc.org specifically -
      # its own comments even document that this host uses a URS OAuth2
      # redirect - but with no netrc entry for that exact host, curl's
      # netrc=1L had no matching credentials to offer at any point in the
      # redirect chain, so the request bottomed out on URS's plain login/
      # authorize HTML page (a normal HTTP 200 - not an error status - which
      # is why it wasn't caught as a failed request). This is the actual root
      # cause of "Final URL: https://urs.earthdata.nasa.gov/oauth/authorize?..."
      # showing up instead of the real directory listing for CMC above.
      sprintf("machine daacdata.apps.nsidc.org login %s password %s", user, pass),
      sprintf("machine data.asf.alaska.edu login %s password %s", user, pass),
      sprintf("machine api.daac.asf.alaska.edu login %s password %s", user, pass),
      sprintf("machine search.asf.alaska.edu login %s password %s", user, pass)
    )
    writeLines(lines, netrc_path)
    try(Sys.chmod(netrc_path, "600"), silent = TRUE)
    invisible(TRUE)
  }
  
  # FIX (point 3): lightweight, ONCE-PER-RUN check of whether the Earthdata
  # username/password themselves are valid, independent of any specific
  # DAAC/app. Used to tell apart the two failure modes that both look
  # identical from acquire_nsidc_https_dir() alone (a 200 response that's
  # actually the URS login page): (a) EARTHDATA_USER/PASS are simply wrong,
  # vs (b) the login is fine but this particular app/host
  # (daacdata.apps.nsidc.org) hasn't been authorized for this account yet -
  # a one-time manual step in a browser that no amount of retrying or
  # cookie-jar handling from R can complete on its own. Hits URS's own
  # profile endpoint, which requires nothing beyond a valid login (no
  # per-dataset authorization), so a login-page redirect THERE specifically
  # means the credentials are the problem.
  .earthdata_login_check_cache <- new.env(parent = emptyenv())
  check_earthdata_login <- function(user = Sys.getenv("EARTHDATA_USER", unset = ""),
                                    pass = Sys.getenv("EARTHDATA_PASS", unset = "")) {
    cache_key <- paste0(user, ":", nchar(pass))  # never cache/print the password itself
    if (!is.null(.earthdata_login_check_cache[[cache_key]])) {
      return(.earthdata_login_check_cache[[cache_key]])
    }
    result <- if (!nzchar(user) || !nzchar(pass)) {
      "no_credentials"
    } else {
      resp <- tryCatch(
        httr::GET("https://urs.earthdata.nasa.gov/profile",
                  httr::config(netrc = 1L, followlocation = TRUE), timeout(30)),
        error = function(e) NULL
      )
      if (is.null(resp) || httr::http_error(resp)) {
        "network_error"
      } else {
        body_check <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
        is_login_page <- grepl("<form[^>]*login|Earthdata Login", body_check, ignore.case = TRUE)
        if (is_login_page) "bad_credentials" else "ok"
      }
    }
    .earthdata_login_check_cache[[cache_key]] <- result
    result
  }
  
  download_file_auth <- function(url, destfile,
                                 user = Sys.getenv("EARTHDATA_USER", unset = ""),
                                 pass = Sys.getenv("EARTHDATA_PASS", unset = "")) {
    ensure_dir(dirname(destfile))
    handle <- curl::new_handle(followlocation = TRUE)
    if (nzchar(user) && nzchar(pass) || file.exists(path.expand("~/.netrc"))) {
      curl::handle_setopt(handle, netrc = 1L)
    }
    curl::curl_download(url, destfile = destfile, handle = handle, quiet = TRUE)
    invisible(destfile)
  }
  
  zenodo_record_files <- function(record_id) {
    rec_url <- sprintf("https://zenodo.org/api/records/%s", record_id)
    rec <- tryCatch(
      jsonlite::fromJSON(rec_url, flatten = TRUE),
      error = function(e) stop("Failed to query Zenodo record ", record_id, ": ", conditionMessage(e), call. = FALSE)
    )
    
    files_tbl <- rec$files
    if (is.null(files_tbl) || (is.data.frame(files_tbl) && nrow(files_tbl) == 0)) {
      return(data.frame(filename = character(0), url = character(0), stringsAsFactors = FALSE))
    }
    
    if (!is.data.frame(files_tbl)) files_tbl <- as.data.frame(files_tbl, stringsAsFactors = FALSE)
    
    file_name_col <- intersect(c("filename", "key", "name"), names(files_tbl))
    if (length(file_name_col) == 0) file_name_col <- names(files_tbl)[1]
    
    url_col <- intersect(c("links.download", "download", "url", "links.self"), names(files_tbl))
    if (length(url_col) == 0) {
      url_col <- grep("download", names(files_tbl), value = TRUE)
    }
    if (length(url_col) == 0) {
      # Fallback to the standard Zenodo content endpoint if the API response
      # does not expose a direct download link in the flattened record.
      urls <- sprintf(
        "https://zenodo.org/api/records/%s/files/%s/content",
        record_id, utils::URLencode(files_tbl[[file_name_col[1]]], reserved = TRUE)
      )
    } else {
      urls <- files_tbl[[url_col[1]]]
    }
    
    data.frame(
      filename = as.character(files_tbl[[file_name_col[1]]]),
      url = as.character(urls),
      stringsAsFactors = FALSE
    )
  }
  
  download_zenodo_record <- function(record_id, out_dir, pattern = NULL, overwrite = FALSE) {
    ensure_dir(out_dir)
    files <- zenodo_record_files(record_id)
    if (nrow(files) == 0) {
      warning("Zenodo record ", record_id, " returned no files.")
      return(invisible(character(0)))
    }
    if (!is.null(pattern)) {
      keep <- grepl(pattern, files$filename, ignore.case = TRUE)
      files <- files[keep, , drop = FALSE]
    }
    if (nrow(files) == 0) {
      warning("No Zenodo files matched pattern '", pattern, "' for record ", record_id, ".")
      return(invisible(character(0)))
    }
    
    out_files <- character(0)
    for (i in seq_len(nrow(files))) {
      fn <- files$filename[i]
      url <- files$url[i]
      dest <- file.path(out_dir, fn)
      if (file.exists(dest) && !overwrite) {
        message("  [OK] ", fn, " already exists - skipping.")
        out_files <- c(out_files, dest)
        next
      }
      message("  [DL] ", fn)
      download_file_auth(url, dest)
      out_files <- c(out_files, dest)
    }
    invisible(out_files)
  }
  
  # ----------------------------------------------------------------------------
  # Direct-HTTPS NSIDC DAAC acquisition
  # ----------------------------------------------------------------------------
  # For legacy NSIDC/ECS collections that exist as a CMR collection stub (so
  # verify_cmr_collection() finds them and confirms the version) but were
  # never granule-indexed into CMR's search API - cmr_download_direct()'s
  # granule query then correctly returns 0 granules every time, which looks
  # identical to "this bbox/time window has no data" but actually means
  # "this collection was never searchable that way." Confirmed for both
  # NSIDC-0447 (CMC) and NSIDC-0668 (CanSISE): their own NSIDC landing pages
  # each document a "Programmatic Access Guide for Data on daacdata.apps" -
  # i.e. direct HTTPS file-system access - rather than a CMR granule search
  # API, which is exactly the same category of problem already documented
  # and worked around below for NSIDC-0595/GlobSnow.
  #
  # This lists the actual directory (an Apache/nginx-style autoindex HTML
  # page) and downloads whatever matches, rather than hardcoding a per-year
  # filename/version pattern - NSIDC reprocesses individual years under
  # different version suffixes (e.g. "v01.2" vs "v01.3"), so a guessed
  # pattern would silently miss some years without ever raising an error.
  acquire_nsidc_https_dir <- function(dir_url, out_dir, file_pattern = "\\.(tif|tiff|zip|nc|txt)$",
                                      start_year = NULL, end_year = NULL,
                                      user = Sys.getenv("EARTHDATA_USER", unset = ""),
                                      pass = Sys.getenv("EARTHDATA_PASS", unset = ""),
                                      query_timeout_sec = 60, query_retries = 4) {
    ensure_dir(out_dir)
    
    # FAIL LOUD AND EARLY if credentials are missing, rather than letting the
    # request go out with no .netrc entries at all and land on the URS login
    # page - which was previously indistinguishable, at a glance, from a
    # working-but-empty listing (both looked like "silence, then nothing
    # downloaded"). This is exactly the failure mode that hit CMC/NSIDC-0447
    # above: no clear signal at the point of failure about WHY it fell
    # through to urs.earthdata.nasa.gov/oauth/authorize.
    if (!nzchar(user) || !nzchar(pass)) {
      message(
        "  [SKIP] EARTHDATA_USER/EARTHDATA_PASS are not set (or empty) - ",
        "daacdata.apps.nsidc.org requires an authenticated Earthdata Login ",
        "session, so this request cannot succeed without them. Add both to ",
        "your .env file (create a free account at ",
        "https://urs.earthdata.nasa.gov/users/new if you don't have one) and re-run."
      )
      warning(
        "acquire_nsidc_https_dir(): EARTHDATA_USER/EARTHDATA_PASS not set - ",
        "cannot authenticate to daacdata.apps.nsidc.org.", call. = FALSE
      )
      return(invisible(character(0)))
    }
    write_earthdata_netrc(user, pass)
    
    # NSIDC's own documented access pattern for daacdata.apps.nsidc.org uses
    # wget with an explicit cookie jar (--load-cookies/--save-cookies
    # ~/.urs_cookies --keep-session-cookies --auth-no-challenge=on), not bare
    # .netrc alone - see "Programmatic Access Guide for Data on daacdata.apps"
    # (https://nsidc.org/data/user-resources/help-center/programmatic-access-guide-data-daacdataapps).
    # daacdata.apps.nsidc.org authenticates via a URS (Earthdata Login) OAuth2
    # redirect, and that redirect/cookie handshake is exactly what netrc=1L
    # alone does not reliably complete - the symptom is a 200 response that is
    # actually the URS login/authorize HTML rather than the real directory
    # listing (see the looks_like_login check below). A persistent cookie jar
    # (httr::config(cookies = ...) via a shared handle) mirrors wget's
    # --keep-session-cookies behaviour across the redirect chain.
    cookie_jar <- file.path(tempdir(), "nsidc_urs_cookies.txt")
    # A stale cookie jar left over from a PRIOR failed/half-completed auth
    # attempt (e.g. an earlier run with the wrong password, before it was
    # fixed in .env) can itself keep redirecting back to the login page even
    # after the .netrc fix above - curl will happily reuse a cookie file that
    # never actually completed the OAuth handshake. Start clean on every call.
    if (file.exists(cookie_jar)) unlink(cookie_jar)
    
    # FIX (point 3): a login-page response is a credential/authorization
    # problem, not a transient network blip - it will not resolve itself
    # between attempts a few seconds apart. Give it a short, separate retry
    # budget (in case it genuinely is just a slow cookie handshake, as the
    # original comment allowed for) instead of burning the full
    # query_retries budget - and its escalating Sys.sleep() backoff, which is
    # there for real network flakiness - against a wall that a 5th or 10th
    # identical attempt would not get past either.
    login_page_retry_cap <- min(query_retries, 2)
    saw_login_page <- FALSE
    resp <- NULL
    for (attempt in seq_len(query_retries)) {
      resp <- tryCatch(
        httr::GET(
          dir_url,
          httr::config(netrc = 1L, followlocation = TRUE, cookiejar = cookie_jar, cookiefile = cookie_jar),
          timeout(query_timeout_sec)
        ),
        error = function(e) {
          message(sprintf("  NSIDC directory listing attempt %d/%d failed: %s",
                          attempt, query_retries, conditionMessage(e)))
          NULL
        }
      )
      
      # BUGFIX: a redirect to the URS login/authorize page is a normal HTTP
      # 200, so the old `!httr::http_error(resp)` check alone treated it as a
      # SUCCESSFUL response and `break`-ed out of the retry loop immediately -
      # exactly what produced "Final URL: https://urs.earthdata.nasa.gov/
      # oauth/authorize?..." above with no further diagnosis. Checking for the
      # login-page signature HERE, inside the loop, means a login-page
      # response is now treated as a failed attempt, but capped at
      # login_page_retry_cap attempts (not the full query_retries) since
      # retrying further would not change the outcome - see above.
      if (!is.null(resp) && !httr::http_error(resp)) {
        body_check <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
        is_login_page <- grepl("urs\\.earthdata\\.nasa\\.gov|Earthdata Login|<form[^>]*login", body_check, ignore.case = TRUE)
        if (!is_login_page) break
        saw_login_page <- TRUE
        message(sprintf(
          "  NSIDC directory listing attempt %d/%d landed on the Earthdata Login page (HTTP 200, but not real content) - %s",
          attempt, login_page_retry_cap,
          if (attempt < login_page_retry_cap) "retrying once (session/cookie may just need another handshake)..." else "giving up - this is a credential/authorization wall, not a network blip, so further retries would not help. Running a one-time login check for a definitive diagnosis..."
        ))
        if (attempt >= login_page_retry_cap) { resp <- NULL; break }
      } else if (!is.null(resp)) {
        message(sprintf(
          "  NSIDC directory listing attempt %d/%d returned HTTP %d - %s",
          attempt, query_retries, httr::status_code(resp),
          if (attempt < query_retries) "retrying..." else "giving up."
        ))
      }
      resp <- NULL
      if (attempt < query_retries) Sys.sleep(10 * attempt)
    }
    if (is.null(resp)) {
      if (saw_login_page) {
        # Definitive diagnosis instead of a hedge between "wrong password" and
        # "needs one-time app authorization" - check_earthdata_login() hits a
        # dataset-agnostic URS endpoint, so its result actually distinguishes
        # the two.
        login_status <- check_earthdata_login(user, pass)
        diagnosis <- switch(login_status,
                            bad_credentials = paste0(
                              "DIAGNOSIS: your EARTHDATA_USER/EARTHDATA_PASS themselves are being rejected by ",
                              "Earthdata Login directly (confirmed via a login-only check independent of this ",
                              "dataset) - fix the credentials in .env."
                            ),
                            ok = paste0(
                              "DIAGNOSIS: EARTHDATA_USER/EARTHDATA_PASS are valid and work against Earthdata Login ",
                              "itself - so this is almost certainly the one-time per-application authorization gate, ",
                              "not a credential problem. Visit ", dir_url, " in a browser while logged in to ",
                              "urs.earthdata.nasa.gov and approve/authorize the application if prompted, then re-run."
                            ),
                            network_error = "DIAGNOSIS: could not even reach urs.earthdata.nasa.gov to check credentials - this may be a general network/connectivity issue rather than anything dataset-specific.",
                            "DIAGNOSIS: could not run the login check (no credentials supplied)."
        )
        message("  [SKIP] Could not list NSIDC directory ", dir_url, " - ", diagnosis)
        warning("Could not list NSIDC directory ", dir_url, " - ", diagnosis, call. = FALSE)
      } else {
        message(
          "  [SKIP] Could not list NSIDC directory ", dir_url, " after ", query_retries,
          " attempts - every attempt either errored or returned a non-200 status ",
          "(not a login-page redirect, so this looks like a network/server issue rather ",
          "than a credential problem)."
        )
        warning("Could not list NSIDC directory ", dir_url, " after ", query_retries, " attempts.", call. = FALSE)
      }
      return(invisible(character(0)))
    }
    
    message("Final URL: ", resp$url)
    html <- httr::content(resp, as = "text", encoding = "UTF-8")
    looks_like_listing <- grepl("Index of|Parent Directory|<title>Index", html, ignore.case = TRUE)
    looks_like_login <- grepl("urs\\.earthdata\\.nasa\\.gov|Earthdata Login|<form[^>]*login", html, ignore.case = TRUE)
    if (!looks_like_listing || looks_like_login) {
      warning(
        "NSIDC directory listing at ", dir_url, " returned HTTP ", httr::status_code(resp),
        " but the response does not look like a real directory index (",
        if (looks_like_login) "it looks like an Earthdata Login page" else "unexpected HTML",
        "). Treating this as a failed listing rather than a legitimate empty directory."
      )
      return(invisible(character(0)))
    }
    
    hrefs <- unique(unlist(regmatches(html, gregexpr('href="[^"]+"', html, perl = TRUE))))
    hrefs <- sub('^href="', '', hrefs)
    hrefs <- sub('"$', '', hrefs)
    candidate_files <- hrefs[grepl(file_pattern, hrefs, ignore.case = TRUE)]
    candidate_files <- candidate_files[!grepl("^\\.\\.?/?$|^\\?|^https?://(?!.*\\bdaacdata\\b)", candidate_files, perl = TRUE)]
    
    if (length(candidate_files) == 0) {
      warning("NSIDC directory listing at ", dir_url, " returned no files matching '", file_pattern, "'.")
      return(invisible(character(0)))
    }
    
    if (!is.null(start_year) && !is.null(end_year)) {
      yr <- suppressWarnings(as.integer(regmatches(
        candidate_files, regexpr("(19|20)\\d{2}", candidate_files)
      )))
      # Keep undated files too (e.g. a mask or readme swept up by the file
      # pattern) rather than dropping them outright - only files that DO have
      # a parseable year and fall outside [start_year, end_year] are excluded.
      keep <- is.na(yr) | (yr >= start_year & yr <= end_year)
      candidate_files <- candidate_files[keep]
    }
    if (length(candidate_files) == 0) {
      warning("NSIDC directory listing at ", dir_url, " had files, but none fell within [",
              start_year, ", ", end_year, "].")
      return(invisible(character(0)))
    }
    
    base <- sub("/[^/]*$", "/", dir_url)
    out <- character(0)
    n_failed <- 0
    for (f in candidate_files) {
      full_url <- if (grepl("^https?://", f)) f else paste0(base, f)
      dest <- file.path(out_dir, basename(f))
      if (!file.exists(dest)) {
        message("  [DL] ", basename(dest))
        ok <- tryCatch({ download_file_auth(full_url, dest); TRUE }, error = function(e) {
          message("    download failed for ", basename(dest), ": ", conditionMessage(e))
          FALSE
        })
        if (!ok) { n_failed <- n_failed + 1; next }
      }
      out <- c(out, dest)
    }
    message(sprintf("  [OK] %d file(s) on disk from %s (%d failed).", length(out), dir_url, n_failed))
    invisible(out)
  }
  
  # Study window from local ERA5-Land SWE, if available
  # ----------------------------------------------------------------------------
  ERA5_LAND_SWE_FILE <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
  study_start <- as.Date("1981-01-01")
  study_end   <- Sys.Date()
  
  if (file.exists(ERA5_LAND_SWE_FILE)) {
    message("[OK] Found local ERA5-Land SWE file: ", ERA5_LAND_SWE_FILE)
    era5_r <- tryCatch(terra::rast(ERA5_LAND_SWE_FILE), error = function(e) NULL)
    if (!is.null(era5_r)) {
      era5_dates <- tryCatch({
        tt <- terra::time(era5_r)
        if (!all(is.na(tt))) as.Date(tt) else NULL
      }, error = function(e) NULL)
      if (is.null(era5_dates) || all(is.na(era5_dates))) {
        era5_dates <- tryCatch({
          nc <- ncdf4::nc_open(ERA5_LAND_SWE_FILE)
          on.exit(ncdf4::nc_close(nc), add = TRUE)
          if ("valid_time" %in% names(nc$var)) {
            tv <- ncdf4::ncvar_get(nc, "valid_time")
            as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC"))
          } else if ("time" %in% names(nc$dim)) {
            tv <- ncdf4::ncvar_get(nc, "time")
            as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC"))
          } else {
            NULL
          }
        }, error = function(e) NULL)
      }
      if (!is.null(era5_dates) && any(!is.na(era5_dates))) {
        study_start <- min(era5_dates, na.rm = TRUE)
        study_end   <- max(era5_dates, na.rm = TRUE)
      }
    }
  } else {
    message("[WARN] Local ERA5-Land SWE file not found; using default study window.")
  }
  message(sprintf("[OK] Study window: %s to %s", study_start, study_end))
  
  # ----------------------------------------------------------------------------
  # Basin geometry (reused by all acquisition sections)
  # ----------------------------------------------------------------------------
  read_kmz_or_kml <- function(path) {
    if (!grepl("\\.kmz$", path, ignore.case = TRUE)) {
      return(sf::st_read(path, quiet = TRUE))
    }
    kml_inside <- utils::unzip(path, list = TRUE)$Name
    kml_inside <- kml_inside[grepl("\\.kml$", kml_inside, ignore.case = TRUE)][1]
    if (is.na(kml_inside) || !nzchar(kml_inside)) stop("No KML found inside ", path)
    vsi_path <- file.path("/vsizip", path, kml_inside)
    result <- tryCatch(sf::st_read(vsi_path, quiet = TRUE), error = function(e) NULL)
    if (!is.null(result)) return(result)
    tmp_dir <- file.path(tempdir(), "kmz_extract")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    utils::unzip(path, files = kml_inside, exdir = tmp_dir)
    sf::st_read(file.path(tmp_dir, kml_inside), quiet = TRUE)
  }
  
  BASIN_KMZ <- "Spatial/Nechako_Basin.kmz"
  basin_sf <- if (exists("basin", inherits = TRUE) && inherits(basin, "sf")) {
    basin
  } else {
    read_kmz_or_kml(BASIN_KMZ)
  }
  basin_wgs84 <- suppressWarnings(sf::st_transform(basin_sf, 4326))
  basin_bbox <- sf::st_bbox(basin_wgs84)
  bbox_vec <- c(xmin = as.numeric(basin_bbox["xmin"]),
                ymin = as.numeric(basin_bbox["ymin"]),
                xmax = as.numeric(basin_bbox["xmax"]),
                ymax = as.numeric(basin_bbox["ymax"]))
  bbox_str <- sprintf("%s,%s,%s,%s", bbox_vec["xmin"], bbox_vec["ymin"], bbox_vec["xmax"], bbox_vec["ymax"])
  
  # ----------------------------------------------------------------------------
  # 1) CanSWE acquisition + merge into daily.csv
  # ----------------------------------------------------------------------------
  acquire_canswe <- function(record_id = 5889352, out_dir = "snow_data/canswe",
                             local_dir = Sys.getenv("CANSWE_LOCAL_DIR", unset = "snow_data/CanSWE-CanEEN_1928-2025_v8")) {
    message("\n===== CanSWE acquisition =====")
    ensure_dir(out_dir)
    
    # ---- Prefer an already-downloaded local copy over hitting Zenodo at all ----
    # CanSWE is distributed as one big versioned bundle (e.g.
    # "CanSWE-CanEEN_1928-2025_v8"), not something worth re-fetching from Zenodo
    # every run once you already have it - and the Zenodo API call in
    # zenodo_record_files() has no retry/timeout around it (see the CMR/ASF
    # calls elsewhere in this script for that pattern), so a transient DNS/
    # network hiccup there (e.g. "Could not resolve hostname") throws and can
    # take down the rest of Part B with it. Checking local_dir first avoids the
    # network call entirely when the data's already on disk. local_dir can be
    # a specific .csv file, or a folder (searched recursively for .csv files) -
    # matching how the downloaded CanSWE bundle actually looks on disk (a
    # folder containing the .csv, a .nc, and the Lisez-moi/ReadMe PDFs).
    csvs <- character(0)
    if (nzchar(local_dir) && file.exists(local_dir)) {
      csvs <- if (dir.exists(local_dir)) {
        list.files(local_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
      } else if (grepl("\\.csv$", local_dir, ignore.case = TRUE)) {
        local_dir
      } else {
        character(0)
      }
    }
    
    if (length(csvs) > 0) {
      message(
        "  [OK] Using local CanSWE data at '", local_dir, "' (", length(csvs),
        " CSV file(s)) - skipping the Zenodo download for record ", record_id, "."
      )
    } else {
      message(
        "  No local CanSWE data found at '", local_dir, "' - falling back to ",
        "downloading record ", record_id, " from Zenodo. (Set CANSWE_LOCAL_DIR, ",
        "or pass local_dir=, to point this at an existing local copy instead.)"
      )
      
      files <- zenodo_record_files(record_id)
      if (nrow(files) == 0) {
        warning("CanSWE Zenodo record returned no files.")
        return(invisible(NULL))
      }
      
      data_files <- files[grepl("\\.(csv|zip|gz)$", files$filename, ignore.case = TRUE), , drop = FALSE]
      if (nrow(data_files) == 0) data_files <- files
      
      downloaded <- character(0)
      for (i in seq_len(nrow(data_files))) {
        fn <- data_files$filename[i]
        url <- data_files$url[i]
        dest <- file.path(out_dir, fn)
        if (!file.exists(dest)) {
          message("  [DL] ", fn)
          download_file_auth(url, dest)
        } else {
          message("  [OK] ", fn, " already exists.")
        }
        downloaded <- c(downloaded, dest)
        
        if (grepl("\\.zip$", fn, ignore.case = TRUE)) {
          utils::unzip(dest, exdir = out_dir)
        }
      }
      
      csvs <- list.files(out_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
    }
    
    if (length(csvs) == 0) {
      stop("No CSV found locally or after downloading CanSWE record ", record_id, ".")
    }
    
    parse_one <- function(csv_path) {
      df <- tryCatch(
        readr::read_csv(csv_path, col_types = readr::cols(.default = "c"), guess_max = Inf, show_col_types = FALSE),
        error = function(e) NULL
      )
      if (is.null(df) || nrow(df) == 0) return(NULL)
      
      orig <- names(df)
      nm <- tolower(gsub("[^a-z0-9]+", "_", orig))
      
      names(df) <- nm
      
      sel_col <- function(nms, pattern) {
        hit <- nms[which(grepl(pattern, nms, ignore.case = TRUE))[1]]
        if (length(hit) == 0 || is.na(hit) || !nzchar(hit)) return(NA_character_)
        hit
      }
      
      # preferred names
      # NOTE: CanSWE's own CSV/NetCDF uses short CF-style variable names, not
      # spelled-out ones - confirmed in Vionnet et al. (2021, ESSD) and the
      # NSIDC/Zenodo description: "snw" = SWE [kg m-2], "snd" = snow depth [m],
      # "den" = bulk density [kg m-3]. The original "swe"/"depth" patterns below
      # never matched "snw"/"snd", so swe_col came back NA, out$swe_mm was all
      # NA, and every row was silently dropped by the is.finite(swe_mm) filter
      # further down - which is exactly the "CSVs found but none could be
      # parsed" warning this produced. "snw"/"snd" are added as exact-name
      # alternatives (anchored with ^...$) so this doesn't start accidentally
      # matching unrelated columns on some other dataset.
      station_col <- sel_col(names(df), "station|site|pillow|id")
      lon_col     <- sel_col(names(df), "lon|long")
      lat_col     <- sel_col(names(df), "lat")
      date_col    <- sel_col(names(df), "^date$|date|time|timestamp|datetime")
      swe_col     <- sel_col(names(df), "swe|^snw$")
      sd_col      <- sel_col(names(df), "snow.*depth|^sd$|depth|^snd$")
      
      if (any(is.na(c(station_col, lon_col, lat_col, date_col)))) {
        message(
          "  [skip] ", basename(csv_path), " - could not identify all required ",
          "columns (station/lon/lat/date) among: ", paste(orig, collapse = ", ")
        )
        return(NULL)
      }
      if (is.na(swe_col)) {
        message(
          "  [warn] ", basename(csv_path), " - no SWE-like column found (looked ",
          "for 'swe' or 'snw') among: ", paste(orig, collapse = ", "),
          " - every row from this file will be dropped."
        )
      }
      
      out <- data.frame(
        station_id = paste0("CanSWE_", normalize_id(df[[station_col]])),
        lon = suppressWarnings(as.numeric(df[[lon_col]])),
        lat = suppressWarnings(as.numeric(df[[lat_col]])),
        date = maybe_parse_date(df[[date_col]]),
        swe_mm = if (!is.na(swe_col)) suppressWarnings(as.numeric(df[[swe_col]])) else NA_real_,
        stringsAsFactors = FALSE
      )
      
      # If SWE is missing but depth exists, keep the record only if no SWE is present.
      if (all(is.na(out$swe_mm)) && !is.na(sd_col)) {
        warning("CSV ", basename(csv_path), " has no SWE column; depth exists but is not converted automatically.")
      }
      
      out <- out[is.finite(out$lon) & is.finite(out$lat) & !is.na(out$date) & is.finite(out$swe_mm), , drop = FALSE]
      out
    }
    
    daily_list <- lapply(csvs, parse_one)
    daily_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, daily_list)
    if (length(daily_list) == 0) {
      warning("CanSWE CSVs were found, but none could be parsed into station/date/SWE records.")
      return(invisible(NULL))
    }
    
    canswe_daily <- dplyr::bind_rows(daily_list)
    canswe_daily <- canswe_daily %>%
      dplyr::distinct(station_id, lon, lat, date, swe_mm, .keep_all = TRUE) %>%
      dplyr::arrange(station_id, date)
    canswe_pull_time <- Sys.time()
    canswe_daily <- canswe_daily %>% dplyr::mutate(canswe_pull_time = canswe_pull_time)
    write_csv(canswe_daily, file.path(out_dir, "CanSWE_daily.csv"))
    message(sprintf("  [OK] Parsed CanSWE records: %d rows, %d stations.",
                    nrow(canswe_daily), dplyr::n_distinct(canswe_daily$station_id)))
    invisible(canswe_daily)
  }
  
  canswe_daily <- acquire_canswe()
  
  # FIX (point 1): daily.csv is NOT overwritten with CanSWE data here anymore.
  # acquire_canswe() (above) already writes the full CanSWE dataset to its own
  # file, snow_data/canswe/CanSWE_daily.csv - unfiltered, national scope (all
  # ~11,000+ CanSWE stations), because CanSWE has no basin filter applied
  # anywhere in this script. Previously this block then bound that entire
  # national dataset onto daily.csv and wrote the result back over
  # DAILY_CSV_PATH - so the file Part A carefully spatial/topo-screened down
  # to a handful of Nechako stations silently ballooned to 11,000+ stations,
  # >99% of them nowhere near the basin, with no warning.
  #
  # H8SWEI_Snow_2_regionalization_bias.R's merge_station_networks() already
  # reads daily.csv AND CanSWE_daily.csv as two separate inputs, and does the
  # real work properly there: it filters EVERY station (from either file) to
  # ones that resolve to a valid (non-NA) cell in the basin's own SWE grid
  # template, and dedupes CanSWE rows by station_id+date identity (not by
  # exact-value match). That is the correct place for CanSWE to be spatially
  # scoped and merged - doing it again here duplicated that work with none of
  # its safeguards, and is exactly the fragile "double-merge" scenario
  # Nechako_Snow_Bridge.R's STEP 1 (verify_canswe_dedup) already checks for
  # and warns about ("remove the CanSWE merge-into-daily.csv step ... and let
  # H8 be the ONLY place merging" - see that script's own comments).
  #
  # daily.csv therefore now stays exactly what Part A produced: ASWS stations
  # screened against the Nechako basin/buffer. CanSWE_daily.csv is still
  # written in full (H8 needs the complete file to do its own filtering), and
  # a preview/audit copy of what a combined view would look like is still
  # written to a clearly-separate sidecar file, in case it's useful for
  # inspection - it is NOT read by any downstream script.
  if (!file.exists(DAILY_CSV_PATH)) {
    warning("daily.csv from Part A was not found; nothing to preview alongside CanSWE.")
  } else if (!is.null(canswe_daily) && nrow(canswe_daily) > 0) {
    asws_daily <- readr::read_csv(DAILY_CSV_PATH, show_col_types = FALSE)
    if (!all(c("station_id", "lon", "lat", "date", "swe_mm") %in% names(asws_daily))) {
      stop("Existing daily.csv does not contain the required columns.")
    }
    asws_daily <- asws_daily %>%
      dplyr::mutate(date = as.Date(date), station_id = as.character(station_id)) %>%
      dplyr::mutate(source = "ASWS")
    canswe_preview <- canswe_daily %>% dplyr::mutate(source = "CanSWE")
    combined_preview <- dplyr::bind_rows(asws_daily, canswe_preview) %>%
      dplyr::arrange(station_id, date, source)
    write_csv(combined_preview, "snow_data/stations_daily_with_source_PREVIEW_ONLY.csv")
    message(sprintf(
      "[OK] daily.csv left as ASWS-only (%d rows, %d station(s)) - NOT merged with CanSWE's %d rows/%d stations. A combined preview (for inspection only, not read by H8 or Bridge) was written to snow_data/stations_daily_with_source_PREVIEW_ONLY.csv. H8SWEI_Snow_2's merge_station_networks() performs the real, spatially-filtered CanSWE merge at run time.",
      nrow(asws_daily), dplyr::n_distinct(asws_daily$station_id),
      nrow(canswe_daily), dplyr::n_distinct(canswe_daily$station_id)
    ))
  } else {
    message("[OK] CanSWE unavailable this run; daily.csv left as ASWS-only (unchanged from Part A).")
  }
  
  # ----------------------------------------------------------------------------
  # 2) Sentinel-1 SAR acquisition (ASF) + optional RCM/C-SNOW manifests
  # ----------------------------------------------------------------------------
  acquire_sar_manifest <- function(manifest_csv, out_dir) {
    if (!file.exists(manifest_csv)) return(invisible(character(0)))
    ensure_dir(out_dir)
    man <- readr::read_csv(manifest_csv, show_col_types = FALSE)
    if (!("url" %in% names(man))) stop("SAR manifest ", manifest_csv, " must contain a 'url' column.")
    if (!("file_name" %in% names(man))) man$file_name <- basename(sub("\\?.*$", "", man$url))
    out <- character(0)
    for (i in seq_len(nrow(man))) {
      dest <- file.path(out_dir, man$file_name[i])
      if (!file.exists(dest)) {
        message("  [DL] ", basename(dest))
        download_file_auth(man$url[i], dest)
      }
      out <- c(out, dest)
    }
    invisible(out)
  }
  
  # ----------------------------------------------------------------------------
  # Copernicus Data Space Ecosystem (CDSE) - fallback Sentinel-1 source
  # ----------------------------------------------------------------------------
  # ASF and CDSE are independent infrastructure serving the same underlying
  # Sentinel-1 archive, so an ASF-side outage (the 504s seen in prior runs)
  # doesn't affect CDSE at all. This is used ONLY as a fallback when ASF's
  # search fails outright or returns nothing - not tried first - since the
  # existing ASF path is already what the rest of this pipeline expects
  # (same filenames, same out_dir).
  #
  # Needs CDSE_CLIENT_ID / CDSE_CLIENT_SECRET in .env - an OAuth client
  # created under an existing Copernicus Data Space Ecosystem / Sentinel Hub
  # account (https://shapps.dataspace.copernicus.eu/dashboard/#/account/settings
  # -> "OAuth clients" -> Create). This uses the OAuth2 client_credentials
  # grant, which is what a Sentinel Hub OAuth client (client_id starting with
  # "sh-...") actually authenticates with - NOT a username/password login.
  #
  # UPDATE: this used to send CDS_USER/CDS_PASS with grant_type = "password"
  # against identity.dataspace.copernicus.eu - that endpoint expects the
  # email+password of an account registered directly at
  # https://dataspace.copernicus.eu/, which is a THIRD, separate credential
  # set from both (a) an ECMWF/Climate Data Store login (cds.climate.
  # copernicus.eu - a different service entirely, despite the similar name)
  # and (b) a Sentinel Hub OAuth client's id+secret (what's actually
  # available on the account dashboard screenshot this was debugged from).
  # Sending either of the wrong credential pairs to the password-grant
  # endpoint reliably 401s - it's not a network/script problem, it's the
  # wrong credentials for that specific grant type. Switched to
  # client_credentials so this matches the OAuth client that's actually easy
  # to obtain from the CDSE/Sentinel Hub dashboard.
  #
  # If those aren't set, this fallback just warns and returns no files,
  # rather than erroring out.
  
  # CDSE access tokens are short-lived (~10 minutes) - see
  # https://documentation.dataspace.copernicus.eu/APIs/Token.html - so this
  # is called fresh at the start of a download batch and re-called partway
  # through a long one, rather than once for the whole acquisition.
  get_cdse_token <- function(client_id = Sys.getenv("CDSE_CLIENT_ID", unset = ""),
                             client_secret = Sys.getenv("CDSE_CLIENT_SECRET", unset = ""),
                             timeout_sec = 30) {
    if (!nzchar(client_id) || !nzchar(client_secret)) {
      # NOTE: this used to `return(NULL)` here with no message at all - the
      # caller (acquire_sentinel1_cdse) DOES explain this case, but only via
      # warning(), which R buffers/defers by default (printed as a summary at
      # the very end of the script, or not at all if the script errors out or
      # is interrupted first) rather than printing immediately like message()/
      # cat() do. The net effect in a prior run's log was the
      # "===== Sentinel-1 (Copernicus Data Space Ecosystem) fallback ====="
      # header (a message(), which DOES print immediately) followed by total
      # silence - no confirmation of success, failure, or reason - because the
      # one line that explained why was sitting in R's warning buffer instead
      # of the console. Printing this immediately here removes the silent gap
      # regardless of what happens later in the script or whether warnings are
      # ever flushed.
      message(
        "  CDSE fallback: CDSE_CLIENT_ID/CDSE_CLIENT_SECRET not set in .env ",
        "(or empty) - skipping. Create an OAuth client at ",
        "https://shapps.dataspace.copernicus.eu/dashboard/#/account/settings ",
        "(under 'OAuth clients') and set both to enable this fallback."
      )
      return(NULL)
    }
    resp <- tryCatch(
      httr::POST(
        "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
        body = list(client_id = client_id, client_secret = client_secret,
                    grant_type = "client_credentials"),
        encode = "form", timeout(timeout_sec)
      ),
      error = function(e) {
        message("  CDSE token request failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(resp) || httr::http_error(resp)) {
      if (!is.null(resp)) {
        message("  CDSE token request returned HTTP ", httr::status_code(resp), " - ",
                "check CDSE_CLIENT_ID/CDSE_CLIENT_SECRET in .env (the client secret ",
                "is only shown once when the OAuth client is created/regenerated - ",
                "regenerate it on the dashboard if it's been lost).")
      }
      return(NULL)
    }
    tok <- tryCatch(httr::content(resp, as = "parsed")$access_token, error = function(e) NULL)
    if (is.null(tok) || !nzchar(tok)) return(NULL)
    tok
  }
  
  # Pages through CDSE's OData catalogue (via @odata.nextLink) for Sentinel-1
  # GRD products intersecting bbox within [start, end]. Retries each page on
  # timeout/5xx the same way cmr_query_page()/the ASF search above do.
  query_cdse_sentinel1 <- function(bbox, start, end, top = 1000,
                                   query_timeout_sec = 60, query_retries = 5) {
    wkt <- sprintf(
      "POLYGON((%1$f %2$f, %3$f %2$f, %3$f %4$f, %1$f %4$f, %1$f %2$f))",
      bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]
    )
    filter_str <- sprintf(
      paste0(
        "Collection/Name eq 'SENTINEL-1' and contains(Name,'GRD') and ",
        "OData.CSC.Intersects(area=geography'SRID=4326;%s') and ",
        "ContentDate/Start gt %sT00:00:00.000Z and ContentDate/Start lt %sT23:59:59.999Z"
      ),
      wkt, start, end
    )
    next_url <- paste0(
      "https://catalogue.dataspace.copernicus.eu/odata/v1/Products",
      "?$filter=", utils::URLencode(filter_str, reserved = TRUE),
      "&$top=", top, "&$orderby=ContentDate/Start%20asc"
    )
    
    all_rows <- list()
    page <- 1
    repeat {
      resp <- NULL
      for (attempt in seq_len(query_retries)) {
        resp <- tryCatch(
          httr::GET(next_url, timeout(query_timeout_sec)),
          error = function(e) {
            message(sprintf("  CDSE catalogue query (page %d) attempt %d/%d failed: %s",
                            page, attempt, query_retries, conditionMessage(e)))
            NULL
          }
        )
        if (!is.null(resp) && !httr::http_error(resp)) break
        if (!is.null(resp)) {
          message(sprintf(
            "  CDSE catalogue query (page %d) attempt %d/%d returned HTTP %d - %s",
            page, attempt, query_retries, httr::status_code(resp),
            if (attempt < query_retries) "retrying..." else "giving up."
          ))
        }
        resp <- NULL
        if (attempt < query_retries) Sys.sleep(10 * attempt)
      }
      if (is.null(resp)) {
        message("  [SKIP] CDSE catalogue query failed after ", query_retries, " attempts on page ", page, " - giving up on this page.")
        break
      }
      
      parsed <- tryCatch(httr::content(resp, as = "parsed"), error = function(e) NULL)
      vals <- parsed$value
      if (is.null(vals) || length(vals) == 0) break
      
      rows <- lapply(vals, function(v) data.frame(id = v$Id, name = v$Name, stringsAsFactors = FALSE))
      all_rows <- c(all_rows, rows)
      
      next_url <- parsed[["@odata.nextLink"]]
      if (is.null(next_url) || !nzchar(next_url)) break
      page <- page + 1
    }
    
    if (length(all_rows) == 0) return(data.frame(id = character(0), name = character(0)))
    dplyr::bind_rows(all_rows)
  }
  
  # Downloads a single product via CDSE's "zipper" endpoint. Returns TRUE/FALSE
  # rather than throwing, so the caller can move on to the next product instead
  # of one bad download stopping the whole batch.
  download_cdse_product <- function(id, name, dest, token, timeout_sec = 900) {
    url <- sprintf("https://zipper.dataspace.copernicus.eu/odata/v1/Products(%s)/$value", id)
    resp <- tryCatch(
      httr::GET(
        url,
        httr::add_headers(Authorization = paste("Bearer", token)),
        httr::write_disk(dest, overwrite = TRUE),
        timeout(timeout_sec)
      ),
      error = function(e) {
        message("    CDSE download failed for ", name, ": ", conditionMessage(e))
        NULL
      }
    )
    ok <- !is.null(resp) && !httr::http_error(resp)
    if (!ok) {
      if (!is.null(resp)) message("    CDSE download HTTP ", httr::status_code(resp), " for ", name)
      if (file.exists(dest)) unlink(dest)
    }
    ok
  }
  
  # Full CDSE fallback: token -> catalogue query -> per-product download.
  SENTINEL1_LAUNCH_DATE <- as.Date("2014-04-03")
  
  # Writes into the SAME out_dir as acquire_sentinel1_asf() so downstream
  # code doesn't need to know which source a given file came from.
  acquire_sentinel1_cdse <- function(out_dir = "snow_data/sar/sentinel1",
                                     start = study_start, end = study_end,
                                     max_products = 5000) {
    message("\n===== Sentinel-1 (Copernicus Data Space Ecosystem) fallback =====")
    ensure_dir(out_dir)
    
    if (start < SENTINEL1_LAUNCH_DATE) {
      message(sprintf(
        "  [INFO] Requested Sentinel-1 start %s predates mission launch; clamping to %s.",
        format(start, "%Y-%m-%d"), format(SENTINEL1_LAUNCH_DATE, "%Y-%m-%d")
      ))
      start <- SENTINEL1_LAUNCH_DATE
    }
    
    token <- get_cdse_token()
    if (is.null(token)) {
      # get_cdse_token() already printed WHY (missing credentials, failed
      # request, bad HTTP status) immediately via message() - this just closes
      # out the section with a clear, immediately-visible outcome instead of
      # a deferred warning() that could otherwise leave this section trailing
      # off with no visible result in the log.
      message("  [SKIP] Sentinel-1 CDSE fallback: no data acquired (no valid access token - see reason above).")
      return(invisible(character(0)))
    }
    
    products <- query_cdse_sentinel1(bbox_vec, format(start, "%Y-%m-%d"), format(end, "%Y-%m-%d"),
                                     top = max_products)
    if (nrow(products) == 0) {
      message("  [SKIP] Sentinel-1 CDSE fallback: catalogue query returned 0 Sentinel-1 GRD products for this bbox/time window.")
      return(invisible(character(0)))
    }
    message(sprintf("  [OK] CDSE catalogue found %d Sentinel-1 GRD product(s).", nrow(products)))
    
    MIN_VALID_BYTES <- 50 * 1024
    is_valid <- function(f) file.exists(f) && file.info(f)$size >= MIN_VALID_BYTES
    
    out <- character(0)
    n_failed <- 0
    token_issued <- Sys.time()
    
    for (i in seq_len(nrow(products))) {
      dest <- file.path(out_dir, paste0(products$name[i], ".zip"))
      if (is_valid(dest)) {
        out <- c(out, dest)
        next
      }
      # Refresh proactively at ~9 minutes rather than waiting for a 401
      # partway through what can be a long batch of large SAR downloads.
      if (difftime(Sys.time(), token_issued, units = "secs") > 540) {
        refreshed <- get_cdse_token()
        if (!is.null(refreshed)) { token <- refreshed; token_issued <- Sys.time() }
      }
      message(sprintf("  [DL %d/%d] %s", i, nrow(products), basename(dest)))
      ok <- download_cdse_product(products$id[i], products$name[i], dest, token)
      if (ok && is_valid(dest)) {
        out <- c(out, dest)
      } else {
        n_failed <- n_failed + 1
      }
    }
    message(sprintf("  [OK] CDSE fallback: %d file(s) on disk, %d failed.", length(out), n_failed))
    invisible(out)
  }
  
  acquire_sentinel1_asf <- function(out_dir = "snow_data/sar/sentinel1",
                                    start = study_start, end = study_end,
                                    max_results = 5000) {
    message("\n===== Sentinel-1 (ASF) acquisition =====")
    ensure_dir(out_dir)
    
    if (start < SENTINEL1_LAUNCH_DATE) {
      message(sprintf(
        "  [INFO] Requested Sentinel-1 start %s predates mission launch; clamping to %s.",
        format(start, "%Y-%m-%d"), format(SENTINEL1_LAUNCH_DATE, "%Y-%m-%d")
      ))
      start <- SENTINEL1_LAUNCH_DATE
    }
    
    q <- list(
      dataset = "SENTINEL-1",
      bbox = bbox_str,
      start = format(start, "%Y-%m-%d"),
      end = format(end, "%Y-%m-%d"),
      processingLevel = "GRD_HD",
      output = "download",
      maxResults = as.integer(max_results)
    )
    
    # Retry with backoff: a bare single-attempt GET here previously meant a
    # transient timeout/5xx from the ASF search endpoint (observed: HTTP 504,
    # a server-side timeout - not evidence of "no data") skipped Sentinel-1
    # acquisition for the whole run. Mirrors the retry pattern already used in
    # cmr_query_page() for the same class of problem.
    search_retries <- 5
    resp <- NULL
    for (attempt in seq_len(search_retries)) {
      resp <- tryCatch(
        httr::GET("https://api.daac.asf.alaska.edu/services/search/param", query = q, timeout(300)),
        error = function(e) {
          message(sprintf("  ASF search attempt %d/%d failed: %s", attempt, search_retries, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(resp) && !httr::http_error(resp)) break
      if (!is.null(resp)) {
        message(sprintf(
          "  ASF search attempt %d/%d returned HTTP %d (often transient, e.g. a server-side timeout) - %s",
          attempt, search_retries, httr::status_code(resp),
          if (attempt < search_retries) "retrying..." else "giving up."
        ))
      }
      resp <- NULL
      if (attempt < search_retries) Sys.sleep(15 * attempt)
    }
    if (is.null(resp)) {
      warning(
        "ASF search request failed after ", search_retries, " attempts - ",
        "falling back to the Copernicus Data Space Ecosystem (independent ",
        "infrastructure, unaffected by an ASF-side outage)."
      )
      return(acquire_sentinel1_cdse(out_dir = out_dir, start = start, end = end))
    }
    
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    urls <- unique(unlist(regmatches(txt, gregexpr("https?://[^[:space:]\"']+", txt, perl = TRUE))))
    urls <- urls[grepl("\\.(zip|SAFE|tif|tiff)$|download|data", urls, ignore.case = TRUE)]
    if (length(urls) == 0) {
      warning(
        "ASF search returned no downloadable URLs - falling back to the ",
        "Copernicus Data Space Ecosystem."
      )
      return(acquire_sentinel1_cdse(out_dir = out_dir, start = start, end = end))
    }
    
    out <- character(0)
    for (u in urls) {
      fn <- basename(sub("\\?.*$", "", u))
      if (!nzchar(fn) || identical(fn, "/")) fn <- paste0("asf_", sprintf("%08d", as.integer(abs(sum(charToRaw(u))) %% 1e8)), ".dat")
      dest <- file.path(out_dir, fn)
      if (!file.exists(dest)) {
        message("  [DL] ", fn)
        download_file_auth(u, dest)
      }
      out <- c(out, dest)
    }
    message(sprintf("  [OK] Sentinel-1 acquisitions stored in %s (%d files).", out_dir, length(out)))
    invisible(out)
  }
  
  # Optional local/manifest-based support for RCM / C-SNOW archives.
  acquire_rcm_or_csnow <- function(provider_name, manifest_csv = NULL, local_dir = NULL) {
    out_dir <- file.path("snow_data", "sar", tolower(provider_name))
    ensure_dir(out_dir)
    if (!is.null(manifest_csv) && file.exists(manifest_csv)) {
      message(sprintf("\n===== %s manifest acquisition =====", provider_name))
      acquire_sar_manifest(manifest_csv, out_dir)
    }
    if (!is.null(local_dir) && dir.exists(local_dir)) {
      message(sprintf("\n===== %s local ingest =====", provider_name))
      files <- list.files(local_dir, full.names = TRUE, recursive = TRUE,
                          pattern = "\\.(zip|SAFE|tif|tiff|nc)$", ignore.case = TRUE)
      for (f in files) {
        dest <- file.path(out_dir, basename(f))
        if (!file.exists(dest)) file.copy(f, dest, overwrite = FALSE)
      }
      message(sprintf("  [OK] Copied %d %s file(s) into %s.", length(files), provider_name, out_dir))
    }
    invisible(list.files(out_dir, full.names = TRUE))
  }
  
  # ----------------------------------------------------------------------------
  # C-SNOW local ingestion (Lievens et al. Sentinel-1 snow depth)
  # ----------------------------------------------------------------------------
  # Helper: find the actual snow-depth variable name inside one .nc file.
  # The variable name INSIDE these NetCDFs isn't confirmed here - only the file
  # NAMES ("snd_YYYYMMDD.nc") are known from the local directory listing.
  # Rather than hardcode a guessed variable name (risking a silent all-NA
  # extraction if it's wrong), this opens one file, lists its actual data
  # variables (excluding obvious dimension/coordinate variables), and picks the
  # most plausible snow-depth-like one - same "detect, don't guess" pattern
  # already used above for the WFS station-ID column (Part A) and CanSWE's
  # swe/snd columns (Part B).
  detect_csnow_varname <- function(nc_path) {
    nc <- tryCatch(ncdf4::nc_open(nc_path), error = function(e) NULL)
    if (is.null(nc)) return(NA_character_)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    
    all_vars <- names(nc$var)
    dim_like <- c("lon", "longitude", "lat", "latitude", "x", "y", "time",
                  "crs", "spatial_ref", "transverse_mercator", "crs_wkt")
    candidates <- all_vars[!tolower(all_vars) %in% dim_like]
    if (length(candidates) == 0) return(NA_character_)
    
    preferred <- candidates[grepl("^snd$|snow.*depth|^sd$|^depth$", candidates, ignore.case = TRUE)]
    if (length(preferred) > 0) return(preferred[1])
    # Fall back to whichever real (non-dimension) variable the file has - most
    # single-product NetCDFs like this one carry only one data variable anyway.
    candidates[1]
  }
  
  # ----------------------------------------------------------------------------
  # ROBUST C-SNOW READER - replaces the old plain-vs-subds terra::rast()
  # ordering "fix", which turned out not to be the actual bug.
  # ----------------------------------------------------------------------------
  # DIAGNOSIS (from inspecting a sample file, snd_20170807.nc, directly with
  # h5py/ncdf4 rather than guessing from the warning text alone):
  #   - lat/lon in this archive ARE proper 1-D NetCDF4 dimension-scale
  #     variables (not 2-D auxiliary coordinates on a projected grid, as
  #     previously assumed from the warning wording alone).
  #   - The grid IS regular: cellsize computed from (max-min)/(n-1) is
  #     ~0.0049603 deg in both lat and lon, matching to 5+ significant
  #     figures. But the raw lat/lon vectors carry float64-vs-float32
  #     rounding jitter (successive diffs vary by ~2e-6 deg, i.e. 4 distinct
  #     rounded diff values instead of exactly 1) - almost certainly an
  #     artifact of whatever regridding/export pipeline produced these files.
  #   - terra's classic NetCDF reader checks dimension spacing against a
  #     tight tolerance and this jitter fails that check, so it treats the
  #     file as irregular/multidimensional and skips the array entirely:
  #       "[rast] skipped multidimensional array: /snd (lat is not regularly spaced)"
  #     This happens for BOTH terra::rast(f) and terra::rast(f, subds=...) -
  #     they hit the same regularity check - which is why reordering plain
  #     vs subds (the previous "FIX (point 2)") did not actually help; it
  #     was solving the wrong problem. subds=/mdim was never the issue here.
  #
  # FIX: bypass terra's NetCDF dimension-regularity check entirely by
  # reading the raw lat/lon/data arrays with ncdf4 and constructing the
  # SpatRaster manually from the array plus an extent computed from
  # min/max lat/lon (not from the noisy per-step diffs). This is safe
  # specifically because we've confirmed the underlying grid is regular to
  # ~1e-6 deg - i.e. we are working around a false-positive irregularity
  # check, not actually handling an irregular/curvilinear grid. If a future
  # file in this archive turns out to be genuinely irregular (real
  # curvilinear lat/lon, non-constant cellsize beyond float jitter), this
  # will produce a silently misaligned raster rather than an error - hence
  # the tolerance check below, which stops loudly instead of guessing.
  read_csnow_raster <- function(nc_path, varname) {
    # Try the fast paths first (works fine for any file that doesn't hit the
    # jitter issue above) before paying for a manual ncdf4 read.
    # UPDATE: wrapped in suppressWarnings() - on this archive, EVERY file hits
    # terra's dimension-regularity check and prints
    #   "[rast] skipped multidimensional array: /flag (lat is not regularly spaced)"
    #   "skipped multidimensional array: /snd (lat is not regularly spaced)"
    # twice per file (once here, once again from the per-tile loop call), so
    # across 1000+ files that's 2000+ near-duplicate warning lines - all of
    # them already expected and already handled by the manual ncdf4 fallback
    # below (see the DIAGNOSIS comment above), not indicating anything new
    # wrong. tryCatch(..., error=) alone can't catch these (they're warnings,
    # not errors - that's WHY the plain-read result is silently NULL instead
    # of throwing) - suppressWarnings() is what actually silences them. This
    # does NOT hide a real failure: if BOTH attempts still return NULL, this
    # function falls through to the manual fallback exactly as before, and
    # that fallback still messages [OK]/[SKIP] per file either way.
    r <- suppressWarnings(tryCatch(terra::rast(nc_path), error = function(e) NULL))
    if (is.null(r)) r <- suppressWarnings(tryCatch(terra::rast(nc_path, subds = varname), error = function(e) NULL))
    if (!is.null(r)) return(r)
    
    # Manual fallback: read raw arrays and build the SpatRaster ourselves.
    nc <- tryCatch(ncdf4::nc_open(nc_path), error = function(e) NULL)
    if (is.null(nc)) return(NULL)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    
    lat <- tryCatch(ncdf4::ncvar_get(nc, "lat"), error = function(e) NULL)
    lon <- tryCatch(ncdf4::ncvar_get(nc, "lon"), error = function(e) NULL)
    dat <- tryCatch(ncdf4::ncvar_get(nc, varname), error = function(e) NULL)
    if (is.null(lat) || is.null(lon) || is.null(dat)) return(NULL)
    if (!is.matrix(dat) || length(dim(dat)) != 2) return(NULL)
    
    # Confirm the grid is regular enough to treat as a simple raster (see
    # DIAGNOSIS above) rather than silently misaligning a genuinely
    # irregular grid. Tolerance is generous relative to the ~2e-6 deg
    # jitter actually observed, tight relative to the ~0.00496 deg cellsize.
    lat_diffs <- diff(lat); lon_diffs <- diff(lon)
    lat_cellsize <- (max(lat) - min(lat)) / (length(lat) - 1)
    lon_cellsize <- (max(lon) - min(lon)) / (length(lon) - 1)
    REGULARITY_TOL <- 1e-4  # degrees; observed jitter was ~2e-6, cellsize ~0.005
    if (any(abs(abs(lat_diffs) - lat_cellsize) > REGULARITY_TOL) ||
        any(abs(abs(lon_diffs) - lon_cellsize) > REGULARITY_TOL)) {
      message(
        "  [SKIP] '", basename(nc_path), "': lat/lon spacing is irregular ",
        "beyond float-precision jitter (exceeds ", REGULARITY_TOL, " deg tolerance) - ",
        "this looks like a genuinely non-regular/curvilinear grid, not just the ",
        "rounding jitter this fallback is designed to work around. Refusing to ",
        "guess at an extent for it rather than silently misaligning the raster."
      )
      return(NULL)
    }
    
    # Build the extent from min/max (not from the noisy per-step diffs).
    ext_r <- terra::ext(min(lon) - lon_cellsize / 2, max(lon) + lon_cellsize / 2,
                        min(lat) - lat_cellsize / 2, max(lat) + lat_cellsize / 2)
    # Don't assume a fixed dim order from ncvar_get() - verify against the
    # actual array shape instead (ncdf4's returned order matches the
    # variable's declared dimension order, which isn't guaranteed to be
    # [lon, lat] vs [lat, lon] the same way for every product/exporter).
    # terra::rast(matrix=) wants rows = lat (nrow = length(lat)), cols = lon.
    if (identical(dim(dat), c(length(lon), length(lat)))) {
      mat <- t(dat)   # was [lon, lat] -> transpose to [lat, lon]
    } else if (identical(dim(dat), c(length(lat), length(lon)))) {
      mat <- dat       # already [lat, lon]
    } else {
      message(
        "  [SKIP] '", basename(nc_path), "': data array dims ", 
        paste(dim(dat), collapse = " x "), " don't match length(lat)=", length(lat),
        " / length(lon)=", length(lon), " in either order - refusing to guess orientation."
      )
      return(NULL)
    }
    # terra expects row 1 = northernmost row (north-up); flip if lat is ascending.
    if (lat[1] < lat[length(lat)]) mat <- mat[nrow(mat):1, , drop = FALSE]
    
    fill_val <- tryCatch(ncdf4::ncatt_get(nc, varname, "_FillValue")$value, error = function(e) NA)
    if (!is.na(fill_val)) mat[mat == fill_val] <- NA
    
    r <- tryCatch(
      terra::rast(mat, extent = ext_r, crs = "EPSG:4326"),
      error = function(e) NULL
    )
    if (!is.null(r)) {
      message("  [OK] '", basename(nc_path), "' opened via manual ncdf4 fallback (terra's regularity check false-positived on float jitter in lat/lon).")
    }
    r
  }
  
  # Ingests the local daily C-SNOW archive (one NetCDF per day, one band each)
  # rather than downloading/reprocessing raw Sentinel-1: extracts (1) a
  # station-level time series at every station_id/lon/lat already known to
  # this pipeline (from daily.csv, i.e. ASWS+CanSWE) - this is exactly what's
  # needed to validate C-SNOW against ground truth you already have, before
  # trusting it further - and (2) a basin-wide daily mean + valid-pixel
  # fraction, which is a cheap way to see how much of the basin C-SNOW
  # actually covers on a given date (C-band SAR snow-depth retrieval is masked
  # out entirely once snow is wet, so coverage should be expected to drop off
  # sharply in spring regardless of anything being wrong).
  acquire_csnow_local <- function(source_dir = CSNOW_LOCAL_DIR,
                                  out_dir = file.path("snow_data", "csnow"),
                                  start = study_start, end = study_end,
                                  station_lookup = NULL,
                                  force_reprocess = isTRUE(as.logical(Sys.getenv("CSNOW_FORCE_REPROCESS", "FALSE")))) {
    message("\n===== C-SNOW (local, Lievens Sentinel-1 snow depth) ingestion =====")
    ensure_dir(out_dir)
    
    if (!dir.exists(source_dir)) {
      # NOTE: warning() alone is not enough here - R buffers/defers warnings
      # by default (see the options(warn = 1) fix near the top of this
      # script), so on some setups this reason would only ever surface as a
      # batch at the very end of the run, well after the "[SKIP]" outcome was
      # already decided elsewhere. message() below prints immediately,
      # exactly where the failure happens, regardless of the warn option or
      # whether the run is later interrupted; warning() is kept alongside it
      # so this still shows up in warnings() and any automated log scraping
      # that greps for "Warning message".
      message(
        "  [SKIP] C-SNOW local directory not found: '", source_dir, "'. ",
        "Set CSNOW_LOCAL_DIR (env var, or edit the default above) if the ",
        "folder lives somewhere else."
      )
      warning(
        "C-SNOW local directory not found: '", source_dir, "' - skipping. ",
        "Set CSNOW_LOCAL_DIR (env var, or edit the default above) if the ",
        "folder lives somewhere else.",
        call. = FALSE
      )
      return(invisible(NULL))
    }
    
    nc_files <- list.files(source_dir, pattern = "^snd_[0-9]{8}\\.nc$",
                           full.names = TRUE, recursive = TRUE)
    if (length(nc_files) == 0) {
      # Same reasoning as above: print immediately via message() rather than
      # relying solely on a deferred warning(). Also report what actually IS
      # in source_dir (up to a handful of names) so a naming-pattern mismatch
      # (e.g. files present but not matching "snd_YYYYMMDD.nc" exactly - a
      # different case, an extra suffix, a subfolder one level deeper than
      # expected) is visible right away instead of requiring a separate manual
      # dir() call to diagnose.
      any_files <- list.files(source_dir, full.names = FALSE, recursive = TRUE)
      sample_note <- if (length(any_files) == 0) {
        "The directory exists but appears to be completely empty."
      } else {
        sprintf(
          "The directory is not empty (%d file(s) total) but none match 'snd_YYYYMMDD.nc'. First few entries found: %s",
          length(any_files), paste(utils::head(any_files, 5), collapse = ", ")
        )
      }
      message(
        "  [SKIP] No 'snd_YYYYMMDD.nc' files found under '", source_dir, "'. ", sample_note
      )
      warning(
        "No 'snd_YYYYMMDD.nc' files found under '", source_dir, "' - skipping. ",
        "Check the folder actually contains the daily NetCDFs (not a further ",
        "subfolder) and that the naming matches this pattern.",
        call. = FALSE
      )
      return(invisible(NULL))
    }
    
    file_dates <- as.Date(sub("^snd_([0-9]{8})\\.nc$", "\\1", basename(nc_files)), format = "%Y%m%d")
    keep <- !is.na(file_dates) & file_dates >= start & file_dates <= end
    n_out_of_window <- sum(!keep, na.rm = TRUE)
    nc_files <- nc_files[keep]
    file_dates <- file_dates[keep]
    
    if (length(nc_files) == 0) {
      message(
        "  [SKIP] Found C-SNOW files under '", source_dir, "', but none fall within the ",
        "study window ", start, " to ", end, " (", n_out_of_window, " file(s) excluded)."
      )
      warning(
        "Found C-SNOW files under '", source_dir, "', but none fall within the ",
        "study window ", start, " to ", end, " (", n_out_of_window, " file(s) excluded).",
        call. = FALSE
      )
      return(invisible(NULL))
    }
    n_files <- length(nc_files)
    message(sprintf(
      "  [OK] Found %d C-SNOW file(s) in the study window (%d outside it, excluded).",
      n_files, n_out_of_window
    ))
    
    station_csv_path <- file.path(out_dir, "CSNOW_station_extract.csv")
    basin_csv_path   <- file.path(out_dir, "CSNOW_basin_summary.csv")
    
    # ------------------------------------------------------------------------
    # RESUMABILITY: skip dates already processed and persisted in a prior run
    # ------------------------------------------------------------------------
    # Every other acquisition function in this script (CanSWE, Sentinel-1,
    # CMC/CanSISE, OpenLandMap) already skips whatever it finds already on
    # disk from a previous run. This loop used to be the one exception: it
    # held ALL of its output in memory (station_rows/basin_rows) and wrote it
    # out only ONCE, after every single file had been processed. If the run
    # was interrupted or paused partway through - the exact situation this
    # section exists to fix - nothing was saved, and the next run started
    # again at file 1 and reprocessed every date from scratch, including ones
    # that had already finished.
    #
    # CSNOW_basin_summary.csv is written for every file that opens
    # successfully, regardless of whether station-level extraction produced
    # anything, so its `date` column is the authoritative "already processed"
    # record. Results are now ALSO written incrementally as the loop runs
    # (see flush_results() below), so an interruption only loses progress
    # back to the last flush (at most FLUSH_EVERY files), not the whole run.
    #
    # Set force_reprocess = TRUE (or env var CSNOW_FORCE_REPROCESS=TRUE) to
    # ignore existing output and reprocess everything from scratch instead -
    # e.g. after changing the extraction logic itself, or after
    # station_lookup has grown to include stations that earlier "already
    # done" dates never saw. This deletes the existing CSVs first so
    # reprocessed dates don't end up duplicated alongside the old rows.
    if (force_reprocess) {
      message("  [INFO] force_reprocess is TRUE - deleting any existing C-SNOW output and reprocessing every file in the study window.")
      if (file.exists(station_csv_path)) unlink(station_csv_path)
      if (file.exists(basin_csv_path)) unlink(basin_csv_path)
    }
    
    already_done_dates <- as.Date(character(0))
    if (!force_reprocess && file.exists(basin_csv_path)) {
      prior_basin <- tryCatch(
        readr::read_csv(basin_csv_path, show_col_types = FALSE),
        error = function(e) NULL
      )
      if (!is.null(prior_basin) && "date" %in% names(prior_basin)) {
        already_done_dates <- unique(as.Date(prior_basin$date))
        message(sprintf(
          "  [OK] %d date(s) already processed in a previous run (found in '%s') - these will be skipped. Set force_reprocess=TRUE (or env var CSNOW_FORCE_REPROCESS=TRUE) to redo them.",
          length(already_done_dates), basin_csv_path
        ))
      } else if (is.null(prior_basin)) {
        warning(
          "'", basin_csv_path, "' exists but could not be read (possibly ",
          "truncated by an earlier hard interruption) - treating it as if no ",
          "prior progress exists. If this is wrong, inspect the file, or set ",
          "force_reprocess=TRUE to rebuild it cleanly."
        )
      }
    }
    
    skip_mask <- file_dates %in% already_done_dates
    n_skip_already_done <- sum(skip_mask)
    nc_files   <- nc_files[!skip_mask]
    file_dates <- file_dates[!skip_mask]
    n_files_todo <- length(nc_files)
    
    if (n_skip_already_done > 0) {
      message(sprintf(
        "  [OK] Skipping %d/%d file(s) already processed in a previous run - %d file(s) left to process this run.",
        n_skip_already_done, n_files, n_files_todo
      ))
    }
    
    # ---- Station lookup: where to extract point values from ----
    # Reads fresh from daily.csv (written by Part A, updated by the CanSWE
    # merge above) rather than assuming Part A's in-memory objects
    # (stations_to_use etc.) still exist - RUN_PART_A may have been FALSE this
    # run, but daily.csv already carries exactly the station_id/lon/lat needed.
    if (is.null(station_lookup)) {
      if (file.exists(DAILY_CSV_PATH)) {
        station_lookup <- readr::read_csv(DAILY_CSV_PATH, show_col_types = FALSE) %>%
          dplyr::distinct(station_id, lon, lat)
      } else {
        warning(
          "No station_lookup supplied and 'daily.csv' does not exist - C-SNOW ",
          "will still be summarized basin-wide below, but no station-level ",
          "extraction will be produced. Run Part A (or pass station_lookup=) ",
          "to get station-level C-SNOW values for validating it against the ",
          "ground-truth SWE/depth record you already have."
        )
        station_lookup <- data.frame(station_id = character(0), lon = numeric(0), lat = numeric(0))
      }
    }
    
    n_read_failed <- 0
    n_processed_this_run <- 0
    
    if (n_files_todo == 0) {
      message("  [OK] Every file in the study window was already processed in a previous run - nothing new to do this run.")
    } else {
      
      # ---- Confirm the data variable name from the FIRST file only ----
      # Assumed consistent across the whole archive (same product, presumably
      # produced by the same export pipeline for every file).
      varname <- detect_csnow_varname(nc_files[1])
      if (is.na(varname)) {
        stop(
          "Could not identify a data variable inside '", basename(nc_files[1]), "' - ",
          "inspect it directly (ncdf4::nc_open(path) -> names(nc$var)) and adjust ",
          "detect_csnow_varname()'s pattern if this auto-detection guessed wrong."
        )
      }
      message("  [OK] Using NetCDF variable '", varname, "' as C-SNOW snow depth.")
      
      station_pts <- NULL
      if (nrow(station_lookup) > 0) {
        station_pts <- sf::st_as_sf(station_lookup, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
      }
      
      # ------------------------------------------------------------------------
      # PERFORMANCE: probe the first (remaining) file's CRS/extent once
      # ------------------------------------------------------------------------
      # This local C-SNOW archive is one tile/product exported the same way for
      # every day (same CRS, same extent) - re-deriving the reprojected basin
      # polygon AND the reprojected station points from scratch on every single
      # file, and running crop()/mask()/global() as separate full-raster
      # passes, added avoidable per-file overhead. This probes file 1 up front
      # so the basin/station geometry can be reprojected ONCE and reused, and
      # so out-of-tile stations can be dropped before the loop even starts (see
      # below). If a later file's CRS is ever found to differ mid-archive, the
      # cache is rebuilt automatically for that file inside the loop.
      # UPDATE: the plain-read-first/subds-fallback ordering below used to be
      # credited as "the fix" for this file's read failure, but inspecting a
      # sample file (snd_20170807.nc) directly showed that was the wrong
      # diagnosis - lat/lon here are proper 1-D dimension-scale variables on
      # a genuinely regular grid, just with float-precision jitter in the
      # stored coordinate values that trips terra's regularity check
      # regardless of plain vs subds. Both attempts return NULL either way.
      # read_csnow_raster() (defined above, next to detect_csnow_varname())
      # keeps this pair as its fast path, then falls back to a manual
      # ncdf4-based read that bypasses the regularity check entirely once
      # it's confirmed the jitter is small relative to cellsize. See that
      # function's DIAGNOSIS comment for the full explanation.
      probe_r <- read_csnow_raster(nc_files[1], varname)
      if (is.null(probe_r)) {
        stop(
          "Could not open the first C-SNOW file ('", basename(nc_files[1]), "') to ",
          "establish its CRS/extent - cannot proceed with C-SNOW ingestion."
        )
      }
      probe_crs <- terra::crs(probe_r)
      probe_ext <- terra::ext(probe_r)
      
      # ------------------------------------------------------------------------
      # PERFORMANCE (correctness-neutral): drop stations outside the tile
      # ------------------------------------------------------------------------
      # terra::extract() at a point outside a raster's extent can only ever
      # return NA - it is not possible for such a station to contribute a real
      # value. daily.csv (the source of station_lookup) now includes every
      # CanSWE station nationwide, not just ones near the Nechako basin (see
      # the CanSWE merge above), so most stations here can sit far outside
      # this C-SNOW tile. This only removes stations that were guaranteed
      # all-NA anyway; it never discards a station that could have produced a
      # real value, and it has no effect on the basin-wide summary below.
      if (!is.null(station_pts) && nrow(station_pts) > 0) {
        pts_native <- tryCatch(sf::st_transform(station_pts, probe_crs), error = function(e) NULL)
        if (!is.null(pts_native)) {
          xy <- sf::st_coordinates(pts_native)
          in_extent <- xy[, "X"] >= probe_ext$xmin & xy[, "X"] <= probe_ext$xmax &
            xy[, "Y"] >= probe_ext$ymin & xy[, "Y"] <= probe_ext$ymax
          n_dropped <- sum(!in_extent, na.rm = TRUE)
          if (n_dropped > 0) {
            message(sprintf(
              "  [OK] Excluding %d/%d station(s) that fall outside the C-SNOW tile's extent from station-level extraction (guaranteed NA there either way; the basin-wide summary below is unaffected).",
              n_dropped, nrow(station_lookup)
            ))
          }
          station_lookup <- station_lookup[in_extent, , drop = FALSE]
          station_pts    <- station_pts[in_extent, , drop = FALSE]
        }
      }
      
      # ------------------------------------------------------------------------
      # PERFORMANCE: cache the reprojected basin polygon + station points
      # ------------------------------------------------------------------------
      cached_crs <- probe_crs
      basin_v_cache <- terra::vect(sf::st_transform(basin_wgs84, probe_crs))
      station_v_cache <- if (!is.null(station_pts) && nrow(station_pts) > 0) {
        terra::vect(sf::st_transform(station_pts, probe_crs))
      } else {
        NULL
      }
      
      refresh_cache_if_needed <- function(r_crs) {
        if (!identical(r_crs, cached_crs)) {
          message("  [INFO] CRS differs from the cached one for this file - rebuilding cached basin/station geometry.")
          basin_v_cache <<- terra::vect(sf::st_transform(basin_wgs84, r_crs))
          station_v_cache <<- if (!is.null(station_pts) && nrow(station_pts) > 0) {
            terra::vect(sf::st_transform(station_pts, r_crs))
          } else {
            NULL
          }
          cached_crs <<- r_crs
        }
        invisible(NULL)
      }
      
      # ------------------------------------------------------------------------
      # RESUMABILITY: flush buffered results to disk periodically
      # ------------------------------------------------------------------------
      # Rather than accumulating every result in memory and writing once at
      # the very end (the previous behaviour, and the reason an interrupted
      # run lost everything), results are buffered only since the last flush
      # and appended to the output CSVs every FLUSH_EVERY files - the same
      # cadence as progress reporting. This bounds how much work an unclean
      # interruption (crash, kill, laptop sleep) can lose to at most
      # FLUSH_EVERY files, instead of the entire run.
      station_csv_exists <- file.exists(station_csv_path)
      basin_csv_exists    <- file.exists(basin_csv_path)
      
      flush_results <- function(buf_station, buf_basin) {
        s_df <- dplyr::bind_rows(Filter(Negate(is.null), buf_station))
        b_df <- dplyr::bind_rows(Filter(Negate(is.null), buf_basin))
        if (nrow(s_df) > 0) {
          readr::write_csv(s_df, station_csv_path, append = station_csv_exists, col_names = !station_csv_exists)
          station_csv_exists <<- TRUE
        }
        if (nrow(b_df) > 0) {
          readr::write_csv(b_df, basin_csv_path, append = basin_csv_exists, col_names = !basin_csv_exists)
          basin_csv_exists <<- TRUE
        }
        invisible(NULL)
      }
      
      buf_station_rows <- list()
      buf_basin_rows <- list()
      FLUSH_EVERY <- 50  # more frequent than the old "write once at the end" so progress is actually saved
      
      t_start <- Sys.time()
      
      for (i in seq_len(n_files_todo)) {
        f <- nc_files[i]
        d <- file_dates[i]
        
        # UPDATE: uses the same read_csnow_raster() helper as the probe read
        # above - see that call site's comment and the function's own
        # DIAGNOSIS block for why the old plain/subds-only ordering wasn't
        # actually sufficient (both fail on this archive's float-jitter lat/lon;
        # the fix needed is the ncdf4 manual-build fallback inside the helper).
        r <- read_csnow_raster(f, varname)
        if (is.null(r)) {
          # NOTE: a file that fails to read is deliberately NOT recorded as
          # "done" anywhere (it never reaches basin_rows below) - so it will
          # be retried automatically on the next run, in case the failure was
          # transient (e.g. a file still being synced) rather than permanent.
          n_read_failed <- n_read_failed + 1
          next
        }
        if (terra::nlyr(r) > 1) r <- r[[1]]
        # Pin the layer to a known name (same reasoning as the extract(dem, ...)
        # note above for the DEM: extract() names its output column after the
        # layer's own name, which - unlike positional indexing - stays correct
        # regardless of what the source NetCDF happened to call its variable).
        names(r) <- "csnow_snd_m"
        
        refresh_cache_if_needed(terra::crs(r))
        
        # Station-level extraction: uses the CACHED, already-reprojected point
        # vector built once above, instead of re-transforming station_pts and
        # rebuilding a SpatVector on every iteration.
        if (!is.null(station_v_cache)) {
          vals <- tryCatch(
            terra::extract(r, station_v_cache)$csnow_snd_m,
            error = function(e) rep(NA_real_, nrow(station_lookup))
          )
          buf_station_rows[[length(buf_station_rows) + 1]] <- data.frame(
            station_id = station_lookup$station_id,
            lon = station_lookup$lon,
            lat = station_lookup$lat,
            date = d,
            csnow_snd_m = vals
          )
        }
        
        # Basin-wide daily summary: mean value + fraction of basin pixels that
        # are actually valid (not masked-out) on this date.
        # PERFORMANCE: crop+mask are done in a single pass via terra::crop()'s
        # own mask= argument (terra >= 1.6) instead of a separate crop() then
        # mask() call, and the two previous terra::global() calls - one of
        # which built an entire extra is.na() raster layer just to count NA
        # cells - are replaced by a single terra::values() pull, with the
        # mean and the valid-pixel fraction both computed from that one
        # in-memory vector in base R.
        # NOTE: tryCatch()'s expr evaluates in THIS frame (not a fresh one), so
        # <<- inside it would skip right past the local basin_mean_val/frac_valid
        # and assign into the enclosing environment instead - a real bug, not
        # just a style nit. Returning a list from expr/the error handler and
        # unpacking it afterward avoids that entirely.
        basin_stats <- tryCatch({
          cropped <- terra::crop(r, basin_v_cache, mask = TRUE)
          v <- terra::values(cropped, mat = FALSE)
          list(
            mean_val   = if (length(v) > 0) mean(v, na.rm = TRUE) else NA_real_,
            frac_valid = if (length(v) > 0) mean(!is.na(v)) else NA_real_
          )
        }, error = function(e) list(mean_val = NA_real_, frac_valid = NA_real_))
        
        buf_basin_rows[[length(buf_basin_rows) + 1]] <- data.frame(
          date = d,
          csnow_basin_mean_m = basin_stats$mean_val,
          frac_basin_valid = basin_stats$frac_valid
        )
        
        n_processed_this_run <- n_processed_this_run + 1
        
        # Flush + report progress together: every FLUSH_EVERY files, and always
        # on the very last file so nothing is left unflushed at the end of a
        # normal (non-interrupted) run.
        if (i %% FLUSH_EVERY == 0 || i == n_files_todo) {
          flush_results(buf_station_rows, buf_basin_rows)
          buf_station_rows <- list()
          buf_basin_rows <- list()
          
          elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
          rate <- if (elapsed_sec > 0) i / elapsed_sec else NA_real_
          eta_sec <- if (!is.na(rate) && rate > 0) (n_files_todo - i) / rate else NA_real_
          message(sprintf(
            "  ...processed %d/%d C-SNOW file(s) this run (%.2f files/sec, elapsed %.0fs, ETA ~%s) - progress saved to disk.",
            i, n_files_todo,
            if (is.na(rate)) 0 else rate,
            elapsed_sec,
            if (is.na(eta_sec)) "unknown" else sprintf("%.0fs", eta_sec)
          ))
        }
      }
      
      if (n_read_failed > 0) {
        message(sprintf(
          "  [WARN] %d/%d file(s) this run could not be read/opened - skipped (corrupt or truncated download?). These dates were NOT recorded as done and will be retried automatically on the next run.",
          n_read_failed, n_files_todo
        ))
      }
    }
    
    # ------------------------------------------------------------------------
    # Read back the FULL accumulated output (this run + all prior runs)
    # ------------------------------------------------------------------------
    # Since results are now written incrementally as the loop runs (or were
    # already fully written by a prior run, if nothing was left to do this
    # time), the CSVs on disk - not the in-memory buffers from this run alone -
    # are the authoritative full dataset. Reading them back keeps this
    # function's return value and summary messages accurate regardless of
    # how many runs it took to get here.
    station_out <- if (file.exists(station_csv_path)) {
      tryCatch(
        readr::read_csv(station_csv_path, show_col_types = FALSE) %>% dplyr::mutate(date = as.Date(date)),
        error = function(e) data.frame()
      )
    } else {
      data.frame()
    }
    basin_out <- if (file.exists(basin_csv_path)) {
      tryCatch(
        readr::read_csv(basin_csv_path, show_col_types = FALSE) %>%
          dplyr::mutate(date = as.Date(date)) %>% dplyr::arrange(date),
        error = function(e) data.frame()
      )
    } else {
      data.frame()
    }
    
    if (nrow(station_out) > 0) {
      n_na <- sum(is.na(station_out$csnow_snd_m))
      message(sprintf(
        "  [OK] '%s' now holds %d row(s) total (%d station(s), %d NA - likely stations outside C-SNOW's valid/unmasked area on that date) across all runs to date.",
        station_csv_path, nrow(station_out), dplyr::n_distinct(station_out$station_id), n_na
      ))
      message(
        "  [NEXT STEP] Join this against daily.csv on station_id/date to see how ",
        "well C-SNOW agrees with your ground-truth SWE/depth record BEFORE relying ",
        "on it further - that comparison is not done automatically by this script."
      )
    }
    if (nrow(basin_out) > 0) {
      mean_cov <- round(100 * mean(basin_out$frac_basin_valid, na.rm = TRUE), 1)
      message(sprintf(
        "  [OK] '%s' now holds %d daily basin-mean row(s) total across all runs to date. Average basin pixel coverage: %.1f%% valid.",
        basin_csv_path, nrow(basin_out), mean_cov
      ))
      message(
        "  [INFO] Low average coverage here is expected in spring/summer (C-SNOW ",
        "masks out wet-snow/snow-free conditions), but if coverage is persistently ",
        "low even in mid-winter, that's a sign the West_US_Canada tile/algorithm ",
        "may not resolve the Nechako's relatively low-relief, forested terrain well - ",
        "worth checking against the station-level comparison above before relying on it."
      )
    }
    
    invisible(list(
      station = station_out, basin = basin_out,
      n_files = n_files,
      n_skipped_already_done = n_skip_already_done,
      n_processed_this_run = n_processed_this_run,
      n_read_failed = n_read_failed
    ))
  }
  
  cmr_download_direct <- function(short_name, version, bbox, start, end,
                                  out_dir, user, pass, page_size = 200,
                                  query_timeout_sec = 60, query_retries = 5,
                                  include_months = 1:12,
                                  n_parallel = 8, file_timeout_sec = 120,
                                  file_retries = 3) {
    
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    
    # ---- HELPER: one CMR metadata request, with a short timeout + retries ----
    # The bare GET() this replaces had no timeout at all, so a stalled
    # connection sat until curl's own low-speed-limit fired ~600s later
    # ("Operation too slow. Less than 1 bytes/sec transferred the last 600
    # seconds") and, since nothing caught it, that error crashed the entire
    # script instead of just this one page request. This wraps each page in
    # httr::timeout() (so a stall gives up in query_timeout_sec, not 600s) and
    # retries with exponential backoff before giving up for good.
    cmr_query_page <- function(page_num) {
      for (attempt in seq_len(query_retries)) {
        resp <- tryCatch(
          GET(
            "https://cmr.earthdata.nasa.gov/search/granules.json",
            query = list(
              short_name   = short_name,
              version      = version,
              bounding_box = paste(bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"], sep = ","),
              temporal     = paste0(start, "T00:00:00Z,", end, "T23:59:59Z"),
              page_size    = page_size,
              page_num     = page_num
            ),
            timeout(query_timeout_sec)
          ),
          error = function(e) {
            cat(sprintf("  CMR query (page %d) attempt %d/%d failed: %s\n",
                        page_num, attempt, query_retries, conditionMessage(e)))
            NULL
          }
        )
        
        if (!is.null(resp)) {
          ok <- tryCatch({ stop_for_status(resp); TRUE },
                         error = function(e) {
                           cat(sprintf("  CMR query (page %d) attempt %d/%d returned HTTP error: %s\n",
                                       page_num, attempt, query_retries, conditionMessage(e)))
                           FALSE
                         })
          if (ok) return(resp)
        }
        
        if (attempt < query_retries) {
          wait_sec <- 10 * attempt  # 10s, 20s, 30s, 40s ... back off between retries
          cat(sprintf("  Retrying page %d in %ds...\n", page_num, wait_sec))
          Sys.sleep(wait_sec)
        }
      }
      stop(sprintf(
        "CMR query for page %d failed after %d attempts (network timeout or repeated HTTP error). ",
        page_num, query_retries),
        "Check https://status.earthdata.nasa.gov for an outage, or re-run the script later - ",
        "already-downloaded granule files are skipped on the next attempt.")
    }
    
    page_num <- 1
    all_granules <- list()
    
    cat("Querying NASA CMR API for granule metadata...\n")
    flush.console()
    
    repeat {
      resp <- cmr_query_page(page_num)
      entries <- content(resp, as = "parsed")$feed$entry
      if (length(entries) == 0) break
      all_granules <- c(all_granules, entries)
      if (length(entries) < page_size) break
      page_num <- page_num + 1
    }
    
    # SEASON_MONTHS filter (applied here, after metadata is fetched but before
    # any actual file download - the metadata query itself is cheap JSON and
    # covers the full [start,end] range regardless, but the granule DOWNLOADS
    # below are the slow part, so filtering the granule list first is what
    # actually saves time).
    if (length(include_months) < 12) {
      granule_month <- vapply(all_granules, function(g) {
        ts <- g$time_start %||% NA_character_
        if (is.na(ts)) return(NA_integer_)
        as.integer(format(as.Date(substr(ts, 1, 10)), "%m"))
      }, integer(1))
      n_before <- length(all_granules)
      all_granules <- all_granules[!is.na(granule_month) & granule_month %in% include_months]
      cat(sprintf(
        "  SEASON_MONTHS filter: kept %d/%d granules in months %s.\n",
        length(all_granules), n_before, paste(month.abb[sort(include_months)], collapse = ", ")
      ))
    }
    
    total_granules <- length(all_granules)
    cat("  ", short_name, " v", version, " - found ", total_granules, " granules total.\n", sep = "")
    flush.console()
    
    MIN_VALID_BYTES <- 50 * 1024
    
    # ---- Build a flat table of every granule's data URL + destination path ----
    granule_rows <- lapply(all_granules, function(g) {
      is_data_link <- sapply(g$links, function(l) isTRUE(grepl("data#$", l$rel)))
      if (!any(is_data_link)) return(NULL)
      data_url  <- g$links[[which(is_data_link)[1]]]$href
      file_name <- basename(data_url)
      data.frame(data_url = data_url, file_name = file_name,
                 dest_file = file.path(out_dir, file_name),
                 stringsAsFactors = FALSE)
    })
    granules_df <- do.call(rbind, granule_rows)
    
    if (is.null(granules_df) || nrow(granules_df) == 0) {
      cat("No downloadable data links found among the granules returned by CMR.\n")
      return(invisible(list(granules = all_granules, total = total_granules, failed = 0L)))
    }
    
    is_valid <- function(f) file.exists(f) && file.info(f)$size >= MIN_VALID_BYTES
    
    # ---- Skip whatever's already on disk and valid (same check as before) ----
    already_ok <- vapply(granules_df$dest_file, is_valid, logical(1))
    cat(sprintf("%d/%d granules already on disk and valid - %d left to download.\n",
                sum(already_ok), nrow(granules_df), sum(!already_ok)))
    flush.console()
    todo <- granules_df[!already_ok, , drop = FALSE]
    
    n_failed <- 0
    n_done <- 0
    
    if (nrow(todo) > 0) {
      cat(sprintf("Downloading %d file(s) using %d parallel connections...\n\n",
                  nrow(todo), n_parallel))
      flush.console()
      
      for (batch_start in seq(1, nrow(todo), by = n_parallel)) {
        batch <- todo[batch_start:min(batch_start + n_parallel - 1, nrow(todo)), , drop = FALSE]
        
        for (attempt in seq_len(file_retries)) {
          reauthenticate()
          
          res <- tryCatch(
            curl::multi_download(
              urls      = batch$data_url,
              destfiles = batch$dest_file,
              resume    = TRUE,
              timeout   = file_timeout_sec,
              progress  = FALSE,
              netrc = 1L, followlocation = TRUE
            ),
            error = function(e) {
              cat(sprintf("  Batch download error (attempt %d/%d): %s\n",
                          attempt, file_retries, conditionMessage(e)))
              NULL
            }
          )
          
          success_col <- if (!is.null(res) && "success" %in% names(res)) res$success else rep(NA, nrow(batch))
          status_col  <- if (!is.null(res) && "status_code" %in% names(res)) res$status_code else rep(NA, nrow(batch))
          
          done <- vapply(seq_len(nrow(batch)), function(i) {
            status_ok <- is.na(status_col[i]) || status_col[i] == 200
            ok <- (is.na(success_col[i]) || isTRUE(success_col[i])) && status_ok &&
              is_valid(batch$dest_file[i])
            if (!ok && file.exists(batch$dest_file[i])) unlink(batch$dest_file[i])
            ok
          }, logical(1))
          
          n_done <- n_done + sum(done)
          batch <- batch[!done, , drop = FALSE]
          if (nrow(batch) == 0) break
          if (attempt < file_retries) Sys.sleep(10 * attempt)
        }
        
        if (nrow(batch) > 0) {
          n_failed <- n_failed + nrow(batch)
          cat(sprintf("  %d file(s) still failing after %d attempts:\n", nrow(batch), file_retries))
          cat(paste0("    ", batch$file_name, collapse = "\n"), "\n")
        }
        
        cat(sprintf("  ...%d/%d downloaded so far, %d failed.\n", n_done, nrow(todo), n_failed))
        flush.console()
      }
    }
    
    cat(sprintf("\nDone. %d/%d granules failed validation after retries and were not saved.\n",
                n_failed, total_granules))
    invisible(list(granules = all_granules, total = total_granules, failed = n_failed))
  }
  
  # ----------------------------------------------------------------------------
  # 3) NSIDC/CMR validation products: CMC, GlobSnow, CanSISE
  # ----------------------------------------------------------------------------
  # CMR's granule search silently returns an EMPTY list (not an error) whenever
  # short_name+version doesn't resolve to a real collection - which looks
  # identical in the logs to "this bbox/time window genuinely has no data".
  # Legacy NSIDC/ECS collections (which is what NSIDC-0447/0595/0668 all are)
  # store their version zero-padded in the actual granule metadata - e.g. the
  # public GranuleUR pattern for these older holdings is "SC:NSIDC-0079.004:...",
  # "SC:NSIDC-0531.001:...", i.e. "004"/"001", not the bare "1"/"2" this script
  # was passing. This check queries CMR's collections.json for the short_name
  # BEFORE running the (slow) granule search, so a bad version string stops
  # immediately with the real version list instead of quietly producing "0
  # granules found" that looks like a data-coverage problem.
  verify_cmr_collection <- function(short_name, version, query_timeout_sec = 30) {
    resp <- tryCatch(
      GET("https://cmr.earthdata.nasa.gov/search/collections.json",
          query = list(short_name = short_name), timeout(query_timeout_sec)),
      error = function(e) NULL
    )
    if (is.null(resp) || httr::http_error(resp)) {
      warning(
        "Could not verify collection '", short_name, "' against CMR (network ",
        "issue) - proceeding without verification; a silent 0-granule result ",
        "may indicate a wrong short_name/version rather than no data."
      )
      return(invisible(NULL))
    }
    entries <- tryCatch(httr::content(resp, as = "parsed")$feed$entry, error = function(e) NULL)
    if (is.null(entries) || length(entries) == 0) {
      stop(
        "CMR has NO collection at all matching short_name = '", short_name, "'. ",
        "Check https://cmr.earthdata.nasa.gov/search/collections.json?short_name=",
        short_name, " directly before continuing - the short_name itself may ",
        "be wrong, or this dataset may not be catalogued in CMR at all (e.g. ",
        "NSIDC-0595/GlobSnow is a landing-page-only record - NSIDC states it ",
        "does not archive that data; get it directly from https://www.globsnow.info/swe/ instead).",
        call. = FALSE
      )
    }
    # BUGFIX: this function used to reference `available_versions` here without
    # ever assigning it from `entries` - every call fell through to R's own
    # "object 'available_versions' not found" error, which tryCatch() in the
    # caller (cmr_download_wrapper) then reported as if it were a genuine CMR
    # failure. That masked the real check below entirely: NSIDC-0447 and
    # NSIDC-0668 were never actually tested for a version mismatch, they just
    # errored on this undefined variable every time. `entries` is the parsed
    # Atom `feed$entry` list from CMR's collections.json; each entry's version
    # is exposed under `version_id` in that schema.
    available_versions <- vapply(
      entries,
      function(e) if (!is.null(e$version_id)) as.character(e$version_id) else NA_character_,
      character(1)
    )
    available_versions <- available_versions[!is.na(available_versions)]
    if (length(available_versions) == 0) {
      stop(
        "CMR collection '", short_name, "' was found, but no entry exposed a ",
        "'version_id' field - the collections.json schema may have changed. ",
        "Inspect https://cmr.earthdata.nasa.gov/search/collections.json?short_name=",
        short_name, " manually to find the correct field name.",
        call. = FALSE
      )
    }
    # Compare after stripping leading zeros, not with an exact string match:
    # the whole point of this check is that legacy NSIDC/ECS collections are
    # INCONSISTENT about zero-padding ("004"/"001" for some collections,
    # bare "1" for others per the collections.json response itself, as
    # confirmed for NSIDC-0447 in practice) - an exact match on our guessed
    # padding therefore stops on a false negative just as easily as it
    # catches a real mismatch. Normalizing and returning whichever actual
    # version string CMR has on file means the caller no longer has to
    # guess the padding correctly to get past this check.
    norm_version <- function(v) sub("^0+(?=[0-9])", "", v, perl = TRUE)
    match_idx <- which(norm_version(available_versions) == norm_version(version))
    if (length(match_idx) == 0) {
      stop(
        "CMR collection '", short_name, "' exists, but not with version = '",
        version, "'. Version(s) actually available: ",
        paste(unique(available_versions), collapse = ", "), ". Legacy NSIDC/ECS ",
        "collections commonly store this zero-padded (e.g. '001' rather than ",
        "'1') - try one of the versions listed above.",
        call. = FALSE
      )
    }
    invisible(available_versions[match_idx[1]])
  }
  
  cmr_download_wrapper <- function(short_name, version, out_dir, start = study_start, end = study_end) {
    message(sprintf("\n===== NSIDC/CMR acquisition: %s v%s =====", short_name, version))
    ensure_dir(out_dir)
    # Use whatever version string CMR actually has on file (verify_cmr_collection
    # now normalizes zero-padding and returns the real one) rather than the
    # possibly mis-padded guess passed in, since cmr_download_direct's granule
    # query below needs the exact version CMR recognizes.
    resolved_version <- verify_cmr_collection(short_name, version)
    if (!is.null(resolved_version)) version <- resolved_version
    write_earthdata_netrc(Sys.getenv("EARTHDATA_USER", unset = ""), Sys.getenv("EARTHDATA_PASS", unset = ""))
    cmr_download_direct(
      short_name = short_name,
      version = version,
      bbox = bbox_vec,
      start = format(start, "%Y-%m-%d"),
      end = format(end, "%Y-%m-%d"),
      out_dir = out_dir,
      user = Sys.getenv("EARTHDATA_USER", unset = ""),
      pass = Sys.getenv("EARTHDATA_PASS", unset = ""),
      include_months = 1:12
    )
  }
  
  # ----------------------------------------------------------------------------
  # 4) OpenLandMap snow-cover climatology (interannual P05/P95 envelope)
  # ----------------------------------------------------------------------------
  # METHODOLOGY NOTE (revised): OpenLandMap/Hengl (2021), record 6011200,
  # distributes THREE quantile products per calendar year-month:
  #   clm_snow.cover_esa.modis_p05_..._YYYY.MM_...tif
  #   clm_snow.cover_esa.modis_p50_..._YYYY.MM_...tif
  #   clm_snow.cover_esa.modis_p95_..._YYYY.MM_...tif
  # Per the Zenodo record description, these p05/p50/p95 values are quantiles
  # of the daily (or weekly) ESA CCI snow-cover-fraction values WITHIN that
  # single year-month - i.e. each file already describes within-month
  # variability for one specific year, not variability across years.
  #
  # The previous version of this function fetched only the p05 and p95 files
  # for each year and averaged them across 2000-2019 with terra::app(mean),
  # producing "the average within-month 5th percentile" and "the average
  # within-month 95th percentile" - a central-tendency climatology of a
  # within-month statistic. That is NOT the same object as an interannual
  # P05-P95 envelope, which is what H8SWEI_Snow.R's validate_scf_timing()
  # actually needs: a band that captures how much a calendar month's
  # snow-cover frequency varies from YEAR TO YEAR, against which one
  # particular year's corrected SCF can be judged as normal or anomalous.
  # Averaging pre-computed within-month percentiles across years smooths
  # over exactly the interannual spread that envelope is supposed to capture,
  # and does so unevenly across months (this function historically averaged
  # anywhere from 11 to 19 of the 20 available years per month, per the
  # download log, so the smoothing bias itself varies by month).
  #
  # FIX: use the p50 (within-month median) file for each year as that year's
  # representative monthly SCF value - the closest single-valued analogue to
  # "this month's snow-cover fraction in this particular year" - then compute
  # the TRUE empirical interannual 5th/95th percentile across the 20 yearly
  # p50 values, pixel-by-pixel. This is computed directly from raw per-year
  # data rather than by averaging someone else's already-reduced percentile
  # product, so it is the methodologically appropriate statistic for an
  # interannual validation envelope.
  acquire_openlandmap_climatology <- function(record_id = 6011200,
                                              basin = basin_wgs84,
                                              out_dir = "snow_data/openlandmap_snow_climatology",
                                              years = 2000:2019) {
    message("\n===== OpenLandMap snow-climatology acquisition =====")
    ensure_dir(out_dir)
    
    rec <- tryCatch(jsonlite::fromJSON(sprintf("https://zenodo.org/api/records/%s", record_id), flatten = TRUE),
                    error = function(e) NULL)
    if (is.null(rec) || is.null(rec$files)) {
      warning("OpenLandMap Zenodo record could not be queried; skipping.")
      return(invisible(NULL))
    }
    
    files <- rec$files
    if (!is.data.frame(files)) files <- as.data.frame(files, stringsAsFactors = FALSE)
    file_name_col <- intersect(c("filename", "key", "name"), names(files))
    if (length(file_name_col) == 0) stop("Could not identify filename column in OpenLandMap record.")
    url_col <- intersect(c("links.download", "url", "download"), names(files))
    if (length(url_col) == 0) url_col <- grep("download", names(files), value = TRUE)
    
    if (length(url_col) == 0) {
      files$url <- sprintf(
        "https://zenodo.org/api/records/%s/files/%s/content",
        record_id, utils::URLencode(files[[file_name_col[1]]], reserved = TRUE)
      )
      url_col <- "url"
    }
    
    # Real filenames from this record look like:
    #   clm_snow.cover_esa.modis_p05_1km_s0..0cm_2000.03_epsg4326_v1.tif
    # i.e. a single "_YYYY.MM_" token (digits on both sides of the dot), not a
    # "_MM.mon_" token with a 3-letter month name - so year and month are
    # extracted together from one capture, not two separate patterns.
    fname_col <- files[[file_name_col[1]]]
    ym_matches <- regmatches(fname_col, regexec("_(20\\d{2})\\.(\\d{2})_", fname_col))
    files$year  <- suppressWarnings(as.integer(vapply(ym_matches, function(m) if (length(m) >= 3) m[2] else NA_character_, character(1))))
    files$month <- suppressWarnings(as.integer(vapply(ym_matches, function(m) if (length(m) >= 3) m[3] else NA_character_, character(1))))
    # This record also has a P.90/SD monthly-climatology family already
    # pre-averaged over 2000-2012 (e.g. "..._apr_p.90_..._2000..2012_v2.0.tif")
    # - those don't match "_YYYY.MM_" at all (they use "_2000..2012_" and a
    # 3-letter month elsewhere in the name) and are correctly excluded here.
    # We now pull "p50" (within-month median) rather than "p05"/"p95", since
    # this function computes its own interannual percentiles from the yearly
    # p50 values instead of averaging OpenLandMap's within-month percentiles.
    files$qtag  <- ifelse(grepl("p50", fname_col, ignore.case = TRUE), "p50", NA_character_)
    
    files <- files[!is.na(files$month) & !is.na(files$qtag) & !is.na(files$year), , drop = FALSE]
    if (nrow(files) == 0) {
      warning("No OpenLandMap monthly p50 COGs were recognized in the Zenodo record.")
      return(invisible(NULL))
    }
    
    basin_ext <- terra::ext(terra::vect(sf::st_transform(basin, 4326)))
    mon_lab_all <- tolower(format(as.Date(paste0("2020-", 1:12, "-01")), "%b"))
    
    # Local cache for the cropped per-year-month source tiles (p50), so a
    # re-run (e.g. after an interruption partway through the 12-month loop,
    # or a re-run after this function's own P05/P95 outputs were deleted to
    # rebuild them) reuses whatever's already on disk instead of re-streaming
    # every year x month combination from Zenodo again. Only the FINAL P05/P95
    # rasters were being checked for pre-existence before; the raw p50 inputs
    # that feed them had no local cache at all.
    cache_dir <- file.path(out_dir, "_cache_p50")
    ensure_dir(cache_dir)
    
    # ------------------------------------------------------------------------
    # CIRCUIT BREAKER: skip the wasted retry/backoff once a month is confirmed
    # streaming-broken
    # ------------------------------------------------------------------------
    # Observed behaviour: /vsicurl/ streaming fails for EVERY year within a
    # given month (e.g. all 20 years of "07"), each one burning through 4
    # attempts x (5+10+15)s = 30s of Sys.sleep() before the full-file fallback
    # finally kicks in and succeeds. That's ~30s x 20 years = 10 minutes of
    # pure waiting per broken month, on top of the actual download time - and
    # this record has up to 12 months x 20 years = 240 combinations total.
    # Once STREAK_THRESHOLD consecutive years in the SAME month have hit the
    # streaming failure, this remembers that month as streaming-broken (in
    # streaming_dead_months, keyed by month number) and every subsequent year
    # in that month skips straight to the full-file download - no more
    # 4-attempt/backoff cycle wasted on a path that's already been shown not
    # to work for this month. A NEW month starts with a clean slate (a broken
    # July says nothing about August), and if a month's very first year
    # succeeds via streaming, its streak is reset to 0 rather than latching a
    # single early failure permanently.
    streaming_dead_months <- new.env(parent = emptyenv())
    streaming_fail_streak <- new.env(parent = emptyenv())
    STREAK_THRESHOLD <- 2  # consecutive streaming failures within one month before skipping the retry cycle for its remaining years
    
    fetch_cropped_month <- function(qtag, year, month, max_retries = 4) {
      cache_file <- file.path(cache_dir, sprintf("%s_%d_%02d.tif", qtag, year, month))
      if (file.exists(cache_file)) {
        r <- tryCatch(terra::rast(cache_file), error = function(e) NULL)
        if (!is.null(r)) return(r)
        # Cached file exists but failed to read (e.g. truncated from a prior
        # interrupted run) - fall through and re-fetch/re-cache it below.
        unlink(cache_file)
      }
      row <- files[files$qtag == qtag & files$year == year & files$month == month, , drop = FALSE]
      if (nrow(row) == 0) return(NULL)
      url <- row[[url_col[1]]][1]
      
      month_key <- sprintf("%s_%02d", qtag, month)
      skip_streaming <- isTRUE(get0(month_key, envir = streaming_dead_months, ifnotfound = FALSE))
      
      # Probe the URL before asking GDAL/terra to parse it. This distinguishes
      # a permanent missing-file condition (404) from a transient server-side
      # failure, and avoids retrying cases that can never succeed.
      probe <- tryCatch(
        httr::HEAD(url, timeout(20)),
        error = function(e) NULL
      )
      if (!is.null(probe)) {
        sc <- httr::status_code(probe)
        if (sc == 404) {
          message(sprintf(
            "    [MISS] %s %d.%02d - Zenodo returned HTTP 404 for the file URL; this file is absent under the expected name, so skipping retries.",
            qtag, year, month
          ))
          return(NULL)
        }
        if (sc >= 500) {
          message(sprintf(
            "    [retry-note] %s %d.%02d - Zenodo returned HTTP %d before parsing; proceeding with normal open/retry logic.",
            qtag, year, month, sc
          ))
        }
      }
      
      # Retry with backoff: this used to be a single bare terra::rast(/vsicurl/...)
      # attempt, so a transient Zenodo-side hiccup (rate limiting, a brief
      # outage) was indistinguishable from "this file doesn't exist" - the
      # whole month got silently abandoned (see the "No OpenLandMap p50
      # rasters found for month NN" warning below) even though every other
      # download path in this script (CMR, ASF, CDSE) already retries
      # transient failures before giving up. A 100%-failure-rate across every
      # year for a given month, starting at the very first live request of a
      # run, is the signature of a transient Zenodo-side issue, not 20
      # simultaneously-broken files - mirrors that same retry/backoff pattern.
      #
      # skip_streaming (set above from streaming_dead_months) short-circuits
      # this loop entirely once STREAK_THRESHOLD consecutive years in THIS
      # month have already exhausted it and needed the full-file fallback -
      # at that point retrying streaming again for the next year is not
      # gathering new information, just re-paying the same ~30s cost this
      # month has already proven it will not avoid.
      r <- NULL
      if (skip_streaming) {
        message(sprintf(
          "    [skip-streaming] %s %d.%02d - month %02d was already confirmed streaming-broken (%d+ consecutive prior years needed the full-file fallback) - going straight to full-file download.",
          qtag, year, month, month, STREAK_THRESHOLD
        ))
      } else {
        for (attempt in seq_len(max_retries)) {
          r <- tryCatch(terra::rast(paste0("/vsicurl/", url)), error = function(e) NULL)
          if (!is.null(r)) break
          if (attempt < max_retries) {
            wait_sec <- 5 * attempt
            message(sprintf(
              "    [retry] %s %d.%02d attempt %d/%d failed to open - retrying in %ds...",
              qtag, year, month, attempt, max_retries, wait_sec
            ))
            Sys.sleep(wait_sec)
          } else {
            if (!is.null(probe) && httr::status_code(probe) >= 500) {
              message(sprintf(
                "    [retry-note] %s %d.%02d - /vsicurl/ streaming failed after %d attempts (Zenodo returned HTTP %d; transient server-side failure). Falling back to a full-file download before giving up.",
                qtag, year, month, max_retries, httr::status_code(probe)
              ))
            } else {
              message(sprintf(
                "    [retry-note] %s %d.%02d - /vsicurl/ streaming failed after %d attempts (streamed response was not a readable GeoTIFF/COG - often a range-request/CDN quirk, not a missing file). Falling back to a full-file download before giving up.",
                qtag, year, month, max_retries
              ))
            }
          }
        }
        
        # Update the circuit breaker based on whether STREAMING (not the
        # eventual full-file fallback) succeeded this year, so a month that
        # keeps failing streaming gets flagged, while a month where streaming
        # is intermittently fine again has its streak cleared.
        prior_streak <- get0(month_key, envir = streaming_fail_streak, ifnotfound = 0L)
        if (is.null(r)) {
          new_streak <- prior_streak + 1L
          assign(month_key, new_streak, envir = streaming_fail_streak)
          if (new_streak >= STREAK_THRESHOLD) {
            assign(month_key, TRUE, envir = streaming_dead_months)
          }
        } else {
          assign(month_key, 0L, envir = streaming_fail_streak)
        }
      }
      
      # FALLBACK: full-file download instead of /vsicurl/ range-request
      # streaming. The "[FAIL] ... not a readable GeoTIFF/COG" pattern is
      # consistent with Zenodo's CDN occasionally mishandling the partial-
      # content (HTTP Range) requests GDAL issues under /vsicurl/, rather
      # than the object itself being gone (a genuine 404 was already caught
      # and returned earlier, before any of this runs). Pulling the whole
      # file with a plain HTTP GET sidesteps range-request handling
      # entirely, so it recovers exactly the cases where streaming fails for
      # that reason. NOTE: this is still the same underlying Zenodo object
      # via the same URL, not an independent second source - Zenodo is the
      # sole distribution point for this record - so it fixes transient/CDN-
      # shaped failures but will not help with a genuinely missing/renamed
      # file (those already return NULL above via the 404 check).
      if (is.null(r)) {
        tmp_full <- file.path(cache_dir, sprintf("_dlfull_%s_%d_%02d.tif", qtag, year, month))
        dl_ok <- tryCatch({
          download_file_auth(url, tmp_full)
          file.exists(tmp_full) && file.info(tmp_full)$size > 0
        }, error = function(e) FALSE)
        if (isTRUE(dl_ok)) {
          r <- tryCatch(terra::rast(tmp_full), error = function(e) NULL)
          if (!is.null(r)) {
            message(sprintf(
              "    [OK] %s %d.%02d - recovered via full-file download fallback (/vsicurl/ streaming had failed).",
              qtag, year, month
            ))
          } else {
            unlink(tmp_full)
          }
        }
      }
      
      if (is.null(r)) {
        message(sprintf(
          "    [FAIL] %s %d.%02d - could not open after %d streaming attempts AND a full-file download fallback; skipping this year/month (check for a renamed/moved object or an HTML/error body at the URL).",
          qtag, year, month, max_retries
        ))
        return(NULL)
      }
      
      r <- tryCatch(terra::crop(r, basin_ext), error = function(e) r)
      if (isTRUE("SpatRaster" %in% class(r))) {
        if (!is.null(basin)) {
          try(r <- terra::mask(r, terra::vect(sf::st_transform(basin, 4326))), silent = TRUE)
        }
      }
      try(terra::writeRaster(r, cache_file, overwrite = TRUE), silent = TRUE)
      # If the full-file-download fallback was used above, its raw temp copy
      # is now redundant once cache_file has been written (terra::rast() on
      # a file path is normally lazy/memory-mapped, so this cleanup is
      # deferred until after the writeRaster() call just above, not placed
      # right after the download itself) - remove it so repeated fallback
      # hits across a full run (up to 240 year-month combinations) don't
      # leave duplicate raw copies sitting in cache_dir alongside the
      # already-cropped/masked cache_file.
      if (exists("tmp_full") && file.exists(tmp_full) && file.exists(cache_file)) {
        unlink(tmp_full)
      }
      # Small pause after every LIVE fetch (cache hits return early above and
      # never reach here) so a full month's worth of requests (up to 20 in a
      # row, one per year) doesn't hammer Zenodo's API in a tight loop - a
      # likely contributor to the rate-limiting/outage pattern this retry
      # logic is otherwise just working around after the fact.
      Sys.sleep(0.5)
      r
    }
    
    cat(sprintf("  [Section 4] Building interannual P05-P95 SCF envelope from yearly p50 (%d-%d)\n",
                min(years), max(years)))
    # Tracks months where fewer years came back than were requested, so the
    # caller (and ultimately the final inventory) can report OpenLandMap as
    # genuinely incomplete rather than a blanket "OK" - previously a
    # month-level shortfall was only ever visible as a warning() buried in
    # the console log, and the caller (acquire_openlandmap_climatology()'s
    # own return value) had no way to know it happened at all.
    incomplete_months <- character(0)
    for (m in 1:12) {
      out_file_p05 <- file.path(out_dir, sprintf("openlandmap_scf_p05_%s.tif", mon_lab_all[m]))
      out_file_p95 <- file.path(out_dir, sprintf("openlandmap_scf_p95_%s.tif", mon_lab_all[m]))
      if (file.exists(out_file_p05) && file.exists(out_file_p95)) {
        message("    [OK] ", basename(out_file_p05), " / ", basename(out_file_p95), " already exist - skipping.")
        next
      }
      year_rasters <- list()
      for (yr in years) {
        r <- tryCatch(fetch_cropped_month("p50", yr, m), error = function(e) NULL)
        if (!is.null(r)) year_rasters[[length(year_rasters) + 1]] <- r
      }
      if (length(year_rasters) < length(years)) {
        incomplete_months <- c(incomplete_months, sprintf(
          "%s (%d/%d years)", mon_lab_all[m], length(year_rasters), length(years)
        ))
      }
      if (length(year_rasters) == 0) {
        warning(sprintf("No OpenLandMap p50 rasters found for month %02d.", m))
        next
      }
      template <- year_rasters[[1]]
      aligned <- lapply(year_rasters, function(r) {
        if (terra::compareGeom(r, template, stopOnError = FALSE)) r else terra::resample(r, template, method = "near")
      })
      stacked <- terra::rast(aligned)
      # True empirical interannual percentiles, computed pixel-by-pixel across
      # the yearly p50 (within-month median) values - this is the actual
      # year-to-year distribution for this calendar month, not an average of
      # a within-month statistic.
      p05_month <- terra::app(stacked, fun = function(x) stats::quantile(x, probs = 0.05, na.rm = TRUE, names = FALSE))
      p95_month <- terra::app(stacked, fun = function(x) stats::quantile(x, probs = 0.95, na.rm = TRUE, names = FALSE))
      names(p05_month) <- sprintf("p05_%s_interannual", mon_lab_all[m])
      names(p95_month) <- sprintf("p95_%s_interannual", mon_lab_all[m])
      terra::writeRaster(p05_month, out_file_p05, overwrite = TRUE)
      terra::writeRaster(p95_month, out_file_p95, overwrite = TRUE)
      message(sprintf("    [OK] Wrote %s and %s (%d year(s) of p50 used to derive interannual P05/P95).",
                      basename(out_file_p05), basename(out_file_p95), length(year_rasters)))
    }
    if (length(incomplete_months) > 0) {
      message(sprintf(
        "  [WARN] %d/12 month(s) have fewer years than requested in their interannual P05/P95 (retries were exhausted for the missing years - see [FAIL] messages above for which ones and why): %s",
        length(incomplete_months), paste(incomplete_months, collapse = ", ")
      ))
    }
    message(sprintf("  [OK] OpenLandMap interannual climatology envelope written to %s", out_dir))
    invisible(list(out_dir = out_dir, incomplete_months = incomplete_months))
  }
  
  # ----------------------------------------------------------------------------
  # GlobSnow (FMI) acquisition - SWE (monthly, bias-corrected) + SE (snow extent)
  # ----------------------------------------------------------------------------
  # GlobSnow is distributed directly by FMI (globsnow.info), not through
  # NSIDC/CMR - NSIDC-0595 is a landing-page-only record ("NSIDC does not
  # archive these data"), so cmr_download_wrapper()/CMR can never return
  # anything for it. FMI's own directory listings ARE plain HTTPS autoindex
  # pages (no Earthdata Login / cookie handshake needed, unlike NSIDC's
  # daacdata.apps - see acquire_nsidc_https_dir() above for that case), so
  # this lists them directly rather than guessing per-year/per-month filenames.
  #
  # Two distinct GlobSnow products are acquired here, matching the two file
  # H8SWEI_Snow.R actually looks for:
  #   1) SWE - already downloaded locally (per the local snow_data/globsnow/swe
  #      folder: "YYYYMM_northern_hemisphere_monthly_biascorrected_swe_0.25grid.nc",
  #      one file per year-month, GlobSnow v3.0 monthly bias-corrected product).
  #      This only fills in whatever months are MISSING locally - it does not
  #      re-download months you already have.
  #   2) SE (snow extent / fractional snow cover) - this is the "which data
  #      should I get" gap: H8 wants a single combined-time NetCDF stack at
  #      monthly_data_direct/snow_cover_monthly.nc (see scf_file further up in
  #      H8SWEI_Snow.R), which needs to be the MONTHLY-AGGREGATED product
  #      (MFSC = Monthly Aggregated Fractional Snow Cover), NOT:
  #        - D4SC/D4SC_ql: daily 4-CLASS (categorical, not fractional) browse
  #          images (the .jpg/.png quicklooks in the archive_v2.1/<year>/D4SC_ql/
  #          folder shown in your screenshot - these are for visual browsing
  #          only, not gridded data for analysis)
  #        - DFSC/DFSC_ql: DAILY fractional snow cover - wrong temporal
  #          resolution for a "monthly" file
  #        - WFSC: WEEKLY aggregated fractional snow cover - also wrong
  #          temporal resolution
  #      MFSC is the one true monthly, NetCDF (not image), fractional (not
  #      4-class) product - the direct analogue of the SWE product's monthly
  #      bias-corrected file, and the only one of the five that H8's monthly
  #      SCF stack actually calls for.
  #      GlobSnow SE archive layout confirmed from the GlobSnow-2 Product User
  #      Guide: one catalogue per product version (e.g. archive_v2.1/), one
  #      sub-catalogue per year, containing DFSC / D4SC / WFSC / MFSC (+ their
  #      _ql quicklook siblings). Filenames follow
  #      "GlobSnow_SE_MFSC_L3B_NH_YYYYMM_v2.1_1.nc" (confirmed against the
  #      D4SC_ql daily naming pattern visible in your screenshot,
  #      "GlobSnow_SE_4CL_L3A_NH_YYYYMMDD_v2.1_1.jpg/.png" - MFSC follows the
  #      same "_v2.1_1" suffix convention, just monthly and NetCDF instead of
  #      daily JPG/PNG).
  
  GLOBSNOW_SE_ARCHIVE_VERSION <- "v2.1"
  # BUGFIX: the previous URL - ".../swe/archive_v3.0/L3B_monthly_biascorrected/NetCDF4/" -
  # 404ed on every attempt (see the log: 4/4 retries returned HTTP 404). The
  # actual directory name is "L3B_monthly_biascorrected_SWE" (note: "_SWE"
  # suffix, not a separate "/NetCDF4/" subfolder underneath a differently-named
  # parent) - confirmed against the exact URL cited independently in three
  # published sources: the GlobSnow v3.0 data-descriptor paper (Luojus et al.
  # 2021, Scientific Data, 10.1038/s41597-021-00939-2), the Nature paper this
  # dataset was produced for (Pulliainen et al. 2020, 10.1038/s41586-020-2258-0),
  # and a 2022 JGR-Atmospheres paper's data-availability statement - all three
  # give "https://www.globsnow.info/swe/archive_v3.0/L3B_monthly_biascorrected_SWE/"
  # verbatim. NetCDF and HDF files for a given product live together in that
  # same folder (per the GlobSnow SWE product guide's own catalogue listing),
  # not sorted into a format-specific subfolder - filter by extension in the
  # directory listing instead of assuming a NetCDF4/ path segment exists.
  GLOBSNOW_SWE_BASE_URL <- "https://www.globsnow.info/swe/archive_v3.0/L3B_monthly_biascorrected_SWE/"
  GLOBSNOW_SE_BASE_URL  <- sprintf("https://www.globsnow.info/se/archive_%s/", GLOBSNOW_SE_ARCHIVE_VERSION)
  
  # Lists an HTTPS autoindex directory (plain Apache/nginx-style listing, no
  # auth) and returns the href values. Shares the same "is this actually a
  # directory listing" caution as acquire_nsidc_https_dir(), minus the
  # cookie-jar/URS handshake machinery that's specific to NSIDC's
  # daacdata.apps login-gated system - globsnow.info serves these directly.
  list_https_dir <- function(dir_url, query_timeout_sec = 60, query_retries = 4) {
    resp <- NULL
    for (attempt in seq_len(query_retries)) {
      resp <- tryCatch(
        httr::GET(dir_url, timeout(query_timeout_sec)),
        error = function(e) {
          message(sprintf("  Directory listing attempt %d/%d failed for %s: %s",
                          attempt, query_retries, dir_url, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(resp) && !httr::http_error(resp)) break
      
      # BUGFIX: a 404 means "this path does not exist" - it is a PERMANENT
      # condition, not a transient server hiccup, so retrying it query_retries
      # times with exponential backoff (10s, 20s, 30s = ~60s total, as seen in
      # the log for the dead GlobSnow SWE URL) only delays reporting a failure
      # that was already certain on attempt 1. Every other status code (5xx,
      # timeouts, connection resets) genuinely can be transient and keeps
      # retrying as before; only 404 short-circuits.
      if (!is.null(resp) && httr::status_code(resp) == 404) {
        message(sprintf(
          "  [SKIP] Directory listing for %s returned HTTP 404 (path does not exist) - not retrying, since a 404 will not change on its own. If this URL used to work, the site may have reorganized its layout; verify the current path in a browser before re-running.",
          dir_url
        ))
        return(character(0))
      }
      
      if (!is.null(resp)) {
        message(sprintf(
          "  Directory listing attempt %d/%d returned HTTP %d for %s - %s",
          attempt, query_retries, httr::status_code(resp), dir_url,
          if (attempt < query_retries) "retrying..." else "giving up."
        ))
      }
      resp <- NULL
      if (attempt < query_retries) Sys.sleep(10 * attempt)
    }
    if (is.null(resp)) {
      message(sprintf(
        "  [SKIP] Could not list directory %s after %d attempt(s) - see HTTP status/error messages above.",
        dir_url, query_retries
      ))
      return(character(0))
    }
    
    html <- httr::content(resp, as = "text", encoding = "UTF-8")
    looks_like_listing <- grepl("Index of|Parent Directory|<title>Index", html, ignore.case = TRUE)
    if (!looks_like_listing) {
      message(sprintf(
        "  [SKIP] URL %s returned HTTP %d but does not look like a directory listing - treating as empty.",
        dir_url, httr::status_code(resp)
      ))
      warning("URL ", dir_url, " returned HTTP ", httr::status_code(resp),
              " but does not look like a directory listing - treating as empty.",
              call. = FALSE)
      return(character(0))
    }
    hrefs <- unique(unlist(regmatches(html, gregexpr('href="[^"]+"', html, perl = TRUE))))
    hrefs <- sub('^href="', '', hrefs)
    hrefs <- sub('"$', '', hrefs)
    hrefs[!grepl("^\\.\\.?/?$|^\\?|^https?://", hrefs, perl = TRUE)]
  }
  
  # ------------------------------------------------------------------
  # Shared helper: given a remote filename, generate all plausible local
  # cache paths that would count as "already have this file" (covers .nc,
  # .zip, .gz, and combinations like .nc.zip).
  # ------------------------------------------------------------------
  # NOTE: this used to be defined ONLY as a nested/local function inside
  # acquire_globsnow_swe() - fine for that function alone, but
  # acquire_globsnow_se_mfsc() is a separate sibling function (not nested
  # inside acquire_globsnow_swe()) and could not see it, which would have
  # thrown "could not find function 'generate_cache_candidates'" the moment
  # acquire_globsnow_se_mfsc() tried to use it. Hoisted here to the shared
  # Part B helper scope (alongside ensure_dir(), download_file_auth(), etc.)
  # so both GlobSnow acquisition functions can call it.
  generate_cache_candidates <- function(fn, cache_dir) {
    base <- fn
    base <- sub("\\.gz$", "", base, ignore.case = TRUE)
    base <- sub("\\.(nc|zip)$", "", base, ignore.case = TRUE)
    base <- sub("\\.nc$", "", base, ignore.case = TRUE)
    
    candidates <- c(
      file.path(cache_dir, fn),
      file.path(cache_dir, paste0(base, ".nc")),
      file.path(cache_dir, paste0(base, ".zip")),
      file.path(cache_dir, paste0(base, ".nc.zip")),
      file.path(cache_dir, paste0(base, ".nc.gz")),
      file.path(cache_dir, paste0(base, ".gz"))
    )
    unique(candidates)
  }
  
  # ------------------------------------------------------------------
  # Shared helper: resolve a file path to a plain .nc for the stack builder.
  # ------------------------------------------------------------------
  # Hoisted out of acquire_globsnow_swe() for the same reason as
  # generate_cache_candidates() above - shared by both GlobSnow acquisition
  # functions instead of being visible to only one of them.
  resolve_to_nc <- function(f, cache_dir) {
    if (!file.exists(f) || file.info(f)$size == 0) return(NA_character_)
    
    if (grepl("\\.nc$", f, ignore.case = TRUE) && !grepl("\\.(zip|gz)$", f, ignore.case = TRUE)) {
      return(f)
    }
    
    base <- sub("\\.(zip|gz)$", "", f, ignore.case = TRUE)
    target_nc <- if (grepl("\\.nc$", base, ignore.case = TRUE)) {
      base
    } else {
      paste0(base, ".nc")
    }
    
    if (file.exists(target_nc) && file.info(target_nc)$size > 0) {
      return(target_nc)
    }
    
    if (grepl("\\.zip$", f, ignore.case = TRUE)) {
      message(sprintf("    [EXTRACT] %s -> %s", basename(f), basename(target_nc)))
      tryCatch({
        utils::unzip(f, exdir = cache_dir)
        if (!file.exists(target_nc)) {
          extracted_ncs <- list.files(cache_dir, pattern = "\\.nc$", full.names = TRUE)
          base_pattern <- sub("\\.nc$", "", basename(target_nc), ignore.case = TRUE)
          matches <- extracted_ncs[grepl(base_pattern, basename(extracted_ncs), ignore.case = TRUE)]
          if (length(matches) > 0) {
            if (matches[1] != target_nc) file.rename(matches[1], target_nc)
          }
        }
      }, error = function(e) {
        message(sprintf("    [WARN] Failed to unzip %s: %s", basename(f), conditionMessage(e)))
      })
    } else if (grepl("\\.gz$", f, ignore.case = TRUE)) {
      message(sprintf("    [DECOMPRESS] %s -> %s", basename(f), basename(target_nc)))
      tryCatch({
        con_in <- gzfile(f, "rb")
        con_out <- file(target_nc, "wb")
        # 8MB read buffer instead of the previous 64KB - these GlobSnow
        # monthly files are only ~150-300KB compressed (a few MB at most
        # decompressed), so a 64KB chunk meant several R-level readBin()/
        # writeBin() round-trips per file for no reason; an 8MB buffer
        # finishes a typical file in a single iteration while still looping
        # safely (no truncation risk) if a future file is ever much larger.
        while (length(chunk <- readBin(con_in, "raw", 8L * 1024L * 1024L)) > 0) {
          writeBin(chunk, con_out)
        }
        close(con_in)
        close(con_out)
      }, error = function(e) {
        message(sprintf("    [WARN] Failed to decompress %s: %s", basename(f), conditionMessage(e)))
      })
    }
    
    if (file.exists(target_nc) && file.info(target_nc)$size > 0) {
      return(target_nc)
    }
    NA_character_
  }
  
  # ------------------------------------------------------------------
  # Shared helper: resolve MANY files to .nc IN PARALLEL.
  # ------------------------------------------------------------------
  # Both GlobSnow acquisition functions used to call resolve_to_nc() once per
  # file inside a plain for() loop - fine for a handful of files, but with
  # 100+ independent .gz files to decompress (the common case the first time
  # this script runs against a fully-populated local archive), that meant
  # 100+ sequential gzfile()/readBin()/writeBin() round-trips, none of them
  # overlapping, even though decompressing one file has zero dependency on
  # any other. This fans that work out across a PSOCK cluster (works on
  # Windows, unlike parallel::mclapply's fork-only approach, which is
  # important since this pipeline runs on a Windows machine per its D:/ paths
  # elsewhere in this script) instead of one process doing it all serially.
  #
  # resolve_to_nc itself is passed in and shipped to each worker as a plain
  # function argument (rather than via clusterExport() + guessing the right
  # environment to export it from) - this only works because resolve_to_nc
  # is self-contained (it only calls base/utils functions - file.exists,
  # gzfile, readBin, unzip, etc. - already available on any fresh R worker),
  # so serializing the function object itself is all that's needed.
  #
  # Falls back to a plain sequential loop - identical to the previous
  # behaviour - whenever there are too few files to be worth cluster startup
  # overhead, only 1 usable core is detected, or cluster setup itself fails
  # for any reason (e.g. a locked-down machine that disallows spawning child
  # R processes). Number of workers is capped via GLOBSNOW_DECOMPRESS_CORES
  # (env var; "1" forces sequential); default is (detected cores - 1) so the
  # machine keeps one core free, with an overall ceiling of 8 (little benefit
  # beyond that for this many small files, and it keeps memory/handle usage
  # bounded).
  resolve_to_nc_batch <- function(files, cache_dir, resolve_fun = resolve_to_nc,
                                  min_files_for_parallel = 4,
                                  max_workers = suppressWarnings(as.integer(Sys.getenv("GLOBSNOW_DECOMPRESS_CORES", "0")))) {
    if (length(files) == 0) return(character(0))
    
    # ---- Cheap, sequential, in-this-process pre-check: is this file ----
    # ---- ALREADY a usable, non-empty .nc with nothing to extract? ----
    # This mirrors resolve_to_nc()'s own first branch (a plain, existing,
    # non-zero-size .nc file returns itself immediately) plus its
    # already-resolved-target check for .zip/.gz inputs - but done here,
    # up front, for the WHOLE batch in one pass, so the file list handed
    # to the parallel path (if any) never includes files with nothing to do.
    is_already_resolved <- function(f) {
      if (!file.exists(f) || file.info(f)$size == 0) return(NA)  # missing/empty - let resolve_fun handle/report it
      
      if (grepl("\\.nc$", f, ignore.case = TRUE) && !grepl("\\.(zip|gz)$", f, ignore.case = TRUE)) {
        return(f)  # already a plain .nc - resolve_to_nc()'s first branch would just return this
      }
      
      base <- sub("\\.(zip|gz)$", "", f, ignore.case = TRUE)
      target_nc <- if (grepl("\\.nc$", base, ignore.case = TRUE)) base else paste0(base, ".nc")
      if (file.exists(target_nc) && file.info(target_nc)$size > 0) {
        return(target_nc)  # already decompressed/extracted in a previous run
      }
      NA_character_  # genuinely needs decompression/extraction this run
    }
    
    precheck <- vapply(files, function(f) {
      r <- is_already_resolved(f)
      if (is.na(r)) NA_character_ else r
    }, character(1))
    
    already_done_idx <- !is.na(precheck)
    n_already <- sum(already_done_idx)
    if (n_already > 0) {
      message(sprintf(
        "  [OK] %d/%d file(s) already resolved to .nc with no work needed this run (skipped before any cluster/sequential dispatch).",
        n_already, length(files)
      ))
    }
    
    files_todo <- files[!already_done_idx]
    out <- character(length(files))
    out[already_done_idx] <- precheck[already_done_idx]
    
    if (length(files_todo) == 0) {
      return(out)  # every file was already resolved - no cluster ever started
    }
    
    # ---- From here down: identical logic to before, but only for the ----
    # ---- subset of files that actually need decompression/extraction. ----
    if (is.na(max_workers) || max_workers <= 0) {
      max_workers <- tryCatch({
        dc <- parallel::detectCores(logical = TRUE)
        if (is.na(dc) || dc < 1) 1L else max(1L, dc - 1L)
      }, error = function(e) 1L)
    }
    max_workers <- min(max_workers, 8L, length(files_todo))
    
    run_sequential <- function(todo) {
      res_out <- character(length(todo))
      for (i in seq_along(todo)) {
        res <- tryCatch(resolve_fun(todo[i], cache_dir), error = function(e) NA_character_)
        res_out[i] <- if (is.null(res)) NA_character_ else res
        if (is.na(res_out[i])) {
          message(sprintf("    [WARN] Could not resolve %s to .nc - skipping from stack.", basename(todo[i])))
        }
      }
      res_out
    }
    
    if (length(files_todo) < min_files_for_parallel || max_workers <= 1) {
      resolved_todo <- run_sequential(files_todo)
    } else {
      cl <- tryCatch(parallel::makeCluster(max_workers, type = "PSOCK"), error = function(e) NULL)
      if (is.null(cl)) {
        message("  [INFO] Could not start a parallel cluster for decompression/extraction - falling back to sequential.")
        resolved_todo <- run_sequential(files_todo)
      } else {
        on.exit(parallel::stopCluster(cl), add = TRUE)
        message(sprintf(
          "  [OK] Decompressing/extracting %d file(s) needing real work across %d parallel worker(s) (set GLOBSNOW_DECOMPRESS_CORES=1 to force sequential)...",
          length(files_todo), max_workers
        ))
        results <- tryCatch(
          parallel::parLapply(cl, files_todo, function(f, resolve_fun, cache_dir) {
            tryCatch(resolve_fun(f, cache_dir), error = function(e) NA_character_)
          }, resolve_fun = resolve_fun, cache_dir = cache_dir),
          error = function(e) {
            message("  [WARN] Parallel decompression failed (", conditionMessage(e), ") - falling back to sequential.")
            NULL
          }
        )
        if (is.null(results)) {
          resolved_todo <- run_sequential(files_todo)
        } else {
          resolved_todo <- vapply(results, function(x) if (is.null(x) || is.na(x)) NA_character_ else x, character(1))
          for (i in seq_along(resolved_todo)) {
            if (is.na(resolved_todo[i])) {
              message(sprintf("    [WARN] Could not resolve %s to .nc - skipping from stack.", basename(files_todo[i])))
            }
          }
        }
      }
    }
    
    out[!already_done_idx] <- resolved_todo
    out
  }
  # ---- 1) GlobSnow SWE (monthly, bias-corrected) - fill in missing months only ----
  # NOTE: this function used to be mistakenly named/defined as
  # "acquire_globsnow_se_mfsc" (a copy-paste of the SE/MFSC acquisition below,
  # pointed at the wrong base URL/layout) - it is the SWE fetcher referenced
  # by the "acquire_globsnow_swe()" call further down, and is fixed here to
  # actually target the SWE archive: GLOBSNOW_SWE_BASE_URL (no per-year
  # "MFSC/" subfolder - the SWE archive is one flat NetCDF4/ directory), and
  # the "YYYYMM_northern_hemisphere_monthly_biascorrected_swe_0.25grid.nc"
  # naming convention documented above. The local-cache-first logic
  # (generate_cache_candidates()/resolve_to_nc(), matching .nc/.zip/.gz,
  # now shared helpers defined just above) is otherwise unchanged from the
  # original body - it was already correct, just mislabeled and pointed at
  # the wrong product.
  acquire_globsnow_swe <- function(
    cache_dir = file.path("snow_data", "globsnow", "swe"),
    out_stack_path = file.path("monthly_data_direct", "globsnow_swe_monthly.nc"),
    start = study_start,
    end = study_end,
    check_remote = as.logical(Sys.getenv("GLOBSNOW_SWE_CHECK_REMOTE", "FALSE"))) {
    message("\n===== GlobSnow SWE (monthly, bias-corrected) acquisition =====")
    ensure_dir(cache_dir)
    ensure_dir(dirname(out_stack_path))
    if (is.na(check_remote)) check_remote <- FALSE
    
    # GlobSnow v3.0 monthly bias-corrected SWE covers a much longer record
    # than the SE/MFSC product (which is clamped to 2003-2012 below) - no
    # equivalent date clamp is applied here beyond the overall study window.
    start <- as.Date(start)
    end   <- as.Date(end)
    if (start > end) {
      message("  [OK] Requested period is outside the study window.")
      return(invisible(NULL))
    }
    
    yrs <- seq(as.integer(format(start, "%Y")), as.integer(format(end, "%Y")))
    
    downloaded <- character(0)
    n_ok <- 0
    n_404 <- 0
    n_failed <- 0
    n_already <- 0
    
    # ---- Local scan FIRST, independent of whether the remote listing works ----
    # This used to only check "is this specific filename (from the remote
    # listing) already in cache_dir" - so if the remote listing itself failed
    # (e.g. GLOBSNOW_SWE_BASE_URL 404ing because FMI reorganized/retired that
    # exact path), data_hrefs came back empty, the per-file loop below never
    # ran at all, and files already sitting in cache_dir were never even
    # looked at - "0 already on disk" even with dozens of real local files
    # present. This scans cache_dir directly for files matching the SWE
    # naming convention, so already-downloaded months are found and used
    # regardless of whether globsnow.info is reachable this run.
    local_swe_pattern <- "^[0-9]{6}.*\\.(nc|zip|gz)$"
    local_files <- list.files(cache_dir, pattern = local_swe_pattern, full.names = TRUE)
    if (length(local_files) > 0) {
      local_ym <- suppressWarnings(as.integer(substr(basename(local_files), 1, 6)))
      start_ym0 <- as.integer(format(start, "%Y%m"))
      end_ym0   <- as.integer(format(end, "%Y%m"))
      keep_local <- !is.na(local_ym) & local_ym >= start_ym0 & local_ym <= end_ym0
      local_files_in_window <- local_files[keep_local]
      n_already <- n_already + length(local_files_in_window)
      downloaded <- c(downloaded, local_files_in_window)
      message(sprintf(
        "  [OK] Found %d local GlobSnow SWE file(s) already in '%s' within the study window (of %d total matching the naming pattern there).",
        length(local_files_in_window), cache_dir, length(local_files)
      ))
    }
    
    # The SWE archive (unlike SE/MFSC) is a single flat directory covering the
    # whole record - "YYYYMM_..._swe_0.25grid.nc" per month, not one
    # subfolder per year - so this lists it once rather than looping a
    # per-year "MFSC/"-style subfolder that doesn't exist for this product.
    dir_url <- GLOBSNOW_SWE_BASE_URL
    
    # ---- Decide whether the remote directory listing is even worth trying ----
    # Previously this was queried UNCONDITIONALLY, regardless of what the local
    # scan just found - so a run with every month already on disk still burned
    # 4 retries of DNS/HTTP backoff against globsnow.info before falling back
    # to the local files anyway (up to ~100s wasted, worse still if the host
    # simply doesn't resolve). Mirrors the local-found-so-skip-network pattern
    # already used for CMC/NSIDC-0447 above: if local files already cover the
    # requested window, skip the network call entirely by default. Set
    # GLOBSNOW_SWE_CHECK_REMOTE=TRUE (env var) to still query globsnow.info for
    # newly published months even when local files already exist - useful the
    # first time you run this after GlobSnow might have added a new month, or
    # periodically if you want to check for updates.
    if (n_already > 0 && !check_remote) {
      message(sprintf(
        "  [OK] %d local file(s) already cover the study window and GLOBSNOW_SWE_CHECK_REMOTE is not TRUE - skipping the globsnow.info directory listing entirely. Set GLOBSNOW_SWE_CHECK_REMOTE=TRUE (env var) to also check for newly published months.",
        n_already
      ))
      data_hrefs <- character(0)
    } else {
      if (n_already == 0) {
        message("  No local GlobSnow SWE files found - querying globsnow.info for the study window.")
      }
      hrefs <- list_https_dir(dir_url)
      data_hrefs <- hrefs[grepl("\\.(nc|gz|zip)$", hrefs, ignore.case = TRUE)]
      
      if (length(data_hrefs) == 0) {
        message(
          "  [INFO] Remote directory listing at ", dir_url, " returned no files ",
          "(could be a genuine empty listing, a changed URL, or the site being ",
          "unreachable - see any HTTP status/retry messages above). Continuing ",
          "with whatever was already found locally above."
        )
      } else {
        # Restrict to files whose YYYYMM prefix falls within [start, end] -
        # mirrors the year-filtering the SE/MFSC loop got "for free" from its
        # per-year subfolders, now done by filename instead since everything
        # lives in one directory here.
        fn_ym <- suppressWarnings(as.integer(substr(data_hrefs, 1, 6)))
        start_ym <- as.integer(format(start, "%Y%m"))
        end_ym   <- as.integer(format(end, "%Y%m"))
        in_window <- !is.na(fn_ym) & fn_ym >= start_ym & fn_ym <= end_ym
        # Keep undated matches too (e.g. a readme swept up by the file pattern)
        # rather than silently dropping anything with an unparseable prefix.
        data_hrefs <- data_hrefs[in_window | is.na(fn_ym)]
      }
      
      if (length(data_hrefs) == 0) {
        message(sprintf(
          "  [INFO] No GlobSnow SWE files found for %d-%d in the remote archive listing.",
          min(yrs), max(yrs)
        ))
      }
    }
    
    for (fn in data_hrefs) {
      candidates <- generate_cache_candidates(fn, cache_dir)
      existing_file <- NULL
      for (candidate in candidates) {
        if (file.exists(candidate) && file.info(candidate)$size > 0) {
          existing_file <- candidate
          break
        }
      }
      
      if (!is.null(existing_file)) {
        n_already <- n_already + 1
        downloaded <- c(downloaded, existing_file)
        next
      }
      
      dest <- file.path(cache_dir, fn)
      url <- paste0(dir_url, fn)
      
      probe <- tryCatch(httr::HEAD(url, timeout(20)), error = function(e) NULL)
      if (!is.null(probe) && httr::status_code(probe) == 404) {
        n_404 <- n_404 + 1
        next
      }
      
      ok <- tryCatch({
        download_file_auth(url, dest, user = "", pass = "")
        file.exists(dest) && file.info(dest)$size > 0
      }, error = function(e) {
        message(sprintf("    [FAIL] %s: %s", fn, conditionMessage(e)))
        if (file.exists(dest)) unlink(dest)
        FALSE
      })
      
      if (ok) {
        n_ok <- n_ok + 1
        message("  [DL] ", fn)
        downloaded <- c(downloaded, dest)
      } else {
        n_failed <- n_failed + 1
      }
    }
    
    # De-duplicate: a file the local scan already found could also appear in
    # the remote listing's loop above (both branches check the same cache
    # candidates), which would otherwise double-count it in n_already and
    # leave a duplicate path in `downloaded`.
    dup_paths <- duplicated(downloaded)
    if (any(dup_paths)) {
      n_already <- n_already - sum(dup_paths)
      downloaded <- downloaded[!dup_paths]
    }
    
    message(sprintf(
      "  [OK] GlobSnow SWE: %d already on disk (including local .zip/.gz files), %d newly downloaded, %d not published (404), %d failed.",
      n_already, n_ok, n_404, n_failed
    ))
    
    if (length(downloaded) == 0) {
      warning("No GlobSnow SWE files available locally or via download.")
      return(invisible(NULL))
    }
    
    resolved_nc_files <- resolve_to_nc_batch(downloaded, cache_dir)
    resolved_nc_files <- unique(resolved_nc_files[!is.na(resolved_nc_files)])
    downloaded <- resolved_nc_files
    
    if (length(downloaded) == 0) {
      warning("No GlobSnow SWE files could be resolved to .nc for stacking.")
      return(invisible(NULL))
    }
    
    needs_rebuild <- TRUE
    if (file.exists(out_stack_path)) {
      existing_stack <- tryCatch(terra::rast(out_stack_path), error = function(e) NULL)
      if (!is.null(existing_stack) && terra::nlyr(existing_stack) >= length(downloaded)) {
        needs_rebuild <- FALSE
        message(sprintf("  [OK] '%s' already contains %d layer(s).", out_stack_path, terra::nlyr(existing_stack)))
      }
    }
    
    if (needs_rebuild) {
      ym <- regmatches(basename(downloaded), regexpr("(19|20)\\d{4}", basename(downloaded)))
      ord <- order(ym)
      downloaded <- downloaded[ord]
      ym <- ym[ord]
      dates <- as.Date(paste0(ym, "01"), format = "%Y%m%d")
      
      rasters <- lapply(downloaded, function(f) tryCatch(terra::rast(f), error = function(e) NULL))
      ok_idx <- !vapply(rasters, is.null, logical(1))
      if (sum(ok_idx) == 0) { warning("None of the GlobSnow SWE files could be opened."); return(invisible(NULL)) }
      rasters <- rasters[ok_idx]
      dates <- dates[ok_idx]
      rasters <- lapply(rasters, function(r) if (terra::nlyr(r) > 1) r[[1]] else r)
      
      template <- rasters[[1]]
      rasters <- lapply(rasters, function(r) {
        if (terra::compareGeom(r, template, stopOnError = FALSE)) r else terra::resample(r, template, method = "near")
      })
      stack <- terra::rast(rasters)
      names(stack) <- format(dates, "%Y-%m")
      terra::time(stack) <- dates
      terra::writeCDF(stack, out_stack_path, varname = "SWE", overwrite = TRUE)
      message(sprintf("  [OK] Wrote %s (%d monthly layers, %s to %s).",
                      out_stack_path, terra::nlyr(stack),
                      format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m")))
    }
    
    invisible(list(cache_files = downloaded, out_stack_path = out_stack_path))
  }
  
  # ---- 2) GlobSnow SE / MFSC (monthly fractional snow cover) ----
  # Downloads whatever MFSC files are missing locally, then stitches them into
  # a single time-stacked NetCDF at monthly_data_direct/snow_cover_monthly.nc -
  # the exact path/shape H8SWEI_Snow.R's scf_file expects (one multi-band
  # stack with a time dimension, matching how snow_depth_water_equivalent_monthly.nc
  # is already structured), rather than one file per month.
  #
  # NOTE: this replaces an earlier, buggy duplicate definition of
  # acquire_globsnow_se_mfsc() that (a) had no local-cache check at all - it
  # queried globsnow.info's directory listing directly on every call, even
  # when the SE/MFSC files were already sitting in cache_dir - and (b) only
  # matched "\\.nc$" in that listing, silently ignoring any already-downloaded
  # .zip/.gz copies. That combination is what produced "No MFSC files found"
  # for every year (2003-2012) followed by "No GlobSnow SE/MFSC files
  # available locally or via download" even when local SE data existed:
  # local files were never checked, and the remote listing itself came back
  # empty for this run. This version restores the local-cache-first check
  # (generate_cache_candidates()/resolve_to_nc(), matching .nc/.zip/.gz) that
  # acquire_globsnow_swe() above already uses, applied to the SE/MFSC
  # per-year "<year>/MFSC/" archive layout instead of the SWE product's flat
  # directory.
  acquire_globsnow_se_mfsc <- function(
    cache_dir = file.path("snow_data", "globsnow", "se_mfsc"),
    out_stack_path = file.path("monthly_data_direct", "snow_cover_monthly.nc"),
    start = study_start,
    end = study_end,
    check_remote = as.logical(Sys.getenv("GLOBSNOW_SE_CHECK_REMOTE", "FALSE"))) {
    message("\n===== GlobSnow SE - MFSC (Monthly Aggregated Fractional Snow Cover) acquisition =====")
    ensure_dir(cache_dir)
    ensure_dir(dirname(out_stack_path))
    if (is.na(check_remote)) check_remote <- FALSE
    
    FIRST_MFSC_DATE <- as.Date("2003-01-01")
    LAST_MFSC_DATE  <- as.Date("2012-04-30")
    start <- max(as.Date(start), FIRST_MFSC_DATE)
    end   <- min(as.Date(end), LAST_MFSC_DATE)
    if (start > end) {
      message("  [OK] Requested period is outside the available GlobSnow SE archive.")
      return(invisible(NULL))
    }
    
    yrs <- seq(as.integer(format(start, "%Y")), as.integer(format(end, "%Y")))
    
    downloaded <- character(0)
    n_ok <- 0
    n_404 <- 0
    n_failed <- 0
    n_already <- 0
    
    # ---- Local scan FIRST, independent of whether globsnow.info is reachable ----
    # Mirrors the local-scan-first logic in acquire_globsnow_swe() above (and
    # the pre-existing CMC/NSIDC-0447 pattern) - this function previously had
    # NO bulk local check at all: it went straight to querying globsnow.info
    # per year, and only ever consulted the local cache for a filename it had
    # already gotten back from a successful remote listing. If the remote
    # listing failed outright (e.g. a DNS resolution failure), already-
    # downloaded local SE/MFSC files were never even looked at, and every
    # single year in [start, end] paid its own retry/backoff cost for nothing.
    # Filenames follow "GlobSnow_SE_MFSC_L3B_NH_YYYYMM_v2.1_1.nc" - the year is
    # extracted the same way the stack-building step further down already
    # does (first "(19|20)dddd" match anywhere in the filename), rather than
    # assuming it's at a fixed position, since it isn't at the start of the
    # name here the way it is for the SWE product.
    local_files <- list.files(cache_dir, pattern = "\\.(nc|zip|gz)$", full.names = TRUE, recursive = TRUE)
    if (length(local_files) > 0) {
      local_ym_chr <- regmatches(basename(local_files), regexpr("(19|20)\\d{4}", basename(local_files)))
      local_year <- suppressWarnings(as.integer(substr(local_ym_chr, 1, 4)))
      keep_local <- !is.na(local_year) &
        local_year >= as.integer(format(start, "%Y")) &
        local_year <= as.integer(format(end, "%Y"))
      local_files_in_window <- local_files[keep_local]
      n_already <- n_already + length(local_files_in_window)
      downloaded <- c(downloaded, local_files_in_window)
      message(sprintf(
        "  [OK] Found %d local GlobSnow SE/MFSC file(s) already in '%s' within the study window (of %d total matching the naming pattern there).",
        length(local_files_in_window), cache_dir, length(local_files)
      ))
    }
    
    # ---- Decide whether the per-year remote listing loop is even worth trying ----
    # Same reasoning as acquire_globsnow_swe() above: skip globsnow.info
    # entirely when local files already cover the window and the caller
    # hasn't asked to check for anything new. Without this, an unreachable
    # globsnow.info meant up to 10 separate years (2003-2012), each paying its
    # own 4-attempt DNS-retry/backoff cost, before falling back to local files
    # that were sitting there the whole time. Set GLOBSNOW_SE_CHECK_REMOTE=TRUE
    # (env var) to still query globsnow.info even when local files exist.
    if (n_already > 0 && !check_remote) {
      message(sprintf(
        "  [OK] %d local file(s) already cover the study window and GLOBSNOW_SE_CHECK_REMOTE is not TRUE - skipping the globsnow.info directory listing entirely (all %d year(s): %d-%d). Set GLOBSNOW_SE_CHECK_REMOTE=TRUE (env var) to also check for newly published months.",
        n_already, length(yrs), min(yrs), max(yrs)
      ))
    } else {
      if (n_already == 0) {
        message("  No local GlobSnow SE/MFSC files found - querying globsnow.info for each year in the study window.")
      }
      for (yr in yrs) {
        dir_url <- sprintf("%s%d/MFSC/", GLOBSNOW_SE_BASE_URL, yr)
        hrefs <- list_https_dir(dir_url)
        
        data_hrefs <- hrefs[grepl("\\.(nc|gz|zip)$", hrefs, ignore.case = TRUE)]
        
        if (length(data_hrefs) == 0) {
          message(sprintf("  [INFO] No MFSC files found for %d.", yr))
          next
        }
        
        for (fn in data_hrefs) {
          # ---- Check the local cache FIRST, before touching the network ----
          # generate_cache_candidates()/resolve_to_nc() are defined once,
          # above, inside acquire_globsnow_swe() - shared here rather than
          # redefined, since both functions cache/resolve GlobSnow NetCDFs
          # the same way.
          candidates <- generate_cache_candidates(fn, cache_dir)
          existing_file <- NULL
          for (candidate in candidates) {
            if (file.exists(candidate) && file.info(candidate)$size > 0) {
              existing_file <- candidate
              break
            }
          }
          
          if (!is.null(existing_file)) {
            n_already <- n_already + 1
            downloaded <- c(downloaded, existing_file)
            next
          }
          
          dest <- file.path(cache_dir, fn)
          url <- paste0(dir_url, fn)
          
          probe <- tryCatch(httr::HEAD(url, timeout(20)), error = function(e) NULL)
          if (!is.null(probe) && httr::status_code(probe) == 404) {
            n_404 <- n_404 + 1
            next
          }
          
          ok <- tryCatch({
            download_file_auth(url, dest, user = "", pass = "")
            file.exists(dest) && file.info(dest)$size > 0
          }, error = function(e) {
            message(sprintf("    [FAIL] %s: %s", fn, conditionMessage(e)))
            if (file.exists(dest)) unlink(dest)
            FALSE
          })
          
          if (ok) {
            n_ok <- n_ok + 1
            message("  [DL] ", fn)
            downloaded <- c(downloaded, dest)
          } else {
            n_failed <- n_failed + 1
          }
        }
      }
    }
    
    # De-duplicate: a file the local pre-scan already found could also appear
    # again via the per-year loop's own local-cache check above (possible when
    # check_remote = TRUE alongside existing local files), which would
    # otherwise double-count it in n_already and leave a duplicate path in
    # `downloaded`. Mirrors the same de-duplication in acquire_globsnow_swe().
    dup_paths <- duplicated(downloaded)
    if (any(dup_paths)) {
      n_already <- n_already - sum(dup_paths)
      downloaded <- downloaded[!dup_paths]
    }
    
    message(sprintf(
      "  [OK] GlobSnow SE/MFSC: %d already on disk (including local .zip/.gz files), %d newly downloaded, %d not published (404), %d failed.",
      n_already, n_ok, n_404, n_failed
    ))
    
    if (length(downloaded) == 0) {
      warning("No GlobSnow SE/MFSC files available locally or via download.")
      return(invisible(NULL))
    }
    
    resolved_nc_files <- resolve_to_nc_batch(downloaded, cache_dir)
    resolved_nc_files <- unique(resolved_nc_files[!is.na(resolved_nc_files)])
    downloaded <- resolved_nc_files
    
    if (length(downloaded) == 0) {
      warning("No GlobSnow SE/MFSC files could be resolved to .nc for stacking.")
      return(invisible(NULL))
    }
    
    needs_rebuild <- TRUE
    if (file.exists(out_stack_path)) {
      existing_stack <- tryCatch(terra::rast(out_stack_path), error = function(e) NULL)
      if (!is.null(existing_stack) && terra::nlyr(existing_stack) >= length(downloaded)) {
        needs_rebuild <- FALSE
        message(sprintf("  [OK] '%s' already contains %d layer(s).", out_stack_path, terra::nlyr(existing_stack)))
      }
    }
    
    if (needs_rebuild) {
      ym <- regmatches(basename(downloaded), regexpr("(19|20)\\d{4}", basename(downloaded)))
      ord <- order(ym)
      downloaded <- downloaded[ord]
      ym <- ym[ord]
      dates <- as.Date(paste0(ym, "01"), format = "%Y%m%d")
      
      rasters <- lapply(downloaded, function(f) tryCatch(terra::rast(f), error = function(e) NULL))
      ok_idx <- !vapply(rasters, is.null, logical(1))
      if (sum(ok_idx) == 0) { warning("None of the MFSC files could be opened."); return(invisible(NULL)) }
      rasters <- rasters[ok_idx]
      dates <- dates[ok_idx]
      rasters <- lapply(rasters, function(r) if (terra::nlyr(r) > 1) r[[1]] else r)
      
      template <- rasters[[1]]
      rasters <- lapply(rasters, function(r) {
        if (terra::compareGeom(r, template, stopOnError = FALSE)) r else terra::resample(r, template, method = "near")
      })
      stack <- terra::rast(rasters)
      names(stack) <- format(dates, "%Y-%m")
      terra::time(stack) <- dates
      terra::writeCDF(stack, out_stack_path, varname = "SCF", overwrite = TRUE)
      message(sprintf("  [OK] Wrote %s (%d monthly layers, %s to %s).",
                      out_stack_path, terra::nlyr(stack),
                      format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m")))
    }
    
    invisible(list(cache_files = downloaded, out_stack_path = out_stack_path))
  }
  # ----------------------------------------------------------------------------
  # Execute acquisitions
  # ----------------------------------------------------------------------------
  # Sentinel-1 automatic acquisition; RCM supported through optional local
  # manifests or directories because that archive is not uniformly exposed
  # through a single public API. C-SNOW is now ingested directly from the
  # local daily-NetCDF archive via acquire_csnow_local() (defined above),
  # rather than through the generic manifest/local-copy helper.
  # NOTE: the inventory below used to hardcode Sentinel-1's status to
  # "SKIP - see log" unconditionally, regardless of whether this call actually
  # produced files (it does, via ASF or the CDSE fallback, return a character
  # vector of the files acquired) - so even a fully successful Sentinel-1 run
  # would have been misreported as skipped. Capturing the result here so the
  # final inventory can report what actually happened.
  if (RUN_SENTINEL1_RAW_ASF) {
    sentinel1_files <- acquire_sentinel1_asf()
  } else {
    message(
      "\n===== Sentinel-1 (ASF) raw acquisition SKIPPED =====\n",
      "  RUN_SENTINEL1_RAW_ASF is FALSE - using the local C-SNOW product ",
      "(already-processed Sentinel-1 snow depth) instead of downloading raw, ",
      "un-terrain-corrected Sentinel-1 GRD scenes and reprocessing them locally. ",
      "Set RUN_SENTINEL1_RAW_ASF=TRUE (env var) if C-SNOW's coverage/validation ",
      "diagnostics below show it isn't usable for the Nechako basin and you need ",
      "to build your own SAR-based retrieval instead."
    )
    sentinel1_files <- character(0)
  }
  
  acquire_rcm_or_csnow(
    provider_name = "RCM",
    manifest_csv  = Sys.getenv("RCM_MANIFEST_CSV", unset = ""),
    local_dir     = Sys.getenv("RCM_LOCAL_DIR", unset = "")
  )
  
  csnow_result <- tryCatch(
    acquire_csnow_local(),
    error = function(e) {
      message("[SKIP] C-SNOW local ingestion failed: ", conditionMessage(e))
      NULL
    }
  )
  
  # Validation products from NSIDC/CMR.
  cmc_dir      <- file.path("snow_data", "cmc")
  globsnow_dir <- file.path("snow_data", "globsnow")
  cansise_dir  <- file.path("snow_data", "cansise")
  ensure_dir(cmc_dir); ensure_dir(globsnow_dir); ensure_dir(cansise_dir)
  
  # UPDATE: cmr_download_wrapper() (CMR granule search) previously returned
  # "0 granules found" for BOTH of these collections even after the
  # zero-padded-version fix confirmed the short_name/version are correct
  # (verify_cmr_collection() finds NSIDC-0447 v1 and NSIDC-0668 v2 with no
  # mismatch). That's because neither collection is actually granule-indexed
  # in CMR at all - both datasets' own NSIDC landing pages document a
  # "Programmatic Access Guide for Data on daacdata.apps" (direct HTTPS
  # file-system access) rather than a CMR-searchable granule collection. This
  # is the exact same category of problem already identified and worked
  # around below for NSIDC-0595/GlobSnow ("NSIDC does not archive this
  # dataset") - CMR was never going to return anything here regardless of
  # short_name/version, so cmr_download_wrapper() has been replaced with
  # direct HTTPS acquisition via acquire_nsidc_https_dir() (defined near the
  # top of Part B, alongside the other shared download helpers).
  # Each acquisition is wrapped in tryCatch: NSIDC-0447 and NSIDC-0668 remain
  # independent of each other and of everything below them (GlobSnow's
  # message, OpenLandMap, the final inventory write) - a caught failure here
  # is reported and the run continues; check the message above the "[SKIP]"
  # line for the reason.
  cmr_ok <- c(NSIDC_0447 = FALSE, NSIDC_0668 = FALSE)
  
  # NSIDC-0447 (CMC): this used to target Snow_Depth/Snow_Depth_Daily_Values/
  # GeoTIFF/ (one GeoTIFF per year, e.g. cmc_sdepth_dly_2020_v01.2.tif) - but
  # that directory holds daily SNOW DEPTH (cm), not SWE. CMC/NSIDC-0447 does
  # NOT produce a daily SWE product at all: per the dataset's User Guide
  # (Sec. 1.5.4), SWE is only ever estimated at MONTHLY resolution, for
  # Oct-Jun of each year (Jul-Sep is skipped for lack of snow-density obs),
  # by applying a fixed monthly snow-density lookup table to the monthly
  # mean snow depth. So a "daily SWE" download target was never reachable -
  # the only genuine CMC SWE source is SWE/SWE_Monthly_Averages/, which
  # (conveniently) ships as a SINGLE file covering the entire 1998-2020
  # record: cmc_swe_mly_1998to2020_v01.2.zip. That's also monthly, matching
  # every other series in this pipeline, so no daily-to-monthly aggregation
  # step is needed downstream either.
  #
  # Files already sitting in cmc_dir (e.g. manually downloaded from
  # https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0447_CMC_snow_depth_v01/
  # via a browser, as direct-HTTPS listing from R has been unreliable here)
  # are used as-is and NOT re-downloaded. Only if no matching file is found
  # locally does this fall back to listing SWE_Monthly_Averages/ over HTTPS.
  message("\n===== CMC Monthly SWE (NSIDC-0447) =====")
  cmc_swe_pattern <- "^cmc_swe_mly_[0-9]{4}to[0-9]{4}_v[0-9.]+\\.zip$"
  cmr_ok["NSIDC_0447"] <- tryCatch({
    local_swe <- list.files(cmc_dir, pattern = cmc_swe_pattern, full.names = TRUE)
    if (length(local_swe) > 0) {
      message("  [OK] Found manually-downloaded CMC monthly SWE file already in ",
              cmc_dir, " - skipping download:")
      for (f in local_swe) message("    ", basename(f))
      TRUE
    } else {
      got <- acquire_nsidc_https_dir(
        dir_url    = "https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0447_CMC_snow_depth_v01/SWE/SWE_Monthly_Averages/",
        out_dir    = cmc_dir,
        file_pattern = cmc_swe_pattern,
        # This product ships as one file for the whole record (not one file
        # per year), so there's nothing to filter by year - NULL/NULL makes
        # acquire_nsidc_https_dir() skip its year filter and take whatever
        # matches file_pattern in the listed directory.
        start_year = NULL,
        end_year   = NULL
      )
      length(got) > 0
    }
  }, error = function(e) {
    message("[SKIP] NSIDC-0447 (CMC) monthly SWE acquisition failed: ", conditionMessage(e))
    FALSE
  })
  
  # NSIDC-0668 (CanSISE): same access pattern as NSIDC-0447 above (a
  # daacdata.apps.nsidc.org HTTPS file system, not a CMR-searchable granule
  # collection - confirmed via its landing page's own "Programmatic Access
  # Guide for Data on daacdata.apps" section), but the exact
  # "pub/DATASETS/nsidc0668_..." subdirectory name isn't published anywhere
  # that can be read without an Earthdata-authenticated browser session, so
  # hand-typing it risks the same "guessed a URL, could silently 404 or hit
  # the wrong dataset" problem this script already avoids elsewhere.
  #
  # FIX: discover the real daacdata.apps.nsidc.org base path from CMR's own
  # collection metadata instead of guessing it - the same "detect, don't
  # guess" approach already used for the WFS station-ID column (Part A) and
  # CanSWE's swe/snd column names (Part B) above. CMR's UMM-JSON collection
  # record exposes a RelatedUrls list, and for legacy NSIDC/ECS collections
  # like this one that list includes the actual "GET DATA" HTTPS file-system
  # URL (the type/subtype CMR uses for this varies by collection, so this
  # matches on the URL's own shape - it lives under daacdata.apps.nsidc.org -
  # rather than assuming one specific type/subtype label). If CMR's metadata
  # doesn't expose it for some reason, this falls back to requiring
  # CANSISE_HTTPS_DIR exactly as before, so there is still a manual escape
  # hatch, but manual guessing is no longer the ONLY path.
  discover_nsidc_https_base <- function(short_name, query_timeout_sec = 30) {
    resp <- tryCatch(
      httr::GET("https://cmr.earthdata.nasa.gov/search/collections.umm_json",
                query = list(short_name = short_name), timeout(query_timeout_sec)),
      error = function(e) NULL
    )
    if (is.null(resp) || httr::http_error(resp)) return(NA_character_)
    parsed <- tryCatch(httr::content(resp, as = "parsed"), error = function(e) NULL)
    items <- parsed$items
    if (is.null(items) || length(items) == 0) return(NA_character_)
    
    urls <- character(0)
    for (item in items) {
      related <- item$umm$RelatedUrls
      if (is.null(related)) next
      for (ru in related) {
        u <- ru$URL %||% NA_character_
        if (!is.na(u) && grepl("daacdata\\.apps\\.nsidc\\.org", u, ignore.case = TRUE)) {
          urls <- c(urls, u)
        }
      }
    }
    if (length(urls) == 0) return(NA_character_)
    
    # Prefer a URL that already looks like a directory (ends in "/"); a
    # collection-level RelatedUrls entry is usually the top-level dataset
    # folder rather than a specific file, but normalize either way by
    # trimming back to its parent directory if it looks like a file path.
    u <- urls[1]
    if (!grepl("/$", u)) u <- sub("/[^/]*$", "/", u)
    u
  }
  
  message("\n===== CanSISE Ensemble SWE (NSIDC-0668) =====")
  # NOTE: CanSISE is OPTIONAL - neither H8SWEI_Snow_1/2 nor
  # Nechako_Snow_Bridge.R reference CanSISE/NSIDC-0668 anywhere. A failure or
  # skip here does not block or degrade any downstream step; it's logged in
  # detail below purely so a real, fixable auth problem doesn't go unnoticed,
  # not because the rest of the pipeline is waiting on it.
  message("  [INFO] CanSISE is optional - not required by H8SWEI_Snow.R or Nechako_Snow_Bridge.R. A skip here does not block the rest of the pipeline.")
  cansise_https_dir <- Sys.getenv("CANSISE_HTTPS_DIR", unset = "")
  if (!nzchar(cansise_https_dir)) {
    message("  Attempting to discover the daacdata.apps.nsidc.org path for NSIDC-0668 from CMR metadata...")
    discovered <- discover_nsidc_https_base("CanSISE_Snow")
    if (is.na(discovered)) discovered <- discover_nsidc_https_base("NSIDC-0668")
    if (!is.na(discovered)) {
      cansise_https_dir <- discovered
      message("  [OK] Discovered CanSISE HTTPS directory from CMR: ", cansise_https_dir)
    }
  }
  if (!nzchar(cansise_https_dir)) {
    message(
      "  [SKIP] Could not auto-discover the daacdata.apps.nsidc.org path for ",
      "NSIDC-0668 from CMR metadata, and CANSISE_HTTPS_DIR is not set in .env. ",
      "To enable manually: (1) visit https://nsidc.org/data/nsidc-0668/versions/2 ",
      "and open its 'Programmatic Access Guide for Data on daacdata.apps' ",
      "section to get the exact daacdata.apps.nsidc.org/pub/DATASETS/... ",
      "directory URL (same pattern as NSIDC-0447 above); (2) set ",
      "CANSISE_HTTPS_DIR to that URL (including a trailing /) in .env; ",
      "(3) re-run. Skipping automatic acquisition until then."
    )
    cmr_ok["NSIDC_0668"] <- FALSE
  } else {
    cmr_ok["NSIDC_0668"] <- tryCatch({
      got <- acquire_nsidc_https_dir(
        dir_url    = cansise_https_dir,
        out_dir    = cansise_dir,
        file_pattern = "\\.(nc|tif|zip)$",
        start_year = as.integer(format(max(study_start, as.Date("1981-01-01")), "%Y")),
        end_year   = as.integer(format(min(study_end, as.Date("2010-12-31")), "%Y"))
      )
      length(got) > 0
    }, error = function(e) {
      message("[SKIP] NSIDC-0668 (CanSISE) direct-HTTPS acquisition failed: ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(cmr_ok["NSIDC_0668"])) {
      message(
        "  [NOTE] If this discovered/used directory turns out to be wrong ",
        "(e.g. an unexpected file naming pattern in the downloaded files), ",
        "set CANSISE_HTTPS_DIR explicitly in .env to override auto-discovery ",
        "next run - see https://nsidc.org/data/nsidc-0668/versions/2 for the ",
        "authoritative path."
      )
    }
  }
  
  # NSIDC-0595 (ESA GlobSnow SWE) is a landing-page-only record at NSIDC - its
  # own dataset page states "NSIDC does not archive these data. Please visit
  # the data set website for data access and citation guidance." There is
  # therefore nothing for CMR to return here regardless of short_name/version,
  # which is why this returned 0 granules alongside the other two. The actual
  # GlobSnow files (both SWE and SE/snow-extent) are distributed directly by
  # FMI, not through NSIDC/CMR - acquire_globsnow_swe() and
  # acquire_globsnow_se_mfsc() (defined above, alongside the other shared
  # acquisition functions) fetch directly from globsnow.info instead.
  #
  # Both functions now skip the globsnow.info directory listing entirely by
  # default whenever local files already cover the requested study window -
  # see the "Decide whether ... is even worth trying" blocks inside each
  # function. Previously they queried globsnow.info unconditionally on every
  # run, so a fully-populated local cache still paid the full DNS/HTTP
  # retry-and-backoff cost (worse for SE/MFSC, which repeats that cost once
  # per year in its archive). Set these to TRUE (env var, or edit the default
  # below) to also check globsnow.info for newly published months even when
  # local files already exist.
  GLOBSNOW_SWE_CHECK_REMOTE <- as.logical(Sys.getenv("GLOBSNOW_SWE_CHECK_REMOTE", "FALSE"))
  GLOBSNOW_SE_CHECK_REMOTE  <- as.logical(Sys.getenv("GLOBSNOW_SE_CHECK_REMOTE", "FALSE"))
  if (is.na(GLOBSNOW_SWE_CHECK_REMOTE)) GLOBSNOW_SWE_CHECK_REMOTE <- FALSE
  if (is.na(GLOBSNOW_SE_CHECK_REMOTE))  GLOBSNOW_SE_CHECK_REMOTE  <- FALSE
  cat(sprintf(
    "[OK] GLOBSNOW_SWE_CHECK_REMOTE = %s, GLOBSNOW_SE_CHECK_REMOTE = %s (whether to query globsnow.info even when local files already cover the study window)\n",
    GLOBSNOW_SWE_CHECK_REMOTE, GLOBSNOW_SE_CHECK_REMOTE
  ))
  
  globsnow_swe_files <- tryCatch(
    acquire_globsnow_swe(cache_dir = file.path(globsnow_dir, "swe"), check_remote = GLOBSNOW_SWE_CHECK_REMOTE),
    error = function(e) {
      message("[SKIP] GlobSnow SWE acquisition failed: ", conditionMessage(e))
      character(0)
    }
  )
  globsnow_se_result <- tryCatch(
    acquire_globsnow_se_mfsc(cache_dir = file.path(globsnow_dir, "se_mfsc"), check_remote = GLOBSNOW_SE_CHECK_REMOTE),
    error = function(e) {
      message("[SKIP] GlobSnow SE/MFSC acquisition failed: ", conditionMessage(e))
      NULL
    }
  )
  
  openlandmap_result <- tryCatch({
    acquire_openlandmap_climatology()
  }, error = function(e) {
    message("[SKIP] OpenLandMap climatology acquisition failed: ", conditionMessage(e))
    NULL
  })
  # Status now reflects actual per-month completeness (see the
  # incomplete_months tracking added inside acquire_openlandmap_climatology()
  # above), not just "did the function throw" - a function that returned
  # normally but left one or more months short on years (e.g. retries
  # exhausted after a Zenodo-side issue) is now reported as PARTIAL instead
  # of a blanket OK, which is what the inventory showed even when months
  # 07-09 (and potentially more) were silently incomplete.
  openlandmap_status <- if (is.null(openlandmap_result)) {
    "SKIP - see log"
  } else if (length(openlandmap_result$incomplete_months) > 0) {
    sprintf("PARTIAL - %d/12 month(s) incomplete, see log", length(openlandmap_result$incomplete_months))
  } else {
    "OK"
  }
  
  # C-SNOW status: reflects whether the local ingestion actually produced
  # rows, distinct from "the function was never even asked to run." Now
  # resumable across runs (see acquire_csnow_local()'s RESUMABILITY sections),
  # so this reports the cumulative total on disk (across every run to date)
  # alongside what THIS run specifically did - newly processed, skipped
  # because already done previously, and failed to read this run.
  csnow_status <- if (is.null(csnow_result)) {
    "SKIP - see log"
  } else if (isTRUE(nrow(csnow_result$basin) > 0 || nrow(csnow_result$station) > 0)) {
    sprintf(
      "OK - %d date(s) ingested total (all runs); this run: %d newly processed, %d skipped (already done), %d failed to read",
      nrow(csnow_result$basin), csnow_result$n_processed_this_run,
      csnow_result$n_skipped_already_done, csnow_result$n_read_failed
    )
  } else {
    "SKIP - see log"
  }
  sentinel1_raw_status <- if (!RUN_SENTINEL1_RAW_ASF) {
    "SKIPPED BY DESIGN - using local C-SNOW instead (set RUN_SENTINEL1_RAW_ASF=TRUE to re-enable)"
  } else if (length(sentinel1_files) > 0) {
    "OK"
  } else {
    "SKIP - see log"
  }
  
  # Final inventory - status column reflects what actually succeeded this run
  # (previously a mid-pipeline stop() meant this file was never written at
  # all on a failed run, so there was no record of what got skipped and why;
  # now every section reports its own status even if an earlier one failed).
  inv <- data.frame(
    source = c("ASWS+CanSWE daily.csv", "Sentinel-1 ASF (raw)", "RCM manifest/local", "C-SNOW (local NetCDF)",
               "CMC monthly SWE", "GlobSnow", "CanSISE", "OpenLandMap"),
    path = c(DAILY_CSV_PATH, "snow_data/sar/sentinel1", "snow_data/sar/rcm", "snow_data/csnow",
             cmc_dir, globsnow_dir, cansise_dir, "snow_data/openlandmap_snow_climatology"),
    status = c("OK",
               sentinel1_raw_status,
               "n/a (manifest/local only)",
               csnow_status,
               ifelse(cmr_ok["NSIDC_0447"], "OK", "SKIP - see log"),
               "SKIP - not archived at NSIDC, see message above",
               ifelse(cmr_ok["NSIDC_0668"], "OK", "SKIP (optional - not required by H8/Bridge) - see log"),
               openlandmap_status),
    stringsAsFactors = FALSE
  )
  write.csv(inv, file.path("snow_data", "revised_acquisition_inventory.csv"), row.names = FALSE)
  cat(sprintf(
    "[OK] Part B (revised data acquisition) complete. Acquisition inventory written to %s.\n",
    normalizePath(file.path("snow_data", "revised_acquisition_inventory.csv"))
  ))
  cat("[NEXT STEP] Run Nechako_Snow_Bridge.R next (before H8SWEI_Snow.R) to build sar_depth_monthly.nc\n")
  cat("            and the cross-validation basin-mean CSVs from what was just downloaded.\n")
  message("\n[OK] Revised acquisition pipeline complete.")
}