# ==============================================================================
# Nechako_Snow_Bridge.R
# ==============================================================================
# RUN ORDER:  Nechako_Snow_Data_download.R  ->  Nechako_Snow_Bridge.R (this)  ->  H8SWEI_Snow.R
#
# WHY THIS EXISTS
# ------------------------------------------------------------------------------
# H8SWEI_Snow.R and Nechako_Snow_Data_download.R were audited output-by-output
# against each other. Most of the hand-off works correctly as-is:
#
#   PRODUCED BY DOWNLOAD SCRIPT                          CONSUMED BY H8 AS
#   daily.csv (ASWS)                                      STATION_OBS_PATH        [OK]
#   snow_data/canswe/CanSWE_daily.csv                     CANSWE_STATION_OBS_PATH [OK]
#   Spatial/nechako_dem_mosaic.tif                        DEM_PATH                [OK]
#   snow_data/openlandmap_snow_climatology/*.tif          SCF_CLIMATOLOGY_DIR     [OK]
#
# But two things H8 is written to consume were never actually produced by the
# download script, even though comments in BOTH scripts refer to them as if
# they exist ("Part D", "Part C", "SEASON_MONTHS"):
#
#   H8 EXPECTS (via ...)          WHY IT'S MISSING FROM THE DOWNLOAD SCRIPT
#   sar_depth_monthly.nc          Only Section named "C-SNOW" exists, and it
#   (SAR_DEPTH_PATH)               writes daily station/basin CSVs, not a
#                                  monthly gridded NetCDF stack. Nothing
#                                  converts the raw daily snd_YYYYMMDD.nc
#                                  tiles it already downloaded into that form.
#
#   snow_data/cross_validation/    "Part D" is referenced in comments in both
#   *_monthly_basinmean.csv        scripts but does not exist as code anywhere.
#   (CROSS_VALIDATION_DIR)         "Part C" (Crocus-ERA5) likewise doesn't
#                                  exist despite ncdf4/ggplot2 being loaded
#                                  for it.
#
# Neither gap is fatal - H8 degrades gracefully and just logs a fallback/skip
# for each - but it means two of the accuracy-enhancing checks you likely
# intend to run are currently silent no-ops, and the download script's own
# C-SNOW section literally tells you to do a comparison it doesn't do itself:
#   "[NEXT STEP] Join this against daily.csv on station_id/date to see how
#    well C-SNOW agrees with your ground-truth SWE/depth record BEFORE
#    relying on it further - that comparison is not done automatically by
#    this script."
#
# This script closes those gaps using ONLY data the download script already
# put on disk - no new external downloads, no Earthdata/CDSE credentials
# needed here. Specifically it:
#
#   1. Converts the C-SNOW daily Sentinel-1 snow-depth NetCDFs
#      (CSNOW_LOCAL_DIR, e.g. snow_data/csnow or wherever you pointed it) into
#      sar_depth_monthly.nc - a monthly-mean gridded stack on the ERA5-Land
#      grid, which is the exact file/format H8's SAR_DEPTH_PATH expects.
#
#   2. Builds snow_data/cross_validation/*_monthly_basinmean.csv from every
#      independent gridded/point product actually available on disk:
#        - cmc_swe_monthly_basinmean.csv         (from snow_data/cmc, extracted from the downloaded zip)
#        - csnow_sar_monthly_basinmean.csv       (from the C-SNOW ingestion,
#                                                  reusing CSNOW_basin_summary.csv)
#      AU_DySno, MODIS-VIIRS, HLS and Crocus-ERA5 are NOT produced here - none
#      of those were actually downloaded by the acquisition script, so
#      fabricating basin-mean CSVs for them would be worse than leaving H8's
#      existing "skip" behaviour alone. If you later add a real downloader
#      for any of those, just have it write another
#      snow_data/cross_validation/<name>_monthly_basinmean.csv with columns
#      (date, <name>) and H8 will pick it up with no changes on its side.
#
#   3. Directly performs the C-SNOW-vs-ground-truth comparison the download
#      script's own message says isn't done automatically, using
#      CSNOW_station_extract.csv against daily.csv.
#
#   4. Re-checks the ASWS/CanSWE dedup assumption described above: confirms
#      that daily.csv's CanSWE-sourced rows really are exact duplicates of
#      CanSWE_daily.csv's rows once both pass through H8's own
#      standardize_station_obs() logic (which is what makes H8's unique()
#      dedup work). This currently holds ONLY because both files happen to
#      use the same "station_id" column name - it is not guaranteed by
#      design, so this is re-verified here rather than assumed.
#
# WHAT THIS SCRIPT DELIBERATELY DOES NOT DO
# ------------------------------------------------------------------------------
# - It does not re-download anything. Every input here must already exist
#   from a completed run of Nechako_Snow_Data_download.R.
# - It does not touch daily.csv, CanSWE_daily.csv, the DEM, or the OpenLandMap
#   climatology rasters - those hand-offs already work and are left alone.
# - It does not invent Crocus-ERA5, AU_DySno, MODIS-VIIRS, or HLS products.
#   Building real downloaders for those is future work, not something that
#   can be faked from what's already on disk.
# ==============================================================================

suppressWarnings({
  old_warn_opt <- getOption("warn")
  if (!is.null(old_warn_opt) && old_warn_opt >= 2) {
    message("[OK] Downgrading options(warn = 2) -> options(warn = 1) for this script's run (same reasoning as Nechako_Snow_Data_download.R).")
    options(warn = 1)
  }
})

library(terra)
library(ncdf4)
library(dplyr)
library(readr)
library(lubridate)

# ==============================================================================
# WORKING DIRECTORY
# ==============================================================================
# MUST be the exact same folder Nechako_Snow_Data_download.R and H8SWEI_Snow.R
# both use - this is where daily.csv, CanSWE_daily.csv, the DEM, and
# snow_data/ all already live from the download run.
setwd("D:/Nechako_Drought/Nechako")

# ==============================================================================
# CONFIGURATION - mirror the equivalent settings in the other two scripts
# ==============================================================================
DEM_PATH               <- "Spatial/nechako_dem_mosaic.tif"      # same as both scripts
STATION_OBS_PATH       <- "snow_data/daily.csv"                            # same as H8SWEI_Snow.R
CANSWE_STATION_OBS_PATH <- "snow_data/canswe/CanSWE_daily.csv"   # NOTE: see path check below
SAR_DEPTH_PATH_OUT     <- "snow_data/sar_depth_monthly.nc"                 # same as H8's SAR_DEPTH_PATH
CROSS_VALIDATION_DIR   <- "snow_data/cross_validation"           # same as H8's CROSS_VALIDATION_DIR
# Same file the download script already treats as the true ERA5-Land SWE
# series/study-window source (see its "Study window from local ERA5-Land
# SWE" section) - used below as the actual "ERA5-Land grid" template for
# sar_depth_monthly.nc, rather than the DEM standing in for it.
ERA5_LAND_SWE_FILE     <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
# Safety cap on how many cells an alignment template is allowed to have
# before Step 2 will coarsen it first. The Nechako DEM is ~490 million cells
# at ~30m resolution (built for per-station elevation/aspect screening in
# Part A) - resampling C-SNOW tiles onto that grid and stacking a month's
# worth of them is what crashed a prior run with "std::bad_alloc". 2 million
# cells is already generous for a basin-mean monthly comparison.
SAR_DEPTH_TEMPLATE_MAX_CELLS <- 2e6
CSNOW_DIR              <- file.path("snow_data", "csnow")        # written by the download script
CMC_DIR                <- file.path("snow_data", "cmc")          # written by the download script

# BUGFIX: this was referenced (STEP 1's mismatch branch) but never actually
# defined anywhere in the script - if verify_canswe_dedup() ever found a
# mismatch, that branch threw "object 'STOP_ON_CANSWE_MISMATCH' not found"
# instead of doing either documented behaviour (stop, or continue anyway).
# Defaults to TRUE (fail closed) to match the "refusing to proceed" wording
# already in that stop() message; set to FALSE to log the mismatch and
# continue instead.
STOP_ON_CANSWE_MISMATCH <- as.logical(Sys.getenv("STOP_ON_CANSWE_MISMATCH", "TRUE"))
if (is.na(STOP_ON_CANSWE_MISMATCH)) STOP_ON_CANSWE_MISMATCH <- TRUE

# Plausible bulk-density range (kg/m3) for the C-SNOW-vs-station comparison
# in compare_csnow_to_ground_truth() below. Used two ways: (1) as a hard
# GATE - any station whose median implied density falls outside this range
# is flagged and excluded from the correlation-only comparison entirely,
# rather than just noted in a log line a person has to go read; and (2) as
# a REGIONAL PRIOR - stations that pass the gate contribute to a single
# basin-wide median density, which is then used to convert C-SNOW depth (m)
# into an approximate SWE (mm) for an actual mm-bias number, clearly labeled
# as density-dependent/uncertain rather than a real SWE-vs-SWE validation.
# 100-400 kg/m3 spans fresh dry snow (~50-150) through settled/wet spring
# snowpack (~300-450) for a mid-latitude seasonal snowpack; values far
# outside this are far more likely to reflect a geolocation/timing mismatch
# or a masked C-SNOW pixel at that station than a real snow density.
CSNOW_PLAUSIBLE_DENSITY_MIN_KGM3 <- 100
CSNOW_PLAUSIBLE_DENSITY_MAX_KGM3 <- 400

# Minimum retained (gate-passing) stations contributing to a SEASON's median
# implied density before that season's prior is trusted at full weight (see
# compare_csnow_to_ground_truth() below). Below this, the season's density is
# empirical-Bayes shrunk toward the basin-wide, all-seasons prior, weighted by
# n_season_stations / this value (capped at 1) - same shrinkage style already
# used for sparse-month Tier 1 SWE bias correction in H8SWEI_Snow.R's
# fit_basin_bias_correction(). Snow density varies systematically over a
# season (fresh/dry early winter, settled/wet spring) - a single basin-wide
# scalar density can't capture that, but splitting by season with no
# safeguard risks noisy medians from 1-2 stations. This gets the seasonal
# signal without the small-sample instability.
CSNOW_MIN_STATIONS_PER_SEASON <- 3

BRIDGE_DIAG_DIR <- file.path("snow_data", "bridge_diagnostics")
if (!dir.exists(BRIDGE_DIAG_DIR)) dir.create(BRIDGE_DIAG_DIR, recursive = TRUE)
BRIDGE_LOG <- character(0)
log_msg <- function(...) {
  m <- sprintf(...)
  cat(m, "\n")
  BRIDGE_LOG <<- c(BRIDGE_LOG, m)
}

log_msg("============================================================")
log_msg("NECHAKO SNOW BRIDGE - filling gaps between download script and H8SWEI_Snow.R")
log_msg("============================================================")

# ------------------------------------------------------------------------------
# PATH CHECK: CanSWE_daily.csv location
# ------------------------------------------------------------------------------
# H8SWEI_Snow.R's own CANSWE_STATION_OBS_PATH default is "CanSWE_daily.csv" in
# the working directory root, but the download script's acquire_canswe()
# writes it to file.path(out_dir, "CanSWE_daily.csv") where out_dir defaults
# to "snow_data/canswe" - i.e. "snow_data/canswe/CanSWE_daily.csv", NOT the
# root. If you haven't already changed CANSWE_STATION_OBS_PATH inside
# H8SWEI_Snow.R to match, H8 silently proceeds with ASWS-only bias correction
# (BIAS_CORRECT still runs off daily.csv alone) and CanSWE never actually gets
# used, without erroring - exactly the kind of silent degrade this bridge is
# meant to catch.
canswe_root    <- "CanSWE_daily.csv"
canswe_nested  <- "snow_data/canswe/CanSWE_daily.csv"
if (file.exists(canswe_nested) && !file.exists(canswe_root)) {
  log_msg("[CHECK] CanSWE_daily.csv found at '%s' but NOT at the root ('%s').", canswe_nested, canswe_root)
  log_msg("        H8SWEI_Snow.R's default CANSWE_STATION_OBS_PATH is \"CanSWE_daily.csv\" (root).")
  log_msg("        Fix ONE of the following before running H8SWEI_Snow.R, or CanSWE will be silently skipped:")
  log_msg("          (a) In H8SWEI_Snow.R, change CANSWE_STATION_OBS_PATH to \"%s\", or", canswe_nested)
  log_msg("          (b) copy/symlink the file to the root: file.copy(\"%s\", \"%s\")", canswe_nested, canswe_root)
  # Do it here so a bridge run + H8 run in sequence "just works" without manual steps,
  # while still logging loudly above so it's not a silent fix.
  file.copy(canswe_nested, canswe_root, overwrite = FALSE)
  log_msg("        [FIXED THIS RUN] Copied '%s' -> '%s'.", canswe_nested, canswe_root)
} else if (file.exists(canswe_root)) {
  log_msg("[OK] CanSWE_daily.csv already present at the root - no path fix needed.")
} else {
  log_msg("[WARN] CanSWE_daily.csv not found at either the root or '%s'.", canswe_nested)
  log_msg("       CanSWE densification will be unavailable to H8's bias correction this run.")
}

# ==============================================================================
# STEP 1 - VERIFY THE daily.csv / CanSWE_daily.csv DEDUP ASSUMPTION
# ==============================================================================
# Nechako_Snow_Data_download.R Part B merges CanSWE INTO daily.csv itself
# (so daily.csv already contains CanSWE_* rows) AND writes CanSWE_daily.csv
# separately. H8SWEI_Snow.R then reads BOTH files and rbind()s + unique()s
# them in merge_station_networks(). That only avoids double-counting CanSWE
# because both reads happen to produce identical station_id/lon/lat/date/
# swe_mm tuples for CanSWE rows - which holds today only because both files
# use the literal column name "station_id" and identical values reach it.
# This is fragile (works by naming coincidence, not by design), so verify it
# explicitly rather than assuming it silently continues to hold.
log_msg("\n----- STEP 1: Verifying CanSWE double-merge is actually being deduplicated -----")

verify_canswe_dedup <- function(daily_path = "snow_data/daily.csv", canswe_path = "snow_data/canswe/CanSWE_daily.csv") {
  if (!file.exists(daily_path) || !file.exists(canswe_path)) {
    log_msg("[SKIP] Need both %s and %s to verify dedup - one or both missing.", daily_path, canswe_path)
    return(invisible(NULL))
  }
  daily <- tryCatch(read_csv(daily_path, show_col_types = FALSE), error = function(e) NULL)
  canswe <- tryCatch(read_csv(canswe_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(daily) || is.null(canswe)) {
    log_msg("[SKIP] Could not read one of the two files.")
    return(invisible(NULL))
  }
  
  daily_canswe_rows <- daily %>% filter(grepl("^CanSWE_", station_id))
  if (nrow(daily_canswe_rows) == 0) {
    log_msg("[OK] daily.csv contains no CanSWE_* rows (CanSWE merge into daily.csv did not happen, or produced 0 rows) - nothing to dedupe, no double-count risk this run.")
    return(invisible(NULL))
  }
  
  n_daily_canswe <- nrow(daily_canswe_rows)
  n_canswe_file  <- nrow(canswe)
  
  # Reproduce H8's own standardize_station_obs()/unique() path at a small
  # scale: same station_id/lon/lat/date/swe_mm tuple test it relies on.
  key_cols <- c("station_id", "lon", "lat", "date", "swe_mm")
  daily_keys  <- daily_canswe_rows %>% mutate(date = as.character(as.Date(date))) %>% select(all_of(key_cols)) %>% distinct()
  canswe_keys <- canswe            %>% mutate(date = as.character(as.Date(date))) %>% select(all_of(key_cols)) %>% distinct()
  
  matched <- dplyr::intersect(daily_keys, canswe_keys)
  n_matched <- nrow(matched)
  n_daily_only  <- nrow(dplyr::setdiff(daily_keys, canswe_keys))
  n_canswe_only <- nrow(dplyr::setdiff(canswe_keys, daily_keys))
  
  log_msg("  daily.csv CanSWE-sourced rows:      %d (%d distinct keys)", n_daily_canswe, nrow(daily_keys))
  log_msg("  CanSWE_daily.csv rows:              %d (%d distinct keys)", n_canswe_file, nrow(canswe_keys))
  log_msg("  Exact-match keys (will be deduped):  %d", n_matched)
  log_msg("  Keys only in daily.csv's copy:        %d", n_daily_only)
  log_msg("  Keys only in CanSWE_daily.csv's copy: %d", n_canswe_only)
  
  if (n_daily_only == 0 && n_canswe_only == 0) {
    log_msg("[OK] Every CanSWE row in daily.csv exactly matches a row in CanSWE_daily.csv - H8's unique() in merge_station_networks() WILL correctly collapse them to one copy each. No double-counting.")
  } else {
    log_msg("[!!] MISMATCH: %d row(s) in daily.csv and %d row(s) in CanSWE_daily.csv do NOT have an exact match on the other side.", n_daily_only, n_canswe_only)
    log_msg("     This means H8's unique()-based dedup will NOT fully collapse these - some CanSWE")
    log_msg("     observations will be counted TWICE in bias correction / RFA pooling (once via")
    log_msg("     daily.csv, once via CanSWE_daily.csv), silently inflating the influence of CanSWE")
    log_msg("     stations relative to true-independent ASWS stations. Likely causes: re-running Part")
    log_msg("     A/B at different times so the two files no longer reflect the same CanSWE vintage,")
    log_msg("     or floating-point formatting drift in swe_mm/lon/lat across the two write_csv() calls.")
    log_msg("     FIX: re-run Nechako_Snow_Data_download.R Part A+B together in one session so both")
    log_msg("     files are written from the same in-memory canswe_daily object, or remove the CanSWE")
    log_msg("     merge-into-daily.csv step (lines ~1532-1540) and let H8 be the ONLY place merging.")
    diag_path <- file.path(BRIDGE_DIAG_DIR, "canswe_dedup_mismatch.csv")
    write_csv(dplyr::bind_rows(
      dplyr::setdiff(daily_keys, canswe_keys)  %>% mutate(only_in = "daily.csv"),
      dplyr::setdiff(canswe_keys, daily_keys)  %>% mutate(only_in = "CanSWE_daily.csv")
    ), diag_path)
    log_msg("     Mismatched rows written to %s for inspection.", diag_path)
    if (n_daily_only > 0 || n_canswe_only > 0) {
      log_msg("[!!] ...")
      if (STOP_ON_CANSWE_MISMATCH) {
        stop("CanSWE dedup mismatch detected - refusing to proceed with a run that would double-count CanSWE stations. Re-run download script Part A+B together, or set STOP_ON_CANSWE_MISMATCH <- FALSE to proceed anyway.")
      }
    }
    invisible(list(n_matched = n_matched, n_daily_only = n_daily_only, n_canswe_only = n_canswe_only))
  }
}

dedup_check <- verify_canswe_dedup(daily_path = "snow_data/daily.csv", canswe_path = "snow_data/canswe/CanSWE_daily.csv")

# ==============================================================================
# STEP 2 - BUILD sar_depth_monthly.nc FROM C-SNOW DAILY NETCDFS
# ==============================================================================
# H8's SAR_DEPTH_PATH ("sar_depth_monthly.nc") is read as an "optional monthly
# SAR-depth / change-detection stack". The download script's C-SNOW ingestion
# only reads the raw daily snd_YYYYMMDD.nc tiles into point/basin-mean CSVs -
# the tiles themselves are still sitting in CSNOW_LOCAL_DIR. This step
# reprojects/resamples each daily tile onto the SAME grid as the DEM (a stand-
# in for "the ERA5-Land grid" - swap in an actual ERA5-Land raster path below
# if you have one handy; the DEM is used here only as a template because
# H8SWEI_Snow.R already guarantees it's on the Nechako extent/resolution the
# rest of the pipeline works with) and writes a monthly-mean stack.
log_msg("\n----- STEP 2: Building sar_depth_monthly.nc from C-SNOW daily tiles -----")

# BUGFIX: this used to default to CSNOW_DIR ("snow_data/csnow") when the
# CSNOW_LOCAL_DIR environment variable wasn't set - but CSNOW_DIR is where
# Nechako_Snow_Data_download.R writes its OWN derived CSVs (CSNOW_basin_
# summary.csv etc.), not where the raw snd_YYYYMMDD.nc tiles live. The
# download script's own default for the raw tiles is
# "snow_data/C-SNOW/West_US_Canada" (confirmed against its logged
# "CSNOW_LOCAL_DIR = ..." line and a direct folder listing) - but that
# default only ever lived inside the OTHER script, so this script's fallback
# never matched it unless CSNOW_LOCAL_DIR happened to also be set as a real
# OS environment variable (it wasn't). The result: this step silently
# reported "no tiles found" every run even with 1000+ real tiles on disk,
# because it was only ever looking in the CSV folder.
CSNOW_LOCAL_DIR <- Sys.getenv("CSNOW_LOCAL_DIR", unset = "snow_data/C-SNOW/West_US_Canada")
# Prefer the RAW source directory (the actual snd_YYYYMMDD.nc tiles) over the
# CSV-output directory - CSNOW_LOCAL_DIR (now correctly defaulted above)
# points at the raw NetCDF source, while CSNOW_DIR ("snow_data/csnow") is
# where the download script writes ITS OWN derived CSVs. Check both,
# preferring whichever actually has the raw .nc tiles.
find_csnow_source_dir <- function(candidates, search_root = "snow_data") {
  for (d in candidates) {
    if (dir.exists(d)) {
      hits <- list.files(d, pattern = "^snd_[0-9]{8}\\.nc$", full.names = TRUE, recursive = TRUE)
      if (length(hits) > 0) return(d)
    }
  }
  # LAST-RESORT FALLBACK: none of the fixed candidate paths panned out (e.g.
  # the folder layout changes again in the future) - rather than silently
  # giving up the way the old default did, search recursively under
  # search_root for ANY folder containing snd_YYYYMMDD.nc files and use the
  # first one found. This is a genuine search, not another guessed constant,
  # so it survives future path changes without needing another manual fix.
  if (dir.exists(search_root)) {
    all_hits <- list.files(search_root, pattern = "^snd_[0-9]{8}\\.nc$", full.names = TRUE, recursive = TRUE)
    if (length(all_hits) > 0) {
      found_dir <- dirname(all_hits[1])
      log_msg("[INFO] C-SNOW tiles not found at any expected path - located via recursive search under '%s' instead: '%s' (%d file(s)).",
              search_root, found_dir, length(all_hits))
      return(found_dir)
    }
  }
  NA_character_
}
csnow_src_dir <- find_csnow_source_dir(unique(c(CSNOW_LOCAL_DIR, CSNOW_DIR, "snow_data/csnow")))

build_sar_depth_monthly <- function(source_dir, template_path, out_path,
                                    era5_land_path = ERA5_LAND_SWE_FILE,
                                    max_template_cells = SAR_DEPTH_TEMPLATE_MAX_CELLS) {
  if (is.na(source_dir)) {
    log_msg("[SKIP] No 'snd_YYYYMMDD.nc' tiles found under any of the checked C-SNOW directories.")
    log_msg("       %s will not be (re)built this run; H8 will fall back to non-SAR predictors.", out_path)
    return(invisible(NULL))
  }
  
  # BUGFIX: this used to be `template <- rast(template_path)` where
  # template_path was always the DEM (~490 million cells at ~30m
  # resolution - built for per-station elevation/aspect screening in Part A,
  # not for use as a resample target). terra::resample(day_tile, template,
  # ...) makes an OUTPUT the size of `template`, so every single daily
  # C-SNOW tile (originally only ~1-2MB) was being blown up into a
  # ~490-million-cell raster; stacking ~30 of those for one month's mean is
  # what crashed a prior run with "Error: std::bad_alloc". This now prefers
  # the actual ERA5-Land grid - which is what this file is documented to
  # align to in the first place, and is naturally coarse (a handful of
  # thousand cells for the basin, not hundreds of millions) - and only falls
  # back to the DEM, coarsened to a safe cell count, if the ERA5-Land file
  # isn't available.
  get_template <- function() {
    if (file.exists(era5_land_path)) {
      r <- tryCatch(terra::rast(era5_land_path), error = function(e) NULL)
      if (!is.null(r)) {
        if (terra::nlyr(r) > 1) r <- r[[1]]
        log_msg("  [OK] Using the ERA5-Land grid ('%s', %.0f cell(s)) as the SAR-depth alignment template.",
                era5_land_path, terra::ncell(r))
        return(r)
      }
      log_msg("  [WARN] Found '%s' but could not open it with terra - falling back to the DEM.", era5_land_path)
    } else {
      log_msg("  [INFO] ERA5-Land file '%s' not found - falling back to the DEM as an alignment template.", era5_land_path)
    }
    
    if (!file.exists(template_path)) {
      log_msg("[SKIP] Fallback template '%s' (DEM) also not found - cannot build %s.", template_path, out_path)
      return(NULL)
    }
    dem_r <- tryCatch(terra::rast(template_path), error = function(e) NULL)
    if (is.null(dem_r)) {
      log_msg("[SKIP] Could not open fallback template '%s' - cannot build %s.", template_path, out_path)
      return(NULL)
    }
    if (terra::nlyr(dem_r) > 1) dem_r <- dem_r[[1]]
    
    n_cells <- terra::ncell(dem_r)
    if (n_cells > max_template_cells) {
      fact <- ceiling(sqrt(n_cells / max_template_cells))
      log_msg("  [WARN] DEM fallback template has %.0f cells - this is what previously caused std::bad_alloc when used as a resample target. Aggregating by factor %d (~%.0f cells after) before use.",
              n_cells, fact, n_cells / fact^2)
      dem_r <- tryCatch(terra::aggregate(dem_r, fact = fact, fun = "mean", na.rm = TRUE), error = function(e) dem_r)
    }
    dem_r
  }
  
  template <- get_template()
  if (is.null(template)) return(invisible(NULL))
  
  nc_files <- list.files(source_dir, pattern = "^snd_[0-9]{8}\\.nc$", full.names = TRUE, recursive = TRUE)
  file_dates <- as.Date(sub("^snd_([0-9]{8})\\.nc$", "\\1", basename(nc_files)), format = "%Y%m%d")
  ok <- !is.na(file_dates)
  nc_files <- nc_files[ok]; file_dates <- file_dates[ok]
  if (length(nc_files) == 0) {
    log_msg("[SKIP] No validly-dated C-SNOW tiles found in '%s'.", source_dir)
    return(invisible(NULL))
  }
  log_msg("  Found %d C-SNOW daily tile(s) in '%s' spanning %s to %s.",
          length(nc_files), source_dir, min(file_dates), max(file_dates))
  
  # Detect the data variable the same way the download script does, from the
  # first file, then assume consistency across the archive (mirrors
  # detect_csnow_varname()'s own documented assumption).
  detect_varname <- function(nc_path) {
    nc <- tryCatch(ncdf4::nc_open(nc_path), error = function(e) NULL)
    if (is.null(nc)) return(NA_character_)
    on.exit(ncdf4::nc_close(nc))
    vn <- names(nc$var)
    # prefer an obviously snow-depth-like name, else just take the first
    # non-dimension variable, mirroring the download script's own fallback
    # behaviour (terra::rast(f, subds=varname) -> terra::rast(f) fallback).
    hit <- vn[grepl("snd|depth|snow", vn, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
    if (length(vn) > 0) return(vn[1])
    NA_character_
  }
  varname <- detect_varname(nc_files[1])
  if (is.na(varname)) {
    log_msg("[SKIP] Could not detect a data variable in '%s' - aborting SAR-depth build.", basename(nc_files[1]))
    return(invisible(NULL))
  }
  log_msg("  Using NetCDF variable '%s' as C-SNOW snow depth.", varname)
  
  yearmon <- format(file_dates, "%Y-%m")
  months_present <- sort(unique(yearmon))
  log_msg("  Aggregating into %d calendar-month(s) of monthly means.", length(months_present))
  
  monthly_layers <- list()
  n_failed <- 0
  n_total_tiles <- length(nc_files)
  n_processed <- 0
  t_start <- Sys.time()
  PROGRESS_EVERY <- 25  # tiles between progress lines - was 0 (silent) before, which is what made a slow run look hung
  
  for (ym in months_present) {
    idx <- which(yearmon == ym)
    day_rasters <- list()
    for (i in idx) {
      # BUGFIX: this used to try `subds = varname` first on every single
      # tile. detect_varname() already confirmed these are classic
      # (single-variable) NetCDFs, not multidim/group files, so passing
      # `subds` routed every one of 1000+ reads through terra's multidim API,
      # which doesn't apply here and issued "[rast] cannot open this file
      # with the multidim API" as a WARNING (not an error) on every tile - a
      # real run logged 50+ of these. tryCatch(..., error=) can't catch a
      # warning, so `r` still ended up populated via terra's own internal
      # fallback and the run completed correctly, but noisily. Try the plain
      # classic read first (fast, no warning for this file type) and only
      # fall back to `subds` for files that genuinely need it.
      r <- tryCatch(terra::rast(nc_files[i]), error = function(e) NULL)
      if (is.null(r)) r <- tryCatch(terra::rast(nc_files[i], subds = varname), error = function(e) NULL)
      if (is.null(r)) { n_failed <- n_failed + 1; next }
      if (terra::nlyr(r) > 1) r <- r[[1]]
      
      # PERFORMANCE FIX: terra::project() rebuilds a full coordinate-
      # transform pipeline and rewarps the whole tile - much more expensive
      # than resample() alone. It was previously called on EVERY one of up
      # to ~1000+ daily tiles even when the tile is already in the
      # template's CRS (the common case here: C-SNOW and ERA5-Land/DEM are
      # both geographic WGS84), making it pure overhead. Only reproject
      # when the CRS genuinely differs; resample() alone handles the
      # grid/resolution alignment either way.
      r_al <- tryCatch({
        if (isTRUE(terra::same.crs(r, template))) {
          # If CRS matches, we only need to align the grid geometries
          terra::resample(r, template, method = "bilinear")
        } else {
          # If CRS differs, warp straight to the final template grid in one step
          terra::project(r, template, method = "bilinear")
        }
      }, error = function(e) NULL)
      if (!is.null(r_al)) day_rasters[[length(day_rasters) + 1]] <- r_al
      
      # PROGRESS: this loop can run for a long time on a full C-SNOW
      # archive (1000+ tiles, each individually reprojected/resampled) -
      # previously nothing was printed between "Aggregating into N
      # month(s)..." and the very end, so a slow-but-working run was
      # indistinguishable from a hung one. Report periodically instead.
      n_processed <- n_processed + 1
      if (n_processed %% PROGRESS_EVERY == 0 || n_processed == n_total_tiles) {
        elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
        rate <- if (elapsed_sec > 0) n_processed / elapsed_sec else NA_real_
        eta_sec <- if (!is.na(rate) && rate > 0) (n_total_tiles - n_processed) / rate else NA_real_
        log_msg("  ...processed %d/%d C-SNOW tile(s) (%.2f tiles/sec, elapsed %.0fs, ETA ~%s).",
                n_processed, n_total_tiles,
                if (is.na(rate)) 0 else rate, elapsed_sec,
                if (is.na(eta_sec)) "unknown" else sprintf("%.0fs", eta_sec))
        flush.console()
      }
    }
    if (length(day_rasters) == 0) {
      log_msg("  [WARN] %s: no readable/alignable tiles - month skipped.", ym)
      next
    }
    stacked <- terra::rast(day_rasters)
    month_mean <- terra::app(stacked, fun = mean, na.rm = TRUE)
    names(month_mean) <- ym
    monthly_layers[[ym]] <- month_mean
  }
  
  if (n_failed > 0) {
    log_msg("  [WARN] %d individual daily tile(s) failed to read/align and were excluded from their month's mean.", n_failed)
  }
  if (length(monthly_layers) == 0) {
    log_msg("[SKIP] No month produced a usable mean - %s not written.", out_path)
    return(invisible(NULL))
  }
  
  full_stack <- terra::rast(monthly_layers)
  time_vals <- as.Date(paste0(names(monthly_layers), "-01"))
  terra::time(full_stack) <- time_vals
  
  terra::writeCDF(full_stack, out_path, varname = "sar_snow_depth_m", unit = "m",
                  longname = "C-SNOW (Lievens Sentinel-1) monthly mean snow depth, aligned to Nechako grid",
                  overwrite = TRUE)
  log_msg("[OK] Wrote %s: %d monthly layer(s), %s to %s.",
          out_path, terra::nlyr(full_stack), min(time_vals), max(time_vals))
  invisible(full_stack)
}

sar_depth_result <- build_sar_depth_monthly(csnow_src_dir, DEM_PATH, SAR_DEPTH_PATH_OUT)

# ==============================================================================
# STEP 3 - BUILD snow_data/cross_validation/*_monthly_basinmean.csv
# ==============================================================================
# H8's validate_against_cross_products() looks for any
# snow_data/cross_validation/*_monthly_basinmean.csv with columns (date,
# <productname>) and reports correlation + bias vs. the bias-corrected
# ERA5-Land basin mean. This builds that file for every product actually on
# disk from the download run: CMC snow depth (GeoTIFFs) and C-SNOW SAR depth
# (reusing the same monthly aggregation as Step 2, at basin-mean resolution
# instead of gridded). AU_DySno / MODIS-VIIRS / HLS / Crocus-ERA5 are NOT
# fabricated here since the download script never actually acquires them.
log_msg("\n----- STEP 3: Building cross-validation basin-mean CSVs -----")
if (!dir.exists(CROSS_VALIDATION_DIR)) dir.create(CROSS_VALIDATION_DIR, recursive = TRUE)

# --- 3a. CMC monthly SWE -> cmc_swe_monthly_basinmean.csv ---
# BUGFIX (was two separate bugs):
#  (1) The download script only ever downloads the single combined-record
#      zip "cmc_swe_mly_<startyr>to<endyr>_v0X.X.zip" (see cmc_swe_pattern in
#      Nechako_Snow_Data_download.R) and never extracts it - a folder listing
#      confirms cmc_dir holds only .zip files, so the old "list *.tif" scan
#      here always found 0 files and silently skipped every run. This now
#      extracts the matching zip itself, once, before scanning for data files.
#  (2) The old code also assumed the wrong internal structure entirely: it
#      was written for one GeoTIFF PER YEAR with one band PER DAY (a layout
#      Nechako_Snow_Data_download.R's comments describe abandoning in favour
#      of the single combined MONTHLY file). It also labelled the output
#      "cmc_snowdepth", but the file this pipeline actually acquires is SWE
#      (SWE/SWE_Monthly_Averages/), not snow depth - per the download
#      script's own documentation. That mislabelling meant H8's
#      validate_against_cross_products() would (if this had ever produced
#      output) treat it as depth and skip the mm-bias check it could
#      legitimately run for a real SWE-vs-SWE comparison (is_swe_like <-
#      grepl("swe", varname, ...)).
#
# This version: extracts the zip if needed, reads the resulting file's own
# band/time metadata to get real dates (does NOT assume the daily-per-year
# structure the old code guessed at), and outputs a "cmc_swe_mm" column so
# the mm-bias check H8 already supports actually runs.
CMC_SWE_ZIP_PATTERN <- "^cmc_swe_mly_[0-9]{4}to[0-9]{4}_v[0-9.]+\\.zip$"  # same pattern as the download script's cmc_swe_pattern - deliberately excludes the _clim_ climatology zips and the _sdepth_ (depth, not SWE) zips also seen sitting in this folder

# ------------------------------------------------------------------------
# CMC SWE flat-file format - CONFIRMED against the NSIDC-0447 User Guide
# ("Canadian Meteorological Centre (CMC) Daily Snow Depth Analysis Data,
# Version 1", section 1.5.4 "Snow Water Equivalent (SWE)"):
#   write(60,*) iyear, imonth, nconf
#   do i = 1, 706
#     write(60, '(706(1x,f6.1))') (swe(i,j), j = 1, 706)
#   end do
# repeated once per month present in the file. A real run's diagnostic
# preview showed the file's first header line as "1998 10 4631" - this is
# NOT a grid dimension (the old guess); it is iyear=1998, imonth=10 (SWE
# starts in October per the User Guide), nconf=4631 (a per-month count of
# land/sea-confusion points excluded from that month's analysis - a
# diagnostic value, not something this reader needs). The grid itself is
# always a fixed 706 x 706, not the header's third number.
# Values are -999.0 where SWE could not be calculated (ocean, Greenland, or
# undefined Sturm snow-class for that cell).
# ------------------------------------------------------------------------
CMC_GRID_N <- 706L

# Grid Specification (User Guide section 1.3.2, "Projection and Grid
# Description"): North Polar Stereographic, ni=nj=706, bottom-left grid
# point=(1,1), grid rotation -10 deg (baked into lon_0 below), pole
# position=(353,353), resolution at 60N = 23812.5 m.
CMC_CRS <- "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=10 +x_0=0 +y_0=0 +R=6371200 +units=m +no_defs=True"
CMC_RES_M <- 23812.5
CMC_POLE_I <- 353.5  # row index of the pole (documented "pole position")
CMC_POLE_J <- 353.5  # column index of the pole

# Canonical CMC matrix accessor - ALWAYS use this (not raw m[...] indexing)
# to pull a SWE value out of a block matrix produced by
# read_cmc_swe_flatfile(). Convention confirmed by CMC_INDEXING_DIAGNOSTIC_v2.R
# against 3 independent, widely-separated, off-diagonal points (Prince George
# BC, Tofino BC, Churchill MB) in a deep-winter month - the earlier
# convention m[(707-j), i] passed only a single near-diagonal test point
# (Prince George) and failed at Churchill; this one passed all three.
# Row index i is used as-is (no flip); only the column index j is flipped.
#
# IMPORTANT: i and j must be PAIRED element-wise (i[1] with j[1], i[2] with
# j[2], ...), matching the (i,j) rows of a lookup table like basin_cells.
# Uses cbind() to build a proper 2-column index matrix - m[i_vec, j_vec]
# with two vectors of length >1 would instead return the outer-product
# submatrix (every i against every j), which is NOT what basin-cell lookups
# need and would silently produce a wrong-shaped/wrong-content result.
cmc_get <- function(m, i, j, n = CMC_GRID_N) {
  m[cbind(i, (n + 1L) - j)]
}

# Cell-CENTER coordinate of grid index k along one axis, from the
# documented "bottom left grid point = (1,1)" + pole-position convention.
.cmc_center <- function(k, pole, res) (k - pole) * res

build_cmc_extent <- function(n = CMC_GRID_N, res = CMC_RES_M,
                             pole_i = CMC_POLE_I, pole_j = CMC_POLE_J) {
  x_min <- .cmc_center(1, pole_j, res) - res / 2
  x_max <- .cmc_center(n, pole_j, res) + res / 2
  y_min <- .cmc_center(1, pole_i, res) - res / 2
  y_max <- .cmc_center(n, pole_i, res) + res / 2
  terra::ext(x_min, x_max, y_min, y_max)
}

# Reads the documented format directly as a flat token stream - "3 tokens
# (iyear, imonth, nconf), then 706*706 value tokens" repeated per month, in
# file order. This sidesteps line-by-line/fixed-width parsing entirely (no
# header keywords to look for, unlike an ESRI ASCII grid) and gives EXACT
# dates straight from the file - no guessing/regex on filenames needed,
# unlike detect_band_dates() below (which stays in use for genuine rasters).
read_cmc_swe_flatfile <- function(f, n = CMC_GRID_N) {
  toks <- scan(f, what = numeric(), quiet = TRUE)
  block_size <- 3 + n * n
  n_blocks <- length(toks) %/% block_size
  leftover <- length(toks) %% block_size
  if (leftover != 0) {
    log_msg("  [WARN] '%s' has %d leftover token(s) after %d complete month-block(s) of %d tokens each - file may be truncated; only complete blocks are used.",
            basename(f), leftover, n_blocks, block_size)
  }
  if (n_blocks == 0) {
    log_msg("  [SKIP] '%s' has too few tokens for even one %dx%d month block - not the expected NSIDC-0447 SWE layout.", basename(f), n, n)
    return(NULL)
  }
  
  out <- vector("list", n_blocks)
  for (b in seq_len(n_blocks)) {
    start <- (b - 1) * block_size
    iyear  <- toks[start + 1]
    imonth <- toks[start + 2]
    vals   <- toks[(start + 4):(start + block_size)]
    # NO vertical flip. Empirically confirmed (CMC_INDEXING_DIAGNOSTIC_v2.R,
    # 3 independent widely-separated points: Prince George BC, Tofino BC,
    # Churchill MB, deep-winter month) that the matrix must be kept exactly
    # as tokens appear in the file - row index i used as-is, NOT flipped.
    # The column index j DOES need to be flipped at extraction time (see
    # cmc_get() below) - that correction happens where values are read out
    # of m, not here at read time, so this matrix is intentionally left in
    # raw file order.
    m <- matrix(vals, nrow = n, ncol = n, byrow = TRUE)
    out[[b]] <- list(date = as.Date(sprintf("%04d-%02d-01", iyear, imonth)), m = m)
  }
  dates <- do.call(c, lapply(out, `[[`, "date"))
  log_msg("  Parsed %d month-block(s) from '%s' (%s to %s), %dx%d grid, per the NSIDC-0447 documented SWE flat-file format.",
          n_blocks, basename(f), min(dates), max(dates), n, n)
  out
}

ensure_cmc_extracted <- function(cmc_dir) {
  zips <- list.files(cmc_dir, pattern = CMC_SWE_ZIP_PATTERN, full.names = TRUE)
  if (length(zips) == 0) {
    log_msg("[SKIP] No file matching '%s' found in '%s' (only the combined-record SWE zip is the real SWE product - _clim_ and _sdepth_ files are a different product/resolution).", CMC_SWE_ZIP_PATTERN, cmc_dir)
    return(character(0))
  }
  if (length(zips) > 1) {
    log_msg("[WARN] %d file(s) matched the CMC SWE pattern in '%s' - using all of them, but check this is expected (usually there should be exactly one combined-record file).", length(zips), cmc_dir)
  }
  
  extracted_dir <- file.path(cmc_dir, "_extracted")
  if (!dir.exists(extracted_dir)) dir.create(extracted_dir, recursive = TRUE)
  
  # BUGFIX: this used to only accept ".tif"/".tiff"/".nc" as the extension of
  # an extracted payload. A real run showed the zip DOES extract successfully
  # ("[EXTRACT] cmc_swe_mly_1998to2020_v01.2.zip -> ...") but every call still
  # logged "No usable CMC SWE data file after checking for/extracting the
  # zip." A direct list.files() on the extraction folder showed why:
  #   snow_data/cmc/_extracted/cmc_swe_mly_1998to2020_v01.2.txt
  # NSIDC-0447's combined-record monthly product ships as a flat text/ASCII
  # grid, not a GeoTIFF or NetCDF, so the old regex silently matched zero
  # files - extraction "succeeded" but the function still returned
  # character(0), and build_cmc_basinmean() reported a misleading skip that
  # made it look like nothing had been extracted at all. Accept the wider set
  # of extensions this product is actually known to use.
  DATA_FILE_EXT_PATTERN <- "\\.(tif|tiff|nc|txt|asc|dat)$"
  data_files <- character(0)
  for (z in zips) {
    stem <- sub("\\.zip$", "", basename(z), ignore.case = TRUE)
    already <- list.files(extracted_dir, pattern = paste0("^", stem, DATA_FILE_EXT_PATTERN),
                          full.names = TRUE, ignore.case = TRUE)
    if (length(already) > 0) {
      log_msg("  [OK] '%s' already extracted -> %s.", basename(z), paste(basename(already), collapse = ", "))
      data_files <- c(data_files, already)
      next
    }
    log_msg("  [EXTRACT] %s -> %s", basename(z), extracted_dir)
    ok <- tryCatch({ utils::unzip(z, exdir = extracted_dir); TRUE },
                   error = function(e) { log_msg("    [FAIL] Could not unzip %s: %s", basename(z), conditionMessage(e)); FALSE })
    if (!ok) next
    # Permissive on purpose: some archives rename/reformat the payload on
    # extraction, so this doesn't require the extracted name to match `stem`
    # - it just picks up whatever data file(s) landed in extracted_dir.
    newly <- list.files(extracted_dir, pattern = DATA_FILE_EXT_PATTERN, full.names = TRUE, ignore.case = TRUE)
    data_files <- unique(c(data_files, newly))
  }
  if (length(data_files) == 0 && length(zips) > 0) {
    log_msg("  [WARN] Zip(s) extracted but nothing matched '%s' in '%s' - list.files(\"%s\") to see what actually landed there and extend DATA_FILE_EXT_PATTERN if it's a format not listed yet.",
            DATA_FILE_EXT_PATTERN, extracted_dir, extracted_dir)
  }
  data_files
}

# Detect real per-band dates rather than assuming a structure. Tries, in
# order: (1) terra's own time() metadata, (2) for NetCDF sources, the raw
# time/valid_time variable via ncdf4, (3) date-like tokens in band names.
# Returns a Date vector of length nlyr(r), or NULL if none of these worked -
# callers must NOT invent dates when this returns NULL (better to skip the
# cross-validation file for this product than silently mislabel the months).
detect_band_dates <- function(r, f) {
  n <- terra::nlyr(r)
  
  t1 <- tryCatch(terra::time(r), error = function(e) NULL)
  if (!is.null(t1) && length(t1) == n && !all(is.na(t1))) {
    return(as.Date(t1))
  }
  
  if (grepl("\\.nc$", f, ignore.case = TRUE)) {
    t2 <- tryCatch({
      nc <- ncdf4::nc_open(f)
      on.exit(ncdf4::nc_close(nc), add = TRUE)
      time_var <- intersect(c("time", "valid_time", "Time"), names(nc$var))
      if (length(time_var) == 0 && "time" %in% names(nc$dim)) time_var <- "time"
      if (length(time_var) == 0) return(NULL)
      tv <- ncdf4::ncvar_get(nc, time_var[1])
      units_str <- tryCatch(ncdf4::ncatt_get(nc, time_var[1], "units")$value, error = function(e) NA_character_)
      origin_match <- regmatches(units_str, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", units_str))
      origin <- if (length(origin_match) > 0 && nzchar(origin_match)) as.Date(origin_match) else as.Date("1970-01-01")
      if (grepl("day", units_str, ignore.case = TRUE)) as.Date(origin + tv)
      else if (grepl("hour", units_str, ignore.case = TRUE)) as.Date(origin + tv / 24)
      else if (grepl("second", units_str, ignore.case = TRUE)) as.Date(origin + tv / 86400)
      else NULL
    }, error = function(e) NULL)
    if (!is.null(t2) && length(t2) == n) return(t2)
  }
  
  nm <- names(r)
  ym <- regmatches(nm, regexpr("(19|20)[0-9]{2}[-_.]?(0[1-9]|1[0-2])", nm))
  if (length(ym) == n && all(nzchar(ym))) {
    parsed <- as.Date(paste0(gsub("[-_.]", "", ym), "01"), format = "%Y%m%d")
    if (!any(is.na(parsed))) return(parsed)
  }
  
  NULL
}

# ==============================================================================
# CMC GEOREFERENCING - build lookup table of (i,j) -> (lon,lat)
# ==============================================================================
# Rather than deriving an extent formula from grid specs, build an explicit
# (i,j)->(x,y) lookup table once, then use points-in-polygon directly.
# This avoids silent errors from subtle off-by-one/half-cell mistakes in
# extent or pole-position formulas.
build_cmc_georef <- function(n = CMC_GRID_N, res = CMC_RES_M,
                             pole_i = CMC_POLE_I, pole_j = CMC_POLE_J,
                             crs = CMC_CRS) {
  # Build (i,j)->(x,y) lookup for every grid cell
  i_idx <- 1:n
  j_idx <- 1:n
  grid_df <- expand.grid(i = i_idx, j = j_idx)
  
  # Cell-center coordinates in the CMC polar-stereographic grid
  grid_df$x <- (grid_df$j - pole_j) * res
  grid_df$y <- (grid_df$i - pole_i) * res
  
  # Project to geographic WGS84
  pts <- terra::vect(as.matrix(grid_df[, c("x", "y")]), type = "points", crs = crs)
  pts_geo <- terra::project(pts, "EPSG:4326")
  coords_geo <- terra::geom(pts_geo)[, c("x", "y")]
  
  grid_df$lon <- coords_geo[, "x"]
  grid_df$lat <- coords_geo[, "y"]
  
  # Sanity check: the pole (353.5, 353.5) should be exactly at (0, 90N)
  pole_row <- which(grid_df$i == CMC_POLE_I & grid_df$j == CMC_POLE_J)
  if (length(pole_row) > 0) {
    pole_check <- grid_df[pole_row, ]
    # North Pole is at (any_lon, 90N). Check latitude is close to 90.
    if (!is.na(pole_check$lat) && abs(pole_check$lat - 90) > 0.5) {
      log_msg("[WARN] CMC pole position (i=%g, j=%g) projects to lat=%g instead of ~90N - georeferencing may be incorrect.",
              CMC_POLE_I, CMC_POLE_J, pole_check$lat)
      return(NULL)
    }
  }
  
  list(table = grid_df, crs_from = crs, crs_to = "EPSG:4326")
}

# Sanity-check a block's values against THREE known ground-truth points, not
# one. A single-point check (previously: Prince George only) is not enough -
# CMC_INDEXING_DIAGNOSTIC_v2.R found that a WRONG orientation (flipI_ij)
# still produced a plausible-looking positive value at Prince George alone
# (126.8 mm), and was only caught by adding Tofino (should be low/coastal)
# and, decisively, Churchill MB - a point placed off-diagonal from Prince
# George (large delta-j / longitude separation, small delta-i / latitude
# separation) specifically so that an i/j-swapped or partially-flipped
# orientation can't pass by coincidence the way it can near the diagonal.
#
# THRESHOLD NOTE: an earlier version of this check required val > 20 mm at
# Prince George and Churchill, based on a single January 1999 sample that
# happened to be a locally deep-snow month. A real run then failed this
# check with Prince George = 17.4 mm - a perfectly plausible early/moderate
# winter SWE value, not a nodata or wrong-cell artifact (Tofino and
# Churchill were correctly ordered around it: 1.1 < 17.4 < 46.9). Any given
# month/year can legitimately have thin-to-moderate snowpack, so a fixed
# magnitude floor is the wrong test. The magnitude floor here is now just
# "not near-zero/nodata" (> 1 mm) - the REAL discriminating signal is the
# relative ORDERING (coastal < interior < subarctic-interior), which only
# holds under the correct orientation and is what actually caught the bug
# during diagnosis.
#
# IMPORTANT: call this against a DEEP-WINTER block (Dec/Jan/Feb), not
# whatever month happens to be first in the file. October (this file's
# first block) is low-snow/shoulder-season and multiple wrong orientations
# can look "plausible" there too - it's the deep-winter check that actually
# discriminates, via ordering, not magnitude alone.
verify_cmc_georef <- function(georef, blk_m, blk_date = NULL) {
  
  find_cell <- function(lon, lat) {
    d <- (georef$table$lon - lon)^2 + (georef$table$lat - lat)^2
    georef$table[which.min(d), ]
  }
  is_plausible_snow <- function(val, min_ok = 1, max_ok = 2000) {
    !is.na(val) && val > min_ok && val <= max_ok
  }
  
  pg <- find_cell(-124.5, 53.0)     # Prince George, BC - interior
  to <- find_cell(-125.9, 49.15)    # Tofino, BC - coastal, should stay low even in winter
  ch <- find_cell(-94.2, 58.77)     # Churchill, MB - off-diagonal from PG, subarctic interior
  
  val_pg <- cmc_get(blk_m, pg$i, pg$j)
  val_to <- cmc_get(blk_m, to$i, to$j)
  val_ch <- cmc_get(blk_m, ch$i, ch$j)
  
  if (!is.null(blk_date) && !(lubridate::month(blk_date) %in% c(12, 1, 2))) {
    log_msg("[WARN] CMC georeferencing check is running against a non-deep-winter month (%s) - a shoulder-season month can look plausible under a WRONG orientation too. Prefer verifying against a Dec/Jan/Feb block if one is available.",
            blk_date)
  }
  
  checks <- c(
    pg_plausible = is_plausible_snow(val_pg),
    to_plausible = is_plausible_snow(val_to, min_ok = -0.01),  # coastal can be genuinely ~0
    to_lt_pg     = (!is.na(val_pg) && !is.na(val_to) && val_to < val_pg),
    ch_plausible = is_plausible_snow(val_ch)
  )
  
  log_msg("  CMC georef check - Prince George: %.1f mm | Tofino: %.1f mm | Churchill: %.1f mm",
          val_pg, val_to, val_ch)
  
  if (!all(checks)) {
    log_msg("[WARN] CMC georeferencing check FAILED: %s - this may indicate wrong indexing or an incorrect extent formula.",
            paste(names(checks)[!checks], collapse = ", "))
    return(FALSE)
  }
  
  log_msg("  CMC georeferencing verified: Prince George (%.1f mm), Tofino (%.1f mm, < PG as expected), Churchill (%.1f mm) all plausible.",
          val_pg, val_to, val_ch)
  TRUE
}

build_cmc_basinmean <- function(cmc_dir, basin_boundary_path, out_dir) {
  if (!dir.exists(cmc_dir)) {
    log_msg("[SKIP] CMC directory '%s' not found.", cmc_dir)
    return(invisible(NULL))
  }
  if (!file.exists(basin_boundary_path)) {
    log_msg("[SKIP] Basin boundary '%s' not found - cannot compute a basin mean for CMC.", basin_boundary_path)
    return(invisible(NULL))
  }
  
  data_files <- ensure_cmc_extracted(cmc_dir)
  if (length(data_files) == 0) {
    log_msg("[SKIP] No usable CMC SWE data file after checking for/extracting the zip.")
    return(invisible(NULL))
  }
  log_msg("  Found %d extracted CMC SWE data file(s).", length(data_files))
  
  basin_v <- tryCatch(terra::vect(basin_boundary_path), error = function(e) NULL)
  if (is.null(basin_v)) {
    log_msg("[SKIP] Could not read basin boundary '%s'.", basin_boundary_path)
    return(invisible(NULL))
  }
  
  # Verify that the basin vector actually has valid geometry (not a silent GDAL parse failure)
  basin_extent <- terra::ext(basin_v)
  if (!is.finite(basin_extent[1]) || !is.finite(basin_extent[2]) || 
      !is.finite(basin_extent[3]) || !is.finite(basin_extent[4])) {
    log_msg("[STOP] Basin boundary '%s' has an invalid/degenerate extent (%g, %g, %g, %g) - this indicates GDAL failed to parse the geometry (likely a KMZ truncation issue).",
            basin_boundary_path, basin_extent[1], basin_extent[2], basin_extent[3], basin_extent[4])
    log_msg("       Recommended fix: convert '%s' to a GeoJSON or Shapefile using ogr2ogr or QGIS, then point basin_boundary_path at that file instead.",
            basin_boundary_path)
    return(invisible(NULL))
  }
  
  # Build georeferencing lookup table (i,j)->(lon,lat) once
  georef <- build_cmc_georef()
  if (is.null(georef)) {
    log_msg("[SKIP] Cannot build CMC SWE basin-mean without valid georeferencing (see message above).")
    return(invisible(NULL))
  }
  
  # Find which grid cells fall inside the basin polygon using points-in-polygon test
  ref_pts <- terra::vect(as.matrix(georef$table[, c("lon", "lat")]), type = "points", crs = "EPSG:4326")
  in_basin <- terra::is.related(ref_pts, basin_v, "intersects")
  basin_cells <- georef$table[in_basin, ]
  log_msg("  %d CMC grid cell(s) fall inside the basin polygon (out of %d total).",
          nrow(basin_cells), nrow(georef$table))
  
  if (nrow(basin_cells) == 0) {
    log_msg("[SKIP] No CMC grid cells intersect the basin polygon - cannot compute a basin mean.")
    log_msg("       This would indicate a real problem (wrong basin file, or still-wrong georeferencing) -")
    log_msg("       investigate before assuming CMC coverage genuinely excludes the Nechako.")
    return(invisible(NULL))
  }
  
  # One-time sanity check, run against the first real block we process
  # below (needs actual SWE values to check against nodata).
  georef_verified <- FALSE
  
  rows <- list()
  for (f in data_files) {
    # Only process flat-file format (txt/asc/dat) - not raster GeoTIFF/NetCDF
    # which would be from a different product structure
    if (!grepl("\\.(txt|asc|dat)$", f, ignore.case = TRUE)) {
      log_msg("  [SKIP] '%s' is not a recognized CMC flat-file extension - skipping.", basename(f))
      next
    }
    
    blocks <- tryCatch(read_cmc_swe_flatfile(f), error = function(e) {
      log_msg("  [WARN] Could not parse '%s': %s", basename(f), conditionMessage(e))
      NULL
    })
    if (is.null(blocks) || length(blocks) == 0) next
    
    # Run the georeferencing sanity check against a DEEP-WINTER block
    # (Dec/Jan/Feb) if one is available in this file, rather than whichever
    # block happens to be first (often October - shoulder-season, and
    # ambiguous: a WRONG orientation can still look "plausible" there, which
    # is exactly what let an earlier bug go undetected). Falls back to the
    # first block if no deep-winter block is present in this file.
    if (!georef_verified) {
      blk_months <- sapply(blocks, function(b) lubridate::month(b$date))
      winter_idx <- which(blk_months %in% c(12, 1, 2))
      verify_blk <- if (length(winter_idx) > 0) blocks[[winter_idx[1]]] else blocks[[1]]
      ok <- verify_cmc_georef(georef, verify_blk$m, blk_date = verify_blk$date)
      if (!ok) {
        log_msg("[STOP] Refusing to compute a CMC SWE basin mean - georeferencing failed its sanity check.")
        return(invisible(NULL))
      }
      georef_verified <- TRUE
    }
    
    n_blocks_used <- 0
    for (blk in blocks) {
      # CORRECTED INDEXING: row index i used as-is, column index j flipped -
      # i.e. m[i, (n+1-j)], applied via the canonical cmc_get() accessor.
      # This superseded an earlier m[(707-j), i] convention that had passed
      # only a single near-diagonal test point (Prince George) and was later
      # found to fail at an off-diagonal point (Churchill, MB) once a
      # 3-point check was added - see CMC_INDEXING_DIAGNOSTIC_v2.R and
      # verify_cmc_georef() above for the full verification.
      vals <- cmc_get(blk$m, basin_cells$i, basin_cells$j)
      vals[vals <= -998] <- NA
      
      # Track validity: if all basin cells are nodata for this month, skip it
      n_valid <- sum(!is.na(vals))
      if (n_valid == 0) {
        log_msg("  [WARN] %s: 0/%d basin cells have valid (non-nodata) data - month skipped.",
                blk$date, nrow(basin_cells))
        next
      }
      
      mean_val <- mean(vals, na.rm = TRUE)
      rows[[length(rows) + 1]] <- data.frame(date = blk$date, cmc_swe_mm = mean_val,
                                             n_valid_cells = n_valid, n_basin_cells = nrow(basin_cells))
      n_blocks_used <- n_blocks_used + 1
    }
    
    log_msg("  [OK] Parsed '%s': %d/%d month(s) produced a usable basin-mean value.",
            basename(f), n_blocks_used, length(blocks))
  }
  
  if (length(rows) == 0) {
    log_msg("[SKIP] No CMC file produced usable dated basin-mean values.")
    return(invisible(NULL))
  }
  
  # CMC flat-file data is already MONTHLY - group only in case two source
  # files ever overlap the same month, rather than to aggregate daily data up.
  all_rows <- dplyr::bind_rows(rows)
  monthly_df <- all_rows %>%
    dplyr::filter(is.finite(cmc_swe_mm)) %>%
    dplyr::mutate(yearmon = format(date, "%Y-%m-01")) %>%
    dplyr::group_by(yearmon) %>%
    dplyr::summarise(cmc_swe_mm = mean(cmc_swe_mm, na.rm = TRUE),
                     mean_n_valid_cells = mean(n_valid_cells), .groups = "drop") %>%
    dplyr::transmute(date = as.Date(yearmon), cmc_swe_mm, mean_n_valid_cells)
  
  if (nrow(monthly_df) == 0) {
    log_msg("[SKIP] All computed basin-mean values were non-finite after filtering - not writing output.")
    return(invisible(NULL))
  }
  
  # Sanity check: SWE values should be plausible
  plausible_max <- max(monthly_df$cmc_swe_mm, na.rm = TRUE)
  if (is.finite(plausible_max) && plausible_max < 5) {
    log_msg("[WARN] CMC basin-mean values look implausibly small for mm SWE (max = %.3f) - verify units.", plausible_max)
  }
  
  out_path <- file.path(out_dir, "cmc_swe_monthly_basinmean.csv")
  write_csv(monthly_df %>% dplyr::select(date, cmc_swe_mm), out_path)
  log_msg("[OK] Wrote %s: %d month(s), %s to %s (avg %.0f/%d basin cells valid per month).",
          out_path, nrow(monthly_df), min(monthly_df$date), max(monthly_df$date),
          mean(monthly_df$mean_n_valid_cells), nrow(basin_cells))
  log_msg("     Georeferencing verified via direct NSIDC lat/lon table lookup, indexing empirically confirmed at 3 independent, off-diagonal points (Prince George BC, Tofino BC, Churchill MB) against a deep-winter block.")
  invisible(monthly_df)
}

# basin boundary: re-use the same KMZ path convention as the download script's
# Part A. Adjust here if your basin file lives somewhere else.
basin_boundary_path <- "Spatial/Nechako_Basin.kmz"
cmc_basinmean <- build_cmc_basinmean(CMC_DIR, basin_boundary_path, CROSS_VALIDATION_DIR)

# --- 3b. C-SNOW SAR depth -> csnow_sar_monthly_basinmean.csv ---
build_csnow_basinmean <- function(csnow_dir, out_dir) {
  basin_summary_path <- file.path(csnow_dir, "CSNOW_basin_summary.csv")
  if (!file.exists(basin_summary_path)) {
    log_msg("[SKIP] '%s' not found - run the download script's C-SNOW ingestion first.", basin_summary_path)
    return(invisible(NULL))
  }
  daily_df <- tryCatch(read_csv(basin_summary_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(daily_df) || !all(c("date", "csnow_basin_mean_m") %in% names(daily_df))) {
    log_msg("[SKIP] '%s' missing expected columns (date, csnow_basin_mean_m).", basin_summary_path)
    return(invisible(NULL))
  }
  monthly_df <- daily_df %>%
    dplyr::mutate(date = as.Date(date), yearmon = format(date, "%Y-%m-01")) %>%
    dplyr::filter(is.finite(csnow_basin_mean_m)) %>%
    dplyr::group_by(yearmon) %>%
    dplyr::summarise(csnow_sar_depth = mean(csnow_basin_mean_m, na.rm = TRUE), .groups = "drop") %>%
    dplyr::transmute(date = as.Date(yearmon), csnow_sar_depth)
  
  out_path <- file.path(out_dir, "csnow_sar_monthly_basinmean.csv")
  write_csv(monthly_df, out_path)
  log_msg("[OK] Wrote %s: %d month(s). Units are metres of SAR-derived snow depth (also correlation-only vs. ERA5-Land SWE, not mm-bias).",
          out_path, nrow(monthly_df))
  invisible(monthly_df)
}

csnow_basinmean <- build_csnow_basinmean(CSNOW_DIR, CROSS_VALIDATION_DIR)

# ==============================================================================
# STEP 4 - C-SNOW vs. GROUND-TRUTH (daily.csv) COMPARISON
# ==============================================================================
# The download script's own C-SNOW section says: "[NEXT STEP] Join this
# against daily.csv on station_id/date to see how well C-SNOW agrees with
# your ground-truth SWE/depth record BEFORE relying on it further - that
# comparison is not done automatically by this script." This does exactly
# that comparison and writes a small report so you have a number (not just a
# vibe) before trusting sar_depth_monthly.nc as a bias-correction predictor.
log_msg("\n----- STEP 4: C-SNOW vs. daily.csv ground-truth comparison -----")

compare_csnow_to_ground_truth <- function(csnow_dir, daily_path, out_dir,
                                          density_min_kgm3 = CSNOW_PLAUSIBLE_DENSITY_MIN_KGM3,
                                          density_max_kgm3 = CSNOW_PLAUSIBLE_DENSITY_MAX_KGM3,
                                          min_stations_per_season = CSNOW_MIN_STATIONS_PER_SEASON) {
  csnow_station_path <- file.path(csnow_dir, "CSNOW_station_extract.csv")
  if (!file.exists(csnow_station_path) || !file.exists(daily_path)) {
    log_msg("[SKIP] Need both '%s' and '%s' for this comparison.", csnow_station_path, daily_path)
    return(invisible(NULL))
  }
  csnow_df <- tryCatch(read_csv(csnow_station_path,col_types = cols(
      station_id  = col_character(),
      lon         = col_double(),
      lat         = col_double(),
      date        = col_date(),
      csnow_snd_m = col_double()   # pin this explicitly — don't let guess_max decide
    )
  ), error = function(e) NULL)
  daily_df <- tryCatch(read_csv(daily_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(csnow_df) || is.null(daily_df)) {
    log_msg("[SKIP] Could not read one of the two files.")
    return(invisible(NULL))
  }
  
  joined <- csnow_df %>%
    dplyr::mutate(date = as.Date(date)) %>%
    dplyr::inner_join(
      daily_df %>% dplyr::mutate(date = as.Date(date)) %>% dplyr::select(station_id, date, swe_mm),
      by = c("station_id", "date")
    ) %>%
    dplyr::filter(is.finite(csnow_snd_m), is.finite(swe_mm), csnow_snd_m > 0)
  if (nrow(joined) == 0) {
    log_msg("[ERROR] 0 overlapping station-date rows. CSNOW range: %s to %s | daily.csv range: %s to %s",
            min(csnow_df$date, na.rm=TRUE), max(csnow_df$date, na.rm=TRUE),
            min(daily_df$date, na.rm=TRUE), max(daily_df$date, na.rm=TRUE))
  }
  if (nrow(joined) < 10) {
    log_msg("[SKIP] Only %d overlapping station-date row(s) between C-SNOW and daily.csv - too few for a meaningful comparison.", nrow(joined))
    return(invisible(NULL))
  }
  
  # Meteorological season, used below to let the implied bulk density vary
  # seasonally (fresh/dry snow accumulates early winter; settled/wet snowpack
  # by spring) instead of a single basin-wide, all-season scalar.
  season_of_month <- function(m) dplyr::case_when(
    m %in% c(12, 1, 2)  ~ "Winter",
    m %in% c(3, 4, 5)   ~ "Spring",
    m %in% c(6, 7, 8)   ~ "Summer",
    m %in% c(9, 10, 11) ~ "Fall",
    TRUE ~ NA_character_
  )
  joined$season <- season_of_month(as.integer(format(joined$date, "%m")))
  
  # C-SNOW gives snow DEPTH (m); daily.csv gives SWE (mm). These are not
  # directly comparable in the same units, so this reports per-station
  # CORRELATION (unit-independent) and an implied bulk density
  # (swe_mm / 1000) / csnow_snd_m as a sanity range, rather than a bias in mm.
  # Kept UNGROUPED BY SEASON here deliberately - this block is the density-
  # PLAUSIBILITY GATE, and splitting each station's sample four ways before
  # gating would make legitimate stations more likely to fail the n>=10 floor
  # and get excluded outright. Seasonal density is computed separately below,
  # only for stations that already passed this whole-record gate.
  by_station <- joined %>%
    dplyr::group_by(station_id) %>%
    dplyr::summarise(
      n = dplyr::n(),
      r = suppressWarnings(cor(csnow_snd_m, swe_mm, use = "complete.obs")),
      median_implied_density_kgm3 = suppressWarnings(stats::median((swe_mm / 1000) / csnow_snd_m * 1000, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(n >= 10)
  
  if (nrow(by_station) == 0) {
    log_msg("[SKIP] No station had >=10 overlapping observations.")
    return(invisible(NULL))
  }
  
  # ------------------------------------------------------------------------
  # DENSITY GATE (was previously only a logged sanity comment - this now
  # actually flags AND excludes stations whose implied density falls outside
  # the plausible physical range, rather than leaving that judgment to
  # whoever happens to read the log). A station outside [density_min_kgm3,
  # density_max_kgm3] most likely reflects a geolocation/timing mismatch or a
  # masked C-SNOW pixel at that location, per the reasoning already in this
  # script's log messages - so it is excluded from BOTH the reported
  # correlation summary AND the basin-wide prior density used below, rather
  # than silently diluting either one.
  # ------------------------------------------------------------------------
  by_station <- by_station %>%
    dplyr::mutate(
      density_flag = dplyr::case_when(
        !is.finite(median_implied_density_kgm3) ~ "non_finite",
        median_implied_density_kgm3 < density_min_kgm3 ~ "implausible_low",
        median_implied_density_kgm3 > density_max_kgm3 ~ "implausible_high",
        TRUE ~ "plausible"
      ),
      excluded_from_summary = density_flag != "plausible"
    )
  
  n_flagged <- sum(by_station$excluded_from_summary)
  if (n_flagged > 0) {
    flagged_ids <- by_station$station_id[by_station$excluded_from_summary]
    log_msg("[FLAG] %d/%d station(s) excluded from the summary/prior-density calculation - implied density outside the plausible %d-%d kg/m3 range (likely geolocation/timing mismatch or a masked C-SNOW pixel, not a real density): %s",
            n_flagged, nrow(by_station), density_min_kgm3, density_max_kgm3, paste(flagged_ids, collapse = ", "))
  }
  
  by_station_kept <- by_station %>% dplyr::filter(!excluded_from_summary)
  
  out_path <- file.path(out_dir, "csnow_vs_groundtruth_by_station.csv")
  write_csv(by_station, out_path)  # full table, including flagged/excluded rows, for audit
  
  if (nrow(by_station_kept) == 0) {
    log_msg("[SKIP] Every station was flagged as implausible-density - no station-level summary or mm-bias conversion computed. See %s for the full (all-flagged) table.", out_path)
    return(invisible(by_station))
  }
  
  log_msg("  Compared %d station(s) total (%d retained after the density gate), %d overlapping station-date rows.",
          nrow(by_station), nrow(by_station_kept), nrow(joined))
  log_msg("  Correlation (C-SNOW depth vs. ground-truth SWE), retained stations only: median r = %.2f (range %.2f to %.2f).",
          stats::median(by_station_kept$r, na.rm = TRUE), min(by_station_kept$r, na.rm = TRUE), max(by_station_kept$r, na.rm = TRUE))
  
  # ------------------------------------------------------------------------
  # BASIN-WIDE, ALL-SEASON PRIOR DENSITY - kept as the shrinkage TARGET for
  # sparse seasons below, and as the fallback when a row's season can't be
  # determined (shouldn't happen, but keeps this from ever producing NA).
  # ------------------------------------------------------------------------
  basin_prior_density_kgm3 <- stats::median(by_station_kept$median_implied_density_kgm3, na.rm = TRUE)
  log_msg("  Basin-wide, all-season prior density (median across %d retained station(s)): %.0f kg/m3.",
          nrow(by_station_kept), basin_prior_density_kgm3)
  
  # ------------------------------------------------------------------------
  # SEASONAL PRIOR DENSITY (replaces the single all-season scalar used below).
  # Same gate-passing stations, now split by season - excluded stations stay
  # excluded here too; the whole-record gate decision above is not
  # re-litigated per season. Each season's median is empirical-Bayes shrunk
  # toward the all-season prior, weighted by how many gate-passing stations
  # actually have data in that season (n_stations / min_stations_per_season,
  # capped at 1) - same shrinkage style as H8SWEI_Snow.R's Tier 1 SWE bias
  # correction (fit_basin_bias_correction()), so a season with only 1-2
  # contributing stations doesn't hand a noisy median straight to every
  # pixel-month; a season with none at all falls back fully to the
  # all-season prior.
  # ------------------------------------------------------------------------
  joined_kept <- joined %>% dplyr::filter(station_id %in% by_station_kept$station_id)
  
  by_station_season <- joined_kept %>%
    dplyr::filter(!is.na(season)) %>%
    dplyr::group_by(station_id, season) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_implied_density_kgm3 = suppressWarnings(stats::median((swe_mm / 1000) / csnow_snd_m * 1000, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(n >= 5, is.finite(median_implied_density_kgm3))  # lighter floor than the whole-record
  # gate above - that gate already ran
  
  all_seasons <- c("Winter", "Spring", "Summer", "Fall")
  season_prior <- data.frame(season = all_seasons, stringsAsFactors = FALSE)
  season_prior$n_stations <- vapply(all_seasons, function(s) sum(by_station_season$season == s), integer(1))
  season_prior$raw_median_density_kgm3 <- vapply(all_seasons, function(s) {
    vals <- by_station_season$median_implied_density_kgm3[by_station_season$season == s]
    if (length(vals) == 0) NA_real_ else stats::median(vals, na.rm = TRUE)
  }, numeric(1))
  season_prior$weight <- pmin(1, season_prior$n_stations / min_stations_per_season)
  season_prior$density_kgm3 <- ifelse(
    is.finite(season_prior$raw_median_density_kgm3),
    season_prior$weight * season_prior$raw_median_density_kgm3 + (1 - season_prior$weight) * basin_prior_density_kgm3,
    basin_prior_density_kgm3   # no gate-passing station had ANY data this season - fall back fully
  )
  
  log_msg("  Seasonal prior density (empirical-Bayes shrunk toward the %.0f kg/m3 all-season prior where a season has < %d contributing station(s)):",
          basin_prior_density_kgm3, min_stations_per_season)
  for (i in seq_len(nrow(season_prior))) {
    log_msg("    %-7s: %.0f kg/m3  (n_stations = %d, shrinkage weight = %.2f)",
            season_prior$season[i], season_prior$density_kgm3[i], season_prior$n_stations[i], season_prior$weight[i])
  }
  
  density_by_season <- setNames(season_prior$density_kgm3, season_prior$season)
  
  out_path_season_density <- file.path(out_dir, "csnow_vs_groundtruth_density_by_season.csv")
  write_csv(season_prior, out_path_season_density)
  log_msg("  Seasonal density prior table written to %s.", out_path_season_density)
  
  # ------------------------------------------------------------------------
  # APPROXIMATE mm-BIAS, now using each row's OWN SEASON's (shrunk) density
  # rather than a single basin-wide, all-season scalar. Still NOT a real
  # SWE-vs-SWE validation - it is an approximate, explicitly-labeled
  # conversion, at the cost of whatever error the assumed density (now
  # seasonally-varying, still not elevation/aspect-varying) introduces.
  # ------------------------------------------------------------------------
  joined_kept$density_used_kgm3 <- ifelse(
    !is.na(joined_kept$season) & joined_kept$season %in% names(density_by_season),
    density_by_season[joined_kept$season],
    basin_prior_density_kgm3
  )
  joined_kept$csnow_implied_swe_mm <- joined_kept$csnow_snd_m * joined_kept$density_used_kgm3
  approx_bias_mm <- mean(joined_kept$swe_mm - joined_kept$csnow_implied_swe_mm, na.rm = TRUE)
  approx_rmse_mm <- sqrt(mean((joined_kept$swe_mm - joined_kept$csnow_implied_swe_mm)^2, na.rm = TRUE))
  
  log_msg("  [APPROXIMATE, density-dependent] Using the SEASONAL prior density to convert C-SNOW depth to an implied SWE:")
  log_msg("    mean bias (ground-truth SWE minus C-SNOW-implied SWE) = %.1f mm; RMSE = %.1f mm, n = %d station-date rows (all seasons pooled).",
          approx_bias_mm, approx_rmse_mm, nrow(joined_kept))
  
  bias_by_season <- joined_kept %>%
    dplyr::filter(!is.na(season)) %>%
    dplyr::group_by(season) %>%
    dplyr::summarise(
      n = dplyr::n(),
      density_used_kgm3 = dplyr::first(density_used_kgm3),
      mean_bias_mm = mean(swe_mm - csnow_implied_swe_mm, na.rm = TRUE),
      rmse_mm = sqrt(mean((swe_mm - csnow_implied_swe_mm)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  log_msg("  Per-season mm-bias (same conversion, broken out - this is the actual payoff of the seasonal split):")
  for (i in seq_len(nrow(bias_by_season))) {
    log_msg("    %-7s: mean bias = %6.1f mm | RMSE = %6.1f mm | density = %.0f kg/m3 | n = %d",
            bias_by_season$season[i], bias_by_season$mean_bias_mm[i], bias_by_season$rmse_mm[i],
            bias_by_season$density_used_kgm3[i], bias_by_season$n[i])
  }
  log_msg("  This mm-bias is UNCERTAIN and density-dependent (density now varies by SEASON but still not by")
  log_msg("  elevation/aspect within a season) - treat it as an approximate sanity number alongside the")
  log_msg("  correlation above, not as a substitute for a real SWE-vs-SWE cross-validation product like CMC.")
  
  approx_bias_path <- file.path(out_dir, "csnow_vs_groundtruth_approx_mm_bias.csv")
  write_csv(
    data.frame(
      n_stations_total = nrow(by_station),
      n_stations_retained = nrow(by_station_kept),
      n_stations_flagged = n_flagged,
      basin_prior_density_kgm3 = basin_prior_density_kgm3,
      approx_mean_bias_mm = approx_bias_mm,
      approx_rmse_mm = approx_rmse_mm,
      n_station_date_rows = nrow(joined_kept)
    ),
    approx_bias_path
  )
  log_msg("  Approximate mm-bias summary (all-seasons-pooled) written to %s.", approx_bias_path)
  
  approx_bias_season_path <- file.path(out_dir, "csnow_vs_groundtruth_approx_mm_bias_by_season.csv")
  write_csv(bias_by_season, approx_bias_season_path)
  log_msg("  Approximate mm-bias summary (by season) written to %s.", approx_bias_season_path)
  
  invisible(by_station)
}

csnow_vs_truth <- compare_csnow_to_ground_truth(CSNOW_DIR, STATION_OBS_PATH, BRIDGE_DIAG_DIR)

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
log_msg("\n============================================================")
log_msg("BRIDGE COMPLETE")
log_msg("============================================================")
log_msg("Outputs that H8SWEI_Snow.R will now be able to find (if produced above):")
log_msg("  %s : %s", SAR_DEPTH_PATH_OUT, ifelse(file.exists(SAR_DEPTH_PATH_OUT), "PRESENT", "still missing - see STEP 2 log above"))
log_msg("  %s/cmc_swe_monthly_basinmean.csv : %s", CROSS_VALIDATION_DIR,
        ifelse(file.exists(file.path(CROSS_VALIDATION_DIR, "cmc_swe_monthly_basinmean.csv")), "PRESENT", "still missing - see STEP 3a log above"))
log_msg("  %s/csnow_sar_monthly_basinmean.csv : %s", CROSS_VALIDATION_DIR,
        ifelse(file.exists(file.path(CROSS_VALIDATION_DIR, "csnow_sar_monthly_basinmean.csv")), "PRESENT", "still missing - see STEP 3b log above"))
log_msg("  CanSWE_daily.csv at root (H8's default CANSWE_STATION_OBS_PATH) : %s",
        ifelse(file.exists("CanSWE_daily.csv"), "PRESENT", "still missing"))
log_msg("\nStill NOT produced (would need new downloaders, not just re-arranging what's on disk):")
log_msg("  AU_DySno, MODIS-VIIRS, HLS, Crocus-ERA5 monthly basin-mean CSVs")
log_msg("  (H8 will continue to skip these cross-validation checks until real data exists for them)")

writeLines(BRIDGE_LOG, file.path(BRIDGE_DIAG_DIR, sprintf("bridge_run_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S"))))
cat(sprintf("\n[OK] Full log written to %s\n", file.path(BRIDGE_DIAG_DIR, "bridge_run_log_*.txt")))
cat("\n--- READY FOR H8SWEI_Snow.R ---\n")