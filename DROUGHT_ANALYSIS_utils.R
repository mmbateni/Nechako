####################################################################################
# DROUGHT_ANALYSIS_utils.R  ·  SHARED UTILITIES FOR ALL NECHAKO DROUGHT SCRIPTS
# ─────────────────────────────────────────────────────────────────────────────────
# Sourced by:
#   6trend_test_ALL.R              — drought index trend analysis
#   4pr_pet_trends.r             — Pr / PET / Temperature trend analysis
#   (and all other w* scripts)
#
# SECTIONS
#   A.  Project configuration (paths, constants, thresholds)
#   B.  Basin-averaged CSV loader / writer / package bootstrap
#   C.  NetCDF file finders
#   D.  NetCDF date parsing
#   E.  Basin spatial helpers
#   F.  Drought event detection
#   F2. Spatial helpers (rasterise / smooth / plot)
#   G.  ggplot helpers
#   H.  Excel summary export
#   I.  w9 shared-state loader
#   J.  [SHARED STATS]  Functions used by BOTH 6trend_test_ALL.R and
#       4pr_pet_trends.r — moved here to eliminate duplication:
#         J1. Utility helpers  (check_min_value_threshold,
#                               calculate_variance_with_ties)
#         J2. PELT changepoint  (detect_regime_shift_pelt)
#         J3. VC-MK per-series  (.mk_vc_single)
#         J4. TFPW-MK per-series (.mk_tfpw_single)
#         J5. Spectral peak detection per-series (.spectral_single)
#         J6. Omnibus single-series wrapper (mk_tfpw_spectral_for_series)
#             — used directly by 4pr_pet_trends.r
#         J7. Vectorised pixel-matrix wrappers  (vectorized_mann_kendall,
#             vectorized_mann_kendall_tfpw, vectorized_regime_shift_pelt,
#             vectorized_spectral_peaks)
#             — used directly by 6trend_test_ALL.R
####################################################################################

## A. PROJECT CONFIGURATION ───────────────────────────────────────────
WD_PATH <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")
BASIN_SHP       <- file.path(WD_PATH, "Spatial/nechakoBound_dissolve.kmz")

## Private helper: unzip a KMZ file into a self-cleaning temp directory and
## call FUN(kml_path), returning FUN's result.  The temp directory is deleted
## by on.exit() when this function returns — so FUN must finish reading the
## KML before .read_kmz() exits (which it always does, since FUN is called
## synchronously inside the function body).
.read_kmz <- function(kmz_path, FUN) {
  tmp <- tempfile()
  dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(kmz_path, exdir = tmp)
  kml <- list.files(tmp, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)
  if (!length(kml))
    stop(".read_kmz: no .kml file found inside ", basename(kmz_path))
  FUN(kml[1])
}

## Helper: load basin as a single dissolved SpatVector (terra).
## terra::vect() cannot read KMZ directly on Windows (missing GDAL KMZ driver).
## .read_kmz() handles the unzip/cleanup; terra::vect() reads the inner KML.
load_basin_vect <- function(path = BASIN_SHP) {
  if (!file.exists(path)) stop("Basin file not found: ", path)
  v <- if (tolower(tools::file_ext(path)) == "kmz")
    .read_kmz(path, terra::vect)
  else
    terra::vect(path)
  if (nrow(v) > 1L) v <- terra::aggregate(v)
  v
}

## Index-specific directories
SPI_SEAS_DIR      <- file.path(WD_PATH, "spi_results_seasonal/")
SPEI_SEAS_DIR     <- file.path(WD_PATH, "spei_results_seasonal/")
SPEI_THW_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal_thw/")   # Thornthwaite PET
SWEI_SEAS_DIR     <- file.path(WD_PATH, "swei_results_seasonal/")

## Output directories
TREND_DIR      <- file.path(WD_PATH,  "temporal_drought/")
BASIN_PLOT_DIR <- file.path(WD_PATH,  "basin_averaged_plots/")
POINT_PLOT_DIR <- file.path(WD_PATH,  "point_timeseries_plots/")
GIF_DIR        <- file.path(TREND_DIR, "drought_gifs/")
CACHE_DIR      <- file.path(WD_PATH,  "temporal_drought/cache/")
TELE_DIR       <- file.path(WD_PATH,  "teleconnections")

## Figure output constants (used by save_figure() and all w* scripts)
FIG_DPI        <- 300L   # PNG resolution
FIG_WIDTH_WIDE <- 14     # wide panel figures (inches)
FIG_WIDTH_STD  <- 10     # standard single-panel figures (inches)
FIG_HEIGHT_STD <-  8     # standard figure height (inches)

## CRS and thresholds
EQUAL_AREA_CRS    <- "EPSG:3005"
DROUGHT_ONSET     <- -0.5
DROUGHT_END       <- -0.5


## Timescales
SPI_SCALES          <- c(1, 2, 3, 6, 12)
SPEI_SCALES         <- c(1, 2, 3, 6, 12)
SWEI_SCALE          <- 3  # Single timescale only
INDICES             <- c("spi", "spei", "spei_thw", "swei")
TIMESCALES_STANDARD <- SPI_SCALES          # used by some codes for time series plots
TIMESCALES_SPATIAL  <- SPI_SCALES          # used by spme codes for spatial figures
MONTH_NAMES    <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## A2. PRECIPITATION UNIT DETECTION & SAFE SCALING ───────────────────
#' Two companion functions:
#'   detect_precip_unit()   — inspect one layer and return the conversion factor
#'   apply_precip_scaling() — apply that factor + days_in_month safely (use this
#'                            in callers instead of writing the if/multiply inline)
#'
#' Detect the unit of a precipitation (or PET) raster by sampling a summer
#' layer and comparing the basin-mean value against known physical ranges.
#'
#' ERA5-Land CDS "monthly averaged reanalysis" stores both total_precipitation
#' and potential_evapotranspiration as mean DAILY RATES in metres/day (m/day).
#' Standard pre-processing multiplies by 1000 to get mm/day, then by
#' days_in_month to get mm/month.  When the × 1000 step is accidentally
#' skipped the values are 1000× too small, which is the root cause of the
#' WB_PM range bug (Issue 1).
#'
#' The function performs a fast, non-destructive test on one summer layer
#' (July preferred; falls back to the first available layer) and returns the
#' scalar conversion factor needed to bring the raster to mm/day:
#'   • 1000  — raw m/day   (ERA5-Land straight from CDS, no prior scaling)
#'   •    1  — already mm/day  (standard pre-processed input)
#'   •   NA  — ambiguous; value printed with a warning so the caller can decide
#'
#' The thresholds are calibrated for a cool boreal/montane basin (Nechako):
#'   < 0.02  m/day  → m/day   (summer mean ~0.001–0.005 m/day → need × 1000)
#'   0.02–20 mm/day → mm/day  (summer mean ~0.5–5 mm/day for Nechako July)
#'   > 20           → mm/month already (do NOT multiply by days_in_month again)
#'
#' @param raster_obj  SpatRaster loaded with terra::rast().
#' @param dates_vec   Date vector of length nlyr(raster_obj) giving each layer
#'                    date.  If NULL the function falls back to layer 1.
#' @param var_label   Character label printed in diagnostic messages
#'                    (e.g. "Precip", "PET").
#' @param july_only   Logical.  If TRUE (default) sample the first July layer;
#'                    otherwise sample the first available layer.
#'
#' @return Named list:
#'   $factor   numeric — multiply raster by this to get mm/day
#'             (1000, 1, or NA if ambiguous)
#'   $unit_in  character — detected input unit label
#'   $mean_val numeric — basin-mean of the sampled layer (original units)
detect_precip_unit <- function(raster_obj,
                               dates_vec  = NULL,
                               var_label  = "variable",
                               july_only  = TRUE) {
  
  ## ── Pick a representative layer ────────────────────────────────────────────
  layer_idx <- NA_integer_
  if (!is.null(dates_vec) && length(dates_vec) == terra::nlyr(raster_obj)) {
    mon_nums <- as.integer(format(dates_vec, "%m"))
    if (july_only) {
      july_hits <- which(mon_nums == 7L)
      layer_idx <- if (length(july_hits)) july_hits[1L] else which(!is.na(mon_nums))[1L]
    } else {
      layer_idx <- which(!is.na(mon_nums))[1L]
    }
  }
  if (is.na(layer_idx) || layer_idx < 1L) layer_idx <- 1L
  
  ## ── Sample basin-mean of that layer ────────────────────────────────────────
  vals <- terra::values(raster_obj[[layer_idx]], mat = FALSE)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    warning(sprintf(
      "detect_precip_unit [%s]: all values NA in layer %d — cannot detect unit.",
      var_label, layer_idx))
    return(list(factor = NA_real_, unit_in = "unknown", mean_val = NA_real_))
  }
  mean_val <- mean(vals, na.rm = TRUE)
  
  ## ── Classify ───────────────────────────────────────────────────────────────
  # Boundary 0.02: highest plausible m/day daily rate for a wet boreal basin
  # in summer is ~0.01 m/day; lowest plausible mm/day daily rate is ~0.2 mm/day.
  # Boundary 20: max plausible mm/day for any month; above that it must be mm/month.
  if (mean_val < 0.02) {
    factor  <- 1000
    unit_in <- "m/day"
  } else if (mean_val <= 20) {
    factor  <- 1
    unit_in <- "mm/day"
  } else {
    factor  <- NA_real_   # already mm/month — caller must not re-scale by days
    unit_in <- "mm/month"
    warning(sprintf(
      "detect_precip_unit [%s]: mean = %.2f looks like mm/month already. ",
      var_label, mean_val),
      "Do NOT multiply by days_in_month again.")
  }
  
  cat(sprintf(
    "  [detect_precip_unit] %s layer %d: mean = %.5f → detected '%s' → factor = %s\n",
    var_label, layer_idx, mean_val, unit_in,
    if (is.na(factor)) "NA (mm/month — skip scaling)" else as.character(factor)))
  
  list(factor = factor, unit_in = unit_in, mean_val = mean_val)
}

#' Apply the unit conversion returned by detect_precip_unit() and then scale
#' by days_in_month — in the correct order and with a hard stop if the input
#' is already in mm/month (factor = NA), which would otherwise cause a silent
#' ×28–31 over-inflation if the days_in_month multiply ran unconditionally.
#'
#' Replace the unsafe inline pattern:
#'   if (!is.na(precip_unit$factor) && precip_unit$factor != 1)
#'     precip <- precip * precip_unit$factor
#'   precip <- precip * days_in_month          # ← bug: always ran, even for NA
#'
#' With:
#'   precip <- apply_precip_scaling(precip, precip_unit, days_in_month)
#'
#' @param raster_obj    SpatRaster to scale (modified in-place via terra copy-on-modify).
#' @param precip_unit   List returned by detect_precip_unit() — must have $factor.
#' @param days_in_month Integer vector of length nlyr(raster_obj) giving the number
#'                      of days in each layer's month (used for mm/day → mm/month).
#' @return Scaled SpatRaster in mm/month.
apply_precip_scaling <- function(raster_obj, precip_unit, days_in_month) {
  
  factor <- precip_unit$factor
  
  ## ── Guard: mm/month input must not be re-scaled by days ────────────────────
  if (is.na(factor)) {
    stop(
      "apply_precip_scaling: detect_precip_unit returned factor = NA, which means\n",
      "  the input raster appears to be in mm/month already (mean > 20).\n",
      "  Multiplying by days_in_month would inflate every value by ~28-31x.\n",
      "  Review your input NetCDF: if it really is mm/month, return it as-is\n",
      "  without any further scaling.  If it is mm/day, re-check the detection\n",
      "  thresholds in detect_precip_unit() for your study region.")
  }
  
  ## ── Step 1: m/day → mm/day (only when factor = 1000) ──────────────────────
  if (factor != 1) {
    cat(sprintf("  [apply_precip_scaling] Applying × %.0f (%s → mm/day)\n",
                factor, precip_unit$unit_in))
    raster_obj <- raster_obj * factor
  }
  
  ## ── Step 2: mm/day → mm/month ─────────────────────────────────────────────
  raster_obj <- raster_obj * days_in_month
  cat("  [apply_precip_scaling] Days-in-month scaling applied → mm/month\n")
  
  raster_obj
}

## B. BASIN-AVERAGED CSV LOADER ─────────────────────────────────────
#' Load a pre-computed basin-averaged drought index CSV and return a tidy
#' long-format data.frame(date, value).
#'
#' AUTHORITATIVE SOURCE: these CSVs are produced by 6trend_test_ALL.R via
#' save_basin_avg_from_pixels() (see Section B2 below).  The methodology is:
#'   (1) compute the drought index at every basin pixel from per-pixel inputs,
#'   (2) take the area-weighted spatial mean of those pixel-level index values.
#' This is the standard hydrometeorological approach: area-weighted averaging
#' of standardised index values that were calibrated on the same spatial basis.
#' Run 6trend_test_ALL.R before basin_timeseries.R.
#'
#' File naming convention: {index}_{scale:02d}_basin_averaged_by_month.csv
#' File structure:
#'   Year | Jan | Feb | Mar | … | Dec
#'   1950 | -1.88 |  NA | …
#'
#' @param index_type  "spi", "spei", or "swei"
#' @param scale       Integer timescale (e.g. 1, 3, 6, 12)
#' @param data_dir    Directory containing the CSV; if NULL, the index-specific
#'                    seasonal directory from section A is used automatically.
#' @return data.frame with columns \code{date} (Date) and \code{value} (numeric),
#'         sorted by date with NA rows removed.  Returns NULL with a warning if
#'         the file is not found (usually means w3 has not been run yet).
load_basin_avg_csv <- function(index_type, scale, data_dir = NULL) {
  dir_map <- list(spi      = SPI_SEAS_DIR,
                  spei     = SPEI_SEAS_DIR,
                  spei_thw = SPEI_THW_SEAS_DIR,
                  swei     = SWEI_SEAS_DIR)
  dir     <- if (!is.null(data_dir)) data_dir else {
    d <- dir_map[[tolower(index_type)]]
    if (is.null(d)) stop("Unknown index_type: ", index_type)
    d
  }
  f <- file.path(dir, sprintf("%s_%02d_basin_averaged_by_month.csv",
                              tolower(index_type), as.integer(scale)))
  if (!file.exists(f)) {
    warning(sprintf("Basin-averaged CSV not found: %s", f))
    return(NULL)
  }
  df <- as.data.frame(data.table::fread(f))  # faster than read.csv for wide CSVs
  if (!"Year" %in% names(df))
    stop("Expected 'Year' column in ", basename(f))
  mon_cols <- intersect(MONTH_NAMES, names(df))
  if (!length(mon_cols))
    stop("No month name columns found in ", basename(f),
         ".  Expected: ", paste(MONTH_NAMES, collapse = ", "))
  long <- do.call(rbind, lapply(seq_along(mon_cols), function(mi) {
    data.frame(
      date  = as.Date(paste(df$Year, mi, "01", sep = "-")),
      value = as.numeric(df[[mon_cols[mi]]]),
      stringsAsFactors = FALSE)
  }))
  long <- long[order(long$date), ]
  long <- long[!is.na(long$value) & is.finite(long$value), ]
  rownames(long) <- NULL
  long
}

## B2. BASIN-AVERAGED CSV WRITER ────────────────────────────────────
#' Compute an area-weighted basin-average time series from a pixel × month
#' matrix and write the result as the standard basin_averaged_by_month.csv
#' expected by load_basin_avg_csv().
#'
#' Called by 6trend_test_ALL.R::process_index() after the pixel time-series
#' matrix has been assembled.  The area weights are the equal-area cell sizes
#' extracted from the raster template (already in EPSG:3005).
#'
#' Methodology:
#'   basin_avg[t] = Σ( index_value[pixel, t] × cell_area[pixel] )
#'                / Σ( cell_area[pixel] )          (sum over non-NA pixels)
#'
#' @param ts_matrix   Numeric matrix (n_valid_pixels × n_months).  Rows
#'                    correspond to valid basin pixels; columns to calendar
#'                    months in chronological order.
#' @param area_weights Numeric vector of length n_valid_pixels; cell areas in
#'                    m² in the equal-area CRS (from terra::cellSize on the
#'                    BC Albers projected raster template).
#' @param dates       Date vector of length n_months giving the first day of
#'                    each month represented by a column of ts_matrix.
#' @param index_type  "spi", "spei", "swei", or any other index key.
#' @param scale       Integer timescale (e.g. 1, 3, 6, 12).
#' @param out_dir     Directory to write the CSV.  Should be the index-specific
#'                    seasonal results directory (SPI_SEAS_DIR etc.) so that
#'                    load_basin_avg_csv() can find it without specifying data_dir.
#' @return Invisible data.frame (Year × month-name columns) as written to disk.
save_basin_avg_from_pixels <- function(ts_matrix, area_weights, dates,
                                       index_type, scale, out_dir) {
  ## ── Input validation ───────────────────────────────────────────────
  if (!is.matrix(ts_matrix))
    stop("ts_matrix must be a numeric matrix (pixels × months)")
  if (nrow(ts_matrix) != length(area_weights))
    stop(sprintf(
      "ts_matrix has %d rows but area_weights has %d elements",
      nrow(ts_matrix), length(area_weights)))
  if (ncol(ts_matrix) != length(dates))
    stop(sprintf(
      "ts_matrix has %d columns but dates has %d elements",
      ncol(ts_matrix), length(dates)))
  
  ## ── Guard: replace zero / negative / NA areas with column median ───
  bad_w <- is.na(area_weights) | !is.finite(area_weights) | area_weights <= 0
  if (any(bad_w)) {
    med_w <- median(area_weights[!bad_w], na.rm = TRUE)
    area_weights[bad_w] <- if (is.finite(med_w) && med_w > 0) med_w else 1
    warning(sprintf(
      "save_basin_avg_from_pixels: %d bad area weight(s) replaced with %.0f m²",
      sum(bad_w), area_weights[which(bad_w)[1]]))
  }
  
  ## ── Area-weighted mean for every time step ─────────────────────────
  basin_avg <- vapply(seq_len(ncol(ts_matrix)), function(t) {
    v  <- ts_matrix[, t]
    ok <- !is.na(v) & is.finite(v)
    if (!any(ok)) return(NA_real_)
    sum(v[ok] * area_weights[ok]) / sum(area_weights[ok])
  }, numeric(1L))
  
  ## ── Reshape to Year × Month wide format ────────────────────────────
  yrs <- as.integer(format(dates, "%Y"))
  mos <- as.integer(format(dates, "%m"))
  all_years <- sort(unique(yrs))
  
  df <- data.frame(Year = all_years, stringsAsFactors = FALSE)
  for (mi in seq_along(MONTH_NAMES)) {
    hits    <- which(mos == mi)
    row_idx <- match(yrs[hits], all_years)        # vectorised: all hits at once
    valid   <- !is.na(row_idx)
    col_vals <- rep(NA_real_, length(all_years))
    col_vals[row_idx[valid]] <- basin_avg[hits[valid]]
    df[[MONTH_NAMES[mi]]] <- col_vals
  }
  
  ## ── Write CSV ───────────────────────────────────────────────────────
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(
    out_dir,
    sprintf("%s_%02d_basin_averaged_by_month.csv",
            tolower(index_type), as.integer(scale)))
  utils::write.csv(df, out_file, row.names = FALSE)
  cat(sprintf(
    "  ✓ Basin-avg CSV written (area-weighted, %d pixels, %d months): %s\n",
    nrow(ts_matrix), ncol(ts_matrix), basename(out_file)))
  invisible(df)
}

## B3. PACKAGE BOOTSTRAP ──────────────────────────────────────────────
utils_load_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org/")
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
  cat("✓ Packages loaded:", paste(pkgs, collapse = ", "), "\n")
}

####################################################################################
# C. NETCDF FILE FINDERS ────────────────────────────────────────────
####################################################################################
find_seasonal_nc_files <- function(data_dir, index_type, scale) {
  idx <- tolower(index_type)
  pattern <- sprintf("^%s_%02d_month\\d{2}_[A-Za-z]+\\.nc$", idx, scale)
  files <- list.files(data_dir, pattern = pattern,
                      full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  if (!length(files)) {
    warning(sprintf("No seasonal files found for %s-%d in: %s",
                    toupper(index_type), scale, data_dir))
    return(character(0))
  }
  mon_nums <- as.integer(regmatches(basename(files),
                                    regexpr("(?<=month)\\d{2}", basename(files), perl = TRUE)))
  files[order(mon_nums)]
}

find_swei_seasonal_files <- function(data_dir, index_type, scale) {
  pattern <- sprintf("^swei_%02d_month\\d{2}_[A-Za-z]+\\.nc$", scale)
  files <- list.files(data_dir, pattern = pattern,
                      full.names = TRUE, ignore.case = TRUE, recursive = TRUE)
  if (!length(files)) {
    warning(sprintf("No SWEI seasonal files found for scale %d in: %s", scale, data_dir))
    return(character(0))
  }
  mon_nums <- as.integer(regmatches(basename(files),
                                    regexpr("(?<=month)\\d{2}", basename(files), perl = TRUE)))
  files[order(mon_nums)]
}

## D. NETCDF DATE PARSING ────────────────────────────────────────────
extract_dates_from_nc <- function(nc_file, n_layers = NULL) {
  r_temp <- terra::rast(nc_file)
  n_layers <- if (is.null(n_layers)) terra::nlyr(r_temp) else n_layers
  dates  <- tryCatch({
    d  <- terra::time(r_temp)
    if (!is.null(d) && length(d) > 0) {
      d_conv  <- tryCatch(as.Date(d), error = function(e) as.Date(d, origin = "1970-01-01"))
      if (all(!is.na(d_conv)) &&
          min(d_conv) >= as.Date("1900-01-01") &&
          max(d_conv) <= as.Date("2030-01-01")) {
        return(.pad_dates(d_conv, n_layers))
      }
    }
    NULL
  }, error = function(e) NULL)
  if (is.null(dates)) {
    cat("  [dates] Using fallback monthly sequence (Jan 1950)\n")
    d <- seq(as.Date("1950-01-01"), by = "month", length.out = n_layers)
    if (tail(d, 1) > as.Date("2025-12-31")) {
      end <- as.Date("2025-12-01")
      start <- end - lubridate::months(n_layers - 1)
      d <- seq(start, by = "month", length.out = n_layers)
    }
    d
  } else dates
}

.pad_dates <- function(d, n) {
  if (length(d) == n) return(d)
  if (length(d) > n) return(d[1:n])
  extra <- seq(tail(d, 1) + lubridate::months(1), by = "month", length.out = n - length(d))
  c(d, extra)
}

## E. BASIN SPATIAL HELPERS ──────────────────────────────────────────
precompute_basin_geometry  <- function(basin_shp_path, raster_file_path, crs = EQUAL_AREA_CRS) {
  basin_v  <- terra::vect(basin_shp_path)
  r_tmpl  <- terra::rast(raster_file_path)[[1]]
  r_proj  <- terra::project(r_tmpl, crs)
  b_proj  <- terra::project(basin_v, crs)
  total_area  <- sum(terra::expanse(b_proj, unit = "m"))
  cell_areas  <- terra::cellSize(r_proj, unit = "m")
  areas_masked  <- terra::mask(cell_areas, b_proj)
  basin_rast  <- terra::rasterize(b_proj, r_proj, cover = TRUE)
  area_v  <- as.vector(terra::values(areas_masked, na.rm = FALSE))
  cover_v  <- as.vector(terra::values(basin_rast, na.rm = FALSE))
  eff_areas  <- area_v * cover_v
  valid_idx  <- which(!is.na(cover_v) & cover_v > 0)
  cat(sprintf("  Basin area: %.1f km²  |  Valid cells: %d\n",
              total_area / 1e6, length(valid_idx)))
  list(total_basin_area = total_area,
       effective_areas = eff_areas,
       valid_cell_idx = valid_idx,
       equal_area_crs = crs,
       basin_shp_path = basin_shp_path,
       raster_file_path = raster_file_path)
}

compute_layer_statistics  <- function(layer_idx, raster_file, basin_geom, threshold = DROUGHT_THRESHOLD) {
  r  <- terra::rast(raster_file)[[layer_idx]]
  r  <- terra::project(r, basin_geom$equal_area_crs)
  v  <- as.vector(terra::values(r, na.rm = FALSE))
  vc  <- v[basin_geom$valid_cell_idx]
  ac  <- basin_geom$effective_areas[basin_geom$valid_cell_idx]
  ok  <- !is.na(vc) & is.finite(vc)
  if (!any(ok)) return(list(mean_value = NA_real_, valid_fraction = 0,
                            drought_fraction = NA_real_, drought_area = 0,
                            n_cells_valid = 0, n_cells_drought = 0))
  valid_area  <- sum(ac[ok])
  mean_val  <- sum(vc[ok] * ac[ok]) / valid_area
  valid_frac  <- valid_area / basin_geom$total_basin_area
  drought_ok  <- ok & (vc <= threshold)
  drought_area  <- sum(ac[drought_ok])
  list(mean_value = mean_val,
       valid_fraction = valid_frac,
       drought_fraction = drought_area / basin_geom$total_basin_area,
       drought_area = drought_area,
       n_cells_valid = sum(ok),
       n_cells_drought = sum(drought_ok))
}
## E2. Area-weighted pixel-matrix column means ──────────────────────────────
#' Compute area-weighted column means of a pixel × time matrix.
#'
#' @param mat     Numeric matrix [n_all_cells × n_time].
#'                NA for non-basin or missing cells.
#' @param weights Numeric vector [n_all_cells].
#'                Set to 0 (not NA) for non-basin cells.
#' @return Numeric vector of length n_time.
area_weighted_colmeans <- function(mat, weights) {
  apply(mat, 2L, function(col) {
    ok <- !is.na(col) & is.finite(col) & weights > 0
    if (!any(ok)) return(NA_real_)
    sum(col[ok] * weights[ok]) / sum(weights[ok])
  })
}
## F. DROUGHT EVENT DETECTION ────────────────────────────────────────
detect_drought_events <- function(df,
                                  onset_threshold       = -0.5,
                                  termination_threshold = -0.5,
                                  min_duration          = 1,
                                  severe_thr            = -1.5,
                                  extreme_thr           = -2.0) {
  stopifnot(all(c("date", "value") %in% names(df)))
  df    <- df[order(df$date), ]
  vals  <- df$value
  dates <- df$date
  
  # Canonical empty frame — returned when there are no events
  empty_events <- data.frame(
    event_id        = integer(),
    start_date      = as.Date(character()),
    end_date        = as.Date(character()),
    duration_months = integer(),
    min_value       = numeric(),
    severity        = character(),
    stringsAsFactors = FALSE)
  
  classify_severity <- function(min_v) {
    if      (min_v <= extreme_thr) "D3: extremely dry"
    else if (min_v <= severe_thr)  "D2: severely dry"
    else                           "D1: moderately dry"
  }
  
  event_list <- list()
  close_event <- function(s, e) {
    dur <- e - s + 1L
    if (dur < min_duration) return()
    mv <- min(vals[s:e], na.rm = TRUE)
    event_list[[length(event_list) + 1L]] <<- data.frame(
      event_id        = length(event_list) + 1L,
      start_date      = dates[s],
      end_date        = dates[e],
      duration_months = dur,
      min_value       = mv,
      severity        = classify_severity(mv),
      stringsAsFactors = FALSE)
  }
  
  in_drought  <- FALSE
  start_idx   <- NA_integer_
  for (i in seq_along(vals)) {
    if (!in_drought && !is.na(vals[i]) && vals[i] < onset_threshold) {
      in_drought  <- TRUE
      start_idx   <- i
    } else if (in_drought && !is.na(vals[i]) && vals[i] >= termination_threshold) {
      close_event(start_idx, i - 1L)
      in_drought  <- FALSE
    }
  }
  if (in_drought) close_event(start_idx, length(vals))
  
  if (!length(event_list)) return(empty_events)
  do.call(rbind, event_list)
}

## F2. SPATIAL HELPERS (used by 4_trends_visualization.R) ──────────
#' Load basin shapefile and reproject to equal-area CRS.
#' KMZ files are handled via the shared .read_kmz() helper.
load_basin <- function(shp_path = BASIN_SHP, crs = EQUAL_AREA_CRS) {
  if (!file.exists(shp_path)) stop("Basin file not found: ", shp_path)
  b <- if (tolower(tools::file_ext(shp_path)) == "kmz")
    .read_kmz(shp_path, function(kml) sf::st_read(kml, quiet = TRUE))
  else
    sf::st_read(shp_path, quiet = TRUE)
  if (nrow(b) > 1L) b <- sf::st_as_sf(sf::st_union(b))
  if (sf::st_crs(b)$input != crs) b <- sf::st_transform(b, crs)
  cat(sprintf("✓ Basin loaded (%d polygon(s) → dissolved, CRS: %s)\n", nrow(b), crs))
  b
}

#' Clip a data.table of lon/lat points to basin using sf::st_intersects
clip_to_basin  <- function(dt, basin_sf, crs = EQUAL_AREA_CRS) {
  basin_v    <- terra::vect(basin_sf)
  basin_ext  <- terra::ext(basin_v)
  dt_pre     <- dt[lon >= basin_ext[1] & lon <= basin_ext[2] &
                     lat >= basin_ext[3] & lat <= basin_ext[4]]
  if (nrow(dt_pre) == 0) { warning("No points within basin extent"); return(dt_pre) }
  pts_sf  <- sf::st_as_sf(terra::vect(dt_pre, geom = c("lon", "lat"), crs = crs))
  idx     <- sapply(sf::st_intersects(pts_sf, basin_sf, sparse = TRUE), length) > 0
  dt_pre[idx, ]
}

#' Create a raster template aligned to basin extent
create_raster_template  <- function(data_dt, basin_sf) {
  data_dt  <- data_dt[!is.na(lon) & !is.na(lat)]
  if (nrow(data_dt) == 0) return(NULL)
  u_lons  <- unique(data_dt$lon); u_lats  <- unique(data_dt$lat)
  if (length(u_lons) >= 2 && length(u_lats) >= 2) {
    dx   <- median(diff(sort(u_lons)), na.rm = TRUE)
    dy   <- median(diff(sort(u_lats)), na.rm = TRUE)
    res  <- min(dx, dy, na.rm = TRUE)
  } else {
    ext  <- terra::ext(terra::vect(basin_sf))
    res  <- min(ext[2]-ext[1], ext[4]-ext[3]) / 50
  }
  if (!is.finite(res) || res <= 0) res  <- 5000
  terra::rast(ext = terra::ext(terra::vect(basin_sf)),
              resolution = res, crs = EQUAL_AREA_CRS)
}

#' Rasterise a column from a data.table onto a template and mask to basin
create_raster_from_points <- function(dt, template, value_col, basin_sf) {
  if (is.null(template)) return(NULL)
  empty <- function() {
    r <- terra::rast(template); terra::values(r) <- NA
    terra::mask(r, terra::vect(basin_sf))
  }
  if (nrow(dt) == 0 || !value_col %in% names(dt)) return(empty())
  dt <- dt[!is.na(lon) & !is.na(lat)]
  if (nrow(dt) == 0 || all(is.na(dt[[value_col]]))) return(empty())
  pts <- terra::vect(dt, geom = c("lon","lat"), crs = terra::crs(template))
  r   <- terra::rasterize(pts, template, field = value_col, touches = TRUE)
  terra::mask(r, terra::vect(basin_sf))
}

#' 3×3 focal mean smoother (passes NA through gracefully)
smooth_raster <- function(r) {
  if (is.null(r) || all(is.na(terra::values(r)))) return(r)
  terra::focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
}

#' Plot a raster with basin overlay; handles all-NA gracefully.
#' Relies on basin_sf_global being set in the calling script.
plot_raster_clean <- function(r, main, zlim = NULL, col, breaks = NULL,
                              legend = TRUE, categorical = FALSE, legend_title = "") {
  if (is.null(r) || all(is.na(terra::values(r)))) {
    plot.new(); text(0.5, 0.5, "No valid data", cex = 1.2, col = "red")
    return(invisible())
  }
  if (!categorical) r <- smooth_raster(r)
  plg <- list(cex = 0.9)
  if (nchar(legend_title)) plg$title <- legend_title
  if (!is.null(breaks)) {
    terra::plot(r, main = main, cex.main = 1.0, col = col, breaks = breaks,
                axes = FALSE, box = FALSE, legend = legend, plg = plg, colNA = NA)
  } else {
    terra::plot(r, main = main, cex.main = 1.0, col = col, zlim = zlim,
                axes = FALSE, box = FALSE, legend = legend, plg = plg, colNA = NA)
  }
  plot(sf::st_geometry(basin_sf_global),
       add = TRUE, col = NA, border = "black", lwd = 2.5)
}

## G. GGPLOT HELPERS ─────────────────────────────────────────────────
## Helper: make one annotation_custom rect band (avoids Date-scale numeric warning)
.band_rect <- function(ymin, ymax, fill, alpha = 0.10) {
  ggplot2::annotation_custom(
    grob = grid::rectGrob(
      gp = grid::gpar(fill = scales::alpha(fill, alpha), col = NA)),
    xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax)
}

drought_band_layers <- function() {
  list(
    .band_rect(-Inf,   -2.0,     "#7B0025"),   # D3: extremely dry
    .band_rect(-2.0,   -1.5,     "#D73027"),   # D2: severely dry
    .band_rect(-1.5,   -1.0,     "#F46D43"),   # D1: moderately dry
    .band_rect(-1.0,    1.0,     "#2E7D32"),   # N0: normal
    .band_rect( 1.0,    1.5,     "#90EE90"),   # W1: moderately wet
    .band_rect( 1.5,    2.0,     "#66BB6A"),   # W2: severely wet
    .band_rect( 2.0,    Inf,     "#4575B4"),   # W3: extremely wet
    ggplot2::geom_hline(yintercept = c(-2, -1.5, -1, 1, 1.5, 2),
                        linetype = "dotted", color = "gray50", linewidth = 0.3)
  )
}

shared_ts_theme <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

index_colours <- c(spi = "#e41a1c", spei = "#377eb8", swei = "#1E90FF")

calc_dynamic_ylim <- function(v, base = c(-4, 4), pad_frac = 0.10) {
  v <- v[is.finite(v)]
  if (!length(v)) return(base)
  r <- range(v)
  span <- diff(r); if (!is.finite(span) || span == 0) span <- 1
  c(min(base[1], r[1] - pad_frac * span),
    max(base[2], r[2] + pad_frac * span))
}

safe_pdf <- function(filename, width = 12, height = 8) {
  while (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
  tryCatch({
    grDevices::pdf(filename, width = width, height = height, family = "Helvetica")
    TRUE
  }, error = function(e) {
    cat("ERROR: Cannot open PDF '", filename, "': ", e$message, "\n")
    FALSE
  })
}

#' Unified figure-saving helper — replaces the scattered mix of pdf()/dev.off()
#' and ggsave() calls across the w* scripts.
save_figure <- function(plot_obj = NULL,
                        stem,
                        width  = FIG_WIDTH_WIDE,
                        height = FIG_HEIGHT_STD,
                        dpi    = FIG_DPI) {
  pdf_path <- paste0(stem, ".pdf")
  png_path <- paste0(stem, ".png")
  
  if (is.null(plot_obj)) {
    safe_pdf(pdf_path, width = width, height = height)
    return(invisible(NULL))
  }
  
  tryCatch(
    ggplot2::ggsave(pdf_path, plot_obj,
                    width = width, height = height, units = "in", device = "pdf"),
    error = function(e) cat(sprintf("  ⚠ PDF save failed (%s): %s\n",
                                    basename(pdf_path), e$message))
  )
  tryCatch(
    ggplot2::ggsave(png_path, plot_obj,
                    width = width, height = height, units = "in",
                    dpi = dpi, device = "png"),
    error = function(e) cat(sprintf("  ⚠ PNG save failed (%s): %s\n",
                                    basename(png_path), e$message))
  )
  cat(sprintf("  ✓ Figure saved: %s (.pdf + .png)\n", basename(stem)))
  invisible(NULL)
}

## H. EXCEL SUMMARY EXPORT ───────────────────────────────────────────
export_summary_excel <- function(timescales, index_types, ts_loader_fn, output_file) {
  wb <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  stats_df <- data.frame(Timescale = character(), Index = character(),
                         Mean = numeric(), Median = numeric(),
                         StdDev = numeric(), Min = numeric(), Max = numeric(),
                         Drought_Events_2020_2025 = integer(),
                         stringsAsFactors = FALSE)
  for (idx in index_types) {
    for (sc in timescales) {
      tryCatch({
        df  <- ts_loader_fn(idx, sc)
        df_rec  <- df[df$date >= as.Date("2020-01-01") & df$date <= as.Date("2025-12-31"), ]
        ev_count  <- nrow(detect_drought_events(df_rec))
        stats_df  <- rbind(stats_df, data.frame(
          Timescale = sprintf("%02d", sc),
          Index = toupper(idx),
          Mean = mean(df$value, na.rm = TRUE),
          Median = median(df$value, na.rm = TRUE),
          StdDev = sd(df$value, na.rm = TRUE),
          Min = min(df$value, na.rm = TRUE),
          Max = max(df$value, na.rm = TRUE),
          Drought_Events_2020_2025 = ev_count,
          stringsAsFactors = FALSE))
      }, error = function(e)
        cat(sprintf("  ⚠ Stats for %s-%02d: %s\n", toupper(idx), sc, e$message)))
    }
  }
  openxlsx::addWorksheet(wb, "Summary_Statistics")
  openxlsx::writeData(wb, "Summary_Statistics", stats_df)
  openxlsx::addStyle(wb, "Summary_Statistics", hdr_style,
                     rows = 1, cols = seq_len(ncol(stats_df)))
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  cat(sprintf("✓ Excel summary saved: %s\n", basename(output_file)))
  invisible(wb)
}

## I. W9 SHARED-STATE LOADER ─────────────────────────────────────────
read_w9_state <- function(path) {
  if (!file.exists(path))
    stop("w9 shared-state RDS not found:\n  ", path,
         "\n  Run w9_atmospheric_diagnostics.R first.")
  
  st <- readRDS(path)
  
  required_keys <- c(
    "base_date_index", "z500_ridge_ts", "slp_nwbc_ts", "sst_nepac_ts",
    "map_layers",
    "anomaly_nc_z500", "anomaly_nc_slp", "anomaly_nc_sst",
    "WD_PATH", "DATA_DIR", "OUT_DIR", "SPI_SEAS_DIR", "SPEI_SEAS_DIR",
    "START_YEAR", "END_YEAR", "CLIM_START", "CLIM_END",
    "RECENT_START", "RECENT_END",
    "DROUGHT_FOCUS_START", "DROUGHT_FOCUS_END",
    "RIDGE_LON_MIN", "RIDGE_LON_MAX", "RIDGE_LAT_MIN", "RIDGE_LAT_MAX",
    "SST_MEAN_LON_MIN", "SST_MEAN_LON_MAX", "SST_MEAN_LAT_MIN", "SST_MEAN_LAT_MAX",
    "EOF_LON_MIN", "EOF_LON_MAX", "EOF_LAT_MIN", "EOF_LAT_MAX",
    "FIGURE_DPI", "FIGURE_WIDTH_WIDE", "FIGURE_WIDTH_STD", "MONTH_LABELS"
  )
  
  missing_keys <- setdiff(required_keys, names(st))
  if (length(missing_keys))
    stop("w9 shared-state is missing ", length(missing_keys), " required key(s):\n",
         "  ", paste(missing_keys, collapse = ", "), "\n",
         "  The RDS schema may have changed since w9 last ran.\n",
         "  Re-run w9_atmospheric_diagnostics.R to regenerate the RDS.")
  
  cat(sprintf("✓ w9 shared state loaded (%d keys): %s\n",
              length(st), basename(path)))
  st
}

####################################################################################
# J. SHARED STATISTICAL FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────────
# These were previously duplicated between 6trend_test_ALL.R (vectorised pixel
# versions) and 4pr_pet_trends.r (single-series version).  They are now defined
# once here and sourced by both scripts.
#
# The public API is intentionally backward-compatible:
#   • mk_tfpw_spectral_for_series() — identical signature to the function that
#     was previously embedded in 4pr_pet_trends.r.  That script calls it
#     unchanged.
#   • vectorized_mann_kendall(), vectorized_mann_kendall_tfpw(),
#     vectorized_regime_shift_pelt(), vectorized_spectral_peaks() — identical
#     signatures to the functions previously in 6trend_test_ALL.R.  That script
#     now omits its local copies and relies on these definitions.
#   • detect_regime_shift_pelt() — single-series PELT helper, previously in
#     4pr_pet_trends.r under that name.  "trend_test_ALL.R" used a different but equivalent
#     vectorised wrapper; both now share this single-series core.
####################################################################################

# ── J0. Utility operator ────────────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# ─────────────────────────────────────────────────────────────────────────────────
# J1. PRIMITIVE HELPERS
# ─────────────────────────────────────────────────────────────────────────────────

#' Check whether too many values in a series are at or below a floor threshold.
#' Returns a list(exceeds_threshold, pct_min_vals).
#' Used by mk_tfpw_spectral_for_series() (formerly embedded in 4pr_pet_trends.r).
check_min_value_threshold <- function(ts_clean,
                                      min_val_threshold   = 0.01,
                                      max_pct             = 50,
                                      is_precip           = FALSE,
                                      max_min_value_pct_precip = 50,
                                      max_min_value_pct_pet    = 50) {
  n_min_vals   <- sum(ts_clean <= min_val_threshold, na.rm = TRUE)
  pct_min_vals <- n_min_vals / length(ts_clean) * 100
  max_allowed  <- if (is_precip) max_min_value_pct_precip else max_min_value_pct_pet
  list(exceeds_threshold = pct_min_vals > max_allowed, pct_min_vals = pct_min_vals)
}

#' Variance of Mann-Kendall S statistic corrected for ties.
#' Used by mk_tfpw_for_series() (formerly embedded in 4pr_pet_trends.r).
calculate_variance_with_ties <- function(S, n, x) {
  tie_table  <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  var_s      <- n * (n - 1) * (2 * n + 5) / 18
  if (length(tie_counts))
    var_s <- var_s - sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
  var_s
}

# ─────────────────────────────────────────────────────────────────────────────────
# J2. PELT CHANGEPOINT — single series
# ─────────────────────────────────────────────────────────────────────────────────
#' Detect the first mean-shift regime in a univariate time series using PELT/BIC.
#' Returns a named list with fields:
#'   changepoint_detected, changepoint_position, n_changepoints,
#'   first_changepoint_year, mean_before, mean_after, magnitude_shift.
#'
#' Previously defined as detect_regime_shift_pelt() in 4pr_pet_trends.r;
#' 6trend_test_ALL.R used an inline equivalent inside vectorized_regime_shift_pelt().
#' Both now share this function.
detect_regime_shift_pelt <- function(ts_vec, min_obs = 20, start_year = NULL) {
  ts_clean <- na.omit(ts_vec)
  n <- length(ts_clean)
  empty <- list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                n_changepoints = 0L, first_changepoint_year = NA_integer_,
                mean_before = NA_real_, mean_after = NA_real_, magnitude_shift = NA_real_)
  if (n < min_obs) return(empty)
  tryCatch({
    cpt      <- changepoint::cpt.mean(ts_clean, method = "PELT", penalty = "BIC")
    cpts_det <- changepoint::cpts(cpt)
    if (!length(cpts_det)) return(empty)
    fc <- cpts_det[1]
    list(changepoint_detected   = TRUE,
         changepoint_position   = fc,
         n_changepoints         = length(cpts_det),
         first_changepoint_year = as.integer(
           if (!is.null(start_year)) start_year + fc - 1L else fc),
         mean_before            = mean(ts_clean[seq_len(fc)], na.rm = TRUE),
         mean_after             = mean(ts_clean[(fc + 1L):n],  na.rm = TRUE),
         magnitude_shift        = mean(ts_clean[(fc + 1L):n],  na.rm = TRUE) -
           mean(ts_clean[seq_len(fc)],   na.rm = TRUE))
  }, error = function(e) empty)
}

# ─────────────────────────────────────────────────────────────────────────────────
# J3. VC-MK — single series (core used by both the vectorised wrapper and
#     mk_tfpw_spectral_for_series)
# ─────────────────────────────────────────────────────────────────────────────────
#' Hamed-Rao (1998) variance-corrected MK + Sen's slope for a single vector.
#' Returns list(tau, pval, slope, rho1, vc_corrected, filtered).
#' Internal helper; not normally called directly by analysis scripts.
.mk_vc_single <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 12L || var(x, na.rm = TRUE) < 1e-6)
    return(list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
                rho1 = NA_real_, vc_corrected = FALSE, filtered = TRUE))
  tryCatch({
    tau_val <- as.numeric(Kendall::MannKendall(x)$tau)
    hr      <- modifiedmk::mmkh(x)
    rho1    <- tryCatch({
      sl    <- as.numeric(trend::sens.slope(x)$estimates)
      det   <- x - sl * seq_along(x)
      acf(det, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]
    }, error = function(e) NA_real_)
    list(tau          = tau_val,
         pval         = as.numeric(hr["new P-value"]),
         slope        = as.numeric(hr["Sen's slope"]),
         rho1         = rho1,
         vc_corrected = !is.na(rho1) && abs(rho1) > 0.1,
         filtered     = FALSE)
  }, error = function(e)
    list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
         rho1 = NA_real_, vc_corrected = FALSE, filtered = TRUE))
}

# ─────────────────────────────────────────────────────────────────────────────────
# J4. TFPW-MK — single series
# ─────────────────────────────────────────────────────────────────────────────────
#' TFPW (Trend-Free Pre-Whitening, Yue et al. 2002) MK for a single vector.
#' Sen's slope is taken from the ORIGINAL series (not the pre-whitened one).
#' Returns list(tau, pval, slope, rho1, tfpw_applied, filtered).
#' Internal helper; not normally called directly by analysis scripts.
.mk_tfpw_single <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 12L || var(x, na.rm = TRUE) < 1e-6)
    return(list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
                rho1 = NA_real_, tfpw_applied = FALSE, filtered = TRUE))
  tryCatch({
    # Step 1-2: Sen's slope on original series
    beta <- as.numeric(trend::sens.slope(x)$estimates)
    if (!is.finite(beta))
      return(list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
                  rho1 = NA_real_, tfpw_applied = FALSE, filtered = TRUE))
    
    # Step 3: Detrend → lag-1 autocorrelation
    t_seq <- seq_len(n)
    x_det <- x - beta * t_seq
    r1    <- cor(x_det[-n], x_det[-1], use = "complete.obs")
    if (is.na(r1)) r1 <- 0
    r1 <- max(-0.99, min(0.99, r1))
    
    # Step 3b: Pre-whiten + blend original trend back in
    x_pw    <- x_det[-1] - r1 * x_det[-n]
    x_final <- x_pw + beta * t_seq[-1]
    if (length(x_final) < 10L)
      return(list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
                  rho1 = r1, tfpw_applied = FALSE, filtered = TRUE))
    
    mk <- Kendall::MannKendall(x_final)
    list(tau          = as.numeric(mk$tau),
         pval         = as.numeric(mk$sl),
         slope        = beta,    # ← original-series Sen's slope (Yue et al. 2002)
         rho1         = r1,
         tfpw_applied = abs(r1) > 0.1,
         filtered     = FALSE)
  }, error = function(e)
    list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
         rho1 = NA_real_, tfpw_applied = FALSE, filtered = TRUE))
}

# ─────────────────────────────────────────────────────────────────────────────────
# J5. SPECTRAL PEAK DETECTION — single series
# ─────────────────────────────────────────────────────────────────────────────────
#' AR(1) red-noise spectral test for a single detrended time series.
#' Counts contiguous groups of periodogram bins above the chi-squared threshold
#' at periods >= min_period_months (default 24 months).
#' Returns list(n_peaks, dominant_period, conf).
#' Internal helper; not normally called directly by analysis scripts.
.spectral_single <- function(x,
                             min_period_months = 24L,
                             sig_level         = 0.95,
                             n_sim             = 500L,
                             alpha             = 0.05,
                             conf_cache_env    = NULL) {
  x <- na.omit(x)
  n <- length(x)
  empty <- list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_)
  if (n < 48L) return(empty)
  
  # Demean and linearly detrend
  t_s  <- seq_len(n)
  cf   <- lm.fit(cbind(1, t_s), x)$coefficients
  x_dt <- x - (cf[1] + cf[2] * t_s)
  
  tryCatch({
    # ── Method 1: AR(1) red-noise via smoothed periodogram (used by 6trend_test_ALL) ──────
    spec_r <- stats::spectrum(x_dt, spans = c(3L, 5L), taper = 0.1,
                              plot = FALSE, na.action = na.omit)
    freq   <- spec_r$freq
    spower <- spec_r$spec
    period <- 1.0 / freq
    keep   <- is.finite(period) & period >= min_period_months
    if (!any(keep)) return(empty)
    
    f_k <- freq[keep]; p_k <- spower[keep]
    r1  <- tryCatch(stats::cor(x_dt[-n], x_dt[-1], use = "complete.obs"),
                    error = function(e) 0.0)
    if (is.na(r1) || !is.finite(r1)) r1 <- 0.0
    r1  <- max(-0.99, min(0.99, r1))
    rn_raw <- (1.0 - r1^2) / (1.0 - 2.0 * r1 * cos(2.0 * pi * f_k) + r1^2)
    rn_scl <- rn_raw * (mean(p_k) / mean(rn_raw))
    dof    <- if (!is.null(spec_r$df)) spec_r$df else 2.0
    thr_ar1 <- rn_scl * stats::qchisq(sig_level, dof) / dof
    
    above_ar1 <- p_k > thr_ar1
    n_peaks_ar1 <- if (any(above_ar1)) as.integer(sum(rle(above_ar1)$values)) else 0L
    dom_period  <- if (n_peaks_ar1 > 0L) {
      pk_idx <- which(above_ar1)
      1.0 / f_k[pk_idx[which.max(p_k[pk_idx])]]
    } else NA_real_
    
    # ── Method 2: Monte-Carlo white-noise envelope (used by 4pr_pet_trends.r) ─
    # Also compute this for the omnibus wrapper so it can serve both scripts.
    conf_limit <- NA_real_
    if (!is.null(conf_cache_env)) {
      key <- paste0("n_", n)
      if (!exists(key, envir = conf_cache_env, inherits = FALSE)) {
        half        <- floor(n / 2)
        max_spectra <- vapply(seq_len(n_sim), function(jj) {
          r  <- rnorm(n, 0, 1)
          fr <- fft(r - mean(r))
          max(Mod(fr[seq_len(half)])^2 / n, na.rm = TRUE)
        }, numeric(1))
        assign(key, stats::quantile(max_spectra, 1 - alpha, na.rm = TRUE),
               envir = conf_cache_env)
      }
      conf_limit <- get(key, envir = conf_cache_env, inherits = FALSE)
      
      # Re-check with MC envelope (used by 4pr_pet_trends.r path)
      sd_d    <- stats::sd(x_dt); if (!is.finite(sd_d) || sd_d == 0) sd_d <- 1
      dstd    <- (x_dt - mean(x_dt)) / sd_d
      half    <- floor(n / 2)
      ff      <- fft(dstd)
      sp_mc   <- Mod(ff[seq_len(half)])^2 / n
      freqs_mc <- seq(0, 0.5, length.out = half)
      pk_mc   <- which(sp_mc > conf_limit)
      pp_mc   <- if (length(pk_mc)) {
        pf <- freqs_mc[pk_mc]; pp <- ifelse(pf > 0, 1 / pf, NA_real_)
        pp[order(sp_mc[pk_mc], decreasing = TRUE)]
      } else numeric(0)
      return(list(n_peaks        = length(pp_mc),
                  dominant_period = if (length(pp_mc)) pp_mc[1] else NA_real_,
                  conf            = conf_limit))
    }
    
    list(n_peaks         = n_peaks_ar1,
         dominant_period = dom_period,
         conf            = NA_real_)
    
  }, error = function(e) empty)
}

# ─────────────────────────────────────────────────────────────────────────────────
# J6. OMNIBUS SINGLE-SERIES WRAPPER
#     mk_tfpw_spectral_for_series()
#     Used directly by 4pr_pet_trends.r (and its visualization script reads
#     outputs produced by it).  Signature is identical to the version that
#     was previously embedded in that script.
# ─────────────────────────────────────────────────────────────────────────────────
#' Full VC-MK + TFPW-MK + spectral + PELT analysis for a single time series.
#'
#' @param ts_vec           Numeric vector (may contain NAs).
#' @param is_precip        Logical; TRUE → precipitation quality filters.
#' @param alpha            Significance level (default 0.05).
#' @param max_tie_pct      Maximum % ties before filtering (default 50).
#' @param n_sim_spectral   Monte-Carlo replicates for spectral null (default 500).
#' @param conf_cache_env   Environment for caching MC spectral thresholds.
#' @param start_year       First year of ts_vec (for PELT year labelling).
#' @param min_nonzero_count  Minimum count of nonzero values (PET monthly filter).
#' @param skip_min_filter  Logical; bypass the min-value quality filter.
#' @param max_min_value_pct_precip  Max % floor values for precipitation.
#' @param max_min_value_pct_pet     Max % floor values for PET.
#' @param min_positive_value        Floor value below which a value is "zero".
#' @param min_obs_changepoint       Minimum obs for PELT (default 20).
#'
#' @return Named list with elements \code{vc}, \code{tf}, \code{spec}, \code{cpt}.
mk_tfpw_spectral_for_series <- function(ts_vec,
                                        is_precip            = FALSE,
                                        alpha                = 0.05,
                                        max_tie_pct          = 50,
                                        n_sim_spectral       = 500L,
                                        conf_cache_env       = NULL,
                                        start_year           = NULL,
                                        min_nonzero_count    = NULL,
                                        skip_min_filter      = FALSE,
                                        max_min_value_pct_precip = 50,
                                        max_min_value_pct_pet    = 50,
                                        min_positive_value   = 0.01,
                                        min_obs_changepoint  = 20L) {
  ts_clean <- na.omit(ts_vec)
  n        <- length(ts_clean)
  
  na_result <- function(reason, n_ties = 0, pct_ties = 0, pct_min = 0) {
    list(
      vc  = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_,
                 varS = NA_real_, n = n, rho1 = NA_real_, vc_corrected = FALSE,
                 n_ties = n_ties, percent_ties = pct_ties, n_min_vals = 0L,
                 percent_min_vals = pct_min, tau_b_adjusted = FALSE,
                 filtered = TRUE, reason = reason),
      tf  = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_,
                 varS = NA_real_, n = n, rho1 = NA_real_, tfpw_applied = FALSE,
                 n_ties = n_ties, percent_ties = pct_ties, n_min_vals = 0L,
                 percent_min_vals = pct_min, tau_b_adjusted = FALSE,
                 filtered = TRUE, reason = reason),
      spec = list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_),
      cpt  = list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                  n_changepoints = 0L, first_changepoint_year = NA_integer_,
                  mean_before = NA_real_, mean_after = NA_real_,
                  magnitude_shift = NA_real_))
  }
  
  if (n < 10L) return(na_result("low_n"))
  
  # ── Quality filters ──────────────────────────────────────────────────────────
  if (skip_min_filter) {
    min_check <- list(exceeds_threshold = FALSE, pct_min_vals = 0)
  } else if (!is.null(min_nonzero_count)) {
    n_nonzero    <- sum(ts_clean > min_positive_value, na.rm = TRUE)
    pct_min_vals <- (1 - n_nonzero / n) * 100
    if (n_nonzero < min_nonzero_count)
      return(na_result("insufficient_nonzero", pct_min = pct_min_vals))
    min_check <- list(exceeds_threshold = FALSE, pct_min_vals = pct_min_vals)
  } else {
    min_val_threshold <- if (is_precip) 0.0 else min_positive_value
    min_check <- check_min_value_threshold(
      ts_clean, min_val_threshold,
      max_pct = if (is_precip) max_min_value_pct_precip else max_min_value_pct_pet,
      is_precip = is_precip,
      max_min_value_pct_precip = max_min_value_pct_precip,
      max_min_value_pct_pet    = max_min_value_pct_pet)
    if (min_check$exceeds_threshold)
      return(na_result("excessive_min_vals", pct_min = min_check$pct_min_vals))
  }
  
  n_unique     <- length(unique(ts_clean))
  n_ties       <- n - n_unique
  percent_ties <- n_ties / n * 100
  if (percent_ties > max_tie_pct)
    return(na_result("excessive_ties", n_ties, percent_ties,
                     min_check$pct_min_vals))
  
  tau_b_adjusted <- (percent_ties > 5)
  
  # ── Sen's slope (on original series, shared by both MK variants) ────────────
  sen_slope <- trend::sens.slope(ts_clean)$estimates
  if (is.na(sen_slope) || is.infinite(sen_slope) || length(sen_slope) == 0)
    return(na_result("sens_slope_fail", n_ties, percent_ties,
                     min_check$pct_min_vals))
  
  # ── Detrend → lag-1 autocorrelation ─────────────────────────────────────────
  time_index <- seq_len(n)
  trend_line <- sen_slope * time_index
  detrended  <- ts_clean - trend_line
  acf_result <- tryCatch(
    acf(detrended, lag.max = 1, plot = FALSE, na.action = na.pass),
    error = function(e) NULL, warning = function(w) NULL)
  rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA_real_
  
  # ── VC-MK (Hamed-Rao 1998) ───────────────────────────────────────────────────
  # Use mmkh() — the same function used by .mk_vc_single() (J3) and
  # vectorized_mann_kendall() (J7a) — so all three code paths produce
  # consistent H-R98 variance-corrected p-values.
  # mk3() is a different modified-MK variant and must NOT be used here.
  vc_corrected <- FALSE
  p_vc <- NA_real_; S_vc <- NA_real_; tau_vc <- NA_real_; varS_final <- NA_real_
  vc_res <- tryCatch(modifiedmk::mmkh(ts_clean), error = function(e) NULL)
  if (!is.null(vc_res)) {
    tau_vc     <- as.numeric(vc_res["Kendall's tau"])
    p_vc       <- as.numeric(vc_res["new P-value"])
    S_vc       <- as.numeric(vc_res["S"])
    varS_final <- as.numeric(vc_res["variance of S"])
    if (!is.na(rho1) && abs(rho1) > 0.1) vc_corrected <- TRUE
  }
  vc_list <- list(tau = tau_vc, p = p_vc, sl = sen_slope, S = S_vc,
                  varS = varS_final, n = n, rho1 = rho1,
                  vc_corrected = vc_corrected,
                  n_ties = n_ties, percent_ties = percent_ties,
                  n_min_vals = 0L, percent_min_vals = min_check$pct_min_vals,
                  tau_b_adjusted = tau_b_adjusted, filtered = FALSE, reason = "none")
  
  # ── TFPW-MK ─────────────────────────────────────────────────────────────────
  tfpw_applied <- FALSE
  if (!is.na(rho1) && abs(rho1) > 0.1) {
    pw    <- numeric(n); pw[1] <- detrended[1]
    for (j in 2:n) pw[j] <- detrended[j] - rho1 * detrended[j - 1]
    corrected    <- pw + trend_line
    tfpw_applied <- TRUE
  } else {
    corrected <- ts_clean
  }
  s_mat_tf <- sign(outer(corrected, corrected, `-`))
  S_tf     <- sum(s_mat_tf[upper.tri(s_mat_tf)], na.rm = TRUE)
  varS_tf  <- calculate_variance_with_ties(S_tf, n, corrected)
  n_pairs  <- n * (n - 1) / 2
  tau_tf   <- S_tf / n_pairs
  p_tf     <- if (varS_tf <= 0) NA_real_ else 2 * pnorm(-abs(S_tf / sqrt(varS_tf)))
  tf_list  <- list(tau = tau_tf, p = p_tf, sl = sen_slope, S = S_tf,
                   varS = varS_tf, n = n, rho1 = rho1,
                   tfpw_applied = tfpw_applied,
                   n_ties = n_ties, percent_ties = percent_ties,
                   n_min_vals = 0L, percent_min_vals = min_check$pct_min_vals,
                   tau_b_adjusted = tau_b_adjusted, filtered = FALSE, reason = "none")
  
  # ── Spectral ────────────────────────────────────────────────────────────────
  spec_result <- .spectral_single(ts_clean,
                                  min_period_months = 24L,
                                  sig_level         = 0.95,
                                  n_sim             = n_sim_spectral,
                                  alpha             = alpha,
                                  conf_cache_env    = conf_cache_env)
  
  # ── PELT ────────────────────────────────────────────────────────────────────
  cpt_result <- detect_regime_shift_pelt(ts_clean,
                                         min_obs    = min_obs_changepoint,
                                         start_year = start_year)
  
  list(vc = vc_list, tf = tf_list, spec = spec_result, cpt = cpt_result)
}

# ─────────────────────────────────────────────────────────────────────────────────
# J7. VECTORISED PIXEL-MATRIX WRAPPERS
#     Used directly by 6trend_test_ALL.R.  Signatures are identical to the
#     functions previously defined in that script.
#     The parallel back-end variables N_CORES, is_windows, cl are expected to
#     be set in the calling script (6trend_test_ALL.R) before these are called.
# ─────────────────────────────────────────────────────────────────────────────────

## J7a. Hamed-Rao (1998) Variance-Corrected MK + Sen's slope ──────────────────
vectorized_mann_kendall <- function(ts_matrix) {
  n_pix <- nrow(ts_matrix)
  
  process_pixel <- function(i) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 12L || var(x, na.rm = TRUE) < 1e-6)
      return(list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
                  filtered = TRUE))
    tryCatch({
      tau_val <- as.numeric(Kendall::MannKendall(x)$tau)
      hr      <- modifiedmk::mmkh(x)
      list(tau      = tau_val,
           pval     = as.numeric(hr["new P-value"]),
           slope    = as.numeric(hr["Sen's slope"]),
           filtered = FALSE)
    }, error = function(e)
      list(tau = NA_real_, pval = NA_real_, slope = NA_real_, filtered = TRUE))
  }
  
  do_par <- N_CORES > 1L && n_pix > 100L
  res <- if (do_par && is_windows) {
    parallel::parLapply(cl, seq_len(n_pix), process_pixel)
  } else if (do_par) {
    parallel::mclapply(seq_len(n_pix), process_pixel, mc.cores = N_CORES)
  } else {
    lapply(seq_len(n_pix), process_pixel)
  }
  
  data.table::data.table(
    tau_vc      = sapply(res, `[[`, "tau"),
    p_value_vc  = sapply(res, `[[`, "pval"),
    sl_vc       = sapply(res, `[[`, "slope"),
    filtered_vc = sapply(res, `[[`, "filtered")
  )
}

## J7b. TFPW Mann-Kendall (Yue et al. 2002) ───────────────────────────────────
vectorized_mann_kendall_tfpw <- function(ts_matrix) {
  n_pix <- nrow(ts_matrix)
  
  process_pixel_tfpw <- function(i) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 12L || var(x, na.rm = TRUE) < 1e-6)
      return(list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
                  filtered = TRUE))
    tryCatch({
      res <- .mk_tfpw_single(x)
      list(tau      = res$tau,
           pval     = res$pval,
           slope    = res$slope,
           filtered = res$filtered)
    }, error = function(e)
      list(tau = NA_real_, pval = NA_real_, slope = NA_real_, filtered = TRUE))
  }
  
  do_par <- N_CORES > 1L && n_pix > 100L
  res <- if (do_par && is_windows) {
    parallel::parLapply(cl, seq_len(n_pix), process_pixel_tfpw)
  } else if (do_par) {
    parallel::mclapply(seq_len(n_pix), process_pixel_tfpw, mc.cores = N_CORES)
  } else {
    lapply(seq_len(n_pix), process_pixel_tfpw)
  }
  
  data.table::data.table(
    tau_tfpw      = sapply(res, `[[`, "tau"),
    p_value_tfpw  = sapply(res, `[[`, "pval"),
    sl_tfpw       = sapply(res, `[[`, "slope"),
    filtered_tfpw = sapply(res, `[[`, "filtered")
  )
}

## J7c. PELT regime-shift — pixel matrix ──────────────────────────────────────
#' Applies detect_regime_shift_pelt() to annual averages of each pixel's
#' monthly time series.
#' NOTE: regime_shift_year is kept as the primary column name so that
#'       4_trends_visualization.R (Fig 3, Panel B) works without modification.
vectorized_regime_shift_pelt <- function(ts_matrix, years, min_obs = 20L) {
  n_pix <- nrow(ts_matrix)
  n_yrs <- length(years)
  
  out <- data.table::data.table(
    regime_shift_year     = rep(NA_real_,  n_pix),
    regime_shift_detected = rep(FALSE,     n_pix),
    n_changepoints        = rep(0L,        n_pix),
    mean_before_shift     = rep(NA_real_,  n_pix),
    mean_after_shift      = rep(NA_real_,  n_pix),
    magnitude_shift       = rep(NA_real_,  n_pix)
  )
  
  for (i in seq_len(n_pix)) {
    row <- ts_matrix[i, ]
    
    ann <- vapply(seq_len(n_yrs), function(k) {
      cols <- ((k - 1L) * 12L + 1L):(k * 12L)
      vals <- row[cols[cols <= length(row)]]
      if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
    }, numeric(1L))
    
    ok <- !is.na(ann)
    if (sum(ok) < min_obs) next
    ann_clean <- ann[ok]
    yrs_clean <- years[ok]
    
    res <- detect_regime_shift_pelt(ann_clean,
                                    min_obs    = min_obs,
                                    start_year = NULL)   # year looked up below
    if (!res$changepoint_detected) next
    
    fc <- res$changepoint_position
    data.table::set(out, i, "regime_shift_year",     as.numeric(yrs_clean[fc]))
    data.table::set(out, i, "regime_shift_detected", TRUE)
    data.table::set(out, i, "n_changepoints",        res$n_changepoints)
    data.table::set(out, i, "mean_before_shift",     res$mean_before)
    data.table::set(out, i, "mean_after_shift",      res$mean_after)
    data.table::set(out, i, "magnitude_shift",       res$magnitude_shift)
  }
  out
}

## J7d. Spectral peak detection — pixel matrix ─────────────────────────────────
vectorized_spectral_peaks <- function(ts_matrix,
                                      min_period_months = 24L,
                                      sig_level         = 0.95) {
  n_pix <- nrow(ts_matrix)
  
  process_pixel_spec <- function(i) {
    x <- ts_matrix[i, ]
    res <- .spectral_single(x,
                            min_period_months = min_period_months,
                            sig_level         = sig_level,
                            conf_cache_env    = NULL)   # AR1 path, no MC cache needed
    res$n_peaks
  }
  
  do_par <- N_CORES > 1L && n_pix > 100L
  res <- if (do_par && is_windows) {
    parallel::parLapply(cl, seq_len(n_pix), process_pixel_spec)
  } else if (do_par) {
    parallel::mclapply(seq_len(n_pix), process_pixel_spec, mc.cores = N_CORES)
  } else {
    lapply(seq_len(n_pix), process_pixel_spec)
  }
  
  data.table::data.table(n_spectral_peaks = as.integer(unlist(res)))
}

#source("utils_teleconnection_addon.R")
cat("✓ DROUGHT_ANALYSIS_utils.R loaded\n")