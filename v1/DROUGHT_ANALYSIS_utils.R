####################################################################################
# DROUGHT_ANALYSIS_utils_v2.R  ??  SHARED UTILITIES FOR ALL NECHAKO DROUGHT SCRIPTS
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# Sourced by:
#   6trend_test_ALL_v2.R              ??? drought index trend analysis
#   3SPEI_ERALand_v2.R                ??? SPI/SPEI-PM/SPEI-Thw computation
#   4pr_pet_trends_v2.r               ??? Pr / PET / Temperature trend analysis
#   9event_ranking.R               ??? drought event ranking and return periods
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
#   J.  [SHARED STATS]  Functions used by BOTH 6trend_test_ALL_v2.R and
#       4pr_pet_trends_v2.r ??? moved here to eliminate duplication:
#         J1. Utility helpers  (check_min_value_threshold,
#                               calculate_variance_with_ties)
#         J2. PELT changepoint  (detect_regime_shift_pelt)
#         J3. VC-MK per-series  (.mk_vc_single)
#         J4. TFPW-MK per-series (.mk_tfpw_single)
#         J5. Spectral peak detection per-series (.spectral_single)
#         J6. Omnibus single-series wrapper (mk_tfpw_spectral_for_series)
#             ??? used directly by 4pr_pet_trends_v2.r
#         J7. Vectorised pixel-matrix wrappers  (vectorized_mann_kendall,
#             vectorized_mann_kendall_tfpw, vectorized_regime_shift_pelt,
#             vectorized_spectral_peaks)
#             ??? used directly by 6trend_test_ALL_v2.R
####################################################################################

## A. PROJECT CONFIGURATION ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
WD_PATH <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")
BASIN_SHP       <- file.path(WD_PATH, "Spatial/Nechako_shps")

## Helper: load the authoritative Nechako basin as a single dissolved
## SpatVector from the supplied subbasin Shapefiles.  The source directory is
## intentionally used instead of a legacy KMZ/KML boundary.
load_basin_vect <- function(path = BASIN_SHP) {
  if (!dir.exists(path))
    stop("Nechako subbasin Shapefile directory not found: ", path)
  nechako_subbasin_union(load_nechako_subbasins(path))
}

## Index-specific directories
SPI_SEAS_DIR      <- file.path(WD_PATH, "spi_results_seasonal/")
SPEI_SEAS_DIR     <- file.path(WD_PATH, "spei_results_seasonal/")
SPEI_ERA5_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal_era5land/")   # ERA5-Land native PEV outputs
SWEI_SEAS_DIR     <- file.path(WD_PATH, "swei_results_seasonal/")
SSPI_SEAS_DIR     <- file.path(WD_PATH, "sspi_results_seasonal/")

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
## Drought-event thresholds (USER SETTINGS; used by every script).
## Onset: index < DROUGHT_ONSET. Termination: index >= DROUGHT_END.
## DROUGHT_END must be >= DROUGHT_ONSET; equal values = closed hysteresis band.
## Minimum event durations per scale: EVENT_MIN_DURATION_BY_SCALE (section F).
DROUGHT_ONSET     <- -0.5
DROUGHT_END       <- -0.5

## Shared basin-mask / basin-mean contract (applies to every script)
##   * Basin geometry : dissolved union of the 7 authoritative subbasin
##                      Shapefiles (load_basin_vect() == nechako_subbasin_union()).
##   * Basin mask     : native-grid cells with any fractional coverage (> 0).
##   * Basin mean     : sum(x * w) / sum(w over finite cells), where
##                      w = cell area (m2) x fractional coverage of the cell.
##   * Reporting rule : a mean is reported only when the finite cells hold at
##                      least BASIN_MIN_VALID_AREA_COVERAGE of the region's total
##                      weight (fraction of area, 0-1); otherwise NA.
## The same rule is used for the whole basin and for every subbasin.
BASIN_MIN_VALID_AREA_COVERAGE <- 0.20


# ==============================================================================
# CANONICAL NATIVE ERA5-LAND GRID CONTRACT
# ==============================================================================
COMMON_GRID_FILE <- file.path(WD_PATH, "era5land_common_analysis_grid.rds")

save_common_analysis_grid <- function(template, basin = NULL,
                                      source = "ERA5-Land native precipitation grid") {
  if (!inherits(template, "SpatRaster")) stop("template must be a SpatRaster.")
  template <- template[[1]]
  if (is.na(terra::crs(template)) || !nzchar(terra::crs(template)))
    stop("Canonical template has no CRS.")
  if (is.null(basin)) basin <- load_basin_vect(BASIN_SHP)
  if (!inherits(basin, "SpatVector")) basin <- terra::vect(basin)
  if (!terra::same.crs(basin, template)) basin <- terra::project(basin, terra::crs(template))
  cover <- terra::rasterize(basin, template, cover = TRUE)
  cover_v <- as.numeric(terra::values(cover, mat = FALSE))
  cover_v[!is.finite(cover_v)] <- 0
  cover_v <- pmax(0, pmin(1, cover_v))
  mask_r <- terra::setValues(template, as.numeric(cover_v > 0))
  coverage_r <- terra::setValues(template, cover_v)
  grid <- list(
    raster_template = terra::wrap(template), basin_mask = terra::wrap(mask_r),
    basin_coverage = terra::wrap(coverage_r), basin_pixels = sum(cover_v > 0),
    total_pixels = terra::ncell(template), crs = terra::crs(template),
    resolution = terra::res(template), extent = terra::ext(template),
    origin = terra::origin(template), nrow = terra::nrow(template),
    ncol = terra::ncol(template), source = source, created = Sys.time()
  )
  saveRDS(grid, COMMON_GRID_FILE, compress = TRUE)
  grid
}

load_common_analysis_grid <- function(required = TRUE) {
  if (!file.exists(COMMON_GRID_FILE)) {
    if (required) stop("Canonical ERA5-Land grid not found: ", COMMON_GRID_FILE)
    return(NULL)
  }
  g <- readRDS(COMMON_GRID_FILE)
  if (!is.list(g) || is.null(g$raster_template) || is.null(g$basin_mask))
    stop("Invalid canonical ERA5-Land grid file: ", COMMON_GRID_FILE)
  g$raster_template <- terra::unwrap(g$raster_template)
  g$basin_mask <- terra::unwrap(g$basin_mask)
  if (!is.null(g$basin_coverage)) g$basin_coverage <- terra::unwrap(g$basin_coverage)
  g
}

validate_common_analysis_grid <- function(x, grid = NULL, label = "raster") {
  if (!inherits(x, "SpatRaster")) stop(label, " must be a SpatRaster.")
  if (is.null(grid)) grid <- load_common_analysis_grid()
  geom_ok <- terra::compareGeom(
    x[[1]],
    grid$raster_template,
    stopOnError = FALSE,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE
  )
  orig_ok <- isTRUE(all.equal(
    terra::origin(x[[1]]),
    terra::origin(grid$raster_template),
    tolerance = 1e-6
  ))
  if (!isTRUE(geom_ok) || !isTRUE(orig_ok))
    stop(label, " does not match the canonical native ERA5-Land grid.")
  invisible(TRUE)
}

## Timescales ??? current manuscript source of truth
MANUSCRIPT_SCALES    <- c(1L, 3L, 6L, 12L, 24L, 36L)
EVENT_RANKING_SCALES <- MANUSCRIPT_SCALES
SPATIAL_PLOT_SCALES  <- c(1L, 3L, 6L, 12L, 24L, 36L)

SPI_SCALES          <- MANUSCRIPT_SCALES
SPEI_SCALES         <- MANUSCRIPT_SCALES
SWEI_SCALE          <- 3L
INDICES             <- c("spi", "spei", "spei_era5", "swei")
TIMESCALES_STANDARD <- MANUSCRIPT_SCALES
TIMESCALES_SPATIAL  <- SPATIAL_PLOT_SCALES
MONTH_NAMES         <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

## A2. VALIDATION HELPERS
validate_scales <- function(observed, expected = MANUSCRIPT_SCALES,
                            allow_extra = TRUE, label = "scales") {
  observed <- sort(unique(as.integer(observed)))
  expected <- sort(unique(as.integer(expected)))
  if (anyNA(observed)) stop(label, " contains invalid scale values.")
  missing <- setdiff(expected, observed)
  extra <- setdiff(observed, expected)
  if (length(missing))
    stop(label, " missing required scale(s): ", paste(missing, collapse = ", "))
  if (length(extra) && !allow_extra)
    stop(label, " contains non-manuscript scale(s): ", paste(extra, collapse = ", "))
  if (length(extra))
    warning(label, " contains additional scale(s): ", paste(extra, collapse = ", "),
            call. = FALSE)
  invisible(TRUE)
}

validate_monthly_record <- function(dates, start_date = NULL, end_date = NULL,
                                    label = "monthly record",
                                    require_complete = TRUE) {
  dates <- as.Date(dates)
  if (!length(dates)) stop(label, " has no dates.")
  if (anyNA(dates)) stop(label, " contains NA/unparseable dates.")
  if (anyDuplicated(dates)) stop(label, " contains duplicated month(s).")
  if (is.unsorted(dates)) stop(label, " is not chronologically sorted.")
  if (any(format(dates, "%d") != "01"))
    stop(label, " contains dates that are not month starts.")
  
  expected_internal <- seq(min(dates), max(dates), by = "month")
  missing_internal <- setdiff(expected_internal, dates)
  if (require_complete && length(missing_internal))
    stop(label, " has missing calendar month(s): ",
         paste(missing_internal, collapse = ", "))
  
  if (!is.null(start_date) || !is.null(end_date)) {
    if (is.null(start_date) || is.null(end_date))
      stop("Supply start_date and end_date together.")
    expected <- seq(as.Date(start_date), as.Date(end_date), by = "month")
    missing <- setdiff(expected, dates)
    extra <- setdiff(dates, expected)
    if (require_complete && length(missing))
      stop(label, " incomplete for requested period; missing: ",
           paste(missing, collapse = ", "))
    if (length(extra))
      stop(label, " contains dates outside requested period.")
  }
  invisible(TRUE)
}

validate_standardized_index <- function(x, label = "standardized index",
                                        absolute_limit = 10) {
  if (!is.numeric(x)) stop(label, " must be numeric.")
  xf <- x[is.finite(x)]
  if (!length(xf)) stop(label, " contains no finite values.")
  if (any(abs(xf) > absolute_limit))
    stop(label, " contains value(s) outside +/-", absolute_limit, ".")
  invisible(TRUE)
}

validate_event_catalog <- function(events, min_duration = 1L,
                                   label = "event catalog") {
  if (!is.data.frame(events)) stop(label, " must be a data.frame.")
  if (!nrow(events)) return(invisible(TRUE))
  req <- c("start_date", "end_date", "duration_months", "mean_int", "severity")
  miss <- setdiff(req, names(events))
  if (length(miss)) stop(label, " missing: ", paste(miss, collapse = ", "))
  if (any(as.Date(events$end_date) < as.Date(events$start_date)))
    stop(label, " contains end_date < start_date.")
  if (any(events$duration_months < as.integer(min_duration), na.rm = TRUE))
    stop(label, " contains event(s) shorter than min_duration.")
  if (any(!is.finite(events$mean_int)) || any(events$mean_int < 0))
    stop(label, " contains invalid intensity.")
  if (any(!is.finite(events$severity)) || any(events$severity < 0))
    stop(label, " contains invalid severity.")
  invisible(TRUE)
}

## B. BASIN-AVERAGED CSV LOADER ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
#' Load a pre-computed basin-averaged drought index CSV and return a tidy
#' long-format data.frame(date, value).
#'
#' AUTHORITATIVE SOURCE: these CSVs are produced by 6trend_test_ALL_v2.R via
#' save_basin_avg_from_pixels() (see Section B2 below).  The methodology is:
#'   (1) compute the drought index at every basin pixel from per-pixel inputs,
#'   (2) take the area-weighted spatial mean of those pixel-level index values.
#' This is the standard hydrometeorological approach: area-weighted averaging
#' of standardised index values that were calibrated on the same spatial basis.
#' Run 6trend_test_ALL_v2.R before basin_timeseries.R.
#'
#' File naming convention: {index}_{scale:02d}_basin_averaged_by_month.csv
#'   where index may be spi, spei (FAO-56 PM), spei_era5 (ERA5-Land native PEV), or swei
#' File structure:
#'   Year | Jan | Feb | Mar | ??? | Dec
#'   1950 | -1.88 |  NA | ???
#'
#' @param index_type  "spi", "spei", "spei_era5", or "swei"
#' @param scale       Integer timescale (e.g. 1, 3, 6, 12)
#' @param data_dir    Directory containing the CSV; if NULL, the index-specific
#'                    seasonal directory from section A is used automatically.
#' @return data.frame with columns \code{date} (Date) and \code{value} (numeric),
#'         sorted by date with missing values retained as NA. Returns NULL with a warning if
#'         the file is not found (usually means w3 has not been run yet).
load_basin_avg_csv <- function(index_type, scale, data_dir = NULL) {
  dir_map <- list(spi      = SPI_SEAS_DIR,
                  spei     = SPEI_SEAS_DIR,
                  spei_era5 = SPEI_ERA5_SEAS_DIR,
                  swei     = SWEI_SEAS_DIR,
                  sspi     = SSPI_SEAS_DIR)
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
  long$value[!is.finite(long$value)] <- NA_real_
  rownames(long) <- NULL
  long
}

## B2. BASIN-AVERAGED CSV WRITER ????????????????????????????????????????????????????????????????????????????????????????????????????????????
#' Compute an area-weighted basin-average time series from a pixel ?? month
#' matrix and write the result as the standard basin_averaged_by_month.csv
#' expected by load_basin_avg_csv().
#'
#' Called by 6trend_test_ALL_v2.R::process_index() after the pixel time-series
#' matrix has been assembled.  The area weights are the equal-area cell sizes
#' extracted from the raster template (already in EPSG:3005).
#'
#' Methodology:
#'   basin_avg[t] = ??( index_value[pixel, t] ?? cell_area[pixel] )
#'                / ??( cell_area[pixel] )          (sum over non-NA pixels)
#'
#' @param ts_matrix   Numeric matrix (n_valid_pixels ?? n_months).  Rows
#'                    correspond to valid basin pixels; columns to calendar
#'                    months in chronological order.
#' @param area_weights Numeric vector of length n_valid_pixels; cell areas in
#'                    m?? in the equal-area CRS (from terra::cellSize on the
#'                    BC Albers projected raster template).
#' @param dates       Date vector of length n_months giving the first day of
#'                    each month represented by a column of ts_matrix.
#' @param index_type  "spi", "spei", "swei", or any other index key.
#' @param scale       Integer timescale (e.g. 1, 3, 6, 12).
#' @param out_dir     Directory to write the CSV.  Should be the index-specific
#'                    seasonal results directory (SPI_SEAS_DIR etc.) so that
#'                    load_basin_avg_csv() can find it without specifying data_dir.
#' @return Invisible data.frame (Year ?? month-name columns) as written to disk.
save_basin_avg_from_pixels <- function(ts_matrix, area_weights, dates,
                                       index_type, scale, out_dir,
                                       min_coverage = BASIN_MIN_VALID_AREA_COVERAGE,
                                       total_weight = NULL) {
  ## area_weights must be cell area x fractional basin coverage (see
  ## basin_fractional_weights()); total_weight is the whole-basin weight
  ## (defaults to the sum of the supplied weights, i.e. rows = all basin cells).
  ## ?????? Input validation ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  if (!is.matrix(ts_matrix))
    stop("ts_matrix must be a numeric matrix (pixels ?? months)")
  if (nrow(ts_matrix) != length(area_weights))
    stop(sprintf(
      "ts_matrix has %d rows but area_weights has %d elements",
      nrow(ts_matrix), length(area_weights)))
  if (ncol(ts_matrix) != length(dates))
    stop(sprintf(
      "ts_matrix has %d columns but dates has %d elements",
      ncol(ts_matrix), length(dates)))
  
  ## Guard: fractional weights are >= 0 by construction. A zero weight means the
  ## cell is outside the basin and is excluded; it is never replaced by a
  ## median area (that would silently add outside-basin cells to the mean).
  bad_w <- !is.finite(area_weights) | area_weights < 0
  if (any(bad_w))
    stop(sprintf("save_basin_avg_from_pixels: %d non-finite or negative area weight(s).",
                 sum(bad_w)))
  
  ## ?????? Fractional-coverage area-weighted mean with the shared coverage rule ??????
  basin_avg <- regional_weighted_mean(ts_matrix, area_weights,
                                      min_coverage = min_coverage,
                                      total_weight = total_weight)
  
  ## ?????? Reshape to Year ?? Month wide format ????????????????????????????????????????????????????????????????????????????????????
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
  
  ## ?????? Write CSV ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(
    out_dir,
    sprintf("%s_%02d_basin_averaged_by_month.csv",
            tolower(index_type), as.integer(scale)))
  utils::write.csv(df, out_file, row.names = FALSE)
  cat(sprintf(
    "  ??? Basin-avg CSV written (area-weighted, %d pixels, %d months): %s\n",
    nrow(ts_matrix), ncol(ts_matrix), basename(out_file)))
  invisible(df)
}

## B3. PACKAGE BOOTSTRAP ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
utils_load_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org/")
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
  cat("??? Packages loaded:", paste(pkgs, collapse = ", "), "\n")
}

####################################################################################
# C. NETCDF FILE FINDERS ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
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

## D. NETCDF DATE PARSING ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
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

## E. BASIN SPATIAL HELPERS ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
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
  cat(sprintf("  Basin area: %.1f km??  |  Valid cells: %d\n",
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
## E2. Area-weighted pixel-matrix column means ??????????????????????????????????????????????????????????????????????????????????????????
#' Compute area-weighted column means of a pixel ?? time matrix.
#'
#' @param mat     Numeric matrix [n_all_cells ?? n_time].
#'                NA for non-basin or missing cells.
#' @param weights Numeric vector [n_all_cells].
#'                Set to 0 (not NA) for non-basin cells.
#' @return Numeric vector of length n_time.
area_weighted_colmeans <- function(mat, weights,
                                   min_coverage = BASIN_MIN_VALID_AREA_COVERAGE) {
  regional_weighted_mean(mat, weights, min_coverage = min_coverage)
}
## F. DROUGHT EVENT DETECTION — SINGLE FRAMEWORK-WIDE DEFINITION ##############
## Every script (SSI, TSDI/NTSDI, trend tests, basin time series, event
## ranking) identifies drought events with the functions in this section.
##
##   Method       : threshold-level run method with hysteresis
##   Onset        : first month with x_t <  DROUGHT_ONSET
##   Termination  : first month with x_t >= DROUGHT_END (event ends the month before)
##   Missing month: terminates an active event (gaps are never bridged)
##   Min duration : scale-dependent, EVENT_MIN_DURATION_BY_SCALE (below);
##                  shorter runs are not events
##   Intensity    : I = mean(DROUGHT_END - x_t) over the event months (>= 0)
##   Severity     : S = I x D = sum(DROUGHT_END - x_t)
##
## DROUGHT_ONSET and DROUGHT_END are set in section A (both -0.5 for now, so
## the hysteresis band is closed and the method equals the single-threshold
## run method). Raising DROUGHT_END (e.g. to 0) opens the band; nothing else
## needs to change.
##
## USER SETTING: minimum event duration (months) per accumulation scale.
## Edit the numbers here; every script reads them through
## get_event_min_duration(). Non-accumulated monthly series (station SSI,
## SSPI-1, TSDI/NTSDI) use the scale-1 entry.
EVENT_MIN_DURATION_BY_SCALE <- c(
  `1`  =  2L,
  `3`  =  3L,
  `6`  =  4L,
  `12` =  6L,
  `24` = 12L,
  `36` = 12L
)

.validate_event_settings <- function() {
  t <- EVENT_MIN_DURATION_BY_SCALE
  if (is.null(names(t)) || anyNA(suppressWarnings(as.integer(names(t)))) ||
      anyDuplicated(names(t)) || any(!is.finite(t)) || any(t < 1) || any(t != round(t)))
    stop("EVENT_MIN_DURATION_BY_SCALE must be named by scale with positive whole-month values.")
  if (!is.numeric(DROUGHT_ONSET) || !is.numeric(DROUGHT_END) ||
      length(DROUGHT_ONSET) != 1L || length(DROUGHT_END) != 1L ||
      !is.finite(DROUGHT_ONSET) || !is.finite(DROUGHT_END))
    stop("DROUGHT_ONSET and DROUGHT_END must be single finite numbers.")
  if (DROUGHT_END < DROUGHT_ONSET)
    stop("DROUGHT_END (termination) must be >= DROUGHT_ONSET (onset).")
  invisible(TRUE)
}
.validate_event_settings()

#' Minimum event duration (months) for an accumulation scale.
get_event_min_duration <- function(scale, table = EVENT_MIN_DURATION_BY_SCALE) {
  if (length(scale) != 1L || !is.finite(suppressWarnings(as.numeric(scale))))
    stop("get_event_min_duration(): scale must be a single number.")
  key <- as.character(as.integer(scale))
  v <- table[key]
  if (is.na(v))
    stop(sprintf(paste0("No minimum event duration is defined for scale %s. ",
                        "Add it to EVENT_MIN_DURATION_BY_SCALE in DROUGHT_ANALYSIS_utils.R."), key))
  as.integer(unname(v))
}

.resolve_event_min_duration <- function(scale = NULL, min_duration = NULL) {
  if (!is.null(min_duration)) {
    md <- as.integer(min_duration)
    if (length(md) != 1L || is.na(md) || md < 1L) stop("min_duration must be a positive integer.")
    return(md)
  }
  if (!is.null(scale)) return(get_event_min_duration(scale))
  stop("Supply the index `scale` (preferred) or an explicit `min_duration`.")
}

#' Core run detector on a regular (gap-free calendar) monthly vector.
#' Returns one row per event: start_index, end_index, duration_months,
#' right_censored (still active at the end of the vector) and
#' gap_terminated (closed by a missing month).
drought_event_runs <- function(x, min_duration,
                               onset_threshold       = DROUGHT_ONSET,
                               termination_threshold = DROUGHT_END) {
  x <- as.numeric(x); n <- length(x)
  min_duration <- .resolve_event_min_duration(min_duration = min_duration)
  if (termination_threshold < onset_threshold)
    stop("termination_threshold must be >= onset_threshold.")
  st <- en <- integer(n); cz <- gp <- logical(n); k <- 0L; in_d <- FALSE; s <- 0L
  for (j in seq_len(n)) {
    v <- x[j]
    if (!is.finite(v)) {
      if (in_d && (j - s) >= min_duration) {
        k <- k + 1L; st[k] <- s; en[k] <- j - 1L; gp[k] <- TRUE
      }
      in_d <- FALSE
    } else if (!in_d) {
      if (v < onset_threshold) { in_d <- TRUE; s <- j }
    } else if (v >= termination_threshold) {
      if ((j - s) >= min_duration) { k <- k + 1L; st[k] <- s; en[k] <- j - 1L }
      in_d <- FALSE
    }
  }
  if (in_d && (n - s + 1L) >= min_duration) {
    k <- k + 1L; st[k] <- s; en[k] <- n; cz[k] <- TRUE
  }
  i <- seq_len(k)
  data.frame(start_index = st[i], end_index = en[i],
             duration_months = en[i] - st[i] + 1L,
             right_censored = cz[i], gap_terminated = gp[i])
}

#' Events with intensity / severity metrics (and dates when supplied).
drought_event_table <- function(x, dates = NULL, min_duration,
                                onset_threshold       = DROUGHT_ONSET,
                                termination_threshold = DROUGHT_END) {
  x <- as.numeric(x)
  runs <- drought_event_runs(x, min_duration, onset_threshold, termination_threshold)
  k <- nrow(runs)
  m <- if (k) t(vapply(seq_len(k), function(i) {
    seg <- x[runs$start_index[i]:runs$end_index[i]]
    d <- termination_threshold - seg
    c(min(seg), mean(d), max(d), sum(d))
  }, numeric(4))) else matrix(numeric(0), 0L, 4L)
  runs$min_value <- m[, 1]
  runs$mean_int  <- m[, 2]
  runs$max_int   <- m[, 3]
  runs$severity  <- m[, 4]
  if (!is.null(dates)) {
    dates <- as.Date(dates)
    if (length(dates) != length(x)) stop("drought_event_table(): dates/x length mismatch.")
    runs$start_date <- dates[runs$start_index]
    runs$end_date   <- dates[runs$end_index]
  }
  attr(runs, "min_duration") <- as.integer(min_duration)
  runs
}

#' Place a (possibly gappy) monthly series on a complete monthly calendar.
#' Missing months become NA and therefore terminate events.
regularize_monthly_series <- function(dates, values) {
  m <- as.Date(format(as.Date(dates), "%Y-%m-01"))
  ok <- !is.na(m)
  if (!any(ok)) stop("regularize_monthly_series(): no valid dates.")
  if (anyDuplicated(m[ok])) stop("regularize_monthly_series(): duplicated calendar months.")
  grid <- seq(min(m[ok]), max(m[ok]), by = "month")
  v <- rep(NA_real_, length(grid))
  v[match(m[ok], grid)] <- as.numeric(values)[ok]
  list(dates = grid, values = v, input_pos = match(m, grid))
}

#' Events for a dated monthly series of any spacing (e.g. seasonal records).
#' attr(result, "in_event") flags, for every INPUT row, whether that month
#' belongs to a detected event (used for figure shading).
identify_drought_events_ts <- function(dates, values, scale = NULL, min_duration = NULL,
                                       onset_threshold       = DROUGHT_ONSET,
                                       termination_threshold = DROUGHT_END) {
  md  <- .resolve_event_min_duration(scale, min_duration)
  reg <- regularize_monthly_series(dates, values)
  ev  <- drought_event_table(reg$values, reg$dates, md, onset_threshold, termination_threshold)
  in_ev <- logical(length(reg$values))
  for (i in seq_len(nrow(ev))) in_ev[ev$start_index[i]:ev$end_index[i]] <- TRUE
  pos <- reg$input_pos
  attr(ev, "in_event") <- !is.na(pos) & in_ev[ifelse(is.na(pos), 1L, pos)]
  attr(ev, "min_duration") <- md
  ev
}

#' Data-frame interface (date, value) used by the SSI, basin time-series and
#' Excel-summary code. Requires a complete monthly calendar.
detect_drought_events <- function(df,
                                  onset_threshold       = DROUGHT_ONSET,
                                  termination_threshold = DROUGHT_END,
                                  min_duration          = NULL,
                                  scale                 = NULL,
                                  severe_thr            = -1.5,
                                  extreme_thr           = -2.0) {
  stopifnot(all(c("date", "value") %in% names(df)))
  md <- .resolve_event_min_duration(scale, min_duration)
  df <- df[order(df$date), ]
  validate_monthly_record(df$date, label = "detect_drought_events input",
                          require_complete = TRUE)
  ev <- drought_event_table(df$value, df$date, md, onset_threshold, termination_threshold)
  classify_severity <- function(min_v)
    ifelse(min_v <= extreme_thr, "D3: extremely dry",
           ifelse(min_v <= severe_thr, "D2: severely dry", "D1: moderately dry"))
  out <- data.frame(
    event_id        = seq_len(nrow(ev)),
    start_date      = ev$start_date,
    end_date        = ev$end_date,
    duration_months = ev$duration_months,
    min_value       = ev$min_value,
    severity_class  = if (nrow(ev)) classify_severity(ev$min_value) else character(),
    mean_int        = ev$mean_int,
    max_int         = ev$max_int,
    severity        = ev$severity,
    right_censored  = ev$right_censored,
    gap_terminated  = ev$gap_terminated,
    stringsAsFactors = FALSE)
  attr(out, "min_duration") <- md
  out
}

## F2. SPATIAL HELPERS (used by 4_trends_visualization.R) ??????????????????????????????
#' Load the authoritative Nechako basin from the supplied subbasin
#' Shapefiles and transform the dissolved union to the requested CRS.
load_basin <- function(shp_path = BASIN_SHP, crs = EQUAL_AREA_CRS) {
  if (!dir.exists(shp_path))
    stop("Nechako subbasin Shapefile directory not found: ", shp_path)
  b <- terra::project(
    nechako_subbasin_union(load_nechako_subbasins(shp_path)),
    crs
  )
  b <- sf::st_as_sf(b)
  cat(sprintf("Nechako basin loaded from %d Shapefiles (dissolved, CRS: %s)\n",
              NECHAKO_SUBBASIN_EXPECTED_N, crs))
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

#' 3??3 focal mean smoother (passes NA through gracefully)
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

## G. GGPLOT HELPERS ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
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

index_colours <- c(spi = "#e41a1c", spei = "#377eb8", spei_era5 = "#d62728", swei = "#1E90FF")

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

#' Unified figure-saving helper ??? replaces the scattered mix of pdf()/dev.off()
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
    error = function(e) cat(sprintf("  ??? PDF save failed (%s): %s\n",
                                    basename(pdf_path), e$message))
  )
  tryCatch(
    ggplot2::ggsave(png_path, plot_obj,
                    width = width, height = height, units = "in",
                    dpi = dpi, device = "png"),
    error = function(e) cat(sprintf("  ??? PNG save failed (%s): %s\n",
                                    basename(png_path), e$message))
  )
  cat(sprintf("  ??? Figure saved: %s (.pdf + .png)\n", basename(stem)))
  invisible(NULL)
}

## H. EXCEL SUMMARY EXPORT ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
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
        # Events detected on the full record, counted if they overlap 2020-2025.
        ev_all    <- identify_drought_events_ts(df$date, df$value, scale = sc)
        ev_count  <- sum(ev_all$end_date >= as.Date("2020-01-01") &
                           ev_all$start_date <= as.Date("2025-12-31"))
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
        cat(sprintf("  ??? Stats for %s-%02d: %s\n", toupper(idx), sc, e$message)))
    }
  }
  openxlsx::addWorksheet(wb, "Summary_Statistics")
  openxlsx::writeData(wb, "Summary_Statistics", stats_df)
  openxlsx::addStyle(wb, "Summary_Statistics", hdr_style,
                     rows = 1, cols = seq_len(ncol(stats_df)))
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  cat(sprintf("??? Excel summary saved: %s\n", basename(output_file)))
  invisible(wb)
}

## I. W9 SHARED-STATE LOADER ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
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
  
  cat(sprintf("??? w9 shared state loaded (%d keys): %s\n",
              length(st), basename(path)))
  st
}

####################################################################################
# J. SHARED STATISTICAL FUNCTIONS
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# These were previously duplicated between 6trend_test_ALL_v2.R (vectorised pixel
# versions) and 4pr_pet_trends_v2.r (single-series version).  They are now defined
# once here and sourced by both scripts.
#
# The public API is intentionally backward-compatible:
#   ??? mk_tfpw_spectral_for_series() ??? identical signature to the function that
#     was previously embedded in 4pr_pet_trends_v2.r.  That script calls it
#     unchanged.
#   ??? vectorized_mann_kendall(), vectorized_mann_kendall_tfpw(),
#     vectorized_regime_shift_pelt(), vectorized_spectral_peaks() ??? identical
#     signatures to the functions previously in 6trend_test_ALL_v2.R.  That script
#     now omits its local copies and relies on these definitions.
#   ??? detect_regime_shift_pelt() ??? single-series PELT helper, previously in
#     4pr_pet_trends_v2.r under that name.  "trend_test_ALL.R" used a different but equivalent
#     vectorised wrapper; both now share this single-series core.
####################################################################################

# ?????? J0. Utility operator ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.atomic(a)) {
    if (length(a) > 0 && !is.na(a[1])) return(a) else return(b)
  }
  # Non-atomic objects (SpatRaster, PackedSpatRaster, list, S4, data.frame, etc.)
  # can't be NA-checked via a[1] - existence (non-NULL) is enough.
  a
}

# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# J1. PRIMITIVE HELPERS
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

#' Check whether too many values in a series are at or below a floor threshold.
#' Returns a list(exceeds_threshold, pct_min_vals).
#' Used by mk_tfpw_spectral_for_series() (formerly embedded in 4pr_pet_trends_v2.r).
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
#' Used by mk_tfpw_for_series() (formerly embedded in 4pr_pet_trends_v2.r).
calculate_variance_with_ties <- function(S, n, x) {
  tie_table  <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  var_s      <- n * (n - 1) * (2 * n + 5) / 18
  if (length(tie_counts))
    var_s <- var_s - sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
  var_s
}

# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# J2. PELT CHANGEPOINT ??? single series
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#' Detect the first mean-shift regime in a univariate time series using PELT/BIC.
#' Returns a named list with fields:
#'   changepoint_detected, changepoint_position, n_changepoints,
#'   first_changepoint_year, mean_before, mean_after, magnitude_shift.
#'
#' Previously defined as detect_regime_shift_pelt() in 4pr_pet_trends_v2.r;
#' 6trend_test_ALL_v2.R used an inline equivalent inside vectorized_regime_shift_pelt().
#' Both now share this function.
detect_regime_shift_pelt <- function(ts_vec, min_obs = 20, start_year = NULL,
                                     minseglen = 5L, standardize = TRUE) {
  # changepoint::cpt.mean() with the default Normal test statistic ASSUMES UNIT
  # VARIANCE. On series whose sd is >> 1 (annual precipitation ~80 mm, PET_ERA5
  # ~300 mm) the BIC penalty is then negligible relative to the cost, and
  # spurious breaks are placed at the series edge (first_changepoint_year 1950/51).
  # standardize = TRUE z-scores the series before fitting; reported means and the
  # shift magnitude are still computed on the ORIGINAL values.
  # minseglen forbids segments shorter than `minseglen` observations, so a break
  # cannot sit within `minseglen` years of either end.
  ts_clean <- as.numeric(stats::na.omit(ts_vec))
  n <- length(ts_clean)
  empty <- list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                n_changepoints = 0L, first_changepoint_year = NA_integer_,
                mean_before = NA_real_, mean_after = NA_real_, magnitude_shift = NA_real_)
  if (n < min_obs || n < 2L * minseglen) return(empty)
  tryCatch({
    x_fit <- ts_clean
    if (isTRUE(standardize)) {
      s <- stats::sd(ts_clean)
      if (!is.finite(s) || s <= 0) return(empty)
      x_fit <- (ts_clean - mean(ts_clean)) / s
    }
    cpt      <- changepoint::cpt.mean(x_fit, method = "PELT", penalty = "BIC",
                                      minseglen = minseglen)
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

# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# J3. VC-MK ??? single series (core used by both the vectorised wrapper and
#     mk_tfpw_spectral_for_series)
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#' Parse modifiedmk::mmkh() output by NAME PATTERN (punctuation differs between
#' fields: "new P-value", "old P.value", "new.variance", "N/N*"). Returns NA for
#' any field that is missing so callers can decide how to degrade.
.mmkh_fields <- function(z) {
  nm <- names(z)
  getv <- function(pattern) {
    k <- grep(pattern, nm, ignore.case = TRUE)
    if (length(k)) suppressWarnings(as.numeric(z[[k[1]]])) else NA_real_
  }
  out <- list(tau      = getv("^Tau$"),
              p_new    = getv("^new[ .]?P"),
              p_old    = getv("^old[ .]?P"),
              var_new  = getv("^new[ .]?var"),
              var_old  = getv("^old[ .]?var"),
              n_over_nstar = getv("^N/N"))
  # Fall back to the variance ratio if the N/N* field is absent.
  if (!is.finite(out$n_over_nstar) && is.finite(out$var_new) &&
      is.finite(out$var_old) && out$var_old > 0)
    out$n_over_nstar <- out$var_new / out$var_old
  out
}

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
  }, error = function(e) {
    
    message("VC error: ", conditionMessage(e))
    
    list(
      tau = NA_real_,
      pval = NA_real_,
      slope = NA_real_,
      filtered = TRUE
    )
  })
}

# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# J4. TFPW-MK ??? single series
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
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
    
    # Step 3: Detrend ??? lag-1 autocorrelation
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
         slope        = beta,    # ??? original-series Sen's slope (Yue et al. 2002)
         rho1         = r1,
         tfpw_applied = abs(r1) > 0.1,
         filtered     = FALSE)
  }, error = function(e)
    list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
         rho1 = NA_real_, tfpw_applied = FALSE, filtered = TRUE))
}

# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# J5. SPECTRAL PEAK DETECTION ??? single series
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#' AR(1) red-noise spectral test for a single detrended time series.
#' Counts contiguous groups of periodogram bins above the chi-squared threshold
#' at periods >= min_period_months (default 24 months).
#' Returns list(n_peaks, dominant_period, conf).
#' Internal helper; not normally called directly by analysis scripts.
#'
#' [BUGFIX, see inline comments in Method 2 below] Both the AR(1) branch
#' (Method 1) and the Monte-Carlo white-noise branch (Method 2, used whenever
#' conf_cache_env is supplied -- i.e. by 4pr_pet_trends_v2.r) now restrict
#' testing to periods >= min_period_months. Previously only Method 1 applied
#' this restriction; Method 2 tested the full frequency range, which inflated
#' "significant peak" detections to near-universal rates unrelated to genuine
#' multi-year drought periodicity.
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
    # ?????? Method 1: AR(1) red-noise via smoothed periodogram (used by 6trend_test_ALL) ??????????????????
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
    
    # ?????? Method 2: Monte-Carlo white-noise envelope (used by 4pr_pet_trends_v2.r) ???
    # Also compute this for the omnibus wrapper so it can serve both scripts.
    #
    # [BUGFIX] This branch previously computed max_spectra (the null) and then
    # tested sp_mc (the observed spectrum) over the FULL frequency range -- it
    # never applied the min_period_months >= 24 drought-band restriction that
    # Method 1 (AR1) applies above. That means "significant peak" was really
    # testing for ANY significant frequency, including sub-annual/seasonal-
    # residual power, not multi-year drought periodicity. This is almost
    # certainly why n_spectral_peaks was flagged "significant" at ~100% of
    # pixels across every index/scale (see 6trend_test_ALL_v2.R QA): the test
    # was not restricted to the periods it was documented to test.
    #
    # Fix: restrict BOTH the null (max_spectra) and the observed periodogram
    # (sp_mc) to the same drought-band frequency bins (period >= min_period_months)
    # before comparing. Using the MAX over the restricted band as the null
    # statistic (rather than testing each bin independently) keeps this a
    # single family-wise threshold across all bins in the band, so the true
    # per-pixel false-positive rate stays at ~alpha regardless of how many
    # bins happen to fall inside the band.
    conf_limit <- NA_real_
    if (!is.null(conf_cache_env)) {
      half      <- floor(n / 2)
      freqs_mc  <- seq(0, 0.5, length.out = half)
      period_mc <- ifelse(freqs_mc > 0, 1.0 / freqs_mc, Inf)
      band_mc   <- is.finite(period_mc) & period_mc >= min_period_months
      
      if (!any(band_mc)) {
        return(list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_))
      }
      
      # Cache key now includes the band definition so thresholds computed
      # under different min_period_months settings never collide.
      key <- paste0("n_", n, "_band", min_period_months)
      if (!exists(key, envir = conf_cache_env, inherits = FALSE)) {
        max_spectra <- vapply(seq_len(n_sim), function(jj) {
          r  <- rnorm(n, 0, 1)
          fr <- fft(r - mean(r))
          sp <- Mod(fr[seq_len(half)])^2 / n
          max(sp[band_mc], na.rm = TRUE)   # global max WITHIN the drought band only
        }, numeric(1))
        assign(key, stats::quantile(max_spectra, 1 - alpha, na.rm = TRUE),
               envir = conf_cache_env)
      }
      conf_limit <- get(key, envir = conf_cache_env, inherits = FALSE)
      
      # Re-check with MC envelope (used by 4pr_pet_trends_v2.r path)
      sd_d    <- stats::sd(x_dt); if (!is.finite(sd_d) || sd_d == 0) sd_d <- 1
      dstd    <- (x_dt - mean(x_dt)) / sd_d
      ff      <- fft(dstd)
      sp_mc   <- Mod(ff[seq_len(half)])^2 / n
      # Only consider bins within the drought band -- this is the actual fix.
      pk_mc   <- which(band_mc & (sp_mc > conf_limit))
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

# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# J6. OMNIBUS SINGLE-SERIES WRAPPER
#     mk_tfpw_spectral_for_series()
#     Used directly by 4pr_pet_trends_v2.r (and its visualization script reads
#     outputs produced by it).  Signature is identical to the version that
#     was previously embedded in that script.
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
#' Full VC-MK + TFPW-MK + spectral + PELT analysis for a single time series.
#'
#' @param ts_vec           Numeric vector (may contain NAs).
#' @param is_precip        Logical; TRUE ??? precipitation quality filters.
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
#' @param also_check_min_pct   If TRUE and min_nonzero_count is given, ALSO reject
#'   series whose floor-value share exceeds max_min_value_pct_* (default FALSE keeps
#'   the historical behaviour, where min_nonzero_count replaced the % check).
#' @param cap_vc_variance_shrinkage  If TRUE (default), the Hamed-Rao correction may
#'   only INFLATE the MK variance (n/n* >= 1). When it would shrink it (negative
#'   autocorrelation) or produce a negative variance (NaN p), the uncorrected MK
#'   p-value is used and the row is tagged in `reason`.
#' @param pelt_minseglen       Minimum PELT segment length (default 5 obs).
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
                                        min_obs_changepoint  = 20L,
                                        also_check_min_pct   = FALSE,
                                        cap_vc_variance_shrinkage = TRUE,
                                        pelt_minseglen       = 5L) {
  ts_clean <- na.omit(ts_vec)
  n        <- length(ts_clean)
  
  na_result <- function(reason, n_ties = 0, pct_ties = 0, pct_min = 0) {
    list(
      vc  = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_,
                 varS = NA_real_, n = n, rho1 = NA_real_, vc_corrected = FALSE,
                 n_ties = n_ties, percent_ties = pct_ties, n_min_vals = 0L,
                 percent_min_vals = pct_min, tau_b_adjusted = FALSE,
                 filtered = TRUE, reason = reason,  mmkh_status = "not_run"),
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
  
  # ?????? Quality filters ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  if (skip_min_filter) {
    min_check <- list(exceeds_threshold = FALSE, pct_min_vals = 0)
  } else if (!is.null(min_nonzero_count)) {
    n_nonzero    <- sum(ts_clean > min_positive_value, na.rm = TRUE)
    pct_min_vals <- (1 - n_nonzero / n) * 100
    if (n_nonzero < min_nonzero_count)
      return(na_result("insufficient_nonzero", pct_min = pct_min_vals))
    if (isTRUE(also_check_min_pct) &&
        pct_min_vals > (if (is_precip) max_min_value_pct_precip else max_min_value_pct_pet))
      return(na_result("excessive_min_vals", pct_min = pct_min_vals))
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
  
  # ?????? Sen's slope (on original series, shared by both MK variants) ????????????????????????????????????
  sen_slope <- trend::sens.slope(ts_clean)$estimates
  if (is.na(sen_slope) || is.infinite(sen_slope) || length(sen_slope) == 0)
    return(na_result("sens_slope_fail", n_ties, percent_ties,
                     min_check$pct_min_vals))
  
  # ?????? Detrend ??? lag-1 autocorrelation ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  time_index <- seq_len(n)
  trend_line <- sen_slope * time_index
  detrended  <- ts_clean - trend_line
  acf_result <- tryCatch(
    acf(detrended, lag.max = 1, plot = FALSE, na.action = na.pass),
    error = function(e) NULL, warning = function(w) NULL)
  rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA_real_
  
  # ?????? VC-MK (Hamed-Rao 1998) ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  # Use mmkh() ??? the same function used by .mk_vc_single() (J3) and
  # vectorized_mann_kendall() (J7a) ??? so all three code paths produce
  # consistent H-R98 variance-corrected p-values.
  # mk3() is a different modified-MK variant and must NOT be used here.
  vc_corrected <- FALSE
  p_vc <- NA_real_; S_vc <- NA_real_; tau_vc <- NA_real_; varS_final <- NA_real_
  vc_failed <- FALSE; vc_reason <- "none"; vc_capped <- FALSE
  n_over_nstar <- NA_real_; p_vc_uncorrected <- NA_real_
  vc_res <- tryCatch(suppressWarnings(modifiedmk::mmkh(ts_clean)), error = function(e) NULL)
  if (!is.null(vc_res)) {
    mm <- .mmkh_fields(vc_res)          # pattern-matched: immune to name punctuation
    S_vc         <- NA_real_            # mmkh() doesn't return S directly
    tau_vc       <- mm$tau
    varS_final   <- mm$var_new
    n_over_nstar <- mm$n_over_nstar
    p_corr       <- mm$p_new
    p_unc        <- mm$p_old
    if (!is.finite(p_unc))
      p_unc <- tryCatch(as.numeric(Kendall::MannKendall(ts_clean)$sl), error = function(e) NA_real_)
    if (!is.finite(tau_vc))
      tau_vc <- tryCatch(as.numeric(Kendall::MannKendall(ts_clean)$tau), error = function(e) NA_real_)
    p_vc_uncorrected <- p_unc
    
    var_ok  <- is.finite(varS_final) && varS_final > 0 && is.finite(p_corr)
    shrunk  <- is.finite(n_over_nstar) && n_over_nstar < 1
    ## Hamed-Rao n/n* < 1 (negative autocorrelation) SHRINKS Var(S) and yields
    ## anti-conservative p-values (e.g. tau = 0.06 with p = 0.0007). n/n* <= 0 gives a
    ## negative variance and sqrt(VS) = NaN. In both cases fall back to the
    ## uncorrected MK p-value instead of discarding or over-trusting the series.
    ## A VC failure must NOT early-return: TFPW / spectral / PELT do not depend on it.
    if (var_ok && !(isTRUE(cap_vc_variance_shrinkage) && shrunk)) {
      p_vc         <- p_corr
      vc_corrected <- is.finite(n_over_nstar) && n_over_nstar > 1
    } else if (isTRUE(cap_vc_variance_shrinkage) && is.finite(p_unc)) {
      p_vc       <- p_unc
      varS_final <- mm$var_old
      vc_capped  <- TRUE
      vc_reason  <- if (var_ok) "vc_shrinkage_capped_uncorrected" else "vc_invalid_variance_uncorrected"
    } else {
      vc_failed <- TRUE
      vc_reason <- if (!is.finite(varS_final) || varS_final <= 0) "vc_invalid_variance" else "vc_nonfinite_p"
    }
    if (!vc_failed && !is.finite(p_vc)) { vc_failed <- TRUE; vc_reason <- "vc_nonfinite_p" }
    vc_status <- if (vc_failed) "invalid" else if (vc_capped) "capped" else "ok"
  } else {
    vc_failed <- TRUE; vc_reason <- "mmkh_error"
    vc_status <- "mmkh_error"
  }
  n_min_obs <- as.integer(round(min_check$pct_min_vals / 100 * n))
  vc_list <- list(tau = if (vc_failed) NA_real_ else tau_vc,
                  p   = if (vc_failed) NA_real_ else p_vc,
                  sl = sen_slope, S = S_vc,
                  varS = varS_final, n = n, rho1 = rho1,
                  vc_corrected = vc_corrected,
                  n_ties = n_ties, percent_ties = percent_ties,
                  n_min_vals = n_min_obs, percent_min_vals = min_check$pct_min_vals,
                  tau_b_adjusted = tau_b_adjusted,
                  filtered = vc_failed, reason = vc_reason, mmkh_status = vc_status,
                  vc_n_over_nstar = n_over_nstar, vc_p_uncorrected = p_vc_uncorrected,
                  vc_capped = vc_capped)
  # execution falls through unconditionally into TFPW-MK / spectral / PELT below
  
  # ?????? TFPW-MK ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  tfpw_applied <- FALSE
  if (!is.na(rho1) && abs(rho1) > 0.1) {
    pw    <- numeric(n); pw[1] <- detrended[1]
    for (j in 2:n) pw[j] <- detrended[j] - rho1 * detrended[j - 1]
    corrected    <- pw + trend_line
    tfpw_applied <- TRUE
  } else {
    corrected <- ts_clean
  }
  # outer(X, Y, FUN)[i, j] = FUN(X[i], Y[j]); the upper triangle is i < j, so
  # FUN(a, b) = b - a gives (later - earlier), the standard MK S. The previous
  # `-` gave (earlier - later), i.e. S and tau with the WRONG SIGN (p unaffected).
  s_mat_tf <- sign(outer(corrected, corrected, function(a, b) b - a))
  S_tf     <- sum(s_mat_tf[upper.tri(s_mat_tf)], na.rm = TRUE)
  varS_tf  <- calculate_variance_with_ties(S_tf, n, corrected)
  n_pairs  <- n * (n - 1) / 2
  tau_tf   <- S_tf / n_pairs
  p_tf     <- if (varS_tf <= 0) NA_real_ else 2 * pnorm(-abs(S_tf / sqrt(varS_tf)))
  tf_list  <- list(tau = tau_tf, p = p_tf, sl = sen_slope, S = S_tf,
                   varS = varS_tf, n = n, rho1 = rho1,
                   tfpw_applied = tfpw_applied,
                   n_ties = n_ties, percent_ties = percent_ties,
                   n_min_vals = n_min_obs, percent_min_vals = min_check$pct_min_vals,
                   tau_b_adjusted = tau_b_adjusted, filtered = FALSE, reason = "none")
  
  # ?????? Spectral ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  spec_result <- .spectral_single(ts_clean,
                                  min_period_months = 24L,
                                  sig_level         = 0.95,
                                  n_sim             = n_sim_spectral,
                                  alpha             = alpha,
                                  conf_cache_env    = conf_cache_env)
  
  # ?????? PELT ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  cpt_result <- detect_regime_shift_pelt(ts_clean,
                                         min_obs    = min_obs_changepoint,
                                         start_year = start_year,
                                         minseglen  = pelt_minseglen)
  
  list(vc = vc_list, tf = tf_list, spec = spec_result, cpt = cpt_result)
}

# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# J7. VECTORISED PIXEL-MATRIX WRAPPERS
#     Used directly by 6trend_test_ALL_v2.R.  Signatures are identical to the
#     functions previously defined in that script.
#     The parallel back-end variables N_CORES, is_windows, cl are expected to
#     be set in the calling script (6trend_test_ALL_v2.R) before these are called.
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????

## J7a. Hamed-Rao (1998) Variance-Corrected MK + Sen's slope ??????????????????????????????????????????????????????
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

## J7b. TFPW Mann-Kendall (Yue et al. 2002) ?????????????????????????????????????????????????????????????????????????????????????????????????????????
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

## J7c. PELT regime-shift ??? pixel matrix ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????
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

## J7d. Spectral peak detection ??? pixel matrix ???????????????????????????????????????????????????????????????????????????????????????????????????
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
cat("??? DROUGHT_ANALYSIS_utils_v2.R loaded\n")
################################################################################
# K. NECHAKO REGIONAL SPATIAL FRAMEWORK
#
# Authoritative spatial source: the 7 supplied ESRI Shapefiles under
#   Spatial/Nechako_shps/<SubbasinName>.shp
# Source CRS: WGS 84 geographic coordinates (EPSG:4326; CRS84-equivalent files
# are accepted). The whole Nechako basin is the dissolved union of these 7
# supplied subbasin polygons. No KMZ/KML boundary is used.
################################################################################

NECHAKO_SUBBASIN_DIR <- file.path(WD_PATH, "Spatial", "Nechako_shps")
NECHAKO_SUBBASIN_FILES <- c(
  "Chilako.shp",
  "Endako.shp",
  "Lower_Nechako.shp",
  "Lower_Stuart.shp",
  "Stellako.shp",
  "Upper_Nechako.shp",
  "Upper_Stuart.shp"
)
NECHAKO_SUBBASIN_EXPECTED_N <- length(NECHAKO_SUBBASIN_FILES)
NECHAKO_SUBBASIN_SOURCE_CRS <- "EPSG:4326"
NECHAKO_REGION_BASIN_ID <- "Nechako_Basin"
NECHAKO_REGION_OUTPUT_DIR <- file.path(TREND_DIR, "regional_timeseries")

sanitize_region_id <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) stop("Cannot construct a non-empty region_id from: ", x)
  x
}

load_nechako_subbasins <- function(subbasin_dir = NECHAKO_SUBBASIN_DIR,
                                   expected_n = NECHAKO_SUBBASIN_EXPECTED_N,
                                   native_crs = NECHAKO_SUBBASIN_SOURCE_CRS) {
  if (!dir.exists(subbasin_dir))
    stop("Nechako subbasin Shapefile directory not found: ", subbasin_dir)
  
  expected_files <- NECHAKO_SUBBASIN_FILES
  if (length(expected_files) != expected_n)
    stop("Internal configuration error: expected_n does not match the authoritative Shapefile list.")
  
  shp <- file.path(subbasin_dir, expected_files)
  missing_files <- shp[!file.exists(shp)]
  if (length(missing_files)) {
    stop(
      "Missing authoritative Nechako subbasin Shapefile(s): ",
      paste(basename(missing_files), collapse = ", "),
      "\nExpected directory: ", subbasin_dir
    )
  }
  
  shp <- normalizePath(shp, winslash = "/", mustWork = TRUE)
  if (length(shp) != expected_n)
    stop("Expected exactly ", expected_n, " authoritative subbasin Shapefiles; found ", length(shp))
  
  regions <- lapply(shp, function(f) {
    v <- tryCatch(
      terra::vect(f),
      error = function(e)
        stop("Failed to read subbasin Shapefile ", f, ": ", conditionMessage(e))
    )
    if (nrow(v) < 1L)
      stop("Subbasin Shapefile contains no features: ", f)
    if (is.na(terra::crs(v)) || !nzchar(terra::crs(v)))
      stop("Subbasin has no CRS: ", f)
    if (nrow(v) > 1L)
      v <- terra::aggregate(v)
    terra::makeValid(v)
  })
  
  ids <- vapply(
    shp,
    function(f) sanitize_region_id(tools::file_path_sans_ext(basename(f))),
    character(1)
  )
  if (anyDuplicated(ids))
    stop("Duplicate subbasin region_id(s): ", paste(ids[duplicated(ids)], collapse = ", "))
  
  # The supplied .prj files are WGS 84 geographic CRS. Lower_Stuart uses a
  # CRS84-equivalent WKT; both are transformable to EPSG:4326 and therefore
  # represent the same geographic coordinate system for this workflow.
  for (i in seq_along(regions)) {
    cr <- terra::crs(regions[[i]])
    if (is.na(cr) || !nzchar(cr))
      stop("Subbasin has an invalid CRS: ", shp[i])
    tryCatch(
      terra::project(regions[[i]], native_crs),
      error = function(e)
        stop("Subbasin CRS cannot be transformed to ", native_crs,
             ": ", shp[i], " | ", conditionMessage(e))
    )
  }
  
  names(regions) <- ids
  attr(regions, "source_files") <- shp
  attr(regions, "native_crs") <- native_crs
  class(regions) <- c("nechako_subbasins", class(regions))
  regions
}

nechako_subbasin_union <- function(subbasins = NULL) {
  if (is.null(subbasins)) subbasins <- load_nechako_subbasins()
  if (!length(subbasins)) stop("No subbasins supplied for union.")
  
  # The 7 supplied subbasins are the authoritative partition of the Nechako
  # Basin. Build the whole-basin geometry by dissolving those polygons as a
  # single geometry.  Do not use attribute-based aggregation: attributes are
  # irrelevant to the basin union and can prevent a true dissolve.
  v <- do.call(rbind, unname(subbasins))
  v <- terra::makeValid(v)
  u <- terra::aggregate(v, dissolve = TRUE)
  u <- terra::makeValid(u)
  u
}

load_nechako_regions <- function(template = NULL, subbasins = NULL,
                                 include_whole_basin = TRUE,
                                 analysis_crs = NULL) {
  if (is.null(subbasins)) subbasins <- load_nechako_subbasins()
  if (length(subbasins) != NECHAKO_SUBBASIN_EXPECTED_N)
    stop("Expected exactly ", NECHAKO_SUBBASIN_EXPECTED_N,
         " Nechako subbasins; found ", length(subbasins), ".")
  if (is.null(analysis_crs) && !is.null(template)) analysis_crs <- terra::crs(template)
  
  # Construct the authoritative basin from the same native subbasin objects
  # before any reprojection. This guarantees that basin and subbasin definitions
  # share exactly the same source geometry.
  basin <- if (include_whole_basin) nechako_subbasin_union(subbasins) else NULL
  
  regions <- list()
  if (include_whole_basin)
    regions[[NECHAKO_REGION_BASIN_ID]] <- basin
  for (id in names(subbasins)) regions[[id]] <- subbasins[[id]]
  
  if (!is.null(analysis_crs)) {
    regions <- lapply(regions, function(v) {
      if (is.na(terra::crs(v)) || !nzchar(terra::crs(v)))
        stop("Region has no valid CRS and cannot be transformed.")
      terra::project(v, analysis_crs)
    })
  }
  
  regions
}

region_catalog <- function(regions = NULL) {
  if (is.null(regions)) regions <- load_nechako_regions()
  data.table::rbindlist(lapply(names(regions), function(id) {
    v <- regions[[id]]
    data.table::data.table(
      region_id = id,
      region_name = gsub("_", " ", id, fixed = TRUE),
      area_km2 = as.numeric(sum(terra::expanse(v, unit = "km"))),
      n_features = nrow(v)
    )
  }), fill = TRUE)
}

validate_nechako_regions <- function(regions = NULL, tolerance = 1e-5,
                                     require_partition = TRUE) {
  if (is.null(regions)) regions <- load_nechako_regions()
  
  ids <- names(regions)
  expected <- c(NECHAKO_REGION_BASIN_ID,
                names(load_nechako_subbasins()))
  if (!setequal(ids, expected))
    stop("Regional definition mismatch. Expected: ", paste(expected, collapse = ", "))
  
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0)
    stop("tolerance must be one positive finite numeric value.")
  
  # Geometry validation must be performed in an equal-area CRS, independently
  # of the raster CRS supplied by the calling analysis script. The original
  # utility assumed that the incoming regions were already in an appropriate
  # CRS; 2PET supplies regions transformed to the ERA5-Land raster CRS, so that
  # assumption was false and made the area/overlap QA fragile.
  subs <- regions[setdiff(ids, NECHAKO_REGION_BASIN_ID)]
  if (length(subs) != NECHAKO_SUBBASIN_EXPECTED_N)
    stop("Expected exactly ", NECHAKO_SUBBASIN_EXPECTED_N,
         " subbasins; found ", length(subs), ".")
  
  # Avoid constructing a synthetic raster merely to test CRS equality.  The
  # target CRS is explicit and terra::project() is safe when source and target
  # CRS are identical, so always project to the declared equal-area CRS.
  subs_aea <- lapply(subs, function(v) {
    if (is.na(terra::crs(v)) || !nzchar(terra::crs(v)))
      stop("Subbasin has no valid CRS for geometric validation.")
    terra::makeValid(terra::project(v, EQUAL_AREA_CRS))
  })
  basin_aea <- {
    v <- regions[[NECHAKO_REGION_BASIN_ID]]
    if (is.na(terra::crs(v)) || !nzchar(terra::crs(v)))
      stop("Whole-basin region has no valid CRS for geometric validation.")
    terra::makeValid(terra::project(v, EQUAL_AREA_CRS))
  }
  
  areas <- vapply(
    subs_aea,
    function(v) sum(terra::expanse(v, unit = "m"), na.rm = TRUE),
    numeric(1)
  )
  union_area <- sum(terra::expanse(basin_aea, unit = "m"), na.rm = TRUE)
  if (!is.finite(union_area) || union_area <= 0)
    stop("Whole-basin union has invalid area.")
  if (any(!is.finite(areas)) || any(areas <= 0))
    stop("One or more subbasins has invalid/non-positive area.")
  
  # Independently reconstruct the authoritative union from the 7 supplied
  # subbasins. This is the reference geometry for the partition test; it avoids
  # making the result depend on how the caller happened to construct the basin
  # object after reprojection.
  sub_union_aea <- terra::aggregate(
    do.call(rbind, unname(subs_aea)),
    dissolve = TRUE
  )
  sub_union_aea <- terra::makeValid(sub_union_aea)
  sub_union_area <- sum(
    terra::expanse(sub_union_aea, unit = "m"),
    na.rm = TRUE
  )
  
  if (!is.finite(sub_union_area) || sub_union_area <= 0)
    stop("Dissolved 7-subbasin union has invalid area.")
  
  # First test: the supplied whole-basin object must agree with the union of
  # the 7 authoritative subbasins. This is distinct from the partition test.
  basin_union_rel_error <- abs(sub_union_area - union_area) / union_area
  if (basin_union_rel_error > tolerance) {
    warning(
      sprintf(
        "Whole-basin geometry differs from the dissolved 7-subbasin union by %.6f%%.",
        100 * basin_union_rel_error
      ),
      call. = FALSE
    )
  }
  
  if (require_partition) {
    # For a mutually exclusive partition, the sum of component areas must equal
    # the area of their dissolved union. Because the union is constructed from
    # exactly these 7 components, this test directly detects overlaps or gaps
    # without relying on the caller's projected basin object.
    rel_partition_error <- abs(sum(areas) - sub_union_area) / sub_union_area
    if (rel_partition_error > tolerance) {
      warning(
        sprintf(
          "The 7 subbasins do not form a mutually exclusive partition: sum-of-areas differs from dissolved union by %.6f%%.",
          100 * rel_partition_error
        ),
        call. = FALSE
      )
    }
  }
  
  # Pairwise overlap test. terra::intersect() emits a benign warning when two
  # mutually exclusive polygons have no intersection. That warning is NOT an
  # overlap. Suppress only that implementation warning and calculate overlap
  # area from polygonal intersection geometries that actually exist.
  overlap_area <- 0
  if (length(subs_aea) > 1L) {
    for (i in seq_len(length(subs_aea) - 1L)) {
      for (j in seq.int(i + 1L, length(subs_aea))) {
        inter <- tryCatch(
          suppressWarnings(terra::intersect(subs_aea[[i]], subs_aea[[j]])),
          error = function(e) {
            stop(
              "Failed pairwise subbasin intersection for ",
              names(subs_aea)[i], " vs ", names(subs_aea)[j],
              ": ", conditionMessage(e)
            )
          }
        )
        
        if (is.null(inter) || nrow(inter) == 0L)
          next
        
        # Shared boundaries can produce line/point intersections. These have
        # zero area and therefore do not constitute polygon overlap.
        poly_idx <- which(terra::geomtype(inter) %in% c("polygons", "polygon"))
        if (!length(poly_idx)) next
        
        a <- sum(
          terra::expanse(inter[poly_idx], unit = "m"),
          na.rm = TRUE
        )
        if (is.finite(a) && a > 0) {
          overlap_area <- overlap_area + a
        }
      }
    }
  }
  
  if (overlap_area > tolerance * sub_union_area) {
    warning(
      sprintf(
        "Detected %.6f km² of positive-area pairwise subbasin overlap.",
        overlap_area / 1e6
      ),
      call. = FALSE
    )
  }
  
  invisible(region_catalog(regions))
}

region_fractional_weights <- function(template, regions = NULL, equal_area_crs = EQUAL_AREA_CRS) {
  if (!inherits(template, "SpatRaster")) stop("template must be a terra SpatRaster.")
  if (terra::nlyr(template) != 1L) template <- template[[1]]
  if (is.null(regions)) regions <- load_nechako_regions(template = template)
  
  cell_area <- terra::cellSize(template, unit = "m")
  cell_area_v <- as.numeric(terra::values(cell_area, mat = FALSE))
  
  out <- lapply(names(regions), function(id) {
    v <- terra::project(regions[[id]], terra::crs(template))
    cover <- terra::rasterize(v, template, cover = TRUE)
    cover_v <- as.numeric(terra::values(cover, mat = FALSE))
    cover_v[!is.finite(cover_v)] <- 0
    cover_v <- pmax(0, pmin(1, cover_v))
    weights <- cell_area_v * cover_v
    list(
      region_id = id,
      region_name = gsub("_", " ", id, fixed = TRUE),
      coverage = cover_v,
      weights = weights,
      area_m2 = sum(weights, na.rm = TRUE),
      n_cells = sum(cover_v > 0, na.rm = TRUE),
      n_partial = sum(cover_v > 0 & cover_v < 1, na.rm = TRUE)
    )
  })
  names(out) <- names(regions)
  attr(out, "template_crs") <- terra::crs(template)
  attr(out, "cell_area_m2") <- cell_area_v
  out
}

## Fractional-coverage area-weighted mean (shared basin/regional mean rule).
##   ts_matrix    : cells in rows, time in columns
##   area_weights : cell area x fractional coverage, one per row (0 = outside)
##   min_coverage : minimum FRACTION (0-1) of the region's total weight that
##                  must carry finite values before a mean is reported
##   total_weight : the region's total weight; defaults to the sum of the
##                  positive supplied weights. Pass the full-region total when
##                  the rows are only a subset of the region's cells.
## A time step with no finite cell is always NA (also when min_coverage = 0).
.region_total_weight <- function(area_weights, total_weight = NULL) {
  if (is.null(total_weight)) {
    ok <- is.finite(area_weights) & area_weights > 0
    total_weight <- sum(area_weights[ok])
  }
  as.numeric(total_weight)
}

regional_valid_coverage <- function(ts_matrix, area_weights, total_weight = NULL) {
  if (!is.matrix(ts_matrix)) ts_matrix <- as.matrix(ts_matrix)
  if (length(area_weights) != nrow(ts_matrix)) stop("area_weights length must equal nrow(ts_matrix).")
  area_weights <- as.numeric(area_weights)
  tw <- .region_total_weight(area_weights, total_weight)
  if (!is.finite(tw) || tw <= 0) return(rep(NA_real_, ncol(ts_matrix)))
  w_ok <- is.finite(area_weights) & area_weights > 0
  vapply(seq_len(ncol(ts_matrix)), function(j) {
    ok <- w_ok & is.finite(ts_matrix[, j])
    sum(area_weights[ok]) / tw
  }, numeric(1))
}

regional_weighted_mean <- function(ts_matrix, area_weights,
                                   min_coverage = BASIN_MIN_VALID_AREA_COVERAGE,
                                   total_weight = NULL) {
  if (!is.matrix(ts_matrix)) stop("ts_matrix must be a matrix with pixels in rows and time in columns.")
  if (length(area_weights) != nrow(ts_matrix)) stop("area_weights length must equal nrow(ts_matrix).")
  if (!is.numeric(min_coverage) || length(min_coverage) != 1L || !is.finite(min_coverage) ||
      min_coverage < 0 || min_coverage > 1)
    stop("min_coverage must be a single fraction between 0 and 1.")
  area_weights <- as.numeric(area_weights)
  tw <- .region_total_weight(area_weights, total_weight)
  if (!is.finite(tw) || tw <= 0) return(rep(NA_real_, ncol(ts_matrix)))
  w_ok <- is.finite(area_weights) & area_weights > 0
  vapply(seq_len(ncol(ts_matrix)), function(j) {
    x <- ts_matrix[, j]
    ok <- w_ok & is.finite(x)
    sw <- sum(area_weights[ok])
    if (!is.finite(sw) || sw <= 0 || (sw / tw) < min_coverage) return(NA_real_)
    sum(x[ok] * area_weights[ok]) / sw
  }, numeric(1))
}

## Whole-basin fractional-coverage weights on a native raster grid.
## Basin = dissolved union of the 7 subbasins; mask = coverage > 0;
## weights = cell area (m2) x coverage (0 outside the basin).
basin_fractional_weights <- function(template, basin = NULL) {
  if (!inherits(template, "SpatRaster")) stop("template must be a terra SpatRaster.")
  template <- template[[1]]
  if (is.na(terra::crs(template)) || !nzchar(terra::crs(template)))
    stop("basin_fractional_weights(): template has no CRS.")
  if (is.null(basin)) basin <- nechako_subbasin_union(load_nechako_subbasins())
  if (!inherits(basin, "SpatVector")) basin <- terra::vect(basin)
  if (!terra::same.crs(basin, template)) basin <- terra::project(basin, terra::crs(template))
  cover_v <- as.numeric(terra::values(terra::rasterize(basin, template, cover = TRUE), mat = FALSE))
  cover_v[!is.finite(cover_v)] <- 0
  if (any(cover_v < -1e-6 | cover_v > 1 + 1e-6))
    stop("Basin coverage outside [0, 1]; refusing to clamp invalid coverage.")
  cover_v <- pmax(0, pmin(1, cover_v))
  cell_area_v <- as.numeric(terra::values(terra::cellSize(template, unit = "m"), mat = FALSE))
  w <- cell_area_v * cover_v
  w[!is.finite(w)] <- 0
  list(
    coverage     = cover_v,
    cell_area_m2 = cell_area_v,
    weights      = w,
    mask         = cover_v > 0,
    mask_raster  = terra::setValues(template, as.numeric(cover_v > 0)),
    area_m2      = sum(w),
    n_cells      = sum(cover_v > 0),
    n_partial    = sum(cover_v > 0 & cover_v < 1)
  )
}

## Basin mean of a cells x time matrix using basin_fractional_weights() output.
##   cells = NULL : rows are all ncell(template) cells in raster order
##   cells = idx  : rows correspond to raster cells `idx` (e.g. valid basin cells)
## Coverage is always measured against the WHOLE-basin weight (bw$area_m2).
basin_weighted_mean <- function(ts_matrix, bw, cells = NULL,
                                min_coverage = BASIN_MIN_VALID_AREA_COVERAGE) {
  if (!is.matrix(ts_matrix)) ts_matrix <- as.matrix(ts_matrix)
  w <- if (is.null(cells)) bw$weights else bw$weights[cells]
  regional_weighted_mean(ts_matrix, w, min_coverage = min_coverage, total_weight = bw$area_m2)
}

basin_valid_coverage <- function(ts_matrix, bw, cells = NULL) {
  if (!is.matrix(ts_matrix)) ts_matrix <- as.matrix(ts_matrix)
  w <- if (is.null(cells)) bw$weights else bw$weights[cells]
  regional_valid_coverage(ts_matrix, w, total_weight = bw$area_m2)
}

## Map lon/lat rows (e.g. from per-month index CSVs) to basin weights.
basin_weights_for_xy <- function(lon, lat, template, bw) {
  cells <- terra::cellFromXY(template[[1]], cbind(as.numeric(lon), as.numeric(lat)))
  w <- rep(0, length(cells))
  ok <- is.finite(cells)
  w[ok] <- bw$weights[cells[ok]]
  w
}

save_regional_index_csv <- function(ts_matrix, region_weights, dates, index_type,
                                    scale, out_dir = NECHAKO_REGION_OUTPUT_DIR) {
  if (!is.matrix(ts_matrix)) stop("ts_matrix must be a matrix.")
  dates <- as.Date(dates)
  validate_monthly_record(dates, label = "regional index dates", require_complete = TRUE)
  if (ncol(ts_matrix) != length(dates)) stop("ts_matrix/date length mismatch.")
  
  rows <- lapply(names(region_weights), function(id) {
    rw <- region_weights[[id]]
    if (length(rw$weights) != nrow(ts_matrix))
      stop("Region ", id, " weight length does not match ts_matrix rows.")
    value <- regional_weighted_mean(ts_matrix, rw$weights,
                                    min_coverage = BASIN_MIN_VALID_AREA_COVERAGE)
    data.frame(
      region_id = id,
      region_name = rw$region_name,
      date = dates,
      value = value,
      coverage_area_m2 = rw$area_m2,
      n_cells = rw$n_cells,
      n_partial_cells = rw$n_partial,
      stringsAsFactors = FALSE
    )
  })
  out <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(out_dir, sprintf("%s_%02d_regional_averaged_by_month.csv",
                                  tolower(index_type), as.integer(scale)))
  data.table::fwrite(out, f)
  invisible(out)
}

load_regional_index_csv <- function(index_type, scale,
                                    file = NULL,
                                    out_dir = NECHAKO_REGION_OUTPUT_DIR) {
  if (is.null(file))
    file <- file.path(out_dir, sprintf("%s_%02d_regional_averaged_by_month.csv",
                                       tolower(index_type), as.integer(scale)))
  if (!file.exists(file)) return(NULL)
  x <- data.table::fread(file)
  req <- c("region_id", "region_name", "date", "value")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Regional index file missing: ", paste(miss, collapse = ", "))
  x[, date := as.Date(date)]
  x[, value := as.numeric(value)]
  x
}

regional_wide_from_long <- function(x, value_col = "value") {
  x <- data.table::as.data.table(x)
  if (!all(c("region_id", "date", value_col) %in% names(x)))
    stop("regional_wide_from_long requires region_id, date and ", value_col)
  x[, year := as.integer(format(as.Date(date), "%Y"))]
  x[, month := as.integer(format(as.Date(date), "%m"))]
  x[, .(value = get(value_col)), by = .(region_id, year, month)]
}

save_region_catalog <- function(regions = NULL, out_dir = NECHAKO_REGION_OUTPUT_DIR) {
  if (is.null(regions)) regions <- load_nechako_regions()
  validate_nechako_regions(regions)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(region_catalog(regions), file.path(out_dir, "region_catalog.csv"))
  invisible(region_catalog(regions))
}