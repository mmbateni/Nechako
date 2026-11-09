####################################################################################
# DROUGHT_ANALYSIS_utils.R  --  SHARED UTILITIES FOR ALL NECHAKO DROUGHT SCRIPTS
#
# FRAMEWORK SCRIPTS (numbers = file-name prefixes; use these in all comments)
#   0   Monthly_ERALand_download.R           ERA5-Land monthly/daily inputs
#   1   1SPI_ERALand.R                       SPI
#   2b  2b_PET_ERALand.R                     FAO-56 PM PET + ERA5-Land PEV
#   3   3SPEI_ERALand.R                      SPEI-PM and SPEI-ERA5
#   4   4snow_drought_indices_ERALand.R      SWEI-3 and SSPI-1
#   5   5SRI_SDI_ERALand.R                   SRI and SDI (ERA5-Land runoff)
#   6   6_StreamFlow_Index_WSCLocal.R        station SSI (WSC gauges)
#   7   7_WBInex_NTSDI_GeoLakes_LocalRS.r    TSDI / NTSDI (lake and reservoir storage)
#   8   8_SSI_SSMI_ERALand.R                 soil moisture SSMI / SSI-SM (optional producer)
#   9   9SGI.R                               groundwater SGI (optional producer)
#   10  10MSDI.r                             multivariate MSDI (optional producer)
#   11  11pr_pet_trends.r                    P / PET / PEV / T / AET / SWE trends
#   12  12pr_pet_trends_visualization.r      figures for Script 11
#   13  13trend_test_ALL.R                   drought-index trend tests
#   14  14basin_timeseries.R                 basin time series and summaries
#   15  15trends_visualization.R             figures for Script 13
#   16  16event_ranking.R                    drought-event ranking and return periods
#
# EXECUTION ORDER: 0 -> 1 -> 2b -> 3 -> 4 -> 5 -> 6 -> 7 -> 8 (optional) -> 9/10 (optional) -> 11 -> 12 -> 13 -> 14 -> 15 -> 16
# Every framework script starts with rm(list = ls()) and WD_PATH <- getwd();
# run each one from the project root.
#
# COMPANION ANALYSES (run after the framework; prerequisites in brackets)
#   20a 20a_dynamic_thermodynamic_decomp.R   legacy SPEI_PM vs SPEI_ERA5 diagnostics [3, 13]
#   20b 20b_PET_bias_nonstationarity.R       PET_PM - ERA5-Land PEV difference      [2b, 11]
#   22  22PET_formulation_sensitivity.R      SPEI_PM vs SPEI_ERA5 agreement metrics  [3, 13]
#   23  23meteorological_controls.R          meteorological controls of |dSPEI|      [2b, 22]
#   24  24circulation_context.R              boundary IVT, MFC, precipitation        [0]
#   25  25PET_shapley_decomposition.R        Shapley allocation of PM PET anomalies  [2b, 11, 22]
#   26  26_climatology_water_budget.R        regional climatology and water budget   [11]
# COMPANION ORDER: 22 -> 23 -> 25; 20a, 20b, 24 and 26 are independent of each other.
# Companion scripts set the project root with
#   setwd(Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/"))
# before sourcing this file, and use the shared settings and helpers in
# Section N (record and recent periods, warm season, date parsing, regional
# CSV reader, bootstrap indices, HAC trend test).
# Scripts named elsewhere in comments with other numbers (e.g. w9_*) are
# external analyses and are NOT part of this framework.
#
# SECTIONS
#   A.  Project configuration (paths, constants, thresholds)
#   B.  Basin-averaged CSV loader / writer / package bootstrap
#   C.  NetCDF file finders
#   D.  NetCDF date parsing
#   E.  Basin spatial helpers (shared with K2 regional framework)
#   F.  Drought event detection (single framework-wide definition)
#   F2. Spatial helpers (rasterise / smooth / plot)
#   G.  ggplot helpers
#   J.  [SHARED STATS]  Functions used by BOTH 13trend_test_ALL.R and
#       11pr_pet_trends.r --- moved here to eliminate duplication:
#         J1. Utility helpers  (check_min_value_threshold,
#                               calculate_variance_with_ties)
#         J2. PELT changepoint  (detect_regime_shift_pelt)
#         J3. Shared VC-MK guard rails (vc_mmkh)
#         J4. TFPW-MK per-series (.mk_tfpw_single)
#         J5. Spectral peak detection per-series (.spectral_single)
#         J6. Omnibus single-series wrapper (mk_tfpw_spectral_for_series)
#             --- used directly by 11pr_pet_trends.r
#         J7. Vectorised pixel-matrix wrappers  (vectorized_mann_kendall,
#             vectorized_mann_kendall_tfpw, vectorized_regime_shift_pelt,
#             vectorized_spectral_peaks)
#             --- used directly by 13trend_test_ALL.R
#   K.  Zero-inflated gamma standardization, Stagge et al. (2015) (SPI, SSPI, SRI)
#   L.  Monthly water-depth units: metadata + magnitude check and conversion
#   M.  Distribution selection by L-moment-ratio (tau4) distance (SPEI, SSI)
#   K2. Nechako regional spatial framework (7 subbasins + whole basin)
#   N.  Companion-analysis helpers (scripts 20a-26)
#   O. Basin series written by the producing scripts.
####################################################################################

## A. PROJECT CONFIGURATION ----------------------------------------
## Project root = the working directory of the calling script (framework rule:
## every script sets WD_PATH <- getwd() before sourcing this file).
## Framework scripts set WD_PATH before sourcing; companion scripts only setwd().
.wd_now <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!exists("WD_PATH", inherits = TRUE)) {
  WD_PATH <- .wd_now
} else if (!identical(normalizePath(WD_PATH, winslash = "/", mustWork = FALSE), .wd_now)) {
  stop("WD_PATH (", WD_PATH, ") differs from the working directory (", .wd_now,
       "). Run every script from the project root.")
}
rm(.wd_now)
BASIN_SHP       <- file.path(WD_PATH, "Spatial/Nechako_shps")

## Helper: load the authoritative Nechako basin as a single dissolved
## SpatVector from the supplied subbasin Shapefiles.  The source directory is
## intentionally used instead of a legacy KMZ/KML boundary.
## Single authoritative basin loader. Subbasins and their union are cached per
## directory, so the 7 Shapefiles are read once per R session.
.NECHAKO_GEOM_CACHE <- new.env(parent = emptyenv())
load_basin_vect <- function(path = BASIN_SHP) {
  key <- paste0("union|", normalizePath(path, winslash = "/", mustWork = FALSE))
  if (!is.null(.NECHAKO_GEOM_CACHE[[key]])) return(.NECHAKO_GEOM_CACHE[[key]])
  if (!dir.exists(path))
    stop("Nechako subbasin Shapefile directory not found: ", path)
  u <- nechako_subbasin_union(load_nechako_subbasins(path))
  assign(key, u, envir = .NECHAKO_GEOM_CACHE)
  u
}

## Index-specific directories
SPI_SEAS_DIR      <- file.path(WD_PATH, "spi_results_seasonal/")
SPEI_SEAS_DIR     <- file.path(WD_PATH, "spei_results_seasonal/")
SPEI_ERA5_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal_era5land/")   # ERA5-Land native PEV outputs
SWEI_SEAS_DIR     <- file.path(WD_PATH, "swei_results_seasonal/")
SSPI_SEAS_DIR     <- file.path(WD_PATH, "sspi_results_seasonal/")

## Output directories
TREND_DIR      <- file.path(WD_PATH,  "temporal_drought/")
BASIN_PLOT_DIR <- file.path(TREND_DIR, "basin_averaged_plots")  # authoritative Script 14 output directory
POINT_PLOT_DIR <- file.path(WD_PATH,  "point_timeseries_plots/")
CACHE_DIR      <- file.path(WD_PATH,  "temporal_drought/cache/")

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

## Shared analysis periods for the companion scripts (20a-26). Scripts alias
## these locally (e.g. ANALYSIS_START <- NECHAKO_RECORD_START) so every
## companion analysis uses the same record, recent period and warm season.
NECHAKO_RECORD_START       <- as.Date("1950-01-01")
NECHAKO_RECORD_END         <- as.Date("2025-12-01")
NECHAKO_RECENT_START       <- as.Date("2022-01-01")
NECHAKO_RECENT_END         <- as.Date("2025-12-01")
NECHAKO_WARM_SEASON_MONTHS <- 5:9   # May-September (MJJAS)

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


## Timescales --- current manuscript source of truth
MANUSCRIPT_SCALES    <- c(1L, 3L, 6L, 12L, 24L, 36L)
EVENT_RANKING_SCALES <- MANUSCRIPT_SCALES
SPATIAL_PLOT_SCALES  <- c(1L, 3L, 6L, 12L, 24L, 36L)

SPI_SCALES          <- MANUSCRIPT_SCALES
SPEI_SCALES         <- MANUSCRIPT_SCALES
SPEI_ERA5_SCALES    <- MANUSCRIPT_SCALES
SRI_SCALES          <- MANUSCRIPT_SCALES
SWEI_SCALE          <- 3L   
SSPI_SCALE          <- 1L   
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

## A3. SHARED ANALYSIS CONSISTENCY CONTRACTS (Scripts 11-16) ---------------
## Statistical settings, annual-series validation, unit QA and output isolation.
## These are defined here so Scripts 11-16 need no companion configuration file.
ANALYSIS_MIN_HOMOG_YEARS <- 30L
ANALYSIS_HOMOG_MC <- 1000L
ANALYSIS_RECENT_START <- NECHAKO_RECENT_START
ANALYSIS_RECENT_END <- seq(NECHAKO_RECENT_END, by = "month", length.out = 2L)[2L] - 1L
ANALYSIS_PERIOD_START <- as.integer(format(NECHAKO_RECORD_START, "%Y"))
ANALYSIS_PERIOD_END <- as.integer(format(NECHAKO_RECORD_END, "%Y"))
ANALYSIS_MANUSCRIPT_SCALES <- as.integer(MANUSCRIPT_SCALES)
## Shared trend / homogeneity settings for Scripts 11-16 -----------------------
ANALYSIS_STEP_YEAR    <- 1979L  # exploratory fixed-date sensitivity split; NOT independently pre-specified
ANALYSIS_MIN_TREND_N  <- 12L    # minimum finite observations for any MK / Sen / TFPW test
ANALYSIS_TIE_MODE_MAX <- 0.50   # filter a series when its largest tied group exceeds this share
ANALYSIS_PERIOD_LABEL       <- sprintf("%d\u2013%d", ANALYSIS_PERIOD_START, ANALYSIS_PERIOD_END)
ANALYSIS_PERIOD_LABEL_ASCII <- sprintf("%d-%d",      ANALYSIS_PERIOD_START, ANALYSIS_PERIOD_END)

## Pixel-to-basin correlation mapped by Script 14 (Fig4c) and Script 15 (Part 8):
## a column of the {idx}_{sc}_pixel_basin_correlation.csv files written by Script 13.
ANALYSIS_PIXEL_BASIN_CORR <- "spearman_rho" #or "pearson_r"   
ANALYSIS_PIXEL_BASIN_CORR_LABEL <- switch(ANALYSIS_PIXEL_BASIN_CORR,
                                          pearson_r = "Pearson r", spearman_rho = "Spearman \u03c1",
                                          stop("ANALYSIS_PIXEL_BASIN_CORR must be 'pearson_r' or 'spearman_rho'."))
ANALYSIS_PIXEL_BASIN_CORR_TAG <- switch(ANALYSIS_PIXEL_BASIN_CORR,
                                        pearson_r = "Pearson_r", spearman_rho = "Spearman_rho")

## Share of the largest group of tied values (shared tie rule for Scripts 11 and 13).
analysis_modal_tie_fraction <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(tabulate(match(x, unique(x)))) / length(x)
}
# Both batteries remove nonfinite observations and retain their ORIGINAL years.
# Serial dependence across a missing calendar year is not silently treated as
# adjacent by AR(1)-GLS; consumers must use actual years in the correlation index.
analysis_homog_clean <- function(x, years, min_n = ANALYSIS_MIN_HOMOG_YEARS) {
  stopifnot(length(x) == length(years))
  ok <- is.finite(x) & is.finite(years)
  x <- as.numeric(x[ok]); years <- as.integer(years[ok])
  if (length(years) && (anyDuplicated(years) || is.unsorted(years, strictly = TRUE)))
    stop("Homogeneity years must be strictly increasing and unique")
  list(x = x, years = years, enough = length(x) >= min_n)
}

# trend::snh.test/br.test report a Monte Carlo p-value, not the exceedance
# count. When it is on the expected b/B grid, reconstruct b and apply the
# (b+1)/(B+1) estimator. Off-grid values are retained except that zero is
# bounded away from zero; an exact plus-one correction is then not identifiable.
analysis_mc_p <- function(returned_p, B = ANALYSIS_HOMOG_MC) {
  if (length(returned_p) != 1L || !is.finite(returned_p)) return(NA_real_)
  if (B < 1L || B != as.integer(B)) stop("Invalid Monte Carlo replicate count")
  p <- as.numeric(returned_p)
  if (p < 0 || p > 1) stop("Monte Carlo p-value outside [0, 1]")
  b <- round(p * B)
  if (abs(p * B - b) <= 1e-7) return((b + 1) / (B + 1))
  max(p, 1 / (B + 1))
}

analysis_check_output_isolation <- function(dir1, dir2) {
  norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  a <- norm(dir1); b <- norm(dir2)
  if (identical(a,b)) stop("Output collision: climate and drought trend roots coincide: ", a)
  invisible(TRUE)
}

analysis_unit_values <- function(r) {
  u <- unique(trimws(as.character(terra::units(r))))
  unique(tolower(u[!is.na(u) & nzchar(u)]))
}
analysis_verify_pm_control_units <- function(r, variable) {
  # Reject absent or unexpected units; guessing scales would silently corrupt
  # trends. User must correct the upstream export instead.
  acceptable <- list(Temperature = c("c", "degc", "degree_celsius", "degrees_celsius", "degree celsius", "celsius", "°c"),
                     NetRadiation = c("mj/m2/day", "mj m-2 day-1", "mj m^-2 day^-1", "mj m-2 d-1", "mj/m^2/day"),
                     VPD = c("kpa"), Wind2m = c("m/s", "m s-1", "m s^-1"),
                     SurfacePressure = c("kpa"))
  if (!variable %in% names(acceptable)) stop("Unknown PM control: ", variable)
  u <- analysis_unit_values(r)
  if (length(u) != 1L || !u %in% acceptable[[variable]])
    stop("[PM UNIT] ", variable, ": expected ", paste(acceptable[[variable]], collapse = " or "),
         "; found ", if (length(u)) paste(u, collapse = "; ") else "<missing>",
         ". Fix upstream units; conversion is not inferred.")
  invisible(u)
}
analysis_negative_qa <- function(r, variable, unit, months = NULL, tol = NA_real_) {
  # One source for pre-clamp and post-clamp negative-value statistics.
  # When months is NULL, retain the original compact output contract.
  # Ordinary post-clamp QA remains disk-backed; only the detailed, month-
  # stratified audit needs materialized per-layer values.
  if (is.null(months)) {
    n <- sum(terra::global(!is.na(r), "sum", na.rm = TRUE)[, 1], na.rm = TRUE)
    k <- sum(terra::global(r < 0, "sum", na.rm = TRUE)[, 1], na.rm = TRUE)
    lo <- if (k) min(terra::global(terra::ifel(r < 0, r, NA_real_),
                                   "min", na.rm = TRUE)[, 1], na.rm = TRUE) else NA_real_
    return(data.frame(variable = variable, unit = unit, n_values = n,
                      n_negative = k, pct_negative = 100 * k / max(n, 1),
                      most_negative = lo))
  }
  v <- terra::values(r, mat = FALSE)
  ok <- is.finite(v)
  neg <- which(ok & v < 0)
  n_ok <- sum(ok)
  out <- data.frame(
    variable = variable, unit = unit, n_values = n_ok,
    n_negative = length(neg), pct_negative = 100 * length(neg) / max(n_ok, 1L),
    most_negative = if (length(neg)) min(v[neg]) else NA_real_)
  if (length(months) != terra::nlyr(r))
    stop("analysis_negative_qa: months must match the raster's layer count")
  layer_month <- if (length(neg)) months[(neg - 1L) %/% terra::ncell(r) + 1L] else integer(0)
  out$p01_negative <- if (length(neg)) unname(stats::quantile(v[neg], 0.01)) else NA_real_
  out$tolerance <- tol
  out$n_beyond_tolerance <- if (is.finite(tol)) sum(v[neg] < -tol) else NA_integer_
  out$pct_negative_in_JunOct <- if (length(neg)) 100 * mean(layer_month %in% 6:10) else NA_real_
  out
}

# Producer-specific depth QA. Scripts 1/2b supply depth rates in mm/day,
# whereas resolve_monthly_depth_unit() handles ambiguous MONTHLY NetCDF units;
# these separate unit contracts must not silently convert each other.
analysis_validate_daily_depth <- function(r, label, months, log_fun,
                                          expected_min = 0.01, expected_max = 20) {
  july_idx <- which(as.integer(months) == 7L)[1L]
  if (!is.finite(july_idx))
    stop(sprintf("[%s] No July layer available for unit QA.", label))
  v <- terra::values(r[[july_idx]], mat = FALSE)
  v <- v[is.finite(v)]
  if (!length(v)) stop(sprintf("[%s] July layer has no finite values.", label))
  q <- stats::quantile(v, probs = c(0.01, 0.50, 0.99),
                       na.rm = TRUE, names = FALSE)
  log_fun(sprintf("  [UNIT:%s] July q01/median/q99 = %.5g / %.5g / %.5g",
                  label, q[1], q[2], q[3]))
  if (q[2] < expected_min)
    stop(sprintf("[%s] Median July %.6g below expected mm/day range; check upstream NetCDF (possible m/day).",
                 label, q[2]))
  if (q[2] > expected_max)
    stop(sprintf("[%s] July median %.6g exceeds configured maximum %.6g mm/day; check upstream NetCDF units.",
                 label, q[2], expected_max))
  log_fun(sprintf("  [UNIT:%s] Magnitude consistent with configured mm/day QA range.", label))
  log_fun(sprintf("  [UNIT:%s] No 1000x conversion applied; monthly totals use calendar month lengths.", label))
  invisible(q)
}

# Canonical P vs PM-PET monthly climatology for Scripts 12 and 14.
# Input is Script 11's basin-area weighted monthly means, in mm/month.
# Pairwise finite months keep both variables on the same year/month support.
analysis_monthly_p_pet_climatology <- function(monthly = NULL) {
  if (is.null(monthly)) {
    f <- file.path(WD_PATH, "trend_analysis_pr_pet", "analysis_metadata.rds")
    if (!file.exists(f)) stop("Missing Script 11 basin-monthly metadata: ", f)
    monthly <- readRDS(f)$basin_avg_monthly
  }
  required <- c("date", "precip_mm_month", "pet_mm_month")
  if (is.null(monthly) || !all(required %in% names(monthly)))
    stop("Script 11 basin-monthly metadata lacks: ", paste(setdiff(required, names(monthly)), collapse = ", "))
  d <- data.frame(date = as.Date(monthly$date), P = as.numeric(monthly$precip_mm_month),
                  PET_PM = as.numeric(monthly$pet_mm_month))
  d <- d[!is.na(d$date) & is.finite(d$P) & is.finite(d$PET_PM), , drop = FALSE]
  if (!nrow(d)) stop("No simultaneous finite basin P and PM-PET months.")
  d$month <- as.integer(format(d$date, "%m"))
  d$year <- as.integer(format(d$date, "%Y"))
  summary_one <- function(m) {
    x <- d[d$month == m, , drop = FALSE]
    c(P_Mean = if (nrow(x)) mean(x$P) else NA_real_,
      P_SD = if (nrow(x) > 1L) stats::sd(x$P) else NA_real_,
      PET_Mean = if (nrow(x)) mean(x$PET_PM) else NA_real_,
      PET_SD = if (nrow(x) > 1L) stats::sd(x$PET_PM) else NA_real_,
      n_yrs = length(unique(x$year)))
  }
  out <- data.frame(month = 1:12, do.call(rbind, lapply(1:12, summary_one)))
  if (any(out$n_yrs == 0L)) stop("Climatology is missing a calendar month.")
  out$Month <- factor(month.abb[out$month], levels = month.abb)
  out$Balance <- out$P_Mean - out$PET_Mean
  attr(out, "years") <- range(d$year)
  out
}

## B. BASIN-AVERAGED CSV LOADER ----------------------------------------
#' Load a pre-computed basin-averaged drought index CSV and return a tidy
#' long-format data.frame(date, value).
#'
#' AUTHORITATIVE SOURCE: each index script writes them via write_basin_series().
#'   The methodology is:
#'   (1) compute the drought index at every basin pixel from per-pixel inputs,
#'   (2) take the area-weighted spatial mean of those pixel-level index values.
#' This is the standard hydrometeorological approach: area-weighted averaging
#' of standardised index values that were calibrated on the same spatial basis.
#' Run 13trend_test_ALL.R before 14basin_timeseries.R.
#'
#' File naming convention: {index}_{scale:02d}_basin_averaged_by_month.csv
#'   where index may be spi, spei (FAO-56 PM), spei_era5 (ERA5-Land native PEV), or swei
#' File structure:
#'   Year | Jan | Feb | Mar | --- | Dec
#'   1950 | -1.88 |  NA | ---
#'
#' @param index_type  "spi", "spei", "spei_era5", or "swei"
#' @param scale       Integer timescale (e.g. 1, 3, 6, 12)
#' @param data_dir    Directory containing the CSV; if NULL, the index-specific
#'                    seasonal directory from section A is used automatically.
#' @return data.frame with columns \code{date} (Date) and \code{value} (numeric),
#'         sorted by date with missing values retained as NA. Returns NULL with a warning if
#'         the file is not found (the corresponding index producer may not have run).
load_basin_avg_csv <- function(index_type, scale, data_dir = NULL) {
  dir_map <- list(spi      = SPI_SEAS_DIR,
                  spei     = SPEI_SEAS_DIR,
                  spei_era5 = SPEI_ERA5_SEAS_DIR,
                  swei     = SWEI_SEAS_DIR,
                  sspi     = SSPI_SEAS_DIR,
                  sri  = file.path(WD_PATH, "sri_results_seasonal", "runoff"),
                  ssmi = file.path(WD_PATH, "ssmi_results_seasonal"),
                  ssi_sm = file.path(WD_PATH, "ssi_sm_results_seasonal"))
  dir <- if (!is.null(data_dir)) data_dir else {
    idx <- tolower(index_type)
    d <- dir_map[[idx]]
    if (is.null(d)) stop("Unknown index_type: ", index_type)
    if (idx %in% c("ssmi", "ssi_sm")) {
      cfg_file <- file.path(TREND_DIR, "script13_analysis_config.rds")
      if (!file.exists(cfg_file))
        stop("Soil series requires ", cfg_file,
             " to resolve the active SOIL_LAYER; pass data_dir explicitly otherwise.")
      cfg <- readRDS(cfg_file)
      layer <- cfg$soil_layer
      if (length(layer) != 1L || is.na(layer) || !grepl("^L[0-9]+_[0-9]+$", layer))
        stop("Invalid soil_layer in ", cfg_file)
      d <- file.path(d, layer)
    }
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
  if (length(mon_cols) != 12L)
    warning(basename(f), " is missing month column(s): ",
            paste(setdiff(MONTH_NAMES, mon_cols), collapse = ", "), call. = FALSE)
  long <- do.call(rbind, lapply(mon_cols, function(mc) {
    data.frame(
      date  = as.Date(sprintf("%04d-%02d-01", as.integer(df$Year),
                              match(mc, MONTH_NAMES))),
      value = as.numeric(df[[mc]]),
      stringsAsFactors = FALSE)
  }))
  long <- long[order(long$date), ]
  long$value[!is.finite(long$value)] <- NA_real_
  rownames(long) <- NULL
  long
}

## B2. BASIN-AVERAGED CSV WRITER ----------------------------------------
#' Compute an area-weighted basin-average time series from a pixel x month
#' matrix and write the result as the standard basin_averaged_by_month.csv
#' expected by load_basin_avg_csv().
#'
#' Optional writer for index-producing scripts after the pixel time-series
#' matrix has been assembled. Script 13 reads producer-written basin series;
#' it does not call this helper. Area weights follow the producer's native-grid
#' fractional coverage and cell-area contract; the raster is not reprojected.
#'
#' Methodology:
#'   basin_avg[t] = sum( index_value[pixel, t] x cell_area[pixel] )
#'                / sum( cell_area[pixel] )          (sum over non-NA pixels)
#'
#' @param ts_matrix   Numeric matrix (n_valid_pixels x n_months).  Rows
#'                    correspond to valid basin pixels; columns to calendar
#'                    months in chronological order.
#' @param area_weights Numeric vector of length n_valid_pixels; cell areas in
#'                    m2 from native-grid cell area times fractional coverage.
#' @param dates       Date vector of length n_months giving the first day of
#'                    each month represented by a column of ts_matrix.
#' @param index_type  "spi", "spei", "swei", or any other index key.
#' @param scale       Integer timescale (e.g. 1, 3, 6, 12).
#' @param out_dir     Directory to write the CSV.  Should be the index-specific
#'                    seasonal results directory (SPI_SEAS_DIR etc.) so that
#'                    load_basin_avg_csv() can find it without specifying data_dir.
#' @return Invisible data.frame (Year x month-name columns) as written to disk.
save_basin_avg_from_pixels <- function(ts_matrix, area_weights, dates,
                                       index_type, scale, out_dir,
                                       min_coverage = BASIN_MIN_VALID_AREA_COVERAGE,
                                       total_weight = NULL) {
  ## area_weights must be cell area x fractional basin coverage (see
  ## basin_fractional_weights()); total_weight is the whole-basin weight
  ## (defaults to the sum of the supplied weights, i.e. rows = all basin cells).
  ## --- Input validation ----------------------------------------
  if (!is.matrix(ts_matrix))
    stop("ts_matrix must be a numeric matrix (pixels x months)")
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
  
  ## --- Fractional-coverage area-weighted mean with the shared coverage rule ---
  basin_avg <- regional_weighted_mean(ts_matrix, area_weights,
                                      min_coverage = min_coverage,
                                      total_weight = total_weight)
  
  ## --- Reshape to Year x Month wide format ----------------------------------------
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
  
  ## --- Write CSV ----------------------------------------
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(
    out_dir,
    sprintf("%s_%02d_basin_averaged_by_month.csv",
            tolower(index_type), as.integer(scale)))
  utils::write.csv(df, out_file, row.names = FALSE)
  cat(sprintf(
    "  [OK] Basin-avg CSV written (area-weighted, %d pixels, %d months): %s\n",
    nrow(ts_matrix), ncol(ts_matrix), basename(out_file)))
  invisible(df)
}

## B3. PACKAGE BOOTSTRAP ----------------------------------------
utils_load_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org/")
  }
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
  cat("[OK] Packages loaded:", paste(pkgs, collapse = ", "), "\n")
}

####################################################################################
# C. NETCDF FILE FINDERS ----------------------------------------
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

## D. NETCDF DATE PARSING ----------------------------------------
extract_dates_from_nc <- function(nc_file, n_layers = NULL) {
  r_temp <- terra::rast(nc_file)
  n_layers <- if (is.null(n_layers)) terra::nlyr(r_temp) else n_layers
  # Seasonal NetCDFs hold one layer PER YEAR at the filename's calendar month.
  month_match <- regmatches(basename(nc_file),
                            regexec("_month([0-9]{2})_", basename(nc_file),
                                    ignore.case = TRUE))[[1L]]
  seasonal_month <- if (length(month_match)) as.integer(month_match[2L]) else NA_integer_
  if (!is.na(seasonal_month) && !(seasonal_month %in% 1:12))
    stop("Invalid calendar month in NetCDF filename: ", basename(nc_file))
  
  dates <- tryCatch({
    d <- terra::time(r_temp)
    if (!is.null(d) && length(d) > 0L) {
      d_conv <- tryCatch(as.Date(d), error = function(e) as.Date(d, origin = "1970-01-01"))
      if (all(!is.na(d_conv)) &&
          min(d_conv) >= as.Date("1900-01-01") &&
          max(d_conv) <= as.Date("2030-01-01")) {
        d_conv
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  if (!is.na(seasonal_month) && !is.null(dates)) {
    year <- as.integer(format(dates, "%Y"))
    month <- as.integer(format(dates, "%m"))
    if (length(dates) != n_layers || anyNA(dates) ||
        any(month != seasonal_month) || is.unsorted(year, strictly = TRUE))
      stop("NetCDF time axis does not match annual month", sprintf("%02d", seasonal_month),
           " file: ", basename(nc_file))
  }
  if (!is.null(dates)) {
    if (is.na(seasonal_month)) return(.pad_dates(dates, n_layers))
    return(dates)
  }
  
  if (!is.na(seasonal_month)) {
    # Without a time coordinate, shorter seasonal stacks have an unknown year
    # origin and must not be silently assigned to 1950 onward.
    expected_years <- ANALYSIS_PERIOD_END - ANALYSIS_PERIOD_START + 1L
    if (n_layers != expected_years)
      stop("Cannot infer years for ", basename(nc_file), ": ", n_layers,
           " annual layers without a valid time axis (expected ", expected_years,
           " for ", ANALYSIS_PERIOD_START, "-", ANALYSIS_PERIOD_END, ").")
    cat(sprintf("  [dates] Falling back to annual month %02d (%d-%d): %s\n",
                seasonal_month, ANALYSIS_PERIOD_START, ANALYSIS_PERIOD_END,
                basename(nc_file)))
    return(seq(as.Date(sprintf("%04d-%02d-01", ANALYSIS_PERIOD_START,
                               seasonal_month)), by = "year", length.out = n_layers))
  }
  
  cat("  [dates] Using fallback monthly sequence (Jan 1950)\n")
  d <- seq(NECHAKO_RECORD_START, by = "month", length.out = n_layers)
  if (tail(d, 1L) > NECHAKO_RECORD_END) {
    end <- NECHAKO_RECORD_END
    start <- seq(end, by = "-1 month", length.out = n_layers)[n_layers]
    d <- seq(start, by = "month", length.out = n_layers)
  }
  d
}

.pad_dates <- function(d, n) {
  if (length(d) == n) return(d)
  if (length(d) > n) return(d[1:n])
  extra <- seq(tail(d, 1), by = "month", length.out = n - length(d) + 1L)[-1L]
  c(d, extra)
}

## E. BASIN SPATIAL HELPERS ----------------------------------------
# Native cell-area x fractional basin-coverage utilities are defined in section K2.
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

## K. ZERO-INFLATED GAMMA STANDARDIZATION — Stagge et al. (2015) ###############
## Used by SPI (Script 1), SSPI (Script 4) and SRI (Script 5).
##   p0      = n_zero / (n + 1)                 probability mass of zeros
##   p(x>0)  = p0 + (1 - p0) * G(x)             G = gamma CDF, L-moment fit
##   p(x=0)  = (n_zero + 1) / (2 * (n + 1))     centre of mass of the zero block
## Assigning the centre of mass (instead of the full p0) to zeros avoids
## pushing dry months towards near-normal values. With no zeros in the
## reference sample the method reduces to the plain L-moment gamma fit.
stagge_zero_probabilities <- function(n_zero, n_total) {
  if (!is.finite(n_total) || n_total < 1) stop("n_total must be >= 1.")
  list(p0     = n_zero / (n_total + 1),
       p_zero = (n_zero + 1) / (2 * (n_total + 1)))
}

#' Non-exceedance probabilities for `obs` from a zero-inflated gamma fitted
#' (L-moments, lmomco) to the positive values of `ref`. Returns NULL when the
#' fit is not possible (too few positive values, failed fit, non-finite CDF),
#' so the caller can use its empirical fallback. Non-finite obs stay NA.
zero_inflated_gamma_prob <- function(obs, ref, min_pos = 10L) {
  ref <- ref[is.finite(ref)]
  nr  <- length(ref)
  if (nr < 1L) return(NULL)
  pos <- ref[ref > 0]
  if (length(pos) < min_pos) return(NULL)
  par <- tryCatch(lmomco::pargam(lmomco::lmoms(pos)), error = function(e) NULL)
  if (is.null(par) || any(!is.finite(par$para)) ||
      (!is.null(par$ifail) && any(par$ifail != 0))) return(NULL)
  sp  <- stagge_zero_probabilities(sum(ref <= 0), nr)
  p   <- rep(NA_real_, length(obs))
  fin <- is.finite(obs)
  ip  <- which(fin & obs > 0)
  if (length(ip)) {
    Fg <- tryCatch(as.numeric(lmomco::cdfgam(obs[ip], par)), error = function(e) NULL)
    if (is.null(Fg) || length(Fg) != length(ip) || any(!is.finite(Fg))) return(NULL)
    p[ip] <- sp$p0 + (1 - sp$p0) * Fg
  }
  p[fin & obs <= 0] <- sp$p_zero
  p
}


## L. MONTHLY WATER-DEPTH UNITS (precipitation, PET, PEV) ######################
## Units are taken from the NetCDF `units` attribute AND inferred from
## magnitude. The run stops if the two disagree, so a producer switching
## between mm/day and mm/month can never be silently double-converted.
## Inference: the record-mean monthly value is turned into an annual total
## under each candidate unit; exactly one candidate must fall inside the
## plausible range for the variable.
MONTHLY_DEPTH_UNITS <- c("mm/day", "m/day", "mm/month", "m/month")
.DEPTH_ANNUAL_FACTOR <- c(`mm/day` = 365.25, `m/day` = 365250, `mm/month` = 12, `m/month` = 12000)

## USER SETTING: plausible basin/domain-mean ANNUAL totals (mm/yr).
UNIT_PLAUSIBLE_ANNUAL_MM <- list(
  precipitation = c(100, 5000),
  pet           = c(150, 3000),    # FAO-56 Penman-Monteith reference ET
  pev           = c(150, 15000)    # ERA5-Land potential evaporation (very large here)
)

#' Map a free-text units attribute to one of MONTHLY_DEPTH_UNITS (or NA).
#' Bare "m"/"mm" without a time denominator are ambiguous and return NA.
normalize_depth_unit <- function(u) {
  if (is.null(u) || !length(u) || is.na(u[1]) || !nzchar(trimws(u[1]))) return(NA_character_)
  s <- gsub("\\s+", " ", tolower(trimws(as.character(u[1]))))
  lead <- sub("^([a-z]+).*$", "\\1", s)
  depth <- if (lead %in% c("mm", "millimeter", "millimetre", "millimeters", "millimetres")) "mm" else
    if (lead %in% c("m", "meter", "metre", "meters", "metres")) "m" else NA_character_
  if (is.na(depth)) return(NA_character_)
  per_day   <- grepl("(/|per )\\s*(day|d)\\b|\\bday-1\\b|\\bd-1\\b|24.?hour|\\b24 ?h\\b", s, perl = TRUE)
  per_month <- grepl("(/|per )\\s*(month|mon)\\b|\\bmonth-1\\b|\\bmon-1\\b", s, perl = TRUE)
  if (per_day && !per_month) return(paste0(depth, "/day"))
  if (per_month && !per_day) return(paste0(depth, "/month"))
  NA_character_
}

#' Read the units attribute of the (single) data variable of a NetCDF file.
read_netcdf_units <- function(file, var = NULL) {
  if (is.null(file) || !file.exists(file)) return(NA_character_)
  nc <- ncdf4::nc_open(file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  if (is.null(var)) {
    vars <- names(nc$var)[vapply(nc$var, function(v) v$ndims > 0, logical(1))]
    if (length(vars) != 1L) return(NA_character_)
    var <- vars[[1L]]
  }
  att <- ncdf4::ncatt_get(nc, var, "units")
  if (isTRUE(att$hasatt)) trimws(as.character(att$value)) else NA_character_
}

#' Infer the unit of a monthly depth raster from its magnitude.
infer_monthly_depth_unit <- function(r, kind, cells = NULL) {
  rng <- UNIT_PLAUSIBLE_ANNUAL_MM[[kind]]
  if (is.null(rng)) stop("Unknown kind '", kind, "'; add it to UNIT_PLAUSIBLE_ANNUAL_MM.")
  if (is.null(cells)) {
    lm <- terra::global(r, "mean", na.rm = TRUE)[, 1]
  } else {
    v  <- terra::values(r, mat = TRUE)[cells, , drop = FALSE]
    lm <- colMeans(v, na.rm = TRUE)
  }
  m <- mean(lm[is.finite(lm)])
  totals <- m * .DEPTH_ANNUAL_FACTOR
  ok <- names(totals)[is.finite(totals) & totals >= rng[1] & totals <= rng[2]]
  list(unit = if (length(ok) == 1L) ok else NA_character_, candidates = ok,
       mean_value = m, annual_totals = totals, range = rng)
}

#' Decide the unit from metadata + magnitude; stop on any conflict.
resolve_monthly_depth_unit <- function(r, kind, file = NULL, var = NULL, cells = NULL,
                                       label = kind) {
  meta_raw <- read_netcdf_units(file, var)
  if (is.na(meta_raw)) {
    tu <- tryCatch(terra::units(r)[1], error = function(e) NA_character_)
    if (length(tu) && !is.na(tu) && nzchar(tu)) meta_raw <- tu
  }
  meta <- normalize_depth_unit(meta_raw)
  inf  <- infer_monthly_depth_unit(r, kind, cells)
  cat(sprintf("\n===== UNIT CHECK: %s =====\n", label))
  cat(sprintf("  Metadata units       : %s -> %s\n",
              ifelse(is.na(meta_raw), "<missing>", meta_raw), ifelse(is.na(meta), "<unrecognised>", meta)))
  cat(sprintf("  Record-mean value    : %.6g per layer\n", inf$mean_value))
  cat(sprintf("  Implied annual total : %s (plausible %g-%g mm/yr)\n",
              paste(sprintf("%s=%.0f", names(inf$annual_totals), inf$annual_totals), collapse = ", "),
              inf$range[1], inf$range[2]))
  cat(sprintf("  Magnitude implies    : %s\n",
              if (length(inf$candidates)) paste(inf$candidates, collapse = " or ") else "<no plausible unit>"))
  fail <- function(msg) stop(sprintf("Unit check failed for %s: %s", label, msg), call. = FALSE)
  if (!length(inf$candidates))
    fail(paste0("no candidate unit gives a plausible annual total. Check the variable, ",
                "or adjust UNIT_PLAUSIBLE_ANNUAL_MM in DROUGHT_ANALYSIS_utils.R."))
  if (!is.na(meta)) {
    if (!meta %in% inf$candidates)
      fail(sprintf("metadata says %s but the magnitude implies %s.",
                   meta, paste(inf$candidates, collapse = " or ")))
    unit <- meta
    src  <- if (length(inf$candidates) == 1L) "metadata and magnitude agree" else
      "metadata (magnitude consistent but not unique)"
  } else {
    if (is.na(inf$unit))
      fail(sprintf("metadata units are missing/unrecognised and the magnitude is ambiguous (%s).",
                   paste(inf$candidates, collapse = " or ")))
    unit <- inf$unit
    src  <- "inferred from magnitude (metadata missing or unrecognised)"
  }
  cat(sprintf("  Unit used            : %s  [%s]\n", unit, src))
  invisible(list(unit = unit, source = src, metadata_raw = meta_raw, metadata_unit = meta,
                 inferred_unit = inf$unit, mean_value = inf$mean_value,
                 annual_totals = inf$annual_totals))
}

days_in_month_vec <- function(dates) {
  f <- as.Date(format(as.Date(dates), "%Y-%m-01"))
  nxt <- as.Date(vapply(f, function(d) seq(d, by = "month", length.out = 2L)[2L], numeric(1)),
                 origin = "1970-01-01")
  as.integer(nxt - f)
}

#' Convert a monthly depth raster between units; `dates` give calendar days.
convert_monthly_depth <- function(r, from, to = c("mm/month", "mm/day"), dates) {
  to <- match.arg(to)
  if (!from %in% MONTHLY_DEPTH_UNITS) stop("Unsupported unit: ", from)
  if (terra::nlyr(r) != length(dates)) stop("convert_monthly_depth(): dates/layers mismatch.")
  days <- days_in_month_vec(dates)
  if (to == "mm/month") {
    switch(from, `mm/month` = r, `m/month` = r * 1000,
           `mm/day` = r * days, `m/day` = r * 1000 * days)
  } else {
    switch(from, `mm/day` = r, `m/day` = r * 1000,
           `mm/month` = r / days, `m/month` = r * 1000 / days)
  }
}


## M. DISTRIBUTION SELECTION BY L-MOMENT-RATIO DISTANCE ########################
## Used by SPEI (Script 3) and SSI (Script 6).
## All candidates are 3-parameter distributions fitted by L-moments (lmomco).
## Such a fit reproduces the sample lambda1, lambda2 AND tau3 exactly, so tau3
## cannot discriminate between candidates; the criterion is the L-kurtosis
## distance |tau4_sample - tau4_fitted| (L-moment ratio diagram; Hosking &
## Wallis 1997). tau3 is still checked: a candidate whose fit does not
## reproduce the sample tau3 (outside the family's attainable range) is rejected.
LMOM_TAU3_TOLERANCE <- 1e-3

#' Fit one candidate and return list(distance, par); distance = Inf on failure.
lmom_ratio_distance <- function(lmom_sample, dist, tau3_tol = LMOM_TAU3_TOLERANCE) {
  bad <- list(distance = Inf, par = NULL)
  par <- tryCatch(lmomco::lmom2par(lmom_sample, type = dist), error = function(e) NULL)
  if (is.null(par) || any(!is.finite(par$para)) ||
      (!is.null(par$ifail) && any(par$ifail != 0))) return(bad)
  th <- tryCatch(lmomco::par2lmom(par), error = function(e) NULL)
  if (is.null(th)) return(bad)
  t3s <- lmom_sample$ratios[3]; t4s <- lmom_sample$ratios[4]
  t3t <- th$ratios[3];          t4t <- th$ratios[4]
  if (!all(is.finite(c(t3s, t4s, t3t, t4t)))) return(bad)
  if (abs(t3s - t3t) > tau3_tol) return(bad)
  list(distance = abs(t4s - t4t), par = par)
}

#' Rank candidates (lmomco type codes, e.g. c("glo","pe3","gev")) for `x`.
#' Returns a data.frame ordered by distance (ties keep the candidate order).
lmom_rank_distributions <- function(x, candidates, tau3_tol = LMOM_TAU3_TOLERANCE) {
  x <- x[is.finite(x)]
  lm <- tryCatch(lmomco::lmoms(x, nmom = 4), error = function(e) NULL)
  if (is.null(lm) || !isTRUE(lmomco::are.lmom.valid(lm)))
    return(data.frame(dist = candidates, distance = Inf, stringsAsFactors = FALSE))
  fits <- lapply(candidates, function(d) lmom_ratio_distance(lm, d, tau3_tol))
  out <- data.frame(dist = candidates,
                    distance = vapply(fits, `[[`, numeric(1), "distance"),
                    stringsAsFactors = FALSE)
  ord <- order(out$distance, seq_along(candidates))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "fits") <- stats::setNames(lapply(fits[ord], `[[`, "par"), candidates[ord])
  out
}

## F2. SPATIAL HELPERS (used by 13trends_visualization.R) ----------------------------------------
#' Load the authoritative Nechako basin from the supplied subbasin
#' Shapefiles and transform the dissolved union to the requested CRS.
load_basin <- function(shp_path = BASIN_SHP, crs = EQUAL_AREA_CRS) {
  if (!dir.exists(shp_path))
    stop("Nechako subbasin Shapefile directory not found: ", shp_path)
  b <- sf::st_as_sf(terra::project(load_basin_vect(shp_path), crs))
  cat(sprintf("Nechako basin loaded from %d Shapefiles (dissolved, CRS: %s)\n",
              NECHAKO_SUBBASIN_EXPECTED_N, crs))
  b
}

# Native-grid rasterization of producer-selected coordinate rows. No resampling,
# polygon centre test, or reprojection; checked duplicate cells prevent overwrite.
analysis_rasterize_table <- function(d, template, value_col,
                                     x_col = "lon", y_col = "lat") {
  if (!all(c(x_col, y_col, value_col) %in% names(d)))
    stop("Missing native rasterization field: ", value_col)
  r <- template[[1]]
  terra::values(r) <- NA_real_
  xy <- as.matrix(as.data.frame(d)[, c(x_col, y_col), drop = FALSE])
  cell <- terra::cellFromXY(r, xy)
  v <- as.numeric(d[[value_col]])
  ok <- !is.na(cell) & cell >= 1L & cell <= terra::ncell(r) & is.finite(v)
  if (!any(ok)) return(r)
  cell <- cell[ok]
  if (anyDuplicated(cell)) stop("Duplicate native-grid cells in ", value_col)
  vals <- terra::values(r, mat = FALSE)
  vals[cell] <- v[ok]
  terra::values(r) <- vals
  r
}

#' Optional 3x3 focal mean; do not fill NA centre cells outside the tested domain.
smooth_raster <- function(r) {
  if (is.null(r) || all(is.na(terra::values(r)))) return(r)
  terra::focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE,
               na.policy = "omit")
}

#' Plot native per-pixel statistics without smoothing unless explicitly requested.
#' Relies on basin_sf_global being set in the calling script.
plot_raster_clean <- function(r, main, zlim = NULL, col, breaks = NULL,
                              legend = TRUE, categorical = FALSE, legend_title = "",
                              smooth = FALSE) {
  if (is.null(r) || all(is.na(terra::values(r)))) {
    plot.new(); text(0.5, 0.5, "No valid data", cex = 1.2, col = "red")
    return(invisible())
  }
  if (isTRUE(smooth) && !categorical) r <- smooth_raster(r)
  plg <- list(cex = 0.9)
  if (nchar(legend_title)) plg$title <- legend_title
  if (!is.null(breaks)) {
    terra::plot(r, main = main, cex.main = 1.0, col = col, breaks = breaks,
                axes = FALSE, box = FALSE, legend = legend, plg = plg, colNA = NA)
  } else {
    terra::plot(r, main = main, cex.main = 1.0, col = col, zlim = zlim,
                axes = FALSE, box = FALSE, legend = legend, plg = plg, colNA = NA)
  }
  basin_overlay <- sf::st_transform(basin_sf_global, crs = terra::crs(r))
  plot(sf::st_geometry(basin_overlay),
       add = TRUE, col = NA, border = "black", lwd = 2.5)
}

## G. GGPLOT HELPERS ----------------------------------------
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

index_colours <- c(spi = "#e41a1c", spei = "#377eb8", spei_era5 = "#d62728",
                   swei = "#1E90FF", sspi = "#9467bd", sri = "#a65628",
                   sdi = "#f781bf", ssmi = "#33a02c", ssi_sm = "#b15928")

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
    if (isTRUE(capabilities("cairo")))
      grDevices::cairo_pdf(filename, width = width, height = height)
    else
      grDevices::pdf(filename, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
    TRUE
  }, error = function(e) {
    cat("ERROR: Cannot open PDF '", filename, "': ", e$message, "\n")
    FALSE
  })
}

#' Unified figure-saving helper --- replaces the scattered mix of pdf()/dev.off()
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
                    width = width, height = height, units = "in",
                    device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf),
    error = function(e) cat(sprintf("  [WARN] PDF save failed (%s): %s\n",
                                    basename(pdf_path), e$message))
  )
  tryCatch(
    ggplot2::ggsave(png_path, plot_obj,
                    width = width, height = height, units = "in",
                    dpi = dpi, device = "png"),
    error = function(e) cat(sprintf("  [WARN] PNG save failed (%s): %s\n",
                                    basename(png_path), e$message))
  )
  cat(sprintf("  [OK] Figure saved: %s (.pdf + .png)\n", basename(stem)))
  invisible(NULL)
}

# J. SHARED STATISTICAL FUNCTIONS

.analysis_homog_try <- function(expr) {
  tryCatch(suppressWarnings(expr), error = function(e) NA_real_)
}

# Common one-row battery. Retains all Script 11 columns and adds the aliases
# expected by Script 13, without changing fixed-step or AR(1)-GLS inference.
analysis_homogeneity_battery <- function(x, years, step_year,
                                         m_sim = ANALYSIS_HOMOG_MC,
                                         min_unique = 5L,
                                         run_sim_tests = TRUE) {
  h <- analysis_homogeneity_tests(x, years, step_year = step_year,
                                  m_sim = m_sim, run_sim_tests = run_sim_tests,
                                  min_unique = min_unique)
  m <- analysis_step_trend_model(x, years, step_year = step_year)
  out <- data.table::as.data.table(cbind(h, m))
  out[, `:=`(n_years = as.integer(n), first_year = as.integer(year_min),
             last_year = as.integer(year_max), diff_at_step = diff_after_minus_before)]
  out
}

# Basin summaries from either pipeline have the same ordered set of columns.
# FDR applies only to the rows passed by the caller (separate hypothesis families);
# the legacy pettitt_p_fdr_basin_all column is retained as a compatibility alias.
analysis_homogeneity_basin_summary <- function(x) {
  d <- data.table::as.data.table(x)
  if (!"pettitt_p_bh_across_series" %in% names(d)) {
    if ("pettitt_p_fdr_basin_all" %in% names(d)) {
      d[, pettitt_p_bh_across_series := pettitt_p_fdr_basin_all]
    } else {
      d[, pettitt_p_bh_across_series := stats::p.adjust(pettitt_p, "BH")]
    }
  }
  if (!"pettitt_p_fdr_basin_all" %in% names(d))
    d[, pettitt_p_fdr_basin_all := pettitt_p_bh_across_series]
  if (!"step_p_bh_across_series" %in% names(d))
    d[, step_p_bh_across_series := stats::p.adjust(step_p, method = "BH")]
  columns <- c(
    "region_id", "region_name", "series", "index_type", "scale", "season", "family", "kind",
    "n", "n_years", "year_min", "year_max", "first_year", "last_year",
    "pettitt_p", "pettitt_last_year_before", "pettitt_first_year_after",
    "pettitt_p_fdr", "pettitt_p_fdr_basin_all", "pettitt_p_bh_across_series",
    "snht_p", "snht_first_year_after", "buishand_range_p",
    "n_before", "n_after", "mean_before", "mean_after",
    "diff_after_minus_before", "diff_at_step", "pct_change", "cohens_d",
    "welch_p", "mannwhitney_p", "variance_fligner_p",
    "trend_no_step", "trend_no_step_p", "trend_with_step", "trend_with_step_p",
    "step_coef", "step_p", "step_p_fdr", "step_p_bh_across_series",
    "delta_AIC_step_vs_nostep", "ar1_phi_with_step", "break_near_step_year")
  for (nm in setdiff(columns, names(d))) d[, (nm) := NA]
  d[, ..columns]
}

analysis_homogeneity_tests <- function(x, yrs, step_year, m_sim = ANALYSIS_HOMOG_MC,
                                       run_sim_tests = TRUE, min_unique = 5L) {
  clean <- analysis_homog_clean(x, yrs)
  x <- clean$x; yrs <- clean$years; n <- length(x)
  o <- list(n = n, year_min = if (n) min(yrs) else NA_integer_,
            year_max = if (n) max(yrs) else NA_integer_,
            pettitt_p = NA_real_, pettitt_last_year_before = NA_integer_,
            pettitt_first_year_after = NA_integer_,
            snht_p = NA_real_, snht_first_year_after = NA_integer_,
            buishand_range_p = NA_real_,
            n_before = NA_integer_, n_after = NA_integer_,
            mean_before = NA_real_, mean_after = NA_real_, diff_after_minus_before = NA_real_,
            pct_change = NA_real_, cohens_d = NA_real_,
            welch_p = NA_real_, mannwhitney_p = NA_real_, variance_fligner_p = NA_real_)
  if (n < ANALYSIS_MIN_HOMOG_YEARS || length(unique(x)) < min_unique) return(data.table::as.data.table(o))
  
  pt <- tryCatch(analysis_pettitt(x), error = function(e) NULL)
  if (!is.null(pt)) {
    K <- as.integer(pt$estimate)[1]
    o$pettitt_p <- pt$p.value
    o$pettitt_last_year_before  <- yrs[K]
    o$pettitt_first_year_after  <- yrs[min(K + 1L, n)]
  }
  if (run_sim_tests) {
    sn <- tryCatch(trend::snh.test(x, m = m_sim), error = function(e) NULL)
    if (!is.null(sn)) { o$snht_p <- analysis_mc_p(sn$p.value, m_sim)
    o$snht_first_year_after <- yrs[min(as.integer(sn$estimate)[1] + 1L, n)] }
    br <- tryCatch(trend::br.test(x, m = m_sim), error = function(e) NULL)
    if (!is.null(br)) o$buishand_range_p <- analysis_mc_p(br$p.value, m_sim)
  }
  b <- x[yrs <  step_year]; a <- x[yrs >= step_year]
  o$n_before <- length(b); o$n_after <- length(a)
  if (length(b) >= 5L && length(a) >= 5L) {
    o$mean_before <- mean(b); o$mean_after <- mean(a)
    o$diff_after_minus_before <- o$mean_after - o$mean_before
    o$pct_change <- if (isTRUE(abs(o$mean_before) > 0)) 100 * o$diff_after_minus_before / abs(o$mean_before) else NA_real_
    sp <- sqrt(((length(b) - 1) * stats::var(b) + (length(a) - 1) * stats::var(a)) /
                 (length(a) + length(b) - 2))
    o$cohens_d <- if (is.finite(sp) && sp > 0) o$diff_after_minus_before / sp else NA_real_
    o$welch_p        <- .analysis_homog_try(stats::t.test(a, b)$p.value)
    o$mannwhitney_p  <- .analysis_homog_try(stats::wilcox.test(a, b, exact = FALSE)$p.value)
    o$variance_fligner_p <- .analysis_homog_try(stats::fligner.test(list(b, a))$p.value)
  }
  data.table::as.data.table(o)
}

# Shared H1/H2/H3 homogeneity panels for both climate-variable and drought-index
# pipelines. The caller supplies its own homogeneity directory and saves the
# returned plots under its existing output filenames. H4/H5/H6 remain
# pipeline-specific because their native grids and annual diagnostics differ.
analysis_homogeneity_plots <- function(homog_dir, step_year, basin_id,
                                       series_levels = NULL, break_window = 3L) {
  rd <- function(file) {
    path <- file.path(homog_dir, file)
    if (file.exists(path)) data.table::fread(path) else NULL
  }
  basin <- rd("homogeneity_break_summary_basin.csv")
  regions <- rd("step_tests_by_region_all.csv")
  pixels <- rd("pixel_step_summary.csv")
  if (is.null(basin) || !nrow(basin))
    return(list(basin = basin, regions = regions, pixels = pixels,
                plots = list(basin = NULL, regions = NULL, pixels = NULL)))
  classes <- c(sprintf("Step at %d favoured (BH; exploratory)", step_year),
               "Pettitt break only (BH)", "No BH-supported break")
  colours <- stats::setNames(c("#b2182b", "#ef8a62", "grey60"), classes)
  classify <- function(dAIC, step_p_bh, pettitt_p_bh) {
    data.table::fifelse(is.finite(dAIC) & dAIC < -2 &
                          is.finite(step_p_bh) & step_p_bh < 0.05, classes[1],
                        data.table::fifelse(is.finite(pettitt_p_bh) & pettitt_p_bh < 0.05,
                                            classes[2], classes[3]))
  }
  if (is.null(series_levels)) series_levels <- unique(basin$series)
  series_levels <- c(series_levels, setdiff(unique(basin$series), series_levels))
  b <- data.table::copy(basin)
  if (!"family" %in% names(b)) b[, family := "index"]
  if (!"kind" %in% names(b)) b[, kind := "raw"]
  if (!"pettitt_p_bh_across_series" %in% names(b))
    b[, pettitt_p_bh_across_series := stats::p.adjust(pettitt_p, method = "BH")]
  if (!"step_p_bh_across_series" %in% names(b))
    b[, step_p_bh_across_series := stats::p.adjust(step_p, method = "BH")]
  b[, `:=`(class = factor(classify(delta_AIC_step_vs_nostep, step_p_bh_across_series,
                                   pettitt_p_bh_across_series), levels = classes),
           series = factor(series, levels = rev(series_levels)))]
  yr0 <- if (all(c("first_year", "last_year") %in% names(b)))
    c(b$first_year, b$last_year) else c(b$year_min, b$year_max)
  yr0 <- yr0[is.finite(yr0)]
  if (!length(yr0)) yr0 <- c(step_year - 30, step_year + 30)
  yr <- range(yr0)
  b[, lab := sprintf("Pettitt BH %.3g | step BH %.3g | dAIC %s",
                     pettitt_p_bh_across_series, step_p_bh_across_series,
                     ifelse(is.finite(delta_AIC_step_vs_nostep),
                            sprintf("%+.1f", delta_AIC_step_vs_nostep), "NA"))]
  p1 <- ggplot2::ggplot(b, ggplot2::aes(pettitt_first_year_after, series)) +
    ggplot2::annotate("rect", xmin = step_year - break_window,
                      xmax = step_year + break_window, ymin = -Inf, ymax = Inf,
                      fill = "grey85", alpha = 0.5) +
    ggplot2::geom_vline(xintercept = step_year, linetype = "dashed", colour = "grey35") +
    ggplot2::geom_vline(xintercept = mean(yr), linetype = "dotted", colour = "grey35") +
    ggplot2::geom_point(ggplot2::aes(colour = class, size = -log10(pmax(pettitt_p, 1e-10))), na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(x = yr[2] + 2, label = lab), hjust = 0, size = 2.3, colour = "grey30") +
    ggplot2::facet_grid(family ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE, name = NULL) +
    ggplot2::scale_size_continuous(range = c(1.2, 4.5), name = "-log10 nominal Pettitt p") +
    ggplot2::scale_x_continuous(limits = c(yr[1], yr[2] + 36)) +
    ggplot2::labs(title = "Break year and step-vs-trend classification — basin mean",
                  subtitle = sprintf("Dashed = exploratory fixed %d; step needs dAIC < -2 and BH-adjusted step p < 0.05 across basin series.", step_year),
                  x = "First year of the later regime", y = NULL) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom", strip.text.y = ggplot2::element_text(angle = 0))
  p2 <- NULL
  r <- regions
  if (!is.null(r) && nrow(r)) {
    r <- data.table::copy(r)
    # Regional sensitivity: BH across regions within each series, NOT basin-series BH.
    if (!"step_p_fdr" %in% names(r))
      r[, step_p_fdr := stats::p.adjust(step_p, method = "BH"), by = series]
    if (!"pettitt_p_fdr" %in% names(r))
      r[, pettitt_p_fdr := stats::p.adjust(pettitt_p, method = "BH"), by = series]
    r[, `:=`(class = factor(classify(delta_AIC_step_vs_nostep, step_p_fdr,
                                     pettitt_p_fdr), levels = classes),
             region = factor(region_id, levels = c(intersect(basin_id, unique(region_id)),
                                                   setdiff(sort(unique(region_id)), basin_id))),
             series = factor(series, levels = rev(series_levels)))]
    r[, lab := ifelse(class == classes[3] | !is.finite(pettitt_first_year_after),
                      "", as.character(pettitt_first_year_after))]
    p2 <- ggplot2::ggplot(r, ggplot2::aes(region, series, fill = class)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(ggplot2::aes(label = lab), size = 2.3) +
      ggplot2::scale_fill_manual(values = colours, drop = FALSE, name = NULL) +
      ggplot2::theme_minimal(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                     legend.position = "bottom") +
      ggplot2::labs(title = "Step-vs-trend classification by region",
                    subtitle = "BH adjustment within each series across regions; exploratory fixed-date step.",
                    x = NULL, y = NULL)
  }
  p3 <- NULL
  p <- pixels
  if (!is.null(p) && nrow(p)) {
    p <- data.table::copy(p)[order(pct_pettitt_fdr05)]
    p[, series := factor(series, levels = unique(series))]
    if (!"pct_welch_fdr05" %in% names(p)) p[, pct_welch_fdr05 := NA_real_]
    p[, lab := ifelse(is.finite(median_break_year_fdr_sig),
                      sprintf("median break %.0f", median_break_year_fdr_sig), "")]
    p3 <- ggplot2::ggplot(p, ggplot2::aes(pct_pettitt_fdr05, series)) +
      ggplot2::geom_col(fill = "#4575b4", width = 0.7, na.rm = TRUE) +
      ggplot2::geom_text(ggplot2::aes(label = lab), hjust = -0.05, size = 2.5) +
      ggplot2::scale_x_continuous(limits = c(0, 130), breaks = seq(0, 100, 25)) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(title = "Share of basin pixels with a Pettitt break (BH-FDR < 0.05)",
                    subtitle = "A midpoint break can reflect a sustained trend rather than an abrupt step.",
                    x = "% of tested pixels", y = NULL)
  }
  list(basin = b, regions = r, pixels = p,
       plots = list(basin = p1, regions = p2, pixels = p3))
}

# Rank-based Pettitt statistic and its asymptotic p-value (single implementation).
# The estimated break is AFTER observation tau; U has length n-1. Non-finite
# inputs are rejected to avoid silently changing the year-to-value alignment.
analysis_pettitt <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3L || any(!is.finite(x)))
    return(list(U = numeric(0), K = NA_real_, tau = NA_integer_, n = n,
                p = NA_real_, p.value = NA_real_, estimate = NA_integer_))
  tt <- seq_len(n - 1L)
  U <- 2 * cumsum(rank(x, ties.method = "average"))[tt] - tt * (n + 1L)
  K <- max(abs(U))
  tau <- which.max(abs(U))
  p <- min(1, 2 * exp(-6 * K^2 / (n^3 + n^2)))
  list(U = U, K = K, tau = tau, n = n, p = p, p.value = p,
       estimate = as.integer(tau))
}

# ------------------------------------------------------------------------------
# Trend with and without a step term, AR(1) errors (nlme::gls, ML so AIC is comparable)
# ------------------------------------------------------------------------------
analysis_step_trend_model <- function(x, yrs, step_year) {
  d <- data.frame(y = as.numeric(x), t = as.integer(yrs) - min(as.integer(yrs), na.rm = TRUE),
                  step = as.integer(yrs >= step_year))
  d <- d[is.finite(d$y), , drop = FALSE]
  na <- data.table::data.table(trend_no_step = NA_real_, trend_no_step_p = NA_real_,
                               trend_with_step = NA_real_, trend_with_step_p = NA_real_,
                               step_coef = NA_real_, step_p = NA_real_, delta_AIC_step_vs_nostep = NA_real_,
                               ar1_phi_with_step = NA_real_)
  if (nrow(d) < ANALYSIS_MIN_HOMOG_YEARS || length(unique(d$step)) < 2L) return(na)
  fit <- function(f) tryCatch(
    suppressWarnings(nlme::gls(f, data = d, correlation = nlme::corAR1(form = ~ t),
                               method = "ML")), error = function(e) NULL)
  m0 <- fit(y ~ t); m1 <- fit(y ~ t + step)
  if (is.null(m0) || is.null(m1)) return(na)
  tt0 <- summary(m0)$tTable; tt1 <- summary(m1)$tTable
  data.table::data.table(trend_no_step = tt0["t", "Value"], trend_no_step_p = tt0["t", "p-value"],
                         trend_with_step = tt1["t", "Value"], trend_with_step_p = tt1["t", "p-value"],
                         step_coef = tt1["step", "Value"], step_p = tt1["step", "p-value"],
                         delta_AIC_step_vs_nostep = AIC(m1) - AIC(m0),   # negative = step model preferred
                         ar1_phi_with_step = tryCatch(as.numeric(coef(m1$modelStruct$corStruct,
                                                                      unconstrained = FALSE)),
                                                      error = function(e) NA_real_))
}


# ----------------------------------------
# These were previously duplicated between 13trend_test_ALL.R (vectorised pixel
# versions) and 11pr_pet_trends.r (single-series version).  They are now defined
# once here and sourced by both scripts.
#
# The public API is intentionally backward-compatible:
#   --- mk_tfpw_spectral_for_series() --- identical signature to the function that
#     was previously embedded in 11pr_pet_trends.r.  That script calls it
#     unchanged.
#   --- vectorized_mann_kendall(), vectorized_mann_kendall_tfpw(),
#     vectorized_regime_shift_pelt(), vectorized_spectral_peaks() --- identical
#     signatures to the functions previously in 13trend_test_ALL.R.  That script
#     now omits its local copies and relies on these definitions.
#   --- detect_regime_shift_pelt() --- single-series PELT helper, previously in
#     11pr_pet_trends.r under that name.  "trend_test_ALL.R" used a different but equivalent
#     vectorised wrapper; both now share this single-series core.
####################################################################################

# --- J0. Utility operator ----------------------------------------
`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.atomic(a)) {
    if (length(a) > 0 && !is.na(a[1])) return(a) else return(b)
  }
  # Non-atomic objects (SpatRaster, PackedSpatRaster, list, S4, data.frame, etc.)
  # can't be NA-checked via a[1] - existence (non-NULL) is enough.
  a
}

# ----------------------------------------
# J1. PRIMITIVE HELPERS
# ----------------------------------------

#' Check whether too many values in a series are at or below a floor threshold.
#' Returns a list(exceeds_threshold, pct_min_vals).
#' Used by mk_tfpw_spectral_for_series() (formerly embedded in 11pr_pet_trends.r).
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
#' Used by mk_tfpw_for_series() (formerly embedded in 11pr_pet_trends.r).
calculate_variance_with_ties <- function(S, n, x) {
  tie_table  <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  var_s      <- n * (n - 1) * (2 * n + 5) / 18
  if (length(tie_counts))
    var_s <- var_s - sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
  var_s
}

# ----------------------------------------
# J2. PELT CHANGEPOINT --- single series
# ----------------------------------------
#' Detect the first mean-shift regime in a univariate time series using PELT/BIC.
#' Returns a named list with fields:
#'   changepoint_detected, changepoint_position, n_changepoints,
#'   first_changepoint_year (first year AFTER the break), mean_before,
#'   mean_after, magnitude_shift. The original field name is retained for CSV compatibility.
#'
#' Previously defined as detect_regime_shift_pelt() in 11pr_pet_trends.r;
#' 13trend_test_ALL.R used an inline equivalent inside vectorized_regime_shift_pelt().
#' Both now share this function.
detect_regime_shift_pelt <- function(ts_vec, min_obs = 20, start_year = NULL,
                                     minseglen = 5L, standardize = TRUE) {
  # changepoint::cpt.mean() with the default Normal test statistic ASSUMES UNIT
  # VARIANCE. On series whose sd is >> 1 (annual precipitation ~80 mm, PET_ERA5
  # ~300 mm) the BIC penalty is then negligible relative to the cost, and
  # spurious breaks are placed at the series edge (e.g., first_changepoint_year 1951/52).
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
    fc <- as.integer(cpts_det[1L])  # index of the last observation BEFORE the break
    if (!is.finite(fc) || fc < 1L || fc >= n) return(empty)
    list(changepoint_detected   = TRUE,
         changepoint_position   = fc,
         n_changepoints         = length(cpts_det),
         first_changepoint_year = as.integer(
           if (length(start_year) == 1L && is.finite(start_year)) start_year + fc else fc + 1L),
         mean_before            = mean(ts_clean[seq_len(fc)], na.rm = TRUE),
         mean_after             = mean(ts_clean[(fc + 1L):n],  na.rm = TRUE),
         magnitude_shift        = mean(ts_clean[(fc + 1L):n],  na.rm = TRUE) -
           mean(ts_clean[seq_len(fc)],   na.rm = TRUE))
  }, error = function(e) empty)
}

# ----------------------------------------
# J3. VC-MK --- single series (core used by both the vectorised wrapper and
#     mk_tfpw_spectral_for_series)
# ----------------------------------------
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

## J3b. Hamed-Rao p-value with the framework guard rails (all VC-MK paths).
##   * Input names are stripped: mmkh() builds its result with c(), so a named
##     input turns "new P-value" into "new P-value.<name>" and exact lookups fail.
##   * n/N* < 1 (Hamed-Rao would SHRINK Var(S)): uncorrected MK p ("capped").
##   * n/N* < 1 is capped even when the corrected variance is nonpositive.
##   * All other Var(S) <= 0 or non-finite p: NA by default ("invalid_variance").
##     The uncorrected p is NOT a safe fallback for accumulated (k >= 3) series,
##     whose strong positive autocorrelation makes it far too small; use
##     invalid = "uncorrected" only for explicitly documented sensitivity runs.
vc_mmkh <- function(x, invalid = c("na", "uncorrected"), cap_shrinkage = TRUE) {
  invalid <- match.arg(invalid)
  x <- as.numeric(x); x <- x[is.finite(x)]
  out <- list(p = NA_real_, slope = NA_real_, tau = NA_real_, n_over_nstar = NA_real_,
              p_uncorrected = NA_real_, varS = NA_real_, varS_old = NA_real_,
              status = "filtered")
  if (length(x) < ANALYSIS_MIN_TREND_N || stats::var(x) < 1e-6) return(out)
  hr <- tryCatch(suppressWarnings(modifiedmk::mmkh(x)), error = function(e) NULL)
  if (is.null(hr)) { out$status <- "mmkh_error"; return(out) }
  mm <- .mmkh_fields(hr)
  sen_col <- grep("^Sen", names(hr))[1L]
  if (length(sen_col) && is.finite(sen_col))
    out$slope <- suppressWarnings(as.numeric(hr[[sen_col]]))
  out$tau <- mm$tau; out$n_over_nstar <- mm$n_over_nstar
  out$p_uncorrected <- mm$p_old
  if (!is.finite(out$p_uncorrected))
    out$p_uncorrected <- tryCatch(as.numeric(Kendall::MannKendall(x)$sl),
                                  error = function(e) NA_real_)
  if (!is.finite(out$tau))
    out$tau <- tryCatch(as.numeric(Kendall::MannKendall(x)$tau),
                        error = function(e) NA_real_)
  out$varS <- mm$var_new; out$varS_old <- mm$var_old
  valid <- is.finite(mm$var_new) && mm$var_new > 0 && is.finite(mm$p_new)
  shrink <- is.finite(mm$n_over_nstar) && mm$n_over_nstar < 1
  # Cap all finite n/n* < 1 cases before testing the corrected variance:
  # extreme negative serial correlation can produce var_new <= 0 and NaN p_new.
  if (isTRUE(cap_shrinkage) && shrink && is.finite(out$p_uncorrected)) {
    out$p <- out$p_uncorrected; out$varS <- mm$var_old
    out$status <- "capped_uncorrected"
  } else if (valid) {
    out$p <- mm$p_new; out$status <- "corrected"
  } else if (invalid == "uncorrected") {
    out$p <- out$p_uncorrected; out$varS <- mm$var_old
    out$status <- "invalid_variance_uncorrected"
  } else {
    out$status <- "invalid_variance"
  }
  out
}
# Canonical shared Yue et al. (2002) TFPW definition for J4 and J6:
# retain the first observation; prewhiten detrended residuals only if |rho1| > 0.1;
# retrend; compute conventional MK S and tie-adjusted variance on the result.
.analysis_tfpw_transform <- function(x, slope, rho1, retain_first = TRUE,
                                     rho_threshold = 0.1) {
  n <- length(x)
  t <- seq_len(n)
  detrended <- x - slope * t
  apply_pw <- is.finite(rho1) && abs(rho1) > rho_threshold
  if (!isTRUE(retain_first)) {
    r <- if (apply_pw) rho1 else 0
    return(detrended[-1L] - r * detrended[-n] + slope * t[-1L])
  }
  if (!apply_pw) return(x)
  c(detrended[1L], detrended[-1L] - rho1 * detrended[-n]) + slope * t
}

# Kendall S with tie-aware Fenwick prefix counts: O(n log n) time, O(n) space.
# Avoids allocating an n x n difference matrix for full-monthly (n ~ 900) pixels.
.analysis_kendall_S <- function(x) {
  n <- length(x)
  ranks <- match(x, sort(unique(x)))
  bit <- integer(max(ranks))
  prefix <- function(k) {
    total <- 0L
    while (k > 0L) {
      total <- total + bit[k]
      k <- k - bitwAnd(k, -k)
    }
    total
  }
  S <- 0
  for (i in seq_len(n)) {
    r <- ranks[i]
    less <- prefix(r - 1L)
    equal <- prefix(r) - less
    S <- S + 2 * less + equal - (i - 1L)
    while (r <= length(bit)) {
      bit[r] <- bit[r] + 1L
      r <- r + bitwAnd(r, -r)
    }
  }
  S
}

.analysis_tfpw_mk <- function(x, slope, rho1) {
  corrected <- .analysis_tfpw_transform(x, slope, rho1)
  n <- length(corrected)
  S <- .analysis_kendall_S(corrected)
  varS <- calculate_variance_with_ties(S, n, corrected)
  list(tau = S / (n * (n - 1) / 2),
       p = if (is.finite(varS) && varS > 0) 2 * stats::pnorm(-abs(S / sqrt(varS))) else NA_real_,
       S = S, varS = varS,
       applied = is.finite(rho1) && abs(rho1) > 0.1)
}

# ----------------------------------------
# J4. TFPW-MK --- single series (same inferential estimator as J6)
# ----------------------------------------
#' Original-series Sen slope and shared TFPW MK inference; tau is tau-a.
.mk_tfpw_single <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < ANALYSIS_MIN_TREND_N || !is.finite(stats::var(x)) || stats::var(x) < 1e-6)
    return(list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
                rho1 = NA_real_, tfpw_applied = FALSE, filtered = TRUE))
  tryCatch({
    beta <- as.numeric(trend::sens.slope(x)$estimates)
    if (!is.finite(beta)) stop("Invalid Sen slope")
    xd <- x - beta * seq_len(n)
    ac <- tryCatch(stats::acf(xd, lag.max = 1L, plot = FALSE)$acf[2L],
                   error = function(e) NA_real_)
    rho1 <- if (is.finite(ac)) as.numeric(ac) else NA_real_
    mk <- .analysis_tfpw_mk(x, beta, rho1)
    list(tau = mk$tau, pval = mk$p, slope = beta, rho1 = rho1,
         tfpw_applied = mk$applied, filtered = !is.finite(mk$p))
  }, error = function(e)
    list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
         rho1 = NA_real_, tfpw_applied = FALSE, filtered = TRUE))
}

# ----------------------------------------
# J5. SPECTRAL PEAK DETECTION --- single series
# ----------------------------------------
#' Monte Carlo AR(1) red-noise spectral test; same implementation for J5 and J7d.
#' The observed and simulated spectra use identical detrending, taper and
#' smoothing. A max-over-frequency threshold controls family-wise Type I error
#' within periods >= min_period_months; contiguous significant bins count once.
#' dominant_period is returned in YEARS for either monthly or annual sampling.
#' conf is the variance-normalized maximum-spectral-power threshold.
.spectral_single <- function(x,
                             min_period_months = 24L,
                             sig_level = 0.95,
                             n_sim = 500L,
                             alpha = 0.05,
                             conf_cache_env = NULL,
                             sample_interval_months = 1L) {
  if (length(sample_interval_months) != 1L || !is.finite(sample_interval_months) ||
      sample_interval_months <= 0) stop("sample_interval_months must be positive")
  if (length(n_sim) != 1L || !is.finite(n_sim) || n_sim < 1L)
    stop("n_sim must be a positive integer")
  empty <- list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_)
  x <- as.numeric(x)
  finite <- which(is.finite(x))
  if (!length(finite)) return(empty)
  x <- x[min(finite):max(finite)]
  # Never collapse internal gaps: that would change the sampling interval.
  if (any(!is.finite(x)) || length(x) < 48L) return(empty)
  n <- length(x)
  x_dt <- stats::residuals(stats::lm(x ~ seq_len(n)))
  vv <- stats::var(x_dt)
  if (!is.finite(vv) || vv <= 0) return(empty)
  x_std <- x_dt / sqrt(vv)
  sp <- tryCatch(stats::spectrum(x_std, spans = c(3L, 5L), taper = 0,
                                 detrend = FALSE, demean = FALSE, plot = FALSE),
                 error = function(e) NULL)
  if (is.null(sp)) return(empty)
  band <- is.finite(sp$freq) & sp$freq > 0 &
    sample_interval_months / sp$freq >= min_period_months
  if (!any(band)) return(empty)
  rho <- suppressWarnings(stats::cor(x_std[-n], x_std[-1L]))
  if (!is.finite(rho)) rho <- 0
  rho_bin <- max(-0.95, min(0.95, round(rho / 0.05) * 0.05))
  env <- if (is.environment(conf_cache_env)) conf_cache_env else SPECTRAL_RN_ENV_CACHE
  key <- paste(n, n_sim, min_period_months, sample_interval_months,
               sprintf("%.3f", sig_level), sprintf("%.2f", rho_bin),
               paste(RNGkind(), collapse = "_"), sep = "_")
  if (!exists(key, envir = env, inherits = FALSE)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(42L)
    maxima <- vapply(seq_len(as.integer(n_sim)), function(i) {
      innov <- stats::rnorm(n + 500L, sd = sqrt(max(1e-8, 1 - rho_bin^2)))
      rn <- utils::tail(as.numeric(stats::filter(innov, filter = rho_bin,
                                                 method = "recursive")), n)
      residual <- stats::residuals(stats::lm(rn ~ seq_len(n)))
      v <- stats::var(residual)
      if (!is.finite(v) || v <= 0) return(NA_real_)
      null_sp <- stats::spectrum(residual / sqrt(v), spans = c(3L, 5L),
                                 taper = 0, detrend = FALSE, demean = FALSE,
                                 plot = FALSE)
      max(null_sp$spec[band], na.rm = TRUE)
    }, numeric(1))
    maxima <- maxima[is.finite(maxima)]
    if (!length(maxima)) return(empty)
    assign(key, as.numeric(stats::quantile(maxima, probs = sig_level,
                                           names = FALSE)), envir = env)
  }
  threshold <- get(key, envir = env, inherits = FALSE)
  above <- band & (sp$spec > threshold)
  n_peaks <- if (any(above)) as.integer(sum(rle(above)$values)) else 0L
  dom <- if (n_peaks) {
    ii <- which(above)
    sample_interval_months / (12 * sp$freq[ii[which.max(sp$spec[ii])]])
  } else NA_real_
  list(n_peaks = n_peaks, dominant_period = as.numeric(dom), conf = threshold)
}

# ----------------------------------------
# J6. OMNIBUS SINGLE-SERIES WRAPPER
#     mk_tfpw_spectral_for_series()
#     Used directly by 11pr_pet_trends.r (and its visualization script reads
#     outputs produced by it).  Signature is identical to the version that
#     was previously embedded in that script.
# ----------------------------------------
#' Full VC-MK + TFPW-MK + spectral + PELT analysis for a single time series.
#'
#' @param ts_vec           Numeric vector (may contain NAs).
#' @param is_precip        Logical; TRUE --- precipitation quality filters.
#' @param alpha            Significance level (default 0.05).
#' @param max_tie_pct      Maximum % ties before filtering (default 50).
#' @param n_sim_spectral   Monte-Carlo replicates for spectral null (default 500).
#' @param sample_interval_months  Months between observations (12 for annual or
#'   fixed-calendar-month series; 1 for continuous monthly series).
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
#'   only INFLATE the MK variance (n/n* >= 1). Shrinkage falls back to the
#'   uncorrected MK p-value even if the corrected variance is nonpositive.
#'   For other invalid variances this wrapper always returns NA (it calls
#'   vc_mmkh(invalid = "na")); sensitivity analyses must call vc_mmkh()
#'   separately with invalid = "uncorrected" and document that choice.
#' @param pelt_minseglen       Minimum PELT segment length (default 5 obs).
#'
#' @return Named list with elements \code{vc}, \code{tf}, \code{spec}, \code{cpt}.
#' Filtered (NA) result in the structure returned by mk_tfpw_spectral_for_series().
.trend_na_result <- function(n, reason, n_ties = 0, pct_ties = 0, pct_min = 0) {
  list(
    vc  = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_,
               varS = NA_real_, n = n, rho1 = NA_real_, vc_corrected = FALSE,
               n_ties = n_ties, percent_ties = pct_ties, n_min_vals = 0L,
               percent_min_vals = pct_min,
               filtered = TRUE, reason = reason, mmkh_status = "not_run"),
    tf  = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_,
               varS = NA_real_, n = n, rho1 = NA_real_, tfpw_applied = FALSE,
               n_ties = n_ties, percent_ties = pct_ties, n_min_vals = 0L,
               percent_min_vals = pct_min,
               filtered = TRUE, reason = reason),
    spec = list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_),
    cpt  = list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                n_changepoints = 0L, first_changepoint_year = NA_integer_,
                mean_before = NA_real_, mean_after = NA_real_,
                magnitude_shift = NA_real_))
}

#' @param max_modal_tie_frac  If finite, also filter series whose largest tied group
#'   exceeds this share (the rule shared with Script 13; ANALYSIS_TIE_MODE_MAX).
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
                                        pelt_minseglen       = 5L,
                                        sample_interval_months = 1L,
                                        max_modal_tie_frac   = NA_real_) {
  ts_clean <- as.numeric(na.omit(ts_vec))
  n        <- length(ts_clean)
  na_result <- function(reason, ...) .trend_na_result(n, reason, ...)
  
  if (n < ANALYSIS_MIN_TREND_N) return(na_result("low_n"))
  
  # --- Quality filters ---------------------------------------------------------
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
  modal_excess <- is.finite(max_modal_tie_frac) &&
    analysis_modal_tie_fraction(ts_clean) > max_modal_tie_frac
  if (percent_ties > max_tie_pct || modal_excess)
    return(na_result("excessive_ties", n_ties, percent_ties, min_check$pct_min_vals))
  
  # --- Sen's slope (original series, shared by both MK variants) ---------------
  sen_slope <- trend::sens.slope(ts_clean)$estimates
  if (is.na(sen_slope) || is.infinite(sen_slope) || length(sen_slope) == 0)
    return(na_result("sens_slope_fail", n_ties, percent_ties, min_check$pct_min_vals))
  
  # --- Detrend -> lag-1 autocorrelation ------------------------------------------
  time_index <- seq_len(n)
  detrended  <- ts_clean - sen_slope * time_index
  acf_result <- tryCatch(
    acf(detrended, lag.max = 1, plot = FALSE, na.action = na.pass),
    error = function(e) NULL, warning = function(w) NULL)
  rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA_real_
  
  # --- VC-MK (Hamed-Rao 1998) via vc_mmkh(): same code path as J3 / J7a ---------
  vc_res <- vc_mmkh(ts_clean,
                    invalid = "na",  # Script 13 contract: invalid VC variance is not significant
                    cap_shrinkage = cap_vc_variance_shrinkage)
  p_vc         <- vc_res$p
  n_over_nstar <- vc_res$n_over_nstar
  vc_failed    <- !is.finite(p_vc)
  vc_capped    <- vc_res$status %in% c("capped_uncorrected", "invalid_variance_uncorrected")
  vc_corrected <- identical(vc_res$status, "corrected") && isTRUE(n_over_nstar > 1)
  vc_reason <- switch(vc_res$status,
                      corrected = "none",
                      capped_uncorrected = "vc_shrinkage_capped_uncorrected",
                      invalid_variance_uncorrected = "vc_invalid_variance_uncorrected",
                      invalid_variance = "vc_invalid_variance",
                      mmkh_error = "mmkh_error",
                      filtered = "vc_filtered",
                      "vc_nonfinite_p")
  vc_status <- if (vc_failed) if (vc_res$status == "mmkh_error") "mmkh_error" else "invalid" else
    if (vc_capped) "capped" else "ok"
  n_min_obs <- as.integer(round(min_check$pct_min_vals / 100 * n))
  vc_list <- list(tau = if (vc_failed) NA_real_ else vc_res$tau,
                  p   = if (vc_failed) NA_real_ else p_vc,
                  sl = sen_slope, S = NA_real_,   # mmkh() does not report raw S
                  varS = vc_res$varS, n = n, rho1 = rho1,
                  vc_corrected = vc_corrected,
                  n_ties = n_ties, percent_ties = percent_ties,
                  n_min_vals = n_min_obs, percent_min_vals = min_check$pct_min_vals,
                  filtered = vc_failed, reason = vc_reason, mmkh_status = vc_status,
                  vc_n_over_nstar = n_over_nstar, vc_p_uncorrected = vc_res$p_uncorrected,
                  vc_capped = vc_capped)
  
  # --- TFPW-MK: identical estimator to J4 -----------------------------------------
  tf_mk <- .analysis_tfpw_mk(ts_clean, sen_slope, rho1)
  tf_list <- list(tau = tf_mk$tau, p = tf_mk$p, sl = sen_slope, S = tf_mk$S,
                  varS = tf_mk$varS, n = n, rho1 = rho1,
                  tfpw_applied = tf_mk$applied,
                  n_ties = n_ties, percent_ties = percent_ties,
                  n_min_vals = n_min_obs, percent_min_vals = min_check$pct_min_vals,
                  filtered = FALSE, reason = "none")
  
  spec_result <- .spectral_single(ts_clean,
                                  min_period_months = 24L,
                                  sig_level         = 0.95,
                                  n_sim             = n_sim_spectral,
                                  alpha             = alpha,
                                  conf_cache_env    = conf_cache_env,
                                  sample_interval_months = sample_interval_months)
  cpt_result <- detect_regime_shift_pelt(ts_clean,
                                         min_obs    = min_obs_changepoint,
                                         start_year = start_year,
                                         minseglen  = pelt_minseglen)
  list(vc = vc_list, tf = tf_list, spec = spec_result, cpt = cpt_result)
}

# ----------------------------------------
# J7. VECTORISED PIXEL-MATRIX WRAPPERS
#     Used directly by 13trend_test_ALL.R.  Signatures are identical to the
#     functions previously defined in that script.
#     The parallel back-end variables N_CORES, is_windows, cl are expected to
#     be set in the calling script (13trend_test_ALL.R) before these are called.
# ----------------------------------------

## J7a. Hamed-Rao (1998) Variance-Corrected MK + Sen's slope (guard rails: vc_mmkh)
vectorized_mann_kendall <- function(ts_matrix, invalid = "na") {
  n_pix <- nrow(ts_matrix)
  process_pixel <- function(i) {
    x <- as.numeric(ts_matrix[i, ]); x <- x[is.finite(x)]
    v <- vc_mmkh(x, invalid)
    tau <- v$tau  # identical VC tau definition in both pipelines (modifiedmk::mmkh)
    list(tau = tau, pval = v$p, slope = v$slope, filtered = !is.finite(v$p),
         status = v$status, nn = v$n_over_nstar)
  }
  do_par <- N_CORES > 1L && n_pix > 100L   # parallel workers must source this file
  res <- if (do_par && is_windows) parallel::parLapply(cl, seq_len(n_pix), process_pixel) else
    if (do_par) parallel::mclapply(seq_len(n_pix), process_pixel, mc.cores = N_CORES) else
      lapply(seq_len(n_pix), process_pixel)
  data.table::data.table(
    tau_vc          = vapply(res, `[[`, numeric(1),   "tau"),
    p_value_vc      = vapply(res, `[[`, numeric(1),   "pval"),
    sl_vc           = vapply(res, `[[`, numeric(1),   "slope"),
    filtered_vc     = vapply(res, `[[`, logical(1),   "filtered"),
    vc_status       = vapply(res, `[[`, character(1), "status"),
    vc_n_over_nstar = vapply(res, `[[`, numeric(1),   "nn"))
}
## J7b. TFPW Mann-Kendall (Yue et al. 2002) ----------------------------------------
vectorized_mann_kendall_tfpw <- function(ts_matrix) {
  n_pix <- nrow(ts_matrix)
  
  process_pixel_tfpw <- function(i) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < ANALYSIS_MIN_TREND_N || var(x, na.rm = TRUE) < 1e-6)
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

## J7c. PELT regime-shift --- pixel matrix ----------------------------------------
#' Applies detect_regime_shift_pelt() to annual averages of each pixel's
#' monthly time series.
#' NOTE: regime_shift_year is the first year AFTER the break, consistent
#'       with Pettitt's pettitt_first_year_after; name kept for CSV compatibility.
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
    if (is.na(fc) || fc >= length(yrs_clean)) next
    data.table::set(out, i, "regime_shift_year",     as.numeric(yrs_clean[fc + 1L]))
    data.table::set(out, i, "regime_shift_detected", TRUE)
    data.table::set(out, i, "n_changepoints",        res$n_changepoints)
    data.table::set(out, i, "mean_before_shift",     res$mean_before)
    data.table::set(out, i, "mean_after_shift",      res$mean_after)
    data.table::set(out, i, "magnitude_shift",       res$magnitude_shift)
  }
  out
}

## J7d. Spectral peak detection --- pixel matrix (same red-noise MC as J5)
SPECTRAL_RN_ENV_CACHE <- new.env(parent = emptyenv())
vectorized_spectral_peaks <- function(ts_matrix, min_period_months = 24L,
                                      sig_level = 0.95, n_sim_spectral = 500L) {
  n_pix <- nrow(ts_matrix)
  out <- rep(NA_integer_, n_pix)
  for (i in seq_len(n_pix)) {
    x <- as.numeric(ts_matrix[i, ])
    ok <- which(is.finite(x))
    if (!length(ok)) next
    x <- x[min(ok):max(ok)]
    if (length(x) < 48L || any(!is.finite(x)) ||
        !is.finite(stats::var(x)) || stats::var(x) <= 0) next
    ss <- .spectral_single(x, min_period_months = min_period_months,
                           sig_level = sig_level, n_sim = n_sim_spectral,
                           conf_cache_env = SPECTRAL_RN_ENV_CACHE,
                           sample_interval_months = 1L)
    out[i] <- ss$n_peaks
  }
  data.table::data.table(n_spectral_peaks = out)
}


#source("utils_teleconnection_addon.R")
################################################################################
# K2. NECHAKO REGIONAL SPATIAL FRAMEWORK
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
  cache_key <- paste0("subbasins|", normalizePath(subbasin_dir, winslash = "/", mustWork = FALSE))
  if (!is.null(.NECHAKO_GEOM_CACHE[[cache_key]])) return(.NECHAKO_GEOM_CACHE[[cache_key]])
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
  assign(cache_key, regions, envir = .NECHAKO_GEOM_CACHE)
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
        "Detected %.6f km2 of positive-area pairwise subbasin overlap.",
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
  # Section O's vectorised implementation is the one numerical kernel for
  # basin and regional means; the public wrapper is retained for Script 11
  # and upstream producer scripts.
  .bs_mean_cov(ts_matrix, area_weights, tw, min_cov = min_coverage)$value
}

## Whole-basin fractional-coverage weights on a native raster grid.
## Basin = dissolved union of the 7 subbasins; mask = coverage > 0;
## weights = cell area (m2) x coverage (0 outside the basin).
basin_fractional_weights <- function(template, basin = NULL) {
  if (!inherits(template, "SpatRaster")) stop("template must be a terra SpatRaster.")
  template <- template[[1]]
  if (is.na(terra::crs(template)) || !nzchar(terra::crs(template)))
    stop("basin_fractional_weights(): template has no CRS.")
  if (is.null(basin)) basin <- load_basin_vect()
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

## Canonical native-grid basin mask for Scripts 1, 5 and 8.
## A mask value of 1 means positive fractional basin coverage; 0 means outside.
## Return a packed terra raster because callers use terra::unwrap().
## The fractional-area weights remain defined by basin_fractional_weights().
save_common_analysis_grid <- function(template, basin = NULL, source = NULL) {
  if (!inherits(template, "SpatRaster"))
    stop("save_common_analysis_grid(): template must be a terra SpatRaster.")
  template <- template[[1L]]
  bw <- basin_fractional_weights(template, basin)
  if (!is.finite(bw$n_cells) || bw$n_cells < 1L)
    stop("save_common_analysis_grid(): basin has no covered native-grid cells.")
  if (!isTRUE(terra::compareGeom(bw$mask_raster, template,
                                 stopOnError = FALSE, crs = TRUE,
                                 ext = TRUE, rowcol = TRUE, res = TRUE)))
    stop("save_common_analysis_grid(): basin mask geometry differs from template.")
  list(
    basin_mask   = terra::wrap(bw$mask_raster),
    basin_pixels = as.integer(bw$n_cells),
    source       = source
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

################################################################################
# N. COMPANION-ANALYSIS HELPERS (scripts 20a-26)
#
# One definition of each helper the companion scripts share, so that date
# parsing, seasons, regional inputs, bootstrap resampling and trend inference
# behave identically everywhere. Periods are set in Section A
# (NECHAKO_RECORD_*, NECHAKO_RECENT_*, NECHAKO_WARM_SEASON_MONTHS).
################################################################################

## Meteorological season label only; a label does not define a season year.
## Use season_year() when grouping, so Dec(Y-1), Jan(Y), Feb(Y) share year Y.
season_of <- function(m) {
  m <- as.integer(m)
  out <- rep(NA_character_, length(m))
  out[m %in% c(12L, 1L, 2L)] <- "DJF"
  out[m %in% 3:5]            <- "MAM"
  out[m %in% 6:8]            <- "JJA"
  out[m %in% 9:11]           <- "SON"
  out
}

## Contiguous meteorological-season year, named for its Jan-Feb winter year.
## Complexity: O(n) time and O(n) storage; NA years/months remain NA.
season_year <- function(year, month) {
  year <- as.integer(year); month <- as.integer(month)
  n <- max(length(year), length(month))
  if (!length(year) || !length(month) ||
      !(length(year) %in% c(1L, n)) || !(length(month) %in% c(1L, n)))
    stop("year and month must be equal-length or scalar/non-empty")
  if (any(!is.na(month) & !(month %in% 1:12)))
    stop("month must lie in 1:12 or be NA")
  result <- rep_len(year, n) + as.integer(rep_len(month, n) == 12L)
  result[is.na(rep_len(year, n)) | is.na(rep_len(month, n))] <- NA_integer_
  result
}

## Parse monthly dates written as YYYY-MM or YYYY-MM-DD (character, factor,
## Date or IDate) and return the first day of each month. Anything else stops
## with the offending values, instead of silently becoming NA.
normalise_month_date <- function(x, label = "date") {
  if (inherits(x, "Date")) {
    out <- as.Date(x)
  } else {
    x <- trimws(as.character(x))
    if (!length(x)) stop(label, " is empty.")
    is_ym  <- grepl("^\\d{4}-\\d{2}$", x)
    is_ymd <- grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
    bad <- !(is_ym | is_ymd)
    if (any(bad))
      stop(label, " contains unsupported date format(s): ",
           paste(head(unique(x[bad]), 10L), collapse = ", "),
           ". Expected YYYY-MM or YYYY-MM-DD.")
    out <- as.Date(ifelse(is_ym, paste0(x, "-01"), x), format = "%Y-%m-%d")
  }
  if (anyNA(out))
    stop(label, " contains missing or invalid calendar date(s).")
  as.Date(format(out, "%Y-%m-01"))
}

## Region ids in the order returned by load_nechako_regions(), derived from
## the configuration only (no Shapefile I/O): whole basin first, then the 7
## subbasins.
nechako_region_ids <- function() {
  c(NECHAKO_REGION_BASIN_ID,
    vapply(tools::file_path_sans_ext(NECHAKO_SUBBASIN_FILES),
           sanitize_region_id, character(1), USE.NAMES = FALSE))
}

## Display name used for a region id in every companion output.
nechako_region_name <- function(region_id) {
  gsub("_", " ", as.character(region_id), fixed = TRUE)
}

## Stop unless `ids` is exactly the authoritative region set.
validate_region_ids <- function(ids, label = "regional input") {
  ids <- unique(as.character(ids))
  expected <- nechako_region_ids()
  missing <- setdiff(expected, ids)
  extra <- setdiff(ids, expected)
  if (length(missing) || length(extra))
    stop(label, " does not contain the authoritative Nechako region set.",
         if (length(missing)) paste0("\n  Missing: ", paste(missing, collapse = ", ")),
         if (length(extra)) paste0("\n  Unexpected: ", paste(extra, collapse = ", ")),
         "\n  Expected: ", paste(expected, collapse = ", "))
  invisible(TRUE)
}

## Read a long-format regional monthly CSV (region_id, date, <value_cols>).
##   * region_id is the join key; region_name is rebuilt from it so names are
##     identical in every companion output;
##   * dates are normalised to month starts and clipped to [start, end];
##   * duplicate region-months are collapsed only when they agree (a row of
##     NAs next to a finite row is redundant); conflicting finite values stop;
##   * the full region set and, optionally, a complete monthly record per
##     region are required.
read_regional_monthly_csv <- function(file, value_cols, label = basename(file),
                                      start = NECHAKO_RECORD_START,
                                      end = NECHAKO_RECORD_END,
                                      require_complete = TRUE,
                                      require_regions = TRUE) {
  if (!file.exists(file)) stop(label, " not found: ", file)
  x <- data.table::fread(file)
  req <- c("region_id", "date", value_cols)
  miss <- setdiff(req, names(x))
  if (length(miss))
    stop(label, " is missing column(s): ", paste(miss, collapse = ", "), " | ", file)
  x <- x[, req, with = FALSE]
  x[, region_id := trimws(as.character(region_id))]
  x[, date := normalise_month_date(date, paste(label, "date"))]
  for (v in value_cols) data.table::set(x, j = v, value = as.numeric(x[[v]]))
  x <- x[date >= as.Date(start) & date <= as.Date(end)]
  if (!nrow(x)) stop(label, " has no rows between ", start, " and ", end, ".")
  
  if (anyDuplicated(x, by = c("region_id", "date"))) {
    n_before <- nrow(x)
    x <- x[, {
      res <- lapply(value_cols, function(v) {
        vals <- unique(get(v)[is.finite(get(v))])
        if (length(vals) > 1L)
          stop(label, ": conflicting duplicate values of ", v, " for region ",
               region_id[1], " at ", format(date[1], "%Y-%m"), ".")
        if (length(vals)) vals else NA_real_
      })
      names(res) <- value_cols
      res
    }, by = .(region_id, date)]
    warning(label, ": collapsed ", n_before - nrow(x),
            " redundant duplicate region-month row(s).", call. = FALSE)
  }
  
  if (require_regions) validate_region_ids(x$region_id, label)
  if (require_complete) {
    for (rid in unique(x$region_id))
      validate_monthly_record(x$date[x$region_id == rid], start, end,
                              label = paste0(label, " [", rid, "]"),
                              require_complete = TRUE)
  }
  x[, region_name := nechako_region_name(region_id)]
  data.table::setcolorder(x, c("region_id", "region_name", "date"))
  data.table::setorder(x, region_id, date)
  x[]
}

## Circular moving-block bootstrap indices for a contiguous monthly series.
circular_block_indices <- function(n, block_length = 12L) {
  n <- as.integer(n)
  block_length <- as.integer(block_length)
  if (n <= 1L) return(seq_len(n))
  if (block_length < 1L) stop("block_length must be >= 1.")
  block_length <- min(block_length, n)
  starts <- sample.int(n, ceiling(n / block_length), replace = TRUE)
  idx <- unlist(lapply(starts, function(st)
    ((st - 1L + seq_len(block_length) - 1L) %% n) + 1L), use.names = FALSE)
  idx[seq_len(n)]
}

## Year-cluster bootstrap indices: whole years are resampled with
## replacement. Used for seasonal or filtered subsets, which are not
## contiguous monthly series. `years` gives the cluster of each row.
year_cluster_indices <- function(years) {
  groups <- split(seq_along(years), years)
  unlist(groups[sample.int(length(groups), length(groups), replace = TRUE)],
         use.names = FALSE)
}

## Bootstrap design shared by the companion scripts: contiguous periods use
## the circular monthly block bootstrap; filtered subsets use the year-cluster
## bootstrap and need at least `min_years` distinct years.
bootstrap_design <- function(contiguous, years, min_years = 8L) {
  if (contiguous) return("circular_monthly_block")
  if (length(unique(years[is.finite(years)])) < min_years)
    return("year_cluster_insufficient_years")
  "year_cluster"
}

bootstrap_indices <- function(design, n, years, block_length = 12L) {
  switch(design,
         circular_monthly_block = circular_block_indices(n, block_length),
         year_cluster = year_cluster_indices(years),
         stop("No bootstrap indices for design: ", design))
}

## OLS trend of y on x with Newey-West (HAC) standard errors; t-based p-value
## and 95% CI on the residual degrees of freedom. Same inference as
## lmtest::coeftest(fit, vcov. = sandwich::NeweyWest(...)).
hac_trend_test <- function(y, x, lag = 3L) {
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]; x <- x[ok]
  n <- length(y)
  empty <- data.frame(n = n, slope = NA_real_, hac_se = NA_real_,
                      hac_t = NA_real_, hac_df = NA_integer_, hac_p = NA_real_,
                      hac_ci95_low = NA_real_, hac_ci95_high = NA_real_,
                      ols_p = NA_real_, hac_lag = as.integer(lag))
  if (n < 5L || stats::sd(x) == 0) return(empty)
  if (!requireNamespace("sandwich", quietly = TRUE))
    stop("Package 'sandwich' is required for HAC trend inference.")
  fit <- stats::lm(y ~ x)
  V <- sandwich::NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE)
  b <- unname(stats::coef(fit)[2])
  se <- sqrt(V[2, 2])
  df <- fit$df.residual
  tval <- b / se
  q <- stats::qt(0.975, df)
  data.frame(n = n, slope = b, hac_se = se, hac_t = tval, hac_df = df,
             hac_p = 2 * stats::pt(abs(tval), df, lower.tail = FALSE),
             hac_ci95_low = b - q * se, hac_ci95_high = b + q * se,
             ols_p = summary(fit)$coefficients[2, 4], hac_lag = as.integer(lag))
}
# ==============================================================================
# BASIN_SERIES_utils
BASIN_SERIES_REGIONAL_DIR <- NECHAKO_REGION_OUTPUT_DIR
BASIN_SERIES_MANIFEST_DIR <- file.path(WD_PATH, "basin_series_manifest")
BASIN_SERIES_TOL          <- 1e-6   # verify_basin_series(): allowed |stored - recomputed|

.BS_CACHE <- new.env(parent = emptyenv())

.bs_grid_key <- function(r) paste(
  terra::nrow(r), terra::ncol(r), paste(signif(as.vector(terra::ext(r)), 12), collapse = ","),
  trimws(terra::crs(r)), paste(signif(terra::origin(r), 12), collapse = ","), sep = "|")

# Basin and regional weights on one native grid (cached per grid).
.bs_weights <- function(template) {
  template <- template[[1]]
  key <- .bs_grid_key(template)
  if (!is.null(.BS_CACHE[[key]])) return(.BS_CACHE[[key]])
  basin <- load_basin_vect()
  if (!terra::same.crs(basin, template)) basin <- terra::project(basin, terra::crs(template))
  bw <- basin_fractional_weights(template, basin)
  
  if (length(bw$weights) != terra::ncell(template))
    stop("basin_fractional_weights() must return one weight per native cell.")
  regions <- load_nechako_regions(template = template)
  validate_nechako_regions(regions)
  out <- list(basin = bw, regions = region_fractional_weights(template, regions),
              catalog = region_catalog(regions))
  assign(key, out, envir = .BS_CACHE)
  out
}

# Weighted mean and valid-area coverage of each column (time step).
.bs_mean_cov <- function(mat, w, total, min_cov = BASIN_MIN_VALID_AREA_COVERAGE) {
  n <- ncol(mat)
  if (!is.finite(total) || total <= 0) return(list(value = rep(NA_real_, n), coverage = rep(NA_real_, n)))
  w <- as.numeric(w); w[!is.finite(w) | w < 0] <- 0
  f <- is.finite(mat); m0 <- mat; m0[!f] <- 0
  va <- colSums(f * w)                       # valid weight per time step
  coverage <- va / total
  value <- ifelse(va > 0 & coverage >= min_cov, colSums(m0 * w) / va, NA_real_)
  list(value = as.numeric(value), coverage = as.numeric(coverage))
}

basin_series_regional_file <- function(index, scale, resolution = "monthly", variant = NULL) {
  dir <- if (is.null(variant)) BASIN_SERIES_REGIONAL_DIR else file.path(BASIN_SERIES_REGIONAL_DIR, variant)
  file.path(dir, sprintf("%s_%02d_regional_averaged_by_%s.csv", tolower(index), as.integer(scale),
                         if (identical(resolution, "annual")) "year" else "month"))
}

# mat: cells x time. cells = native cell numbers of the rows (NULL = all cells in order).
compute_basin_series <- function(mat, template, cells = NULL) {
  template <- template[[1]]; n_all <- terra::ncell(template)
  mat <- as.matrix(mat); mat[!is.finite(mat)] <- NA_real_
  if (is.null(cells)) {
    if (nrow(mat) != n_all)
      stop("compute_basin_series: ", nrow(mat), " rows but ", n_all, " native cells; pass 'cells'.")
    cells <- seq_len(n_all)
  } else if (length(cells) != nrow(mat)) stop("compute_basin_series: length(cells) != nrow(mat).")
  W <- .bs_weights(template)
  regional <- lapply(names(W$regions), function(rid) {
    rw <- W$regions[[rid]]; w_all <- as.numeric(rw$weights)
    tot <- sum(w_all[is.finite(w_all) & w_all > 0])
    s <- .bs_mean_cov(mat, w_all[cells], tot)
    list(region_id = rid, region_name = as.character(rw$region_name[1]), value = s$value,
         coverage = s$coverage, coverage_area_m2 = tot, n_cells = rw$n_cells, n_partial_cells = rw$n_partial)
  })
  list(basin = .bs_mean_cov(mat, W$basin$weights[cells], W$basin$area_m2),
       regional = regional, basin_meta = W$basin, catalog = W$catalog)
}

# Called by each index script once per index x scale, right after the pixel values exist.
write_basin_series <- function(mat, dates, template, index, scale, out_dir, cells = NULL,
                               resolution = c("monthly", "annual"), variant = NULL,
                               producer = NA_character_) {
  resolution <- match.arg(resolution)
  index <- tolower(index); scale <- as.integer(scale); dates <- as.Date(dates)
  if (length(dates) != ncol(as.matrix(mat))) stop("write_basin_series: length(dates) != ncol(mat).")
  s <- compute_basin_series(mat, template, cells)
  tag <- sprintf("%s_%02d", index, scale)
  reg_file <- basin_series_regional_file(index, scale, resolution, variant)
  for (d in c(out_dir, dirname(reg_file), BASIN_SERIES_MANIFEST_DIR))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  cat_file <- file.path(BASIN_SERIES_REGIONAL_DIR, "region_catalog.csv")
  if (!file.exists(cat_file)) data.table::fwrite(s$catalog, cat_file)
  
  reg <- data.table::rbindlist(lapply(s$regional, function(r) data.table::data.table(
    region_id = r$region_id, region_name = r$region_name, date = dates, value = r$value,
    valid_area_coverage_fraction = r$coverage, coverage_area_m2 = r$coverage_area_m2,
    n_cells = r$n_cells, n_partial_cells = r$n_partial_cells)))
  reg[, `:=`(index_type = index, scale = scale)]
  
  if (resolution == "monthly") {
    yr <- as.integer(format(dates, "%Y")); mo <- as.integer(format(dates, "%m")); uy <- sort(unique(yr))
    to_wide <- function(v) {
      w <- data.frame(Year = uy)
      for (m in 1:12) { col <- rep(NA_real_, length(uy)); h <- which(mo == m); col[match(yr[h], uy)] <- v[h]
      w[[month.abb[m]]] <- col }
      w
    }
    files <- c(file.path(out_dir, sprintf("%s_basin_averaged_by_month.csv", tag)),
               file.path(out_dir, sprintf("%s_basin_valid_area_coverage_by_month.csv", tag)), reg_file)
    utils::write.csv(to_wide(s$basin$value), files[1], row.names = FALSE, na = "")
    utils::write.csv(to_wide(s$basin$coverage), files[2], row.names = FALSE, na = "")
    data.table::fwrite(reg, reg_file)
  } else {
    hy <- as.integer(format(dates, "%Y"))
    files <- c(file.path(out_dir, sprintf("%s_basin_averaged_by_year.csv", tag)), reg_file)
    data.table::fwrite(data.table::data.table(hydro_year_end = hy, value = s$basin$value,
                                              valid_area_coverage_fraction = s$basin$coverage), files[1])
    reg[, hydro_year_end := as.integer(format(date, "%Y"))]
    data.table::fwrite(reg[, .(region_id, region_name, hydro_year_end, value, valid_area_coverage_fraction,
                               index_type, scale)], reg_file)
  }
  
  v <- s$basin$value; sdv <- stats::sd(v, na.rm = TRUE)
  meta <- data.table::data.table(
    index_type = index, scale = scale, variant = if (is.null(variant)) NA_character_ else variant,
    resolution = resolution, producer = producer, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    first = format(min(dates)), last = format(max(dates)), n_steps = length(dates), n_finite = sum(is.finite(v)),
    basin_mean = mean(v, na.rm = TRUE), basin_sd = sdv,
    n_cells = s$basin_meta$n_cells, n_partial = s$basin_meta$n_partial,
    basin_area_km2 = s$basin_meta$area_m2 / 1e6, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE,
    out_dir = out_dir, files = paste(basename(files), collapse = "; "))
  data.table::fwrite(meta, file.path(BASIN_SERIES_MANIFEST_DIR,
                                     paste0(tag, if (!is.null(variant)) paste0("_", variant), ".csv")))
  cat(sprintf("  [basin series] %s%s: %d steps | basin SD %.2f%s\n", tag,
              if (!is.null(variant)) paste0(" (", variant, ")") else "", length(dates), sdv,
              if (is.finite(sdv) && sdv < 0.9) "  (below 1: averaging standardized pixels shrinks variance)" else ""))
  invisible(list(basin = s$basin, regional = reg, manifest = meta))
}

# Reader. times = Date vector (monthly) or integer hydrological years (annual); when given,
# the series are aligned to it and a stale / incomplete file stops the run.
read_basin_series <- function(index, scale, out_dir, times = NULL,
                              resolution = c("monthly", "annual"), variant = NULL) {
  resolution <- match.arg(resolution)
  tag <- sprintf("%s_%02d", tolower(index), as.integer(scale))
  need <- function(f) {
    if (!file.exists(f)) stop("Basin series not found: ", f, "\n  Re-run the script that produces ",
                              toupper(index), "; it writes its own basin series.", call. = FALSE)
    f
  }
  if (resolution == "monthly") {
    rd <- function(f) {
      w <- utils::read.csv(need(f), check.names = FALSE); yrs <- as.integer(w$Year)
      list(time = as.Date(sprintf("%d-%02d-01", rep(yrs, each = 12L), rep(1:12, times = length(yrs)))),
           v = as.numeric(t(as.matrix(w[, month.abb]))))
    }
    b  <- rd(file.path(out_dir, sprintf("%s_basin_averaged_by_month.csv", tag)))
    cv <- rd(file.path(out_dir, sprintf("%s_basin_valid_area_coverage_by_month.csv", tag)))
    basin <- data.table::data.table(time = b$time, value = b$v, coverage = cv$v[match(b$time, cv$time)])
    reg <- data.table::fread(need(basin_series_regional_file(index, scale, "monthly", variant)))
    reg[, date := as.Date(date)]; tcol <- "date"
  } else {
    b <- data.table::fread(need(file.path(out_dir, sprintf("%s_basin_averaged_by_year.csv", tag))))
    basin <- data.table::data.table(time = as.integer(b$hydro_year_end), value = as.numeric(b$value),
                                    coverage = as.numeric(b$valid_area_coverage_fraction))
    reg <- data.table::fread(need(basin_series_regional_file(index, scale, "annual", variant)))
    reg[, hydro_year_end := as.integer(hydro_year_end)]; tcol <- "hydro_year_end"
  }
  if (!is.null(times)) {
    i <- match(times, basin$time)
    if (anyNA(i)) stop(sprintf("%s basin series does not cover the %d requested time steps (stale file?).",
                               tag, length(times)), call. = FALSE)
    basin <- basin[i]
    reg <- data.table::rbindlist(lapply(split(reg, by = "region_id"), function(d) {
      j <- match(times, d[[tcol]])
      if (anyNA(j)) stop(sprintf("%s regional series for %s does not cover the requested period (stale file?).",
                                 tag, d$region_id[1]), call. = FALSE)
      d[j]
    }))
  }
  list(basin = basin, regional = reg)
}

# Recompute from the pixel matrix and stop if the stored series differs (stale or edited file).
verify_basin_series <- function(stored, mat, template, cells = NULL, label = "", tol = BASIN_SERIES_TOL) {
  s <- compute_basin_series(mat, template, cells)
  cmp <- function(a, b, what) {
    if (length(a) != length(b) || any(is.finite(a) != is.finite(b)))
      stop(label, ": ", what, " differs from the producer's file (missing-value pattern). ",
           "Re-run the producing script.", call. = FALSE)
    d <- suppressWarnings(max(abs(a - b), na.rm = TRUE))
    if (is.finite(d) && d > tol)
      stop(sprintf("%s: %s differs from the producer's file (max |diff| %.3g). Re-run the producing script.",
                   label, what, d), call. = FALSE)
  }
  cmp(stored$basin$value, s$basin$value, "basin mean")
  for (r in s$regional) cmp(stored$regional[region_id == r$region_id]$value, r$value, paste("region", r$region_id))
  invisible(TRUE)
}

# All manifests in one table (provenance / freshness checks in Scripts 10 and 14).
basin_series_manifest <- function() {
  f <- list.files(BASIN_SERIES_MANIFEST_DIR, pattern = "\\.csv$", full.names = TRUE)
  if (!length(f)) return(data.table::data.table())
  data.table::rbindlist(lapply(f, data.table::fread), fill = TRUE)
}

cat("[OK] DROUGHT_ANALYSIS_utils.R loaded\n")