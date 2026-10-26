####################################################################################
# 4pr_pet_trends.r  - Trend Analysis: Climate, PM-PET Controls, PET, AET, Snow, Runoff
# ---------------------------------------------------------------------------------
# Performs comprehensive VC-MK + TFPW-MK + Spectral + Changepoint trend analysis
# for ALL BASIN PIXELS and SPECIFIC POINTS for climate variables, PM-PET controls,
# PET-method comparison diagnostics:
#   - Precipitation   (ERA5-Land, mm/month -> mm/year, full-year sum)
#   - PET             (Penman-Monteith, mm/month -> warm-season sum Apr-Oct)
#   - PET_ERA5        (ERA5-Land PEV, mm/month - warm-season sum Apr-Oct)
#   - AET             (ERA5-Land total evaporation, mm/month -> calendar-year sum)
#   - Temperature     (ERA5-Land 2 m Tair, K -> degC, annual mean + 12 calendar months)
#   - SWE, SWE_Peak   (ERA5-Land snow water equivalent)
#   - Runoff, Snowfall, SnowDepth, SnowDensity (ERA5-Land)
#
#
# OUTPUTS (binary .rds + CSV):
#   all_results.rds              pixel + basin trend statistics for all configured variables
#   analysis_metadata.rds        time series, extents, parameters, QA results
#   clipped_template.rds         native ERA5-Land precipitation raster template for visualization
#   summary_statistics.csv       per-variable significance counts
#   point_trend_stats.csv        per-point stats (annual + monthly; all configured variables)
#   point_monthly_timeseries.csv raw monthly values for all configured variables at each point
#   pet_source_qa.csv, swe_input_qa.csv
#   data_processing.log
#   regional/*                   regional monthly/annual means (incl. PM controls)
#   swe_homogeneity/*            SWE step-change / homogeneity tables, pixel scan, GeoTIFFs, figure
#   homogeneity/*                multi-variable homogeneity tables, pixel scans, GeoTIFFs, figure
#   negative_clamp_qa.csv        negatives before clamping (count, most negative, share in Jun-Oct)
#   root_cause_filter_reasons.csv, vc_correction_tags.csv, root_cause_summary.csv
#   basin_boundary.rds, native_basin_mask.rds
####################################################################################

# rm(list = ls())  # commented out: safe to source from a master script
# gc()             # commented out: safe to source from a master script

library(ncdf4)
library(terra)
library(data.table)
library(Kendall)
library(zoo)
library(sf)
library(changepoint)
library(trend)
library(modifiedmk)

# -- Shared statistical functions (VC-MK, TFPW-MK, PELT, spectral) --
# check_min_value_threshold(), calculate_variance_with_ties(),
# detect_regime_shift_pelt(), and mk_tfpw_spectral_for_series() are defined in
# DROUGHT_ANALYSIS_utils.R (Section J).
# ================= USER / ENV =================
WD_PATH <- Sys.getenv(
  "NECHAKO_WD",
  unset = "D:/Nechako_Drought/Nechako/"
)
if (!dir.exists(WD_PATH)) {
  stop("Working directory not found: ", WD_PATH)
}
setwd(WD_PATH)
source(file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R"))

out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ===== LOGGING SETUP =====
LOG_FILE <- file.path(out_dir, "data_processing.log")
log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, "   ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, "   ", msg))
}
cat(
  paste0(
    "\n============================================================\n",
    "Trend Analysis Data Processing - Started\n",
    "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "============================================================\n"
  ),
  file = LOG_FILE,
  append = TRUE
)
##
library(future)
library(future.apply)
num_cores <- min(
  max(1L, parallel::detectCores(logical = TRUE) - 1L),
  8L
)

future::plan(
  future::multisession,
  workers = num_cores
)

log_event(paste("Using", num_cores, "future multisession workers"))

####################################################################################
# -- SPECIFIC POINTS CONFIG
####################################################################################
SPECIFIC_PTS <- data.frame(
  x = c(-124.5, -123.0, -125.0),
  y = c(  54.0,   53.8,   54.2)
)

####################################################################################
# -- PARAMETERS
####################################################################################
alpha                    <- 0.05
n_sim_spectral           <- 500
max_tie_percent          <- 95
max_min_value_pct_precip <- 50
max_min_value_pct_pet    <- 50
min_positive_value       <- 0.01
min_obs_changepoint      <- 20
pelt_minseglen           <- 5L     # no PELT break within 5 yr of either end
# Monthly PET: also reject series whose floor-value (0.01 mm) share exceeds
# max_min_value_pct_pet.
PET_MONTHLY_ALSO_PCT_FILTER <- TRUE

warm_months_pet          <- 4:10
min_nonzero_pet_monthly  <- 10

# ===== INPUT PATHS =====
precip_path   <- "monthly_data_direct/total_precipitation_monthly.nc"
pet_path      <- "monthly_data_direct/potential_evapotranspiration_monthly.nc"
# ERA5-Land PEV PET - produced by 2b_PET_ERALand.R
pet_era5_path  <- "monthly_data_direct/potential_evaporation_era5land_monthly.nc"
# ERA5-Land total evaporation (e), used as AET.
AET_NC_CANDIDATES <- c(
  "monthly_data_direct/total_evaporation_monthly.nc",        # current download
  "monthly_data_direct/total_evapotranspiration_monthly.nc"  # older name
)
aet_path <- AET_NC_CANDIDATES[file.exists(AET_NC_CANDIDATES)][1]
if (is.na(aet_path)) aet_path <- AET_NC_CANDIDATES[1]

TAIR_NC_CANDIDATES <- c(
  "monthly_data_direct/2m_temperature_monthly.nc",
  "monthly_data_direct/t2m_monthly.nc",
  "monthly_data_direct/air_temperature_2m_monthly.nc",
  "monthly_data_direct/ERA5Land_T2m_monthly.nc"
)

# PM-PET meteorological-control fields exported by 2b_PET_ERALand.R.
PM_CONTROL_NC <- c(
  Temperature = "monthly_data_direct/pm_control_2m_temperature_monthly.nc",
  NetRadiation = "monthly_data_direct/pm_control_net_radiation_monthly.nc",
  VPD = "monthly_data_direct/pm_control_vpd_monthly.nc",
  Wind2m = "monthly_data_direct/pm_control_wind2m_monthly.nc",
  SurfacePressure = "monthly_data_direct/pm_control_surface_pressure_monthly.nc"
)

PM_CONTROL_LABELS <- c(
  Temperature = "PM_Temperature",
  NetRadiation = "PM_NetRadiation",
  VPD = "PM_VPD",
  Wind2m = "PM_Wind2m",
  SurfacePressure = "PM_SurfacePressure"
)

PM_CONTROL_UNITS <- c(
  Temperature = "degC",
  NetRadiation = "MJ/m2/day",
  VPD = "kPa",
  Wind2m = "m/s",
  SurfacePressure = "kPa"
)

# ===== SWE (ERA5-Land snow water equivalent): TREND + HOMOGENEITY CONFIG =====
swe_path          <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
SWE_UNITS         <- "m"          # FALLBACK only: used when the file carries no units attribute
RUN_SWE           <- file.exists(swe_path)
SWE_ANNUAL_MONTHS <- 3:4          # annual SWE metric ("SWE") = mean of Mar-Apr monthly-mean SWE
SWE_PEAK_MONTHS   <- 1:6          # "SWE_Peak" = maximum monthly-mean SWE over Jan-Jun
SWE_STEP_YEAR     <- 1979L        # a-priori break hypothesis: FIRST year of the later regime
SWE_ACCUM_MONTHS  <- c(11L, 12L, 1L, 2L, 3L)  # Nov(Y-1)-Mar(Y)
N_SIM_HOMOG       <- 10000L       # Monte-Carlo size for SNHT / Buishand p-values (SWE block)
swe_out_dir       <- file.path(out_dir, "swe_homogeneity")

# ===== ADDITIONAL ERA5-LAND VARIABLES: RUNOFF, SNOWFALL, SNOW DEPTH, SNOW DENSITY (v8) =====
SNOW_PRESENT_DEPTH_CM <- 1
EXTRA_VARS <- list(
  Runoff      = list(path = "monthly_data_direct/runoff_monthly.nc",       nc_var = "runoff",
                     conv = "accum_m_per_day", unit = "mm/month",
                     col_month = "runoff_mm_month",         col_year = "runoff_mm_year",
                     annual = "sum",    is_precip = TRUE,  skip_min_filter = FALSE),
  Snowfall    = list(path = "monthly_data_direct/snowfall_monthly.nc",     nc_var = "snowfall",
                     conv = "accum_m_per_day", unit = "mm w.e./month",
                     col_month = "snowfall_mm_month",       col_year = "snowfall_mm_year",
                     annual = "sum",    is_precip = TRUE,  skip_min_filter = FALSE),
  SnowDepth   = list(path = "monthly_data_direct/snow_depth_monthly.nc",   nc_var = "snow_depth",
                     conv = "m_to_cm", unit = "cm",
                     col_month = "snow_depth_cm_month",     col_year = "snow_depth_marapr_cm",
                     annual = "marapr", is_precip = TRUE,  skip_min_filter = FALSE),
  SnowDensity = list(path = "monthly_data_direct/snow_density_monthly.nc", nc_var = "snow_density",
                     conv = "none", unit = "kg m-3",
                     col_month = "snow_density_kgm3_month", col_year = "snow_density_marapr_kgm3",
                     annual = "marapr", is_precip = FALSE, skip_min_filter = TRUE)
)
.extra_missing <- names(EXTRA_VARS)[!vapply(EXTRA_VARS, function(v) file.exists(v$path), logical(1))]
if (length(.extra_missing))
  log_event(paste("  [EXTRA] Not found, skipped:", paste(.extra_missing, collapse = ", ")))
EXTRA_VARS <- EXTRA_VARS[setdiff(names(EXTRA_VARS), .extra_missing)]
log_event(paste("  [EXTRA] Variables to process:",
                if (length(EXTRA_VARS)) paste(names(EXTRA_VARS), collapse = ", ") else "none"))

# ===== MULTI-VARIABLE HOMOGENEITY CONFIG (v9) =====
homog_out_dir      <- file.path(out_dir, "homogeneity")
N_SIM_HOMOG_ALL    <- 5000L   # Monte-Carlo size for SNHT / Buishand in the multi-variable block
HOMOG_BREAK_WINDOW <- 3L      # breaks within +/- this many years of SWE_STEP_YEAR are flagged

log_event(paste("Using", num_cores, "cores for parallel processing"))

####################################################################################
# -- BASIN BOUNDARY - AUTHORITATIVE 7-SUBBASIN UNION
####################################################################################
log_event("Loading authoritative Nechako basin boundary from the 7 subbasin Shapefiles...")
regional_subbasins <- load_nechako_subbasins()
validate_nechako_regions(load_nechako_regions(subbasins=regional_subbasins))
basin_boundary <- sf::st_as_sf(nechako_subbasin_union(regional_subbasins))
basin_boundary <- sf::st_transform(basin_boundary, "EPSG:4326")
log_event(sprintf("Loaded dissolved union of %d subbasins; CRS=%s",
                  length(regional_subbasins), sf::st_crs(basin_boundary)$input))
saveRDS(basin_boundary, file.path(out_dir, "basin_boundary.rds"))

####################################################################################
# -- HELPER FUNCTIONS
####################################################################################
aggregate_to_annual_fast <- function(monthly_matrix, years, method = "sum") {
  y_levels      <- unique(years)
  yfac          <- match(years, y_levels)
  count_by_year <- rowsum((!is.na(monthly_matrix)) * 1L, group = yfac, reorder = FALSE)
  m0            <- monthly_matrix; m0[is.na(m0)] <- 0
  sum_by_year   <- rowsum(m0, group = yfac, reorder = FALSE)
  res           <- if (identical(method, "sum")) sum_by_year else
    sum_by_year / pmax(count_by_year, 1)
  res[count_by_year < 6] <- NA_real_
  rownames(res) <- y_levels
  res
}

aggregate_to_annual_warmseason_fast <- function(monthly_matrix, years, months_vec,
                                                warm_months = 4:10, method = "sum") {
  warm_idx    <- which(months_vec %in% warm_months)
  if (!length(warm_idx)) stop("No warm-season months found in months_vec.")
  warm_matrix <- monthly_matrix[warm_idx, , drop = FALSE]
  warm_years  <- years[warm_idx]
  y_levels    <- unique(warm_years)
  yfac        <- match(warm_years, y_levels)
  count       <- rowsum((!is.na(warm_matrix)) * 1L, group = yfac, reorder = FALSE)
  m0          <- warm_matrix; m0[is.na(m0)] <- 0
  sums        <- rowsum(m0, group = yfac, reorder = FALSE)
  res         <- if (identical(method, "sum")) sums else sums / pmax(count, 1)
  res[count < ceiling(length(warm_months) * 0.5)] <- NA_real_
  rownames(res) <- y_levels
  res
}

compute_month_index <- function(months) split(seq_along(months), months)

# ------------------------------------------------------------------------------
# Calendar verification (v9). Every monthly input must have exactly the precipitation
# calendar. A readable time axis must match by year-month; an unreadable one is accepted
# only when the layer count matches, and this is logged.
# ------------------------------------------------------------------------------
verify_monthly_calendar <- function(r, expected_dates, label) {
  n <- terra::nlyr(r)
  if (n != length(expected_dates))
    stop(sprintf("[%s] %d layers, but the precipitation record has %d months (%s to %s).",
                 label, n, length(expected_dates),
                 format(min(expected_dates), "%Y-%m"), format(max(expected_dates), "%Y-%m")))
  d <- tryCatch({
    td <- terra::time(r)
    if (is.null(td) || length(td) != n || anyNA(td)) NULL else as.Date(format(as.Date(td), "%Y-%m-01"))
  }, error = function(e) NULL)
  if (is.null(d)) {
    log_event(sprintf("  [CAL:%s] no readable time axis; layer count matches, using the precipitation calendar (%s to %s).",
                      label, format(min(expected_dates), "%Y-%m"), format(max(expected_dates), "%Y-%m")))
    return(invisible(expected_dates))
  }
  if (!all(d == expected_dates)) {
    bad <- which(d != expected_dates)[1]
    stop(sprintf("[%s] time axis differs from the precipitation calendar at layer %d: %s vs %s.",
                 label, bad, format(d[bad], "%Y-%m"), format(expected_dates[bad], "%Y-%m")))
  }
  log_event(sprintf("  [CAL:%s] %d months, %s to %s - matches the precipitation calendar.",
                    label, n, format(min(d), "%Y-%m"), format(max(d), "%Y-%m")))
  invisible(d)
}
# ------------------------------------------------------------------------------
# Negative-value QA before clamping (v10). Logs how many values are negative, the most
# negative value, how many exceed `tol` (not plausible as round-off), and the share in
# Jun-Oct; then sets negatives to 0. Results go to negative_clamp_qa.csv.
# ------------------------------------------------------------------------------
NEG_CLAMP_QA  <- list()
NEG_CLAMP_TOL <- c(SWE = 1, SnowDepth = 0.1, Runoff = 0.1, Snowfall = 0.1)   # in the converted units
clamp_negatives <- function(r, label, unit, tol = NEG_CLAMP_TOL[[label]]) {
  v   <- terra::values(r, mat = FALSE)            # cells x layers, layer-major
  nc  <- terra::ncell(r)
  ok  <- is.finite(v)
  neg <- which(ok & v < 0)
  mo  <- months[(neg - 1L) %/% nc + 1L]
  qa <- data.table::data.table(
    variable = label, unit = unit, n_values = sum(ok), n_negative = length(neg),
    pct_negative = 100 * length(neg) / max(sum(ok), 1L),
    most_negative = if (length(neg)) min(v[neg]) else NA_real_,
    p01_negative  = if (length(neg)) unname(stats::quantile(v[neg], 0.01)) else NA_real_,
    tolerance = tol, n_beyond_tolerance = sum(v[neg] < -tol),
    pct_negative_in_JunOct = if (length(neg)) 100 * mean(mo %in% 6:10) else NA_real_)
  NEG_CLAMP_QA[[label]] <<- qa
  if (!length(neg)) { log_event(sprintf("  [NEG:%s] no negative values", label)); return(r) }
  log_event(sprintf("  [NEG:%s] %d of %d values negative (%.1f%%); most negative %.3g %s; %d beyond -%g %s; %.0f%% of negatives in Jun-Oct",
                    label, qa$n_negative, qa$n_values, qa$pct_negative, qa$most_negative, unit,
                    qa$n_beyond_tolerance, tol, unit, qa$pct_negative_in_JunOct))
  if (qa$n_beyond_tolerance > 0)
    log_event(sprintf("  [NEG:%s] WARNING: %d values are below -%g %s - not round-off; inspect the source file.",
                      label, qa$n_beyond_tolerance, tol, unit))
  if (length(neg)) r <- terra::ifel(r < 0, 0, r)
  r
}
####################################################################################
# SWE HELPERS (annual metrics, season aggregation, homogeneity / step-change tests)
####################################################################################
# Annual metric from a monthly matrix (time x series). `fun` = "mean" or "max" over
# `use_months`; a year is NA unless >= `min_months` finite months are present.
annual_metric_fast <- function(mat, years, months, use_months, fun = c("mean", "max"),
                               min_months = length(use_months)) {
  fun <- match.arg(fun)
  if (!is.matrix(mat)) mat <- matrix(as.numeric(mat), ncol = 1L)
  y_all <- sort(unique(years))
  out   <- matrix(NA_real_, length(y_all), ncol(mat), dimnames = list(y_all, NULL))
  idx   <- which(months %in% use_months)
  for (k in seq_along(y_all)) {
    rows <- idx[years[idx] == y_all[k]]
    if (!length(rows)) next
    blk  <- mat[rows, , drop = FALSE]
    n_ok <- colSums(is.finite(blk))
    v <- if (fun == "mean") {
      colSums(blk, na.rm = TRUE) / pmax(n_ok, 1L)
    } else {
      suppressWarnings(apply(blk, 2L, max, na.rm = TRUE))
    }
    v[n_ok < min_months] <- NA_real_
    v[!is.finite(v)]     <- NA_real_
    out[k, ] <- v
  }
  out
}

# Season aggregation that crosses the calendar-year boundary (e.g. Nov(Y-1)-Mar(Y)).
season_by_year <- function(x, years, months, season_months, wrap_from = 11L,
                           fun = c("sum", "mean")) {
  fun  <- match.arg(fun)
  keep <- months %in% season_months
  sy   <- (years + as.integer(months >= wrap_from))[keep]
  agg  <- vapply(split(x[keep], sy), function(v) {
    if (sum(is.finite(v)) < length(season_months)) return(NA_real_)
    if (fun == "sum") sum(v) else mean(v)
  }, numeric(1))
  yrs_out <- sort(unique(years))
  out <- agg[as.character(yrs_out)]
  names(out) <- yrs_out
  out
}

.p_floor <- function(p, m) if (is.finite(p)) max(p, 1 / m) else NA_real_
.try_na  <- function(expr, default = NA_real_) tryCatch(suppressWarnings(expr), error = function(e) default)

# ------------------------------------------------------------------------------
# Single-series homogeneity / step-change battery.
#   * Unknown-date tests: Pettitt, SNHT, Buishand range.
#   * A-priori fixed-date comparison at `step_year` (first year of the later regime):
#     Welch t, Mann-Whitney, Cohen's d, Fligner-Killeen variance test.
# ------------------------------------------------------------------------------
homogeneity_tests <- function(x, yrs, step_year = SWE_STEP_YEAR, m_sim = N_SIM_HOMOG,
                              run_sim_tests = TRUE) {
  ok <- is.finite(x); x <- as.numeric(x[ok]); yrs <- as.integer(yrs[ok]); n <- length(x)
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
  if (n < 20L || length(unique(x)) < 5L) return(as.data.table(o))
  
  pt <- tryCatch(trend::pettitt.test(x), error = function(e) NULL)
  if (!is.null(pt)) {
    K <- as.integer(pt$estimate)[1]
    o$pettitt_p <- pt$p.value
    o$pettitt_last_year_before  <- yrs[K]
    o$pettitt_first_year_after  <- yrs[min(K + 1L, n)]
  }
  if (run_sim_tests) {
    sn <- tryCatch(trend::snh.test(x, m = m_sim), error = function(e) NULL)
    if (!is.null(sn)) { o$snht_p <- .p_floor(sn$p.value, m_sim)
    o$snht_first_year_after <- yrs[min(as.integer(sn$estimate)[1] + 1L, n)] }
    br <- tryCatch(trend::br.test(x, m = m_sim), error = function(e) NULL)
    if (!is.null(br)) o$buishand_range_p <- .p_floor(br$p.value, m_sim)
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
    o$welch_p        <- .try_na(stats::t.test(a, b)$p.value)
    o$mannwhitney_p  <- .try_na(stats::wilcox.test(a, b, exact = FALSE)$p.value)
    o$variance_fligner_p <- .try_na(stats::fligner.test(list(b, a))$p.value)
  }
  as.data.table(o)
}

# ------------------------------------------------------------------------------
# Trend with and without a step term, AR(1) errors (nlme::gls, ML so AIC is comparable)
# ------------------------------------------------------------------------------
step_trend_model <- function(x, yrs, step_year = SWE_STEP_YEAR) {
  d <- data.frame(y = as.numeric(x), t = as.integer(yrs) - min(as.integer(yrs), na.rm = TRUE),
                  step = as.integer(yrs >= step_year))
  d <- d[is.finite(d$y), , drop = FALSE]
  na <- data.table(trend_no_step = NA_real_, trend_no_step_p = NA_real_,
                   trend_with_step = NA_real_, trend_with_step_p = NA_real_,
                   step_coef = NA_real_, step_p = NA_real_, delta_AIC_step_vs_nostep = NA_real_,
                   ar1_phi_with_step = NA_real_)
  if (nrow(d) < 20L || length(unique(d$step)) < 2L) return(na)
  fit <- function(f) tryCatch(
    suppressWarnings(nlme::gls(f, data = d, correlation = nlme::corAR1(form = ~ t),
                               method = "ML")), error = function(e) NULL)
  m0 <- fit(y ~ t); m1 <- fit(y ~ t + step)
  if (is.null(m0) || is.null(m1)) return(na)
  tt0 <- summary(m0)$tTable; tt1 <- summary(m1)$tTable
  data.table(trend_no_step = tt0["t", "Value"], trend_no_step_p = tt0["t", "p-value"],
             trend_with_step = tt1["t", "Value"], trend_with_step_p = tt1["t", "p-value"],
             step_coef = tt1["step", "Value"], step_p = tt1["step", "p-value"],
             delta_AIC_step_vs_nostep = AIC(m1) - AIC(m0),   # negative = step model preferred
             ar1_phi_with_step = tryCatch(as.numeric(coef(m1$modelStruct$corStruct,
                                                          unconstrained = FALSE)),
                                          error = function(e) NA_real_))
}

# ------------------------------------------------------------------------------
# VC-MK + TFPW-MK (the pipeline's own run_trend_bundle) on sub-periods.
# ------------------------------------------------------------------------------
subperiod_trends <- function(x, yrs, series_label, region_id = NA_character_,
                             periods, is_precip = TRUE, skip_min_filter = FALSE) {
  rows <- lapply(names(periods), function(pn) {
    a <- periods[[pn]][1]; b <- periods[[pn]][2]
    idx <- yrs >= a & yrs <= b
    xs  <- as.numeric(x[idx])
    if (sum(is.finite(xs)) < 10L) {
      return(data.table(region_id = region_id, series = series_label, period = pn,
                        year_start = a, year_end = b, n = sum(is.finite(xs)),
                        mean_value = NA_real_, sen_slope_per_yr = NA_real_,
                        sen_pct_of_mean_per_decade = NA_real_,
                        tau_vc = NA_real_, p_vc = NA_real_, filtered_vc = TRUE,
                        tau_tfpw = NA_real_, p_tfpw = NA_real_, filtered_tfpw = TRUE))
    }
    r <- run_trend_bundle(xs, is_precip = is_precip, alpha = alpha,
                          max_tie_percent = max_tie_percent, n_sim_spectral = n_sim_spectral,
                          conf_env = new.env(parent = emptyenv()), start_year = as.integer(a),
                          skip_min_filter = skip_min_filter)
    mu <- mean(xs, na.rm = TRUE)
    data.table(region_id = region_id, series = series_label, period = pn,
               year_start = a, year_end = b, n = sum(is.finite(xs)), mean_value = mu,
               sen_slope_per_yr = r$tf$sl,
               sen_pct_of_mean_per_decade = if (isTRUE(abs(mu) > 0)) 100 * 10 * r$tf$sl / abs(mu) else NA_real_,
               tau_vc = r$vc$tau, p_vc = r$vc$p, filtered_vc = r$vc$filtered,
               tau_tfpw = r$tf$tau, p_tfpw = r$tf$p, filtered_tfpw = r$tf$filtered)
  })
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

# ------------------------------------------------------------------------------
# SWEI-style reference sensitivity for an annual series (Mar-Apr SWE proxy).
# ------------------------------------------------------------------------------
.gringorten_z_vec <- function(x, ref, eps = 1e-6) {
  ref <- sort(ref[is.finite(ref)]); out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (length(ref) < 10L || !any(ok)) return(out)
  n_lt <- findInterval(x[ok], ref, left.open = TRUE); n_le <- findInterval(x[ok], ref)
  p <- ((n_lt + n_le + 1) / 2 - 0.44) / (length(ref) + 0.12)
  out[ok] <- stats::qnorm(pmin(pmax(p, eps), 1 - eps)); out
}
reference_sensitivity <- function(x, yrs, region_id, series_label,
                                  refs = list(whole_record = c(1950L, 2025L),
                                              post_step = c(SWE_STEP_YEAR, 2025L),
                                              clim_1991_2020 = c(1991L, 2020L)),
                                  eras = list(pre_step = c(1950L, SWE_STEP_YEAR - 1L),
                                              early_post = c(SWE_STEP_YEAR, 1990L),
                                              clim_1991_2020 = c(1991L, 2020L),
                                              recent = c(2021L, 2025L))) {
  rbindlist(lapply(names(refs), function(rn) {
    ref <- x[yrs >= refs[[rn]][1] & yrs <= refs[[rn]][2]]
    z_g <- .gringorten_z_vec(x, ref)
    mu <- mean(ref, na.rm = TRUE); sdv <- stats::sd(ref, na.rm = TRUE)
    z_p <- (x - mu) / sdv
    rbindlist(lapply(names(eras), function(en) {
      i <- yrs >= eras[[en]][1] & yrs <= eras[[en]][2]
      data.table(region_id = region_id, series = series_label, reference = rn,
                 ref_start = refs[[rn]][1], ref_end = refs[[rn]][2],
                 era = en, era_start = eras[[en]][1], era_end = eras[[en]][2],
                 n_years = sum(is.finite(x[i])),
                 mean_z_gringorten = mean(z_g[i], na.rm = TRUE),
                 mean_z_parametric = mean(z_p[i], na.rm = TRUE),
                 pct_years_below_minus1 = 100 * mean(z_g[i] < -1, na.rm = TRUE))
    }))
  }))
}

####################################################################################
# -- ROOT-CAUSE DIAGNOSTICS / QUALITY SCREEN
####################################################################################
diag_store <- new.env(parent = emptyenv())
diag_store$qa_tables <- list()

series_quality_profile <- function(x, floor_value = NA_real_) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) {
    return(list(
      n = 0L, sd = NA_real_, uniq = 0L, zero_frac = NA_real_,
      floor_frac = NA_real_, min = NA_real_, median = NA_real_, max = NA_real_
    ))
  }
  list(
    n = n,
    sd = if (n >= 2L) stats::sd(x) else NA_real_,
    uniq = length(unique(x)),
    zero_frac = mean(x == 0),
    floor_frac = if (is.finite(floor_value)) mean(x <= floor_value) else NA_real_,
    min = min(x),
    median = stats::median(x),
    max = max(x)
  )
}

log_series_quality <- function(x, label, floor_value = NA_real_) {
  prof <- series_quality_profile(x, floor_value = floor_value)
  floor_txt <- if (is.finite(prof$floor_frac)) sprintf("%.2f%%", 100 * prof$floor_frac) else "NA"
  log_event(sprintf(
    "  [QA:%s] n=%d | sd=%s | unique=%d | zero_frac=%s | floor_frac=%s | min=%.4g | median=%.4g | max=%.4g",
    label,
    prof$n,
    ifelse(is.na(prof$sd), "NA", format(signif(prof$sd, 4), scientific = TRUE)),
    prof$uniq,
    ifelse(is.na(prof$zero_frac), "NA", sprintf("%.2f%%", 100 * prof$zero_frac)),
    floor_txt,
    prof$min, prof$median, prof$max
  ))
  invisible(prof)
}

matrix_quality_profile <- function(mat, label, floor_value = NA_real_) {
  if (is.null(mat) || !length(mat)) {
    log_event(sprintf("  [QA:%s] matrix empty", label))
    return(invisible(NULL))
  }
  
  n_time   <- nrow(mat)
  n_series <- ncol(mat)
  
  sd_vals <- vapply(seq_len(n_series), function(i) {
    x <- mat[, i]
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    stats::sd(x)
  }, numeric(1))
  
  uniq_vals <- vapply(seq_len(n_series), function(i) {
    x <- mat[, i]
    length(unique(x[is.finite(x)]))
  }, integer(1))
  
  zero_frac <- vapply(seq_len(n_series), function(i) {
    x <- mat[, i]
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    mean(x == 0)
  }, numeric(1))
  
  floor_frac <- if (is.finite(floor_value)) {
    vapply(seq_len(n_series), function(i) {
      x <- mat[, i]
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      mean(x <= floor_value)
    }, numeric(1))
  } else {
    rep(NA_real_, n_series)
  }
  
  out <- data.table(
    label = label,
    n_time = n_time,
    n_series = n_series,
    pct_na = 100 * sum(!is.finite(mat)) / length(mat),
    median_sd = stats::median(sd_vals, na.rm = TRUE),
    n_zero_sd = sum(sd_vals <= 1e-6, na.rm = TRUE),
    n_low_sd = sum(sd_vals < 0.01, na.rm = TRUE),
    median_unique = stats::median(uniq_vals, na.rm = TRUE),
    n_very_tied = sum(uniq_vals <= 3L, na.rm = TRUE),
    median_zero_frac = stats::median(zero_frac, na.rm = TRUE),
    median_floor_frac = stats::median(floor_frac, na.rm = TRUE),
    n_floor_dominant = sum(floor_frac > 0.5, na.rm = TRUE)
  )
  
  floor_txt <- if (is.finite(floor_value)) sprintf(" | floor>%.3g=%d", floor_value, out$n_floor_dominant) else ""
  log_event(sprintf(
    "  [QA:%s] series=%d time=%d | NA=%.2f%% | median sd=%s | sd<1e-6=%d | sd<0.01=%d | median unique=%s | uniq<=3=%d | median zero=%.2f%%%s",
    label,
    n_series, n_time, out$pct_na,
    ifelse(is.na(out$median_sd), "NA", format(signif(out$median_sd, 4), scientific = TRUE)),
    out$n_zero_sd, out$n_low_sd,
    ifelse(is.na(out$median_unique), "NA", format(out$median_unique)),
    out$n_very_tied,
    100 * out$median_zero_frac,
    floor_txt
  ))
  
  diag_store$qa_tables[[length(diag_store$qa_tables) + 1L]] <- out
  invisible(out)
}

log_result_root_cause_summary <- function(results_dt) {
  if (is.null(results_dt) || !nrow(results_dt)) {
    log_event("  [DIAG] No results table available for root-cause summary.")
    return(invisible(NULL))
  }
  log_event(sprintf("  [DIAG] non-finite p-values among retained series: VC=%d | TFPW=%d",
                    results_dt[!filtered_vc & !is.finite(p_value_vc), .N],
                    results_dt[!filtered_tfpw & !is.finite(p_value_tfpw), .N]))
  
  # 1. Series actually REMOVED (filtered = TRUE), by reason
  removed <- rbindlist(list(
    results_dt[filtered_vc %in% TRUE, .N, by = .(variable, period, is_basin_average, reason = filter_reason_vc)][, method := "VC"],
    results_dt[filtered_tfpw %in% TRUE, .N, by = .(variable, period, is_basin_average, reason = filter_reason_tfpw)][, method := "TFPW"]))
  setnames(removed, "N", "n")
  fwrite(removed[order(-n)], file.path(out_dir, "root_cause_filter_reasons.csv"))
  if (nrow(removed)) {
    log_event("  [DIAG] series REMOVED from trend testing (top 8):")
    top <- removed[order(-n)][seq_len(min(8L, .N))]
    for (i in seq_len(nrow(top)))
      log_event(sprintf("    [%s/%s/%s] %s = %d", top$variable[i], top$period[i], top$method[i], top$reason[i], top$n[i]))
  } else log_event("  [DIAG] no series removed.")
  
  # 2. Series RETAINED but carrying a VC tag (variance correction not applied)
  #    vc_shrinkage_capped_uncorrected : Hamed-Rao would shrink the variance (n/n* < 1);
  #                                      the uncorrected variance was kept (conservative)
  #    vc_invalid_variance_uncorrected : corrected variance undefined; uncorrected kept
  tags <- results_dt[!(filtered_vc %in% TRUE) & !is.na(filter_reason_vc) & filter_reason_vc != "none",
                     .N, by = .(variable, period, is_basin_average, tag = filter_reason_vc)]
  setnames(tags, "N", "n")
  fwrite(tags[order(-n)], file.path(out_dir, "vc_correction_tags.csv"))
  if (nrow(tags)) log_event(sprintf(
    "  [DIAG] %d retained series carry a VC tag (uncorrected variance kept - informational, not a removal); see vc_correction_tags.csv",
    sum(tags$n)))
  
  diag_summary <- results_dt[, .(
    n_total = .N,
    n_filtered_vc = sum(filtered_vc, na.rm = TRUE),
    n_filtered_tfpw = sum(filtered_tfpw, na.rm = TRUE),
    n_vc_tagged_uncorrected = sum(!(filtered_vc %in% TRUE) & !is.na(filter_reason_vc) & filter_reason_vc != "none"),
    n_nonfinite_p_vc = sum(!filtered_vc & !is.finite(p_value_vc), na.rm = TRUE),
    n_nonfinite_p_tfpw = sum(!filtered_tfpw & !is.finite(p_value_tfpw), na.rm = TRUE),
    n_vc_corrected = sum(vc_corrected %in% TRUE, na.rm = TRUE),
    n_tfpw_applied = sum(tfpw_applied %in% TRUE, na.rm = TRUE)
  ), by = .(variable, period, is_basin_average)]
  fwrite(diag_summary, file.path(out_dir, "root_cause_summary.csv"))
  invisible(diag_summary)
}
####################################################################################
# -- DATA LOADING & PREPROCESSING
####################################################################################
log_event("Loading Precipitation, PET (PM and ERA5-Land PEV), AET, SWE and Temperature data...")

precip_full <- rast(precip_path)
pet_full    <- rast(pet_path)

if (!file.exists(pet_era5_path))
  stop("ERA5-Land PEV PET NetCDF not found: ", pet_era5_path,
       "\nRun 2b_PET_ERALand.R first.")
pet_era5_full <- rast(pet_era5_path)
log_event(sprintf("  [PET_ERA5] Loaded: %d layers from %s", nlyr(pet_era5_full), basename(pet_era5_path)))

# ERA5-Land total evaporation (e) -> AET
if (!file.exists(aet_path)) {
  stop(
    "ERA5-Land total evaporation/AET NetCDF not found. Tried:\n  ",
    paste(AET_NC_CANDIDATES, collapse = "\n  ")
  )
}
aet_full <- rast(aet_path)
log_event(sprintf("  [AET] Loaded: %d layers from %s", nlyr(aet_full), basename(aet_path)))

# ERA5-Land snow water equivalent (sd). Units are read from the file (v9).
swe_full <- NULL
swe_qa   <- list(file = swe_path, run = RUN_SWE)
if (RUN_SWE) {
  swe_full <- rast(swe_path)
  log_event(sprintf("  [SWE] Loaded: %d layers from %s", nlyr(swe_full), basename(swe_path)))
  .u <- unique(terra::units(swe_full)); .u <- .u[nzchar(.u)]
  swe_qa$units_attribute <- if (length(.u)) paste(.u, collapse = "; ") else NA_character_
  .detected <- if (!length(.u)) NA_character_ else
    if (any(grepl("^mm( |$)", .u))) "mm" else
      if (any(grepl("^m( |$)", .u))) "m" else "unknown"
  if (identical(.detected, "unknown"))
    stop("[SWE] Unrecognised units attribute '", swe_qa$units_attribute,
         "'. Expected metres or millimetres of water equivalent.")
  if (is.na(.detected)) {
    log_event(sprintf("  [SWE] WARNING: no units attribute in the file; using configured SWE_UNITS = '%s'.", SWE_UNITS))
  } else if (.detected != SWE_UNITS) {
    log_event(sprintf("  [SWE] WARNING: file declares '%s' but SWE_UNITS = '%s'; using the file's units.",
                      .detected, SWE_UNITS))
    SWE_UNITS <- .detected
  } else {
    log_event(sprintf("  [SWE] Units from file: '%s'.", swe_qa$units_attribute))
  }
  swe_qa$units_used <- SWE_UNITS
} else {
  log_event(paste("  [SWE] NOT FOUND:", swe_path, "- SWE trend/homogeneity analysis will be skipped."))
}

# -- Temperature --
tair_nc <- NULL
for (cand in TAIR_NC_CANDIDATES) {
  if (file.exists(cand)) { tair_nc <- cand; break }
}
if (is.null(tair_nc))
  stop("ERA5-Land 2 m temperature NetCDF not found. Tried:\n  ",
       paste(TAIR_NC_CANDIDATES, collapse="\n  "))
log_event(paste("  [temp] Using:", basename(tair_nc)))
tair_full <- rast(tair_nc)

# -- Date extraction: the precipitation record defines the calendar --
n_layers   <- nlyr(precip_full)
dates_full <- tryCatch({
  d <- as.Date(terra::time(precip_full))
  if (is.null(d) || length(d) != n_layers || all(is.na(d))) NULL else d
}, error=function(e) NULL)
if (is.null(dates_full)) {
  dates_full <- seq(as.Date("1950-01-01"), by="month", length.out=n_layers)
  log_event(paste("Generated", n_layers, "monthly dates from 1950"))
}
dates  <- as.Date(format(dates_full, "%Y-%m-01"))
years  <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))
if (any(diff(years * 12L + months) != 1L))
  stop("The precipitation calendar is not a continuous monthly sequence.")
log_event(sprintf("  [CAL:Precipitation] %d months, %s to %s (reference calendar).",
                  length(dates), format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m")))

# -- Calendar verification for every other input (v9) --
verify_monthly_calendar(pet_full,      dates, "PET_PM")
verify_monthly_calendar(pet_era5_full, dates, "PET_ERA5")
verify_monthly_calendar(aet_full,      dates, "AET")
verify_monthly_calendar(tair_full,     dates, "2-m Temperature")
if (RUN_SWE) verify_monthly_calendar(swe_full, dates, "SWE")

# Temperature uses the precipitation calendar (verified above); one month index for all runs.
dates_t  <- dates
years_t  <- years
months_t <- months

# ------------------------------------------------------------------------------
# PM-PET METEOROLOGICAL CONTROLS
# ------------------------------------------------------------------------------
log_event("Loading PM-PET meteorological-control fields...")
pm_control_full <- lapply(names(PM_CONTROL_NC), function(nm) {
  path <- PM_CONTROL_NC[[nm]]
  if (!file.exists(path)) {
    stop(
      "Required PM-PET meteorological-control NetCDF not found: ", path,
      "\nRun 2b_PET_ERALand.R / the PM-control export step first."
    )
  }
  r <- terra::rast(path)
  log_event(sprintf("  [PM control] %s: %d layers from %s",
                    nm, terra::nlyr(r), basename(path)))
  r
})
names(pm_control_full) <- names(PM_CONTROL_NC)

pm_control_dates <- lapply(names(pm_control_full), function(nm)
  verify_monthly_calendar(pm_control_full[[nm]], dates, paste0("PM control ", nm)))
names(pm_control_dates) <- names(pm_control_full)

# -- Unit conversions / input-unit QA --
# The upstream ERA5-Land precipitation product used by Scripts 1/3 is already
# a mean daily precipitation rate in mm/day. The PM PET product and the
# ERA5-Land native PEV product written by Script 2b are likewise in mm/day.
classify_daily_depth <- function(r, label, expected_min = 0.01, expected_max = 20) {
  july_idx <- which(months == 7L)[1L]
  if (!is.finite(july_idx))
    stop(sprintf("[%s] No July layer available for unit QA.", label))
  
  v <- terra::values(r[[july_idx]], mat = FALSE)
  v <- v[is.finite(v)]
  if (!length(v))
    stop(sprintf("[%s] July layer contains no finite values; cannot validate units.", label))
  
  q <- stats::quantile(v, probs = c(0.01, 0.50, 0.99), na.rm = TRUE, names = FALSE)
  log_event(sprintf("  [UNIT:%s] July q01/median/q99 = %.5g / %.5g / %.5g",
                    label, q[1], q[2], q[3]))
  
  if (q[2] < expected_min) {
    stop(sprintf(
      "[%s] Median July value %.6g is below the expected mm/day range. The file may still be in m/day; do not silently rescale it. Inspect the upstream NetCDF.",
      label, q[2]))
  }
  if (q[2] > expected_max) {
    log_event(sprintf(
      "  [UNIT:%s] July median within the configured QA range (%.3g-%.3g mm/day).",
      label, expected_min, expected_max))
  } else {
    log_event(sprintf(
      "  [UNIT:%s] Magnitude consistent with the configured mm/day QA range.",
      label))
  }
  
  log_event(sprintf("  [UNIT:%s] No 1000x conversion applied; monthly totals will use calendar month length.", label))
  invisible(q)
}

classify_daily_depth(precip_full,  "Precipitation")
classify_daily_depth(pet_full,     "PET_PM")
# expected_max raised to 50 mm/day only so that PEV passes; its magnitude is not physical
classify_daily_depth(pet_era5_full, "PET_ERA5", expected_max = 50)

# ------------------------------------------------------------------------------
# AET unit/sign conversion
#   AET [mm/month] = sign * e [m/day] * 1000 * days_in_month
# The sign convention is detected from the data (ECMWF: evaporation negative).
# ------------------------------------------------------------------------------
aet_july_idx <- which(months == 7L)[1L]
if (is.na(aet_july_idx)) stop("No July layer available for AET unit/sign QA.")
aet_july_vals <- terra::values(aet_full[[aet_july_idx]], mat=FALSE)
aet_july_vals <- aet_july_vals[is.finite(aet_july_vals)]
if (!length(aet_july_vals)) stop("AET July layer contains no finite values.")
aet_q <- stats::quantile(aet_july_vals, c(0.01, 0.50, 0.99), na.rm=TRUE, names=FALSE)
log_event(sprintf("  [UNIT:AET] Native July e q01/median/q99 = %.6g / %.6g / %.6g m/day",
                  aet_q[1], aet_q[2], aet_q[3]))
if (abs(aet_q[2]) > 0.05)
  stop(sprintf("[AET] July median |e| = %.4g looks like mm/day, not m/day. Inspect %s before rescaling.",
               abs(aet_q[2]), basename(aet_path)))

aet_positive_count <- sum(terra::global(aet_full > 0, "sum", na.rm=TRUE)[,1], na.rm=TRUE)
aet_negative_count <- sum(terra::global(aet_full < 0, "sum", na.rm=TRUE)[,1], na.rm=TRUE)
log_event(sprintf("  [AET] Native sign counts: negative=%d, positive=%d",
                  aet_negative_count, aet_positive_count))

aet_sign <- if (aet_negative_count >= aet_positive_count) -1 else 1
if (aet_sign == -1) {
  log_event("  [AET] ECMWF sign convention detected (evaporation negative).")
  if (aet_positive_count > 0L) {
    aet_full <- terra::ifel(aet_full > 0, 0, aet_full)
    log_event(sprintf("  [AET] Set %d positive (condensation) cells to zero before sign conversion.",
                      aet_positive_count))
  }
} else {
  log_event("  [AET] WARNING: evaporation is stored as POSITIVE in this file; no sign flip applied.")
  if (aet_negative_count > 0L) {
    aet_full <- terra::ifel(aet_full < 0, 0, aet_full)
    log_event(sprintf("  [AET] Set %d negative cells to zero.", aet_negative_count))
  }
}

month_days_base  <- c(31,28,31,30,31,30,31,31,30,31,30,31)
is_leap <- function(yr) (yr%%4==0) & (yr%%100!=0 | yr%%400==0)
days_in_month <- month_days_base[months]
days_in_month[is_leap(years) & months==2] <- 29

log_event("Converting Pr, PET (PM), and PET (ERA5) from mm/day to mm/month using actual calendar month length...")
for (i in seq_len(nlyr(precip_full)))    precip_full[[i]]   <- precip_full[[i]]   * days_in_month[i]
for (i in seq_len(nlyr(pet_full)))       pet_full[[i]]      <- pet_full[[i]]      * days_in_month[i]
for (i in seq_len(nlyr(pet_era5_full)))  pet_era5_full[[i]]  <- pet_era5_full[[i]]  * days_in_month[i]
for (i in seq_len(nlyr(aet_full))) {
  aet_full[[i]] <- aet_sign * aet_full[[i]] * 1000 * days_in_month[i]
}

aet_negative_after_conversion <- sum(
  terra::global(aet_full < 0, "sum", na.rm=TRUE)[,1],
  na.rm=TRUE
)
if (aet_negative_after_conversion > 0L) {
  stop("Negative AET remains after ERA5-Land e sign conversion: ",
       aet_negative_after_conversion, " cells.")
}
log_event("  [AET] Converted native e (m/day) to positive AET (mm/month).")

# Temperature: Kelvin -> degC
tair_sample <- as.numeric(terra::global(tair_full[[1]], "mean", na.rm=TRUE))
if (!is.na(tair_sample) && tair_sample > 200) {
  log_event("  [temp] Converting K -> degC...")
  tair_full <- tair_full - 273.15
}
log_event(sprintf("  [temp] Date range: %s to %s (precipitation calendar)", min(dates_t), max(dates_t)))

####################################################################################
# NATIVE ERA5-LAND GRID + BASIN CLIPPING
####################################################################################
# Spatial contract:
#   1. The monthly ERA5-Land precipitation NetCDF defines the native analysis grid.
#   2. All other inputs must already match that native precipitation geometry.
#   3. No raster is reprojected, resampled, aligned, or interpolated.
#   4. Only the basin vector is transformed to the native raster CRS.
#   5. The basin mask is rasterized directly onto the native precipitation grid.
####################################################################################
log_event("Using the native ERA5-Land precipitation grid as the sole spatial template...")

analysis_template <- precip_full[[1]]

# Origins are compared modulo one cell (+0.05 and -0.05 describe the same lattice when a
# grid edge lies exactly on a half-cell coordinate such as -127.85). If the only difference
# is floating-point noise in the stored corners, the template extent is copied onto x
# (a relabelling of identical cells, not a resampling). Returns x invisibly.
assert_native_grid <- function(x, template, label, snap_tol = 1e-6) {
  if (!inherits(x, "SpatRaster")) stop(label, " is not a SpatRaster.")
  if (terra::nlyr(x) < 1L) stop(label, " contains no raster layers.")
  
  rs <- terra::res(template)
  ex <- as.vector(terra::ext(x)); et <- as.vector(terra::ext(template))
  dims_ok <- terra::nrow(x) == terra::nrow(template) && terra::ncol(x) == terra::ncol(template)
  res_ok  <- all(abs(terra::res(x) - rs) / rs < snap_tol)
  ext_dev <- max(abs(ex - et) / rep(rs, each = 2))
  crs_ok  <- isTRUE(tryCatch(terra::same.crs(x, template), error = function(e) FALSE))
  d_org   <- (terra::origin(x) - terra::origin(template)) / rs
  org_ok  <- all(abs(d_org - round(d_org)) < snap_tol)
  
  if (dims_ok && res_ok && crs_ok && org_ok && ext_dev < snap_tol) {
    if (ext_dev > 0) terra::ext(x) <- terra::ext(template)
    return(invisible(x))
  }
  
  stop(sprintf(paste0(
    "%s does not match the native ERA5-Land precipitation grid.\n",
    "  rows x cols : %d x %d (template %d x %d)\n",
    "  resolution  : %.12g x %.12g (template %.12g x %.12g)\n",
    "  extent      : %s (template %s); max deviation %.3g cells\n",
    "  origin      : %s (template %s)\n",
    "  same CRS    : %s\n",
    "Script 4 does not reproject or resample; fix the upstream NetCDF."),
    label, terra::nrow(x), terra::ncol(x), terra::nrow(template), terra::ncol(template),
    terra::res(x)[1], terra::res(x)[2], rs[1], rs[2],
    paste(signif(ex, 12), collapse = ", "), paste(signif(et, 12), collapse = ", "), ext_dev,
    paste(signif(terra::origin(x), 12), collapse = ", "),
    paste(signif(terra::origin(template), 12), collapse = ", "), crs_ok))
}

precip_full   <- assert_native_grid(precip_full,   analysis_template, "Precipitation")
pet_full      <- assert_native_grid(pet_full,      analysis_template, "FAO-56 PM PET")
pet_era5_full <- assert_native_grid(pet_era5_full, analysis_template, "ERA5-Land native PEV")
aet_full      <- assert_native_grid(aet_full,      analysis_template, "ERA5-Land total evaporation / AET")
tair_full     <- assert_native_grid(tair_full,     analysis_template, "2-m Temperature")
if (RUN_SWE) swe_full <- assert_native_grid(swe_full, analysis_template, "ERA5-Land SWE")

for (nm in names(pm_control_full)) {
  pm_control_full[[nm]] <- assert_native_grid(
    pm_control_full[[nm]], analysis_template,
    paste0("PM-PET meteorological control: ", nm))
}

if (is.na(terra::crs(analysis_template)) ||
    !nzchar(terra::crs(analysis_template))) {
  stop("Native ERA5-Land precipitation grid has no CRS.")
}

log_event(sprintf(
  "Native ERA5-Land grid: %d x %d = %d cells | res %.10f x %.10f | CRS=%s",
  terra::nrow(analysis_template), terra::ncol(analysis_template),
  terra::ncell(analysis_template),
  terra::res(analysis_template)[1], terra::res(analysis_template)[2],
  terra::crs(analysis_template)
))

if (!terra::same.crs(terra::vect(basin_boundary), analysis_template)) {
  basin_boundary <- sf::st_transform(
    basin_boundary,
    sf::st_crs(terra::crs(analysis_template))
  )
}

log_event("Rasterizing the authoritative Nechako basin boundary on the native ERA5-Land grid...")

basin_coverage <- terra::rasterize(
  basin_boundary,
  analysis_template,
  cover = TRUE
)

basin_coverage_values <- as.numeric(
  terra::values(basin_coverage, mat = FALSE)
)
basin_coverage_values[!is.finite(basin_coverage_values)] <- 0
basin_coverage_values <- pmax(0, pmin(1, basin_coverage_values))

mask_r <- terra::setValues(
  analysis_template,
  as.numeric(basin_coverage_values > 0)
)

mask_geom_ok <- terra::compareGeom(
  mask_r,
  analysis_template,
  stopOnError = FALSE,
  crs = TRUE,
  ext = TRUE,
  rowcol = TRUE,
  res = TRUE
)
if (!mask_geom_ok) {
  stop("Native ERA5-Land basin mask geometry does not match the precipitation grid.")
}

precip_clipped   <- terra::mask(precip_full,   mask_r, maskvalues = 0, updatevalue = NA)
pet_clipped      <- terra::mask(pet_full,      mask_r, maskvalues = 0, updatevalue = NA)
pet_era5_clipped <- terra::mask(pet_era5_full, mask_r, maskvalues = 0, updatevalue = NA)
aet_clipped      <- terra::mask(aet_full,      mask_r, maskvalues = 0, updatevalue = NA)
tair_clipped     <- terra::mask(tair_full,     mask_r, maskvalues = 0, updatevalue = NA)

if (RUN_SWE) {
  swe_clipped <- terra::mask(swe_full, mask_r, maskvalues = 0, updatevalue = NA)
  if (SWE_UNITS == "m") swe_clipped <- swe_clipped * 1000          # m w.e. -> mm w.e.
  swe_clipped <- clamp_negatives(swe_clipped, "SWE", "mm w.e.")
  swe_qa$n_negative_clamped <- NEG_CLAMP_QA$SWE$n_negative
  swe_qa$most_negative_mm   <- NEG_CLAMP_QA$SWE$most_negative
  .swe_basin_cells <- terra::values(mask_r, mat = FALSE) > 0
  log_event(sprintf("  [SWE] Clipped to basin; %d of %d basin cells have SWE data (units: mm w.e.); %d negative values clamped to 0",
                    sum(!is.na(terra::values(swe_clipped[[1]], mat = FALSE)) & .swe_basin_cells),
                    sum(.swe_basin_cells), swe_qa$n_negative_clamped))
  
  # Magnitude QA (v9): basin-cell median of the March climatology, in mm w.e.
  .mar <- which(months == 3L)
  .v <- terra::values(terra::app(swe_clipped[[.mar]], mean, na.rm = TRUE), mat = FALSE)
  .v <- .v[is.finite(.v)]
  if (!length(.v)) stop("[SWE] March SWE is empty over the basin.")
  swe_qa$march_median_mm <- stats::median(.v)
  swe_qa$march_p99_mm    <- unname(stats::quantile(.v, 0.99))
  log_event(sprintf("  [SWE] March climatology over basin cells: median %.1f mm, 99th percentile %.1f mm",
                    swe_qa$march_median_mm, swe_qa$march_p99_mm))
  if (swe_qa$march_median_mm < 1)
    stop("[SWE] Median March SWE < 1 mm w.e.: the file is probably already in mm while metres were assumed, or the field is empty.")
  if (swe_qa$march_median_mm > 5000)
    stop("[SWE] Median March SWE > 5000 mm w.e.: values were probably scaled twice (mm treated as m).")
}

pm_control_clipped <- lapply(pm_control_full, function(r) {
  terra::mask(r, mask_r, maskvalues = 0, updatevalue = NA)
})
names(pm_control_clipped) <- names(pm_control_full)

# ------------------------------------------------------------------------------
# Additional variables: load, check calendar + grid, clip to basin, convert units
# ------------------------------------------------------------------------------
extra_clipped <- list()
for (nm in names(EXTRA_VARS)) {
  cfg <- EXTRA_VARS[[nm]]
  r <- tryCatch(terra::rast(cfg$path, subds = cfg$nc_var), error = function(e) terra::rast(cfg$path))
  verify_monthly_calendar(r, dates, nm)
  r <- assert_native_grid(r, analysis_template, paste("ERA5-Land", nm))
  
  u <- unique(terra::units(r)); u <- u[nzchar(u)]
  log_event(sprintf("  [%s] %d layers | native units: %s", nm, terra::nlyr(r),
                    if (length(u)) paste(u, collapse = "; ") else "not stored"))
  if (cfg$conv != "none" && length(u) && !any(grepl("^m( |$)", u)))
    stop(sprintf("[%s] expected units in metres, found '%s'; refusing to rescale.", nm, paste(u, collapse = "; ")))
  
  r <- terra::mask(r, mask_r, maskvalues = 0, updatevalue = NA)
  if (cfg$conv == "accum_m_per_day") {
    for (i in seq_len(terra::nlyr(r))) r[[i]] <- r[[i]] * (1000 * days_in_month[i])
  } else if (cfg$conv == "m_to_cm") {
    r <- r * 100
  }
  if (cfg$conv != "none") r <- clamp_negatives(r, nm, cfg$unit)
  log_event(sprintf("  [%s] Clipped to basin and converted to %s", nm, cfg$unit))
  extra_clipped[[nm]] <- r
}

# Snow density is undefined without snow.
if ("SnowDensity" %in% names(extra_clipped)) {
  if ("SnowDepth" %in% names(extra_clipped)) {
    extra_clipped$SnowDensity <- terra::ifel(extra_clipped$SnowDepth < SNOW_PRESENT_DEPTH_CM, NA,
                                             extra_clipped$SnowDensity)
    log_event(sprintf("  [SnowDensity] masked where snow depth < %g cm", SNOW_PRESENT_DEPTH_CM))
  } else if (RUN_SWE) {
    extra_clipped$SnowDensity <- terra::ifel(swe_clipped < 1, NA, extra_clipped$SnowDensity)
    log_event("  [SnowDensity] masked where SWE < 1 mm (snow-depth file not available)")
  } else {
    log_event("  [SnowDensity] WARNING: nothing to mask snow-free cells with; values include snow-free months")
  }
}
fwrite(rbindlist(NEG_CLAMP_QA), file.path(out_dir, "negative_clamp_qa.csv"))
# ------------------------------------------------------------------------------
# SWE input QA, cell level (v9): spatial pattern vs snow depth (catches a flipped or
# shifted file) and implied bulk density = 100 * SWE[mm] / depth[cm] (catches unit errors).
# ------------------------------------------------------------------------------
if (RUN_SWE && "SnowDepth" %in% names(extra_clipped)) {
  .mar <- which(months == 3L)
  .dep <- terra::values(terra::app(extra_clipped$SnowDepth[[.mar]], mean, na.rm = TRUE), mat = FALSE)  # cm
  .swe <- terra::values(terra::app(swe_clipped[[.mar]], mean, na.rm = TRUE), mat = FALSE)              # mm
  .ok  <- is.finite(.dep) & is.finite(.swe)
  swe_qa$depth_swe_march_spatial_r <- if (sum(.ok) > 20L) stats::cor(.dep[.ok], .swe[.ok]) else NA_real_
  .deep <- .ok & .dep >= 10
  swe_qa$implied_density_march_median <- if (sum(.deep) >= 10L) stats::median(100 * .swe[.deep] / .dep[.deep]) else NA_real_
  log_event(sprintf("  [SWE QA] March spatial correlation with snow depth r = %.3f | implied bulk density (cells with depth >= 10 cm) = %.0f kg m-3",
                    swe_qa$depth_swe_march_spatial_r, swe_qa$implied_density_march_median))
  if (!is.finite(swe_qa$depth_swe_march_spatial_r) || swe_qa$depth_swe_march_spatial_r < 0.7)
    stop("[SWE QA] SWE and snow-depth March patterns disagree (r = ", signif(swe_qa$depth_swe_march_spatial_r, 3),
         "); one of the files may be flipped or shifted.")
  .d <- swe_qa$implied_density_march_median
  if (is.finite(.d) && (.d < 30 || .d > 1000))
    stop("[SWE QA] Implied bulk density ", round(.d), " kg m-3 is physically impossible: SWE or snow-depth units are wrong.")
  if (is.finite(.d) && (.d < 100 || .d > 550))
    log_event("  [SWE QA] WARNING: implied bulk density outside the usual 100-550 kg m-3 range; check SWE and snow-depth units.")
  if ("SnowDensity" %in% names(extra_clipped)) {
    .den <- terra::values(terra::app(extra_clipped$SnowDensity[[.mar]], mean, na.rm = TRUE), mat = FALSE)
    .sel <- .deep & is.finite(.den)
    swe_qa$file_density_march_median  <- if (sum(.sel) >= 10L) stats::median(.den[.sel]) else NA_real_
    swe_qa$implied_over_file_density  <- .d / swe_qa$file_density_march_median
    log_event(sprintf("  [SWE QA] Snow-density file (same cells) median %.0f kg m-3; implied / file = %.2f",
                      swe_qa$file_density_march_median, swe_qa$implied_over_file_density))
    if (is.finite(swe_qa$implied_over_file_density) &&
        (swe_qa$implied_over_file_density < 0.7 || swe_qa$implied_over_file_density > 1.3))
      log_event("  [SWE QA] WARNING: implied and stored snow density differ by more than 30%.")
  }
}

# ------------------------------------------------------------------------------
# PET zero floor (v9: BEFORE any averaging, so basin, regional, pixel and point PET
# are the same series). Exact zeros become min_positive_value.
# ------------------------------------------------------------------------------
floor_pet_zeros <- function(r, label) {
  nz <- sum(terra::global(r == 0, "sum", na.rm = TRUE)[, 1], na.rm = TRUE)
  nm <- names(r)
  r  <- terra::ifel(r == 0, min_positive_value, r)
  names(r) <- nm
  log_event(sprintf("  Replaced %d PET zeros with %.3f mm (%s), before basin averaging",
                    nz, min_positive_value, label))
  r
}
pet_clipped      <- floor_pet_zeros(pet_clipped,      "PET_PM")
pet_era5_clipped <- floor_pet_zeros(pet_era5_clipped, "PET_ERA5")

clipped_template <- terra::rast(analysis_template, nlyrs = 1)
terra::values(clipped_template) <- NA_real_
saveRDS(clipped_template, file.path(out_dir, "clipped_template.rds"))

analysis_mask <- as.logical(
  terra::values(mask_r, mat = FALSE)
)

n_bbox_cells   <- terra::ncell(analysis_template)
valid_mask     <- !is.na(terra::values(precip_clipped[[1]], mat = FALSE))
n_basin_pixels <- sum(valid_mask)
reduction_pct  <- 100 * (1 - n_basin_pixels / n_bbox_cells)

if (length(analysis_mask) != n_bbox_cells) {
  stop("Native basin mask length does not match native grid cell count.")
}

if (!identical(as.integer(n_basin_pixels),
               as.integer(sum(analysis_mask)))) {
  stop("Native basin mask/data mismatch: valid data cells do not equal native basin pixels.")
}

native_mask <- list(
  template = clipped_template,
  basin_mask = analysis_mask,
  basin_pixels = n_basin_pixels,
  ncell = n_bbox_cells,
  crs = terra::crs(analysis_template),
  extent = as.vector(terra::ext(analysis_template)),
  resolution = terra::res(analysis_template),
  source = "Native ERA5-Land precipitation grid",
  created = Sys.time()
)
saveRDS(native_mask,
        file.path(out_dir, "native_basin_mask.rds"),
        compress = TRUE)

log_event(sprintf(
  "NATIVE ERA5-LAND BASIN GRID: %d cells | %d basin pixels (%.1f%% outside mask)",
  n_bbox_cells, n_basin_pixels, reduction_pct
))

if (n_basin_pixels == 0L) {
  stop("CRITICAL: No valid cells after native-grid basin clipping.")
}

####################################################################################
# -- BASIN AVERAGE TIME SERIES
####################################################################################
log_event(sprintf("Computing basin-averaged monthly time series (7-subbasin union; fractional-coverage area weights; min valid coverage %.2f)...", BASIN_MIN_VALID_AREA_COVERAGE))
whole_region_weights <- region_fractional_weights(analysis_template,
                                                  load_nechako_regions(template=analysis_template))[[NECHAKO_REGION_BASIN_ID]]$weights
whole_region_weights <- whole_region_weights[valid_mask]
precip_monthly_avg    <- regional_weighted_mean(terra::values(precip_clipped)[valid_mask,,drop=FALSE], whole_region_weights, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE)
pet_monthly_avg       <- regional_weighted_mean(terra::values(pet_clipped)[valid_mask,,drop=FALSE], whole_region_weights, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE)
pet_era5_monthly_avg  <- regional_weighted_mean(terra::values(pet_era5_clipped)[valid_mask,,drop=FALSE], whole_region_weights, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE)
aet_monthly_avg       <- regional_weighted_mean(terra::values(aet_clipped)[valid_mask,,drop=FALSE], whole_region_weights, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE)
tair_monthly_avg      <- regional_weighted_mean(terra::values(tair_clipped)[valid_mask,,drop=FALSE], whole_region_weights, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE)

pm_control_monthly_avg <- lapply(pm_control_clipped, function(r) {
  regional_weighted_mean(
    terra::values(r)[valid_mask, , drop = FALSE],
    whole_region_weights,
    min_coverage = BASIN_MIN_VALID_AREA_COVERAGE
  )
})
names(pm_control_monthly_avg) <- names(pm_control_clipped)

# Annual aggregation
precip_annual_avg_matrix   <- aggregate_to_annual_fast(
  matrix(precip_monthly_avg, ncol=1), years, method="sum")
pet_annual_avg_matrix      <- aggregate_to_annual_warmseason_fast(
  matrix(pet_monthly_avg, ncol=1), years, months_vec=months,
  warm_months=warm_months_pet, method="sum")
pet_era5_annual_avg_matrix  <- aggregate_to_annual_warmseason_fast(
  matrix(pet_era5_monthly_avg, ncol=1), years, months_vec=months,
  warm_months=warm_months_pet, method="sum")
tair_annual_avg_matrix     <- aggregate_to_annual_fast(
  matrix(tair_monthly_avg, ncol=1), years_t, method="mean")
aet_annual_avg_matrix       <- aggregate_to_annual_fast(
  matrix(aet_monthly_avg, ncol=1), years, method="sum")

pm_control_annual_avg_matrix <- lapply(pm_control_monthly_avg, function(x) {
  aggregate_to_annual_fast(matrix(x, ncol = 1), years, method = "mean")
})
names(pm_control_annual_avg_matrix) <- names(pm_control_monthly_avg)

precip_annual_avg    <- precip_annual_avg_matrix[,1]
pet_annual_avg       <- pet_annual_avg_matrix[,1]
pet_era5_annual_avg   <- pet_era5_annual_avg_matrix[,1]
tair_annual_avg      <- tair_annual_avg_matrix[,1]
aet_annual_avg       <- aet_annual_avg_matrix[,1]
pm_control_annual_avg <- lapply(pm_control_annual_avg_matrix, function(m) m[,1])

annual_years         <- as.integer(rownames(precip_annual_avg_matrix))
tair_annual_years    <- as.integer(rownames(tair_annual_avg_matrix))

if (RUN_SWE) {
  swe_monthly_avg  <- regional_weighted_mean(terra::values(swe_clipped)[valid_mask,,drop=FALSE], whole_region_weights, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE)
  swe_marapr_basin <- annual_metric_fast(swe_monthly_avg, years, months, SWE_ANNUAL_MONTHS, "mean")[,1]
  swe_peak_basin   <- annual_metric_fast(swe_monthly_avg, years, months, SWE_PEAK_MONTHS,  "max")[,1]
} else {
  swe_monthly_avg  <- rep(NA_real_, length(dates))
  swe_marapr_basin <- swe_peak_basin <- rep(NA_real_, length(annual_years))
}

.extra_annual <- function(mat, cfg) {
  if (!is.matrix(mat)) mat <- matrix(as.numeric(mat), ncol = 1L)
  if (identical(cfg$annual, "sum")) aggregate_to_annual_fast(mat, years, method = "sum")
  else annual_metric_fast(mat, years, months, SWE_ANNUAL_MONTHS, "mean")
}
extra_monthly_avg <- lapply(extra_clipped, function(r)
  regional_weighted_mean(terra::values(r)[valid_mask, , drop = FALSE], whole_region_weights, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE))
extra_annual_avg <- lapply(names(extra_monthly_avg), function(nm)
  .extra_annual(extra_monthly_avg[[nm]], EXTRA_VARS[[nm]])[, 1])
names(extra_annual_avg) <- names(extra_monthly_avg)

if ("Snowfall" %in% names(extra_annual_avg)) {
  .r_sf <- mean(extra_annual_avg$Snowfall, na.rm = TRUE) / mean(precip_annual_avg, na.rm = TRUE)
  log_event(sprintf("  [QA] basin snowfall / precipitation = %.2f", .r_sf))
  if (!is.finite(.r_sf) || .r_sf > 1.05 || .r_sf < 0.01)
    stop("Snowfall/precipitation ratio is implausible (", signif(.r_sf, 3), "): check units of both files.")
}
if ("Runoff" %in% names(extra_annual_avg)) {
  .r_ro <- mean(extra_annual_avg$Runoff, na.rm = TRUE) / mean(precip_annual_avg, na.rm = TRUE)
  log_event(sprintf("  [QA] basin runoff / precipitation = %.2f", .r_ro))
  if (!is.finite(.r_ro) || .r_ro > 1.5 || .r_ro < 0.01)
    warning("Runoff/precipitation ratio ", signif(.r_ro, 3), " is outside 0.01-1.5: check units.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# SWE input QA, basin level (v9): peak SWE cannot greatly exceed the cold-season water input.
# ------------------------------------------------------------------------------
if (RUN_SWE) {
  .cold <- c(10L, 11L, 12L, 1L, 2L, 3L, 4L)
  .p_cold <- season_by_year(precip_monthly_avg, years, months, .cold, wrap_from = 10L, fun = "sum")
  swe_qa$peak_swe_over_OctApr_precip <- mean(swe_peak_basin, na.rm = TRUE) / mean(.p_cold, na.rm = TRUE)
  log_event(sprintf("  [SWE QA] basin mean peak SWE / mean Oct-Apr precipitation = %.2f",
                    swe_qa$peak_swe_over_OctApr_precip))
  if (!is.finite(swe_qa$peak_swe_over_OctApr_precip) ||
      swe_qa$peak_swe_over_OctApr_precip > 2 || swe_qa$peak_swe_over_OctApr_precip < 0.02)
    stop("[SWE QA] Peak SWE / Oct-Apr precipitation = ", signif(swe_qa$peak_swe_over_OctApr_precip, 3),
         " is impossible: SWE or precipitation units are wrong.")
  if (swe_qa$peak_swe_over_OctApr_precip > 1.2)
    log_event("  [SWE QA] WARNING: peak SWE exceeds 120% of Oct-Apr precipitation.")
  if ("Snowfall" %in% names(extra_monthly_avg)) {
    .s_cold <- season_by_year(extra_monthly_avg$Snowfall, years, months, .cold, wrap_from = 10L, fun = "sum")
    swe_qa$peak_swe_over_OctApr_snowfall <- mean(swe_peak_basin, na.rm = TRUE) / mean(.s_cold, na.rm = TRUE)
    log_event(sprintf("  [SWE QA] basin mean peak SWE / mean Oct-Apr snowfall = %.2f",
                      swe_qa$peak_swe_over_OctApr_snowfall))
    if (!is.finite(swe_qa$peak_swe_over_OctApr_snowfall) ||
        swe_qa$peak_swe_over_OctApr_snowfall > 3 || swe_qa$peak_swe_over_OctApr_snowfall < 0.02)
      stop("[SWE QA] Peak SWE / Oct-Apr snowfall = ", signif(swe_qa$peak_swe_over_OctApr_snowfall, 3),
           " is impossible: SWE or snowfall units are wrong.")
    if (swe_qa$peak_swe_over_OctApr_snowfall > 1.5)
      log_event("  [SWE QA] WARNING: peak SWE exceeds 150% of Oct-Apr snowfall.")
  }
}
swe_qa_dt <- as.data.table(lapply(swe_qa, function(v) if (is.null(v) || !length(v)) NA else v[1]))
fwrite(swe_qa_dt, file.path(out_dir, "swe_input_qa.csv"))

# Pixel matrices (time x pixels)
valid_cell_indices <- which(valid_mask)
valid_xy           <- xyFromCell(precip_clipped, valid_cell_indices)
precip_matrix      <- t(values(precip_clipped)[valid_mask,    , drop=FALSE])
pet_matrix         <- t(values(pet_clipped)[valid_mask,       , drop=FALSE])
pet_era5_matrix    <- t(values(pet_era5_clipped)[valid_mask,  , drop=FALSE])
aet_matrix         <- t(values(aet_clipped)[valid_mask,       , drop=FALSE])
tair_matrix        <- t(values(tair_clipped)[valid_mask,      , drop=FALSE])
aet_annual_pix     <- aggregate_to_annual_fast(aet_matrix, years, method = "sum")   # v9: AET trends
swe_marapr_pix <- swe_peak_pix <- NULL
if (RUN_SWE) {
  swe_matrix     <- t(values(swe_clipped)[valid_mask, , drop=FALSE])
  swe_marapr_pix <- annual_metric_fast(swe_matrix, years, months, SWE_ANNUAL_MONTHS, "mean")
  swe_peak_pix   <- annual_metric_fast(swe_matrix, years, months, SWE_PEAK_MONTHS,  "max")
  if (!identical(rownames(swe_marapr_pix), as.character(annual_years)))
    stop("SWE annual metric years do not match the precipitation annual years.")
  log_event(sprintf("  [SWE matrix] %d x %d | annual Mar-Apr matrix %d x %d",
                    nrow(swe_matrix), ncol(swe_matrix), nrow(swe_marapr_pix), ncol(swe_marapr_pix)))
}

extra_matrix <- lapply(extra_clipped, function(r) t(values(r)[valid_mask, , drop = FALSE]))
extra_annual_pix <- lapply(names(extra_matrix), function(nm) {
  a <- .extra_annual(extra_matrix[[nm]], EXTRA_VARS[[nm]])
  if (!identical(rownames(a), as.character(annual_years)))
    stop(sprintf("[%s] annual years do not match the precipitation annual years.", nm))
  a
})
names(extra_annual_pix) <- names(extra_matrix)
for (nm in names(extra_matrix))
  log_event(sprintf("  [%s matrix] %d x %d | annual matrix %d x %d", nm,
                    nrow(extra_matrix[[nm]]), ncol(extra_matrix[[nm]]),
                    nrow(extra_annual_pix[[nm]]), ncol(extra_annual_pix[[nm]])))

pm_control_matrices <- lapply(pm_control_clipped, function(r) {
  t(values(r)[valid_mask, , drop = FALSE])
})
names(pm_control_matrices) <- names(pm_control_clipped)
.pm_t_maxdiff <- max(abs(pm_control_matrices$Temperature - tair_matrix), na.rm = TRUE)
PM_T_IS_TAIR  <- is.finite(.pm_t_maxdiff) && .pm_t_maxdiff < 1e-3
log_event(sprintf("  [PM control] Temperature vs 2-m temperature input: max |difference| = %.3g degC -> %s",
                  .pm_t_maxdiff, if (PM_T_IS_TAIR) "same field; PM_Temperature results are copied from Temperature"
                  else "different fields; both are computed"))
.relabel <- function(dt, new) {
  dt <- data.table::copy(dt); dt[, variable := new]
  if ("variable_name" %in% names(dt)) dt[, variable_name := new]
  dt
}
coords_dt <- data.table(
  space_idx=seq_len(n_basin_pixels), x=valid_xy[,1], y=valid_xy[,2],
  is_basin_average=FALSE)
log_event(sprintf("  Matrices: Pr/PET %d x %d | PET_ERA5 %d x %d | AET %d x %d | Tair %d x %d",
                  nrow(precip_matrix), ncol(precip_matrix),
                  nrow(pet_era5_matrix), ncol(pet_era5_matrix),
                  nrow(aet_matrix), ncol(aet_matrix),
                  nrow(tair_matrix), ncol(tair_matrix)))
for (nm in names(pm_control_matrices)) {
  mm <- pm_control_matrices[[nm]]
  log_event(sprintf("  [PM control matrix] %s: %d x %d", nm, nrow(mm), ncol(mm)))
}

# One month index for every pixel and basin run (v9).
month_index_list   <- compute_month_index(months)
month_index_list_t <- month_index_list

####################################################################################
# -- HELPER: BUILD RESULT DATA.TABLE ROW SET (annual + 12 monthly)
####################################################################################
unpack_results <- function(res_list, which = c("vc","tf")) {
  which <- match.arg(which)
  r <- lapply(res_list, `[[`, which)
  data.table(
    tau       = vapply(r, `[[`, numeric(1),   "tau"),
    p         = vapply(r, `[[`, numeric(1),   "p"),
    sl        = vapply(r, `[[`, numeric(1),   "sl"),
    S         = vapply(r, `[[`, numeric(1),   "S"),
    varS      = vapply(r, `[[`, numeric(1),   "varS"),
    n         = vapply(r, `[[`, numeric(1),   "n"),
    rho1      = vapply(r, `[[`, numeric(1),   "rho1"),
    n_ties    = vapply(r, `[[`, numeric(1),   "n_ties"),
    pct_ties  = vapply(r, `[[`, numeric(1),   "percent_ties"),
    n_min     = vapply(r, `[[`, numeric(1),   "n_min_vals"),
    pct_min   = vapply(r, `[[`, numeric(1),   "percent_min_vals"),
    tau_b_adj = vapply(r, `[[`, logical(1),   "tau_b_adjusted"),
    filtered  = vapply(r, `[[`, logical(1),   "filtered"),
    reason    = vapply(r, `[[`, character(1), "reason"),
    extra     = if (which=="vc")
      vapply(r, function(z) z$vc_corrected,  logical(1)) else
        vapply(r, function(z) z$tfpw_applied,  logical(1))
  )
}

build_annual_dt <- function(res_annual, coords_dt, var_name) {
  vc_a   <- unpack_results(res_annual, "vc")
  tf_a   <- unpack_results(res_annual, "tf")
  spec_a <- data.table(
    n_spectral_peaks    = vapply(lapply(res_annual,`[[`,"spec"), `[[`, integer(1), "n_peaks"),
    dominant_period     = vapply(lapply(res_annual,`[[`,"spec"), `[[`, numeric(1), "dominant_period"),
    spectral_confidence = vapply(lapply(res_annual,`[[`,"spec"), `[[`, numeric(1), "conf"))
  cpt_a  <- data.table(
    changepoint_detected   = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, logical(1),   "changepoint_detected"),
    changepoint_position   = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, integer(1),   "changepoint_position"),
    n_changepoints         = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, integer(1),   "n_changepoints"),
    first_changepoint_year = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, integer(1),   "first_changepoint_year"),
    mean_before_shift      = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, numeric(1),   "mean_before"),
    mean_after_shift       = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, numeric(1),   "mean_after"),
    magnitude_shift        = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, numeric(1),   "magnitude_shift"))
  cbind(coords_dt,
        variable="annual", period="annual", month=NA_integer_, variable_name=var_name,
        data.table(tau_vc=vc_a$tau, p_value_vc=vc_a$p, sl_vc=vc_a$sl,
                   vc_corrected=vc_a$extra, n_ties_vc=vc_a$n_ties,
                   percent_ties_vc=vc_a$pct_ties, n_min_vals_vc=vc_a$n_min,
                   percent_min_vals_vc=vc_a$pct_min, tau_b_adjusted_vc=vc_a$tau_b_adj,
                   filtered_vc=vc_a$filtered, filter_reason_vc=vc_a$reason),
        data.table(tau_tfpw=tf_a$tau, p_value_tfpw=tf_a$p, sl_tfpw=tf_a$sl,
                   tfpw_applied=tf_a$extra, n_ties_tfpw=tf_a$n_ties,
                   percent_ties_tfpw=tf_a$pct_ties, n_min_vals_tfpw=tf_a$n_min,
                   percent_min_vals_tfpw=tf_a$pct_min, tau_b_adjusted_tfpw=tf_a$tau_b_adj,
                   filtered_tfpw=tf_a$filtered, filter_reason_tfpw=tf_a$reason),
        n=vc_a$n, rho1_vc=vc_a$rho1, rho1_tfpw=vc_a$rho1,
        spec_a, cpt_a,
        same_significance=(vc_a$p < alpha) == (tf_a$p < alpha),
        same_direction   =sign(vc_a$tau) == sign(tf_a$tau))
}

####################################################################################
# ROBUST TREND-TEST WRAPPER
####################################################################################
.run_modifiedmk_tfpw <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NULL)
  z <- tryCatch(suppressWarnings(modifiedmk::tfpwmk(x)), error = function(e) NULL)
  if (is.null(z)) return(NULL)
  nm <- names(z)
  getv <- function(pattern) {
    k <- grep(pattern, nm, ignore.case = TRUE)
    if (length(k)) unname(z[[k[1]]]) else NA_real_
  }
  tau <- getv("^Tau$")
  p   <- getv("^P-?value$")
  sl  <- getv("^Old Sen")
  if (!is.finite(tau) || !is.finite(p) || !is.finite(sl)) return(NULL)
  ties <- length(x) - length(unique(x))
  list(
    tau = tau, p = p, sl = sl,
    S = getv("^S$"), varS = getv("^Var\\(S\\)$"),
    n = length(x), rho1 = NA_real_, tfpw_applied = TRUE,
    n_ties = ties, percent_ties = 100 * ties / length(x),
    n_min_vals = 0, percent_min_vals = 0, tau_b_adjusted = FALSE,
    filtered = FALSE, reason = "none"
  )
}

# Start-up self-test: fail loudly if the VC or TFPW code paths cannot parse the package
# output or violate the variance-cap rule.
local({
  set.seed(42)
  x_test <- as.numeric(stats::filter(rnorm(76), 0.4, method = "recursive")) +
    0.03 * seq_len(76)
  chk <- mk_tfpw_spectral_for_series(x_test, FALSE, alpha, max_tie_percent, 20L,
                                     new.env(parent = emptyenv()), start_year = 1950L)
  t_shared <- chk$tf
  t_pkg <- .run_modifiedmk_tfpw(x_test)
  mm <- .mmkh_fields(modifiedmk::mmkh(x_test))
  if (isTRUE(chk$vc$filtered) || !is.finite(chk$vc$p) ||
      isTRUE(t_shared$filtered) || !is.finite(t_shared$p) || is.null(t_pkg) ||
      !all(is.finite(unlist(mm[c("tau", "p_new", "p_old", "var_new", "var_old")])))) {
    stop("VC/TFPW self-test failed. VC reason=", chk$vc$reason,
         "; TFPW reason=", t_shared$reason,
         "; names(mmkh) = ", paste(names(modifiedmk::mmkh(x_test)), collapse = " | "),
         "; names(tfpwmk) = ", paste(names(modifiedmk::tfpwmk(x_test)), collapse = " | "))
  }
  x_neg <- as.numeric(stats::filter(rnorm(76), -0.6, method = "recursive")) +
    0.01 * seq_len(76)
  c2 <- mk_tfpw_spectral_for_series(x_neg, FALSE, alpha, max_tie_percent, 20L,
                                    new.env(parent = emptyenv()), start_year = 1950L)
  if (isTRUE(c2$vc$filtered) || (is.finite(c2$vc$p) && is.finite(c2$vc$vc_p_uncorrected) &&
                                 c2$vc$p < c2$vc$vc_p_uncorrected - 1e-9))
    stop("VC self-test: corrected p is smaller than the uncorrected p (cap not applied).")
  set.seed(1)
  x_big <- 3000 + 300 * rnorm(76)
  cp <- detect_regime_shift_pelt(x_big, min_obs = min_obs_changepoint,
                                 start_year = 1950L, minseglen = pelt_minseglen)
  log_event(sprintf(
    "[SELFTEST] VC p=%.4f (n/n*=%.2f, %s) | neg-AR VC p=%.4f vs uncorrected %.4f (%s) | shared TFPW p=%.4f tau=%.3f sl=%.4f | package TFPW QA=%s | PELT on noise (sd=300): %s",
    chk$vc$p, chk$vc$vc_n_over_nstar, chk$vc$reason,
    c2$vc$p, c2$vc$vc_p_uncorrected, c2$vc$reason,
    t_shared$p, t_shared$tau, t_shared$sl,
    if (is.null(t_pkg)) "FAIL" else "available",
    if (cp$changepoint_detected) paste("break", cp$first_changepoint_year) else "no break"))
})

.internal_gap <- function(x) {
  f <- which(is.finite(x))
  length(f) > 0L && any(!is.finite(x[f[1]:f[length(f)]]))
}

run_trend_bundle <- function(x, is_precip, alpha, max_tie_percent,
                             n_sim_spectral, conf_env,
                             start_year = NA_integer_,
                             min_nonzero_count = NULL,
                             skip_min_filter = FALSE) {
  x0  <- as.numeric(x)
  res <- mk_tfpw_spectral_for_series(
    x0, is_precip, alpha, max_tie_percent, n_sim_spectral, conf_env,
    start_year = start_year,
    min_nonzero_count = min_nonzero_count,
    skip_min_filter = skip_min_filter,
    max_min_value_pct_precip = max_min_value_pct_precip,
    max_min_value_pct_pet = max_min_value_pct_pet,
    min_positive_value = min_positive_value,
    min_obs_changepoint = min_obs_changepoint,
    also_check_min_pct = PET_MONTHLY_ALSO_PCT_FILTER,
    pelt_minseglen = pelt_minseglen
  )
  # v10: missing years inside the record -> both tests filtered with the same reason
  # (MK, Sen, TFPW, spectral and PELT all assume an evenly spaced series).
  if (.internal_gap(x0) && !isTRUE(res$tf$filtered)) {
    for (k in c("vc", "tf")) {
      res[[k]]$filtered <- TRUE; res[[k]]$reason <- "internal_gaps"
      res[[k]]$tau <- NA_real_; res[[k]]$p <- NA_real_; res[[k]]$sl <- NA_real_
    }
    res$spec$n_peaks <- 0L; res$spec$dominant_period <- NA_real_; res$spec$conf <- NA_real_
    res$cpt$changepoint_detected <- FALSE; res$cpt$n_changepoints <- 0L
    res$cpt$changepoint_position <- NA_integer_; res$cpt$first_changepoint_year <- NA_integer_
    res$cpt$mean_before <- NA_real_; res$cpt$mean_after <- NA_real_; res$cpt$magnitude_shift <- NA_real_
  }
  res
}

####################################################################################
# -- BASIN AVERAGE TREND PROCESSING
####################################################################################
process_basin_series <- function(ts_monthly, ts_annual, var_name,
                                 is_precip, skip_min_filter = FALSE,
                                 m_index_list = month_index_list,
                                 annual_only = FALSE) {
  conf_env <- new.env(parent=emptyenv())
  log_series_quality(ts_annual, sprintf("%s basin annual", var_name),
                     floor_value = if (is_precip) 0 else if (skip_min_filter) NA_real_ else min_positive_value)
  log_series_quality(ts_monthly, sprintf("%s basin monthly", var_name),
                     floor_value = if (is_precip) 0 else if (skip_min_filter) NA_real_ else min_positive_value)
  res_ann  <- run_trend_bundle(ts_annual, is_precip, alpha, max_tie_percent,
                               n_sim_spectral, conf_env,
                               start_year      = min(annual_years),
                               skip_min_filter = skip_min_filter)
  vc_a  <- res_ann$vc;  tf_a  <- res_ann$tf
  spec_a <- res_ann$spec; cpt_a <- res_ann$cpt
  
  annual_dt <- data.table(
    space_idx=NA_integer_, x=NA_real_, y=NA_real_,
    variable=var_name, period="annual", month=NA_integer_, variable_name=var_name,
    tau_vc=vc_a$tau, p_value_vc=vc_a$p, sl_vc=vc_a$sl, vc_corrected=vc_a$vc_corrected,
    n_ties_vc=vc_a$n_ties, percent_ties_vc=vc_a$percent_ties,
    n_min_vals_vc=vc_a$n_min_vals, percent_min_vals_vc=vc_a$percent_min_vals,
    tau_b_adjusted_vc=vc_a$tau_b_adjusted, filtered_vc=vc_a$filtered,
    filter_reason_vc=vc_a$reason,
    tau_tfpw=tf_a$tau, p_value_tfpw=tf_a$p, sl_tfpw=tf_a$sl, tfpw_applied=tf_a$tfpw_applied,
    n_ties_tfpw=tf_a$n_ties, percent_ties_tfpw=tf_a$percent_ties,
    n_min_vals_tfpw=tf_a$n_min_vals, percent_min_vals_tfpw=tf_a$percent_min_vals,
    tau_b_adjusted_tfpw=tf_a$tau_b_adjusted, filtered_tfpw=tf_a$filtered,
    filter_reason_tfpw=tf_a$reason,
    n=NA_integer_, rho1_vc=NA_real_, rho1_tfpw=NA_real_,
    n_spectral_peaks=spec_a$n_peaks, dominant_period=spec_a$dominant_period,
    spectral_confidence=spec_a$conf,
    changepoint_detected=cpt_a$changepoint_detected,
    changepoint_position=cpt_a$changepoint_position,
    n_changepoints=cpt_a$n_changepoints,
    first_changepoint_year=cpt_a$first_changepoint_year,
    mean_before_shift=cpt_a$mean_before, mean_after_shift=cpt_a$mean_after,
    magnitude_shift=cpt_a$magnitude_shift,
    same_significance=(vc_a$p < alpha)==(tf_a$p < alpha),
    same_direction=sign(vc_a$tau)==sign(tf_a$tau),
    is_basin_average=TRUE)
  
  if (annual_only) return(annual_dt)
  
  monthly_dts <- lapply(1:12, function(m) {
    mi <- m_index_list[[as.character(m)]]
    mn_count <- if (!skip_min_filter && !is_precip) min_nonzero_pet_monthly else NULL
    res_m <- if (length(mi) && !all(is.na(ts_monthly[mi])))
      run_trend_bundle(ts_monthly[mi], is_precip, alpha, max_tie_percent,
                       n_sim_spectral, conf_env,
                       min_nonzero_count = mn_count,
                       skip_min_filter   = skip_min_filter)
    else list(vc=list(tau=NA_real_,p=NA_real_,sl=NA_real_,S=NA_real_,varS=NA_real_,
                      n=0,rho1=NA_real_,vc_corrected=FALSE,n_ties=0,percent_ties=0,
                      n_min_vals=0,percent_min_vals=0,tau_b_adjusted=FALSE,
                      filtered=TRUE,reason="no_data"),
              tf=list(tau=NA_real_,p=NA_real_,sl=NA_real_,S=NA_real_,varS=NA_real_,
                      n=0,rho1=NA_real_,tfpw_applied=FALSE,n_ties=0,percent_ties=0,
                      n_min_vals=0,percent_min_vals=0,tau_b_adjusted=FALSE,
                      filtered=TRUE,reason="no_data"),
              spec=list(n_peaks=0L,dominant_period=NA_real_,conf=NA_real_),
              cpt=list(changepoint_detected=FALSE,changepoint_position=NA_integer_,
                       n_changepoints=0L,first_changepoint_year=NA_integer_,
                       mean_before=NA_real_,mean_after=NA_real_,magnitude_shift=NA_real_))
    vc_m <- res_m$vc; tf_m <- res_m$tf; sp_m <- res_m$spec; cp_m <- res_m$cpt
    data.table(
      space_idx=NA_integer_,x=NA_real_,y=NA_real_,
      variable=var_name, period="monthly", month=m, variable_name=var_name,
      tau_vc=vc_m$tau, p_value_vc=vc_m$p, sl_vc=vc_m$sl, vc_corrected=vc_m$vc_corrected,
      n_ties_vc=vc_m$n_ties, percent_ties_vc=vc_m$percent_ties,
      n_min_vals_vc=vc_m$n_min_vals, percent_min_vals_vc=vc_m$percent_min_vals,
      tau_b_adjusted_vc=vc_m$tau_b_adjusted, filtered_vc=vc_m$filtered,
      filter_reason_vc=vc_m$reason,
      tau_tfpw=tf_m$tau, p_value_tfpw=tf_m$p, sl_tfpw=tf_m$sl, tfpw_applied=tf_m$tfpw_applied,
      n_ties_tfpw=tf_m$n_ties, percent_ties_tfpw=tf_m$percent_ties,
      n_min_vals_tfpw=tf_m$n_min_vals, percent_min_vals_tfpw=tf_m$percent_min_vals,
      tau_b_adjusted_tfpw=tf_m$tau_b_adjusted, filtered_tfpw=tf_m$filtered,
      filter_reason_tfpw=tf_m$reason,
      n=NA_integer_, rho1_vc=NA_real_, rho1_tfpw=NA_real_,
      n_spectral_peaks=sp_m$n_peaks, dominant_period=sp_m$dominant_period,
      spectral_confidence=sp_m$conf,
      changepoint_detected=cp_m$changepoint_detected,
      changepoint_position=cp_m$changepoint_position,
      n_changepoints=cp_m$n_changepoints,
      first_changepoint_year=cp_m$first_changepoint_year,
      mean_before_shift=cp_m$mean_before, mean_after_shift=cp_m$mean_after,
      magnitude_shift=cp_m$magnitude_shift,
      same_significance=(vc_m$p < alpha)==(tf_m$p < alpha),
      same_direction=sign(vc_m$tau)==sign(tf_m$tau),
      is_basin_average=TRUE)
  })
  rbindlist(c(list(annual_dt), monthly_dts))
}

log_event("Processing basin-averaged series...")
log_event("  [DIAG] Basin-average preflight: reporting variance / tie density before trend testing.")
# ------------------------------------------------------------------------------
# PET-source QA (NON-FATAL diagnostic)
# ------------------------------------------------------------------------------
PET_RATIO_REVIEW_RANGE <- c(0.5, 2.0)

.pm_mean_mm  <- mean(pet_monthly_avg,      na.rm = TRUE)
.pev_mean_mm <- mean(pet_era5_monthly_avg, na.rm = TRUE)
.pet_pm_era5_ratio <- .pev_mean_mm / .pm_mean_mm

log_event(sprintf(
  "[QA] Basin monthly PET_PM = %.3f mm; PET_ERA5 = %.3f mm; PEV/PM ratio = %.2f ",
  .pm_mean_mm, .pev_mean_mm, .pet_pm_era5_ratio
))

if (!is.finite(.pet_pm_era5_ratio)) {
  stop("PET_PM / PET_ERA5 basin means are missing or zero; PEV/PM ratio is non-finite. Inspect the upstream PET NetCDFs.")
}

pet_era5_qa_status <- "PASS"
if (.pet_pm_era5_ratio > PET_RATIO_REVIEW_RANGE[2] ||
    .pet_pm_era5_ratio < PET_RATIO_REVIEW_RANGE[1]) {
  pet_era5_qa_status <- "REVIEW"
  log_event(sprintf(
    "WARNING: PEV/PM = %.2f is outside the diagnostic %.1f-%.1f range. PET_ERA5 is retained UNSCALED as a labelled formulation-sensitivity product (not rescaled, not rejected).",
    .pet_pm_era5_ratio, PET_RATIO_REVIEW_RANGE[1], PET_RATIO_REVIEW_RANGE[2]
  ))
}

pet_era5_qa <- list(
  status              = pet_era5_qa_status,
  pev_pm_ratio        = .pet_pm_era5_ratio,
  pm_mean_mm_month    = .pm_mean_mm,
  pev_mean_mm_month   = .pev_mean_mm,
  review_range        = PET_RATIO_REVIEW_RANGE,
  fatal               = FALSE,
  provenance_check    = "2026-09-29 July-1950: CDS hourly 22.05 / CDS standard monthly 21.94 / archive 21.94 mm/day",
  interpretation      = "ERA5-Land PEV; analysed unscaled as PET-formulation sensitivity"
)
fwrite(
  data.table(
    status            = pet_era5_qa$status,
    pev_pm_ratio      = pet_era5_qa$pev_pm_ratio,
    pm_mean_mm_month  = pet_era5_qa$pm_mean_mm_month,
    pev_mean_mm_month = pet_era5_qa$pev_mean_mm_month,
    review_min        = PET_RATIO_REVIEW_RANGE[1],
    review_max        = PET_RATIO_REVIEW_RANGE[2],
    provenance_check  = pet_era5_qa$provenance_check
  ),
  file.path(out_dir, "pet_source_qa.csv")
)
log_series_quality(precip_monthly_avg,   "Basin Pr monthly", 0)
log_series_quality(pet_monthly_avg,      "Basin PET monthly", min_positive_value)
log_series_quality(pet_era5_monthly_avg, "Basin PET_ERA5 monthly", min_positive_value)
log_series_quality(aet_monthly_avg,      "Basin AET monthly", 0)
log_series_quality(tair_monthly_avg,     "Basin Tair monthly", NA_real_)
log_series_quality(precip_annual_avg,    "Basin Pr annual", 0)
log_series_quality(pet_annual_avg,       "Basin PET annual", min_positive_value)
log_series_quality(pet_era5_annual_avg,  "Basin PET_ERA5 annual", min_positive_value)
log_series_quality(aet_annual_avg,       "Basin AET annual", 0)
log_series_quality(tair_annual_avg,      "Basin Tair annual", NA_real_)
basin_precip_results    <- process_basin_series(precip_monthly_avg,  precip_annual_avg,
                                                "Precipitation", is_precip=TRUE)
basin_pet_results       <- process_basin_series(pet_monthly_avg,     pet_annual_avg,
                                                "PET", is_precip=FALSE)
basin_pet_era5_results  <- process_basin_series(pet_era5_monthly_avg, pet_era5_annual_avg,
                                                "PET_ERA5", is_precip=FALSE)
basin_aet_results       <- process_basin_series(aet_monthly_avg,     aet_annual_avg,       # v9
                                                "AET", is_precip=TRUE)
basin_tair_results      <- process_basin_series(tair_monthly_avg,    tair_annual_avg,
                                                "Temperature", is_precip=FALSE,
                                                skip_min_filter=TRUE)

pm_control_basin_results <- lapply(names(pm_control_monthly_avg), function(nm) {
  if (nm == "Temperature" && PM_T_IS_TAIR) return(.relabel(basin_tair_results, PM_CONTROL_LABELS[[nm]]))
  process_basin_series(
    pm_control_monthly_avg[[nm]],
    pm_control_annual_avg[[nm]],
    PM_CONTROL_LABELS[[nm]],
    is_precip = FALSE,
    skip_min_filter = TRUE
  )
})
names(pm_control_basin_results) <- names(pm_control_monthly_avg)

basin_swe_results <- basin_swe_peak_results <- NULL
if (RUN_SWE) {
  log_series_quality(swe_monthly_avg,  "Basin SWE monthly", 0)
  log_series_quality(swe_marapr_basin, "Basin SWE Mar-Apr annual", 0)
  basin_swe_results      <- process_basin_series(swe_monthly_avg, swe_marapr_basin, "SWE", is_precip=TRUE)
  basin_swe_peak_results <- process_basin_series(swe_monthly_avg, swe_peak_basin, "SWE_Peak",
                                                 is_precip=TRUE, annual_only=TRUE)
}

basin_extra_results <- lapply(names(extra_monthly_avg), function(nm) {
  cfg <- EXTRA_VARS[[nm]]
  process_basin_series(extra_monthly_avg[[nm]], extra_annual_avg[[nm]], nm,
                       is_precip = cfg$is_precip, skip_min_filter = cfg$skip_min_filter)
})
names(basin_extra_results) <- names(extra_monthly_avg)

####################################################################################
# -- PIXEL-LEVEL PROCESSING
####################################################################################
process_variable_final <- function(data_matrix, var_name, coords_dt,
                                   is_precip=FALSE, skip_min_filter=FALSE,
                                   yrs=years, mos=months,
                                   annual_override=NULL, run_monthly=TRUE) {
  log_event(paste("Processing", var_name, "pixels..."))
  # v9: pixel and basin runs share one month index; refuse a different calendar.
  if (!identical(as.integer(mos), as.integer(months)) || !identical(as.integer(yrs), as.integer(years)))
    stop(sprintf("[%s] calendar passed to process_variable_final() differs from the precipitation calendar.", var_name))
  if (!is.null(annual_override)) {
    log_event("  Using supplied annual metric matrix...")
    annual_matrix <- annual_override
  } else if (!is_precip && !skip_min_filter) {
    log_event("  Aggregating PET to warm-season annual sum (Apr-Oct)...")
    annual_matrix <- aggregate_to_annual_warmseason_fast(
      data_matrix, yrs, months_vec=mos, warm_months=warm_months_pet, method="sum")
  } else if (skip_min_filter) {
    log_event("  Aggregating to full-year annual mean...")
    annual_matrix <- aggregate_to_annual_fast(data_matrix, yrs, method="mean")
  } else {
    log_event("  Aggregating to annual sum...")
    annual_matrix <- aggregate_to_annual_fast(data_matrix, yrs, method="sum")
  }
  
  conf_env          <- new.env(parent=emptyenv())
  start_year_annual <- min(as.integer(rownames(annual_matrix)))
  
  log_event("  Running VC + TFPW + Spectral + Changepoint (annual, parallel)...")
  matrix_quality_profile(annual_matrix, sprintf("%s annual aggregated", var_name),
                         floor_value = if (is_precip) 0 else if (skip_min_filter) NA_real_ else min_positive_value)
  res_annual <- future_lapply(
    seq_len(ncol(annual_matrix)),
    function(i) run_trend_bundle(
      annual_matrix[,i], is_precip, alpha, max_tie_percent,
      n_sim_spectral, conf_env,
      start_year      = start_year_annual,
      skip_min_filter = skip_min_filter),
    future.seed=TRUE)
  
  annual_dt <- build_annual_dt(res_annual, coords_dt, var_name)
  annual_dt[, variable := var_name]
  
  if (!run_monthly) return(annual_dt)
  log_event("  Processing 12 calendar months (parallel)...")
  monthly_results_list <- vector("list", 12L)
  for (m in 1:12) {
    mi_key  <- as.character(m)
    mi_use  <- month_index_list[[mi_key]]     # same index as the basin runs
    monthly_sub <- if (length(mi_use)) data_matrix[mi_use, , drop=FALSE] else
      matrix(NA_real_, 0, ncol(data_matrix))
    mn_count  <- if (!skip_min_filter && !is_precip) min_nonzero_pet_monthly else NULL
    matrix_quality_profile(monthly_sub, sprintf("%s month %02d aggregated", var_name, m),
                           floor_value = if (is_precip) 0 else if (skip_min_filter) NA_real_ else min_positive_value)
    
    res_m <- future_lapply(
      seq_len(ncol(monthly_sub)),
      function(i) run_trend_bundle(
        monthly_sub[,i], is_precip, alpha, max_tie_percent,
        n_sim_spectral, conf_env,
        min_nonzero_count = mn_count,
        skip_min_filter   = skip_min_filter),
      future.seed=TRUE)
    
    vc_m  <- unpack_results(res_m, "vc")
    tf_m  <- unpack_results(res_m, "tf")
    spec_m <- data.table(
      n_spectral_peaks    = vapply(lapply(res_m,`[[`,"spec"),`[[`,integer(1),"n_peaks"),
      dominant_period     = vapply(lapply(res_m,`[[`,"spec"),`[[`,numeric(1),"dominant_period"),
      spectral_confidence = vapply(lapply(res_m,`[[`,"spec"),`[[`,numeric(1),"conf"))
    cpt_m  <- data.table(
      changepoint_detected   = vapply(lapply(res_m,`[[`,"cpt"),`[[`,logical(1),"changepoint_detected"),
      changepoint_position   = vapply(lapply(res_m,`[[`,"cpt"),`[[`,integer(1),"changepoint_position"),
      n_changepoints         = vapply(lapply(res_m,`[[`,"cpt"),`[[`,integer(1),"n_changepoints"),
      first_changepoint_year = vapply(lapply(res_m,`[[`,"cpt"),`[[`,integer(1),"first_changepoint_year"),
      mean_before_shift      = vapply(lapply(res_m,`[[`,"cpt"),`[[`,numeric(1),"mean_before"),
      mean_after_shift       = vapply(lapply(res_m,`[[`,"cpt"),`[[`,numeric(1),"mean_after"),
      magnitude_shift        = vapply(lapply(res_m,`[[`,"cpt"),`[[`,numeric(1),"magnitude_shift"))
    
    monthly_results_list[[m]] <- cbind(
      coords_dt, variable=var_name, period="monthly", month=m, variable_name=var_name,
      data.table(tau_vc=vc_m$tau, p_value_vc=vc_m$p, sl_vc=vc_m$sl,
                 vc_corrected=vc_m$extra, n_ties_vc=vc_m$n_ties,
                 percent_ties_vc=vc_m$pct_ties, n_min_vals_vc=vc_m$n_min,
                 percent_min_vals_vc=vc_m$pct_min, tau_b_adjusted_vc=vc_m$tau_b_adj,
                 filtered_vc=vc_m$filtered, filter_reason_vc=vc_m$reason),
      data.table(tau_tfpw=tf_m$tau, p_value_tfpw=tf_m$p, sl_tfpw=tf_m$sl,
                 tfpw_applied=tf_m$extra, n_ties_tfpw=tf_m$n_ties,
                 percent_ties_tfpw=tf_m$pct_ties, n_min_vals_tfpw=tf_m$n_min,
                 percent_min_vals_tfpw=tf_m$pct_min, tau_b_adjusted_tfpw=tf_m$tau_b_adj,
                 filtered_tfpw=tf_m$filtered, filter_reason_tfpw=tf_m$reason),
      n=NA_integer_, rho1_vc=NA_real_, rho1_tfpw=NA_real_,
      spec_m, cpt_m,
      same_significance=(vc_m$p < alpha)==(tf_m$p < alpha),
      same_direction=sign(vc_m$tau)==sign(tf_m$tau))
  }
  monthly_results <- rbindlist(monthly_results_list, use.names=TRUE, fill=TRUE)
  rbindlist(list(annual_dt, monthly_results), use.names=TRUE, fill=TRUE)
}

log_event("=== STARTING PRECIPITATION ANALYSIS ===")
precip_results   <- process_variable_final(precip_matrix,  "Precipitation", coords_dt,
                                           is_precip=TRUE)
log_event("=== STARTING PET (Penman-Monteith) ANALYSIS ===")
pet_results      <- process_variable_final(pet_matrix,     "PET",       coords_dt,
                                           is_precip=FALSE)
log_event("=== STARTING PET (ERA5-Land PEV) TREND ANALYSIS ===")
pet_era5_results <- process_variable_final(pet_era5_matrix, "PET_ERA5", coords_dt,
                                           is_precip=FALSE)
log_event("=== STARTING AET (ERA5-Land total evaporation) TREND ANALYSIS (annual = calendar-year total) ===")
aet_results      <- process_variable_final(aet_matrix,     "AET",       coords_dt,
                                           is_precip=TRUE, annual_override=aet_annual_pix)
log_event("=== STARTING TEMPERATURE ANALYSIS ===")
tair_results     <- process_variable_final(tair_matrix,    "Temperature", coords_dt,
                                           is_precip=FALSE, skip_min_filter=TRUE)

pm_control_results <- lapply(names(pm_control_matrices), function(nm) {
  if (nm == "Temperature" && PM_T_IS_TAIR) return(.relabel(tair_results, PM_CONTROL_LABELS[[nm]]))
  process_variable_final(
    pm_control_matrices[[nm]],
    PM_CONTROL_LABELS[[nm]],
    coords_dt,
    is_precip = FALSE,
    skip_min_filter = TRUE
  )
})
names(pm_control_results) <- names(pm_control_matrices)

swe_results <- swe_peak_results <- NULL
if (RUN_SWE) {
  log_event("=== STARTING SWE (ERA5-Land) TREND ANALYSIS (annual = Mar-Apr mean; + 12 months) ===")
  swe_results <- process_variable_final(swe_matrix, "SWE", coords_dt, is_precip=TRUE,
                                        annual_override=swe_marapr_pix)
  log_event("=== STARTING SWE_Peak (annual maximum monthly SWE, Jan-Jun) TREND ANALYSIS ===")
  swe_peak_results <- process_variable_final(swe_matrix, "SWE_Peak", coords_dt, is_precip=TRUE,
                                             annual_override=swe_peak_pix, run_monthly=FALSE)
}

extra_results <- lapply(names(extra_matrix), function(nm) {
  cfg <- EXTRA_VARS[[nm]]
  log_event(sprintf("=== STARTING %s TREND ANALYSIS (annual = %s) ===", nm,
                    if (cfg$annual == "sum") "calendar-year total" else "Mar-Apr mean"))
  process_variable_final(extra_matrix[[nm]], nm, coords_dt,
                         is_precip = cfg$is_precip, skip_min_filter = cfg$skip_min_filter,
                         annual_override = extra_annual_pix[[nm]])
})
names(extra_results) <- names(extra_matrix)

####################################################################################
# -- SPECIFIC POINT EXTRACTION AND TREND ANALYSIS
####################################################################################
log_event("=== SPECIFIC POINT ANALYSIS (all configured variables) ===")

pts_wgs84  <- vect(SPECIFIC_PTS, geom=c("x","y"), crs="EPSG:4326")
pts_proj   <- project(pts_wgs84, terra::crs(analysis_template))
n_pts      <- nrow(pts_proj)
cat(sprintf("  Extracting at %d specific point(s)...\n", n_pts))

extract_point_monthly <- function(rast_clipped, pts, dates_vec, var_name, unit_label="") {
  rows <- vector("list", nrow(pts))
  for (i in seq_len(nrow(pts))) {
    pt_v  <- pts[i, ]
    vals  <- tryCatch(as.numeric(terra::extract(rast_clipped, pt_v)[1, -1]),
                      error=function(e) rep(NA_real_, nlyr(rast_clipped)))
    n_use <- min(length(dates_vec), length(vals))
    rows[[i]] <- data.frame(
      point_id  = i,
      lon_wgs84 = unname(crds(pts_wgs84)[i, 1]),
      lat_wgs84 = unname(crds(pts_wgs84)[i, 2]),
      date      = dates_vec[seq_len(n_use)],
      year      = as.integer(format(dates_vec[seq_len(n_use)], "%Y")),
      month     = as.integer(format(dates_vec[seq_len(n_use)], "%m")),
      value     = vals[seq_len(n_use)],
      variable  = var_name,
      unit      = unit_label,
      stringsAsFactors = FALSE
    )
  }
  rbindlist(rows)
}

pt_pr       <- extract_point_monthly(precip_clipped,   pts_proj, dates, "Precipitation", "mm/month")
pt_pet      <- extract_point_monthly(pet_clipped,      pts_proj, dates, "PET",           "mm/month")
pt_pet_era5 <- extract_point_monthly(pet_era5_clipped, pts_proj, dates, "PET_ERA5",      "mm/month")
pt_aet      <- extract_point_monthly(aet_clipped,      pts_proj, dates, "AET",           "mm/month")
pt_tair     <- extract_point_monthly(tair_clipped,     pts_proj, dates, "Temperature",   "degC")
pt_swe      <- if (RUN_SWE) extract_point_monthly(swe_clipped, pts_proj, dates, "SWE", "mm w.e.") else NULL

pt_pm_controls <- lapply(names(pm_control_clipped), function(nm) {
  extract_point_monthly(
    pm_control_clipped[[nm]],
    pts_proj,
    dates,
    PM_CONTROL_LABELS[[nm]],
    PM_CONTROL_UNITS[[nm]]
  )
})
names(pt_pm_controls) <- names(pm_control_clipped)

pt_extra <- lapply(names(extra_clipped), function(nm)
  extract_point_monthly(extra_clipped[[nm]], pts_proj, dates, nm, EXTRA_VARS[[nm]]$unit))
names(pt_extra) <- names(extra_clipped)

point_monthly_ts <- rbindlist(
  c(list(pt_pr, pt_pet, pt_pet_era5, pt_aet, pt_tair, pt_swe), pt_pm_controls, pt_extra),
  use.names=TRUE, fill=TRUE
)
fwrite(point_monthly_ts, file.path(out_dir, "point_monthly_timeseries.csv"))
log_event(sprintf("[OK] Saved point_monthly_timeseries.csv  (%d rows)", nrow(point_monthly_ts)))

run_point_trends <- function(pt_df, var_name, is_precip,
                             skip_min_filter = FALSE,
                             warm_season_annual = FALSE,
                             annual_metric = NULL,
                             annual_only = FALSE) {
  results <- list()
  for (pid in unique(pt_df$point_id)) {
    sub <- pt_df[point_id == pid]
    lon <- sub$lon_wgs84[1]; lat <- sub$lat_wgs84[1]
    conf_env <- new.env(parent=emptyenv())
    
    if (!is.null(annual_metric)) {
      ann_mat <- annual_metric(sub$value, sub$year, sub$month)
    } else if (warm_season_annual) {
      ann_mat <- aggregate_to_annual_warmseason_fast(
        matrix(sub$value, ncol=1), sub$year, months_vec=sub$month,
        warm_months=warm_months_pet, method="sum")
    } else if (skip_min_filter) {
      ann_mat <- aggregate_to_annual_fast(matrix(sub$value, ncol=1), sub$year, method="mean")
    } else {
      ann_mat <- aggregate_to_annual_fast(matrix(sub$value, ncol=1), sub$year, method="sum")
    }
    ann_vec  <- ann_mat[, 1]
    ann_yrs  <- as.integer(rownames(ann_mat))
    
    res_ann <- run_trend_bundle(
      ann_vec, is_precip, alpha, max_tie_percent, n_sim_spectral, conf_env,
      start_year=min(ann_yrs, na.rm=TRUE), skip_min_filter=skip_min_filter)
    vc_a <- res_ann$vc; tf_a <- res_ann$tf; cpt_a <- res_ann$cpt
    
    results[[length(results)+1]] <- data.frame(
      point_id=pid, lon_wgs84=lon, lat_wgs84=lat,
      variable=var_name, period="annual", month=NA_integer_, month_abb="Annual",
      n_obs=vc_a$n,
      year_min=min(ann_yrs, na.rm=TRUE), year_max=max(ann_yrs, na.rm=TRUE),
      tau_vc=vc_a$tau, p_value_vc=vc_a$p, sl_vc=vc_a$sl,
      filtered_vc=vc_a$filtered, filter_reason_vc=vc_a$reason,
      tau_tfpw=tf_a$tau, p_value_tfpw=tf_a$p, sl_tfpw=tf_a$sl,
      filtered_tfpw=tf_a$filtered, filter_reason_tfpw=tf_a$reason,
      rho1=vc_a$rho1,
      changepoint_detected=cpt_a$changepoint_detected,
      first_changepoint_year=cpt_a$first_changepoint_year,
      magnitude_shift=cpt_a$magnitude_shift,
      sig_vc  =!is.na(vc_a$p) && vc_a$p < alpha,
      sig_tfpw=!is.na(tf_a$p) && tf_a$p < alpha,
      stringsAsFactors=FALSE)
    
    if (annual_only) next
    for (m in 1:12) {
      sub_m <- sub[month == m]
      if (nrow(sub_m) < 10) next
      mn_count <- if (!skip_min_filter && !is_precip) min_nonzero_pet_monthly else NULL
      res_m <- run_trend_bundle(
        sub_m$value, is_precip, alpha, max_tie_percent, n_sim_spectral, conf_env,
        min_nonzero_count=mn_count, skip_min_filter=skip_min_filter)
      vc_m <- res_m$vc; tf_m <- res_m$tf; cpt_m <- res_m$cpt
      
      results[[length(results)+1]] <- data.frame(
        point_id=pid, lon_wgs84=lon, lat_wgs84=lat,
        variable=var_name, period="monthly", month=m, month_abb=month.abb[m],
        n_obs=vc_m$n,
        year_min=min(sub_m$year, na.rm=TRUE), year_max=max(sub_m$year, na.rm=TRUE),
        tau_vc=vc_m$tau, p_value_vc=vc_m$p, sl_vc=vc_m$sl,
        filtered_vc=vc_m$filtered, filter_reason_vc=vc_m$reason,
        tau_tfpw=tf_m$tau, p_value_tfpw=tf_m$p, sl_tfpw=tf_m$sl,
        filtered_tfpw=tf_m$filtered, filter_reason_tfpw=tf_m$reason,
        rho1=vc_m$rho1,
        changepoint_detected=cpt_m$changepoint_detected,
        first_changepoint_year=cpt_m$first_changepoint_year,
        magnitude_shift=cpt_m$magnitude_shift,
        sig_vc  =!is.na(vc_m$p) && vc_m$p < alpha,
        sig_tfpw=!is.na(tf_m$p) && tf_m$p < alpha,
        stringsAsFactors=FALSE)
    }
  }
  rbindlist(results, fill=TRUE)
}

pt_pr_trends       <- run_point_trends(pt_pr,       "Precipitation", is_precip=TRUE)
pt_pet_trends      <- run_point_trends(pt_pet,      "PET",           is_precip=FALSE,
                                       warm_season_annual=TRUE)
pt_pet_era5_trends <- run_point_trends(pt_pet_era5, "PET_ERA5",      is_precip=FALSE,
                                       warm_season_annual=TRUE)
pt_aet_trends      <- run_point_trends(pt_aet,      "AET",           is_precip=TRUE)          # v9
pt_tair_trends     <- run_point_trends(pt_tair,     "Temperature",   is_precip=FALSE,
                                       skip_min_filter=TRUE)
pt_swe_trends      <- if (RUN_SWE) run_point_trends(pt_swe, "SWE", is_precip=TRUE,
                                                    annual_metric=function(v, yr, mo) annual_metric_fast(v, yr, mo, SWE_ANNUAL_MONTHS, "mean")) else NULL
pt_swe_peak_trends <- if (RUN_SWE) run_point_trends(pt_swe, "SWE_Peak", is_precip=TRUE,                # v9
                                                    annual_metric=function(v, yr, mo) annual_metric_fast(v, yr, mo, SWE_PEAK_MONTHS, "max"),
                                                    annual_only=TRUE) else NULL

pm_control_point_trends <- lapply(names(pt_pm_controls), function(nm) {
  if (nm == "Temperature" && PM_T_IS_TAIR) return(.relabel(pt_tair_trends, PM_CONTROL_LABELS[[nm]]))
  run_point_trends(
    pt_pm_controls[[nm]],
    PM_CONTROL_LABELS[[nm]],
    is_precip = FALSE,
    skip_min_filter = TRUE
  )
})
names(pm_control_point_trends) <- names(pt_pm_controls)

pt_extra_trends <- lapply(names(pt_extra), function(nm) {
  cfg <- EXTRA_VARS[[nm]]
  run_point_trends(pt_extra[[nm]], nm, is_precip = cfg$is_precip, skip_min_filter = cfg$skip_min_filter,
                   annual_metric = if (cfg$annual == "marapr")
                     function(v, yr, mo) annual_metric_fast(v, yr, mo, SWE_ANNUAL_MONTHS, "mean") else NULL)
})

point_trend_stats <- rbindlist(
  c(list(pt_pr_trends, pt_pet_trends, pt_pet_era5_trends, pt_aet_trends, pt_tair_trends,
         pt_swe_trends, pt_swe_peak_trends),
    pm_control_point_trends, pt_extra_trends),
  use.names=TRUE, fill=TRUE)
point_trend_stats[, tfpw_method := "Yue2002_shared_J4"]
fwrite(point_trend_stats, file.path(out_dir, "point_trend_stats.csv"))
log_event(sprintf("[OK] Saved point_trend_stats.csv  (%d rows x %d cols)",
                  nrow(point_trend_stats), ncol(point_trend_stats)))

cat("\n-- Specific Point Trend Summary --\n")
smry_pts <- point_trend_stats[, .(
  n_periods  = .N,
  sig_vc_n   = sum(sig_vc,   na.rm=TRUE),
  sig_tfpw_n = sum(sig_tfpw, na.rm=TRUE)
), by=.(point_id, variable)]
print(as.data.frame(smry_pts))

####################################################################################
# -- COMBINE ALL RESULTS & SAVE
####################################################################################
log_event("Combining all pixel-level and basin results (all configured variables)...")
all_results <- rbindlist(
  c(
    list(
      precip_results, pet_results, pet_era5_results, aet_results, tair_results,
      basin_precip_results, basin_pet_results, basin_pet_era5_results, basin_aet_results,
      basin_tair_results
    ),
    pm_control_results,
    pm_control_basin_results,
    Filter(Negate(is.null), list(swe_results, swe_peak_results, basin_swe_results, basin_swe_peak_results)),
    extra_results,
    basin_extra_results
  ),
  use.names=TRUE, fill=TRUE)
all_results[, tfpw_method := "Yue2002_shared_J4"]

saveRDS(all_results, file.path(out_dir, "all_results.rds"), compress="gzip")

saveRDS(list(
  precip_results         = precip_results,
  pet_results            = pet_results,
  pet_era5_results       = pet_era5_results,
  aet_results            = aet_results,
  tair_results           = tair_results,
  basin_precip_results   = basin_precip_results,
  basin_pet_results      = basin_pet_results,
  basin_pet_era5_results = basin_pet_era5_results,
  basin_aet_results      = basin_aet_results,
  basin_tair_results     = basin_tair_results,
  swe_results            = swe_results,
  swe_peak_results       = swe_peak_results,
  basin_swe_results      = basin_swe_results,
  basin_swe_peak_results = basin_swe_peak_results,
  swe_config = list(file = swe_path, units_in = SWE_UNITS, annual_months = SWE_ANNUAL_MONTHS,
                    peak_months = SWE_PEAK_MONTHS, step_year = SWE_STEP_YEAR, run = RUN_SWE),
  swe_input_qa           = swe_qa,
  negative_clamp_qa      = rbindlist(NEG_CLAMP_QA),
  aet_config             = list(file = aet_path, sign_applied = aet_sign),
  basin_pixels          = n_basin_pixels,
  bbox_cells            = n_bbox_cells,
  reduction_pct         = reduction_pct,
  processing_date       = Sys.time(),
  pet_era5_qa           = pet_era5_qa,
  parameters = list(alpha=alpha, max_tie_percent=max_tie_percent,
                    tfpw_method = "Yue2002_shared_J4",
                    tfpw_reference_implementation = "modifiedmk::tfpwmk (QA only)",
                    max_min_value_pct_precip=max_min_value_pct_precip,
                    max_min_value_pct_pet=max_min_value_pct_pet,
                    pet_zero_floor = "applied before basin averaging (v9)",
                    pm_control_files=PM_CONTROL_NC,
                    pm_control_units=PM_CONTROL_UNITS),
  clipped_extent = list(xmin=xmin(clipped_template), xmax=xmax(clipped_template),
                        ymin=ymin(clipped_template), ymax=ymax(clipped_template),
                        nrows=nrow(clipped_template), ncols=ncol(clipped_template),
                        res=res(clipped_template)),
  basin_avg_monthly = data.frame(
    date               = dates,
    precip_mm_month    = precip_monthly_avg,
    pet_mm_month       = pet_monthly_avg,
    pet_era5_mm_month  = pet_era5_monthly_avg,
    aet_mm_month       = aet_monthly_avg,
    swe_mm_month       = swe_monthly_avg,
    tair_degC_month    = tair_monthly_avg,
    PM_Temperature     = pm_control_monthly_avg$Temperature,
    PM_NetRadiation    = pm_control_monthly_avg$NetRadiation,
    PM_VPD             = pm_control_monthly_avg$VPD,
    PM_Wind2m          = pm_control_monthly_avg$Wind2m,
    PM_SurfacePressure = pm_control_monthly_avg$SurfacePressure
  ),
  basin_avg_annual = data.frame(
    year               = annual_years,
    precip_mm_year     = precip_annual_avg,
    pet_mm_year        = pet_annual_avg,
    pet_era5_mm_year   = pet_era5_annual_avg,
    aet_mm_year        = aet_annual_avg,
    swe_marapr_mm_year = swe_marapr_basin,
    swe_peak_mm_year   = swe_peak_basin,
    tair_degC_year     = tair_annual_avg,
    PM_Temperature     = pm_control_annual_avg$Temperature,
    PM_NetRadiation    = pm_control_annual_avg$NetRadiation,
    PM_VPD             = pm_control_annual_avg$VPD,
    PM_Wind2m          = pm_control_annual_avg$Wind2m,
    PM_SurfacePressure = pm_control_annual_avg$SurfacePressure
  ),
  specific_pts = SPECIFIC_PTS
), file.path(out_dir, "analysis_metadata.rds"))

summary_stats <- all_results[, .(
  n_total            = .N,
  n_valid_vc         = sum(!filtered_vc & !is.na(p_value_vc)),
  n_significant_vc   = sum(!filtered_vc & p_value_vc < alpha, na.rm=TRUE),
  n_valid_tfpw       = sum(!filtered_tfpw & !is.na(p_value_tfpw)),
  n_significant_tfpw = sum(!filtered_tfpw & p_value_tfpw < alpha, na.rm=TRUE)
), by=.(variable, period, month, is_basin_average)]
fwrite(summary_stats, file.path(out_dir, "summary_statistics.csv"))

log_event("------------------------------------------------------------")
log_event("DIAGNOSTIC SUMMARY: root-cause screen for MK / TFPW warnings.")
diag_summary <- log_result_root_cause_summary(all_results)
if (!is.null(diag_summary)) {
  print(diag_summary)
}
log_event("------------------------------------------------------------")
log_event("TREND PROCESSING COMPLETE (Pr + PET_PM + PET_ERA5 + AET + Temperature + PM controls + SWE + runoff/snow)")
cat("\n[OK] TREND PROCESSING COMPLETED\n")
cat(sprintf("  Basin pixels: %d  |  Specific points: %d\n", n_basin_pixels, n_pts))
cat(sprintf("  Variables: Precipitation, PET (PM), PET (ERA5-Land PEV), AET, Temperature, PM T2m, Rn, VPD, 2-m wind, surface pressure%s%s\n",
            if (RUN_SWE) ", SWE, SWE_Peak" else "",
            if (length(EXTRA_VARS)) paste0(", ", paste(names(EXTRA_VARS), collapse = ", ")) else ""))
cat(sprintf("  Outputs: %s\n", normalizePath(out_dir)))

future::plan(future::sequential)
####################################################################################
# REGIONAL EXTENSION - 7 SUBBASINS + WHOLE NECHAKO BASIN
####################################################################################
log_event("Building regional monthly/annual series for 7 subbasins + whole basin...")

regional_regions <- load_nechako_regions(template = analysis_template)
validate_nechako_regions(regional_regions)
regional_weights_all <- region_fractional_weights(analysis_template, regional_regions)

regional_weights <- lapply(regional_weights_all, function(x) {
  x$weights <- x$weights[valid_mask]
  x$coverage <- x$coverage[valid_mask]
  x
})

regional_mean_from_matrix <- function(mat, dates_vec, variable_name) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (nrow(mat) != length(dates_vec))
    stop(sprintf(
      "Regional climate extension: %s has %d time rows but %d dates.",
      variable_name, nrow(mat), length(dates_vec)))
  
  out <- lapply(names(regional_weights), function(id) {
    rw <- regional_weights[[id]]
    w <- rw$weights
    if (length(w) != ncol(mat))
      stop(sprintf(
        "Regional climate extension: %s/%s has %d pixel columns but %d weights.",
        variable_name, id, ncol(mat), length(w)))
    
    data.table(
      region_id = id,
      region_name = rw$region_name,
      date = as.Date(dates_vec),
      value = regional_weighted_mean(t(mat), w, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE, total_weight = rw$area_m2),
      variable = variable_name,
      area_m2 = rw$area_m2,
      n_cells = rw$n_cells,
      n_partial_cells = rw$n_partial
    )
  })
  rbindlist(out, use.names = TRUE, fill = TRUE)
}

if (ncol(tair_matrix) != ncol(precip_matrix))
  stop("Regional climate extension: temperature and precipitation pixel counts differ.")
for (.m in list(precip_matrix, pet_matrix, pet_era5_matrix, aet_matrix, tair_matrix))
  if (nrow(.m) != length(dates)) stop("Regional climate extension: matrix/date count mismatch.")

.keep_cols <- c("region_id", "region_name", "date", "area_m2", "n_cells", "n_partial_cells")
.part <- function(rr, col) { out <- rr[, c(.keep_cols, "value"), with = FALSE]; data.table::setnames(out, "value", col); out }

regional_pr  <- regional_mean_from_matrix(precip_matrix,   dates, "precip_mm_month")
regional_pm  <- regional_mean_from_matrix(pet_matrix,      dates, "pet_mm_month")
regional_e5  <- regional_mean_from_matrix(pet_era5_matrix, dates, "pet_era5_mm_month")
regional_aet <- regional_mean_from_matrix(aet_matrix,      dates, "aet_mm_month")
regional_t   <- regional_mean_from_matrix(tair_matrix,     dates, "tair_degC_month")
regional_swe_part <- if (RUN_SWE) list(.part(regional_mean_from_matrix(swe_matrix, dates, "swe_mm_month"), "swe_mm_month")) else list()

regional_extra_parts <- lapply(names(extra_matrix), function(nm)
  .part(regional_mean_from_matrix(extra_matrix[[nm]], dates, nm), EXTRA_VARS[[nm]]$col_month))

# v9: PM controls in the regional tables (column names = PM_CONTROL_LABELS)
regional_pm_control_parts <- lapply(names(pm_control_matrices), function(nm) {
  lab <- PM_CONTROL_LABELS[[nm]]
  .part(regional_mean_from_matrix(pm_control_matrices[[nm]], dates, lab), lab)
})

regional_monthly <- Reduce(function(x, y) merge(x, y, by = .keep_cols, all = TRUE, sort = TRUE),
                           c(list(.part(regional_pr,  "precip_mm_month"),
                                  .part(regional_pm,  "pet_mm_month"),
                                  .part(regional_e5,  "pet_era5_mm_month"),
                                  .part(regional_aet, "aet_mm_month"),
                                  .part(regional_t,   "tair_degC_month")),
                             regional_swe_part, regional_extra_parts, regional_pm_control_parts))
if (!"swe_mm_month" %in% names(regional_monthly)) regional_monthly[, swe_mm_month := NA_real_]
setorder(regional_monthly, region_id, date)

regional_monthly[, year := as.integer(format(date, "%Y"))]
regional_monthly[, month := as.integer(format(date, "%m"))]
regional_monthly[, season := fifelse(month %in% c(12L,1L,2L), "DJF",
                                     fifelse(month %in% 3:5, "MAM",
                                             fifelse(month %in% 6:8, "JJA", "SON")))]
regional_monthly[, pet_pm_fullyear_component := pet_mm_month]
regional_monthly[, pet_era5_fullyear_component := pet_era5_mm_month]

regional_annual <- regional_monthly[, {
  pet_warm <- pet_mm_month[month %in% warm_months_pet]
  e5_warm  <- pet_era5_mm_month[month %in% warm_months_pet]
  list(
    precip_mm_year = if (sum(is.finite(precip_mm_month)) >= 6L) sum(precip_mm_month, na.rm=TRUE) else NA_real_,
    pet_mm_year = if (sum(is.finite(pet_warm)) >= ceiling(length(warm_months_pet)/2)) sum(pet_warm, na.rm=TRUE) else NA_real_,
    pet_era5_mm_year = if (sum(is.finite(e5_warm)) >= ceiling(length(warm_months_pet)/2)) sum(e5_warm, na.rm=TRUE) else NA_real_,
    pet_mm_year_full = if (sum(is.finite(pet_mm_month)) >= 6L) sum(pet_mm_month, na.rm=TRUE) else NA_real_,
    pet_era5_mm_year_full = if (sum(is.finite(pet_era5_mm_month)) >= 6L) sum(pet_era5_mm_month, na.rm=TRUE) else NA_real_,
    aet_mm_year = if (sum(is.finite(aet_mm_month)) >= 6L) sum(aet_mm_month, na.rm=TRUE) else NA_real_,
    tair_degC_year = if (sum(is.finite(tair_degC_month)) >= 6L) mean(tair_degC_month, na.rm=TRUE) else NA_real_,
    swe_marapr_mm = { v <- swe_mm_month[month %in% SWE_ANNUAL_MONTHS]; if (sum(is.finite(v)) == length(SWE_ANNUAL_MONTHS)) mean(v) else NA_real_ },
    swe_peak_mm   = { v <- swe_mm_month[month %in% SWE_PEAK_MONTHS];   if (sum(is.finite(v)) == length(SWE_PEAK_MONTHS))   max(v)  else NA_real_ },
    n_months = .N
  )
}, by=.(region_id, region_name, year, area_m2, n_cells, n_partial_cells)]

for (nm in names(extra_matrix)) {
  cm <- EXTRA_VARS[[nm]]$col_month; cy <- EXTRA_VARS[[nm]]$col_year
  ann <- if (EXTRA_VARS[[nm]]$annual == "sum") {
    regional_monthly[, .(v = { x <- get(cm); if (sum(is.finite(x)) >= 6L) sum(x, na.rm = TRUE) else NA_real_ }),
                     by = .(region_id, year)]
  } else {
    regional_monthly[month %in% SWE_ANNUAL_MONTHS,
                     .(v = { x <- get(cm); if (sum(is.finite(x)) == length(SWE_ANNUAL_MONTHS)) mean(x) else NA_real_ }),
                     by = .(region_id, year)]
  }
  data.table::setnames(ann, "v", cy)
  regional_annual <- merge(regional_annual, ann, by = c("region_id", "year"), all.x = TRUE)
}

# v9: PM controls - calendar-year mean and Apr-Oct mean per region
for (nm in names(pm_control_matrices)) {
  lab <- PM_CONTROL_LABELS[[nm]]
  ann <- regional_monthly[, .(
    v_year = { x <- get(lab); if (sum(is.finite(x)) >= 6L) mean(x, na.rm = TRUE) else NA_real_ },
    v_ao   = { x <- get(lab)[month %in% warm_months_pet]
    if (sum(is.finite(x)) >= ceiling(length(warm_months_pet) / 2)) mean(x, na.rm = TRUE) else NA_real_ }),
    by = .(region_id, year)]
  data.table::setnames(ann, c("v_year", "v_ao"), c(paste0(lab, "_year"), paste0(lab, "_AprOct")))
  regional_annual <- merge(regional_annual, ann, by = c("region_id", "year"), all.x = TRUE)
}

REGIONAL_DIR <- file.path(out_dir, "regional")
dir.create(REGIONAL_DIR, recursive=TRUE, showWarnings=FALSE)
fwrite(regional_monthly, file.path(REGIONAL_DIR, "regional_avg_monthly.csv"))
fwrite(regional_annual, file.path(REGIONAL_DIR, "regional_avg_annual.csv"))
fwrite(region_catalog(regional_regions), file.path(REGIONAL_DIR, "region_catalog.csv"))

metadata_path <- file.path(out_dir, "analysis_metadata.rds")
metadata <- readRDS(metadata_path)
metadata$regional_avg_monthly <- as.data.frame(regional_monthly)
metadata$regional_avg_annual  <- as.data.frame(regional_annual)
metadata$regional_regions <- region_catalog(regional_regions)
metadata$regional_source <- list(
  shapefile_directory = NECHAKO_SUBBASIN_DIR,
  native_crs = NECHAKO_SUBBASIN_SOURCE_CRS,
  analysis_grid_crs = terra::crs(analysis_template),
  whole_basin_definition = "dissolved union of 7 supplied subbasins",
  weighting = "fractional polygon coverage x equal-area cell area"
)
for (nm in names(extra_monthly_avg)) {
  metadata$basin_avg_monthly[[EXTRA_VARS[[nm]]$col_month]] <- extra_monthly_avg[[nm]]
  metadata$basin_avg_annual[[EXTRA_VARS[[nm]]$col_year]]   <- extra_annual_avg[[nm]]
}
metadata$extra_vars_config     <- lapply(EXTRA_VARS, `[`, c("path", "unit", "annual", "col_month", "col_year"))
metadata$snow_present_depth_cm <- SNOW_PRESENT_DEPTH_CM
saveRDS(metadata, metadata_path, compress=TRUE)

regional_pm_summary <- regional_monthly[, .(
  date = format(date, "%Y-%m"), year, month,
  region_id, region_name,
  mean_pet = pet_mm_month / lubridate::days_in_month(date)
)]
regional_e5_summary <- regional_monthly[, .(
  date = format(date, "%Y-%m"), year, month,
  region_id, region_name,
  mean_pet = pet_era5_mm_month / lubridate::days_in_month(date)
)]
fwrite(regional_pm_summary, file.path(REGIONAL_DIR, "regional_PM_PET_monthly_summary.csv"))
fwrite(regional_e5_summary, file.path(REGIONAL_DIR, "regional_ERA5Land_PEV_monthly_summary.csv"))

log_event(sprintf("REGIONAL PRODUCTS: %d regions | %d monthly rows | %d annual rows",
                  uniqueN(regional_monthly$region_id), nrow(regional_monthly), nrow(regional_annual)))
log_event("Regional extension complete.")

####################################################################################
# SWE HOMOGENEITY / STEP-CHANGE INVESTIGATION  (basin + 7 subbasins + every pixel)
####################################################################################
if (RUN_SWE) {
  log_event("=== SWE HOMOGENEITY / STEP-CHANGE DIAGNOSTICS ===")
  dir.create(swe_out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(20260101)
  
  homog_long <- list()
  for (id in unique(regional_monthly$region_id)) {
    rm_ <- regional_monthly[region_id == id][order(date)]
    yrs_all <- sort(unique(rm_$year))
    swe_ma <- annual_metric_fast(rm_$swe_mm_month, rm_$year, rm_$month, SWE_ANNUAL_MONTHS, "mean")[, 1]
    swe_pk <- annual_metric_fast(rm_$swe_mm_month, rm_$year, rm_$month, SWE_PEAK_MONTHS, "max")[, 1]
    p_acc  <- season_by_year(rm_$precip_mm_month,  rm_$year, rm_$month, SWE_ACCUM_MONTHS, fun = "sum")
    t_acc  <- season_by_year(rm_$tair_degC_month,  rm_$year, rm_$month, SWE_ACCUM_MONTHS, fun = "mean")
    
    resid <- rep(NA_real_, length(yrs_all)); r2 <- NA_real_
    d <- data.frame(y = as.numeric(swe_ma), p = as.numeric(p_acc), t = as.numeric(t_acc))
    okf <- is.finite(d$y) & is.finite(d$p) & is.finite(d$t)
    if (sum(okf) >= 30L && stats::sd(d$y[okf]) > 0) {
      fit <- stats::lm(y ~ p + t, data = d[okf, ])
      resid[okf] <- stats::residuals(fit); r2 <- summary(fit)$r.squared
    }
    homog_long[[id]] <- rbindlist(list(
      data.table(region_id = id, year = yrs_all, series = "SWE_MarApr",          value = as.numeric(swe_ma)),
      data.table(region_id = id, year = yrs_all, series = "SWE_Peak",            value = as.numeric(swe_pk)),
      data.table(region_id = id, year = yrs_all, series = "Precip_NovMar",       value = as.numeric(p_acc)),
      data.table(region_id = id, year = yrs_all, series = "Temp_NovMar",         value = as.numeric(t_acc)),
      data.table(region_id = id, year = yrs_all, series = "SWE_MarApr_resid_PT", value = resid)
    ))
    homog_long[[id]][series == "SWE_MarApr_resid_PT", fit_r2 := r2]
  }
  homog_long <- rbindlist(homog_long, fill = TRUE)
  fwrite(homog_long, file.path(swe_out_dir, "swe_homogeneity_annual_series.csv"))
  
  log_event("  [SWE-HOMOG] Running Pettitt / SNHT / Buishand / fixed-date / step+trend tests...")
  step_tbl <- rbindlist(lapply(split(homog_long, by = c("region_id", "series"), drop = TRUE), function(dd) {
    cbind(data.table(region_id = dd$region_id[1], series = dd$series[1]),
          homogeneity_tests(dd$value, dd$year),
          step_trend_model(dd$value, dd$year))
  }), fill = TRUE)
  step_tbl[, pettitt_p_fdr := stats::p.adjust(pettitt_p, "BH"), by = series]
  step_tbl[, step_p_fdr    := stats::p.adjust(step_p,    "BH"), by = series]
  fwrite(step_tbl, file.path(swe_out_dir, "swe_step_tests_by_region.csv"))
  
  log_event("  [SWE-HOMOG] Sub-period VC-MK / TFPW-MK trends (full, pre-step, post-step, 1991-2020)...")
  periods <- list(full = c(1950L, 2025L), pre_step = c(1950L, SWE_STEP_YEAR - 1L),
                  post_step = c(SWE_STEP_YEAR, 2025L), clim_1991_2020 = c(1991L, 2020L))
  sub_tbl <- rbindlist(lapply(split(homog_long[series %in% c("SWE_MarApr", "SWE_Peak", "Precip_NovMar", "Temp_NovMar")],
                                    by = c("region_id", "series"), drop = TRUE), function(dd) {
                                      is_t <- dd$series[1] == "Temp_NovMar"
                                      subperiod_trends(dd$value, dd$year, dd$series[1], dd$region_id[1], periods,
                                                       is_precip = !is_t, skip_min_filter = is_t)
                                    }), fill = TRUE)
  fwrite(sub_tbl, file.path(swe_out_dir, "swe_subperiod_trends.csv"))
  
  era_def <- list(`1950-1978` = c(1950L, 1978L), `1979-1990 (post-hoc window)` = c(1979L, 1990L),
                  `1979-2025` = c(1979L, 2025L), `1991-2020` = c(1991L, 2020L), `2021-2025` = c(2021L, 2025L))
  era_tbl <- rbindlist(lapply(names(era_def), function(en) {
    homog_long[year >= era_def[[en]][1] & year <= era_def[[en]][2],
               .(era = en, n_years = sum(is.finite(value)), mean_value = mean(value, na.rm = TRUE),
                 sd_value = stats::sd(value, na.rm = TRUE)), by = .(region_id, series)]
  }))
  fwrite(era_tbl, file.path(swe_out_dir, "swe_period_means.csv"))
  
  ref_tbl <- rbindlist(lapply(split(homog_long[series == "SWE_MarApr"], by = "region_id"), function(dd)
    reference_sensitivity(dd$value, dd$year, dd$region_id[1], "SWE_MarApr")))
  fwrite(ref_tbl, file.path(swe_out_dir, "swe_reference_sensitivity.csv"))
  
  log_event("  [SWE-HOMOG] Pixel-level Pettitt + fixed-date scan of Mar-Apr SWE...")
  yrs_pix <- as.integer(rownames(swe_marapr_pix))
  pix_rows <- lapply(seq_len(ncol(swe_marapr_pix)), function(i) {
    x <- swe_marapr_pix[, i]; ok <- is.finite(x)
    na <- data.table(n = sum(ok), mean_marapr_mm = NA_real_, mean_before = NA_real_, mean_after = NA_real_,
                     diff_mm = NA_real_, pct_change = NA_real_, welch_p = NA_real_, mannwhitney_p = NA_real_,
                     pettitt_p = NA_real_, pettitt_first_year_after = NA_integer_,
                     sen_full = NA_real_, sen_pre = NA_real_, sen_post = NA_real_)
    if (sum(ok) < 30L || stats::sd(x[ok]) <= 0) return(na)
    xx <- x[ok]; yy <- yrs_pix[ok]
    b <- xx[yy < SWE_STEP_YEAR]; a <- xx[yy >= SWE_STEP_YEAR]
    pt <- tryCatch(trend::pettitt.test(xx), error = function(e) NULL)
    sen <- function(v) .try_na(as.numeric(trend::sens.slope(v)$estimates))
    data.table(n = sum(ok), mean_marapr_mm = mean(xx), mean_before = mean(b), mean_after = mean(a),
               diff_mm = mean(a) - mean(b),
               pct_change = if (mean(b) > 0) 100 * (mean(a) - mean(b)) / mean(b) else NA_real_,
               welch_p = .try_na(stats::t.test(a, b)$p.value),
               mannwhitney_p = .try_na(stats::wilcox.test(a, b, exact = FALSE)$p.value),
               pettitt_p = if (is.null(pt)) NA_real_ else pt$p.value,
               pettitt_first_year_after = if (is.null(pt)) NA_integer_ else yy[min(as.integer(pt$estimate)[1] + 1L, length(yy))],
               sen_full = sen(xx), sen_pre = sen(b), sen_post = sen(a))
  })
  pix_tbl <- cbind(coords_dt[, .(space_idx, x, y)], rbindlist(pix_rows))
  pix_tbl[, `:=`(pettitt_p_fdr = stats::p.adjust(pettitt_p, "BH"),
                 welch_p_fdr   = stats::p.adjust(welch_p, "BH"))]
  fwrite(pix_tbl, file.path(swe_out_dir, "swe_pixel_step_tests.csv"))
  
  tryCatch({
    for (nm in c("pct_change", "pettitt_first_year_after")) {
      r <- terra::rast(analysis_template); terra::values(r) <- NA_real_
      r[valid_cell_indices] <- as.numeric(pix_tbl[[nm]])
      names(r) <- paste0("swe_marapr_", nm)
      terra::writeRaster(r, file.path(swe_out_dir, paste0("swe_marapr_", nm, ".tif")), overwrite = TRUE)
    }
  }, error = function(e) log_event(paste("  [SWE-HOMOG] raster export skipped:", conditionMessage(e))))
  
  gb <- function(s, col) step_tbl[region_id == NECHAKO_REGION_BASIN_ID & series == s, get(col)]
  eb <- function(en) era_tbl[region_id == NECHAKO_REGION_BASIN_ID & series == "SWE_MarApr" & era == en, mean_value]
  log_event(sprintf("  [SWE-HOMOG] Basin Mar-Apr SWE: mean 1950-78 = %.1f mm | 1979-90 = %.1f mm | 1979-2025 = %.1f mm | 1991-2020 = %.1f mm",
                    eb("1950-1978"), eb("1979-1990 (post-hoc window)"), eb("1979-2025"), eb("1991-2020")))
  log_event(sprintf("  [SWE-HOMOG] Basin Mar-Apr SWE: Pettitt p=%.4g (last yr before break %s) | SNHT p=%.4g | Buishand p=%.4g",
                    gb("SWE_MarApr", "pettitt_p"), gb("SWE_MarApr", "pettitt_last_year_before"),
                    gb("SWE_MarApr", "snht_p"), gb("SWE_MarApr", "buishand_range_p")))
  log_event(sprintf("  [SWE-HOMOG] Fixed split at %d: change=%.1f mm (%.1f%%), Cohen d=%.2f, Welch p=%.4g, Mann-Whitney p=%.4g",
                    SWE_STEP_YEAR, gb("SWE_MarApr", "diff_after_minus_before"), gb("SWE_MarApr", "pct_change"),
                    gb("SWE_MarApr", "cohens_d"), gb("SWE_MarApr", "welch_p"), gb("SWE_MarApr", "mannwhitney_p")))
  log_event(sprintf("  [SWE-HOMOG] Trend (mm/yr) without step = %.3f (p=%.4g); with step term = %.3f (p=%.4g); step = %.1f mm (p=%.4g); dAIC = %.1f (negative favours step)",
                    gb("SWE_MarApr", "trend_no_step"), gb("SWE_MarApr", "trend_no_step_p"),
                    gb("SWE_MarApr", "trend_with_step"), gb("SWE_MarApr", "trend_with_step_p"),
                    gb("SWE_MarApr", "step_coef"), gb("SWE_MarApr", "step_p"),
                    gb("SWE_MarApr", "delta_AIC_step_vs_nostep")))
  log_event(sprintf("  [SWE-HOMOG] Forcing check (basin): Nov-Mar precip Pettitt p=%.4g, change %.1f%% | Nov-Mar T Pettitt p=%.4g, change %.2f degC | SWE residual after P,T: Pettitt p=%.4g, step change %.1f mm (Welch p=%.4g)",
                    gb("Precip_NovMar", "pettitt_p"), gb("Precip_NovMar", "pct_change"),
                    gb("Temp_NovMar", "pettitt_p"), gb("Temp_NovMar", "diff_after_minus_before"),
                    gb("SWE_MarApr_resid_PT", "pettitt_p"), gb("SWE_MarApr_resid_PT", "diff_after_minus_before"),
                    gb("SWE_MarApr_resid_PT", "welch_p")))
  n_pix <- sum(is.finite(pix_tbl$pettitt_p))
  sp <- tryCatch(suppressWarnings(stats::cor.test(pix_tbl$pct_change, pix_tbl$mean_marapr_mm,
                                                  method = "spearman", exact = FALSE)), error = function(e) NULL)
  log_event(sprintf("  [SWE-HOMOG] Pixels: %d tested | Pettitt p<%.2f: %d (FDR-adjusted: %d) | negative 1979 step: %.1f%% | median change %.1f%% | Spearman(pct change, mean SWE) rho=%s",
                    n_pix, alpha, sum(pix_tbl$pettitt_p < alpha, na.rm = TRUE),
                    sum(pix_tbl$pettitt_p_fdr < alpha, na.rm = TRUE),
                    100 * mean(pix_tbl$diff_mm < 0, na.rm = TRUE), stats::median(pix_tbl$pct_change, na.rm = TRUE),
                    if (is.null(sp)) "NA" else sprintf("%.2f (p=%.3g)", sp$estimate, sp$p.value)))
  by_dec <- pix_tbl[is.finite(pettitt_p) & pettitt_p < alpha, .N, by = .(decade = 10L * (pettitt_first_year_after %/% 10L))][order(decade)]
  if (nrow(by_dec)) log_event(paste("  [SWE-HOMOG] Significant-pixel break decades:",
                                    paste(sprintf("%ds=%d", by_dec$decade, by_dec$N), collapse = ", ")))
  
  tryCatch({
    bb <- homog_long[region_id == NECHAKO_REGION_BASIN_ID]
    panels <- list(list("SWE_MarApr", "Mar-Apr SWE (mm)"), list("Precip_NovMar", "Nov-Mar precipitation (mm)"),
                   list("Temp_NovMar", "Nov-Mar temperature (degC)"), list("SWE_MarApr_resid_PT", "Mar-Apr SWE residual after P, T (mm)"))
    png(file.path(swe_out_dir, "swe_homogeneity_basin.png"), width = 1500, height = 1700, res = 150)
    op <- par(mfrow = c(4, 1), mar = c(3, 4.5, 2.2, 1))
    for (pn in panels) {
      d <- bb[series == pn[[1]]]
      plot(d$year, d$value, type = "b", pch = 16, cex = 0.6, xlab = "", ylab = pn[[2]],
           main = sprintf("%s  |  Pettitt p=%.3g", pn[[1]], gb(pn[[1]], "pettitt_p")))
      mb <- mean(d$value[d$year < SWE_STEP_YEAR], na.rm = TRUE); ma <- mean(d$value[d$year >= SWE_STEP_YEAR], na.rm = TRUE)
      segments(min(d$year), mb, SWE_STEP_YEAR - 0.5, mb, col = "firebrick", lwd = 2)
      segments(SWE_STEP_YEAR - 0.5, ma, max(d$year), ma, col = "steelblue", lwd = 2)
      abline(v = SWE_STEP_YEAR - 0.5, lty = 2, col = "grey40")
      pk <- gb(pn[[1]], "pettitt_last_year_before")
      if (is.finite(pk)) abline(v = pk + 0.5, lty = 3, col = "darkorange", lwd = 2)
      if (pn[[1]] == "SWE_MarApr") legend("topright", bty = "n", cex = 0.8, lty = c(1, 1, 2, 3), col = c("firebrick", "steelblue", "grey40", "darkorange"),
                                          legend = c("mean before step year", "mean from step year", "a-priori step year", "Pettitt break"))
    }
    par(op); invisible(dev.off())
  }, error = function(e) { try(dev.off(), silent = TRUE); log_event(paste("  [SWE-HOMOG] figure skipped:", conditionMessage(e))) })
  
  log_event(sprintf("  [SWE-HOMOG] Outputs written to %s", normalizePath(swe_out_dir, winslash = "/", mustWork = FALSE)))
} else {
  log_event("SWE file not found - SWE trend and homogeneity diagnostics skipped.")
}

####################################################################################
# MULTI-VARIABLE HOMOGENEITY / STEP-CHANGE DIAGNOSTICS (v9)
#
# Series tested (basin + 7 subbasins), each with the same annual metric as its trend test:
#   raw      PET_AprOct, PET_ERA5_AprOct, AET_annual, Runoff_annual, Snowfall_annual,
#            SnowDepth_MarApr, SnowDensity_MarApr, PM_<control>_annual, PM_<control>_AprOct
#   ratio    RunoffCoef, SnowfallFraction, AET_over_P, AET_over_PET_AprOct, PEV_over_PM_AprOct
#   residual Runoff_resid_P (~ P + P[t-1]), Snowfall_resid_PT (~ P + T), AET_resid_P_PET
#            (~ P + PET), PET_resid_PM and PET_ERA5_resid_PM (~ Apr-Oct PM controls)
# Tests: Pettitt, SNHT, Buishand range (unknown date); Welch / Mann-Whitney / Cohen d /
# Fligner at SWE_STEP_YEAR (a-priori); AR(1) GLS trend with/without a step; VC-MK/TFPW-MK on
# sub-periods (raw + ratio series). BH-FDR across regions within each series and across all
# basin series. Pixel-level Pettitt + fixed-date scans with BH-FDR and GeoTIFFs.
#
# Interpretation:
#   * All series are ERA5-Land products sharing the same atmospheric forcing. A break that
#     appears in a raw series AND its ratio/residual series points to the land-surface
#     product; a break that disappears in the ratio/residual series is inherited from the
#     forcing (or is a real climate shift). Separating forcing artefacts from real change
#     needs independent observations (WSC gauges, snow courses, stations).
#   * Pettitt / SNHT / Buishand assume independent years; positive autocorrelation inflates
#     false detections. The AR(1) step+trend model is the autocorrelation-aware check.
#   * A trend can trigger the unknown-date tests (break near the record midpoint); compare
#     with delta_AIC_step_vs_nostep and trend_with_step before calling a break.
####################################################################################
log_event("=== MULTI-VARIABLE HOMOGENEITY / STEP-CHANGE DIAGNOSTICS ===")
dir.create(homog_out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260102)
bid <- NECHAKO_REGION_BASIN_ID

.reg_annual <- function(rm_, col, use_months, fun = c("sum", "mean")) {
  fun <- match.arg(fun)
  yrs_all <- sort(unique(rm_$year))
  if (!col %in% names(rm_) || !any(is.finite(rm_[[col]]))) return(rep(NA_real_, length(yrs_all)))
  m <- annual_metric_fast(rm_[[col]], rm_$year, rm_$month, use_months, "mean")[, 1]
  if (fun == "sum") m * length(use_months) else m
}
.ratio <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b > 0, a / b, NA_real_)
.lag1  <- function(v) c(NA_real_, head(v, -1L))
.resid_on <- function(y, X) {
  d <- data.frame(y = as.numeric(y), X)
  ok <- stats::complete.cases(d)
  out <- rep(NA_real_, length(y)); r2 <- NA_real_
  if (sum(ok) >= 30L && stats::sd(d$y[ok]) > 0) {
    fit <- stats::lm(y ~ ., data = d[ok, , drop = FALSE])
    out[ok] <- stats::residuals(fit); r2 <- summary(fit)$r.squared
  }
  list(resid = out, r2 = r2)
}
.hs <- function(col, use_months, fun, family, is_precip, skip, pct_ok = TRUE)
  list(col = col, months = use_months, fun = fun, family = family,
       is_precip = is_precip, skip = skip, pct_ok = pct_ok)

HOMOG_SERIES <- list(
  PET_AprOct         = .hs("pet_mm_month",            warm_months_pet,   "sum",  "PET",         FALSE, FALSE),
  PET_ERA5_AprOct    = .hs("pet_era5_mm_month",       warm_months_pet,   "sum",  "PET",         FALSE, FALSE),
  AET_annual         = .hs("aet_mm_month",            1:12,              "sum",  "AET",         TRUE,  FALSE),
  Runoff_annual      = .hs("runoff_mm_month",         1:12,              "sum",  "Runoff",      TRUE,  FALSE),
  Snowfall_annual    = .hs("snowfall_mm_month",       1:12,              "sum",  "Snowfall",    TRUE,  FALSE),
  SnowDepth_MarApr   = .hs("snow_depth_cm_month",     SWE_ANNUAL_MONTHS, "mean", "SnowDepth",   TRUE,  FALSE),
  SnowDensity_MarApr = .hs("snow_density_kgm3_month", SWE_ANNUAL_MONTHS, "mean", "SnowDensity", FALSE, TRUE))
for (lab in unname(PM_CONTROL_LABELS)) {
  pct_ok <- lab != "PM_Temperature"   # a percent change of a temperature is meaningless
  HOMOG_SERIES[[paste0(lab, "_annual")]] <- .hs(lab, 1:12,            "mean", "PM_control", FALSE, TRUE, pct_ok)
  HOMOG_SERIES[[paste0(lab, "_AprOct")]] <- .hs(lab, warm_months_pet, "mean", "PM_control", FALSE, TRUE, pct_ok)
}
HOMOG_SERIES <- HOMOG_SERIES[vapply(HOMOG_SERIES, function(s)
  s$col %in% names(regional_monthly) && any(is.finite(regional_monthly[[s$col]])), logical(1))]
log_event(sprintf("  [HOMOG] %d raw series configured: %s", length(HOMOG_SERIES),
                  paste(names(HOMOG_SERIES), collapse = ", ")))

# ---- 1. Annual series per region (raw, ratio, residual) ----------------------------
homog_rows <- list()
for (id in unique(regional_monthly$region_id)) {
  rm_ <- regional_monthly[region_id == id][order(date)]
  yrs_all <- sort(unique(rm_$year))
  add <- function(series, value, family, kind, is_precip, skip, pct_ok, r2 = NA_real_) {
    homog_rows[[length(homog_rows) + 1L]] <<- data.table(
      region_id = id, year = yrs_all, series = series, family = family, kind = kind,
      value = as.numeric(value), is_precip = is_precip, skip_min_filter = skip,
      pct_ok = pct_ok, fit_r2 = r2)
  }
  raw <- list()
  for (sn in names(HOMOG_SERIES)) {
    s <- HOMOG_SERIES[[sn]]
    raw[[sn]] <- .reg_annual(rm_, s$col, s$months, s$fun)
    add(sn, raw[[sn]], s$family, "raw", s$is_precip, s$skip, s$pct_ok)
  }
  P  <- .reg_annual(rm_, "precip_mm_month", 1:12, "sum")
  Tm <- .reg_annual(rm_, "tair_degC_month", 1:12, "mean")
  
  if (!is.null(raw$Runoff_annual)) {
    add("RunoffCoef", .ratio(raw$Runoff_annual, P), "Runoff", "ratio", FALSE, TRUE, TRUE)
    rr <- .resid_on(raw$Runoff_annual, data.frame(P = P, P_lag1 = .lag1(P)))
    add("Runoff_resid_P", rr$resid, "Runoff", "residual", FALSE, TRUE, FALSE, rr$r2)
  }
  if (!is.null(raw$Snowfall_annual)) {
    add("SnowfallFraction", .ratio(raw$Snowfall_annual, P), "Snowfall", "ratio", FALSE, TRUE, TRUE)
    rr <- .resid_on(raw$Snowfall_annual, data.frame(P = P, T = Tm))
    add("Snowfall_resid_PT", rr$resid, "Snowfall", "residual", FALSE, TRUE, FALSE, rr$r2)
  }
  if (!is.null(raw$AET_annual)) {
    add("AET_over_P", .ratio(raw$AET_annual, P), "AET", "ratio", FALSE, TRUE, TRUE)
    if (!is.null(raw$PET_AprOct)) {
      aet_ao <- .reg_annual(rm_, "aet_mm_month", warm_months_pet, "sum")
      add("AET_over_PET_AprOct", .ratio(aet_ao, raw$PET_AprOct), "AET", "ratio", FALSE, TRUE, TRUE)
      rr <- .resid_on(raw$AET_annual, data.frame(P = P, PET = raw$PET_AprOct))
      add("AET_resid_P_PET", rr$resid, "AET", "residual", FALSE, TRUE, FALSE, rr$r2)
    }
  }
  if (all(c("PET_AprOct", "PET_ERA5_AprOct") %in% names(raw)))
    add("PEV_over_PM_AprOct", .ratio(raw$PET_ERA5_AprOct, raw$PET_AprOct), "PET", "ratio", FALSE, TRUE, TRUE)
  pm_ao <- intersect(paste0(unname(PM_CONTROL_LABELS), "_AprOct"), names(raw))
  if (length(pm_ao) >= 2L) {
    X <- as.data.frame(raw[pm_ao])
    for (tgt in intersect(c("PET_AprOct", "PET_ERA5_AprOct"), names(raw))) {
      rr <- .resid_on(raw[[tgt]], X)
      add(paste0(sub("_AprOct$", "", tgt), "_resid_PM"), rr$resid, "PET", "residual", FALSE, TRUE, FALSE, rr$r2)
    }
  }
}
homog_all <- rbindlist(homog_rows, fill = TRUE)
fwrite(homog_all, file.path(homog_out_dir, "homogeneity_annual_series_all.csv"))
log_event(sprintf("  [HOMOG] %d series x %d regions built (raw %d, ratio %d, residual %d)",
                  uniqueN(homog_all$series), uniqueN(homog_all$region_id),
                  uniqueN(homog_all[kind == "raw", series]), uniqueN(homog_all[kind == "ratio", series]),
                  uniqueN(homog_all[kind == "residual", series])))

# ---- 2. Step-change battery + step/trend model ---------------------------------------
log_event("  [HOMOG] Pettitt / SNHT / Buishand / fixed-date / AR(1) step+trend for every region x series...")
step_all <- rbindlist(lapply(split(homog_all, by = c("region_id", "series"), drop = TRUE), function(dd) {
  cbind(dd[1, .(region_id, series, family, kind, pct_ok)],
        homogeneity_tests(dd$value, dd$year, m_sim = N_SIM_HOMOG_ALL),
        step_trend_model(dd$value, dd$year))
}), fill = TRUE)
step_all[pct_ok == FALSE, pct_change := NA_real_]
step_all[, pct_ok := NULL]
step_all[, pettitt_p_fdr := stats::p.adjust(pettitt_p, "BH"), by = series]          # across regions
step_all[, step_p_fdr    := stats::p.adjust(step_p,    "BH"), by = series]
step_all[, pettitt_p_fdr_basin_all := NA_real_]
step_all[region_id == bid, pettitt_p_fdr_basin_all := stats::p.adjust(pettitt_p, "BH")]   # across basin series
step_all[, break_near_step_year := is.finite(pettitt_first_year_after) &
           abs(pettitt_first_year_after - SWE_STEP_YEAR) <= HOMOG_BREAK_WINDOW]
fwrite(step_all, file.path(homog_out_dir, "step_tests_by_region_all.csv"))

# ---- 3. Sub-period trends (raw + ratio series) ---------------------------------------
log_event("  [HOMOG] Sub-period VC-MK / TFPW-MK trends (raw and ratio series)...")
periods_h <- list(full = c(min(annual_years), max(annual_years)),
                  pre_step = c(min(annual_years), SWE_STEP_YEAR - 1L),
                  post_step = c(SWE_STEP_YEAR, max(annual_years)),
                  clim_1991_2020 = c(1991L, 2020L))
sub_all <- rbindlist(lapply(split(homog_all[kind %in% c("raw", "ratio")], by = c("region_id", "series"), drop = TRUE),
                            function(dd) subperiod_trends(dd$value, dd$year, dd$series[1], dd$region_id[1], periods_h,
                                                          is_precip = dd$is_precip[1],
                                                          skip_min_filter = dd$skip_min_filter[1])), fill = TRUE)
fwrite(sub_all, file.path(homog_out_dir, "subperiod_trends_all.csv"))

# ---- 4. Period means -------------------------------------------------------------------
era_def_h <- list(`1950-1978` = c(1950L, SWE_STEP_YEAR - 1L), `1979-1990 (post-hoc window)` = c(SWE_STEP_YEAR, 1990L),
                  `1979-2025` = c(SWE_STEP_YEAR, max(annual_years)), `1991-2020` = c(1991L, 2020L),
                  `2021-2025` = c(2021L, max(annual_years)))
era_all <- rbindlist(lapply(names(era_def_h), function(en) {
  homog_all[year >= era_def_h[[en]][1] & year <= era_def_h[[en]][2],
            .(era = en, n_years = sum(is.finite(value)), mean_value = mean(value, na.rm = TRUE),
              sd_value = stats::sd(value, na.rm = TRUE)), by = .(region_id, series, kind)]
}))
fwrite(era_all, file.path(homog_out_dir, "period_means_all.csv"))

# ---- 5. Cross-variable break summary (basin), including the SWE series ------------------
brk_cols <- c("region_id", "series", "family", "kind", "n", "pettitt_p", "pettitt_p_fdr",
              "pettitt_p_fdr_basin_all", "pettitt_first_year_after", "snht_p", "snht_first_year_after",
              "buishand_range_p", "diff_after_minus_before", "pct_change", "cohens_d", "welch_p",
              "trend_no_step", "trend_with_step", "step_coef", "step_p", "delta_AIC_step_vs_nostep",
              "break_near_step_year")
brk <- step_all[region_id == bid, intersect(brk_cols, names(step_all)), with = FALSE]
if (RUN_SWE && exists("step_tbl")) {
  s <- data.table::copy(step_tbl[region_id == bid])
  s[, `:=`(family = "SWE", kind = fifelse(grepl("resid", series), "residual", "raw"),
           pettitt_p_fdr_basin_all = NA_real_,
           break_near_step_year = is.finite(pettitt_first_year_after) &
             abs(pettitt_first_year_after - SWE_STEP_YEAR) <= HOMOG_BREAK_WINDOW)]
  brk <- rbind(brk, s[, intersect(brk_cols, names(s)), with = FALSE], fill = TRUE)
}
setorder(brk, family, kind, series)
fwrite(brk, file.path(homog_out_dir, "homogeneity_break_summary_basin.csv"))
.sig <- step_all[region_id == bid & is.finite(pettitt_p_fdr_basin_all) & pettitt_p_fdr_basin_all < alpha]
log_event(sprintf("  [HOMOG] Basin: %d of %d series have a Pettitt break (BH across basin series, p < %.2f); %d within +/-%d yr of %d%s",
                  nrow(.sig), nrow(step_all[region_id == bid]), alpha,
                  sum(.sig$break_near_step_year), HOMOG_BREAK_WINDOW, SWE_STEP_YEAR,
                  if (nrow(.sig)) paste0(": ", paste(sprintf("%s (%s)", .sig$series, .sig$pettitt_first_year_after), collapse = ", ")) else ""))

# ---- 6. Pixel-level scans ---------------------------------------------------------------
log_event("  [HOMOG] Pixel-level Pettitt + fixed-date scans...")
.pixel_step_scan <- function(ann_mat, positive = TRUE) {
  yrs <- as.integer(rownames(ann_mat))
  rbindlist(lapply(seq_len(ncol(ann_mat)), function(i) {
    x <- ann_mat[, i]; ok <- is.finite(x)
    na <- data.table(n = sum(ok), mean_value = NA_real_, mean_before = NA_real_, mean_after = NA_real_,
                     diff = NA_real_, pct_change = NA_real_, welch_p = NA_real_, mannwhitney_p = NA_real_,
                     pettitt_p = NA_real_, pettitt_first_year_after = NA_integer_, sen_full = NA_real_)
    if (sum(ok) < 30L || !isTRUE(stats::sd(x[ok]) > 0)) return(na)
    xx <- x[ok]; yy <- yrs[ok]
    b <- xx[yy < SWE_STEP_YEAR]; a <- xx[yy >= SWE_STEP_YEAR]
    if (length(b) < 5L || length(a) < 5L) return(na)
    pt <- tryCatch(trend::pettitt.test(xx), error = function(e) NULL)
    data.table(n = sum(ok), mean_value = mean(xx), mean_before = mean(b), mean_after = mean(a),
               diff = mean(a) - mean(b),
               pct_change = if (positive && mean(b) > 0) 100 * (mean(a) - mean(b)) / mean(b) else NA_real_,
               welch_p = .try_na(stats::t.test(a, b)$p.value),
               mannwhitney_p = .try_na(stats::wilcox.test(a, b, exact = FALSE)$p.value),
               pettitt_p = if (is.null(pt)) NA_real_ else pt$p.value,
               pettitt_first_year_after = if (is.null(pt)) NA_integer_ else
                 yy[min(as.integer(pt$estimate)[1] + 1L, length(yy))],
               sen_full = .try_na(as.numeric(trend::sens.slope(xx)$estimates)))
  }))
}

pix_mats <- list(
  PET_AprOct      = list(m = aggregate_to_annual_warmseason_fast(pet_matrix, years, months, warm_months_pet, "sum"), positive = TRUE),
  PET_ERA5_AprOct = list(m = aggregate_to_annual_warmseason_fast(pet_era5_matrix, years, months, warm_months_pet, "sum"), positive = TRUE),
  AET_annual      = list(m = aet_annual_pix, positive = TRUE))
for (nm in names(extra_annual_pix))
  pix_mats[[paste0(nm, if (EXTRA_VARS[[nm]]$annual == "sum") "_annual" else "_MarApr")]] <-
  list(m = extra_annual_pix[[nm]], positive = TRUE)
for (nm in names(pm_control_matrices)) {
  lab <- PM_CONTROL_LABELS[[nm]]; pos <- lab != "PM_Temperature"
  pix_mats[[paste0(lab, "_annual")]] <- list(m = aggregate_to_annual_fast(pm_control_matrices[[nm]], years, "mean"), positive = pos)
  pix_mats[[paste0(lab, "_AprOct")]] <- list(m = aggregate_to_annual_warmseason_fast(pm_control_matrices[[nm]], years, months,
                                                                                     warm_months_pet, "mean"), positive = pos)
}

pix_dir <- file.path(homog_out_dir, "pixel")
dir.create(pix_dir, recursive = TRUE, showWarnings = FALSE)
pix_summary <- rbindlist(lapply(names(pix_mats), function(lab) {
  tb <- cbind(coords_dt[, .(space_idx, x, y)], .pixel_step_scan(pix_mats[[lab]]$m, pix_mats[[lab]]$positive))
  tb[, `:=`(pettitt_p_fdr = stats::p.adjust(pettitt_p, "BH"), welch_p_fdr = stats::p.adjust(welch_p, "BH"))]
  fwrite(tb, file.path(pix_dir, sprintf("pixel_step_tests_%s.csv", lab)))
  tryCatch({
    lyr <- c("diff", "pct_change", "pettitt_first_year_after", "pettitt_sig_fdr")
    vals <- cbind(tb$diff, tb$pct_change, tb$pettitt_first_year_after,
                  as.numeric(is.finite(tb$pettitt_p_fdr) & tb$pettitt_p_fdr < alpha))
    mm <- matrix(NA_real_, terra::ncell(analysis_template), length(lyr))
    mm[valid_cell_indices, ] <- vals
    r <- terra::rast(analysis_template, nlyrs = length(lyr))
    terra::values(r) <- mm
    names(r) <- paste0(lab, "_", lyr)
    terra::writeRaster(r, file.path(pix_dir, sprintf("pixel_step_%s.tif", lab)), overwrite = TRUE)
  }, error = function(e) log_event(sprintf("  [HOMOG] raster export skipped for %s: %s", lab, conditionMessage(e))))
  sig <- tb[is.finite(pettitt_p_fdr) & pettitt_p_fdr < alpha]
  out <- data.table(
    series = lab, n_pixels_tested = sum(is.finite(tb$pettitt_p)),
    pct_pettitt_p05   = 100 * mean(tb$pettitt_p < alpha, na.rm = TRUE),
    pct_pettitt_fdr05 = 100 * mean(tb$pettitt_p_fdr < alpha, na.rm = TRUE),
    pct_negative_step = 100 * mean(tb$diff < 0, na.rm = TRUE),
    median_diff = stats::median(tb$diff, na.rm = TRUE),
    median_pct_change = suppressWarnings(stats::median(tb$pct_change, na.rm = TRUE)),
    median_break_year_fdr_sig = if (nrow(sig)) stats::median(sig$pettitt_first_year_after) else NA_real_,
    pct_fdr_sig_breaks_near_step_year = if (nrow(sig)) 100 * mean(abs(sig$pettitt_first_year_after - SWE_STEP_YEAR) <= HOMOG_BREAK_WINDOW) else NA_real_)
  log_event(sprintf("    %-28s tested %4d | Pettitt FDR<%.2f: %5.1f%% | negative step %5.1f%% | median FDR-sig break year %s",
                    lab, out$n_pixels_tested, alpha, out$pct_pettitt_fdr05, out$pct_negative_step,
                    if (is.finite(out$median_break_year_fdr_sig)) sprintf("%.0f", out$median_break_year_fdr_sig) else "-"))
  out
}), fill = TRUE)
fwrite(pix_summary, file.path(homog_out_dir, "pixel_step_summary.csv"))

# ---- 7. Diagnostic figure (basin, all series) ----------------------------------------
tryCatch({
  bb  <- homog_all[region_id == bid & is.finite(value)]
  ser <- unique(bb$series); nser <- length(ser)
  nc  <- 4L; nr <- ceiling(nser / nc)
  png(file.path(homog_out_dir, "homogeneity_basin_all_series.png"), width = 2400, height = 420 * nr, res = 150)
  op <- par(mfrow = c(nr, nc), mar = c(2.5, 4, 2.4, 0.8))
  for (s in ser) {
    d  <- bb[series == s][order(year)]
    st <- step_all[region_id == bid & series == s]
    plot(d$year, d$value, type = "b", pch = 16, cex = 0.45, xlab = "", ylab = "",
         main = sprintf("%s | Pettitt p=%.3g", s, st$pettitt_p), cex.main = 0.85)
    mb <- mean(d$value[d$year < SWE_STEP_YEAR], na.rm = TRUE); ma <- mean(d$value[d$year >= SWE_STEP_YEAR], na.rm = TRUE)
    segments(min(d$year), mb, SWE_STEP_YEAR - 0.5, mb, col = "firebrick", lwd = 2)
    segments(SWE_STEP_YEAR - 0.5, ma, max(d$year), ma, col = "steelblue", lwd = 2)
    abline(v = SWE_STEP_YEAR - 0.5, lty = 2, col = "grey40")
    pk <- st$pettitt_first_year_after
    if (length(pk) && is.finite(pk)) abline(v = pk - 0.5, lty = 3, col = "darkorange", lwd = 2)
  }
  par(op); invisible(dev.off())
}, error = function(e) { try(dev.off(), silent = TRUE); log_event(paste("  [HOMOG] figure skipped:", conditionMessage(e))) })

metadata <- readRDS(metadata_path)
metadata$homogeneity <- list(dir = homog_out_dir, step_year = SWE_STEP_YEAR, n_sim = N_SIM_HOMOG_ALL,
                             break_window = HOMOG_BREAK_WINDOW, series = unique(homog_all[, .(series, family, kind)]))
saveRDS(metadata, metadata_path, compress = TRUE)

log_event(sprintf("  [HOMOG] Outputs written to %s", normalizePath(homog_out_dir, winslash = "/", mustWork = FALSE)))
log_event("DATA PROCESSING COMPLETE")
cat("\n[OK] DATA PROCESSING COMPLETED SUCCESSFULLY\n")