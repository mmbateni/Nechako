# ==============================================================================
#   13trend_test_ALL.R 
#   SPATIAL TREND, EVENT AND HOMOGENEITY ANALYSIS FOR ALL DROUGHT INDICES
# ==============================================================================
# INDICES : SPI, SPEI-PM and SRI (1, 3, 6, 12, 24, 36 months), SWEI-3,
#           SSPI-1, SDI (hydrological-year periods k1-k4, Oct start), and the
#           soil-moisture indices SSMI and SSI-SM (1 and 3 months; layer
#           SOIL_LAYER, written by 8_SSI_SSMI_ERALand.R).
#
# GRID CONTRACT (unchanged): every index x scale uses the native grid of its own
# NetCDF files; the Nechako polygon is rasterized on that grid; fractional coverage
# x cell area gives the weights. No canonical-grid RDS, no resampling/reprojection.
#
# STATISTICAL DESIGN
#   Standardized indices have no seasonal cycle in the mean, but a k-month
#   accumulation makes consecutive monthly values share k-1 months: lag-1
#   autocorrelation ~ (k-1)/k (SPI-12 ~0.92, SPI-36 ~0.97). The same calendar month
#   in consecutive years shares max(0, k-12) months (~0 for k <= 12).
#   A. Full monthly series (n ~ 912): Hamed-Rao VC-MK = primary. TFPW-MK reported
#      but its AR(1) prewhitening mis-specifies the MA(k-1) dependence for k >= 3
#      (REJECTED for inference at k >= 3; values retained only as a model-failure
#      diagnostic). Basin series: detrended-null block-bootstrap MK, block =
#      max(12, 2k) months, as an assumption-light check.
#   B. Calendar-month series (n ~ 76 per month) uses the same monthly analysis as Script 11:
#      VC-MK + TFPW-MK per pixel, BH-FDR within index x scale x month.
#   C. Annual drought-month fraction (months < DROUGHT_ONSET) uses complete
#      12-month calendar years for non-snow indices; SWEI uses Jan-Apr and SSPI
#      uses Dec-Mar (winter labelled by its ending year), with all 4 months required.
#      VC-MK per pixel; >12-month accumulation scales are exploratory only.
#   D. Events: ONE framework-wide definition from DROUGHT_ANALYSIS_utils.R
#      (threshold run method; hysteresis only if DROUGHT_END > DROUGHT_ONSET;
#      when equal, there is no hysteresis. Onset < DROUGHT_ONSET;
#      termination >= DROUGHT_END; minimum duration comes from
#      EVENT_MIN_DURATION_BY_SCALE). Event counts,
#      duration classes, the runs test and the basin event catalogue use these
#      events. PELT (on annual means) and the red-noise spectra use the raw index
#      series; basin drought extent is the monthly area share below DROUGHT_ONSET
#      (no minimum duration).
#   E. Basin and regional series: A + B (no gap splicing). The whole-basin row of
#      the regional tables reuses the basin tests (same series; not retested).
#   F. Shared annual-resolution homogeneity battery in utils (Pettitt, SNHT, Buishand, Welch at
#      the exploratory fixed-step year, AR(1)-GLS trend vs trend+step, sub-period trends,
#      pixel Pettitt scans with BH-FDR) ONLY on non-overlapping annual-resolution
#      series (HOMOG_SPECS): these tests assume independent observations.
#   G. Series QA (same no-gap principle and tie rule as Script 11):
#      full-monthly n >= 120; annual/calendar-month n >= ANALYSIS_MIN_HOMOG_YEARS;
#      filter largest tied group > TIE_MODE_MAX (= ANALYSIS_TIE_MODE_MAX, the same
#      rule Script 11 applies). Sub-period and Script 11 tests use the shared
#      minimum of ANALYSIS_MIN_TREND_N observations.
#      Filtered series: NA statistics + reason code. SWEI gaps are never spliced.
#   H. Point trends at Script 11's specific points.
#
# OUTPUTS (TREND_DIR)
#   {idx}_{sc}_results.csv                 full-series pixel results (+ events, QA, dfreq)
#   {idx}_{sc}_calmonth_results.csv        calendar-month pixel results
#   sdi_{P}_results.csv                    SDI pixel results (annual)
#   basin_trends/ regional_trends/ regional_timeseries/ homogeneity/ points/ qa/
#   duration-class CSVs, D13+ NetCDFs, basin_extent/, script13_run_status.csv,
#   script13_analysis_config.rds
#
# EXECUTION ORDER: index scripts 1, 2b, 3, 4, 5 (+8 soil; +6/7 optional) -> 11 -> 12 -> 13 -> 14 -> 15 -> 16
# ==============================================================================


# ==============================================================================
#   0. ENVIRONMENT AND USER SWITCHES
# ==============================================================================
# WD_PATH / UTILS_FILE are set from getwd() at the top of this script.
rm(list = ls())

# ---- Working directory (framework rule) -------------------------------------
# Run every script from the project root: set the working directory there
# BEFORE running (setwd(), or open the RStudio project). All paths below are
# derived from getwd(); no path is hard-coded.
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")

source(UTILS_FILE)
# Shared analysis contracts are defined in DROUGHT_ANALYSIS_utils.R.
required_utils <- c("ANALYSIS_STEP_YEAR", "ANALYSIS_TIE_MODE_MAX", "ANALYSIS_MIN_TREND_N",
                    "analysis_modal_tie_fraction", ".analysis_kendall_S",
                    "SPI_SCALES", "SPEI_SCALES", "SRI_SCALES",
                    "SWEI_SCALE", "SSPI_SCALE")
missing_utils <- required_utils[!vapply(required_utils, exists, logical(1), inherits = TRUE)]
if (length(missing_utils))
  stop("DROUGHT_ANALYSIS_utils.R lacks: ", paste(missing_utils, collapse = ", "),
       "\nApply the utils patches (U1-U7) before running this script.")
analysis_check_output_isolation(file.path(WD_PATH, "trend_analysis_pr_pet"), TREND_DIR)

utils_load_packages(c("terra", "data.table", "Kendall", "trend", "modifiedmk", "lubridate",
                      "changepoint", "zoo", "ggplot2", "dplyr", "patchwork", "nlme"))

# The shared MK wrappers in DROUGHT_ANALYSIS_utils.R can use fork workers on
# Unix and a pre-initialized PSOCK cluster on Windows. Keep index processing
# serial to avoid concurrent writes to the same result files.
.detected_cores_13 <- suppressWarnings(parallel::detectCores(logical = TRUE))
N_CORES <- if (length(.detected_cores_13) != 1L ||
               !is.finite(.detected_cores_13)) 1L else
                 max(1L, min(as.integer(.detected_cores_13) - 1L, 8L))
is_windows <- identical(.Platform$OS.type, "windows")
cl <- NULL
if (N_CORES > 1L && is_windows) {
  # On Windows, helper functions must exist in each worker session. If cluster
  # initialization fails, proceed serially rather than producing partial CSVs.
  cl <- tryCatch(parallel::makeCluster(N_CORES), error = function(e) {
    warning("Script 13 parallel cluster unavailable: ", conditionMessage(e),
            "; falling back to serial execution.", call. = FALSE)
    NULL
  })
  if (!is.null(cl)) {
    initialized <- tryCatch({
      parallel::clusterCall(cl, function(wd, uf) {
        setwd(wd)
        assign("WD_PATH", wd, envir = .GlobalEnv)
        assign("UTILS_FILE", uf, envir = .GlobalEnv)
        source(uf, local = FALSE)
        invisible(TRUE)
      }, WD_PATH, UTILS_FILE)
      TRUE
    }, error = function(e) {
      warning("Script 13 worker utility initialization failed: ", conditionMessage(e),
              "; falling back to serial execution.", call. = FALSE)
      FALSE
    })
    if (!initialized) {
      try(parallel::stopCluster(cl), silent = TRUE)
      cl <- NULL
      N_CORES <- 1L
    }
  } else N_CORES <- 1L
}
cat(sprintf("  Script 13 MK workers: %d (%s)\n", N_CORES,
            if (N_CORES == 1L) "serial" else if (is_windows) "PSOCK" else "fork"))


RUN_CALENDAR_MONTH  <- TRUE     # B
RUN_DROUGHT_FREQ    <- TRUE     # C
RUN_HOMOGENEITY     <- TRUE     # F
RUN_BLOCK_BOOTSTRAP <- TRUE     # A (basin series)
VERIFY_BASIN_SERIES <- TRUE     # recompute from pixels and stop if a producer's file differs
SAVE_TS_RDS         <- TRUE     # pixel time series for 16event_ranking.R 
MBB_B               <- 1000L    # block-bootstrap replicates
HOMOG_MC            <- ANALYSIS_HOMOG_MC    # Monte Carlo replicates for SNHT / Buishand p-values
MIN_N_MONTHLY       <- 120L     # min finite months for a full-series test
MIN_N_YEARS         <- ANALYSIS_MIN_HOMOG_YEARS      # >=30: more conservative than Script 11 (>=12); do not conflate tested domains
TIE_MODE_MAX        <- ANALYSIS_TIE_MODE_MAX   # shared with Script 11 (utils)
DFREQ_MIN_MONTHS    <- 12L      # non-snow: complete calendar years ONLY (fixed denominator)
DFREQ_SNOW_MONTHS   <- list(swei = 1:4, sspi = c(12L, 1L, 2L, 3L))
DFREQ_DEFINITION_VERSION <- "full_year_12_snow_4_end_year_v1"
# Separate full-record snow-index trends from physically distinct seasons.
# Fixed calendar windows, not an inferred snow mask. Each estimate uses exactly
# one non-overlapping seasonal mean per year; incomplete seasons are NA.
SNOW_BASIN_WINDOWS <- list(
  swei = list(accumulation = 1:4, melt = 5:9),
  sspi = list(accumulation = c(12L, 1L, 2L, 3L), melt = 5:9))
SNOW_SEASON_MBB_YEARS <- 5L
MIN_REGIONAL_VALID_AREA_COVERAGE <- BASIN_MIN_VALID_AREA_COVERAGE   # framework-wide rule (0.20), see utils

BASIN_TREND_DIR    <- file.path(TREND_DIR, "basin_trends")
REGIONAL_TREND_DIR <- file.path(TREND_DIR, "regional_trends")
REGIONAL_INDEX_DIR <- file.path(TREND_DIR, "regional_timeseries")
HOMOG_DIR          <- file.path(TREND_DIR, "homogeneity")
POINT_DIR          <- file.path(TREND_DIR, "points")
QA_DIR             <- file.path(TREND_DIR, "qa")
for (d in c(TREND_DIR, BASIN_TREND_DIR, REGIONAL_TREND_DIR, REGIONAL_INDEX_DIR, HOMOG_DIR,
            file.path(HOMOG_DIR, "pixel"), POINT_DIR, QA_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Specific points come from Script 11's metadata. The step year is ANALYSIS_STEP_YEAR
# (utils); Script 11's metadata is only checked, so a stale Script 11 run stops here.
S4_META_FILE <- file.path(WD_PATH, "trend_analysis_pr_pet", "analysis_metadata.rds")
S4_META <- if (file.exists(S4_META_FILE)) tryCatch(readRDS(S4_META_FILE), error = function(e) NULL) else NULL
STEP_YEAR <- ANALYSIS_STEP_YEAR
local({
  sy <- tryCatch(S4_META$step_year %||% S4_META$swe_config$step_year, error = function(e) NULL)
  if (length(sy) == 1L && is.finite(sy) && as.integer(sy) != STEP_YEAR)
    stop("Script 11 outputs were produced with step year ", sy, " but ANALYSIS_STEP_YEAR is ",
         STEP_YEAR, ". Re-run 11pr_pet_trends.r.")
})

.read_specific_points <- function(meta) {
  p <- tryCatch(meta$specific_pts, error = function(e) NULL)
  if (is.null(p) && exists("SPECIFIC_PTS")) p <- get("SPECIFIC_PTS")
  if (is.null(p)) return(NULL)
  p <- tryCatch(data.table::as.data.table(p), error = function(e) NULL)
  if (is.null(p) || !nrow(p)) return(NULL)
  nm <- tolower(names(p))
  lonc <- names(p)[match(TRUE, nm %in% c("lon_wgs84", "lon", "longitude", "x"))]
  latc <- names(p)[match(TRUE, nm %in% c("lat_wgs84", "lat", "latitude", "y"))]
  if (is.na(lonc) || is.na(latc)) return(NULL)
  idc <- names(p)[match(TRUE, nm %in% c("point_id", "id"))]
  data.table::data.table(point_id = if (is.na(idc)) seq_len(nrow(p)) else as.integer(p[[idc]]),
                         lon = as.numeric(p[[lonc]]), lat = as.numeric(p[[latc]]))
}
SPECIFIC_PTS_DT <- .read_specific_points(S4_META)
BASIN_RID <- if (exists("NECHAKO_REGION_BASIN_ID")) NECHAKO_REGION_BASIN_ID else "Nechako_Basin"

cat("\n============================================================\n",
    "SCRIPT 13 — NATIVE ERA5-LAND GRID | step year ", STEP_YEAR,
    " | specific points: ", if (is.null(SPECIFIC_PTS_DT)) 0L else nrow(SPECIFIC_PTS_DT),
    "\n============================================================\n", sep = "")


# ==============================================================================
#   1. NATIVE-GRID HELPERS AND CACHES
# ==============================================================================
assert_same_native_grid <- function(x, template, label = "Raster") {
  if (!inherits(x, "SpatRaster") || !inherits(template, "SpatRaster")) stop(label, ": not a SpatRaster.")
  geom_ok <- terra::compareGeom(x[[1]], template[[1]], stopOnError = FALSE,
                                crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)
  orig_ok <- tryCatch(isTRUE(all.equal(terra::origin(x[[1]]), terra::origin(template[[1]]), tolerance = 1e-6)),
                      error = function(e) FALSE)
  if (!isTRUE(geom_ok) || !isTRUE(orig_ok))
    stop(label, " does not match the native grid of the first file. Script never resamples.")
  invisible(TRUE)
}

BASIN_GRID_CACHE <- new.env(parent = emptyenv())
get_basin_grid_cached <- function(r_tmpl) {
  key <- .bs_grid_key(r_tmpl)
  if (!is.null(BASIN_GRID_CACHE[[key]])) return(BASIN_GRID_CACHE[[key]])
  bw <- basin_fractional_weights(r_tmpl)
  if (length(bw$coverage) != terra::ncell(r_tmpl)) stop("Basin weight/grid mismatch")
  out <- list(cover_all = bw$coverage, area_all = bw$cell_area_m2)
  BASIN_GRID_CACHE[[key]] <- out
  out
}

# ==============================================================================
#   2. INDICES, SCALES, DIRECTORIES, HOMOGENEITY SERIES
# ==============================================================================
# Non-soil products must match the producer precipitation raster. Load only once.
PRECIP_TEMPLATE_13 <- terra::rast(file.path(WD_PATH, "monthly_data_direct",
                                            "total_precipitation_monthly.nc"))[[1L]]
# SPI_SCALES, SPEI_SCALES, SRI_SCALES, SWEI_SCALE and SSPI_SCALE
# are defined once in DROUGHT_ANALYSIS_utils.R.
SDI_PERIODS <- c(3L, 6L, 9L, 12L)     # cumulative months from October: k1..k4
# Water years are labelled by their ENDING year (WY Y = Oct Y-1 to Sep Y), so k1
# (Oct-Dec) of the last calendar year carries the next label (Oct-Dec 2025 = WY 2026).
# TRUE keeps only water years ending by ANALYSIS_PERIOD_END: k1-k4 then share
# WY 1951-2025 and every label stays inside the record. FALSE keeps WY 2026 for k1.
SDI_LIMIT_TO_RECORD_END <- TRUE
SWEI_SEAS_DIR      <- file.path(WD_PATH, "swei_results_seasonal")
SSPI_SEAS_DIR      <- file.path(WD_PATH, "sspi_results_seasonal")
SRI_SEAS_DIR       <- file.path(WD_PATH, "sri_results_seasonal", "runoff")   # total runoff (ro)
SDI_DIR            <- file.path(WD_PATH, "sdi_results", "runoff")

# Soil-moisture indices from 8_SSI_SSMI_ERALand.R. Index key "ssi_sm" (not "ssi")
# because "ssi" is the streamflow index of Script 6 / Script 12.
RUN_SOIL_MOISTURE <- TRUE
SOIL_LAYER        <- "L1_3"          # "L1_3" (0-100 cm) or "L1_4" (0-289 cm)
SSMI_SEAS_DIR     <- file.path(WD_PATH, "ssmi_results_seasonal", SOIL_LAYER)
SSI_SM_SEAS_DIR   <- file.path(WD_PATH, "ssi_sm_results_seasonal", SOIL_LAYER)
BASIN_SERIES_DIRS <- list(spi = SPI_SEAS_DIR, spei = SPEI_SEAS_DIR,
                          swei = SWEI_SEAS_DIR, sspi = SSPI_SEAS_DIR, sri = SRI_SEAS_DIR,
                          ssmi = SSMI_SEAS_DIR, ssi_sm = SSI_SM_SEAS_DIR)
# Scales without a framework minimum event duration (EVENT_MIN_DURATION_BY_SCALE)
# are dropped here rather than failing inside process_index(); currently scale 2
# is unsupported, so only soil scales 1 and 3 are analyzed.
.soil_scale_ok <- function(s) !is.null(tryCatch(get_event_min_duration(s), error = function(e) NULL))
SSMI_SCALES <- SSI_SM_SCALES <- Filter(.soil_scale_ok, c(1L, 2L, 3L))
if (RUN_SOIL_MOISTURE && !identical(SSMI_SCALES, c(1L, 2L, 3L)))
  cat("  NOTE: soil-moisture scales limited to", paste(SSMI_SCALES, collapse = ", "),
      "(no EVENT_MIN_DURATION_BY_SCALE entry for the others)\n")

INDEX_LABEL  <- c(spi = "SPI", spei = "SPEI-PM", swei = "SWEI",
                  sspi = "SSPI", sri = "SRI", sdi = "SDI", ssmi = "SSMI", ssi_sm = "SSI-SM")
INDEX_FAMILY <- c(spi = "meteorological", spei = "meteorological",
                  swei = "snow", sspi = "snow", sri = "runoff", sdi = "runoff",
                  ssmi = "soil_moisture", ssi_sm = "soil_moisture")

# Non-overlapping annual-resolution series for the homogeneity battery.
.hs <- function(idx, sc, months, suffix, desc)
  data.table::data.table(index_type = idx, scale = as.integer(sc), months = list(as.integer(months)),
                         suffix = suffix, desc = desc)
HOMOG_SPECS <- data.table::rbindlist(c(
  unlist(lapply(c("spi", "spei", "sri"), function(i) list(
    .hs(i, 12, 12, "Dec", "Jan-Dec (calendar year)"),
    .hs(i, 12,  9, "Sep", "Oct-Sep (water year)"),
    .hs(i,  6,  9, "Sep", "Apr-Sep (warm season)"),
    .hs(i,  3,  8, "Aug", "Jun-Aug (summer)"))), recursive = FALSE),
  list(.hs("swei", 3, 4, "Apr", "Feb-Apr snowpack"),
       .hs("sspi", 1, c(11, 12, 1, 2, 3), "NovMar", "Nov-Mar mean"),
       .hs("sdi", 12, 0, "WY", "Oct-Sep water-year runoff volume"))))

# Minimum event durations come from get_event_min_duration() in
# DROUGHT_ANALYSIS_utils.R (EVENT_MIN_DURATION_BY_SCALE); edit them there.


# ==============================================================================
#   3. SERIES QA AND TREND HELPERS
# ==============================================================================
.bh <- function(p) {
  out <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (sum(ok) > 1L) out[ok] <- stats::p.adjust(p[ok], method = "BH")
  out
}
.trim_finite <- function(x) { f <- which(is.finite(x)); if (!length(f)) x[0] else x[min(f):max(f)] }
.modal_frac  <- analysis_modal_tie_fraction
.qa_reason <- function(x, min_n) {
  xf <- x[is.finite(x)]
  if (length(xf) < min_n) return("too_few_values")
  if (anyNA(.trim_finite(x))) return("internal_gaps")
  if (length(unique(xf)) < 3L) return("constant")
  if (.modal_frac(xf) > TIE_MODE_MAX) return("ties")
  "ok"
}
.qa_rows         <- function(mat, min_n) apply(mat, 1L, .qa_reason, min_n = min_n)
.modal_frac_rows <- function(mat) apply(mat, 1L, .modal_frac)
.rho1_one <- function(x) {
  x <- .trim_finite(x); n <- length(x)
  if (n < 10L || anyNA(x)) return(NA_real_)
  t <- seq_len(n) - (n + 1) / 2; xc <- x - mean(x)
  e <- xc - t * sum(t * xc) / sum(t^2)
  suppressWarnings(stats::cor(e[-1L], e[-n]))
}
.rho1_rows <- function(mat) apply(mat, 1L, .rho1_one)
.mk_na <- function(n) data.table::data.table(tau_vc = rep(NA_real_, n), p_value_vc = NA_real_, sl_vc = NA_real_,
                                             vc_status = NA_character_, vc_n_over_nstar = NA_real_,
                                             tau_tfpw = NA_real_, p_value_tfpw = NA_real_, sl_tfpw = NA_real_)
MK_REQUIRED <- c("tau_vc", "p_value_vc", "tau_tfpw", "p_value_tfpw", "sl_tfpw")

# Pixel MK battery (shared VC-HR98 and TFPW implementations) on the rows that pass QA.
.mk_pixel <- function(mat, valid) {
  n <- nrow(mat); vi <- which(valid)
  if (!length(vi)) return(.mk_na(n))
  sub <- mat[vi, , drop = FALSE]
  r <- cbind(data.table::as.data.table(vectorized_mann_kendall(sub)),
             data.table::as.data.table(vectorized_mann_kendall_tfpw(sub)))
  data.table::setnames(r, make.unique(names(r)))
  miss <- setdiff(MK_REQUIRED, names(r))
  if (length(miss)) stop("Shared MK utilities did not return column(s): ", paste(miss, collapse = ", "))
  r[match(seq_len(n), vi)]
}

# Single-series VC-HR98 (basin / regional / sub-period).
.mk_series <- function(x, min_n) {
  reason <- .qa_reason(x, min_n)
  out <- data.table::data.table(n = sum(is.finite(x)), qa = reason, tau = NA_real_, p_vc = NA_real_,
                                sen = NA_real_, n_over_neff = NA_real_, n_over_neff_raw = NA_real_,
                                p_mk_plain = NA_real_, rho1 = NA_real_, vc_status = NA_character_)
  if (reason != "ok") return(out)
  v  <- .trim_finite(x)
  vc <- vc_mmkh(v)
  # Shared vc_mmkh() may expose the uncapped variance ratio. Any ratio below
  # one must use the PLAIN MK variance (N/N* = 1), never the deflated VC variance.
  # Recompute the plain p-value if the utility does not expose it; a missing
  # p_uncorrected must never silently leave the anti-conservative p in place.
  raw_ratio <- as.numeric(vc$n_over_nstar)
  if (length(raw_ratio) != 1L) stop("vc_mmkh() returned an invalid N/N* ratio")
  is_capped <- is.finite(raw_ratio) && raw_ratio < 1
  corrected_ratio <- if (is_capped) 1 else raw_ratio
  plain_p <- suppressWarnings(as.numeric(vc$p_uncorrected))
  if (length(plain_p) != 1L || !is.finite(plain_p))
    plain_p <- tryCatch(as.numeric(Kendall::MannKendall(v)$sl), error = function(e) NA_real_)
  if (length(plain_p) != 1L || !is.finite(plain_p)) plain_p <- NA_real_
  if (is_capped && !is.finite(plain_p))
    stop("Cannot cap VC variance: uncorrected MK p-value is unavailable")
  vc_p <- if (is_capped) plain_p else as.numeric(vc$p)
  if (length(vc_p) != 1L) stop("vc_mmkh() returned an invalid VC p-value")
  status_to_report <- if (is_capped) "capped_uncorrected" else as.character(vc$status)
  out[, `:=`(tau = vc$tau, p_vc = vc_p, sen = vc$slope,
             n_over_neff = corrected_ratio, n_over_neff_raw = raw_ratio,
             p_mk_plain = plain_p, vc_status = status_to_report, rho1 = .rho1_one(v))]
  out
}
# Columns returned by .mk_series(); used to reuse basin results in the regional tables.
MK_SERIES_COLS <- c("n", "qa", "tau", "p_vc", "sen", "n_over_neff", "n_over_neff_raw",
                    "p_mk_plain", "rho1", "vc_status")

.tfpw_na <- function() data.table::data.table(tau_tfpw = NA_real_, p_tfpw = NA_real_, sen_tfpw = NA_real_)
.tfpw_series <- function(v) {
  r <- tryCatch(data.table::as.data.table(vectorized_mann_kendall_tfpw(matrix(v, nrow = 1L))), error = function(e) NULL)
  if (is.null(r) || !all(c("tau_tfpw", "p_value_tfpw", "sl_tfpw") %in% names(r))) return(.tfpw_na())
  data.table::data.table(tau_tfpw = r$tau_tfpw[1], p_tfpw = r$p_value_tfpw[1], sen_tfpw = r$sl_tfpw[1])
}

# Moving-block-bootstrap MK: resampling blocks keeps the short-range dependence and destroys any trend.
# MK S = concordant minus discordant ordered pairs, via the shared Fenwick-tree
# implementation in utils (.analysis_kendall_S): O(n log n), ties contribute zero.
.mk_S <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  as.numeric(.analysis_kendall_S(x))
}

# Detrend before resampling, so bootstrap samples approximate the stationary
# no-trend null while preserving within-block serial dependence. Observed S
# remains that of the ORIGINAL series. This is a sensitivity p-value, not
# the primary VC-HR98 / BH-FDR inference.
.mbb_mk_p <- function(x, L, B, seed) {
  n <- length(x)
  S0 <- .mk_S(x)
  if (!is.finite(S0) || n < 3L * L) return(c(S = S0, p = NA_real_))
  x_null <- as.numeric(stats::residuals(stats::lm(x ~ seq_along(x))))
  starts <- seq_len(n - L + 1L)
  nb <- ceiling(n / L)
  set.seed(seed)
  Sb <- vapply(seq_len(B), function(b) {
    s <- sample(starts, nb, replace = TRUE)
    ii <- as.vector(outer(0:(L - 1L), s, "+"))[seq_len(n)]
    .mk_S(x_null[ii])
  }, numeric(1))
  c(S = S0, p = (1 + sum(abs(Sb) >= abs(S0), na.rm = TRUE)) /
      (1 + sum(is.finite(Sb))))
}

.basin_full_row <- function(x, idx, sc, series, unit, L) {
  r <- .mk_series(x, if (series == "full_monthly") MIN_N_MONTHLY else MIN_N_YEARS)
  ok <- r$qa == "ok"; v <- .trim_finite(x)
  do_mbb <- ok && RUN_BLOCK_BOOTSTRAP
  mb <- if (do_mbb) .mbb_mk_p(v, L, MBB_B, seed = 7919L * as.integer(sc) + sum(utf8ToInt(idx))) else c(S = NA_real_, p = NA_real_)
  out <- cbind(data.table::data.table(index_type = idx, scale = as.integer(sc), series = series, slope_unit = unit), r,
               if (ok) .tfpw_series(v) else .tfpw_na(),
               data.table::data.table(mbb_S = unname(mb["S"]), mbb_p = unname(mb["p"]),
                                      mbb_block = if (do_mbb) as.integer(L) else NA_integer_,
                                      mbb_B = if (do_mbb) MBB_B else NA_integer_))
  out[, `:=`(vc_significant = is.finite(p_vc) & p_vc < 0.05,
             mbb_significant = is.finite(mbb_p) & mbb_p < 0.05)]
  # Soil-moisture basin trends require BOTH VC and block-bootstrap support.
  # Other indices retain VC as primary; bootstrap remains a robustness check.
  is_soil <- idx %in% c("ssmi", "ssi_sm") && identical(series, "full_monthly")
  out[, `:=`(vc_and_mbb_confirmed = vc_significant & mbb_significant,
             robust_significant = if (is_soil) vc_significant & mbb_significant else vc_significant,
             inference_note = if (is_soil && !is.finite(mbb_p)) "soil: bootstrap unavailable; not confirmed" else
               if (is_soil && vc_significant && !mbb_significant) "soil: VC-only; NOT robustly significant" else
                 if (is_soil) "soil: significance requires VC and bootstrap" else
                   "VC-primary significance; bootstrap is a separate sensitivity (robust field is VC-only here)")]
  out
}

ACC <- new.env()
for (s in c("basin_full", "basin_calmonth", "basin_monthly", "basin_annual",
            "snow_basin_season", "snow_basin_annual", "regional_full",
            "regional_calmonth", "homog_region", "homog_subperiod", "homog_series", "homog_pixel",
            "points", "point_ts", "qa")) ACC[[s]] <- list()
.acc <- function(slot, x) ACC[[slot]][[length(ACC[[slot]]) + 1L]] <- x
.record_qa <- function(idx, sc, analysis, reason) {
  tb <- table(reason)
  .acc("qa", data.table::data.table(index_type = idx, scale = as.integer(sc), analysis = analysis,
                                    qa = names(tb), n = as.integer(tb)))
}

# Uniform annual drought fractions: the denominator is fixed for every pixel
# and every year within each index. Annual windows missing even one required
# month remain NA, not a fraction over a changing subset of months. This avoids
# snow-season/full-year mixtures and partial accumulation years.
.annual_drought_fraction <- function(ts_mat, dates, index_type) {
  months <- as.integer(format(as.Date(dates), "%m"))
  years <- as.integer(format(as.Date(dates), "%Y"))
  requested <- DFREQ_SNOW_MONTHS[[index_type]]
  is_snow <- !is.null(requested)
  if (!is_snow) requested <- seq_len(DFREQ_MIN_MONTHS)
  # For SSPI, Dec of year y belongs to winter ending in year y+1.
  if (identical(index_type, "sspi")) years[months == 12L] <- years[months == 12L] + 1L
  in_season <- months %in% requested
  annual_years <- sort(unique(years[in_season]))
  n_expected <- length(requested)
  nr <- nrow(ts_mat)
  fraction <- matrix(NA_real_, nrow = nr, ncol = length(annual_years),
                     dimnames = list(NULL, as.character(annual_years)))
  for (j in seq_along(annual_years)) {
    at <- which(in_season & years == annual_years[j])
    # Guard missing/duplicated calendar months, including truncated winters.
    if (length(at) != n_expected || !setequal(months[at], requested) ||
        anyDuplicated(months[at])) next
    values <- ts_mat[, at, drop = FALSE]
    good <- rowSums(is.finite(values)) == n_expected
    if (any(good)) fraction[good, j] <-
      rowSums(values[good, , drop = FALSE] < DROUGHT_ONSET) / n_expected
  }
  definition <- if (identical(index_type, "swei")) "Jan-Apr; 4/4 months; calendar year" else
    if (identical(index_type, "sspi")) "Dec-Mar; 4/4 months; winter end year" else
      "Jan-Dec; 12/12 months; calendar year"
  list(fraction = fraction, years = annual_years, n_expected = n_expected,
       definition = definition)
}

# Non-overlapping, complete-year snow-window means (one value per winter or
# melt season). For SSPI, December is assigned to the winter ending next year.
# We deliberately do NOT concatenate separated months and call them adjacent
# observations in the monthly serial-correlation model.
.snow_basin_annual_means <- function(values, dates, index_type, month_window) {
  stopifnot(length(values) == length(dates), !anyDuplicated(dates))
  mon <- as.integer(format(dates, "%m"))
  yr <- as.integer(format(dates, "%Y"))
  if (identical(index_type, "sspi") && 12L %in% month_window)
    yr[mon == 12L] <- yr[mon == 12L] + 1L
  dat <- data.table::data.table(year = yr, month = mon, value = as.numeric(values))
  dat <- dat[month %in% month_window]
  if (!nrow(dat)) return(data.table::data.table(year = integer(), value = numeric()))
  ans <- dat[, .(value = if (.N == length(month_window) &&
                             !anyDuplicated(month) && setequal(month, month_window) &&
                             all(is.finite(value))) mean(value) else NA_real_), by = year]
  data.table::setorder(ans, year)
  ans
}

.snow_basin_season_test <- function(annual, idx, sc, window, label) {
  r <- .mk_series(annual$value, MIN_N_YEARS)
  ok <- identical(r$qa[[1L]], "ok")
  v <- .trim_finite(annual$value)
  mb <- if (ok && RUN_BLOCK_BOOTSTRAP)
    .mbb_mk_p(v, SNOW_SEASON_MBB_YEARS, MBB_B,
              seed = 7919L * as.integer(sc) + sum(utf8ToInt(idx)) +
                if (identical(label, "accumulation")) 177L else 319L)
  else c(S = NA_real_, p = NA_real_)
  cbind(data.table::data.table(index_type = idx, scale = as.integer(sc),
                               season = label, months = paste(month.abb[window], collapse = ","),
                               season_year_rule = if (identical(idx, "sspi") && 12L %in% window)
                                 "December belongs to next year" else "calendar year",
                               n_complete_years = sum(is.finite(annual$value)),
                               first_complete_year = if (any(is.finite(annual$value)))
                                 min(annual$year[is.finite(annual$value)]) else NA_integer_,
                               last_complete_year = if (any(is.finite(annual$value)))
                                 max(annual$year[is.finite(annual$value)]) else NA_integer_),
        r, data.table::data.table(mbb_p = unname(mb["p"]),
                                  mbb_block_years = SNOW_SEASON_MBB_YEARS,
                                  vc_significant = is.finite(r$p_vc) && r$p_vc < 0.05,
                                  mbb_significant = is.finite(mb["p"]) && mb["p"] < 0.05,
                                  # This is a secondary diagnostic; do not retrospectively select
                                  # a season or label it a confirmatory snow-only trend.
                                  inference_scope = "exploratory complete-season basin annual mean; VC primary, MBB sensitivity"))
}

# ==============================================================================
#   4. DROUGHT EVENTS, RUNS TEST, SPECTRAL PEAKS (unchanged statistics)
# ==============================================================================
# Per-pixel event statistics using the shared framework event definition.
vectorized_event_stats <- function(ts_matrix, min_duration,
                                   onset_thr = DROUGHT_ONSET, end_thr = DROUGHT_END) {
  n_pix <- nrow(ts_matrix)
  ne <- d12 <- d36 <- d712 <- d13 <- d4 <- integer(n_pix)
  md <- mx <- mI <- mS <- rep(NA_real_, n_pix)
  for (i in seq_len(n_pix)) {
    ev <- drought_event_table(ts_matrix[i, ], min_duration = min_duration,
                              onset_threshold = onset_thr, termination_threshold = end_thr)
    if (!nrow(ev)) next
    durs <- ev$duration_months; I_ev <- ev$mean_int; S_ev <- ev$severity
    ne[i] <- length(durs); md[i] <- mean(durs); mx[i] <- max(ev$max_int, na.rm = TRUE)
    d12[i] <- sum(durs <= 2L); d36[i] <- sum(durs >= 3L & durs <= 6L); d712[i] <- sum(durs >= 7L & durs <= 12L)
    d13[i] <- sum(durs >= 13L); d4[i] <- sum(durs >= 4L)
    mI[i] <- mean(I_ev, na.rm = TRUE); mS[i] <- mean(S_ev, na.rm = TRUE)
  }
  data.table::data.table(n_events = ne, mean_duration = md, max_intensity = mx, n_events_D12 = d12,
                         n_events_D36 = d36, n_events_D712 = d712, n_events_D13p = d13, n_D4p = d4,
                         mean_I = mI, mean_S = mS)
}

vectorized_runs_test <- function(ts_matrix, min_duration,
                                 onset_thr = DROUGHT_ONSET, end_thr = DROUGHT_END) {
  n_pix <- nrow(ts_matrix)
  pv_all <- rep(NA_real_, n_pix); clus_all <- rep(NA_character_, n_pix); filt <- rep(TRUE, n_pix)
  for (i in seq_len(n_pix)) {
    x <- .trim_finite(as.numeric(ts_matrix[i, ]))
    if (anyNA(x) || length(x) < 10L) next      # never join non-adjacent months
    # Binary event-membership series from the shared framework event definition.
    b <- integer(length(x))
    ev <- drought_event_runs(x, min_duration, onset_thr, end_thr)
    for (k in seq_len(nrow(ev))) b[ev$start_index[k]:ev$end_index[k]] <- 1L
    n <- length(b); n1 <- sum(b); n0 <- n - n1
    if (n1 < 5L || n0 < 5L) next
    nr <- length(rle(b)$lengths); er <- 2 * n0 * n1 / n + 1
    vr <- (2 * n0 * n1 * (2 * n0 * n1 - n)) / (n^2 * (n - 1))
    if (vr <= 0) next
    z <- (nr - er) / sqrt(vr)
    pv_all[i] <- 2 * stats::pnorm(-abs(z))
    clus_all[i] <- if (nr < er) "clustered" else if (nr > er) "dispersed" else "random"
    filt[i] <- FALSE
  }
  res <- data.table::data.table(p_value_runs = pv_all, clustering = clus_all, filtered_runs = filt)
  res[, p_value_runs_fdr := .bh(p_value_runs)]
  ok <- is.finite(res$p_value_runs)
  res$clustering[ok & (!is.finite(res$p_value_runs_fdr) | res$p_value_runs_fdr >= 0.05)] <- "not significant"
  res$clustering[!ok] <- NA_character_
  res
}

n_sim_spectral <- 100L
# Red-noise Monte Carlo spectral detection is implemented by utils J7d.
# For k >= 3, the MA(k-1) accumulation overlap violates both the runs-test
# independence null and an AR(1) red-noise spectral null. These are descriptive
# diagnostics ONLY, not evidence of clustering or spectral periodicity.
# The same caution applies to the spatial runs test above.
# ==============================================================================
#   5. BASIN DROUGHT EXTENT AND EVENT CATALOGUE (unchanged)
# ==============================================================================
compute_basin_extent_and_event_catalog <- function(ts_matrix, n_years, index_label, out_dir, lon_vec, lat_vec,
                                                   onset_thr = DROUGHT_ONSET, end_thr = DROUGHT_END,
                                                   start_year = 1950L, win_years = 30L,
                                                   area_weights = NULL, min_duration) {
  stopifnot(length(lon_vec) == nrow(ts_matrix), length(lat_vec) == nrow(ts_matrix))
  use_w <- !is.null(area_weights) && length(area_weights) == nrow(ts_matrix) &&
    all(is.finite(area_weights)) && all(area_weights > 0)
  if (!is.null(area_weights) && !use_w) warning("area_weights invalid; equal weighting used.")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ms <- seq_len(n_years * 12L)
  yr_seq <- start_year + floor((ms - 1L) / 12L); mo_seq <- ((ms - 1L) %% 12L) + 1L
  n_col <- min(length(ms), ncol(ts_matrix))
  tm <- ts_matrix[, seq_len(n_col), drop = FALSE]; ok <- !is.na(tm)
  w <- if (use_w) as.numeric(area_weights) else rep(1, nrow(tm))
  dr <- ok & (tm < onset_thr)
  pct_area <- ifelse(colSums(ok) > 0, 100 * colSums(dr * w) / colSums(ok * w), NA_real_)
  extent_ts <- data.frame(date = as.Date(sprintf("%04d-%02d-01", yr_seq, mo_seq))[seq_len(n_col)],
                          year = yr_seq[seq_len(n_col)], month = mo_seq[seq_len(n_col)], pct_area = pct_area)
  utils::write.csv(extent_ts, file.path(out_dir, sprintf("%s_basin_extent_monthly.csv", index_label)), row.names = FALSE)
  extent_annual <- extent_ts |> dplyr::group_by(year) |>
    dplyr::summarise(n_valid_months = sum(is.finite(pct_area)),
                     mean_pct = if (sum(is.finite(pct_area)) == 12L)
                       mean(pct_area) else NA_real_, .groups = "drop")
  # Only 12/12 months count. Do not splice gaps between the first and last
  # eligible year; partial accumulation years must not enter the MK test.
  utils::write.csv(extent_annual, file.path(out_dir, sprintf("%s_basin_extent_annual.csv", index_label)),
                   row.names = FALSE)
  valid_years <- which(is.finite(extent_annual$mean_pct))
  xc <- if (length(valid_years))
    extent_annual$mean_pct[min(valid_years):max(valid_years)] else numeric()
  mk_ext <- list(tau = NA_real_, p_value = NA_real_, sens_slope = NA_real_,
                 n_complete_years = length(valid_years),
                 n_incomplete_years = sum(extent_annual$n_valid_months < 12L),
                 qa = if (!length(valid_years)) "no_complete_years" else if (anyNA(xc))
                   "internal_year_gaps" else "ok")
  if (identical(mk_ext$qa, "ok") && length(xc) >= MIN_N_YEARS &&
      length(unique(xc)) >= 5L) {
    hr <- tryCatch(vc_mmkh(xc), error = function(e) NULL)
    if (!is.null(hr)) mk_ext[c("tau", "p_value", "sens_slope")] <-
        list(as.numeric(hr$tau), as.numeric(hr$p), as.numeric(hr$slope))
  }
  mk_ext_df <- as.data.frame(mk_ext)
  utils::write.csv(mk_ext_df, file.path(out_dir, sprintf("%s_basin_extent_MK.csv", index_label)), row.names = FALSE)
  if (nrow(extent_annual) >= win_years)
    utils::write.csv(data.frame(year = extent_annual$year,
                                rolling_mean = zoo::rollmean(extent_annual$mean_pct, k = win_years, fill = NA, align = "right")),
                     file.path(out_dir, sprintf("%s_basin_extent_%dyr_rolling.csv", index_label, win_years)), row.names = FALSE)
  p_ext <- ggplot2::ggplot(extent_annual, ggplot2::aes(year, mean_pct)) + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(title = sprintf("%s — Annual Mean Basin-Wide Drought Extent", index_label),
                  subtitle = sprintf("Hamed–Rao MK: tau = %.3f, p = %.4g", mk_ext$tau, mk_ext$p_value),
                  x = "Year", y = "% Basin Area in Drought") +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file.path(out_dir, sprintf("%s_basin_extent_timeseries.png", index_label)), p_ext, width = 12, height = 5, dpi = 150)
  ev_list <- lapply(seq_len(nrow(ts_matrix)), function(i) {
    ev <- drought_event_table(ts_matrix[i, ], min_duration = min_duration,
                              onset_threshold = onset_thr, termination_threshold = end_thr)
    if (!nrow(ev)) return(NULL)
    data.frame(pixel = i, lon = lon_vec[i], lat = lat_vec[i], ev)
  })
  events_df <- do.call(rbind, ev_list)
  if (!is.null(events_df) && nrow(events_df)) {
    events_df$duration_class <- cut(events_df$duration_months, breaks = c(0, 2, 6, 12, Inf),
                                    labels = c("D1-2", "D3-6", "D7-12", "D13+"), right = TRUE)
    events_df$index <- index_label; events_df$minimum_duration_months <- as.integer(min_duration)
    utils::write.csv(events_df, file.path(out_dir, sprintf("%s_basin_event_catalog.csv", index_label)), row.names = FALSE)
  }
  invisible(list(extent_ts = extent_ts, mk_extent = mk_ext_df, events_df = events_df))
}


# ==============================================================================
#   6. HOMOGENEITY BATTERY (Script 11 equivalent, annual-resolution series only)
# ==============================================================================
.season_matrix <- function(mat, dates, months) {
  yr <- as.integer(format(dates, "%Y")); mo <- as.integer(format(dates, "%m"))
  wrap <- any(months >= 10L) && any(months <= 3L)
  sy <- yr + as.integer(wrap & mo >= 10L)
  keep <- mo %in% months
  yrs <- sort(unique(sy[keep]))
  out <- matrix(NA_real_, nrow(mat), length(yrs))
  for (j in seq_along(yrs)) {
    cols <- which(keep & sy == yrs[j])
    if (length(cols) == length(months)) out[, j] <- rowMeans(mat[, cols, drop = FALSE])   # NA if any month missing
  }
  list(mat = out, years = yrs)
}


# The shared utils routine is the only homogeneity battery used by Scripts 11 and 13.
.homog_battery <- function(v, yrs, sy)
  analysis_homogeneity_battery(v, yrs, step_year = sy,
                               m_sim = HOMOG_MC, min_unique = 5L)

.subperiod_trends <- function(v, yrs, sy) {
  win <- list(full = range(yrs), before_step = c(min(yrs), sy - 1L), from_step = c(sy, max(yrs)),
              wmo_1991_2020 = c(1991L, 2020L))
  data.table::rbindlist(lapply(names(win), function(w) {
    k <- yrs >= win[[w]][1] & yrs <= win[[w]][2]
    r <- .mk_series(v[k], ANALYSIS_MIN_TREND_N)
    cbind(data.table::data.table(window = w, year_from = win[[w]][1], year_to = win[[w]][2]),
          r[, .(n, qa, tau, p_vc, sen)])
  }))
}

.pixel_step <- function(v, yrs, sy) {
  na <- data.table::data.table(n_years = sum(is.finite(v)), pettitt_p = NA_real_, pettitt_first_year_after = NA_real_,
                               mean_before = NA_real_, diff_at_step = NA_real_, welch_p = NA_real_)
  ok <- which(is.finite(v)); if (length(ok) < MIN_N_YEARS) return(na)
  i <- min(ok):max(ok); if (anyNA(v[i])) return(na)
  v <- v[i]; y <- yrs[i]; if (length(unique(v)) < 3L) return(na)
  pc <- analysis_pettitt(v); b <- v[y < sy]; a <- v[y >= sy]
  data.table::data.table(n_years = length(v), pettitt_p = pc$p, pettitt_first_year_after = y[pc$tau] + 1,
                         mean_before = if (length(b)) mean(b) else NA_real_,
                         diff_at_step = if (length(a) && length(b)) mean(a) - mean(b) else NA_real_,
                         welch_p = if (length(a) >= 5L && length(b) >= 5L)
                           tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_) else NA_real_)
}

# Scan independent pixels in contiguous blocks. The order of the output rows is
# identical to lapply(seq_len(nrow(mat)), ...). The only expensive work on each
# worker is a deterministic Pettitt/Welch test; no Monte Carlo is parallelized.
.pixel_scan_parallel <- function(mat, yrs, step_year) {
  n <- nrow(mat)
  if (n == 0L) return(list())
  if (N_CORES <= 1L) return(lapply(seq_len(n), function(i)
    .pixel_step(mat[i, ], yrs, step_year)))
  block_size <- max(1L, ceiling(n / (N_CORES * 4L)))
  ii <- seq_len(n)
  blocks <- split(ii, ceiling(ii / block_size))
  work <- function(ids) lapply(ids, function(i) .pixel_step(mat[i, ], yrs, step_year))
  if (is_windows) {
    if (is.null(cl)) stop("Parallel pixel scan requires an initialized PSOCK cluster")
    parallel::clusterExport(cl, c(".pixel_step", "MIN_N_YEARS"), envir = .GlobalEnv)
    out <- parallel::parLapply(cl, blocks, function(ids, mat, yrs, step_year) {
      lapply(ids, function(i) .pixel_step(mat[i, ], yrs, step_year))
    }, mat = mat, yrs = yrs, step_year = step_year)
  } else {
    out <- parallel::mclapply(blocks, work,
                              mc.cores = min(N_CORES, length(blocks)),
                              mc.preschedule = TRUE)
  }
  unlist(out, recursive = FALSE)
}

# time_info: list(annual = FALSE, dates = <monthly dates>) or list(annual = TRUE, years = <years>)
run_homog_spec <- function(spec, pix_mat, lon, lat, reg_vals, reg_names, time_info) {
  sid <- paste0(INDEX_LABEL[[spec$index_type]], "-", spec$scale, "_", spec$suffix)
  months <- spec$months[[1]]
  season <- function(m) if (isTRUE(time_info$annual)) list(mat = m, years = time_info$years)
  else .season_matrix(m, time_info$dates, months)
  meta <- data.table::data.table(series = sid, index_type = spec$index_type, scale = spec$scale, season = spec$desc)
  for (rid in names(reg_vals)) {
    s <- season(matrix(reg_vals[[rid]], nrow = 1L)); v <- as.numeric(s$mat[1, ]); y <- s$years
    .acc("homog_region", cbind(data.table::data.table(region_id = rid, region_name = reg_names[[rid]]), meta,
                               .homog_battery(v, y, STEP_YEAR)))
    .acc("homog_subperiod", cbind(data.table::data.table(region_id = rid), meta, .subperiod_trends(v, y, STEP_YEAR)))
    .acc("homog_series", data.table::data.table(region_id = rid, series = sid, year = y, value = v))
  }
  s <- season(pix_mat); y <- s$years
  px <- data.table::rbindlist(.pixel_scan_parallel(s$mat, y, STEP_YEAR))
  px <- cbind(data.table::data.table(lon = lon, lat = lat), px)
  px[, `:=`(pettitt_p_fdr = .bh(pettitt_p), welch_p_fdr = .bh(welch_p))]
  data.table::fwrite(px, file.path(HOMOG_DIR, "pixel", sprintf("pixel_step_tests_%s.csv", sid)))
  tested <- is.finite(px$pettitt_p); sig <- tested & is.finite(px$pettitt_p_fdr) & px$pettitt_p_fdr < 0.05
  .acc("homog_pixel", cbind(meta, data.table::data.table(
    n_pixels = sum(tested), pct_pettitt_fdr05 = if (any(tested)) 100 * mean(sig[tested]) else NA_real_,
    median_break_year_fdr_sig = if (any(sig)) stats::median(px$pettitt_first_year_after[sig]) else NA_real_,
    pct_welch_fdr05 = if (any(is.finite(px$welch_p))) 100 * mean((px$welch_p_fdr < 0.05)[is.finite(px$welch_p)]) else NA_real_,
    median_diff_at_step = stats::median(px$diff_at_step, na.rm = TRUE))))
  cat(sprintf("    homogeneity %-18s pixels tested %4d | Pettitt FDR-sig %5.1f%%\n", sid, sum(tested),
              if (any(tested)) 100 * mean(sig[tested]) else NA_real_))
}


# ==============================================================================
#   7. TIME-SERIES MATRIX BUILDER (unchanged contract)
# ==============================================================================
build_index_matrix <- function(files, stacks, n_pix) {
  n_files <- length(stacks); n_l <- terra::nlyr(stacks[[1]])
  if (n_files == 1L && n_l > 12L) {
    dates <- extract_dates_from_nc(files[1], n_l)
    m <- as.matrix(terra::values(stacks[[1]], mat = TRUE, na.rm = FALSE))
    if (!identical(dim(m), c(as.integer(n_pix), as.integer(n_l)))) stop("Unexpected single-file matrix dimensions.")
    years <- sort(unique(as.integer(format(dates, "%Y"))))
  } else {
    if (n_files != 12L) stop("Expected 12 calendar-month NetCDF files; found ", n_files, ".")
    mo <- suppressWarnings(as.integer(sub(".*month([0-9]{2}).*", "\\1", basename(files), ignore.case = TRUE)))
    if (anyNA(mo) || !identical(sort(mo), 1:12)) stop("Seasonal files must contain one file per month 01-12: ",
                                                      paste(basename(files), collapse = ", "))
    o <- order(mo); files <- files[o]; stacks <- stacks[o]
    d1 <- extract_dates_from_nc(files[1], terra::nlyr(stacks[[1]]))
    years <- as.integer(format(d1, "%Y")); ny <- length(years)
    m <- matrix(NA_real_, n_pix, ny * 12L)
    for (k in 1:12) {
      v <- as.matrix(terra::values(stacks[[k]], mat = TRUE, na.rm = FALSE))
      dk <- extract_dates_from_nc(files[k], terra::nlyr(stacks[[k]]))
      if (length(dk) != ny || anyNA(dk) || any(as.integer(format(dk, "%m")) != k) ||
          !identical(as.integer(format(dk, "%Y")), years))
        stop("Date/year alignment failure in ", basename(files[k]))
      if (!identical(dim(v), c(as.integer(n_pix), as.integer(ny)))) stop("Unexpected dimensions in ", basename(files[k]))
      m[, (seq_len(ny) - 1L) * 12L + k] <- v
    }
    dates <- as.Date(sprintf("%04d-%02d-01", rep(years, each = 12L), rep(1:12, times = ny)))
  }
  if (ncol(m) != length(dates)) stop("Matrix/date mismatch.")
  validate_monthly_record(dates, label = "6trend time-series matrix calendar", require_complete = TRUE)
  list(ts_mat = m, dates = dates, years = years)
}

.point_rows <- function(r_tmpl, cell_vec) {
  if (is.null(SPECIFIC_PTS_DT) || !nrow(SPECIFIC_PTS_DT)) return(NULL)
  v <- terra::vect(as.data.frame(SPECIFIC_PTS_DT), geom = c("lon", "lat"), crs = "EPSG:4326")
  if (!terra::is.lonlat(r_tmpl)) v <- terra::project(v, terra::crs(r_tmpl))
  cells <- terra::cellFromXY(r_tmpl, terra::crds(v))
  out <- data.table::data.table(point_id = SPECIFIC_PTS_DT$point_id, lon_wgs84 = SPECIFIC_PTS_DT$lon,
                                lat_wgs84 = SPECIFIC_PTS_DT$lat, cell = cells, row = match(cells, cell_vec))
  if (any(is.na(out$row))) cat("    ⚠ point(s) outside the valid basin pixels:", paste(out$point_id[is.na(out$row)], collapse = ", "), "\n")
  out[!is.na(row)]
}


# ==============================================================================
#   8. MAIN PROCESSING FUNCTION — MONTHLY INDICES
# ==============================================================================
process_index <- function(index_type, sc, find_fn, seas_dir) {
  idx <- tolower(index_type); sc <- as.integer(sc)
  index_label <- sprintf("%s-%02d", INDEX_LABEL[[idx]], sc)
  min_dur <- get_event_min_duration(sc)
  # Inclusive wall-clock timings locate the expensive stages without changing
  # any statistics, random-number draws, or result-file schemas.
  .time_stage <- local({
    t_all <- proc.time()[["elapsed"]]
    t_prev <- t_all
    last <- NULL
    function(next_stage = NULL) {
      now <- proc.time()[["elapsed"]]
      if (!is.null(last))
        cat(sprintf("  [TIME] %s: %.1f s\n", last, now - t_prev))
      if (is.null(next_stage))
        cat(sprintf("  [TIME] %s total: %.1f s\n", index_label, now - t_all))
      last <<- next_stage
      t_prev <<- now
      invisible(NULL)
    }
  })
  on.exit(.time_stage(NULL), add = TRUE)
  cat(sprintf("\n╔══ %s ══╗  (event minimum duration %d months)\n", index_label, min_dur))
  
  # ---- [1/9] files and native grid ------------------------------------------------------------
  .time_stage("[1/9] Files/native grid")
  cat("[1/9] Files and native grid...\n")
  files <- find_fn(seas_dir, index_type, sc)
  if (!length(files)) stop(sprintf("No monthly NetCDF files for %s in %s", index_label, seas_dir))
  stacks <- lapply(files, function(f) tryCatch(terra::rast(f), error = function(e) {
    cat("  ⚠ Failed: ", basename(f), " | ", conditionMessage(e), "\n", sep = ""); NULL }))
  bad <- vapply(stacks, is.null, logical(1))
  if (any(bad)) stop(sprintf("Failed to load %d NetCDF(s) for %s: %s", sum(bad), index_label, paste(basename(files[bad]), collapse = ", ")))
  r_tmpl <- stacks[[1]][[1]]
  if (!(idx %in% c("ssmi", "ssi_sm"))) {
    # Script 15 rasterizes non-soil CSV cells on the ERA5-Land precipitation grid.
    # Fail here rather than silently mislocating spatial statistics downstream.
    assert_same_native_grid(r_tmpl, PRECIP_TEMPLATE_13, paste0(index_label, " vs precipitation grid"))
  }
  if (is.na(terra::crs(r_tmpl)) || !nzchar(terra::crs(r_tmpl))) stop("Native index grid has no CRS: ", basename(files[1]))
  for (k in seq_along(stacks)) assert_same_native_grid(stacks[[k]], r_tmpl, sprintf("%s file %d (%s)", index_label, k, basename(files[k])))
  cat(sprintf("  ✓ %d files | %d x %d cells | res %.6f x %.6f\n", length(files), terra::nrow(r_tmpl),
              terra::ncol(r_tmpl), terra::res(r_tmpl)[1], terra::res(r_tmpl)[2]))
  xy_all <- terra::xyFromCell(r_tmpl, seq_len(terra::ncell(r_tmpl)))
  
  # ---- [2/9] time-series matrix ---------------------------------------------------------------
  .time_stage("[2/9] Time-series matrix")
  cat("[2/9] Time-series matrix...\n")
  tsm <- build_index_matrix(files, stacks, terra::ncell(r_tmpl))
  ts_mat <- tsm$ts_mat; dates_ts <- tsm$dates; years <- tsm$years; n_years <- length(years)
  rm(stacks, tsm); invisible(gc(verbose = FALSE))
  ts_mat[!is.finite(ts_mat)] <- NA_real_
  mo_ts <- as.integer(format(dates_ts, "%m")); yr_ts <- as.integer(format(dates_ts, "%Y"))
  cat(sprintf("  ✓ %d months (%s to %s) | finite %.1f%%\n", ncol(ts_mat), min(dates_ts), max(dates_ts), 100 * mean(is.finite(ts_mat))))
  
  bg <- get_basin_grid_cached(r_tmpl)
  valid_pix <- (rowSums(is.finite(ts_mat)) > 0) & bg$cover_all > 0
  if (!any(valid_pix)) stop("No finite-data pixels inside the basin for ", index_label)
  eff_area_all <- bg$area_all * bg$cover_all
  area_valid <- eff_area_all[valid_pix]
  if (any(!is.finite(area_valid) | area_valid <= 0)) stop("Invalid effective area for in-basin pixels: ", index_label)
  ts_mat <- ts_mat[valid_pix, , drop = FALSE]
  cell_vec <- which(valid_pix); lon_vec <- xy_all[valid_pix, 1]; lat_vec <- xy_all[valid_pix, 2]; n_pix <- length(cell_vec)
  cat(sprintf("  ✓ Basin pixels: %d (partial coverage %d)\n", n_pix, sum(bg$cover_all[valid_pix] < 1)))
  
  # ---- [3/9] regional and basin series + their trends -----------------------------------------
  .time_stage("[3/9] Basin/regional trends")
  cat("[3/9] Regional and basin series and trends...\n")
  variant <- if (idx %in% c("ssmi", "ssi_sm")) SOIL_LAYER else NULL
  basin_csv_dir <- BASIN_SERIES_DIRS[[idx]]
  if (is.null(basin_csv_dir)) stop("No basin-series folder configured for ", idx)
  stored <- read_basin_series(idx, sc, basin_csv_dir, times = dates_ts, variant = variant)
  if (VERIFY_BASIN_SERIES)
    verify_basin_series(stored, ts_mat, r_tmpl, cells = cell_vec, label = index_label)
  regional_index <- stored$regional
  basin_avg_vec  <- as.numeric(stored$basin$value)
  basin_stats    <- list(value = basin_avg_vec, coverage = as.numeric(stored$basin$coverage))
  .acc("basin_monthly", data.table::data.table(index_type = idx, scale = sc, date = dates_ts, value = basin_avg_vec,
                                               valid_area_coverage_fraction = basin_stats$coverage))
  
  reg_split <- split(regional_index[order(region_id, date)], by = "region_id")
  reg_vals  <- lapply(reg_split, function(d) d$value)
  reg_names <- vapply(reg_split, function(d) as.character(d$region_name[1]), character(1))
  if (!BASIN_RID %in% names(reg_vals)) { reg_vals[[BASIN_RID]] <- basin_avg_vec; reg_names[BASIN_RID] <- "Nechako Basin" }
  
  bfull <- .basin_full_row(basin_avg_vec, idx, sc, "full_monthly", "per month", max(12L, 2L * sc))
  .acc("basin_full", bfull)
  if (idx %in% names(SNOW_BASIN_WINDOWS)) {
    for (season in names(SNOW_BASIN_WINDOWS[[idx]])) {
      win <- SNOW_BASIN_WINDOWS[[idx]][[season]]
      ann <- .snow_basin_annual_means(basin_avg_vec, dates_ts, idx, win)
      .acc("snow_basin_annual", data.table::data.table(
        index_type = idx, scale = as.integer(sc), season = season,
        year = ann$year, value = ann$value))
      sr <- .snow_basin_season_test(ann, idx, sc, win, season)
      .acc("snow_basin_season", sr)
      cat(sprintf("  ✓ %s (%s) basin: %d complete years | VC p=%.3g | block-bootstrap p=%.3g | QA=%s (exploratory)\n",
                  season, paste(month.abb[win], collapse = ","),
                  sr$n_complete_years, sr$p_vc, sr$mbb_p, sr$qa))
    }
  }
  bcm <- data.table::rbindlist(lapply(1:12, function(m) {
    x <- basin_avg_vec[mo_ts == m]; r <- .mk_series(x, MIN_N_YEARS)
    cbind(data.table::data.table(index_type = idx, scale = sc, month = m, slope_unit = "per year"), r,
          if (r$qa == "ok") .tfpw_series(.trim_finite(x)) else .tfpw_na())
  }))
  .acc("basin_calmonth", bcm)
  
  # Regional tables. The whole-basin row reuses the basin tests above: it is the
  # same series (verified against the pixels), so it is not tested a second time.
  for (rid in names(reg_vals)) {
    rhead <- data.table::data.table(region_id = rid, region_name = reg_names[[rid]],
                                    index_type = idx, scale = sc)
    if (rid == BASIN_RID) {
      full_r <- bfull[, ..MK_SERIES_COLS]
      cm_r   <- bcm[, c("month", MK_SERIES_COLS), with = FALSE]
    } else {
      x <- reg_vals[[rid]]
      full_r <- .mk_series(x, MIN_N_MONTHLY)
      cm_r <- data.table::rbindlist(lapply(1:12, function(m)
        cbind(data.table::data.table(month = m), .mk_series(x[mo_ts == m], MIN_N_YEARS))))
    }
    .acc("regional_full", cbind(rhead, series = "full_monthly", slope_unit = "per month", full_r))
    .acc("regional_calmonth", cbind(rhead[rep(1L, nrow(cm_r))], slope_unit = "per year", cm_r))
  }
  cat(sprintf("  ✓ Basin full series: tau %+.3f | p VC %.3g | p TFPW %.3g [diagnostic only if k>=3] | p block-bootstrap %.3g | N/N* %.2f (raw %.2f) | primary-rule %s | VC+MBB %s\n",
              bfull$tau, bfull$p_vc, bfull$p_tfpw, bfull$mbb_p, bfull$n_over_neff,
              bfull$n_over_neff_raw, if (bfull$robust_significant) "YES" else "NO",
              if (bfull$vc_and_mbb_confirmed) "YES" else "NO"))
  
  # ---- [4/9] full-series pixel trends ---------------------------------------------------------
  .time_stage("[4/9] Full-series pixel MK")
  cat("[4/9] Full-series pixel trends (VC-HR98 primary; TFPW k>=3 REJECTED diagnostic)...\n")
  results <- data.table::data.table(cell = cell_vec, lon = lon_vec, lat = lat_vec, n_months = rowSums(is.finite(ts_mat)))
  results[, qa_full := .qa_rows(ts_mat, MIN_N_MONTHLY)]
  results[, `:=`(modal_tie_fraction = .modal_frac_rows(ts_mat), rho1_detrended = .rho1_rows(ts_mat))]
  results <- cbind(results, .mk_pixel(ts_mat, results$qa_full == "ok"))
  results[, `:=`(tfpw_method = "Yue2002_shared_J4_J6",
                 tfpw_interpretation = if (sc >= 3L) "REJECTED: AR1 prewhitening invalid for overlapping monthly k>=3" else "diagnostic: k=1") ]
  .record_qa(idx, sc, "full_monthly", results$qa_full)
  .record_qa(idx, sc, "vc_status_full_monthly",data.table::fifelse(is.na(results$vc_status), "not_tested", results$vc_status))
  cat(sprintf("  ✓ tested %d / %d pixels | mean detrended lag-1 rho %.2f\n",
              sum(results$qa_full == "ok"), n_pix, mean(results$rho1_detrended, na.rm = TRUE)))
  
  # ---- [5/9] events, PELT, runs, spectra, drought-month fraction ------------------------------
  .time_stage("[5/9] Events/PELT/runs/spectra")
  cat("[5/9] Events, PELT, runs test, spectra, drought-month fraction...\n")
  results <- cbind(results, vectorized_event_stats(ts_mat, min_dur, DROUGHT_ONSET, DROUGHT_END))
  cat(sprintf("  ✓ events (onset < %.2f, end >= %.2f, min %d months): mean %.1f per pixel (D13+ max %d)\n",
              DROUGHT_ONSET, DROUGHT_END, min_dur, mean(results$n_events), max(results$n_events_D13p)))
  results <- cbind(results, vectorized_regime_shift_pelt(ts_mat, years))
  results <- cbind(results, vectorized_runs_test(ts_mat, min_dur, DROUGHT_ONSET, DROUGHT_END))
  spectral_t0 <- proc.time()[["elapsed"]]
  results <- cbind(results, vectorized_spectral_peaks(ts_mat, n_sim_spectral = n_sim_spectral))
  cat(sprintf("  [TIME] Spectral Monte Carlo (B=%d): %.1f s\n",
              n_sim_spectral, proc.time()[["elapsed"]] - spectral_t0))
  
  if (RUN_DROUGHT_FREQ) {
    annual <- .annual_drought_fraction(ts_mat, dates_ts, idx)
    fr <- annual$fraction
    q <- .qa_rows(fr, MIN_N_YEARS)
    mk <- .mk_pixel(fr, q == "ok")
    avg <- rowMeans(fr, na.rm = TRUE)
    avg[!is.finite(avg)] <- NA_real_
    results[, `:=`(dfreq_mean = avg, dfreq_qa = q, dfreq_tau_vc = mk$tau_vc,
                   dfreq_p_vc = mk$p_value_vc, dfreq_sl_vc = mk$sl_vc,
                   dfreq_n_complete_years = rowSums(is.finite(fr)),
                   dfreq_months_per_year = annual$n_expected,
                   dfreq_definition = annual$definition)]
    results[, dfreq_p_fdr_vc := .bh(dfreq_p_vc)]
    .record_qa(idx, sc, "drought_month_fraction", q)
    cat(sprintf("  ✓ Annual fraction (%s): VC tested %d/%d pixels; ties %d; BH family n=%d%s\n",
                annual$definition, sum(is.finite(results$dfreq_p_vc)), n_pix,
                sum(q == "ties"), sum(is.finite(results$dfreq_p_vc)),
                if (sc > 12L) "; exploratory, not in Fig 9" else ""))
  }
  results[, `:=`(p_fdr_vc = .bh(p_value_vc), p_fdr_tfpw = .bh(p_value_tfpw))]
  # No per-pixel block-bootstrap test is performed. VC BH discoveries cannot
  # inherit basin bootstrap confirmation (different statistical units).
  results[, `:=`(pixel_inference_scope = "pixel VC-HR98 BH-FDR only; no pixel block-bootstrap",
                 basin_full_vc_p = bfull$p_vc,
                 basin_full_mbb_p = bfull$mbb_p,
                 basin_full_robust_significant = bfull$robust_significant,
                 basin_full_vc_and_mbb_confirmed = bfull$vc_and_mbb_confirmed)]
  if (idx %in% c("ssmi", "ssi_sm"))
    results[, pixel_inference_scope := "SOIL: VC-only pixel discoveries; basin VC+MBB NOT confirmed if robust=FALSE"]
  if (idx == "swei")
    results[, dfreq_inference_note := "Jan-Apr 4-month fractions; filter if tied; no trend estimate when qa=ties"]
  
  # ---- [6/9] calendar-month pixel trends ------------------------------------------------------
  cm_all <- NULL
  cm_file <- file.path(TREND_DIR, sprintf("%s_%02d_calmonth_results.csv", idx, sc))
  if (RUN_CALENDAR_MONTH) {
    .time_stage("[6/9] Calendar-month MK")
    cat("[6/9] Calendar-month pixel trends (VC-HR98 + TFPW, BH-FDR within month)...\n")
    cm_all <- data.table::rbindlist(lapply(1:12, function(m) {
      sub <- ts_mat[, mo_ts == m, drop = FALSE]; q <- .qa_rows(sub, MIN_N_YEARS)
      d <- cbind(data.table::data.table(cell = cell_vec, lon = lon_vec, lat = lat_vec, month = m,
                                        n_years = rowSums(is.finite(sub)), qa = q), .mk_pixel(sub, q == "ok"))
      d[, `:=`(p_fdr_vc = .bh(p_value_vc), p_fdr_tfpw = .bh(p_value_tfpw))]
      d
    }), use.names = TRUE, fill = TRUE)
    cm_all[, `:=`(index_type = idx, timescale = sc, slope_unit = "per year",
                  pixel_inference_scope = "calendar-month pixel VC-HR98 BH-FDR only; no bootstrap")]
    if (idx %in% c("ssmi", "ssi_sm"))
      cm_all[, pixel_inference_scope := "SOIL: calendar-month pixel VC-only; not confirmed by basin bootstrap"]
    if (idx %in% c("swei", "sspi"))
      cm_all[, snow_interpretation := ifelse(month %in% 5:9,
                                             "melt or low-snow season: standardized index may reflect snow presence",
                                             "fixed calendar month: raw SWE thresholds and snow presence not verified here")]
    .record_qa(idx, sc, "calendar_month", cm_all$qa)
    .record_qa(idx, sc, "vc_status_calendar_month",data.table::fifelse(is.na(cm_all$vc_status), "not_tested", cm_all$vc_status))
    data.table::fwrite(cm_all, cm_file)
    s <- cm_all[, .(tested = sum(qa == "ok"), fdr = sum(p_fdr_vc < 0.05, na.rm = TRUE)), by = month]
    cat("  ✓ VC FDR-significant pixels by month: ",
        paste(sprintf("%s %d/%d", month.abb[s$month], s$fdr, s$tested), collapse = " | "), "\n", sep = "")
  }
  
  # ---- [7/9] homogeneity on annual-resolution series ------------------------------------------
  if (RUN_HOMOGENEITY) {
    specs <- HOMOG_SPECS[index_type == idx & scale == sc]
    if (nrow(specs)) {
      .time_stage("[7/9] Homogeneity")
      cat("[7/9] Homogeneity battery...\n")
      for (q in seq_len(nrow(specs)))
        run_homog_spec(specs[q], ts_mat, lon_vec, lat_vec, reg_vals, reg_names, list(annual = FALSE, dates = dates_ts))
    }
  }
  
  # ---- [8/9] basin CSVs, extent, points -------------------------------------------------------
  .time_stage("[8/9] Extent/points")
  cat("[8/9] Basin CSVs, extent, specific points...\n")
  basin_avg_path <- file.path(basin_csv_dir, sprintf("%s_%02d_basin_averaged_by_month.csv", idx, sc))
  basin_cov_path <- file.path(basin_csv_dir, sprintf("%s_%02d_basin_valid_area_coverage_by_month.csv", idx, sc))
  
  extent_out <- file.path(TREND_DIR, "basin_extent")
  compute_basin_extent_and_event_catalog(ts_mat, n_years, index_label, extent_out, lon_vec, lat_vec,
                                         start_year = min(years), area_weights = area_valid, min_duration = min_dur)
  
  pr <- .point_rows(r_tmpl, cell_vec)
  if (!is.null(pr) && nrow(pr)) {
    keep_cols <- intersect(c("cell", "lon", "lat", "n_months", "qa_full", "rho1_detrended", "tau_vc", "p_value_vc",
                             "sl_vc", "tau_tfpw", "p_value_tfpw", "sl_tfpw", "p_fdr_vc", "p_fdr_tfpw",
                             "dfreq_qa", "dfreq_n_complete_years", "dfreq_months_per_year", "dfreq_definition",
                             "dfreq_tau_vc", "dfreq_p_vc", "dfreq_sl_vc", "dfreq_p_fdr_vc"), names(results))
    cm_cols <- c("month", "n_years", "qa", "tau_vc", "p_value_vc", "sl_vc", "tau_tfpw", "p_value_tfpw", "sl_tfpw",
                 "p_fdr_vc", "p_fdr_tfpw")
    for (j in seq_len(nrow(pr))) {
      rr <- pr$row[j]
      base <- data.table::data.table(point_id = pr$point_id[j], lon_wgs84 = pr$lon_wgs84[j], lat_wgs84 = pr$lat_wgs84[j],
                                     index_type = idx, scale = sc)
      .acc("points", cbind(base, period = "full_monthly", results[rr, ..keep_cols]))
      if (!is.null(cm_all)) {
        cc <- results$cell[rr]
        cmr <- cm_all[cell == cc, intersect(cm_cols, names(cm_all)), with = FALSE]
        .acc("points", cbind(base[rep(1L, nrow(cmr))], period = "calendar_month", cmr))
      }
      .acc("point_ts", data.table::data.table(point_id = pr$point_id[j], index_type = idx, scale = sc,
                                              date = dates_ts, value = ts_mat[rr, ]))
    }
  }
  
  # ---- [9/9] write pixel results --------------------------------------------------------------
  .time_stage("[9/9] Save pixel outputs")
  cat("[9/9] Writing pixel outputs...\n")
  out_file <- file.path(TREND_DIR, sprintf("%s_%02d_results.csv", idx, sc))
  data.table::fwrite(results, out_file, showProgress = FALSE)
  if (idx %in% c("ssmi", "ssi_sm")) {
    # Explicitly labelled VC-only soil pixel map. No pixel-scale bootstrap has
    # been run; the basin confirmation shown in the subtitle is a distinct test.
    .soil_plot <- data.table::copy(results)
    .soil_plot[, vc_pixel_result := factor(ifelse(is.finite(p_fdr_vc) & p_fdr_vc < 0.05,
                                                  "Pixel VC BH-FDR < 0.05", "Not VC BH-FDR significant / filtered"),
                                           levels = c("Not VC BH-FDR significant / filtered", "Pixel VC BH-FDR < 0.05"))]
    .soil_fig <- ggplot2::ggplot(.soil_plot, ggplot2::aes(lon, lat, colour = vc_pixel_result)) +
      ggplot2::geom_point(shape = 15, size = 1.1) + ggplot2::coord_fixed() +
      ggplot2::scale_colour_manual(values = c("grey80", "#b2182b"), name = "Pixel VC-only") +
      ggplot2::labs(
        title = sprintf("%s: pixel VC-HR98 only (BH-FDR)", index_label),
        subtitle = sprintf("Basin VC p=%.3g; basin bootstrap p=%.3g; basin dual-confirmed: %s. Pixel tests have NO bootstrap confirmation.",
                           bfull$p_vc, bfull$mbb_p, if (bfull$robust_significant) "YES" else "NO"),
        x = "Longitude", y = "Latitude") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(legend.position = "bottom")
    ggplot2::ggsave(file.path(TREND_DIR, sprintf("%s_%02d_pixel_VC_ONLY.png", idx, sc)),
                    .soil_fig, width = 8, height = 5.3, dpi = 220)
  }
  cat(sprintf("  ✓ %s | VC%s FDR-sig %d / %d | TFPW raw BH flags (NOT valid if k>=3) %d / %d\n", basename(out_file),
              if (idx %in% c("ssmi", "ssi_sm")) "-ONLY (soil pixels; basin confirmation separate)" else "",
              sum(results$p_fdr_vc < 0.05, na.rm = TRUE), sum(is.finite(results$p_value_vc)),
              sum(results$p_fdr_tfpw < 0.05, na.rm = TRUE), sum(is.finite(results$p_value_tfpw))))
  
  dc_cols <- c("lon", "lat", "n_events_D12", "n_events_D36", "n_events_D712", "n_events_D13p", "n_D4p", "mean_I", "mean_S")
  data.table::fwrite(results[, ..dc_cols], file.path(TREND_DIR, sprintf("%s_%02d_duration_class_map.csv", idx, sc)))
  
  write_d13p_nc <- function(values_vec) {
    r_out <- r_tmpl; all_vals <- rep(NA_real_, terra::ncell(r_tmpl)); all_vals[valid_pix] <- values_vec
    terra::values(r_out) <- all_vals
    names(r_out) <- sprintf("%s_%02d_n_D13p_longterm_droughts", idx, sc)
    terra::writeCDF(r_out, file.path(TREND_DIR, sprintf("%s_%02d_D13p_pixel_frequency.nc", idx, sc)), overwrite = TRUE)
  }
  write_d13p_nc(results$n_events_D13p)
  
  req <- c(out_file, basin_avg_path, basin_cov_path,
           file.path(TREND_DIR, sprintf("%s_%02d_duration_class_map.csv", idx, sc)),
           file.path(TREND_DIR, sprintf("%s_%02d_D13p_pixel_frequency.nc", idx, sc)),
           basin_series_regional_file(idx, sc, "monthly", variant),
           file.path(extent_out, sprintf("%s_basin_extent_monthly.csv", index_label)),
           if (RUN_CALENDAR_MONTH) cm_file)
  miss <- req[!file.exists(req) | is.na(file.info(req)$size) | file.info(req)$size <= 0]
  if (length(miss)) stop("Required output validation failed for ", index_label, ":\n", paste(miss, collapse = "\n"))
  if (SAVE_TS_RDS)
    saveRDS(list(ts_mat = ts_mat, lon = lon_vec, lat = lat_vec, cell = cell_vec, dates = dates_ts,
                 record_start_year = min(years)),
            file.path(TREND_DIR, sprintf("%s_%02d_timeseries.rds", idx, sc)), compress = "gzip")
  rm(ts_mat, results, cm_all); invisible(gc(verbose = FALSE))
  invisible(TRUE)
}


# ==============================================================================
#   9. SDI — HYDROLOGICAL-YEAR SERIES
# ==============================================================================
# One value per hydrological year and pixel (Oct start). Same QA + VC-HR98 + TFPW + BH-FDR,
# basin/regional trends, homogeneity (k4 = water year) and points. Monthly event / runs /
# spectral statistics do not apply. 'scale' = cumulative months from October (3, 6, 9, 12).
process_sdi <- function(P) {
  P <- as.integer(P); k <- match(P, SDI_PERIODS)
  f <- list.files(SDI_DIR, pattern = sprintf("^sdi_k%d_.*\\.nc$", k), full.names = TRUE)
  if (length(f) != 1L) stop(sprintf("Expected one sdi_k%d_*.nc in %s; found %d", k, SDI_DIR, length(f)))
  cat(sprintf("\n╔══ SDI k%d (%d months from October) ══╗  %s\n", k, P, basename(f)))
  r <- terra::rast(f); r_tmpl <- r[[1]]
  assert_same_native_grid(r_tmpl, PRECIP_TEMPLATE_13, paste0("SDI k", k, " vs precipitation grid"))
  hy_end <- as.Date(terra::time(r))
  if (anyNA(hy_end)) stop("SDI layers carry no time stamps: ", basename(f))
  source_years <- as.integer(format(hy_end, "%Y"))
  if (k == 1L && !all(format(hy_end, "%m") == "12"))
    stop("SDI k1 must have December timestamps before water-year relabelling.")
  yrs <- source_years + as.integer(k == 1L)          # water year = ending year
  if (any(diff(yrs) != 1L)) stop("SDI k", k, " is not a run of consecutive hydrological years.")
  smat <- as.matrix(terra::values(r, mat = TRUE)); smat[!is.finite(smat)] <- NA_real_
  bg <- get_basin_grid_cached(r_tmpl)
  valid <- bg$cover_all > 0 & rowSums(is.finite(smat)) > 0
  if (!any(valid)) stop("No valid basin pixels for SDI k", k)
  xy <- terra::xyFromCell(r_tmpl, which(valid)); smat <- smat[valid, , drop = FALSE]; cell_vec <- which(valid)
  
  # Producer series: read by the unmodified timestamps and verified on the FULL record,
  # then relabelled to water years. Trimming happens only after verification.
  stored <- read_basin_series("sdi", P, SDI_DIR, times = source_years, resolution = "annual")
  if (VERIFY_BASIN_SERIES)
    verify_basin_series(stored, smat, r_tmpl, cells = cell_vec, label = sprintf("SDI-%02d", P))
  if ("hydro_year_end" %in% names(stored$regional))
    stored$regional[, hydro_year_end := yrs[match(hydro_year_end, source_years)]]
  
  keep_y <- if (SDI_LIMIT_TO_RECORD_END) yrs <= ANALYSIS_PERIOD_END else rep(TRUE, length(yrs))
  if (!any(keep_y)) stop("No SDI k", k, " water year ends by ", ANALYSIS_PERIOD_END)
  n_drop <- sum(!keep_y)
  smat <- smat[, keep_y, drop = FALSE]; yrs <- yrs[keep_y]; hy_end <- hy_end[keep_y]
  basin_value <- stored$basin$value[keep_y]; basin_cov <- stored$basin$coverage[keep_y]
  stored$regional <- stored$regional[hydro_year_end %in% yrs]
  cal_start <- sprintf("Oct %d", min(yrs) - 1L)
  cal_end   <- sprintf("%s %d", month.abb[((8L + P) %% 12L) + 1L], max(yrs) - 1L + (8L + P) %/% 12L)
  
  res <- data.table::data.table(cell = cell_vec, lon = xy[, 1], lat = xy[, 2], n_years = rowSums(is.finite(smat)))
  res[, `:=`(qa_full = .qa_rows(smat, MIN_N_YEARS), modal_tie_fraction = .modal_frac_rows(smat), rho1_detrended = .rho1_rows(smat))]
  res <- cbind(res, .mk_pixel(smat, res$qa_full == "ok"))
  res[, `:=`(tfpw_method = "Yue2002_shared_J4_J6", p_fdr_vc = .bh(p_value_vc), p_fdr_tfpw = .bh(p_value_tfpw),
             n_years_sdi_lt_m1 = rowSums(smat < -1, na.rm = TRUE), n_years_sdi_lt_m1p5 = rowSums(smat < -1.5, na.rm = TRUE),
             driest_year = yrs[apply(smat, 1L, function(z) if (all(is.na(z))) NA_integer_ else which.min(z))],
             first_year = min(yrs), last_year = max(yrs),
             calendar_start = cal_start, calendar_end = cal_end,
             water_year_convention = "WY Y = Oct Y-1 to Sep Y (labelled by ending year)",
             period_label = sub("^sdi_k\\d+_(.*)\\.nc$", "\\1", basename(f)), slope_unit = "per hydrological year")]
  .record_qa("sdi", P, "annual_hydro_year", res$qa_full)
  .record_qa("sdi", P, "vc_status_annual_hydro_year", data.table::fifelse(is.na(res$vc_status), "not_tested", res$vc_status))
  
  bdt <- data.table::data.table(index_type = "sdi", scale = P, hydro_year_end = yrs,
                                period_end_date = hy_end, value = basin_value,
                                valid_area_coverage_fraction = basin_cov)
  .acc("basin_annual", bdt)
  bsdi <- .basin_full_row(bdt$value, "sdi", P, "annual_hydro_year", "per hydrological year", 5L)
  .acc("basin_full", bsdi)
  
  reg_split <- split(stored$regional[order(region_id, hydro_year_end)], by = "region_id")
  reg_vals  <- lapply(reg_split, function(d) d$value)
  reg_names <- vapply(reg_split, function(d) as.character(d$region_name[1]), character(1))
  if (!BASIN_RID %in% names(reg_vals)) { reg_vals[[BASIN_RID]] <- bdt$value; reg_names[BASIN_RID] <- "Nechako Basin" }
  
  for (rid in names(reg_vals))
    .acc("regional_full", cbind(data.table::data.table(region_id = rid, region_name = reg_names[[rid]], index_type = "sdi",
                                                       scale = P, series = "annual_hydro_year", slope_unit = "per hydrological year"),
                                if (rid == BASIN_RID) bsdi[, ..MK_SERIES_COLS] else .mk_series(reg_vals[[rid]], MIN_N_YEARS)))
  
  if (RUN_HOMOGENEITY && P == 12L)
    run_homog_spec(HOMOG_SPECS[index_type == "sdi"][1], smat, xy[, 1], xy[, 2], reg_vals, reg_names, list(annual = TRUE, years = yrs))
  
  pr <- .point_rows(r_tmpl, cell_vec)
  if (!is.null(pr) && nrow(pr)) for (j in seq_len(nrow(pr))) {
    rr <- pr$row[j]
    .acc("points", cbind(data.table::data.table(point_id = pr$point_id[j], lon_wgs84 = pr$lon_wgs84[j], lat_wgs84 = pr$lat_wgs84[j],
                                                index_type = "sdi", scale = P, period = "annual_hydro_year"),
                         res[rr, .(cell, lon, lat, n_years, qa_full, tau_vc, p_value_vc, sl_vc, tau_tfpw, p_value_tfpw, sl_tfpw, p_fdr_vc, p_fdr_tfpw)]))
    .acc("point_ts", data.table::data.table(point_id = pr$point_id[j], index_type = "sdi", scale = P, date = hy_end, water_year = yrs, value = smat[rr, ]))
  }
  out_file <- file.path(TREND_DIR, sprintf("sdi_%02d_results.csv", P))
  data.table::fwrite(res, out_file, showProgress = FALSE)
  cat(sprintf("  ✓ %s | %d pixels x %d water years (WY %d-%d = %s to %s)%s | VC FDR-sig %d | TFPW FDR-sig %d\n",
              basename(out_file), nrow(res), length(yrs), min(yrs), max(yrs), cal_start, cal_end,
              if (n_drop) sprintf(" | %d WY after %d dropped", n_drop, ANALYSIS_PERIOD_END) else "",
              sum(res$p_fdr_vc < 0.05, na.rm = TRUE), sum(res$p_fdr_tfpw < 0.05, na.rm = TRUE)))
  invisible(TRUE)
}

# ==============================================================================
#   10. RUN ALL INDICES
# ==============================================================================
cat("\n  Indices : SPI, SPEI-PM, SRI (", paste(SPI_SCALES, collapse = ", "), "), SWEI-3, SSPI-1, SDI k1-k4",
    if (RUN_SOIL_MOISTURE) sprintf(", SSMI / SSI-SM %s (%s)", paste(SSMI_SCALES, collapse = ", "), SOIL_LAYER) else "",
    "\n", sep = "")
cat(sprintf("  Events  : onset < %.2f, termination >= %.2f, min duration (months) by scale: %s\n",
            DROUGHT_ONSET, DROUGHT_END,
            paste(sprintf("%s=%d", names(EVENT_MIN_DURATION_BY_SCALE), EVENT_MIN_DURATION_BY_SCALE), collapse = ", ")))
total_start <- Sys.time()
RUN_STATUS <- data.table::data.table(index_type = character(), scale = integer(), success = logical(), message = character())
.status <- function(idx, sc, ok, msg) RUN_STATUS <<- data.table::rbindlist(list(RUN_STATUS,
                                                                                data.table::data.table(index_type = idx, scale = as.integer(sc), success = ok, message = msg)), use.names = TRUE)

run_index_safe <- function(index_type, scales, find_fn, seas_dir) {
  for (sc in scales) {
    msg <- "completed"
    ok <- tryCatch({ process_index(index_type, sc, find_fn, seas_dir); TRUE },
                   error = function(e) { msg <<- conditionMessage(e)
                   cat(sprintf("\n  ❌ %s-%02d failed: %s\n", INDEX_LABEL[[tolower(index_type)]], sc, msg)); FALSE })
    .status(tolower(index_type), sc, ok, msg)
  }
}
.missing_dir <- function(idx, scales, dir, hint) {
  cat(sprintf("\n── %s skipped: directory not found: %s ──\n", INDEX_LABEL[[tolower(idx)]], dir))
  for (sc in scales) .status(idx, sc, FALSE, paste("directory not found;", hint))
}

cat("\n── SPI ──\n");      run_index_safe("spi", SPI_SCALES, find_seasonal_nc_files, SPI_SEAS_DIR)
cat("\n── SPEI (PM) ──\n"); run_index_safe("spei", SPEI_SCALES, find_seasonal_nc_files, SPEI_SEAS_DIR)
if (dir.exists(SWEI_SEAS_DIR)) { cat("\n── SWEI ──\n"); run_index_safe("swei", SWEI_SCALE, find_swei_seasonal_files, SWEI_SEAS_DIR)
} else .missing_dir("swei", SWEI_SCALE, SWEI_SEAS_DIR, "run 4snow_drought_indices_ERALand.R")
if (dir.exists(SSPI_SEAS_DIR)) { cat("\n── SSPI ──\n"); run_index_safe("sspi", SSPI_SCALE, find_seasonal_nc_files, SSPI_SEAS_DIR)
} else .missing_dir("sspi", SSPI_SCALE, SSPI_SEAS_DIR, "run 4snow_drought_indices_ERALand.R")
if (dir.exists(SRI_SEAS_DIR)) { cat("\n── SRI (total runoff) ──\n"); run_index_safe("sri", SRI_SCALES, find_seasonal_nc_files, SRI_SEAS_DIR)
} else .missing_dir("sri", SRI_SCALES, SRI_SEAS_DIR, "run 5SRI_SDI_ERALand.R")
if (dir.exists(SDI_DIR)) {
  cat("\n── SDI (hydrological year) ──\n")
  for (P in SDI_PERIODS) {
    msg <- "completed"
    ok <- tryCatch({ process_sdi(P); TRUE }, error = function(e) { msg <<- conditionMessage(e)
    cat(sprintf("\n  ❌ SDI-%02d failed: %s\n", P, msg)); FALSE })
    .status("sdi", P, ok, msg)
  }
} else .missing_dir("sdi", SDI_PERIODS, SDI_DIR, "run 5SRI_SDI_ERALand.R")
if (RUN_SOIL_MOISTURE) {
  if (dir.exists(SSMI_SEAS_DIR)) { cat("\n── SSMI (soil moisture, ", SOIL_LAYER, ") ──\n", sep = "")
    run_index_safe("ssmi", SSMI_SCALES, find_seasonal_nc_files, SSMI_SEAS_DIR)
  } else .missing_dir("ssmi", SSMI_SCALES, SSMI_SEAS_DIR, "run 8_SSI_SSMI_ERALand.R")
  if (dir.exists(SSI_SM_SEAS_DIR)) { cat("\n── SSI-SM (soil moisture, ", SOIL_LAYER, ") ──\n", sep = "")
    run_index_safe("ssi_sm", SSI_SM_SCALES, find_seasonal_nc_files, SSI_SM_SEAS_DIR)
  } else .missing_dir("ssi_sm", SSI_SM_SCALES, SSI_SM_SEAS_DIR, "run 8_SSI_SSMI_ERALand.R")
}


# The MK worker pool is no longer needed for the combined output tables.
# Close it before the correlation export and restore serial defaults.
if (!is.null(cl)) {
  parallel::stopCluster(cl)
  cl <- NULL
}
N_CORES <- 1L

# ==============================================================================
#   11. COMBINED BASIN / REGIONAL / HOMOGENEITY / POINT / QA TABLES
# ==============================================================================
cat("\n══ Writing combined tables ══\n")
.bind <- function(slot) data.table::rbindlist(ACC[[slot]], use.names = TRUE, fill = TRUE)
.save <- function(dt, path) if (!is.null(dt) && nrow(dt)) { data.table::fwrite(dt, path); cat("  ✓", sub(paste0(TREND_DIR, "/"), "", path, fixed = TRUE), "\n") }

bc <- .bind("basin_calmonth"); if (nrow(bc)) bc[, p_vc_bh_months := .bh(p_vc), by = .(index_type, scale)]
rc <- .bind("regional_calmonth"); if (nrow(rc)) rc[, p_vc_bh_months := .bh(p_vc), by = .(region_id, index_type, scale)]
.save(.bind("basin_full"),    file.path(BASIN_TREND_DIR, "basin_trend_full_series.csv"))
.save(.bind("snow_basin_season"), file.path(BASIN_TREND_DIR, "snow_basin_season_trends.csv"))
.save(.bind("snow_basin_annual"), file.path(BASIN_TREND_DIR, "snow_basin_annual_season_means.csv"))
.save(bc,                     file.path(BASIN_TREND_DIR, "basin_trend_calendar_month.csv"))
.save(.bind("basin_monthly"), file.path(BASIN_TREND_DIR, "basin_index_monthly_all.csv"))
.save(.bind("basin_annual"),  file.path(BASIN_TREND_DIR, "basin_sdi_annual_all.csv"))
.save(.bind("regional_full"), file.path(REGIONAL_TREND_DIR, "regional_trend_full_series.csv"))
.save(rc,                     file.path(REGIONAL_TREND_DIR, "regional_trend_calendar_month.csv"))

hr <- .bind("homog_region")
if (nrow(hr)) {
  hr[, family := unname(INDEX_FAMILY[index_type])]
  .save(hr, file.path(HOMOG_DIR, "step_tests_by_region_all.csv"))
  hb <- data.table::copy(hr[region_id == BASIN_RID])
  hb[, `:=`(kind = "index", pettitt_p_bh_across_series = .bh(pettitt_p), step_p_bh_across_series = .bh(step_p))]
  .save(analysis_homogeneity_basin_summary(hb), file.path(HOMOG_DIR, "homogeneity_break_summary_basin.csv"))
}
.save(.bind("homog_subperiod"), file.path(HOMOG_DIR, "subperiod_trends_all.csv"))
.save(.bind("homog_series"),    file.path(HOMOG_DIR, "homogeneity_annual_series_all.csv"))
.save(.bind("homog_pixel"),     file.path(HOMOG_DIR, "pixel_step_summary.csv"))
.save(.bind("points"),          file.path(POINT_DIR, "point_trend_stats_drought.csv"))
.save(.bind("point_ts"),        file.path(POINT_DIR, "point_index_timeseries.csv"))
qa <- .bind("qa")
if (nrow(qa)) {
  .save(qa, file.path(QA_DIR, "filter_reasons_by_index_scale.csv"))
  .save(qa[qa != "ok", .(n = sum(n)), by = .(analysis, qa)][order(analysis, -n)], file.path(QA_DIR, "root_cause_filter_reasons.csv"))
}
saveRDS(list(step_year = STEP_YEAR, min_n_monthly = MIN_N_MONTHLY, min_n_years = MIN_N_YEARS, tie_mode_max = TIE_MODE_MAX,
             mbb_B = MBB_B, mbb_block_rule = "max(12, 2k) months; 5 years for SDI", homog_mc = HOMOG_MC,
             soil_layer = SOIL_LAYER, soil_scales = as.integer(SSMI_SCALES), run_soil_moisture = RUN_SOIL_MOISTURE,
             homog_specs = HOMOG_SPECS[, .(index_type, scale, months = vapply(months, paste, "", collapse = ","), suffix, desc)],
             sdi_water_year_convention = "WY Y = Oct Y-1 to Sep Y", sdi_limit_to_record_end = SDI_LIMIT_TO_RECORD_END,
             dfreq_min_months = DFREQ_MIN_MONTHS, dfreq_snow_months = DFREQ_SNOW_MONTHS,
             dfreq_definition_version = DFREQ_DEFINITION_VERSION,
             dfreq_requires_complete_seasons = TRUE,
             snow_basin_windows = SNOW_BASIN_WINDOWS,
             snow_season_mbb_block_years = SNOW_SEASON_MBB_YEARS,
             soil_pixel_inference = "VC-HR98 BH-FDR pixel only; basin soil requires VC and MBB",
             basin_region_id = BASIN_RID, created = Sys.time()),
        file.path(TREND_DIR, "script13_analysis_config.rds"))


# ==============================================================================
#   12. RUN QA
# ==============================================================================
expected_runs <- data.table::rbindlist(list(
  data.table::data.table(index_type = "spi", scale = SPI_SCALES),
  data.table::data.table(index_type = "spei", scale = SPEI_SCALES),
  data.table::data.table(index_type = "swei", scale = SWEI_SCALE),
  data.table::data.table(index_type = "sspi", scale = SSPI_SCALE),
  data.table::data.table(index_type = "sri", scale = SRI_SCALES),
  data.table::data.table(index_type = "sdi", scale = SDI_PERIODS),
  if (RUN_SOIL_MOISTURE) data.table::data.table(index_type = "ssmi", scale = SSMI_SCALES),
  if (RUN_SOIL_MOISTURE) data.table::data.table(index_type = "ssi_sm", scale = SSI_SM_SCALES)))
run_status_path <- file.path(TREND_DIR, "script13_run_status.csv")
data.table::fwrite(RUN_STATUS, run_status_path)
missing_runs <- expected_runs[!RUN_STATUS, on = .(index_type, scale)]
failed_runs  <- RUN_STATUS[success == FALSE]
cat("\n══════════════════════════════════════════════════════\nSCRIPT 13 EXECUTION QA\n",
    "══════════════════════════════════════════════════════\n", sep = "")
cat(sprintf("  Expected %d | successful %d | failed %d | missing %d\n  Run-status file: %s\n",
            nrow(expected_runs), sum(RUN_STATUS$success), nrow(failed_runs), nrow(missing_runs), run_status_path))
if (nrow(failed_runs)) print(failed_runs)
if (nrow(missing_runs)) print(missing_runs)
run_status_ok <- !nrow(failed_runs) && !nrow(missing_runs) && nrow(RUN_STATUS) == nrow(expected_runs)
if (run_status_ok){
  cat(sprintf("  ✓ All %d expected index x scale combinations completed.\n", nrow(expected_runs)))
} else {
  cat("  ⚠ Run-status QA failed; correlation exports will still be attempted before the final stop.\n")
}

# ==============================================================================
#   13. PIXEL-TO-BASIN CORRELATION FILES (read by Scripts 14 and 15)
# ==============================================================================
# Pearson and Spearman correlation of each pixel's index series with the
# producer's basin series (verified against the pixels in step [3/9]).
# Writes {idx}_{sc}_pixel_basin_correlation.csv into each producer folder.
cat("\n══ Pixel-to-basin correlation files ══\n")
if (!SAVE_TS_RDS) stop("Pixel-to-basin correlations need SAVE_TS_RDS = TRUE.")
CORR_MIN_PAIRS <- 60L
CORR_REQUIRED  <- c("spi", "spei")
CORR_SCALES <- list(spi = SPI_SCALES, spei = SPEI_SCALES,
                    swei = SWEI_SCALE, sspi = SSPI_SCALE, sri = SRI_SCALES)

pixel_basin_corr <- function(ts_mat, basin, min_pairs) {
  out <- matrix(NA_real_, nrow(ts_mat), 3L,
                dimnames = list(NULL, c("pearson_r", "spearman_rho", "n_pairs")))
  bf <- is.finite(basin)
  for (i in seq_len(nrow(ts_mat))) {
    ok <- bf & is.finite(ts_mat[i, ])
    out[i, 3L] <- sum(ok)
    if (sum(ok) < min_pairs) next
    x <- ts_mat[i, ok]; y <- basin[ok]
    if (stats::sd(x) == 0 || stats::sd(y) == 0) next
    out[i, 1L] <- stats::cor(x, y)
    out[i, 2L] <- stats::cor(x, y, method = "spearman")
  }
  out
}

CORR_STATUS <- list()
for (idx in names(CORR_SCALES)) for (sc in CORR_SCALES[[idx]]) {
  lab <- sprintf("%s-%02d", INDEX_LABEL[[idx]], sc)
  rds <- file.path(TREND_DIR, sprintf("%s_%02d_timeseries.rds", idx, sc))
  out <- file.path(BASIN_SERIES_DIRS[[idx]], sprintf("%s_%02d_pixel_basin_correlation.csv", idx, sc))
  CORR_STATUS[[length(CORR_STATUS) + 1L]] <- tryCatch({
    # A failed index run may leave a previous run's RDS on disk. Never publish
    # correlations from stale data when run-status QA reports this pair failed.
    current_run <- RUN_STATUS[index_type == idx & scale == sc]
    if (nrow(current_run) != 1L || !isTRUE(current_run$success[[1L]]))
      stop("Index x scale did not complete in this run; stale cached series rejected")
    tsd <- readRDS(rds)
    if (ncol(tsd$ts_mat) != length(tsd$dates)) stop("matrix/date mismatch in ", basename(rds))
    dates_corr <- as.Date(tsd$dates)
    month_corr <- as.integer(format(dates_corr, "%m"))
    snow_months <- switch(idx, swei = 1:4, sspi = c(12L, 1L, 2L, 3L), 1:12)
    use_month <- month_corr %in% snow_months
    pixel_results <- data.table::fread(file.path(TREND_DIR, sprintf("%s_%02d_results.csv", idx, sc)),
                                       select = c("cell", "qa_full"))
    if (anyDuplicated(pixel_results$cell) || !setequal(pixel_results$cell, tsd$cell))
      stop("Trend-QA cells do not match pixel time-series cache: ", idx, "-", sc)
    qa_full <- pixel_results$qa_full[match(tsd$cell, pixel_results$cell)]
    if (idx %in% c("swei", "sspi")) {
      # One common, pre-specified seasonal window for ALL pixels; calendar
      # summer gaps cannot alter a pixel's selection of correlation months.
      ts_snow <- tsd$ts_mat[, use_month, drop = FALSE]
      qa_snow <- .qa_rows(ts_snow, MIN_N_YEARS * length(snow_months))
      eligible <- qa_snow == "ok"
      qa_basis <- "supported_snow_season"
    } else {
      qa_snow <- rep(NA_character_, length(qa_full))
      eligible <- qa_full == "ok"
      qa_basis <- "full_monthly"
    }
    basin <- as.numeric(read_basin_series(idx, sc, BASIN_SERIES_DIRS[[idx]], times = dates_corr)$basin$value)
    if (length(basin) != length(dates_corr)) stop("Basin correlation calendar mismatch: ", idx)
    cc <- pixel_basin_corr(tsd$ts_mat[, use_month, drop = FALSE], basin[use_month], CORR_MIN_PAIRS)
    cc[!eligible, 1:2] <- NA_real_
    dt <- data.table::data.table(cell = as.integer(tsd$cell), x = tsd$lon, y = tsd$lat,
                                 pearson_r = cc[, 1L], spearman_rho = cc[, 2L],
                                 n_pairs = as.integer(cc[, 3L]), qa_full = qa_full,
                                 qa_snow_season = qa_snow, qa_basis = qa_basis,
                                 correlation_eligible = eligible,
                                 months_used = paste(month.abb[sort(unique(snow_months))], collapse = ","))
    data.table::fwrite(dt, out)
    v <- dt[is.finite(spearman_rho)]
    cat(sprintf("  ✓ %-12s %4d pixels (%4d valid) | median r %.3f | median rho %.3f\n",
                lab, nrow(dt), nrow(v), stats::median(v$pearson_r), stats::median(v$spearman_rho)))
    data.table::data.table(index_type = idx, scale = sc, ok = TRUE, n_pixels = nrow(dt), n_valid = nrow(v),
                           median_pearson = stats::median(v$pearson_r),
                           median_spearman = stats::median(v$spearman_rho), message = basename(out))
  }, error = function(e) {
    cat(sprintf("  ❌ %-12s %s\n", lab, conditionMessage(e)))
    data.table::data.table(index_type = idx, scale = sc, ok = FALSE, message = conditionMessage(e))
  })
}
CORR_STATUS <- data.table::rbindlist(CORR_STATUS, fill = TRUE)
data.table::fwrite(CORR_STATUS, file.path(QA_DIR, "pixel_basin_correlation_status.csv"))
CORR_BAD <- CORR_STATUS[ok == FALSE & index_type %in% CORR_REQUIRED]
cat(sprintf("  ✓ %d of %d correlation files written.\n", sum(CORR_STATUS$ok), nrow(CORR_STATUS)))
if (!run_status_ok || nrow(CORR_BAD)) {
  issues <- c(
    if (!run_status_ok) sprintf("incomplete index x scale runs (see %s)", run_status_path),
    if (nrow(CORR_BAD)) paste0("required SPI/SPEI correlation exports failed: ",
                               paste(sprintf("%s-%02d", INDEX_LABEL[CORR_BAD$index_type], CORR_BAD$scale), collapse = ", "))
  )
  stop("SCRIPT 13 FAILED QA: ", paste(issues, collapse = "; "))
}

# ==============================================================================
#   14. SUMMARY
# ==============================================================================
cat("\n══════════════════════════════════════════════════════\nSCRIPT 13 SUMMARY\n",
    "══════════════════════════════════════════════════════\n",
    "A full monthly series : VC-HR98 primary; TFPW k>=3 REJECTED (diagnostic only); basin block-bootstrap MK\n",
    "B calendar-month      : VC-HR98 primary; TFPW retained as diagnostic, BH-FDR within month\n",
    "C drought-month frac. : VC-HR98; 12/12 calendar or 4/4 snow-season months; >12 exploratory\n",
    "D events / PELT / runs / red-noise spectra / basin extent\n",
    "E basin + regional    : A + B for every region; snow fixed-season annual basin sensitivity\n",
    sprintf("F homogeneity         : %d annual-resolution series; step year %d\n", nrow(HOMOG_SPECS), STEP_YEAR),
    sprintf("G QA                  : min %d months / %d years; ties > %.0f%%; internal gaps never spliced\n",
            MIN_N_MONTHLY, MIN_N_YEARS, 100 * TIE_MODE_MAX),
    sprintf("Total elapsed time    : %.1f minutes\n", as.numeric(difftime(Sys.time(), total_start, units = "mins"))), sep = "")