# ==============================================================================
#   11trend_test_ALL.R 
#   SPATIAL TREND, EVENT AND HOMOGENEITY ANALYSIS FOR ALL DROUGHT INDICES
# ==============================================================================
# INDICES : SPI, SPEI-PM, SPEI-ERA5 and SRI (1, 3, 6, 12, 24, 36 months), SWEI-3,
#           SSPI-1, SDI (hydrological-year periods k1-k4, Oct start), and the
#           soil-moisture indices SSMI and SSI-SM (1, 2, 3 months; layer
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
#      (sensitivity only). Basin series: moving-block-bootstrap MK, block =
#      max(12, 2k) months, as an assumption-light check.
#   B. Calendar-month series (n ~ 76 per month) = Script 9's (9pr_pet_trends.r) monthly analysis:
#      VC-MK + TFPW-MK per pixel, BH-FDR within index x scale x month.
#   C. Annual drought-month fraction (months < DROUGHT_ONSET per year) = Script 9's
#      annual-aggregate analogue: VC-MK per pixel.
#   D. Events: ONE framework-wide definition from DROUGHT_ANALYSIS_utils.R
#      (hysteresis run method, onset < DROUGHT_ONSET, termination >= DROUGHT_END,
#      scale-dependent minimum duration EVENT_MIN_DURATION_BY_SCALE). PELT, runs
#      test, red-noise spectra and basin drought extent use the same events.
#   E. Basin and regional series: A + B (no gap splicing).
#   F. Homogeneity battery from Script 9 (Pettitt, SNHT, Buishand range, Welch at
#      the a-priori step year, AR(1)-GLS trend vs trend+step, sub-period trends,
#      pixel Pettitt scans with BH-FDR) ONLY on non-overlapping annual-resolution
#      series (HOMOG_SPECS): these tests assume independent observations.
#   G. Series QA (Script 9 rules): too few values, internal gaps, constant,
#      ties (largest tied group > TIE_MODE_MAX). Filtered series: NA statistics +
#      reason code. SWEI summer gaps are never spliced.
#   H. Point trends at Script 9's specific points.
#
# OUTPUTS (TREND_DIR)
#   {idx}_{sc}_results.csv                 full-series pixel results (+ events, QA, dfreq)
#   {idx}_{sc}_calmonth_results.csv        calendar-month pixel results
#   sdi_{P}_results.csv                    SDI pixel results (annual)
#   basin_trends/ regional_trends/ regional_timeseries/ homogeneity/ points/ qa/
#   duration-class CSVs, D13+ NetCDFs, basin_extent/, script6_run_status.csv,
#   script6_analysis_config.rds
#
# EXECUTION ORDER: index scripts 1, 2b, 3, 4, 5, 8 (and 6, 7) -> 11 (this script) -> 12 -> 13 -> 14
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
utils_load_packages(c("terra", "data.table", "Kendall", "trend", "modifiedmk", "lubridate",
                      "changepoint", "zoo", "ggplot2", "dplyr", "patchwork", "nlme"))

# The shared MK / spectral wrappers reference these; Script 11 runs serially.
N_CORES <- 1L
is_windows <- identical(.Platform$OS.type, "windows")
cl <- NULL


RUN_CALENDAR_MONTH  <- TRUE     # B
RUN_DROUGHT_FREQ    <- TRUE     # C
RUN_HOMOGENEITY     <- TRUE     # F
RUN_BLOCK_BOOTSTRAP <- TRUE     # A (basin series)
MBB_B               <- 1000L    # block-bootstrap replicates
HOMOG_MC            <- 2000L    # Monte Carlo replicates for SNHT / Buishand p-values
MIN_N_MONTHLY       <- 120L     # min finite months for a full-series test
MIN_N_YEARS         <- 30L      # min finite years for calendar-month / annual tests
TIE_MODE_MAX        <- 0.50     # filter when the largest tied group exceeds this share
DFREQ_MIN_MONTHS    <- 6L       # min finite months in a year for the drought-month fraction
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

# Step year and specific points are taken from Script 9 so both pipelines agree.
S4_META_FILE <- file.path(WD_PATH, "trend_analysis_pr_pet", "analysis_metadata.rds")
S4_META <- if (file.exists(S4_META_FILE)) tryCatch(readRDS(S4_META_FILE), error = function(e) NULL) else NULL
STEP_YEAR <- local({
  sy <- tryCatch(S4_META$swe_config$step_year, error = function(e) NULL)
  if (length(sy) == 1L && is.finite(sy)) as.integer(sy) else 1979L
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
    "SCRIPT 6 v2 — NATIVE ERA5-LAND GRID | step year ", STEP_YEAR,
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
    stop(label, " does not match the native grid of the first file. Script 11 never resamples.")
  invisible(TRUE)
}

native_grid_key <- function(r) paste(
  terra::nrow(r), terra::ncol(r), paste(signif(as.vector(terra::ext(r)), 12), collapse = ","),
  trimws(terra::crs(r)), paste(signif(terra::origin(r), 12), collapse = ","), sep = "|")

REGIONAL_WEIGHT_CACHE <- new.env(parent = emptyenv())
get_regional_weights_cached <- function(r_tmpl) {
  key <- native_grid_key(r_tmpl)
  if (exists(key, envir = REGIONAL_WEIGHT_CACHE, inherits = FALSE)) return(get(key, envir = REGIONAL_WEIGHT_CACHE))
  regions <- load_nechako_regions(template = r_tmpl)
  validate_nechako_regions(regions)
  rw <- region_fractional_weights(r_tmpl, regions)
  assign(key, rw, envir = REGIONAL_WEIGHT_CACHE)
  rw
}

BASIN_GRID_CACHE <- new.env(parent = emptyenv())
get_basin_grid_cached <- function(r_tmpl) {
  key <- native_grid_key(r_tmpl)
  if (exists(key, envir = BASIN_GRID_CACHE, inherits = FALSE)) return(get(key, envir = BASIN_GRID_CACHE))
  basin_v <- nechako_subbasin_union(load_nechako_subbasins())   # 7-subbasin union, as in Script 3
  if (!terra::same.crs(basin_v, r_tmpl)) basin_v <- terra::project(basin_v, terra::crs(r_tmpl))
  cover_all <- as.numeric(terra::values(terra::rasterize(basin_v, r_tmpl, cover = TRUE), mat = FALSE))
  cover_all[!is.finite(cover_all)] <- 0
  if (any(cover_all < -1e-6 | cover_all > 1 + 1e-6)) stop("Basin coverage outside [0, 1].")
  cover_all <- pmax(0, pmin(1, cover_all))
  area_all <- as.numeric(terra::values(terra::cellSize(r_tmpl, unit = "m"), na.rm = FALSE))
  a_r <- sum(area_all * cover_all, na.rm = TRUE)
  a_v <- tryCatch(sum(as.numeric(terra::expanse(basin_v, unit = "m", transform = TRUE))), error = function(e) NA_real_)
  if (is.finite(a_v) && a_v > 0) {
    rel <- abs(a_r - a_v) / a_v
    cat(sprintf("  Basin area check: raster %.1f km² | polygon %.1f km² | diff %.2f%%\n", a_r / 1e6, a_v / 1e6, 100 * rel))
    if (rel > 0.05) stop(sprintf("Rasterized basin area differs from polygon area by %.2f%%.", 100 * rel))
  }
  out <- list(cover_all = cover_all, area_all = area_all)
  assign(key, out, envir = BASIN_GRID_CACHE)
  out
}

regional_mean_with_coverage <- function(ts_matrix, area_weights, total_region_area_m2,
                                        min_coverage = MIN_REGIONAL_VALID_AREA_COVERAGE) {
  if (!is.matrix(ts_matrix) || length(area_weights) != nrow(ts_matrix)) stop("regional_mean_with_coverage: dimension mismatch.")
  w <- as.numeric(area_weights)
  if (!is.finite(total_region_area_m2) || total_region_area_m2 <= 0)
    return(list(value = rep(NA_real_, ncol(ts_matrix)), coverage = rep(NA_real_, ncol(ts_matrix))))
  coverage <- value <- rep(NA_real_, ncol(ts_matrix))
  for (j in seq_len(ncol(ts_matrix))) {
    ok <- is.finite(ts_matrix[, j]) & is.finite(w) & w > 0
    va <- sum(w[ok]); coverage[j] <- va / total_region_area_m2
    if (va > 0 && coverage[j] >= min_coverage) value[j] <- sum(ts_matrix[ok, j] * w[ok]) / va
  }
  list(value = value, coverage = coverage)
}


# ==============================================================================
#   2. INDICES, SCALES, DIRECTORIES, HOMOGENEITY SERIES
# ==============================================================================
SPI_SCALES <- SPEI_SCALES <- SPEI_ERA5_SCALES <- SRI_SCALES <- c(1L, 3L, 6L, 12L, 24L, 36L)
SWEI_SCALES <- 3L
SSPI_SCALES <- 1L
SDI_PERIODS <- c(3L, 6L, 9L, 12L)     # cumulative months from October: k1..k4

SPEI_ERA5_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal_era5land")
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
# Scales without a framework minimum event duration (EVENT_MIN_DURATION_BY_SCALE)
# are dropped here rather than failing inside process_index().
.soil_scale_ok <- function(s) !is.null(tryCatch(get_event_min_duration(s), error = function(e) NULL))
SSMI_SCALES <- SSI_SM_SCALES <- Filter(.soil_scale_ok, c(1L, 2L, 3L))
if (RUN_SOIL_MOISTURE && !identical(SSMI_SCALES, c(1L, 2L, 3L)))
  cat("  NOTE: soil-moisture scales limited to", paste(SSMI_SCALES, collapse = ", "),
      "(no EVENT_MIN_DURATION_BY_SCALE entry for the others)\n")

INDEX_LABEL  <- c(spi = "SPI", spei = "SPEI-PM", spei_era5 = "SPEI-ERA5", swei = "SWEI",
                  sspi = "SSPI", sri = "SRI", sdi = "SDI", ssmi = "SSMI", ssi_sm = "SSI-SM")
INDEX_FAMILY <- c(spi = "meteorological", spei = "meteorological", spei_era5 = "meteorological",
                  swei = "snow", sspi = "snow", sri = "runoff", sdi = "runoff",
                  ssmi = "soil_moisture", ssi_sm = "soil_moisture")

# Non-overlapping annual-resolution series for the homogeneity battery.
.hs <- function(idx, sc, months, suffix, desc)
  data.table::data.table(index_type = idx, scale = as.integer(sc), months = list(as.integer(months)),
                         suffix = suffix, desc = desc)
HOMOG_SPECS <- data.table::rbindlist(c(
  unlist(lapply(c("spi", "spei", "spei_era5", "sri"), function(i) list(
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
.modal_frac  <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(tabulate(match(x, unique(x)))) / length(x) }
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
.sen_rows <- function(mat) apply(mat, 1L, function(x) {
  t <- which(is.finite(x)); if (length(t) < 3L) return(NA_real_)
  v <- x[t]; d <- outer(v, v, "-"); dt <- outer(t, t, "-"); k <- lower.tri(d)
  stats::median(d[k] / dt[k])
})

.mk_na <- function(n) data.table::data.table(tau_vc = rep(NA_real_, n), p_value_vc = NA_real_, sl_vc = NA_real_,
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
  if (!"sl_vc" %in% names(r)) r[, sl_vc := .sen_rows(sub)]
  r[match(seq_len(n), vi)]
}

# Single-series VC-HR98 (basin / regional / sub-period).
.mk_series <- function(x, min_n) {
  reason <- .qa_reason(x, min_n)
  out <- data.table::data.table(n = sum(is.finite(x)), qa = reason, tau = NA_real_, p_vc = NA_real_,
                                sen = NA_real_, n_over_neff = NA_real_, p_mk_plain = NA_real_, rho1 = NA_real_)
  if (reason != "ok") return(out)
  v <- .trim_finite(x)
  hr <- tryCatch(modifiedmk::mmkh(v), error = function(e) NULL)
  if (!is.null(hr)) out[, `:=`(tau = unname(hr["Tau"]), p_vc = unname(hr["new P-value"]),
                               sen = unname(hr["Sen's slope"]), n_over_neff = unname(hr["N/N*"]),
                               p_mk_plain = unname(hr["old P.value"]))]
  out[, rho1 := .rho1_one(v)]
  out
}
.tfpw_na <- function() data.table::data.table(tau_tfpw = NA_real_, p_tfpw = NA_real_, sen_tfpw = NA_real_)
.tfpw_series <- function(v) {
  r <- tryCatch(data.table::as.data.table(vectorized_mann_kendall_tfpw(matrix(v, nrow = 1L))), error = function(e) NULL)
  if (is.null(r) || !all(c("tau_tfpw", "p_value_tfpw", "sl_tfpw") %in% names(r))) return(.tfpw_na())
  data.table::data.table(tau_tfpw = r$tau_tfpw[1], p_tfpw = r$p_value_tfpw[1], sen_tfpw = r$sl_tfpw[1])
}

# Moving-block-bootstrap MK: resampling blocks keeps the short-range dependence and destroys any trend.
.mk_S <- function(x) {
  x <- x[is.finite(x)]; if (length(x) < 3L) return(NA_real_)
  d <- outer(x, x, "-"); sum(sign(d[lower.tri(d)]))
}
.mbb_mk_p <- function(x, L, B, seed) {
  n <- length(x); S0 <- .mk_S(x)
  if (!is.finite(S0) || n < 3L * L) return(c(S = S0, p = NA_real_))
  starts <- seq_len(n - L + 1L); nb <- ceiling(n / L)
  set.seed(seed)
  Sb <- vapply(seq_len(B), function(b) {
    s <- sample(starts, nb, replace = TRUE)
    .mk_S(x[as.vector(outer(0:(L - 1L), s, "+"))][seq_len(n)])
  }, numeric(1))
  c(S = S0, p = (1 + sum(abs(Sb) >= abs(S0), na.rm = TRUE)) / (1 + sum(is.finite(Sb))))
}

.basin_full_row <- function(x, idx, sc, series, unit, L) {
  r <- .mk_series(x, if (series == "full_monthly") MIN_N_MONTHLY else MIN_N_YEARS)
  ok <- r$qa == "ok"; v <- .trim_finite(x)
  do_mbb <- ok && RUN_BLOCK_BOOTSTRAP
  mb <- if (do_mbb) .mbb_mk_p(v, L, MBB_B, seed = 7919L * as.integer(sc) + sum(utf8ToInt(idx))) else c(S = NA_real_, p = NA_real_)
  cbind(data.table::data.table(index_type = idx, scale = as.integer(sc), series = series, slope_unit = unit), r,
        if (ok) .tfpw_series(v) else .tfpw_na(),
        data.table::data.table(mbb_S = unname(mb["S"]), mbb_p = unname(mb["p"]),
                               mbb_block = if (do_mbb) as.integer(L) else NA_integer_,
                               mbb_B = if (do_mbb) MBB_B else NA_integer_))
}

ACC <- new.env()
for (s in c("basin_full", "basin_calmonth", "basin_monthly", "basin_annual", "regional_full",
            "regional_calmonth", "homog_region", "homog_subperiod", "homog_series", "homog_pixel",
            "points", "point_ts", "qa")) ACC[[s]] <- list()
.acc <- function(slot, x) ACC[[slot]][[length(ACC[[slot]]) + 1L]] <- x
.record_qa <- function(idx, sc, analysis, reason) {
  tb <- table(reason)
  .acc("qa", data.table::data.table(index_type = idx, scale = as.integer(sc), analysis = analysis,
                                    qa = names(tb), n = as.integer(tb)))
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

n_sim_spectral <- 500L
SPECTRAL_RN_ENV_CACHE <- new.env(parent = emptyenv())
vectorized_spectral_peaks <- function(ts_matrix, min_period_months = 24L, sig_level = 0.95, n_sim_spectral = 500L) {
  n_pix <- nrow(ts_matrix); out <- rep(NA_integer_, n_pix)
  get_thr <- function(n_len, band, rho1) {
    rb <- max(-0.95, min(0.95, round(rho1 / 0.05) * 0.05))
    key <- paste0(n_len, "_r", sprintf("%.2f", rb))
    if (!is.null(SPECTRAL_RN_ENV_CACHE[[key]])) return(SPECTRAL_RN_ENV_CACHE[[key]])
    had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had) old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({ if (had) assign(".Random.seed", old, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv) }, add = TRUE)
    set.seed(42L)
    mx <- vapply(seq_len(n_sim_spectral), function(s) {
      innov <- stats::rnorm(n_len + 500L, sd = sqrt(max(1e-8, 1 - rb^2)))
      rn <- utils::tail(as.numeric(stats::filter(innov, filter = rb, method = "recursive")), n_len)
      rn <- stats::residuals(stats::lm(rn ~ seq_along(rn)))
      vv <- stats::var(rn); if (!is.finite(vv) || vv <= 0) return(NA_real_)
      sp <- stats::spectrum(rn / sqrt(vv), spans = c(3L, 5L), taper = 0, detrend = FALSE, plot = FALSE)
      max(sp$spec[band], na.rm = TRUE)
    }, numeric(1))
    mx <- mx[is.finite(mx)]
    if (!length(mx)) stop("Red-noise simulation produced no valid spectral maxima.")
    thr <- stats::quantile(mx, probs = sig_level, names = FALSE)
    SPECTRAL_RN_ENV_CACHE[[key]] <- thr
    thr
  }
  for (i in seq_len(n_pix)) {
    x <- .trim_finite(as.numeric(ts_matrix[i, ]))
    if (anyNA(x) || length(x) < 48L) next
    x <- stats::residuals(stats::lm(I(x - mean(x)) ~ seq_along(x)))
    pv <- stats::var(x); if (!is.finite(pv) || pv <= 0) next
    rho1 <- tryCatch(stats::cor(x[-length(x)], x[-1]), error = function(e) 0)
    if (!is.finite(rho1)) rho1 <- 0
    sp <- tryCatch(stats::spectrum(x, spans = c(3L, 5L), taper = 0, detrend = FALSE, plot = FALSE), error = function(e) NULL)
    if (is.null(sp)) next
    band <- (1 / sp$freq) >= min_period_months
    if (!any(band)) next
    sig_bins <- band & (sp$spec > get_thr(length(x), band, max(-0.95, min(0.95, rho1))) * pv)
    if (any(sig_bins)) { r <- rle(sig_bins); out[i] <- sum(r$values) }
  }
  data.table::data.table(n_spectral_peaks = out)
}


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
    dplyr::summarise(mean_pct = mean(pct_area, na.rm = TRUE), .groups = "drop")
  xc <- extent_annual$mean_pct[is.finite(extent_annual$mean_pct)]
  mk_ext <- list(tau = NA_real_, p_value = NA_real_, sens_slope = NA_real_)
  if (length(xc) >= 10L) {
    hr <- tryCatch(modifiedmk::mmkh(xc), error = function(e) NULL)
    tau_val <- as.numeric(Kendall::MannKendall(xc)$tau)
    mk_ext <- if (!is.null(hr)) list(tau = tau_val, p_value = unname(hr["new P-value"]), sens_slope = unname(hr["Sen's slope"]))
    else list(tau = tau_val, p_value = as.numeric(Kendall::MannKendall(xc)$sl), sens_slope = NA_real_)
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
#   6. HOMOGENEITY BATTERY (Script 9 equivalent, annual-resolution series only)
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

.pettitt_fast <- function(x) {
  n <- length(x); r <- rank(x); t <- seq_len(n - 1L)
  U <- 2 * cumsum(r)[t] - t * (n + 1); K <- max(abs(U))
  list(tau = which.max(abs(U)), p = min(1, 2 * exp(-6 * K^2 / (n^3 + n^2))))
}

.gls_ar1 <- function(f, d) tryCatch(suppressWarnings(
  nlme::gls(f, data = d, correlation = nlme::corAR1(form = ~ tt), method = "ML")), error = function(e) NULL)

.homog_na <- function(n = NA_integer_) data.table::data.table(
  n_years = as.integer(n), first_year = NA_integer_, last_year = NA_integer_, pettitt_p = NA_real_,
  pettitt_first_year_after = NA_real_, snht_p = NA_real_, snht_first_year_after = NA_real_,
  buishand_range_p = NA_real_, mean_before = NA_real_, mean_after = NA_real_, diff_at_step = NA_real_,
  welch_p = NA_real_, trend_no_step = NA_real_, trend_no_step_p = NA_real_, trend_with_step = NA_real_,
  trend_with_step_p = NA_real_, step_coef = NA_real_, step_p = NA_real_,
  delta_AIC_step_vs_nostep = NA_real_, ar1_phi_with_step = NA_real_)

.homog_battery <- function(v, yrs, sy) {
  ok <- which(is.finite(v)); if (!length(ok)) return(.homog_na(0L))
  i <- min(ok):max(ok); v <- v[i]; y <- yrs[i]
  if (anyNA(v) || length(v) < MIN_N_YEARS || length(unique(v)) < 3L) return(.homog_na(sum(is.finite(v))))
  out <- .homog_na(length(v)); out[, `:=`(first_year = min(y), last_year = max(y))]
  pt <- tryCatch(trend::pettitt.test(v), error = function(e) NULL)
  if (!is.null(pt)) { k <- as.integer(pt$estimate)[1]; out[, `:=`(pettitt_p = pt$p.value, pettitt_first_year_after = y[k] + 1)] }
  sn <- tryCatch(trend::snh.test(v, m = HOMOG_MC), error = function(e) NULL)
  if (!is.null(sn)) { k <- as.integer(sn$estimate)[1]; out[, `:=`(snht_p = sn$p.value, snht_first_year_after = y[k] + 1)] }
  br <- tryCatch(trend::br.test(v, m = HOMOG_MC), error = function(e) NULL)
  if (!is.null(br)) out[, buishand_range_p := br$p.value]
  b <- v[y < sy]; a <- v[y >= sy]
  if (length(a) >= 5L && length(b) >= 5L)
    out[, `:=`(mean_before = mean(b), mean_after = mean(a), diff_at_step = mean(a) - mean(b),
               welch_p = tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_))]
  d <- data.frame(v = v, tt = y - min(y), s = as.integer(y >= sy))
  g0 <- .gls_ar1(v ~ tt, d)
  if (!is.null(g0)) { t0 <- summary(g0)$tTable; out[, `:=`(trend_no_step = t0["tt", "Value"], trend_no_step_p = t0["tt", "p-value"])] }
  if (length(unique(d$s)) == 2L) {
    g1 <- .gls_ar1(v ~ tt + s, d)
    if (!is.null(g1)) {
      t1 <- summary(g1)$tTable
      out[, `:=`(trend_with_step = t1["tt", "Value"], trend_with_step_p = t1["tt", "p-value"],
                 step_coef = t1["s", "Value"], step_p = t1["s", "p-value"],
                 ar1_phi_with_step = as.numeric(stats::coef(g1$modelStruct$corStruct, unconstrained = FALSE)))]
      if (!is.null(g0)) out[, delta_AIC_step_vs_nostep := stats::AIC(g1) - stats::AIC(g0)]
    }
  }
  out
}

.subperiod_trends <- function(v, yrs, sy) {
  win <- list(full = range(yrs), before_step = c(min(yrs), sy - 1L), from_step = c(sy, max(yrs)),
              wmo_1991_2020 = c(1991L, 2020L))
  data.table::rbindlist(lapply(names(win), function(w) {
    k <- yrs >= win[[w]][1] & yrs <= win[[w]][2]
    r <- .mk_series(v[k], 10L)
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
  pc <- .pettitt_fast(v); b <- v[y < sy]; a <- v[y >= sy]
  data.table::data.table(n_years = length(v), pettitt_p = pc$p, pettitt_first_year_after = y[pc$tau] + 1,
                         mean_before = if (length(b)) mean(b) else NA_real_,
                         diff_at_step = if (length(a) && length(b)) mean(a) - mean(b) else NA_real_,
                         welch_p = if (length(a) >= 5L && length(b) >= 5L)
                           tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_) else NA_real_)
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
  px <- data.table::rbindlist(lapply(seq_len(nrow(s$mat)), function(i) .pixel_step(s$mat[i, ], y, STEP_YEAR)))
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
  index_label <- sprintf("%s-%02d", toupper(index_type), sc)
  min_dur <- get_event_min_duration(sc)
  cat(sprintf("\n╔══ %s ══╗  (event minimum duration %d months)\n", index_label, min_dur))
  
  # ---- [1/9] files and native grid ------------------------------------------------------------
  cat("[1/9] Files and native grid...\n")
  files <- find_fn(seas_dir, index_type, sc)
  if (!length(files)) stop(sprintf("No monthly NetCDF files for %s in %s", index_label, seas_dir))
  stacks <- lapply(files, function(f) tryCatch(terra::rast(f), error = function(e) {
    cat("  ⚠ Failed: ", basename(f), " | ", conditionMessage(e), "\n", sep = ""); NULL }))
  bad <- vapply(stacks, is.null, logical(1))
  if (any(bad)) stop(sprintf("Failed to load %d NetCDF(s) for %s: %s", sum(bad), index_label, paste(basename(files[bad]), collapse = ", ")))
  r_tmpl <- stacks[[1]][[1]]
  if (is.na(terra::crs(r_tmpl)) || !nzchar(terra::crs(r_tmpl))) stop("Native index grid has no CRS: ", basename(files[1]))
  for (k in seq_along(stacks)) assert_same_native_grid(stacks[[k]], r_tmpl, sprintf("%s file %d (%s)", index_label, k, basename(files[k])))
  cat(sprintf("  ✓ %d files | %d x %d cells | res %.6f x %.6f\n", length(files), terra::nrow(r_tmpl),
              terra::ncol(r_tmpl), terra::res(r_tmpl)[1], terra::res(r_tmpl)[2]))
  xy_all <- terra::xyFromCell(r_tmpl, seq_len(terra::ncell(r_tmpl)))
  
  # ---- [2/9] time-series matrix ---------------------------------------------------------------
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
  cat("[3/9] Regional and basin series and trends...\n")
  rw_all <- get_regional_weights_cached(r_tmpl)
  regional_index <- data.table::rbindlist(lapply(names(rw_all), function(rid) {
    rw <- rw_all[[rid]]; w_all <- as.numeric(rw$weights)
    tot <- sum(w_all[is.finite(w_all) & w_all > 0], na.rm = TRUE)
    st <- regional_mean_with_coverage(ts_mat, w_all[valid_pix], tot)
    data.table::data.table(region_id = rid, region_name = rw$region_name[1], date = dates_ts, value = st$value,
                           valid_area_coverage_fraction = st$coverage, coverage_area_m2 = tot,
                           n_cells = rw$n_cells, n_partial_cells = rw$n_partial)
  }), use.names = TRUE, fill = TRUE)
  regional_index[, `:=`(index_type = idx, scale = sc)]
  cat_file <- file.path(REGIONAL_INDEX_DIR, "region_catalog.csv")
  if (!file.exists(cat_file)) data.table::fwrite(region_catalog(load_nechako_regions(template = r_tmpl)), cat_file)
  data.table::fwrite(regional_index, file.path(REGIONAL_INDEX_DIR, sprintf("%s_%02d_regional_averaged_by_month.csv", idx, sc)))
  
  basin_stats <- regional_mean_with_coverage(ts_mat, area_valid, sum(eff_area_all, na.rm = TRUE))
  basin_avg_vec <- as.numeric(basin_stats$value)
  .acc("basin_monthly", data.table::data.table(index_type = idx, scale = sc, date = dates_ts, value = basin_avg_vec,
                                               valid_area_coverage_fraction = basin_stats$coverage))
  
  reg_split <- split(regional_index[order(region_id, date)], by = "region_id")
  reg_vals  <- lapply(reg_split, function(d) d$value)
  reg_names <- vapply(reg_split, function(d) as.character(d$region_name[1]), character(1))
  if (!BASIN_RID %in% names(reg_vals)) { reg_vals[[BASIN_RID]] <- basin_avg_vec; reg_names[BASIN_RID] <- "Nechako Basin" }
  
  for (rid in names(reg_vals)) {
    x <- reg_vals[[rid]]
    .acc("regional_full", cbind(data.table::data.table(region_id = rid, region_name = reg_names[[rid]], index_type = idx,
                                                       scale = sc, series = "full_monthly", slope_unit = "per month"),
                                .mk_series(x, MIN_N_MONTHLY)))
    .acc("regional_calmonth", data.table::rbindlist(lapply(1:12, function(m)
      cbind(data.table::data.table(region_id = rid, region_name = reg_names[[rid]], index_type = idx, scale = sc,
                                   month = m, slope_unit = "per year"), .mk_series(x[mo_ts == m], MIN_N_YEARS)))))
  }
  bfull <- .basin_full_row(basin_avg_vec, idx, sc, "full_monthly", "per month", max(12L, 2L * sc))
  .acc("basin_full", bfull)
  .acc("basin_calmonth", data.table::rbindlist(lapply(1:12, function(m) {
    x <- basin_avg_vec[mo_ts == m]; r <- .mk_series(x, MIN_N_YEARS)
    cbind(data.table::data.table(index_type = idx, scale = sc, month = m, slope_unit = "per year"), r,
          if (r$qa == "ok") .tfpw_series(.trim_finite(x)) else .tfpw_na())
  })))
  cat(sprintf("  ✓ Basin full series: tau %+.3f | p VC %.3g | p TFPW %.3g | p block-bootstrap %.3g | N/N* %.1f\n",
              bfull$tau, bfull$p_vc, bfull$p_tfpw, bfull$mbb_p, bfull$n_over_neff))
  
  # ---- [4/9] full-series pixel trends ---------------------------------------------------------
  cat("[4/9] Full-series pixel trends (VC-HR98 primary; TFPW sensitivity)...\n")
  results <- data.table::data.table(cell = cell_vec, lon = lon_vec, lat = lat_vec, n_months = rowSums(is.finite(ts_mat)))
  results[, qa_full := .qa_rows(ts_mat, MIN_N_MONTHLY)]
  results[, `:=`(modal_tie_fraction = .modal_frac_rows(ts_mat), rho1_detrended = .rho1_rows(ts_mat))]
  results <- cbind(results, .mk_pixel(ts_mat, results$qa_full == "ok"))
  results[, tfpw_method := "Yue2002_shared_J4"]
  .record_qa(idx, sc, "full_monthly", results$qa_full)
  cat(sprintf("  ✓ tested %d / %d pixels | mean detrended lag-1 rho %.2f\n",
              sum(results$qa_full == "ok"), n_pix, mean(results$rho1_detrended, na.rm = TRUE)))
  
  # ---- [5/9] events, PELT, runs, spectra, drought-month fraction ------------------------------
  cat("[5/9] Events, PELT, runs test, spectra, drought-month fraction...\n")
  results <- cbind(results, vectorized_event_stats(ts_mat, min_dur, DROUGHT_ONSET, DROUGHT_END))
  cat(sprintf("  ✓ events (onset < %.2f, end >= %.2f, min %d months): mean %.1f per pixel (D13+ max %d)\n",
              DROUGHT_ONSET, DROUGHT_END, min_dur, mean(results$n_events), max(results$n_events_D13p)))
  results <- cbind(results, vectorized_regime_shift_pelt(ts_mat, years))
  results <- cbind(results, vectorized_runs_test(ts_mat, min_dur, DROUGHT_ONSET, DROUGHT_END))
  results <- cbind(results, vectorized_spectral_peaks(ts_mat, n_sim_spectral = n_sim_spectral))
  
  if (RUN_DROUGHT_FREQ) {
    uy <- sort(unique(yr_ts))
    fr <- matrix(vapply(uy, function(y) {
      sub <- ts_mat[, yr_ts == y, drop = FALSE]; nf <- rowSums(is.finite(sub))
      ifelse(nf >= DFREQ_MIN_MONTHS, rowSums(sub < DROUGHT_ONSET, na.rm = TRUE) / nf, NA_real_)
    }, numeric(n_pix)), nrow = n_pix)
    q <- .qa_rows(fr, MIN_N_YEARS); mk <- .mk_pixel(fr, q == "ok")
    results[, `:=`(dfreq_mean = rowMeans(fr, na.rm = TRUE), dfreq_qa = q, dfreq_tau_vc = mk$tau_vc,
                   dfreq_p_vc = mk$p_value_vc, dfreq_sl_vc = mk$sl_vc)]
    results[, dfreq_p_fdr_vc := .bh(dfreq_p_vc)]
    .record_qa(idx, sc, "drought_month_fraction", q)
  }
  results[, `:=`(p_fdr_vc = .bh(p_value_vc), p_fdr_tfpw = .bh(p_value_tfpw))]
  results[, `:=`(p_fdr_vc_raw = p_fdr_vc, p_fdr_tfpw_raw = p_fdr_tfpw)]
  
  # ---- [6/9] calendar-month pixel trends ------------------------------------------------------
  cm_all <- NULL
  cm_file <- file.path(TREND_DIR, sprintf("%s_%02d_calmonth_results.csv", idx, sc))
  if (RUN_CALENDAR_MONTH) {
    cat("[6/9] Calendar-month pixel trends (VC-HR98 + TFPW, BH-FDR within month)...\n")
    cm_all <- data.table::rbindlist(lapply(1:12, function(m) {
      sub <- ts_mat[, mo_ts == m, drop = FALSE]; q <- .qa_rows(sub, MIN_N_YEARS)
      d <- cbind(data.table::data.table(cell = cell_vec, lon = lon_vec, lat = lat_vec, month = m,
                                        n_years = rowSums(is.finite(sub)), qa = q), .mk_pixel(sub, q == "ok"))
      d[, `:=`(p_fdr_vc = .bh(p_value_vc), p_fdr_tfpw = .bh(p_value_tfpw))]
      d
    }), use.names = TRUE, fill = TRUE)
    cm_all[, `:=`(index_type = idx, timescale = sc, slope_unit = "per year")]
    .record_qa(idx, sc, "calendar_month", cm_all$qa)
    data.table::fwrite(cm_all, cm_file)
    s <- cm_all[, .(tested = sum(qa == "ok"), fdr = sum(p_fdr_vc < 0.05, na.rm = TRUE)), by = month]
    cat("  ✓ VC FDR-significant pixels by month: ",
        paste(sprintf("%s %d/%d", month.abb[s$month], s$fdr, s$tested), collapse = " | "), "\n", sep = "")
  }
  
  # ---- [7/9] homogeneity on annual-resolution series ------------------------------------------
  if (RUN_HOMOGENEITY) {
    specs <- HOMOG_SPECS[index_type == idx & scale == sc]
    if (nrow(specs)) {
      cat("[7/9] Homogeneity battery...\n")
      for (q in seq_len(nrow(specs)))
        run_homog_spec(specs[q], ts_mat, lon_vec, lat_vec, reg_vals, reg_names, list(annual = FALSE, dates = dates_ts))
    }
  }
  
  # ---- [8/9] basin CSVs, extent, points -------------------------------------------------------
  cat("[8/9] Basin CSVs, extent, specific points...\n")
  seas_dir_map <- list(spi = SPI_SEAS_DIR, spei = SPEI_SEAS_DIR, spei_era5 = SPEI_ERA5_SEAS_DIR,
                       swei = SWEI_SEAS_DIR, sspi = SSPI_SEAS_DIR, sri = SRI_SEAS_DIR,
                       ssmi = SSMI_SEAS_DIR, ssi_sm = SSI_SM_SEAS_DIR)
  basin_csv_dir <- seas_dir_map[[idx]]
  if (is.null(basin_csv_dir) || !nzchar(basin_csv_dir)) basin_csv_dir <- file.path(WD_PATH, paste0(idx, "_results_seasonal"))
  dir.create(basin_csv_dir, showWarnings = FALSE, recursive = TRUE)
  uy <- sort(unique(yr_ts))
  to_wide <- function(values) {
    w <- data.frame(Year = uy)
    for (m in 1:12) { col <- rep(NA_real_, length(uy)); hit <- which(mo_ts == m); col[match(yr_ts[hit], uy)] <- values[hit]; w[[month.abb[m]]] <- col }
    w
  }
  basin_avg_path <- file.path(basin_csv_dir, sprintf("%s_%02d_basin_averaged_by_month.csv", idx, sc))
  basin_cov_path <- file.path(basin_csv_dir, sprintf("%s_%02d_basin_valid_area_coverage_by_month.csv", idx, sc))
  utils::write.csv(to_wide(basin_avg_vec), basin_avg_path, row.names = FALSE, na = "")
  utils::write.csv(to_wide(basin_stats$coverage), basin_cov_path, row.names = FALSE, na = "")
  
  extent_out <- file.path(TREND_DIR, "basin_extent")
  compute_basin_extent_and_event_catalog(ts_mat, n_years, index_label, extent_out, lon_vec, lat_vec,
                                         start_year = min(years), area_weights = area_valid, min_duration = min_dur)
  
  pr <- .point_rows(r_tmpl, cell_vec)
  if (!is.null(pr) && nrow(pr)) {
    keep_cols <- intersect(c("cell", "lon", "lat", "n_months", "qa_full", "rho1_detrended", "tau_vc", "p_value_vc",
                             "sl_vc", "tau_tfpw", "p_value_tfpw", "sl_tfpw", "p_fdr_vc", "p_fdr_tfpw",
                             "dfreq_tau_vc", "dfreq_p_vc", "dfreq_sl_vc"), names(results))
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
  cat("[9/9] Writing pixel outputs...\n")
  out_file <- file.path(TREND_DIR, sprintf("%s_%02d_results.csv", idx, sc))
  data.table::fwrite(results, out_file, showProgress = FALSE)
  cat(sprintf("  ✓ %s | VC FDR-sig %d / %d | TFPW FDR-sig %d / %d\n", basename(out_file),
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
           file.path(REGIONAL_INDEX_DIR, sprintf("%s_%02d_regional_averaged_by_month.csv", idx, sc)),
           file.path(extent_out, sprintf("%s_basin_extent_monthly.csv", index_label)),
           if (RUN_CALENDAR_MONTH) cm_file)
  miss <- req[!file.exists(req) | is.na(file.info(req)$size) | file.info(req)$size <= 0]
  if (length(miss)) stop("Required output validation failed for ", index_label, ":\n", paste(miss, collapse = "\n"))
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
  hy_end <- as.Date(terra::time(r))
  if (anyNA(hy_end)) stop("SDI layers carry no time stamps: ", basename(f))
  yrs <- as.integer(format(hy_end, "%Y"))
  if (any(diff(yrs) != 1L)) stop("SDI k", k, " is not a run of consecutive hydrological years.")
  smat <- as.matrix(terra::values(r, mat = TRUE)); smat[!is.finite(smat)] <- NA_real_
  bg <- get_basin_grid_cached(r_tmpl)
  valid <- bg$cover_all > 0 & rowSums(is.finite(smat)) > 0
  if (!any(valid)) stop("No valid basin pixels for SDI k", k)
  xy <- terra::xyFromCell(r_tmpl, which(valid)); smat <- smat[valid, , drop = FALSE]; cell_vec <- which(valid)
  w_all <- bg$area_all * bg$cover_all
  
  res <- data.table::data.table(cell = cell_vec, lon = xy[, 1], lat = xy[, 2], n_years = rowSums(is.finite(smat)))
  res[, `:=`(qa_full = .qa_rows(smat, MIN_N_YEARS), modal_tie_fraction = .modal_frac_rows(smat), rho1_detrended = .rho1_rows(smat))]
  res <- cbind(res, .mk_pixel(smat, res$qa_full == "ok"))
  res[, `:=`(tfpw_method = "Yue2002_shared_J4", p_fdr_vc = .bh(p_value_vc), p_fdr_tfpw = .bh(p_value_tfpw),
             n_years_sdi_lt_m1 = rowSums(smat < -1, na.rm = TRUE), n_years_sdi_lt_m1p5 = rowSums(smat < -1.5, na.rm = TRUE),
             driest_year = yrs[apply(smat, 1L, function(z) if (all(is.na(z))) NA_integer_ else which.min(z))],
             first_year = min(yrs), last_year = max(yrs),
             period_label = sub("^sdi_k\\d+_(.*)\\.nc$", "\\1", basename(f)), slope_unit = "per hydrological year")]
  .record_qa("sdi", P, "annual_hydro_year", res$qa_full)
  
  bs <- regional_mean_with_coverage(smat, w_all[valid], sum(w_all, na.rm = TRUE))
  bdt <- data.table::data.table(index_type = "sdi", scale = P, hydro_year_end = yrs, value = bs$value,
                                valid_area_coverage_fraction = bs$coverage)
  data.table::fwrite(bdt[, .(hydro_year_end, value, valid_area_coverage_fraction)],
                     file.path(SDI_DIR, sprintf("sdi_%02d_basin_averaged_by_year.csv", P)))
  .acc("basin_annual", bdt)
  .acc("basin_full", .basin_full_row(bs$value, "sdi", P, "annual_hydro_year", "per hydrological year", 5L))
  
  rw_all <- get_regional_weights_cached(r_tmpl)
  reg_vals <- list(); reg_names <- character(0)
  for (rid in names(rw_all)) {
    wa <- as.numeric(rw_all[[rid]]$weights)
    reg_vals[[rid]] <- regional_mean_with_coverage(smat, wa[valid], sum(wa[is.finite(wa) & wa > 0], na.rm = TRUE))$value
    reg_names[rid] <- as.character(rw_all[[rid]]$region_name[1])
  }
  if (!BASIN_RID %in% names(reg_vals)) { reg_vals[[BASIN_RID]] <- bs$value; reg_names[BASIN_RID] <- "Nechako Basin" }
  data.table::fwrite(data.table::rbindlist(lapply(names(reg_vals), function(rid)
    data.table::data.table(region_id = rid, region_name = reg_names[[rid]], hydro_year_end = yrs, value = reg_vals[[rid]],
                           index_type = "sdi", scale = P))),
    file.path(REGIONAL_INDEX_DIR, sprintf("sdi_%02d_regional_averaged_by_year.csv", P)))
  for (rid in names(reg_vals))
    .acc("regional_full", cbind(data.table::data.table(region_id = rid, region_name = reg_names[[rid]], index_type = "sdi",
                                                       scale = P, series = "annual_hydro_year", slope_unit = "per hydrological year"),
                                .mk_series(reg_vals[[rid]], MIN_N_YEARS)))
  
  if (RUN_HOMOGENEITY && P == 12L)
    run_homog_spec(HOMOG_SPECS[index_type == "sdi"][1], smat, xy[, 1], xy[, 2], reg_vals, reg_names, list(annual = TRUE, years = yrs))
  
  pr <- .point_rows(r_tmpl, cell_vec)
  if (!is.null(pr) && nrow(pr)) for (j in seq_len(nrow(pr))) {
    rr <- pr$row[j]
    .acc("points", cbind(data.table::data.table(point_id = pr$point_id[j], lon_wgs84 = pr$lon_wgs84[j], lat_wgs84 = pr$lat_wgs84[j],
                                                index_type = "sdi", scale = P, period = "annual_hydro_year"),
                         res[rr, .(cell, lon, lat, n_years, qa_full, tau_vc, p_value_vc, sl_vc, tau_tfpw, p_value_tfpw, sl_tfpw, p_fdr_vc, p_fdr_tfpw)]))
    .acc("point_ts", data.table::data.table(point_id = pr$point_id[j], index_type = "sdi", scale = P, date = hy_end, value = smat[rr, ]))
  }
  out_file <- file.path(TREND_DIR, sprintf("sdi_%02d_results.csv", P))
  data.table::fwrite(res, out_file, showProgress = FALSE)
  cat(sprintf("  ✓ %s | %d pixels x %d years (%d-%d) | VC FDR-sig %d | TFPW FDR-sig %d\n", basename(out_file), nrow(res),
              length(yrs), min(yrs), max(yrs), sum(res$p_fdr_vc < 0.05, na.rm = TRUE), sum(res$p_fdr_tfpw < 0.05, na.rm = TRUE)))
  invisible(TRUE)
}


# ==============================================================================
#   10. RUN ALL INDICES
# ==============================================================================
cat("\n  Indices : SPI, SPEI-PM, SPEI-ERA5, SRI (", paste(SPI_SCALES, collapse = ", "), "), SWEI-3, SSPI-1, SDI k1-k4",
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
                   cat(sprintf("\n  ❌ %s-%02d failed: %s\n", toupper(index_type), sc, msg)); FALSE })
    .status(tolower(index_type), sc, ok, msg)
  }
}
.missing_dir <- function(idx, scales, dir, hint) {
  cat(sprintf("\n── %s skipped: directory not found: %s ──\n", toupper(idx), dir))
  for (sc in scales) .status(idx, sc, FALSE, paste("directory not found;", hint))
}

find_spei_era5_nc_safe <- function(data_dir, index_type, scale)
  tryCatch(suppressWarnings(find_seasonal_nc_files(data_dir, "spei", scale)), error = function(e) character(0))

cat("\n── SPI ──\n");      run_index_safe("spi", SPI_SCALES, find_seasonal_nc_files, SPI_SEAS_DIR)
cat("\n── SPEI (PM) ──\n"); run_index_safe("spei", SPEI_SCALES, find_seasonal_nc_files, SPEI_SEAS_DIR)
if (dir.exists(SPEI_ERA5_SEAS_DIR)) { cat("\n── SPEI (ERA5-Land PEV) ──\n")
  run_index_safe("spei_era5", SPEI_ERA5_SCALES, find_spei_era5_nc_safe, SPEI_ERA5_SEAS_DIR)
} else .missing_dir("spei_era5", SPEI_ERA5_SCALES, SPEI_ERA5_SEAS_DIR, "run 3SPEI_ERALand.R")
if (dir.exists(SWEI_SEAS_DIR)) { cat("\n── SWEI ──\n"); run_index_safe("swei", SWEI_SCALES, find_swei_seasonal_files, SWEI_SEAS_DIR)
} else .missing_dir("swei", SWEI_SCALES, SWEI_SEAS_DIR, "run 4snow_drought_indices_ERALand.R")
if (dir.exists(SSPI_SEAS_DIR)) { cat("\n── SSPI ──\n"); run_index_safe("sspi", SSPI_SCALES, find_seasonal_nc_files, SSPI_SEAS_DIR)
} else .missing_dir("sspi", SSPI_SCALES, SSPI_SEAS_DIR, "run 4snow_drought_indices_ERALand.R")
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


# ==============================================================================
#   11. COMBINED BASIN / REGIONAL / HOMOGENEITY / POINT / QA TABLES
# ==============================================================================
cat("\n══ Writing combined tables ══\n")
.bind <- function(slot) data.table::rbindlist(ACC[[slot]], use.names = TRUE, fill = TRUE)
.save <- function(dt, path) if (!is.null(dt) && nrow(dt)) { data.table::fwrite(dt, path); cat("  ✓", sub(paste0(TREND_DIR, "/"), "", path, fixed = TRUE), "\n") }

bc <- .bind("basin_calmonth"); if (nrow(bc)) bc[, p_vc_bh_months := .bh(p_vc), by = .(index_type, scale)]
rc <- .bind("regional_calmonth"); if (nrow(rc)) rc[, p_vc_bh_months := .bh(p_vc), by = .(region_id, index_type, scale)]
.save(.bind("basin_full"),    file.path(BASIN_TREND_DIR, "basin_trend_full_series.csv"))
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
  .save(hb, file.path(HOMOG_DIR, "homogeneity_break_summary_basin.csv"))
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
             homog_specs = HOMOG_SPECS[, .(index_type, scale, months = vapply(months, paste, "", collapse = ","), suffix, desc)],
             dfreq_min_months = DFREQ_MIN_MONTHS, basin_region_id = BASIN_RID, created = Sys.time()),
        file.path(TREND_DIR, "script6_analysis_config.rds"))


# ==============================================================================
#   12. RUN QA
# ==============================================================================
expected_runs <- data.table::rbindlist(list(
  data.table::data.table(index_type = "spi", scale = SPI_SCALES),
  data.table::data.table(index_type = "spei", scale = SPEI_SCALES),
  data.table::data.table(index_type = "spei_era5", scale = SPEI_ERA5_SCALES),
  data.table::data.table(index_type = "swei", scale = SWEI_SCALES),
  data.table::data.table(index_type = "sspi", scale = SSPI_SCALES),
  data.table::data.table(index_type = "sri", scale = SRI_SCALES),
  data.table::data.table(index_type = "sdi", scale = SDI_PERIODS),
  if (RUN_SOIL_MOISTURE) data.table::data.table(index_type = "ssmi", scale = SSMI_SCALES),
  if (RUN_SOIL_MOISTURE) data.table::data.table(index_type = "ssi_sm", scale = SSI_SM_SCALES)))
run_status_path <- file.path(TREND_DIR, "script6_run_status.csv")
data.table::fwrite(RUN_STATUS, run_status_path)
missing_runs <- expected_runs[!RUN_STATUS, on = .(index_type, scale)]
failed_runs  <- RUN_STATUS[success == FALSE]
cat("\n══════════════════════════════════════════════════════\nSCRIPT 6 EXECUTION QA\n",
    "══════════════════════════════════════════════════════\n", sep = "")
cat(sprintf("  Expected %d | successful %d | failed %d | missing %d\n  Run-status file: %s\n",
            nrow(expected_runs), sum(RUN_STATUS$success), nrow(failed_runs), nrow(missing_runs), run_status_path))
if (nrow(failed_runs)) print(failed_runs)
if (nrow(missing_runs)) print(missing_runs)
if (nrow(failed_runs) || nrow(missing_runs) || nrow(RUN_STATUS) != nrow(expected_runs))
  stop("SCRIPT 6 FAILED QA: not all ", nrow(expected_runs), " index x scale combinations completed. Inspect ", run_status_path)
cat(sprintf("  ✓ All %d expected index x scale combinations completed.\n", nrow(expected_runs)))


# ==============================================================================
#   13. SUMMARY
# ==============================================================================
cat("\n══════════════════════════════════════════════════════\nSCRIPT 6 v2 SUMMARY\n",
    "══════════════════════════════════════════════════════\n",
    "A full monthly series : VC-HR98 (primary), TFPW (sensitivity for k >= 3), basin block-bootstrap MK\n",
    "B calendar-month      : VC-HR98 + TFPW per pixel, BH-FDR within month\n",
    "C drought-month frac. : VC-HR98 per pixel\n",
    "D events / PELT / runs / red-noise spectra / basin extent\n",
    "E basin + regional    : A + B for every region\n",
    sprintf("F homogeneity         : %d annual-resolution series; step year %d\n", nrow(HOMOG_SPECS), STEP_YEAR),
    sprintf("G QA                  : min %d months / %d years; ties > %.0f%%; internal gaps never spliced\n",
            MIN_N_MONTHLY, MIN_N_YEARS, 100 * TIE_MODE_MAX),
    sprintf("Total elapsed time    : %.1f minutes\n", as.numeric(difftime(Sys.time(), total_start, units = "mins"))), sep = "")