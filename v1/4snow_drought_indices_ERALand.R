##############################################################################
# snow_drought_indices.R
#
# SIMPLE PIPELINE:  monthly ERA5-Land SWE  ->  two snow drought indices per pixel
#
#   1. SWEI    Standardized SWE Index  (Huning & AghaKouchak 2020)
#              k-month (default 3) right-aligned SWE sum -> for each calendar
#              month, nonparametric Gringorten plotting position -> N(0,1).
#
#   2. SSPI-1  Standardized SnowPack Index (zero-inflated gamma, after
#              Stagge et al. 2015). Monthly SWE -> for each calendar month,
#              P(zero) + gamma (L-moment fit) on the positive reference values
#              -> N(0,1).
#
# CALIBRATION: both SWEI and SSPI are standardized against the WHOLE
#   1950-01 to 2025-12 record (every year of each calendar month enters the
#   reference distribution). No 1991-2020 sub-period is used.
#
# INPUT  : one NetCDF of monthly ERA5-Land snow water equivalent
#          (monthly_data_direct/snow_depth_water_equivalent_monthly.nc)
# OUTPUT : snow_drought_results/ full-record SWEI/SSPI NetCDF stacks and QA products
#          swei_results_seasonal/ 12 month-of-year NetCDF files + basin-average CSV
#          sspi_results_seasonal/ 12 month-of-year NetCDF files + basin-average CSV
#          temporal_drought/regional_timeseries/ area-weighted basin/subbasin series
#
# TWO VALID-AREA-COVERAGE CASES (this version)
#   Regional (basin + every subbasin) area-weighted means are reported for two
#   values of MIN_VALID_AREA_COVERAGE, the minimum fraction of a region's area that
#   must hold a finite index value before a monthly mean is reported:
#     case "custom": MIN_VALID_AREA_COVERAGE (= BASIN_MIN_VALID_AREA_COVERAGE from
#                    DROUGHT_ANALYSIS_utils.R, i.e. 0.20; framework-wide rule)
#     case "zero"  : 0  (report whenever at least one finite cell exists)
#   - "custom" case keeps the ORIGINAL file names, so Script 6 / downstream code
#     is unaffected:  <index>_basin_averaged_by_month.csv,
#                     <index>_<scale>_regional_averaged_by_month.csv
#   - "zero" case adds the suffix _mincov0 to the same two files.
#   - NEW: <index>_<scale>_regional_averaged_by_month_coverage_comparison.csv
#          (both cases side by side, per region and month)
#   - NEW: snow_drought_regional_mean_summary.csv (written to REGIONAL_DIR):
#          mean index per region x case x calendar month (and all months).
#
# Framework contract: native precipitation grid; shared basin/subbasin definitions;
# area-weighted regional means; monthly seasonal NetCDF naming used by Script 6.
#
# Negative values = drier / lower snowpack than normal, positive = wetter.
#
# WHAT IS DELIBERATELY *NOT* IN HERE (compared with the 6-script pipeline):
#   - station bias correction, Sentinel-1 / C-SNOW depth fusion, temporal-
#     transfer ML, regional frequency analysis (pooled regional gamma),
#     trend screening and the extra validation diagnostics.
#   - the snow-cover-fraction (>= 5 %) domain mask. Only SWE is used here, so
#     the SWEI domain mask is the SWE-derived one the original script already
#     uses as its fallback: a pixel-month is scored only if k-month SWE is
#     non-zero in >= 75 % of years.
#
# Packages: terra, lmomco   (ncdf4 only needed if the NetCDF has no readable time)
##############################################################################

## ============================ CONFIGURATION ================================
# Follow the shared Nechako framework path contract used by Scripts 1 and 3.
WD_PATH <- Sys.getenv("NECHAKO_WD", unset = "D:/Nechako_Drought/Nechako/")
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH, call. = FALSE)
setwd(WD_PATH)

UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE)) stop("Shared utility file not found: ", UTILS_FILE, call. = FALSE)
source(UTILS_FILE)
utils_load_packages(c("terra", "ncdf4", "lmomco", "data.table"))
suppressPackageStartupMessages(library(data.table))   # dt[, ...] syntax used by the regional-mean summary

SWE_FILE <- file.path(WD_PATH, "monthly_data_direct", "snow_depth_water_equivalent_monthly.nc")
NATIVE_PRECIP_FILE <- file.path(WD_PATH, "monthly_data_direct", "total_precipitation_monthly.nc")
OUT_DIR <- file.path(WD_PATH, "snow_drought_results")
SWEI_DIR <- SWEI_SEAS_DIR
SSPI_DIR <- file.path(WD_PATH, "sspi_results_seasonal")
REGIONAL_DIR <- NECHAKO_REGION_OUTPUT_DIR

ANALYSIS_START <- as.Date("1950-01-01")
ANALYSIS_END   <- as.Date("2025-12-01")

# Calibration (reference) period for BOTH indices: the whole 1950-2025 record.
CALIB_START <- ANALYSIS_START
CALIB_END   <- ANALYSIS_END

SWEI_SCALE     <- 3L
SSPI_SCALE     <- 1L
SWEI_REFERENCE <- "whole_record"  # 1950-2025; documented in manifest
SSPI_REFERENCE <- "whole_record"  # 1950-2025; documented in manifest
SWE_UNITS      <- "m"             # ERA5-Land SWE is metres of water equivalent

# Minimum fraction of a region's area with a valid (finite) index value for its
# area-weighted mean to be reported. Cases compared in the regional outputs:
#   custom = MIN_VALID_AREA_COVERAGE (framework rule, 0.20)   |   zero = 0 (sensitivity only)
# The framework-wide value lives in DROUGHT_ANALYSIS_utils.R so every script
# (SPEI, SRI/SDI, trend tests, event ranking, ...) applies the same rule.
MIN_VALID_AREA_COVERAGE <- BASIN_MIN_VALID_AREA_COVERAGE
COVERAGE_ZERO           <- 0

MIN_NONZERO_FRACTION <- 0.75
MIN_GAMMA_N          <- 15L
MIN_POSITIVE_N       <- 10L
PROB_EPS             <- 1e-6
N_CORES              <- 1L

if (!identical(SWEI_REFERENCE, "whole_record") || !identical(SSPI_REFERENCE, "whole_record"))
  stop("SWEI_REFERENCE and SSPI_REFERENCE must both be 'whole_record' (1950-2025).", call. = FALSE)
if (!SWE_UNITS %in% c("m", "mm"))
  stop("SWE_UNITS must be 'm' or 'mm'. ERA5-Land defaults to metres.", call. = FALSE)
if (SWEI_SCALE < 1L || SSPI_SCALE < 1L)
  stop("Index accumulation scales must be positive integers.", call. = FALSE)
if (N_CORES < 1L) stop("N_CORES must be >= 1.", call. = FALSE)
if (!is.numeric(MIN_VALID_AREA_COVERAGE) || length(MIN_VALID_AREA_COVERAGE) != 1L ||
    !is.finite(MIN_VALID_AREA_COVERAGE) || MIN_VALID_AREA_COVERAGE < 0 || MIN_VALID_AREA_COVERAGE > 1)
  stop("MIN_VALID_AREA_COVERAGE must be a single number between 0 and 1.", call. = FALSE)

OUT_DIRS <- c(OUT_DIR, SWEI_DIR, SSPI_DIR, REGIONAL_DIR)
for (d in OUT_DIRS) dir.create(d, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
terra::terraOptions(progress = 1)
## ===========================================================================

## ============================== FUNCTIONS ==================================

# ---- time axis: use terra's, else read the NetCDF "time"/"valid_time" variable ----
get_dates <- function(r, file) {
  tt <- terra::time(r)
  if (length(tt) == terra::nlyr(r) && !anyNA(tt)) return(as.Date(tt))
  if (requireNamespace("ncdf4", quietly = TRUE)) {
    nc <- ncdf4::nc_open(file); on.exit(ncdf4::nc_close(nc))
    tn <- intersect(c("valid_time", "time", "date"), names(nc$dim))
    if (length(tn)) {
      tv <- ncdf4::ncvar_get(nc, tn[1])
      un <- ncdf4::ncatt_get(nc, tn[1], "units")$value
      if (is.character(un)) {
        m <- regmatches(un, regexec("^(seconds|hours|days) since (.+)$", un))[[1]]
        if (length(m) == 3) {
          d <- tryCatch(as.Date(as.POSIXct(m[3], tz = "UTC") +
                                  tv * c(seconds = 1, hours = 3600, days = 86400)[[m[2]]]),
                        error = function(e) NULL)
          if (!is.null(d) && length(d) == terra::nlyr(r)) return(d)
        }
      }
    }
  }
  warning("No usable time axis in NetCDF; assuming monthly layers starting 1950-01-01.")
  seq(as.Date("1950-01-01"), by = "month", length.out = terra::nlyr(r))
}

# ---- k-month RIGHT-aligned rolling sum over matrix columns (NA propagates) ----
roll_sum_right <- function(M, k) {
  if (k == 1L) return(M)
  n <- ncol(M); acc <- M
  for (s in seq_len(k - 1L))
    acc <- acc + cbind(matrix(NA_real_, nrow(M), s), M[, seq_len(n - s), drop = FALSE])
  acc
}

# ---- empirical Gringorten score of x against a reference sample ----
# mid-rank handles ties (zeros stay tied: no random jitter -> reproducible)
gringorten_z <- function(x, ref, eps = 1e-6) {
  ref <- sort(ref[is.finite(ref)])
  out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  if (length(ref) < 3L || !any(ok)) return(out)
  n_lt <- findInterval(x[ok], ref, left.open = TRUE)
  n_le <- findInterval(x[ok], ref)
  p    <- ((n_lt + n_le + 1) / 2 - 0.44) / (length(ref) + 0.12)
  out[ok] <- qnorm(pmin(pmax(p, eps), 1 - eps))
  out
}

# ---- SWEI for ONE pixel: v = k-month SWE sums (length = n months) ----
swei_pixel <- function(v, mon, ref_mask, min_nonzero, eps) {
  z <- rep(NA_real_, length(v))
  if (!any(is.finite(v))) return(z)
  for (m in 1:12) {
    im  <- which(mon == m)
    ref <- v[im[ref_mask[im]]]; ref <- ref[is.finite(ref)]
    if (length(ref) < 3L || mean(ref > 0) < min_nonzero) next   # not enough snow: leave NA
    z[im] <- gringorten_z(v[im], ref, eps)
  }
  z
}

# ---- SSPI for ONE pixel ----
# method codes: 1 = zero-inflated gamma | 2 = empirical (Gringorten) | 3 = undefined reference (NA)
sspi_pixel <- function(v, mon, ref_mask, min_n, min_pos, eps) {
  n <- length(v)
  z <- rep(NA_real_, n); meth <- rep(NA_integer_, n)
  if (!any(is.finite(v))) return(list(z = z, method = meth))
  for (m in 1:12) {
    im  <- which(mon == m)
    obs <- v[im]; fin <- is.finite(obs)
    if (!any(fin)) next
    ref <- obs[ref_mask[im]]; ref <- ref[is.finite(ref)]
    nr  <- length(ref)
    
    # no usable reference (too short, or constant e.g. always zero in summer)
    if (nr < 3L || !is.finite(var(ref)) || var(ref) < .Machine$double.eps) {
      meth[im[fin]] <- 3L; next
    }
    # short or almost-constant reference -> empirical
    if (nr < min_n || var(ref) < 0.01) {
      z[im] <- gringorten_z(obs, ref, eps); meth[im[fin]] <- 2L; next
    }
    
    # zero-inflated gamma
    n0 <- sum(ref <= 0)
    p0 <- n0 / (nr + 1)                    # mass for the mixture
    pz <- (n0 + 1) / (2 * (nr + 1))        # probability assigned to an observed zero
    pos <- ref[ref > 0]
    done <- FALSE
    if (length(pos) >= min_pos) {
      par <- try(lmomco::pargam(lmomco::lmoms(pos)), silent = TRUE)
      ok  <- !inherits(par, "try-error") && !is.null(par) &&
        (is.null(par$ifail) || (length(par$ifail) == 1L && par$ifail == 0))
      if (ok) {
        p  <- rep(NA_real_, length(obs))
        ip <- which(fin & obs > 0)
        if (length(ip)) {
          Fg <- try(as.numeric(lmomco::cdfgam(obs[ip], par)), silent = TRUE)
          if (!inherits(Fg, "try-error") && length(Fg) == length(ip) && all(is.finite(Fg)))
            p[ip] <- p0 + (1 - p0) * Fg
        }
        p[fin & obs <= 0] <- pz
        if (all(is.finite(p[fin]))) {
          z[im[fin]] <- qnorm(pmin(pmax(p[fin], eps), 1 - eps))
          meth[im[fin]] <- 1L; done <- TRUE
        }
      }
    }
    if (!done) {                           # gamma failed -> empirical fallback
      z[im] <- gringorten_z(obs, ref, eps); meth[im[fin]] <- 2L
    }
  }
  list(z = z, method = meth)
}

# ---- apply FUN to every row of X (serial or PSOCK-parallel) ----
map_rows <- function(X, FUN, ...) {
  rows <- lapply(seq_len(nrow(X)), function(i) X[i, ])
  if (N_CORES > 1L) {
    cl <- parallel::makeCluster(N_CORES); on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("gringorten_z", "swei_pixel", "sspi_pixel"), envir = globalenv())
    parallel::parLapply(cl, rows, FUN, ...)
  } else lapply(rows, FUN, ...)
}

## ================================ 1. LOAD ===================================
cat("\n== 1. Load ERA5-Land SWE and validate the framework grid ==\n")
if (!file.exists(SWE_FILE)) stop("SWE file not found: ", SWE_FILE, call. = FALSE)
if (!file.exists(NATIVE_PRECIP_FILE))
  stop("Native ERA5-Land precipitation grid not found: ", NATIVE_PRECIP_FILE,
       "\nRun Monthly_ERALand_download.R first.", call. = FALSE)

swe <- terra::rast(SWE_FILE)
native_precip <- terra::rast(NATIVE_PRECIP_FILE)
if (!inherits(swe, "SpatRaster") || terra::nlyr(swe) < 1L)
  stop("Could not read SWE NetCDF as a non-empty SpatRaster.", call. = FALSE)
if (!inherits(native_precip, "SpatRaster") || terra::nlyr(native_precip) < 1L)
  stop("Could not read precipitation NetCDF as a non-empty SpatRaster.", call. = FALSE)
analysis_template <- native_precip[[1L]]
if (is.na(terra::crs(analysis_template)) || !nzchar(terra::crs(analysis_template)))
  stop("Native precipitation grid has no CRS.", call. = FALSE)
if (is.na(terra::crs(swe)) || !nzchar(terra::crs(swe)))
  stop("SWE grid has no CRS.", call. = FALSE)

# Scripts 1, 3, 6 and 7 preserve the native precipitation geometry. Do not rotate,
# crop, resample, or reproject SWE here: a geometry mismatch is a data-contract
# failure and must be corrected in the upstream downloader instead.
geom_ok <- terra::compareGeom(swe[[1L]], analysis_template, stopOnError = FALSE,
                              crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)
origin_ok <- isTRUE(all.equal(terra::origin(swe[[1L]]),
                              terra::origin(analysis_template), tolerance = 1e-6))
if (!isTRUE(geom_ok) || !origin_ok) {
  stop("SWE NetCDF does not match the native ERA5-Land precipitation grid.\n",
       "Script 8 does not resample, rotate, or crop raster products.\n",
       "Check monthly_data_direct/snow_depth_water_equivalent_monthly.nc.",
       call. = FALSE)
}

# Extract and normalize monthly dates.
dates <- get_dates(swe, SWE_FILE)
dates <- as.Date(format(dates, "%Y-%m-01"))
if (length(dates) != terra::nlyr(swe) || anyNA(dates))
  stop("SWE time axis is missing or does not match the raster layer count.", call. = FALSE)
if (is.unsorted(dates)) {
  o <- order(dates)
  swe <- swe[[o]]
  dates <- dates[o]
}
ym <- as.integer(format(dates, "%Y")) * 12L + as.integer(format(dates, "%m"))
if (anyDuplicated(ym))
  stop("Duplicate months in the SWE file. Keep one layer per calendar month.", call. = FALSE)
if (length(dates) > 1L && any(diff(ym) != 1L))
  stop("SWE series is not a continuous monthly record; rolling sums would be wrong.", call. = FALSE)
if (!identical(dates, seq(ANALYSIS_START, ANALYSIS_END, by = "month"))) {
  stop("SWE time axis must cover the framework record 1950-01 through 2025-12 ",
       "(912 monthly layers); found ", format(min(dates)), " through ", format(max(dates)),
       " with ", length(dates), " layers.", call. = FALSE)
}
validate_monthly_record(dates, ANALYSIS_START, ANALYSIS_END,
                        label = "SWE monthly record", require_complete = TRUE)
terra::time(swe) <- dates
nt <- length(dates)
cat(sprintf("   Native grid: %d rows x %d columns x %d months (%s to %s)\n",
            terra::nrow(swe), terra::ncol(swe), nt,
            format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m")))

# Authoritative whole-basin geometry: dissolved union of the 7 subbasins
# (identical construction to 3SPEI_ERALand.R).
basin <- nechako_subbasin_union(load_nechako_subbasins())
if (!terra::same.crs(basin, analysis_template))
  basin <- terra::project(basin, terra::crs(analysis_template))
coverage <- terra::rasterize(basin, analysis_template, cover = TRUE)
coverage_values <- as.numeric(terra::values(coverage, mat = FALSE))
coverage_values[!is.finite(coverage_values)] <- 0
if (any(coverage_values < -1e-6 | coverage_values > 1 + 1e-6, na.rm = TRUE)) {
  stop("Basin coverage raster contains values outside [0, 1]. Refusing to clamp invalid coverage.", call. = FALSE)
}
coverage_values <- pmax(0, pmin(1, coverage_values))
cell_area_check <- as.numeric(terra::values(terra::cellSize(analysis_template, unit = "m"), mat = FALSE))
raster_basin_area_m2 <- sum(cell_area_check * coverage_values, na.rm = TRUE)
polygon_basin_area_m2 <- tryCatch(
  sum(as.numeric(terra::expanse(basin, unit = "m", transform = TRUE)), na.rm = TRUE),
  error = function(e) NA_real_
)
if (is.finite(polygon_basin_area_m2) && polygon_basin_area_m2 > 0) {
  area_rel_diff <- abs(raster_basin_area_m2 - polygon_basin_area_m2) / polygon_basin_area_m2
  cat(sprintf("   Basin area check: raster %.3f km²; polygon %.3f km²; difference %.2f%%\n",
              raster_basin_area_m2 / 1e6, polygon_basin_area_m2 / 1e6, 100 * area_rel_diff))
  if (area_rel_diff > 0.05) {
    stop(sprintf("Rasterized basin area differs from polygon area by %.2f%% (>5%%). Check native-grid coverage.",
                 100 * area_rel_diff), call. = FALSE)
  }
} else {
  warning("Could not calculate polygon area for the basin coverage sanity check.")
}
basin_mask <- terra::setValues(analysis_template, as.numeric(coverage_values > 0))
swe <- terra::mask(swe, basin_mask, maskvalues = 0, updatevalue = NA)
if (sum(coverage_values > 0) == 0L)
  stop("Authoritative Nechako basin mask contains zero cells.", call. = FALSE)
cat(sprintf("   Basin mask: %d native cells; fractional coverage retained for means.\n",
            sum(coverage_values > 0)))

# ERA5-Land snow water equivalent is stored in metres of water equivalent.
# Convert explicitly to mm before calculating rolling sums and distributions.
if (SWE_UNITS == "m") {
  swe <- swe * 1000
  cat("   Converted SWE m -> mm using the configured ERA5-Land units.\n")
}
X <- terra::values(swe, mat = TRUE)                            # cells x months
X[is.finite(X) & X < 0] <- 0
valid_rows <- which(rowSums(is.finite(X)) > 0)
if (!length(valid_rows)) stop("No valid SWE pixels after masking.", call. = FALSE)
cat(sprintf("   %d valid pixels\n", length(valid_rows)))

mon      <- as.integer(format(dates, "%m"))
# Whole-record calibration mask (1950-01 .. 2025-12): every month is in the
# reference sample for both SWEI and SSPI.
ref_mask <- dates >= CALIB_START & dates <= CALIB_END
if (!all(ref_mask))
  stop("Calibration mask does not cover the whole 1950-2025 record.", call. = FALSE)
cat(sprintf("   Calibration period (SWEI and SSPI): %s to %s (%d months)\n",
            format(CALIB_START, "%Y-%m"), format(CALIB_END, "%Y-%m"), sum(ref_mask)))
Xv <- X[valid_rows, , drop = FALSE]

## ================================ 2. SWEI ===================================
cat(sprintf("\n== 2. SWEI-%d (Gringorten, %s) ==\n", SWEI_SCALE, SWEI_REFERENCE))
swei_ref <- ref_mask
res <- map_rows(roll_sum_right(Xv, SWEI_SCALE), swei_pixel,
                mon = mon, ref_mask = swei_ref, min_nonzero = MIN_NONZERO_FRACTION, eps = PROB_EPS)
SWEI <- matrix(NA_real_, nrow(X), nt); SWEI[valid_rows, ] <- do.call(rbind, res)

## ================================ 3. SSPI ===================================
cat(sprintf("\n== 3. SSPI-%d (zero-inflated gamma, reference %s to %s) ==\n",
            SSPI_SCALE, format(CALIB_START), format(CALIB_END)))
Xs        <- roll_sum_right(Xv, SSPI_SCALE)
# Do not apply a whole-year zero fraction filter here: seasonal snow pixels
# commonly have zero SWE for much of the year but valid winter distributions.
# sspi_pixel() decides independently for each calendar month whether its
# reference distribution is estimable and returns NA/method 3 when undefined.
res  <- map_rows(Xs, sspi_pixel, mon = mon, ref_mask = ref_mask,
                 min_n = MIN_GAMMA_N, min_pos = MIN_POSITIVE_N, eps = PROB_EPS)
SSPI <- METH <- matrix(NA_real_, nrow(X), nt)
SSPI[valid_rows, ] <- do.call(rbind, lapply(res, `[[`, "z"))
METH[valid_rows, ] <- do.call(rbind, lapply(res, `[[`, "method"))

## ================================ 4. OUTPUT =================================
cat("\n== 4. Write framework-compatible NetCDF and regional-average outputs ==\n")

month_names <- tolower(month.abb)

write_full_stack <- function(mat, path, varname, longname, unit = "standardized_index") {
  r <- terra::rast(analysis_template, nlyrs = nt)
  terra::values(r) <- mat
  terra::time(r) <- dates
  names(r) <- format(dates, "%Y-%m")
  terra::writeCDF(r, path, varname = varname, longname = longname,
                  unit = unit, missval = -9999, overwrite = TRUE)
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0)
    stop("NetCDF write failed or produced an empty file: ", path, call. = FALSE)
  invisible(path)
}

# Preserve convenient full-record stacks, while also writing the standard
# 12-file seasonal layout consumed by find_seasonal_nc_files() in Script 6.
write_seasonal_files <- function(mat, index_name, index_dir, varname, longname) {
  for (m in 1:12) {
    idx <- which(as.integer(format(dates, "%m")) == m)
    if (!length(idx)) stop("No observations for calendar month ", m, call. = FALSE)
    r <- terra::rast(analysis_template, nlyrs = length(idx))
    terra::values(r) <- mat[, idx, drop = FALSE]
    terra::time(r) <- dates[idx]
    names(r) <- format(dates[idx], "%Y-%m")
    path <- file.path(index_dir,
                      sprintf("%s_month%02d_%s.nc", index_name, m, month_names[m]))
    terra::writeCDF(r, path, varname = varname, longname = longname,
                    unit = "standardized_index", missval = -9999, overwrite = TRUE)
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0)
      stop("Seasonal NetCDF write failed or produced an empty file: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

swei_name <- sprintf("swei_%02d", SWEI_SCALE)
sspi_name <- sprintf("sspi_%02d", SSPI_SCALE)
write_full_stack(SWEI, file.path(OUT_DIR, paste0(swei_name, ".nc")),
                 "swei", sprintf("Standardized SWE Index (%d-month)", SWEI_SCALE))
write_full_stack(SSPI, file.path(OUT_DIR, paste0(sspi_name, ".nc")),
                 "sspi", sprintf("Standardized SnowPack Index (%d-month)", SSPI_SCALE))
write_full_stack(METH, file.path(OUT_DIR, paste0(sspi_name, "_method.nc")),
                 "sspi_method", "SSPI method: 1 = zero-inflated gamma, 2 = empirical, 3 = undefined reference", "code")
write_seasonal_files(SWEI, swei_name, SWEI_DIR, "swei",
                     sprintf("Standardized SWE Index (%d-month)", SWEI_SCALE))
write_seasonal_files(SSPI, sspi_name, SSPI_DIR, "sspi",
                     sprintf("Standardized SnowPack Index (%d-month)", SSPI_SCALE))

# Area-weighted regional means use the shared native-grid fractional coverage
# and cell-area weights. Regions include the whole basin and every supplied
# subbasin; no unweighted cell mean is used for framework-facing averages.
regions <- load_nechako_regions(template = analysis_template)
validate_nechako_regions(regions)
region_weights <- region_fractional_weights(analysis_template, regions)

# Renormalize weights over finite cells only when at least `min_coverage` of the
# region's total area has a valid index value. This avoids summer means based on
# a few surviving high-elevation cells being presented as whole-region means.
# Two cases are reported (see header): the user-set MIN_VALID_AREA_COVERAGE
# ("custom", defined in CONFIGURATION) and 0 ("zero": any region-month with at
# least one finite cell gets a mean).
# regional_valid_coverage() and regional_weighted_mean() come from
# DROUGHT_ANALYSIS_utils.R so this script uses exactly the same estimator as
# every other script. A month with no finite cell is always NA.
regional_weighted_mean_covered <- function(mat, weights, min_coverage) {
  regional_weighted_mean(mat, weights, min_coverage = min_coverage)
}

# The two cases: name, threshold, and file-name suffix. The "custom" case keeps
# the original file names (suffix ""), "zero" adds "_mincov0".
COVERAGE_CASES <- list(
  custom = list(min_coverage = MIN_VALID_AREA_COVERAGE, suffix = ""),
  zero   = list(min_coverage = COVERAGE_ZERO,           suffix = "_mincov0")
)
cat(sprintf("   Valid-area coverage cases: custom = %.2f | zero = %.2f\n",
            COVERAGE_CASES$custom$min_coverage, COVERAGE_CASES$zero$min_coverage))
if (isTRUE(all.equal(COVERAGE_CASES$custom$min_coverage, COVERAGE_CASES$zero$min_coverage)))
  message("   NOTE: MIN_VALID_AREA_COVERAGE equals 0, so both cases are identical.")

# Year x month wide table (Jan..Dec) from a monthly vector aligned with `dates`.
by_month_wide <- function(v) {
  yrs <- as.integer(format(dates, "%Y"))
  mos <- as.integer(format(dates, "%m"))
  years <- sort(unique(yrs))
  wide <- data.frame(Year = years)
  for (m in 1:12) {
    vv <- rep(NA_real_, length(years))
    hit <- which(mos == m)
    vv[match(yrs[hit], years)] <- v[hit]
    wide[[month.abb[m]]] <- vv
  }
  wide
}

# Whole-basin Year x month CSV for one coverage case. The valid-area-coverage CSV
# does not depend on the threshold, so it is written once (primary case only).
write_basin_average <- function(mat, index_name, index_dir, min_coverage,
                                suffix = "", write_coverage = (suffix == "")) {
  rw <- region_weights[[NECHAKO_REGION_BASIN_ID]]
  value <- regional_weighted_mean_covered(mat, rw$weights, min_coverage)
  coverage <- regional_valid_coverage(mat, rw$weights)
  utils::write.csv(by_month_wide(value),
                   file.path(index_dir, sprintf("%s_basin_averaged_by_month%s.csv", index_name, suffix)),
                   row.names = FALSE, na = "")
  if (write_coverage)
    utils::write.csv(by_month_wide(coverage),
                     file.path(index_dir, sprintf("%s_basin_valid_area_coverage_by_month.csv", index_name)),
                     row.names = FALSE, na = "")
  invisible(data.frame(date = dates, value = value, valid_area_coverage_fraction = coverage))
}

# Basin + every subbasin, monthly, for one coverage case (original long format).
regional_average_table <- function(mat, index_name, scale, min_coverage) {
  data.table::rbindlist(lapply(names(region_weights), function(id) {
    rw <- region_weights[[id]]
    data.table::data.table(
      index_type = index_name,
      scale = as.integer(scale),
      region_id = id,
      region_name = rw$region_name,
      date = dates,
      value = regional_weighted_mean_covered(mat, rw$weights, min_coverage),
      valid_area_coverage_fraction = regional_valid_coverage(mat, rw$weights),
      coverage_area_m2 = rw$area_m2,
      n_cells = rw$n_cells,
      n_partial_cells = rw$n_partial
    )
  }), use.names = TRUE)
}

# Writes: (a) the original-format regional CSV for each case, and (b) one
# comparison CSV with both cases side by side. Returns a list of the two tables.
write_regional_average <- function(mat, index_name, scale) {
  tabs <- lapply(COVERAGE_CASES, function(cs)
    regional_average_table(mat, index_name, scale, cs$min_coverage))
  for (cn in names(COVERAGE_CASES))
    data.table::fwrite(tabs[[cn]],
                       file.path(REGIONAL_DIR, sprintf("%s_%02d_regional_averaged_by_month%s.csv",
                                                       index_name, scale, COVERAGE_CASES[[cn]]$suffix)),
                       na = "NA")
  cmp <- data.table::copy(tabs$custom)[, c("value", "coverage_area_m2", "n_cells", "n_partial_cells") := NULL]
  cmp[, min_coverage_custom := COVERAGE_CASES$custom$min_coverage]
  cmp[, value_custom := tabs$custom$value]
  cmp[, value_cov0   := tabs$zero$value]
  cmp[, reported_custom := is.finite(value_custom)]
  cmp[, reported_cov0   := is.finite(value_cov0)]
  data.table::fwrite(cmp,
                     file.path(REGIONAL_DIR, sprintf("%s_%02d_regional_averaged_by_month_coverage_comparison.csv",
                                                     index_name, scale)),
                     na = "NA")
  invisible(tabs)
}

basin_swei <- write_basin_average(SWEI, swei_name, SWEI_DIR, COVERAGE_CASES$custom$min_coverage, "")
basin_sspi <- write_basin_average(SSPI, sspi_name, SSPI_DIR, COVERAGE_CASES$custom$min_coverage, "")
basin_swei_cov0 <- write_basin_average(SWEI, swei_name, SWEI_DIR, COVERAGE_CASES$zero$min_coverage, "_mincov0")
basin_sspi_cov0 <- write_basin_average(SSPI, sspi_name, SSPI_DIR, COVERAGE_CASES$zero$min_coverage, "_mincov0")
regional_swei <- write_regional_average(SWEI, "swei", SWEI_SCALE)
regional_sspi <- write_regional_average(SSPI, "sspi", SSPI_SCALE)

# ---- Average values per region x coverage case x calendar month ----
# mean_value             : mean of the reported monthly area-weighted index values
# n_months_reported      : months in which a regional mean was reported (finite)
# mean_valid_area_coverage: mean fraction of the region's area with a valid index
# period "all" pools every month of the record; Jan..Dec are calendar months.
regional_mean_summary <- data.table::rbindlist(lapply(
  list(regional_swei, regional_sspi), function(tabs) {
    data.table::rbindlist(lapply(names(COVERAGE_CASES), function(cn) {
      tab <- data.table::copy(tabs[[cn]])
      tab[, mo := as.integer(format(date, "%m"))]
      summ <- function(d) list(
        n_months_total = nrow(d),
        n_months_reported = sum(is.finite(d$value)),
        mean_value = if (any(is.finite(d$value))) mean(d$value[is.finite(d$value)]) else NA_real_,
        mean_valid_area_coverage = mean(d$valid_area_coverage_fraction, na.rm = TRUE))
      by_mo  <- tab[, c(list(period = month.abb[mo[1]]), summ(.SD)),
                    by = .(index_type, scale, region_id, region_name, mo),
                    .SDcols = c("value", "valid_area_coverage_fraction")]
      by_all <- tab[, c(list(period = "all", mo = 0L), summ(.SD)),
                    by = .(index_type, scale, region_id, region_name),
                    .SDcols = c("value", "valid_area_coverage_fraction")]
      out <- data.table::rbindlist(list(by_all, by_mo), use.names = TRUE)
      out[, `:=`(coverage_case = cn,
                 min_valid_area_coverage = COVERAGE_CASES[[cn]]$min_coverage)]
      out
    }), use.names = TRUE)
  }), use.names = TRUE)
data.table::setorder(regional_mean_summary, index_type, region_id, mo, coverage_case)
regional_mean_summary[, mo := NULL]
data.table::setcolorder(regional_mean_summary,
                        c("index_type", "scale", "region_id", "region_name", "coverage_case",
                          "min_valid_area_coverage", "period", "n_months_total",
                          "n_months_reported", "mean_value", "mean_valid_area_coverage"))
data.table::fwrite(regional_mean_summary,
                   file.path(REGIONAL_DIR, "snow_drought_regional_mean_summary.csv"), na = "NA")

# Keep the earlier simple two-index basin chronology as a convenience product.
basin_df <- data.frame(
  date = dates,
  SWEI = basin_swei$value,
  SSPI = basin_sspi$value,
  n_valid_swei = colSums(is.finite(SWEI[valid_rows, , drop = FALSE])),
  n_valid_sspi = colSums(is.finite(SSPI[valid_rows, , drop = FALSE]))
)
utils::write.csv(basin_df, file.path(OUT_DIR, "basin_mean_indices.csv"),
                 row.names = FALSE, na = "")

png(file.path(OUT_DIR, "basin_mean_indices.png"), width = 1600, height = 900, res = 150)
op <- par(mfrow = c(2, 1), mar = c(3, 4, 2.5, 1))
for (nm in c("SWEI", "SSPI")) {
  z <- basin_df[[nm]]
  plot(dates, z, type = "h", ylim = c(-3.5, 3.5), xlab = "", ylab = "z-score",
       col = ifelse(is.na(z) | z >= 0, "steelblue", "firebrick"),
       main = paste("Area-weighted Nechako basin mean", nm))
  abline(h = c(-2, -1, 0, 1, 2), lty = c(3, 3, 1, 3, 3), col = "grey50")
}
par(op)
invisible(dev.off())

# Machine-readable manifest lets downstream scripts verify exact products.
manifest <- data.table::data.table(
  index = c("SWEI", "SSPI"),
  scale = c(SWEI_SCALE, SSPI_SCALE),
  output_directory = c(normalizePath(SWEI_DIR, winslash = "/", mustWork = FALSE),
                       normalizePath(SSPI_DIR, winslash = "/", mustWork = FALSE)),
  seasonal_netcdf_files = 12L,
  first_date = min(dates),
  last_date = max(dates),
  n_months = nt,
  swei_reference = "1950-2025 whole record",
  sspi_reference = "1950-2025 whole record",
  regional_min_valid_area_coverage = COVERAGE_CASES$custom$min_coverage,
  regional_min_valid_area_coverage_zero_case = COVERAGE_CASES$zero$min_coverage,
  era5land_swe_caveat = "Pre-1979 SWE is less observationally constrained; no glacier mask is applied.",
  spatial_grid = "native ERA5-Land precipitation grid"
)
data.table::fwrite(manifest, file.path(OUT_DIR, "snow_drought_output_manifest.csv"))

## ================================ 5. SUMMARY ================================
share <- function(M, cond) 100 * mean(cond(M[valid_rows, , drop = FALSE]), na.rm = TRUE)
cat(sprintf("\n   %s: %.1f %% of pixel-months scored | %.1f %% below -1 | %.1f %% below -2\n",
            swei_name, share(SWEI, is.finite), share(SWEI, function(m) ifelse(is.finite(m), m < -1, NA)),
            share(SWEI, function(m) ifelse(is.finite(m), m < -2, NA))))
cat(sprintf("   %s: %.1f %% of pixel-months scored | %.1f %% below -1 | %.1f %% below -2\n",
            sspi_name, share(SSPI, is.finite), share(SSPI, function(m) ifelse(is.finite(m), m < -1, NA)),
            share(SSPI, function(m) ifelse(is.finite(m), m < -2, NA))))
mt <- table(METH[valid_rows, ][is.finite(METH[valid_rows, ])])
cat("   SSPI methods (pixel-months): ",
    paste(sprintf("%s=%d", c("1=gamma", "2=empirical", "3=undefined")[as.integer(names(mt))], mt), collapse = ", "), "\n")
required_outputs <- c(
  file.path(OUT_DIR, paste0(swei_name, ".nc")),
  file.path(OUT_DIR, paste0(sspi_name, ".nc")),
  file.path(OUT_DIR, "basin_mean_indices.csv"),
  file.path(OUT_DIR, "snow_drought_output_manifest.csv"),
  unlist(lapply(1:12, function(m) file.path(SWEI_DIR,
                                            sprintf("%s_month%02d_%s.nc", swei_name, m, month_names[m]))), use.names = FALSE),
  unlist(lapply(1:12, function(m) file.path(SSPI_DIR,
                                            sprintf("%s_month%02d_%s.nc", sspi_name, m, month_names[m]))), use.names = FALSE),
  file.path(SWEI_DIR, sprintf("%s_basin_averaged_by_month.csv", swei_name)),
  file.path(SSPI_DIR, sprintf("%s_basin_averaged_by_month.csv", sspi_name)),
  file.path(SWEI_DIR, sprintf("%s_basin_valid_area_coverage_by_month.csv", swei_name)),
  file.path(SSPI_DIR, sprintf("%s_basin_valid_area_coverage_by_month.csv", sspi_name)),
  # two-coverage-case products
  file.path(SWEI_DIR, sprintf("%s_basin_averaged_by_month_mincov0.csv", swei_name)),
  file.path(SSPI_DIR, sprintf("%s_basin_averaged_by_month_mincov0.csv", sspi_name)),
  unlist(lapply(c(sprintf("swei_%02d", SWEI_SCALE), sprintf("sspi_%02d", SSPI_SCALE)), function(nm)
    file.path(REGIONAL_DIR, sprintf("%s_regional_averaged_by_month%s.csv", nm,
                                    c("", "_mincov0", "_coverage_comparison")))), use.names = FALSE),
  file.path(REGIONAL_DIR, "snow_drought_regional_mean_summary.csv")
)
missing_outputs <- required_outputs[!file.exists(required_outputs) |
                                      is.na(file.info(required_outputs)$size) |
                                      file.info(required_outputs)$size <= 0]
if (length(missing_outputs))
  stop("Final output validation failed; missing/empty product(s):\n",
       paste(missing_outputs, collapse = "\n"), call. = FALSE)

cat("\n   Area-weighted mean index, all months (basin and subbasins):\n")
tab_all <- data.table::dcast(regional_mean_summary[period == "all"],
                             index_type + region_id + region_name ~ coverage_case,
                             value.var = c("mean_value", "n_months_reported"))
print(tab_all, digits = 3, row.names = FALSE)
cat(sprintf("   (custom = MIN_VALID_AREA_COVERAGE %.2f | zero = 0)\n", COVERAGE_CASES$custom$min_coverage))

cat("\nFramework output validation passed.\n")
cat(sprintf("  SWEI seasonal files : %s\n", normalizePath(SWEI_DIR, winslash = "/", mustWork = FALSE)))
cat(sprintf("  SSPI seasonal files : %s\n", normalizePath(SSPI_DIR, winslash = "/", mustWork = FALSE)))
cat(sprintf("  Regional averages   : %s\n", normalizePath(REGIONAL_DIR, winslash = "/", mustWork = FALSE)))
cat(sprintf("  Mean summary        : %s\n", file.path(normalizePath(REGIONAL_DIR, winslash = "/", mustWork = FALSE),
                                                      "snow_drought_regional_mean_summary.csv")))
cat(sprintf("Done in %.1f min.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))