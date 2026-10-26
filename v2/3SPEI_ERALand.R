##############################################
# SPEI CALCULATION FOR NECHAKO RIVER BASIN
#
# Current manuscript framing:
#   - SPEI_PM  = observed Penman???Monteith SPEI
#   - SPEI_ERA5 = observed SPEI using ERA5-Land native PEV
#   - SPEI_PET = SPEI_PM - SPEI_ERA5 is a PET-formulation sensitivity residual
#     and is not interpreted as an additive physical decomposition.
#
# This script computes the two observed SPEI branches only.
# The detrended/counterfactual PET branch has been removed.
##############################################


# ---- Libraries ----
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

library(terra)
library(lubridate)
library(zoo)
library(writexl)
library(lmomco)
library(parallel)

# ---- Paths ----
out_dir     <- "spei_results_seasonal"
out_dir_era5 <- "spei_results_seasonal_era5land"
source(UTILS_FILE)

if (!dir.exists(out_dir))     dir.create(out_dir,     recursive = TRUE)
if (!dir.exists(out_dir_era5)) dir.create(out_dir_era5, recursive = TRUE)

terraOptions(progress = 0)

# ---- Load authoritative whole-basin boundary ----
# The whole basin is the dissolved union of the 7 supplied subbasins.
basin <- nechako_subbasin_union(load_nechako_subbasins())
cat(sprintf("Basin boundary loaded from 7 subbasins (CRS: %s)\n", terra::crs(basin)))
##############################################
# ---- Load Data ----
cat("\n===== LOADING INPUT DATA =====\n")
precip <- rast("monthly_data_direct/total_precipitation_monthly.nc")
# ---- PRECIPITATION UNITS ----
# Shared contract (DROUGHT_ANALYSIS_utils.R, resolve_monthly_depth_unit()):
# the NetCDF units attribute and the data magnitude must agree, otherwise the
# run STOPS. Each input is converted once to mm/day here; the water balance is
# then multiplied by actual calendar days exactly once (below), so a product
# delivered as mm/month can no longer be double-counted.
precip_dates_raw <- as.Date(time(precip))
if (is.null(precip_dates_raw) || all(is.na(precip_dates_raw))) {
  precip_dates_raw <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(precip))
}
precip_unit_info <- resolve_monthly_depth_unit(
  precip, kind = "precipitation",
  file = "monthly_data_direct/total_precipitation_monthly.nc",
  label = "precipitation (SPEI input)")
precip <- convert_monthly_depth(precip, from = precip_unit_info$unit,
                                to = "mm/day", dates = precip_dates_raw)

pet    <- rast("monthly_data_direct/potential_evapotranspiration_monthly.nc")
# ---- PET (FAO-56 PM) UNITS: same shared contract; converted once to mm/day ----
pet_dates_raw <- as.Date(time(pet))
if (is.null(pet_dates_raw) || all(is.na(pet_dates_raw))) {
  pet_dates_raw <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(pet))
}
pet_unit_info <- resolve_monthly_depth_unit(
  pet, kind = "pet",
  file = "monthly_data_direct/potential_evapotranspiration_monthly.nc",
  label = "PET FAO-56 PM (SPEI input)")
pet <- convert_monthly_depth(pet, from = pet_unit_info$unit,
                             to = "mm/day", dates = pet_dates_raw)

# ERA5-Land native potential evaporation (PEV), already converted to mm/day by 2b.
pet_era5_file <- "monthly_data_direct/potential_evaporation_era5land_monthly.nc"
if (!file.exists(pet_era5_file))
  stop(paste("ERA5-Land PEV file not found:", pet_era5_file,
             "\nRun 2b_PET_ERALand.R first."))
pet_era5 <- rast(pet_era5_file)
cat("ERA5-Land native PEV loaded\n")
pev_dates_raw <- as.Date(time(pet_era5))
if (is.null(pev_dates_raw) || all(is.na(pev_dates_raw))) {
  pev_dates_raw <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(pet_era5))
}
pev_unit_info <- resolve_monthly_depth_unit(
  pet_era5, kind = "pev", file = pet_era5_file,
  label = "ERA5-Land PEV (SPEI input)")
pet_era5 <- convert_monthly_depth(pet_era5, from = pev_unit_info$unit,
                                  to = "mm/day", dates = pev_dates_raw)

##############################################
# NATIVE ERA5-LAND GRID CONTRACT
#
# Framework rule:
#   Script 3 uses the actual spatial grid carried by the native ERA5-Land
#   precipitation NetCDF. No canonical/projected/resampled raster grid is
#   introduced and no raster product is silently interpolated onto another
#   grid.
#
# The precipitation raster is therefore the native spatial template for the
# SPEI calculation. The basin vector may be transformed to the raster CRS,
# but the ERA5-Land raster geometry itself is never changed.
##############################################
analysis_template <- precip[[1]]

assert_native_grid <- function(x, template, label) {
  if (!inherits(x, "SpatRaster")) {
    stop(label, " is not a SpatRaster.")
  }
  
  # NOTE: terra >=1.8 dropped the `orig` argument from compareGeom(). Passing
  # orig = TRUE previously fell into `...` (treated as "additional SpatRaster
  # to compare") and raised:
  #   Error: [compareGeom] all additional arguments must be a SpatRaster
  # Origin equality is checked separately below with terra::origin().
  geom_ok <- terra::compareGeom(
    x[[1]],
    template,
    stopOnError = FALSE,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE
  )
  
  orig_ok <- isTRUE(all.equal(
    terra::origin(x[[1]]),
    terra::origin(template),
    tolerance = 1e-6
  ))
  
  if (!geom_ok || !orig_ok) {
    stop(
      label,
      " does not match the native ERA5-Land precipitation grid.\n",
      "Script 3 does not reproject or resample raster inputs. Check the upstream NetCDF producer."
    )
  }
  
  invisible(TRUE)
}

assert_native_grid(precip, analysis_template, "Precipitation")
assert_native_grid(pet, analysis_template, "FAO-56 PM PET")
assert_native_grid(pet_era5, analysis_template, "ERA5-Land native PEV")

if (is.na(terra::crs(analysis_template)) ||
    !nzchar(terra::crs(analysis_template))) {
  stop("Native ERA5-Land precipitation grid has no CRS.")
}

if (!terra::same.crs(basin, analysis_template)) {
  basin <- terra::project(basin, terra::crs(analysis_template))
}

cat(sprintf(
  "Native ERA5-Land grid: %d x %d = %d cells | res %.10f x %.10f | CRS=%s | Time steps: %d\n",
  terra::ncol(analysis_template), terra::nrow(analysis_template),
  terra::ncell(analysis_template),
  terra::res(analysis_template)[1], terra::res(analysis_template)[2],
  terra::crs(analysis_template), terra::nlyr(precip)
))

# ------------------------------------------------------------------------------
# Build the basin mask directly on the native ERA5-Land precipitation grid.
# No stored raster mask from another grid is used. `cover = TRUE` preserves
# partial-cell intersection information; cells with any basin coverage are
# retained, matching the existing framework's basin inclusion rule.
# ------------------------------------------------------------------------------
basin_coverage <- terra::rasterize(
  basin,
  analysis_template,
  cover = TRUE
)

basin_coverage_values <- as.numeric(
  terra::values(basin_coverage, mat = FALSE)
)
basin_coverage_values[!is.finite(basin_coverage_values)] <- 0
basin_coverage_values <- pmax(0, pmin(1, basin_coverage_values))

basin_mask_raster <- terra::setValues(
  analysis_template,
  as.numeric(basin_coverage_values > 0)
)

basin_mask_geom_ok <- terra::compareGeom(
  basin_mask_raster,
  analysis_template,
  stopOnError = FALSE,
  crs = TRUE,
  ext = TRUE,
  rowcol = TRUE,
  res = TRUE
)

basin_mask_orig_ok <- isTRUE(all.equal(
  terra::origin(basin_mask_raster),
  terra::origin(analysis_template),
  tolerance = 1e-6
))

if (!basin_mask_geom_ok || !basin_mask_orig_ok) {
  stop("Native ERA5-Land basin mask geometry does not match the precipitation grid.")
}

# Basin masking is performed directly on the native ERA5-Land grid.
cat("\n===== MASKING INPUTS TO NATIVE ERA5-LAND BASIN GRID =====\n")
precip <- terra::mask(precip, basin_mask_raster, maskvalues = 0, updatevalue = NA)
pet <- terra::mask(pet, basin_mask_raster, maskvalues = 0, updatevalue = NA)
pet_era5 <- terra::mask(pet_era5, basin_mask_raster, maskvalues = 0, updatevalue = NA)
basin_mask <- as.logical(terra::values(basin_mask_raster, mat = FALSE))
basin_pixels <- sum(basin_mask, na.rm = TRUE)
total_pixels <- terra::ncell(analysis_template)

cat(sprintf(
  "Basin masking complete: %d pixels (%.1f%% of native raster)\n",
  basin_pixels,
  100 * basin_pixels / total_pixels
))

if (basin_pixels == 0L) {
  stop("CRITICAL: Native ERA5-Land basin mask produced zero valid pixels.")
}

# ------------------------------------------------------------------------------
# Fractional-coverage basin weights (shared framework rule, see
# DROUGHT_ANALYSIS_utils.R): w = cell area x fractional coverage of the
# 7-subbasin union; basin means are reported only when the finite cells hold
# >= BASIN_MIN_VALID_AREA_COVERAGE (0.20) of the total basin weight.
# ------------------------------------------------------------------------------
basin_bw <- basin_fractional_weights(analysis_template, basin)
if (length(basin_bw$mask) != length(basin_mask) || any(basin_bw$mask != basin_mask)) {
  stop("Fractional-weight basin mask differs from the SPEI basin mask.")
}
cat(sprintf(
  "Basin weights: %d cells (%d partial) | weighted area %.1f km2 | min valid coverage %.2f\n",
  basin_bw$n_cells, basin_bw$n_partial, basin_bw$area_m2 / 1e6,
  BASIN_MIN_VALID_AREA_COVERAGE
))

# Dates
# ==============================================================================
# Time axis
# Normalize all monthly timestamps to the first day of each month.
# ERA5-Land monthly NetCDFs may decode to the 2nd day depending on the calendar.
# ==============================================================================
dates <- as.Date(time(precip))

if (is.null(dates) || all(is.na(dates)) || length(dates) != nlyr(precip)) {
  
  cat("Time extraction failed - reconstructing from 1950-01...\n")
  
  dates <- seq(
    as.Date("1950-01-01"),
    by = "month",
    length.out = nlyr(precip)
  )
  
} else {
  
  # Force all monthly dates to YYYY-MM-01
  dates <- as.Date(format(dates, "%Y-%m-01"))
  
}

terra::time(precip)  <- dates
terra::time(pet)     <- dates
terra::time(pet_era5) <- dates
first_of_month    <- as.Date(format(dates, "%Y-%m-01"))
first_next_month  <- seq(first_of_month[1], by = "month", length.out = length(dates) + 1)[-1]
days_in_month     <- as.integer(first_next_month - first_of_month)
cat(sprintf("Time period: %s to %s (%d months)\n",
            min(dates), max(dates), length(dates)))

# Basin mask is already applied above. Keep the full projected raster
# geometry unchanged; outside-basin cells remain NA.
cat(sprintf(
  "Final SPEI grid: %d x %d = %d cells; basin cells = %d\n",
  nrow(precip),
  ncol(precip),
  ncell(precip),
  basin_pixels
))

# Days-in-month vector
month_nums    <- as.integer(format(dates, "%m"))
year_nums     <- as.integer(format(dates, "%Y"))
days_in_month <- as.integer(first_next_month - first_of_month)

month_numbers <- month_nums   # alias used inside variance_aware_spei
month_names   <- c("Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec")

# ==============================================================================
#   WATER BALANCE ??? PENMAN-MONTEITH  (WB_pm = P - PET_PM)
# ==============================================================================
cat("\n===== WATER BALANCE (PM) =====\n")
wb <- precip - pet
cat("Converting WB_pm from mm/day to mm/month...\n")
for (i in seq_len(nlyr(wb))) wb[[i]] <- wb[[i]] * days_in_month[i]
terra::time(wb) <- dates
cat(sprintf("WB range: %d-%d days per month\n", min(days_in_month), max(days_in_month)))

wb_mat <- values(wb, mat = TRUE)
cat(sprintf(
  "SPEI processing matrix: %d raster cells x %d months; basin cells = %d\n",
  nrow(wb_mat),
  ncol(wb_mat),
  basin_pixels
))

# Zero diagnostics for WB_pm
cat("\n===== WB_PM ZERO DIAGNOSTICS =====\n")
wb_basin <- wb_mat[basin_mask, , drop = FALSE]
total_v  <- length(wb_basin) - sum(is.na(wb_basin))
zero_c   <- sum(wb_basin == 0, na.rm = TRUE)
zero_pct <- 100 * zero_c / total_v
cat(sprintf("Total values: %d | Zeros: %d (%.4f%%)\n", total_v, zero_c, zero_pct))
if (zero_pct > 1) cat("WARNING: High proportion of zeros detected!\n")

cat("\nZero counts by month (WB_pm):\n")
for (m in 1:12) {
  idx <- which(month_numbers == m)
  if (length(idx) == 0) next
  mv  <- wb_basin[, idx, drop = FALSE]
  tot <- length(mv) - sum(is.na(mv))
  zer <- sum(mv == 0, na.rm = TRUE)
  cat(sprintf("  %s: %d / %d (%.2f%%)\n", month_names[m], zer, tot, 100 * zer / tot))
}
cat("=========================================\n\n")

# ==============================================================================
#   WATER BALANCE - ERA5-LAND NATIVE PEV  (WB_era5 = P - PET_ERA5)
# ==============================================================================
cat("\n===== WATER BALANCE (ERA5-Land PEV) =====\n")
wb_era5 <- precip - pet_era5
cat("Converting WB_era5 from mm/day to mm/month...\n")
for (i in seq_len(nlyr(wb_era5))) wb_era5[[i]] <- wb_era5[[i]] * days_in_month[i]
terra::time(wb_era5) <- dates

wb_era5_mat   <- values(wb_era5, mat = TRUE)
wb_era5_basin <- wb_era5_mat[basin_mask, , drop = FALSE]
cat(sprintf("Basin pixels for WB_era5: %d\n", nrow(wb_era5_basin)))

total_vt  <- length(wb_era5_basin) - sum(is.na(wb_era5_basin))
zero_ct   <- sum(wb_era5_basin == 0, na.rm = TRUE)
cat(sprintf("WB_era5 zeros: %d / %d (%.4f%%)\n", zero_ct, total_vt,
            100 * zero_ct / total_vt))
if (100 * zero_ct / total_vt > 1) cat("WARNING: High zeros in WB_era5!\n")

# NOTE: WB_era5 = P - PET_ERA5, where PET_ERA5 is the ERA5-Land native
# potential evaporation (PEV). PEV is a surface-energy-balance product, not a
# temperature-based PET: it is generally non-zero in winter and is not set to
# zero below 0 C. In this basin it is several times larger than FAO-56 PM PET
# (see the 2b_PET_ERALand.R diagnostics), so WB_era5 is dominated by PEV.
cat("\nZero counts by month (WB_era5):\n")
for (m in 1:12) {
  idx <- which(month_numbers == m)
  if (length(idx) == 0) next
  mv  <- wb_era5_basin[, idx, drop = FALSE]
  tot <- length(mv) - sum(is.na(mv))
  zer <- sum(mv == 0, na.rm = TRUE)
  cat(sprintf("  %s: %d / %d (%.2f%%)\n", month_names[m], zer, tot, 100 * zer / tot))
}
cat("=========================================\n\n")

# ==============================================================================
#   EXPORT: WB_basin_average_monthly.csv
#
#   Diagnostic export for external companion analyses (e.g. a PET-bias /
#   non-stationarity check); no script in this framework requires it.
#   Written to spei_results_seasonal/ (= out_dir).
#
#   Columns:
#     date           character  "YYYY-MM-DD"  (first of each month)
#     year           integer
#     month          integer    1???12
#     days_in_month  integer    28/29/30/31  (actual, leap-year aware)
#     wb_pm_mm_month numeric    basin-mean WB_PM  = P ??? PET_PM  [mm/month]
#     wb_era5_mm_month numeric   basin-mean WB_ERA5 = P ??? PET_ERA5 [mm/month]
#
#     wb_pm_valid_area_coverage / wb_era5_valid_area_coverage  numeric  0-1
#
#   Basin mean: fractional-coverage area-weighted mean (cell area x coverage
#   of the 7-subbasin union), reported only when finite cells hold >= 0.20 of
#   the total basin weight (BASIN_MIN_VALID_AREA_COVERAGE). This is the same
#   estimator used for every basin series in the framework.
#
#   NOTE: Both wb_basin (PM) and wb_era5_basin (ERA5) are already in mm/month
#   at this point ??? the ?? days_in_month scaling was applied to the raster
#   layers on lines 275 and 309 before values() was called.  No additional
#   unit conversion is needed here.
# ==============================================================================
cat("\n===== EXPORT: WB_basin_average_monthly.csv =====\n")

basin_cells_idx    <- which(basin_mask)
wb_pm_basin_mean   <- basin_weighted_mean(wb_basin,      basin_bw, cells = basin_cells_idx)
wb_era5_basin_mean <- basin_weighted_mean(wb_era5_basin, basin_bw, cells = basin_cells_idx)
wb_pm_coverage     <- basin_valid_coverage(wb_basin,      basin_bw, cells = basin_cells_idx)
wb_era5_coverage   <- basin_valid_coverage(wb_era5_basin, basin_bw, cells = basin_cells_idx)

wb_export <- data.frame(
  date            = format(dates, "%Y-%m-%d"),
  year            = year_nums,
  month           = month_nums,
  days_in_month   = days_in_month,
  wb_pm_mm_month  = round(wb_pm_basin_mean,  4L),
  wb_era5_mm_month = round(wb_era5_basin_mean, 4L),
  wb_pm_valid_area_coverage   = round(wb_pm_coverage,   4L),
  wb_era5_valid_area_coverage = round(wb_era5_coverage, 4L),
  stringsAsFactors = FALSE
)

wb_export_path <- file.path(out_dir, "WB_basin_average_monthly.csv")
write.csv(wb_export, wb_export_path, row.names = FALSE, quote = FALSE)

cat(sprintf("  Written : %s\n", wb_export_path))
cat(sprintf("  Rows    : %d  (%s to %s)\n",
            nrow(wb_export),
            wb_export$date[1],
            wb_export$date[nrow(wb_export)]))
cat(sprintf("  WB_PM   range: %.1f to %.1f mm/month\n",
            min(wb_export$wb_pm_mm_month,  na.rm = TRUE),
            max(wb_export$wb_pm_mm_month,  na.rm = TRUE)))
cat(sprintf("  WB_ERA5  range: %.1f to %.1f mm/month\n",
            min(wb_export$wb_era5_mm_month, na.rm = TRUE),
            max(wb_export$wb_era5_mm_month, na.rm = TRUE)))
cat("  (Diagnostic export; not required by other framework scripts.)\n")
cat("=================================================\n\n")

# ==============================================================================

# ==============================================================================
#   PARAMETRIC FITTING HELPER
# ==============================================================================
# Candidate 3-parameter distributions (lmomco type codes). GLO is listed first
# (the standard SPEI choice) so it wins exact ties.
SPEI_CANDIDATE_DISTS <- c("glo", "pe3", "gev")

# Best distribution = smallest L-moment-ratio distance |tau4_sample - tau4_fit|
# (lmom_rank_distributions() in DROUGHT_ANALYSIS_utils.R). All candidates are
# 3-parameter L-moment fits, which reproduce the sample tau3 exactly, so tau4
# is the discriminating ratio. If the winner's CDF is not finite for every
# value, the next-ranked candidate is used; NULL if none works.
try_parametric <- function(values_m) {
  rk <- lmom_rank_distributions(values_m, SPEI_CANDIDATE_DISTS)
  fits <- attr(rk, "fits")
  for (k in seq_len(nrow(rk))) {
    if (!is.finite(rk$distance[k])) next
    par <- fits[[rk$dist[k]]]
    p_fitted <- tryCatch(lmomco::par2cdf(values_m, par), error = function(e) NULL)
    if (!is.null(p_fitted) && length(p_fitted) == length(values_m) && all(is.finite(p_fitted)))
      return(p_fitted)
  }
  NULL
}

# ==============================================================================
#   VARIANCE-AWARE SPEI
# ==============================================================================
variance_aware_spei <- function(x, month_numbers, scale = 1) {
  x_agg <- if (scale > 1) {
    xa <- zoo::rollapply(x, scale, sum, align = "right", fill = NA, na.rm = FALSE)
    if (length(xa) < length(x)) c(rep(NA_real_, length(x) - length(xa)), xa) else xa
  } else x
  
  n      <- length(x_agg)
  z      <- rep(NA_real_, n)
  method <- rep(NA_integer_, n)   # 1=Parametric  2=ZeroVar  NA=Failed
  
  for (m in 1:12) {
    idx <- which(month_numbers == m & is.finite(x_agg))
    if (length(idx) < 5) next
    values_m <- x_agg[idx]
    v0       <- var(values_m, na.rm = TRUE)
    
    if (!is.finite(v0) || v0 < .Machine$double.eps) {
      z[idx] <- 0; method[idx] <- 2; next
    }
    p <- try_parametric(values_m)
    if (is.null(p)) {
      # L-moments invalid or all distributions failed ??? return NA
      method[idx] <- NA_integer_
    } else {
      method[idx] <- 1
    }
    if (!is.null(p)) z[idx] <- qnorm(pmax(pmin(p, 1 - 1e-6), 1e-6))
  }
  list(z = z, method = method)
}

# ==============================================================================
#   run_spei_loop()
#   Runs the full SPEI calculation for one WB input, writing all outputs to
#   target_dir. Called identically for SPEI_PM and SPEI_ERA5.
#
#   Arguments:
#     wb_input          : full raster matrix [n_pixels x n_months] water balance
#                         (mm/month); cells outside the basin are all NA
#     target_dir        : output directory (string)
#     pet_label         : label appended to filenames/messages ("_PM" or "_ERA5")
#     cl_in             : parallel cluster (already running)
#     wb_template       : SpatRaster for spatial reconstruction (wb[[1]] or wb_era5[[1]])
#     summary_file_path : path for the combined summary text file
# ==============================================================================
run_spei_loop <- function(wb_input, target_dir, pet_label,
                          cl_in, wb_template, summary_file_path) {
  
  # Enforce that every SPEI branch uses the actual native ERA5-Land raster
  # geometry established by the precipitation NetCDF.
  assert_native_grid(wb_template, analysis_template,
                     paste0("SPEI", pet_label, " template"))
  
  if (terra::ncell(wb_template) != terra::ncell(analysis_template)) {
    stop(
      "SPEI template cell count does not match the native ERA5-Land grid: ",
      terra::ncell(wb_template), " vs ",
      terra::ncell(analysis_template), "."
    )
  }
  
  if (nrow(wb_input) != terra::ncell(analysis_template)) {
    stop(
      "SPEI input matrix has ", nrow(wb_input),
      " rows but the native ERA5-Land grid has ",
      terra::ncell(analysis_template), " cells."
    )
  }
  
  
  cat(sprintf("SPEI%s SUMMARY (Variance-Aware)\nGenerated: %s\n%s\n\n",
              pet_label, Sys.time(), strrep("=", 60)),
      file = summary_file_path)
  
  scales <- c(1, 3, 6, 12, 24, 36)
  
  # ?????? FIX (part 1): export wb_input ONCE before the loop.
  # wb_input is constant across all scales; only `scale` changes per iteration.
  clusterExport(cl_in, varlist = "wb_input", envir = environment())
  
  for (scale in scales) {
    cat(sprintf("\n===== SPEI%s-%d (Variance-Aware Calculation) =====\n",
                pet_label, scale))
    
    # Export only the scalar that changes each iteration
    clusterExport(cl_in, varlist = "scale", envir = environment())
    
    # ?????? FIX (part 2): strip the worker function's closure.
    # parLapply serializes the function AND its enclosing environment to each
    # worker. If the function is defined as an anonymous lambda inside
    # run_spei_loop, its enclosing env is run_spei_loop's local frame, which
    # contains wb_input (~6 MB). Setting environment() to globalenv() makes
    # R serialize only a reference to globalenv instead of copying the 6 MB
    # object with every parLapply call. Workers resolve wb_input, scale, and
    # month_numbers from their own globalenv (populated by clusterExport).
    spei_worker <- function(i) variance_aware_spei(wb_input[i, ], month_numbers, scale)
    environment(spei_worker) <- globalenv()
    
    t0         <- Sys.time()
    pixel_list <- parLapply(cl_in, seq_len(nrow(wb_input)), spei_worker)
    cat(sprintf("  Parallel processing done (%.1f min)\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    
    res_basin    <- do.call(rbind, lapply(pixel_list, `[[`, "z"))
    method_basin <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
    storage.mode(method_basin) <- "integer"
    
    na_rate_basin <- 100 * mean(is.na(res_basin[basin_mask, , drop = FALSE]))
    cat(sprintf("NA rate (basin pixels): %.3f%%\n", na_rate_basin))
    cat("Variance threshold: v < 2.2e-16 -> SPEI=0; v >= 2.2e-16 -> Parametric (GLO/PE3/GEV, smallest L-moment-ratio tau4 distance)\n")
    
    # Parametric fit diagnostic (first valid basin pixel)
    cat("  Parametric fit check (representative basin pixel):\n")
    representative_pixel <- which(basin_mask)[1L]
    
    if (is.na(representative_pixel)) {
      cat("    No valid basin pixel available for diagnostic.\n")
    } else {
      x_rep <- wb_input[representative_pixel, ]
      
      for (m_d in 1:12) {
        idx_d <- which(month_numbers == m_d & is.finite(x_rep))
        if (length(idx_d) < 5) next
        
        vals_d <- if (scale > 1) {
          xa <- zoo::rollapply(
            x_rep,
            scale,
            sum,
            align = "right",
            fill = NA,
            na.rm = FALSE
          )
          xa <- xa[idx_d]
        } else {
          x_rep[idx_d]
        }
        
        vals_d <- vals_d[is.finite(vals_d)]
        if (length(vals_d) < 5) next
        
        v_d <- var(vals_d, na.rm = TRUE)
        
        if (!is.finite(v_d) || v_d < .Machine$double.eps) {
          cat(sprintf(
            "    %s: zero variance -> SPEI=0\n",
            month_names[m_d]
          ))
          next
        }
        
        lmom_d <- tryCatch(
          lmoms(vals_d, nmom = 4),
          error = function(e) NULL
        )
        
        if (is.null(lmom_d) || !are.lmom.valid(lmom_d)) {
          cat(sprintf(
            "    %s: L-moments invalid -> NA\n",
            month_names[m_d]
          ))
          next
        }
        
        rk_d <- lmom_rank_distributions(vals_d, SPEI_CANDIDATE_DISTS)
        
        cat(sprintf(
          "    %s: %s (|dtau4| = %s) parametric\n",
          month_names[m_d],
          toupper(rk_d$dist[1]),
          paste(sprintf("%s %.4f", toupper(rk_d$dist), rk_d$distance), collapse = ", ")
        ))
      }
    }
    
    # The SPEI result is already spatially complete because wb_input contains
    # every cell of the masked precipitation grid. Cells outside the basin
    # remain NA throughout the calculation.
    res_full <- res_basin
    
    # Spatial method map: most frequent method across valid months for each
    # raster cell. Outside-basin cells remain NA.
    method_full <- apply(
      method_basin,
      1,
      function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(NA_integer_)
        tab <- table(x)
        as.integer(names(tab)[which.max(tab)])
      }
    )
    
    # ------------------------------------------------------------------
    # OUTPUT 1: CSV (one per calendar month)
    # ------------------------------------------------------------------
    cat("  -> Saving CSV files...\n")
    for (m in 1:12) {
      idx_m <- which(month_numbers == m)
      if (length(idx_m) == 0) next
      df         <- as.data.frame(res_full[, idx_m, drop = FALSE])
      colnames(df) <- format(dates[idx_m], "%Y")
      coords     <- xyFromCell(wb_template, 1:nrow(df))
      df         <- cbind(lon = coords[, 1], lat = coords[, 2], df)
      write.csv(df,
                file.path(target_dir,
                          sprintf("spei_%02d_month%02d_%s.csv", scale, m, month_names[m])),
                row.names = FALSE, na = "")
    }
    cat(sprintf("  Saved 12 CSV files for SPEI%s-%d\n", pet_label, scale))
    
    # ------------------------------------------------------------------
    # OUTPUT 2: Excel workbook
    # ------------------------------------------------------------------
    cat("  -> Saving Excel workbook...\n")
    excel_data <- list()
    for (m in 1:12) {
      idx_m <- which(month_numbers == m)
      if (length(idx_m) == 0) next
      df         <- as.data.frame(res_full[, idx_m, drop = FALSE])
      colnames(df) <- format(dates[idx_m], "%Y")
      coords     <- xyFromCell(wb_template, 1:nrow(df))
      excel_data[[month_names[m]]] <- cbind(lon = coords[, 1], lat = coords[, 2], df)
    }
    xlsx_file <- file.path(target_dir, sprintf("spei_%02d_all_months.xlsx", scale))
    if (file.exists(xlsx_file)) file.remove(xlsx_file)
    Sys.sleep(0.3)
    write_xlsx(excel_data, xlsx_file)
    cat(sprintf("  Saved Excel: %s\n", xlsx_file))
    
    # ------------------------------------------------------------------
    # OUTPUT 3: Distribution map (PNG)
    # ------------------------------------------------------------------
    cat("  -> Creating distribution map...\n")
    spei_dist_raster <- rast(wb_template)
    values(spei_dist_raster) <- method_full
    mapping <- data.frame(code  = c(1, 2),
                          name  = c("Parametric", "Zero Var"),
                          color = c("#4575b4", "#91bfdb"),
                          stringsAsFactors = FALSE)
    png_file <- file.path(target_dir,
                          sprintf("spei_%02d_distribution_map.png", scale))
    png(png_file, width = 1200, height = 800, res = 150)
    plot(spei_dist_raster, col = mapping$color, breaks = c(0.5, 1.5, 2.5),
         legend = FALSE,
         main = sprintf("SPEI%s-%d: Parametric Method (GLO/PE3/GEV)", pet_label, scale),
         axes = FALSE, box = FALSE)
    if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
    valid_pixels <- sum(!is.na(method_full))
    dist_counts  <- table(factor(method_full[!is.na(method_full)],
                                 levels = 1:2, labels = mapping$name))
    legend("bottomright",
           legend = sprintf("%s (%d)", names(dist_counts), dist_counts),
           fill = mapping$color, title = "Method", cex = 0.85,
           bg = "white", bty = "o", border = "gray50")
    mtext(sprintf("Valid pixels: %d | Failed/NA: %d",
                  valid_pixels,
                  sum(is.na(method_full[basin_mask]))),
          side = 1, line = 0.5, cex = 0.75, col = "gray30")
    dev.off()
    cat(sprintf("  Saved distribution map: %s\n", png_file))
    
    # ------------------------------------------------------------------
    # OUTPUT 4: Summary statistics (append to text file)
    # ------------------------------------------------------------------
    cat("  -> Appending summary statistics...\n")
    summary_text <- capture.output({
      cat(sprintf("\n---------- SPEI%s-%d SUMMARY ----------\n", pet_label, scale))
      cat(sprintf("Date: %s\n\n", Sys.time()))
      cat(sprintf("Grid: %d x %d | Basin pixels: %d\n",
                  ncol(wb_template), nrow(wb_template), sum(basin_mask)))
      cat(sprintf("Period: %s to %s (%d months)\n",
                  min(dates), max(dates), length(dates)))
      cat("\nVariance Handling:\n")
      cat("  Case 1 zero var  (v < 2.2e-16): SPEI = 0\n")
      cat("  Case 2 suff var  (v >= 2.2e-16): Parametric (GLO/PE3/GEV, smallest L-moment-ratio tau4 distance); NA if fitting fails\n")
      cat(sprintf("\nNA rate (basin): %.3f%%  ", na_rate_basin))
      cat(ifelse(na_rate_basin < 0.5, "[MSPEI-ready]\n", "[CAUTION: >0.5% NAs]\n"))
      cat(sprintf("\nMethod distribution:\n"))
      cat(sprintf("  Parametric:  %d (%.1f%%)\n",
                  sum(method_full == 1, na.rm = TRUE),
                  100 * sum(method_full == 1, na.rm = TRUE) / valid_pixels))
      cat(sprintf("  Zero var:    %d (%.1f%%)\n",
                  sum(method_full == 2, na.rm = TRUE),
                  100 * sum(method_full == 2, na.rm = TRUE) / valid_pixels))
      cat(sprintf("  Failed/NA:   %d (%.1f%%)\n",
                  sum(is.na(method_full[basin_mask])),
                  100 * sum(is.na(method_full[basin_mask])) / valid_pixels))
      cat("\nCase counts by calendar month (basin pixels):\n")
      cat(sprintf("  %-5s  %8s  %6s  %8s  %6s\n",
                  "Month", "Param(n)", "Par(%)", "ZeroVar(n)", "ZV(%)"))
      for (mm in 1:12) {
        idx_mm   <- which(month_numbers == mm)
        meth_mm  <- as.integer(
          method_basin[
            basin_mask,
            idx_mm[idx_mm <= ncol(method_basin)],
            drop = FALSE
          ]
        )
        n_tot    <- length(meth_mm)
        if (n_tot == 0) next
        n1 <- sum(meth_mm == 1L, na.rm = TRUE)
        n2 <- sum(meth_mm == 2L, na.rm = TRUE)
        cat(sprintf("  %-5s  %8d  %5.1f%%  %8d  %5.1f%%\n",
                    month_names[mm],
                    n1, 100*n1/n_tot, n2, 100*n2/n_tot))
      }
      representative_pixel <- which(basin_mask)[1L]
      first_idx <- head(which(!is.na(res_full[representative_pixel, ])), 12)
      if (length(first_idx) >= 12) {
        fs <- res_full[basin_mask, first_idx, drop = FALSE]
        cat( "\nDrought/Wet frequency (first 12 months):\n ")
        cat(sprintf( "  W3 Extremely wet (SPEI ???  2.0): %.1f%%\n ",100 * mean(fs >=  2.0, na.rm = TRUE)))
        cat(sprintf( "  W2 Severely wet  (1.5 ??? SPEI < 2): %.1f%%\n ",100 * mean(fs >=  1.5 & fs < 2.0, na.rm = TRUE)))
        cat(sprintf( "  W1 Moderately wet (1 ??? SPEI < 1.5): %.1f%%\n ",100 * mean(fs >=  1.0 & fs < 1.5, na.rm = TRUE)))
        cat(sprintf( "  N0 Normal         (-1 < SPEI < 1):  %.1f%%\n ",100 * mean(fs > -1.0 & fs < 1.0, na.rm = TRUE)))
        cat(sprintf( "  D1 Moderately dry (-1.5 < SPEI ??? -1): %.1f%%\n ",100 * mean(fs <= -1.0 & fs > -1.5, na.rm = TRUE)))
        cat(sprintf( "  D2 Severely dry   (-2 < SPEI ??? -1.5): %.1f%%\n ",100 * mean(fs <= -1.5 & fs > -2.0, na.rm = TRUE)))
        cat(sprintf( "  D3 Extremely dry  (SPEI ??? -2.0): %.1f%%\n ",100 * mean(fs <= -2.0, na.rm = TRUE)))
      }
      cat("\n")
    })
    cat(summary_text, file = summary_file_path, append = TRUE, sep = "\n")
    
    # ------------------------------------------------------------------
    # OUTPUT 5: NetCDF (with retry logic for Windows file locks)
    # ------------------------------------------------------------------
    cat("  -> Saving NetCDF files...\n")
    MAX_RETRIES <- 3
    nc_failed   <- c()
    
    for (m in 1:12) {
      idx_m <- which(month_numbers == m)
      if (length(idx_m) == 0) next
      spei_rast <- rast(wb_template)
      spei_rast <- rep(spei_rast, length(idx_m))
      values(spei_rast) <- res_full[, idx_m, drop = FALSE]
      terra::time(spei_rast) <- dates[idx_m]
      
      nc_file <- file.path(target_dir,
                           sprintf("spei_%02d_month%02d_%s.nc",
                                   scale, m, month_names[m]))
      
      # Pre-emptive delete (releases Windows file lock before overwrite)
      if (file.exists(nc_file)) {
        del_ok <- tryCatch({ file.remove(nc_file); TRUE }, error = function(e) FALSE)
        if (!del_ok) {
          cat(sprintf("    WARNING: Cannot delete locked file %s -- skipped\n",
                      basename(nc_file)))
          cat("      Close it in QGIS/Panoply/ArcGIS and re-run.\n")
          nc_failed <- c(nc_failed, nc_file); next
        }
        Sys.sleep(0.1)
      }
      
      write_ok <- FALSE
      for (attempt in seq_len(MAX_RETRIES)) {
        ok <- tryCatch({
          suppressMessages(
            writeCDF(spei_rast, nc_file,
                     varname  = "spei",
                     longname = sprintf(
                       "Standardized Precipitation Evapotranspiration Index (SPEI%s-%d)",
                       pet_label, scale),
                     unit     = "standardized_index",
                     missval  = -9999,
                     overwrite = TRUE))
          # Verify the written NetCDF immediately. This prevents a spatial
          # geometry mismatch from propagating into downstream analyses.
          written_nc <- terra::rast(nc_file)
          
          if (!terra::compareGeom(
            written_nc,
            wb_template,
            stopOnError = FALSE,
            crs = TRUE,
            ext = TRUE,
            rowcol = TRUE,
            res = TRUE
          )) {
            stop(
              "Written SPEI NetCDF geometry does not match wb_template: ",
              nc_file
            )
          }
          
          if (terra::ncell(written_nc) != terra::ncell(analysis_template)) {
            stop(
              "Written SPEI NetCDF has ",
              terra::ncell(written_nc),
              " cells; native ERA5-Land grid has ",
              terra::ncell(analysis_template),
              "."
            )
          }
          
          TRUE
        }, error = function(e) {
          cat(sprintf("    WARNING: Attempt %d/%d failed: %s\n",
                      attempt, MAX_RETRIES, conditionMessage(e)))
          FALSE
        })
        if (isTRUE(ok)) { write_ok <- TRUE; break }
        if (attempt < MAX_RETRIES) Sys.sleep(1.5 * attempt)
      }
      if (!write_ok) {
        nc_failed <- c(nc_failed, nc_file)
        cat(sprintf("    SKIPPED (all %d attempts failed): %s\n",
                    MAX_RETRIES, basename(nc_file)))
      }
    }
    
    n_written <- 12 - length(nc_failed)
    if (length(nc_failed) == 0) {
      cat(sprintf("  Saved 12 NetCDF files for SPEI%s-%d\n", pet_label, scale))
    } else {
      cat(sprintf("  WARNING: Saved %d/12 NetCDF for SPEI%s-%d\n",
                  n_written, pet_label, scale))
      for (f in nc_failed) cat(sprintf("    Failed: %s\n", f))
    }
    
  }  # end for(scale in scales)
}  # end run_spei_loop()

# ==============================================================================
#   PARALLEL CLUSTER (shared by both SPEI runs)
# ==============================================================================
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\nStarting parallel cluster with %d cores\n", n_cores))

# ?????? FIX: redirect temp files to D: to avoid C: disk-full crash ?????????????????????????????????????????????
# makeCluster() writes per-worker log files to tempdir(), which sits on C:.
# When C: is full this raises "No space left on device" before a single worker
# starts.  Pointing TMPDIR to the project drive resolves it completely.
# outfile = "" echoes worker stdout to the console and avoids writing Rout
# log files to the temp dir on Windows.
.old_tmpdir <- Sys.getenv("TMPDIR")
.new_tmpdir <- file.path(dirname(WD_PATH), "tmp_r_cluster")   # sibling of the project root (same drive)
dir.create(.new_tmpdir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = .new_tmpdir)
cat(sprintf("  Temp dir redirected to: %s\n", .new_tmpdir))

cl <- makeCluster(n_cores, outfile = "")

clusterExport(cl, varlist = c("variance_aware_spei", "try_parametric",
                              "lmom_rank_distributions", "lmom_ratio_distance",
                              "LMOM_TAU3_TOLERANCE", "SPEI_CANDIDATE_DISTS",
                              "month_numbers"),
              envir = environment())
clusterEvalQ(cl, { library(zoo); library(lmomco) })

# ==============================================================================
#   RUN 1: SPEI_PM  (Penman-Monteith PET, full physics)
# ==============================================================================
cat("\n##############################################################\n")
cat("  RUN 1 OF 2: SPEI_PM  (Penman-Monteith PET, full physics)\n")
cat("  Output -> spei_results_seasonal/\n")
cat("##############################################################\n")

run_spei_loop(
  wb_input          = wb_mat,
  target_dir        = out_dir,
  pet_label         = "_PM",
  cl_in             = cl,
  wb_template       = wb[[1]],
  summary_file_path = file.path(out_dir, "spei_all_scales_summary.txt")
)

# ==============================================================================
#   RUN 2: SPEI_ERA5  (ERA5-Land native potential evaporation, PEV)
# ==============================================================================
cat("\n##############################################################\n")
cat("  RUN 2 OF 2: SPEI_ERA5 (ERA5-Land native potential evaporation, PEV)\n")
cat("  Output -> spei_results_seasonal_era5land/\n")
cat("  NOTE: PEV is an energy-balance product (not temperature-based) and is\n")
cat("  generally non-zero in winter; WB_era5 = P - PEV is dominated by PEV here.\n")
cat("##############################################################\n")

run_spei_loop(
  wb_input          = wb_era5_mat,
  target_dir        = out_dir_era5,
  pet_label         = "_ERA5",
  cl_in             = cl,
  wb_template       = wb_era5[[1]],
  summary_file_path = file.path(out_dir_era5, "spei_era5_all_scales_summary.txt")
)

# ==============================================================================

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\nParallel cluster stopped\n")

# Restore original TMPDIR (clean up after temp redirect)
if (nchar(.old_tmpdir) > 0) {
  Sys.setenv(TMPDIR = .old_tmpdir)
} else {
  Sys.unsetenv("TMPDIR")
}
cat(sprintf("  Temp dir restored to: %s\n",
            if (nchar(.old_tmpdir) > 0) .old_tmpdir else "(system default)"))

cat("\n============================================================\n")
cat("SPEI CALCULATION COMPLETE -- OBSERVED BRANCH (2 PET VARIANTS)\n")
cat("============================================================\n")
cat(sprintf("\nSPEI_PM   output : %s\n", normalizePath(out_dir)))
cat(sprintf("SPEI_ERA5  output : %s\n", normalizePath(out_dir_era5)))


# ==============================================================================
# ------------------------------------------------------------------------------
# Basin-average CSV export for the observed SPEI branches
# These files are consumed by downstream observed-framework analyses.
# ------------------------------------------------------------------------------

## ?????? Helper: read all 12 per-month CSVs for one scale ?? directory and ??????????????????????????????
## assemble a chronological data.frame(date, value) of basin-mean SPEI values.
##
## Arguments:
##   scale    : integer 1, 3, 6, 12, 24 or 36
##   dir      : path to the directory written by run_spei_loop()
##   mn_names : character(12) ??? month abbreviations matching file naming
##   all_dates: Date vector for the full record (length = total months)
##
## Returns a data.frame(date, value) sorted chronologically, one row per month.
assemble_basin_avg <- function(scale, dir, mn_names, all_dates) {
  
  # Pre-allocate output aligned to all_dates
  out_value <- rep(NA_real_, length(all_dates))
  
  for (m in seq_len(12L)) {
    csv_f <- file.path(dir,
                       sprintf("spei_%02d_month%02d_%s.csv", scale, m, mn_names[m]))
    if (!file.exists(csv_f)) {
      cat(sprintf("    WARNING: missing %s ??? skipping month %d\n",
                  basename(csv_f), m))
      next
    }
    df_m <- tryCatch(read.csv(csv_f, stringsAsFactors = FALSE, check.names = FALSE),
                     error = function(e) {
                       cat(sprintf("    ERROR reading %s: %s\n",
                                   basename(csv_f), e$message))
                       NULL
                     })
    if (is.null(df_m)) next
    
    # Identify year columns (4-digit numeric names, exclude lon/lat)
    yr_cols <- names(df_m)[grepl("^[0-9]{4}$", names(df_m))]
    if (!length(yr_cols)) next
    
    # Fractional-coverage area-weighted basin mean (shared 0.20 coverage rule).
    # run_spei_loop() writes one row per native raster cell in cell order; if a
    # file does not follow that layout, rows are mapped to cells via lon/lat.
    pixel_mat <- as.matrix(df_m[, yr_cols, drop = FALSE])
    row_w <- if (nrow(pixel_mat) == length(basin_bw$weights)) {
      basin_bw$weights
    } else {
      basin_weights_for_xy(df_m$lon, df_m$lat, analysis_template, basin_bw)
    }
    basin_means <- regional_weighted_mean(pixel_mat, row_w,
                                          min_coverage = BASIN_MIN_VALID_AREA_COVERAGE,
                                          total_weight = basin_bw$area_m2)
    names(basin_means) <- yr_cols
    
    # Map to all_dates positions (month == m)
    idx_m   <- which(format(all_dates, "%m") == sprintf("%02d", m))
    yr_here <- as.integer(format(all_dates[idx_m], "%Y"))
    
    for (k in seq_along(idx_m)) {
      yr_k <- as.character(yr_here[k])
      if (yr_k %in% names(basin_means)) {
        out_value[idx_m[k]] <- basin_means[[yr_k]]
      }
    }
  }
  
  data.frame(date  = format(all_dates, "%Y-%m-%d"),
             value = round(out_value, 6L),
             stringsAsFactors = FALSE)
}

## ?????? (A + B) Assemble and (C) write basin-average CSVs for scales 1,3,6,12,24,36 ???????????????????????????

EXPORT_SCALES <- c(1, 3, 6, 12, 24, 36)

# Output directory for the new ERA5 basin-average CSVs (same folder as the
# per-month CSVs so everything related to SPEI_ERA5 lives in one place; the
# files are also copied to the working directory so 12basin_timeseries.R can
# find them without path changes).
era5_avg_dir <- out_dir_era5   # "spei_results_seasonal_era5land"

basin_avg_pm   <- list()   # PM   basin-average data.frames (one per scale)
basin_avg_era5 <- list()   # ERA5 basin-average data.frames (one per scale)

for (sc in EXPORT_SCALES) {
  
  cat(sprintf("\n  Scale %d: assembling ERA5 basin-average series...\n", sc))
  
  # --- ERA5 ---
  df_era5 <- assemble_basin_avg(sc, out_dir_era5, month_names, dates)
  basin_avg_era5[[as.character(sc)]] <- df_era5
  
  out_csv_era5 <- file.path(era5_avg_dir,
                            sprintf("spei_era5_%02d_basin_average.csv", sc))
  write.csv(df_era5, out_csv_era5, row.names = FALSE, quote = FALSE)
  cat(sprintf("    Written: %s\n", out_csv_era5))
  
  # --- PM ---
  cat(sprintf("  Scale %d: assembling PM  basin-average series...\n", sc))
  df_pm <- assemble_basin_avg(sc, out_dir, month_names, dates)
  basin_avg_pm[[as.character(sc)]] <- df_pm
  
  out_csv_pm <- file.path(out_dir, sprintf("spei_pm_%02d_basin_average.csv", sc))
  write.csv(df_pm, out_csv_pm, row.names = FALSE, quote = FALSE)
  cat(sprintf("    Written: %s\n", out_csv_pm))
}

cat("\n  All SPEI_PM and SPEI_ERA5 basin-average CSVs written.\n")

cat("\nObserved-framework outputs complete.\n")
cat("Basin-average CSVs written for SPEI_PM and SPEI_ERA5.\n")
cat("============================================================\n")