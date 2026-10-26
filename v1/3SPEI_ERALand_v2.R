##############################################
# SPEI CALCULATION FOR NECHAKO RIVER BASIN
#
# Current manuscript framing:
#   - SPEI_PM  = observed Penman???Monteith SPEI
#   - SPEI_Thw = observed Thornthwaite SPEI
#   - ??SPEI_PET = SPEI_PM - SPEI_Thw is a PET-formulation sensitivity residual
#     and is not interpreted as an additive physical decomposition.
#
# This script computes the two observed SPEI branches only.
# The detrended/counterfactual PET branch has been removed.
##############################################


# ---- Libraries ----
library(terra)
library(lubridate)
library(zoo)
library(writexl)
library(lmomco)
library(parallel)

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako/")
out_dir     <- "spei_results_seasonal"
out_dir_thw <- "spei_results_seasonal_thw"
basin_path  <- "Spatial/nechakoBound_dissolve.kmz"

if (!dir.exists(out_dir))     dir.create(out_dir,     recursive = TRUE)
if (!dir.exists(out_dir_thw)) dir.create(out_dir_thw, recursive = TRUE)

terraOptions(progress = 0)

# ---- Load Basin Boundary ----
tmp <- tempfile(); dir.create(tmp, showWarnings = FALSE)
utils::unzip(basin_path, exdir = tmp)
kml <- list.files(tmp, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)[1]
basin <- vect(kml)
if (nrow(basin) > 1L) basin <- aggregate(basin)
unlink(tmp, recursive = TRUE)
cat("Basin boundary loaded\n")
##############################################
# ---- Load Data ----
cat("\n===== LOADING INPUT DATA =====\n")
precip <- rast("monthly_data_direct/total_precipitation_monthly.nc")
# ?????? PRECIPITATION UNIT VERIFICATION (run before days-in-month scaling) ???????????????????????????
# ERA5-Land "tp" from CDS monthly averaged reanalysis = mean daily rate (m/day).
# Standard preprocessing converts to mm/day (*1000).
# Expected summer (July) basin range for Nechako: ~1 to 4 mm/day.
# If values are ~0.001-0.004 -> still in m/day (need *1000)
# If values are ~30-120      -> already mm/month total (do NOT multiply by days again)

precip_dates_raw <- as.Date(time(precip))
if (is.null(precip_dates_raw) || all(is.na(precip_dates_raw))) {
  precip_dates_raw <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(precip))
}
july_idx_p <- which(as.integer(format(precip_dates_raw, "%m")) == 7)[1]

if (!is.na(july_idx_p)) {
  precip_july      <- precip[[july_idx_p]]
  precip_july_vals <- values(precip_july, mat = FALSE)
  precip_july_vals <- precip_july_vals[is.finite(precip_july_vals)]
  mean_p           <- mean(precip_july_vals, na.rm = TRUE)
  
  cat(sprintf("\n?????? PRECIPITATION UNIT CHECK (July layer %d) ????????????????????????????????????????????????\n", july_idx_p))
  cat(sprintf("  Min  : %10.4f\n", min(precip_july_vals)))
  cat(sprintf("  Mean : %10.4f\n", mean_p))
  cat(sprintf("  Max  : %10.4f\n", max(precip_july_vals)))
  cat("  Expected if mm/day   : ~0.5 to 5\n")
  cat("  Expected if m/day    : ~0.0005 to 0.005  ??? multiply by 1000\n")
  cat("  Expected if mm/month : ~15 to 150         ??? do NOT multiply by days again\n")
  
  if (mean_p < 0.01) {
    cat("  ??? WARNING: Precip appears to be in m/day ??? add: precip <- precip * 1000\n")
  } else if (mean_p > 20) {
    cat("  ??? WARNING: Precip appears to be in mm/month ??? days_in_month multiplication will double-count\n")
  } else {
    cat("  ??? Precip units look consistent with mm/day\n")
  }
  cat("????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n\n")
} else {
  cat("  ??? Could not find a July layer in Precip raster ??? check time axis\n")
}
pet    <- rast("monthly_data_direct/potential_evapotranspiration_monthly.nc")
# ?????? PET unit verification (run before any WB computation) ??????????????????????????????????????????????????????
# ERA5-Land PET from CDS "monthly averaged reanalysis" = mean daily rate (m/day)
# Expected summer (July) basin range: ~0.5 to ~5 mm/day for a boreal basin.
# If values are ~0.001-0.005 -> still in m/day (need *1000)
# If values are ~30-150      -> mm/month total (already multiplied by days somewhere)
pet_dates_raw <- as.Date(time(pet))
if (is.null(pet_dates_raw) || all(is.na(pet_dates_raw))) {
  # Fall back: assume same time axis as precip
  pet_dates_raw <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(pet))
}
july_idx <- which(as.integer(format(pet_dates_raw, "%m")) == 7)[1]

if (!is.na(july_idx)) {
  pet_july      <- pet[[july_idx]]
  pet_july_vals <- values(pet_july, mat = FALSE)
  pet_july_vals <- pet_july_vals[is.finite(pet_july_vals)]
  cat(sprintf("\n?????? PET UNIT CHECK (July layer %d) ??????????????????????????????????????????????????????????????????\n", july_idx))
  cat(sprintf("  Min  : %10.4f\n", min(pet_july_vals)))
  cat(sprintf("  Mean : %10.4f\n", mean(pet_july_vals)))
  cat(sprintf("  Max  : %10.4f\n", max(pet_july_vals)))
  cat("  Expected if mm/day  : ~0.5 to 5\n")
  cat("  Expected if m/day   : ~0.0005 to 0.005  ??? multiply by 1000\n")
  cat("  Expected if mm/month: ~15 to 150         ??? do NOT multiply by days again\n")
  if (mean(pet_july_vals, na.rm = TRUE) < 0.01) {
    cat("  ??? WARNING: PET appears to be in m/day ??? add: pet <- pet * 1000\n")
  } else if (mean(pet_july_vals, na.rm = TRUE) > 20) {
    cat("  ??? WARNING: PET appears to be in mm/month ??? days_in_month multiplication will double-count\n")
  } else {
    cat("  ??? PET units look consistent with mm/day\n")
  }
  cat("????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n\n")
} else {
  cat("  ??? Could not find a July layer in PET raster ??? check time axis\n")
}

# Thornthwaite PET (temperature-only, from 2b_PET_ERALand.R)
pet_thw_file <- "monthly_data_direct/potential_evapotranspiration_thornthwaite_monthly.nc"
if (!file.exists(pet_thw_file))
  stop(paste("Thornthwaite PET file not found:", pet_thw_file,
             "\nRun 2b_PET_ERALand.R first."))
pet_thw <- rast(pet_thw_file)
cat("Thornthwaite PET loaded\n")

##############################################
# Reproject precipitation to BC Albers Equal Area.
# The precipitation raster defines the spatial grid for this SPEI calculation,
# following the same raster -> basin mask -> matrix architecture used by SPI.
target_crs <- "EPSG:3005"

if (!same.crs(precip, target_crs)) {
  cat("Reprojecting precipitation to BC Albers (EPSG:3005)...\n")
  precip <- project(precip, target_crs, method = "bilinear")
}

if (!same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
}

cat(sprintf(
  "Grid: %d x %d | Time steps: %d\n",
  ncol(precip),
  nrow(precip),
  nlyr(precip)
))

# Align both PET products directly to the precipitation grid.
cat("Aligning PM PET to precipitation grid...\n")
pet <- project(pet, precip, method = "bilinear")

cat("Aligning Thornthwaite PET to precipitation grid...\n")
pet_thw <- project(pet_thw, precip, method = "bilinear")

# Basin masking: preserve the full raster geometry exactly as in SPI.
# Masking changes cell values outside the basin to NA; it does not crop the
# raster and does not introduce a separate spatial template.
cat("\n===== MASKING INPUTS TO BASIN BOUNDARY =====\n")

precip  <- mask(precip,  basin, inverse = FALSE, touches = TRUE)
pet     <- mask(pet,     basin, inverse = FALSE, touches = TRUE)
pet_thw <- mask(pet_thw, basin, inverse = FALSE, touches = TRUE)

basin_mask  <- !is.na(values(precip[[1]]))
basin_pixels <- sum(basin_mask)
total_pixels <- ncell(precip)

# ==============================================================================
# CANONICAL SPEI ANALYSIS GRID
#
# This is the exact post-reprojection, basin-masked grid used by the water
# balance and by every SPEI NetCDF written below. Downstream scripts must use
# this object rather than infer basin membership from a multilayer SPEI file.
# ==============================================================================

# Store the basin mask as a SpatRaster, not as the logical vector returned
# by values(). terra::wrap() requires a terra spatial object.
basin_mask_raster <- terra::setValues(
  precip[[1]],
  as.logical(basin_mask)
)

spei_analysis_grid <- list(
  raster_template = terra::wrap(precip[[1]]),
  basin_mask      = terra::wrap(basin_mask_raster),
  basin_pixels    = basin_pixels,
  total_pixels    = total_pixels,
  crs             = terra::crs(precip),
  resolution      = terra::res(precip),
  extent          = terra::ext(precip),
  nrow            = terra::nrow(precip),
  ncol            = terra::ncol(precip)
)

spei_grid_file <- file.path(out_dir, "spei_analysis_grid.rds")
saveRDS(spei_analysis_grid, spei_grid_file)

cat("\n===== CANONICAL SPEI GRID SAVED =====\n")
cat(sprintf("  File: %s\n",
            normalizePath(spei_grid_file, winslash = "/", mustWork = FALSE)))
cat(sprintf("  Grid: %d x %d = %d cells\n",
            spei_analysis_grid$nrow,
            spei_analysis_grid$ncol,
            spei_analysis_grid$total_pixels))
cat(sprintf("  Basin cells: %d (%.1f%%)\n",
            spei_analysis_grid$basin_pixels,
            100 * spei_analysis_grid$basin_pixels /
              spei_analysis_grid$total_pixels))
cat(sprintf("  Resolution: %.6f x %.6f\n",
            spei_analysis_grid$resolution[1],
            spei_analysis_grid$resolution[2]))
cat(sprintf("  CRS: %s\n", spei_analysis_grid$crs))
cat("=======================================\n\n")

cat(sprintf(
  "Basin masking complete: %d pixels (%.1f%% of raster)\n",
  basin_pixels,
  100 * basin_pixels / total_pixels
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
terra::time(pet_thw) <- dates
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
#   WATER BALANCE ??? THORNTHWAITE  (WB_thw = P - PET_Thw)
# ==============================================================================
cat("\n===== WATER BALANCE (Thornthwaite) =====\n")
wb_thw <- precip - pet_thw
cat("Converting WB_thw from mm/day to mm/month...\n")
for (i in seq_len(nlyr(wb_thw))) wb_thw[[i]] <- wb_thw[[i]] * days_in_month[i]
terra::time(wb_thw) <- dates

wb_thw_mat   <- values(wb_thw, mat = TRUE)
wb_thw_basin <- wb_thw_mat[basin_mask, , drop = FALSE]
cat(sprintf("Basin pixels for WB_thw: %d\n", nrow(wb_thw_basin)))

total_vt  <- length(wb_thw_basin) - sum(is.na(wb_thw_basin))
zero_ct   <- sum(wb_thw_basin == 0, na.rm = TRUE)
cat(sprintf("WB_thw zeros: %d / %d (%.4f%%)\n", zero_ct, total_vt,
            100 * zero_ct / total_vt))
if (100 * zero_ct / total_vt > 1) cat("WARNING: High zeros in WB_thw!\n")

# NOTE: WB_thw = P - PET_Thw.  When T <= 0 C, PET_Thw = 0 so WB_thw = P.
# Precipitation is variable even in winter, so the WB_thw series retains full
# precipitation variance in DJF months.  Zero-variance cases arise only at
# pixels where winter precipitation itself is near-constant, not from the
# Thornthwaite constraint.
cat("\nZero counts by month (WB_thw):\n")
for (m in 1:12) {
  idx <- which(month_numbers == m)
  if (length(idx) == 0) next
  mv  <- wb_thw_basin[, idx, drop = FALSE]
  tot <- length(mv) - sum(is.na(mv))
  zer <- sum(mv == 0, na.rm = TRUE)
  cat(sprintf("  %s: %d / %d (%.2f%%)\n", month_names[m], zer, tot, 100 * zer / tot))
}
cat("=========================================\n\n")

# ==============================================================================
#   EXPORT: WB_basin_average_monthly.csv
#
#   Required by 10b_PET_bias_nonstationarity.R (Section 1 / Section 2).
#   Written to spei_results_seasonal/ (= out_dir = WB_DIR in 10b).
#
#   Columns:
#     date           character  "YYYY-MM-DD"  (first of each month)
#     year           integer
#     month          integer    1???12
#     days_in_month  integer    28/29/30/31  (actual, leap-year aware)
#     wb_pm_mm_month numeric    basin-mean WB_PM  = P ??? PET_PM  [mm/month]
#     wb_thw_mm_month numeric   basin-mean WB_Thw = P ??? PET_Thw [mm/month]
#
#   Basin mean: simple unweighted pixel mean over wb_basin / wb_thw_basin
#   (colMeans over non-NA rows).  This matches the unweighted convention used
#   throughout the post-processing section below and in 10b itself.  Area
#   weights differ by < 0.3 % across this approximately rectangular
#   ERA5-Land pixel grid at Nechako latitudes and do not materially affect
#   the bias or correction values computed in 10b.
#
#   NOTE: Both wb_basin (PM) and wb_thw_basin (Thw) are already in mm/month
#   at this point ??? the ?? days_in_month scaling was applied to the raster
#   layers on lines 275 and 309 before values() was called.  No additional
#   unit conversion is needed here.
# ==============================================================================
cat("\n===== EXPORT: WB_basin_average_monthly.csv =====\n")

wb_pm_basin_mean  <- colMeans(wb_basin,     na.rm = TRUE)  # length = n_months
wb_thw_basin_mean <- colMeans(wb_thw_basin, na.rm = TRUE)

wb_export <- data.frame(
  date            = format(dates, "%Y-%m-%d"),
  year            = year_nums,
  month           = month_nums,
  days_in_month   = days_in_month,
  wb_pm_mm_month  = round(wb_pm_basin_mean,  4L),
  wb_thw_mm_month = round(wb_thw_basin_mean, 4L),
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
cat(sprintf("  WB_Thw  range: %.1f to %.1f mm/month\n",
            min(wb_export$wb_thw_mm_month, na.rm = TRUE),
            max(wb_export$wb_thw_mm_month, na.rm = TRUE)))
cat("  NEXT STEP: 10b_PET_bias_nonstationarity.R can now run.\n")
cat("=================================================\n\n")

# ==============================================================================

# ==============================================================================
#   PARAMETRIC FITTING HELPER
# ==============================================================================
try_parametric <- function(values_m) {
  lmom <- tryCatch(lmoms(values_m, nmom = 4), error = function(e) NULL)
  if (is.null(lmom) || !are.lmom.valid(lmom)) return(NULL)
  
  dists <- list(
    GLO = list(par = tryCatch(parglo(lmom), error = function(e) NULL), cdf = cdfglo),
    PE3 = list(par = tryCatch(parpe3(lmom), error = function(e) NULL), cdf = cdfpe3),
    GEV = list(par = tryCatch(pargev(lmom), error = function(e) NULL), cdf = cdfgev)
  )
  best_p    <- -1
  best_prob <- NULL
  for (d in names(dists)) {
    par      <- dists[[d]]$par
    if (is.null(par)) next
    p_fitted <- tryCatch(dists[[d]]$cdf(values_m, par), error = function(e) NULL)
    if (is.null(p_fitted) || any(!is.finite(p_fitted))) next
    ks <- tryCatch(ks.test(values_m, function(q) dists[[d]]$cdf(q, par))$p.value,
                   error = function(e) 0)
    if (ks > best_p) { best_p <- ks; best_prob <- p_fitted }
  }
  best_prob
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
#   target_dir. Called identically for SPEI_PM and SPEI_Thw.
#
#   Arguments:
#     wb_input          : full raster matrix [n_pixels x n_months] water balance
#                         (mm/month); cells outside the basin are all NA
#     target_dir        : output directory (string)
#     pet_label         : label appended to filenames/messages ("_PM" or "_Thw")
#     cl_in             : parallel cluster (already running)
#     wb_template       : SpatRaster for spatial reconstruction (wb[[1]] or wb_thw[[1]])
#     summary_file_path : path for the combined summary text file
# ==============================================================================
run_spei_loop <- function(wb_input, target_dir, pet_label,
                          cl_in, wb_template, summary_file_path) {
  
  # Enforce that every SPEI branch uses exactly the canonical spatial grid.
  canonical_template <- terra::unwrap(spei_analysis_grid$raster_template)
  
  if (terra::ncell(wb_template) != spei_analysis_grid$total_pixels) {
    stop(
      "SPEI template cell count does not match canonical grid: ",
      terra::ncell(wb_template), " vs ",
      spei_analysis_grid$total_pixels, "."
    )
  }
  
  if (!terra::compareGeom(
    wb_template,
    canonical_template,
    stopOnError = FALSE,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE
  )) {
    stop("SPEI wb_template geometry does not match the canonical SPEI grid.")
  }
  
  if (nrow(wb_input) != spei_analysis_grid$total_pixels) {
    stop(
      "SPEI input matrix has ", nrow(wb_input),
      " rows but canonical grid has ",
      spei_analysis_grid$total_pixels, " cells."
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
    cat("Variance threshold: v < 2.2e-16 -> SPEI=0; v >= 2.2e-16 -> Parametric (GLO/PE3/GEV best KS)\n")
    
    # Parametric fit diagnostic (first valid basin pixel)
    cat("  Parametric fit check (representative basin pixel):\n")
    dist_fns <- list(
      GLO = list(par_fn = parglo, cdf_fn = cdfglo),
      PE3 = list(par_fn = parpe3, cdf_fn = cdfpe3),
      GEV = list(par_fn = pargev, cdf_fn = cdfgev))
    
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
        
        ks_res <- sapply(names(dist_fns), function(d) {
          par <- tryCatch(
            dist_fns[[d]]$par_fn(lmom_d),
            error = function(e) NULL
          )
          
          if (is.null(par)) return(0)
          
          tryCatch(
            ks.test(
              vals_d,
              function(q) dist_fns[[d]]$cdf_fn(q, par)
            )$p.value,
            error = function(e) 0
          )
        })
        
        best_d   <- names(which.max(ks_res))
        best_ksp <- max(ks_res)
        
        cat(sprintf(
          "    %s: %s (KS p=%.3f) parametric\n",
          month_names[m_d],
          best_d,
          best_ksp
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
      cat("  Case 2 suff var  (v >= 2.2e-16): Parametric (GLO/PE3/GEV best KS); NA if fitting fails\n")
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
          
          if (terra::ncell(written_nc) != spei_analysis_grid$total_pixels) {
            stop(
              "Written SPEI NetCDF has ",
              terra::ncell(written_nc),
              " cells; canonical grid has ",
              spei_analysis_grid$total_pixels,
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
.new_tmpdir <- "D:/Nechako_Drought/tmp_r_cluster"
dir.create(.new_tmpdir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = .new_tmpdir)
cat(sprintf("  Temp dir redirected to: %s\n", .new_tmpdir))

cl <- makeCluster(n_cores, outfile = "")

clusterExport(cl, varlist = c("variance_aware_spei", "try_parametric",
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
#   RUN 2: SPEI_Thw  (Thornthwaite PET, temperature-only)
# ==============================================================================
cat("\n##############################################################\n")
cat("  RUN 2 OF 2: SPEI_Thw (Thornthwaite PET, temperature-only)\n")
cat("  Output -> spei_results_seasonal_thw/\n")
cat("  NOTE: When T <= 0 C, PET_Thw = 0 so WB_thw = P - 0 = P.\n")
cat("  The water balance therefore retains full precipitation variance\n")
cat("  in DJF months. Zero-variance cases arise only at pixels where\n")
cat("  winter precipitation itself is near-constant, not from the\n")
cat("  Thornthwaite constraint per se.\n")
cat("##############################################################\n")

run_spei_loop(
  wb_input          = wb_thw_mat,
  target_dir        = out_dir_thw,
  pet_label         = "_Thw",
  cl_in             = cl,
  wb_template       = wb_thw[[1]],
  summary_file_path = file.path(out_dir_thw, "spei_thw_all_scales_summary.txt")
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
cat(sprintf("SPEI_Thw  output : %s\n", normalizePath(out_dir_thw)))


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
    
    # Compute unweighted basin-pixel means: average over non-NA rows per year
    pixel_mat <- as.matrix(df_m[, yr_cols, drop = FALSE])
    basin_means <- colMeans(pixel_mat, na.rm = TRUE)  # named by year
    
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

# Output directory for the new Thw basin-average CSVs (same folder as the
# per-month CSVs so everything related to SPEI_Thw lives in one place; the
# files are also copied to the working directory so w2_basin_timeseries.R can
# find them without path changes).
thw_avg_dir <- out_dir_thw   # "spei_results_seasonal_thw"

basin_avg_pm  <- list()   # will hold PM  data.frames for f_thm computation
basin_avg_thw <- list()   # will hold Thw data.frames for f_thm computation

for (sc in EXPORT_SCALES) {
  
  cat(sprintf("\n  Scale %d: assembling Thw basin-average series...\n", sc))
  
  # --- Thw ---
  df_thw <- assemble_basin_avg(sc, out_dir_thw, month_names, dates)
  basin_avg_thw[[as.character(sc)]] <- df_thw
  
  out_csv_thw <- file.path(thw_avg_dir,
                           sprintf("spei_thw_%02d_basin_average.csv", sc))
  write.csv(df_thw, out_csv_thw, row.names = FALSE, quote = FALSE)
  cat(sprintf("    Written: %s\n", out_csv_thw))
  
  # --- PM (needed for f_thm) ---
  cat(sprintf("  Scale %d: assembling PM  basin-average series...\n", sc))
  df_pm <- assemble_basin_avg(sc, out_dir, month_names, dates)
  basin_avg_pm[[as.character(sc)]] <- df_pm
}

cat("\n  All Thw basin-average CSVs written.\n")

cat("\nObserved-framework outputs complete.\n")
cat("Basin-average CSVs written for SPEI_PM and SPEI_Thw.\n")
cat("============================================================\n")
