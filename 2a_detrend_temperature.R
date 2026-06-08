##############################################
# SCRIPT 2a — MONTH-SPECIFIC TEMPERATURE DETRENDING
# Nechako River Basin Drought Analysis Framework
# Counterfactual Branch (proposed by S. Dery, JHM revision, June 2026)
#
# PURPOSE:
#   Applies month-specific ordinary-least-squares linear detrending to the
#   ERA5-Land 2 m air temperature (T2m) and 2 m dewpoint temperature (Tdew)
#   fields over the full 1950–2025 record.  The detrended fields are saved as
#   NetCDF files that feed into 2b_PET_ERALand.R, which reads them alongside
#   the observed temperatures and computes both observed and counterfactual
#   PET in a single run.
#
# METHOD: Month-Specific Linear Detrending
#   For each calendar month m ∈ {1 … 12} and each basin pixel, an OLS
#   linear regression of temperature vs. year is fitted over 1950–2025:
#
#       T(year, m) = a_m + b_m × year + ε_m
#
#   The detrended series is:
#
#       T_det(year, m) = T(year, m) − b_m × year
#                      ≡ T(year, m) − [T_fitted(year,m) − ā_m]
#
#   i.e. the linear trend is removed while the long-term monthly mean ā_m
#   is retained.  This preserves the natural inter-annual variance and the
#   seasonal amplitude; only the monotonic warming signal is removed.
#
#   The identical procedure is applied to Tdew so that
#
#       Tdew_det(year, m) ≤ T_det(year, m)   for all pixels and years
#
#   which guarantees VPD > 0 and RH < 100% — physical requirements for the
#   Penman-Monteith equation.  Any remaining Tdew_det > T_det instances
#   (can occur at pixels where the OLS fits produce a crossing) are clamped
#   to T_det as a physics guard.
#
# INPUTS (from monthly_data_direct/):
#   2m_temperature_monthly.nc             ERA5-Land T2m (K), 1950-01 to 2025-12
#   2m_dewpoint_temperature_monthly.nc    ERA5-Land Tdew (K), same time axis
#
# OUTPUTS (to monthly_data_direct/):
#   2m_temperature_detrended_monthly.nc   Detrended T2m (K), same grid/time
#   2m_dewpoint_detrended_monthly.nc      Detrended Tdew (K), same grid/time
#
# DIAGNOSTICS (to monthly_data_direct/detrend_diagnostics/):
#   slope_T2m_month{mm}.tif               Pixel-level b_m for T2m  (°C/year)
#   slope_Tdew_month{mm}.tif              Pixel-level b_m for Tdew (°C/year)
#   rsq_T2m_month{mm}.tif                 Pixel-level R² for T2m OLS fits
#   rsq_Tdew_month{mm}.tif                Pixel-level R² for Tdew OLS fits
#   detrend_diagnostics.pdf               Multi-panel diagnostic figures:
#       A. Basin-average obs vs. detrended T2m — 12-panel monthly overlay
#       B. Spatial warming rate maps — b_m (°C/decade) for each month
#       C. R² distributions — histograms of OLS goodness-of-fit per month
#       D. Residual autocorrelation check — lag-1 ACF of detrended residuals
#       E. Physics check — VPD preservation before/after detrending
#       F. QC summary table
#   detrend_summary.csv                   Per-month basin-mean statistics
#
# EXECUTION ORDER:
#   Step 1 → 2a_detrend_temperature.R  ← THIS SCRIPT
#   Step 2 → 2b_PET_ERALand.R   (reads observed + detrended T; computes both
#                                  observed PET and counterfactual PET in one run)
#   Step 3 → 3SPEI_ERALand.R Runs 1–4  (SPI/SPEI from observed + detrended PET)
#
# The numeric prefixes reflect execution order: 2a runs before 2b.
#
# NOTES:
#   1. Outputs are in Kelvin (same unit as ERA5-Land inputs) so that
#      2b_PET_ERALand.R Section 5c can load them without any unit
#      conversion — it subtracts 273.15 exactly as it does for observed T2m.
#   2. No setwd() call — caller controls the working directory.
#   3. safe_writeCDF() helper mirrors 2b_PET_ERALand.R to prevent Windows
#      file-lock errors when .nc files are open in QGIS / Panoply.
#   4. terraOptions(memfrac) is set conservatively; increase on machines
#      with > 32 GB RAM.
##############################################

library(terra)
library(ncdf4)
library(lubridate)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

MONTHLY_DIR  <- "D:/Nechako_Drought/Nechako/monthly_data_direct"
DIAG_DIR     <- file.path(MONTHLY_DIR, "detrend_diagnostics")
LOG_FILE     <- file.path(DIAG_DIR, "detrend_temperature.log")

# Create output directories
if (!dir.exists(DIAG_DIR)) dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

terraOptions(memfrac = 0.7, progress = 0)

# --- Logging helper ---
log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line      <- paste0(timestamp, " | ", msg, "\n")
  cat(line, file = LOG_FILE, append = file.exists(LOG_FILE))
  message(paste0(timestamp, " | ", msg))
}

cat("Month-Specific Temperature Detrending Log\n", file = LOG_FILE)
cat(paste("Processing started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)

# ==============================================================================
# safe_writeCDF: writeCDF wrapper with pre-emptive delete + retry.
# Mirrors implementation in 2b_PET_ERALand.R (Section 1).
# ==============================================================================
safe_writeCDF <- function(raster_obj, filepath, varname, longname, unit,
                          overwrite = TRUE, max_retries = 3) {
  if (file.exists(filepath)) {
    del_ok <- tryCatch({ file.remove(filepath); TRUE }, error = function(e) FALSE)
    if (!del_ok) {
      log_event(sprintf("  WARNING: Cannot delete locked file: %s", basename(filepath)))
      log_event("    Close it in QGIS/Panoply/ArcGIS and re-run.")
      return(invisible(FALSE))
    }
    Sys.sleep(0.15)
  }
  for (attempt in seq_len(max_retries)) {
    ok <- tryCatch({
      suppressMessages(
        writeCDF(raster_obj, filepath, varname = varname,
                 longname = longname, unit = unit, overwrite = overwrite))
      TRUE
    }, error = function(e) {
      log_event(sprintf("  WARNING: Write attempt %d/%d failed: %s",
                        attempt, max_retries, conditionMessage(e)))
      FALSE
    })
    if (isTRUE(ok)) return(invisible(TRUE))
    if (attempt < max_retries) Sys.sleep(1.5 * attempt)
  }
  log_event(sprintf("  ERROR: All %d attempts failed for: %s",
                    max_retries, basename(filepath)))
  invisible(FALSE)
}

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================
log_event("Loading ERA5-Land temperature fields...")

f_t2m  <- file.path(MONTHLY_DIR, "2m_temperature_monthly.nc")
f_tdew <- file.path(MONTHLY_DIR, "2m_dewpoint_temperature_monthly.nc")

tryCatch({
  t2m_b  <- rast(f_t2m)
  tdew_b <- rast(f_tdew)
  log_event("Successfully loaded T2m and Tdew.")
}, error = function(e) stop("Input file loading failed: ", e$message))

# --- Robust date extraction (mirrors Script 2) ---
extract_time_robust <- function(raster_obj, file_path, var_name) {
  time_info <- time(raster_obj)
  if ((inherits(time_info, "Date") || inherits(time_info, "POSIXt")) &&
      !all(is.na(time_info))) {
    log_event(sprintf("  Time extracted via terra::time() for %s", var_name))
    return(as.Date(time_info))
  }
  log_event(sprintf("  Trying ncdf4 for %s...", var_name))
  tryCatch({
    nc         <- nc_open(file_path)
    time_data  <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")$value
    m          <- regexpr("since\\s+([0-9]{4}-[0-9]{2}-[0-9]{2})", time_units)
    if (m > 0) {
      origin <- as.Date(sub("since\\s+", "", regmatches(time_units, m)))
      tv     <- if (grepl("hours", time_units, ignore.case = TRUE))
        origin + time_data / 24 else origin + time_data
      nc_close(nc)
      log_event(sprintf("  Time extracted via ncdf4 for %s", var_name))
      return(as.Date(tv))
    }
    nc_close(nc)
  }, error = function(e)
    log_event(sprintf("  ncdf4 failed for %s: %s", var_name, conditionMessage(e))))
  n  <- nlyr(raster_obj)
  tv <- seq.Date(as.Date("1950-01-01"), by = "month", length.out = n)
  log_event(sprintf("  Reconstructed time for %s (%d layers from 1950-01)", var_name, n))
  tv
}

dates    <- extract_time_robust(t2m_b, f_t2m, "T2m")
n_months <- nlyr(t2m_b)

log_event(paste("Total time steps:", n_months))
log_event(sprintf("Time range: %s to %s", format(min(dates)), format(max(dates))))

years_all   <- as.integer(format(dates, "%Y"))
months_all  <- as.integer(format(dates, "%m"))

# Verify Tdew has matching dimensions
if (nlyr(tdew_b) != n_months) {
  stop(sprintf("Tdew has %d layers but T2m has %d — inputs must match.",
               nlyr(tdew_b), n_months))
}

# Quick statistics on raw inputs
s <- global(t2m_b - 273.15,  c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Input T2m  (C): %.2f to %.2f (mean %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))
s <- global(tdew_b - 273.15, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Input Tdew (C): %.2f to %.2f (mean %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))

# ==============================================================================
# 3. MONTH-SPECIFIC OLS DETRENDING
#
# For each calendar month m:
#   (a) Extract the sub-stack of T2m layers belonging to month m  (n_years layers)
#   (b) Regress each pixel time series on year: lm(T_pixel ~ year)
#       b_m = pixel-wise OLS slope  (°C/year in K/year here, same thing)
#   (c) Detrended series: T_det = T - b_m * year
#       (retains long-term mean ā_m; removes linear trend only)
#   (d) Same procedure for Tdew
#
# Implementation strategy:
#   Use terra::lapp() pixel-wise for each month sub-stack.
#   The regressor (year vector) is encoded as an auxiliary argument.
#   Outputs are placed back into the full 912-layer stacks at the correct
#   time positions.
# ==============================================================================
log_event("==========================================")
log_event("Starting month-specific OLS detrending...")
log_event("==========================================")

# Pre-allocate output stacks (same type and extent as inputs, all NaN initially)
t2m_det  <- t2m_b  * NA
tdew_det <- tdew_b * NA

# Per-month diagnostic containers
slope_t2m_list  <- vector("list", 12)
slope_tdew_list <- vector("list", 12)
rsq_t2m_list    <- vector("list", 12)
rsq_tdew_list   <- vector("list", 12)

# Basin-average summary table
detrend_summary <- data.frame(
  month          = integer(12),
  month_name     = character(12),
  n_years        = integer(12),
  slope_T2m_mean = numeric(12),   # basin-mean b_m  (K/year)
  slope_T2m_sd   = numeric(12),
  slope_T2m_max  = numeric(12),
  slope_Tdew_mean= numeric(12),
  slope_Tdew_sd  = numeric(12),
  rsq_T2m_mean   = numeric(12),
  rsq_Tdew_mean  = numeric(12),
  T2m_obs_mean   = numeric(12),   # basin-mean raw T (°C) for context
  T2m_det_mean   = numeric(12),   # basin-mean detrended T (°C)
  stringsAsFactors = FALSE
)

for (m in 1:12) {
  
  log_event(sprintf("  Processing month %2d (%s)...", m, month.name[m]))
  
  # ----- Index of all time steps belonging to this calendar month -----
  idx_m  <- which(months_all == m)
  yrs_m  <- years_all[idx_m]           # year values for regression x-axis
  n_yrs  <- length(idx_m)
  
  if (n_yrs < 10L) {
    log_event(sprintf("    WARNING: Only %d years for month %d — skipping (need >= 10).", n_yrs, m))
    next
  }
  
  # ----- Sub-stacks (one layer per year within this calendar month) -----
  T_sub    <- t2m_b[[idx_m]]
  Tdew_sub <- tdew_b[[idx_m]]
  
  # ----- OLS fit per pixel using terra::regress() -----
  # terra::regress(y, x) fits y_pixel ~ x for each pixel in the stack y.
  # x must be a single numeric vector of length nlyr(y).
  # Returns a SpatRaster with layers: intercept (a_m), slope (b_m), R-squared.
  # NOTE: terra >= 1.7 supports regress(); older versions need lapp() workaround.
  # We use regress() as it is vectorized, GPU-friendly, and handles NA pixels.
  
  tryCatch({
    fit_t2m  <- regress(T_sub,    yrs_m, formula = "y~x")
    fit_tdew <- regress(Tdew_sub, yrs_m, formula = "y~x")
  }, error = function(e) {
    stop(sprintf("terra::regress() failed for month %d: %s", m, conditionMessage(e)))
  })
  
  # Use positional indexing to avoid terra version name mismatches
  a_t2m    <- fit_t2m[[1]]
  b_t2m    <- fit_t2m[[2]]
  
  a_tdew   <- fit_tdew[[1]]
  b_tdew   <- fit_tdew[[2]]
  # R-squared extraction/calculation
  # Many terra versions do not output R2 by default. If missing, we compute it mathematically.
  if (nlyr(fit_t2m) >= 3) {
    r2_t2m <- fit_t2m[[3]]
  } else {
    r2_t2m <- (b_t2m^2) * var(yrs_m,na.rm = TRUE) / app(T_sub, fun = var, na.rm = TRUE)
  }
  
  if (nlyr(fit_tdew) >= 3) {
    r2_tdew <- fit_tdew[[3]]
  } else {
    r2_tdew <- (b_tdew^2) * var(yrs_m) / app(Tdew_sub, fun = var, na.rm = TRUE)
  }
  # Store diagnostic rasters
  slope_t2m_list[[m]]  <- b_t2m
  slope_tdew_list[[m]] <- b_tdew
  rsq_t2m_list[[m]]    <- r2_t2m
  rsq_tdew_list[[m]]   <- r2_tdew
  
  # ----- Compute detrended series -----
  # T_det(year, m) = T(year, m) - b_m * year
  # This retains the long-term mean ā_m = a_m + b_m * mean(years_m) while
  # removing the trend b_m * (year - mean(years_m)) relative to mean(years_m).
  # Equivalently: T_det = T - (T_fitted - ā_m) = T - b_m * year + b_m * mean(year)
  # We use the simpler form T_det = T - b_m * year which shifts the mean of the
  # series to a_m (the fitted intercept at year=0) rather than ā_m.
  # To retain the actual long-term mean, we add back mean(T_fitted):
  #   mean_fitted = a_m + b_m * mean(years_m)
  # so: T_det = T - b_m * year + b_m * mean(years_m)
  #           = T - b_m * (year - mean(years_m))
  # This is the standard "remove slope, preserve mean" detrend.
  mean_yr <- mean(yrs_m)
  
  for (i in seq_along(idx_m)) {
    yr_i                  <- yrs_m[i]
    t2m_det[[idx_m[i]]]  <- T_sub[[i]]    - b_t2m  * (yr_i - mean_yr)
    tdew_det[[idx_m[i]]] <- Tdew_sub[[i]] - b_tdew * (yr_i - mean_yr)
  }
  
  # ----- Per-month diagnostics -----
  s_b_t2m  <- global(b_t2m,  c("mean", "sd", "max"), na.rm = TRUE)
  s_b_tdew  <- global(b_tdew, c("mean", "sd", "max"),  na.rm = TRUE)
  s_r2_t2m <- global(r2_t2m, "mean",                  na.rm = TRUE)
  s_r2_tdew<- global(r2_tdew,"mean",                  na.rm = TRUE)
  
  # Convert C/year to C/decade for logging
  log_event(sprintf("    T2m  slope: mean=%+.4f  sd=%.4f  max=%+.4f C/year  (= %+.3f C/decade)",
                    s_b_t2m$mean, s_b_t2m$sd, s_b_t2m$max,
                    s_b_t2m$mean * 10))
  log_event(sprintf("    Tdew slope: mean=%+.4f C/year  (= %+.3f C/decade)",
                    s_b_tdew$mean, s_b_tdew$mean * 10))
  log_event(sprintf("    R²: T2m mean=%.3f,  Tdew mean=%.3f",
                    s_r2_t2m[1, 1], s_r2_tdew[1, 1]))
  
  # Fill summary table
  s_obs  <- global(T_sub,              "mean", na.rm = TRUE)
  s_det  <- global(t2m_det[[idx_m]],   "mean", na.rm = TRUE)
  
  detrend_summary[m, ] <- list(
    month           = m,
    month_name      = month.name[m],
    n_years         = n_yrs,
    slope_T2m_mean  = s_b_t2m$mean,
    slope_T2m_sd    = s_b_t2m$sd,
    slope_T2m_max   = s_b_t2m$max,
    slope_Tdew_mean = s_b_tdew$mean,
    slope_Tdew_sd   = s_b_tdew$sd,
    rsq_T2m_mean    = s_r2_t2m[1, 1],
    rsq_Tdew_mean   = s_r2_tdew[1, 1],
    T2m_obs_mean    = mean(s_obs[, 1], na.rm = TRUE) - 273.15,
    T2m_det_mean    = mean(s_det[, 1], na.rm = TRUE) - 273.15
  )
  
  # Save per-month slope and R² GeoTIFFs
  writeRaster(b_t2m,  file.path(DIAG_DIR, sprintf("slope_T2m_month%02d.tif",  m)), overwrite = TRUE)
  writeRaster(b_tdew, file.path(DIAG_DIR, sprintf("slope_Tdew_month%02d.tif", m)), overwrite = TRUE)
  writeRaster(r2_t2m, file.path(DIAG_DIR, sprintf("rsq_T2m_month%02d.tif",   m)), overwrite = TRUE)
  writeRaster(r2_tdew,file.path(DIAG_DIR, sprintf("rsq_Tdew_month%02d.tif",  m)), overwrite = TRUE)
  
}  # end month loop

log_event("Month-specific OLS detrending complete.")

# ==============================================================================
# 4. PHYSICS CHECK: Tdew_det <= T_det
# ==============================================================================
log_event("==========================================")
log_event("Physics check: Tdew_det <= T_det ...")

vpd_viol <- global(tdew_det > t2m_det, "sum", na.rm = TRUE)
n_viol   <- sum(vpd_viol[, 1], na.rm = TRUE)

if (n_viol > 0L) {
  log_event(sprintf("  Clamping %d Tdew_det > T_det instances to preserve VPD > 0.", n_viol))
  tdew_det <- ifel(tdew_det > t2m_det, t2m_det, tdew_det)
  n_viol_after <- sum(global(tdew_det > t2m_det, "sum", na.rm = TRUE)[, 1], na.rm = TRUE)
  log_event(sprintf("  Violations after clamp: %d", n_viol_after))
} else {
  log_event("  PASS: Tdew_det <= T_det everywhere; no clamping required.")
}

# ==============================================================================
# 5. ASSIGN TIME AXIS AND LAYER NAMES
# ==============================================================================
terra::time(t2m_det)  <- dates
terra::time(tdew_det) <- dates

names(t2m_det)  <- sprintf("T2m_DET_%04d_%02d",  years_all, months_all)
names(tdew_det) <- sprintf("Tdew_DET_%04d_%02d", years_all, months_all)

# ==============================================================================
# 6. GLOBAL STATISTICS ON DETRENDED OUTPUTS
# ==============================================================================
log_event("==========================================")
log_event("Global statistics on detrended outputs:")

s_obs  <- global(t2m_b   - 273.15, c("min", "max", "mean"), na.rm = TRUE)
s_det  <- global(t2m_det - 273.15, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  T2m  obs  (°C): %.2f to %.2f (mean %.2f)",
                  min(s_obs$min), max(s_obs$max), mean(s_obs$mean)))
log_event(sprintf("  T2m  det  (°C): %.2f to %.2f (mean %.2f)",
                  min(s_det$min), max(s_det$max), mean(s_det$mean)))

s_obs  <- global(tdew_b   - 273.15, c("min", "max", "mean"), na.rm = TRUE)
s_det  <- global(tdew_det - 273.15, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  Tdew obs  (°C): %.2f to %.2f (mean %.2f)",
                  min(s_obs$min), max(s_obs$max), mean(s_obs$mean)))
log_event(sprintf("  Tdew det  (°C): %.2f to %.2f (mean %.2f)",
                  min(s_det$min), max(s_det$max), mean(s_det$mean)))

# ==============================================================================
# 7. WRITE OUTPUT NetCDF FILES
# ==============================================================================
log_event("==========================================")
log_event("Writing detrended NetCDF outputs...")

f_out_t2m  <- file.path(MONTHLY_DIR, "2m_temperature_detrended_monthly.nc")
f_out_tdew <- file.path(MONTHLY_DIR, "2m_dewpoint_detrended_monthly.nc")

safe_writeCDF(
  t2m_det, f_out_t2m,
  varname  = "t2m_detrended",
  longname = "2m Air Temperature (month-specific linearly detrended, K)",
  unit     = "K")
log_event(paste("Saved:", basename(f_out_t2m)))

safe_writeCDF(
  tdew_det, f_out_tdew,
  varname  = "d2m_detrended",
  longname = "2m Dewpoint Temperature (month-specific linearly detrended, K)",
  unit     = "K")
log_event(paste("Saved:", basename(f_out_tdew)))

# ==============================================================================
# 8. WRITE SUMMARY CSV
# ==============================================================================
write.csv(detrend_summary,
          file.path(DIAG_DIR, "detrend_summary.csv"),
          row.names = FALSE)
log_event("Summary CSV saved: detrend_summary.csv")

# ==============================================================================
# 9. DIAGNOSTIC PLOTS
# ==============================================================================
log_event("==========================================")
log_event("Generating diagnostic plots...")

# Compute basin-average monthly time series for obs and detrended T2m
basin_t2m_obs <- global(t2m_b   - 273.15, "mean", na.rm = TRUE)[, 1]
basin_t2m_det <- global(t2m_det - 273.15, "mean", na.rm = TRUE)[, 1]

# Helper to open PDF safely (mirrors Script 2 pattern)
tryCatch(
  pdf(file.path(DIAG_DIR, "detrend_diagnostics.pdf"), width = 14, height = 10),
  error = function(e) stop("Cannot write PDF: close it in your viewer first."))

# ---------------------------------------------------------------------------
# Helper: render a plot function into the open PDF AND save a PNG
# ---------------------------------------------------------------------------
save_plot <- function(name, plot_fn, width = 14, height = 10, res = 150) {
  plot_fn()
  png_path <- file.path(DIAG_DIR, paste0(name, ".png"))
  tryCatch({
    png(png_path, width = width, height = height, units = "in", res = res)
    plot_fn()
    dev.off()
    log_event(sprintf("    PNG saved: %s", basename(png_path)))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    log_event(sprintf("    WARNING: PNG failed for %s — %s", name, conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# Panel A: Basin-average observed vs. detrended T2m — 12-panel monthly overlay
# ---------------------------------------------------------------------------
log_event("... A: Basin-average obs vs. detrended T2m (12 months)")
save_plot("A_obs_vs_det_T2m_monthly", function() {
  
  par(mfrow = c(3, 4), mar = c(3.5, 4, 2.5, 1), oma = c(0, 0, 3, 0))
  for (m in 1:12) {
    idx_m  <- which(months_all == m)
    yrs_m  <- years_all[idx_m]
    t_obs  <- basin_t2m_obs[idx_m]
    t_det  <- basin_t2m_det[idx_m]
    
    y_range <- range(c(t_obs, t_det), na.rm = TRUE)
    y_pad   <- diff(y_range) * 0.12
    y_lim   <- c(y_range[1] - y_pad, y_range[2] + y_pad)
    
    plot(yrs_m, t_obs, type = "l", col = "steelblue", lwd = 1.5,
         ylim = y_lim, xlab = "", ylab = "T2m (°C)",
         main = month.name[m], cex.main = 1.1, cex.lab = 0.9, cex.axis = 0.85)
    lines(yrs_m, t_det, col = "firebrick2", lwd = 1.5, lty = 2)
    
    # Add OLS trend line for observed (visual confirmation of what was removed)
    lm_obs <- lm(t_obs ~ yrs_m)
    abline(lm_obs, col = "steelblue", lty = 3, lwd = 1.2)
    
    # Slope annotation (K/decade)
    b_decade <- coef(lm_obs)[2] * 10
    legend("topleft",
           legend = sprintf("%+.3f C/dec", b_decade),
           bty = "n", cex = 0.8, text.col = "steelblue")
    grid(lty = "dotted", col = "grey88")
  }
  mtext("Basin-Average T2m: Observed (blue) vs. Detrended (red dashed)",
        outer = TRUE, cex = 1.2, font = 2)
  legend("bottom",
         legend = c("Observed T2m", "Detrended T2m", "OLS trend (observed)"),
         col    = c("steelblue", "firebrick2", "steelblue"),
         lty    = c(1, 2, 3), lwd    = 2, bty = "n", horiz = TRUE, cex = 0.85,
         inset  = -0.08, xpd = NA)
})

# ---------------------------------------------------------------------------
# Panel B: Spatial warming rate maps — b_m (K/decade) for each calendar month
# ---------------------------------------------------------------------------
log_event("... B: Spatial warming rate maps (K/decade)")
save_plot("B_warming_rate_maps_T2m", function() {
  
  par(mfrow = c(3, 4), mar = c(1.5, 1.5, 2.5, 3.5), oma = c(0, 0, 3, 0))
  all_vals <- unlist(lapply(slope_t2m_list,
                            function(r) if (!is.null(r)) as.vector(values(r)) else NULL))
  sym_lim  <- max(abs(quantile(all_vals * 10, c(0.01, 0.99), na.rm = TRUE)))
  brks     <- seq(-sym_lim, sym_lim, length.out = 101)
  pal      <- colorRampPalette(c("#313695", "#74add1", "#ffffbf", "#f46d43", "#a50026"))(100)
  
  for (m in 1:12) {
    if (is.null(slope_t2m_list[[m]])) next
    plot(slope_t2m_list[[m]] * 10,
         col     = pal,
         breaks  = brks,
         legend  = TRUE,
         main    = month.name[m],
         cex.main = 1.0,
         axes    = FALSE)
    box()
  }
  mtext("T2m Warming Rate (C/decade, 1950–2025) — OLS slope per calendar month",
        outer = TRUE, cex = 1.1, font = 2)
})

# ---------------------------------------------------------------------------
# Panel C: R² distribution histograms — OLS goodness-of-fit per month
# ---------------------------------------------------------------------------
log_event("... C: R² distribution histograms")
save_plot("C_R2_distributions", function() {
  
  par(mfrow = c(3, 4), mar = c(3.5, 4, 2.5, 1), oma = c(0, 0, 3, 0))
  for (m in 1:12) {
    if (is.null(rsq_t2m_list[[m]])) next
    r2_vals <- as.vector(values(rsq_t2m_list[[m]]))
    r2_vals <- r2_vals[!is.na(r2_vals)]
    hist(r2_vals, breaks = 20, col = "#4393c3", border = "white",
         main = month.name[m], xlab = "R²", xlim = c(0, 1),
         cex.main = 1.1, cex.axis = 0.85, cex.lab = 0.9)
    abline(v = median(r2_vals, na.rm = TRUE), col = "red", lwd = 1.5)
    legend("topright",
           legend = sprintf("med=%.2f", median(r2_vals, na.rm = TRUE)),
           bty = "n", cex = 0.8, text.col = "red")
    grid(lty = "dotted", col = "grey88")
  }
  mtext("OLS R² per Calendar Month — T2m Linear Trend Fit (1950–2025)",
        outer = TRUE, cex = 1.1, font = 2)
})

# ---------------------------------------------------------------------------
# Panel D: Residual lag-1 autocorrelation check
#   Verifies detrended residuals are approximately white noise.
#   Uses basin-average time series for each calendar month.
# ---------------------------------------------------------------------------
log_event("... D: Residual autocorrelation check")
save_plot("D_residual_acf", function() {
  
  par(mfrow = c(3, 4), mar = c(3.5, 4, 2.5, 1), oma = c(0, 0, 3, 0))
  for (m in 1:12) {
    idx_m   <- which(months_all == m)
    yrs_m   <- years_all[idx_m]
    t_obs   <- basin_t2m_obs[idx_m]
    resid_m <- residuals(lm(t_obs ~ yrs_m))
    n_lag   <- min(15L, floor(length(resid_m) / 3))
    
    tryCatch({
      acf(resid_m, lag.max = n_lag, main = month.name[m],
          col = "steelblue", lwd = 2,
          cex.main = 1.0, cex.axis = 0.85, cex.lab = 0.9,
          ylab = "ACF", ci.col = "red")
      lag1 <- cor(resid_m[-length(resid_m)], resid_m[-1], use = "complete.obs")
      legend("topright",
             legend = sprintf("lag-1=%.2f", lag1),
             bty = "n", cex = 0.8)
    }, error = function(e)
      plot(1, type = "n", main = month.name[m], xlab = "", ylab = ""))
  }
  mtext("Residual ACF — T2m Detrended Basin-Average (white noise = good)",
        outer = TRUE, cex = 1.1, font = 2)
})

# ---------------------------------------------------------------------------
# Panel E: VPD preservation check
#   Plots basin-average VPD = es(T) - ea(Tdew) before and after detrending
#   for all 912 time steps (should always be positive after clamping).
# ---------------------------------------------------------------------------
log_event("... E: VPD preservation plot")
save_plot("E_VPD_preservation", function() {
  
  par(mfrow = c(1, 2), mar = c(4, 5, 3, 1.5), oma = c(0, 0, 2, 0))
  
  # Basin-average VPD from observed T / Tdew  (kPa)
  T_obs_C    <- global(t2m_b   - 273.15, "mean", na.rm = TRUE)[, 1]
  Tdew_obs_C <- global(tdew_b  - 273.15, "mean", na.rm = TRUE)[, 1]
  es_obs     <- 0.6108 * exp((17.27 * T_obs_C)    / (T_obs_C    + 237.3))
  ea_obs     <- 0.6108 * exp((17.27 * Tdew_obs_C) / (Tdew_obs_C + 237.3))
  vpd_obs    <- es_obs - ea_obs
  
  T_det_C    <- global(t2m_det  - 273.15, "mean", na.rm = TRUE)[, 1]
  Tdew_det_C <- global(tdew_det - 273.15, "mean", na.rm = TRUE)[, 1]
  es_det     <- 0.6108 * exp((17.27 * T_det_C)    / (T_det_C    + 237.3))
  ea_det     <- 0.6108 * exp((17.27 * Tdew_det_C) / (Tdew_det_C + 237.3))
  vpd_det    <- es_det - ea_det
  
  # Left panel: observed VPD time series
  plot(dates, vpd_obs, type = "l", col = "steelblue", lwd = 1.2,
       xlab = "Date", ylab = "Basin-Mean VPD (kPa)",
       main = "VPD — Observed T/Tdew", ylim = c(0, max(c(vpd_obs, vpd_det))))
  abline(h = 0, col = "red", lty = 2)
  grid(lty = "dotted", col = "grey88")
  legend("topleft", legend = sprintf("Min = %.4f kPa", min(vpd_obs, na.rm = TRUE)),
         bty = "n", cex = 0.9)
  
  # Right panel: detrended VPD time series
  plot(dates, vpd_det, type = "l", col = "firebrick2", lwd = 1.2,
       xlab = "Date", ylab = "Basin-Mean VPD (kPa)",
       main = "VPD — Detrended T/Tdew", ylim = c(0, max(c(vpd_obs, vpd_det))))
  abline(h = 0, col = "red", lty = 2)
  grid(lty = "dotted", col = "grey88")
  legend("topleft", legend = sprintf("Min = %.4f kPa", min(vpd_det, na.rm = TRUE)),
         bty = "n", cex = 0.9)
  
  mtext("VPD Preservation Check — VPD must remain ≥ 0 everywhere",
        outer = TRUE, cex = 1.0, font = 2)
})

# ---------------------------------------------------------------------------
# Panel F: Per-month detrending summary table
# ---------------------------------------------------------------------------
log_event("... F: QC summary table")
save_plot("F_detrend_QC_table", function() {
  
  par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  title("Detrending QC Summary — T2m OLS Linear Fit (1950–2025)", cex.main = 1.3)
  
  # Column headers
  cols <- c("Month", "N years", "Slope mean\n(C/decade)", "Slope max\n(C/decade)",
            "R² mean", "T_obs mean (°C)", "T_det mean (°C)")
  x_pos <- c(0.03, 0.17, 0.32, 0.49, 0.63, 0.76, 0.90)
  y0  <- 0.93
  ys  <- 0.065
  
  for (j in seq_along(cols))
    text(x_pos[j], y0, cols[j], adj = 0, cex = 0.72, font = 2)
  lines(c(0.02, 0.98), c(y0 - 0.015, y0 - 0.015), lwd = 1.5)
  
  for (i in 1:12) {
    y_row <- y0 - i * ys
    vals  <- c(
      detrend_summary$month_name[i],
      detrend_summary$n_years[i],
      sprintf("%+.3f", detrend_summary$slope_T2m_mean[i] * 10),
      sprintf("%+.3f", detrend_summary$slope_T2m_max[i]  * 10),
      sprintf("%.3f",  detrend_summary$rsq_T2m_mean[i]),
      sprintf("%.2f",  detrend_summary$T2m_obs_mean[i]),
      sprintf("%.2f",  detrend_summary$T2m_det_mean[i])
    )
    bg_col <- if (i %% 2 == 0) "grey96" else "white"
    rect(0.01, y_row - ys * 0.45, 0.99, y_row + ys * 0.45,
         col = bg_col, border = NA)
    for (j in seq_along(x_pos))
      text(x_pos[j], y_row, vals[j], adj = 0, cex = 0.72)
  }
})

dev.off()
log_event("Diagnostic PDF saved: detrend_diagnostics.pdf")

# ==============================================================================
# 10. FINISH
# ==============================================================================
log_event("==========================================")

# Final global summary
s_t2m_det  <- global(t2m_det  - 273.15, c("min","max","mean"), na.rm = TRUE)
s_tdew_det <- global(tdew_det - 273.15, c("min","max","mean"), na.rm = TRUE)

log_event(sprintf("Detrended T2m  (°C): %.2f to %.2f  mean=%.2f",
                  min(s_t2m_det$min), max(s_t2m_det$max), mean(s_t2m_det$mean)))
log_event(sprintf("Detrended Tdew (°C): %.2f to %.2f  mean=%.2f",
                  min(s_tdew_det$min), max(s_tdew_det$max), mean(s_tdew_det$mean)))
log_event(paste("Total months processed:", n_months))
log_event("Outputs:")
log_event(paste0("  ", basename(f_out_t2m)))
log_event(paste0("  ", basename(f_out_tdew)))
log_event(sprintf("  Diagnostic GeoTIFFs and PDF in: %s/", DIAG_DIR))
log_event("PROCESSING COMPLETE")
log_event(paste("Finished:", Sys.time()))