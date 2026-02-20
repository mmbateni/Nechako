####################################################################################
# Trend Analysis Data Processing with Basin Averages
# Purpose: Perform comprehensive trend analysis for pixels AND basin averages
# Output: Binary data files including basin-averaged time series and statistics
####################################################################################
rm(list = ls())  # Clear workspace
gc()  # Garbage collect
library(ncdf4)
library(terra)
library(data.table)
library(Kendall)
library(future.apply)
library(zoo)
library(sf)
library(changepoint)  # NEW: For PELT changepoint detection

# ================= USER / ENV =================
setwd("D:/Nechako_Drought/Nechako/")

# ===== LOGGING SETUP =====
out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
LOG_FILE <- file.path(out_dir, "data_processing.log")
log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, "  ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, "  ", msg))
}
cat("Trend Analysis Data Processing - Started\n", file = LOG_FILE)
cat(paste("Timestamp:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat("==========================================\n", file = LOG_FILE, append = TRUE)

# ===== STEP 1: LOAD BASIN BOUNDARY =====
log_event("Searching for Nechako Basin boundary shapefile...")
basin_boundary <- NULL
basin_files <- c(
  "nechako_basin.shp", "Nechako_Basin.shp", "basin_boundary.shp",
  "../nechako_basin.shp", "data/nechako_basin.shp",
  "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
)
for (bf in basin_files) {
  if (file.exists(bf)) {
    tryCatch({
      basin_boundary <- st_read(bf, quiet = TRUE)
      target_crs <- "EPSG:3005" # BC Albers
      basin_boundary <- st_transform(basin_boundary, target_crs)
      log_event(paste("✓ Loaded Nechako Basin boundary from:", bf))
      log_event(paste("  Basin area:", round(as.numeric(st_area(basin_boundary))/1e6, 2), "km²"))
      break
    }, error = function(e) {
      log_event(paste("  Error loading", bf, ":", e$message))
    })
  }
}
if (is.null(basin_boundary)) {
  stop("CRITICAL: Nechako Basin boundary NOT FOUND. Required for clipping.\nPlease place shapefile in working directory.")
}
# Save basin boundary for visualization script
saveRDS(basin_boundary, file.path(out_dir, "basin_boundary.rds"))
log_event("Basin boundary saved for visualization script")

# ===== INPUT PATHS =====
precip_path <- "monthly_data_direct/total_precipitation_monthly.nc"
pet_path    <- "monthly_data_direct/potential_evapotranspiration_monthly.nc"

# ===== PARAMETERS =====
alpha <- 0.05
n_sim_spectral <- 500
max_tie_percent <- 50
max_min_value_pct_precip <- 80
max_min_value_pct_pet    <- 50
min_positive_value <- 0.01
min_obs_changepoint <- 20  # NEW: Minimum observations for PELT changepoint detection

# Parallel setup
num_cores <- min(parallel::detectCores() - 1, 8)
if (num_cores < 1) num_cores <- 1
plan(multisession, workers = num_cores)
log_event(paste("Using", num_cores, "cores for parallel processing"))

# ===== HELPER FUNCTIONS =====
check_min_value_threshold <- function(ts_clean, min_val_threshold = 0.01, max_pct = 50,
                                      var_name = "Unknown", is_precip = FALSE) {
  n_total <- length(ts_clean)
  n_min_vals <- sum(ts_clean <= min_val_threshold, na.rm = TRUE)
  pct_min_vals <- (n_min_vals / n_total) * 100
  max_allowed <- if (is_precip) max_min_value_pct_precip else max_min_value_pct_pet
  exceeds_threshold <- pct_min_vals > max_allowed
  list(exceeds_threshold = exceeds_threshold, pct_min_vals = pct_min_vals)
}

calculate_sens_slope_manual <- function(x) {
  n <- length(x)
  if (n < 2) return(NA_real_)
  slopes <- numeric()
  for (i in 1:(n-1)) {
    xi <- x[i]
    if (is.na(xi)) next
    for (j in (i+1):n) {
      xj <- x[j]
      if (!is.na(xj)) {
        slopes <- c(slopes, (xj - xi) / (j - i))
      }
    }
  }
  if (!length(slopes)) return(NA_real_)
  median(slopes, na.rm = TRUE)
}

calculate_variance_with_ties <- function(S, n, x) {
  tie_table  <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  var_s <- n * (n - 1) * (2 * n + 5) / 18
  if (length(tie_counts)) {
    tie_adjustment <- sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
    var_s <- var_s - tie_adjustment
  }
  var_s
}

# ===== NEW: PELT CHANGEPOINT DETECTION FUNCTION =====
# Detects regime shifts in mean level using PELT algorithm
detect_regime_shift_pelt <- function(ts_vec, min_obs = 20, start_year = NULL) {
  ts_clean <- na.omit(ts_vec)
  n <- length(ts_clean)
  
  # Require minimum observations for reliable detection
  if (n < min_obs) {
    return(list(
      changepoint_detected = FALSE,
      changepoint_position = NA_integer_,
      n_changepoints = 0L,
      first_changepoint_year = NA_integer_,
      mean_before = NA_real_,
      mean_after = NA_real_,
      magnitude_shift = NA_real_
    ))
  }
  
  result <- tryCatch({
    # PELT algorithm with BIC penalty (conservative, avoids over-detection)
    cpt <- cpt.mean(ts_clean, method = "PELT", penalty = "BIC")
    cpts_detected <- cpts(cpt)
    
    if (length(cpts_detected) > 0) {
      # Get first changepoint position
      first_cpt <- cpts_detected[1]
      
      # Calculate means before and after first changepoint
      mean_before <- mean(ts_clean[1:first_cpt], na.rm = TRUE)
      mean_after <- mean(ts_clean[(first_cpt+1):n], na.rm = TRUE)
      magnitude <- mean_after - mean_before
      
      # Convert position to actual year if start_year provided
      cpt_year <- if (!is.null(start_year)) {
        start_year + first_cpt - 1
      } else {
        first_cpt
      }
      
      list(
        changepoint_detected = TRUE,
        changepoint_position = first_cpt,
        n_changepoints = length(cpts_detected),
        first_changepoint_year = as.integer(cpt_year),
        mean_before = mean_before,
        mean_after = mean_after,
        magnitude_shift = magnitude
      )
    } else {
      list(
        changepoint_detected = FALSE,
        changepoint_position = NA_integer_,
        n_changepoints = 0L,
        first_changepoint_year = NA_integer_,
        mean_before = NA_real_,
        mean_after = NA_real_,
        magnitude_shift = NA_real_
      )
    }
  }, error = function(e) {
    # Return NA structure if changepoint detection fails
    list(
      changepoint_detected = FALSE,
      changepoint_position = NA_integer_,
      n_changepoints = 0L,
      first_changepoint_year = NA_integer_,
      mean_before = NA_real_,
      mean_after = NA_real_,
      magnitude_shift = NA_real_
    )
  })
  
  return(result)
}

aggregate_to_annual_fast <- function(monthly_matrix, years, method = "sum") {
  y_levels <- unique(years)
  yfac <- match(years, y_levels)
  count_by_year <- rowsum((!is.na(monthly_matrix)) * 1L, group = yfac, reorder = FALSE)
  m0 <- monthly_matrix
  m0[is.na(m0)] <- 0
  sum_by_year <- rowsum(m0, group = yfac, reorder = FALSE)
  if (identical(method, "sum")) {
    res <- sum_by_year
  } else {
    res <- sum_by_year / pmax(count_by_year, 1)
  }
  res[count_by_year < 6] <- NA_real_
  rownames(res) <- y_levels
  res
}

compute_month_index <- function(months) {
  split(seq_along(months), months)
}

mk_tfpw_spectral_for_series <- function(ts_vec, is_precip, alpha, max_tie_pct,
                                        n_sim_spectral, conf_cache_env, start_year = NULL) {
  ts_clean <- na.omit(ts_vec)
  n <- length(ts_clean)
  
  if (n < 10) {
    return(list(
      vc = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                n = n, rho1 = NA_real_, vc_corrected = FALSE, n_ties = 0, percent_ties = 0,
                n_min_vals = 0, percent_min_vals = 0, tau_b_adjusted = FALSE,
                filtered = TRUE, reason = "low_n"),
      tf = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                n = n, rho1 = NA_real_, tfpw_applied = FALSE, n_ties = 0, percent_ties = 0,
                n_min_vals = 0, percent_min_vals = 0, tau_b_adjusted = FALSE,
                filtered = TRUE, reason = "low_n"),
      spec = list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_),
      cpt = list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                 n_changepoints = 0L, first_changepoint_year = NA_integer_,
                 mean_before = NA_real_, mean_after = NA_real_, magnitude_shift = NA_real_)
    ))
  }
  
  min_val_threshold <- if (is_precip) 0.0 else min_positive_value
  min_check <- check_min_value_threshold(
    ts_clean, min_val_threshold,
    max_pct = if (is_precip) max_min_value_pct_precip else max_min_value_pct_pet,
    is_precip = is_precip
  )
  if (min_check$exceeds_threshold) {
    return(list(
      vc = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                n = n, rho1 = NA_real_, vc_corrected = FALSE, n_ties = 0, percent_ties = 0,
                n_min_vals = 0, percent_min_vals = min_check$pct_min_vals, tau_b_adjusted = FALSE,
                filtered = TRUE, reason = "excessive_min_vals"),
      tf = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                n = n, rho1 = NA_real_, tfpw_applied = FALSE, n_ties = 0, percent_ties = 0,
                n_min_vals = 0, percent_min_vals = min_check$pct_min_vals, tau_b_adjusted = FALSE,
                filtered = TRUE, reason = "excessive_min_vals"),
      spec = list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_),
      cpt = list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                 n_changepoints = 0L, first_changepoint_year = NA_integer_,
                 mean_before = NA_real_, mean_after = NA_real_, magnitude_shift = NA_real_)
    ))
  }
  
  n_unique <- length(unique(ts_clean))
  n_ties <- n - n_unique
  percent_ties <- (n_ties / n) * 100
  if (percent_ties > max_tie_pct) {
    return(list(
      vc = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                n = n, rho1 = NA_real_, vc_corrected = FALSE, n_ties = n_ties, percent_ties = percent_ties,
                n_min_vals = 0, percent_min_vals = min_check$pct_min_vals, tau_b_adjusted = FALSE,
                filtered = TRUE, reason = "excessive_ties"),
      tf = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                n = n, rho1 = NA_real_, tfpw_applied = FALSE, n_ties = n_ties, percent_ties = percent_ties,
                n_min_vals = 0, percent_min_vals = min_check$pct_min_vals, tau_b_adjusted = FALSE,
                filtered = TRUE, reason = "excessive_ties"),
      spec = list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_),
      cpt = list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                 n_changepoints = 0L, first_changepoint_year = NA_integer_,
                 mean_before = NA_real_, mean_after = NA_real_, magnitude_shift = NA_real_)
    ))
  }
  
  sen_slope <- calculate_sens_slope_manual(ts_clean)
  if (is.na(sen_slope) || is.infinite(sen_slope)) {
    return(list(
      vc = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                n = n, rho1 = NA_real_, vc_corrected = FALSE, n_ties = n_ties, percent_ties = percent_ties,
                n_min_vals = 0, percent_min_vals = min_check$pct_min_vals, tau_b_adjusted = FALSE,
                filtered = TRUE, reason = "sens_slope_fail"),
      tf = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                n = n, rho1 = NA_real_, tfpw_applied = FALSE, n_ties = n_ties, percent_ties = percent_ties,
                n_min_vals = 0, percent_min_vals = min_check$pct_min_vals, tau_b_adjusted = FALSE,
                filtered = TRUE, reason = "sens_slope_fail"),
      spec = list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_),
      cpt = list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                 n_changepoints = 0L, first_changepoint_year = NA_integer_,
                 mean_before = NA_real_, mean_after = NA_real_, magnitude_shift = NA_real_)
    ))
  }
  
  time_index <- seq_len(n)
  trend_line <- sen_slope * time_index
  detrended <- ts_clean - trend_line
  
  acf_result <- tryCatch(acf(detrended, lag.max = 1, plot = FALSE, na.action = na.pass),
                         error = function(e) NULL, warning = function(w) NULL)
  rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA_real_
  
  # Variance-corrected MK (VC)
  s_mat <- sign(outer(ts_clean, ts_clean, `-`))
  S_vc <- sum(s_mat[upper.tri(s_mat)], na.rm = TRUE)
  varS_taub <- calculate_variance_with_ties(S_vc, n, ts_clean)
  n_pairs <- n * (n - 1) / 2
  tau_vc <- S_vc / n_pairs
  tau_b_adjusted <- (percent_ties > 5)
  
  vc_corrected <- FALSE
  varS_final <- varS_taub
  if (!is.na(rho1) && abs(rho1) > 0.1) {
    correction_factor <- 1 + (2 * rho1 * (n - 1 - 2 * (n - 1) * rho1 + 3 * rho1 * rho1)) /
      ((n - 1) * (1 - rho1) * (1 - rho1))
    varS_final <- varS_taub * correction_factor
    vc_corrected <- TRUE
  }
  p_vc <- if (varS_final <= 0) NA_real_ else {
    z_stat <- S_vc / sqrt(varS_final)
    2 * pnorm(-abs(z_stat))
  }
  
  vc_list <- list(
    tau = tau_vc, p = p_vc, sl = sen_slope, S = S_vc, varS = varS_final,
    n = n, rho1 = rho1, vc_corrected = vc_corrected,
    n_ties = n_ties, percent_ties = percent_ties,
    n_min_vals = 0, percent_min_vals = min_check$pct_min_vals,
    tau_b_adjusted = tau_b_adjusted, filtered = FALSE, reason = "none"
  )
  
  # TFPW MK
  tfpw_applied <- FALSE
  if (!is.na(rho1) && abs(rho1) > 0.1) {
    pw <- numeric(n)
    pw[1] <- detrended[1]
    for (j in 2:n) pw[j] <- detrended[j] - rho1 * detrended[j-1]
    corrected <- pw + trend_line
    tfpw_applied <- TRUE
  } else {
    corrected <- ts_clean
  }
  
  s_mat_tf <- sign(outer(corrected, corrected, `-`))
  S_tf <- sum(s_mat_tf[upper.tri(s_mat_tf)], na.rm = TRUE)
  varS_tf <- calculate_variance_with_ties(S_tf, n, corrected)
  tau_tf <- S_tf / n_pairs
  p_tf <- if (varS_tf <= 0) NA_real_ else {
    z_stat <- S_tf / sqrt(varS_tf)
    2 * pnorm(-abs(z_stat))
  }
  
  tf_list <- list(
    tau = tau_tf, p = p_tf, sl = sen_slope, S = S_tf, varS = varS_tf,
    n = n, rho1 = rho1, tfpw_applied = tfpw_applied,
    n_ties = n_ties, percent_ties = percent_ties,
    n_min_vals = 0, percent_min_vals = min_check$pct_min_vals,
    tau_b_adjusted = tau_b_adjusted, filtered = FALSE, reason = "none"
  )
  
  # Spectral analysis
  dstd <- detrended - mean(detrended)
  sd_d <- stats::sd(dstd)
  if (!is.finite(sd_d) || sd_d == 0) sd_d <- 1
  dstd <- dstd / sd_d
  
  key <- paste0("n_", n)
  if (!exists(key, envir = conf_cache_env, inherits = FALSE)) {
    max_spectra <- numeric(n_sim_spectral)
    half <- floor(n / 2)
    for (jj in seq_len(n_sim_spectral)) {
      r <- rnorm(n, 0, 1)
      fr <- fft(r - mean(r))
      sp <- Mod(fr[1:half])^2 / n
      max_spectra[jj] <- max(sp, na.rm = TRUE)
    }
    assign(key, stats::quantile(max_spectra, 1 - alpha, na.rm = TRUE), envir = conf_cache_env)
  }
  conf_limit <- get(key, envir = conf_cache_env, inherits = FALSE)
  
  half <- floor(n / 2)
  ff <- fft(dstd)
  spectral_density <- Mod(ff[1:half])^2 / n
  freqs <- seq(0, 0.5, length.out = half)
  significant <- spectral_density > conf_limit
  peak_idx <- which(significant & !is.na(significant))
  peak_periods <- if (length(peak_idx)) {
    pf <- freqs[peak_idx]
    pp <- ifelse(pf > 0, 1/pf, NA_real_)
    ord <- order(spectral_density[peak_idx], decreasing = TRUE)
    pp[ord]
  } else numeric(0)
  
  spec_list <- list(
    n_peaks = length(peak_periods),
    dominant_period = if (length(peak_periods)) peak_periods[1] else NA_real_,
    conf = conf_limit
  )
  
  # NEW: Apply PELT changepoint detection to the time series
  cpt_result <- detect_regime_shift_pelt(ts_clean, min_obs = min_obs_changepoint, start_year = start_year)
  
  list(vc = vc_list, tf = tf_list, spec = spec_list, cpt = cpt_result)
}

# ===== DATA LOADING & PREPROCESSING =====
log_event("Loading precipitation and PET data...")
precip_full <- rast(precip_path)
pet_full    <- rast(pet_path)



# ===== EXTRACT TIME AND COMPUTE DAYS PER MONTH =====
log_event("Extracting time dimension for full rasters...")
dates_full <- NULL
tryCatch({
  dates_full <- terra::time(precip_full)
  if (!is.null(dates_full) && length(dates_full) > 0) {
    log_event("Time extracted using terra::time()")
  }
}, error = function(e) {})
if (is.null(dates_full) || length(dates_full) == 0 || all(is.na(dates_full))) {
  n_layers <- nlyr(precip_full)
  start_year <- 1950  # Changed from 1980 to match actual data range
  dates_full <- seq(as.Date(paste0(start_year, "-01-01")), by = "month", length.out = n_layers)
  log_event(paste("Generated", n_layers, "monthly dates from", start_year))
}
years_full  <- as.integer(format(dates_full, "%Y"))
months_full <- as.integer(format(dates_full, "%m"))

log_event("Computing days per month for each layer...")
month_days_base <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
is_leap <- function(year) (year %% 4 == 0) & (year %% 100 != 0 | year %% 400 == 0)
days_in_month <- month_days_base[months_full]
leap_years <- is_leap(years_full)
days_in_month[leap_years & months_full == 2] <- 29

# Convert from meters to millimeters (mm/day)
precip_full <- precip_full * 1000

# Now convert from mm/day to monthly totals using days per month

log_event("Converting precipitation from mm/day to mm/month...")
for (i in seq_len(nlyr(precip_full))) {
  precip_full[[i]] <- precip_full[[i]] * days_in_month[i]
}

# For PET, since ti is in mm/day monthly totals (mm/month) do the same:
log_event("Converting PET from mm/day to mm/month...")
for (i in seq_len(nlyr(pet_full))) {
  pet_full[[i]] <- pet_full[[i]] * days_in_month[i]
}
cat("\n=== Pr AFTER × 1000 * N.days CONVERSION (mm/month) ===\n")
cat("Precipitation min:", min(values(precip_full), na.rm=TRUE), "\n")
cat("Precipitation max:", max(values(precip_full), na.rm=TRUE), "\n")
cat("Precipitation mean:", mean(values(precip_full), na.rm=TRUE), "\n")
cat("=== PET AFTER × N.days CONVERSION (mm/month) ===\n")
pet_min    <- min(global(pet_full,    "min", na.rm = TRUE))
pet_max    <- max(global(pet_full,    "max", na.rm = TRUE))
pet_mean    <- mean(global(pet_full,    "max", na.rm = TRUE))
cat("PET – min:", pet_min, " max:", pet_max, "\n")
cat("PET Mean:",pet_mean, "\n")

original_template <- rast(precip_full, nlyrs = 1)
saveRDS(original_template, file.path(out_dir, "original_template.rds"))
log_event(sprintf("✓ Saved ORIGINAL raster template (full extent: %d x %d cells, res: %.1f m)",
                  nrow(original_template), ncol(original_template), res(original_template)[1]))

log_event("Reprojecting to BC Albers (EPSG:3005)...")
precip_full <- project(precip_full, "EPSG:3005", method = "bilinear")
pet_full    <- project(pet_full, "EPSG:3005", method = "bilinear")

log_event("Clipping to Nechako Basin extent for processing...")
basin_extent <- ext(basin_boundary)
precip_clipped <- crop(precip_full, basin_extent)
pet_clipped    <- crop(pet_full, basin_extent)
precip_clipped <- mask(precip_clipped, vect(basin_boundary),touches = TRUE)
pet_clipped    <- mask(pet_clipped, vect(basin_boundary),touches = TRUE)

n_bbox_cells <- ncell(precip_clipped)
first_layer_vals <- values(precip_clipped[[1]], mat = FALSE)
valid_mask <- !is.na(first_layer_vals)
n_basin_pixels <- sum(valid_mask)
reduction_pct <- 100 * (1 - n_basin_pixels / n_bbox_cells)
log_event(sprintf("✓ BASIN CLIPPING: %d bbox cells → %d basin pixels (%.1f%% reduction)",
                  n_bbox_cells, n_basin_pixels, reduction_pct))
if (n_basin_pixels == 0) stop("CRITICAL ERROR: No valid cells after basin clipping. Check CRS alignment.")

log_event("Extracting time dimension...")
dates <- NULL
tryCatch({
  dates_temp <- terra::time(precip_clipped)
  if (!is.null(dates_temp) && length(dates_temp) > 0) {
    dates <- dates_temp
    log_event("Time extracted using terra::time()")
  }
}, error = function(e) {})
if (is.null(dates) || length(dates) == 0 || all(is.na(dates))) {
  n_layers <- nlyr(precip_clipped)
  start_year <- 1980
  dates <- seq(as.Date(paste0(start_year, "-01-01")), by = "month", length.out = n_layers)
  log_event(paste("Generated", n_layers, "monthly dates from", start_year))
}
years  <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))

# After extracting dates, years, months
log_event("Computing days per month for each layer...")
month_days_base <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Leap year test function
is_leap <- function(year) {
  (year %% 4 == 0) & (year %% 100 != 0 | year %% 400 == 0)
}

# Create a vector of days for each layer
days_in_month <- month_days_base[months]
leap_years <- is_leap(years)
days_in_month[leap_years & months == 2] <- 29   # adjust February in leap years

# ===== BASIN AVERAGE TIME SERIES COMPUTATION =====
# log_event("Converting units from m to mm...")
# precip_monthly_avg <- precip_monthly_avg * 1000
# pet_monthly_avg <- pet_monthly_avg * 1000
log_event("Computing basin-averaged monthly time series...")

precip_monthly_avg <- as.vector(global(precip_clipped, fun = "mean", na.rm = TRUE)[, 1])
pet_monthly_avg <- as.vector(global(pet_clipped, fun = "mean", na.rm = TRUE)[, 1])


if (length(precip_monthly_avg) != length(dates)) stop("Length mismatch in basin average Pr")
if (length(pet_monthly_avg) != length(dates)) stop("Length mismatch in basin average PET")

log_event(sprintf("  Basin-averaged monthly Pr: %d values, range: %.1f to %.1f mm/month", 
                  length(precip_monthly_avg), 
                  min(precip_monthly_avg, na.rm = TRUE), 
                  max(precip_monthly_avg, na.rm = TRUE)))
log_event(sprintf("  Basin-averaged monthly PET: %d values, range: %.1f to %.1f mm/month", 
                  length(pet_monthly_avg), 
                  min(pet_monthly_avg, na.rm = TRUE), 
                  max(pet_monthly_avg, na.rm = TRUE)))

# Annual aggregation for basin averages
log_event("Aggregating basin averages to annual...")
precip_annual_avg_matrix <- aggregate_to_annual_fast(matrix(precip_monthly_avg, ncol = 1), years, method = "sum")
pet_annual_avg_matrix <- aggregate_to_annual_fast(matrix(pet_monthly_avg, ncol = 1), years, method = "mean")
precip_annual_avg <- precip_annual_avg_matrix[, 1]
pet_annual_avg <- pet_annual_avg_matrix[, 1]
annual_years <- as.integer(rownames(precip_annual_avg_matrix))

log_event(sprintf("  Basin-averaged annual Pr: %d values, range: %.1f to %.1f mm/year", 
                  length(precip_annual_avg), 
                  min(precip_annual_avg, na.rm = TRUE), 
                  max(precip_annual_avg, na.rm = TRUE)))
log_event(sprintf("  Basin-averaged annual PET: %d values, range: %.1f to %.1f mm/year", 
                  length(pet_annual_avg), 
                  min(pet_annual_avg, na.rm = TRUE), 
                  max(pet_annual_avg, na.rm = TRUE)))

# Pre-process PET: Replace zeros
log_event(paste("Replacing PET zeros with", min_positive_value, "mm..."))
pet_vals <- values(pet_clipped)
n_valid_data <- sum(!is.na(pet_vals))
n_zeros_before <- sum(pet_vals == 0, na.rm = TRUE)
pet_vals[pet_vals == 0] <- min_positive_value
values(pet_clipped) <- pet_vals
log_event(paste("  Replaced", n_zeros_before, "zeros (", round(n_zeros_before/n_valid_data*100, 2), "%)"))

# Reshape data to matrices (time × basin_pixels)
log_event("Reshaping data for processing...")
n_time <- nlyr(precip_clipped)
valid_cell_indices <- which(valid_mask)
valid_xy <- xyFromCell(precip_clipped, valid_cell_indices)

precip_vals <- values(precip_clipped)
pet_vals    <- values(pet_clipped)
precip_matrix <- t(precip_vals[valid_mask, , drop = FALSE])
pet_matrix    <- t(pet_vals[valid_mask,    , drop = FALSE])

coords_dt <- data.table(
  space_idx = seq_len(n_basin_pixels),
  x = valid_xy[, 1],
  y = valid_xy[, 2],
  is_basin_average = FALSE  # Flag for pixel data
)
log_event(paste("Data reshaped to", n_time, "×", n_basin_pixels, "matrix"))

month_index_list <- compute_month_index(months)

# ===== PROCESS BASIN AVERAGES =====
log_event("Processing basin-averaged time series through trend tests...")

process_basin_series <- function(ts_monthly, ts_annual, var_name, is_precip) {
  conf_cache_env_basin <- new.env(parent = emptyenv())
  
  # Annual analysis
  res_annual <- mk_tfpw_spectral_for_series(ts_annual, is_precip, alpha, max_tie_percent,
                                            n_sim_spectral, conf_cache_env_basin, 
                                            start_year = min(annual_years))  # NEW: Pass start year
  
  # Monthly analysis (by calendar month)
  monthly_results <- vector("list", 12)
  for (m in 1:12) {
    mi <- month_index_list[[as.character(m)]]
    if (length(mi) == 0 || all(is.na(ts_monthly[mi]))) {
      monthly_results[[m]] <- list(
        vc = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                  n = 0, rho1 = NA_real_, vc_corrected = FALSE, n_ties = 0, percent_ties = 0,
                  n_min_vals = 0, percent_min_vals = 0, tau_b_adjusted = FALSE,
                  filtered = TRUE, reason = "no_data"),
        tf = list(tau = NA_real_, p = NA_real_, sl = NA_real_, S = NA_real_, varS = NA_real_,
                  n = 0, rho1 = NA_real_, tfpw_applied = FALSE, n_ties = 0, percent_ties = 0,
                  n_min_vals = 0, percent_min_vals = 0, tau_b_adjusted = FALSE,
                  filtered = TRUE, reason = "no_data"),
        spec = list(n_peaks = 0L, dominant_period = NA_real_, conf = NA_real_),
        cpt = list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                   n_changepoints = 0L, first_changepoint_year = NA_integer_,
                   mean_before = NA_real_, mean_after = NA_real_, magnitude_shift = NA_real_)
      )
    } else {
      ts_vec <- ts_monthly[mi]
      # Monthly changepoint detection not recommended - see documentation
      monthly_results[[m]] <- mk_tfpw_spectral_for_series(ts_vec, is_precip, alpha, max_tie_percent,
                                                          n_sim_spectral, conf_cache_env_basin,
                                                          start_year = NULL)  # No year conversion for monthly
    }
  }
  
  # Create annual result row
  vc_a <- res_annual$vc
  tf_a <- res_annual$tf
  spec_a <- res_annual$spec
  cpt_a <- res_annual$cpt  # NEW: Extract changepoint results
  annual_dt <- data.table(
    space_idx = -1L,
    x = NA_real_, y = NA_real_,
    variable = var_name,
    period = "annual",
    month = NA_integer_,
    tau_vc = vc_a$tau, p_value_vc = vc_a$p, sl_vc = vc_a$sl,
    vc_corrected = vc_a$vc_corrected, n_ties_vc = vc_a$n_ties,
    percent_ties_vc = vc_a$percent_ties, n_min_vals_vc = vc_a$n_min_vals,
    percent_min_vals_vc = vc_a$percent_min_vals, tau_b_adjusted_vc = vc_a$tau_b_adjusted,
    filtered_vc = vc_a$filtered, filter_reason_vc = vc_a$reason,
    tau_tfpw = tf_a$tau, p_value_tfpw = tf_a$p, sl_tfpw = tf_a$sl,
    tfpw_applied = tf_a$tfpw_applied, n_ties_tfpw = tf_a$n_ties,
    percent_ties_tfpw = tf_a$percent_ties, n_min_vals_tfpw = tf_a$n_min_vals,
    percent_min_vals_tfpw = tf_a$percent_min_vals, tau_b_adjusted_tfpw = tf_a$tau_b_adjusted,
    filtered_tfpw = tf_a$filtered, filter_reason_tfpw = tf_a$reason,
    n = vc_a$n, rho1 = vc_a$rho1,
    n_spectral_peaks = spec_a$n_peaks,
    dominant_period = spec_a$dominant_period,
    spectral_confidence = spec_a$conf,
    # NEW: Add changepoint columns
    changepoint_detected = cpt_a$changepoint_detected,
    changepoint_position = cpt_a$changepoint_position,
    n_changepoints = cpt_a$n_changepoints,
    first_changepoint_year = cpt_a$first_changepoint_year,
    mean_before_shift = cpt_a$mean_before,
    mean_after_shift = cpt_a$mean_after,
    magnitude_shift = cpt_a$magnitude_shift,
    same_significance = (vc_a$p < alpha) == (tf_a$p < alpha),
    same_direction = sign(vc_a$tau) == sign(tf_a$tau),
    is_basin_average = TRUE
  )
  
  # Create monthly result rows
  monthly_dts <- list()
  for (m in 1:12) {
    vc_m <- monthly_results[[m]]$vc
    tf_m <- monthly_results[[m]]$tf
    spec_m <- monthly_results[[m]]$spec
    cpt_m <- monthly_results[[m]]$cpt  # NEW: Extract changepoint results
    monthly_dts[[m]] <- data.table(
      space_idx = -1L,
      x = NA_real_, y = NA_real_,
      variable = var_name,
      period = "monthly",
      month = m,
      tau_vc = vc_m$tau, p_value_vc = vc_m$p, sl_vc = vc_m$sl,
      vc_corrected = FALSE,
      n_ties_vc = vc_m$n_ties, percent_ties_vc = vc_m$percent_ties,
      n_min_vals_vc = vc_m$n_min_vals, percent_min_vals_vc = vc_m$percent_min_vals,
      tau_b_adjusted_vc = vc_m$tau_b_adjusted, filtered_vc = vc_m$filtered,
      filter_reason_vc = vc_m$reason,
      tau_tfpw = tf_m$tau, p_value_tfpw = tf_m$p, sl_tfpw = tf_m$sl,
      tfpw_applied = FALSE,
      n_ties_tfpw = tf_m$n_ties, percent_ties_tfpw = tf_m$percent_ties,
      n_min_vals_tfpw = tf_m$n_min_vals, percent_min_vals_tfpw = tf_m$percent_min_vals,
      tau_b_adjusted_tfpw = tf_m$tau_b_adjusted, filtered_tfpw = tf_m$filtered,
      filter_reason_tfpw = tf_m$reason,
      n = NA_integer_, rho1 = NA_real_,
      n_spectral_peaks = spec_m$n_peaks,
      dominant_period = spec_m$dominant_period,
      spectral_confidence = spec_m$conf,
      # NEW: Add changepoint columns (usually NA for monthly - not recommended)
      changepoint_detected = cpt_m$changepoint_detected,
      changepoint_position = cpt_m$changepoint_position,
      n_changepoints = cpt_m$n_changepoints,
      first_changepoint_year = cpt_m$first_changepoint_year,
      mean_before_shift = cpt_m$mean_before,
      mean_after_shift = cpt_m$mean_after,
      magnitude_shift = cpt_m$magnitude_shift,
      same_significance = (vc_m$p < alpha) == (tf_m$p < alpha),
      same_direction = sign(vc_m$tau) == sign(tf_m$tau),
      is_basin_average = TRUE
    )
  }
  monthly_dt <- rbindlist(monthly_dts)
  
  rbindlist(list(annual_dt, monthly_dt))
}

basin_precip_results <- process_basin_series(precip_monthly_avg, precip_annual_avg, "Precipitation", is_precip = TRUE)
basin_pet_results <- process_basin_series(pet_monthly_avg, pet_annual_avg, "PET", is_precip = FALSE)

# ===== MAIN PROCESSING (PIXELS) =====
process_variable_final <- function(data_matrix, var_name, coords_dt, is_precip = FALSE) {
  log_event(paste("Processing", var_name, "..."))
  
  agg_method <- if (is_precip) "sum" else "mean"
  log_event(paste("  Aggregating to annual", agg_method, "..."))
  
  annual_matrix <- aggregate_to_annual_fast(data_matrix, years, method = agg_method)
  
  conf_cache_env <- new.env(parent = emptyenv())
  
  log_event("  Running VC + TFPW MK + Spectral + Changepoint (annual, parallel)...")
  future.seed <- TRUE
  start_year_annual <- min(as.integer(rownames(annual_matrix)))  # NEW: Get first year
  res_annual <- future_lapply(
    seq_len(ncol(annual_matrix)),
    function(i) {
      mk_tfpw_spectral_for_series(
        annual_matrix[, i], is_precip, alpha, max_tie_percent,
        n_sim_spectral, conf_cache_env, start_year = start_year_annual  # NEW: Pass start year
      )
    },
    future.seed = TRUE
  )
  
  unpack_mk <- function(which = c("vc", "tf")) {
    which <- match.arg(which)
    r <- lapply(res_annual, `[[`, which)
    
    DT <- data.table(
      tau      = vapply(r, function(z) z$tau,      numeric(1)),
      p        = vapply(r, function(z) z$p,        numeric(1)),
      sl       = vapply(r, function(z) z$sl,       numeric(1)),
      S        = vapply(r, function(z) z$S,        numeric(1)),
      varS     = vapply(r, function(z) z$varS,     numeric(1)),
      n        = vapply(r, function(z) z$n,        numeric(1)),
      rho1     = vapply(r, function(z) z$rho1,     numeric(1)),
      n_ties   = vapply(r, function(z) z$n_ties,   numeric(1)),
      pct_ties = vapply(r, function(z) z$percent_ties, numeric(1)),
      n_min    = vapply(r, function(z) z$n_min_vals,   numeric(1)),
      pct_min  = vapply(r, function(z) z$percent_min_vals, numeric(1)),
      tau_b_adj= vapply(r, function(z) z$tau_b_adjusted,   logical(1)),
      filtered = vapply(r, function(z) z$filtered,         logical(1)),
      reason   = vapply(r, function(z) z$reason,           character(1))
    )
    
    if (which == "vc") {
      DT[, vc_corrected := vapply(r, function(z) z$vc_corrected, logical(1))]
    } else {
      DT[, tfpw_applied := vapply(r, function(z) z$tfpw_applied, logical(1))]
    }
    DT
  }
  vc_annual   <- unpack_mk("vc")
  tfpw_annual <- unpack_mk("tf")
  spec_annual <- {
    r <- lapply(res_annual, `[[`, "spec")
    data.table(
      n_spectral_peaks = vapply(r, `[[`, integer(1), "n_peaks"),
      dominant_period  = vapply(r, `[[`, numeric(1),  "dominant_period"),
      spectral_confidence = vapply(r, `[[`, numeric(1), "conf")
    )
  }
  
  # NEW: Extract changepoint results
  cpt_annual <- {
    r <- lapply(res_annual, `[[`, "cpt")
    data.table(
      changepoint_detected = vapply(r, `[[`, logical(1), "changepoint_detected"),
      changepoint_position = vapply(r, `[[`, integer(1), "changepoint_position"),
      n_changepoints = vapply(r, `[[`, integer(1), "n_changepoints"),
      first_changepoint_year = vapply(r, `[[`, integer(1), "first_changepoint_year"),
      mean_before_shift = vapply(r, `[[`, numeric(1), "mean_before"),
      mean_after_shift = vapply(r, `[[`, numeric(1), "mean_after"),
      magnitude_shift = vapply(r, `[[`, numeric(1), "magnitude_shift")
    )
  }
  
  annual_results <- cbind(
    coords_dt,
    variable = var_name,
    period = "annual",
    month = NA_integer_,
    data.table(
      tau_vc = vc_annual$tau, p_value_vc = vc_annual$p, sl_vc = vc_annual$sl,
      vc_corrected = vc_annual$vc_corrected, n_ties_vc = vc_annual$n_ties,
      percent_ties_vc = vc_annual$pct_ties, n_min_vals_vc = vc_annual$n_min,
      percent_min_vals_vc = vc_annual$pct_min, tau_b_adjusted_vc = vc_annual$tau_b_adj,
      filtered_vc = vc_annual$filtered, filter_reason_vc = vc_annual$reason
    ),
    data.table(
      tau_tfpw = tfpw_annual$tau, p_value_tfpw = tfpw_annual$p, sl_tfpw = tfpw_annual$sl,
      tfpw_applied = tfpw_annual$tfpw_applied, n_ties_tfpw = tfpw_annual$n_ties,
      percent_ties_tfpw = tfpw_annual$pct_ties, n_min_vals_tfpw = tfpw_annual$n_min,
      percent_min_vals_tfpw = tfpw_annual$pct_min, tau_b_adjusted_tfpw = tfpw_annual$tau_b_adj,
      filtered_tfpw = tfpw_annual$filtered, filter_reason_tfpw = tfpw_annual$reason
    ),
    n = vc_annual$n,
    rho1 = vc_annual$rho1,
    spec_annual,
    cpt_annual,  # NEW: Add changepoint columns
    same_significance = (vc_annual$p < alpha) == (tfpw_annual$p < alpha),
    same_direction    = sign(vc_annual$tau) == sign(tfpw_annual$tau)
  )
  
  log_event("  Processing monthly data (12 calendar months, parallel)...")
  monthly_results_list <- vector("list", 12L)
  for (m in 1:12) {
    mi <- month_index_list[[as.character(m)]]
    monthly_subset <- if (length(mi)) data_matrix[mi, , drop = FALSE] else matrix(NA_real_, 0, ncol(data_matrix))
    
    res_month <- future_lapply(
      seq_len(ncol(monthly_subset)),
      function(i) {
        mk_tfpw_spectral_for_series(
          monthly_subset[, i], is_precip, alpha, max_tie_percent,
          n_sim_spectral, conf_cache_env
        )
      },
      future.seed = TRUE
    )
    
    unpack_month <- function(which = c("vc","tf")) {
      which <- match.arg(which)
      r <- lapply(res_month, `[[`, which)
      data.table(
        tau   = vapply(r, `[[`, numeric(1), "tau"),
        p     = vapply(r, `[[`, numeric(1), "p"),
        sl    = vapply(r, `[[`, numeric(1), "sl"),
        n_ties= vapply(r, `[[`, numeric(1), "n_ties"),
        pct_ties = vapply(r, `[[`, numeric(1), "percent_ties"),
        n_min = vapply(r, `[[`, numeric(1), "n_min_vals"),
        pct_min = vapply(r, `[[`, numeric(1), "percent_min_vals"),
        tau_b_adj = vapply(r, `[[`, logical(1), "tau_b_adjusted"),
        filtered  = vapply(r, `[[`, logical(1), "filtered"),
        reason    = vapply(r, `[[`, character(1), "reason")
      )
    }
    vc_m   <- unpack_month("vc")
    tf_m   <- unpack_month("tf")
    spec_m <- {
      r <- lapply(res_month, `[[`, "spec")
      data.table(
        n_spectral_peaks = vapply(r, `[[`, integer(1), "n_peaks"),
        dominant_period  = vapply(r, `[[`, numeric(1),  "dominant_period"),
        spectral_confidence = vapply(r, `[[`, numeric(1), "conf")
      )
    }
    # NEW: Extract changepoint results for monthly data
    cpt_m <- {
      r <- lapply(res_month, `[[`, "cpt")
      data.table(
        changepoint_detected = vapply(r, `[[`, logical(1), "changepoint_detected"),
        changepoint_position = vapply(r, `[[`, integer(1), "changepoint_position"),
        n_changepoints = vapply(r, `[[`, integer(1), "n_changepoints"),
        first_changepoint_year = vapply(r, `[[`, integer(1), "first_changepoint_year"),
        mean_before_shift = vapply(r, `[[`, numeric(1), "mean_before"),
        mean_after_shift = vapply(r, `[[`, numeric(1), "mean_after"),
        magnitude_shift = vapply(r, `[[`, numeric(1), "magnitude_shift")
      )
    }
    
    monthly_results_list[[m]] <- cbind(
      coords_dt,
      variable = var_name,
      period = "monthly",
      month = m,
      data.table(
        tau_vc = vc_m$tau, p_value_vc = vc_m$p, sl_vc = vc_m$sl,
        vc_corrected = FALSE,
        n_ties_vc = vc_m$n_ties, percent_ties_vc = vc_m$pct_ties,
        n_min_vals_vc = vc_m$n_min, percent_min_vals_vc = vc_m$pct_min,
        tau_b_adjusted_vc = vc_m$tau_b_adj, filtered_vc = vc_m$filtered,
        filter_reason_vc = vc_m$reason
      ),
      data.table(
        tau_tfpw = tf_m$tau, p_value_tfpw = tf_m$p, sl_tfpw = tf_m$sl,
        tfpw_applied = FALSE,
        n_ties_tfpw = tf_m$n_ties, percent_ties_tfpw = tf_m$pct_ties,
        n_min_vals_tfpw = tf_m$n_min, percent_min_vals_tfpw = tf_m$pct_min,
        tau_b_adjusted_tfpw = tf_m$tau_b_adj, filtered_tfpw = tf_m$filtered,
        filter_reason_tfpw = tf_m$reason
      ),
      n = NA_integer_, rho1 = NA_real_,
      spec_m,
      cpt_m,  # NEW: Add changepoint columns
      same_significance = (vc_m$p < alpha) == (tf_m$p < alpha),
      same_direction    = sign(vc_m$tau) == sign(tf_m$tau)
    )
  }
  monthly_results <- rbindlist(monthly_results_list, use.names = TRUE, fill = TRUE)
  
  all_results <- rbindlist(list(annual_results, monthly_results), use.names = TRUE, fill = TRUE)
  all_results
}

log_event("=== STARTING PRECIPITATION ANALYSIS ===")
precip_results <- process_variable_final(precip_matrix, "Precipitation", coords_dt, is_precip = TRUE)

log_event("=== STARTING PET ANALYSIS ===")
pet_results <- process_variable_final(pet_matrix, "PET", coords_dt, is_precip = FALSE)

# ===== COMBINE PIXEL AND BASIN RESULTS =====
log_event("Combining pixel-level and basin-averaged results...")
all_results <- rbindlist(list(precip_results, pet_results, basin_precip_results, basin_pet_results), 
                         use.names = TRUE, fill = TRUE)

# ===== EFFICIENT STORAGE =====
log_event("Saving results in RDS format (base R, no external packages)...")
saveRDS(all_results, file.path(out_dir, "all_results.rds"), compress = "gzip")

# Save comprehensive metadata including raw time series
saveRDS(list(
  precip_results = precip_results,
  pet_results = pet_results,
  basin_precip_results = basin_precip_results,
  basin_pet_results = basin_pet_results,
  basin_pixels = n_basin_pixels,
  bbox_cells = n_bbox_cells,
  reduction_pct = reduction_pct,
  processing_date = Sys.time(),
  parameters = list(
    alpha = alpha,
    max_tie_percent = max_tie_percent,
    max_min_value_pct_precip = max_min_value_pct_precip,
    max_min_value_pct_pet = max_min_value_pct_pet
  ),
  original_extent = list(
    xmin = xmin(original_template),
    xmax = xmax(original_template),
    ymin = ymin(original_template),
    ymax = ymax(original_template),
    nrows = nrow(original_template),
    ncols = ncol(original_template),
    res = res(original_template)
  ),
  # CRITICAL: Store raw basin-averaged time series for visualization
  basin_avg_monthly = data.frame(
    date = dates,
    precip_mm_month = precip_monthly_avg,
    pet_mm_month = pet_monthly_avg
  ),
  basin_avg_annual = data.frame(
    year = annual_years,
    precip_mm_year = precip_annual_avg,
    pet_mm_year = pet_annual_avg
  )
), file.path(out_dir, "analysis_metadata.rds"))

# Save summary statistics
summary_stats <- all_results[, .(
  n_total = .N,
  n_valid_vc = sum(!filtered_vc & !is.na(p_value_vc)),
  n_significant_vc = sum(!filtered_vc & p_value_vc < 0.05, na.rm = TRUE),
  n_valid_tfpw = sum(!filtered_tfpw & !is.na(p_value_tfpw)),
  n_significant_tfpw = sum(!filtered_tfpw & p_value_tfpw < 0.05, na.rm = TRUE)
), by = .(variable, period, month, is_basin_average)]
fwrite(summary_stats, file.path(out_dir, "summary_statistics.csv"))

log_event("==========================================")
log_event("DATA PROCESSING COMPLETE")
log_event("==========================================")
log_event(sprintf("Results saved to: %s", out_dir))
log_event(" - all_results.rds (compressed base R format)")
log_event(" - analysis_metadata.rds (includes basin-averaged time series)")
log_event(" - summary_statistics.csv (quick overview)")
log_event(" - original_template.rds (FULL raster extent for proper outputs)")
log_event(" - basin_boundary.rds (for mapping)")
log_event("==========================================")
log_event("NEW FEATURE: PELT Changepoint Detection")
log_event(" - Applied to ANNUAL precipitation and PET time series")
log_event(" - Detects regime shifts in mean level (e.g., 1976-77 PDO shift)")
log_event(" - New columns: changepoint_detected, first_changepoint_year, magnitude_shift")
log_event(" - Monthly changepoints included but NOT RECOMMENDED (see documentation)")
log_event("==========================================")

cat("\n✓ DATA PROCESSING COMPLETED SUCCESSFULLY\n")
cat(sprintf("  Basin pixels processed: %d (from %d bbox cells)\n", n_basin_pixels, n_bbox_cells))
cat(sprintf("  Basin-averaged time series computed and analyzed\n"))
cat(sprintf("  PELT changepoint detection applied to annual data\n"))
cat(sprintf("  Results stored in: %s\n", out_dir))

plan(sequential)