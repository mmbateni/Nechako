####################################################################################
# Trend Analysis Data Processing
# Purpose: Perform comprehensive trend analysis and save results efficiently
# Output: Binary data files for fast reloading by visualization script
####################################################################################
####################################################################################
library(ncdf4)
library(terra)
library(data.table)
library(Kendall)
library(future.apply)
library(zoo)
library(sf)

setwd("D:/Nechako_Drought/Nechako/")

# ===== LOGGING SETUP =====
out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
LOG_FILE <- file.path(out_dir, "data_processing.log")

log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, " | ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, " | ", msg))
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
      target_crs <- "EPSG:3005"  # BC Albers
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
max_min_value_pct_pet <- 50
min_positive_value <- 0.01
num_cores <- min(parallel::detectCores() - 1, 8)
if (num_cores < 1) num_cores <- 1
plan(multisession, workers = num_cores)
log_event(paste("Using", num_cores, "cores for parallel processing"))

# ===== HELPER FUNCTIONS =====
check_min_value_threshold <- function(ts_clean, min_val_threshold = 0.01, max_pct = 50, var_name = "Unknown", is_precip = FALSE) {
  n_total <- length(ts_clean)
  n_min_vals <- sum(ts_clean <= min_val_threshold, na.rm = TRUE)
  pct_min_vals <- (n_min_vals / n_total) * 100
  max_allowed <- if (is_precip) max_min_value_pct_precip else max_min_value_pct_pet
  exceeds_threshold <- pct_min_vals > max_allowed
  return(list(exceeds_threshold = exceeds_threshold, pct_min_vals = pct_min_vals))
}

calculate_sens_slope_manual <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  slopes <- numeric()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (!is.na(x[i]) && !is.na(x[j])) {
        slope <- (x[j] - x[i]) / (j - i)
        slopes <- c(slopes, slope)
      }
    }
  }
  if (length(slopes) == 0) return(NA)
  return(median(slopes, na.rm = TRUE))
}

calculate_variance_with_ties <- function(S, n, x) {
  tie_table <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  if (length(tie_counts) == 0) {
    var_s <- n * (n - 1) * (2 * n + 5) / 18
  } else {
    var_s <- n * (n - 1) * (2 * n + 5) / 18
    tie_adjustment <- sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
    var_s <- var_s - tie_adjustment
  }
  return(var_s)
}

# ===== DATA LOADING & PREPROCESSING =====
log_event("Loading precipitation and PET data...")
precip_full <- rast(precip_path)
pet_full <- rast(pet_path)

# CRITICAL FIX: Save ORIGINAL raster template BEFORE clipping (for proper output extent)
original_template <- rast(precip_full, nlyrs = 1)
saveRDS(original_template, file.path(out_dir, "original_template.rds"))
log_event(sprintf("✓ Saved ORIGINAL raster template (full extent: %d x %d cells, res: %.1f m)",
                  nrow(original_template), ncol(original_template), res(original_template)[1]))

# Reproject to BC Albers
log_event("Reprojecting to BC Albers (EPSG:3005)...")
precip_full <- project(precip_full, "EPSG:3005", method = "bilinear")
pet_full <- project(pet_full, "EPSG:3005", method = "bilinear")

# Clip to basin extent FOR PROCESSING ONLY (speed optimization)
log_event("Clipping to Nechako Basin extent for processing...")
basin_extent <- ext(basin_boundary)
precip_clipped <- crop(precip_full, basin_extent)
pet_clipped <- crop(pet_full, basin_extent)
precip_clipped <- mask(precip_clipped, vect(basin_boundary))
pet_clipped <- mask(pet_clipped, vect(basin_boundary))

# Diagnostic: Cell reduction
n_bbox_cells <- ncell(precip_clipped)
first_layer_vals <- values(precip_clipped[[1]], mat = FALSE)
valid_mask <- !is.na(first_layer_vals)
n_basin_pixels <- sum(valid_mask)
reduction_pct <- 100 * (1 - n_basin_pixels / n_bbox_cells)
log_event(sprintf("✓ BASIN CLIPPING: %d bbox cells → %d basin pixels (%.1f%% reduction)",
                  n_bbox_cells, n_basin_pixels, reduction_pct))

if (n_basin_pixels == 0) {
  stop("CRITICAL ERROR: No valid cells after basin clipping. Check CRS alignment.")
}

# Extract time dimension
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
years <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))

# Unit conversion (m → mm)
log_event("Converting units from m to mm...")
precip_clipped <- precip_clipped * 1000
pet_clipped <- pet_clipped * 1000

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
pet_vals <- values(pet_clipped)
precip_matrix <- t(precip_vals[valid_mask, ])
pet_matrix <- t(pet_vals[valid_mask, ])

coords_dt <- data.table(
  space_idx = 1:n_basin_pixels,
  x = valid_xy[, 1],
  y = valid_xy[, 2]
)
log_event(paste("Data reshaped to", n_time, "×", n_basin_pixels, "matrix"))

# ===== AGGREGATION FUNCTIONS =====
aggregate_to_annual <- function(monthly_matrix, years, method = "sum") {
  n_years <- length(unique(years))
  n_space <- ncol(monthly_matrix)
  annual_matrix <- matrix(NA, nrow = n_years, ncol = n_space)
  for (i in 1:n_years) {
    year_indices <- which(years == unique(years)[i])
    if (method == "sum") {
      year_data <- colSums(monthly_matrix[year_indices, , drop = FALSE], na.rm = TRUE)
    } else {
      year_data <- colMeans(monthly_matrix[year_indices, , drop = FALSE], na.rm = TRUE)
    }
    n_valid <- colSums(!is.na(monthly_matrix[year_indices, , drop = FALSE]))
    year_data[n_valid < 6] <- NA
    annual_matrix[i, ] <- year_data
  }
  return(annual_matrix)
}

extract_monthly_subset <- function(monthly_matrix, months, target_month) {
  month_indices <- which(months == target_month)
  return(monthly_matrix[month_indices, , drop = FALSE])
}

# ===== TRENDS ANALYSIS FUNCTIONS (same as original code) =====
modified_mann_kendall_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50,
                                       var_name = "Unknown", is_precip = FALSE) {
  n_time <- nrow(ts_matrix)
  n_basin_pixels <- ncol(ts_matrix)
  results <- data.frame(
    tau = numeric(n_basin_pixels),
    p.value = numeric(n_basin_pixels),
    sl = numeric(n_basin_pixels),
    S = numeric(n_basin_pixels),
    varS = numeric(n_basin_pixels),
    n = integer(n_basin_pixels),
    rho1 = numeric(n_basin_pixels),
    vc_corrected = logical(n_basin_pixels),
    n_ties = integer(n_basin_pixels),
    percent_ties = numeric(n_basin_pixels),
    n_min_vals = integer(n_basin_pixels),
    percent_min_vals = numeric(n_basin_pixels),
    tau_b_adjusted = logical(n_basin_pixels),
    filtered = logical(n_basin_pixels),
    filter_reason = character(n_basin_pixels),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_basin_pixels) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    if (n < 10) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, 0, 0, 0, 0, FALSE, TRUE, "low_n")
      next
    }
    
    min_val_threshold <- if (is_precip) 0.0 else min_positive_value
    min_check <- check_min_value_threshold(ts_clean, min_val_threshold, 
                                           max_pct = if(is_precip) max_min_value_pct_precip else max_min_value_pct_pet,
                                           var_name = var_name, is_precip = is_precip)
    if (min_check$exceeds_threshold) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, 0, 0, 0, min_check$pct_min_vals, 
                           FALSE, TRUE, "excessive_min_vals")
      next
    }
    
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    if (percent_ties > max_tie_pct) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, 0, 0, 
                           FALSE, TRUE, "excessive_ties")
      next
    }
    
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, 0, 0, 
                           FALSE, TRUE, "sens_slope_fail")
      next
    }
    
    S <- 0
    for (j in 1:(n-1)) {
      for (k in (j+1):n) {
        S <- S + sign(ts_clean[k] - ts_clean[j])
      }
    }
    
    varS_taub <- calculate_variance_with_ties(S, n, ts_clean)
    tau_b_adjusted <- (percent_ties > 5)
    n_pairs <- n * (n - 1) / 2
    tau <- S / n_pairs
    
    vc_corrected <- FALSE
    varS_final <- varS_taub
    rho1 <- NA
    acf_result <- tryCatch({
      acf(ts_clean, lag.max = 1, plot = FALSE, na.action = na.pass)
    }, error = function(e) NULL, warning = function(w) NULL)
    rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA
    
    if (!is.na(rho1) && abs(rho1) > 0.1) {
      correction_factor <- 1 + (2 * rho1 * (n - 1 - 2 * (n - 1) * rho1 + 3 * rho1 * rho1)) /
        ((n - 1) * (1 - rho1) * (1 - rho1))
      varS_final <- varS_taub * correction_factor
      vc_corrected <- TRUE
    }
    
    if (varS_final <= 0) {
      p_value <- NA
    } else {
      z_stat <- S / sqrt(varS_final)
      p_value <- 2 * pnorm(-abs(z_stat))
    }
    
    results[i, ] <- list(
      tau = tau, p.value = p_value, sl = sen_slope, S = S, varS = varS_final, n = n,
      rho1 = rho1, vc_corrected = vc_corrected, n_ties = n_ties, percent_ties = percent_ties,
      n_min_vals = 0, percent_min_vals = min_check$pct_min_vals,
      tau_b_adjusted = tau_b_adjusted, filtered = FALSE, filter_reason = "none"
    )
  }
  return(results)
}

perform_tfpw_mk_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50,
                                 var_name = "Unknown", is_precip = FALSE) {
  n_time <- nrow(ts_matrix)
  n_basin_pixels <- ncol(ts_matrix)
  results <- data.frame(
    tau = numeric(n_basin_pixels),
    p.value = numeric(n_basin_pixels),
    sl = numeric(n_basin_pixels),
    S = numeric(n_basin_pixels),
    varS = numeric(n_basin_pixels),
    n = integer(n_basin_pixels),
    rho1 = numeric(n_basin_pixels),
    tfpw_applied = logical(n_basin_pixels),
    n_ties = integer(n_basin_pixels),
    percent_ties = numeric(n_basin_pixels),
    n_min_vals = integer(n_basin_pixels),
    percent_min_vals = numeric(n_basin_pixels),
    tau_b_adjusted = logical(n_basin_pixels),
    filtered = logical(n_basin_pixels),
    filter_reason = character(n_basin_pixels),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_basin_pixels) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    if (n < 10) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, 0, 0, 0, 0, FALSE, TRUE, "low_n")
      next
    }
    
    min_val_threshold <- if (is_precip) 0.0 else min_positive_value
    min_check <- check_min_value_threshold(ts_clean, min_val_threshold, 
                                           max_pct = if(is_precip) max_min_value_pct_precip else max_min_value_pct_pet,
                                           var_name = var_name, is_precip = is_precip)
    if (min_check$exceeds_threshold) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, 0, 0, 0, min_check$pct_min_vals, 
                           FALSE, TRUE, "excessive_min_vals")
      next
    }
    
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    if (percent_ties > max_tie_pct) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, 0, 0, 
                           FALSE, TRUE, "excessive_ties")
      next
    }
    
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, 0, 0, 
                           FALSE, TRUE, "sens_slope_fail")
      next
    }
    
    time_index <- 1:n
    trend_line <- sen_slope * time_index
    detrended_series <- ts_clean - trend_line
    
    acf_result <- tryCatch({
      acf(detrended_series, lag.max = 1, plot = FALSE, na.action = na.pass)
    }, error = function(e) NULL, warning = function(w) NULL)
    rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA
    
    tfpw_applied <- FALSE
    if (!is.na(rho1) && abs(rho1) > 0.1) {
      prewhitened_detrended <- numeric(n)
      prewhitened_detrended[1] <- detrended_series[1]
      for (j in 2:n) {
        prewhitened_detrended[j] <- detrended_series[j] - rho1 * detrended_series[j-1]
      }
      corrected_series <- prewhitened_detrended + trend_line
      tfpw_applied <- TRUE
    } else {
      corrected_series <- ts_clean
    }
    
    S <- 0
    for (j in 1:(n-1)) {
      for (k in (j+1):n) {
        S <- S + sign(corrected_series[k] - corrected_series[j])
      }
    }
    
    varS_taub <- calculate_variance_with_ties(S, n, corrected_series)
    tau_b_adjusted <- (percent_ties > 5)
    n_pairs <- n * (n - 1) / 2
    tau <- S / n_pairs
    
    if (varS_taub <= 0) {
      p_value <- NA
    } else {
      z_stat <- S / sqrt(varS_taub)
      p_value <- 2 * pnorm(-abs(z_stat))
    }
    
    results[i, ] <- list(
      tau = tau, p.value = p_value, sl = sen_slope, S = S, varS = varS_taub, n = n,
      rho1 = rho1, tfpw_applied = tfpw_applied, n_ties = n_ties, percent_ties = percent_ties,
      n_min_vals = 0, percent_min_vals = min_check$pct_min_vals,
      tau_b_adjusted = tau_b_adjusted, filtered = FALSE, filter_reason = "none"
    )
  }
  return(results)
}

perform_spectral_analysis_vectorized <- function(ts_matrix, n_sim = 500, alpha = 0.05) {
  n_time <- nrow(ts_matrix)
  n_space <- ncol(ts_matrix)
  results <- vector("list", n_space)
  
  for (i in 1:n_space) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    if (n < 20) {
      results[[i]] <- list(n_peaks = 0, dominant_period = NA, confidence_limit = NA)
      next
    }
    
    sen_slope <- tryCatch(sens.slope(ts_clean)$estimates, error = function(e) 0)
    detrended_series <- ts_clean - sen_slope * (1:n)
    fft_result <- fft(detrended_series - mean(detrended_series))
    spectral_density <- Mod(fft_result[1:(n/2)])^2 / n
    frequencies <- seq(0, 0.5, length.out = n/2)
    
    max_spectra <- numeric(n_sim)
    for (j in 1:n_sim) {
      random_series <- rnorm(n, mean = mean(detrended_series), sd = sd(detrended_series))
      fft_rand <- fft(random_series - mean(random_series))
      spectral_rand <- Mod(fft_rand[1:(n/2)])^2 / n
      max_spectra[j] <- max(spectral_rand, na.rm = TRUE)
    }
    
    conf_limit <- quantile(max_spectra, 1 - alpha, na.rm = TRUE)
    significant_peaks <- spectral_density > conf_limit
    peak_indices <- which(significant_peaks & !is.na(significant_peaks))
    peak_frequencies <- frequencies[peak_indices]
    peak_periods <- ifelse(peak_frequencies > 0, 1/peak_frequencies, NA)
    
    if (length(peak_periods) > 0) {
      sorted_indices <- order(spectral_density[peak_indices], decreasing = TRUE)
      peak_periods <- peak_periods[sorted_indices]
    }
    
    results[[i]] <- list(
      n_peaks = length(peak_periods),
      dominant_period = if (length(peak_periods) > 0) peak_periods[1] else NA,
      confidence_limit = conf_limit
    )
  }
  return(results)
}

# ===== MAIN PROCESSING FUNCTION =====
process_variable_final <- function(data_matrix, var_name, coords_dt, is_precip = FALSE) {
  log_event(paste("Processing", var_name, "..."))
  agg_method <- if (is_precip) "sum" else "mean"
  
  # Annual processing
  log_event(paste("  Aggregating to annual", agg_method, "..."))
  annual_matrix <- aggregate_to_annual(data_matrix, years, method = agg_method)
  
  log_event("  Running VC Mann-Kendall...")
  vc_annual <- modified_mann_kendall_taub(annual_matrix, alpha, max_tie_percent,
                                          var_name = var_name, is_precip = is_precip)
  
  log_event("  Running TFPW Mann-Kendall...")
  tfpw_annual <- perform_tfpw_mk_taub(annual_matrix, alpha, max_tie_percent,
                                      var_name = var_name, is_precip = is_precip)
  
  log_event("  Running spectral analysis...")
  spectral_annual <- perform_spectral_analysis_vectorized(annual_matrix, n_sim_spectral)
  spectral_df <- data.table(
    n_spectral_peaks = sapply(spectral_annual, function(x) x$n_peaks),
    dominant_period = sapply(spectral_annual, function(x) x$dominant_period),
    spectral_confidence = sapply(spectral_annual, function(x) x$confidence_limit)
  )
  
  annual_results <- cbind(
    coords_dt,
    variable = var_name,
    period = "annual",
    month = NA_integer_,
    setDT(vc_annual)[, .(tau_vc = tau, p_value_vc = p.value, sl_vc = sl,
                         vc_corrected = vc_corrected, n_ties_vc = n_ties,
                         percent_ties_vc = percent_ties, n_min_vals_vc = n_min_vals,
                         percent_min_vals_vc = percent_min_vals,
                         tau_b_adjusted_vc = tau_b_adjusted,
                         filtered_vc = filtered, filter_reason_vc = filter_reason)],
    setDT(tfpw_annual)[, .(tau_tfpw = tau, p_value_tfpw = p.value, sl_tfpw = sl,
                           tfpw_applied = tfpw_applied, n_ties_tfpw = n_ties,
                           percent_ties_tfpw = percent_ties, n_min_vals_tfpw = n_min_vals,
                           percent_min_vals_tfpw = percent_min_vals,
                           tau_b_adjusted_tfpw = tau_b_adjusted,
                           filtered_tfpw = filtered, filter_reason_tfpw = filter_reason)],
    n = vc_annual$n,
    rho1 = vc_annual$rho1,
    spectral_df,
    same_significance = (vc_annual$p.value < alpha) == (tfpw_annual$p.value < alpha),
    same_direction = sign(vc_annual$tau) == sign(tfpw_annual$tau)
  )
  
  # Monthly processing
  log_event("  Processing monthly data (12 calendar months)...")
  monthly_results_list <- vector("list", 12)
  for (m in 1:12) {
    log_event(paste("    Month", m, "..."))
    monthly_subset <- extract_monthly_subset(data_matrix, months, m)
    
    vc_monthly <- modified_mann_kendall_taub(monthly_subset, alpha, max_tie_percent,
                                             var_name = paste(var_name, "month", m),
                                             is_precip = is_precip)
    tfpw_monthly <- perform_tfpw_mk_taub(monthly_subset, alpha, max_tie_percent,
                                         var_name = paste(var_name, "month", m),
                                         is_precip = is_precip)
    spectral_monthly <- perform_spectral_analysis_vectorized(monthly_subset, n_sim_spectral)
    
    spectral_df_m <- data.table(
      n_spectral_peaks = sapply(spectral_monthly, function(x) x$n_peaks),
      dominant_period = sapply(spectral_monthly, function(x) x$dominant_period),
      spectral_confidence = sapply(spectral_monthly, function(x) x$confidence_limit)
    )
    
    monthly_results_list[[m]] <- cbind(
      coords_dt,
      variable = var_name,
      period = "monthly",
      month = m,
      setDT(vc_monthly)[, .(tau_vc = tau, p_value_vc = p.value, sl_vc = sl,
                            vc_corrected = vc_corrected, n_ties_vc = n_ties,
                            percent_ties_vc = percent_ties, n_min_vals_vc = n_min_vals,
                            percent_min_vals_vc = percent_min_vals,
                            tau_b_adjusted_vc = tau_b_adjusted,
                            filtered_vc = filtered, filter_reason_vc = filter_reason)],
      setDT(tfpw_monthly)[, .(tau_tfpw = tau, p_value_tfpw = p.value, sl_tfpw = sl,
                              tfpw_applied = tfpw_applied, n_ties_tfpw = n_ties,
                              percent_ties_tfpw = percent_ties, n_min_vals_tfpw = n_min_vals,
                              percent_min_vals_tfpw = percent_min_vals,
                              tau_b_adjusted_tfpw = tau_b_adjusted,
                              filtered_tfpw = filtered, filter_reason_tfpw = filter_reason)],
      n = vc_monthly$n,
      rho1 = vc_monthly$rho1,
      spectral_df_m,
      same_significance = (vc_monthly$p.value < alpha) == (tfpw_monthly$p.value < alpha),
      same_direction = sign(vc_monthly$tau) == sign(tfpw_monthly$tau)
    )
  }
  
  monthly_results <- rbindlist(monthly_results_list)
  all_results <- rbindlist(list(annual_results, monthly_results))
  return(all_results)
}

# ===== EXECUTE PROCESSING =====
log_event("=== STARTING PRECIPITATION ANALYSIS ===")
precip_results <- process_variable_final(precip_matrix, "Precipitation", coords_dt, is_precip = TRUE)
log_event("=== STARTING PET ANALYSIS ===")
pet_results <- process_variable_final(pet_matrix, "PET", coords_dt, is_precip = FALSE)

# ===== EFFICIENT STORAGE (NO fst REQUIRED) =====
log_event("Saving results in RDS format (base R, no external packages)...")
all_results <- rbindlist(list(precip_results, pet_results))

# Save full results with compression
saveRDS(all_results, file.path(out_dir, "all_results.rds"), compress = "gzip")
saveRDS(list(
  precip_results = precip_results,
  pet_results = pet_results,
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
  )
), file.path(out_dir, "analysis_metadata.rds"))

# Save summary statistics
summary_stats <- all_results[, .(
  n_total = .N,
  n_valid_vc = sum(!filtered_vc & !is.na(p_value_vc)),
  n_significant_vc = sum(!filtered_vc & p_value_vc < 0.05, na.rm = TRUE),
  n_valid_tfpw = sum(!filtered_tfpw & !is.na(p_value_tfpw)),
  n_significant_tfpw = sum(!filtered_tfpw & p_value_tfpw < 0.05, na.rm = TRUE)
), by = .(variable, period, month)]

fwrite(summary_stats, file.path(out_dir, "summary_statistics.csv"))

log_event("==========================================")
log_event("DATA PROCESSING COMPLETE")
log_event("==========================================")
log_event(sprintf("Results saved to: %s", out_dir))
log_event("  - all_results.rds (compressed base R format)")
log_event("  - analysis_metadata.rds (full context)")
log_event("  - summary_statistics.csv (quick overview)")
log_event("  - original_template.rds (FULL raster extent for proper outputs)")
log_event("  - basin_boundary.rds (for mapping)")
log_event("==========================================")
cat("\n✓ DATA PROCESSING COMPLETED SUCCESSFULLY\n")
cat(sprintf("  Basin pixels processed: %d (from %d bbox cells)\n", n_basin_pixels, n_bbox_cells))
cat(sprintf("  Processing time reduction: %.1f%%\n", reduction_pct))
cat(sprintf("  Results stored in: %s\n", out_dir))
plan(sequential)