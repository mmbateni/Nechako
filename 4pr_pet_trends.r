####################################################################################
# Comprehensive Trend Analysis
# Literature-Based Tie Handling + Data Filtering + Kendall Tau-b Adjustment
# Based on: Kendall (1945, 1976), Hirsch et al. (1982), Helsel (2005)
####################################################################################

library(ncdf4)
library(terra)
library(data.table)
library(Kendall)
library(future.apply)
library(zoo)

setwd("D:/Nechako_Drought/Nechako/")

# Input files
precip_path <- "monthly_data_direct/total_precipitation_monthly.nc"
pet_path    <- "monthly_data_direct/potential_evapotranspiration_monthly.nc"

# Output directory
out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Log file
LOG_FILE <- file.path(out_dir, "comprehensive_analysis_final.log")
cat("Comprehensive Trend Analysis - FINAL VERSION\n", file = LOG_FILE)
cat("Literature-Based Tie Corrections Implemented\n", file = LOG_FILE, append = TRUE)
cat(paste("Analysis started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat("==========================================\n", file = LOG_FILE, append = TRUE)

log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, " | ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, " | ", msg))
}

# Parameters
alpha <- 0.05    # significance level
n_sim_spectral <- 500  # Monte Carlo simulations for spectral analysis
max_tie_percent <- 50  # Threshold for filtering problematic grid points
min_positive_value <- 0.01  # Replacement value for zeros in PET data

# Set up parallel processing
num_cores <- min(parallel::detectCores() - 1, 8)
if (num_cores < 1) num_cores <- 1
plan(multisession, workers = num_cores)
log_event(paste("Using", num_cores, "cores for parallel processing"))

# ---- KENDALL TAU-B VARIANCE ADJUSTMENT FOR TIES ----
# Based on Kendall (1945, 1976) and implemented following literature standards

calculate_variance_with_ties <- function(S, n, x) {
  # Calculate variance of S accounting for ties in x
  # Based on Kendall (1976) formula for tied ranks
  
  # Count tied groups
  tie_table <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  
  if (length(tie_counts) == 0) {
    # No ties - use standard variance formula
    var_s <- n * (n - 1) * (2 * n + 5) / 18
  } else {
    # Ties present - use adjusted variance (Kendall's tau-b)
    # Var(S) = [n(n-1)(2n+5) - Σt(t-1)(2t+5)] / 18
    n0 <- n * (n - 1) / 2  # Total number of pairs
    
    # Base variance
    var_s <- n * (n - 1) * (2 * n + 5) / 18
    
    # Subtract tie adjustment
    tie_adjustment <- sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
    var_s <- var_s - tie_adjustment
  }
  
  return(var_s)
}

# ---- MODIFIED MANN-KENDALL WITH KENDALL TAU-B VARIANCE ----

# ==============================================================================
# FIXED VERSION OF modified_mann_kendall_taub FUNCTION
# Replace lines 49-194 in your original script with this
# ==============================================================================

# ---- MANUAL SEN'S SLOPE CALCULATION ----
calculate_sens_slope_manual <- function(x) {
  # Calculate Sen's slope: median of all pairwise slopes
  # More robust than sens.slope() from Kendall package
  n <- length(x)
  if (n < 2) return(NA)
  
  slopes <- numeric()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      slope <- (x[j] - x[i]) / (j - i)
      slopes <- c(slopes, slope)
    }
  }
  
  return(median(slopes, na.rm = TRUE))
}


# ---- KENDALL TAU-B VARIANCE ADJUSTMENT FOR TIES ----
calculate_variance_with_ties <- function(S, n, x) {
  # Calculate variance of S accounting for ties in x
  # Based on Kendall (1976) formula for tied ranks
  
  tie_table <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  
  if (length(tie_counts) == 0) {
    # No ties - use standard variance formula
    var_s <- n * (n - 1) * (2 * n + 5) / 18
  } else {
    # Ties present - use adjusted variance (Kendall's tau-b)
    var_s <- n * (n - 1) * (2 * n + 5) / 18
    tie_adjustment <- sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
    var_s <- var_s - tie_adjustment
  }
  
  return(var_s)
}

# ---- MODIFIED MANN-KENDALL WITH KENDALL TAU-B VARIANCE ----
# ---- MODIFIED MANN-KENDALL WITH KENDALL TAU-B VARIANCE ----
modified_mann_kendall_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50, 
                                       apply_vc = TRUE) {
  n_time <- nrow(ts_matrix)
  n_space <- ncol(ts_matrix)
  
  results <- data.frame(
    tau = numeric(n_space),
    p.value = numeric(n_space),
    sl = numeric(n_space),
    S = numeric(n_space),
    varS = numeric(n_space),
    n = integer(n_space),
    rho1 = numeric(n_space),
    vc_corrected = logical(n_space),
    n_ties = integer(n_space),
    percent_ties = numeric(n_space),
    tau_b_adjusted = logical(n_space),
    filtered = logical(n_space)
  )
  
  # Counter for diagnostics
  n_filtered_low_n <- 0
  n_filtered_high_ties <- 0
  n_filtered_sens_fail <- 0
  n_success <- 0
  
  for (i in 1:n_space) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    # Filter: insufficient data
    if (n < 10) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, 0, 0, FALSE, TRUE)
      n_filtered_low_n <- n_filtered_low_n + 1
      next
    }
    
    # Calculate tie statistics
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    
    # Filter: excessive ties
    if (percent_ties > max_tie_pct) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      n_filtered_high_ties <- n_filtered_high_ties + 1
      next
    }
    
    # ===== FIXED: Calculate Sen's slope manually =====
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    
    # Only filter if Sen's slope genuinely fails (very rare)
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      n_filtered_sens_fail <- n_filtered_sens_fail + 1
      next
    }
    # ===== END FIX =====
    
    # Calculate Mann-Kendall S statistic manually
    S <- 0
    for (j in 1:(n-1)) {
      for (k in (j+1):n) {
        S <- S + sign(ts_clean[k] - ts_clean[j])
      }
    }
    
    # Calculate variance with Kendall's tau-b adjustment for ties
    varS_taub <- calculate_variance_with_ties(S, n, ts_clean)
    tau_b_adjusted <- (percent_ties > 5)
    
    # Calculate Kendall's tau
    n_pairs <- n * (n - 1) / 2
    tau <- S / n_pairs
    
    # Apply Variance Correction for autocorrelation if requested
    vc_corrected <- FALSE
    varS_final <- varS_taub
    rho1 <- NA
    
    if (apply_vc) {
      acf_result <- tryCatch({
        acf(ts_clean, lag.max = 1, plot = FALSE, na.action = na.pass)
      }, error = function(e) NULL, warning = function(w) NULL)
      
      rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA
      
      if (!is.na(rho1) && abs(rho1) > 0.1) {
        # Hamed and Rao (1998) variance correction
        correction_factor <- 1 + (2 * rho1 * (n - 1 - 2 * (n - 1) * rho1 + 3 * rho1 * rho1)) / 
          ((n - 1) * (1 - rho1) * (1 - rho1))
        
        varS_final <- varS_taub * correction_factor
        vc_corrected <- TRUE
      }
    }
    
    # Calculate p-value
    if (varS_final <= 0) {
      p_value <- NA
    } else {
      z_stat <- S / sqrt(varS_final)
      p_value <- 2 * pnorm(-abs(z_stat))
    }
    
    results[i, ] <- list(
      tau = tau,
      p.value = p_value,
      sl = sen_slope,  # FIXED: using manual calculation
      S = S,
      varS = varS_final,
      n = n,
      rho1 = rho1,
      vc_corrected = vc_corrected,
      n_ties = n_ties,
      percent_ties = percent_ties,
      tau_b_adjusted = tau_b_adjusted,
      filtered = FALSE  # Not filtered!
    )
    
    n_success <- n_success + 1
  }
  
  # Print diagnostic summary
  cat("\n=== VC Mann-Kendall Summary ===\n")
  cat("Total points:", n_space, "\n")
  cat("Successfully processed:", n_success, 
      sprintf("(%.1f%%)\n", 100 * n_success / n_space))
  cat("Filtered - low n (<10):", n_filtered_low_n, "\n")
  cat("Filtered - high ties (>50%):", n_filtered_high_ties, "\n")
  cat("Filtered - Sen's slope fail:", n_filtered_sens_fail, "\n")
  cat("===============================\n\n")
  
  return(results)
}

# ==============================================================================

# ---- TFPW MANN-KENDALL WITH TAU-B ----

# ==============================================================================
# FIXED VERSION OF perform_tfpw_mk_taub FUNCTION
# Replace lines 196-316 in your original script with this
# ==============================================================================

perform_tfpw_mk_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50) {
  n_time <- nrow(ts_matrix)
  n_space <- ncol(ts_matrix)
  
  results <- data.frame(
    tau = numeric(n_space),
    p.value = numeric(n_space),
    sl = numeric(n_space),
    S = numeric(n_space),
    varS = numeric(n_space),
    n = integer(n_space),
    rho1 = numeric(n_space),
    tfpw_applied = logical(n_space),
    n_ties = integer(n_space),
    percent_ties = numeric(n_space),
    tau_b_adjusted = logical(n_space),
    filtered = logical(n_space)
  )
  
  # Counter for diagnostics
  n_filtered_low_n <- 0
  n_filtered_high_ties <- 0
  n_filtered_sens_fail <- 0
  n_success <- 0
  
  for (i in 1:n_space) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    # Filter: insufficient data
    if (n < 10) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, 0, 0, FALSE, TRUE)
      n_filtered_low_n <- n_filtered_low_n + 1
      next
    }
    
    # Calculate tie statistics
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    
    # Filter: excessive ties
    if (percent_ties > max_tie_pct) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      n_filtered_high_ties <- n_filtered_high_ties + 1
      next
    }
    
    # ===== FIXED: Step 1 - Calculate Sen's slope manually =====
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      n_filtered_sens_fail <- n_filtered_sens_fail + 1
      next
    }
    # ===== END FIX =====
    
    # Step 2: Detrend
    time_index <- 1:n
    trend_line <- sen_slope * time_index
    detrended_series <- ts_clean - trend_line
    
    # Step 3: Calculate autocorrelation on detrended series
    acf_result <- tryCatch({
      acf(detrended_series, lag.max = 1, plot = FALSE, na.action = na.pass)
    }, error = function(e) NULL, warning = function(w) NULL)
    rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA
    
    # Step 4: Pre-whitening if needed
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
    
    # Step 5: Calculate Mann-Kendall with tau-b variance
    S <- 0
    for (j in 1:(n-1)) {
      for (k in (j+1):n) {
        S <- S + sign(corrected_series[k] - corrected_series[j])
      }
    }
    
    # Calculate variance with tie adjustment
    varS_taub <- calculate_variance_with_ties(S, n, corrected_series)
    tau_b_adjusted <- (percent_ties > 5)
    
    # Calculate tau
    n_pairs <- n * (n - 1) / 2
    tau <- S / n_pairs
    
    # Calculate p-value
    if (varS_taub <= 0) {
      p_value <- NA
    } else {
      z_stat <- S / sqrt(varS_taub)
      p_value <- 2 * pnorm(-abs(z_stat))
    }
    
    results[i, ] <- list(
      tau = tau,
      p.value = p_value,
      sl = sen_slope,  # FIXED: using manual calculation
      S = S,
      varS = varS_taub,
      n = n,
      rho1 = rho1,
      tfpw_applied = tfpw_applied,
      n_ties = n_ties,
      percent_ties = percent_ties,
      tau_b_adjusted = tau_b_adjusted,
      filtered = FALSE  # Not filtered!
    )
    
    n_success <- n_success + 1
  }
  
  # Print diagnostic summary
  cat("\n=== TFPW Mann-Kendall Summary ===\n")
  cat("Total points:", n_space, "\n")
  cat("Successfully processed:", n_success, 
      sprintf("(%.1f%%)\n", 100 * n_success / n_space))
  cat("Filtered - low n (<10):", n_filtered_low_n, "\n")
  cat("Filtered - high ties (>50%):", n_filtered_high_ties, "\n")
  cat("Filtered - Sen's slope fail:", n_filtered_sens_fail, "\n")
  cat("=================================\n\n")
  
  return(results)
}


# ---- SPECTRAL ANALYSIS  ----

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

# ---- LOAD AND PREPARE DATA ----

log_event("Loading precipitation and PET data...")
precip <- rast(precip_path)
pet <- rast(pet_path)

# AFTER loading precip/pet rasters, BEFORE any calculations:
target_crs <- "EPSG:3005"  # BC Albers Equal Area

# Reproject rasters to BC Albers
if (!same.crs(precip, target_crs)) {
  cat("✓ Reprojecting to BC Albers (EPSG:3005) for area-accurate calculations...\n")
  precip <- project(precip, target_crs, method = "bilinear")
  pet <- project(pet, target_crs, method = "bilinear")
}

# # Reproject basin boundary to match
# if (!same.crs(basin, target_crs)) {
#   basin <- project(basin, target_crs)
# }
# Get time dimension
log_event("Extracting time dimension...")
dates <- NULL
tryCatch({
  dates_temp <- terra::time(precip)
  if (!is.null(dates_temp) && length(dates_temp) > 0) {
    dates <- dates_temp
    log_event("Time extracted using terra::time()")
  }
}, error = function(e) {})

if (is.null(dates) || length(dates) == 0 || all(is.na(dates))) {
  n_layers <- nlyr(precip)
  start_year <- 1980
  dates <- seq(as.Date(paste0(start_year, "-01-01")), by = "month", length.out = n_layers)
  log_event(paste("Generated", n_layers, "monthly dates starting from", start_year))
}

years <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))

lon_coords <- xFromCol(precip, 1:ncol(precip))
lat_coords <- yFromRow(precip, 1:nrow(precip))

unique_years <- unique(years)
n_years <- length(unique_years)

log_event(paste("Processing", n_years, "years from", min(unique_years), "to", max(unique_years)))
log_event(paste("Grid dimensions:", nrow(precip), "x", ncol(precip), "=", nrow(precip) * ncol(precip), "cells"))

# Convert units (m to mm)
log_event("Converting units from m to mm...")
precip <- precip * 1000
pet <- pet * 1000

# PRE-PROCESS PET: Replace zeros with minimum positive value
log_event(paste("Pre-processing PET: replacing zeros with", min_positive_value, "mm..."))
pet_vals <- values(pet)
n_zeros_before <- sum(pet_vals == 0, na.rm = TRUE)
pet_vals[pet_vals == 0] <- min_positive_value
values(pet) <- pet_vals
n_zeros_after <- sum(values(pet) == 0, na.rm = TRUE)
log_event(paste("  Replaced", n_zeros_before, "zeros (", 
                round(n_zeros_before/length(pet_vals)*100, 2), "% of data)"))

# Reshape data
log_event("Reshaping data for vectorized processing...")
n_rows <- nrow(precip)
n_cols <- ncol(precip)
n_time <- nlyr(precip)
n_space <- n_rows * n_cols

precip_matrix <- matrix(values(precip), nrow = n_time, ncol = n_space, byrow = FALSE)
pet_matrix <- matrix(values(pet), nrow = n_time, ncol = n_space, byrow = FALSE)

coords_dt <- data.table(
  space_idx = 1:n_space,
  lat = rep(lat_coords, each = n_cols),
  lon = rep(lon_coords, times = n_rows)
)

log_event(paste("Data reshaped to", n_time, "x", n_space, "matrix"))

# ---- AGGREGATION FUNCTIONS ----

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

# ---- MAIN PROCESSING FUNCTION ----

process_variable_final <- function(data_matrix, var_name, coords_dt) {
  log_event(paste("Processing", var_name, "with literature-based methods..."))
  
  agg_method <- if (var_name == "Precipitation") "sum" else "mean"
  
  # Annual processing
  log_event(paste("  Aggregating to annual", agg_method, "..."))
  annual_matrix <- aggregate_to_annual(data_matrix, years, method = agg_method)
  
  log_event("  Running VC Mann-Kendall with Kendall tau-b variance adjustment...")
  vc_annual <- modified_mann_kendall_taub(annual_matrix, alpha, max_tie_percent, apply_vc = TRUE)
  
  log_event("  Running TFPW Mann-Kendall with Kendall tau-b variance adjustment...")
  tfpw_annual <- perform_tfpw_mk_taub(annual_matrix, alpha, max_tie_percent)
  
  log_event("  Running spectral analysis...")
  spectral_annual <- perform_spectral_analysis_vectorized(annual_matrix, n_sim_spectral)
  
  spectral_df <- data.table(
    n_spectral_peaks = sapply(spectral_annual, function(x) x$n_peaks),
    dominant_period = sapply(spectral_annual, function(x) x$dominant_period),
    spectral_confidence = sapply(spectral_annual, function(x) x$confidence_limit)
  )
  
  # Combine annual results
  annual_results <- cbind(
    coords_dt,
    variable = var_name,
    period = "annual",
    month = NA,
    setDT(vc_annual)[, .(tau_vc = tau, p_value_vc = p.value, sl_vc = sl, 
                         vc_corrected = vc_corrected, n_ties_vc = n_ties, 
                         percent_ties_vc = percent_ties, tau_b_adjusted_vc = tau_b_adjusted,
                         filtered_vc = filtered)],
    setDT(tfpw_annual)[, .(tau_tfpw = tau, p_value_tfpw = p.value, sl_tfpw = sl, 
                           tfpw_applied = tfpw_applied, n_ties_tfpw = n_ties, 
                           percent_ties_tfpw = percent_ties, tau_b_adjusted_tfpw = tau_b_adjusted,
                           filtered_tfpw = filtered)],
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
    log_event(paste("    Processing month", m, "..."))
    monthly_subset <- extract_monthly_subset(data_matrix, months, m)
    
    vc_monthly <- modified_mann_kendall_taub(monthly_subset, alpha, max_tie_percent, apply_vc = TRUE)
    tfpw_monthly <- perform_tfpw_mk_taub(monthly_subset, alpha, max_tie_percent)
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
                            percent_ties_vc = percent_ties, tau_b_adjusted_vc = tau_b_adjusted,
                            filtered_vc = filtered)],
      setDT(tfpw_monthly)[, .(tau_tfpw = tau, p_value_tfpw = p.value, sl_tfpw = sl, 
                              tfpw_applied = tfpw_applied, n_ties_tfpw = n_ties, 
                              percent_ties_tfpw = percent_ties, tau_b_adjusted_tfpw = tau_b_adjusted,
                              filtered_tfpw = filtered)],
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

# ---- PROCESS BOTH VARIABLES ----

precip_results <- process_variable_final(precip_matrix, "Precipitation", coords_dt)
precip_file <- file.path(out_dir, "precipitation_final_results.csv")
fwrite(precip_results, precip_file)
log_event(paste("Precipitation results saved:", precip_file))

pet_results <- process_variable_final(pet_matrix, "PET", coords_dt)
pet_file <- file.path(out_dir, "pet_final_results.csv")
fwrite(pet_results, pet_file)
log_event(paste("PET results saved:", pet_file))

# Combine results
all_results <- rbindlist(list(precip_results, pet_results))
combined_file <- file.path(out_dir, "combined_final_results.csv")
fwrite(all_results, combined_file)
log_event(paste("Combined results saved:", combined_file))

# ---- SUMMARY STATISTICS ----

log_event("Creating comprehensive summary statistics...")

comprehensive_summary <- all_results[, .(
  total_points = .N,
  
  # VC method statistics
  valid_vc = sum(!is.na(tau_vc)),
  filtered_vc = sum(filtered_vc, na.rm = TRUE),
  significant_vc = sum(p_value_vc < alpha, na.rm = TRUE),
  percent_significant_vc = mean(p_value_vc < alpha, na.rm = TRUE) * 100,
  mean_tau_vc = mean(tau_vc, na.rm = TRUE),
  mean_sl_vc = mean(sl_vc, na.rm = TRUE),
  percent_vc_corrected = mean(vc_corrected, na.rm = TRUE) * 100,
  percent_tau_b_adjusted_vc = mean(tau_b_adjusted_vc, na.rm = TRUE) * 100,
  mean_ties_vc = mean(percent_ties_vc, na.rm = TRUE),
  
  # TFPW method statistics
  valid_tfpw = sum(!is.na(tau_tfpw)),
  filtered_tfpw = sum(filtered_tfpw, na.rm = TRUE),
  significant_tfpw = sum(p_value_tfpw < alpha, na.rm = TRUE),
  percent_significant_tfpw = mean(p_value_tfpw < alpha, na.rm = TRUE) * 100,
  mean_tau_tfpw = mean(tau_tfpw, na.rm = TRUE),
  mean_sl_tfpw = mean(sl_tfpw, na.rm = TRUE),
  percent_tfpw_applied = mean(tfpw_applied, na.rm = TRUE) * 100,
  percent_tau_b_adjusted_tfpw = mean(tau_b_adjusted_tfpw, na.rm = TRUE) * 100,
  mean_ties_tfpw = mean(percent_ties_tfpw, na.rm = TRUE),
  
  # Agreement between methods
  agreement_significance = mean(same_significance, na.rm = TRUE) * 100,
  agreement_direction = mean(same_direction, na.rm = TRUE) * 100,
  
  # Spectral analysis
  points_with_peaks = sum(n_spectral_peaks > 0, na.rm = TRUE),
  percent_with_peaks = mean(n_spectral_peaks > 0, na.rm = TRUE) * 100,
  mean_dominant_period = mean(dominant_period, na.rm = TRUE),
  
  # Autocorrelation
  mean_rho1 = mean(rho1, na.rm = TRUE)
  
), by = .(variable, period, month)]

summary_file <- file.path(out_dir, "summary_statistics_final.csv")
fwrite(comprehensive_summary, summary_file)
log_event(paste("Summary statistics saved:", summary_file))

# ---- CREATE SPATIAL RASTERS ----

log_event("Creating spatial rasters...")

for (var in c("Precipitation", "PET")) {
  for (method in c("vc", "tfpw")) {
    
    p_col <- paste0("p_value_", method)
    tau_col <- paste0("tau_", method)
    filtered_col <- paste0("filtered_", method)
    
    # Annual
    annual_data <- all_results[variable == var & period == "annual"]
    annual_data[, sig_dir := fifelse(get(p_col) < 0.05 & !get(filtered_col), 
                                     sign(get(tau_col)), 0)]
    
    r_sig <- rast(precip, nlyrs = 1)
    values(r_sig) <- NA
    xy <- cbind(annual_data$lon, annual_data$lat)
    cell_indices <- cellFromXY(r_sig, xy)
    values(r_sig)[cell_indices] <- annual_data$sig_dir
    
    writeRaster(r_sig, 
                file.path(out_dir, sprintf("%s_annual_%s_significant.tif", 
                                           tolower(var), method)), 
                overwrite = TRUE)
    
    # Monthly
    for (m in 1:12) {
      monthly_data <- all_results[variable == var & period == "monthly" & month == m]
      monthly_data[, sig_dir := fifelse(get(p_col) < 0.05 & !get(filtered_col), 
                                        sign(get(tau_col)), 0)]
      
      r_sig_m <- rast(precip, nlyrs = 1)
      values(r_sig_m) <- NA
      xy_m <- cbind(monthly_data$lon, monthly_data$lat)
      cell_indices_m <- cellFromXY(r_sig_m, xy_m)
      values(r_sig_m)[cell_indices_m] <- monthly_data$sig_dir
      
      writeRaster(r_sig_m, 
                  file.path(out_dir, sprintf("%s_month%02d_%s_significant.tif", 
                                             tolower(var), m, method)), 
                  overwrite = TRUE)
    }
  }
}

log_event("Raster creation complete!")

# ---- PDF PLOTS ----

log_event("Creating PDF plots...")

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

sym_limits <- function(x, probs = c(0.02, 0.98)) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(-1, 1))
  q <- quantile(x, probs = probs, na.rm = TRUE)
  m <- max(abs(q), na.rm = TRUE)
  if (!is.finite(m) || m == 0) m <- max(abs(x), na.rm = TRUE)
  if (!is.finite(m) || m == 0) m <- 1
  c(-m, m)
}

create_comprehensive_pdf <- function(results_dt, var_name, out_pdf) {
  
  lim_sl_vc <- sym_limits(results_dt[variable == var_name & !filtered_vc]$sl_vc)
  lim_sl_tfpw <- sym_limits(results_dt[variable == var_name & !filtered_tfpw]$sl_tfpw)
  lim_tau_vc <- sym_limits(results_dt[variable == var_name & !filtered_vc]$tau_vc)
  lim_tau_tfpw <- sym_limits(results_dt[variable == var_name & !filtered_tfpw]$tau_tfpw)
  
  col_div <- hcl.colors(101, "RdBu", rev = TRUE)
  
  pdf(out_pdf, width = 14, height = 10, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  
  plot_comparison_page <- function(period_name, month_num = NA) {
    
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 6), oma = c(0, 0, 2, 0))
    
    if (period_name == "annual") {
      subset_data <- results_dt[variable == var_name & period == period_name]
      page_title <- paste(var_name, "- Annual Trends (VC vs TFPW) - Literature-Based")
    } else {
      subset_data <- results_dt[variable == var_name & period == period_name & month == month_num]
      page_title <- paste0(var_name, " - ", month_names[month_num], " Trends (VC vs TFPW)")
    }
    
    if (nrow(subset_data) == 0) {
      plot.new()
      title(paste(page_title, "\n(No data)"))
      return(invisible(NULL))
    }
    
    # Create rasters
    r_sl_vc <- rast(precip, nlyrs = 1); values(r_sl_vc) <- NA
    r_tau_vc <- rast(precip, nlyrs = 1); values(r_tau_vc) <- NA
    r_sig_vc <- rast(precip, nlyrs = 1); values(r_sig_vc) <- NA
    
    r_sl_tfpw <- rast(precip, nlyrs = 1); values(r_sl_tfpw) <- NA
    r_tau_tfpw <- rast(precip, nlyrs = 1); values(r_tau_tfpw) <- NA
    r_sig_tfpw <- rast(precip, nlyrs = 1); values(r_sig_tfpw) <- NA
    
    xy <- cbind(subset_data$lon, subset_data$lat)
    cells <- cellFromXY(r_sl_vc, xy)
    
    values(r_sl_vc)[cells] <- subset_data$sl_vc
    values(r_tau_vc)[cells] <- subset_data$tau_vc
    values(r_sig_vc)[cells] <- fifelse(subset_data$p_value_vc < 0.05 & !subset_data$filtered_vc, 
                                       sign(subset_data$tau_vc), 0)
    
    values(r_sl_tfpw)[cells] <- subset_data$sl_tfpw
    values(r_tau_tfpw)[cells] <- subset_data$tau_tfpw
    values(r_sig_tfpw)[cells] <- fifelse(subset_data$p_value_tfpw < 0.05 & !subset_data$filtered_tfpw, 
                                         sign(subset_data$tau_tfpw), 0)
    
    # Plot VC
    plot(r_sl_vc, main = "VC: Sen slope", col = col_div, range = lim_sl_vc)
    plot(r_tau_vc, main = "VC: Kendall tau", col = col_div, range = lim_tau_vc)
    plot(r_sig_vc, main = "VC: Significant (p<0.05)", 
         col = c("blue", "grey85", "red"), breaks = c(-1.5, -0.5, 0.5, 1.5))
    
    # Plot TFPW
    plot(r_sl_tfpw, main = "TFPW: Sen slope", col = col_div, range = lim_sl_tfpw)
    plot(r_tau_tfpw, main = "TFPW: Kendall tau", col = col_div, range = lim_tau_tfpw)
    plot(r_sig_tfpw, main = "TFPW: Significant (p<0.05)", 
         col = c("blue", "grey85", "red"), breaks = c(-1.5, -0.5, 0.5, 1.5))
    
    mtext(page_title, outer = TRUE, cex = 1.2, font = 2)
  }
  
  # Annual
  plot_comparison_page("annual")
  
  # Monthly
  for (m in 1:12) {
    plot_comparison_page("monthly", m)
  }
  
  invisible(TRUE)
}

precip_pdf <- file.path(out_dir, "precip_final_maps.pdf")
pet_pdf <- file.path(out_dir, "pet_final_maps.pdf")

create_comprehensive_pdf(precip_results, "Precipitation", precip_pdf)
log_event(paste("Precipitation PDF saved:", precip_pdf))

create_comprehensive_pdf(pet_results, "PET", pet_pdf)
log_event(paste("PET PDF saved:", pet_pdf))

# ---- CLEANUP AND FINAL REPORT ----

plan(sequential)

log_event("==========================================")
log_event("COMPREHENSIVE ANALYSIS COMPLETE - FINAL VERSION")
log_event("==========================================")
log_event("METHODS IMPLEMENTED:")
log_event("1. Kendall's tau-b variance adjustment for ties (Kendall 1945, 1976)")
log_event("2. Hamed-Rao variance correction for autocorrelation (1998)")
log_event("3. Trend-Free Pre-Whitening (Yue et al. 2002)")
log_event("4. Grid point filtering (>50% ties excluded)")
log_event("5. PET zero replacement (0 → 0.01 mm)")
log_event("==========================================")
log_event(paste("Files saved in:", out_dir))
log_event(paste("  - Precipitation results:", precip_file))
log_event(paste("  - PET results:", pet_file))
log_event(paste("  - Combined results:", combined_file))
log_event(paste("  - Summary statistics:", summary_file))
log_event("==========================================")

cat("\n✓ FINAL ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("Literature-based methods with tie corrections implemented.\n")
cat("Check", LOG_FILE, "for full details.\n")
cat("\nKey improvements:\n")
cat("  • Kendall's tau-b variance for ties\n")
cat("  • PET zeros replaced with 0.01 mm\n")
cat("  • Grid points with >50% ties filtered\n")
cat("  • No IFAULT=12 errors expected\n")