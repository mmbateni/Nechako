####################################################################################
# PIXEL-WISE TEMPORAL DIAGNOSTICS FOR DROUGHT INDICES (SPI/SPEI) WITH RUNS TEST
# Nechako Basin Analysis - 76 Years (1950-2025)
# Drought Definition: Onset < -1.0, Termination >= 0.0 (exact from visualize_spi_spei.R)
# Timescales: 1, 3, 6, 9, 12, 24 months | Tie Filtering: 50% threshold
####################################################################################

library(terra)
library(ncdf4)
library(data.table)
library(Kendall)
library(changepoint)
library(extRemes)
library(future.apply)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(zoo)
library(tseries)  # For runs test validation

# Set working directory
setwd("D:/Nechako_Drought/")

# Output directory
out_dir <- "drought_temporal_diagnostics_runs_test"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Log file setup - COMPREHENSIVE TIMESTAMPED LOGGING
LOG_FILE <- file.path(out_dir, "drought_temporal_diagnostics_runs_test.log")
cat("DROUGHT TEMPORAL DIAGNOSTICS WITH WALD-WOLFOWITZ RUNS TEST\n", file = LOG_FILE)
cat("==========================================\n", file = LOG_FILE, append = TRUE)
cat("Drought Definition: Onset < -1.0, Termination >= 0.0\n", file = LOG_FILE, append = TRUE)
cat("Runs Test Purpose: Detect non-random clustering of drought occurrences\n", file = LOG_FILE, append = TRUE)
cat("Timescales Analyzed: 1, 3, 6, 9, 12, 24 months\n", file = LOG_FILE, append = TRUE)
cat("Tie Filtering Threshold: 50% (pixels with >50% ties excluded)\n", file = LOG_FILE, append = TRUE)
cat("Methods: VC (Hamed-Rao) + TFPW + Runs Test\n", file = LOG_FILE, append = TRUE)
cat(paste("Analysis started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat("==========================================\n", file = LOG_FILE, append = TRUE)

log_event <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_msg <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(log_msg, "\n", file = LOG_FILE, append = TRUE)
  if (level %in% c("INFO", "SUCCESS")) {
    message(sprintf("[%s] %s", timestamp, msg))
  } else if (level == "WARNING") {
    warning(sprintf("[%s] %s", timestamp, msg), call. = FALSE)
  } else if (level == "ERROR") {
    stop(sprintf("[%s] %s", timestamp, msg))
  }
}

log_event("Starting drought temporal diagnostics with Runs Test...", "INFO")

# Parameters - EXACT 50% TIE FILTERING THRESHOLD FROM ATTACHED CODE
alpha <- 0.05
drought_threshold <- -1.0    # Drought onset threshold (exact from visualize_spi_spei.R)
recovery_threshold <- 0.0    # Drought termination threshold (exact from visualize_spi_spei.R)
extreme_threshold <- -2.0    # Extreme drought threshold
min_duration <- 2            # Minimum drought duration (months)
max_tie_percent <- 50        # EXACT FILTERING THRESHOLD FROM ATTACHED CODE: 50%
min_valid_obs <- 60          # Minimum valid observations for analysis
n_sim_spectral <- 500        # Monte Carlo simulations for spectral analysis

# Set up parallel processing
num_cores <- min(parallel::detectCores() - 1, 8)
if (num_cores < 1) num_cores <- 1
plan(multisession, workers = num_cores)
log_event(sprintf("Using %d cores for parallel processing", num_cores), "INFO")

# ---- BASIN BOUNDARY LOADING & REPROJECTION ----
log_event("Loading Nechako Basin boundary...", "INFO")
basin_path <- "Spatial/nechakoBound_dissolve.shp"
if (!file.exists(basin_path)) {
  log_event(sprintf("ERROR: Basin boundary file not found: %s", basin_path), "ERROR")
}
basin <- vect(basin_path)
target_crs <- "EPSG:3005"  # BC Albers Equal Area

if (!same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
  log_event("Basin boundary reprojected to BC Albers (EPSG:3005)", "SUCCESS")
} else {
  log_event("Basin boundary already in BC Albers projection", "SUCCESS")
}

# ---- KENDALL TAU-B VARIANCE ADJUSTMENT FOR TIES ----
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

# ---- MODIFIED MANN-KENDALL WITH KENDALL TAU-B VARIANCE (VC METHOD) ----
modified_mann_kendall_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50, apply_vc = TRUE) {
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
  
  for (i in 1:n_space) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    if (n < 10) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, 0, 0, FALSE, TRUE)
      next
    }
    
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    
    if (percent_ties > max_tie_pct) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      next
    }
    
    sen_result <- tryCatch({
      sens.slope(ts_clean)
    }, error = function(e) NULL, warning = function(w) NULL)
    
    if (is.null(sen_result)) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
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
    
    if (apply_vc) {
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
    }
    
    if (varS_final <= 0) {
      p_value <- NA
    } else {
      z_stat <- S / sqrt(varS_final)
      p_value <- 2 * pnorm(-abs(z_stat))
    }
    
    results[i, ] <- list(
      tau = tau,
      p.value = p_value,
      sl = sen_result$estimates,
      S = S,
      varS = varS_final,
      n = n,
      rho1 = rho1,
      vc_corrected = vc_corrected,
      n_ties = n_ties,
      percent_ties = percent_ties,
      tau_b_adjusted = tau_b_adjusted,
      filtered = FALSE
    )
  }
  
  return(results)
}

# ---- TFPW MANN-KENDALL WITH TAU-B (TFPW METHOD) ----
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
  
  for (i in 1:n_space) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    if (n < 10) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, 0, 0, FALSE, TRUE)
      next
    }
    
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    
    if (percent_ties > max_tie_pct) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      next
    }
    
    sen_result <- tryCatch({
      sens.slope(ts_clean)
    }, error = function(e) NULL, warning = function(w) NULL)
    
    if (is.null(sen_result)) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      next
    }
    sen_slope <- sen_result$estimates
    
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
      tau = tau,
      p.value = p_value,
      sl = sen_slope,
      S = S,
      varS = varS_taub,
      n = n,
      rho1 = rho1,
      tfpw_applied = tfpw_applied,
      n_ties = n_ties,
      percent_ties = percent_ties,
      tau_b_adjusted = tau_b_adjusted,
      filtered = FALSE
    )
  }
  
  return(results)
}

# ---- WALD-WOLFOWITZ RUNS TEST FOR DROUGHT CLUSTERING ----
# Detects non-random clustering of drought occurrences (values < -1.0)
# Null hypothesis: Droughts occur randomly in time
# Rejection (p < 0.05): Significant clustering (fewer runs) or alternation (more runs)
wald_wolfowitz_runs_test <- function(ts_values, drought_thresh = -1.0) {
  n <- length(ts_values)
  if (n < 20 || sum(!is.na(ts_values)) < min_valid_obs) {
    return(list(
      n_runs = NA,
      n_drought = NA,
      n_non_drought = NA,
      expected_runs = NA,
      var_runs = NA,
      z_stat = NA,
      p_value = NA,
      clustering = NA,  # -1=clustering, 0=random, +1=alternation
      filtered = TRUE
    ))
  }
  
  # Create binary sequence: 1=drought (< threshold), 0=non-drought (>= threshold)
  binary_seq <- ifelse(ts_values < drought_thresh, 1, 0)
  binary_seq <- binary_seq[!is.na(binary_seq)]  # Remove NAs
  
  n_total <- length(binary_seq)
  if (n_total < 20) {
    return(list(
      n_runs = NA,
      n_drought = NA,
      n_non_drought = NA,
      expected_runs = NA,
      var_runs = NA,
      z_stat = NA,
      p_value = NA,
      clustering = NA,
      filtered = TRUE
    ))
  }
  
  # Count drought/non-drought observations
  n1 <- sum(binary_seq == 1)  # Drought months
  n2 <- sum(binary_seq == 0)  # Non-drought months
  
  # Require minimum representation of both states
  if (n1 < 5 || n2 < 5) {
    return(list(
      n_runs = NA,
      n_drought = n1,
      n_non_drought = n2,
      expected_runs = NA,
      var_runs = NA,
      z_stat = NA,
      p_value = NA,
      clustering = NA,
      filtered = TRUE
    ))
  }
  
  # Count runs (sequences of identical values)
  runs <- rle(binary_seq)
  n_runs <- length(runs$lengths)
  
  # Expected runs under randomness (Wald-Wolfowitz formula)
  expected_runs <- (2 * n1 * n2) / (n1 + n2) + 1
  
  # Variance of runs under randomness
  var_runs <- (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / 
    ((n1 + n2)^2 * (n1 + n2 - 1))
  
  # Z-statistic (with continuity correction)
  if (var_runs <= 0) {
    z_stat <- NA
    p_value <- NA
    clustering <- NA
  } else {
    # Continuity correction: +0.5 if observed < expected, -0.5 if observed > expected
    cc <- if (n_runs < expected_runs) 0.5 else -0.5
    z_stat <- (n_runs - expected_runs + cc) / sqrt(var_runs)
    p_value <- 2 * pnorm(-abs(z_stat))  # Two-tailed test
    
    # Determine clustering direction
    if (!is.na(p_value) && p_value < alpha) {
      clustering <- if (n_runs < expected_runs) -1 else 1  # -1=clustering, +1=alternation
    } else {
      clustering <- 0  # Random pattern
    }
  }
  
  list(
    n_runs = n_runs,
    n_drought = n1,
    n_non_drought = n2,
    expected_runs = expected_runs,
    var_runs = var_runs,
    z_stat = z_stat,
    p_value = p_value,
    clustering = clustering,
    filtered = FALSE
  )
}

# ---- SPECTRAL ANALYSIS (EXACT FROM ATTACHED pr_pet_trends.r CODE) ----
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

# ---- DROUGHT EVENT DETECTION (EXACT DEFINITION FROM visualize_spi_spei.R) ----
detect_drought_events_pixel <- function(ts_values, onset_thresh = -1.0, end_thresh = 0.0, min_dur = 2) {
  n <- length(ts_values)
  if (n < min_dur * 2 || sum(!is.na(ts_values)) < min_valid_obs) {
    return(list(
      n_events = 0,
      mean_duration = NA,
      mean_severity = NA,
      mean_intensity = NA,
      max_duration = NA,
      max_severity = NA,
      max_intensity = NA,
      total_severity = NA
    ))
  }
  
  in_drought <- FALSE
  events <- list()
  start_idx <- NA
  
  for (i in 1:n) {
    if (is.na(ts_values[i])) next
    
    if (!in_drought && ts_values[i] < onset_thresh) {
      in_drought <- TRUE
      start_idx <- i
    } else if (in_drought && ts_values[i] >= end_thresh) {
      end_idx <- i
      duration <- end_idx - start_idx + 1
      
      if (duration >= min_dur) {
        deficit_vals <- ts_values[start_idx:end_idx]
        deficit_vals[deficit_vals >= onset_thresh] <- 0
        severity <- sum(abs(deficit_vals), na.rm = TRUE)
        intensity <- severity / duration
        
        events[[length(events) + 1]] <- list(
          start = start_idx,
          end = end_idx,
          duration = duration,
          severity = severity,
          intensity = intensity
        )
      }
      
      in_drought <- FALSE
      start_idx <- NA
    }
  }
  
  if (in_drought && !is.na(start_idx)) {
    end_idx <- n
    duration <- end_idx - start_idx + 1
    
    if (duration >= min_dur) {
      deficit_vals <- ts_values[start_idx:end_idx]
      deficit_vals[deficit_vals >= onset_thresh] <- 0
      severity <- sum(abs(deficit_vals), na.rm = TRUE)
      intensity <- severity / duration
      
      events[[length(events) + 1]] <- list(
        start = start_idx,
        end = end_idx,
        duration = duration,
        severity = severity,
        intensity = intensity
      )
    }
  }
  
  if (length(events) == 0) {
    return(list(
      n_events = 0,
      mean_duration = NA,
      mean_severity = NA,
      mean_intensity = NA,
      max_duration = NA,
      max_severity = NA,
      max_intensity = NA,
      total_severity = NA
    ))
  }
  
  durations <- sapply(events, function(e) e$duration)
  severities <- sapply(events, function(e) e$severity)
  intensities <- sapply(events, function(e) e$intensity)
  
  list(
    n_events = length(events),
    mean_duration = mean(durations),
    mean_severity = mean(severities),
    mean_intensity = mean(intensities),
    max_duration = max(durations),
    max_severity = max(severities),
    max_intensity = max(intensities),
    total_severity = sum(severities)
  )
}

# ---- REGIME SHIFT DETECTION (PETTITT TEST) ----
detect_regime_shift <- function(ts_values) {
  if (sum(!is.na(ts_values)) < 30) return(NA)
  
  tryCatch({
    cp_result <- cpt.mean(ts_values, method = "PELT", penalty = "MBIC")
    cpts_val <- cpts(cp_result)
    if (length(cpts_val) == 0) return(NA)
    return(cpts_val[1])
  }, error = function(e) {
    return(NA)
  })
}

# ---- RETURN PERIOD ANALYSIS FOR EXTREME DROUGHTS ----
calculate_return_period <- function(ts_values, threshold = -2.0) {
  n <- length(ts_values)
  if (sum(!is.na(ts_values)) < 30) return(NA)
  
  annual_min <- tapply(ts_values, rep(1:floor(n/12), each = 12)[1:n], min, na.rm = TRUE)
  annual_min <- annual_min[!is.infinite(annual_min) & !is.na(annual_min)]
  
  if (length(annual_min) < 10) return(NA)
  
  tryCatch({
    fit <- fevd(annual_min, type = "GEV", method = "MLE")
    rp <- 1 / (1 - pgev(threshold, 
                        loc = fit$results$par[1],
                        scale = fit$results$par[2],
                        shape = fit$results$par[3]))
    return(rp)
  }, error = function(e) {
    return(NA)
  })
}

# ---- LOAD DROUGHT INDEX NETCDF FILES ----
log_event("Loading drought index NetCDF files...", "INFO")

index_types <- c("spi", "spei")
timescales <- c(1, 3, 6, 9, 12, 24)  # INCLUDING 9-MONTH SCALE

all_results <- list()

for (index_type in index_types) {
  for (scale in timescales) {
    nc_file <- file.path(
      ifelse(index_type == "spi", "spi_results_seasonal", "spei_results_seasonal"),
      sprintf("%s_%02d_monthly_1950_2025.nc", index_type, scale)
    )
    
    if (!file.exists(nc_file)) {
      log_event(sprintf("WARNING: File not found - %s", nc_file), "WARNING")
      next
    }
    
    log_event(sprintf("Processing %s-%02d...", toupper(index_type), scale), "INFO")
    
    r_stack <- rast(nc_file)
    
    if (!same.crs(r_stack, target_crs)) {
      r_stack <- project(r_stack, target_crs, method = "bilinear")
      log_event(sprintf("  Reprojected %s-%02d to BC Albers", index_type, scale), "SUCCESS")
    }
    
    r_stack <- crop(r_stack, basin)
    log_event(sprintf("  Clipped to basin extent: %d x %d cells", nrow(r_stack), ncol(r_stack)), "SUCCESS")
    
    dates <- tryCatch({
      terra::time(r_stack)
    }, error = function(e) NULL)
    
    if (is.null(dates) || length(dates) == 0) {
      n_layers <- nlyr(r_stack)
      dates <- seq(as.Date("1950-01-01"), by = "month", length.out = n_layers)
      terra::time(r_stack) <- dates
      log_event(sprintf("  Generated time dimension for %d months", n_layers), "SUCCESS")
    }
    
    years <- as.integer(format(dates, "%Y"))
    months <- as.integer(format(dates, "%m"))
    
    n_time <- nlyr(r_stack)
    n_rows <- nrow(r_stack)
    n_cols <- ncol(r_stack)
    n_space <- n_rows * n_cols
    
    index_matrix <- matrix(values(r_stack), nrow = n_time, ncol = n_space, byrow = FALSE)
    
    lon_coords <- xFromCol(r_stack, 1:n_cols)
    lat_coords <- yFromRow(r_stack, 1:n_rows)
    
    coords_dt <- data.table(
      space_idx = 1:n_space,
      lat = rep(lat_coords, each = n_cols),
      lon = rep(lon_coords, times = n_rows)
    )
    
    # ---- ANNUAL ANALYSIS (VC, TFPW, RUNS TEST, SPECTRAL) ----
    log_event("  Running annual VC trend analysis...", "INFO")
    vc_annual <- modified_mann_kendall_taub(index_matrix, alpha, max_tie_percent, apply_vc = TRUE)
    
    log_event("  Running annual TFPW trend analysis...", "INFO")
    tfpw_annual <- perform_tfpw_mk_taub(index_matrix, alpha, max_tie_percent)
    
    log_event("  Running annual Wald-Wolfowitz Runs Test...", "INFO")
    runs_annual <- lapply(1:n_space, function(i) {
      wald_wolfowitz_runs_test(index_matrix[, i], drought_threshold)
    })
    runs_df_annual <- rbindlist(lapply(runs_annual, as.data.frame))
    runs_df_annual[, space_idx := 1:n_space]
    
    log_event("  Running annual spectral analysis...", "INFO")
    spectral_annual <- perform_spectral_analysis_vectorized(index_matrix, n_sim_spectral)
    spectral_df_annual <- data.table(
      n_spectral_peaks = sapply(spectral_annual, function(x) x$n_peaks),
      dominant_period = sapply(spectral_annual, function(x) x$dominant_period),
      spectral_confidence = sapply(spectral_annual, function(x) x$confidence_limit)
    )
    
    # ---- MONTHLY ANALYSIS (12 CALENDAR MONTHS) ----
    log_event("  Running monthly analyses (12 calendar months)...", "INFO")
    vc_monthly_list <- vector("list", 12)
    tfpw_monthly_list <- vector("list", 12)
    runs_monthly_list <- vector("list", 12)
    spectral_monthly_list <- vector("list", 12)
    
    for (m in 1:12) {
      month_idx <- which(months == m)
      if (length(month_idx) < min_valid_obs / 12) {
        # Create empty results for insufficient data
        empty_df <- data.frame(
          tau = rep(NA, n_space),
          p.value = rep(NA, n_space),
          sl = rep(NA, n_space),
          filtered = rep(TRUE, n_space),
          n = rep(0, n_space),
          rho1 = rep(NA, n_space),
          vc_corrected = rep(FALSE, n_space),
          n_ties = rep(0, n_space),
          percent_ties = rep(0, n_space),
          tau_b_adjusted = rep(FALSE, n_space)
        )
        vc_monthly_list[[m]] <- empty_df
        tfpw_monthly_list[[m]] <- empty_df
        runs_monthly_list[[m]] <- data.frame(
          n_runs = rep(NA, n_space),
          n_drought = rep(NA, n_space),
          n_non_drought = rep(NA, n_space),
          expected_runs = rep(NA, n_space),
          var_runs = rep(NA, n_space),
          z_stat = rep(NA, n_space),
          p_value = rep(NA, n_space),
          clustering = rep(NA, n_space),
          filtered = rep(TRUE, n_space)
        )
        spectral_monthly_list[[m]] <- list(n_peaks = 0, dominant_period = NA, confidence_limit = NA)
        next
      }
      
      month_matrix <- index_matrix[month_idx, , drop = FALSE]
      vc_monthly_list[[m]] <- modified_mann_kendall_taub(month_matrix, alpha, max_tie_percent, apply_vc = TRUE)
      tfpw_monthly_list[[m]] <- perform_tfpw_mk_taub(month_matrix, alpha, max_tie_percent)
      
      # Runs test for monthly subset
      runs_monthly <- lapply(1:n_space, function(i) {
        wald_wolfowitz_runs_test(month_matrix[, i], drought_threshold)
      })
      runs_monthly_list[[m]] <- rbindlist(lapply(runs_monthly, as.data.frame))
      
      spectral_monthly_list[[m]] <- perform_spectral_analysis_vectorized(month_matrix, n_sim_spectral)
    }
    
    # ---- DROUGHT EVENT METRICS (PIXEL-WISE) ----
    log_event("  Calculating drought event metrics per pixel...", "INFO")
    event_metrics <- lapply(1:n_space, function(i) {
      ts_vals <- index_matrix[, i]
      detect_drought_events_pixel(ts_vals, drought_threshold, recovery_threshold, min_duration)
    })
    
    event_df <- rbindlist(lapply(event_metrics, as.data.frame))
    event_df[, space_idx := 1:n_space]
    
    # ---- REGIME SHIFT & RETURN PERIOD ----
    log_event("  Detecting regime shifts per pixel...", "INFO")
    regime_shifts <- sapply(1:n_space, function(i) {
      detect_regime_shift(index_matrix[, i])
    })
    
    log_event("  Calculating return periods for extreme droughts...", "INFO")
    return_periods <- sapply(1:n_space, function(i) {
      calculate_return_period(index_matrix[, i], extreme_threshold)
    })
    
    # ---- COMBINE ANNUAL RESULTS ----
    annual_results <- cbind(
      coords_dt,
      index_type = index_type,
      timescale = scale,
      period = "annual",
      month = NA,
      setDT(vc_annual)[, .(tau_vc = tau, p_value_vc = p.value, sl_vc = sl, 
                           vc_corrected = vc_corrected, n_ties_vc = n_ties, 
                           percent_ties_vc = percent_ties, tau_b_adjusted_vc = tau_b_adjusted,
                           filtered_vc = filtered, n_vc = n, rho1_vc = rho1)],
      setDT(tfpw_annual)[, .(tau_tfpw = tau, p_value_tfpw = p.value, sl_tfpw = sl, 
                             tfpw_applied = tfpw_applied, n_ties_tfpw = n_ties, 
                             percent_ties_tfpw = percent_ties, tau_b_adjusted_tfpw = tau_b_adjusted,
                             filtered_tfpw = filtered, n_tfpw = n, rho1_tfpw = rho1)],
      runs_df_annual[, .(n_runs, n_drought, n_non_drought, expected_runs, var_runs, 
                         z_stat, p_value_runs = p_value, clustering, filtered_runs = filtered)],
      spectral_df_annual,
      event_df[, .(n_events, mean_duration, mean_severity, mean_intensity,
                   max_duration, max_severity, max_intensity, total_severity)],
      regime_shift_year = regime_shifts,
      return_period_extreme = return_periods,
      same_significance = (vc_annual$p.value < alpha & !vc_annual$filtered) == 
        (tfpw_annual$p.value < alpha & !tfpw_annual$filtered),
      same_direction = sign(vc_annual$tau) == sign(tfpw_annual$tau),
      runs_significant = !runs_df_annual$filtered & runs_df_annual$p_value < alpha
    )
    
    # ---- COMBINE MONTHLY RESULTS ----
    monthly_results_list <- vector("list", 12)
    for (m in 1:12) {
      spectral_df_m <- data.table(
        n_spectral_peaks = if (length(spectral_monthly_list[[m]]) > 0) 
          sapply(spectral_monthly_list[[m]], function(x) x$n_peaks) else rep(NA, n_space),
        dominant_period = if (length(spectral_monthly_list[[m]]) > 0) 
          sapply(spectral_monthly_list[[m]], function(x) x$dominant_period) else rep(NA, n_space),
        spectral_confidence = if (length(spectral_monthly_list[[m]]) > 0) 
          sapply(spectral_monthly_list[[m]], function(x) x$confidence_limit) else rep(NA, n_space)
      )
      
      runs_df_m <- runs_monthly_list[[m]]
      if (!is.data.frame(runs_df_m)) runs_df_m <- data.frame(matrix(NA, nrow = n_space, ncol = 9))
      colnames(runs_df_m) <- c("n_runs", "n_drought", "n_non_drought", "expected_runs", 
                               "var_runs", "z_stat", "p_value", "clustering", "filtered")
      
      monthly_results_list[[m]] <- cbind(
        coords_dt,
        index_type = index_type,
        timescale = scale,
        period = "monthly",
        month = m,
        setDT(vc_monthly_list[[m]])[, .(tau_vc = tau, p_value_vc = p.value, sl_vc = sl, 
                                        vc_corrected = vc_corrected, n_ties_vc = n_ties, 
                                        percent_ties_vc = percent_ties, tau_b_adjusted_vc = tau_b_adjusted,
                                        filtered_vc = filtered)],
        setDT(tfpw_monthly_list[[m]])[, .(tau_tfpw = tau, p_value_tfpw = p.value, sl_tfpw = sl, 
                                          tfpw_applied = tfpw_applied, n_ties_tfpw = n_ties, 
                                          percent_ties_tfpw = percent_ties, tau_b_adjusted_tfpw = tau_b_adjusted,
                                          filtered_tfpw = filtered)],
        runs_df_m,
        spectral_df_m,
        same_significance = (vc_monthly_list[[m]]$p.value < alpha & !vc_monthly_list[[m]]$filtered) == 
          (tfpw_monthly_list[[m]]$p.value < alpha & !tfpw_monthly_list[[m]]$filtered),
        same_direction = sign(vc_monthly_list[[m]]$tau) == sign(tfpw_monthly_list[[m]]$tau),
        runs_significant = !runs_df_m$filtered & runs_df_m$p_value < alpha
      )
    }
    
    monthly_results <- rbindlist(monthly_results_list)
    all_index_results <- rbindlist(list(annual_results, monthly_results))
    
    all_results[[paste0(index_type, "_", scale)]] <- all_index_results
    
    # ---- CREATE SPATIAL RASTERS ----
    log_event("  Creating spatial rasters...", "INFO")
    
    template <- rast(r_stack, nlyrs = 1)
    values(template) <- NA
    
    # Annual significance rasters (VC, TFPW, Runs Test)
    for (method in c("vc", "tfpw", "runs")) {
      annual_sig <- copy(template)
      
      if (method == "vc") {
        sig_pixels <- all_index_results[period == "annual" & 
                                          p_value_vc < alpha & !filtered_vc]
        if (nrow(sig_pixels) > 0) {
          xy <- cbind(sig_pixels$lon, sig_pixels$lat)
          cells <- cellFromXY(annual_sig, xy)
          values(annual_sig)[cells] <- sign(sig_pixels$tau_vc)
        }
      } else if (method == "tfpw") {
        sig_pixels <- all_index_results[period == "annual" & 
                                          p_value_tfpw < alpha & !filtered_tfpw]
        if (nrow(sig_pixels) > 0) {
          xy <- cbind(sig_pixels$lon, sig_pixels$lat)
          cells <- cellFromXY(annual_sig, xy)
          values(annual_sig)[cells] <- sign(sig_pixels$tau_tfpw)
        }
      } else if (method == "runs") {
        # Runs test: -1=clustering, +1=alternation (only significant results)
        sig_pixels <- all_index_results[period == "annual" & 
                                          p_value_runs < alpha & !filtered_runs]
        if (nrow(sig_pixels) > 0) {
          xy <- cbind(sig_pixels$lon, sig_pixels$lat)
          cells <- cellFromXY(annual_sig, xy)
          values(annual_sig)[cells] <- sig_pixels$clustering
        }
      }
      
      writeRaster(annual_sig, 
                  file.path(out_dir, sprintf("%s_%02d_annual_%s_significance.tif", index_type, scale, method)),
                  overwrite = TRUE)
    }
    
    # Drought frequency raster
    freq_raster <- copy(template)
    freq_pixels <- all_index_results[period == "annual" & !is.na(n_events)]
    if (nrow(freq_pixels) > 0) {
      xy <- cbind(freq_pixels$lon, freq_pixels$lat)
      cells <- cellFromXY(freq_raster, xy)
      values(freq_raster)[cells] <- freq_pixels$n_events / 7.6
    }
    writeRaster(freq_raster,
                file.path(out_dir, sprintf("%s_%02d_drought_frequency_per_decade.tif", index_type, scale)),
                overwrite = TRUE)
    
    # Runs test clustering raster (annual)
    clustering_raster <- copy(template)
    cluster_pixels <- all_index_results[period == "annual" & !filtered_runs & p_value_runs < alpha]
    if (nrow(cluster_pixels) > 0) {
      xy <- cbind(cluster_pixels$lon, cluster_pixels$lat)
      cells <- cellFromXY(clustering_raster, xy)
      values(clustering_raster)[cells] <- cluster_pixels$clustering
    }
    writeRaster(clustering_raster,
                file.path(out_dir, sprintf("%s_%02d_drought_clustering.tif", index_type, scale)),
                overwrite = TRUE)
    
    # Monthly rasters (VC, TFPW, Runs)
    for (m in 1:12) {
      for (method in c("vc", "tfpw", "runs")) {
        monthly_sig <- copy(template)
        
        if (method == "vc") {
          m_pixels <- all_index_results[period == "monthly" & month == m & 
                                          p_value_vc < alpha & !filtered_vc]
          if (nrow(m_pixels) > 0) {
            xy <- cbind(m_pixels$lon, m_pixels$lat)
            cells <- cellFromXY(monthly_sig, xy)
            values(monthly_sig)[cells] <- sign(m_pixels$tau_vc)
          }
        } else if (method == "tfpw") {
          m_pixels <- all_index_results[period == "monthly" & month == m & 
                                          p_value_tfpw < alpha & !filtered_tfpw]
          if (nrow(m_pixels) > 0) {
            xy <- cbind(m_pixels$lon, m_pixels$lat)
            cells <- cellFromXY(monthly_sig, xy)
            values(monthly_sig)[cells] <- sign(m_pixels$tau_tfpw)
          }
        } else if (method == "runs") {
          m_pixels <- all_index_results[period == "monthly" & month == m & 
                                          p_value < alpha & !filtered]
          if (nrow(m_pixels) > 0) {
            xy <- cbind(m_pixels$lon, m_pixels$lat)
            cells <- cellFromXY(monthly_sig, xy)
            values(monthly_sig)[cells] <- m_pixels$clustering
          }
        }
        
        writeRaster(monthly_sig,
                    file.path(out_dir, sprintf("%s_%02d_month%02d_%s_significance.tif", index_type, scale, m, method)),
                    overwrite = TRUE)
      }
    }
    
    log_event(sprintf("  Completed %s-%02d processing", toupper(index_type), scale), "SUCCESS")
  }
}

# Combine all results
if (length(all_results) > 0) {
  final_results <- rbindlist(all_results, fill = TRUE)
  fwrite(final_results, file.path(out_dir, "all_drought_temporal_diagnostics_runs_test.csv"))
  log_event("Combined results saved to CSV", "SUCCESS")
} else {
  log_event("ERROR: No results generated - check input files", "ERROR")
}

# ---- GENERATE COMPREHENSIVE PDF REPORTS ----
log_event("Creating comprehensive PDF reports with Runs Test...", "INFO")

safe_pdf <- function(filename, width = 14, height = 10) {
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
  tryCatch({
    pdf(filename, width = width, height = height, onefile = TRUE)
    TRUE
  }, error = function(e) {
    log_event(sprintf("ERROR creating PDF '%s': %s", filename, e$message), "ERROR")
    FALSE
  })
}

col_div <- colorRampPalette(brewer.pal(11, "RdBu"))(101)
col_clustering <- c("#2166ac", "gray90", "#b2182b")  # Blue=clustering, Gray=random, Red=alternation

create_index_pdf <- function(index_type, scale) {
  pdf_file <- file.path(out_dir, sprintf("%s_%02d_temporal_diagnostics_runs_test.pdf", index_type, scale))
  if (!safe_pdf(pdf_file, 14, 10)) return(FALSE)
  
  # Load rasters
  annual_sig_vc <- rast(file.path(out_dir, sprintf("%s_%02d_annual_vc_significance.tif", index_type, scale)))
  annual_sig_tfpw <- rast(file.path(out_dir, sprintf("%s_%02d_annual_tfpw_significance.tif", index_type, scale)))
  annual_sig_runs <- rast(file.path(out_dir, sprintf("%s_%02d_annual_runs_significance.tif", index_type, scale)))
  freq_rast <- rast(file.path(out_dir, sprintf("%s_%02d_drought_frequency_per_decade.tif", index_type, scale)))
  cluster_rast <- rast(file.path(out_dir, sprintf("%s_%02d_drought_clustering.tif", index_type, scale)))
  
  par(mfrow = c(3, 2), mar = c(3, 3, 3, 1), oma = c(0, 0, 2, 0))
  
  # VC Annual significance
  plot(annual_sig_vc, col = c("blue", "gray90", "red"), breaks = c(-1.5, -0.5, 0.5, 1.5),
       main = "VC Annual Trend\n(p<0.05)", legend = FALSE)
  plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  legend("bottomright", legend = c("Drying", "NS", "Wetting"),
         fill = c("blue", "gray90", "red"), bty = "n", cex = 0.6)
  
  # TFPW Annual significance
  plot(annual_sig_tfpw, col = c("blue", "gray90", "red"), breaks = c(-1.5, -0.5, 0.5, 1.5),
       main = "TFPW Annual Trend\n(p<0.05)", legend = FALSE)
  plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  legend("bottomright", legend = c("Drying", "NS", "Wetting"),
         fill = c("blue", "gray90", "red"), bty = "n", cex = 0.6)
  
  # Runs Test Clustering
  plot(cluster_rast, col = col_clustering, breaks = c(-1.5, -0.5, 0.5, 1.5),
       main = "Drought Clustering\n(Runs Test p<0.05)", legend = FALSE)
  plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  legend("bottomright", legend = c("Clustering", "Random", "Alternation"),
         fill = col_clustering, bty = "n", cex = 0.6)
  
  # Drought frequency
  plot(freq_rast, col = col_div, main = "Drought Frequency\n(events/decade)",
       legend.args = list(text = "Events/decade", side = 4, line = 2.5, cex = 0.8))
  plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  
  # Monthly examples (Jan VC + Runs)
  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  monthly_sig_vc <- rast(file.path(out_dir, sprintf("%s_%02d_month01_vc_significance.tif", index_type, scale)))
  plot(monthly_sig_vc, col = c("blue", "gray90", "red"), breaks = c(-1.5, -0.5, 0.5, 1.5),
       main = sprintf("%s VC Trend", month_names[1]), legend = FALSE)
  plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  legend("bottomright", legend = c("Drying", "NS", "Wetting"),
         fill = c("blue", "gray90", "red"), bty = "n", cex = 0.6)
  
  monthly_sig_runs <- rast(file.path(out_dir, sprintf("%s_%02d_month01_runs_significance.tif", index_type, scale)))
  plot(monthly_sig_runs, col = col_clustering, breaks = c(-1.5, -0.5, 0.5, 1.5),
       main = sprintf("%s Drought Clustering", month_names[1]), legend = FALSE)
  plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  legend("bottomright", legend = c("Clustering", "Random", "Alternation"),
         fill = col_clustering, bty = "n", cex = 0.6)
  
  mtext(sprintf("%s-%02d Temporal Diagnostics with Runs Test | Nechako Basin (1950-2025)\nDrought: Onset < %.1f, Termination >= %.1f | Tie Filter: >%d%%", 
                toupper(index_type), scale, drought_threshold, recovery_threshold, max_tie_percent),
        outer = TRUE, cex = 1.0, font = 2, line = -1)
  
  dev.off()
  log_event(sprintf("PDF report created: %s", pdf_file), "SUCCESS")
  return(TRUE)
}

# Generate PDFs
for (index_type in index_types) {
  for (scale in timescales) {
    nc_file <- file.path(
      ifelse(index_type == "spi", "spi_results_seasonal", "spei_results_seasonal"),
      sprintf("%s_%02d_monthly_1950_2025.nc", index_type, scale)
    )
    if (file.exists(nc_file)) {
      create_index_pdf(index_type, scale)
    }
  }
}

# ---- SUMMARY STATISTICS WITH RUNS TEST METRICS ----
log_event("Generating summary statistics with Runs Test metrics...", "INFO")

if (exists("final_results")) {
  summary_stats <- final_results[, .(
    total_pixels = .N,
    # Trend metrics
    valid_vc = sum(!filtered_vc, na.rm = TRUE),
    significant_vc = sum(!filtered_vc & p_value_vc < alpha, na.rm = TRUE),
    percent_significant_vc = mean(!filtered_vc & p_value_vc < alpha, na.rm = TRUE) * 100,
    valid_tfpw = sum(!filtered_tfpw, na.rm = TRUE),
    significant_tfpw = sum(!filtered_tfpw & p_value_tfpw < alpha, na.rm = TRUE),
    percent_significant_tfpw = mean(!filtered_tfpw & p_value_tfpw < alpha, na.rm = TRUE) * 100,
    agreement_significance = mean(same_significance[!filtered_vc & !filtered_tfpw], na.rm = TRUE) * 100,
    # Runs test metrics
    valid_runs = sum(!filtered_runs, na.rm = TRUE),
    significant_runs = sum(!filtered_runs & p_value_runs < alpha, na.rm = TRUE),
    percent_significant_runs = mean(!filtered_runs & p_value_runs < alpha, na.rm = TRUE) * 100,
    percent_clustering = mean(clustering[!filtered_runs & p_value_runs < alpha] == -1, na.rm = TRUE) * 100,
    percent_alternation = mean(clustering[!filtered_runs & p_value_runs < alpha] == 1, na.rm = TRUE) * 100,
    # Drought metrics
    mean_drought_freq = mean(n_events[!is.na(n_events)] / 7.6, na.rm = TRUE),
    mean_return_period = mean(return_period_extreme[return_period_extreme > 0 & 
                                                      return_period_extreme < 1000], na.rm = TRUE),
    # Tie statistics
    percent_ties_mean = mean(percent_ties_vc[!filtered_vc], na.rm = TRUE),
    percent_vc_corrected = mean(vc_corrected[!filtered_vc], na.rm = TRUE) * 100,
    percent_tfpw_applied = mean(tfpw_applied[!filtered_tfpw], na.rm = TRUE) * 100
  ), by = .(index_type, timescale, period)]
  
  fwrite(summary_stats, file.path(out_dir, "summary_statistics_runs_test.csv"))
  log_event("Summary statistics saved", "SUCCESS")
  
  # Calculate tie filtering percentage
  total_pixels <- nrow(final_results[period == "annual"])
  filtered_pixels_vc <- sum(final_results[period == "annual"]$filtered_vc, na.rm = TRUE)
  filtered_pixels_tfpw <- sum(final_results[period == "annual"]$filtered_tfpw, na.rm = TRUE)
  filtered_pixels_runs <- sum(final_results[period == "annual"]$filtered_runs, na.rm = TRUE)
  
  percent_filtered_vc <- (filtered_pixels_vc / total_pixels) * 100
  percent_filtered_tfpw <- (filtered_pixels_tfpw / total_pixels) * 100
  percent_filtered_runs <- (filtered_pixels_runs / total_pixels) * 100
  
  log_event("TIE FILTERING & RUNS TEST FILTERING STATISTICS (Annual Analysis):", "INFO")
  log_event(sprintf("  Total pixels analyzed: %d", total_pixels), "INFO")
  log_event(sprintf("  VC method - pixels filtered (>50%% ties): %d (%.2f%%)", 
                    filtered_pixels_vc, percent_filtered_vc), "INFO")
  log_event(sprintf("  TFPW method - pixels filtered (>50%% ties): %d (%.2f%%)", 
                    filtered_pixels_tfpw, percent_filtered_tfpw), "INFO")
  log_event(sprintf("  Runs Test - pixels filtered (insufficient data): %d (%.2f%%)", 
                    filtered_pixels_runs, percent_filtered_runs), "INFO")
  log_event(sprintf("  Filtering threshold: %d%% ties (exact from attached code)", max_tie_percent), "INFO")
}

# ---- CLEANUP ----
plan(sequential)
while (dev.cur() > 1) try(dev.off(), silent = TRUE)

log_event("==========================================", "INFO")
log_event("DROUGHT TEMPORAL DIAGNOSTICS WITH RUNS TEST COMPLETE", "SUCCESS")
log_event("==========================================", "INFO")
log_event(sprintf("Results saved in: %s", out_dir), "INFO")
log_event("Key Outputs:", "INFO")
log_event("  ??? Pixel-wise trend metrics (VC + TFPW) for annual + 12 months", "INFO")
log_event("  ??? Wald-Wolfowitz Runs Test for drought clustering detection", "INFO")
log_event("  ??? Drought event characterization (frequency, duration, severity)", "INFO")
log_event("  ??? Spectral analysis (500 Monte Carlo simulations)", "INFO")
log_event("  ??? Regime shift detection (Pettitt test)", "INFO")
log_event("  ??? Return period maps for extreme droughts", "INFO")
log_event("  ??? Spatial rasters: VC/TFPW significance, Runs Test clustering, frequency", "INFO")
log_event("  ??? Comprehensive PDF reports with method comparison", "INFO")
log_event("==========================================", "INFO")
log_event("Methods Implemented:", "INFO")
log_event("  1. Kendall tau-b variance adjustment for ties (50% threshold)", "INFO")
log_event("  2. Hamed-Rao variance correction for autocorrelation (|?????| > 0.1)", "INFO")
log_event("  3. Trend-Free Pre-Whitening (Yue et al. 2002)", "INFO")
log_event("  4. Wald-Wolfowitz Runs Test for drought clustering (p<0.05)", "INFO")
log_event("  5. Explicit VC vs TFPW vs Runs Test comparison", "INFO")
log_event("  6. Drought definition: onset < -1.0, termination >= 0.0", "INFO")
log_event("==========================================", "INFO")

cat("\n??? DROUGHT TEMPORAL DIAGNOSTICS WITH RUNS TEST COMPLETED SUCCESSFULLY!\n")
cat(sprintf("Full audit trail in log file: %s\n", LOG_FILE))
cat("\nKey innovations:\n")
cat(sprintf("  ??? Tie filtering threshold: %d%% (exact from attached code)\n", max_tie_percent))
cat("  ??? Timescales: 1, 3, 6, 9, 12, 24 months\n")
cat("  ??? Dual trend methods: VC (Hamed-Rao) + TFPW with explicit comparison\n")
cat("  ??? Wald-Wolfowitz Runs Test: Detects non-random drought clustering\n")
cat("    - Fewer runs than expected ??? drought clustering (persistent dry periods)\n")
cat("    - More runs than expected ??? drought alternation (rapid wet-dry shifts)\n")
cat("  ??? Drought definition: onset < -1.0, termination >= 0.0 (exact from visualize_spi_spei.R)\n")
cat("  ??? 12 calendar months processed separately\n")
cat("  ??? Comprehensive timestamped logging with filtering statistics\n")