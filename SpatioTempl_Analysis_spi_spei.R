####################################################################################
# PIXEL-WISE TEMPORAL DIAGNOSTICS FOR DROUGHT INDICES (SPI/SPEI) WITH RUNS TEST
# Nechako Basin Analysis - 76 Years (1950-2025)
# Drought Definition: Onset < -1.0, Termination >= 0.0
# SPI/SPEI Timescales: 1, 3, 6, 9, 12, 24 months | Tie Filtering: 50% threshold
# DATA STRUCTURE: Monthly composites in CSV format (spi_01_month01_Jan.csv = all Jan values 1950-2025)
# RECONSTRUCTION: Interleaves monthly composites into chronological series
# Trend Detection: It uses two variations of the Mann-Kendall test—Variance Correction (VC) and 
# Trend-Free Pre-Whitening (TFPW)—to identify significant long-term increases or decreases in
# drought intensity while accounting for data "ties" and serial correlation.
####################################################################################
####################################################################################

library(terra)
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
library(tseries)


# Output directory
setwd("D:/Nechako_Drought/")
out_dir <- "temporal_spi_spei"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Log file setup
LOG_FILE <- file.path(out_dir, "drought_temporal_diagnostics_runs_test.log")
cat("DROUGHT TEMPORAL DIAGNOSTICS WITH WALD-WOLFOWITZ RUNS TEST\n", file = LOG_FILE)
cat("==========================================\n", file = LOG_FILE, append = TRUE)
cat("FIXED VERSION - Manual Sen's slope calculation\n", file = LOG_FILE, append = TRUE)
cat("Drought Definition: Onset < -0.52, Termination >= -0\n", file = LOG_FILE, append = TRUE)
cat("Runs Test Purpose: Detect non-random clustering of drought occurrences\n", file = LOG_FILE, append = TRUE)
cat("Timescales Analyzed: 1, 3, 6, 9, 12, 24 months\n", file = LOG_FILE, append = TRUE)
cat("Tie Filtering Threshold: 50% (pixels with >50% ties excluded)\n", file = LOG_FILE, append = TRUE)
cat("Methods: VC (Hamed-Rao) + TFPW + Runs Test\n", file = LOG_FILE, append = TRUE)
cat("Data Structure: Monthly composites (CSV) reconstructed into chronological series\n", file = LOG_FILE, append = TRUE)
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

log_event("Starting drought temporal diagnostics with Runs Test (CSV input)...", "INFO")

# Parameters
alpha <- 0.05
drought_threshold <- -0.52
recovery_threshold <- 0
extreme_threshold <- -2.0
min_duration <- 2
max_tie_percent <- 50
min_valid_obs <- 60
n_sim_spectral <- 500
n_years_total <- 76

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
target_crs <- "EPSG:3005"

if (!same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
  log_event("Basin boundary reprojected to BC Albers (EPSG:3005)", "SUCCESS")
} else {
  log_event("Basin boundary already in BC Albers projection", "SUCCESS")
}

# ===============================================================================
# FIX 1: MANUAL SEN'S SLOPE CALCULATION
# ===============================================================================
calculate_sens_slope_manual <- function(x) {
  # Calculate Sen's slope manually: median of all pairwise slopes
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

# ===============================================================================
# FIX 2: MODIFIED MANN-KENDALL WITH MANUAL SEN'S SLOPE (VC METHOD)
# ===============================================================================
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
    
    # ===== FIXED: Use manual Sen's slope calculation =====
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      next
    }
    # ===== END FIX =====
    
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
      sl = sen_slope,  # FIXED: using manual calculation
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

# ===============================================================================
# FIX 3: TFPW MANN-KENDALL WITH MANUAL SEN'S SLOPE (TFPW METHOD)
# ===============================================================================
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
    
    # ===== FIXED: Use manual Sen's slope calculation =====
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      results[i, ] <- list(NA, NA, NA, NA, NA, n, NA, FALSE, n_ties, percent_ties, FALSE, TRUE)
      next
    }
    # ===== END FIX =====
    
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
      sl = sen_slope,  # FIXED: using manual calculation
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
      clustering = NA,
      filtered = TRUE
    ))
  }
  
  binary_seq <- ifelse(ts_values < drought_thresh, 1, 0)
  binary_seq <- binary_seq[!is.na(binary_seq)]
  
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
  
  n1 <- sum(binary_seq == 1)
  n2 <- sum(binary_seq == 0)
  
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
  
  runs <- rle(binary_seq)
  n_runs <- length(runs$lengths)
  
  expected_runs <- (2 * n1 * n2) / (n1 + n2) + 1
  var_runs <- (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / 
    ((n1 + n2)^2 * (n1 + n2 - 1))
  
  if (var_runs <= 0) {
    z_stat <- NA
    p_value <- NA
    clustering <- NA
  } else {
    cc <- if (n_runs < expected_runs) 0.5 else -0.5
    z_stat <- (n_runs - expected_runs + cc) / sqrt(var_runs)
    p_value <- 2 * pnorm(-abs(z_stat))
    
    if (!is.na(p_value) && p_value < alpha) {
      clustering <- if (n_runs < expected_runs) -1 else 1
    } else {
      clustering <- 0
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

# ---- SPECTRAL ANALYSIS ----
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
    
    # FIXED: Use manual Sen's slope
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    if (is.na(sen_slope)) sen_slope <- 0
    
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

# ---- DROUGHT EVENT DETECTION ----
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

# ================================================================================
# CSV FILE READER & RECONSTRUCTION
# ================================================================================
read_monthly_csv_and_reconstruct <- function(index_type, scale, subfolder) {
  search_path <- file.path(getwd(), subfolder)
  
  if (!dir.exists(search_path)) {
    log_event(sprintf("ERROR: Subfolder not found: %s", search_path), "ERROR")
    return(NULL)
  }
  
  pattern <- sprintf("%s_0?%d_month\\d+.*\\.csv$", index_type, scale)
  
  monthly_files <- list.files(
    path = search_path,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(monthly_files) < 12) {
    log_event(sprintf("DEBUG: Only found %d files in %s matching pattern %s", 
                      length(monthly_files), subfolder, pattern), "WARNING")
    return(NULL)
  }
  
  month_nums <- as.integer(gsub(".*month(\\d+).*", "\\1", basename(monthly_files), ignore.case = TRUE))
  
  ord <- order(month_nums)
  monthly_files <- monthly_files[ord]
  
  log_event(sprintf("Reconstructing %s-%02d from %s", toupper(index_type), scale, subfolder), "INFO")
  
  df_template <- fread(monthly_files[1], header = FALSE, na.strings = c("", "NA", ","), fill = TRUE)
  n_years <- ncol(df_template) - 2
  n_pixels <- nrow(df_template)
  
  coords <- df_template[, 1:2, with = FALSE]
  setnames(coords, c("x", "y"))
  
  full_matrix <- matrix(NA, nrow = n_years * 12, ncol = n_pixels)
  
  for (m in 1:12) {
    df_month <- fread(monthly_files[m], header = FALSE, na.strings = c("", "NA", ","), fill = TRUE)
    year_values <- as.matrix(df_month[, 3:ncol(df_month), with = FALSE])
    for (yr in 1:n_years) {
      full_matrix[(yr - 1) * 12 + m, ] <- year_values[, yr]
    }
  }
  
  list(index_matrix = full_matrix, coords = coords, n_years = n_years, n_pixels = n_pixels)
}

# ================================================================================
# MAIN PROCESSING LOOP
# ================================================================================
log_event("Searching for monthly composite CSV files...", "INFO")

index_types <- c("spi", "spei")
timescales <- c(1, 3, 6, 9, 12, 24)

all_results <- list()

for (index_type in index_types) {
  subfolder_name <- paste0(index_type, "_results_seasonal")
  for (scale in timescales) {
    reconstruction <- read_monthly_csv_and_reconstruct(index_type, scale, subfolder_name)
    
    if (is.null(reconstruction)) {
      log_event(sprintf("Skipping %s-%02d due to insufficient files", index_type, scale), "WARNING")
      next
    }
    
    index_matrix <- reconstruction$index_matrix
    coords <- reconstruction$coords
    n_years <- reconstruction$n_years
    n_pixels <- reconstruction$n_pixels
    
    # Generate date sequences
    n_months <- nrow(index_matrix)
    years_full <- rep(1:n_years, each = 12)[1:n_months]
    months_full <- rep(1:12, times = n_years)[1:n_months]
    
    log_event(sprintf("Processing %s-%02d with reconstructed time series...", toupper(index_type), scale), "INFO")
    
    coords_dt <- data.table(
      space_idx = 1:n_pixels,
      lon = coords$x,
      lat = coords$y
    )
    
    # ---- ANNUAL ANALYSIS ----
    log_event("  Running annual VC trend analysis...", "INFO")
    vc_annual <- modified_mann_kendall_taub(index_matrix, alpha, max_tie_percent, apply_vc = TRUE)
    
    log_event("  Running annual TFPW trend analysis...", "INFO")
    tfpw_annual <- perform_tfpw_mk_taub(index_matrix, alpha, max_tie_percent)
    
    log_event("  Running annual Wald-Wolfowitz Runs Test...", "INFO")
    runs_annual <- lapply(1:n_pixels, function(i) {
      wald_wolfowitz_runs_test(index_matrix[, i], drought_threshold)
    })
    runs_df_annual <- rbindlist(lapply(runs_annual, as.data.frame))
    runs_df_annual[, space_idx := 1:n_pixels]
    
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
      month_idx <- which(months_full == m)
      if (length(month_idx) < min_valid_obs / 12) {
        # ===== FIX 4: COMPLETE EMPTY DATAFRAME STRUCTURE =====
        empty_df <- data.frame(
          tau = rep(NA, n_pixels),
          p.value = rep(NA, n_pixels),
          sl = rep(NA, n_pixels),
          S = rep(NA, n_pixels),
          varS = rep(NA, n_pixels),
          n = rep(0, n_pixels),
          rho1 = rep(NA, n_pixels),
          vc_corrected = rep(FALSE, n_pixels),
          tfpw_applied = rep(FALSE, n_pixels),  # ADDED THIS LINE
          n_ties = rep(0, n_pixels),
          percent_ties = rep(0, n_pixels),
          tau_b_adjusted = rep(FALSE, n_pixels),
          filtered = rep(TRUE, n_pixels)
        )
        # ===== END FIX =====
        
        vc_monthly_list[[m]] <- empty_df
        tfpw_monthly_list[[m]] <- empty_df
        runs_monthly_list[[m]] <- data.frame(
          n_runs = rep(NA, n_pixels),
          n_drought = rep(NA, n_pixels),
          n_non_drought = rep(NA, n_pixels),
          expected_runs = rep(NA, n_pixels),
          var_runs = rep(NA, n_pixels),
          z_stat = rep(NA, n_pixels),
          p_value = rep(NA, n_pixels),
          clustering = rep(NA, n_pixels),
          filtered = rep(TRUE, n_pixels)
        )
        spectral_monthly_list[[m]] <- list(n_peaks = 0, dominant_period = NA, confidence_limit = NA)
        next
      }
      
      month_matrix <- index_matrix[month_idx, , drop = FALSE]
      vc_monthly_list[[m]] <- modified_mann_kendall_taub(month_matrix, alpha, max_tie_percent, apply_vc = TRUE)
      tfpw_monthly_list[[m]] <- perform_tfpw_mk_taub(month_matrix, alpha, max_tie_percent)
      
      runs_monthly <- lapply(1:n_pixels, function(i) {
        wald_wolfowitz_runs_test(month_matrix[, i], drought_threshold)
      })
      runs_monthly_list[[m]] <- rbindlist(lapply(runs_monthly, as.data.frame))
      
      spectral_monthly_list[[m]] <- perform_spectral_analysis_vectorized(month_matrix, n_sim_spectral)
    }
    
    # ---- DROUGHT EVENT METRICS ----
    log_event("  Calculating drought event metrics per pixel...", "INFO")
    event_metrics <- lapply(1:n_pixels, function(i) {
      ts_vals <- index_matrix[, i]
      detect_drought_events_pixel(ts_vals, drought_threshold, recovery_threshold, min_duration)
    })
    
    event_df <- rbindlist(lapply(event_metrics, as.data.frame))
    event_df[, space_idx := 1:n_pixels]
    
    # ---- REGIME SHIFT & RETURN PERIOD ----
    log_event("  Detecting regime shifts per pixel...", "INFO")
    regime_shifts <- sapply(1:n_pixels, function(i) {
      detect_regime_shift(index_matrix[, i])
    })
    
    log_event("  Calculating return periods for extreme droughts...", "INFO")
    return_periods <- sapply(1:n_pixels, function(i) {
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
      
      safe_extract <- function(list_obj, field) {
        sapply(list_obj, function(x) {
          if (is.list(x) && field %in% names(x)) {
            return(x[[field]])
          } else {
            return(NA)
          }
        })
      }
      
      spectral_df_m <- data.table(
        n_spectral_peaks = if (!is.null(spectral_monthly_list[[m]])) 
          safe_extract(spectral_monthly_list[[m]], "n_peaks") 
        else rep(NA, n_pixels),
        dominant_period = if (!is.null(spectral_monthly_list[[m]])) 
          safe_extract(spectral_monthly_list[[m]], "dominant_period") 
        else rep(NA, n_pixels),
        spectral_confidence = if (!is.null(spectral_monthly_list[[m]])) 
          safe_extract(spectral_monthly_list[[m]], "confidence_limit") 
        else rep(NA, n_pixels)
      )
      
      runs_df_m <- runs_monthly_list[[m]]
      if (!is.data.frame(runs_df_m)) {
        runs_df_m <- data.table(matrix(NA, nrow = n_pixels, ncol = 9))
        colnames(runs_df_m) <- c("n_runs", "n_drought", "n_non_drought", "expected_runs", 
                                 "var_runs", "z_stat", "p_value", "clustering", "filtered")
      } else {
        setDT(runs_df_m)
      }
      
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
        runs_df_m[, .(n_runs, n_drought, n_non_drought, expected_runs, var_runs, 
                      z_stat, p_value_runs = p_value, clustering, filtered_runs = filtered)],
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
    
    template <- rast(ext(basin), nrows = 200, ncols = 200, crs = target_crs)
    template <- rasterize(basin, template, values = NA)
    
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
    
    log_event(sprintf("  Completed %s-%02d processing", toupper(index_type), scale), "SUCCESS")
  }
}

# Combine all results
if (length(all_results) > 0) {
  final_results <- rbindlist(all_results, fill = TRUE)
  fwrite(final_results, file.path(out_dir, "all_drought_temporal_diagnostics_runs_test.csv"))
  log_event("Combined results saved to CSV", "SUCCESS")
  
  # Summary statistics
  n_pixels_total <- nrow(final_results[period == "annual"])
  n_sig_vc <- sum(final_results[period == "annual" & !filtered_vc & p_value_vc < alpha], na.rm = TRUE)
  n_sig_tfpw <- sum(final_results[period == "annual" & !filtered_tfpw & p_value_tfpw < alpha], na.rm = TRUE)
  n_sig_runs <- sum(final_results[period == "annual" & !filtered_runs & p_value_runs < alpha], na.rm = TRUE)
  
  cat("\n=== ANALYSIS SUMMARY ===\n", file = LOG_FILE, append = TRUE)
  cat(sprintf("Total pixels analyzed: %d\n", n_pixels_total), file = LOG_FILE, append = TRUE)
  cat(sprintf("Significant VC trends (p<0.05): %d (%.1f%%)\n", n_sig_vc, n_sig_vc/n_pixels_total*100), file = LOG_FILE, append = TRUE)
  cat(sprintf("Significant TFPW trends (p<0.05): %d (%.1f%%)\n", n_sig_tfpw, n_sig_tfpw/n_pixels_total*100), file = LOG_FILE, append = TRUE)
  cat(sprintf("Significant Runs Test (p<0.05): %d (%.1f%%)\n", n_sig_runs, n_sig_runs/n_pixels_total*100), file = LOG_FILE, append = TRUE)
  cat(sprintf("Results saved to: %s\n", file.path(out_dir, "all_drought_temporal_diagnostics_runs_test.csv")), file = LOG_FILE, append = TRUE)
  cat("==========================================\n", file = LOG_FILE, append = TRUE)
  
  log_event(sprintf("Analysis complete! Results saved to %s", out_dir), "SUCCESS")
} else {
  log_event("ERROR: No results generated - check input files and data structure", "ERROR")
}

plan(sequential)
cat("\n✓ ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("All fixes applied:\n")
cat("  • Manual Sen's slope calculation (robust replacement for sens.slope())\n")
cat("  • Fixed empty_df structure to include all required columns\n")
cat("  • Spectral analysis updated with manual Sen's slope\n")
cat(sprintf("Check %s for full details.\n", LOG_FILE))