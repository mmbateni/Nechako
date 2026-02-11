####################################################################################
# PIXEL-WISE TEMPORAL DIAGNOSTICS FOR DROUGHT INDICES (SPI/SPEI) WITH RUNS TEST
# Nechako Basin Analysis - 76 Years (1950-2025)
# 1. Added diagnostic plots and statistics for VC vs TFPW discrepancy
# 2. Added diagnostics for runs test dispersion dominance
# 3. Created publication-quality maps with proper legends
# 4. Added comparative visualizations and summary reports
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
library(viridis)
library(scales)

# Output directory
setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "temporal_spi_spei"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Create diagnostic subdirectories
diag_dir <- file.path(out_dir, "diagnostics")
maps_dir <- file.path(out_dir, "maps_with_legends")
reports_dir <- file.path(out_dir, "reports")
for (d in c(diag_dir, maps_dir, reports_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Log file setup
LOG_FILE <- file.path(out_dir, "drought_temporal_diagnostics.log")
cat("ENHANCED DROUGHT TEMPORAL DIAGNOSTICS WITH COMPREHENSIVE DIAGNOSTICS\n", file = LOG_FILE)
cat("======================================================================\n", file = LOG_FILE, append = TRUE)
cat(paste("Analysis started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat("======================================================================\n", file = LOG_FILE, append = TRUE)

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

log_event("Starting enhanced drought temporal diagnostics...", "INFO")

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
num_cores <- parallel::detectCores(logical = FALSE) - 1
if (is.na(num_cores) || num_cores < 1) num_cores <- 1
options(future.globals.maxSize = 2000 * 1024^2)
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
# CORE FUNCTIONS
# ===============================================================================

calculate_sens_slope_manual <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  d_mat <- outer(x, x, "-")
  t_mat <- outer(1:n, 1:n, "-")
  upper.tri.inds <- which(upper.tri(d_mat))
  slopes <- d_mat[upper.tri.inds] / t_mat[upper.tri.inds]
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

modified_mann_kendall_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50, apply_vc = TRUE) {
  n_space <- ncol(ts_matrix)
  
  results_list <- future_lapply(1:n_space, function(i) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    if (n < 10) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA, 
                  vc_corrected=FALSE, n_ties=0, percent_ties=0, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    
    if (percent_ties > max_tie_pct) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA, 
                  vc_corrected=FALSE, n_ties=n_ties, percent_ties=percent_ties, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA, 
                  vc_corrected=FALSE, n_ties=n_ties, percent_ties=percent_ties, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    diff_mat <- outer(ts_clean, ts_clean, "-")
    S <- sum(sign(diff_mat[upper.tri(diff_mat)]))
    
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
    
    return(list(
      tau = tau, p.value = p_value, sl = sen_slope, S = S, varS = varS_final,
      n = n, rho1 = rho1, vc_corrected = vc_corrected, n_ties = n_ties,
      percent_ties = percent_ties, tau_b_adjusted = tau_b_adjusted, filtered = FALSE
    ))
  }, future.seed = TRUE)
  
  return(rbindlist(results_list))
}

perform_tfpw_mk_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50) {
  n_space <- ncol(ts_matrix)
  
  results_list <- future_lapply(1:n_space, function(i) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    if (n < 10) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA, 
                  tfpw_applied=FALSE, n_ties=0, percent_ties=0, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    
    if (percent_ties > max_tie_pct) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA, 
                  tfpw_applied=FALSE, n_ties=n_ties, percent_ties=percent_ties, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA, 
                  tfpw_applied=FALSE, n_ties=n_ties, percent_ties=percent_ties, tau_b_adjusted=FALSE, filtered=TRUE))
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
      prewhitened_detrended[2:n] <- detrended_series[2:n] - rho1 * detrended_series[1:(n-1)]
      corrected_series <- prewhitened_detrended + trend_line
      tfpw_applied <- TRUE
    } else {
      corrected_series <- ts_clean
    }
    
    diff_mat <- outer(corrected_series, corrected_series, "-")
    S <- sum(sign(diff_mat[upper.tri(diff_mat)]))
    
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
    
    return(list(
      tau = tau, p.value = p_value, sl = sen_slope, S = S, varS = varS_taub,
      n = n, rho1 = rho1, tfpw_applied = tfpw_applied, n_ties = n_ties,
      percent_ties = percent_ties, tau_b_adjusted = tau_b_adjusted, filtered = FALSE
    ))
  }, future.seed = TRUE)
  
  return(rbindlist(results_list))
}

wald_wolfowitz_runs_test <- function(ts_values, drought_thresh = -1.0) {
  n <- length(ts_values)
  if (n < 20 || sum(!is.na(ts_values)) < min_valid_obs) {
    return(list(n_runs=NA, n_drought=NA, n_non_drought=NA, expected_runs=NA, 
                var_runs=NA, z_stat=NA, p_value=NA, clustering=NA, filtered=TRUE))
  }
  
  binary_seq <- ifelse(ts_values < drought_thresh, 1, 0)
  binary_seq <- binary_seq[!is.na(binary_seq)]
  
  n_total <- length(binary_seq)
  if (n_total < 20) return(list(n_runs=NA, n_drought=NA, n_non_drought=NA, expected_runs=NA, 
                                var_runs=NA, z_stat=NA, p_value=NA, clustering=NA, filtered=TRUE))
  
  n1 <- sum(binary_seq == 1)
  n2 <- sum(binary_seq == 0)
  
  if (n1 < 5 || n2 < 5) {
    return(list(n_runs=NA, n_drought=n1, n_non_drought=n2, expected_runs=NA, 
                var_runs=NA, z_stat=NA, p_value=NA, clustering=NA, filtered=TRUE))
  }
  
  runs <- rle(binary_seq)
  n_runs <- length(runs$lengths)
  
  expected_runs <- (2 * n1 * n2) / (n1 + n2) + 1
  var_runs <- (2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)) / ((n1 + n2)^2 * (n1 + n2 - 1))
  
  if (var_runs <= 0) {
    z_stat <- NA; p_value <- NA; clustering <- NA
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
  
  list(n_runs=n_runs, n_drought=n1, n_non_drought=n2, expected_runs=expected_runs, 
       var_runs=var_runs, z_stat=z_stat, p_value=p_value, clustering=clustering, filtered=FALSE)
}

perform_spectral_analysis_vectorized <- function(ts_matrix, n_sim = 500, alpha = 0.05) {
  n_space <- ncol(ts_matrix)
  
  results_list <- future_lapply(1:n_space, function(i) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    if (n < 20) {
      return(list(n_peaks = 0, dominant_period = NA, confidence_limit = NA))
    }
    
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    if (is.na(sen_slope)) sen_slope <- 0
    
    detrended_series <- ts_clean - sen_slope * (1:n)
    detrended_series <- detrended_series - mean(detrended_series)
    
    fft_result <- fft(detrended_series)
    spectral_density <- Mod(fft_result[1:(n/2)])^2 / n
    
    frequencies <- seq(0, 0.5, length.out = n/2)
    
    random_matrix <- matrix(rnorm(n * n_sim, mean = 0, sd = sd(detrended_series)), nrow = n)
    fft_rand <- mvfft(random_matrix)
    spectral_rand <- Mod(fft_rand[1:(n/2), ])^2 / n
    max_spectra <- apply(spectral_rand, 2, max, na.rm = TRUE)
    conf_limit <- quantile(max_spectra, 1 - alpha, na.rm = TRUE)
    
    significant_peaks <- spectral_density > conf_limit
    peak_indices <- which(significant_peaks & !is.na(significant_peaks))
    
    peak_frequencies <- frequencies[peak_indices]
    peak_periods <- ifelse(peak_frequencies > 0, 1/peak_frequencies, NA)
    
    if (length(peak_periods) > 0) {
      sorted_indices <- order(spectral_density[peak_indices], decreasing = TRUE)
      peak_periods <- peak_periods[sorted_indices]
    }
    
    return(list(
      n_peaks = length(peak_periods),
      dominant_period = if (length(peak_periods) > 0) peak_periods[1] else NA,
      confidence_limit = conf_limit
    ))
  }, future.seed = TRUE)
  
  return(results_list)
}

detect_drought_events_pixel <- function(ts_values, onset_thresh = -1.0, end_thresh = 0.0, min_dur = 2) {
  n <- length(ts_values)
  if (n < min_dur * 2 || sum(!is.na(ts_values)) < min_valid_obs) {
    return(list(n_events=0, mean_duration=NA, mean_severity=NA, mean_intensity=NA, 
                max_duration=NA, max_severity=NA, max_intensity=NA, total_severity=NA))
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
          start = start_idx, end = end_idx, duration = duration,
          severity = severity, intensity = intensity
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
        start = start_idx, end = end_idx, duration = duration,
        severity = severity, intensity = intensity
      )
    }
  }
  
  if (length(events) == 0) {
    return(list(n_events=0, mean_duration=NA, mean_severity=NA, mean_intensity=NA, 
                max_duration=NA, max_severity=NA, max_intensity=NA, total_severity=NA))
  }
  
  durations <- sapply(events, function(e) e$duration)
  severities <- sapply(events, function(e) e$severity)
  intensities <- sapply(events, function(e) e$intensity)
  
  list(
    n_events = length(events),
    mean_duration = mean(durations), mean_severity = mean(severities),
    mean_intensity = mean(intensities), max_duration = max(durations),
    max_severity = max(severities), max_intensity = max(intensities),
    total_severity = sum(severities)
  )
}

detect_regime_shift <- function(ts_values) {
  if (sum(!is.na(ts_values)) < 30) return(NA)
  tryCatch({
    cp_result <- cpt.mean(ts_values, method = "PELT", penalty = "MBIC")
    cpts_val <- cpts(cp_result)
    if (length(cpts_val) == 0) return(NA)
    return(cpts_val[1])
  }, error = function(e) NA)
}

calculate_return_period <- function(ts_values, threshold = -2.0) {
  n <- length(ts_values)
  if (sum(!is.na(ts_values)) < 30) return(NA)
  
  annual_min <- tapply(ts_values, rep(1:floor(n/12), each = 12)[1:n], min, na.rm = TRUE)
  annual_min <- annual_min[!is.infinite(annual_min) & !is.na(annual_min)]
  
  if (length(annual_min) < 10) return(NA)
  
  tryCatch({
    fit <- fevd(annual_min, type = "GEV", method = "MLE")
    rp <- 1 / (1 - pgev(threshold, loc = fit$results$par[1],
                        scale = fit$results$par[2], shape = fit$results$par[3]))
    return(rp)
  }, error = function(e) NA)
}

read_monthly_csv_and_reconstruct <- function(index_type, scale, subfolder) {
  search_path <- file.path(getwd(), subfolder)
  
  if (!dir.exists(search_path)) {
    log_event(sprintf("ERROR: Subfolder not found: %s", search_path), "ERROR")
    return(NULL)
  }
  
  pattern <- sprintf("%s_0?%d_month\\d+.*\\.csv$", index_type, scale)  
  monthly_files <- list.files(
    path = search_path, pattern = pattern, full.names = TRUE, ignore.case = TRUE
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

# ===============================================================================
# NEW: ENHANCED VISUALIZATION FUNCTIONS
# ===============================================================================

create_map_with_legend <- function(raster_data, basin_boundary, output_file, 
                                   title, legend_labels, color_palette,
                                   subtitle = NULL) {
  
  png(output_file, width = 10, height = 8, units = "in", res = 300)
  
  par(mar = c(1, 1, 3, 8), xpd = FALSE)
  
  # Plot raster
  plot(raster_data, col = color_palette, axes = FALSE, box = FALSE,
       legend = FALSE, main = title, cex.main = 1.5, font.main = 2)
  
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = 0.5, cex = 0.9, col = "gray30")
  }
  
  # Add basin boundary
  plot(basin_boundary, add = TRUE, border = "black", lwd = 1.5)
  
  # Add scale bar
  sbar(10000, xy = "bottomleft", type = "bar", below = "meters", divs = 4)
  
  # Add north arrow
  north(xy = "topright", type = 1, cex = 1.2)
  
  # Create custom legend
  par(xpd = TRUE)
  legend_x <- par("usr")[2] + (par("usr")[2] - par("usr")[1]) * 0.02
  legend_y <- mean(par("usr")[3:4])
  
  n_colors <- length(color_palette)
  legend_height <- (par("usr")[4] - par("usr")[3]) * 0.5
  color_height <- legend_height / n_colors
  
  for (i in 1:n_colors) {
    y_bottom <- legend_y - legend_height/2 + (i-1) * color_height
    y_top <- y_bottom + color_height
    
    rect(legend_x, y_bottom, legend_x + (par("usr")[2] - par("usr")[1]) * 0.05, 
         y_top, col = color_palette[i], border = "black", lwd = 0.5)
  }
  
  # Add legend text
  text_x <- legend_x + (par("usr")[2] - par("usr")[1]) * 0.06
  for (i in 1:length(legend_labels)) {
    y_pos <- legend_y - legend_height/2 + (i - 0.5) * color_height
    text(text_x, y_pos, legend_labels[i], pos = 4, cex = 0.9)
  }
  
  # Add legend title
  text(legend_x, legend_y + legend_height/2 + color_height/2, 
       "Classification", pos = 4, cex = 1, font = 2)
  
  # Add metadata box
  metadata_text <- sprintf("Analysis: %s\nPeriod: 1950-2025 (76 years)\nα = %.2f", 
                           title, alpha)
  mtext(metadata_text, side = 1, line = -2, adj = 0.02, cex = 0.7, col = "gray40")
  
  dev.off()
  
  log_event(sprintf("Created map with legend: %s", basename(output_file)), "SUCCESS")
}

# ===============================================================================
# NEW: DIAGNOSTIC FUNCTIONS
# ===============================================================================

# DIAGNOSTIC 1: VC vs TFPW Comparison Analysis
compare_vc_tfpw <- function(vc_results, tfpw_results, index_name, output_dir) {
  
  log_event(sprintf("Running VC vs TFPW diagnostic for %s...", index_name), "INFO")
  
  # Create copies and rename columns to match the function's expected logic
  vc_dt <- copy(vc_results)
  setnames(vc_dt, 
           old = c("filtered", "rho1", "p.value"), 
           new = c("filtered_vc", "rho1_vc", "p_value_vc"), 
           skip_absent = TRUE)
  
  tfpw_dt <- copy(tfpw_results)
  setnames(tfpw_dt, 
           old = c("filtered", "rho1", "p.value"), 
           new = c("filtered_tfpw", "rho1_tfpw", "p_value_tfpw"), 
           skip_absent = TRUE)
  
  # Extract non-filtered results
  vc_valid <- vc_dt[filtered_vc == FALSE] 
  tfpw_valid <- tfpw_dt[filtered_tfpw == FALSE]
  
  # Diagnostic 1: Autocorrelation distribution
  png(file.path(output_dir, sprintf("%s_diagnostic_autocorrelation.png", index_name)),
      width = 12, height = 8, units = "in", res = 300)
  par(mfrow = c(2, 2))
  
  hist(vc_valid$rho1_vc, breaks = 50, main = "VC Method: Autocorrelation Distribution",
       xlab = "Lag-1 Autocorrelation (ρ₁)", col = "skyblue", border = "white")
  abline(v = c(-0.1, 0.1), col = "red", lty = 2, lwd = 2)
  text(0, par("usr")[4] * 0.9, sprintf("Mean ρ₁ = %.3f\nMedian ρ₁ = %.3f", 
                                       mean(vc_valid$rho1_vc, na.rm = TRUE),
                                       median(vc_valid$rho1_vc, na.rm = TRUE)), pos = 4)
  
  hist(tfpw_valid$rho1_tfpw, breaks = 50, main = "TFPW Method: Autocorrelation Distribution",
       xlab = "Lag-1 Autocorrelation (ρ₁) - Detrended", col = "lightcoral", border = "white")
  abline(v = c(-0.1, 0.1), col = "red", lty = 2, lwd = 2)
  text(0, par("usr")[4] * 0.9, sprintf("Mean ρ₁ = %.3f\nMedian ρ₁ = %.3f", 
                                       mean(tfpw_valid$rho1_tfpw, na.rm = TRUE),
                                       median(tfpw_valid$rho1_tfpw, na.rm = TRUE)), pos = 4)
  
  # Diagnostic 2: P-value comparison
  plot(vc_valid$p_value_vc, tfpw_valid$p_value_tfpw, 
       pch = 16, col = scales::alpha("black", 0.3), cex = 0.5,
       xlab = "VC Method p-value", ylab = "TFPW Method p-value",
       main = "P-value Comparison: VC vs TFPW")
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(h = 0.05, v = 0.05, col = "blue", lty = 3)
  
  # Count quadrants
  q1 <- sum(vc_valid$p_value_vc < 0.05 & tfpw_valid$p_value_tfpw < 0.05, na.rm = TRUE)
  q2 <- sum(vc_valid$p_value_vc >= 0.05 & tfpw_valid$p_value_tfpw < 0.05, na.rm = TRUE)
  q3 <- sum(vc_valid$p_value_vc < 0.05 & tfpw_valid$p_value_tfpw >= 0.05, na.rm = TRUE)
  q4 <- sum(vc_valid$p_value_vc >= 0.05 & tfpw_valid$p_value_tfpw >= 0.05, na.rm = TRUE)
  
  legend("topright", 
         legend = c(sprintf("Both significant: %d", q1),
                    sprintf("TFPW only: %d", q2),
                    sprintf("VC only: %d", q3),
                    sprintf("Neither: %d", q4)),
         cex = 0.8, bg = "white")
  
  # Diagnostic 3: Variance inflation
  vc_corrected_idx <- which(vc_valid$vc_corrected)
  if (length(vc_corrected_idx) > 0) {
    # Calculate variance inflation factor
    # This would require storing varS before and after correction
    # For now, plot autocorrelation vs significance
    
    plot(vc_valid$rho1_vc, -log10(vc_valid$p_value_vc), 
         pch = 16, col = alpha("darkblue", 0.3), cex = 0.6,
         xlab = "Autocorrelation (ρ₁)", ylab = "-log10(p-value)",
         main = "VC Method: Autocorrelation Impact on Significance")
    abline(h = -log10(0.05), col = "red", lty = 2, lwd = 2)
    abline(v = c(-0.1, 0.1), col = "orange", lty = 3)
    
    text(0.5, par("usr")[4] * 0.9, 
         sprintf("Pixels with |ρ₁| > 0.1: %d (%.1f%%)\nVC correction applied: %d",
                 sum(abs(vc_valid$rho1_vc) > 0.1, na.rm = TRUE),
                 sum(abs(vc_valid$rho1_vc) > 0.1, na.rm = TRUE) / nrow(vc_valid) * 100,
                 sum(vc_valid$vc_corrected, na.rm = TRUE)),
         pos = 2, cex = 0.8)
  }
  
  dev.off()
  
  # Create diagnostic report
  report_file <- file.path(output_dir, sprintf("%s_VC_vs_TFPW_diagnostic_report.txt", index_name))
  
  cat("=" , rep("=", 78), "=\n", sep = "", file = report_file)
  cat("DIAGNOSTIC REPORT: VC vs TFPW Comparison\n", file = report_file, append = TRUE)
  cat("Index:", index_name, "\n", file = report_file, append = TRUE)
  cat("=" , rep("=", 78), "=\n\n", sep = "", file = report_file, append = TRUE)
  
  cat("1. AUTOCORRELATION STATISTICS\n", file = report_file, append = TRUE)
  cat("   VC Method (original series):\n", file = report_file, append = TRUE)
  cat(sprintf("     Mean ρ₁: %.4f\n", mean(vc_valid$rho1_vc, na.rm = TRUE)), file = report_file, append = TRUE)
  cat(sprintf("     Median ρ₁: %.4f\n", median(vc_valid$rho1_vc, na.rm = TRUE)), file = report_file, append = TRUE)
  cat(sprintf("     SD ρ₁: %.4f\n", sd(vc_valid$rho1_vc, na.rm = TRUE)), file = report_file, append = TRUE)
  cat(sprintf("     |ρ₁| > 0.1: %d pixels (%.1f%%)\n", 
              sum(abs(vc_valid$rho1_vc) > 0.1, na.rm = TRUE),
              sum(abs(vc_valid$rho1_vc) > 0.1, na.rm = TRUE) / nrow(vc_valid) * 100),
      file = report_file, append = TRUE)
  
  cat("\n   TFPW Method (detrended series):\n", file = report_file, append = TRUE)
  cat(sprintf("     Mean ρ₁: %.4f\n", mean(tfpw_valid$rho1_tfpw, na.rm = TRUE)), file = report_file, append = TRUE)
  cat(sprintf("     Median ρ₁: %.4f\n", median(tfpw_valid$rho1_tfpw, na.rm = TRUE)), file = report_file, append = TRUE)
  cat(sprintf("     SD ρ₁: %.4f\n", sd(tfpw_valid$rho1_tfpw, na.rm = TRUE)), file = report_file, append = TRUE)
  cat(sprintf("     |ρ₁| > 0.1: %d pixels (%.1f%%)\n", 
              sum(abs(tfpw_valid$rho1_tfpw) > 0.1, na.rm = TRUE),
              sum(abs(tfpw_valid$rho1_tfpw) > 0.1, na.rm = TRUE) / nrow(tfpw_valid) * 100),
      file = report_file, append = TRUE)
  
  cat("\n2. SIGNIFICANCE DETECTION\n", file = report_file, append = TRUE)
  cat(sprintf("   VC Method: %d pixels significant (%.2f%%)\n",
              sum(vc_valid$p_value_vc < 0.05, na.rm = TRUE),
              sum(vc_valid$p_value_vc < 0.05, na.rm = TRUE) / nrow(vc_valid) * 100),
      file = report_file, append = TRUE)
  cat(sprintf("   TFPW Method: %d pixels significant (%.2f%%)\n",
              sum(tfpw_valid$p_value_tfpw < 0.05, na.rm = TRUE),
              sum(tfpw_valid$p_value_tfpw < 0.05, na.rm = TRUE) / nrow(tfpw_valid) * 100),
      file = report_file, append = TRUE)
  
  cat("\n3. AGREEMENT ANALYSIS\n", file = report_file, append = TRUE)
  cat(sprintf("   Both methods agree (both sig or both not sig): %d pixels (%.1f%%)\n",
              q1 + q4, (q1 + q4) / nrow(vc_valid) * 100),
      file = report_file, append = TRUE)
  cat(sprintf("   Disagreement: %d pixels (%.1f%%)\n",
              q2 + q3, (q2 + q3) / nrow(vc_valid) * 100),
      file = report_file, append = TRUE)
  cat(sprintf("     - TFPW detects but VC doesn't: %d\n", q2), file = report_file, append = TRUE)
  cat(sprintf("     - VC detects but TFPW doesn't: %d\n", q3), file = report_file, append = TRUE)
  
  cat("\n4. POSSIBLE CAUSES OF DISCREPANCY\n", file = report_file, append = TRUE)
  
  if (mean(abs(vc_valid$rho1_vc), na.rm = TRUE) > 0.3) {
    cat("   ⚠ HIGH AUTOCORRELATION DETECTED (mean |ρ₁| > 0.3)\n", file = report_file, append = TRUE)
    cat("     - VC method may be over-inflating variance\n", file = report_file, append = TRUE)
    cat("     - Consider: autocorrelation may be partly due to trends\n", file = report_file, append = TRUE)
  }
  
  if (sum(tfpw_valid$p_value_tfpw < 0.05, na.rm = TRUE) > 
      sum(vc_valid$p_value_vc < 0.05, na.rm = TRUE) * 10) {
    cat("   ⚠ TFPW DETECTS 10× MORE TRENDS THAN VC\n", file = report_file, append = TRUE)
    cat("     - VC variance correction may be too conservative\n", file = report_file, append = TRUE)
    cat("     - TFPW pre-whitening removes trend-induced autocorrelation\n", file = report_file, append = TRUE)
    cat("     - Recommendation: Use TFPW for trend detection in this dataset\n", file = report_file, append = TRUE)
  }
  
  cat("\n5. RECOMMENDATION\n", file = report_file, append = TRUE)
  agreement_pct <- (q1 + q4) / nrow(vc_valid) * 100
  
  if (agreement_pct > 90) {
    cat("   ✓ Methods show high agreement (>90%). Either method is suitable.\n", 
        file = report_file, append = TRUE)
  } else if (q2 > q1 * 2) {
    cat("   → Use TFPW METHOD (preferred)\n", file = report_file, append = TRUE)
    cat("     Reason: VC appears over-conservative due to trend-induced autocorrelation\n", 
        file = report_file, append = TRUE)
  } else if (q3 > q1 * 2) {
    cat("   → Use VC METHOD (preferred)\n", file = report_file, append = TRUE)
    cat("     Reason: TFPW may be removing legitimate autocorrelation structure\n", 
        file = report_file, append = TRUE)
  } else {
    cat("   → Report both methods with caveats\n", file = report_file, append = TRUE)
    cat("     Reason: Substantial disagreement suggests sensitivity to assumptions\n", 
        file = report_file, append = TRUE)
  }
  
  cat("\n" , rep("=", 80), "\n", sep = "", file = report_file, append = TRUE)
  
  log_event(sprintf("VC vs TFPW diagnostic complete: %s", basename(report_file)), "SUCCESS")
}

# DIAGNOSTIC 2: Runs Test Dispersion Analysis
analyze_runs_dispersion <- function(runs_results, index_matrix, index_name, output_dir, 
                                    coords, drought_threshold) {
  
  log_event(sprintf("Running Runs Test diagnostic for %s...", index_name), "INFO")
  
  # Create copy and rename columns to match function logic
  runs_dt <- copy(runs_results)
  setnames(runs_dt, "filtered", "filtered_runs", skip_absent = TRUE)
  
  # Extract non-filtered results
  runs_valid <- runs_dt[filtered_runs == FALSE]
  
  # Diagnostic plots
  png(file.path(output_dir, sprintf("%s_diagnostic_runs_test.png", index_name)),
      width = 14, height = 10, units = "in", res = 300)
  par(mfrow = c(3, 2))
  
  # 1. Distribution of observed vs expected runs
  runs_ratio <- runs_valid$n_runs / runs_valid$expected_runs
  hist(runs_ratio, breaks = 50, main = "Runs Ratio Distribution",
       xlab = "Observed Runs / Expected Runs", col = "lightblue", border = "white")
  abline(v = 1, col = "red", lwd = 2, lty = 2)
  text(par("usr")[2] * 0.7, par("usr")[4] * 0.9,
       sprintf("Mean ratio: %.3f\nMedian ratio: %.3f\n< 1 (clustering): %d\n> 1 (dispersion): %d",
               mean(runs_ratio, na.rm = TRUE),
               median(runs_ratio, na.rm = TRUE),
               sum(runs_ratio < 1, na.rm = TRUE),
               sum(runs_ratio > 1, na.rm = TRUE)),
       pos = 2, cex = 0.8)
  
  # 2. Drought proportion vs runs
  drought_prop <- runs_valid$n_drought / (runs_valid$n_drought + runs_valid$n_non_drought)
  plot(drought_prop, runs_ratio, pch = 16, col = alpha("darkgreen", 0.3),
       xlab = "Proportion of Time in Drought", ylab = "Runs Ratio",
       main = "Drought Frequency vs Temporal Pattern")
  abline(h = 1, col = "red", lty = 2)
  
  # 3. Z-statistic distribution
  hist(runs_valid$z_stat, breaks = 50, main = "Runs Test Z-Statistic Distribution",
       xlab = "Z-statistic", col = "salmon", border = "white")
  abline(v = c(-1.96, 1.96), col = "blue", lty = 2, lwd = 2)
  text(par("usr")[2] * 0.6, par("usr")[4] * 0.9,
       sprintf("Mean Z: %.3f\nZ > 1.96: %d (%.1f%%)\nZ < -1.96: %d (%.1f%%)",
               mean(runs_valid$z_stat, na.rm = TRUE),
               sum(runs_valid$z_stat > 1.96, na.rm = TRUE),
               sum(runs_valid$z_stat > 1.96, na.rm = TRUE) / nrow(runs_valid) * 100,
               sum(runs_valid$z_stat < -1.96, na.rm = TRUE),
               sum(runs_valid$z_stat < -1.96, na.rm = TRUE) / nrow(runs_valid) * 100),
       pos = 2, cex = 0.8)
  
  # 4. P-value distribution
  hist(runs_valid$p_value, breaks = 50, main = "Runs Test P-value Distribution",
       xlab = "P-value", col = "lightgreen", border = "white", xlim = c(0, 1))
  abline(v = 0.05, col = "red", lty = 2, lwd = 2)
  text(0.5, par("usr")[4] * 0.9,
       sprintf("p < 0.05: %d (%.1f%%)",
               sum(runs_valid$p_value < 0.05, na.rm = TRUE),
               sum(runs_valid$p_value < 0.05, na.rm = TRUE) / nrow(runs_valid) * 100),
       pos = 4, cex = 0.9)
  
  # 5. Clustering classification
  clust_counts <- table(runs_valid$clustering)
  barplot(clust_counts, main = "Temporal Pattern Classification",
          names.arg = c("Not Significant", "Clustering", "Dispersion"),
          col = c("gray", "red", "blue"), border = "white",
          ylab = "Number of Pixels")
  text(1:3, clust_counts + max(clust_counts) * 0.05,
       sprintf("n=%d\n(%.1f%%)", clust_counts, clust_counts/sum(clust_counts)*100),
       pos = 3)
  
  # 6. Example time series
  # Select representative pixels
  set.seed(123)
  clust_idx <- which(runs_valid$clustering == -1)
  disp_idx <- which(runs_valid$clustering == 1)
  
  if (length(clust_idx) > 0 && length(disp_idx) > 0) {
    example_clust <- sample(clust_idx, min(1, length(clust_idx)))
    example_disp <- sample(disp_idx, min(1, length(disp_idx)))
    
    ts_clust <- index_matrix[, example_clust]
    ts_disp <- index_matrix[, example_disp]
    
    plot(ts_clust, type = "l", col = "red", lwd = 1.5,
         main = "Example Time Series",
         xlab = "Month", ylab = "Drought Index",
         ylim = range(c(ts_clust, ts_disp), na.rm = TRUE))
    lines(ts_disp, col = "blue", lwd = 1.5)
    abline(h = drought_threshold, col = "black", lty = 2)
    legend("topright", 
           legend = c("Clustering", "Dispersion", "Drought threshold"),
           col = c("red", "blue", "black"),
           lty = c(1, 1, 2), lwd = c(1.5, 1.5, 1),
           cex = 0.8, bg = "white")
  }
  
  dev.off()
  
  # Create diagnostic report
  report_file <- file.path(output_dir, sprintf("%s_runs_test_diagnostic_report.txt", index_name))
  
  cat("=", rep("=", 78), "=\n", sep = "", file = report_file)
  cat("DIAGNOSTIC REPORT: Runs Test Dispersion Analysis\n", file = report_file, append = TRUE)
  cat("Index:", index_name, "\n", file = report_file, append = TRUE)
  cat("Drought Threshold:", drought_threshold, "\n", file = report_file, append = TRUE)
  cat("=", rep("=", 78), "=\n\n", sep = "", file = report_file, append = TRUE)
  
  cat("1. OVERALL STATISTICS\n", file = report_file, append = TRUE)
  cat(sprintf("   Total pixels analyzed: %d\n", nrow(runs_valid)), file = report_file, append = TRUE)
  cat(sprintf("   Pixels with significant pattern (p<0.05): %d (%.1f%%)\n",
              sum(runs_valid$p_value < 0.05, na.rm = TRUE),
              sum(runs_valid$p_value < 0.05, na.rm = TRUE) / nrow(runs_valid) * 100),
      file = report_file, append = TRUE)
  
  cat("\n2. TEMPORAL PATTERN CLASSIFICATION\n", file = report_file, append = TRUE)
  cat(sprintf("   Clustering (observed < expected runs): %d (%.1f%%)\n",
              sum(runs_valid$clustering == -1, na.rm = TRUE),
              sum(runs_valid$clustering == -1, na.rm = TRUE) / nrow(runs_valid) * 100),
      file = report_file, append = TRUE)
  cat(sprintf("   Dispersion (observed > expected runs): %d (%.1f%%)\n",
              sum(runs_valid$clustering == 1, na.rm = TRUE),
              sum(runs_valid$clustering == 1, na.rm = TRUE) / nrow(runs_valid) * 100),
      file = report_file, append = TRUE)
  cat(sprintf("   No significant pattern: %d (%.1f%%)\n",
              sum(runs_valid$clustering == 0, na.rm = TRUE),
              sum(runs_valid$clustering == 0, na.rm = TRUE) / nrow(runs_valid) * 100),
      file = report_file, append = TRUE)
  
  cat("\n3. RUNS STATISTICS\n", file = report_file, append = TRUE)
  cat(sprintf("   Mean observed runs: %.2f\n", mean(runs_valid$n_runs, na.rm = TRUE)), 
      file = report_file, append = TRUE)
  cat(sprintf("   Mean expected runs: %.2f\n", mean(runs_valid$expected_runs, na.rm = TRUE)), 
      file = report_file, append = TRUE)
  cat(sprintf("   Mean runs ratio (obs/exp): %.3f\n", mean(runs_ratio, na.rm = TRUE)), 
      file = report_file, append = TRUE)
  cat(sprintf("   Median runs ratio: %.3f\n", median(runs_ratio, na.rm = TRUE)), 
      file = report_file, append = TRUE)
  
  cat("\n4. DROUGHT CHARACTERISTICS\n", file = report_file, append = TRUE)
  cat(sprintf("   Mean drought proportion: %.3f\n", mean(drought_prop, na.rm = TRUE)), 
      file = report_file, append = TRUE)
  cat(sprintf("   Mean drought months per pixel: %.1f\n", 
              mean(runs_valid$n_drought, na.rm = TRUE)), 
      file = report_file, append = TRUE)
  cat(sprintf("   Mean non-drought months per pixel: %.1f\n", 
              mean(runs_valid$n_non_drought, na.rm = TRUE)), 
      file = report_file, append = TRUE)
  
  cat("\n5. INTERPRETATION OF HIGH DISPERSION\n", file = report_file, append = TRUE)
  
  disp_pct <- sum(runs_valid$clustering == 1, na.rm = TRUE) / nrow(runs_valid) * 100
  
  if (disp_pct > 90) {
    cat("   ⚠ EXTREMELY HIGH DISPERSION DETECTED (>90%)\n\n", file = report_file, append = TRUE)
    cat("   Possible causes:\n", file = report_file, append = TRUE)
    cat("   a) Drought threshold may be too sensitive:\n", file = report_file, append = TRUE)
    cat(sprintf("      - Current threshold: %.2f\n", drought_threshold), file = report_file, append = TRUE)
    cat("      - This may capture normal monthly variability rather than droughts\n", 
        file = report_file, append = TRUE)
    cat("      - Consider stricter threshold (e.g., -1.0 or -1.5)\n", 
        file = report_file, append = TRUE)
    
    cat("\n   b) Monthly timescale artifacts:\n", file = report_file, append = TRUE)
    cat("      - Monthly data shows natural wet/dry alternation\n", file = report_file, append = TRUE)
    cat("      - Seasonal cycles create apparent 'dispersion'\n", file = report_file, append = TRUE)
    cat("      - Consider: aggregate to seasonal or annual timescales\n", 
        file = report_file, append = TRUE)
    
    cat("\n   c) Climate regime characteristics:\n", file = report_file, append = TRUE)
    cat("      - Region may have highly variable precipitation\n", file = report_file, append = TRUE)
    cat("      - Frequent transitions between wet/dry states\n", file = report_file, append = TRUE)
    cat("      - This could be a real climatic feature\n", file = report_file, append = TRUE)
    
  } else if (disp_pct > 70) {
    cat("   ⚠ HIGH DISPERSION DETECTED (>70%)\n", file = report_file, append = TRUE)
    cat("   - Basin shows frequent wet/dry alternation\n", file = report_file, append = TRUE)
    cat("   - Limited multi-month drought persistence\n", file = report_file, append = TRUE)
  } else if (disp_pct < 30) {
    cat("   ✓ BALANCED PATTERN DISTRIBUTION\n", file = report_file, append = TRUE)
    cat("   - Mix of clustering and dispersion is typical\n", file = report_file, append = TRUE)
  }
  
  cat("\n6. RECOMMENDATIONS\n", file = report_file, append = TRUE)
  
  if (disp_pct > 90) {
    cat("   1. Re-run analysis with drought threshold = -1.0\n", file = report_file, append = TRUE)
    cat("   2. Test multiple thresholds: -0.8, -1.0, -1.5, -2.0\n", file = report_file, append = TRUE)
    cat("   3. Examine seasonal patterns separately\n", file = report_file, append = TRUE)
    cat("   4. Calculate runs test on longer timescales (3-month, 6-month)\n", 
        file = report_file, append = TRUE)
    cat("   5. Compare with historical drought records\n", file = report_file, append = TRUE)
  }
  
  cat("\n", rep("=", 80), "\n", sep = "", file = report_file, append = TRUE)
  
  log_event(sprintf("Runs test diagnostic complete: %s", basename(report_file)), "SUCCESS")
}

# ===============================================================================
# MAIN PROCESSING LOOP WITH ENHANCED OUTPUTS
# ===============================================================================

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
    log_event("  Running annual VC trend analysis (Parallel)...", "INFO")
    vc_annual <- modified_mann_kendall_taub(index_matrix, alpha, max_tie_percent, apply_vc = TRUE)
    
    log_event("  Running annual TFPW trend analysis (Parallel)...", "INFO")
    tfpw_annual <- perform_tfpw_mk_taub(index_matrix, alpha, max_tie_percent)
    
    log_event("  Running annual Wald-Wolfowitz Runs Test (Parallel)...", "INFO")
    runs_annual <- future_lapply(1:n_pixels, function(i) {
      wald_wolfowitz_runs_test(index_matrix[, i], drought_threshold)
    }, future.seed = TRUE)
    runs_df_annual <- rbindlist(lapply(runs_annual, as.data.frame))
    runs_df_annual[, space_idx := 1:n_pixels]
    
    log_event("  Running annual spectral analysis (Parallel)...", "INFO")
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
        empty_df <- data.frame(
          tau = rep(NA, n_pixels), p.value = rep(NA, n_pixels), sl = rep(NA, n_pixels),
          S = rep(NA, n_pixels), varS = rep(NA, n_pixels), n = rep(0, n_pixels),
          rho1 = rep(NA, n_pixels), vc_corrected = rep(FALSE, n_pixels),
          tfpw_applied = rep(FALSE, n_pixels), n_ties = rep(0, n_pixels),
          percent_ties = rep(0, n_pixels), tau_b_adjusted = rep(FALSE, n_pixels),
          filtered = rep(TRUE, n_pixels)
        )
        
        vc_monthly_list[[m]] <- empty_df
        tfpw_monthly_list[[m]] <- empty_df
        runs_monthly_list[[m]] <- data.frame(
          n_runs = rep(NA, n_pixels), n_drought = rep(NA, n_pixels),
          n_non_drought = rep(NA, n_pixels), expected_runs = rep(NA, n_pixels),
          var_runs = rep(NA, n_pixels), z_stat = rep(NA, n_pixels),
          p_value = rep(NA, n_pixels), clustering = rep(NA, n_pixels),
          filtered = rep(TRUE, n_pixels)
        )
        spectral_monthly_list[[m]] <- list(n_peaks = 0, dominant_period = NA, confidence_limit = NA)
        next
      }
      
      month_matrix <- index_matrix[month_idx, , drop = FALSE]
      vc_monthly_list[[m]] <- modified_mann_kendall_taub(month_matrix, alpha, max_tie_percent, apply_vc = TRUE)
      tfpw_monthly_list[[m]] <- perform_tfpw_mk_taub(month_matrix, alpha, max_tie_percent)
      
      runs_monthly <- future_lapply(1:n_pixels, function(i) {
        wald_wolfowitz_runs_test(month_matrix[, i], drought_threshold)
      }, future.seed = TRUE)
      runs_monthly_list[[m]] <- rbindlist(lapply(runs_monthly, as.data.frame))
      
      spectral_monthly_list[[m]] <- perform_spectral_analysis_vectorized(month_matrix, n_sim_spectral)
    }
    
    # ---- DROUGHT EVENT METRICS ----
    log_event("  Calculating drought event metrics per pixel (Parallel)...", "INFO")
    event_metrics <- future_lapply(1:n_pixels, function(i) {
      ts_vals <- index_matrix[, i]
      detect_drought_events_pixel(ts_vals, drought_threshold, recovery_threshold, min_duration)
    }, future.seed = TRUE)
    
    event_df <- rbindlist(lapply(event_metrics, as.data.frame))
    event_df[, space_idx := 1:n_pixels]
    
    # ---- REGIME SHIFT & RETURN PERIOD ----
    log_event("  Detecting regime shifts & return periods (Parallel)...", "INFO")
    regime_shifts <- future_sapply(1:n_pixels, function(i) {
      detect_regime_shift(index_matrix[, i])
    }, future.seed = TRUE)
    
    return_periods <- future_sapply(1:n_pixels, function(i) {
      calculate_return_period(index_matrix[, i], extreme_threshold)
    }, future.seed = TRUE)
    
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
    all_index_results <- rbindlist(list(annual_results, monthly_results), fill = TRUE)
    
    all_results[[paste0(index_type, "_", scale)]] <- all_index_results
    
    # ===============================================================================
    # NEW: RUN DIAGNOSTICS
    # ===============================================================================
    
    index_name <- sprintf("%s_%02d", toupper(index_type), scale)
    
    # Diagnostic 1: VC vs TFPW comparison
    compare_vc_tfpw(vc_annual, tfpw_annual, index_name, diag_dir)
    
    # Diagnostic 2: Runs test dispersion analysis
    analyze_runs_dispersion(runs_df_annual, index_matrix, index_name, diag_dir, 
                            coords, drought_threshold)
    
    # ===============================================================================
    # NEW: CREATE ENHANCED MAPS WITH PROPER LEGENDS
    # ===============================================================================
    
    log_event("  Creating enhanced spatial rasters with legends...", "INFO")
    
    template <- rast(ext(basin), nrows = 200, ncols = 200, crs = target_crs)
    template <- rasterize(basin, template, values = NA)
    
    # Function to populate raster
    populate_raster <- function(raster_template, data_table, value_col, filter_col = NULL) {
      result_rast <- copy(raster_template)
      
      if (!is.null(filter_col)) {
        valid_data <- data_table[!get(filter_col) & !is.na(get(value_col))]
      } else {
        valid_data <- data_table[!is.na(get(value_col))]
      }
      
      if (nrow(valid_data) > 0) {
        xy <- cbind(valid_data$lon, valid_data$lat)
        cells <- cellFromXY(result_rast, xy)
        valid_mask <- !is.na(cells)
        
        if (any(valid_mask)) {
          values(result_rast)[cells[valid_mask]] <- valid_data[[value_col]][valid_mask]
        }
      }
      
      return(result_rast)
    }
    
    # 1. Runs test significance map
    runs_sig_rast <- populate_raster(template, annual_results, "clustering", "filtered_runs")
    runs_sig_rast <- runs_sig_rast * (annual_results$p_value_runs[match(cells(runs_sig_rast), 
                                                                        cellFromXY(template, 
                                                                                   cbind(annual_results$lon, annual_results$lat)))] < alpha)
    runs_sig_rast[runs_sig_rast == 0] <- NA
    
    writeRaster(runs_sig_rast, 
                file.path(out_dir, sprintf("%s_%02d_annual_runs_significance.tif", index_type, scale)),
                overwrite = TRUE)
    
    create_map_with_legend(
      raster_data = runs_sig_rast,
      basin_boundary = basin,
      output_file = file.path(maps_dir, sprintf("%s_%02d_runs_test_map.png", index_type, scale)),
      title = sprintf("%s-%02d: Wald-Wolfowitz Runs Test", toupper(index_type), scale),
      subtitle = sprintf("Temporal Clustering Analysis (p < %.2f) | Threshold: %.2f", alpha, drought_threshold),
      legend_labels = c("Significant Clustering\n(Too few runs)", 
                        "Significant Dispersion\n(Too many runs)"),
      color_palette = c("#d62728", "#1f77b4")  # red, blue
    )
    
    # 2. VC trend significance map
    vc_sig_rast <- populate_raster(template, annual_results, "tau_vc", "filtered_vc")
    vc_sig_rast <- vc_sig_rast * (annual_results$p_value_vc[match(cells(vc_sig_rast), 
                                                                  cellFromXY(template, 
                                                                             cbind(annual_results$lon, annual_results$lat)))] < alpha)
    vc_sig_rast[vc_sig_rast == 0] <- NA
    vc_sig_rast <- sign(vc_sig_rast)
    
    writeRaster(vc_sig_rast,
                file.path(out_dir, sprintf("%s_%02d_annual_vc_significance.tif", index_type, scale)),
                overwrite = TRUE)
    
    if (sum(!is.na(values(vc_sig_rast))) > 0) {
      create_map_with_legend(
        raster_data = vc_sig_rast,
        basin_boundary = basin,
        output_file = file.path(maps_dir, sprintf("%s_%02d_vc_trends_map.png", index_type, scale)),
        title = sprintf("%s-%02d: Mann-Kendall (VC Method)", toupper(index_type), scale),
        subtitle = sprintf("Variance-Corrected Trend Analysis (p < %.2f)", alpha),
        legend_labels = c("Significant Drying Trend", "Significant Wetting Trend"),
        color_palette = c("#8B4513", "#228B22")  # brown, green
      )
    }
    
    # 3. TFPW trend significance map
    tfpw_sig_rast <- populate_raster(template, annual_results, "tau_tfpw", "filtered_tfpw")
    tfpw_sig_rast <- tfpw_sig_rast * (annual_results$p_value_tfpw[match(cells(tfpw_sig_rast), 
                                                                        cellFromXY(template, 
                                                                                   cbind(annual_results$lon, annual_results$lat)))] < alpha)
    tfpw_sig_rast[tfpw_sig_rast == 0] <- NA
    tfpw_sig_rast <- sign(tfpw_sig_rast)
    
    writeRaster(tfpw_sig_rast,
                file.path(out_dir, sprintf("%s_%02d_annual_tfpw_significance.tif", index_type, scale)),
                overwrite = TRUE)
    
    if (sum(!is.na(values(tfpw_sig_rast))) > 0) {
      create_map_with_legend(
        raster_data = tfpw_sig_rast,
        basin_boundary = basin,
        output_file = file.path(maps_dir, sprintf("%s_%02d_tfpw_trends_map.png", index_type, scale)),
        title = sprintf("%s-%02d: Mann-Kendall (TFPW Method)", toupper(index_type), scale),
        subtitle = sprintf("Trend-Free Pre-Whitened Trend Analysis (p < %.2f)", alpha),
        legend_labels = c("Significant Drying Trend", "Significant Wetting Trend"),
        color_palette = c("#8B4513", "#228B22")  # brown, green
      )
    }
    
    # 4. Drought frequency map
    freq_raster <- populate_raster(template, annual_results, "n_events")
    freq_raster <- freq_raster / 7.6  # Convert to events per decade
    
    writeRaster(freq_raster,
                file.path(out_dir, sprintf("%s_%02d_drought_frequency_per_decade.tif", index_type, scale)),
                overwrite = TRUE)
    
    if (sum(!is.na(values(freq_raster))) > 0) {
      # Create categorical version for better visualization
      freq_cat <- freq_raster
      values(freq_cat) <- cut(values(freq_raster), 
                              breaks = c(0, 2, 4, 6, 8, Inf),
                              labels = FALSE)
      
      create_map_with_legend(
        raster_data = freq_cat,
        basin_boundary = basin,
        output_file = file.path(maps_dir, sprintf("%s_%02d_drought_frequency_map.png", index_type, scale)),
        title = sprintf("%s-%02d: Drought Frequency", toupper(index_type), scale),
        subtitle = "Number of drought events per decade",
        legend_labels = c("0-2 events", "2-4 events", "4-6 events", "6-8 events", ">8 events"),
        color_palette = c("#ffffcc", "#fee391", "#fec44f", "#fe9929", "#d95f0e")
      )
    }
    
    log_event(sprintf("  Completed %s-%02d processing with diagnostics", toupper(index_type), scale), "SUCCESS")
  }
}

# Combine all results
if (length(all_results) > 0) {
  final_results <- rbindlist(all_results, fill = TRUE)
  fwrite(final_results, file.path(out_dir, "all_drought_temporal_diagnostics.csv"))
  log_event("Combined results saved to CSV", "SUCCESS")
  
  # Summary statistics
  n_pixels_total <- nrow(final_results[period == "annual"])
  n_sig_vc <- nrow(final_results[period == "annual" & !filtered_vc & p_value_vc < alpha, ])
  n_sig_tfpw <- nrow(final_results[period == "annual" & !filtered_tfpw & p_value_tfpw < alpha, ])
  n_sig_runs <- nrow(final_results[period == "annual" & !filtered_runs & p_value_runs < alpha, ])
  
  # Create comprehensive summary report
  summary_file <- file.path(reports_dir, "ANALYSIS_SUMMARY.txt")
  
  cat("=" , rep("=", 78), "=\n", sep = "", file = summary_file)
  cat("ENHANCED DROUGHT TEMPORAL DIAGNOSTICS - ANALYSIS SUMMARY\n", file = summary_file, append = TRUE)
  cat("=" , rep("=", 78), "=\n", sep = "", file = summary_file, append = TRUE)
  cat(paste("Analysis completed:", Sys.time(), "\n"), file = summary_file, append = TRUE)
  cat(paste("Total processing time:", difftime(Sys.time(), as.POSIXct("2026-02-10 14:42:58"), units = "mins"), "minutes\n"), 
      file = summary_file, append = TRUE)
  cat("=" , rep("=", 78), "=\n\n", sep = "", file = summary_file, append = TRUE)
  
  cat("OVERALL STATISTICS\n", file = summary_file, append = TRUE)
  cat(sprintf("Total pixels analyzed: %d\n", n_pixels_total), file = summary_file, append = TRUE)
  cat(sprintf("Indices processed: %d (SPI + SPEI across 6 timescales)\n", 
              length(unique(final_results$timescale)) * 2), file = summary_file, append = TRUE)
  cat(sprintf("Total time series: %d\n", n_pixels_total), file = summary_file, append = TRUE)
  
  cat("\nTREND DETECTION RESULTS (Annual Analysis)\n", file = summary_file, append = TRUE)
  cat(sprintf("  Variance Correction (VC) method:\n"), file = summary_file, append = TRUE)
  cat(sprintf("    Significant trends (p<0.05): %d pixels (%.2f%%)\n", 
              n_sig_vc, n_sig_vc/n_pixels_total*100), file = summary_file, append = TRUE)
  
  cat(sprintf("\n  Trend-Free Pre-Whitening (TFPW) method:\n"), file = summary_file, append = TRUE)
  cat(sprintf("    Significant trends (p<0.05): %d pixels (%.2f%%)\n", 
              n_sig_tfpw, n_sig_tfpw/n_pixels_total*100), file = summary_file, append = TRUE)
  
  cat(sprintf("\n  ⚠ DISCREPANCY: TFPW detects %.1fx more trends than VC\n", 
              n_sig_tfpw / max(n_sig_vc, 1)), file = summary_file, append = TRUE)
  cat("    → See diagnostic reports in /diagnostics/ for detailed analysis\n", 
      file = summary_file, append = TRUE)
  
  cat("\nTEMPORAL PATTERN RESULTS (Runs Test)\n", file = summary_file, append = TRUE)
  cat(sprintf("  Significant temporal patterns (p<0.05): %d pixels (%.2f%%)\n", 
              n_sig_runs, n_sig_runs/n_pixels_total*100), file = summary_file, append = TRUE)
  
  n_clustering <- nrow(final_results[period == "annual" & clustering == -1 & p_value_runs < alpha])
  n_dispersion <- nrow(final_results[period == "annual" & clustering == 1 & p_value_runs < alpha])
  
  cat(sprintf("    Clustering: %d pixels (%.1f%% of significant)\n", 
              n_clustering, n_clustering/max(n_sig_runs, 1)*100), file = summary_file, append = TRUE)
  cat(sprintf("    Dispersion: %d pixels (%.1f%% of significant)\n", 
              n_dispersion, n_dispersion/max(n_sig_runs, 1)*100), file = summary_file, append = TRUE)
  
  if (n_dispersion / max(n_sig_runs, 1) > 0.9) {
    cat("\n  ⚠ EXTREMELY HIGH DISPERSION DETECTED (>90%)\n", file = summary_file, append = TRUE)
    cat("    → See diagnostic reports for possible causes and recommendations\n", 
        file = summary_file, append = TRUE)
  }
  
  cat("\nOUTPUT FILES CREATED\n", file = summary_file, append = TRUE)
  cat(sprintf("  Main results CSV: %s\n", 
              file.path(out_dir, "all_drought_temporal_diagnostics.csv")), 
      file = summary_file, append = TRUE)
  cat(sprintf("  Diagnostic reports: %s/*.txt\n", diag_dir), file = summary_file, append = TRUE)
  cat(sprintf("  Diagnostic plots: %s/*.png\n", diag_dir), file = summary_file, append = TRUE)
  cat(sprintf("  Maps with legends: %s/*.png\n", maps_dir), file = summary_file, append = TRUE)
  cat(sprintf("  GeoTIFF rasters: %s/*.tif\n", out_dir), file = summary_file, append = TRUE)
  
  cat("\n" , rep("=", 80), "\n", sep = "", file = summary_file, append = TRUE)
  cat("✓ ENHANCED ANALYSIS COMPLETED SUCCESSFULLY!\n", file = summary_file, append = TRUE)
  cat("=" , rep("=", 80), "\n", sep = "", file = summary_file, append = TRUE)
  
  log_event(sprintf("Analysis complete! Results saved to %s", out_dir), "SUCCESS")
  log_event(sprintf("Summary report: %s", summary_file), "SUCCESS")
  
} else {
  log_event("ERROR: No results generated - check input files and data structure", "ERROR")
}

plan(sequential)
cat("\n✓ ENHANCED ANALYSIS WITH COMPREHENSIVE DIAGNOSTICS COMPLETED!\n")
cat(sprintf("\nCheck these directories:\n"))
cat(sprintf("  - Diagnostics: %s\n", diag_dir))
cat(sprintf("  - Maps: %s\n", maps_dir))
cat(sprintf("  - Reports: %s\n", reports_dir))