# Required libraries
library(ncdf4)
library(terra)
library(tidyverse)
library(Kendall)  # For Mann-Kendall test and Sen's slope
library(foreach)
library(doParallel)
library(broom)

# Set up parallel processing
num_cores <- detectCores() - 1
if (num_cores < 1) num_cores <- 1
registerDoParallel(cores = num_cores)
cat("Using", num_cores, "cores for parallel processing\n")

# Function to perform Trend-Free Pre-Whitening (TFPW)
perform_tfpw <- function(time_series, alpha = 0.05) {
  # Remove NA values
  ts_clean <- na.omit(time_series)
  n <- length(ts_clean)
  
  if (n < 10) {
    return(list(
      corrected_series = time_series,
      original_series = time_series,
      trend_removed = FALSE,
      rho = NA,
      p_value_rho = NA,
      sen_slope = NA,
      trend_line = rep(NA, n)
    ))
  }
  
  # Step 1: Calculate Sen's slope (non-parametric trend estimator)
  sen_result <- sens.slope(ts_clean)
  sen_slope <- sen_result$tau
  
  # Step 2: Remove trend to get detrended series
  time_index <- 1:n
  trend_line <- sen_slope * time_index
  detrended_series <- ts_clean - trend_line
  
  # Step 3: Calculate lag-1 autocorrelation on detrended series
  acf_result <- tryCatch({
    acf(detrended_series, lag.max = 1, plot = FALSE, na.action = na.pass)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(acf_result)) {
    rho <- NA
    p_value_rho <- NA
  } else {
    rho <- acf_result$acf[2]
    
    # Calculate variance of autocorrelation
    rho_var <- (1 + 2 * sum(acf(detrended_series, lag.max = 10, plot = FALSE)$acf[-1]^2, na.rm = TRUE)) / n
    
    # Test significance of autocorrelation
    z_rho <- rho / sqrt(rho_var)
    p_value_rho <- 2 * pnorm(-abs(z_rho))
  }
  
  # Step 4: Apply pre-whitening if autocorrelation is significant
  if (!is.na(rho) && !is.na(p_value_rho) && p_value_rho < alpha && abs(rho) > 0.1) {
    # Pre-whiten the detrended series
    prewhitened_detrended <- numeric(n)
    prewhitened_detrended[1] <- detrended_series[1]  # First value unchanged
    
    for (i in 2:n) {
      prewhitened_detrended[i] <- detrended_series[i] - rho * detrended_series[i-1]
    }
    
    # Step 5: Add back the original trend
    corrected_series <- prewhitened_detrended + trend_line
    
    trend_removed <- TRUE
  } else {
    # No significant autocorrelation, use original series
    corrected_series <- ts_clean
    trend_removed <- FALSE
  }
  
  return(list(
    corrected_series = corrected_series,
    original_series = ts_clean,
    trend_removed = trend_removed,
    rho = rho,
    p_value_rho = p_value_rho,
    sen_slope = sen_slope,
    trend_line = trend_line
  ))
}

# Function to perform spectral analysis with Monte Carlo confidence limits
perform_spectral_analysis <- function(time_series, n_sim = 1000, alpha = 0.05) {
  ts_clean <- na.omit(time_series)
  n <- length(ts_clean)
  
  if (n < 20) {
    return(list(
      frequencies = NA,
      spectral_density = NA,
      significant_peaks = NA,
      confidence_limits = NA,
      peak_frequencies = NA,
      peak_periods = NA
    ))
  }
  
  # Perform FFT
  fft_result <- fft(ts_clean - mean(ts_clean))
  spectral_density <- Mod(fft_result[1:(n/2)])^2 / n
  
  # Frequency vector
  frequencies <- seq(0, 0.5, length.out = n/2)
  
  # Monte Carlo simulation for confidence limits
  max_spectra <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Generate random series with same length and variance
    random_series <- rnorm(n, mean = mean(ts_clean), sd = sd(ts_clean))
    fft_rand <- fft(random_series - mean(random_series))
    spectral_rand <- Mod(fft_rand[1:(n/2)])^2 / n
    max_spectra[i] <- max(spectral_rand, na.rm = TRUE)
  }
  
  # Calculate confidence limits (95th percentile)
  conf_limit <- quantile(max_spectra, 1 - alpha, na.rm = TRUE)
  
  # Identify significant peaks
  significant_peaks <- spectral_density > conf_limit
  peak_indices <- which(significant_peaks & !is.na(significant_peaks))
  
  peak_frequencies <- frequencies[peak_indices]
  peak_periods <- ifelse(peak_frequencies > 0, 1/peak_frequencies, NA)
  
  return(list(
    frequencies = frequencies,
    spectral_density = spectral_density,
    significant_peaks = significant_peaks,
    confidence_limits = conf_limit,
    peak_frequencies = peak_frequencies,
    peak_periods = peak_periods,
    n_simulations = n_sim
  ))
}

# Function to process a single grid point for annual data
process_annual_point <- function(data_array, lat_idx, lon_idx, n_years) {
  tryCatch({
    # Extract time series for this grid point
    ts_monthly <- data_array[, lat_idx, lon_idx]
    
    # Aggregate to annual values (mean for PET, sum for precipitation - you can modify this)
    ts_annual <- numeric(n_years)
    
    for (i in 1:n_years) {
      year_data <- ts_monthly[((i-1)*12 + 1):(i*12)]
      if (sum(!is.na(year_data)) >= 6) {  # At least 6 months of data
        ts_annual[i] <- mean(year_data, na.rm = TRUE)  # Use mean for PET, sum for precipitation
      } else {
        ts_annual[i] <- NA
      }
    }
    
    # Apply TFPW
    tfpw_result <- perform_tfpw(ts_annual)
    
    # Perform Mann-Kendall test on pre-whitened series
    mk_result <- if (!all(is.na(tfpw_result$corrected_series))) {
      MannKendall(tfpw_result$corrected_series)
    } else {
      list(tau = NA, sl = NA, S = NA, varS = NA)
    }
    
    # Perform spectral analysis on pre-whitened series
    spectral_result <- if (!all(is.na(tfpw_result$corrected_series))) {
      perform_spectral_analysis(tfpw_result$corrected_series)
    } else {
      list(peak_periods = NA, significant_peaks = NA)
    }
    
    # Count significant spectral peaks
    n_peaks <- if (!all(is.na(spectral_result$significant_peaks))) {
      sum(spectral_result$significant_peaks, na.rm = TRUE)
    } else {
      0
    }
    
    return(data.frame(
      lat = lat_coords[lat_idx],
      lon = lon_coords[lon_idx],
      period_type = "annual",
      month = NA,
      tau = ifelse(!is.null(mk_result$tau), mk_result$tau, NA),
      p_value_mk = ifelse(!is.null(mk_result$sl), mk_result$sl, NA),
      n = length(na.omit(ts_annual)),
      autocorrelation_corrected = tfpw_result$trend_removed,
      rho = tfpw_result$rho,
      p_value_rho = tfpw_result$p_value_rho,
      sen_slope = tfpw_result$sen_slope,
      n_spectral_peaks = n_peaks,
      dominant_period = if (length(spectral_result$peak_periods) > 0) {
        if (!all(is.na(spectral_result$peak_periods))) {
          spectral_result$peak_periods[1]
        } else NA
      } else NA,
      spectral_confidence = if (!is.null(spectral_result$confidence_limits)) {
        spectral_result$confidence_limits
      } else NA
    ))
  }, error = function(e) {
    return(data.frame(
      lat = lat_coords[lat_idx],
      lon = lon_coords[lon_idx],
      period_type = "annual",
      month = NA,
      tau = NA,
      p_value_mk = NA,
      n = 0,
      autocorrelation_corrected = FALSE,
      rho = NA,
      p_value_rho = NA,
      sen_slope = NA,
      n_spectral_peaks = 0,
      dominant_period = NA,
      spectral_confidence = NA
    ))
  })
}

# Function to process a single grid point for a specific calendar month
process_monthly_point <- function(data_array, lat_idx, lon_idx, month_num, n_years) {
  tryCatch({
    # Extract time series for this grid point
    ts_monthly <- data_array[, lat_idx, lon_idx]
    
    # Extract data for the specific calendar month across all years
    month_indices <- seq(month_num, length(ts_monthly), by = 12)
    ts_specific_month <- ts_monthly[month_indices]
    
    # Apply TFPW
    tfpw_result <- perform_tfpw(ts_specific_month)
    
    # Perform Mann-Kendall test on pre-whitened series
    mk_result <- if (!all(is.na(tfpw_result$corrected_series))) {
      MannKendall(tfpw_result$corrected_series)
    } else {
      list(tau = NA, sl = NA, S = NA, varS = NA)
    }
    
    # Perform spectral analysis on pre-whitened series
    spectral_result <- if (!all(is.na(tfpw_result$corrected_series))) {
      perform_spectral_analysis(tfpw_result$corrected_series)
    } else {
      list(peak_periods = NA, significant_peaks = NA)
    }
    
    # Count significant spectral peaks
    n_peaks <- if (!all(is.na(spectral_result$significant_peaks))) {
      sum(spectral_result$significant_peaks, na.rm = TRUE)
    } else {
      0
    }
    
    return(data.frame(
      lat = lat_coords[lat_idx],
      lon = lon_coords[lon_idx],
      period_type = "monthly",
      month = month_num,
      tau = ifelse(!is.null(mk_result$tau), mk_result$tau, NA),
      p_value_mk = ifelse(!is.null(mk_result$sl), mk_result$sl, NA),
      n = length(na.omit(ts_specific_month)),
      autocorrelation_corrected = tfpw_result$trend_removed,
      rho = tfpw_result$rho,
      p_value_rho = tfpw_result$p_value_rho,
      sen_slope = tfpw_result$sen_slope,
      n_spectral_peaks = n_peaks,
      dominant_period = if (length(spectral_result$peak_periods) > 0) {
        if (!all(is.na(spectral_result$peak_periods))) {
          spectral_result$peak_periods[1]
        } else NA
      } else NA,
      spectral_confidence = if (!is.null(spectral_result$confidence_limits)) {
        spectral_result$confidence_limits
      } else NA
    ))
  }, error = function(e) {
    return(data.frame(
      lat = lat_coords[lat_idx],
      lon = lon_coords[lon_idx],
      period_type = "monthly",
      month = month_num,
      tau = NA,
      p_value_mk = NA,
      n = 0,
      autocorrelation_corrected = FALSE,
      rho = NA,
      p_value_rho = NA,
      sen_slope = NA,
      n_spectral_peaks = 0,
      dominant_period = NA,
      spectral_confidence = NA
    ))
  })
}

# Main function to process the entire NetCDF file
process_combined_analysis <- function(nc_file_path, var_name, is_precipitation = FALSE, 
                                      output_prefix = "combined_analysis", n_sim_spectral = 500) {
  cat("Reading NetCDF file:", nc_file_path, "\n")
  
  # Read the NetCDF file using terra
  r <- rast(nc_file_path, varname = var_name)
  cat("Data dimensions - nlyr:", nlyr(r), "nrow:", nrow(r), "ncol:", ncol(r), "\n")
  
  # Get coordinates
  lon_coords <<- xFromCol(r, 1:ncol(r))
  lat_coords <<- yFromRow(r, 1:nrow(r))
  
  # Get time information
  time_info <- r$times
  
  if (is.null(time_info)) {
    # Alternative way to get time dimension
    nc_file <- nc_open(nc_file_path)
    time_var <- ncvar_get(nc_file, "time")
    time_units <- ncatt_get(nc_file, "time", "units")$value
    nc_close(nc_file)
    
    # Convert time to POSIXct
    if (grepl("since", time_units)) {
      origin_time <- strsplit(time_units, " ")[[1]][3]
      time_info <- as.POSIXct((time_var * 86400) + as.numeric(as.POSIXct(origin_time, tz = "UTC")), 
                              origin = "1970-01-01", tz = "UTC")
    } else {
      time_info <- as.POSIXct(time_var, origin = "1970-01-01", tz = "UTC")
    }
  }
  
  # Extract years
  years <- as.numeric(format(time_info, "%Y"))
  unique_years <- unique(years)
  n_years <- length(unique_years)
  
  cat("Processing", n_years, "years of data from", min(unique_years), "to", max(unique_years), "\n")
  
  # Extract the data array
  data_array <- values(r)
  dim(data_array) <- c(nrow(r), ncol(r), nlyr(r))
  data_array <- aperm(data_array, c(3, 1, 2))  # Time, Lat, Lon
  
  cat("Starting annual data processing...\n")
  # Process annual data for all grid points
  annual_results <- foreach(i = 1:nrow(r), .combine = 'rbind', .packages = c('Kendall')) %dopar% {
    results_i <- foreach(j = 1:ncol(r), .combine = 'rbind') %do% {
      process_annual_point(data_array, i, j, n_years)
    }
    results_i
  }
  
  cat("Starting monthly data processing (all 12 calendar months)...\n")
  # Process monthly data for all grid points and all 12 months
  monthly_results <- foreach(month_num = 1:12, .combine = 'rbind', .packages = c('Kendall')) %dopar% {
    cat("Processing month", month_num, "\n")
    results_month <- foreach(i = 1:nrow(r), .combine = 'rbind') %do% {
      results_i <- foreach(j = 1:ncol(r), .combine = 'rbind') %do% {
        process_monthly_point(data_array, i, j, month_num, n_years)
      }
      results_i
    }
    results_month
  }
  
  # Combine results
  all_results <- bind_rows(annual_results, monthly_results)
  
  # Save results
  output_file <- paste0(output_prefix, "_combined_results.csv")
  write.csv(all_results, output_file, row.names = FALSE)
  cat("Combined results saved to:", output_file, "\n")
  
  # Create summary statistics
  summary_stats <- all_results %>%
    group_by(period_type, month) %>%
    summarise(
      total_points = n(),
      valid_points = sum(!is.na(tau), na.rm = TRUE),
      percent_valid = mean(!is.na(tau), na.rm = TRUE) * 100,
      significant_trends = sum(p_value_mk < 0.05, na.rm = TRUE),
      percent_significant_trends = mean(p_value_mk < 0.05, na.rm = TRUE) * 100,
      mean_tau = mean(tau, na.rm = TRUE),
      mean_sen_slope = mean(sen_slope, na.rm = TRUE),
      points_with_peaks = sum(n_spectral_peaks > 0, na.rm = TRUE),
      percent_with_peaks = mean(n_spectral_peaks > 0, na.rm = TRUE) * 100,
      mean_dominant_period = mean(dominant_period, na.rm = TRUE),
      .groups = 'drop'
    )
  
  summary_file <- paste0(output_prefix, "_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  cat("Summary statistics saved to:", summary_file, "\n")
  
  # Create spatial outputs
  if (!is.null(annual_results) && nrow(na.omit(annual_results)) > 0) {
    # Annual significant trends raster
    annual_trend <- annual_results %>%
      filter(!is.na(p_value_mk)) %>%
      mutate(significant_trend = as.numeric(p_value_mk < 0.05))
    
    if (nrow(annual_trend) > 0) {
      r_trend <- rast(r, nlyrs = 1)
      values(r_trend) <- NA
      xy <- cbind(annual_trend$lon, annual_trend$lat)
      cell_indices <- cellFromXY(r_trend, xy)
      values(r_trend)[cell_indices] <- annual_trend$significant_trend
      
      trend_file <- paste0(output_prefix, "_annual_significant_trends.tif")
      writeRaster(r_trend, trend_file, overwrite = TRUE)
      cat("Annual significant trends raster saved to:", trend_file, "\n")
    }
    
    # Annual spectral peaks raster
    annual_peaks <- annual_results %>%
      filter(!is.na(n_spectral_peaks)) %>%
      mutate(has_peaks = as.numeric(n_spectral_peaks > 0))
    
    if (nrow(annual_peaks) > 0) {
      r_peaks <- rast(r, nlyrs = 1)
      values(r_peaks) <- NA
      xy <- cbind(annual_peaks$lon, annual_peaks$lat)
      cell_indices <- cellFromXY(r_peaks, xy)
      values(r_peaks)[cell_indices] <- annual_peaks$has_peaks
      
      peaks_file <- paste0(output_prefix, "_annual_spectral_peaks.tif")
      writeRaster(r_peaks, peaks_file, overwrite = TRUE)
      cat("Annual spectral peaks raster saved to:", peaks_file, "\n")
    }
  }
  
  return(list(
    full_results = all_results,
    summary_stats = summary_stats
  ))
}

# Example usage - modify these parameters for your data
nc_file_path <- "your_data.nc"  # Replace with your NetCDF file path
var_name <- "pr"  # Replace with your variable name
is_precipitation <- TRUE  # TRUE for precipitation (sum annual), FALSE for PET (mean annual)
output_prefix <- "combined_tfpw_spectral"
n_sim_spectral <- 500  # Number of Monte Carlo simulations for spectral confidence limits

# Run the combined analysis
results <- process_combined_analysis(nc_file_path, var_name, is_precipitation, 
                                     output_prefix, n_sim_spectral)

# Stop parallel processing
stopImplicitCluster()

cat("Combined TFPW + Mann-Kendall + Spectral analysis completed successfully!\n")