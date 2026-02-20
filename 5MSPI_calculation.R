##############################################
# MSPI & MSPEI Calculation - BAZRAFSHAN 2014 METHOD 
# Computes both indices with .nc and .xlsx outputs
# PCA on SPI/SPEI scales AT SAME CALENDAR TIME
##############################################
library(terra)
library(stringr)
library(openxlsx)  # For Excel output

# Utility: Format seconds to HH:MM:SS
format_time <- function(seconds) {
  h <- floor(seconds / 3600)
  m <- floor((seconds %% 3600) / 60)
  s <- round(seconds %% 60)
  if (h > 0) sprintf("%02dh %02dm %02ds", h, m, s)
  else if (m > 0) sprintf("%02dm %02ds", m, s)
  else sprintf("%02ds", s)
}

#==============================================================================
# Load and Validate Basin Boundary
#==============================================================================
load_basin_boundary <- function(basin_path) {
  cat(sprintf("\n[%s] STEP 1: Loading basin boundary...\n",
              format(Sys.time(), "%H:%M:%S")))
  
  if (!file.exists(basin_path)) {
    stop(sprintf("ERROR: Basin boundary not found: %s",
                 normalizePath(basin_path)))
  }
  
  basin <- vect(basin_path)
  cat(sprintf("OK: Basin loaded: %s\n", basename(basin_path)))
  cat(sprintf("  -> Geometry: %s | Features: %d | CRS: %s\n",
              geomtype(basin), nrow(basin), crs(basin, describe = TRUE)$name))
  
  # Show basin extent + area
  cat(sprintf("  -> Basin extent: xmin=%.2f, xmax=%.2f, ymin=%.2f, ymax=%.2f\n",
              xmin(basin), xmax(basin), ymin(basin), ymax(basin)))
  
  if (geomtype(basin) == "polygons") {
    area_km2 <- expanse(basin) / 1e6
    cat(sprintf("  -> Basin area: %.1f km2\n", area_km2))
    
    # Show sample coordinates
    coords_sample <- crds(basin[1])
    cat(sprintf("  -> Sample boundary coordinates (first 5 vertices):\n"))
    for (i in 1:min(5, nrow(coords_sample))) {
      cat(sprintf("       (%.4f, %.4f)\n", coords_sample[i,1], coords_sample[i,2]))
    }
  }
  
  return(basin)
}

#==============================================================================
# Create Basin Mask
#==============================================================================
create_basin_mask <- function(ref_rast, basin) {
  cat(sprintf("\n[%s] STEP 2: Creating basin mask...\n",
              format(Sys.time(), "%H:%M:%S")))
  
  if (!same.crs(ref_rast, basin)) {
    cat(sprintf("  -> Reprojecting basin from %s to raster CRS...\n", 
                crs(basin, describe = TRUE)$name))
    basin <- project(basin, crs(ref_rast))
  }
  
  mask <- rasterize(basin, ref_rast, field = 1, background = NA, touches = TRUE)
  basin_pixels <- which(!is.na(values(mask)))
  n_basin <- length(basin_pixels)
  
  cat(sprintf("OK: Basin mask created\n"))
  cat(sprintf("  -> Basin pixels in raster domain: %s\n", format(n_basin, big.mark = ",")))
  cat(sprintf("  -> Basin coverage: %.2f%% of raster cells\n", 
              100 * n_basin / ncell(ref_rast)))
  
  return(list(
    mask = mask,
    basin_pixels = basin_pixels,
    n_basin = n_basin
  ))
}

#==============================================================================
# Reconstruct Full Time Series (improved: reads actual time dimension)
#==============================================================================
reconstruct_full_timeseries <- function(input_dir, index_type = c("SPI", "SPEI"),
                                        start_year = 1950, end_year = 2025) {
  cat(sprintf("\n[%s] STEP 3: Reconstructing %s time series...\n", 
              format(Sys.time(), "%H:%M:%S"), toupper(index_type)))
  
  index_type <- match.arg(index_type)
  
  pattern <- if (index_type == "SPI") 
    "spi_\\d{2}_month\\d{2}_[A-Za-z]{3}\\.nc$" 
  else 
    "spei_\\d{2}_month\\d{2}_[A-Za-z]{3}\\.nc$"
  
  all_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  all_files <- all_files[grepl("\\.nc$", all_files)]
  
  if (length(all_files) == 0) {
    stop(sprintf("ERROR: No %s files found matching pattern!", index_type))
  }
  
  cat(sprintf("  -> Found %d files matching pattern\n", length(all_files)))
  
  # Parse file metadata
  file_info <- data.frame(
    filepath = all_files,
    scale = as.integer(str_extract(basename(all_files), "(?<=spei_|spi_)\\d{2}")),
    month_num = as.integer(str_extract(basename(all_files), "(?<=month)\\d{2}")),
    stringsAsFactors = FALSE
  )
  file_info <- file_info[!is.na(file_info$scale) & !is.na(file_info$month_num), ]
  
  scales <- sort(unique(file_info$scale))
  months <- sort(unique(file_info$month_num))
  cat(sprintf("  -> Detected scales: %s\n", paste(scales, collapse = ", ")))
  cat(sprintf("  -> Detected months: %s\n", paste(months, collapse = ", ")))
  
  # Load reference raster (first file)
  ref_rast <- rast(file_info$filepath[1])
  n_pixels <- ncell(ref_rast)
  cat(sprintf("  -> Reference raster: %d layers, %dx%d cells (%s pixels)\n",
              nlyr(ref_rast), ncol(ref_rast), nrow(ref_rast), 
              format(n_pixels, big.mark = ",")))
  
  # Build global timeline from start_year to end_year
  # Each file covers ONE calendar month across ALL years (layers = years).
  # File spi_XX_monthMM_Mon.nc: layer k (1-indexed) = year (start_year + k - 1), month MM.
  # Global index formula: (year - start_year) * 12 + MM
  global_dates <- seq(as.Date(sprintf("%d-01-01", start_year)), 
                      by = "month", length.out = (end_year - start_year + 1) * 12)
  total_months <- length(global_dates)
  
  cat(sprintf("  -> Global timeline: %s to %s (%d months)\n",
              format(min(global_dates), "%Y-%m"), 
              format(max(global_dates), "%Y-%m"), total_months))
  
  # Initialize array
  ts_array <- array(NA, dim = c(n_pixels, total_months, length(scales)))
  dimnames(ts_array) <- list(NULL, NULL, sprintf("scale_%02d", scales))
  
  cat(sprintf("  -> Target array dimensions: %s pixels x %d months x %d scales\n",
              format(n_pixels, big.mark = ","), total_months, length(scales)))
  
  # Number of years spanned
  n_years <- end_year - start_year + 1
  
  # Populate array
  cat("  -> Populating time series array...\n")
  pb <- txtProgressBar(min = 0, max = length(scales) * 12, style = 3)
  counter <- 0
  
  for (sc in scales) {
    scale_idx <- which(scales == sc)
    for (mon in 1:12) {
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
      if (Sys.info()["sysname"] == "Windows") flush.console()
      
      idx <- which(file_info$scale == sc & file_info$month_num == mon)
      if (length(idx) == 0) next
      
      r <- rast(file_info$filepath[idx[1]])
      n_layers_actual <- nlyr(r)
      
      # Each layer k corresponds to year (start_year + k - 1) for this month.
      # Map directly: global_time_idx = (year - start_year) * 12 + mon
      for (layer_idx in 1:n_layers_actual) {
        year_for_layer <- start_year + layer_idx - 1
        if (year_for_layer > end_year) next
        global_time_idx <- (year_for_layer - start_year) * 12 + mon
        if (global_time_idx < 1 || global_time_idx > total_months) next
        
        vals <- values(r[[layer_idx]], mat = FALSE)
        if (length(vals) == n_pixels) {
          ts_array[, global_time_idx, scale_idx] <- vals
        }
      }
    }
  }
  close(pb)
  
  # Validate completeness
  total_valid <- sum(!is.na(ts_array))
  total_cells <- prod(dim(ts_array))
  valid_pct_total <- 100 * total_valid / total_cells
  cat(sprintf("\n  -> Data completeness: %.2f%% (%s / %s values)\n",
              valid_pct_total,
              format(total_valid, big.mark = ","),
              format(total_cells, big.mark = ",")))
  
  if (valid_pct_total < 1.0) {
    stop(sprintf("\nERROR: Time series array is %.2f%% NA values!\n", 
                 100 - valid_pct_total))
  }
  
  return(list(
    ts_array = ts_array,
    dates = global_dates,
    ref_rast = ref_rast,
    scales = scales
  ))
}

#==============================================================================
# Validate Basin Data
#==============================================================================
validate_basin_data <- function(ts_array, basin_pixels) {
  cat(sprintf("\n[%s] STEP 4: Validating basin data quality...\n", format(Sys.time(), "%H:%M:%S")))
  
  n_basin <- length(basin_pixels)
  cat(sprintf("  -> Analyzing %s basin pixels...\n", format(n_basin, big.mark = ",")))
  
  valid_counts <- rowSums(!is.na(ts_array[basin_pixels, , , drop = FALSE]), na.rm = TRUE)
  total_per_pixel <- prod(dim(ts_array)[2:3])
  basin_data_quality <- valid_counts / total_per_pixel
  
  cat(sprintf("  -> Data quality stats:\n"))
  cat(sprintf("      Mean: %.1f%% | Median: %.1f%% | Min: %.1f%% | Max: %.1f%%\n",
              100 * mean(basin_data_quality),
              100 * median(basin_data_quality),
              100 * min(basin_data_quality),
              100 * max(basin_data_quality)))
  
  sufficient_data <- basin_data_quality >= 0.05
  basin_pixels_filtered <- basin_pixels[sufficient_data]
  n_filtered <- length(basin_pixels_filtered)
  
  cat(sprintf("\nOK: Proceeding with %s valid basin pixels (%.1f%% of basin)\n",
              format(n_filtered, big.mark = ","), 100*n_filtered/n_basin))
  
  if (n_filtered == 0) {
    stop("ERROR: NO valid pixels in basin!")
  }
  
  return(list(
    basin_pixels_valid = basin_pixels_filtered,
    n_valid = n_filtered,
    data_quality = basin_data_quality[sufficient_data]
  ))
}

#==============================================================================
# MSPI/MSPEI COMPUTATION - BAZRAFSHAN 2014 METHOD
# (improved missing‑data handling: drop scales with >50% missing)
#==============================================================================
# FUNCTION: Compute MSPI with fixed loadings
#==============================================================================
# MSPI/MSPEI COMPUTATION - BAZRAFSHAN 2014 METHOD (
# Each month-year uses its OWN 12 SPI/SPEI scales
# First valid value: December 1950 (index 12)
#==============================================================================
#==============================================================================
# MSPI/MSPEI COMPUTATION - BAZRAFSHAN 2014 METHOD (CORRECTED)
# FIXED: ONE PCA from ALL month-years, first valid = December 1950 (index 12)
#==============================================================================
#==============================================================================
# MSPI/MSPEI COMPUTATION - BAZRAFSHAN 2014 METHOD (CORRECTED)
# FIXED: Requires ALL 12 scales (1-12), first valid = December 1950 (index 12)
#==============================================================================
compute_mspi_bazrafshan <- function(index_tensor, dates, basin_pixels, scales) {
  n_pixels_total <- dim(index_tensor)[1]
  n_time <- dim(index_tensor)[2]
  n_scales <- length(scales)
  
  mspi_matrix <- matrix(NA, nrow = n_pixels_total, ncol = n_time)
  var_explained_matrix <- matrix(NA, nrow = n_pixels_total, ncol = 1)
  
  n_basin <- length(basin_pixels)
  
  cat(sprintf("\n[%s] STEP 5: Computing MSPI/MSPEI for %s pixels...\n", 
              format(Sys.time(), "%H:%M:%S"), format(n_basin, big.mark = ",")))
  cat("  -> Method: ONE PCA from ALL month-years (Bazrafshan 2014)\n")
  cat("  -> REQUIRED: All 12 scales (1-12) must be present\n")
  cat("  -> First valid value: December 1950 (index 12)\n")
  
  # === Filter to only scales 1-12 for PCA ===
  pca_scales_idx <- which(scales >= 1 & scales <= 12)
  pca_scales <- scales[pca_scales_idx]
  n_scales_pca <- length(pca_scales_idx)
  
  if (n_scales_pca != 12) {
    warning(sprintf("Expected 12 scales (1-12), found %d", n_scales_pca))
  }
  
  pb <- txtProgressBar(min = 0, max = n_basin, style = 3)
  stats <- list(successful = 0, failed = 0)
  
  for (i in 1:n_basin) {
    setTxtProgressBar(pb, i)
    if (Sys.info()["sysname"] == "Windows") flush.console()
    
    px <- basin_pixels[i]
    pixel_data <- matrix(index_tensor[px, , ], nrow = n_time, ncol = n_scales)
    
    # === Use ONLY scales 1-12 ===
    pixel_data_pca <- pixel_data[, pca_scales_idx, drop = FALSE]
    
    # === STEP 1: Build training matrix - REQUIRE ALL 12 SCALES ===
    # CRITICAL FIX: Require ALL 12 scales present (not 50%)
    valid_rows <- rowSums(!is.na(pixel_data_pca)) == n_scales_pca
    
    if (sum(valid_rows) < 12) next
    
    training_data <- pixel_data_pca[valid_rows, , drop = FALSE]
    
    # Center (don't scale per Bazrafshan)
    training_centered <- scale(training_data, center = TRUE, scale = FALSE)
    training_complete <- na.omit(training_centered)
    
    if (nrow(training_complete) < 12) next
    
    # === STEP 2: Run ONE PCA to get loadings ===
    pca_result <- tryCatch({
      prcomp(training_complete, center = FALSE, scale. = FALSE)
    }, error = function(e) NULL)
    
    if (is.null(pca_result) || length(pca_result$sdev) == 0) next
    
    pc1_loadings <- pca_result$rotation[, 1]
    training_means <- colMeans(training_data)
    
    total_var <- sum(pca_result$sdev^2)
    var_explained_matrix[px, 1] <- pca_result$sdev[1]^2 / total_var
    
    # === STEP 3: Apply loadings - REQUIRE ALL 12 SCALES for each month ===
    for (t in 1:n_time) {
      month_vals <- pixel_data_pca[t, ]
      # CRITICAL: Require ALL 12 scales present (not 50%)
      if (sum(!is.na(month_vals)) == n_scales_pca) {
        centered_vals <- month_vals - training_means
        mspi_matrix[px, t] <- sum(centered_vals * pc1_loadings, na.rm = TRUE)
      }
    }
    
    # === STEP 4: Standardize final MSPI series ===
    mspi_vals <- mspi_matrix[px, ]
    valid_idx <- !is.na(mspi_vals)
    
    if (sum(valid_idx) > 1) {
      mspi_std <- (mspi_vals[valid_idx] - mean(mspi_vals[valid_idx])) / sd(mspi_vals[valid_idx])
      mspi_matrix[px, valid_idx] <- mspi_std
    }
    
    if (sum(!is.na(mspi_matrix[px, ])) > 0) stats$successful <- stats$successful + 1
    else stats$failed <- stats$failed + 1
  }
  close(pb)
  
  cat(sprintf("\nOK: Completed %d/%d pixels successful (%.1f%%)\n",
              stats$successful, n_basin, 100 * stats$successful / n_basin))
  
  return(list(
    mspi = mspi_matrix,
    variance_explained = var_explained_matrix,
    stats = stats
  ))
}
#==============================================================================
export_basin_timeseries_excel <- function(mspi_rast, basin_mask, dates, 
                                          output_file, index_name = "MSPI") {
  cat(sprintf("\n[%s] Exporting basin time series to Excel...\n", 
              format(Sys.time(), "%H:%M:%S")))
  
  # Extract basin pixels
  basin_pixels <- which(!is.na(values(basin_mask)))
  n_basin <- length(basin_pixels)
  n_time <- nlyr(mspi_rast)
  
  cat(sprintf("  -> Extracting %d time steps for %s basin pixels...\n", 
              n_time, format(n_basin, big.mark = ",")))
  
  # Extract values for basin pixels
  basin_values <- values(mspi_rast)[basin_pixels, , drop = FALSE]
  
  # --- Missingness report ---
  all_na_steps <- which(colSums(!is.na(basin_values)) == 0)
  n_all_na <- length(all_na_steps)
  pct_all_na <- 100 * n_all_na / n_time
  cat(sprintf("  -> Missingness report: %d / %d time steps (%.1f%%) have ALL basin pixels NA\n",
              n_all_na, n_time, pct_all_na))
  if (n_all_na > 0) {
    cat(sprintf("     Example all‑NA step indices: %s\n", 
                paste(head(all_na_steps, 5), collapse = ", ")))
  }
  
  # Compute basin-wide statistics (handle all‑NA columns gracefully)
  basin_mean <- colMeans(basin_values, na.rm = TRUE)
  basin_median <- apply(basin_values, 2, median, na.rm = TRUE)
  # For min, max, sd: return NA if all values are NA
  basin_min <- apply(basin_values, 2, function(x) if(all(is.na(x))) NA else min(x, na.rm = TRUE))
  basin_max <- apply(basin_values, 2, function(x) if(all(is.na(x))) NA else max(x, na.rm = TRUE))
  basin_sd  <- apply(basin_values, 2, function(x) if(all(is.na(x))) NA else sd(x, na.rm = TRUE))
  
  # Count valid pixels per time step
  valid_count <- colSums(!is.na(basin_values))
  
  # Create data frame
  df <- data.frame(
    Date = format(dates, "%Y-%m"),
    Year = as.numeric(format(dates, "%Y")),
    Month = as.numeric(format(dates, "%m")),
    Basin_Mean = round(basin_mean, 3),
    Basin_Median = round(basin_median, 3),
    Basin_Min = round(basin_min, 3),
    Basin_Max = round(basin_max, 3),
    Basin_SD = round(basin_sd, 3),
    Valid_Pixels = valid_count,
    Total_Pixels = n_basin,
    Coverage_Pct = round(100 * valid_count / n_basin, 1)
  )
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add time series sheet
  addWorksheet(wb, "Basin_Timeseries")
  writeData(wb, "Basin_Timeseries", df, startRow = 1)
  
  # Format header
  headerStyle <- createStyle(
    fontSize = 11,
    fontColour = "#FFFFFF",
    halign = "center",
    fgFill = "#4F81BD",
    border = "TopBottomLeftRight",
    textDecoration = "bold"
  )
  addStyle(wb, "Basin_Timeseries", headerStyle, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  
  # Add metadata sheet
  addWorksheet(wb, "Metadata")
  metadata <- data.frame(
    Parameter = c("Index Type", "Basin Area (pixels)", "Time Period", 
                  "Start Date", "End Date", "Total Time Steps",
                  "Computation Method", "Reference"),
    Value = c(index_name, n_basin, 
              sprintf("%s to %s", format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m")),
              format(min(dates), "%Y-%m-%d"), format(max(dates), "%Y-%m-%d"), n_time,
              "Bazrafshan 2014 PCA method",
              "Bazrafshan et al. (2014) Water Resources Management")
  )
  writeData(wb, "Metadata", metadata)
  addStyle(wb, "Metadata", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
  
  # Save workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat(sprintf("OK: Excel file saved: %s\n", basename(output_file)))
  cat(sprintf("  -> Sheets: Basin_Timeseries, Metadata\n"))
  cat(sprintf("  -> Time steps: %d | Basin pixels: %s\n", n_time, format(n_basin, big.mark = ",")))
  
  return(invisible(df))
}
#==============================================================================
# Diagnostic: MSPI/MSPEI Missingness Report
#==============================================================================
diagnose_mspi_missingness <- function(mspi_matrix, basin_pixels, dates, index_name) {
  # Subset to basin pixels only
  basin_mspi <- mspi_matrix[basin_pixels, , drop = FALSE]
  n_time <- ncol(basin_mspi)
  n_pixels <- nrow(basin_mspi)
  
  # Valid pixels per time step
  valid_per_time <- colSums(!is.na(basin_mspi))
  first_valid <- which(valid_per_time > 0)[1]
  last_valid <- tail(which(valid_per_time > 0), 1)
  
  cat(sprintf("\n[%s] MISSINGNESS DIAGNOSTIC for %s:\n", format(Sys.time(), "%H:%M:%S"), index_name))
  cat(sprintf("  Basin pixels: %d\n", n_pixels))
  cat(sprintf("  Time steps: %d\n", n_time))
  cat(sprintf("  First time step with any data: %s (index %d)\n", 
              ifelse(is.na(first_valid), "none", format(dates[first_valid], "%Y-%m")), first_valid))
  cat(sprintf("  Last time step with any data: %s (index %d)\n", 
              ifelse(is.na(last_valid), "none", format(dates[last_valid], "%Y-%m")), last_valid))
  
  # Time steps with zero valid pixels
  zero_steps <- which(valid_per_time == 0)
  n_zero <- length(zero_steps)
  cat(sprintf("  Time steps with zero valid pixels: %d (%.1f%%)\n", n_zero, 100 * n_zero / n_time))
  if (n_zero > 0 && n_zero <= 20) {
    cat("    Indices: ", paste(zero_steps, collapse = ", "), "\n")
    cat("    Dates: ", paste(format(dates[zero_steps], "%Y-%m"), collapse = ", "), "\n")
  } else if (n_zero > 20) {
    cat("    First 10 zero-step indices: ", paste(head(zero_steps, 10), collapse = ", "), "\n")
    cat("    Corresponding dates: ", paste(format(dates[head(zero_steps, 10)], "%Y-%m"), collapse = ", "), "\n")
  }
  
  # Data availability by month-of-year
  month_of_year <- as.numeric(format(dates, "%m"))
  years <- as.numeric(format(dates, "%Y"))
  cat("\n  Data availability by month-of-year:\n")
  for (mon in 1:12) {
    time_idx <- which(month_of_year == mon)
    valid_in_month <- sapply(time_idx, function(t) valid_per_time[t] > 0)
    n_years_with_data <- sum(valid_in_month)
    total_years <- length(time_idx)
    cat(sprintf("    Month %02d: %d/%d years have data (%.1f%%)\n", 
                mon, n_years_with_data, total_years, 100 * n_years_with_data / total_years))
    if (n_years_with_data < total_years) {
      missing_years <- years[time_idx][!valid_in_month]
      cat(sprintf("      Missing years: %s\n", paste(missing_years, collapse = ", ")))
    }
  }
  
  invisible(list(valid_per_time = valid_per_time, zero_steps = zero_steps))
}
#==============================================================================
# Detect Drought Events from Basin-Averaged Time Series
# (Authoritative definition shared with visualization script)
# Onset: value < onset_threshold; Termination: value >= termination_threshold
# Severity based on minimum value during event
#==============================================================================
detect_drought_events <- function(df, onset_threshold = -1.0,
                                  termination_threshold = 0.0,
                                  min_duration = 2,
                                  severe_threshold = -1.3,
                                  extreme_threshold = -1.6,
                                  exceptional_threshold = -2.0) {
  df     <- df[order(df$Date_formatted), ]
  values <- df$Basin_Mean
  dates  <- df$Date_formatted
  
  events <- data.frame(
    event_id       = integer(),
    start_date     = as.Date(character()),
    end_date       = as.Date(character()),
    duration_months = integer(),
    min_value      = numeric(),
    severity       = character(),
    stringsAsFactors = FALSE
  )
  
  event_id <- 0
  i <- 1
  
  while (i <= length(values)) {
    if (!is.na(values[i]) && values[i] < onset_threshold) {
      event_id  <- event_id + 1
      start_idx <- i
      min_val   <- values[i]
      
      j <- i + 1
      while (j <= length(values) && !is.na(values[j]) &&
             values[j] < termination_threshold) {
        min_val <- min(min_val, values[j], na.rm = TRUE)
        j <- j + 1
      }
      
      end_idx  <- j - 1
      duration <- end_idx - start_idx + 1
      
      if (duration >= min_duration) {
        severity <- if (min_val < exceptional_threshold) "Exceptional" else
          if (min_val < extreme_threshold)     "Extreme"     else
            if (min_val < severe_threshold)      "Severe"      else
              "Moderate"
        
        events <- rbind(events, data.frame(
          event_id        = event_id,
          start_date      = dates[start_idx],
          end_date        = dates[end_idx],
          duration_months = duration,
          min_value       = round(min_val, 3),
          severity        = severity,
          stringsAsFactors = FALSE
        ))
      }
      i <- end_idx + 1
    } else {
      i <- i + 1
    }
  }
  
  return(events)
}

#==============================================================================
# Export Combined Summary Statistics Excel
# Consolidates basin stats + drought events for both indices into one workbook.
# Called once at the end of compute_both_indices() so downstream users do not
# need to run the visualization script to obtain the drought catalogue.
#==============================================================================
export_summary_excel <- function(mspi_df, mspei_df, output_file,
                                 onset_threshold     = -1.0,
                                 termination_threshold = 0.0,
                                 min_duration        = 2) {
  cat(sprintf("\n[%s] Exporting summary statistics Excel...\n",
              format(Sys.time(), "%H:%M:%S")))
  
  # --- Add Date_formatted column (Date object) if not already present.
  #     export_basin_timeseries_excel() stores dates as "YYYY-MM" strings in
  #     the 'Date' column; detect_drought_events() needs actual Date objects
  #     in 'Date_formatted'. ---
  if (!"Date_formatted" %in% names(mspi_df)) {
    mspi_df$Date_formatted  <- as.Date(paste0(mspi_df$Date,  "-01"))
  }
  if (!"Date_formatted" %in% names(mspei_df)) {
    mspei_df$Date_formatted <- as.Date(paste0(mspei_df$Date, "-01"))
  }
  
  # --- Drought detection ---
  mspi_droughts_full   <- detect_drought_events(mspi_df,  onset_threshold,
                                                termination_threshold, min_duration)
  mspei_droughts_full  <- detect_drought_events(mspei_df, onset_threshold,
                                                termination_threshold, min_duration)
  
  mspi_recent  <- mspi_df [mspi_df $Date_formatted >= as.Date("2020-01-01"), ]
  mspei_recent <- mspei_df[mspei_df$Date_formatted >= as.Date("2020-01-01"), ]
  mspi_droughts_recent  <- detect_drought_events(mspi_recent,  onset_threshold,
                                                 termination_threshold, min_duration)
  mspei_droughts_recent <- detect_drought_events(mspei_recent, onset_threshold,
                                                 termination_threshold, min_duration)
  
  # --- Summary row ---
  stats_data <- data.frame(
    Index  = c("MSPI", "MSPEI"),
    Mean   = round(c(mean(mspi_df$Basin_Mean,  na.rm = TRUE),
                     mean(mspei_df$Basin_Mean, na.rm = TRUE)), 4),
    Median = round(c(median(mspi_df$Basin_Mean,  na.rm = TRUE),
                     median(mspei_df$Basin_Mean, na.rm = TRUE)), 4),
    StdDev = round(c(sd(mspi_df$Basin_Mean,  na.rm = TRUE),
                     sd(mspei_df$Basin_Mean, na.rm = TRUE)), 4),
    Min    = round(c(min(mspi_df$Basin_Mean,  na.rm = TRUE),
                     min(mspei_df$Basin_Mean, na.rm = TRUE)), 4),
    Max    = round(c(max(mspi_df$Basin_Mean,  na.rm = TRUE),
                     max(mspei_df$Basin_Mean, na.rm = TRUE)), 4),
    Drought_Events_Full_Period  = c(nrow(mspi_droughts_full),
                                    nrow(mspei_droughts_full)),
    Drought_Events_2020_Present = c(nrow(mspi_droughts_recent),
                                    nrow(mspei_droughts_recent))
  )
  
  headerStyle <- createStyle(
    fontSize = 11, fontColour = "#FFFFFF", halign = "center",
    fgFill = "#4F81BD", border = "TopBottomLeftRight", textDecoration = "bold"
  )
  
  wb <- createWorkbook()
  
  # Sheet 1: summary
  addWorksheet(wb, "Summary_Statistics")
  writeData(wb, "Summary_Statistics", stats_data, rowNames = FALSE)
  addStyle(wb, "Summary_Statistics", headerStyle,
           rows = 1, cols = 1:ncol(stats_data), gridExpand = TRUE)
  
  # Sheets 2-5: full-period and recent drought catalogues
  add_drought_sheet <- function(wb, sheet_name, df_events) {
    if (nrow(df_events) > 0) {
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, df_events)
      addStyle(wb, sheet_name, headerStyle,
               rows = 1, cols = 1:ncol(df_events), gridExpand = TRUE)
    }
  }
  
  add_drought_sheet(wb, "MSPI_Droughts_Full",        mspi_droughts_full)
  add_drought_sheet(wb, "MSPEI_Droughts_Full",       mspei_droughts_full)
  add_drought_sheet(wb, "MSPI_Droughts_2020_Present",  mspi_droughts_recent)
  add_drought_sheet(wb, "MSPEI_Droughts_2020_Present", mspei_droughts_recent)
  
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat(sprintf("OK: Summary Excel saved: %s\n", basename(output_file)))
  cat(sprintf("  -> MSPI  full-period drought events: %d\n", nrow(mspi_droughts_full)))
  cat(sprintf("  -> MSPEI full-period drought events: %d\n", nrow(mspei_droughts_full)))
  
  return(invisible(list(
    mspi_droughts_full   = mspi_droughts_full,
    mspei_droughts_full  = mspei_droughts_full,
    mspi_droughts_recent  = mspi_droughts_recent,
    mspei_droughts_recent = mspei_droughts_recent
  )))
}

#==============================================================================
# Main Function - Compute Both MSPI and MSPEI
#==============================================================================
compute_both_indices <- function(start_year = 1950,
                                 end_year = 2025,
                                 base_dir = ".",
                                 basin_path = NULL) {
  
  overall_start <- Sys.time()
  
  cat(sprintf("\n%s\n", paste(rep("=", 80), collapse = "")))
  cat(sprintf("MSPI & MSPEI CALCULATION - BAZRAFSHAN 2014 METHOD\n"))
  cat(sprintf("   Period: %d-%d | Basin: %s\n", 
              start_year, end_year, 
              ifelse(is.null(basin_path), "N/A", basename(basin_path))))
  cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
  
  if (is.null(basin_path) || !file.exists(basin_path)) {
    stop("ERROR: basin_path is required and must exist!")
  }
  
  # Load basin boundary (shared for both indices)
  basin <- load_basin_boundary(basin_path)
  
  results <- list()
  
  # ===== PROCESS MSPI =====
  cat(sprintf("\n%s\n", paste(rep("-", 80), collapse = "")))
  cat(sprintf("COMPUTING MSPI (Standardized Precipitation Index)\n"))
  cat(sprintf("%s\n", paste(rep("-", 80), collapse = "")))
  
  spi_start <- Sys.time()
  
  # Reconstruct SPI time series
  spi_data <- reconstruct_full_timeseries(
    input_dir = file.path(base_dir, "spi_results_seasonal"),
    index_type = "SPI",
    start_year = start_year,
    end_year = end_year
  )
  
  # Create mask
  mask_info <- create_basin_mask(spi_data$ref_rast, basin)
  
  # Validate data
  validation_spi <- validate_basin_data(spi_data$ts_array, mask_info$basin_pixels)
  
  # Compute MSPI
  mspi_result <- compute_mspi_bazrafshan(
    index_tensor = spi_data$ts_array,
    dates = spi_data$dates,
    basin_pixels = validation_spi$basin_pixels_valid,
    scales = spi_data$scales
  )
  
  # Diagnostic: MSPI missingness
  diagnose_mspi_missingness(mspi_result$mspi, 
                            validation_spi$basin_pixels_valid, 
                            spi_data$dates, 
                            "MSPI")
  # Assemble raster
  cat(sprintf("\n[%s] Assembling MSPI raster stack...\n", format(Sys.time(), "%H:%M:%S")))
  mspi_rast <- rast(spi_data$ref_rast[[1]])
  mspi_rast <- rep(mspi_rast, length(spi_data$dates))
  values(mspi_rast) <- mspi_result$mspi
  terra::time(mspi_rast) <- spi_data$dates
  
  # Save MSPI outputs
  output_dir_spi <- file.path(base_dir, "mspi_results")
  if (!dir.exists(output_dir_spi)) dir.create(output_dir_spi, recursive = TRUE)
  
  cat(sprintf("\n[%s] Saving MSPI outputs...\n", format(Sys.time(), "%H:%M:%S")))
  
  # NetCDF output
  nc_file_spi <- file.path(output_dir_spi, sprintf("mspi_monthly_%d-%d.nc", start_year, end_year))
  writeCDF(mspi_rast, nc_file_spi,
           varname = "mspi",
           longname = "Multivariate Standardized Precipitation Index",
           unit = "standardized_index",
           missval = -9999,
           overwrite = TRUE)
  cat(sprintf("OK: MSPI NetCDF saved: %s\n", basename(nc_file_spi)))
  
  # Excel output
  xlsx_file_spi <- file.path(output_dir_spi, sprintf("mspi_basin_timeseries_%d-%d.xlsx", 
                                                     start_year, end_year))
  mspi_df <- export_basin_timeseries_excel(
    mspi_rast = mspi_rast,
    basin_mask = mask_info$mask,
    dates = spi_data$dates,
    output_file = xlsx_file_spi,
    index_name = "MSPI"
  )
  
  # Variance explained maps
  var_rast_spi <- rast(spi_data$ref_rast[[1]])
  band_vals <- mspi_result$variance_explained[, 1] * 100
  band_vals[is.na(band_vals)] <- -9999
  values(var_rast_spi) <- band_vals
  names(var_rast_spi) <- "PC1_variance_explained"
  var_file_spi <- file.path(output_dir_spi, "mspi_variance_explained_pc1.nc")
  writeCDF(var_rast_spi, var_file_spi,
           varname  = "PC1_variance_explained",
           longname = "MSPI PC1 Variance Explained (%)",
           unit     = "percent",
           missval  = -9999,
           overwrite = TRUE)
  cat(sprintf("OK: MSPI variance map saved: %s\n", basename(var_file_spi)))
  
  spi_elapsed <- difftime(Sys.time(), spi_start, units = "secs")
  cat(sprintf("\nMSPI processing completed in %s\n", format_time(as.numeric(spi_elapsed))))
  
  results$mspi <- list(
    raster = mspi_rast,
    timeseries = mspi_df,
    variance = var_rast_spi,  
    stats = mspi_result$stats
  )
  
  # ===== PROCESS MSPEI =====
  cat(sprintf("\n%s\n", paste(rep("-", 80), collapse = "")))
  cat(sprintf("COMPUTING MSPEI (Standardized Precipitation-Evapotranspiration Index)\n"))
  cat(sprintf("%s\n", paste(rep("-", 80), collapse = "")))
  
  spei_start <- Sys.time()
  
  # Reconstruct SPEI time series
  spei_data <- reconstruct_full_timeseries(
    input_dir = file.path(base_dir, "spei_results_seasonal"),
    index_type = "SPEI",
    start_year = start_year,
    end_year = end_year
  )
  
  # Use same mask (already created)
  validation_spei <- validate_basin_data(spei_data$ts_array, mask_info$basin_pixels)
  
  # Compute MSPEI
  mspei_result <- compute_mspi_bazrafshan(
    index_tensor = spei_data$ts_array,
    dates = spei_data$dates,
    basin_pixels = validation_spei$basin_pixels_valid,
    scales = spei_data$scales
  )
  
  # Diagnostic: MSPEI missingness
  diagnose_mspi_missingness(mspei_result$mspi, 
                            validation_spei$basin_pixels_valid, 
                            spei_data$dates, 
                            "MSPEI")
  # Assemble raster
  cat(sprintf("\n[%s] Assembling MSPEI raster stack...\n", format(Sys.time(), "%H:%M:%S")))
  mspei_rast <- rast(spei_data$ref_rast[[1]])
  mspei_rast <- rep(mspei_rast, length(spei_data$dates))
  values(mspei_rast) <- mspei_result$mspi
  terra::time(mspei_rast) <- spei_data$dates
  
  # Save MSPEI outputs
  output_dir_spei <- file.path(base_dir, "mspei_results")
  if (!dir.exists(output_dir_spei)) dir.create(output_dir_spei, recursive = TRUE)
  
  cat(sprintf("\n[%s] Saving MSPEI outputs...\n", format(Sys.time(), "%H:%M:%S")))
  
  # NetCDF output
  nc_file_spei <- file.path(output_dir_spei, sprintf("mspei_monthly_%d-%d.nc", start_year, end_year))
  writeCDF(mspei_rast, nc_file_spei,
           varname = "mspei",
           longname = "Multivariate Standardized Precipitation-Evapotranspiration Index",
           unit = "standardized_index",
           missval = -9999,
           overwrite = TRUE)
  cat(sprintf("OK: MSPEI NetCDF saved: %s\n", basename(nc_file_spei)))
  
  # Excel output
  xlsx_file_spei <- file.path(output_dir_spei, sprintf("mspei_basin_timeseries_%d-%d.xlsx", 
                                                       start_year, end_year))
  mspei_df <- export_basin_timeseries_excel(
    mspi_rast = mspei_rast,
    basin_mask = mask_info$mask,
    dates = spei_data$dates,
    output_file = xlsx_file_spei,
    index_name = "MSPEI"
  )
  
  # Variance explained maps
  var_rast_spei <- rast(spei_data$ref_rast[[1]])
  band_vals <- mspei_result$variance_explained[, 1] * 100
  band_vals[is.na(band_vals)] <- -9999
  values(var_rast_spei) <- band_vals
  names(var_rast_spei) <- "PC1_variance_explained"
  var_file_spei <- file.path(output_dir_spei, "mspei_variance_explained_pc1.nc")
  writeCDF(var_rast_spei, var_file_spei,
           varname  = "PC1_variance_explained",
           longname = "MSPEI PC1 Variance Explained (%)",
           unit     = "percent",
           missval  = -9999,
           overwrite = TRUE)
  cat(sprintf("OK: MSPEI variance map saved: %s\n", basename(var_file_spei)))
  
  spei_elapsed <- difftime(Sys.time(), spei_start, units = "secs")
  cat(sprintf("\nMSPEI processing completed in %s\n", format_time(as.numeric(spei_elapsed))))
  
  results$mspei <- list(
    raster = mspei_rast,
    timeseries = mspei_df,
    variance = var_rast_spei, 
    stats = mspei_result$stats
  )
  
  # ===== COMBINED SUMMARY STATISTICS & DROUGHT CATALOGUE =====
  cat(sprintf("\n%s\n", paste(rep("-", 80), collapse = "")))
  cat(sprintf("EXPORTING COMBINED SUMMARY STATISTICS\n"))
  cat(sprintf("%s\n", paste(rep("-", 80), collapse = "")))
  
  summary_file <- file.path(base_dir,
                            sprintf("MSPI_MSPEI_Summary_Statistics_%d-%d.xlsx",
                                    start_year, end_year))
  summary_results <- export_summary_excel(
    mspi_df    = mspi_df,
    mspei_df   = mspei_df,
    output_file = summary_file
  )
  
  results$summary <- summary_results
  
  # ===== FINAL SUMMARY =====
  total_elapsed <- difftime(Sys.time(), overall_start, units = "secs")
  
  cat(sprintf("\n%s\n", paste(rep("=", 80), collapse = "")))
  cat(sprintf("PROCESSING COMPLETE - TOTAL TIME: %s\n", format_time(as.numeric(total_elapsed))))
  cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
  cat(sprintf("OUTPUT SUMMARY:\n"))
  cat(sprintf("\nMSPI:\n"))
  cat(sprintf("   NetCDF:    %s\n", basename(nc_file_spi)))
  cat(sprintf("   Excel:     %s\n", basename(xlsx_file_spi)))
  cat(sprintf("   Variance:  mspi_variance_explained_pc1.nc\n"))
  cat(sprintf("   Time:      %s\n", format_time(as.numeric(spi_elapsed))))
  cat(sprintf("\nMSPEI:\n"))
  cat(sprintf("   NetCDF:    %s\n", basename(nc_file_spei)))
  cat(sprintf("   Excel:     %s\n", basename(xlsx_file_spei)))
  cat(sprintf("   Variance:  mspei_variance_explained_pc1.nc\n"))
  cat(sprintf("   Time:      %s\n", format_time(as.numeric(spei_elapsed))))
  cat(sprintf("\nCombined Summary:\n"))
  cat(sprintf("   Excel:     %s\n", basename(summary_file)))
  cat(sprintf("\nBasin: %s (%s pixels)\n", basename(basin_path), 
              format(mask_info$n_basin, big.mark = ",")))
  cat(sprintf("%s\n\n", paste(rep("=", 80), collapse = "")))
  
  return(invisible(results))
}

#==============================================================================
# EXECUTION
#==============================================================================
if (interactive()) {
  basin_path <- "Spatial/nechakoBound_dissolve.shp"
  
  if (!file.exists(basin_path)) {
    stop(sprintf("ERROR: Basin file not found: %s", basin_path))
  }
  
  # Compute both MSPI and MSPEI
  results <- compute_both_indices(
    start_year = 1950,
    end_year = 2025,
    base_dir = getwd(),
    basin_path = basin_path
  )
  
  cat("\nACCESSING RESULTS:\n")
  cat("  results$mspi$raster          - MSPI raster stack\n")
  cat("  results$mspi$timeseries      - MSPI basin time series (data.frame)\n")
  cat("  results$mspei$raster         - MSPEI raster stack\n")
  cat("  results$mspei$timeseries     - MSPEI basin time series (data.frame)\n")
  cat("  results$summary$mspi_droughts_full   - Full-period MSPI drought catalogue\n")
  cat("  results$summary$mspei_droughts_full  - Full-period MSPEI drought catalogue\n")
}