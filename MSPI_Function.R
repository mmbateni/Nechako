##############################################
# MSPI & MSPEI Calculation - BAZRAFSHAN 2014 METHOD 
# Computes both indices with .nc and .xlsx outputs
# PCA on SPI/SPEI scales AT SAME CALENDAR TIME
##############################################
library(terra)
library(ncdf4)
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
  
  # AFTER loading precip/pet rasters, BEFORE any calculations:
  target_crs <- "EPSG:3005"  # BC Albers Equal Area
  
  # Reproject basin boundary to match
  if (!same.crs(basin, target_crs)) {
    basin <- project(basin, target_crs)
  }
  
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
# Reconstruct Full Time Series
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
  
  # Load reference raster
  ref_rast <- rast(file_info$filepath[1])
  n_pixels <- ncell(ref_rast)
  cat(sprintf("  -> Reference raster: %d layers, %dx%d cells (%s pixels)\n",
              nlyr(ref_rast), ncol(ref_rast), nrow(ref_rast), 
              format(n_pixels, big.mark = ",")))
  
  # Time dimensions
  n_years <- end_year - start_year + 1
  total_months <- n_years * 12
  
  # Initialize array
  ts_array <- array(NA, dim = c(n_pixels, total_months, length(scales)))
  dimnames(ts_array) <- list(NULL, NULL, sprintf("scale_%02d", scales))
  
  cat(sprintf("  -> Target array dimensions: %s pixels x %d months x %d scales\n",
              format(n_pixels, big.mark = ","), total_months, length(scales)))
  
  # Reconstruct time series
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
      
      for (layer_idx in 1:n_layers_actual) {
        year_offset <- layer_idx - 1
        time_idx <- year_offset * 12 + mon
        
        if (time_idx <= total_months) {
          vals <- values(r[[layer_idx]], mat = FALSE)
          if (length(vals) == n_pixels) {
            ts_array[, time_idx, scale_idx] <- vals
          }
        }
      }
    }
  }
  close(pb)
  
  dates <- seq(as.Date(sprintf("%d-01-01", start_year)), by = "month", length.out = total_months)
  cat(sprintf("\nOK: Time series reconstructed: %s to %s (%d months)\n",
              format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m"), total_months))
  
  # Validate completeness
  total_valid <- sum(!is.na(ts_array))
  total_cells <- prod(dim(ts_array))
  valid_pct_total <- 100 * total_valid / total_cells
  
  cat(sprintf("  -> Data completeness: %.2f%% (%s / %s values)\n",
              valid_pct_total,
              format(total_valid, big.mark = ","),
              format(total_cells, big.mark = ",")))
  
  if (valid_pct_total < 1.0) {
    stop(sprintf("\nERROR: Time series array is %.2f%% NA values!\n", 
                 100 - valid_pct_total))
  }
  
  return(list(
    ts_array = ts_array,
    dates = dates,
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
# 
# HOW THIS WORKS:
# 1. For each calendar month (Jan-Dec), extracts all SPI/SPEI scales ending at
#    the same calendar time (e.g., Dec 1950: [SPI-1(Dec), SPI-2(Dec), ..., SPI-12(Dec)])
# 2. Builds [years × scales] matrices for each month-of-year (12 separate PCAs)
# 3. Performs PCA on centered (NOT scaled) data to extract PC1
# 4. Standardizes PC1 by month: MSPI = PC1 / SD_monthly (Equation 3, Bazrafshan 2014)
# 5. Result integrates multi-scale drought information into single time series
#==============================================================================
compute_mspi_bazrafshan <- function(index_tensor, dates, basin_pixels, scales) {
  n_pixels_total <- dim(index_tensor)[1]
  n_time <- dim(index_tensor)[2]
  n_scales <- length(scales)
  
  mspi_matrix <- matrix(NA, nrow = n_pixels_total, ncol = n_time)
  var_explained_matrix <- matrix(NA, nrow = n_pixels_total, ncol = 12)
  
  n_basin <- length(basin_pixels)
  month_of_year <- as.numeric(format(dates, "%m"))
  
  cat(sprintf("\n[%s] STEP 5: Computing MSPI/MSPEI for %s pixels...\n", 
              format(Sys.time(), "%H:%M:%S"), format(n_basin, big.mark = ",")))
  cat("  -> Method: PCA on scales AT SAME CALENDAR TIME\n")
  
  pb <- txtProgressBar(min = 0, max = n_basin, style = 3)
  stats <- list(successful = 0, failed = 0)
  
  for (i in 1:n_basin) {
    setTxtProgressBar(pb, i)
    if (Sys.info()["sysname"] == "Windows") flush.console()
    
    px <- basin_pixels[i]
    
    # Extract full time series: [time x scales]
    pixel_data <- matrix(index_tensor[px, , ], nrow = n_time, ncol = n_scales)
    
    # Process each month-of-year separately (seasonal PCA)
    for (mon in 1:12) {
      # Get time indices for this month-of-year
      time_idx <- which(month_of_year == mon)
      
      # Extract data matrix: rows=years, cols=scales (ALL ENDING at same calendar month)
      month_data <- pixel_data[time_idx, , drop = FALSE]
      
      # Remove rows with insufficient data
      valid_rows <- rowSums(!is.na(month_data)) >= (n_scales * 0.5)
      if (sum(valid_rows) < 3) next
      
      clean_data <- month_data[valid_rows, , drop = FALSE]
      
      # Center but DON'T scale (per Bazrafshan methodology)
      clean_centered <- scale(clean_data, center = TRUE, scale = FALSE)
      
      # Remove remaining NA values
      clean_complete <- na.omit(clean_centered)
      if (nrow(clean_complete) < 3 || ncol(clean_complete) < 2) next
      
      # Run PCA
      pca_result <- tryCatch({
        prcomp(clean_complete, center = FALSE, scale. = FALSE)
      }, error = function(e) NULL)
      
      if (is.null(pca_result) || length(pca_result$sdev) == 0) next
      
      # Extract PC1 scores
      pc1_scores <- pca_result$x[, 1]
      
      # Assign to MSPI matrix
      valid_time_idx <- time_idx[valid_rows][!is.na(rowSums(clean_centered))]
      if (length(pc1_scores) == length(valid_time_idx)) {
        mspi_matrix[px, valid_time_idx] <- pc1_scores
        
        # Store variance explained
        total_var <- sum(pca_result$sdev^2)
        var_explained_matrix[px, mon] <- pca_result$sdev[1]^2 / total_var
      }
    }
    
    # Standardize final MSPI series by month (Eq. 3 in Bazrafshan)
    for (mon in 1:12) {
      time_idx <- which(month_of_year == mon)
      mspi_vals <- mspi_matrix[px, time_idx]
      valid_idx <- !is.na(mspi_vals)
      
      if (sum(valid_idx) > 1) {
        mspi_std <- (mspi_vals[valid_idx] - mean(mspi_vals[valid_idx])) / sd(mspi_vals[valid_idx])
        mspi_matrix[px, time_idx[valid_idx]] <- mspi_std
      }
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
# Export Basin-Averaged Time Series to Excel
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
  
  # Compute basin-wide statistics
  basin_mean <- colMeans(basin_values, na.rm = TRUE)
  basin_median <- apply(basin_values, 2, median, na.rm = TRUE)
  basin_min <- apply(basin_values, 2, min, na.rm = TRUE)
  basin_max <- apply(basin_values, 2, max, na.rm = TRUE)
  basin_sd <- apply(basin_values, 2, sd, na.rm = TRUE)
  
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
  var_monthly_rast_spi <- rast(spi_data$ref_rast[[1]], nlyrs = 12)
  for (mon in 1:12) {
    band_vals <- mspi_result$variance_explained[, mon] * 100
    band_vals[is.na(band_vals)] <- -9999
    values(var_monthly_rast_spi[[mon]]) <- band_vals
  }
  names(var_monthly_rast_spi) <- sprintf("PC1_variance_month_%02d", 1:12)
  
  var_file_spi <- file.path(output_dir_spi, "mspi_variance_explained_monthly_12band.tif")
  writeRaster(var_monthly_rast_spi, var_file_spi, datatype = "FLT4S", 
              overwrite = TRUE, NAflag = -9999)
  cat(sprintf("OK: MSPI variance map saved: %s\n", basename(var_file_spi)))
  
  spi_elapsed <- difftime(Sys.time(), spi_start, units = "secs")
  cat(sprintf("\nMSPI processing completed in %s\n", format_time(as.numeric(spi_elapsed))))
  
  results$mspi <- list(
    raster = mspi_rast,
    timeseries = mspi_df,
    variance = var_monthly_rast_spi,
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
  var_monthly_rast_spei <- rast(spei_data$ref_rast[[1]], nlyrs = 12)
  for (mon in 1:12) {
    band_vals <- mspei_result$variance_explained[, mon] * 100
    band_vals[is.na(band_vals)] <- -9999
    values(var_monthly_rast_spei[[mon]]) <- band_vals
  }
  names(var_monthly_rast_spei) <- sprintf("PC1_variance_month_%02d", 1:12)
  
  var_file_spei <- file.path(output_dir_spei, "mspei_variance_explained_monthly_12band.tif")
  writeRaster(var_monthly_rast_spei, var_file_spei, datatype = "FLT4S", 
              overwrite = TRUE, NAflag = -9999)
  cat(sprintf("OK: MSPEI variance map saved: %s\n", basename(var_file_spei)))
  
  spei_elapsed <- difftime(Sys.time(), spei_start, units = "secs")
  cat(sprintf("\nMSPEI processing completed in %s\n", format_time(as.numeric(spei_elapsed))))
  
  results$mspei <- list(
    raster = mspei_rast,
    timeseries = mspei_df,
    variance = var_monthly_rast_spei,
    stats = mspei_result$stats
  )
  
  # ===== FINAL SUMMARY =====
  total_elapsed <- difftime(Sys.time(), overall_start, units = "secs")
  
  cat(sprintf("\n%s\n", paste(rep("=", 80), collapse = "")))
  cat(sprintf("PROCESSING COMPLETE - TOTAL TIME: %s\n", format_time(as.numeric(total_elapsed))))
  cat(sprintf("%s\n", paste(rep("=", 80), collapse = "")))
  cat(sprintf("OUTPUT SUMMARY:\n"))
  cat(sprintf("\nMSPI:\n"))
  cat(sprintf("   NetCDF: %s\n", basename(nc_file_spi)))
  cat(sprintf("   Excel:  %s\n", basename(xlsx_file_spi)))
  cat(sprintf("   Time:   %s\n", format_time(as.numeric(spi_elapsed))))
  cat(sprintf("\nMSPEI:\n"))
  cat(sprintf("   NetCDF: %s\n", basename(nc_file_spei)))
  cat(sprintf("   Excel:  %s\n", basename(xlsx_file_spei)))
  cat(sprintf("   Time:   %s\n", format_time(as.numeric(spei_elapsed))))
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
  cat("  results$mspi$raster       - MSPI raster stack\n")
  cat("  results$mspi$timeseries   - MSPI basin time series (data.frame)\n")
  cat("  results$mspei$raster      - MSPEI raster stack\n")
  cat("  results$mspei$timeseries  - MSPEI basin time series (data.frame)\n")
}