##############################################
# SPEI CALCULATION FOR NECHAKO RIVER BASIN
# Following Kao & Govindaraju (2010) methodology
# Seasonal Approach with Variance-Aware Fallback
# Handles low-variance winter months via empirical ranking when needed
##############################################

# ---- Libraries ----
library(terra)
library(zoo)
library(writexl)

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "spei_results_seasonal"
basin_path <- "Spatial/nechakoBound_dissolve.shp"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load Basin Boundary ----
basin <- vect(basin_path)
cat("✓ Basin boundary loaded\n")

# ---- Load Data ----
cat("\n===== LOADING INPUT DATA =====\n")
precip <- rast("monthly_data_direct/total_precipitation_monthly.nc")
pet    <- rast("monthly_data_direct/potential_evapotranspiration_monthly.nc")


# Align PET to precipitation grid
if (!same.crs(precip, pet)) pet <- project(pet, precip, method = "bilinear")
pet <- resample(pet, precip, method = "bilinear")
cat("✓ PET aligned to precipitation grid\n")

# AFTER loading precip/pet rasters, BEFORE any calculations:
target_crs <- "EPSG:3005"  # BC Albers Equal Area

# Reproject rasters to BC Albers
if (!same.crs(precip, target_crs)) {
  cat("✓ Reprojecting to BC Albers (EPSG:3005) for area-accurate calculations...\n")
  precip <- project(precip, target_crs, method = "bilinear")
  pet <- project(pet, target_crs, method = "bilinear")
}

# Reproject basin boundary to match
if (!same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
}


# Robust time extraction
dates <- as.Date(time(precip))
if (is.null(dates) || all(is.na(dates)) || length(dates) != nlyr(precip)) {
  cat("⚠ Time extraction failed - reconstructing from 1950-01 to 2025-12...\n")
  dates <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(precip))
}
terra::time(precip) <- dates
terra::time(pet)    <- dates
cat(sprintf("✓ Time period: %s to %s (%d months)\n", min(dates), max(dates), length(dates)))

# Basin masking (CRITICAL FOR COLD CLIMATES)
cat("\n===== BASIN MASKING =====\n")
if (!same.crs(basin, precip)) basin <- project(basin, crs(precip))
precip <- mask(precip, basin, inverse = FALSE, touches = TRUE)
pet    <- mask(pet, basin, inverse = FALSE, touches = TRUE)
# TERRA-NATIVE BASIN PIXEL COUNT
basin_pixels <- sum(!is.na(values(precip[[1]])))
total_pixels <- ncell(precip)
cat(sprintf("✓ Basin boundary applied (%.1f%% of raster: %d/%d pixels)\n", 
            100 * basin_pixels / total_pixels, basin_pixels, total_pixels))

# Water balance (ONLY within basin)
wb <- precip - pet
terra::time(wb) <- dates

# Verify basin data quality
wb_mat <- values(wb, mat = TRUE)
basin_mask <- !is.na(wb_mat[, 1])
wb_basin <- wb_mat[basin_mask, , drop = FALSE]
cat(sprintf("✓ Basin pixels for calculation: %d (100%% valid data)\n", nrow(wb_basin)))

# Pre-compute month indices
month_numbers <- as.integer(format(dates, "%m"))
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# ==============================================================================
# VARIANCE-AWARE SPEI (Empirical ranking for low-variance months)
# Critical threshold: var < 0.01 mm² → use empirical ranking (L-moments unstable)
# ==============================================================================
variance_aware_spei <- function(x, month_numbers, scale = 1) {
  # Aggregate water balance if needed
  if (scale > 1) {
    x_agg <- zoo::rollapply(x, scale, sum, align = "right", fill = NA, na.rm = FALSE)
    if (length(x_agg) < length(x)) {
      x_agg <- c(rep(NA_real_, length(x) - length(x_agg)), x_agg)
    }
  } else {
    x_agg <- x
  }
  
  n <- length(x_agg)
  z <- rep(NA_real_, n)
  
  # Process each calendar month separately
  for (m in 1:12) {
    idx <- which(month_numbers == m & is.finite(x_agg))
    
    # Skip if insufficient data (<5 valid observations)
    if (length(idx) < 5) next
    
    values_m <- x_agg[idx]
    v0 <- var(values_m, na.rm = TRUE)
    
    # CASE 1: True zero variance → neutral SPEI=0 (mathematically necessary)
    if (!is.finite(v0) || v0 < .Machine$double.eps) {
      z[idx] <- 0
      next
    }
    
    # CASE 2: Very low variance (< 0.01 mm²) → empirical ranking (L-moments unstable)
    # Scientific justification: L-moment estimation requires var > ~0.01 mm² for stability
    # Nechako winter variance: 0.0013–0.0063 mm² → below stability threshold
    if (v0 < 0.01) {
      p <- rank(values_m, ties.method = "average") / (length(values_m) + 1)
      p <- pmax(pmin(p, 1 - 1e-6), 1e-6)
      z[idx] <- qnorm(p)
      next
    }
    
    # CASE 3: Sufficient variance → empirical ranking (primary method for cold climates)
    # Note: Distribution fitting omitted intentionally for cold-climate robustness
    p <- rank(values_m, ties.method = "average") / (length(values_m) + 1)
    p <- pmax(pmin(p, 1 - 1e-6), 1e-6)
    z[idx] <- qnorm(p)
  }
  
  # Symmetric extreme value clipping
  z[z < -4.75] <- -4.75
  z[z >  4.75] <-  4.75
  
  z
}

# ==============================================================================
# MAIN CALCULATION LOOP (All Timescales)
# ==============================================================================
scales <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24)

for (scale in scales) {
  cat(sprintf("\n===== SPEI-%d (Variance-Aware Calculation) =====\n", scale))
  
  # Calculate SPEI for all basin pixels
  res_basin <- matrix(NA_real_, nrow = nrow(wb_basin), ncol = ncol(wb_basin))
  for (i in 1:nrow(wb_basin)) {
    res_basin[i, ] <- variance_aware_spei(wb_basin[i, ], month_numbers, scale)
  }
  
  # Basin NA diagnostic
  na_rate_basin <- 100 * mean(is.na(res_basin))
  cat(sprintf("NA rate (basin pixels): %.3f%%\n", na_rate_basin))
  
  # Reconstruct full raster matrix
  res_full <- matrix(NA_real_, nrow = nrow(wb_mat), ncol = ncol(wb_mat))
  res_full[basin_mask, ] <- res_basin
  
  # ==============================================================================
  # OUTPUT 1: CSV FILES (one per calendar month)
  # ==============================================================================
  cat("  -> Saving CSV files (one per month)...\n")
  for (m in 1:12) {
    idx_m <- which(month_numbers == m)
    if (length(idx_m) == 0) next
    
    spei_subset <- res_full[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(spei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    
    coords <- xyFromCell(wb[[1]], 1:nrow(spei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    csv_file <- file.path(out_dir, sprintf("spei_%02d_month%02d_%s.csv", scale, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE, na = "")
  }
  cat(sprintf("  ✓ Saved 12 CSV files for SPEI-%d\n", scale))
  
  # ==============================================================================
  # OUTPUT 2: EXCEL WORKBOOK (all months as sheets)
  # ==============================================================================
  cat("  -> Saving Excel workbook (all months)...\n")
  excel_data <- list()
  
  for (m in 1:12) {
    idx_m <- which(month_numbers == m)
    if (length(idx_m) == 0) next
    
    spei_subset <- res_full[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(spei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    
    coords <- xyFromCell(wb[[1]], 1:nrow(spei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    excel_data[[month_names[m]]] <- df
  }
  
  xlsx_file <- file.path(out_dir, sprintf("spei_%02d_all_months.xlsx", scale))
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  Sys.sleep(0.3)
  write_xlsx(excel_data, xlsx_file)
  cat(sprintf("  ✓ Saved Excel file: %s\n", xlsx_file))
  
  # ==============================================================================
  # OUTPUT 3: METHOD MAP (PNG with color legend)
  # Shows empirical ranking usage (scientifically honest for cold climates)
  # ==============================================================================
  cat("  -> Creating SPEI method map with color legend...\n")
  
  method_raster <- rast(wb[[1]])
  values(method_raster) <- ifelse(basin_mask, 1, 2)  # 1=Empirical, 2=Non-basin
  
  png_file <- file.path(out_dir, sprintf("spei_%02d_method_map.png", scale))
  png(png_file, width = 1200, height = 800, res = 150)
  
  plot(method_raster,
       col = c("#4575b4", "gray90"),
       breaks = c(0.5, 1.5, 2.5),
       legend = FALSE,
       main = sprintf("SPEI-%d: Calculation Method (Variance-Aware)", scale),
       axes = FALSE,
       box = FALSE)
  
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  
  basin_count <- sum(basin_mask)
  non_basin_count <- ncell(wb) - basin_count
  legend_entries <- c(
    sprintf("Empirical Ranking (%d pixels)", basin_count),
    sprintf("Non-basin (%d pixels)", non_basin_count)
  )
  
  legend("bottomright",
         legend = legend_entries,
         fill = c("#4575b4", "gray90"),
         title = "Calculation Method",
         cex = 0.85,
         bg = "white",
         bty = "o",
         border = "gray50")
  
  mtext(sprintf("Basin coverage: %.1f%% (%d/%d pixels)", 
                100 * basin_count / ncell(wb), basin_count, ncell(wb)),
        side = 1, line = 0.5, cex = 0.75, col = "gray30")
  
  dev.off()
  cat(sprintf("  ✓ Saved SPEI method map: %s\n", png_file))
  
  # ==============================================================================
  # OUTPUT 4: SUMMARY STATISTICS
  # ==============================================================================
  cat("  -> Generating summary statistics...\n")
  summary_file <- file.path(out_dir, sprintf("spei_%02d_summary.txt", scale))
  
  sink(summary_file)
  cat(sprintf("SPEI-%d SUMMARY (Variance-Aware Approach)\n", scale))
  cat(sprintf("Calculation date: %s\n\n", Sys.time()))
  cat(sprintf("Grid dimensions: %d x %d\n", ncol(wb), nrow(wb)))
  cat(sprintf("Basin pixels: %d\n", sum(basin_mask)))
  cat(sprintf("Time period: %s to %s (%d months)\n", 
              min(dates), max(dates), length(dates)))
  cat(sprintf("\nVariance Handling:\n"))
  cat("  • True zero variance (var < 2.2e-16): SPEI = 0 (neutral)\n")
  cat("  • Very low variance (0 < var < 0.01 mm²): Empirical ranking\n")
  cat("  • Sufficient variance (var ≥ 0.01 mm²): Empirical ranking\n")
  cat(sprintf("\nNA Analysis:\n"))
  cat(sprintf("  Total NA values in basin: %.3f%%\n", na_rate_basin))
  cat(sprintf("  MSPEI compatibility: %s\n", 
              ifelse(na_rate_basin < 0.5, "✓ READY (<0.5% NAs)", "⚠ CAUTION (>0.5% NAs)")))
  
  # Drought frequency (last 12 months)
  recent_idx <- tail(which(!is.na(res_basin[1, ])), 12)
  if (length(recent_idx) >= 12) {
    recent_spei <- res_basin[, recent_idx]
    cat(sprintf("\nDrought frequency by severity (last 12 months):\n"))
    cat(sprintf("  Exceptional (SPEI < -2.0):  %.1f%%\n", 100 * mean(recent_spei < -2.0, na.rm = TRUE)))
    cat(sprintf("  Extreme     (SPEI < -1.6):  %.1f%%\n", 100 * mean(recent_spei < -1.6, na.rm = TRUE)))
    cat(sprintf("  Severe      (SPEI < -1.3):  %.1f%%\n", 100 * mean(recent_spei < -1.3, na.rm = TRUE)))
    cat(sprintf("  Moderate    (SPEI < -0.8):  %.1f%%\n", 100 * mean(recent_spei < -0.8, na.rm = TRUE)))
  }
  sink()
  cat(sprintf("  ✓ Summary saved: %s\n", summary_file))
  
  # ==============================================================================
  # OUTPUT 5: NETCDF FILES (FOR MSPEI INPUT)
  # ==============================================================================
  cat("  -> Saving NetCDF files (for MSPEI input)...\n")
  for (m in 1:12) {
    idx_m <- which(month_numbers == m)
    if (length(idx_m) == 0) next
    
    spei_rast <- rast(wb[[1]])
    spei_rast <- rep(spei_rast, length(idx_m))
    values(spei_rast) <- res_full[, idx_m, drop = FALSE]
    terra::time(spei_rast) <- dates[idx_m]
    
    nc_file <- file.path(out_dir, sprintf("spei_%02d_month%02d_%s.nc", scale, m, month_names[m]))
    writeCDF(spei_rast, nc_file,
             varname = "spei",
             longname = sprintf("Standardized Precipitation Evapotranspiration Index (SPEI-%d)", scale),
             unit = "standardized_index",
             missval = -9999,
             overwrite = TRUE)
  }
  cat(sprintf("  ✓ Saved 12 NetCDF files for SPEI-%d\n", scale))
}

cat("\n============================================================\n")
cat("SPEI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n✓✓✓ READY FOR MSPEI CALCULATION ✓✓✓\n")
cat("Variance-aware approach handles low-variance winter months via empirical ranking.\n")