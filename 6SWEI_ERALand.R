##############################################
# CORRECTED SWEI CALCULATION WITH SEASONAL APPROACH
# 1-month output resolution with 3-month SWE smoothing + 3-month SCF mask
# Following Huning & AghaKouchak (2020) and Hajivand Paydari et al. (2025)
##############################################

library(terra)
library(parallel)
library(ncdf4)
library(zoo)
library(writexl)

# Configuration
setwd("D:/Nechako_Drought/Nechako")
out_dir <- "swei_seasonal"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load basin boundary
basin_path <- "Spatial/nechakoBound_dissolve.shp"
basin <- if (file.exists(basin_path)) vect(basin_path) else NULL
if (is.null(basin)) stop("Basin boundary not found: ", basin_path)

# Load SWE and SCF data
swe_file <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
scf_file <- "monthly_data_direct/snow_cover_monthly.nc"

if (!file.exists(swe_file)) stop("SWE file not found: ", swe_file)
if (!file.exists(scf_file)) stop("SCF file not found: ", scf_file)

cat("Reading SWE and SCF data...\n")
swe <- rast(swe_file)
scf_stack <- rast(scf_file)

# Time extraction function
extract_time_dimension <- function(raster_obj, file_path) {
  t <- time(raster_obj)
  if (!all(is.na(t))) return(as.Date(t))
  
  nc <- nc_open(file_path)
  on.exit(nc_close(nc))
  
  if ("valid_time" %in% names(nc$var)) {
    tv <- ncvar_get(nc, "valid_time")
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  if ("time" %in% names(nc$dim)) {
    tv <- ncvar_get(nc, "time")
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  
  n <- nlyr(raster_obj)
  start <- as.Date("1950-01-01")
  seq(start, by = "month", length.out = n)
}

dates_swe <- extract_time_dimension(swe, swe_file)
dates_scf <- extract_time_dimension(scf_stack, scf_file)

terra::time(swe) <- dates_swe
terra::time(scf_stack) <- dates_scf

# Convert SCF to fraction if in percent
scf_max <- global(scf_stack[[1]], "max", na.rm = TRUE)$max
if (!is.na(scf_max) && scf_max > 1.5) scf_stack <- scf_stack / 100

# Reproject to common grid if needed
if (!compareGeom(scf_stack[[1]], swe[[1]], stopOnError = FALSE)) {
  swe <- resample(swe, scf_stack, method = "bilinear")
}

# Reproject basin to match CRS
if (!is.null(basin) && !same.crs(basin, crs(scf_stack))) {
  basin <- project(basin, crs(scf_stack))
}

# Mask SWE to basin
swe <- mask(swe, basin, inverse = FALSE, touches = TRUE)

# Convert SWE to mm if needed (meters -> mm)
swe_mean_sample <- global(swe[[1]], "mean", na.rm = TRUE)$mean
if (!is.na(swe_mean_sample) && swe_mean_sample < 10) {
  swe <- swe * 1000
  cat("✓ Converted SWE from meters to mm\n")
}

# Prepare matrices
swe_matrix <- values(swe, mat = TRUE)
scf_matrix <- values(scf_stack, mat = TRUE)
dates <- dates_swe
months_all <- as.integer(format(dates, "%m"))
years_all <- as.integer(format(dates, "%Y"))

cat(sprintf("SWE matrix dimensions: %d pixels x %d time steps\n", nrow(swe_matrix), ncol(swe_matrix)))
cat(sprintf("Time period: %s to %s\n", min(dates), max(dates)))

# ============================================================================
# STEP 1: Apply 3-month moving average to RAW SWE values (centered window)
# This is PREPROCESSING - not the drought timescale itself
# ============================================================================
cat("\n===== STEP 1: Applying 3-month centered moving average to RAW SWE =====\n")
swe_smoothed <- matrix(NA, nrow = nrow(swe_matrix), ncol = ncol(swe_matrix))

for (i in 1:nrow(swe_matrix)) {
  # Centered 3-month window (Huning & AghaKouchak 2020 approach)
  swe_smoothed[i, ] <- zoo::rollapply(swe_matrix[i, ], 
                                      width = 3, 
                                      FUN = mean, 
                                      align = "center", 
                                      fill = NA, 
                                      na.rm = TRUE)
}
cat("✓ 3-month moving average applied to SWE values\n")

# ============================================================================
# STEP 2: LOAD SCF MASK FOR k=3 (CRITICAL CORRECTION)
# Since we're smoothing SWE with 3 months, mask must reflect 3-month snow regimes
# ============================================================================
cat("\n===== STEP 2: Loading SCF mask for k=3 (3-month climatology) =====\n")
scf_diag_dir <- "scf_timescale_diagnostics"

# CORRECTED: Use k=3 for mask/threshold selection (not k=1)
target_k <- 3  # Must match smoothing window size

# Load diagnostics to get optimal threshold for k=3
scf_diag_csv <- file.path(scf_diag_dir, "scf_timescale_final_recommendations.csv")
scf_threshold <- 0.10  # Default fallback

if (file.exists(scf_diag_csv)) {
  diag_df <- read.csv(scf_diag_csv, stringsAsFactors = FALSE)
  if ("timescale" %in% names(diag_df) && "chosen_threshold" %in% names(diag_df)) {
    k3_row <- diag_df[as.character(diag_df$timescale) == as.character(target_k), ]
    if (nrow(k3_row) > 0) {
      scf_threshold <- as.numeric(k3_row$chosen_threshold[1])
      cat(sprintf("✓ Loaded optimal SCF threshold for k=%d: %.4f\n", target_k, scf_threshold))
    } else {
      cat(sprintf("Warning: No k=%d entry in diagnostics. Using default threshold %.4f\n", 
                  target_k, scf_threshold))
    }
  }
}

# Load pre-computed k=3 mask (MUST match target_k)
mask_file_pattern <- sprintf("scf_mask_k%d_*.nc", target_k)
mask_candidates <- list.files(scf_diag_dir, 
                              pattern = glob2rx(mask_file_pattern), 
                              full.names = TRUE)
if (length(mask_candidates) == 0) {
  stop(sprintf("CRITICAL ERROR: No SCF mask found for k=%d in %s", target_k, scf_diag_dir))
}

# Select best match (closest to our threshold)
best_mask <- mask_candidates[[1]]
if (length(mask_candidates) > 1) {
  # Find mask with threshold closest to our selected value
  min_diff <- Inf
  for (f in mask_candidates) {
    # Extract threshold from filename (e.g., "scf_mask_k03_threshold_0.1000.nc")
    m <- regexec("threshold_([0-9.]+)\\.nc$", basename(f))
    if (length(m[[1]]) > 1) {
      th_in_file <- as.numeric(regmatches(basename(f), m)[[1]][2])
      diff <- abs(th_in_file - scf_threshold)
      if (diff < min_diff) {
        min_diff <- diff
        best_mask <- f
      }
    }
  }
}

cat(sprintf("✓ Loading SCF mask: %s\n", basename(best_mask)))
scf_mask_stack <- rast(best_mask)
scf_mask_list <- vector("list", 12)

# Extract monthly masks (12 layers = 12 calendar months)
for (m in 1:12) {
  if (m <= nlyr(scf_mask_stack)) {
    mask_vals <- values(scf_mask_stack[[m]])
    scf_mask_list[[m]] <- !is.na(mask_vals) & (mask_vals == 1)
  } else {
    scf_mask_list[[m]] <- rep(FALSE, ncell(swe))
  }
}
cat(sprintf("✓ Loaded k=%d SCF masks for 12 calendar months\n", target_k))

# ============================================================================
# STEP 3: Seasonal SWEI calculation function (Kao & Govindaraju approach)
# Applied to 3-month smoothed SWE values
# ============================================================================
gringorten_swei_seasonal <- function(v_smoothed, dates_vec, scf_mask_list, cell_index, eps = 1e-6) {
  # Defensive checks
  if (length(v_smoothed) == 0 || all(is.na(v_smoothed)) || is.null(v_smoothed)) {
    return(rep(NA_real_, length(dates_vec)))
  }
  
  v_clean <- v_smoothed
  v_clean[!is.finite(v_clean)] <- NA_real_
  
  n <- length(v_clean)
  z <- rep(NA_real_, n)
  mon <- as.integer(format(dates_vec, "%m"))
  
  # Process each calendar month separately (seasonal grouping)
  for (m in 1:12) {
    # Get indices for this calendar month across all years
    idx <- which(mon == m & is.finite(v_clean))
    
    # Skip if insufficient data
    if (length(idx) < 20) next
    
    # Apply SCF mask for this calendar month (k=3 mask)
    if (is.null(scf_mask_list) || length(scf_mask_list) < m) next
    mask_vec <- scf_mask_list[[m]]
    if (length(mask_vec) < cell_index || !mask_vec[cell_index]) next
    
    # Extract values for this calendar month
    samp <- v_clean[idx]
    
    # Replace zeros with tiny positive value (0.001% of smallest non-zero)
    # Following Hajivand Paydari et al. 2025 methodology
    zero_idx <- which(samp == 0 & !is.na(samp))
    if (length(zero_idx) > 0) {
      nonzero_vals <- samp[samp > 0 & !is.na(samp)]
      if (length(nonzero_vals) > 0) {
        tiny <- min(nonzero_vals) * 1e-5  # 0.001% = 1e-5
        samp[zero_idx] <- tiny
      }
    }
    
    # Remove remaining NAs
    valid_idx <- which(!is.na(samp))
    if (length(valid_idx) < 3) next
    
    samp_valid <- samp[valid_idx]
    idx_valid <- idx[valid_idx]
    
    # Compute Gringorten probabilities (non-parametric)
    r <- rank(samp_valid, ties.method = "average")
    N <- length(samp_valid)
    p_val <- (r - 0.44) / (N + 0.12)  # Gringorten plotting position
    
    # Clip probabilities to avoid Inf in qnorm
    p_val <- pmax(pmin(p_val, 1 - eps), eps)
    
    # Transform to standard normal (SWEI)
    z_val <- qnorm(p_val)
    
    # Assign back to correct positions
    z[idx_valid] <- z_val
  }
  
  # Clip extreme values (±4.75 corresponds to ~1:10,000 return period)
  z[!is.na(z) & z < -4.75] <- -4.75
  z[!is.na(z) & z >  4.75] <-  4.75
  
  return(z)
}

# ============================================================================
# STEP 4: Parallel computation of seasonal SWEI on smoothed data
# Output resolution: 1-month (k=1) but based on 3-month smoothed SWE
# ============================================================================
cat("\n===== STEP 3: Computing Seasonal SWEI (1-month output from 3-month smoothed SWE) =====\n")
n_cores <- max(1, detectCores() - 1)
cat(sprintf("Using %d CPU cores for parallel computation\n", n_cores))
eps <- 1e-6
cl <- makeCluster(n_cores)
clusterExport(cl, c("swe_smoothed", "dates", "gringorten_swei_seasonal", 
                    "scf_mask_list", "eps"), envir = environment())
clusterEvalQ(cl, { library(zoo) })

eps <- 1e-6
swei_results <- parLapply(cl, 1:nrow(swe_smoothed), function(i) {
  tryCatch({
    gringorten_swei_seasonal(swe_smoothed[i, ], dates, scf_mask_list, i, eps = eps)
  }, error = function(e) {
    rep(NA_real_, length(dates))
  })
})

stopCluster(cl)

# Combine results
swei_indices <- do.call(rbind, swei_results)

# Verify dimensions
if (is.null(dim(swei_indices)) || nrow(swei_indices) != nrow(swe_matrix)) {
  stop("SWEI calculation failed: dimension mismatch")
}

cat(sprintf("✓ Computed seasonal SWEI for %d pixels\n", nrow(swei_indices)))

# ============================================================================
# STEP 5: Basin-averaged time series and outputs
# ============================================================================
cat("\n===== STEP 4: Generating outputs =====\n")

# Compute basin average (masked by original SWE presence)
# Compute basin average (masked by original SWE presence)
basin_mask <- !is.na(swe_matrix[, 1])
if (sum(basin_mask) == 0) stop("No valid basin pixels found!")

swei_basin_avg <- apply(swei_indices[basin_mask, , drop = FALSE], 2, mean, na.rm = TRUE)

# Create time series dataframe
swei_ts <- data.frame(
  date = dates,
  year = years_all,
  month = months_all,
  swei_basin_avg = swei_basin_avg,
  swe_raw_mm = colMeans(swe_matrix[basin_mask, , drop = FALSE], na.rm = TRUE),
  swe_smoothed_mm = colMeans(swe_smoothed[basin_mask, , drop = FALSE], na.rm = TRUE)
)

# Save basin time series
csv_file <- file.path(out_dir, "swei_k01_smooth3_timeseries.csv")
write.csv(swei_ts, csv_file, row.names = FALSE)
cat(sprintf("✓ Saved basin-averaged time series: %s\n", basename(csv_file)))

# Save full time series NetCDF
# Create a new raster stack with clean time dimension
swei_full_rast <- rast(swe[[1]])                     # template from first layer
swei_full_rast <- rep(swei_full_rast, ncol(swei_indices))  # create empty layers
values(swei_full_rast) <- swei_indices
terra::time(swei_full_rast) <- dates                  # assign sorted dates

nc_full_file <- file.path(out_dir, "swei_k01_smooth3_full_timeseries.nc")
writeCDF(swei_full_rast, nc_full_file,
         varname = "swei",
         longname = "Seasonal SWEI (1-month output, 3-month SWE smoothing, k=3 SCF mask)",
         unit = "standardized_index",
         missval = -9999,
         overwrite = TRUE)
cat(sprintf("✓ Saved full time series NetCDF: %s\n", basename(nc_full_file)))

# ============================================================================
# STEP 6: Summary statistics
# ============================================================================
summary_file <- file.path(out_dir, "swei_k01_smooth3_summary.txt")
sink(summary_file)
cat("SEASONAL SWEI SUMMARY (1-MONTH OUTPUT WITH 3-MONTH SWE SMOOTHING)\n")
cat(sprintf("Calculation date: %s\n\n", Sys.time()))
cat(sprintf("Methodology:\n"))
cat("  • 3-month centered moving average applied to RAW SWE values (preprocessing)\n")
cat("  • SCF mask derived from 3-month climatology (k=3) - CRITICAL ALIGNMENT\n")
cat("  • Seasonal grouping by calendar month (Kao & Govindaraju 2010 approach)\n")
cat("  • Zero SWE values replaced with 0.001% of smallest non-zero value\n")
cat("  • Gringorten plotting position for non-parametric probability estimation\n")
cat("  • Transformation to standard normal distribution (mean=0, sd=1)\n")
cat("  • Output resolution: 1-month (k=1) but based on smoothed values\n\n")
cat(sprintf("Grid dimensions: %d x %d pixels\n", ncol(swe), nrow(swe)))
cat(sprintf("Valid basin pixels: %d\n", sum(basin_mask)))
cat(sprintf("Time period: %s to %s (%d months)\n\n", min(dates), max(dates), length(dates)))
cat(sprintf("SCF mask configuration:\n"))
cat(sprintf("  • Timescale: k=%d months\n", target_k))
cat(sprintf("  • Threshold: %.4f\n", scf_threshold))
cat(sprintf("  • Source: %s\n\n", basename(best_mask)))

# Drought frequency statistics
cat("Drought frequency by severity (entire period):\n")
cat(sprintf("  Exceptional (SWEI < -2.0): %.2f%%\n", 
            100 * mean(swei_indices[basin_mask, ] < -2.0, na.rm = TRUE)))
cat(sprintf("  Extreme     (SWEI < -1.6): %.2f%%\n", 
            100 * mean(swei_indices[basin_mask, ] < -1.6, na.rm = TRUE)))
cat(sprintf("  Severe      (SWEI < -1.3): %.2f%%\n", 
            100 * mean(swei_indices[basin_mask, ] < -1.3, na.rm = TRUE)))
cat(sprintf("  Moderate    (SWEI < -0.8): %.2f%%\n", 
            100 * mean(swei_indices[basin_mask, ] < -0.8, na.rm = TRUE)))
cat(sprintf("  Near Normal (-0.5 to +0.5): %.2f%%\n", 
            100 * mean(abs(swei_indices[basin_mask, ]) <= 0.5, na.rm = TRUE)))

# NA analysis
na_rate <- 100 * mean(is.na(swei_indices[basin_mask, ]))
cat(sprintf("\nNA rate in basin pixels: %.3f%%\n", na_rate))
sink()
cat(sprintf("✓ Saved summary statistics: %s\n", basename(summary_file)))

cat("\n============================================================\n")
cat("SEASONAL SWEI CALCULATION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Output directory: %s\n", normalizePath(out_dir)))
cat("\nKey methodological alignment:\n")
cat("  ✓ 1-month output resolution (k=1 timescale)\n")
cat("  ✓ 3-month moving average applied to RAW SWE (preprocessing step)\n")
cat("  ✓ SCF mask derived from 3-month climatology (k=3) - PROPER ALIGNMENT\n")
cat("  ✓ Seasonal grouping by calendar month (accounts for snow seasonality)\n")
cat("  ✓ Non-parametric Gringorten approach (distribution-free)\n")
cat("\nCRITICAL NOTE: The SCF mask timescale (k=3) matches the SWE smoothing\n")
cat("window size, ensuring methodological consistency with Huning & AghaKouchak (2020)\n")