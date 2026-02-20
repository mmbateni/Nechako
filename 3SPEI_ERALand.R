# ##############################################
# SPEI CALCULATION FOR NECHAKO RIVER BASIN
# Following Kao & Govindaraju (2010) methodology
# Seasonal Approach with Variance-Aware Fallback
# Case 1: True zero variance         → SPEI = 0
# Case 2: Low variance (< 0.001 mm²) → Gringorten empirical ranking
# Case 3: Sufficient variance        → Parametric (GLO/PE3/GEV via L-moments),
#                                     Gringorten fallback if KS test fails
# ##############################################

# ---- Libraries ----
library(terra)
library(zoo)
library(writexl)
library(lmomco)   # Parametric distribution fitting via L-moments
library(parallel) # Parallel processing

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "spei_results_seasonal"
basin_path <- "Spatial/nechakoBound_dissolve.shp"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Suppress terra progress bars for clean console output
terraOptions(progress = 0)

# ---- Load Basin Boundary ----
basin <- vect(basin_path)
cat("✓ Basin boundary loaded\n")

# ---- Load Data ----
cat("\n===== LOADING INPUT DATA =====\n")
precip <- rast("monthly_data_direct/total_precipitation_monthly.nc")
pet    <- rast("monthly_data_direct/potential_evapotranspiration_monthly.nc")
precip <- precip * 1000  # m → mm

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
basin_pixels <- global(precip[[1]], "notNA")$notNA
total_pixels <- ncell(precip)                 # Should be 787056
cat(sprintf("✓ Basin boundary applied (%.1f%% of raster: %d/%d pixels)\n",
            100 * basin_pixels / total_pixels, basin_pixels, total_pixels))

# Water balance with days-in-month correction
cat("\n===== WATER BALANCE CALCULATION =====\n")
wb <- precip - pet

# Get number of days in each month
month_nums <- as.integer(format(dates, "%m"))
year_nums  <- as.integer(format(dates, "%Y"))
days_lookup <- c(31,28,31,30,31,30,31,31,30,31,30,31)
days_in_month <- days_lookup[month_nums]

# Adjust for leap years
leap_years <- (year_nums %% 4 == 0 & year_nums %% 100 != 0) | (year_nums %% 400 == 0)
days_in_month[month_nums == 2 & leap_years] <- 29

cat("✓ Converting from mm/day to mm/month (multiplying by days in month)...\n")
# FIX: Corrected indexing wb[[i]]
for (i in 1:nlyr(wb)) {
  wb[[i]] <- wb[[i]] * days_in_month[i]
}
terra::time(wb) <- dates
cat(sprintf("✓ Water balance now in mm/month (range: %d-%d days per month)\n",
            min(days_in_month), max(days_in_month)))

# Verify basin data quality
wb_mat <- values(wb, mat = TRUE)
basin_mask <- !is.na(wb_mat)
wb_basin <- wb_mat
cat(sprintf("✓ Basin pixels for calculation: %d (100%% valid data)\n", nrow(wb_basin)))

# ==============================================================================
# ADDED: BASIN-AVERAGED WATER BALANCE SETUP
# ==============================================================================
wb_basin_avg <- colMeans(wb_basin, na.rm = TRUE)
basin_avg_spei_results <- list()
cat("✓ Basin-averaged water balance time series created\n")

# ==============================================================================
# WATER BALANCE ZERO DIAGNOSTICS
# ==============================================================================
cat("\n===== WATER BALANCE ZERO DIAGNOSTICS =====\n")
total_vals <- length(wb_basin) - sum(is.na(wb_basin))
zero_count <- sum(wb_basin == 0, na.rm = TRUE)
zero_pct <- 100 * zero_count / total_vals
cat(sprintf("Total basin water balance values: %d\n", total_vals))
cat(sprintf("Zero values count: %d (%.4f%%)\n", zero_count, zero_pct))
if (zero_pct > 1) cat("⚠ High proportion of zeros detected!\n")

# Per-month zero counts
month_numbers   <- as.integer(format(dates, "%m"))
month_names_vec <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
cat("\nZero counts by month:\n")
for (m in 1:12) {
  month_idx  <- which(month_numbers == m)
  if (length(month_idx) == 0) next
  month_vals  <- wb_basin[, month_idx]
  total_month <- length(month_vals) - sum(is.na(month_vals))
  zero_month  <- sum(month_vals == 0, na.rm = TRUE)
  pct_month   <- 100 * zero_month / total_month
  cat(sprintf("  %s: %d / %d (%.2f%%)\n", month_names_vec[m], zero_month, total_month, pct_month))
}
cat("=========================================\n\n")

# ==============================================================================
# PREPARE COMBINED SUMMARY FILE (single file for all scales)
# ==============================================================================
summary_file <- file.path(out_dir, "spei_all_scales_summary.txt")
cat("SPEI SUMMARY FOR ALL TIMESCALES (Variance-Aware Approach)\n",
    paste("Generated:", Sys.time()), "\n",
    "============================================================\n\n",
    file = summary_file, sep = "")

# ==============================================================================
# PARAMETRIC FITTING HELPER
# ==============================================================================
try_parametric <- function(values_m) {
  lmom <- tryCatch(lmoms(values_m, nmom = 4), error = function(e) NULL)
  if (is.null(lmom) || !are.lmom.valid(lmom)) return(NULL)
  
  dists <- list(
    GLO = list(par = tryCatch(parglo(lmom), error = function(e) NULL),
               cdf = cdfglo),
    PE3 = list(par = tryCatch(parpe3(lmom), error = function(e) NULL),
               cdf = cdfpe3),
    GEV = list(par = tryCatch(pargev(lmom), error = function(e) NULL),
               cdf = cdfgev)
  )
  
  best_p    <- -1
  best_prob <- NULL
  best_dist <- NULL
  
  for (d in names(dists)) {
    # FIX: Corrected list indexing dists[[d]]
    par <- dists[[d]]$par
    if (is.null(par)) next
    
    p_fitted <- tryCatch(dists[[d]]$cdf(values_m, par), error = function(e) NULL)
    if (is.null(p_fitted) || any(!is.finite(p_fitted))) next
    
    # KS goodness-of-fit test
    ks <- tryCatch(ks.test(values_m, function(q) dists[[d]]$cdf(q, par))$p.value,
                   error = function(e) 0)
    
    if (ks > 0.05 && ks > best_p) {   # satisfactory AND better than previous
      best_p    <- ks
      best_prob <- p_fitted
      best_dist <- d
    }
  }
  # Return list(prob, dist) if a distribution passed KS test, else NULL
  if (is.null(best_prob)) NULL else list(prob = best_prob, dist = best_dist)
}

# ==============================================================================
# VARIANCE-AWARE SPEI (Parametric-first with Gringorten fallback)
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
  method <- rep(NA_integer_, n)  # Track method: 1=GLO, 2=PE3, 3=GEV, 4=Gringorten, 5=Zero Var
  dist_code <- c(GLO = 1L, PE3 = 2L, GEV = 3L)
  
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
      method[idx] <- 5L  # Zero variance
      next
    }
    
    # CASE 2: Very low variance (< 0.001 mm²) → Gringorten directly (L-moments unstable)
    if (v0 < 0.001) {
      p <- (rank(values_m, ties.method = "average") - 0.44) / (length(values_m) + 0.12)
      p <- pmax(pmin(p, 1 - 1e-6), 1e-6)
      z[idx] <- qnorm(p)
      method[idx] <- 4L  # Gringorten (low variance)
      next
    }
    
    # CASE 3: Sufficient variance → try parametric (GLO, PE3, GEV); fall back to Gringorten
    fit <- try_parametric(values_m)
    if (is.null(fit)) {
      # Fallback: Gringorten empirical ranking
      p <- (rank(values_m, ties.method = "average") - 0.44) / (length(values_m) + 0.12)
      method[idx] <- 4L  # Gringorten fallback
    } else {
      p <- fit$prob
      method[idx] <- dist_code[fit$dist] # 1=GLO, 2=PE3, 3=GEV
    }
    p <- pmax(pmin(p, 1 - 1e-6), 1e-6)
    z[idx] <- qnorm(p)
  }
  list(z = z, method = method)  # Return both Z-scores and method
}

# ==============================================================================
# MAIN CALCULATION LOOP (All Timescales) - PARALLEL VERSION
# ==============================================================================
# ---- Set up parallel cluster ----
n_cores <- max(1L, parallel::detectCores() - 1L)   # leave 1 core for the OS
cat(sprintf("\n✓ Starting parallel cluster with %d cores\n", n_cores))
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)  # fix Windows temp-dir issue
cl <- makeCluster(n_cores)

# Export everything the worker function needs
clusterExport(cl, varlist = c("variance_aware_spei", "try_parametric",
                              "wb_basin", "month_numbers"),
              envir = environment())
clusterEvalQ(cl, { library(zoo); library(lmomco) })

scales <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24)

for (scale_idx in seq_along(scales)) {
  scale <- scales[scale_idx]
  cat(sprintf("\n===== SPEI-%d (Variance-Aware Calculation) =====\n", scale))
  
  # ==============================================================================
  # ADDED: CALCULATE BASIN-AVERAGED SPEI
  # ==============================================================================
  avg_spei_out <- variance_aware_spei(wb_basin_avg, month_numbers, scale)
  basin_avg_spei_results[[as.character(scale)]] <- avg_spei_out$z
  
  # Calculate SPEI for all basin pixels - PARALLEL
  clusterExport(cl, varlist = "scale", envir = environment())
  start_time <- Sys.time()
  
  # FIX: Pass row vector wb_basin[i, ] to function
  pixel_list <- parLapply(cl, seq_len(nrow(wb_basin)), function(i) {
    variance_aware_spei(wb_basin[i, ], month_numbers, scale)
  })
  
  end_time <- Sys.time()
  cat(sprintf("  Parallel processing done (%.1f min)\n",
              as.numeric(difftime(end_time, start_time, units = "mins"))))
  
  # Extract results
  res_basin     <- do.call(rbind, lapply(pixel_list, `[[`, "z"))
  method_basin  <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
  
  # REMOVED ORPHANED DIAGNOSTIC BLOCK (idx_diag undefined)
  
  # Reconstruct full raster matrices
  res_full <- matrix(NA_real_, nrow = nrow(wb_mat), ncol = ncol(wb_mat))
  res_full <- res_basin
  
  method_full <- rep(NA_integer_, nrow(wb_mat))
  method_first <- apply(method_basin, 1, function(x) {
    vals <- x[!is.na(x)]
    if (length(vals) == 0) return(NA_integer_)
    as.integer(names(which.max(table(vals))))
  })
  method_full <- method_first
  
  # ==============================================================================
  # OUTPUT 1: CSV FILES (one per calendar month)
  # ==============================================================================
  cat("  -> Saving CSV files (one per month)...\n")
  for (m in 1:12) {
    idx_m <- which(month_numbers == m)
    if (length(idx_m) == 0) next
    spei_subset  <- res_full[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    df           <- as.data.frame(spei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    coords       <- xyFromCell(wb, 1:nrow(spei_subset))
    df           <- cbind(lon = coords[, 1], lat = coords[, 2], df)
    csv_file     <- file.path(out_dir, sprintf("spei_%02d_month%02d_%s.csv",
                                               scale, m, month_names_vec[m]))
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
    spei_subset <- res_full[, idx_m, drop=FALSE]
    dates_subset <- dates[idx_m]
    df <- as.data.frame(spei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    # FIX: Corrected xyFromCell syntax
    coords <- xyFromCell(wb, 1:nrow(spei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    excel_data[[month_names_vec[m]]] <- df
  }
  xlsx_file <- file.path(out_dir, sprintf("spei_%02d_all_months.xlsx", scale))
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  Sys.sleep(0.3)
  write_xlsx(excel_data, xlsx_file)
  cat(sprintf("  ✓ Saved Excel file: %s\n", xlsx_file))
  
  # ==============================================================================
  # OUTPUT 3: DISTRIBUTION MAP (PNG)
  # ==============================================================================
  cat("  -> Creating SPEI distribution map with color legend...\n")
  # FIX: Corrected rast(wb) syntax
  spei_dist_raster  <- rast(wb[[1]])
  values(spei_dist_raster) <- method_full
  spei_dist_mapping  <- data.frame(
    code  = c(1, 2, 3, 4, 5),
    name  = c("GLO", "PE3", "GEV", "Gringorten", "Zero Var"),
    color = c("#1a9641", "#a6d96a", "#fdae61", "#d7191c", "#b0b0b0"),
    stringsAsFactors = FALSE
  )
  valid_pixels  <- sum(!is.na(method_full))
  empirical_pct <- 100 * sum(method_full == 4, na.rm = TRUE) / valid_pixels
  png_file <- file.path(out_dir, sprintf("spei_%02d_distribution_map.png", scale))
  png(png_file, width = 2000, height = 1000, res = 150)
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  par(mar = c(3, 2, 3, 1))
  plot(spei_dist_raster,
       col    = spei_dist_mapping$color,
       breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5),
       legend = FALSE,
       main   = sprintf("SPEI-%d: Dominant Parametric Distribution (modal across months)", scale),
       axes   = FALSE, box = FALSE)
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  mtext(sprintf("Valid pixels: %d | Gringorten (empirical): %.1f%%", valid_pixels, empirical_pct),
        side = 1, line = 1, cex = 0.75, col = "gray30")
  par(mar = c(3, 1, 3, 2))
  plot.new()
  dist_counts    <- table(factor(method_full,
                                 levels = spei_dist_mapping$code,
                                 labels = spei_dist_mapping$name))
  legend_entries <- sprintf("%s\n(%d px)", names(dist_counts), as.integer(dist_counts))
  legend("center",
         legend = legend_entries,
         fill   = spei_dist_mapping$color,
         title  = "Distribution\n(modal)",
         cex    = 0.90,
         y.intersp = 1.6,
         bg     = "white",
         bty    = "o",
         border = "gray50",
         xpd    = TRUE)
  dev.off()
  layout(1)   # reset layout
  cat(sprintf("  ✓ Saved SPEI distribution map: %s\n", png_file))
  
  # ==============================================================================
  # OUTPUT 3b: PER-MONTH METHOD MAPS (12 PNGs per scale)
  # ==============================================================================
  cat("  -> Creating per-month method maps (12 PNGs)...\n")
  for (m in 1:12) {
    idx_m <- which(month_numbers == m)
    if (length(idx_m) == 0) next
    method_month_basin <- method_basin[, idx_m, drop=FALSE]
    method_month_vec <- apply(method_month_basin, 1, function(x) {
      v <- x[!is.na(x)]
      if (length(v) == 0) NA_integer_ else v[1] # Take first available method for month
    })
    method_month_full <- rep(NA_integer_, nrow(wb_mat))
    method_month_full <- method_month_vec
    method_month_rast <- rast(wb[[1]])
    values(method_month_rast) <- method_month_full
    m_counts <- table(factor(method_month_full,
                             levels = spei_dist_mapping$code,
                             labels = spei_dist_mapping$name))
    present   <- as.integer(m_counts) > 0
    m_colors  <- spei_dist_mapping$color
    m_entries <- sprintf("%s\n(%d px)", names(m_counts), as.integer(m_counts))
    m_breaks  <- c(spei_dist_mapping$code - 0.5,
                   max(spei_dist_mapping$code) + 0.5)
    png_method <- file.path(out_dir,
                            sprintf("spei_%02d_method_map_month%02d_%s.png",
                                    scale, m, month_names_vec[m]))
    png(png_method, width = 2000, height = 1000, res = 150)
    layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
    par(mar = c(3, 2, 3, 1))
    plot(method_month_rast,
         col    = m_colors,
         breaks = m_breaks,
         legend = FALSE,
         main   = sprintf("SPEI-%d: Fitting Method — %s",
                          scale, month_names_vec[m]),
         axes   = FALSE, box = FALSE)
    if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
    valid_m   <- sum(!is.na(method_month_full))
    empir_m   <- 100 * sum(method_month_full == 4, na.rm = TRUE) / max(valid_m, 1)
    mtext(sprintf("Valid pixels: %d | Gringorten: %.1f%%", valid_m, empir_m),
          side = 1, line = 1, cex = 0.75, col = "gray30")
    par(mar = c(3, 1, 3, 2))
    plot.new()
    legend("center",
           legend  = m_entries,
           fill    = m_colors,
           title   = "Distribution",
           cex     = 0.90,
           y.intersp = 1.6,
           bg      = "white",
           bty     = "o",
           border  = "gray50",
           xpd     = TRUE)
    dev.off()
    layout(1)
  }
  cat(sprintf("  ✓ Saved 12 per-month method maps for SPEI-%d\n", scale))
  
  # ==============================================================================
  # OUTPUT 4: SUMMARY STATISTICS (APPENDED TO SINGLE FILE)
  # ==============================================================================
  cat("  -> Appending summary statistics to combined file...\n")
  
  # Calculate NA rate for summary
  na_rate_basin <- sum(is.na(wb_basin)) / length(wb_basin)
  
  summary_text  <- capture.output({
    cat(sprintf("\n---------- SPEI-%d SUMMARY ----------\n", scale))
    cat(sprintf("Calculation date: %s\n\n", Sys.time()))
    cat(sprintf("Grid dimensions: %d x %d\n", ncol(wb), nrow(wb)))
    cat(sprintf("Basin pixels: %d\n", sum(!is.na(wb_mat[, 1]))))
    cat(sprintf("Time period: %s to %s (%d months)\n",
                min(dates), max(dates), length(dates)))
    cat(sprintf("\nVariance Handling:\n"))
    cat("  • True zero variance  (var < 2.2e-16):      SPEI = 0 (neutral)\n")
    cat("  • Very low variance   (0 < var < 0.001 mm²): Gringorten empirical ranking\n")
    cat("  • Sufficient variance (var ≥ 0.001 mm²):    Parametric (best of GLO/PE3/GEV),\n")
    cat("                                                Gringorten fallback if KS p ≤ 0.05\n")
    cat(sprintf("\nNA Analysis:\n"))
    cat(sprintf("  Total NA values in basin: %.3f%%\n", na_rate_basin))
    cat(sprintf("  MSPEI compatibility: %s\n",
                ifelse(na_rate_basin < 0.5, "✓ READY (<0.5% NAs)", "⚠ CAUTION (>0.5% NAs)")))
    cat(sprintf("\nMethod Distribution (modal across calendar months per pixel):\n"))
    cat(sprintf("  GLO (Log-Logistic):  %d pixels (%.1f%%)\n",
                sum(method_full == 1, na.rm = TRUE),
                100 * sum(method_full == 1, na.rm = TRUE) / valid_pixels))
    cat(sprintf("  PE3 (Pearson III):   %d pixels (%.1f%%)\n",
                sum(method_full == 2, na.rm = TRUE),
                100 * sum(method_full == 2, na.rm = TRUE) / valid_pixels))
    cat(sprintf("  GEV:                 %d pixels (%.1f%%)\n",
                sum(method_full == 3, na.rm = TRUE),
                100 * sum(method_full == 3, na.rm = TRUE) / valid_pixels))
    cat(sprintf("  Gringorten fallback: %d pixels (%.1f%%)\n",
                sum(method_full == 4, na.rm = TRUE),
                100 * sum(method_full == 4, na.rm = TRUE) / valid_pixels))
    cat(sprintf("  Zero variance:       %d pixels (%.1f%%)\n",
                sum(method_full == 5, na.rm = TRUE),
                100 * sum(method_full == 5, na.rm = TRUE) / valid_pixels))
    
    # Find first 12 non-NA indices for drought freq check
    first_idx  <- head(which(!is.na(res_basin)), 12)
    if (length(first_idx) >= 12) {
      # FIX: Corrected indexing for first_spei
      first_spei <- res_basin[first_idx]
      cat(sprintf("\nDrought frequency by severity (first 12 months):\n"))
      cat(sprintf("  Exceptional (SPEI < -2.0):  %.1f%%\n", 100 * mean(first_spei < -2.0, na.rm = TRUE)))
      cat(sprintf("  Extreme     (SPEI < -1.6):  %.1f%%\n", 100 * mean(first_spei < -1.6, na.rm = TRUE)))
      cat(sprintf("  Severe      (SPEI < -1.3):  %.1f%%\n", 100 * mean(first_spei < -1.3, na.rm = TRUE)))
      cat(sprintf("  Moderate    (SPEI < -0.8):  %.1f%%\n", 100 * mean(first_spei < -0.8, na.rm = TRUE)))
    }
    cat("\n")
  })
  cat(summary_text, file = summary_file, append = TRUE, sep = "\n")
  
  # ==============================================================================
  # OUTPUT 5: NETCDF FILES (FOR MSPEI INPUT)
  # ==============================================================================
  # ==============================================================================
  # ==============================================================================
  cat("  -> Saving NetCDF files (for MSPEI input)...\n")
  
  # Create basin mask for value reconstruction
  basin_cells <- which(!is.na(values(wb)[,1]))  # Get cell indices of basin pixels
  
  for (m in 1:12) {
    idx_m <- which(month_numbers == m)
    if (length(idx_m) == 0) next
    
    # Create raster template with correct number of layers
    spei_rast <- rep(rast(wb[[1]]), length(idx_m))  # single-layer template × length(idx_m) layers
    values(spei_rast) <- res_full[, idx_m]    
    terra::time(spei_rast) <- dates[idx_m]
    
    nc_file <- file.path(out_dir, sprintf("spei_%02d_month%02d_%s.nc", 
                                          scale, m, month_names_vec[m]))
    suppressMessages(
      writeCDF(spei_rast, nc_file,
               varname = "spei",
               longname = sprintf("Standardized Precipitation Evapotranspiration Index (SPEI-%d)", scale),
               unit = "standardized_index",
               missval = -9999,
               overwrite = TRUE)
    )
  }
  cat(sprintf("  ✓ Saved 12 NetCDF files for SPEI-%d\n", scale))
}

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\n✓ Parallel cluster stopped\n")

# ==============================================================================
# ADDED: SAVE BASIN-AVERAGED SPEI TO CSV
# ==============================================================================
cat("\n===== SAVING BASIN-AVERAGED SPEI =====\n")

months_all <- as.integer(format(dates, "%m"))
years_all  <- as.integer(format(dates, "%Y"))
month_names_save <- c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec")

for (scale in scales) {
  spei_series <- basin_avg_spei_results[[as.character(scale)]]
  
  yrs <- sort(unique(years_all))
  df_out <- data.frame(Year = yrs)
  
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    yr_m  <- years_all[idx_m]
    val_m <- spei_series[idx_m]
    col   <- setNames(val_m, yr_m)
    df_out[[month_names_save[m]]] <- col[as.character(yrs)]
  }
  
  csv_file <- file.path(out_dir,
                        sprintf("spei_%02d_basin_averaged_by_month.csv", scale))
  write.csv(df_out, csv_file, row.names = FALSE, na = "")
  cat(sprintf("✓ Saved basin-averaged SPEI-%d (12 monthly series) to: %s\n", scale, csv_file))
}

cat("\n============================================================\n")
cat("SPEI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat(sprintf("\nCombined summary file: %s\n", summary_file))
cat("Variance-aware approach: Parametric (GLO/PE3/GEV) with Gringorten fallback.\n")
