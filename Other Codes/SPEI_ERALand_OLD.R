##############################################
# SPEI calculation (terra + SPEI) for Nechako River Basin
# OPTIMIZED VERSION with Progress Bars
# ERA5-Land MODA monthly mean → monthly totals → SPI
# Distribution Selection: RMSE vs empirical CDF
# + Probability clipping to avoid ±Inf, with reporting
#Not using any GIS file
##############################################

library(terra)
library(SPEI)
library(parallel)
library(zoo)
library(ncdf4)
library(writexl)
library(RColorBrewer)
library(pbapply)  # For progress bars with parallel processing

if (!requireNamespace("lmomco", quietly = TRUE)) {
  stop("Package 'lmomco' is required.")
}

setwd("D:/Nechako_Drought/")
out_dir <- "spei_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

diag_dir <- file.path(out_dir, "diagnostics")
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

cat("\n========================================\n")
cat("SPEI CALCULATION - OPTIMIZED VERSION\n")
cat("========================================\n\n")

precip_path <- "monthly_data_direct/total_precipitation_monthly.nc"
pet_path    <- "monthly_data_direct/potential_evapotranspiration_monthly.nc"

cat("Loading raster data...\n")
precip <- rast(precip_path)
pet    <- rast(pet_path)

cat("Checking CRS and extent alignment...\n")
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
if (exists("basin", inherits = FALSE) && !is.null(basin)) {
  basin <- tryCatch({
    if (!same.crs(basin, target_crs)) project(basin, target_crs) else basin
  }, error = function(e) {
    warning("Skipping basin projection: ", conditionMessage(e))
    basin
  })
}


extract_time_robust <- function(raster_obj, file_path, var_name) {
  
  cat(sprintf("Extracting time for %s...\n", var_name))
  
  # 1. Try terra time()
  tvals <- try(terra::time(raster_obj), silent = TRUE)
  if (!inherits(tvals, "try-error") && !all(is.na(tvals))) {
    cat(" ✓ Time extracted from terra::time()\n")
    return(as.Date(tvals))
  }
  
  # 2. Try NetCDF dimension or variable
  tryCatch({
    
    nc <- nc_open(file_path)
    on.exit(nc_close(nc))
    
    # Check for time as variable
    if ("time" %in% names(nc$var)) {
      time_data  <- ncvar_get(nc, "time")
      time_units <- ncatt_get(nc, "time", "units")$value
    } else if ("time" %in% names(nc$dim)) {
      time_data  <- nc$dim$time$vals
      time_units <- nc$dim$time$units
    } else {
      stop("No time variable or dimension found")
    }
    
    origin_str <- sub(".*since ", "", time_units)
    origin_date <- as.Date(origin_str)
    
    if (grepl("hours", time_units, ignore.case = TRUE)) {
      tvals <- origin_date + time_data / 24
    } else {
      tvals <- origin_date + time_data
    }
    
    cat(" ✓ Time extracted from NetCDF\n")
    return(as.Date(tvals))
    
  }, error = function(e) {
    cat(" ⚠ NetCDF time not found, reconstructing sequence\n")
  })
  
  # 3. Fallback: reconstruct monthly sequence
  n_layers <- nlyr(raster_obj)
  start_date <- as.Date("1950-01-01")
  tvals <- seq(start_date, by = "month", length.out = n_layers)
  
  cat(" ✓ Time reconstructed from layer count\n")
  return(tvals)
}

dates <- extract_time_robust(precip, precip_path, "Precip")

cat("\nCalculating water balance (P - PET)...\n")
wb <- precip - pet

# Check lengths match
stopifnot(length(dates) == terra::nlyr(wb))

# Correct setter
terra::time(wb) <- dates

# Extract values and clip extremes
cat("Extracting and clipping values...\n")
wb_vals <- terra::values(wb)
wb_vals[wb_vals >  300] <-  300
wb_vals[wb_vals < -300] <- -300
values(wb) <- wb_vals
cat("  ✓ Values clipped to [-300, 300]\n\n")

# OPTIMIZATION: Pre-compile rolling sum function using Rcpp if available
roll_sum_right <- function(x, k) {
  # Optimized version using filter for speed
  if (all(is.na(x))) return(rep(NA, length(x)))
  zoo::rollapply(x, k, sum, align = "right", fill = NA, partial = FALSE)
}

clip_prob <- function(p, eps = 1e-6)
  pmax(pmin(p, 1 - eps), eps)

# OPTIMIZATION: Cache L-moments calculation
spei_modified_nonseasonal <- function(x_agg, dist, eps = 1e-6) {
  idx <- which(is.finite(x_agg))
  if (length(idx) < 10) return(rep(NA_real_, length(x_agg)))
  samp <- x_agg[idx]
  
  # Calculate L-moments once
  lm <- lmomco::lmoms(samp)
  
  # Use tryCatch to handle fitting failures gracefully
  par <- tryCatch({
    switch(
      dist,
      "log-Logistic" = lmomco::parglo(lm),
      "GEV"          = lmomco::pargev(lm),
      "PearsonIII"   = lmomco::parpe3(lm)
    )
  }, error = function(e) NULL)
  
  if (is.null(par)) return(rep(NA_real_, length(x_agg)))
  
  cdf <- tryCatch({
    switch(
      dist,
      "log-Logistic" = lmomco::cdfglo(x_agg[idx], par),
      "GEV"          = lmomco::cdfgev(x_agg[idx], par),
      "PearsonIII"   = lmomco::cdfpe3(x_agg[idx], par)
    )
  }, error = function(e) NULL)
  
  if (is.null(cdf)) return(rep(NA_real_, length(x_agg)))
  
  z <- rep(NA_real_, length(x_agg))
  z[idx] <- qnorm(clip_prob(cdf, eps))
  z
}

# OPTIMIZATION: Vectorized ECDF calculation
spei_fun_best_nonseasonal <- function(v, scale, dates, eps_clip) {
  x_agg <- roll_sum_right(v, scale)
  
  # Skip if too many NAs
  if (sum(is.finite(x_agg)) < 30) {
    return(list(index = rep(NA_real_, length(x_agg)), 
                best_dist = NA, 
                rmse = NA))
  }
  
  dists <- c("log-Logistic", "GEV", "PearsonIII")
  res   <- lapply(dists, spei_modified_nonseasonal, x_agg = x_agg)
  
  # Pre-calculate ECDF once
  idx_finite <- is.finite(x_agg)
  if (sum(idx_finite) >= 10) {
    ecdf_fun <- ecdf(x_agg[idx_finite])
    emp_cdf <- ecdf_fun(x_agg[idx_finite])
    
    rmses <- sapply(res, function(z) {
      idx <- is.finite(z) & idx_finite
      if (sum(idx) < 10) return(Inf)
      sqrt(mean((pnorm(z[idx]) - emp_cdf)^2))
    })
  } else {
    rmses <- rep(Inf, length(dists))
  }
  
  best <- which.min(rmses)
  list(index = res[[best]], best_dist = dists[best], rmse = rmses[best])
}

wb_matrix <- values(wb, mat = TRUE)
cat(sprintf("Processing %d pixels across %d time steps\n", 
            nrow(wb_matrix), ncol(wb_matrix)))
cat(sprintf("Total data points: %s\n\n", 
            format(nrow(wb_matrix) * ncol(wb_matrix), big.mark = ",")))

scales <- c(1, 3, 6, 9, 12, 24)
plot_scales  <- c(1, 3, 6, 12)
plot_indices <- c(24, 241, 611, 906)  # Starting from Dec 1951 (index 24)

# Pre-create cluster once for better efficiency
n_cores <- max(1, detectCores() - 1)
cat(sprintf("Using %d CPU cores for parallel processing\n\n", n_cores))

cat("========================================\n")
cat("STARTING SPEI CALCULATIONS\n")
cat("========================================\n\n")

total_start_time <- Sys.time()
spei_nc_files <- list()  # NEW: store output NetCDF paths by scale
for (scale_idx in seq_along(scales)) {
  
  scale <- scales[scale_idx]
  scale_start_time <- Sys.time()
  
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("SCALE %d of %d: SPEI-%d months\n", 
              scale_idx, length(scales), scale))
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  
  # Create cluster
  cat(sprintf("Initializing cluster with %d cores...\n", n_cores))
  cl <- makeCluster(n_cores)
  
  clusterExport(
    cl,
    c("spei_fun_best_nonseasonal", "wb_matrix", "scale",
      "dates", "clip_prob", "roll_sum_right",
      "spei_modified_nonseasonal"),
    envir = environment()
  )
  clusterEvalQ(cl, library(lmomco))
  clusterEvalQ(cl, library(zoo))
  
  cat("Processing pixels with progress bar...\n")
  
  # Use pbapply for progress bar
  results_list <- pbapply(
    wb_matrix, 1,
    spei_fun_best_nonseasonal,
    scale = scale,
    dates = dates,
    eps_clip = 1e-6,
    cl = cl
  )
  
  stopCluster(cl)
  
  cat("Assembling results...\n")
  final_indices <- sapply(results_list, `[[`, "index")
  if (nrow(final_indices) != nrow(wb_matrix))
    final_indices <- t(final_indices)
  
  final_indices <- qnorm(
    clip_prob(pnorm(final_indices), 1e-6)
  )
  # ==============================================================================
  # WRITE SPEI NETCDF (one per scale) 
  # ==============================================================================
  cat("Writing SPEI NetCDF file...\n")
  
  spei_r <- rast(wb)  # template: same grid/time length as wb
  
  for (i in 1:nlyr(spei_r)) {
    values(spei_r[[i]]) <- final_indices[, i]
  }
  
  terra::time(spei_r) <- dates
  names(spei_r) <- sprintf("SPEI_%02d_%s", scale, format(dates, "%Y%m"))
  
  out_nc <- file.path(out_dir, sprintf("spei_%02dmonth.nc", scale))
  
  writeCDF(
    spei_r, out_nc,
    overwrite = TRUE,
    varname = "spei",
    unit = "unitless",
    prec = "float",
    missval = -9999
  )
  
  spei_nc_files[[as.character(scale)]] <- out_nc
  cat(sprintf(" ✓ Saved: %s\n", basename(out_nc)))
  
  # ==============================================================================
  # Create combined selected-layer PNG from multiple SPEI NetCDFs
  # Produces 1 PNG with 3 scales x 4 dates (12 panels)
  # ==============================================================================
  create_selected_layers_plot_group_from_nc <- function(nc_files, scales, diag_dir,
                                                        plot_indices = c(24, 241, 611, 906),
                                                        zlim = c(-2.5, 2.5),
                                                        prefix = "SPEI") {
    
    # Load rasters
    r_list <- lapply(scales, function(sc) {
      f <- nc_files[[as.character(sc)]]
      if (is.null(f) || !file.exists(f)) stop(sprintf("Missing NetCDF for %s-%d: %s", prefix, sc, f))
      r <- rast(f)
      # convert -9999 to NA
      r[r <= -9990] <- NA
      r
    })
    
    # Time axis from first raster
    dates <- terra::time(r_list[[1]])
    if (is.null(dates) || all(is.na(dates))) stop("No time axis found in NetCDF (terra::time is NA).")
    
    # Build 12-layer stack in row-major order: scale1(4 dates), scale2(4 dates), scale3(4 dates)
    layers <- list()
    nms <- character(0)
    
    for (sc_i in seq_along(scales)) {
      sc <- scales[sc_i]
      r  <- r_list[[sc_i]]
      
      for (idx in plot_indices) {
        layers[[length(layers) + 1]] <- r[[idx]]
        nms <- c(nms, sprintf("%s-%d | %s", prefix, sc, format(dates[idx], "%b %Y")))
      }
    }
    
    stack_r <- do.call(c, layers)
    names(stack_r) <- nms
    
    cols <- colorRampPalette(c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"))(100)
    
    png_file <- file.path(
      diag_dir,
      sprintf("%s_%s_selected_layers.png", prefix, paste(sprintf("%02d", scales), collapse = "_"))
    )
    
    png(png_file, width = 4200, height = 3000, res = 200)
    
    terra::plot(
      stack_r,
      nc   = length(plot_indices),  # 4 columns
      col  = cols,
      zlim = zlim,
      axes = FALSE,
      plg  = list(x = "right"),
      mar  = c(2, 2, 2, 6)
    )
    
    mtext(
      sprintf("Selected %s layers | Scales: %s", prefix, paste(scales, collapse = ", ")),
      side = 3, line = 1, cex = 1.2, font = 2
    )
    
    dev.off()
    cat(sprintf(" ✓ Combined selected-layers plot saved: %s\n", basename(png_file)))
    
    invisible(png_file)
  }
  ### >>> MAP PLOTTING SECTION <<<
  
  # if (scale %in% plot_scales) {
  #   
  #   cat("Generating selected layers maps...\n")
  #   
  #   # Reversed color palette: Blue (dry/low) -> White -> Red (wet/high)
  #   cols <- colorRampPalette(c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"))(100)
  #   zlim <- c(-2.5, 2.5)
  #   
  #   png(
  #     filename = file.path(
  #       diag_dir,
  #       sprintf("SPEI_%02d_selected_layers.png", scale)
  #     ),
  #     width = 1200, height = 3000, res = 200
  #   )
  #   
  #   layout(matrix(1:8, nrow = 4, byrow = TRUE),
  #          widths = c(4, 0.6))
  #   
  #   for (i in seq_along(plot_indices)) {
  #     
  #     idx <- plot_indices[i]
  #     r <- rast(wb[[1]])
  #     values(r) <- final_indices[, idx]
  #     
  #     par(mar = c(4, 4, 2, 1))
  #     plot(
  #       r,
  #       col = cols,
  #       zlim = zlim,
  #       main = sprintf(
  #         "SPEI-%d: %s (Layer %d)",
  #         scale,
  #         format(dates[idx], "%B %Y"),
  #         idx
  #       ),
  #       axes = TRUE,
  #       box = FALSE,
  #       legend = FALSE
  #     )
  #     
  #     par(mar = c(4, 1, 2, 3))
  #     image(
  #       x = 1,
  #       y = seq(zlim[1], zlim[2], length.out = length(cols)),
  #       z = matrix(seq(zlim[1], zlim[2], length.out = length(cols)), nrow = 1),
  #       col = cols,
  #       axes = FALSE
  #     )
  #     axis(4, las = 1, cex.axis = 0.7)
  #     box()
  #   }
  #   
  #   layout(1)
  #   dev.off()
  #   
  #   cat(sprintf("  ✓ Selected layers plot saved: SPEI_%02d_selected_layers.png\n", scale))
  # }
  # 
  # ==============================================================================
  # CREATE DISTRIBUTION VISUALIZATION (like SPI)
  # ==============================================================================
  cat("Creating distribution visualization...\n")
  
  # Extract distribution selection from results
  selected_dists <- sapply(results_list, `[[`, "best_dist")
  rmse_values <- sapply(results_list, `[[`, "rmse")
  
  # Create distribution raster
  dist_raster <- rast(wb[[1]])
  dist_codes <- match(selected_dists, c("log-Logistic", "GEV", "PearsonIII"))
  values(dist_raster) <- dist_codes
  
  # Create RMSE raster
  rmse_raster <- rast(wb[[1]])
  values(rmse_raster) <- rmse_values
  
  # Create distribution map PNG
  png_file <- file.path(diag_dir, sprintf("spei_%02dmonth_distribution_map.png", scale))
  png(png_file, width = 1600, height = 1200, res = 150)
  
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(0.65, 0.35))
  
  # Get extent and convert to matrices
  dist_mat <- as.matrix(dist_raster, wide = TRUE)
  rmse_mat <- as.matrix(rmse_raster, wide = TRUE)
  e <- ext(dist_raster)
  
  # 1. Distribution map
  par(mar = c(4, 4, 4, 2))
  image(x = seq(e[1], e[2], length.out = ncol(dist_mat)),
        y = seq(e[3], e[4], length.out = nrow(dist_mat)),
        z = t(dist_mat[nrow(dist_mat):1, ]),
        col = c("#2E86AB", "#A23B72", "#E6A23C"),  # log-Logistic, GEV, PearsonIII
        main = sprintf("Best Distribution Selection for SPEI-%d", scale),
        xlab = "Longitude", ylab = "Latitude",
        asp = 1)
  box()
  
  legend("topright", 
         legend = c("log-Logistic", "GEV", "Pearson III"),
         fill = c("#2E86AB", "#A23B72", "#E6A23C"),
         bty = "n", cex = 1.1)
  
  # 2. RMSE map - REVERSED PALETTE (warm colors for high RMSE)
  par(mar = c(4, 4, 3, 2))
  valid_rmse <- rmse_values[!is.na(rmse_values) & is.finite(rmse_values)]
  rmse_pal <- colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#FDDBC7", "#D6604D"))(100)
  
  # Normalize RMSE
  rmse_zlim <- c(0, quantile(valid_rmse, 0.95, na.rm = TRUE))
  rmse_mat_normalized <- (rmse_mat - rmse_zlim[1]) / (rmse_zlim[2] - rmse_zlim[1])
  rmse_mat_normalized[rmse_mat_normalized < 0] <- 0
  rmse_mat_normalized[rmse_mat_normalized > 1] <- 1
  
  image(x = seq(e[1], e[2], length.out = ncol(rmse_mat)),
        y = seq(e[3], e[4], length.out = nrow(rmse_mat)),
        z = t(rmse_mat_normalized[nrow(rmse_mat):1, ]),
        col = rmse_pal,
        main = "RMSE of Selected Distribution\n(Lower = Better Fit)",
        xlab = "Longitude", ylab = "Latitude",
        asp = 1)
  box()
  
  legend("topright", 
         legend = c("Excellent (Low RMSE)", "Good", "Poor (High RMSE)"),
         fill = c("#2166AC", "#F7F7F7", "#D6604D"),
         bty = "n", cex = 0.9)
  
  # 3. Statistical summary
  par(mar = c(1, 1, 1, 1))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
  
  # Calculate statistics
  total_cells <- length(selected_dists)
  loglog_cells <- sum(selected_dists == "log-Logistic", na.rm = TRUE)
  gev_cells <- sum(selected_dists == "GEV", na.rm = TRUE)
  pe3_cells <- sum(selected_dists == "PearsonIII", na.rm = TRUE)
  
  loglog_pct <- round(100 * loglog_cells / total_cells, 1)
  gev_pct <- round(100 * gev_cells / total_cells, 1)
  pe3_pct <- round(100 * pe3_cells / total_cells, 1)
  
  loglog_rmse <- mean(rmse_values[selected_dists == "log-Logistic"], na.rm = TRUE)
  gev_rmse <- mean(rmse_values[selected_dists == "GEV"], na.rm = TRUE)
  pe3_rmse <- mean(rmse_values[selected_dists == "PearsonIII"], na.rm = TRUE)
  
  summary_text <- sprintf(
    "DISTRIBUTION SELECTION SUMMARY (SPEI-%d)\n\nTotal grid cells: %d\n\nlog-Logistic: %d cells (%.1f%%)\n  Mean RMSE: %.4f\n\nGEV: %d cells (%.1f%%)\n  Mean RMSE: %.4f\n\nPearson III: %d cells (%.1f%%)\n  Mean RMSE: %.4f",
    scale, total_cells,
    loglog_cells, loglog_pct, loglog_rmse,
    gev_cells, gev_pct, gev_rmse,
    pe3_cells, pe3_pct, pe3_rmse
  )
  
  text(0.05, 0.95, summary_text, pos = 4, cex = 1.0, family = "mono", adj = c(0, 1))
  
  dev.off()
  cat(sprintf("  ✓ Distribution map saved: spei_%02dmonth_distribution_map.png\n", scale))
  
  # Create text report
  report_file <- file.path(diag_dir, sprintf("spei_%02dmonth_distribution_report.txt", scale))
  cat(
    "=====================================================================\n",
    sprintf("BEST DISTRIBUTION SELECTION REPORT (SPEI-%d)\n", scale),
    "=====================================================================\n\n",
    sprintf("Total Grid Cells: %d\n\n", total_cells),
    
    "DISTRIBUTION SELECTION FREQUENCY:\n",
    sprintf("  log-Logistic:  %5d cells (%6.2f%%) - RMSE: %.5f\n", loglog_cells, loglog_pct, loglog_rmse),
    sprintf("  GEV:           %5d cells (%6.2f%%) - RMSE: %.5f\n", gev_cells, gev_pct, gev_rmse),
    sprintf("  Pearson III:   %5d cells (%6.2f%%) - RMSE: %.5f\n\n", pe3_cells, pe3_pct, pe3_rmse),
    
    "=====================================================================\n",
    file = report_file
  )
  cat(sprintf("  ✓ Distribution report saved: spei_%02dmonth_distribution_report.txt\n", scale))
  
  ### <<< END MAP PLOTTING SECTION <<<
  
  scale_end_time <- Sys.time()
  scale_duration <- difftime(scale_end_time, scale_start_time, units = "mins")
  
  cat(sprintf("\n✓ SPEI-%d completed in %.1f minutes\n", scale, scale_duration))
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"))
  
  rm(results_list, final_indices)
  gc(verbose = FALSE)
}
# ==============================================================================
# AFTER ALL SCALES: CREATE 2 COMBINED PNGs (3 scales each)
# ==============================================================================
cat("\nCreating combined selected-layers PNGs...\n")

create_selected_layers_plot_group_from_nc(
  nc_files = spei_nc_files,
  scales   = c(1, 3, 6),
  diag_dir = diag_dir,
  plot_indices = c(24, 241, 611, 906),
  prefix = "SPEI"
)

create_selected_layers_plot_group_from_nc(
  nc_files = spei_nc_files,
  scales   = c(9, 12, 24),
  diag_dir = diag_dir,
  plot_indices = c(24, 241, 611, 906),
  prefix = "SPEI"
)
total_end_time <- Sys.time()
total_duration <- difftime(total_end_time, total_start_time, units = "mins")

cat("\n========================================\n")
cat("SPEI CALCULATION COMPLETE!\n")
cat("========================================\n")
cat(sprintf("Total processing time: %.1f minutes (%.1f hours)\n", 
            total_duration, total_duration / 60))
cat(sprintf("\n✓ All SPEI NetCDF files saved to: %s\n", normalizePath(out_dir)))
cat(sprintf("✓ All diagnostics, maps, and reports saved to: %s\n", diag_dir))
cat("========================================\n\n")