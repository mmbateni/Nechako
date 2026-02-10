##############################################
# SPI calculation (terra + SPEI) for Nechako River Basin
# OPTIMIZED VERSION with Progress Bars & Enhanced Diagnostics
# ERA5-Land MODA monthly mean → monthly totals → SPI
# Candidate distributions: Gamma + Weibull ONLY 
# Gamma via SPEI::spi (ub-pwm)
# Weibull via custom non-seasonal L-moments (lmomco)
# Selection: RMSE vs empirical CDF
##############################################

# ---- Libraries ----
library(terra)
library(SPEI)
library(parallel)
library(ncdf4)
library(zoo)
library(lmomco)
library(pbapply)  # For progress bars with parallel processing

# Optional packages for enhanced visualization
has_brewer <- requireNamespace("RColorBrewer", quietly = TRUE)
if (has_brewer) library(RColorBrewer)

# ---- Paths ----
setwd("D:/Nechako_Drought/")
out_dir <- "spi_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

diag_dir <- file.path(out_dir, "diagnostics")
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

cat("\n========================================\n")
cat("SPI CALCULATION - OPTIMIZED VERSION\n")
cat("========================================\n\n")

# ---- Input file ----
precip_file <- "monthly_data_direct/total_precipitation_monthly.nc"
if (!file.exists(precip_file)) stop("Input file not found")

cat("Loading precipitation data...\n")
full_raster <- rast(precip_file)

# ---- Time extraction ----
extract_time_dimension <- function(raster_obj, file_path) {
  cat("Extracting time dimension...\n")
  
  t <- time(raster_obj)
  if (!all(is.na(t))) {
    cat(" ✓ Time extracted from terra::time()\n")
    return(as.Date(t))
  }
  
  nc <- nc_open(file_path)
  if ("valid_time" %in% names(nc$var)) {
    tv <- ncvar_get(nc, "valid_time")
    nc_close(nc)
    cat(" ✓ Time extracted from NetCDF valid_time\n")
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  nc_close(nc)
  
  cat(" ⚠ Time not found, reconstructing sequence\n")
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
if (exists("basin", inherits = FALSE) && !is.null(basin)) {
  basin <- tryCatch({
    if (!same.crs(basin, target_crs)) project(basin, target_crs) else basin
  }, error = function(e) {
    warning("Skipping basin projection: ", conditionMessage(e))
    basin
  })
}

cat(sprintf("\n✓ Grid dimensions: %d x %d\n", ncol(precip), nrow(precip)))
cat(sprintf("✓ Time steps: %d (from %s to %s)\n", 
            nlyr(precip), 
            format(min(dates), "%Y-%m"),
            format(max(dates), "%Y-%m")))

# ---- MODA fix: monthly mean → monthly total ----
cat("\nConverting monthly mean → monthly total...\n")
first_of_month <- as.Date(format(dates, "%Y-%m-01"))
first_next_month <- seq(first_of_month[1], by = "month", length.out = length(dates) + 1)[-1]
days_in_month <- as.integer(first_next_month - first_of_month)

precip <- precip * days_in_month
cat("  ✓ Days-in-month scaling applied\n")

# ---- Units: m → mm ----
precip <- precip * 1000
mean_precip <- global(precip[[1]], "mean", na.rm = TRUE)$mean
cat(sprintf("  ✓ Converted to mm (example layer 1 mean: %.2f mm/month)\n\n", mean_precip))

# ---- Prepare matrix ----
cat("Extracting raster values to matrix...\n")
precip_matrix <- values(precip, mat = TRUE)
cat(sprintf("  ✓ Processing %d pixels across %d time steps\n", 
            nrow(precip_matrix), ncol(precip_matrix)))
cat(sprintf("  ✓ Total data points: %s\n\n", 
            format(nrow(precip_matrix) * ncol(precip_matrix), big.mark = ",")))

# ==============================================================================
# OPTIMIZED Helper Functions
# ==============================================================================
roll_sum_right <- function(x, k) {
  if (all(is.na(x))) return(rep(NA, length(x)))
  zoo::rollapply(x, width = k, FUN = sum, align = "right", fill = NA_real_, partial = FALSE)
}

clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

# OPTIMIZED: Weibull SPI with better error handling
spi_weibull_custom <- function(v, scale, dates_vec, eps = 1e-6) {
  v[!is.finite(v)] <- NA_real_
  x_agg <- roll_sum_right(v, scale)
  
  n <- length(v)
  z <- rep(NA_real_, n)
  
  idx <- which(is.finite(x_agg))
  if (length(idx) < 10) return(z)
  
  samp <- x_agg[idx]
  p0 <- mean(samp <= 0, na.rm = TRUE)
  samp_pos <- samp[samp > 0]
  
  if (length(samp_pos) < 10) return(z)
  
  # Enhanced error handling
  lm <- try(lmomco::lmoms(samp_pos), silent = TRUE)
  if (inherits(lm, "try-error") || is.null(lm)) return(z)
  
  par <- try(lmomco::parwei(lm), silent = TRUE)
  if (inherits(par, "try-error") || is.null(par)) return(z)
  
  x_m <- x_agg[idx]
  p_m <- rep(NA_real_, length(x_m))
  
  pos <- is.finite(x_m) & x_m > 0
  if (any(pos)) {
    Fw <- try(lmomco::cdfwei(x_m[pos], par), silent = TRUE)
    if (!inherits(Fw, "try-error")) {
      p_m[pos] <- p0 + (1 - p0) * as.numeric(Fw)
    }
  }
  
  zero <- is.finite(x_m) & x_m <= 0
  if (any(zero)) p_m[zero] <- p0
  
  p_m <- clip_prob(p_m, eps = eps)
  z[idx] <- qnorm(p_m)
  return(z)
}

# OPTIMIZED: SPI best distribution with vectorized ECDF
spi_fun_best <- function(v, scale, dates_vec, eps_clip = 1e-6) {
  # Early exit for invalid data
  if (all(is.na(v))) {
    return(list(index = rep(NA_real_, length(v)), best_dist = NA_character_, rmse = NA_real_))
  }
  
  v[!is.finite(v)] <- NA_real_
  x_agg <- roll_sum_right(v, scale)
  
  # Early exit if insufficient data
  if (sum(is.finite(x_agg)) < 20) {
    return(list(index = rep(NA_real_, length(v)), best_dist = NA_character_, rmse = NA_real_))
  }
  
  ts_v <- ts(v, start = 1, frequency = 1)
  
  dists <- c("Gamma", "Weibull")
  results <- list()
  rmse_vals <- setNames(rep(Inf, length(dists)), dists)
  
  # Gamma via SPEI::spi
  res_g <- try(suppressWarnings(
    SPEI::spi(ts_v, scale = scale, distribution = "Gamma", fit = "ub-pwm", na.rm = FALSE)
  ), silent = TRUE)
  
  if (!inherits(res_g, "try-error") && !is.null(res_g$fitted)) {
    z <- as.numeric(res_g$fitted)
    z[!is.finite(z)] <- NA_real_
    results[["Gamma"]] <- z
  }
  
  # Weibull custom
  results[["Weibull"]] <- spi_weibull_custom(v, scale, dates_vec, eps = eps_clip)
  
  # OPTIMIZATION: Pre-calculate ECDF once for RMSE calculation
  valid_idx <- is.finite(x_agg)
  if (sum(valid_idx) > 10) {
    ecdf_fun <- ecdf(x_agg[valid_idx])
    empirical_cdf <- ecdf_fun(x_agg[valid_idx])
    
    for (d in dists) {
      z <- results[[d]]
      if (is.null(z)) next
      
      valid_z_idx <- is.finite(z) & valid_idx
      if (sum(valid_z_idx) <= 10) next
      
      theoretical_cdf <- pnorm(z[valid_z_idx])
      rmse_vals[d] <- sqrt(mean((theoretical_cdf - empirical_cdf[valid_idx[valid_z_idx]])^2))
    }
  }
  
  best_dist <- names(which.min(rmse_vals))
  best_rmse <- rmse_vals[best_dist]
  
  if (!is.finite(best_rmse)) {
    return(list(index = rep(NA_real_, length(v)), best_dist = NA_character_, rmse = NA_real_))
  }
  
  list(index = results[[best_dist]], best_dist = best_dist, rmse = best_rmse)
}

# ==============================================================================
# ENHANCED VISUALIZATION FUNCTION FOR DISTRIBUTION SELECTION
# ==============================================================================
create_distribution_visualization <- function(selected_dists, rmse_values, template_raster, scale, diag_dir) {
  # Create distribution raster (1=Gamma, 2=Weibull)
  dist_raster <- rast(template_raster[[1]])
  dist_codes <- match(selected_dists, c("Gamma", "Weibull"))
  values(dist_raster) <- dist_codes
  
  # Create RMSE raster
  rmse_raster <- rast(template_raster[[1]])
  values(rmse_raster) <- rmse_values
  
  # NOTE: Raw distribution TIF is NOT saved (as requested)
  
  # Create styled PNG visualization with legend and context
  png_file <- file.path(diag_dir, sprintf("spi_%02dmonth_distribution_map.png", scale))
  png(png_file, width = 1600, height = 1200, res = 150)
  
  # Layout: Distribution map (top), RMSE map (bottom left), Statistics (bottom right)
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(0.65, 0.35))
  
  # Get raster extent and convert to matrices
  dist_mat <- as.matrix(dist_raster, wide = TRUE)
  rmse_mat <- as.matrix(rmse_raster, wide = TRUE)
  e <- ext(dist_raster)
  
  # 1. Distribution map with color legend
  par(mar = c(4, 4, 4, 2))
  image(x = seq(e[1], e[2], length.out = ncol(dist_mat)),
        y = seq(e[3], e[4], length.out = nrow(dist_mat)),
        z = t(dist_mat[nrow(dist_mat):1, ]),
        col = c("#2E86AB", "#A23B72"),  # Blue for Gamma, Magenta for Weibull
        main = sprintf("Best Distribution Selection for SPI-%d", scale),
        xlab = "Longitude", ylab = "Latitude",
        asp = 1)
  box()
  
  # Add legend manually for clarity
  legend("topright", 
         legend = c("Gamma Distribution", "Weibull Distribution"),
         fill = c("#2E86AB", "#A23B72"),
         bty = "n", cex = 1.1)
  
  # 2. RMSE map (fit quality) - REVERSED PALETTE (warm colors for high RMSE)
  par(mar = c(4, 4, 3, 2))
  valid_rmse <- rmse_values[!is.na(rmse_values) & is.finite(rmse_values)]
  # Reversed: Blue (low RMSE/good) -> Yellow -> Red (high RMSE/poor)
  rmse_pal <- colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#FDDBC7", "#D6604D"))(100)
  
  # Normalize RMSE for color mapping
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
  
  # 3. Statistical summary panel
  par(mar = c(1, 1, 1, 1))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
  
  # Calculate statistics
  total_cells <- length(selected_dists)
  gamma_cells <- sum(selected_dists == "Gamma", na.rm = TRUE)
  weibull_cells <- sum(selected_dists == "Weibull", na.rm = TRUE)
  gamma_pct <- round(100 * gamma_cells / total_cells, 1)
  weibull_pct <- round(100 * weibull_cells / total_cells, 1)
  
  gamma_rmse <- mean(rmse_values[selected_dists == "Gamma"], na.rm = TRUE)
  weibull_rmse <- mean(rmse_values[selected_dists == "Weibull"], na.rm = TRUE)
  
  # Create summary text
  summary_text <- sprintf(
    "DISTRIBUTION SELECTION SUMMARY (SPI-%d)\n\nTotal grid cells: %d\n\nGamma selected: %d cells (%.1f%%)\n  • Mean RMSE: %.4f\n\nWeibull selected: %d cells (%.1f%%)\n  • Mean RMSE: %.4f\n\nInterpretation:\n• Blue areas: Gamma distribution\n  provided better fit to precipitation\n  data (typically for moderate climates)\n\n• Magenta areas: Weibull distribution\n  better captured extreme values\n  (typically in drier regions)",
    scale,
    total_cells,
    gamma_cells, gamma_pct, gamma_rmse,
    weibull_cells, weibull_pct, weibull_rmse
  )
  
  text(0.05, 0.95, summary_text, pos = 4, cex = 1.1, family = "mono", adj = c(0, 1))
  
  dev.off()
  cat(sprintf("  ✓ Distribution map saved: %s\n", png_file))
  
  # Create comprehensive text report
  report_file <- file.path(diag_dir, sprintf("spi_%02dmonth_distribution_report.txt", scale))
  cat(
    "=====================================================================\n",
    sprintf("BEST DISTRIBUTION SELECTION REPORT (SPI-%d)\n", scale),
    "=====================================================================\n\n",
    sprintf("Analysis Period: %s to %s\n", format(min(dates), "%Y-%m"), format(max(dates), "%Y-%m")),
    sprintf("Total Grid Cells Analyzed: %d\n", total_cells),
    sprintf("Valid Cells (with sufficient data): %d\n\n", sum(!is.na(selected_dists))),
    
    "DISTRIBUTION SELECTION FREQUENCY:\n",
    sprintf("  Gamma Distribution:    %5d cells (%6.2f%%)\n", gamma_cells, gamma_pct),
    sprintf("  Weibull Distribution:  %5d cells (%6.2f%%)\n\n", weibull_cells, weibull_pct),
    
    "FIT QUALITY (RMSE - Lower is Better):\n",
    sprintf("  Gamma   (mean RMSE):   %.5f\n", gamma_rmse),
    sprintf("  Weibull (mean RMSE):   %.5f\n\n", weibull_rmse),
    
    "METHODOLOGY:\n",
    "  Distributions were selected based on RMSE comparison between\n",
    "  theoretical CDF (from fitted distribution) and empirical CDF\n",
    "  of aggregated precipitation data. Lower RMSE indicates better\n",
    "  representation of the actual precipitation distribution.\n\n",
    
    "INTERPRETATION GUIDANCE:\n",
    "  • Gamma distribution typically fits well in regions with\n",
    "    moderate precipitation regimes and fewer zero-rain days.\n",
    "  • Weibull distribution often better captures precipitation\n",
    "    patterns with higher skewness or more extreme events.\n",
    "  • Spatial patterns in distribution selection may reflect\n",
    "    climatic gradients across the study area.\n\n",
    
    "FILES GENERATED:\n",
    sprintf("  • Visual summary map:   spi_%02dmonth_distribution_map.png\n", scale),
    sprintf("  • This report:          spi_%02dmonth_distribution_report.txt\n", scale),
    "=====================================================================\n",
    file = report_file
  )
  cat(sprintf("  ✓ Distribution report saved: %s\n", report_file))
  
  # Return objects for potential further use
  list(distribution_raster = dist_raster, rmse_raster = rmse_raster)
}

# ==============================================================================
# FUNCTION TO CREATE SELECTED LAYERS VISUALIZATION
# ==============================================================================
# create_selected_layers_plot <- function(spi_raster, dates, scale, diag_dir) {
#   # Select 4 representative time steps starting from December 1951 (index 24)
#   n_layers <- nlyr(spi_raster)
#   plot_indices <- round(seq(24, n_layers, length.out = min(4, n_layers)))
#   
#   # Reversed color palette: Blue (low/dry) -> White -> Red (high/wet)
#   cols <- colorRampPalette(c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"))(100)
#   zlim <- c(-2.5, 2.5)
#   
#   png_file <- file.path(diag_dir, sprintf("SPI_%02d_selected_layers.png", scale))
#   png(png_file, width = 1200, height = 3000, res = 200)
#   
#   layout(matrix(1:8, nrow = 4, byrow = TRUE), widths = c(4, 0.6))
#   
#   for (i in seq_along(plot_indices)) {
#     idx <- plot_indices[i]
#     r <- spi_raster[[idx]]
#     
#     par(mar = c(4, 4, 2, 1))
#     plot(r, col = cols, zlim = zlim,
#          main = sprintf("SPI-%d: %s (Layer %d)", scale, format(dates[idx], "%B %Y"), idx),
#          axes = TRUE, legend = FALSE, box = FALSE)
#     
#     par(mar = c(4, 1, 2, 3))
#     image(x = 1, y = seq(zlim[1], zlim[2], length.out = length(cols)),
#           z = matrix(seq(zlim[1], zlim[2], length.out = length(cols)), nrow = 1),
#           col = cols, axes = FALSE)
#     axis(4, las = 1, cex.axis = 0.7)
#     box()
#   }
#   
#   layout(1)
#   dev.off()
#   cat(sprintf("  ✓ Selected layers plot saved: %s\n", basename(png_file)))
# }
# ==============================================================================
# FUNCTION: CREATE COMBINED SELECTED-LAYERS PLOT FOR MULTIPLE SCALES
# ==============================================================================
# ==============================================================================
# BETTER: Combined selected-layers plot using terra::plot() on a stacked raster
# Produces a clean n_scales x n_times panel plot + one legend.
# ==============================================================================
create_selected_layers_plot_group <- function(nc_files, scales, diag_dir,
                                              start_index = 24, n_time_panels = 4,
                                              zlim = c(-2.5, 2.5)) {
  
  # Load rasters
  spi_list <- lapply(scales, function(sc) {
    f <- nc_files[[as.character(sc)]]
    if (is.null(f) || !file.exists(f)) stop(sprintf("Missing NetCDF for SPI-%d: %s", sc, f))
    r <- rast(f)
    
    # Convert -9999 to NA if present
    r[r <= -9990] <- NA
    
    r
  })
  
  # Use time from first raster
  dates <- terra::time(spi_list[[1]])
  if (is.null(dates) || all(is.na(dates))) {
    stop("Time dimension missing in SPI NetCDF. terra::time(rast(nc)) returned NA.")
  }
  
  n_layers <- nlyr(spi_list[[1]])
  plot_indices <- round(seq(start_index, n_layers, length.out = min(n_time_panels, n_layers)))
  
  # Build a 12-layer stack: (SPI-1 at 4 dates, SPI-3 at 4 dates, SPI-6 at 4 dates)
  layers <- list()
  layer_names <- character(0)
  
  for (i in seq_along(scales)) {
    sc <- scales[i]
    r  <- spi_list[[i]]
    
    for (idx in plot_indices) {
      lyr <- r[[idx]]
      layers[[length(layers) + 1]] <- lyr
      layer_names <- c(layer_names, sprintf("SPI-%d | %s", sc, format(dates[idx], "%b %Y")))
    }
  }
  
  stack_r <- do.call(c, layers)
  names(stack_r) <- layer_names
  
  # Color palette (same style as your original selected-layer plots) [1](https://gounbc-my.sharepoint.com/personal/bateni_unbc_ca/Documents/Microsoft%20Copilot%20Chat%20Files/SPI_ERALand_OLD.txt)
  cols <- colorRampPalette(c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"))(100)
  
  # Output
  png_file <- file.path(
    diag_dir,
    sprintf("SPI_%s_selected_layers.png", paste(sprintf("%02d", scales), collapse = "_"))
  )
  
  # Device size tuned for 3 rows x 4 cols
  png(png_file, width = 4200, height = 3000, res = 200)
  
  # terra will lay these out properly with nc=4
  terra::plot(
    stack_r,
    nc   = length(plot_indices),      # 4 columns
    col  = cols,
    zlim = zlim,
    axes = FALSE,
    plg  = list(x = "right"),         # legend on the right
    mar  = c(2, 2, 2, 6)              # room for legend
  )
  
  mtext(
    sprintf("Selected SPI layers | Scales: %s", paste(scales, collapse = ", ")),
    side = 3, line = 1, outer = FALSE, cex = 1.2, font = 2
  )
  
  dev.off()
  cat(sprintf(" ✓ Combined selected-layers plot saved: %s\n", basename(png_file)))
  
  invisible(png_file)
}
# ==============================================================================
# MAIN LOOP WITH PROGRESS TRACKING
# ==============================================================================
n_cores <- max(1, detectCores() - 1)
cat(sprintf("Using %d CPU cores for parallel processing\n\n", n_cores))

timescales <- c(1, 3, 6, 9, 12, 24)
eps_clip <- 1e-6

distribution_summary_all <- list()

cat("========================================\n")
cat("STARTING SPI CALCULATIONS\n")
cat("========================================\n\n")

total_start_time <- Sys.time()
spi_nc_files <- list()   # NEW: store output NetCDF paths by scale
for (scale_idx in seq_along(timescales)) {
  
  sc <- timescales[scale_idx]
  scale_start_time <- Sys.time()
  
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("SCALE %d of %d: SPI-%d months (Gamma vs Weibull)\n", 
              scale_idx, length(timescales), sc))
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  
  cat(sprintf("Initializing cluster with %d cores...\n", n_cores))
  cl <- makeCluster(n_cores)
  clusterExport(
    cl,
    c("precip_matrix", "dates", "spi_fun_best", "spi_weibull_custom",
      "roll_sum_right", "clip_prob", "sc", "eps_clip"),
    envir = environment()
  )
  clusterEvalQ(cl, {
    library(SPEI)
    library(zoo)
    library(lmomco)
  })
  
  cat("Processing pixels with progress bar...\n")
  
  # Use pbapply for progress bar
  spi_results <- pbapply(precip_matrix, 1, function(v) {
    spi_fun_best(v, sc, dates, eps_clip = eps_clip)
  }, cl = cl)
  
  stopCluster(cl)
  
  cat("Assembling results...\n")
  spi_indices <- sapply(spi_results, function(x) x$index)
  selected_dists <- sapply(spi_results, function(x) x$best_dist)
  rmse_values <- sapply(spi_results, function(x) x$rmse)
  
  if (nrow(spi_indices) != nrow(precip_matrix)) spi_indices <- t(spi_indices)
  
  cat("\nDistribution Selection Summary:\n")
  dist_summary <- table(selected_dists, useNA = "ifany")
  print(dist_summary)
  
  distribution_summary_all[[as.character(sc)]] <- list(
    scale = sc,
    frequency = dist_summary,
    mean_rmse = tapply(rmse_values, selected_dists, mean, na.rm = TRUE)
  )
  
  # ==============================================================================
  # Probability clipping for FINAL SPI matrix + reporting
  # ==============================================================================
  cat("\nApplying probability clipping for safety...\n")
  
  z <- spi_indices
  z[!is.finite(z)] <- NA_real_
  
  p <- pnorm(z)
  p_clip <- clip_prob(p, eps = eps_clip)
  z_clip <- qnorm(p_clip)
  z_clip[is.na(z)] <- NA_real_
  
  changed <- which(is.finite(z) & is.finite(z_clip) & (z != z_clip), arr.ind = TRUE)
  if (nrow(changed) > 0) {
    cat(sprintf("  → Clipped %d index values (eps=%g)\n", nrow(changed), eps_clip))
    cells <- changed[, "row"]
    tcols <- changed[, "col"]
    coords <- terra::xyFromCell(precip[[1]], cells)
    
    clip_report <- data.frame(
      timescale = sc,
      cell = cells,
      x = coords[,1],
      y = coords[,2],
      time_col = tcols,
      date = dates[tcols],
      spi_raw = z[changed],
      spi_clipped = z_clip[changed],
      best_dist = selected_dists[cells],
      rmse = rmse_values[cells],
      stringsAsFactors = FALSE
    )
    
    out_csv <- file.path(diag_dir, sprintf("spi_%02d_prob_clipping_report.csv", sc))
    write.csv(clip_report, out_csv, row.names = FALSE)
    cat(sprintf("  ✓ Clipping report saved: %s\n", basename(out_csv)))
  } else {
    cat("  → No values required clipping adjustment\n")
  }
  
  spi_indices <- z_clip
  rm(z, p, p_clip, z_clip)
  
  # ==============================================================================
  # Write SPI raster safely
  # ==============================================================================
  cat("Writing SPI NetCDF file...\n")
  spi_r <- rast(precip)
  for (i in 1:nlyr(spi_r)) {
    values(spi_r[[i]]) <- spi_indices[, i]
  }
  
  terra::time(spi_r) <- dates
  names(spi_r) <- sprintf("SPI_%02d_%s", sc, format(dates, "%Y%m"))
  
  out_nc <- file.path(out_dir, sprintf("spi_%02dmonth.nc", sc))
  writeCDF(
    spi_r, out_nc,
    overwrite = TRUE,
    varname = "spi",
    unit = "unitless",
    prec = "float",
    missval = -9999
  )
  spi_nc_files[[as.character(sc)]] <- out_nc   # NEW: record the NetCDF path
  cat(sprintf("  ✓ Saved: %s\n", basename(out_nc)))
  
  # ==============================================================================
  # CREATE ENHANCED DISTRIBUTION VISUALIZATION
  # ==============================================================================
  cat("Creating distribution visualization...\n")
  vis_results <- create_distribution_visualization(
    selected_dists = selected_dists,
    rmse_values = rmse_values,
    template_raster = precip,
    scale = sc,
    diag_dir = diag_dir
  )
  
  # ==============================================================================
  # CREATE SELECTED LAYERS VISUALIZATION
  # ==============================================================================
  cat("Creating selected layers visualization...\n")
  #create_selected_layers_plot(spi_r, dates, sc, diag_dir)
  
  scale_end_time <- Sys.time()
  scale_duration <- difftime(scale_end_time, scale_start_time, units = "mins")
  
  cat(sprintf("\n✓ SPI-%d completed in %.1f minutes\n", sc, scale_duration))
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"))
  
  gc(verbose = FALSE)
}
# ==============================================================================
# AFTER ALL SCALES: CREATE 2 COMBINED PNGs
# ==============================================================================
cat("\nCreating combined selected-layers PNGs...\n")


create_selected_layers_plot_group(spi_nc_files, c(1, 3, 6),  diag_dir)
create_selected_layers_plot_group(spi_nc_files, c(9, 12, 24), diag_dir)

total_end_time <- Sys.time()
total_duration <- difftime(total_end_time, total_start_time, units = "mins")

cat("\n========================================\n")
cat("SPI CALCULATION COMPLETE!\n")
cat("========================================\n")
cat(sprintf("Total processing time: %.1f minutes (%.1f hours)\n", 
            total_duration, total_duration / 60))
cat(sprintf("\n✓ All SPI NetCDF files saved to: %s\n", normalizePath(out_dir)))
cat(sprintf("✓ All diagnostics, maps, and reports saved to: %s\n", diag_dir))
cat("========================================\n\n")
