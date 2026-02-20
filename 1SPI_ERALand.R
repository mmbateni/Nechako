##############################################
# SPI CALCULATION - PARALLEL VERSION
# Uses parallel::parLapply - Works on Windows and Linux
##############################################

# ---- Libraries ----
library(terra)
library(SPEI)
library(ncdf4)
library(zoo)
library(writexl)
library(lmomco)
library(parallel) # Parallel processing

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")
out_dir <- "spi_results_seasonal"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Load Basin Boundary ----
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
cat("VARIANCE-AWARE SPI CALCULATION (PARALLEL MODE)\n")
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

# ---- Reproject to BC Albers ----
target_crs <- "EPSG:3005"

if (!same.crs(precip, target_crs)) {
  cat("✓ Reprojecting to BC Albers (EPSG:3005)...\n")
  precip <- project(precip, target_crs, method = "bilinear")
}

if (!same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
}

cat(sprintf("Grid: %d x %d | Time steps: %d\n", ncol(precip), nrow(precip), nlyr(precip)))
cat(sprintf("Time period: %s to %s\n", min(dates), max(dates)))

# ---- Basin masking ----
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
n_pixels <- nrow(precip_matrix)
cat(sprintf("Processing %d basin pixels...\n", n_pixels))

# ==============================================================================
# Helper Functions
# ==============================================================================
roll_sum_right <- function(x, k) {
  if (length(x) < k) return(rep(NA_real_, length(x)))
  zoo::rollapply(x, width = k, FUN = sum, align = "right", fill = NA_real_, partial = FALSE)
}

clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

# ==============================================================================
# VARIANCE-AWARE SPI FUNCTION
# ==============================================================================
variance_aware_spi <- function(v, scale, dates_vec, eps = 1e-6) {
  if (length(v) == 0 || all(is.na(v)) || is.null(v)) {
    return(list(spi = rep(NA_real_, length(dates_vec)), 
                method = rep(NA_integer_, length(dates_vec))))
  }
  
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  
  x_agg <- if (scale == 1) v_clean else roll_sum_right(v_clean, scale)
  
  if (length(x_agg) != length(dates_vec)) {
    x_agg <- c(rep(NA_real_, length(dates_vec) - length(x_agg)), x_agg[1:length(dates_vec)])
  }
  
  n <- length(x_agg)
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n)
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    idx <- which(mon == m & is.finite(x_agg))
    if (length(idx) < 5) next
    
    samp <- x_agg[idx]
    samp_var <- var(samp, na.rm = TRUE)
    
    # CASE 1: Zero variance
    if (!is.finite(samp_var) || samp_var < .Machine$double.eps) {
      z[idx] <- 0
      method_used[idx] <- 3
      next
    }
    
    # CASE 2: Low variance
    if (samp_var < 0.01) {
      p <- (rank(samp, ties.method = "average") - 0.44) / (length(samp) + 0.12)
      p <- clip_prob(p, eps = eps)
      z[idx] <- qnorm(p)
      method_used[idx] <- 2
      next
    }
    
    # CASE 3: Gamma fitting
    p0 <- mean(samp <= 0, na.rm = TRUE)
    samp_pos <- samp[samp > 0]
    
    if (length(samp_pos) >= 10) {
      lm <- try(lmomco::lmoms(samp_pos), silent = TRUE)
      if (!inherits(lm, "try-error") && !is.null(lm)) {
        par <- try(lmomco::pargam(lm), silent = TRUE)
        if (!inherits(par, "try-error") && !is.null(par)) {
          x_m <- x_agg[idx]
          p_m <- rep(NA_real_, length(x_m))
          
          pos_idx <- which(is.finite(x_m) & x_m > 0)
          if (length(pos_idx) > 0) {
            Fg <- try(lmomco::cdfgam(x_m[pos_idx], par), silent = TRUE)
            if (!inherits(Fg, "try-error")) {
              p_m[pos_idx] <- p0 + (1 - p0) * as.numeric(Fg)
            }
          }
          
          zero_idx <- which(is.finite(x_m) & x_m <= 0)
          if (length(zero_idx) > 0) p_m[zero_idx] <- p0
          
          p_m <- clip_prob(p_m, eps = eps)
          z[idx] <- qnorm(p_m)
          method_used[idx] <- 1
          next
        }
      }
    }
    
    # FALLBACK: Empirical
    p <- (rank(samp, ties.method = "average") - 0.44) / (length(samp) + 0.12)
    p <- clip_prob(p, eps = eps)
    z[idx] <- qnorm(p)
    method_used[idx] <- 2
  }
  
  finite_idx <- is.finite(z)
  z[finite_idx & z < -4.75] <- -4.75
  z[finite_idx & z >  4.75] <-  4.75
  
  list(spi = z, method = method_used)
}

# ==============================================================================
# MAIN CALCULATION LOOP - PARALLEL
# ==============================================================================
timescales <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24)

# ---- Set up parallel cluster (once, reused across all timescales) ----
n_cores <- max(1L, detectCores() - 1L)   # leave 1 core for the OS
cat(sprintf("\n✓ Starting parallel cluster with %d cores\n", n_cores))
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)  # fix Windows temp-dir issue
cl <- makeCluster(n_cores)

# Export data and functions needed by worker processes
eps_clip <- 1e-6
clusterExport(cl, varlist = c("variance_aware_spi", "roll_sum_right", "clip_prob",
                              "precip_matrix", "dates", "eps_clip"),
              envir = environment())
clusterEvalQ(cl, { library(zoo); library(lmomco) })
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# ---- Summary file setup ----
summary_combined_file <- file.path(out_dir, "spi_all_timescales_summary.txt")
all_summaries <- list()

for (sc in timescales) {
  cat(sprintf("\n===== SPI-%d (Sequential Calculation) =====\n", sc))
  cat(sprintf("Processing %d pixels in parallel... ", n_pixels))
  
  start_time <- Sys.time()
  
  # ---- PARALLEL PROCESSING ----
  clusterExport(cl, varlist = "sc", envir = environment())  # sc changes each iteration
  
  pixel_list <- parLapply(cl, seq_len(n_pixels), function(i) {
    tryCatch(
      variance_aware_spi(precip_matrix[i, ], sc, dates, eps = eps_clip),
      error = function(e) list(spi = rep(NA_real_, length(dates)),
                               method = rep(NA_integer_, length(dates)))
    )
  })
  
  # Convert to matrix
  spi_indices  <- do.call(rbind, lapply(pixel_list, `[[`, "spi"))
  method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
  
  end_time <- Sys.time()
  cat(sprintf("done (%.1f min)\n", as.numeric(difftime(end_time, start_time, units = "mins"))))
  
  if (is.null(dim(spi_indices)) || nrow(spi_indices) != nrow(precip_matrix)) {
    stop(sprintf("SPI-%d calculation failed: unexpected result dimensions", sc))
  }
  
  z <- spi_indices
  z[!is.finite(z)] <- NA_real_
  spi_indices <- z
  
  # Diagnostics
  clipped_low  <- sum(z < -4.7, na.rm = TRUE)
  clipped_high <- sum(z >  4.7, na.rm = TRUE)
  cat(sprintf("  Clipping: %d dry, %d wet\n", clipped_low, clipped_high))
  
  basin_mask_spi <- !is.na(precip_matrix[, 1])
  na_rate_basin <- 100 * mean(is.na(spi_indices[basin_mask_spi, ]))
  cat(sprintf("  NA rate (basin): %.3f%%\n", na_rate_basin))
  
  # Method distribution
  gamma_pct <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
  empirical_pct <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
  zero_var_pct <- 100 * sum(method_matrix == 3, na.rm = TRUE) / sum(!is.na(method_matrix))
  
  cat(sprintf("  Methods: Gamma=%.1f%%, Empirical=%.1f%%, Zero-var=%.1f%%\n", 
              gamma_pct, empirical_pct, zero_var_pct))
  
  # Store summary
  all_summaries[[as.character(sc)]] <- list(
    sc = sc,
    date = Sys.time(),
    na_rate = na_rate_basin,
    gamma_pct = gamma_pct,
    empirical_pct = empirical_pct,
    zero_var_pct = zero_var_pct,
    spi_indices = spi_indices,
    method_matrix = method_matrix
  )
  
  # ---- OUTPUTS ----
  months_all <- as.numeric(format(dates, "%m"))
  
  # CSV files
  cat("  -> Saving CSV files...\n")
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
  
  # Excel workbook
  cat("  -> Saving Excel workbook...\n")
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
  write_xlsx(excel_data, xlsx_file)
  
  # Distribution map
  cat("  -> Creating distribution map...\n")
  
  spi_dist_mapping <- data.frame(
    code = c(1, 2, 3),
    name = c("Gamma", "Empirical", "Zero-Variance"),
    color = c("#4575b4", "#d73027", "#91bfdb"),
    stringsAsFactors = FALSE
  )
  
  spi_dist_raster <- rast(precip[[1]])
  dist_sim <- apply(method_matrix, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_integer_)
    as.integer(names(sort(table(x), decreasing = TRUE))[1])
  })
  dist_sim[is.na(values(precip[[1]]))] <- NA
  values(spi_dist_raster) <- dist_sim
  
  png_file_spi <- file.path(out_dir, sprintf("spi_%02d_distribution_map.png", sc))
  png(png_file_spi, width = 1800, height = 1000, res = 150)
  
  # Split canvas: 80% map | 20% legend
  layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
  
  # ── Panel 1: Map ──────────────────────────────────────────
  plot(spi_dist_raster,
       col    = spi_dist_mapping$color,
       breaks = c(0.5, 1.5, 2.5, 3.5),
       legend = FALSE,
       main   = sprintf("SPI-%d: Variance-Aware Method", sc),
       axes   = FALSE, box = FALSE)
  
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  
  valid_pixels <- sum(!is.na(values(precip[[1]])))
  mtext(sprintf("Valid pixels: %d | Empirical: %.1f%%", valid_pixels, empirical_pct),
        side = 1, line = 1, cex = 0.75, col = "gray30")
  
  # ── Panel 2: Legend (empty plot, then draw legend into it) ─
  par(mar = c(0, 0, 0, 0))
  plot.new()   # blank panel
  
  dist_counts    <- table(factor(dist_sim[!is.na(dist_sim)], levels = 1:3,
                                 labels = spi_dist_mapping$name))
  legend_entries <- sprintf("%s (%d pixels)", names(dist_counts), dist_counts)
  
  legend("center",
         legend = legend_entries,
         fill   = spi_dist_mapping$color,
         title  = "Method",
         cex    = 0.95,
         bg     = "white",
         bty    = "o",
         border = "gray50")
  
  dev.off()
  layout(1)   # reset layout
  
  # NetCDF files
  cat("  -> Saving NetCDF files...\n")
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
  
  cat(sprintf("✓ SPI-%d complete\n", sc))
  gc()
}

# ==============================================================================
# WRITE SUMMARY FILE AT END
# ==============================================================================
cat("\n===== WRITING SUMMARY FILE =====\n")

cat("============================================================\n", 
    file = summary_combined_file)
cat("COMBINED SPI SUMMARY FOR ALL TIMESCALES\n", 
    file = summary_combined_file)
cat("Variance-Aware Approach (Gamma/Empirical fallback)\n", 
    file = summary_combined_file)
cat("============================================================\n\n", 
    file = summary_combined_file)

for (sc in timescales) {
  s <- all_summaries[[as.character(sc)]]
  
  cat(sprintf("\n\n========== SPI-%d ==========\n", s$sc), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("Calculation date: %s\n", s$date), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("Grid dimensions: %d x %d\n", ncol(precip), nrow(precip)), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("Basin pixels: %d\n", basin_pixels), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("Time period: %s to %s (%d months)\n", 
              min(dates), max(dates), length(dates)), 
      file = summary_combined_file, append = TRUE)
  
  cat("\nVariance Handling:\n", 
      file = summary_combined_file, append = TRUE)
  cat("  • True zero variance: SPI = 0 (neutral)\n", 
      file = summary_combined_file, append = TRUE)
  cat("  • Very low variance (< 0.01 mm²): Empirical ranking\n", 
      file = summary_combined_file, append = TRUE)
  cat("  • Sufficient variance (≥ 0.01 mm²): Gamma fitting\n", 
      file = summary_combined_file, append = TRUE)
  
  cat("\nMethod Distribution:\n", 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Gamma fitting: %.1f%%\n", s$gamma_pct), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Empirical ranking: %.1f%%\n", s$empirical_pct), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Zero-variance: %.1f%%\n", s$zero_var_pct), 
      file = summary_combined_file, append = TRUE)
  
  cat("\nNA Analysis:\n", 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Total NA values in basin: %.3f%%\n", s$na_rate), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  MSPEI compatibility: %s\n", 
              ifelse(s$na_rate < 0.5, "✓ READY (<0.5% NAs)", "⚠ CAUTION (>0.5% NAs)")), 
      file = summary_combined_file, append = TRUE)
  
  first_idx <- 1:min(12, ncol(s$spi_indices))
  first_spi <- s$spi_indices[, first_idx, drop = FALSE]
  
  cat("\nDrought frequency (first 12 months):\n", 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Exceptional (SPI < -2.0): %.1f%%\n", 
              100 * mean(first_spi < -2.0, na.rm = TRUE)), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Extreme     (SPI < -1.6): %.1f%%\n", 
              100 * mean(first_spi < -1.6, na.rm = TRUE)), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Severe      (SPI < -1.3): %.1f%%\n", 
              100 * mean(first_spi < -1.3, na.rm = TRUE)), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Moderate    (SPI < -0.8): %.1f%%\n", 
              100 * mean(first_spi < -0.8, na.rm = TRUE)), 
      file = summary_combined_file, append = TRUE)
}

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\n✓ Parallel cluster stopped\n")

cat("\n============================================================\n")
cat("SPI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat(sprintf("Combined summary file: %s\n", summary_combined_file))
cat("\n✓✓✓ READY FOR MSPEI CALCULATION ✓✓✓\n")