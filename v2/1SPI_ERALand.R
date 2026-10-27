##############################################
# SPI CALCULATION - PARALLEL VERSION
# Shared-path version aligned with the observed FAO-56 PM + ERA5-Land PEV framework
##############################################

# ---- Libraries ----
rm(list = ls())

# ---- Working directory (framework rule) -------------------------------------
# Run every script from the project root: set the working directory there
# BEFORE running (setwd(), or open the RStudio project). All paths below are
# derived from getwd(); no path is hard-coded.
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")

library(terra)
library(SPEI)
library(ncdf4)
library(zoo)
library(writexl)
library(lmomco)
library(parallel) # Parallel processing

# ---- Paths / shared utilities ----
source(UTILS_FILE)

out_dir <- file.path(WD_PATH, "spi_results_seasonal")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Whole basin = dissolved union of the 7 authoritative subbasins
# (identical construction to 3SPEI_ERALand.R).
basin <- nechako_subbasin_union(load_nechako_subbasins())
cat("Basin boundary loaded (union of 7 subbasins)\n")

# ---- Input file ----
precip_file <- file.path(WD_PATH, "monthly_data_direct", "total_precipitation_monthly.nc")
if (!file.exists(precip_file)) stop("Input file not found: ", precip_file)


cat("\n============================================================\n")
cat("SPI CALCULATION (PARALLEL MODE)\n")
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

# -----------------------------------------------------------------------------
# PRECIPITATION UNIT / TEMPORAL-REPRESENTATION VERIFICATION
# -----------------------------------------------------------------------------
# Units are resolved by resolve_monthly_depth_unit() (DROUGHT_ANALYSIS_utils.R):
# the NetCDF `units` attribute AND the magnitude of the data are both checked.
# If they disagree the run stops; if the attribute is missing/unrecognised the
# unit is inferred from magnitude. Values are then converted ONCE to monthly
# totals (mm/month) with actual calendar days, so no double conversion is
# possible whichever of mm/day, m/day, mm/month or m/month the producer writes.
# -----------------------------------------------------------------------------
precip_dates_raw <- as.Date(time(precip))
if (is.null(precip_dates_raw) || all(is.na(precip_dates_raw))) {
  precip_dates_raw <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(precip))
}

precip_unit_info <- resolve_monthly_depth_unit(
  precip, kind = "precipitation", file = precip_file,
  label = "precipitation (SPI input)"
)
precip <- convert_monthly_depth(precip, from = precip_unit_info$unit,
                                to = "mm/month", dates = precip_dates_raw)
cat(sprintf("  Conversion applied: %s -> mm/month using actual calendar days.\n",
            precip_unit_info$unit))

july_idx_p <- which(as.integer(format(precip_dates_raw, "%m")) == 7L)[1L]
if (is.finite(july_idx_p)) {
  july_vals <- terra::values(precip[[july_idx_p]], mat = FALSE)
  july_vals <- july_vals[is.finite(july_vals)]
  cat(sprintf("  July converted mean       : %.3f mm/month\n",
              mean(july_vals, na.rm = TRUE)))
}
cat("===== PRECIPITATION UNIT CHECK PASSED =====\n\n")
# ---- Establish the canonical native ERA5-Land grid ----
target_crs <- terra::crs(precip)
if (is.na(target_crs) || !nzchar(target_crs)) stop("ERA5-Land precipitation raster has no CRS.")
basin <- terra::project(basin, target_crs)
if (!terra::same.crs(precip, basin)) stop("Basin and precipitation CRS are incompatible after transformation.")

common_grid <- save_common_analysis_grid(
  precip[[1]], basin = basin,
  source = "1SPI_ERALand.R canonical native ERA5-Land precipitation grid"
)
common_mask <- terra::unwrap(common_grid$basin_mask)

cat("Native ERA5-Land precipitation CRS:
")
print(target_crs)
cat(sprintf("Canonical grid: %d x %d = %d cells | res %.10f x %.10f | Time steps: %d
",
            ncol(precip), nrow(precip), ncell(precip),
            res(precip)[1], res(precip)[2], nlyr(precip)))
cat(sprintf("Time period: %s to %s
", min(dates), max(dates)))
cat(sprintf("Canonical basin cells: %d
", common_grid$basin_pixels))

cat("
===== MASKING PRECIPITATION TO CANONICAL BASIN MASK =====
")
precip <- terra::mask(precip, common_mask, maskvalues = 0, updatevalue = NA)
basin_pixels <- sum(!is.na(terra::values(precip[[1]])))
total_pixels <- terra::ncell(precip)
cat(sprintf("✓ Basin masking complete: %d pixels (%.1f%% of raster)
",
            basin_pixels, 100 * basin_pixels / total_pixels))
if (basin_pixels == 0L) stop("CRITICAL: Canonical basin mask produced zero valid ERA5-Land pixels.
")

spi_basin_bw <- basin_fractional_weights(precip[[1]], basin)
cat(sprintf("Example basin mean (layer 1, fractional-coverage weighted): %.2f mm/month
",
            basin_weighted_mean(terra::values(precip[[1]], mat = TRUE), spi_basin_bw)))

# ---- Prepare matrix: process only cells with valid basin precipitation ----
precip_values_all <- terra::values(precip, mat = TRUE)
valid_cells <- which(!is.na(terra::values(precip[[1]])))
if (!length(valid_cells))
  stop("No valid basin pixels available for SPI calculation.")
precip_matrix <- precip_values_all[valid_cells, , drop = FALSE]
n_pixels <- nrow(precip_matrix)
cat(sprintf("Processing %d valid basin pixels...\n", n_pixels))

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
    
    # CASE 3: zero-inflated gamma (L-moments) with Stagge et al. (2015)
    # centre-of-mass zero probability: p0 = n0/(n+1) for the mixture and
    # (n0+1)/(2(n+1)) for observed zeros. Shared with SSPI and SRI.
    p_m <- zero_inflated_gamma_prob(samp, samp, min_pos = 10L)
    if (!is.null(p_m) && all(is.finite(p_m))) {
      p_m <- clip_prob(p_m, eps = eps)
      z[idx] <- qnorm(p_m)
      method_used[idx] <- 1
      next
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
timescales <- c(1, 3, 6, 12, 24, 36)

# ---- Set up parallel cluster (once, reused across all timescales) ----
n_cores <- max(1L, detectCores() - 1L)   # leave 1 core for the OS
cat(sprintf("\n??? Starting parallel cluster with %d cores\n", n_cores))
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)  # fix Windows temp-dir issue
cl <- makeCluster(n_cores)
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
eps_clip <- 1e-6


# Export data and functions needed by worker processes
clusterExport(cl, varlist = c("variance_aware_spi", "roll_sum_right", "clip_prob",
                              "zero_inflated_gamma_prob", "stagge_zero_probabilities",
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
      {
        ans <- variance_aware_spi(precip_matrix[i, ], sc, dates, eps = eps_clip)
        ans$error <- NULL
        ans
      },
      error = function(e) list(
        spi = rep(NA_real_, length(dates)),
        method = rep(NA_integer_, length(dates)),
        error = conditionMessage(e)
      )
    )
  })
  
  worker_errors <- vapply(
    pixel_list,
    function(x) !is.null(x$error) && nzchar(x$error),
    logical(1)
  )
  if (any(worker_errors)) {
    first_error <- pixel_list[[which(worker_errors)[1L]]]$error
    stop(sprintf(
      "SPI-%d calculation failed in %d/%d pixels. First worker error: %s",
      sc, sum(worker_errors), length(worker_errors), first_error
    ))
  }
  
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
  clipped_low  <- sum(z <= -4.75, na.rm = TRUE)
  clipped_high <- sum(z >=  4.75, na.rm = TRUE)
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
  mon_all <- as.integer(format(dates, "%m"))
  cat(sprintf("  Case counts by calendar month (basin pixels, SPI-%d):\n", sc))
  cat(sprintf("  %-5s  %8s  %6s  %8s  %6s  %8s  %6s\n",
              "Month", "Case1(n)", "C1(%)", "Case2(n)", "C2(%)", "Case3(n)", "C3(%)"))
  for (mm in 1:12) {
    idx_mm  <- which(mon_all == mm)
    meth_mm <- as.integer(method_matrix[, idx_mm])
    meth_mm <- meth_mm[!is.na(meth_mm)]
    n_tot   <- length(meth_mm)
    if (n_tot == 0) next
    n1 <- sum(meth_mm == 1L); n2 <- sum(meth_mm == 2L); n3 <- sum(meth_mm == 3L)
    cat(sprintf("  %-5s  %8d  %5.1f%%  %8d  %5.1f%%  %8d  %5.1f%%\n",
                month.abb[mm],
                n1, 100*n1/n_tot, n2, 100*n2/n_tot, n3, 100*n3/n_tot))
  }
  
  # Store summary
  cbm_rows <- lapply(1:12, function(mm) {
    idx_mm  <- which(mon_all == mm)
    meth_mm <- as.integer(method_matrix[, idx_mm])
    meth_mm <- meth_mm[!is.na(meth_mm)]
    n_tot   <- length(meth_mm)
    if (n_tot == 0) return(NULL)
    n1 <- sum(meth_mm==1L); n2 <- sum(meth_mm==2L); n3 <- sum(meth_mm==3L)
    data.frame(month=month.abb[mm], n1=n1, pct1=100*n1/n_tot,
               n2=n2, pct2=100*n2/n_tot, n3=n3, pct3=100*n3/n_tot)
  })
  cbm <- do.call(rbind, Filter(Negate(is.null), cbm_rows))
  
  all_summaries[[as.character(sc)]] <- list(
    sc = sc,
    date = Sys.time(),
    na_rate = na_rate_basin,
    gamma_pct = gamma_pct,
    empirical_pct = empirical_pct,
    zero_var_pct = zero_var_pct,
    case_by_month = cbm,
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
    
    coords <- terra::xyFromCell(precip[[1]], valid_cells)
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
    
    coords <- terra::xyFromCell(precip[[1]], valid_cells)
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
  
  spi_dist_raster <- terra::rast(precip[[1]])
  # Use the dominant fitting method across the full time record rather than
  # an arbitrary 24th observation. Reconstruct the full raster from the
  # valid-cell subset used by the SPI calculation.
  dist_sim_valid <- apply(method_matrix, 1, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    tab <- table(x)
    as.numeric(names(tab)[which.max(tab)])
  })
  dist_sim <- rep(NA_real_, terra::ncell(precip))
  dist_sim[valid_cells] <- dist_sim_valid
  terra::values(spi_dist_raster) <- dist_sim
  
  png_file_spi <- file.path(out_dir, sprintf("spi_%02d_distribution_map.png", sc))
  png(png_file_spi, width = 1200, height = 800, res = 150)
  
  plot(spi_dist_raster,
       col = spi_dist_mapping$color,
       breaks = c(0.5, 1.5, 2.5, 3.5),
       legend = FALSE,
       main = sprintf("SPI-%d", sc),
       axes = FALSE, box = FALSE)
  
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  
  dist_counts <- table(factor(dist_sim[!is.na(dist_sim)], levels = 1:3, 
                              labels = spi_dist_mapping$name))
  legend_entries <- sprintf("%s (%d pixels)", names(dist_counts), dist_counts)
  
  legend("bottomright", legend = legend_entries, fill = spi_dist_mapping$color,
         title = "Method", cex = 0.85, bg = "white", bty = "o", border = "gray50")
  
  valid_pixels <- length(valid_cells)
  mtext(sprintf("Valid pixels: %d | Empirical: %.1f%%", valid_pixels, empirical_pct),
        side = 1, line = 0.5, cex = 0.75, col = "gray30")
  
  dev.off()
  
  # NetCDF files
  cat("  -> Saving NetCDF files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    
    spi_rast <- rast(precip[[1]])
    spi_rast <- rep(spi_rast, length(idx_m))
    spi_values_full <- matrix(NA_real_, nrow = terra::ncell(precip), ncol = length(idx_m))
    spi_values_full[valid_cells, ] <- spi_indices[, idx_m, drop = FALSE]
    terra::values(spi_rast) <- spi_values_full
    terra::time(spi_rast) <- dates[idx_m]
    
    nc_file <- file.path(out_dir, sprintf("spi_%02d_month%02d_%s.nc", sc, m, month_names[m]))
    writeCDF(spi_rast, nc_file,
             varname = "spi",
             longname = sprintf("Standardized Precipitation Index (SPI-%d)", sc),
             unit = "standardized_index",
             missval = -9999,
             overwrite = TRUE)
  }
  
  cat(sprintf("??? SPI-%d complete\n", sc))
  gc()
}

# ==============================================================================
# WRITE SUMMARY FILE AT END
# ==============================================================================
cat("\n===== WRITING SUMMARY FILE =====\n")

cat("============================================================\n", 
    file = summary_combined_file)
cat("COMBINED SPI SUMMARY FOR ALL TIMESCALES\n", 
    file = summary_combined_file, append = TRUE)
cat("Zero-inflated gamma (L-moments; Stagge et al. 2015 zero handling) with empirical fallback for low-variance months\n", 
    file = summary_combined_file, append = TRUE)
cat("============================================================\n\n", 
    file = summary_combined_file, append = TRUE)

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
  cat("  ??? True zero variance: SPI = 0 (neutral)\n", 
      file = summary_combined_file, append = TRUE)
  cat("  ??? Very low variance (< 0.01 mm??): Empirical ranking\n", 
      file = summary_combined_file, append = TRUE)
  cat("  ??? Sufficient variance (??? 0.01 mm??): zero-inflated gamma fitting (Stagge et al. 2015 zero probability)\n", 
      file = summary_combined_file, append = TRUE)
  
  cat("\nMethod Distribution:\n", 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Gamma fitting: %.1f%%\n", s$gamma_pct), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Empirical ranking: %.1f%%\n", s$empirical_pct), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Zero-variance: %.1f%%\n", s$zero_var_pct), 
      file = summary_combined_file, append = TRUE)
  if (!is.null(s$case_by_month)) {
    cat("\nCase counts by calendar month:\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  %-5s  %8s  %6s  %8s  %6s  %8s  %6s\n",
                "Month","Case1(n)","C1(%)","Case2(n)","C2(%)","Case3(n)","C3(%)"),
        file = summary_combined_file, append = TRUE)
    for (i in seq_len(nrow(s$case_by_month))) {
      r <- s$case_by_month[i, ]
      cat(sprintf("  %-5s  %8d  %5.1f%%  %8d  %5.1f%%  %8d  %5.1f%%\n",
                  r$month, r$n1, r$pct1, r$n2, r$pct2, r$n3, r$pct3),
          file = summary_combined_file, append = TRUE)
    }
  }
  
  cat("\nNA Analysis:\n", 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Total NA values in basin: %.3f%%\n", s$na_rate), 
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  MSPEI compatibility: %s\n", 
              ifelse(s$na_rate < 0.5, "??? READY (<0.5% NAs)", "??? CAUTION (>0.5% NAs)")), 
      file = summary_combined_file, append = TRUE)
  
  # First 12 months that carry SPI values: for SPI-k the first k-1 months are
  # NA by construction (accumulation warm-up), so start at the first valid month.
  valid_cols <- which(colSums(is.finite(s$spi_indices)) > 0)
  first_idx  <- utils::head(valid_cols, 12L)
  first_spi  <- s$spi_indices[, first_idx, drop = FALSE]
  cat(sprintf("\nDrought/Wet frequency (first 12 valid months, %s to %s):\n ",
              if (length(first_idx)) format(dates[first_idx[1]], "%Y-%m") else "NA",
              if (length(first_idx)) format(dates[utils::tail(first_idx, 1)], "%Y-%m") else "NA"),
      file = summary_combined_file, append = TRUE)
  cat(sprintf( "  W3 Extremely wet (SPI ???  2.0): %.1f%%\n ",
               100 * mean(first_spi >=  2.0, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf( "  W2 Severely wet  (1.5 ??? SPI < 2): %.1f%%\n ",
               100 * mean(first_spi >=  1.5 & first_spi < 2.0, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf( "  W1 Moderately wet (1 ??? SPI < 1.5): %.1f%%\n ",
               100 * mean(first_spi >=  1.0 & first_spi < 1.5, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf( "  N0 Normal         (-1 < SPI < 1):  %.1f%%\n ",
               100 * mean(first_spi > -1.0 & first_spi < 1.0, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf( "  D1 Moderately dry (-1.5 < SPI ??? -1): %.1f%%\n ",
               100 * mean(first_spi <= -1.0 & first_spi > -1.5, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf( "  D2 Severely dry   (-2 < SPI ??? -1.5): %.1f%%\n ",
               100 * mean(first_spi <= -1.5 & first_spi > -2.0, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf( "  D3 Extremely dry  (SPI ??? -2.0): %.1f%%\n ",
               100 * mean(first_spi <= -2.0, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
}

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\n??? Parallel cluster stopped\n")

cat("\n============================================================\n")
cat("SPI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat(sprintf("Combined summary file: %s\n", summary_combined_file))
cat("\n????????? READY FOR MSPEI CALCULATION ?????????\n")