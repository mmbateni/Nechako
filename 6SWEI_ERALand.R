##############################################
# SWEI CALCULATION - PARALLEL VERSION (SPI-STYLE OUTPUTS)
# Uses parallel::parLapply - Works on Windows and Linux
# Output naming matches SPI_ERALand.R convention
##############################################
# ---- Libraries ----
  library(terra)
library(ncdf4)
library(zoo)
library(writexl)
library(parallel) # Parallel processing
# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")
out_dir <- "swei_results_seasonal"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# ---- Load Basin Boundary ----
basin_path <- "Spatial/nechakoBound_dissolve.shp"
basin <- if (file.exists(basin_path)) vect(basin_path) else NULL
if (!is.null(basin)) {
  cat("✓ Basin boundary loaded\n")
} else {
  stop("Basin boundary not found: ", basin_path)
}
# ---- Input files ----
swe_file <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
scf_file <- "monthly_data_direct/snow_cover_monthly.nc"
if (!file.exists(swe_file)) stop("SWE file not found: ", swe_file)
if (!file.exists(scf_file)) stop("SCF file not found: ", scf_file)

cat("\n============================================================\n")
cat("SEASONAL SWEI CALCULATION (PARALLEL MODE)\n")
cat("============================================================\n\n")

cat("===== READING NETCDF FILES =====\n")
swe <- rast(swe_file)
scf_stack <- rast(scf_file)

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
  if ("time" %in% names(nc$dim)) {
    tv <- ncvar_get(nc, "time")
    nc_close(nc)
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  nc_close(nc)
  n <- nlyr(raster_obj)
  start <- as.Date("1950-01-01")
  seq(start, by = "month", length.out = n)
}
dates_swe <- extract_time_dimension(swe, swe_file)
dates_scf <- extract_time_dimension(scf_stack, scf_file)
terra::time(swe) <- dates_swe
terra::time(scf_stack) <- dates_scf
dates <- dates_swe

# ---- Reproject to common grid ----
if (!compareGeom(scf_stack[[1]], swe[[1]], stopOnError = FALSE)) {
  cat("✓ Resampling SWE to match SCF grid...\n")
  swe <- resample(swe, scf_stack, method = "bilinear")
}
target_crs <- crs(scf_stack)
if (!is.null(basin) && !same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
}

# ---- Basin masking ----
cat("\n===== MASKING SWE TO BASIN BOUNDARY =====\n")
swe <- mask(swe, basin, inverse = FALSE, touches = TRUE)
basin_pixels <- global(swe[[1]], "notNA")$notNA
total_pixels <- ncell(swe)
cat(sprintf("✓ Basin masking complete: %d pixels (%.1f%% of raster)\n",
            basin_pixels, 100 * basin_pixels / total_pixels))

# ---- Convert SWE to mm if needed ----
swe_mean_sample <- global(swe[[1]], "mean", na.rm = TRUE)$mean
if (!is.na(swe_mean_sample) && swe_mean_sample < 10) {
  swe <- swe * 1000
  cat("✓ Converted SWE from meters to mm\n")
}

# ---- Convert SCF to fraction if needed ----
scf_max <- global(scf_stack[[1]], "max", na.rm = TRUE)$max
if (!is.na(scf_max) && scf_max > 1.5) {
  scf_stack <- scf_stack / 100
  cat("✓ Converted SCF from percent to fraction\n")
}

# ---- Prepare matrices ----
swe_matrix <- values(swe, mat = TRUE)
dates <- dates_swe
months_all <- as.integer(format(dates, "%m"))
years_all <- as.integer(format(dates, "%Y"))
n_pixels <- nrow(swe_matrix)
cat(sprintf("Processing %d basin pixels...\n", n_pixels))

# ==============================================================================
#   STEP 1: Apply 3-month moving average to RAW SWE values (PREPROCESSING)
# ==============================================================================
cat("\n===== STEP 1: Applying 3-month centered moving average to RAW SWE =====\n")
swe_smoothed <- matrix(NA, nrow = nrow(swe_matrix), ncol = ncol(swe_matrix))
for (i in 1:nrow(swe_matrix)) {
  swe_smoothed[i, ] <- zoo::rollapply(swe_matrix[i, ],
                                      width = 3,
                                      FUN = mean,
                                      align = "center",
                                      fill = NA,
                                      na.rm = TRUE)
}
cat("✓ 3-month moving average applied to SWE values\n")

# ==============================================================================
#   STEP 2: LOAD SCF MASK FOR k=3
# ==============================================================================
cat("\n===== STEP 2: Loading SCF mask for k=3 =====\n")
scf_diag_dir <- "scf_timescale_diagnostics"
target_k <- 3
scf_threshold <- 0.10
scf_diag_csv <- file.path(scf_diag_dir, "scf_timescale_final_recommendations.csv")

if (file.exists(scf_diag_csv)) {
  diag_df <- read.csv(scf_diag_csv, stringsAsFactors = FALSE)
  if ("timescale" %in% names(diag_df) && "chosen_threshold" %in% names(diag_df)) {
    k3_row <- diag_df[as.character(diag_df$timescale) == as.character(target_k), ]
    if (nrow(k3_row) > 0) {
      scf_threshold <- as.numeric(k3_row$chosen_threshold[1])
    }
  }
}

mask_file_pattern <- sprintf("scf_mask_k%d_*.nc", target_k)
mask_candidates <- list.files(scf_diag_dir, pattern = glob2rx(mask_file_pattern), full.names = TRUE)
if (length(mask_candidates) == 0) {
  stop(sprintf("CRITICAL ERROR: No SCF mask found for k=%d in %s", target_k, scf_diag_dir))
}
best_mask <- mask_candidates[[1]]
if (length(mask_candidates) > 1) {
  min_diff <- Inf
  for (f in mask_candidates) {
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
for (m in 1:12) {
  if (m <= nlyr(scf_mask_stack)) {
    mask_vals <- values(scf_mask_stack[[m]])
    scf_mask_list[[m]] <- !is.na(mask_vals) & (mask_vals == 1)
  } else {
    scf_mask_list[[m]] <- rep(FALSE, ncell(swe))
  }
}

# ==============================================================================
#   BASIN-LEVEL SCF MASK
#   For the basin-averaged SWEI path, determine which calendar months are
#   "snowy" at the basin scale.  A month is included when the majority of
#   basin pixels that have a valid SCF mask are flagged as snow-covered.
# ==============================================================================
cat("\n===== COMPUTING BASIN-LEVEL SCF MASK =====\n")
basin_scf_mask <- logical(12)
for (m in 1:12) {
  if (m <= nlyr(scf_mask_stack)) {
    mv <- as.vector(values(scf_mask_stack[[m]]))  # binary: 1 = include
    mv_valid <- mv[!is.na(mv)]
    # Month is "snowy enough" at basin scale when ≥50% of basin pixels are flagged
    basin_scf_mask[m] <- length(mv_valid) > 0 && mean(mv_valid, na.rm = TRUE) >= 0.5
  }
}
cat(sprintf("✓ Basin SCF mask: months included = %s\n",
            paste(month.abb[which(basin_scf_mask)], collapse = ", ")))

# ==============================================================================
#   BASIN-AVERAGED SWE SETUP
#   Correct order (matches SPI/SPEI approach):
#     1. Average RAW SWE across basin pixels   → basin time series
#     2. Apply 3-month rolling average          → smoothed basin series
#   (The old approach averaged the already-smoothed per-pixel matrix, which
#    gives the same numbers but couples the spatial and temporal smoothing
#    steps and makes the basin series dependent on per-pixel NA patterns.)
# ==============================================================================
cat("\n===== COMPUTING BASIN-AVERAGED SWE (raw → average → smooth) =====\n")
swe_basin_avg_raw      <- colMeans(swe_matrix, na.rm = TRUE)
swe_basin_avg_smoothed <- zoo::rollapply(swe_basin_avg_raw,
                                         width = 3,
                                         FUN   = mean,
                                         align = "center",
                                         fill  = NA,
                                         na.rm = TRUE)
cat(sprintf("✓ Basin-averaged smoothed SWE: %d time steps, %.1f%% non-NA\n",
            length(swe_basin_avg_smoothed),
            100 * mean(!is.na(swe_basin_avg_smoothed))))
basin_avg_swei_results <- list()

# ==============================================================================
#   Helper Functions
# ==============================================================================
clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

# ==============================================================================
#   SEASONAL SWEI FUNCTION (Gringorten)
# ==============================================================================
gringorten_swei_seasonal <- function(v_smoothed, dates_vec, scf_mask_list,
                                     cell_index = NULL,
                                     basin_scf_mask = NULL,
                                     eps = 1e-6) {
  # cell_index    : integer  → per-pixel mask lookup in scf_mask_list
  # basin_scf_mask: logical(12) → basin-level mask, one value per calendar month
  # If both are NULL the mask step is skipped entirely (not recommended)
  if (length(v_smoothed) == 0 || all(is.na(v_smoothed)) || is.null(v_smoothed)) {
    return(list(swei = rep(NA_real_, length(dates_vec)), method = rep(NA_integer_, length(dates_vec))))
  }
  v_clean <- v_smoothed
  v_clean[!is.finite(v_clean)] <- NA_real_
  n <- length(v_clean)
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n) # 1=Gringorten, 2=NA/Masked
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    idx <- which(mon == m & is.finite(v_clean))
    if (length(idx) < 20) next
    
    if (is.null(scf_mask_list) || length(scf_mask_list) < m) next
    # Check mask for this pixel (if cell_index provided)
    if (!is.null(cell_index)) {
      # Per-pixel path: check whether this cell is flagged as snowy for month m
      mask_vec <- scf_mask_list[[m]]
      if (length(mask_vec) < cell_index || !mask_vec[cell_index]) next
    } else if (!is.null(basin_scf_mask)) {
      # Basin-averaged path: use the pre-computed basin-level mask
      if (length(basin_scf_mask) < m || !basin_scf_mask[m]) next
    }
    # If both are NULL: no mask applied (original fallback; not recommended)
    
    samp <- v_clean[idx]
    zero_idx <- which(samp == 0 & !is.na(samp))
    if (length(zero_idx) > 0) {
      nonzero_vals <- samp[samp > 0 & !is.na(samp)]
      if (length(nonzero_vals) > 0) {
        tiny <- min(nonzero_vals) * 1e-5
        samp[zero_idx] <- tiny
      }
    }
    valid_idx <- which(!is.na(samp))
    if (length(valid_idx) < 3) next
    
    samp_valid <- samp[valid_idx]
    idx_valid <- idx[valid_idx]
    
    r <- rank(samp_valid, ties.method = "average")
    N <- length(samp_valid)
    p_val <- (r - 0.44) / (N + 0.12)
    p_val <- clip_prob(p_val, eps = eps)
    z_val <- qnorm(p_val)
    
    z[idx_valid] <- z_val
    method_used[idx_valid] <- 1
  }
  z[z < -4.75] <- -4.75
  z[z > 4.75] <- 4.75
  list(swei = z, method = method_used)
}

# ==============================================================================
#   MAIN CALCULATION LOOP - PARALLEL (Matches SPI Structure)
# ==============================================================================
timescales <- c(1) # SWEI typically 1-month output, but structured like SPI/SPEI/...
# ---- Set up parallel cluster ----
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\n✓ Starting parallel cluster with %d cores\n", n_cores))
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
cl <- makeCluster(n_cores)

clusterExport(cl, varlist = c("gringorten_swei_seasonal", "clip_prob",
                              "swe_smoothed", "dates", "scf_mask_list",
                              "basin_scf_mask"),
              envir = environment())
clusterEvalQ(cl, { library(zoo) })

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
# ---- Summary file setup ----
summary_combined_file <- file.path(out_dir, "swei_all_timescales_summary.txt")
all_summaries <- list()

for (sc in timescales) {
  cat(sprintf("\n===== SWEI-%d (Sequential Calculation) =====\n", sc))
  # ==============================================================================
  #   CALCULATE BASIN-AVERAGED SWEI
  # ==============================================================================
  avg_swei_out <- tryCatch(
    gringorten_swei_seasonal(swe_basin_avg_smoothed, dates, scf_mask_list,
                             cell_index = NULL, basin_scf_mask = basin_scf_mask,
                             eps = 1e-6),
    error = function(e) list(swei = rep(NA_real_, length(dates)), method = rep(NA_integer_, length(dates)))
  )
  basin_avg_swei_results[[as.character(sc)]] <- avg_swei_out$swei
  
  cat(sprintf("Processing %d pixels in parallel... ", n_pixels))
  start_time <- Sys.time()
  # ---- PARALLEL PROCESSING ----
  pixel_list <- parLapply(cl, seq_len(n_pixels), function(i) {
    tryCatch(
      gringorten_swei_seasonal(swe_smoothed[i, ], dates, scf_mask_list, cell_index = i, eps = 1e-6),
      error = function(e) list(swei = rep(NA_real_, length(dates)),
                               method = rep(NA_integer_, length(dates)))
    )
  })
  
  # Convert to matrix
  swei_indices <- do.call(rbind, lapply(pixel_list, `[[`, "swei"))
  method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
  
  # Diagnostics
  clipped_low   <- sum(swei_indices < -4.7, na.rm = TRUE)
  clipped_high  <- sum(swei_indices >  4.7, na.rm = TRUE)
  cat(sprintf("  Clipping: %d dry, %d wet\n", clipped_low, clipped_high))
  basin_mask_swei  <- !is.na(swe_smoothed)
  na_rate_basin  <- 100 * mean(is.na(swei_indices))
  cat(sprintf("  NA rate (basin): %.3f%%\n", na_rate_basin))
  
  # Method distribution
  gringorten_pct  <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
  masked_pct      <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
  
  cat(sprintf("  Methods: Gringorten=%.1f%%, Masked/NA=%.1f%%\n",
              gringorten_pct, masked_pct))
  
  # Store summary
  all_summaries[[as.character(sc)]] <- list(
    sc = sc,
    date = Sys.time(),
    na_rate = na_rate_basin,
    gringorten_pct = gringorten_pct,
    masked_pct = masked_pct,
    swei_indices = swei_indices,
    method_matrix = method_matrix
  )
  
  # ---- OUTPUTS (MATCHING SPI NAMING) ----
  months_all <- as.numeric(format(dates, "%m"))
  # CSV files
  cat("  -> Saving CSV files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    
    coords <- xyFromCell(swe, 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    csv_file <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.csv", sc, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE)
  }
  
  # Excel workbook
  cat("  -> Saving Excel workbook...\n")
  excel_data <- list()
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    
    coords <- xyFromCell(swe, 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    excel_data[[as.character(m)]] <- df
  }
  xlsx_file <- file.path(out_dir, sprintf("swei_%02d_all_months.xlsx", sc))
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  write_xlsx(excel_data, xlsx_file)
  
  # Distribution map
  cat("  -> Creating distribution map...\n")
  swei_dist_mapping <- data.frame(
    code = c(1, 2),
    name = c("Gringorten", "Masked/NA"),
    color = c("#4575b4", "#91bfdb"),
    stringsAsFactors = FALSE
  )
  swei_dist_raster <- rast(swe[[1]])
  dist_sim <- apply(method_matrix, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_integer_)
    as.integer(names(sort(table(x), decreasing = TRUE))[1])
  })
  dist_sim[is.na(dist_sim)] <- NA
  values(swei_dist_raster) <- dist_sim
  png_file_swei <- file.path(out_dir, sprintf("swei_%02d_distribution_map.png", sc))
  png(png_file_swei, width = 1800, height = 1000, res = 150)
  layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
  # ── Panel 1: Map ──────────────────────────────────────────
  plot(swei_dist_raster,
       col    = swei_dist_mapping$color,
       breaks = c(0.5, 1.5, 2.5),
       legend = FALSE,
       main   = sprintf("SWEI-%d: Seasonal Method", sc),
       axes   = FALSE, box = FALSE)
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  valid_pixels <- sum(!is.na(values(swe)))
  mtext(sprintf("Valid pixels: %d | Gringorten: %.1f%%", valid_pixels, gringorten_pct),
        side = 1, line = 1, cex = 0.75, col = "gray30")
  # ── Panel 2: Legend ─
  par(mar = c(0, 0, 0, 0))
  plot.new()
  dist_counts <- table(factor(dist_sim, levels = 1:2,
                              labels = swei_dist_mapping$name))
  legend_entries <- sprintf("%s (%d pixels)", names(dist_counts), dist_counts)
  legend("center",
         legend = legend_entries,
         fill   = swei_dist_mapping$color,
         title  = "Method",
         cex    = 0.95,
         bg     = "white",
         bty    = "o",
         border = "gray50")
  dev.off()
  layout(1)
  
  # NetCDF files (Monthly split like SPI)
  cat("  -> Saving NetCDF files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    swei_rast  <- rep(rast(swe[[1]]), length(idx_m))
    sub_matrix  <- swei_indices[, idx_m, drop = FALSE]
    if (nrow(sub_matrix) != ncell(swei_rast)) {
      stop(sprintf("Dimension mismatch: Matrix %d rows vs Raster %d cells.",
                   nrow(sub_matrix), ncell(swei_rast)))
    }
    values(swei_rast)  <- sub_matrix
    terra::time(swei_rast)  <- dates[idx_m]
    nc_file  <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.nc", sc, m, month_names[m]))
    writeCDF(swei_rast, nc_file,
             varname = "swei",
             longname = sprintf("Seasonal SWEI (SWEI-%d)", sc),
             unit = "standardized_index",
             missval = -9999,
             overwrite = TRUE)
  }
  cat(sprintf("✓ SWEI-%d complete\n", sc))
  gc()
}

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\n✓ Parallel cluster stopped\n")

# ==============================================================================
#   SAVE BASIN-AVERAGED SWEI TO CSV (Wide Format like SPI)
# ==============================================================================
cat("\n===== SAVING BASIN-AVERAGED SWEI =====\n")
months_all  <- as.integer(format(dates, "%m"))
years_all   <- as.integer(format(dates, "%Y"))
month_names_save  <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
for (sc in timescales) {
  swei_series <- basin_avg_swei_results[[as.character(sc)]]
  yrs <- sort(unique(years_all))
  df_out <- data.frame(Year = yrs)
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    yr_m  <- years_all[idx_m]
    val_m <- swei_series[idx_m]
    col   <- setNames(val_m, yr_m)
    df_out[[month_names_save[m]]] <- col[as.character(yrs)]
  }
  csv_file <- file.path(out_dir,
                        sprintf("swei_%02d_basin_averaged_by_month.csv", sc))
  write.csv(df_out, csv_file, row.names = FALSE, na = "")
  cat(sprintf("✓ Saved basin-averaged SWEI-%d (12 monthly series) to: %s\n", sc, csv_file))
}

# ==============================================================================
#   WRITE SUMMARY FILE AT END (Matches SPI Layout)
# ==============================================================================
cat("\n===== WRITING SUMMARY FILE =====\n")
cat(
  "============================================================\n",
  "COMBINED SWEI SUMMARY FOR ALL TIMESCALES\n",
  "Seasonal Approach (3-month Smoothed SWE + k=3 SCF Mask)\n",
  "============================================================\n\n",
  file = summary_combined_file, sep = ""
)
for (sc in timescales) {
  s  <- all_summaries[[as.character(sc)]]
  cat(sprintf("\n\n========== SWEI-%d ==========\n", s$sc),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("Calculation date: %s\n", s$date),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("Grid dimensions: %d x %d\n", ncol(swe), nrow(swe)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("Basin pixels: %d\n", basin_pixels),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("Time period: %s to %s (%d months)\n",
              min(dates), max(dates), length(dates)),
      file = summary_combined_file, append = TRUE)
  cat("\nMethodology:\n",
      file = summary_combined_file, append = TRUE)
  cat("  • 3-month centered moving average applied to RAW SWE\n",
      file = summary_combined_file, append = TRUE)
  cat("  • SCF mask derived from 3-month climatology (k=3)\n",
      file = summary_combined_file, append = TRUE)
  cat("  • Gringorten plotting position (non-parametric)\n",
      file = summary_combined_file, append = TRUE)
  cat("\nMethod Distribution:\n",
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Gringorten fitting: %.1f%%\n", s$gringorten_pct),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Masked/NA: %.1f%%\n", s$masked_pct),
      file = summary_combined_file, append = TRUE)
  cat("\nNA Analysis:\n",
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Total NA values in basin: %.3f%%\n", s$na_rate),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Compatibility: %s\n",
              ifelse(s$na_rate < 0.5, "✓ READY (<0.5% NAs)", "⚠ CAUTION (>0.5% NAs)")),
      file = summary_combined_file, append = TRUE)
  first_idx  <- 1:min(12, ncol(s$swei_indices))
  first_swei  <- s$swei_indices[, first_idx, drop = FALSE]
  cat("\nDrought frequency (first 12 months):\n",
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Exceptional (SWEI < -2.0): %.1f%%\n",
              100 * mean(first_swei < -2.0, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Extreme     (SWEI < -1.6): %.1f%%\n",
              100 * mean(first_swei < -1.6, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Severe      (SWEI < -1.3): %.1f%%\n",
              100 * mean(first_swei < -1.3, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
  cat(sprintf("  Moderate    (SWEI < -0.8): %.1f%%\n",
              100 * mean(first_swei < -0.8, na.rm = TRUE)),
      file = summary_combined_file, append = TRUE)
}
cat("\n============================================================\n")
cat("SWEI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat(sprintf("Combined summary file: %s\n", summary_combined_file))
cat("\n✓✓✓ READY FOR ANALYSIS ✓✓✓\n")