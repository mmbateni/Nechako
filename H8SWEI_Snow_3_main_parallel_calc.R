##############################################
# 03_main_parallel_calc.R
# STAGE 3 of 4 — Main parallel calculation loop
#
# THIS IS ALMOST CERTAINLY THE SINGLE MOST TIME-CONSUMING AND MOST
# ERROR-PRONE STAGE (cluster worker crashes, memory pressure, per-pixel
# file-write failures). It is isolated here so that a crash never forces
# you to redo the (slow) bias-correction/regionalization stage or the
# (fast but fragile) data-loading stage.
#
# WHAT THIS SCRIPT DOES:
#   - Sets up the parallel cluster (makeCluster)
#   - Runs parLapply over all pixels to compute SWEI/SSPI-1
#   - Writes per-pixel/per-month .nc, .csv, and .xlsx outputs directly to
#     out_dir as it goes (these are NOT re-saved into the .rds checkpoint -
#     they are already safely on disk in out_dir)
#
# INPUT:  stage_checkpoints/stage2.rds  (from 02_regionalization_bias.R)
# OUTPUT: stage_checkpoints/stage3.rds  (small - just the basin-averaged
#         series + summary stats that 04 needs; the large gridded outputs
#         are already written to out_dir by this script directly)
##############################################

STAGE_DIR <- "stage_checkpoints"
stage2_path <- file.path(STAGE_DIR, "stage2.rds")
if (!file.exists(stage2_path)) {
  stop(sprintf(
    "Cannot find %s. Run 01_setup_and_data.R and 02_regionalization_bias.R first.",
    stage2_path
  ), call. = FALSE)
}

cat("============================================================\n")
cat("STAGE 3 - loading checkpoint from script 02\n")
cat("============================================================\n")
stage2_data <- readRDS(stage2_path)
list2env(stage2_data, envir = environment())
rm(stage2_data)

# ---- Rehydrate SpatRaster/SpatVector objects from terra's own files ----
# script 02 deliberately EXCLUDES swe/scf_stack/basin from stage2.rds (their
# terra external C++ pointers cannot survive saveRDS()/readRDS() - see the
# "external pointer is not valid" crash this replaced) and instead persists
# them as stage2_swe.tif / stage2_scf_stack.tif / stage2_basin.gpkg. Restore
# them here the same way script 02 restores script 01's versions. Done AFTER
# list2env() so these are the last word - nothing here gets clobbered by the
# (deliberately swe/scf_stack/basin-free) RDS payload.
swe <- terra::rast(file.path(STAGE_DIR, "stage2_swe.tif"))
terra::time(swe) <- dates_swe
scf_stack_path <- file.path(STAGE_DIR, "stage2_scf_stack.tif")
scf_stack <- if (file.exists(scf_stack_path)) {
  s <- terra::rast(scf_stack_path)
  terra::time(s) <- dates_scf
  s
} else NULL
basin <- terra::vect(file.path(STAGE_DIR, "stage2_basin.gpkg"))

cat(sprintf("[OK] Restored %d object(s) from %s\n", length(ls()), stage2_path))

# ==============================================================================
#   MAIN CALCULATION LOOP - PARALLEL 
# ==============================================================================
timescales <- if (scf_method == "2") c(1) else c(3)

# ---- Set up parallel execution safely ----
# Stage 2 has already shown that Windows PSOCK startup can fail/hang on this
# machine. Do not unconditionally request detectCores()-1 workers. Cap the
# worker count, catch startup failures, and fall back to sequential lapply().
if (!exists("ENABLE_PSOCK_PARALLEL", inherits = FALSE)) ENABLE_PSOCK_PARALLEL <- TRUE
if (!exists("MAX_PSOCK_WORKERS", inherits = FALSE)) MAX_PSOCK_WORKERS <- 4L

detected_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
n_cores <- max(1L, min(as.integer(MAX_PSOCK_WORKERS), as.integer(detected_cores) - 1L))

cat(sprintf("\n- CPU detection: %d logical core(s); Stage 3 worker cap: %d\n",
            detected_cores, MAX_PSOCK_WORKERS))

cl <- NULL
USE_PARALLEL_STAGE3 <- FALSE

if (isTRUE(ENABLE_PSOCK_PARALLEL) && n_cores >= 2L) {
  cat(sprintf("- Attempting Windows PSOCK cluster with %d worker(s)...\n", n_cores))
  dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
  cl <- tryCatch(
    parallel::makePSOCKcluster(
      n_cores,
      setup_timeout = 30,
      setup_strategy = "parallel",
      outfile = ""
    ),
    error = function(e) {
      warning(sprintf("Stage 3 PSOCK cluster creation failed: %s\nFalling back to sequential execution.",
                      conditionMessage(e)), call. = FALSE)
      NULL
    }
  )
  if (!is.null(cl)) {
    USE_PARALLEL_STAGE3 <- TRUE
    parallel::clusterSetRNGStream(cl, iseed = BASIN_SWEI_RANDOM_SEED)
    cat(sprintf("[OK] Stage 3 PSOCK cluster started with %d worker(s).\n", length(cl)))
  }
}

if (!USE_PARALLEL_STAGE3) {
  cl <- NULL
  cat("[NOTE] Stage 3 will run sequentially.\n",
      "       This avoids the Windows PSOCK socket problem but will be slower.\n", sep = "")
  set.seed(BASIN_SWEI_RANDOM_SEED)
}

if (USE_PARALLEL_STAGE3) {
  if (scf_method == "2" && REGIONALIZE_SSPI) {
    parallel::clusterExport(cl, varlist = c(
      "monthly_sspi", "monthly_sspi_regional", "clip_prob",
      "swe_matrix", "dates", "ref_idx", "region_id", "index_values",
      "C_matrix", "region_par", "discordant_mat", "MIN_PIXEL_GAMMA_SAMPLE"
    ), envir = environment())
    parallel::clusterEvalQ(cl, { library(zoo); library(lmomco) })
  } else if (scf_method == "2") {
    parallel::clusterExport(cl, varlist = c(
      "monthly_sspi", "clip_prob", "swe_matrix", "dates", "ref_idx",
      "MIN_PIXEL_GAMMA_SAMPLE"
    ), envir = environment())
    parallel::clusterEvalQ(cl, { library(zoo); library(lmomco) })
  } else {
    parallel::clusterExport(cl, varlist = c(
      "gringorten_swei_seasonal", "clip_prob", "dates", "scf_mask_list",
      "basin_scf_mask"
    ), envir = environment())
    parallel::clusterEvalQ(cl, { library(zoo) })
  }
}

stage3_lapply <- function(X, FUN) {
  if (USE_PARALLEL_STAGE3) parallel::parLapply(cl, X, FUN) else lapply(X, FUN)
}

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# ---- Summary file setup ----
summary_combined_file <- file.path(out_dir,
                                   if (scf_method == "2") "sspi_summary.txt" else "swei_all_timescales_summary.txt")
all_summaries <- list()

for (sc in timescales) {
  
  if (scf_method == "2") {
    # ============================================================================
    #   OPTION 2: SSPI-1 CALCULATION
    # ============================================================================
    cat(sprintf("\n===== SSPI-%d (Monthly, Parallel) =====\n", sc))
    cat(sprintf("Processing %d pixels in parallel (%s)... ", n_pixels,
                if (REGIONALIZE_SSPI) "regional pooled gamma" else "pixel-by-pixel gamma"))
    start_time <- Sys.time()
    
    if (REGIONALIZE_SSPI) {
      pixel_list <- stage3_lapply(seq_len(n_pixels), function(i) {
        tryCatch(
          monthly_sspi_regional(swe_matrix[i, ], dates, ref_idx, pixel_idx = i,
                                region_vec = region_id, index_values = index_values,
                                C = C_matrix, region_par = region_par,
                                discordant_mat = discordant_mat, eps = 1e-6),
          error = function(e) list(sspi = rep(NA_real_, length(dates)),
                                   method = rep(NA_integer_, length(dates)))
        )
      })
    } else {
      pixel_list <- stage3_lapply(seq_len(n_pixels), function(i) {
        tryCatch(
          monthly_sspi(swe_matrix[i, ], dates, ref_idx, eps = 1e-6),
          error = function(e) list(sspi = rep(NA_real_, length(dates)),
                                   method = rep(NA_integer_, length(dates)))
        )
      })
    }
    
    sspi_indices  <- do.call(rbind, lapply(pixel_list, `[[`, "sspi"))
    method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
    end_time <- Sys.time()
    cat(sprintf("done (%.1f min)\n", as.numeric(difftime(end_time, start_time, units = "mins"))))
    
    if (is.null(dim(sspi_indices)) || nrow(sspi_indices) != nrow(swe_matrix))
      stop("SSPI-1 calculation failed: unexpected result dimensions")
    sspi_indices[!is.finite(sspi_indices)] <- NA_real_
    
    clipped_low  <- sum(sspi_indices < -4.9, na.rm = TRUE)
    clipped_high <- sum(sspi_indices >  4.9, na.rm = TRUE)
    cat(sprintf("  Clipping: %d dry, %d wet\n", clipped_low, clipped_high))
    
    basin_mask_sspi <- !is.na(swe_matrix[, 1])
    na_rate_basin   <- 100 * mean(is.na(sspi_indices[basin_mask_sspi, ]))
    cat(sprintf("  NA rate (basin): %.3f%%\n", na_rate_basin))
    basin_avg_swei_results[[as.character(sc)]] <- colMeans(sspi_indices[basin_mask_sspi, , drop = FALSE], na.rm = TRUE)
    
    if (REGIONALIZE_SSPI) {
      n_valid_m   <- sum(!is.na(method_matrix))
      reg_gamma_pct  <- 100 * sum(method_matrix == 1, na.rm = TRUE) / n_valid_m
      reg_emp_pct    <- 100 * sum(method_matrix == 2, na.rm = TRUE) / n_valid_m
      local_gamma_pct <- 100 * sum(method_matrix == 3, na.rm = TRUE) / n_valid_m
      local_emp_pct   <- 100 * sum(method_matrix == 4, na.rm = TRUE) / n_valid_m
      local_zv_pct    <- 100 * sum(method_matrix == 5, na.rm = TRUE) / n_valid_m
      cat(sprintf("  Methods: Regional-Gamma=%.1f%%, Regional-Empirical=%.1f%%, Local-Gamma(fallback)=%.1f%%, Local-Empirical(fallback)=%.1f%%, Local-ZeroVar(fallback)=%.1f%%\n",
                  reg_gamma_pct, reg_emp_pct, local_gamma_pct, local_emp_pct, local_zv_pct))
      gamma_pct <- reg_gamma_pct; empirical_pct <- reg_emp_pct
      zero_var_pct <- local_zv_pct; excluded_pct <- local_gamma_pct + local_emp_pct
      all_summaries[[as.character(sc)]] <- list(
        sc = sc, date = Sys.time(), na_rate = na_rate_basin,
        gamma_pct = reg_gamma_pct, empirical_pct = reg_emp_pct,
        zero_var_pct = local_zv_pct, excluded_pct = local_gamma_pct + local_emp_pct,
        regional_gamma_pct = reg_gamma_pct, regional_empirical_pct = reg_emp_pct,
        local_gamma_pct = local_gamma_pct, local_empirical_pct = local_emp_pct,
        local_zerovar_pct = local_zv_pct,
        sspi_indices = sspi_indices, method_matrix = method_matrix
      )
    } else {
      gamma_pct      <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
      empirical_pct  <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
      zero_var_pct   <- 100 * sum(method_matrix == 3, na.rm = TRUE) / sum(!is.na(method_matrix))
      excluded_pct   <- 0
      cat(sprintf("  Methods: Gamma=%.1f%%, Empirical=%.1f%%, Zero-var=%.1f%%\n",
                  gamma_pct, empirical_pct, zero_var_pct))
      all_summaries[[as.character(sc)]] <- list(
        sc = sc, date = Sys.time(), na_rate = na_rate_basin,
        gamma_pct = gamma_pct, empirical_pct = empirical_pct,
        zero_var_pct = zero_var_pct, excluded_pct = excluded_pct,
        sspi_indices = sspi_indices, method_matrix = method_matrix
      )
    }
    
    months_all_sc <- as.numeric(format(dates, "%m"))
    
    cat("  -> Saving CSV files...\n")
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sub  <- sspi_indices[, idx_m, drop = FALSE]
      df   <- as.data.frame(sub)
      colnames(df) <- format(dates[idx_m], "%Y")
      coords <- xyFromCell(swe[[1]], 1:nrow(sub))
      df <- cbind(lon = coords[,1], lat = coords[,2], df)
      write.csv(df, file.path(out_dir, sprintf("sspi_%02d_month%02d_%s.csv", sc, m, month_names[m])), row.names = FALSE)
    }
    
    cat("  -> Saving Excel workbook...\n")
    excel_data <- list()
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sub  <- sspi_indices[, idx_m, drop = FALSE]
      df   <- as.data.frame(sub)
      colnames(df) <- format(dates[idx_m], "%Y")
      coords <- xyFromCell(swe[[1]], 1:nrow(sub))
      excel_data[[month_names[m]]] <- cbind(lon = coords[,1], lat = coords[,2], df)
    }
    xlsx_file <- file.path(out_dir, sprintf("sspi_%02d_all_months.xlsx", sc))
    if (file.exists(xlsx_file)) file.remove(xlsx_file)
    write_xlsx(excel_data, xlsx_file)
    
    cat("  -> Creating distribution map...\n")
    if (REGIONALIZE_SSPI) {
      # Code 2 "Regional Empirical" is retained only for backward-compatible
      # numbering; it is never assigned (no regional-empirical fallback is
      # used - see monthly_sspi_regional). It will always show 0 pixels below.
      method_mapping <- data.frame(
        code  = c(1, 2, 3, 4, 5),
        name  = c("Regional Gamma", "Regional Empirical (unused)", "Local Gamma (fallback)",
                  "Local Empirical (fallback)", "Local Zero-Var (fallback)"),
        color = c("#4575b4", "#74add1", "#fdae61", "#d73027", "#91bfdb"),
        stringsAsFactors = FALSE
      )
    } else {
      method_mapping <- data.frame(
        code  = c(1, 2, 3, 4),
        name  = c("Gamma", "Empirical", "Zero-Variance", "Excluded"),
        color = c("#4575b4", "#d73027", "#91bfdb", "#cccccc"),
        stringsAsFactors = FALSE
      )
    }
    method_raster <- rast(swe[[1]])
    dist_vals <- method_matrix[, 1]
    dist_vals[is.na(values(swe[[1]]))] <- NA
    values(method_raster) <- dist_vals
    n_cat <- nrow(method_mapping)
    png_file <- file.path(out_dir, sprintf("sspi_%02d_method_map.png", sc))
    png(png_file, width = 1200, height = 800, res = 150)
    plot(method_raster, col = method_mapping$color, breaks = seq(0.5, n_cat + 0.5, by = 1),
         legend = FALSE, main = sprintf("SSPI-%d: Fitting Method Distribution", sc),
         axes = FALSE, box = FALSE)
    if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
    dist_counts <- table(factor(dist_vals[!is.na(dist_vals)], levels = 1:n_cat,
                                labels = method_mapping$name))
    legend("bottomright", legend = sprintf("%s (%d pixels)", names(dist_counts), dist_counts),
           fill = method_mapping$color, title = "Method", cex = 0.8,
           bg = "white", bty = "o", border = "gray50")
    valid_pix <- sum(!is.na(values(swe[[1]])))
    mtext(sprintf("Valid pixels: %d | Gamma: %.1f%%", valid_pix, gamma_pct),
          side = 1, line = 0.5, cex = 0.7, col = "gray30")
    dev.off()
    
    cat("  -> Saving NetCDF files...\n")
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sspi_rast <- rep(rast(swe[[1]]), length(idx_m))
      values(sspi_rast) <- sspi_indices[, idx_m, drop = FALSE]
      terra::time(sspi_rast) <- dates[idx_m]
      writeCDF(sspi_rast, file.path(out_dir,
                                    sprintf("sspi_%02d_month%02d_%s.nc", sc, m, month_names[m])),
               varname = "sspi",
               longname = sprintf("Standardized SnowPack Index (SSPI-%d, monthly)", sc),
               unit = "standardized_index", missval = -9999, overwrite = TRUE)
    }
    cat(sprintf("[OK] SSPI-%d complete\n", sc))
    gc()
    next  
  }
  
  # ============================================================================
  #   OPTION 1: SWEI CALCULATION
  # ============================================================================
  cat(sprintf("\n===== SWEI-%d (%d-month window, Sequential Calculation) =====\n", sc, sc))
  
  swe_basin_avg_smoothed <- zoo::rollapply(swe_basin_avg_raw,
                                           width = sc,
                                           FUN   = sum,
                                           align = "right",
                                           fill  = NA,
                                           na.rm = FALSE)
  cat(sprintf("[OK] Basin-averaged %d-month smoothed SWE: %d steps, %.1f%% non-NA\n",
              sc, length(swe_basin_avg_smoothed),
              100 * mean(!is.na(swe_basin_avg_smoothed))))
  
  swe_smoothed <- matrix(NA, nrow = nrow(swe_matrix), ncol = ncol(swe_matrix))
  for (i in 1:nrow(swe_matrix)) {
    swe_smoothed[i, ] <- zoo::rollapply(swe_matrix[i, ],
                                        width = sc,
                                        FUN = sum,
                                        align = "right",
                                        fill = NA,
                                        na.rm = FALSE) 
  }
  cat(sprintf("[OK] %d-month rolling sum applied to SWE (SWEI preprocessing)\n", sc))
  
  if (USE_PARALLEL_STAGE3) {
    parallel::clusterExport(cl, varlist = "swe_smoothed", envir = environment())
  }
  
  avg_swei_out <- tryCatch(
    gringorten_swei_seasonal(swe_basin_avg_smoothed, dates, scf_mask_list,
                             cell_index = NULL, basin_scf_mask = basin_scf_mask,
                             eps = 1e-6),
    error = function(e) list(swei = rep(NA_real_, length(dates)), method = rep(NA_integer_, length(dates)))
  )
  basin_avg_swei_results[[as.character(sc)]] <- avg_swei_out$swei
  
  cat(sprintf("Processing %d pixels in parallel... ", n_pixels))
  start_time <- Sys.time()
  
  pixel_list <- stage3_lapply(seq_len(n_pixels), function(i) {
    tryCatch(
      gringorten_swei_seasonal(swe_smoothed[i, ], dates, scf_mask_list, cell_index = i, eps = 1e-6),
      error = function(e) list(swei = rep(NA_real_, length(dates)),
                               method = rep(NA_integer_, length(dates)))
    )
  })
  
  swei_indices <- do.call(rbind, lapply(pixel_list, `[[`, "swei"))
  method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
  
  extreme_dry  <- sum(swei_indices < -3.0, na.rm = TRUE)
  extreme_wet  <- sum(swei_indices >  3.0, na.rm = TRUE)
  cat(sprintf("  Extreme values (|SWEI| > 3): %d dry, %d wet (retained, no clipping)\n",
              extreme_dry, extreme_wet))
  
  basin_pix_rows      <- which(!is.na(swe_matrix[, 1]))
  
  # Protect against empty vectors generating NaN
  if (length(basin_pix_rows) == 0) {
    na_rate_grid <- NA
  } else {
    na_rate_grid <- 100 * mean(is.na(swei_indices[basin_pix_rows, , drop = FALSE]))
  }
  
  basin_mean_complete <- 100 * mean(!is.na(basin_avg_swei_results[[as.character(sc)]]))
  
  cat(sprintf("  Basin mean series:  %.1f%% complete (%.1f%% NA)\n",
              basin_mean_complete, 100 - basin_mean_complete))
  cat(sprintf("  Per-pixel gridded:  %.3f%% NA  (basin pixels only)\n", na_rate_grid))
  
  grid_status <- if (is.na(na_rate_grid)) {
    "UNKNOWN (0 valid basin pixels)" 
  } else if (na_rate_grid < 0.5) {
    "READY  (<0.5% NAs)"
  } else {
    "CAUTION (>0.5% NAs)"
  }
  
  cat(sprintf("  Gridded compatibility: %s\n", grid_status))
  
  if (!is.na(na_rate_grid) && na_rate_grid >= 0.5) {
    cat("  -  High per-pixel NA is expected: the 75%%-nonzero-years criterion\n")
    cat("     and SCF mask are applied pixel-by-pixel, which is much stricter\n")
    cat("     than the basin-mean aggregate implies. The basin mean series is\n")
    cat("     suitable for time-series analysis; use gridded files for spatial\n")
    cat("     analysis with awareness of actual spatial coverage.\n")
  }
  
  gringorten_pct  <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
  masked_pct      <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
  
  cat(sprintf("  Methods: Gringorten=%.1f%%, Masked/NA=%.1f%%\n",
              gringorten_pct, masked_pct))
  
  all_summaries[[as.character(sc)]] <- list(
    sc = sc, date = Sys.time(), na_rate = na_rate_grid, 
    basin_mean_complete = basin_mean_complete,
    gringorten_pct = gringorten_pct, masked_pct = masked_pct,
    swei_indices = swei_indices, method_matrix = method_matrix
  )
  
  months_all <- as.numeric(format(dates, "%m"))
  cat("  -> Saving CSV files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    csv_file <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.csv", sc, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE)
  }
  
  cat("  -> Saving Excel workbook...\n")
  excel_data <- list()
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    excel_data[[as.character(m)]] <- df
  }
  xlsx_file <- file.path(out_dir, sprintf("swei_%02d_all_months.xlsx", sc))
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  write_xlsx(excel_data, xlsx_file)
  
  cat("  -> Creating distribution map...\n")
  swei_dist_mapping <- data.frame(
    code = c(1, 2), name = c("Gringorten", "Masked/NA"),
    color = c("#4575b4", "#91bfdb"), stringsAsFactors = FALSE
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
  plot(swei_dist_raster, col = swei_dist_mapping$color, breaks = c(0.5, 1.5, 2.5),
       legend = FALSE, main = sprintf("SWEI-%d: Seasonal Method", sc), axes = FALSE, box = FALSE)
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  valid_pixels <- sum(!is.na(values(swe[[1]])))
  mtext(sprintf("Valid pixels: %d | Gringorten: %.1f%%", valid_pixels, gringorten_pct),
        side = 1, line = 1, cex = 0.75, col = "gray30")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  dist_counts <- table(factor(dist_sim, levels = 1:2, labels = swei_dist_mapping$name))
  legend_entries <- sprintf("%s (%d pixels)", names(dist_counts), dist_counts)
  legend("center", legend = legend_entries, fill = swei_dist_mapping$color,
         title = "Method", cex = 0.95, bg = "white", bty = "o", border = "gray50")
  dev.off()
  layout(1)
  
  cat("  -> Saving NetCDF files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
    swei_rast  <- rep(rast(swe[[1]]), length(idx_m))
    sub_matrix  <- swei_indices[, idx_m, drop = FALSE]
    if (nrow(sub_matrix) != ncell(swei_rast)) {
      stop(sprintf("Dimension mismatch: Matrix %d rows vs Raster %d cells.",
                   nrow(sub_matrix), ncell(swei_rast)))
    }
    values(swei_rast)  <- sub_matrix
    terra::time(swei_rast)  <- dates[idx_m]
    nc_file  <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.nc", sc, m, month_names[m]))
    writeCDF(swei_rast, nc_file, varname = "swei",
             longname = sprintf("Seasonal SWEI (SWEI-%d)", sc),
             unit = "standardized_index", missval = -9999, overwrite = TRUE)
  }
  cat(sprintf("[OK] SWEI-%d complete\n", sc))
  gc()
}

# ---- Shut down parallel cluster ----
if (USE_PARALLEL_STAGE3 && !is.null(cl)) {
  try(parallel::stopCluster(cl), silent = TRUE)
  cl <- NULL
  cat("\n- Parallel cluster stopped\n")
} else {
  cat("\n- Stage 3 sequential execution complete\n")
}
gc()


##############################################
# ---- STAGE 3 CHECKPOINT: save everything for script 04 ----
##############################################
# NOTE: the large per-pixel gridded outputs (.nc/.csv/.xlsx, one set per
# month) have ALREADY been written directly to out_dir by the loop above -
# they are not duplicated into this .rds. What's saved here is the
# workspace needed for script 04's basin-averaged CSVs, time-series plots,
# and the final summary/reporting text files (all_summaries and
# basin_avg_swei_results are the main - and still comparatively small -
# objects; at basin-grid scale this whole checkpoint is normally just a
# few MB to a few tens of MB).
cat("\n============================================================\n")
cat("STAGE 3 COMPLETE - saving workspace checkpoint for script 04\n")
cat("============================================================\n")

# ---- Persist SpatRaster/SpatVector objects via terra's own I/O ----
# Same reason as scripts 01/02's checkpoints: saveRDS() cannot round-trip
# terra's external C++ pointers ("external pointer is not valid" at save
# time). swe/scf_stack/basin were rehydrated from stage2's .tif/.gpkg files
# at the top of this script and are unchanged here, so simply re-copy those
# same files forward as this stage's checkpoint. EXCLUDE them from the RDS
# payload below so the RDS write itself doesn't touch a live pointer.
file.copy(file.path(STAGE_DIR, "stage2_swe.tif"),
          file.path(STAGE_DIR, "stage3_swe.tif"), overwrite = TRUE)
if (!is.null(scf_stack)) {
  file.copy(file.path(STAGE_DIR, "stage2_scf_stack.tif"),
            file.path(STAGE_DIR, "stage3_scf_stack.tif"), overwrite = TRUE)
}
if (!is.null(basin)) {
  file.copy(file.path(STAGE_DIR, "stage2_basin.gpkg"),
            file.path(STAGE_DIR, "stage3_basin.gpkg"), overwrite = TRUE)
}

# Auto-detect every remaining terra Spat* object before serialization.
# This prevents a newly introduced SpatRaster/SpatVector from silently entering
# the RDS payload and reproducing the external-pointer failure seen in Stage 2.
spat_leftover <- Filter(function(nm) {
  obj <- tryCatch(get(nm, envir = environment(), inherits = FALSE), error = function(e) NULL)
  inherits(obj, c("SpatRaster", "SpatVector", "SpatRasterDataset", "SpatRasterCollection"))
}, ls(all.names = TRUE, envir = environment()))
if (length(spat_leftover) > 0L) {
  warning("Excluding leftover terra Spat* object(s) from stage3.rds: ",
          paste(spat_leftover, collapse = ", "), call. = FALSE)
}

stage3_objects <- setdiff(
  ls(all.names = TRUE, envir = environment()),
  c("STAGE_DIR", "stage2_path", "swe", "scf_stack", "scf_stack_path", "basin",
    "spat_leftover")
)
stage3_objects <- setdiff(stage3_objects, spat_leftover)
stage3_path <- file.path(STAGE_DIR, "stage3.rds")
saveRDS(mget(stage3_objects, envir = environment()), stage3_path, compress = "xz")

# Verify the checkpoint immediately while any serialization error is still
# attributable to this stage.
stage3_verify <- readRDS(stage3_path)
if (!is.list(stage3_verify) || length(stage3_verify) != length(stage3_objects)) {
  stop("stage3.rds verification failed after saveRDS().", call. = FALSE)
}
rm(stage3_verify)

cat(sprintf("[OK] Saved %d object(s) to: %s\n", length(stage3_objects), normalizePath(stage3_path)))
cat(sprintf("[OK] File size: %.1f MB\n", file.size(stage3_path) / 1024^2))
cat("\nNext: run 04_plots_and_summary.R\n")
cat("(Once 04 has run successfully and you've confirmed the outputs in\n")
cat(" out_dir look right, it's safe to delete stage1.rds, stage2.rds, AND\n")
cat(" stage3.rds - all the real deliverables are the files already written\n")
cat(" to out_dir, not these checkpoints.)\n")