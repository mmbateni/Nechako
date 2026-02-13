# SWEI CALCULATION WITH TIMESCALE-SPECIFIC SCF MASKS AND GRINGORTEN PLOTTING POSITION
# Merged modifications: read SCF, compute k-month SCF climatology, apply per-timescale SCF threshold masks,
# replace zeros with tiny positive values before ranking, compute SWEI using Gringorten plotting position.
# Requirements: terra, parallel, ncdf4, zoo, writexl
# This script reads ERA5???Land SWE and SCF NetCDFs, reprojects and masks SWE
# to your basin, and ensures units and time vectors are consistent.
# It computes centered k???month SCF climatologies and builds a separate binary mask
# for each center month using the user???specified SCF threshold 
# for that timescale. For every grid cell and timescale the code 
# forms k???month aggregated SWE (right???aligned), replaces zero SWE values 
# with a tiny positive value derived from the smallest nonzero sample, 
# and then ranks the masked samples by calendar month. 
# It converts Gringorten plotting positions to standard normal variates
# to produce the SWEI time series, clips extreme z values, 
# and reports simple diagnostics such as clipping counts and basin NA rate.
# Processing is parallelized across cells and repeated for each timescale,
# and the script writes per???month CSVs, Excel workbooks, NetCDF SWEI files,
# and the SCF masks used. The result is a set of timescale???specific SWEI products
# and masks that reflect the different SCF thresholds you chose.
library(terra)
library(parallel)
library(ncdf4)
library(zoo)
library(writexl)

# ---------------------------
# User inputs / file paths
# ---------------------------
setwd("D:/Nechako_Drought/Nechako")
out_dir <- "swei_results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

basin_path <- "Spatial/nechakoBound_dissolve.shp"
basin <- if (file.exists(basin_path)) vect(basin_path) else NULL
if (is.null(basin)) stop("Basin boundary not found: ", basin_path)

swe_file <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
scf_file <- "monthly_data_direct/snow_cover_monthly.nc"
if (!file.exists(swe_file)) stop("SWE file not found: ", swe_file)
if (!file.exists(scf_file)) stop("SCF file not found: ", scf_file)

# ---------------------------
# Read SWE and SCF
# ---------------------------
cat("Reading SWE and SCF...\n")
swe <- rast(swe_file)
scf_stack <- rast(scf_file)

# Robust time extraction
extract_time_dimension <- function(raster_obj, file_path) {
  t <- time(raster_obj)
  if (!all(is.na(t))) return(as.Date(t))
  nc <- nc_open(file_path)
  on.exit(nc_close(nc))
  # try common names
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
scf_max <- global(scf_stack, "max", na.rm = TRUE)$max
if (length(scf_max) > 1) scf_max <- max(scf_max, na.rm = TRUE)
if (!is.na(scf_max) && scf_max > 1.5) scf_stack <- scf_stack / 100

# Reproject SWE to SCF grid if needed
if (!compareGeom(scf_stack[[1]], swe[[1]], stopOnError = FALSE)) {
  swe <- resample(swe, scf_stack, method = "bilinear")
}

# Reproject basin to SCF grid CRS if needed
if (!is.null(basin) && !same.crs(basin, crs(scf_stack))) basin <- project(basin, crs(scf_stack))

# Mask SWE to basin early
swe <- mask(swe, basin, inverse = FALSE, touches = TRUE)

# Convert SWE units if needed (meters -> mm)
swe_mean_sample <- global(swe[[1]], "mean", na.rm = TRUE)$mean
if (!is.na(swe_mean_sample) && swe_mean_sample < 10) swe <- swe * 1000

# Prepare SWE matrix (cells x time)
swe_matrix <- values(swe, mat = TRUE)
seq(start, by = "month", length.out = n)
}

dates_swe <- extract_time_dimension(swe, swe_file)
dates_scf <- extract_time_dimension(scf_stack, scf_file)
terra::time(swe) <- dates_swe
terra::time(scf_stack) <- dates_scf

# Convert SCF to fraction if in percent
scf_max <- global(scf_stack, "max", na.rm = TRUE)$max
if (length(scf_max) > 1) scf_max <- max(scf_max, na.rm = TRUE)
if (!is.na(scf_max) && scf_max > 1.5) scf_stack <- scf_stack / 100

# Reproject SWE to SCF grid if needed
if (!compareGeom(scf_stack[[1]], swe[[1]], stopOnError = FALSE)) {
  swe <- resample(swe, scf_stack, method = "bilinear")
}

# Reproject basin to SCF grid CRS if needed
if (!is.null(basin) && !same.crs(basin, crs(scf_stack))) basin <- project(basin, crs(scf_stack))

# Mask SWE to basin early
swe <- mask(swe, basin, inverse = FALSE, touches = TRUE)

# Convert SWE units if needed (meters -> mm)
swe_mean_sample <- global(swe[[1]], "mean", na.rm = TRUE)$mean
if (!is.na(swe_mean_sample) && swe_mean_sample < 10) swe <- swe * 1000

# Prepare SWE matrix (cells x time)
swe_matrix <- values(swe, mat = TRUE)
dates <- dates_swe
months_all <- as.integer(format(dates, "%m"))

# ---------------------------
# User-defined per-timescale SCF thresholds
# Edit these values as desired; keys are timescales in months
# ---------------------------
scf_thresholds <- c(
  "1"  = 0.05,
  "2"  = 0.05,
  "3"  = 0dates <- dates_swe
  months_all <- as.integer(format(dates, "%m"))
  
  # ---------------------------
  # User-defined per-timescale SCF thresholds
  # Edit these values as desired; keys are timescales in months
  # ---------------------------
  scf_thresholds <- c(
    "1"  = 0.05,
    "2"  = 0.05,
    "3"  = 0.075,
    "4"  = 0.08,
    "5"  = 0.09,
    "6"  = 0.10,
    "7"  = 0.11,
    "8"  = 0.12,
    "9"  = 0.13,
    "10" = 0.14,
    "11" = 0.145,
    "12" = 0.15
  )
  
  # ---------------------------
  # Helper functions
  # ---------------------------
  # Right.075,
  "4"  = 0.08,
  "5"  = 0.09,
  "6"  = 0.10,
  "7"  = 0.11,
  "8"  = 0.12,
  "9"  = 0.13,
  "10" = 0.14,
  "11" = 0.145,
  "12" = 0.15
)

# ---------------------------
# Helper functions
# ---------------------------
# Right-aligned rolling-aligned rolling sum (keeps original behavior)
roll_sum_right <- function(x, k) {
  if (length(x) < k) return(rep(NA_real_, length(x)))
  zoo::rollapply(x, width = k, FUN = sum, align = "right", fill = NA_real_, partial = FALSE)
}

# Centered rolling sum (optional; not used here but available)
centered_roll_sum <- function(x, k) {
  sum (keeps original behavior)
  roll_sum_right <- function(x, k) {
    if (length(x) < k) return(rep(NA_real_, length(x)))
    zoo::rollapply(x, width = k, FUN = sum, align = "right", fill = NA_real_, partial = FALSE)
  }
  
  # Centered rolling sum (optional; not used here but available)
  centered_roll_sum <- function(x, k) {
    if (k == 1) return  if (k == 1) return(x)
    zoo::rollapply(x, width = k, FUN = sum, align = "center", fill = NA_real_, partial = FALSE)
  }
  
  clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)
  
  # Compute k-month SCF climatology per center month (returns list of 12 rasters)
  compute_kmonth_scf_climatology <- function(scf_rast, dates, k) {
    (x)
    zoo::rollapply(x, width = k, FUN = sum, align = "center", fill = NA_real_, partial = FALSE)
  }
  
  clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)
  
  # Compute k-month SCF climatology per center month (returns list of 12 rasters)
  compute_kmonth_scf_climatology <- function(scf_rast, dates, k) {
    months <- as.integer(format(dates, "%m"))
    centers <- 1:12
    scf_clim <- vector("list", 12)
    for (m in centers) {
      # build centered window months for width k (wrap-around)
      half <- floor((k-1)/2)
      months_window <- ((m - half - 1):(m + (k - 1 - half) - 1)) %% 12 + 1
      idx <- which(months %in% months_window)
      if (length(idx) == 0) {
        scf_clim[[m]] <- rast(scf_rast[[1]])
        values(scf_clim[[m]]) <- NA
      } else {
        scf_clim[[m]] <- mean(scf_rast[[idx]], na.rm = TRUE)
      }
    }
    scf_clim
  }
  
  # ---------------------------
  # Precompute SCF climatologies for all timescales we will compute
  # ---------------------------
  timescales <- as.integer(names(scf_thresholds))
  scf_clims_all <- list()
  for (k in unique(timescales)) {
    scf_clims_all[[  months <- as.integer(format(dates, "%m"))
                     centers <- 1:12
                     scf_clim <- vector("list", 12)
                     for (m in centers) {
                       # build centered window months for width k (wrap-around)
                       half <- floor((k-1)/2)
                       months_window <- ((m - half - 1):(m + (k - 1 - half) - 1)) %% 12 + 1
                       idx <- which(months %in% months_window)
                       if (length(idx) == 0) {
                         scf_clim[[m]] <- rast(scf_rast[[1]])
                         values(scf_clim[[m]]) <- NA
                       } else {
                         scf_clim[[m]] <- mean(scf_rast[[idx]], na.rm = TRUE)
                       }
                     }
                     scf_clim
  }
  
  # ---------------------------
  # Precompute SCF climatologies for all timescales we will compute
  # ---------------------------
  timescales <- as.integer(names(scf_thresholds))
  scf_clims_all <- list()
  for (k in unique(timescales)) {
    scf_clims_all[[as.character(k)]] <- compute_kmonth_scf_climatology(scf_stack, dates_scf, k)
  }
  
  # ---------------------------
  # Modified Gringorten SWEIas.character(k)]] <- compute_kmonth_scf_climatology(scf_stack, dates_scf, k)
}

# ---------------------------
# Modified Gringorten SWEI function that applies per-center-month mask and zero-replacement
# Inputs:
#   v: numeric vector of SWE time series for a single cell
#   scale: integer timescale (months)
#   dates_vec: vector of dates
#   mask_months: list of 12 logical vectors (length = ncell) indicating which cells to keep for each center month
#   cell_index: integer index of the cell being processed
# ---------------------------
gringorten_swei_masked <- function(v, scale, dates_vec, mask_months, cell_index, eps = 1e-6) {
  if (length(v) == 0 || all(is.na(v)) || is.null(v)) return(rep(NA_real_, length(dates_vec)))
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  # aggregate (right-aligned to preserve original behavior)
  x_agg <- if (scale == 1) v_clean else roll_sum_right(v_clean, scale)
  if (length(x_agg) != length(dates_vec)) {
    x_agg function that applies per-center-month mask and zero-replacement
    # Inputs:
    #   v: numeric vector of SWE time series for a single cell
    #   scale: integer timescale (months)
    #   dates_vec: vector of dates
    #   mask_months: list of 12 logical vectors (length = ncell) indicating which cells to keep for each center month
    #   cell_index: integer index of the cell being processed
    # ---------------------------
    gringorten_swei_masked <- function(v, scale, dates_vec, mask_months, cell_index, eps = 1e-6) {
      if (length(v) == 0 || all(is.na(v)) || is.null(v)) return(rep(NA_real_, length(dates_vec)))
      v_clean <- v
      v_clean[!is.finite(v_clean)] <- NA_real_
      # aggregate (right-aligned to preserve original behavior)
      x_agg <- if (scale == 1) v_clean else roll_sum_right(v_clean, scale)
      if (length(x_agg) != length(dates_vec)) {
        x_agg <- c(rep(NA_real_, length(dates_vec) - length(x_agg)), x_agg[1:length(dates_vec)])
        <- c(rep(NA_real_, length(dates_vec) - length(x_agg)), x_agg[1:length(dates_vec)])
      }
      n <- length(x_agg)
      z <- rep(NA_real_, n)
      mon <- as.integer(format(dates_vec, "%m"))
      # For each calendar month, only include indices where mask_months[[m]][cell_index] is TRUE
      for (m in 1:12) {
        # indices for this calendar month
        idx_all <- which(mon == m & is.finite(x_agg))
        if (length(idx_all) < 3) next
        # apply mask for this cell and month
        keep_mask <- TRUE
        if (!is.null(mask_months) && length(mask  }
      n <- length(x_agg)
      z <- rep(NA_real_, n)
      mon <- as.integer(format(dates_vec, "%m"))
      # For each calendar month, only include indices where mask_months[[m]][cell_index] is TRUE
      for (m in 1:12) {
        # indices for this calendar month
        idx_all <- which(mon == m & is.finite(x_agg))
        if (length(idx_all) < 3) next
        # apply mask for this cell and month
        keep_mask <- TRUE
        if (!is.null(mask_months) && length(mask_months) >= m) {
          mask_vec <- mask_months[[m]]
          if (length(mask_vec) >= cell_index) keep_mask_months) >= m) {
          mask_vec <- mask_months[[m]]
          if (length(mask_vec) >= cell_index) keep_mask <- isTRUE(mask_vec[cell_index]) else keep_mask <- FALSE
        }
        if (!keep_mask) next
        idx <- idx_all
        samp <- x_agg[idx]
        # Replace zeros with tiny positive value per SSRN: 0.001% of smallest non-zero
        if (any(!is.na(samp) & samp == 0)) {
          nonzero_vals <- samp[!is.na(samp) & samp > 0]
          if (length(nonzero_vals) > 0) {
            tiny <- min(nonzero_vals, na.rm = TRUE) * 1e-5  # 0.001% = 1e-5
            samp[samp <- isTRUE(mask_vec[cell_index]) else keep_mask <- FALSE
          }
          if (!keep_mask) next
          idx <- idx_all
          samp <- x_agg[idx]
          # Replace zeros with tiny positive value per SSRN: 0.001% of smallest non-zero
          if (any(!is.na(samp) & samp == 0)) {
            nonzero_vals <- samp[!is.na(samp) & samp > 0]
            if (length(nonzero_vals) > 0) {
              tiny <- min(nonzero_vals, na.rm = TRUE) * 1e-5  # 0.001% = 1e-5
              samp[samp == 0] <- tiny
            } else {
              # if no nonzero values, leave zeros == 0] <- tiny
            } else {
              # if no nonzero values, leave zeros as-is (will rank ties)
              samp[samp == 0] <- 0
            }
          }
          # compute ranks and Gringorten p
          N <- length(samp)
          ranks <- rank(samp, ties.method = "average")
          p <- (ranks - 0.44) / (N + 0.12)
          p <- clip_prob(p, eps = eps)
          z[idx] <- qnorm(p)
        }
        # clip extremes
        finite_idx <- is.finite(z)
        z[finite_idx & z < -4.75] <- -4.75
        z[finite_idx & as-is (will rank ties)
          samp[samp == 0] <- 0
        }
      }
      # compute ranks and Gringorten p
      N <- length(samp)
      ranks <- rank(samp, ties.method = "average")
      p <- (ranks - 0.44) / (N + 0.12)
      p <- clip_prob(p, eps = eps)
      z[idx] <- qnorm(p)
    }
    # clip extremes
    finite_idx <- is.finite(z)
    z[finite_idx & z < -4.75] <- -4.75
    z[finite_idx & z >  4.75] <-  4.75
    z
  }
  
  # ---------------------------
  # Main z >  4.75] <-  4.75
  z
}

# ---------------------------
# Main SWEI calculation loop with per-timescale SCF masks applied
# ---------------------------
n_cores <- max(1, detectCores() - 1)
cat(sprintf("Using %d CPU cores\n", n_cores))

month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
eps_clip <- 1e-6

for (sc in sort(unique(timescales))) {
  sc_str <- as.character(sc)
  cat(sprintf("\n=== SWEI-%d (applying SCF threshold = %.3f) ===\n", sc, scf_threshold SWEI calculation loop with per-timescale SCF masks applied
              # ---------------------------
              n_cores <- max(1, detectCores() - 1)
              cat(sprintf("Using %d CPU cores\n", n_cores))
              
              month_names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
              eps_clip <- 1e-6
              
              for (sc in sort(unique(timescales))) {
                sc_str <- as.character(sc)
                cat(sprintf("\n=== SWEI-%d (applying SCF threshold = %.3f) ===\n", sc, scf_thresholds[sc_str]))
                # build mask_months for this timescale: list of 12 logical vectors (length = ncell)
                scf_clim_list <- scf_clims_all[[sc_str]]
                threshold <- scf_thresholds[sc_str]
                mask_months <- vector("list", 12)
                for (m in 1:12) {
                  r <- scf_clim_list[[m]]
                  if (is.null(r)) {
                    mask_months[[m]] <- rep(FALSE, ncell(swe))
                  } else {
                    vals <- values(r)
                    mask_months[[m]] <- !is.na(vals) & (vals >= threshold)
                    s[sc_str]))
  # build mask_months for this timescale: list of 12 logical vectors (length = ncell)
  scf_clim_list <- scf_clims_all[[sc_str]]
  threshold <- scf_thresholds[sc_str]
  mask_months <- vector("list", 12)
  for (m in 1:12) {
    r <- scf_clim_list[[m]]
    if (is.null(r)) {
      mask_months[[m]] <- rep(FALSE, ncell(swe))
    } else {
      vals <- values(r)
      mask_months[[m]] <- !is.na(vals) & (vals >= threshold)
    }
  }
  # Export }
                  }
                  # Export variables to cluster
                  cl <- makeCluster(n_cores)
                  clusterExport(cl, c("swe_matrix", "dates", "gringorten_swei_masked", "roll_sum_right", "clip_prob", "sc", "eps_clip", "mask_months"), envir = environment())
                  clusterEvalQ(cl, { library(zoo) })
                  # compute SWEI per cell (parallel)
                  swei_results <- parLapply(cl, 1:nrow(swe_matrix), function(i) {
                    tryCatch({
                      gringorten_swei_masked(swe_matrix[i, ], sc, dates, mask_months, i variables to cluster
                                             cl <- makeCluster(n_cores)
                                             clusterExport(cl, c("swe_matrix", "dates", "gringorten_swei_masked", "roll_sum_right", "clip_prob", ", eps = eps_clip)
    }, error = function(e) {
      rep(NA_real_, length(dates))
    })
  })
  stopCluster(cl)
  swei_indices <- do.call(rbind, swei_results)
  if (is.null(dim(swei_indices)) || nrow(swei_indices) != nrow(swe_matrix)) stop(sprintf("SWEI-%d calculation failed", sc))
  # Clip and re-normalize probabilities to avoid numerical issues
  z <- swei_indices
  z[!is.finite(z)] <- NA_real_
  p <- pnorm(z)
  p_clip <- clip_prob(p, eps = eps_clip)
  z_clip <- qnorm(p_clip)
  z_clip[is.na(z)] <- NA_real_
  swei_indices <- z_clip
  # Diagnostics
  clipped_low  <- sum(z < -4.7, na.rm = TRUE)
  clipped_high <- sum(z >  4.7, na.rm = TRUE)
  cat(sprintf("  Clipping: low=%d high=%d ratio=%.2f\n", clipped_low, clipped_high, clipped_low / max(1, clipped_high)))
  basin_mask_swei <- !is.na(swe_matrix[, 1])
  na_rate_basin <- 100 * mean(is.na(swei_indices[basin_mask_swei, ]))
  cat(sprintf("  NA rate (basin pixels): %.3f%%\n", na_rate_basin))
  # Outputs: CSV per month, Excel workbook, NetCDF per month (same as original)
  months_all <- as.numeric(format(dates, "%m"))
  # CSV files
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    csv_file <- file.path(out_dir, sprintf("swei_%02d_k%02d_month%02d_%s.csv", sc, sc, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE)
  }
  # Excel workbook
  excel_data <- list()
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
d high=%d ratio=%.2f\n", clipped_low, clipped_high, clipped_low / max(1, clipped_high)))
                                             basin_mask_swei <- !is.na(swe_matrix[, 1])
                                             na_rate_basin <- 100 * mean(is.na(swei_indices[basin_mask_swei, ]))
                                             cat(sprintf("  NA rate (basin pixels): %.3f%%\n", na_rate_basin))
                                             # Outputs: CSV per month, Excel workbook, NetCDF per month (same as original)
                                             months_all <- as.numeric(format(dates, "%m"))
                                             # CSV files
                                             for (m in 1:12) {
                                               idx_m <- which(months_all == m)
                                               if (length(idx_m) == 0) next
                                               swei_subset <- swei_indices[, idx_m, drop = FALSE]
                                               dates_subset <- dates[idx_m]
                                               df <- as.data.frame(swei_subset)
                                               colnames(df) <- format(dates_subset, "%Y")
                                               coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
                                               df <- cbind(lon = coords[,1], lat = coords[,2], df)
                                               csv_file <- file.path(out_dir, sprintf("swei_%02d_k%02d_month%02d_%s.csv", sc, sc, m, month_names[m]))
                                               write.csv(df, csv_file, row.names = FALSE)
                                             }
                                             # Excel workbook
                                             excel_data <- list()
                                             for (m in 1:12) {
                                               idx_m <- which(months_all == m)
                                               if (length(idx_m) == 0) next
                                               swei_subset <- swei_indices[, idx_m, drop = FALSE]
                                               dates_subset <- dates[idx_m]
                                               df <- as.data.frame(swei_subset)
                                               colnames(df) <- format(dates_subset, "%Y")
                                               coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
                                               df <- cbind(lon = coords[,1], lat = coords[,2], df)
                                               excel_data[[month_names[m]]] <- df
                                             }
                                             xlsx_file <- file.path(out_dir, sprintf("swei_k%02d_all_months.xlsx", sc))
                                             if (file.exists(xlsx_file)) file.remove(xlsx_file)
                                             write_xlsx(ex    excel_data[[month_names[m]]] <- df
                    }
                    xlsx_file <- file.path(out_dir, sprintf("swei_k%02d_all_months.xlsx", sc))
                    if (file.exists(xlsx_file)) file.remove(xlsx_file)
                    write_xlsx(excel_data, xlsx_file)
                    # Save NetCDFcel_data, xlsx_file)
                    # Save NetCDF per month
                    for (m in 1:12) {
                      idx_m <- which(months_all == m)
                      if (length(idx_m) == 0) next
                      swei_rast <- rast(swe[[1]])
                      swei_rast <- rep(swei_rast, length(idx_m))
                      values(swei_rast) <- swei_indices[, idx_m, drop = FALSE]
                      terra::time(swei_rast) <- dates[idx_m]
                      nc_file <- file.path(out_dir, sprintf("swei_k%02d_month%02d_%s.nc", sc, m, month_names[m]))
                      writeCDF(swei_rast, nc_file,
                               varname = "swei",
                               longname = sprintf("SWEI k=%d months (SCF threshold=%.3f)", sc, threshold),
                               unit = "standardized_index",
                               missval = -9999,
                               overwrite = TRUE)
                    }
                    # Save mask used for this timescale (12 layers)
                    mask_stack <- rast(scf_stack[[1]])
                    masks <- lapply(mask_months, function(m) {
                      r <- rast(scf_stack[[1]])
                      per month
                      for (m in 1:12) {
                        idx_m <- which(months_all == m)
                        if (length(idx_m) == 0) next
                        swei_rast <- rast(swe[[1]])
                        swei_rast <- rep(swei_rast, length(idx_m))
                        values(swei_rast) <- swei_indices[, idx_m, drop = FALSE]
                        terra::time(swei_rast) <- dates[idx_m]
                        nc_file <- file.path(out_dir, sprintf("swei_k%02d_month%02d_%s.nc", sc, m, month_names[m]))
                        writeCDF(swei_rast, nc_file,
                                 varname = "swei",
                                 longname = sprintf("SWEI k=%d months (SCF threshold=%.3f)", sc, threshold),
                                 unit = "standardized_index",
                                 missval = -9999,
                                 overwrite = TRUE)
                      }
                      # Save mask used for this timescale (12 layers)
                      mask_stack <- rast(scf_stack[[1]])
                      masks <- lapply(mask_months, function(m) {
                        r <- rast(scf_stack[[1]])
                        vals <- rep(NA_real_, ncell(r))
                        vals <- rep(NA_real_, ncell(r))
                        vals[which(!is.na(values(scf_stack[[1]])))] <- as.numeric(m)
                        values(r) <- vals
                        r
                      })
                      mask_stack <- rast(masks)
                      names(mask_stack) <- paste0("mask_k", sc, "_m", sprintf("%02d", 1:12))
                      writeCDF(mask_stack, file.path(out_dir, sprintf("scf_mask_k%02d_threshold_%.3f.nc", sc, threshold)),
                               varname = "scf_mask", longname = sprintf("SCF mask k=%d threshold %.3f", sc, threshold), unit = "logical", overwrite = TRUE)
                      gc()
                    }
                    
                    cat("\nSWEI calculation    vals[which(!is.na(values(scf_stack[[1]])))] <- as.numeric(m)
    values(r) <- vals
    r
  })
  mask_stack <- rast(masks)
  names(mask_stack) <- paste0("mask_k", sc, "_m", sprintf("%02d", 1:12))
  writeCDF(mask_stack, file.path(out_dir, sprintf("scf_mask_k%02d_threshold_%.3f.nc", sc, threshold)),
           varname = "scf_mask", longname = sprintf("SCF mask k=%d threshold %.3f", sc, threshold), unit = "logical", overwrite = TRUE)
  gc()
}

cat("\nSWEI calculation with timescale-specific SCF masks complete with timescale-specific SCF masks complete. Outputs in:", normalizePath(out_dir), "\n")
```. Outputs in:", normalizePath(out_dir), "\n")
                    