####################################################################################
# PIXEL-WISE TEMPORAL DIAGNOSTICS FOR DROUGHT INDICES (SPI/SPEI) - OPTIMIZED v3.2+
# Nechako Basin Analysis - 76 Years (1950-2025)
# MODIFIED TO PRODUCE BOTH TIFS (GIS) AND PDFS (WINDOWS-FRIENDLY)
# NOTE (UPDATED):
#   • No PDFs are written to `maps_with_legends` (PNGs only)
#   • Every TIF saved in `temporal_spi_spei` gets a PDF sibling with matching basename
#   • PNG map titles are drawn in outer margin (no overlap)
#   • Legend label clipping removed in PNGs
#   • (NEW) PDF significance maps use a dedicated legend panel (no trimming)
#   • (NEW) Drought frequency PDF maps have no bottom footer (no date/descriptions)
####################################################################################

library(terra)
library(data.table)
library(Kendall)
library(changepoint)
library(extRemes)
library(future.apply)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(zoo)
library(tseries)
library(viridis)
library(scales)
library(stringr)
library(gstat) # IDW
library(sf)    # sf objects for IDW data & grid

# ------------------------------------------------------------------------------
# Output directory
# ------------------------------------------------------------------------------
setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "temporal_spi_spei"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Diagnostic subdirectories
maps_dir <- file.path(out_dir, "maps_with_legends")
for (d in c(maps_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Log file setup
# ------------------------------------------------------------------------------
LOG_FILE <- file.path(out_dir, "drought_temporal_diagnostics_optimized.log")
cat("OPTIMIZED DROUGHT TEMPORAL DIAGNOSTICS (v3.2+)\n", file = LOG_FILE)
cat("======================================================================\n", file = LOG_FILE, append = TRUE)
cat(paste("Analysis started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat("======================================================================\n", file = LOG_FILE, append = TRUE)

log_event <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_msg <- sprintf("[%s] [%s] %s", timestamp, level, msg)
  cat(log_msg, "\n", file = LOG_FILE, append = TRUE)
  if (level %in% c("INFO", "SUCCESS")) {
    message(sprintf("[%s] %s", timestamp, msg))
  } else if (level == "WARNING") {
    warning(sprintf("[%s] %s", timestamp, msg), call. = FALSE)
  } else if (level == "ERROR") {
    stop(sprintf("[%s] %s", timestamp, msg))
  }
}
log_event("Starting drought temporal diagnostics (v3.2+)...", "INFO")

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
alpha <- 0.05
drought_threshold <- -0.52
recovery_threshold <- 0
extreme_threshold <- -2.0
min_duration <- 2
max_tie_percent <- 50
min_valid_obs <- 60 # Critical pre-filtering threshold
n_years_total <- 76

# ------------------------------------------------------------------------------
# Parallel setup
# ------------------------------------------------------------------------------
num_cores <- min(parallel::detectCores(logical = FALSE) - 1, 6)
if (is.na(num_cores) || num_cores < 1) num_cores <- 1
options(future.globals.maxSize = 1500 * 1024^2)
plan(multisession, workers = num_cores)
log_event(sprintf("Base parallel configuration: %d cores", num_cores), "INFO")

# ------------------------------------------------------------------------------
# Basin boundary load/reproject
# ------------------------------------------------------------------------------
log_event("Loading Nechako Basin boundary...", "INFO")
basin_path <- "Spatial/nechakoBound_dissolve.shp"
if (!file.exists(basin_path)) {
  log_event(sprintf("ERROR: Basin boundary file not found: %s", basin_path), "ERROR")
}
basin <- vect(basin_path)
target_crs <- "EPSG:3005"
if (!same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
  log_event("Basin boundary reprojected to BC Albers (EPSG:3005)", "SUCCESS")
} else {
  log_event("Basin boundary already in BC Albers projection", "SUCCESS")
}

# ==============================================================================
# CORE STATISTICAL FUNCTIONS
# ==============================================================================
calculate_sens_slope_manual <- function(x) {
  n <- length(x)
  if (n < 2) return(NA)
  d_mat <- outer(x, x, "-")
  t_mat <- outer(1:n, 1:n, "-")
  upper.tri.inds <- which(upper.tri(d_mat))
  slopes <- d_mat[upper.tri.inds] / t_mat[upper.tri.inds]
  return(median(slopes, na.rm = TRUE))
}

calculate_variance_with_ties <- function(S, n, x) {
  tie_table <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  if (length(tie_counts) == 0) {
    var_s <- n * (n - 1) * (2 * n + 5) / 18
  } else {
    var_s <- n * (n - 1) * (2 * n + 5) / 18
    tie_adjustment <- sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
    var_s <- var_s - tie_adjustment
  }
  return(var_s)
}

# ============================================================
# Mann–Kendall (variance-corrected) with detrended ρ1
# Matches codeR4 logic: Sen–detrend -> ρ1 on residuals -> VC
# ============================================================
modified_mann_kendall_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50, apply_vc = TRUE) {
  
  n_space <- ncol(ts_matrix)
  
  results_list <- future.apply::future_lapply(seq_len(n_space), function(i) {
    
    ts_clean <- stats::na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    if (n < 10) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  vc_corrected=FALSE, n_ties=0, percent_ties=0,
                  tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    # Ties screen
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    if (percent_ties > max_tie_pct) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  vc_corrected=FALSE, n_ties=n_ties, percent_ties=percent_ties,
                  tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    # Sen slope on original series (for reporting AND detrending)
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  vc_corrected=FALSE, n_ties=n_ties, percent_ties=percent_ties,
                  tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    # Detrend using Sen slope, then estimate rho1 on residuals (not raw)
    t_idx <- seq_len(n)
    detrended <- ts_clean - (sen_slope * t_idx)
    acf_result <- tryCatch(stats::acf(detrended, lag.max = 1, plot = FALSE, na.action = stats::na.pass),
                           error = function(e) NULL, warning = function(w) NULL)
    rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA_real_
    
    # MK core on ORIGINAL (not detrended) series
    diff_mat <- outer(ts_clean, ts_clean, "-")
    S <- sum(sign(diff_mat[upper.tri(diff_mat)]), na.rm = TRUE)
    varS_taub <- calculate_variance_with_ties(S, n, ts_clean)
    n_pairs <- n * (n - 1) / 2
    tau_val <- S / n_pairs
    tau_b_adjusted <- (percent_ties > 5)
    
    # Variance correction if requested and |rho1|>0.1
    vc_corrected <- FALSE
    varS_final <- varS_taub
    if (apply_vc && !is.na(rho1) && abs(rho1) > 0.1) {
      # same correction structure as codeR4
      correction_factor <- 1 + (2 * rho1 * (n - 1 - 2 * (n - 1) * rho1 + 3 * rho1 * rho1)) /
        ((n - 1) * (1 - rho1) * (1 - rho1))
      varS_final <- varS_taub * correction_factor
      vc_corrected <- TRUE
    }
    
    # p-value
    p_value <- if (varS_final <= 0) NA_real_ else {
      z_stat <- S / sqrt(varS_final)
      2 * stats::pnorm(-abs(z_stat))
    }
    
    list(
      tau = tau_val, p.value = p_value, sl = sen_slope, S = S, varS = varS_final,
      n = n, rho1 = rho1, vc_corrected = vc_corrected, n_ties = n_ties,
      percent_ties = percent_ties, tau_b_adjusted = tau_b_adjusted, filtered = FALSE
    )
  }, future.seed = TRUE)
  
  data.table::rbindlist(results_list)
}
# ============================================================
# FULL TFPW (Yue–Pilon style) MK test
# Sen–detrend -> estimate rho1 on residuals -> AR(1) pre-whiten ->
# add trend back -> MK on corrected series
# ============================================================
perform_tfpw_mk_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50) {
  
  n_space <- ncol(ts_matrix)
  
  results_list <- future.apply::future_lapply(seq_len(n_space), function(i) {
    
    ts_clean <- stats::na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    
    if (n < 10) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  tfpw_applied=FALSE, n_ties=0, percent_ties=0,
                  tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    # Ties screen
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    if (percent_ties > max_tie_pct) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  tfpw_applied=FALSE, n_ties=n_ties, percent_ties=percent_ties,
                  tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    # Sen slope on original series
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  tfpw_applied=FALSE, n_ties=n_ties, percent_ties=percent_ties,
                  tau_b_adjusted=FALSE, filtered=TRUE))
    }
    
    # Detrend -> rho1 on detrended
    t_idx <- seq_len(n)
    detrended <- ts_clean - (sen_slope * t_idx)
    acf_result <- tryCatch(stats::acf(detrended, lag.max = 1, plot = FALSE, na.action = stats::na.pass),
                           error = function(e) NULL, warning = function(w) NULL)
    rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA_real_
    
    # If |rho1|>0.1, pre-whiten; otherwise use original
    tfpw_applied <- !is.na(rho1) && abs(rho1) > 0.1
    if (tfpw_applied) {
      # AR(1) pre-whitening on detrended series
      pw <- detrended
      if (n >= 2) {
        pw[2:n] <- detrended[2:n] - rho1 * detrended[1:(n - 1)]
      }
      # Reintroduce trend
      corrected <- pw + (sen_slope * t_idx)
    } else {
      corrected <- ts_clean
    }
    
    # MK on corrected series
    diff_mat <- outer(corrected, corrected, "-")
    S <- sum(sign(diff_mat[upper.tri(diff_mat)]), na.rm = TRUE)
    varS_taub <- calculate_variance_with_ties(S, n, corrected)
    n_pairs <- n * (n - 1) / 2
    tau_val <- S / n_pairs
    tau_b_adjusted <- (percent_ties > 5)
    
    p_value <- if (varS_taub <= 0) NA_real_ else {
      z_stat <- S / sqrt(varS_taub)
      2 * stats::pnorm(-abs(z_stat))
    }
    
    list(
      tau = tau_val, p.value = p_value, sl = sen_slope, S = S, varS = varS_taub,
      n = n, rho1 = rho1, tfpw_applied = tfpw_applied, n_ties = n_ties,
      percent_ties = percent_ties, tau_b_adjusted = tau_b_adjusted, filtered = FALSE
    )
  }, future.seed = TRUE)
  
  data.table::rbindlist(results_list)
}

wald_wolfowitz_runs_test <- function(ts_vec, threshold) {
  ts_clean <- na.omit(ts_vec)
  n <- length(ts_clean)
  if (n < 5) {
    return(list(n_runs=NA, n_drought=NA, n_non_drought=NA, expected_runs=NA,
                var_runs=NA, z_stat=NA, p_value=NA, clustering=NA, filtered=TRUE))
  }
  binary <- ifelse(ts_clean <= threshold, 1, 0)
  n1 <- sum(binary == 1)
  n2 <- sum(binary == 0)
  if (n1 == 0 || n2 == 0) {
    return(list(n_runs=NA, n_drought=n1, n_non_drought=n2, expected_runs=NA,
                var_runs=NA, z_stat=NA, p_value=NA, clustering=NA, filtered=TRUE))
  }
  runs <- 1
  for (i in 2:n) if (binary[i] != binary[i - 1]) runs <- runs + 1
  expected <- (2 * n1 * n2) / n + 1
  variance <- (2 * n1 * n2 * (2 * n1 * n2 - n)) / (n^2 * (n - 1))
  if (variance <= 0) {
    return(list(n_runs=runs, n_drought=n1, n_non_drought=n2, expected_runs=expected,
                var_runs=variance, z_stat=NA, p_value=NA, clustering=NA, filtered=TRUE))
  }
  z_stat <- (runs - expected) / sqrt(variance)
  p_value <- 2 * pnorm(-abs(z_stat))
  clustering <- ifelse(runs < expected, "clustering", "dispersion")
  return(list(n_runs=runs, n_drought=n1, n_non_drought=n2, expected_runs=expected,
              var_runs=variance, z_stat=z_stat, p_value=p_value, clustering=clustering, filtered=FALSE))
}

detect_drought_events_pixel <- function(ts_vec, drought_thresh, recovery_thresh, min_dur) {
  ts_clean <- na.omit(ts_vec)
  n <- length(ts_clean)
  if (n < min_dur) {
    return(list(n_events=0, mean_duration=NA, mean_severity=NA, mean_intensity=NA,
                max_duration=NA, max_severity=NA, max_intensity=NA, total_severity=0))
  }
  in_drought <- ts_clean <= drought_thresh
  events <- list()
  start_idx <- NULL
  for (i in 1:n) {
    if (in_drought[i] && is.null(start_idx)) {
      start_idx <- i
    } else if (!in_drought[i] && !is.null(start_idx)) {
      end_idx <- i - 1
      duration <- end_idx - start_idx + 1
      if (duration >= min_dur) {
        severity <- sum(ts_clean[start_idx:end_idx])
        intensity <- severity / duration
        events[[length(events) + 1]] <- list(duration=duration, severity=severity, intensity=intensity)
      }
      start_idx <- NULL
    }
  }
  if (!is.null(start_idx)) {
    duration <- n - start_idx + 1
    if (duration >= min_dur) {
      severity <- sum(ts_clean[start_idx:n])
      intensity <- severity / duration
      events[[length(events) + 1]] <- list(duration=duration, severity=severity, intensity=intensity)
    }
  }
  if (length(events) == 0) {
    return(list(n_events=0, mean_duration=NA, mean_severity=NA, mean_intensity=NA,
                max_duration=NA, max_severity=NA, max_intensity=NA, total_severity=0))
  }
  durations  <- sapply(events, function(e) e$duration)
  severities <- sapply(events, function(e) e$severity)
  intensities<- sapply(events, function(e) e$intensity)
  return(list(
    n_events = length(events),
    mean_duration = mean(durations),
    mean_severity = mean(severities),
    mean_intensity = mean(intensities),
    max_duration = max(durations),
    max_severity = max(severities),
    max_intensity = max(intensities),
    total_severity = sum(severities)
  ))
}

detect_regime_shift <- function(ts_vec) {
  ts_clean <- na.omit(ts_vec)
  if (length(ts_clean) < 20) return(NA)
  result <- tryCatch({
    cpt <- cpt.mean(ts_clean, method = "PELT", penalty = "BIC")
    if (length(cpts(cpt)) > 0) cpts(cpt)[1] else NA
  }, error = function(e) NA)
  return(result)
}

calculate_return_period <- function(ts_vec, extreme_thresh = -2, conf = 0.90) {
  ts_clean <- na.omit(ts_vec)
  if (length(ts_clean) < 10) return(c(rl = NA_real_, lower = NA_real_, upper = NA_real_))
  
  extremes <- ts_clean[ts_clean <= extreme_thresh]
  if (length(extremes) < 5) return(c(rl = NA_real_, lower = NA_real_, upper = NA_real_))
  
  y <- -extremes
  
  out <- tryCatch({
    fit <- extRemes::fevd(y, type = "GEV")
    rl  <- extRemes::return.level(fit, return.period = 100, do.ci = TRUE, conf = conf)
    c(rl = -as.numeric(rl["Estimate"]), 
      lower = -as.numeric(rl["Lower CI"]), 
      upper = -as.numeric(rl["Upper CI"]))
  }, error = function(e) c(rl = NA_real_, lower = NA_real_, upper = NA_real_))
  
  out
}


read_monthly_csv_and_reconstruct <- function(index_type, scale, subfolder) {
  search_path <- file.path(getwd(), subfolder)
  if (!dir.exists(search_path)) {
    log_event(sprintf("WARNING: Directory not found: %s", search_path), "WARNING")
    return(NULL)
  }
  pattern <- sprintf("%s_0?%d_month\\d+.*\\.csv$", index_type, scale)
  monthly_files <- list.files(
    path = search_path, pattern = pattern, full.names = TRUE, ignore.case = TRUE
  )
  if (length(monthly_files) < 12) {
    log_event(sprintf("DEBUG: Only found %d files in %s matching pattern %s",
                      length(monthly_files), subfolder, pattern), "WARNING")
    return(NULL)
  }
  month_nums <- as.integer(gsub(".*month(\\d+).*", "\\1", basename(monthly_files), ignore.case = TRUE))
  ord <- order(month_nums)
  monthly_files <- monthly_files[ord]
  
  log_event(sprintf("Reconstructing %s-%02d from %s", toupper(index_type), scale, subfolder), "INFO")
  
  # Read first file to determine structure
  df_template <- fread(monthly_files[1], header = FALSE, na.strings = c("", "NA", "nan", "NaN", "."), fill = TRUE)
  n_years <- ncol(df_template) - 2
  n_pixels <- nrow(df_template)
  
  # CRITICAL FIX: Explicitly extract and convert coordinates to numeric
  lon_col <- as.numeric(as.character(df_template[[1]]))
  lat_col <- as.numeric(as.character(df_template[[2]]))
  
  # Validate coordinate conversion
  if (sum(is.na(lon_col)) > 0.1 * length(lon_col) || sum(is.na(lat_col)) > 0.1 * length(lat_col)) {
    log_event(sprintf("ERROR: >10%% of coordinates failed numeric conversion in %s", basename(monthly_files[1])), "ERROR")
    return(NULL)
  }
  
  coords <- data.table(lon = lon_col, lat = lat_col)
  
  # Build full time series matrix
  full_matrix <- matrix(NA, nrow = n_years * 12, ncol = n_pixels)
  for (m in 1:12) {
    df_month <- fread(monthly_files[m], header = FALSE, na.strings = c("", "NA", "nan", "NaN", "."), fill = TRUE)
    # Defensive check: ensure consistent dimensions
    if ((ncol(df_month) - 2) != n_years || nrow(df_month) != n_pixels) {
      log_event(sprintf("WARNING: Dimension mismatch in month %d file", m), "WARNING")
      next
    }
    year_values <- as.matrix(df_month[, 3:ncol(df_month), with = FALSE])
    for (yr in 1:n_years) {
      full_matrix[(yr - 1) * 12 + m, ] <- year_values[, yr]
    }
  }
  
  list(index_matrix = full_matrix, coords = coords, n_years = n_years, n_pixels = n_pixels)
}

# ==============================================================================
# CRITICAL: Filter to pixels with sufficient non-NA observations
# ==============================================================================
filter_basin_pixels <- function(index_matrix, coords, min_valid_obs = 60) {
  valid_obs <- colSums(!is.na(index_matrix))
  basin_pixel_idx <- which(valid_obs >= min_valid_obs)
  if (length(basin_pixel_idx) == 0) {
    return(list(index_matrix = matrix(NA, nrow = nrow(index_matrix), ncol = 0),
                coords = coords[0, ],
                n_pixels = 0,
                reduction_pct = 100))
  }
  index_matrix_basin <- index_matrix[, basin_pixel_idx, drop = FALSE]
  coords_basin <- coords[basin_pixel_idx, , drop = FALSE]
  n_pixels_orig <- ncol(index_matrix)
  n_pixels_basin <- length(basin_pixel_idx)
  reduction_pct <- 100 * (1 - n_pixels_basin / n_pixels_orig)
  list(index_matrix = index_matrix_basin,
       coords = coords_basin,
       n_pixels = n_pixels_basin,
       reduction_pct = reduction_pct)
}

# ==============================================================================
# OPTIMIZED INTERPOLATION WITH CRS HANDLING & BASIN MASKING (PATCHED)
# ==============================================================================
create_interpolated_raster_optimized <- function(coords_dt, value_col, raster_template,
                                                 basin_boundary, source_crs = NULL) {
  # Strong filtering: require non-NA value and non-NA lon/lat
  valid_data <- coords_dt[!is.na(get(value_col)) & !is.na(lon) & !is.na(lat)]
  if (nrow(valid_data) == 0) {
    log_event(sprintf("No valid data (or coords) for '%s' -> returning empty raster", value_col), "WARNING")
    result_rast <- rast(raster_template); values(result_rast) <- NA; return(result_rast)
  }
  
  coords_clean <- data.frame(
    x = as.numeric(valid_data$lon),
    y = as.numeric(valid_data$lat),
    value = as.numeric(valid_data[[value_col]])
  )
  coords_clean <- coords_clean[!is.na(coords_clean$x) & !is.na(coords_clean$y) & !is.na(coords_clean$value), ]
  if (nrow(coords_clean) < 3) {
    log_event(sprintf("Insufficient points (%d) for interpolation -> empty raster", nrow(coords_clean)), "WARNING")
    result_rast <- rast(raster_template); values(result_rast) <- NA; return(result_rast)
  }
  
  # Auto-detect source CRS if not supplied: degrees vs meters
  if (is.null(source_crs)) {
    looks_like_deg <- (max(abs(coords_clean$x)) <= 180) && (max(abs(coords_clean$y)) <= 90)
    source_crs <- if (looks_like_deg) "EPSG:4326" else "EPSG:3005"
  }
  
  # Build SpatVector and project to template CRS
  pts <- vect(coords_clean, geom = c("x","y"), crs = source_crs)
  if (!terra::same.crs(pts, raster_template)) {
    pts <- project(pts, crs(raster_template))
  }
  
  # Guard: zero points after projection (CRS mismatch / invalid coords)
  if (nrow(pts) == 0) {
    log_event("No valid geometries after projection — CRS mismatch suspected -> empty raster", "WARNING")
    result_rast <- rast(raster_template); values(result_rast) <- NA; return(result_rast)
  }
  
  # Try IDW using sf objects for both data and grid (gstat-friendly)
  result_rast <- tryCatch({
    # Points as sf
    pts_vals <- data.frame(value = values(pts)[,1])
    pts_xy <- as.data.frame(crds(pts))
    names(pts_xy) <- c("x","y")
    pts_sf <- sf::st_as_sf(
      cbind(pts_vals, pts_xy),
      coords = c("x","y"),
      crs = sf::st_crs(terra::crs(pts, proj = TRUE))
    )
    # Grid centers as sf points
    grid_pts <- terra::as.points(raster_template)
    grid_xy <- as.data.frame(terra::crds(grid_pts))
    names(grid_xy) <- c("x","y")
    grid_sf <- sf::st_as_sf(grid_xy, coords = c("x","y"),
                            crs = sf::st_crs(terra::crs(raster_template, proj = TRUE)))
    # IDW prediction
    pred_sf <- gstat::idw(value ~ 1, pts_sf, newdata = grid_sf, nmax = 12)
    pred_vec <- terra::vect(pred_sf) # back to SpatVector
    terra::rasterize(pred_vec, raster_template, field = "var1.pred")
  }, error = function(e) {
    log_event(sprintf("IDW failed (%s) -> using nearest neighbor fallback", e$message), "WARNING")
    terra::rasterize(pts, raster_template, field = "value", touches = TRUE)
  })
  
  # Mask to basin
  result_rast <- mask(result_rast, basin_boundary)
  
  # Coverage report
  valid_cells <- sum(!is.na(values(result_rast)))
  coverage_pct <- valid_cells / ncell(result_rast) * 100
  log_event(sprintf("Raster created: %d valid cells (%.1f%% coverage)", valid_cells, coverage_pct), "SUCCESS")
  
  return(result_rast)
}

# ==============================================================================
# HELPER: WRAP LEGEND LABELS
# ==============================================================================
wrap_legend_label <- function(label, max_width = 25) {
  if (nchar(label) <= max_width) return(label)
  words <- str_split(label, "\\s+")[[1]]
  lines <- character()
  current_line <- ""
  for (word in words) {
    test_line <- if (current_line == "") word else paste(current_line, word)
    if (nchar(test_line) > max_width && current_line != "") {
      lines <- c(lines, current_line)
      current_line <- word
    } else {
      current_line <- test_line
    }
  }
  if (current_line != "") lines <- c(lines, current_line)
  paste(lines, collapse = "\n")
}

# ==============================================================================
# MAP CREATION (continuous legend) -- PNGs, titles in outer margin; no clipping
# ==============================================================================
create_map_with_continuous_legend <- function(raster_data, basin_boundary, output_file,
                                              title, subtitle = NULL,
                                              color_palette = NULL, n_colors = 9,
                                              legend_title = "Value") {
  valid_vals  <- values(raster_data)
  valid_cells <- sum(!is.na(valid_vals))
  if (valid_cells == 0) {
    log_event(sprintf("Warning: No valid data for %s", basename(output_file)), "WARNING")
    return(FALSE)
  }
  data_range <- range(valid_vals, na.rm = TRUE)
  log_event(sprintf("Creating map: %d valid cells, range [%.3f, %.3f]",
                    valid_cells, data_range[1], data_range[2]), "INFO")
  
  if (is.null(color_palette)) {
    color_palette <- rev(brewer.pal(n_colors, "RdYlBu"))
  }
  
  png(output_file, width = 16, height = 9, units = "in", res = 300)
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.83, 0.17))
  
  # Reserve outer margin for titles; small inner margins so no overlap
  par(oma = c(2.2, 2, 4.5, 2))     # outer margins (bottom, left, top, right)
  par(mar = c(3, 2, 3, 3))         # inner margins
  
  # Panel 1 – map
  plot(raster_data, col = color_palette, axes = FALSE, box = FALSE, legend = FALSE, main = "")
  plot(basin_boundary, add = TRUE, border = "black", lwd = 2.5)
  tryCatch({
    sbar(10000, xy = "bottomleft", type = "bar", below = "meters", divs = 4, lwd = 2)
    north(xy = "topleft", type = 1, cex = 1.5)
  }, error = function(e) NULL)
  
  metadata_text <- sprintf("Period: 1950-2025 (76 years)  \u03B1 = %.2f", alpha)
  mtext(metadata_text, side = 1, line = 1.0, adj = 0.05, cex = 0.95, col = "gray40")
  
  # Panel 2 – legend
  par(mar = c(3, 1, 3, 2))
  plot.new()
  legend_x_left  <- 0.1; legend_x_right <- 0.45
  legend_y_bottom <- 0.15; legend_y_top <- 0.85
  n_grad <- length(color_palette)
  y_seq <- seq(legend_y_bottom, legend_y_top, length.out = n_grad + 1)
  for (i in 1:n_grad) {
    rect(legend_x_left, y_seq[i], legend_x_right, y_seq[i + 1], col = color_palette[i], border = NA, xpd = TRUE)
  }
  rect(legend_x_left, legend_y_bottom, legend_x_right, legend_y_top, col = NA, border = "black", lwd = 2)
  
  # Legend labels (no clipping)
  n_labels <- 5
  label_vals <- seq(data_range[1], data_range[2], length.out = n_labels)
  label_y_pos <- seq(legend_y_bottom, legend_y_top, length.out = n_labels)
  old_xpd <- par(xpd = NA)
  for (i in 1:n_labels) {
    text(legend_x_right + 0.08, label_y_pos[i], sprintf("%.2f", label_vals[i]), adj = 0, cex = 1.25, font = 1)
  }
  par(xpd = old_xpd)
  text((legend_x_left + legend_x_right) / 2, legend_y_top + 0.09, legend_title, cex = 1.6, font = 2, col = "black")
  
  # Titles in outer margin (prevents any overlap with map)
  mtext(title, side = 3, outer = TRUE, line = 2.6, cex = 1.8, font = 2)
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, outer = TRUE, line = 1.2, cex = 1.15, font = 1, col = "gray30")
  }
  
  dev.off()
  log_event(sprintf("Created continuous legend map: %s", basename(output_file)), "SUCCESS")
  return(TRUE)
}

# ==============================================================================
# MAP CREATION (categorical legend) -- PNGs, titles in outer margin; no clipping
# ==============================================================================
create_map_with_categorical_legend <- function(raster_data, basin_boundary, output_file,
                                               title, subtitle = NULL,
                                               legend_labels, color_palette) {
  valid_vals  <- values(raster_data)
  valid_cells <- sum(!is.na(valid_vals))
  if (valid_cells == 0) {
    log_event(sprintf("Warning: No valid data for %s", basename(output_file)), "WARNING")
    return(FALSE)
  }
  log_event(sprintf("Creating categorical map: %d valid cells", valid_cells), "INFO")
  
  png(output_file, width = 16, height = 9, units = "in", res = 300)
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.83, 0.17))
  par(oma = c(2.2, 2, 4.5, 2))
  par(mar = c(3, 2, 3, 3))
  
  # Panel 1 – map
  plot(raster_data, col = color_palette, axes = FALSE, box = FALSE, legend = FALSE, main = "")
  plot(basin_boundary, add = TRUE, border = "black", lwd = 2.5)
  tryCatch({
    sbar(10000, xy = "bottomleft", type = "bar", below = "meters", divs = 4, lwd = 2)
    north(xy = "topleft", type = 1, cex = 1.5)
  }, error = function(e) NULL)
  
  metadata_text <- sprintf("Period: 1950-2025 (76 years)  \u03B1 = %.2f", alpha)
  mtext(metadata_text, side = 1, line = 1.0, adj = 0.05, cex = 0.95, col = "gray40")
  
  # Panel 2 – legend
  par(mar = c(3, 1, 3, 2))
  plot.new()
  n_cats <- length(legend_labels)
  cat_height <- 0.14
  total_height <- n_cats * cat_height
  y_start <- 0.5 - total_height / 2
  legend_x_left <- 0.1; legend_x_right <- 0.45
  for (i in 1:n_cats) {
    y_b <- y_start + (i - 1) * cat_height
    y_t <- y_b + cat_height
    rect(legend_x_left, y_b, legend_x_right, y_t, col = color_palette[i], border = "black", lwd = 2)
  }
  
  # Labels (no clipping)
  old_xpd <- par(xpd = NA)
  for (i in 1:n_cats) {
    y_pos <- y_start + (i - 0.5) * cat_height
    wrapped_label <- wrap_legend_label(legend_labels[i], max_width = 22)
    text(legend_x_right + 0.08, y_pos, wrapped_label, adj = 0, cex = 1.15, font = 1)
  }
  par(xpd = old_xpd)
  
  text((legend_x_left + legend_x_right) / 2, y_start + total_height + 0.09,
       "Classification", cex = 1.6, font = 2)
  
  # Titles in outer margin
  mtext(title, side = 3, outer = TRUE, line = 2.6, cex = 1.8, font = 2)
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, outer = TRUE, line = 1.2, cex = 1.15, font = 1, col = "gray30")
  }
  
  dev.off()
  log_event(sprintf("Created categorical legend map: %s", basename(output_file)), "SUCCESS")
  return(TRUE)
}

# ==============================================================================
# PDF EXPORT - CONTINUOUS
#    (NEW) add_footer switch -> FALSE for drought frequency PDFs (no bottom text)
# ==============================================================================
create_map_pdf_continuous <- function(raster_data, basin_boundary, output_file,
                                      title, subtitle = NULL,
                                      color_palette = NULL, n_colors = 9,
                                      legend_title = "Value",
                                      add_footer = TRUE) {
  valid_vals <- values(raster_data)
  valid_cells <- sum(!is.na(valid_vals))
  if (valid_cells == 0) {
    log_event(sprintf("Warning: No valid data for PDF %s", basename(output_file)), "WARNING")
    return(FALSE)
  }
  data_range <- range(valid_vals, na.rm = TRUE)
  if (is.null(color_palette)) {
    color_palette <- rev(brewer.pal(n_colors, "RdYlBu"))
  }
  
  # A4 landscape (11.7 x 8.3 inches)
  pdf(output_file, width = 11.7, height = 8.3, family = "Helvetica")
  par(mar = c(4.5, 4.5, 3.5, 2.5)) # Bottom, left, top, right
  
  # Raster plot
  plot(raster_data, col = color_palette, axes = TRUE, box = TRUE,
       main = "", xlab = "Easting (m)", ylab = "Northing (m)")
  
  # Title & subtitle
  title(main = title, line = 2.8, cex.main = 1.7, font.main = 2)
  if (!is.null(subtitle)) {
    title(main = subtitle, line = 1.3, cex.main = 1.1, font.main = 1, col.main = "gray30")
  }
  
  # Basin boundary
  plot(basin_boundary, add = TRUE, border = "black", lwd = 2.5)
  
  # Continuous legend (right side)
  legend("right",
         legend = sprintf("%.2f", seq(data_range[1], data_range[2], length.out = 5)),
         fill = color_palette,
         bty = "n", cex = 1.1, title = legend_title)
  
  # Footer (omit when add_footer = FALSE)
  if (isTRUE(add_footer)) {
    metadata_text <- sprintf("Nechako Basin  \nPeriod: 1950-2025 (76 years)  \n\u03B1 = %.2f",
                             alpha)
    mtext(metadata_text, side = 1, line = 2.8, cex = 0.85, col = "gray40")
  }
  
  dev.off()
  log_event(sprintf("Created PDF map: %s", basename(output_file)), "SUCCESS")
  return(TRUE)
}

# ==============================================================================
# PDF EXPORT - CATEGORICAL (Runs/VC/TFPW)
#    (NEW) Use two-panel layout with custom legend panel to avoid trimming
# ==============================================================================
create_map_pdf_categorical <- function(raster_data, basin_boundary, output_file,
                                       title, subtitle = NULL,
                                       legend_labels, color_palette) {
  valid_vals <- values(raster_data)
  valid_cells <- sum(!is.na(valid_vals))
  if (valid_cells == 0) {
    log_event(sprintf("Warning: No valid data for PDF %s", basename(output_file)), "WARNING")
    return(FALSE)
  }
  
  pdf(output_file, width = 11.7, height = 8.3, family = "Helvetica")
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.83, 0.17))
  
  # Panel 1 — map
  par(mar = c(4.5, 4.5, 3.5, 1.0))
  plot(raster_data, col = color_palette, axes = TRUE, box = TRUE,
       main = "", xlab = "Easting (m)", ylab = "Northing (m)")
  title(main = title, line = 2.8, cex.main = 1.7, font.main = 2)
  if (!is.null(subtitle)) {
    title(main = subtitle, line = 1.3, cex.main = 1.1, font.main = 1, col.main = "gray30")
  }
  plot(basin_boundary, add = TRUE, border = "black", lwd = 2.5)
  
  # Panel 2 — legend (custom, no trimming)
  par(mar = c(4.5, 1.0, 3.5, 2.5))
  plot.new()
  n_cats <- length(legend_labels)
  cat_height <- 0.14
  total_height <- n_cats * cat_height
  y_start <- 0.5 - total_height / 2
  legend_x_left <- 0.1; legend_x_right <- 0.45
  for (i in 1:n_cats) {
    y_b <- y_start + (i - 1) * cat_height
    y_t <- y_b + cat_height
    rect(legend_x_left, y_b, legend_x_right, y_t, col = color_palette[i], border = "black", lwd = 2)
  }
  old_xpd <- par(xpd = NA)
  for (i in 1:n_cats) {
    y_pos <- y_start + (i - 0.5) * cat_height
    wrapped_label <- wrap_legend_label(legend_labels[i], max_width = 22)
    text(legend_x_right + 0.08, y_pos, wrapped_label, adj = 0, cex = 1.15)
  }
  par(xpd = old_xpd)
  text((legend_x_left + legend_x_right) / 2, y_start + total_height + 0.09,
       "Classification", cex = 1.6, font = 2)
  
  dev.off()
  log_event(sprintf("Created categorical PDF map: %s", basename(output_file)), "SUCCESS")
  return(TRUE)
}

# ==============================================================================
# MAIN PROCESSING LOOP
# ==============================================================================
log_event("Searching for monthly composite CSV files...", "INFO")
index_types <- c("spi", "spei")
timescales <- c(1, 3, 6, 9, 12, 24)

all_results <- list()
total_original_pixels <- 0
total_basin_pixels <- 0
start_time <- Sys.time()

for (index_type in index_types) {
  subfolder_name <- paste0(index_type, "_results_seasonal")
  for (scale_idx in seq_along(timescales)) {
    scale <- timescales[scale_idx]
    
    reconstruction <- read_monthly_csv_and_reconstruct(index_type, scale, subfolder_name)
    if (is.null(reconstruction)) {
      log_event(sprintf("Skipping %s-%02d due to insufficient files", index_type, scale), "WARNING")
      next
    }
    
    # Track original pixel count
    total_original_pixels <- total_original_pixels + reconstruction$n_pixels
    
    # Decide source CRS once per dataset (auto-detect from provided coords)
    coords_all <- reconstruction$coords
    looks_like_deg_all <- (max(abs(coords_all$lon), na.rm=TRUE) <= 180) &&
      (max(abs(coords_all$lat), na.rm=TRUE) <= 90)
    this_source_crs <- if (looks_like_deg_all) "EPSG:4326" else "EPSG:3005"
    
    # Filter to pixels with sufficient valid observations BEFORE stats
    log_event("Filtering to basin pixels with sufficient valid observations...", "INFO")
    filtered <- filter_basin_pixels(reconstruction$index_matrix, reconstruction$coords, min_valid_obs)
    if (filtered$n_pixels == 0) {
      log_event(sprintf("Skipping %s-%02d: No pixels with >=%d valid observations",
                        index_type, scale, min_valid_obs), "WARNING")
      next
    }
    
    # Track basin pixel count
    total_basin_pixels <- total_basin_pixels + filtered$n_pixels
    log_event(sprintf("Pixel reduction: %d → %d pixels (%.1f%% reduction)",
                      reconstruction$n_pixels, filtered$n_pixels, filtered$reduction_pct), "SUCCESS")
    
    # Basin-only data
    index_matrix <- filtered$index_matrix
    coords <- filtered$coords
    n_pixels <- filtered$n_pixels
    n_months <- nrow(index_matrix)
    
    # Adjust workers based on basin pixels
    optimal_workers <- min(num_cores, max(1, floor(n_pixels / 50)))
    if (optimal_workers < num_cores) {
      log_event(sprintf("Adjusting workers: %d → %d (based on %d basin pixels)",
                        num_cores, optimal_workers, n_pixels), "INFO")
      plan(multisession, workers = optimal_workers)
    }
    
    log_event(sprintf("Processing %s-%02d (%d basin pixels, %d months)...",
                      toupper(index_type), scale, n_pixels, n_months), "INFO")
    
    coords_dt <- data.table(
      space_idx = 1:n_pixels,
      lon = coords$lon,
      lat = coords$lat
    )
    
    # ---- Full Monthly analyses ----
    log_event("  Running annual VC trend analysis...", "INFO")
    vc_annual <- modified_mann_kendall_taub(index_matrix, alpha, max_tie_percent, apply_vc = TRUE)
    
    log_event("  Running annual TFPW trend analysis...", "INFO")
    tfpw_annual <- perform_tfpw_mk_taub(index_matrix, alpha, max_tie_percent)
    
    log_event("  Running annual Wald-Wolfowitz Runs Test...", "INFO")
    runs_annual <- future_lapply(1:n_pixels, function(i) {
      wald_wolfowitz_runs_test(index_matrix[, i], drought_threshold)
    }, future.seed = TRUE)
    runs_df_annual <- rbindlist(lapply(runs_annual, as.data.frame))
    runs_df_annual[, space_idx := 1:n_pixels]
    
    log_event("  Calculating drought event metrics...", "INFO")
    event_metrics <- future_lapply(1:n_pixels, function(i) {
      ts_vals <- index_matrix[, i]
      detect_drought_events_pixel(ts_vals, drought_threshold, recovery_threshold, min_duration)
    }, future.seed = TRUE)
    event_df <- rbindlist(lapply(event_metrics, as.data.frame))
    event_df[, space_idx := 1:n_pixels]
    
    log_event("  Detecting regime shifts & return periods...", "INFO")
    regime_shifts <- future_sapply(1:n_pixels, function(i) {
      detect_regime_shift(index_matrix[, i])
    }, future.seed = TRUE)
    return_periods <- future_sapply(1:n_pixels, function(i) {
      calculate_return_period(index_matrix[, i], extreme_threshold)
    }, future.seed = TRUE)
    
    # Combine
    annual_results <- cbind(
      coords_dt,
      index_type = index_type,
      timescale = scale,
      period = "annual",
      month = NA,
      setDT(vc_annual)[, .(tau_vc = tau, p_value_vc = p.value, sl_vc = sl,
                           vc_corrected = vc_corrected, n_ties_vc = n_ties,
                           percent_ties_vc = percent_ties, tau_b_adjusted_vc = tau_b_adjusted,
                           filtered_vc = filtered, n_vc = n, rho1_vc = rho1)],
      setDT(tfpw_annual)[, .(tau_tfpw = tau, p_value_tfpw = p.value, sl_tfpw = sl,
                             tfpw_applied = tfpw_applied, n_ties_tfpw = n_ties,
                             percent_ties_tfpw = percent_ties, tau_b_adjusted_tfpw = tau_b_adjusted,
                             filtered_tfpw = filtered, n_tfpw = n, rho1_tfpw = rho1)],
      runs_df_annual[, .(n_runs, n_drought, n_non_drought, expected_runs, var_runs,
                         z_stat, p_value_runs = p_value, clustering, filtered_runs = filtered)],
      event_df[, .(n_events, mean_duration, mean_severity, mean_intensity,
                   max_duration, max_severity, max_intensity, total_severity)],
      regime_shift_year = regime_shifts,
      return_period_extreme = return_periods
    )
    
    # Derived flags
    annual_results[, runs_significant := !filtered_runs & p_value_runs < alpha]
    annual_results[, same_significance := (p_value_vc < alpha & !filtered_vc) ==
                     (p_value_tfpw < alpha & !filtered_tfpw)]
    annual_results[, same_direction := sign(tau_vc) == sign(tau_tfpw)]
    
    # Template (coarse) for interpolation
    template <- rast(ext(basin), nrows = 100, ncols = 100, crs = target_crs)
    
    # --------------------------------------------------------------------------
    # Create rasters & maps (TIFs + PDFs in out_dir; PNGs in maps_with_legends)
    # --------------------------------------------------------------------------
    log_event("  Creating spatial rasters and maps...", "INFO")
    
    # 1) Runs test significance map
    log_event("   • Creating runs test map...", "INFO")
    annual_results[, clustering_numeric := fifelse(clustering == "clustering", -1,
                                                   fifelse(clustering == "dispersion", 1, NA_real_))]
    annual_results_sig_runs <- annual_results[runs_significant == TRUE & !is.na(clustering_numeric)]
    if (nrow(annual_results_sig_runs) > 0) {
      runs_rast <- create_interpolated_raster_optimized(
        annual_results_sig_runs, "clustering_numeric", template, basin, source_crs = this_source_crs
      )
      
      # --- TIF EXPORT ---
      tif_path <- file.path(out_dir, sprintf("%s_%02d_annual_runs_significance.tif", index_type, scale))
      if (file.exists(tif_path)) file.remove(tif_path)
      writeRaster(runs_rast, tif_path, overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
      tif_check <- rast(tif_path)
      valid_pct <- 100 * sum(!is.na(values(tif_check))) / ncell(tif_check)
      log_event(sprintf("TIF validation: %.1f%% valid cells in %s", valid_pct, basename(tif_path)), "INFO")
      
      # --- PDF EXPORT (to out_dir) ---
      if (sum(!is.na(values(runs_rast))) > 0) {
        pdf_path <- file.path(out_dir, sprintf("%s_%02d_annual_runs_significance.pdf", index_type, scale))
        if (file.exists(pdf_path)) file.remove(pdf_path)
        create_map_pdf_categorical(
          runs_rast, basin, pdf_path,
          sprintf("%s-%02d: Wald-Wolfowitz Runs Test", toupper(index_type), scale),
          sprintf("Temporal Clustering (p < %.2f) \n Threshold: %.2f", alpha, drought_threshold),
          c("Clustering\n(Too few runs)", "Dispersion\n(Too many runs)"),
          c("#d62728", "#1f77b4")
        )
      }
      
      # --- PNG EXPORT ---
      if (sum(!is.na(values(runs_rast))) > 0) {
        create_map_with_categorical_legend(
          runs_rast, basin,
          file.path(maps_dir, sprintf("%s_%02d_runs_test_map.png", index_type, scale)),
          sprintf("%s-%02d: Wald-Wolfowitz Runs Test", toupper(index_type), scale),
          sprintf("Temporal Clustering Analysis (p < %.2f) \n Threshold: %.2f", alpha, drought_threshold),
          c("Clustering\n(Too few runs)", "Dispersion\n(Too many runs)"),
          c("#d62728", "#1f77b4")
        )
      }
    } else {
      log_event("   • No significant runs - skipping map", "WARNING")
    }
    
    # 2) VC trend significance map
    log_event("   • Creating VC trends map...", "INFO")
    annual_results[, tau_vc_sign := sign(tau_vc)]
    annual_results_sig_vc <- annual_results[p_value_vc < alpha & filtered_vc == FALSE]
    if (nrow(annual_results_sig_vc) > 0) {
      vc_rast <- create_interpolated_raster_optimized(
        annual_results_sig_vc, "tau_vc_sign", template, basin, source_crs = this_source_crs
      )
      
      # --- TIF EXPORT ---
      tif_path <- file.path(out_dir, sprintf("%s_%02d_annual_vc_significance.tif", index_type, scale))
      if (file.exists(tif_path)) file.remove(tif_path)
      writeRaster(vc_rast, tif_path, overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
      tif_check <- rast(tif_path)
      valid_pct <- 100 * sum(!is.na(values(tif_check))) / ncell(tif_check)
      log_event(sprintf("TIF validation: %.1f%% valid cells in %s", valid_pct, basename(tif_path)), "INFO")
      
      # --- PDF EXPORT (to out_dir) ---
      if (sum(!is.na(values(vc_rast))) > 0) {
        pdf_path <- file.path(out_dir, sprintf("%s_%02d_annual_vc_significance.pdf", index_type, scale))
        if (file.exists(pdf_path)) file.remove(pdf_path)
        create_map_pdf_categorical(
          vc_rast, basin, pdf_path,
          sprintf("%s-%02d: Mann-Kendall (VC Method)", toupper(index_type), scale),
          sprintf("Variance-Corrected Trend Analysis (p < %.2f)", alpha),
          c("Drying Trend", "Wetting Trend"),
          c("#8B4513", "#228B22")
        )
      }
      
      # --- PNG EXPORT ---
      if (sum(!is.na(values(vc_rast))) > 0) {
        create_map_with_categorical_legend(
          vc_rast, basin,
          file.path(maps_dir, sprintf("%s_%02d_vc_trends_map.png", index_type, scale)),
          sprintf("%s-%02d: Mann-Kendall (VC Method)", toupper(index_type), scale),
          sprintf("Variance-Corrected Trend Analysis (p < %.2f)", alpha),
          c("Drying Trend", "Wetting Trend"),
          c("#8B4513", "#228B22")
        )
      }
    } else {
      log_event("   • No significant VC trends - skipping map", "WARNING")
    }
    
    # 3) TFPW trend significance map
    log_event("   • Creating TFPW trends map...", "INFO")
    annual_results[, tau_tfpw_sign := sign(tau_tfpw)]
    annual_results_sig_tfpw <- annual_results[p_value_tfpw < alpha & filtered_tfpw == FALSE]
    if (nrow(annual_results_sig_tfpw) > 0) {
      tfpw_rast <- create_interpolated_raster_optimized(
        annual_results_sig_tfpw, "tau_tfpw_sign", template, basin, source_crs = this_source_crs
      )
      
      # --- TIF EXPORT ---
      tif_path <- file.path(out_dir, sprintf("%s_%02d_annual_tfpw_significance.tif", index_type, scale))
      if (file.exists(tif_path)) file.remove(tif_path)
      writeRaster(tfpw_rast, tif_path, overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
      tif_check <- rast(tif_path)
      valid_pct <- 100 * sum(!is.na(values(tif_check))) / ncell(tif_check)
      log_event(sprintf("TIF validation: %.1f%% valid cells in %s", valid_pct, basename(tif_path)), "INFO")
      
      # --- PDF EXPORT (to out_dir) ---
      if (sum(!is.na(values(tfpw_rast))) > 0) {
        pdf_path <- file.path(out_dir, sprintf("%s_%02d_annual_tfpw_significance.pdf", index_type, scale))
        if (file.exists(pdf_path)) file.remove(pdf_path)
        create_map_pdf_categorical(
          tfpw_rast, basin, pdf_path,
          sprintf("%s-%02d: Mann-Kendall (TFPW Method)", toupper(index_type), scale),
          sprintf("Trend-Free Pre-Whitened Analysis (p < %.2f)", alpha),
          c("Drying Trend", "Wetting Trend"),
          c("#8B4513", "#228B22")
        )
      }
      
      # --- PNG EXPORT ---
      if (sum(!is.na(values(tfpw_rast))) > 0) {
        create_map_with_categorical_legend(
          tfpw_rast, basin,
          file.path(maps_dir, sprintf("%s_%02d_tfpw_trends_map.png", index_type, scale)),
          sprintf("%s-%02d: Mann-Kendall (TFPW Method)", toupper(index_type), scale),
          sprintf("Trend-Free Pre-Whitened Analysis (p < %.2f)", alpha),
          c("Drying Trend", "Wetting Trend"),
          c("#8B4513", "#228B22")
        )
      }
    } else {
      log_event("   • No significant TFPW trends - skipping map", "WARNING")
    }
    
    # 4) Drought frequency map (continuous)
    log_event("   • Creating drought frequency map...", "INFO")
    annual_results[, freq_per_decade := n_events / 7.6]
    freq_rast <- create_interpolated_raster_optimized(
      annual_results, "freq_per_decade", template, basin, source_crs = this_source_crs
    )
    
    # --- TIF EXPORT ---
    tif_path <- file.path(out_dir, sprintf("%s_%02d_drought_frequency_per_decade.tif", index_type, scale))
    if (file.exists(tif_path)) file.remove(tif_path)
    writeRaster(freq_rast, tif_path, overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
    tif_check <- rast(tif_path)
    valid_pct <- 100 * sum(!is.na(values(tif_check))) / ncell(tif_check)
    log_event(sprintf("TIF validation: %.1f%% valid cells in %s", valid_pct, basename(tif_path)), "INFO")
    
    # --- PDF EXPORT (to out_dir) — NO FOOTER per request ---
    if (sum(!is.na(values(freq_rast))) > 0) {
      pdf_path <- file.path(out_dir, sprintf("%s_%02d_drought_frequency_per_decade.pdf", index_type, scale))
      if (file.exists(pdf_path)) file.remove(pdf_path)
      create_map_pdf_continuous(
        freq_rast, basin, pdf_path,
        sprintf("%s-%02d: Drought Frequency", toupper(index_type), scale),
        "Number of drought events per decade",
        color_palette = c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c",
                          "#fc4e2a", "#e31a1c", "#bd0026", "#800026"),
        legend_title = "Events/decade",
        add_footer = FALSE   # <<< remove bottom descriptions/date
      )
    }
    
    # --- PNG EXPORT ---
    if (sum(!is.na(values(freq_rast))) > 0) {
      create_map_with_continuous_legend(
        freq_rast, basin,
        file.path(maps_dir, sprintf("%s_%02d_drought_frequency_map.png", index_type, scale)),
        sprintf("%s-%02d: Drought Frequency", toupper(index_type), scale),
        "Number of drought events per decade",
        color_palette = c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c",
                          "#fc4e2a", "#e31a1c", "#bd0026", "#800026"),
        legend_title = "Events/decade"
      )
    }
    
    # Store results
    all_results[[paste0(index_type, "_", scale)]] <- annual_results
    
    # Restore original worker count
    if (optimal_workers != num_cores) {
      plan(multisession, workers = num_cores)
    }
    
    # Progress tracking
    completed_idx <- (which(index_types == index_type) - 1) * length(timescales) + scale_idx
    total_tasks <- length(index_types) * length(timescales)
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    remaining <- (elapsed / completed_idx) * (total_tasks - completed_idx)
    log_event(sprintf("   Completed %s-%02d (%d/%d)\n   Basin pixels: %d\n   Elapsed: %.1f min\n   Remaining: %.1f min",
                      toupper(index_type), scale, completed_idx, total_tasks,
                      n_pixels, as.numeric(elapsed), as.numeric(remaining)), "SUCCESS")
  }
}

# ==============================================================================
# OPTIONAL CLEANUP: Remove any legacy PDFs from maps_with_legends
# ==============================================================================
old_pdfs <- list.files(maps_dir, pattern = "\\.pdf$", full.names = TRUE)
if (length(old_pdfs)) {
  log_event(sprintf("Removing %d legacy PDF(s) from maps_with_legends...", length(old_pdfs)), "INFO")
  file.remove(old_pdfs)
}

# ==============================================================================
# Ensure every TIF in temporal_spi_spei has a PDF sibling (preview)
# ==============================================================================
log_event("Ensuring PDF previews for all TIFs in temporal_spi_spei...", "INFO")
tif_files <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
for (tp in tif_files) {
  pdf_p <- sub("\\.tif$", ".pdf", tp)
  if (!file.exists(pdf_p)) {
    log_event(sprintf("Creating PDF preview for %s", basename(tp)), "INFO")
    rtry <- try({
      rr <- rast(tp)
      pdf(pdf_p, width = 11.7, height = 8.3, family = "Helvetica")
      par(mar = c(4.5, 4.5, 3.5, 2.5))
      plot(rr, main = tools::file_path_sans_ext(basename(tp)), axes = TRUE, box = TRUE)
      plot(basin, add = TRUE, border = "black", lwd = 2.0)
      dev.off()
    }, silent = TRUE)
    if (inherits(rtry, "try-error")) {
      log_event(sprintf("Failed to build PDF for %s", basename(tp)), "WARNING")
      if (file.exists(pdf_p)) try(file.remove(pdf_p), silent = TRUE)
    }
  }
}

# ==============================================================================
# SAVE FINAL RESULTS + SUMMARY
# ==============================================================================
log_event("Saving consolidated results...", "INFO")
if (length(all_results) > 0) {
  final_results <- rbindlist(all_results, fill = TRUE)
  fwrite(final_results, file.path(out_dir, "all_temporal_diagnostics_optimized.csv"))
  
  # Summary with optimization metrics
  total_pixels <- nrow(final_results)
  sig_vc <- sum(final_results$p_value_vc < alpha & !final_results$filtered_vc, na.rm = TRUE)
  sig_tfpw <- sum(final_results$p_value_tfpw < alpha & !final_results$filtered_tfpw, na.rm = TRUE)
  sig_runs <- sum(final_results$runs_significant, na.rm = TRUE)
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  reduction_pct <- 100 * (1 - total_basin_pixels / total_original_pixels)
  estimated_original_time <- if (reduction_pct < 100) total_time / (1 - reduction_pct/100) else total_time
  
  cat("\n=== OPTIMIZATION SUMMARY ===\n", file = LOG_FILE, append = TRUE)
  cat(sprintf("Original pixels processed: %d\n", total_original_pixels), file = LOG_FILE, append = TRUE)
  cat(sprintf("Basin pixels processed: %d\n", total_basin_pixels), file = LOG_FILE, append = TRUE)
  cat(sprintf("Spatial reduction: %.1f%%\n", reduction_pct), file = LOG_FILE, append = TRUE)
  cat(sprintf("Total processing time: %.1f minutes\n", total_time), file = LOG_FILE, append = TRUE)
  cat(sprintf("Estimated time w/o opt: %.1f minutes\n", estimated_original_time), file = LOG_FILE, append = TRUE)
  cat(sprintf("Time saved: %.1f minutes (%.1f%% reduction)\n",
              estimated_original_time - total_time,
              100 * (estimated_original_time - total_time) / estimated_original_time),
      file = LOG_FILE, append = TRUE)
  
  cat("\n=== ANALYSIS RESULTS ===\n", file = LOG_FILE, append = TRUE)
  cat(sprintf("Total basin pixels analyzed: %d\n", total_pixels), file = LOG_FILE, append = TRUE)
  cat(sprintf("VC significant trends: %d (%.1f%%)\n", sig_vc, sig_vc/total_pixels*100), file = LOG_FILE, append = TRUE)
  cat(sprintf("TFPW significant trends: %d (%.1f%%)\n", sig_tfpw, sig_tfpw/total_pixels*100), file = LOG_FILE, append = TRUE)
  cat(sprintf("Significant runs patterns: %d (%.1f%%)\n", sig_runs, sig_runs/total_pixels*100), file = LOG_FILE, append = TRUE)
  cat("========================\n\n", file = LOG_FILE, append = TRUE)
  
  log_event(sprintf("SUCCESS! %.1f%% spatial reduction achieved", reduction_pct), "SUCCESS")
  log_event(sprintf("Total runtime: %.1f minutes (saved %.1f minutes vs non-optimized)",
                    total_time, estimated_original_time - total_time), "SUCCESS")
} else {
  log_event("No results generated", "WARNING")
}

plan(sequential)
log_event("Analysis completed successfully", "SUCCESS")
cat("\n✅ ANALYSIS COMPLETE - All TIFs and PDF siblings in temporal_spi_spei; PNGs in maps_with_legends only!\n",
    file = LOG_FILE, append = TRUE)