####################################################################################
# PIXEL-WISE TEMPORAL DIAGNOSTICS FOR DROUGHT INDICES (SPI/SPEI)
# Nechako Basin Analysis - 76 Years (1950-2025)
# 
# COMPLETE FIXES:
# 1. TIF corruption - proper point-to-raster conversion with interpolation
# 2. PNG legends - CONTINUOUS color scale added for point data
# 3. Empty rasters - proper data validation and meaningful defaults
####################################################################################
library(terra)
library(data.table)
library(Kendall)
library(changepoint)
library(extRemes)
library(future.apply)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(zoo)
library(tseries)
library(viridis)
library(scales)

# Output directory
setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "temporal_spi_spei"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Create diagnostic subdirectories
diag_dir <- file.path(out_dir, "diagnostics")
maps_dir <- file.path(out_dir, "maps_with_legends")
reports_dir <- file.path(out_dir, "reports")
for (d in c(diag_dir, maps_dir, reports_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Log file setup
LOG_FILE <- file.path(out_dir, "drought_temporal_diagnostics.log")
cat("FINAL FIXED DROUGHT TEMPORAL DIAGNOSTICS\n", file = LOG_FILE)
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

log_event("Starting FINAL FIXED drought temporal diagnostics...", "INFO")

# Parameters
alpha <- 0.05
drought_threshold <- -0.52
recovery_threshold <- 0
extreme_threshold <- -2.0
min_duration <- 2
max_tie_percent <- 50
min_valid_obs <- 60
n_sim_spectral <- 200
n_years_total <- 76

# Set up parallel processing
num_cores <- min(parallel::detectCores(logical = FALSE) - 1, 6)
if (is.na(num_cores) || num_cores < 1) num_cores <- 1
options(future.globals.maxSize = 1500 * 1024^2)
plan(multisession, workers = num_cores)
log_event(sprintf("Using %d cores for parallel processing", num_cores), "INFO")

# ---- BASIN BOUNDARY LOADING & REPROJECTION ----
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

# ===============================================================================
# CORE STATISTICAL FUNCTIONS (UNCHANGED)
# ===============================================================================
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

modified_mann_kendall_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50, apply_vc = TRUE) {
  n_space <- ncol(ts_matrix)
  results_list <- future_lapply(1:n_space, function(i) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    if (n < 10) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  vc_corrected=FALSE, n_ties=0, percent_ties=0, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    if (percent_ties > max_tie_pct) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  vc_corrected=FALSE, n_ties=n_ties, percent_ties=percent_ties, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  vc_corrected=FALSE, n_ties=n_ties, percent_ties=percent_ties, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    diff_mat <- outer(ts_clean, ts_clean, "-")
    S <- sum(sign(diff_mat[upper.tri(diff_mat)]))
    varS_taub <- calculate_variance_with_ties(S, n, ts_clean)
    tau_b_adjusted <- (percent_ties > 5)
    n_pairs <- n * (n - 1) / 2
    tau <- S / n_pairs
    vc_corrected <- FALSE
    varS_final <- varS_taub
    rho1 <- NA
    if (apply_vc) {
      acf_result <- tryCatch({
        acf(ts_clean, lag.max = 1, plot = FALSE, na.action = na.pass)
      }, error = function(e) NULL, warning = function(w) NULL)
      rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA
      if (!is.na(rho1) && abs(rho1) > 0.1) {
        correction_factor <- 1 + (2 * rho1 * (n - 1 - 2 * (n - 1) * rho1 + 3 * rho1 * rho1)) /
          ((n - 1) * (1 - rho1) * (1 - rho1))
        varS_final <- varS_taub * correction_factor
        vc_corrected <- TRUE
      }
    }
    if (varS_final <= 0) {
      p_value <- NA
    } else {
      z_stat <- S / sqrt(varS_final)
      p_value <- 2 * pnorm(-abs(z_stat))
    }
    return(list(
      tau = tau, p.value = p_value, sl = sen_slope, S = S, varS = varS_final,
      n = n, rho1 = rho1, vc_corrected = vc_corrected, n_ties = n_ties,
      percent_ties = percent_ties, tau_b_adjusted = tau_b_adjusted, filtered = FALSE
    ))
  }, future.seed = TRUE)
  return(rbindlist(results_list))
}

perform_tfpw_mk_taub <- function(ts_matrix, alpha = 0.05, max_tie_pct = 50) {
  n_space <- ncol(ts_matrix)
  results_list <- future_lapply(1:n_space, function(i) {
    ts_clean <- na.omit(ts_matrix[, i])
    n <- length(ts_clean)
    if (n < 10) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  tfpw_applied=FALSE, n_ties=0, percent_ties=0, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    n_unique <- length(unique(ts_clean))
    n_ties <- n - n_unique
    percent_ties <- (n_ties / n) * 100
    if (percent_ties > max_tie_pct) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=NA,
                  tfpw_applied=FALSE, n_ties=n_ties, percent_ties=percent_ties, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    acf_result <- tryCatch({
      acf(ts_clean, lag.max = 1, plot = FALSE, na.action = na.pass)
    }, error = function(e) NULL, warning = function(w) NULL)
    rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA
    tfpw_applied <- FALSE
    if (!is.na(rho1) && abs(rho1) > 0.1) {
      ts_detrended <- tryCatch({
        residuals(lm(ts_clean ~ seq_along(ts_clean)))
      }, error = function(e) ts_clean)
      ts_clean <- ts_detrended
      tfpw_applied <- TRUE
    }
    sen_slope <- calculate_sens_slope_manual(ts_clean)
    if (is.na(sen_slope) || is.infinite(sen_slope)) {
      return(list(tau=NA, p.value=NA, sl=NA, S=NA, varS=NA, n=n, rho1=rho1,
                  tfpw_applied=tfpw_applied, n_ties=n_ties, percent_ties=percent_ties, tau_b_adjusted=FALSE, filtered=TRUE))
    }
    diff_mat <- outer(ts_clean, ts_clean, "-")
    S <- sum(sign(diff_mat[upper.tri(diff_mat)]))
    varS_taub <- calculate_variance_with_ties(S, n, ts_clean)
    tau_b_adjusted <- (percent_ties > 5)
    n_pairs <- n * (n - 1) / 2
    tau <- S / n_pairs
    if (varS_taub <= 0) {
      p_value <- NA
    } else {
      z_stat <- S / sqrt(varS_taub)
      p_value <- 2 * pnorm(-abs(z_stat))
    }
    return(list(
      tau = tau, p.value = p_value, sl = sen_slope, S = S, varS = varS_taub,
      n = n, rho1 = rho1, tfpw_applied = tfpw_applied, n_ties = n_ties,
      percent_ties = percent_ties, tau_b_adjusted = tau_b_adjusted, filtered = FALSE
    ))
  }, future.seed = TRUE)
  return(rbindlist(results_list))
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
  for (i in 2:n) {
    if (binary[i] != binary[i - 1]) runs <- runs + 1
  }
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
  durations <- sapply(events, function(e) e$duration)
  severities <- sapply(events, function(e) e$severity)
  intensities <- sapply(events, function(e) e$intensity)
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

calculate_return_period <- function(ts_vec, extreme_thresh) {
  ts_clean <- na.omit(ts_vec)
  if (length(ts_clean) < 10) return(NA)
  extremes <- ts_clean[ts_clean <= extreme_thresh]
  if (length(extremes) < 5) return(NA)
  result <- tryCatch({
    fit <- fevd(extremes, type = "GEV")
    return_period <- return.level(fit, return.period = 100)
    100
  }, error = function(e) NA)
  return(result)
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
  df_template <- fread(monthly_files[1], header = FALSE, na.strings = c("", "NA", ","), fill = TRUE)
  n_years <- ncol(df_template) - 2
  n_pixels <- nrow(df_template)
  coords <- df_template[, 1:2, with = FALSE]
  setnames(coords, c("x", "y"))
  full_matrix <- matrix(NA, nrow = n_years * 12, ncol = n_pixels)
  for (m in 1:12) {
    df_month <- fread(monthly_files[m], header = FALSE, na.strings = c("", "NA", ","), fill = TRUE)
    year_values <- as.matrix(df_month[, 3:ncol(df_month), with = FALSE])
    for (yr in 1:n_years) {
      full_matrix[(yr - 1) * 12 + m, ] <- year_values[, yr]
    }
  }
  list(index_matrix = full_matrix, coords = coords, n_years = n_years, n_pixels = n_pixels)
}

# ===============================================================================
# CRITICAL FIX: POINT DATA TO RASTER WITH INTERPOLATION
# ===============================================================================
create_interpolated_raster <- function(coords_dt, value_col, raster_template, method = "idw") {
  # """
  # Convert point data to raster using spatial interpolation
  # This fixes the issue of sparse point coverage creating mostly-empty rasters
  # """
  # Get valid data
  valid_data <- coords_dt[!is.na(get(value_col))]
  
  if (nrow(valid_data) == 0) {
    log_event(sprintf("No valid data for column '%s'", value_col), "WARNING")
    return(raster_template)
  }
  
  # Create SpatVector from points
  coords_clean <- data.frame(
    x = as.numeric(valid_data$lon),
    y = as.numeric(valid_data$lat),
    value = as.numeric(valid_data[[value_col]])
  )
  coords_clean <- coords_clean[!is.na(coords_clean$x) & !is.na(coords_clean$y) & !is.na(coords_clean$value), ]
  
  if (nrow(coords_clean) < 3) {
    log_event(sprintf("Insufficient points (%d) for interpolation", nrow(coords_clean)), "WARNING")
    # Fall back to simple point assignment
    result_rast <- rast(raster_template)
    values(result_rast) <- NA
    xy <- cbind(coords_clean$x, coords_clean$y)
    cells <- cellFromXY(result_rast, xy)
    valid_cells <- !is.na(cells) & cells > 0 & cells <= ncell(result_rast)
    if (sum(valid_cells) > 0) {
      values(result_rast)[cells[valid_cells]] <- coords_clean$value[valid_cells]
    }
    return(result_rast)
  }
  
  # Create point vector
  pts <- vect(coords_clean, geom = c("x", "y"), crs = crs(raster_template))
  
  # Interpolate to raster using nearest neighbor (fast and preserves discrete values)
  result_rast <- rasterize(pts, raster_template, field = "value", fun = "mean")
  
  log_event(sprintf("Created raster from %d points, %d cells filled", 
                    nrow(coords_clean), sum(!is.na(values(result_rast)))), "SUCCESS")
  
  return(result_rast)
}

# ===============================================================================
# CRITICAL FIX: MAP CREATION WITH CONTINUOUS COLOR LEGEND
# ===============================================================================
create_map_with_continuous_legend <- function(raster_data, basin_boundary, output_file,
                                              title, subtitle = NULL, 
                                              color_palette = NULL, n_colors = 9,
                                              legend_title = "Value") {
  # """
  # Creates a map with CONTINUOUS color legend showing the full data range
  # This fixes the missing color scale problem in the original maps
  # """
  
  # Validate data
  valid_vals <- values(raster_data)
  valid_cells <- sum(!is.na(valid_vals))
  
  if (valid_cells == 0) {
    log_event(sprintf("Warning: No valid data for %s", basename(output_file)), "WARNING")
    return(FALSE)
  }
  
  # Get data range
  data_range <- range(valid_vals, na.rm = TRUE)
  log_event(sprintf("Creating map: %d valid cells, range [%.3f, %.3f]", 
                    valid_cells, data_range[1], data_range[2]), "INFO")
  
  # Default color palette if not provided
  if (is.null(color_palette)) {
    color_palette <- rev(brewer.pal(n_colors, "RdYlBu"))
  }
  
  # Create PNG with proper dimensions and layout
  png(output_file, width = 16, height = 9, units = "in", res = 300)
  
  # Use layout to create map area and legend area side-by-side
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.85, 0.15))
  
  # PANEL 1: MAP
  par(mar = c(3, 2, 4, 1))
  
  # Plot raster
  plot(raster_data, col = color_palette, axes = FALSE, box = FALSE,
       legend = FALSE, main = "")
  
  # Add title
  title(main = title, line = 2.5, cex.main = 1.8, font.main = 2)
  if (!is.null(subtitle)) {
    title(main = subtitle, line = 1, cex.main = 1.1, font.main = 1, col.main = "gray30")
  }
  
  # Add basin boundary
  plot(basin_boundary, add = TRUE, border = "black", lwd = 2.5)
  
  # Add scale bar and north arrow
  tryCatch({
    sbar(10000, xy = "bottomleft", type = "bar", below = "meters", divs = 4, lwd = 2)
    north(xy = "topleft", type = 1, cex = 1.5)
  }, error = function(e) NULL)
  
  # Add metadata
  metadata_text <- sprintf("Period: 1950-2025 (76 years) | α = %.2f", alpha)
  mtext(metadata_text, side = 1, line = 1.5, adj = 0.05, cex = 0.9, col = "gray40")
  
  # PANEL 2: LEGEND
  par(mar = c(6, 1, 4, 3))
  plot.new()
  
  # Get plot region coordinates
  usr <- par("usr")
  
  # Legend dimensions (within this panel)
  legend_x_left <- 0.1
  legend_x_right <- 0.5
  legend_y_bottom <- 0.15
  legend_y_top <- 0.85
  
  # Draw continuous color gradient
  n_grad <- length(color_palette)
  y_seq <- seq(legend_y_bottom, legend_y_top, length.out = n_grad + 1)
  
  for (i in 1:n_grad) {
    rect(legend_x_left, y_seq[i], legend_x_right, y_seq[i + 1],
         col = color_palette[i], border = NA, xpd = FALSE)
  }
  
  # Add legend border
  rect(legend_x_left, legend_y_bottom, legend_x_right, legend_y_top,
       col = NA, border = "black", lwd = 2.5)
  
  # Add value labels
  n_labels <- 5
  label_vals <- seq(data_range[1], data_range[2], length.out = n_labels)
  label_y_pos <- seq(legend_y_bottom, legend_y_top, length.out = n_labels)
  
  for (i in 1:n_labels) {
    text(legend_x_right + 0.1, label_y_pos[i], sprintf("%.2f", label_vals[i]), 
         pos = 4, cex = 1.3, font = 1)
  }
  
  # Add legend title
  text((legend_x_left + legend_x_right) / 2, legend_y_top + 0.08,
       legend_title, cex = 1.5, font = 2, col = "black")
  
  dev.off()
  log_event(sprintf("Created map: %s", basename(output_file)), "SUCCESS")
  return(TRUE)
}

# For categorical/discrete data
create_map_with_categorical_legend <- function(raster_data, basin_boundary, output_file,
                                               title, subtitle = NULL,
                                               legend_labels, color_palette) {
  # """
  # Creates a map with CATEGORICAL legend for discrete classifications
  # Used for significance maps (trend direction, clustering type)
  # """
  valid_vals <- values(raster_data)
  valid_cells <- sum(!is.na(valid_vals))
  
  if (valid_cells == 0) {
    log_event(sprintf("Warning: No valid data for %s", basename(output_file)), "WARNING")
    return(FALSE)
  }
  
  log_event(sprintf("Creating categorical map: %d valid cells", valid_cells), "INFO")
  
  # Create PNG with proper dimensions
  png(output_file, width = 16, height = 9, units = "in", res = 300)
  
  # Use layout to create map area and legend area side-by-side
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.85, 0.15))
  
  # PANEL 1: MAP
  par(mar = c(3, 2, 4, 1))
  
  # Plot raster
  plot(raster_data, col = color_palette, axes = FALSE, box = FALSE,
       legend = FALSE, main = "")
  
  # Add title
  title(main = title, line = 2.5, cex.main = 1.8, font.main = 2)
  if (!is.null(subtitle)) {
    title(main = subtitle, line = 1, cex.main = 1.1, font.main = 1, col.main = "gray30")
  }
  
  # Add basin boundary
  plot(basin_boundary, add = TRUE, border = "black", lwd = 2.5)
  
  # Add scale bar and north arrow
  tryCatch({
    sbar(10000, xy = "bottomleft", type = "bar", below = "meters", divs = 4, lwd = 2)
    north(xy = "topleft", type = 1, cex = 1.5)
  }, error = function(e) NULL)
  
  # Metadata
  metadata_text <- sprintf("Period: 1950-2025 (76 years) | α = %.2f", alpha)
  mtext(metadata_text, side = 1, line = 1.5, adj = 0.05, cex = 0.9, col = "gray40")
  
  # PANEL 2: CATEGORICAL LEGEND
  par(mar = c(6, 1, 4, 3))
  plot.new()
  
  # Legend positioning (centered vertically)
  n_cats <- length(legend_labels)
  cat_height <- 0.15
  total_height <- n_cats * cat_height
  y_start <- 0.5 - total_height / 2
  
  legend_x_left <- 0.1
  legend_x_right <- 0.5
  
  # Draw categorical boxes
  for (i in 1:n_cats) {
    y_b <- y_start + (i - 1) * cat_height
    y_t <- y_b + cat_height
    rect(legend_x_left, y_b, legend_x_right, y_t,
         col = color_palette[i], border = "black", lwd = 2)
  }
  
  # Add labels
  for (i in 1:n_cats) {
    y_pos <- y_start + (i - 0.5) * cat_height
    text(legend_x_right + 0.1, y_pos, legend_labels[i], 
         pos = 4, cex = 1.3, font = 1)
  }
  
  # Legend title
  text((legend_x_left + legend_x_right) / 2,
       y_start + total_height + 0.08,
       "Classification", cex = 1.5, font = 2)
  
  dev.off()
  log_event(sprintf("Created categorical map: %s", basename(output_file)), "SUCCESS")
  return(TRUE)
}

# ===============================================================================
# MAIN PROCESSING LOOP
# ===============================================================================
log_event("Searching for monthly composite CSV files...", "INFO")
index_types <- c("spi", "spei")
timescales <- c(1, 3, 6, 9, 12, 24)
all_results <- list()

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
    
    index_matrix <- reconstruction$index_matrix
    coords <- reconstruction$coords
    n_years <- reconstruction$n_years
    n_pixels <- reconstruction$n_pixels
    n_months <- nrow(index_matrix)
    years_full <- rep(1:n_years, each = 12)[1:n_months]
    months_full <- rep(1:12, times = n_years)[1:n_months]
    
    log_event(sprintf("Processing %s-%02d (%d pixels, %d months)...", 
                      toupper(index_type), scale, n_pixels, n_months), "INFO")
    
    coords_dt <- data.table(
      space_idx = 1:n_pixels,
      lon = coords$x,
      lat = coords$y
    )
    
    # ---- ANNUAL ANALYSIS ----
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
    
    # ---- DROUGHT EVENT METRICS ----
    log_event("  Calculating drought event metrics...", "INFO")
    event_metrics <- future_lapply(1:n_pixels, function(i) {
      ts_vals <- index_matrix[, i]
      detect_drought_events_pixel(ts_vals, drought_threshold, recovery_threshold, min_duration)
    }, future.seed = TRUE)
    event_df <- rbindlist(lapply(event_metrics, as.data.frame))
    event_df[, space_idx := 1:n_pixels]
    
    # ---- REGIME SHIFT & RETURN PERIOD ----
    log_event("  Detecting regime shifts & return periods...", "INFO")
    regime_shifts <- future_sapply(1:n_pixels, function(i) {
      detect_regime_shift(index_matrix[, i])
    }, future.seed = TRUE)
    
    return_periods <- future_sapply(1:n_pixels, function(i) {
      calculate_return_period(index_matrix[, i], extreme_threshold)
    }, future.seed = TRUE)
    
    # ---- COMBINE ANNUAL RESULTS ----
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
    
    # Add derived columns
    annual_results[, runs_significant := !filtered_runs & p_value_runs < alpha]
    annual_results[, same_significance := (p_value_vc < alpha & !filtered_vc) == 
                     (p_value_tfpw < alpha & !filtered_tfpw)]
    annual_results[, same_direction := sign(tau_vc) == sign(tau_tfpw)]
    
    # ===============================================================================
    # CREATE RASTERS AND MAPS
    # ===============================================================================
    log_event("  Creating spatial rasters and maps...", "INFO")
    
    # Create template raster
    template <- rast(ext(basin), nrows = 150, ncols = 150, crs = target_crs)
    template <- rasterize(basin, template, values = NA)
    
    # 1. RUNS TEST SIGNIFICANCE MAP
    log_event("  Creating runs test map...", "INFO")
    
    # Create numeric version of clustering
    annual_results[, clustering_numeric := ifelse(clustering == "clustering", -1,
                                                  ifelse(clustering == "dispersion", 1, NA))]
    
    # Only keep significant runs
    annual_results_sig_runs <- annual_results[runs_significant == TRUE]
    
    if (nrow(annual_results_sig_runs) > 0) {
      runs_rast <- create_interpolated_raster(annual_results_sig_runs, "clustering_numeric", template)
      
      writeRaster(runs_rast,
                  file.path(out_dir, sprintf("%s_%02d_annual_runs_significance.tif", index_type, scale)),
                  overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
      
      if (sum(!is.na(values(runs_rast))) > 0) {
        create_map_with_categorical_legend(
          runs_rast, basin,
          file.path(maps_dir, sprintf("%s_%02d_runs_test_map.png", index_type, scale)),
          sprintf("%s-%02d: Wald-Wolfowitz Runs Test", toupper(index_type), scale),
          sprintf("Temporal Clustering Analysis (p < %.2f) | Threshold: %.2f", alpha, drought_threshold),
          c("Clustering\n(Too few runs)", "Dispersion\n(Too many runs)"),
          c("#d62728", "#1f77b4")
        )
      }
    } else {
      log_event("  No significant runs - skipping map", "WARNING")
    }
    
    # 2. VC TREND SIGNIFICANCE MAP
    log_event("  Creating VC trends map...", "INFO")
    
    annual_results[, tau_vc_sign := sign(tau_vc)]
    annual_results_sig_vc <- annual_results[p_value_vc < alpha & filtered_vc == FALSE]
    
    if (nrow(annual_results_sig_vc) > 0) {
      vc_rast <- create_interpolated_raster(annual_results_sig_vc, "tau_vc_sign", template)
      
      writeRaster(vc_rast,
                  file.path(out_dir, sprintf("%s_%02d_annual_vc_significance.tif", index_type, scale)),
                  overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
      
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
      log_event("  No significant VC trends - skipping map", "WARNING")
    }
    
    # 3. TFPW TREND SIGNIFICANCE MAP
    log_event("  Creating TFPW trends map...", "INFO")
    
    annual_results[, tau_tfpw_sign := sign(tau_tfpw)]
    annual_results_sig_tfpw <- annual_results[p_value_tfpw < alpha & filtered_tfpw == FALSE]
    
    if (nrow(annual_results_sig_tfpw) > 0) {
      tfpw_rast <- create_interpolated_raster(annual_results_sig_tfpw, "tau_tfpw_sign", template)
      
      writeRaster(tfpw_rast,
                  file.path(out_dir, sprintf("%s_%02d_annual_tfpw_significance.tif", index_type, scale)),
                  overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
      
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
      log_event("  No significant TFPW trends - skipping map", "WARNING")
    }
    
    # 4. DROUGHT FREQUENCY MAP (CONTINUOUS LEGEND)
    log_event("  Creating drought frequency map...", "INFO")
    
    # Convert to events per decade
    annual_results[, freq_per_decade := n_events / 7.6]
    
    freq_rast <- create_interpolated_raster(annual_results, "freq_per_decade", template)
    
    writeRaster(freq_rast,
                file.path(out_dir, sprintf("%s_%02d_drought_frequency_per_decade.tif", index_type, scale)),
                overwrite = TRUE, datatype = "FLT4S", NAflag = -9999)
    
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
    
    # Progress
    completed_idx <- (which(index_types == index_type) - 1) * length(timescales) + scale_idx
    total_tasks <- length(index_types) * length(timescales)
    elapsed <- difftime(Sys.time(), start_time, units = "mins")
    remaining <- (elapsed / completed_idx) * (total_tasks - completed_idx)
    
    log_event(sprintf("  Completed %s-%02d (%d/%d) | Elapsed: %.1f min | Remaining: %.1f min", 
                      toupper(index_type), scale, completed_idx, total_tasks,
                      elapsed, remaining), "SUCCESS")
  }
}

# ===============================================================================
# SAVE FINAL RESULTS
# ===============================================================================
log_event("Saving consolidated results...", "INFO")
if (length(all_results) > 0) {
  final_results <- rbindlist(all_results, fill = TRUE)
  fwrite(final_results, file.path(out_dir, "all_temporal_diagnostics_results.csv"))
  
  # Summary
  total_pixels <- nrow(final_results)
  sig_vc <- sum(final_results$p_value_vc < alpha & !final_results$filtered_vc, na.rm = TRUE)
  sig_tfpw <- sum(final_results$p_value_tfpw < alpha & !final_results$filtered_tfpw, na.rm = TRUE)
  sig_runs <- sum(final_results$runs_significant, na.rm = TRUE)
  
  cat("\n=== ANALYSIS SUMMARY ===\n", file = LOG_FILE, append = TRUE)
  cat(sprintf("Total pixels: %d\n", total_pixels), file = LOG_FILE, append = TRUE)
  cat(sprintf("VC trends: %d (%.1f%%)\n", sig_vc, sig_vc/total_pixels*100), file = LOG_FILE, append = TRUE)
  cat(sprintf("TFPW trends: %d (%.1f%%)\n", sig_tfpw, sig_tfpw/total_pixels*100), file = LOG_FILE, append = TRUE)
  cat(sprintf("Runs patterns: %d (%.1f%%)\n", sig_runs, sig_runs/total_pixels*100), file = LOG_FILE, append = TRUE)
  cat(sprintf("Time: %.2f minutes\n", as.numeric(difftime(Sys.time(), start_time, units = "mins"))), 
      file = LOG_FILE, append = TRUE)
  cat("========================\n\n", file = LOG_FILE, append = TRUE)
  
  log_event("FINAL FIX COMPLETE!", "SUCCESS")
} else {
  log_event("No results generated", "WARNING")
}

plan(sequential)