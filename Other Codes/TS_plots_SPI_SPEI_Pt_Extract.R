# ============================================================================
# SCRIPT: Create Time Series Plots with Trend Analysis 
# ============================================================================
# This script creates 4-panel plots showing SPI and SPEI at 1, 3, 6, and 12 month timescales
# PLUS trend analysis with explicit handling and visualization of data gaps

# Load required libraries
library(terra)       # For raster processing
library(ggplot2)     # For plotting
library(dplyr)       # For data manipulation
library(tidyr)       # For data reshaping
library(patchwork)   # For combining plots
library(lubridate)   # For date handling
library(ncdf4)       # For NetCDF handling
library(trend)       # For Mann-Kendall and Sen's slope tests
library(Kendall)     # For Mann-Kendall trend test

# ===================== CRITICAL UNDERSTANDING: DATA GAPS IN YOUR PLOTS =====================
#
# YOU ARE SEEING MISSING DATA AFTER 2020 IN SPI-12 AND SPI-6 PLOTS
#
# This is NOT normal behavior and suggests one of these issues:
#
# 1. **REFERENCE PERIOD LIMITATION**:
# - Many SPI implementations use a fixed reference period (e.g., 1981-2010)
# - Data OUTSIDE this reference period may not be properly standardized
# - Check: What reference period was used in your SPI calculation?
#
# 2. **INCOMPLETE INPUT DATA**:
# - Your precipitation or PET data may actually end around 2020
# - Even if NetCDF claims to have 2024 data, actual values might be NA
# - Check: Open your NetCDF files and look at the actual data values for recent years
#
# 3. **DISTRIBUTION FITTING FAILURE**:
# - SPI requires fitting a probability distribution to the data
# - If insufficient data exists for recent periods, fitting fails → NA values
# - This happens if using a moving/rolling reference window
#
# 4. **SOFTWARE-SPECIFIC BEHAVIOR**:
# - What software/package calculated your SPI? (CDO, Python SPEI, R package?)
# - Some implementations have known issues with edge effects
#
# DIAGNOSTIC STEPS YOU SHOULD TAKE:
#
# A. Check your NetCDF files directly:
# ```R
# nc <- nc_open("spi_results/spi_12month.nc")
# print(nc)
# time_vals <- ncvar_get(nc, "time")
# spi_vals <- ncvar_get(nc, "spi") # or whatever the variable name is
# # Check if last ~50 time steps have actual data or all NAs
# ```
#
# B. Check your INPUT precipitation data:
# - Does it actually extend to 2024?
# - Are there gaps in recent years?
#
# C. Review your SPI calculation code/parameters:
# - What reference period was specified?
# - What distribution fitting method was used?
#
# THIS SCRIPT NOW INCLUDES:
# - Explicit gap visualization in ALL plots (including trend plots)
# - Detailed diagnostic messages about data availability
# - Warnings when recent data is systematically missing
# ============================================================================

# ===================== CONFIGURATION =====================
setwd("D:/Nechako_Drought/")

spi_data_dir <- "spi_results/"
spei_data_dir <- "spei_results/"
basin_shp_path <- "Spatial/nechakoBound_dissolve.shp"
output_dir <- "point_timeseries_plots_with_trends_and_diagnostics/"
dir.create(output_dir, showWarnings = FALSE)

spi_files <- c(
  "spi_01month.nc",
  "spi_03month.nc",
  "spi_06month.nc",
  "spi_12month.nc"
)

spei_files <- c(
  "spei_01month.nc",
  "spei_03month.nc",
  "spei_06month.nc",
  "spei_12month.nc"
)

use_specific_points <- FALSE
specific_points <- data.frame(
  x = c(-123.5, -124.2, -123.8),
  y = c(54.2, 53.8, 54.5)
)

n_random_points <- 5

# ===================== TREND ANALYSIS CONFIGURATION =====================
alpha <- 0.05
trend_start_year <- NULL
trend_end_year <- NULL

# =============================================================================
# NEW: Helper to prevent fake "gaps" caused by ylim() dropping values < -4 or > 4
# =============================================================================
# coord_cartesian() zooms without discarding data (unlike ylim()).
# This function expands beyond [-4,4] only when needed, with padding.
calc_dynamic_ylim <- function(v, base = c(-4, 4), pad_frac = 0.10) {
  v <- v[is.finite(v)]
  if (length(v) == 0) return(base)
  
  r <- range(v)
  span <- diff(r)
  if (!is.finite(span) || span == 0) span <- 1
  
  y_min <- min(base[1], r[1] - pad_frac * span)
  y_max <- max(base[2], r[2] + pad_frac * span)
  
  c(y_min, y_max)
}

# ===================== FUNCTION: EXTRACT DATES FROM NETCDF =====================
extract_dates_from_nc <- function(nc_file_path) {
  cat(" Extracting dates from:", basename(nc_file_path), "\n")
  r_temp <- rast(nc_file_path)
  n_layers <- nlyr(r_temp)
  cat(" Number of layers in file:", n_layers, "\n")
  
  tryCatch({
    dates <- time(r_temp)
    if (!is.null(dates) && length(dates) > 0) {
      cat(" Using terra time() method\n")
      dates_converted <- tryCatch({
        as.Date(dates)
      }, error = function(e) {
        as.Date(dates, origin = "1970-01-01")
      })
      
      if (all(!is.na(dates_converted)) &&
          all(dates_converted >= as.Date("1900-01-01")) &&
          all(dates_converted <= as.Date("2025-12-31"))) {
        
        if (length(dates_converted) > n_layers) {
          dates_converted <- dates_converted[1:n_layers]
          cat(" Trimmed dates to match layer count\n")
        } else if (length(dates_converted) < n_layers) {
          cat(" WARNING: Fewer dates than layers. Generating missing dates.\n")
          last_date <- tail(dates_converted, 1)
          missing_dates <- seq(last_date + months(1), by = "month", length.out = n_layers - length(dates_converted))
          dates_converted <- c(dates_converted, missing_dates)
        }
        
        cat(" Date range:", format(min(dates_converted), "%Y-%m-%d"), "to",
            format(max(dates_converted), "%Y-%m-%d"), "\n")
        return(dates_converted)
      }
    }
    
    cat(" Trying ncdf4 method\n")
    nc <- nc_open(nc_file_path)
    
    if ("time" %in% names(nc$var)) {
      time_vals <- ncvar_get(nc, "time")
      time_units <- ncatt_get(nc, "time", "units")$value
      
      if (grepl("since", time_units)) {
        origin_str <- sub(".*since ", "", time_units)
        origin_date <- as.Date(origin_str)
        
        if (grepl("^days", time_units)) {
          dates <- origin_date + time_vals
        } else if (grepl("^months", time_units)) {
          dates <- seq(origin_date, by = "month", length.out = min(length(time_vals), n_layers))
        }
        
        nc_close(nc)
        
        if (all(dates >= as.Date("1900-01-01")) &&
            all(dates <= as.Date("2025-12-31"))) {
          cat(" Successfully extracted dates using ncdf4\n")
          if (length(dates) > n_layers) {
            dates <- dates[1:n_layers]
          } else if (length(dates) < n_layers) {
            cat(" WARNING: Fewer dates than layers. Generating missing dates.\n")
            last_date <- tail(dates, 1)
            missing_dates <- seq(last_date + months(1), by = "month", length.out = n_layers - length(dates))
            dates <- c(dates, missing_dates)
          }
          
          cat(" Date range:", format(min(dates), "%Y-%m-%d"), "to",
              format(max(dates), "%Y-%m-%d"), "\n")
          return(as.Date(dates))
        }
      }
    }
    
    nc_close(nc)
  }, error = function(e) {
    cat(" Date extraction failed:", e$message, "\n")
  })
  
  cat(" Using fallback date generation\n")
  start_date <- as.Date("1954-01-01")
  dates <- seq(start_date, by = "month", length.out = n_layers)
  
  end_date <- dates[length(dates)]
  if (end_date > as.Date("2024-12-31")) {
    cat(" WARNING: Generated dates exceed Dec 2024. Adjusting...\n")
    end_target <- as.Date("2024-12-01")
    start_date <- end_target - months(n_layers - 1)
    dates <- seq(start_date, by = "month", length.out = n_layers)
  }
  
  cat(" Generated dates from", format(dates[1], "%Y-%m-%d"),
      "to", format(dates[length(dates)], "%Y-%m-%d"), "\n")
  return(dates)
}

# ===================== FUNCTION: GET GRID POINTS INSIDE BASIN =====================
get_grid_points_inside_basin <- function(basin_shp, raster_template) {
  cat("Creating raster mask of basin...\n")
  mask_raster <- rast(raster_template)
  values(mask_raster) <- 1
  mask_raster <- mask(mask_raster, basin_shp)
  
  mask_df <- as.data.frame(mask_raster, xy = TRUE, na.rm = TRUE)
  if (ncol(mask_df) >= 2) {
    colnames(mask_df)[1:2] <- c("x", "y")
  } else {
    stop("Error: Could not extract coordinates from raster mask")
  }
  
  valid_points <- mask_df[, c("x", "y")]
  cat("Found", nrow(valid_points), "grid points inside basin\n")
  return(valid_points)
}

# ===================== FUNCTION: CHECK IF POINTS ARE INSIDE BASIN =====================
check_points_inside_basin <- function(points_df, basin_shp, raster_template) {
  pts <- vect(points_df[, c("x", "y")], geom = c("x", "y"), crs = crs(raster_template))
  inside <- is.inside(pts, basin_shp)
  
  valid_points <- points_df[inside, ]
  if (nrow(valid_points) == 0) {
    stop("Error: None of the specified points are inside the basin!")
  } else if (nrow(valid_points) < nrow(points_df)) {
    cat("Warning:", nrow(points_df) - nrow(valid_points), "points were outside the basin and removed\n")
  }
  
  return(valid_points)
}

# ===================== FUNCTION: GET SELECTED POINTS =====================
get_selected_points <- function(basin_shp, raster_template, use_specific = FALSE,
                                specific_pts = NULL, n_random = 5) {
  if (use_specific) {
    if (is.null(specific_pts)) {
      stop("Error: use_specific is TRUE but no specific_pts provided")
    }
    cat("Using", nrow(specific_pts), "specific point(s):\n")
    print(specific_pts)
    selected_points <- check_points_inside_basin(specific_pts, basin_shp, raster_template)
    cat("Valid points inside basin:", nrow(selected_points), "\n")
  } else {
    cat("Selecting", n_random, "random points from basin grid...\n")
    all_points <- get_grid_points_inside_basin(basin_shp, raster_template)
    
    if (nrow(all_points) < n_random) {
      cat("Warning: Requested", n_random, "points but only", nrow(all_points), "available. Using all.\n")
      selected_points <- all_points
    } else {
      set.seed(42)
      sample_indices <- sample(1:nrow(all_points), n_random)
      selected_points <- all_points[sample_indices, ]
    }
  }
  return(selected_points)
}

# ===================== FUNCTION: EXTRACT TIME SERIES AT POINT =====================
extract_timeseries_at_point <- function(raster_stack, x, y, index_name, nc_file_path) {
  cat(" Extracting", index_name, "...\n")
  point_vect <- vect(data.frame(x = x, y = y), geom = c("x", "y"), crs = crs(raster_stack))
  
  # CRITICAL: Use terra::extract explicitly to avoid namespace conflict with dplyr
  values <- terra::extract(raster_stack, point_vect)
  values_numeric <- as.numeric(values[1, -1])
  
  dates <- extract_dates_from_nc(nc_file_path)
  
  if (length(dates) != length(values_numeric)) {
    cat(" WARNING: Date and value length mismatch for", index_name, "\n")
    cat(" Dates:", length(dates), "Values:", length(values_numeric), "\n")
    min_length <- min(length(dates), length(values_numeric))
    dates <- dates[1:min_length]
    values_numeric <- values_numeric[1:min_length]
  }
  
  data.frame(
    Date = dates,
    Value = values_numeric,
    Index = index_name,
    stringsAsFactors = FALSE
  )
}

# ===================== NEW FUNCTION: DIAGNOSE DATA GAPS =====================
diagnose_data_gaps <- function(ts_data, index_name) {
  ts_with_na <- ts_data %>%
    mutate(
      Year = year(Date),
      is_missing = is.na(Value)
    )
  
  n_total <- nrow(ts_with_na)
  n_missing <- sum(ts_with_na$is_missing)
  pct_missing <- round(100 * n_missing / n_total, 1)
  
  recent_years <- ts_with_na %>% filter(Year >= 2020)
  if (nrow(recent_years) > 0) {
    n_recent_missing <- sum(recent_years$is_missing)
    pct_recent_missing <- round(100 * n_recent_missing / nrow(recent_years), 1)
    if (pct_recent_missing > 50) {
      cat(" ⚠️ WARNING: ", index_name, " has ", pct_recent_missing,
          "% missing data after 2020!\n", sep = "")
      cat(" This suggests issues with reference period or input data\n")
    }
  }
  
  valid_data <- ts_with_na %>% filter(!is_missing)
  if (nrow(valid_data) > 0) {
    last_valid <- max(valid_data$Date)
    cat(" Last valid value:", format(last_valid, "%Y-%m-%d"), "\n")
    
    data_end <- max(ts_with_na$Date)
    months_gap <- interval(last_valid, data_end) %/% months(1)
    if (months_gap > 6) {
      cat(" ⚠️ GAP: ", months_gap, " months of missing data at end of series\n", sep = "")
    }
  }
  
  return(list(
    n_total = n_total,
    n_missing = n_missing,
    pct_missing = pct_missing,
    last_valid = if (nrow(valid_data) > 0) last_valid else NA
  ))
}

# ===================== NEW FUNCTION: CALCULATE TREND STATISTICS WITH GAP HANDLING =====================
calculate_trend_statistics <- function(ts_data, index_name, start_year = NULL, end_year = NULL) {
  if (!is.null(start_year)) {
    ts_data <- ts_data %>% filter(year(Date) >= start_year)
  }
  if (!is.null(end_year)) {
    ts_data <- ts_data %>% filter(year(Date) <= end_year)
  }
  
  ts_original <- ts_data
  ts_clean <- ts_data %>% filter(!is.na(Value))
  
  if (nrow(ts_clean) < 10) {
    cat(" Insufficient data for trend analysis of", index_name, "\n")
    return(NULL)
  }
  
  completeness <- nrow(ts_clean) / nrow(ts_original) * 100
  if (completeness < 70) {
    cat(" ⚠️ WARNING:", index_name, "only has", round(completeness, 1),
        "% complete data. Trend may be unreliable.\n")
  }
  
  ts_clean$time_index <- 1:nrow(ts_clean)
  
  lm_model <- lm(Value ~ time_index, data = ts_clean)
  lm_slope <- coef(lm_model)[2]
  lm_pvalue <- summary(lm_model)$coefficients[2, 4]
  lm_rsquared <- summary(lm_model)$r.squared
  
  sens_result <- tryCatch({
    sens.slope(ts_clean$Value)
  }, error = function(e) {
    cat(" Error calculating Sen's slope for", index_name, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(sens_result)) {
    return(NULL)
  }
  
  sens_slope <- sens_result$estimates
  sens_pvalue <- sens_result$p.value
  
  mk_result <- tryCatch({
    MannKendall(ts_clean$Value)
  }, error = function(e) {
    cat(" Error in Mann-Kendall test for", index_name, ":", e$message, "\n")
    return(NULL)
  })
  
  mk_tau <- if (!is.null(mk_result)) mk_result$tau else NA
  mk_pvalue <- if (!is.null(mk_result)) mk_result$sl else NA
  
  results <- data.frame(
    Index = index_name,
    N_obs = nrow(ts_clean),
    N_total = nrow(ts_original),
    Completeness_pct = round(completeness, 1),
    Period_start = format(min(ts_clean$Date), "%Y-%m"),
    Period_end = format(max(ts_clean$Date), "%Y-%m"),
    LR_slope = lm_slope,
    LR_pvalue = lm_pvalue,
    LR_significant = lm_pvalue < alpha,
    LR_rsquared = lm_rsquared,
    Sens_slope = sens_slope,
    Sens_pvalue = sens_pvalue,
    Sens_significant = sens_pvalue < alpha,
    MK_tau = mk_tau,
    MK_pvalue = mk_pvalue,
    MK_significant = mk_pvalue < alpha,
    stringsAsFactors = FALSE
  )
  
  ts_clean$LR_fitted <- predict(lm_model)
  ts_clean$Sens_fitted <- median(ts_clean$Value) +
    sens_slope * (ts_clean$time_index - median(ts_clean$time_index))
  
  return(list(
    statistics = results,
    data_clean = ts_clean,
    data_original = ts_original
  ))
}

# ===================== NEW FUNCTION: CREATE TREND PLOT WITH GAPS SHOWN =====================
create_trend_plot <- function(trend_result, index_name) {
  if (is.null(trend_result)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste(index_name, "- Insufficient data"),
                      size = 4, color = "red") +
             theme_void())
  }
  
  ts_original <- trend_result$data_original
  ts_clean <- trend_result$data_clean
  stats <- trend_result$statistics
  
  completeness_note <- if (stats$Completeness_pct < 80)
    sprintf("\n⚠️ Only %.0f%% complete", stats$Completeness_pct) else ""
  
  lr_label <- sprintf("Linear: slope=%.4f (p=%.3f)%s",
                      stats$LR_slope, stats$LR_pvalue,
                      ifelse(stats$LR_significant, "*", ""))
  
  sens_label <- sprintf("Sen's: slope=%.4f (p=%.3f)%s",
                        stats$Sens_slope, stats$Sens_pvalue,
                        ifelse(stats$Sens_significant, "*", ""))
  
  mk_label <- sprintf("M-K: τ=%.3f (p=%.3f)%s%s",
                      stats$MK_tau, stats$MK_pvalue,
                      ifelse(stats$MK_significant, "*", ""),
                      completeness_note)
  
  p <- ggplot(ts_original, aes(x = Date, y = Value)) +
    geom_line(color = "blue", alpha = 0.6, na.rm = FALSE) +
    geom_point(color = "blue", size = 0.5, alpha = 0.4, na.rm = TRUE) +
    geom_line(data = ts_clean, aes(y = LR_fitted), color = "red", linewidth = 1, alpha = 0.7) +
    geom_line(data = ts_clean, aes(y = Sens_fitted), color = "darkgreen", linewidth = 1,
              linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    annotate("text",
             x = min(ts_original$Date, na.rm = TRUE),
             y = max(ts_original$Value, na.rm = TRUE),
             label = paste(lr_label, sens_label, mk_label, sep = "\n"),
             hjust = 0, vjust = 1, size = 2.5, color = "black",
             fontface = "plain") +
    labs(
      title = index_name,
      x = "Date",
      y = "Index Value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# ===================== FUNCTION: CREATE 4-PANEL PLOT FOR 1 & 12 MONTHS =====================
create_four_panel_plot_1_12 <- function(point_data, point_coords, point_id, output_file) {
  expected_indices <- c("SPI1", "SPI12", "SPEI1", "SPEI12")
  indices <- unique(point_data$Index)
  plots <- list()
  
  for (index in expected_indices) {
    if (index %in% indices) {
      index_data <- point_data %>% filter(Index == index)
      index_data$Value[is.infinite(index_data$Value)] <- NA
      
      # NEW: Dynamic y-limits that DO NOT drop data (prevents fake gaps)
      ylims <- calc_dynamic_ylim(index_data$Value)
      
      p <- ggplot(index_data, aes(x = Date, y = Value)) +
        geom_line(color = "blue", linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2,
                 fill = "#8B0000", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = -1.5,
                 fill = "#FF0000", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.5, ymax = -1,
                 fill = "#FF8C00", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = -0.5,
                 fill = "#FFA500", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1,
                 fill = "#90EE90", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1.5,
                 fill = "#00FF00", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2,
                 fill = "#008000", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf,
                 fill = "#006400", alpha = 0.1) +
        geom_hline(yintercept = c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2),
                   linetype = "dotted", color = "gray50", linewidth = 0.3) +
        labs(
          title = index,
          x = ifelse(index == "SPEI12", "Date", ""),
          y = "Index Value"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          panel.grid.minor = element_blank(),
          panel.spacing = grid::unit(0.5, "lines")
        ) +
        # FIX: coord_cartesian zooms without dropping data (unlike ylim)
        coord_cartesian(ylim = ylims)
      
      plots[[index]] <- p
    } else {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste(index, "data not available"),
                 size = 5, color = "red") +
        theme_void()
      plots[[index]] <- p
    }
  }
  
  combined_plot <- (plots[["SPI1"]] | plots[["SPI12"]]) /
    (plots[["SPEI1"]] | plots[["SPEI12"]]) +
    plot_annotation(
      title = sprintf("Point #%d: (%.4f°N, %.4f°W) - 1 & 12 Month Timescales",
                      point_id, point_coords$y, abs(point_coords$x)),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10)
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(output_file, combined_plot, width = 12, height = 10, dpi = 150)
  cat(" Saved plot to:", output_file, "\n")
}

# ===================== FUNCTION: CREATE 4-PANEL PLOT FOR 3 & 6 MONTHS =====================
create_four_panel_plot_3_6 <- function(point_data, point_coords, point_id, output_file) {
  expected_indices <- c("SPI3", "SPI6", "SPEI3", "SPEI6")
  indices <- unique(point_data$Index)
  plots <- list()
  
  for (index in expected_indices) {
    if (index %in% indices) {
      index_data <- point_data %>% filter(Index == index)
      index_data$Value[is.infinite(index_data$Value)] <- NA
      
      # NEW: Dynamic y-limits that DO NOT drop data (prevents fake gaps)
      ylims <- calc_dynamic_ylim(index_data$Value)
      
      p <- ggplot(index_data, aes(x = Date, y = Value)) +
        geom_line(color = "blue", linewidth = 1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2,
                 fill = "#8B0000", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2, ymax = -1.5,
                 fill = "#FF0000", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.5, ymax = -1,
                 fill = "#FF8C00", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = -0.5,
                 fill = "#FFA500", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1,
                 fill = "#90EE90", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1.5,
                 fill = "#00FF00", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2,
                 fill = "#008000", alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf,
                 fill = "#006400", alpha = 0.1) +
        geom_hline(yintercept = c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2),
                   linetype = "dotted", color = "gray50", linewidth = 0.3) +
        labs(
          title = index,
          x = ifelse(index == "SPEI6", "Date", ""),
          y = "Index Value"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          panel.grid.minor = element_blank(),
          panel.spacing = grid::unit(0.5, "lines")
        ) +
        # FIX: coord_cartesian zooms without dropping data (unlike ylim)
        coord_cartesian(ylim = ylims)
      
      plots[[index]] <- p
    } else {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = paste(index, "data not available"),
                 size = 5, color = "red") +
        theme_void()
      plots[[index]] <- p
    }
  }
  
  combined_plot <- (plots[["SPI3"]] | plots[["SPI6"]]) /
    (plots[["SPEI3"]] | plots[["SPEI6"]]) +
    plot_annotation(
      title = sprintf("Point #%d: (%.4f°N, %.4f°W) - 3 & 6 Month Timescales",
                      point_id, point_coords$y, abs(point_coords$x)),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10)
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(output_file, combined_plot, width = 12, height = 10, dpi = 150)
  cat(" Saved plot to:", output_file, "\n")
}

# ===================== NEW FUNCTION: CREATE TREND ANALYSIS PLOT =====================
create_trend_analysis_plot <- function(all_trends, point_coords, point_id, output_file) {
  plots <- list()
  indices <- c("SPI3", "SPI6", "SPEI3", "SPEI6", "SPI1", "SPI12", "SPEI1", "SPEI12")
  
  for (index in indices) {
    if (index %in% names(all_trends)) {
      plots[[index]] <- create_trend_plot(all_trends[[index]], index)
    } else {
      plots[[index]] <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste(index, "- No data"),
                 size = 4, color = "red") +
        theme_void()
    }
  }
  
  combined_plot <- (plots[["SPI1"]] | plots[["SPI3"]] | plots[["SPI6"]] | plots[["SPI12"]]) /
    (plots[["SPEI1"]] | plots[["SPEI3"]] | plots[["SPEI6"]] | plots[["SPEI12"]]) +
    plot_annotation(
      title = sprintf("Trend Analysis - Point #%d: (%.4f°N, %.4f°W)\nRed=Linear Regression, Green=Sen's Slope, *=Significant (α=%.2f)",
                      point_id, point_coords$y, abs(point_coords$x), alpha),
      subtitle = "Gaps in data are shown. Trend lines fitted to available data only.",
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "darkred"),
        plot.margin = margin(10, 10, 10, 10)
      )
    )
  
  ggsave(output_file, combined_plot, width = 16, height = 8, dpi = 150)
  cat(" Saved trend plot to:", output_file, "\n")
}

# ===================== MAIN PROCESSING =====================
cat("Loading basin shapefile...\n")
basin_shp <- vect(basin_shp_path)

cat("Loading template raster...\n")
template_file <- file.path(spi_data_dir, spi_files[1])
template_raster <- rast(template_file)

if (!same.crs(basin_shp, template_raster)) {
  cat("Reprojecting basin to match raster CRS...\n")
  basin_shp <- project(basin_shp, crs(template_raster))
}

selected_points <- get_selected_points(
  basin_shp,
  template_raster,
  use_specific = use_specific_points,
  specific_pts = specific_points,
  n_random = n_random_points
)

selected_points$point_id <- 1:nrow(selected_points)
all_trend_stats <- list()

# Process each selected point
for (i in 1:nrow(selected_points)) {
  point <- selected_points[i, ]
  cat("\n========================================\n")
  cat("Processing point", i, "of", nrow(selected_points), ":",
      sprintf("(%.4f, %.4f)", point$x, point$y), "\n")
  cat("========================================\n")
  
  cat(" Loading raster data...\n")
  spi1_file <- file.path(spi_data_dir, spi_files[1])
  spi3_file <- file.path(spi_data_dir, spi_files[2])
  spi6_file <- file.path(spi_data_dir, spi_files[3])
  spi12_file <- file.path(spi_data_dir, spi_files[4])
  
  spi1_stack <- rast(spi1_file)
  spi3_stack <- rast(spi3_file)
  spi6_stack <- rast(spi6_file)
  spi12_stack <- rast(spi12_file)
  
  spei1_file <- file.path(spei_data_dir, spei_files[1])
  spei3_file <- file.path(spei_data_dir, spei_files[2])
  spei6_file <- file.path(spei_data_dir, spei_files[3])
  spei12_file <- file.path(spei_data_dir, spei_files[4])
  
  spei1_stack <- rast(spei1_file)
  spei3_stack <- rast(spei3_file)
  spei6_stack <- rast(spei6_file)
  spei12_stack <- rast(spei12_file)
  
  cat(" Extracting time series...\n")
  spi1_ts <- extract_timeseries_at_point(spi1_stack, point$x, point$y, "SPI1", spi1_file)
  spi3_ts <- extract_timeseries_at_point(spi3_stack, point$x, point$y, "SPI3", spi3_file)
  spi6_ts <- extract_timeseries_at_point(spi6_stack, point$x, point$y, "SPI6", spi6_file)
  spi12_ts <- extract_timeseries_at_point(spi12_stack, point$x, point$y, "SPI12", spi12_file)
  
  spei1_ts <- extract_timeseries_at_point(spei1_stack, point$x, point$y, "SPEI1", spei1_file)
  spei3_ts <- extract_timeseries_at_point(spei3_stack, point$x, point$y, "SPEI3", spei3_file)
  spei6_ts <- extract_timeseries_at_point(spei6_stack, point$x, point$y, "SPEI6", spei6_file)
  spei12_ts <- extract_timeseries_at_point(spei12_stack, point$x, point$y, "SPEI12", spei12_file)
  
  # ==================== DATA GAP DIAGNOSTICS ====================
  cat("\n DATA GAP DIAGNOSTICS:\n")
  cat(" =====================\n")
  for (idx_name in c("SPI1", "SPI3", "SPI6", "SPI12", "SPEI1", "SPEI3", "SPEI6", "SPEI12")) {
    ts_data <- switch(idx_name,
                      "SPI1" = spi1_ts, "SPI3" = spi3_ts, "SPI6" = spi6_ts, "SPI12" = spi12_ts,
                      "SPEI1" = spei1_ts, "SPEI3" = spei3_ts, "SPEI6" = spei6_ts, "SPEI12" = spei12_ts)
    cat(" ", idx_name, ":\n", sep = "")
    diagnose_data_gaps(ts_data, idx_name)
  }
  
  # Combine time series
  all_ts_1_12 <- bind_rows(spi1_ts, spi12_ts, spei1_ts, spei12_ts)
  all_ts_3_6  <- bind_rows(spi3_ts, spi6_ts, spei3_ts, spei6_ts)
  
  # ==================== TREND ANALYSIS ====================
  cat("\n Performing trend analysis...\n")
  all_trends <- list()
  point_trend_stats <- list()
  
  for (idx in c("SPI1", "SPI3", "SPI6", "SPI12", "SPEI1", "SPEI3", "SPEI6", "SPEI12")) {
    ts_data <- switch(idx,
                      "SPI1" = spi1_ts, "SPI3" = spi3_ts, "SPI6" = spi6_ts, "SPI12" = spi12_ts,
                      "SPEI1" = spei1_ts, "SPEI3" = spei3_ts, "SPEI6" = spei6_ts, "SPEI12" = spei12_ts)
    
    cat(" Analyzing trends for", idx, "...\n")
    trend_result <- calculate_trend_statistics(ts_data, idx, trend_start_year, trend_end_year)
    
    if (!is.null(trend_result)) {
      all_trends[[idx]] <- trend_result
      point_trend_stats[[idx]] <- trend_result$statistics
    }
  }
  
  if (length(point_trend_stats) > 0) {
    point_stats_df <- bind_rows(point_trend_stats) %>%
      mutate(Point_ID = point$point_id,
             Point_X = point$x,
             Point_Y = point$y)
    all_trend_stats[[i]] <- point_stats_df
  }
  
  # Create output filenames
  output_file_1_12 <- file.path(output_dir,
                                sprintf("point_%d_%.4f_%.4f_timeseries_1_12months.png",
                                        point$point_id, point$y, abs(point$x)))
  output_file_3_6 <- file.path(output_dir,
                               sprintf("point_%d_%.4f_%.4f_timeseries_3_6months.png",
                                       point$point_id, point$y, abs(point$x)))
  output_file_trends <- file.path(output_dir,
                                  sprintf("point_%d_%.4f_%.4f_trend_analysis.png",
                                          point$point_id, point$y, abs(point$x)))
  
  # Create and save the plots
  create_four_panel_plot_1_12(all_ts_1_12, point, point$point_id, output_file_1_12)
  create_four_panel_plot_3_6(all_ts_3_6, point, point$point_id, output_file_3_6)
  create_trend_analysis_plot(all_trends, point, point$point_id, output_file_trends)
  
  # Clean up memory
  rm(spi1_stack, spi3_stack, spi6_stack, spi12_stack)
  rm(spei1_stack, spei3_stack, spei6_stack, spei12_stack)
  rm(spi1_ts, spi3_ts, spi6_ts, spi12_ts)
  rm(spei1_ts, spei3_ts, spei6_ts, spei12_ts)
  rm(all_trends)
  gc()
}

# ==================== SAVE TREND STATISTICS TO CSV ====================
if (length(all_trend_stats) > 0) {
  cat("\n========================================\n")
  cat("Saving trend statistics...\n")
  cat("========================================\n")
  
  all_stats_combined <- bind_rows(all_trend_stats)
  output_csv <- file.path(output_dir, "trend_statistics_summary.csv")
  write.csv(all_stats_combined, output_csv, row.names = FALSE)
  cat("Trend statistics saved to:", output_csv, "\n")
  
  cat("\nTrend Analysis Summary:\n")
  cat("=======================\n")
  for (idx in unique(all_stats_combined$Index)) {
    idx_data <- all_stats_combined %>% filter(Index == idx)
    n_sig_lr <- sum(idx_data$LR_significant, na.rm = TRUE)
    n_sig_sens <- sum(idx_data$Sens_significant, na.rm = TRUE)
    n_sig_mk <- sum(idx_data$MK_significant, na.rm = TRUE)
    n_total <- nrow(idx_data)
    avg_completeness <- mean(idx_data$Completeness_pct, na.rm = TRUE)
    
    cat(sprintf("\n%s (avg %.1f%% complete data):\n", idx, avg_completeness))
    cat(sprintf(" Linear Regression: %d/%d points show significant trend\n", n_sig_lr, n_total))
    cat(sprintf(" Sen's Slope: %d/%d points show significant trend\n", n_sig_sens, n_total))
    cat(sprintf(" Mann-Kendall: %d/%d points show significant trend\n", n_sig_mk, n_total))
  }
}

cat("\n==========================================\n")
cat("Analysis complete!\n")
cat("All plots and statistics saved in:", output_dir, "\n")
cat("\n⚠️ IMPORTANT NOTES:\n")
cat(" - Data gaps are now visible in ALL plots\n")
cat(" - Check console output for gap diagnostics\n")
cat(" - Trend statistics include data completeness info\n")
cat(" - Investigate why SPI-12/SPI-6 have gaps after 2020\n")
cat("==========================================\n")