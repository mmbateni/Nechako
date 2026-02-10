# ============================================================================
# SCRIPT: Basin-Mean Drought Index Time Series Analysis (PARALLEL)
# ============================================================================
# Combines:
#  - The main extraction workflow (parallel-safe terra, caching, layer chunking)
# ===================== LIBRARIES =====================
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(lubridate)
library(ncdf4)
library(sf)
library(future.apply)
library(data.table)

# ===================== CONFIGURATION =====================
setwd("D:/Nechako_Drought/")

n_cores <- max(1, parallel::detectCores() - 1)
cat("Using", n_cores, "cores for parallel processing\n")

# IMPORTANT: Use 'multisession' for Windows, 'multicore' for Linux/Mac
plan(multisession, workers = n_cores)

spi_data_dir  <- "spi_results/"
spei_data_dir <- "spei_results/"
basin_shp_path <- "Spatial/nechakoBound_dissolve.shp"

output_dir <- "basin_mean_timeseries_plots/"
cache_dir  <- "cache/"
dir.create(output_dir, showWarnings = FALSE)
dir.create(cache_dir, showWarnings = FALSE)

EQUAL_AREA_CRS <- "EPSG:3005"
DROUGHT_THRESHOLD <- -1.0
CHUNK_SIZE <- 100

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

# ===================== CACHE MANAGEMENT =====================
get_cache_path <- function(index_name, metric = "timeseries") {
  file.path(cache_dir, paste0(index_name, "_", metric, ".rds"))
}

load_from_cache <- function(index_name, metric = "timeseries") {
  cache_path <- get_cache_path(index_name, metric)
  if (file.exists(cache_path)) {
    cat(" Loading from cache:", basename(cache_path), "\n")
    return(readRDS(cache_path))
  }
  return(NULL)
}

save_to_cache <- function(data, index_name, metric = "timeseries") {
  cache_path <- get_cache_path(index_name, metric)
  saveRDS(data, cache_path)
  cat(" Saved to cache:", basename(cache_path), "\n")
}

# ===================== EXTRACT DATES =====================
extract_dates_from_nc <- function(nc_file_path) {
  cache_key <- paste0("dates_", basename(nc_file_path))
  cached <- load_from_cache(cache_key, "dates")
  if (!is.null(cached)) return(cached)
  
  cat(" Extracting dates from:", basename(nc_file_path), "\n")
  r_temp <- rast(nc_file_path)
  n_layers <- nlyr(r_temp)
  cat(" Number of layers in file:", n_layers, "\n")
  
  dates <- tryCatch({
    # 1) Try terra time()
    dates_raw <- time(r_temp)
    if (!is.null(dates_raw) && length(dates_raw) > 0) {
      cat(" Using terra time() method\n")
      dates_converted <- tryCatch(
        as.Date(dates_raw),
        error = function(e) as.Date(dates_raw, origin = "1970-01-01")
      )
      
      if (all(!is.na(dates_converted)) &&
          all(dates_converted >= as.Date("1920-01-01")) &&
          all(dates_converted <= as.Date("2030-01-01"))) {
        
        # Match length
        if (length(dates_converted) != n_layers) {
          dates_converted <- dates_converted[1:min(length(dates_converted), n_layers)]
          if (length(dates_converted) < n_layers) {
            last_date <- tail(dates_converted, 1)
            missing_dates <- seq(last_date + months(1), by = "month",
                                 length.out = n_layers - length(dates_converted))
            dates_converted <- c(dates_converted, missing_dates)
          }
        }
        return(dates_converted)
      }
    }
    
    # 2) Try ncdf4 time units
    nc <- nc_open(nc_file_path)
    on.exit(nc_close(nc), add = TRUE)
    
    if ("time" %in% names(nc$var)) {
      time_vals  <- ncvar_get(nc, "time")
      time_units <- ncatt_get(nc, "time", "units")$value
      
      if (grepl("since", time_units)) {
        origin_str  <- sub(".*since ", "", time_units)
        origin_date <- as.Date(origin_str)
        
        if (grepl("^days", time_units)) {
          dates <- origin_date + time_vals
        } else if (grepl("^months", time_units)) {
          dates <- seq(origin_date, by = "month",
                       length.out = min(length(time_vals), n_layers))
        }
        
        if (all(dates >= as.Date("1920-01-01")) &&
            all(dates <= as.Date("2030-01-01"))) {
          cat(" Successfully extracted dates using ncdf4\n")
          
          if (length(dates) != n_layers) {
            dates <- dates[1:min(length(dates), n_layers)]
            if (length(dates) < n_layers) {
              last_date <- tail(dates, 1)
              missing_dates <- seq(last_date + months(1), by = "month",
                                   length.out = n_layers - length(dates))
              dates <- c(dates, missing_dates)
            }
          }
          return(as.Date(dates))
        }
      }
    }
    
    # 3) Fallback generation
    cat(" Using fallback date generation\n")
    start_date <- as.Date("1954-01-01")
    dates <- seq(start_date, by = "month", length.out = n_layers)
    
    cat(" Generated dates from", format(dates[1], "%Y-%m-%d"),
        "to", format(tail(dates, 1), "%Y-%m-%d"), "\n")
    return(dates)
    
  }, error = function(e) {
    cat(" Date extraction error:", e$message, "\n")
    dates <- seq(as.Date("1954-01-01"), by = "month", length.out = n_layers)
    return(dates)
  })
  
  save_to_cache(dates, cache_key, "dates")
  return(dates)
}

# ===================== PRECOMPUTE BASIN GEOMETRY =====================
# Returns ONLY serializable objects (no terra objects) so workers won't fail.
precompute_basin_geometry <- function(basin_shp_path, raster_file_path,
                                      equal_area_crs = "EPSG:3005") {
  cat(" Precomputing basin geometry...\n")
  
  basin_polygon <- vect(basin_shp_path)
  raster_template <- rast(raster_file_path)[[1]]
  
  raster_proj <- project(raster_template, equal_area_crs)
  basin_proj  <- project(basin_polygon, equal_area_crs)
  
  total_basin_area <- expanse(basin_proj, unit = "m")
  
  cell_areas <- cellSize(raster_proj, unit = "m")
  areas_masked <- mask(cell_areas, basin_proj)
  
  basin_rast <- rasterize(basin_proj, raster_proj, cover = TRUE)
  
  area_vals     <- as.vector(values(areas_masked, dataframe = FALSE, na.rm = FALSE))
  coverage_vals <- as.vector(values(basin_rast, dataframe = FALSE, na.rm = FALSE))
  
  effective_areas <- area_vals * coverage_vals
  valid_cell_idx  <- which(!is.na(coverage_vals) & coverage_vals > 0)
  
  cat(" Total basin area:", round(total_basin_area / 1e6, 2), "km²\n")
  cat(" Valid cells:", length(valid_cell_idx), "\n")
  
  return(list(
    total_basin_area = total_basin_area,
    effective_areas  = effective_areas,
    coverage_vals    = coverage_vals,
    valid_cell_idx   = valid_cell_idx,
    equal_area_crs   = equal_area_crs,
    basin_shp_path   = basin_shp_path,
    raster_file_path = raster_file_path
  ))
}

# ===================== COMPUTE LAYER STATISTICS =====================
# Reconstruct terra objects inside each worker (parallel-safe).
compute_layer_statistics <- function(layer_idx, raster_file_path, basin_geom,
                                     drought_threshold = -1.0) {
  raster_stack <- rast(raster_file_path)
  raster_layer <- raster_stack[[layer_idx]]
  
  raster_template <- raster_stack[[1]]
  raster_proj <- project(raster_template, basin_geom$equal_area_crs)
  raster_layer_proj <- project(raster_layer, raster_proj)
  
  raster_vals <- as.vector(values(raster_layer_proj, dataframe = FALSE, na.rm = FALSE))
  
  valid_cells  <- basin_geom$valid_cell_idx
  cell_values  <- raster_vals[valid_cells]
  cell_areas   <- basin_geom$effective_areas[valid_cells]
  
  data_valid_idx <- !is.na(cell_values) & is.finite(cell_values)
  
  n_cells_valid <- sum(data_valid_idx)
  valid_area <- sum(cell_areas[data_valid_idx], na.rm = TRUE)
  
  if (n_cells_valid == 0 || valid_area == 0) {
    return(list(
      mean_value = NA_real_,
      valid_fraction = 0,
      drought_fraction = NA_real_,
      drought_area = 0,
      n_cells_valid = 0,
      n_cells_drought = 0
    ))
  }
  
  weighted_sum <- sum(cell_values[data_valid_idx] * cell_areas[data_valid_idx])
  mean_value <- weighted_sum / valid_area
  valid_fraction <- valid_area / basin_geom$total_basin_area
  
  drought_idx <- data_valid_idx & (cell_values <= drought_threshold)
  drought_area <- sum(cell_areas[drought_idx], na.rm = TRUE)
  drought_fraction <- drought_area / basin_geom$total_basin_area
  
  return(list(
    mean_value = mean_value,
    valid_fraction = valid_fraction,
    drought_fraction = drought_fraction,
    drought_area = drought_area,
    n_cells_valid = n_cells_valid,
    n_cells_drought = sum(drought_idx)
  ))
}

# ===================== PARALLEL BASIN MEAN TIME SERIES =====================
extract_basin_mean_timeseries_parallel <- function(raster_file_path, basin_shp_path,
                                                   index_name,
                                                   equal_area_crs = "EPSG:3005",
                                                   drought_threshold = -1.0,
                                                   use_cache = TRUE) {
  cat(" Extracting basin-mean time series for", index_name, "\n")
  
  if (use_cache) {
    cached_data <- load_from_cache(index_name, "timeseries")
    if (!is.null(cached_data)) return(cached_data)
  }
  
  dates <- extract_dates_from_nc(raster_file_path)
  
  raster_stack <- rast(raster_file_path)
  n_layers <- nlyr(raster_stack)
  rm(raster_stack)
  gc()
  
  cat(" Processing", n_layers, "time steps\n")
  
  basin_geom <- precompute_basin_geometry(basin_shp_path, raster_file_path, equal_area_crs)
  
  cat(" Processing layers in parallel...\n")
  n_chunks <- ceiling(n_layers / CHUNK_SIZE)
  all_results <- list()
  
  for (chunk_idx in 1:n_chunks) {
    start_idx <- (chunk_idx - 1) * CHUNK_SIZE + 1
    end_idx   <- min(chunk_idx * CHUNK_SIZE, n_layers)
    chunk_layers <- start_idx:end_idx
    
    cat(sprintf(" Chunk %d/%d (layers %d-%d)...\n", chunk_idx, n_chunks, start_idx, end_idx))
    
    chunk_results <- future_lapply(chunk_layers, function(i) {
      compute_layer_statistics(i, raster_file_path, basin_geom, drought_threshold)
    }, future.seed = TRUE)
    
    all_results <- c(all_results, chunk_results)
    gc()
  }
  
  results_dt <- data.table(
    Date = dates[1:n_layers],
    Mean_Value = sapply(all_results, function(x) x$mean_value),
    Valid_Fraction = sapply(all_results, function(x) x$valid_fraction),
    Drought_Fraction = sapply(all_results, function(x) x$drought_fraction),
    Drought_Valid_Fraction = sapply(all_results, function(x) x$valid_fraction),
    Index = index_name
  )
  
  results <- as.data.frame(results_dt)
  
  cat(" Mean value range:",
      round(min(results$Mean_Value, na.rm = TRUE), 2), "to",
      round(max(results$Mean_Value, na.rm = TRUE), 2), "\n")
  cat(" Average valid fraction:",
      round(mean(results$Valid_Fraction, na.rm = TRUE) * 100, 1), "%\n")
  cat(" Average drought extent:",
      round(mean(results$Drought_Fraction, na.rm = TRUE) * 100, 1), "%\n")
  
  if (use_cache) save_to_cache(results, index_name, "timeseries")
  return(results)
}

# ============================================================================
# ===================== IMPROVED PLOTTING FUNCTION ===========================
# ============================================================================
# Dynamic y-limits based on actual data (no hard cutoff at -4 to 4).
create_basin_mean_plot_improved <- function(basin_data, output_file, title_suffix = "") {
  
  indices <- unique(basin_data$Index)
  plots_mean <- list()
  plots_drought <- list()
  
  shared_theme <- theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank()
    )
  
  drought_bands <- list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2,   fill = "#8B0000", alpha = 0.10),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2,   ymax = -1.5, fill = "#FF0000", alpha = 0.10),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.5, ymax = -1,   fill = "#FF8C00", alpha = 0.10),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1,   ymax = -0.5, fill = "#FFA500", alpha = 0.10),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5,  ymax = 1,    fill = "#90EE90", alpha = 0.10),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1,    ymax = 1.5,  fill = "#00FF00", alpha = 0.10),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5,  ymax = 2,    fill = "#008000", alpha = 0.10),
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2,    ymax = Inf,  fill = "#006400", alpha = 0.10)
  )
  
  reference_lines <- geom_hline(
    yintercept = c(-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2),
    linetype = "dotted", color = "gray50", linewidth = 0.3
  )
  
  for (index in indices) {
    index_data <- basin_data %>% filter(Index == index)
    index_data$Mean_Value[is.infinite(index_data$Mean_Value)] <- NA
    
    # Dynamic y-limits based on actual data range
    data_min <- suppressWarnings(min(index_data$Mean_Value, na.rm = TRUE))
    data_max <- suppressWarnings(max(index_data$Mean_Value, na.rm = TRUE))
    
    # Handle edge case: all NA
    if (!is.finite(data_min) || !is.finite(data_max)) {
      data_min <- -4
      data_max <-  4
    }
    
    y_range <- data_max - data_min
    if (!is.finite(y_range) || y_range == 0) y_range <- 1
    
    y_min <- min(-4, data_min - 0.1 * y_range)
    y_max <- max( 4, data_max + 0.1 * y_range)
    
    if (data_min < -4 || data_max > 4) {
      cat(sprintf("NOTE: %s has extreme values (%.2f to %.2f) - adjusting y-limits\n",
                  index, data_min, data_max))
    }
    
    p_mean <- ggplot(index_data, aes(x = Date, y = Mean_Value)) +
      drought_bands +
      reference_lines +
      geom_line(color = "blue", linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
      geom_hline(yintercept = DROUGHT_THRESHOLD, linetype = "dashed",
                 color = "red", linewidth = 0.8) +
      labs(
        title = paste(index, "- Basin-Mean (Area-Weighted)"),
        x = "",
        y = "Index Value"
      ) +
      shared_theme +
      coord_cartesian(ylim = c(y_min, y_max))  # better than ylim() (doesn't drop data)
    
    p_drought <- ggplot(index_data, aes(x = Date, y = Drought_Fraction * 100)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,  ymax = 25,  fill = "#FFA500", alpha = 0.05) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 25, ymax = 50,  fill = "#FF6B6B", alpha = 0.05) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 50, ymax = 75,  fill = "#CC0000", alpha = 0.05) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 75, ymax = 100, fill = "#8B0000", alpha = 0.05) +
      geom_hline(yintercept = c(25, 50, 75), linetype = "dotted",
                 color = "gray50", linewidth = 0.3) +
      geom_area(fill = "#FF6B6B", alpha = 0.6) +
      geom_line(color = "#CC0000", linewidth = 0.8) +
      labs(
        title = paste(index, "- Drought Extent (% of Basin Area)"),
        x = "Date",
        y = "Basin Area in Drought (%)"
      ) +
      shared_theme +
      coord_cartesian(ylim = c(0, 100))
    
    plots_mean[[index]] <- p_mean
    plots_drought[[index]] <- p_drought
  }
  
  # Combine plots
  if (length(indices) == 4) {
    combined_plot <-
      (plots_mean[[indices[1]]] | plots_drought[[indices[1]]]) /
      (plots_mean[[indices[2]]] | plots_drought[[indices[2]]]) /
      (plots_mean[[indices[3]]] | plots_drought[[indices[3]]]) /
      (plots_mean[[indices[4]]] | plots_drought[[indices[4]]]) +
      plot_annotation(
        title = paste("Nechako Basin Drought Analysis", title_suffix),
        subtitle = sprintf(
          "Area-weighted basin means (left) and drought extent (right)\nDrought threshold: %.1f | Equal-area CRS: %s\nY-limits adjusted to data",
          DROUGHT_THRESHOLD, EQUAL_AREA_CRS
        ),
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10)
        )
      )
  } else {
    combined_plot <- wrap_plots(c(plots_mean, plots_drought), ncol = 2) +
      plot_annotation(
        title = paste("Nechako Basin Drought Analysis", title_suffix),
        subtitle = sprintf(
          "Area-weighted basin means (left) and drought extent (right)\nDrought threshold: %.1f | Equal-area CRS: %s\nY-limits adjusted to data",
          DROUGHT_THRESHOLD, EQUAL_AREA_CRS
        ),
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5),
          plot.margin = margin(10, 10, 10, 10)
        )
      )
  }
  
  ggsave(output_file, combined_plot, width = 16, height = 14, dpi = 150)
  cat(" Saved plot to:", output_file, "\n")
}

# ===================== MAIN PROCESSING =====================
cat("==============================================\n")
cat("Basin-Mean Drought Index Time Series Analysis\n")
cat("(PARALLEL SAFE + IMPROVED PLOTS)\n")
cat("==============================================\n\n")

cat("Loading basin shapefile...\n")
basin_shp <- vect(basin_shp_path)
cat("Basin CRS:", crs(basin_shp), "\n")
basin_proj <- project(basin_shp, EQUAL_AREA_CRS)
basin_area_km2 <- expanse(basin_proj, unit = "km")
cat("Basin area:", round(basin_area_km2, 2), "km²\n\n")

# ===================== PROCESS SPI DATA =====================
cat("Processing SPI data...\n")
cat("======================\n")
spi_data_list <- list()

for (i in seq_along(spi_files)) {
  file_path <- file.path(spi_data_dir, spi_files[i])
  timescale <- as.integer(gsub(".*_(\\d+)month.*", "\\1", spi_files[i]))
  index_name <- paste0("SPI", timescale)
  
  cat("\nProcessing", index_name, "\n")
  spi_ts <- extract_basin_mean_timeseries_parallel(
    file_path, basin_shp_path, index_name,
    EQUAL_AREA_CRS, DROUGHT_THRESHOLD
  )
  spi_data_list[[index_name]] <- spi_ts
  gc()
}

all_spi_data <- bind_rows(spi_data_list)

# ===================== PROCESS SPEI DATA =====================
cat("\n\nProcessing SPEI data...\n")
cat("=======================\n")
spei_data_list <- list()

for (i in seq_along(spei_files)) {
  file_path <- file.path(spei_data_dir, spei_files[i])
  timescale <- as.integer(gsub(".*_(\\d+)month.*", "\\1", spei_files[i]))
  index_name <- paste0("SPEI", timescale)
  
  cat("\nProcessing", index_name, "\n")
  spei_ts <- extract_basin_mean_timeseries_parallel(
    file_path, basin_shp_path, index_name,
    EQUAL_AREA_CRS, DROUGHT_THRESHOLD
  )
  spei_data_list[[index_name]] <- spei_ts
  gc()
}

all_spei_data <- bind_rows(spei_data_list)

# ===================== CREATE IMPROVED PLOTS =====================
cat("\n\nCreating IMPROVED plots...\n")
cat("==========================\n")

all_data <- bind_rows(all_spi_data, all_spei_data)

# Save combined data (optional but useful)
output_csv <- file.path(output_dir, "basin_mean_timeseries_all.csv")
fwrite(all_data, output_csv)
cat("Saved combined data to:", output_csv, "\n\n")

# 1 & 12 month panel
data_1_12 <- all_data %>%
  filter(Index %in% c("SPI1", "SPI12", "SPEI1", "SPEI12"))

output_1_12 <- file.path(output_dir, "basin_mean_1_12_months_IMPROVED.png")
create_basin_mean_plot_improved(data_1_12, output_1_12, "- 1 & 12 Month Timescales")

# 3 & 6 month panel
data_3_6 <- all_data %>%
  filter(Index %in% c("SPI3", "SPI6", "SPEI3", "SPEI6"))

output_3_6 <- file.path(output_dir, "basin_mean_3_6_months_IMPROVED.png")
create_basin_mean_plot_improved(data_3_6, output_3_6, "- 3 & 6 Month Timescales")

cat("\nImproved plots saved!\n")

# ===================== CLEAN UP =====================
plan(sequential)
cat("\n==============================================\n")
cat("Processing complete!\n")
cat("Improved plots are in:\n")
cat(" ", output_dir, "\n")
cat("==============================================\n")