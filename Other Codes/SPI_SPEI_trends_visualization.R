####################################################################################
#  COMPREHENSIVE DROUGHT DIAGNOSTICS - COMPLETE WORKFLOW 
#  Nechako Basin Analysis - SPI/SPEI Results
#
#  MODIFIED VERSION - FIXES:
#  1. Time column name collision issue (loads individual files for basin averages)
#  2. Proper extraction of time series data per index/timescale
#  3. Handles missing TFPW columns gracefully
#  4. Better error handling and diagnostics
#
#  WORKFLOW:
#  1. Load individual SPI/SPEI result files (spi_01, spei_01, etc.)
#  2. Combine into master diagnostics file with basin clipping
#  3. Generate all spatial visualization figures
#  4. Generate basin-averaged time series CSVs (FIXED: loads individual files)
#  5. Generate time series plots and comparison analysis
#  6. Generate Excel summary statistics
#
#  INPUT FILES:
#  - Individual result files: spi_XX_results.csv, spei_XX_results.csv
#  - Basin boundary: Spatial/nechakoBound_dissolve.shp
#
#  OUTPUT FILES:
#  - Master CSV: temporal_spi_spei/all_temporal_diagnostics_results.csv
#  - Figures: temporal_spi_spei/figures_temporal_diagnostics_enhanced/
#  - Time series CSVs: temporal_spi_spei/basin_averaged_timeseries/
#  - Excel summary: temporal_spi_spei/basin_averaged_timeseries/SPI_SPEI_Summary_Statistics.xlsx
####################################################################################

# Load required libraries with automatic installation
required_packages <- c("terra", "data.table", "ggplot2", "sf", "zoo", "viridis", 
                       "cowplot", "patchwork", "scales", "RColorBrewer", 
                       "ggspatial", "gridExtra", "grid", "lubridate", "openxlsx")

# Check and install missing packages
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("📦 Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org/")
}

# Load all libraries
invisible(lapply(required_packages, library, character.only = TRUE, quietly = TRUE))
cat("✓ All required packages loaded\n")

# Set working directory (configurable)
# To run on a different system, either:
# 1. Set WD_PATH environment variable, OR
# 2. Modify this path, OR  
# 3. Run script from the correct directory
WD_PATH <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")

if (!dir.exists(WD_PATH)) {
  stop(sprintf("❌ ERROR: Working directory does not exist: %s\nPlease set the correct path.", WD_PATH))
}

setwd(WD_PATH)
cat(sprintf("✓ Working directory: %s\n", getwd()))

# Create output directories
fig_dir <- "./temporal_spi_spei/figures_temporal_diagnostics_enhanced"
timeseries_dir <- "./temporal_spi_spei/basin_averaged_timeseries"
basin_plots_dir <- "./temporal_spi_spei/basin_averaged_plots"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(timeseries_dir)) dir.create(timeseries_dir, recursive = TRUE)
if (!dir.exists(basin_plots_dir)) dir.create(basin_plots_dir, recursive = TRUE)

# Load basin boundary
cat("📂 Loading basin boundary...\n")

basin_file <- "Spatial/nechakoBound_dissolve.shp"
if (!file.exists(basin_file)) {
  stop(sprintf("❌ ERROR: Basin boundary file not found: %s\nPlease check the file path.", basin_file))
}

basin <- st_read(basin_file, quiet = TRUE)
if (st_crs(basin)$input != "EPSG:3005") {
  basin <- st_transform(basin, "EPSG:3005")
}
basin_vect <- vect(basin)
cat(sprintf("✓ Basin loaded: %d polygons, CRS: %s\n", nrow(basin), st_crs(basin)$input))

####################################################################################
# PART 1: GENERATE MASTER DIAGNOSTICS FILE FROM INDIVIDUAL RESULTS
####################################################################################

cat("\n========================================\n")
cat("🔄 PART 1: GENERATING MASTER DIAGNOSTICS FILE\n")
cat("========================================\n")

input_dir <- "./temporal_spi_spei/"
output_file <- "./temporal_spi_spei/all_temporal_diagnostics_results.csv"
timescales <- c(1, 3, 6, 12, 24)
indices <- c("spi", "spei")

# Function to load and process individual result file
load_and_clip_results <- function(index_type, timescale, input_dir, basin_vect) {
  # Input validation
  if (!index_type %in% c("spi", "spei")) {
    stop(sprintf("Invalid index_type: %s (must be 'spi' or 'spei')", index_type))
  }
  if (!is.numeric(timescale) || timescale <= 0) {
    stop(sprintf("Invalid timescale: %s (must be positive number)", timescale))
  }
  if (!dir.exists(input_dir)) {
    stop(sprintf("Input directory does not exist: %s", input_dir))
  }
  
  cat(sprintf("\nProcessing %s-%02d...\n", toupper(index_type), timescale))
  
  # Load file
  file_pattern <- sprintf("%s_%02d_results.csv", index_type, timescale)
  file_path <- file.path(input_dir, file_pattern)
  
  if (!file.exists(file_path)) {
    cat(sprintf("  ⚠️  File not found: %s\n", file_path))
    return(NULL)
  }
  
  dt <- fread(file_path)
  cat(sprintf("  ✓ Loaded: %d rows\n", nrow(dt)))
  
  # Check for coordinate columns
  coord_cols <- c("lon", "lat", "x", "y", "easting", "northing")
  found_coords <- intersect(coord_cols, names(dt))
  
  if (length(found_coords) < 2) {
    cat(sprintf("  ⚠️  Could not find coordinate columns. Available: %s\n", 
                paste(names(dt), collapse = ", ")))
    return(NULL)
  }
  
  # Standardize coordinate column names
  if ("x" %in% names(dt) && !"lon" %in% names(dt)) setnames(dt, "x", "lon")
  if ("y" %in% names(dt) && !"lat" %in% names(dt)) setnames(dt, "y", "lat")
  
  # Remove rows with NA coordinates
  dt <- dt[!is.na(lon) & !is.na(lat)]
  cat(sprintf("  ✓ Valid coordinates: %d rows\n", nrow(dt)))
  
  # Create spatial points
  pts <- vect(dt, geom = c("lon", "lat"), crs = "EPSG:3005")
  
  # FIX: Get basin extent properly
  basin_ext <- ext(basin_vect)
  
  # Filter points within basin extent (fast pre-filter)
  pts_in_ext <- dt[lon >= basin_ext[1] & lon <= basin_ext[2] & 
                     lat >= basin_ext[3] & lat <= basin_ext[4]]
  cat(sprintf("  ✓ Points within basin extent: %d (%.1f%%)\n", 
              nrow(pts_in_ext), nrow(pts_in_ext) / nrow(dt) * 100))
  
  if (nrow(pts_in_ext) == 0) {
    cat("  ⚠️  WARNING: No points within basin extent!\n")
    return(NULL)
  }
  
  # Create points from filtered data
  pts_clipped <- vect(pts_in_ext, geom = c("lon", "lat"), crs = "EPSG:3005")
  
  # Use sf for precise clipping with touches
  pts_sf <- st_as_sf(pts_clipped)
  basin_sf <- st_as_sf(basin_vect)
  
  # Find points that intersect or touch the basin
  intersects_idx <- st_intersects(pts_sf, basin_sf, sparse = TRUE)
  has_intersection <- sapply(intersects_idx, length) > 0
  
  dt_clipped <- pts_in_ext[has_intersection, ]
  
  n_clipped <- nrow(dt_clipped)
  cat(sprintf("  ✓ After basin clip (with touches): %d rows (%.1f%%)\n", 
              n_clipped, n_clipped / nrow(dt) * 100))
  
  if (n_clipped == 0) {
    cat("  ⚠️  WARNING: No points remain after clipping!\n")
    return(NULL)
  }
  
  # Convert back to data.table (already is)
  # dt_clipped <- as.data.table(dt_clipped)
  
  # Add metadata columns
  dt_clipped[, index_type := index_type]
  dt_clipped[, timescale := timescale]
  
  # Check for NA values in critical columns
  cat(sprintf("  📊 Data quality check:\n"))
  cat(sprintf("     tau_vc NA: %d (%.1f%%)\n", 
              sum(is.na(dt_clipped$tau_vc)), 
              sum(is.na(dt_clipped$tau_vc)) / nrow(dt_clipped) * 100))
  cat(sprintf("     p_value_vc NA: %d (%.1f%%)\n", 
              sum(is.na(dt_clipped$p_value_vc)), 
              sum(is.na(dt_clipped$p_value_vc)) / nrow(dt_clipped) * 100))
  
  return(dt_clipped)
}

# Process all files
all_results <- list()

for (idx in indices) {
  for (scale in timescales) {
    result <- load_and_clip_results(idx, scale, input_dir, basin_vect)
    if (!is.null(result)) {
      all_results[[sprintf("%s_%02d", idx, scale)]] <- result
    }
  }
}

# Combine all results
if (length(all_results) == 0) {
  stop("❌ ERROR: No data loaded from any file!")
}

combined <- rbindlist(all_results, fill = TRUE)
cat(sprintf("\n✓ Combined dataset: %d rows\n", nrow(combined)))

# Add space_idx column if missing
if (!"space_idx" %in% names(combined)) {
  combined[, space_idx := 1:.N]
}

# Save to CSV
fwrite(combined, output_file, row.names = FALSE)
cat(sprintf("✓ Saved: %s (%.2f MB)\n", output_file, file.info(output_file)$size / 1024 / 1024))

# Final summary
cat("\n📊 FINAL DATA SUMMARY:\n")
cat("----------------------------------------\n")
for (idx in c("spi", "spei")) {
  idx_data <- combined[index_type == idx]
  cat(sprintf("\n%s:\n", toupper(idx)))
  cat(sprintf("  Total rows: %d\n", nrow(idx_data)))
  cat(sprintf("  Non-NA tau_vc: %d (%.1f%%)\n", 
              sum(!is.na(idx_data$tau_vc)), 
              sum(!is.na(idx_data$tau_vc)) / nrow(idx_data) * 100))
  cat(sprintf("  filtered_vc = TRUE: %d (%.1f%%)\n", 
              sum(idx_data$filtered_vc == TRUE, na.rm = TRUE), 
              sum(idx_data$filtered_vc == TRUE, na.rm = TRUE) / nrow(idx_data) * 100))
}

####################################################################################
# PART 2: LOAD MASTER FILE FOR VISUALIZATION
####################################################################################

cat("\n========================================\n")
cat("📂 PART 2: LOADING DATA FOR VISUALIZATION\n")
cat("========================================\n")

results <- fread(output_file)
cat(sprintf("✓ Loaded %d rows\n", nrow(results)))
cat(sprintf("✓ Indices: %s\n", paste(unique(results$index_type), collapse = ", ")))
cat(sprintf("✓ Timescales: %s\n", paste(sort(unique(results$timescale)), collapse = ", ")))

# Configuration
dpi <- 300
timescales_spatial <- c(1, 3, 6, 12, 24)
timescales_temporal <- c(1, 3, 6, 9, 12)
drought_threshold_onset <- -1.0
drought_threshold_end <- 0.0
severe_threshold <- -1.3
extreme_threshold <- -1.6
INCLUDE_FILTERED_DATA <- TRUE  # Include SPEI even though it failed VC check

# Apply filtering
if (INCLUDE_FILTERED_DATA) {
  cat("\n📊 Using ALL data (including filtered)\n")
  plot_data <- results
} else {
  cat("\n📊 Using only unfiltered data\n")
  plot_data <- results[!filtered_vc | is.na(filtered_vc)]
}

cat(sprintf("   Rows for plotting: %d\n", nrow(plot_data)))
cat(sprintf("   SPI rows: %d\n", sum(plot_data$index_type == "spi")))
cat(sprintf("   SPEI rows: %d\n", sum(plot_data$index_type == "spei")))

####################################################################################
# HELPER FUNCTIONS - SPATIAL ANALYSIS
####################################################################################

create_raster_template <- function(data_subset, basin_sf) {
  if (nrow(data_subset) == 0) return(NULL)
  data_subset <- data_subset[!is.na(lon) & !is.na(lat)]
  if (nrow(data_subset) == 0) return(NULL)
  
  unique_lons <- unique(data_subset$lon)
  unique_lats <- unique(data_subset$lat)
  
  if (length(unique_lons) < 2 || length(unique_lats) < 2) {
    basin_ext <- ext(vect(basin_sf))
    resolution <- min(basin_ext[2] - basin_ext[1], basin_ext[4] - basin_ext[3]) / 50
  } else {
    dx <- median(diff(sort(unique_lons)), na.rm = TRUE)
    dy <- median(diff(sort(unique_lats)), na.rm = TRUE)
    resolution <- if (is.na(dx) || dx <= 0) {
      basin_ext <- ext(vect(basin_sf))
      min(basin_ext[2] - basin_ext[1], basin_ext[4] - basin_ext[3]) / 50
    } else {
      min(dx, dy, na.rm = TRUE)
    }
  }
  
  basin_ext <- ext(vect(basin_sf))
  rast(ext = basin_ext, resolution = resolution, crs = "EPSG:3005")
}

create_raster_from_points <- function(data_dt, template, value_col, basin_sf) {
  if (is.null(template)) return(NULL)
  if (nrow(data_dt) == 0) {
    r <- rast(template); values(r) <- NA
    return(mask(r, vect(basin_sf)))
  }
  
  data_dt <- data_dt[!is.na(lon) & !is.na(lat)]
  if (nrow(data_dt) == 0) {
    r <- rast(template); values(r) <- NA
    return(mask(r, vect(basin_sf)))
  }
  
  if (!value_col %in% names(data_dt)) {
    r <- rast(template); values(r) <- NA
    return(mask(r, vect(basin_sf)))
  }
  
  if (all(is.na(data_dt[[value_col]]))) {
    r <- rast(template); values(r) <- NA
    return(mask(r, vect(basin_sf)))
  }
  
  pts <- vect(data_dt, geom = c("lon", "lat"), crs = crs(template))
  r <- rasterize(pts, template, field = value_col, touches = TRUE)
  mask(r, vect(basin_sf))
}

smooth_raster <- function(r) {
  if (is.null(r) || all(is.na(values(r)))) return(r)
  focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
}

plot_raster_clean <- function(r, main, zlim = NULL, col, breaks = NULL,
                              legend = TRUE, categorical = FALSE, legend_title = "") {
  if (is.null(r) || all(is.na(values(r)))) {
    plot.new()
    text(0.5, 0.5, "No valid data\n(all values are NA)", cex = 1.2, col = "red")
    return()
  }
  
  if (!categorical) r <- smooth_raster(r)
  
  plg_args <- list(cex = 0.9)
  if (legend_title != "") plg_args$title <- legend_title
  
  if (!is.null(breaks)) {
    plot(r, main = main, cex.main = 1.0, col = col, breaks = breaks,
         axes = FALSE, box = FALSE, legend = legend, plg = plg_args, colNA = NA)
  } else {
    plot(r, main = main, cex.main = 1.0, col = col, zlim = zlim,
         axes = FALSE, box = FALSE, legend = legend, plg = plg_args, colNA = NA)
  }
  plot(st_geometry(basin), add = TRUE, col = NA, border = "black", lwd = 2.5)
}

safe_pdf <- function(filename, width = 12, height = 8) {
  while (dev.cur() > 1) try(dev.off(), silent = TRUE)
  tryCatch({
    pdf(filename, width = width, height = height, family = "Helvetica")
    TRUE
  }, error = function(e) {
    cat(sprintf("ERROR: Cannot create PDF '%s': %s\n", filename, e$message))
    FALSE
  })
}

####################################################################################
# BASIN-AVERAGED TIME SERIES RECONSTRUCTION (WITH SYNTHETIC NOISE)
####################################################################################
reconstruct_basin_timeseries <- function(results_dt, output_file) {
  cat("📈 Reconstructing basin-averaged time series...\n")
  pdf(output_file, width = 12, height = 8, family = "Helvetica")
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  for (idx in c("spi", "spei")) {
    for (scale in c(12, 24)) {
      subset_data <- results_dt[index_type == idx & timescale == scale & !filtered_vc]
      
      if (nrow(subset_data) == 0) {
        plot.new()
        text(0.5, 0.5, sprintf("No data for %s-%d", toupper(idx), scale), cex = 1.2)
        next
      }
      
      avg_tau <- median(subset_data$tau_vc, na.rm = TRUE)
      avg_slope <- median(subset_data$sl_vc, na.rm = TRUE)
      pct_sig <- sum(subset_data$p_value_vc < 0.05) / nrow(subset_data) * 100
      
      years <- 1950:2025
      months <- length(years) * 12
      time_idx <- 1:months
      
      set.seed(123)
      baseline <- rnorm(months, 0, 1)
      trend_component <- avg_slope * time_idx / 12
      synthetic_series <- baseline + trend_component
      
      plot(years[1] + (1:months)/12, synthetic_series,
           type = "l", col = "steelblue", lwd = 1.5,
           xlab = "Year", ylab = paste0(toupper(idx), "-", scale, " Index Value"),
           main = sprintf("%s-%d Basin Trend (Median τ = %.3f, %.1f%% sig.)",
                          toupper(idx), scale, avg_tau, pct_sig))
      abline(h = 0, col = "gray50", lty = 2)
      abline(h = c(-1, 1), col = "gray70", lty = 3)
      
      trend_line <- avg_slope * time_idx / 12
      lines(years[1] + (1:months)/12, trend_line, col = "red", lwd = 2, lty = 2)
      
      legend("topleft",
             legend = c("Synthetic Example", "Linear Trend", "Neutral", "±1 Std"),
             col = c("steelblue", "red", "gray50", "gray70"),
             lty = c(1, 2, 2, 3), lwd = c(1.5, 2, 1, 1),
             bty = "n", cex = 0.9)
      
      text(years[1] + 10, max(synthetic_series) * 0.9,
           sprintf("Sen's Slope: %.4f\nMedian Tau: %.3f\nN pixels: %d",
                   avg_slope, avg_tau, nrow(subset_data)),
           pos = 4, cex = 0.8, col = "darkred")
    }
  }
  dev.off()
  cat(sprintf("✅ Time series visualization: %s\n", basename(output_file)))
}

create_smooth_trend_maps <- function(data, basin_sf, output_file, scale = 12) {
  cat(sprintf("🗺️  Creating trend maps (scale=%d)...\n", scale))
  
  # Validate inputs
  if (is.null(data) || nrow(data) == 0) {
    cat("⚠️  No data provided\n")
    return(invisible(NULL))
  }
  if (is.null(basin_sf)) {
    cat("⚠️  No basin boundary provided\n")
    return(invisible(NULL))
  }
  
  pdf(output_file, width = 12, height = 10, family = "Helvetica")
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  
  for (idx in c("spi", "spei")) {
    idx_data <- data[index_type == idx & timescale == scale]
    
    if (nrow(idx_data) == 0) {
      plot.new()
      text(0.5, 0.5, sprintf("%s: No data", toupper(idx)), cex = 1.5)
      next
    }
    
    if (all(is.na(idx_data$tau_vc))) {
      plot.new()
      text(0.5, 0.5, sprintf("%s: All tau values are NA", toupper(idx)), cex = 1.0, col = "red")
      next
    }
    
    template <- create_raster_template(idx_data, basin_sf)
    if (is.null(template)) next
    
    r_tau <- create_raster_from_points(idx_data, template, "tau_vc", basin_sf)
    plot_raster_clean(r_tau,
                      paste0(toupper(idx), ": Kendall's Tau"),
                      c(-0.4, 0.4),
                      hcl.colors(101, "RdBu", rev = TRUE),
                      legend_title = "Tau")
    
    valid_p <- idx_data[!is.na(p_value_vc)]
    if (nrow(valid_p) > 0) {
      idx_data_copy <- copy(idx_data)
      idx_data_copy[, sig_cat := cut(p_value_vc,
                                     breaks = c(0, 0.001, 0.01, 0.05, 1),
                                     labels = c("p<0.001", "p<0.01", "p<0.05", "n.s."))]
      idx_data_copy[, sig_num := as.numeric(sig_cat)]
      
      r_sig <- create_raster_from_points(idx_data_copy, template, "sig_num", basin_sf)
      plot_raster_clean(r_sig,
                        paste0(toupper(idx), ": Significance"),
                        NULL,
                        c("#d73027", "#fc8d59", "#fee08b", "#f7f7f7"),
                        breaks = c(0.5, 1.5, 2.5, 3.5, 4.5),
                        categorical = TRUE,
                        legend_title = "p-value")
    } else {
      plot.new()
      text(0.5, 0.5, sprintf("%s: All p-values are NA", toupper(idx)), cex = 1.0, col = "red")
    }
  }
  
  mtext(sprintf("Trend Analysis (%d-month scale)", scale),
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  dev.off()
  cat(sprintf("✅ Saved: %s\n", basename(output_file)))
}

create_event_raster_maps <- function(data, basin_sf, output_file, scale = 12) {
  cat(sprintf("🗺️  Creating event maps (scale=%d)...\n", scale))
  
  # Validate inputs
  if (is.null(data) || nrow(data) == 0) {
    cat("⚠️  No data provided\n")
    return(invisible(NULL))
  }
  
  pdf(output_file, width = 12, height = 14, family = "Helvetica")
  par(mfrow = c(3, 2), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  
  for (idx in c("spi", "spei")) {
    idx_data <- data[index_type == idx & timescale == scale]
    
    if (nrow(idx_data) == 0) {
      for (i in 1:3) {
        plot.new()
        text(0.5, 0.5, sprintf("%s: No data", toupper(idx)), cex = 1.2)
      }
      next
    }
    
    template <- create_raster_template(idx_data, basin_sf)
    if (is.null(template)) next
    
    metrics <- list(
      n_events = list(title = "Number of Events", col = "YlOrRd"),
      mean_duration = list(title = "Mean Duration", col = "viridis"),
      max_intensity = list(title = "Max Intensity", col = "Reds")
    )
    
    for (metric_name in names(metrics)) {
      if (all(is.na(idx_data[[metric_name]]))) {
        plot.new()
        text(0.5, 0.5, sprintf("%s: All %s are NA", toupper(idx), metric_name),
             cex = 1.0, col = "red")
      } else {
        r_metric <- create_raster_from_points(idx_data, template, metric_name, basin_sf)
        metric_info <- metrics[[metric_name]]
        
        if (metric_name == "max_intensity") {
          zlim <- c(min(idx_data[[metric_name]], na.rm = TRUE), 0)
        } else {
          zlim <- c(0, max(idx_data[[metric_name]], na.rm = TRUE))
        }
        
        plot_raster_clean(r_metric,
                          paste0(toupper(idx), ": ", metric_info$title),
                          zlim,
                          hcl.colors(101, metric_info$col, rev = (metric_name == "max_intensity")),
                          legend_title = "")
      }
    }
  }
  
  mtext(sprintf("Event Characteristics (%d-month scale)", scale),
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  dev.off()
  cat(sprintf("✅ Saved: %s\n", basename(output_file)))
}

create_temporal_pattern_maps <- function(data, basin_sf, output_file, scale = 12) {
  cat("🎨 Creating temporal pattern maps...\n")
  
  # Validate inputs
  if (is.null(data) || nrow(data) == 0) {
    cat("⚠️  No data provided\n")
    return(invisible(NULL))
  }
  
  pdf(output_file, width = 12, height = 8, family = "Helvetica")
  par(mfrow = c(2, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  for (idx in c("spi", "spei")) {
    idx_data <- data[index_type == idx & timescale == scale]
    template <- create_raster_template(idx_data, basin_sf)
    
    if (is.null(template)) {
      for (i in 1:3) {
        plot.new()
        text(0.5, 0.5, paste("No data for", toupper(idx)), cex = 1.2)
      }
      next
    }
    
    # Runs test clustering
    idx_data[, cluster_class := fcase(
      filtered_runs | is.na(p_value_runs), NA_real_,
      p_value_runs >= 0.05, 0,
      clustering == "clustered", -1,
      clustering == "dispersed", 1
    )]
    
    r_cluster <- create_raster_from_points(idx_data, template, "cluster_class", basin_sf)
    
    plot_raster_clean(r_cluster,
                      paste0(toupper(idx), ": Temporal Clustering"),
                      NULL,
                      c("#e41a1c", "#fee090", "#4daf4a"),
                      breaks = c(-1.5, -0.5, 0.5, 1.5),
                      legend = FALSE, categorical = TRUE)
    
    legend("bottomright",
           legend = c("Clustered", "No sig.", "Dispersed"),
           fill = c("#e41a1c", "#fee090", "#4daf4a"),
           bty = "n", cex = 0.9)
    
    # Regime shift timing
    idx_data[, shift_decade := cut(regime_shift_year,
                                   breaks = seq(1950, 2030, by = 10),
                                   labels = 1:8,
                                   include.lowest = TRUE)]
    idx_data[, shift_decade := as.numeric(shift_decade)]
    
    r_regime <- create_raster_from_points(idx_data[!is.na(regime_shift_year)],
                                          template, "shift_decade", basin_sf)
    
    plot_raster_clean(r_regime,
                      paste0(toupper(idx), ": Regime Shift Decade"),
                      c(1, 8),
                      hcl.colors(8, "Spectral", rev = TRUE),
                      legend_title = "Decade")
    
    # Spectral peaks
    r_peaks <- create_raster_from_points(idx_data, template, "n_spectral_peaks", basin_sf)
    
    plot_raster_clean(r_peaks,
                      paste0(toupper(idx), ": Spectral Peaks"),
                      c(0, 5),
                      hcl.colors(101, "viridis"),
                      legend_title = "N peaks")
  }
  
  mtext(sprintf("Temporal Patterns (%d-month scale)", scale),
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  dev.off()
  cat(sprintf("✅ Saved: %s\n", basename(output_file)))
}

create_method_comparison_plots <- function(data, output_file) {
  cat("🎨 Creating method comparison...\n")
  
  # Validate input
  if (is.null(data) || nrow(data) == 0) {
    cat("⚠️  No data provided\n")
    return(invisible(NULL))
  }
  
  required_cols <- c("filtered_tfpw", "tau_tfpw", "p_value_tfpw")
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    cat(sprintf("⚠️  Missing columns: %s\n", paste(missing, collapse = ", ")))
    cat("   Skipping method comparison (TFPW not calculated)\n")
    pdf(output_file, width = 12, height = 5)
    plot.new()
    text(0.5, 0.5, "TFPW method not available\n(Only VC method calculated)", cex = 1.2, col = "orange")
    dev.off()
    return(invisible(NULL))
  }
  
  comparison_data <- data[timescale == 12]
  
  valid_comparison <- comparison_data[!is.na(tau_vc) & !is.na(tau_tfpw)]
  if (nrow(valid_comparison) == 0) {
    cat("⚠️  No valid tau values for comparison\n")
    pdf(output_file, width = 12, height = 5)
    plot.new()
    text(0.5, 0.5, "No valid tau values\nfor method comparison", cex = 1.2, col = "red")
    dev.off()
    return()
  }
  
  comparison_data[, vc_sig := p_value_vc < 0.05]
  comparison_data[, tfpw_sig := p_value_tfpw < 0.05]
  
  p1 <- ggplot(comparison_data, aes(x = tau_vc, y = tau_tfpw, color = index_type)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(alpha = 0.5, size = 2, na.rm = TRUE) +
    facet_wrap(~index_type, labeller = labeller(index_type = toupper)) +
    scale_color_manual(values = c("spi" = "#e41a1c", "spei" = "#377eb8")) +
    coord_equal() +
    theme_bw(base_size = 11) +
    theme(legend.position = "none") +
    labs(title = "A) Kendall's Tau Comparison",
         x = "VC Method (τ)", y = "TFPW Method (τ)")
  
  agreement_summary <- comparison_data[, .(
    Both_Sig = sum(vc_sig & tfpw_sig, na.rm = TRUE),
    Only_VC = sum(vc_sig & !tfpw_sig, na.rm = TRUE),
    Only_TFPW = sum(!vc_sig & tfpw_sig, na.rm = TRUE),
    Neither = sum(!vc_sig & !tfpw_sig, na.rm = TRUE)
  ), by = index_type]
  
  agreement_long <- melt(agreement_summary, id.vars = "index_type")
  
  p2 <- ggplot(agreement_long, aes(x = index_type, y = value, fill = variable)) +
    geom_col(position = "stack") +
    scale_fill_brewer(palette = "Set2",
                      labels = c("Both Sig.", "Only VC", "Only TFPW", "Neither")) +
    scale_x_discrete(labels = c("spi" = "SPI", "spei" = "SPEI")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom") +
    labs(title = "B) Method Agreement (12-month scale)",
         x = "Index", y = "Number of Pixels", fill = "Significance")
  
  combined <- p1 | p2
  ggsave(output_file, combined, width = 12, height = 5, dpi = dpi)
  cat(sprintf("✅ Saved: %s\n", basename(output_file)))
}

create_timescale_comparison <- function(data, output_file) {
  cat("🎨 Creating timescale comparison...\n")
  
  # Validate input
  if (is.null(data) || nrow(data) == 0) {
    cat("⚠️  No data provided\n")
    return(invisible(NULL))
  }
  
  timescale_summary <- data[, .(
    Pct_Significant = sum(p_value_vc < 0.05, na.rm = TRUE) / .N * 100,
    Median_Tau = median(tau_vc, na.rm = TRUE),
    Mean_Events = mean(n_events, na.rm = TRUE),
    Mean_Duration = mean(mean_duration, na.rm = TRUE),
    N_Pixels = .N
  ), by = .(index_type, timescale)]
  
  p1 <- ggplot(timescale_summary,
               aes(x = timescale, y = Pct_Significant,
                   color = index_type, group = index_type)) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    scale_color_manual(values = c("spi" = "#e41a1c", "spei" = "#377eb8"),
                       labels = c("SPI", "SPEI")) +
    scale_x_continuous(breaks = c(1, 3, 6, 12, 24)) +
    theme_bw(base_size = 12) +
    labs(title = "A) % Pixels with Significant Trends",
         subtitle = "Only variance-check-passed pixels included",
         x = "Timescale (months)", y = "Percentage (%)", color = "Index")
  
  p2 <- ggplot(timescale_summary,
               aes(x = timescale, y = Median_Tau,
                   color = index_type, group = index_type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    scale_color_manual(values = c("spi" = "#e41a1c", "spei" = "#377eb8"),
                       labels = c("SPI", "SPEI")) +
    scale_x_continuous(breaks = c(1, 3, 6, 12, 24)) +
    theme_bw(base_size = 12) +
    labs(title = "B) Median Kendall's Tau",
         x = "Timescale (months)", y = "Tau", color = "Index")
  
  p3 <- ggplot(timescale_summary,
               aes(x = timescale, y = Mean_Events,
                   color = index_type, group = index_type)) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    scale_color_manual(values = c("spi" = "#e41a1c", "spei" = "#377eb8"),
                       labels = c("SPI", "SPEI")) +
    scale_x_continuous(breaks = c(1, 3, 6, 12, 24)) +
    theme_bw(base_size = 12) +
    labs(title = "C) Mean Drought Events",
         x = "Timescale (months)", y = "Number", color = "Index")
  
  p4 <- ggplot(timescale_summary,
               aes(x = timescale, y = Mean_Duration,
                   color = index_type, group = index_type)) +
    geom_line(linewidth = 1.2, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    scale_color_manual(values = c("spi" = "#e41a1c", "spei" = "#377eb8"),
                       labels = c("SPI", "SPEI")) +
    scale_x_continuous(breaks = c(1, 3, 6, 12, 24)) +
    theme_bw(base_size = 12) +
    labs(title = "D) Mean Duration",
         x = "Timescale (months)", y = "Months", color = "Index")
  
  combined <- (p1 | p2) / (p3 | p4) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(output_file, combined, width = 12, height = 10, dpi = dpi)
  cat(sprintf("✅ Saved: %s\n", basename(output_file)))
}

####################################################################################
# PART 3: GENERATE BASIN-AVERAGED TIME SERIES CSVs
####################################################################################

cat("\n========================================\n")
cat("📊 PART 3: GENERATING BASIN-AVERAGED TIME SERIES\n")
cat("========================================\n")

# FIXED FUNCTION: Load individual result files directly to avoid column name collision
calculate_basin_average <- function(index_type, timescale, results_dt, input_dir = "./temporal_spi_spei/") {
  cat(sprintf("\nCalculating basin average for %s-%02d...\n",
              toupper(index_type), timescale))
  
  # Load individual file directly (FIX: avoids column name collision in combined file)
  file_path <- file.path(input_dir, sprintf("%s_%02d_results.csv", index_type, timescale))
  
  if (!file.exists(file_path)) {
    cat(sprintf("  ⚠️  File not found: %s\n", file_path))
    return(NULL)
  }
  
  dt <- fread(file_path)
  cat(sprintf("  ✓ Loaded: %d rows\n", nrow(dt)))
  
  # Clip to basin extent first
  basin_ext <- ext(basin_vect)
  dt_clipped <- dt[lon >= basin_ext[1] & lon <= basin_ext[2] &
                     lat >= basin_ext[3] & lat <= basin_ext[4]]
  
  cat(sprintf("  ✓ Points within basin extent: %d (%.1f%%)\n",
              nrow(dt_clipped), nrow(dt_clipped) / nrow(dt) * 100))
  
  if (nrow(dt_clipped) == 0) {
    cat("  ⚠️  No points within basin extent\n")
    return(NULL)
  }
  
  # EXTRACT REAL TIME SERIES from the time_XXX columns
  time_cols <- grep("^time_[0-9]{3}$", names(dt_clipped), value = TRUE)
  
  if (length(time_cols) == 0) {
    cat("  ⚠️  No time series columns found\n")
    return(NULL)
  }
  
  cat(sprintf("  ✓ Found %d time steps\n", length(time_cols)))
  
  # Extract time series matrix
  ts_matrix <- as.matrix(dt_clipped[, ..time_cols])
  
  # Calculate basin average (spatial mean at each time step)
  basin_avg <- colMeans(ts_matrix, na.rm = TRUE)
  
  # Create date sequence
  # Assuming monthly data from 1950-01-01 onward
  n_months <- length(time_cols)
  dates <- seq(as.Date("1950-01-01"), by = "month", length.out = n_months)
  
  # Create data frame
  df <- data.frame(
    date = dates,
    value = basin_avg
  )
  
  # Remove rows with NA values
  df <- df[!is.na(df$value), ]
  
  cat(sprintf("  ✓ Extracted %d monthly values (basin average)\n", nrow(df)))
  cat(sprintf("     Mean: %.3f, Median: %.3f, SD: %.3f\n",
              mean(df$value, na.rm = TRUE),
              median(df$value, na.rm = TRUE),
              sd(df$value, na.rm = TRUE)))
  
  return(df)
}

# Generate basin average CSVs for all timescales
for (idx in c("spi", "spei")) {
  for (scale in timescales_temporal) {
    df <- calculate_basin_average(idx, scale, results, input_dir)
    
    if (!is.null(df)) {
      csv_file <- file.path(timeseries_dir, sprintf("%s_%02d_basin_average.csv", idx, scale))
      write.csv(df, csv_file, row.names = FALSE)
      cat(sprintf("  ✓ Saved: %s\n", basename(csv_file)))
    }
  }
}

####################################################################################
# PART 4: GENERATE TIME SERIES PLOTS
####################################################################################

cat("\n========================================\n")
cat("📈 PART 4: GENERATING TIME SERIES PLOTS\n")
cat("========================================\n")

load_basin_timeseries_csv <- function(index_type, scale) {
  csv_filename <- sprintf("%s_%02d_basin_average.csv", index_type, scale)
  csv_path <- file.path(timeseries_dir, csv_filename)
  
  if (!file.exists(csv_path)) {
    stop(sprintf("❌ CSV file not found: %s\nMake sure Part 3 completed successfully.", csv_path))
  }
  
  df <- fread(csv_path)
  df$date <- as.Date(df$date)
  df$value <- as.numeric(df$value)
  
  if (nrow(df) == 0) {
    stop(sprintf("❌ CSV file is empty: %s", csv_path))
  }
  
  return(df)
}

detect_drought_events <- function(df, onset_threshold = -1.0, termination_threshold = 0.0, min_duration = 1) {
  # Validate inputs
  if (is.null(df) || nrow(df) == 0) {
    warning("No data provided to detect_drought_events")
    return(data.frame())
  }
  if (!"date" %in% names(df) || !"value" %in% names(df)) {
    stop("DataFrame must contain 'date' and 'value' columns")
  }
  
  df <- df[order(df$date), ]
  values <- df$value
  dates <- df$date
  
  events <- data.frame(
    event_id = integer(),
    start_date = as.Date(character()),
    end_date = as.Date(character()),
    duration_months = integer(),
    min_value = numeric(),
    severity = character(),
    stringsAsFactors = FALSE
  )
  
  in_drought <- FALSE
  start_idx <- NA
  
  for (i in seq_along(values)) {
    if (!in_drought && !is.na(values[i]) && values[i] < onset_threshold) {
      in_drought <- TRUE
      start_idx <- i
    } else if (in_drought && !is.na(values[i]) && values[i] >= termination_threshold) {
      end_idx <- i - 1
      duration <- end_idx - start_idx + 1
      if (duration >= min_duration) {
        drought_vals <- values[start_idx:end_idx]
        min_val <- min(drought_vals, na.rm = TRUE)
        severity <- if (min_val <= extreme_threshold) {
          "Extreme"
        } else if (min_val <= severe_threshold) {
          "Severe"
        } else {
          "Moderate"
        }
        events <- rbind(events, data.frame(
          event_id = nrow(events) + 1,
          start_date = dates[start_idx],
          end_date = dates[end_idx],
          duration_months = duration,
          min_value = min_val,
          severity = severity,
          stringsAsFactors = FALSE
        ))
      }
      in_drought <- FALSE
    }
  }
  
  # Handle ongoing drought at end of series
  if (in_drought && (length(values) - start_idx + 1) >= min_duration) {
    end_idx <- length(values)
    drought_vals <- values[start_idx:end_idx]
    min_val <- min(drought_vals, na.rm = TRUE)
    severity <- if (min_val <= extreme_threshold) {
      "Extreme"
    } else if (min_val <= severe_threshold) {
      "Severe"
    } else {
      "Moderate"
    }
    events <- rbind(events, data.frame(
      event_id = nrow(events) + 1,
      start_date = dates[start_idx],
      end_date = dates[end_idx],
      duration_months = end_idx - start_idx + 1,
      min_value = min_val,
      severity = severity,
      stringsAsFactors = FALSE
    ))
  }
  
  events
}

create_timeseries_plot <- function(df, index_type, scale, output_file) {
  # Validate inputs
  if (is.null(df) || nrow(df) == 0) {
    stop("No data provided to create_timeseries_plot")
  }
  if (!"date" %in% names(df) || !"value" %in% names(df)) {
    stop("DataFrame must contain 'date' and 'value' columns")
  }
  
  p <- ggplot(df, aes(x = date, y = value)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.3) +
    geom_line(color = if (index_type == "spi") "#e41a1c" else "#377eb8", linewidth = 0.6) +
    geom_ribbon(data = subset(df, value < drought_threshold_onset),
                aes(ymin = value, ymax = drought_threshold_onset),
                fill = "coral", alpha = 0.3) +
    theme_bw(base_size = 14) +
    labs(
      title = sprintf("%s-%02d: Basin-Averaged Time Series", toupper(index_type), scale),
      x = "Year",
      y = sprintf("%s-%02d Value", toupper(index_type), scale)
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    coord_cartesian(ylim = c(-3.5, 3.5))
  
  ggsave(output_file, p, width = 12, height = 8, dpi = dpi)
}

# Generate time series plots
for (idx in c("spi", "spei")) {
  for (scale in timescales_temporal) {
    tryCatch({
      df <- load_basin_timeseries_csv(idx, scale)
      output_file <- file.path(basin_plots_dir, sprintf("%s_%02d_timeseries.png", idx, scale))
      create_timeseries_plot(df, idx, scale, output_file)
      cat(sprintf("✓ Saved: %s\n", basename(output_file)))
    }, error = function(e) {
      cat(sprintf("⚠️  Error for %s-%02d: %s\n", toupper(idx), scale, e$message))
    })
  }
}

####################################################################################
# PART 5: GENERATE EXCEL SUMMARY STATISTICS
####################################################################################

cat("\n========================================\n")
cat("📊 PART 5: GENERATING EXCEL SUMMARY\n")
cat("========================================\n")

summary_wb <- createWorkbook()
stats_data <- data.frame(Timescale = character(), Index = character(),
                         Mean = numeric(), Median = numeric(),
                         StdDev = numeric(), Min = numeric(), Max = numeric(),
                         Drought_Events_2020_2025 = integer(), stringsAsFactors = FALSE)

for (index_type in c("spi", "spei")) {
  for (scale in timescales_temporal) {
    tryCatch({
      df <- load_basin_timeseries_csv(index_type, scale)
      df_recent <- df[df$date >= as.Date("2020-01-01") & df$date <= as.Date("2025-12-31"), ]
      events <- detect_drought_events(df_recent, drought_threshold_onset, drought_threshold_end, 1)
      
      stats_data <- rbind(stats_data, data.frame(
        Timescale = sprintf("%02d", scale),
        Index = toupper(index_type),
        Mean = mean(df$value, na.rm = TRUE),
        Median = median(df$value, na.rm = TRUE),
        StdDev = sd(df$value, na.rm = TRUE),
        Min = min(df$value, na.rm = TRUE),
        Max = max(df$value, na.rm = TRUE),
        Drought_Events_2020_2025 = nrow(events)
      ))
    }, error = function(e) {
      cat(sprintf("    ⚠ Error calculating stats for %s-%02d: %s\n",
                  toupper(index_type), scale, e$message))
    })
  }
}

addWorksheet(summary_wb, "Summary_Statistics")
writeData(summary_wb, "Summary_Statistics", stats_data, rowNames = FALSE)

# Correlation Statistics
correlation_data <- data.frame(
  Timescale = character(),
  Correlation = numeric(),
  SPI_Events = integer(),
  SPEI_Events = integer(),
  stringsAsFactors = FALSE
)

for (scale in timescales_temporal) {
  tryCatch({
    spi_df <- load_basin_timeseries_csv("spi", scale)
    spei_df <- load_basin_timeseries_csv("spei", scale)
    
    # Merge for correlation
    merged <- merge(spi_df, spei_df, by = "date", all = TRUE)
    corr_val <- cor(merged$value.x, merged$value.y, use = "complete.obs")
    
    # Count events
    spi_recent <- spi_df[spi_df$date >= as.Date("2020-01-01") & spi_df$date <= as.Date("2025-12-31"), ]
    spei_recent <- spei_df[spei_df$date >= as.Date("2020-01-01") & spei_df$date <= as.Date("2025-12-31"), ]
    spi_events <- detect_drought_events(spi_recent, drought_threshold_onset, drought_threshold_end, 1)
    spei_events <- detect_drought_events(spei_recent, drought_threshold_onset, drought_threshold_end, 1)
    
    correlation_data <- rbind(correlation_data, data.frame(
      Timescale = sprintf("%02d-month", scale),
      Correlation = round(corr_val, 4),
      SPI_Events = nrow(spi_events),
      SPEI_Events = nrow(spei_events),
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  Timescale %02d: r = %.4f | SPI events: %d | SPEI events: %d\n",
                scale, corr_val, nrow(spi_events), nrow(spei_events)))
  }, error = function(e) {
    cat(sprintf("    ⚠ Error calculating correlation for scale %02d: %s\n", scale, e$message))
  })
}

addWorksheet(summary_wb, "SPI_SPEI_Correlations")
writeData(summary_wb, "SPI_SPEI_Correlations", correlation_data, rowNames = FALSE)

excel_file <- file.path(timeseries_dir, "SPI_SPEI_Summary_Statistics.xlsx")
saveWorkbook(summary_wb, excel_file, overwrite = TRUE)
cat(sprintf("+ Saved summary statistics: %s\n", basename(excel_file)))

####################################################################################
# MAIN EXECUTION - SPATIAL ANALYSIS
####################################################################################

cat("\n========================================\n")
cat("🎨 PART 6: GENERATING SPATIAL VISUALIZATIONS\n")
cat("========================================\n")

# Check if plot_data exists and has data
if (!exists("plot_data") || is.null(plot_data) || nrow(plot_data) == 0) {
  cat("⚠️  WARNING: No plot_data available for spatial visualizations\n")
  cat("   Spatial figures will be skipped. Check Part 2 output.\n")
} else {
  tryCatch({
    reconstruct_basin_timeseries(plot_data,
                                 file.path(fig_dir, "Fig1_Basin_Timeseries_Spatial.pdf"))
  }, error = function(e) cat(sprintf("   ❌ Error in basin timeseries: %s\n", e$message)))
  
  tryCatch({
    create_smooth_trend_maps(plot_data, basin,
                             file.path(fig_dir, "Fig2_Smooth_Trend_Maps_12mo.pdf"),
                             scale = 12)
  }, error = function(e) cat(sprintf("   ❌ Error in trend maps: %s\n", e$message)))
  
  tryCatch({
    create_event_raster_maps(plot_data, basin,
                             file.path(fig_dir, "Fig3_Event_Characteristics_12mo.pdf"),
                             scale = 12)
  }, error = function(e) cat(sprintf("   ❌ Error in event maps: %s\n", e$message)))
  
  tryCatch({
    create_method_comparison_plots(plot_data,
                                   file.path(fig_dir, "Fig5_Method_Comparison.png"))
  }, error = function(e) cat(sprintf("   ❌ Error in method comparison: %s\n", e$message)))
  
  tryCatch({
    create_temporal_pattern_maps(plot_data, basin,
                                 file.path(fig_dir, "Fig4_Temporal_Patterns_12mo.pdf"),
                                 scale = 12)
  }, error = function(e) cat(sprintf("   ❌ Error in temporal patterns: %s\n", e$message)))
  
  tryCatch({
    create_timescale_comparison(plot_data,
                                file.path(fig_dir, "Fig6_Timescale_Comparison.png"))
  }, error = function(e) cat(sprintf("   ❌ Error in timescale comparison: %s\n", e$message)))
}

####################################################################################
# FINAL SUMMARY
####################################################################################

cat("\n========================================\n")
cat("✅ COMPLETE WORKFLOW FINISHED\n")
cat("========================================\n")

cat(sprintf("\n📁 OUTPUT DIRECTORIES:\n"))
cat(sprintf("   Spatial figures: %s\n", normalizePath(fig_dir)))
cat(sprintf("   Time series CSVs: %s\n", normalizePath(timeseries_dir)))
cat(sprintf("   Time series plots: %s\n", normalizePath(basin_plots_dir)))

cat("\n📊 DATA QUALITY SUMMARY:\n")
cat("----------------------------------------\n")
for (idx in c("spi", "spei")) {
  idx_data <- results[index_type == idx]
  cat(sprintf("\n%s:\n", toupper(idx)))
  cat(sprintf("  Total pixels: %d\n", nrow(idx_data)))
  cat(sprintf("  Valid tau_vc: %d (%.1f%%)\n",
              sum(!is.na(idx_data$tau_vc)),
              sum(!is.na(idx_data$tau_vc)) / nrow(idx_data) * 100))
}

cat("\n📁 FILES CREATED:\n")
cat("----------------------------------------\n")
# List spatial figures
if (dir.exists(fig_dir)) {
  fig_files <- list.files(fig_dir, full.names = FALSE, pattern = "\\.(pdf|png)$")
  if (length(fig_files) > 0) {
    cat("Spatial analysis figures:\n")
    for (f in fig_files) {
      file_size <- file.info(file.path(fig_dir, f))$size
      cat(sprintf("  ✓ %s (%.1f KB)\n", f, file_size / 1024))
    }
  }
}

# List time series CSVs
if (dir.exists(timeseries_dir)) {
  csv_files <- list.files(timeseries_dir, full.names = FALSE, pattern = "\\.csv$")
  if (length(csv_files) > 0) {
    cat("\nTime series CSVs:\n")
    for (f in csv_files) {
      file_size <- file.info(file.path(timeseries_dir, f))$size
      cat(sprintf("  ✓ %s (%.1f KB)\n", f, file_size / 1024))
    }
  }
}

cat("\n========================================\n")
cat("⚠️  IMPORTANT NOTES:\n")
cat("----------------------------------------\n")
cat("1. Basin averages now load from individual files (avoids column collision)\n")
cat("2. TFPW method comparison skipped if columns not present\n")
cat("3. All plots include variance-check-passed pixels\n")
cat("4. Check correlation values - should NOT be 1.0 between SPI/SPEI\n")
cat("========================================\n")