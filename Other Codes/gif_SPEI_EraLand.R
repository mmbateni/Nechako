# ============================================================================
# Create Animated GIFs for SPEI Indices - ENHANCED AESTHETICS VERSION
# ============================================================================
# This script creates GIF animations for SPEI time series with basin overlay.
#
# ============================================================================

# --------------------- Libraries ---------------------
library(terra)      # raster/vector processing
library(magick)     # GIF creation
library(ggplot2)    # plotting
library(lubridate)  # date arithmetic
library(scales)     # for squish function

# --------------------- Configuration ---------------------
# SPEI timescales to visualize
spei_scales <- c(3, 6, 12)

# Working directory
setwd("D:/Nechako_Drought/")

# Inputs
spei_data_dir   <- "spei_results/"
basin_shp_path  <- "Spatial/nechakoBound_dissolve.shp"

# Outputs
output_dir <- "spei_gifs/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------- Helper: Extract dates from NetCDF ---------------------
extract_dates_from_nc <- function(nc_file, n_layers = NULL) {
  if (is.null(n_layers)) {
    n_layers <- nlyr(rast(nc_file))
  }
  
  # Method 1: terra time() if present
  try({
    r <- rast(nc_file)
    dates <- time(r)
    if (!is.null(dates) && length(dates) == n_layers) {
      dates_converted <- tryCatch(as.Date(dates), error = function(e) as.Date(dates, origin = "1970-01-01"))
      if (all(!is.na(dates_converted)) &&
          min(dates_converted) >= as.Date("1900-01-01") &&
          max(dates_converted) <= as.Date("2025-01-01")) {
        cat(" Extracted dates from NetCDF time dimension (terra::time)\n")
        return(dates_converted)
      }
    }
  }, silent = TRUE)
  
  # Method 2: read time variable and units via ncdf4 (if available)
  if (requireNamespace("ncdf4", quietly = TRUE)) {
    try({
      nc <- ncdf4::nc_open(nc_file)
      on.exit(ncdf4::nc_close(nc), add = TRUE)
      
      if ("time" %in% names(nc$dim) || "time" %in% names(nc$var)) {
        # time as a variable is more common; try var first, then dim
        time_vals  <- tryCatch(ncdf4::ncvar_get(nc, "time"), error = function(e) nc$dim$time$vals)
        time_units <- tryCatch(ncdf4::ncatt_get(nc, "time", "units")$value, error = function(e) NA_character_)
        
        if (!is.na(time_units) && grepl("since", time_units, ignore.case = TRUE)) {
          origin_str  <- sub(".*since\\s+", "", time_units)
          origin_date <- as.Date(origin_str)
          
          if (grepl("^days", time_units, ignore.case = TRUE)) {
            dates <- origin_date + as.numeric(time_vals)
          } else if (grepl("^months", time_units, ignore.case = TRUE)) {
            # month-aware addition
            dates <- origin_date %m+% months(as.numeric(time_vals))
          } else if (grepl("^hours", time_units, ignore.case = TRUE)) {
            dates <- as.Date(as.POSIXct(origin_date) + as.numeric(time_vals) * 3600)
          } else {
            dates <- NULL
          }
          
          if (!is.null(dates) && length(dates) == n_layers &&
              min(dates) >= as.Date("1900-01-01") &&
              max(dates) <= as.Date("2025-01-01")) {
            cat(" Extracted dates from NetCDF time units (ncdf4)\n")
            return(as.Date(dates))
          }
        }
      }
    }, silent = TRUE)
  }
  
  # Method 3: fallback (monthly sequence)
  cat(" Generating sequential monthly dates starting from Jan 1954 (fallback)\n")
  start_date <- as.Date("1954-01-01")
  dates <- seq(start_date, by = "month", length.out = n_layers)
  
  # If end exceeds Dec 2024, anchor to Dec 2024 and go backwards month-wise
  if (dates[length(dates)] > as.Date("2024-12-31")) {
    end_date <- as.Date("2024-12-01")
    dates <- seq(end_date, by = "-1 month", length.out = n_layers)
    dates <- rev(dates)
  }
  
  return(dates)
}

# --------------------- Helper: Create a single SPEI frame ---------------------
create_spei_frame <- function(spei_raster, time_index, basin_shp, scale, date_label = NULL) {
  spei_layer <- spei_raster[[time_index]]
  
  # Match CRS
  if (!same.crs(basin_shp, spei_layer)) {
    basin_shp <- project(basin_shp, crs(spei_layer))
    cat(" Reprojected basin to match SPEI CRS\n")
  }
  
  # Raster -> data frame for ggplot
  spei_df <- as.data.frame(spei_layer, xy = TRUE, na.rm = FALSE)
  colnames(spei_df) <- c("x", "y", "spei_value")
  
  # Keep NA values (so they map to na.value color), only drop missing coords
  spei_df <- spei_df[!is.na(spei_df$x) & !is.na(spei_df$y), ]
  
  # Basin -> coordinate data frame for polygons/lines (terra::geom is reliable)
  basin_geom <- terra::geom(basin_shp)
  basin_df   <- as.data.frame(basin_geom)
  
  # Ensure coordinate column names are x/y
  if (!all(c("x","y") %in% names(basin_df))) {
    names(basin_df)[1:2] <- c("x","y")
  }
  
  # Group by available geometry columns so polygons draw correctly
  group_cols <- intersect(c("geom","part","hole"), names(basin_df))
  if (length(group_cols) == 0) {
    basin_df$group <- 1
  } else {
    basin_df$group <- interaction(basin_df[, group_cols], drop = TRUE)
  }
  
  
  ext <- ext(spei_layer)
  
  p <- ggplot() +
    geom_raster(data = spei_df, aes(x = x, y = y, fill = spei_value)) +
    geom_path(data = basin_df,
              aes(x = x, y = y, group = group),
              color = "darkgreen", linewidth = 1.2) +
    scale_fill_gradientn(
      colors = c("#8B0000", "#FF0000", "#FFA500", "#FFFF00", "#ADFF2F", "#00FA9A", "#1E90FF"),
      limits = c(-4, 4),  # EXPANDED from -3,3 to capture extreme values
      breaks = seq(-4, 4, 1),
      name = "SPEI Value",
      na.value = "grey90",
      oob = scales::squish  # CRITICAL: squish values outside limits instead of dropping them
    ) +
    labs(
      title = paste0("SPEI-", scale, " Drought Index"),
      subtitle = if (!is.null(date_label)) date_label else "",
      x = "X", y = "Y",
      caption = "Negative values = drought conditions"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5, lineheight = 1.2),
      plot.caption = element_text(size = 10, hjust = 0.5, color = "gray40"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "grey90", linewidth = 0.2),
      legend.position = "right"
    ) +
    coord_equal(xlim = c(ext[1], ext[2]), ylim = c(ext[3], ext[4]))
  
  return(p)
}

# ============================================================================
# Main processing
# ============================================================================

basin_shp <- vect(basin_shp_path)

for (scale in spei_scales) {
  cat("Creating GIF for SPEI-", scale, "...\n", sep = "")
  
  possible_files <- c(
    paste0(spei_data_dir, "spei_", sprintf("%02d", scale), "month.nc"),
    paste0(spei_data_dir, "spei_", scale, "month.nc"),
    paste0(spei_data_dir, "spei_", scale, ".nc"),
    paste0(spei_data_dir, "SPEI_", scale, ".nc")
  )
  
  spei_file <- possible_files[file.exists(possible_files)][1]
  if (is.na(spei_file)) {
    warning("Could not find SPEI-", scale, " file. Tried:\n ", paste(possible_files, collapse = "\n "))
    next
  }
  
  cat(" Loading: ", spei_file, "\n", sep = "")
  spei_data <- rast(spei_file)
  n_layers <- nlyr(spei_data)
  
  cat(" SPI CRS: ", as.character(crs(spei_data)), "\n", sep = "")
  cat(" Basin CRS: ", as.character(crs(basin_shp)), "\n", sep = "")
  
  dates <- extract_dates_from_nc(spei_file, n_layers)
  
  cat(" Date range: ", format(dates[1], "%b %Y"), " to ", format(dates[length(dates)], "%b %Y"), "\n", sep = "")
  cat(" Total frames: ", n_layers, "\n", sep = "")
  
  start_date <- as.Date("2009-01-01")
  valid_dates <- dates >= start_date
  start_index <- 1
  end_index <- n_layers
  
  if (any(valid_dates)) {
    start_index <- which(valid_dates)[1]
    cat(" Starting animation from: ", format(dates[start_index], "%b %Y"), " (frame ", start_index, ")\n", sep = "")
    n_layers <- sum(valid_dates)
    end_index <- length(dates)
  } else {
    cat(" Warning: No data from 2009 onwards. Using all available data.\n")
  }
  
  frames <- list()
  frame_counter <- 0
  
  for (i in start_index:end_index) {
    frame_counter <- frame_counter + 1
    if (frame_counter %% 20 == 0 || frame_counter == 1) {
      cat(" Processing frame ", frame_counter, " of ", n_layers, "\n", sep = "")
    }
    
    date_label <- NULL
    if (!is.null(dates) && i <= length(dates)) {
      date_label <- paste("Beginning Month:", format(dates[i], "%b %Y"))
    }
    
    frame_plot <- create_spei_frame(spei_data, i, basin_shp, scale, date_label)
    
    temp_file <- tempfile(fileext = ".png")
    ggsave(temp_file, frame_plot, width = 10, height = 7, dpi = 100)
    
    frames[[frame_counter]] <- image_read(temp_file)
    
    file.remove(temp_file)
  }
  
  cat(" Creating GIF...\n")
  
  gif <- image_animate(image_join(frames), fps = 4, optimize = TRUE)
  gif_file <- paste0(output_dir, "SPEI_", scale, "_animation.gif")
  image_write(gif, gif_file)
  
  cat(" GIF saved as: ", gif_file, "\n", sep = "")
  cat(" File size: ", round(file.size(gif_file) / 1024^2, 2), " MB\n", sep = "")
  
  # Optional yearly (June) summary
  if (length(frames) > 24 && !is.null(dates)) {
    dates_subset <- dates[start_index:end_index]
    june_idx <- which(as.integer(format(dates_subset, "%m")) == 6)
    
    if (length(june_idx) > 1) {
      cat(" Creating yearly summary (", length(june_idx), " June frames)...\n", sep = "")
      yearly_frames <- frames[june_idx]
      yearly_gif <- image_animate(image_join(yearly_frames), fps = 1, optimize = TRUE)
      yearly_file <- paste0(output_dir, "SPEI_", scale, "_yearly_animation.gif")
      image_write(yearly_gif, yearly_file)
      cat(" Yearly summary GIF saved as: ", yearly_file, "\n", sep = "")
    }
  }
  
  rm(spei_data, frames, gif)
  gc()
}

cat("\n==========================================\n")
cat("SPEI GIF creation complete!\n")
cat("All animations saved in: ", output_dir, "\n", sep = "")

# ============================================================================
# Optional: Create composite comparison GIF
# ============================================================================
cat("\nCreating composite comparison GIF...\n")

composite_scales <- c(3, 6, 12)
composite_files  <- list()

for (scale in composite_scales) {
  possible_files <- c(
    paste0(spei_data_dir, "spei_", sprintf("%02d", scale), "month.nc"),
    paste0(spei_data_dir, "spei_", scale, "month.nc"),
    paste0(spei_data_dir, "spei_", scale, ".nc"),
    paste0(spei_data_dir, "SPEI_", scale, ".nc")
  )
  
  f <- possible_files[file.exists(possible_files)][1]
  if (!is.na(f)) composite_files[[as.character(scale)]] <- f
}

if (length(composite_files) >= 2) {
  layer_counts <- sapply(composite_files, function(f) nlyr(rast(f)))
  min_layers   <- min(layer_counts)
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  library(patchwork)
  
  composite_data  <- lapply(composite_files, rast)
  composite_dates <- mapply(extract_dates_from_nc, composite_files, n_layers = layer_counts, SIMPLIFY = FALSE)
  
  composite_frames <- vector("list", length = min_layers)
  frame_counter <- 0
  
  for (i in 1:min_layers) {
    frame_counter <- frame_counter + 1
    
    if (frame_counter %% 10 == 0 || frame_counter == 1) {
      cat(" Processing composite frame ", frame_counter, " of ", min_layers, "\n", sep = "")
    }
    
    plot_list <- list()
    for (scale in names(composite_files)) {
      date_label <- paste0("SPEI-", scale, "\n", format(composite_dates[[scale]][i], "%b %Y"))
      plot_list[[scale]] <- create_spei_frame(composite_data[[scale]], i, basin_shp, scale, date_label)
    }
    
    composite_plot <- wrap_plots(plot_list, ncol = 2) +
      plot_annotation(
        title = "SPEI Drought Index Comparison",
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )
    
    temp_file <- tempfile(fileext = ".png")
    ggsave(temp_file, composite_plot, width = 14, height = 9, dpi = 80)
    composite_frames[[frame_counter]] <- image_read(temp_file)
    file.remove(temp_file)
  }
  
  if (length(composite_frames) > 0) {
    cat(" Creating composite GIF...\n")
    composite_gif  <- image_animate(image_join(composite_frames), fps = 4, optimize = TRUE)
    composite_file <- paste0(output_dir, "SPEI_comparison_animation.gif")
    image_write(composite_gif, composite_file)
    cat(" Composite GIF saved as: ", composite_file, "\n", sep = "")
  }
  
} else {
  cat(" Not enough SPEI files found for composite (need at least 2)\n")
}

# ============================================================================
# Drought event summary (fast basin-average series)
# ============================================================================
cat("\nCreating drought event summary...\n")

primary_scale <- 12
primary_candidates <- c(
  paste0(spei_data_dir, "spei_", sprintf("%02d", primary_scale), "month.nc"),
  paste0(spei_data_dir, "spei_", primary_scale, "month.nc"),
  paste0(spei_data_dir, "spei_", primary_scale, ".nc"),
  paste0(spei_data_dir, "SPEI_", primary_scale, ".nc")
)
primary_file <- primary_candidates[file.exists(primary_candidates)][1]

if (!is.na(primary_file)) {
  spei_data <- rast(primary_file)
  dates <- extract_dates_from_nc(primary_file, nlyr(spei_data))
  
  # Reproject basin to match raster CRS
  basin_shp_matched <- basin_shp
  if (!same.crs(basin_shp_matched, spei_data)) {
    cat(" Reprojecting basin to match SPEI CRS for drought analysis...\n")
    basin_shp_matched <- project(basin_shp, crs(spei_data))
  }
  
  # Quick overlap check
  overlap_ok <- TRUE
  tryCatch({ crop(spei_data[[1]], basin_shp_matched) }, error = function(e) { overlap_ok <<- FALSE })
  
  if (!overlap_ok) {
    cat(" WARNING: Extents do not overlap. Cannot create drought summary.\n")
  } else {
    # Fast basin mean for all layers (single polygon => one row)
    cat(" Calculating basin-average SPEI for ", nlyr(spei_data), " time steps (fast extract)...\n", sep = "")
    ex <- terra::extract(spei_data, basin_shp_matched, fun = mean, na.rm = TRUE)
    basin_series <- as.numeric(ex[1, -1])  # drop ID column
    
    # Identify drought events (SPEI < -0.5 for >= 3 consecutive months)
    drought_mask  <- basin_series < -0.5
    drought_events <- rle(drought_mask)
    
    drought_summary <- data.frame()
    start_idx <- 1
    for (j in seq_along(drought_events$lengths)) {
      if (isTRUE(drought_events$values[j]) && drought_events$lengths[j] >= 3) {
        end_idx <- start_idx + drought_events$lengths[j] - 1
        drought_summary <- rbind(
          drought_summary,
          data.frame(
            Start_Date      = format(dates[start_idx], "%Y-%m-%d"),
            End_Date        = format(dates[end_idx], "%Y-%m-%d"),
            Duration_Months = drought_events$lengths[j],
            Min_SPEI        = round(min(basin_series[start_idx:end_idx], na.rm = TRUE), 3),
            Mean_SPEI       = round(mean(basin_series[start_idx:end_idx], na.rm = TRUE), 3)
          )
        )
      }
      start_idx <- start_idx + drought_events$lengths[j]
    }
    
    if (nrow(drought_summary) > 0) {
      summary_file <- paste0(output_dir, "SPEI_", primary_scale, "_drought_events.csv")
      write.csv(drought_summary, summary_file, row.names = FALSE)
      cat(" Drought event summary saved as: ", summary_file, "\n", sep = "")
      cat("\nDrought Events Detected (SPEI-", primary_scale, "):\n", sep = "")
      print(drought_summary)
    } else {
      cat(" No significant drought events detected.\n")
    }
    
    # Time-series plot
    ts_plot <- ggplot(data.frame(Date = dates, SPEI = basin_series), aes(x = Date, y = SPEI)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = -0.5, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = -1, linetype = "dashed", color = "red") +
      labs(
        title = paste0("Basin-Average SPEI-", primary_scale, " Time Series"),
        x = "Date", y = "SPEI Value",
        caption = "Dashed lines: 0 (normal), -0.5 (moderate drought), -1 (severe drought)"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.caption = element_text(size = 9, color = "gray40"))
    
    ts_file <- paste0(output_dir, "SPEI_", primary_scale, "_timeseries.png")
    ggsave(ts_file, ts_plot, width = 12, height = 6, dpi = 150)
    cat(" Time series plot saved as: ", ts_file, "\n", sep = "")
  }
  
} else {
  cat(" SPEI file not found for drought summary. Tried:\n ", paste(primary_candidates, collapse = "\n "), "\n", sep = "")
}

cat("\nAll SPEI processing completed successfully!\n")