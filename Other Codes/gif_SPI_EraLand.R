# ============================================================================
# Create Animated GIFs for SPI Indices - ENHANCED AESTHETICS VERSION 
# ============================================================================
# This script creates GIF animations for SPEI time series with basin overlay.
#
# ============================================================================
library(terra)
library(magick)
library(ggplot2)
library(viridis)
library(lubridate)
library(scales)  # for squish function

# ===================== CONFIGURATION =====================
spi_scales <- c(3, 6, 12)

setwd("D:/Nechako_Drought/")
spi_data_dir <- "spi_results/"
basin_shp_path <- "Spatial/nechakoBound_dissolve.shp"
output_dir <- "spi_gifs/"
dir.create(output_dir, showWarnings = FALSE)

# ===================== FUNCTION: CREATE SINGLE SPI FRAME =====================
create_spi_frame <- function(spi_raster, time_index, basin_shp, scale, date_label = NULL) {
  # Extract the specific time layer
  spi_layer <- spi_raster[[time_index]]
  
  # Reproject basin to match SPI raster CRS if needed
  if (!same.crs(basin_shp, spi_layer)) {
    basin_shp <- project(basin_shp, crs(spi_layer))
    cat("    Reprojected basin to match SPI CRS\n")
  }
  
  # Convert to data frame for ggplot
  spi_df <- as.data.frame(spi_layer, xy = TRUE, na.rm = FALSE)
  colnames(spi_df) <- c("x", "y", "spi_value")
  
  # Keep rows with valid coordinates (keep NA values for proper rendering)
  spi_df <- spi_df[!is.na(spi_df$x) & !is.na(spi_df$y), ]
  
  # Convert basin shape to data frame for ggplot
  basin_coords <- geom(basin_shp)
  basin_df <- as.data.frame(basin_coords)
  
  # Get extent for proper axis labels
  ext <- ext(spi_layer)
  
  # Create the plot with SPEI-style aesthetics
  p <- ggplot() +
    # Plot SPI raster
    geom_raster(data = spi_df, aes(x = x, y = y, fill = spi_value)) +
    
    # Add basin boundaries (matching SPEI style)
    geom_path(data = basin_df, 
              aes(x = x, y = y, group = interaction(geom, part)),
              color = "darkgreen", linewidth = 1.2) +
    
    # Color scale for SPI - EXPANDED LIMITS and SQUISH to capture extreme values
    scale_fill_gradientn(
      colors = c("#8B0000", "#FF0000", "#FFA500", "#FFFF00", "#ADFF2F", "#00FA9A", "#1E90FF"),
      limits = c(-4, 4),  # EXPANDED from -3,3 to capture extreme drought values
      breaks = seq(-4, 4, 1),
      name = "SPI Value",
      na.value = "grey90",
      oob = scales::squish  # CRITICAL: squish values outside limits instead of dropping them
    ) +
    
    # Title and labels (matching SPEI style)
    labs(
      title = paste0("SPI-", scale, " Drought Index"),
      subtitle = if (!is.null(date_label)) date_label else "",
      x = "Longitude",
      y = "Latitude",
      caption = "Negative values = Drought conditions"
    ) +
    
    # Theme customization (matching SPEI style)
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5, lineheight = 1.2),
      plot.caption = element_text(size = 10, hjust = 0.5, color = "gray40"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "grey90", linewidth = 0.2)
    ) +
    
    # Coordinate system with proper limits
    coord_equal(xlim = c(ext[1], ext[2]), ylim = c(ext[3], ext[4])) +
    
    # Legend position
    theme(legend.position = "right")
  
  return(p)
}

# ===================== FUNCTION: EXTRACT DATES FROM NETCDF =====================
extract_dates_from_nc <- function(nc_file, n_layers = NULL) {
  if (is.null(n_layers)) {
    n_layers <- nlyr(rast(nc_file))
  }
  
  # Try multiple methods to extract dates
  tryCatch({
    # Method 1: Try using terra's time()
    r <- rast(nc_file)
    dates <- time(r)
    
    if (!is.null(dates) && length(dates) > 0 && length(dates) == n_layers) {
      # Try to convert to Date
      dates_converted <- tryCatch({
        as.Date(dates)
      }, error = function(e) {
        as.Date(dates, origin = "1970-01-01")
      })
      
      # Check if dates are reasonable (between 1900 and 2030)
      if (all(!is.na(dates_converted)) && 
          all(dates_converted >= as.Date("1900-01-01")) &&
          all(dates_converted <= as.Date("2030-01-01"))) {
        cat("  Extracted dates from NetCDF time dimension\n")
        return(dates_converted)
      }
    }
    
    # Method 2: Try reading time variable directly from NetCDF
    if (requireNamespace("ncdf4", quietly = TRUE)) {
      library(ncdf4)
      nc <- nc_open(nc_file)
      
      # Check for time variable
      if ("time" %in% names(nc$dim)) {
        time_vals <- ncvar_get(nc, "time")
        time_units <- ncatt_get(nc, "time", "units")$value
        
        # Parse units (e.g., "days since 1900-01-01" or "months since 1900-01-01")
        if (grepl("since", time_units)) {
          origin_str <- sub(".*since ", "", time_units)
          origin_date <- as.Date(origin_str)
          
          if (grepl("^days", time_units)) {
            dates <- origin_date + time_vals
          } else if (grepl("^months", time_units)) {
            # For months, add months properly
            dates <- origin_date + months(time_vals - 1)
          }
          
          nc_close(nc)
          
          # Validate dates
          if (all(dates >= as.Date("1900-01-01")) && 
              all(dates <= as.Date("2030-01-01"))) {
            cat("  Extracted dates from NetCDF units attribute\n")
            return(as.Date(dates))
          }
        }
      }
      nc_close(nc)
    }
    
  }, error = function(e) {
    cat("  Date extraction failed:", e$message, "\n")
  })
  
  # Method 3: Generate sequential monthly dates
  cat("  Generating sequential dates starting from Jan 1950\n")
  start_date <- as.Date("1950-01-01")
  dates <- seq(start_date, by = "month", length.out = n_layers)
  
  # Validate that we don't exceed Dec 2025
  end_date <- dates[length(dates)]
  if (end_date > as.Date("2025-12-31")) {
    cat("  WARNING: Generated dates exceed Dec 2025. Adjusting start date...\n")
    all_dates <- seq(as.Date("1950-01-01"), as.Date("2025-12-31"), by = "month")
    if (length(all_dates) >= n_layers) {
      dates <- tail(all_dates, n_layers)
    } 
  }
  
  return(dates)
}

# ===================== MAIN PROCESSING =====================
# Load basin shapefile
basin_shp <- vect(basin_shp_path)

# Process each SPI scale
for (scale in spi_scales) {
  cat("Creating GIF for SPI-", scale, "...\n", sep = "")
  
  # Try different possible filenames (including leading zeros)
  possible_files <- c(
    paste0(spi_data_dir, "spi_", sprintf("%02d", scale), "month.nc"),
    paste0(spi_data_dir, "spi_", scale, "month.nc"),
    paste0(spi_data_dir, "spi_", scale, ".nc"),
    paste0(spi_data_dir, "SPI_", scale, ".nc")
  )
  
  spi_file <- NULL
  for (f in possible_files) {
    if (file.exists(f)) {
      spi_file <- f
      break
    }
  }
  
  if (is.null(spi_file)) {
    warning("Could not find SPI-", scale, " file. Tried:\n  ",
            paste(possible_files, collapse = "\n  "))
    next
  }
  
  cat("  Loading:", spi_file, "\n")
  spi_data <- rast(spi_file)
  n_layers <- nlyr(spi_data)
  
  # Print CRS info
  cat("  SPI CRS:", as.character(crs(spi_data)), "\n")
  cat("  Basin CRS:", as.character(crs(basin_shp)), "\n")
  
  # Extract dates
  dates <- extract_dates_from_nc(spi_file, n_layers)
  
  cat("  Date range:", format(dates[1], "%b %Y"), "to", 
      format(dates[length(dates)], "%b %Y"), "\n")
  cat("  Total frames:", n_layers, "\n")
  
  # Filter to start from Jan 2009
  start_date <- as.Date("2009-01-01")
  valid_dates <- dates >= start_date
  if (any(valid_dates)) {
    start_index <- which(valid_dates)[1]
    cat("  Starting animation from:", format(dates[start_index], "%b %Y"), 
        " (frame", start_index, ")\n")
    n_layers <- sum(valid_dates)
  } else {
    cat("  Warning: No data from 2009 onwards. Using all available data.\n")
    start_index <- 1
  }
  
  # Create frames
  frames <- list()
  frame_counter <- 0
  
  for (i in start_index:min(nlyr(spi_data), start_index + n_layers - 1)) {
    frame_counter <- frame_counter + 1
    if (frame_counter %% 20 == 0 || frame_counter == 1) {
      cat("  Processing frame", frame_counter, "of", n_layers, "\n")
    }
    
    # Create date label in format "Beginning Month: Jan 1950"
    date_label <- NULL
    if (!is.null(dates) && i <= length(dates)) {
      date_label <- paste("Beginning Month:", format(dates[i], "%b %Y"))
    } else {
      date_label <- paste("Beginning Month: Jan 2009")
    }
    
    # Create the frame
    frame_plot <- create_spi_frame(spi_data, i, basin_shp, scale, date_label)
    
    # Save frame as temporary PNG (matching SPEI dimensions)
    temp_file <- tempfile(fileext = ".png")
    ggsave(temp_file, frame_plot, width = 10, height = 7, dpi = 100)
    
    # Read the PNG into magick image
    frames[[frame_counter]] <- image_read(temp_file)
    
    # Clean up temp file
    file.remove(temp_file)
  }
  
  cat("  Creating GIF...\n")
  
  # Combine frames into animated GIF
  gif <- image_animate(
    image_join(frames),
    fps = 4,
    optimize = TRUE
  )
  
  # Save the GIF
  gif_file <- paste0(output_dir, "SPI_", scale, "_animation.gif")
  image_write(gif, gif_file)
  
  cat("  GIF saved as:", gif_file, "\n")
  cat("  File size:", round(file.size(gif_file) / 1024^2, 2), "MB\n")
  
  # Optional: Create yearly summary if many frames
  if (length(frames) > 24 && !is.null(dates)) {
    # Get June indices (month 6) for frames from start_index onwards
    dates_subset <- dates[start_index:length(dates)]
    june_indices_subset <- which(as.numeric(format(dates_subset, "%m")) == 6)
    
    if (length(june_indices_subset) > 1) {
      cat("  Creating yearly summary (", length(june_indices_subset), " June frames)...\n")
      yearly_frames <- frames[june_indices_subset]
      yearly_gif <- image_animate(image_join(yearly_frames), fps = 1, optimize = TRUE)
      yearly_file <- paste0(output_dir, "SPI_", scale, "_yearly_animation.gif")
      image_write(yearly_gif, yearly_file)
      cat("  Yearly summary GIF saved as:", yearly_file, "\n")
    }
  }
  
  # Clean up
  rm(spi_data, frames, gif)
  gc()
}

cat("\n==========================================\n")
cat("SPI GIF creation complete!\n")
cat("All animations saved in:", output_dir, "\n")

# ===================== OPTIONAL: CREATE COMPOSITE GIF =====================
cat("\nCreating composite comparison GIF...\n")

composite_scales <- c(3, 6, 12)
composite_files <- list()

# Find existing files for composite
for (scale in composite_scales) {
  possible_files <- c(
    paste0(spi_data_dir, "spi_", sprintf("%02d", scale), "month.nc"),
    paste0(spi_data_dir, "spi_", scale, "month.nc"),
    paste0(spi_data_dir, "spi_", scale, ".nc"),
    paste0(spi_data_dir, "SPI_", scale, ".nc")
  )
  
  for (f in possible_files) {
    if (file.exists(f)) {
      composite_files[[as.character(scale)]] <- f
      break
    }
  }
}

if (length(composite_files) >= 2) {
  # Get minimum number of layers
  layer_counts <- sapply(composite_files, function(f) nlyr(rast(f)))
  min_layers <- min(layer_counts)
  
  # Install patchwork if needed
  if (!require("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
    library(patchwork)
  }
  
  composite_frames <- list()
  
  # Load all datasets and extract dates
  composite_data <- list()
  composite_dates <- list()
  
  for (scale in names(composite_files)) {
    composite_data[[scale]] <- rast(composite_files[[scale]])
    composite_dates[[scale]] <- extract_dates_from_nc(composite_files[[scale]], 
                                                      nlyr(composite_data[[scale]]))
  }
  
  frame_counter <- 0
  for (i in 1:min_layers) {
    frame_counter <- frame_counter + 1
    if (frame_counter %% 10 == 0 || frame_counter == 1) {
      cat("  Processing composite frame", frame_counter, "of", min_layers, "\n")
    }
    
    plot_list <- list()
    
    for (scale in names(composite_files)) {
      # Get date label
      date_label <- NULL
      if (!is.null(composite_dates[[scale]]) && i <= length(composite_dates[[scale]])) {
        date_label <- paste0("SPI-", scale, " | ", format(composite_dates[[scale]][i], "%b %Y"))
      } else {
        date_label <- paste0("SPI-", scale, " | Jan 1950")
      }
      
      p <- create_spi_frame(composite_data[[scale]], i, basin_shp, scale, date_label)
      
      plot_list[[scale]] <- p
    }
    
    composite_plot <- wrap_plots(plot_list, ncol = 2) +
      plot_annotation(
        title = "SPI Drought Index Comparison",
        theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
      )
    
    temp_file <- tempfile(fileext = ".png")
    ggsave(temp_file, composite_plot, width = 14, height = 9, dpi = 80)
    composite_frames[[frame_counter]] <- image_read(temp_file)
    file.remove(temp_file)
  }
  
  if (length(composite_frames) > 0) {
    cat("  Creating composite GIF...\n")
    composite_gif <- image_animate(image_join(composite_frames), fps = 4, optimize = TRUE)
    composite_file <- paste0(output_dir, "SPI_comparison_animation.gif")
    image_write(composite_gif, composite_file)
    cat("  Composite GIF saved as:", composite_file, "\n")
  }
} else {
  cat("  Not enough SPI files found for composite (need at least 2)\n")
}

cat("\nAll SPI animations completed successfully!\n")