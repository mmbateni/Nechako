# ============================================================================
# Create Animated GIFs for SPI Indices
# ============================================================================
# This script creates GIF animations for SPI timeseries with basin overlay

# Load required libraries
library(terra)      # For raster processing
library(magick)     # For GIF creation
library(ggplot2)    # For plotting (alternative approach)
library(viridis)    # For color palettes

# ===================== CONFIGURATION =====================
# Define SPI timescales to visualize (adjust based on your needs)
spi_scales <- c(3, 6, 12)  # You can add 1, 9, 24 if needed

# Path to SPI data (from Script 2 output)
spi_data_dir <- "spi_results/"

# Path to basin shapefile
basin_shp_path <- "Spatial/nechakoBound_dissolve.shp"

# Output directory for GIFs
output_dir <- "spi_gifs/"
dir.create(output_dir, showWarnings = FALSE)

# ===================== FUNCTION: CREATE SINGLE SPI FRAME =====================
create_spi_frame <- function(spi_raster, time_index, basin_shp, scale) {
  # Creates a single frame for the GIF animation
  
  # Extract the specific time layer
  spi_layer <- spi_raster[[time_index]]
  layer_date <- time(spi_layer)
  
  # Convert to data frame for ggplot
  spi_df <- as.data.frame(spi_layer, xy = TRUE)
  colnames(spi_df) <- c("x", "y", "spi_value")
  
  # Convert basin shape to data frame for ggplot
  basin_df <- as.data.frame(basin_shp, geom = "XY")
  
  # Create the plot
  p <- ggplot() +
    # Plot SPI raster
    geom_raster(data = spi_df, aes(x = x, y = y, fill = spi_value)) +
    
    # Add basin boundaries
    geom_polygon(data = basin_df, aes(x = x, y = y, group = part),
                 fill = NA, color = "black", linewidth = 1) +
    
    # Color scale for SPI (standard drought colors)
    scale_fill_gradientn(
      colors = c("#8B0000", "#FF0000", "#FFFF00", "#FFFFFF", "#00FFFF", "#0000FF", "#00008B"),
      limits = c(-3, 3),
      breaks = seq(-3, 3, 1),
      name = "SPI Value",
      na.value = "grey90"
    ) +
    
    # Title and labels
    labs(
      title = paste0("SPI-", scale, " Drought Index"),
      subtitle = format(layer_date, "%B %Y"),
      x = "Longitude",
      y = "Latitude"
    ) +
    
    # Theme customization
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "grey90", linewidth = 0.2)
    ) +
    
    # Coordinate system
    coord_sf(crs = crs(spi_layer), default = TRUE) +
    
    # Legend position
    theme(legend.position = "right")
  
  return(p)
}

# ===================== MAIN PROCESSING =====================
# Load basin shapefile
basin_shp <- vect(basin_shp_path)

# Process each SPI scale
for (scale in spi_scales) {
  cat("Creating GIF for SPI-", scale, "...\n", sep = "")
  
  # Load SPI data for this scale
  spi_file <- paste0(spi_data_dir, "spi_", scale, "month.nc")
  
  if (!file.exists(spi_file)) {
    warning("File not found: ", spi_file, ". Skipping SPI-", scale)
    next
  }
  
  spi_data <- rast(spi_file)
  dates <- time(spi_data)
  
  # Create a list to store individual frames
  frames <- list()
  
  # Generate each frame
  for (i in 1:nlyr(spi_data)) {
    cat("\r  Processing frame", i, "of", nlyr(spi_data))
    
    # Create the frame
    frame_plot <- create_spi_frame(spi_data, i, basin_shp, scale)
    
    # Save frame as temporary PNG
    temp_file <- tempfile(fileext = ".png")
    ggsave(temp_file, frame_plot, width = 10, height = 7, dpi = 100)
    
    # Read the PNG into magick image
    frames[[i]] <- image_read(temp_file)
    
    # Clean up temp file
    file.remove(temp_file)
  }
  
  cat("\n  Creating GIF...\n")
  
  # Combine frames into animated GIF
  gif <- image_animate(
    image_join(frames),
    fps = 2,  # Frames per second (adjust as needed)
    optimize = TRUE
  )
  
  # Save the GIF
  gif_file <- paste0(output_dir, "SPI_", scale, "_animation.gif")
  image_write(gif, gif_file)
  
  cat("  GIF saved as:", gif_file, "\n")
  
  # Optional: Also create a faster/smaller preview GIF
  if (nlyr(spi_data) > 24) {  # If many time steps, create a yearly summary
    # Select one frame per year (e.g., June of each year)
    yearly_indices <- which(month(dates) == 6)  # June frames
    
    if (length(yearly_indices) > 1) {
      yearly_frames <- frames[yearly_indices]
      yearly_gif <- image_animate(image_join(yearly_frames), fps = 1, optimize = TRUE)
      yearly_file <- paste0(output_dir, "SPI_", scale, "_yearly_animation.gif")
      image_write(yearly_gif, yearly_file)
      cat("  Yearly summary GIF saved as:", yearly_file, "\n")
    }
  }
}

cat("\n==========================================\n")
cat("SPI GIF creation complete!\n")
cat("All animations saved in:", output_dir, "\n")

# ===================== OPTIONAL: CREATE COMPOSITE GIF =====================
# Create a side-by-side comparison of different SPI scales
cat("\nCreating composite comparison GIF...\n")

# Select a subset of scales for composite
composite_scales <- c(3, 6, 12)  # Ensure these files exist

composite_frames <- list()

# Get common time steps (minimum number of layers across all scales)
layer_counts <- sapply(composite_scales, function(s) {
  file <- paste0(spi_data_dir, "spi_", s, "month.nc")
  if (file.exists(file)) nlyr(rast(file)) else 0
})
min_layers <- min(layer_counts[layer_counts > 0])

if (min_layers > 0) {
  for (i in 1:min_layers) {
    cat("\r  Processing composite frame", i, "of", min_layers)
    
    # Create a composite plot with multiple panels
    plot_list <- list()
    
    for (scale in composite_scales) {
      spi_file <- paste0(spi_data_dir, "spi_", scale, "month.nc")
      if (!file.exists(spi_file)) next
      
      spi_data <- rast(spi_file)
      frame_date <- time(spi_data[[i]])
      
      # Create individual plot for this scale
      p <- create_spi_frame(spi_data, i, basin_shp, scale) +
        labs(subtitle = paste0("SPI-", scale, ": ", format(frame_date, "%b %Y")))
      
      plot_list[[as.character(scale)]] <- p
    }
    
    # Combine plots (requires patchwork package)
    if (length(plot_list) > 1) {
      # Install patchwork if not available
      if (!require("patchwork", quietly = TRUE)) {
        install.packages("patchwork")
        library(patchwork)
      }
      
      composite_plot <- wrap_plots(plot_list, ncol = 2) +
        plot_annotation(
          title = "SPI Drought Index Comparison",
          theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
        )
      
      # Save temporary file
      temp_file <- tempfile(fileext = ".png")
      ggsave(temp_file, composite_plot, width = 16, height = 10, dpi = 100)
      composite_frames[[i]] <- image_read(temp_file)
      file.remove(temp_file)
    }
  }
  
  # Create composite GIF
  if (length(composite_frames) > 0) {
    composite_gif <- image_animate(image_join(composite_frames), fps = 1, optimize = TRUE)
    composite_file <- paste0(output_dir, "SPI_comparison_animation.gif")
    image_write(composite_gif, composite_file)
    cat("\n  Composite GIF saved as:", composite_file, "\n")
  }
}

cat("\nAll SPI animations completed successfully!\n")