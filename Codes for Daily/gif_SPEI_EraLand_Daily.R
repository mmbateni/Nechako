# ============================================================================
# SCRIPT 5: Create Animated GIFs for SPEI Indices
# ============================================================================
# This script creates GIF animations for SPEI timeseries with basin overlay

# Load required libraries
library(terra)      # For raster processing
library(magick)     # For GIF creation
library(ggplot2)    # For plotting
library(viridis)    # For color palettes

# ===================== CONFIGURATION =====================
# Define SPEI timescales to visualize
spei_scales <- c(3, 6, 12)  # You can add 1, 9, 24 if needed

# Path to SPEI data (from Script 3 output)
spei_data_dir <- "spei_results/"

# Path to basin shapefile
basin_shp_path <- "Spatial/nechakoBound_dissolve.shp"

# Output directory for GIFs
output_dir <- "spei_gifs/"
dir.create(output_dir, showWarnings = FALSE)

# ===================== FUNCTION: CREATE SINGLE SPEI FRAME =====================
create_spei_frame <- function(spei_raster, time_index, basin_shp, scale) {
  # Creates a single frame for the SPEI GIF animation
  
  # Extract the specific time layer
  spei_layer <- spei_raster[[time_index]]
  layer_date <- time(spei_layer)
  
  # Convert to data frame for ggplot
  spei_df <- as.data.frame(spei_layer, xy = TRUE)
  colnames(spei_df) <- c("x", "y", "spei_value")
  
  # Convert basin shape to data frame for ggplot
  basin_df <- as.data.frame(basin_shp, geom = "XY")
  
  # Create the plot
  p <- ggplot() +
    # Plot SPEI raster
    geom_raster(data = spei_df, aes(x = x, y = y, fill = spei_value)) +
    
    # Add basin boundaries
    geom_polygon(data = basin_df, aes(x = x, y = y, group = part),
                 fill = NA, color = "darkgreen", linewidth = 1.2) +
    
    # Color scale for SPEI (standard drought colors)
    scale_fill_gradientn(
      colors = c("#8B0000", "#FF0000", "#FFA500", "#FFFF00", "#ADFF2F", "#00FA9A", "#1E90FF"),
      limits = c(-3, 3),
      breaks = seq(-3, 3, 1),
      name = "SPEI Value",
      na.value = "grey90"
    ) +
    
    # Title and labels
    labs(
      title = paste0("SPEI-", scale, " Drought Index"),
      subtitle = paste(format(layer_date, "%B %Y"), 
                       "\nNechako River Basin"),
      x = "Longitude",
      y = "Latitude",
      caption = "Negative values = Drought conditions"
    ) +
    
    # Theme customization
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
    
    # Coordinate system
    coord_sf(crs = crs(spei_layer), default = TRUE) +
    
    # Legend position
    theme(legend.position = "right")
  
  return(p)
}

# ===================== MAIN PROCESSING =====================
# Load basin shapefile
basin_shp <- vect(basin_shp_path)

# Process each SPEI scale
for (scale in spei_scales) {
  cat("Creating GIF for SPEI-", scale, "...\n", sep = "")
  
  # Load SPEI data for this scale
  spei_file <- paste0(spei_data_dir, "spei_", scale, "month.nc")
  
  if (!file.exists(spei_file)) {
    warning("File not found: ", spei_file, ". Skipping SPEI-", scale)
    next
  }
  
  spei_data <- rast(spei_file)
  dates <- time(spei_data)
  
  # Create a list to store individual frames
  frames <- list()
  
  # Generate each frame
  for (i in 1:nlyr(spei_data)) {
    cat("\r  Processing frame", i, "of", nlyr(spei_data))
    
    # Create the frame
    frame_plot <- create_spei_frame(spei_data, i, basin_shp, scale)
    
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
    fps = 2,  # Frames per second
    optimize = TRUE
  )
  
  # Save the GIF
  gif_file <- paste0(output_dir, "SPEI_", scale, "_animation.gif")
  image_write(gif, gif_file)
  
  cat("  GIF saved as:", gif_file, "\n")
  
  # Create a progress bar version showing drought development
  cat("  Creating progress bar visualization...\n")
  
  # Extract SPEI values averaged over the basin for each time step
  basin_avg <- numeric(nlyr(spei_data))
  for (i in 1:nlyr(spei_data)) {
    # Crop and mask SPEI data to basin
    spei_cropped <- crop(spei_data[[i]], basin_shp)
    spei_masked <- mask(spei_cropped, basin_shp)
    
    # Calculate mean SPEI value within basin
    basin_avg[i] <- global(spei_masked, "mean", na.rm = TRUE)[1, 1]
  }
  
  # Create a time series plot of basin-average SPEI
  ts_data <- data.frame(
    Date = dates,
    SPEI = basin_avg
  )
  
  # Add drought severity categories
  ts_data$Category <- cut(ts_data$SPEI,
                          breaks = c(-Inf, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, Inf),
                          labels = c("Extreme Drought", "Severe Drought", "Moderate Drought", 
                                     "Mild Drought", "Normal", "Mild Wet", "Moderate Wet", 
                                     "Severe Wet", "Extreme Wet")
  )
  
  # Save the time series plot
  ts_plot_file <- paste0(output_dir, "SPEI_", scale, "_timeseries.png")
  png(ts_plot_file, width = 12, height = 6, units = "in", res = 150)
  
  plot(dates, basin_avg, type = "l", lwd = 2, col = "blue",
       main = paste0("Nechako River Basin: SPEI-", scale, " Time Series"),
       xlab = "Date", ylab = "SPEI Value",
       ylim = c(-3, 3))
  
  # Add reference lines
  abline(h = 0, lty = 2, col = "gray")
  abline(h = c(-1, -2, 1, 2), lty = 3, col = "gray")
  
  # Add drought/wetness regions
  rect(min(dates), -3, max(dates), -2, col = rgb(1, 0, 0, 0.2), border = NA)
  rect(min(dates), -2, max(dates), -1, col = rgb(1, 0.5, 0, 0.2), border = NA)
  rect(min(dates), 1, max(dates), 2, col = rgb(0, 0.5, 1, 0.2), border = NA)
  rect(min(dates), 2, max(dates), 3, col = rgb(0, 0, 1, 0.2), border = NA)
  
  # Redraw the line on top
  lines(dates, basin_avg, lwd = 2, col = "blue")
  
  # Add legend
  legend("bottomleft", 
         legend = c("Extreme Drought", "Severe Drought", "Moderate Drought", 
                    "Normal", "Moderate Wet", "Severe Wet"),
         fill = c(rgb(1, 0, 0, 0.5), rgb(1, 0.5, 0, 0.5), rgb(1, 1, 0, 0.5),
                  "white", rgb(0.5, 0.5, 1, 0.5), rgb(0, 0, 1, 0.5)),
         border = NA, bg = "white", cex = 0.8)
  
  dev.off()
  cat("  Time series plot saved as:", ts_plot_file, "\n")
}

cat("\n==========================================\n")
cat("SPEI GIF creation complete!\n")
cat("All animations saved in:", output_dir, "\n")

# ===================== CREATE DROUGHT EVENT SUMMARY =====================
cat("\nCreating drought event summary...\n")

# Identify drought events for the primary scale (e.g., SPEI-12)
primary_scale <- 12
spei_file <- paste0(spei_data_dir, "spei_", primary_scale, "month.nc")

if (file.exists(spei_file)) {
  spei_data <- rast(spei_file)
  dates <- time(spei_data)
  
  # Calculate basin-average SPEI for entire period
  basin_series <- numeric(nlyr(spei_data))
  for (i in 1:nlyr(spei_data)) {
    spei_cropped <- crop(spei_data[[i]], basin_shp)
    spei_masked <- mask(spei_cropped, basin_shp)
    basin_series[i] <- global(spei_masked, "mean", na.rm = TRUE)[1, 1]
  }
  
  # Identify drought events (SPEI < -0.5 for at least 3 months)
  drought_mask <- basin_series < -0.5
  
  # Find consecutive drought periods
  drought_events <- rle(drought_mask)
  
  # Summarize drought events
  drought_summary <- data.frame()
  start_idx <- 1
  
  for (j in 1:length(drought_events$lengths)) {
    if (drought_events$values[j] && drought_events$lengths[j] >= 3) {
      end_idx <- start_idx + drought_events$lengths[j] - 1
      
      event <- data.frame(
        Start_Date = dates[start_idx],
        End_Date = dates[end_idx],
        Duration_Months = drought_events$lengths[j],
        Min_SPEI = min(basin_series[start_idx:end_idx], na.rm = TRUE),
        Mean_SPEI = mean(basin_series[start_idx:end_idx], na.rm = TRUE)
      )
      
      drought_summary <- rbind(drought_summary, event)
    }
    start_idx <- start_idx + drought_events$lengths[j]
  }
  
  # Save drought summary
  if (nrow(drought_summary) > 0) {
    summary_file <- paste0(output_dir, "SPEI_", primary_scale, "_drought_events.csv")
    write.csv(drought_summary, summary_file, row.names = FALSE)
    cat("  Drought event summary saved as:", summary_file, "\n")
    
    # Print summary to console
    cat("\nDrought Events Detected (SPEI-", primary_scale, "):\n", sep = "")
    print(drought_summary)
  } else {
    cat("  No significant drought events detected.\n")
  }
}

cat("\nAll SPEI processing completed successfully!\n")