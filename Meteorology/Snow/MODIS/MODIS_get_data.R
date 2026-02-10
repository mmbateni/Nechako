# -------------------------------
# 1. Load Required Libraries
# -------------------------------
install.packages(c("terra", "sf", "lubridate", "dplyr", "zoo", "ggplot2", "httr", "xml2", "jsonlite"))
library(terra)
library(sf)
library(lubridate)
library(dplyr)
library(zoo)
library(ggplot2)
library(httr)
library(xml2)
library(jsonlite)

# -------------------------------
# 2. Setup Directory Structure
# -------------------------------
base_dir <- "D:/Nechako_Drought"
if (!dir.exists(file.path(base_dir, "MODIS_raw"))) dir.create(file.path(base_dir, "MODIS_raw"), recursive = TRUE)
if (!dir.exists(file.path(base_dir, "MODIS_processed"))) dir.create(file.path(base_dir, "MODIS_processed"), recursive = TRUE)

# -------------------------------
# 3. Manual Download Instructions (Most Reliable Method)
# -------------------------------
cat("=== MANUAL DOWNLOAD REQUIRED ===\n")
cat("1. Go to NASA Earthdata Search: https://search.earthdata.nasa.gov/\n")
cat("2. Search for: 'MOD10A1.061' (Terra Snow Cover Daily L3 Global 500m)\n")
cat("3. Set date range: January 1, 2020 to December 31, 2020\n")
cat("4. Draw bounding box for Nechako River Basin:\n")
cat("   - West: -128.5, East: -122\n")
cat("   - South: 51.9, North: 56.5\n")
cat("5. Download HDF files to: ", file.path(base_dir, "MODIS_raw"), "\n")
cat("6. After downloading, run the processing code below\n")

# -------------------------------
# 4. Alternative: Semi-Automated Download (If You Have Earthdata Account)
# -------------------------------
download_modis_data <- function() {
  cat("NOTE: This requires an Earthdata account (free at https://urs.earthdata.nasa.gov)\n")
  username <- readline("Enter your Earthdata username: ")
  password <- readline("Enter your Earthdata password: ", ask = TRUE)
  
  # MOD10A1.061 collection ID
  collection_id <- "C2045979734-LPDAAC_ECS"
  
  # Date range
  start_date <- "2020-01-01"
  end_date <- "2020-12-31"
  
  # Bounding box for Nechako Basin
  bbox <- "-128.5,51.9,-122,56.5"
  
  # Build CMR query URL
  cmr_url <- sprintf(
    "https://cmr.earthdata.nasa.gov/search/granules.json?collection_id=%s&temporal=%s,%s&bounding_box=%s&page_size=2000",
    collection_id, start_date, end_date, bbox
  )
  
  cat("Searching for MODIS files...\n")
  response <- GET(cmr_url, authenticate(username, password))
  
  if (status_code(response) != 200) {
    cat("Error: Could not retrieve file list. Status code:", status_code(response), "\n")
    return(NULL)
  }
  
  # Parse JSON response
  json_data <- content(response, "parsed")
  granules <- json_data$feed$entry
  
  cat(sprintf("Found %d MODIS files to download\n", length(granules)))
  
  # Download files
  download_results <- list()
  for (i in 1:length(granules)) {
    granule <- granules[[i]]
    url <- granule$link[[1]]$href
    filename <- basename(url)
    output_path <- file.path(base_dir, "MODIS_raw", filename)
    
    cat(sprintf("Downloading [%d/%d]: %s\n", i, length(granules), filename))
    
    # Download with progress
    download_result <- tryCatch({
      response <- GET(url, authenticate(username, password), 
                      write_disk(output_path, overwrite = TRUE),
                      progress())
      if (status_code(response) == 200) {
        "Success"
      } else {
        sprintf("Failed (status %d)", status_code(response))
      }
    }, error = function(e) {
      sprintf("Error: %s", e$message)
    })
    
    download_results[[filename]] <- download_result
  }
  
  return(download_results)
}

# Uncomment to run download (requires Earthdata account)
# download_results <- download_modis_data()

# -------------------------------
# 5. Process Downloaded HDF Files (Run After Download)
# -------------------------------
process_modis_files <- function() {
  # List all HDF files in raw directory
  hdf_files <- list.files(file.path(base_dir, "MODIS_raw"), 
                          pattern = "MOD10A1.*\\.hdf$", 
                          full.names = TRUE, recursive = TRUE)
  
  if (length(hdf_files) == 0) {
    cat("No HDF files found! Please download data first.\n")
    cat("Download from: https://search.earthdata.nasa.gov/\n")
    return(NULL)
  }
  
  cat(sprintf("Found %d HDF files to process\n", length(hdf_files)))
  
  # Study area extent
  nechako_ext <- ext(c(xmin = -128.5, xmax = -122, ymin = 51.9, ymax = 56.5))
  
  # Container for results
  snow_data <- data.frame()
  
  # Process each file
  for (i in 1:length(hdf_files)) {
    hdf_file <- hdf_files[i]
    filename <- basename(hdf_file)
    
    cat(sprintf("Processing [%d/%d]: %s\n", i, length(hdf_files), filename))
    
    tryCatch({
      # Extract date from filename (format: MOD10A1.A2020001.h25v04.061.2020003123456.hdf)
      date_match <- regmatches(filename, regexpr("A[0-9]{7}", filename))
      if (length(date_match) > 0) {
        julian_date <- substr(date_match, 2, 8)  # Remove 'A' prefix
        year <- as.numeric(substr(julian_date, 1, 4))
        day_of_year <- as.numeric(substr(julian_date, 5, 7))
        date <- as.Date(sprintf("%d-%03d", year, day_of_year), format = "%Y-%j")
      } else {
        cat(sprintf("Warning: Could not extract date from filename: %s\n", filename))
        next
      }
      
      # Read NDSI_Snow_Cover band (typically band 1 in MOD10A1)
      r <- tryCatch({
        rast(hdf_file, subds = "NDSI_Snow_Cover")
      }, error = function(e) {
        cat(sprintf("Warning: Could not read NDSI_Snow_Cover from %s: %s\n", filename, e$message))
        # Try alternative subdataset names
        try(rast(hdf_file, subds = "1"), silent = TRUE)
      })
      
      if (is.null(r) || class(r)[1] != "SpatRaster") {
        cat(sprintf("Warning: Could not read raster data from %s\n", filename))
        next
      }
      
      # Project to geographic coordinates if needed
      if (crs(r) != "EPSG:4326") {
        r <- project(r, "EPSG:4326")
      }
      
      # Crop to Nechako basin
      r_cropped <- crop(r, nechako_ext)
      
      # Calculate mean snow cover (values 0-100, where 200+ are special values)
      # MOD10A1 uses: 0-100 = snow cover %, 200 = cloud, 201 = no decision, 210+ = other
      snow_values <- values(r_cropped)
      valid_snow <- snow_values[snow_values >= 0 & snow_values <= 100]
      
      if (length(valid_snow) > 0) {
        mean_snow <- mean(valid_snow, na.rm = TRUE)
        
        # Save processed raster
        output_file <- file.path(base_dir, "MODIS_processed", 
                                 sprintf("snow_cover_%s.tif", format(date, "%Y%m%d")))
        writeRaster(r_cropped, output_file, overwrite = TRUE)
        
        # Add to results
        snow_data <- rbind(snow_data, data.frame(
          date = date,
          mean_snow_cover = mean_snow,
          file = output_file
        ))
        
        cat(sprintf("  Date: %s, Mean Snow Cover: %.2f%%\n", format(date, "%Y-%m-%d"), mean_snow))
      } else {
        cat(sprintf("  Warning: No valid snow cover values found for %s\n", format(date, "%Y-%m-%d")))
      }
      
    }, error = function(e) {
      cat(sprintf("  Error processing %s: %s\n", filename, e$message))
    })
    
    # Small delay to avoid overwhelming system
    Sys.sleep(0.1)
  }
  
  if (nrow(snow_data) == 0) {
    cat("No data processed successfully!\n")
    return(NULL)
  }
  
  # Sort by date
  snow_data <- snow_data[order(snow_data$date), ]
  
  # Save daily data
  write.csv(snow_data, file.path(base_dir, "MODIS_processed", "daily_snow_cover_2020.csv"), row.names = FALSE)
  
  return(snow_data)
}

# Run processing
snow_data <- process_modis_files()

# -------------------------------
# 6. Apply 7-Day Temporal Filter & Monthly Aggregation
# -------------------------------
if (!is.null(snow_data) && nrow(snow_data) > 7) {
  # Apply 7-day moving average (cloud-aware)
  snow_data$filtered_snow <- NA
  
  for (i in 4:(nrow(snow_data) - 3)) {
    window_vals <- snow_data$mean_snow_cover[(i-3):(i+3)]
    if (sum(!is.na(window_vals)) > 0) {
      snow_data$filtered_snow[i] <- mean(window_vals, na.rm = TRUE)
    }
  }
  
  # Remove NA values
  filtered_data <- snow_data[!is.na(snow_data$filtered_snow), ]
  
  # Monthly aggregation
  filtered_data$month <- format(filtered_data$date, "%Y-%m")
  monthly_data <- filtered_data %>%
    group_by(month) %>%
    summarise(
      mean_monthly_snow = mean(filtered_snow, na.rm = TRUE),
      mid_month_date = as.Date(paste0(month, "-15")),
      n_days = n()
    ) %>%
    ungroup()
  
  # -------------------------------
  # 7. Plot Results
  # -------------------------------
  ggplot(monthly_data, aes(x = mid_month_date, y = mean_monthly_snow)) +
    geom_line(color = "blue", linewidth = 1.5, na.rm = TRUE) +
    geom_point(color = "red", size = 2, na.rm = TRUE) +
    labs(
      title = "Monthly Mean Snow Cover in Nechako River Basin (2020)",
      subtitle = "MODIS MOD10A1.061 | 7-day temporal filter applied",
      x = "Date", 
      y = "Mean Snow Cover (%)",
      caption = "Data source: NASA MODIS via Earthdata Search"
    ) +
    scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "gray50"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Save plot
  ggsave(file.path(base_dir, "MODIS_processed", "snow_cover_trend_2020.png"), 
         width = 12, height = 8, dpi = 300)
  
  # Save monthly data
  write.csv(monthly_data, file.path(base_dir, "MODIS_processed", "monthly_snow_cover_2020.csv"), row.names = FALSE)
  
  cat("Analysis complete! Results saved to:", file.path(base_dir, "MODIS_processed"), "\n")
  
} else {
  cat("Not enough data for analysis. Please ensure you have downloaded MODIS files first.\n")
}