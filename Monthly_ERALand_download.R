# ============================================================================
# Download ERA5-Land Monthly Values
# ============================================================================
# Load required libraries
library(ecmwfr)  # For CDS API access
library(terra)   # For raster processing
library(lubridate)
library(ncdf4)

# Set working directory (ensure this path exists)
setwd("D:/Nechako_Drought/")

# ===================== CONFIGURATION =====================
# Set credentials 
wf_set_key(
  user =  "mehdi.bateni@iusspavia.it",
  key = "44ace8e1-8ec8-4246-a44e-f30409c90c5b"
)

# Study Area Bounding Box (British Columbia)
bbox <- c(-127.8, 52.9, -122.7, 56.2)  # xmin, ymin, xmax, ymax

# Time Period (adjust as needed)
start_date <- "1950-01-01"
end_date <- "2025-12-31"

# Variables needed for SPI and SPEI
variables <- c(
  "snow_cover",
  "snow_depth_water_equivalent"
  # "surface_solar_radiation_downwards",
  # "10m_u_component_of_wind",
  # "10m_v_component_of_wind",
  # "geopotential",
  # "surface_pressure",
  # "total_precipitation",
  # "potential_evaporation",
  # "2m_temperature",
  # "2m_dewpoint_temperature",
  # "skin_reservoir_content",
  # "volumetric_soil_water_layer_1",
  # "volumetric_soil_water_layer_2",
  # "volumetric_soil_water_layer_3",
  # "volumetric_soil_water_layer_4"
)

# Create output directories
dir.create("monthly_data_direct", showWarnings = FALSE)
dir.create("downloads", showWarnings = FALSE)

# ===================== FUNCTION: DOWNLOAD ALL DATA FOR ONE VARIABLE =====================
download_era5_variable <- function(variable, start_year, end_year) {
  # Downloads all monthly data for one variable across multiple years
  
  # Generate year and month vectors
  years <- as.character(start_year:end_year)
  months <- sprintf("%02d", 1:12)
  
  # Construct request for monthly means dataset
  request <- list(
    dataset_short_name = "reanalysis-era5-land-monthly-means",
    product_type = "monthly_averaged_reanalysis",
    variable = variable,
    year = years,  # All years at once
    month = months,  # All months at once
    time = "00:00",
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(bbox[4], bbox[1], bbox[2], bbox[3])  # North, West, South, East
  )
  
  # Submit request to CDS
  file <- wf_request(
    user = "mehdi.bateni@iusspavia.it",
    request = request,
    transfer = TRUE,
    path = "downloads/",
    time_out = 7200  # Increased timeout for large downloads (2 hours)
  )
  
  return(file)
}

# ===================== MAIN PROCESSING =====================
cat("Starting download of ERA5-Land monthly means data...\n")
cat("Time period:", start_date, "to", end_date, "\n\n")

start_year <- year(as.Date(start_date))
end_year <- year(as.Date(end_date))

# Wrap entire processing in tryCatch to ensure safe cleanup only on success
tryCatch({
  # Process each variable
  for (var in variables) {
    cat("Processing variable:", var, "\n")
    cat("  Downloading data for", start_year, "-", end_year, "... ")
    
    # Download all data for this variable
    nc_file <- download_era5_variable(var, start_year, end_year)
    
    cat("Downloaded.\n")
    cat("  Reading NetCDF file... ")
    
    # Read the NetCDF file as SpatRaster
    result <- rast(nc_file)
    
    cat("Done.\n")
    cat("  Processing", nlyr(result), "layers... ")
    
    # --------------------------------------------------------------------------
    # NEW BLOCK: Check and Report NA Values
    # --------------------------------------------------------------------------
    cat("\n  --- NA Check ---\n")
    
    # Count NAs across the entire raster stack
    total_cells <- ncell(result) * nlyr(result)
    na_count <- global(result, fun = "isNA", na.rm = FALSE) 
    total_na <- sum(na_count$isNA)
    
    # Calculate percentage
    na_percent <- (total_na / total_cells) * 100
    
    # Report to console
    cat(sprintf("  Total Cells: %d\n", total_cells))
    cat(sprintf("  Total NAs:   %d\n", total_na))
    cat(sprintf("  NA Percent:  %.4f%%\n", na_percent))
    
    if (total_na > 0) {
      cat("  WARNING: Missing values detected.\n")
    } else {
      cat("  STATUS: No missing values.\n")
    }
    cat("  ----------------\n")
    # --------------------------------------------------------------------------
    
    # Apply unit conversions if needed
    if (var == "total_precipitation") {
      # Convert from meters to millimeters (1 m = 1000 mm)
      result <- result * 1000
      cat("  Converted precipitation from m to mm.\n")
    }
    
    if (grepl("temperature", var)) {
      # Convert from Kelvin to Celsius if desired
      # result <- result - 273.15
      # cat("  Converted temperature from K to °C.\n")
    }
    
    # Save to file
    output_file <- paste0("monthly_data_direct/", var, "_monthly.nc")
    cat("  Saving to:", output_file, "... ")
    writeCDF(result, output_file, overwrite = TRUE)
    cat("Done.\n")
    
    # Clean up downloaded file
    file.remove(nc_file)
    
    cat("  Completed", var, "\n\n")
  }
  
  cat("All downloads complete!\n")
  
  # ===================== CLEANUP: DELETE DOWNLOADS FOLDER =====================
  cat("\n--- Cleanup Phase ---\n")
  downloads_path <- "downloads"
  
  if (dir.exists(downloads_path)) {
    # Verify directory is empty (all files should have been removed during processing)
    if (length(list.files(downloads_path)) == 0) {
      unlink(downloads_path, recursive = TRUE)
      cat(sprintf("✓ Successfully deleted folder: %s\n", file.path(getwd(), downloads_path)))
    } else {
      warning(sprintf("Downloads folder not empty (%d files remaining). Skipping deletion for safety.",
                      length(list.files(downloads_path))))
    }
  } else {
    cat(sprintf("✓ Downloads folder already removed or never created.\n"))
  }
  
  cat("\nScript execution completed successfully.\n")
  
}, error = function(e) {
  cat("\n❌ ERROR: Processing failed!\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("\n⚠️  Downloads folder NOT deleted for debugging purposes.\n")
  cat("You can manually inspect 'D:/Nechako_Drought/downloads/' before deletion.\n")
  stop(e)
})