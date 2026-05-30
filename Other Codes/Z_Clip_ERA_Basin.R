# Load required packages
library(terra)
library(sf)
library(httr)
library(curl)

# 1. Download the Nechako river basin shapefile from GitHub
# First, create a temporary directory to store the shapefile components
temp_dir <- tempdir()
setwd(temp_dir)

# The GitHub repository likely contains multiple shapefile components (.shp, .shx, .dbf, .prj)
# We need to download all of them. Common naming patterns in this repository might be:
shapefile_name <- "Nechako_basin"  # This might need to be adjusted based on actual filename

# Download shapefile components - you may need to adjust these URLs based on actual files
download.file("https://raw.githubusercontent.com/mmbateni/Nechako/master/Spatial/Nechako_basin.shp", 
              destfile = "Nechako_basin.shp", mode = "wb")
download.file("https://raw.githubusercontent.com/mmbateni/Nechako/master/Spatial/Nechako_basin.shx", 
              destfile = "Nechako_basin.shx", mode = "wb")
download.file("https://raw.githubusercontent.com/mmbateni/Nechako/master/Spatial/Nechako_basin.dbf", 
              destfile = "Nechako_basin.dbf", mode = "wb")
download.file("https://raw.githubusercontent.com/mmbateni/Nechako/master/Spatial/Nechako_basin.prj", 
              destfile = "Nechako_basin.prj", mode = "wb")

# 2. Read the shapefile
if (file.exists("Nechako_basin.shp")) {
  nechako_basin <- st_read("Nechako_basin.shp")
  message("Successfully loaded Nechako river basin shapefile")
} else {
  # Alternative approach - try other common filenames
  possible_names <- c("Nechako", "basin", "Nechako_river", "Nechako_watershed", "NechakoBasin")
  for (name in possible_names) {
    if (file.exists(paste0(name, ".shp"))) {
      nechako_basin <- st_read(paste0(name, ".shp"))
      message(paste("Successfully loaded shapefile with name:", name))
      break
    }
  }
}

# 3. Get the bounding box of the Nechako basin
# Convert to the same CRS as your ERA5 data (likely WGS84)
nechako_basin <- st_transform(nechako_basin, 4326)
nechako_bbox <- st_bbox(nechako_basin)

# Create a SpatExtent object for terra
nechako_extent <- ext(nechako_bbox["xmin"], nechako_bbox["xmax"], 
                      nechako_bbox["ymin"], nechako_bbox["ymax"])

message("Nechako basin bounding box:")
print(nechako_extent)

# 4. Read and crop the ERA5-Land data to the Nechako basin extent
# Set your file path here
era5_file <- "era5_land_snow_depth_merged.nc"

# Read the NetCDF file
era5_rast <- rast(era5_file)

# Crop to the Nechako basin extent
message("Cropping ERA5 data to Nechako river basin...")
era5_nechako <- crop(era5_rast, nechako_extent)

# 5. Optional: Mask the data to the exact basin boundary (not just bounding box)
message("Masking data to exact basin boundary...")
era5_nechako_masked <- mask(era5_nechako, vect(nechako_basin))

# 6. Save the cropped data (optional)
output_file <- "era5_nechako_basin_snow_depth.nc"
message(paste("Saving cropped data to:", output_file))
writeCDF(era5_nechako_masked, output_file, overwrite = TRUE)

# 7. For visualization - extract a specific time slice (e.g., 2020-01-01)
target_time <- "2020-01-01"
matching_layers <- grep(target_time, names(era5_nechako_masked))
if (length(matching_layers) > 0) {
  jan_2020_data <- era5_nechako_masked[[matching_layers[1]]]
  
  # Plot the data
  plot(jan_2020_data, 
       main = paste("Snow Depth - Nechako Basin", target_time),
       col = hcl.colors(100, "viridis"))
} else {
  message("No data found for the specified date. Showing first time slice instead.")
  plot(era5_nechako_masked[[1]], 
       main = "Snow Depth - Nechako Basin (First Time Slice)",
       col = hcl.colors(100, "viridis"))
}

# 8. Alternative approach if GitHub download fails - manual download
# If the GitHub download fails, you can manually download the shapefile and use:
# nechako_basin <- st_read("path/to/your/downloaded/shapefile.shp")