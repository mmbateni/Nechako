# ---- 1. Setup & Libraries ----
library(ecmwfr)   # Pure R interface to CDS (No Python setup needed)
library(terra)    # Memory-safe raster handling
library(lubridate)# Robust date handling (Handles leap years correctly)

# User Credentials (REPLACE WITH YOUR DETAILS)
user_uid <- "mehdi.bateni@iusspavia.it"
user_key <- "f0b6de88-c627-4f98-98b4-22453fb77cac"

# Register the key
wf_set_key(user = user_uid, key = user_key)

# Settings
dataset <- "reanalysis-era5-land"
variable_name <- "snowfall" 

# Directory Setup
base_download_directory <- "D:/Nechako_Drought/ERA5-Land download-N and B"
output_directory <- file.path(base_download_directory, variable_name)
if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

# Time & Area Settings
start_year <- 1950
end_year <- 2024
months_of_interest <- 1:12 # Full year (Python code used 1:10, change if needed)
area_bbox <- c(56.5, -128.5, 51.9, -122.0) # [North, West, South, East]
times <- sprintf("%02d:00", 0:23) # "00:00" to "23:00"

# ---- 2. Download Loop (Improved) ----
for (year in start_year:end_year) {
  for (month_num in months_of_interest) {
    
    # IMPROVEMENT: Use lubridate for exact days in month (handles leap years)
    current_date <- make_date(year, month_num, 1)
    days_in_current_month <- days_in_month(current_date)
    days_list <- sprintf("%02d", 1:days_in_current_month)
    month_str <- sprintf("%02d", month_num)
    
    # Define filename
    file_name <- sprintf("era5_%s_%d_%s.nc", variable_name, year, month_str)
    full_path <- file.path(output_directory, file_name)
    
    if (file.exists(full_path)) {
      message(paste("Skipping:", file_name, "(Already exists)"))
      next
    }
    
    message(paste("Requesting:", year, "-", month_str))
    
    request <- list(
      product_type = "reanalysis",
      format = "netcdf",
      variable = variable_name,
      year = as.character(year),
      month = month_str,
      day = days_list,
      time = times,
      area = area_bbox,
      dataset_short_name = dataset,
      target = file_name
    )
    
    # Robust error handling
    tryCatch({
      wf_request(request = request, user = user_uid, path = output_directory, verbose = FALSE)
      message(paste("Saved:", file_name))
    }, error = function(e) {
      message(paste("FAILED:", year, "-", month_str, "Error:", e$message))
    })
  }
}

message("Downloads complete. Starting Merge...")

# ---- 3. Smart Merge (New Hybrid Approach) ----
merged_filename <- file.path(output_directory, paste0("era5_land_", variable_name, "_merged.nc"))
nc_files <- list.files(output_directory, pattern = paste0("era5_", variable_name, "_.*\\.nc$"), full.names = TRUE)

if (length(nc_files) == 0) stop("No files found to merge.")

# Check if CDO is available (The "Fast" option from the other script)
cdo_path <- Sys.which("cdo")

if (nzchar(cdo_path)) {
  # OPTION A: CDO (Fastest, best for Linux/Power Users)
  message("CDO found! Using high-speed merge...")
  
  # On Windows, command line length is limited. We write file list to a text file.
  file_list_txt <- file.path(output_directory, "filelist.txt")
  writeLines(nc_files, file_list_txt)
  
  # Execute CDO merge
  # Syntax: cdo mergetime [files] output
  # Using system2 to call the external tool
  tryCatch({
    # Note: 'cat' is used to pass the list, varies by OS. 
    # For Windows simplicity, we just pass the glob pattern if possible, 
    # or loop. Here we try a direct wildcard merge which CDO supports.
    wildcard_path <- file.path(output_directory, paste0("era5_", variable_name, "_*.nc"))
    
    system2("cdo", args = c("mergetime", wildcard_path, merged_filename))
    message(paste("Merge Successful via CDO:", merged_filename))
    
  }, error = function(e) {
    message("CDO merge failed, falling back to Terra...")
  })
  
} else {
  # OPTION B: Terra (Memory-Safe R fallback)
  # This is superior to 'abind' because it does NOT load 74 years into RAM.
  message("CDO not found. Using Terra (Memory-safe R merge)...")
  
  # Create a virtual raster stack (lazy loading)
  r_stack <- terra::rast(nc_files)
  
  # Sort by time just in case file system order is wrong
  # (Terra usually respects the time index inside the NetCDF, but sorting files helps)
  r_stack <- r_stack[[order(time(r_stack))]]
  
  message(paste("Merging", length(nc_files), "files... this may take a few minutes."))
  
  terra::writeCDF(r_stack, merged_filename, 
                  overwrite = TRUE, 
                  varname = variable_name,
                  compression = 5) # Added compression to save disk space
  
  message(paste("Merge Successful via Terra:", merged_filename))
}