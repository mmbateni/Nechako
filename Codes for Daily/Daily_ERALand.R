library(ecmwfr)
library(terra)
library(foreach)
library(doParallel)

# --- Configuration -----------------------------------------------------------
user_email <- "mehdi.bateni@iusspavia.it"
pat_token  <- "44ace8e1-8ec8-4246-a44e-f30409c90c5b"
wf_set_key(user = user_email, key = pat_token)

setwd("D:/Nechako_Drought/")
dir.create("monthly_data", showWarnings = FALSE)
dir.create("temp_monthly", showWarnings = FALSE)

# Area: North, West, South, East
bbox  <- c(56.2, -127.8, 52.9, -122.7)
years <- 1950:2024

# --- Fixed Variable Configuration --------------------------------------------
var_cfg <- list(
  "total_precipitation" = list(
    dataset      = "derived-era5-land-daily-statistics",
    daily_stat   = "daily_sum",
    monthly_fun  = "sum",
    convert      = function(r) r * 1000,  # Convert m to mm
    needs_ptype  = FALSE
  ),
  "2m_temperature" = list(
    dataset      = "derived-era5-land-daily-statistics",
    daily_stat   = "daily_mean",
    monthly_fun  = "mean",
    convert      = function(r) r - 273.15,  # Convert K to °C
    needs_ptype  = FALSE
  ),
  "potential_evaporation" = list(
    dataset      = "derived-era5-land-daily-statistics",
    daily_stat   = "daily_sum",
    monthly_fun  = "sum",
    convert      = function(r) r * 1000,  # Convert m to mm
    needs_ptype  = FALSE
  )
)

# --- Helper function to get days in month (handles leap years) -----------------
get_days_in_month <- function(year, month) {
  days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  # Handle leap year for February
  if (month == 2 && ((year %% 4 == 0 && year %% 100 != 0) || year %% 400 == 0)) {
    return(29)
  }
  return(days_in_month[month])
}

# --- Parallel Setup ----------------------------------------------------------
num_cores <- 3
cl <- makeCluster(num_cores)
registerDoParallel(cl)

for (var_name in names(var_cfg)) {
  cfg <- var_cfg[[var_name]]
  cat("\n--- Processing:", var_name, "---\n")
  
  foreach(yr = years, .packages = c("ecmwfr", "terra")) %dopar% {
    # Register key for each worker
    ecmwfr::wf_set_key(user = user_email, key = pat_token)
    
    out_nc <- file.path("monthly_data", sprintf("%s_monthly_%d.nc", var_name, yr))
    if (file.exists(out_nc)) {
      cat(sprintf("Skipping %s %d - already exists\n", var_name, yr))
      return(NULL)
    }
    
    monthly_rasters <- list()
    
    for (m in 1:12) {
      month_str <- sprintf("%02d", m)
      days_in_month <- get_days_in_month(yr, m)
      days_str <- sprintf("%02d", 1:days_in_month)
      
      cat(sprintf("Processing %s %d-%02d (%d days)\n", var_name, yr, m, days_in_month))
      
      # Request structure matching Python sample
      req <- list(
        dataset_short_name = cfg$dataset,
        variable           = var_name,
        year               = as.character(yr),
        month              = month_str,
        day                = days_str,
        daily_statistic    = cfg$daily_stat,
        time_zone          = "utc+00:00",
        frequency          = "1_hourly"
      )
      
      tryCatch({
        # Download monthly data
        temp_file <- file.path("temp_monthly", sprintf("%s_%d_%02d.nc", var_name, yr, m))
        
        if (!file.exists(temp_file)) {
          result <- wf_request(
            request  = req,
            user     = user_email,
            transfer = TRUE,
            path     = temp_file,
            time_out = 3600
          )
          cat(sprintf("Downloaded %s %d-%02d\n", var_name, yr, m))
        } else {
          result <- temp_file
          cat(sprintf("Using cached file for %s %d-%02d\n", var_name, yr, m))
        }
        
        if (file.exists(result) && file.size(result) > 0) {
          # Process the data
          r <- rast(result)
          
          # Convert values (mm for precip/PET, °C for temp)
          r_converted <- cfg$convert(r)
          
          # Calculate monthly statistic
          if (cfg$monthly_fun == "sum") {
            monthly_result <- sum(r_converted, na.rm = TRUE)
          } else if (cfg$monthly_fun == "mean") {
            monthly_result <- mean(r_converted, na.rm = TRUE)
          }
          
          monthly_rasters[[m]] <- monthly_result
          
          # Cleanup temporary file
          #file.remove(result)
        } else {
          warning(sprintf("File %s is empty or doesn't exist", result))
        }
        
      }, error = function(e) {
        error_msg <- paste("ERROR:", yr, month_str, var_name, conditionMessage(e))
        cat(error_msg, "\n")
        write(paste(Sys.time(), "|", error_msg, sep = " "),
              file = "error_log.txt", append = TRUE)
      })
    }
    
    # Combine monthly results into yearly dataset
    if (length(monthly_rasters) == 12) {
      yearly_raster <- rast(monthly_rasters)
      
      # Set proper time dimension
      dates <- seq(as.Date(paste0(yr, "-01-01")), as.Date(paste0(yr, "-12-01")), by = "month")
      time(yearly_raster) <- dates
      
      # Write to NetCDF
      writeCDF(yearly_raster, out_nc, overwrite = TRUE)
      cat(sprintf("Successfully saved %s for year %d\n", var_name, yr))
    } else {
      warning(sprintf("Incomplete data for %s %d - only %d months processed", 
                      var_name, yr, length(monthly_rasters)))
    }
    
    return(NULL)
  }
}

stopCluster(cl)
cat("\n--- Process Finished ---\n")