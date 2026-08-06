# ============================================================================
# Unified ERA5-Land + ERA5 Atmospheric Circulation/SST Monthly Download
#
# Merged from Monthly_ERALand_download.R and Monthly_ERA_atm_download.R.
# ============================================================================

# ── Libraries ────────────────────────────────────────────────────────────────
library(ecmwfr)
library(terra)
library(lubridate)
library(ncdf4)
library(dotenv)

# ── Environment and working directory ───────────────────────────────────────
dotenv::load_dot_env(".env")

WD_PATH <- Sys.getenv("WD_PATH")
CDS_USER <- Sys.getenv("CDS_USER")
CDS_KEY <- Sys.getenv("CDS_KEY")

if (!nzchar(WD_PATH)) stop("WD_PATH is not defined in .env or the environment.")
if (!nzchar(CDS_USER)) stop("CDS_USER is not defined in .env or the environment.")
if (!nzchar(CDS_KEY)) stop("CDS_KEY is not defined in .env or the environment.")
if (!dir.exists(WD_PATH)) stop("WD_PATH does not exist: ", WD_PATH)

setwd(WD_PATH)

# ── CDS credentials ──────────────────────────────────────────────────────────
wf_set_key(user = CDS_USER, key = CDS_KEY)

# ── Shared output directories ────────────────────────────────────────────────
dir.create("monthly_data_direct", showWarnings = FALSE, recursive = TRUE)
dir.create("downloads", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# PART 1 — ERA5-LAND MONTHLY VALUES
# ============================================================================

# Study Area Bounding Box (British Columbia)
bbox <- c(-127.8, 52.9, -122.7, 56.2)  # xmin, ymin, xmax, ymax

# Time Period (adjust as needed)
start_date <- "1950-01-01"
end_date <- "2025-12-31"

# Variables needed for SPI and SPEI
variables <- c(
  #"snow_cover",
  #"snow_depth_water_equivalent",      # SWE proxy (m water equiv.) — Basin Memory
  # "surface_solar_radiation_downwards",
  # "10m_u_component_of_wind",
  # "10m_v_component_of_wind",
  # "geopotential",
  # "surface_pressure",
  "total_precipitation"
  # "potential_evaporation"
  #"2m_temperature",
  #"2m_dewpoint_temperature",
  # "skin_reservoir_content",
  #"volumetric_soil_water_layer_1",    # 0–7 cm   — near-surface, responds to AR events
  #"volumetric_soil_water_layer_2"     # 7–28 cm  — shallow memory, most drought-relevant
  # "volumetric_soil_water_layer_3",  # 28–100 cm  — medium-term memory  [add if needed]
  # "volumetric_soil_water_layer_4"  , # 100–289 cm — deep storage         [add if needed]
  # "total_evaporation"
)

# Create output directories
dir.create("monthly_data_direct", showWarnings = FALSE)
dir.create("downloads", showWarnings = FALSE)

# ===================== FUNCTION: DOWNLOAD ALL DATA FOR ONE VARIABLE =====================
download_era5_land_variable <- function(variable, start_year, end_year) {
  # Downloads all monthly data for one variable across multiple years
  
  # Generate year and month vectors
  years <- as.character(start_year:end_year)
  months <- sprintf("%02d", 1:12)
  
  # Construct request for monthly means dataset
  request <- list(
    dataset_short_name = "reanalysis-era5-land-monthly-means",
    product_type = if (variable %in% c("total_precipitation", "potential_evaporation",
                                       "total_evaporation",
                                       "surface_solar_radiation_downwards"))
      "monthly_averaged_reanalysis_by_hour_of_day"
    else
      "monthly_averaged_reanalysis",
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
    nc_file <- download_era5_land_variable(var, start_year, end_year)
    
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
    
    # ==========================================================================
    # UNIT CONVERSION + EXPLICIT NETCDF UNITS
    # ==========================================================================
    
    # Keep track of the physical units that will actually be stored in the
    # derived NetCDF. This is critical because raster arithmetic can remove or
    # leave stale NetCDF metadata.
    output_units <- NULL
    
    if (var == "total_precipitation") {
      
      # ERA5-Land monthly-mean accumulated precipitation from this CDS request
      # is interpreted as mean daily accumulation in metres of water equivalent.
      #
      # Convert depth units only:
      #
      #   m/day -> mm/day
      #
      # DO NOT multiply by the number of days in the month here.
      # Monthly totals are calculated downstream.
      result <- result * 1000
      
      output_units <- "mm/day"
      
      cat("  Converted precipitation from m/day to mm/day.\n")
      cat("  NOTE: values remain daily rates; they are NOT mm/month.\n")
    }
    
    if (var == "total_evaporation") {
      
      # Keep ERA5-Land total evaporation in its native representation.
      # Downstream ET processing handles:
      #   1. ECMWF evaporation sign convention
      #   2. m -> mm conversion
      #   3. multiplication by days in each calendar month
      output_units <- "m/day"
      
      cat("  Total evaporation retained in native m/day units.\n")
      cat("  ERA5-Land evaporation sign convention retained.\n")
    }
    
    if (var == "potential_evaporation") {
      
      # Keep native ERA5-Land potential evaporation representation.
      output_units <- "m/day"
      
      cat("  Potential evaporation retained in native m/day units.\n")
    }
    
    if (grepl("temperature", var)) {
      
      # No conversion currently applied.
      output_units <- "K"
      
      # To convert in the future:
      # result <- result - 273.15
      # output_units <- "degC"
    }
    
    if (var == "snow_depth_water_equivalent") {
      
      # Native ERA5-Land SWE.
      output_units <- "m"
      
      cat(
        "  SWE: units are metres water equivalent ",
        "(no conversion applied).\n",
        sep = ""
      )
    }
    
    if (grepl("volumetric_soil_water", var)) {
      
      # Native ERA5-Land volumetric soil-water fraction.
      output_units <- "m3/m3"
      
      cat(
        "  Soil moisture: units are m3/m3 ",
        "(no conversion applied).\n",
        sep = ""
      )
    }
    
    
    # ==========================================================================
    # SAVE DERIVED NETCDF
    # ==========================================================================
    
    output_file <- paste0(
      "monthly_data_direct/",
      var,
      "_monthly.nc"
    )
    
    cat("  Saving to:", output_file, "... ")
    
    terra::writeCDF(
      result,
      output_file,
      overwrite = TRUE
    )
    
    cat("Done.\n")
    
    
    # ==========================================================================
    # WRITE EXPLICIT NETCDF METADATA
    #
    # Do this AFTER writeCDF(). The NetCDF is reopened in write mode and the
    # physical units of the values actually stored in the file are explicitly
    # recorded.
    # ==========================================================================
    
    if (!is.null(output_units)) {
      
      nc <- ncdf4::nc_open(
        output_file,
        write = TRUE
      )
      
      all_nc_vars <- names(nc$var)
      
      # terra::writeCDF() adds a dimensionless grid-mapping variable
      # (typically named "crs") alongside the actual data variable to store
      # the CRS/projection. That variable has zero dimensions, whereas the
      # real data variable is indexed by lon/lat/time. Filter it out instead
      # of assuming there is exactly one variable in the file.
      nc_vars <- all_nc_vars[
        sapply(all_nc_vars, function(v) nc$var[[v]]$ndims > 0)
      ]
      
      if (length(nc_vars) != 1L) {
        
        ncdf4::nc_close(nc)
        
        stop(
          "Expected exactly one data variable (excluding CRS/grid-mapping ",
          "variables) in ",
          output_file,
          ", but found: ",
          paste(nc_vars, collapse = ", "),
          " (all variables in file: ",
          paste(all_nc_vars, collapse = ", "),
          ")"
        )
      }
      
      nc_var <- nc_vars[[1L]]
      
      # Physical units of the values stored in this derived file.
      ncdf4::ncatt_put(
        nc,
        nc_var,
        "units",
        output_units
      )
      
      # Record the project processing convention.
      if (var == "total_precipitation") {
        
        ncdf4::ncatt_put(
          nc,
          nc_var,
          "project_conversion",
          paste(
            "ERA5-Land precipitation converted from m/day to mm/day",
            "by multiplying by 1000.",
            "Values are not monthly totals."
          )
        )
        
        ncdf4::ncatt_put(
          nc,
          nc_var,
          "temporal_representation",
          "monthly mean of daily accumulated precipitation"
        )
      }
      
      if (var %in% c(
        "total_evaporation",
        "potential_evaporation"
      )) {
        
        ncdf4::ncatt_put(
          nc,
          nc_var,
          "project_conversion",
          "Native ERA5-Land values retained; no unit conversion applied."
        )
        
        ncdf4::ncatt_put(
          nc,
          nc_var,
          "temporal_representation",
          "monthly mean of daily accumulation"
        )
      }
      
      # Global provenance attributes.
      ncdf4::ncatt_put(
        nc,
        0,
        "source_dataset",
        "ERA5-Land monthly means"
      )
      
      ncdf4::ncatt_put(
        nc,
        0,
        "project_output_units",
        output_units
      )
      
      ncdf4::nc_sync(nc)
      ncdf4::nc_close(nc)
      
      
      # ========================================================================
      # VERIFY THAT THE UNITS ATTRIBUTE WAS ACTUALLY WRITTEN
      # ========================================================================
      
      nc_check <- ncdf4::nc_open(output_file)
      
      all_check_vars <- names(nc_check$var)
      
      check_vars <- all_check_vars[
        sapply(all_check_vars, function(v) nc_check$var[[v]]$ndims > 0)
      ]
      
      if (length(check_vars) != 1L) {
        
        ncdf4::nc_close(nc_check)
        
        stop(
          "Metadata verification failed for ",
          output_file,
          ": unexpected number of data variables (excluding CRS/grid-",
          "mapping variables). Found: ",
          paste(check_vars, collapse = ", "),
          " (all variables in file: ",
          paste(all_check_vars, collapse = ", "),
          ")"
        )
      }
      
      check_var <- check_vars[[1L]]
      
      unit_att <- ncdf4::ncatt_get(
        nc_check,
        check_var,
        "units"
      )
      
      ncdf4::nc_close(nc_check)
      
      if (
        !isTRUE(unit_att$hasatt) ||
        !identical(
          trimws(as.character(unit_att$value)),
          output_units
        )
      ) {
        
        actual_units <- if (isTRUE(unit_att$hasatt)) {
          as.character(unit_att$value)
        } else {
          "<missing>"
        }
        
        stop(
          "NetCDF unit verification failed for ",
          output_file,
          ". Expected units='",
          output_units,
          "', but found units='",
          actual_units,
          "'."
        )
      }
      
      cat(
        sprintf(
          "  ✓ NetCDF metadata verified: units = %s\n",
          output_units
        )
      )
    }
    
    # Clean up downloaded file
    #
    # terra keeps a live GDAL file handle open on `result`'s source file
    # (nc_file) even after `result` has been reassigned by arithmetic like
    # `result * 1000`. On Windows that handle isn't released until R runs
    # garbage collection, which can make file.remove() silently fail with
    # the file still locked. Drop the raster and force gc() first, then
    # verify the removal actually happened instead of assuming it did.
    rm(result)
    invisible(gc())
    
    removed <- file.remove(nc_file)
    
    if (!isTRUE(removed)) {
      # One retry in case the handle release was still in flight.
      invisible(gc())
      Sys.sleep(1)
      removed <- file.remove(nc_file)
    }
    
    if (!isTRUE(removed)) {
      warning(sprintf(
        "Could not remove downloaded file '%s' (still locked). It will remain in 'downloads/' for manual cleanup.",
        nc_file
      ))
    }
    
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

# ============================================================================
# PART 2 — ERA5 ATMOSPHERIC CIRCULATION + SST MONTHLY VALUES
# ============================================================================

# PART 1 removes downloads/ when it is empty. Recreate it for atmospheric data.
dir.create("downloads", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
#  DOWNLOAD CONTROL
#  Set each flag to TRUE to download & process, FALSE to skip entirely.
#  Example: you already have z500 and SLP on disk from a previous run —
#  set z500 = FALSE, slp = FALSE, sst = TRUE to only fetch SST.
# ============================================================================
RUN <- list(
  #z500 = TRUE,    # 500 hPa Geopotential Height
  #slp  = TRUE,    # Mean Sea Level Pressure
  sst  = TRUE     # Sea Surface Temperature (NE Pacific)
)

# ============================================================================
#  CONFIGURATION
# ============================================================================

# ── Spatial domains ──────────────────────────────────────────────────────────
# ERA5 area format for all API requests: c(North, West, South, East)

# Atmospheric dynamics — captures PNA action centres, Aleutian Low,
# Gulf of Alaska ridge/trough, full synoptic systems driving Nechako drought.
# Justified by TRARE ERA5 map extents ~45-65N, 170-110W (Hurley et al. 2025)
AREA_ATM <- c(70, -180, 40, -100)

# NE Pacific SST — PDO pattern, Marine Heatwave / Blob region, California
# Current. Excludes tropical Pacific (ENSO better represented by ONI directly).
AREA_SST <- c(65, -180, 20, -110)

# ── Time period ──────────────────────────────────────────────────────────────
START_DATE <- "1950-01-01"
END_DATE   <- "2025-12-31"

start_year <- year(as.Date(START_DATE))
end_year   <- year(as.Date(END_DATE))

# ── Climatology reference period ─────────────────────────────────────────────
# WMO standard 1991-2020 — consistent with ECCC normals and SPEI baseline (w1)
CLIM_START <- 1991
CLIM_END   <- 2020

stopifnot(
  "CLIM_START must be >= start_year" = CLIM_START >= start_year,
  "CLIM_END must be <= end_year"     = CLIM_END   <= end_year
)

# ── Variable definitions ──────────────────────────────────────────────────────
atm_variables <- list(
  
  z500 = list(
    run               = RUN$z500,
    name              = "geopotential",
    dataset           = "reanalysis-era5-pressure-levels-monthly-means",
    pressure_level    = "500",
    area              = AREA_ATM,
    target            = paste0("z500_monthly_", start_year, "_", end_year, ".nc"),
    output_name       = "z500_monthly.nc",
    conversion        = "geopotential_to_height",   # m2/s2 -> m
    calculate_anomaly = TRUE,
    description       = "500 hPa Geopotential Height"
  ),
  
  slp = list(
    run               = RUN$slp,
    name              = "mean_sea_level_pressure",
    dataset           = "reanalysis-era5-single-levels-monthly-means",
    pressure_level    = NULL,          # SLP has no pressure dimension
    area              = AREA_ATM,
    target            = paste0("slp_monthly_", start_year, "_", end_year, ".nc"),
    output_name       = "slp_monthly.nc",
    conversion        = "Pa_to_hPa",
    calculate_anomaly = FALSE,         # raw hPa is standard for SLP composites
    description       = "Mean Sea Level Pressure"
  ),
  
  sst = list(
    run               = RUN$sst,
    name              = "sea_surface_temperature",
    dataset           = "reanalysis-era5-single-levels-monthly-means",
    pressure_level    = NULL,          # SST has no pressure dimension
    area              = AREA_SST,      # tighter NE Pacific domain
    target            = paste0("sst_monthly_", start_year, "_", end_year, ".nc"),
    output_name       = "sst_monthly.nc",
    conversion        = "K_to_C",      # ERA5 SST is delivered in Kelvin
    calculate_anomaly = TRUE,
    description       = "Sea Surface Temperature (NE Pacific)"
  )
)

# ============================================================================
#  FUNCTION: BUILD AND SUBMIT CDS REQUEST
# ============================================================================
download_era5_atmos_variable <- function(variable_info, start_yr, end_yr) {
  #
  # Handles both pressure-level and single-level variables.
  # 'area' in variable_info is already c(N, W, S, E) — passed directly to API.
  #
  years  <- as.character(start_yr:end_yr)
  months <- sprintf("%02d", 1:12)
  
  request <- list(
    dataset_short_name = variable_info$dataset,
    product_type       = "monthly_averaged_reanalysis",
    variable           = variable_info$name,
    year               = years,
    month              = months,
    time               = "00:00",
    data_format        = "netcdf",
    download_format    = "unarchived",
    area               = variable_info$area,    # c(N, W, S, E)
    target             = variable_info$target
  )
  
  # Pressure-levels dataset requires an explicit pressure_level field
  if (!is.null(variable_info$pressure_level)) {
    request$pressure_level <- variable_info$pressure_level
  }
  
  cat(sprintf("  Submitting CDS request → target: %s\n", variable_info$target))
  
  file_path <- wf_request(
    user     = CDS_USER,
    request  = request,
    transfer = TRUE,
    path     = "downloads/",
    time_out = 7200
  )
  
  return(file_path)
}

# ============================================================================
#  FUNCTION: NA DIAGNOSTIC REPORT
# ============================================================================
report_na <- function(r, label = "") {
  total_cells <- ncell(r) * nlyr(r)
  na_vals     <- sum(global(r, fun = "isNA", na.rm = FALSE)$isNA)
  na_pct      <- (na_vals / total_cells) * 100
  lbl         <- if (nchar(label) > 0) paste0(" — ", label) else ""
  cat(sprintf("  [NA Check%s]\n", lbl))
  cat(sprintf("    Total cells : %d\n",     total_cells))
  cat(sprintf("    NA count    : %d\n",     na_vals))
  cat(sprintf("    NA percent  : %.4f%%\n", na_pct))
  if (na_vals > 0) cat("    ⚠ WARNING: Missing values detected.\n")
  else             cat("    ✓ No missing values.\n")
}

# ============================================================================
#  FUNCTION: FULL-PERIOD ANOMALY
#  Climatology = mean across ALL downloaded years for each calendar month.
#  Non-standard but useful for internal comparisons and flexible exploration.
# ============================================================================
calculate_anomaly_full_period <- function(raster_data, n_years) {
  n_layers     <- nlyr(raster_data)
  anomaly_list <- vector("list", 12)
  
  cat(sprintf("    Climatology base: all %d years in record\n", n_years))
  
  for (m in 1:12) {
    idx         <- seq(m, n_layers, by = 12)
    month_data  <- subset(raster_data, idx)
    climatology <- app(month_data, mean, na.rm = TRUE)
    anomaly_list[[m]] <- month_data - climatology
  }
  
  combined    <- do.call(c, anomaly_list)
  layer_order <- order(c(sapply(1:12, function(m) seq(m, n_layers, by = 12))))
  combined    <- combined[[layer_order]]
  names(combined) <- names(raster_data)
  return(combined)
}

# ============================================================================
#  FUNCTION: FIXED-PERIOD ANOMALY (1991-2020 WMO standard)
#  Climatology = mean of reference period only; applied to ALL years.
#  Reproducible and comparable to ECCC normals and published literature.
# ============================================================================
calculate_anomaly_fixed_period <- function(raster_data,
                                           all_years,
                                           clim_start = 1991,
                                           clim_end   = 2020) {
  n_layers     <- nlyr(raster_data)
  clim_mask    <- all_years >= clim_start & all_years <= clim_end
  n_clim_years <- sum(clim_mask)
  
  cat(sprintf("    Climatology base: %d-%d (%d years)\n",
              clim_start, clim_end, n_clim_years))
  
  if (n_clim_years < 10)
    warning("Fewer than 10 climatology years — results may be unreliable.")
  
  anomaly_list <- vector("list", 12)
  
  for (m in 1:12) {
    all_idx     <- seq(m, n_layers, by = 12)
    clim_idx    <- all_idx[clim_mask]
    climatology <- app(subset(raster_data, clim_idx), mean, na.rm = TRUE)
    anomaly_list[[m]] <- subset(raster_data, all_idx) - climatology
  }
  
  combined    <- do.call(c, anomaly_list)
  layer_order <- order(c(sapply(1:12, function(m) seq(m, n_layers, by = 12))))
  combined    <- combined[[layer_order]]
  names(combined) <- names(raster_data)
  return(combined)
}

# ============================================================================
#  HELPER: RUN BOTH ANOMALIES AND SAVE TO DISK
#  Called for any variable with calculate_anomaly = TRUE
# ============================================================================
run_anomalies_and_save <- function(result, var_info, all_years_vec) {
  
  n_years_total <- length(all_years_vec)
  yr0           <- all_years_vec[1]
  yr1           <- tail(all_years_vec, 1)
  base_name     <- sub("\\.nc$", "", var_info$output_name)
  
  # ── Full-period anomaly ───────────────────────────────────────────────────
  cat(sprintf("\n  [Anomaly 1/2] Full-period (%d–%d)\n", yr0, yr1))
  anom_full <- calculate_anomaly_full_period(result, n_years_total)
  report_na(anom_full, label = "full-period anomaly")
  
  anom_full_file <- file.path(
    "monthly_data_direct",
    sprintf("%s_anomaly_full_%d_%d.nc", base_name, yr0, yr1)
  )
  cat(sprintf("  Saving → %s... ", anom_full_file))
  writeCDF(anom_full, anom_full_file, overwrite = TRUE)
  cat("Done.\n")
  
  # ── Fixed-period anomaly (WMO 1991-2020) ─────────────────────────────────
  cat(sprintf("\n  [Anomaly 2/2] Fixed-period (%d–%d WMO standard)\n",
              CLIM_START, CLIM_END))
  anom_fixed <- calculate_anomaly_fixed_period(
    raster_data = result,
    all_years   = all_years_vec,
    clim_start  = CLIM_START,
    clim_end    = CLIM_END
  )
  report_na(anom_fixed, label = "fixed-period anomaly")
  
  anom_fixed_file <- file.path(
    "monthly_data_direct",
    sprintf("%s_anomaly_clim%d_%d.nc", base_name, CLIM_START, CLIM_END)
  )
  cat(sprintf("  Saving → %s... ", anom_fixed_file))
  writeCDF(anom_fixed, anom_fixed_file, overwrite = TRUE)
  cat("Done.\n")
  
  invisible(list(full = anom_full_file, fixed = anom_fixed_file))
}

# ============================================================================
#  MAIN PROCESSING LOOP
# ============================================================================
cat("============================================================\n")
cat(" ERA5 Download — Atmospheric Circulation + NE Pacific SST\n")
cat("============================================================\n")
cat(sprintf(" Period      : %d – %d\n", start_year, end_year))
cat(sprintf(" AREA_ATM    : N=%d W=%d S=%d E=%d\n",
            AREA_ATM[1], AREA_ATM[2], AREA_ATM[3], AREA_ATM[4]))
cat(sprintf(" AREA_SST    : N=%d W=%d S=%d E=%d\n",
            AREA_SST[1], AREA_SST[2], AREA_SST[3], AREA_SST[4]))
cat(sprintf(" Clim ref    : %d – %d (WMO)\n", CLIM_START, CLIM_END))
cat("------------------------------------------------------------\n")
cat(" Download plan:\n")
for (vname in names(atm_variables)) {
  v <- atm_variables[[vname]]
  status <- if (isTRUE(v$run)) "WILL DOWNLOAD" else "SKIP"
  cat(sprintf("   %-5s  %-38s  %s\n",
              toupper(vname), v$description, status))
}
cat("============================================================\n\n")

all_years_vec <- start_year:end_year
output_log    <- list()

tryCatch({
  
  for (vname in names(atm_variables)) {
    
    var_info <- atm_variables[[vname]]
    
    # ── SKIP check ───────────────────────────────────────────────────────────
    if (!isTRUE(var_info$run)) {
      cat(sprintf("⏭  SKIPPING %s (%s)  [RUN$%s = FALSE]\n\n",
                  toupper(vname), var_info$description, vname))
      next
    }
    
    cat(sprintf("▶  Processing: %s — %s\n", toupper(vname), var_info$description))
    cat("------------------------------------------------------------\n")
    cat(sprintf("   Dataset : %s\n", var_info$dataset))
    if (!is.null(var_info$pressure_level))
      cat(sprintf("   Level   : %s hPa\n", var_info$pressure_level))
    else
      cat("   Level   : single-level (surface)\n")
    cat(sprintf("   Domain  : N=%d W=%d S=%d E=%d\n",
                var_info$area[1], var_info$area[2],
                var_info$area[3], var_info$area[4]))
    cat(sprintf("   Anomaly : %s\n",
                ifelse(isTRUE(var_info$calculate_anomaly),
                       "YES — full-period + fixed 1991-2020", "NO")))
    cat("------------------------------------------------------------\n")
    
    # ── 1. Download ──────────────────────────────────────────────────────────
    nc_file <- download_era5_atmos_variable(var_info, start_year, end_year)
    cat(sprintf("  ✓ Downloaded: %s\n", nc_file))
    
    # ── 2. Read ───────────────────────────────────────────────────────────────
    cat("  Reading NetCDF... ")
    result <- rast(nc_file)
    cat(sprintf("Done. (%d layers, %d x %d cells)\n",
                nlyr(result), nrow(result), ncol(result)))
    
    # ── 3. NA check on raw download ──────────────────────────────────────────
    report_na(result, label = "raw download")
    
    # ── 4. Unit conversions ───────────────────────────────────────────────────
    if (var_info$conversion == "geopotential_to_height") {
      cat("  Converting geopotential (m2/s2) to geopotential height (m)... ")
      result <- result / 9.80665
      cat("Done.\n")
      
    } else if (var_info$conversion == "Pa_to_hPa") {
      cat("  Converting SLP: Pa to hPa... ")
      result <- result / 100
      cat("Done.\n")
      
    } else if (var_info$conversion == "K_to_C") {
      # ERA5 SST is in Kelvin. Convert to Celsius for readability.
      # Anomaly values are identical in K or C (constant offset cancels),
      # but absolute SST maps are far more interpretable in Celsius.
      cat("  Converting SST: K to degrees C (subtract 273.15)... ")
      result <- result - 273.15
      cat("Done.\n")
    }
    
    # ── 5. Save raw converted data ────────────────────────────────────────────
    raw_file <- file.path("monthly_data_direct", var_info$output_name)
    cat(sprintf("  Saving raw converted data to %s... ", raw_file))
    writeCDF(result, raw_file, overwrite = TRUE)
    cat("Done.\n")
    output_log[[vname]] <- list(raw = raw_file)
    
    # ── 6. Anomaly calculations ───────────────────────────────────────────────
    if (isTRUE(var_info$calculate_anomaly)) {
      anom_files <- run_anomalies_and_save(result, var_info, all_years_vec)
      output_log[[vname]]$anom_full  <- anom_files$full
      output_log[[vname]]$anom_fixed <- anom_files$fixed
    }
    
    # ── 7. Clean up server download file ─────────────────────────────────────
    file.remove(nc_file)
    cat(sprintf("\n  ✓ Completed: %s\n\n", var_info$description))
  }
  
  # ── Cleanup downloads folder ──────────────────────────────────────────────
  cat("------------------------------------------------------------\n")
  cat(" Cleanup\n")
  if (dir.exists("downloads")) {
    remaining <- length(list.files("downloads"))
    if (remaining == 0) {
      unlink("downloads", recursive = TRUE)
      cat("✓ Deleted empty downloads folder.\n")
    } else {
      warning(sprintf(
        "downloads/ still has %d file(s) — skipping deletion for safety.",
        remaining))
    }
  }
  
  # ── Final summary ─────────────────────────────────────────────────────────
  cat("\n============================================================\n")
  cat(" Output files written to: monthly_data_direct/\n")
  cat("============================================================\n")
  for (vname in names(output_log)) {
    fl <- output_log[[vname]]
    cat(sprintf("\n  %s:\n", toupper(vname)))
    cat(sprintf("    Raw             : %s\n", basename(fl$raw)))
    if (!is.null(fl$anom_full))
      cat(sprintf("    Full anomaly    : %s\n", basename(fl$anom_full)))
    if (!is.null(fl$anom_fixed))
      cat(sprintf("    Fixed anomaly   : %s\n", basename(fl$anom_fixed)))
  }
  cat("============================================================\n")
  
}, error = function(e) {
  cat("\n ERROR: Processing failed!\n")
  cat("  Message:", conditionMessage(e), "\n")
  cat("\n  downloads/ folder preserved for debugging.\n")
  stop(e)
})