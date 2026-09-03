# ==============================================================================
# ERA5-Land download and low-storage preprocessing pipeline (merged version)
#
# Purpose
#   1) Preserve the existing monthly ERA5-Land state-variable archive and native PEV product.
#   2) Preserve the optional CDS hourly-download branch for FAO-56 wind/radiation.
#   3) Convert locally available daily GeoTIFF archives (for example from GEE)
#      into annual NetCDF forcing files that are digestible by 2b_PET_ERALand_v2.R.
#   4) Preserve the circulation-context downloader for Z500 and IVT.
#
# Local GeoTIFF contract
#   Expected source files in monthly_data_direct/ (or GEE_DAILY_TIF_DIR):
#     - u2_mean_daily_YYYY.tif
#     - Rs_daily_YYYY.tif
#   Each file must contain one band per day, with date-stamped band descriptions
#   such as u2_mean_19500101 or Rs_19500101.
#
# Output contract
#   daily_eto_forcing/u2_mean_daily_YYYY.nc
#   daily_eto_forcing/Rs_daily_YYYY.nc
#
# These NetCDF files preserve a valid daily time axis and the layer order needed
# by 2b_PET_ERALand_v2.R. All retained gridded data are NetCDF.
# ==============================================================================

library(ecmwfr)
library(terra)
library(lubridate)
library(ncdf4)
library(dotenv)

dotenv::load_dot_env(".env")
setwd(Sys.getenv("WD_PATH"))

# ============================== CONFIGURATION ================================

wf_set_key(
  user = Sys.getenv("CDS_USER"),
  key  = Sys.getenv("CDS_KEY")
)

cds_user <- Sys.getenv("CDS_USER")
if (!nzchar(cds_user)) {
  stop("CDS_USER is missing from .env")
}

# British Columbia / Nechako-region bounding box
# c(xmin, ymin, xmax, ymax)
bbox <- c(-127.8, 52.9, -122.7, 56.2)

start_date <- as.Date("1950-01-01")
end_date   <- as.Date("2025-12-31")

# Existing monthly state variables retained for basin-memory / climate analyses.
# Wind components are intentionally omitted from the monthly archive because the
# FAO-56 forcing branch now derives them from daily products before PET uses them.
monthly_variables <- c(
  "snow_depth_water_equivalent",
  "surface_pressure",
  "2m_temperature",
  "2m_dewpoint_temperature",
  "volumetric_soil_water_layer_1",
  "volumetric_soil_water_layer_2",
  "potential_evaporation"
)

# Minimum hourly forcing set required by the FAO-56 framework.
hourly_variables <- c(
  "10m_u_component_of_wind",
  "10m_v_component_of_wind",
  "surface_solar_radiation_downwards"
)

# 10-m -> 2-m FAO-56 logarithmic wind-height conversion factor.
wind_10m_to_2m_factor <- 4.87 / log(67.8 * 10 - 5.42)

# Output locations.
monthly_dir  <- "monthly_data_direct"
daily_dir    <- "daily_eto_forcing"
download_dir <- "downloads"
carry_dir    <- file.path(download_dir, "carry")
gee_tif_dir  <- normalizePath(
  Sys.getenv("GEE_DAILY_TIF_DIR", file.path(getwd(), "monthly_data_direct")),
  winslash = "/", mustWork = FALSE
)

dir.create(monthly_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(daily_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(carry_dir,    showWarnings = FALSE, recursive = TRUE)

# NetCDF compression level for retained files.
nc_compression <- 9L

# ================================ UTILITIES ==================================

month_sequence <- function(start_date, end_date) {
  seq(
    floor_date(start_date, unit = "month"),
    floor_date(end_date, unit = "month"),
    by = "month"
  )
}

month_days <- function(month_date) {
  sprintf("%02d", seq_len(days_in_month(month_date)))
}

all_hours <- sprintf("%02d:00", 0:23)

safe_remove <- function(path) {
  if (length(path) == 1L && !is.na(path) && file.exists(path)) {
    ok <- file.remove(path)
    if (!ok) warning("Could not delete temporary file: ", path)
  }
  invisible(TRUE)
}

check_time_axis <- function(x, label) {
  tt <- time(x)
  if (is.null(tt) || length(tt) != nlyr(x) || anyNA(tt)) {
    stop("Missing/invalid NetCDF time axis for ", label)
  }
  as.POSIXct(tt, tz = "UTC")
}

report_na <- function(x, label) {
  cat("\n  --- NA Check:", label, "---\n")
  na_by_layer <- global(is.na(x), "sum", na.rm = TRUE)
  total_na <- sum(na_by_layer[, 1], na.rm = TRUE)
  total_cells <- as.double(ncell(x)) * nlyr(x)
  na_percent <- if (total_cells > 0) 100 * total_na / total_cells else NA_real_
  cat(sprintf("  Total cells: %.0f\n", total_cells))
  cat(sprintf("  Total NAs:   %.0f\n", total_na))
  cat(sprintf("  NA percent:  %.6f%%\n", na_percent))
  if (total_na > 0) warning(label, " contains missing values.")
  invisible(total_na)
}

# Robustly extract a variable from a NetCDF that may expose variables either as
# subdatasets or as a single raster with variable-prefixed layer names.
read_nc_variable <- function(nc_file, candidates) {
  candidates <- unique(tolower(candidates))
  
  s <- try(sds(nc_file), silent = TRUE)
  if (!inherits(s, "try-error") && length(s) > 0L) {
    s_names <- tolower(names(s))
    idx <- which(s_names %in% candidates)
    if (length(idx) == 0L) {
      idx <- which(vapply(
        s_names,
        function(z) any(vapply(candidates, function(cn) grepl(cn, z, fixed = TRUE), logical(1))),
        logical(1)
      ))
    }
    if (length(idx) > 0L) return(s[[idx[1]]])
  }
  
  r <- rast(nc_file)
  layer_names <- tolower(names(r))
  idx <- which(vapply(
    layer_names,
    function(z) any(vapply(candidates, function(cn) grepl(cn, z, fixed = TRUE), logical(1))),
    logical(1)
  ))
  
  if (length(idx) == 0L) {
    stop(
      "Could not identify variable in NetCDF. Candidates: ",
      paste(candidates, collapse = ", "),
      ". Available layers/subdatasets should be inspected in: ", nc_file
    )
  }
  
  r[[idx]]
}

# Append a SpatRaster to an annual NetCDF by rewriting the compact daily annual
# file. At this study-domain size, the annual daily product is small; this avoids
# retaining any hourly data and keeps the implementation robust in terra.
nc_time_origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")

spatraster_to_nc_array <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  nl <- nlyr(x)
  
  if (nr <= 0L || nc <= 0L || nl <= 0L) {
    stop("Cannot write an empty SpatRaster to NetCDF.")
  }
  
  values_mat <- terra::values(x, mat = TRUE)
  if (is.null(dim(values_mat))) {
    values_mat <- matrix(values_mat, ncol = 1L)
  }
  
  if (nrow(values_mat) != nr * nc) {
    stop("Unexpected raster value matrix shape while preparing NetCDF output.")
  }
  
  arr <- array(NA_real_, dim = c(nc, nr, nl))
  for (k in seq_len(nl)) {
    layer_mat <- matrix(values_mat[, k], nrow = nr, ncol = nc, byrow = TRUE)
    arr[, , k] <- t(layer_mat)
  }
  
  arr
}

write_daily_netcdf_from_spatraster <- function(x, dates, out_file, variable_name,
                                               long_name, units, overwrite = TRUE) {
  if (nlyr(x) != length(dates)) {
    stop("Layer/date mismatch while writing ", variable_name)
  }
  if (nlyr(x) == 0L) return(invisible(NULL))
  
  dates <- as.Date(dates)
  if (anyNA(dates)) {
    stop("Invalid dates while writing ", variable_name)
  }
  
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  
  tmp_file <- paste0(out_file, ".tmp.nc")
  if (file.exists(tmp_file)) {
    safe_remove(tmp_file)
  }
  
  if (overwrite && file.exists(out_file)) {
    safe_remove(out_file)
  }
  
  nr <- nrow(x)
  nc <- ncol(x)
  xext <- terra::ext(x)
  xres <- (xext[2] - xext[1]) / nc
  yres <- (xext[4] - xext[3]) / nr
  
  lon <- seq(xext[1] + xres / 2, xext[2] - xres / 2, by = xres)
  lat <- seq(xext[4] - yres / 2, xext[3] + yres / 2, by = -yres)
  time_vals <- as.numeric(difftime(as.POSIXct(dates, tz = "UTC"), nc_time_origin,
                                   units = "days"))
  
  data_arr <- spatraster_to_nc_array(x)
  fill_value <- -9999.0
  data_arr[is.na(data_arr)] <- fill_value
  
  lon_dim <- ncdim_def("lon", "degrees_east", lon, create_dimvar = TRUE)
  lat_dim <- ncdim_def("lat", "degrees_north", lat, create_dimvar = TRUE)
  time_dim <- ncdim_def(
    "time",
    paste0("days since ", format(nc_time_origin, "%Y-%m-%d %H:%M:%S")),
    time_vals,
    create_dimvar = TRUE,
    unlim = FALSE
  )
  
  var_def <- ncvar_def(
    variable_name,
    units,
    list(lon_dim, lat_dim, time_dim),
    missval = fill_value,
    longname = long_name,
    prec = "float",
    compression = nc_compression
  )
  
  nc <- nc_create(tmp_file, vars = list(var_def), force_v4 = TRUE)
  on.exit({
    try(nc_close(nc), silent = TRUE)
    if (file.exists(tmp_file)) {
      file.remove(tmp_file)
    }
  }, add = TRUE)
  
  ncvar_put(nc, var_def, data_arr)
  ncatt_put(nc, 0, "Conventions", "CF-1.8")
  ncatt_put(nc, 0, "title", long_name)
  ncatt_put(nc, 0, "source", "Monthly_ERALand_download_v2")
  ncatt_put(nc, 0, "history", paste("Created", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
  ncatt_put(nc, "lon", "axis", "X")
  ncatt_put(nc, "lon", "long_name", "longitude")
  ncatt_put(nc, "lat", "axis", "Y")
  ncatt_put(nc, "lat", "long_name", "latitude")
  ncatt_put(nc, "time", "axis", "T")
  ncatt_put(nc, "time", "calendar", "standard")
  ncatt_put(nc, "time", "long_name", "time")
  
  nc_close(nc)
  on.exit(NULL, add = TRUE)
  
  if (!file.rename(tmp_file, out_file)) {
    stop("Failed to finalize annual NetCDF: ", out_file)
  }
  
  invisible(out_file)
}

append_daily_to_year <- function(x, dates, variable_name, long_name, units) {
  if (nlyr(x) != length(dates)) {
    stop("Layer/date mismatch while writing ", variable_name)
  }
  if (nlyr(x) == 0L) return(invisible(NULL))
  
  years <- unique(year(dates))
  
  for (yy in years) {
    keep <- which(year(dates) == yy)
    x_new <- x[[keep]]
    d_new <- as.Date(dates[keep])
    out_file <- file.path(daily_dir, sprintf("%s_daily_%04d.nc", variable_name, yy))
    
    if (file.exists(out_file)) {
      x_old <- rast(out_file)
      t_old <- check_time_axis(x_old, out_file)
      d_old <- as.Date(t_old, tz = "UTC")
      
      x_all <- c(x_old, x_new)
      d_all <- c(d_old, d_new)
      
      keep_unique <- !duplicated(d_all, fromLast = TRUE)
      x_all <- x_all[[keep_unique]]
      d_all <- d_all[keep_unique]
      
      ord <- order(d_all)
      x_all <- x_all[[ord]]
      d_all <- d_all[ord]
    } else {
      x_all <- x_new
      d_all <- d_new
    }
    
    write_daily_netcdf_from_spatraster(
      x_all,
      d_all,
      out_file,
      variable_name = variable_name,
      long_name = long_name,
      units = units,
      overwrite = TRUE
    )
  }
  
  invisible(NULL)
}

# ========================== MONTHLY STATE DOWNLOAD ============================

download_monthly_variable <- function(variable, start_year, end_year) {
  # ECMWF documents an accumulated-field error in the standard
  # monthly_averaged_reanalysis product for September 2022-February 2024.
  # For potential evaporation (pev), use the corrected synoptic monthly mean
  # at 00 UTC throughout the record. This avoids mixing two definitions of the
  # same accumulated variable in the drought period of interest.
  product_type <- if (identical(variable, "potential_evaporation")) {
    "monthly_averaged_reanalysis_by_hour_of_day"
  } else {
    "monthly_averaged_reanalysis"
  }

  request <- list(
    dataset_short_name = "reanalysis-era5-land-monthly-means",
    product_type = product_type,
    variable = variable,
    year = as.character(start_year:end_year),
    month = sprintf("%02d", 1:12),
    time = "00:00",
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(bbox[4], bbox[1], bbox[2], bbox[3])
  )
  
  wf_request(
    user = cds_user,
    request = request,
    transfer = TRUE,
    path = download_dir,
    time_out = 7200
  )
}

process_monthly_state_variables <- function() {
  start_year <- year(start_date)
  end_year <- year(end_date)
  
  cat("\n============================================================\n")
  cat("PHASE 1: ERA5-Land monthly state variables\n")
  cat("============================================================\n")
  
  for (var in monthly_variables) {
    output_file <- file.path(monthly_dir, paste0(var, "_monthly.nc"))
    
    if (file.exists(output_file)) {
      cat("Skipping existing monthly file:", output_file, "\n")
      next
    }
    
    cat("\nProcessing monthly variable:", var, "\n")
    nc_file <- download_monthly_variable(var, start_year, end_year)
    
    result <- rast(nc_file)
    report_na(result, var)
    
    if (var == "snow_depth_water_equivalent") {
      cat("  SWE retained in metres water equivalent.\n")
    } else if (grepl("volumetric_soil_water", var, fixed = TRUE)) {
      cat("  Soil moisture retained in m3 m-3.\n")
    } else if (var == "surface_pressure") {
      cat("  Surface pressure retained in Pa.\n")
    } else if (grepl("temperature", var, fixed = TRUE)) {
      cat("  Temperature retained in K.\n")
    } else if (identical(var, "potential_evaporation")) {
      cat("  Potential evaporation retained from ERA5-Land PEV (m/day in monthly means); 2b converts to mm/day.\n")
    }
    
    writeCDF(
      result,
      output_file,
      overwrite = TRUE,
      compression = nc_compression
    )
    
    safe_remove(nc_file)
    rm(result)
    gc(verbose = FALSE)
    
    cat("  Saved:", output_file, "\n")
  }
}

# ============================ HOURLY DOWNLOAD ================================

download_hourly_month <- function(month_date) {
  request <- list(
    dataset_short_name = "reanalysis-era5-land",
    variable = hourly_variables,
    year = format(month_date, "%Y"),
    month = format(month_date, "%m"),
    day = month_days(month_date),
    time = all_hours,
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(bbox[4], bbox[1], bbox[2], bbox[3])
  )
  
  wf_request(
    user = cds_user,
    request = request,
    transfer = TRUE,
    path = download_dir,
    time_out = 7200
  )
}

# ============================== DAILY WIND ===================================

process_daily_wind <- function(nc_file, month_date) {
  u10 <- read_nc_variable(
    nc_file,
    c("10m_u_component_of_wind", "u10", "10u")
  )
  v10 <- read_nc_variable(
    nc_file,
    c("10m_v_component_of_wind", "v10", "10v")
  )
  
  tu <- check_time_axis(u10, "u10")
  tv <- check_time_axis(v10, "v10")
  
  if (length(tu) != length(tv) || any(tu != tv)) {
    stop("u10 and v10 time axes are not identical for ", format(month_date, "%Y-%m"))
  }
  
  # Hourly scalar wind magnitude MUST be calculated before temporal averaging.
  ws10 <- sqrt(u10^2 + v10^2)
  u2_hourly <- ws10 * wind_10m_to_2m_factor
  dates_hourly <- as.Date(tu, tz = "UTC")
  
  u2_daily <- tapp(u2_hourly, index = dates_hourly, fun = mean, na.rm = TRUE)
  daily_dates <- sort(unique(dates_hourly))
  
  # Respect user-requested overall date interval.
  keep <- daily_dates >= start_date & daily_dates <= end_date
  u2_daily <- u2_daily[[which(keep)]]
  daily_dates <- daily_dates[keep]
  
  report_na(u2_daily, paste0("u2 daily ", format(month_date, "%Y-%m")))
  
  append_daily_to_year(
    u2_daily,
    daily_dates,
    variable_name = "u2_mean",
    long_name = "Daily mean 2-m wind speed derived from hourly ERA5-Land 10-m u/v components",
    units = "m s-1"
  )
  
  rm(u10, v10, ws10, u2_hourly, u2_daily)
  gc(verbose = FALSE)
}

# ============================= DAILY SOLAR ===================================

# Save the current month's final daily Rs temporarily because its required
# accumulation appears at 00 UTC on the first day of the NEXT month.
save_rs_carry <- function(x, date_value) {
  carry_file <- file.path(carry_dir, sprintf("Rs_carry_%s.nc", format(date_value, "%Y%m%d")))
  write_daily_netcdf_from_spatraster(
    x,
    as.Date(date_value),
    carry_file,
    variable_name = "Rs",
    long_name = "Daily incoming solar radiation",
    units = "MJ m-2 day-1",
    overwrite = TRUE
  )
  carry_file
}

flush_previous_rs_carry <- function(ssrd, ssrd_times, month_date) {
  # The first 00 UTC layer of this month is the 24-h accumulation for the final
  # UTC day of the previous month.
  first_midnight <- as.POSIXct(floor_date(month_date, "month"), tz = "UTC")
  idx <- which(ssrd_times == first_midnight)
  
  if (length(idx) != 1L) {
    stop("Expected exactly one 00 UTC radiation layer at ", first_midnight)
  }
  
  previous_date <- as.Date(first_midnight, tz = "UTC") - 1L
  
  if (previous_date < start_date || previous_date > end_date) {
    return(invisible(NULL))
  }
  
  # J m-2 accumulated over the previous UTC day -> MJ m-2 day-1.
  rs_prev <- ssrd[[idx]] / 1e6
  
  append_daily_to_year(
    rs_prev,
    previous_date,
    variable_name = "Rs",
    long_name = "Daily incoming solar radiation from ERA5-Land SSRD",
    units = "MJ m-2 day-1"
  )
  
  # If an earlier carry marker exists, it is now complete and can be removed.
  old_carry <- file.path(carry_dir, sprintf("Rs_carry_%s.nc", format(previous_date, "%Y%m%d")))
  safe_remove(old_carry)
  
  rm(rs_prev)
  gc(verbose = FALSE)
  invisible(NULL)
}

process_daily_radiation <- function(nc_file, month_date) {
  ssrd <- read_nc_variable(
    nc_file,
    c("surface_solar_radiation_downwards", "ssrd")
  )
  tt <- check_time_axis(ssrd, "ssrd")
  
  # First complete the previous month's final day using this month's 00 UTC.
  flush_previous_rs_carry(ssrd, tt, month_date)
  
  # For dates within this month, timestamp D 00:00 belongs to D-1. Therefore,
  # 00 UTC values on days 2..last provide daily totals for days 1..last-1.
  is_midnight <- format(tt, "%H:%M:%S", tz = "UTC") == "00:00:00"
  idx_midnight <- which(is_midnight)
  
  if (length(idx_midnight) == 0L) {
    stop("No 00 UTC SSRD layers found for ", format(month_date, "%Y-%m"))
  }
  
  rs_dates <- as.Date(tt[idx_midnight], tz = "UTC") - 1L
  in_current_month <- floor_date(rs_dates, "month") == floor_date(month_date, "month")
  
  idx_use <- idx_midnight[in_current_month]
  rs_dates <- rs_dates[in_current_month]
  
  if (length(idx_use) > 0L) {
    rs_daily <- ssrd[[idx_use]] / 1e6
    
    keep <- rs_dates >= start_date & rs_dates <= end_date
    if (any(keep)) {
      rs_daily <- rs_daily[[which(keep)]]
      rs_dates_keep <- rs_dates[keep]
      
      report_na(rs_daily, paste0("Rs daily ", format(month_date, "%Y-%m")))
      
      append_daily_to_year(
        rs_daily,
        rs_dates_keep,
        variable_name = "Rs",
        long_name = "Daily incoming solar radiation from ERA5-Land SSRD",
        units = "MJ m-2 day-1"
      )
    }
    
    rm(rs_daily)
  }
  
  # Create a tiny marker/carry file for the month's final day. The actual value
  # will be supplied by next month's first 00 UTC layer. This file is deliberately
  # small and makes interrupted runs explicit/recoverable.
  last_date <- ceiling_date(month_date, "month") - days(1)
  if (last_date >= start_date && last_date <= end_date) {
    template <- ssrd[[1]]
    values(template) <- NA_real_
    save_rs_carry(template, as.Date(last_date))
    rm(template)
  }
  
  rm(ssrd)
  gc(verbose = FALSE)
}

# Final study day needs the next day's 00 UTC accumulation. Download only one
# radiation time step after end_date, rather than an entire additional month.
finalize_last_rs_day <- function() {
  if (end_date >= Sys.Date()) {
    warning("end_date is not safely historical; final SSRD accumulation may be unavailable.")
  }
  
  carry_file <- file.path(carry_dir, sprintf("Rs_carry_%s.nc", format(end_date, "%Y%m%d")))
  if (!file.exists(carry_file)) return(invisible(NULL))
  
  next_day <- end_date + 1L
  cat("\nFinalizing Rs for", as.character(end_date), "using", as.character(next_day), "00 UTC...\n")
  
  request <- list(
    dataset_short_name = "reanalysis-era5-land",
    variable = "surface_solar_radiation_downwards",
    year = format(next_day, "%Y"),
    month = format(next_day, "%m"),
    day = format(next_day, "%d"),
    time = "00:00",
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(bbox[4], bbox[1], bbox[2], bbox[3])
  )
  
  nc_file <- wf_request(
    user = cds_user,
    request = request,
    transfer = TRUE,
    path = download_dir,
    time_out = 7200
  )
  
  ssrd <- read_nc_variable(nc_file, c("surface_solar_radiation_downwards", "ssrd"))
  tt <- check_time_axis(ssrd, "final ssrd")
  
  expected <- as.POSIXct(next_day, tz = "UTC")
  idx <- which(tt == expected)
  if (length(idx) != 1L) stop("Could not identify final 00 UTC SSRD layer.")
  
  rs_last <- ssrd[[idx]] / 1e6
  
  append_daily_to_year(
    rs_last,
    end_date,
    variable_name = "Rs",
    long_name = "Daily incoming solar radiation from ERA5-Land SSRD",
    units = "MJ m-2 day-1"
  )
  
  safe_remove(carry_file)
  safe_remove(nc_file)
  rm(ssrd, rs_last)
  gc(verbose = FALSE)
}

# ========================= HOURLY -> DAILY DRIVER =============================


process_hourly_forcing <- function() {
  months <- month_sequence(start_date, end_date)
  
  cat("\n============================================================\n")
  cat("PHASE 2: Hourly u10/v10/SSRD -> daily u2/Rs\n")
  cat("============================================================\n")
  cat(sprintf("10-m -> 2-m wind factor: %.8f\n", wind_10m_to_2m_factor))
  
  for (month_date in months) {
    label <- format(month_date, "%Y-%m")
    cat("\n------------------------------------------------------------\n")
    cat("Processing hourly month:", label, "\n")
    
    nc_file <- NULL
    
    tryCatch({
      nc_file <- download_hourly_month(month_date)
      cat("  Downloaded:", nc_file, "\n")
      
      process_daily_wind(nc_file, month_date)
      process_daily_radiation(nc_file, month_date)
      
      # Delete hourly NetCDF ONLY after both daily products succeed.
      safe_remove(nc_file)
      nc_file <- NULL
      gc(verbose = FALSE)
      
      cat("  Completed:", label, "\n")
    }, error = function(e) {
      cat("\nERROR while processing ", label, ": ", conditionMessage(e), "\n", sep = "")
      if (!is.null(nc_file) && file.exists(nc_file)) {
        cat("Temporary NetCDF retained for debugging:", nc_file, "\n")
      }
      stop(e)
    })
  }
  
  finalize_last_rs_day()
}


# ==================== LOCAL GEE DAILY TIFF -> NETCDF ===========================
# This branch converts annual daily GeoTIFF stacks (for example GEE exports)
# into annual NetCDF files with a proper daily time axis.  It is the preferred
# path when locally exported daily GeoTIFFs already exist in monthly_data_direct/.
# ==============================================================================
is_leap_year <- function(yr) {
  (yr %% 4L == 0L) & ((yr %% 100L) != 0L | (yr %% 400L) == 0L)
}

find_local_daily_tifs <- function(variable_name, source_dir = gee_tif_dir) {
  if (!dir.exists(source_dir)) return(character(0))
  pat <- paste0("^", variable_name, "_daily_[0-9]{4}\\.tif(f)?$")
  list.files(source_dir, pattern = pat, full.names = TRUE)
}

# ------------------------------------------------------------------------------
# Determine whether the complete local GeoTIFF archive already exists.
# If TRUE, the hourly ERA5-Land download branch can be skipped safely.
# ------------------------------------------------------------------------------
has_complete_local_daily_tifs <- function(source_dir = gee_tif_dir,
                                          start_year = year(start_date),
                                          end_year   = year(end_date)) {
  
  yrs <- start_year:end_year
  
  u_expected <- file.path(
    source_dir,
    sprintf("u2_mean_daily_%04d.tif", yrs)
  )
  
  r_expected <- file.path(
    source_dir,
    sprintf("Rs_daily_%04d.tif", yrs)
  )
  
  missing_u <- u_expected[!file.exists(u_expected)]
  missing_r <- r_expected[!file.exists(r_expected)]
  
  if (length(missing_u) == 0 && length(missing_r) == 0) {
    message("\nUsing local daily GeoTIFF archive.")
    message("Hourly ERA5-Land download will be skipped.")
    return(TRUE)
  }
  
  message("\nLocal GeoTIFF archive is incomplete.")
  if (length(missing_u))
    message("Missing u2 files: ", length(missing_u))
  if (length(missing_r))
    message("Missing Rs files: ", length(missing_r))
  
  FALSE
}
extract_dates_from_band_names <- function(r, tif_file, variable_name) {
  band_names <- names(r)
  if (is.null(band_names) || !length(band_names)) {
    stop("No layer names found in GeoTIFF: ", tif_file)
  }
  
  # Prefer date tokens embedded in band descriptions/names.
  date_tokens <- regmatches(band_names, regexpr("(19|20)[0-9]{6}", band_names))
  has_dates <- lengths(date_tokens) > 0L
  if (!all(has_dates)) {
    stop(
      "Could not extract YYYYMMDD dates from all band names in: ", tif_file,
      "\nExample band names: ", paste(head(band_names, 3), collapse = ", ")
    )
  }
  
  dates <- as.Date(vapply(date_tokens, `[`, character(1), 1L, USE.NAMES = FALSE), format = "%Y%m%d")
  if (anyNA(dates)) {
    stop("Unparseable band dates in GeoTIFF: ", tif_file)
  }
  
  # Sort defensively in case layer order is not chronological.
  ord <- order(dates)
  r <- r[[ord]]
  dates <- dates[ord]
  
  # Enforce a single calendar year per archive.
  tif_year <- suppressWarnings(as.integer(sub(".*_daily_([0-9]{4})\\.tif(f)?$", "\\1", basename(tif_file))))
  if (is.na(tif_year)) {
    stop("Could not parse year from filename: ", basename(tif_file))
  }
  if (length(unique(as.integer(format(dates, "%Y")))) != 1L ||
      unique(as.integer(format(dates, "%Y"))) != tif_year) {
    stop("Band dates do not match file year for: ", basename(tif_file))
  }
  
  expected <- seq(
    as.Date(sprintf("%04d-01-01", tif_year)),
    as.Date(sprintf("%04d-12-31", tif_year)),
    by = "day"
  )
  # Compare by numeric value (days since epoch) rather than identical(): Date
  # vectors built via seq.Date() vs. as.Date(character) are not guaranteed to
  # share the same internal storage mode (integer vs. double) across R
  # versions/platforms, which makes identical() an unreliable equality check
  # here even when every date genuinely matches.
  if (length(dates) != length(expected) ||
      !isTRUE(all(as.numeric(dates) == as.numeric(expected)))) {
    missing_dates <- setdiff(expected, dates)
    extra_dates <- setdiff(dates, expected)
    stop(
      "Daily GeoTIFF does not contain exactly the expected daily calendar for ",
      tif_year, ".\nMissing: ",
      paste(head(missing_dates, 10), collapse = ", "),
      if (length(extra_dates)) paste0("\nExtra: ", paste(head(extra_dates, 10), collapse = ", ")) else ""
    )
  }
  
  list(raster = r, dates = dates, year = tif_year)
}

write_daily_tif_archive_to_nc <- function(tif_file, variable_name, long_name, units,
                                          output_dir = daily_dir) {
  cat("\nConverting local GeoTIFF -> NetCDF:\n  ", basename(tif_file), "\n", sep = "")
  
  x <- rast(tif_file)
  parsed <- extract_dates_from_band_names(x, tif_file, variable_name)
  x <- parsed$raster
  dates <- parsed$dates
  tif_year <- parsed$year
  
  out_file <- file.path(output_dir, sprintf("%s_daily_%04d.nc", variable_name, tif_year))
  
  if (file.exists(out_file)) {
    cat("Skipping existing NetCDF:", basename(out_file), "\n")
    return(invisible(out_file))
  }
  
  write_daily_netcdf_from_spatraster(
    x,
    dates,
    out_file,
    variable_name = variable_name,
    long_name = long_name,
    units = units,
    overwrite = TRUE
  )
  
  # Validate immediately after write.
  x_out <- rast(out_file)
  tt <- check_time_axis(x_out, out_file)
  dd <- as.Date(tt, tz = "UTC")
  if (!identical(dd, dates)) {
    stop(
      "Written NetCDF time axis mismatch for ", basename(out_file),
      "\nExpected first/last dates: ", as.character(min(dates)), " / ", as.character(max(dates)),
      "\nObserved first/last dates: ", as.character(min(dd)), " / ", as.character(max(dd))
    )
  }
  
  cat("  Saved:", out_file, "\n")
  invisible(out_file)
}

process_local_daily_tifs <- function(source_dir = gee_tif_dir) {
  
  cat("\n============================================================\n")
  cat("PHASE 2: Local daily GeoTIFF archives -> annual NetCDF forcing\n")
  cat("============================================================\n")
  cat("GeoTIFF source directory:", source_dir, "\n")
  
  u_files <- find_local_daily_tifs("u2_mean", source_dir)
  r_files <- find_local_daily_tifs("Rs", source_dir)
  
  if (!length(u_files) && !length(r_files)) {
    stop(
      "No local daily GeoTIFFs found in: ", source_dir,
      "\nExpected files such as u2_mean_daily_1950.tif and Rs_daily_1950.tif."
    )
  }
  
  if (length(u_files)) {
    cat("\n-- Processing daily wind archives --\n")
    for (f in sort(u_files)) {
      write_daily_tif_archive_to_nc(
        f,
        variable_name = "u2_mean",
        long_name = "Daily mean 2-m wind speed derived from local GeoTIFF archive",
        units = "m s-1"
      )
    }
  } else {
    warning("No u2_mean_daily_YYYY.tif files found in ", source_dir)
  }
  
  if (length(r_files)) {
    cat("\n-- Processing daily radiation archives --\n")
    for (f in sort(r_files)) {
      write_daily_tif_archive_to_nc(
        f,
        variable_name = "Rs",
        long_name = "Daily incoming solar radiation derived from local GeoTIFF archive",
        units = "MJ m-2 day-1"
      )
    }
  } else {
    warning("No Rs_daily_YYYY.tif files found in ", source_dir)
  }
  
  invisible(TRUE)
}

# Preserve the original CDS hourly branch as a fallback / legacy option.
# The local GeoTIFF conversion branch is the preferred path for the current
# workflow because 2b_PET_ERALand_v2.R consumes annual daily NetCDF forcing.
# ==============================================================================
# ============================== VALIDATION ===================================

validate_daily_archive <- function(variable_name) {
  files <- list.files(
    daily_dir,
    pattern = paste0("^", variable_name, "_daily_[0-9]{4}\\.nc$"),
    full.names = TRUE
  )
  
  if (length(files) == 0L) {
    stop("No daily output files found for ", variable_name)
  }
  
  all_dates <- do.call(c, lapply(files, function(f) {
    x <- rast(f)
    as.Date(check_time_axis(x, f), tz = "UTC")
  }))
  
  all_dates <- sort(unique(all_dates))
  expected <- seq(start_date, end_date, by = "day")
  missing_dates <- setdiff(expected, all_dates)
  extra_dates <- setdiff(all_dates, expected)
  
  cat("\nValidation:", variable_name, "\n")
  cat("  Expected days:", length(expected), "\n")
  cat("  Stored days:  ", length(all_dates), "\n")
  cat("  Missing days: ", length(missing_dates), "\n")
  cat("  Extra days:   ", length(extra_dates), "\n")
  
  if (length(missing_dates) > 0L) {
    warning(
      variable_name, " has missing dates. First missing dates: ",
      paste(head(missing_dates, 10), collapse = ", ")
    )
  }
  if (length(extra_dates) > 0L) {
    warning(variable_name, " has dates outside requested interval.")
  }
  
  invisible(list(missing = missing_dates, extra = extra_dates))
}

# ======================= PHASE 3: CIRCULATION-CONTEXT FIELDS ==================
# MANUSCRIPT ALIGNMENT (Methods 3g, "Circulation Context of the 2022-2025
# Drought"): the large-scale atmospheric context is evaluated using 500-hPa
# geopotential height (ridging indicator) and a moisture-transport diagnostic
# (integrated vapor transport / vertically integrated moisture flux). Neither
# field is part of ERA5-Land (a land-surface reanalysis), so this phase downloads
# the required fields from ERA5 proper.
#
# Output contract required by 14circulation_context_v2.R:
#   monthly_data_direct/geopotential_500hPa_monthly.nc
#   monthly_data_direct/ivt_components_monthly.nc
#
# The circulation domain is deliberately wider than the Nechako basin because
# Methods 3g concerns large-scale circulation over the North Pacific / western
# North America sector.
# ============================================================================== 

circ_bbox <- c(-170, 25, -100, 65)  # xmin, ymin, xmax, ymax
G_GRAV_CONST <- 9.80665

expected_circulation_dates <- function() {
  seq(
    floor_date(start_date, unit = "month"),
    floor_date(end_date, unit = "month"),
    by = "month"
  )
}

download_z500_monthly <- function(start_year, end_year) {
  # ERA5 pressure-level monthly means. CDS stores geopotential in m2 s-2;
  # Phase 3 converts it to geopotential height in metres before retention.
  request <- list(
    dataset_short_name = "reanalysis-era5-pressure-levels-monthly-means",
    product_type = "monthly_averaged_reanalysis",
    variable = "geopotential",
    pressure_level = "500",
    year = as.character(start_year:end_year),
    month = sprintf("%02d", 1:12),
    time = "00:00",
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(circ_bbox[4], circ_bbox[1], circ_bbox[2], circ_bbox[3])
  )
  
  wf_request(
    user = cds_user,
    request = request,
    transfer = TRUE,
    path = download_dir,
    time_out = 7200
  )
}


find_existing_ivt_download <- function() {
  files <- list.files(
    download_dir,
    pattern = "^ecmwfr_.*\\.nc$",
    full.names = TRUE
  )
  if (!length(files)) return(NA_character_)
  
  hits <- vapply(files, function(f) {
    ok <- FALSE
    nc <- try(ncdf4::nc_open(f), silent = TRUE)
    if (!inherits(nc, "try-error")) {
      on.exit(ncdf4::nc_close(nc), add = TRUE)
      vn <- tolower(names(nc$var))
      ok <- any(grepl("vertical_integral_of_eastward|eastward", vn)) &&
        any(grepl("vertical_integral_of_northward|northward", vn))
    }
    ok
  }, logical(1))
  
  candidates <- files[hits]
  if (!length(candidates)) return(NA_character_)
  
  # If an interrupted Phase 3 left multiple CDS files, use the newest matching
  # IVT file; all files are still retained until successful finalization.
  candidates[which.max(file.info(candidates)$mtime)]
}

download_ivt_components_monthly <- function(start_year, end_year) {
  # ERA5 single-level monthly means of vertically integrated eastward and
  # northward water-vapour flux. IVT magnitude is calculated downstream as
  # sqrt(U^2 + V^2) in 14circulation_context_v2.R.
  request <- list(
    dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
    product_type = "monthly_averaged_reanalysis",
    variable = c(
      "vertical_integral_of_eastward_water_vapour_flux",
      "vertical_integral_of_northward_water_vapour_flux"
    ),
    year = as.character(start_year:end_year),
    month = sprintf("%02d", 1:12),
    time = "00:00",
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(circ_bbox[4], circ_bbox[1], circ_bbox[2], circ_bbox[3])
  )
  
  wf_request(
    user = cds_user,
    request = request,
    transfer = TRUE,
    path = download_dir,
    time_out = 7200
  )
}

validate_circulation_record <- function(x, label) {
  tt <- check_time_axis(x, label)
  dates <- as.Date(tt, tz = "UTC")
  expected <- expected_circulation_dates()
  
  if (length(dates) != length(expected) ||
      !isTRUE(all(as.numeric(dates) == as.numeric(expected)))) {
    missing_dates <- setdiff(expected, dates)
    extra_dates <- setdiff(dates, expected)
    stop(
      "Invalid monthly time axis for ", label,
      ". Expected ", length(expected), " months from ",
      as.character(min(expected)), " to ", as.character(max(expected)),
      ". Missing: ", paste(head(missing_dates, 12), collapse = ", "),
      if (length(extra_dates))
        paste0(". Extra: ", paste(head(extra_dates, 12), collapse = ", "))
      else
        ""
    )
  }
  
  cat("  ", label, ": ", length(dates), " months; ",
      as.character(min(dates)), " to ", as.character(max(dates)), "\n", sep = "")
  invisible(dates)
}


# ------------------------------------------------------------------------------
# ERA5 monthly NetCDF time-coordinate reader.
# terra::sds() can expose ERA5 multi-variable NetCDF subdatasets without attaching
# the CF time coordinate to each returned SpatRaster. Phase 3 therefore reads the
# authoritative time coordinate directly with ncdf4 and assigns it explicitly.
# ------------------------------------------------------------------------------
read_netcdf_time_axis <- function(nc_file, expected_n = NULL) {
  nc <- ncdf4::nc_open(nc_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  dim_names <- names(nc$dim)
  if (!length(dim_names)) {
    stop("No dimensions found in NetCDF: ", nc_file)
  }
  
  time_idx <- which(tolower(dim_names) %in% c("time", "valid_time"))
  if (!length(time_idx)) {
    time_idx <- which(vapply(
      nc$dim,
      function(d) grepl("since", tolower(d$units %||% ""), fixed = TRUE),
      logical(1)
    ))
  }
  if (!length(time_idx)) {
    stop("Could not identify a CF time dimension in NetCDF: ", nc_file)
  }
  
  d <- nc$dim[[time_idx[1]]]
  vals <- as.numeric(d$vals)
  units <- as.character(d$units)
  
  if (!length(vals) || !nzchar(units)) {
    stop("NetCDF time dimension has no values/units: ", nc_file)
  }
  
  m <- regexec(
    "^\\s*([A-Za-z]+)\\s+since\\s+(.+?)\\s*$",
    units,
    ignore.case = TRUE
  )
  parts <- regmatches(units, m)[[1]]
  if (length(parts) != 3L) {
    stop("Unsupported CF time units in ", nc_file, ": ", units)
  }
  
  unit_name <- tolower(parts[2])
  origin_txt <- trimws(parts[3])
  
  multiplier <- switch(
    unit_name,
    second = 1,
    seconds = 1,
    sec = 1,
    secs = 1,
    minute = 60,
    minutes = 60,
    min = 60,
    mins = 60,
    hour = 3600,
    hours = 3600,
    hr = 3600,
    hrs = 3600,
    day = 86400,
    days = 86400,
    stop("Unsupported CF time unit '", unit_name, "' in ", nc_file)
  )
  
  origin <- as.POSIXct(origin_txt, tz = "UTC")
  if (is.na(origin)) {
    origin <- as.POSIXct(
      sub("\\.0+$", "", origin_txt),
      tz = "UTC"
    )
  }
  if (is.na(origin)) {
    stop("Could not parse NetCDF time origin '", origin_txt,
         "' in ", nc_file)
  }
  
  tt <- origin + vals * multiplier
  dates <- as.Date(tt, tz = "UTC")
  
  if (!is.null(expected_n) && length(dates) != expected_n) {
    stop(
      "NetCDF time-axis length (", length(dates),
      ") does not match raster layer count (", expected_n, "): ", nc_file
    )
  }
  
  if (anyNA(dates)) {
    stop("NetCDF time axis contains unparseable dates: ", nc_file)
  }
  
  dates
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) y else x
}

# ------------------------------------------------------------------------------
# Identify the two IVT variables from the CDS NetCDF using the NetCDF variable
# metadata, then attach the CF time coordinate explicitly to the terra rasters.
# ------------------------------------------------------------------------------
extract_ivt_components <- function(nc_file) {
  nc <- ncdf4::nc_open(nc_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  var_names <- names(nc$var)
  if (!length(var_names)) {
    stop("No variables found in IVT NetCDF: ", nc_file)
  }
  
  lower_names <- tolower(var_names)
  
  east_idx <- which(grepl("vertical_integral_of_eastward|eastward", lower_names))
  north_idx <- which(grepl("vertical_integral_of_northward|northward", lower_names))
  
  # Fall back to the first two non-coordinate variables only if the CDS variable
  # names are unavailable/changed. The request order is eastward, northward.
  if (!length(east_idx) || !length(north_idx)) {
    east_idx <- integer(0)
    north_idx <- integer(0)
    noncoord <- which(!lower_names %in% c("time", "valid_time", "latitude",
                                          "longitude", "lat", "lon"))
    if (length(noncoord) >= 2L) {
      east_idx <- noncoord[1]
      north_idx <- noncoord[2]
    }
  }
  
  if (!length(east_idx) || !length(north_idx)) {
    stop(
      "Could not identify eastward/northward IVT variables in: ", nc_file,
      "\nAvailable variables: ", paste(var_names, collapse = ", ")
    )
  }
  
  east_name <- var_names[east_idx[1]]
  north_name <- var_names[north_idx[1]]
  
  s <- try(terra::sds(nc_file), silent = TRUE)
  if (inherits(s, "try-error") || length(s) < 2L) {
    stop(
      "Downloaded IVT NetCDF does not expose both variables as terra subdatasets: ",
      nc_file
    )
  }
  
  s_names <- tolower(names(s))
  east_sds <- which(grepl("eastward|east", s_names))
  north_sds <- which(grepl("northward|north", s_names))
  
  if (!length(east_sds) || !length(north_sds)) {
    # CDS request order is eastward followed by northward.
    east_sds <- 1L
    north_sds <- 2L
  }
  
  east <- s[[east_sds[1]]]
  north <- s[[north_sds[1]]]
  
  dates <- read_netcdf_time_axis(nc_file, expected_n = terra::nlyr(east))
  
  if (terra::nlyr(north) != length(dates)) {
    stop("Northward IVT layer count differs from NetCDF time axis.")
  }
  
  terra::time(east) <- as.POSIXct(dates, tz = "UTC")
  terra::time(north) <- as.POSIXct(dates, tz = "UTC")
  
  attr(east, "nc_variable_name") <- east_name
  attr(north, "nc_variable_name") <- north_name
  
  list(
    east = east,
    north = north,
    dates = dates,
    east_name = east_name,
    north_name = north_name
  )
}

# ------------------------------------------------------------------------------
# Write the two IVT components to a compact, CF-compliant NetCDF with one shared
# explicit time dimension. Writing one time slice at a time avoids creating
# multi-hundred-MB in-memory 3-D arrays for the full 1950-2025 domain.
# ------------------------------------------------------------------------------
write_ivt_components_netcdf <- function(east, north, dates, out_file,
                                        chunk_layers = 12L,
                                        compression = 4L) {
  if (terra::nlyr(east) != length(dates) ||
      terra::nlyr(north) != length(dates)) {
    stop("IVT layer/date mismatch while writing: ", out_file)
  }
  
  if (terra::nrow(east) != terra::nrow(north) ||
      terra::ncol(east) != terra::ncol(north)) {
    stop("Eastward and northward IVT row/column dimensions differ.")
  }
  
  # terra::ext() returns an S4 SpatExtent object. Do not coerce it with
  # as.numeric(); that is not a supported coercion in current terra versions.
  east_ext <- terra::ext(east)
  north_ext <- terra::ext(north)
  
  east_ext_values <- c(
    terra::xmin(east), terra::xmax(east),
    terra::ymin(east), terra::ymax(east)
  )
  north_ext_values <- c(
    terra::xmin(north), terra::xmax(north),
    terra::ymin(north), terra::ymax(north)
  )
  
  if (!isTRUE(all.equal(
    east_ext_values,
    north_ext_values,
    tolerance = 1e-8
  ))) {
    stop("Eastward and northward IVT extents differ.")
  }
  
  east_res <- terra::res(east)
  north_res <- terra::res(north)
  if (!isTRUE(all.equal(
    as.numeric(east_res),
    as.numeric(north_res),
    tolerance = 1e-8
  ))) {
    stop("Eastward and northward IVT resolutions differ.")
  }
  
  east_origin <- terra::origin(east)
  north_origin <- terra::origin(north)
  if (!isTRUE(all.equal(
    as.numeric(east_origin),
    as.numeric(north_origin),
    tolerance = 1e-8
  ))) {
    stop("Eastward and northward IVT origins differ.")
  }
  
  east_crs <- terra::crs(east, proj = TRUE)
  north_crs <- terra::crs(north, proj = TRUE)
  if (!identical(east_crs, north_crs)) {
    stop("Eastward and northward IVT coordinate reference systems differ.")
  }
  
  dates <- as.Date(dates)
  if (anyNA(dates) || anyDuplicated(dates)) {
    stop("IVT dates contain NA values or duplicates.")
  }
  
  expected_dates <- expected_circulation_dates()
  if (length(dates) != length(expected_dates) ||
      !isTRUE(all(as.numeric(dates) == as.numeric(expected_dates)))) {
    stop("IVT dates are not the complete expected chronological monthly sequence.")
  }
  
  nlat <- terra::nrow(east)
  nlon <- terra::ncol(east)
  ntime <- length(dates)
  
  ext <- east_ext
  xres <- (ext[2] - ext[1]) / nlon
  yres <- (ext[4] - ext[3]) / nlat
  
  lon <- seq(ext[1] + xres / 2, ext[2] - xres / 2, length.out = nlon)
  lat <- seq(ext[4] - yres / 2, ext[3] + yres / 2, length.out = nlat)
  
  time_vals <- as.numeric(
    difftime(
      as.POSIXct(dates, tz = "UTC"),
      nc_time_origin,
      units = "days"
    )
  )
  
  tmp_file <- paste0(out_file, ".tmp.nc")
  safe_remove(tmp_file)
  safe_remove(out_file)
  
  chunk_layers <- as.integer(chunk_layers)
  compression <- as.integer(compression)
  
  if (length(chunk_layers) != 1L || is.na(chunk_layers) || chunk_layers < 1L) {
    stop("chunk_layers must be a positive integer.")
  }
  if (length(compression) != 1L || is.na(compression) ||
      compression < 0L || compression > 9L) {
    stop("compression must be an integer from 0 to 9.")
  }
  
  lon_dim <- ncdf4::ncdim_def(
    "lon", "degrees_east", lon, create_dimvar = TRUE
  )
  lat_dim <- ncdf4::ncdim_def(
    "lat", "degrees_north", lat, create_dimvar = TRUE
  )
  time_dim <- ncdf4::ncdim_def(
    "time",
    "days since 1970-01-01 00:00:00",
    time_vals,
    create_dimvar = TRUE,
    unlim = FALSE
  )
  
  fill_value <- -9999.0
  
  east_var <- ncdf4::ncvar_def(
    "ivt_eastward",
    "kg m-1 s-1",
    list(lon_dim, lat_dim, time_dim),
    missval = fill_value,
    longname = "Vertically integrated eastward water vapour flux",
    prec = "float",
    compression = compression
  )
  
  north_var <- ncdf4::ncvar_def(
    "ivt_northward",
    "kg m-1 s-1",
    list(lon_dim, lat_dim, time_dim),
    missval = fill_value,
    longname = "Vertically integrated northward water vapour flux",
    prec = "float",
    compression = compression
  )
  
  ncdf <- ncdf4::nc_create(
    tmp_file,
    vars = list(east_var, north_var),
    force_v4 = TRUE
  )
  
  closed <- FALSE
  on.exit({
    if (!closed) try(ncdf4::nc_close(ncdf), silent = TRUE)
    if (file.exists(tmp_file)) safe_remove(tmp_file)
  }, add = TRUE)
  
  ncdf4::ncatt_put(ncdf, 0, "Conventions", "CF-1.8")
  ncdf4::ncatt_put(
    ncdf, 0, "title",
    "ERA5 vertically integrated water vapour transport components"
  )
  ncdf4::ncatt_put(
    ncdf, 0, "source",
    "ERA5; processed by Monthly_ERALand_download_v2"
  )
  ncdf4::ncatt_put(
    ncdf, 0, "history",
    paste("Created", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  )
  
  ncdf4::ncatt_put(ncdf, "lon", "axis", "X")
  ncdf4::ncatt_put(ncdf, "lon", "long_name", "longitude")
  ncdf4::ncatt_put(ncdf, "lat", "axis", "Y")
  ncdf4::ncatt_put(ncdf, "lat", "long_name", "latitude")
  ncdf4::ncatt_put(ncdf, "time", "axis", "T")
  ncdf4::ncatt_put(ncdf, "time", "calendar", "standard")
  ncdf4::ncatt_put(ncdf, "time", "long_name", "time")
  
  # The previous implementation performed 912 terra reads per component and
  # 1,824 compressed ncvar_put() calls. With compression=9 this can appear to
  # stall because every monthly write invokes the compressor separately.
  # Twelve monthly layers are therefore read and written as one block.
  chunk_layers <- min(chunk_layers, ntime)
  n_chunks <- ceiling(ntime / chunk_layers)
  
  cat("  IVT writer:\n")
  cat("    grid: ", nlon, " x ", nlat, "\n", sep = "")
  cat("    months: ", ntime, "\n", sep = "")
  cat("    write chunks: ", n_chunks, " (", chunk_layers,
      " months/chunk)\n", sep = "")
  cat("    NetCDF compression: ", compression, "\n", sep = "")
  flush.console()
  
  for (chunk_id in seq_len(n_chunks)) {
    k1 <- (chunk_id - 1L) * chunk_layers + 1L
    k2 <- min(chunk_id * chunk_layers, ntime)
    kk <- k1:k2
    nk <- length(kk)
    
    east_values <- terra::values(east[[kk]], mat = TRUE)
    north_values <- terra::values(north[[kk]], mat = TRUE)
    
    if (!is.matrix(east_values) ||
        nrow(east_values) != nlat * nlon ||
        ncol(east_values) != nk) {
      stop("Unexpected eastward IVT block shape in chunk ", chunk_id)
    }
    
    if (!is.matrix(north_values) ||
        nrow(north_values) != nlat * nlon ||
        ncol(north_values) != nk) {
      stop("Unexpected northward IVT block shape in chunk ", chunk_id)
    }
    
    east_arr <- array(NA_real_, dim = c(nlon, nlat, nk))
    north_arr <- array(NA_real_, dim = c(nlon, nlat, nk))
    
    for (j in seq_len(nk)) {
      east_mat <- matrix(
        east_values[, j],
        nrow = nlat,
        ncol = nlon,
        byrow = TRUE
      )
      north_mat <- matrix(
        north_values[, j],
        nrow = nlat,
        ncol = nlon,
        byrow = TRUE
      )
      
      east_arr[, , j] <- t(east_mat)
      north_arr[, , j] <- t(north_mat)
    }
    
    east_arr[is.na(east_arr)] <- fill_value
    north_arr[is.na(north_arr)] <- fill_value
    
    ncdf4::ncvar_put(
      ncdf, east_var, east_arr,
      start = c(1L, 1L, k1),
      count = c(nlon, nlat, nk)
    )
    ncdf4::ncvar_put(
      ncdf, north_var, north_arr,
      start = c(1L, 1L, k1),
      count = c(nlon, nlat, nk)
    )
    
    cat(
      sprintf(
        "    chunk %d/%d: %s to %s\n",
        chunk_id, n_chunks,
        format(dates[k1], "%Y-%m"),
        format(dates[k2], "%Y-%m")
      )
    )
    flush.console()
    
    rm(east_values, north_values, east_arr, north_arr)
    gc(verbose = FALSE)
  }
  
  ncdf4::nc_close(ncdf)
  closed <- TRUE
  
  if (!file.rename(tmp_file, out_file)) {
    stop("Failed to finalize IVT NetCDF: ", out_file)
  }
  
  cat("  IVT NetCDF finalized:", out_file, "\n")
  invisible(out_file)
}

process_circulation_context_fields <- function() {
  start_year <- year(start_date)
  end_year <- year(end_date)
  
  cat("\n============================================================\n")
  cat("PHASE 3: ERA5 circulation-context fields (Methods 3g)\n")
  cat("============================================================\n")
  cat("Circulation domain [xmin ymin xmax ymax]: ",
      paste(circ_bbox, collapse = " "), "\n", sep = "")
  
  z_out <- file.path(monthly_dir, "geopotential_500hPa_monthly.nc")
  if (!file.exists(z_out)) {
    cat("\nDownloading 500-hPa geopotential height...\n")
    nc_file <- download_z500_monthly(start_year, end_year)
    cat("  Downloaded:", nc_file, "\n")
    
    z_raw <- terra::rast(nc_file)
    z_dates <- validate_circulation_record(z_raw, "downloaded Z500")
    
    # ERA5 pressure-level geopotential is m2 s-2. Convert to geopotential
    # height in metres using the standard gravitational acceleration.
    z <- z_raw / G_GRAV_CONST
    terra::time(z) <- terra::time(z_raw)
    
    terra::writeCDF(
      z,
      z_out,
      overwrite = TRUE,
      compression = nc_compression
    )
    
    safe_remove(nc_file)
    rm(z_raw, z)
    gc(verbose = FALSE)
    cat("  Saved:", z_out, "\n")
  } else {
    cat("Skipping existing:", z_out, "\n")
    z_existing <- terra::rast(z_out)
    validate_circulation_record(z_existing, "existing Z500")
    rm(z_existing)
  }
  
  ivt_out <- file.path(monthly_dir, "ivt_components_monthly.nc")
  if (!file.exists(ivt_out)) {
    cat("\nDownloading vertically integrated moisture-flux components...\n")
    nc_file <- find_existing_ivt_download()
    if (is.na(nc_file)) {
      nc_file <- download_ivt_components_monthly(start_year, end_year)
      cat("  Downloaded:", nc_file, "\n")
    } else {
      cat("  Reusing existing interrupted IVT download:", nc_file, "\n")
    }
    
    ivt <- extract_ivt_components(nc_file)
    east_dates <- ivt$dates
    north_dates <- as.Date(terra::time(ivt$north), tz = "UTC")
    
    if (!identical(east_dates, north_dates)) {
      stop("Downloaded eastward and northward IVT time axes differ.")
    }
    
    validate_circulation_record(ivt$east, "downloaded eastward IVT")
    validate_circulation_record(ivt$north, "downloaded northward IVT")
    
    cat("  Writing terra-compatible two-variable IVT NetCDF...\n")
    write_ivt_components_netcdf(
      ivt$east,
      ivt$north,
      east_dates,
      ivt_out
    )
    
    safe_remove(nc_file)
    
    # Re-open and validate the project output exactly as 14circulation_context_v2.R
    # will consume it.
    ivt_check <- try(terra::sds(ivt_out), silent = TRUE)
    if (inherits(ivt_check, "try-error") || length(ivt_check) < 2L) {
      stop(
        "Written IVT NetCDF does not expose two subdatasets: ",
        ivt_out
      )
    }
    
    nm_check <- tolower(names(ivt_check))
    east_check_idx <- which(grepl("eastward|east", nm_check))
    north_check_idx <- which(grepl("northward|north", nm_check))
    
    if (!length(east_check_idx) || !length(north_check_idx)) {
      east_check_idx <- which(grepl("ivt_eastward|east", nm_check))
      north_check_idx <- which(grepl("ivt_northward|north", nm_check))
    }
    
    if (!length(east_check_idx) || !length(north_check_idx)) {
      stop(
        "Could not identify eastward/northward variables in written IVT NetCDF: ",
        ivt_out
      )
    }
    
    east_check <- ivt_check[[east_check_idx[1]]]
    north_check <- ivt_check[[north_check_idx[1]]]
    
    # terra should now recover the explicit CF time dimension. Assigning the
    # dates here is only a defensive check; the saved NetCDF itself contains
    # the time coordinate and requires no R-session metadata.
    east_check_dates <- read_netcdf_time_axis(
      ivt_out,
      expected_n = terra::nlyr(east_check)
    )
    north_check_dates <- read_netcdf_time_axis(
      ivt_out,
      expected_n = terra::nlyr(north_check)
    )
    
    validate_circulation_record(east_check, "retained eastward IVT")
    validate_circulation_record(north_check, "retained northward IVT")
    
    if (!identical(east_check_dates, north_check_dates)) {
      stop("Retained eastward and northward IVT time axes differ.")
    }
    
    rm(
      ivt, ivt_check, east_check, north_check,
      east_check_dates, north_check_dates
    )
    gc(verbose = FALSE)
    
    cat("  Saved:", ivt_out, "\n")
    cat("  IVT magnitude = sqrt(eastward_flux^2 + northward_flux^2); compute in 14circulation_context_v2.R\n")
  } else {
    cat("Skipping existing:", ivt_out, "\n")
    
    ivt_existing <- try(terra::sds(ivt_out), silent = TRUE)
    if (inherits(ivt_existing, "try-error") || length(ivt_existing) < 2L) {
      stop(
        "Existing IVT output does not expose two component subdatasets: ",
        ivt_out,
        ". Delete it and rerun Phase 3."
      )
    }
    
    nm <- tolower(names(ivt_existing))
    east_idx <- which(grepl("eastward|east", nm))
    north_idx <- which(grepl("northward|north", nm))
    
    if (!length(east_idx) || !length(north_idx)) {
      stop(
        "Existing IVT output does not expose identifiable eastward/northward ",
        "variables: ", ivt_out, ". Delete it and rerun Phase 3."
      )
    }
    
    east_existing <- ivt_existing[[east_idx[1]]]
    north_existing <- ivt_existing[[north_idx[1]]]
    
    # Validate the NetCDF coordinate directly rather than relying only on terra's
    # interpretation of the CDS metadata.
    existing_dates <- read_netcdf_time_axis(
      ivt_out,
      expected_n = terra::nlyr(east_existing)
    )
    
    if (terra::nlyr(north_existing) != length(existing_dates)) {
      stop("Existing northward IVT layer count differs from time axis.")
    }
    
    validate_circulation_record(east_existing, "existing eastward IVT")
    validate_circulation_record(north_existing, "existing northward IVT")
    
    rm(ivt_existing, east_existing, north_existing, existing_dates)
    gc(verbose = FALSE)
  }
  
  cat("\nPhase 3 completed and circulation inputs validated.\n")
  cat("  Z500: ", z_out, "\n", sep = "")
  cat("  IVT : ", ivt_out, "\n", sep = "")
  invisible(TRUE)
}

# ================================ MAIN =======================================
# The local GeoTIFF/hourly decision affects only the FAO-56 daily forcing branch.
# Phase 3 is independent and is explicitly enabled below.
process_monthly_state_variables()

if (has_complete_local_daily_tifs()) {
  process_local_daily_tifs()
} else {
  process_hourly_forcing()
}

validate_daily_archive("u2_mean")
validate_daily_archive("Rs")

# Required for Methods 3g / 14circulation_context_v2.R.
RUN_CIRCULATION_CONTEXT <- TRUE

if (RUN_CIRCULATION_CONTEXT) {
  process_circulation_context_fields()
}

cat("\n============================================================\n")
cat("ERA5-Land / ERA5 processing completed successfully.\n")
cat("Daily FAO-56 forcing:\n")
cat("  u2_mean : m s-1\n")
cat("  Rs      : MJ m-2 day-1\n")
cat("Circulation context:\n")
cat("  Z500    : 500-hPa geopotential height [m]\n")
cat("  IVT     : eastward/northward moisture-flux components [kg m-1 s-1]\n")
cat("============================================================\n")