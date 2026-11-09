# ==============================================================================
# Unified ERA5-Land / ERA5 Nechako download pipeline
#
# One DAILY forcing workflow ONLY: CDS ERA5-Land -> monthly_data_direct/*.nc.
# No GEE GeoTIFF conversion and no legacy month-by-month hourly alternative.
#
# PHASE 1: monthly ERA5-Land state/hydrological variables, including native PEV.
# PHASE 2: CDS ERA5-Land wind/solar, batched and resumable, on exact native grid.
# PHASE 3: optional ERA5 Z500 and IVT circulation context (existing implementation).
#
# Daily output location (point 2b_PET_ERALand.R here as well if needed):
#   monthly_data_direct/u2_mean_daily_YYYY.nc  (u2_mean, m s-1)
#   monthly_data_direct/Rs_daily_YYYY.nc       (Rs, MJ m-2 day-1)
#
# u2 = daily mean of sqrt(u10^2 + v10^2) * 4.87/log(67.8*10 - 5.42).
# Rs(D) = ERA5-Land ssrd at 00 UTC of D+1 / 1e6; first 1950-01-01 00
# stamp is not used, as it would describe 1949-12-31.
#
# Existing annual NetCDFs are never silently labelled as CDS outputs. The
# optional reuse mode preserves legacy files only after date/grid/variable checks;
# their actual source and wind calculation cannot be established by those checks.
# To rebuild all annual files, back up the archive and set OVERWRITE_DAILY <- TRUE.
# A replacement is validated before the previous file is moved.
#
# Run from the project root (folder containing .env).
# ==============================================================================

rm(list = ls())
WD_PATH <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(WD_PATH, ".env"))) {
  stop(".env not found in working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running.")
}

library(ecmwfr)
library(terra)
library(lubridate)
library(ncdf4)
library(dotenv)
if (!requireNamespace("curl", quietly = TRUE))
  stop("Package 'curl' is required for CDS connectivity checking.")

dotenv::load_dot_env(file.path(WD_PATH, ".env"))
cds_user <- Sys.getenv("CDS_USER")
cds_key  <- Sys.getenv("CDS_KEY")
if (!nzchar(cds_user) || !nzchar(cds_key))
  stop("CDS_USER and/or CDS_KEY is missing from .env")
wf_set_key(user = cds_user, key = cds_key)

# ============================= CONFIGURATION =================================
RUN_MONTHLY     <- TRUE
RUN_DAILY       <- TRUE
RUN_CIRCULATION <- FALSE   # TRUE to also run the preserved ERA5 Z500/IVT phase
REGENERATE_PEV  <- FALSE   # TRUE explicitly refreshes existing, verified PEV

# Both daily forcing and monthly products are stored in monthly_data_direct/.
# Keep existing annual files unchanged by default; generate only missing years.
# Legacy annual files may have no CDS provenance attributes. Setting TRUE allows
# their reuse ONLY if variable name, complete daily dates and grid are validated.
# A reused legacy file is *not* certified as CDS-derived or methodologically
# identical to new outputs. Set FALSE to enforce strict CDS-only provenance.
REUSE_EXISTING_UNVERIFIED_DAILY <- TRUE
# Set TRUE only after backing up your existing archive: rebuilds all years.
OVERWRITE_DAILY <- FALSE
CLEAN_RAW_AFTER_SUCCESS <- TRUE

bbox <- c(-127.8, 52.9, -122.7, 56.2)  # xmin, ymin, xmax, ymax
start_date <- as.Date("1950-01-01")
end_date   <- as.Date("2025-12-31")
YEARS <- seq.int(as.integer(format(start_date, "%Y")),
                 as.integer(format(end_date, "%Y")))
if (format(start_date, "%m-%d") != "01-01" ||
    format(end_date, "%m-%d") != "12-31")
  stop("Daily CDS annual workflow requires complete calendar years.")

monthly_variables <- c(
  "snow_depth", "snow_density", "snowfall", "skin_reservoir_content",
  "volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2",
  "volumetric_soil_water_layer_3", "volumetric_soil_water_layer_4",
  "runoff", "sub_surface_runoff", "surface_runoff",
  "total_evaporation", "potential_evaporation"
)
runoff_variables <- c("runoff", "sub_surface_runoff", "surface_runoff")
corrected_accumulation_variables <- c(
  "total_evaporation", "potential_evaporation", "snowfall", runoff_variables
)

# Wind is a sampled scalar mean, NOT the magnitude of the mean wind vector.
WORKERS <- 8L                    # lower to 4 if CDS rejects batched submissions
WIND_HOUR_STEP <- 3L             # 1 = exact hourly mean, 3/6 = sampled daily mean
WIND_MONTHS_CANDIDATES <- c(12L, 6L, 3L, 1L)
SSRD_YEARS_PER_REQUEST <- 10L
MAX_ATTEMPTS_PER_GROUP <- 5L
MAX_WAIT_MIN <- 180
WIND_METHOD <- "mean_hourly_speed"
stopifnot(WIND_METHOD %in% c("mean_hourly_speed", "vector_mean"),
          WIND_HOUR_STEP %in% c(1L, 2L, 3L, 4L, 6L, 8L, 12L, 24L),
          24L %% WIND_HOUR_STEP == 0L,
          all(12L %% WIND_MONTHS_CANDIDATES == 0L), WORKERS >= 1L)
WIND_HOURS <- seq.int(0L, 23L, by = WIND_HOUR_STEP)
wind_10m_to_2m_factor <- 4.87 / log(67.8 * 10 - 5.42)
nc_compression <- 9L

# Optional quality comparison is read-only; this is NOT another forcing source.
RUN_COMPARISON <- FALSE
COMPARISON_YEAR <- 2020L
REFERENCE_DIR <- file.path(WD_PATH, "daily_eto_forcing_gee_reference")

monthly_dir <- file.path(WD_PATH, "monthly_data_direct")
daily_dir <- monthly_dir  # IMPORTANT: annual Rs/u2 outputs live alongside monthly NetCDFs
download_dir <- file.path(WD_PATH, "downloads")
RAW_DIR <- file.path(download_dir, "raw_cds_daily")
for (d in c(monthly_dir, daily_dir, download_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# =============================== SHARED HELPERS ==============================
month_sequence <- function(start_date, end_date) {
  seq(lubridate::floor_date(start_date, "month"),
      lubridate::floor_date(end_date, "month"), by = "month")
}

safe_remove <- function(path) {
  if (length(path) == 1L && !is.na(path) && file.exists(path)) {
    if (!file.remove(path)) warning("Could not delete file: ", path)
  }
  invisible(TRUE)
}

check_time_axis <- function(x, label) {
  tt <- terra::time(x)
  if (is.null(tt) || length(tt) != terra::nlyr(x) || anyNA(tt))
    stop("Missing/invalid NetCDF time axis for ", label)
  as.POSIXct(tt, tz = "UTC")
}

report_na <- function(x, label) {
  cat("\n  --- NA Check: ", label, " ---\n", sep = "")
  na_by_layer <- terra::global(is.na(x), "sum", na.rm = TRUE)
  total_na <- sum(na_by_layer[, 1], na.rm = TRUE)
  total_cells <- as.double(terra::ncell(x)) * terra::nlyr(x)
  cat(sprintf("  Total cells: %.0f\n  Total NAs: %.0f\n  NA percent: %.6f%%\n",
              total_cells, total_na, 100 * total_na / total_cells))
  if (total_na > 0) warning(label, " contains missing values.")
  invisible(total_na)
}

read_nc_variable <- function(nc_file, candidates) {
  candidates <- unique(tolower(candidates))
  s <- try(terra::sds(nc_file), silent = TRUE)
  if (!inherits(s, "try-error") && length(s) > 0L) {
    s_names <- tolower(names(s))
    idx <- which(vapply(s_names, function(z)
      any(vapply(candidates, function(cn) grepl(cn, z, fixed = TRUE), logical(1))),
      logical(1)))
    if (length(idx)) return(s[[idx[1]]])
    if (length(s) == 1L) return(s[[1]])
  }
  r <- terra::rast(nc_file)
  nms <- tolower(names(r))
  idx <- which(vapply(nms, function(z)
    any(vapply(candidates, function(cn) grepl(cn, z, fixed = TRUE), logical(1))),
    logical(1)))
  if (length(idx)) return(r[[idx]])
  stop("Could not identify variable ", paste(candidates, collapse = "/"),
       " in ", nc_file)
}

spatraster_to_nc_array <- function(x) {
  nr <- terra::nrow(x); nc <- terra::ncol(x); nl <- terra::nlyr(x)
  if (nr < 1L || nc < 1L || nl < 1L) stop("Cannot write empty SpatRaster.")
  m <- terra::values(x, mat = TRUE)
  if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
  if (nrow(m) != nr * nc || ncol(m) != nl) stop("Unexpected raster matrix shape")
  a <- array(NA_real_, dim = c(nc, nr, nl))
  for (k in seq_len(nl))
    a[, , k] <- t(matrix(m[, k], nrow = nr, ncol = nc, byrow = TRUE))
  a
}

nc_time_origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")

# Annual writer: staged NetCDF -> validate -> install. Never remove an old output
# until the replacement is complete. A locked old file is left untouched.
write_daily_netcdf_from_spatraster <- function(x, dates, out_file, variable_name,
                                               long_name, units, overwrite = FALSE) {
  dates <- as.Date(dates)
  if (terra::nlyr(x) != length(dates) || anyNA(dates) ||
      anyDuplicated(dates) || is.unsorted(dates))
    stop("Invalid date/layer mapping for ", variable_name)
  if (!overwrite && file.exists(out_file))
    stop("Output exists and overwrite=FALSE: ", out_file)
  if (!terra::is.lonlat(x)) stop("Expected longitude/latitude NetCDF grid: ", out_file)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  stage <- paste0(out_file, ".staging.nc")
  backup <- paste0(out_file, ".replacement-backup.nc")
  if (file.exists(stage) || file.exists(backup))
    stop("Staging/backup file already exists; inspect before retry: ", out_file)
  cleanup_stage <- TRUE
  nc_is_open <- FALSE
  on.exit({
    if (nc_is_open) try(ncdf4::nc_close(nc), silent = TRUE)
    if (cleanup_stage && file.exists(stage)) safe_remove(stage)
  }, add = TRUE)
  
  nr <- terra::nrow(x); ncol_x <- terra::ncol(x)
  e <- as.vector(terra::ext(x)) # xmin, xmax, ymin, ymax
  xres <- (e[2] - e[1]) / ncol_x
  yres <- (e[4] - e[3]) / nr
  lon <- e[1] + (seq_len(ncol_x) - 0.5) * xres
  lat <- e[4] - (seq_len(nr) - 0.5) * yres
  time_vals <- as.numeric(difftime(as.POSIXct(dates, tz = "UTC"),
                                   nc_time_origin, units = "days"))
  arr <- spatraster_to_nc_array(x)
  fill_value <- -9999.0
  arr[!is.finite(arr)] <- fill_value
  lon_dim <- ncdf4::ncdim_def("lon", "degrees_east", lon)
  lat_dim <- ncdf4::ncdim_def("lat", "degrees_north", lat)
  time_dim <- ncdf4::ncdim_def("time", "days since 1970-01-01 00:00:00",
                               time_vals, unlim = FALSE)
  var_def <- ncdf4::ncvar_def(variable_name, units,
                              list(lon_dim, lat_dim, time_dim),
                              missval = fill_value, longname = long_name,
                              prec = "float", compression = nc_compression)
  nc <- ncdf4::nc_create(stage, vars = list(var_def), force_v4 = TRUE)
  nc_is_open <- TRUE
  ncdf4::ncvar_put(nc, var_def, arr)
  ncdf4::ncatt_put(nc, 0, "Conventions", "CF-1.8")
  ncdf4::ncatt_put(nc, 0, "title", long_name)
  ncdf4::ncatt_put(nc, 0, "source", "ERA5-Land CDS unified daily forcing pipeline")
  ncdf4::ncatt_put(nc, 0, "source_dataset", "reanalysis-era5-land")
  ncdf4::ncatt_put(nc, 0, "source_workflow", "cds_daily_v5_unified")
  ncdf4::ncatt_put(nc, 0, "history", paste("Created", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
  ncdf4::ncatt_put(nc, "lon", "axis", "X")
  ncdf4::ncatt_put(nc, "lat", "axis", "Y")
  ncdf4::ncatt_put(nc, "time", "axis", "T")
  ncdf4::ncatt_put(nc, "time", "calendar", "standard")
  if (identical(variable_name, "u2_mean")) {
    ncdf4::ncatt_put(nc, 0, "wind_method", WIND_METHOD)
    ncdf4::ncatt_put(nc, 0, "wind_hour_step", as.integer(WIND_HOUR_STEP))
  } else if (identical(variable_name, "Rs")) {
    ncdf4::ncatt_put(nc, 0, "radiation_daily_alignment", "D+1 00:00 UTC accumulation for day D")
  }
  ncdf4::nc_close(nc)
  nc_is_open <- FALSE
  
  check <- ncdf4::nc_open(stage)
  stage_valid <- tryCatch({
    identical(as.integer(check$dim$time$len), as.integer(length(dates))) &&
      identical(as.numeric(check$dim$time$vals), time_vals) &&
      variable_name %in% names(check$var) &&
      identical(as.integer(check$dim$lon$len), as.integer(ncol_x)) &&
      identical(as.integer(check$dim$lat$len), as.integer(nr))
  }, error = function(e) FALSE)
  ncdf4::nc_close(check)
  if (!stage_valid) stop("Staged NetCDF failed time/grid/variable validation: ", stage)
  cleanup_stage <- FALSE  # validated replacement should survive a Windows lock
  if (file.exists(out_file)) {
    if (!overwrite) stop("Output exists: ", out_file)
    if (!file.rename(out_file, backup))
      stop("Cannot replace locked output: ", out_file,
           "\nValidated replacement remains at ", stage)
    if (!file.rename(stage, out_file)) {
      if (!file.rename(backup, out_file))
        stop("Failed installation AND rollback. Restore backup at ", backup)
      stop("Replacement install failed; old output restored. New file at ", stage)
    }
    if (!file.remove(backup)) warning("Old-file backup retained: ", backup)
  } else if (!file.rename(stage, out_file)) {
    stop("Failed to install validated NetCDF: ", stage)
  }
  invisible(out_file)
}

# ========================== MONTHLY STATE DOWNLOAD ============================

download_monthly_variable <- function(variable, start_year, end_year) {
  # ECMWF documents an accumulated-field error in the standard
  # monthly_averaged_reanalysis product for September 2022-February 2024.
  # For potential evaporation (PEV), use the corrected synoptic monthly mean
  # at 00 UTC throughout the record.
  #
  # IMPORTANT: Do NOT add a `step` field to this CDS request. The
  # monthly_averaged_reanalysis_by_hour_of_day catalogue request is specified
  # by product_type + time; the 00 UTC product carries the relevant accumulated
  # hydrological field. Script 2b documents the resulting 24-hour accumulation
  # interpretation and performs the m -> mm conversion.
  # IMPORTANT ECMWF data-quality handling:
  # Accumulated hydrological variables in the standard
  # `monthly_averaged_reanalysis` product were affected by a CDS/ERA5-Land
  # data-copy problem for September 2022 through February 2024. ECMWF
  # identifies `monthly_averaged_reanalysis_by_hour_of_day` as the recovery
  # product. Total evaporation (`e`) is an accumulated hydrological variable,
  # so it MUST use the corrected by-hour-of-day product, just as PEV does.
  #
  # We deliberately use the corrected product for the ENTIRE 1950-2025 record
  # rather than mixing products across time. This keeps the AET source
  # homogeneous and avoids introducing a product-boundary artefact.
  # (The list of affected variables, `corrected_accumulation_variables`, is
  # defined once in the CONFIGURATION section and now also covers runoff.)
  
  product_type <- if (variable %in% corrected_accumulation_variables) {
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

# ------------------------------------------------------------------------------
# Safely install a regenerated monthly product. The old target is never removed
# before the new product has been downloaded and successfully written. This is
# critical on Windows, where QGIS/Panoply/ArcGIS may keep a NetCDF file locked.
# ------------------------------------------------------------------------------
install_monthly_product <- function(result, output_file, variable) {
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  tmp_output <- file.path(
    dirname(output_file),
    paste0(".", tools::file_path_sans_ext(basename(output_file)),
           "_new.nc")
  )
  
  safe_remove(tmp_output)
  
  cat("  Writing replacement NetCDF:", tmp_output, "\n")
  writeCDF(
    result,
    tmp_output,
    varname = variable,
    longname = variable,
    overwrite = TRUE,
    compression = nc_compression
  )
  
  if (!file.exists(tmp_output)) {
    stop("Replacement NetCDF was not created: ", tmp_output)
  }
  
  # Validate the replacement before touching the existing target.
  check <- try(rast(tmp_output), silent = TRUE)
  if (inherits(check, "try-error")) {
    safe_remove(tmp_output)
    stop("Replacement NetCDF cannot be reopened: ", tmp_output)
  }
  
  if (nlyr(check) != nlyr(result)) {
    n_bad <- nlyr(check)
    n_exp <- nlyr(result)
    rm(check)
    safe_remove(tmp_output)
    stop(
      "Replacement NetCDF layer count mismatch for ", variable,
      ": expected ", n_exp, ", got ", n_bad
    )
  }
  
  replacement_time <- try(check_time_axis(check, tmp_output), silent = TRUE)
  expected_time <- try(check_time_axis(result, variable), silent = TRUE)
  if (inherits(replacement_time, "try-error") ||
      inherits(expected_time, "try-error") ||
      length(replacement_time) != length(expected_time) ||
      !isTRUE(all(as.numeric(as.Date(replacement_time)) ==
                  as.numeric(as.Date(expected_time))))) {
    rm(check)
    safe_remove(tmp_output)
    stop("Replacement NetCDF time axis validation failed for ", variable)
  }
  rm(check, replacement_time, expected_time)
  gc(verbose = FALSE)
  
  # Install safely: NEVER delete the old file before the new one is in place.
  # On any installation failure, restore the old file from its backup.
  backup <- paste0(output_file, ".replacement-backup.nc")
  if (file.exists(output_file)) {
    if (file.exists(backup))
      stop("Monthly replacement backup already exists; inspect before retry: ", backup)
    if (!file.rename(output_file, backup))
      stop("Cannot move existing monthly NetCDF to a safety backup (possibly locked): ",
           output_file, "\nValidated replacement remains at: ", tmp_output)
    if (!file.rename(tmp_output, output_file)) {
      if (!file.rename(backup, output_file))
        stop("Monthly installation and rollback both failed. Restore: ", backup,
             "\nValidated replacement remains at: ", tmp_output)
      stop("Monthly installation failed; old file restored: ", output_file,
           "\nValidated replacement remains at: ", tmp_output)
    }
    if (!file.remove(backup)) warning("Old monthly backup retained: ", backup)
  } else if (!file.rename(tmp_output, output_file)) {
    stop("Failed to install monthly NetCDF: ", output_file,
         "\nValidated output remains at: ", tmp_output)
  }
  
  cat("  Saved:", output_file, "\n")
  invisible(output_file)
}
resolve_nc_data_var <- function(nc, variable) {
  short_names <- c(
    snowfall = "sf", potential_evaporation = "pev", total_evaporation = "e",
    runoff = "ro", surface_runoff = "sro", sub_surface_runoff = "ssro"
  )
  coord_like <- c("time", "valid_time", "latitude", "longitude", "lat", "lon",
                  "crs", "expver", "number", "spatial_ref", "time_bnds")
  vn <- names(nc$var)
  data_vars <- vn[!tolower(vn) %in% coord_like]
  
  hit <- data_vars[tolower(data_vars) %in% tolower(c(variable, short_names[variable]))]
  if (length(hit)) return(hit[1])
  if (length(data_vars) == 1L) return(data_vars[1])   # single data variable: unambiguous
  stop("Could not identify '", variable, "'. Available variables: ",
       paste(vn, collapse = ", "))
}
add_accumulation_source_metadata <- function(output_file, variable) {
  nc <- ncdf4::nc_open(output_file, write = TRUE)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  target_name <- resolve_nc_data_var(nc, variable)
  ncdf4::ncatt_put(
    nc, 0, "source_product",
    "monthly_averaged_reanalysis_by_hour_of_day"
  )
  ncdf4::ncatt_put(nc, 0, "source_time", "00:00 UTC")
  ncdf4::ncatt_put(
    nc, 0, "accumulation_interpretation",
    "24-hour accumulation ending at 00 UTC; monthly mean of synoptic accumulations"
  )
  ncdf4::ncatt_put(
    nc, 0, "data_quality_handling",
    "By-hour-of-day product used to avoid the documented standard monthly accumulated-variable data-copy issue for September 2022-February 2024"
  )
  ncdf4::ncatt_put(
    nc, 0, "native_units_interpretation",
    if (variable %in% c("potential_evaporation", "total_evaporation"))
      "signed metres of water equivalent per 24-hour accumulation"
    else "non-negative metres of water equivalent per 24-hour accumulation"
  )
  
  if (identical(variable, "potential_evaporation")) {
    ncdf4::ncatt_put(
      nc, 0, "downstream_conversion",
      "PEV_mm_day = -PEV_native_m * 1000; PEV is potential/open-water evaporation and is not AET"
    )
    ncdf4::ncatt_put(
      nc, target_name, "long_name",
      "ERA5-Land native potential evaporation; 00 UTC 24-hour accumulation"
    )
  } else if (identical(variable, "total_evaporation")) {
    ncdf4::ncatt_put(
      nc, 0, "downstream_use",
      "Source for actual evapotranspiration (AET) preprocessing; do not substitute PEV"
    )
    ncdf4::ncatt_put(
      nc, 0, "downstream_conversion",
      "AET_mm_day = -total_evaporation_native_m * 1000; verify returned CDS variable metadata before monthly aggregation"
    )
    ncdf4::ncatt_put(
      nc, target_name, "long_name",
      "ERA5-Land total evaporation; 00 UTC 24-hour accumulation; source for AET"
    )
  }
  
  if (variable %in% runoff_variables) {
    ncdf4::ncatt_put(
      nc, 0, "downstream_conversion",
      paste0(variable, "_mm_day = ", variable, "_native_m * 1000; runoff is a non-negative accumulated water depth")
    )
    ncdf4::ncatt_put(
      nc, target_name, "long_name",
      paste0("ERA5-Land ", gsub("_", " ", variable), "; 00 UTC 24-hour accumulation")
    )
  }
  
  ncdf4::ncatt_put(
    nc, target_name, "units",
    "m per 24-hour accumulation"
  )
  
  invisible(TRUE)
}

# TRUE if the retained NetCDF already carries the corrected-product provenance
# attributes written by add_accumulation_source_metadata().
has_provenance <- function(f) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  attr <- ncdf4::ncatt_get(nc, 0, "source_product")
  isTRUE(attr$hasatt) &&
    identical(as.character(attr$value), "monthly_averaged_reanalysis_by_hour_of_day")
}

process_monthly_state_variables <- function() {
  start_year <- year(start_date)
  end_year <- year(end_date)
  
  cat("\n============================================================\n")
  cat("PHASE 1: ERA5-Land monthly state variables\n")
  cat("============================================================\n")
  
  for (var in monthly_variables) {
    output_file <- file.path(monthly_dir, paste0(var, "_monthly.nc"))
    force_regenerate <- (identical(var, "potential_evaporation") && REGENERATE_PEV) ||
      (file.exists(output_file) && var %in% corrected_accumulation_variables &&
         !has_provenance(output_file))
    
    if (file.exists(output_file) && !force_regenerate) {
      cat("Skipping verified/existing monthly file:", output_file, "\n")
      next
    }
    if (force_regenerate && file.exists(output_file) &&
        var %in% corrected_accumulation_variables && !has_provenance(output_file)) {
      cat("Unverified monthly accumulation will be regenerated from CDS: ",
          output_file, "\n", sep = "")
    }
    if (force_regenerate && file.exists(output_file) &&
        identical(var, "potential_evaporation")) {
      cat(
        "\nRegenerating existing PEV file using corrected synoptic product: ",
        output_file, "\n", sep = ""
      )
      cat(
        "  Existing file is NOT deleted before the replacement has been fully ",
        "downloaded and validated.\n"
      )
    }
    
    cat("\nProcessing monthly variable:", var, "\n")
    nc_file <- download_monthly_variable(var, start_year, end_year)
    
    result <- rast(nc_file)
    report_na(result, var)
    
    # Monthly ERA5-Land accumulated evaporation variables are retained in their
    # native signed water-depth representation. The sign/unit conversion belongs
    # to the downstream preprocessing script, not this downloader.
    #
    # For both PEV and total evaporation, the by-hour-of-day product at 00 UTC
    # provides the synoptic accumulation required to avoid the documented
    # standard-monthly-product data-copy issue.
    #
    # IMPORTANT: total_evaporation is the source for downstream AET; it is NOT
    # reconstructed from PEV and PEV is never used as an AET surrogate.
    if (var == "snow_depth") {
      cat("  Snow depth retained in metres.\n")
    } else if (var == "snow_density") {
      cat("  Snow density retained in kg m-3.\n")
    } else if (var == "snowfall") {
      cat("  Snowfall source: monthly_averaged_reanalysis_by_hour_of_day (00:00 UTC).\n")
      cat("  Snowfall retained in metres of water equivalent per 24-hour accumulation.\n")
    } else if (var == "snow_depth_water_equivalent") {
      cat("  SWE retained in metres water equivalent.\n")
    } else if (grepl("volumetric_soil_water", var, fixed = TRUE)) {
      cat("  Soil moisture retained in m3 m-3.\n")
    } else if (identical(var, "skin_reservoir_content")) {
      cat("  Skin reservoir content retained in m of water equivalent.\n")
    } else if (var %in% runoff_variables) {
      cat("  Runoff source: monthly_averaged_reanalysis_by_hour_of_day.\n")
      cat("  Runoff time: 00:00 UTC.\n")
      cat("  Runoff interpretation: 24-hour accumulation ending at 00 UTC.\n")
      cat("  Runoff native units retained: m per 24-hour accumulation.\n")
      cat("  Downstream conversion: mm/day = value * 1000.\n")
    } else if (var == "surface_pressure") {
      cat("  Surface pressure retained in Pa.\n")
    } else if (grepl("temperature", var, fixed = TRUE)) {
      cat("  Temperature retained in K.\n")
    } else if (identical(var, "potential_evaporation")) {
      cat("  PEV source: monthly_averaged_reanalysis_by_hour_of_day.\n")
      cat("  PEV time: 00:00 UTC.\n")
      cat("  PEV interpretation: 24-hour accumulation ending at 00 UTC.\n")
      cat("  PEV native units retained: m per 24-hour accumulation.\n")
      cat("  2b conversion: signed native metres -> positive mm/day by -value * 1000.\n")
    } else if (identical(var, "total_evaporation")) {
      cat("  Total evaporation source: monthly_averaged_reanalysis_by_hour_of_day.\n")
      cat("  Total evaporation time: 00:00 UTC.\n")
      cat("  Total evaporation interpretation: 24-hour accumulation ending at 00 UTC.\n")
      cat("  Native signed water-depth units are retained for downstream AET preprocessing.\n")
      cat("  Downstream AET conversion must reverse the ERA5 evaporation sign and convert m -> mm.\n")
      cat("  PEV is NOT used as an AET surrogate.\n")
    }
    
    install_monthly_product(
      result,
      output_file,
      var
    )
    
    if (var %in% corrected_accumulation_variables) {
      add_accumulation_source_metadata(output_file, var)
      cat("  Corrected accumulated-variable provenance metadata written to retained NetCDF.\n")
    }
    
    safe_remove(nc_file)
    rm(result)
    gc(verbose = FALSE)
  }
}


# ==================== PHASE 2: ONE CDS DAILY FORCING PATH =====================
# The CDS workflow below supersedes both the older hourly-month loop and the
# GEE GeoTIFF conversion path; it writes ONLY daily_dir/ annual NetCDF files.
# ---- REQUEST BUILDERS -------------------------------------------------------
mk_request <- function(vars, years, months, days, hours, target) {
  list(
    dataset_short_name = "reanalysis-era5-land",
    variable           = vars,
    year               = as.character(years),
    month              = sprintf("%02d", months),
    day                = sprintf("%02d", days),
    time               = sprintf("%02d:00", hours),
    data_format        = "netcdf",
    download_format    = "unarchived",
    area               = AREA,
    target             = target
  )
}

make_chunks <- function(years, n) unname(split(years, ceiling(seq_along(years) / n)))
ssrd_chunks <- make_chunks(YEARS, SSRD_YEARS_PER_REQUEST)

wind_blocks <- function(k) {
  starts <- seq(1, 12, by = k)
  unlist(lapply(YEARS, function(y) lapply(starts, function(m0)
    list(year = y, months = m0:(m0 + k - 1)))), recursive = FALSE)
}
wind_file <- function(b) file.path(RAW_DIR, sprintf("u10v10_h%02d_%04d_%02d-%02d.nc",
                                                    WIND_HOUR_STEP, b$year, min(b$months), max(b$months)))
ssrd_file <- function(ch) file.path(RAW_DIR, sprintf("ssrd00_%04d_%04d.nc", min(ch), max(ch)))
SSRD_FINAL_FILE <- file.path(RAW_DIR, sprintf("ssrd00_%04d_jan01.nc", max(YEARS) + 1L))
SSRD_VAR <- "surface_solar_radiation_downwards"

ssrd_requests <- function() c(
  lapply(ssrd_chunks, function(ch)
    mk_request(SSRD_VAR, ch, 1:12, 1:31, 0, basename(ssrd_file(ch)))),
  list(mk_request(SSRD_VAR, max(YEARS) + 1L, 1, 1, 0, basename(SSRD_FINAL_FILE)))
)
wind_requests <- function(k) lapply(wind_blocks(k), function(b)
  mk_request(c("10m_u_component_of_wind", "10m_v_component_of_wind"),
             b$year, b$months, 1:31, WIND_HOURS, basename(wind_file(b))))

WIND_K <- NA_integer_  # discovered once in process_daily_cds_forcing()

# ---- TIME AXIS FROM THE REQUEST ---------------------------------------------
SMOKE_REQ <- NULL   # initialized after TEMPLATE and AREA are loaded
ALL_REQS  <- list()
register_requests <- function(reqs) {
  names(reqs) <- vapply(reqs, `[[`, "", "target")
  ALL_REQS <<- c(ALL_REQS[setdiff(names(ALL_REQS), names(reqs))], reqs)
}
# Requests are registered only inside process_daily_cds_forcing(), after AREA
# is initialized from the actual project template.

# ERA5-Land starts on 1950-01-01 and CDS omits some stamps on that first day
# (00:00 for accumulated ssrd, and also for the 10-m wind). Missing stamps are
# tolerated ONLY before this instant; a gap anywhere else still stops the script.
ERA5L_TOLERANCE_END <- as.POSIXct("1950-01-02 00:00:00", tz = "UTC")

expected_stamps <- function(req) {
  yrs <- as.integer(req$year); mos <- as.integer(req$month); dys <- as.integer(req$day)
  hrs <- as.integer(substr(req$time, 1, 2))
  d <- seq(as.Date(sprintf("%04d-01-01", min(yrs))), as.Date(sprintf("%04d-12-31", max(yrs))), by = "day")
  d <- d[as.integer(format(d, "%Y")) %in% yrs & as.integer(format(d, "%m")) %in% mos &
           as.integer(format(d, "%d")) %in% dys]
  as.POSIXct(paste(rep(format(d, "%Y-%m-%d"), each = length(hrs)),
                   sprintf("%02d:00:00", rep(hrs, times = length(d)))), tz = "UTC")
}

# Read the time coordinate straight from the NetCDF (CDS uses "valid_time"; older files "time").
# Independent of terra, which returns NA times for these subdatasets.
nc_time_axis <- function(f) {
  nc <- ncdf4::nc_open(f)
  on.exit(ncdf4::nc_close(nc))
  nm <- intersect(c("valid_time", "time"), c(names(nc$dim), names(nc$var)))
  if (!length(nm)) return(NULL)
  nm <- nm[1]
  if (nm %in% names(nc$dim) && length(nc$dim[[nm]]$vals) > 0L) {
    vals  <- as.numeric(nc$dim[[nm]]$vals)
    units <- nc$dim[[nm]]$units
  } else {
    vals  <- as.numeric(ncdf4::ncvar_get(nc, nm))
    units <- ncdf4::ncatt_get(nc, nm, "units")$value
  }
  m <- regmatches(units, regexec("^\\s*(seconds|minutes|hours|days)\\s+since\\s+(.+)$",
                                 units, ignore.case = TRUE))[[1]]
  if (length(m) != 3L) return(NULL)
  mult   <- c(seconds = 1, minutes = 60, hours = 3600, days = 86400)[[tolower(m[2])]]
  origin <- as.POSIXct(substr(trimws(m[3]), 1, 19), tz = "UTC",
                       tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"))
  as.POSIXct(round(as.numeric(origin) + vals * mult), origin = "1970-01-01", tz = "UTC")
}

get_axis <- function(r, f) {
  req <- ALL_REQS[[basename(f)]]
  if (is.null(req)) stop("No registered CDS request for ", basename(f))
  expected <- expected_stamps(req)
  observed <- tryCatch(nc_time_axis(f), error = function(e) NULL)
  if (is.null(observed) || length(observed) != terra::nlyr(r) || anyNA(observed))
    stop(basename(f), ": missing/invalid NetCDF time axis.")
  if (is.unsorted(observed) || anyDuplicated(observed))
    stop(basename(f), ": timestamps are duplicated or not sorted.")
  exp_s <- round(as.numeric(expected))
  obs_s <- round(as.numeric(observed))
  missing <- expected[!(exp_s %in% obs_s)]
  extra <- observed[!(obs_s %in% exp_s)]
  if (length(extra) || any(missing >= ERA5L_TOLERANCE_END))
    stop(basename(f), ": time axis differs from request. Missing: ",
         paste(format(head(missing, 5), "%Y-%m-%d %H:%M"), collapse = ", "),
         " | Unexpected: ", paste(format(head(extra, 5)), collapse = ", "))
  if (length(missing))
    message(basename(f), ": tolerated ", length(missing),
            " unavailable ERA5-Land stamp(s) on 1950-01-01.")
  observed
}

# ---- NETWORK HELPERS --------------------------------------------------------
cds_reachable <- function() {
  ip <- tryCatch(curl::nslookup("cds.climate.copernicus.eu", error = FALSE),
                 error = function(e) NULL)
  !is.null(ip)
}
wait_for_cds <- function(max_wait_min = MAX_WAIT_MIN) {
  t0 <- Sys.time()
  while (!cds_reachable()) {
    if (as.numeric(difftime(Sys.time(), t0, units = "mins")) > max_wait_min)
      stop("No connection to CDS for ", max_wait_min,
           " min. Fix the network and re-run (finished files are kept).", call. = FALSE)
    message(format(Sys.time(), "%H:%M:%S"), "  CDS unreachable; retrying in 60 s...")
    Sys.sleep(60)
  }
  invisible(TRUE)
}

# ---- 0. SMOKE TEST ----------------------------------------------------------
smoke_test <- function() {
  f <- file.path(RAW_DIR, "_smoke_ssrd_2020-01.nc")
  if (!raw_request_valid(SMOKE_REQ)) {
    wait_for_cds()
    message("Smoke test: one small CDS request (31 fields)...")
    res <- tryCatch(
      wf_request(user = cds_user, request = mk_request(SSRD_VAR, 2020, 1, 1:31, 0, basename(f)),
                 transfer = TRUE, path = RAW_DIR, time_out = 1800),
      error = function(e) e)
    if (inherits(res, "error"))
      stop("CDS smoke test failed: ", conditionMessage(res),
           if (grepl("licen|terms|403", conditionMessage(res), ignore.case = TRUE))
             "\nAccept the ERA5-Land licence on the dataset page of the CDS website, then re-run." else "",
           call. = FALSE)
    got <- as.character(unlist(res))[1]
    if (!file.exists(f) && !is.na(got) && file.exists(got)) file.rename(got, f)
    if (!file.exists(f)) stop("CDS smoke test: no file returned.", call. = FALSE)
  }
  r  <- read_nc_variable(f, c(SSRD_VAR, "ssrd"))
  invisible(get_axis(r, f))
  if (nlyr(r) != 31L) stop("Smoke test: expected 31 layers, got ", nlyr(r), call. = FALSE)
  if (!isTRUE(terra::compareGeom(r[[1]], TEMPLATE, stopOnError = FALSE, crs = FALSE,
                                 ext = TRUE, rowcol = TRUE, res = TRUE)))
    stop("Smoke test: CDS subset grid differs from the project template (",
         nrow(r), " x ", ncol(r), " vs ", nrow(TEMPLATE), " x ", ncol(TEMPLATE), ").", call. = FALSE)
  message("Smoke test passed: credentials, licence, variable name, time axis and grid are OK.")
}
# ---- 1. DOWNLOAD (with connectivity wait + retry) ----------------------------
# Check completed raw files by their actual CDS time coordinates. This prevents
# a partial/interrupted NetCDF from masquerading as a successful download.
raw_request_valid <- function(req) {
  f <- file.path(RAW_DIR, req$target)
  if (!file.exists(f)) return(FALSE)
  tt <- tryCatch(nc_time_axis(f), error = function(e) NULL)
  expected <- expected_stamps(req)
  good <- !is.null(tt) && length(tt) > 0L && !anyNA(tt) &&
    !is.unsorted(tt) && !anyDuplicated(tt) &&
    all(round(as.numeric(tt)) %in% round(as.numeric(expected))) &&
    !any(expected[!(round(as.numeric(expected)) %in% round(as.numeric(tt)))] >=
           ERA5L_TOLERANCE_END)
  if (!good) {
    warning("Partial/invalid raw CDS NetCDF will be retried: ", f)
    if (!file.remove(f)) stop("Could not remove invalid raw NetCDF: ", f)
    return(FALSE)
  }
  TRUE
}

run_requests <- function(reqs) {
  pending <- function(g) Filter(function(r) !raw_request_valid(r), g)
  todo <- pending(reqs)
  if (!length(todo)) return(invisible(character(0)))
  groups <- split(todo, ceiling(seq_along(todo) / WORKERS))
  message(sprintf("Submitting %d CDS request(s) in %d group(s) of up to %d...",
                  length(todo), length(groups), WORKERS))
  for (gi in seq_along(groups)) {
    g <- groups[[gi]]
    for (attempt in seq_len(MAX_ATTEMPTS_PER_GROUP)) {
      g <- pending(g)
      if (!length(g)) break
      wait_for_cds()
      message(sprintf("Group %d/%d, attempt %d: %d request(s)",
                      gi, length(groups), attempt, length(g)))
      res <- tryCatch(
        wf_request_batch(g, workers = length(g), user = cds_user, path = RAW_DIR,
                         time_out = 6 * 3600),
        error = function(e) { message("Group reported: ", conditionMessage(e)); NULL })
      tg  <- vapply(g, `[[`, "", "target")
      got <- if (is.null(res)) character(0) else as.character(unlist(res))
      if (length(got) == length(tg))
        for (i in seq_along(tg))
          if (!file.exists(file.path(RAW_DIR, tg[i])) && file.exists(got[i]))
            file.rename(got[i], file.path(RAW_DIR, tg[i]))
      if (!length(pending(g))) break
      Sys.sleep(30)
    }
  }
  tg_all <- vapply(todo, `[[`, "", "target")
  is_good <- vapply(todo, raw_request_valid, logical(1))
  invisible(tg_all[!is_good])
}

# 1a. Largest wind request CDS accepts (remembered in a small text file)
WIND_K_FILE <- file.path(RAW_DIR, sprintf("_wind_months_per_request_h%02d.txt", WIND_HOUR_STEP))
probe_wind_size <- function() {
  for (k in WIND_MONTHS_CANDIDATES) {
    wait_for_cds()
    res <- tryCatch(wf_request(user = cds_user, request = wind_requests(k)[[1]],
                               transfer = FALSE, path = RAW_DIR),
                    error = function(e) e)
    if (!inherits(res, "error")) {
      try(res$delete(), silent = TRUE)
      message(sprintf("Wind requests of %d month(s) are accepted by CDS.", k))
      return(k)
    }
    msg <- conditionMessage(res)
    if (grepl("too large|cost limit", msg, ignore.case = TRUE)) {
      message(sprintf("Wind requests of %d month(s) exceed the CDS cost limit; trying a smaller size...", k))
      next
    }
    stop("Wind probe request failed for a reason other than size: ", msg, call. = FALSE)
  }
  stop("Even the smallest wind request exceeds the CDS cost limit. Increase WIND_HOUR_STEP (e.g. 6).",
       call. = FALSE)
}
# ---- 2. BUILD ANNUAL NetCDFs --------------------------------------------------
read_checked <- function(f, cands) {
  r <- read_nc_variable(f, cands)
  stamps <- get_axis(r, f)
  if (!isTRUE(terra::compareGeom(r[[1]], TEMPLATE, stopOnError = FALSE, crs = FALSE,
                                 ext = TRUE, rowcol = TRUE, res = TRUE)))
    stop(basename(f), ": grid differs from the project template.", call. = FALSE)
  list(r = r, stamps = stamps)
}

# Daily 2-m wind for one year (all blocks of that year), for a given method
daily_wind_year <- function(year, method = WIND_METHOD) {
  blocks <- Filter(function(b) b$year == year, wind_blocks(WIND_K))
  out <- lapply(blocks, function(b) {
    f <- wind_file(b)
    if (!file.exists(f)) stop("Missing raw wind file: ", basename(f), call. = FALSE)
    u <- read_checked(f, c("u10", "10m_u_component_of_wind"))
    v <- read_checked(f, c("v10", "10m_v_component_of_wind"))
    if (nlyr(u$r) != nlyr(v$r) || any(u$stamps != v$stamps))
      stop(basename(f), ": u10 and v10 time axes differ.", call. = FALSE)
    day <- format(u$stamps, "%Y-%m-%d")
    idx <- as.integer(factor(day))
    d <- if (method == "mean_hourly_speed") {
      terra::tapp(sqrt(u$r^2 + v$r^2), idx, "mean")
    } else {
      sqrt(terra::tapp(u$r, idx, "mean")^2 + terra::tapp(v$r, idx, "mean")^2)
    }
    list(r = d * wind_10m_to_2m_factor, dates = as.Date(sort(unique(day))))
  })
  r     <- do.call(c, lapply(out, `[[`, "r"))
  dates <- do.call(c, lapply(out, `[[`, "dates"))
  n_exp <- as.integer(as.Date(sprintf("%d-12-31", year)) - as.Date(sprintf("%d-01-01", year))) + 1L
  if (length(dates) != n_exp || nlyr(r) != n_exp)
    stop("Year ", year, ": expected ", n_exp, " daily wind layers, got ", nlyr(r), call. = FALSE)
  list(r = r, dates = dates)
}

build_u2 <- function(year) {
  f_out <- file.path(daily_dir, sprintf("u2_mean_daily_%d.nc", year))
  w <- daily_wind_year(year)
  expected <- seq(as.Date(sprintf("%04d-01-01", year)),
                  as.Date(sprintf("%04d-12-31", year)), by = "day")
  if (!identical(as.numeric(w$dates), as.numeric(expected)))
    stop("Year ", year, ": wind dates do not match the full calendar.")
  write_daily_netcdf_from_spatraster(
    w$r, w$dates, f_out, "u2_mean",
    sprintf("Daily mean 2-m wind speed (%s; every %d hour(s))",
            WIND_METHOD, WIND_HOUR_STEP),
    "m s-1", overwrite = OVERWRITE_DAILY)
  invisible("written")
}

# Read only one calendar year's required 00 UTC SSRD stamps rather than holding
# the entire 1950-2026 raster in memory. CDS files may span 10 years.
ssrd_requests_for_year <- function(year) {
  required_years <- c(year, year + 1L)
  reqs <- ssrd_requests()
  Filter(function(r) any(as.integer(r$year) %in% required_years), reqs)
}

build_rs <- function(year) {
  f_out <- file.path(daily_dir, sprintf("Rs_daily_%d.nc", year))
  days <- seq(as.Date(sprintf("%04d-01-01", year)),
              as.Date(sprintf("%04d-12-31", year)), by = "day")
  want <- as.numeric(as.POSIXct(paste(format(days + 1L), "00:00:00"), tz = "UTC"))
  reqs <- ssrd_requests_for_year(year)
  files <- vapply(reqs, function(r) file.path(RAW_DIR, r$target), character(1))
  parts <- lapply(files, function(f) read_checked(f, c(SSRD_VAR, "ssrd")))
  stamps <- unlist(lapply(parts, function(p) as.numeric(p$stamps)), use.names = FALSE)
  idx <- match(want, stamps)
  if (anyNA(idx))
    stop("Year ", year, ": missing required D+1 00:00 UTC SSRD stamps: ",
         paste(format(head(days[is.na(idx)], 5)), collapse = ", "))
  r <- do.call(c, lapply(parts, `[[`, "r"))
  rs <- terra::clamp(r[[idx]] / 1e6, lower = 0, values = TRUE)
  write_daily_netcdf_from_spatraster(
    rs, days, f_out, "Rs",
    "Daily surface solar radiation downwards (D+1 00 UTC, ERA5-Land)",
    "MJ m-2 day-1", overwrite = OVERWRITE_DAILY)
  invisible("written")
}

# ---- 3. VALIDATION, SAFE RESUME AND SINGLE DRIVER -----------------------------
# All existing outputs require complete dates and an exact grid. A file lacking
# CDS provenance can be retained only by the explicit REUSE_EXISTING_UNVERIFIED_DAILY
# option, and its original source/methodology remains unknown.
inspect_daily_output <- function(year, variable, ignore_overwrite = FALSE) {
  f <- file.path(daily_dir, sprintf("%s_daily_%d.nc", variable, year))
  if (!file.exists(f)) return(FALSE)
  validate <- tryCatch({
    nc <- ncdf4::nc_open(f)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    workflow_att <- ncdf4::ncatt_get(nc, 0, "source_workflow")
    is_cds_workflow <- isTRUE(workflow_att$hasatt) &&
      identical(as.character(workflow_att$value), "cds_daily_v5_unified")
    if (!is_cds_workflow && !REUSE_EXISTING_UNVERIFIED_DAILY)
      stop("unverified source_workflow (possibly an older GEE output)")
    if (!variable %in% names(nc$var)) stop("missing ", variable, " data variable")
    # Only unified CDS products have comparable wind_method metadata. Existing
    # legacy wind products are kept as-is; their calculation is not certified.
    if (identical(variable, "u2_mean") && is_cds_workflow) {
      method <- ncdf4::ncatt_get(nc, 0, "wind_method")
      step <- ncdf4::ncatt_get(nc, 0, "wind_hour_step")
      if (!isTRUE(method$hasatt) || !identical(as.character(method$value), WIND_METHOD) ||
          !isTRUE(step$hasatt) || as.integer(step$value) != WIND_HOUR_STEP)
        stop("wind method or sample interval has changed")
    }
    expected <- seq(as.Date(sprintf("%04d-01-01", year)),
                    as.Date(sprintf("%04d-12-31", year)), by = "day")
    tt <- nc_time_axis(f)
    if (is.null(tt) || length(tt) != length(expected) || anyNA(tt) ||
        !identical(as.numeric(as.Date(tt, tz = "UTC")), as.numeric(expected)))
      stop("incomplete or invalid daily time axis")
    r <- read_nc_variable(f, variable)
    if (terra::nlyr(r) != length(expected) ||
        !isTRUE(terra::compareGeom(r[[1]], TEMPLATE,
                                   stopOnError = FALSE, crs = FALSE,
                                   ext = TRUE, rowcol = TRUE, res = TRUE)))
      stop("layer count or grid does not match template")
    TRUE
  }, error = function(e) {
    structure(FALSE, reason = conditionMessage(e))
  })
  if (!isTRUE(validate) && (!OVERWRITE_DAILY || ignore_overwrite))
    stop("Existing daily forcing is not valid for safe reuse: ", f,
         "\nReason: ", attr(validate, "reason"),
         "\nBack up or move that file and rerun to regenerate it from CDS.",
         call. = FALSE)
  isTRUE(validate) && (!OVERWRITE_DAILY || ignore_overwrite)
}

check_outputs <- function() {
  for (y in YEARS) {
    for (variable in c("u2_mean", "Rs")) {
      f <- file.path(daily_dir, sprintf("%s_daily_%d.nc", variable, y))
      if (!file.exists(f)) stop("Missing daily forcing: ", f)
      # Post-write verification, independent of OVERWRITE_DAILY.
      ok <- inspect_daily_output(y, variable, ignore_overwrite = TRUE)
      if (!ok) stop("Daily output failed validation: ", f)
    }
  }
  message("Validated daily dates, variables and native grid for ", length(YEARS),
          " annual u2 and ", length(YEARS), " annual Rs files.")
  if (REUSE_EXISTING_UNVERIFIED_DAILY)
    message("Provenance note: any reused legacy files remain unverified as CDS-derived ",
            "and have NOT been relabelled.")
  invisible(TRUE)
}

process_daily_cds_forcing <- function() {
  cat("\n============================================================\n")
  cat("PHASE 2: CDS ERA5-Land -> annual u2/Rs (single daily workflow)\n")
  cat("============================================================\n")
  template_file <- file.path(monthly_dir, "total_precipitation_monthly.nc")
  if (!file.exists(template_file))
    stop("Daily forcing requires the project grid template: ", template_file,
         "\nCreate/provide the existing total_precipitation_monthly.nc first.")
  TEMPLATE <<- terra::rast(template_file)[[1]]
  e <- as.vector(terra::ext(TEMPLATE))  # xmin, xmax, ymin, ymax
  AREA <<- round(c(e[4], e[1], e[3], e[2]), 4)  # N W S E
  message(sprintf("Template %d x %d | CDS area (N,W,S,E) = %s",
                  terra::nrow(TEMPLATE), terra::ncol(TEMPLATE),
                  paste(AREA, collapse = ", ")))
  SMOKE_REQ <<- mk_request(SSRD_VAR, 2020, 1, 1:31, 0, "_smoke_ssrd_2020-01.nc")
  ALL_REQS <<- list()
  register_requests(c(ssrd_requests(), list(SMOKE_REQ)))
  
  # Preflight before ANY large CDS request, to detect a legacy archive safely.
  if (!OVERWRITE_DAILY) {
    status_u <- vapply(YEARS, inspect_daily_output, logical(1), variable = "u2_mean")
    status_r <- vapply(YEARS, inspect_daily_output, logical(1), variable = "Rs")
  } else {
    status_u <- status_r <- rep(FALSE, length(YEARS))
  }
  if (all(status_u) && all(status_r)) {
    message("Daily forcing archive is structurally complete. No download needed; ",
            "existing legacy files, if any, remain unverified as CDS outputs.")
    check_outputs()
    if (RUN_COMPARISON) compare_to_gee_reference()
    return(invisible(TRUE))
  }
  
  dir.create(RAW_DIR, showWarnings = FALSE, recursive = TRUE)
  smoke_test()
  if (file.exists(WIND_K_FILE)) {
    previous_k <- suppressWarnings(as.integer(readLines(WIND_K_FILE, n = 1L)))
    if (length(previous_k) != 1L || is.na(previous_k) ||
        !(previous_k %in% WIND_MONTHS_CANDIDATES))
      stop("Invalid saved wind chunk size: ", WIND_K_FILE)
    WIND_K <<- previous_k
  } else {
    WIND_K <<- probe_wind_size()
    writeLines(as.character(WIND_K), WIND_K_FILE)
  }
  register_requests(wind_requests(WIND_K))
  message("Wind chunks: ", WIND_K, " month(s); wind sampling: every ",
          WIND_HOUR_STEP, " hour(s).")
  
  for (i in seq_along(YEARS)) {
    year <- YEARS[i]
    if (status_u[i] && status_r[i]) {
      message(year, ": existing daily products validated for reuse; skipping")
      next
    }
    reqs <- list()
    if (!status_u[i]) {
      reqs <- c(reqs, Filter(function(b) b$year == year, wind_requests(WIND_K)))
    }
    if (!status_r[i]) reqs <- c(reqs, ssrd_requests_for_year(year))
    # Deduplicate multi-year SSRD chunks selected for adjacent years.
    names(reqs) <- vapply(reqs, `[[`, "", "target")
    reqs <- reqs[!duplicated(names(reqs))]
    missing <- run_requests(reqs)
    if (length(missing))
      stop("CDS downloads incomplete for year ", year, ": ",
           paste(missing, collapse = ", "),
           "\nRerun to resume; completed request files were kept.")
    if (!status_u[i]) build_u2(year)
    if (!status_r[i]) build_rs(year)
    # Require both products to be valid before cleaning any per-year wind raw.
    if (!all(vapply(c("u2_mean", "Rs"), function(v) {
      f <- file.path(daily_dir, sprintf("%s_daily_%d.nc", v, year))
      file.exists(f)
    }, logical(1)))) stop("Missing annual output after year ", year)
    if (CLEAN_RAW_AFTER_SUCCESS && !status_u[i]) {
      for (b in Filter(function(b) b$year == year, wind_blocks(WIND_K)))
        safe_remove(wind_file(b))
    }
    message(year, ": CDS daily forcing completed")
    invisible(gc(verbose = FALSE))
  }
  check_outputs()
  if (RUN_COMPARISON) compare_to_gee_reference()
  if (CLEAN_RAW_AFTER_SUCCESS) {
    for (req in c(ssrd_requests(), wind_requests(WIND_K), list(SMOKE_REQ)))
      safe_remove(file.path(RAW_DIR, req$target))
    message("Successfully cleaned temporary CDS raw files; annual NetCDFs retained.")
  }
  cat("CDS daily output: ", normalizePath(daily_dir, winslash = "/"), "\n", sep = "")
  invisible(TRUE)
}


# ---- OPTIONAL READ-ONLY COMPARISON WITH GEE REFERENCE ARCHIVE ------------------------------------
find_old_file <- function(year, key_regex) {
  if (!dir.exists(REFERENCE_DIR)) return(NA_character_)
  f <- list.files(REFERENCE_DIR, pattern = sprintf("%d.*\\.(nc|tif)$", year), full.names = TRUE)
  f <- f[grepl(key_regex, basename(f), ignore.case = TRUE)]
  if (length(f)) f[1] else NA_character_
}

read_old <- function(f, cands) {
  if (grepl("\\.nc$", f, ignore.case = TRUE))
    tryCatch(read_nc_variable(f, cands), error = function(e) terra::rast(f))
  else terra::rast(f)
}

compare_var <- function(label, f_new, new_cands, f_old, old_cands) {
  rn <- read_nc_variable(f_new, new_cands)
  ro <- read_old(f_old, old_cands)
  if (nlyr(rn) != nlyr(ro)) {
    message(sprintf("%s: layer counts differ (CDS %d vs GEE %d); comparison skipped.",
                    label, nlyr(rn), nlyr(ro)))
    return(invisible(NULL))
  }
  if (!nzchar(terra::crs(rn))) terra::crs(rn) <- "EPSG:4326"
  if (!nzchar(terra::crs(ro))) terra::crs(ro) <- "EPSG:4326"
  if (!isTRUE(terra::compareGeom(rn, ro, stopOnError = FALSE, crs = FALSE)))
    ro <- terra::resample(ro, rn, method = "bilinear")
  
  a <- as.vector(terra::values(rn)); b <- as.vector(terra::values(ro))
  ok <- is.finite(a) & is.finite(b)
  a <- a[ok]; b <- b[ok]
  
  # Domain-mean daily series: correlation at lags -1, 0, +1 day exposes date offsets
  sa <- terra::global(rn, "mean", na.rm = TRUE)[, 1]
  sb <- terra::global(ro, "mean", na.rm = TRUE)[, 1]
  n  <- length(sa)
  lag_cor <- c(
    "-1" = cor(sa[-1], sb[-n], use = "complete.obs"),
    "0"  = cor(sa, sb, use = "complete.obs"),
    "+1" = cor(sa[-n], sb[-1], use = "complete.obs"))
  
  message(sprintf(
    "%s %d | mean CDS %.3f, GEE %.3f | bias %.3f | RMSE %.3f | r %.4f | ratio %.3f",
    label, COMPARISON_YEAR, mean(a), mean(b), mean(a - b), sqrt(mean((a - b)^2)),
    cor(a, b), mean(a) / mean(b)))
  message(sprintf("   domain-mean series correlation by lag (days, CDS vs GEE): %s",
                  paste(sprintf("%s: %.4f", names(lag_cor), lag_cor), collapse = " | ")))
  if (which.max(lag_cor) != 2L)
    message("   NOTE: best correlation is not at lag 0; the two archives may define days differently.")
  invisible(list(bias = mean(a - b), rmse = sqrt(mean((a - b)^2)), r = cor(a, b), lag_cor = lag_cor))
}


compare_to_gee_reference <- function() {
  if (normalizePath(REFERENCE_DIR, mustWork = FALSE) ==
      normalizePath(daily_dir, mustWork = FALSE))
    stop("REFERENCE_DIR must differ from daily_dir to avoid comparing an archive to itself.")
  if (!dir.exists(REFERENCE_DIR)) {
    warning("Optional GEE reference directory does not exist: ", REFERENCE_DIR)
    return(invisible(FALSE))
  }
  y <- COMPARISON_YEAR
  f_old_u2 <- find_old_file(y, "u2|wind")
  f_old_rs <- find_old_file(y, "^rs|_rs|rad|ssrd")
  if (is.na(f_old_u2)) message("Comparison: no old u2 file for ", y, " in ", REFERENCE_DIR, "; skipped.")
  else compare_var("u2 (m s-1)", file.path(daily_dir, sprintf("u2_mean_daily_%d.nc", y)), "u2_mean",
                   f_old_u2, c("u2_mean", "u2", "wind"))
  if (is.na(f_old_rs)) message("Comparison: no old Rs file for ", y, " in ", REFERENCE_DIR, "; skipped.")
  else compare_var("Rs (MJ m-2 d-1)", file.path(daily_dir, sprintf("Rs_daily_%d.nc", y)), "Rs",
                   f_old_rs, c("rs", "ssrd", "radiation"))
}

# ======================= PHASE 3: CIRCULATION-CONTEXT FIELDS ==================
# MANUSCRIPT ALIGNMENT (Methods 3g, "Circulation Context of the 2022-2025
# Drought"): the large-scale atmospheric context is evaluated using 500-hPa
# geopotential height (ridging indicator) and a moisture-transport diagnostic
# (integrated vapor transport / vertically integrated moisture flux). Neither
# field is part of ERA5-Land (a land-surface reanalysis), so this phase downloads
# the required fields from ERA5 proper.
#
# Output contract required by the external circulation-context script
# (14circulation_context_v2.R; NOT part of this framework and unrelated to 14event_ranking.R):
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
  # sqrt(U^2 + V^2) in the external circulation-context script.
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
    "ERA5; processed by Monthly_ERALand_download.R"
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
    
    # Re-open and validate the project output exactly as the external circulation-context script
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
    cat("  IVT magnitude = sqrt(eastward_flux^2 + northward_flux^2); compute in the external circulation-context script\n")
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


# =================================== MAIN =====================================
# Preflight before downloading: inspect the archive in monthly_data_direct/.
# Reuse is explicitly opt-in above, and never adds CDS provenance metadata.
preflight_legacy_daily_archive <- function() {
  if (!RUN_DAILY || OVERWRITE_DAILY) return(invisible(TRUE))
  f <- list.files(daily_dir,
                  pattern = "^(u2_mean|Rs)_daily_[0-9]{4}\\.nc$",
                  full.names = TRUE)
  unverified <- character(0)
  for (p in f) {
    verified <- tryCatch({
      n <- ncdf4::nc_open(p)
      a <- ncdf4::ncatt_get(n, 0, "source_workflow")
      ncdf4::nc_close(n)
      isTRUE(a$hasatt) &&
        identical(as.character(a$value), "cds_daily_v5_unified")
    }, error = function(e) FALSE)
    if (!verified) unverified <- c(unverified, p)
  }
  if (length(unverified)) {
    if (!REUSE_EXISTING_UNVERIFIED_DAILY)
      stop(length(unverified), " existing daily file(s) lack CDS provenance, first: ",
           unverified[1],
           "\nBack up legacy files and set OVERWRITE_DAILY <- TRUE to rebuild, ",
           "or explicitly set REUSE_EXISTING_UNVERIFIED_DAILY <- TRUE for ",
           "date/grid/variable-validated reuse.")
    warning(length(unverified), " daily NetCDF(s) in monthly_data_direct/ lack ",
            "unified CDS provenance (first: ", basename(unverified[1]), "). ",
            "Existing files will be preserved IF date/variable/grid checks pass. ",
            "Their provenance and calculation method are NOT verified.")
  }
  message("Daily output directory: ", normalizePath(daily_dir, winslash = "/"))
  invisible(TRUE)
}
preflight_legacy_daily_archive()
if (RUN_DAILY && OVERWRITE_DAILY)
  message("OVERWRITE_DAILY=TRUE: existing annual files will be replaced only after staging and validation.")

# Only a single CDS-based daily workflow is runnable from this script.
# Existing monthly/circulation workflows are distinct scientific products,
# not competing ways of generating the daily u2/Rs archive.
if (RUN_MONTHLY) process_monthly_state_variables()
if (RUN_DAILY) process_daily_cds_forcing()
if (RUN_CIRCULATION) process_circulation_context_fields()
cat("\nPipeline finished.\n")
