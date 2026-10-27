##############################################
# SCRIPT 2b  — PET CALCULATION
# FAO-56 Penman-Monteith PET + ERA5-Land native potential evaporation (PEV)
#
# EXECUTION ORDER: Run after Monthly_ERALand_download.R.
#   The downloader supplies the atmospheric fields required for the independent
#   FAO-56 PM calculation and the ERA5-Land native potential_evaporation (PEV)
#   product. This script produces both PET products in one reproducible stage.
#
# Scientific distinction:
#   PET_PM   = independently calculated FAO-56 Penman-Monteith reference ET.
#   PET_ERA5 = ERA5-Land native potential evaporation (PEV), documented by ECMWF
#              as open-water/pan evaporation. It is not relabelled as PM ET.
#
# No Thornthwaite or detrended/counterfactual PET branch is generated.
#
# MODIFICATIONS vs. the previous script:
#   1. Wind and solar radiation are read from annual daily NetCDF forcing and
#      aggregated to monthly means for FAO-56 PM PET.
#   2. The previous Thornthwaite calculation is removed completely.
#   3. ERA5-Land native potential_evaporation is loaded from the corrected
#      by-hour-of-day monthly product at 00 UTC. The native signed field is
#      interpreted as a 24-hour accumulation and converted to positive mm/day.
#   4. PM-vs-ERA5-Land diagnostics replace all PM-vs-Thornthwaite diagnostics.
#   5. Missing dewpoint is no longer replaced by an artificial T-2°C value.
#   6. Daily forcing archives are checked for duplicate dates before aggregation.
#   7. NetCDF writes are now fatal if the output file is not successfully created.
#   8. PM and ERA5-Land PEV remain explicitly mm/day; 3SPEI_ERALand.R converts
#      both water-demand products to mm/month using the actual calendar month length.
#   9. Daily forcing archives are hard-validated for duplicated calendar dates.
#  10. The PM algebra is centralized in a canonical helper for exact downstream
#      reproduction.
#
# MODIFICATIONS in this revision (v4, QA pass):
#  11. Fixed a double-execution bug: save_plot() draws each figure twice (PDF
#      + PNG), so the Sen's-slope/Mann-Kendall calculation and its log_event()
#      call in panel E (and the F2e trend fits) were being run and logged
#      twice. These are now computed once, outside the plot closures.
#  12. Corrupted "???" placeholder characters (originally arrows/dashes) in
#      comments and log strings replaced with plain ASCII.
#  13. FAO-56 soil-heat-flux equation citation corrected: the 0.14 x
#      (T_i - T_{i-1}) form used here is Eq.44, not Eq.45; Eq.43 is the
#      0.07 x (T_{i+1} - T_{i-1}) form (not used, since T_{i+1} is unknown).
#  14. Removed the unused top-level `gamma` variable; the canonical PM helper
#      already computes gamma_local internally.
#  15. Added explicit logging of (a) how many cell-months had negative raw
#      VPD before flooring to zero, and (b) how many cell-months triggered
#      the last-resort Hargreaves-style Rs fallback -- both were previously
#      silent.
#  16. Warnings raised by validate_nechako_regions() are now captured and
#      written to the processing log instead of only being visible via
#      warnings() in an interactive session.
#  17. PEV source metadata are validated against the attributes actually written
#      by Monthly_ERALand_download.R.
#  18. The PEV metadata validator does not require a nonexistent CDS `step`
#      request field; the 24-hour accumulation semantics are validated from
#      the retained NetCDF provenance metadata.
#  19. PEV_FILE is resolved from MONTHLY_DIR so execution does not depend on
#      the caller's current working directory.
#  20. PEV conversion is explicitly documented as signed metres per 24-hour
#      accumulation -> positive mm/day; no calendar-month divisor is applied.
#  21. PEV-vs-PM magnitude and energy comparisons are diagnostic (non-fatal);
#      PEV is an open-water/pan product and is not required to equal FAO-56 PM ET0.
#
# MODIFICATIONS in this revision (v5):
#  22. Verified 2026-09-29 (July 1950): the PEV archive matches CDS hourly
#      ERA5-Land (00 UTC) and CDS standard monthly means (21.94 vs 22.05 / 21.94
#      mm/day). The ~7x PEV/PM ratio is therefore a property of ERA5-Land PEV,
#      not a download, sign or unit error. The REVIEW message now says so.
#  23. QC-table label corrected: soil-heat-flux coefficient is FAO-56 Eq.44
#
# MODIFICATIONS in this revision (v6, basin-mean contract):
#  25. FAO-56 Eq. 12 / Eq. 39 temperature basis: saturation vapour pressure is
#      es = [e0(Tmax) + e0(Tmin)] / 2 and net long-wave radiation uses
#      [(Tmax,K)^4 + (Tmin,K)^4] / 2 when a monthly Tmax/Tmin file that holds
#      MONTHLY MEANS OF DAILY Tmax/Tmin is supplied (TMAXMIN_FILE). The file is
#      validated first (units, calendar, grid, Tmax >= Tmin, mean diurnal range,
#      midpoint vs T2m). Monthly ABSOLUTE extremes (min/max of all hours in the
#      month) are rejected: they inflate es far beyond FAO-56. With
#      PM_TEMPERATURE_BASIS = "auto" a rejected/missing file falls back to the
#      mean-T2m formulation with a logged warning; "tmax_tmin" stops instead.
#      The basis actually used is written to the log and to the PET NetCDF.
#  24. All "basin" summary series (PM PET, ERA5-Land PEV, meteorological
#      controls, PEV QA means, elevation) are now restricted to the basin mask
#      (union of the 7 subbasins, cells with coverage > 0) and averaged with
#      fractional-coverage area weights (cell area x coverage). A monthly mean
#      is reported only when finite cells hold >= BASIN_MIN_VALID_AREA_COVERAGE
#      (0.20) of the total basin weight. Previously these were unweighted
#      means over the full 52 x 34 bounding box. min/max/sd columns are now
#      computed over basin-mask cells. Gridded NetCDF outputs are unchanged
#      (full native grid; downstream scripts apply the same mask).
#      (label previously read Eq.45).
##############################################
rm(list = ls())

# ---- Working directory (framework rule) -------------------------------------
# Run every script from the project root: set the working directory there
# BEFORE running (setwd(), or open the RStudio project). All paths below are
# derived from getwd(); no path is hard-coded.
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")

library(terra)
library(ncdf4)
library(lubridate)

# --- 1. CONFIGURATION ---
# The project root is the working directory (WD_PATH <- getwd(), set at the
# top of this script). All downstream paths are derived from it.

MONTHLY_DIR   <- file.path(WD_PATH, "monthly_data_direct")
DAILY_DIR     <- file.path(WD_PATH, "daily_eto_forcing")
PET_CALCS_DIR <- file.path(WD_PATH, "pet_calcs")

# Canonical project utility (shared native-grid/CRS and regional helpers).
source(UTILS_FILE)

# OUTPUT LOCATION POLICY
# -----------------------------------------------------------------------------
# NetCDF rasters  -> MONTHLY_DIR  (consumed by 3SPEI_ERALand.R and 9pr_pet_trends.r)
# Summary CSVs    -> MONTHLY_DIR  (consumed by 12basin_timeseries.R Part 4d;
#                                 also available to external companion analyses
#                                 that are not part of this framework)
# Auxiliary products -> PET_CALCS_DIR
# Diagnostic PDFs -> PET_CALCS_DIR/diagnostics/
# Processing log  -> PET_CALCS_DIR/
# -----------------------------------------------------------------------------


# Create output directories if they do not yet exist
if (!dir.exists(MONTHLY_DIR))   dir.create(MONTHLY_DIR,   recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(DAILY_DIR))     stop("DAILY_DIR does not exist: ", DAILY_DIR, "\nRun Monthly_ERALand_download.R first.")
if (!dir.exists(PET_CALCS_DIR)) dir.create(PET_CALCS_DIR, recursive = TRUE, showWarnings = FALSE)

G_GRAV   <- 9.80665
K_RS     <- 0.16   # Interior coefficient for Nechako Basin

# ---- FAO-56 temperature basis for es and net long-wave radiation (USER SETTINGS)
# "auto"      : use Tmax/Tmin if TMAXMIN_FILE exists and passes validation,
#               otherwise the mean-T2m formulation (logged warning)
# "tmax_tmin" : require a valid Tmax/Tmin file; stop otherwise
# "tmean"     : always use monthly mean T2m (previous behaviour)
PM_TEMPERATURE_BASIS <- "auto"
TMAXMIN_FILE         <- file.path(MONTHLY_DIR, "tmin_tmax_monthly_1950_2025.nc")
TMIN_VAR             <- "tmin"
TMAX_VAR             <- "tmax"
# The file must hold MONTHLY MEANS OF DAILY Tmax / Tmin (FAO-56). Checks that
# separate these from monthly absolute extremes of hourly values:
TMAXMIN_MAX_MEAN_DTR_K      <- 18   # record-mean diurnal range; daily-extreme means are ~5-15 K
TMAXMIN_MAX_MIDPOINT_BIAS_K <- 2    # |mean((Tmax+Tmin)/2 - T2m)|; daily-extreme means are within ~1-2 K
stopifnot(PM_TEMPERATURE_BASIS %in% c("auto", "tmax_tmin", "tmean"))
LOG_FILE <- file.path(PET_CALCS_DIR, "monthly_pet_processing.log")

terraOptions(memfrac = 0.8, progress = 0)

cat("Monthly PET Processing Log (Vectorized)\n", file = LOG_FILE)
cat(paste("Processing started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)

log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, " | ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, " | ", msg))
}

# ==============================================================================
# safe_writeCDF: writeCDF wrapper with pre-emptive delete + retry.
# Prevents "Permission denied" when .nc file is open in QGIS/Panoply (Windows).
# ==============================================================================
safe_writeCDF <- function(raster_obj, filepath, varname, longname, unit,
                          overwrite = TRUE, max_retries = 3) {
  if (file.exists(filepath)) {
    del_ok <- tryCatch({ file.remove(filepath); TRUE }, error = function(e) FALSE)
    if (!del_ok) {
      log_event(sprintf("  ERROR: Cannot delete locked file: %s", basename(filepath)))
      log_event("    Close it in QGIS/Panoply/ArcGIS and re-run.")
      return(invisible(FALSE))
    }
    Sys.sleep(0.15)
  }
  for (attempt in seq_len(max_retries)) {
    ok <- tryCatch({
      suppressMessages(
        writeCDF(raster_obj, filepath, varname = varname,
                 longname = longname, unit = unit, overwrite = overwrite))
      TRUE
    }, error = function(e) {
      log_event(sprintf("  WARNING: Write attempt %d/%d failed: %s",
                        attempt, max_retries, conditionMessage(e)))
      FALSE
    })
    if (isTRUE(ok) && file.exists(filepath)) {
      log_event(sprintf("  NetCDF write successful: %s", basename(filepath)))
      return(invisible(TRUE))
    }
    if (file.exists(filepath)) {
      try(file.remove(filepath), silent = TRUE)
    }
    if (attempt < max_retries) Sys.sleep(1.5 * attempt)
  }
  log_event(sprintf("  ERROR: All %d attempts failed for: %s",
                    max_retries, basename(filepath)))
  invisible(FALSE)
}

# ==============================================================================
# terra compatibility shim
# ==============================================================================
# Some installed terra builds expose time() getters but not the replacement
# method time<-() for SpatRaster.  The workflow only needs stable monthly layer
# ordering, so we retain the calendar vector separately and, where possible,
# stamp it into layer names for provenance.
set_monthly_layer_names <- function(raster_obj, dates, variable_name = NULL) {
  dates <- as.Date(dates)
  if (length(dates) == 0L || anyNA(dates)) return(raster_obj)
  if (!is.null(variable_name) && nzchar(variable_name)) {
    names(raster_obj) <- sprintf("%s_%04d_%02d", variable_name,
                                 as.integer(format(dates, "%Y")),
                                 as.integer(format(dates, "%m")))
  } else {
    names(raster_obj) <- sprintf("%04d_%02d", as.integer(format(dates, "%Y")),
                                 as.integer(format(dates, "%m")))
  }
  attr(raster_obj, "monthly_dates") <- dates
  raster_obj
}

# --- 2. LOAD DATA ---
log_event("Loading ERA5-Land data files...")

tryCatch({
  # The ERA5-Land precipitation grid is the project-wide spatial analysis grid.
  # Scripts 1 and 3 use this raster natively; all PET/meteorological inputs are
  # aligned to it here before any PET calculation is performed.
  precip_grid <- rast(file.path(MONTHLY_DIR, "total_precipitation_monthly.nc"))
  t2m_b <- rast(file.path(MONTHLY_DIR, "2m_temperature_monthly.nc"))
  d2m_b <- rast(file.path(MONTHLY_DIR, "2m_dewpoint_temperature_monthly.nc"))
  sp_b  <- rast(file.path(MONTHLY_DIR, "surface_pressure_monthly.nc"))
  log_event("Successfully loaded monthly precipitation, temperature, dewpoint, and pressure inputs.")
}, error = function(e) stop("Input file loading failed: ", e$message))

# Exact common grid: the precipitation NetCDF used by Script 1 (1SPI_ERALand.R).
# Do not invent a second spatial template in the PET pipeline.
PET_ANALYSIS_TEMPLATE <- precip_grid[[1]]
align_to_common_grid <- function(x, label) {
  if (terra::compareGeom(x[[1]], PET_ANALYSIS_TEMPLATE,
                         stopOnError = FALSE,
                         crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)) {
    log_event(sprintf("  [GRID] %s already matches the common ERA5-Land grid.", label))
    return(x)
  }
  log_event(sprintf("  [GRID] Aligning %s to the common ERA5-Land precipitation grid.", label))
  terra::project(x, PET_ANALYSIS_TEMPLATE, method = "bilinear")
}

t2m_b <- align_to_common_grid(t2m_b, "2-m temperature")
d2m_b <- align_to_common_grid(d2m_b, "2-m dewpoint temperature")
sp_b  <- align_to_common_grid(sp_b,  "surface pressure")

log_event(sprintf("[GRID] Common ERA5-Land analysis grid: %d x %d = %d cells | res %.6f x %.6f | CRS=%s",
                  terra::nrow(PET_ANALYSIS_TEMPLATE), terra::ncol(PET_ANALYSIS_TEMPLATE),
                  terra::ncell(PET_ANALYSIS_TEMPLATE),
                  terra::res(PET_ANALYSIS_TEMPLATE)[1], terra::res(PET_ANALYSIS_TEMPLATE)[2],
                  terra::crs(PET_ANALYSIS_TEMPLATE)))

# ------------------------------------------------------------------------------
# Basin mask and fractional-coverage weights (shared framework rule).
# Basin = dissolved union of the 7 subbasins (as in 3SPEI_ERALand.R).
# ------------------------------------------------------------------------------
PET_BASIN_BW <- basin_fractional_weights(
  PET_ANALYSIS_TEMPLATE,
  nechako_subbasin_union(load_nechako_subbasins())
)
log_event(sprintf(
  "[BASIN] %d basin cells (%d partial) | weighted area %.1f km2 | min valid coverage %.2f",
  PET_BASIN_BW$n_cells, PET_BASIN_BW$n_partial, PET_BASIN_BW$area_m2 / 1e6,
  BASIN_MIN_VALID_AREA_COVERAGE
))

# Per-layer basin statistics: weighted mean (shared 0.20 coverage rule) plus
# min / max / sd over basin-mask cells. Column names match terra::global().
basin_layer_summary <- function(r, bw = PET_BASIN_BW) {
  v <- terra::values(r, mat = TRUE)
  if (nrow(v) != length(bw$weights))
    stop("basin_layer_summary(): raster cell count does not match the basin weights.")
  vb <- v[bw$mask, , drop = FALSE]
  col_stat <- function(f) apply(vb, 2L, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) f(x) else NA_real_
  })
  data.frame(
    mean = basin_weighted_mean(v, bw),
    min  = col_stat(min),
    max  = col_stat(max),
    sd   = col_stat(function(x) if (length(x) > 1L) stats::sd(x) else NA_real_)
  )
}

# Elevation is required only for FAO-56 clear-sky radiation Rso.  The revised
# downloader does not retain a geopotential time series.  If an existing
# geopotential_monthly.nc is available, use it; otherwise estimate a static
# elevation field from long-term mean surface pressure using the inverse of
# FAO-56 Eq. 7.  This avoids adding another large monthly archive solely for Rso.
z_file <- file.path(MONTHLY_DIR, "geopotential_monthly.nc")
if (file.exists(z_file)) {
  log_event("Using existing geopotential_monthly.nc for elevation.")
  z_b <- rast(z_file)
  h_m <- z_b[[1]] / G_GRAV
  h_m <- resample(h_m, t2m_b[[1]], method = "bilinear")
  rm(z_b)
} else {
  log_event("geopotential_monthly.nc not found; estimating static elevation from climatological surface pressure.")
  p_mean_kpa <- mean(sp_b, na.rm = TRUE) / 1000
  h_m <- (293 / 0.0065) * (1 - (p_mean_kpa / 101.3)^(1 / 5.26))
  h_m <- clamp(h_m, lower = -500, upper = 6000, values = TRUE)
  rm(p_mean_kpa)
}
avg_elevation <- basin_weighted_mean(terra::values(h_m, mat = TRUE), PET_BASIN_BW)[1]
log_event(paste("Average basin elevation used for Rso:", round(avg_elevation, 1), "m"))

# --- 3. ROBUST DATE HANDLING ---
log_event("Processing time information...")

extract_time_robust <- function(raster_obj, file_path, var_name) {
  time_info <- time(raster_obj)
  if ((inherits(time_info, "Date") || inherits(time_info, "POSIXt")) &&
      !all(is.na(time_info))) {
    log_event(sprintf("  Time extracted via terra::time() for %s", var_name))
    return(as.Date(time_info))
  }
  log_event(sprintf("  Trying ncdf4 for %s...", var_name))
  tryCatch({
    nc         <- nc_open(file_path)
    time_data  <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")$value
    m          <- regexpr("since\\s+([0-9]{4}-[0-9]{2}-[0-9]{2})", time_units)
    if (m > 0) {
      origin  <- as.Date(sub("since\\s+", "", regmatches(time_units, m)))
      tv      <- if (grepl("hours", time_units, ignore.case = TRUE))
        origin + time_data / 24 else origin + time_data
      nc_close(nc)
      log_event(sprintf("  Time extracted via ncdf4 for %s", var_name))
      return(as.Date(tv))
    }
    nc_close(nc)
  }, error = function(e)
    log_event(sprintf("  ncdf4 failed for %s: %s", var_name, conditionMessage(e))))
  n  <- nlyr(raster_obj)
  tv <- seq.Date(as.Date("1950-01-01"), by = "month", length.out = n)
  log_event(sprintf("  Reconstructed time for %s (%d layers from 1950-01)", var_name, n))
  tv
}

dates    <- extract_time_robust(t2m_b, file.path(MONTHLY_DIR, "2m_temperature_monthly.nc"), "Temperature")
n_months <- nlyr(t2m_b)
log_event(paste("Processing", n_months, "time steps."))
log_event(sprintf("Time range: %s to %s", format(min(dates)), format(max(dates))))

month_indices <- as.numeric(format(dates, "%m"))

# ------------------------------------------------------------------------------
# Validate monthly state-variable calendars before any arithmetic.  The PM
# equation assumes that T2m, dewpoint, and surface pressure refer to the same
# calendar month in the same layer position.
# ------------------------------------------------------------------------------
d2m_dates <- extract_time_robust(
  d2m_b,
  file.path(MONTHLY_DIR, "2m_dewpoint_temperature_monthly.nc"),
  "Dewpoint temperature"
)
sp_dates <- extract_time_robust(
  sp_b,
  file.path(MONTHLY_DIR, "surface_pressure_monthly.nc"),
  "Surface pressure"
)

normalize_month_dates <- function(x) as.Date(format(as.Date(x), "%Y-%m-01"))

dates_norm    <- normalize_month_dates(dates)
d2m_dates_norm <- normalize_month_dates(d2m_dates)
sp_dates_norm  <- normalize_month_dates(sp_dates)

if (nlyr(d2m_b) != n_months || nlyr(sp_b) != n_months) {
  stop(
    "Monthly state-variable layer mismatch: T2m=", n_months,
    ", dewpoint=", nlyr(d2m_b),
    ", surface pressure=", nlyr(sp_b), "."
  )
}

if (!identical(dates_norm, d2m_dates_norm) ||
    !identical(dates_norm, sp_dates_norm)) {
  stop(
    "Monthly state-variable calendars are not identical. ",
    "T2m, dewpoint, and surface pressure must have the same YYYY-MM sequence."
  )
}

dates <- dates_norm
month_indices <- as.integer(format(dates, "%m"))

# ------------------------------------------------------------------------------
# ERA5-Land native potential evaporation (PEV)
# ------------------------------------------------------------------------------
# ERA5-Land provides `potential_evaporation` (pev) as a ready-to-use
# potential-evaporation product. It is scientifically distinct from the
# independently calculated FAO-56 Penman-Monteith reference ET below: ECMWF
# documents ERA5-Land PEV as an open-water/pan evaporation quantity. We retain
# that distinction explicitly in filenames, labels, and downstream analyses.
#
# The downloader retrieves PEV from:
#   product_type = monthly_averaged_reanalysis_by_hour_of_day
#   time         = 00:00 UTC
#
# The retained native field is documented as a 24-hour accumulation ending at
# 00 UTC, represented as a monthly mean of daily accumulations. Script 2b then
# converts the signed native water depth from metres to positive mm/day.
PEV_FILE <- file.path(MONTHLY_DIR, "potential_evaporation_monthly.nc")
if (!file.exists(PEV_FILE)) {
  stop(
    "ERA5-Land native potential evaporation file not found: ", PEV_FILE,
    "\nRun Monthly_ERALand_download.R first."
  )
}

log_event("Loading ERA5-Land native potential evaporation (PEV)...")

# Validate provenance before reading values. This prevents an older or
# differently constructed PEV archive from silently entering the PET workflow.
# The validator uses the metadata names actually written by the downloader.
validate_pev_source_metadata <- function(file_path) {
  
  nc <- tryCatch(
    ncdf4::nc_open(file_path),
    error = function(e) NULL
  )
  
  if (is.null(nc)) {
    stop(
      "Cannot open ERA5-Land PEV NetCDF for metadata validation: ",
      file_path
    )
  }
  
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  get_global_att <- function(name) {
    a <- tryCatch(
      ncdf4::ncatt_get(nc, 0, name),
      error = function(e) NULL
    )
    
    if (is.null(a) || !isTRUE(a$hasatt)) {
      return(NA_character_)
    }
    
    as.character(a$value)
  }
  
  product_type <- get_global_att("pev_source_product")
  source_time  <- get_global_att("pev_source_time")
  accumulation <- get_global_att("pev_accumulation_interpretation")
  native_units <- get_global_att("pev_native_units")
  sign_conv    <- get_global_att("pev_sign_convention")
  downstream   <- get_global_att("pev_downstream_conversion")
  
  log_event(sprintf(
    paste0(
      "PEV source metadata: product=%s; time=%s; ",
      "native_units=%s; accumulation=%s"
    ),
    product_type,
    source_time,
    native_units,
    accumulation
  ))
  
  if (!identical(
    product_type,
    "monthly_averaged_reanalysis_by_hour_of_day"
  )) {
    stop(
      "PEV source product metadata are missing or incorrect. Expected ",
      "'monthly_averaged_reanalysis_by_hour_of_day'; got: ",
      product_type,
      ". Re-run Monthly_ERALand_download.R."
    )
  }
  
  if (!identical(source_time, "00:00 UTC")) {
    stop(
      "PEV source time metadata are missing or incorrect. Expected ",
      "'00:00 UTC'; got: ",
      source_time,
      ". Re-run Monthly_ERALand_download.R."
    )
  }
  
  if (!identical(
    native_units,
    "m per 24-hour accumulation"
  )) {
    stop(
      "PEV native-units metadata are missing or incorrect. Expected ",
      "'m per 24-hour accumulation'; got: ",
      native_units,
      ". Re-run Monthly_ERALand_download.R."
    )
  }
  
  if (!identical(
    accumulation,
    "24-hour accumulation ending at 00 UTC; monthly mean of daily accumulations"
  )) {
    stop(
      "PEV accumulation metadata are missing or inconsistent. Got: ",
      accumulation,
      ". Re-run Monthly_ERALand_download.R."
    )
  }
  
  if (!is.character(sign_conv) ||
      !grepl("ERA5 evaporation is negative", sign_conv, fixed = TRUE)) {
    stop(
      "PEV sign-convention metadata are missing or inconsistent. ",
      "Re-run Monthly_ERALand_download.R."
    )
  }
  
  if (!identical(
    downstream,
    "PEV_mm_day = -PEV_native_m * 1000"
  )) {
    stop(
      "PEV downstream-conversion metadata are missing or inconsistent. ",
      "Expected 'PEV_mm_day = -PEV_native_m * 1000'; got: ",
      downstream,
      ". Re-run Monthly_ERALand_download.R."
    )
  }
  
  log_event(
    paste(
      "PEV source metadata validation PASSED:",
      "by-hour-of-day product;",
      "00:00 UTC;",
      "24-hour accumulation;",
      "native metres;",
      "conversion = -native * 1000."
    )
  )
  
  invisible(TRUE)
}

validate_pev_source_metadata(PEV_FILE)

pev_raw <- rast(PEV_FILE)
pev_dates <- extract_time_robust(
  pev_raw,
  PEV_FILE,
  "potential_evaporation"
)
pev_dates <- as.Date(format(pev_dates, "%Y-%m-01"))
pev_idx <- match(dates, pev_dates)

if (anyNA(pev_idx)) {
  stop(
    "ERA5-Land PEV does not cover the complete PET analysis period. Missing: ",
    paste(format(dates[is.na(pev_idx)], "%Y-%m"), collapse = ", ")
  )
}

pev_raw <- pev_raw[[pev_idx]]

# ERA5-Land PEV native-grid handling
# ----------------------------------
# The downloaded native PEV diagnostic establishes:
#   * 52 x 34 cells and 912 monthly layers
#   * WGS 84 / CRS84
#   * identical native geometry to the ERA5-Land temperature grid
#   * native values predominantly <= 0, with a small positive tail
#
# The native positive tail is the source of the negative values produced by
# the ECMWF evaporation-sign conversion below.  Positive values are therefore
# explicitly treated as zero BEFORE conversion.  This is not a tolerance-based
# repair of the converted product and does not alter valid negative PEV values.
#
# No projection/resampling is performed when the native geometry already
# matches t2m_b.  A conditional fallback remains for future archives whose
# geometry genuinely differs.

raw_pev_stats <- global(
  pev_raw,
  c("min", "max", "mean", "sd"),
  na.rm = TRUE
)

raw_negative_count <- sum(
  global(pev_raw < 0, "sum", na.rm = TRUE)[, 1],
  na.rm = TRUE
)

raw_positive_count <- sum(
  global(pev_raw > 0, "sum", na.rm = TRUE)[, 1],
  na.rm = TRUE
)

raw_zero_count <- sum(
  global(pev_raw == 0, "sum", na.rm = TRUE)[, 1],
  na.rm = TRUE
)

raw_max_positive_by_layer <- global(
  ifel(pev_raw > 0, pev_raw, NA),
  "max",
  na.rm = TRUE
)[, 1]

raw_max_positive <- if (any(is.finite(raw_max_positive_by_layer))) {
  max(raw_max_positive_by_layer[is.finite(raw_max_positive_by_layer)])
} else {
  NA_real_
}

log_event(sprintf(
  "Raw ERA5-Land PEV: %.6f to %.6f m (mean across layers: %.6f)",
  min(raw_pev_stats$min, na.rm = TRUE),
  max(raw_pev_stats$max, na.rm = TRUE),
  mean(raw_pev_stats$mean, na.rm = TRUE)
))

log_event(sprintf(
  "Raw PEV sign counts across all layers: negative=%d, positive=%d, zero=%d",
  raw_negative_count,
  raw_positive_count,
  raw_zero_count
))

if (is.finite(raw_max_positive)) {
  log_event(sprintf(
    "Maximum positive native PEV: %.9f m (%.6f mm/day)",
    raw_max_positive,
    raw_max_positive * 1000
  ))
}

# The observed native positive tail is a sign-inconsistent value for the
# ECMWF evaporation convention.  Set ONLY positive native values to zero.
if (raw_positive_count > 0) {
  pev_raw <- ifel(pev_raw > 0, 0, pev_raw)
  log_event(sprintf(
    "Corrected %d positive native PEV cells to zero before sign conversion.",
    raw_positive_count
  ))
}

# ERA5-Land PEV is retained from the monthly averaged reanalysis by-hour-of-day
# product (CDS stream=mnth) at 00:00 UTC. Per ECMWF's ERA5-Land accumulation
# convention, CDS validity time 00:00 UTC maps to MARS step=24, i.e. the native
# field is a 24-hour accumulated water depth ending at 00 UTC, expressed in
# metres; the monthly field is the mean of those daily 24-hour accumulations.
# (https://confluence.ecmwf.int/display/CKB/ERA5-Land:+data+documentation,
# sections "Accumulations" and "Monthly means".)
#
# ERA5-Land uses the ECMWF evaporation sign convention, with evaporation
# negative. Convert to positive water depth in mm/day. No calendar-month
# scaling is applied.
pev_era5 <- -pev_raw * 1000

# Hard validation: after native positive-tail correction, negative evaporative
# demand must not exist. Do not silently clamp unexpected residual negatives.
negative_after_conversion <- global(
  pev_era5 < 0,
  "sum",
  na.rm = TRUE
)[, 1]

negative_after_conversion <- sum(
  negative_after_conversion,
  na.rm = TRUE
)

if (negative_after_conversion > 0) {
  most_negative <- global(
    pev_era5,
    "min",
    na.rm = TRUE
  )[, 1]
  most_negative <- min(
    most_negative[is.finite(most_negative)],
    na.rm = TRUE
  )
  stop(
    "Negative ERA5-Land PEV remains after native sign conversion: ",
    negative_after_conversion,
    " cells; minimum = ",
    sprintf("%.9f mm/day", most_negative),
    ". Inspect the native PEV input."
  )
}

# Native PEV and T2m are already on the same ERA5-Land grid. Only align if a
# future input archive actually differs in CRS or geometry.
if (!terra::compareGeom(pev_era5[[1]], PET_ANALYSIS_TEMPLATE,
                        stopOnError = FALSE, crs = TRUE, ext = TRUE,
                        rowcol = TRUE, res = TRUE)) {
  log_event("PEV geometry differs from the common ERA5-Land grid; projecting to common grid...")
  pev_era5 <- terra::project(pev_era5, PET_ANALYSIS_TEMPLATE, method = "bilinear")
}

# Final hard validation after optional spatial alignment.
final_negative_count <- sum(
  global(pev_era5 < 0, "sum", na.rm = TRUE)[, 1],
  na.rm = TRUE
)

if (final_negative_count > 0) {
  most_negative <- global(pev_era5, "min", na.rm = TRUE)[, 1]
  most_negative <- min(most_negative[is.finite(most_negative)], na.rm = TRUE)
  stop(
    "Negative values remain in final ERA5-Land PEV after spatial alignment: ",
    final_negative_count,
    " cells; minimum = ",
    sprintf("%.9f mm/day", most_negative),
    "."
  )
}

pev_era5 <- set_monthly_layer_names(pev_era5, dates, "PET_ERA5LAND")

pev_stats <- global(
  pev_era5,
  c("min", "max", "mean", "sd"),
  na.rm = TRUE
)

log_event(sprintf(
  "Final ERA5-Land PEV: %.3f to %.3f mm/day (mean: %.3f)",
  min(pev_stats$min, na.rm = TRUE),
  max(pev_stats$max, na.rm = TRUE),
  mean(pev_stats$mean, na.rm = TRUE)
))

# ------------------------------------------------------------------------------
# 3.25. DAILY FORCING -> MONTHLY MEANS
# ------------------------------------------------------------------------------
# The revised downloader provides daily 2-m wind-speed magnitude and daily
# shortwave radiation.  PM remains a monthly formulation here: the daily scalar
# forcing is aggregated to monthly mean daily values before entering the
# monthly FAO-56 equation.
aggregate_daily_archive_to_monthly <- function(variable_name,
                                               expected_dates,
                                               template) {
  years_needed <- sort(unique(as.integer(format(expected_dates, "%Y"))))
  monthly_parts <- vector("list", length(years_needed))
  monthly_dates <- vector("list", length(years_needed))
  
  for (j in seq_along(years_needed)) {
    yr <- years_needed[j]
    f <- file.path(DAILY_DIR, sprintf("%s_daily_%04d.nc", variable_name, yr))
    
    if (!file.exists(f)) {
      stop("Missing revised daily forcing file: ", f)
    }
    
    x <- terra::rast(f)
    tt <- terra::time(x)
    
    if (!(inherits(tt, "Date") || inherits(tt, "POSIXt")) || all(is.na(tt))) {
      stop("No valid NetCDF time axis in: ", f)
    }
    
    dd <- as.Date(tt)
    
    # Hard temporal invariant: exactly one layer per calendar date.
    if (anyDuplicated(dd)) {
      dup_dates <- unique(dd[duplicated(dd)])
      stop(
        variable_name,
        " contains duplicated daily date(s) in ",
        yr,
        ": ",
        paste(format(dup_dates), collapse = ", ")
      )
    }
    
    keep <- which(as.integer(format(dd, "%Y")) == yr)
    if (!length(keep)) {
      stop("No dates for year ", yr, " in ", f)
    }
    
    x <- x[[keep]]
    dd <- dd[keep]
    
    expected_daily <- seq(
      as.Date(sprintf("%04d-01-01", yr)),
      as.Date(sprintf("%04d-12-31", yr)),
      by = "day"
    )
    
    expected_daily <- expected_daily[
      expected_daily >= min(expected_dates) &
        expected_daily <= (max(expected_dates) + months(1) - days(1))
    ]
    
    missing_daily <- setdiff(expected_daily, dd)
    if (length(missing_daily)) {
      stop(
        variable_name,
        " has missing daily dates in ",
        yr,
        ": ",
        paste(head(missing_daily, 10), collapse = ", ")
      )
    }
    
    month_index <- as.integer(format(dd, "%m"))
    xm <- terra::tapp(x, index = month_index, fun = mean, na.rm = TRUE)
    md <- as.Date(paste0(sort(unique(format(dd, "%Y-%m"))), "-01"))
    
    if (terra::nlyr(xm) != length(md)) {
      stop(
        variable_name,
        ": monthly aggregation produced ",
        terra::nlyr(xm),
        " layers but ",
        length(md),
        " monthly dates."
      )
    }
    
    monthly_parts[[j]] <- xm
    monthly_dates[[j]] <- md
    rm(x, xm)
    gc(verbose = FALSE)
  }
  
  out <- do.call(c, monthly_parts)
  od <- do.call(c, monthly_dates)
  idx <- match(expected_dates, od)
  
  if (anyNA(idx)) {
    stop(
      variable_name,
      " monthly aggregation is missing: ",
      paste(format(expected_dates[is.na(idx)], "%Y-%m"), collapse = ", ")
    )
  }
  
  out <- out[[idx]]
  
  if (!terra::compareGeom(out[[1]], template[[1]], stopOnError = FALSE)) {
    log_event(sprintf(
      "  Resampling %s monthly forcing to temperature grid.",
      variable_name
    ))
    out <- terra::resample(out, template[[1]], method = "bilinear")
  }
  
  out
}

log_event("Aggregating revised daily u2 and Rs archives to monthly means...")
u2 <- aggregate_daily_archive_to_monthly("u2_mean", dates, t2m_b)
Rs <- aggregate_daily_archive_to_monthly("Rs", dates, t2m_b)
names(u2) <- sprintf(
  "u2_%04d_%02d",
  as.integer(format(dates, "%Y")),
  as.integer(format(dates, "%m"))
)
names(Rs) <- sprintf(
  "Rs_%04d_%02d",
  as.integer(format(dates, "%Y")),
  as.integer(format(dates, "%m"))
)

# --- 3.5. DAYS IN MONTH ---
log_event("Calculating days in each month...")
first_of_month <- as.Date(format(dates, "%Y-%m-01"))
first_next_month <- seq(
  first_of_month[1],
  by = "month",
  length.out = n_months + 1
)[-1]
days_in_month <- as.integer(first_next_month - first_of_month)
log_event(sprintf(
  "Days range: %d to %d",
  min(days_in_month),
  max(days_in_month)
))

# --- 4. PRE-CALCULATIONS ---
log_event("Calculating Ra (Extraterrestrial Radiation)...")
lat_r <- init(t2m_b, "y")
lat_rad <- lat_r * pi / 180

mid_month_days <- c(16, 46, 75, 105, 135, 166, 196, 227, 258, 288, 319, 349)
ra_list <- vector("list", 12)

for (m in 1:12) {
  J <- mid_month_days[m]
  delta <- 0.409 * sin(2 * pi * J / 365 - 1.39)
  dr <- 1 + 0.033 * cos(2 * pi * J / 365)
  tan_prod <- ifel(
    tan(lat_rad) * tan(delta) > 1,
    1,
    ifel(
      tan(lat_rad) * tan(delta) < -1,
      -1,
      tan(lat_rad) * tan(delta)
    )
  )
  ws <- acos(-tan_prod)
  Ra_val <- (24 * 60 / pi) * 0.0820 * dr *
    (ws * sin(lat_rad) * sin(delta) +
       cos(lat_rad) * cos(delta) * sin(ws))
  ra_list[[m]] <- ifel(Ra_val < 0, 0, Ra_val)
}

ra_clim <- rast(ra_list)
Ra_stack <- ra_clim[[month_indices]]

log_event("Ra climatology statistics:")
for (m in 1:12) {
  s_ra <- global(ra_clim[[m]], c("min", "max", "mean"), na.rm = TRUE)
  log_event(sprintf(
    "  Month %2d: min=%.2f, max=%.2f, mean=%.2f MJ/m2/day",
    m,
    s_ra$min,
    s_ra$max,
    s_ra$mean
  ))
}
log_event("Ra stack generated.")

# --- 5. PENMAN-MONTEITH INPUT FIELDS ---
log_event("Preparing FAO-56 Penman-Monteith input fields...")

# A. Temperature and pressure
# NOTE: psychrometric constant (gamma) is intentionally NOT computed here.
# It is computed once, internally, inside the canonical fa056_pm_et0() helper
# below (as gamma_local) so there is a single authoritative place where the
# PM algebra lives, per modification #10 in the header.
T_mean <- t2m_b - 273.15
P_kpa <- sp_b / 1000

s <- global(T_mean, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf(
  "Temperature range: %.2f to %.2f degC (mean: %.2f)",
  min(s$min),
  max(s$max),
  mean(s$mean)
))

# A2. Tmax / Tmin for FAO-56 Eq. 12 (es) and Eq. 39 (net long-wave radiation)
load_validated_tmax_tmin <- function() {
  reject <- function(msg) {
    if (PM_TEMPERATURE_BASIS == "tmax_tmin")
      stop("PM_TEMPERATURE_BASIS = 'tmax_tmin' but the Tmax/Tmin input was rejected: ", msg, call. = FALSE)
    log_event(paste0("WARNING: Tmax/Tmin input not used (", msg,
                     "). es and Rnl fall back to monthly mean T2m."))
    NULL
  }
  if (PM_TEMPERATURE_BASIS == "tmean") {
    log_event("PM temperature basis: mean T2m (PM_TEMPERATURE_BASIS = 'tmean').")
    return(NULL)
  }
  if (!file.exists(TMAXMIN_FILE)) return(reject(paste("file not found:", TMAXMIN_FILE)))
  tn <- tryCatch(terra::rast(TMAXMIN_FILE, subds = TMIN_VAR), error = function(e) NULL)
  tx <- tryCatch(terra::rast(TMAXMIN_FILE, subds = TMAX_VAR), error = function(e) NULL)
  if (is.null(tn) || is.null(tx)) return(reject(sprintf("variables '%s'/'%s' not found", TMIN_VAR, TMAX_VAR)))
  if (terra::nlyr(tn) != terra::nlyr(tx)) return(reject("tmin and tmax have different layer counts"))
  
  # Calendar: same YYYY-MM sequence as T2m (subset/reorder by month, never shift)
  td <- as.Date(format(as.Date(terra::time(tn)), "%Y-%m-01"))
  if (anyNA(td) || anyDuplicated(td)) return(reject("missing or duplicated time stamps"))
  idx <- match(dates, td)
  if (anyNA(idx)) return(reject(sprintf("%d T2m months are missing from the Tmax/Tmin file", sum(is.na(idx)))))
  tn <- tn[[idx]]; tx <- tx[[idx]]
  
  # Units: attribute if present, otherwise magnitude (K vs degC)
  u <- tolower(trimws(c(terra::units(tn)[1], terra::units(tx)[1])))
  mean_raw <- mean(terra::global(tn[[seq_len(min(24L, terra::nlyr(tn)))]], "mean", na.rm = TRUE)[, 1])
  is_kelvin <- if (all(u %in% c("k", "kelvin"))) TRUE else
    if (all(u %in% c("degc", "c", "celsius", "deg_c", "degrees_celsius"))) FALSE else mean_raw > 150
  if (is_kelvin) { tn <- tn - 273.15; tx <- tx - 273.15 }
  log_event(sprintf("  Tmax/Tmin units: %s (%s)", if (is_kelvin) "K -> degC" else "degC",
                    if (all(nzchar(u)) && !anyNA(u)) "metadata" else "inferred from magnitude"))
  
  tn <- align_to_common_grid(tn, "monthly Tmin")
  tx <- align_to_common_grid(tx, "monthly Tmax")
  
  # Physical / statistical validation
  dtr <- tx - tn
  n_neg <- sum(terra::global(dtr < -0.05, "sum", na.rm = TRUE)[, 1], na.rm = TRUE)
  if (n_neg > 0) return(reject(sprintf("%d cell-months with Tmax < Tmin", n_neg)))
  mean_dtr <- mean(terra::global(dtr, "mean", na.rm = TRUE)[, 1], na.rm = TRUE)
  mid_bias <- mean(terra::global((tx + tn) / 2 - T_mean, "mean", na.rm = TRUE)[, 1], na.rm = TRUE)
  mon <- as.integer(format(dates, "%m"))
  dtr_l <- terra::global(dtr, "mean", na.rm = TRUE)[, 1]
  log_event(sprintf(paste0("  Tmax/Tmin checks: record-mean DTR %.1f K (limit %.0f) | ",
                           "mean (Tmax+Tmin)/2 - T2m %+.2f K (limit +/-%.1f) | DJF DTR %.1f K vs JJA %.1f K"),
                    mean_dtr, TMAXMIN_MAX_MEAN_DTR_K, mid_bias, TMAXMIN_MAX_MIDPOINT_BIAS_K,
                    mean(dtr_l[mon %in% c(12, 1, 2)], na.rm = TRUE), mean(dtr_l[mon %in% 6:8], na.rm = TRUE)))
  if (!is.finite(mean_dtr) || mean_dtr > TMAXMIN_MAX_MEAN_DTR_K)
    return(reject(sprintf(paste0("record-mean diurnal range %.1f K exceeds %.0f K; the file looks like monthly ",
                                 "ABSOLUTE extremes, not monthly means of daily Tmax/Tmin"),
                          mean_dtr, TMAXMIN_MAX_MEAN_DTR_K)))
  if (!is.finite(mid_bias) || abs(mid_bias) > TMAXMIN_MAX_MIDPOINT_BIAS_K)
    return(reject(sprintf("mean (Tmax+Tmin)/2 differs from T2m by %+.2f K", mid_bias)))
  n_na <- sum(terra::global(is.na(tn) & !is.na(T_mean), "sum", na.rm = TRUE)[, 1], na.rm = TRUE)
  if (n_na > 0) log_event(sprintf("  NOTE: %d cell-months lack Tmax/Tmin; those cells use mean T2m.", n_na))
  log_event("PM temperature basis: FAO-56 Tmax/Tmin (es Eq. 12, Rnl Eq. 39).")
  list(tmin = tn, tmax = tx, mean_dtr = mean_dtr, midpoint_bias = mid_bias)
}
TXN <- load_validated_tmax_tmin()
PM_TEMPERATURE_BASIS_USED <- if (is.null(TXN)) "mean T2m" else "FAO-56 Tmax/Tmin"

# B. Vapour pressure
# No artificial dewpoint substitution is used. Missing ERA5-Land dewpoint
# therefore propagates as NA into ea/VPD and ultimately into PM PET.
Tdew <- d2m_b - 273.15
missing_dewpoint <- sum(global(is.na(Tdew), "sum", na.rm = FALSE)[, 1], na.rm = TRUE)
if (missing_dewpoint > 0) {
  log_event(sprintf(
    "WARNING: %d dewpoint cells are missing; no artificial dewpoint replacement is applied.",
    missing_dewpoint
  ))
}

e0_kpa <- function(T) 0.6108 * exp((17.27 * T) / (T + 237.3))   # FAO-56 Eq. 11
es_tmean <- e0_kpa(T_mean)
es <- if (is.null(TXN)) es_tmean else {
  # FAO-56 Eq. 12; any cell-month without Tmax/Tmin keeps the mean-T2m value
  es_tx <- (e0_kpa(TXN$tmax) + e0_kpa(TXN$tmin)) / 2
  terra::cover(es_tx, es_tmean)
}
ea <- e0_kpa(Tdew)
if (!is.null(TXN)) {
  r_es <- mean(terra::global(es / es_tmean, "mean", na.rm = TRUE)[, 1], na.rm = TRUE)
  log_event(sprintf("  es (Tmax/Tmin) / es (mean T2m): record mean %.3f", r_es))
}
vpd_raw <- es - ea

# Cells where Tdew > T2m give a physically invalid negative VPD (supersaturation
# artifact in the reanalysis fields). These are floored to zero, but the count
# is logged explicitly so a systematic dewpoint/temperature inconsistency does
# not pass through silently.
n_negative_vpd <- sum(global(vpd_raw < 0, "sum", na.rm = TRUE)[, 1], na.rm = TRUE)
if (n_negative_vpd > 0) {
  log_event(sprintf(
    "WARNING: %d cell-months had negative raw VPD (Tdew > T2m); floored to 0.",
    n_negative_vpd
  ))
}
vpd <- clamp(vpd_raw, lower = 0, upper = Inf, values = TRUE)

s <- global(vpd, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf(
  "VPD range: %.3f to %.3f kPa (mean: %.3f)",
  min(s$min),
  max(s$max),
  mean(s$mean)
))

# C. Wind speed
log_event("Using monthly mean of daily scalar 2-m wind speed forcing...")
s <- global(u2, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf(
  "Wind speed (2m) range: %.2f to %.2f m/s (mean: %.2f)",
  min(s$min),
  max(s$max),
  mean(s$mean)
))

# D. Net radiation
log_event("Calculating FAO-56 net radiation from revised daily Rs forcing...")
s <- global(Rs, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf(
  "Rs range: %.2f to %.2f MJ/m2/day (mean: %.2f)",
  min(s$min),
  max(s$max),
  mean(s$mean)
))
if (mean(s$mean) < 8 || mean(s$mean) > 20) {
  log_event("WARNING: long-term Rs mean outside broad expected 8-20 MJ/m2/day range.")
}

# Missing-cell fallback is retained only as a last-resort guard. Missing dates
# have already been rejected above. This does not change the monthly time basis.
# The fallback assumes a fixed 10 degC diurnal temperature range (Hargreaves-
# style) rather than the actual Tmax-Tmin, so it is deliberately only a rare
# gap-filler; the count of cells actually using it is logged for audit.
n_missing_rs <- sum(global(is.na(Rs), "sum", na.rm = TRUE)[, 1], na.rm = TRUE)
if (n_missing_rs > 0) {
  log_event(sprintf(
    "WARNING: %d cell-months had missing Rs; filled via Hargreaves-style fallback (K_RS * sqrt(10) * Ra).",
    n_missing_rs
  ))
}
Rs <- if (is.null(TXN)) ifel(is.na(Rs), K_RS * sqrt(10) * Ra_stack, Rs) else
  ifel(is.na(Rs), K_RS * sqrt(terra::cover(terra::clamp(TXN$tmax - TXN$tmin, lower = 0), 10 + 0 * Ra_stack)) * Ra_stack, Rs)
Rso <- (0.75 + 2e-5 * h_m) * Ra_stack
Rns <- (1 - 0.23) * Rs

Rs_Rso <- clamp(Rs / Rso, lower = 0.3, upper = 1.0)
f_cloud <- 1.35 * Rs_Rso - 0.35
# FAO-56 Eq. 39: mean of (Tmax,K)^4 and (Tmin,K)^4 when Tmax/Tmin are available
T4_tmean <- (T_mean + 273.16)^4
T4 <- if (is.null(TXN)) T4_tmean else
  terra::cover(((TXN$tmax + 273.16)^4 + (TXN$tmin + 273.16)^4) / 2, T4_tmean)
Rnl <- 4.903e-9 *
  T4 *
  (0.34 - 0.14 * sqrt(ea)) *
  f_cloud
Rn <- Rns - Rnl

s <- global(Rn, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf(
  "Rn range (before G): %.2f to %.2f MJ/m2/day (mean: %.2f)",
  min(s$min),
  max(s$max),
  mean(s$mean)
))
log_event(sprintf(
  "Layers with negative Rn pixels: %d of %d",
  sum(global(Rn < 0, "sum", na.rm = TRUE)[, 1] > 0),
  nlyr(Rn)
))

# ===============================================================================
# D2. METEOROLOGICAL-CONTROLS BASIN-MEAN EXPORT  [NEW]
# ===============================================================================
# MANUSCRIPT ALIGNMENT (Methods 3f, "Meteorological Controls Associated with
# PET-Formulation Sensitivity"): an external meteorological-controls analysis
# (a companion script, NOT part of this framework) needs basin-mean
# monthly air temperature, net radiation, vapor-pressure deficit, and 2-m wind
# speed to compute standardized calendar-month anomalies
# Z_X,t = (X_t - Xbar_m) / sigma_X,m. All four fields already exist as
# SpatRasters at this point in the script (T_mean, Rn, VPD, u2), so they
# are aggregated to a basin-mean CSV here rather than recomputed later.
#
# Basin mean uses the SAME fractional-coverage area-weighted estimator (basin
# mask = 7-subbasin union, minimum valid coverage 0.20) as the PM / ERA5-Land
# PEV summary CSVs below.
#
# Net radiation is exported BEFORE the soil-heat-flux subtraction (Rn, not
# Rn - G): Rn is the meteorological-control term used by Methods 3f; G is a
# FAO-56 bookkeeping term for the reference-crop surface.
log_event("... Meteorological-controls basin-mean export (Methods 3f)")

met_controls_summary <- data.frame(
  date = format(dates, "%Y-%m"),
  year = as.integer(format(dates, "%Y")),
  month = as.integer(format(dates, "%m")),
  t2m_degC = basin_weighted_mean(terra::values(T_mean, mat = TRUE), PET_BASIN_BW),
  rn_MJ_m2_d = basin_weighted_mean(terra::values(Rn, mat = TRUE), PET_BASIN_BW),
  vpd_kPa = basin_weighted_mean(terra::values(vpd, mat = TRUE), PET_BASIN_BW),
  wind2m_ms = basin_weighted_mean(terra::values(u2, mat = TRUE), PET_BASIN_BW),
  pressure_kPa = basin_weighted_mean(terra::values(P_kpa, mat = TRUE), PET_BASIN_BW)
)

write.csv(
  met_controls_summary,
  file.path(MONTHLY_DIR, "meteorological_controls_monthly_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# D3. GRIDDED PM-PET METEOROLOGICAL-CONTROL EXPORTS
# ------------------------------------------------------------------------------
# These are the gridded monthly fields actually used by the FAO-56 PM
# calculation and are exported so Script 9 (9pr_pet_trends.r) can perform pixel-level trend
# analysis on the physical PM inputs without treating PM-PET itself as a panel.
# The fields are written on the same native ERA5-Land analysis grid and retain
# the monthly 1950-01 to 2025-12 time axis.
#
# PM-input fields: 2-m air temperature, net radiation, VPD, 2-m wind speed,
# and surface pressure.
# ------------------------------------------------------------------------------
log_event("... Writing gridded PM-PET meteorological-control NetCDFs")

control_nc <- list(
  temperature = file.path(MONTHLY_DIR, "pm_control_2m_temperature_monthly.nc"),
  net_radiation = file.path(MONTHLY_DIR, "pm_control_net_radiation_monthly.nc"),
  vpd = file.path(MONTHLY_DIR, "pm_control_vpd_monthly.nc"),
  wind2m = file.path(MONTHLY_DIR, "pm_control_wind2m_monthly.nc"),
  surface_pressure = file.path(MONTHLY_DIR, "pm_control_surface_pressure_monthly.nc")
)

control_write_ok <- c(
  temperature = safe_writeCDF(
    T_mean, control_nc$temperature,
    varname = "t2m",
    longname = "2-m air temperature used by FAO-56 Penman-Monteith PET",
    unit = "degC"
  ),
  net_radiation = safe_writeCDF(
    Rn, control_nc$net_radiation,
    varname = "net_radiation",
    longname = "Net radiation used by FAO-56 Penman-Monteith PET",
    unit = "MJ m-2 d-1"
  ),
  vpd = safe_writeCDF(
    vpd, control_nc$vpd,
    varname = "vpd",
    longname = "Vapour pressure deficit used by FAO-56 Penman-Monteith PET",
    unit = "kPa"
  ),
  wind2m = safe_writeCDF(
    u2, control_nc$wind2m,
    varname = "wind2m",
    longname = "2-m wind speed used by FAO-56 Penman-Monteith PET",
    unit = "m s-1"
  ),
  surface_pressure = safe_writeCDF(
    P_kpa, control_nc$surface_pressure,
    varname = "surface_pressure",
    longname = "Surface pressure used by FAO-56 Penman-Monteith PET",
    unit = "kPa"
  )
)

if (!all(control_write_ok)) {
  stop(
    "Failed to write one or more gridded PM-PET meteorological-control NetCDFs: ",
    paste(names(control_write_ok)[!control_write_ok], collapse = ", ")
  )
}

# Regional meteorological-control export. The 7 source shapefiles are in
# WGS 84 (EPSG:4326); they are transformed to the native CRS of the current raster
# grid before fractional rasterisation.
#
# The canonical grid-aware utility was already sourced during configuration.
# Do not source the legacy utility here: doing so can replace the shared
# regional/grid functions with an incompatible version.
regional_template <- T_mean[[1]]
regional_regions <- load_nechako_regions(template = regional_template)

# Region validation can raise warnings (e.g. boundary/CRS mismatches for the
# 7 source shapefiles). These are captured and written to LOG_FILE instead
# of being left only in the interactive session's warning buffer, since the
# regional CSV below is consumed downstream and any alignment issue should be
# auditable from the log alone.
withCallingHandlers(
  validate_nechako_regions(regional_regions),
  warning = function(w) {
    log_event(paste("REGION VALIDATION WARNING:", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
)

regional_rw <- region_fractional_weights(regional_template, regional_regions)
regional_met_rows <- lapply(names(regional_rw), function(rid) {
  rw <- regional_rw[[rid]]
  data.frame(
    region_id = rid,
    region_name = rw$region_name,
    date = format(dates, "%Y-%m"),
    year = as.integer(format(dates, "%Y")),
    month = as.integer(format(dates, "%m")),
    t2m_degC = regional_weighted_mean(terra::values(T_mean, mat = TRUE), rw$weights,
                                      min_coverage = BASIN_MIN_VALID_AREA_COVERAGE),
    rn_MJ_m2_d = regional_weighted_mean(terra::values(Rn, mat = TRUE), rw$weights,
                                        min_coverage = BASIN_MIN_VALID_AREA_COVERAGE),
    vpd_kPa = regional_weighted_mean(terra::values(vpd, mat = TRUE), rw$weights,
                                     min_coverage = BASIN_MIN_VALID_AREA_COVERAGE),
    wind2m_ms = regional_weighted_mean(terra::values(u2, mat = TRUE), rw$weights,
                                       min_coverage = BASIN_MIN_VALID_AREA_COVERAGE),
    pressure_kPa = regional_weighted_mean(terra::values(P_kpa, mat = TRUE), rw$weights,
                                          min_coverage = BASIN_MIN_VALID_AREA_COVERAGE),
    area_m2 = rw$area_m2,
    n_cells = rw$n_cells,
    n_partial_cells = rw$n_partial,
    stringsAsFactors = FALSE
  )
})
regional_met <- data.table::rbindlist(
  regional_met_rows,
  use.names = TRUE,
  fill = TRUE
)
regional_met_dir <- file.path(MONTHLY_DIR, "regional")
dir.create(regional_met_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(
  regional_met,
  file.path(regional_met_dir, "regional_meteorological_controls_monthly_summary.csv")
)

log_event(sprintf(
  "  Meteorological-controls summary CSV -> %s/meteorological_controls_monthly_summary.csv",
  basename(MONTHLY_DIR)
))
log_event(sprintf(
  "  Basin means -- T: %.2f degC, Rn: %.2f MJ/m2/day, VPD: %.3f kPa, u2: %.2f m/s",
  mean(met_controls_summary$t2m_degC, na.rm = TRUE),
  mean(met_controls_summary$rn_MJ_m2_d, na.rm = TRUE),
  mean(met_controls_summary$vpd_kPa, na.rm = TRUE),
  mean(met_controls_summary$wind2m_ms, na.rm = TRUE)
))

# E. Soil Heat Flux (G) -- FAO-56 Eq.44 with cold-region constraint
log_event("... Soil Heat Flux (G)")
# FAO-56 monthly soil heat flux: Eq.44 gives G = 0.14 x (T_i - T_{i-1}), used
# when the following month's mean temperature (T_{i+1}) is not yet known.
# Coefficient 0.14 is CORRECT here.
# (0.07 belongs to Eq.43, which uses T_{i+1} - T_{i-1} -- a different formula
#  that requires the next month's temperature and is not used in this script.)
if (nlyr(T_mean) > 1) {
  t_prev <- c(T_mean[[1]] * 0, T_mean[[1:(nlyr(T_mean) - 1)]])
  G      <- 0.14 * (T_mean - t_prev)   # FAO-56 Eq.44
  G[[1]] <- 0                           # Jan 1950: no prior month
  log_event("  G = 0.14 x (T_i - T_{i-1})  [FAO-56 Eq.44]")
  log_event("  G[[1]] = 0 (January 1950)")
  
  # COLD-REGION CONSTRAINT
  # Frozen soil / snow cover: thermal conductivity drops ~3x, actual G rarely
  # exceeds 15-25% of Rn (Langer et al. 2011; Katul et al. 2011).
  # Cap |G| at 0.3 x |Rn| when T < 0 degC.
  G_cap <- 0.3 * abs(Rn)
  G     <- ifel(T_mean < 0,
                ifel(G >  G_cap,  G_cap,
                     ifel(G < -G_cap, -G_cap, G)),
                G)
  log_event("  Cold-region constraint applied: |G| <= 0.3 x |Rn| when T < 0 degC")
} else {
  G <- T_mean * 0
}

s <- global(G, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("G range (after constraint): %.2f to %.2f MJ/m2/day (mean: %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))

# F. Final ET0
# Canonical FAO-56 PM helper.  The cold-region G constraint is applied before
# this function is called; the helper itself contains only the PM algebra.
fa056_pm_et0 <- function(T, P_kpa, Rn, G, u2, VPD) {
  es_local <- 0.6108 * exp((17.27 * T) / (T + 237.3))
  delta_local <- (4098 * es_local) / (T + 237.3)^2
  gamma_local <- 0.000665 * P_kpa
  
  num <- 0.408 * delta_local * (Rn - G) +
    gamma_local * (900 / (T + 273)) * u2 * VPD
  den <- delta_local + gamma_local * (1 + 0.34 * u2)
  
  num / den
}

log_event("... Final ET0 calculation using canonical FAO-56 PM helper")
et0 <- fa056_pm_et0(
  T = T_mean,
  P_kpa = P_kpa,
  Rn = Rn,
  G = G,
  u2 = u2,
  VPD = vpd
)

neg_count        <- global(et0 < 0, "sum", na.rm = TRUE)
total_pixels     <- ncell(et0[[1]])
problematic_layers <- sum(neg_count / total_pixels > 0.05)

log_event("=== NEGATIVE ET0 DIAGNOSTIC ===")
for (i in 1:min(12, nlyr(et0)))
  log_event(sprintf("  Layer %d: %.1f%% pixels negative",
                    i, neg_count[i, 1] / total_pixels * 100))
log_event(sprintf("Layers with >5%% negative pixels: %d of %d",
                  problematic_layers, nlyr(et0)))
if (problematic_layers > nlyr(et0) * 0.1) {
  log_event("WARNING: >10% of layers have significant negative values!")
}

et0 <- ifel(et0 < 0, 0, et0)

names(et0) <- sprintf("PET_%04d_%02d",
                      as.numeric(format(dates, "%Y")),
                      as.numeric(format(dates, "%m")))

s <- global(et0, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Final ET0 range: %.3f to %.3f mm/day (mean: %.3f)",
                  min(s$min), max(s$max), mean(s$mean)))

# ==============================================================================
# --- 5b. ERA5-LAND NATIVE POTENTIAL EVAPORATION (PEV) ---
log_event("==========================================")
log_event("ERA5-LAND NATIVE PEV SUMMARY")
log_event("==========================================")

pev_global_stats <- basin_layer_summary(pev_era5)
pet_era5_summary <- data.frame(
  date     = format(dates, "%Y-%m"),
  year     = as.integer(format(dates, "%Y")),
  month    = as.integer(format(dates, "%m")),
  mean_pet = pev_global_stats$mean,
  min_pet  = pev_global_stats$min,
  max_pet  = pev_global_stats$max,
  std_dev  = pev_global_stats$sd
)

# ------------------------------------------------------------------------------
# 5b-QA. PHYSICAL PLAUSIBILITY CHECK: PEV VS. AVAILABLE ENERGY
# ------------------------------------------------------------------------------
# validate_pev_source_metadata() above confirms only that the archive's
# self-reported provenance ATTRIBUTES match the strings the downloader is
# expected to have written. That is a metadata self-consistency check, not
# a verification that the retrieved NUMERIC field is physically correct -- a
# downloader defect could reproduce the expected attributes while still
# writing the wrong data (wrong step, wrong parameter, stale cache, etc.).
# The check below is independent of that metadata and tests the numbers
# directly against an energy-balance bound.
#
# Long-term mean evaporation -- open-water/pan or otherwise -- cannot greatly
# exceed the energy available to vaporize water. Priestley-Taylor bounds
# wet-surface evaporation at alpha * (Rn / lambda), with alpha ~1.26
# empirically; PT_ALPHA_CEILING = 3.0 is a deliberately generous ceiling
# (~2x the empirical PT value) so that only a genuine, large anomaly trips
# it. Rn here is pre-soil-heat-flux (Rn, not Rn - G), consistent with its
# use in the Section D2 meteorological-controls export above.
LAMBDA_MJ_KG     <- 2.45   # FAO-56 latent heat of vaporization, MJ/kg
PT_ALPHA_CEILING <- 3.0    # generous multiple of Rn/lambda; true PT alpha ~1.26

basin_pev_mean     <- mean(pev_global_stats$mean, na.rm = TRUE)
basin_rn_mean      <- mean(basin_weighted_mean(terra::values(Rn, mat = TRUE), PET_BASIN_BW), na.rm = TRUE)
basin_pm_mean      <- mean(basin_weighted_mean(terra::values(et0, mat = TRUE), PET_BASIN_BW), na.rm = TRUE)
basin_energy_bound <- (basin_rn_mean / LAMBDA_MJ_KG) * PT_ALPHA_CEILING
basin_pev_pm_ratio <- basin_pev_mean / basin_pm_mean

log_event(sprintf(
  paste0(
    "PEV PHYSICAL PLAUSIBILITY CHECK: basin-mean ERA5-Land PEV = %.3f mm/day; ",
    "basin-mean Rn = %.3f MJ/m2/day -> energy-limited ceiling (PT alpha=%.1f) ",
    "= %.3f mm/day; basin-mean PM ET0 = %.3f mm/day (PEV:PM ratio = %.2f)."
  ),
  basin_pev_mean, basin_rn_mean, PT_ALPHA_CEILING, basin_energy_bound,
  basin_pm_mean, basin_pev_pm_ratio
))

# IMPORTANT INTERPRETATION:
# ERA5-Land PEV is an open-water/pan potential-evaporation product, whereas
# PM ET0 is FAO-56 reference-crop evapotranspiration. Their magnitudes are NOT
# required to agree, so PEV/PM is a diagnostic only and must not be used as a
# hard rejection criterion.
#
# The previous revision incorrectly stopped whenever PEV/PM differed by >2x or
# PEV exceeded a Priestley-Taylor-style energy ceiling. That conflated two
# different evaporation definitions and made a validly read PEV archive
# impossible to propagate even when its provenance and sign/unit conversion
# were correct.
#
# We retain the checks as NON-FATAL diagnostics. They:
#   1. report the PEV/PM ratio;
#   2. report the generous energy-bound comparison;
#   3. flag the archive for review when either diagnostic is exceeded;
#   4. never rescale PEV to force agreement with PM.
#
# The actual hard safeguards remain the provenance, time-axis, sign, unit,
# non-negativity, geometry, and finite-value checks above. If the upstream
# PEV download is known to be wrong, the correct remedy is to repair the
# downloader/archive rather than alter PEV numerically in this script.

pev_physical_qa_status <- "PASS"
pev_physical_qa_reasons <- character(0)

if (!is.finite(basin_pev_pm_ratio)) {
  pev_physical_qa_status <- "REVIEW"
  pev_physical_qa_reasons <- c(
    pev_physical_qa_reasons,
    "PEV/PM ratio is non-finite"
  )
} else if (basin_pev_pm_ratio > 2 || basin_pev_pm_ratio < 0.5) {
  pev_physical_qa_status <- "REVIEW"
  pev_physical_qa_reasons <- c(
    pev_physical_qa_reasons,
    sprintf("PEV/PM ratio %.2fx is outside the diagnostic 0.5-2.0 range",
            basin_pev_pm_ratio)
  )
}

if (is.finite(basin_pev_mean) &&
    is.finite(basin_energy_bound) &&
    basin_pev_mean > basin_energy_bound) {
  pev_physical_qa_status <- "REVIEW"
  pev_physical_qa_reasons <- c(
    pev_physical_qa_reasons,
    sprintf(
      "PEV %.3f mm/day exceeds the deliberately generous energy diagnostic %.3f mm/day",
      basin_pev_mean, basin_energy_bound
    )
  )
}

if (pev_physical_qa_status == "PASS") {
  log_event(
    "PEV physical QA status: PASS. PEV/PM and energy diagnostics are within the configured review ranges."
  )
} else {
  log_event(
    paste0(
      "WARNING: PEV physical QA status = REVIEW. ",
      paste(pev_physical_qa_reasons, collapse = "; "),
      ". PEV is NOT rescaled or rejected because ERA5-Land PEV and FAO-56 PM ET0 ",
      "are different physical products. A 2026-09-29 spot check (July 1950) matched ",
      "the archive to CDS hourly and standard monthly ERA5-Land PEV, so a high ratio ",
      "alone does not indicate a download or conversion error; re-verify if the ",
      "archive is regenerated."
    )
  )
}

ERA5_PET_NC <- file.path(MONTHLY_DIR, "potential_evaporation_era5land_monthly.nc")
if (!isTRUE(safe_writeCDF(
  pev_era5, ERA5_PET_NC,
  varname = "pet_era5land",
  longname = "ERA5-Land native potential evaporation (PEV; open-water/pan evaporation)",
  unit = "mm/day"
))) {
  stop("Failed to write ERA5-Land PEV NetCDF: ", ERA5_PET_NC)
}

write.csv(pet_era5_summary,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_ERA5Land_monthly_summary.csv"),
          row.names = FALSE)

monthly_stats_era5 <- aggregate(
  pet_era5_summary[, c("mean_pet", "min_pet", "max_pet", "std_dev")],
  by = list(month = pet_era5_summary$month), FUN = mean, na.rm = TRUE)
monthly_stats_era5$month_name <- month.name[monthly_stats_era5$month]
write.csv(monthly_stats_era5,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_ERA5Land_calendar_month_stats.csv"),
          row.names = FALSE)

# --- 6. STATISTICS & OUTPUT (Penman-Monteith) ---
log_event("Generating PM statistics and writing files...")

global_stats <- basin_layer_summary(et0)
pet_summary  <- data.frame(
  date     = format(dates, "%Y-%m"),
  year     = as.numeric(format(dates, "%Y")),
  month    = as.numeric(format(dates, "%m")),
  mean_pet = global_stats$mean,
  min_pet  = global_stats$min,
  max_pet  = global_stats$max,
  std_dev  = global_stats$sd
)

log_event("Writing PM PET NetCDF...")
varnames(et0)  <- "pet"
longnames(et0) <- "Reference Evapotranspiration (FAO-56 Penman-Monteith)"
units(et0)     <- "mm/day"
nc_file        <- file.path(MONTHLY_DIR, "potential_evapotranspiration_monthly.nc")
if (!isTRUE(safe_writeCDF(
  et0, nc_file,
  varname  = "pet",
  longname = sprintf("Reference Evapotranspiration (FAO-56 Penman-Monteith; es and Rnl from %s)",
                     PM_TEMPERATURE_BASIS_USED),
  unit     = "mm/day"
))) {
  stop("Failed to write PM PET NetCDF: ", nc_file)
}

# Written to MONTHLY_DIR so 12basin_timeseries.R (Part 4d) and any external
# companion analysis (not part of this framework) can locate it.
write.csv(pet_summary,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_monthly_summary.csv"),
          row.names = FALSE)
log_event(sprintf("  PM summary CSV -> %s/ERA5Land_Nechako_PET_monthly_summary.csv",
                  basename(MONTHLY_DIR)))

log_event("Calculating calendar month climatology...")
monthly_stats <- aggregate(
  pet_summary[, c("mean_pet", "min_pet", "max_pet", "std_dev")],
  by = list(month = pet_summary$month), FUN = mean, na.rm = TRUE)
monthly_stats$month_name <- month.name[monthly_stats$month]
write.csv(monthly_stats,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_calendar_month_stats.csv"),
          row.names = FALSE)
log_event(sprintf("  PM calendar-month stats CSV -> %s/ERA5Land_Nechako_PET_calendar_month_stats.csv",
                  basename(MONTHLY_DIR)))

# --- 7. DIAGNOSTIC PLOTS ---
log_event("Generating diagnostic plots...")
graphics.off()
diag_dir <- file.path(PET_CALCS_DIR, "diagnostics")
if (!dir.exists(diag_dir)) dir.create(diag_dir, showWarnings = FALSE)

# Average days per calendar month (using 28.25 for Feb to account for leap years)
days_per_month <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# ---------------------------------------------------------------------------
# Helper: render a plot function into the currently open PDF device AND save
# an individual PNG file.
#   name     : filename stem (no extension) for the PNG
#   plot_fn  : zero-argument function containing all plotting calls
#   width/height : inches (shared by PDF page and PNG canvas)
#   res      : PNG resolution in dpi
# ---------------------------------------------------------------------------
save_plot <- function(name, plot_fn, width = 12, height = 8, res = 150) {
  # 1. Draw into the already-open PDF device
  plot_fn()
  # 2. Also write an individual PNG
  png_path <- file.path(diag_dir, paste0(name, ".png"))
  tryCatch({
    png(png_path, width = width, height = height, units = "in", res = res)
    plot_fn()
    dev.off()
    log_event(sprintf("    PNG saved: %s", basename(png_path)))
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    log_event(sprintf("    WARNING: PNG failed for %s -- %s", name, conditionMessage(e)))
  })
}

tryCatch(
  pdf(file.path(diag_dir, "PET_diagnostics.pdf"), width = 12, height = 8),
  error = function(e) stop("Cannot write PDF: close it in your viewer first."))

# A. PM time series
log_event("... A: PM time series")
save_plot("A_PM_timeseries", function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
  plot(dates, pet_summary$mean_pet, type = "l", col = "blue", lwd = 2,
       xlab = "Date", ylab = "PET (mm/day)",
       main = "Basin-Average PET Time Series (Penman-Monteith)",
       cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, font.main = 2)
  grid(lty = "dotted", col = "grey88")
  if (nrow(pet_summary) >= 12) {
    ma_12 <- stats::filter(pet_summary$mean_pet, rep(1/12, 12), sides = 2)
    lines(dates, ma_12, col = "red", lwd = 2)
    legend("topleft", legend = c("Monthly PET", "12-month MA"),
           col = c("blue", "red"), lwd = 2, bty = "n", cex = 1.2)
  }
})

# B. Seasonal cycle
log_event("... B: Seasonal cycle boxplot")
save_plot("B_seasonal_boxplot", function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
  boxplot(mean_pet ~ month, data = pet_summary,
          xlab = "Month", ylab = "PET (mm/day)",
          main = "Seasonal Distribution of PET (PM)",
          names = month.abb, col = "lightblue", border = "darkblue",
          cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, font.main = 2)
  grid(lty = "dotted", col = "grey88")
})

# C. Range validation
log_event("... C: Range validation")
anomalies <- pet_summary[pet_summary$mean_pet < 0 | pet_summary$mean_pet > 8, ]
if (nrow(anomalies) > 0) {
  log_event(paste("WARNING:", nrow(anomalies), "months PET outside 0-8 mm/day"))
  write.csv(anomalies, file.path(diag_dir, "PET_anomalies.csv"), row.names = FALSE)
} else {
  log_event("All PM PET values within expected range (0-8 mm/day)")
}

# D. Monthly climatology (PM only, mm/day)
log_event("... D: Monthly climatology (PM)")
save_plot("D_monthly_climatology_PM", function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
  plot(monthly_stats$month, monthly_stats$mean_pet, type = "b", pch = 19, col = "darkblue",
       xlab = "Month", ylab = "Mean PET (mm/day)",
       main = "Mean PM PET Climatology by Calendar Month",
       xaxt = "n", ylim = c(0, max(monthly_stats$mean_pet, na.rm = TRUE) * 1.1),
       cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, font.main = 2)
  axis(1, at = 1:12, labels = month.abb, cex.axis = 1.2)
  arrows(monthly_stats$month,
         monthly_stats$mean_pet - monthly_stats$std_dev,
         monthly_stats$month,
         monthly_stats$mean_pet + monthly_stats$std_dev,
         angle = 90, code = 3, length = 0.05, col = "gray50")
  grid(lty = "dotted", col = "grey88")
})

# E. Inter-annual variability with Sen's slope
log_event("... E: Inter-annual variability")
annual_mean <- aggregate(mean_pet ~ year, data = pet_summary, FUN = mean, na.rm = TRUE)

# Trend statistics (linear regression, Theil-Sen slope, Mann-Kendall test) are
# computed ONCE here, outside the plot closure. save_plot() draws plot_fn()
# twice -- once to the open PDF device and once to a PNG device -- so any
# computation or log_event() call placed inside plot_fn was previously
# executed and logged twice (harmless here since the arithmetic is
# deterministic, but wasteful and confusing in the log). The plot closure
# below only re-draws using these precomputed values; it does not recompute
# or re-log anything.
has_trend <- nrow(annual_mean) > 2
if (has_trend) {
  trend_lm   <- lm(mean_pet ~ year, data = annual_mean)
  n          <- nrow(annual_mean)
  slopes_mat <- outer(seq_len(n), seq_len(n), function(i, j)
    ifelse(j > i,
           (annual_mean$mean_pet[j] - annual_mean$mean_pet[i]) /
             (annual_mean$year[j]   - annual_mean$year[i]), NA))
  sen_slope     <- median(slopes_mat, na.rm = TRUE)
  sen_intercept <- median(annual_mean$mean_pet - sen_slope * annual_mean$year, na.rm = TRUE)
  S     <- sum(sign(outer(annual_mean$mean_pet, annual_mean$mean_pet, "-")[lower.tri(matrix(0, n, n))]))
  var_S <- n * (n - 1) * (2 * n + 5) / 18
  Z     <- ifelse(S > 0, (S - 1) / sqrt(var_S), ifelse(S < 0, (S + 1) / sqrt(var_S), 0))
  mk_p  <- 2 * pnorm(-abs(Z))
  log_event(sprintf("Sen's slope: %.6f mm/day/year (MK p=%.3f)", sen_slope, mk_p))
}

save_plot("E_interannual_trend_PM", function() {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  plot(annual_mean$year, annual_mean$mean_pet, type = "b", pch = 19, col = "darkgreen",
       xlab = "Year", ylab = "Annual Mean PET (mm/day)", main = "Annual Mean PET Trend (PM)")
  grid()
  if (has_trend) {
    abline(trend_lm, col = "red", lwd = 2, lty = 2)
    abline(a = sen_intercept, b = sen_slope, col = "blue", lwd = 2)
    legend("topleft",
           legend = c(sprintf("LR: %.4f mm/day/yr (p=%.3f)",
                              coef(trend_lm)[2], summary(trend_lm)$coefficients[2, 4]),
                      sprintf("Sen: %.4f mm/day/yr (MK p=%.3f)", sen_slope, mk_p)),
           col = c("red", "blue"), lwd = 2, lty = c(2, 1), bty = "n")
  }
})

# F. PM sample spatial maps
log_event("... F: PM spatial maps")
sample_months <- c(1, 256, 619, 910)
save_plot("F_PM_spatial_maps", function() {
  if (nlyr(et0) >= 12) {
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))
    for (idx in sample_months)
      plot(et0[[idx]], col = rev(hcl.colors(100, "YlOrRd")),
           main = sprintf("PET_PM %s", names(et0)[idx]))
  }
})

# F2. ERA5-Land PEV vs. Penman-Monteith comparison
log_event("... F2: PM vs. ERA5-Land native PEV comparison")

save_plot("F2a_timeseries_PM_vs_ERA5Land", function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
  y_lim <- range(c(pet_summary$mean_pet, pet_era5_summary$mean_pet), na.rm = TRUE)
  plot(dates, pet_summary$mean_pet, type = "l", col = "blue", lwd = 2,
       ylim = y_lim, xlab = "Date", ylab = "PET / PEV (mm/day)",
       main = "Basin-Average PET: FAO-56 PM vs. ERA5-Land PEV",
       cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, font.main = 2)
  lines(dates, pet_era5_summary$mean_pet, col = "darkorange", lwd = 2, lty = 2)
  legend("topright", legend = c("FAO-56 Penman-Monteith", "ERA5-Land PEV (pan evaporation)"),
         col = c("blue", "darkorange"), lwd = 2, lty = c(1, 2), bty = "n", cex = 1.05)
  grid(lty = "dotted", col = "grey88")
})

pet_clim_pm <- monthly_stats$mean_pet * days_per_month
pet_clim_era5 <- monthly_stats_era5$mean_pet * days_per_month
sd_pm <- monthly_stats$std_dev * days_per_month
sd_era5 <- monthly_stats_era5$std_dev * days_per_month
save_plot("F2b_seasonal_climatology_PM_vs_ERA5Land", function() {
  par(mfrow = c(1, 1), mar = c(5, 6, 4, 2))
  y_max <- max(c(pet_clim_pm, pet_clim_era5), na.rm = TRUE) * 1.15
  plot(1:12, pet_clim_pm, type = "b", pch = 19, col = "blue", lwd = 2.5,
       xlab = "Month", ylab = "Mean PET / PEV (mm/month)",
       main = "Seasonal Climatology: FAO-56 PM vs. ERA5-Land PEV",
       xaxt = "n", ylim = c(0, y_max), cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.3, font.main = 2)
  axis(1, at = 1:12, labels = month.abb, cex.axis = 1.3)
  lines(1:12, pet_clim_era5, type = "b", pch = 17, col = "darkorange", lwd = 2.5, lty = 2)
  suppressWarnings({
    arrows(1:12, pet_clim_pm - sd_pm, 1:12, pet_clim_pm + sd_pm, angle = 90, code = 3, length = 0.05, col = "blue")
    arrows(1:12, pet_clim_era5 - sd_era5, 1:12, pet_clim_era5 + sd_era5, angle = 90, code = 3, length = 0.05, col = "darkorange")
  })
  legend("topleft", legend = c("FAO-56 Penman-Monteith", "ERA5-Land PEV (pan evaporation)"),
         col = c("blue", "darkorange"), lwd = 2.5, pch = c(19, 17), lty = c(1, 2), bty = "n", cex = 1.0)
  grid(lty = "dotted", col = "grey88")
})

delta_pet_monthly <- (monthly_stats$mean_pet - monthly_stats_era5$mean_pet) * days_per_month
save_plot("F2c_monthly_difference_PM_minus_ERA5Land", function() {
  par(mfrow = c(1, 1), mar = c(4, 5, 3, 1))
  barplot(delta_pet_monthly, names.arg = month.abb,
          col = ifelse(delta_pet_monthly > 0, "#e31a1c", "#1f78b4"), border = "white",
          xlab = "Month", ylab = "PM minus ERA5-Land PEV (mm/month)",
          main = "Monthly PET-source difference")
  abline(h = 0, lwd = 1.5)
  grid(nx = NA, ny = NULL)
})

annual_pm_era5 <- merge(
  aggregate(mean_pet ~ year, data = pet_summary, FUN = mean, na.rm = TRUE),
  aggregate(mean_pet ~ year, data = pet_era5_summary, FUN = mean, na.rm = TRUE),
  by = "year", suffixes = c("_pm", "_era5land"))
r_pm_era5 <- cor(annual_pm_era5$mean_pet_pm, annual_pm_era5$mean_pet_era5land, use = "complete.obs")
save_plot("F2d_annual_scatter_PM_vs_ERA5Land", function() {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  plot(annual_pm_era5$mean_pet_era5land, annual_pm_era5$mean_pet_pm,
       pch = 19, col = "#4393c3", cex = 0.9,
       xlab = "ERA5-Land PEV (mm/day)", ylab = "FAO-56 PM PET (mm/day)",
       main = sprintf("Annual Mean PET/PEV: PM vs. ERA5-Land (r = %.3f)", r_pm_era5))
  abline(0, 1, col = "black", lty = 2, lwd = 1.5)
  abline(lm(mean_pet_pm ~ mean_pet_era5land, data = annual_pm_era5), col = "red", lwd = 2)
  grid()
})

# Fitted once, outside the plot closure, for the same reason as panel E above
# (save_plot() draws plot_fn() twice, once per output device).
lm_pm <- lm(mean_pet_pm ~ year, data = annual_pm_era5)
lm_e5 <- lm(mean_pet_era5land ~ year, data = annual_pm_era5)

save_plot("F2e_annual_trends_PM_vs_ERA5Land", function() {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  y_r <- range(c(annual_pm_era5$mean_pet_pm, annual_pm_era5$mean_pet_era5land), na.rm = TRUE)
  plot(annual_pm_era5$year, annual_pm_era5$mean_pet_pm, type = "b", pch = 19, col = "blue", lwd = 1.5,
       ylim = y_r, xlab = "Year", ylab = "Annual Mean PET / PEV (mm/day)",
       main = "Annual PET/PEV Trends: PM vs. ERA5-Land")
  lines(annual_pm_era5$year, annual_pm_era5$mean_pet_era5land, type = "b", pch = 17, col = "darkorange", lwd = 1.5)
  abline(lm_pm, col = "blue", lty = 2, lwd = 2)
  abline(lm_e5, col = "darkorange", lty = 2, lwd = 2)
  legend("topleft", legend = c(sprintf("PM: %+.4f mm/day/yr (p=%.3f)", coef(lm_pm)[2], summary(lm_pm)$coefficients[2,4]),
                               sprintf("ERA5-Land PEV: %+.4f mm/day/yr (p=%.3f)", coef(lm_e5)[2], summary(lm_e5)$coefficients[2,4])),
         col = c("blue", "darkorange"), pch = c(19, 17), lty = 1, lwd = 2, bty = "n")
  grid()
})

save_plot("F2f_ERA5Land_PEV_spatial_maps", function() {
  if (nlyr(pev_era5) >= 12) {
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))
    for (idx in sample_months) if (idx <= nlyr(pev_era5))
      plot(pev_era5[[idx]], col = rev(hcl.colors(100, "YlOrRd")), main = sprintf("PET_ERA5Land %s", names(pev_era5)[idx]))
  }
})

# G. Distribution histogram (PM)
log_event("... G: Distribution histogram")
save_plot("G_distribution_histogram_PM", function() {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  hist(pet_summary$mean_pet, breaks = 30, col = "skyblue", border = "white",
       xlab = "PET (mm/day)", main = "Distribution of Basin-Mean PET (PM)", freq = FALSE)
  lines(density(pet_summary$mean_pet, na.rm = TRUE), col = "red", lwd = 2)
  curve(dnorm(x, mean(pet_summary$mean_pet, na.rm = TRUE),
              sd(pet_summary$mean_pet,   na.rm = TRUE)),
        add = TRUE, col = "blue", lwd = 2, lty = 2)
  legend("topright", legend = c("Kernel density", "Normal fit"),
         col = c("red", "blue"), lwd = 2, lty = c(1, 2), bty = "n")
})

# H. QC summary table
log_event("... H: QC summary table")
qc <- data.frame(
  Metric = c("Total Months",
             "PM Mean PET (mm/day)", "PM Min", "PM Max", "PM Std Dev",
             "ERA5-Land PEV Mean (mm/day)", "ERA5-Land PEV Max",
             "PM: % months PET > 8",
             "G coefficient (FAO-56 Eq.44)",
             "Cold-region G constraint",
             "Rs conversion"),
  Value = c(
    nrow(pet_summary),
    round(mean(pet_summary$mean_pet,     na.rm = TRUE), 3),
    round(min(pet_summary$min_pet,       na.rm = TRUE), 3),
    round(max(pet_summary$max_pet,       na.rm = TRUE), 3),
    round(mean(pet_summary$std_dev,      na.rm = TRUE), 3),
    round(mean(pet_era5_summary$mean_pet, na.rm = TRUE), 3),
    round(max(pet_era5_summary$max_pet,   na.rm = TRUE), 3),
    round(sum(pet_summary$mean_pet > 8,  na.rm = TRUE) / nrow(pet_summary) * 100, 2),
    "0.14 x (T_i - T_{i-1})",
    "|G| <= 0.3 x |Rn| when T < 0 C",
    "monthly mean of Rs_daily_YYYY.nc [MJ/m2/day]"
  )
)
save_plot("H_QC_summary_table", function() {
  par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
  title("Quality Control Summary", cex.main = 1.4)
  y0 <- 0.93; ys <- 0.075
  for (i in seq_len(nrow(qc))) {
    text(0.04, y0 - (i - 1) * ys, qc$Metric[i], adj = 0, cex = 0.9, font = 2)
    text(0.58, y0 - (i - 1) * ys, qc$Value[i],  adj = 0, cex = 0.9)
  }
  lines(c(0.04, 0.96), c(y0 + 0.02, y0 + 0.02), lwd = 2)
})

dev.off()
log_event("Diagnostic PDF saved.")
log_event(sprintf("Individual PNGs saved in: %s/", diag_dir))

# ==============================================================================
# --- 5c. COUNTERFACTUAL / DETRENDED PET BRANCH REMOVED ---
# The manuscript workflow now uses observed FAO-56 PM PET and observed
# ERA5-Land native PEV. No detrended PET product is generated in this script.

# --- 8. FINISH ---
log_event("==========================================")
log_event(paste("Total months processed:", n_months))
log_event(paste("Penman-Monteith mean PET:",
                round(mean(pet_summary$mean_pet,     na.rm = TRUE), 2), "mm/day"))
log_event(paste("ERA5-Land PEV mean:       ",
                round(mean(pet_era5_summary$mean_pet, na.rm = TRUE), 2), "mm/day"))
log_event(paste("PM layers with >5% negative pixels:", problematic_layers))
log_event("Outputs written to MONTHLY_DIR (monthly_data_direct/):")
log_event("  potential_evapotranspiration_monthly.nc")
log_event("  potential_evaporation_era5land_monthly.nc")
log_event("  ERA5Land_Nechako_PET_monthly_summary.csv")
log_event("  ERA5Land_Nechako_PET_calendar_month_stats.csv")
log_event("  ERA5Land_Nechako_PET_ERA5Land_monthly_summary.csv")
log_event("  ERA5Land_Nechako_PET_ERA5Land_calendar_month_stats.csv")
log_event("  [detrended/counterfactual PET branch removed from production workflow]")
log_event("Diagnostics written to PET_CALCS_DIR (pet_calcs/):")
log_event("  diagnostics/PET_diagnostics.pdf  + individual PNGs")
log_event("  [no Thornthwaite auxiliary rasters]")
log_event("  monthly_pet_processing.log")
log_event("PROCESSING COMPLETE")
log_event(paste("Finished:", Sys.time()))