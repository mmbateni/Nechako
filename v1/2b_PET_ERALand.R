##############################################
# SCRIPT 2b ??? PET CALCULATION
# FAO-56 Penman-Monteith PET + ERA5-Land native potential evaporation (PEV)
#
# EXECUTION ORDER: Run after Monthly_ERALand_download_v2.R.
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
#   3. ERA5-Land native potential_evaporation is loaded and converted from m/day
#      to mm/day as the second observed PET/PEV product.
#   4. PM-vs-ERA5-Land diagnostics replace all PM-vs-Thornthwaite diagnostics.
##############################################
library(terra)
library(ncdf4)
library(lubridate)

# --- 1. CONFIGURATION ---
# Working directory remains: D:/Nechako_Drought/Nechako/
# (do NOT call setwd ??? caller controls the working directory)

MONTHLY_DIR   <- "D:/Nechako_Drought/Nechako/monthly_data_direct"  # monthly ERA5-Land state variables
DAILY_DIR     <- "D:/Nechako_Drought/Nechako/daily_eto_forcing"    # revised daily u2 + Rs archives
PET_CALCS_DIR <- "D:/Nechako_Drought/Nechako/pet_calcs"            # diagnostics, aux rasters, log

# OUTPUT LOCATION POLICY
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# NetCDF rasters  ??? MONTHLY_DIR  (consumed by 3SPEI_ERALand.R, 10a, 10b)
# Summary CSVs    ??? MONTHLY_DIR  (consumed by 7basin_timeseries.R Part 4d
#                                 and 10b_PET_bias_nonstationarity.R)
# Auxiliary products -> PET_CALCS_DIR
# Diagnostic PDFs ??? PET_CALCS_DIR/diagnostics/
# Processing log  ??? PET_CALCS_DIR/
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????


# Create output directories if they do not yet exist
if (!dir.exists(MONTHLY_DIR))   dir.create(MONTHLY_DIR,   recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(DAILY_DIR))     stop("DAILY_DIR does not exist: ", DAILY_DIR, "\nRun Monthly_ERALand_download_v2.R first.")
if (!dir.exists(PET_CALCS_DIR)) dir.create(PET_CALCS_DIR, recursive = TRUE, showWarnings = FALSE)

G_GRAV   <- 9.80665
K_RS     <- 0.16   # Interior coefficient for Nechako Basin
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
      log_event(sprintf("  WARNING: Cannot delete locked file: %s", basename(filepath)))
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
    if (isTRUE(ok)) return(invisible(TRUE))
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
  t2m_b <- rast(file.path(MONTHLY_DIR, "2m_temperature_monthly.nc"))
  d2m_b <- rast(file.path(MONTHLY_DIR, "2m_dewpoint_temperature_monthly.nc"))
  sp_b  <- rast(file.path(MONTHLY_DIR, "surface_pressure_monthly.nc"))
  log_event("Successfully loaded monthly temperature, dewpoint, and pressure inputs.")
}, error = function(e) stop("Input file loading failed: ", e$message))

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
avg_elevation <- global(h_m, "mean", na.rm = TRUE)[1, 1]
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
# ERA5-Land native potential evaporation (PEV)
# ------------------------------------------------------------------------------
# ERA5-Land provides `potential_evaporation` (pev) as a ready-to-use
# potential-evaporation product. It is scientifically distinct from the
# independently calculated FAO-56 Penman-Monteith reference ET below: ECMWF
# documents ERA5-Land PEV as an open-water/pan-evaporation quantity. We retain
# that distinction explicitly in filenames, labels, and downstream analyses.
PEV_FILE <- file.path(MONTHLY_DIR, "potential_evaporation_monthly.nc")
if (!file.exists(PEV_FILE)) {
  stop("ERA5-Land native potential evaporation file not found: ", PEV_FILE,
       "\nRun Monthly_ERALand_download_v2.R first.")
}

log_event("Loading ERA5-Land native potential evaporation (PEV)...")
pev_raw <- rast(PEV_FILE)
pev_dates <- extract_time_robust(pev_raw, PEV_FILE, "potential_evaporation")
pev_dates <- as.Date(format(pev_dates, "%Y-%m-01"))
pev_idx <- match(dates, pev_dates)
if (anyNA(pev_idx)) {
  stop("ERA5-Land PEV does not cover the complete PET analysis period. Missing: ",
       paste(format(dates[is.na(pev_idx)], "%Y-%m"), collapse = ", "))
}
pev_raw <- pev_raw[[pev_idx]]
if (!compareGeom(pev_raw[[1]], t2m_b[[1]], stopOnError = FALSE)) {
  log_event("Resampling ERA5-Land PEV to the temperature grid...")
  pev_raw <- project(pev_raw, t2m_b[[1]], method = "bilinear")
}
# Monthly averaged ERA5-Land hydrological accumulations are expressed as
# metres per day in the monthly-mean product; convert to mm/day.
pev_era5 <- pev_raw * 1000
# PEV is a non-negative evaporative-demand quantity; negative values indicate
# malformed/incorrectly decoded input and are not silently converted to positive.
if (global(pev_era5 < 0, "sum", na.rm = TRUE)[1,1] > 0)
  stop("Negative values detected in ERA5-Land PEV after unit conversion; inspect the downloaded NetCDF.")
pev_era5 <- set_monthly_layer_names(pev_era5, dates, "PET_ERA5LAND")
pev_stats <- global(pev_era5, c("min", "max", "mean", "sd"), na.rm = TRUE)
log_event(sprintf("ERA5-Land PEV: %.3f to %.3f mm/day (mean: %.3f)",
                  min(pev_stats$min), max(pev_stats$max), mean(pev_stats$mean)))


# --- 3.25 REVISED DAILY FORCING -> MONTHLY MEANS ---
# Monthly PM still operates on one layer per calendar month, but wind and Rs are
# now derived from physically correct daily products.  Wind magnitude was formed
# from hourly u/v before daily averaging in Monthly_ERALand_download_modified.R;
# here we only average those daily scalar values within each calendar month.
aggregate_daily_archive_to_monthly <- function(variable_name,
                                               expected_dates,
                                               template) {
  
  years_needed <- sort(unique(as.integer(format(expected_dates, "%Y"))))
  monthly_parts <- vector("list", length(years_needed))
  monthly_dates <- vector("list", length(years_needed))
  
  for (j in seq_along(years_needed)) {
    
    yr <- years_needed[j]
    
    f <- file.path(
      DAILY_DIR,
      sprintf("%s_daily_%04d.nc", variable_name, yr)
    )
    
    if (!file.exists(f))
      stop("Missing revised daily forcing file: ", f)
    
    x <- terra::rast(f)
    
    tt <- terra::time(x)
    
    if (!(inherits(tt, "Date") || inherits(tt, "POSIXt")) ||
        all(is.na(tt))) {
      stop("No valid NetCDF time axis in: ", f)
    }
    
    dd <- as.Date(tt)
    
    keep <- which(as.integer(format(dd, "%Y")) == yr)
    
    if (!length(keep))
      stop("No dates for year ", yr, " in ", f)
    
    x  <- x[[keep]]
    dd <- dd[keep]
    
    # ------------------------------------------------------------
    # Daily completeness check
    # ------------------------------------------------------------
    
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
    
    # ------------------------------------------------------------
    # Aggregate daily -> monthly
    # ------------------------------------------------------------
    
    month_index <- as.integer(format(dd, "%m"))
    
    xm <- terra::tapp(
      x,
      index = month_index,
      fun = mean,
      na.rm = TRUE
    )
    
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
  
  # ------------------------------------------------------------
  # Combine annual pieces
  # ------------------------------------------------------------
  
  out <- do.call(c, monthly_parts)
  
  od <- do.call(c, monthly_dates)
  
  idx <- match(expected_dates, od)
  
  if (anyNA(idx)) {
    stop(
      variable_name,
      " monthly aggregation is missing: ",
      paste(
        format(expected_dates[is.na(idx)], "%Y-%m"),
        collapse = ", "
      )
    )
  }
  
  out <- out[[idx]]
  
  # ------------------------------------------------------------
  # Geometry check
  # ------------------------------------------------------------
  
  if (!terra::compareGeom(
    out[[1]],
    template[[1]],
    stopOnError = FALSE
  )) {
    
    log_event(sprintf(
      "  Resampling %s monthly forcing to temperature grid.",
      variable_name
    ))
    
    out <- terra::resample(
      out,
      template[[1]],
      method = "bilinear"
    )
  }
  
  out
}
log_event("Aggregating revised daily u2 and Rs archives to monthly means...")
u2 <- aggregate_daily_archive_to_monthly("u2_mean", dates, t2m_b)
Rs <- aggregate_daily_archive_to_monthly("Rs",      dates, t2m_b)
names(u2) <- sprintf("u2_%04d_%02d", as.integer(format(dates, "%Y")), as.integer(format(dates, "%m")))
names(Rs) <- sprintf("Rs_%04d_%02d", as.integer(format(dates, "%Y")), as.integer(format(dates, "%m")))

# --- 3.5 DAYS IN MONTH ---
log_event("Calculating days in each month...")
first_of_month   <- as.Date(format(dates, "%Y-%m-01"))
first_next_month <- seq(first_of_month[1], by = "month", length.out = n_months + 1)[-1]
days_in_month    <- as.integer(first_next_month - first_of_month)
log_event(sprintf("Days range: %d to %d", min(days_in_month), max(days_in_month)))

# --- 4. PRE-CALCULATIONS ---

# A. Extraterrestrial Radiation (Ra) - vectorised
log_event("Calculating Ra (Extraterrestrial Radiation)...")
lat_r   <- init(t2m_b, "y")
lat_rad <- lat_r * pi / 180

mid_month_days <- c(16, 46, 75, 105, 135, 166, 196, 227, 258, 288, 319, 349)
ra_list <- list()

for (m in 1:12) {
  J        <- mid_month_days[m]
  delta    <- 0.409 * sin(2 * pi * J / 365 - 1.39)
  dr       <- 1 + 0.033 * cos(2 * pi * J / 365)
  tan_prod <- ifel(tan(lat_rad) * tan(delta) >  1,  1,
                   ifel(tan(lat_rad) * tan(delta) < -1, -1,
                        tan(lat_rad) * tan(delta)))
  ws       <- acos(-tan_prod)
  Ra_val   <- (24 * 60 / pi) * 0.0820 * dr *
    (ws * sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * sin(ws))
  ra_list[[m]] <- ifel(Ra_val < 0, 0, Ra_val)
}
ra_clim  <- rast(ra_list)
Ra_stack <- ra_clim[[month_indices]]

log_event("Ra climatology statistics:")
for (m in 1:12) {
  s <- global(ra_clim[[m]], c("min", "max", "mean"), na.rm = TRUE)
  log_event(sprintf("  Month %2d: min=%.2f, max=%.2f, mean=%.2f MJ/m2/day",
                    m, s$min, s$max, s$mean))
}
log_event("Ra stack generated.")

# --- 5. PENMAN-MONTEITH PET CALCULATION ---
log_event("Executing vectorized Penman-Monteith calculation...")

# A. Temperature & Pressure
log_event("... Temperature and Pressure")
T_mean <- t2m_b - 273.15
P_kpa  <- sp_b  / 1000
gamma  <- 0.000665 * P_kpa

s <- global(T_mean, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Temperature range: %.2f to %.2f degC (mean: %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))

# B. Vapor Pressure
log_event("... Vapor Pressure")
es   <- 0.6108 * exp((17.27 * T_mean) / (T_mean + 237.3))
Tdew <- ifel(is.na(d2m_b - 273.15), T_mean - 2, d2m_b - 273.15)
ea   <- 0.6108 * exp((17.27 * Tdew) / (Tdew + 237.3))
s    <- global(es - ea, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("VPD range: %.3f to %.3f kPa (mean: %.3f)",
                  min(s$min), max(s$max), mean(s$mean)))

# C. Wind Speed
log_event("... Wind Speed from revised daily u2 archive")
s <- global(u2, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Wind speed (2m; monthly mean of daily scalar wind) range: %.2f to %.2f m/s (mean: %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))

# D. Solar Radiation
log_event("... Net Radiation from revised daily Rs archive")
# Rs is already MJ m-2 day-1 in Rs_daily_YYYY.nc.  Monthly aggregation above
# computes the mean daily solar energy for each calendar month; no unit conversion
# and no division by days_in_month is required here.
s <- global(Rs, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  Rs range: %.2f to %.2f MJ/m2/day (mean across layers: %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))
if (mean(s$mean) < 8 || mean(s$mean) > 20) {
  log_event("  WARNING: long-term Rs mean outside broad expected 8-20 MJ/m2/day range.")
}

# Missing-cell fallback is retained only as a last-resort guard.  It does not
# replace missing dates, which are rejected during archive validation above.
Rs <- ifel(is.na(Rs), K_RS * sqrt(10) * Ra_stack, Rs)
Rso <- (0.75 + 2e-5 * h_m) * Ra_stack
Rns <- (1 - 0.23) * Rs

Rs_Rso  <- clamp(Rs / Rso, lower = 0.3, upper = 1.0)
f_cloud <- 1.35 * Rs_Rso - 0.35
Rnl     <- 4.903e-9 * (T_mean + 273.15)^4 * (0.34 - 0.14 * sqrt(ea)) * f_cloud
Rn      <- Rns - Rnl

s <- global(Rn, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Rn range (before G): %.2f to %.2f MJ/m2/day (mean: %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))
log_event(sprintf("Layers with negative Rn pixels: %d of %d",
                  sum(global(Rn < 0, "sum", na.rm = TRUE) > 0), nlyr(Rn)))

# ==============================================================================
# D2. METEOROLOGICAL-CONTROLS BASIN-MEAN EXPORT  [NEW]
# ==============================================================================
# MANUSCRIPT ALIGNMENT (Methods 3f, "Meteorological Controls Associated with
# PET-Formulation Sensitivity"): 13meteorological_controls.R needs basin-mean
# monthly air temperature, net radiation, vapor-pressure deficit, and 2-m wind
# speed to compute standardized calendar-month anomalies
# Z_X,t = (X_t - Xbar_m) / sigma_X,m. All four fields already exist as
# SpatRasters at this point in the script (T_mean, Rn, es - ea, u2), so they
# are aggregated to a basin-mean CSV here rather than recomputed later.
#
# Basin mean uses the SAME unweighted global() mean already used for the PM/
# ERA5-Land PEV summary CSVs below, so this file is
# directly comparable to ERA5Land_Nechako_PET_monthly_summary.csv. This is a
# grid of ~0.1 deg cells over a single basin; 3SPEI_ERALand.R documents the
# same unweighted-vs-area-weighted distinction for WB_basin_average_monthly.csv
# (area-weighting is reserved for standardized-index basin averages computed
# from equal-area-projected pixel grids in 6trend_test_ALL.R).
#
# Net radiation is exported BEFORE the soil-heat-flux subtraction (Rn, not
# Rn - G): Rn is the term Methods 3f actually names ("net radiation"); G is a
# FAO-56 bookkeeping term for the reference-crop surface, not a meteorological
# control in its own right.
log_event("... Meteorological-controls basin-mean export (Methods 3f)")

vpd <- es - ea

met_controls_summary <- data.frame(
  date        = format(dates, "%Y-%m"),
  year        = as.numeric(format(dates, "%Y")),
  month       = as.numeric(format(dates, "%m")),
  t2m_degC    = global(T_mean, "mean", na.rm = TRUE)[, 1],
  rn_MJ_m2_d  = global(Rn,     "mean", na.rm = TRUE)[, 1],
  vpd_kPa     = global(vpd,    "mean", na.rm = TRUE)[, 1],
  wind2m_ms   = global(u2,     "mean", na.rm = TRUE)[, 1],
  pressure_kPa = global(P_kpa, "mean", na.rm = TRUE)[, 1]
)

write.csv(met_controls_summary,
          file.path(MONTHLY_DIR, "meteorological_controls_monthly_summary.csv"),
          row.names = FALSE)
log_event(sprintf("  Meteorological-controls summary CSV -> %s/meteorological_controls_monthly_summary.csv",
                  basename(MONTHLY_DIR)))
log_event(sprintf("  Basin means -- T: %.2f degC, Rn: %.2f MJ/m2/day, VPD: %.3f kPa, u2: %.2f m/s",
                  mean(met_controls_summary$t2m_degC, na.rm = TRUE),
                  mean(met_controls_summary$rn_MJ_m2_d, na.rm = TRUE),
                  mean(met_controls_summary$vpd_kPa, na.rm = TRUE),
                  mean(met_controls_summary$wind2m_ms, na.rm = TRUE)))

# E. Soil Heat Flux (G) -- FAO-56 Eq.45 with cold-region constraint
log_event("... Soil Heat Flux (G)")
# Eq.45: G = 0.14 x (T_i - T_{i-1}).  Coefficient 0.14 is CORRECT here.
# (0.07 belongs to Eq.44 which uses T_{i+1} - T_{i-1} -- a different formula.)
if (nlyr(T_mean) > 1) {
  t_prev <- c(T_mean[[1]] * 0, T_mean[[1:(nlyr(T_mean) - 1)]])
  G      <- 0.14 * (T_mean - t_prev)   # FAO-56 Eq.45
  G[[1]] <- 0                           # Jan 1950: no prior month
  log_event("  G = 0.14 x (T_i - T_{i-1})  [FAO-56 Eq.45]")
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
log_event("... Final ET0 calculation")
delta_slope <- (4098 * es) / (T_mean + 237.3)^2
num         <- 0.408 * delta_slope * (Rn - G) +
  gamma * (900 / (T_mean + 273)) * u2 * (es - ea)
den         <- delta_slope + gamma * (1 + 0.34 * u2)
et0         <- num / den

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

pev_global_stats <- global(pev_era5, c("mean", "min", "max", "sd"), na.rm = TRUE)
pet_era5_summary <- data.frame(
  date     = format(dates, "%Y-%m"),
  year     = as.integer(format(dates, "%Y")),
  month    = as.integer(format(dates, "%m")),
  mean_pet = pev_global_stats$mean,
  min_pet  = pev_global_stats$min,
  max_pet  = pev_global_stats$max,
  std_dev  = pev_global_stats$sd
)

ERA5_PET_NC <- file.path(MONTHLY_DIR, "potential_evaporation_era5land_monthly.nc")
safe_writeCDF(pev_era5, ERA5_PET_NC,
              varname = "pet_era5land",
              longname = "ERA5-Land native potential evaporation (PEV; open-water/pan evaporation)",
              unit = "mm/day")

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

global_stats <- global(et0, c("mean", "min", "max", "sd"), na.rm = TRUE)
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
safe_writeCDF(et0, nc_file,
              varname  = "pet",
              longname = "Reference Evapotranspiration (FAO-56 Penman-Monteith)",
              unit     = "mm/day")

# Written to MONTHLY_DIR so 7basin_timeseries.R (Part 4d) and
# 10b_PET_bias_nonstationarity.R can locate it without path changes.
write.csv(pet_summary,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_monthly_summary.csv"),
          row.names = FALSE)
log_event(sprintf("  PM summary CSV ??? %s/ERA5Land_Nechako_PET_monthly_summary.csv",
                  basename(MONTHLY_DIR)))

log_event("Calculating calendar month climatology...")
monthly_stats <- aggregate(
  pet_summary[, c("mean_pet", "min_pet", "max_pet", "std_dev")],
  by = list(month = pet_summary$month), FUN = mean, na.rm = TRUE)
monthly_stats$month_name <- month.name[monthly_stats$month]
write.csv(monthly_stats,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_calendar_month_stats.csv"),
          row.names = FALSE)
log_event(sprintf("  PM calendar-month stats CSV ??? %s/ERA5Land_Nechako_PET_calendar_month_stats.csv",
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
    log_event(sprintf("    WARNING: PNG failed for %s ??? %s", name, conditionMessage(e)))
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
save_plot("E_interannual_trend_PM", function() {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  plot(annual_mean$year, annual_mean$mean_pet, type = "b", pch = 19, col = "darkgreen",
       xlab = "Year", ylab = "Annual Mean PET (mm/day)", main = "Annual Mean PET Trend (PM)")
  grid()
  if (nrow(annual_mean) > 2) {
    trend_lm   <- lm(mean_pet ~ year, data = annual_mean)
    abline(trend_lm, col = "red", lwd = 2, lty = 2)
    n <- nrow(annual_mean)
    slopes_mat <- outer(seq_len(n), seq_len(n), function(i, j)
      ifelse(j > i,
             (annual_mean$mean_pet[j] - annual_mean$mean_pet[i]) /
               (annual_mean$year[j]   - annual_mean$year[i]), NA))
    sen_slope     <- median(slopes_mat, na.rm = TRUE)
    sen_intercept <- median(annual_mean$mean_pet - sen_slope * annual_mean$year, na.rm = TRUE)
    abline(a = sen_intercept, b = sen_slope, col = "blue", lwd = 2)
    S     <- sum(sign(outer(annual_mean$mean_pet, annual_mean$mean_pet, "-")[lower.tri(matrix(0, n, n))]))
    var_S <- n * (n - 1) * (2 * n + 5) / 18
    Z     <- ifelse(S > 0, (S - 1) / sqrt(var_S), ifelse(S < 0, (S + 1) / sqrt(var_S), 0))
    mk_p  <- 2 * pnorm(-abs(Z))
    log_event(sprintf("Sen's slope: %.6f mm/day/year (MK p=%.3f)", sen_slope, mk_p))
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

save_plot("F2e_annual_trends_PM_vs_ERA5Land", function() {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  y_r <- range(c(annual_pm_era5$mean_pet_pm, annual_pm_era5$mean_pet_era5land), na.rm = TRUE)
  plot(annual_pm_era5$year, annual_pm_era5$mean_pet_pm, type = "b", pch = 19, col = "blue", lwd = 1.5,
       ylim = y_r, xlab = "Year", ylab = "Annual Mean PET / PEV (mm/day)",
       main = "Annual PET/PEV Trends: PM vs. ERA5-Land")
  lines(annual_pm_era5$year, annual_pm_era5$mean_pet_era5land, type = "b", pch = 17, col = "darkorange", lwd = 1.5)
  lm_pm <- lm(mean_pet_pm ~ year, data = annual_pm_era5)
  lm_e5 <- lm(mean_pet_era5land ~ year, data = annual_pm_era5)
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
             "G coefficient (FAO-56 Eq.45)",
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