##############################################
# SCRIPT 2b — PET CALCULATION
# Following FAO Penman-Monteith monthly equation
#
# EXECUTION ORDER: Run AFTER 2a_detrend_temperature.R.
#   2a_detrend_temperature.R produces the detrended T2m and Tdew NetCDFs;
#   this script then reads both observed and detrended temperature fields and
#   computes observed PET (Sections 1–6) and counterfactual PET (Section 5c)
#   in a single run — no two-pass execution required.
#
# MODIFICATIONS vs. original:
#   1. Rs conversion fixed: ERA5-Land monthly-averaged SSRD is already
#      J/m2/day (mean daily) -- only divide by 1e6 (J->MJ), NOT days_in_month.
#   2. G coefficient: 0.14 retained (correct for FAO-56 Eq.45 T_i - T_{i-1}).
#      Cold-region constraint: |G| capped at 0.3 x |Rn| when T < 0 degC
#      (prevents over-estimated G under frozen/snow-covered soil from making
#       Rn-G strongly negative in winter months).
#   3. Thornthwaite PET (temperature-only) added as Section 5b.
#   4. safe_writeCDF() helper: retry/delete prevents Windows file-lock errors.
#   5. Comparison diagnostic plots (PM vs. Thornthwaite) added.
#   6. [Counterfactual Branch] Section 5c: PET recomputed from month-specific
#      linearly detrended T2m and Tdew (outputs of 2a_detrend_temperature.R).
#      Produces: potential_evapotranspiration_PM_detrended.nc
#                potential_evapotranspiration_Thw_detrended.nc
#      Wind, radiation, and pressure remain at observed values so that only
#      the temperature-mediated warming effect on PET is isolated.
#      RUN_DETRENDED_BRANCH (Section 1) controls whether Section 5c executes.
##############################################
library(terra)
library(ncdf4)
library(lubridate)

# --- 1. CONFIGURATION ---
# Working directory remains: D:/Nechako_Drought/Nechako/
# (do NOT call setwd — caller controls the working directory)

MONTHLY_DIR   <- "D:/Nechako_Drought/Nechako/monthly_data_direct"  # NetCDF outputs + summary CSVs
PET_CALCS_DIR <- "D:/Nechako_Drought/Nechako/pet_calcs"            # diagnostics, aux rasters, log

# OUTPUT LOCATION POLICY
# ─────────────────────────────────────────────────────────────────────────────
# NetCDF rasters  → MONTHLY_DIR  (consumed by 3SPEI_ERALand.R, 10a, 10b)
# Summary CSVs    → MONTHLY_DIR  (consumed by 7basin_timeseries.R Part 4d
#                                 and 10b_PET_bias_nonstationarity.R)
# Auxiliary tifs  → PET_CALCS_DIR  (Thornthwaite I, a, N rasters)
# Diagnostic PDFs → PET_CALCS_DIR/diagnostics/
# Processing log  → PET_CALCS_DIR/
# ─────────────────────────────────────────────────────────────────────────────

# Counterfactual branch — requires 2a_detrend_temperature.R to have been run first.
# Set TRUE  : Section 5c computes counterfactual PET from detrended T2m/Tdew.
# Set FALSE : Section 5c is skipped; only observed PET (Sections 1–6) is produced.
# Default is TRUE because 2a_detrend_temperature.R is step 1 in this pipeline.
RUN_DETRENDED_BRANCH <- TRUE

# Create output directories if they do not yet exist
if (!dir.exists(MONTHLY_DIR))   dir.create(MONTHLY_DIR,   recursive = TRUE, showWarnings = FALSE)
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

# --- 2. LOAD DATA ---
log_event("Loading ERA5-Land data files...")

tryCatch({
  t2m_b  <- rast(file.path(MONTHLY_DIR, "2m_temperature_monthly.nc"))
  d2m_b  <- rast(file.path(MONTHLY_DIR, "2m_dewpoint_temperature_monthly.nc"))
  u10_b  <- rast(file.path(MONTHLY_DIR, "10m_u_component_of_wind_monthly.nc"))
  v10_b  <- rast(file.path(MONTHLY_DIR, "10m_v_component_of_wind_monthly.nc"))
  ssrd_b <- rast(file.path(MONTHLY_DIR, "surface_solar_radiation_downwards_monthly.nc"))
  sp_b   <- rast(file.path(MONTHLY_DIR, "surface_pressure_monthly.nc"))
  z_b    <- rast(file.path(MONTHLY_DIR, "geopotential_monthly.nc"))
  log_event("Successfully loaded all inputs.")
}, error = function(e) stop("Input file loading failed: ", e$message))

h_m           <- z_b[[1]] / G_GRAV
h_m           <- resample(h_m, t2m_b[[1]], method = "bilinear")  # <-- to avoid error "extents do not match"
avg_elevation <- global(h_m, "mean", na.rm = TRUE)[1, 1]
log_event(paste("Average basin elevation:", round(avg_elevation, 1), "m"))

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
log_event("... Wind Speed")
u2 <- ifel(is.na(sqrt(u10_b^2 + v10_b^2) * 4.87 / log(67.8 * 10 - 5.42)), 2.0,
           sqrt(u10_b^2 + v10_b^2) * 4.87 / log(67.8 * 10 - 5.42))
s  <- global(u2, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Wind speed (2m) range: %.2f to %.2f m/s (mean: %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))

# D. Solar Radiation
log_event("... Net Radiation")
# FIX: ERA5-Land monthly-averaged SSRD is already J/m2/day (mean daily).
# Divide by 1e6 only (J->MJ). Do NOT divide by days_in_month.
log_event("  Rs = ssrd_b / 1e6  (monthly-averaged product already per day)")
Rs <- ssrd_b / 1e6

s <- global(Rs, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  Rs range: %.2f to %.2f MJ/m2/day (annual mean: %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))
if (mean(s$mean) < 8 || mean(s$mean) > 20) {
  log_event("  WARNING: Rs annual mean outside expected 8-20 MJ/m2/day!")
} else {
  log_event("  Rs in expected range for ~54 degN (8-20 MJ/m2/day)")
}

log_event(sprintf("  January Rs: %.2f, July Rs: %.2f MJ/m2/day",
                  global(Rs[[1]], "mean", na.rm = TRUE)[1, 1],
                  global(Rs[[7]], "mean", na.rm = TRUE)[1, 1]))

Rs  <- ifel(is.na(Rs), K_RS * sqrt(10) * Ra_stack, Rs)
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
# --- 5b. THORNTHWAITE PET (Temperature-Only) ---
# Pure thermodynamic signal for Dynamic vs. Thermodynamic decomposition.
# Formula:  PET_thw = 16 x (10T/I)^a x (N/12) x (d/30)   [mm/day]
#           PET_thw = 0  when T <= 0 degC
# I, a computed per pixel from 1991-2020 climatological T.
# N reuses Ra-loop sunset hour angle geometry (no duplication).
# ==============================================================================
log_event("==========================================")
log_event("Starting THORNTHWAITE PET (temperature-only)...")
log_event("==========================================")

# 5b-1. Climatological monthly mean T (1991-2020) per pixel
log_event("... 1991-2020 climatological T per pixel")
clim_idx_thw <- which(as.numeric(format(dates, "%Y")) >= 1991 &
                        as.numeric(format(dates, "%Y")) <= 2020)

T_clim_12 <- rast(lapply(1:12, function(m) {
  idx_m <- clim_idx_thw[month_indices[clim_idx_thw] == m]
  if (length(idx_m) == 0) { warning(sprintf("No clim months for %d", m)); return(T_mean[[1]] * NA) }
  mean(T_mean[[idx_m]], na.rm = TRUE)
}))
names(T_clim_12) <- month.abb

s <- global(T_clim_12, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  Climatological T: %.2f to %.2f degC (mean: %.2f)",
                  min(s$min), max(s$max), mean(s$mean)))

# 5b-2. Annual heat index I (per pixel)
log_event("... Heat index I per pixel")
# NOTE: use T_clim_12[[1]] * 0 (NOT rast(...)) to keep spatial values
I_rast <- T_clim_12[[1]] * 0
for (m in 1:12)
  I_rast <- I_rast + (ifel(T_clim_12[[m]] <= 0, 0, T_clim_12[[m]]) / 5)^1.514
names(I_rast) <- "heat_index_I"
I_rast        <- ifel(I_rast <= 0, NA_real_, I_rast)  # guard Inf for all-cold pixels

s <- global(I_rast, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  I: min=%.2f, max=%.2f, mean=%.2f", s$min, s$max, s$mean))

# 5b-3. Exponent a (per pixel)
log_event("... Exponent a per pixel")
a_rast <- 6.75e-7 * I_rast^3 - 7.71e-5 * I_rast^2 + 1.792e-2 * I_rast + 0.49239
names(a_rast) <- "exponent_a"
s <- global(a_rast, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  a: min=%.4f, max=%.4f, mean=%.4f", s$min, s$max, s$mean))

# 5b-4. Monthly daylight hours N_m (reuse Ra-loop geometry)
log_event("... Daylight hours N_m from Ra geometry")
N_list <- list()
for (m in 1:12) {
  J      <- mid_month_days[m]
  delta  <- 0.409 * sin(2 * pi * J / 365 - 1.39)
  tp     <- ifel(tan(lat_rad) * tan(delta) >  1,  1,
                 ifel(tan(lat_rad) * tan(delta) < -1, -1,
                      tan(lat_rad) * tan(delta)))
  N_list[[m]] <- (24 / pi) * acos(-tp)
}
N_clim_12 <- rast(N_list)
names(N_clim_12) <- month.abb
s <- global(N_clim_12, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("  N: DJF mean=%.1f hr, JJA mean=%.1f hr",
                  mean(s$mean[c(1, 2, 12)]), mean(s$mean[6:8])))

# 5b-5. Apply formula to all 912 time steps
log_event(sprintf("... Applying Thornthwaite to all %d time steps", n_months))
et0_thw <- rast(lapply(seq_len(n_months), function(i) {
  m   <- month_indices[i]
  d_m <- days_in_month[i]
  T_i <- T_mean[[i]]
  ratio     <- ifel(T_i <= 0, 0, ifel((10 * T_i) / I_rast < 0, 0, (10 * T_i) / I_rast))
  pet_unadj <- ifel(T_i <= 0, 0, ifel(is.na(I_rast), 0, 16 * ratio^a_rast))
  pet_month <- ifel(pet_unadj * (N_clim_12[[m]] / 12) * (d_m / 30) < 0, 0,
                    pet_unadj * (N_clim_12[[m]] / 12) * (d_m / 30))
  pet_month / d_m   # mm/month -> mm/day
}))
terra::time(et0_thw) <- dates
names(et0_thw) <- sprintf("PET_THW_%04d_%02d",
                          as.numeric(format(dates, "%Y")),
                          as.numeric(format(dates, "%m")))

s <- global(et0_thw, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Thornthwaite PET: %.3f to %.3f mm/day (mean: %.3f)",
                  min(s$min), max(s$max), mean(s$mean)))
log_event(sprintf("  Jan=%.3f mm/day, Jul=%.3f mm/day",
                  global(et0_thw[[1]], "mean", na.rm = TRUE)[1, 1],
                  global(et0_thw[[7]], "mean", na.rm = TRUE)[1, 1]))
if (global(et0_thw[[7]], "mean", na.rm = TRUE)[1, 1] >
    global(et0_thw[[1]], "mean", na.rm = TRUE)[1, 1]) {
  log_event("  Seasonal cycle correct (July > January)")
} else {
  log_event("  WARNING: July <= January -- check latitude or T data!")
}

# 5b-6. Save Thornthwaite NetCDF
log_event("Writing Thornthwaite PET NetCDF...")
nc_file_thw <- file.path(MONTHLY_DIR, "potential_evapotranspiration_thornthwaite_monthly.nc")
safe_writeCDF(et0_thw, nc_file_thw,
              varname  = "pet_thw",
              longname = "Reference Evapotranspiration (Thornthwaite 1948, temperature-only)",
              unit     = "mm/day")
log_event(sprintf("Thornthwaite PET saved: %s", nc_file_thw))

# 5b-7. Save auxiliary rasters
writeRaster(I_rast,    file.path(PET_CALCS_DIR, "thornthwaite_heat_index_I.tif"),  overwrite = TRUE)
writeRaster(a_rast,    file.path(PET_CALCS_DIR, "thornthwaite_exponent_a.tif"),    overwrite = TRUE)
writeRaster(N_clim_12, file.path(PET_CALCS_DIR, "thornthwaite_daylight_N.tif"),    overwrite = TRUE)
log_event("Thornthwaite auxiliary rasters saved (I, a, N)")

# 5b-8. Basin-summary statistics for Thornthwaite
global_stats_thw <- global(et0_thw, c("mean", "min", "max", "sd"), na.rm = TRUE)
pet_thw_summary  <- data.frame(
  date     = format(dates, "%Y-%m"),
  year     = as.numeric(format(dates, "%Y")),
  month    = as.numeric(format(dates, "%m")),
  mean_pet = global_stats_thw$mean,
  min_pet  = global_stats_thw$min,
  max_pet  = global_stats_thw$max,
  std_dev  = global_stats_thw$sd
)
# Written to MONTHLY_DIR so 7basin_timeseries.R (Part 4d) and
# 10b_PET_bias_nonstationarity.R can locate it without path changes.
write.csv(pet_thw_summary,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_Thornthwaite_monthly_summary.csv"),
          row.names = FALSE)
log_event(sprintf("  Thornthwaite summary CSV → %s/ERA5Land_Nechako_PET_Thornthwaite_monthly_summary.csv",
                  basename(MONTHLY_DIR)))

monthly_stats_thw <- aggregate(
  pet_thw_summary[, c("mean_pet", "min_pet", "max_pet", "std_dev")],
  by = list(month = pet_thw_summary$month), FUN = mean, na.rm = TRUE)
monthly_stats_thw$month_name <- month.name[monthly_stats_thw$month]
write.csv(monthly_stats_thw,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_Thornthwaite_calendar_month_stats.csv"),
          row.names = FALSE)
log_event(sprintf("  Thornthwaite calendar-month stats CSV → %s/ERA5Land_Nechako_PET_Thornthwaite_calendar_month_stats.csv",
                  basename(MONTHLY_DIR)))

log_event("==========================================")
log_event("THORNTHWAITE PET COMPLETE")
log_event(sprintf("  Mean: %.3f mm/day", mean(pet_thw_summary$mean_pet, na.rm = TRUE)))
log_event("==========================================")

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
log_event(sprintf("  PM summary CSV → %s/ERA5Land_Nechako_PET_monthly_summary.csv",
                  basename(MONTHLY_DIR)))

log_event("Calculating calendar month climatology...")
monthly_stats <- aggregate(
  pet_summary[, c("mean_pet", "min_pet", "max_pet", "std_dev")],
  by = list(month = pet_summary$month), FUN = mean, na.rm = TRUE)
monthly_stats$month_name <- month.name[monthly_stats$month]
write.csv(monthly_stats,
          file.path(MONTHLY_DIR, "ERA5Land_Nechako_PET_calendar_month_stats.csv"),
          row.names = FALSE)
log_event(sprintf("  PM calendar-month stats CSV → %s/ERA5Land_Nechako_PET_calendar_month_stats.csv",
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
    log_event(sprintf("    WARNING: PNG failed for %s — %s", name, conditionMessage(e)))
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

# F2. Thornthwaite vs. Penman-Monteith comparison
# Scientific interpretation:
#   Large PM-Thw summer bias  => radiation + wind also drive PET (not purely T)
#   Small PM-Thw summer bias  => T dominates PET => Thornthwaite is good proxy
log_event("... F2: PM vs. Thornthwaite comparison")

# F2a. Time series overlay
save_plot("F2a_timeseries_PM_vs_Thw", function() {
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
  y_lim <- range(c(pet_summary$mean_pet, pet_thw_summary$mean_pet), na.rm = TRUE)
  plot(dates, pet_summary$mean_pet, type = "l", col = "blue", lwd = 2,
       ylim = y_lim, xlab = "Date", ylab = "PET (mm/day)",
       main = "Basin-Average PET: Penman-Monteith vs. Thornthwaite",
       cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, font.main = 2)
  lines(dates, pet_thw_summary$mean_pet, col = "darkorange", lwd = 2, lty = 2)
  legend("topright",
         legend = c("Penman-Monteith", "Thornthwaite (T-only)"),
         col = c("blue", "darkorange"), lwd = 2, lty = c(1, 2), bty = "n", cex = 1.2)
  grid(lty = "dotted", col = "grey88")
})

# F2b. Seasonal climatology with error bars — UNIT FIX: convert mm/day -> mm/month
log_event("... F2b: Seasonal PET climatology PM vs Thw (mm/month)")
pet_clim_pm  <- monthly_stats$mean_pet     * days_per_month
pet_clim_thw <- monthly_stats_thw$mean_pet * days_per_month
sd_pm        <- monthly_stats$std_dev      * days_per_month
sd_thw       <- monthly_stats_thw$std_dev  * days_per_month

save_plot("F2b_seasonal_climatology_PM_vs_Thw", function() {
  # ── Increased font sizes for readability ──────────────────────────────────
  CEX_MAIN <- 1.65   # plot title
  CEX_LAB  <- 1.60   # axis labels (x / y)
  CEX_AXIS <- 1.50   # tick labels (both x month names and y numeric values)
  CEX_LEG  <- 1.25   # legend text
  LWD_LINE <- 2.5    # line width
  
  par(mfrow = c(1, 1), mar = c(5, 6, 4, 2))
  y_max <- max(c(pet_clim_pm, pet_clim_thw), na.rm = TRUE) * 1.15
  
  # ── Base plot: Penman-Monteith (blue solid, filled circles) ───────────────
  # cex.axis here controls the y-axis tick numbers; x-axis is redrawn below.
  plot(1:12, pet_clim_pm,
       type = "b", pch = 19, col = "blue", lwd = LWD_LINE,
       xlab = "Month", ylab = "Mean PET (mm/month)",
       main = "Seasonal PET Climatology: PM vs. Thornthwaite",
       xaxt = "n", ylim = c(0, y_max),
       cex.main = CEX_MAIN, cex.lab = CEX_LAB, cex.axis = CEX_AXIS, font.main = 2)
  axis(1, at = 1:12, labels = month.abb, cex.axis = CEX_AXIS)
  
  # ── Overlay: Thornthwaite (orange dashed, filled triangles) ──────────────
  lines(1:12, pet_clim_thw,
        type = "b", pch = 17, col = "darkorange", lwd = LWD_LINE, lty = 2)
  
  # ── Error bars (±1 SD) ────────────────────────────────────────────────────
  suppressWarnings({
    arrows(1:12, pet_clim_pm  - sd_pm,  1:12, pet_clim_pm  + sd_pm,
           angle = 90, code = 3, length = 0.05, col = "blue",       lwd = 1.5)
    arrows(1:12, pet_clim_thw - sd_thw, 1:12, pet_clim_thw + sd_thw,
           angle = 90, code = 3, length = 0.05, col = "darkorange", lwd = 1.5)
  })
  
  # ── Legend (label fixed: "Penman-Monteith" not "Penman-Monteith (ERA5)") ─
  legend("topleft",
         legend = c("Penman-Monteith", "Thornthwaite (T-only)"),
         col    = c("blue", "darkorange"),
         lwd    = LWD_LINE,
         pch    = c(19L, 17L),
         lty    = c(1L, 2L),
         bty    = "n",
         cex    = CEX_LEG)
  
  # ── Light dotted reference grid ───────────────────────────────────────────
  grid(lty = "dotted", col = "grey88", lwd = 0.8)
})

# F2c. Bias bar chart — UNIT FIX: convert mm/day -> mm/month
log_event("... F2c: Monthly PET bias PM minus Thw (mm/month)")
bias_monthly <- (monthly_stats$mean_pet - monthly_stats_thw$mean_pet) * days_per_month

save_plot("F2c_monthly_bias_PM_minus_Thw", function() {
  par(mfrow = c(1, 1), mar = c(4, 5, 3, 1))
  barplot(bias_monthly, names.arg = month.abb,
          col    = ifelse(bias_monthly > 0, "#e31a1c", "#1f78b4"),
          border = "white",
          xlab = "Month", ylab = "PET bias (mm/month)",
          main = "PM minus Thornthwaite PET\n(red = PM higher [radiation/wind contribution])")
  abline(h = 0, lwd = 1.5)
  grid(nx = NA, ny = NULL)
  legend("topleft",
         legend = c("PM > Thornthwaite", "Thornthwaite > PM"),
         fill = c("#e31a1c", "#1f78b4"), bty = "n")
})

# F2d. Annual scatter
annual_pm_thw <- merge(
  aggregate(mean_pet ~ year, data = pet_summary,     FUN = mean, na.rm = TRUE),
  aggregate(mean_pet ~ year, data = pet_thw_summary, FUN = mean, na.rm = TRUE),
  by = "year", suffixes = c("_pm", "_thw"))
r_pm_thw <- cor(annual_pm_thw$mean_pet_pm, annual_pm_thw$mean_pet_thw,
                use = "complete.obs")
save_plot("F2d_annual_scatter_PM_vs_Thw", function() {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  plot(annual_pm_thw$mean_pet_thw, annual_pm_thw$mean_pet_pm,
       pch = 19, col = "#4393c3", cex = 0.9,
       xlab = "Thornthwaite PET (mm/day)", ylab = "Penman-Monteith PET (mm/day)",
       main = sprintf("Annual Mean PET: PM vs. Thornthwaite  (r = %.3f)", r_pm_thw))
  abline(0, 1, col = "black", lty = 2, lwd = 1.5)
  abline(lm(mean_pet_pm ~ mean_pet_thw, data = annual_pm_thw), col = "red", lwd = 2)
  legend("topleft", legend = c("1:1 line", "Linear fit"),
         col = c("black", "red"), lty = c(2, 1), lwd = 2, bty = "n")
  grid()
})

# F2e. Annual trend comparison
save_plot("F2e_annual_trends_PM_vs_Thw", function() {
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
  y_r <- range(c(annual_pm_thw$mean_pet_pm, annual_pm_thw$mean_pet_thw), na.rm = TRUE)
  plot(annual_pm_thw$year, annual_pm_thw$mean_pet_pm, type = "b",
       pch = 19, col = "blue", lwd = 1.5, ylim = y_r,
       xlab = "Year", ylab = "Annual Mean PET (mm/day)",
       main = "Annual PET Trends: PM vs. Thornthwaite")
  lines(annual_pm_thw$year, annual_pm_thw$mean_pet_thw,
        type = "b", pch = 17, col = "darkorange", lwd = 1.5)
  lm_pm  <- lm(mean_pet_pm  ~ year, data = annual_pm_thw)
  lm_thw <- lm(mean_pet_thw ~ year, data = annual_pm_thw)
  abline(lm_pm,  col = "blue",       lty = 2, lwd = 2)
  abline(lm_thw, col = "darkorange", lty = 2, lwd = 2)
  legend("topleft",
         legend = c(sprintf("PM: %+.4f mm/day/yr (p=%.3f)",
                            coef(lm_pm)[2],  summary(lm_pm)$coefficients[2, 4]),
                    sprintf("Thornthwaite: %+.4f mm/day/yr (p=%.3f)",
                            coef(lm_thw)[2], summary(lm_thw)$coefficients[2, 4])),
         col = c("blue", "darkorange"), pch = c(19, 17), lty = 1, lwd = 2, bty = "n")
  grid()
})

# F2f. Thornthwaite sample spatial maps
save_plot("F2f_Thw_spatial_maps", function() {
  if (nlyr(et0_thw) >= 12) {
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))
    for (idx in sample_months)
      if (idx <= nlyr(et0_thw))
        plot(et0_thw[[idx]], col = rev(hcl.colors(100, "YlOrRd")),
             main = sprintf("PET_Thw %s", names(et0_thw)[idx]))
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
             "Thw Mean PET (mm/day)", "Thw Max",
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
    round(mean(pet_thw_summary$mean_pet, na.rm = TRUE), 3),
    round(max(pet_thw_summary$max_pet,   na.rm = TRUE), 3),
    round(sum(pet_summary$mean_pet > 8,  na.rm = TRUE) / nrow(pet_summary) * 100, 2),
    "0.14 x (T_i - T_{i-1})",
    "|G| <= 0.3 x |Rn| when T < 0 C",
    "ssrd/1e6  [J/m2/day -> MJ/m2/day]"
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
# --- 5c. COUNTERFACTUAL PET FROM DETRENDED T2m / Tdew  [NEW] ---
# Methodology: identical PM + Thornthwaite calculation as Sections 5 and 5b,
# but T2m and Tdew are replaced by their month-specific linearly detrended
# counterparts produced by 2a_detrend_temperature.R.  All other forcing
# (wind, radiation, pressure) remains at observed values so that only the
# temperature-mediated warming signal is removed from PET.
#
# Outputs:
#   monthly_data_direct/potential_evapotranspiration_PM_detrended.nc
#   monthly_data_direct/potential_evapotranspiration_Thw_detrended.nc
# ==============================================================================
if (RUN_DETRENDED_BRANCH) {
  
  log_event("==========================================")
  log_event("Section 5c: COUNTERFACTUAL PET (detrended T2m + Tdew)")
  log_event("==========================================")
  
  # ── 5c-0. Check prerequisites ───────────────────────────────────────────────
  f_t2m_det  <- file.path(MONTHLY_DIR, "2m_temperature_detrended_monthly.nc")
  f_tdew_det <- file.path(MONTHLY_DIR, "2m_dewpoint_detrended_monthly.nc")
  
  if (!file.exists(f_t2m_det) || !file.exists(f_tdew_det)) {
    log_event("  WARNING: Detrended NC files not found — skipping Section 5c.")
    log_event(sprintf("  Expected: %s", f_t2m_det))
    log_event(sprintf("  Expected: %s", f_tdew_det))
    log_event("  Run 2a_detrend_temperature.R first, then re-run this script.")
  } else {
    
    # ── 5c-1. Load detrended temperature fields ────────────────────────────────
    log_event("  Loading detrended T2m and Tdew...")
    t2m_det_raw  <- rast(f_t2m_det)
    tdew_det_raw <- rast(f_tdew_det)
    
    # Align spatial grid to the already-loaded observed T2m (t2m_b)
    if (!same.crs(t2m_det_raw, t2m_b) || !compareGeom(t2m_det_raw, t2m_b, stopOnError = FALSE)) {
      t2m_det_raw  <- resample(t2m_det_raw,  t2m_b, method = "bilinear")
      tdew_det_raw <- resample(tdew_det_raw, t2m_b, method = "bilinear")
    }
    
    # 2b outputs are in Kelvin (same units as ERA5-Land originals) — convert to °C
    T_det  <- t2m_det_raw  - 273.15
    Td_det <- tdew_det_raw - 273.15
    
    # Physics guard: Tdew must not exceed T (RH <= 100%)
    n_viol <- global(Td_det > T_det, "sum", na.rm = TRUE)
    n_viol_total <- sum(n_viol[, 1], na.rm = TRUE)
    if (n_viol_total > 0L) {
      log_event(sprintf("  Clamping %d Tdew_det > T_det instances to preserve VPD > 0.",
                        n_viol_total))
      Td_det <- ifel(Td_det > T_det, T_det, Td_det)
    } else {
      log_event("  Physics check passed: Tdew_det <= T_det everywhere.")
    }
    
    # Verify time dimension matches observed data
    n_det <- nlyr(T_det)
    if (n_det != n_months) {
      log_event(sprintf("  WARNING: Detrended stack has %d layers; observed has %d.",
                        n_det, n_months))
      log_event("  Truncating/padding to match observed time axis.")
      if (n_det > n_months) {
        T_det  <- T_det[[1:n_months]]
        Td_det <- Td_det[[1:n_months]]
      }
    }
    
    s <- global(T_det,  c("min", "max", "mean"), na.rm = TRUE)
    log_event(sprintf("  Detrended T:    %.2f to %.2f degC (mean %.2f)",
                      min(s$min), max(s$max), mean(s$mean)))
    s <- global(Td_det, c("min", "max", "mean"), na.rm = TRUE)
    log_event(sprintf("  Detrended Tdew: %.2f to %.2f degC (mean %.2f)",
                      min(s$min), max(s$max), mean(s$mean)))
    
    # ── 5c-2. Penman-Monteith with detrended T (observed Rn, G uses det T) ────
    log_event("  Computing PM PET (detrended T)...")
    
    # Saturation and actual vapour pressures from detrended temperatures
    es_det    <- 0.6108 * exp((17.27 * T_det)  / (T_det  + 237.3))
    ea_det    <- 0.6108 * exp((17.27 * Td_det) / (Td_det + 237.3))
    # Slope of saturation VP curve at detrended T
    delta_det <- (4098 * es_det) / (T_det + 237.3)^2
    
    # Ground heat flux using detrended T sequence (FAO-56 Eq.45)
    if (n_months > 1L) {
      t_prev_det  <- c(T_det[[1]] * 0, T_det[[1:(n_months - 1)]])
      G_det       <- 0.14 * (T_det - t_prev_det)
      G_det[[1]]  <- G_det[[1]] * 0   # first layer: no prior month
      # Cold-region constraint with detrended T
      G_cap_det   <- 0.3 * abs(Rn)    # Rn is unchanged (observed)
      G_det       <- ifel(T_det < 0,
                          ifel(G_det >  G_cap_det,  G_cap_det,
                               ifel(G_det < -G_cap_det, -G_cap_det, G_det)),
                          G_det)
    } else {
      G_det <- T_det * 0
    }
    
    # PM numerator & denominator  (u2, gamma, Rn are all unchanged/observed)
    num_det <- 0.408 * delta_det * (Rn - G_det) +
      gamma * (900 / (T_det + 273)) * u2 * (es_det - ea_det)
    den_det <- delta_det + gamma * (1 + 0.34 * u2)
    et0_det <- ifel(num_det / den_det < 0, 0, num_det / den_det)
    
    terra::time(et0_det) <- dates
    names(et0_det) <- sprintf("PET_DET_%04d_%02d",
                              as.numeric(format(dates, "%Y")),
                              as.numeric(format(dates, "%m")))
    
    s <- global(et0_det, c("min", "max", "mean"), na.rm = TRUE)
    log_event(sprintf("  PM PET (det): %.3f to %.3f mm/day (mean %.3f)",
                      min(s$min), max(s$max), mean(s$mean)))
    
    # ── 5c-3. Thornthwaite with detrended T ──────────────────────────────────
    log_event("  Computing Thornthwaite PET (detrended T)...")
    
    # 1991-2020 climatological mean T from the DETRENDED series (per pixel)
    clim_idx_det <- which(as.numeric(format(dates, "%Y")) >= 1991 &
                            as.numeric(format(dates, "%Y")) <= 2020)
    T_clim_det_12 <- rast(lapply(1:12, function(m) {
      idx_m <- clim_idx_det[month_indices[clim_idx_det] == m]
      if (length(idx_m) == 0) return(T_det[[1]] * NA)
      mean(T_det[[idx_m]], na.rm = TRUE)
    }))
    names(T_clim_det_12) <- month.abb
    
    # Heat index I from detrended climatology
    I_det <- T_clim_det_12[[1]] * 0
    for (m in 1:12)
      I_det <- I_det + (ifel(T_clim_det_12[[m]] <= 0, 0, T_clim_det_12[[m]]) / 5)^1.514
    I_det <- ifel(I_det <= 0, NA_real_, I_det)
    
    # Exponent a from detrended I
    a_det <- 6.75e-7 * I_det^3 - 7.71e-5 * I_det^2 + 1.792e-2 * I_det + 0.49239
    
    # Apply Thornthwaite to all time steps (reuse N_clim_12 from Section 5b-4)
    et0_thw_det <- rast(lapply(seq_len(n_months), function(i) {
      m_i   <- month_indices[i]
      d_m   <- days_in_month[i]
      T_i   <- T_det[[i]]
      ratio     <- ifel(T_i <= 0, 0,
                        ifel((10 * T_i) / I_det < 0, 0, (10 * T_i) / I_det))
      pet_unadj <- ifel(T_i <= 0, 0,
                        ifel(is.na(I_det), 0, 16 * ratio^a_det))
      pet_month <- ifel(pet_unadj * (N_clim_12[[m_i]] / 12) * (d_m / 30) < 0, 0,
                        pet_unadj * (N_clim_12[[m_i]] / 12) * (d_m / 30))
      pet_month / d_m   # mm/month -> mm/day
    }))
    terra::time(et0_thw_det) <- dates
    names(et0_thw_det) <- sprintf("PET_THW_DET_%04d_%02d",
                                  as.numeric(format(dates, "%Y")),
                                  as.numeric(format(dates, "%m")))
    
    s <- global(et0_thw_det, c("min", "max", "mean"), na.rm = TRUE)
    log_event(sprintf("  Thw PET (det): %.3f to %.3f mm/day (mean %.3f)",
                      min(s$min), max(s$max), mean(s$mean)))
    
    # ── 5c-4. Save NetCDFs ────────────────────────────────────────────────────
    log_event("  Writing detrended PET NetCDF files...")
    
    safe_writeCDF(
      et0_det,
      file.path(MONTHLY_DIR, "potential_evapotranspiration_PM_detrended.nc"),
      varname  = "pet_det",
      longname = "Reference Evapotranspiration (FAO-56 Penman-Monteith, detrended T2m+Tdew)",
      unit     = "mm/day")
    log_event("  Saved: potential_evapotranspiration_PM_detrended.nc")
    
    safe_writeCDF(
      et0_thw_det,
      file.path(MONTHLY_DIR, "potential_evapotranspiration_Thw_detrended.nc"),
      varname  = "pet_thw_det",
      longname = "Reference Evapotranspiration (Thornthwaite 1948, detrended T2m)",
      unit     = "mm/day")
    log_event("  Saved: potential_evapotranspiration_Thw_detrended.nc")
    
    # ── 5c-5. Quick comparison diagnostics (obs vs. detrended) ───────────────
    log_event("  Generating obs-vs-detrended comparison diagnostics...")
    
    # Basin-average annual mean PET for each branch
    ann_pm_obs  <- aggregate(global(et0,         "mean", na.rm = TRUE)[, 1],
                             by = list(year = as.numeric(format(dates, "%Y"))), FUN = mean)
    ann_pm_det  <- aggregate(global(et0_det,     "mean", na.rm = TRUE)[, 1],
                             by = list(year = as.numeric(format(dates, "%Y"))), FUN = mean)
    ann_thw_obs <- aggregate(global(et0_thw,     "mean", na.rm = TRUE)[, 1],
                             by = list(year = as.numeric(format(dates, "%Y"))), FUN = mean)
    ann_thw_det <- aggregate(global(et0_thw_det, "mean", na.rm = TRUE)[, 1],
                             by = list(year = as.numeric(format(dates, "%Y"))), FUN = mean)
    
    # Monthly climatology (12-value vectors)
    clim_pm_obs  <- tapply(global(et0,         "mean", na.rm = TRUE)[, 1], month_indices, mean)
    clim_pm_det  <- tapply(global(et0_det,     "mean", na.rm = TRUE)[, 1], month_indices, mean)
    clim_thw_obs <- tapply(global(et0_thw,     "mean", na.rm = TRUE)[, 1], month_indices, mean)
    clim_thw_det <- tapply(global(et0_thw_det, "mean", na.rm = TRUE)[, 1], month_indices, mean)
    
    det_diag_dir <- file.path(PET_CALCS_DIR, "detrended_diagnostics")
    if (!dir.exists(det_diag_dir)) dir.create(det_diag_dir, showWarnings = FALSE)
    
    tryCatch({
      pdf(file.path(det_diag_dir, "PET_obs_vs_detrended.pdf"), width = 14, height = 10)
      par(mfrow = c(2, 2), mar = c(4, 4.5, 3.5, 1.5), oma = c(0, 0, 2, 0))
      
      # Panel 1: PM seasonal climatology (mm/day)
      y_lim <- range(c(clim_pm_obs, clim_pm_det), na.rm = TRUE) * c(0.9, 1.15)
      plot(1:12, clim_pm_obs, type = "b", pch = 19, col = "blue", lwd = 2,
           xaxt = "n", ylim = y_lim,
           xlab = "Month", ylab = "PET (mm/day)",
           main = "PM: Observed vs. Detrended", font.main = 2)
      lines(1:12, clim_pm_det, type = "b", pch = 17, col = "firebrick2", lwd = 2, lty = 2)
      axis(1, at = 1:12, labels = month.abb)
      legend("topleft",
             legend = c("Observed T", "Detrended T"),
             col = c("blue", "firebrick2"), pch = c(19, 17), lty = c(1, 2), lwd = 2, bty = "n")
      grid(lty = "dotted", col = "grey88")
      
      # Panel 2: Thornthwaite seasonal climatology (mm/day)
      y_lim2 <- range(c(clim_thw_obs, clim_thw_det), na.rm = TRUE) * c(0.9, 1.15)
      plot(1:12, clim_thw_obs, type = "b", pch = 19, col = "darkorange", lwd = 2,
           xaxt = "n", ylim = y_lim2,
           xlab = "Month", ylab = "PET (mm/day)",
           main = "Thornthwaite: Observed vs. Detrended", font.main = 2)
      lines(1:12, clim_thw_det, type = "b", pch = 17, col = "purple3", lwd = 2, lty = 2)
      axis(1, at = 1:12, labels = month.abb)
      legend("topleft",
             legend = c("Observed T", "Detrended T"),
             col = c("darkorange", "purple3"), pch = c(19, 17), lty = c(1, 2), lwd = 2, bty = "n")
      grid(lty = "dotted", col = "grey88")
      
      # Panel 3: Monthly PM warming signal removed (mm/day)
      delta_pm <- clim_pm_obs - clim_pm_det
      barplot(delta_pm, names.arg = month.abb,
              col    = ifelse(delta_pm >= 0, "#d73027", "#4575b4"),
              border = "white",
              xlab = "Month", ylab = "ΔPET (mm/day)",
              main = "PM Warming Signal Removed\n(obs − detrended; red = warming raises PET)",
              font.main = 2)
      abline(h = 0, lwd = 1.5)
      
      # Panel 4: Annual PM trends (obs vs. detrended)
      yrs <- ann_pm_obs$year
      y_r <- range(c(ann_pm_obs$x, ann_pm_det$x), na.rm = TRUE)
      plot(yrs, ann_pm_obs$x, type = "l", col = "blue", lwd = 1.5,
           ylim = y_r, xlab = "Year", ylab = "Annual Mean PET (mm/day)",
           main = "PM Annual PET: Observed vs. Detrended", font.main = 2)
      lines(yrs, ann_pm_det$x, col = "firebrick2", lwd = 1.5, lty = 2)
      lm_obs <- lm(x ~ year, data = ann_pm_obs)
      lm_det <- lm(x ~ year, data = ann_pm_det)
      abline(lm_obs, col = "blue",       lty = 3, lwd = 2)
      abline(lm_det, col = "firebrick2", lty = 3, lwd = 2)
      legend("topleft",
             legend = c(
               sprintf("Obs: %+.4f mm/day/yr", coef(lm_obs)[2]),
               sprintf("Det: %+.4f mm/day/yr", coef(lm_det)[2])),
             col = c("blue", "firebrick2"), lwd = 2, lty = c(1, 2), bty = "n")
      grid(lty = "dotted", col = "grey88")
      
      mtext("PET: Observed vs. Detrended (counterfactual branch, Script 2b — Section 5c)",
            outer = TRUE, cex = 1.1, font = 2)
      dev.off()
      log_event(sprintf("  Diagnostic PDF saved: %s/PET_obs_vs_detrended.pdf", det_diag_dir))
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      log_event(sprintf("  WARNING: Diagnostic PDF failed — %s", conditionMessage(e)))
    })
    
    log_event("==========================================")
    log_event("Section 5c COMPLETE")
    log_event(sprintf("  PM  PET (det) mean: %.3f mm/day",
                      mean(global(et0_det,     "mean", na.rm = TRUE)[, 1], na.rm = TRUE)))
    log_event(sprintf("  Thw PET (det) mean: %.3f mm/day",
                      mean(global(et0_thw_det, "mean", na.rm = TRUE)[, 1], na.rm = TRUE)))
    log_event("==========================================")
    
  }  # end detrended-files-exist block
}  # end RUN_DETRENDED_BRANCH

# --- 8. FINISH ---
log_event("==========================================")
log_event(paste("Total months processed:", n_months))
log_event(paste("Penman-Monteith mean PET:",
                round(mean(pet_summary$mean_pet,     na.rm = TRUE), 2), "mm/day"))
log_event(paste("Thornthwaite mean PET:   ",
                round(mean(pet_thw_summary$mean_pet, na.rm = TRUE), 2), "mm/day"))
log_event(paste("PM layers with >5% negative pixels:", problematic_layers))
log_event("Outputs written to MONTHLY_DIR (monthly_data_direct/):")
log_event("  potential_evapotranspiration_monthly.nc")
log_event("  potential_evapotranspiration_thornthwaite_monthly.nc")
log_event("  ERA5Land_Nechako_PET_monthly_summary.csv")
log_event("  ERA5Land_Nechako_PET_calendar_month_stats.csv")
log_event("  ERA5Land_Nechako_PET_Thornthwaite_monthly_summary.csv")
log_event("  ERA5Land_Nechako_PET_Thornthwaite_calendar_month_stats.csv")
if (RUN_DETRENDED_BRANCH &&
    file.exists(file.path(MONTHLY_DIR, "potential_evapotranspiration_PM_detrended.nc"))) {
  log_event("  potential_evapotranspiration_PM_detrended.nc  [counterfactual]")
  log_event("  potential_evapotranspiration_Thw_detrended.nc [counterfactual]")
} else if (RUN_DETRENDED_BRANCH) {
  log_event("  [counterfactual outputs SKIPPED — detrended NC inputs were missing]")
} else {
  log_event("  [counterfactual branch not run — RUN_DETRENDED_BRANCH = FALSE]")
}
log_event("Diagnostics written to PET_CALCS_DIR (pet_calcs/):")
log_event("  diagnostics/PET_diagnostics.pdf  + individual PNGs")
log_event("  thornthwaite_heat_index_I.tif, thornthwaite_exponent_a.tif")
log_event("  thornthwaite_daylight_N.tif")
log_event("  monthly_pet_processing.log")
log_event("PROCESSING COMPLETE")
log_event(paste("Finished:", Sys.time()))