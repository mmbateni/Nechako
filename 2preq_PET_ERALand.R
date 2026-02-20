##############################################
# PET CALCULATION 
# Following FAO Penman-Monteith monthly equation
# MODIFIED: SSRD conversion + Soil heat flux fixes
##############################################
library(terra)
library(ncdf4)

# --- 1. CONFIGURATION ---
# Set working directory
setwd("D:/Nechako_Drought/Nechako/monthly_data_direct")

# Constants
G_GRAV <- 9.80665
K_RS <- 0.16  # Interior coefficient for Nechako Basin
LOG_FILE <- "monthly_pet_processing.log"

# Terra Options (Optimization for large files)
terraOptions(memfrac = 0.8, progress = 0) # Use 80% RAM, suppress progress bars for speed

# Initialize Log
cat("Monthly PET Processing Log (Vectorized)\n", file = LOG_FILE)
cat(paste("Processing started:", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)

log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, " | ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, " | ", msg))
}

# --- 2. LOAD DATA ---
log_event("Loading ERA5-Land data files...")

tryCatch({
  t2m_b  <- rast("2m_temperature_monthly.nc")
  d2m_b  <- rast("2m_dewpoint_temperature_monthly.nc")
  u10_b  <- rast("10m_u_component_of_wind_monthly.nc")
  v10_b  <- rast("10m_v_component_of_wind_monthly.nc")
  ssrd_b <- rast("surface_solar_radiation_downwards_monthly.nc")
  sp_b   <- rast("surface_pressure_monthly.nc")
  z_b    <- rast("geopotential_monthly.nc")
  
  log_event(paste("Successfully loaded all inputs."))
}, error = function(e) {
  stop("Input file loading failed: ", e$message)
})

# Static Elevation Calculation
h_m <- z_b[[1]] / G_GRAV 
avg_elevation <- global(h_m, "mean", na.rm = TRUE)[1,1]
log_event(paste("Average basin elevation:", round(avg_elevation, 1), "m"))


# --- 3. ROBUST DATE HANDLING WITH MULTIPLE FALLBACKS ---
log_event("Processing time information...")

extract_time_robust <- function(raster_obj, file_path, var_name) {
  # Method 1: Try terra's time()
  time_info <- time(raster_obj)
  
  if (inherits(time_info, "Date") || inherits(time_info, "POSIXt")) {
    if (!all(is.na(time_info))) {
      log_event(sprintf("  ✓ Time extracted via terra::time() for %s", var_name))
      return(as.Date(time_info))
    }
  }
  
  # Method 2: Try ncdf4
  log_event(sprintf("  Trying ncdf4 for %s...", var_name))
  tryCatch({
    nc <- nc_open(file_path)
    
    if ("time" %in% names(nc$var)) {
      time_data <- ncvar_get(nc, "time")
      time_units <- ncatt_get(nc, "time", "units")$value
      
      origin_match <- regexpr("since\\s+([0-9]{4}-[0-9]{2}-[0-9]{2})", time_units)
      if (origin_match > 0) {
        origin_str <- sub("since\\s+", "", regmatches(time_units, origin_match))
        origin_date <- as.Date(origin_str)
        
        if (grepl("hours", time_units, ignore.case=TRUE)) {
          time_vals <- origin_date + time_data / 24
        } else if (grepl("days", time_units, ignore.case=TRUE)) {
          time_vals <- origin_date + time_data
        } else {
          time_vals <- origin_date + time_data
        }
        
        nc_close(nc)
        log_event(sprintf("  ✓ Time extracted via ncdf4 for %s", var_name))
        return(as.Date(time_vals))
      }
    }
    nc_close(nc)
  }, error = function(e) {
    log_event(sprintf("  ncdf4 failed for %s: %s", var_name, conditionMessage(e)))
  })
  
  # Method 3: Reconstruct from layer count
  n_layers <- nlyr(raster_obj)
  log_event(sprintf("  Reconstructing time for %s (%d layers)...", var_name, n_layers))
  
  if (n_layers == 912) {
    start_date <- as.Date("1950-01-01")
    log_event("    Assuming 1950-01 to 2025-12")
  } else {
    start_date <- as.Date("1950-01-01")
    log_event(sprintf("    Assuming 1950-01 start, %d months", n_layers))
  }
  
  time_vals <- seq.Date(from=start_date, by="month", length.out=n_layers)
  log_event("  ✓ Time sequence reconstructed")
  return(time_vals)
}

# Extract time for temperature (use as reference)
t2m_file <- "2m_temperature_monthly.nc"
dates <- extract_time_robust(t2m_b, t2m_file, "Temperature")

n_months <- nlyr(t2m_b)
log_event(paste("Processing", n_months, "time steps."))
log_event(sprintf("Time range: %s to %s", format(min(dates)), format(max(dates))))


# Create Month Indices Vector
month_indices <- as.numeric(format(dates, "%m"))

# --- 3.5 CALCULATE DAYS IN MONTH (CRITICAL FOR SSRD CONVERSION) ---
log_event("Calculating days in each month for SSRD conversion...")

first_of_month <- as.Date(format(dates, "%Y-%m-01"))
first_next_month <- seq(first_of_month[1], by = "month", length.out = n_months + 1)[-1]
days_in_month <- as.integer(first_next_month - first_of_month)

log_event(sprintf("Days in month range: %d to %d", min(days_in_month), max(days_in_month)))
log_event(sprintf("January days: %d, July days: %d", days_in_month[1], days_in_month[7]))

# --- 4. VECTORIZED PRE-CALCULATIONS ---

# A. Extraterrestrial Radiation (Ra) - Vectorized
log_event("Calculating Ra (Extraterrestrial Radiation)...")

lat_r <- init(t2m_b, "y")
lat_rad <- lat_r * pi / 180

ra_list <- list()
mid_month_days <- c(16, 46, 75, 105, 135, 166, 196, 227, 258, 288, 319, 349)

for (m in 1:12) {
  J <- mid_month_days[m]
  delta <- 0.409 * sin(2 * pi * J / 365 - 1.39)
  dr <- 1 + 0.033 * cos(2 * pi * J / 365)
  
  tan_prod <- tan(lat_rad) * tan(delta)
  # Clamp values
  tan_prod <- ifel(tan_prod > 1, 1, ifel(tan_prod < -1, -1, tan_prod))
  ws <- acos(-tan_prod)
  
  # Ra Formula
  Ra_val <- (24 * 60 / pi) * 0.0820 * dr * (
    ws * sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * sin(ws)
  )
  ra_list[[m]] <- ifel(Ra_val < 0, 0, Ra_val)
}

# Convert to 12-layer stack
ra_clim <- rast(ra_list)

# Check Ra values
log_event("Ra climatology statistics:")
for (m in 1:12) {
  ra_stats <- global(ra_clim[[m]], c("min", "max", "mean"), na.rm = TRUE)
  log_event(sprintf("  Month %2d: min=%.2f, max=%.2f, mean=%.2f MJ/m²/day", 
                    m, ra_stats$min, ra_stats$max, ra_stats$mean))
}

# Expand Ra to match full time series
Ra_stack <- ra_clim[[month_indices]]
log_event("Ra stack generated.")

# --- 5. VECTORIZED PET CALCULATION ---
log_event("Executing vectorized Penman-Monteith calculation...")

# A. Temperature & Pressure
log_event("... Temperature and Pressure")
T_mean <- t2m_b - 273.15      # Kelvin -> Celsius
P_kpa  <- sp_b / 1000         # Pa -> kPa
gamma  <- 0.000665 * P_kpa    # Psychrometric constant

# Check temperature range
temp_stats <- global(T_mean, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Temperature range: %.2f to %.2f °C (mean: %.2f)", 
                  min(temp_stats$min), max(temp_stats$max), mean(temp_stats$mean)))

# B. Vapor Pressure
log_event("... Vapor Pressure")
# Saturation Vapor Pressure (es)
es <- 0.6108 * exp((17.27 * T_mean) / (T_mean + 237.3))

# Actual Vapor Pressure (ea)
Tdew <- d2m_b - 273.15
Tdew <- ifel(is.na(Tdew), T_mean - 2, Tdew)
ea <- 0.6108 * exp((17.27 * Tdew) / (Tdew + 237.3))

# Check VPD
vpd_stats <- global(es - ea, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("VPD range: %.3f to %.3f kPa (mean: %.3f)", 
                  min(vpd_stats$min), max(vpd_stats$max), mean(vpd_stats$mean)))

# C. Wind Speed
log_event("... Wind Speed")
u10_speed <- sqrt(u10_b^2 + v10_b^2)
u2 <- u10_speed * 4.87 / log(67.8 * 10 - 5.42)
u2 <- ifel(is.na(u2), 2.0, u2)

wind_stats <- global(u2, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Wind speed (2m) range: %.2f to %.2f m/s (mean: %.2f)", 
                  min(wind_stats$min), max(wind_stats$max), mean(wind_stats$mean)))

# D. Radiation
log_event("... Net Radiation")
# Solar Radiation (Rs) - CORRECTED FOR MONTHLY ACCUMULATED DATA
# CRITICAL FIX: ERA5-Land monthly SSRD is ACCUMULATED over the month
# Must divide by days in month to get daily average for Penman-Monteith

log_event("... Converting SSRD from J/m²/month to MJ/m²/day")
log_event("  WARNING: ERA5-Land SSRD is monthly accumulated, not daily average")
log_event("  Applying conversion: (J/m²/month) / 1e6 / days_in_month")

# Convert: J/m²/month → MJ/m²/day (unit conversion AND temporal scaling)
Rs <- ssrd_b
for (i in 1:n_months) {
  Rs[[i]] <- (ssrd_b[[i]] / 1e6) / days_in_month[i]
}

# Check if conversion makes sense
rs_stats <- global(Rs, c("min", "max", "mean"), na.rm=TRUE)

log_event("SSRD conversion check:")
log_event(sprintf("  Rs range: %.2f to %.2f MJ/m²/day (mean: %.2f)",
                  min(rs_stats$min), max(rs_stats$max), mean(rs_stats$mean)))

# Expected daily Rs values for Nechako Basin latitude (~54°N):
# Winter: 2-5 MJ/m²/day, Summer: 20-30 MJ/m²/day, Annual mean: 12-15 MJ/m²/day
expected_mean <- mean(rs_stats$mean)
if (expected_mean < 8 || expected_mean > 20) {
  log_event("  ⚠️ WARNING: Rs annual mean outside expected range (8-20 MJ/m²/day)!")
} else {
  log_event("  ✓ Rs values in expected range for Nechako Basin")
}

# Check seasonal pattern
jan_rs <- global(Rs[[1]], "mean", na.rm=TRUE)[1,1]
jul_rs <- global(Rs[[7]], "mean", na.rm=TRUE)[1,1]
log_event(sprintf("  January Rs: %.2f, July Rs: %.2f MJ/m²/day", jan_rs, jul_rs))

if (jan_rs < 2 || jan_rs > 8) {
  log_event("    ⚠️ January Rs outside expected 2-8 range")
}
if (jul_rs < 18 || jul_rs > 35) {
  log_event("    ⚠️ July Rs outside expected 18-35 range")
}

# Handle missing data with Hargreaves estimate
Rs <- ifel(is.na(Rs), K_RS * sqrt(10) * Ra_stack, Rs)

# Check Rs
rs_stats <- global(Rs, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Rs range: %.2f to %.2f MJ/m²/day (mean: %.2f)", 
                  min(rs_stats$min), max(rs_stats$max), mean(rs_stats$mean)))

# Clear-sky Radiation (Rso)
Rso <- (0.75 + 2e-5 * h_m) * Ra_stack

# Net Shortwave (Rns)
Rns <- (1 - 0.23) * Rs

# Net Longwave (Rnl)
Rs_Rso <- Rs / Rso
Rs_Rso <- clamp(Rs_Rso, lower=0.3, upper=1.0)
f_cloud <- 1.35 * Rs_Rso - 0.35

sigma <- 4.903e-9
T_kelvin_4 <- (T_mean + 273.15)^4
Rnl <- sigma * T_kelvin_4 * (0.34 - 0.14 * sqrt(ea)) * f_cloud

# Net Radiation (Rn)
Rn <- Rns - Rnl

# Check Rn before constraining
rn_stats <- global(Rn, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Rn range (before constraint): %.2f to %.2f MJ/m²/day (mean: %.2f)", 
                  min(rn_stats$min), max(rn_stats$max), mean(rn_stats$mean)))

# Count negative Rn values
negative_rn <- global(Rn < 0, "sum", na.rm = TRUE)
log_event(sprintf("Layers with negative Rn pixels: %d of %d", 
                  sum(negative_rn > 0), nlyr(Rn)))

# E. Soil Heat Flux (G) - CORRECTED FOR MONTHLY TIMESTEP
log_event("... Soil Heat Flux (G)")
# CORRECTED: Use 0.07 coefficient for MONTHLY (not 0.14 which is for WEEKLY)
# FAO-56 Chapter 4, p.89: G_month = 0.07 × (T_month_i - T_month_i-1)
if (nlyr(T_mean) > 1) {
  # Create lagged temperature: first month gets 0, rest get previous month
  t_prev <- c(T_mean[[1]] * 0, T_mean[[1:(nlyr(T_mean)-1)]])
  G <- 0.14 * (T_mean - t_prev)  # CHANGED FROM 0.14 TO 0.07
  
  # Reset G to 0 for January of each year (no previous month in same year)
  january_indices <- which(month_indices == 1)
  G[[1]] <- 0   # Jan 1950 only: no prior month available
  log_event("  ✓ G reset to 0 for January 1950")
} else {
  G <- T_mean * 0
}

g_stats <- global(G, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("G range: %.2f to %.2f MJ/m²/day (mean: %.2f)", 
                  min(g_stats$min), max(g_stats$max), mean(g_stats$mean)))

# F. Final ET0 Calculation
log_event("... Final ET0 calculation")

delta_slope <- (4098 * es) / (T_mean + 237.3)^2

num <- 0.408 * delta_slope * (Rn - G) + gamma * (900 / (T_mean + 273)) * u2 * (es - ea)
den <- delta_slope + gamma * (1 + 0.34 * u2)

et0 <- num / den

# Check for negative values BEFORE clipping
negative_mask <- et0 < 0
negative_count <- global(negative_mask, "sum", na.rm = TRUE)
total_pixels <- ncell(et0[[1]])

log_event("=== NEGATIVE ET0 DIAGNOSTIC ===")
for (i in 1:min(12, nlyr(et0))) {
  pct_negative <- (negative_count[i,1] / total_pixels) * 100
  log_event(sprintf("  Layer %d (%s): %.1f%% pixels negative", 
                    i, names(et0)[i], pct_negative))
}

# Count layers with >5% negative values
problematic_layers <- sum(negative_count / total_pixels > 0.05)
log_event(sprintf("Layers with >5%% negative pixels: %d of %d", 
                  problematic_layers, nlyr(et0)))

if (problematic_layers > nlyr(et0) * 0.1) {
  log_event("WARNING: >10% of layers have significant negative values - check calculation!")
}

# Physical constraints
et0 <- ifel(et0 < 0, 0, et0)

# Set names
names(et0) <- sprintf("PET_%04d_%02d", as.numeric(format(dates, "%Y")), 
                      as.numeric(format(dates, "%m")))

# Final statistics after clipping
et0_stats <- global(et0, c("min", "max", "mean"), na.rm = TRUE)
log_event(sprintf("Final ET0 range: %.3f to %.3f mm/day (mean: %.3f)", 
                  min(et0_stats$min), max(et0_stats$max), mean(et0_stats$mean)))

# --- 6. STATISTICS & OUTPUT ---
log_event("Calculation complete. Generating statistics and writing files...")

# Calculate Global Stats
global_stats <- global(et0, c("mean", "min", "max", "sd"), na.rm=TRUE)

pet_summary <- data.frame(
  date = format(dates, "%Y-%m"),
  year = as.numeric(format(dates, "%Y")),
  month = as.numeric(format(dates, "%m")),
  mean_pet = global_stats$mean,
  min_pet = global_stats$min,
  max_pet = global_stats$max,
  std_dev = global_stats$sd
)

# NetCDF Output
log_event("Writing NetCDF...")
varnames(et0) <- "pet"
longnames(et0) <- "Reference Evapotranspiration (FAO-56 Penman-Monteith)"
units(et0) <- "mm/day"
nc_file <- "potential_evapotranspiration_monthly.nc"
writeCDF(et0, nc_file, varname = "pet", overwrite = TRUE, 
         unit = "mm/day", 
         longname = "Reference Evapotranspiration (FAO-56 Penman-Monteith)")

# CSV Output
csv_file <- "ERA5Land_Nechako_PET_monthly_summary.csv"
write.csv(pet_summary, csv_file, row.names = FALSE)

# Calendar Month Stats
log_event("Calculating calendar month climatology...")
monthly_stats <- aggregate(pet_summary[, c("mean_pet", "min_pet", "max_pet", "std_dev")], 
                           by = list(month = pet_summary$month), 
                           FUN = mean, na.rm = TRUE)
monthly_stats$month_name <- month.name[monthly_stats$month]
write.csv(monthly_stats, "ERA5Land_Nechako_PET_calendar_month_stats.csv", row.names=FALSE)

# --- 7. DIAGNOSTIC PLOTS & VALIDATION ---
log_event("Generating diagnostic plots...")

# Close any existing plot devices to prevent sharing violations
graphics.off() 

# Create diagnostics directory
diag_dir <- "diagnostics"
if (!dir.exists(diag_dir)) dir.create(diag_dir, showWarnings = FALSE)

# Set up PDF for all diagnostic plots
pdf_path <- file.path(diag_dir, "PET_diagnostics.pdf")

# Attempt to open PDF - if this fails, the file is likely open in Foxit
tryCatch({
  pdf(pdf_path, width = 12, height = 8)
}, error = function(e) {
  stop("CANNOT WRITE PDF: Please close 'PET_diagnostics.pdf' in your PDF reader and run again.")
})

# A. TIME SERIES PLOT
log_event("... Time series plot")
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(dates, pet_summary$mean_pet, type = "l", col = "blue", lwd = 2,
     xlab = "Date", ylab = "PET (mm/day)",
     main = "Basin-Average PET Time Series")
grid()
if (nrow(pet_summary) >= 12) {
  ma_12 <- stats::filter(pet_summary$mean_pet, rep(1/12, 12), sides = 2)
  lines(dates, ma_12, col = "red", lwd = 2)
  legend("topleft", legend = c("Monthly PET", "12-month MA"), 
         col = c("blue", "red"), lwd = 2, bty = "n")
}

# B. SEASONAL CYCLE
log_event("... Seasonal cycle")
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
boxplot(mean_pet ~ month, data = pet_summary,
        xlab = "Month", ylab = "PET (mm/day)",
        main = "Seasonal Distribution of PET",
        names = month.abb,
        col = "lightblue", border = "darkblue")
grid()

# C. EXPECTED RANGE CHECK
log_event("... Range validation")
lower_bound <- 0
upper_bound <- 8
anomalies <- pet_summary[pet_summary$mean_pet < lower_bound | 
                           pet_summary$mean_pet > upper_bound, ]

if (nrow(anomalies) > 0) {
  log_event(paste("WARNING:", nrow(anomalies), "months with PET outside expected range"))
  write.csv(anomalies, "diagnostics/PET_anomalies.csv", row.names = FALSE)
} else {
  log_event("All PET values within expected range (0-8 mm/day)")
}

# D. MONTHLY CLIMATOLOGY
log_event("... Monthly climatology comparison")
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(monthly_stats$month, monthly_stats$mean_pet, type = "b", pch = 19, col = "darkblue",
     xlab = "Month", ylab = "Mean PET (mm/day)",
     main = "Mean PET Climatology by Calendar Month",
     xaxt = "n", ylim = c(0, max(monthly_stats$mean_pet, na.rm = TRUE) * 1.1))
axis(1, at = 1:12, labels = month.abb)
grid()

arrows(monthly_stats$month, 
       monthly_stats$mean_pet - monthly_stats$std_dev,
       monthly_stats$month, 
       monthly_stats$mean_pet + monthly_stats$std_dev,
       angle = 90, code = 3, length = 0.05, col = "gray50")

# E. INTER-ANNUAL VARIABILITY (WITH BOTH LINEAR REGRESSION AND SEN'S SLOPE)
log_event("... Inter-annual variability")
annual_mean <- aggregate(mean_pet ~ year, data = pet_summary, FUN = mean, na.rm = TRUE)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
plot(annual_mean$year, annual_mean$mean_pet, type = "b", pch = 19, col = "darkgreen",
     xlab = "Year", ylab = "Annual Mean PET (mm/day)",
     main = "Annual Mean PET Trend")
grid()

if (nrow(annual_mean) > 2) {
  # Linear Regression Trend
  trend_lm <- lm(mean_pet ~ year, data = annual_mean)
  abline(trend_lm, col = "red", lwd = 2, lty = 2)
  trend_slope_lm <- coef(trend_lm)[2]
  trend_pval_lm <- summary(trend_lm)$coefficients[2, 4]
  
  # Sen's Slope Estimator
  log_event("... Calculating Sen's slope")
  n <- nrow(annual_mean)
  slopes <- c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      slope_ij <- (annual_mean$mean_pet[j] - annual_mean$mean_pet[i]) / 
        (annual_mean$year[j] - annual_mean$year[i])
      slopes <- c(slopes, slope_ij)
    }
  }
  sen_slope <- median(slopes, na.rm = TRUE)
  
  # Draw Sen's slope line
  sen_intercept <- median(annual_mean$mean_pet - sen_slope * annual_mean$year, na.rm = TRUE)
  abline(a = sen_intercept, b = sen_slope, col = "blue", lwd = 2, lty = 1)
  
  # Mann-Kendall test for significance (simple implementation)
  S <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      S <- S + sign(annual_mean$mean_pet[j] - annual_mean$mean_pet[i])
    }
  }
  var_S <- n * (n - 1) * (2 * n + 5) / 18
  Z <- ifelse(S > 0, (S - 1) / sqrt(var_S), ifelse(S < 0, (S + 1) / sqrt(var_S), 0))
  mk_pval <- 2 * pnorm(-abs(Z))
  
  log_event(sprintf("Sen's slope: %.6f mm/day/year (MK p=%.3f)", sen_slope, mk_pval))
  
  legend("topleft", 
         legend = c(
           sprintf("Linear Regression: %.4f mm/day/year (p=%.3f)", 
                   trend_slope_lm, trend_pval_lm),
           sprintf("Sen's Slope: %.4f mm/day/year (p=%.3f)", 
                   sen_slope, mk_pval)
         ),
         col = c("red", "blue"),
         lwd = 2,
         lty = c(2, 1),
         bty = "n")
}
# F. SAMPLE SPATIAL MAPS (WITH CORRECTED COLOR PALETTE)
log_event("... Sample monthly maps")
sample_months <- c(1, 256, 619, 910)  # Jan, Apr, Jul, Oct of some years
if (nlyr(et0) >= 12) {
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))
  for (idx in sample_months) {
    # REVISED color palette: red for high values, yellow for low values
    plot(et0[[idx]], col = rev(hcl.colors(100, "YlOrRd")),
         main = sprintf("PET %s", names(et0)[idx]))
  }
}

# G. DISTRIBUTION HISTOGRAM
log_event("... PET distribution histogram")
par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
hist(pet_summary$mean_pet, breaks = 30, col = "skyblue", border = "white",
     xlab = "PET (mm/day)", main = "Distribution of Basin-Mean PET",
     freq = FALSE)
lines(density(pet_summary$mean_pet, na.rm = TRUE), col = "red", lwd = 2)
mean_pet <- mean(pet_summary$mean_pet, na.rm = TRUE)
sd_pet <- sd(pet_summary$mean_pet, na.rm = TRUE)
curve(dnorm(x, mean_pet, sd_pet), add = TRUE, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("Actual", "Normal fit"), 
       col = c("red", "blue"), lwd = 2, lty = c(1, 2), bty = "n")

# H. QUALITY CONTROL SUMMARY
log_event("... Quality control summary")
qc_summary <- data.frame(
  Metric = c("Total Months", "Mean PET (mm/day)", "Median PET (mm/day)", 
             "Min PET (mm/day)", "Max PET (mm/day)", "Std Dev (mm/day)",
             "% Months with PET < 0", "% Months with PET > 8",
             "% Missing Values"),
  Value = c(
    nrow(pet_summary),
    round(mean(pet_summary$mean_pet, na.rm = TRUE), 3),
    round(median(pet_summary$mean_pet, na.rm = TRUE), 3),
    round(min(pet_summary$mean_pet, na.rm = TRUE), 3),
    round(max(pet_summary$mean_pet, na.rm = TRUE), 3),
    round(sd(pet_summary$mean_pet, na.rm = TRUE), 3),
    round(sum(pet_summary$mean_pet < 0, na.rm = TRUE) / nrow(pet_summary) * 100, 2),
    round(sum(pet_summary$mean_pet > 8, na.rm = TRUE) / nrow(pet_summary) * 100, 2),
    round(sum(is.na(pet_summary$mean_pet)) / nrow(pet_summary) * 100, 2)
  )
)

par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
     xlim = c(0, 1), ylim = c(0, 1))
title("Quality Control Summary", cex.main = 1.5)

y_start <- 0.85
y_step <- 0.08
for (i in 1:nrow(qc_summary)) {
  text(0.05, y_start - (i-1) * y_step, qc_summary$Metric[i], 
       adj = 0, cex = 1.1, font = 2)
  text(0.65, y_start - (i-1) * y_step, qc_summary$Value[i], 
       adj = 0, cex = 1.1)
}
lines(c(0.05, 0.95), c(y_start + 0.02, y_start + 0.02), lwd = 2)

# CLOSE THE PDF DEVICE (Critical to prevent sharing violations)
dev.off()
log_event("Main diagnostic PDF closed and saved.")

log_event("Diagnostics complete. Check 'diagnostics/' folder.")
log_event("Main diagnostic PDF: diagnostics/PET_diagnostics.pdf")

# --- 8. FINISH ---
log_event("==========================================")
log_event(paste("Total months processed:", n_months))
log_event(paste("Overall mean PET:", round(mean(pet_summary$mean_pet, na.rm=TRUE), 2), "mm/day"))
log_event(paste("Months with >5% negative pixels:", problematic_layers))
log_event("PROCESSING COMPLETE")
log_event(paste("Processing finished:", Sys.time()))