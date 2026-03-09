# ============================================================================
#  Download ERA5 Atmospheric Circulation + SST Variables
#
#  Variables:
#    [1] 500 hPa Geopotential Height  (pressure-levels, AREA_ATM)
#    [2] Mean Sea Level Pressure      (single-levels,   AREA_ATM)
#    [3] Sea Surface Temperature      (single-levels,   AREA_SST)
#
#  Anomalies computed for [1] and [3]:
#    (a) Full-period  : value minus climatology over ALL downloaded years
#    (b) Fixed-period : value minus 1991-2020 WMO standard climatology
#
#  Domains:
#    AREA_ATM  N=70, W=-180, S=40,  E=-100  (500 hPa + SLP)
#    AREA_SST  N=65, W=-180, S=20,  E=-110  (SST NE Pacific)
#
#  Justified by:
#    - TRARE ERA5 synoptic map extents (Hurley et al. 2025)
#    - Richards-Thomas et al. (2024) ERA5/ERA5-Land usage
#    - PDO/Marine Heatwave domain literature
# ============================================================================

# ── Libraries ────────────────────────────────────────────────────────────────
library(ecmwfr)
library(terra)
library(lubridate)
library(ncdf4)

# ── Working directory ────────────────────────────────────────────────────────
setwd("D:/Nechako_Drought/")

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

# ── CDS credentials ──────────────────────────────────────────────────────────
CDS_USER <- "mehdi.bateni@iusspavia.it"   # Numeric UID from cds.climate.copernicus.eu
CDS_KEY  <- "44ace8e1-8ec8-4246-a44e-f30409c90c5b" 

wf_set_key(user = CDS_USER, key = CDS_KEY)

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

# ── Output directories ────────────────────────────────────────────────────────
dir.create("monthly_data_direct", showWarnings = FALSE)
dir.create("downloads",           showWarnings = FALSE)

# ============================================================================
#  FUNCTION: BUILD AND SUBMIT CDS REQUEST
# ============================================================================
download_era5_variable <- function(variable_info, start_yr, end_yr) {
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
    nc_file <- download_era5_variable(var_info, start_year, end_year)
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