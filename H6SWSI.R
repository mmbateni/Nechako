# ============================================================================
# SURFACE-WATER SUPPLY INDEX (SWSI) — NECHAKO BASIN
# DUAL-TRACK IMPLEMENTATION  |  Garen (1993) Revised — Z-score formulation
#
# ┌──────────────────────────────────────────────────────────────────────────┐
# │ TRACK 1 — H-SWSI  Hydrological Drought Index                            │
# │   Objective : Isolate climate-driven anomalies (unregulated basin)       │
# │   Baseline  : Pre-regulation naturalized record (pre-1952 preferred)     │
# │   Components: Qnat · LakeStorage · Precipitation · SWE                  │
# │   Weights   : Equal (0.25 each) or proportional to hydroclimate input   │
# ├──────────────────────────────────────────────────────────────────────────┤
# │ TRACK 2 — O-SWSI  Operational Supply Drought Index                      │
# │   Objective : Measure actual system stress & management vulnerability    │
# │   Baseline  : Post-regulation (1984+ constrained by GloLakes start)     │
# │   Components: Active Reservoir Storage · SWE · Unregulated Inflow (Qin) │
# │   Weights   : Heavily toward Storage (0.50) + SWE (0.35)                │
# └──────────────────────────────────────────────────────────────────────────┘
#
# Mathematical Core (Garen 1993 Revised):
#   Pi  = non-exceedance probability [%] from historical baseline CDF
#   Zi  = Φ⁻¹(Pi/100)   — standard normal variate (Z-score)
#   SWSI = Σ(Wi · Zi)    — weighted sum; capped at ±4.2
#
# Station Notes:
#   08JC001  Nechako R at Vanderhoof     | regulated + naturalized  (Qnat & Qin)
#   08JC002  Nechako R at Isle Pierre   | regulated + naturalized  (Qnat & Qin PRIMARY)
#   08JB002  STELLAKO RIVER AT GLENANNAN    | unaffected by regulation 1929-2024 
#   08JA015  Stuart Lake outlet          | unregulated 1976-2023    [supplementary Qnat]
#
# Natural Lakes for H-SWSI LakeNat component (GloLakes v2):
#   Stuart Lake  ~825 km²  54.9°N 124.1°W  — Nechako basin
#   Fraser Lake  ~176 km²  54.0°N 124.8°W  — Nechako basin
#
# References:
#   Garen, D.C. (1993). Revised Surface-Water Supply Index for Western
#     United States. J. Water Resour. Plan. Manage. 119(4):437-454.
#     https://doi.org/10.1061/(ASCE)0733-9496(1993)119:4(437)
#   Hou, J. et al. (2024). GloLakes v2. Earth Syst. Sci. Data.
#
# Data Sources:
#   (1) GloLakes v2      — Lake / reservoir storage 1984-present
#   (2) ERA5-Land        — SWE & Precipitation (monthly_data_direct/)
#   (3) WSC              — Naturalized & regulated discharge
#   (4) BC RFC           — Snow pillow SWE (Hydrology/SWE_data/)
# ============================================================================

rm(list = ls())

# ============================================================================
# SECTION 1: PACKAGES
# ============================================================================
packages_needed <- c("tidyverse", "lubridate", "zoo", "ncdf4",
                     "lmomco", "terra", "scales", "patchwork")
for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# SECTION 2: FILE PATHS & SETTINGS
# ============================================================================
setwd("D:/Nechako_Drought/Nechako")

# ── Input Directories ────────────────────────────────────────────────────────
INPUT_GLOLAKES_DIR   <- "Lakes/GloLakes/nechako_extracted"
INPUT_ERA5_DIR       <- "monthly_data_direct"
INPUT_NATURAL_DIR    <- "Hydrology/data_retrievalWaterSurveyofCanada/naturalized_flows"
INPUT_WSC_DIR        <- "Hydrology/data_retrievalWaterSurveyofCanada/data_downloads_geomet_api"
INPUT_SNOWPILLOW_DIR <- "Hydrology/SWE_data"

# ── Output Directories ───────────────────────────────────────────────────────
MAIN_OUTPUT_DIR <- "swsi_results_dual_track"
for (sub in c("", "processed_data", "era5_extracted", "diagnostics",
              "hswsi", "hswsi/components", "hswsi/droughts",
              "oswsi", "oswsi/components", "oswsi/droughts",
              "comparison")) {
  dir.create(file.path(MAIN_OUTPUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

# ── Basin / Domain ───────────────────────────────────────────────────────────
BASIN_KMZ_PATH   <- "Spatial/nechakoBound_dissolve.kmz"
BASIN_LAT_MIN    <- 52.5
BASIN_LAT_MAX    <- 55.5
BASIN_LON_MIN    <- -127.0
BASIN_LON_MAX    <- -120.5
BASIN_AREA_KM2   <- 47000       # Nechako above Isle Pierre

# ── Analysis Parameters ──────────────────────────────────────────────────────
WATER_YEAR_START_MONTH <- 10L
ACCUMULATION_MONTHS    <- c(10, 11, 12, 1, 2, 3)   # Oct–Mar
MELT_MONTHS            <- c(4, 5, 6, 7, 8, 9)       # Apr–Sep
MIN_YEARS_FOR_DIST     <- 10    # minimum baseline years for CDF fitting
MIN_COMPLETENESS_PCT   <- 70
MAX_FILL_DAYS          <- 14
MIN_DAYS_PER_MONTH     <- 15

# ── SWSI Mathematical Scale ──────────────────────────────────────────────────
# Garen (1993): ±4.2 corresponds to non-exceedance probabilities of ~0.001/99.999
SWSI_CAP          <- 4.2

# ── Drought Triggers (Garen 1993, p.444) ────────────────────────────────────
TRIGGER_WATCH     <- -1.0
TRIGGER_EMERGENCY <- -2.0

# ── Regulation Epoch ─────────────────────────────────────────────────────────
# Kenney Dam closed October 1952; inter-basin diversion via Kemano tunnel began.
REGULATION_START  <- as.Date("1952-10-01")

# ============================================================================
# TRACK 1 — H-SWSI PARAMETERS
# ============================================================================
# Baseline: pre-regulation naturalized data (pre-1952) is preferred.
# Falls back to the full naturalized record when <MIN_YEARS_FOR_DIST pre-1952
# months are available (e.g. for GloLakes lake storage, which starts 1984,
# the entire record is used as the reference distribution).
HSWSI_BASELINE_END <- REGULATION_START

# Component weights — must sum to 1.
# Default: equal weighting (0.25 each).
# Adjust proportionally if area-weighted hydroclimate contributions are known.
HSWSI_WEIGHTS <- c(
  Qnat    = 0.25,   # Naturalized streamflow   (WSC 08JC002 primary)
  LakeNat = 0.25,   # Aggregate natural lake storage (GloLakes)
  Precip  = 0.25,   # Basin-averaged YTD precipitation (ERA5)
  SWE     = 0.25    # Snow water equivalent, bias-corrected (ERA5 + BC RFC)
)

# Stations used for H-SWSI Qnat component.
# 08JC002 (primary): best naturalized record, unaffected by Nautley confluence.
# 08JA015 (supplementary): Stuart Lake outlet — naturally unregulated 1976-2023.
HSWSI_QNAT_PRIMARY     <- "08JC002"
HSWSI_QNAT_SUPPLEMENT  <- "08JA015"   # loaded if available; not blended automatically

# ============================================================================
# TRACK 2 — O-SWSI PARAMETERS
# ============================================================================
# Baseline: post-regulation data only.
# 1984-01-01 is used as the common start (GloLakes constraint).
# Using a common baseline start across all three components ensures the
# standardisation reflects the post-Kemano-diversion operating regime.
OSWSI_BASELINE_START <- as.Date("1984-01-01")

# Component weights — must sum to 1.
# Storage and SWE are the primary buffers; Qin is secondary.
OSWSI_WEIGHTS <- c(
  Storage = 0.50,   # Reservoir active storage (Kenney/Nechako Reservoir)
  SWE     = 0.35,   # Snowpack — primary spring refill signal
  Qin     = 0.15    # Unregulated inflow proxy (08JC002 naturalized)
)

# Station used as proxy for unregulated reservoir inflow (Qin).
# 08JC002 naturalized flow represents the hydrological input that would
# have reached the reservoir under natural conditions.
OSWSI_QIN_STATION <- "08JC002"

# ── Reservoir Storage Adjustments (O-SWSI) ───────────────────────────────────
# Active storage = Total GloLakes storage
#                − Dead storage (volume below minimum operating level)
#                − Environmental flow reserve (legal obligation)
#
# *** UPDATE THESE FROM BC HYDRO / RIO TINTO OPERATING LICENCE ***
KENNEY_DEAD_STORAGE_MM3  <- 2500    # [Mm³]  Minimum operating level (dead storage)
ENV_FLOW_RESERVE_MM3     <- 100     # [Mm³]  Legally reserved environmental flow volume

# Kemano inter-basin diversion subtraction from Qin.
# Set to 0 (default) to retain raw naturalized flow as Qin without adjustment.
# Set to the long-term monthly-mean diversion volume to isolate purely local supply.
KEMANO_BASELINE_DIV_MM3_MONTH <- 0  # [Mm³/month]  0 = no adjustment

# ============================================================================
# NATURAL LAKE CATALOGUE  (H-SWSI LakeNat component)
# ============================================================================
NATURAL_LAKES <- data.frame(
  name             = c("Stuart Lake",  "Fraser Lake"),
  area_km2         = c(825L,           176L),
  lat              = c(54.9,           54.0),
  lon              = c(-124.1,         -124.8),
  in_nechako_basin = c(TRUE,           TRUE),
  stringsAsFactors = FALSE
)

# ── Optional: Explicit GloLakes file-path overrides ──────────────────────────
# Leave as NULL to use auto-detection (name search inside INPUT_GLOLAKES_DIR).
# Example: GLOLAKES_LAKE_FILES[["Stuart Lake"]] <- "Lakes/GloLakes/.../lake_xyz.csv"
GLOLAKES_LAKE_FILES <- list(
  "Stuart Lake"      = NULL,
  "Fraser Lake"      = NULL,
  "Kenney Reservoir" = NULL    # NULL → uses BEST_ESTIMATE_basin_total_Mm3.csv
)

# ── Station Routing Table ────────────────────────────────────────────────────
STATIONS <- data.frame(
  StationID       = c("08JC001",                    "08JC002",
                      "08JB002",                    "08JA015"),
  Name            = c("Nechako R at Vanderhoof",    "Nechako R above Nautley R",
                      "Nechako R at Isle Pierre",   "Stuart Lake outlet"),
  Has_NatFlow     = c(TRUE,                          TRUE,  FALSE, FALSE),
  Use_HSWSI_Qnat  = c(TRUE,                          TRUE,  FALSE, TRUE),
  Use_OSWSI_Qin   = c(FALSE,                         TRUE,  FALSE, FALSE),
  Is_Regulated    = c(TRUE,                          TRUE,  TRUE,  FALSE),
  Notes           = c("Regulated + naturalized",
                      "Regulated + naturalized; primary station",
                      "Regulated ONLY 1929-2024; downstream of reservoir",
                      "Unregulated 1976-2023; Stuart Lake outlet"),
  stringsAsFactors = FALSE
)

# ── Plot Focus Window ─────────────────────────────────────────────────────────
FOCUS_START <- as.Date("2022-01-01")
FOCUS_END   <- as.Date("2025-12-31")

# ============================================================================
# SECTION 3: BASIN BOUNDARY — Nechako polygon from KMZ
# ============================================================================
load_basin_boundary <- function(kmz_path = BASIN_KMZ_PATH) {
  cat("\n[Basin] Loading Nechako boundary from KMZ...\n")
  
  if (!file.exists(kmz_path))
    stop("Basin KMZ not found: ", kmz_path,
         "\n  Expected at: Spatial/nechakoBound_dissolve.kmz")
  
  tmp <- tempfile()
  dir.create(tmp, showWarnings = FALSE)
  utils::unzip(kmz_path, exdir = tmp)
  
  kml_files <- list.files(tmp, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)
  if (length(kml_files) == 0) stop("No .kml found inside ", kmz_path)
  
  v <- terra::vect(kml_files[1])
  if (nrow(v) > 1L) v <- terra::aggregate(v)
  unlink(tmp, recursive = TRUE)
  
  area_km2 <- as.numeric(terra::expanse(terra::project(v, "EPSG:3005"), unit = "km"))
  cat(sprintf("  Polygon area: %.0f km²  |  BASIN_AREA_KM2 constant: %.0f km²\n",
              area_km2, BASIN_AREA_KM2))
  return(v)
}

# ============================================================================
# SECTION 4: GLOLAKES STORAGE — RESERVOIR (O-SWSI) + NATURAL LAKES (H-SWSI)
# ============================================================================

# ── 4a. Reservoir / basin-total storage (O-SWSI) ────────────────────────────
load_reservoir_storage <- function() {
  cat("\n[GloLakes | O-SWSI] Loading reservoir storage (best-estimate)...\n")
  
  # Allow explicit override path for Kenney Reservoir
  f <- if (!is.null(GLOLAKES_LAKE_FILES[["Kenney Reservoir"]]) &&
           file.exists(GLOLAKES_LAKE_FILES[["Kenney Reservoir"]]))
    GLOLAKES_LAKE_FILES[["Kenney Reservoir"]]
  else
    file.path(INPUT_GLOLAKES_DIR, "BEST_ESTIMATE_basin_total_Mm3.csv")
  
  if (!file.exists(f))
    stop("GloLakes best-estimate file not found: ", f,
         "\nRun GloLakes_Nechako_Extraction_v2.R first.")
  
  raw <- read.csv(f, stringsAsFactors = FALSE)
  
  reservoir <- raw %>%
    mutate(
      date       = as.Date(date),
      storage    = as.numeric(total_storage_Mm3),
      year       = year(date),
      month      = month(date),
      water_year = if_else(month >= WATER_YEAR_START_MONTH, year, year - 1L),
      n_lakes    = as.numeric(n_lakes_contributing),
      products   = product
    ) %>%
    filter(!is.na(storage), storage >= 0, n_lakes >= 1) %>%
    dplyr::select(date, storage, year, month, water_year, n_lakes, products)
  
  cat(sprintf("  Loaded: %d months  (%s to %s)\n",
              nrow(reservoir), min(reservoir$date), max(reservoir$date)))
  cat(sprintf("  Storage range: %.0f – %.0f Mm³  |  mean %.0f Mm³\n",
              min(reservoir$storage), max(reservoir$storage),
              mean(reservoir$storage, na.rm = TRUE)))
  
  if (nrow(reservoir) > 12) {
    p <- ggplot(reservoir, aes(x = date, y = storage)) +
      geom_line(colour = "#1f78b4", linewidth = 0.8) +
      labs(title    = "Nechako Basin — GloLakes Reservoir Storage",
           subtitle = "Source: GloLakes v2 (Hou et al. 2024) | O-SWSI component",
           x = NULL, y = "Storage (Mm³)") +
      theme_bw()
    ggsave(file.path(MAIN_OUTPUT_DIR, "diagnostics", "plot_reservoir_storage.png"),
           p, width = 12, height = 5, dpi = 300)
  }
  return(reservoir)
}

# ── 4b. Individual natural lake storage (H-SWSI LakeNat component) ──────────
#
# Searches INPUT_GLOLAKES_DIR for per-lake CSV files matching each lake name.
# Auto-detection uses partial name matching (e.g. "Stuart" for "Stuart Lake").
# Override auto-detection via GLOLAKES_LAKE_FILES[["Stuart Lake"]] = "path/to/file.csv".
#
# Expected CSV schema (column names are case-insensitive):
#   date        — YYYY-MM-DD or parseable date string
#   storage     — lake volume in Mm³   (or 'volume', 'level')
#
load_natural_lake_storage <- function(lake_catalogue = NATURAL_LAKES,
                                      glolakes_dir   = INPUT_GLOLAKES_DIR,
                                      explicit_files = GLOLAKES_LAKE_FILES) {
  cat("\n[GloLakes | H-SWSI] Loading individual natural lake storage...\n")
  
  all_csv <- list.files(glolakes_dir, pattern = "\\.csv$",
                        full.names = TRUE, recursive = FALSE)
  
  lake_results <- list()
  
  for (i in seq_len(nrow(lake_catalogue))) {
    lk      <- lake_catalogue[i, ]
    lk_name <- lk$name
    
    # Basin membership flag
    basin_flag <- if (lk$in_nechako_basin) "Nechako basin" else
      "⚠ Skeena drainage — regional indicator"
    cat(sprintf("  %-14s  (%.1f°N, %.1f°W, %.0f km²)  [%s]\n",
                lk_name, lk$lat, abs(lk$lon), lk$area_km2, basin_flag))
    
    # Resolve file path: explicit override > name match
    fname <- explicit_files[[lk_name]]
    if (!is.null(fname) && file.exists(fname)) {
      cat(sprintf("    Explicit path: %s\n", basename(fname)))
    } else {
      key        <- sub(" Lake$| Lake$", "", lk_name, ignore.case = TRUE)
      candidates <- all_csv[
        grepl(key, basename(all_csv), ignore.case = TRUE) &
          !grepl("basin_total|BEST_ESTIMATE|reservoir", basename(all_csv),
                 ignore.case = TRUE)
      ]
      if (length(candidates) > 0) {
        fname <- candidates[1]
        cat(sprintf("    Auto-detected: %s\n", basename(fname)))
      } else {
        cat(sprintf("    Not found in %s — lake will be absent from LakeNat\n",
                    basename(glolakes_dir)))
        next
      }
    }
    
    raw <- tryCatch(read.csv(fname, stringsAsFactors = FALSE),
                    error = function(e) { warning(e$message); NULL })
    if (is.null(raw) || nrow(raw) == 0) next
    
    colnames(raw) <- tolower(trimws(colnames(raw)))
    date_col <- grep("^date",             colnames(raw), value = TRUE)[1]
    stor_col <- grep("storage|volume|level", colnames(raw), value = TRUE,
                     ignore.case = TRUE)[1]
    
    if (is.na(date_col) || is.na(stor_col)) {
      warning(sprintf("Cannot identify date/storage columns in %s", basename(fname)))
      next
    }
    
    df <- tibble(
      date        = as.Date(raw[[date_col]]),
      storage_Mm3 = suppressWarnings(as.numeric(raw[[stor_col]])),
      lake_name   = lk_name,
      area_km2    = lk$area_km2,
      in_nechako  = lk$in_nechako_basin
    ) %>%
      filter(!is.na(date), !is.na(storage_Mm3), storage_Mm3 >= 0) %>%
      mutate(year  = year(date),
             month = month(date),
             date  = as.Date(paste0(year, "-", sprintf("%02d", month), "-01"))) %>%
      group_by(date, year, month, lake_name, area_km2, in_nechako) %>%
      summarise(storage_Mm3 = mean(storage_Mm3, na.rm = TRUE), .groups = "drop")
    
    lake_results[[lk_name]] <- df
    cat(sprintf("    Loaded: %d months  (%s – %s)  |  %.1f – %.1f Mm³\n",
                nrow(df), min(df$date), max(df$date),
                min(df$storage_Mm3), max(df$storage_Mm3)))
  }
  
  if (length(lake_results) == 0) {
    warning("[H-SWSI] No individual natural lake files loaded. LakeNat component will be absent.")
    return(NULL)
  }
  
  bind_rows(lake_results)
}

# ── 4c. Aggregate individual lakes → single monthly basin indicator ──────────
aggregate_natural_lakes <- function(lake_df) {
  if (is.null(lake_df) || nrow(lake_df) == 0) return(NULL)
  
  lake_df %>%
    group_by(date,
             year  = year(date),
             month = month(date)) %>%
    summarise(
      nat_lake_total_Mm3 = sum(storage_Mm3, na.rm = TRUE),
      n_lakes            = n(),
      lakes_included     = paste(lake_name, collapse = " + "),
      .groups            = "drop"
    ) %>%
    mutate(water_year = if_else(month >= WATER_YEAR_START_MONTH, year, year - 1L))
}

# ============================================================================
# SECTION 5: ERA5-LAND CLIMATE DATA
# ============================================================================
ERA5_EPOCH <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC")

.era5_decode_dates <- function(r, nc_path) {
  t <- terra::time(r)
  if (!all(is.na(t))) return(as.Date(t))
  
  nc         <- ncdf4::nc_open(nc_path)
  time_raw   <- ncdf4::ncvar_get(nc, "valid_time")
  time_units <- ncdf4::ncatt_get(nc, "valid_time", "units")$value
  ncdf4::nc_close(nc)
  
  cat(sprintf("    time units: %s\n", time_units))
  cat(sprintf("    time_raw range: %.0f to %.0f\n", min(time_raw), max(time_raw)))
  
  if (grepl("seconds since 1970", time_units, ignore.case = TRUE))
    as.Date(as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC"))
  else if (grepl("hours since 1900", time_units, ignore.case = TRUE))
    as.Date(ERA5_EPOCH + as.difftime(time_raw, units = "hours"))
  else {
    cat(sprintf("    WARNING: unrecognised time units '%s' — assuming Unix seconds\n",
                time_units))
    as.Date(as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC"))
  }
}

load_era5_basin_mean <- function(nc_path, var_name, basin_vect = NULL) {
  cat(sprintf("  ERA5: reading '%s' from %s\n", var_name, basename(nc_path)))
  
  r <- terra::rast(nc_path)
  if (terra::ext(r)$xmax > 180) {
    r <- terra::rotate(r)
    cat(sprintf("    Longitude rotated 0-360 → -180/180\n"))
  }
  
  dates <- .era5_decode_dates(r, nc_path)
  dates <- as.Date(format(dates, "%Y-%m-01"))
  cat(sprintf("    Dates: %s to %s\n", min(dates), max(dates)))
  
  if (!is.null(basin_vect)) {
    basin_proj <- if (!terra::same.crs(basin_vect, terra::crs(r)))
      terra::project(basin_vect, terra::crs(r)) else basin_vect
    r <- terra::crop(r, basin_proj, snap = "out")
    r <- terra::mask(r, basin_proj, touches = TRUE)
    n_basin <- sum(!is.na(terra::values(r[[1]])))
    cat(sprintf("    Basin pixels after KMZ mask: %d\n", n_basin))
    if (n_basin == 0)
      stop(sprintf("No ERA5 cells overlap the basin polygon for '%s'", var_name))
  } else {
    cat("    WARNING: no basin polygon — using full raster extent\n")
  }
  
  vals       <- terra::values(r, mat = TRUE)
  vals[vals >= 9e36] <- NA_real_
  cell_coords <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  lat_w       <- cos(cell_coords[, 2L] * pi / 180)
  
  basin_mean <- vapply(seq_len(ncol(vals)), function(t) {
    v  <- vals[, t]; ok <- !is.na(v)
    if (sum(ok) == 0L) return(NA_real_)
    weighted.mean(v[ok], lat_w[ok])
  }, numeric(1L))
  
  out <- tibble(
    date       = dates,
    year       = year(dates),
    month      = month(dates),
    water_year = if_else(month(dates) >= WATER_YEAR_START_MONTH,
                         year(dates), year(dates) - 1L)
  )
  out[[var_name]] <- basin_mean
  cat(sprintf("    Non-NA values: %d / %d\n", sum(!is.na(basin_mean)), length(basin_mean)))
  return(out)
}

load_era5_all <- function(era5_dir = INPUT_ERA5_DIR, basin_vect = NULL) {
  cat("\n[ERA5-Land] Loading basin-averaged monthly climate data...\n")
  
  f_snowc <- file.path(era5_dir, "snow_cover_monthly.nc")
  f_sd    <- file.path(era5_dir, "snow_depth_water_equivalent_monthly.nc")
  f_tp    <- file.path(era5_dir, "total_precipitation_monthly.nc")
  
  for (f in c(f_snowc, f_sd, f_tp))
    if (!file.exists(f)) stop("ERA5 file missing: ", f)
  
  df_snowc <- load_era5_basin_mean(f_snowc, "snowc", basin_vect)
  df_sd    <- load_era5_basin_mean(f_sd,    "sd",    basin_vect)
  df_tp    <- load_era5_basin_mean(f_tp,    "tp",    basin_vect)
  
  cat(sprintf("  Date ranges — snowc: %s to %s | sd: %s to %s | tp: %s to %s\n",
              min(df_snowc$date), max(df_snowc$date),
              min(df_sd$date),    max(df_sd$date),
              min(df_tp$date),    max(df_tp$date)))
  
  era5 <- df_snowc %>%
    inner_join(df_sd %>% dplyr::select(date, sd), by = "date") %>%
    inner_join(df_tp %>% dplyr::select(date, tp), by = "date") %>%
    mutate(
      swe_mm         = sd * 1000,
      swe_Mm3        = sd * BASIN_AREA_KM2 * 1e6 / 1e6,
      precip_mm      = tp * 1000,
      precip_Mm3     = tp * BASIN_AREA_KM2 * 1e6 / 1e6,
      snow_cover_pct = snowc * 100
    ) %>%
    arrange(date) %>%
    group_by(water_year) %>%
    mutate(precip_ytd_Mm3 = cumsum(precip_Mm3),
           precip_ytd_mm  = cumsum(precip_mm)) %>%
    ungroup()
  
  cat(sprintf("  ERA5 merged: %d months (%s to %s)\n",
              nrow(era5), min(era5$date), max(era5$date)))
  return(era5)
}

# ============================================================================
# SECTION 6: SNOW PILLOW BIAS CORRECTION (ERA5 SWE)
# ============================================================================
load_all_snow_pillows <- function(pillow_dir = INPUT_SNOWPILLOW_DIR) {
  if (!dir.exists(pillow_dir)) {
    cat("  Snow pillow directory not found — skipping bias correction\n")
    return(NULL)
  }
  cat("\n[Snow Pillows] Loading BC RFC snow pillow SWE data...\n")
  
  all_files <- list.files(pillow_dir, full.names = TRUE, pattern = "\\.csv$")
  if (length(all_files) == 0) {
    cat("  No CSV files found — skipping\n"); return(NULL)
  }
  
  all_loaded <- list()
  for (f in all_files) {
    raw <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(raw) || nrow(raw) == 0) next
    colnames(raw) <- trimws(colnames(raw))
    col_start <- grep("start|Start|Date", colnames(raw), value = TRUE)[1]
    col_swe   <- grep("SW|swe|SWE|Average", colnames(raw), value = TRUE,
                      ignore.case = TRUE)[1]
    if (is.na(col_start) || is.na(col_swe)) next
    
    df <- data.frame(
      date_start = tryCatch(
        lubridate::parse_date_time(raw[[col_start]],
                                   orders = c("d/m/Y H:M","m/d/Y H:M","Y-m-d H:M","Y-m-d"),
                                   tz = "UTC"),
        error = function(e) NA),
      swe_mm = suppressWarnings(as.numeric(raw[[col_swe]])),
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(date_start), !is.na(swe_mm), swe_mm > 0) %>%
      mutate(year = year(date_start), month = month(date_start))
    
    if (nrow(df) > 0) {
      monthly <- df %>%
        group_by(year, month) %>%
        summarise(swe_mm_mean = mean(swe_mm, na.rm = TRUE),
                  n_readings = n(), .groups = "drop") %>%
        filter(n_readings >= 5) %>%
        mutate(date = as.Date(paste0(year, "-", sprintf("%02d", month), "-01")))
      if (nrow(monthly) > 0) all_loaded[[basename(f)]] <- monthly
    }
  }
  
  if (length(all_loaded) == 0) return(NULL)
  combined   <- bind_rows(all_loaded)
  basin_mean <- combined %>%
    group_by(year, month, date) %>%
    summarise(swe_pillow_mm    = mean(swe_mm_mean, na.rm = TRUE),
              n_stations_avail = n(), .groups = "drop")
  cat(sprintf("  Basin-mean pillow SWE: %d months\n", nrow(basin_mean)))
  list(by_station = combined, basin_mean = basin_mean)
}

bias_correct_era5_swe <- function(era5_data, pillow_data) {
  if (is.null(pillow_data)) {
    cat("  No pillow data — ERA5 SWE used without bias correction\n")
    era5_data$swe_corrected_mm  <- era5_data$swe_mm
    era5_data$swe_corrected_Mm3 <- era5_data$swe_Mm3
    era5_data$swe_bias_factor   <- 1.0
    return(era5_data)
  }
  cat("\n[SWE Bias Correction] Estimating ERA5 SWE bias from snow pillows...\n")
  
  paired <- era5_data %>%
    dplyr::select(year, month, swe_mm, swe_Mm3) %>%
    inner_join(pillow_data$basin_mean %>% dplyr::select(year, month, swe_pillow_mm),
               by = c("year", "month")) %>%
    filter(swe_mm > 1.0, swe_pillow_mm > 1.0)
  
  if (nrow(paired) < MIN_YEARS_FOR_DIST) {
    cat("  Insufficient overlap — using raw ERA5 SWE\n")
    era5_data$swe_corrected_mm  <- era5_data$swe_mm
    era5_data$swe_corrected_Mm3 <- era5_data$swe_Mm3
    era5_data$swe_bias_factor   <- 1.0
    return(era5_data)
  }
  
  bf <- paired %>%
    group_by(month) %>%
    summarise(bias_factor = median(swe_pillow_mm / swe_mm, na.rm = TRUE),
              n_pairs = n(), .groups = "drop") %>%
    mutate(bias_factor_used = if_else(n_pairs >= MIN_YEARS_FOR_DIST,
                                      bias_factor, 1.0))
  
  era5_corrected <- era5_data %>%
    left_join(bf[, c("month","bias_factor_used")], by = "month") %>%
    mutate(
      bias_factor_used  = replace_na(bias_factor_used, 1.0),
      swe_corrected_mm  = swe_mm  * bias_factor_used,
      swe_corrected_Mm3 = swe_Mm3 * bias_factor_used,
      swe_bias_factor   = bias_factor_used
    )
  
  write.csv(bf, file.path(MAIN_OUTPUT_DIR, "era5_extracted", "swe_bias_factors.csv"),
            row.names = FALSE)
  return(era5_corrected)
}

# ============================================================================
# SECTION 7: WSC DISCHARGE DATA
# ============================================================================
load_discharge_data <- function(station_id, data_path = INPUT_WSC_DIR) {
  patterns <- c(
    file.path(data_path, paste0(station_id, "_WSC_DISCHARGE.csv")),
    file.path(data_path, paste0(station_id, "_discharge.csv")),
    file.path(data_path, paste0(station_id, ".csv"))
  )
  raw <- NULL
  for (p in patterns) { if (file.exists(p)) { raw <- read.csv(p, stringsAsFactors=FALSE); break } }
  if (is.null(raw)) { warning(sprintf("No discharge file for %s", station_id)); return(NULL) }
  
  colnames(raw) <- tolower(colnames(raw))
  date_col  <- grep("date",              colnames(raw), value=TRUE)[1]
  disch_col <- grep("discharge|flow|q_|value", colnames(raw), value=TRUE)[1]
  if (is.na(date_col) || is.na(disch_col)) {
    warning(sprintf("Cannot identify date/discharge columns for %s", station_id))
    return(NULL)
  }
  raw$date      <- as.Date(raw[[date_col]])
  raw$discharge <- suppressWarnings(as.numeric(raw[[disch_col]]))
  raw           <- raw[!is.na(raw$date) & !is.na(raw$discharge), ]
  raw$year      <- year(raw$date)
  raw$month     <- month(raw$date)
  raw$station_id <- station_id
  raw$flow_type  <- "observed"
  cat(sprintf("    Observed: %d records (%d-%d)\n",
              nrow(raw), min(raw$year), max(raw$year)))
  return(raw)
}

load_naturalized_flow <- function(station_id, nat_path = INPUT_NATURAL_DIR) {
  fname <- file.path(nat_path, paste0("dat_natural_", station_id, ".csv"))
  if (!file.exists(fname)) {
    warning(sprintf("Naturalized flow file not found: %s", fname))
    return(NULL)
  }
  raw <- read.csv(fname, stringsAsFactors=FALSE)
  colnames(raw) <- tolower(trimws(colnames(raw)))
  data <- raw %>%
    mutate(date      = as.Date(date),
           discharge = suppressWarnings(as.numeric(value))) %>%
    filter(!is.na(date), tolower(symbol) == "s", !is.na(discharge), discharge >= 0) %>%
    mutate(year       = year(date),
           month      = month(date),
           station_id = station_id,
           flow_type  = "naturalized") %>%
    dplyr::select(date, discharge, year, month, station_id, flow_type)
  cat(sprintf("    Naturalized: %d records (%d-%d)\n",
              nrow(data), min(data$year), max(data$year)))
  return(data)
}

fill_missing_data <- function(data, max_fill = MAX_FILL_DAYS) {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  full_d  <- seq(min(data$date), max(data$date), by = "day")
  cdata   <- merge(data.frame(date = full_d), data, by = "date", all.x = TRUE)
  cdata   <- cdata[order(cdata$date), ]
  q       <- cdata$discharge; n <- length(q)
  runs    <- rle(is.na(q))
  ends    <- cumsum(runs$lengths)
  starts  <- ends - runs$lengths + 1
  filled  <- 0L
  
  for (k in seq_along(runs$lengths)) {
    if (!runs$values[k]) next
    rs <- starts[k]; re <- ends[k]; rl <- runs$lengths[k]
    if (rl > max_fill) next
    li <- rs - 1L; ri <- re + 1L
    if (li < 1L || ri > n || is.na(q[li]) || is.na(q[ri])) next
    ctx <- max(1L, li - 15L):min(n, ri + 15L)
    obs <- !is.na(q[ctx])
    if (sum(obs) < 2) next
    sf       <- splinefun(ctx[obs], q[ctx[obs]], method = "monoH.FC")
    q[rs:re] <- pmax(0, sf(rs:re))
    filled   <- filled + rl
  }
  
  cdata$discharge <- q
  cdata$year      <- year(cdata$date)
  cdata$month     <- month(cdata$date)
  if ("station_id" %in% colnames(data)) cdata$station_id <- data$station_id[1]
  if ("flow_type"  %in% colnames(data)) cdata$flow_type  <- data$flow_type[1]
  cat(sprintf("      Filled %d values (PCHIP)\n", filled))
  return(cdata)
}

aggregate_to_monthly_volume <- function(daily_data) {
  if (is.null(daily_data) || nrow(daily_data) == 0) return(NULL)
  daily_data %>%
    group_by(year, month) %>%
    summarise(n_valid        = sum(!is.na(discharge)),
              discharge_mean = mean(discharge, na.rm = TRUE),
              .groups        = "drop") %>%
    filter(n_valid >= MIN_DAYS_PER_MONTH) %>%
    mutate(
      days_in_mo = days_in_month(as.Date(paste0(year, "-", sprintf("%02d", month), "-01"))),
      volume_Mm3 = discharge_mean * days_in_mo * 86400 / 1e6,
      date       = as.Date(paste0(year, "-", sprintf("%02d", month), "-15")),
      water_year = if_else(month >= WATER_YEAR_START_MONTH, year, year - 1L)
    )
}

# ============================================================================
# SECTION 8: LEGACY FORECAST MODELS (optional; not used in dual-track pipeline)
# ============================================================================
# Retained for reference. The H-SWSI and O-SWSI tracks use OBSERVED naturalized
# streamflow (Qnat) and OBSERVED reservoir storage, so regression forecasts are
# not required for index computation. Uncomment and call build_forecast_models()
# separately if seasonal outlook products are also needed.
#
# build_forecast_models <- function(era5_data, monthly_vol_primary) {
#   ...  [original Section 8 code unchanged]
# }

# ============================================================================
# SECTION 9: GAREN (1993) REVISED — Z-SCORE MATHEMATICAL CORE
# ============================================================================
#
# The original Garen (1993) formula was:  SWSI = (P − 50) / 12
# where P is the Weibull rank percentile of the composite variable.
#
# The REVISED formulation (now standard in the literature) uses:
#
#   Step 1:  Pi  = non-exceedance probability [%] from the historical CDF
#   Step 2:  Zi  = Φ⁻¹(Pi / 100)   — standard normal variate
#   Step 3:  SWSI = Σ(Wi · Zi)      — weighted sum across components
#   Step 4:  Clip to ±SWSI_CAP (±4.2)
#
# Benefits over original:
#   • Each component Zi is individually standardised on a common normal scale.
#   • Weights (Wi) allow component importance to differ by track.
#   • Baseline period can differ per component (pre- vs post-regulation).
#   • Missing components handled by renormalising remaining weights.

# ── 9a. Weibull plotting-position estimator ──────────────────────────────────
compute_weibull_prob <- function(hist_vals, x_new) {
  v <- sort(hist_vals[!is.na(hist_vals)])
  n <- length(v)
  if (n == 0 || is.na(x_new)) return(NA_real_)
  if (length(unique(v)) < 2)  return(50)           # neutral value when no variance
  P <- seq_len(n) / (n + 1) * 100
  pmax(0.01, pmin(99.99, approx(v, P, xout = x_new, rule = 2, method = "linear")$y))
}

# ── 9b. L-moment parametric CDF fitting ─────────────────────────────────────
fit_parametric_cdf <- function(hist_vals) {
  v <- hist_vals[!is.na(hist_vals) & hist_vals > 0]
  if (length(v) < MIN_YEARS_FOR_DIST || var(v) < 1e-6) return(NULL)
  
  lmom <- tryCatch(lmomco::lmoms(v, nmom = 4), error = function(e) NULL)
  if (is.null(lmom)) return(NULL)
  
  fits <- lapply(c("gamma", "lnorm", "norm"), function(d) {
    tryCatch({
      if (d == "gamma") {
        par <- lmomco::pargam(lmom);   cdf <- function(x) lmomco::cdfgam(x,   par)
      } else if (d == "lnorm") {
        par <- lmomco::parlnorm(lmom); cdf <- function(x) lmomco::cdflnorm(x, par)
      } else {
        par <- lmomco::parnor(lmom);   cdf <- function(x) lmomco::cdfnor(x,   par)
      }
      list(dist = d, cdf = cdf)
    }, error = function(e) list(dist = d, cdf = NULL))
  })
  
  best <- fits[[which.min(sapply(fits, function(x) if (is.null(x$cdf)) Inf else 0))]]
  if (is.null(best$cdf)) return(NULL)
  list(best_dist = best$dist, cdf_fn = best$cdf)
}

# ── 9c. Probability → Z-score  (REVISED CORE) ───────────────────────────────
prob_to_zscore <- function(P_pct, cap = SWSI_CAP) {
  # P_pct : non-exceedance probability in percent [0, 100]
  # Returns Zi = Φ⁻¹(P_pct/100), capped to ±cap to match the SWSI scale.
  P_pct <- pmax(0.01, pmin(99.99, P_pct))   # guard ±Inf from qnorm
  pmax(-cap, pmin(cap, qnorm(P_pct / 100)))
}

# ── 9d. Per-component Z-score computation ────────────────────────────────────
# KEY: baseline_rows
#   NULL             → use all rows in df for that month as the historical CDF
#   integer vector   → restrict calibration CDF to those row indices only
#                      (enables separate pre- / post-regulation baselines)
#
# Z-scores are computed for ALL rows in df, evaluated against the baseline CDF.
#
compute_component_zscore <- function(df, value_col, label, baseline_rows = NULL) {
  df[[paste0("P_",    label)]] <- NA_real_
  df[[paste0("Z_",    label)]] <- NA_real_
  df[[paste0("dist_", label)]] <- NA_character_
  
  for (m in 1:12) {
    row_m     <- which(df$month == m)
    hist_rows <- if (!is.null(baseline_rows)) intersect(row_m, baseline_rows) else row_m
    hist_vals <- df[[value_col]][hist_rows]
    
    n_hist <- sum(!is.na(hist_vals))
    if (n_hist < MIN_YEARS_FOR_DIST) {
      message(sprintf("  [Z: %s] Month %d — only %d baseline values; skipping",
                      label, m, n_hist))
      next
    }
    
    param_fit <- fit_parametric_cdf(hist_vals)
    
    for (r in row_m) {
      v_r <- df[[value_col]][r]
      if (is.na(v_r)) next
      
      P_emp <- compute_weibull_prob(hist_vals, v_r)
      
      if (!is.null(param_fit)) {
        P_use <- tryCatch(
          pmax(0.01, pmin(99.99, param_fit$cdf_fn(v_r) * 100)),
          error = function(e) P_emp)
        dist_used <- param_fit$best_dist
      } else {
        P_use     <- P_emp
        dist_used <- "empirical"
      }
      
      df[[paste0("P_",    label)]][r] <- P_use
      df[[paste0("Z_",    label)]][r] <- prob_to_zscore(P_use)
      df[[paste0("dist_", label)]][r] <- dist_used
    }
  }
  return(df)
}

# ── 9e. Weighted SWSI aggregation ────────────────────────────────────────────
# When one or more component Z-scores are NA for a given row, the remaining
# available weights are renormalised so the index still sums to 1.
# A minimum of min_components non-NA components is required; otherwise NA.
#
compute_weighted_swsi <- function(df, z_cols, weights, label = "swsi",
                                  min_components = 2L) {
  w <- weights / sum(weights, na.rm = TRUE)
  
  result <- vapply(seq_len(nrow(df)), function(i) {
    avail <- vapply(z_cols, function(col) {
      col %in% colnames(df) && !is.na(df[[col]][i])
    }, logical(1))
    
    if (sum(avail) < min_components) return(NA_real_)
    
    w_i <- w[avail] / sum(w[avail])                          # renormalise
    z_i <- vapply(z_cols[avail], function(col) df[[col]][i], numeric(1))
    pmax(-SWSI_CAP, pmin(SWSI_CAP, sum(z_i * w_i)))
  }, numeric(1))
  
  df[[label]] <- result
  return(df)
}

# ── 9f. Shared drought classification ────────────────────────────────────────
classify_swsi <- function(x) {
  cut(x,
      breaks = c(-Inf, -4.0, -3.0, -2.0, -1.0, 0.0, 2.0, Inf),
      labels = c("Extreme drought", "Severe drought", "Moderate drought",
                 "Drought watch", "Below normal", "Near normal", "Abundant supply"),
      right  = TRUE)
}

# ============================================================================
# SECTION 10: H-SWSI — HYDROLOGICAL DROUGHT INDEX
# ============================================================================
#
#  Components  : Qnat · LakeNat · Precip · SWE
#  Baseline    : Pre-regulation naturalized record
#    – Qnat    : pre-REGULATION_START rows (if ≥ MIN_YEARS_FOR_DIST), else full record
#    – LakeNat : full GloLakes record (satellite era starts 1984; no pre-1952 data)
#    – Precip  : full ERA5 record (ERA5 starts 1940; no pre-1952 monthly data)
#    – SWE     : full ERA5 record (same rationale)
#  Weights     : HSWSI_WEIGHTS (default 0.25 each)
#  Excluded    : Reservoir storage, regulated downstream flow

compute_hswsi <- function(era5_data,
                          monthly_vol_qnat,
                          nat_lake_agg    = NULL,
                          weights         = HSWSI_WEIGHTS,
                          baseline_end    = HSWSI_BASELINE_END) {
  
  sep <- strrep("-", 62)
  cat(sprintf("\n%s\n[TRACK 1]  H-SWSI — Hydrological Drought Index\n%s\n", sep, sep))
  cat("  Objective  : Isolate climate-driven anomalies (unregulated basin)\n")
  cat(sprintf("  Formula    : SWSI = Σ(Wi·Zi),  Zi = Φ⁻¹(Pi),  cap = ±%.1f\n", SWSI_CAP))
  cat(sprintf("  Baseline   : Naturalized flow pre-%s; lakes/climate full record\n",
              baseline_end))
  cat(sprintf("  Weights    : Qnat=%.2f  LakeNat=%.2f  Precip=%.2f  SWE=%.2f\n",
              weights["Qnat"], weights["LakeNat"], weights["Precip"], weights["SWE"]))
  cat("  Excluded   : Reservoir storage, regulated downstream flow\n")
  
  # ── Assemble dataset ──────────────────────────────────────────────────────
  base <- monthly_vol_qnat %>%
    dplyr::select(year, month, water_year, date, volume_Mm3) %>%
    rename(Qnat_Mm3 = volume_Mm3) %>%
    left_join(
      era5_data %>% dplyr::select(year, month, swe_corrected_Mm3, precip_ytd_Mm3),
      by = c("year", "month")
    )
  
  # ── Natural lake storage (Stuart / Fraser) ───────────────────────
  has_lake <- !is.null(nat_lake_agg) && nrow(nat_lake_agg) > 0
  if (has_lake) {
    base <- base %>%
      left_join(nat_lake_agg %>% dplyr::select(year, month, nat_lake_total_Mm3),
                by = c("year", "month"))
    cat(sprintf("\n  LakeNat overlap: %d months (%s – %s)\n",
                sum(!is.na(base$nat_lake_total_Mm3)),
                min(nat_lake_agg$date), max(nat_lake_agg$date)))
  } else {
    base$nat_lake_total_Mm3 <- NA_real_
    cat("\n  LakeNat : ABSENT — weight redistributed to remaining components\n")
  }
  
  # ── Qnat baseline: pre-regulation rows if sufficient ─────────────────────
  qnat_pre_rows <- which(!is.na(base$Qnat_Mm3) & base$date < baseline_end)
  if (length(qnat_pre_rows) >= MIN_YEARS_FOR_DIST) {
    qnat_baseline_rows <- qnat_pre_rows
    cat(sprintf("  Qnat baseline : pre-%d  (%d months)\n",
                year(baseline_end), length(qnat_pre_rows)))
  } else {
    qnat_baseline_rows <- NULL
    cat(sprintf("  Qnat baseline : full record (only %d pre-regulation months)\n",
                length(qnat_pre_rows)))
  }
  cat("  LakeNat / Precip / SWE baseline : full available record\n\n")
  
  # ── Component Z-scores ────────────────────────────────────────────────────
  cat("  [1/4] Qnat   — naturalized streamflow\n")
  base <- compute_component_zscore(base, "Qnat_Mm3",          "Qnat",    qnat_baseline_rows)
  
  if (has_lake && !all(is.na(base$nat_lake_total_Mm3))) {
    cat("  [2/4] LakeNat — aggregate natural lake storage\n")
    base <- compute_component_zscore(base, "nat_lake_total_Mm3", "LakeNat", NULL)
  } else {
    base <- base %>%
      mutate(P_LakeNat = NA_real_, Z_LakeNat = NA_real_, dist_LakeNat = NA_character_)
  }
  
  cat("  [3/4] Precip  — basin-averaged YTD precipitation\n")
  base <- compute_component_zscore(base, "precip_ytd_Mm3",    "Precip",  NULL)
  
  cat("  [4/4] SWE     — snow water equivalent (bias-corrected)\n")
  base <- compute_component_zscore(base, "swe_corrected_Mm3", "SWE",     NULL)
  
  # ── Weighted SWSI ──────────────────────────────────────────────────────────
  z_cols <- c("Z_Qnat", "Z_LakeNat", "Z_Precip", "Z_SWE")
  w_vec  <- unname(weights[c("Qnat", "LakeNat", "Precip", "SWE")])
  base   <- compute_weighted_swsi(base, z_cols, w_vec, label = "hswsi")
  
  base$drought_category <- classify_swsi(base$hswsi)
  base$track            <- "H-SWSI"
  
  cat(sprintf("\n  ✓ H-SWSI  : %d valid months\n", sum(!is.na(base$hswsi))))
  return(base)
}

# ============================================================================
# SECTION 11: O-SWSI — OPERATIONAL SUPPLY DROUGHT INDEX
# ============================================================================
#
#  Components  : Active Storage (R) · SWE · Unregulated Inflow (Qin)
#  Baseline    : Post-regulation data (OSWSI_BASELINE_START = 1984-01-01)
#    – All three components share the same 1984+ baseline so that
#      standardisation reflects the post-Kemano-diversion operating regime.
#  Weights     : OSWSI_WEIGHTS (Storage=0.50, SWE=0.35, Qin=0.15)
#  Excluded    : Raw basin precipitation, regulated downstream flow
#
#  Reservoir adjustments (isolate usable active storage):
#    Active_Storage = Total_GloLakes − dead_storage − env_flow_reserve
#    Qin (optional) = Naturalized flow − baseline_Kemano_diversion

# ── 11a. Reservoir adjustment ─────────────────────────────────────────────────
adjust_reservoir_storage <- function(reservoir_data,
                                     dead_storage     = KENNEY_DEAD_STORAGE_MM3,
                                     env_flow_reserve = ENV_FLOW_RESERVE_MM3) {
  cat(sprintf("\n[O-SWSI] Reservoir storage adjustment:\n"))
  cat(sprintf("  Total storage\n"))
  cat(sprintf("    − dead storage       %.0f Mm³  (minimum operating level)\n", dead_storage))
  cat(sprintf("    − env-flow reserve   %.0f Mm³  (legally reserved volume)\n", env_flow_reserve))
  cat(sprintf("  = active_storage_Mm3  (clipped to ≥ 0)\n"))
  
  reservoir_data %>%
    mutate(
      dead_storage_Mm3   = dead_storage,
      env_reserve_Mm3    = env_flow_reserve,
      active_storage_Mm3 = pmax(0, storage - dead_storage - env_flow_reserve)
    )
}

# ── 11b. O-SWSI computation ───────────────────────────────────────────────────
compute_oswsi <- function(era5_data,
                          reservoir_data_adj,
                          monthly_vol_qin       = NULL,
                          weights               = OSWSI_WEIGHTS,
                          baseline_start        = OSWSI_BASELINE_START,
                          kemano_div_mm3_month  = KEMANO_BASELINE_DIV_MM3_MONTH) {
  
  sep <- strrep("-", 62)
  cat(sprintf("\n%s\n[TRACK 2]  O-SWSI — Operational Supply Drought Index\n%s\n", sep, sep))
  cat("  Objective  : Measure actual system stress & water management vulnerability\n")
  cat(sprintf("  Formula    : SWSI = Σ(Wi·Zi),  Zi = Φ⁻¹(Pi),  cap = ±%.1f\n", SWSI_CAP))
  cat(sprintf("  Baseline   : Post-regulation data from %s\n", baseline_start))
  cat(sprintf("  Weights    : Storage=%.2f  SWE=%.2f  Qin=%.2f\n",
              weights["Storage"], weights["SWE"], weights["Qin"]))
  cat("  Excluded   : Raw basin precipitation, regulated downstream flow\n")
  
  # ── Assemble dataset ──────────────────────────────────────────────────────
  base <- reservoir_data_adj %>%
    dplyr::select(date, year, month, water_year,
                  storage, active_storage_Mm3) %>%
    left_join(
      era5_data %>% dplyr::select(year, month, swe_corrected_Mm3),
      by = c("year", "month")
    )
  
  # ── Unregulated inflow (Qin) ─────────────────────────────────────────────
  # 08JC002 naturalized flow is used as the proxy for the pre-diversion inflow
  # that would have entered the reservoir under natural conditions.
  # Optionally subtract the long-term baseline Kemano diversion (KEMANO_BASELINE_DIV_MM3_MONTH)
  # to isolate purely local supply; default = 0 (no adjustment).
  if (!is.null(monthly_vol_qin)) {
    base <- base %>%
      left_join(
        monthly_vol_qin %>% dplyr::select(year, month, volume_Mm3) %>%
          rename(Qin_raw_Mm3 = volume_Mm3),
        by = c("year", "month")
      ) %>%
      mutate(Qin_Mm3 = pmax(0, Qin_raw_Mm3 - kemano_div_mm3_month))
    cat(sprintf("\n  Qin source  : %s naturalized  |  Kemano deduction: %.0f Mm³/month\n",
                OSWSI_QIN_STATION, kemano_div_mm3_month))
    cat(sprintf("  Qin overlap : %d months\n", sum(!is.na(base$Qin_Mm3))))
  } else {
    base$Qin_raw_Mm3 <- NA_real_
    base$Qin_Mm3     <- NA_real_
    cat("\n  Qin : ABSENT — Storage + SWE weights redistributed\n")
  }
  
  # ── Post-regulation baseline rows ────────────────────────────────────────
  baseline_rows <- which(base$date >= baseline_start)
  cat(sprintf("  Baseline    : %d rows (%s to %s)\n\n",
              length(baseline_rows), baseline_start, max(base$date)))
  
  # ── Component Z-scores ────────────────────────────────────────────────────
  cat("  [1/3] Active Reservoir Storage\n")
  base <- compute_component_zscore(base, "active_storage_Mm3", "Storage", baseline_rows)
  
  cat("  [2/3] SWE — snow water equivalent (bias-corrected)\n")
  base <- compute_component_zscore(base, "swe_corrected_Mm3",  "SWE",     baseline_rows)
  
  if (!all(is.na(base$Qin_Mm3))) {
    cat("  [3/3] Qin — unregulated reservoir inflow\n")
    base <- compute_component_zscore(base, "Qin_Mm3", "Qin", baseline_rows)
  } else {
    base <- base %>%
      mutate(P_Qin = NA_real_, Z_Qin = NA_real_, dist_Qin = NA_character_)
  }
  
  # ── Weighted SWSI ──────────────────────────────────────────────────────────
  z_cols <- c("Z_Storage", "Z_SWE", "Z_Qin")
  w_vec  <- unname(weights[c("Storage", "SWE", "Qin")])
  base   <- compute_weighted_swsi(base, z_cols, w_vec, label = "oswsi")
  
  base$drought_category <- classify_swsi(base$oswsi)
  base$track            <- "O-SWSI"
  
  cat(sprintf("\n  ✓ O-SWSI  : %d valid months\n", sum(!is.na(base$oswsi))))
  return(base)
}

# ============================================================================
# SECTION 12: DROUGHT EVENT IDENTIFICATION (DUAL-TRACK)
# ============================================================================
# Works for both tracks: pass index_col = "hswsi" or "oswsi".
# Identifies runs below TRIGGER_WATCH and summarises their severity,
# duration, and primary driver (component with lowest Z-score minimum).

identify_droughts <- function(swsi_df, index_col = "hswsi",
                              track_label = "H-SWSI") {
  if (is.null(swsi_df) || !index_col %in% colnames(swsi_df)) return(NULL)
  
  safe_min <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_real_ else min(x) }
  
  data <- swsi_df %>%
    arrange(date) %>%
    mutate(swsi_val   = .data[[index_col]],
           in_drought = !is.na(swsi_val) & swsi_val <= TRIGGER_WATCH)
  data$event_id <- with(rle(data$in_drought), rep(seq_along(lengths), lengths))
  
  # Identify Z-score driver columns present in this track
  z_driver_cols <- intersect(
    c("Z_Qnat", "Z_LakeNat", "Z_Precip", "Z_SWE",   # H-SWSI
      "Z_Storage", "Z_Qin"),                           # O-SWSI
    colnames(data)
  )
  
  event_info <- data %>%
    filter(in_drought) %>%
    group_by(event_id) %>%
    summarise(
      start_date       = min(date),
      end_date         = max(date),
      duration_months  = n(),
      start_year       = first(year),
      water_year_start = first(water_year),
      swsi_min         = safe_min(swsi_val),
      swsi_mean        = mean(swsi_val, na.rm = TRUE),
      across(all_of(z_driver_cols),
             list(min = ~ safe_min(.x)),
             .names = "min_{.col}"),
      .groups = "drop"
    ) %>%
    filter(duration_months >= 1) %>%
    mutate(
      track         = track_label,
      peak_category = dplyr::case_when(
        swsi_min <= -4.0 ~ "Extreme drought",
        swsi_min <= -3.0 ~ "Severe drought",
        swsi_min <= -2.0 ~ "Moderate drought",
        TRUE             ~ "Drought watch"
      ),
      primary_driver = {
        min_z_df <- dplyr::select(., starts_with("min_Z_"))
        sapply(seq_len(nrow(min_z_df)), function(i) {
          row  <- as.numeric(min_z_df[i, ])
          if (all(is.na(row))) return("Unknown")
          col  <- sub("^min_Z_", "", colnames(min_z_df)[which.min(row)])
          switch(col,
                 Qnat    = "Naturalized streamflow deficit",
                 LakeNat = "Natural lake storage deficit",
                 Precip  = "Precipitation deficit",
                 SWE     = "Snowpack deficit",
                 Storage = "Reservoir storage deficit",
                 Qin     = "Inflow deficit",
                 col)
        })
      }
    )
  
  n_emergency <- sum(event_info$swsi_min <= TRIGGER_EMERGENCY, na.rm = TRUE)
  cat(sprintf("  [%s] Drought events: %d  |  Emergency (≤%.1f): %d\n",
              track_label, nrow(event_info), TRIGGER_EMERGENCY, n_emergency))
  
  list(monthly_data   = data,
       drought_events  = event_info,
       n_events        = nrow(event_info),
       n_emergency     = n_emergency)
}

# ============================================================================
# SECTION 13: MAIN PROCESSING PIPELINE
# ============================================================================
cat("\n", strrep("=", 70), "\n",
    "DUAL-TRACK SWSI  —  Nechako Basin (Garen 1993 Revised)\n",
    "Track 1: H-SWSI (Hydrological)   |   Track 2: O-SWSI (Operational)\n",
    strrep("=", 70), "\n", sep = "")

# ── A. Basin boundary ─────────────────────────────────────────────────────────
basin_vect <- load_basin_boundary(BASIN_KMZ_PATH)
terra::crs(basin_vect, describe = TRUE)

# ── B. Reservoir storage (O-SWSI primary component) ──────────────────────────
reservoir_data_raw <- load_reservoir_storage()

# ── C. Natural lake storage (H-SWSI LakeNat component) ───────────────────────
nat_lake_df  <- load_natural_lake_storage()
nat_lake_agg <- aggregate_natural_lakes(nat_lake_df)

if (!is.null(nat_lake_agg)) {
  cat(sprintf("  Lake aggregate: %d months (%s – %s) | lakes: %s\n",
              nrow(nat_lake_agg), min(nat_lake_agg$date), max(nat_lake_agg$date),
              nat_lake_agg$lakes_included[1]))
}

# ── D. ERA5-Land climate (SWE + precipitation, masked to Nechako basin) ───────
era5_raw    <- load_era5_all(INPUT_ERA5_DIR, basin_vect = basin_vect)
pillow_data <- load_all_snow_pillows(INPUT_SNOWPILLOW_DIR)
era5_data   <- bias_correct_era5_swe(era5_raw, pillow_data)
write.csv(era5_data,
          file.path(MAIN_OUTPUT_DIR, "era5_extracted", "era5_basin_monthly_corrected.csv"),
          row.names = FALSE)

# ── E. WSC streamflow — load all stations assigned to either track ─────────────
cat("\n[WSC] Loading discharge / naturalized flow data...\n")
all_monthly_vol <- list()

for (i in seq_len(nrow(STATIONS))) {
  sid  <- STATIONS$StationID[i]
  snam <- STATIONS$Name[i]
  
  if (!STATIONS$Use_HSWSI_Qnat[i] && !STATIONS$Use_OSWSI_Qin[i]) next
  
  cat(sprintf("\n  %s  —  %s\n", sid, snam))
  if (!is.na(STATIONS$Notes[i])) cat(sprintf("    Note: %s\n", STATIONS$Notes[i]))
  
  # Prefer naturalized flow for all stations used in either track
  daily_data <- if (STATIONS$Has_NatFlow[i])
    load_naturalized_flow(sid, INPUT_NATURAL_DIR)
  else
    load_discharge_data(sid, INPUT_WSC_DIR)
  
  if (is.null(daily_data)) {
    cat(sprintf("    ✗ No data for %s\n", sid))
    next
  }
  
  data_filled <- fill_missing_data(daily_data)
  monthly_vol <- aggregate_to_monthly_volume(data_filled)
  if (!is.null(monthly_vol)) all_monthly_vol[[sid]] <- monthly_vol
}

# Verify primary stations
if (!HSWSI_QNAT_PRIMARY %in% names(all_monthly_vol))
  stop(sprintf("Primary Qnat station %s not loaded — cannot compute H-SWSI",
               HSWSI_QNAT_PRIMARY))

if (!OSWSI_QIN_STATION %in% names(all_monthly_vol))
  warning(sprintf("Qin station %s not loaded — O-SWSI Qin component will be absent",
                  OSWSI_QIN_STATION))

# ── F. Reservoir storage adjustment (O-SWSI) ─────────────────────────────────
reservoir_data_adj <- adjust_reservoir_storage(reservoir_data_raw)
write.csv(reservoir_data_adj,
          file.path(MAIN_OUTPUT_DIR, "processed_data", "reservoir_storage_adjusted.csv"),
          row.names = FALSE)

# ── G. Compute H-SWSI (Track 1) ──────────────────────────────────────────────
hswsi_result <- compute_hswsi(
  era5_data        = era5_data,
  monthly_vol_qnat = all_monthly_vol[[HSWSI_QNAT_PRIMARY]],
  nat_lake_agg     = nat_lake_agg,
  weights          = HSWSI_WEIGHTS,
  baseline_end     = HSWSI_BASELINE_END
)

hswsi_droughts <- identify_droughts(hswsi_result, "hswsi", "H-SWSI")

write.csv(hswsi_result,
          file.path(MAIN_OUTPUT_DIR, "hswsi", "hswsi_timeseries.csv"),
          row.names = FALSE)
if (!is.null(hswsi_droughts$drought_events) &&
    nrow(hswsi_droughts$drought_events) > 0)
  write.csv(hswsi_droughts$drought_events,
            file.path(MAIN_OUTPUT_DIR, "hswsi", "droughts", "hswsi_drought_events.csv"),
            row.names = FALSE)

# ── H. Compute O-SWSI (Track 2) ──────────────────────────────────────────────
qin_monthly <- if (OSWSI_QIN_STATION %in% names(all_monthly_vol))
  all_monthly_vol[[OSWSI_QIN_STATION]] else NULL

oswsi_result <- compute_oswsi(
  era5_data          = era5_data,
  reservoir_data_adj = reservoir_data_adj,
  monthly_vol_qin    = qin_monthly,
  weights            = OSWSI_WEIGHTS,
  baseline_start     = OSWSI_BASELINE_START
)

oswsi_droughts <- identify_droughts(oswsi_result, "oswsi", "O-SWSI")

write.csv(oswsi_result,
          file.path(MAIN_OUTPUT_DIR, "oswsi", "oswsi_timeseries.csv"),
          row.names = FALSE)
if (!is.null(oswsi_droughts$drought_events) &&
    nrow(oswsi_droughts$drought_events) > 0)
  write.csv(oswsi_droughts$drought_events,
            file.path(MAIN_OUTPUT_DIR, "oswsi", "droughts", "oswsi_drought_events.csv"),
            row.names = FALSE)

# ── I. Wide-format monthly matrices ──────────────────────────────────────────
make_wide <- function(df, index_col) {
  df %>%
    mutate(Year  = year(date),
           Month = month(date, label = TRUE, abbr = TRUE)) %>%
    dplyr::select(Year, Month, value = all_of(index_col)) %>%
    tidyr::pivot_wider(names_from = Month, values_from = value) %>%
    arrange(Year)
}
write.csv(make_wide(hswsi_result, "hswsi"),
          file.path(MAIN_OUTPUT_DIR, "hswsi", "hswsi_wide_by_month.csv"), row.names = FALSE)
write.csv(make_wide(oswsi_result, "oswsi"),
          file.path(MAIN_OUTPUT_DIR, "oswsi", "oswsi_wide_by_month.csv"), row.names = FALSE)
cat("✓ Wide-format monthly CSV files saved\n")

# ============================================================================
# SECTION 14: VISUALISATIONS
# ============================================================================

# ── Shared ribbon-plot builder ────────────────────────────────────────────────
swsi_ribbon_plot <- function(df, index_col, title_str, subtitle_str,
                             focus_start = NULL, focus_end = NULL,
                             date_break  = "5 years") {
  df$val <- df[[index_col]]
  df     <- df[!is.na(df$val), ]
  
  p <- ggplot(df, aes(x = date, y = val)) +
    geom_ribbon(aes(ymin = pmin(val, 0), ymax = 0),
                fill = "#d73027", alpha = 0.35) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(val, 0)),
                fill = "#4575b4", alpha = 0.35) +
    geom_line(linewidth = 0.7, colour = "grey20") +
    geom_hline(yintercept = c(0, -1, -2, -3, -4),
               linetype   = c("solid","dashed","dashed","dashed","dashed"),
               colour     = c("grey50","#fc8d59","#d73027","#a50026","#4c0000"),
               linewidth  = 0.4) +
    annotate("text",
             x      = max(df$date) - 180,
             y      = c(-0.85, -1.85, -2.85, -3.85),
             label  = c("Watch", "Emergency", "Severe", "Extreme"),
             hjust  = 1, size = 2.6,
             colour = c("#fc8d59","#d73027","#a50026","#4c0000")) +
    scale_x_date(date_breaks = date_break, date_labels = "%Y") +
    scale_y_continuous(limits = c(-SWSI_CAP - 0.15, SWSI_CAP + 0.15),
                       breaks = seq(-4, 4, 1)) +
    labs(title    = title_str,
         subtitle = subtitle_str,
         x        = NULL,
         y        = sprintf("SWSI  (Zi = Φ⁻¹(Pi),  cap ±%.1f)", SWSI_CAP)) +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  
  if (!is.null(focus_start) && !is.null(focus_end))
    p <- p + annotate("rect",
                      xmin = focus_start, xmax = focus_end,
                      ymin = -Inf,         ymax = Inf,
                      fill = "grey80",      alpha = 0.40)
  p
}

# ── H-SWSI: full time series ──────────────────────────────────────────────────
p_h_full <- swsi_ribbon_plot(
  hswsi_result, "hswsi",
  title_str    = "H-SWSI  —  Hydrological Drought Index  |  Nechako Basin",
  subtitle_str = sprintf("Components: Qnat·LakeStorage·Precip·SWE  |  Station: %s  |  Grey: 2022–2025",
                         HSWSI_QNAT_PRIMARY),
  focus_start = FOCUS_START, focus_end = FOCUS_END
)
ggsave(file.path(MAIN_OUTPUT_DIR, "hswsi", "hswsi_timeseries_full.png"),
       p_h_full, width = 14, height = 5, dpi = 300)
cat("✓ H-SWSI full time series plot saved\n")

# ── O-SWSI: full time series ──────────────────────────────────────────────────
p_o_full <- swsi_ribbon_plot(
  oswsi_result, "oswsi",
  title_str    = "O-SWSI  —  Operational Supply Drought Index  |  Nechako Basin",
  subtitle_str = sprintf("Components: Active Storage·SWE·Qin  |  Baseline: %s+  |  Grey: 2022–2025",
                         year(OSWSI_BASELINE_START)),
  focus_start = FOCUS_START, focus_end = FOCUS_END
)
ggsave(file.path(MAIN_OUTPUT_DIR, "oswsi", "oswsi_timeseries_full.png"),
       p_o_full, width = 14, height = 5, dpi = 300)
cat("✓ O-SWSI full time series plot saved\n")

# ── H-SWSI: component Z-score panels ─────────────────────────────────────────
h_comp_cols <- c("Z_Qnat", "Z_LakeNat", "Z_Precip", "Z_SWE")
h_comp_labels <- c("Naturalized Flow (Qnat)", "Natural Lake Storage",
                   "YTD Precipitation",       "SWE")
names(h_comp_labels) <- h_comp_cols

h_long <- hswsi_result %>%
  dplyr::select(date, all_of(h_comp_cols), hswsi) %>%
  filter(!is.na(hswsi)) %>%
  tidyr::pivot_longer(c(all_of(h_comp_cols), hswsi),
                      names_to = "component", values_to = "Z") %>%
  mutate(component = recode(component,
                            Z_Qnat    = "Naturalized Flow",
                            Z_LakeNat = "Lake Storage",
                            Z_Precip  = "Precipitation",
                            Z_SWE     = "SWE",
                            hswsi     = "★ H-SWSI (composite)"),
         component = factor(component,
                            levels = c("Naturalized Flow","Lake Storage","Precipitation","SWE","★ H-SWSI (composite)")))

p_h_comp <- ggplot(h_long, aes(x = date, y = Z, colour = component)) +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = 0,  linetype = "dashed", colour = "grey50",  linewidth = 0.35) +
  geom_hline(yintercept = -1, linetype = "dashed", colour = "#fc8d59", linewidth = 0.35) +
  geom_hline(yintercept = -2, linetype = "dashed", colour = "#d73027", linewidth = 0.35) +
  geom_line(linewidth = 0.6) +
  scale_colour_brewer(palette = "Set1", guide = "none") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  labs(title    = "H-SWSI Component Z-scores  —  Nechako Basin",
       subtitle = "Each panel: Zi = Φ⁻¹(Pi)  |  Bottom panel: composite weighted sum",
       x = NULL, y = "Z-score") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(MAIN_OUTPUT_DIR, "hswsi", "components", "hswsi_component_zscores.png"),
       p_h_comp, width = 14, height = 11, dpi = 300)
cat("✓ H-SWSI component Z-score panels saved\n")

# ── O-SWSI: component Z-score panels ─────────────────────────────────────────
o_long <- oswsi_result %>%
  dplyr::select(date, Z_Storage, Z_SWE, Z_Qin, oswsi) %>%
  filter(!is.na(oswsi)) %>%
  tidyr::pivot_longer(c(Z_Storage, Z_SWE, Z_Qin, oswsi),
                      names_to = "component", values_to = "Z") %>%
  mutate(component = recode(component,
                            Z_Storage = "Active Reservoir Storage",
                            Z_SWE     = "SWE",
                            Z_Qin     = "Unregulated Inflow (Qin)",
                            oswsi     = "★ O-SWSI (composite)"),
         component = factor(component,
                            levels = c("Active Reservoir Storage","SWE","Unregulated Inflow (Qin)","★ O-SWSI (composite)")))

p_o_comp <- ggplot(o_long, aes(x = date, y = Z, colour = component)) +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  geom_hline(yintercept = 0,  linetype = "dashed", colour = "grey50",  linewidth = 0.35) +
  geom_hline(yintercept = -1, linetype = "dashed", colour = "#fc8d59", linewidth = 0.35) +
  geom_hline(yintercept = -2, linetype = "dashed", colour = "#d73027", linewidth = 0.35) +
  geom_line(linewidth = 0.6) +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  labs(title    = "O-SWSI Component Z-scores  —  Nechako Basin",
       subtitle = "Each panel: Zi = Φ⁻¹(Pi)  |  Bottom panel: composite weighted sum",
       x = NULL, y = "Z-score") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))
ggsave(file.path(MAIN_OUTPUT_DIR, "oswsi", "components", "oswsi_component_zscores.png"),
       p_o_comp, width = 14, height = 10, dpi = 300)
cat("✓ O-SWSI component Z-score panels saved\n")

# ── Track comparison: H-SWSI vs O-SWSI ───────────────────────────────────────
compare_df <- inner_join(
  hswsi_result %>% dplyr::select(date, hswsi),
  oswsi_result %>% dplyr::select(date, oswsi),
  by = "date"
) %>%
  filter(!is.na(hswsi) | !is.na(oswsi)) %>%
  tidyr::pivot_longer(c(hswsi, oswsi), names_to = "track", values_to = "swsi") %>%
  mutate(track = recode(track,
                        hswsi = "H-SWSI  (Climate-driven)",
                        oswsi = "O-SWSI  (Operational supply)"))

p_compare <- ggplot(compare_df, aes(x = date, y = swsi, colour = track)) +
  annotate("rect",
           xmin = FOCUS_START, xmax = FOCUS_END,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.40) +
  geom_hline(yintercept = c(0, -1, -2, -3, -4),
             linetype = c("solid","dashed","dashed","dashed","dashed"),
             colour   = c("grey50","#fc8d59","#d73027","#a50026","#4c0000"),
             linewidth = 0.4) +
  geom_line(linewidth = 0.80, alpha = 0.90, na.rm = TRUE) +
  scale_colour_manual(values = c("H-SWSI  (Climate-driven)"    = "#2c7bb6",
                                 "O-SWSI  (Operational supply)" = "#d7191c")) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(limits = c(-SWSI_CAP - 0.15, SWSI_CAP + 0.15),
                     breaks = seq(-4, 4, 1)) +
  labs(title    = "H-SWSI vs O-SWSI  —  Nechako Basin",
       subtitle = "Blue: hydrological (unregulated)  |  Red: operational (managed)  |  Grey: 2022–2025",
       x        = NULL, y = "SWSI", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title      = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())
ggsave(file.path(MAIN_OUTPUT_DIR, "comparison", "hswsi_vs_oswsi_full.png"),
       p_compare, width = 14, height = 5, dpi = 300)
cat("✓ Track comparison plot (H-SWSI vs O-SWSI) saved\n")

# ── Focus period: 2022–2025 ───────────────────────────────────────────────────
focus_df <- compare_df %>%
  filter(date >= FOCUS_START, date <= FOCUS_END, !is.na(swsi))

if (nrow(focus_df) > 0) {
  p_focus <- ggplot(focus_df, aes(x = date, y = swsi, colour = track)) +
    geom_hline(yintercept = c(0, -1, -2, -3, -4),
               linetype = c("solid","dashed","dashed","dashed","dashed"),
               colour   = c("grey50","#fc8d59","#d73027","#a50026","#4c0000"),
               linewidth = 0.4) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c("H-SWSI  (Climate-driven)"    = "#2c7bb6",
                                   "O-SWSI  (Operational supply)" = "#d7191c")) +
    scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
    scale_y_continuous(limits = c(-SWSI_CAP - 0.15, SWSI_CAP + 0.15),
                       breaks = seq(-4, 4, 1)) +
    labs(title    = "H-SWSI vs O-SWSI  —  2022–2025 Drought Focus  |  Nechako Basin",
         subtitle = "Garen (1993) Revised Z-score  |  Blue = climate signal  |  Red = operational stress",
         x        = NULL, y = "SWSI", colour = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title      = element_text(face = "bold"),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
  ggsave(file.path(MAIN_OUTPUT_DIR, "comparison", "focus_2022_2025_both_tracks.png"),
         p_focus, width = 14, height = 6, dpi = 300)
  cat("✓ 2022–2025 focus comparison plot saved\n")
}

# ── Diagnostic: natural lake storage time series ──────────────────────────────
if (!is.null(nat_lake_df)) {
  p_lakes <- ggplot(nat_lake_df, aes(x = date, y = storage_Mm3, colour = lake_name)) +
    geom_line(linewidth = 0.8) +
    scale_colour_brewer(palette = "Set2") +
    labs(title    = "Natural Lake Storage  —  GloLakes v2 (H-SWSI LakeNat component)",
         subtitle = "Stuart Lake (Nechako) + Fraser Lake (Nechako)",
         x = NULL, y = "Storage (Mm³)", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(MAIN_OUTPUT_DIR, "hswsi", "components", "natural_lake_storage.png"),
         p_lakes, width = 12, height = 5, dpi = 300)
  cat("✓ Natural lake storage diagnostic plot saved\n")
}

# ── Diagnostic: reservoir active storage vs total ─────────────────────────────
p_res_adj <- ggplot(reservoir_data_adj, aes(x = date)) +
  geom_line(aes(y = storage,            colour = "Total (GloLakes)"),    linewidth = 0.7) +
  geom_line(aes(y = active_storage_Mm3, colour = "Active (adjusted)"),  linewidth = 0.8) +
  geom_hline(yintercept = KENNEY_DEAD_STORAGE_MM3 + ENV_FLOW_RESERVE_MM3,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  annotate("text",
           x      = min(reservoir_data_adj$date) + 100,
           y      = KENNEY_DEAD_STORAGE_MM3 + ENV_FLOW_RESERVE_MM3 + 50,
           label  = sprintf("Dead storage + env reserve = %.0f Mm³",
                            KENNEY_DEAD_STORAGE_MM3 + ENV_FLOW_RESERVE_MM3),
           size   = 3, hjust = 0, colour = "grey40") +
  scale_colour_manual(values = c("Total (GloLakes)" = "#1f78b4",
                                 "Active (adjusted)" = "#e31a1c")) +
  labs(title    = "Kenney Reservoir Storage  —  Total vs Active (O-SWSI)",
       subtitle = sprintf("Active = Total − %.0f Mm³ dead − %.0f Mm³ env-flow reserve",
                          KENNEY_DEAD_STORAGE_MM3, ENV_FLOW_RESERVE_MM3),
       x = NULL, y = "Storage (Mm³)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(MAIN_OUTPUT_DIR, "oswsi", "components", "reservoir_total_vs_active.png"),
       p_res_adj, width = 12, height = 5, dpi = 300)
cat("✓ Reservoir active-storage diagnostic plot saved\n")

# ============================================================================
# SECTION 15: SAVE MASTER RDS & SUMMARY REPORT
# ============================================================================
all_results <- list(
  hswsi           = hswsi_result,
  oswsi           = oswsi_result,
  hswsi_droughts  = hswsi_droughts,
  oswsi_droughts  = oswsi_droughts,
  era5_data       = era5_data,
  reservoir_data  = reservoir_data_adj,
  nat_lake_data   = nat_lake_df,
  nat_lake_agg    = nat_lake_agg,
  stations        = STATIONS,
  metadata        = list(
    index_type        = "SWSI — Dual-Track (H-SWSI + O-SWSI)",
    reference         = "Garen (1993) J.WRPM 119(4):437-454  |  Z-score revision",
    basin             = "Nechako River Basin, BC, Canada",
    formula           = "SWSI = Σ(Wi·Zi),  Zi = Φ⁻¹(Pi/100),  capped ±4.2",
    hswsi_components  = c("Qnat","LakeNat","Precip","SWE"),
    hswsi_weights     = HSWSI_WEIGHTS,
    hswsi_baseline    = sprintf("Pre-regulation through %s (Qnat); full record (lake/climate)",
                                HSWSI_BASELINE_END),
    oswsi_components  = c("Active_Storage","SWE","Qin"),
    oswsi_weights     = OSWSI_WEIGHTS,
    oswsi_baseline    = sprintf("Post-regulation from %s", OSWSI_BASELINE_START),
    reservoir_adj_Mm3 = sprintf("dead=%.0f, env_flow=%.0f",
                                KENNEY_DEAD_STORAGE_MM3, ENV_FLOW_RESERVE_MM3),
    reservoir_source  = "GloLakes v2 (Hou et al. 2024)",
    era5_variables    = c("sd→SWE", "tp→Precip"),
    natural_lakes     = paste(NATURAL_LAKES$name, collapse=", "),
    qnat_station      = HSWSI_QNAT_PRIMARY,
    qin_station       = OSWSI_QIN_STATION,
    processed_date    = Sys.Date(),
    script_version    = "4.0_DualTrack_Zscore"
  )
)
saveRDS(all_results, file.path(MAIN_OUTPUT_DIR, "all_results_dual_track.rds"))
cat("✓ Master RDS saved: all_results_dual_track.rds\n")

# ── Summary report ────────────────────────────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 70)))
cat("DUAL-TRACK SWSI SUMMARY  —  Nechako Basin\n")
cat(sprintf("%s\n", strrep("=", 70)))
cat(sprintf("Formula     : SWSI = Σ(Wi·Zi),  Zi = Φ⁻¹(Pi/100),  scale ±%.1f\n", SWSI_CAP))
cat("Reference   : Garen (1993) J.WRPM 119(4):437-454  [Z-score revision]\n\n")

cat("─── TRACK 1  H-SWSI  (Hydrological Drought Index) ────────────────────\n")
cat(sprintf("  Objective  : Isolate climate-driven anomalies (unregulated basin)\n"))
cat(sprintf("  Qnat source: %s (naturalized)\n", HSWSI_QNAT_PRIMARY))
cat(sprintf("  Lakes      : %s\n", paste(NATURAL_LAKES$name, collapse=", ")))
cat(sprintf("  Weights    : Qnat=%.2f  LakeNat=%.2f  Precip=%.2f  SWE=%.2f\n",
            HSWSI_WEIGHTS["Qnat"], HSWSI_WEIGHTS["LakeNat"],
            HSWSI_WEIGHTS["Precip"], HSWSI_WEIGHTS["SWE"]))
if (sum(!is.na(hswsi_result$hswsi)) > 0)
  cat(sprintf("  Period     : %s – %s\n",
              min(hswsi_result$date[!is.na(hswsi_result$hswsi)]),
              max(hswsi_result$date[!is.na(hswsi_result$hswsi)])))
cat(sprintf("  Droughts   : %d events  |  Emergency (≤%.1f): %d\n",
            hswsi_droughts$n_events, TRIGGER_EMERGENCY, hswsi_droughts$n_emergency))

cat("\n─── TRACK 2  O-SWSI  (Operational Supply Drought Index) ──────────────\n")
cat(sprintf("  Objective  : Measure actual system stress & mgmt vulnerability\n"))
cat(sprintf("  Baseline   : Post-regulation from %s\n", OSWSI_BASELINE_START))
cat(sprintf("  Weights    : Storage=%.2f  SWE=%.2f  Qin=%.2f\n",
            OSWSI_WEIGHTS["Storage"], OSWSI_WEIGHTS["SWE"], OSWSI_WEIGHTS["Qin"]))
cat(sprintf("  Reservoir  : GloLakes v2  |  dead=%.0f Mm³  env-flow=%.0f Mm³\n",
            KENNEY_DEAD_STORAGE_MM3, ENV_FLOW_RESERVE_MM3))
cat(sprintf("  Qin source : %s naturalized (Kemano deduction: %.0f Mm³/mo)\n",
            OSWSI_QIN_STATION, KEMANO_BASELINE_DIV_MM3_MONTH))
if (sum(!is.na(oswsi_result$oswsi)) > 0)
  cat(sprintf("  Period     : %s – %s\n",
              min(oswsi_result$date[!is.na(oswsi_result$oswsi)]),
              max(oswsi_result$date[!is.na(oswsi_result$oswsi)])))
cat(sprintf("  Droughts   : %d events  |  Emergency (≤%.1f): %d\n",
            oswsi_droughts$n_events, TRIGGER_EMERGENCY, oswsi_droughts$n_emergency))

cat(sprintf("\n%s\n", strrep("=", 70)))
cat("SWSI Drought Scale:\n")
cat("   ≥ +2.0          Abundant supply\n")
cat("   0.0 to +2.0     Near normal\n")
cat("  -1.0 to  0.0     Below normal\n")
cat("  -2.0 to -1.0     Drought watch        (trigger ≤ -1.0)\n")
cat("  -3.0 to -2.0     Moderate drought     (trigger ≤ -2.0)\n")
cat("  -4.0 to -3.0     Severe drought\n")
cat("  < -4.0           Extreme drought\n")
cat(sprintf("%s\n", strrep("=", 70)))
cat(sprintf("\nAll outputs → %s/\n", MAIN_OUTPUT_DIR))
cat("  hswsi/             H-SWSI time series, component Z-scores, events\n")
cat("  oswsi/             O-SWSI time series, component Z-scores, events\n")
cat("  comparison/        H-SWSI vs O-SWSI overlay plots\n")
cat("  processed_data/    Adjusted reservoir storage CSV\n")
cat("  era5_extracted/    ERA5 basin-mean monthly CSV + SWE bias factors\n")
cat("  all_results_dual_track.rds  — master list for downstream analysis\n")
cat(sprintf("%s\n", strrep("=", 70)))