# ============================================================================
#  w9_atmospheric_diagnostics.R
#  Nechako Watershed Drought — ERA5 ANOMALY PRODUCER
#
#  ROLE (post-refactor):
#    This script is the SOLE producer of ERA5 anomaly NetCDF files and the
#    shared atmospheric state object consumed by w10a and w10c.  It does NOT
#    produce any figures or composite analyses — all that logic lives in w10a.
#
#  What this script does:
#    1. Defines all shared configuration (paths, climatology window, domain
#       boxes, drought focus period, EOF analysis domain).
#    2. Provides compute_anomaly() and load_or_compute_anomaly() which are
#       used here and re-exported via the shared-state RDS for reference.
#    3. Loads raw ERA5 monthly fields (Z500, SLP, SST).
#    4. Computes 1991-2020 anomalies for each field.
#    5. Saves the three anomaly NetCDFs to DATA_DIR so that w10a / w10c can
#       load them without repeating the expensive computation.
#    6. Builds the date spine, area-averaged time series, and map layers
#       (objects that are needed by multiple sections of w10a).
#    7. Saves a shared-state RDS (w9_shared_state.rds) to OUT_DIR that w10a
#       and w10c read with a single readRDS() call to restore all objects.
#
#  OUTPUTS written to DATA_DIR:
#    z500_monthly_anomaly_clim{CLIM_START}_{CLIM_END}.nc
#    slp_monthly_anomaly_clim{CLIM_START}_{CLIM_END}.nc
#    sst_monthly_anomaly_clim{CLIM_START}_{CLIM_END}.nc
#
#  OUTPUTS written to OUT_DIR:
#    w9_shared_state.rds   — list of shared objects (see Section 5)
#
#  Run BEFORE: w10a_further_atm_diag_era5_composites.R
#              w10c_eof_pca_circulation_indices.R
#  Dependencies: terra, lubridate, sf, rnaturalearth, rnaturalearthdata
#
#  CHANGE LOG
#  ----------
#  [w10c addition] Section 1: added EOF_LON_MIN/MAX, EOF_LAT_MIN/MAX constants.
#  [w10c addition] Section 5: added four EOF domain scalars to w9_shared list.
#  [w10c addition] Section 5: updated RDS header comment to document new keys.
#  All other code is unchanged from the previous version.
# ============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(lubridate)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

# ============================================================================
#  SECTION 1: CONFIGURATION
#  All constants used downstream by w10a and w10c are stored in the
#  shared-state RDS at the end of this script, so neither downstream script
#  ever needs to redeclare them.
# ============================================================================
cat("============================================================\n")
cat(" SECTION 1: Configuration\n")
cat("============================================================\n")

WD_PATH  <- "D:/Nechako_Drought/Nechako/"
DATA_DIR <- "D:/Nechako_Drought/Nechako/monthly_data_direct"
OUT_DIR  <- file.path(WD_PATH, "atmospheric_diagnostics")

setwd(WD_PATH)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Index CSV directories (must match DROUGHT_ANALYSIS_utils.R) ──────────────
SPI_SEAS_DIR  <- file.path(WD_PATH, "spi_results_seasonal/")
SPEI_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal/")

# Scales to attempt (scale 9 is skipped gracefully if no file found)
SPI_SCALES_TRY  <- c(1, 3, 6, 9, 12)
SPEI_SCALES_TRY <- c(1, 3, 6, 9, 12)

# ── Analysis period ──────────────────────────────────────────────────────────
START_YEAR <- 1950
END_YEAR   <- 2025

# ── Climatology reference window (WMO 1991-2020) ────────────────────────────
CLIM_START <- 1991
CLIM_END   <- 2020

# ── Drought event focus window ───────────────────────────────────────────────
RECENT_START        <- 2022
RECENT_END          <- END_YEAR
DROUGHT_FOCUS_START <- 2022L
DROUGHT_FOCUS_END   <- 2025L

# ── Spatial averaging boxes ──────────────────────────────────────────────────
# NW Canada ridge box (Gulf of Alaska action centre, PNA pattern ~55 N, 120 W)
RIDGE_LON_MIN <- -140; RIDGE_LON_MAX <- -110
RIDGE_LAT_MIN <-   50; RIDGE_LAT_MAX <-   65

# NE Pacific SST box (PDO warm tongue; excludes tropical Pacific to avoid
# aliasing ENSO variance that is handled separately via ONI in w7-w8)
SST_MEAN_LON_MIN <- -160; SST_MEAN_LON_MAX <- -130
SST_MEAN_LAT_MIN <-   40; SST_MEAN_LAT_MAX <-   55

# ── Plot settings (passed through shared state so w10a stays consistent) ─────
ANOMALY_PALETTE   <- "RdBu"
SST_PALETTE       <- "RdBu"
FIGURE_DPI        <- 200
FIGURE_WIDTH_WIDE <- 16
FIGURE_WIDTH_STD  <- 12
MONTH_LABELS      <- month.abb

# ── North Pacific EOF/PCA domain (0–360 longitude convention) ────────────────
# Used exclusively by w10c_eof_pca_circulation_indices.R.  Stored here so that
# w10c inherits them from the shared-state RDS without redeclaring.
#
# Domain choice: 140°E–120°W (= 240° in 0-360), 20°N–65°N.
# Rationale: spans the full North Pacific, capturing the PDO warm-tongue,
# the PNA ridge-trough couplet over the Gulf of Alaska, and ENSO-forced
# teleconnections that affect BC precipitation.
#
# IMPORTANT — longitude convention:
#   ERA5 rasters loaded by terra are stored in the 0-360 convention internally
#   (confirmed: DOM_XMIN = 100, DOM_XMAX = 250 in Section 6 below).
#   terra::ext(140, 240, 20, 65) therefore works WITHOUT any terra::rotate()
#   call.  Wrapping to -180:+180 is applied only for ggplot2 display in w10c
#   via a wrap_lon() helper (x > 180  →  x - 360).
EOF_LON_MIN <- 140L   # 140°E  in 0-360 convention
EOF_LON_MAX <- 240L   # 240°E  = 120°W in 0-360 convention
EOF_LAT_MIN <-  20L   # 20°N
EOF_LAT_MAX <-  65L   # 65°N

cat("  Configuration loaded.\n")
cat(sprintf("  EOF domain (0-360): lon [%d, %d], lat [%d, %d]\n",
            EOF_LON_MIN, EOF_LON_MAX, EOF_LAT_MIN, EOF_LAT_MAX))
cat("\n")

# ============================================================================
#  SECTION 2: ANOMALY COMPUTATION FUNCTIONS
#  compute_anomaly()         — vectorised raster anomaly computation
#  load_or_compute_anomaly() — cache-aware wrapper (skips expensive recompute
#                              if the anomaly NetCDF already exists on disk)
# ============================================================================
cat("============================================================\n")
cat(" SECTION 2: Anomaly computation functions\n")
cat("============================================================\n")

#' Subtract the 1991-2020 climatological monthly mean from every layer of a
#' SpatRaster.  terra broadcasts the single-layer climatology mean across all
#' layers of the same calendar month in one vectorised operation, avoiding
#' the O(n_years) layer-by-layer loop of older implementations.
#'
#' Date metadata is extracted from terra::time().  If the time axis is absent,
#' uses a non-CF calendar, or produces dates outside the range 1900-2030, the
#' function falls back to a positional monthly sequence starting Jan 1950 (the
#' ERA5 start year used throughout this project) and prints a warning so the
#' user can verify the assumption.  The fallback covers the common case of
#' ERA5 monthly files that were saved without explicit CF time metadata.
#'
#' @param r           SpatRaster.  terra::time() should ideally be set, but
#'                    a positional fallback is used automatically if not.
#' @param clim_start  First year of the reference climatology (default 1991).
#' @param clim_end    Last  year of the reference climatology (default 2020).
#' @return SpatRaster with the same number of layers as r, terra::time()
#'         preserved (or set to the fallback sequence), values = anomalies.
compute_anomaly <- function(r, clim_start = 1991, clim_end = 2020) {
  n_lyr <- nlyr(r)
  
  # ── Robust date extraction with positional fallback ───────────────────────
  # terra::time() can return NULL, a numeric 0-epoch vector, or dates outside
  # any plausible range when the NetCDF time axis is missing or uses an
  # unrecognised calendar (common for pre-processed ERA5 monthly files saved
  # without explicit CF time metadata).  In that case we reconstruct the date
  # sequence positionally: monthly from Jan 1950, which is the ERA5 start year
  # used throughout this project.  The fallback is announced so the user can
  # verify the assumption is correct for their data.
  dates_r <- tryCatch({
    d <- as.Date(time(r))
    if (is.null(d) || length(d) != n_lyr || all(is.na(d)) ||
        min(d, na.rm = TRUE) < as.Date("1900-01-01") ||
        max(d, na.rm = TRUE) > as.Date("2030-01-01")) {
      NULL   # signal: use fallback
    } else {
      d
    }
  }, error = function(e) NULL)
  
  if (is.null(dates_r)) {
    cat(sprintf(paste0(
      "  [warning] terra::time() did not return valid dates for this raster\n",
      "            (%d layers).  Assigning positional monthly sequence:\n",
      "            Jan 1950 + %d months.\n"), n_lyr, n_lyr))
    dates_r <- seq(as.Date("1950-01-01"), by = "month", length.out = n_lyr)
  }
  
  month_r  <- as.integer(format(dates_r, "%m"))
  year_r   <- as.integer(format(dates_r, "%Y"))
  clim_idx <- which(year_r >= clim_start & year_r <= clim_end)
  
  # ── One climatology mean per calendar month ───────────────────────────────
  anom_parts <- vector("list", 12)
  for (m in 1:12) {
    m_all  <- which(month_r == m)
    m_clim <- intersect(m_all, clim_idx)
    if (!length(m_all) || !length(m_clim)) next
    clim_m          <- app(r[[m_clim]], mean, na.rm = TRUE)  # 1-layer mean
    anom_parts[[m]] <- list(idx = m_all, anom = r[[m_all]] - clim_m)
    # terra broadcasts clim_m (1 layer) across all layers of r[[m_all]]
  }
  
  # ── Guard: at least one month must have been processed ───────────────────
  valid <- !vapply(anom_parts, is.null, logical(1))
  if (!any(valid))
    stop(
      "compute_anomaly: no calendar months could be processed.\n",
      "  Likely cause: the raster time axis could not be parsed and the\n",
      "  positional fallback date sequence does not overlap the climatology\n",
      "  window (", clim_start, "-", clim_end, ").\n",
      "  Check START_YEAR, CLIM_START, CLIM_END and the number of raster layers."
    )
  
  # ── Reassemble in original chronological order ────────────────────────────
  # use.names = FALSE + explicit as.integer() guarantees a plain integer
  # vector.  Without this, unlist() on an empty list returns NULL (not a
  # vector), and unlist() on named-integer elements can return a named
  # integer — both cause order() to throw "argument 1 is not a vector".
  all_idx  <- as.integer(unlist(lapply(anom_parts[valid], `[[`, "idx"),
                                use.names = FALSE))
  all_anom <- do.call(c, lapply(anom_parts[valid], `[[`, "anom"))
  result   <- all_anom[[order(all_idx)]]
  terra::time(result) <- dates_r
  result
}

#' Load a pre-computed anomaly NetCDF if it exists; otherwise compute from the
#' raw ERA5 file, apply optional unit conversion, and save to anomaly_nc so
#' that subsequent runs skip the computation entirely.
#'
#' @param raw_nc       Path to the raw ERA5 monthly NetCDF.
#' @param anomaly_nc   Path where the anomaly NetCDF is/will be stored.
#' @param unit_divisor Divide raw values by this before computing anomalies
#'                     (e.g. 9.80665 to convert geopotential J/kg → m;
#'                           100     to convert Pa → hPa).
#' @param unit_offset  Add this constant after dividing
#'                     (e.g. -273.15 to convert K → °C).
#' @param clim_start   First year of reference climatology.
#' @param clim_end     Last  year of reference climatology.
#' @return SpatRaster of anomalies with terra::time() set.
load_or_compute_anomaly <- function(raw_nc, anomaly_nc,
                                    unit_divisor = 1,
                                    unit_offset  = 0,
                                    clim_start   = CLIM_START,
                                    clim_end     = CLIM_END) {
  if (file.exists(anomaly_nc)) {
    cat(sprintf("  [cache] Reusing existing anomaly: %s\n", basename(anomaly_nc)))
    return(rast(anomaly_nc))
  }
  if (!file.exists(raw_nc))
    stop(sprintf("Raw ERA5 file not found: %s", raw_nc))
  cat(sprintf("  [compute] Loading raw file: %s\n", basename(raw_nc)))
  r <- rast(raw_nc)
  if (unit_divisor != 1) r <- r / unit_divisor
  if (unit_offset  != 0) r <- r + unit_offset
  cat(sprintf("  [compute] Computing %d-%d anomalies (%d layers)... ",
              clim_start, clim_end, nlyr(r)))
  anom <- compute_anomaly(r, clim_start, clim_end)
  cat("Done.\n")
  cat(sprintf("  [save]    Writing: %s\n", basename(anomaly_nc)))
  terra::writeCDF(anom, anomaly_nc, overwrite = TRUE)
  anom
}

cat("  Functions defined.\n\n")

# ============================================================================
#  SECTION 3: LOAD / COMPUTE ERA5 ANOMALIES
#  Canonical file names used by both w9 and w10a.  If the anomaly file
#  already exists from a previous run it is loaded directly; otherwise it is
#  computed from the raw ERA5 file and saved for future reuse.
# ============================================================================
cat("============================================================\n")
cat(" SECTION 3: Loading / computing ERA5 anomalies\n")
cat("============================================================\n")

# ── Canonical anomaly file paths (used by w10a and w10c as well) ──────────────
z500_nc_path <- file.path(DATA_DIR,
                          sprintf("z500_monthly_anomaly_clim%d_%d.nc", CLIM_START, CLIM_END))
slp_nc_path  <- file.path(DATA_DIR,
                          sprintf("slp_monthly_anomaly_clim%d_%d.nc",  CLIM_START, CLIM_END))
sst_nc_path  <- file.path(DATA_DIR,
                          sprintf("sst_monthly_anomaly_clim%d_%d.nc",  CLIM_START, CLIM_END))

# ── Raw ERA5 file paths ───────────────────────────────────────────────────────
# The raw file for Z500 stores geopotential in J/kg; divide by g=9.80665 for
# geopotential height in metres.
# The raw file for SLP stores pressure in Pa; divide by 100 for hPa.
# The raw file for SST stores temperature in K; subtract 273.15 for °C.
# Alternative raw file names are tried in order to accommodate ERA5 naming
# conventions from different download workflows.
find_raw <- function(...) {
  candidates <- c(...)
  found <- candidates[file.exists(candidates)]
  if (!length(found))
    stop("None of these raw ERA5 files were found:\n  ",
         paste(candidates, collapse = "\n  "))
  found[1]
}

z500_raw <- find_raw(
  file.path(DATA_DIR, "geopotential_500hPa_monthly.nc"),
  file.path(DATA_DIR, "z500_monthly.nc"))

slp_raw_path <- find_raw(
  file.path(DATA_DIR, "mean_sea_level_pressure_monthly.nc"),
  file.path(DATA_DIR, "slp_monthly.nc"))

sst_raw <- find_raw(
  file.path(DATA_DIR, "sea_surface_temperature_monthly.nc"),
  file.path(DATA_DIR, "sst_monthly.nc"))

# ── Detect SLP unit: Pa (needs /100) or already hPa? ─────────────────────────
# ERA5 downloads via the CDS API deliver Pa.  Some pre-processed local files
# are already in hPa.  A quick range check decides which case applies.
slp_unit_divisor <- local({
  r_tmp  <- rast(slp_raw_path)[[1]]
  v_med  <- median(values(r_tmp), na.rm = TRUE)
  if (v_med > 10000) 100 else 1  # Pa → hPa if median > 10000
})
cat(sprintf("  SLP unit divisor: %g  (median raw value = %.0f)\n",
            slp_unit_divisor,
            median(values(rast(slp_raw_path)[[1]]), na.rm = TRUE)))

# ── Compute / load all three anomaly fields ───────────────────────────────────
cat("\n  Z500:\n")
z500_anom <- load_or_compute_anomaly(
  raw_nc       = z500_raw,
  anomaly_nc   = z500_nc_path,
  unit_divisor = 9.80665)

cat("\n  SLP:\n")
slp_anom  <- load_or_compute_anomaly(
  raw_nc       = slp_raw_path,
  anomaly_nc   = slp_nc_path,
  unit_divisor = slp_unit_divisor)

cat("\n  SST:\n")
sst_anom  <- load_or_compute_anomaly(
  raw_nc      = sst_raw,
  anomaly_nc  = sst_nc_path,
  unit_offset = -273.15)

cat(sprintf("\n  Anomaly layers — Z500: %d | SLP: %d | SST: %d\n",
            nlyr(z500_anom), nlyr(slp_anom), nlyr(sst_anom)))

# ============================================================================
#  SECTION 4: BUILD SHARED ATMOSPHERIC STATE OBJECTS
#  These objects are computationally inexpensive to create but are needed by
#  multiple sections of w10a.  Building them here (once) and persisting them
#  avoids duplicating logic across scripts.
# ============================================================================
cat("============================================================\n")
cat(" SECTION 4: Building shared state objects\n")
cat("============================================================\n")

# ── Date spine aligned to ERA5 raster layers ─────────────────────────────────
all_dates <- seq.Date(as.Date(sprintf("%d-01-01", START_YEAR)),
                      as.Date(sprintf("%d-12-01", END_YEAR)),
                      by = "month")
n_total <- length(all_dates)
if (nlyr(z500_anom) != n_total)
  warning(sprintf("Z500 has %d layers but date spine expects %d.",
                  nlyr(z500_anom), n_total))

base_date_index <- data.frame(
  layer  = seq_len(n_total),
  date   = all_dates,
  year   = as.integer(format(all_dates, "%Y")),
  month  = as.integer(format(all_dates, "%m")),
  season = ifelse(as.integer(format(all_dates, "%m")) %in% c(12, 1, 2),  "DJF",
                  ifelse(as.integer(format(all_dates, "%m")) %in% c(3,  4, 5),  "MAM",
                         ifelse(as.integer(format(all_dates, "%m")) %in% c(6,  7, 8),  "JJA",
                                "SON"))),
  stringsAsFactors = FALSE
)
cat(sprintf("  Date spine: %s to %s (%d months)\n",
            min(all_dates), max(all_dates), n_total))

# ── Area-averaged time series (ridge box and NE Pacific SST box) ──────────────
# Extracted once here; w10a uses them throughout Sections 4, 4b, and 5
# without touching the full rasters again.
ridge_ext <- terra::ext(RIDGE_LON_MIN, RIDGE_LON_MAX, RIDGE_LAT_MIN, RIDGE_LAT_MAX)
sst_ext   <- terra::ext(SST_MEAN_LON_MIN, SST_MEAN_LON_MAX,
                        SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX)

z500_ridge_ts <- terra::global(terra::crop(z500_anom, ridge_ext),
                               "mean", na.rm = TRUE)[, 1]
slp_nwbc_ts   <- terra::global(terra::crop(slp_anom,  ridge_ext),
                               "mean", na.rm = TRUE)[, 1]
sst_nepac_ts  <- terra::global(terra::crop(sst_anom,  sst_ext),
                               "mean", na.rm = TRUE)[, 1]
cat("  Area-averaged time series extracted.\n")

# ── Map layers ────────────────────────────────────────────────────────────────
cat("  Loading map layers... ")
LAKES_SHP <- "D:/Nechako_Drought/Nechako/Spatial/ne_50m_lakes/ne_50m_lakes.shp"
borders <- ne_countries(scale = "medium", returnclass = "sf")
lakes   <- tryCatch(
  sf::st_read(LAKES_SHP, quiet = TRUE),
  error = function(e) {
    warning(sprintf(
      "[w9] Local lakes shapefile not found: %s\n  Continuing without lake polygons.",
      LAKES_SHP))
    sf::st_sf(geometry = sf::st_sfc(crs = 4326))
  })
map_layers <- list(borders = borders, lakes = lakes)
cat("Done.\n\n")

# ============================================================================
#  SECTION 5: SAVE SHARED STATE RDS
#  w10a and w10c read this single file to restore all shared objects without
#  any ERA5 computation, heavy I/O, or redeclaration of configuration
#  constants.
#
#  The RDS contains:
#    base_date_index  — full date spine data.frame (n_months × 5)
#    z500_ridge_ts    — numeric vector (n_months): Z500 anom ridge-box mean
#    slp_nwbc_ts      — numeric vector (n_months): SLP anom ridge-box mean
#    sst_nepac_ts     — numeric vector (n_months): SST anom NE-Pac box mean
#    map_layers       — list(borders = sf, lakes = sf)
#    anomaly_nc_z500  — character: full path to z500 anomaly NetCDF
#    anomaly_nc_slp   — character: full path to slp anomaly NetCDF
#    anomaly_nc_sst   — character: full path to sst anomaly NetCDF
#    + all scalar configuration constants (START_YEAR ... MONTH_LABELS)
#    EOF_LON_MIN      — integer: western bound of EOF domain (0-360, °E)  [w10c]
#    EOF_LON_MAX      — integer: eastern bound of EOF domain (0-360, °E)  [w10c]
#    EOF_LAT_MIN      — integer: southern bound of EOF domain (°N)        [w10c]
#    EOF_LAT_MAX      — integer: northern bound of EOF domain (°N)        [w10c]
# ============================================================================
cat("============================================================\n")
cat(" SECTION 5: Saving shared state RDS\n")
cat("============================================================\n")

w9_shared <- list(
  # Shared data objects
  base_date_index  = base_date_index,
  z500_ridge_ts    = z500_ridge_ts,
  slp_nwbc_ts      = slp_nwbc_ts,
  sst_nepac_ts     = sst_nepac_ts,
  map_layers       = map_layers,
  # Anomaly file paths (w10a and w10c load these directly)
  anomaly_nc_z500  = z500_nc_path,
  anomaly_nc_slp   = slp_nc_path,
  anomaly_nc_sst   = sst_nc_path,
  # Configuration scalars
  WD_PATH          = WD_PATH,
  DATA_DIR         = DATA_DIR,
  OUT_DIR          = OUT_DIR,
  SPI_SEAS_DIR     = SPI_SEAS_DIR,
  SPEI_SEAS_DIR    = SPEI_SEAS_DIR,
  SPI_SCALES_TRY   = SPI_SCALES_TRY,
  SPEI_SCALES_TRY  = SPEI_SCALES_TRY,
  START_YEAR       = START_YEAR,
  END_YEAR         = END_YEAR,
  CLIM_START       = CLIM_START,
  CLIM_END         = CLIM_END,
  RECENT_START     = RECENT_START,
  RECENT_END       = RECENT_END,
  DROUGHT_FOCUS_START = DROUGHT_FOCUS_START,
  DROUGHT_FOCUS_END   = DROUGHT_FOCUS_END,
  RIDGE_LON_MIN    = RIDGE_LON_MIN,
  RIDGE_LON_MAX    = RIDGE_LON_MAX,
  RIDGE_LAT_MIN    = RIDGE_LAT_MIN,
  RIDGE_LAT_MAX    = RIDGE_LAT_MAX,
  SST_MEAN_LON_MIN = SST_MEAN_LON_MIN,
  SST_MEAN_LON_MAX = SST_MEAN_LON_MAX,
  SST_MEAN_LAT_MIN = SST_MEAN_LAT_MIN,
  SST_MEAN_LAT_MAX = SST_MEAN_LAT_MAX,
  ANOMALY_PALETTE  = ANOMALY_PALETTE,
  SST_PALETTE      = SST_PALETTE,
  FIGURE_DPI       = FIGURE_DPI,
  FIGURE_WIDTH_WIDE = FIGURE_WIDTH_WIDE,
  FIGURE_WIDTH_STD  = FIGURE_WIDTH_STD,
  MONTH_LABELS     = MONTH_LABELS,
  # ── EOF analysis domain (used by w10c) ─────────────────────────────────────
  # 0-360 longitude convention — matches how terra stores ERA5 rasters.
  # terra::ext(EOF_LON_MIN, EOF_LON_MAX, EOF_LAT_MIN, EOF_LAT_MAX) works
  # directly without terra::rotate().  See w10c wrap_lon() for display.
  EOF_LON_MIN      = EOF_LON_MIN,
  EOF_LON_MAX      = EOF_LON_MAX,
  EOF_LAT_MIN      = EOF_LAT_MIN,
  EOF_LAT_MAX      = EOF_LAT_MAX
)

rds_path <- file.path(OUT_DIR, "w9_shared_state.rds")
saveRDS(w9_shared, rds_path)
cat(sprintf("  Shared state saved: %s\n", rds_path))

# ============================================================================
#  SECTION 6: MANUSCRIPT FIGURE 7
#  Three-panel composite anomaly map (Z500 | SLP | SST) for the 2022–2025
#  drought period, plotted over the North Pacific–North America domain.
#
#  Panel design:
#   (a) Z500 anomaly (m)  — diverging RdBu palette; stippling for months
#       where the anomaly is statistically significant vs the 1950–2025
#       background (pointwise two-sided t-test, α = 0.05); ridge-box outline.
#   (b) SLP anomaly (hPa) — same palette; ridge-box outline for context.
#   (c) SST anomaly (°C)  — same palette; annotated PDO warm-tongue box;
#       contour lines at −0.5 and +0.5 °C to guide the eye.
#
#  Domain: 120°E–90°W (100°E–70°W as buffer), 15°N–75°N.
#  The composite is the mean over all months in DROUGHT_FOCUS_START–
#  DROUGHT_FOCUS_END (inclusive); the significance test uses the full
#  1950–CLIM_END monthly background.
#
#  Outputs (OUT_DIR):
#    Fig7_MS_AtmosphericComposite_2022_2025.pdf   (vector, 7.5 × 8.5 in)
#    Fig7_MS_AtmosphericComposite_2022_2025.png   (300 DPI)
# ============================================================================
cat("============================================================\n")
cat(" SECTION 6: Manuscript Figure 7 — atmospheric composite\n")
cat("============================================================\n")

tryCatch({
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(scales)
  })
  
  # ── domain ────────────────────────────────────────────────────────────────
  DOM_XMIN <- 100; DOM_XMAX <- 250   # 100°E–110°W in 0–360 convention
  DOM_YMIN <-  15; DOM_YMAX <-  75
  # In -180/+180 convention for sf / ggplot2:
  DOM_XMIN_WGS <- -260 + 360; DOM_XMAX_WGS <- 270 - 360  # handled via wrap below
  
  # ── layer index for the drought focus period ───────────────────────────────
  focus_idx <- which(base_date_index$year >= DROUGHT_FOCUS_START &
                       base_date_index$year <= DROUGHT_FOCUS_END)
  cat(sprintf("  Composite: %d months (%d-%d)\n",
              length(focus_idx), DROUGHT_FOCUS_START, DROUGHT_FOCUS_END))
  
  # ── helper: compute composite mean + pointwise t-test ─────────────────────
  make_composite <- function(anom_rast, focus_lyr, label) {
    n_focus <- length(focus_lyr)
    n_total <- terra::nlyr(anom_rast)
    
    # composite mean (focus period)
    comp_mean <- terra::app(anom_rast[[focus_lyr]], mean, na.rm = TRUE)
    
    # background mean and SD (full record)
    bg_mean <- terra::app(anom_rast, mean, na.rm = TRUE)
    bg_sd   <- terra::app(anom_rast, sd,   na.rm = TRUE)
    
    # pointwise t-statistic: is focus mean significantly ≠ background mean?
    # t = (comp_mean - bg_mean) / (bg_sd / sqrt(n_focus))
    t_stat  <- (comp_mean - bg_mean) /
      (bg_sd / sqrt(as.numeric(n_focus)))
    t_crit  <- qt(0.975, df = n_focus - 1)
    sig_rast <- abs(t_stat) > t_crit   # TRUE = significant at α=0.05
    
    cat(sprintf("  %s: composite range [%.2f, %.2f]\n", label,
                min(terra::values(comp_mean), na.rm = TRUE),
                max(terra::values(comp_mean), na.rm = TRUE)))
    
    list(mean = comp_mean, sig = sig_rast, t = t_stat)
  }
  
  z500_comp <- make_composite(z500_anom, focus_idx, "Z500")
  slp_comp  <- make_composite(slp_anom,  focus_idx, "SLP")
  sst_comp  <- make_composite(sst_anom,  focus_idx, "SST")
  
  # ── helper: raster → long data frame, WGS84 coords ─────────────────────────
  rast_to_df <- function(r, val_name = "value") {
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(df)[3] <- val_name
    # wrap longitudes from [0,360] to [-180,180] if needed
    if (max(df$x, na.rm = TRUE) > 181)
      df$x <- ifelse(df$x > 180, df$x - 360, df$x)
    df
  }
  
  # composite mean data frames
  df_z500 <- rast_to_df(z500_comp$mean, "Z500_anom")
  df_slp  <- rast_to_df(slp_comp$mean,  "SLP_anom")
  df_sst  <- rast_to_df(sst_comp$mean,  "SST_anom")
  
  # significance points (keep only TRUE pixels; low-res: every 3rd pixel)
  sig_to_pts <- function(sig_rast, stride = 3) {
    df <- as.data.frame(sig_rast, xy = TRUE, na.rm = FALSE)
    names(df)[3] <- "sig"
    df <- df[!is.na(df$sig) & df$sig == 1, ]
    if (max(df$x, na.rm = TRUE) > 181)
      df$x <- ifelse(df$x > 180, df$x - 360, df$x)
    # thin to avoid overplotting
    df[seq(1, nrow(df), by = stride), ]
  }
  pts_z500 <- sig_to_pts(z500_comp$sig)
  pts_slp  <- sig_to_pts(slp_comp$sig)
  
  # SST contour data (for ±0.5°C lines)
  sst_for_contour <- sst_comp$mean
  df_sst_cont <- rast_to_df(sst_for_contour, "SST_anom")
  
  # ── shared map infrastructure ─────────────────────────────────────────────
  # crop borders and lakes to domain
  DOMAIN_RECT <- sf::st_bbox(c(xmin = -260, xmax = -70,
                               ymin = DOM_YMIN, ymax = DOM_YMAX),
                             crs = sf::st_crs(4326))
  crop_sf <- function(layer) {
    tryCatch(
      sf::st_crop(layer, DOMAIN_RECT),
      error   = function(e) layer,
      warning = function(w) suppressWarnings(sf::st_crop(layer, DOMAIN_RECT)))
  }
  borders_crop <- crop_sf(map_layers$borders)
  lakes_crop   <- tryCatch(crop_sf(map_layers$lakes),
                           error = function(e) NULL)
  
  # ── shared theme ──────────────────────────────────────────────────────────
  theme_atm <- ggplot2::theme_void(base_size = 8.5) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 8.5, face = "bold",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 7, colour = "grey35",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 2)),
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.40, "cm"),
      legend.text      = ggplot2::element_text(size = 6.5),
      legend.title     = ggplot2::element_text(size = 7, face = "bold"),
      panel.border     = ggplot2::element_rect(colour = "grey40",
                                               fill = NA, linewidth = 0.35),
      plot.margin      = ggplot2::margin(3, 5, 3, 5)
    )
  
  # ridge-box rectangle (used on Z500 and SLP panels)
  ridge_rect <- data.frame(
    xmin = RIDGE_LON_MIN, xmax = RIDGE_LON_MAX,
    ymin = RIDGE_LAT_MIN, ymax = RIDGE_LAT_MAX
  )
  # PDO warm-tongue box (used on SST panel)
  pdo_rect <- data.frame(
    xmin = SST_MEAN_LON_MIN, xmax = SST_MEAN_LON_MAX,
    ymin = SST_MEAN_LAT_MIN, ymax = SST_MEAN_LAT_MAX
  )
  
  # colour limits (symmetric, clipped)
  clamp <- function(x, pct = 0.98) {
    q <- quantile(abs(x), pct, na.rm = TRUE)
    c(-q, q)
  }
  lim_z500 <- clamp(df_z500$Z500_anom)
  lim_slp  <- clamp(df_slp$SLP_anom)
  lim_sst  <- clamp(df_sst$SST_anom)
  
  # ── geom_raster + map layers helper ──────────────────────────────────────
  base_map_layers <- function() {
    list(
      ggplot2::geom_sf(data        = borders_crop,
                       fill        = NA,
                       colour      = "grey30",
                       linewidth   = 0.25,
                       inherit.aes = FALSE),
      if (!is.null(lakes_crop) && nrow(lakes_crop) > 0)
        ggplot2::geom_sf(data        = lakes_crop,
                         fill        = "white",
                         colour      = "grey50",
                         linewidth   = 0.15,
                         inherit.aes = FALSE)
      else NULL
    )
  }
  
  map_coord <- function() {
    ggplot2::coord_sf(xlim   = c(-260, -70),
                      ylim   = c(DOM_YMIN, DOM_YMAX),
                      expand = FALSE,
                      crs    = sf::st_crs(4326))
  }
  
  # ════════════════════════════════════════════════════════════════════════════
  # PANEL (a): Z500 anomaly
  # ════════════════════════════════════════════════════════════════════════════
  p_z500 <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = df_z500,
      ggplot2::aes(x = x, y = y, fill = Z500_anom)) +
    ggplot2::scale_fill_distiller(
      palette  = "RdBu",
      direction = -1,
      limits   = lim_z500,
      oob      = scales::squish,
      name     = "Z500 anomaly (m)",
      guide    = ggplot2::guide_colorbar(
        barwidth       = 7, barheight = 0.45,
        title.position = "top")) +
    
    # stippling: significant pixels
    {if (nrow(pts_z500) > 0)
      ggplot2::geom_point(
        data = pts_z500,
        ggplot2::aes(x = x, y = y),
        shape = 20, size = 0.4,
        colour = "grey15", alpha = 0.65,
        inherit.aes = FALSE)
      else NULL} +
    
    # ridge action-centre box
    ggplot2::geom_rect(
      data = ridge_rect,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, colour = "#d73027", linewidth = 0.55,
      linetype = "dashed", inherit.aes = FALSE) +
    ggplot2::annotate("text",
                      x = RIDGE_LON_MIN + 1, y = RIDGE_LAT_MAX + 1.5,
                      label = "Ridge\ncentre", size = 2.0,
                      hjust = 0, colour = "#d73027", fontface = "bold") +
    
    base_map_layers() +
    map_coord() +
    
    ggplot2::labs(
      title    = sprintf("(a)  Z500 anomaly \u2014 %d\u2013%d composite",
                         DROUGHT_FOCUS_START, DROUGHT_FOCUS_END),
      subtitle = "vs WMO 1991\u20132020 climatology. Stippling: p < 0.05 (pointwise t-test)") +
    theme_atm
  
  # ════════════════════════════════════════════════════════════════════════════
  # PANEL (b): SLP anomaly
  # ════════════════════════════════════════════════════════════════════════════
  p_slp <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = df_slp,
      ggplot2::aes(x = x, y = y, fill = SLP_anom)) +
    ggplot2::scale_fill_distiller(
      palette   = "RdBu",
      direction = -1,
      limits    = lim_slp,
      oob       = scales::squish,
      name      = "SLP anomaly (hPa)",
      guide     = ggplot2::guide_colorbar(
        barwidth = 7, barheight = 0.45,
        title.position = "top")) +
    
    # stippling
    {if (nrow(pts_slp) > 0)
      ggplot2::geom_point(
        data = pts_slp,
        ggplot2::aes(x = x, y = y),
        shape = 20, size = 0.4,
        colour = "grey15", alpha = 0.65,
        inherit.aes = FALSE)
      else NULL} +
    
    ggplot2::geom_rect(
      data = ridge_rect,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, colour = "#d73027", linewidth = 0.55,
      linetype = "dashed", inherit.aes = FALSE) +
    
    base_map_layers() +
    map_coord() +
    
    ggplot2::labs(
      title    = "(b)  SLP anomaly",
      subtitle = "vs WMO 1991\u20132020 climatology. Dashed box = ridge centre (Z500 panel a)") +
    theme_atm
  
  # ════════════════════════════════════════════════════════════════════════════
  # PANEL (c): SST anomaly
  # ════════════════════════════════════════════════════════════════════════════
  p_sst <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = df_sst,
      ggplot2::aes(x = x, y = y, fill = SST_anom)) +
    ggplot2::scale_fill_distiller(
      palette   = "RdBu",
      direction = -1,
      limits    = lim_sst,
      oob       = scales::squish,
      name      = "SST anomaly (\u00b0C)",
      guide     = ggplot2::guide_colorbar(
        barwidth = 7, barheight = 0.45,
        title.position = "top")) +
    
    # ±0.5°C contour lines
    ggplot2::geom_contour(
      data    = df_sst_cont,
      ggplot2::aes(x = x, y = y, z = SST_anom),
      breaks  = c(-0.5, 0.5),
      colour  = "grey20",
      linewidth = 0.35,
      linetype  = "solid",
      inherit.aes = FALSE) +
    
    # PDO warm-tongue box
    ggplot2::geom_rect(
      data = pdo_rect,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = NA, colour = "#2166ac", linewidth = 0.55,
      linetype = "dashed", inherit.aes = FALSE) +
    ggplot2::annotate("text",
                      x = SST_MEAN_LON_MIN + 1, y = SST_MEAN_LAT_MAX + 1.5,
                      label = "PDO warm\ntongue", size = 2.0,
                      hjust = 0, colour = "#2166ac", fontface = "bold") +
    
    base_map_layers() +
    map_coord() +
    
    ggplot2::labs(
      title    = "(c)  SST anomaly",
      subtitle = "vs WMO 1991\u20132020. Contours at \u00b10.5\u00b0C. Dashed box = PDO warm-tongue region.") +
    theme_atm
  
  # ── assemble ─────────────────────────────────────────────────────────────
  fig7 <- p_z500 / p_slp / p_sst +
    patchwork::plot_layout(guides = "keep", heights = c(1, 1, 1)) +
    patchwork::plot_annotation(
      title    = sprintf(
        "ERA5 atmospheric composite anomalies \u2014 Nechako drought period %d\u2013%d",
        DROUGHT_FOCUS_START, DROUGHT_FOCUS_END),
      subtitle = paste0(
        "Domain: 100\u00b0E\u2013110\u00b0W, 15\u00b0N\u201375\u00b0N.  ",
        "Anomalies relative to WMO 1991\u20132020 climatology.  ",
        "Stippling on Z500 and SLP: p < 0.05 (two-sided pointwise t-test vs 1950\u2013",
        CLIM_END, " background)."),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 10, face = "bold",
                                              hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey35",
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 4))
      )
    )
  
  # ── save ─────────────────────────────────────────────────────────────────
  fig7_pdf <- file.path(OUT_DIR, "Fig7_MS_AtmosphericComposite_2022_2025.pdf")
  fig7_png <- file.path(OUT_DIR, "Fig7_MS_AtmosphericComposite_2022_2025.png")
  
  tryCatch({
    ggplot2::ggsave(fig7_pdf, fig7,
                    width = 7.5, height = 8.5, units = "in", device = "pdf")
    cat(sprintf("  \u2713 Fig7 (PDF): %s\n", basename(fig7_pdf)))
  }, error = function(e) cat(sprintf("  \u26a0 Fig7 PDF: %s\n", e$message)))
  
  tryCatch({
    ggplot2::ggsave(fig7_png, fig7,
                    width = 7.5, height = 8.5, units = "in",
                    dpi = 300, device = "png")
    cat(sprintf("  \u2713 Fig7 (PNG): %s\n", basename(fig7_png)))
  }, error = function(e) cat(sprintf("  \u26a0 Fig7 PNG: %s\n", e$message)))
  
}, error = function(e) {
  cat(sprintf("  \u26a0 Section 6 failed: %s\n", e$message))
  cat("  Check that Sections 3\u20135 completed successfully before re-running.\n")
})

# ============================================================================
#  DONE
# ============================================================================
cat("============================================================\n")
cat(" w9 COMPLETE\n")
cat("============================================================\n")
cat(sprintf("  Anomaly NetCDFs written to: %s\n", DATA_DIR))
cat(sprintf("    %s\n", basename(z500_nc_path)))
cat(sprintf("    %s\n", basename(slp_nc_path)))
cat(sprintf("    %s\n", basename(sst_nc_path)))
cat(sprintf("  Shared state RDS: %s\n", rds_path))
cat("  NEW manuscript figure:\n")
cat("    Fig7_MS_AtmosphericComposite_2022_2025.pdf/.png\n")
cat("  Next steps:\n")
cat("    w10a_further_atm_diag_era5_composites.R  (composite maps)\n")
cat("    w10b_further_atm_diag_index_skill.R      (index skill, independent)\n")
cat("    w10c_eof_pca_circulation_indices.R        (EOF/PCA, reads w9 RDS)\n")
cat("============================================================\n")