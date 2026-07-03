##############################################
# SWEI / SSPI-1 CALCULATION - PARALLEL VERSION (SPI-STYLE OUTPUTS)
# Uses parallel::parLapply - Works on Windows and Linux
# Output naming matches SPI_ERALand.R convention
#
# METHOD (user-selected at runtime):
#   1 = Huning & AghaKouchak (2020) fixed 5% SCF threshold (SWEI)
#       Right-aligned 3-month climatology, matches original paper.
#   2 = SSPI-1: Standardized SnowPack Index (monthly, no SCF mask)
#       Zero-inflated gamma distribution, following Stagge et al. (2015) and JRC EDO (2020).
#       Output to: sspi_results_monthly/
#
# Other methodological alignments with Huning & AghaKouchak (2020):
#   - Zero-SWE perturbation: random uniform in (0, min_nonzero)  [method 1]
#   - Minimum data: >= 75% of years per calendar month must be nonzero [method 1]
#   - No clipping of SWEI values  [method 1]
#
# ------------------------------------------------------------------------------
# REGIONAL FREQUENCY ANALYSIS (RFA) EXTENSIONS  [method 2 / SSPI-1 only]
# ------------------------------------------------------------------------------
# Adapted from the index-flood RFA framework of Hosking & Wallis (1993, 1997),
# as applied to ERA5-Land snow variables by Fontana et al. (2026, Weather and
# Climate Extremes). Rather than fitting a zero-inflated gamma independently to
# each pixel's short, often zero-inflated reference-period record, pixels are:
#
#   1. REGIONALIZED into homogeneous sub-regions (elevation band x aspect class,
#      or a spatial fallback if no DEM is supplied) - see build_regions().
#   2. NORMALIZED by a per-pixel, per-calendar-month index value (median or mean
#      of that pixel's reference-period SWE for the month) - Eq. (1)-(2) analogue.
#   3. SCREENED for discordancy using Hosking & Wallis L-moment discordancy (Di)
#      and regional heterogeneity (H0, Monte-Carlo on a fitted kappa distribution)
#      before pooling - see discordancy_test() / heterogeneity_test(). Region/
#      months with H0 >= H0_REGION_FLAG_THRESHOLD (or an untested/failed H0)
#      are EXCLUDED from pooling entirely, not just flagged.
#   4. OPTIONALLY BIAS-CORRECTED against independent station SWE observations -
#      see the "BASIN-WIDE BIAS CORRECTION" section below and
#      fit_basin_bias_correction(). Applied ONCE, basin-wide, directly to
#      swe_matrix, before regionalization/gamma fitting even begins, so it
#      benefits both the regionalized and pixel-by-pixel paths identically.
#      (The original per-region quantile-mapping fit_bias_correction() is kept
#      in the script, defined but unused - see Section F of the RFA
#      construction block for why it was replaced.)
#   5. POOLED within each homogeneous (H0-passing) region/month to fit a single
#      regional gamma distribution (lmomco, L-moments), tested for goodness-of-
#      fit with a Kolmogorov-Smirnov test - see fit_regional_gamma() /
#      gof_test_gamma(). A region/month that fails GoF, or was excluded at the
#      H0 step, falls straight through to per-pixel local fitting
#      (monthly_sspi) - there is deliberately NO regional-empirical CDF
#      fallback between "regional gamma" and "local".
#
# This stabilizes the gamma-parameter estimates for short/zero-inflated pixel
# records (the regional sample size is far larger than any single pixel's record)
# while keeping each pixel's own zero-probability and index value, so the
# standardized index remains pixel-specific in its mean state.
#
# Toggle with REGIONALIZE_SSPI (TRUE/FALSE) in the CONFIGURATION block below.
# Set REGIONALIZE_SSPI <- FALSE to reproduce the original pixel-by-pixel logic.
# ------------------------------------------------------------------------------
#
# RUNNING THE SCRIPT:
#   Interactive (RStudio / R console): source("8SWEI_Snow.R")
#   Command line with argument:        Rscript 8SWEI_Snow.R 1
#   Batch / non-interactive fallback:  edit the DEFAULT_SCF_METHOD line below
##############################################

# ---- Libraries ----
library(terra)
library(ncdf4)
library(zoo)
library(writexl)
library(parallel)
library(lmomco)   # required for option 2 (SSPI-1 zero-inflated gamma)
# also used for the regional L-moment / kappa-distribution
# machinery (discordancy, heterogeneity, regional gamma fit)
# stats::ks.test (base R) is used for the gamma goodness-of-fit test.

# ---- USER SELECTION: SCF DOMAIN METHOD ----
cat("\n============================================================\n")
cat("SCF DOMAIN DEFINITION: SELECT METHOD\n")
cat("============================================================\n")
cat("  1 = Huning & AghaKouchak (2020) ??? fixed 5% SCF threshold\n")
cat("                          (standard approach from the original global paper)\n")
cat("  2 = SSPI-1 (Standardized SnowPack Index, monthly)\n")
cat("                          (zero-inflated gamma, no SCF mask)\n\n")

# ---- To use in batch/scheduled jobs, set DEFAULT_SCF_METHOD to "1" or "2"
DEFAULT_SCF_METHOD <- "1"

scf_method <- ""

# Priority 1: command-line argument (e.g. Rscript script.R 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && trimws(args[1]) %in% c("1", "2")) {
  scf_method <- trimws(args[1])
  cat(sprintf("??? Method set from command-line argument: %s\n", scf_method))
  
  # Priority 2: RStudio GUI dialog (works with Run button AND source())
} else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable()) {
  raw_input <- rstudioapi::showPrompt(
    title   = "SCF Domain Method",
    message = paste0(
      "Enter method number (1 or 2):\n\n",
      "  1 = Huning & AghaKouchak (2020)  [default]\n",
      "  2 = SSPI-1 (zero-inflated gamma)\n"
    ),
    default = DEFAULT_SCF_METHOD
  )
  if (is.null(raw_input)) {
    # User clicked Cancel ??? use the default
    scf_method <- DEFAULT_SCF_METHOD
    cat(sprintf("??? Dialog cancelled: using default Method %s\n", scf_method))
  } else {
    scf_method <- trimws(raw_input)
  }
  while (!scf_method %in% c("1", "2")) {
    raw_input <- rstudioapi::showPrompt(
      title   = "Invalid input ??? try again",
      message = "Please enter 1 or 2:",
      default = DEFAULT_SCF_METHOD
    )
    if (is.null(raw_input)) {
      scf_method <- DEFAULT_SCF_METHOD
      break
    }
    scf_method <- trimws(raw_input)
  }
  cat(sprintf("??? Method set from RStudio dialog: %s\n", scf_method))
  
  # Priority 3: plain readline() fallback (R console, non-RStudio interactive)
} else if (interactive()) {
  while (!scf_method %in% c("1", "2")) {
    scf_method <- trimws(readline(prompt = "Enter 1 or 2: "))
    if (!scf_method %in% c("1", "2")) cat("  Invalid input. Please enter 1 or 2.\n")
  }
  
  # Priority 4: non-interactive fallback default (batch / Rscript without args)
} else {
  scf_method <- DEFAULT_SCF_METHOD
  cat(sprintf("??? Non-interactive session detected: using default Method %s\n", scf_method))
  cat("  (Edit DEFAULT_SCF_METHOD at the top of the script to change this.)\n")
}

if (scf_method == "1") {
  cat("\n??? Method selected: Huning & AghaKouchak (2020) fixed 5% SCF threshold\n")
} else {
  cat("\n??? Method selected: SSPI-1 (zero-inflated gamma, monthly)\n")
}

scf_method_label <- if (scf_method == "1") "Huning & AghaKouchak (2020)" else "SSPI-1 (zero-inflated gamma, no SCF mask)"

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")
out_dir <- if (scf_method == "2") "sspi_results_monthly" else "swei_results_seasonal"
TIMESTAMP_OUTPUT_SUBDIR <- TRUE  # set TRUE to give every run its own folder
if (TIMESTAMP_OUTPUT_SUBDIR) {
  out_dir <- file.path(out_dir, format(Sys.time(), "%Y%m%d_%H%M%S"))
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# NOTE: run_settings.txt is written further below, after the RFA CONFIGURATION
# block, because it reports on REGIONALIZE_SSPI / BIAS_CORRECT / STATION_OBS_PATH,
# which aren't defined until that block runs.
# ==============================================================================
#   REGIONAL FREQUENCY ANALYSIS (RFA) CONFIGURATION  [method 2 / SSPI-1 only]
# ==============================================================================
# Set REGIONALIZE_SSPI <- FALSE to fall back to the original pixel-by-pixel
# SSPI-1 fitting (kept fully intact as a fallback path throughout the script).
REGIONALIZE_SSPI <- FALSE

# --- Regionalization (homogeneous sub-region) settings ---
# Optional DEM raster used to build elevation-band x aspect-class regions,
# analogous to the paper's A1-A6 avalanche-occurrence regions. If this file
# does not exist, the script automatically falls back to a spatial k-means
# regionalization on pixel coordinates (still pools pixels, but without a
# physical elevation/aspect basis) and prints a warning.

# NOTE: the DEM on disk is an ESRI GRID export (folder "nechako-dem"
# containing hdr.adf / w001001.adf / etc., produced via ArcGIS Extract By
# Mask) - not a GeoTIFF, and the folder is named with a hyphen, not an
# underscore. The previous path ("Spatial/nechako_dem.tif") never existed,
# which silently sent every run down the unweighted k-means spatial
# fallback. terra::rast() reads ESRI GRID directories directly when pointed
# at the folder itself (no file extension needed), so just correcting the
# path is sufficient - no change to build_regions() is required.
DEM_PATH          <- "Spatial/nechako-dem"
N_ELEV_BANDS      <- 4      # number of elevation quantile bands
USE_ASPECT        <- TRUE   # split each elevation band further by aspect
# (Restored after testing: USE_ASPECT <- FALSE roughly tripled mean H0
#  (12.9 -> 31.7) and pushed the flag rate to 94% - aspect is capturing real
#  solar-loading/melt-timing differences between slopes, not adding noise.
#  N_ELEV_BANDS 4 -> 3 made no material difference to heterogeneity, so left
#  at 4 for a somewhat finer elevation basis.)
N_ASPECT_CLASSES  <- 4      # N / E / S / W (flat slopes pooled into nearest band only)
N_REGIONS_FALLBACK <- 6     # number of regions for the k-means spatial fallback
MIN_REGION_PIXELS <- 8      # regions smaller than this are merged into their
# nearest neighbouring region before pooling
# (Hosking & Wallis 1997 recommend >= ~7 sites)

# --- Index-value normalization (Eq. 1-2 analogue) ---
INDEX_VALUE_STAT <- "median"  # "median" (robust to zero-inflation skew) or "mean"

# --- Discordancy / heterogeneity screening (Hosking & Wallis 1993, 1997) ---
DISCORDANCY_SCREENING <- TRUE
HETEROGENEITY_NSIM    <- 10   # Monte-Carlo simulations for H0 (paper uses 5000;
# reduced here for tractability at pixel-grid scale -
# H0 is stable to within ~0.05-0.1 at Nsim=10 for
# regions of this size; increase if runtime allows)
H0_REGION_FLAG_THRESHOLD <- 1  # H0 >= 1 -> "possibly heterogeneous" flag (paper convention)

# --- Goodness-of-fit screening on the pooled regional gamma fit ---
GOF_ALPHA <- 0.05   # KS-test significance level; region/months failing this
# fall back to a regional empirical (Gringorten-style) CDF

# --- Bias correction against independent station SWE observations ---
# Optional CSV with (at minimum) columns: station_id, lon, lat, date, swe_mm
# Dates should bracket the SSPI reference period. If this file is absent,
# bias correction is automatically skipped (BIAS_CORRECT is forced to FALSE)
# and a console message documents that ERA5-Land is being used uncorrected.
#
# NOTE: correction is applied basin-wide (single scaling factor per calendar
# month, pooled across all stations), NOT per-region. With only a handful of
# stations spanning a narrow elevation range, per-region quantile mapping
# (the original design) had far too little data per region/month to be
# trustworthy - see fit_basin_bias_correction() and the "BASIN-WIDE BIAS
# CORRECTION" section below.
STATION_OBS_PATH <- "Spatial/station_swe_obs.csv"
BIAS_CORRECT     <- TRUE
BIAS_CORRECT_CLIP <- c(0.5, 2.0)  # sanity clip on correction coefficients C,
# mirroring the paper's flag of |B(TR)|>0.5
# as the threshold for "correction needed"

cat(sprintf("\n??? RFA regionalization for SSPI-1: %s\n",
            if (scf_method == "2" && REGIONALIZE_SSPI) "ENABLED" else "disabled"))

writeLines(
  c(sprintf("run_timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("scf_method: %s (%s)", scf_method, scf_method_label),
    sprintf("REGIONALIZE_SSPI: %s", if (scf_method == "2") REGIONALIZE_SSPI else NA),
    sprintf("BIAS_CORRECT: %s", if (scf_method == "2") BIAS_CORRECT else NA),
    sprintf("station_file_found: %s", if (scf_method == "2") file.exists(STATION_OBS_PATH) else NA)),
  file.path(out_dir, "run_settings.txt")
)


# ---- Load Basin Boundary ----
basin_path <- "Spatial/nechakoBound_dissolve.kmz"
basin <- if (file.exists(basin_path)) {
  tmp <- tempfile(); dir.create(tmp, showWarnings = FALSE)
  utils::unzip(basin_path, exdir = tmp)
  kml <- list.files(tmp, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)[1]
  v <- vect(kml); if (nrow(v) > 1L) v <- aggregate(v)
  unlink(tmp, recursive = TRUE); v
} else NULL
if (!is.null(basin)) {
  cat("??? Basin boundary loaded\n")
} else {
  stop("Basin boundary not found: ", basin_path)
}

# ---- Input files ----
swe_file <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
scf_file <- "monthly_data_direct/snow_cover_monthly.nc"
if (!file.exists(swe_file)) stop("SWE file not found: ", swe_file)
if (scf_method == "1" && !file.exists(scf_file))
  stop("SCF file not found (required for method 1): ", scf_file)

cat("\n============================================================\n")
cat(if (scf_method == "2") "SSPI-1 CALCULATION (PARALLEL MODE)\n"
    else "SEASONAL SWEI CALCULATION (PARALLEL MODE)\n")
cat("============================================================\n\n")

cat("===== READING NETCDF FILES =====\n")
swe <- rast(swe_file)
if (scf_method == "1") {
  scf_stack <- rast(scf_file)
} else {
  scf_stack <- NULL  # not needed for SSPI-1
}

# ---- Time extraction ----
extract_time_dimension <- function(raster_obj, file_path) {
  t <- time(raster_obj)
  if (!all(is.na(t))) return(as.Date(t))
  nc <- nc_open(file_path)
  if ("valid_time" %in% names(nc$var)) {
    tv <- ncvar_get(nc, "valid_time")
    nc_close(nc)
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  if ("time" %in% names(nc$dim)) {
    tv <- ncvar_get(nc, "time")
    nc_close(nc)
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  nc_close(nc)
  n <- nlyr(raster_obj)
  start <- as.Date("1950-01-01")
  seq(start, by = "month", length.out = n)
}
dates_swe <- extract_time_dimension(swe, swe_file)
if (!is.null(scf_stack)) {
  dates_scf <- extract_time_dimension(scf_stack, scf_file)
  terra::time(scf_stack) <- dates_scf
} else {
  dates_scf <- NULL
}
terra::time(swe) <- dates_swe
dates <- dates_swe

# ---- Align Spatial Extents and Reproject ----
cat("\n===== SPATIAL ALIGNMENT =====\n")

if (!is.null(scf_stack)) {
  # 1. Fix Longitude Mismatches (0-360 vs -180-180)
  swe_ext <- ext(swe)
  scf_ext <- ext(scf_stack)
  
  if (scf_ext[2] > 180 && swe_ext[2] <= 180) {
    cat("??? Rotating SCF from 0-360 to -180-180 to match SWE...\n")
    scf_stack <- terra::rotate(scf_stack)
  } else if (swe_ext[2] > 180 && scf_ext[2] <= 180) {
    cat("??? Rotating SWE from 0-360 to -180-180 to match SCF...\n")
    swe <- terra::rotate(swe)
  }
}

# 2. Crop early to save memory and isolate the basin
if (!is.null(basin)) {
  cat("??? Cropping rasters to basin extent...\n")
  swe <- crop(swe, project(basin, crs(swe)))
  if (!is.null(scf_stack)) {
    scf_stack <- crop(scf_stack, project(basin, crs(scf_stack)))
  }
}

# 3. Resample to common grid
if (!is.null(scf_stack) && !compareGeom(scf_stack[[1]], swe[[1]], stopOnError = FALSE)) {
  cat("??? Resampling SWE to match SCF grid...\n")
  swe <- resample(swe, scf_stack, method = "bilinear")
}

target_crs <- if (!is.null(scf_stack)) crs(scf_stack) else crs(swe)
if (!is.null(basin) && !same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
}
# ---- Basin masking ----
cat("\n===== MASKING SWE TO BASIN BOUNDARY =====\n")
swe <- mask(swe, basin, inverse = FALSE, touches = TRUE)
basin_pixels <- global(swe[[1]], "notNA")$notNA
total_pixels <- ncell(swe)
cat(sprintf("??? Basin masking complete: %d pixels (%.1f%% of raster)\n",
            basin_pixels, 100 * basin_pixels / total_pixels))

# Safety Check against NaN calculation
if (is.na(basin_pixels) || basin_pixels == 0) {
  stop("FATAL ERROR: Basin masking resulted in 0 valid pixels. The SWE and SCF rasters may not overlap the basin boundary properly after resampling.")
}

# ---- Convert SWE to mm if needed ----
swe_mean_sample <- global(swe[[1]], "mean", na.rm = TRUE)$mean
if (!is.na(swe_mean_sample) && swe_mean_sample < 10) {
  swe <- swe * 1000
  cat("??? Converted SWE from meters to mm\n")
}

# ---- Convert SCF to fraction if needed ----
if (!is.null(scf_stack)) {
  scf_max <- global(scf_stack[[1]], "max", na.rm = TRUE)$max
  if (!is.na(scf_max) && scf_max > 1.5) {
    scf_stack <- scf_stack / 100
    cat("??? Converted SCF from percent to fraction\n")
  }
}

# ---- Prepare matrices ----
swe_matrix <- values(swe, mat = TRUE)
dates <- dates_swe
months_all <- as.integer(format(dates, "%m"))
years_all <- as.integer(format(dates, "%Y"))
n_pixels <- nrow(swe_matrix)
cat(sprintf("Processing %d pixels of the rectangle covering basin ...\n", n_pixels))

# ==============================================================================
#   USER CHOICE: SCF DOMAIN DEFINITION METHOD
# ==============================================================================
cat(sprintf("\n??? Proceeding with Method %s: %s\n", scf_method, scf_method_label))

cat("\n===== STEP 2: Building SCF mask =====\n")
target_k <- 3

if (scf_method == "2") {
  
  # ============================================================================
  #   BRANCH 2 ??? SSPI-1 (Zero-inflated Gamma, Monthly, No SCF Mask)
  # ============================================================================
  cat("\n??? SSPI-1 selected: skipping SCF mask build (not applicable).\n")
  cat("  Zero-inflated gamma fitting will be applied per calendar month.\n")
  scf_mask_list  <- NULL
  basin_scf_mask <- NULL
  
  # Reference period for SSPI-1
  sspi_ref_start <- min(dates) 
  sspi_ref_end   <- max(dates) 
  ref_idx        <- which(dates >= sspi_ref_start & dates <= sspi_ref_end)
  if (length(ref_idx) < 400) {
    warning(sprintf("SSPI-1 reference period has only %d months (recommended: >=480)", length(ref_idx)))
  }
  cat(sprintf("??? SSPI-1 reference period: %s to %s (%d months)\n",
              dates[ref_idx[1]], dates[ref_idx[length(ref_idx)]], length(ref_idx)))
  
  # Zero-inflation screening
  cat("\n===== ZERO-INFLATION SCREENING (SSPI-1) =====\n")
  ref_swe_sspi   <- swe_matrix[, ref_idx, drop = FALSE]
  zero_prop_sspi <- rowMeans(ref_swe_sspi <= 0, na.rm = TRUE)
  invalid_pix_sspi <- which(zero_prop_sspi > 0.50)
  valid_pix_sspi   <- which(zero_prop_sspi <= 0.50)
  cat(sprintf("??? Valid pixels (<=50%% zero SWE in ref period): %d (%.1f%%)\n",
              length(valid_pix_sspi), 100 * length(valid_pix_sspi) / n_pixels))
  if (length(invalid_pix_sspi) > 0) {
    cat(sprintf("  Excluded pixels (>50%% zero SWE): %d\n", length(invalid_pix_sspi)))
    swe_matrix[invalid_pix_sspi, ] <- NA
  }
  # ---- Diagnostic: is this screen actually doing anything, or is the
  # ~48% "valid pixel" figure just the basin boundary mask? Distinguishing
  # these matters: if the screen removes 0 interior basin pixels, the 50%
  # threshold is currently a no-op for this basin/period, not evidence the
  # domain is being aggressively filtered.
  basin_rows_prescreen <- which(is.finite(zero_prop_sspi))  # rows with real basin data
  invalid_within_basin <- intersect(invalid_pix_sspi, basin_rows_prescreen)
  
  cat("\n  --- Zero-inflation screen coverage check ---\n")
  cat(sprintf("  Rectangle pixels (grid cells, incl. outside basin): %d\n", n_pixels))
  cat(sprintf("  Basin-boundary pixels (from earlier masking step):  %d (%.1f%% of rectangle)\n",
              basin_pixels, 100 * basin_pixels / n_pixels))
  cat(sprintf("  Basin pixels with computable zero-fraction:         %d%s\n",
              length(basin_rows_prescreen),
              if (length(basin_rows_prescreen) != basin_pixels)
                sprintf("  [!] differs from basin_pixels (%d) - check masking/cropping alignment", basin_pixels)
              else "  (matches basin_pixels - masking consistent)"))
  cat(sprintf("  Of those, excluded here for >50%% zero SWE:         %d (%.1f%% of basin)\n",
              length(invalid_within_basin),
              100 * length(invalid_within_basin) / max(1, length(basin_rows_prescreen))))
  
  if (length(invalid_within_basin) == 0) {
    cat("  -> The zero-inflation screen removed ZERO interior basin pixels.\n")
    cat("     The overall valid-pixel percentage reported above is driven entirely\n")
    cat("     by the basin boundary mask, not by this screen. For this basin/period\n")
    cat("     the 50% threshold is currently a no-op - it remains in place as a\n")
    cat("     safeguard for other basins or reference periods where zero-inflation\n")
    cat("     could plausibly exceed 50% (e.g. warmer/lower-elevation domains).\n")
  } else {
    cat(sprintf("  -> The zero-inflation screen actively excluded %d interior basin\n",
                length(invalid_within_basin)))
    cat("     pixel(s) beyond what the basin boundary mask alone removed. Review\n")
    cat("     these if the count is unexpectedly high relative to basin size.\n")
  }
  cat("  ---------------------------------------------\n")
} else if (scf_method == "1") {
  
  # ============================================================================
  #   BRANCH 1 ??? HUNING & AGHAKOUCHAK (2020) fixed 5% SCF threshold
  # ============================================================================
  scf_threshold <- 0.05
  cat(sprintf("??? H&A (2020): building SCF mask with fixed %.0f%% threshold\n",
              100 * scf_threshold))
  
  # Mask SCF stack to the basin so pixel indices align exactly with SWE
  scf_basin      <- mask(scf_stack, basin, inverse = FALSE, touches = TRUE)
  scf_matrix_ha  <- values(scf_basin, mat = TRUE)
  months_scf_ha  <- as.integer(format(dates_scf, "%m"))
  
  scf_mask_list <- vector("list", 12)
  for (m in 1:12) {
    win_months        <- ((c(m - 3L, m - 2L, m - 1L) + 120L) %% 12L) + 1L
    idx_win           <- which(months_scf_ha %in% win_months)
    if (length(idx_win) == 0) {
      scf_mask_list[[m]] <- rep(FALSE, ncell(swe))
    } else {
      scf_clim           <- rowMeans(scf_matrix_ha[, idx_win, drop = FALSE], na.rm = TRUE)
      scf_mask_list[[m]] <- !is.na(scf_clim) & (scf_clim >= scf_threshold)
    }
    cat(sprintf("  Month %2d (%s): %d pixels pass SCF >= %.0f%%\n",
                m, month.abb[m],
                sum(scf_mask_list[[m]], na.rm = TRUE),
                100 * scf_threshold))
  }
  rm(scf_basin, scf_matrix_ha)
  
  cat("\n===== COMPUTING BASIN-LEVEL SCF MASK (H&A 2020) =====\n")
  basin_pixel_idx <- which(!is.na(swe_matrix[, 1]))
  basin_scf_mask  <- logical(12)
  for (m in 1:12) {
    mv_basin          <- scf_mask_list[[m]][basin_pixel_idx]
    basin_scf_mask[m] <- length(mv_basin) > 0 && mean(mv_basin, na.rm = TRUE) >= 0.5
  }
}

if (scf_method == "1") {
  cat(sprintf("??? Basin SCF mask: months included = %s\n",
              paste(month.abb[which(basin_scf_mask)], collapse = ", ")))
}

cat("\n===== COMPUTING BASIN-AVERAGED SWE =====\n")
if (scf_method == "1") {
  swe_basin_avg_raw <- colMeans(swe_matrix, na.rm = TRUE)
  cat(sprintf("??? Basin-averaged raw SWE: %d time steps\n", length(swe_basin_avg_raw)))
} else {
  swe_basin_avg_raw <- NULL
  cat("??? Basin-averaged SWE: skipped (not used by SSPI-1)\n")
}

swe_basin_avg_smoothed <- NULL
basin_avg_swei_results <- list()

# ==============================================================================
#   Helper Functions
# ==============================================================================
clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

# ==============================================================================
#   REGIONAL FREQUENCY ANALYSIS (RFA) HELPER FUNCTIONS  [SSPI-1 / method 2]
# ==============================================================================
# Implements the index-flood regionalization workflow described in
# Fontana et al. (2026), adapted from snow-depth extremes (block maxima, GEV)
# to monthly SWE values (full distribution, zero-inflated gamma). The unit of
# pooling here is "pixel-month within a homogeneous region" rather than
# "station-year within a homogeneous region", but the L-moment machinery
# (index value normalization, discordancy Di, heterogeneity H0, regional
# fitting, goodness-of-fit) is otherwise a direct analogue.

## ---- 1. Regionalization: elevation x aspect bands (DEM) or spatial fallback ----
build_regions <- function(swe_template, basin, dem_path,
                          n_elev_bands = 4, use_aspect = TRUE,
                          n_aspect_classes = 4, n_regions_fallback = 6,
                          min_region_pixels = 8) {
  
  n_pix <- ncell(swe_template)
  coords <- xyFromCell(swe_template, 1:n_pix)
  valid_mask <- !is.na(values(swe_template))
  
  region_id  <- rep(NA_integer_, n_pix)
  elevation  <- rep(NA_real_, n_pix)
  aspect_deg <- rep(NA_real_, n_pix)
  method_used <- "none"
  
  if (!is.null(dem_path) && file.exists(dem_path)) {
    cat(sprintf("??? Building regions from DEM: %s\n", dem_path))
    dem <- rast(dem_path)
    if (!is.null(basin)) dem <- crop(dem, project(basin, crs(dem)))
    if (!compareGeom(dem, swe_template, stopOnError = FALSE)) {
      dem <- resample(dem, swe_template, method = "bilinear")
    }
    elev_vals <- values(dem)[, 1]
    elevation[valid_mask] <- elev_vals[valid_mask]
    
    # Elevation quantile bands (so each band gets a comparable sample size)
    elev_breaks <- quantile(elevation[valid_mask], probs = seq(0, 1, length.out = n_elev_bands + 1),
                            na.rm = TRUE, type = 7)
    elev_breaks[1] <- -Inf; elev_breaks[length(elev_breaks)] <- Inf
    elev_band <- cut(elevation, breaks = unique(elev_breaks), labels = FALSE, include.lowest = TRUE)
    
    if (use_aspect) {
      asp <- terrain(dem, v = "aspect", unit = "degrees")
      asp_vals <- values(asp)[, 1]
      aspect_deg[valid_mask] <- asp_vals[valid_mask]
      # 4-class compass binning (N/E/S/W), flat (-1 aspect from terra) -> class 0
      aspect_class <- rep(NA_integer_, n_pix)
      aspect_class[valid_mask] <- {
        a <- aspect_deg[valid_mask]
        cls <- floor(((a + 45) %% 360) / 90) + 1   # 1=N,2=E,3=S,4=W
        cls[is.na(a) | a < 0] <- 0                  # flat / undefined
        as.integer(cls)
      }
      region_id <- as.integer(elev_band) * 10L + aspect_class
    } else {
      region_id <- as.integer(elev_band)
    }
    method_used <- "dem_elevation_aspect"
    
  } else {
    cat(sprintf("??? DEM not found at '%s' - falling back to spatial k-means regionalization.\n", dem_path))
    cat("  (Regions will not have a physical elevation/aspect basis. Provide DEM_PATH\n")
    cat("   for an elevation-band x aspect regionalization analogous to the paper's A1-A6.)\n")
    valid_idx <- which(valid_mask)
    if (length(valid_idx) >= n_regions_fallback) {
      set.seed(42)
      km <- kmeans(coords[valid_idx, , drop = FALSE], centers = n_regions_fallback, nstart = 10)
      region_id[valid_idx] <- km$cluster
    } else {
      region_id[valid_idx] <- 1L
    }
    method_used <- "kmeans_spatial_fallback"
  }
  
  # ---- Merge undersized regions into their nearest (by centroid) neighbour ----
  region_id[!valid_mask] <- NA_integer_
  tab <- table(region_id[valid_mask])
  small_regions <- as.integer(names(tab)[tab < min_region_pixels])
  if (length(small_regions) > 0 && length(tab) > 1) {
    centroids <- sapply(as.integer(names(tab)), function(r) {
      idx <- which(region_id == r)
      colMeans(coords[idx, , drop = FALSE])
    })
    centroids <- t(centroids)
    rownames(centroids) <- names(tab)
    for (r in small_regions) {
      others <- setdiff(as.integer(names(tab)), small_regions)
      if (length(others) == 0) next
      d <- sqrt(rowSums(sweep(centroids[as.character(others), , drop = FALSE], 2,
                              centroids[as.character(r), ], "-")^2))
      nearest <- others[which.min(d)]
      region_id[region_id == r & !is.na(region_id)] <- nearest
    }
    cat(sprintf("??? Merged %d undersized region(s) (< %d pixels) into nearest neighbour.\n",
                length(small_regions), min_region_pixels))
  }
  
  region_id <- match(region_id, sort(unique(region_id[valid_mask])))  # renumber 1..K, NA stays NA
  n_regions <- length(unique(region_id[valid_mask]))
  cat(sprintf("??? Regionalization complete (%s): %d homogeneous sub-regions, %d valid pixels.\n",
              method_used, n_regions, sum(valid_mask)))
  
  list(region_id = region_id, elevation = elevation, aspect_deg = aspect_deg,
       n_regions = n_regions, method = method_used)
}

## ---- 2. Index value (Eq. 1-2 analogue): per-pixel, per-calendar-month scale ----
compute_index_values <- function(swe_matrix, dates_vec, ref_idx, stat = "median") {
  mon <- as.integer(format(dates_vec, "%m"))
  n_pix <- nrow(swe_matrix)
  idx_val <- matrix(NA_real_, nrow = n_pix, ncol = 12)
  fn <- if (stat == "mean") mean else median
  for (m in 1:12) {
    cols <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
    if (length(cols) == 0) next
    sub <- swe_matrix[, cols, drop = FALSE]
    # Index value computed from the FULL reference-period record for that
    # calendar month (mirrors mu_i = mean of yearly maxima in Eq. 1; here the
    # "extreme" is replaced by the central tendency of the monthly record).
    idx_val[, m] <- apply(sub, 1, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 3) return(NA_real_)
      fn(x)
    })
  }
  idx_val[idx_val <= 0 | !is.finite(idx_val)] <- NA_real_
  idx_val
}

## ---- 3. L-moment ratios per pixel-month (for discordancy / heterogeneity) ----
lmom_ratios_site <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 5) return(c(t2 = NA_real_, t3 = NA_real_, t4 = NA_real_, n = length(x)))
  lm <- try(lmomco::lmoms(x), silent = TRUE)
  if (inherits(lm, "try-error") || is.null(lm) || any(!is.finite(lm$ratios[2:4]))) {
    return(c(t2 = NA_real_, t3 = NA_real_, t4 = NA_real_, n = length(x)))
  }
  c(t2 = lm$ratios[2], t3 = lm$ratios[3], t4 = lm$ratios[4], n = length(x))
}

## ---- 4. Discordancy measure Di (Hosking & Wallis 1993, Eq. 3-5) ----
# Standard critical values of Di (Hosking & Wallis 1997, Table 3)
.Dcr_table <- c(`5` = 1.333, `6` = 1.648, `7` = 1.917, `8` = 2.140, `9` = 2.329,
                `10` = 2.491, `11` = 2.632, `12` = 2.757, `13` = 2.869, `14` = 2.971,
                `15` = 3.000)
get_Dcr <- function(N) {
  if (N < 5) return(Inf)   # too few sites for a meaningful discordancy test - never flag
  key <- as.character(min(N, 15))
  unname(.Dcr_table[key])
}

discordancy_test <- function(lmom_mat) {
  # lmom_mat: matrix with columns t2, t3, t4 (one row per site/pixel)
  ok <- which(rowSums(is.finite(lmom_mat[, c("t2", "t3", "t4"), drop = FALSE])) == 3)
  N <- length(ok)
  Di <- rep(NA_real_, nrow(lmom_mat))
  if (N < 4) {
    # S becomes singular for N=3 and Di is non-informative for N<=4
    # (exactly as noted in the paper); skip the test, keep all sites.
    return(list(Di = Di, Dcr = get_Dcr(N), discordant = rep(FALSE, nrow(lmom_mat)), N = N))
  }
  U <- lmom_mat[ok, c("t2", "t3", "t4"), drop = FALSE]
  ubar <- colMeans(U)
  # Hosking & Wallis (1993) Eq. (5): S = (N-1)^-1 * sum_i (u_i-ubar)(u_i-ubar)^T
  # which is exactly the sample covariance matrix.
  S <- stats::cov(U)
  Sinv <- try(solve(S), silent = TRUE)
  if (inherits(Sinv, "try-error")) {
    return(list(Di = Di, Dcr = get_Dcr(N), discordant = rep(FALSE, nrow(lmom_mat)), N = N))
  }
  Dcr <- get_Dcr(N)
  for (k in seq_along(ok)) {
    d <- as.numeric(U[k, ] - ubar)
    Di[ok[k]] <- (1 / 3) * as.numeric(t(d) %*% Sinv %*% d)
  }
  list(Di = Di, Dcr = Dcr, discordant = !is.na(Di) & Di >= Dcr, N = N)
}

## ---- 5. Heterogeneity measure H0 (Hosking & Wallis 1993/1997, Eq. 6) ----
# Monte-Carlo simulation from a regional 4-parameter kappa distribution fitted
# to the (discordancy-screened) pooled normalized sample. Falls back to a
# gamma-based simulation if the kappa fit fails (common with heavily
# zero-inflated / short pooled samples).
heterogeneity_test <- function(site_samples, Nsim = 500, seed = 123) {
  # site_samples: list of numeric vectors (normalized, nonzero, finite values),
  # one per site/pixel within the candidate region for one calendar month.
  site_samples <- site_samples[sapply(site_samples, length) >= 5]
  Nsites <- length(site_samples)
  if (Nsites < 2) return(list(H0 = NA_real_, V0 = NA_real_, status = "too_few_sites"))
  
  pooled <- unlist(site_samples)
  reg_lm <- try(lmomco::lmoms(pooled), silent = TRUE)
  if (inherits(reg_lm, "try-error")) return(list(H0 = NA_real_, V0 = NA_real_, status = "lmom_failed"))
  
  par_obj <- try(lmomco::parkap(reg_lm), silent = TRUE)
  sim_fn <- NULL
  if (!inherits(par_obj, "try-error") && !is.null(par_obj) && isTRUE(par_obj$ifail == 0)) {
    sim_fn <- function(n) lmomco::quakap(stats::runif(n), par_obj)
    dist_used <- "kappa"
  } else {
    gam_par <- try(lmomco::pargam(reg_lm), silent = TRUE)
    if (inherits(gam_par, "try-error") || is.null(gam_par)) {
      return(list(H0 = NA_real_, V0 = NA_real_, status = "fit_failed"))
    }
    sim_fn <- function(n) lmomco::quagam(stats::runif(n), gam_par)
    dist_used <- "gamma"
  }
  
  site_n <- sapply(site_samples, length)
  lcv_at_site <- function(samples_list) {
    sapply(samples_list, function(x) {
      lm <- try(lmomco::lmoms(x), silent = TRUE)
      if (inherits(lm, "try-error") || is.null(lm)) return(NA_real_)
      lm$ratios[2]
    })
  }
  V0 <- function(lcv, n) {
    w <- n / sum(n)
    sqrt(sum(w * (lcv - sum(w * lcv, na.rm = TRUE))^2, na.rm = TRUE))
  }
  
  lcv_obs <- lcv_at_site(site_samples)
  V0_obs <- V0(lcv_obs, site_n)
  
  set.seed(seed)
  V0_sim <- numeric(Nsim)
  for (s in seq_len(Nsim)) {
    sim_samples <- lapply(site_n, function(n) pmax(sim_fn(n), 1e-6))
    lcv_sim <- lcv_at_site(sim_samples)
    V0_sim[s] <- V0(lcv_sim, site_n)
  }
  V0_sim <- V0_sim[is.finite(V0_sim)]
  if (length(V0_sim) < 10) return(list(H0 = NA_real_, V0 = V0_obs, status = "sim_failed"))
  
  H0 <- (V0_obs - mean(V0_sim)) / stats::sd(V0_sim)
  list(H0 = H0, V0 = V0_obs, status = paste0("ok_", dist_used))
}

## ---- 6. Pooled regional gamma fit + goodness-of-fit (KS test) ----
fit_regional_gamma <- function(pooled_sample) {
  pooled_sample <- pooled_sample[is.finite(pooled_sample) & pooled_sample > 0]
  if (length(pooled_sample) < 15) return(list(par = NULL, status = "insufficient_data"))
  lm <- try(lmomco::lmoms(pooled_sample), silent = TRUE)
  if (inherits(lm, "try-error")) return(list(par = NULL, status = "lmom_failed"))
  par_obj <- try(lmomco::pargam(lm), silent = TRUE)
  if (inherits(par_obj, "try-error") || is.null(par_obj)) return(list(par = NULL, status = "fit_failed"))
  list(par = par_obj, status = "ok")
}

gof_test_gamma <- function(pooled_sample, par_obj, alpha = 0.05) {
  pooled_sample <- pooled_sample[is.finite(pooled_sample) & pooled_sample > 0]
  if (is.null(par_obj) || length(pooled_sample) < 15) {
    return(list(D = NA_real_, p_value = NA_real_, pass = FALSE))
  }
  cdf_fun <- function(x) as.numeric(lmomco::cdfgam(x, par_obj))
  ks <- try(stats::ks.test(pooled_sample, cdf_fun), silent = TRUE)
  if (inherits(ks, "try-error")) return(list(D = NA_real_, p_value = NA_real_, pass = FALSE))
  list(D = unname(ks$statistic), p_value = ks$p.value, pass = ks$p.value >= alpha)
}

## ---- 7. Bias correction vs. independent station SWE obs (Eq. 9-10 analogue) ----
# Simple parametric scaling: x*_corr = C * x*_E5L, with C estimated as the
# ratio of station-to-ERA5-Land dimensionless quantiles, pooled across a
# probability grid (more stable than a single quantile). One C per
# region x calendar-month, applied to raw pixel SWE before normalization.
#
# NOTE: kept for reference / future use with a denser station network, but
# NOT called anywhere in this script anymore - see fit_basin_bias_correction()
# below and Section F of the RFA construction block for why it was replaced
# by a single basin-wide monthly correction.
fit_bias_correction <- function(station_obs, swe, region_id, swe_matrix, dates_vec,
                                ref_idx, index_values, clip_range = c(0.5, 2.0)) {
  mon <- as.integer(format(dates_vec, "%m"))
  n_regions <- max(region_id, na.rm = TRUE)
  C <- matrix(1, nrow = n_regions, ncol = 12)
  details <- data.frame()
  
  if (is.null(station_obs) || nrow(station_obs) == 0) {
    cat("??? No station observations supplied: bias correction coefficients set to 1 (no correction).\n")
    return(list(C = C, details = details))
  }
  
  station_obs$date <- as.Date(station_obs$date)
  pts <- vect(station_obs, geom = c("lon", "lat"), crs = crs(swe))
  cellnum <- cells(swe[[1]], pts)[, "cell"]
  station_obs$cell <- cellnum
  station_obs$region <- region_id[cellnum]
  station_obs <- station_obs[is.finite(station_obs$region), ]
  
  if (nrow(station_obs) == 0) {
    cat("??? Station observations did not intersect any valid basin pixel: skipping bias correction.\n")
    return(list(C = C, details = details))
  }
  
  probs <- seq(0.1, 0.9, by = 0.1)
  for (r in sort(unique(station_obs$region))) {
    for (m in 1:12) {
      sub <- station_obs[station_obs$region == r &
                           as.integer(format(station_obs$date, "%m")) == m &
                           station_obs$date >= dates_vec[ref_idx[1]] &
                           station_obs$date <= dates_vec[ref_idx[length(ref_idx)]], ]
      if (nrow(sub) < 10) next  # not enough station-months to estimate a quantile ratio
      
      # station dimensionless values: normalize each station by ITS OWN
      # reference-period median/mean for that month (paper's index-value logic
      # applied to the station record itself).
      station_idx_val <- stats::median(sub$swe_mm, na.rm = TRUE)
      if (!is.finite(station_idx_val) || station_idx_val <= 0) next
      station_star <- sub$swe_mm / station_idx_val
      q_station <- stats::quantile(station_star, probs, na.rm = TRUE, type = 7)
      
      # ERA5-Land dimensionless values pooled over all pixels in the region for
      # this month, reference period only.
      pix_in_region <- which(region_id == r)
      cols <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
      if (length(pix_in_region) == 0 || length(cols) == 0) next
      era_vals <- swe_matrix[pix_in_region, cols, drop = FALSE]
      era_idx  <- index_values[pix_in_region, m]
      era_star <- sweep(era_vals, 1, era_idx, "/")
      era_star <- era_star[is.finite(era_star) & era_star > 0]
      if (length(era_star) < 30) next
      q_era <- stats::quantile(era_star, probs, na.rm = TRUE, type = 7)
      
      ratio <- q_station / q_era
      ratio <- ratio[is.finite(ratio) & ratio > 0]
      if (length(ratio) < 3) next
      c_val <- stats::median(ratio)
      c_val <- max(min(c_val, clip_range[2]), clip_range[1])
      C[r, m] <- c_val
      details <- rbind(details, data.frame(region = r, month = m, n_station = nrow(sub),
                                           n_era5land = length(era_star), C = c_val))
    }
  }
  cat(sprintf("??? Bias correction: %d region/month combinations corrected from station data.\n",
              nrow(details)))
  list(C = C, details = details)
}

## ---- 7b. Basin-wide monthly bias correction against station SWE obs ----
prep_station_file <- function(path, station_id, lon, lat, good_grades = c(0, 1)) {
  raw <- read.csv(path, stringsAsFactors = FALSE)
  raw$date <- as.Date(raw$`Start.of.Interval..UTC.`, format = "%m/%d/%Y")  # confirm actual format
  raw <- raw[raw$Grade.Code %in% good_grades & is.finite(raw$SW.Daily.Average..mm.), ]
  
  raw$year_month <- format(raw$date, "%Y-%m")
  monthly <- aggregate(SW.Daily.Average..mm. ~ year_month, data = raw, FUN = mean)
  monthly$date <- as.Date(paste0(monthly$year_month, "-01"))
  
  data.frame(station_id = station_id, lon = lon, lat = lat,
             date = monthly$date, swe_mm = monthly$SW.Daily.Average..mm.)
}
# Single multiplicative scaling factor per calendar month (<=12 numbers
# total), pooled across ALL stations and applied uniformly to every basin
# pixel for that month. For each station, uses the nearest ERA5-Land pixel;
# for each calendar month, pools the ratio (station SWE mm) / (ERA5-Land SWE
# mm at that pixel) across all stations and overlapping years, and takes the
# median as the correction factor. Deliberately NOT quantile mapping (too few
# points per probability bin) and NOT spatially interpolated (a handful of
# stations can't support kriging/IDW without inventing spatial structure the
# data can't back up) - this is the standard low-variance "delta method"
# correction appropriate for a scarcely-observed basin.
fit_basin_bias_correction <- function(station_obs, swe_template, swe_matrix,
                                      dates_vec, clip_range = c(0.5, 2.0),
                                      min_ratios_per_month = 5) {
  factor_by_month <- rep(1, 12)
  details <- data.frame()
  
  if (is.null(station_obs) || nrow(station_obs) == 0) {
    cat("??? No station observations supplied: basin-wide bias correction set to 1 (no correction).\n")
    return(list(factor = factor_by_month, details = details))
  }
  
  station_obs$date <- as.Date(station_obs$date)
  pts <- vect(station_obs, geom = c("lon", "lat"), crs = crs(swe_template))
  
  # Nearest cell to each station (works whether or not the station falls
  # exactly inside a basin-masked NA-free pixel).
  cellnum <- cells(swe_template, pts)[, "cell"]
  station_obs$cell <- cellnum
  station_obs <- station_obs[is.finite(station_obs$cell), ]
  if (nrow(station_obs) == 0) {
    cat("??? Stations did not resolve to any raster cell: skipping basin-wide bias correction.\n")
    return(list(factor = factor_by_month, details = details))
  }
  
  mon <- as.integer(format(dates_vec, "%m"))
  yr  <- as.integer(format(dates_vec, "%Y"))
  station_obs$s_mon <- as.integer(format(station_obs$date, "%m"))
  station_obs$s_yr  <- as.integer(format(station_obs$date, "%Y"))
  
  for (m in 1:12) {
    sub_m <- station_obs[station_obs$s_mon == m, ]
    if (nrow(sub_m) == 0) next
    
    ratios_m <- c()
    for (i in seq_len(nrow(sub_m))) {
      cix <- sub_m$cell[i]
      if (is.na(cix) || cix < 1 || cix > nrow(swe_matrix)) next
      col <- which(mon == m & yr == sub_m$s_yr[i])
      if (length(col) != 1) next               # need exactly one matching month
      era_val <- swe_matrix[cix, col]
      sta_val <- sub_m$swe_mm[i]
      if (!is.finite(era_val) || era_val <= 0) next
      if (!is.finite(sta_val) || sta_val <= 0) next
      ratios_m <- c(ratios_m, sta_val / era_val)
    }
    ratios_m <- ratios_m[is.finite(ratios_m) & ratios_m > 0]
    
    if (length(ratios_m) < min_ratios_per_month) {
      details <- rbind(details, data.frame(month = m, n_ratios = length(ratios_m),
                                           factor = 1, applied = FALSE))
      next
    }
    c_val <- stats::median(ratios_m)
    c_val <- max(min(c_val, clip_range[2]), clip_range[1])
    factor_by_month[m] <- c_val
    details <- rbind(details, data.frame(month = m, n_ratios = length(ratios_m),
                                         factor = c_val, applied = TRUE))
  }
  
  n_applied <- sum(details$applied, na.rm = TRUE)
  cat(sprintf("??? Basin-wide bias correction: %d/12 months corrected from %d station(s) (need >=%d ratios/month).\n",
              n_applied, length(unique(station_obs$station_id)), min_ratios_per_month))
  if (n_applied > 0) {
    applied_rows <- details[details$applied, ]
    cat("  Monthly factors: ",
        paste(sprintf("%s=%.2f (n=%d)", month.abb[applied_rows$month],
                      applied_rows$factor, applied_rows$n_ratios), collapse = ", "),
        "\n", sep = "")
  }
  if (any(!details$applied)) {
    skipped <- details$month[!details$applied]
    cat(sprintf("  Skipped (insufficient overlapping station-months): %s\n",
                paste(month.abb[skipped], collapse = ", ")))
  }
  
  list(factor = factor_by_month, details = details)
}

gringorten_swei_seasonal <- function(v_smoothed, dates_vec, scf_mask_list,
                                     cell_index = NULL,
                                     basin_scf_mask = NULL,
                                     eps = 1e-6) {
  if (length(v_smoothed) == 0 || all(is.na(v_smoothed)) || is.null(v_smoothed)) {
    return(list(swei = rep(NA_real_, length(dates_vec)), method = rep(NA_integer_, length(dates_vec))))
  }
  v_clean <- v_smoothed
  v_clean[!is.finite(v_clean)] <- NA_real_
  n <- length(v_clean)
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n) 
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    all_m_idx <- which(mon == m)
    if (length(all_m_idx) == 0) next
    nonzero_frac <- sum(v_clean[all_m_idx] > 0, na.rm = TRUE) / length(all_m_idx)
    if (nonzero_frac < 0.75) next
    
    idx <- which(mon == m & is.finite(v_clean))
    if (length(idx) == 0) next
    
    if (is.null(scf_mask_list) || length(scf_mask_list) < m) next
    if (!is.null(cell_index)) {
      mask_vec <- scf_mask_list[[m]]
      if (length(mask_vec) < cell_index || !mask_vec[cell_index]) next
    } else if (!is.null(basin_scf_mask)) {
      if (length(basin_scf_mask) < m || !basin_scf_mask[m]) next
    }
    
    samp <- v_clean[idx]
    zero_idx <- which(samp == 0 & !is.na(samp))
    if (length(zero_idx) > 0) {
      nonzero_vals <- samp[samp > 0 & !is.na(samp)]
      if (length(nonzero_vals) > 0) {
        min_nonzero <- min(nonzero_vals)
        samp[zero_idx] <- runif(length(zero_idx), min = 0, max = min_nonzero)
      }
    }
    valid_idx <- which(!is.na(samp))
    if (length(valid_idx) < 3) next
    
    samp_valid <- samp[valid_idx]
    idx_valid <- idx[valid_idx]
    
    r <- rank(samp_valid, ties.method = "average")
    N <- length(samp_valid)
    p_val <- (r - 0.44) / (N + 0.12)
    p_val <- clip_prob(p_val, eps = eps)
    z_val <- qnorm(p_val)
    
    z[idx_valid] <- z_val
    method_used[idx_valid] <- 1
  }
  list(swei = z, method = method_used)
}

monthly_sspi <- function(v, dates_vec, ref_idx, eps = 1e-6) {
  n <- length(v)
  if (n == 0 || all(is.na(v)) || is.null(v)) {
    return(list(sspi = rep(NA_real_, n), method = rep(NA_integer_, n)))
  }
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n) 
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    ref_mon_idx <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
    ref_samp <- v_clean[ref_mon_idx]
    ref_samp <- ref_samp[is.finite(ref_samp)]
    obs_mon_idx <- which(mon == m)
    
    if (length(ref_samp) < 10) {
      if (length(obs_mon_idx) >= 3) {
        obs_vals <- v_clean[obs_mon_idx]
        fin_idx  <- is.finite(obs_vals)
        n_fin    <- sum(fin_idx)
        if (n_fin >= 3) {
          p <- (rank(obs_vals[fin_idx], ties.method = "average") - 0.44) / (n_fin + 0.12)
          p <- clip_prob(p, eps = eps)
          z[obs_mon_idx[fin_idx]] <- qnorm(p)
          method_used[obs_mon_idx[fin_idx]] <- 2
        }
      }
      next
    }
    ref_var <- var(ref_samp, na.rm = TRUE)
    if (!is.finite(ref_var) || ref_var < .Machine$double.eps) {
      z[obs_mon_idx] <- 0
      method_used[obs_mon_idx] <- 3
      next
    }
    if (ref_var < 0.01) {
      obs_vals <- v_clean[obs_mon_idx]
      fin_idx <- is.finite(obs_vals)
      if (sum(fin_idx) >= 3) {
        p <- (rank(obs_vals[fin_idx], ties.method = "average") - 0.44) /
          (sum(fin_idx) + 0.12)
        p <- clip_prob(p, eps = eps)
        z[obs_mon_idx[fin_idx]] <- qnorm(p)
        method_used[obs_mon_idx[fin_idx]] <- 2
      }
      next
    }
    p0 <- mean(ref_samp <= 0, na.rm = TRUE)
    ref_pos <- ref_samp[ref_samp > 0]
    if (length(ref_pos) >= 10) {
      lm_obj <- try(lmomco::lmoms(ref_pos), silent = TRUE)
      if (!inherits(lm_obj, "try-error") && !is.null(lm_obj)) {
        par_obj <- try(lmomco::pargam(lm_obj), silent = TRUE)
        if (!inherits(par_obj, "try-error") && !is.null(par_obj)) {
          obs_vals <- v_clean[obs_mon_idx]
          p_m <- rep(NA_real_, length(obs_vals))
          pos_idx <- which(is.finite(obs_vals) & obs_vals > 0)
          if (length(pos_idx) > 0) {
            Fg <- try(lmomco::cdfgam(obs_vals[pos_idx], par_obj), silent = TRUE)
            if (!inherits(Fg, "try-error")) p_m[pos_idx] <- p0 + (1 - p0) * as.numeric(Fg)
          }
          zero_idx2 <- which(is.finite(obs_vals) & obs_vals <= 0)
          if (length(zero_idx2) > 0) p_m[zero_idx2] <- p0
          p_m <- clip_prob(p_m, eps = eps)
          z[obs_mon_idx] <- qnorm(p_m)
          method_used[obs_mon_idx] <- 1
          next
        }
      }
    }
    all_samp <- c(ref_samp, v_clean[obs_mon_idx])
    all_samp <- all_samp[is.finite(all_samp)]
    obs_vals <- v_clean[obs_mon_idx]
    fin_obs <- is.finite(obs_vals)
    if (length(all_samp) >= 10 && sum(fin_obs) > 0) {
      combined <- c(all_samp, obs_vals[fin_obs])
      r_all <- rank(combined, ties.method = "average")
      r_obs <- r_all[(length(all_samp) + 1):(length(all_samp) + sum(fin_obs))]
      p <- (r_obs - 0.44) / (length(combined) + 0.12)
      p <- clip_prob(p, eps = eps)
      z[obs_mon_idx[fin_obs]] <- qnorm(p)
      method_used[obs_mon_idx[fin_obs]] <- 2
    }
  }
  fin_z <- is.finite(z)
  z[fin_z & z < -5] <- -5
  z[fin_z & z >  5] <-  5
  list(sspi = z, method = method_used)
}

## ---- 8. Regionalized SSPI-1 (uses pooled regional gamma / empirical fit) ----
# Method codes: 1=regional gamma (GoF pass), 2=regional pooled empirical
# (GoF fail, still regional), 3=local gamma fallback (no usable region),
# 4=local empirical fallback, 5=local zero-variance fallback.
# Falls back lazily to the original monthly_sspi() (per-pixel) whenever the
# regional fit is unusable for a given pixel/month (no region assigned,
# missing/zero index value, discordant pixel, or insufficient pooled sample).
monthly_sspi_regional <- function(v, dates_vec, ref_idx, pixel_idx,
                                  region_vec, index_values, C, region_par,
                                  discordant_mat, eps = 1e-6) {
  n <- length(v)
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n)
  mon <- as.integer(format(dates_vec, "%m"))
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  
  region_here <- if (pixel_idx <= length(region_vec)) region_vec[pixel_idx] else NA_integer_
  old_result <- NULL  # lazily computed pixel-local fallback (only if actually needed)
  
  for (m in 1:12) {
    obs_mon_idx <- which(mon == m)
    if (length(obs_mon_idx) == 0) next
    
    mu <- if (!is.na(region_here)) index_values[pixel_idx, m] else NA_real_
    disc_flag <- if (!is.null(discordant_mat)) isTRUE(discordant_mat[pixel_idx, m]) else FALSE
    par_info <- if (!is.na(region_here) && !is.null(region_par[[region_here]])) {
      region_par[[region_here]][[m]]
    } else NULL
    # pool_ok now additionally requires gof_pass == TRUE. Region/months that
    # failed the KS goodness-of-fit, OR were excluded from pooling upstream
    # for failing the H0 heterogeneity screen (region_par[[r]][[m]]$status ==
    # "heterogeneous_skip", set in the RFA construction block), have
    # par = NULL / gof_pass = FALSE and so never qualify here. There is
    # deliberately no regional-empirical fallback: a region/month either uses
    # the pooled regional gamma, or falls straight through to per-pixel local
    # fitting below (monthly_sspi) - see prior discussion on why a pooled
    # empirical CDF was dropped as an intermediate option.
    pool_ok <- !is.null(par_info) && isTRUE(par_info$gof_pass) && !is.null(par_info$par) &&
      !is.null(par_info$pooled_sample) && length(par_info$pooled_sample) >= 15
    use_regional <- !is.na(region_here) && is.finite(mu) && mu > 0 && !disc_flag && pool_ok
    
    if (use_regional) {
      c_val <- C[region_here, m]
      mu_corr <- mu * c_val
      x_star <- (v_clean[obs_mon_idx] * c_val) / mu_corr
      
      # zero-probability stays pixel-specific (paper keeps index values site
      # specific even though the dimensionless shape is shared regionally)
      ref_mon_idx <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
      ref_samp <- v_clean[ref_mon_idx]
      ref_samp <- ref_samp[is.finite(ref_samp)]
      p0 <- if (length(ref_samp) >= 5) mean(ref_samp <= 0, na.rm = TRUE) else mean(v_clean[obs_mon_idx] <= 0, na.rm = TRUE)
      if (!is.finite(p0)) p0 <- 0
      
      p_m <- rep(NA_real_, length(obs_mon_idx))
      pos_idx <- which(is.finite(x_star) & x_star > 0)
      
      # ---- Pooled regional gamma (passed KS goodness-of-fit AND H0 screen) ----
      if (length(pos_idx) > 0) {
        Fg <- try(lmomco::cdfgam(x_star[pos_idx], par_info$par), silent = TRUE)
        if (!inherits(Fg, "try-error")) p_m[pos_idx] <- p0 + (1 - p0) * as.numeric(Fg)
      }
      method_code <- 1L
      zero_idx2 <- which(is.finite(x_star) & x_star <= 0)
      if (length(zero_idx2) > 0) p_m[zero_idx2] <- p0
      p_m <- clip_prob(p_m, eps = eps)
      z[obs_mon_idx] <- qnorm(p_m)
      method_used[obs_mon_idx] <- method_code
      next
    }
    
    # ---- Fallback: no usable regional fit for this pixel/month ----
    if (is.null(old_result)) {
      old_result <- monthly_sspi(v, dates_vec, ref_idx, eps = eps)
    }
    z[obs_mon_idx] <- old_result$sspi[obs_mon_idx]
    old_codes <- old_result$method[obs_mon_idx]   # 1=gamma,2=empirical,3=zero-var (local)
    method_used[obs_mon_idx] <- ifelse(is.na(old_codes), NA_integer_, old_codes + 2L)
  }
  
  fin_z <- is.finite(z)
  z[fin_z & z < -5] <- -5
  z[fin_z & z >  5] <-  5
  list(sspi = z, method = method_used)
}

# ==============================================================================
#   BASIN-WIDE BIAS CORRECTION AGAINST STATION SWE OBSERVATIONS  [method 2 only]
# ==============================================================================
# Runs ONCE, here - after fit_basin_bias_correction() has been defined above,
# and BEFORE the RFA construction block / main calculation loop use swe_matrix
# for anything. Applying it at this point (rather than deep inside the
# REGIONALIZE_SSPI==TRUE branch) means BOTH the regionalized path
# (monthly_sspi_regional, via build_regions/compute_index_values/the pooled
# gamma fits) and the pixel-by-pixel path (monthly_sspi) see the same
# corrected swe_matrix, regardless of REGIONALIZE_SSPI.
#
# This is the ONLY place station correction is applied - see Section F inside
# the RFA construction block below, which intentionally no-ops the old
# per-region fit_bias_correction() to avoid double-correcting.
if (scf_method == "2") {
  cat("\n===== BASIN-WIDE BIAS CORRECTION vs STATION OBSERVATIONS =====\n")
  
  if (BIAS_CORRECT && file.exists(STATION_OBS_PATH)) {
    station_obs_raw <- tryCatch(read.csv(STATION_OBS_PATH, stringsAsFactors = FALSE),
                                error = function(e) NULL)
    bbc <- fit_basin_bias_correction(station_obs_raw, swe[[1]], swe_matrix, dates,
                                     clip_range = BIAS_CORRECT_CLIP)
    basin_bias_factor <- bbc$factor
    basin_bias_log     <- bbc$details
    
    months_vec_bc <- as.integer(format(dates, "%m"))
    for (m in 1:12) {
      if (basin_bias_factor[m] != 1) {
        cols_m <- which(months_vec_bc == m)
        swe_matrix[, cols_m] <- swe_matrix[, cols_m] * basin_bias_factor[m]
      }
    }
    
    bc_diag_dir <- file.path(out_dir, "rfa_diagnostics")
    if (!dir.exists(bc_diag_dir)) dir.create(bc_diag_dir, recursive = TRUE)
    if (nrow(basin_bias_log) > 0) {
      write.csv(basin_bias_log, file.path(bc_diag_dir, "basin_wide_bias_correction_by_month.csv"),
                row.names = FALSE)
    }
  } else {
    basin_bias_factor <- rep(1, 12)
    basin_bias_log <- data.frame()
    if (isTRUE(BIAS_CORRECT)) {
      cat(sprintf("??? BIAS_CORRECT = TRUE but no station file found at '%s'.\n", STATION_OBS_PATH))
      cat("  Proceeding WITHOUT bias correction (factor = 1 for every month). ERA5-Land\n")
      cat("  snow variables carry known systematic biases (Kouki et al. 2023; Blau et al.\n")
      cat("  2024) - supply station SWE observations at STATION_OBS_PATH (columns:\n")
      cat("  station_id, lon, lat, date, swe_mm) to enable this correction.\n")
    } else {
      cat("??? BIAS_CORRECT = FALSE: skipping basin-wide bias correction.\n")
    }
  }
}

# ==============================================================================
#   REGIONAL FREQUENCY ANALYSIS (RFA) CONSTRUCTION  [SSPI-1 / method 2 only]
# ==============================================================================
# Builds everything monthly_sspi_regional() needs: homogeneous regions, per-
# pixel index values, discordancy screening, heterogeneity diagnostics, and
# pooled regional gamma fits with goodness-of-fit screening. All of this runs
# ONCE (not per-pixel-in-parallel) and the resulting small objects are
# exported to the cluster workers. Station bias correction has already been
# applied basin-wide, upstream, directly to swe_matrix (see above) - so this
# block operates on already-corrected SWE.
if (scf_method == "2" && REGIONALIZE_SSPI) {
  
  cat("\n============================================================\n")
  cat("REGIONAL FREQUENCY ANALYSIS (RFA) CONSTRUCTION FOR SSPI-1\n")
  cat("============================================================\n")
  
  ## ---- A. Regionalization (elevation x aspect, or spatial k-means fallback) ----
  regions_built <- build_regions(swe[[1]], basin, DEM_PATH,
                                 n_elev_bands = N_ELEV_BANDS, use_aspect = USE_ASPECT,
                                 n_aspect_classes = N_ASPECT_CLASSES,
                                 n_regions_fallback = N_REGIONS_FALLBACK,
                                 min_region_pixels = MIN_REGION_PIXELS)
  region_id <- regions_built$region_id
  region_id[invalid_pix_sspi] <- NA_integer_   # zero-inflation exclusions stay excluded
  n_regions <- length(unique(region_id[!is.na(region_id)]))
  
  ## ---- B. Per-pixel, per-calendar-month index values (Eq. 1-2 analogue) ----
  cat("\n===== COMPUTING PER-PIXEL INDEX VALUES =====\n")
  index_values <- compute_index_values(swe_matrix, dates, ref_idx, stat = INDEX_VALUE_STAT)
  cat(sprintf("??? Index values (%s of reference-period monthly SWE): %.1f%% pixel-months valid\n",
              INDEX_VALUE_STAT, 100 * mean(is.finite(index_values))))
  
  ## ---- C. L-moment ratios per pixel-month, on normalized reference samples ----
  cat("\n===== COMPUTING L-MOMENT RATIOS PER PIXEL-MONTH =====\n")
  mon_ref <- as.integer(format(dates[ref_idx], "%m"))
  lmom_array <- array(NA_real_, dim = c(n_pixels, 12, 3),
                      dimnames = list(NULL, NULL, c("t2", "t3", "t4")))
  for (m in 1:12) {
    cols <- ref_idx[mon_ref == m]
    if (length(cols) == 0) next
    sub  <- swe_matrix[, cols, drop = FALSE]
    idxv <- index_values[, m]
    elig <- which(!is.na(region_id) & is.finite(idxv) & idxv > 0)
    for (i in elig) {
      x_star <- sub[i, ]
      x_star <- x_star[is.finite(x_star)] / idxv[i]
      lr <- lmom_ratios_site(x_star)
      lmom_array[i, m, ] <- lr[c("t2", "t3", "t4")]
    }
  }
  cat("??? L-moment ratios computed for all eligible pixel-months.\n")
  
  ## ---- D. Discordancy screening (Hosking & Wallis Di, Eq. 3-5) ----
  cat("\n===== DISCORDANCY SCREENING (Hosking & Wallis D_i) =====\n")
  discordant_mat <- matrix(FALSE, n_pixels, 12)
  discordancy_log <- data.frame()
  if (DISCORDANCY_SCREENING) {
    for (r in seq_len(n_regions)) {
      pix_r <- which(region_id == r)
      if (length(pix_r) < 4) next  # Di non-informative for N<=4 (S singular at N=3)
      for (m in 1:12) {
        lm_mat <- cbind(t2 = lmom_array[pix_r, m, "t2"],
                        t3 = lmom_array[pix_r, m, "t3"],
                        t4 = lmom_array[pix_r, m, "t4"])
        dt <- discordancy_test(lm_mat)
        if (dt$N >= 4 && any(dt$discordant)) {
          flagged <- pix_r[which(dt$discordant)]
          discordant_mat[flagged, m] <- TRUE
          discordancy_log <- rbind(discordancy_log,
                                   data.frame(region = r, month = m, pixel = flagged,
                                              Di = dt$Di[dt$discordant], Dcr = dt$Dcr))
        }
      }
    }
  }
  cat(sprintf("??? Discordancy screening: %d pixel-months flagged discordant (excluded from pooling).\n",
              sum(discordant_mat)))
  
  ## ---- E. Heterogeneity test (Hosking & Wallis H0, Monte Carlo) ----
  cat("\n===== HETEROGENEITY TEST (Hosking & Wallis H0) =====\n")
  
  # Each region/month H0 test is completely independent of every other one
  # (own pooled sample, own kappa/gamma fit, own Nsim simulations), so build
  # the task list up front and farm it out across cores rather than looping
  # serially. This is the dominant cost in the RFA construction block.
  het_tasks <- list()
  for (r in seq_len(n_regions)) {
    pix_r <- which(region_id == r)
    if (length(pix_r) < 2) next
    for (m in 1:12) {
      keep <- pix_r[!discordant_mat[pix_r, m]]
      if (length(keep) < 2) next
      het_tasks[[length(het_tasks) + 1]] <- list(region = r, month = m, keep = keep)
    }
  }
  cat(sprintf("??? %d region-month heterogeneity tests queued (Nsim=%d each)\n",
              length(het_tasks), HETEROGENEITY_NSIM))
  
  if (length(het_tasks) > 0) {
    n_cores_het <- max(1L, detectCores() - 1L)
    cat(sprintf("??? Running heterogeneity tests in parallel (%d cores)\n", n_cores_het))
    cl_het <- makeCluster(n_cores_het)
    clusterSetRNGStream(cl_het, iseed = 40)
    clusterExport(cl_het, varlist = c("heterogeneity_test", "swe_matrix", "index_values",
                                      "ref_idx", "mon_ref", "HETEROGENEITY_NSIM"),
                  envir = environment())
    clusterEvalQ(cl_het, library(lmomco))
    
    # Worker defined at top level so it serializes against each node's own
    # .GlobalEnv (same mechanism clusterExport() relies on) and can see the
    # exported objects above.
    .het_worker <- function(task) {
      cols <- ref_idx[mon_ref == task$month]
      site_samples <- lapply(task$keep, function(i) {
        x <- swe_matrix[i, cols] / index_values[i, task$month]
        x[is.finite(x) & x > 0]
      })
      ht <- tryCatch(
        heterogeneity_test(site_samples, Nsim = HETEROGENEITY_NSIM),
        error = function(e) list(H0 = NA_real_, V0 = NA_real_, status = "error")
      )
      data.frame(region = task$region, month = task$month, H0 = ht$H0,
                 status = ht$status, n_sites = length(task$keep))
    }
    
    # ---- Manual load-balanced dispatch with a live progress readout ----
    # parLapply()/clusterApplyLB() block silently until the WHOLE batch is
    # done, which is exactly what produced the opaque stall. This reproduces
    # the same load-balanced scheduling using the low-level sendCall() /
    # recvOneResult() primitives that those functions are themselves built
    # on, so the master can print a line every time one (region, month) test
    # finishes instead of waiting on all 240 at once.
    n_tasks    <- length(het_tasks)
    n_workers  <- length(cl_het)
    results    <- vector("list", n_tasks)
    start_time <- Sys.time()
    
    n_initial <- min(n_workers, n_tasks)
    for (i in seq_len(n_initial)) {
      parallel:::sendCall(cl_het[[i]], .het_worker, list(het_tasks[[i]]), tag = i)
    }
    next_task <- n_initial + 1L
    completed <- 0L
    
    tryCatch({
      while (completed < n_tasks) {
        res <- parallel:::recvOneResult(cl_het)
        results[[res$tag]] <- res$value
        completed <- completed + 1L
        
        if (next_task <= n_tasks) {
          parallel:::sendCall(cl_het[[res$node]], .het_worker, list(het_tasks[[next_task]]), tag = next_task)
          next_task <- next_task + 1L
        }
        
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        rate    <- completed / elapsed
        eta_sec <- if (is.finite(rate) && rate > 0) (n_tasks - completed) / rate else NA_real_
        cat(sprintf("\r  Progress: %4d/%4d (%5.1f%%) | elapsed %5.0fs | est. remaining %5.0fs   ",
                    completed, n_tasks, 100 * completed / n_tasks, elapsed,
                    ifelse(is.na(eta_sec), 0, eta_sec)))
        flush.console()
      }
    }, error = function(e) { try(stopCluster(cl_het), silent = TRUE); stop(e) })
    cat("\n")
    
    stopCluster(cl_het)
    heterogeneity_log <- do.call(rbind, results)
  } else {
    heterogeneity_log <- data.frame()
  }
  
  
  # ---- E.1 Turn the H0 flag into an actual pooling decision ----
  # Previously heterogeneity_log was written to disk and printed as a
  # diagnostic but never actually consulted when fitting the regional gamma
  # (every region/month pooled regardless of H0). Region/months flagged
  # H0 >= H0_REGION_FLAG_THRESHOLD are now excluded from regional pooling
  # entirely and fall through to per-pixel local fitting in
  # monthly_sspi_regional() instead - see heterogeneous_mat below and its use
  # in section G. NA/failed H0 tests (lmom_failed, fit_failed, sim_failed,
  # too_few_sites) are treated as heterogeneous too, since "unknown" is not
  # the same as "known homogeneous" and pooling on an untested sample
  # defeats the purpose of the screen.
  heterogeneous_mat <- matrix(TRUE, n_regions, 12)  # default: don't pool
  if (nrow(heterogeneity_log) > 0) {
    for (k in seq_len(nrow(heterogeneity_log))) {
      r <- heterogeneity_log$region[k]; m <- heterogeneity_log$month[k]
      h <- heterogeneity_log$H0[k]
      heterogeneous_mat[r, m] <- !(is.finite(h) && h < H0_REGION_FLAG_THRESHOLD)
    }
    n_flag <- sum(heterogeneity_log$H0 >= H0_REGION_FLAG_THRESHOLD, na.rm = TRUE)
    n_untested <- sum(!is.finite(heterogeneity_log$H0))
    cat(sprintf("??? Heterogeneity test: %d/%d region-months flagged H0 >= %.0f (possibly/definitely heterogeneous).\n",
                n_flag, nrow(heterogeneity_log), H0_REGION_FLAG_THRESHOLD))
    cat(sprintf("  %d additional region-months had no valid H0 (untested) and are also excluded from pooling.\n",
                n_untested))
    cat("  Flagged/untested region-months are EXCLUDED from regional pooling and fall back\n")
    cat("  to per-pixel local fitting (monthly_sspi) for that region/month - no regional-\n")
    cat("  empirical fallback is used. See rfa_diagnostics/ for the full H0 breakdown.\n")
  }
  
  ## ---- F. Bias correction vs. independent station SWE observations ----
  # Station correction already applied basin-wide (once) upstream, directly
  # to swe_matrix, before this RFA block runs - see the "BASIN-WIDE BIAS
  # CORRECTION" section above. Re-running fit_bias_correction() here would
  # double-correct every pixel, so it is intentionally NOT called. C_matrix
  # is fixed at identity so monthly_sspi_regional()'s `mu_corr <- mu * c_val`
  # step remains a no-op. fit_bias_correction() is left defined earlier in
  # the script (unused) in case a denser station network later makes
  # per-region quantile mapping viable.
  C_matrix <- matrix(1, nrow = max(n_regions, 1), ncol = 12)
  bias_correction_log <- data.frame()
  cat("??? Per-region station bias correction: skipped (basin-wide correction already applied upstream; see basin_bias_log).\n")
  
  ## ---- G. Pooled regional gamma fit + KS goodness-of-fit per region/month ----
  cat("\n===== FITTING POOLED REGIONAL GAMMA DISTRIBUTIONS =====\n")
  region_par <- vector("list", n_regions)
  gof_log <- data.frame()
  for (r in seq_len(n_regions)) {
    region_par[[r]] <- vector("list", 12)
    pix_r <- which(region_id == r)
    for (m in 1:12) {
      keep <- pix_r[!discordant_mat[pix_r, m]]
      cols <- ref_idx[mon_ref == m]
      if (length(keep) == 0 || length(cols) == 0) {
        region_par[[r]][[m]] <- list(par = NULL, pooled_sample = numeric(0),
                                     gof_pass = FALSE, status = "no_data")
        # Previously silent (no gof_log row) - logging it so gof_log always
        # has n_regions*12 rows and "missing" region-months don't have to be
        # inferred by diffing row counts against n_regions*12.
        gof_log <- rbind(gof_log, data.frame(region = r, month = m, n_pool = 0L,
                                             fit_status = "no_data", ks_D = NA_real_,
                                             ks_p = NA_real_, gamma_used = FALSE))
        next
      }
      if (isTRUE(heterogeneous_mat[r, m])) {
        # Region/month failed (or never passed) the H0 heterogeneity screen -
        # do NOT pool. Leaving par = NULL and pooled_sample empty means
        # pool_ok is FALSE in monthly_sspi_regional(), so every pixel in this
        # region/month falls through to per-pixel local fitting instead of a
        # regional gamma OR a regional-empirical CDF.
        region_par[[r]][[m]] <- list(par = NULL, pooled_sample = numeric(0),
                                     gof_pass = FALSE, status = "heterogeneous_skip")
        # n_pool here is the pooled SAMPLE size the region/month WOULD have
        # had (matching the units of every other row's n_pool), not the site
        # count - a prior version of this log logged length(keep) (site
        # count, e.g. 36) here instead of the actual pooled sample size (e.g.
        # ~450-3000), which made heterogeneous_skip rows look like they had
        # far less data than they actually did.
        c_val_diag <- C_matrix[r, m]
        sub_diag  <- swe_matrix[keep, cols, drop = FALSE] * c_val_diag
        idxv_diag <- index_values[keep, m] * c_val_diag
        x_star_diag <- sweep(sub_diag, 1, idxv_diag, "/")
        n_pool_diag <- sum(is.finite(x_star_diag) & x_star_diag > 0)
        gof_log <- rbind(gof_log, data.frame(region = r, month = m, n_pool = n_pool_diag,
                                             fit_status = "heterogeneous_skip", ks_D = NA_real_,
                                             ks_p = NA_real_, gamma_used = FALSE))
        next
      }
      c_val <- C_matrix[r, m]
      sub  <- swe_matrix[keep, cols, drop = FALSE] * c_val
      idxv <- index_values[keep, m] * c_val
      x_star <- sweep(sub, 1, idxv, "/")
      x_star <- x_star[is.finite(x_star) & x_star > 0]
      
      fit <- fit_regional_gamma(x_star)
      gof <- if (!is.null(fit$par)) {
        gof_test_gamma(x_star, fit$par, alpha = GOF_ALPHA)
      } else list(D = NA_real_, p_value = NA_real_, pass = FALSE)
      
      region_par[[r]][[m]] <- list(par = fit$par, pooled_sample = x_star,
                                   gof_pass = isTRUE(gof$pass), status = fit$status,
                                   n_pool = length(x_star))
      gof_log <- rbind(gof_log, data.frame(region = r, month = m, n_pool = length(x_star),
                                           fit_status = fit$status, ks_D = gof$D,
                                           ks_p = gof$p_value, gamma_used = isTRUE(gof$pass)))
    }
  }
  if (nrow(gof_log) > 0) {
    n_gamma <- sum(gof_log$gamma_used, na.rm = TRUE)
    n_het_skip <- sum(gof_log$fit_status == "heterogeneous_skip", na.rm = TRUE)
    cat(sprintf("??? Regional gamma fitting: %d/%d region-months pass KS goodness-of-fit (alpha=%.2f)\n",
                n_gamma, nrow(gof_log), GOF_ALPHA))
    cat(sprintf("  Region-months excluded from pooling for failing the H0 heterogeneity screen: %d\n",
                n_het_skip))
    # No regional-empirical fallback: every region-month that fails GoF or the
    # H0 screen falls straight through to per-pixel local fitting instead
    # (see monthly_sspi_regional / "Regional Empirical" is always 0% below).
    cat(sprintf("  Region-months failing GoF (but passing H0) that fall back to LOCAL fitting: %d\n",
                sum(!gof_log$gamma_used & gof_log$fit_status != "heterogeneous_skip" &
                      gof_log$n_pool >= 15, na.rm = TRUE)))
  }
  
  ## ---- H. Diagnostic outputs (mirrors paper's Table 3/6 + Fig. 6) ----
  cat("\n===== WRITING RFA DIAGNOSTIC OUTPUTS =====\n")
  rfa_dir <- file.path(out_dir, "rfa_diagnostics")
  if (!dir.exists(rfa_dir)) dir.create(rfa_dir, recursive = TRUE)
  if (nrow(discordancy_log) > 0)
    write.csv(discordancy_log, file.path(rfa_dir, "discordancy_flagged_pixels.csv"), row.names = FALSE)
  if (nrow(heterogeneity_log) > 0)
    write.csv(heterogeneity_log, file.path(rfa_dir, "heterogeneity_H0_by_region_month.csv"), row.names = FALSE)
  if (nrow(gof_log) > 0)
    write.csv(gof_log, file.path(rfa_dir, "gamma_goodness_of_fit_by_region_month.csv"), row.names = FALSE)
  if (nrow(bias_correction_log) > 0)
    write.csv(bias_correction_log, file.path(rfa_dir, "bias_correction_coefficients.csv"), row.names = FALSE)
  
  region_raster <- rast(swe[[1]])
  values(region_raster) <- region_id
  png(file.path(rfa_dir, "regions_map.png"), width = 1200, height = 900, res = 150)
  plot(region_raster, main = sprintf("SSPI-1 homogeneous regions (%s)", regions_built$method),
       col = grDevices::hcl.colors(max(n_regions, 1), "Dark 3"), axes = FALSE, box = FALSE)
  if (!is.null(basin)) plot(basin, add = TRUE, border = "black", lwd = 1.5)
  dev.off()
  cat(sprintf("??? RFA diagnostics written to: %s\n", normalizePath(rfa_dir)))
  cat("============================================================\n")
  
} else if (scf_method == "2" && !REGIONALIZE_SSPI) {
  cat("\n??? REGIONALIZE_SSPI = FALSE: using original pixel-by-pixel SSPI-1 fitting.\n")
  region_id <- NULL; index_values <- NULL; C_matrix <- NULL
  region_par <- NULL; discordant_mat <- NULL
}

# ==============================================================================
#   MAIN CALCULATION LOOP - PARALLEL 
# ==============================================================================
timescales <- if (scf_method == "2") c(1) else c(3)

# ---- Set up parallel cluster ----
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\n??? Starting parallel cluster with %d cores\n", n_cores))
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, iseed = 40)  
if (scf_method == "2" && REGIONALIZE_SSPI) {
  clusterExport(cl, varlist = c("monthly_sspi", "monthly_sspi_regional", "clip_prob",
                                "swe_matrix", "dates", "ref_idx",
                                "region_id", "index_values", "C_matrix",
                                "region_par", "discordant_mat"),
                envir = environment())
  clusterEvalQ(cl, { library(zoo); library(lmomco) })
} else if (scf_method == "2") {
  clusterExport(cl, varlist = c("monthly_sspi", "clip_prob",
                                "swe_matrix", "dates", "ref_idx"),
                envir = environment())
  clusterEvalQ(cl, { library(zoo); library(lmomco) })
} else {
  clusterExport(cl, varlist = c("gringorten_swei_seasonal", "clip_prob",
                                "dates", "scf_mask_list",
                                "basin_scf_mask"),
                envir = environment())
  clusterEvalQ(cl, { library(zoo) })
}

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# ---- Summary file setup ----
summary_combined_file <- file.path(out_dir,
                                   if (scf_method == "2") "sspi_summary.txt" else "swei_all_timescales_summary.txt")
all_summaries <- list()

for (sc in timescales) {
  
  if (scf_method == "2") {
    # ============================================================================
    #   OPTION 2: SSPI-1 CALCULATION
    # ============================================================================
    cat(sprintf("\n===== SSPI-%d (Monthly, Parallel) =====\n", sc))
    cat(sprintf("Processing %d pixels in parallel (%s)... ", n_pixels,
                if (REGIONALIZE_SSPI) "regional pooled gamma" else "pixel-by-pixel gamma"))
    start_time <- Sys.time()
    
    if (REGIONALIZE_SSPI) {
      pixel_list <- parLapply(cl, seq_len(n_pixels), function(i) {
        tryCatch(
          monthly_sspi_regional(swe_matrix[i, ], dates, ref_idx, pixel_idx = i,
                                region_vec = region_id, index_values = index_values,
                                C = C_matrix, region_par = region_par,
                                discordant_mat = discordant_mat, eps = 1e-6),
          error = function(e) list(sspi = rep(NA_real_, length(dates)),
                                   method = rep(NA_integer_, length(dates)))
        )
      })
    } else {
      pixel_list <- parLapply(cl, seq_len(n_pixels), function(i) {
        tryCatch(
          monthly_sspi(swe_matrix[i, ], dates, ref_idx, eps = 1e-6),
          error = function(e) list(sspi = rep(NA_real_, length(dates)),
                                   method = rep(NA_integer_, length(dates)))
        )
      })
    }
    
    sspi_indices  <- do.call(rbind, lapply(pixel_list, `[[`, "sspi"))
    method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
    end_time <- Sys.time()
    cat(sprintf("done (%.1f min)\n", as.numeric(difftime(end_time, start_time, units = "mins"))))
    
    if (is.null(dim(sspi_indices)) || nrow(sspi_indices) != nrow(swe_matrix))
      stop("SSPI-1 calculation failed: unexpected result dimensions")
    sspi_indices[!is.finite(sspi_indices)] <- NA_real_
    
    clipped_low  <- sum(sspi_indices < -4.9, na.rm = TRUE)
    clipped_high <- sum(sspi_indices >  4.9, na.rm = TRUE)
    cat(sprintf("  Clipping: %d dry, %d wet\n", clipped_low, clipped_high))
    
    basin_mask_sspi <- !is.na(swe_matrix[, 1])
    na_rate_basin   <- 100 * mean(is.na(sspi_indices[basin_mask_sspi, ]))
    cat(sprintf("  NA rate (basin): %.3f%%\n", na_rate_basin))
    basin_avg_swei_results[[as.character(sc)]] <- colMeans(sspi_indices[basin_mask_sspi, , drop = FALSE], na.rm = TRUE)
    
    if (REGIONALIZE_SSPI) {
      n_valid_m   <- sum(!is.na(method_matrix))
      reg_gamma_pct  <- 100 * sum(method_matrix == 1, na.rm = TRUE) / n_valid_m
      reg_emp_pct    <- 100 * sum(method_matrix == 2, na.rm = TRUE) / n_valid_m
      local_gamma_pct <- 100 * sum(method_matrix == 3, na.rm = TRUE) / n_valid_m
      local_emp_pct   <- 100 * sum(method_matrix == 4, na.rm = TRUE) / n_valid_m
      local_zv_pct    <- 100 * sum(method_matrix == 5, na.rm = TRUE) / n_valid_m
      cat(sprintf("  Methods: Regional-Gamma=%.1f%%, Regional-Empirical=%.1f%%, Local-Gamma(fallback)=%.1f%%, Local-Empirical(fallback)=%.1f%%, Local-ZeroVar(fallback)=%.1f%%\n",
                  reg_gamma_pct, reg_emp_pct, local_gamma_pct, local_emp_pct, local_zv_pct))
      gamma_pct <- reg_gamma_pct; empirical_pct <- reg_emp_pct
      zero_var_pct <- local_zv_pct; excluded_pct <- local_gamma_pct + local_emp_pct
      all_summaries[[as.character(sc)]] <- list(
        sc = sc, date = Sys.time(), na_rate = na_rate_basin,
        gamma_pct = reg_gamma_pct, empirical_pct = reg_emp_pct,
        zero_var_pct = local_zv_pct, excluded_pct = local_gamma_pct + local_emp_pct,
        regional_gamma_pct = reg_gamma_pct, regional_empirical_pct = reg_emp_pct,
        local_gamma_pct = local_gamma_pct, local_empirical_pct = local_emp_pct,
        local_zerovar_pct = local_zv_pct,
        sspi_indices = sspi_indices, method_matrix = method_matrix
      )
    } else {
      gamma_pct      <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
      empirical_pct  <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
      zero_var_pct   <- 100 * sum(method_matrix == 3, na.rm = TRUE) / sum(!is.na(method_matrix))
      excluded_pct   <- 0
      cat(sprintf("  Methods: Gamma=%.1f%%, Empirical=%.1f%%, Zero-var=%.1f%%\n",
                  gamma_pct, empirical_pct, zero_var_pct))
      all_summaries[[as.character(sc)]] <- list(
        sc = sc, date = Sys.time(), na_rate = na_rate_basin,
        gamma_pct = gamma_pct, empirical_pct = empirical_pct,
        zero_var_pct = zero_var_pct, excluded_pct = excluded_pct,
        sspi_indices = sspi_indices, method_matrix = method_matrix
      )
    }
    
    months_all_sc <- as.numeric(format(dates, "%m"))
    
    cat("  -> Saving CSV files...\n")
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sub  <- sspi_indices[, idx_m, drop = FALSE]
      df   <- as.data.frame(sub)
      colnames(df) <- format(dates[idx_m], "%Y")
      coords <- xyFromCell(swe[[1]], 1:nrow(sub))
      df <- cbind(lon = coords[,1], lat = coords[,2], df)
      write.csv(df, file.path(out_dir, sprintf("sspi_%02d_month%02d_%s.csv", sc, m, month_names[m])), row.names = FALSE)
    }
    
    cat("  -> Saving Excel workbook...\n")
    excel_data <- list()
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sub  <- sspi_indices[, idx_m, drop = FALSE]
      df   <- as.data.frame(sub)
      colnames(df) <- format(dates[idx_m], "%Y")
      coords <- xyFromCell(swe[[1]], 1:nrow(sub))
      excel_data[[month_names[m]]] <- cbind(lon = coords[,1], lat = coords[,2], df)
    }
    xlsx_file <- file.path(out_dir, sprintf("sspi_%02d_all_months.xlsx", sc))
    if (file.exists(xlsx_file)) file.remove(xlsx_file)
    write_xlsx(excel_data, xlsx_file)
    
    cat("  -> Creating distribution map...\n")
    if (REGIONALIZE_SSPI) {
      # Code 2 "Regional Empirical" is retained only for backward-compatible
      # numbering; it is never assigned (no regional-empirical fallback is
      # used - see monthly_sspi_regional). It will always show 0 pixels below.
      method_mapping <- data.frame(
        code  = c(1, 2, 3, 4, 5),
        name  = c("Regional Gamma", "Regional Empirical (unused)", "Local Gamma (fallback)",
                  "Local Empirical (fallback)", "Local Zero-Var (fallback)"),
        color = c("#4575b4", "#74add1", "#fdae61", "#d73027", "#91bfdb"),
        stringsAsFactors = FALSE
      )
    } else {
      method_mapping <- data.frame(
        code  = c(1, 2, 3, 4),
        name  = c("Gamma", "Empirical", "Zero-Variance", "Excluded"),
        color = c("#4575b4", "#d73027", "#91bfdb", "#cccccc"),
        stringsAsFactors = FALSE
      )
    }
    method_raster <- rast(swe[[1]])
    dist_vals <- method_matrix[, 1]
    dist_vals[is.na(values(swe[[1]]))] <- NA
    values(method_raster) <- dist_vals
    n_cat <- nrow(method_mapping)
    png_file <- file.path(out_dir, sprintf("sspi_%02d_method_map.png", sc))
    png(png_file, width = 1200, height = 800, res = 150)
    plot(method_raster, col = method_mapping$color, breaks = seq(0.5, n_cat + 0.5, by = 1),
         legend = FALSE, main = sprintf("SSPI-%d: Fitting Method Distribution", sc),
         axes = FALSE, box = FALSE)
    if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
    dist_counts <- table(factor(dist_vals[!is.na(dist_vals)], levels = 1:n_cat,
                                labels = method_mapping$name))
    legend("bottomright", legend = sprintf("%s (%d pixels)", names(dist_counts), dist_counts),
           fill = method_mapping$color, title = "Method", cex = 0.8,
           bg = "white", bty = "o", border = "gray50")
    valid_pix <- sum(!is.na(values(swe[[1]])))
    mtext(sprintf("Valid pixels: %d | Gamma: %.1f%%", valid_pix, gamma_pct),
          side = 1, line = 0.5, cex = 0.7, col = "gray30")
    dev.off()
    
    cat("  -> Saving NetCDF files...\n")
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sspi_rast <- rep(rast(swe[[1]]), length(idx_m))
      values(sspi_rast) <- sspi_indices[, idx_m, drop = FALSE]
      terra::time(sspi_rast) <- dates[idx_m]
      writeCDF(sspi_rast, file.path(out_dir,
                                    sprintf("sspi_%02d_month%02d_%s.nc", sc, m, month_names[m])),
               varname = "sspi",
               longname = sprintf("Standardized SnowPack Index (SSPI-%d, monthly)", sc),
               unit = "standardized_index", missval = -9999, overwrite = TRUE)
    }
    cat(sprintf("??? SSPI-%d complete\n", sc))
    gc()
    next  
  }
  
  # ============================================================================
  #   OPTION 1: SWEI CALCULATION
  # ============================================================================
  cat(sprintf("\n===== SWEI-%d (%d-month window, Sequential Calculation) =====\n", sc, sc))
  
  swe_basin_avg_smoothed <- zoo::rollapply(swe_basin_avg_raw,
                                           width = sc,
                                           FUN   = sum,
                                           align = "right",
                                           fill  = NA,
                                           na.rm = FALSE)
  cat(sprintf("??? Basin-averaged %d-month smoothed SWE: %d steps, %.1f%% non-NA\n",
              sc, length(swe_basin_avg_smoothed),
              100 * mean(!is.na(swe_basin_avg_smoothed))))
  
  swe_smoothed <- matrix(NA, nrow = nrow(swe_matrix), ncol = ncol(swe_matrix))
  for (i in 1:nrow(swe_matrix)) {
    swe_smoothed[i, ] <- zoo::rollapply(swe_matrix[i, ],
                                        width = sc,
                                        FUN = sum,
                                        align = "right",
                                        fill = NA,
                                        na.rm = FALSE) 
  }
  cat(sprintf("??? %d-month rolling sum applied to SWE (SWEI preprocessing)\n", sc))
  
  clusterExport(cl, varlist = c("swe_smoothed"), envir = environment())
  
  avg_swei_out <- tryCatch(
    gringorten_swei_seasonal(swe_basin_avg_smoothed, dates, scf_mask_list,
                             cell_index = NULL, basin_scf_mask = basin_scf_mask,
                             eps = 1e-6),
    error = function(e) list(swei = rep(NA_real_, length(dates)), method = rep(NA_integer_, length(dates)))
  )
  basin_avg_swei_results[[as.character(sc)]] <- avg_swei_out$swei
  
  cat(sprintf("Processing %d pixels in parallel... ", n_pixels))
  start_time <- Sys.time()
  
  pixel_list <- parLapply(cl, seq_len(n_pixels), function(i) {
    tryCatch(
      gringorten_swei_seasonal(swe_smoothed[i, ], dates, scf_mask_list, cell_index = i, eps = 1e-6),
      error = function(e) list(swei = rep(NA_real_, length(dates)),
                               method = rep(NA_integer_, length(dates)))
    )
  })
  
  swei_indices <- do.call(rbind, lapply(pixel_list, `[[`, "swei"))
  method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
  
  extreme_dry  <- sum(swei_indices < -3.0, na.rm = TRUE)
  extreme_wet  <- sum(swei_indices >  3.0, na.rm = TRUE)
  cat(sprintf("  Extreme values (|SWEI| > 3): %d dry, %d wet (retained, no clipping)\n",
              extreme_dry, extreme_wet))
  
  basin_pix_rows      <- which(!is.na(swe_matrix[, 1]))
  
  # Protect against empty vectors generating NaN
  if (length(basin_pix_rows) == 0) {
    na_rate_grid <- NA
  } else {
    na_rate_grid <- 100 * mean(is.na(swei_indices[basin_pix_rows, , drop = FALSE]))
  }
  
  basin_mean_complete <- 100 * mean(!is.na(basin_avg_swei_results[[as.character(sc)]]))
  
  cat(sprintf("  Basin mean series:  %.1f%% complete (%.1f%% NA)\n",
              basin_mean_complete, 100 - basin_mean_complete))
  cat(sprintf("  Per-pixel gridded:  %.3f%% NA  (basin pixels only)\n", na_rate_grid))
  
  grid_status <- if (is.na(na_rate_grid)) {
    "UNKNOWN (0 valid basin pixels)" 
  } else if (na_rate_grid < 0.5) {
    "READY  (<0.5% NAs)"
  } else {
    "CAUTION (>0.5% NAs)"
  }
  
  cat(sprintf("  Gridded compatibility: %s\n", grid_status))
  
  if (!is.na(na_rate_grid) && na_rate_grid >= 0.5) {
    cat("  ???  High per-pixel NA is expected: the 75%%-nonzero-years criterion\n")
    cat("     and SCF mask are applied pixel-by-pixel, which is much stricter\n")
    cat("     than the basin-mean aggregate implies. The basin mean series is\n")
    cat("     suitable for time-series analysis; use gridded files for spatial\n")
    cat("     analysis with awareness of actual spatial coverage.\n")
  }
  
  gringorten_pct  <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
  masked_pct      <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
  
  cat(sprintf("  Methods: Gringorten=%.1f%%, Masked/NA=%.1f%%\n",
              gringorten_pct, masked_pct))
  
  all_summaries[[as.character(sc)]] <- list(
    sc = sc, date = Sys.time(), na_rate = na_rate_grid, 
    basin_mean_complete = basin_mean_complete,
    gringorten_pct = gringorten_pct, masked_pct = masked_pct,
    swei_indices = swei_indices, method_matrix = method_matrix
  )
  
  months_all <- as.numeric(format(dates, "%m"))
  cat("  -> Saving CSV files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    csv_file <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.csv", sc, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE)
  }
  
  cat("  -> Saving Excel workbook...\n")
  excel_data <- list()
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    excel_data[[as.character(m)]] <- df
  }
  xlsx_file <- file.path(out_dir, sprintf("swei_%02d_all_months.xlsx", sc))
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  write_xlsx(excel_data, xlsx_file)
  
  cat("  -> Creating distribution map...\n")
  swei_dist_mapping <- data.frame(
    code = c(1, 2), name = c("Gringorten", "Masked/NA"),
    color = c("#4575b4", "#91bfdb"), stringsAsFactors = FALSE
  )
  swei_dist_raster <- rast(swe[[1]])
  dist_sim <- apply(method_matrix, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_integer_)
    as.integer(names(sort(table(x), decreasing = TRUE))[1])
  })
  dist_sim[is.na(dist_sim)] <- NA
  values(swei_dist_raster) <- dist_sim
  png_file_swei <- file.path(out_dir, sprintf("swei_%02d_distribution_map.png", sc))
  png(png_file_swei, width = 1800, height = 1000, res = 150)
  layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
  plot(swei_dist_raster, col = swei_dist_mapping$color, breaks = c(0.5, 1.5, 2.5),
       legend = FALSE, main = sprintf("SWEI-%d: Seasonal Method", sc), axes = FALSE, box = FALSE)
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  valid_pixels <- sum(!is.na(values(swe[[1]])))
  mtext(sprintf("Valid pixels: %d | Gringorten: %.1f%%", valid_pixels, gringorten_pct),
        side = 1, line = 1, cex = 0.75, col = "gray30")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  dist_counts <- table(factor(dist_sim, levels = 1:2, labels = swei_dist_mapping$name))
  legend_entries <- sprintf("%s (%d pixels)", names(dist_counts), dist_counts)
  legend("center", legend = legend_entries, fill = swei_dist_mapping$color,
         title = "Method", cex = 0.95, bg = "white", bty = "o", border = "gray50")
  dev.off()
  layout(1)
  
  cat("  -> Saving NetCDF files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
    swei_rast  <- rep(rast(swe[[1]]), length(idx_m))
    sub_matrix  <- swei_indices[, idx_m, drop = FALSE]
    if (nrow(sub_matrix) != ncell(swei_rast)) {
      stop(sprintf("Dimension mismatch: Matrix %d rows vs Raster %d cells.",
                   nrow(sub_matrix), ncell(swei_rast)))
    }
    values(swei_rast)  <- sub_matrix
    terra::time(swei_rast)  <- dates[idx_m]
    nc_file  <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.nc", sc, m, month_names[m]))
    writeCDF(swei_rast, nc_file, varname = "swei",
             longname = sprintf("Seasonal SWEI (SWEI-%d)", sc),
             unit = "standardized_index", missval = -9999, overwrite = TRUE)
  }
  cat(sprintf("??? SWEI-%d complete\n", sc))
  gc()
}

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\n??? Parallel cluster stopped\n")

# ==============================================================================
#   SAVE BASIN-AVERAGED SWEI TO CSV 
# ==============================================================================
index_lbl <- if (scf_method == "2") "sspi" else "swei"
cat(sprintf("\n===== SAVING BASIN-AVERAGED %s =====\n", toupper(index_lbl)))
months_all  <- as.integer(format(dates, "%m"))
years_all   <- as.integer(format(dates, "%Y"))
month_names_save  <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
for (sc in timescales) {
  basin_series <- basin_avg_swei_results[[as.character(sc)]]
  yrs <- sort(unique(years_all))
  df_out <- data.frame(Year = yrs)
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    yr_m  <- years_all[idx_m]
    val_m <- basin_series[idx_m]
    col   <- setNames(val_m, yr_m)
    df_out[[month_names_save[m]]] <- col[as.character(yrs)]
  }
  csv_file <- file.path(out_dir, sprintf("%s_%02d_basin_averaged_by_month.csv", index_lbl, sc))
  write.csv(df_out, csv_file, row.names = FALSE, na = "")
  cat(sprintf("??? Saved basin-averaged %s-%d (12 monthly series) to: %s\n", toupper(index_lbl), sc, csv_file))
}

# ==============================================================================
#   WRITE SUMMARY FILE AT END 
# ==============================================================================
cat("\n===== WRITING SUMMARY FILE =====\n")

if (scf_method == "2") {
  cat(
    "============================================================\n",
    "MONTHLY SSPI-1 SUMMARY - NECHAKO BASIN\n",
    "ERA5-Land Snow Water Equivalent\n",
    "Method 2: Zero-Inflated Gamma (no SCF mask)\n",
    sprintf("Regional pooling (index-flood RFA): %s\n", if (REGIONALIZE_SSPI) "ENABLED" else "disabled"),
    "============================================================\n\n",
    file = summary_combined_file, sep = ""
  )
  for (sc in timescales) {
    s <- all_summaries[[as.character(sc)]]
    cat(sprintf("\n========== SSPI-%d ==========\n", s$sc), file = summary_combined_file, append = TRUE)
    cat(sprintf("Calculation date: %s\n", s$date), file = summary_combined_file, append = TRUE)
    cat(sprintf("Grid dimensions: %d x %d\n", ncol(swe), nrow(swe)), file = summary_combined_file, append = TRUE)
    cat(sprintf("Basin pixels (total): %d\n", basin_pixels), file = summary_combined_file, append = TRUE)
    cat(sprintf("Time period: %s to %s (%d months)\n", min(dates), max(dates), length(dates)), file = summary_combined_file, append = TRUE)
    cat(sprintf("Reference period: %s to %s\n", sspi_ref_start, sspi_ref_end), file = summary_combined_file, append = TRUE)
    cat("\nZero-Inflation Screening:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Excluded pixels (>50%% zero SWE in ref period): %d (%.1f%%)\n", length(invalid_pix_sspi), 100 * length(invalid_pix_sspi) / n_pixels), file = summary_combined_file, append = TRUE)
    cat("\nMethodology:\n", file = summary_combined_file, append = TRUE)
    if (REGIONALIZE_SSPI) {
      cat(sprintf("  - Regionalization: %s, %d homogeneous sub-regions (index-flood RFA, Fontana et al. 2026)\n",
                  regions_built$method, n_regions), file = summary_combined_file, append = TRUE)
      cat(sprintf("  - Per-pixel index value: %s of reference-period monthly SWE\n", INDEX_VALUE_STAT),
          file = summary_combined_file, append = TRUE)
      cat(sprintf("  - Discordancy screening (Hosking & Wallis D_i): %s, %d pixel-months excluded from pooling\n",
                  ifelse(DISCORDANCY_SCREENING, "ON", "OFF"), sum(discordant_mat)),
          file = summary_combined_file, append = TRUE)
      if (nrow(heterogeneity_log) > 0) {
        cat(sprintf("  - Heterogeneity (H0, Monte Carlo, Nsim=%d): %d/%d region-months flagged H0>=%.0f\n",
                    HETEROGENEITY_NSIM, sum(heterogeneity_log$H0 >= H0_REGION_FLAG_THRESHOLD, na.rm = TRUE),
                    nrow(heterogeneity_log), H0_REGION_FLAG_THRESHOLD),
            file = summary_combined_file, append = TRUE)
      }
      if (nrow(gof_log) > 0) {
        cat(sprintf("  - Pooled regional gamma GoF (KS, alpha=%.2f): %d/%d region-months pass\n",
                    GOF_ALPHA, sum(gof_log$gamma_used, na.rm = TRUE), nrow(gof_log)),
            file = summary_combined_file, append = TRUE)
      }
      cat(sprintf("  - Bias correction vs. station obs: %s (basin-wide, %d/12 months corrected; see basin_bias_log)\n",
                  ifelse(BIAS_CORRECT && file.exists(STATION_OBS_PATH), "ENABLED", "disabled (no station file)"),
                  sum(basin_bias_log$applied, na.rm = TRUE)),
          file = summary_combined_file, append = TRUE)
      cat("  - Per-pixel zero-probability retained (only the gamma shape/scale is pooled regionally)\n",
          file = summary_combined_file, append = TRUE)
      cat("  - Region/months failing GoF, OR failing the H0 heterogeneity screen, are excluded\n    from pooling entirely and fall back to per-pixel local fitting (no regional-\n    empirical CDF fallback is used)\n",
          file = summary_combined_file, append = TRUE)
      cat("  - Pixels with no usable region/index value fall back to the original per-pixel fit\n",
          file = summary_combined_file, append = TRUE)
      cat(sprintf("  - Diagnostics written to: %s\n", file.path(out_dir, "rfa_diagnostics")),
          file = summary_combined_file, append = TRUE)
    } else {
      cat("  - No SCF domain mask applied (SSPI-1 is valid for all months)\n  - Calendar-month stratified fitting\n  - Zero-inflated gamma (center of probability mass)\n  - L-moments (lmomco) for gamma parameter estimation\n  - Empirical fallback when gamma fit fails or variance too low\n  - SSPI values clipped at +/-5\n", file = summary_combined_file, append = TRUE)
    }
    if (REGIONALIZE_SSPI) {
      cat("\nMethod Distribution:\n", file = summary_combined_file, append = TRUE)
      cat(sprintf("  Regional gamma:             %.1f%%\n", s$regional_gamma_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Regional empirical (unused, no fallback CDF): %.1f%%\n", s$regional_empirical_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Local gamma (fallback):      %.1f%%\n", s$local_gamma_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Local empirical (fallback):  %.1f%%\n", s$local_empirical_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Local zero-variance (fallback): %.1f%%\n", s$local_zerovar_pct), file = summary_combined_file, append = TRUE)
    } else {
      cat("\nMethod Distribution:\n", file = summary_combined_file, append = TRUE)
      cat(sprintf("  Gamma fitting:   %.1f%%\n", s$gamma_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Empirical:       %.1f%%\n", s$empirical_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Zero-variance:   %.1f%%\n", s$zero_var_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Excluded:        %.1f%%\n", s$excluded_pct), file = summary_combined_file, append = TRUE)
    }
    cat("\nNA Analysis:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Total NA values in basin: %.3f%%\n", s$na_rate), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Compatibility: %s\n", ifelse(s$na_rate < 0.5, "READY (<0.5% NAs)", "CAUTION (>0.5% NAs)")), file = summary_combined_file, append = TRUE)
    fi <- 1:min(12, ncol(s$sspi_indices))
    fs <- s$sspi_indices[, fi, drop = FALSE]
    cat("\nSnow drought frequency (first 12 months):\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional deficit (SSPI < -2.0): %.1f%%\n", 100 * mean(fs < -2.0, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Extreme deficit     (SSPI < -1.5): %.1f%%\n", 100 * mean(fs < -1.5, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Moderate deficit    (SSPI < -1.0): %.1f%%\n", 100 * mean(fs < -1.0, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional surplus (SSPI >  2.0): %.1f%%\n", 100 * mean(fs >  2.0, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
  }
} else {
  cat(
    "============================================================\n",
    "COMBINED SWEI SUMMARY FOR ALL TIMESCALES\n",
    "Seasonal Approach: SWE rolling sum = timescale months (Option B)\n",
    "SCF mask window fixed at k=3 per Huning & AghaKouchak (2020)\n",
    "============================================================\n\n",
    file = summary_combined_file, sep = ""
  )
  for (sc in timescales) {
    s  <- all_summaries[[as.character(sc)]]
    cat(sprintf("\n\n========== SWEI-%d ==========\n", s$sc), file = summary_combined_file, append = TRUE)
    cat(sprintf("Calculation date: %s\n", s$date), file = summary_combined_file, append = TRUE)
    cat(sprintf("Grid dimensions: %d x %d\n", ncol(swe), nrow(swe)), file = summary_combined_file, append = TRUE)
    cat(sprintf("Basin pixels: %d\n", basin_pixels), file = summary_combined_file, append = TRUE)
    cat(sprintf("Time period: %s to %s (%d months)\n", min(dates), max(dates), length(dates)), file = summary_combined_file, append = TRUE)
    cat("\nMethodology:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  - %d-month rolling sum applied to RAW SWE (window = sc = %d)\n", sc, sc), file = summary_combined_file, append = TRUE)
    cat("  - SCF mask window fixed at k=3 (intentionally independent of SWE window,\n    per Huning & AghaKouchak 2020 ??? see target_k comment in script)\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  - SCF domain method: %s\n", scf_method_label), file = summary_combined_file, append = TRUE)
    cat(sprintf("    -> Fixed %.0f%% SCF threshold, right-aligned 3-month climatology\n", 100 * 0.05), file = summary_combined_file, append = TRUE)
    cat("  - Gringorten plotting position (non-parametric)\n  - Zero-SWE perturbation: random uniform in (0, min_nonzero) per Huning & AghaKouchak (2020)\n  - Minimum data requirement: >= 75% of years for a calendar month must have nonzero SWE\n  - No clipping applied to SWEI values (full distribution retained)\n", file = summary_combined_file, append = TRUE)
    cat("\nMethod Distribution:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Gringorten fitting: %.1f%%\n", s$gringorten_pct), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Masked/NA: %.1f%%\n", s$masked_pct), file = summary_combined_file, append = TRUE)
    cat("\nNA Analysis:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Basin mean series:     %.1f%% complete (%.1f%% NA)\n", s$basin_mean_complete, 100 - s$basin_mean_complete), file = summary_combined_file, append = TRUE)
    
    # Safely print na_rate even if NA
    if (is.na(s$na_rate)) {
      cat("  Per-pixel gridded:     UNKNOWN% NA  (0 basin pixels)\n", file = summary_combined_file, append = TRUE)
      cat("  Gridded compatibility: UNKNOWN (0 valid basin pixels)\n", file = summary_combined_file, append = TRUE)
    } else {
      cat(sprintf("  Per-pixel gridded:     %.3f%% NA  (basin pixels only)\n", s$na_rate), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Gridded compatibility: %s\n", ifelse(s$na_rate < 0.5, "READY  (<0.5% NAs)", "CAUTION (>0.5% NAs)")), file = summary_combined_file, append = TRUE)
    }
    
    if (!is.na(s$na_rate) && s$na_rate >= 0.5) {
      cat("  NOTE: The large gap between basin-mean completeness and per-pixel NA\n  is expected ??? spatial averaging masks individual pixel gaps, while the\n  75%%-nonzero-years criterion and SCF mask applied pixel-by-pixel are\n  much stricter. The basin mean series is suitable for time-series analysis;\n  use gridded files for spatial analysis with awareness of spatial coverage.\n", file = summary_combined_file, append = TRUE)
    }
    first_idx  <- 1:min(12, ncol(s$swei_indices))
    first_swei  <- s$swei_indices[, first_idx, drop = FALSE]
    cat("\nDrought frequency (first 12 months):\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional (SWEI < -2.0): %.1f%%\n", 100 * mean(first_swei < -2.0, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Extreme     (SWEI < -1.6): %.1f%%\n", 100 * mean(first_swei < -1.6, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Severe      (SWEI < -1.3): %.1f%%\n", 100 * mean(first_swei < -1.3, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Moderate    (SWEI < -0.8): %.1f%%\n", 100 * mean(first_swei < -0.8, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
  }
}

cat("\n============================================================\n")
cat(if (scf_method == "2") "SSPI-1 CALCULATION COMPLETE!\n" else "SWEI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat(sprintf("Summary file: %s\n", summary_combined_file))
cat("\n--- READY FOR ANALYSIS ---\n")
