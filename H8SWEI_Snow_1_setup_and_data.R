##############################################
# 01_setup_and_data.R
# STAGE 1 of 4 — Setup, configuration, and data loading
#
# This is the first of four scripts split from the original H8SWEI_Snow.R
# monolith so that a crash/error in a later, slower stage never forces you
# to re-run this (fast but fragile — typo'd paths, missing files) stage.
#
# WHAT THIS SCRIPT DOES:
#   - Loads libraries, resolves the SCF method (1 or 2), runs preflight checks
#   - Loads the basin boundary + ERA5-Land SWE/SCF NetCDFs
#   - Reprojects/resamples/masks to the basin, builds swe_matrix + dates
#   - Builds the SCF domain mask (method 1) or SSPI-1 zero-inflation screen (method 2)
#
# OUTPUT:
#   stage1.rds  — every object in the workspace at the end of this stage,
#                 written to STAGE_DIR (see below). This is the ONLY thing
#                 script 02 needs from this script.
#
# RUN ORDER: 01 -> 02 -> 03 -> 04 (each script loads the previous stage's .rds)
##############################################

# ---- Where staged .rds checkpoint files live ----
# Change this if you want the checkpoints somewhere other than the working
# directory. Each stage script creates this if it doesn't exist.
STAGE_DIR <- "stage_checkpoints"
if (!dir.exists(STAGE_DIR)) dir.create(STAGE_DIR, recursive = TRUE)

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
#       [Nechako SAF pipeline] Run this method (2) BEFORE A03_saf_core_pipeline.R.
#       Files are written as sspi_01_month{MM}_{Mon}.nc, which
#       assemble_index_raster(prefix="sspi", scale_tag="01") in
#       A03_saf_core_pipeline.R reads directly -- no renaming needed.
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

# ---- Total run timer ----
# Captured first, before anything else runs (including library loading), so
# the total elapsed time printed at the very end of the script reflects the
# full wall-clock cost of the run, not just the main calculation loop.
SCRIPT_START_TIME <- Sys.time()

# ---- Libraries ----
library(terra)
library(ncdf4)
library(zoo)
library(writexl)
library(parallel)
library(trend)
library(data.table)
library(lmomco)   # required for option 2 (SSPI-1 zero-inflated gamma)
library(ggplot2)  # basin-averaged time series plots (see TIME SERIES PLOTS section)
# also used for the regional L-moment / kappa-distribution
# machinery (discordancy, heterogeneity, regional gamma fit)
# stats::ks.test (base R) is used for the gamma goodness-of-fit test.

# ---- USER SELECTION: SCF DOMAIN METHOD ----
cat("\n============================================================\n")
cat("SCF DOMAIN DEFINITION: SELECT METHOD\n")
cat("============================================================\n")
cat("  1 = Huning & AghaKouchak (2020) - fixed 5% SCF threshold\n")
cat("                          (standard approach from the original global paper)\n")
cat("  2 = SSPI-1 (Standardized SnowPack Index, monthly)\n")
cat("                          (zero-inflated gamma, no SCF mask)\n\n")

# ---- To use in batch/scheduled jobs, set DEFAULT_SCF_METHOD to "1" or "2"
# NOTE: changed from "1" to "2" so that non-interactive / batch sourcing
# (e.g. as part of the automated Nechako SAF pipeline, which now consumes
# SSPI-1 as an additional input alongside SPI-1/SPEI-1/SPI-3/SPEI-3) produces
# SSPI-1 (sspi_results_monthly/) by default. Set back to "1" if you need the
# original SWEI (Huning & AghaKouchak 2020) output instead.
DEFAULT_SCF_METHOD <- "2"

scf_method <- ""

# Priority 1: command-line argument (e.g. Rscript script.R 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && trimws(args[1]) %in% c("1", "2")) {
  scf_method <- trimws(args[1])
  cat(sprintf("[OK] Method set from command-line argument: %s\n", scf_method))
  
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
    # User clicked Cancel - use the default
    scf_method <- DEFAULT_SCF_METHOD
    cat(sprintf("[OK] Dialog cancelled: using default Method %s\n", scf_method))
  } else {
    scf_method <- trimws(raw_input)
  }
  while (!scf_method %in% c("1", "2")) {
    raw_input <- rstudioapi::showPrompt(
      title   = "Invalid input - try again",
      message = "Please enter 1 or 2:",
      default = DEFAULT_SCF_METHOD
    )
    if (is.null(raw_input)) {
      scf_method <- DEFAULT_SCF_METHOD
      break
    }
    scf_method <- trimws(raw_input)
  }
  cat(sprintf("[OK] Method set from RStudio dialog: %s\n", scf_method))
  
  # Priority 3: plain readline() fallback (R console, non-RStudio interactive)
} else if (interactive()) {
  while (!scf_method %in% c("1", "2")) {
    scf_method <- trimws(readline(prompt = "Enter 1 or 2: "))
    if (!scf_method %in% c("1", "2")) cat("  Invalid input. Please enter 1 or 2.\n")
  }
  
  # Priority 4: non-interactive fallback default (batch / Rscript without args)
} else {
  scf_method <- DEFAULT_SCF_METHOD
  cat(sprintf("[OK] Non-interactive session detected: using default Method %s\n", scf_method))
  cat("  (Edit DEFAULT_SCF_METHOD at the top of the script to change this.)\n")
}

if (scf_method == "1") {
  cat("\n- Method selected: Huning & AghaKouchak (2020) fixed 5% SCF threshold\n")
} else {
  cat("\n- Method selected: SSPI-1 (zero-inflated gamma, monthly)\n")
}

scf_method_label <- if (scf_method == "1") "Huning & AghaKouchak (2020)" else "SSPI-1 (zero-inflated gamma, no SCF mask)"

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")

# ==============================================================================
# PREFLIGHT CHECK (merged from former preflight_check.R)
# ==============================================================================
# Confirms every file this script actually reads is present in THIS working
# directory before spending time on spatial alignment, masking, and the
# parallel cluster - so a missing/typo'd path fails in seconds with a clear
# message instead of partway into a long run, or silently degrading to a
# fallback method (see log_pipeline_fallback() further below).
{
  preflight_required <- list(
    "Basin boundary (hard requirement - script stops without it)" =
      "Spatial/Nechako_Basin.kmz",
    "ERA5-Land SWE (hard requirement - NOT produced by Nechako_Snow_Data_download.R; must come from a separate ERA5-Land acquisition step)" =
      "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
  )
  
  preflight_optional <- list(
    "Station SWE observations (Nechako_Snow_Data_download.R Part A output - enables bias correction; falls back to uncorrected ERA5-Land if missing)" =
      "daily.csv",
    "DEM (enables elevation x aspect regionalization; falls back to spatial k-means if missing)" =
      "Spatial/nechako_dem_mosaic.tif"
  )
  
  cat("============================================================\n")
  cat("PREFLIGHT CHECK\n")
  cat("Working directory:", normalizePath("."), "\n")
  cat("============================================================\n\n")
  
  cat("--- Required ---\n")
  preflight_missing <- character(0)
  for (label in names(preflight_required)) {
    path <- preflight_required[[label]]
    ok <- file.exists(path)
    cat(sprintf("  [%s] %-70s %s\n", if (ok) "OK" else "MISSING", path, if (ok) "" else paste0("- ", label)))
    if (!ok) preflight_missing <- c(preflight_missing, path)
  }
  
  cat("\n--- Optional (script will run without these, but with a fallback) ---\n")
  for (label in names(preflight_optional)) {
    path <- preflight_optional[[label]]
    ok <- file.exists(path)
    cat(sprintf("  [%s] %-30s %s\n", if (ok) "OK" else "will fall back", path, paste0("- ", label)))
  }
  cat("\n")
  
  if (length(preflight_missing) > 0) {
    cat("============================================================\n")
    cat("PREFLIGHT FAILED - missing required input(s):\n")
    for (p in preflight_missing) cat("  -", p, "\n")
    cat("This script will stop() on these itself, but failing here saves the\n")
    cat("time it would otherwise spend on spatial alignment first.\n")
    cat("============================================================\n")
    stop("Preflight check failed - see missing required input(s) above.", call. = FALSE)
  } else {
    cat("[OK] All required inputs present. Continuing with the SWEI/SSPI-1 run.\n\n")
  }
  
  rm(preflight_required, preflight_optional, preflight_missing, label, path, ok)
}
# ---- Enforce pipeline run order: download script -> Bridge -> H8 ----
# sar_depth_monthly.nc and snow_data/cross_validation/*_monthly_basinmean.csv
# are NOT produced by Nechako_Snow_Data_download.R (see header note above) -
# they only exist after Nechako_Snow_Bridge.R has run. Both are optional
# accuracy enhancements (H8 degrades gracefully without them), so this is a
# warning, not a stop - but it's better to know NOW that they're missing
# because the Bridge script was skipped, rather than discover it later as an
# unexplained "SAR_DEPTH_PATH not found" fallback deep in the log.
bridge_outputs_missing <- !file.exists("sar_depth_monthly.nc") &&
  !dir.exists("snow_data/cross_validation")
if (bridge_outputs_missing) {
  cat("[NOTE] Neither sar_depth_monthly.nc nor snow_data/cross_validation/ exist yet.\n")
  cat("       These are built by Nechako_Snow_Bridge.R, not by Nechako_Snow_Data_download.R.\n")
  cat("       If you intend to use SAR-augmented bias correction and/or the independent-\n")
  cat("       product cross-validation, run Nechako_Snow_Bridge.R now, then re-source this script.\n")
  cat("       Proceeding without them - both fall back gracefully (see PIPELINE_FALLBACK_ALERTS).\n")
}
# ==============================================================================
# END PREFLIGHT CHECK
# ==============================================================================

out_dir <- if (scf_method == "2") "sspi_Snow_results_monthly" else "swei_results_seasonal"
TIMESTAMP_OUTPUT_SUBDIR <- FALSE  # set TRUE to give every run its own timestamped folder
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
REGIONALIZE_SSPI <- TRUE

# --- Regionalization (homogeneous sub-region) settings ---
# Optional DEM raster used to build elevation-band x aspect-class regions,
# analogous to the paper's A1-A6 avalanche-occurrence regions. If this file
# does not exist, the script automatically falls back to a spatial k-means
# regionalization on pixel coordinates (still pools pixels, but without a
# physical elevation/aspect basis) and prints a warning.

# NOTE: Nechako_Snow_Data_download.R Part A now writes its validated,
# mosaicked DEM to Spatial/nechako_dem_mosaic.tif after building it (whether
# from the NRCan CDEM tiles, the ESRI Grid fallback, or an elevatr download) -
# see the "Persist the DEM to a fixed path" step in that script. This DEM_PATH
# must be that SAME file (file.exists() is a literal string match - a
# previous version of this comment claimed the two scripts already agreed on
# a path, but they didn't: Part A never wrote ANY DEM file to disk at all, so
# this always silently fell back to the unweighted k-means regionalization
# below, regardless of what DEM_PATH said).
#
# It's deliberately the FULL 50 km-buffer mosaic Part A builds (to screen
# neighbouring-basin stations too), not a basin-only crop - build_regions()
# below already does `dem <- crop(dem, project(basin, crs(dem)))` internally
# before building elevation x aspect regions, so handing it the wider raster
# is exactly right and no second, basin-only DEM needs to be produced.
DEM_PATH          <- "Spatial/nechako_dem_mosaic.tif"
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
HETEROGENEITY_NSIM    <- 1000  # Monte-Carlo simulations for H0. 200 was too noisy near H0=1;
# 1000 is a compromise between Monte-Carlo stability and runtime. The script
# records Monte-Carlo uncertainty and uses a conservative decision band around
# H0=1 rather than treating a single noisy estimate as an exact threshold.
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
STATION_OBS_PATH <- "snow_data/daily.csv"  # produced by Nechako_Snow_Data_download.R Part A
# (spatially/topographically screened ASWS stations, columns:
# station_id, lon, lat, date, swe_mm)
BIAS_CORRECT     <- TRUE
BIAS_CORRECT_CLIP <- c(0.5, 2.0)  # sanity clip on correction coefficients C,
# mirroring the paper's flag of |B(TR)|>0.5
# as the threshold for "correction needed"

# --- Strict-mode switches: turn a silent degrade into a loud one (or a hard stop) ---
# By default this pipeline is permissive: if DEM_PATH or STATION_OBS_PATH can't
# be found, it falls back (unweighted spatial k-means regions / no bias
# correction) rather than dying, which is convenient for exploratory runs
# where those optional inputs genuinely aren't available yet. The risk is
# that a TYPO'd path degrades results quietly - the previous version of this
# script only ever logged this to the console via cat(), which is easy to
# miss in a long batch log.
#
# Two things now change that, regardless of these flags:
#   1. Every fallback is logged via log_pipeline_fallback() below, which
#      prints a hard-to-miss console banner, raises a real R warning()
#      (so it also shows up in R's end-of-run "Warning messages" list and in
#      any log-scraping that greps for "Warning"), and is written to
#      <out_dir>/PIPELINE_FALLBACK_ALERTS.txt - a dedicated file that always
#      gets created, even when empty, so "the file exists and says nothing
#      happened" is itself a positive confirmation, not just an absence.
#   2. run_settings.txt (written just below) now records whether the DEM and
#      station files were actually found, at the very top of the run, before
#      any silent degrade could happen.
#
# Set either flag below to TRUE to additionally turn that specific fallback
# into a hard stop() - use this for a "production" run where you specifically
# expect the DEM and/or station bias correction to be used, and would rather
# the script fail loudly than complete on the wrong inputs.
REQUIRE_DEM         <- FALSE  # TRUE = stop() if DEM_PATH doesn't resolve (instead of k-means fallback)
REQUIRE_STATION_OBS <- FALSE  # TRUE = stop() if STATION_OBS_PATH doesn't resolve and BIAS_CORRECT is TRUE

# Accumulates every fallback message triggered during the run (see
# log_pipeline_fallback()) so it can be written out as a single file at the
# end, in addition to being surfaced immediately when it happens.
PIPELINE_FALLBACK_ALERTS <- character(0)

# ---- Central alert helper: call this at every point where the pipeline is
# about to silently substitute a lesser data source / method for the one the
# user configured. Never call cat() alone for this kind of event again - use
# this instead, so the visibility behaviour (console banner + warning() +
# alerts file + optional hard stop) stays consistent everywhere it's used.
log_pipeline_fallback <- function(msg, hard_stop = FALSE) {
  banner <- paste0(
    "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
    "!!  FALLBACK TRIGGERED - RESULTS MAY BE DEGRADED\n",
    "!!  ", msg, "\n",
    "!!  Logged to: ", file.path(out_dir, "PIPELINE_FALLBACK_ALERTS.txt"), "\n",
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  )
  cat(banner)
  PIPELINE_FALLBACK_ALERTS <<- c(PIPELINE_FALLBACK_ALERTS,
                                 sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
  if (isTRUE(hard_stop)) {
    stop(msg, call. = FALSE)  # hard stop takes precedence - no point also raising a mere warning
  } else {
    warning(msg, call. = FALSE)  # ensures it also surfaces in R's end-of-run warning summary
  }
  invisible(NULL)
}

# --- Tiered bias correction (see Nechako SWE bias-correction framework) ---
# Extends the single basin-wide monthly ratio above into a hierarchy that
# mirrors Zahmatkesh et al. (2019), King et al. (2020), and Kanda & Fletcher
# (2025): Tier 1 = basin-wide mean/median bias subtraction (the block above,
# now tested in BOTH additive and multiplicative form); Tier 2 = regionalized
# CDF/quantile-matching using the SAME elevation x aspect regions built for
# SSPI pooling (build_regions()); Tier 3 = multiple linear regression using
# elevation, slope, aspect, latitude, longitude, and day-of-year (sin/cos) as
# predictors, following the Kanda & Fletcher (2025) predictor set. All three
# tiers are cross-validated with leave-one-station-out (LOSO) CV and the
# tier with the lowest CV RMSE is applied. A full (non-locally-trained)
# Random Forest tier is intentionally NOT included here: with only a
# handful of Nechako stations, training RF locally risks exactly the
# overfitting/poor-spatial-transfer behaviour both ML papers warn about.
# The framework's recommended path for an RF tier is to train it on the
# public CanSWE/NorSWE network for this ecozone and layer this script's
# local correction on top as a residual step - that training happens
# offline (outside this script) and the result would be supplied the same
# way STATION_OBS_PATH is supplied now.
#
# NOTE: this tiered framework (Tier 0 QA screen, Tier 1 additive/multiplicative
# with empirical-Bayes shrinkage, Tier 2 regionalized CDF-matching, Tier 3 MLR,
# LOSO-CV tier selection, and the 3-way validation protocol) is now WIRED INTO
# the "BASIN-WIDE BIAS CORRECTION" block below - it used to be defined but
# dead code (only Tier 1 ever ran). See rfa_diagnostics/ for the per-tier
# fit details and rfa_diagnostics/bias_correction_validation.csv for the
# random k-fold / spatial-block-LOSO / temporal-transfer scores.
BIAS_TIER_MODE   <- "auto"  # "tier1_only" (old behaviour), "tier2", "tier2_sar", "tier3", "auto" (CV-selects best)
MIN_STATION_MONTHS_TIER2 <- 8   # per region x month, mirrors MIN_RATIOS_PER_MONTH logic
MIN_STATION_MONTHS_TIER3 <- 15  # pooled across all regions/months for the MLR fit
LOSO_CV_MIN_STATIONS <- 4       # need at least this many distinct stations to run LOSO CV at all;
# below this, BIAS_TIER_MODE is forced to "tier1_only" since a meaningful
# held-out test isn't possible with fewer stations than that
#
# --- LOSO-CV held-out set size cap (station-count reduction) ---
# CONTEXT: the underlying prediction grid is only 790 valid 0.1 deg ERA5-Land
# pixels grouped into 17 homogeneous elevation x aspect regions
# (build_regions()); the station network merged in merge_station_networks()
# (ASWS + CanSWE) can nonetheless contain thousands of distinct station_id
# values (observed on the real Nechako run: 2958). Running true
# leave-one-STATION-out CV over every one of those IDs refits Tier 1/2/3 (each
# itself a nontrivial fit - Tier 2 alone refits 204 region-months) once per
# station, i.e. ~2958 full refits, which is not a station-count problem so
# much as a validation-DESIGN mismatch: many of those station IDs share the
# SAME 0.1 deg grid cell (see `cells(swe_template, pts)[, "cell"]` above), so
# holding out station A and station B in the same cell tests nothing new
# after the first - and the correction itself is only ~17 region-month
# curves (Tier 2) or one basin-wide MLR surface (Tier 3), not 2958 free
# parameters. The "spatial-block CV" logic Kanda & Fletcher (2025) recommend
# is about testing transfer to unseen LOCATIONS, and with 790 pixels / 17
# regions, a much smaller held-out set already tests that fully.
#
# FIX (see build_loso_holdout_stations() below): before LOSO-CV runs, the
# station list actually iterated is reduced by (1) collapsing to one
# representative station per OCCUPIED grid cell (co-located stations are
# redundant for spatial-transfer testing), then (2) stratified sampling
# across the 17 build_regions() regions so every region keeps at least
# min_per_region held-out sites, capped in total at LOSO_CV_MAX_STATIONS.
# Tier 1/2/3 themselves are still FIT on the FULL, un-subsampled station
# network in every fold (train_obs = station_obs[station_obs$station_id !=
# s, ]) - only the set of stations actually left OUT and scored is reduced.
# Set LOSO_CV_MAX_STATIONS <- Inf to restore the original exhaustive
# behaviour (every distinct station individually held out).
LOSO_CV_MAX_STATIONS   <- 40  # cap on distinct stations actually held out/scored in LOSO-CV
LOSO_CV_MIN_PER_REGION <- 2   # floor: try to keep at least this many held-out stations per region
# (so a small region is never silently unvalidated by the cap/subsampling)
#
# NOTE: these lines above are the ONLY tuning knobs that actually exist
# for tiered bias correction. There is no separate ENABLE_LOSO_CV flag - LOSO
# CV always runs automatically whenever BIAS_TIER_MODE == "auto" AND the
# station file has >= LOSO_CV_MIN_STATIONS distinct stations (see the
# "Only N distinct station(s)..." message in the BASIN-WIDE BIAS CORRECTION
# section if it doesn't). There is also no ENFORCE_MONOTONICITY flag -
# monotonicity rejection is unconditional inside fit_poly_qm() (Tier 2) and
# cannot be turned off. And Tier 2's per-region-month minimum is named
# MIN_STATION_MONTHS_TIER2 (matched station-MONTH pairs, not distinct
# stations) - dropping it to 3 would let very small, noisy samples through;
# lower it cautiously (e.g. to 5-6) only if you deliberately want more Tier-2
# coverage at the cost of shakier per-region fits, and check
# rfa_diagnostics/bias_correction_validation.csv afterwards to confirm it
# actually helped instead of just adding noise.

# --- Independent snow-cover-timing validation (OpenLandMap climatology) ---
# Cropped by Nechako_Snow_Data_download.R Section 4 to
# snow_data/openlandmap_snow_climatology/. Purely a VALIDATION input (it's
# snow-cover fraction, not SWE, so it never touches swe_matrix) - it checks
# whether the (bias-corrected) ERA5-Land-implied snow season timing lines up
# with a long, independent satellite record. See the bias-correction
# framework document, section 4.
SCF_CLIMATOLOGY_DIR <- "snow_data/openlandmap_snow_climatology"
SCF_CLIMATOLOGY_VALIDATE <- TRUE
SNOW_COVER_SWE_THRESHOLD_MM <- 4  # SWE >= this = "snow covered", per Bresson et al. (2026)

# --- Method-1 domain-mask fallback (when optical snow-cover is unavailable) ---
SWE_DOMAIN_NONZERO_FRACTION <- 0.75  # pixel-month is "in domain" if SWE is nonzero in >= 75% of reference years
SWE_DOMAIN_USE_FALLBACK    <- TRUE   # allow SWE-derived fallback for Method 1 when snow_cover_monthly.nc is missing

# --- Optional Tier 1.5 / Tier 2 SAR hooks (used only when files exist) ---
CANSWE_STATION_OBS_PATH <- "./snow_data/canswe/CanSWE_daily.csv"  # optional densification file; same canonical columns as daily.csv after standardization
SAR_DEPTH_PATH          <- "snow_data/sar_depth_monthly.nc"  # optional monthly SAR-depth / change-detection stack
USE_CANSWE_NETWORK      <- TRUE
USE_SAR_PROXY           <- TRUE


# ---- Reproducibility / Reference-period configuration ----
REFERENCE_PERIOD_MODE <- "fixed"   # "fixed" or "whole_record"
REFERENCE_START <- as.Date("1991-01-01")
REFERENCE_END   <- as.Date("2020-12-31")

MIN_PIXEL_GAMMA_SAMPLE    <- 15
MIN_REGIONAL_GAMMA_SAMPLE <- 25
MIN_BIAS_FORM_SELECTION_N <- 6
BASIN_SWEI_RANDOM_SEED    <- 40
set.seed(BASIN_SWEI_RANDOM_SEED)

# Parallel startup policy. On Windows, PSOCK uses localhost TCP sockets. If
# Rscript.exe is blocked by firewall/antivirus, cluster startup fails and the
# script falls back to serial execution. Set FALSE only after PSOCK is known to
# work on this machine. This avoids repeatedly attempting a broken backend.
ENABLE_PSOCK_PARALLEL <- TRUE
MAX_PSOCK_WORKERS <- 8L

# H0 Monte-Carlo uncertainty guard. A region-month is pooled only when the
# upper confidence bound is below the Hosking-Wallis H0=1 threshold.
H0_MC_CONF_LEVEL <- 0.95

# IMPORTANT METHODOLOGICAL NOTE:
# Hosking-Wallis RFA assumes independent sites. ERA5-Land pixels are spatially
# autocorrelated pseudo-sites, therefore Di/H0/KS diagnostics should be
# interpreted as descriptive rather than exact inferential statistics.
# Nearest-pixel station matching also introduces representativeness error in
# complex terrain because point observations and grid-cell means may differ.

cat(sprintf("\n- RFA regionalization for SSPI-1: %s\n",
            if (scf_method == "2" && REGIONALIZE_SSPI) "ENABLED" else "disabled"))

# NOTE: run_settings.txt is written further below (after "Input files"), not
# here - this used to be right at this point in the script, but it reads
# scf_file and effectively reports on DOMAIN_MASK_SOURCE, and NEITHER of
# those exists yet this early: scf_file isn't assigned until the "Input
# files" section below, and DOMAIN_MASK_SOURCE is only ever set (see the
# Method-1 SCF-fallback check right after scf_file) once we know whether the
# domain mask actually came from SCF or had to fall back to SWE. Writing the
# log here unconditionally threw "object 'scf_file' not found" (a hard
# Error, not a graceful fallback) on every run, before a single output file
# was produced - see run_settings.txt writeLines() call below for where this
# now actually happens.

# ---- Load Basin Boundary ----
basin_path <- "Spatial/Nechako_Basin.kmz"
basin <- if (file.exists(basin_path)) {
  tmp <- tempfile(); dir.create(tmp, showWarnings = FALSE)
  utils::unzip(basin_path, exdir = tmp)
  kml <- list.files(tmp, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)[1]
  v <- vect(kml); if (nrow(v) > 1L) v <- aggregate(v)
  unlink(tmp, recursive = TRUE); v
} else NULL
if (!is.null(basin)) {
  cat("[OK] Basin boundary loaded\n")
} else {
  stop("Basin boundary not found: ", basin_path)
}

# ---- Input files ----
swe_file <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
scf_file <- "monthly_data_direct/snow_cover_monthly.nc"
if (!file.exists(swe_file)) stop("SWE file not found: ", swe_file)

# DOMAIN_MASK_SOURCE records which domain-mask source Method 1 actually used
# this run (only meaningful for scf_method == "1" - Method 2/SSPI-1 doesn't
# use an SCF domain mask at all). Previously referenced in run_settings.txt
# via exists()-guarded NA fallback but never actually assigned anywhere in
# the script, so that line always silently printed NA even when a real
# fallback had occurred - defined here, at the one place that actually
# decides it, so the log means something.
if (scf_method == "1" && !file.exists(scf_file)) {
  DOMAIN_MASK_SOURCE <- "SWE-derived fallback (SCF file missing)"
  log_pipeline_fallback(
    sprintf("Method 1 selected, but %s is missing - deriving the domain mask directly from ERA5-Land SWE instead.", scf_file)
  )
} else if (scf_method == "1") {
  DOMAIN_MASK_SOURCE <- sprintf("SCF (%s)", scf_file)
} else {
  DOMAIN_MASK_SOURCE <- "not applicable (Method 2 / SSPI-1 does not use an SCF domain mask)"
}

# run_settings.txt is written HERE, not earlier in the script, because it
# reports on scf_file and DOMAIN_MASK_SOURCE, both of which are only just
# now defined (see the comment left at the old, premature call site above -
# writing it any earlier either hard-crashed on a not-yet-existing scf_file
# or silently reported DOMAIN_MASK_SOURCE as NA regardless of what actually
# happened).
writeLines(
  c(sprintf("run_timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("scf_method: %s (%s)", scf_method, scf_method_label),
    sprintf("DOMAIN_MASK_SOURCE: %s", DOMAIN_MASK_SOURCE),
    sprintf("SCF_FILE_FOUND (%s): %s", scf_file, file.exists(scf_file)),
    sprintf("REGIONALIZE_SSPI: %s", if (scf_method == "2") REGIONALIZE_SSPI else NA),
    sprintf("BIAS_CORRECT: %s", if (scf_method == "2") BIAS_CORRECT else NA),
    sprintf("station_file_found (%s): %s", STATION_OBS_PATH,
            if (scf_method == "2") file.exists(STATION_OBS_PATH) else NA),
    sprintf("canswe_file_found (%s): %s", CANSWE_STATION_OBS_PATH, file.exists(CANSWE_STATION_OBS_PATH)),
    sprintf("sar_proxy_found (%s): %s", SAR_DEPTH_PATH, file.exists(SAR_DEPTH_PATH)),
    sprintf("dem_file_found (%s): %s", DEM_PATH, file.exists(DEM_PATH)),
    sprintf("REQUIRE_DEM: %s", REQUIRE_DEM),
    sprintf("REQUIRE_STATION_OBS: %s", REQUIRE_STATION_OBS),
    "--- If a found line above is FALSE and the corresponding requirement flag is FALSE,",
    "    this run silently fell back to a lesser method - check",
    "    PIPELINE_FALLBACK_ALERTS.txt in this same folder for details. ---"),
  file.path(out_dir, "run_settings.txt")
)

cat("\n============================================================\n")
cat(if (scf_method == "2") "SSPI-1 CALCULATION (PARALLEL MODE)\n"
    else "SEASONAL SWEI CALCULATION (PARALLEL MODE)\n")
cat("============================================================\n\n")

cat("===== READING NETCDF FILES =====\n")
swe <- rast(swe_file)
if (scf_method == "1" && file.exists(scf_file)) {
  scf_stack <- rast(scf_file)
} else {
  scf_stack <- NULL  # not needed for SSPI-1; or unavailable, in which case Method 1 falls back to SWE-derived domain masking
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
    cat("[OK] Rotating SCF from 0-360 to -180-180 to match SWE...\n")
    scf_stack <- terra::rotate(scf_stack)
  } else if (swe_ext[2] > 180 && scf_ext[2] <= 180) {
    cat("[OK] Rotating SWE from 0-360 to -180-180 to match SCF...\n")
    swe <- terra::rotate(swe)
  }
}

# 2. Crop early to save memory and isolate the basin
if (!is.null(basin)) {
  cat("[OK] Cropping rasters to basin extent...\n")
  swe <- crop(swe, project(basin, crs(swe)))
  if (!is.null(scf_stack)) {
    scf_stack <- crop(scf_stack, project(basin, crs(scf_stack)))
  }
}

# 3. Resample to common grid
if (!is.null(scf_stack) && !compareGeom(scf_stack[[1]], swe[[1]], stopOnError = FALSE)) {
  cat("[OK] Resampling SWE to match SCF grid...\n")
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
cat(sprintf("[OK] Basin masking complete: %d pixels (%.1f%% of raster)\n",
            basin_pixels, 100 * basin_pixels / total_pixels))

# Safety Check against NaN calculation
if (is.na(basin_pixels) || basin_pixels == 0) {
  stop("FATAL ERROR: Basin masking resulted in 0 valid pixels. The SWE and SCF rasters may not overlap the basin boundary properly after resampling.")
}

# ---- Convert SWE to mm if needed ----
# Prefer NetCDF units; fall back to whole-stack maximum if unavailable.
#
# BUGFIX: swe is a multi-layer SpatRaster (one layer per month - 912 layers
# for a multi-decade monthly record). terra::global(swe, "max", ...) returns
# ONE ROW PER LAYER, not a single whole-stack value, despite the comment
# above calling it a "whole-stack maximum" - so swe_mean_sample was a
# 912-element vector, and `swe_mean_sample < 10` inside `&&` (which requires
# a length-1 logical) threw "'length = 912' in coercion to logical(1)".
# Collapsing with max() over the per-layer maxima gives the actual
# whole-stack maximum the unit-detection logic was always meant to check.
swe_layer_maxima <- global(swe, "max", na.rm = TRUE)$max
swe_mean_sample <- suppressWarnings(max(swe_layer_maxima, na.rm = TRUE))
if (is.finite(swe_mean_sample) && swe_mean_sample < 10) {
  swe <- swe * 1000
  cat("[OK] Converted SWE from meters to mm\n")
}

# ---- Convert SCF to fraction if needed ----
if (!is.null(scf_stack)) {
  # Same fix as above: scf_stack is also multi-layer, so global(..., "max")
  # returns one row per layer - collapse to a single whole-stack scalar
  # before the comparison.
  scf_layer_maxima <- global(scf_stack, "max", na.rm = TRUE)$max
  scf_max <- suppressWarnings(max(scf_layer_maxima, na.rm = TRUE))
  if (is.finite(scf_max) && scf_max > 1.5) {
    scf_stack <- scf_stack / 100
    cat("[OK] Converted SCF from percent to fraction\n")
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
cat(sprintf("\n- Proceeding with Method %s: %s\n", scf_method, scf_method_label))

cat("\n===== STEP 2: Building SCF mask =====\n")
target_k <- 3

if (scf_method == "2") {
  
  # ============================================================================
  #   BRANCH 2 - SSPI-1 (Zero-inflated Gamma, Monthly, No SCF Mask)
  # ============================================================================
  cat("\n- SSPI-1 selected: skipping SCF mask build (not applicable).\n")
  cat("  Zero-inflated gamma fitting will be applied per calendar month.\n")
  scf_mask_list  <- NULL
  basin_scf_mask <- NULL
  
  # Reference period for SSPI-1
  if (REFERENCE_PERIOD_MODE == "fixed") {
    sspi_ref_start <- REFERENCE_START
    sspi_ref_end   <- REFERENCE_END
  } else {
    sspi_ref_start <- min(dates)
    sspi_ref_end   <- max(dates)
  }
  ref_idx <- which(dates >= sspi_ref_start & dates <= sspi_ref_end)
  if (length(ref_idx) < 400) {
    warning(sprintf("SSPI-1 reference period has only %d months (recommended: >=480)", length(ref_idx)))
  }
  cat(sprintf("[OK] SSPI-1 reference period: %s to %s (%d months)\n",
              dates[ref_idx[1]], dates[ref_idx[length(ref_idx)]], length(ref_idx)))
  write(sprintf(
    "\nSSPI_REFERENCE_PERIOD_NOTE: %s to %s (n=%d) is both FIT to and SCORED against the same gamma distribution (standard SPI/SWEI practice, not a bug) - severity values inside this window are not out-of-sample. See drought_frequency_ref_vs_postref.csv if produced.",
    format(sspi_ref_start), format(sspi_ref_end), length(ref_idx)),
    file = file.path(out_dir, "run_settings.txt"), append = TRUE)
  
  # Zero-inflation screening
  cat("\n===== ZERO-INFLATION SCREENING (SSPI-1) =====\n")
  ref_swe_sspi   <- swe_matrix[, ref_idx, drop = FALSE]
  zero_prop_sspi <- rowMeans(ref_swe_sspi <= 0, na.rm = TRUE)
  invalid_pix_sspi <- which(zero_prop_sspi > 0.50)
  valid_pix_sspi   <- which(zero_prop_sspi <= 0.50)
  cat(sprintf("[OK] Valid pixels (<=50%% zero SWE in ref period): %d (%.1f%%)\n",
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
  #   BRANCH 1 - HUNING & AGHAKOUCHAK (2020) fixed 5% SCF threshold
  #            or SWE-derived domain mask fallback when optical SCF is absent
  # ============================================================================
  if (!is.null(scf_stack)) {
    domain_mask_source <- "SCF"
    scf_threshold <- 0.05
    cat(sprintf("[OK] H&A (2020): building SCF mask with fixed %.0f%% threshold\n",
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
    rm(scf_basin, scf_matrix_ha,v)
    
    cat("\n===== COMPUTING BASIN-LEVEL SCF MASK (H&A 2020) =====\n")
    basin_pixel_idx <- which(!is.na(swe_matrix[, 1]))
    basin_scf_mask  <- logical(12)
    for (m in 1:12) {
      mv_basin          <- scf_mask_list[[m]][basin_pixel_idx]
      basin_scf_mask[m] <- length(mv_basin) > 0 && mean(mv_basin, na.rm = TRUE) >= 0.5
    }
    cat(sprintf("[OK] Basin SCF mask: months included = %s\n",
                paste(month.abb[which(basin_scf_mask)], collapse = ", ")))
    
  } else {
    domain_mask_source <- "SWE_FALLBACK"
    if (isTRUE(SWE_DOMAIN_USE_FALLBACK)) {
      log_pipeline_fallback(sprintf(
        "Method 1 is running without snow_cover_monthly.nc - deriving the domain mask from ERA5-Land SWE persistence (SWE > 0 in >= %.0f%% of reference years, by calendar month).",
        100 * SWE_DOMAIN_NONZERO_FRACTION
      ))
      scf_mask_list <- vector("list", 12)
      basin_pixel_idx <- which(!is.na(swe_matrix[, 1]))
      for (m in 1:12) {
        cols_m <- which(months_all == m)
        if (length(cols_m) == 0) {
          scf_mask_list[[m]] <- rep(FALSE, ncell(swe))
          cat(sprintf("  Month %2d (%s): 0 pixels pass SWE-derived domain mask\n", m, month.abb[m]))
          next
        }
        swe_m <- swe_matrix[, cols_m, drop = FALSE]
        nonzero_year_fraction <- rowMeans(swe_m > 0, na.rm = TRUE)
        scf_mask_list[[m]] <- is.finite(nonzero_year_fraction) &
          (nonzero_year_fraction >= SWE_DOMAIN_NONZERO_FRACTION)
        cat(sprintf("  Month %2d (%s): %d pixels pass SWE-derived domain mask (>= %.0f%% nonzero years)\n",
                    m, month.abb[m],
                    sum(scf_mask_list[[m]], na.rm = TRUE),
                    100 * SWE_DOMAIN_NONZERO_FRACTION))
      }
      
      cat("\n===== COMPUTING BASIN-LEVEL DOMAIN MASK (SWE-derived fallback) =====\n")
      basin_scf_mask <- logical(12)
      for (m in 1:12) {
        mv_basin <- scf_mask_list[[m]][basin_pixel_idx]
        basin_scf_mask[m] <- length(mv_basin) > 0 && mean(mv_basin, na.rm = TRUE) >= 0.5
      }
      cat(sprintf("[OK] Basin domain mask (SWE fallback): months included = %s\n",
                  paste(month.abb[which(basin_scf_mask)], collapse = ", ")))
    } else {
      stop("Method 1 requested but no SCF file is available and SWE fallback is disabled.")
    }
  }
}

if (scf_method == "1" && !exists("domain_mask_source", inherits = FALSE)) {
  domain_mask_source <- if (!is.null(scf_stack)) "SCF" else "SWE_FALLBACK"
}


cat("\n===== COMPUTING BASIN-AVERAGED SWE =====\n")
if (scf_method == "1") {
  swe_basin_avg_raw <- colMeans(swe_matrix, na.rm = TRUE)
  cat(sprintf("[OK] Basin-averaged raw SWE: %d time steps\n", length(swe_basin_avg_raw)))
} else {
  swe_basin_avg_raw <- NULL
  cat("[OK] Basin-averaged SWE: skipped (not used by SSPI-1)\n")
}

swe_basin_avg_smoothed <- NULL
basin_avg_swei_results <- list()

# Before calling gringorten_swei_seasonal() for basin-average SWEI:
# old_seed <- if (exists(".Random.seed")) .Random.seed else NULL
# set.seed(BASIN_SWEI_RANDOM_SEED)
# ... gringorten_swei_seasonal(...)
# if (!is.null(old_seed)) .Random.seed <- old_seed


##############################################
# ---- STAGE 1 CHECKPOINT: save everything for script 02 ----
##############################################


# ---- Persist SpatRaster/SpatVector objects via terra's own I/O ----
# saveRDS() cannot round-trip these - their external C++ pointers do not
# survive serialization (see the "external pointer is not valid" crash this
# replaced). Write them with writeRaster()/writeVector() instead, and EXCLUDE
# them from the RDS payload below so nothing overwrites the rehydrated
# versions when script 02 unpacks the RDS.
terra::writeRaster(swe, file.path(STAGE_DIR, "stage1_swe.tif"), overwrite = TRUE)
if (!is.null(scf_stack)) {
  terra::writeRaster(scf_stack, file.path(STAGE_DIR, "stage1_scf_stack.tif"), overwrite = TRUE)
}
if (!is.null(basin)) {
  terra::writeVector(basin, file.path(STAGE_DIR, "stage1_basin.gpkg"), overwrite = TRUE)
}

# Defensive terra-pointer sweep: never let any Spat* object enter saveRDS().
# This catches future temporary rasters/vectors as well as the known core
# objects. Recursive (not just a top-level inherits() check) so a Spat*
# object buried inside a list or data.frame column - which the old
# top-level-only check would silently miss and let flow straight into
# saveRDS() - gets caught too.
has_nested_spat <- function(x, depth = 0L, max_depth = 6L) {
  if (depth > max_depth) return(FALSE)
  if (inherits(x, c("SpatRaster", "SpatVector",
                    "SpatRasterDataset", "SpatRasterCollection"))) return(TRUE)
  if (is.list(x)) return(any(vapply(x, has_nested_spat, logical(1),
                                    depth = depth + 1L, max_depth = max_depth)))
  FALSE
}
all_stage1_names <- ls(all.names = TRUE, envir = environment())
spat_leftover <- Filter(function(nm) {
  obj <- tryCatch(get(nm, envir = environment(), inherits = FALSE), error = function(e) NULL)
  !is.null(obj) && has_nested_spat(obj)
}, all_stage1_names)
if (length(spat_leftover) > 0L) {
  cat(sprintf("[CHECKPOINT] Excluding %d live terra Spat* object(s) from stage1.rds: %s\n",
              length(spat_leftover), paste(spat_leftover, collapse = ", ")))
}

stage1_objects <- setdiff(all_stage1_names,
                          unique(c("STAGE_DIR", "stage1_path", "has_nested_spat",
                                   "all_stage1_names", "spat_leftover", spat_leftover)))
stage1_path <- file.path(STAGE_DIR, "stage1.rds")
saveRDS(mget(stage1_objects, envir = environment(), inherits = FALSE), stage1_path, compress = "xz")

# Verify before releasing terra pointers, then remove all live Spat* objects so
# RStudio Environment refresh / workspace autosave cannot inspect stale C++ pointers.
stage1_verify <- readRDS(stage1_path)
if (!is.list(stage1_verify) || length(stage1_verify) != length(stage1_objects)) {
  stop("stage1.rds verification failed: object count/type mismatch after saveRDS().", call. = FALSE)
}
rm(stage1_verify)
if (length(spat_leftover) > 0L) rm(list = spat_leftover, envir = environment())
gc(verbose = FALSE)
gc(verbose = FALSE)  # second pass: catches wrapper objects whose finalizer
# only actually freed the underlying pointer during the
# first pass, before RStudio's Environment pane refresh
# or .RData autosave gets a chance to touch them

cat(sprintf("[OK] Saved %d object(s) to: %s\n", length(stage1_objects), normalizePath(stage1_path)))
cat(sprintf("[OK] File size: %.1f MB\n", file.size(stage1_path) / 1024^2))
cat("\nNext: run 02_regionalization_bias.R\n")