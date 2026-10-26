####################################################################################
# 12basin_timeseries.R  ·  BASIN-AVERAGED TIME SERIES — OBSERVED FRAMEWORK
# ─────────────────────────────────────────────────────────────────────────────────
# SCOPE  (authoritative home for all basin-averaged outputs):
#   PART 1 – Load basin-averaged CSVs produced by 11trend_test_ALL.R for
#             SPI, SPEI-PM, SPEI-ERA5, SWEI, SSPI and optionally MSPI / MSPEI.
#             SPEI-PM and SPEI-ERA5 are parallel observed formulations;
#             no detrended/counterfactual branch is used.
#   PART 2 – Write per-index flat CSVs to basin_averaged_timeseries/
#             (re-exports the authoritative *_basin_averaged_by_month.csv files
#             into a plain date/value layout for downstream use).
#   PART 3 – Combined CSV (basin_mean_all_indices.csv — all indices long format).
#   PART 4 – Time series plots (one publication-quality PNG per index × scale).
#   PART 4d– Seasonality bar chart: mean monthly P vs PET(PM) (12 values each),
#             1950–2025 average; climatic water balance panel.
#             Requires ERA5Land_Nechako_PET_monthly_summary.csv from
#             2b_PET_ERALand.R.  Output: FigS_Seasonality_P_PET.pdf/.png.
#   PART 4e– Full historical 18-panel chronology + recent formulation comparison.
#   PART 4e-C2– Manuscript 12-panel SPI + SPEI-PM chronology across all six scales.
#   PART 4e-D– Supplementary Fig. S1: 1- and 6-month historical chronology.
#             Fig2_MS_AllIndices_3x6.pdf/.png (18-panel: SPI, SPEI-PM, SPEI-ERA5 × 1, 3, 6, 12, 24, 36 months).
#             Fig2S_MS_Recent_AllScales_3x2.pdf/.png (2022–2025 SPI + SPEI-PM comparison across all six scales).
#             The historical figure is the primary Results chronology for Section 4.2.
#             Section 4.2 of the manuscript.
#   PART 5 – Excel summary: descriptive statistics + drought event counts
#             (2020-2025) + pairwise correlation matrix (Drought_Summary.xlsx).
#   PART 6 – Optional support-aware exports for station SSI and surface-water
#             TSDI / NTSDI. These are kept separate from basin-mean indices.
#
# PREREQUISITE: run 11trend_test_ALL.R first.
# NOTE: 13trends_visualization.R handles ONLY spatial figures (Parts 1-2 of 8).
#       All basin-averaged time-series outputs live here.
#
# METHODOLOGY NOTE:
#   Step 1 — SPI/SPEI are computed per-pixel by 1/3SPI/SPEI_ERALand.R /
#             3SPEI_ERALand.R; results live in seasonal NetCDF files.
#   Step 2 — 11trend_test_ALL.R area-weights per-pixel values with
#             terra::cellSize() on the common ERA5-Land grid and writes
#             {index}_{scale:02d}_basin_averaged_by_month.csv.
#   Step 3 — This script reads those authoritative CSVs.
#             Do NOT recompute averages here.
#
# ─────────────────────────────────────────────────────────────────────────────────
# FIX LOG (this revision):
#   The SPEI-ERA5 branches of load_basin_index() and load_any_ts() referenced an
#   undefined constant SPEI_ERA5_BASIN_DIR and parsed a wide-format
#   "spei_era5_%02d_basin_average.csv" file. That code path is inconsistent with
#   this script's own documented contract above (SPEI-ERA5 comes from
#   11trend_test_ALL.R via the same *_basin_averaged_by_month.csv convention as
#   SPI / SPEI-PM) and is almost certainly leftover from before the pipeline was
#   consolidated. Both functions now route "spei_era5" through the same
#   load_basin_avg_csv() call used for "spi" and "spei", removing the dependency
#   on the missing constant.
#
#   If load_basin_avg_csv("spei_era5", sc) throws a file-not-found error when you
#   re-run this, it means SPEI-ERA5 genuinely lives in a separate wide-format
#   file/directory rather than the unified convention — in that case the old
#   branch needs to be restored, but pointed at the *real* directory (please
#   supply that path).
####################################################################################
# WD_PATH / UTILS_FILE are set from getwd() at the top of this script.
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

source(UTILS_FILE)

# Explicit scales shared with Script 11; scale 2 is not a manuscript output.
SPI_SCALES <- c(1L, 3L, 6L, 12L, 24L, 36L)
SPEI_SCALES <- c(1L, 3L, 6L, 12L, 24L, 36L)
SPEI_ERA5_SCALES <- c(1L, 3L, 6L, 12L, 24L, 36L)
MANUSCRIPT_SCALES <- SPI_SCALES
SWEI_SCALE <- 3L
SSPI_SCALE <- 1L
SWEI_SEAS_DIR <- file.path(WD_PATH, "swei_results_seasonal")
SSPI_SEAS_DIR <- file.path(WD_PATH, "sspi_results_seasonal")
DROUGHT_ONSET <- if (exists("DROUGHT_ONSET", inherits = TRUE)) DROUGHT_ONSET else -0.5

utils_load_packages(c("terra", "ggplot2", "lubridate", "dplyr",
                      "data.table", "patchwork", "openxlsx", "zoo",
                      "RColorBrewer", "scales", "sf"))


## ── Output directories ──────────────────────────────────────────────────────────
BASIN_PLOT_DIR <- file.path(WD_PATH, "temporal_drought", "basin_averaged_plots")
BASIN_TS_DIR   <- file.path(WD_PATH, "temporal_drought", "basin_averaged_timeseries")
for (d in c(BASIN_PLOT_DIR, BASIN_TS_DIR, CACHE_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# NATIVE ERA5-LAND GRID CONTRACT
#
# Script 12 uses the original/native ERA5-Land precipitation grid.
#
# The first layer of:
#   monthly_data_direct/total_precipitation_monthly.nc
#
# is the sole raster template. No canonical/project-wide analysis grid is
# loaded, generated, or required.
#
# Raster products are never reprojected/resampled merely to satisfy Script 12.
# Only vector geometries may be transformed to the native raster CRS.
# ==============================================================================

native_precip_file <- file.path(
  WD_PATH,
  "monthly_data_direct",
  "total_precipitation_monthly.nc"
)

if (!file.exists(native_precip_file)) {
  stop(
    "Native ERA5-Land precipitation NetCDF not found: ",
    native_precip_file,
    "\nRun Monthly_ERALand_download.R first."
  )
}

native_precip <- terra::rast(native_precip_file)

if (!inherits(native_precip, "SpatRaster")) {
  stop(
    "Failed to load native ERA5-Land precipitation NetCDF as a SpatRaster."
  )
}

if (terra::nlyr(native_precip) < 1L) {
  stop(
    "Native ERA5-Land precipitation NetCDF contains no raster layers."
  )
}

# Sole authoritative spatial template.
analysis_template <- native_precip[[1]]

if (is.na(terra::crs(analysis_template)) ||
    !nzchar(terra::crs(analysis_template))) {
  stop(
    "Native ERA5-Land precipitation grid has no CRS."
  )
}

cat(sprintf(
  "[GRID] Native ERA5-Land precipitation grid: %d x %d = %d cells | ",
  terra::nrow(analysis_template),
  terra::ncol(analysis_template),
  terra::ncell(analysis_template)
))

cat(sprintf(
  "resolution %.10f x %.10f | CRS=%s\n",
  terra::res(analysis_template)[1],
  terra::res(analysis_template)[2],
  terra::crs(analysis_template)
))

# ------------------------------------------------------------------------------
# Native-grid geometry validator
# ------------------------------------------------------------------------------
assert_native_grid <- function(x, template, label) {
  
  if (!inherits(x, "SpatRaster")) {
    stop(label, " is not a SpatRaster.")
  }
  
  if (terra::nlyr(x) < 1L) {
    stop(label, " contains no raster layers.")
  }
  
  geom_ok <- terra::compareGeom(
    x[[1]],
    template,
    stopOnError = FALSE,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE
  )
  
  orig_ok <- isTRUE(all.equal(
    terra::origin(x[[1]]),
    terra::origin(template),
    tolerance = 1e-6
  ))
  
  if (!geom_ok || !orig_ok) {
    stop(
      label,
      " does not match the native ERA5-Land precipitation grid.\n",
      "Script 12 does not reproject or resample raster inputs.\n",
      "Check the upstream ERA5-Land NetCDF files."
    )
  }
  
  invisible(TRUE)
}

## ── Basin geometry and fractional-coverage weights ──────────────────────────────
## Basin = dissolved union of the 7 subbasins (as in 3SPEI_ERALand.R). Any basin
## mean computed in this script uses cell area x fractional coverage and the
## shared BASIN_MIN_VALID_AREA_COVERAGE (0.20) reporting rule.
basin_shp <- nechako_subbasin_union(load_nechako_subbasins())
if (!terra::same.crs(basin_shp, analysis_template)) {
  basin_shp <- terra::project(
    basin_shp,
    terra::crs(analysis_template)
  )
}
S12_BASIN_BW <- basin_fractional_weights(analysis_template, basin_shp)
cat(sprintf(
  "✓ Basin: %.1f km² (polygon) | %d cells, %d partial | weighted %.1f km² | min valid coverage %.2f\n",
  sum(terra::expanse(basin_shp, unit = "km")),
  S12_BASIN_BW$n_cells, S12_BASIN_BW$n_partial, S12_BASIN_BW$area_m2 / 1e6,
  BASIN_MIN_VALID_AREA_COVERAGE
))

## Fractional-coverage basin mean of a CSV pixel matrix whose rows carry lon/lat.
s12_basin_mean_xy <- function(mat, lon, lat) {
  w <- basin_weights_for_xy(lon, lat, analysis_template, S12_BASIN_BW)
  regional_weighted_mean(as.matrix(mat), w,
                         min_coverage = BASIN_MIN_VALID_AREA_COVERAGE,
                         total_weight = S12_BASIN_BW$area_m2)
}

## ── MSPI / MSPEI single scale (mirror of 8 constant) ───────────────────────────
MSPI_MSPEI_SCALE <- 1L
if (!"spei_era5" %in% names(index_colours)) index_colours["spei_era5"] <- "#d62728"
index_colours["swei"] <- "#1E90FF"
index_colours["sspi"] <- "#9467bd"

################################################################################
# HELPERS
################################################################################

# Base-R month arithmetic avoids the non-exported lubridate::months symbol.
# Inputs/outputs are Date values representing the first day of a month.
next_month_start <- function(x) {
  x <- as.Date(x)
  if (!length(x) || all(is.na(x))) return(as.Date(NA))
  month_start <- as.Date(format(max(x, na.rm = TRUE), "%Y-%m-01"))
  seq.Date(month_start, by = "month", length.out = 2L)[2L]
}

month_lengths <- function(x) {
  month_start <- as.Date(format(as.Date(x), "%Y-%m-01"))
  next_start <- as.Date(format(month_start + 32, "%Y-%m-01"))
  as.integer(next_start - month_start)
}

load_supported_basin_csv <- function(index_type, scale) {
  idx <- tolower(index_type)
  if (idx == "sspi") {
    return(load_basin_avg_csv(idx, scale, data_dir = SSPI_SEAS_DIR))
  }
  load_basin_avg_csv(idx, scale)
}

## Load the authoritative basin-average CSV and return a data.frame
## (Date, Index, Mean_Value, Drought_Fraction) ready for plotting.
##
## FIXED: spi / spei / spei_era5 are all produced by 11trend_test_ALL.R in the
## same *_basin_averaged_by_month.csv convention, so all three now go through
## the same load_basin_avg_csv() call. The old spei_era5-specific branch that
## referenced the undefined SPEI_ERA5_BASIN_DIR constant has been removed.
load_basin_index <- function(index_type, scale) {
  label <- sprintf("%s%d", toupper(index_type), as.integer(scale))
  
  df <- tryCatch(
    load_supported_basin_csv(index_type, as.integer(scale)),
    error = function(e) {
      stop(sprintf("%s: %s", label, conditionMessage(e)), call. = FALSE)
    }
  )
  if (is.null(df) || !nrow(df))
    stop("Required basin-average product is empty for ", label, call. = FALSE)
  
  data.frame(
    Date             = df$date,
    Index            = label,
    Mean_Value       = as.numeric(df$value),
    Drought_Fraction = as.numeric(as.numeric(df$value) < DROUGHT_ONSET),
    stringsAsFactors = FALSE
  )
}

## Load MSPI / MSPEI time series from the Excel file written by their own scripts.
## Returns data.frame(date, value) or NULL if the file is absent.
load_mspi_mspei_ts <- function(index_type) {
  results_dir <- file.path(WD_PATH,
                           if (index_type == "mspi") "mspi_results" else "mspei_results")
  pat  <- sprintf("%s_basin_timeseries.*\\.xlsx$", index_type)
  fls  <- list.files(results_dir, pattern = pat, full.names = TRUE)
  if (!length(fls)) {
    cat(sprintf("  ⚠ No Excel timeseries found for %s in %s\n",
                toupper(index_type), results_dir))
    return(NULL)
  }
  tryCatch({
    df       <- openxlsx::read.xlsx(fls[1], sheet = "Basin_Timeseries")
    df$date  <- as.Date(paste0(df$Date, "-01"))
    df$value <- df$Basin_Mean
    df[, c("date", "value")]
  }, error = function(e) {
    cat(sprintf("  ⚠ Could not read %s Excel: %s\n", toupper(index_type), e$message))
    NULL
  })
}

## Unified loader — returns data.frame(date, value) for any supported index.
## Used by the CSV writer, plot builder, and Excel export.
##
## FIXED: same change as load_basin_index() — spei_era5 now goes through
## load_basin_avg_csv() instead of a separate wide-format CSV parser that
## depended on the undefined SPEI_ERA5_BASIN_DIR constant.
load_any_ts <- function(idx, sc) {
  if (idx %in% c("mspi", "mspei")) {
    df <- load_mspi_mspei_ts(idx)
    return(if (!is.null(df)) as.data.frame(df) else NULL)
  }
  df <- tryCatch(load_supported_basin_csv(idx, sc), error = function(e) NULL)
  if (!is.null(df)) as.data.frame(df) else NULL
}

################################################################################
# PART 1 – LOAD ALL INDICES
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 1: Loading basin-averaged indices\n")
cat("══════════════════════════════════════\n")

all_data <- list()

cat("\n── SPI ──\n")
for (sc in SPI_SCALES) {
  ts <- tryCatch(load_basin_index("spi", sc), error = function(e) {
    cat(sprintf("  ⚠ SPI-%02d unavailable: %s\n", sc, e$message)); NULL
  })
  if (!is.null(ts)) {
    all_data[[sprintf("SPI%d", sc)]] <- ts
    cat(sprintf("  ✓ SPI%d: %d time steps\n", sc, nrow(ts)))
  }
}

cat("\n── SPEI-PM (Penman-Monteith) ──\n")
for (sc in SPEI_SCALES) {
  ts <- tryCatch(load_basin_index("spei", sc), error = function(e) {
    cat(sprintf("  ⚠ SPEI-PM-%02d unavailable: %s\n", sc, e$message)); NULL
  })
  if (!is.null(ts)) {
    ts$Index <- sprintf("SPEI_PM%d", sc)
    all_data[[sprintf("SPEI_PM%d", sc)]] <- ts
    cat(sprintf("  ✓ SPEI-PM-%d: %d time steps\n", sc, nrow(ts)))
  }
}

cat("\n── SPEI-ERA5 (ERA5-Land native PEV) ──\n")
for (sc in SPEI_ERA5_SCALES) {
  ts <- tryCatch(load_basin_index("spei_era5", sc), error = function(e) {
    cat(sprintf("  ⚠ SPEI-ERA5-%02d unavailable: %s\n", sc, e$message)); NULL
  })
  if (!is.null(ts)) {
    ts$Index <- sprintf("SPEI_ERA5%d", sc)
    all_data[[sprintf("SPEI_ERA5%d", sc)]] <- ts
    cat(sprintf("  ✓ SPEI-ERA5-%d: %d time steps\n", sc, nrow(ts)))
  }
}

cat("\n── SWEI (snow water equivalent drought; scale 3) ──\n")
ts <- tryCatch(load_basin_index("swei", SWEI_SCALE), error = function(e) {
  cat(sprintf("  ⚠ SWEI-%02d unavailable: %s\n", SWEI_SCALE, e$message)); NULL
})
if (!is.null(ts)) {
  ts$Index <- sprintf("SWEI%d", SWEI_SCALE)
  all_data[[sprintf("SWEI%d", SWEI_SCALE)]] <- ts
  cat(sprintf("  ✓ SWEI-%02d: %d time steps\n", SWEI_SCALE, nrow(ts)))
}

cat("\n── SSPI (snowpack drought; scale 1) ──\n")
ts <- tryCatch(load_basin_index("sspi", SSPI_SCALE), error = function(e) {
  cat(sprintf("  ⚠ SSPI-%02d unavailable: %s\n", SSPI_SCALE, e$message)); NULL
})
if (!is.null(ts)) {
  ts$Index <- sprintf("SSPI%d", SSPI_SCALE)
  all_data[[sprintf("SSPI%d", SSPI_SCALE)]] <- ts
  cat(sprintf("  ✓ SSPI-%02d: %d time steps\n", SSPI_SCALE, nrow(ts)))
}

cat("\n── MSPI / MSPEI (optional) ──\n")
for (multi_idx in c("mspi", "mspei")) {
  df_multi <- load_mspi_mspei_ts(multi_idx)
  if (!is.null(df_multi) && nrow(df_multi)) {
    label    <- sprintf("%s%d", toupper(multi_idx), MSPI_MSPEI_SCALE)
    ts_multi <- data.frame(
      Date             = df_multi$date,
      Index            = label,
      Mean_Value       = df_multi$value,
      Drought_Fraction = as.numeric(df_multi$value < DROUGHT_ONSET),
      stringsAsFactors = FALSE)
    all_data[[label]] <- ts_multi
    cat(sprintf("  ✓ %s: %d time steps\n", label, nrow(ts_multi)))
  }
}

if (!length(all_data))
  stop("No basin-averaged drought-index data loaded. Run Script 11 (11trend_test_ALL.R) and Script 4 (4snow_drought_indices_ERALand.R) first, and verify their output paths.")

################################################################################
# PART 2 – PER-INDEX FLAT CSVs  (date / value layout)
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 2: Writing per-index flat CSVs\n")
cat("══════════════════════════════════════════════\n")

write_basin_ts_csv <- function(idx, sc) {
  df <- load_any_ts(idx, sc)
  if (is.null(df) || !nrow(df)) {
    cat(sprintf("  ⚠ %s-%02d: no data available\n", toupper(idx), sc))
    return(invisible(NULL))
  }
  out <- file.path(BASIN_TS_DIR,
                   sprintf("%s_%02d_basin_average.csv", idx, sc))
  utils::write.csv(df, out, row.names = FALSE)
  cat(sprintf("  ✓ %s-%02d: %d months → %s\n",
              toupper(idx), sc, nrow(df), basename(out)))
  invisible(df)
}

cat("\n  SPI, SPEI-PM, SPEI-ERA5, SWEI & SSPI...\n")
for (idx in c("spi", "spei", "spei_era5")) {
  scales <- if (idx == "spi") SPI_SCALES else if (idx == "spei_era5") SPEI_ERA5_SCALES else SPEI_SCALES
  for (sc in scales) write_basin_ts_csv(idx, sc)
}
write_basin_ts_csv("swei", SWEI_SCALE)
write_basin_ts_csv("sspi", SSPI_SCALE)

cat("\n  MSPI & MSPEI (if available)...\n")
for (idx in c("mspi", "mspei")) write_basin_ts_csv(idx, MSPI_MSPEI_SCALE)

################################################################################
# PART 3 – COMBINED LONG-FORMAT CSV
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 3: Combined long-format CSV\n")
cat("══════════════════════════════════════════════\n")

all_basin_data <- dplyr::bind_rows(all_data)
combined_csv   <- file.path(BASIN_PLOT_DIR, "basin_mean_all_indices.csv")
write.csv(all_basin_data, combined_csv, row.names = FALSE)
cat(sprintf("  ✓ %d rows → %s\n", nrow(all_basin_data), basename(combined_csv)))

################################################################################
# PART 4 – TIME SERIES PLOTS  (one PNG per index × scale)
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 4: Time series plots\n")
cat("══════════════════════════════════════════════\n")

## Helper: read back a per-index flat CSV written in Part 2.
load_flat_csv <- function(idx, sc) {
  f <- file.path(BASIN_TS_DIR,
                 sprintf("%s_%02d_basin_average.csv", idx, sc))
  if (!file.exists(f)) stop("CSV not found: ", f)
  df       <- data.table::fread(f)
  df$date  <- as.Date(df$date)
  df$value <- as.numeric(df$value)
  as.data.frame(df)
}

## Publication-quality basin-average time series plot.
## Drought-fill ribbon, categorical band shading, event count subtitle.
make_basin_ts_plot <- function(df, idx, sc) {
  idx_lc    <- tolower(idx)
  idx_label <- if (idx_lc == "spei") "SPEI-PM" else
    if (idx_lc == "spei_era5") "SPEI-ERA5" else toupper(idx)
  idx_color <- index_colours[idx_lc]
  if (is.na(idx_color)) idx_color <- "steelblue"
  
  ## Drought events: shared framework definition (DROUGHT_ONSET / DROUGHT_END,
  ## scale-dependent minimum duration). Shading marks event months only.
  min_dur <- get_event_min_duration(sc)
  ev_tab <- tryCatch(
    identify_drought_events_ts(df$date, df$value, scale = sc),
    error = function(e) NULL)
  n_ev <- if (is.null(ev_tab)) NA_integer_ else nrow(ev_tab)
  df$in_event <- if (is.null(ev_tab)) FALSE else attr(ev_tab, "in_event")
  
  ## Dynamic y-axis: expand 10% beyond the observed min/max,
  ## but never narrower than ±3.5 (the standard SPI display range).
  y_min  <- min(df$value, na.rm = TRUE)
  y_max  <- max(df$value, na.rm = TRUE)
  y_lo   <- min(-3.5, floor(y_min * 1.10))   # at least -3.5; expand below if needed
  y_hi   <- max( 3.5, ceiling(y_max * 1.10)) # at least +3.5; expand above if needed
  
  ggplot2::ggplot(df, ggplot2::aes(date, value)) +
    drought_band_layers() +
    ggplot2::geom_hline(yintercept = 0,
                        color = "black", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1, 1),
                        linetype = "dashed", color = "gray50", linewidth = 0.3) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = ifelse(in_event, value, 0),
        ymax = 0),
      fill = "#c0392b", alpha = 0.35) +
    ggplot2::geom_line(color = idx_color, linewidth = 0.6,na.rm=TRUE) +
    shared_ts_theme(14) +
    ggplot2::scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2025-12-31"), by = "5 years"), date_labels = "%Y") +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s-%02d  Basin-Averaged Time Series (Area-Weighted)",
                         idx_label, sc),
      subtitle = sprintf("%d drought events  (onset < %.2f, end >= %.2f, min %d months)",
                         n_ev, DROUGHT_ONSET, DROUGHT_END, min_dur),
      x = "Year",
      y = sprintf("%s-%02d Index Value", idx_label, sc))
}

cat("\n  SPI, SPEI, SWEI & SSPI...\n")
for (idx in c("spi", "spei", "spei_era5", "swei", "sspi")) {
  scales <- switch(idx,
                   spi = SPI_SCALES,
                   spei = SPEI_SCALES,
                   spei_era5 = SPEI_ERA5_SCALES,
                   swei = SWEI_SCALE,
                   sspi = SSPI_SCALE
  )
  for (sc in scales) {
    tryCatch({
      df <- load_flat_csv(idx, sc)
      out <- file.path(BASIN_PLOT_DIR, sprintf("%s_%02d_timeseries.png", idx, sc))
      ggplot2::ggsave(out, make_basin_ts_plot(df, idx, sc),
                      width = 12, height = 6, dpi = 150)
      cat(sprintf("  ✓ %s-%02d\n", toupper(idx), sc))
    }, error = function(e) {
      cat(sprintf("  ⚠ %s-%02d: %s\n", toupper(idx), sc, e$message))
    })
  }
}

cat("\n  MSPI & MSPEI (if available)...\n")
for (idx in c("mspi", "mspei")) {
  tryCatch({
    df  <- load_flat_csv(idx, MSPI_MSPEI_SCALE)
    out <- file.path(BASIN_PLOT_DIR,
                     sprintf("%s_%02d_timeseries.png", idx, MSPI_MSPEI_SCALE))
    ggplot2::ggsave(out, make_basin_ts_plot(df, idx, MSPI_MSPEI_SCALE),
                    width = 12, height = 6, dpi = 150)
    cat(sprintf("  ✓ %s-%02d\n", toupper(idx), MSPI_MSPEI_SCALE))
  }, error = function(e)
    cat(sprintf("  ⚠ %s-%02d: %s\n", toupper(idx), MSPI_MSPEI_SCALE, e$message)))
}

################################################################################
# PART 4b – MANUSCRIPT FIGURE 4  (two PET-method versions + SPEI-2 fix)
#
# ROOT CAUSE OF SPEI-2 MISSING:
#   SPEI_SCALES in DROUGHT_ANALYSIS_utils.R = c(1, 2, 3, 6, 12).  Scale 2 is
#   never processed by 11trend_test_ALL.R, so no basin-averaged CSV exists.
#   However, 3SPEI_ERALand.R writes per-calendar-month CSV files for every
#   scale it computes (including scale 2) inside spei_results_seasonal/ and
#   spei_results_seasonal_era5land/.  load_spei_seasonal_csvs() reads those files
#   directly, computing an area-weighted basin mean on-the-fly.  This bypasses
#
# TWO PET METHOD VERSIONS:
#   Version A — SPEI_PM  : Penman-Monteith PET  (spei_results_seasonal/)
#   Version B — SPEI_ERA5 : ERA5-Land PEV PET      (spei_results_seasonal_era5land/)
#   Each version produces its own PDF + PNG.  An overlay version (both on
#   the same panels) is also saved.
#
# Outputs (all in BASIN_PLOT_DIR):
#   Fig4_MS_SPEI_123_PM_timeseries.pdf/.png       — PM only (3 panels)
#   Fig4_MS_SPEI_123_ERA5_timeseries.pdf/.png      — ERA5-Land PEV only
#   Fig4_MS_SPEI_123_overlay_timeseries.pdf/.png  — both overlaid per panel
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 4b: Manuscript Figure 4 — SPEI-1/2/3 time series (PM + ERA5-Land PEV)\n")
cat("══════════════════════════════════════════════\n")

# ── Directory paths for the two SPEI variants ──────────────────────────────────
SPEI_PM_DIR  <- file.path(WD_PATH, "spei_results_seasonal")
SPEI_ERA5_DIR <- file.path(WD_PATH, "spei_results_seasonal_era5land")
SPI_SEAS_DIR <- file.path(WD_PATH, "spi_results_seasonal")   # used by Part 4c + 4e
SPEI_ERA5_BASIN_DIR <- file.path(WD_PATH, "spei_results_seasonal_era5land")
MONTH_NAMES_3 <- c("Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec")

# ── Direct per-month CSV loader (bypasses 6 for any scale) ───────────────────
# Reads the 12 per-calendar-month CSV files written by 3SPEI_ERALand.R and
# computes the fractional-coverage area-weighted basin mean (0.20 coverage rule).
# Returns data.frame(date, value) in chronological order, or NULL on failure.
load_spei_seasonal_csvs <- function(scale, seas_dir) {
  if (!dir.exists(seas_dir)) return(NULL)
  
  month_rows <- vector("list", 12L)
  
  for (m in seq_len(12L)) {
    csv_file <- file.path(
      seas_dir,
      sprintf("spei_%02d_month%02d_%s.csv", scale, m, MONTH_NAMES_3[m])
    )
    if (!file.exists(csv_file)) next
    
    df <- tryCatch(data.table::fread(csv_file, data.table = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    
    yr_cols <- setdiff(names(df), c("lon", "lat"))
    yr_cols <- yr_cols[grepl("^[0-9]{4}$", yr_cols)]
    if (!length(yr_cols)) next
    
    mat_yr <- as.matrix(df[, yr_cols, drop = FALSE])
    vals <- s12_basin_mean_xy(mat_yr, df$lon, df$lat)
    
    vals[is.nan(vals)] <- NA_real_
    
    month_rows[[m]] <- data.frame(
      date = as.Date(sprintf("%s-%02d-01", yr_cols, m)),
      value = as.numeric(vals),
      stringsAsFactors = FALSE
    )
  }
  
  valid <- Filter(Negate(is.null), month_rows)
  if (!length(valid)) return(NULL)
  
  out <- do.call(rbind, valid)
  out <- out[order(out$date), ]
  out <- out[is.finite(out$value) & !is.na(out$value), ]
  rownames(out) <- NULL
  out
}


# ── Shared panel constants ──────────────────────────────────────────────────────
THR_EVENT       <- -0.5
THR_MOD         <- -1.0
HIGHLIGHT_START <- as.Date("2022-01-01")
HIGHLIGHT_END   <- as.Date("2025-12-31")

FILL_WET  <- "#4393c3"
FILL_MILD <- "#ffffbf"
FILL_MOD  <- "#fdae61"
FILL_SEV  <- "#f46d43"
FILL_EXT  <- "#d73027"

fill_pal <- c(
  "Wet"              = FILL_WET,
  "Mild deficit"     = FILL_MILD,
  "Moderate drought" = FILL_MOD,
  "Severe drought"   = FILL_SEV,
  "Extreme drought"  = FILL_EXT
)

# ── Single-version panel builder ───────────────────────────────────────────────
# pet_type : "PM" or "ERA5" — controls which directory and line colour
# is_first : TRUE = show legend; FALSE = suppress
make_spei_panel_v2 <- function(sc, panel_lab, pet_type = "PM",
                               show_xlab = FALSE, is_first = FALSE) {
  seas_dir  <- if (pet_type == "PM") SPEI_PM_DIR else SPEI_ERA5_DIR
  col_line  <- if (pet_type == "PM") "#1f77b4" else "#d62728"
  pet_label <- if (pet_type == "PM") "SPEI-PM (Penman-Monteith)"
  else                  "SPEI-ERA5 (ERA5-Land PEV)"
  
  # Prefer authoritative basin-average CSVs from 11trend_test_ALL.R.
  # Fall back to per-month seasonal CSVs only when an authoritative file is absent.
  idx_key <- if (pet_type == "PM") "spei" else "spei_era5"
  df <- tryCatch(load_flat_csv(idx_key, sc), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) {
    df <- tryCatch(load_spei_seasonal_csvs(sc, seas_dir), error = function(e) NULL)
  }
  
  if (is.null(df) || !nrow(df))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("No data\nSPEI-%d (%s)", sc, pet_type),
                               size = 4.5, colour = "grey50") +
             ggplot2::theme_void())
  
  df$date  <- as.Date(df$date)
  df$value <- as.numeric(df$value)
  
  df_fill <- df
  df_fill$fill_cat <- factor(dplyr::case_when(
    df_fill$value >  0    ~ "Wet",
    df_fill$value > -0.5  ~ "Mild deficit",
    df_fill$value > -1.0  ~ "Moderate drought",
    df_fill$value > -1.5  ~ "Severe drought",
    TRUE                  ~ "Extreme drought"),
    levels = c("Wet","Mild deficit","Moderate drought",
               "Severe drought","Extreme drought"))
  df_fill$date_end <- c(df_fill$date[-1L], next_month_start(df_fill$date))
  
  ylo <- min(-2.8, min(df$value, na.rm = TRUE) * 1.05)
  yhi <- max( 2.2, max(df$value, na.rm = TRUE) * 1.05)
  
  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
                      xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                      ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.55) +
    ggplot2::geom_rect(
      data = df_fill,
      ggplot2::aes(xmin = date, xmax = date_end,
                   ymin = -Inf, ymax = Inf, fill = fill_cat),
      alpha = 0.30, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(
      values = fill_pal, name = "Drought category", drop = FALSE,
      guide  = if (is_first)
        ggplot2::guide_legend(title.position = "top", nrow = 1,
                              override.aes = list(alpha = 0.7))
      else "none") +
    ggplot2::geom_hline(yintercept = 0,
                        colour = "grey40", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = THR_EVENT,
                        colour = "grey30", linewidth = 0.50, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = THR_MOD,
                        colour = "grey50", linewidth = 0.35, linetype = "dotted") +
    ggplot2::geom_line(
      data = df, ggplot2::aes(x = date, y = value),
      colour = col_line, linewidth = 0.55, inherit.aes = FALSE) +
    ggplot2::annotate("text",
                      x = min(df$date) + 365, y = yhi * 0.88,
                      label = panel_lab, size = 3.2, fontface = "bold",
                      hjust = 0, vjust = 1, colour = "grey10") +
    ggplot2::annotate("text",
                      x = as.Date("1951-01-01"), y = THR_EVENT + 0.08,
                      label = "x₀ = −0.5", size = 2.4,
                      hjust = 0, colour = "grey30") +
    ggplot2::annotate("text",
                      x = HIGHLIGHT_START + 365, y = yhi * 0.97,
                      label = "2022–2025", size = 2.4, hjust = 0.5,
                      colour = "grey30", fontface = "italic") +
    ggplot2::scale_x_date(
      date_breaks = "10 years", date_labels = "%Y",
      expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::scale_y_continuous(
      limits = c(ylo, yhi),
      breaks = c(-2, -1.5, -1, -0.5, 0, 1, 2),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title = sprintf("SPEI-%d  (%s,  area-weighted basin mean)", sc, pet_label),
      x     = if (show_xlab) "Year" else NULL,
      y     = sprintf("SPEI-%d", sc)) +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 8, face = "bold",
                                           hjust = 0,
                                           margin = ggplot2::margin(b = 2)),
      axis.text.x  = if (show_xlab) ggplot2::element_text(size = 7)
      else ggplot2::element_blank(),
      axis.ticks.x = if (show_xlab) ggplot2::element_line()
      else ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      axis.title.y = ggplot2::element_text(size = 7.5, angle = 90),
      axis.title.x = ggplot2::element_text(size = 7.5),
      legend.position = "bottom",
      legend.text  = ggplot2::element_text(size = 6.5),
      legend.key.size = ggplot2::unit(0.35, "cm"),
      legend.title = ggplot2::element_text(size = 7),
      panel.grid   = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(3, 5, 2, 4))
  
  p
}

# ── Overlay panel builder (PM + ERA5-Land PEV on the same axes) ─────────────────
make_spei_overlay_panel <- function(sc, panel_lab, show_xlab = FALSE,
                                    is_first = FALSE) {
  load_one <- function(pet_type) {
    seas_dir <- if (pet_type == "PM") SPEI_PM_DIR else SPEI_ERA5_DIR
    idx_key <- if (pet_type == "PM") "spei" else "spei_era5"
    df <- tryCatch(load_flat_csv(idx_key, sc), error = function(e) NULL)
    if (is.null(df) || !nrow(df))
      df <- tryCatch(load_spei_seasonal_csvs(sc, seas_dir), error = function(e) NULL)
    if (!is.null(df) && nrow(df)) {
      df$date    <- as.Date(df$date)
      df$value   <- as.numeric(df$value)
      df$pet_type <- pet_type
    }
    df
  }
  df_pm  <- load_one("PM")
  df_era5 <- load_one("ERA5")
  
  if (is.null(df_pm) && is.null(df_era5))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("No data\nSPEI-%d", sc),
                               size = 4.5, colour = "grey50") +
             ggplot2::theme_void())
  
  # Use PM for background fill if available, else ERA5
  df_fill_src <- if (!is.null(df_pm)) df_pm else df_era5
  df_fill <- df_fill_src
  df_fill$fill_cat <- factor(dplyr::case_when(
    df_fill$value >  0    ~ "Wet",
    df_fill$value > -0.5  ~ "Mild deficit",
    df_fill$value > -1.0  ~ "Moderate drought",
    df_fill$value > -1.5  ~ "Severe drought",
    TRUE                  ~ "Extreme drought"),
    levels = c("Wet","Mild deficit","Moderate drought",
               "Severe drought","Extreme drought"))
  df_fill$date_end <- c(df_fill$date[-1L], next_month_start(df_fill$date))
  
  all_vals <- c(if (!is.null(df_pm))  df_pm$value,
                if (!is.null(df_era5)) df_era5$value)
  ylo <- min(-2.8, min(all_vals, na.rm = TRUE) * 1.05)
  yhi <- max( 2.2, max(all_vals, na.rm = TRUE) * 1.05)
  
  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
                      xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                      ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.55) +
    ggplot2::geom_rect(
      data = df_fill,
      ggplot2::aes(xmin = date, xmax = date_end,
                   ymin = -Inf, ymax = Inf, fill = fill_cat),
      alpha = 0.25, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(
      values = fill_pal, name = "Drought category", drop = FALSE,
      guide  = if (is_first)
        ggplot2::guide_legend(title.position = "top", nrow = 1,
                              override.aes = list(alpha = 0.6))
      else "none") +
    ggplot2::geom_hline(yintercept = 0,     colour = "grey40", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = THR_EVENT, colour = "grey30",
                        linewidth = 0.50, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = THR_MOD, colour = "grey50",
                        linewidth = 0.35, linetype = "dotted")
  
  if (!is.null(df_pm))
    p <- p + ggplot2::geom_line(
      data = df_pm, ggplot2::aes(x = date, y = value,
                                 colour = "SPEI-PM (Penman-Monteith)"),
      linewidth = 0.65, inherit.aes = FALSE)
  
  if (!is.null(df_era5))
    p <- p + ggplot2::geom_line(
      data = df_era5, ggplot2::aes(x = date, y = value,
                                   colour = "SPEI-ERA5 (ERA5-Land PEV)"),
      linewidth = 0.55, linetype = "solid", alpha = 0.75, inherit.aes = FALSE)
  
  p <- p +
    ggplot2::scale_colour_manual(
      values = c("SPEI-PM (Penman-Monteith)" = "#1f77b4",
                 "SPEI-ERA5 (ERA5-Land PEV)"          = "#d62728"),
      name   = "PET method",
      guide  = if (is_first)
        ggplot2::guide_legend(title.position = "top", nrow = 1,
                              override.aes = list(linewidth = 1.2))
      else "none") +
    ggplot2::annotate("text",
                      x = min(df_fill_src$date) + 365, y = yhi * 0.88,
                      label = panel_lab, size = 3.2, fontface = "bold",
                      hjust = 0, vjust = 1, colour = "grey10") +
    ggplot2::annotate("text",
                      x = as.Date("1951-01-01"), y = THR_EVENT + 0.08,
                      label = "x₀ = −0.5", size = 2.4,
                      hjust = 0, colour = "grey30") +
    ggplot2::annotate("text",
                      x = HIGHLIGHT_START + 365, y = yhi * 0.97,
                      label = "2022–2025", size = 2.4, hjust = 0.5,
                      colour = "grey30", fontface = "italic") +
    ggplot2::scale_x_date(breaks = seq(as.Date("1950-01-01"), as.Date("2025-12-31"), by = "5 years"), date_labels = "%Y", 
                          expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::scale_y_continuous(limits = c(ylo, yhi),
                                breaks = c(-2, -1.5, -1, -0.5, 0, 1, 2),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title = sprintf("SPEI-%d  (area-weighted basin mean)", sc),
      x     = if (show_xlab) "Year" else NULL,
      y     = sprintf("SPEI-%d", sc)) +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 8, face = "bold", hjust = 0,
                                           margin = ggplot2::margin(b = 2)),
      axis.text.x  = if (show_xlab) ggplot2::element_text(size = 7)
      else ggplot2::element_blank(),
      axis.ticks.x = if (show_xlab) ggplot2::element_line()
      else ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      axis.title.y = ggplot2::element_text(size = 7.5, angle = 90),
      axis.title.x = ggplot2::element_text(size = 7.5),
      legend.position = "bottom",
      legend.text  = ggplot2::element_text(size = 6.5),
      legend.key.size = ggplot2::unit(0.35, "cm"),
      legend.title = ggplot2::element_text(size = 7),
      panel.grid   = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(3, 5, 2, 4))
  
  p
}

# ── Helper: assemble + save a 3-panel figure ───────────────────────────────────
save_spei_fig <- function(panels, title_str, subtitle_str, base_name) {
  fig <- panels[[1]] / panels[[2]] / panels[[3]] +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title    = title_str,
      subtitle = subtitle_str,
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40",
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 4)),
        legend.position = "bottom",
        legend.text     = ggplot2::element_text(size = 7.5),
        legend.key.size = ggplot2::unit(0.38, "cm")))
  
  for (ext in c("pdf", "png")) {
    out <- file.path(BASIN_PLOT_DIR, paste0(base_name, ".", ext))
    tryCatch(
      ggplot2::ggsave(out, fig, width = 7.0, height = 9.0, units = "in",
                      dpi = if (ext == "png") 300 else "print",
                      device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png"),
      error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message)))
    cat(sprintf("  ✓ Saved: %s\n", basename(out)))
  }
  invisible(fig)
}

common_subtitle <- paste0(
  "Fill: drought category  |  ",
  "Dashed: event threshold (x₀ = −0.5, 31st percentile)  |  ",
  "Dotted: moderate drought (−1.0)  |  ",
  "Grey band: 2022–2025  |  ",
  "Red segment: persistent sub-threshold anomaly 2023–2025 (SPEI-3)")

tryCatch({
  # ── Version A: PM only ─────────────────────────────────────────────────────
  cat("  Building Version A: SPEI-PM (Penman-Monteith)...\n")
  pm_panels <- list(
    make_spei_panel_v2(1, "(a)", "PM", show_xlab = FALSE, is_first = TRUE),
    make_spei_panel_v2(2, "(b)", "PM", show_xlab = FALSE, is_first = FALSE),
    make_spei_panel_v2(3, "(c)", "PM", show_xlab = TRUE,  is_first = FALSE)
  )
  save_spei_fig(pm_panels,
                "Basin-averaged SPEI-PM (Penman-Monteith PET) — Nechako River Basin (1950–2025)",
                common_subtitle,
                "Fig4_MS_SPEI_123_PM_timeseries")
  cat("  ✓ Version A complete\n")
}, error = function(e) cat(sprintf("  ⚠ Version A failed: %s\n", e$message)))

tryCatch({
  # ── Version B: ERA5-Land PEV only ──────────────────────────────────────────
  cat("  Building Version B: SPEI-ERA5 (ERA5-Land PEV)...\n")
  era5_panels <- list(
    make_spei_panel_v2(1, "(a)", "ERA5", show_xlab = FALSE, is_first = TRUE),
    make_spei_panel_v2(2, "(b)", "ERA5", show_xlab = FALSE, is_first = FALSE),
    make_spei_panel_v2(3, "(c)", "ERA5", show_xlab = TRUE,  is_first = FALSE)
  )
  save_spei_fig(era5_panels,
                "Basin-averaged SPEI-T (ERA5-Land PEV PET) — Nechako River Basin (1950–2025)",
                common_subtitle,
                "Fig4_MS_SPEI_123_ERA5_timeseries")
  cat("  ✓ Version B complete\n")
}, error = function(e) cat(sprintf("  ⚠ Version B failed: %s\n", e$message)))

tryCatch({
  # ── Version C: Overlay (PM + ERA5-Land PEV on same axes) ───────────────────
  cat("  Building Version C: overlay (PM + ERA5-Land PEV)...\n")
  ov_panels <- list(
    make_spei_overlay_panel(1, "(a)", show_xlab = FALSE, is_first = TRUE),
    make_spei_overlay_panel(2, "(b)", show_xlab = FALSE, is_first = FALSE),
    make_spei_overlay_panel(3, "(c)", show_xlab = TRUE,  is_first = FALSE)
  )
  save_spei_fig(ov_panels,
                "Basin-averaged SPEI — Nechako River Basin (1950–2025)",
                paste0(common_subtitle,
                       "  |  Blue = PM;  Red = ERA5-Land PEV"),
                "Fig4_MS_SPEI_123_overlay_timeseries")
  cat("  ✓ Version C complete\n")
}, error = function(e) cat(sprintf("  ⚠ Version C failed: %s\n", e$message)))

################################################################################
# PART 4c – SPATIAL REPRESENTATIVENESS FIGURE  (fig3)
# Output:
#   Fig4c_Pearson_r_SPEI_PM_SPI.pdf/.png
# Layout:
#   Rows 1-2: 2x3 Pearson r maps
# Row 1 (a-c): SPEI_PM scales 1, 3, 6  (Blues palette)
# Row 2 (d-f): SPI     scales 1, 3, 6  (Greens palette)
# Bottom row: 2 density panels side by side
# (g) Pixel Pearson r — Penman-Monteith PET (SPEI-1/2/3 PM)
# (h) Pixel Pearson r — SPI (scales 1, 3, 6)
# All panels use ggplot2 + patchwork (no terra::plot).
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4c: Spatial representativeness — SPEI-PM + SPI Pearson r maps\n")
cat("══════════════════════════════════════════════════════════════\n")

## Helper: detect basin boundary and return sf object in specified CRS
detect_basin_sf <- function(coords_df, target_crs = terra::crs(analysis_template)) {
  tryCatch({
    # Convert SpatVector to sf if needed
    if (inherits(basin_shp, "SpatVector")) {
      basin_sf <- sf::st_as_sf(basin_shp)
    } else {
      basin_sf <- basin_shp
    }
    
    is_projected <- max(abs(coords_df$lon), na.rm = TRUE) > 200
    if (is_projected) {
      basin_sf <- sf::st_transform(basin_sf, "EPSG:3005")
    } else {
      basin_sf <- sf::st_transform(basin_sf, "EPSG:4326")
    }
    basin_sf
  }, error = function(e) {
    cat(sprintf("  ⚠ detect_basin_sf failed: %s\n", e$message))
    NULL
  })
}

## Helper: build one Pearson r map panel ────────────────────────────────
## Transform coordinates to WGS84 (degrees) using terra::crds()
transform_coords_to_degrees <- function(coords_df) {
  # If coordinates are already geographic (degrees), return as-is
  if (max(abs(coords_df$lon), na.rm = TRUE) <= 200) {
    return(data.frame(lon = coords_df$lon, lat = coords_df$lat))
  }
  # Otherwise they are in BC Albers (EPSG:3005) — project to WGS84
  v      <- terra::vect(coords_df, geom = c("lon", "lat"), crs = "EPSG:3005")
  v_4326 <- terra::project(v, "EPSG:4326")
  crds   <- terra::crds(v_4326)
  data.frame(lon = crds[, 1], lat = crds[, 2])
}
#
calc_r_fill_limits <- function(r_list, floor_to = 0.0, ceil_to = 1.0) {
  vals <- unlist(r_list)
  vals <- vals[is.finite(vals)]
  
  if (!length(vals)) {
    return(list(
      limits = c(0.0, 1.0),
      breaks = seq(0.0, 1.0, by = 0.1)
    ))
  }
  
  r_min <- max(floor_to, floor(min(vals, na.rm = TRUE) * 10) / 10)
  r_max <- min(ceil_to,  ceiling(max(vals, na.rm = TRUE) * 10) / 10)
  
  if (!is.finite(r_min) || !is.finite(r_max) || r_max <= r_min) {
    r_min <- 0.0
    r_max <- 1.0
  }
  
  brks <- seq(r_min, r_max, by = 0.10)
  if (length(brks) < 2L) brks <- c(r_min, r_max)
  
  list(limits = c(r_min, r_max), breaks = brks)
}
#
make_map_panel <- function(r_v, coords, basin_sf, title_str, subtitle_str,
                           palette = "Blues",
                           snap_m = 1000,
                           overfill = 1.01,
                           show_legend = TRUE,
                           legend_title = "Pearson r",
                           fill_limits = NULL,
                           fill_breaks = NULL) {
  
  df_points <- data.frame(
    x = coords$lon,
    y = coords$lat,
    r = as.numeric(r_v)
  )
  df_points <- df_points[
    is.finite(df_points$r) &
      is.finite(df_points$x) &
      is.finite(df_points$y), ]
  
  if (!nrow(df_points)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  
  is_geo_coords <- max(abs(df_points$x), na.rm = TRUE) <= 200
  snap <- if (is_geo_coords) 0.01 else snap_m
  df_points$x_r <- if (snap == 0) df_points$x else round(df_points$x / snap) * snap
  df_points$y_r <- if (snap == 0) df_points$y else round(df_points$y / snap) * snap
  
  df_grid <- aggregate(r ~ x_r + y_r, data = df_points, FUN = mean)
  
  if (!nrow(df_grid)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  
  ux <- sort(unique(df_grid$x_r))
  uy <- sort(unique(df_grid$y_r))
  
  dx <- if (length(ux) >= 2L) stats::median(diff(ux), na.rm = TRUE) else snap
  dy <- if (length(uy) >= 2L) stats::median(diff(uy), na.rm = TRUE) else snap
  
  if (!is.finite(dx) || dx <= 0) dx <- snap
  if (!is.finite(dy) || dy <= 0) dy <- snap
  
  x_lim <- range(df_grid$x_r, na.rm = TRUE)
  y_lim <- range(df_grid$y_r, na.rm = TRUE)
  
  x_pad <- max(dx * 0.55, diff(x_lim) * 0.03)
  y_pad <- max(dy * 0.55, diff(y_lim) * 0.03)
  
  # Shared limits/breaks across all comparable panels
  if (is.null(fill_limits)) {
    r_min <- max(0.0, floor(min(df_grid$r, na.rm = TRUE) * 10) / 10)
    r_max <- min(1.0, ceiling(max(df_grid$r, na.rm = TRUE) * 10) / 10)
    
    if (!is.finite(r_min) || !is.finite(r_max) || r_max <= r_min) {
      r_min <- 0.0
      r_max <- 1.0
    }
    
    fill_limits <- c(r_min, r_max)
  }
  
  if (is.null(fill_breaks)) {
    fill_breaks <- seq(fill_limits[1], fill_limits[2], by = 0.10)
    if (length(fill_breaks) < 2L) fill_breaks <- fill_limits
  }
  
  if (palette == "Blues") {
    col_low  <- "#c6dbef"
    col_high <- "#08519c"
  } else if (palette == "Greens") {
    col_low  <- "#bae4bc"
    col_high <- "#005824"
  } else {
    brewer_cols <- RColorBrewer::brewer.pal(9, palette)
    col_low  <- brewer_cols[2]
    col_high <- brewer_cols[9]
  }
  
  ## ── Axis ticks in native plot coordinates ────────────────────────────────────
  ## Target geographic graticule values (no CRS transformation of the plot):
  ##   Longitude: 128°W to 123°W (6 ticks, step 1°), labelled as e.g. "128°W"
  ##   Latitude : 53°N  to 56°N  (4 ticks, step 1°), labelled as e.g. "53°N"
  ##
  ## Strategy: detect coordinate system from value range.
  ##   • Geographic (lon in −200..200): tick positions are the degree values.
  ##   • Projected (BC Albers EPSG:3005): reproject a reference grid of
  ##     lon/lat points to Albers and interpolate the native x/y for each
  ##     desired graticule line, using a linear mapping fitted to the data.
  
  lon_targets <- seq(-128, -123, by = 1)   # 6 ticks
  lat_targets <- seq(53, 56, by = 1)       # 4 ticks
  lon_labels  <- paste0(abs(lon_targets), "\u00b0W")
  lat_labels  <- paste0(lat_targets, "\u00b0N")
  
  is_geo <- max(abs(df_grid$x_r), na.rm = TRUE) <= 200
  
  if (is_geo) {
    ## Geographic coords — tick positions are the degree values directly
    x_breaks <- lon_targets
    y_breaks <- lat_targets
  } else {
    ## Projected coords (BC Albers) — map graticule degrees → native units
    ## Build a small reference grid covering the basin and project it
    ref_lon <- seq(-130, -120, by = 0.5)
    ref_lat <- seq(52,   58,   by = 0.5)
    ref_grid <- expand.grid(lon = ref_lon, lat = ref_lat)
    ref_vect <- terra::vect(ref_grid, geom = c("lon", "lat"), crs = "EPSG:4326")
    ref_alb  <- terra::project(ref_vect, "EPSG:3005")
    ref_crds <- terra::crds(ref_alb)
    ref_df   <- data.frame(
      lon    = ref_grid$lon,
      lat    = ref_grid$lat,
      x_alb  = ref_crds[, 1],
      y_alb  = ref_crds[, 2]
    )
    ## Linear model: native_x ~ lon,  native_y ~ lat
    ## (adequate for a ~5°×5° domain; avoids full reprojection of plot axes)
    lm_x     <- stats::lm(x_alb ~ lon, data = ref_df)
    lm_y     <- stats::lm(y_alb ~ lat, data = ref_df)
    x_breaks <- stats::predict(lm_x, newdata = data.frame(lon = lon_targets))
    y_breaks <- stats::predict(lm_y, newdata = data.frame(lat = lat_targets))
    ## Keep only ticks that fall within the plot extent (with a small margin)
    x_margin <- diff(x_lim) * 0.05
    y_margin <- diff(y_lim) * 0.05
    keep_x   <- x_breaks >= (x_lim[1] - x_margin) & x_breaks <= (x_lim[2] + x_margin)
    keep_y   <- y_breaks >= (y_lim[1] - y_margin) & y_breaks <= (y_lim[2] + y_margin)
    x_breaks  <- x_breaks[keep_x];  lon_labels <- lon_labels[keep_x]
    y_breaks  <- y_breaks[keep_y];  lat_labels <- lat_labels[keep_y]
  }
  
  dx <- stats::median(
    diff(sort(unique(df_grid$x_r))),
    na.rm = TRUE
  )
  
  dy <- stats::median(
    diff(sort(unique(df_grid$y_r))),
    na.rm = TRUE
  )
  
  p_map <- ggplot2::ggplot(
    df_grid,
    ggplot2::aes(x = x_r, y = y_r, fill = r)
  ) +
    ggplot2::geom_tile(
      width  = dx,
      height = dy
    ) +
    ggplot2::scale_fill_gradient(
      low = col_low,
      high = col_high,
      limits = fill_limits,
      breaks = fill_breaks,
      name = legend_title,
      oob = scales::squish,
      na.value = "transparent",
      guide = if (show_legend) {
        ggplot2::guide_colourbar(
          title.position = "top",
          barheight = grid::unit(2.8, "cm"),
          barwidth  = grid::unit(0.35, "cm"),
          frame.colour = "grey50",
          ticks.colour = "grey40"
        )
      } else {
        "none"
      }
    ) +
    ggplot2::coord_cartesian(
      xlim = c(x_lim[1] - x_pad, x_lim[2] + x_pad),
      ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad),
      expand = FALSE
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = x_breaks,
      labels = lon_labels
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      breaks = y_breaks,
      labels = lat_labels
    ) +
    ggplot2::labs(
      title    = title_str,
      subtitle = subtitle_str,
      x        = "Longitude (\u00b0W)",
      y        = "Latitude (\u00b0N)"
    ) +
    theme_map4c +
    ggplot2::theme(
      aspect.ratio  = 1,
      axis.title.x  = ggplot2::element_text(size = 6.5, margin = ggplot2::margin(t = 2)),
      axis.title.y  = ggplot2::element_text(size = 6.5, margin = ggplot2::margin(r = 2)),
      axis.text.x   = ggplot2::element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
      axis.text.y   = ggplot2::element_text(size = 6),
      axis.ticks     = ggplot2::element_line(linewidth = 0.3, colour = "grey40"),
      axis.ticks.length = ggplot2::unit(2, "pt")
    )
  if (!is.null(basin_sf) && inherits(basin_sf, "sf")) {
    p_map <- p_map + ggplot2::geom_sf(data = basin_sf, inherit.aes = FALSE,
                                      fill = NA, colour = "black",
                                      linewidth = 0.35, show.legend = FALSE)
  }
  p_map
}
## ── Part 4c / 4c-2 shared helpers (global scope) ─────────────────────────────
## Defined here — OUTSIDE both tryCatch blocks — so that Part 4c-2 can always
## find them even if Part 4c encountered an error before assigning them.

## Helper: load 12 per-month CSVs for one scale x directory x prefix ──────────
load_pixel_mat <- function(scale, seas_dir, prefix = "spei") {
  month_list   <- vector("list", 12L)
  coords_out   <- NULL
  year_out     <- NULL
  for (m in seq_len(12L)) {
    csv_f   <- file.path(seas_dir,
                         sprintf("%s_%02d_month%02d_%s.csv",
                                 prefix, scale, m, MONTH_NAMES_3[m]))
    if (!file.exists(csv_f)) next
    df_m   <- tryCatch(data.table::fread(csv_f, data.table = FALSE),
                       error = function(e) NULL)
    if (is.null(df_m)) next
    yr_c   <- setdiff(names(df_m), c("lon", "lat"))
    yr_c   <- yr_c[grepl("^[0-9]{4}$", yr_c)]
    if (!length(yr_c)) next
    if (is.null(coords_out)) { coords_out <- df_m[, c("lon", "lat")]; year_out <- yr_c }
    month_list[[m]] <- as.matrix(df_m[, yr_c, drop = FALSE])
  }
  valid <- which(!sapply(month_list, is.null))
  if (!length(valid) || is.null(coords_out)) return(NULL)
  list(mat          = do.call(cbind, month_list[valid]),
       coords       = coords_out,
       year_int     = as.integer(year_out),
       yr_chr       = year_out,
       valid_months = valid,
       month_list   = month_list)
}

## Helper: per-pixel Pearson r with basin mean ─────────────────────────────────
## Basin mean = fractional-coverage area-weighted mean (0.20 coverage rule).
compute_r <- function(obj) {
  mat    <- obj$mat
  coords <- obj$coords
  bavg   <- s12_basin_mean_xy(mat, coords$lon, coords$lat)
  vapply(seq_len(nrow(mat)), function(i) {
    xi <- mat[i, ]; ok <- !is.na(xi) & !is.na(bavg)
    if (sum(ok) < 10L) NA_real_ else cor(xi[ok], bavg[ok])
  }, numeric(1L))
}

## Shared map theme ─────────────────────────────────────────────────────────────
theme_map4c <- ggplot2::theme_bw(base_size = 10.5) +
  ggplot2::theme(
    axis.title        = ggplot2::element_blank(),
    axis.text         = ggplot2::element_text(size = 7),
    axis.text.x       = ggplot2::element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
    axis.text.y       = ggplot2::element_text(size = 7),
    plot.title        = ggplot2::element_text(size = 7.5, face = "bold", hjust = 0.5,
                                              lineheight = 1.1,
                                              margin = ggplot2::margin(t = 2, b = 2)),
    plot.subtitle     = ggplot2::element_text(size = 7, colour = "grey30",
                                              hjust = 0.5,
                                              lineheight = 0.9,
                                              margin = ggplot2::margin(t = 0, b = 2)),
    legend.key.height = ggplot2::unit(0.8, "cm"),
    legend.key.width  = ggplot2::unit(0.25, "cm"),
    legend.text       = ggplot2::element_text(size = 8),
    legend.title      = ggplot2::element_text(size = 8.5),
    legend.position   = "right",
    panel.grid        = ggplot2::element_blank(),
    plot.margin       = ggplot2::margin(5, 3, 1, 3))

## Helper: build one r-density panel ──────────────────────────────────────────
make_density_panel <- function(r_list, combo_names, col_pal,
                               panel_lab, title_str) {
  df  <- do.call(rbind, lapply(seq_along(r_list), function(k) {
    rv  <- r_list[[k]]
    if (is.null(rv)) return(NULL)
    data.frame(r = rv, combo = combo_names[k], stringsAsFactors = FALSE)
  }))
  if (is.null(df) || !nrow(df))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = "No data", size = 4, colour = "grey50") +
             ggplot2::theme_void())
  df      <- df[!is.na(df$r), ]
  df$combo   <- factor(df$combo, levels = combo_names)
  med_df   <- do.call(rbind, lapply(combo_names, function(g) {
    x   <- df$r[df$combo == g]
    data.frame(combo = g, med_r = median(x, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  med_df$combo   <- factor(med_df$combo, levels = combo_names)
  
  ## ── Stagger median labels vertically so they don't collide ────────────────
  ## Rank by median r (ascending); assign vjust levels that stack labels
  ## from the top down: smallest median gets highest vjust (lowest on page),
  ## largest median gets smallest vjust (closest to top edge).
  ## Levels: 1.1  (≈ 1 line below top)
  ##         2.8  (≈ 2 lines below top)
  ##         4.5  (≈ 3 lines below top)
  rk <- rank(med_df$med_r, ties.method = "first")
  vjust_levels <- c(1.1, 2.8, 4.5)
  med_df$vjust_pos <- vjust_levels[rk]
  
  ## ── Nudge xintercepts slightly so overlapping lines stay visually distinct ──
  ## When all 3 medians cluster within ~0.01 r units, solid/dashed/dotted lines
  ## pile on top of each other.  A ±0.003 cosmetic offset (<0.5% of x range)
  ## separates the lines without distorting the displayed label (which still
  ## shows the true med_r value).
  n_sc     <- length(combo_names)
  offsets  <- seq(-0.003 * (n_sc - 1) / 2,
                  0.003 * (n_sc - 1) / 2,
                  length.out = n_sc)
  med_df$xint_nudged <- med_df$med_r +
    offsets[as.integer(med_df$combo)]
  
  ## Line types per scale: solid / dashed / dotted — visually separates
  ## lines that may be very close in x position
  lt_vals <- setNames(
    c("solid", "dashed", "dotted")[seq_along(combo_names)],
    combo_names)
  
  ## ── Per-scale density linewidths: thick→thin so overlapping curves stay ──────
  ## distinguishable even when their x distributions are nearly identical.
  ## Scale 1 = thickest (most prominent); scale 3 = thinnest.
  lw_vals <- setNames(
    c(1.6, 1.0, 0.45)[seq_along(combo_names)],
    combo_names)
  
  ## ── Label x-offset: place text just to the RIGHT of its vertical line ────────
  ## hjust = 0  → left edge of text anchored at x
  ## x = xint_nudged + 0.008  → small gap between line and start of label
  ## This prevents the line from bisecting the numeral.
  LABEL_GAP <- 0.008
  
  ggplot2::ggplot(df, ggplot2::aes(x = r, colour = combo, fill = combo)) +
    ## Map linewidth to combo so each density curve has a distinct stroke weight
    ggplot2::geom_density(ggplot2::aes(linewidth = combo), alpha = 0.20) +
    ggplot2::geom_vline(data = med_df,
                        ggplot2::aes(xintercept = xint_nudged, colour = combo,
                                     linetype = combo),
                        linewidth = 0.75,
                        show.legend = FALSE) +
    ## Text: staggered vertically AND shifted right of the line
    ggplot2::geom_text(data = med_df,
                       ggplot2::aes(x     = xint_nudged + LABEL_GAP,
                                    y     = Inf,
                                    label = sprintf("%.2f", med_r),
                                    colour = combo,
                                    vjust  = vjust_pos),
                       hjust = 0, size = 3.2,
                       show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = col_pal, name = NULL) +
    ggplot2::scale_fill_manual(values   = col_pal, name = NULL) +
    ggplot2::scale_linewidth_manual(values = lw_vals, guide = "none") +
    ggplot2::scale_linetype_manual(values = lt_vals) +
    ggplot2::scale_x_continuous(breaks = seq(0.4, 1.0, by = 0.1)) +
    ggplot2::coord_cartesian(xlim = c(0.4, 1.01)) +
    ## RESERVE ~60% OF TOP PANEL FOR LABELS to guarantee no curve overlap
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.70))) +
    ## NOTE: subtitle deliberately omitted here — a single shared note is added
    ##       once via the patchwork caption in the fig_dens assembly below,
    ##       preventing the identical text from appearing in both panels and
    ##       overflowing into the adjacent panel's space.
    ggplot2::labs(
      title = sprintf("%s  %s", panel_lab, title_str),
      x     = "Pearson r (pixel vs basin mean)",
      y     = "Density") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 8.5, face = "bold", hjust = 0,
                                              margin = ggplot2::margin(t = 3, b = 3)),
      legend.text     = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(0.4, "cm"),
      axis.text       = ggplot2::element_text(size = 9),
      plot.margin     = ggplot2::margin(5, 6, 2, 4))
}

tryCatch({
  ## 1. Compute SPEI PM Pearson r for scales 1, 3, 6 ──────────────────────
  ### 1. Compute SPEI PM Pearson r for scales 1, 3, 6 ──────────────────────
  
  spei_scales <- c(1L, 2L, 3L)
  spei_r_list <- vector("list", 3L)
  spei_obj_list <- vector("list", 3L)
  spei_med_list <- rep(NA_real_, 3L)
  spei_pct80_list <- rep(NA_real_, 3L)
  spei_plots <- vector("list", 3L)
  spei_labels <- c("(a)", "(b)", "(c)")
  
  for (k in seq_along(spei_scales)) {
    sc <- spei_scales[k]
    cat(sprintf(" Loading SPEI-PM-%d...\n", sc))
    
    obj <- tryCatch(
      load_pixel_mat(sc, SPEI_PM_DIR, prefix = "spei"),
      error = function(e) NULL
    )
    
    if (is.null(obj)) {
      cat(" No data\n")
      next
    }
    
    rv <- compute_r(obj)
    
    spei_obj_list[[k]] <- obj
    spei_r_list[[k]] <- rv
    spei_med_list[k] <- median(rv, na.rm = TRUE)
    spei_pct80_list[k] <- 100 * mean(rv >= 0.80, na.rm = TRUE)
  }
  
  # Shared colour scale for all 3 SPEI-PM maps
  spei_fill_cfg <- calc_r_fill_limits(spei_r_list)
  
  # Build SPEI-PM map panels AFTER shared limits are known
  for (k in seq_along(spei_scales)) {
    obj <- spei_obj_list[[k]]
    rv  <- spei_r_list[[k]]
    
    if (is.null(obj) || is.null(rv)) next
    
    bsf <- detect_basin_sf(obj$coords)
    
    spei_plots[[k]] <- make_map_panel(
      rv, obj$coords, bsf,
      title_str = sprintf("%s SPEI-%d (PM)", spei_labels[k], spei_scales[k]),
      subtitle_str = sprintf("median r = %.3f\n%.0f%% \u2265 0.80",
                             spei_med_list[k], spei_pct80_list[k]),
      palette = "Blues",
      show_legend = TRUE,
      legend_title = "Pearson r (SPEI-PM)",
      fill_limits = spei_fill_cfg$limits,
      fill_breaks = spei_fill_cfg$breaks
    )
  }
  
  
  ## 2. Compute SPI Pearson r for scales 1, 3, 6 ──────────────────────────
  
  spi_scales <- c(1L, 2L, 3L)
  spi_r_list <- vector("list", 3L)
  spi_obj_list <- vector("list", 3L)
  spi_med_list <- rep(NA_real_, 3L)
  spi_pct80_list <- rep(NA_real_, 3L)
  spi_plots <- vector("list", 3L)
  spi_labels <- c("(d)", "(e)", "(f)")
  
  for (k in seq_along(spi_scales)) {
    sc <- spi_scales[k]
    cat(sprintf(" Loading SPI-%d...\n", sc))
    
    obj <- tryCatch(
      load_pixel_mat(sc, SPI_SEAS_DIR, prefix = "spi"),
      error = function(e) NULL
    )
    
    if (is.null(obj)) {
      cat(" No data\n")
      next
    }
    
    rv <- compute_r(obj)
    
    spi_obj_list[[k]] <- obj
    spi_r_list[[k]] <- rv
    spi_med_list[k] <- median(rv, na.rm = TRUE)
    spi_pct80_list[k] <- 100 * mean(rv >= 0.80, na.rm = TRUE)
  }
  
  # Shared colour scale for all 3 SPI maps
  spi_fill_cfg <- calc_r_fill_limits(spi_r_list)
  
  # Build SPI map panels AFTER shared limits are known
  for (k in seq_along(spi_scales)) {
    obj <- spi_obj_list[[k]]
    rv  <- spi_r_list[[k]]
    
    if (is.null(obj) || is.null(rv)) next
    
    bsf <- detect_basin_sf(obj$coords)
    
    spi_plots[[k]] <- make_map_panel(
      rv, obj$coords, bsf,
      title_str = sprintf("%s SPI-%d", spi_labels[k], spi_scales[k]),
      subtitle_str = sprintf("median r = %.3f\n%.0f%% ≥ 0.80",
                             spi_med_list[k], spi_pct80_list[k]),
      palette = "Greens",
      show_legend = TRUE,
      legend_title = "Pearson r (SPI)",
      fill_limits = spi_fill_cfg$limits,
      fill_breaks = spi_fill_cfg$breaks
    )
  }
  
  
  any_spei <- any(!sapply(spei_plots, is.null))
  any_spi  <- any(!sapply(spi_plots,  is.null))
  if (!any_spei && !any_spi) { cat("  No data for any combination\n"); stop("no_data") }
  
  ## 3. Density panel (g): SPEI-PM r values for scales 1, 3, 6 ────────────
  col_spei_pm   <- c("SPEI-PM-1" = "#1f77b4",
                     "SPEI-PM-3" = "#6baed6",
                     "SPEI-PM-6" = "#2ca02c")
  pc_pm   <- make_density_panel(
    r_list      = spei_r_list,
    combo_names = names(col_spei_pm),
    col_pal     = col_spei_pm,
    panel_lab   = "(g)",
    title_str   = "Pixel Pearson r — Penman-Monteith PET")
  
  ## 4. Density panel (h): SPI r values for scales 1, 3, 6 ────────────────
  col_spi   <- c("SPI-1" = "#1b7837",
                 "SPI-3" = "#5aae61",
                 "SPI-6" = "#a6d96a")
  pl_spi   <- make_density_panel(
    r_list      = spi_r_list,
    combo_names = names(col_spi),
    col_pal     = col_spi,
    panel_lab   = "(h)",
    title_str   = "Pixel Pearson r — SPI (scales 1, 3, 6)")
  
  ## 5. Assemble figure ────────────────────────────────────────────────────
  top_plots <- c(spei_plots, spi_plots)
  top_plots <- lapply(top_plots, function(p) {
    if (is.null(p))
      ggplot2::ggplot() + ggplot2::theme_void()
    else p
  })
  
  ## Ensure NO legends on map panels (override any defaults)
  fig_maps <- patchwork::wrap_plots(top_plots, nrow = 2, ncol = 3, guides = "collect") &
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical"
    )
  
  ## Density panels keep their legends.
  ## The shared subtitle note ("Vertical lines = median r …") is placed ONCE as
  ## a patchwork caption below both panels — this prevents the identical text
  ## from appearing in each individual panel and overflowing into its neighbour.
  fig_dens <- (pc_pm | pl_spi) +
    patchwork::plot_layout(guides = "keep") +
    patchwork::plot_annotation(
      caption = "Vertical lines = median r per scale  (solid / dashed / dotted)",
      theme   = ggplot2::theme(
        plot.caption = ggplot2::element_text(
          size   = 7.5, colour = "grey40",
          hjust  = 0.5,
          margin = ggplot2::margin(t = 2, b = 2))))
  
  ## Combine: maps (no legends) + density plots (with their legends)
  fig4c  <- fig_maps / fig_dens +
    patchwork::plot_layout(heights = c(2, 0.7)) +
    patchwork::plot_annotation(
      title    = "Pixel-to-basin-mean Pearson r — SPEI 1,2,3 + SPI 1,2,3 — Nechako River Basin",
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5,
                                              margin = ggplot2::margin(b = 4)),
        plot.subtitle = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text     = ggplot2::element_text(size = 10),
        legend.key.size = ggplot2::unit(0.45, "cm")))
  
  ## Save with Climate Dynamics specs
  for (ext in c("pdf", "png")) {
    out_f   <- file.path(BASIN_PLOT_DIR,
                         paste0("Fig4c_Pearson_r_SPEI_PM_SPI.", ext))
    tryCatch(
      ggplot2::ggsave(out_f, fig4c,
                      width  = 6.85,
                      height = 9.21,
                      units  = "in",
                      dpi    = if (ext == "png") 600 else "print",
                      device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png",
                      limitsize = FALSE),
      error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message)))
    cat(sprintf("  ✓ Saved: %s\n", basename(out_f)))
  }
  
  ## Console summary
  cat("\n  Representativeness summary:\n")
  cat(sprintf("  %-14s  median_r  pct_r≥80\n", "Combination"))
  cat("    ", paste(rep("-", 38), collapse = ""), "\n", sep = "")
  for (k in seq_along(spei_scales)) {
    rv   <- spei_r_list[[k]]
    if (is.null(rv)) next
    cat(sprintf("  SPEI-PM-%d       %.3f    %.0f%%\n",
                spei_scales[k], median(rv, na.rm=TRUE),
                100 * mean(rv >= 0.80, na.rm=TRUE)))
  }
  for (k in seq_along(spi_scales)) {
    rv   <- spi_r_list[[k]]
    if (is.null(rv)) next
    cat(sprintf("  SPI-%d           %.3f    %.0f%%\n",
                spi_scales[k], median(rv, na.rm=TRUE),
                100 * mean(rv >= 0.80, na.rm=TRUE)))
  }
}, error = function(e) {
  if (!identical(conditionMessage(e), "no_data"))
    cat(sprintf("  ⚠ Part 4c failed: %s\n", e$message))
})

################################################################################
# PART 4c-2 – SPEI-PM 1/3/6 PM-ONLY SPATIAL REPRESENTATIVENESS FIGURE
# Output:
#   Fig4c_Pearson_r_SPEI136_PM.pdf/.png
# Layout (ERA5-Land PEV maps REMOVED — SPEI-PM panels only):
#   Row 1 (a-c): SPEI-PM scales 1, 3, 6  (Blues palette)
#   Bottom row : single density panel — Pixel Pearson r for SPEI-PM 1, 2, 3
# NOTE: This is a modified version of the PM+ERA5 figure; the 3 SPEI-ERA5
#       maps have been intentionally omitted at the author's request.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4c-2: SPEI123 PM-only spatial representativeness figure\n")
cat("          (ERA5-Land PEV maps excluded)\n")
cat("══════════════════════════════════════════════════════════════\n")

tryCatch({
  
  ## 1. Compute SPEI-PM Pearson r for scales 1, 3, 6 ─────────────────────────
  spei_scales_pm2  <- c(1L, 3L, 6L)
  spei_r_list_pm2  <- vector("list", 3L)
  spei_plots_pm2   <- vector("list", 3L)
  spei_labels_pm2  <- c("(a)", "(b)", "(c)")
  
  for (k in seq_along(spei_scales_pm2)) {
    sc  <- spei_scales_pm2[k]
    cat(sprintf("  Loading SPEI-PM-%d...\n", sc))
    obj <- tryCatch(load_pixel_mat(sc, SPEI_PM_DIR, prefix = "spei"),
                    error = function(e) NULL)
    if (is.null(obj)) { cat(sprintf("    No data for SPEI-PM-%d\n", sc)); next }
    rv        <- compute_r(obj)
    bsf       <- detect_basin_sf(obj$coords)
    med_r     <- median(rv, na.rm = TRUE)
    pct80     <- 100 * mean(rv >= 0.80, na.rm = TRUE)
    spei_r_list_pm2[[k]] <- rv
    spei_plots_pm2[[k]]  <- make_map_panel(
      rv, obj$coords, bsf,
      title_str    = sprintf("%s  SPEI-%d  (Penman-Monteith)", spei_labels_pm2[k], sc),
      subtitle_str = sprintf("median r = %.3f\n%.0f%% \u2265 0.80", med_r, pct80),
      palette      = "Blues")
  }
  
  pm2_fill_cfg <- calc_r_fill_limits(spei_r_list_pm2)
  for (k in seq_along(spei_scales_pm2)) {
    if (is.null(spei_plots_pm2[[k]]) || is.null(spei_r_list_pm2[[k]])) next
    sc <- spei_scales_pm2[k]
    obj <- tryCatch(load_pixel_mat(sc, SPEI_PM_DIR, prefix = "spei"), error = function(e) NULL)
    if (is.null(obj)) next
    spei_plots_pm2[[k]] <- make_map_panel(
      spei_r_list_pm2[[k]], obj$coords, detect_basin_sf(obj$coords),
      title_str = sprintf("%s SPEI-PM-%d", spei_labels_pm2[k], sc),
      subtitle_str = sprintf("median r = %.3f", median(spei_r_list_pm2[[k]], na.rm = TRUE)),
      palette = "Blues", show_legend = (k == 1L),
      fill_limits = pm2_fill_cfg$limits, fill_breaks = pm2_fill_cfg$breaks
    )
  }
  any_pm2 <- any(!sapply(spei_plots_pm2, is.null))
  if (!any_pm2) { cat("  No SPEI-PM data available\n"); stop("no_data") }
  
  ## 2. Single density panel (d): SPEI-PM r values for scales 1, 3, 6 ────────
  col_spei_pm2 <- c("SPEI-1 PM" = "#1f77b4",
                    "SPEI-2 PM" = "#6baed6",
                    "SPEI-3 PM" = "#2ca02c")
  pc_pm2 <- make_density_panel(
    r_list      = spei_r_list_pm2,
    combo_names = names(col_spei_pm2),
    col_pal     = col_spei_pm2,
    panel_lab   = "(d)",
    title_str   = "Pixel Pearson r \u2014 Penman-Monteith PET (scales 1, 3, 6)")
  
  ## 3. Assemble: 3 SPEI-PM maps (top) + density panel (bottom) ───────────────
  top_plots_pm2 <- lapply(spei_plots_pm2, function(p) {
    if (is.null(p)) ggplot2::ggplot() + ggplot2::theme_void() else p
  })
  
  fig_maps_pm2 <- patchwork::wrap_plots(top_plots_pm2, nrow = 1, ncol = 3)
  
  fig4c_pm2 <- fig_maps_pm2 / pc_pm2 +
    patchwork::plot_layout(heights = c(1.6, 0.7), guides = "collect") +
    patchwork::plot_annotation(
      title = paste0("Pixel-to-basin-mean Pearson r \u2014 SPEI-PM 1, 2, 3 \u2014 ",
                     "Nechako River Basin"),
      theme = ggplot2::theme(
        plot.title      = ggplot2::element_text(size = 12, face = "bold",
                                                hjust = 0.5,
                                                margin = ggplot2::margin(b = 4)),
        plot.subtitle   = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text     = ggplot2::element_text(size = 10),
        legend.key.size = ggplot2::unit(0.45, "cm")))
  
  ## 4. Save ───────────────────────────────────────────────────────────────────
  for (ext in c("pdf", "png")) {
    out_pm2 <- file.path(BASIN_PLOT_DIR,
                         paste0("Fig4c_Pearson_r_SPEI136_PM.", ext))
    tryCatch(
      ggplot2::ggsave(out_pm2, fig4c_pm2,
                      width     = 6.85,
                      height    = 6.50,   # shorter: one map row instead of two
                      units     = "in",
                      dpi       = if (ext == "png") 600 else "print",
                      device    = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png",
                      limitsize = FALSE),
      error = function(e) cat(sprintf("  \u26a0 %s: %s\n", ext, e$message)))
    cat(sprintf("  \u2713 Saved: %s\n", basename(out_pm2)))
  }
  
  ## 5. Console summary ────────────────────────────────────────────────────────
  cat("\n  SPEI-PM representativeness summary (PM-only figure):\n")
  cat(sprintf("  %-14s  median_r  pct_r\u226580\n", "Combination"))
  cat("    ", paste(rep("-", 38), collapse = ""), "\n", sep = "")
  for (k in seq_along(spei_scales_pm2)) {
    rv <- spei_r_list_pm2[[k]]
    if (is.null(rv)) next
    cat(sprintf("  SPEI-PM-%d       %.3f    %.0f%%\n",
                spei_scales_pm2[k],
                median(rv, na.rm = TRUE),
                100 * mean(rv >= 0.80, na.rm = TRUE)))
  }
  
}, error = function(e) {
  if (!identical(conditionMessage(e), "no_data"))
    cat(sprintf("  \u26a0 Part 4c-2 failed: %s\n", e$message))
})

cat("\n══ PART 4c-2 complete ══\n")

################################################################################
# PART 4c-3 — SNOW-INDEX SPATIAL REPRESENTATIVENESS
#
# Counterpart to the SPI/SPEI pixel-to-basin Pearson-r maps. SWEI-3 and SSPI-1
# are read from their native seasonal NetCDF files (not the SPI/SPEI CSV loader).
# Each pixel's monthly series is correlated with the grid-wide, latitude-weighted
# mean, matching the existing Part 4c diagnostic definition.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4c-3: Spatial representativeness — SWEI-3 + SSPI-1\n")
cat("══════════════════════════════════════════════════════════════\n")

load_snow_pixel_mat <- function(scale, seas_dir, prefix) {
  month_mats <- vector("list", 12L)
  template <- NULL
  coords_out <- NULL
  n_years <- NULL
  for (m in seq_len(12L)) {
    pat <- sprintf("^%s_%02d_month%02d_[A-Za-z]+\\.nc$", prefix, as.integer(scale), m)
    files <- list.files(seas_dir, pattern = pat, full.names = TRUE, ignore.case = TRUE,
                        recursive = TRUE)
    if (length(files) != 1L) {
      stop(sprintf("Expected one seasonal NetCDF for %s-%02d month %02d; found %d in %s",
                   toupper(prefix), scale, m, length(files), seas_dir))
    }
    r <- terra::rast(files[1L])
    if (is.null(template)) {
      template <- r[[1L]]
      xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
      coords_out <- data.frame(lon = xy[, 1L], lat = xy[, 2L])
    } else if (!isTRUE(terra::compareGeom(r[[1L]], template, stopOnError = FALSE,
                                          crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE))) {
      stop("Snow seasonal files do not share identical native raster geometry: ", basename(files[1L]))
    }
    vals <- as.matrix(terra::values(r, mat = TRUE))
    vals[!is.finite(vals)] <- NA_real_
    if (is.null(n_years)) n_years <- ncol(vals)
    if (ncol(vals) != n_years) {
      stop(sprintf("Inconsistent layer count across %s seasonal files: month %02d has %d, expected %d",
                   toupper(prefix), m, ncol(vals), n_years))
    }
    month_mats[[m]] <- vals
    rm(r, vals)
  }
  list(mat = do.call(cbind, month_mats), coords = coords_out,
       valid_months = seq_len(12L), index = prefix, scale = as.integer(scale))
}

tryCatch({
  snow_specs <- list(
    list(index = "swei", scale = SWEI_SCALE, dir = SWEI_SEAS_DIR,
         label = "SWEI-3", palette = "Blues", colour = "#1E90FF"),
    list(index = "sspi", scale = SSPI_SCALE, dir = SSPI_SEAS_DIR,
         label = "SSPI-1", palette = "Purples", colour = "#9467bd")
  )
  snow_r_list <- vector("list", length(snow_specs))
  names(snow_r_list) <- vapply(snow_specs, `[[`, character(1L), "label")
  snow_obj_list <- vector("list", length(snow_specs))
  snow_map_list <- vector("list", length(snow_specs))
  snow_medians <- rep(NA_real_, length(snow_specs))
  snow_pct80 <- rep(NA_real_, length(snow_specs))
  
  for (k in seq_along(snow_specs)) {
    spec <- snow_specs[[k]]
    cat(sprintf("  Loading %s from native seasonal NetCDFs...\n", spec$label))
    obj <- tryCatch(load_snow_pixel_mat(spec$scale, spec$dir, spec$index),
                    error = function(e) {
                      cat(sprintf("    ⚠ %s unavailable: %s\n", spec$label, e$message))
                      NULL
                    })
    if (is.null(obj)) next
    rv <- tryCatch(compute_r(obj), error = function(e) {
      cat(sprintf("    ⚠ Correlation failed for %s: %s\n", spec$label, e$message))
      NULL
    })
    if (is.null(rv)) next
    snow_obj_list[[k]] <- obj
    snow_r_list[[k]] <- rv
    snow_medians[k] <- median(rv, na.rm = TRUE)
    snow_pct80[k] <- 100 * mean(rv >= 0.80, na.rm = TRUE)
  }
  
  if (!any(vapply(snow_r_list, function(x) !is.null(x) && any(is.finite(x)), logical(1L)))) {
    cat("  ⚠ Neither snow index could produce valid pixel-to-basin correlations; figure skipped.\n")
  } else {
    snow_fill_cfg <- calc_r_fill_limits(snow_r_list)
    for (k in seq_along(snow_specs)) {
      spec <- snow_specs[[k]]
      obj <- snow_obj_list[[k]]
      rv <- snow_r_list[[k]]
      if (is.null(obj) || is.null(rv)) next
      snow_map_list[[k]] <- make_map_panel(
        rv, obj$coords, detect_basin_sf(obj$coords),
        title_str = sprintf("(%s) %s pixel-to-basin correlation", letters[k], spec$label),
        subtitle_str = sprintf("median r = %.3f; %.0f%% ≥ 0.80", snow_medians[k], snow_pct80[k]),
        palette = spec$palette, show_legend = (k == 1L),
        legend_title = "Pearson r", fill_limits = snow_fill_cfg$limits,
        fill_breaks = snow_fill_cfg$breaks
      )
    }
    snow_map_panels <- lapply(snow_map_list, function(p) {
      if (is.null(p)) ggplot2::ggplot() + ggplot2::theme_void() else p
    })
    snow_density <- make_density_panel(
      r_list = snow_r_list,
      combo_names = names(snow_r_list),
      col_pal = setNames(vapply(snow_specs, `[[`, character(1L), "colour"), names(snow_r_list)),
      panel_lab = "(c)",
      title_str = "Pixel Pearson r — snow-drought indices"
    )
    snow_spatial_fig <- patchwork::wrap_plots(snow_map_panels, nrow = 1L, ncol = 2L,
                                              guides = "collect") / snow_density +
      patchwork::plot_layout(heights = c(1.5, 0.7), guides = "collect") +
      patchwork::plot_annotation(
        title = "Pixel-to-basin-mean Pearson correlation — snow-drought indices",
        subtitle = "Native-grid SWEI-3 and SSPI-1; regional means are not recomputed here",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, colour = "grey35"),
          legend.position = "bottom"
        )
      )
    for (ext in c("pdf", "png")) {
      out <- file.path(BASIN_PLOT_DIR, paste0("Fig4c_Pearson_r_SWEI_SSPI.", ext))
      tryCatch({
        ggplot2::ggsave(out, snow_spatial_fig, width = 10.5, height = 7.8, units = "in",
                        dpi = if (ext == "png") 600 else "print",
                        device = if (ext == "pdf") {
                          if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
                        } else "png", limitsize = FALSE, bg = "white")
        cat(sprintf("  ✓ Saved: %s\n", basename(out)))
      }, error = function(e) cat(sprintf("  ⚠ %s: %s\n", basename(out), e$message)))
    }
    cat("  Snow-index representativeness summary:\n")
    for (k in seq_along(snow_specs)) {
      if (is.null(snow_r_list[[k]])) next
      cat(sprintf("    %-8s median r = %.3f; %.0f%% of pixels r ≥ 0.80\n",
                  snow_specs[[k]]$label, snow_medians[k], snow_pct80[k]))
    }
  }
}, error = function(e) {
  cat(sprintf("  ⚠ Snow-index spatial representativeness figure failed: %s\n", e$message))
})

################################################################################
# PART 4d – SEASONALITY FIGURE  (P vs PET(PM) — 12 monthly means)
#
# Two-panel bar chart:
#   Panel (a): mean monthly precipitation P vs PET(PM)  [mm month⁻¹],
#              exactly 12 bars per variable averaged over all 1950–2025 years.
#   Panel (b): climatic water balance P − PET(PM)  [mm month⁻¹].
#
# DATA
#   Precipitation: monthly_data_direct/total_precipitation_monthly.nc
#                  (m/day mean rate) → × days_in_month × 1000 = mm month⁻¹
#                  Same conversion as 1SPI_ERALand.R.
#   PET(PM):       monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv
#                  (mean_pet in mm day⁻¹) → × mean-days-per-month = mm month⁻¹
#                  Written by 2b_PET_ERALand.R.
#
# OUTPUT  FigS_Seasonality_P_PET.pdf + .png  in BASIN_PLOT_DIR
################################################################################
cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("PART 4d: Seasonality figure — P vs PET(PM) (12 monthly means)\n")
cat("═══════════════════════════════════════════════════════════════════\n")

tryCatch({
  
  PRECIP_NC_4d  <- file.path(WD_PATH, "monthly_data_direct",
                             "total_precipitation_monthly.nc")
  PET_SEAS_CSV  <- file.path(WD_PATH, "monthly_data_direct",
                             "ERA5Land_Nechako_PET_monthly_summary.csv")
  
  if (!file.exists(PRECIP_NC_4d))
    stop(sprintf("Precip NetCDF not found: %s", PRECIP_NC_4d))
  if (!file.exists(PET_SEAS_CSV)) {
    cat(sprintf("  ⚠ PET CSV not found: %s\n", PET_SEAS_CSV))
    cat("    Run 2b_PET_ERALand.R first, then re-run this script.\n")
    stop("pet_csv_missing_4d")
  }
  pet_raw_4d       <- data.table::fread(PET_SEAS_CSV, data.table = FALSE)
  pet_raw_4d$date  <- as.Date(paste0(pet_raw_4d$year, "-", 
                                     sprintf("%02d", pet_raw_4d$month), "-01"))
  COL_P_4d   <- "#4393c3"
  COL_PET_4d <- "#d7301f"
  ## ── Precipitation from NetCDF (mm month⁻¹) ──────────────────────────────
  ## Mirrors the conversion in 1SPI_ERALand.R exactly:
  cat("  Loading precipitation from NetCDF...\n")
  pr_r  <- terra::rast(PRECIP_NC_4d)
  pr_dates <- tryCatch({
    t <- terra::time(pr_r)
    if (!all(is.na(t))) as.Date(t)
    else {
      nc  <- ncdf4::nc_open(PRECIP_NC_4d)
      out <- if ("valid_time" %in% names(nc$var)) {
        tv <- ncdf4::ncvar_get(nc, "valid_time")
        as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC"))
      } else {
        seq(as.Date("1950-01-01"), by = "month", length.out = terra::nlyr(pr_r))
      }
      ncdf4::nc_close(nc)
      out
    }
  }, error = function(e)
    seq(as.Date("1950-01-01"), by = "month", length.out = terra::nlyr(pr_r)))
  
  # Use the same native ERA5-Land analysis grid as Scripts 1–3.
  if (!terra::compareGeom(pr_r[[1]], analysis_template, stopOnError = FALSE,
                          crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)) {
    pr_r <- terra::project(pr_r, analysis_template, method = "bilinear")
  }
  # Basin mask = 7-subbasin union, cells with coverage > 0 (shared rule).
  pr_r <- terra::mask(pr_r, S12_BASIN_BW$mask_raster, maskvalues = 0, updatevalue = NA)
  
  # days_in_month from actual dates
  # days_in_month from actual dates
  fom_4d  <- as.Date(format(pr_dates, "%Y-%m-01"))
  fnm_4d  <- seq(fom_4d[1], by = "month", length.out = length(pr_dates) + 1)[-1]
  dim_4d  <- as.integer(fnm_4d - fom_4d)
  
  # ── PRECIPITATION UNIT CHECK (mirrors 1SPI_ERALand.R / 3SPEI_ERALand.R) ──────
  # ERA5-Land "tp" from CDS:
  #   "monthly averaged reanalysis" = mean daily rate → m/day  (raw, needs *1000 then *days)
  #   "monthly aggregated"          = monthly total   → m/month (needs *1000 only)
  #   Already-converted             = mm/day          (needs *days only)
  #   Already-converted + totalled  = mm/month        (needs nothing)
  #
  # Decision thresholds on RAW July basin-mean:
  #   < 0.005   → m/day    → multiply by days_in_month AND 1000   (→ mm/month)
  #   0.005–20  → mm/day   → multiply by days_in_month ONLY        (→ mm/month)
  #   20–500    → mm/month → DO NOT multiply (already monthly total)
  #   > 500     → unknown unit error
  #
  # NOTE: the old threshold was ≤1 for mm/day, but ERA5-Land July basin means
  # can exceed 1 mm/day (Nechako July ≈ 1.6 mm/day).  Raised to ≤20 mm/day
  # to avoid mis-classifying valid mm/day values as mm/month.
  # ─────────────────────────────────────────────────────────────────────────────
  july_idx_4d <- which(as.integer(format(pr_dates, "%m")) == 7L)[1L]
  
  if (!is.na(july_idx_4d)) {
    raw_july_vals <- terra::values(pr_r[[july_idx_4d]], mat = FALSE)
    raw_july_vals <- raw_july_vals[is.finite(raw_july_vals)]
    raw_mean_p    <- mean(raw_july_vals, na.rm = TRUE)
    
    cat(sprintf(
      "\n── 7 Part 4d PRECIPITATION UNIT CHECK (July layer %d, BEFORE conversion) ──\n",
      july_idx_4d))
    cat(sprintf("  Raw Min  : %12.6f\n", min(raw_july_vals)))
    cat(sprintf("  Raw Mean : %12.6f\n", raw_mean_p))
    cat(sprintf("  Raw Max  : %12.6f\n", max(raw_july_vals)))
    cat("  ── Thresholds ──────────────────────────────────────────────────────\n")
    cat("  < 0.005        → m/day        → apply * days_in_month * 1000\n")
    cat("  0.005 to 20    → mm/day       → apply * days_in_month  (NO *1000)\n")
    cat("  20 to 500      → mm/month     → apply * 1              (NO further scaling)\n")
    cat("  > 500          → already large; check source carefully\n")
    cat("  ─────────────────────────────────────────────────────────────────────\n")
    
    if (raw_mean_p < 0.005) {
      precip_unit_4d  <- "m/day"
      precip_scale_4d <- dim_4d * 1000          # → mm/month
      cat(sprintf(
        "  ✓ DETECTED: m/day  →  applying * days_in_month * 1000 (→ mm/month)\n"))
      cat(sprintf(
        "  Expected July result: ~%.0f mm/month\n", raw_mean_p * 31 * 1000))
      
    } else if (raw_mean_p <= 20) {
      precip_unit_4d  <- "mm/day"
      precip_scale_4d <- dim_4d * 1            # → mm/month  (no *1000)
      cat(sprintf(
        "  ✓ DETECTED: mm/day →  applying * days_in_month ONLY (→ mm/month)\n"))
      cat(sprintf(
        "  Expected July result: ~%.0f mm/month\n", raw_mean_p * 31))
      cat("  ⚠  NOTE: old 7 line was  pr_r * dim_4d * 1000  — that would be WRONG here.\n")
      
    } else if (raw_mean_p <= 500) {
      precip_unit_4d  <- "mm/month"
      precip_scale_4d <- rep(1L, length(dim_4d))   # no scaling
      cat(sprintf(
        "  ✓ DETECTED: mm/month → NO further scaling applied\n"))
      cat("  ⚠  NOTE: applying * days_in_month * 1000 here would produce ~1000x inflated values.\n")
      
    } else {
      precip_unit_4d  <- "UNKNOWN (>500)"
      precip_scale_4d <- rep(1L, length(dim_4d))
      cat(sprintf(
        "  ✗ ERROR: raw mean = %.1f — unit unrecognised. No scaling applied.\n",
        raw_mean_p))
      cat("    Check your NetCDF source before continuing.\n")
      warning(sprintf(
        "[7 Part 4d] Precipitation unit unrecognised (July mean = %.1f). Manual check required.",
        raw_mean_p))
    }
    
    # Post-conversion sanity check threshold (July should be 20–200 mm/month for Nechako)
    expected_after <- raw_mean_p * mean(precip_scale_4d[as.integer(format(pr_dates, "%m")) == 7])
    cat(sprintf(
      "  Post-conversion July estimate: ~%.1f mm/month\n", expected_after))
    if (expected_after < 5 || expected_after > 500) {
      warning(sprintf(
        "[7 Part 4d] Post-conversion July precipitation (%.1f mm/month) is outside the\n",
        expected_after),
        "  expected range 5–500 mm/month for Nechako. Re-check unit conversion.")
      cat(sprintf(
        "  ✗ WARNING: %.1f mm/month is outside the expected range [5, 500] — review conversion!\n",
        expected_after))
    } else {
      cat(sprintf(
        "  ✓ Post-conversion value looks physically reasonable (5–500 mm/month range).\n"))
    }
    cat(sprintf("  Detected unit: %s\n", precip_unit_4d))
    cat("────────────────────────────────────────────────────────────────────────\n\n")
    
  } else {
    cat("  ⚠⚠⚠ Part 4d: Could not find a July layer — defaulting to * dim_4d * 1000.\n")
    precip_scale_4d <- dim_4d * 1000
  }
  
  # Apply the auto-detected conversion
  pr_r <- pr_r * precip_scale_4d   # → mm/month (scale factor chosen above)
  
  # Basin mean per layer: fractional-coverage area weights, 0.20 coverage rule
  # (same estimator as the PET(PM) basin CSV written by 2b_PET_ERALand.R).
  pr_vals <- basin_weighted_mean(terra::values(pr_r, mat = TRUE), S12_BASIN_BW)
  prec_df_4d <- data.frame(
    date  = pr_dates,
    value = as.numeric(pr_vals),
    month = as.integer(format(pr_dates, "%m")))
  
  cm_prec_4d <- prec_df_4d |>
    dplyr::group_by(month) |>
    dplyr::summarise(mean_val = mean(value, na.rm = TRUE),
                     sd_val   = sd(value,   na.rm = TRUE),
                     n_yrs    = dplyr::n(),
                     .groups  = "drop") |>
    dplyr::mutate(
      month_lab = factor(month.abb[month], levels = month.abb),
      variable  = "Precipitation (P)")
  
  ## ── PET(PM) from CSV (mm month⁻¹) ────────────────────────────────────────
  pet_raw_4d$days_in_month <- month_lengths(pet_raw_4d$date)
  pet_raw_4d$pet_mm <- pet_raw_4d$mean_pet * pet_raw_4d$days_in_month
  
  cm_pet_4d <- pet_raw_4d |>
    dplyr::group_by(month) |>
    dplyr::summarise(mean_val = mean(pet_mm, na.rm = TRUE),
                     sd_val   = sd(pet_mm,   na.rm = TRUE),
                     n_yrs    = dplyr::n(),
                     .groups  = "drop") |>
    dplyr::mutate(
      month_lab = factor(month.abb[month], levels = month.abb),
      variable  = "PET(PM)")
  
  ann_P   <- round(sum(cm_prec_4d$mean_val))
  ann_PET <- round(sum(cm_pet_4d$mean_val))
  ai      <- round(ann_PET / ann_P, 2)
  nyrs    <- cm_prec_4d$n_yrs[1]
  
  cat(sprintf("  Precip  : %d yrs, mean annual = %d mm yr⁻¹\n", nyrs, ann_P))
  cat(sprintf("  PET(PM) : %d yrs, mean annual = %d mm yr⁻¹\n",
              cm_pet_4d$n_yrs[1], ann_PET))
  cat(sprintf("  Aridity index (PET/P) = %.2f\n", ai))
  
  ## ── Panel (a): grouped bar chart P vs PET(PM) ────────────────────────────
  cm_all_4d          <- rbind(cm_prec_4d, cm_pet_4d)
  cm_all_4d$variable <- factor(cm_all_4d$variable,
                               levels = c("Precipitation (P)", "PET(PM)"))
  
  p4d_a <- ggplot2::ggplot(
    cm_all_4d,
    ggplot2::aes(x = month_lab, y = mean_val,
                 fill = variable, group = variable)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(0.72),
      width = 0.65, colour = "white", linewidth = 0.15) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
      position = ggplot2::position_dodge(0.72),
      width = 0.22, colour = "grey25", linewidth = 0.40) +
    ggplot2::scale_fill_manual(
      values = c("Precipitation (P)" = COL_P_4d, "PET(PM)" = COL_PET_4d),
      name   = NULL) +
    ggplot2::scale_y_continuous(
      name   = expression(mm~month^{-1}),
      expand = ggplot2::expansion(mult = c(0, 0.09))) +
    ggplot2::labs(
      title    = "(a)  Mean monthly precipitation and PET(PM)",
      subtitle = paste0(
        nyrs, "-year mean (1950–2025)  |  Error bars = ±1 SD  |  ",
        "Annual totals:  P = ", ann_P, " mm;  PET(PM) = ", ann_PET, " mm"),
      x = NULL) +
    ggplot2::theme_classic(base_size = 9.5) +
    ggplot2::theme(
      legend.position    = c(0.86, 0.88),
      legend.background  = ggplot2::element_rect(fill = "white", colour = NA),
      legend.text        = ggplot2::element_text(size = 8.5),
      legend.key.size    = ggplot2::unit(0.40, "cm"),
      plot.title         = ggplot2::element_text(size = 9.5, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 7.5, colour = "grey40"),
      axis.text          = ggplot2::element_text(size = 8.5),
      axis.title.y       = ggplot2::element_text(size = 8.5),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      plot.margin        = ggplot2::margin(4, 6, 2, 4))
  
  ## ── Panel (b): water balance P − PET(PM) ─────────────────────────────────
  cm_wb_4d <- dplyr::inner_join(
    cm_prec_4d[, c("month", "month_lab", "mean_val")],
    cm_pet_4d[,  c("month", "mean_val")],
    by = "month", suffix = c("_p", "_pet")) |>
    dplyr::mutate(wb       = mean_val_p - mean_val_pet,
                  fill_col = ifelse(wb >= 0, COL_P_4d, COL_PET_4d))
  
  # x positions of surplus and deficit groups (for annotations outside bars)
  surplus_x  <- mean(which(cm_wb_4d$wb >= 0))
  deficit_x  <- mean(which(cm_wb_4d$wb <  0))
  # y positions placed just outside the zero line so they never overlap a bar
  surplus_y  <-  max(cm_wb_4d$wb[cm_wb_4d$wb >= 0], na.rm = TRUE) * 1.18
  deficit_y  <- -abs(min(cm_wb_4d$wb, na.rm = TRUE)) * 1.14
  
  p4d_b <- ggplot2::ggplot(
    cm_wb_4d,
    ggplot2::aes(x = month_lab, y = wb, fill = fill_col)) +
    ggplot2::geom_col(colour = "white", linewidth = 0.15, width = 0.65) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.55) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(
      name   = expression(P - PET(PM)~~(mm~month^{-1})),
      expand = ggplot2::expansion(mult = c(0.16, 0.16))) +
    # Labels placed ABOVE the tallest surplus bar and BELOW the deepest deficit bar
    ggplot2::annotate("text",
                      x = surplus_x, y = surplus_y,
                      label = "P surplus", colour = COL_P_4d,
                      size = 3.0, fontface = "italic") +
    ggplot2::annotate("text",
                      x = deficit_x, y = deficit_y,
                      label = "PET deficit", colour = COL_PET_4d,
                      size = 3.0, fontface = "italic") +
    ggplot2::labs(
      title    = "(b)  Climatic water balance  P − PET(PM)",
      subtitle = "Blue = monthly surplus (P > PET);  red = monthly deficit (PET > P)",
      x        = "Calendar month") +
    ggplot2::theme_classic(base_size = 9.5) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(size = 9.5, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 7.5, colour = "grey40"),
      axis.text          = ggplot2::element_text(size = 8.5),
      axis.title.y       = ggplot2::element_text(size = 8.5),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      plot.margin        = ggplot2::margin(2, 6, 4, 4))
  
  ## ── Assemble and save ────────────────────────────────────────────────────
  fig_seas_4d <- p4d_a / p4d_b +
    patchwork::plot_annotation(
      title    = paste0(
        "Nechako River Basin — Seasonal cycle of precipitation ",
        "and PET(PM)"),
      subtitle = paste0(
        "Each bar = mean of ", nyrs, " years (1950–2025).  ",
        "Annual totals: P = ", ann_P, " mm;  PET(PM) = ", ann_PET, " mm;  ",
        "aridity index (PET/P) = ", ai, "."),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 11, face = "bold", hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey35",
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 4))))
  
  for (ext in c("pdf", "png")) {
    out_seas <- file.path(BASIN_PLOT_DIR,
                          paste0("FigS_Seasonality_P_PET.", ext))
    tryCatch(
      ggplot2::ggsave(out_seas, fig_seas_4d,
                      width  = 7.2,
                      height = 7.0,
                      units  = "in",
                      dpi    = if (ext == "png") 300 else "print",
                      device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png"),
      error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message)))
    cat(sprintf("  ✓ Saved: %s\n", basename(out_seas)))
  }
  
}, error = function(e) {
  if (!identical(conditionMessage(e), "pet_csv_missing_4d"))
    cat(sprintf("  ⚠ Part 4d failed: %s\n", e$message))
})

################################################################################
################################################################################
# PART 4e – MANUSCRIPT RESULTS FIGURES: FULL HISTORICAL DROUGHT CHRONOLOGY
#
# Purpose:
#   Provide the historical time-series evidence needed to compare the
#   2022–2025 drought with all prior drought episodes across the complete
#   1950–2025 record and across all six manuscript accumulation scales.
#
# Primary manuscript figure:
#   Fig2_MS_AllIndices_3x6.pdf/.png
#       18 panels = 3 index families × 6 accumulation scales
#       Row 1: SPI-1, SPI-3, SPI-6, SPI-12, SPI-24, SPI-36
#       Row 2: SPEI-PM-1, SPEI-PM-3, SPEI-PM-6, SPEI-PM-12, SPEI-PM-24, SPEI-PM-36
#       Row 3: SPEI-ERA5-1, SPEI-ERA5-3, SPEI-ERA5-6, SPEI-ERA5-12, SPEI-ERA5-24, SPEI-ERA5-36
#       Full 1950–2025 chronology; 2022–2025 is highlighted.
#
# Secondary focused figure:
#   Fig2S_MS_Recent_AllScales_3x2.pdf/.png
#       Six panels = accumulation scales 1, 3, 6, 12, 24, 36 months.
#       Each panel overlays SPI and SPEI-PM for 2022–2025. SPEI-ERA5 is
#       intentionally excluded from this manuscript-facing recent figure.
#
# Regional focused figure:
#   Fig2S_MS_Recent_Subbasins_7x6.pdf/.png
#       Seven rows = the seven supplied subbasins; six columns = accumulation
#       scales. Each panel overlays SPI and SPEI-PM for 2022–2025. The
#       whole-basin union and SPEI-ERA5 are intentionally excluded.
#
# Scientific interpretation:
#   These are descriptive basin-average chronologies. The 2022–2025 highlight
#   identifies the focus interval but does not alter the historical reference
#   period. No detrending, counterfactual reconstruction, or PET decomposition
#   is introduced here.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4e: Full historical drought chronology + recent comparison\n")
cat("══════════════════════════════════════════════════════════════\n")

# Manuscript scales are fixed explicitly here so that the figure cannot silently
# inherit a reduced utility-scale vector (e.g. 1, 2, 3, 6, 12) and omit the
# long accumulation periods required for historical persistence comparison.
MANUSCRIPT_TS_SCALES <- c(1L, 3L, 6L, 12L, 24L, 36L)

THR_EVENT_4e       <- -0.5
THR_MOD_4e         <- -1.0
HIGHLIGHT_START_4e <- as.Date("2022-01-01")
HIGHLIGHT_END_4e   <- as.Date("2025-12-31")
RECENT_START_4e    <- as.Date("2022-01-01")
RECENT_END_4e      <- as.Date("2025-12-01")

fill_pal_4e <- c(
  "Wet"              = "#4393c3",
  "Mild deficit"     = "#ffffbf",
  "Moderate drought" = "#fdae61",
  "Severe drought"   = "#f46d43",
  "Extreme drought"  = "#d73027"
)

LINE_COLOURS_4e <- c(
  "SPI"     = "#2ca02c",
  "SPEI-PM" = "#1f77b4",
  "SPEI-ERA5" = "#d62728"
)

# ── SPI per-month CSV fallback loader ─────────────────────────────────────────
load_spi_seasonal_csvs <- function(scale, seas_dir = SPI_SEAS_DIR) {
  if (!dir.exists(seas_dir)) {
    cat(sprintf("    [load_spi_seasonal_csvs] Directory not found: %s\n", seas_dir))
    return(NULL)
  }
  
  month_rows <- vector("list", 12L)
  
  for (m in seq_len(12L)) {
    csv_file <- file.path(
      seas_dir,
      sprintf("spi_%02d_month%02d_%s.csv", scale, m, MONTH_NAMES_3[m])
    )
    
    if (!file.exists(csv_file)) next
    
    tryCatch({
      df <- data.table::fread(csv_file, data.table = FALSE)
      yr_cols <- setdiff(names(df), c("lon", "lat"))
      yr_cols <- yr_cols[grepl("^[0-9]{4}$", yr_cols)]
      if (!length(yr_cols)) next
      
      mat_yr <- as.matrix(df[, yr_cols, drop = FALSE])
      is_geo <- max(abs(df$lon), na.rm = TRUE) <= 200
      
      if (is_geo) {
        w <- cos(df$lat * pi / 180)
        vals <- apply(mat_yr, 2, function(x) {
          ok <- is.finite(x) & is.finite(w)
          if (!any(ok)) return(NA_real_)
          sum(x[ok] * w[ok]) / sum(w[ok])
        })
      } else {
        vals <- apply(mat_yr, 2, function(x) {
          ok <- is.finite(x)
          if (!any(ok)) return(NA_real_)
          mean(x[ok])
        })
      }
      
      vals[is.nan(vals)] <- NA_real_
      
      month_rows[[m]] <- data.frame(
        date  = as.Date(sprintf("%s-%02d-01", yr_cols, m)),
        value = as.numeric(vals),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat(sprintf("    ⚠ Could not read %s: %s\n", basename(csv_file), e$message))
    })
  }
  
  valid <- Filter(Negate(is.null), month_rows)
  if (!length(valid)) return(NULL)
  
  out <- do.call(rbind, valid)
  out <- out[order(out$date), , drop = FALSE]
  out <- out[is.finite(out$value), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ── Unified historical time-series loader ─────────────────────────────────────
# IMPORTANT: this function only reads existing products. It does not recompute
# spatial means from the drought-index NetCDF files unless a required scale is
# absent from the authoritative script-6 CSVs, in which case it uses the same
# per-calendar-month fallback already used elsewhere in this script.
load_manuscript_ts_4e <- function(index_type, scale) {
  index_type <- match.arg(index_type, c("spi", "spei_pm", "spei_era5"))
  
  if (index_type == "spi") {
    df <- tryCatch(load_flat_csv("spi", scale), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) {
      df <- tryCatch(load_spi_seasonal_csvs(scale, SPI_SEAS_DIR),
                     error = function(e) NULL)
    }
  } else if (index_type == "spei_pm") {
    df <- tryCatch(load_flat_csv("spei", scale), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) {
      df <- tryCatch(load_spei_seasonal_csvs(scale, SPEI_PM_DIR),
                     error = function(e) NULL)
    }
  } else {
    # Script 11 is authoritative for basin averages. Read its exact CSV path;
    # do not silently recompute the basin mean from pixel-level seasonal files.
    df <- tryCatch(load_any_ts("spei_era5", scale), error = function(e) NULL)
  }
  
  if (is.null(df) || !nrow(df)) return(NULL)
  
  df$date <- as.Date(df$date)
  df$value <- as.numeric(df$value)
  df <- df[is.finite(df$value) & !is.na(df$date), c("date", "value"), drop = FALSE]
  df <- df[order(df$date), , drop = FALSE]
  rownames(df) <- NULL
  df
}

# ── Historical panel builder: one index × one accumulation scale ──────────────
make_historical_panel_4e <- function(sc, panel_lab, index_type,
                                     show_xlab = FALSE, show_ylabel = TRUE,
                                     is_first = FALSE) {
  df <- load_manuscript_ts_4e(index_type, sc)
  
  index_label <- switch(
    index_type,
    spi      = sprintf("SPI-%d", sc),
    spei_pm  = sprintf("SPEI-PM-%d", sc),
    spei_era5 = sprintf("SPEI-ERA5-%d", sc)
  )
  
  line_col <- unname(LINE_COLOURS_4e[
    switch(index_type, spi = "SPI", spei_pm = "SPEI-PM", spei_era5 = "SPEI-ERA5")
  ])
  
  if (is.null(df) || !nrow(df)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = sprintf("No data\n%s", index_label),
          size = 3.2, colour = "grey50"
        ) +
        ggplot2::theme_void()
    )
  }
  
  # Retain the full historical record. The 2022–2025 band is an annotation,
  # not a data filter.
  df_fill <- df
  df_fill <- df_fill[is.finite(df_fill$value), , drop = FALSE]
  df_fill$fill_cat <- factor(
    dplyr::case_when(
      df_fill$value > 0    ~ "Wet",
      df_fill$value > -0.5 ~ "Mild deficit",
      df_fill$value > -1.0 ~ "Moderate drought",
      df_fill$value > -1.5 ~ "Severe drought",
      TRUE                 ~ "Extreme drought"
    ),
    levels = names(fill_pal_4e)
  )
  df_fill$date_end <- c(df_fill$date[-1L], next_month_start(df_fill$date))
  
  ylo <- min(-3.2, floor(min(df$value, na.rm = TRUE) * 1.08 * 2) / 2)
  yhi <- max( 2.5, ceiling(max(df$value, na.rm = TRUE) * 1.08 * 2) / 2)
  
  p <- ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = HIGHLIGHT_START_4e,
      xmax = HIGHLIGHT_END_4e,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey85",
      alpha = 0.55
    ) +
    ggplot2::geom_rect(
      data = df_fill,
      ggplot2::aes(
        xmin = date,
        xmax = date_end,
        ymin = -Inf,
        ymax = Inf,
        fill = fill_cat
      ),
      alpha = 0.25,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = fill_pal_4e,
      name = "Drought category",
      drop = FALSE,
      guide = if (is_first) {
        ggplot2::guide_legend(
          title.position = "top",
          nrow = 1,
          override.aes = list(alpha = 0.65)
        )
      } else {
        "none"
      }
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey40",
      linewidth = 0.35
    ) +
    ggplot2::geom_hline(
      yintercept = THR_EVENT_4e,
      colour = "grey25",
      linewidth = 0.50,
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = THR_MOD_4e,
      colour = "grey50",
      linewidth = 0.30,
      linetype = "dotted"
    ) +
    ggplot2::geom_line(
      data = df,
      ggplot2::aes(x = date, y = value),
      colour = line_col,
      linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = min(df$date) + 365,
      y = yhi * 0.90,
      label = panel_lab,
      size = 3.0,
      fontface = "bold",
      hjust = 0,
      vjust = 1,
      colour = "grey10"
    ) +
    ggplot2::annotate(
      "text",
      x = as.Date("1951-01-01"),
      y = THR_EVENT_4e + 0.08,
      label = "−0.5",
      size = 2.2,
      hjust = 0,
      colour = "grey30"
    ) +
    ggplot2::annotate(
      "text",
      x = as.Date("2023-12-01"),
      y = yhi * 0.97,
      label = "2022–2025",
      size = 2.1,
      hjust = 0.5,
      colour = "grey30",
      fontface = "italic"
    ) +
    ggplot2::scale_x_date(
      date_breaks = "10 years",
      date_labels = "%Y",
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(ylo, yhi),
      breaks = seq(floor(ylo), ceiling(yhi), by = 1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = index_label,
      x = if (show_xlab) "Year" else NULL,
      y = if (show_ylabel) "Standardized index" else NULL
    ) +
    ggplot2::theme_classic(base_size = 8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 8.2,
        face = "bold",
        hjust = 0,
        margin = ggplot2::margin(b = 2)
      ),
      axis.text.x = if (show_xlab) {
        ggplot2::element_text(size = 6.5)
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (show_xlab) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = ggplot2::element_text(size = 6.3),
      axis.title.y = if (show_ylabel) {
        ggplot2::element_text(size = 7)
      } else {
        ggplot2::element_blank()
      },
      axis.title.x = ggplot2::element_text(size = 7),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 6.2),
      legend.title = ggplot2::element_text(size = 6.5),
      legend.key.size = ggplot2::unit(0.30, "cm"),
      plot.margin = ggplot2::margin(2, 4, 2, 3)
    )
  
  p
}

# ── Save the 18-panel historical chronology ───────────────────────────────────
save_historical_18panel_4e <- function(panels) {
  fig <- patchwork::wrap_plots(
    panels,
    nrow = 3,
    ncol = 6,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = paste0(
        "Basin-averaged drought-index chronology — Nechako River Basin (1950–2025)"
      ),
      subtitle = paste0(
        "Rows: SPI, SPEI-PM, SPEI-ERA5  |  Columns: 1-, 3-, 6-, 12-, 24-, 36-month accumulation scales  |  ",
        "Dashed line: drought-event threshold (x₀ = −0.5)  |  Grey band: 2022–2025"
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 11,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(b = 3)
        ),
        plot.subtitle = ggplot2::element_text(
          size = 7.5,
          colour = "grey40",
          hjust = 0.5,
          margin = ggplot2::margin(b = 5)
        ),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 7),
        legend.title = ggplot2::element_text(size = 7.2),
        legend.key.size = ggplot2::unit(0.35, "cm")
      )
    )
  
  for (ext in c("pdf", "png")) {
    out <- file.path(
      BASIN_PLOT_DIR,
      paste0("Fig2_MS_AllIndices_3x6.", ext)
    )
    tryCatch(
      ggplot2::ggsave(
        out,
        fig,
        width = 16,
        height = 9.5,
        units = "in",
        dpi = if (ext == "png") 300 else "print",
        device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png",
        limitsize = FALSE
      ),
      error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message))
    )
    cat(sprintf("  ✓ Saved: %s\n", basename(out)))
  }
  
  invisible(fig)
}

# ── Build historical 18-panel figure ──────────────────────────────────────────
tryCatch({
  cat("  Building 18-panel historical chronology figure...\n")
  
  index_grid <- rep(c("spi", "spei_pm", "spei_era5"), each = length(MANUSCRIPT_TS_SCALES))
  scale_grid <- rep(MANUSCRIPT_TS_SCALES, times = 3L)
  label_grid <- sprintf("(%s)", letters[seq_along(index_grid)])
  
  panels_18 <- lapply(seq_along(index_grid), function(k) {
    row_idx <- ceiling(k / length(MANUSCRIPT_TS_SCALES))
    col_idx <- ((k - 1L) %% length(MANUSCRIPT_TS_SCALES)) + 1L
    
    make_historical_panel_4e(
      sc = scale_grid[k],
      panel_lab = label_grid[k],
      index_type = index_grid[k],
      show_xlab = row_idx == 3L,
      show_ylabel = col_idx == 1L,
      is_first = k == 1L
    )
  })
  
  save_historical_18panel_4e(panels_18)
  cat("  ✓ 18-panel historical chronology complete\n")
}, error = function(e) {
  cat(sprintf("  ⚠ Historical 18-panel figure failed: %s\n", e$message))
})

# ── Recent six-panel formulation comparison ───────────────────────────────────
# Each panel is one accumulation scale. This figure is deliberately restricted
# to 2022–2025 and overlays the three independently standardized drought
# formulations. It complements, rather than replaces, the full historical
# chronology above.
make_recent_scale_panel_4e <- function(sc, panel_lab, show_xlab = FALSE,
                                       show_ylabel = FALSE, is_first = FALSE) {
  load_one <- function(index_type) {
    df <- load_manuscript_ts_4e(index_type, sc)
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- df[
      df$date >= RECENT_START_4e & df$date <= RECENT_END_4e,
      c("date", "value"),
      drop = FALSE
    ]
    if (!nrow(df)) return(NULL)
    df$index_type <- index_type
    df
  }
  
  # Manuscript-facing recent comparison deliberately contains only SPI and
  # SPEI-PM. SPEI-ERA5 remains available in the historical/diagnostic products
  # but is not plotted in this figure.
  df_all <- dplyr::bind_rows(
    load_one("spi"),
    load_one("spei_pm")
  )
  
  if (!nrow(df_all)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = sprintf("No data\nscale %d", sc),
          size = 3.2, colour = "grey50"
        ) +
        ggplot2::theme_void()
    )
  }
  
  ylo <- min(-3.2, floor(min(df_all$value, na.rm = TRUE) * 1.08 * 2) / 2)
  yhi <- max( 2.5, ceiling(max(df_all$value, na.rm = TRUE) * 1.08 * 2) / 2)
  
  df_all$label <- factor(
    df_all$index_type,
    levels = c("spi", "spei_pm"),
    labels = c("SPI", "SPEI-PM")
  )
  
  p <- ggplot2::ggplot(
    df_all,
    ggplot2::aes(x = date, y = value, colour = label)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey40",
      linewidth = 0.35
    ) +
    ggplot2::geom_hline(
      yintercept = THR_EVENT_4e,
      colour = "grey25",
      linewidth = 0.50,
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = THR_MOD_4e,
      colour = "grey50",
      linewidth = 0.30,
      linetype = "dotted"
    ) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::scale_colour_manual(
      values = LINE_COLOURS_4e,
      name = "Index",
      guide = if (is_first) {
        ggplot2::guide_legend(
          title.position = "top",
          nrow = 1,
          override.aes = list(linewidth = 1.0)
        )
      } else {
        "none"
      }
    ) +
    ggplot2::annotate(
      "text",
      x = as.Date("2022-02-01"),
      y = yhi * 0.92,
      label = panel_lab,
      size = 3.0,
      fontface = "bold",
      hjust = 0,
      colour = "grey10"
    ) +
    ggplot2::annotate(
      "text",
      x = as.Date("2022-02-01"),
      y = THR_EVENT_4e + 0.08,
      label = "−0.5",
      size = 2.1,
      hjust = 0,
      colour = "grey30"
    ) +
    ggplot2::scale_x_date(
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = ggplot2::expansion(add = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(ylo, yhi),
      breaks = seq(floor(ylo), ceiling(yhi), by = 1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = sprintf("%d-month accumulation", sc),
      x = if (show_xlab) "Year" else NULL,
      y = if (show_ylabel) "Standardized index" else NULL
    ) +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 8.5,
        face = "bold",
        hjust = 0,
        margin = ggplot2::margin(b = 2)
      ),
      axis.text.x = if (show_xlab) {
        ggplot2::element_text(size = 6.5)
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (show_xlab) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = ggplot2::element_text(size = 6.3),
      axis.title.y = if (show_ylabel) {
        ggplot2::element_text(size = 7)
      } else {
        ggplot2::element_blank()
      },
      axis.title.x = ggplot2::element_text(size = 7),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 6.5),
      legend.title = ggplot2::element_text(size = 7),
      legend.key.size = ggplot2::unit(0.35, "cm"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 4, 2, 3)
    )
  
  p
}

tryCatch({
  cat("  Building six-panel 2022–2025 formulation-comparison figure...\n")
  
  recent_panels <- lapply(seq_along(MANUSCRIPT_TS_SCALES), function(k) {
    sc <- MANUSCRIPT_TS_SCALES[k]
    make_recent_scale_panel_4e(
      sc = sc,
      panel_lab = sprintf("(%s)", letters[k]),
      show_xlab = TRUE,
      show_ylabel = k == 1L,
      is_first = k == 1L
    )
  })
  
  fig_recent <- patchwork::wrap_plots(
    recent_panels,
    nrow = 2,
    ncol = 3,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = "2022–2025 drought evolution across accumulation scales",
      subtitle = paste0(
        "Each panel overlays SPI and SPEI-PM.  ",
        "SPEI-ERA5 is excluded. Dashed line: drought-event threshold (x₀ = −0.5); ",
        "dotted line: −1.0."
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 11,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(b = 3)
        ),
        plot.subtitle = ggplot2::element_text(
          size = 7.5,
          colour = "grey40",
          hjust = 0.5,
          margin = ggplot2::margin(b = 5)
        ),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 7),
        legend.title = ggplot2::element_text(size = 7.2),
        legend.key.size = ggplot2::unit(0.35, "cm")
      )
    )
  
  for (ext in c("pdf", "png")) {
    out <- file.path(
      BASIN_PLOT_DIR,
      paste0("Fig2S_MS_Recent_AllScales_3x2.", ext)
    )
    tryCatch(
      ggplot2::ggsave(
        out,
        fig_recent,
        width = 12,
        height = 7.5,
        units = "in",
        dpi = if (ext == "png") 300 else "print",
        device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png",
        limitsize = FALSE
      ),
      error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message))
    )
    cat(sprintf("  ✓ Saved: %s\n", basename(out)))
  }
  
  cat("  ✓ Six-panel recent comparison complete\n")
}, error = function(e) {
  cat(sprintf("  ⚠ Recent six-panel figure failed: %s\n", e$message))
})

cat("\n══ PART 4e complete ══\n")

################################################################################
# PART 4e-B — REDUCED 12-PANEL HISTORICAL DROUGHT CHRONOLOGY
#              SCALES 3, 12, 24, AND 36 MONTHS ONLY
#
# Manuscript layout:
#
#                         3 mo          12 mo          24 mo          36 mo
#              SPI        (a)            (b)            (c)            (d)
#              SPEI-PM    (e)            (f)            (g)            (h)
#              SPEI-ERA5   (i)            (j)            (k)            (l)
#
# Implementation:
#   The shared column headers and row labels are placed directly inside the
#   ggplot panels rather than in separate patchwork/grob rows. This avoids the
#   large vertical whitespace produced by text-only layout cells.
#
# Output:
#   Fig2_MS_AllIndices_3x4_Scales3_12_24_36.pdf
#   Fig2_MS_AllIndices_3x4_Scales3_12_24_36.png
################################################################################

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("PART 4e-B: Reduced 12-panel historical drought chronology\n")
cat("          Manuscript layout: 3 index families × 4 accumulation scales\n")
cat("══════════════════════════════════════════════════════════════════\n")

REDUCED_HIST_SCALES <- c(3L, 12L, 24L, 36L)
REDUCED_HIST_INDEXES <- c("spi", "spei_pm", "spei_era5")
REDUCED_HIST_INDEX_LABELS <- c("SPI", "SPEI-PM", "SPEI-ERA5")
REDUCED_HIST_SCALE_LABELS <- paste0(REDUCED_HIST_SCALES, " mo")
REDUCED_HIST_PANEL_LETTERS <- letters[seq_len(12L)]

# Validate that all 12 required historical series are available.
reduced_panel_order <- unlist(
  lapply(
    c("SPI", "SPEI_PM", "SPEI_ERA5"),
    function(idx) paste0(idx, REDUCED_HIST_SCALES)
  ),
  use.names = FALSE
)

reduced_historical_data <- all_data[reduced_panel_order[reduced_panel_order %in% names(all_data)]]

if (!length(reduced_historical_data)) {
  stop(
    "No data available for the reduced historical chronology. ",
    "Expected SPI, SPEI-PM, and SPEI-ERA5 at scales ",
    paste(REDUCED_HIST_SCALES, collapse = ", "), "."
  )
}

missing_reduced <- setdiff(reduced_panel_order, names(reduced_historical_data))
if (length(missing_reduced)) {
  stop(
    "Reduced historical chronology is missing required panel(s): ",
    paste(missing_reduced, collapse = ", ")
  )
}

reduced_historical_data <- all_data[reduced_panel_order]

# Build the 12 panels in row-major order.
#
# IMPORTANT:
#   * Column headers are implemented as the titles of the four top-row plots.
#   * SPI/SPEI-PM/SPEI-ERA5 are implemented as the vertical y-axis titles of
#     the first-column plots only.
#   * "Standardized index" is therefore removed from all 12 panels.
#   * Panel letters remain inside the plotting region.
#
# This avoids separate header/row-label grobs and therefore removes the
# excessive blank space between labels and the actual plotting regions.

reduced_panels <- lapply(seq_along(reduced_panel_order), function(i) {
  
  row_idx <- ceiling(i / length(REDUCED_HIST_SCALES))
  col_idx <- ((i - 1L) %% length(REDUCED_HIST_SCALES)) + 1L
  
  index_type <- REDUCED_HIST_INDEXES[row_idx]
  scale <- REDUCED_HIST_SCALES[col_idx]
  
  panel_letter <- sprintf("(%s)", REDUCED_HIST_PANEL_LETTERS[i])
  
  # Only the first row receives the accumulation-scale header.
  column_header <- if (row_idx == 1L) {
    REDUCED_HIST_SCALE_LABELS[col_idx]
  } else {
    NULL
  }
  
  # Only the first column receives the vertical drought-index label.
  # This replaces the original "Standardized index" y-axis title.
  row_label <- if (col_idx == 1L) {
    REDUCED_HIST_INDEX_LABELS[row_idx]
  } else {
    NULL
  }
  
  p <- make_historical_panel_4e(
    sc = scale,
    panel_lab = panel_letter,
    index_type = index_type,
    show_xlab = row_idx == 3L,
    show_ylabel = FALSE,
    is_first = i == 1L
  )
  
  p +
    ggplot2::labs(
      title = column_header,
      y = row_label
    ) +
    ggplot2::theme(
      # Shared scale headers: large, bold, and immediately adjacent to
      # the plotting region.
      plot.title = if (row_idx == 1L) {
        ggplot2::element_text(
          size = 15,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(
            t = 0,
            r = 0,
            b = 0,
            l = 0
          )
        )
      } else {
        ggplot2::element_blank()
      },
      
      # Vertical drought-index labels only in the first column.
      axis.title.y = if (col_idx == 1L) {
        ggplot2::element_text(
          size = 15,
          face = "bold",
          angle = 90,
          margin = ggplot2::margin(
            t = 0,
            r = 2,
            b = 0,
            l = 0
          )
        )
      } else {
        ggplot2::element_blank()
      },
      
      # Keep numeric y-axis tick labels but remove their textual axis title
      # from every panel except the first-column index-family labels above.
      axis.text.y = ggplot2::element_text(size = 7),
      axis.ticks.y = ggplot2::element_line(),
      
      # Compact x-axis treatment: only the bottom row shows years.
      axis.text.x = if (row_idx == 3L) {
        ggplot2::element_text(size = 7)
      } else {
        ggplot2::element_blank()
      },
      
      axis.ticks.x = if (row_idx == 3L) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      },
      
      axis.title.x = if (row_idx == 3L) {
        ggplot2::element_text(
          size = 8,
          margin = ggplot2::margin(t = 2)
        )
      } else {
        ggplot2::element_blank()
      },
      
      # No artificial vertical padding above/below the plotting region.
      plot.margin = ggplot2::margin(
        t = if (row_idx == 1L) 0 else 0,
        r = 2,
        b = if (row_idx == 3L) 1 else 0,
        l = 1
      ),
      
      # Readable manuscript legend; the guide is collected from the first
      # panel only by make_historical_panel_4e().
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10.5),
      legend.key.size = ggplot2::unit(0.60, "cm")
    )
})

# Direct 4 × 3 patchwork. No separate header row, no separate row-label
# column, and therefore no large whitespace between labels and plots.
reduced_historical_figure <- patchwork::wrap_plots(
  reduced_panels,
  ncol = 4,
  nrow = 3,
  guides = "collect"
) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10.5),
      legend.key.size = ggplot2::unit(0.60, "cm"),
      legend.spacing.x = ggplot2::unit(0.45, "cm"),
      legend.margin = ggplot2::margin(
        t = 2,
        r = 2,
        b = 1,
        l = 2
      ),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
  )

reduced_pdf <- file.path(
  BASIN_PLOT_DIR,
  "Fig2_MS_AllIndices_3x4_Scales3_12_24_36.pdf"
)

reduced_png <- file.path(
  BASIN_PLOT_DIR,
  "Fig2_MS_AllIndices_3x4_Scales3_12_24_36.png"
)

ggplot2::ggsave(
  reduced_pdf,
  reduced_historical_figure,
  width = 14,
  height = 9.5,
  units = "in",
  device = grDevices::cairo_pdf,
  limitsize = FALSE
)

ggplot2::ggsave(
  reduced_png,
  reduced_historical_figure,
  width = 14,
  height = 9.5,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

cat(
  sprintf(
    "  ✓ Reduced 12-panel historical chronology saved:\n    %s\n    %s\n",
    reduced_pdf,
    reduced_png
  )
)



################################################################################
# PART 4e-C — ADDITIONAL 9-PANEL HISTORICAL DROUGHT CHRONOLOGIES
#
# Purpose:
#   Retain the existing 12-panel figure (3, 12, 24, 36 months) and additionally
#   produce two matched 9-panel figures separating short- and long-term scales:
#
#   1. Fig2_MS_AllIndices_3x3_Scales1_3_6
#      Columns = 1, 3, 6 months
#
#   2. Fig2_MS_AllIndices_3x3_Scales12_24_36
#      Columns = 12, 24, 36 months
#
# Layout for both:
#              scale 1       scale 3       scale 6
#   SPI        (a)           (b)           (c)
#   SPEI-PM    (d)           (e)           (f)
#   SPEI-ERA5   (g)           (h)           (i)
#
#   The long-term figure uses the same layout with scales 12, 24, and 36 months.
#   The plotting function, drought classes, thresholds, 2022–2025 highlight,
#   colours, and historical reference period are inherited unchanged from the
#   existing make_historical_panel_4e() implementation.
#
# Output:
#   Fig2_MS_AllIndices_3x3_Scales1_3_6.pdf/.png
#   Fig2_MS_AllIndices_3x3_Scales12_24_36.pdf/.png
################################################################################

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("PART 4e-C: Additional 9-panel historical drought chronologies\n")
cat("          Short-term (1, 3, 6 months) + long-term (12, 24, 36 months)\n")
cat("══════════════════════════════════════════════════════════════════\n")

build_9panel_historical_figure_4e <- function(scales, file_stub) {
  scales <- as.integer(scales)
  if (length(scales) != 3L || any(!scales %in% MANUSCRIPT_TS_SCALES)) {
    stop(
      "build_9panel_historical_figure_4e(): scales must contain exactly three " ,
      "values from MANUSCRIPT_TS_SCALES."
    )
  }
  
  index_keys <- c("spi", "spei_pm", "spei_era5")
  index_labels <- c("SPI", "SPEI-PM", "SPEI-ERA5")
  scale_labels <- paste0(scales, " mo")
  panel_letters <- letters[seq_len(9L)]
  
  panel_order <- unlist(
    lapply(
      c("SPI", "SPEI_PM", "SPEI_ERA5"),
      function(idx) paste0(idx, scales)
    ),
    use.names = FALSE
  )
  
  required_names <- c(
    paste0("SPI", scales),
    paste0("SPEI_PM", scales),
    paste0("SPEI_ERA5", scales)
  )
  
  missing_panels <- setdiff(required_names, names(all_data))
  if (length(missing_panels)) {
    stop(
      "9-panel historical chronology is missing required panel(s): ",
      paste(missing_panels, collapse = ", ")
    )
  }
  
  panels <- lapply(seq_along(panel_order), function(i) {
    row_idx <- ceiling(i / 3L)
    col_idx <- ((i - 1L) %% 3L) + 1L
    
    index_type <- index_keys[row_idx]
    scale <- scales[col_idx]
    
    column_header <- if (row_idx == 1L) scale_labels[col_idx] else NULL
    row_label <- if (col_idx == 1L) index_labels[row_idx] else NULL
    
    p <- make_historical_panel_4e(
      sc = scale,
      panel_lab = sprintf("(%s)", panel_letters[i]),
      index_type = index_type,
      show_xlab = row_idx == 3L,
      show_ylabel = FALSE,
      is_first = i == 1L
    )
    
    p +
      ggplot2::labs(
        title = column_header,
        y = row_label
      ) +
      ggplot2::theme(
        plot.title = if (row_idx == 1L) {
          ggplot2::element_text(
            size = 15,
            face = "bold",
            hjust = 0.5,
            margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
          )
        } else {
          ggplot2::element_blank()
        },
        axis.title.y = if (col_idx == 1L) {
          ggplot2::element_text(
            size = 15,
            face = "bold",
            angle = 90,
            margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0)
          )
        } else {
          ggplot2::element_blank()
        },
        axis.text.y = ggplot2::element_text(size = 8),
        axis.ticks.y = ggplot2::element_line(),
        axis.text.x = if (row_idx == 3L) {
          ggplot2::element_text(size = 8)
        } else {
          ggplot2::element_blank()
        },
        axis.ticks.x = if (row_idx == 3L) {
          ggplot2::element_line()
        } else {
          ggplot2::element_blank()
        },
        axis.title.x = if (row_idx == 3L) {
          ggplot2::element_text(
            size = 9,
            margin = ggplot2::margin(t = 2)
          )
        } else {
          ggplot2::element_blank()
        },
        plot.margin = ggplot2::margin(
          t = 0,
          r = 2,
          b = if (row_idx == 3L) 1 else 0,
          l = 1
        ),
        legend.text = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 10.5),
        legend.key.size = ggplot2::unit(0.60, "cm")
      )
  })
  
  figure <- patchwork::wrap_plots(
    panels,
    ncol = 3,
    nrow = 3,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 10.5),
        legend.key.size = ggplot2::unit(0.60, "cm"),
        legend.spacing.x = ggplot2::unit(0.45, "cm"),
        legend.margin = ggplot2::margin(t = 2, r = 2, b = 1, l = 2),
        plot.margin = ggplot2::margin(0, 0, 0, 0)
      )
    )
  
  pdf_file <- file.path(BASIN_PLOT_DIR, paste0(file_stub, ".pdf"))
  png_file <- file.path(BASIN_PLOT_DIR, paste0(file_stub, ".png"))
  
  ggplot2::ggsave(
    pdf_file,
    figure,
    width = 12,
    height = 9.5,
    units = "in",
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    png_file,
    figure,
    width = 12,
    height = 9.5,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  cat(sprintf(
    "  ✓ Saved 9-panel figure (%s):\n    %s\n    %s\n",
    paste(scales, collapse = ", "),
    pdf_file,
    png_file
  ))
  
  invisible(figure)
}

build_9panel_historical_figure_4e(
  scales = c(1L, 3L, 6L),
  file_stub = "Fig2_MS_AllIndices_3x3_Scales1_3_6"
)

build_9panel_historical_figure_4e(
  scales = c(12L, 24L, 36L),
  file_stub = "Fig2_MS_AllIndices_3x3_Scales12_24_36"
)

cat("  ✓ Both additional 9-panel historical chronologies complete\n")


################################################################################
# PART 4e-C2 — MANUSCRIPT FIGURE: SPI + SPEI-PM ACROSS ALL SIX SCALES
#
# Purpose:
#   Produce the manuscript-focused 12-panel chronology requested for the
#   Results section. It contains only the two drought formulations needed for
#   the precipitation-versus-climatic-water-balance comparison:
#
#       SPI     : 1, 3, 6, 12, 24, 36 months
#       SPEI-PM : 1, 3, 6, 12, 24, 36 months
#
# Layout: 4 rows × 3 columns (12 panels)
#   Row 1: SPI-1, SPI-3, SPI-6
#   Row 2: SPI-12, SPI-24, SPI-36
#   Row 3: SPEI-PM-1, SPEI-PM-3, SPEI-PM-6
#   Row 4: SPEI-PM-12, SPEI-PM-24, SPEI-PM-36
#
# This preserves all six manuscript accumulation scales while keeping the
# figure restricted to SPI and SPEI-PM. The plotting layers, drought-category
# shading, thresholds, 2022–2025 highlighting, and historical reference period
# are inherited directly from make_historical_panel_4e().
#
# Output:
#   Fig3_MS_SPI_SPEIPM_AllScales_4x3.pdf
#   Fig3_MS_SPI_SPEIPM_AllScales_4x3.png
################################################################################

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("PART 4e-C2: Manuscript 12-panel SPI + SPEI-PM chronology\n")
cat("            All six accumulation scales; 4 rows × 3 columns\n")
cat("══════════════════════════════════════════════════════════════════\n")

FIG3_HIST_SCALES <- MANUSCRIPT_TS_SCALES
FIG3_HIST_INDEXES <- c("spi", "spei_pm")
FIG3_HIST_INDEX_LABELS <- c("SPI", "SPEI-PM")
FIG3_HIST_PANEL_LETTERS <- letters[seq_len(12L)]

# Row-major panel specification. The six scales are split across two rows for
# each index so that all six accumulation periods are retained in a compact
# 4 × 3 manuscript layout.
FIG3_HIST_PANEL_SPEC <- data.frame(
  index_type = c(
    rep("spi", 6L),
    rep("spei_pm", 6L)
  ),
  scale = c(
    FIG3_HIST_SCALES,
    FIG3_HIST_SCALES
  ),
  stringsAsFactors = FALSE
)

FIG3_REQUIRED_NAMES <- paste0(
  ifelse(FIG3_HIST_PANEL_SPEC$index_type == "spi", "SPI", "SPEI_PM"),
  FIG3_HIST_PANEL_SPEC$scale
)

FIG3_MISSING <- setdiff(FIG3_REQUIRED_NAMES, names(all_data))
if (length(FIG3_MISSING)) {
  stop(
    "Fig. 3 SPI + SPEI-PM chronology is missing required panel(s): ",
    paste(FIG3_MISSING, collapse = ", ")
  )
}

FIG3_PANELS <- lapply(seq_len(nrow(FIG3_HIST_PANEL_SPEC)), function(i) {
  row_idx <- ceiling(i / 3L)
  col_idx <- ((i - 1L) %% 3L) + 1L
  index_type <- FIG3_HIST_PANEL_SPEC$index_type[i]
  scale <- FIG3_HIST_PANEL_SPEC$scale[i]
  
  # Each panel carries its accumulation-scale header because the six scales
  # occupy two rows for each drought formulation in the 4 × 3 layout.
  column_header <- paste0(scale, " mo")
  
  # Index-family labels appear once, at the first column of each two-row block.
  row_label <- if (col_idx == 1L && row_idx %in% c(1L, 3L)) {
    FIG3_HIST_INDEX_LABELS[match(index_type, FIG3_HIST_INDEXES)]
  } else {
    NULL
  }
  
  p <- make_historical_panel_4e(
    sc = scale,
    panel_lab = sprintf("(%s)", FIG3_HIST_PANEL_LETTERS[i]),
    index_type = index_type,
    show_xlab = row_idx == 4L,
    show_ylabel = FALSE,
    is_first = i == 1L
  )
  
  p +
    ggplot2::labs(
      title = column_header,
      y = row_label
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 13,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
      ),
      axis.title.y = if (col_idx == 1L && row_idx %in% c(1L, 3L)) {
        ggplot2::element_text(
          size = 13,
          face = "bold",
          angle = 90,
          margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0)
        )
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = ggplot2::element_text(size = 7),
      axis.ticks.y = ggplot2::element_line(),
      axis.text.x = if (row_idx == 4L) {
        ggplot2::element_text(size = 7)
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (row_idx == 4L) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      },
      axis.title.x = if (row_idx == 4L) {
        ggplot2::element_text(size = 8, margin = ggplot2::margin(t = 2))
      } else {
        ggplot2::element_blank()
      },
      plot.margin = ggplot2::margin(
        t = 0,
        r = 2,
        b = if (row_idx == 4L) 1 else 0,
        l = 1
      ),
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 9.5),
      legend.key.size = ggplot2::unit(0.55, "cm")
    )
})

FIG3_HISTORICAL <- patchwork::wrap_plots(
  FIG3_PANELS,
  ncol = 3,
  nrow = 4,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = "Basin-average drought-index chronology — SPI and SPEI-PM (1950–2025)",
    subtitle = paste0(
      "Rows 1–2: SPI; rows 3–4: SPEI-PM  |  ",
      "All six accumulation scales: 1, 3, 6, 12, 24, and 36 months  |  ",
      "Dashed line: drought-event threshold (x₀ = −0.5)  |  ",
      "Dotted line: −1.0  |  Grey band: 2022–2025"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 3)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 7.2,
        colour = "grey40",
        hjust = 0.5,
        margin = ggplot2::margin(b = 4)
      ),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 9.5),
      legend.key.size = ggplot2::unit(0.55, "cm"),
      legend.spacing.x = ggplot2::unit(0.4, "cm"),
      legend.margin = ggplot2::margin(t = 2, r = 2, b = 1, l = 2),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
  )

FIG3_PDF <- file.path(
  BASIN_PLOT_DIR,
  "Fig3_MS_SPI_SPEIPM_AllScales_4x3.pdf"
)

FIG3_PNG <- file.path(
  BASIN_PLOT_DIR,
  "Fig3_MS_SPI_SPEIPM_AllScales_4x3.png"
)

tryCatch({
  ggplot2::ggsave(
    FIG3_PDF,
    FIG3_HISTORICAL,
    width = 12,
    height = 12,
    units = "in",
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    FIG3_PNG,
    FIG3_HISTORICAL,
    width = 12,
    height = 12,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  cat(sprintf(
    "  ✓ Manuscript 12-panel SPI + SPEI-PM chronology saved:\n    %s\n    %s\n",
    FIG3_PDF,
    FIG3_PNG
  ))
}, error = function(e) {
  cat(sprintf("  ⚠ Fig. 3 SPI + SPEI-PM chronology failed: %s\n", e$message))
})


################################################################################
# PART 4e-D — SUPPLEMENTARY FIGURE S1
#              1- AND 6-MONTH HISTORICAL DROUGHT CHRONOLOGY
#
# Purpose:
#   Provide the supplementary short-term accumulation-scale chronology
#   corresponding to the manuscript statement that the 1- and 6-month
#   drought-index series show the same broad temporal pattern as the
#   principal 3-, 12-, 24-, and 36-month series.
#
# Layout:
#
#                         1 mo                 6 mo
#              SPI        (a)                 (b)
#              SPEI-PM    (c)                 (d)
#              SPEI-ERA5   (e)                 (f)
#
#   The panel construction, drought-category shading, thresholds,
#   2022–2025 highlighting, colours, and historical period are inherited
#   directly from make_historical_panel_4e(), so S1 remains visually
#   consistent with Fig. 2 and the other historical chronology products.
#
# Output:
#   FigS1_MS_AllIndices_3x2_Scales1_6.pdf
#   FigS1_MS_AllIndices_3x2_Scales1_6.png
################################################################################

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("PART 4e-D: Supplementary Fig. S1 — 1- and 6-month chronology\n")
cat("══════════════════════════════════════════════════════════════════\n")

S1_HIST_SCALES <- c(1L, 6L)
S1_HIST_INDEXES <- c("spi", "spei_pm", "spei_era5")
S1_HIST_INDEX_LABELS <- c("SPI", "SPEI-PM", "SPEI-ERA5")
S1_HIST_SCALE_LABELS <- paste0(S1_HIST_SCALES, " mo")
S1_HIST_PANEL_LETTERS <- letters[seq_len(6L)]

S1_HIST_PANEL_ORDER <- unlist(
  lapply(
    c("SPI", "SPEI_PM", "SPEI_ERA5"),
    function(idx) paste0(idx, S1_HIST_SCALES)
  ),
  use.names = FALSE
)

S1_HISTORICAL_DATA <- all_data[
  names(all_data) %in% S1_HIST_PANEL_ORDER
]

if (!length(S1_HISTORICAL_DATA)) {
  stop(
    "No data available for Supplementary Fig. S1. Expected SPI, SPEI-PM, " ,
    "and SPEI-ERA5 at scales 1 and 6 months."
  )
}

S1_MISSING <- setdiff(S1_HIST_PANEL_ORDER, names(S1_HISTORICAL_DATA))
if (length(S1_MISSING)) {
  stop(
    "Supplementary Fig. S1 is missing required panel(s): ",
    paste(S1_MISSING, collapse = ", ")
  )
}

S1_PANELS <- lapply(seq_along(S1_HIST_PANEL_ORDER), function(i) {
  row_idx <- ceiling(i / length(S1_HIST_SCALES))
  col_idx <- ((i - 1L) %% length(S1_HIST_SCALES)) + 1L
  
  if (row_idx == 1L) {
    index_type <- "spi"
  } else if (row_idx == 2L) {
    index_type <- "spei_pm"
  } else {
    index_type <- "spei_era5"
  }
  
  scale <- S1_HIST_SCALES[col_idx]
  
  p <- make_historical_panel_4e(
    sc = scale,
    panel_lab = sprintf("(%s)", S1_HIST_PANEL_LETTERS[i]),
    index_type = index_type,
    show_xlab = row_idx == 3L,
    show_ylabel = FALSE,
    is_first = i == 1L
  )
  
  p +
    ggplot2::labs(
      title = if (row_idx == 1L) S1_HIST_SCALE_LABELS[col_idx] else NULL,
      y = if (col_idx == 1L) S1_HIST_INDEX_LABELS[row_idx] else NULL
    ) +
    ggplot2::theme(
      plot.title = if (row_idx == 1L) {
        ggplot2::element_text(
          size = 15,
          face = "bold",
          hjust = 0.5,
          margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
        )
      } else {
        ggplot2::element_blank()
      },
      axis.title.y = if (col_idx == 1L) {
        ggplot2::element_text(
          size = 15,
          face = "bold",
          angle = 90,
          margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0)
        )
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = ggplot2::element_text(size = 8),
      axis.ticks.y = ggplot2::element_line(),
      axis.text.x = if (row_idx == 3L) {
        ggplot2::element_text(size = 8)
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (row_idx == 3L) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      },
      axis.title.x = if (row_idx == 3L) {
        ggplot2::element_text(
          size = 9,
          margin = ggplot2::margin(t = 2)
        )
      } else {
        ggplot2::element_blank()
      },
      plot.margin = ggplot2::margin(
        t = 0,
        r = 2,
        b = if (row_idx == 3L) 1 else 0,
        l = 1
      ),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10.5),
      legend.key.size = ggplot2::unit(0.60, "cm")
    )
})

S1_FIGURE <- patchwork::wrap_plots(
  S1_PANELS,
  ncol = 2,
  nrow = 3,
  guides = "collect"
) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10.5),
      legend.key.size = ggplot2::unit(0.60, "cm"),
      legend.spacing.x = ggplot2::unit(0.45, "cm"),
      legend.margin = ggplot2::margin(t = 2, r = 2, b = 1, l = 2),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
  )

S1_PDF <- file.path(
  BASIN_PLOT_DIR,
  "FigS1_MS_AllIndices_3x2_Scales1_6.pdf"
)

S1_PNG <- file.path(
  BASIN_PLOT_DIR,
  "FigS1_MS_AllIndices_3x2_Scales1_6.png"
)

tryCatch({
  ggplot2::ggsave(
    S1_PDF,
    S1_FIGURE,
    width = 12,
    height = 9.5,
    units = "in",
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
  
  ggplot2::ggsave(
    S1_PNG,
    S1_FIGURE,
    width = 12,
    height = 9.5,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
  
  cat(sprintf(
    "  ✓ Supplementary Fig. S1 saved:\n    %s\n    %s\n",
    S1_PDF,
    S1_PNG
  ))
}, error = function(e) {
  cat(sprintf("  ⚠ Supplementary Fig. S1 failed: %s\n", e$message))
})


# PART 5 – EXCEL SUMMARY (statistics + drought event counts + correlations)
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 5: Excel summary\n")
cat("══════════════════════════════════════════════\n")

## Build one stats row for a given index × scale.
make_stats_row <- function(idx, sc) {
  df <- tryCatch(load_any_ts(idx, sc), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  # Events are detected on the FULL record (so events that began before 2020
  # keep their true duration) and counted when they overlap 2020-2025.
  ev_count <- tryCatch({
    ev <- identify_drought_events_ts(df$date, df$value, scale = sc)
    sum(ev$end_date >= as.Date("2020-01-01") & ev$start_date <= as.Date("2025-12-31"))
  }, error = function(e) NA_integer_)
  data.frame(
    Timescale                = sprintf("%02d", sc),
    Index                    = if (idx == "spei") "SPEI-PM" else if (idx == "spei_era5") "SPEI-ERA5" else toupper(idx),
    Mean                     = round(mean(df$value,   na.rm = TRUE), 4),
    Median                   = round(median(df$value, na.rm = TRUE), 4),
    StdDev                   = round(sd(df$value,     na.rm = TRUE), 4),
    Min                      = round(min(df$value,    na.rm = TRUE), 4),
    Max                      = round(max(df$value,    na.rm = TRUE), 4),
    Drought_Events_2020_2025 = ev_count,
    stringsAsFactors = FALSE)
}

export_basin_excel <- function(output_file) {
  wb        <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  
  ## ── Sheet 1: Summary Statistics ─────────────────────────────────────────────
  stats_rows <- list()
  for (idx in c("spi", "spei", "spei_era5", "swei", "sspi")) {
    scales <- if (idx == "spi") SPI_SCALES else if (idx == "spei_era5") SPEI_ERA5_SCALES else if (idx == "swei") SWEI_SCALE else if (idx == "sspi") SSPI_SCALE else SPEI_SCALES
    for (sc in scales) {
      r <- tryCatch(make_stats_row(idx, sc),
                    error = function(e) {
                      cat(sprintf("  ⚠ Stats %s-%02d: %s\n",
                                  toupper(idx), sc, e$message))
                      NULL
                    })
      if (!is.null(r)) stats_rows[[length(stats_rows) + 1L]] <- r
    }
  }
  for (idx in c("mspi", "mspei")) {
    r <- tryCatch(make_stats_row(idx, MSPI_MSPEI_SCALE),
                  error = function(e) NULL)
    if (!is.null(r)) stats_rows[[length(stats_rows) + 1L]] <- r
  }
  
  stats_df <- if (length(stats_rows)) do.call(rbind, stats_rows) else data.frame(Index = character(), Scale = integer())
  openxlsx::addWorksheet(wb, "Summary_Statistics")
  openxlsx::writeData(wb, "Summary_Statistics", stats_df)
  openxlsx::addStyle(wb, "Summary_Statistics", hdr_style,
                     rows = 1, cols = seq_len(ncol(stats_df)))
  cat(sprintf("  ✓ Summary_Statistics sheet: %d rows\n", nrow(stats_df)))
  
  ## ── Sheet 2: Pairwise Correlations ──────────────────────────────────────────
  ## Automatically correlate every available index pair on their common date range.
  corr_df <- tryCatch({
    ## -------------------------------------------------------------------------
    ## Build correlation dataset from ALL manuscript drought indices and scales
    ## -------------------------------------------------------------------------
    series_list <- list()
    
    # SPI
    for (sc in SPI_SCALES) {
      ts <- tryCatch(load_any_ts("spi", sc), error = function(e) NULL)
      if (!is.null(ts))
        series_list[[sprintf("SPI%d", sc)]] <- ts
    }
    
    # SPEI (Penman–Monteith)
    for (sc in SPEI_SCALES) {
      ts <- tryCatch(load_any_ts("spei", sc), error = function(e) NULL)
      if (!is.null(ts))
        series_list[[sprintf("SPEI_PM%d", sc)]] <- ts
    }
    
    # SPEI (ERA5-Land PEV)
    for (sc in SPEI_ERA5_SCALES) {
      ts <- tryCatch(load_any_ts("spei_era5", sc), error = function(e) NULL)
      if (!is.null(ts))
        series_list[[sprintf("SPEI_ERA5%d", sc)]] <- ts
    }
    
    # Snow drought indices: included where authoritative basin averages exist.
    ts <- tryCatch(load_any_ts("swei", SWEI_SCALE), error = function(e) NULL)
    if (!is.null(ts)) series_list[[sprintf("SWEI%d", SWEI_SCALE)]] <- ts
    ts <- tryCatch(load_any_ts("sspi", SSPI_SCALE), error = function(e) NULL)
    if (!is.null(ts)) series_list[[sprintf("SSPI%d", SSPI_SCALE)]] <- ts
    
    # Optional multivariate indices
    ts <- tryCatch(load_any_ts("mspi", MSPI_MSPEI_SCALE), error = function(e) NULL)
    if (!is.null(ts))
      series_list[[sprintf("MSPI%d", MSPI_MSPEI_SCALE)]] <- ts
    
    ts <- tryCatch(load_any_ts("mspei", MSPI_MSPEI_SCALE), error = function(e) NULL)
    if (!is.null(ts))
      series_list[[sprintf("MSPEI%d", MSPI_MSPEI_SCALE)]] <- ts
    
    series_list <- Filter(Negate(is.null), series_list)
    nms <- names(series_list)
    
    if (!length(nms)) {
      data.frame(Comparison = character(), Correlation = numeric(), N_months = integer(), stringsAsFactors = FALSE)
    } else {
      # Pairwise date joins avoid truncating all correlations to the shortest series.
      pairs <- if (length(nms) >= 2L) combn(nms, 2L, simplify = FALSE) else list()
      if (!length(pairs)) {
        data.frame(Comparison = character(), Correlation = numeric(), N_months = integer(), stringsAsFactors = FALSE)
      } else do.call(rbind, lapply(pairs, function(pair) {
        a <- series_list[[pair[1L]]][, c("date", "value"), drop = FALSE]
        b <- series_list[[pair[2L]]][, c("date", "value"), drop = FALSE]
        names(a)[2L] <- "x"; names(b)[2L] <- "y"
        joined <- merge(a, b, by = "date", all = FALSE)
        ok <- is.finite(joined$x) & is.finite(joined$y)
        n_ok <- sum(ok)
        r <- if (n_ok >= 3L && stats::sd(joined$x[ok]) > 0 && stats::sd(joined$y[ok]) > 0) stats::cor(joined$x[ok], joined$y[ok]) else NA_real_
        data.frame(Comparison = sprintf("%s vs %s", pair[1L], pair[2L]), Correlation = round(r, 4), N_months = n_ok, stringsAsFactors = FALSE)
      }))
    }
  }, error = function(e) {
    cat(" ⚠ Correlation sheet skipped:", e$message, "\n")
    data.frame(Comparison = character(), Correlation = numeric(), 
               N_months = integer(), stringsAsFactors = FALSE)
  })
  openxlsx::addWorksheet(wb, "Correlations")
  openxlsx::writeData(wb, "Correlations", corr_df)
  openxlsx::addStyle(wb, "Correlations", hdr_style,
                     rows = 1, cols = seq_len(ncol(corr_df)))
  cat(sprintf("  ✓ Correlations sheet: %d pairs\n", nrow(corr_df)))
  
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  cat(sprintf("✓ Excel summary saved: %s\n", basename(output_file)))
  invisible(wb)
}

excel_out <- file.path(BASIN_PLOT_DIR, "Drought_Summary.xlsx")
export_basin_excel(excel_out)

################################################################################
# FINAL SUMMARY
################################################################################
cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  12basin_timeseries.R  DONE                 ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

cat("Outputs:\n")
cat("  Per-index flat CSVs   : ", normalizePath(BASIN_TS_DIR),   "\n")
cat("  Combined long CSV     : ", normalizePath(BASIN_PLOT_DIR),  "\n")
cat("  Time series PNGs      : ", normalizePath(BASIN_PLOT_DIR),  "\n")
cat("  MS figure SPEI (PM)   : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4_MS_SPEI_123_PM_timeseries.pdf")), "\n")
cat("  MS figure SPEI (ERA5)  : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4_MS_SPEI_123_ERA5_timeseries.pdf")), "\n")
cat("  MS figure SPEI (both) : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4_MS_SPEI_123_overlay_timeseries.pdf")), "\n")
cat("  MS figure historical 18-panel : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig2_MS_AllIndices_3x6.pdf")), "\n")
cat("  Recent 6-panel comparison : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig2S_MS_Recent_AllScales_3x2.pdf")), "\n")
cat("  Spatial representativeness (PDF): ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4c_Pearson_r_SPEI_PM_SPI.pdf")), "\n")
cat("  Spatial representativeness (PNG): ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4c_Pearson_r_SPEI_PM_SPI.png")), "\n")
cat("  SPEI-PM only fig (PDF)          : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4c_Pearson_r_SPEI136_PM.pdf")), "\n")
cat("  SPEI-PM only fig (PNG)          : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4c_Pearson_r_SPEI136_PM.png")), "\n")
cat("  Excel summary         : ", normalizePath(excel_out),       "\n")
cat("  Seasonality P/PET(PM) : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "FigS_Seasonality_P_PET.pdf")), "\n")
cat("\n✓ All basin-averaged output tasks complete!\n")
cat("  Next: run 13trends_visualization.R for spatial figures.\n\n")

################################################################################
# REGIONAL EXTENSION — 7 SUBBASINS + WHOLE NECHAKO BASIN
#
# The original basin outputs remain untouched. This section reads the
# authoritative regional CSVs produced by 11trend_test_ALL.R and creates
# region-aware time-series tables and figures without recomputing spatial means.
################################################################################
REGIONAL_TS_DIR <- file.path(TREND_DIR, "regional_timeseries")
REGIONAL_PLOT_DIR <- file.path(TREND_DIR, "regional_timeseries_plots")
dir.create(REGIONAL_TS_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(REGIONAL_PLOT_DIR, recursive=TRUE, showWarnings=FALSE)

regional_index_types <- c("spi", "spei", "spei_era5", "swei", "sspi")
regional_scales <- MANUSCRIPT_SCALES
regional_files <- list.files(
  REGIONAL_TS_DIR,
  pattern = "_regional_averaged_by_month\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

regional_all <- data.table::rbindlist(lapply(regional_files, function(f) {
  x <- data.table::fread(f)
  if (!all(c("region_id","region_name","date","value") %in% names(x))) return(NULL)
  x[, date := as.Date(date)]
  x[, value := as.numeric(value)]
  x[, source_file := basename(f)]
  x
}), use.names=TRUE, fill=TRUE)

if (nrow(regional_all)) {
  regional_all <- regional_all[
    tolower(index_type) %in% regional_index_types & as.integer(scale) %in% regional_scales
  ]
  setorder(regional_all, region_id, index_type, scale, date)
  fwrite(regional_all, file.path(REGIONAL_TS_DIR, "regional_mean_all_indices.csv"))
  
  ################################################################################
  # REGIONAL RECENT FIGURE — 7 SUBBASINS × 6 SCALES
  #
  # Purpose:
  #   Produce a direct regional analogue of Fig2S_MS_Recent_AllScales_3x2.
  #   Rows are the seven supplied subbasins; columns are the six manuscript
  #   accumulation scales. Every panel overlays SPI and SPEI-PM for 2022–2025.
  #
  # Scientific/data contract:
  #   * Reads the authoritative regional CSVs produced by 11trend_test_ALL.R.
  #   * Does not recompute regional spatial means.
  #   * Excludes the whole-basin dissolved union.
  #   * Excludes SPEI-ERA5 from the figure by design.
  ################################################################################
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("REGIONAL RECENT FIGURE: 7 subbasins × 6 accumulation scales\n")
  cat("                 SPI + SPEI-PM, 2022–2025\n")
  cat("══════════════════════════════════════════════════════════════════\n")
  
  tryCatch({
    if (!exists("regional_all", inherits = TRUE) || !nrow(regional_all)) {
      cat("  ⚠ Regional CSVs are not available; subbasin figure skipped.\n")
    } else {
      regional_recent <- regional_all[
        tolower(index_type) %in% c("spi", "spei") &
          as.integer(scale) %in% MANUSCRIPT_TS_SCALES &
          date >= RECENT_START_4e & date <= RECENT_END_4e
      ]
      
      # The authoritative regional products contain the seven supplied subbasins
      # plus the dissolved whole-basin union. Remove only the whole-basin region.
      whole_id <- NULL
      if (exists("NECHAKO_REGION_BASIN_ID", inherits = TRUE)) {
        whole_id <- get("NECHAKO_REGION_BASIN_ID", inherits = TRUE)
      }
      if (is.null(whole_id) || !length(whole_id)) {
        whole_id_candidates <- unique(
          regional_recent$region_id[
            grepl("whole|basin|union", regional_recent$region_name, ignore.case = TRUE)
          ]
        )
        if (length(whole_id_candidates) == 1L) whole_id <- whole_id_candidates
      }
      
      if (!is.null(whole_id) && length(whole_id)) {
        regional_recent <- regional_recent[region_id != whole_id[1L]]
      }
      
      subbasin_ids <- unique(regional_recent$region_id)
      if (length(subbasin_ids) != 7L) {
        stop(
          "Expected exactly 7 supplied subbasins after excluding the whole-basin union; found ",
          length(subbasin_ids), ". Check regional_timeseries/region_catalog.csv and ",
          "the authoritative load_nechako_regions() definition."
        )
      }
      
      # Stable row order: use the region_catalog written by 11trend_test_ALL.R when
      # available. Otherwise preserve the order of the regional product files.
      region_catalog_file <- file.path(REGIONAL_TS_DIR, "region_catalog.csv")
      if (file.exists(region_catalog_file)) {
        rc <- data.table::fread(region_catalog_file)
        if (all(c("region_id", "region_name") %in% names(rc))) {
          rc_ids <- unique(as.character(rc$region_id))
          subbasin_ids <- c(
            rc_ids[rc_ids %in% subbasin_ids],
            setdiff(subbasin_ids, rc_ids)
          )
        }
      }
      
      regional_recent[, region_id := factor(region_id, levels = subbasin_ids)]
      region_labels <- unique(
        regional_recent[, .(region_id, region_name)]
      )
      region_labels <- region_labels[match(subbasin_ids, as.character(region_id))]
      regional_recent[, region_name_plot := factor(
        region_name,
        levels = region_labels$region_name
      )]
      regional_recent[, scale_plot := factor(
        as.integer(scale),
        levels = MANUSCRIPT_TS_SCALES,
        labels = paste0(MANUSCRIPT_TS_SCALES, "-month")
      )]
      
      regional_recent$Index <- factor(
        regional_recent$index_type,
        levels = c("spi", "spei"),
        labels = c("SPI", "SPEI-PM")
      )
      
      fig_subbasins <- ggplot2::ggplot(
        regional_recent,
        ggplot2::aes(x = date, y = value, colour = Index, group = Index)
      ) +
        ggplot2::geom_hline(
          yintercept = 0,
          colour = "grey40",
          linewidth = 0.25
        ) +
        ggplot2::geom_hline(
          yintercept = DROUGHT_ONSET,
          colour = "grey25",
          linewidth = 0.35,
          linetype = "dashed"
        ) +
        ggplot2::geom_hline(
          yintercept = -1,
          colour = "grey50",
          linewidth = 0.22,
          linetype = "dotted"
        ) +
        ggplot2::geom_line(
          linewidth = 0.48,
          na.rm = TRUE
        ) +
        ggplot2::scale_colour_manual(
          values = LINE_COLOURS_4e[c("SPI", "SPEI-PM")],
          name = "Index",
          guide = ggplot2::guide_legend(
            title.position = "top",
            nrow = 1,
            override.aes = list(linewidth = 1.0)
          )
        ) +
        ggplot2::scale_x_date(
          date_breaks = "2 years",
          date_labels = "%Y",
          expand = ggplot2::expansion(add = c(0, 0))
        ) +
        ggplot2::scale_y_continuous(
          limits = c(-3.2, 2.5),
          breaks = c(-3, -2, -1, 0, 1, 2),
          expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(region_name_plot),
          cols = ggplot2::vars(scale_plot),
          scales = "fixed"
        ) +
        ggplot2::labs(
          title = "2022–2025 drought evolution across Nechako subbasins",
          subtitle = paste0(
            "Each panel overlays SPI and SPEI-PM; rows = seven supplied subbasins; ",
            "columns = accumulation scale. Dashed line: x₀ = −0.5; dotted line: −1.0."
          ),
          x = "Year",
          y = "Standardized index",
          colour = "Index"
        ) +
        ggplot2::theme_classic(base_size = 8) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 11,
            face = "bold",
            hjust = 0.5,
            margin = ggplot2::margin(b = 2)
          ),
          plot.subtitle = ggplot2::element_text(
            size = 7,
            colour = "grey40",
            hjust = 0.5,
            margin = ggplot2::margin(b = 4)
          ),
          strip.background = ggplot2::element_rect(
            fill = "grey95",
            colour = "grey70",
            linewidth = 0.25
          ),
          strip.text.x = ggplot2::element_text(size = 7, face = "bold"),
          strip.text.y = ggplot2::element_text(size = 7, face = "bold"),
          axis.text.x = ggplot2::element_text(size = 5.5),
          axis.text.y = ggplot2::element_text(size = 5.5),
          axis.title.x = ggplot2::element_text(size = 7),
          axis.title.y = ggplot2::element_text(size = 7),
          axis.ticks = ggplot2::element_line(linewidth = 0.25),
          legend.position = "bottom",
          legend.text = ggplot2::element_text(size = 7),
          legend.title = ggplot2::element_text(size = 7),
          legend.key.size = ggplot2::unit(0.35, "cm"),
          # Thin panel borders provide explicit vertical separators between
          # adjacent accumulation-scale panels, matching the visual separation
          # visible in Fig2S_MS_Recent_AllScales_3x2.
          panel.border = ggplot2::element_rect(
            colour = "grey65",
            fill = NA,
            linewidth = 0.35
          ),
          panel.spacing = ggplot2::unit(0.18, "lines"),
          panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(3, 4, 2, 3)
        )
      
      for (ext in c("pdf", "png")) {
        out <- file.path(
          REGIONAL_PLOT_DIR,
          paste0("Fig2S_MS_Recent_Subbasins_7x6.", ext)
        )
        tryCatch(
          ggplot2::ggsave(
            out,
            fig_subbasins,
            width = 16,
            height = 14,
            units = "in",
            dpi = if (ext == "png") 300 else "print",
            device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png",
            limitsize = FALSE
          ),
          error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message))
        )
        cat(sprintf("  ✓ Saved: %s\n", basename(out)))
      }
      
      cat(sprintf(
        "  ✓ Seven-subbasin recent figure complete: %d subbasins × %d scales × 2 indices\n",
        length(subbasin_ids),
        length(MANUSCRIPT_TS_SCALES)
      ))
      
      ################################################################################
      # REGIONAL RECENT FIGURE — 7 SUBBASINS × 3 SELECTED SCALES
      #
      # Purpose:
      #   Produce a compact regional recent-period figure restricted to the
      #   manuscript's 3-, 12-, and 36-month accumulation scales. Each panel
      #   overlays SPI and SPEI-PM for 2022–2025.
      #
      # Layout:
      #   7 rows = the seven supplied subbasins
      #   3 columns = 3-, 12-, and 36-month accumulation scales
      #
      # Output:
      #   Fig2S_MS_Recent_Subbasins_7x3_Scales3_12_36.pdf/.png
      ################################################################################
      
      regional_recent_7x3 <- regional_recent[
        as.integer(scale) %in% c(3L, 12L, 36L)
      ]
      
      if (!nrow(regional_recent_7x3)) {
        stop(
          "No regional SPI/SPEI-PM data are available for the requested " ,
          "3-, 12-, and 36-month scales."
        )
      }
      
      regional_recent_7x3[, scale_plot := factor(
        as.integer(scale),
        levels = c(3L, 12L, 36L),
        labels = c("3-month", "12-month", "36-month")
      )]
      
      fig_subbasins_7x3 <- ggplot2::ggplot(
        regional_recent_7x3,
        ggplot2::aes(x = date, y = value, colour = Index, group = Index)
      ) +
        ggplot2::geom_hline(
          yintercept = 0,
          colour = "grey40",
          linewidth = 0.25
        ) +
        ggplot2::geom_hline(
          yintercept = DROUGHT_ONSET,
          colour = "grey25",
          linewidth = 0.35,
          linetype = "dashed"
        ) +
        ggplot2::geom_hline(
          yintercept = -1,
          colour = "grey50",
          linewidth = 0.22,
          linetype = "dotted"
        ) +
        ggplot2::geom_line(
          linewidth = 0.48,
          na.rm = TRUE
        ) +
        ggplot2::scale_colour_manual(
          values = LINE_COLOURS_4e[c("SPI", "SPEI-PM")],
          name = "Index",
          guide = ggplot2::guide_legend(
            title.position = "top",
            nrow = 1,
            override.aes = list(linewidth = 1.0)
          )
        ) +
        ggplot2::scale_x_date(
          date_breaks = "2 years",
          date_labels = "%Y",
          expand = ggplot2::expansion(add = c(0, 0))
        ) +
        ggplot2::scale_y_continuous(
          limits = c(-3.2, 2.5),
          breaks = c(-3, -2, -1, 0, 1, 2),
          expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(region_name_plot),
          cols = ggplot2::vars(scale_plot),
          scales = "fixed"
        ) +
        ggplot2::labs(
          title = "2022–2025 drought evolution across Nechako subbasins",
          subtitle = paste0(
            "SPI and SPEI-PM at 3-, 12-, and 36-month accumulation scales; ",
            "dashed line: x₀ = −0.5; dotted line: −1.0."
          ),
          x = "Year",
          y = "Standardized index",
          colour = "Index"
        ) +
        ggplot2::theme_classic(base_size = 8) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 11,
            face = "bold",
            hjust = 0.5,
            margin = ggplot2::margin(b = 2)
          ),
          plot.subtitle = ggplot2::element_text(
            size = 7,
            colour = "grey40",
            hjust = 0.5,
            margin = ggplot2::margin(b = 4)
          ),
          strip.background = ggplot2::element_rect(
            fill = "grey95",
            colour = "grey70",
            linewidth = 0.25
          ),
          strip.text.x = ggplot2::element_text(size = 7, face = "bold"),
          strip.text.y = ggplot2::element_text(size = 7, face = "bold"),
          axis.text.x = ggplot2::element_text(size = 5.5),
          axis.text.y = ggplot2::element_text(size = 5.5),
          axis.title.x = ggplot2::element_text(size = 7),
          axis.title.y = ggplot2::element_text(size = 7),
          axis.ticks = ggplot2::element_line(linewidth = 0.25),
          legend.position = "bottom",
          legend.text = ggplot2::element_text(size = 7),
          legend.title = ggplot2::element_text(size = 7),
          legend.key.size = ggplot2::unit(0.35, "cm"),
          # Thin panel borders provide explicit vertical separators between
          # adjacent accumulation-scale panels, matching the visual separation
          # visible in Fig2S_MS_Recent_AllScales_3x2.
          panel.border = ggplot2::element_rect(
            colour = "grey65",
            fill = NA,
            linewidth = 0.35
          ),
          panel.spacing = ggplot2::unit(0.18, "lines"),
          panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(3, 4, 2, 3)
        )
      
      for (ext in c("pdf", "png")) {
        out <- file.path(
          REGIONAL_PLOT_DIR,
          paste0("Fig2S_MS_Recent_Subbasins_7x3_Scales3_12_36.", ext)
        )
        tryCatch(
          ggplot2::ggsave(
            out,
            fig_subbasins_7x3,
            width = 12,
            height = 14,
            units = "in",
            dpi = if (ext == "png") 300 else "print",
            device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png",
            limitsize = FALSE
          ),
          error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message))
        )
        cat(sprintf("  ✓ Saved: %s\n", basename(out)))
      }
      
      cat("  ✓ Seven-subbasin 7×3 recent figure complete: scales 3, 12, 36 months × 2 indices\n")
    }
  }, error = function(e) {
    cat(sprintf("  ⚠ Seven-subbasin recent figure failed: %s\n", e$message))
  })
  
  
  
  ################################################################################
  # REGIONAL RECENT FIGURE — SNOW DROUGHT: 7 SUBBASINS × 2 INDEX/SCALE PANELS
  #
  # SWEI is currently available at scale 3 and SSPI at scale 1. They cannot
  # populate the same six-scale grid as SPI/SPEI, so each gets its own column.
  # Regional means are read from Script 11 outputs; this block does not recompute
  # any spatial averages.
  ################################################################################
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("REGIONAL RECENT FIGURE: snow drought across 7 subbasins\n")
  cat("                         SWEI-3 + SSPI-1, 2022–2025\n")
  cat("══════════════════════════════════════════════════════════════════\n")
  
  tryCatch({
    snow_recent <- regional_all[
      tolower(index_type) %in% c("swei", "sspi") &
        ((tolower(index_type) == "swei" & as.integer(scale) == SWEI_SCALE) |
           (tolower(index_type) == "sspi" & as.integer(scale) == SSPI_SCALE)) &
        date >= RECENT_START_4e & date <= RECENT_END_4e
    ]
    
    if (!nrow(snow_recent)) {
      cat("  ⚠ No regional SWEI-3 / SSPI-1 data found; snow subbasin figure skipped.\n")
    } else {
      # Exclude the dissolved whole-basin union, retaining the seven subbasins.
      whole_id_snow <- NULL
      if (exists("NECHAKO_REGION_BASIN_ID", inherits = TRUE)) {
        whole_id_snow <- get("NECHAKO_REGION_BASIN_ID", inherits = TRUE)
      }
      if (is.null(whole_id_snow) || !length(whole_id_snow)) {
        whole_candidates_snow <- unique(
          snow_recent$region_id[
            grepl("whole|basin|union", snow_recent$region_name, ignore.case = TRUE)
          ]
        )
        if (length(whole_candidates_snow) == 1L) whole_id_snow <- whole_candidates_snow
      }
      if (!is.null(whole_id_snow) && length(whole_id_snow)) {
        snow_recent <- snow_recent[region_id != whole_id_snow[1L]]
      } else {
        # Name-based fallback if the catalog does not expose a whole-basin ID.
        snow_recent <- snow_recent[
          !grepl("whole basin|basin union|entire basin", region_name, ignore.case = TRUE)
        ]
      }
      
      # Require seven unique subbasins; fail with an explicit diagnostic rather
      # than silently presenting a partial regional comparison.
      snow_region_ids <- unique(as.character(snow_recent$region_id))
      if (length(snow_region_ids) != 7L) {
        stop("Expected 7 subbasins for snow indices after excluding the whole-basin union; found ",
             length(snow_region_ids), ". Check Script 11 regional SWEI/SSPI CSVs.")
      }
      
      # Follow the authoritative region catalog order when available.
      snow_catalog <- file.path(REGIONAL_TS_DIR, "region_catalog.csv")
      if (file.exists(snow_catalog)) {
        snow_rc <- data.table::fread(snow_catalog)
        if (all(c("region_id", "region_name") %in% names(snow_rc))) {
          snow_rc_ids <- unique(as.character(snow_rc$region_id))
          snow_region_ids <- c(snow_rc_ids[snow_rc_ids %in% snow_region_ids],
                               setdiff(snow_region_ids, snow_rc_ids))
        }
      }
      snow_recent[, region_id := factor(as.character(region_id), levels = snow_region_ids)]
      snow_region_labels <- unique(snow_recent[, .(region_id, region_name)])
      snow_region_labels <- snow_region_labels[match(snow_region_ids, as.character(region_id))]
      snow_recent[, region_name_plot := factor(region_name,
                                               levels = snow_region_labels$region_name)]
      snow_recent[, IndexPanel := ifelse(tolower(index_type) == "swei", "SWEI-3", "SSPI-1")]
      snow_recent[, IndexPanel := factor(IndexPanel, levels = c("SWEI-3", "SSPI-1"))]
      
      snow_colours <- c("SWEI-3" = "#1E90FF", "SSPI-1" = "#9467bd")
      fig_snow_subbasins <- ggplot2::ggplot(
        snow_recent,
        ggplot2::aes(x = date, y = value, colour = IndexPanel, group = IndexPanel)
      ) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.25) +
        ggplot2::geom_hline(yintercept = DROUGHT_ONSET, colour = "grey25",
                            linewidth = 0.35, linetype = "dashed") +
        ggplot2::geom_hline(yintercept = -1, colour = "grey50", linewidth = 0.22,
                            linetype = "dotted") +
        ggplot2::geom_line(linewidth = 0.48, na.rm = TRUE) +
        ggplot2::scale_colour_manual(values = snow_colours, name = "Index") +
        ggplot2::scale_x_date(date_breaks = "2 years", date_labels = "%Y",
                              expand = ggplot2::expansion(add = c(0, 0))) +
        ggplot2::scale_y_continuous(limits = c(-3.2, 2.5),
                                    breaks = c(-3, -2, -1, 0, 1, 2),
                                    expand = ggplot2::expansion(mult = c(0, 0))) +
        ggplot2::facet_grid(rows = ggplot2::vars(region_name_plot),
                            cols = ggplot2::vars(IndexPanel), scales = "fixed",
                            drop = FALSE) +
        ggplot2::labs(
          title = "2022–2025 snow drought evolution across Nechako subbasins",
          subtitle = paste0("SWEI-3 (3-month) and SSPI-1 (1-month); rows = seven supplied subbasins. ",
                            "Dashed line: x₀ = −0.5; dotted line: −1.0."),
          x = "Year", y = "Standardized index", colour = "Index"
        ) +
        ggplot2::theme_classic(base_size = 8) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5,
                                             margin = ggplot2::margin(b = 2)),
          plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40", hjust = 0.5,
                                                margin = ggplot2::margin(b = 4)),
          strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70",
                                                   linewidth = 0.25),
          strip.text.x = ggplot2::element_text(size = 8, face = "bold"),
          strip.text.y = ggplot2::element_text(size = 7, face = "bold"),
          axis.text.x = ggplot2::element_text(size = 5.5),
          axis.text.y = ggplot2::element_text(size = 5.5),
          axis.title.x = ggplot2::element_text(size = 7),
          axis.title.y = ggplot2::element_text(size = 7),
          axis.ticks = ggplot2::element_line(linewidth = 0.25),
          legend.position = "bottom", legend.text = ggplot2::element_text(size = 7),
          legend.title = ggplot2::element_text(size = 7),
          legend.key.size = ggplot2::unit(0.35, "cm"),
          panel.border = ggplot2::element_rect(colour = "grey65", fill = NA,
                                               linewidth = 0.35),
          panel.spacing = ggplot2::unit(0.18, "lines"),
          panel.grid = ggplot2::element_blank(),
          plot.margin = ggplot2::margin(3, 4, 2, 3)
        )
      
      for (ext in c("pdf", "png")) {
        snow_out <- file.path(REGIONAL_PLOT_DIR,
                              paste0("Fig2S_Recent_Subbasins_SnowDrought_7x2.", ext))
        tryCatch(
          ggplot2::ggsave(
            snow_out, fig_snow_subbasins, width = 10, height = 14, units = "in",
            dpi = if (ext == "png") 300 else "print",
            device = if (ext == "pdf") {
              if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
            } else "png",
            limitsize = FALSE
          ),
          error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message))
        )
        cat(sprintf("  ✓ Saved: %s\n", basename(snow_out)))
      }
      cat(sprintf("  ✓ Snow subbasin figure complete: %d subbasins × 2 index/scale panels\n",
                  length(snow_region_ids)))
    }
  }, error = function(e) {
    cat(sprintf("  ⚠ Snow-drought subbasin figure failed: %s\n", e$message))
  })
  
  
  # One diagnostic/publication plot per region × index × scale.
  make_regional_plot <- function(d, idx, sc, region_name) {
    ggplot2::ggplot(d, ggplot2::aes(date, value)) +
      ggplot2::geom_hline(yintercept = DROUGHT_ONSET, linetype = "dashed", linewidth = 0.35) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.25) +
      ggplot2::geom_line(linewidth = 0.45, na.rm=TRUE) +
      ggplot2::labs(
        title = sprintf("%s — %s-%d", region_name, toupper(idx), as.integer(sc)),
        subtitle = "Area-weighted mean of pixel-level drought-index values; fractional boundary-cell coverage",
        x = NULL, y = "Index value"
      ) +
      ggplot2::theme_bw(base_size=11)
  }
  
  for (rid in unique(regional_all$region_id)) {
    rname <- unique(regional_all$region_name[regional_all$region_id == rid])[1]
    for (idx in regional_index_types) {
      for (sc in regional_scales) {
        d <- regional_all[region_id == rid & tolower(index_type) == idx & scale == sc]
        if (!nrow(d)) next
        ggplot2::ggsave(
          file.path(REGIONAL_PLOT_DIR, sprintf("%s_%s_%02d_timeseries.png", rid, idx, sc)),
          make_regional_plot(d, idx, sc, rname), width=12, height=5.5, dpi=300
        )
      }
    }
  }
  
  # Comparative figure for each scale: regional facets × available indices (SPI/SPEI/SWEI/SSPI).
  for (sc in intersect(regional_scales, unique(regional_all$scale))) {
    d <- regional_all[scale == sc & index_type %in% regional_index_types]
    if (!nrow(d)) next
    d[, Index := factor(
      ifelse(index_type == "spei", "SPEI-PM",
             ifelse(index_type == "spei_era5", "SPEI-ERA5",
                    ifelse(index_type == "swei", "SWEI",
                           ifelse(index_type == "sspi", "SSPI", "SPI")))),
      levels=c("SPI","SPEI-PM","SPEI-ERA5","SWEI","SSPI")
    )]
    p <- ggplot2::ggplot(d, ggplot2::aes(date, value)) +
      ggplot2::geom_hline(yintercept=DROUGHT_ONSET, linetype="dashed", linewidth=0.25) +
      ggplot2::geom_line(linewidth=0.25, na.rm=TRUE) +
      ggplot2::facet_grid(Index ~ region_name, scales="free_y") +
      ggplot2::labs(title=sprintf("Nechako regional drought indices — %d-month scale", sc),
                    x=NULL, y="Index value") +
      ggplot2::theme_bw(base_size=8) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1),
                     strip.text=ggplot2::element_text(size=7))
    ggplot2::ggsave(file.path(REGIONAL_PLOT_DIR, sprintf("regional_all_regions_scale_%02d.png", sc)),
                    p, width=18, height=8, dpi=300)
  }
  
  # Regional descriptive summary.
  regional_summary <- regional_all[, .(
    n_months = sum(is.finite(value)),
    mean_value = mean(value, na.rm=TRUE),
    sd_value = sd(value, na.rm=TRUE),
    min_value = min(value, na.rm=TRUE),
    max_value = max(value, na.rm=TRUE),
    drought_months = sum(value < DROUGHT_ONSET, na.rm=TRUE),
    drought_fraction_pct = 100 * mean(value < DROUGHT_ONSET, na.rm=TRUE)
  ), by=.(region_id, region_name, index_type, scale)]
  fwrite(regional_summary, file.path(REGIONAL_TS_DIR, "regional_index_summary.csv"))
  cat(sprintf("✓ Regional time-series products written for %d regions.\n",
              uniqueN(regional_all$region_id)))
} else {
  cat("⚠ No regional index files found; run 11trend_test_ALL.R first.\n")
}

################################################################################
# PART 6 — OTHER HYDROLOGICAL INDICES WITH DIFFERENT SPATIAL SUPPORT
#
# SSI is station-specific, TSDI is derived from basin-wide lake storage, and
# NTSDI represents the Nechako Reservoir. These are deliberately exported to a
# separate table instead of being labelled as area-weighted basin means.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 6: SSI, TSDI and NTSDI support-aware time-series exports\n")
cat("══════════════════════════════════════════════════════════════\n")

OTHER_TS_DIR <- file.path(WD_PATH, "temporal_drought", "other_index_timeseries")
OTHER_TS_PLOT_DIR <- file.path(OTHER_TS_DIR, "plots")
dir.create(OTHER_TS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OTHER_TS_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
other_ts_rows <- list()

add_other_series <- function(date, value, index, support_type, site_id, site_name,
                             source_label) {
  date <- as.Date(date)
  value <- as.numeric(value)
  if (length(date) != length(value)) return(invisible(NULL))
  d <- data.frame(
    date = date, value = value, index = index,
    support_type = support_type, site_id = site_id,
    site_name = site_name, source = source_label,
    stringsAsFactors = FALSE
  )
  d <- d[!is.na(d$date), , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  other_ts_rows[[length(other_ts_rows) + 1L]] <<- d
  
  safe_id <- gsub("[^A-Za-z0-9_-]+", "_", paste(index, site_id, source_label, sep = "_"))
  utils::write.csv(d, file.path(OTHER_TS_DIR, paste0(safe_id, ".csv")), row.names = FALSE)
  plot_df <- d[is.finite(d$value), , drop = FALSE]
  if (nrow(plot_df)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(date, value)) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
      ggplot2::geom_hline(yintercept = DROUGHT_ONSET, linetype = "dashed", linewidth = 0.4) +
      ggplot2::geom_line(linewidth = 0.5, colour = "#1F77B4", na.rm = TRUE) +
      ggplot2::labs(
        title = sprintf("%s — %s", index, site_name),
        subtitle = sprintf("Support: %s; source: %s", support_type, source_label),
        x = "Year", y = "Standardized index value"
      ) + ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    ggplot2::ggsave(file.path(OTHER_TS_PLOT_DIR, paste0(safe_id, ".png")),
                    p, width = 11, height = 5, dpi = 200, bg = "white")
  }
  invisible(d)
}

# Surface-water storage indices: use normalized (standardized) output columns.
ntsdi_dir <- file.path(WD_PATH, "ntsdi_results")
ntsdi_specs <- list(
  list(file = "Nechako_TSDI_Pipeline1_Output.csv", col = "TSDI_norm",
       index = "TSDI", support = "basin-wide lake-storage aggregate",
       site = "Nechako Basin lakes", source = "7_WBInex_NTSDI_GeoLakes_LocalRS.r"),
  list(file = "Nechako_NTSDI_Pipeline2_Output.csv", col = "NTSDI_norm",
       index = "NTSDI", support = "single reservoir storage record",
       site = "Nechako Reservoir", source = "7_WBInex_NTSDI_GeoLakes_LocalRS.r")
)
for (spec in ntsdi_specs) {
  f <- file.path(ntsdi_dir, spec$file)
  if (!file.exists(f)) {
    cat(sprintf("  ⚠ %s not found; skipped: %s\n", spec$index, f))
    next
  }
  tryCatch({
    x <- data.table::fread(f, data.table = FALSE)
    date_col <- intersect(c("date", "Date", "date_month"), names(x))[1]
    if (is.na(date_col) || !spec$col %in% names(x))
      stop("Expected date column and standardized column ", spec$col)
    add_other_series(x[[date_col]], x[[spec$col]], spec$index,
                     spec$support, spec$index, spec$site, spec$source)
    cat(sprintf("  ✓ %s: %d rows\n", spec$index, nrow(x)))
  }, error = function(e) cat(sprintf("  ⚠ %s skipped: %s\n", spec$index, e$message)))
}

# Streamflow SSI is gauge-specific and is not spatially averaged. Preserve the
# station and processing pipeline (B = naturalized; C = original/observed).
for (pipeline in c("ssi", "ssi_C")) {
  ssi_dir <- file.path(WD_PATH, "streamflow_results", pipeline, "monthly_data")
  rds_files <- list.files(ssi_dir, pattern = "_monthly_ssi\\.rds$", full.names = TRUE)
  if (!length(rds_files)) {
    cat(sprintf("  ⚠ No SSI monthly RDS files found under %s\n", ssi_dir))
    next
  }
  for (f in rds_files) {
    tryCatch({
      x <- readRDS(f)
      if (!all(c("date", "ssi") %in% names(x))) stop("Expected date and ssi columns")
      station_id <- sub("_monthly_ssi\\.rds$", "", basename(f))
      station_name <- station_id
      meta_file <- file.path(WD_PATH, "streamflow_results", "shared",
                             if (pipeline == "ssi_C") "stations_metadata_pipelineC.rds" else "stations_metadata.rds")
      if (file.exists(meta_file)) {
        meta <- tryCatch(readRDS(meta_file), error = function(e) NULL)
        if (!is.null(meta) && all(c("StationID", "StationName") %in% names(meta))) {
          hit <- which(as.character(meta$StationID) == station_id)
          if (length(hit)) station_name <- as.character(meta$StationName[hit[1]])
        }
      }
      pipe_label <- if (pipeline == "ssi_C") "Pipeline C observed/original flow" else "Pipeline B naturalized flow"
      add_other_series(x$date, x$ssi, "SSI", "hydrometric station", station_id,
                       station_name, pipe_label)
      cat(sprintf("  ✓ SSI %s (%s): %d rows\n", station_id, pipeline, nrow(x)))
    }, error = function(e) {
      cat(sprintf("  ⚠ SSI file %s skipped: %s\n", basename(f), e$message))
    })
  }
}

if (length(other_ts_rows)) {
  other_ts_all <- data.table::rbindlist(other_ts_rows, use.names = TRUE, fill = TRUE)
  data.table::fwrite(other_ts_all,
                     file.path(OTHER_TS_DIR, "other_hydrological_indices_long.csv"),
                     na = "NA")
  cat(sprintf("  ✓ Wrote %d rows across %d distinct series to %s\n",
              nrow(other_ts_all),
              length(unique(paste(other_ts_all$index, other_ts_all$site_id,
                                  other_ts_all$source, sep = "|"))),
              OTHER_TS_DIR))
} else {
  cat("  ⚠ No TSDI/NTSDI or station SSI source files were available.\n")
}

################################################################################
# PART 7 — SNOW-DROUGHT COUNTERPARTS FOR SPI/SPEI-ONLY FIGURE FAMILIES
#
# SPI/SPEI manuscript figures use six accumulation scales. Snow products in the
# current framework are SWEI-3 and SSPI-1, so they receive dedicated, equivalent
# chronology figures rather than being forced into incompatible six-scale grids.
# This section writes:
#   Fig2_MS_SnowDrought_Historical_2x1.pdf/.png       (full record)
#   Fig2S_MS_Recent_SnowDrought_2x1.pdf/.png          (2022–2025)
#   FigS_SnowDrought_Subbasins_Historical_7x2.pdf/.png (seven subbasins, full record)
# The recent 7×2 subbasin figure is created in the regional section above.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 7: Snow-drought counterparts for historical/recent figures\n")
cat("══════════════════════════════════════════════════════════════\n")

snow_basin_series <- list(
  "SWEI-3" = tryCatch(load_any_ts("swei", SWEI_SCALE), error = function(e) NULL),
  "SSPI-1" = tryCatch(load_any_ts("sspi", SSPI_SCALE), error = function(e) NULL)
)
snow_basin_series <- Filter(function(x) !is.null(x) && nrow(x) > 0L, snow_basin_series)

if (length(snow_basin_series)) {
  snow_basin_long <- data.table::rbindlist(lapply(names(snow_basin_series), function(lbl) {
    d <- as.data.frame(snow_basin_series[[lbl]])
    if (!all(c("date", "value") %in% names(d))) return(NULL)
    data.frame(date = as.Date(d$date), value = as.numeric(d$value),
               IndexPanel = lbl, stringsAsFactors = FALSE)
  }), use.names = TRUE, fill = TRUE)
  snow_basin_long <- snow_basin_long[is.finite(snow_basin_long$value) & !is.na(snow_basin_long$date), ]
  snow_basin_long$IndexPanel <- factor(snow_basin_long$IndexPanel,
                                       levels = c("SWEI-3", "SSPI-1"))
  
  if (nrow(snow_basin_long)) {
    save_snow_basin_chronology <- function(d, file_stub, title_text, subtitle_text) {
      if (!nrow(d)) {
        cat(sprintf("  ⚠ No data for %s; figure skipped.\n", file_stub))
        return(invisible(NULL))
      }
      finite_y <- d$value[is.finite(d$value)]
      ylo <- min(-3.5, floor(min(finite_y, na.rm = TRUE) * 1.05))
      yhi <- max( 3.5, ceiling(max(finite_y, na.rm = TRUE) * 1.05))
      cols <- c("SWEI-3" = "#1E90FF", "SSPI-1" = "#9467bd")
      g <- ggplot2::ggplot(d, ggplot2::aes(x = date, y = value, colour = IndexPanel)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.3) +
        ggplot2::geom_hline(yintercept = DROUGHT_ONSET, colour = "grey25",
                            linetype = "dashed", linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = -1, colour = "grey55",
                            linetype = "dotted", linewidth = 0.3) +
        ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
        ggplot2::scale_colour_manual(values = cols, name = "Snow index", drop = FALSE) +
        ggplot2::scale_x_date(date_breaks = if (min(d$date) >= as.Date("2020-01-01")) "1 year" else "10 years",
                              date_labels = "%Y", expand = ggplot2::expansion(mult = c(0.005, 0.005)) ) +
        ggplot2::coord_cartesian(ylim = c(ylo, yhi)) +
        ggplot2::facet_wrap(~IndexPanel, ncol = 1, scales = "fixed", drop = FALSE) +
        ggplot2::labs(title = title_text, subtitle = subtitle_text,
                      x = "Year", y = "Standardized drought-index value") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey35", size = 9),
          strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey65"),
          strip.text = ggplot2::element_text(face = "bold", size = 10),
          legend.position = "bottom",
          panel.grid.minor = ggplot2::element_blank()
        )
      for (ext in c("pdf", "png")) {
        out <- file.path(BASIN_PLOT_DIR, paste0(file_stub, ".", ext))
        tryCatch({
          ggplot2::ggsave(out, g, width = 12, height = 7.5, units = "in",
                          dpi = if (ext == "png") 300 else "print",
                          device = if (ext == "pdf") {
                            if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
                          } else "png",
                          limitsize = FALSE, bg = "white")
          cat(sprintf("  ✓ Saved: %s\n", basename(out)))
        }, error = function(e) cat(sprintf("  ⚠ %s: %s\n", basename(out), e$message)))
      }
      invisible(g)
    }
    
    save_snow_basin_chronology(
      snow_basin_long,
      "Fig2_MS_SnowDrought_Historical_2x1",
      "Historical snow-drought chronology for the Nechako Basin",
      "SWEI-3 and SSPI-1; dashed line: drought onset threshold; dotted line: −1"
    )
    snow_recent_basin <- snow_basin_long[
      snow_basin_long$date >= RECENT_START_4e & snow_basin_long$date <= RECENT_END_4e, , drop = FALSE
    ]
    save_snow_basin_chronology(
      snow_recent_basin,
      "Fig2S_MS_Recent_SnowDrought_2x1",
      "Recent snow-drought chronology for the Nechako Basin (2022–2025)",
      "SWEI-3 and SSPI-1; dashed line: drought onset threshold; dotted line: −1"
    )
  } else {
    cat("  ⚠ SWEI/SSPI inputs had no finite dated values; basin snow chronology figures skipped.\n")
  }
} else {
  cat("  ⚠ Neither SWEI-3 nor SSPI-1 basin series is available; snow chronology figures skipped.\n")
}

# Full-record subbasin counterpart to the existing 2022–2025 7×2 snow figure.
# Uses Script 11 regional products only; no basin or regional means are recomputed.
if (exists("regional_all", inherits = TRUE) && nrow(regional_all)) {
  tryCatch({
    snow_hist_regional <- data.table::copy(regional_all)[
      tolower(index_type) %in% c("swei", "sspi") &
        ((tolower(index_type) == "swei" & as.integer(scale) == SWEI_SCALE) |
           (tolower(index_type) == "sspi" & as.integer(scale) == SSPI_SCALE))
    ]
    if (!nrow(snow_hist_regional)) {
      cat("  ⚠ No regional SWEI-3/SSPI-1 history found; historical subbasin figure skipped.\n")
    } else {
      # Remove the whole-basin union, if present, retaining subbasin series only.
      whole_hist_id <- NULL
      if (exists("NECHAKO_REGION_BASIN_ID", inherits = TRUE)) {
        whole_hist_id <- get("NECHAKO_REGION_BASIN_ID", inherits = TRUE)
      }
      if (is.null(whole_hist_id) || !length(whole_hist_id)) {
        whole_hist_candidates <- unique(snow_hist_regional$region_id[
          grepl("whole|basin|union", snow_hist_regional$region_name, ignore.case = TRUE)
        ])
        if (length(whole_hist_candidates) == 1L) whole_hist_id <- whole_hist_candidates
      }
      if (!is.null(whole_hist_id) && length(whole_hist_id)) {
        snow_hist_regional <- snow_hist_regional[region_id != whole_hist_id[1L]]
      } else {
        snow_hist_regional <- snow_hist_regional[
          !grepl("whole basin|basin union|entire basin", region_name, ignore.case = TRUE)
        ]
      }
      hist_region_ids <- unique(as.character(snow_hist_regional$region_id))
      if (length(hist_region_ids) != 7L) {
        stop("Expected 7 subbasins for historical SWEI/SSPI figure; found ",
             length(hist_region_ids), ".")
      }
      hist_catalog_file <- file.path(REGIONAL_TS_DIR, "region_catalog.csv")
      if (file.exists(hist_catalog_file)) {
        hist_rc <- data.table::fread(hist_catalog_file)
        if (all(c("region_id", "region_name") %in% names(hist_rc))) {
          hist_rc_ids <- unique(as.character(hist_rc$region_id))
          hist_region_ids <- c(hist_rc_ids[hist_rc_ids %in% hist_region_ids],
                               setdiff(hist_region_ids, hist_rc_ids))
        }
      }
      snow_hist_regional[, region_id_chr := as.character(region_id)]
      snow_hist_regional[, region_id := factor(region_id_chr, levels = hist_region_ids)]
      hist_region_names <- unique(snow_hist_regional[, .(region_id_chr, region_name)])
      hist_region_names <- hist_region_names[match(hist_region_ids, hist_region_names$region_id_chr)]
      snow_hist_regional[, region_name_plot := factor(region_name,
                                                      levels = hist_region_names$region_name)]
      snow_hist_regional[, IndexPanel := ifelse(tolower(index_type) == "swei", "SWEI-3", "SSPI-1")]
      snow_hist_regional[, IndexPanel := factor(IndexPanel, levels = c("SWEI-3", "SSPI-1"))]
      snow_hist_colours <- c("SWEI-3" = "#1E90FF", "SSPI-1" = "#9467bd")
      p_snow_hist_subbasins <- ggplot2::ggplot(
        snow_hist_regional,
        ggplot2::aes(x = date, y = value, colour = IndexPanel, group = IndexPanel)
      ) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.25) +
        ggplot2::geom_hline(yintercept = DROUGHT_ONSET, colour = "grey25",
                            linewidth = 0.35, linetype = "dashed") +
        ggplot2::geom_hline(yintercept = -1, colour = "grey50", linewidth = 0.22,
                            linetype = "dotted") +
        ggplot2::geom_line(linewidth = 0.35, na.rm = TRUE) +
        ggplot2::scale_colour_manual(values = snow_hist_colours, name = "Snow index", drop = FALSE) +
        ggplot2::scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                              expand = ggplot2::expansion(mult = c(0.005, 0.005))) +
        ggplot2::facet_grid(rows = ggplot2::vars(region_name_plot),
                            cols = ggplot2::vars(IndexPanel), scales = "fixed", drop = FALSE) +
        ggplot2::labs(
          title = "Historical snow-drought evolution across Nechako subbasins",
          subtitle = "SWEI-3 and SSPI-1; dashed line: drought onset threshold; dotted line: −1",
          x = "Year", y = "Standardized index value"
        ) +
        ggplot2::theme_classic(base_size = 8) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 7, colour = "grey40", hjust = 0.5),
          strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70", linewidth = 0.25),
          strip.text.x = ggplot2::element_text(size = 8, face = "bold"),
          strip.text.y = ggplot2::element_text(size = 7, face = "bold"),
          axis.text.x = ggplot2::element_text(size = 5.5),
          axis.text.y = ggplot2::element_text(size = 5.5),
          axis.title = ggplot2::element_text(size = 7),
          axis.ticks = ggplot2::element_line(linewidth = 0.25),
          legend.position = "bottom", legend.text = ggplot2::element_text(size = 7),
          legend.title = ggplot2::element_text(size = 7),
          legend.key.size = ggplot2::unit(0.35, "cm"),
          panel.border = ggplot2::element_rect(colour = "grey65", fill = NA, linewidth = 0.35),
          panel.spacing = ggplot2::unit(0.18, "lines")
        )
      for (ext in c("pdf", "png")) {
        out <- file.path(REGIONAL_PLOT_DIR, paste0("FigS_SnowDrought_Subbasins_Historical_7x2.", ext))
        tryCatch({
          ggplot2::ggsave(out, p_snow_hist_subbasins, width = 16, height = 10.5,
                          units = "in", dpi = if (ext == "png") 300 else "print",
                          device = if (ext == "pdf") {
                            if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
                          } else "png", limitsize = FALSE, bg = "white")
          cat(sprintf("  ✓ Saved: %s\n", basename(out)))
        }, error = function(e) cat(sprintf("  ⚠ %s: %s\n", basename(out), e$message)))
      }
      cat("  ✓ Historical snow subbasin figure complete: 7 subbasins × 2 snow indices.\n")
    }
  }, error = function(e) {
    cat(sprintf("  ⚠ Historical snow subbasin figure failed: %s\n", e$message))
  })
} else {
  cat("  ⚠ Regional products unavailable; historical snow subbasin figure skipped.\n")
}

cat("\nScript 12 completed. Basin-mean indices and non-basin-support indices are kept in separate outputs.\n")

