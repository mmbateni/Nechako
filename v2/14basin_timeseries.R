# Pixel-to-basin correlation figures (Part 4c) map ANALYSIS_PIXEL_BASIN_CORR (utils);
# Script 15 Part 8 maps the same statistic for the supplementary indices and scales.
####################################################################################
# 14basin_timeseries.R  ·  BASIN-AVERAGED TIME SERIES — OBSERVED FRAMEWORK
# ─────────────────────────────────────────────────────────────────────────────────
# SCOPE  (authoritative home for all basin-averaged outputs):
#   PART 1 – Load producer basin-averaged CSVs for SPI, SPEI-PM,
#             SWEI, SSPI, soil SSMI / SSI-SM and optionally MSPI / MSPEI.
#             no detrended/counterfactual branch is used.
#   PART 2 – Use authoritative producer basin series directly (no duplicate flat CSVs).
#   PART 3 – Combined CSV (basin_mean_all_indices.csv — all indices long format).
#   PART 4 – Time series plots (one publication-quality PNG per index × scale).
#   PART 4d– Seasonality bar chart: mean monthly P vs PET(PM) (12 values each),
#             common Script 11 basin-average record; climatic water balance panel.
#             Reads Script 11 basin monthly series through
#             analysis_monthly_p_pet_climatology().
#             Output: FigS_Seasonality_P_PET.pdf/.png.
#   PART 4e– Full historical 12-panel chronology + recent SPI/SPEI-PM comparison.
#   PART 4e-C2– Manuscript 12-panel SPI + SPEI-PM chronology across all six scales.
#   PART 4e-D– Supplementary Fig. S1: 1- and 6-month historical chronology.
#             Fig2_MS_AllIndices_2x6.pdf/.png (12-panel: SPI and SPEI-PM × 1, 3, 6, 12, 24, 36 months).
#             Fig2S_MS_Recent_AllScales_3x2.pdf/.png (the configured recent period SPI + SPEI-PM comparison across all six scales).
#             The historical figure is the primary Results chronology for Section 4.2.
#             Section 4.2 of the manuscript.
#   PART 5 – Excel summary: descriptive statistics + drought event counts
#             (recent window) + pairwise correlation matrix (Drought_Summary.xlsx).
#   PART 6 – Optional support-aware exports for station SSI, surface-water
#             TSDI / NTSDI / SWSI, groundwater SGI (9SGI.R, well-network
#             median) and the multivariate MSDI (10MSDI.r, vine copula).
#             These are kept separate from basin-mean indices.
#
# PREREQUISITE: run 13trend_test_ALL.R first.
# NOTE: 15trends_visualization.R handles spatial, basin, and homogeneity diagnostics.
#       All basin-averaged time-series outputs live here.
#
# METHODOLOGY NOTE:
#   Step 1 — SPI/SPEI are computed per-pixel by 1/3SPI/SPEI_ERALand.R /
#             3SPEI_ERALand.R; results live in seasonal NetCDF files.
#   Step 2 — index producers area-weight each index on its native grid using
#             cell area and fractional basin coverage and write
#             {index}_{scale:02d}_basin_averaged_by_month.csv.
#   Step 3 — This script reads those authoritative CSVs.
#             Do NOT recompute averages here.
#
# WD_PATH / UTILS_FILE are set from getwd() at the top of this script.
rm(list = ls())
SCRIPT14_STARTED <- Sys.time()  # used to reject stale products in the closing inventory

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
# Shared analysis contracts are defined in DROUGHT_ANALYSIS_utils.R.


# Fail early if the framework utility file lacks a required definition. Spatial
# correlations and basin precipitation are read from producer CSVs below; no
# S12_BASIN_BW, s12_basin_mean_xy(), or basin_weighted_mean() is needed.
required_utils <- c(
  "CACHE_DIR", "TREND_DIR", "DROUGHT_END", "index_colours",
  "utils_load_packages", "load_basin_avg_csv", "get_event_min_duration",
  "identify_drought_events_ts", "shared_ts_theme", "days_in_month_vec",
  "load_nechako_subbasins", "nechako_subbasin_union", "load_basin_vect",
  "SPI_SCALES", "SPEI_SCALES", "SWEI_SCALE", "SSPI_SCALE",
  "ANALYSIS_PERIOD_LABEL", "ANALYSIS_RECENT_START", "ANALYSIS_RECENT_END",
  "ANALYSIS_PIXEL_BASIN_CORR", "ANALYSIS_PIXEL_BASIN_CORR_LABEL", "ANALYSIS_PIXEL_BASIN_CORR_TAG"
)
missing_utils <- required_utils[!vapply(required_utils, exists, logical(1),
                                        inherits = TRUE)]
if (length(missing_utils))
  stop("Required framework definitions missing from DROUGHT_ANALYSIS_utils.R: ",
       paste(missing_utils, collapse = ", "),
       "\nApply the utils patches (U1-U7) before running this script.")


# SPI_SCALES, SPEI_SCALES, SWEI_SCALE and SSPI_SCALE come from
# DROUGHT_ANALYSIS_utils.R; MANUSCRIPT_SCALES is supplied by the shared framework.
# Soil-moisture outputs retain their original native grid; no manuscript-scale promotion.
CONFIG13_SOIL <- file.path(TREND_DIR, "script13_analysis_config.rds")
if (!file.exists(CONFIG13_SOIL))
  stop("Missing required Script 13 config: ", CONFIG13_SOIL,
       "; run 13trend_test_ALL.R before Script 14.")
cfg_soil <- readRDS(CONFIG13_SOIL)
if (length(cfg_soil$soil_layer) != 1L || is.na(cfg_soil$soil_layer))
  stop("Script 13 config lacks a valid soil_layer: ", CONFIG13_SOIL)
SOIL_LAYER <- as.character(cfg_soil$soil_layer)
SOIL_SCALES <- as.integer(cfg_soil$soil_scales)
if (anyNA(SOIL_SCALES)) stop("Script 13 config contains invalid soil scales")
SSMI_SEAS_DIR <- file.path(WD_PATH, "ssmi_results_seasonal", SOIL_LAYER)
SSI_SM_SEAS_DIR <- file.path(WD_PATH, "ssi_sm_results_seasonal", SOIL_LAYER)

SWEI_SEAS_DIR <- file.path(WD_PATH, "swei_results_seasonal")
SSPI_SEAS_DIR <- file.path(WD_PATH, "sspi_results_seasonal")
if (!exists("DROUGHT_ONSET", inherits = TRUE))
  stop("DROUGHT_ONSET must be defined by DROUGHT_ANALYSIS_utils.R")

utils_load_packages(c("terra", "ggplot2", "lubridate", "dplyr",
                      "data.table", "patchwork", "openxlsx", "zoo",
                      "RColorBrewer", "scales", "sf"))


## ── Output directories ──────────────────────────────────────────────────────────
# BASIN_PLOT_DIR comes from utils (WD/temporal_drought/basin_averaged_plots).
for (d in c(BASIN_PLOT_DIR, CACHE_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# NATIVE ERA5-LAND GRID CONTRACT
#
# Script 14 uses the original/native ERA5-Land precipitation grid.
#
# The first layer of:
#   monthly_data_direct/total_precipitation_monthly.nc
#
# is the sole raster template. No canonical/project-wide analysis grid is
# loaded, generated, or required.
#
# Raster products are never reprojected/resampled merely to satisfy Script 14.
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

## ── Basin geometry and fractional-coverage weights ──────────────────────────────
## Basin = dissolved union of the 7 subbasins (cached loader in utils). Any basin
## mean computed in this script uses cell area x fractional coverage and the
## shared BASIN_MIN_VALID_AREA_COVERAGE (0.20) reporting rule.
basin_shp <- load_basin_vect()
if (!terra::same.crs(basin_shp, analysis_template)) {
  basin_shp <- terra::project(
    basin_shp,
    terra::crs(analysis_template)
  )
}

## ── MSPI / MSPEI single scale (mirror of 8 constant) ───────────────────────────
MSPI_MSPEI_SCALE <- 1L
# Colours are defined only in DROUGHT_ANALYSIS_utils.R; no local overrides.

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
## Date-safe copy of utils' drought_band_layers(): same seven bands, colours,
## alpha (0.10) and dotted class lines. utils draws the bands with
## annotation_custom(xmin = -Inf, xmax = Inf), whose numeric edges trigger one
## warning each on a Date x-axis (7 bands x 2 = 14 per plot). Here the edges are
## Date-typed -Inf/Inf, so the Date scale accepts them silently.
date_band_layers_14 <- function(alpha = 0.10) {
  bands <- data.frame(
    ymin = c(-Inf, -2.0, -1.5, -1.0, 1.0, 1.5, 2.0),
    ymax = c(-2.0, -1.5, -1.0,  1.0, 1.5, 2.0,  Inf),
    fill = c("#7B0025",   # D3 extremely dry
             "#D73027",   # D2 severely dry
             "#F46D43",   # D1 moderately dry
             "#2E7D32",   # N0 normal
             "#90EE90",   # W1 moderately wet
             "#66BB6A",   # W2 severely wet
             "#4575B4"),  # W3 extremely wet
    stringsAsFactors = FALSE)
  x_lo <- as.Date(-Inf, origin = "1970-01-01")
  x_hi <- as.Date( Inf, origin = "1970-01-01")
  c(lapply(seq_len(nrow(bands)), function(i)
    ggplot2::annotate("rect", xmin = x_lo, xmax = x_hi,
                      ymin = bands$ymin[i], ymax = bands$ymax[i],
                      fill = bands$fill[i], alpha = alpha, colour = NA)),
    list(ggplot2::geom_hline(yintercept = c(-2, -1.5, -1, 1, 1.5, 2),
                             linetype = "dotted", color = "gray50", linewidth = 0.3)))
}
## Plotting-warning log for this run. warnings() is not reliable here: it can show
## another script's leftover warning, and Source() holds back in-script warnings.
WARN_LOG_14 <- character()
log_warnings_14 <- function(expr, label) withCallingHandlers(expr, warning = function(w) {
  WARN_LOG_14 <<- c(WARN_LOG_14, sprintf("%s: %s", label, conditionMessage(w)))
  invokeRestart("muffleWarning")
})

## Save PDF + PNG and report only what was actually written.
save_pdf_png_14 <- function(plot, out_dir, stub, width, height, dpi = 300, bg = "white") {
  written <- character()
  for (ext in c("pdf", "png")) {
    out <- file.path(out_dir, paste0(stub, ".", ext))
    unlink(out, force = TRUE)
    ok <- tryCatch({
      log_warnings_14(label = basename(out), expr = if (ext == "pdf")
        ggplot2::ggsave(out, plot, width = width, height = height, units = "in",
                        limitsize = FALSE, bg = bg,
                        device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf)
        else
          ggplot2::ggsave(out, plot, width = width, height = height, units = "in",
                          limitsize = FALSE, bg = bg, dpi = dpi, device = "png"))
      file.exists(out) && isTRUE(file.info(out)$size > 0)
    }, error = function(e) {
      cat(sprintf("  ⚠ %s: %s\n", basename(out), conditionMessage(e))); FALSE
    })
    if (ok) { cat(sprintf("  ✓ Saved: %s\n", basename(out))); written <- c(written, out) }
    else cat(sprintf("  ⚠ NOT written: %s\n", basename(out)))
  }
  invisible(written)
}
legend_inside_14 <- function(pos) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    ggplot2::theme(legend.position = "inside", legend.position.inside = pos)
  else ggplot2::theme(legend.position = pos)
}

## Basin outline as plain lon/lat paths (avoids coord_sf break/label bug).
basin_path_df_14 <- function(basin_sf) {
  if (is.null(basin_sf) || !inherits(basin_sf, "sf")) return(NULL)
  tryCatch({
    g  <- sf::st_geometry(sf::st_transform(basin_sf, 4326))
    g  <- g[!sf::st_is_empty(g)]
    xy <- sf::st_coordinates(sf::st_cast(g, "MULTIPOLYGON"))
    lc <- grep("^L[0-9]+$", colnames(xy), value = TRUE)
    data.frame(X = xy[, "X"], Y = xy[, "Y"],
               grp = do.call(paste, c(lapply(lc, function(k) xy[, k]), sep = "_")))
  }, error = function(e) { cat("  ⚠ basin outline unavailable:", conditionMessage(e), "\n"); NULL })
}
load_supported_basin_csv <- function(index_type, scale) {
  idx <- tolower(index_type)
  if (idx == "sspi") return(load_basin_avg_csv(idx, scale, data_dir = SSPI_SEAS_DIR))
  if (idx == "ssmi") return(load_basin_avg_csv(idx, scale, data_dir = SSMI_SEAS_DIR))
  if (idx == "ssi_sm") return(load_basin_avg_csv(idx, scale, data_dir = SSI_SM_SEAS_DIR))
  load_basin_avg_csv(idx, scale)
}

## Load the authoritative basin-average CSV and return a data.frame
## (Date, Index, Mean_Value, Basin_Drought_Flag) ready for plotting.
## This is a basin-mean threshold flag (0/1), NOT a fractional drought area.
## Real areal fractions are in Script 13 basin_extent products.
##
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
    Basin_Drought_Flag = as.integer(as.numeric(df$value) < DROUGHT_ONSET),
    stringsAsFactors = FALSE
  )
}

## Load MSPI / MSPEI time series from the Excel file written by their own scripts.
## Returns data.frame(date, value) or NULL if the file is absent.
MULTI_TS_CACHE_14 <- new.env(parent = emptyenv())
load_mspi_mspei_ts <- function(index_type) {
  stopifnot(index_type %in% c("mspi", "mspei"))
  if (exists(index_type, envir = MULTI_TS_CACHE_14, inherits = FALSE))
    return(get(index_type, envir = MULTI_TS_CACHE_14, inherits = FALSE)$data)
  results_dir <- file.path(WD_PATH,
                           if (index_type == "mspi") "mspi_results" else "mspei_results")
  pat  <- sprintf("%s_basin_timeseries.*\\.xlsx$", index_type)
  fls  <- list.files(results_dir, pattern = pat, full.names = TRUE)
  if (!length(fls)) {
    cat(sprintf("  ⚠ No Excel timeseries found for %s in %s\n",
                toupper(index_type), results_dir))
    assign(index_type, list(data = NULL), envir = MULTI_TS_CACHE_14)
    return(NULL)
  }
  result <- tryCatch({
    df       <- openxlsx::read.xlsx(fls[1], sheet = "Basin_Timeseries")
    df$date  <- as.Date(paste0(df$Date, "-01"))
    df$value <- df$Basin_Mean
    df[, c("date", "value")]
  }, error = function(e) {
    cat(sprintf("  ⚠ Could not read %s Excel: %s\n", toupper(index_type), e$message))
    NULL
  })
  assign(index_type, list(data = result), envir = MULTI_TS_CACHE_14)
  result
}

## Unified loader — returns data.frame(date, value) for any supported index.
## Used by the producer-series plot builder and Excel export.
##
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

cat("\n── SSMI / SSI-SM (soil moisture; native scales) ──\n")
for (idx in c("ssmi", "ssi_sm")) {
  for (sc in SOIL_SCALES) {
    ts <- tryCatch(load_basin_index(idx, sc), error = function(e) {
      cat(sprintf("  Soil-moisture %s-%02d unavailable: %s\n", toupper(idx), sc, e$message)); NULL
    })
    if (!is.null(ts)) {
      ts$Index <- sprintf("%s%d", toupper(idx), sc)
      all_data[[ts$Index[1L]]] <- ts
    }
  }
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
      Basin_Drought_Flag = as.integer(df_multi$value < DROUGHT_ONSET),
      stringsAsFactors = FALSE)
    all_data[[label]] <- ts_multi
    cat(sprintf("  ✓ %s: %d time steps\n", label, nrow(ts_multi)))
  }
}

if (!length(all_data))
  stop("No basin-averaged drought-index data loaded. Run Script 13 (13trend_test_ALL.R) and Script 4 (4snow_drought_indices_ERALand.R) first, and verify their output paths.")

################################################################################
# PART 2 – AUTHORITATIVE PRODUCER SERIES (NO RE-EXPORT COPIES)
################################################################################
# The producer CSVs loaded in Part 1 remain the only persisted source.
# Downstream plotting and summary functions call load_any_ts() directly.

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

## Compatibility interface for the rest of Script 14: load the *producer* file,
## not a separately cached or re-exported copy.
load_producer_ts <- function(idx, sc) {
  df <- load_any_ts(idx, sc)
  if (is.null(df) || !nrow(df))
    stop(sprintf("Producer basin series unavailable: %s-%02d", idx, as.integer(sc)))
  df$date <- as.Date(df$date)
  df$value <- as.numeric(df$value)
  as.data.frame(df)
}

## Publication-quality basin-average time series plot.
## Drought-fill ribbon, categorical band shading, event count subtitle.
make_basin_ts_plot <- function(df, idx, sc) {
  idx_lc    <- tolower(idx)
  idx_label <- if (idx_lc == "spei") "SPEI-PM" else toupper(idx)
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
    date_band_layers_14() +
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
    ggplot2::scale_x_date(breaks = seq(NECHAKO_RECORD_START, NECHAKO_RECORD_END, by = "5 years"), date_labels = "%Y") +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s-%02d  Basin-Averaged Time Series (Area-Weighted)",
                         idx_label, sc),
      subtitle = sprintf("%d drought events  (onset < %.2f, end >= %.2f, min %d months)",
                         n_ev, DROUGHT_ONSET, DROUGHT_END, min_dur),
      x = "Year",
      y = sprintf("%s-%02d Index Value", idx_label, sc))
}

cat("\n  SPI, SPEI, SWEI, SSPI, SSMI & SSI-SM...\n")
for (idx in c("spi", "spei", "swei", "sspi", "ssmi", "ssi_sm")) {
  scales <- switch(idx,
                   spi = SPI_SCALES,
                   spei = SPEI_SCALES,
                   swei = SWEI_SCALE,
                   sspi = SSPI_SCALE,
                   ssmi = SOIL_SCALES,
                   ssi_sm = SOIL_SCALES
  )
  for (sc in scales) {
    tryCatch({
      df <- load_producer_ts(idx, sc)
      out <- file.path(BASIN_PLOT_DIR, sprintf("%s_%02d_timeseries.png", idx, sc))
      log_warnings_14(ggplot2::ggsave(out, make_basin_ts_plot(df, idx, sc),
                                      width = 12, height = 6, dpi = 150),
                      sprintf("%s-%02d", toupper(idx), sc))
      cat(sprintf("  ✓ %s-%02d\n", toupper(idx), sc))
    }, error = function(e) {
      cat(sprintf("  ⚠ %s-%02d: %s\n", toupper(idx), sc, e$message))
    })
  }
}

cat("\n  MSPI & MSPEI (if available)...\n")
for (idx in c("mspi", "mspei")) {
  tryCatch({
    df  <- load_producer_ts(idx, MSPI_MSPEI_SCALE)
    out <- file.path(BASIN_PLOT_DIR,
                     sprintf("%s_%02d_timeseries.png", idx, MSPI_MSPEI_SCALE))
    ggplot2::ggsave(out, make_basin_ts_plot(df, idx, MSPI_MSPEI_SCALE),
                    width = 12, height = 6, dpi = 150)
    cat(sprintf("  ✓ %s-%02d\n", toupper(idx), MSPI_MSPEI_SCALE))
  }, error = function(e)
    cat(sprintf("  ⚠ %s-%02d: %s\n", toupper(idx), MSPI_MSPEI_SCALE, e$message)))
}

################################################################################
# PART 4b – MANUSCRIPT FIGURE 4 (SPEI-PM, scales 1/3/6)
#
# SCALE CONTRACT:
#   Manuscript SPEI figures here use scales 1, 3, and 6 only. Scale 2 is not
#   part of the authoritative manuscript index set; basin-average panels
#   load existing CSV products rather than recomputing spatial means.
#
# Output (BASIN_PLOT_DIR): Fig4_MS_SPEI_136_PM_timeseries.pdf/.png
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 4b: Manuscript Figure 4 — SPEI-PM-1/3/6 time series\n")
cat("══════════════════════════════════════════════\n")

# ── Directory paths ─────────────────────────────────────────────────────────────
SPEI_PM_DIR  <- file.path(WD_PATH, "spei_results_seasonal")
SPI_SEAS_DIR <- file.path(WD_PATH, "spi_results_seasonal")   # used by Part 4c + 4e
# ── Shared panel constants ──────────────────────────────────────────────────────
THR_EVENT       <- DROUGHT_ONSET
THR_MOD         <- -1.0
HIGHLIGHT_START <- ANALYSIS_RECENT_START
HIGHLIGHT_END   <- ANALYSIS_RECENT_END
RECENT_LABEL_14 <- sprintf("%d–%d", as.integer(format(HIGHLIGHT_START, "%Y")),
                           as.integer(format(HIGHLIGHT_END, "%Y")))
SNOW_LAB_SWEI <- sprintf("SWEI-%d", SWEI_SCALE)
SNOW_LAB_SSPI <- sprintf("SSPI-%d", SSPI_SCALE)

# Match the seven standardized classes and colours in drought_band_layers()
# (DROUGHT_ANALYSIS_utils.R). The right-closed dry intervals also match the
# boundary convention in detect_drought_events(): x == -2 is extremely dry.
fill_pal <- c(
  "Extremely dry"  = "#7B0025",
  "Severely dry"   = "#D73027",
  "Moderately dry" = "#F46D43",
  "Normal"         = "#2E7D32",
  "Moderately wet" = "#90EE90",
  "Severely wet"   = "#66BB6A",
  "Extremely wet"  = "#4575B4"
)

classify_drought_band_14 <- function(x) {
  cut(x,
      breaks = c(-Inf, -2, -1.5, -1, 1, 1.5, 2, Inf),
      labels = names(fill_pal), right = TRUE, include.lowest = TRUE)
}

# ── SPEI-PM panel builder ─────────────────────────────────────────────────────
# is_first : TRUE = show legend; FALSE = suppress
make_spei_panel_v2 <- function(sc, panel_lab,
                               show_xlab = FALSE, is_first = FALSE) {
  col_line  <- unname(index_colours[["spei"]])
  pet_label <- "SPEI-PM (Penman-Monteith)"
  
  # Read the authoritative producer basin-average CSV directly (no fallback).
  df <- tryCatch(load_producer_ts("spei", sc), error = function(e) NULL)
  
  if (is.null(df) || !nrow(df))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("No data\nSPEI-PM-%d", sc),
                               size = 4.5, colour = "grey50") +
             ggplot2::theme_void())
  
  df$date  <- as.Date(df$date)
  df$value <- as.numeric(df$value)
  
  df_fill <- df
  # End of each month's rectangle is computed on the full monthly record first,
  # so every rectangle spans exactly one month even if a gap exists.
  df_fill$date_end <- df_fill$date + days_in_month_vec(df_fill$date) 
  # Then drop months with no index value (leading accumulation months of
  # SPEI-3/6): no NA-coloured rectangles and no "NA" key in the legend.
  df_fill <- df_fill[is.finite(df_fill$value), , drop = FALSE]
  df_fill$fill_cat <- classify_drought_band_14(df_fill$value)
  
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
        ggplot2::guide_legend(title.position = "top", nrow = 2,
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
      colour = col_line, linewidth = 0.55, inherit.aes = FALSE,
      na.rm = TRUE) +
    ggplot2::annotate("text",
                      x = min(df$date) + 365, y = yhi * 0.88,
                      label = panel_lab, size = 3.2, fontface = "bold",
                      hjust = 0, vjust = 1, colour = "grey10") +
    ggplot2::annotate("text",
                      x = seq(as.Date(sprintf("%04d-01-01", ANALYSIS_PERIOD_START)),
                              by = "year", length.out = 2L)[2L], y = THR_EVENT + 0.08,
                      label = sprintf("x₀ = %.2f", DROUGHT_ONSET), size = 2.4,
                      hjust = 0, colour = "grey30") +
    ggplot2::annotate("text",
                      x = as.Date((as.numeric(HIGHLIGHT_START) + as.numeric(HIGHLIGHT_END)) / 2, origin = "1970-01-01"), y = yhi * 0.97,
                      label = RECENT_LABEL_14, size = 2.4, hjust = 0.5,
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
    unlink(out, force = TRUE)
    tryCatch({
      ggplot2::ggsave(out, fig, width = 7.0, height = 9.0, units = "in",
                      dpi = if (ext == "png") 300 else "print",
                      device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png")
      if (!file.exists(out) || !isTRUE(file.info(out)$size > 0)) stop("Output missing or empty")
      cat(sprintf("  ✓ Saved: %s\n", basename(out)))
    }, error = function(e) cat(sprintf("  ⚠ %s: %s\n", ext, e$message)))
  }
  invisible(fig)
}

common_subtitle <- paste0(
  "Fill: drought category  |  ",
  sprintf("Dashed: event threshold (x₀ = %.2f)  |  ", THR_EVENT),
  "Dotted: −1.0 reference (NOT the event threshold)  |  ",
  paste0("Grey band: ", RECENT_LABEL_14, "  |  "),
  "Recent-period shading shown in every panel")

tryCatch({
  # ── SPEI-PM only ─────────────────────────────────────────────────────────
  cat("  Building SPEI-PM (Penman-Monteith)...\n")
  pm_panels <- list(
    make_spei_panel_v2(1, "(a)", show_xlab = FALSE, is_first = TRUE),
    make_spei_panel_v2(3, "(b)", show_xlab = FALSE, is_first = FALSE),
    make_spei_panel_v2(6, "(c)", show_xlab = TRUE,  is_first = FALSE)
  )
  save_spei_fig(pm_panels,
                sprintf("Basin-averaged SPEI-PM (Penman-Monteith PET) — Nechako River Basin (%s)", ANALYSIS_PERIOD_LABEL),
                common_subtitle,
                "Fig4_MS_SPEI_136_PM_timeseries")
  cat("  ✓ SPEI-PM figure complete\n")
}, error = function(e) cat(sprintf("  ⚠ SPEI-PM figure failed: %s\n", e$message)))

################################################################################
# PART 4c – SPATIAL REPRESENTATIVENESS FIGURE
# Output: Fig4c_<statistic>_SPEI_PM_SPI.pdf/.png (statistic = ANALYSIS_PIXEL_BASIN_CORR)
# Rows 1-2: 2x3 maps -- SPEI-PM scales 1, 3, 6 (Blues) and SPI scales 1, 3, 6 (Greens)
# Bottom row: (g) SPEI-PM and (h) SPI pixel-to-basin correlation densities
# The correlation CSVs are written by 13trend_test_ALL.R (section 13).
# All panels use ggplot2 + patchwork (no terra::plot).
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat(sprintf("PART 4c: Spatial representativeness — SPEI-PM + SPI %s maps\n", ANALYSIS_PIXEL_BASIN_CORR_LABEL))
cat("══════════════════════════════════════════════════════════════\n")

## Helper: detect basin boundary and return sf object in specified CRS
detect_basin_sf <- function(coords_df) {
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

# Pixel-to-basin correlations are descriptive; each pixel contributes to the
# area-weighted basin mean (small positive self-inclusion component).
## Helper: build one correlation map panel ─────────────────────────────
## Plot producer correlation coordinates in their native geographic CRS.
#
# Common signed colour domain for all correlation panels: [-1, 1].
# This is intentional for cross-index comparability, not a data-dependent scale.
make_map_panel <- function(r_v, coords, basin_sf, title_str, subtitle_str,
                           show_legend = TRUE, legend_title = CORR_LAB) {
  df_points <- data.frame(x = coords$lon, y = coords$lat, r = as.numeric(r_v))
  df_points <- df_points[is.finite(df_points$r) & is.finite(df_points$x) & is.finite(df_points$y), ]
  if (!nrow(df_points)) return(ggplot2::ggplot() + ggplot2::theme_void())
  if (any(abs(df_points$x) > 180 | abs(df_points$y) > 90))
    stop("make_map_panel requires geographic longitude/latitude coordinates.")
  
  snap <- 0.01
  df_points$x_r <- round(df_points$x / snap) * snap
  df_points$y_r <- round(df_points$y / snap) * snap
  df_grid <- aggregate(r ~ x_r + y_r, data = df_points, FUN = mean)
  if (!nrow(df_grid)) return(ggplot2::ggplot() + ggplot2::theme_void())
  
  ux <- sort(unique(df_grid$x_r)); uy <- sort(unique(df_grid$y_r))
  dx <- if (length(ux) >= 2L) stats::median(diff(ux)) else snap
  dy <- if (length(uy) >= 2L) stats::median(diff(uy)) else snap
  if (!is.finite(dx) || dx <= 0) dx <- snap
  if (!is.finite(dy) || dy <= 0) dy <- snap
  
  x_lim <- range(df_grid$x_r); y_lim <- range(df_grid$y_r)
  xlim_p <- x_lim + c(-1, 1) * max(dx * 0.55, diff(x_lim) * 0.03)
  ylim_p <- y_lim + c(-1, 1) * max(dy * 0.55, diff(y_lim) * 0.03)
  
  # Only breaks inside the panel (out-of-range breaks caused the original error).
  x_breaks <- seq(-128, -123, by = 1); x_breaks <- x_breaks[x_breaks >= xlim_p[1] & x_breaks <= xlim_p[2]]
  y_breaks <- seq(53, 56, by = 1);     y_breaks <- y_breaks[y_breaks >= ylim_p[1] & y_breaks <= ylim_p[2]]
  
  p_map <- ggplot2::ggplot(df_grid, ggplot2::aes(x = x_r, y = y_r, fill = r)) +
    ggplot2::geom_tile(width = dx, height = dy)
  bpath <- basin_path_df_14(basin_sf)
  if (!is.null(bpath))
    p_map <- p_map + ggplot2::geom_path(data = bpath, ggplot2::aes(X, Y, group = grp),
                                        inherit.aes = FALSE, colour = "black", linewidth = 0.35)
  p_map +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5), name = legend_title,
      oob = scales::squish, na.value = "transparent",
      guide = if (show_legend) ggplot2::guide_colourbar(
        title.position = "top", barheight = grid::unit(2.8, "cm"),
        barwidth = grid::unit(0.35, "cm"), frame.colour = "grey50",
        ticks.colour = "grey40") else "none") +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = x_breaks,
                                labels = function(x) paste0(abs(x), "\u00b0W")) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = y_breaks,
                                labels = function(y) paste0(y, "\u00b0N")) +
    # Plain lon/lat cartesian frame with latitude-corrected aspect; no raster is
    # reprojected or resampled.
    ggplot2::coord_fixed(ratio = 1 / cos(mean(ylim_p) * pi / 180),
                         xlim = xlim_p, ylim = ylim_p, expand = FALSE) +
    ggplot2::labs(title = title_str, subtitle = subtitle_str,
                  x = "Longitude (\u00b0W)", y = "Latitude (\u00b0N)") +
    theme_map4c +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 6.5, margin = ggplot2::margin(t = 2)),
      axis.title.y = ggplot2::element_text(size = 6.5, margin = ggplot2::margin(r = 2)),
      axis.text.x  = ggplot2::element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
      axis.text.y  = ggplot2::element_text(size = 6),
      axis.ticks   = ggplot2::element_line(linewidth = 0.3, colour = "grey40"),
      axis.ticks.length = ggplot2::unit(2, "pt"))
}
## ── Part 4c / 4c-2 shared helpers (global scope) ─────────────────────────────
## Correlations are calculated once by 13trend_test_ALL.R (section 13) against each
## producer's fractional-coverage basin series; this script only reads them.

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

## Helper: build one correlation-density panel ──────────────────────────────────
make_density_panel <- function(r_list, combo_names, col_pal,
                               panel_lab, title_str) {
  available <- vapply(r_list, function(rv) !is.null(rv) && any(is.finite(rv)), logical(1))
  r_list <- r_list[available]
  combo_names <- combo_names[available]
  col_pal <- col_pal[combo_names]
  df  <- do.call(rbind, lapply(seq_along(r_list), function(k) {
    rv <- r_list[[k]]
    data.frame(r = rv, combo = combo_names[k], stringsAsFactors = FALSE)
  }))
  if (is.null(df) || !nrow(df))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = "No data", size = 4, colour = "grey50") +
             ggplot2::theme_void())
  df      <- df[is.finite(df$r), ]
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
  
  ## Keep median lines at the true estimate; linetype separates coincident lines.
  med_df$xint_nudged <- med_df$med_r  # median lines must show actual estimates
  
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
    ggplot2::scale_x_continuous(breaks = seq(-1.0, 1.0, by = 0.2)) +
    ggplot2::coord_cartesian(xlim = c(-1.0, 1.01)) +
    ## RESERVE ~60% OF TOP PANEL FOR LABELS to guarantee no curve overlap
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.70))) +
    ## NOTE: subtitle deliberately omitted here — a single shared note is added
    ##       once via the patchwork caption in the fig_dens assembly below,
    ##       preventing the identical text from appearing in both panels and
    ##       overflowing into the adjacent panel's space.
    ggplot2::labs(
      title = sprintf("%s  %s", panel_lab, title_str),
      x     = paste(CORR_LAB, "(pixel vs basin mean)"),
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

## Correlation statistic shared with Script 15 (set in DROUGHT_ANALYSIS_utils.R).
CORR_LAB <- ANALYSIS_PIXEL_BASIN_CORR_LABEL
CORR_TAG <- ANALYSIS_PIXEL_BASIN_CORR_TAG
load_pixel_corr <- function(idx, sc, dir, stat = ANALYSIS_PIXEL_BASIN_CORR) {
  f <- file.path(dir, sprintf("%s_%02d_pixel_basin_correlation.csv", idx, as.integer(sc)))
  if (!file.exists(f)) stop("Pixel-basin correlation not found: ", f,
                            "\n  Re-run 13trend_test_ALL.R (section 13 writes these files).")
  d <- data.table::fread(f)
  required_cols <- c("x", "y", stat)
  if (!all(required_cols %in% names(d)))
    stop("Invalid pixel-basin correlation CSV: ", f, "; missing column(s): ",
         paste(setdiff(required_cols, names(d)), collapse = ", "))
  list(r = as.numeric(d[[stat]]),
       coords = data.frame(lon = as.numeric(d$x), lat = as.numeric(d$y)))
}

tryCatch({
  ## 1. Read SPEI-PM correlations for scales 1, 3, 6 ─────────────────
  
  spei_scales <- c(1L, 3L, 6L)
  spei_r_list <- vector("list", 3L)
  spei_obj_list <- vector("list", 3L)
  spei_med_list <- rep(NA_real_, 3L)
  spei_pct80_list <- rep(NA_real_, 3L)
  spei_plots <- vector("list", 3L)
  spei_labels <- c("(a)", "(b)", "(c)")
  
  for (k in seq_along(spei_scales)) {
    sc <- spei_scales[k]
    cat(sprintf(" Loading SPEI-PM-%d...\n", sc))
    
    obj <- tryCatch(load_pixel_corr("spei", sc, SPEI_PM_DIR),
                    error = function(e) { cat("   ", conditionMessage(e), "\n"); NULL })
    if (is.null(obj)) next
    rv <- obj$r
    
    spei_obj_list[[k]] <- obj
    spei_r_list[[k]] <- rv
    spei_med_list[k] <- median(rv, na.rm = TRUE)
    spei_pct80_list[k] <- 100 * mean(rv >= 0.80, na.rm = TRUE)
  }
  
  # Shared fixed [-1, 1] correlation colour scale for all panels.
  
  # Build SPEI-PM map panels with the common fixed colour domain
  for (k in seq_along(spei_scales)) {
    obj <- spei_obj_list[[k]]
    rv  <- spei_r_list[[k]]
    
    if (is.null(obj) || is.null(rv)) next
    
    bsf <- detect_basin_sf(obj$coords)
    
    spei_plots[[k]] <- make_map_panel(
      rv, obj$coords, bsf,
      title_str = sprintf("%s SPEI-%d (PM)", spei_labels[k], spei_scales[k]),
      subtitle_str = sprintf("median = %.3f\n%.0f%% \u2265 0.80",
                             spei_med_list[k], spei_pct80_list[k]),
      show_legend = TRUE,
      legend_title = paste(CORR_LAB, "(SPEI-PM)")
    )
  }
  
  
  ## 2. Read SPI correlations for scales 1, 3, 6 ─────────────────────
  
  spi_scales <- c(1L, 3L, 6L)
  spi_r_list <- vector("list", 3L)
  spi_obj_list <- vector("list", 3L)
  spi_med_list <- rep(NA_real_, 3L)
  spi_pct80_list <- rep(NA_real_, 3L)
  spi_plots <- vector("list", 3L)
  spi_labels <- c("(d)", "(e)", "(f)")
  
  for (k in seq_along(spi_scales)) {
    sc <- spi_scales[k]
    cat(sprintf(" Loading SPI-%d...\n", sc))
    
    
    obj <- tryCatch(load_pixel_corr("spi", sc, SPI_SEAS_DIR),
                    error = function(e) { cat("   ", conditionMessage(e), "\n"); NULL })
    if (is.null(obj)) next
    rv <- obj$r
    
    spi_obj_list[[k]] <- obj
    spi_r_list[[k]] <- rv
    spi_med_list[k] <- median(rv, na.rm = TRUE)
    spi_pct80_list[k] <- 100 * mean(rv >= 0.80, na.rm = TRUE)
  }
  
  # Shared fixed [-1, 1] correlation colour scale for all panels.
  
  # Build SPI map panels with the common fixed colour domain
  for (k in seq_along(spi_scales)) {
    obj <- spi_obj_list[[k]]
    rv  <- spi_r_list[[k]]
    
    if (is.null(obj) || is.null(rv)) next
    
    bsf <- detect_basin_sf(obj$coords)
    
    spi_plots[[k]] <- make_map_panel(
      rv, obj$coords, bsf,
      title_str = sprintf("%s SPI-%d", spi_labels[k], spi_scales[k]),
      subtitle_str = sprintf("median = %.3f\n%.0f%% ≥ 0.80",
                             spi_med_list[k], spi_pct80_list[k]),
      show_legend = TRUE,
      legend_title = paste(CORR_LAB, "(SPI)")
    )
  }
  
  
  any_spei <- any(!sapply(spei_plots, is.null))
  any_spi  <- any(!sapply(spi_plots,  is.null))
  if (!any_spei && !any_spi) { cat("  No data for any combination\n"); stop("no_data") }
  
  ## 3. Density panel (g): SPEI-PM values for scales 1, 3, 6 ────────────
  col_spei_pm <- stats::setNames(
    grDevices::colorRampPalette(c("#08306b", unname(index_colours[["spei"]]), "#9ecae1"))(3L),
    c("SPEI-PM-1", "SPEI-PM-3", "SPEI-PM-6"))
  pc_pm   <- make_density_panel(
    r_list      = spei_r_list,
    combo_names = names(col_spei_pm),
    col_pal     = col_spei_pm,
    panel_lab   = "(g)",
    title_str   = paste("Pixel", CORR_LAB, "— Penman-Monteith PET"))
  
  ## 4. Density panel (h): SPI values for scales 1, 3, 6 ────────────────
  col_spi <- stats::setNames(
    grDevices::colorRampPalette(c("#99000d", unname(index_colours[["spi"]]), "#fcbba1"))(3L),
    c("SPI-1", "SPI-3", "SPI-6"))
  pl_spi   <- make_density_panel(
    r_list      = spi_r_list,
    combo_names = names(col_spi),
    col_pal     = col_spi,
    panel_lab   = "(h)",
    title_str   = paste("Pixel", CORR_LAB, "— SPI (scales 1, 3, 6)"))
  
  ## 5. Assemble figure ────────────────────────────────────────────────────
  top_plots <- c(spei_plots, spi_plots)
  top_plots <- lapply(top_plots, function(p) {
    if (is.null(p))
      ggplot2::ggplot() + ggplot2::theme_void()
    else p
  })
  
  ## Map legend is shown on selected panels (override any defaults)
  fig_maps <- patchwork::wrap_plots(top_plots, nrow = 2, ncol = 3, guides = "collect") &
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical"
    )
  
  ## Density panels keep their legends.
  ## The shared subtitle note ("Vertical lines = median …") is placed ONCE as
  ## a patchwork caption below both panels — this prevents the identical text
  ## from appearing in each individual panel and overflowing into its neighbour.
  fig_dens <- (pc_pm | pl_spi) +
    patchwork::plot_layout(guides = "keep") +
    patchwork::plot_annotation(
      caption = "Vertical lines = median per scale  (solid / dashed / dotted)",
      theme   = ggplot2::theme(
        plot.caption = ggplot2::element_text(
          size   = 7.5, colour = "grey40",
          hjust  = 0.5,
          margin = ggplot2::margin(t = 2, b = 2))))
  
  ## Combine: maps (no legends) + density plots (with their legends)
  fig4c  <- fig_maps / fig_dens +
    patchwork::plot_layout(heights = c(2, 0.7)) +
    patchwork::plot_annotation(
      title    = paste("Pixel-to-basin-mean", CORR_LAB, "— SPEI 1,3,6 + SPI 1,3,6 — Nechako River Basin"),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5,
                                              margin = ggplot2::margin(b = 4)),
        plot.subtitle = ggplot2::element_blank(),
        legend.position = "bottom",
        legend.text     = ggplot2::element_text(size = 10),
        legend.key.size = ggplot2::unit(0.45, "cm")))
  
  ## Save 
  for (ext in c("pdf", "png")) {
    out_f   <- file.path(BASIN_PLOT_DIR,
                         sprintf("Fig4c_%s_SPEI_PM_SPI.%s", CORR_TAG, ext))
    unlink(out_f, force = TRUE)
    tryCatch({
      ggplot2::ggsave(out_f, fig4c,
                      width  = 6.85,
                      height = 9.21,
                      units  = "in",
                      dpi    = if (ext == "png") 600 else "print",
                      device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png",
                      limitsize = FALSE)
      if (!file.exists(out_f) || !isTRUE(file.info(out_f)$size > 0)) stop("Output missing or empty")
      cat(sprintf("  ✓ Saved: %s\n", basename(out_f)))
    }, error = function(e) { unlink(out_f, force = TRUE); cat(sprintf("  ⚠ %s: %s\n", ext, e$message)) })
  }
  
  ## Console summary
  cat(sprintf("\n  Representativeness summary (%s):\n", CORR_LAB))
  cat(sprintf("  %-14s  median    pct≥0.80\n", "Combination"))
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
# PART 4c-2 – SPEI-PM 1/3/6 SPATIAL REPRESENTATIVENESS FIGURE (PM only)
# Output: Fig4c_<statistic>_SPEI136_PM.pdf/.png
# Row 1 (a-c): the SPEI-PM maps built in Part 4c; row 2 (d): their density panel.
# Nothing is re-read: the maps and correlations come from Part 4c.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4c-2: SPEI-PM 1/3/6 spatial representativeness (reuses Part 4c)\n")
cat("══════════════════════════════════════════════════════════════\n")

tryCatch({
  if (!exists("spei_plots") || !exists("col_spei_pm") ||
      all(vapply(spei_plots, is.null, logical(1L)))) stop("no_data")
  col_spei_pm2 <- stats::setNames(unname(col_spei_pm), sprintf("SPEI-%d PM", spei_scales))
  pc_pm2 <- make_density_panel(
    r_list = spei_r_list, combo_names = names(col_spei_pm2), col_pal = col_spei_pm2,
    panel_lab = "(d)",
    title_str = paste("Pixel", CORR_LAB, "\u2014 Penman-Monteith PET (scales 1, 3, 6)"))
  top_plots_pm2 <- lapply(spei_plots, function(p)
    if (is.null(p)) ggplot2::ggplot() + ggplot2::theme_void() else p)
  fig4c_pm2 <- patchwork::wrap_plots(top_plots_pm2, nrow = 1, ncol = 3) / pc_pm2 +
    patchwork::plot_layout(heights = c(1.6, 0.7), guides = "collect") +
    patchwork::plot_annotation(
      title = paste("Pixel-to-basin-mean", CORR_LAB, "\u2014 SPEI-PM 1, 3, 6 \u2014 Nechako River Basin"),
      theme = ggplot2::theme(
        plot.title      = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5,
                                                margin = ggplot2::margin(b = 4)),
        legend.position = "bottom",
        legend.text     = ggplot2::element_text(size = 10),
        legend.key.size = ggplot2::unit(0.45, "cm")))
  for (ext in c("pdf", "png")) {
    out_pm2 <- file.path(BASIN_PLOT_DIR, sprintf("Fig4c_%s_SPEI136_PM.%s", CORR_TAG, ext))
    unlink(out_pm2, force = TRUE)
    tryCatch({
      ggplot2::ggsave(out_pm2, fig4c_pm2, width = 6.85, height = 6.50, units = "in",
                      dpi = if (ext == "png") 600 else "print",
                      device = if (ext == "pdf") { if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf } else "png",
                      limitsize = FALSE)
      if (!file.exists(out_pm2) || !isTRUE(file.info(out_pm2)$size > 0)) stop("Output missing or empty")
      cat(sprintf("  \u2713 Saved: %s\n", basename(out_pm2)))
    }, error = function(e) { unlink(out_pm2, force = TRUE); cat(sprintf("  \u26a0 %s: %s\n", ext, e$message)) })
  }
}, error = function(e) {
  if (!identical(conditionMessage(e), "no_data"))
    cat(sprintf("  \u26a0 Part 4c-2 failed: %s\n", e$message))
})

cat("\n══ PART 4c-2 complete ══\n")

################################################################################
# PART 4c-3 — SNOW-INDEX SPATIAL REPRESENTATIVENESS
#
# Counterpart to the SPI/SPEI pixel-to-basin maps. SWEI-3 and SSPI-1 use the
# correlation CSVs written by 13trend_test_ALL.R (section 13) into each index's
# seasonal folder, exactly like SPI/SPEI: Script 13 correlates every basin pixel
# with the producer's basin series. Nothing is resampled or recomputed here.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4c-3: Spatial representativeness — SWEI-3 + SSPI-1\n")
cat("══════════════════════════════════════════════════════════════\n")

# The CSVs contain basin cells only, so the producers' basin-cell support is kept.

tryCatch({
  snow_specs <- list(
    list(index = "swei", scale = SWEI_SCALE, dir = SWEI_SEAS_DIR,
         label = SNOW_LAB_SWEI, colour = "#1E90FF"),
    list(index = "sspi", scale = SSPI_SCALE, dir = SSPI_SEAS_DIR,
         label = SNOW_LAB_SSPI, colour = "#9467bd")
  )
  snow_r_list <- vector("list", length(snow_specs))
  names(snow_r_list) <- vapply(snow_specs, `[[`, character(1L), "label")
  snow_obj_list <- vector("list", length(snow_specs))
  snow_map_list <- vector("list", length(snow_specs))
  snow_medians <- rep(NA_real_, length(snow_specs))
  snow_pct80 <- rep(NA_real_, length(snow_specs))
  
  for (k in seq_along(snow_specs)) {
    spec <- snow_specs[[k]]
    cat(sprintf("  Loading %s correlation CSV (Script 13)...\n", spec$label))
    obj <- tryCatch(load_pixel_corr(spec$index, spec$scale, spec$dir),
                    error = function(e) {
                      cat(sprintf("    ⚠ %s unavailable: %s\n", spec$label, e$message))
                      NULL
                    })
    if (is.null(obj)) next
    rv <- obj$r
    snow_obj_list[[k]] <- obj
    snow_r_list[[k]] <- rv
    snow_medians[k] <- median(rv, na.rm = TRUE)
    snow_pct80[k] <- 100 * mean(rv >= 0.80, na.rm = TRUE)
  }
  
  if (!any(vapply(snow_r_list, function(x) !is.null(x) && any(is.finite(x)), logical(1L)))) {
    cat("  ⚠ Neither snow index could produce valid pixel-to-basin correlations; figure skipped.\n")
  } else {
    for (k in seq_along(snow_specs)) {
      spec <- snow_specs[[k]]
      obj <- snow_obj_list[[k]]
      rv <- snow_r_list[[k]]
      if (is.null(obj) || is.null(rv)) next
      snow_map_list[[k]] <- make_map_panel(
        rv, obj$coords, detect_basin_sf(obj$coords),
        title_str = sprintf("(%s) %s pixel-to-basin correlation", letters[k], spec$label),
        subtitle_str = sprintf("median = %.3f; %.0f%% ≥ 0.80", snow_medians[k], snow_pct80[k]),
        show_legend = (k == 1L),
        legend_title = CORR_LAB
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
      title_str = paste("Pixel", CORR_LAB, "— snow-drought indices")
    )
    snow_spatial_fig <- patchwork::wrap_plots(snow_map_panels, nrow = 1L, ncol = 2L,
                                              guides = "collect") / snow_density +
      patchwork::plot_layout(heights = c(1.5, 0.7), guides = "collect") +
      patchwork::plot_annotation(
        title = paste("Pixel-to-basin-mean", CORR_LAB, "— snow-drought indices"),
        subtitle = paste0("Native-grid ", SNOW_LAB_SWEI, " and ", SNOW_LAB_SSPI,
                          "; availability-dependent snow-pixel support; regional means not recomputed"),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
          plot.subtitle = ggplot2::element_text(size = 9, hjust = 0.5, colour = "grey35"),
          legend.position = "bottom"
        )
      )
    for (ext in c("pdf", "png")) {
      out <- file.path(BASIN_PLOT_DIR, sprintf("Fig4c_%s_SWEI_SSPI.%s", CORR_TAG, ext))
      unlink(out, force = TRUE)
      tryCatch({
        ggplot2::ggsave(out, snow_spatial_fig, width = 10.5, height = 7.8, units = "in",
                        dpi = if (ext == "png") 600 else "print",
                        device = if (ext == "pdf") {
                          if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
                        } else "png", limitsize = FALSE, bg = "white")
        if (!file.exists(out) || !isTRUE(file.info(out)$size > 0)) stop("Output missing or empty")
        cat(sprintf("  ✓ Saved: %s\n", basename(out)))
      }, error = function(e) { unlink(out, force = TRUE); cat(sprintf("  ⚠ %s: %s\n", basename(out), e$message)) })
    }
    cat("  Snow-index representativeness summary:\n")
    for (k in seq_along(snow_specs)) {
      if (is.null(snow_r_list[[k]])) next
      cat(sprintf("    %-8s median = %.3f; %.0f%% of pixels ≥ 0.80\n",
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
#              exactly 12 bars per variable averaged over all complete years.
#   Panel (b): climatic water balance P − PET(PM)  [mm month⁻¹].
#
# DATA
#   analysis_monthly_p_pet_climatology() loads the authoritative Script 11
#   basin monthly P and PM-PET series and computes monthly climatologies
#   from their paired, complete calendar years (mm/month).
#
# OUTPUT  FigS_Seasonality_P_PET.pdf + .png  in BASIN_PLOT_DIR
################################################################################
cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("PART 4d: Seasonality figure — P vs PET(PM) (12 monthly means)\n")
cat("═══════════════════════════════════════════════════════════════════\n")

tryCatch({
  
  # Use the shared Script 11 basin series and computation.
  clim_4d <- analysis_monthly_p_pet_climatology()
  yr_range_4d <- attr(clim_4d, "years")
  COL_P_4d <- "#4393c3"; COL_PET_4d <- "#d7301f"
  cm_prec_4d <- data.frame(month = clim_4d$month, month_lab = clim_4d$Month,
                           mean_val = clim_4d$P_Mean, sd_val = clim_4d$P_SD,
                           n_yrs = clim_4d$n_yrs, variable = "Precipitation (P)")
  cm_pet_4d <- data.frame(month = clim_4d$month, month_lab = clim_4d$Month,
                          mean_val = clim_4d$PET_Mean, sd_val = clim_4d$PET_SD,
                          n_yrs = clim_4d$n_yrs, variable = "PET(PM)")
  ann_P <- round(sum(clim_4d$P_Mean))
  ann_PET <- round(sum(clim_4d$PET_Mean))
  ai <- round(ann_PET / ann_P, 2)
  nyrs <- min(clim_4d$n_yrs)
  cat(sprintf("  Common climatology: %d complete paired years (%d-%d); annual P=%d, PET=%d mm; Budyko dryness index PET/P=%.2f\n",
              nyrs, yr_range_4d[1L], yr_range_4d[2L], ann_P, ann_PET, ai))
  
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
        nyrs, "-year mean (", yr_range_4d[1L], "–", yr_range_4d[2L], ")  |  Error bars = ±1 SD  |  ",
        "Annual totals:  P = ", ann_P, " mm;  PET(PM) = ", ann_PET, " mm"),
      x = NULL) +
    ggplot2::theme_classic(base_size = 9.5) +
    ggplot2::theme(
      legend.background  = ggplot2::element_rect(fill = "white", colour = NA),
      legend.text        = ggplot2::element_text(size = 8.5),
      legend.key.size    = ggplot2::unit(0.40, "cm"),
      plot.title         = ggplot2::element_text(size = 9.5, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 7.5, colour = "grey40"),
      axis.text          = ggplot2::element_text(size = 8.5),
      axis.title.y       = ggplot2::element_text(size = 8.5),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      plot.margin        = ggplot2::margin(4, 6, 2, 4))
  p4d_a <- p4d_a + legend_inside_14(c(0.86, 0.88))
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
        "Each bar = mean of ", nyrs, " years (", yr_range_4d[1L], "–", yr_range_4d[2L], ").  ",
        "Annual totals: P = ", ann_P, " mm;  PET(PM) = ", ann_PET, " mm;  ",
        "Budyko dryness index (PET/P) = ", ai, "."),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 11, face = "bold", hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey35",
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 4))))
  
  save_pdf_png_14(fig_seas_4d, BASIN_PLOT_DIR, "FigS_Seasonality_P_PET", width = 7.2, height = 7.0)
  
}, error = function(e) {
  cat(sprintf("  ⚠ Part 4d failed: %s\n", e$message))
})

################################################################################
################################################################################
# PART 4e – MANUSCRIPT RESULTS FIGURES: FULL HISTORICAL DROUGHT CHRONOLOGY
#
# Purpose:
#   Provide the historical time-series evidence needed to compare the
#   the configured recent period drought with all prior drought episodes across the complete
#   record and across all six manuscript accumulation scales.
#
# Primary manuscript figure (built by build_historical_chronology_4e() in
# Part 4e-B below, as the first entry of CHRONOLOGY_VARIANTS_4E):
#   Fig2_MS_AllIndices_2x6.pdf/.png
#       12 panels = 2 index families × 6 accumulation scales
#       Row 1: SPI-1, SPI-3, SPI-6, SPI-12, SPI-24, SPI-36
#       Row 2: SPEI-PM-1, SPEI-PM-3, SPEI-PM-6, SPEI-PM-12, SPEI-PM-24, SPEI-PM-36
#       Full chronology; the configured recent period is highlighted.
#
# Secondary focused figure:
#   Fig2S_MS_Recent_AllScales_3x2.pdf/.png
#       Six panels = accumulation scales 1, 3, 6, 12, 24, 36 months.
#       Each panel overlays SPI and SPEI-PM for the configured recent period.
#
# Regional focused figure:
#   Fig2S_MS_Recent_Subbasins_7x6.pdf/.png
#       Seven rows = the seven supplied subbasins; six columns = accumulation
#       scales. Each panel overlays SPI and SPEI-PM for the configured recent period. The
#       whole-basin union is excluded.
#
# Scientific interpretation:
#   These are descriptive basin-average chronologies. The the configured recent period highlight
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
MANUSCRIPT_TS_SCALES <- as.integer(MANUSCRIPT_SCALES)

THR_EVENT_4e       <- DROUGHT_ONSET
THR_MOD_4e         <- -1.0
HIGHLIGHT_START_4e <- HIGHLIGHT_START
HIGHLIGHT_END_4e   <- HIGHLIGHT_END
RECENT_START_4e    <- HIGHLIGHT_START
RECENT_END_4e      <- as.Date(format(HIGHLIGHT_END, "%Y-%m-01"))

fill_pal_4e <- fill_pal  # Same standardized classes and colours as Part 4b.

LINE_COLOURS_4e <- c(
  "SPI"       = unname(index_colours[["spi"]]),
  "SPEI-PM"   = unname(index_colours[["spei"]])
)

# ── Unified historical time-series loader ─────────────────────────────────────
# IMPORTANT: this function only reads existing products. It does not recompute
# spatial means from NetCDF files; each requested scale must have its
# existing producer CSV or the panel is explicitly unavailable.
load_manuscript_ts_4e <- function(index_type, scale) {
  index_type <- match.arg(index_type, c("spi", "spei_pm"))
  
  if (index_type == "spi") {
    df <- tryCatch(load_producer_ts("spi", scale), error = function(e) NULL)
  } else if (index_type == "spei_pm") {
    df <- tryCatch(load_producer_ts("spei", scale), error = function(e) NULL)
  }
  if (is.null(df) || !nrow(df))
    stop(sprintf("No basin-average time series for %s scale %d", index_type, scale))
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
    spei_pm  = sprintf("SPEI-PM-%d", sc)
  )
  
  line_col <- unname(LINE_COLOURS_4e[
    switch(index_type, spi = "SPI", spei_pm = "SPEI-PM")
  ])
  
  # Retain the full historical record. The the configured recent period band is an annotation,
  # not a data filter.
  df_fill <- df
  df_fill <- df_fill[is.finite(df_fill$value), , drop = FALSE]
  df_fill$fill_cat <- classify_drought_band_14(df_fill$value)
  df_fill$date_end <- df_fill$date + days_in_month_vec(df_fill$date)  # always one calendar month  
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
          nrow = 2,
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
      x = seq(as.Date(sprintf("%04d-01-01", ANALYSIS_PERIOD_START)),
              by = "year", length.out = 2L)[2L],
      y = THR_EVENT_4e + 0.08,
      label = sprintf("%.2f", DROUGHT_ONSET),
      size = 2.2,
      hjust = 0,
      colour = "grey30"
    ) +
    ggplot2::annotate(
      "text",
      x = as.Date((as.numeric(HIGHLIGHT_START_4e) + as.numeric(HIGHLIGHT_END_4e)) / 2, origin = "1970-01-01"),
      y = yhi * 0.97,
      label = RECENT_LABEL_14,
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

# ── Recent six-panel formulation comparison ───────────────────────────────────
# Each panel is one accumulation scale. This figure is deliberately restricted
# to the configured recent period and overlays SPI with SPEI-PM.
# It complements, rather than replaces, the full historical
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
  
  # Manuscript-facing recent comparison contains SPI and SPEI-PM.
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
      x = seq(as.Date(format(HIGHLIGHT_START, "%Y-%m-01")),
              by = "month", length.out = 2L)[2L],
      y = yhi * 0.92,
      label = panel_lab,
      size = 3.0,
      fontface = "bold",
      hjust = 0,
      colour = "grey10"
    ) +
    ggplot2::annotate(
      "text",
      x = seq(as.Date(format(HIGHLIGHT_START, "%Y-%m-01")),
              by = "month", length.out = 2L)[2L],
      y = THR_EVENT_4e + 0.08,
      label = sprintf("%.2f", DROUGHT_ONSET),
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
  cat(paste0("  Building six-panel ", RECENT_LABEL_14, " formulation-comparison figure...\n"))
  
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
      title = paste0(RECENT_LABEL_14, " drought evolution across accumulation scales"),
      subtitle = paste0(
        "Each panel overlays SPI and SPEI-PM.  ",
        sprintf("Dashed line: drought-event threshold (x₀ = %.2f); ", DROUGHT_ONSET),
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
  
  save_pdf_png_14(fig_recent, BASIN_PLOT_DIR, "Fig2S_MS_Recent_AllScales_3x2",
                  width = 12, height = 7.5)
  
  cat("  ✓ Six-panel recent comparison complete\n")
}, error = function(e) {
  cat(sprintf("  ⚠ Recent six-panel figure failed: %s\n", e$message))
})

cat("\n══ PART 4e complete ══\n")

################################################################################
# PART 4e-B / C / C2 / D – CONFIGURATION-DRIVEN HISTORICAL CHRONOLOGIES
# Every product, including the 12-panel Fig2_MS_AllIndices_2x6, uses
# make_historical_panel_4e(), a single grid composer, and the same nonfatal
# error handling. Missing inputs skip only the affected figure.
################################################################################

build_historical_chronology_4e <- function(index_keys, scales, ncol, file_stub,
                                           width = 12, height = 9.5,
                                           title = NULL, subtitle = NULL,
                                           repeat_headers = FALSE) {
  index_prefix <- c(spi = "SPI", spei_pm = "SPEI_PM")
  index_labels <- c(spi = "SPI", spei_pm = "SPEI-PM")
  if (any(!index_keys %in% names(index_prefix)) ||
      any(!scales %in% MANUSCRIPT_TS_SCALES) || ncol < 1L)
    stop("Invalid index family, scale, or column count for ", file_stub)
  specification <- expand.grid(scale = as.integer(scales), index_type = index_keys,
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # expand.grid varies its FIRST argument fastest, preserving the manuscript
  # row order: each index family by ascending requested accumulation scale.
  specification$key <- paste0(unname(index_prefix[specification$index_type]), specification$scale)
  missing <- setdiff(specification$key, names(all_data))
  if (length(missing)) {
    cat(sprintf("  [SKIP] %s: missing producer series: %s\n",
                file_stub, paste(missing, collapse = ", ")))
    return(invisible(NULL))
  }
  n <- nrow(specification)
  if (n %% ncol != 0L) stop("Number of panels must be divisible by ncol: ", file_stub)
  nrow_grid <- n %/% ncol
  panels <- lapply(seq_len(n), function(i) {
    row_i <- (i - 1L) %/% ncol + 1L
    col_i <- (i - 1L) %% ncol + 1L
    idx <- specification$index_type[i]
    sc <- specification$scale[i]
    row_header <- col_i == 1L && (i == 1L ||
                                    specification$index_type[i - 1L] != idx)
    p <- make_historical_panel_4e(
      sc = sc, panel_lab = sprintf("(%s)", letters[i]),
      index_type = idx, show_xlab = row_i == nrow_grid,
      show_ylabel = FALSE, is_first = i == 1L)
    top_header <- isTRUE(repeat_headers) || row_i == 1L
    p + ggplot2::labs(title = if (top_header) paste0(sc, " mo") else NULL,
                      y = if (row_header) unname(index_labels[idx]) else NULL) +
      ggplot2::theme(
        plot.title = if (top_header)
          ggplot2::element_text(size = 13, face = "bold", hjust = 0.5,
                                margin = ggplot2::margin(b = 1)) else ggplot2::element_blank(),
        axis.title.y = if (row_header)
          ggplot2::element_text(size = 13, face = "bold", angle = 90,
                                margin = ggplot2::margin(r = 2)) else ggplot2::element_blank(),
        axis.text.x = if (row_i == nrow_grid) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
        axis.ticks.x = if (row_i == nrow_grid) ggplot2::element_line() else ggplot2::element_blank(),
        axis.title.x = if (row_i == nrow_grid) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 8),
        plot.margin = ggplot2::margin(0, 2, if (row_i == nrow_grid) 1 else 0, 1),
        legend.text = ggplot2::element_text(size = 9),
        legend.title = ggplot2::element_text(size = 9.5))
  })
  figure <- patchwork::wrap_plots(panels, ncol = ncol, nrow = nrow_grid,
                                  guides = "collect") +
    patchwork::plot_annotation(
      title = title, subtitle = subtitle,
      theme = ggplot2::theme(
        legend.position = "bottom", legend.text = ggplot2::element_text(size = 9),
        legend.title = ggplot2::element_text(size = 9.5),
        plot.title = ggplot2::element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 7.5, hjust = 0.5)))
  for (fmt in c("pdf", "png")) {
    target <- file.path(BASIN_PLOT_DIR, paste0(file_stub, ".", fmt))
    tryCatch({
      if (fmt == "pdf") ggplot2::ggsave(target, plot = figure, width = width,
                                        height = height, units = "in", limitsize = FALSE,
                                        device = if (isTRUE(capabilities("cairo")))
                                          grDevices::cairo_pdf else grDevices::pdf)
      else ggplot2::ggsave(target, plot = figure, width = width, height = height,
                           units = "in", dpi = 300, limitsize = FALSE)
      cat(sprintf("  [OK] %s\n", target))
    }, error = function(e)
      cat(sprintf("  [WARN] Could not write %s: %s\n", target, conditionMessage(e))))
  }
  invisible(figure)
}

# Each variant is independent: one absent series or graphics-device error does
# not interrupt Part 5 or the other manuscript/supplementary layouts.
CHRONOLOGY_VARIANTS_4E <- list(
  list(name = "Fig2_MS_AllIndices_2x6",
       indices = c("spi", "spei_pm"), scales = as.integer(MANUSCRIPT_TS_SCALES),
       ncol = 6L, width = 16, height = 7.0,
       title = sprintf("Basin-averaged drought-index chronology — Nechako River Basin (%s)", ANALYSIS_PERIOD_LABEL),
       subtitle = sprintf(paste0("Rows: SPI, SPEI-PM  |  Columns: 1-, 3-, 6-, 12-, 24-, 36-month scales  |  ",
                                 paste0("Dashed line: event threshold (x₀ = %.2f)  |  Grey band: ", RECENT_LABEL_14)), DROUGHT_ONSET)),
  list(name = "Fig2_MS_AllIndices_2x4_Scales3_12_24_36",
       indices = c("spi", "spei_pm"), scales = c(3L, 12L, 24L, 36L),
       ncol = 4L, width = 14, height = 7.0),
  list(name = "Fig2_MS_AllIndices_2x3_Scales1_3_6",
       indices = c("spi", "spei_pm"), scales = c(1L, 3L, 6L),
       ncol = 3L, width = 12, height = 7.0),
  list(name = "Fig2_MS_AllIndices_2x3_Scales12_24_36",
       indices = c("spi", "spei_pm"), scales = c(12L, 24L, 36L),
       ncol = 3L, width = 12, height = 7.0),
  list(name = "Fig3_MS_SPI_SPEIPM_AllScales_4x3",
       indices = c("spi", "spei_pm"), scales = as.integer(MANUSCRIPT_TS_SCALES),
       ncol = 3L, width = 12, height = 12, repeat_headers = TRUE,
       title = sprintf("Basin-average drought-index chronology — SPI and SPEI-PM (%s)", ANALYSIS_PERIOD_LABEL),
       subtitle = paste0("All six accumulation scales; dashed line = event threshold; grey band = ", RECENT_LABEL_14)),
  list(name = "FigS1_MS_AllIndices_2x2_Scales1_6",
       indices = c("spi", "spei_pm"), scales = c(1L, 6L),
       ncol = 2L, width = 12, height = 7.0)
)
for (cfg in CHRONOLOGY_VARIANTS_4E) {
  cat(sprintf("\nPart 4e: %s\n", cfg$name))
  tryCatch(log_warnings_14(label = cfg$name, expr = build_historical_chronology_4e(
    index_keys = cfg$indices, scales = cfg$scales,
    ncol = cfg$ncol, file_stub = cfg$name,
    width = cfg$width, height = cfg$height,
    title = cfg$title, subtitle = cfg$subtitle,
    repeat_headers = isTRUE(cfg$repeat_headers))),
    error = function(e) cat(sprintf("  [WARN] %s skipped: %s\n", cfg$name, conditionMessage(e))))
}

################################################################################
# PART 5 – EXCEL SUMMARY (statistics + drought event counts + correlations)
################################################################################
# Surface accumulated ggplot warnings from Parts 4/4e; they were previously
# only reported as "50 or more warnings" by R at script exit.
cat("\n[PART 4 WARNING AUDIT] Plotting warnings raised by this run:\n")
if (length(WARN_LOG_14)) {
  tab <- sort(table(WARN_LOG_14), decreasing = TRUE)
  for (k in seq_along(tab)) cat(sprintf("  ⚠ %s  (x%d)\n", names(tab)[k], tab[[k]]))
} else cat("  ✓ none\n")
cat("\n══════════════════════════════════════════════\n")
cat("PART 5: Excel summary\n")
cat("══════════════════════════════════════════════\n")

## Build one stats row for a given index × scale.
RECENT_EVENT_COL <- sprintf("Drought_Events_%s_%s",
                            format(ANALYSIS_RECENT_START, "%Y"), format(ANALYSIS_RECENT_END, "%Y"))
make_stats_row <- function(idx, sc) {
  df <- tryCatch(load_any_ts(idx, sc), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  # Events are detected on the FULL record (so events that began earlier keep their
  # true duration) and counted when they overlap the shared recent window.
  ev_count <- tryCatch({
    ev <- identify_drought_events_ts(df$date, df$value, scale = sc)
    sum(ev$end_date >= ANALYSIS_RECENT_START & ev$start_date <= ANALYSIS_RECENT_END)
  }, error = function(e) NA_integer_)
  out <- data.frame(
    Timescale = sprintf("%02d", sc),
    Index     = if (idx == "spei") "SPEI-PM" else if (idx == "ssi_sm") "SSI-SM" else toupper(idx),
    Mean      = round(mean(df$value,   na.rm = TRUE), 4),
    Median    = round(median(df$value, na.rm = TRUE), 4),
    StdDev    = round(sd(df$value,     na.rm = TRUE), 4),
    Min       = round(min(df$value,    na.rm = TRUE), 4),
    Max       = round(max(df$value,    na.rm = TRUE), 4),
    Drought_Events_Recent = ev_count,
    stringsAsFactors = FALSE)
  names(out)[names(out) == "Drought_Events_Recent"] <- RECENT_EVENT_COL
  out
}

export_basin_excel <- function(output_file) {
  wb        <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  
  ## ── Sheet 1: Summary Statistics ─────────────────────────────────────────────
  stats_rows <- list()
  for (idx in c("spi", "spei", "swei", "sspi", "ssmi", "ssi_sm")) {
    scales <- if (idx == "spi") SPI_SCALES else if (idx == "swei") SWEI_SCALE else if (idx == "sspi") SSPI_SCALE else if (idx %in% c("ssmi", "ssi_sm")) SOIL_SCALES else SPEI_SCALES
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
    
    # Snow drought indices: included where authoritative basin averages exist.
    ts <- tryCatch(load_any_ts("swei", SWEI_SCALE), error = function(e) NULL)
    if (!is.null(ts)) series_list[[sprintf("SWEI%d", SWEI_SCALE)]] <- ts
    ts <- tryCatch(load_any_ts("sspi", SSPI_SCALE), error = function(e) NULL)
    if (!is.null(ts)) series_list[[sprintf("SSPI%d", SSPI_SCALE)]] <- ts
    
    # Soil-moisture indices: match the population in Summary_Statistics.
    for (soil_idx in c("ssmi", "ssi_sm")) for (soil_sc in SOIL_SCALES) {
      soil_ts <- tryCatch(load_any_ts(soil_idx, soil_sc), error = function(e) NULL)
      if (!is.null(soil_ts))
        series_list[[sprintf("%s%d", if (soil_idx == "ssi_sm") "SSI-SM" else "SSMI", soil_sc)]] <- soil_ts
    }
    
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
# REGIONAL EXTENSION — 7 SUBBASINS + WHOLE NECHAKO BASIN
#
# The original basin outputs remain untouched. This section reads the
# authoritative regional CSVs produced by 13trend_test_ALL.R and creates
# region-aware time-series tables and figures without recomputing spatial means.
################################################################################
REGIONAL_TS_DIR <- file.path(TREND_DIR, "regional_timeseries")
REGIONAL_PLOT_DIR <- file.path(TREND_DIR, "regional_timeseries_plots")
dir.create(REGIONAL_TS_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(REGIONAL_PLOT_DIR, recursive=TRUE, showWarnings=FALSE)

# Authoritative dissolved-basin exclusion and stable seven-subbasin order.
# The same contract applies to regional SPI/SPEI and recent/historical snow panels.
select_subbasins_14 <- function(dt, label = "regional figure") {
  dt <- data.table::copy(data.table::as.data.table(dt))
  whole_id <- if (exists("NECHAKO_REGION_BASIN_ID", inherits = TRUE))
    as.character(get("NECHAKO_REGION_BASIN_ID", inherits = TRUE))[1L] else NA_character_
  if (is.na(whole_id) || !nzchar(whole_id)) {
    candidates <- unique(as.character(dt$region_id)[grepl("whole|basin|union", dt$region_name,
                                                          ignore.case = TRUE)])
    if (length(candidates) == 1L) whole_id <- candidates
  }
  if (!is.na(whole_id) && nzchar(whole_id))
    dt <- dt[as.character(region_id) != whole_id]
  else dt <- dt[!grepl("whole basin|basin union|entire basin", region_name, ignore.case = TRUE)]
  ids <- unique(as.character(dt$region_id))
  if (length(ids) != 7L)
    stop(label, ": expected 7 subbasins after excluding whole-basin union; got ", length(ids))
  catalog <- file.path(REGIONAL_TS_DIR, "region_catalog.csv")
  if (file.exists(catalog)) {
    rc <- data.table::fread(catalog)
    if ("region_id" %in% names(rc)) {
      catalog_ids <- unique(as.character(rc$region_id))
      ids <- c(catalog_ids[catalog_ids %in% ids], setdiff(ids, catalog_ids))
    }
  }
  list(data = dt, ids = ids)
}

regional_index_types <- c("spi", "spei", "sri", "swei", "sspi", "ssmi", "ssi_sm")
regional_index_labels <- c(spi = "SPI", spei = "SPEI-PM", sri = "SRI",
                           swei = "SWEI", sspi = "SSPI", ssmi = "SSMI", ssi_sm = "SSI-SM")
regional_scales <- MANUSCRIPT_SCALES
# Non-soil products are immediate children; only the configured soil layer
# is eligible in its own subdirectory. Never read stale alternative layers.
regional_files <- c(
  list.files(REGIONAL_TS_DIR, pattern = "_regional_averaged_by_month\\.csv$",
             full.names = TRUE, recursive = FALSE, ignore.case = TRUE),
  list.files(file.path(REGIONAL_TS_DIR, SOIL_LAYER),
             pattern = "_regional_averaged_by_month\\.csv$",
             full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
)

regional_all <- data.table::rbindlist(lapply(regional_files, function(f) {
  x <- data.table::fread(f)
  is_soil_path <- basename(dirname(f)) == SOIL_LAYER
  if (!"index_type" %in% names(x)) return(NULL)
  if (any(tolower(x$index_type) %in% c("ssmi", "ssi_sm")) != is_soil_path)
    return(NULL)
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
  #   accumulation scales. Every panel overlays SPI and SPEI-PM for the configured recent period.
  #
  # Scientific/data contract:
  #   * Reads the authoritative regional CSVs produced by 13trend_test_ALL.R.
  #   * Does not recompute regional spatial means.
  #   * Excludes the whole-basin dissolved union.
  ################################################################################
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("REGIONAL RECENT FIGURE: 7 subbasins × 6 accumulation scales\n")
  cat(paste0("                 SPI + SPEI-PM, ", RECENT_LABEL_14, "\n"))
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
      
      selected <- select_subbasins_14(regional_recent, "SPI/SPEI recent subbasin figure")
      regional_recent <- selected$data
      subbasin_ids <- selected$ids
      
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
          breaks = c(-3, -2, -1, 0, 1, 2),
          expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::coord_cartesian(ylim = c(-3.2, 2.5)) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(region_name_plot),
          cols = ggplot2::vars(scale_plot),
          scales = "fixed"
        ) +
        ggplot2::labs(
          title = paste0(RECENT_LABEL_14, " drought evolution across Nechako subbasins"),
          subtitle = paste0(
            "Each panel overlays SPI and SPEI-PM; rows = seven supplied subbasins; ",
            sprintf("columns = accumulation scale. Dashed line: x₀ = %.2f; dotted line: −1.0.", DROUGHT_ONSET)
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
      
      save_pdf_png_14(fig_subbasins, REGIONAL_PLOT_DIR, "Fig2S_MS_Recent_Subbasins_7x6",
                      width = 16, height = 14)
      
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
      #   overlays SPI and SPEI-PM for the configured recent period.
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
          breaks = c(-3, -2, -1, 0, 1, 2),
          expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::coord_cartesian(ylim = c(-3.2, 2.5)) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(region_name_plot),
          cols = ggplot2::vars(scale_plot),
          scales = "fixed"
        ) +
        ggplot2::labs(
          title = paste0(RECENT_LABEL_14, " drought evolution across Nechako subbasins"),
          subtitle = paste0(
            "SPI and SPEI-PM at 3-, 12-, and 36-month accumulation scales; ",
            sprintf("dashed line: x₀ = %.2f; dotted line: −1.0.", DROUGHT_ONSET)
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
      
      save_pdf_png_14(fig_subbasins_7x3, REGIONAL_PLOT_DIR,
                      "Fig2S_MS_Recent_Subbasins_7x3_Scales3_12_36",
                      width = 12, height = 14)
      
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
  # Regional means are read from Script 13 outputs; this block does not recompute
  # any spatial averages.
  ################################################################################
  cat("\n══════════════════════════════════════════════════════════════════\n")
  cat("REGIONAL RECENT FIGURE: snow drought across 7 subbasins\n")
  cat(paste0("                         ", SNOW_LAB_SWEI, " + ", SNOW_LAB_SSPI,
             ", ", RECENT_LABEL_14, "\n"))
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
      selected <- select_subbasins_14(snow_recent, "Recent snow subbasin figure")
      snow_recent <- selected$data
      snow_region_ids <- selected$ids
      
      snow_recent[, region_id := factor(as.character(region_id), levels = snow_region_ids)]
      snow_region_labels <- unique(snow_recent[, .(region_id, region_name)])
      snow_region_labels <- snow_region_labels[match(snow_region_ids, as.character(region_id))]
      snow_recent[, region_name_plot := factor(region_name,
                                               levels = snow_region_labels$region_name)]
      snow_recent[, IndexPanel := ifelse(tolower(index_type) == "swei", SNOW_LAB_SWEI, SNOW_LAB_SSPI)]
      snow_recent[, IndexPanel := factor(IndexPanel, levels = c(SNOW_LAB_SWEI, SNOW_LAB_SSPI))]
      
      snow_colours <- stats::setNames(c("#1E90FF", "#9467bd"), c(SNOW_LAB_SWEI, SNOW_LAB_SSPI))
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
        ggplot2::scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2),
                                    expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        ggplot2::coord_cartesian(ylim = c(-3.2, 2.5)) +
        ggplot2::facet_grid(rows = ggplot2::vars(region_name_plot),
                            cols = ggplot2::vars(IndexPanel), scales = "fixed",
                            drop = FALSE) +
        ggplot2::labs(
          title = paste0(RECENT_LABEL_14, " snow drought evolution across Nechako subbasins"),
          subtitle = paste0(SNOW_LAB_SWEI, " and ", SNOW_LAB_SSPI,
                            "; rows = seven supplied subbasins. ",
                            sprintf("Dashed line: x₀ = %.2f; dotted line: −1.0.", DROUGHT_ONSET)),
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
      
      save_pdf_png_14(fig_snow_subbasins, REGIONAL_PLOT_DIR,
                      "Fig2S_Recent_Subbasins_SnowDrought_7x2",
                      width = 10, height = 14)
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
        title = sprintf("%s — %s-%d", region_name, regional_index_labels[[idx]], as.integer(sc)),
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
  
  # Comparative figure for each scale: regional facets × available indices, including soil moisture.
  for (sc in intersect(regional_scales, unique(regional_all$scale))) {
    d <- regional_all[scale == sc & index_type %in% regional_index_types]
    if (!nrow(d)) next
    comparative_index_labels <- regional_index_labels
    unknown_indices <- setdiff(unique(d$index_type), names(comparative_index_labels))
    if (length(unknown_indices)) stop("Unlabelled regional indices: ", paste(unknown_indices, collapse = ", "))
    d[, Index := factor(comparative_index_labels[index_type], levels = unname(comparative_index_labels))]
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
  cat("⚠ No regional index files found; run 13trend_test_ALL.R first.\n")
}

################################################################################
# PART 6 — OTHER HYDROLOGICAL INDICES WITH DIFFERENT SPATIAL SUPPORT
#
# SSI is station-specific, TSDI is derived from basin-wide lake storage,
# NTSDI represents the Nechako Reservoir, SGI is a median over observation
# wells and MSDI is a vine-copula composite of four components (Supply, Snow, Soil, Flow). These
# are deliberately exported to a separate table instead of being labelled as
# area-weighted basin means.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 6: SSI, TSDI, NTSDI, SWSI, SGI and MSDI support-aware time-series exports\n")
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
       site = "Nechako Reservoir", source = "7_WBInex_NTSDI_GeoLakes_LocalRS.r"),
  # Dual-track SWSI (Pipeline 3 of Script 7), written to swsi_results/
  list(dir = file.path(WD_PATH, "swsi_results", "hswsi"), file = "hswsi_timeseries.csv", col = "hswsi",
       index = "H-SWSI", support = "basin composite (gauge + lakes + ERA5-Land basin means)",
       site = "Nechako Basin", source = "7_WBInex_NTSDI_GeoLakes_LocalRS.r (Pipeline 3)"),
  list(dir = file.path(WD_PATH, "swsi_results", "oswsi"), file = "oswsi_timeseries.csv", col = "oswsi",
       index = "O-SWSI", support = "managed-system composite (reservoir + SWE + inflow)",
       site = "Nechako Reservoir system", source = "7_WBInex_NTSDI_GeoLakes_LocalRS.r (Pipeline 3)"),
  # Groundwater SGI (9SGI.R): median of per-well SGI over qualifying wells
  # inside the 7-subbasin union; the contributing well count varies in time.
  list(dir = file.path(WD_PATH, "sgi_results"), file = "sgi_nechako_multisource_basin.csv",
       col = "sgi_median", index = "SGI", support = "groundwater observation-well network median",
       site = "Nechako Basin wells", source = "9SGI.R"),
  # Multivariate MSDI (10MSDI.r v5): regime-based vine-copula composite
  # of Supply, Snow, Soil, and Flow.
  list(dir = file.path(WD_PATH, "msdi_results_N5_regime"), file = "Nechako_MSDI_VineCopula.csv",
       col = "MSDI_seasonal", index = "MSDI",
       support = "regime-based vine-copula composite (Supply, Snow, Soil, Flow)",
       site = "Nechako Basin", source = "10MSDI.r (v5)")
)
for (spec in ntsdi_specs) {
  f <- file.path(if (!is.null(spec$dir)) spec$dir else ntsdi_dir, spec$file)
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
  cat("  ⚠ No TSDI/NTSDI/SWSI/SGI/MSDI or station SSI source files were available.\n")
}

################################################################################
# PART 7 — SNOW-DROUGHT COUNTERPARTS FOR SPI/SPEI-ONLY FIGURE FAMILIES
#
# SPI/SPEI manuscript figures use six accumulation scales. Snow products in the
# current framework are SWEI-3 and SSPI-1, so they receive dedicated, equivalent
# chronology figures rather than being forced into incompatible six-scale grids.
# This section writes:
#   Fig2_MS_SnowDrought_Historical_2x1.pdf/.png       (full record)
#   Fig2S_MS_Recent_SnowDrought_2x1.pdf/.png          (the configured recent period)
#   FigS_SnowDrought_Subbasins_Historical_7x2.pdf/.png (seven subbasins, full record)
# The recent 7×2 subbasin figure is created in the regional section above.
################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 7: Snow-drought counterparts for historical/recent figures\n")
cat("══════════════════════════════════════════════════════════════\n")

snow_basin_series <- stats::setNames(list(
  tryCatch(load_any_ts("swei", SWEI_SCALE), error = function(e) NULL),
  tryCatch(load_any_ts("sspi", SSPI_SCALE), error = function(e) NULL)
), c(SNOW_LAB_SWEI, SNOW_LAB_SSPI))
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
                                       levels = c(SNOW_LAB_SWEI, SNOW_LAB_SSPI))
  
  if (nrow(snow_basin_long)) {
    save_snow_basin_chronology <- function(d, file_stub, title_text, subtitle_text) {
      if (!nrow(d)) {
        cat(sprintf("  ⚠ No data for %s; figure skipped.\n", file_stub))
        return(invisible(NULL))
      }
      finite_y <- d$value[is.finite(d$value)]
      ylo <- min(-3.5, floor(min(finite_y, na.rm = TRUE) * 1.05))
      yhi <- max( 3.5, ceiling(max(finite_y, na.rm = TRUE) * 1.05))
      cols <- stats::setNames(c("#1E90FF", "#9467bd"), c(SNOW_LAB_SWEI, SNOW_LAB_SSPI))
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
      paste0(SNOW_LAB_SWEI, " and ", SNOW_LAB_SSPI, "; dashed line: drought onset threshold; dotted line: −1")
    )
    snow_recent_basin <- snow_basin_long[
      snow_basin_long$date >= RECENT_START_4e & snow_basin_long$date <= RECENT_END_4e, , drop = FALSE
    ]
    save_snow_basin_chronology(
      snow_recent_basin,
      "Fig2S_MS_Recent_SnowDrought_2x1",
      paste0("Recent snow-drought chronology for the Nechako Basin (", RECENT_LABEL_14, ")"),
      paste0(SNOW_LAB_SWEI, " and ", SNOW_LAB_SSPI, "; dashed line: drought onset threshold; dotted line: −1")
    )
  } else {
    cat("  ⚠ SWEI/SSPI inputs had no finite dated values; basin snow chronology figures skipped.\n")
  }
} else {
  cat("  ⚠ Neither SWEI-3 nor SSPI-1 basin series is available; snow chronology figures skipped.\n")
}

# Full-record subbasin counterpart to the existing the configured recent period 7×2 snow figure.
# Uses Script 13 regional products only; no basin or regional means are recomputed.
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
      selected <- select_subbasins_14(snow_hist_regional, "Historical snow subbasin figure")
      snow_hist_regional <- selected$data
      hist_region_ids <- selected$ids
      
      snow_hist_regional[, region_id_chr := as.character(region_id)]
      snow_hist_regional[, region_id := factor(region_id_chr, levels = hist_region_ids)]
      hist_region_names <- unique(snow_hist_regional[, .(region_id_chr, region_name)])
      hist_region_names <- hist_region_names[match(hist_region_ids, hist_region_names$region_id_chr)]
      snow_hist_regional[, region_name_plot := factor(region_name,
                                                      levels = hist_region_names$region_name)]
      snow_hist_regional[, IndexPanel := ifelse(tolower(index_type) == "swei", SNOW_LAB_SWEI, SNOW_LAB_SSPI)]
      snow_hist_regional[, IndexPanel := factor(IndexPanel, levels = c(SNOW_LAB_SWEI, SNOW_LAB_SSPI))]
      snow_hist_colours <- stats::setNames(c("#1E90FF", "#9467bd"), c(SNOW_LAB_SWEI, SNOW_LAB_SSPI))
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
          subtitle = paste0(SNOW_LAB_SWEI, " and ", SNOW_LAB_SSPI, "; dashed line: drought onset threshold; dotted line: −1"),
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

################################################################################
# FINAL SUMMARY (printed after every part has run)
################################################################################
cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  14basin_timeseries.R  DONE                 ║\n")
cat("╚══════════════════════════════════════════════╝\n\nOutputs in ",
    normalizePath(BASIN_PLOT_DIR, winslash = "/", mustWork = FALSE), ":\n", sep = "")
for (f in c("basin_mean_all_indices.csv",
            "Fig4_MS_SPEI_136_PM_timeseries.pdf", "Fig2_MS_AllIndices_2x6.pdf",
            "Fig2S_MS_Recent_AllScales_3x2.pdf",
            sprintf("Fig4c_%s_SPEI_PM_SPI.pdf", CORR_TAG), sprintf("Fig4c_%s_SPEI136_PM.pdf", CORR_TAG),
            sprintf("Fig4c_%s_SWEI_SSPI.pdf", CORR_TAG), "FigS_Seasonality_P_PET.pdf",
            "Drought_Summary.xlsx"))
{
  output_path <- file.path(BASIN_PLOT_DIR, f)
  info <- file.info(output_path)
  this_run <- nrow(info) == 1L && !is.na(info$mtime) &&
    info$mtime >= SCRIPT14_STARTED && !is.na(info$size) && info$size > 0
  cat(sprintf("  %-45s %s\n", f, if (this_run) "written this run" else "NOT written this run (missing, stale or empty)"))
}
cat("  Regional series : ", normalizePath(REGIONAL_TS_DIR,   winslash = "/", mustWork = FALSE), "\n")
cat("  Regional plots  : ", normalizePath(REGIONAL_PLOT_DIR, winslash = "/", mustWork = FALSE), "\n")
cat("  Other indices   : ", normalizePath(OTHER_TS_DIR,      winslash = "/", mustWork = FALSE), "\n")
# Audit Fig4c provenance: this run must have produced BOTH formats for all
# three delegated manuscript figures. Failures are explicit in this manifest.
fig4c_files_14 <- sprintf("Fig4c_%s_%s.%s", CORR_TAG,
                          rep(c("SPEI_PM_SPI", "SPEI136_PM", "SWEI_SSPI"), each = 2L),
                          rep(c("pdf", "png"), times = 3L))
fig4c_paths_14 <- file.path(BASIN_PLOT_DIR, fig4c_files_14)
fig4c_info_14 <- file.info(fig4c_paths_14)
fig4c_manifest_14 <- data.frame(
  filename = fig4c_files_14,
  written_this_run = !is.na(fig4c_info_14$mtime) & fig4c_info_14$mtime >= SCRIPT14_STARTED &
    !is.na(fig4c_info_14$size) & fig4c_info_14$size > 0,
  mtime_utc = format(fig4c_info_14$mtime, tz = "UTC", usetz = TRUE),
  size_bytes = fig4c_info_14$size,
  script14_start_utc = format(SCRIPT14_STARTED, tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE)
utils::write.csv(fig4c_manifest_14, file.path(BASIN_PLOT_DIR, "Fig4c_export_manifest.csv"), row.names = FALSE)
if (!all(fig4c_manifest_14$written_this_run))
  warning("Script 14 Fig4c exports incomplete: ",
          paste(fig4c_manifest_14$filename[!fig4c_manifest_14$written_this_run], collapse = ", "))
cat("\nBasin-mean indices and non-basin-support indices are kept in separate outputs.\n")
cat("Next: run 15trends_visualization.R.\n\n")