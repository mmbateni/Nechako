####################################################################################
# 9pr_pet_trends_visualization.r  (FIXED VERSION)
# Trend maps + temporal change + dual-PET plots + seasonal plots
#
# Reads results produced by 9pr_pet_trends.r (observed PM + ERA5-Land PEV framework).
#
# SUMMARY OF FIXES vs. the original
#   1. subbasin_wgs is now built ONCE at top level and used by BOTH map figures
#      (previously created inside Function 0b but used in Function 0 -> fatal error).
#   2. Every generation step runs inside .run_step(): one failure no longer stops the
#      run or drops you into the debugger (Browse[1]>).
#   3. VC-MK and TFPW statistics are repaired upstream in 9pr_pet_trends.r using
#      the documented modifiedmk implementations; visualization does not alter
#      statistical signs or significance.
#   4. melt() name collision fixed; agreement percentages are now within variable.
#   5. Output files locked by a viewer (Windows) no longer abort a figure: the script
#      tries to overwrite, otherwise writes a time-stamped copy and says so.
#   6. Script is pure ASCII (Unicode via \u escapes); PDFs use cairo_pdf when
#      available so tau / degree / minus characters render correctly.
#   7. Publication maps use symmetric, data-driven colour limits (old fixed
#      c(0, 1.5) clipped negative PET slopes).
#   8. Seasonal plots (Functions 15/16/17) share one helper and drop incomplete
#      seasons (the first DJF had only 2 months but was multiplied by 3).
#   9. Diagnostics retained: tau_tfpw vs Sen-slope concordance, pixel/grid alignment,
#      PET-PM vs PET-ERA5 magnitude check. These are QA checks only.
#  10. Duplicate header / setup blocks removed; Functions 16 and 17 integrated.
#  11. SWE (ERA5-Land snow water equivalent; variables SWE and SWE_Peak from 9pr_pet_trends.r v7)
#      is now treated like the other variables everywhere: publication maps, tile Sen-slope maps
#      (annual + by calendar month), tau comparison, regime-shift / spectral / spatial-pattern maps,
#      method comparison, basin time series, calendar-month plots and slope profile, seasonal
#      trends (Function 18), specific-point plots, and regional plots.
#  12. SWE step-change outputs (Function 19): figures (PNG) and tables (CSV) written to
#      swe_homogeneity/report/ for use in a later report: tests of whether the post-1978 SWE drop is a
#      step, whether Nov-Mar temperature has a step, and whether any SWE step survives removal of
#      precipitation and temperature. tables/Table00_key_numbers.csv lists each number with its source.
#      No narrative text is written. Reads the CSVs written by 9pr_pet_trends.r (v7); needs no pixel data.
#  13. Evidence outputs (Functions 21-23, swe_homogeneity/evidence and swe_homogeneity/link_to_summer):
#      Claim 1 - Nov-Mar temperature change is better described as a trend than a step (AR(1) GLS model
#      comparison, parametric bootstrap, monthly contributions); Claim 2 - precipitation matters about as
#      much as temperature for the SWE change (bootstrap attribution by period, variance shares, decades);
#      Function 23 - does a low-snow spring precede low SPEI (PM and ERA5) in Jul-Sep, from the existing basin
#      CSVs, with overlap flags, detrending and multiple-testing control. Evidence_manifest.csv lists everything.
#  14. Runoff, Snowfall, SnowDepth and SnowDensity (9pr_pet_trends.r v8) are added everywhere: publication
#      maps, a new annual Sen-slope figure (Function 0e), monthly Sen-slope maps (Function 0f),
#      regime-shift / spectral / spatial-pattern maps (now laid out on a grid), method comparison, basin
#      time series, calendar-month plots and slope profile, seasonal trends (Function 18b), point plots
#      and regional plots.
#  15. Downstream of 9pr_pet_trends.r v9-v10: AET (now trend-tested upstream) is handled like the
#      runoff/snow variables everywhere (publication maps, Functions 0e/0f, method comparison, basin and
#      regional series, calendar-month plots, seasonal trends, point plots). Functions 0e, 0f, 18b and the
#      runoff/snow/AET publication maps are now actually called. PM controls get regional figures.
#      Function 24 plots the multi-variable homogeneity outputs (homogeneity/): basin break years with a
#      step-vs-trend classification, the classification by region, pixel break shares and pixel maps.
#      Diagnostic C reports the known PEV magnitude (and the diagnose_pet_era5.R energy check) instead of
#      warning; ERA5-Land PEV gets its own basin time-series panel. Upstream input QA is printed.
#      Function 7 writes summary_statistics_visualization.csv instead of overwriting summary_statistics.csv.
####################################################################################

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

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(zoo)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(data.table)   # attached after terra so data.table verbs win any masking
})

# ---- user switches ---------------------------------------------------------------
# Minimum number of pixels valid for BOTH VC and TFPW to draw the method comparison.
MIN_PIXELS_METHOD_COMPARISON <- 30L

# Shared project constants and basin utilities (load_basin(), load_nechako_regions(), ...)
source(UTILS_FILE, encoding = "UTF-8")

out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

####################################################################################
# DEVICE / OUTPUT HELPERS
####################################################################################
.HAS_CAIRO <- isTRUE(unname(capabilities("cairo")))
.PDF_DEV   <- if (.HAS_CAIRO) grDevices::cairo_pdf else "pdf"
TAU        <- if (.HAS_CAIRO) "\u03c4" else "tau"     # Greek tau only when device can draw it
DEGC       <- "\u00b0C"

# If the target file is locked (open in a viewer) fall back to a time-stamped name.
.resolve_out <- function(path) {
  if (!file.exists(path)) return(path)
  removed <- suppressWarnings(tryCatch(file.remove(path), error = function(e) FALSE))
  if (isTRUE(removed)) return(path)
  alt <- sub("(\\.[^.]+)$",
             paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S"), "\\1"), path)
  cat(sprintf("  [WARN] %s is locked (open in a viewer?) -> writing %s instead\n",
              basename(path), basename(alt)))
  alt
}

# Save a ggplot / patchwork object as PNG and/or PDF, verifying the file exists.
.save_gg <- function(stem, plot, width, height, dpi = 300,
                     formats = c("png", "pdf")) {
  out <- character(0)
  for (fmt in formats) {
    path <- .resolve_out(paste0(stem, ".", fmt))
    ok <- tryCatch({
      withCallingHandlers({
        if (fmt == "pdf") {
          ggplot2::ggsave(path, plot, width = width, height = height, units = "in",
                          device = .PDF_DEV)
        } else {
          ggplot2::ggsave(path, plot, width = width, height = height, units = "in",
                          dpi = dpi)
        }
      }, warning = function(w) {
        cat(sprintf("  [WARN] %s: %s\n", basename(path), conditionMessage(w)))
        invokeRestart("muffleWarning")
      })
      # a device warning ("could not write") leaves no file, so verify on disk
      file.exists(path) && isTRUE(file.info(path)$size > 0)
    }, error = function(e) {
      cat(sprintf("  [WARN] %s not saved: %s\n", basename(path), conditionMessage(e)))
      FALSE
    })
    if (isTRUE(ok)) out <- c(out, basename(path))
  }
  invisible(out)
}

.pdf_safe <- function(pdf_path, width, height, expr_fn) {
  pdf_path <- .resolve_out(pdf_path)
  opened <- FALSE
  tryCatch({
    if (.HAS_CAIRO) grDevices::cairo_pdf(pdf_path, width = width, height = height)
    else            grDevices::pdf(pdf_path, width = width, height = height)
    opened <- TRUE
    expr_fn()
    dev.off()
    opened <- FALSE
    invisible(pdf_path)
  }, error = function(e) {
    if (opened && dev.cur() > 1) dev.off()
    cat(sprintf("  [WARN] PDF failed (%s): %s\n", basename(pdf_path), conditionMessage(e)))
    invisible(NULL)
  })
}

# Run one output step; a failure is reported but never stops the script.
.run_step <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    cat(sprintf("  [FAILED] %s: %s\n", label, conditionMessage(e)))
    invisible(NULL)
  })
}

####################################################################################
# LOAD DATA
####################################################################################
cat("Loading results, basin boundary, and metadata...\n")
all_results    <- readRDS(file.path(out_dir, "all_results.rds"))
basin_boundary <- readRDS(file.path(out_dir, "basin_boundary.rds"))
metadata       <- readRDS(file.path(out_dir, "analysis_metadata.rds"))

data.table::setDT(all_results)
data.table::setnames(all_results, trimws(names(all_results)))
all_results[, variable := trimws(variable)]
all_results[, period   := trimws(period)]

names(metadata$basin_avg_monthly) <- trimws(names(metadata$basin_avg_monthly))
names(metadata$basin_avg_annual)  <- trimws(names(metadata$basin_avg_annual))
basin_avg_monthly <- metadata$basin_avg_monthly
basin_avg_annual  <- metadata$basin_avg_annual
basin_avg_monthly$date <- as.Date(basin_avg_monthly$date)
SPECIFIC_PTS      <- metadata$specific_pts

pt_trend_file <- file.path(out_dir, "point_trend_stats.csv")
pt_ts_file    <- file.path(out_dir, "point_monthly_timeseries.csv")
if (!file.exists(pt_trend_file))
  warning("point_trend_stats.csv not found - skipping specific-point plots")
if (!file.exists(pt_ts_file))
  warning("point_monthly_timeseries.csv not found - skipping specific-point TS plots")
point_trend_stats <- if (file.exists(pt_trend_file)) fread(pt_trend_file) else NULL
point_monthly_ts  <- if (file.exists(pt_ts_file))    fread(pt_ts_file)    else NULL

# ---- authoritative raster template = native ERA5-Land precipitation grid ----------
native_precip_file <- file.path(WD_PATH, "monthly_data_direct",
                                "total_precipitation_monthly.nc")
if (!file.exists(native_precip_file))
  stop("Native ERA5-Land precipitation NetCDF not found: ", native_precip_file)

native_precip <- terra::rast(native_precip_file)
if (!inherits(native_precip, "SpatRaster") || terra::nlyr(native_precip) < 1L)
  stop("Failed to load native ERA5-Land precipitation NetCDF as a SpatRaster.")

analysis_template <- native_precip[[1]]
if (is.na(terra::crs(analysis_template)) || !nzchar(terra::crs(analysis_template)))
  stop("Native ERA5-Land precipitation grid has no CRS.")
template_bc <- analysis_template   # name used throughout the functions below

cat(sprintf("Native ERA5-Land grid: %d x %d, res = %.4f x %.4f deg\n",
            terra::nrow(analysis_template), terra::ncol(analysis_template),
            terra::res(analysis_template)[1], terra::res(analysis_template)[2]))

# Optional QA of the cached template (never used for geometry)
clipped_template_file <- file.path(out_dir, "clipped_template.rds")
if (file.exists(clipped_template_file)) {
  clipped_template <- tryCatch(readRDS(clipped_template_file), error = function(e) NULL)
  if (inherits(clipped_template, "SpatRaster")) {
    same <- tryCatch(
      isTRUE(terra::compareGeom(clipped_template, analysis_template, stopOnError = FALSE,
                                crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)),
      error = function(e) FALSE)
    if (same) cat("clipped_template.rds matches the native grid.\n")
    else cat("WARNING: clipped_template.rds does NOT match the native grid; ignored.\n")
  }
}

basin_boundary <- sf::st_transform(basin_boundary,
                                   sf::st_crs(terra::crs(analysis_template)))
basin_wgs <- tryCatch(sf::st_transform(sf::st_as_sf(basin_boundary), 4326),
                      error = function(e) sf::st_as_sf(basin_boundary))

####################################################################################
# SUB-BASIN BOUNDARIES (built once, used by all map figures)  [FIX #1]
# Only the geometry column is kept, one dissolved feature per sub-basin, so regions
# with different attribute columns can never break the combine step.
####################################################################################
load_subbasins_wgs <- function() {
  if (!exists("load_nechako_regions", mode = "function")) return(NULL)
  regs <- tryCatch(load_nechako_regions(), error = function(e) {
    cat("  [WARN] load_nechako_regions() failed: ", conditionMessage(e), "\n", sep = "")
    NULL
  })
  if (!is.list(regs) || !length(regs)) return(NULL)
  whole_id <- if (exists("NECHAKO_REGION_BASIN_ID", inherits = TRUE))
    get("NECHAKO_REGION_BASIN_ID", inherits = TRUE) else "Nechako_Basin"
  regs <- regs[setdiff(names(regs), whole_id)]
  regs <- regs[vapply(regs, inherits, logical(1), what = "SpatVector")]
  if (!length(regs)) return(NULL)
  geoms <- lapply(regs, function(rv) {
    x <- sf::st_make_valid(sf::st_as_sf(rv))
    x <- sf::st_transform(x, 4326)
    sf::st_geometry(sf::st_union(x))          # one feature per region
  })
  sf::st_sf(region = names(regs), geometry = do.call(c, geoms))
}
subbasin_wgs <- tryCatch(load_subbasins_wgs(), error = function(e) {
  cat("  [WARN] Sub-basin overlay disabled: ", conditionMessage(e), "\n", sep = "")
  NULL
})
cat(sprintf("Sub-basin overlay: %s\n",
            if (is.null(subbasin_wgs)) "not available" else
              sprintf("%d sub-basins", nrow(subbasin_wgs))))

####################################################################################
# SANITY CHECKS
####################################################################################
.vars_present <- unique(all_results$variable)
cat("\nVariables in all_results:", paste(.vars_present, collapse = ", "), "\n")
has_pet_era5 <- "PET_ERA5" %in% .vars_present
if (!has_pet_era5)
  cat("\nWARNING: PET_ERA5 NOT in all_results - dual-PET comparison plots will be skipped.\n")
SWE_VARS_PRESENT <- intersect(c("SWE", "SWE_Peak"), .vars_present)
cat("SWE variables present: ",
    if (length(SWE_VARS_PRESENT)) paste(SWE_VARS_PRESENT, collapse = ", ") else "NONE (re-run 9pr_pet_trends.r v7 for SWE outputs)",
    "\n", sep = "")

PM_CONTROL_VARS <- c("PM_Temperature", "PM_NetRadiation", "PM_VPD",
                     "PM_Wind2m", "PM_SurfacePressure")
PM_CONTROL_LABELS <- c(PM_Temperature = "2-m air temperature",
                       PM_NetRadiation = "Net radiation",
                       PM_VPD = "Vapour pressure deficit",
                       PM_Wind2m = "2-m wind speed",
                       PM_SurfacePressure = "Surface pressure")
PM_CONTROL_UNITS <- c(PM_Temperature = "degC yr-1",
                      PM_NetRadiation = "MJ m-2 d-1 yr-1",
                      PM_VPD = "kPa yr-1",
                      PM_Wind2m = "m s-1 yr-1",
                      PM_SurfacePressure = "kPa yr-1")
PM_CONTROL_PRESENT <- intersect(PM_CONTROL_VARS, .vars_present)
PM_CONTROL_MISSING <- setdiff(PM_CONTROL_VARS, .vars_present)
cat("PM-control variables present: ",
    if (length(PM_CONTROL_PRESENT)) paste(PM_CONTROL_PRESENT, collapse = ", ") else "NONE",
    "\n", sep = "")
if (length(PM_CONTROL_MISSING))
  cat("WARNING: PM-control variables missing: ",
      paste(PM_CONTROL_MISSING, collapse = ", "), "\n", sep = "")

# ---- runoff / snowfall / snow depth / snow density (9pr_pet_trends.r v8 EXTRA_VARS) ----
EXTRA_VAR_LABELS  <- c(AET = "Actual evapotranspiration", Runoff = "Runoff", Snowfall = "Snowfall",
                       SnowDepth = "Snow depth (Mar-Apr mean)", SnowDensity = "Snow density (Mar-Apr mean)")
EXTRA_SLOPE_UNITS <- c(AET = "mm yr-1", Runoff = "mm yr-1", Snowfall = "mm yr-1", SnowDepth = "cm yr-1",
                       SnowDensity = "kg m-3 yr-1")
EXTRA_MONTHLY_SLOPE_UNITS <- c(AET = "mm month-1 yr-1", Runoff = "mm month-1 yr-1", Snowfall = "mm month-1 yr-1",
                               SnowDepth = "cm yr-1", SnowDensity = "kg m-3 yr-1")
EXTRA_MONTH_COLS  <- c(AET = "aet_mm_month", Runoff = "runoff_mm_month", Snowfall = "snowfall_mm_month",
                       SnowDepth = "snow_depth_cm_month", SnowDensity = "snow_density_kgm3_month")
EXTRA_MONTH_UNITS <- c(AET = "mm/month", Runoff = "mm/month", Snowfall = "mm w.e./month", SnowDepth = "cm",
                       SnowDensity = "kg m-3")
EXTRA_COLOURS     <- c(AET = "#66a61e", Runoff = "#2c7fb8", Snowfall = "#7570b3", SnowDepth = "#1b9e77",
                       SnowDensity = "#a6761d")
.short_lab <- function(v) sub(" \\(.*", "", EXTRA_VAR_LABELS[[v]])
EXTRA_VARS_PRESENT <- intersect(names(EXTRA_VAR_LABELS), .vars_present)
cat("Runoff/snow variables present: ",
    if (length(EXTRA_VARS_PRESENT)) paste(EXTRA_VARS_PRESENT, collapse = ", ") else "NONE (re-run 9pr_pet_trends.r v8)",
    "\n", sep = "")
MAP_VARS <- intersect(c("Precipitation", "PET", "PET_ERA5", "Temperature", "SWE", "SWE_Peak",
                        EXTRA_VARS_PRESENT), .vars_present)
.grid_dims <- function(n, max_cols = 5L) c(nr = ceiling(n / max_cols), nc = min(n, max_cols))

n_coords_ok <- sum(!all_results$is_basin_average &
                     !is.na(all_results$x) & !is.na(all_results$y))
cat(sprintf("Pixel rows with valid x/y: %d\n", n_coords_ok))
if (n_coords_ok == 0) stop("All pixel coordinates are NA - re-run the processing script.")

required_cols <- c("x", "y", "variable", "period", "tau_vc", "p_value_vc", "sl_vc",
                   "tau_tfpw", "p_value_tfpw", "sl_tfpw", "filtered_vc",
                   "filtered_tfpw", "is_basin_average")
missing_cols <- setdiff(required_cols, names(all_results))
if (length(missing_cols))
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))

# ---- DIAGNOSTIC A: tau_tfpw sign vs Sen slope sign --------------------------------
.check_tfpw_sign <- function() {
  chk <- all_results[period == "annual" & !is_basin_average &
                       is.finite(tau_tfpw) & is.finite(sl_tfpw) & sl_tfpw != 0,
                     .(n = .N,
                       agree_all      = mean(sign(tau_tfpw) == sign(sl_tfpw)),
                       n_sig          = sum(p_value_tfpw < 0.05, na.rm = TRUE),
                       agree_sig      = { s <- which(p_value_tfpw < 0.05)
                       if (length(s)) mean(sign(tau_tfpw[s]) == sign(sl_tfpw[s])) else NA_real_ },
                       agree_vc_slope = mean(sign(tau_vc) == sign(sl_vc), na.rm = TRUE)),
                     by = variable]
  cat("\nDIAGNOSTIC: sign(tau_tfpw) == sign(Sen slope), all pixels and TFPW-significant pixels\n")
  print(as.data.frame(chk), digits = 3)
  bad <- chk[is.finite(agree_sig) & agree_sig < 0.99]
  if (nrow(bad)) {
    warning("TFPW tau and Sen slope disagree in sign among SIGNIFICANT pixels for: ",
            paste(bad$variable, collapse = ", "), ". Inspect the TFPW implementation.", call. = FALSE)
  } else if (any(chk$agree_all < 0.99)) {
    cat("  Sign disagreements occur only in non-significant (near-zero-trend) pixels: expected for TFPW,\n",
        "  where tau comes from the prewhitened series and the slope from the original series.\n", sep = "")
  }
  invisible(chk)
}
.run_step("tfpw sign diagnostic", .check_tfpw_sign())

# ---- DIAGNOSTIC B: pixel coordinates vs native grid cell centres -----------------
.check_grid_alignment <- function() {
  px <- unique(all_results[!is_basin_average & is.finite(x) & is.finite(y), .(x, y)])
  cells <- terra::cellFromXY(template_bc, as.matrix(px[, .(x, y)]))
  n_out <- sum(is.na(cells))
  keep  <- !is.na(cells)
  ctr   <- terra::xyFromCell(template_bc, cells[keep])
  rs    <- terra::res(template_bc)
  off   <- mean(abs(px$x[keep] - ctr[, 1]) / rs[1] > 0.05 |
                  abs(px$y[keep] - ctr[, 2]) / rs[2] > 0.05)
  dup   <- sum(duplicated(cells[keep]))
  cat(sprintf(paste0("\nDIAGNOSTIC: %d unique pixel coordinates | outside grid: %d | ",
                     "off cell-centre: %.1f%% | sharing a cell: %d\n"),
              nrow(px), n_out, 100 * off, dup))
  if (n_out > 0 || off > 0.01 || dup > 0)
    warning("Pixel coordinates in all_results do not line up with the native ",
            "ERA5-Land grid; maps may be shifted or overwritten.", call. = FALSE)
  invisible(NULL)
}
.run_step("grid alignment diagnostic", .check_grid_alignment())

# ---- DIAGNOSTIC C: PET-PM vs PET-ERA5 magnitude / sign ---------------------------
.check_pet_units <- function() {
  if (!all(c("pet_mm_month", "pet_era5_mm_month") %in% names(basin_avg_monthly)))
    return(invisible(NULL))
  m_pm   <- mean(basin_avg_monthly$pet_mm_month,      na.rm = TRUE)
  m_era5 <- mean(basin_avg_monthly$pet_era5_mm_month, na.rm = TRUE)
  cat(sprintf("\nDIAGNOSTIC: mean monthly PET  PM = %.1f mm | ERA5 PEV = %.1f mm | ratio = %.2f\n",
              m_pm, m_era5, m_era5 / m_pm))
  qa <- metadata$pet_era5_qa
  if (!is.null(qa)) cat(sprintf("  Upstream QA (pet_source_qa): status %s - %s\n", qa$status, qa$interpretation))
  d2f <- file.path(out_dir, "pet_era5_diagnostics", "D2_monthly_climatology.csv")
  if (file.exists(d2f)) {
    j <- data.table::fread(d2f)[month == 7]
    if (nrow(j)) cat(sprintf("  diagnose_pet_era5.R (July): PEV / radiation limit = %.2f | PEV / Penman open water = %.2f%s\n",
                             j$pev_over_rn, j$pev_over_e0,
                             if ("pev_over_ssrd" %in% names(j)) sprintf(" | PEV / all incoming shortwave = %.2f", j$pev_over_ssrd) else ""))
  } else cat("  Run diagnose_pet_era5.R for the energy-bound check of ERA5-Land PEV.\n")
  if (!is.finite(m_era5 / m_pm) || m_era5 < 0)
    warning("PET_ERA5 basin mean is missing, zero or negative - check the upstream PEV file.", call. = FALSE)
  invisible(NULL)
}
.run_step("PET unit diagnostic", .check_pet_units())

# ---- Upstream input QA (9pr_pet_trends.r v9-v10) ----------------------------------
.nz1 <- function(x) if (is.null(x) || !length(x)) NA else x[[1]]
.print_input_qa <- function() {
  nq <- metadata$negative_clamp_qa
  if (!is.null(nq) && nrow(nq)) {
    cat("\nINPUT QA: negative values set to 0 upstream\n")
    tq <- data.table::as.data.table(nq)[, .(variable, unit, n_negative, pct_negative, most_negative,
                                            n_beyond_tolerance, pct_negative_in_JunOct)]
    tq[, most_negative := formatC(most_negative, format = "e", digits = 2)]
    print(as.data.frame(tq), digits = 3)  }
  q <- metadata$swe_input_qa
  if (!is.null(q))
    cat(sprintf("INPUT QA: SWE units '%s' | March median %.0f mm | r(snow depth) = %.3f | implied density %.0f kg m-3 | peak SWE / Oct-Apr P = %.2f\n",
                as.character(.nz1(q$units_used)), as.numeric(.nz1(q$march_median_mm)),
                as.numeric(.nz1(q$depth_swe_march_spatial_r)), as.numeric(.nz1(q$implied_density_march_median)),
                as.numeric(.nz1(q$peak_swe_over_OctApr_precip))))
  invisible(NULL)
}
.run_step("input QA summary", .print_input_qa())

####################################################################################
# RASTER HELPERS
####################################################################################
create_raster_from_table <- function(results_dt, template, value_col) {
  value_col <- trimws(value_col)
  if (!value_col %in% names(results_dt))
    stop(sprintf("Column '%s' not found.", value_col))
  results_dt <- results_dt[!is.na(x) & !is.na(y)]
  if (!nrow(results_dt)) { r <- template; terra::values(r) <- NA; return(r) }
  pts <- terra::vect(as.data.frame(results_dt), geom = c("x", "y"),
                     crs = terra::crs(template))
  r <- terra::rasterize(pts, template, field = value_col, touches = TRUE)
  terra::mask(r, terra::vect(basin_boundary))
}

smooth_raster <- function(r) {
  if (all(is.na(terra::values(r)))) return(r)
  terra::focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
}

# Symmetric, data-driven limits for diverging palettes
.sym_lim <- function(r, fallback = 1) {
  if (is.null(r)) return(c(-fallback, fallback))
  v <- terra::values(r, mat = FALSE)
  v <- v[is.finite(v)]
  z <- if (length(v)) max(abs(v)) else fallback
  if (!is.finite(z) || z <= 0) z <- fallback
  c(-z, z)
}
.pos_lim <- function(r, fallback = 1) {
  if (is.null(r)) return(c(0, fallback))
  v <- terra::values(r, mat = FALSE)
  v <- v[is.finite(v)]
  z <- if (length(v)) max(v) else fallback
  if (!is.finite(z) || z <= 0) z <- fallback
  c(0, z)
}

plot_raster_panel <- function(r, main, zlim = NULL, col, breaks = NULL,
                              legend = TRUE, categorical = FALSE, legend_title = "",
                              smooth = FALSE) {
  if (is.null(r) || all(is.na(terra::values(r)))) {
    plot.new()
    text(0.5, 0.5, "No valid trend data\n(pixels filtered by quality criteria)",
         cex = 1.1, col = "gray40", font = 3)
    return(invisible(NULL))
  }
  vals <- terra::values(r)
  valid_pct <- sum(!is.na(vals)) / length(vals) * 100
  if (isTRUE(smooth) && !categorical) r <- smooth_raster(r)
  plg_args <- list(cex = 0.9)
  if (nzchar(legend_title)) plg_args$title <- legend_title
  if (!is.null(breaks))
    plot(r, main = main, cex.main = 1.0, col = col, breaks = breaks,
         axes = FALSE, box = FALSE, legend = legend, plg = plg_args, colNA = NA)
  else
    plot(r, main = main, cex.main = 1.0, col = col, zlim = zlim,
         axes = FALSE, box = FALSE, legend = legend, plg = plg_args, colNA = NA)
  plot(sf::st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
  if (valid_pct < 10)
    mtext("Note: most pixels filtered by data-quality criteria (ties, min-values, autocorrelation)",
          side = 1, line = 0.5, cex = 0.7, col = "gray40")
  invisible(NULL)
}

prep_annual <- function(var_name) {
  dt <- all_results[variable == var_name & period == "annual" & !is_basin_average]
  cat(sprintf("  %-20s annual pixels: %4d | VC-valid: %4d | TFPW-valid: %4d\n",
              var_name, nrow(dt),
              if (nrow(dt)) sum(!dt$filtered_vc, na.rm = TRUE) else 0L,
              if (nrow(dt)) sum(!dt$filtered_tfpw, na.rm = TRUE) else 0L))
  if (!nrow(dt)) return(NULL)
  dt <- data.table::copy(dt)
  # VC and TFPW are rendered independently: do NOT drop rows by VC validity.
  dt[, combined_vc   := fifelse(p_value_vc   < 0.05, abs(tau_vc),   0)]
  dt[, combined_tfpw := fifelse(p_value_tfpw < 0.05, abs(tau_tfpw), 0)]
  dt[, direction_vc := fifelse(p_value_vc < 0.05 & tau_vc < 0, -1L,
                               fifelse(p_value_vc < 0.05 & tau_vc > 0, 1L, 0L))]
  dt[, direction_tfpw := fifelse(p_value_tfpw < 0.05 & tau_tfpw < 0, -1L,
                                 fifelse(p_value_tfpw < 0.05 & tau_tfpw > 0, 1L, 0L))]
  dt
}

make_rasters <- function(dt) {
  if (is.null(dt)) return(NULL)
  list(
    comb_vc   = create_raster_from_table(dt, template_bc, "combined_vc"),
    comb_tfpw = create_raster_from_table(dt, template_bc, "combined_tfpw"),
    dir_vc    = create_raster_from_table(dt, template_bc, "direction_vc"),
    dir_tfpw  = create_raster_from_table(dt, template_bc, "direction_tfpw"),
    mag_vc    = create_raster_from_table(dt, template_bc, "sl_vc"),
    mag_tfpw  = create_raster_from_table(dt, template_bc, "sl_tfpw"),
    tau_vc    = create_raster_from_table(dt, template_bc, "tau_vc"),
    tau_tfpw  = create_raster_from_table(dt, template_bc, "tau_tfpw")
  )
}

cat("\nPreparing trend metrics...\n")
precip_annual   <- prep_annual("Precipitation")
pet_annual      <- prep_annual("PET")
pet_era5_annual <- if (has_pet_era5) prep_annual("PET_ERA5") else NULL
tair_annual     <- prep_annual("Temperature")
pm_control_annual <- lapply(PM_CONTROL_PRESENT, prep_annual)
names(pm_control_annual) <- PM_CONTROL_PRESENT
swe_annual      <- if ("SWE" %in% .vars_present)      prep_annual("SWE")      else NULL
swe_peak_annual <- if ("SWE_Peak" %in% .vars_present) prep_annual("SWE_Peak") else NULL

cat("Generating rasters (Pr, PET_PM, PET_ERA5, Temperature, PM controls, SWE, runoff/snow)...\n")
r_precip   <- make_rasters(precip_annual)
r_pet      <- make_rasters(pet_annual)
r_pet_era5 <- if (has_pet_era5) make_rasters(pet_era5_annual) else NULL
r_tair     <- make_rasters(tair_annual)
r_pm_control <- lapply(pm_control_annual, make_rasters)
names(r_pm_control) <- PM_CONTROL_PRESENT
r_swe      <- if (!is.null(swe_annual))      make_rasters(swe_annual)      else NULL
r_swe_peak <- if (!is.null(swe_peak_annual)) make_rasters(swe_peak_annual) else NULL
extra_annual_dt <- lapply(EXTRA_VARS_PRESENT, prep_annual); names(extra_annual_dt) <- EXTRA_VARS_PRESENT
r_extra <- lapply(extra_annual_dt, make_rasters);           names(r_extra)         <- EXTRA_VARS_PRESENT

####################################################################################
# SHARED HELPERS FOR THE TILE-STYLE (ggplot) SEN-SLOPE MAPS
####################################################################################
# Dual-test significance (framework rule for the hydroclimatic driver maps):
# BH-FDR is applied separately to the variance-corrected (Hamed-Rao VC-MK) and
# the TFPW-MK p-values within each group. Filled black points = significant
# under BOTH tests; open points = significant under only one of them. The
# colour is the TFPW Sen slope. sig_fdr is kept as an alias of sig_both.
.bh_finite <- function(p) {
  out <- rep(NA_real_, length(p)); ok <- is.finite(p)
  if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = "BH")
  out
}
.add_dual_fdr <- function(d, by = NULL) {
  if (!"filtered_vc"   %in% names(d)) d[, filtered_vc   := FALSE]
  if (!"filtered_tfpw" %in% names(d)) d[, filtered_tfpw := FALSE]
  d[, p_vc_use   := data.table::fifelse(filtered_vc   %in% TRUE, NA_real_, as.numeric(p_value_vc))]
  d[, p_tfpw_use := data.table::fifelse(filtered_tfpw %in% TRUE, NA_real_, as.numeric(p_value_tfpw))]
  if (is.null(by)) {
    d[, `:=`(p_fdr_vc_plot = .bh_finite(p_vc_use), p_fdr_tfpw_plot = .bh_finite(p_tfpw_use))]
  } else {
    d[, `:=`(p_fdr_vc_plot = .bh_finite(p_vc_use), p_fdr_tfpw_plot = .bh_finite(p_tfpw_use)), by = by]
  }
  d[, sig_fdr_vc   := is.finite(p_fdr_vc_plot)   & p_fdr_vc_plot   < 0.05]
  d[, sig_fdr_tfpw := is.finite(p_fdr_tfpw_plot) & p_fdr_tfpw_plot < 0.05]
  d[, sig_both := sig_fdr_vc & sig_fdr_tfpw]
  d[, sig_one  := xor(sig_fdr_vc, sig_fdr_tfpw)]
  d[, sig_fdr  := sig_both]
  d[, c("p_vc_use", "p_tfpw_use") := NULL]
  d
}

# Annual points with dual-test BH-FDR applied separately within each variable.
.prep_tfpw_points <- function(vars) {
  d <- data.table::copy(all_results[period == "annual" & !is_basin_average &
                                      variable %in% vars &
                                      is.finite(sl_tfpw) & is.finite(p_value_tfpw)])
  if (!nrow(d)) return(NULL)
  d <- .add_dual_fdr(d, by = "variable")
  pts <- sf::st_as_sf(as.data.frame(d), coords = c("x", "y"),
                      crs = terra::crs(template_bc), remove = FALSE)
  pts <- sf::st_transform(pts, 4326)
  xy  <- sf::st_coordinates(pts)
  pts$plot_x <- xy[, 1]
  pts$plot_y <- xy[, 2]
  pts
}

.grid_spacing <- function(v) {
  u  <- sort(unique(round(v, 8)))
  dd <- diff(u)
  dd <- dd[is.finite(dd) & dd > 0]
  if (!length(dd)) return(0.05)
  stats::median(dd)
}

# ---- colour palettes for the tile-style Sen-slope maps ---------------------------------
# Moisture variables use a brown <-> blue-green diverging ramp (ColorBrewer BrBG, 9 classes,
# near-white at zero) so that colour always encodes wetter vs drier:
#   .pal_wet_positive : positive slope = wetter -> blue-green ; negative = drier  -> brown
#   .pal_dry_positive : positive slope = drier  -> brown      ; negative = wetter -> blue-green
# Non-moisture variables keep the conventional thermal ramp (red = increase, blue = decrease).
.BRBG9 <- c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5",
            "#C7EAE5", "#80CDC1", "#35978F", "#01665E")
.pal_wet_positive <- function(n = 101) grDevices::colorRampPalette(.BRBG9)(n)
.pal_dry_positive <- function(n = 101) grDevices::colorRampPalette(rev(.BRBG9))(n)
.pal_thermal      <- function(n = 101) grDevices::hcl.colors(n, "RdBu", rev = TRUE)

# Which palette (and legend direction hint) each variable gets. Anything not listed
# falls back to the thermal ramp. Edit here to re-assign a variable.
.TREND_PALETTE <- list(
  Precipitation = list(pal = .pal_wet_positive, hint = "\u2190 drier  |  wetter \u2192"),
  PM_VPD        = list(pal = .pal_dry_positive, hint = "\u2190 wetter  |  drier \u2192"),
  SWE           = list(pal = .pal_wet_positive, hint = "\u2190 less snow  |  more snow \u2192"),
  SWE_Peak      = list(pal = .pal_wet_positive, hint = "\u2190 less snow  |  more snow \u2192"),
  Runoff        = list(pal = .pal_wet_positive, hint = "\u2190 less runoff  |  more runoff \u2192"),
  Snowfall      = list(pal = .pal_wet_positive, hint = "\u2190 less snowfall  |  more snowfall \u2192"),
  SnowDepth     = list(pal = .pal_wet_positive, hint = "\u2190 shallower  |  deeper \u2192"),
  SnowDensity   = list(pal = .pal_thermal,      hint = "\u2190 less dense  |  denser \u2192"),
  AET           = list(pal = .pal_thermal,      hint = "\u2190 less ET  |  more ET \u2192"))
.trend_palette_for <- function(v) {
  x <- .TREND_PALETTE[[v]]
  if (is.null(x)) list(pal = .pal_thermal, hint = NULL) else x
}

# At most `max_n` round-number breaks inside the (symmetric) colour limits, so the tick
# labels of a single colourbar can never run into each other.
.nice_breaks <- function(lim, max_n = 5L) {
  for (n in 6:2) {
    b <- pretty(lim, n = n)
    b <- b[b >= lim[1] - 1e-12 & b <= lim[2] + 1e-12]
    if (length(b) >= 3L && length(b) <= max_n) return(b)
  }
  c(lim[1], 0, lim[2])
}

# One common number of decimals for all labels of a bar (no "0.0005 / 0.0000" mix, no "-0").
.break_labels <- function(x) {
  x0 <- ifelse(is.finite(x) & abs(x) < 1e-12, 0, x)
  dec <- 8L
  for (d in 0:8) {
    if (all(abs(x0 - round(x0, d)) < 1e-9 * pmax(1, abs(x0)), na.rm = TRUE)) { dec <- d; break }
  }
  formatC(x0, format = "f", digits = dec)
}

# Fixed-length horizontal colourbar. Its length is set explicitly (in cm) so that it is
# always narrower than the panel above it; it is version-safe across ggplot2 3.4 / 3.5+.
.cbar_guide <- function(len_cm = 5.5, h_cm = 0.28) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::guide_colourbar(
      theme = ggplot2::theme(
        legend.key.width      = ggplot2::unit(len_cm, "cm"),
        legend.key.height     = ggplot2::unit(h_cm, "cm"),
        legend.title.position = "top"))
  } else {
    ggplot2::guide_colourbar(
      title.position = "top",
      barwidth  = ggplot2::unit(len_cm, "cm"),
      barheight = ggplot2::unit(h_cm, "cm"))
  }
}

# Plain-text unit string ("MJ m-2 d-1 yr-1", "degC yr-1") -> plotmath string with true
# superscripts (MJ m^-2 d^-1 yr^-1) and a proper degree sign. Every base is quoted so
# plotmath never mistakes a unit for a symbol; tokens it cannot parse are kept verbatim.
.unit_plotmath <- function(unit_lab) {
  toks <- strsplit(trimws(unit_lab), "\\s+")[[1]]
  parts <- vapply(toks, function(t) {
    if (t == "degC") return("DEGC_")                  # swapped for the DEGC glyph below
    m <- regmatches(t, regexec("^([A-Za-z]+)(-?[0-9]+)$", t))[[1]]
    if (length(m) == 3L) return(sprintf("'%s'^%s", m[2], m[3]))
    sprintf("'%s'", gsub("['\\\\]", "", t))   # unparsed token: keep as literal text
  }, character(1))
  paste(parts, collapse = "~")
}

# Bold, centred legend title "Sen's slope (<unit with superscripts>)" with an optional
# second line (the drier/wetter hint). Returns an expression, which ggplot2 draws with
# plotmath; bold() is explicit because element_text(face =) is ignored for expressions.
.sen_legend_title <- function(unit_lab, hint = NULL) {
  l1 <- parse(text = sprintf("bold(\"Sen's slope (\"*%s*\")\")", .unit_plotmath(unit_lab)))[[1]]
  l1 <- do.call(substitute, list(l1, list(DEGC_ = DEGC)))   # real degree-sign glyph (DEGC)
  # The hint (contains arrows) is spliced in as an object, not parsed from text, so its
  # Unicode characters are never re-encoded by the locale.
  e <- if (!is.null(hint) && nzchar(hint)) bquote(atop(.(l1), bold(.(hint)))) else l1
  as.expression(e)
}

.trend_tile_panel <- function(s, title, unit_lab, zlim, dx, dy,
                              palette_fn = .pal_thermal, hint = NULL, bar_len_cm = 5.5) {
  s <- sf::st_drop_geometry(s)
  n_total <- sum(is.finite(s$sl_tfpw))
  n_both  <- sum(s$sig_both, na.rm = TRUE)
  n_vc1   <- sum(s$sig_fdr_vc & !s$sig_fdr_tfpw, na.rm = TRUE)
  n_tf1   <- sum(s$sig_fdr_tfpw & !s$sig_fdr_vc, na.rm = TRUE)
  pct     <- if (n_total > 0) 100 * n_both / n_total else NA_real_
  leg_title <- .sen_legend_title(unit_lab, hint)   # true superscripts via plotmath
  p <- ggplot2::ggplot(s, ggplot2::aes(x = plot_x, y = plot_y, fill = sl_tfpw)) +
    ggplot2::geom_tile(width = dx, height = dy) +
    ggplot2::scale_fill_gradientn(
      colours = palette_fn(101),
      limits = zlim, oob = scales::squish,
      breaks = .nice_breaks(zlim), labels = .break_labels,
      name = leg_title,
      guide = .cbar_guide(bar_len_cm))
  sig1 <- s[s$sig_one %in% TRUE, , drop = FALSE]     # one test only: open points
  if (nrow(sig1))
    p <- p + ggplot2::geom_point(data = sig1, ggplot2::aes(x = plot_x, y = plot_y),
                                 inherit.aes = FALSE, shape = 21, size = 0.65,
                                 stroke = 0.25, fill = NA, colour = "black")
  sig <- s[s$sig_both %in% TRUE, , drop = FALSE]     # both VC and TFPW: filled points
  if (nrow(sig))
    p <- p + ggplot2::geom_point(data = sig, ggplot2::aes(x = plot_x, y = plot_y),
                                 inherit.aes = FALSE, shape = 21, size = 0.65,
                                 stroke = 0.15, fill = "black", colour = "black")
  # Internal sub-basin boundaries: geographic context only.
  if (!is.null(subbasin_wgs) && nrow(subbasin_wgs) > 0L)
    p <- p + ggplot2::geom_sf(data = subbasin_wgs, inherit.aes = FALSE, fill = NA,
                              colour = "grey25", linewidth = 0.35, linetype = "dashed")
  p +
    ggplot2::geom_sf(data = basin_wgs, inherit.aes = FALSE, fill = NA,
                     colour = "black", linewidth = 0.65) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("BH-FDR: both %d/%d (%.1f%%) | VC only %d | TFPW only %d",
                                     n_both, n_total, pct, n_vc1, n_tf1),
                  x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 7, hjust = 0.5, colour = "grey30",
                                            margin = ggplot2::margin(b = 2)),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = ggplot2::element_text(size = 7, face = "bold", hjust = 0.5,
                                           lineheight = 1.0),
      legend.text = ggplot2::element_text(size = 6),
      legend.margin = ggplot2::margin(t = 4, b = 2),
      plot.margin = ggplot2::margin(2, 3, 2, 3))
}

.blank_panel <- function(msg = NULL) {
  p <- ggplot2::ggplot() + ggplot2::theme_void()
  if (!is.null(msg))
    p <- p + ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg, size = 4)
  p
}

####################################################################################
# FUNCTION 0: HYDROCLIMATIC TREND MAP (Pr, PET-PM, PET-ERA5, Tair)
####################################################################################
create_supp_hydroclimatic_trend_figure <- function() {
  cat("  Hydroclimatic Sen-slope/FDR map: ")
  vars <- c("Precipitation", "PET", "PET_ERA5", "Temperature")
  labels <- c(Precipitation = "Precipitation", PET = "PET-PM",
              PET_ERA5 = "PET-ERA5-Land PEV", Temperature = "2-m air temperature")
  pts <- .prep_tfpw_points(vars)
  if (is.null(pts)) { cat("no valid annual TFPW trend rows found\n"); return(invisible(NULL)) }
  
  # One symmetric range for the mm/yr variables, one for degC/yr.
  mm_z   <- max(abs(pts$sl_tfpw[pts$variable %in% c("Precipitation", "PET", "PET_ERA5")]),
                na.rm = TRUE)
  temp_z <- max(abs(pts$sl_tfpw[pts$variable == "Temperature"]), na.rm = TRUE)
  if (!is.finite(mm_z)   || mm_z   <= 0) mm_z   <- 1
  if (!is.finite(temp_z) || temp_z <= 0) temp_z <- 0.1
  dx <- .grid_spacing(pts$plot_x)
  dy <- .grid_spacing(pts$plot_y)
  
  panels <- lapply(vars, function(v) {
    s <- pts[pts$variable == v, ]
    if (!nrow(s)) return(.blank_panel(paste(labels[[v]], "\nno valid data")))
    is_t <- v == "Temperature"
    .trend_tile_panel(s, labels[[v]], if (is_t) "degC/yr" else "mm/yr",
                      if (is_t) c(-temp_z, temp_z) else c(-mm_z, mm_z), dx, dy)
  })
  
  fig <- patchwork::wrap_plots(panels, nrow = 2, ncol = 2, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  fig <- fig + patchwork::plot_annotation(
    title = "Spatial hydroclimatic trends \u2014 Nechako River Basin (1950\u20132025)",
    subtitle = paste0("Annual TFPW Mann-Kendall Sen's slope; filled black points = BH-FDR-significant under BOTH VC-MK and TFPW-MK, open points = one test only ",
                      "pixels (alpha = 0.05). Dashed lines = sub-basins.\n",
                      "ERA5-Land PEV exceeds radiation-limited evaporation and is shown as a PET-formulation sensitivity case."),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0.5, colour = "grey30",
                                            margin = ggplot2::margin(b = 5))))
  saved <- .save_gg(file.path(out_dir, "FigSx_Spatial_Hydroclimatic_Trends_SenSlope"),
                    fig, 10.5, 8.8)
  cat(paste(saved, collapse = " / "), "\n")
  invisible(fig)
}

####################################################################################
# FUNCTION 0b: PRECIPITATION + PM-PET METEOROLOGICAL-CONTROL TREND MAP
####################################################################################
create_pm_control_trend_figure <- function() {
  cat("  Precipitation + PM-control Sen-slope/FDR map: ")
  vars <- c("Precipitation", PM_CONTROL_PRESENT)
  var_labels <- c(Precipitation = "Precipitation", PM_CONTROL_LABELS)
  var_units  <- c(Precipitation = "mm yr-1", PM_CONTROL_UNITS)
  
  pts <- .prep_tfpw_points(vars)
  if (is.null(pts)) { cat("no valid annual TFPW trend rows found\n"); return(invisible(NULL)) }
  dx <- .grid_spacing(pts$plot_x)
  dy <- .grid_spacing(pts$plot_y)
  
  panels <- lapply(seq_along(vars), function(i) {
    v <- vars[[i]]
    panel_label <- sprintf("%s) %s", letters[[i]], var_labels[[v]])
    s <- pts[pts$variable == v, ]
    if (!nrow(s)) return(.blank_panel(paste(panel_label, "\nno valid data")))
    z <- max(abs(s$sl_tfpw), na.rm = TRUE)
    if (!is.finite(z) || z <= 0) z <- 1
    # Moisture variables: brown = drying, blue-green = wetting; others: thermal red/blue.
    pal <- .trend_palette_for(v)
    .trend_tile_panel(s, panel_label, var_units[[v]], c(-z, z), dx, dy,
                      palette_fn = pal$pal, hint = pal$hint, bar_len_cm = 5.5)
  })
  if (length(panels) < 6L) {  # keep the 2 x 3 geometry
    panels <- c(panels, replicate(6L - length(panels), .blank_panel(), simplify = FALSE))
  }
  
  # Each panel has its own units/scale, so every panel keeps its own colourbar
  # (guides = "keep"); the bars are fixed at 5.5 cm so they never overlap a neighbour.
  fig <- patchwork::wrap_plots(panels, nrow = 2, ncol = 3, guides = "keep") &
    ggplot2::theme(legend.position = "bottom")
  caption_vars <- paste0(letters[seq_along(vars)], ") ", unname(var_labels[vars]))
  fig <- fig + patchwork::plot_annotation(
    title = "Spatial trends in precipitation and PM-PET meteorological controls \u2014 Nechako River Basin (1950\u20132025)",
    subtitle = paste0("Annual TFPW Mann-Kendall Sen's slope; filled black points = BH-FDR-significant under BOTH VC-MK and TFPW-MK, open points = one test only ",
                      "pixels (alpha = 0.05). Dashed lines = sub-basins. ",
                      "Controls are the fields used by the FAO-56 PM calculation; PM-PET itself is excluded."),
    caption = paste(caption_vars, collapse = "; "),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0.5, colour = "grey30",
                                            margin = ggplot2::margin(b = 5))))
  saved <- .save_gg(file.path(out_dir, "FigSx_Spatial_PM_Meteorological_Controls_TFPW_SenSlope"),
                    fig, 12.5, 8.8)
  cat(paste(saved, collapse = " / "), "\n")
  invisible(fig)
}

####################################################################################
# FUNCTION 1: Publication map per variable (data-driven colour limits)
####################################################################################
create_publication_figure <- function(var_name, method = "vc") {
  rlist <- switch(var_name, PET = r_pet, PET_ERA5 = r_pet_era5,
                  Precipitation = r_precip, Temperature = r_tair,
                  SWE = r_swe, SWE_Peak = r_swe_peak, NULL)
  if (is.null(rlist) && var_name %in% names(r_extra)) rlist <- r_extra[[var_name]]
  if (is.null(rlist)) { cat("  No rasters for", var_name, "- skipping\n"); return(NULL) }
  pdf_path <- file.path(out_dir, sprintf("%s_%s_publication_map.pdf",
                                         gsub("_", "", tolower(var_name)), method))
  r_comb <- rlist[[paste0("comb_", method)]]
  r_dir  <- rlist[[paste0("dir_",  method)]]
  r_mag  <- rlist[[paste0("mag_",  method)]]
  r_tau  <- rlist[[paste0("tau_",  method)]]
  
  is_temp  <- var_name == "Temperature"
  unit_lab <- if (is_temp) paste0(DEGC, "/yr") else
    if (var_name %in% names(EXTRA_SLOPE_UNITS)) sub(" yr-1$", "/yr", EXTRA_SLOPE_UNITS[[var_name]]) else "mm/yr"
  dir_lab  <- if (is_temp) c("Cooling", "Non-sig.", "Warming")
  else         c("Decreasing", "Non-sig.", "Increasing")
  dir_col  <- c("#4575b4", "#f0f0f0", "#d73027")
  rdbu     <- grDevices::hcl.colors(101, "RdBu", rev = TRUE)
  
  result <- .pdf_safe(pdf_path, width = 10, height = 6.67, expr_fn = function() {
    par(mfrow = c(1, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
    if (is_temp) {
      plot_raster_panel(r_tau, paste("Kendall's", TAU), .sym_lim(r_tau, 0.3), rdbu,
                        legend_title = TAU)
    } else {
      plot_raster_panel(r_comb, "Combined trend strength", .pos_lim(r_comb, 0.25),
                        grDevices::hcl.colors(101, "YlOrRd"),
                        legend_title = paste0("|", TAU, "|"))
      mtext("|tau| shown only where p < 0.05; others zero", side = 3, line = 0.1,
            cex = 0.7, col = "gray50")
    }
    plot_raster_panel(r_mag, paste0("Sen's slope (", unit_lab, ")"), .sym_lim(r_mag, 1),
                      rdbu, legend_title = unit_lab)
    plot_raster_panel(r_dir, "Direction (p < 0.05)", NULL, dir_col,
                      breaks = c(-1.5, -0.5, 0.5, 1.5), legend = FALSE, categorical = TRUE)
    legend("bottomright", legend = dir_lab, fill = dir_col, bty = "n", cex = 1.1)
    mtext(sprintf("%s Annual Trends (1950-2025) | %s", var_name,
                  ifelse(method == "vc", "Variance-Corrected MK", "TFPW-Corrected MK")),
          outer = TRUE, cex = 1.1, font = 2, line = 0.5)
  })
  if (!is.null(result)) cat(sprintf("  %s %s map: %s\n", var_name, method, basename(result)))
  invisible(result)
}

####################################################################################
# FUNCTION 2: 4-variable Kendall-tau comparison map (TFPW - VC is ~fully filtered)
####################################################################################
create_comparison_figure <- function() {
  rl <- list(Precipitation = r_precip, PET_PM = r_pet,
             PET_ERA5 = r_pet_era5, Temperature = r_tair, SWE = r_swe)
  lab <- c(Precipitation = paste("Precipitation", TAU), PET_PM = paste("PET_PM", TAU),
           PET_ERA5 = paste("PET_ERA5", TAU), Temperature = paste("Temperature", TAU),
           SWE = paste("SWE (Mar-Apr)", TAU))
  ok <- vapply(rl, function(x)
    !is.null(x) && !is.null(x$tau_tfpw) && any(is.finite(terra::values(x$tau_tfpw))),
    logical(1))
  if (!any(ok)) { cat("  No variables available for comparison map.\n"); return(NULL) }
  rl <- rl[ok]
  zl <- max(vapply(rl, function(x) max(abs(.sym_lim(x$tau_tfpw, 0.3))), numeric(1)))
  n  <- length(rl)
  result <- .pdf_safe(file.path(out_dir, "variables_trend_comparison.pdf"),
                      width = 5 * n, height = 6, expr_fn = function() {
                        par(mfrow = c(1, n), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
                        for (v in names(rl))
                          plot_raster_panel(rl[[v]]$tau_tfpw, lab[[v]], c(-zl, zl),
                                            grDevices::hcl.colors(101, "RdBu", rev = TRUE), legend_title = TAU)
                        ttl <- "Nechako Basin Annual Trends (1950-2025) | TFPW-MK"
                        if (n < 4) ttl <- paste0(ttl, "\nMissing: ", paste(setdiff(names(ok), names(rl)), collapse = ", "))
                        mtext(ttl, outer = TRUE, cex = 1.2, font = 2, line = 0.5)
                      })
  if (!is.null(result)) cat(sprintf("  %d-variable comparison map: %s\n", n, basename(result)))
  invisible(result)
}

####################################################################################
# FUNCTION 3: Regime shift maps
####################################################################################
create_regime_shift_maps <- function() {
  need <- c("changepoint_detected", "first_changepoint_year")
  if (!all(need %in% names(all_results))) {
    cat("  Changepoint columns missing - skipping regime shift maps.\n"); return(NULL)
  }
  vars_here <- MAP_VARS; gd <- .grid_dims(length(vars_here))
  result <- .pdf_safe(file.path(out_dir, "regime_shifts_changepoints.pdf"),
                      width = 5 * gd[["nc"]], height = 5 * gd[["nr"]], expr_fn = function() {
                        par(mfrow = c(gd[["nr"]], gd[["nc"]]), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
                        for (var in vars_here) {
                          vd <- all_results[variable == var & period == "annual" & !is_basin_average &
                                              changepoint_detected == TRUE]
                          if (!nrow(vd)) { plot.new(); text(0.5, 0.5, paste("No regime shifts\n", var), cex = 1.2); next }
                          vd <- data.table::copy(vd)
                          vd[, shift_decade_num := as.numeric(cut(first_changepoint_year,
                                                                  breaks = seq(1950, 2030, by = 10),
                                                                  include.lowest = TRUE))]
                          r_dec <- create_raster_from_table(vd, template_bc, "shift_decade_num")
                          plot_raster_panel(r_dec, paste0(var, ": regime-shift decade"), c(1, 8),
                                            grDevices::hcl.colors(8, "Spectral", rev = TRUE),
                                            legend_title = "1 = 1950s")
                        }
                        mtext("Regime Shift Detection (PELT) - 1950-2025", outer = TRUE, cex = 1.1, font = 2, line = 0.5)
                      })
  if (!is.null(result)) cat(sprintf("  Regime shift map: %s\n", basename(result)))
  invisible(result)
}

####################################################################################
# FUNCTION 4: Spectral maps
####################################################################################
create_spectral_analysis_maps <- function() {
  need <- c("n_spectral_peaks", "dominant_period")
  if (!all(need %in% names(all_results))) {
    cat("  Spectral columns missing - skipping spectral maps.\n"); return(NULL)
  }
  vars_here <- MAP_VARS; gd <- .grid_dims(length(vars_here))
  result <- .pdf_safe(file.path(out_dir, "spectral_analysis.pdf"),
                      width = 5 * gd[["nc"]], height = 5 * gd[["nr"]], expr_fn = function() {
                        par(mfrow = c(gd[["nr"]], gd[["nc"]]), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
                        for (var in vars_here) {
                          vd <- all_results[variable == var & period == "annual" & !is_basin_average &
                                              n_spectral_peaks > 0]
                          if (!nrow(vd)) { plot.new(); text(0.5, 0.5, paste("No spectral peaks\n", var), cex = 1.2); next }
                          r_p <- create_raster_from_table(vd, template_bc, "dominant_period")
                          plot_raster_panel(r_p, paste0(var, ": dominant cycle (yr)"),
                                            range(vd$dominant_period, na.rm = TRUE),
                                            grDevices::hcl.colors(101, "plasma"), legend_title = "Years")
                        }
                        mtext("Spectral Analysis: Dominant Periodicities", outer = TRUE, cex = 1.1, font = 2, line = 0.5)
                      })
  if (!is.null(result)) cat(sprintf("  Spectral map: %s\n", basename(result)))
  invisible(result)
}

####################################################################################
# FUNCTION 5: Spatial pattern analysis
# The old "CV" (|slope| / (|mean slope| + 0.001)) was unit-dependent and not a CV.
# Replaced by the standardised departure of each pixel's TFPW Sen slope from the
# basin mean (z-score), which is unit-free and comparable across variables.
####################################################################################
create_spatial_pattern_analysis <- function() {
  vars_here <- MAP_VARS; gd <- .grid_dims(length(vars_here))
  result <- .pdf_safe(file.path(out_dir, "spatial_patterns.pdf"),
                      width = 5 * gd[["nc"]], height = 5 * gd[["nr"]], expr_fn = function() {
                        par(mfrow = c(gd[["nr"]], gd[["nc"]]), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
                        for (var in vars_here) {
                          vd <- all_results[variable == var & period == "annual" & !is_basin_average &
                                              !filtered_tfpw & is.finite(sl_tfpw)]
                          if (nrow(vd) < 3L) { plot.new(); text(0.5, 0.5, paste("No data\n", var), cex = 1.2); next }
                          vd <- data.table::copy(vd)
                          s <- stats::sd(vd$sl_tfpw)
                          if (!is.finite(s) || s <= 0) { plot.new(); text(0.5, 0.5, paste("No spatial variance\n", var), cex = 1.2); next }
                          vd[, spatial_z := (sl_tfpw - mean(sl_tfpw)) / s]
                          r_z <- create_raster_from_table(vd, template_bc, "spatial_z")
                          plot_raster_panel(r_z, paste0(var, ": slope anomaly (z)"), c(-3, 3),
                                            grDevices::hcl.colors(101, "RdBu", rev = TRUE), legend_title = "z")
                        }
                        mtext("Spatial variability of annual TFPW Sen slopes (z-score vs basin mean)",
                              outer = TRUE, cex = 1.1, font = 2, line = 0.5)
                      })
  if (!is.null(result)) cat(sprintf("  Spatial pattern map: %s\n", basename(result)))
  invisible(result)
}

####################################################################################
# FUNCTION 6: Method comparison (VC vs TFPW)
####################################################################################
create_method_comparison <- function() {
  vars4 <- c("Precipitation", "PET", "PET_ERA5", "Temperature", "SWE", EXTRA_VARS_PRESENT)
  comp <- all_results[period == "annual" & !is_basin_average &
                        !filtered_vc & !filtered_tfpw & variable %in% vars4 &
                        is.finite(tau_vc) & is.finite(tau_tfpw)]
  if (nrow(comp) < MIN_PIXELS_METHOD_COMPARISON) {
    cat(sprintf(paste0("  Only %d pixels are valid for BOTH VC and TFPW (< %d) - ",
                       "method comparison skipped. Most VC results are filtered upstream.\n"),
                nrow(comp), MIN_PIXELS_METHOD_COMPARISON))
    return(NULL)
  }
  comp <- data.table::copy(comp)
  comp[, vc_sig   := p_value_vc   < 0.05]
  comp[, tfpw_sig := p_value_tfpw < 0.05]
  tau_range <- range(c(comp$tau_vc, comp$tau_tfpw), na.rm = TRUE)
  var_cols <- c(Precipitation = "#4575b4", PET = "#d73027",
                PET_ERA5 = "#e08214", Temperature = "#1a9641", SWE = "#01665e", EXTRA_COLOURS)
  p1 <- ggplot(comp, aes(x = tau_vc, y = tau_tfpw, colour = variable)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray40") +
    geom_point(alpha = 0.35, size = 1.0) +
    facet_wrap(~variable, labeller = labeller(variable = c(
      Precipitation = "Precipitation", PET = "PET (PM)",
      PET_ERA5 = "PET (ERA5)", Temperature = "Temperature", SWE = "SWE (Mar-Apr)",
      EXTRA_VAR_LABELS))) +
    scale_colour_manual(values = var_cols) +
    coord_equal(xlim = tau_range, ylim = tau_range) +
    theme_bw(base_size = 10) + theme(legend.position = "none") +
    labs(title = paste0("A) Kendall's ", TAU, ": VC vs TFPW"),
         x = paste0("VC method (", TAU, ")"), y = paste0("TFPW method (", TAU, ")"))
  
  ag <- comp[, .(Both = sum(vc_sig & tfpw_sig, na.rm = TRUE),
                 Only_VC = sum(vc_sig & !tfpw_sig, na.rm = TRUE),
                 Only_TFPW = sum(!vc_sig & tfpw_sig, na.rm = TRUE),
                 Neither = sum(!vc_sig & !tfpw_sig, na.rm = TRUE)), by = variable]
  ag_l <- data.table::melt(ag, id.vars = "variable", variable.name = "category",
                           value.name = "n_pix")                     # [FIX #4]
  ag_l[, pct := 100 * n_pix / sum(n_pix), by = variable]             # % within variable
  ag_l[, category := factor(category, levels = c("Both", "Only_VC", "Only_TFPW", "Neither"))]
  p2 <- ggplot(ag_l, aes(x = category, y = n_pix, fill = category)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.0f%%", pct)), vjust = -0.3, size = 2.8) +
    facet_wrap(~variable, scales = "free_y") +
    scale_fill_manual(values = c(Both = "#1b7837", Only_VC = "#762a83",
                                 Only_TFPW = "#e08214", Neither = "grey70")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom") +
    labs(title = "B) Significance agreement (p < 0.05)", x = NULL, y = "N pixels", fill = NULL)
  
  combined <- p1 / p2 +
    plot_annotation(title = "Mann-Kendall Method Comparison (VC vs TFPW)",
                    theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)))
  saved <- .save_gg(file.path(out_dir, "method_comparison_vc_vs_tfpw"), combined, 14, 10,
                    formats = "png")
  cat("  Method comparison plot saved:", paste(saved, collapse = ", "), "\n")
  invisible(combined)
}

####################################################################################
# FUNCTION 7: Statistics table
####################################################################################
create_statistics_table <- function() {
  has_rho_vc <- "rho1_vc" %in% names(all_results)
  .fmt <- function(v, f) if (length(v) > 0 && any(is.finite(v))) sprintf(f, median(v, na.rm = TRUE)) else "N/A"
  stats_summary <- all_results[period == "annual" & !is_basin_average, .(
    N_Pixels = .N,
    N_Valid_VC = sum(!filtered_vc, na.rm = TRUE),
    N_Valid_TFPW = sum(!filtered_tfpw, na.rm = TRUE),
    Pct_Sig_VC = {
      n <- sum(!filtered_vc, na.rm = TRUE)
      if (n > 0) sprintf("%.1f%%", sum(p_value_vc < 0.05 & !filtered_vc, na.rm = TRUE) / n * 100) else "N/A"
    },
    Pct_Sig_TFPW = {
      n <- sum(!filtered_tfpw, na.rm = TRUE)
      if (n > 0) sprintf("%.1f%%", sum(p_value_tfpw < 0.05 & !filtered_tfpw, na.rm = TRUE) / n * 100) else "N/A"
    },
    Median_Tau_VC     = .fmt(tau_vc[!filtered_vc], "%.3f"),
    Median_Slope_VC   = .fmt(sl_vc[!filtered_vc], "%.4f"),
    Median_Tau_TFPW   = .fmt(tau_tfpw[!filtered_tfpw], "%.3f"),
    Median_Slope_TFPW = .fmt(sl_tfpw[!filtered_tfpw], "%.4f"),
    Mean_Autocorr = if (has_rho_vc) {
      v <- rho1_vc[!filtered_vc]
      if (length(v) > 0) sprintf("%.3f", mean(v, na.rm = TRUE)) else "N/A"
    } else NA_character_
  ), by = .(variable)]
  fwrite(stats_summary, .resolve_out(file.path(out_dir, "summary_statistics_visualization.csv")))
  cat("  Statistics table saved\n")
  print(as.data.frame(stats_summary))
  invisible(stats_summary)
}

####################################################################################
# FUNCTION 8: Basin time series (Pr + PET_PM + PET_ERA5 + Temperature + SWE + runoff/snow)
####################################################################################
create_basin_timeseries_plot <- function() {
  has_tair     <- "tair_degC_month" %in% names(basin_avg_monthly)
  has_pet_era5 <- "pet_era5_mm_month" %in% names(basin_avg_monthly)
  has_swe_ts <- "swe_mm_month" %in% names(basin_avg_monthly) && any(is.finite(basin_avg_monthly$swe_mm_month))
  extra_ts <- EXTRA_MONTH_COLS[EXTRA_MONTH_COLS %in% names(basin_avg_monthly)]
  extra_ts <- extra_ts[vapply(extra_ts, function(cl) any(is.finite(basin_avg_monthly[[cl]])), logical(1))]
  n_panels <- 2L + has_pet_era5 * 1L + has_tair * 1L + has_swe_ts * 1L + length(extra_ts)
  png_path <- .resolve_out(file.path(out_dir, "basin_average_timeseries.png"))
  png(png_path, width = 12, height = 4 * n_panels, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  layout(matrix(seq_len(n_panels), nrow = n_panels), heights = rep(1, n_panels))
  
  d <- basin_avg_monthly$date
  pr <- basin_avg_monthly$precip_mm_month
  par(mar = c(3, 5, 3.5, 2))
  plot(d, pr, type = "l", col = rgb(0.2, 0.4, 0.8, 0.3), lwd = 0.8,
       xlab = "", ylab = "Precipitation (mm/month)",
       main = "Nechako Basin Climate Variables (1950-2025)",
       ylim = c(0, max(pr, na.rm = TRUE) * 1.1), cex.lab = 1.1, cex.main = 1.2)
  lines(d, zoo::rollmean(pr, k = 12, fill = NA, align = "center"), col = "blue", lwd = 2.5)
  legend("topright", legend = c("Pr (monthly)", "Pr 12-mo mean"),
         col = c(rgb(0.2, 0.4, 0.8, 0.3), "blue"), lwd = c(0.8, 2.5), bty = "n", cex = 0.9)
  
  pet_pm <- basin_avg_monthly$pet_mm_month
  plot(d, pet_pm, type = "l", col = rgb(0.8, 0.2, 0.2, 0.25), lwd = 0.6,
       xlab = "", ylab = "PET-PM (mm/month)", ylim = c(0, max(pet_pm, na.rm = TRUE) * 1.15), cex.lab = 1.1,
       main = "Potential Evapotranspiration - FAO-56 Penman-Monteith")
  lines(d, zoo::rollmean(pet_pm, k = 12, fill = NA, align = "center"), col = "red3", lwd = 2.5)
  legend("topright", legend = c("PET_PM (monthly)", "PET_PM 12-mo"),
         col = c(rgb(0.8, 0.2, 0.2, 0.3), "red3"), lwd = c(0.6, 2.5), bty = "n", cex = 0.9)
  if (has_pet_era5) {
    pet_era5 <- basin_avg_monthly$pet_era5_mm_month
    plot(d, pet_era5, type = "l", col = rgb(0.9, 0.5, 0, 0.25), lwd = 0.6, xlab = "",
         ylab = "PEV (mm/month)", ylim = c(0, max(pet_era5, na.rm = TRUE) * 1.15), cex.lab = 1.1,
         main = "ERA5-Land PEV (exceeds radiation-limited evaporation; formulation-sensitivity case)")
    lines(d, zoo::rollmean(pet_era5, k = 12, fill = NA, align = "center"), col = "#e08214", lwd = 2.5)
    legend("topright", legend = c("PEV (monthly)", "PEV 12-mo"),
           col = c(rgb(0.9, 0.5, 0, 0.3), "#e08214"), lwd = c(0.6, 2.5), bty = "n", cex = 0.9)
  }
  
  if (has_tair) {
    tair <- basin_avg_monthly$tair_degC_month
    par(mar = c(4, 5, 2, 2))
    plot(d, tair, type = "l", col = rgb(0.6, 0.2, 0.6, 0.25), lwd = 0.6,
         xlab = "", ylab = paste0("Temperature (", DEGC, ")"), cex.lab = 1.1,
         ylim = range(tair, na.rm = TRUE) + c(-0.5, 0.5),
         main = "Basin-Mean 2-m Air Temperature")
    lines(d, zoo::rollmean(tair, k = 12, fill = NA, align = "center"), col = "purple", lwd = 2.5)
    yr_frac <- as.numeric(format(d, "%Y")) + (as.numeric(format(d, "%m")) - 0.5) / 12
    ok <- !is.na(yr_frac) & !is.na(tair)
    fit_t <- if (sum(ok) >= 10) lm(tair[ok] ~ yr_frac[ok]) else NULL
    if (!is.null(fit_t)) {
      yr_seq <- seq(min(yr_frac[ok]), max(yr_frac[ok]), length.out = 200)
      pred_y <- coef(fit_t)[1] + coef(fit_t)[2] * yr_seq
      # convert fractional year -> Date so the trend line is smooth (not stair-stepped)
      d_seq <- as.Date(paste0(floor(yr_seq), "-01-01")) + (yr_seq - floor(yr_seq)) * 365.25
      lines(d_seq, pred_y, col = "#d73027", lwd = 2, lty = 2)
      mtext(sprintf("OLS: %+.4f %s/yr", coef(fit_t)[2], DEGC),
            side = 3, line = 0, adj = 1, cex = 0.85, col = "#d73027", font = 2)
    }
    legend("topleft", legend = c("Tair (monthly)", "12-mo mean", "OLS trend"),
           col = c(rgb(0.6, 0.2, 0.6, 0.25), "purple", "#d73027"),
           lwd = c(0.6, 2.5, 2), lty = c(1, 1, 2), bty = "n", cex = 0.9)
  }
  if (has_swe_ts) {
    swe <- basin_avg_monthly$swe_mm_month
    par(mar = c(4, 5, 2, 2))
    plot(d, swe, type = "l", col = rgb(0.0, 0.4, 0.4, 0.25), lwd = 0.6,
         xlab = "", ylab = "SWE (mm w.e.)", cex.lab = 1.1,
         ylim = c(0, max(swe, na.rm = TRUE) * 1.1),
         main = "Basin-Mean Snow Water Equivalent (ERA5-Land)")
    lines(d, zoo::rollmean(swe, k = 12, fill = NA, align = "center"), col = "#01665e", lwd = 2.5)
    leg <- c("SWE (monthly)", "SWE 12-mo mean"); lcol <- c(rgb(0.0, 0.4, 0.4, 0.4), "#01665e")
    llwd <- c(0.6, 2.5); lpch <- c(NA, NA)
    if (!is.null(basin_avg_annual) && "swe_marapr_mm_year" %in% names(basin_avg_annual)) {
      points(as.Date(paste0(basin_avg_annual$year, "-04-01")), basin_avg_annual$swe_marapr_mm_year,
             pch = 19, cex = 0.7, col = "#d73027")
      leg <- c(leg, "Mar-Apr mean (annual)"); lcol <- c(lcol, "#d73027"); llwd <- c(llwd, NA); lpch <- c(lpch, 19)
    }
    legend("topright", legend = leg, col = lcol, lwd = llwd, pch = lpch, bty = "n", cex = 0.9)
  }
  # Runoff / snowfall / snow depth / snow density panels (v8)
  .roll12 <- function(v) zoo::rollapply(v, 12, function(z) if (sum(is.finite(z)) >= 6L) mean(z, na.rm = TRUE) else NA_real_,
                                        fill = NA, align = "center")
  for (nm in names(extra_ts)) {
    v <- basin_avg_monthly[[extra_ts[[nm]]]]; col <- EXTRA_COLOURS[[nm]]
    yl <- if (nm == "SnowDensity") range(v, na.rm = TRUE) else c(0, max(v, na.rm = TRUE) * 1.1)
    par(mar = c(4, 5, 2, 2))
    plot(d, v, type = "l", col = grDevices::adjustcolor(col, 0.3), lwd = 0.6, xlab = "", ylim = yl,
         ylab = sprintf("%s (%s)", .short_lab(nm), EXTRA_MONTH_UNITS[[nm]]), cex.lab = 1.1,
         main = sprintf("Basin-mean %s (ERA5-Land)", tolower(.short_lab(nm))))
    lines(d, .roll12(v), col = col, lwd = 2.5)
    legend("topright", legend = c("monthly", "12-mo running mean"),
           col = c(grDevices::adjustcolor(col, 0.5), col), lwd = c(0.6, 2.5), bty = "n", cex = 0.9)
  }
  cat(sprintf("  Basin time series: %s\n", basename(png_path)))
  invisible(png_path)
}

####################################################################################
# FUNCTION 13: CALENDAR-MONTHLY TEMPORAL CHANGE PLOTS
####################################################################################
.build_calmonth_df <- function(dates, monthly_vec) {
  data.frame(year = as.integer(format(dates, "%Y")),
             month = as.integer(format(dates, "%m")),
             value = as.numeric(monthly_vec), stringsAsFactors = FALSE)
}

create_calmonth_temporal_changes <- function(var_label, monthly_vec, dates, y_unit,
                                             color_line, color_ols, pdf_stem) {
  df <- .build_calmonth_df(dates, monthly_vec)
  df <- df[!is.na(df$value), ]
  if (!nrow(df)) { cat(sprintf("  %s: no data - skipping\n", var_label)); return(invisible(NULL)) }
  panels <- lapply(1:12, function(m) {
    sub <- df[df$month == m, ]
    n_obs <- nrow(sub)
    y_lbl <- if (m %in% c(1, 5, 9)) y_unit else NULL
    p <- ggplot2::ggplot(sub, ggplot2::aes(x = year, y = value)) +
      ggplot2::geom_point(colour = color_line, size = 0.9, alpha = 0.65) +
      ggplot2::geom_line(colour = color_line, linewidth = 0.45, alpha = 0.5)
    sig_label <- ""
    if (n_obs >= 10) {
      fit <- tryCatch(lm(value ~ year, data = sub), error = function(e) NULL)
      if (!is.null(fit)) {
        cs <- summary(fit)$coefficients
        if ("year" %in% rownames(cs)) {
          pval <- cs["year", 4]; slope <- coef(fit)[["year"]]
          pred <- data.frame(year = range(sub$year, na.rm = TRUE))
          pred$y <- predict(fit, newdata = pred)
          p <- p + ggplot2::geom_line(data = pred, ggplot2::aes(x = year, y = y),
                                      colour = color_ols, linewidth = 1.1, inherit.aes = FALSE)
          sig_label <- if (!is.na(pval) && pval < 0.05) sprintf("%+.3f/yr *", slope)
          else sprintf("%+.3f/yr", slope)
        }
      }
    }
    p +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = sig_label, hjust = -0.1, vjust = 1.4,
                        size = 2.5, colour = if (grepl("\\*", sig_label)) color_ols else "grey45",
                        fontface = "bold") +
      ggplot2::scale_x_continuous(breaks = seq(1960, 2020, by = 20)) +
      ggplot2::labs(title = month.abb[m], x = NULL, y = y_lbl, subtitle = sprintf("n=%d", n_obs)) +
      ggplot2::theme_classic(base_size = 8.5) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 9, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 7, colour = "grey50", hjust = 0.5),
        axis.title.y = ggplot2::element_text(size = 8),
        axis.text = ggplot2::element_text(size = 7),
        panel.grid.major = ggplot2::element_line(colour = "grey94", linewidth = 0.25),
        plot.margin = ggplot2::margin(2, 3, 2, 3))
  })
  fig <- patchwork::wrap_plots(panels, ncol = 4) +
    patchwork::plot_annotation(
      title = sprintf("Calendar-Monthly Temporal Changes \u2014 %s \u2014 Nechako Basin", var_label),
      subtitle = paste0("Each panel = one calendar month; each dot = one year's value. ",
                        "Solid line = OLS trend. * p < 0.05. Slope in ", y_unit, "/yr."),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35", hjust = 0)))
  saved <- .save_gg(file.path(out_dir, pdf_stem), fig, 13, 9)
  cat(sprintf("  %s: %s\n", var_label, paste(saved, collapse = " / ")))
  invisible(fig)
}

create_all_calmonth_temporal_changes <- function() {
  cat("\n-- Function 13: Calendar-monthly temporal change plots --\n")
  d <- basin_avg_monthly$date
  create_calmonth_temporal_changes("Precipitation", basin_avg_monthly$precip_mm_month, d,
                                   "mm/month", "#4575b4", "#08306b", "calmonth_temporal_Precipitation")
  create_calmonth_temporal_changes("PET \u2014 Penman-Monteith", basin_avg_monthly$pet_mm_month, d,
                                   "mm/month", "#d73027", "#67000d", "calmonth_temporal_PET_PM")
  if ("pet_era5_mm_month" %in% names(basin_avg_monthly))
    create_calmonth_temporal_changes("PET \u2014 ERA5-Land PEV", basin_avg_monthly$pet_era5_mm_month, d,
                                     "mm/month", "#e08214", "#7f2704", "calmonth_temporal_PET_ERA5")
  else cat("  13c: pet_era5_mm_month not found - skipped\n")
  if ("tair_degC_month" %in% names(basin_avg_monthly))
    create_calmonth_temporal_changes("Temperature (2-m)", basin_avg_monthly$tair_degC_month, d,
                                     DEGC, "#6a3d9a", "#3d0072", "calmonth_temporal_Temperature")
  else cat("  13d: tair_degC_month not found - skipped\n")
  if ("swe_mm_month" %in% names(basin_avg_monthly) && any(is.finite(basin_avg_monthly$swe_mm_month)))
    create_calmonth_temporal_changes("Snow water equivalent", basin_avg_monthly$swe_mm_month, d,
                                     "mm w.e.", "#01665e", "#00342e", "calmonth_temporal_SWE")
  else cat("  13f: swe_mm_month not found - skipped\n")
  for (nm in names(EXTRA_MONTH_COLS)) if (EXTRA_MONTH_COLS[[nm]] %in% names(basin_avg_monthly))
    create_calmonth_temporal_changes(.short_lab(nm), basin_avg_monthly[[EXTRA_MONTH_COLS[[nm]]]], d,
                                     EXTRA_MONTH_UNITS[[nm]], EXTRA_COLOURS[[nm]], "black",
                                     paste0("calmonth_temporal_", nm))
  
  cat("  13e: combined per-month slope profile...\n")
  .slope_profile_df <- function(dates, vals, var_label, unit) {
    do.call(rbind, lapply(1:12, function(m) {
      yr <- as.integer(format(dates, "%Y")); mo <- as.integer(format(dates, "%m"))
      sub <- data.frame(year = yr[mo == m], value = vals[mo == m])
      sub <- sub[!is.na(sub$value), ]
      empty <- data.frame(month = m, month_abb = month.abb[m], slope = NA_real_, pval = NA_real_,
                          sig = FALSE, variable = var_label, unit = unit, stringsAsFactors = FALSE)
      if (nrow(sub) < 10) return(empty)
      fit <- tryCatch(lm(value ~ year, data = sub), error = function(e) NULL)
      if (is.null(fit)) return(empty)
      cs <- summary(fit)$coefficients
      data.frame(month = m, month_abb = month.abb[m], slope = unname(cs["year", 1]),
                 pval = unname(cs["year", 4]), sig = unname(cs["year", 4] < 0.05),
                 variable = var_label, unit = unit, stringsAsFactors = FALSE)
    }))
  }
  sp_list <- list(
    .slope_profile_df(d, basin_avg_monthly$precip_mm_month, "Precipitation", "mm/month"),
    .slope_profile_df(d, basin_avg_monthly$pet_mm_month, "PET_PM", "mm/month"))
  if ("pet_era5_mm_month" %in% names(basin_avg_monthly))
    sp_list[[length(sp_list) + 1]] <- .slope_profile_df(d, basin_avg_monthly$pet_era5_mm_month, "PET_ERA5", "mm/month")
  if ("tair_degC_month" %in% names(basin_avg_monthly))
    sp_list[[length(sp_list) + 1]] <- .slope_profile_df(d, basin_avg_monthly$tair_degC_month, "Temperature", DEGC)
  if ("swe_mm_month" %in% names(basin_avg_monthly) && any(is.finite(basin_avg_monthly$swe_mm_month)))
    sp_list[[length(sp_list) + 1]] <- .slope_profile_df(d, basin_avg_monthly$swe_mm_month, "SWE", "mm w.e.")
  for (nm in names(EXTRA_MONTH_COLS)) if (EXTRA_MONTH_COLS[[nm]] %in% names(basin_avg_monthly))
    sp_list[[length(sp_list) + 1]] <- .slope_profile_df(d, basin_avg_monthly[[EXTRA_MONTH_COLS[[nm]]]],
                                                        nm, EXTRA_MONTH_UNITS[[nm]])
  sp_df <- do.call(rbind, sp_list)
  sp_df$month_abb <- factor(sp_df$month_abb, levels = month.abb)
  sp_df$variable  <- factor(sp_df$variable, levels = c("Precipitation", "PET_PM", "PET_ERA5", "Temperature", "SWE",
                                                       names(EXTRA_MONTH_COLS)))
  fwrite(as.data.table(sp_df), .resolve_out(file.path(out_dir, "calmonth_ols_slope_profiles.csv")))
  
  var_pal <- c(Precipitation = "#4575b4", PET_PM = "#d73027", PET_ERA5 = "#e08214",
               Temperature = "#6a3d9a", SWE = "#01665e", EXTRA_COLOURS)
  extra_lab <- setNames(paste0(vapply(names(EXTRA_VAR_LABELS), .short_lab, character(1)), " (",
                               EXTRA_MONTH_UNITS[names(EXTRA_VAR_LABELS)], "/yr)"), names(EXTRA_VAR_LABELS))
  # star offset computed per variable so it sits just beyond each bar
  sp_df$star_y <- ave(sp_df$slope, sp_df$variable,
                      FUN = function(s) s + sign(s) * max(abs(s), na.rm = TRUE) * 0.07)
  sp_sig <- sp_df[which(sp_df$sig), ]
  p13e <- ggplot2::ggplot(sp_df, ggplot2::aes(x = month_abb, y = slope, fill = variable, alpha = sig)) +
    ggplot2::geom_col(position = "dodge", colour = "white", linewidth = 0.2, width = 0.8) +
    ggplot2::geom_point(data = sp_sig, ggplot2::aes(x = month_abb, y = star_y, colour = variable),
                        shape = 8, size = 2.0, inherit.aes = FALSE, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.45) +
    ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y", labeller = ggplot2::labeller(variable = c(
      Precipitation = "Precipitation (mm/month/yr)",
      PET_PM = "PET \u2014 Penman-Monteith (mm/month/yr)",
      PET_ERA5 = "PET \u2014 ERA5-Land PEV (mm/month/yr)",
      Temperature = paste0("Temperature (", DEGC, "/yr)"),
      SWE = "SWE (mm w.e./yr)",
      extra_lab))) +
    ggplot2::scale_fill_manual(values = var_pal, guide = "none") +
    ggplot2::scale_colour_manual(values = var_pal, guide = "none") +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.40), guide = "none") +
    ggplot2::labs(title = "Per-Calendar-Month OLS Slope Profile \u2014 All Variables",
                  subtitle = "Solid bars = significant (p < 0.05); faded = non-significant. * marks p < 0.05.",
                  x = "Calendar month", y = "OLS slope (unit/yr)") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5),
                   strip.text = ggplot2::element_text(face = "bold", size = 9),
                   axis.title.x = ggplot2::element_text(size = 9),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
                   plot.margin = ggplot2::margin(4, 10, 4, 4))
  saved <- .save_gg(file.path(out_dir, "calmonth_slope_profile_all_vars"), p13e, 12,
                    2 * nlevels(droplevels(sp_df$variable)) + 1)
  cat(sprintf("  13e combined slope profile: %s\n", paste(saved, collapse = " / ")))
  cat("-- Function 13 complete.\n")
  invisible(NULL)
}

####################################################################################
# FUNCTION 14: PET_PM vs PET_ERA5 COMPARISON
####################################################################################
create_pet_comparison_plots <- function() {
  if (!has_pet_era5 || !"pet_era5_mm_month" %in% names(basin_avg_monthly)) {
    cat("\n-- PET_ERA5 not present - skipping Function 14.\n"); return(invisible(NULL))
  }
  cat("\n-- Function 14: PET_PM vs PET_ERA5 comparison plots --\n")
  pet_pm   <- basin_avg_monthly$pet_mm_month
  pet_era5 <- basin_avg_monthly$pet_era5_mm_month
  dates <- basin_avg_monthly$date
  yr <- as.integer(format(dates, "%Y")); mo <- as.integer(format(dates, "%m"))
  theme_cmp <- ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey40", hjust = 0),
                   panel.grid.major = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
                   plot.margin = ggplot2::margin(4, 8, 4, 4))
  
  cat("  14a: time series overlay...\n")
  df_ts <- data.frame(date = dates, PET_PM = pet_pm, PET_ERA5 = pet_era5,
                      pm_roll = zoo::rollmean(pet_pm, k = 12, fill = NA, align = "center"),
                      era5_roll = zoo::rollmean(pet_era5, k = 12, fill = NA, align = "center"))
  p14a <- ggplot2::ggplot(df_ts, ggplot2::aes(x = date)) +
    ggplot2::geom_line(ggplot2::aes(y = PET_PM), colour = rgb(0.84, 0.19, 0.15, 0.25), linewidth = 0.45) +
    ggplot2::geom_line(ggplot2::aes(y = PET_ERA5), colour = rgb(0.88, 0.51, 0.08, 0.25), linewidth = 0.45) +
    ggplot2::geom_line(ggplot2::aes(y = pm_roll), colour = "#d73027", linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = era5_roll), colour = "#e08214", linewidth = 1.2,
                       linetype = "dashed", na.rm = TRUE) +
    ggplot2::scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                          expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(title = "Basin-Mean PET \u2014 Penman-Monteith vs ERA5-Land PEV (1950\u20132025)",
                  subtitle = "Thick lines = 12-month rolling mean. Red = PM, orange dashed = ERA5.",
                  x = "Year", y = "PET (mm/month)") + theme_cmp
  .save_gg(file.path(out_dir, "pet_comparison_14a_timeseries"), p14a, 12, 5)
  cat("  14a saved\n")
  
  cat("  14b: per-month slope comparison...\n")
  .fit <- function(x, m) {
    sub <- data.frame(yr = yr[mo == m], v = x[mo == m]); sub <- sub[!is.na(sub$v), ]
    if (nrow(sub) < 10) return(c(NA_real_, NA_real_))
    f <- tryCatch(lm(v ~ yr, data = sub), error = function(e) NULL)
    if (is.null(f)) return(c(NA_real_, NA_real_))
    cs <- summary(f)$coefficients
    c(unname(cs["yr", 1]), unname(cs["yr", 4]))
  }
  slopes_b <- do.call(rbind, lapply(1:12, function(m) {
    spm <- .fit(pet_pm, m); sth <- .fit(pet_era5, m)
    data.frame(month = m, month_abb = month.abb[m],
               slope_PM = spm[1], pval_PM = spm[2], sig_PM = !is.na(spm[2]) & spm[2] < 0.05,
               slope_ERA5 = sth[1], pval_ERA5 = sth[2], sig_ERA5 = !is.na(sth[2]) & sth[2] < 0.05)
  }))
  slopes_b$month_abb <- factor(slopes_b$month_abb, levels = month.abb)
  slopes_long <- rbind(
    data.frame(month_abb = slopes_b$month_abb, slope = slopes_b$slope_PM,
               sig = slopes_b$sig_PM, method_lbl = "Penman-Monteith"),
    data.frame(month_abb = slopes_b$month_abb, slope = slopes_b$slope_ERA5,
               sig = slopes_b$sig_ERA5, method_lbl = "ERA5-Land PEV"))
  p14b <- ggplot2::ggplot(slopes_long, ggplot2::aes(x = month_abb, y = slope, fill = method_lbl, alpha = sig)) +
    ggplot2::geom_col(position = "dodge", colour = "white", linewidth = 0.2, width = 0.75) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = c("Penman-Monteith" = "#d73027", "ERA5-Land PEV" = "#e08214"), name = NULL) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.35), guide = "none") +
    ggplot2::labs(title = "Per-Month OLS Slope: PET_PM vs PET_ERA5",
                  subtitle = "Solid bars = p < 0.05; faded = non-significant (mm/month per year).",
                  x = "Calendar month", y = "OLS slope (mm/month/yr)") +
    theme_cmp + ggplot2::theme(legend.position = "top")
  .save_gg(file.path(out_dir, "pet_comparison_14b_monthly_slopes"), p14b, 11, 5)
  cat("  14b saved\n")
  
  cat("  14c: basin trend maps (TFPW Sen's slope, PM vs ERA5)...\n")
  # VC-MK is essentially fully filtered for PET, so the spatial comparison uses TFPW.
  plot_pet_trend_panel <- function(r, title_text) {
    if (is.null(r)) { plot.new(); text(0.5, 0.5, "No valid trend raster", cex = 1.1, col = "grey40", font = 3); return(invisible(NULL)) }
    plot_raster_panel(r, title_text, zlim = .sym_lim(r, 1),
                      col = grDevices::hcl.colors(101, "RdBu", rev = TRUE), legend_title = "mm/yr")
  }
  res14c <- .pdf_safe(file.path(out_dir, "pet_comparison_14c_trendsmap.pdf"), width = 10, height = 5,
                      expr_fn = function() {
                        par(mfrow = c(1, 2), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
                        plot_pet_trend_panel(if (!is.null(r_pet)) r_pet$mag_tfpw else NULL, "PET-PM: TFPW Sen's slope (mm/yr)")
                        plot_pet_trend_panel(if (!is.null(r_pet_era5)) r_pet_era5$mag_tfpw else NULL, "PET-ERA5 PEV: TFPW Sen's slope (mm/yr)")
                        mtext("PET trend comparison - TFPW Sen's slope | Penman-Monteith vs ERA5-Land PEV",
                              outer = TRUE, cex = 1.0, font = 2, line = 0.5)
                      })
  if (!is.null(res14c)) cat("  14c saved\n")
  
  cat("  14d: monthly PET bias (PM - ERA5)...\n")
  bias_df <- do.call(rbind, lapply(1:12, function(m) {
    mu_pm <- mean(pet_pm[mo == m], na.rm = TRUE); mu_e <- mean(pet_era5[mo == m], na.rm = TRUE)
    data.frame(month = m, month_abb = month.abb[m], bias = mu_pm - mu_e)
  }))
  bias_df$month_abb <- factor(bias_df$month_abb, levels = month.abb)
  bias_df$sign_col  <- ifelse(bias_df$bias > 0, "PM > ERA5", "ERA5 > PM")
  p14d <- ggplot2::ggplot(bias_df, ggplot2::aes(x = month_abb, y = bias, fill = sign_col)) +
    ggplot2::geom_col(colour = "white", linewidth = 0.25) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = c("PM > ERA5" = "#d73027", "ERA5 > PM" = "#e08214"), name = NULL) +
    ggplot2::labs(title = "Monthly Mean PET Bias: Penman-Monteith minus ERA5-Land PEV",
                  subtitle = "Positive = PM exceeds ERA5.", x = "Calendar month", y = "PET bias (mm/month)") +
    theme_cmp + ggplot2::theme(legend.position = "top")
  .save_gg(file.path(out_dir, "pet_comparison_14d_monthly_bias"), p14d, 10, 5)
  cat("  14d saved\n")
  
  cat("  14e: per-month scatter PM vs ERA5...\n")
  sc_df <- data.frame(date = dates, pm = pet_pm, era5 = pet_era5,
                      month_abb = factor(month.abb[mo], levels = month.abb),
                      decade = factor(paste0(floor(yr / 10) * 10, "s")))
  sc_df <- sc_df[!is.na(sc_df$pm) & !is.na(sc_df$era5), ]
  p14e <- ggplot2::ggplot(sc_df, ggplot2::aes(x = era5, y = pm, colour = decade)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(size = 0.8, alpha = 0.5) +
    ggplot2::facet_wrap(~month_abb, ncol = 4, scales = "free") +
    ggplot2::scale_colour_brewer(palette = "RdYlBu", direction = -1, name = "Decade") +
    ggplot2::labs(title = "Monthly PET Scatter: Penman-Monteith (y) vs ERA5-Land PEV (x)",
                  subtitle = "Dashed 1:1 line. Colour = decade.",
                  x = "PET ERA5-Land PEV (mm/month)", y = "PET Penman-Monteith (mm/month)") +
    theme_cmp + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", size = 8),
                               legend.position = "right")
  .save_gg(file.path(out_dir, "pet_comparison_14e_scatter"), p14e, 13, 9)
  cat("  14e saved\n-- Function 14 complete.\n")
  invisible(NULL)
}

####################################################################################
# SEASONAL TRENDS - shared helpers (Functions 15, 16, 17)
# Incomplete seasons are dropped (the first DJF only has Jan-Feb; a x3 total from
# 2 months was biased low).
####################################################################################
.SEASON_LEVELS  <- c("DJF (Dec\u2013Feb)", "MAM (Mar\u2013May)", "JJA (Jun\u2013Aug)", "SON (Sep\u2013Nov)")
.SEASON_COLOURS <- stats::setNames(c("#4472C4", "#70AD47", "#C00000", "#ED7D31"), .SEASON_LEVELS)

# total = TRUE -> seasonal total (mean x 3, for mm/month); FALSE -> seasonal mean (temperature)
.seas_means_generic <- function(df, value_col, year_start, year_end, total = TRUE) {
  dt  <- data.table::as.data.table(df)
  raw <- data.table::data.table(date = as.Date(dt[["date"]]), val = as.numeric(dt[[value_col]]))
  raw <- raw[!is.na(val)]
  raw[, `:=`(year = as.integer(format(date, "%Y")), month = as.integer(format(date, "%m")))]
  raw <- raw[year >= year_start & year <= year_end]
  raw[, season := data.table::fcase(
    month %in% c(12L, 1L, 2L),  .SEASON_LEVELS[1],
    month %in% c(3L, 4L, 5L),   .SEASON_LEVELS[2],
    month %in% c(6L, 7L, 8L),   .SEASON_LEVELS[3],
    month %in% c(9L, 10L, 11L), .SEASON_LEVELS[4])]
  raw[, season_year := data.table::fifelse(month == 12L, year + 1L, year)]  # Dec -> next DJF
  mult <- if (total) 3 else 1
  sm <- raw[, .(seas_mean = mean(val) * mult, n_months = .N), by = .(season_year, season)]
  sm <- sm[n_months == 3L & season_year >= year_start & season_year <= year_end]
  sm[, season := factor(season, levels = .SEASON_LEVELS)]
  sm[]
}

.seas_trend_label <- function(sm, lv, unit) {
  parts <- vapply(lv, function(s) {
    d <- sm[as.character(season) == s]
    short <- sub(" \\(.*", "", s)
    if (nrow(d) < 5L) return(sprintf("%s: n/a", short))
    cs <- summary(lm(seas_mean ~ season_year, data = d))$coefficients
    p  <- cs["season_year", 4]
    sprintf("%s: %+.4f %s%s", short, cs["season_year", 1], unit,
            if (is.finite(p) && p < 0.05) " *" else "")
  }, character(1))
  paste(parts, collapse = "   |   ")
}

.seas_ggplot <- function(sm, lv, title_str, subtitle_str, ylab_str, year_start, year_end,
                         base_size = 16) {
  ggplot2::ggplot(sm, ggplot2::aes(x = season_year, y = seas_mean, colour = season, group = season)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
    ggplot2::geom_point(size = 1.2, alpha = 0.75) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed", linewidth = 1.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50", linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = .SEASON_COLOURS[lv], name = NULL,
                                 guide = ggplot2::guide_legend(direction = "horizontal", nrow = 1,
                                                               override.aes = list(linewidth = 1.4, size = 2))) +
    ggplot2::scale_x_continuous(breaks = seq(year_start, year_end, by = 10),
                                expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::labs(title = title_str, subtitle = subtitle_str, x = "Year", y = ylab_str) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = base_size - 2, colour = "grey30", margin = ggplot2::margin(b = 8)),
      legend.position = "top", legend.text = ggplot2::element_text(size = base_size - 1),
      legend.key.width = ggplot2::unit(1.6, "cm"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1))
}

# Generic seasonal-trend figure for one basin-average column
.make_seasonal_figure <- function(col, what, unit_lab, slope_unit, stem, total,
                                  year_start = 1950L, year_end = 2025L,
                                  basin_label = "Nechako River Basin", base_size = 16) {
  if (is.null(basin_avg_monthly) || !all(c("date", col) %in% names(basin_avg_monthly))) {
    cat(sprintf("  %s: column '%s' not found - skipped\n", what, col)); return(invisible(NULL))
  }
  res <- .seas_means_generic(basin_avg_monthly, col, year_start, year_end, total = total)
  if (!nrow(res)) { cat(sprintf("  %s: no complete seasons - skipped\n", what)); return(invisible(NULL)) }
  p <- .seas_ggplot(res, .SEASON_LEVELS,
                    sprintf("Seasonal %s Trends \u2014 %s (%d\u2013%d)", what, basin_label, year_start, year_end),
                    paste0("Annual seasonal ", if (total) "totals" else "means",
                           " with OLS trend lines (dashed).  * p < 0.05\n",
                           .seas_trend_label(res, .SEASON_LEVELS, slope_unit)),
                    unit_lab, year_start, year_end, base_size)
  saved <- .save_gg(file.path(out_dir, stem), p, 14, 7)
  cat(sprintf("  %s seasonal trends: %s\n", what, paste(saved, collapse = " / ")))
  invisible(list(plot = p, seas = res))
}

####################################################################################
# FUNCTION 15: PET seasonal trends (PM, ERA5, ERA5 - PM difference)
####################################################################################
create_pet_seasonal_trends <- function(year_start = 1950L, year_end = 2025L,
                                       basin_label = "Nechako River Basin") {
  cat("\n-- Function 15: PET seasonal trends (PM, ERA5, difference) --\n")
  r_pm <- .make_seasonal_figure("pet_mm_month", "PET-PM", "Seasonal total PET-PM (mm/season)",
                                "mm/yr", "FigS_PET_Seasonal_Trends", TRUE, year_start, year_end,
                                basin_label, 16)
  r_e5 <- .make_seasonal_figure("pet_era5_mm_month", "PET-ERA5", "Seasonal total PET-ERA5 (mm/season)",
                                "mm/yr", "FigS_PET_ERA5_Seasonal_Trends", TRUE, year_start, year_end,
                                basin_label, 16)
  if (!is.null(r_pm) && !is.null(r_e5)) {
    cat("  15c: PET difference (ERA5 - PM), MAM/JJA/SON...\n")
    a <- data.table::copy(r_pm$seas)[, .(season_year, season, pm = seas_mean)]
    b <- data.table::copy(r_e5$seas)[, .(season_year, season, era5 = seas_mean)]
    sm <- merge(a, b, by = c("season_year", "season"))
    sm[, seas_mean := era5 - pm]
    lv <- .SEASON_LEVELS[-1]                              # drop winter (DJF)
    sm <- sm[as.character(season) %in% lv]
    sm[, season := factor(as.character(season), levels = lv)]
    p <- .seas_ggplot(sm, lv,
                      sprintf("Seasonal PET Difference (ERA5-Land PEV \u2212 Penman\u2013Monteith) \u2014 %s (%d\u2013%d)",
                              basin_label, year_start, year_end),
                      paste0("Seasonal totals of PET-ERA5 minus PET-PM with OLS trend lines (dashed).  * p < 0.05\n",
                             .seas_trend_label(sm, lv, "mm/yr")),
                      "PET-ERA5 minus PET-PM (mm/season)", year_start, year_end, 16)
    saved <- .save_gg(file.path(out_dir, "FigS_PET_Both_Seasonal_Trends"), p, 14, 7)
    cat(sprintf("  15c: %s\n", paste(saved, collapse = " / ")))
  }
  cat("-- Function 15 complete.\n")
  invisible(NULL)
}

####################################################################################
# FUNCTION 16: Precipitation seasonal trends
####################################################################################
create_precip_seasonal_trends <- function() {
  cat("\n-- Function 16: Precipitation seasonal trends --\n")
  r <- .make_seasonal_figure("precip_mm_month", "Precipitation", "Seasonal total precipitation (mm/season)",
                             "mm/yr", "FigS_Precip_Seasonal_Trends", TRUE, base_size = 22)
  if (!is.null(r)) {
    chk <- r$seas[, .(n = .N, mean_mm = round(mean(seas_mean), 1)), by = season]
    print(as.data.frame(chk))
  }
  invisible(r)
}

####################################################################################
# FUNCTION 17: Temperature seasonal trends (seasonal MEAN, not total)
####################################################################################
create_temp_seasonal_trends <- function() {
  cat("\n-- Function 17: Temperature seasonal trends --\n")
  r <- .make_seasonal_figure("tair_degC_month", "Temperature", paste0("Seasonal mean temperature (", DEGC, ")"),
                             paste0(DEGC, "/yr"), "FigS_Temp_Seasonal_Trends", FALSE, base_size = 22)
  if (!is.null(r)) {
    chk <- r$seas[, .(n = .N, mean_degC = round(mean(seas_mean), 2)), by = season]
    print(as.data.frame(chk))
  }
  invisible(r)
}

####################################################################################
# FUNCTION: Seasonal climatology P vs PET (bar chart)
####################################################################################
create_seasonality_p_pet_plot <- function() {
  cat("\n-- Seasonality P vs PET plot --\n")
  df <- data.frame(date = basin_avg_monthly$date, P = basin_avg_monthly$precip_mm_month,
                   PET_PM = basin_avg_monthly$pet_mm_month)
  df$Month <- factor(month.abb[as.integer(format(df$date, "%m"))], levels = month.abb)
  n_years <- length(unique(as.integer(format(df$date, "%Y"))))
  yr_rng  <- range(as.integer(format(df$date, "%Y")))
  
  stats_p   <- aggregate(P ~ Month, data = df, FUN = function(x) c(mean = mean(x), sd = sd(x)))
  stats_pet <- aggregate(PET_PM ~ Month, data = df, FUN = function(x) c(mean = mean(x), sd = sd(x)))
  clim_df <- data.frame(Month = stats_p$Month, P_Mean = stats_p$P[, "mean"], P_SD = stats_p$P[, "sd"],
                        PET_Mean = stats_pet$PET_PM[, "mean"], PET_SD = stats_pet$PET_PM[, "sd"])
  clim_df$Month   <- factor(clim_df$Month, levels = month.abb)
  clim_df$Balance <- clim_df$P_Mean - clim_df$PET_Mean
  ann_P <- sum(clim_df$P_Mean); ann_PET <- sum(clim_df$PET_Mean)
  subtitle_txt <- sprintf("Each bar = mean of %d years (%d\u2013%d). Annual totals: P = %.0f mm; PET(PM) = %.0f mm; aridity index (PET/P) = %.2f.",
                          n_years, yr_rng[1], yr_rng[2], ann_P, ann_PET, ann_PET / ann_P)
  
  plot_data_a <- data.frame(Month = rep(clim_df$Month, 2),
                            Variable = factor(rep(c("Precipitation (P)", "PET(PM)"), each = nrow(clim_df)),
                                              levels = c("Precipitation (P)", "PET(PM)")),
                            Mean = c(clim_df$P_Mean, clim_df$PET_Mean),
                            SD = c(clim_df$P_SD, clim_df$PET_SD))
  p_a <- ggplot(plot_data_a, aes(x = Month, y = Mean, fill = Variable)) +
    geom_col(position = position_dodge(0.8), width = 0.7, colour = "white", linewidth = 0.2) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(0.8),
                  width = 0.2, linewidth = 0.6) +
    scale_fill_manual(values = c("Precipitation (P)" = "#4575b4", "PET(PM)" = "#d73027")) +
    labs(y = "mm/month", title = "(a) Mean monthly precipitation and PET(PM)",
         subtitle = sprintf("%d-year mean (%d\u2013%d) | Error bars = \u00b11 SD", n_years, yr_rng[1], yr_rng[2])) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
          plot.subtitle = element_text(size = 8, colour = "grey40", hjust = 0),
          axis.title.x = element_blank(), legend.position = "top", legend.justification = "right",
          panel.grid.major.y = element_line(colour = "grey90"))
  
  clim_df$BalanceCol <- ifelse(clim_df$Balance > 0, "#4575b4", "#d73027")
  p_b <- ggplot(clim_df, aes(x = Month, y = Balance, fill = BalanceCol)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2, show.legend = FALSE) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey50") +
    scale_fill_identity() +
    labs(x = "Calendar month", y = "P \u2212 PET(PM) (mm/month)",
         title = "(b) Climatic water balance P \u2212 PET(PM)",
         subtitle = "Blue = monthly surplus (P > PET); red = monthly deficit (PET > P)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0),
          plot.subtitle = element_text(size = 8, colour = "grey40", hjust = 0),
          panel.grid.major.y = element_line(colour = "grey90"))
  if (any(clim_df$Balance > 0, na.rm = TRUE))
    p_b <- p_b + annotate("text", x = 6.5, y = max(clim_df$Balance, na.rm = TRUE) * 0.8,
                          label = "P surplus", colour = "#4575b4", size = 4, fontface = "italic")
  
  combined <- (p_a / p_b) +
    plot_annotation(title = "Nechako River Basin \u2014 Seasonal cycle of precipitation and PET(PM)",
                    subtitle = subtitle_txt,
                    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
                                  plot.subtitle = element_text(size = 10, colour = "grey40", hjust = 0)))
  saved <- .save_gg(file.path(out_dir, "FigS_Seasonality_P_PET"), combined, 10, 9, formats = "png")
  cat("  Saved:", paste(saved, collapse = ", "), "\n")
  invisible(combined)
}

####################################################################################
# FUNCTIONS 9-11: Specific-point plots
####################################################################################
create_point_timeseries_plots <- function() {
  if (is.null(point_monthly_ts)) { cat("  point_monthly_timeseries.csv not available - skipping.\n"); return(NULL) }
  cat("Creating specific-point time series plots...\n")
  pm_ts <- data.table::copy(point_monthly_ts)
  pm_ts[, date := as.Date(date)]
  for (pid in unique(pm_ts$point_id)) {
    sub <- pm_ts[point_id == pid]
    lon <- round(sub$lon_wgs84[1], 4); lat <- round(sub$lat_wgs84[1], 4)
    title_str <- sprintf("Point #%d (%.4f\u00b0N, %.4f\u00b0W)", pid, lat, abs(lon))
    make_var_panel <- function(var_name, y_label, color_line, color_smooth) {
      df <- data.table::copy(sub[variable == var_name])
      if (!nrow(df)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste(var_name, "not available")) + theme_void())
      df <- df[order(date)]
      df[, smooth := zoo::rollmean(value, k = min(12, .N), fill = NA, align = "center")]
      ggplot(df, aes(x = date, y = value)) +
        geom_line(colour = color_line, alpha = 0.45, linewidth = 0.5) +
        geom_line(aes(y = smooth), colour = color_smooth, linewidth = 1.1, na.rm = TRUE) +
        labs(title = var_name, x = NULL, y = y_label) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 30, hjust = 1))
    }
    panels <- list(make_var_panel("Precipitation", "mm/month", "steelblue", "blue"),
                   make_var_panel("PET", "mm/month", "tomato", "red3"),
                   make_var_panel("PET_ERA5", "mm/month", "orange", "#e08214"),
                   make_var_panel("Temperature", DEGC, "orchid", "purple4"))
    if (any(sub$variable == "SWE"))
      panels[[length(panels) + 1L]] <- make_var_panel("SWE", "mm w.e.", "#80cdc1", "#01665e")
    for (nm in EXTRA_VARS_PRESENT) if (any(sub$variable == nm))
      panels[[length(panels) + 1L]] <- make_var_panel(nm, EXTRA_MONTH_UNITS[[nm]],
                                                      grDevices::adjustcolor(EXTRA_COLOURS[[nm]], 0.6), EXTRA_COLOURS[[nm]])
    combined <- patchwork::wrap_plots(panels, ncol = 1) +
      plot_annotation(title = title_str, subtitle = "ERA5-Land monthly | coloured line = 12-month rolling mean",
                      theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                                    plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5)))
    stem <- file.path(out_dir, sprintf("point_%d_%.4fN_%.4fW_ts_pr_pet_temp", pid, lat, abs(lon)))
    saved <- .save_gg(stem, combined, 12, 3.25 * length(panels), dpi = 150, formats = "png")
    cat("  Saved:", paste(saved, collapse = ", "), "\n")
  }
  invisible(NULL)
}

create_point_monthly_trend_plots <- function() {
  if (is.null(point_trend_stats)) { cat("  point_trend_stats.csv not available - skipping.\n"); return(NULL) }
  cat("Creating specific-point monthly trend bar charts...\n")
  var_units <- c(Precipitation = "Sen slope (mm/month/yr)", PET = "Sen slope (mm/month/yr)",
                 PET_ERA5 = "Sen slope (mm/month/yr)", Temperature = paste0("Sen slope (", DEGC, "/yr)"),
                 SWE = "Sen slope (mm w.e./yr)")
  var_colors <- c(Precipitation = "#4575b4", PET = "#d73027", PET_ERA5 = "#e08214", Temperature = "#7b2d8b", SWE = "#01665e")
  var_units  <- c(var_units, setNames(paste0("Sen slope (", EXTRA_MONTH_UNITS, "/yr)"), names(EXTRA_MONTH_UNITS)))
  var_colors <- c(var_colors, EXTRA_COLOURS)
  for (pid in unique(point_trend_stats$point_id)) {
    sub <- data.table::copy(point_trend_stats[point_id == pid & period == "monthly"])
    if (!nrow(sub)) next
    lon <- round(sub$lon_wgs84[1], 4); lat <- round(sub$lat_wgs84[1], 4)
    sub[, month_abb := factor(month_abb, levels = month.abb)]
    sub[, sig_label := ifelse(!is.na(sig_vc) & sig_vc, "*", "")]
    sub[, label_y := sl_vc + ifelse(!is.na(sl_vc) & sl_vc >= 0, 1, -1) *
          0.03 * max(abs(sl_vc), na.rm = TRUE), by = variable]      # offset scaled per variable
    pt_vars <- c("Precipitation", "PET", "PET_ERA5", "Temperature", if (any(sub$variable == "SWE")) "SWE",
                 intersect(EXTRA_VARS_PRESENT, unique(sub$variable)))
    panels <- lapply(pt_vars, function(v) {
      d <- sub[variable == v]
      if (!nrow(d)) return(ggplot() + theme_void() + labs(title = paste(v, "- no data")))
      ggplot(d, aes(x = month_abb, y = sl_vc)) +
        geom_bar(stat = "identity", colour = "black", linewidth = 0.3, width = 0.7, fill = var_colors[[v]], na.rm = TRUE) +
        geom_hline(yintercept = 0, linewidth = 0.5) +
        geom_text(aes(y = label_y, label = sig_label), size = 4.5, fontface = "bold", na.rm = TRUE) +
        scale_x_discrete(labels = substr(month.abb, 1, 1)) +
        labs(title = v, x = "Month", y = var_units[[v]]) +
        theme_classic(base_size = 10) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
    })
    combined <- patchwork::wrap_plots(panels, nrow = 1) +
      plot_annotation(title = sprintf("Per-Month Sen's Slope \u2014 Point #%d (%.4f\u00b0N, %.4f\u00b0W)", pid, lat, abs(lon)),
                      subtitle = "* = significant by VC Mann-Kendall (p < 0.05)",
                      theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
                                    plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5)))
    stem <- file.path(out_dir, sprintf("point_%d_%.4fN_%.4fW_monthly_slopes", pid, lat, abs(lon)))
    saved <- .save_gg(stem, combined, 4 * length(panels), 5, dpi = 150, formats = "png")
    cat("  Saved:", paste(saved, collapse = ", "), "\n")
  }
  invisible(NULL)
}

create_point_summary_table <- function() {
  if (is.null(point_trend_stats)) { cat("  point_trend_stats.csv not available.\n"); return(NULL) }
  cat("\n====================================================\n")
  cat(" SPECIFIC POINT TREND STATISTICS SUMMARY\n")
  cat("====================================================\n")
  annual <- point_trend_stats[period == "annual"]
  if (nrow(annual)) {
    cat("Annual trends:\n")
    cols <- intersect(c("point_id", "lat_wgs84", "lon_wgs84", "variable", "n_obs", "year_min", "year_max",
                        "tau_vc", "p_value_vc", "sig_vc", "sl_vc", "tau_tfpw", "p_value_tfpw",
                        "sig_tfpw", "sl_tfpw", "first_changepoint_year"), names(annual))
    tab <- as.data.frame(annual[, ..cols])
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], round, 4)
    print(tab)
  }
  sig_monthly <- point_trend_stats[period == "monthly" & !is.na(sig_vc) & sig_vc == TRUE]
  if (nrow(sig_monthly)) {
    cat("\nSignificant monthly trends (VC-MK p < 0.05):\n")
    cols <- intersect(c("point_id", "variable", "month_abb", "tau_vc", "p_value_vc", "sl_vc", "sig_tfpw"),
                      names(sig_monthly))
    tab <- as.data.frame(sig_monthly[, ..cols])
    num <- vapply(tab, is.numeric, logical(1))
    tab[num] <- lapply(tab[num], round, 5)
    print(tab)
  } else cat("\n (No significant monthly trends at any specific point)\n")
  invisible(annual)
}

####################################################################################
# FUNCTION 12: Temperature dedicated analysis
####################################################################################
create_temperature_dedicated_plots <- function() {
  if (!"tair_degC_month" %in% names(basin_avg_monthly)) {
    cat("  basin_avg_monthly$tair_degC_month not found - skipping Function 12.\n")
    return(invisible(NULL))
  }
  tair <- basin_avg_monthly$tair_degC_month
  dates <- basin_avg_monthly$date
  yr <- as.integer(format(dates, "%Y")); mo <- as.integer(format(dates, "%m"))
  yr_frac <- yr + (mo - 0.5) / 12
  .ols <- function(x, y) { keep <- !is.na(x) & !is.na(y); if (sum(keep) < 5) return(NULL); lm(y[keep] ~ x[keep]) }
  theme_tair <- ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey40", hjust = 0),
                   panel.grid.major = ggplot2::element_line(colour = "grey94", linewidth = 0.3),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(4, 8, 4, 4))
  cat("\n-- Function 12: Temperature dedicated plots --\n")
  tair_roll <- zoo::rollmean(tair, k = 12, fill = NA, align = "center")
  fit_all <- .ols(yr_frac, tair)
  pred_df <- NULL; ann_slope <- "trend not available"
  if (!is.null(fit_all)) {
    ok <- !is.na(tair)
    pred_df <- data.frame(date = dates[ok], tfit = coef(fit_all)[1] + coef(fit_all)[2] * yr_frac[ok])
    ann_slope <- sprintf("%+.4f %s/yr (%+.3f %s/decade)", coef(fit_all)[2], DEGC, coef(fit_all)[2] * 10, DEGC)
  }
  ann_df <- NULL
  if (!is.null(basin_avg_annual) && "tair_degC_year" %in% names(basin_avg_annual))
    ann_df <- data.frame(date = as.Date(paste0(basin_avg_annual$year, "-07-01")), tair = basin_avg_annual$tair_degC_year)
  df_ts <- data.frame(date = dates, tair = tair, troll = tair_roll)
  p12a <- ggplot2::ggplot(df_ts, ggplot2::aes(x = date)) +
    ggplot2::geom_line(ggplot2::aes(y = tair), colour = rgb(0.6, 0.2, 0.6, 0.3), linewidth = 0.45, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = troll), colour = "purple", linewidth = 1.6, na.rm = TRUE)
  if (!is.null(ann_df))
    p12a <- p12a + ggplot2::geom_point(data = ann_df, ggplot2::aes(x = date, y = tair),
                                       colour = rgb(0.4, 0, 0.4, 0.8), size = 1.4, shape = 19)
  if (!is.null(pred_df))
    p12a <- p12a + ggplot2::geom_line(data = pred_df, ggplot2::aes(x = date, y = tfit),
                                      colour = "#d73027", linewidth = 1.2, linetype = "dashed")
  p12a <- p12a +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dotted", linewidth = 0.5) +
    ggplot2::scale_x_date(date_breaks = "10 years", date_labels = "%Y", expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::labs(title = "Basin-Mean 2-m Air Temperature \u2014 Nechako River Basin (1950\u20132025)",
                  subtitle = paste0("Monthly (faded), 12-mo mean (purple), annual mean (dots), OLS (red dashed)\nOLS slope: ", ann_slope),
                  x = "Year", y = paste0("Temperature (", DEGC, ")")) + theme_tair
  .save_gg(file.path(out_dir, "temperature_basin_timeseries"), p12a, 10, 5)
  cat("  12a saved\n")
  
  ann_tair <- tapply(tair, yr, mean, na.rm = TRUE)
  ann_years <- as.integer(names(ann_tair)); ann_tair <- as.numeric(ann_tair)
  clim_idx <- ann_years >= 1991 & ann_years <= 2020
  clim_mean <- if (sum(clim_idx) > 0) mean(ann_tair[clim_idx], na.rm = TRUE) else mean(ann_tair, na.rm = TRUE)
  df_anom <- data.frame(year = ann_years, anomaly = ann_tair - clim_mean)
  df_anom$sign <- ifelse(df_anom$anomaly >= 0, "Warm", "Cool")
  fit_anom <- .ols(df_anom$year, df_anom$anomaly)
  anom_str <- if (!is.null(fit_anom)) sprintf("OLS: %+.4f %s/yr", coef(fit_anom)[2], DEGC) else ""
  p12b <- ggplot2::ggplot(df_anom, ggplot2::aes(x = year, y = anomaly, fill = sign)) +
    ggplot2::geom_bar(stat = "identity", width = 0.85, colour = "white", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = c(Warm = "#d73027", Cool = "#4575b4"), name = NULL) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50", linetype = "dotted", linewidth = 0.5) +
    ggplot2::scale_x_continuous(breaks = seq(1950, 2025, by = 10)) +
    ggplot2::labs(title = "Annual Temperature Anomaly \u2014 Nechako River Basin",
                  subtitle = sprintf("Baseline: WMO 1991\u20132020 (%.2f%s). %s", clim_mean, DEGC, anom_str),
                  x = "Year", y = paste0("Temperature anomaly (", DEGC, ")")) + theme_tair +
    ggplot2::theme(legend.position = "top")
  .save_gg(file.path(out_dir, "temperature_annual_anomaly"), p12b, 10, 5)
  cat("  12b saved\n-- Function 12 complete (12a, 12b).\n")
  invisible(NULL)
}

####################################################################################
# REGIONAL EXTENSION - climate / PET figures for the sub-basins + whole basin
####################################################################################
create_regional_climate_plots <- function() {
  if (is.null(metadata$regional_avg_monthly)) {
    cat("  No regional_avg_monthly in metadata - skipping regional plots.\n"); return(invisible(NULL))
  }
  reg_clim <- data.table::as.data.table(metadata$regional_avg_monthly)
  need <- c("date", "region_id", "region_name", "precip_mm_month", "pet_mm_month",
            "pet_era5_mm_month", "tair_degC_month")
  miss <- setdiff(need, names(reg_clim))
  if (length(miss)) { cat("  Regional table missing columns:", paste(miss, collapse = ", "), "\n"); return(invisible(NULL)) }
  reg_clim[, date := as.Date(date)]
  reg_dir <- file.path(out_dir, "regional_plots")
  dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (rid in unique(reg_clim$region_id)) {
    d <- reg_clim[region_id == rid][order(date)]
    reg_cols <- c("precip_mm_month", "pet_mm_month", "pet_era5_mm_month", "tair_degC_month",
                  if ("swe_mm_month" %in% names(d) && any(is.finite(d$swe_mm_month))) "swe_mm_month")
    reg_labs <- c("Precipitation", "PET-PM", "PET-ERA5", "Temperature", if ("swe_mm_month" %in% reg_cols) "SWE")
    ex <- EXTRA_MONTH_COLS[EXTRA_MONTH_COLS %in% names(d)]
    ex <- ex[vapply(ex, function(cl) any(is.finite(d[[cl]])), logical(1))]
    reg_cols <- c(reg_cols, unname(ex))
    reg_labs <- c(reg_labs, vapply(names(ex), .short_lab, character(1)))
    dlong <- data.table::melt(d[, c("date", reg_cols), with = FALSE],
                              id.vars = "date", variable.name = "var", value.name = "value")
    dlong[, var := factor(var, levels = reg_cols, labels = reg_labs)]
    p <- ggplot2::ggplot(dlong, ggplot2::aes(date, value)) +
      ggplot2::geom_line(linewidth = 0.35, na.rm = TRUE) +
      ggplot2::facet_wrap(~var, ncol = 1, scales = "free_y") +
      ggplot2::labs(title = unique(d$region_name)[1],
                    subtitle = "Regional area-weighted climate series on the native ERA5-Land grid",
                    x = NULL, y = NULL) + ggplot2::theme_bw(base_size = 10)
    .save_gg(file.path(reg_dir, sprintf("%s_climate_timeseries", rid)), p, 12, 2 * length(reg_cols), formats = "png")
  }
  
  if (!is.null(metadata$regional_avg_annual)) {
    reg_ann <- data.table::as.data.table(metadata$regional_avg_annual)
    if (nrow(reg_ann) && all(c("year", "pet_mm_year_full", "region_name") %in% names(reg_ann))) {
      pdat <- reg_ann[is.finite(pet_mm_year_full)]
      p <- ggplot2::ggplot(pdat, ggplot2::aes(year, pet_mm_year_full)) +
        ggplot2::geom_line(linewidth = 0.35) +
        ggplot2::facet_wrap(~region_name, ncol = 3, scales = "free_y") +
        ggplot2::labs(title = "Annual full-year PET-PM by Nechako region",
                      subtitle = "Calendar-year totals (the same Jan-Dec window as the PET trend metric)",
                      x = "Year", y = "PET-PM (mm/year)") + ggplot2::theme_bw(base_size = 9)
      .save_gg(file.path(reg_dir, "regional_annual_PET_fullyear"), p, 14, 12, formats = "png")
    }
  }
  cat("  Regional climate/PET visualization products written.\n")
  invisible(NULL)
}

####################################################################################
# FUNCTION 0c: SWE TREND MAPS (annual Mar-Apr SWE and peak SWE; TFPW Sen slope + BH-FDR)
# SWE is a snow-depth variable: positive slope = more snow (blue-green), negative = less.
####################################################################################
# Symmetric colour limit = 98th percentile of |slope| (a few extreme high-elevation / glacier
# pixels would otherwise stretch the scale and wash out the rest of the basin; values beyond the
# limit are squished to the end colour by .trend_tile_panel).
.robust_lim <- function(x) {
  a <- abs(x[is.finite(x)])
  z <- if (length(a)) stats::quantile(a, 0.98, names = FALSE) else NA_real_
  if (!is.finite(z) || z <= 0) z <- if (length(a) && max(a) > 0) max(a) else 1
  z
}

create_swe_trend_figure <- function() {
  cat("  SWE Sen-slope/FDR map: ")
  vars <- intersect(c("SWE", "SWE_Peak"), .vars_present)
  if (!length(vars)) { cat("SWE not in all_results - skipped\n"); return(invisible(NULL)) }
  labels <- c(SWE = "Mar-Apr mean SWE", SWE_Peak = "Peak monthly SWE (Jan-Jun maximum)")
  pts <- .prep_tfpw_points(vars)
  if (is.null(pts)) { cat("no valid annual TFPW trend rows found\n"); return(invisible(NULL)) }
  dx <- .grid_spacing(pts$plot_x)
  dy <- .grid_spacing(pts$plot_y)
  panels <- lapply(vars, function(v) {
    s <- pts[pts$variable == v, ]
    if (!nrow(s)) return(.blank_panel(paste(labels[[v]], "\nno valid data")))
    z <- .robust_lim(s$sl_tfpw)
    pal <- .trend_palette_for(v)
    .trend_tile_panel(s, labels[[v]], "mm/yr", c(-z, z), dx, dy,
                      palette_fn = pal$pal, hint = pal$hint, bar_len_cm = 5.5)
  })
  fig <- patchwork::wrap_plots(panels, nrow = 1, guides = "keep") &
    ggplot2::theme(legend.position = "bottom")
  fig <- fig + patchwork::plot_annotation(
    title = "Spatial trends in snow water equivalent \u2014 Nechako River Basin (1950\u20132025)",
    subtitle = paste0("Annual TFPW Mann-Kendall Sen's slope (mm w.e. per year); filled black points = BH-FDR-significant under BOTH VC-MK and TFPW-MK, open points = one test only ",
                      "pixels (alpha = 0.05). Colour scale clipped at the 98th percentile of |slope|. Dashed lines = sub-basins."),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0.5, colour = "grey30", margin = ggplot2::margin(b = 5))))
  saved <- .save_gg(file.path(out_dir, "FigSx_Spatial_SWE_Trends_SenSlope"), fig, 5.6 * length(vars), 5.8)
  cat(paste(saved, collapse = " / "), "\n")
  invisible(fig)
}

####################################################################################
# FUNCTION 0d: SWE TREND MAPS BY CALENDAR MONTH (TFPW Sen slope, BH-FDR within month)
# Shows WHEN in the year the pixel-level trends occur. Months with too few snow-covered
# years (Jul-Aug) are filtered upstream and appear as blank panels.
####################################################################################
.prep_tfpw_points_month <- function(var, mm) {
  d <- data.table::copy(all_results[period == "monthly" & !is_basin_average & variable == var &
                                      month == mm & is.finite(sl_tfpw) & is.finite(p_value_tfpw)])
  if (nrow(d) < 5L) return(NULL)
  d <- .add_dual_fdr(d)
  pts <- sf::st_as_sf(as.data.frame(d), coords = c("x", "y"),
                      crs = terra::crs(template_bc), remove = FALSE)
  pts <- sf::st_transform(pts, 4326)
  xy  <- sf::st_coordinates(pts)
  pts$plot_x <- xy[, 1]
  pts$plot_y <- xy[, 2]
  pts
}

create_swe_monthly_trend_maps <- function() {
  cat("  SWE monthly Sen-slope maps: ")
  if (!"SWE" %in% .vars_present) { cat("SWE not in all_results - skipped\n"); return(invisible(NULL)) }
  px <- unique(all_results[!is_basin_average & is.finite(x) & is.finite(y), .(x, y)])
  dx <- .grid_spacing(px$x)
  dy <- .grid_spacing(px$y)
  pal <- .trend_palette_for("SWE")
  panels <- lapply(1:12, function(mm) {
    s <- .prep_tfpw_points_month("SWE", mm)
    if (is.null(s)) return(.blank_panel(paste0(month.name[mm], "\nfiltered (snow-free / too few non-zero years)")))
    z <- .robust_lim(s$sl_tfpw)
    .trend_tile_panel(s, month.name[mm], "mm/yr", c(-z, z), dx, dy,
                      palette_fn = pal$pal, hint = NULL, bar_len_cm = 3.4)
  })
  fig <- patchwork::wrap_plots(panels, nrow = 3, ncol = 4, guides = "keep") &
    ggplot2::theme(legend.position = "bottom")
  fig <- fig + patchwork::plot_annotation(
    title = "SWE trends by calendar month \u2014 Nechako River Basin (1950\u20132025)",
    subtitle = paste0("TFPW Mann-Kendall Sen's slope (mm w.e. per year) of the monthly-mean SWE; each panel has its own ",
                      "colour scale (clipped at the 98th percentile of |slope|). Filled black points = BH-FDR-significant within that month under BOTH VC-MK and TFPW-MK; open points = one test only."),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0.5, colour = "grey30", margin = ggplot2::margin(b = 5))))
  saved <- .save_gg(file.path(out_dir, "FigSx_Spatial_SWE_Monthly_Trends_SenSlope"), fig, 13, 12.5)
  cat(paste(saved, collapse = " / "), "\n")
  invisible(fig)
}

####################################################################################
# FUNCTION 0e: RUNOFF / SNOWFALL / SNOW DEPTH / SNOW DENSITY - annual TFPW Sen slope + BH-FDR
####################################################################################
create_extra_trend_figure <- function() {
  cat("  Runoff/snow Sen-slope/FDR map: ")
  vars <- EXTRA_VARS_PRESENT
  if (!length(vars)) { cat("no runoff/snow variables in all_results - skipped\n"); return(invisible(NULL)) }
  pts <- .prep_tfpw_points(vars)
  if (is.null(pts)) { cat("no valid annual TFPW trend rows\n"); return(invisible(NULL)) }
  dx <- .grid_spacing(pts$plot_x); dy <- .grid_spacing(pts$plot_y)
  panels <- lapply(seq_along(vars), function(i) {
    v <- vars[[i]]; lab <- sprintf("%s) %s", letters[i], EXTRA_VAR_LABELS[[v]])
    s <- pts[pts$variable == v, ]
    if (!nrow(s)) return(.blank_panel(paste(lab, "\nno valid data")))
    z <- .robust_lim(s$sl_tfpw); pal <- .trend_palette_for(v)
    .trend_tile_panel(s, lab, EXTRA_SLOPE_UNITS[[v]], c(-z, z), dx, dy,
                      palette_fn = pal$pal, hint = pal$hint, bar_len_cm = 5.5)
  })
  nc <- min(2L, length(panels))
  fig <- patchwork::wrap_plots(panels, ncol = nc, guides = "keep") & ggplot2::theme(legend.position = "bottom")
  fig <- fig + patchwork::plot_annotation(
    title = "Spatial trends in evapotranspiration, runoff and snowpack \u2014 Nechako River Basin (1950\u20132025)",
    subtitle = paste0("Annual TFPW Mann-Kendall Sen's slope (AET, runoff, snowfall: calendar-year totals; snow depth, density: ",
                      "Mar-Apr means). Filled black points = BH-FDR-significant under BOTH VC-MK and TFPW-MK, open points = one test only (alpha = 0.05). ",
                      "Colour scale clipped at the 98th percentile of |slope|. Dashed lines = sub-basins."),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
                           plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0.5, colour = "grey30")))
  saved <- .save_gg(file.path(out_dir, "FigSx_Spatial_Runoff_Snow_Trends_SenSlope"), fig,
                    5.6 * nc, 5.8 * ceiling(length(panels) / nc))
  cat(paste(saved, collapse = " / "), "\n")
  invisible(fig)
}

####################################################################################
# FUNCTION 0f: monthly Sen-slope maps for one runoff/snow variable (BH-FDR within month)
####################################################################################
create_var_monthly_trend_maps <- function(v) {
  cat(sprintf("  %s monthly Sen-slope maps: ", v))
  if (!v %in% .vars_present) { cat("not in all_results - skipped\n"); return(invisible(NULL)) }
  px <- unique(all_results[!is_basin_average & is.finite(x) & is.finite(y), .(x, y)])
  dx <- .grid_spacing(px$x); dy <- .grid_spacing(px$y)
  pal <- .trend_palette_for(v)
  panels <- lapply(1:12, function(mm) {
    s <- .prep_tfpw_points_month(v, mm)
    if (is.null(s)) return(.blank_panel(paste0(month.name[mm], "\nfiltered / too few valid years")))
    z <- .robust_lim(s$sl_tfpw)
    .trend_tile_panel(s, month.name[mm], EXTRA_MONTHLY_SLOPE_UNITS[[v]], c(-z, z), dx, dy,
                      palette_fn = pal$pal, hint = NULL, bar_len_cm = 3.4)
  })
  fig <- patchwork::wrap_plots(panels, nrow = 3, ncol = 4, guides = "keep") & ggplot2::theme(legend.position = "bottom")
  fig <- fig + patchwork::plot_annotation(
    title = sprintf("%s trends by calendar month \u2014 Nechako River Basin (1950\u20132025)", .short_lab(v)),
    subtitle = "TFPW Mann-Kendall Sen's slope of the monthly series; own colour scale per panel; filled black points = BH-FDR-significant within that month under BOTH VC-MK and TFPW-MK; open points = one test only.",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 13, face = "bold", hjust = 0.5),
                           plot.subtitle = ggplot2::element_text(size = 8.5, hjust = 0.5, colour = "grey30")))
  saved <- .save_gg(file.path(out_dir, sprintf("FigSx_Spatial_%s_Monthly_Trends_SenSlope", v)), fig, 13, 12.5)
  cat(paste(saved, collapse = " / "), "\n")
  invisible(fig)
}

####################################################################################
# FUNCTION 18: SWE seasonal trends (seasonal MEAN SWE, not total)
####################################################################################
create_swe_seasonal_trends <- function() {
  cat("\n-- Function 18: SWE seasonal trends --\n")
  r <- .make_seasonal_figure("swe_mm_month", "SWE", "Seasonal mean SWE (mm w.e.)",
                             "mm/yr", "FigS_SWE_Seasonal_Trends", FALSE, base_size = 22)
  if (!is.null(r)) {
    chk <- r$seas[, .(n = .N, mean_mm = round(mean(seas_mean), 1)), by = season]
    print(as.data.frame(chk))
  }
  invisible(r)
}

####################################################################################
# FUNCTION 18b: runoff / snow seasonal trends
####################################################################################
create_extra_seasonal_trends <- function() {
  cat("\n-- Function 18b: runoff / snow seasonal trends --\n")
  spec <- list(
    AET         = list(total = TRUE,  ylab = "Seasonal total AET (mm/season)",            su = "mm/yr"),
    Runoff      = list(total = TRUE,  ylab = "Seasonal total runoff (mm/season)",         su = "mm/yr"),
    Snowfall    = list(total = TRUE,  ylab = "Seasonal total snowfall (mm w.e./season)",  su = "mm/yr"),
    SnowDepth   = list(total = FALSE, ylab = "Seasonal mean snow depth (cm)",             su = "cm/yr"),
    SnowDensity = list(total = FALSE, ylab = "Seasonal mean snow density (kg m-3)",       su = "kg m-3/yr"))
  for (nm in names(spec))
    .run_step(paste("seasonal", nm),
              .make_seasonal_figure(EXTRA_MONTH_COLS[[nm]], .short_lab(nm), spec[[nm]]$ylab, spec[[nm]]$su,
                                    sprintf("FigS_%s_Seasonal_Trends", nm), spec[[nm]]$total, base_size = 22))
  invisible(NULL)
}

####################################################################################
# SWE STEP-CHANGE OUTPUTS (Function 19): figures and tables for a later report
# Reads the tables written by 9pr_pet_trends.r (trend_analysis_pr_pet/swe_homogeneity/)
# and writes figures (PNG) and tables (CSV) to swe_homogeneity/report/. No narrative text
# is written. The pipeline's own Pettitt / Welch p-values are cross-checked against an
# independent recomputation, and every number needed for the two headline statements
# (Nov-Mar temperature step; SWE residual after removing P and T) is collected in
# tables/Table00_key_numbers.csv with the file and column it comes from.
####################################################################################
SWE_HOMOG_DIR  <- file.path(out_dir, "swe_homogeneity")
SWE_REPORT_DIR <- file.path(SWE_HOMOG_DIR, "report")
MINUS <- "\u2212"

.swe_basin_id <- function()
  if (exists("NECHAKO_REGION_BASIN_ID", inherits = TRUE)) get("NECHAKO_REGION_BASIN_ID", inherits = TRUE) else "Nechako_Basin"
.swe_step_year <- function() {
  sy <- tryCatch(metadata$swe_config$step_year, error = function(e) NULL)
  if (is.null(sy) || length(sy) != 1L || !is.finite(sy)) 1979L else as.integer(sy)
}

.SWE_SERIES_LAB  <- c(Precip_NovMar = "Nov-Mar precipitation", Temp_NovMar = "Nov-Mar temperature",
                      SWE_MarApr = "Mar-Apr SWE", SWE_MarApr_resid_PT = "Mar-Apr SWE residual (P and T removed)")
.SWE_SERIES_UNIT <- c(Precip_NovMar = "mm", Temp_NovMar = DEGC, SWE_MarApr = "mm w.e.", SWE_MarApr_resid_PT = "mm w.e.")
.SWE_SERIES_COL  <- c(Precip_NovMar = "#4575b4", Temp_NovMar = "#6a3d9a", SWE_MarApr = "#01665e", SWE_MarApr_resid_PT = "#8c510a")

# ---- formatting helpers ----------------------------------------------------------------
.fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
# share of the raw step accounted for (over 100% = residual step has the opposite sign)
.fmt_expl <- function(f) if (!is.finite(f)) "NA" else if (f > 100) ">100%" else sprintf("%.0f%%", f)
.fmt_sgn <- function(x, d = 1L) ifelse(is.finite(x), sprintf(paste0("%+.", d, "f"), x), "NA")

# ---- loading -----------------------------------------------------------------------------
.swe_load <- function() {
  f <- function(n) { p <- file.path(SWE_HOMOG_DIR, n); if (file.exists(p)) data.table::fread(p) else NULL }
  list(annual = f("swe_homogeneity_annual_series.csv"), tests = f("swe_step_tests_by_region.csv"),
       era = f("swe_period_means.csv"), sub = f("swe_subperiod_trends.csv"),
       refs = f("swe_reference_sensitivity.csv"), pix = f("swe_pixel_step_tests.csv"))
}
.swe_get_series <- function(annual, rid, sname) {
  d <- annual[region_id == rid & series == sname]
  d <- d[order(year)]
  d[is.finite(value)]
}

# ---- Pettitt (1979) statistic curve ---------------------------------------------------------
# U_t = 2 * sum_{i<=t} rank_i - t (n + 1); K = max |U_t|; the break lies AFTER index tau.
# p ~ 2 exp(-6 K^2 / (n^3 + n^2)) (the approximation used by trend::pettitt.test).
.pettitt_curve <- function(x) {
  n <- length(x); r <- rank(x); tt <- seq_len(n - 1L)
  U <- 2 * cumsum(r)[tt] - tt * (n + 1)
  K <- max(abs(U)); tau <- which.max(abs(U))
  list(U = U, K = K, tau = tau, n = n, p = min(1, 2 * exp(-6 * K^2 / (n^3 + n^2))))
}

# ---- before/after comparison at a given first year of the later regime -----------------------
.split_stats <- function(x, yrs, cut_year, B = 5000L, level = 0.95) {
  ok <- is.finite(x); x <- x[ok]; yrs <- yrs[ok]
  b <- x[yrs < cut_year]; a <- x[yrs >= cut_year]
  nb <- length(b); na_ <- length(a)
  out <- data.table::data.table(cut_year = cut_year, n_before = nb, n_after = na_,
                                mean_before = NA_real_, mean_after = NA_real_, diff = NA_real_, diff_lo = NA_real_, diff_hi = NA_real_,
                                boot_lo = NA_real_, boot_hi = NA_real_, pct = NA_real_, pct_lo = NA_real_, pct_hi = NA_real_,
                                cohens_d = NA_real_, d_lo = NA_real_, d_hi = NA_real_, welch_p = NA_real_, mannwhitney_p = NA_real_,
                                sd_pooled = NA_real_, mde_80 = NA_real_, mde_80_pct = NA_real_)
  if (nb < 5L || na_ < 5L) return(out)
  tt <- tryCatch(stats::t.test(a, b, conf.level = level), error = function(e) NULL)
  sp <- sqrt(((nb - 1) * stats::var(b) + (na_ - 1) * stats::var(a)) / (nb + na_ - 2))
  d  <- if (is.finite(sp) && sp > 0) (mean(a) - mean(b)) / sp else NA_real_
  se_d <- sqrt((nb + na_) / (nb * na_) + d^2 / (2 * (nb + na_)))
  bt <- replicate(B, mean(sample(a, replace = TRUE)) - mean(sample(b, replace = TRUE)))
  q  <- stats::quantile(bt, c((1 - level) / 2, 1 - (1 - level) / 2), names = FALSE)
  mb <- abs(mean(b))
  mde <- (stats::qt(0.975, nb + na_ - 2) + stats::qt(0.80, nb + na_ - 2)) * sp * sqrt(1 / nb + 1 / na_)
  out[, `:=`(mean_before = mean(b), mean_after = mean(a), diff = mean(a) - mean(b),
             diff_lo = if (is.null(tt)) NA_real_ else tt$conf.int[1],
             diff_hi = if (is.null(tt)) NA_real_ else tt$conf.int[2],
             boot_lo = q[1], boot_hi = q[2],
             pct = if (mb > 0) 100 * (mean(a) - mean(b)) / mb else NA_real_,
             pct_lo = if (mb > 0 && !is.null(tt)) 100 * tt$conf.int[1] / mb else NA_real_,
             pct_hi = if (mb > 0 && !is.null(tt)) 100 * tt$conf.int[2] / mb else NA_real_,
             cohens_d = d, d_lo = d - stats::qnorm(0.975) * se_d, d_hi = d + stats::qnorm(0.975) * se_d,
             welch_p = if (is.null(tt)) NA_real_ else tt$p.value,
             mannwhitney_p = tryCatch(suppressWarnings(stats::wilcox.test(a, b, exact = FALSE)$p.value),
                                      error = function(e) NA_real_),
             sd_pooled = sp, mde_80 = mde, mde_80_pct = if (mb > 0) 100 * mde / mb else NA_real_)]
  out
}

# ---- regressions that decide whether temperature/precipitation explain the SWE step ------------
.swe_model_frame <- function(annual, rid) {
  g <- function(s) annual[region_id == rid & series == s][order(year), .(year, v = value)]
  y <- g("SWE_MarApr"); p <- g("Precip_NovMar"); tm <- g("Temp_NovMar")
  if (!nrow(y) || !nrow(p) || !nrow(tm)) return(NULL)
  df <- Reduce(function(a, b) merge(a, b, by = "year"),
               list(data.table::setnames(data.table::copy(y), "v", "y"), data.table::setnames(data.table::copy(p), "v", "P"),
                    data.table::setnames(data.table::copy(tm), "v", "Tm")))
  df <- as.data.frame(df[is.finite(y) & is.finite(P) & is.finite(Tm)])
  if (nrow(df) < 30L) NULL else df
}

# M0: SWE ~ P + T over the full record (coefficients, standardised betas, R2)
.swe_tbl_regression <- function(df, rid) {
  m <- stats::lm(y ~ P + Tm, data = df); cs <- summary(m)$coefficients
  sdy <- stats::sd(df$y)
  data.table::data.table(region_id = rid, term = c("Intercept", "Nov-Mar precipitation (mm)", "Nov-Mar temperature (degC)"),
                         estimate = cs[, 1], std_error = cs[, 2], t_value = cs[, 3], p_value = cs[, 4],
                         std_beta = c(NA_real_, cs[2, 1] * stats::sd(df$P) / sdy, cs[3, 1] * stats::sd(df$Tm) / sdy),
                         r_squared = summary(m)$r.squared, n_years = nrow(df),
                         cor_T_P = stats::cor(df$Tm, df$P))
}

# Nested models: how the SWE 'step' coefficient changes as P and T are added (AR(1) GLS where possible)
.swe_tbl_nested <- function(df, rid, step_year) {
  df$step <- as.integer(df$year >= step_year); df$tt <- df$year - min(df$year)
  forms <- list(`SWE ~ step` = y ~ step, `SWE ~ P + step` = y ~ P + step,
                `SWE ~ T + step` = y ~ Tm + step, `SWE ~ P + T + step` = y ~ P + Tm + step)
  data.table::rbindlist(lapply(names(forms), function(nm) {
    f <- forms[[nm]]
    fit <- tryCatch(suppressWarnings(nlme::gls(f, data = df, correlation = nlme::corAR1(form = ~ tt),
                                               method = "ML")), error = function(e) NULL)
    ols <- stats::lm(f, data = df)
    if (!is.null(fit)) {
      ct <- summary(fit)$tTable; meth <- "GLS-AR(1)"
      est <- ct["step", "Value"]; se <- ct["step", "Std.Error"]; p <- ct["step", "p-value"]
    } else {
      ct <- summary(ols)$coefficients; meth <- "OLS"
      est <- ct["step", 1]; se <- ct["step", 2]; p <- ct["step", 4]
    }
    dfree <- nrow(df) - length(stats::coef(ols))
    data.table::data.table(region_id = rid, model = nm, method = meth, step_year = step_year,
                           step_coef_mm = est, std_error = se, lo95 = est - stats::qt(0.975, dfree) * se,
                           hi95 = est + stats::qt(0.975, dfree) * se, p_value = p,
                           r_squared_ols = summary(ols)$r.squared, n_years = nrow(df),
                           cor_T_step = stats::cor(df$Tm, df$step))
  }))
}

# Out-of-sample test: calibrate SWE ~ P + T on pre-step years only, predict the post-step years
.swe_tbl_oos <- function(df, rid, step_year) {
  pre <- df[df$year < step_year, ]; post <- df[df$year >= step_year, ]
  if (nrow(pre) < 15L || nrow(post) < 10L) return(NULL)
  m <- stats::lm(y ~ P + Tm, data = pre)
  e_pt <- post$y - stats::predict(m, newdata = post)
  e_cl <- post$y - mean(pre$y)
  tt <- tryCatch(stats::t.test(e_pt), error = function(e) NULL)
  data.table::data.table(region_id = rid, step_year = step_year, n_pre = nrow(pre), n_post = nrow(post),
                         r2_pre = summary(m)$r.squared, rmse_pre = sqrt(mean(stats::residuals(m)^2)),
                         mean_obs_post = mean(post$y), mean_pred_post_PT = mean(post$y) - mean(e_pt),
                         bias_PT_mm = mean(e_pt), bias_PT_lo = if (is.null(tt)) NA_real_ else tt$conf.int[1],
                         bias_PT_hi = if (is.null(tt)) NA_real_ else tt$conf.int[2],
                         bias_PT_p = if (is.null(tt)) NA_real_ else tt$p.value,
                         rmse_post_PT = sqrt(mean(e_pt^2)), bias_climatology_mm = mean(e_cl),
                         bias_PT_pct_of_pre_mean = 100 * mean(e_pt) / mean(pre$y),
                         bias_climatology_pct = 100 * mean(e_cl) / mean(pre$y))
}

# ---- figure helpers ---------------------------------------------------------------------------
.swe_theme <- function(base = 10) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = base + 1, hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = base - 2, colour = "grey30", hjust = 0),
                   panel.grid.major = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
                   strip.text = ggplot2::element_text(face = "bold", size = base - 1),
                   legend.position = "bottom")
}
# save PNG into report/figures and return the path relative to the report folder
.swe_save_fig <- function(stem, plot, w, h, dpi = 200) {
  if (is.null(plot)) return(NULL)
  fdir <- file.path(SWE_REPORT_DIR, "figures")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  saved <- .save_gg(file.path(fdir, stem), plot, w, h, dpi = dpi, formats = "png")
  if (!length(saved)) return(NULL)
  file.path("figures", saved[1])
}

# Fig 1/2: annual series + period means + fixed step year + Pettitt break
.swe_fig_series <- function(annual, rid, sname, step_year, title) {
  d <- .swe_get_series(annual, rid, sname)
  if (nrow(d) < 20L) return(NULL)
  pc <- .pettitt_curve(d$value)
  ss <- .split_stats(d$value, d$year, step_year, B = 500L)
  brk <- d$year[pc$tau]
  unit <- .SWE_SERIES_UNIT[[sname]]
  seg <- data.frame(x0 = c(min(d$year), step_year - 0.5), x1 = c(step_year - 0.5, max(d$year)),
                    y = c(ss$mean_before, ss$mean_after),
                    grp = factor(c(sprintf("mean before %d", step_year), sprintf("mean from %d", step_year)),
                                 levels = c(sprintf("mean before %d", step_year), sprintf("mean from %d", step_year))))
  sub1 <- sprintf("Pettitt (searches all years): p = %s, break after %d", .fmt_p(pc$p), brk)
  sub2 <- sprintf("Fixed split at %d: %s %s (%s%%), Welch p = %s", step_year,
                  .fmt_sgn(ss$diff, 2), unit, .fmt_sgn(ss$pct, 1), .fmt_p(ss$welch_p))
  ggplot2::ggplot(d, ggplot2::aes(year, value)) +
    ggplot2::geom_line(colour = .SWE_SERIES_COL[[sname]], linewidth = 0.4, alpha = 0.7) +
    ggplot2::geom_point(colour = .SWE_SERIES_COL[[sname]], size = 1.3) +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = x0, xend = x1, y = y, yend = y, colour = grp),
                          inherit.aes = FALSE, linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = step_year - 0.5, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_vline(xintercept = brk + 0.5, linetype = "dotted", colour = "darkorange", linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = stats::setNames(c("firebrick", "steelblue"), levels(seg$grp)), name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(1950, 2020, by = 10)) +
    ggplot2::labs(title = title, subtitle = paste0(sub1, "\n", sub2,
                                                   "   (dashed = a-priori step year; dotted orange = Pettitt break)"),
                  x = "Year", y = sprintf("%s (%s)", .SWE_SERIES_LAB[[sname]], unit)) +
    .swe_theme()
}

# Fig 3: observed vs P,T-predicted Mar-Apr SWE, and the residual
.swe_fig_obs_fit <- function(annual, rid, step_year) {
  df <- .swe_model_frame(annual, rid)
  if (is.null(df)) return(NULL)
  m <- stats::lm(y ~ P + Tm, data = df)
  df$fit <- stats::fitted(m); df$res <- stats::residuals(m)
  pc <- .pettitt_curve(df$res); ss <- .split_stats(df$res, df$year, step_year, B = 500L)
  r2 <- summary(m)$r.squared
  long <- data.frame(year = rep(df$year, 2), value = c(df$y, df$fit),
                     what = factor(rep(c("Observed (ERA5-Land)", "Predicted from Nov-Mar P and T"), each = nrow(df)),
                                   levels = c("Observed (ERA5-Land)", "Predicted from Nov-Mar P and T")))
  pA <- ggplot2::ggplot(long, ggplot2::aes(year, value, colour = what)) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.1) +
    ggplot2::geom_vline(xintercept = step_year - 0.5, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(values = c("#01665e", "#e08214"), name = NULL) +
    ggplot2::labs(title = "A) Mar-Apr SWE: observed vs predicted from Nov-Mar precipitation and temperature",
                  subtitle = sprintf("SWE ~ P + T fitted over %d-%d (n = %d); R-squared = %.2f", min(df$year), max(df$year), nrow(df), r2),
                  x = NULL, y = "SWE (mm w.e.)") + .swe_theme()
  mu_b <- mean(df$res[df$year < step_year]); mu_a <- mean(df$res[df$year >= step_year])
  seg <- data.frame(x0 = c(min(df$year), step_year - 0.5), x1 = c(step_year - 0.5, max(df$year)), y = c(mu_b, mu_a))
  pB <- ggplot2::ggplot(df, ggplot2::aes(year, res)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_col(fill = "#8c510a", alpha = 0.55, width = 0.8) +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE,
                          colour = "black", linewidth = 1.2) +
    ggplot2::geom_vline(xintercept = step_year - 0.5, linetype = "dashed", colour = "grey40") +
    ggplot2::labs(title = "B) Residual after removing the P and T signal",
                  subtitle = sprintf("Pettitt p = %s | change at %d: %s mm (Welch p = %s). Black lines = period means.",
                                     .fmt_p(pc$p), step_year, .fmt_sgn(ss$diff, 1), .fmt_p(ss$welch_p)),
                  x = "Year", y = "Residual (mm w.e.)") + .swe_theme()
  pA / pB + patchwork::plot_layout(heights = c(1.1, 1))
}

# Fig 4: Pettitt U_t curve for each series (shows the search over every candidate break)
.swe_fig_pettitt <- function(annual, rid, step_year) {
  rows <- lapply(names(.SWE_SERIES_LAB), function(s) {
    d <- if (s == "SWE_MarApr_resid_PT") {
      df <- .swe_model_frame(annual, rid)
      if (is.null(df)) NULL else data.frame(year = df$year, value = stats::residuals(stats::lm(y ~ P + Tm, data = df)))
    } else .swe_get_series(annual, rid, s)
    if (is.null(d) || nrow(d) < 20L) return(NULL)
    pc <- .pettitt_curve(d$value)
    list(curve = data.frame(year = d$year[seq_along(pc$U)], U = pc$U, lab = .SWE_SERIES_LAB[[s]]),
         ann = data.frame(lab = .SWE_SERIES_LAB[[s]], year = d$year[pc$tau], U = pc$U[pc$tau],
                          txt = sprintf("K = %.0f, p %s\nbreak after %d", pc$K, ifelse(pc$p < 0.001, "< 0.001", paste0("= ", .fmt_p(pc$p))), d$year[pc$tau]),
                          pos = pc$U[pc$tau] >= 0, left_half = d$year[pc$tau] < stats::median(d$year)))
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  cv <- data.table::rbindlist(lapply(rows, `[[`, "curve")); an <- data.table::rbindlist(lapply(rows, `[[`, "ann"))
  lv <- unname(.SWE_SERIES_LAB[names(.SWE_SERIES_LAB) %in% names(.SWE_SERIES_LAB)])
  cv[, lab := factor(lab, levels = lv)]; an[, lab := factor(lab, levels = lv)]
  # label in the empty corner on the side away from the peak
  an[, `:=`(xp = ifelse(left_half, Inf, -Inf), yp = ifelse(pos, Inf, -Inf),
            hj = ifelse(left_half, 1.05, -0.05), vj = ifelse(pos, 1.25, -0.4))]
  ggplot2::ggplot(cv, ggplot2::aes(year, U)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_line(linewidth = 0.8, colour = "grey20") +
    ggplot2::geom_vline(xintercept = step_year - 0.5, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(data = an, size = 2.6, colour = "darkorange") +
    ggplot2::geom_text(data = an, ggplot2::aes(x = xp, y = yp, label = txt, hjust = hj, vjust = vj), size = 3.0, inherit.aes = FALSE) +
    ggplot2::facet_wrap(~lab, ncol = 2, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = seq(1950, 2020, by = 20), expand = ggplot2::expansion(mult = 0.03)) +
    ggplot2::labs(title = "Pettitt statistic U(t) across every candidate break year",
                  subtitle = paste0("The peak |U| marks the most probable break (orange); its p-value already accounts for ",
                                    "searching all years. Dashed = a-priori step year."),
                  x = "Year (last year of the earlier regime)", y = "Pettitt U(t)") + .swe_theme()
}

# Fig 5: standardised step (Cohen's d, 95% CI) at the step year, by region and series
.swe_fig_forest <- function(tests, step_year) {
  d <- tests[series %in% names(.SWE_SERIES_LAB) & is.finite(cohens_d) & is.finite(n_before) & is.finite(n_after)]
  if (!nrow(d)) return(NULL)
  d <- data.table::copy(d)
  d[, se_d := sqrt((n_before + n_after) / (n_before * n_after) + cohens_d^2 / (2 * (n_before + n_after)))]
  d[, `:=`(lo = cohens_d - 1.96 * se_d, hi = cohens_d + 1.96 * se_d)]
  bid <- .swe_basin_id()
  d[, region := factor(region_id, levels = rev(c(bid, setdiff(sort(unique(region_id)), bid))))]
  d[, lab := factor(.SWE_SERIES_LAB[series], levels = unname(.SWE_SERIES_LAB))]
  d[, sig := ifelse(is.finite(welch_p) & welch_p < 0.05, "Welch p < 0.05", "not significant")]
  ggplot2::ggplot(d, ggplot2::aes(cohens_d, region, colour = sig)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = lo, xmax = hi), size = 0.35) +
    ggplot2::facet_wrap(~lab, nrow = 1) +
    ggplot2::scale_colour_manual(values = c("Welch p < 0.05" = "#b2182b", "not significant" = "grey40"), name = NULL) +
    ggplot2::labs(title = sprintf("Standardised change after %d (Cohen's d, 95%% CI) by region", step_year),
                  subtitle = "Negative = lower after the step year. Sub-basins are nested/overlapping, so panels are not independent.",
                  x = "Cohen's d", y = NULL) + .swe_theme(9)
}

# Fig 6: how the before/after contrast depends on the chosen step year (basin)
.swe_scan_split <- function(d, cands, min_n = 10L) {
  data.table::rbindlist(lapply(cands, function(cy) {
    if (sum(d$year < cy) < min_n || sum(d$year >= cy) < min_n) return(NULL)
    b <- d$value[d$year < cy]; a <- d$value[d$year >= cy]
    sp <- sqrt(((length(b) - 1) * stats::var(b) + (length(a) - 1) * stats::var(a)) / (length(a) + length(b) - 2))
    data.table::data.table(cut_year = cy, d = (mean(a) - mean(b)) / sp,
                           p = tryCatch(stats::t.test(a, b)$p.value, error = function(e) NA_real_))
  }))
}
.swe_fig_scan <- function(annual, rid, step_year, brk_by_series) {
  parts <- lapply(names(.SWE_SERIES_LAB), function(s) {
    d <- if (s == "SWE_MarApr_resid_PT") {
      df <- .swe_model_frame(annual, rid)
      if (is.null(df)) NULL else data.frame(year = df$year, value = stats::residuals(stats::lm(y ~ P + Tm, data = df)))
    } else .swe_get_series(annual, rid, s)
    if (is.null(d) || nrow(d) < 30L) return(NULL)
    r <- .swe_scan_split(d, seq(1962L, 1995L)); if (is.null(r)) return(NULL)
    r[, lab := .SWE_SERIES_LAB[[s]]]; r
  })
  sc <- data.table::rbindlist(Filter(Negate(is.null), parts))
  if (!nrow(sc)) return(NULL)
  lv <- unname(.SWE_SERIES_LAB)
  long <- data.table::rbindlist(list(
    sc[, .(cut_year, lab, metric = "Cohen's d (after - before)", value = d)],
    sc[, .(cut_year, lab, metric = "-log10(Welch p)", value = -log10(pmax(p, 1e-12)))]))
  long[, lab := factor(lab, levels = lv)]
  long[, metric := factor(metric, levels = c("Cohen's d (after - before)", "-log10(Welch p)"))]
  hl <- data.frame(metric = factor("-log10(Welch p)", levels = levels(long$metric)), yint = -log10(0.05))
  ggplot2::ggplot(long, ggplot2::aes(cut_year, value, colour = lab)) +
    ggplot2::geom_hline(data = hl, ggplot2::aes(yintercept = yint), linetype = "dotted", colour = "grey40") +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = step_year, linetype = "dashed", colour = "grey40") +
    ggplot2::facet_grid(metric ~ lab, scales = "free_y") +
    ggplot2::scale_colour_manual(values = stats::setNames(unname(.SWE_SERIES_COL), lv), guide = "none") +
    ggplot2::labs(title = "Sensitivity of the before/after contrast to the chosen step year",
                  subtitle = paste0("Each point splits the record at that year (first year of the later regime). Dashed = ",
                                    step_year, "; dotted = p = 0.05. Descriptive only: p-values are NOT adjusted for scanning years."),
                  x = "First year of the later regime", y = NULL) + .swe_theme(9)
}

# Fig 7 + Table 9: month-by-month change after the step year (basin; P, T, SWE)
.swe_monthly_change <- function(step_year) {
  b <- basin_avg_monthly
  vars <- list(`Precipitation (mm/month)` = b$precip_mm_month, `Temperature (degC)` = b$tair_degC_month,
               `SWE (mm w.e.)` = b$swe_mm_month)
  yr <- as.integer(format(b$date, "%Y")); mo <- as.integer(format(b$date, "%m"))
  data.table::rbindlist(lapply(names(vars), function(vn) {
    v <- as.numeric(vars[[vn]])
    if (!length(v) || !any(is.finite(v))) return(NULL)
    data.table::rbindlist(lapply(1:12, function(m) {
      a <- v[mo == m & yr >= step_year]; bb <- v[mo == m & yr < step_year]
      a <- a[is.finite(a)]; bb <- bb[is.finite(bb)]
      if (length(a) < 5L || length(bb) < 5L) return(NULL)
      tt <- tryCatch(stats::t.test(a, bb), error = function(e) NULL)
      data.table::data.table(variable = vn, month = m, month_abb = month.abb[m], mean_before = mean(bb), mean_after = mean(a),
                             diff = mean(a) - mean(bb), welch_p = if (is.null(tt)) NA_real_ else tt$p.value)
    }))
  }))
}
.swe_fig_monthly <- function(mc, step_year) {
  if (is.null(mc) || !nrow(mc)) return(NULL)
  mc <- data.table::copy(mc)
  mc[, month_abb := factor(month_abb, levels = month.abb)]
  mc[, variable := factor(variable, levels = unique(variable))]
  mc[, sig := is.finite(welch_p) & welch_p < 0.05]
  ggplot2::ggplot(mc, ggplot2::aes(month_abb, diff, fill = variable, alpha = sig)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(data = mc[sig == TRUE], ggplot2::aes(label = "*"), vjust = ifelse(mc[sig == TRUE]$diff >= 0, -0.1, 1.1),
                       size = 5, show.legend = FALSE, alpha = 1) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4) +
    ggplot2::facet_wrap(~variable, ncol = 1, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("#4575b4", "#6a3d9a", "#01665e"), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.12))) +
    ggplot2::labs(title = sprintf("Change in the monthly mean, %d-%d vs 1950-%d (basin)", step_year, 2025L, step_year - 1L),
                  subtitle = "Solid bars / * = Welch p < 0.05 for that calendar month. Descriptive; no multiplicity adjustment.",
                  x = NULL, y = "Mean after minus mean before") + .swe_theme(9)
}

# Fig 8: pixel-level step vs snow regime (scatter) and its map
.swe_fig_pixels <- function(pix, step_year) {
  if (is.null(pix) || !nrow(pix) || !all(c("pct_change", "mean_marapr_mm", "pettitt_p") %in% names(pix))) return(NULL)
  d <- data.table::copy(pix[is.finite(pct_change) & is.finite(mean_marapr_mm)])
  d[, sig := ifelse(is.finite(pettitt_p_fdr) & pettitt_p_fdr < 0.05, "Pettitt FDR < 0.05", "not significant")]
  pS <- ggplot2::ggplot(d, ggplot2::aes(mean_marapr_mm, pct_change, colour = sig)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_point(size = 1.4, alpha = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c("Pettitt FDR < 0.05" = "#b2182b", "not significant" = "grey55"), name = NULL) +
    ggplot2::labs(title = sprintf("A) Change after %d vs mean Mar-Apr SWE (each pixel)", step_year),
                  x = "Mean Mar-Apr SWE (mm w.e., log scale)", y = "Change after step year (%)") + .swe_theme(9)
  pM <- ggplot2::ggplot(d, ggplot2::aes(x, y, fill = pct_change)) +
    ggplot2::geom_tile(width = .grid_spacing(d$x), height = .grid_spacing(d$y)) +
    ggplot2::scale_fill_gradientn(colours = .pal_wet_positive(101), limits = c(-max(abs(d$pct_change)), max(abs(d$pct_change))),
                                  name = "% change") +
    { if (!is.null(basin_wgs)) ggplot2::geom_sf(data = basin_wgs, inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.6) } +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(title = "B) Map of the change", x = NULL, y = NULL) + ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10), legend.position = "bottom")
  pS | pM
}

# ---- one row of key statistics for a (region, series) -----------------------------------------
.swe_key_row <- function(annual, tests, rid, s, sy, B = 5000L) {
  d <- .swe_get_series(annual, rid, s)
  if (nrow(d) < 20L) return(NULL)
  pc <- .pettitt_curve(d$value)
  brk_last <- d$year[pc$tau]; brk_first <- d$year[pc$tau + 1L]
  at_step <- .split_stats(d$value, d$year, sy, B = B)
  at_brk  <- .split_stats(d$value, d$year, brk_first, B = 200L)
  tr <- tests[region_id == rid & series == s]
  g  <- function(col) if (nrow(tr) && col %in% names(tr)) tr[[col]][1] else NA_real_
  if (s %in% c("Temp_NovMar", "SWE_MarApr_resid_PT")) {   # means near/below zero: a percent change is meaningless
    at_step[, c("pct", "pct_lo", "pct_hi", "mde_80_pct") := NA_real_]; at_brk[, pct := NA_real_]
  }
  cbind(
    data.table::data.table(region_id = rid, series = s, label = .SWE_SERIES_LAB[[s]], unit = .SWE_SERIES_UNIT[[s]],
                           n_years = nrow(d), pettitt_p = pc$p, pettitt_p_pipeline = g("pettitt_p"),
                           pettitt_last_year_before = brk_last, pettitt_first_year_after = brk_first,
                           snht_p = g("snht_p"), snht_first_year_after = g("snht_first_year_after"), buishand_p = g("buishand_range_p"),
                           step_year = sy),
    at_step[, .(n_before, n_after, mean_before, mean_after, diff, diff_lo, diff_hi, boot_lo, boot_hi, pct, pct_lo, pct_hi,
                cohens_d, d_lo, d_hi, welch_p, mannwhitney_p, sd_pooled, mde_80, mde_80_pct)],
    data.table::data.table(welch_p_pipeline = g("welch_p"), diff_at_pettitt_break = at_brk$diff,
                           pct_at_pettitt_break = at_brk$pct))
}

.sig_txt <- function(p) if (!is.finite(p)) "not available" else if (p < 0.05) "significant" else "not significant"
.step_found <- function(r) isTRUE(r$pettitt_p < 0.05) || isTRUE(r$welch_p < 0.05)

####################################################################################
# FUNCTION 19: SWE step-change figures + tables (no narrative report is written)
####################################################################################
create_swe_step_change_outputs <- function() {
  cat("\n-- Function 19: SWE step-change tables and figures --\n")
  dat <- .swe_load()
  if (is.null(dat$annual) || is.null(dat$tests)) {
    cat(sprintf("  SWE homogeneity tables not found in %s - run 9pr_pet_trends.r (v7) first.\n", SWE_HOMOG_DIR))
    return(invisible(NULL))
  }
  annual <- dat$annual; tests <- dat$tests
  dir.create(file.path(SWE_REPORT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
  set.seed(20260101)
  sy  <- .swe_step_year(); bid <- .swe_basin_id()
  rids <- unique(annual$region_id); rids <- c(intersect(bid, rids), setdiff(sort(rids), bid))
  if (!bid %in% rids) { cat("  Basin series not found in the annual table - skipped.\n"); return(invisible(NULL)) }
  save_tbl <- function(dt, stem) {
    data.table::fwrite(dt, .resolve_out(file.path(SWE_REPORT_DIR, "tables", paste0(stem, ".csv"))))
    invisible(dt)
  }
  ser4 <- names(.SWE_SERIES_LAB)
  
  # ---- tables ----------------------------------------------------------------------------------
  cat("  tables...\n")
  key_all <- data.table::rbindlist(lapply(rids, function(r)
    data.table::rbindlist(lapply(ser4, function(s) .swe_key_row(annual, tests, r, s, sy, B = if (r == bid) 5000L else 1000L)))),
    fill = TRUE)
  save_tbl(key_all, "Table08_all_regions_four_series")
  key <- key_all[region_id == bid]
  if (!nrow(key)) { cat("  No basin series long enough for the step tests - skipped.\n"); return(invisible(NULL)) }
  save_tbl(key, "Table01_key_results_basin")
  save_tbl(key[, .(series = label, pettitt_p, pettitt_first_year_after, snht_p, snht_first_year_after, buishand_p)],
           "Table03_break_years_by_test")
  kP <- key[series == "Precip_NovMar"]; kT <- key[series == "Temp_NovMar"]
  kS <- key[series == "SWE_MarApr"];    kR <- key[series == "SWE_MarApr_resid_PT"]
  
  # integrity check: independent recomputation vs the pipeline's own p-values
  chk <- suppressWarnings(key_all[, .(max_abs_diff_pettitt_p = max(abs(pettitt_p - pettitt_p_pipeline), na.rm = TRUE),
                                      max_abs_diff_welch_p = max(abs(welch_p - welch_p_pipeline), na.rm = TRUE),
                                      n_compared = sum(is.finite(pettitt_p_pipeline)))])
  chk[!is.finite(max_abs_diff_pettitt_p), max_abs_diff_pettitt_p := NA_real_]
  chk[!is.finite(max_abs_diff_welch_p), max_abs_diff_welch_p := NA_real_]
  cat(sprintf("  cross-check vs pipeline: max |dPettitt p| = %.2e, max |dWelch p| = %.2e (%d series compared)\n",
              chk$max_abs_diff_pettitt_p, chk$max_abs_diff_welch_p, chk$n_compared))
  
  reg_all <- data.table::rbindlist(lapply(rids, function(r) { df <- .swe_model_frame(annual, r)
  if (is.null(df)) NULL else .swe_tbl_regression(df, r) }))
  nest_all <- data.table::rbindlist(lapply(rids, function(r) { df <- .swe_model_frame(annual, r)
  if (is.null(df)) NULL else .swe_tbl_nested(df, r, sy) }))
  oos_all <- data.table::rbindlist(lapply(rids, function(r) { df <- .swe_model_frame(annual, r)
  if (is.null(df)) NULL else .swe_tbl_oos(df, r, sy) }))
  save_tbl(reg_all, "Table04_regression_SWE_on_P_T"); save_tbl(nest_all, "Table05_step_coefficient_nested_models")
  save_tbl(oos_all, "Table06_out_of_sample_prediction")
  save_tbl(key[, .(label, unit, mean_before, diff, diff_lo, diff_hi, boot_lo, boot_hi, pct, pct_lo, pct_hi, sd_pooled, mde_80, mde_80_pct)],
           "Table02_effect_size_intervals_and_detectable_step")
  sub_b <- if (!is.null(dat$sub)) dat$sub[region_id == bid & series %in% c("Precip_NovMar", "Temp_NovMar", "SWE_MarApr")] else NULL
  if (!is.null(sub_b)) save_tbl(sub_b, "Table10_subperiod_trends_basin")
  era_b <- if (!is.null(dat$era)) dat$era[region_id == bid & series %in% ser4] else NULL
  if (!is.null(era_b)) save_tbl(era_b, "Table11_period_means_basin")
  ref_b <- if (!is.null(dat$refs)) dat$refs[region_id == bid] else NULL
  if (!is.null(ref_b)) save_tbl(ref_b, "Table09_reference_sensitivity_basin")
  mc <- tryCatch(.swe_monthly_change(sy), error = function(e) NULL)
  if (!is.null(mc)) save_tbl(mc, "Table12_monthly_change_basin")
  pix_sum <- NULL
  if (!is.null(dat$pix) && all(c("pct_change", "mean_marapr_mm", "pettitt_p") %in% names(dat$pix))) {
    pp <- data.table::copy(dat$pix[is.finite(pct_change) & is.finite(mean_marapr_mm)])
    pp[, tercile := cut(mean_marapr_mm, stats::quantile(mean_marapr_mm, 0:3 / 3), include.lowest = TRUE,
                        labels = c("lowest third (least snow)", "middle third", "highest third (most snow)"))]
    pix_sum <- rbind(
      pp[, .(group = "all pixels", n_pixels = .N, pct_negative = 100 * mean(diff_mm < 0), median_pct_change = stats::median(pct_change),
             pct_sig_pettitt = 100 * mean(pettitt_p < 0.05, na.rm = TRUE), pct_sig_pettitt_fdr = 100 * mean(pettitt_p_fdr < 0.05, na.rm = TRUE))],
      pp[, .(n_pixels = .N, pct_negative = 100 * mean(diff_mm < 0), median_pct_change = stats::median(pct_change),
             pct_sig_pettitt = 100 * mean(pettitt_p < 0.05, na.rm = TRUE), pct_sig_pettitt_fdr = 100 * mean(pettitt_p_fdr < 0.05, na.rm = TRUE)),
         by = .(group = as.character(tercile))][order(match(group, c("lowest third (least snow)", "middle third", "highest third (most snow)")))])
    save_tbl(pix_sum, "Table07_pixel_summary")
  }
  
  # ---- figures -------------------------------------------------------------------------------------
  cat("  figures...\n")
  fig <- list()
  fig$swe   <- .swe_save_fig("SWE_Fig6_MarApr_SWE_series", .swe_fig_series(annual, bid, "SWE_MarApr", sy, "Mar-Apr SWE (basin mean)"), 9, 4.6)
  fig$p     <- .swe_save_fig("SWE_Fig1_NovMar_precipitation", .swe_fig_series(annual, bid, "Precip_NovMar", sy, "Nov-Mar precipitation (basin mean)"), 9, 4.6)
  fig$t     <- .swe_save_fig("SWE_Fig2_NovMar_temperature", .swe_fig_series(annual, bid, "Temp_NovMar", sy, "Nov-Mar air temperature (basin mean)"), 9, 4.6)
  fig$fit   <- .swe_save_fig("SWE_Fig7_observed_vs_predicted_and_residual", .swe_fig_obs_fit(annual, bid, sy), 9, 8)
  fig$pet   <- .swe_save_fig("SWE_Fig3_Pettitt_statistic_curves", .swe_fig_pettitt(annual, bid, sy), 10, 7)
  fig$forest<- .swe_save_fig("SWE_Fig9_effect_size_by_region", .swe_fig_forest(tests, sy), 13, 4.8)
  fig$scan  <- .swe_save_fig("SWE_Fig4_sensitivity_to_step_year", .swe_fig_scan(annual, bid, sy), 13, 6.5)
  fig$month <- .swe_save_fig("SWE_Fig5_monthly_change", .swe_fig_monthly(mc, sy), 9, 8.5)
  fig$pix   <- .swe_save_fig("SWE_Fig8_pixel_step_vs_snow_regime", .swe_fig_pixels(dat$pix, sy), 12, 5.2)
  # ---- report ----------------------------------------------------------------------------------------
  cat("  key numbers and manifest...\n")
  nb <- nest_all[region_id == bid]; m1 <- nb[model == "SWE ~ step"]; m4 <- nb[model == "SWE ~ P + T + step"]
  m2 <- nb[model == "SWE ~ P + step"]; m3 <- nb[model == "SWE ~ T + step"]
  oo <- oos_all[region_id == bid]; rg <- reg_all[region_id == bid]
  mp <- if (!is.null(mc)) mc else NULL
  signed <- function(x, d = 1L) .fmt_sgn(x, d)
  sw <- if (nrow(kS)) .step_found(kS) else NA; tp <- if (nrow(kT)) .step_found(kT) else NA
  pp_ <- if (nrow(kP)) .step_found(kP) else NA; rs <- if (nrow(kR)) .step_found(kR) else NA
  # Stronger checks (rules of thumb, not formal tests): the residual test alone is partly circular
  exp_frac <- if (nrow(kS) && nrow(kR) && is.finite(kS$diff) && kS$diff != 0) 100 * (1 - kR$diff / kS$diff) else NA_real_
  nest_ok <- nrow(m1) > 0 && nrow(m4) > 0 && is.finite(m4$p_value) &&
    (m4$p_value >= 0.05 || abs(m4$step_coef_mm) < 0.5 * abs(m1$step_coef_mm))
  oos_ok  <- nrow(oo) > 0 && is.finite(oo$bias_PT_mm) &&
    (isTRUE(oo$bias_PT_p >= 0.05) || abs(oo$bias_PT_mm) < 0.5 * abs(oo$bias_climatology_mm))
  corroborated <- isTRUE(nest_ok) && isTRUE(oos_ok) && isTRUE(exp_frac >= 50)
  exp_txt <- if (is.finite(exp_frac)) sprintf("; P and T account for %s of the raw SWE change", .fmt_expl(exp_frac)) else ""
  chk_tbl <- NULL
  if (nrow(kR) && nrow(m1) && nrow(m4) && nrow(oo)) {
    chk_tbl <- data.table::data.table(
      Check = c("Residual shows no significant step (Pettitt and Welch both p >= 0.05)",
                "P and T account for at least 50% of the raw SWE change",
                "Step coefficient with P and T is less than half the no-covariate step, or not significant",
                "Out-of-sample bias with P and T is not significant, or less than half the no-covariate bias"),
      Value = c(sprintf("Pettitt p = %s; Welch p = %s", .fmt_p(kR$pettitt_p), .fmt_p(kR$welch_p)),
                .fmt_expl(exp_frac),
                sprintf("%s vs %s mm (p = %s)", signed(m4$step_coef_mm, 1), signed(m1$step_coef_mm, 1), .fmt_p(m4$p_value)),
                sprintf("%s vs %s mm (p = %s)", signed(oo$bias_PT_mm, 1), signed(oo$bias_climatology_mm, 1), .fmt_p(oo$bias_PT_p))),
      Passes = c(ifelse(isTRUE(!rs), "yes", "no"), ifelse(isTRUE(exp_frac >= 50), "yes", "no"),
                 ifelse(isTRUE(nest_ok), "yes", "no"), ifelse(isTRUE(oos_ok), "yes", "no")))
    save_tbl(chk_tbl, "Table05a_evidence_checklist")
  }
  # ---- key numbers behind the two headline statements -------------------------------------------------
  kn <- function(group, item, value, unit, source) data.table::data.table(group = group, item = item,
                                                                          value = as.numeric(value), unit = unit, source = source)
  T1 <- "Table01_key_results_basin.csv"
  rowsT <- if (nrow(kT)) list(
    kn("Nov-Mar temperature", "Pettitt p (searches all years)", kT$pettitt_p, "", paste0(T1, ": pettitt_p")),
    kn("Nov-Mar temperature", "Pettitt p as written by the pipeline", kT$pettitt_p_pipeline, "", paste0(T1, ": pettitt_p_pipeline")),
    kn("Nov-Mar temperature", "Pettitt break: last year of earlier regime", kT$pettitt_last_year_before, "year", paste0(T1, ": pettitt_last_year_before")),
    kn("Nov-Mar temperature", "SNHT p", kT$snht_p, "", paste0(T1, ": snht_p")),
    kn("Nov-Mar temperature", "Buishand range p", kT$buishand_p, "", paste0(T1, ": buishand_p")),
    kn("Nov-Mar temperature", "A-priori step year (first year of later regime)", sy, "year", paste0(T1, ": step_year")),
    kn("Nov-Mar temperature", "Mean before step year", kT$mean_before, "degC", paste0(T1, ": mean_before")),
    kn("Nov-Mar temperature", "Mean from step year", kT$mean_after, "degC", paste0(T1, ": mean_after")),
    kn("Nov-Mar temperature", "Change at step year (after - before)", kT$diff, "degC", paste0(T1, ": diff")),
    kn("Nov-Mar temperature", "95% CI of change, low", kT$diff_lo, "degC", paste0(T1, ": diff_lo")),
    kn("Nov-Mar temperature", "95% CI of change, high", kT$diff_hi, "degC", paste0(T1, ": diff_hi")),
    kn("Nov-Mar temperature", "Change at the Pettitt-detected break", kT$diff_at_pettitt_break, "degC", paste0(T1, ": diff_at_pettitt_break")),
    kn("Nov-Mar temperature", "Welch p at the step year", kT$welch_p, "", paste0(T1, ": welch_p")),
    kn("Nov-Mar temperature", "Mann-Whitney p at the step year", kT$mannwhitney_p, "", paste0(T1, ": mannwhitney_p")),
    kn("Nov-Mar temperature", "Cohen's d at the step year", kT$cohens_d, "", paste0(T1, ": cohens_d"))) else list()
  rowsR <- if (nrow(kR)) list(
    kn("SWE residual after removing P and T", "Pettitt p (searches all years)", kR$pettitt_p, "", paste0(T1, ": pettitt_p")),
    kn("SWE residual after removing P and T", "Pettitt p as written by the pipeline", kR$pettitt_p_pipeline, "", paste0(T1, ": pettitt_p_pipeline")),
    kn("SWE residual after removing P and T", "Pettitt break: last year of earlier regime", kR$pettitt_last_year_before, "year", paste0(T1, ": pettitt_last_year_before")),
    kn("SWE residual after removing P and T", "Mean residual before step year", kR$mean_before, "mm w.e.", paste0(T1, ": mean_before")),
    kn("SWE residual after removing P and T", "Mean residual from step year", kR$mean_after, "mm w.e.", paste0(T1, ": mean_after")),
    kn("SWE residual after removing P and T", "Residual change at step year (after - before)", kR$diff, "mm w.e.", paste0(T1, ": diff")),
    kn("SWE residual after removing P and T", "95% CI of residual change, low", kR$diff_lo, "mm w.e.", paste0(T1, ": diff_lo")),
    kn("SWE residual after removing P and T", "95% CI of residual change, high", kR$diff_hi, "mm w.e.", paste0(T1, ": diff_hi")),
    kn("SWE residual after removing P and T", "Welch p at the step year", kR$welch_p, "", paste0(T1, ": welch_p")),
    kn("SWE residual after removing P and T", "Mann-Whitney p at the step year", kR$mannwhitney_p, "", paste0(T1, ": mannwhitney_p")),
    kn("SWE residual after removing P and T", "Detectable step with 80% power", kR$mde_80, "mm w.e.", paste0(T1, ": mde_80")),
    kn("SWE residual after removing P and T", "Share of the raw SWE change accounted for by P and T", exp_frac, "%", "1 - residual change / raw SWE change")) else list()
  rowsS <- if (nrow(kS)) list(
    kn("Raw Mar-Apr SWE (for comparison)", "Pettitt p (searches all years)", kS$pettitt_p, "", paste0(T1, ": pettitt_p")),
    kn("Raw Mar-Apr SWE (for comparison)", "Pettitt break: last year of earlier regime", kS$pettitt_last_year_before, "year", paste0(T1, ": pettitt_last_year_before")),
    kn("Raw Mar-Apr SWE (for comparison)", "Change at step year (after - before)", kS$diff, "mm w.e.", paste0(T1, ": diff")),
    kn("Raw Mar-Apr SWE (for comparison)", "Welch p at the step year", kS$welch_p, "", paste0(T1, ": welch_p"))) else list()
  rowsM <- c(
    if (nrow(rg)) list(kn("Regression SWE ~ P + T", "R-squared", rg$r_squared[1], "", "Table04_regression_SWE_on_P_T.csv: r_squared"),
                       kn("Regression SWE ~ P + T", "Years used", rg$n_years[1], "years", "Table04_regression_SWE_on_P_T.csv: n_years")),
    if (nrow(m1)) list(kn("Step coefficient, no covariates", "SWE ~ step", m1$step_coef_mm, "mm w.e.", "Table05_step_coefficient_nested_models.csv: step_coef_mm"),
                       kn("Step coefficient, no covariates", "p", m1$p_value, "", "Table05_step_coefficient_nested_models.csv: p_value")),
    if (nrow(m4)) list(kn("Step coefficient with P and T", "SWE ~ P + T + step", m4$step_coef_mm, "mm w.e.", "Table05_step_coefficient_nested_models.csv: step_coef_mm"),
                       kn("Step coefficient with P and T", "p", m4$p_value, "", "Table05_step_coefficient_nested_models.csv: p_value"),
                       kn("Step coefficient with P and T", "Correlation between temperature and the step term", m4$cor_T_step, "", "Table05_step_coefficient_nested_models.csv: cor_T_step")),
    if (nrow(oo)) list(kn("Out-of-sample (fit before step year)", "Bias with P and T", oo$bias_PT_mm, "mm w.e.", "Table06_out_of_sample_prediction.csv: bias_PT_mm"),
                       kn("Out-of-sample (fit before step year)", "p of that bias", oo$bias_PT_p, "", "Table06_out_of_sample_prediction.csv: bias_PT_p"),
                       kn("Out-of-sample (fit before step year)", "Bias with no covariates", oo$bias_climatology_mm, "mm w.e.", "Table06_out_of_sample_prediction.csv: bias_climatology_mm")))
  key_numbers <- data.table::rbindlist(c(rowsT, rowsR, rowsS, rowsM), use.names = TRUE)
  save_tbl(key_numbers, "Table00_key_numbers")
  
  # ---- manifest of everything written --------------------------------------------------------------------
  tdir <- file.path(SWE_REPORT_DIR, "tables"); fdir <- file.path(SWE_REPORT_DIR, "figures")
  desc <- c(
    Table00_key_numbers = "Every number behind the Nov-Mar temperature step and the SWE-residual statements, with its source file/column",
    Table01_key_results_basin = "Basin: Pettitt/SNHT/Buishand p, break years, before/after means, change with CIs, Welch/Mann-Whitney p, Cohen's d, detectable step",
    Table02_effect_size_intervals_and_detectable_step = "Basin: change with Welch and bootstrap CIs and the smallest step detectable with 80% power",
    Table03_break_years_by_test = "Basin: break year and p-value from each unknown-date test",
    Table04_regression_SWE_on_P_T = "SWE ~ Nov-Mar precipitation + temperature (all regions): coefficients, standardised betas, R-squared",
    Table05_step_coefficient_nested_models = "SWE step coefficient with and without P and T, AR(1) GLS (all regions)",
    Table05a_evidence_checklist = "Basin: four rule-of-thumb checks on whether P and T account for the SWE step",
    Table06_out_of_sample_prediction = "Model fitted before the step year, used to predict later years (all regions)",
    Table07_pixel_summary = "Pixel-level scan of Mar-Apr SWE by snow-regime tercile",
    Table08_all_regions_four_series = "All regions x four series: same statistics as Table01",
    Table09_reference_sensitivity_basin = "Mean standardised SWE by era under three reference periods",
    Table10_subperiod_trends_basin = "Mann-Kendall trends full record / before / after step / 1991-2020",
    Table11_period_means_basin = "Period means (the 1979-1990 window is post hoc, descriptive only)",
    Table12_monthly_change_basin = "Month-by-month change after the step year: precipitation, temperature, SWE",
    SWE_Fig1_NovMar_precipitation = "Nov-Mar precipitation with period means and Pettitt break",
    SWE_Fig2_NovMar_temperature = "Nov-Mar temperature with period means and Pettitt break",
    SWE_Fig3_Pettitt_statistic_curves = "Pettitt U(t) over every candidate break year, four series",
    SWE_Fig4_sensitivity_to_step_year = "Before/after contrast versus chosen step year, four series",
    SWE_Fig5_monthly_change = "Month-by-month change in P, T and SWE after the step year",
    SWE_Fig6_MarApr_SWE_series = "Mar-Apr SWE with period means and Pettitt break",
    SWE_Fig7_observed_vs_predicted_and_residual = "Observed vs P,T-predicted SWE, and the residual with period means",
    SWE_Fig8_pixel_step_vs_snow_regime = "Pixel-level change vs mean SWE, and its map",
    SWE_Fig9_effect_size_by_region = "Cohen's d (95% CI) by region and series")
  have <- c(list.files(tdir, pattern = "\\.csv$"), list.files(fdir, pattern = "\\.png$"))
  man <- data.table::data.table(file = have, kind = ifelse(grepl("\\.csv$", have), "table", "figure"),
                                folder = ifelse(grepl("\\.csv$", have), "tables", "figures"),
                                description = unname(desc[sub("\\.(csv|png)$", "", have)]))
  data.table::fwrite(man[order(kind, file)], .resolve_out(file.path(SWE_REPORT_DIR, "SWE_report_manifest.csv")))
  
  peq <- function(p) ifelse(is.finite(p) & p < 0.001, "p < 0.001", paste0("p = ", .fmt_p(p)))
  # ---- console summary of the two headline results ---------------------------------------------------
  if (nrow(kT)) cat(sprintf("  Nov-Mar temperature: Pettitt %s (break after %d) | change at %d = %s degC (Welch %s)\n",
                            peq(kT$pettitt_p), kT$pettitt_last_year_before, sy, .fmt_sgn(kT$diff, 2), peq(kT$welch_p)))
  if (nrow(kR)) cat(sprintf("  SWE residual (P, T removed): Pettitt %s | change at %d = %s mm (Welch %s) | P,T account for %s of the raw SWE change\n",
                            peq(kR$pettitt_p), sy, .fmt_sgn(kR$diff, 1), peq(kR$welch_p), .fmt_expl(exp_frac)))
  if (!is.null(chk_tbl)) cat(sprintf("  Evidence checklist passed: %d of %d (see Table05a_evidence_checklist.csv)\n", sum(chk_tbl$Passes == "yes"), nrow(chk_tbl)))
  cat(sprintf("  Outputs: %s (%d tables, %d figures; see SWE_report_manifest.csv)\n",
              normalizePath(SWE_REPORT_DIR, winslash = "/", mustWork = FALSE), sum(man$kind == "table"), sum(man$kind == "figure")))
  invisible(list(key = key, key_numbers = key_numbers, nested = nest_all, oos = oos_all, regression = reg_all))
}

####################################################################################
# EVIDENCE OUTPUTS FOR TWO REPORT STATEMENTS (Functions 21-22)
#   Claim 1: the Nov-Mar temperature change is better described as a trend than as a step
#   Claim 2: precipitation matters about as much as temperature for the SWE change
# Reads swe_homogeneity/swe_homogeneity_annual_series.csv and swe_step_tests_by_region.csv
# (written by 9pr_pet_trends.r) and writes tables/figures to swe_homogeneity/evidence/.
# No narrative text is written; tables/Claims_key_numbers.csv lists every number with its source.
####################################################################################
SWE_EVID_DIR      <- file.path(SWE_HOMOG_DIR, "evidence")
N_BOOT_ATTR       <- 2000L   # bootstrap replicates for the attribution (Claim 2)
N_BOOT_MS_BASIN   <- 1000L   # parametric-bootstrap replicates for step-vs-trend, basin
N_BOOT_MS_REGION  <- 200L    # ... and for each sub-basin

.ev_save_tbl <- function(dt, stem) {
  dir.create(file.path(SWE_EVID_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(dt, .resolve_out(file.path(SWE_EVID_DIR, "tables", paste0(stem, ".csv"))))
  invisible(dt)
}
.ev_save_fig <- function(stem, plot, w, h) {
  if (is.null(plot)) return(NULL)
  fdir <- file.path(SWE_EVID_DIR, "figures"); dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  saved <- .save_gg(file.path(fdir, stem), plot, w, h, dpi = 200, formats = "png")
  if (length(saved)) file.path("figures", saved[1]) else NULL
}
.ev_region_levels <- function(ids) { bid <- .swe_basin_id(); c(intersect(bid, ids), setdiff(sort(ids), bid)) }

# Season aggregation across the year boundary (same rule as the pipeline: months >= wrap_from -> next year)
.ev_season <- function(x, years, months, season_months, wrap_from = 11L, fun = mean) {
  keep <- months %in% season_months
  sy <- (years + as.integer(months >= wrap_from))[keep]
  agg <- vapply(split(x[keep], sy), function(v) if (sum(is.finite(v)) == length(season_months)) fun(v) else NA_real_, numeric(1))
  yrs <- sort(unique(years)); out <- agg[as.character(yrs)]; data.frame(year = yrs, v = as.numeric(out))
}

# ------------------------------------------------------------------------------------------------
# Claim 1 - step vs trend
# ------------------------------------------------------------------------------------------------
.ev_gls <- function(f, d) tryCatch(suppressWarnings(nlme::gls(f, data = d, correlation = nlme::corAR1(form = ~ tt),
                                                              method = "ML")), error = function(e) NULL)

# AR(1)-GLS comparison of: constant, linear trend, step at the a-priori year, step at the Pettitt (searched)
# year, and trend + step. The searched step is charged 2 AIC units for having searched the break year.
.ev_model_compare <- function(d, step_year) {
  d <- d[is.finite(d$v), ]; d <- d[order(d$year), ]
  if (nrow(d) < 30L) return(NULL)
  d$tt <- d$year - min(d$year); d$s79 <- as.integer(d$year >= step_year)
  pt <- tryCatch(trend::pettitt.test(d$v), error = function(e) NULL)
  if (is.null(pt)) return(NULL)
  kb <- d$year[as.integer(pt$estimate)[1]] + 1L; d$sb <- as.integer(d$year >= kb)   # [1]: ties in |U| can return 2 positions
  m <- list(constant = .ev_gls(v ~ 1, d), trend = .ev_gls(v ~ tt, d), step_fixed = .ev_gls(v ~ s79, d),
            step_searched = .ev_gls(v ~ sb, d), trend_plus_step = .ev_gls(v ~ tt + s79, d))
  if (any(vapply(m, is.null, logical(1)))) return(NULL)
  aic_raw <- vapply(m, stats::AIC, numeric(1)); aic <- aic_raw; aic[["step_searched"]] <- aic[["step_searched"]] + 2
  dA <- aic - min(aic); w <- exp(-0.5 * dA); w <- w / sum(w)
  pl <- function(a, b) tryCatch(stats::anova(a, b)[["p-value"]][2], error = function(e) NA_real_)
  tt1 <- summary(m$trend)$tTable; tt2 <- summary(m$trend_plus_step)$tTable; ts1 <- summary(m$step_fixed)$tTable
  data.table::data.table(n_years = nrow(d), year_min = min(d$year), year_max = max(d$year),
                         record_midpoint = mean(range(d$year)), pettitt_p = pt$p.value, pettitt_first_year_after = kb,
                         dAIC_constant = dA[["constant"]], dAIC_trend = dA[["trend"]], dAIC_step_fixed = dA[["step_fixed"]],
                         dAIC_step_searched_penalised = dA[["step_searched"]], dAIC_trend_plus_step = dA[["trend_plus_step"]],
                         dAIC_step_searched_unpenalised = aic_raw[["step_searched"]] - min(aic),
                         weight_trend = w[["trend"]], weight_step_fixed = w[["step_fixed"]], weight_step_searched = w[["step_searched"]],
                         weight_trend_plus_step = w[["trend_plus_step"]], weight_constant = w[["constant"]],
                         best_model = names(which.min(dA)),
                         p_step_needed_given_trend = pl(m$trend, m$trend_plus_step), p_trend_needed_given_step = pl(m$step_fixed, m$trend_plus_step),
                         trend_alone_per_yr = tt1["tt", "Value"], trend_alone_p = tt1["tt", "p-value"],
                         step_fixed_alone = ts1["s79", "Value"], step_fixed_alone_p = ts1["s79", "p-value"],
                         trend_in_trend_plus_step = tt2["tt", "Value"], step_in_trend_plus_step = tt2["s79", "Value"],
                         step_in_trend_plus_step_p = tt2["s79", "p-value"],
                         ar1_phi_trend = as.numeric(stats::coef(m$trend$modelStruct$corStruct, unconstrained = FALSE)))
}

.ev_sim_ar1 <- function(mu, phi, sig) {
  n <- length(mu); e <- numeric(n); e[1] <- stats::rnorm(1, 0, sig)
  for (i in 2:n) e[i] <- phi * e[i - 1] + stats::rnorm(1, 0, sig * sqrt(max(1 - phi^2, 1e-8)))
  as.numeric(mu) + e
}
# D = AIC(trend) - AIC(step at the Pettitt year): large D = the step fits better than the trend
.ev_stat_D <- function(y, yrs) {
  dd <- data.frame(v = as.numeric(y), tt = yrs - min(yrs))
  pt <- tryCatch(trend::pettitt.test(dd$v), error = function(e) NULL)
  if (is.null(pt)) return(c(NA_real_, NA_real_, NA_real_))
  kb <- yrs[as.integer(pt$estimate)[1]] + 1L; dd$sb <- as.integer(yrs >= kb)
  mt <- .ev_gls(v ~ tt, dd); ms <- .ev_gls(v ~ sb, dd)
  if (is.null(mt) || is.null(ms)) return(c(NA_real_, NA_real_, NA_real_))
  c(stats::AIC(mt) - stats::AIC(ms), kb, pt$p.value)
}
# Parametric bootstrap: simulate from the fitted trend model and from the fitted step model (both with AR(1)
# noise), apply the SAME procedure (Pettitt search, then AIC comparison) and ask which truth reproduces the data.
.ev_boot_modelsel <- function(d, B, seed = 20260101L) {
  d <- d[is.finite(d$v), ]; d <- d[order(d$year), ]; d$tt <- d$year - min(d$year)
  set.seed(seed)
  obs <- .ev_stat_D(d$v, d$year); if (anyNA(obs)) return(NULL)
  d$sb <- as.integer(d$year >= obs[2])
  mt <- .ev_gls(v ~ tt, d); ms <- .ev_gls(v ~ sb, d); if (is.null(mt) || is.null(ms)) return(NULL)
  gp <- function(m) as.numeric(stats::coef(m$modelStruct$corStruct, unconstrained = FALSE))
  run <- function(m) do.call(rbind, lapply(seq_len(B), function(i) .ev_stat_D(.ev_sim_ar1(stats::fitted(m), gp(m), m$sigma), d$year)))
  St <- run(mt); Ss <- run(ms)
  St <- St[stats::complete.cases(St), , drop = FALSE]; Ss <- Ss[stats::complete.cases(Ss), , drop = FALSE]
  mid <- mean(range(d$year))
  summ <- data.table::data.table(
    n_sims_trend = nrow(St), n_sims_step = nrow(Ss), observed_D = obs[1], observed_break_first_year_after = obs[2],
    observed_pettitt_p = obs[3], record_midpoint = mid,
    p_reject_trend = (1 + sum(St[, 1] >= obs[1])) / (nrow(St) + 1),     # trend truth: P(step looks at least this good)
    p_reject_step  = (1 + sum(Ss[, 1] <= obs[1])) / (nrow(Ss) + 1),     # step truth: P(trend looks at least this good)
    D_median_if_trend = stats::median(St[, 1]), D_median_if_step = stats::median(Ss[, 1]),
    pettitt_p_lt_0.05_if_trend = mean(St[, 3] < 0.05), pettitt_p_as_small_if_trend = mean(St[, 3] <= obs[3]),
    pettitt_p_lt_0.05_if_step = mean(Ss[, 3] < 0.05),
    break_year_median_if_trend = stats::median(St[, 2]), break_year_q25_if_trend = unname(stats::quantile(St[, 2], 0.25)),
    break_year_q75_if_trend = unname(stats::quantile(St[, 2], 0.75)),
    break_within_5yr_of_midpoint_if_trend = mean(abs(St[, 2] - mid) <= 5),
    break_year_median_if_step = stats::median(Ss[, 2]),
    fitted_trend_per_yr = unname(stats::coef(mt)["tt"]), fitted_step_size = unname(stats::coef(ms)["sb"]))
  list(summary = summ, sims_trend = St, sims_step = Ss, observed = obs)
}

.ev_series_set <- function(annual, rid) {
  out <- lapply(c("Temp_NovMar", "Precip_NovMar", "SWE_MarApr"), function(s) {
    d <- .swe_get_series(annual, rid, s); data.frame(year = d$year, v = d$value) })
  names(out) <- c("Temp_NovMar", "Precip_NovMar", "SWE_MarApr"); out
}

# Nov-Mar temperature WITHOUT January (January dominates the winter change): from the monthly means in the metadata
.ev_noJan_series <- function(rid) {
  src <- NULL
  if (!is.null(metadata$regional_avg_monthly) && "tair_degC_month" %in% names(metadata$regional_avg_monthly)) {
    r <- data.table::as.data.table(metadata$regional_avg_monthly); r <- r[region_id == rid]
    if (nrow(r)) src <- data.frame(date = as.Date(r$date), t = r$tair_degC_month)
  }
  if (is.null(src) && rid == .swe_basin_id() && "tair_degC_month" %in% names(basin_avg_monthly))
    src <- data.frame(date = as.Date(basin_avg_monthly$date), t = basin_avg_monthly$tair_degC_month)
  if (is.null(src) || !nrow(src)) return(NULL)
  src <- src[order(src$date), ]
  .ev_season(src$t, as.integer(format(src$date, "%Y")), as.integer(format(src$date, "%m")), c(11L, 12L, 2L, 3L))
}

create_claim1_trend_vs_step_outputs <- function() {
  cat("\n-- Function 21: Claim 1 (temperature: trend vs step) --\n")
  dat <- .swe_load()
  if (is.null(dat$annual)) { cat("  swe_homogeneity_annual_series.csv not found - skipped.\n"); return(invisible(NULL)) }
  annual <- dat$annual; tests <- dat$tests; sy <- .swe_step_year(); bid <- .swe_basin_id()
  rids <- .ev_region_levels(unique(annual$region_id))
  cat("  model comparison (AR(1) GLS)...\n")
  rows <- list(); sets <- list()
  for (r in rids) {
    ss <- .ev_series_set(annual, r)
    nj <- tryCatch(.ev_noJan_series(r), error = function(e) NULL)
    if (!is.null(nj)) ss$Temp_NovMar_noJan <- nj
    sets[[r]] <- ss
    for (s in names(ss)) { x <- .ev_model_compare(ss[[s]], sy)
    if (!is.null(x)) rows[[length(rows) + 1L]] <- cbind(data.table::data.table(region_id = r, series = s), x) }
  }
  mc <- data.table::rbindlist(rows, fill = TRUE)
  if (!nrow(mc)) { cat("  no series long enough - skipped.\n"); return(invisible(NULL)) }
  .ev_save_tbl(mc, "Claim1_Table1_model_comparison")
  
  cat("  parametric bootstrap (this takes a minute)...\n")
  bs <- list(); sims <- list()
  for (r in rids) {
    ss <- sets[[r]]
    for (s in intersect(names(ss), if (r == bid) c("Temp_NovMar", "Precip_NovMar", "SWE_MarApr", "Temp_NovMar_noJan") else "Temp_NovMar")) {
      B <- if (r == bid) (if (s == "Temp_NovMar") N_BOOT_MS_BASIN else max(200L, N_BOOT_MS_BASIN %/% 2L)) else N_BOOT_MS_REGION
      b <- .ev_boot_modelsel(ss[[s]], B)
      if (!is.null(b)) { bs[[length(bs) + 1L]] <- cbind(data.table::data.table(region_id = r, series = s), b$summary); sims[[paste(r, s)]] <- b }
    }
  }
  bst <- data.table::rbindlist(bs, fill = TRUE); .ev_save_tbl(bst, "Claim1_Table2_parametric_bootstrap")
  
  # monthly contributions to the Nov-Mar temperature change (basin)
  mch <- tryCatch(.swe_monthly_change(sy), error = function(e) NULL)
  mtab <- NULL
  if (!is.null(mch) && nrow(mch[grepl("^Temperature", variable)])) {
    tm <- mch[grepl("^Temperature", variable), .(month, month_abb, mean_before, mean_after, change = diff, welch_p)]
    inw <- tm$month %in% c(11, 12, 1, 2, 3)
    tm[, in_NovMar := inw]; tm[, contribution_to_NovMar_mean := ifelse(inw, change / 5, NA_real_)]
    tm[, share_of_NovMar_change_pct := ifelse(inw, 100 * contribution_to_NovMar_mean / sum(contribution_to_NovMar_mean, na.rm = TRUE), NA_real_)]
    yr <- as.integer(format(basin_avg_monthly$date, "%Y")); mo <- as.integer(format(basin_avg_monthly$date, "%m"))
    sl <- do.call(rbind, lapply(1:12, function(m) { z <- data.frame(y = yr[mo == m], v = basin_avg_monthly$tair_degC_month[mo == m]); z <- z[is.finite(z$v), ]
    if (nrow(z) < 20L) return(data.frame(month = m, ols_slope_degC_per_yr = NA_real_, ols_slope_p = NA_real_))
    cs <- summary(stats::lm(v ~ y, data = z))$coefficients; data.frame(month = m, ols_slope_degC_per_yr = cs[2, 1], ols_slope_p = cs[2, 4]) }))
    mtab <- merge(tm, data.table::as.data.table(sl), by = "month"); mtab <- mtab[order(month)]
    .ev_save_tbl(mtab, "Claim1_Table3_monthly_contributions")
  }
  # coefficients when BOTH a trend and a step are allowed (pipeline output): T needs no step, P and SWE do
  if (!is.null(tests)) {
    cols <- intersect(c("region_id", "series", "trend_no_step", "trend_no_step_p", "trend_with_step", "trend_with_step_p",
                        "step_coef", "step_p", "delta_AIC_step_vs_nostep", "ar1_phi_with_step"), names(tests))
    .ev_save_tbl(tests[series %in% c("Temp_NovMar", "Precip_NovMar", "SWE_MarApr"), ..cols][order(series, region_id)], "Claim1_Table4_trend_and_step_coefficients")
  }
  
  cat("  figures...\n")
  # Fig 1: basin Nov-Mar temperature, trend vs steps
  fig1 <- NULL
  d <- sets[[bid]]$Temp_NovMar; mrow <- mc[region_id == bid & series == "Temp_NovMar"]; brow <- bst[region_id == bid & series == "Temp_NovMar"]
  if (!is.null(d) && nrow(mrow)) {
    d <- d[is.finite(d$v), ]; d$tt <- d$year - min(d$year); kb <- mrow$pettitt_first_year_after
    ft <- stats::fitted(.ev_gls(v ~ tt, d)); d$trend <- as.numeric(ft)
    seg <- function(cut, lab) { b <- mean(d$v[d$year < cut]); a <- mean(d$v[d$year >= cut])
    data.frame(x0 = c(min(d$year), cut - 0.5), x1 = c(cut - 0.5, max(d$year)), y = c(b, a), model = lab) }
    sg <- rbind(seg(sy, sprintf("Step at %d (a priori)", sy)), seg(kb, sprintf("Step at %d (Pettitt, searched)", kb)))
    sub <- sprintf("dAIC: trend %.1f | step %d %.1f | step %d (searched, +2) %.1f | trend+step %.1f.  Adding a step to the trend: p = %s; adding a trend to a step: p = %s",
                   mrow$dAIC_trend, sy, mrow$dAIC_step_fixed, kb, mrow$dAIC_step_searched_penalised, mrow$dAIC_trend_plus_step,
                   .fmt_p(mrow$p_step_needed_given_trend), .fmt_p(mrow$p_trend_needed_given_step))
    if (nrow(brow)) sub <- paste0(sub, sprintf("\nBootstrap: pure-step truth rejected at p = %s; pure-trend truth rejected at p = %s",
                                               .fmt_p(brow$p_reject_step), .fmt_p(brow$p_reject_trend)))
    fig1 <- ggplot2::ggplot(d, ggplot2::aes(year, v)) +
      ggplot2::geom_line(colour = "grey60", linewidth = 0.4) + ggplot2::geom_point(colour = "#6a3d9a", size = 1.4) +
      ggplot2::geom_line(ggplot2::aes(y = trend, linetype = "Linear trend"), colour = "#1b7837", linewidth = 1.2) +
      ggplot2::geom_segment(data = sg, ggplot2::aes(x = x0, xend = x1, y = y, yend = y, colour = model), inherit.aes = FALSE, linewidth = 1.1) +
      ggplot2::scale_colour_manual(values = stats::setNames(c("#d73027", "#fc8d59"), unique(sg$model)), name = NULL) +
      ggplot2::scale_linetype_manual(values = c("Linear trend" = "solid"), name = NULL) +
      ggplot2::labs(title = "Nov-Mar temperature (basin): linear trend vs abrupt step", subtitle = sub,
                    x = "Year", y = paste0("Temperature (", DEGC, ")")) + .swe_theme(10)
  }
  p1 <- .ev_save_fig("Claim1_Fig1_trend_vs_step", fig1, 10, 5.6)
  # Fig 2: dAIC by region
  fig2 <- NULL
  l2 <- data.table::melt(mc[series %in% c("Temp_NovMar", "Precip_NovMar", "SWE_MarApr"),
                            .(region_id, series, trend = dAIC_trend, `step (a priori)` = dAIC_step_fixed,
                              `step (searched, +2)` = dAIC_step_searched_penalised, `trend + step` = dAIC_trend_plus_step)],
                         id.vars = c("region_id", "series"), variable.name = "model", value.name = "dAIC")
  if (nrow(l2)) { l2[, region_id := factor(region_id, levels = rev(.ev_region_levels(unique(region_id))))]
    l2[, series := factor(.SWE_SERIES_LAB[as.character(series)], levels = unname(.SWE_SERIES_LAB[c("Temp_NovMar", "Precip_NovMar", "SWE_MarApr")]))]
    fig2 <- ggplot2::ggplot(l2, ggplot2::aes(dAIC, region_id, colour = model, shape = model)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60") + ggplot2::geom_point(size = 2.4, alpha = 0.9) +
      ggplot2::facet_wrap(~series, nrow = 1) +
      ggplot2::scale_colour_manual(values = c(trend = "#1b7837", `step (a priori)` = "#d73027", `step (searched, +2)` = "#fc8d59", `trend + step` = "#4575b4"), name = NULL) +
      ggplot2::scale_shape_manual(values = c(16, 15, 17, 18), name = NULL) +
      ggplot2::labs(title = "Model support by region (AR(1) GLS; dAIC = AIC - best AIC; 0 = best)",
                    subtitle = "Temperature favours the trend; precipitation and SWE favour a step (with or without a trend).", x = "dAIC", y = NULL) + .swe_theme(9) }
  p2 <- .ev_save_fig("Claim1_Fig2_dAIC_by_region", fig2, 13, 4.8)
  # Fig 3: bootstrap
  fig3 <- NULL; bb <- sims[[paste(bid, "Temp_NovMar")]]
  if (!is.null(bb)) {
    dd <- rbind(data.frame(D = bb$sims_trend[, 1], truth = "truth = linear trend"), data.frame(D = bb$sims_step[, 1], truth = "truth = abrupt step"))
    pA <- ggplot2::ggplot(dd, ggplot2::aes(D, fill = truth)) + ggplot2::geom_histogram(bins = 40, alpha = 0.6, position = "identity", colour = NA) +
      ggplot2::geom_vline(xintercept = bb$observed[1], linewidth = 1) +
      ggplot2::scale_fill_manual(values = c("truth = linear trend" = "#1b7837", "truth = abrupt step" = "#d73027"), name = NULL) +
      ggplot2::labs(title = "A) Which truth reproduces the observed fit?", subtitle = "D = AIC(trend) - AIC(step at the Pettitt year); large D = step fits better. Black line = observed.",
                    x = "D", y = "Simulated series") + .swe_theme(9)
    bk <- data.frame(brk = bb$sims_trend[, 2])
    pB <- ggplot2::ggplot(bk, ggplot2::aes(brk)) + ggplot2::geom_histogram(binwidth = 2, fill = "#1b7837", alpha = 0.7, colour = "white") +
      ggplot2::geom_vline(xintercept = bb$observed[2], linewidth = 1) + ggplot2::geom_vline(xintercept = mean(range(d$year)), linetype = "dashed", colour = "grey30") +
      ggplot2::labs(title = "B) Pettitt break year when the truth is a pure trend",
                    subtitle = sprintf("Black = observed (%d); dashed = record midpoint. Pettitt p < 0.05 in %.0f%% of pure-trend simulations.", bb$observed[2], 100 * mean(bb$sims_trend[, 3] < 0.05)),
                    x = "First year of the later regime", y = "Simulated series") + .swe_theme(9)
    fig3 <- pA | pB }
  p3 <- .ev_save_fig("Claim1_Fig3_bootstrap_step_vs_trend", fig3, 13, 4.8)
  # Fig 4: monthly contributions
  fig4 <- NULL
  if (!is.null(mtab)) {
    mtab[, month_abb := factor(month_abb, levels = month.abb)]
    pa <- ggplot2::ggplot(mtab, ggplot2::aes(month_abb, change, fill = in_NovMar)) + ggplot2::geom_col(width = 0.75) + ggplot2::geom_hline(yintercept = 0) +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "#6a3d9a", `FALSE` = "grey70"), labels = c("other months", "Nov-Mar"), name = NULL) +
      ggplot2::labs(title = sprintf("A) Temperature change after %d, by month", sy), x = NULL, y = paste0("Change (", DEGC, ")")) + .swe_theme(9)
    pb <- ggplot2::ggplot(mtab, ggplot2::aes(month_abb, ols_slope_degC_per_yr, fill = in_NovMar)) + ggplot2::geom_col(width = 0.75) + ggplot2::geom_hline(yintercept = 0) +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "#6a3d9a", `FALSE` = "grey70"), guide = "none") +
      ggplot2::labs(title = "B) Linear trend by month", x = NULL, y = paste0("OLS slope (", DEGC, " per year)")) + .swe_theme(9)
    fig4 <- pa / pb }
  p4 <- .ev_save_fig("Claim1_Fig4_monthly_contributions", fig4, 9, 7)
  invisible(list(model_comparison = mc, bootstrap = bst, monthly = mtab, figs = c(p1, p2, p3, p4)))
}

# ------------------------------------------------------------------------------------------------
# Claim 2 - how much of the SWE change is precipitation vs temperature
#   Contribution of X to the SWE change between two periods = beta_X * (mean X after - mean X before),
#   with beta from SWE ~ P + T fitted over the full record. The two contributions plus the residual
#   change add up to the observed change. Uncertainty: stratified (by period) bootstrap of years.
# ------------------------------------------------------------------------------------------------
.ev_attr_once <- function(df, before, after) {
  m <- stats::lm(y ~ P + Tm, data = df); b <- stats::coef(m)
  ib <- df$year %in% before; ia <- df$year %in% after
  dP <- mean(df$P[ia]) - mean(df$P[ib]); dT <- mean(df$Tm[ia]) - mean(df$Tm[ib]); dY <- mean(df$y[ia]) - mean(df$y[ib])
  cP <- unname(b["P"]) * dP; cT <- unname(b["Tm"]) * dT
  c(beta_P = unname(b["P"]), beta_T = unname(b["Tm"]), dP = dP, dT = dT, dSWE = dY, cP = cP, cT = cT, resid = dY - cP - cT)
}
.ev_attr_boot <- function(df, before, after, B) {
  idx_b <- which(df$year %in% before); idx_a <- which(df$year %in% after); idx_o <- setdiff(seq_len(nrow(df)), c(idx_b, idx_a))
  one <- function() { s <- c(sample(idx_b, replace = TRUE), sample(idx_a, replace = TRUE), if (length(idx_o)) sample(idx_o, replace = TRUE))
  d <- df[s, ]; d$year <- c(rep(before[1], length(idx_b)), rep(after[1], length(idx_a)), rep(-1L, length(idx_o)))   # strata labels
  m <- stats::lm(y ~ P + Tm, data = d); b <- stats::coef(m); ib <- d$year == before[1]; ia <- d$year == after[1]
  dP <- mean(d$P[ia]) - mean(d$P[ib]); dT <- mean(d$Tm[ia]) - mean(d$Tm[ib]); dY <- mean(d$y[ia]) - mean(d$y[ib])
  c(cP = unname(b["P"]) * dP, cT = unname(b["Tm"]) * dT, dSWE = dY) }
  M <- t(replicate(B, one())); M
}
.ev_windows <- function(sy) list(
  `1950-78 vs 1979-2025` = list(before = 1950:(sy - 1L), after = sy:2025L),
  `1950-78 vs 1979-1990 (post hoc)` = list(before = 1950:(sy - 1L), after = sy:1990L),
  `1950-78 vs 1991-2020` = list(before = 1950:(sy - 1L), after = 1991:2020),
  `1950-78 vs 2021-2025` = list(before = 1950:(sy - 1L), after = 2021:2025))

# Shapley (LMG) split of R-squared between P and T
.ev_shapley <- function(df) {
  r2 <- function(f) summary(stats::lm(f, data = df))$r.squared
  rp <- r2(y ~ P); rt <- r2(y ~ Tm); rpt <- r2(y ~ P + Tm)
  c(R2_P_alone = rp, R2_T_alone = rt, R2_both = rpt, share_P = 0.5 * (rp + rpt - rt) / rpt, share_T = 0.5 * (rt + rpt - rp) / rpt)
}

create_claim2_attribution_outputs <- function() {
  cat("\n-- Function 22: Claim 2 (precipitation vs temperature) --\n")
  dat <- .swe_load()
  if (is.null(dat$annual)) { cat("  swe_homogeneity_annual_series.csv not found - skipped.\n"); return(invisible(NULL)) }
  annual <- dat$annual; sy <- .swe_step_year(); bid <- .swe_basin_id(); set.seed(20260102)
  rids <- .ev_region_levels(unique(annual$region_id)); wins <- .ev_windows(sy)
  frames <- stats::setNames(lapply(rids, function(r) .swe_model_frame(annual, r)), rids)
  frames <- frames[!vapply(frames, is.null, logical(1))]
  if (!length(frames)) { cat("  no complete P/T/SWE years - skipped.\n"); return(invisible(NULL)) }
  cat("  attribution (stratified bootstrap)...\n")
  rows <- list()
  for (r in names(frames)) {
    df <- frames[[r]]
    for (wn in names(wins)) {
      if (r != bid && wn != names(wins)[1]) next          # sub-basins: main window only
      w <- wins[[wn]]; if (sum(df$year %in% w$before) < 8L || sum(df$year %in% w$after) < 3L) next
      est <- .ev_attr_once(df, w$before, w$after)
      M <- .ev_attr_boot(df, w$before, w$after, N_BOOT_ATTR)
      q <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
      fit <- M[, "cP"] + M[, "cT"]; shP <- 100 * M[, "cP"] / fit; diffPT <- M[, "cP"] - M[, "cT"]
      rows[[length(rows) + 1L]] <- data.table::data.table(region_id = r, window = wn, n_before = sum(df$year %in% w$before), n_after = sum(df$year %in% w$after),
                                                          beta_P_mm_per_mm = est[["beta_P"]], beta_T_mm_per_degC = est[["beta_T"]], dP_mm = est[["dP"]], dT_degC = est[["dT"]],
                                                          dSWE_observed_mm = est[["dSWE"]], dSWE_observed_lo = q(M[, "dSWE"])[1], dSWE_observed_hi = q(M[, "dSWE"])[2],
                                                          contribution_P_mm = est[["cP"]], contribution_P_lo = q(M[, "cP"])[1], contribution_P_hi = q(M[, "cP"])[2],
                                                          contribution_T_mm = est[["cT"]], contribution_T_lo = q(M[, "cT"])[1], contribution_T_hi = q(M[, "cT"])[2],
                                                          residual_change_mm = est[["resid"]], fitted_change_mm = est[["cP"]] + est[["cT"]],
                                                          share_P_of_fitted_pct = 100 * est[["cP"]] / (est[["cP"]] + est[["cT"]]), share_P_lo = q(shP)[1], share_P_hi = q(shP)[2],
                                                          share_T_of_fitted_pct = 100 * est[["cT"]] / (est[["cP"]] + est[["cT"]]),
                                                          contribution_P_minus_T_mm = est[["cP"]] - est[["cT"]], P_minus_T_lo = q(diffPT)[1], P_minus_T_hi = q(diffPT)[2],
                                                          prob_P_at_least_as_large_as_T = mean(abs(M[, "cP"]) >= abs(M[, "cT"]), na.rm = TRUE),
                                                          boot_n = nrow(M))
    }
  }
  at <- data.table::rbindlist(rows, fill = TRUE)
  .ev_save_tbl(at, "Claim2_Table1_attribution_by_window")
  
  # variance (interannual) split + standardised coefficients, with bootstrap CIs
  vs <- data.table::rbindlist(lapply(names(frames), function(r) { df <- frames[[r]]; e <- .ev_shapley(df)
  bs <- t(replicate(N_BOOT_ATTR, .ev_shapley(df[sample(nrow(df), replace = TRUE), ])))
  q <- function(j) stats::quantile(bs[, j], c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  data.table::data.table(region_id = r, n_years = nrow(df), R2_P_alone = e[["R2_P_alone"]], R2_T_alone = e[["R2_T_alone"]], R2_both = e[["R2_both"]],
                         share_of_R2_P = 100 * e[["share_P"]], share_of_R2_P_lo = 100 * q("share_P")[1], share_of_R2_P_hi = 100 * q("share_P")[2],
                         share_of_R2_T = 100 * e[["share_T"]], share_of_R2_T_lo = 100 * q("share_T")[1], share_of_R2_T_hi = 100 * q("share_T")[2],
                         std_beta_P = .swe_tbl_regression(df, r)$std_beta[2], std_beta_T = .swe_tbl_regression(df, r)$std_beta[3]) }))
  .ev_save_tbl(vs, "Claim2_Table2_variance_shares")
  
  # decade-by-decade contributions (basin and regions): anomalies relative to the full-record mean
  dec <- data.table::rbindlist(lapply(names(frames), function(r) { df <- frames[[r]]; m <- stats::lm(y ~ P + Tm, data = df); b <- stats::coef(m)
  x <- data.table::data.table(year = df$year, obs = df$y - mean(df$y), P_part = unname(b["P"]) * (df$P - mean(df$P)),
                              T_part = unname(b["Tm"]) * (df$Tm - mean(df$Tm)))
  x[, resid_part := obs - P_part - T_part]
  x[, decade := paste0(10L * (year %/% 10L), "s")]; x[decade == "2020s", decade := "2020-25"]
  cbind(region_id = r, x[, .(n_years = .N, obs_anomaly_mm = mean(obs), P_part_mm = mean(P_part), T_part_mm = mean(T_part), residual_mm = mean(resid_part)), by = decade]) }))
  .ev_save_tbl(dec, "Claim2_Table3_decadal_contributions")
  
  # trend-based attribution over the whole record (consistent with Claim 1: warming is a trend)
  tr <- data.table::rbindlist(lapply(names(frames), function(r) { df <- frames[[r]]; b <- stats::coef(stats::lm(y ~ P + Tm, data = df)); span <- max(df$year) - min(df$year)
  sl <- function(v) { cs <- summary(stats::lm(v ~ df$year))$coefficients; c(cs[2, 1], cs[2, 4]) }
  sP <- sl(df$P); sT <- sl(df$Tm); sY <- sl(df$y)
  data.table::data.table(region_id = r, span_years = span, P_trend_mm_per_yr = sP[1], P_trend_p = sP[2], T_trend_degC_per_yr = sT[1], T_trend_p = sT[2],
                         SWE_trend_mm_per_yr = sY[1], SWE_trend_p = sY[2], contribution_P_mm = unname(b["P"]) * sP[1] * span, contribution_T_mm = unname(b["Tm"]) * sT[1] * span,
                         SWE_trend_total_mm = sY[1] * span) }))
  tr[, residual_trend_mm := SWE_trend_total_mm - contribution_P_mm - contribution_T_mm]
  .ev_save_tbl(tr, "Claim2_Table4_trend_based_attribution")
  
  cat("  figures...\n")
  lab_r <- function(x) factor(x, levels = rev(.ev_region_levels(unique(x))))
  # Fig 1: basin, fixed split: observed vs P vs T vs residual (with CI)
  fig1 <- NULL; a1 <- at[region_id == bid & window == names(wins)[1]]
  if (nrow(a1)) {
    p1 <- data.frame(part = factor(c("Observed SWE change", "Precipitation part", "Temperature part", "Residual"),
                                   levels = c("Observed SWE change", "Precipitation part", "Temperature part", "Residual")),
                     est = c(a1$dSWE_observed_mm, a1$contribution_P_mm, a1$contribution_T_mm, a1$residual_change_mm),
                     lo = c(a1$dSWE_observed_lo, a1$contribution_P_lo, a1$contribution_T_lo, NA), hi = c(a1$dSWE_observed_hi, a1$contribution_P_hi, a1$contribution_T_hi, NA))
    pA <- ggplot2::ggplot(p1, ggplot2::aes(part, est, fill = part)) + ggplot2::geom_col(width = 0.7) + ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = 0.18, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = c("grey40", "#4575b4", "#d73027", "grey75"), guide = "none") +
      ggplot2::labs(title = sprintf("A) Basin Mar-Apr SWE change, %s", names(wins)[1]),
                    subtitle = sprintf("Precipitation %.0f%% / temperature %.0f%% of the fitted change; P at least as large as T in %.0f%% of bootstrap draws",
                                       a1$share_P_of_fitted_pct, a1$share_T_of_fitted_pct, 100 * a1$prob_P_at_least_as_large_as_T),
                    x = NULL, y = "Change in Mar-Apr SWE (mm w.e.)") + .swe_theme(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
    ar <- at[window == names(wins)[1]]
    lr <- rbind(data.frame(region = ar$region_id, part = "Precipitation", est = ar$contribution_P_mm, lo = ar$contribution_P_lo, hi = ar$contribution_P_hi),
                data.frame(region = ar$region_id, part = "Temperature", est = ar$contribution_T_mm, lo = ar$contribution_T_lo, hi = ar$contribution_T_hi))
    lr$region <- lab_r(lr$region)
    pB <- ggplot2::ggplot(lr, ggplot2::aes(est, region, colour = part)) + ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = lo, xmax = hi), position = ggplot2::position_dodge(width = 0.6), size = 0.35) +
      ggplot2::scale_colour_manual(values = c(Precipitation = "#4575b4", Temperature = "#d73027"), name = NULL) +
      ggplot2::labs(title = "B) Contributions by region (95% bootstrap CI)", x = "Contribution to SWE change (mm w.e.)", y = NULL) + .swe_theme(9)
    fig1 <- pA | pB }
  p1f <- .ev_save_fig("Claim2_Fig1_attribution_basin_and_regions", fig1, 13, 5.2)
  # Fig 2: basin, by comparison window
  fig2 <- NULL; aw <- at[region_id == bid]
  if (nrow(aw) > 1L) {
    lw <- rbind(data.frame(window = aw$window, part = "Precipitation", est = aw$contribution_P_mm, lo = aw$contribution_P_lo, hi = aw$contribution_P_hi),
                data.frame(window = aw$window, part = "Temperature", est = aw$contribution_T_mm, lo = aw$contribution_T_lo, hi = aw$contribution_T_hi))
    lw$window <- factor(lw$window, levels = names(wins))
    fig2 <- ggplot2::ggplot(lw, ggplot2::aes(window, est, fill = part)) + ggplot2::geom_col(position = ggplot2::position_dodge(0.75), width = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), position = ggplot2::position_dodge(0.75), width = 0.2) + ggplot2::geom_hline(yintercept = 0) +
      ggplot2::scale_fill_manual(values = c(Precipitation = "#4575b4", Temperature = "#d73027"), name = NULL) +
      ggplot2::labs(title = "The balance between precipitation and temperature depends on the comparison period",
                    subtitle = "Basin; each bar = beta x (change in the driver) relative to 1950-78; error bars = 95% bootstrap CI",
                    x = NULL, y = "Contribution to SWE change (mm w.e.)") + .swe_theme(9) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1)) }
  p2f <- .ev_save_fig("Claim2_Fig2_attribution_by_window", fig2, 10, 5.2)
  # Fig 3: decadal stack (basin)
  fig3 <- NULL; db <- dec[region_id == bid]
  if (nrow(db)) {
    db[, decade := factor(decade, levels = c("1950s", "1960s", "1970s", "1980s", "1990s", "2000s", "2010s", "2020-25"))]
    ld <- data.table::melt(db[, .(decade, `Precipitation part` = P_part_mm, `Temperature part` = T_part_mm, Residual = residual_mm)], id.vars = "decade", variable.name = "part", value.name = "mm")
    fig3 <- ggplot2::ggplot() + ggplot2::geom_col(data = ld, ggplot2::aes(decade, mm, fill = part), width = 0.75) +
      ggplot2::geom_point(data = db, ggplot2::aes(decade, obs_anomaly_mm), shape = 23, size = 3, fill = "black") + ggplot2::geom_hline(yintercept = 0) +
      ggplot2::scale_fill_manual(values = c("Precipitation part" = "#4575b4", "Temperature part" = "#d73027", Residual = "grey75"), name = NULL) +
      ggplot2::labs(title = "Decade-by-decade: what the Mar-Apr SWE anomaly is made of (basin)",
                    subtitle = "Bars = precipitation, temperature and residual parts of the anomaly (relative to the 1951-2025 mean); diamond = observed anomaly.",
                    x = NULL, y = "Mar-Apr SWE anomaly (mm w.e.)") + .swe_theme(9) }
  p3f <- .ev_save_fig("Claim2_Fig3_decadal_contributions", fig3, 10, 5.2)
  # Fig 4: standardised series with 11-year running means (basin)
  fig4 <- NULL; df <- frames[[bid]]
  if (!is.null(df)) {
    z <- function(v) (v - mean(v)) / stats::sd(v)
    s <- rbind(data.frame(year = df$year, z = z(df$P), what = "Nov-Mar precipitation"), data.frame(year = df$year, z = z(df$Tm), what = "Nov-Mar temperature"),
               data.frame(year = df$year, z = z(df$y), what = "Mar-Apr SWE"))
    s$what <- factor(s$what, levels = c("Nov-Mar precipitation", "Nov-Mar temperature", "Mar-Apr SWE"))
    s <- do.call(rbind, lapply(split(s, s$what), function(g) { g$roll <- as.numeric(zoo::rollmean(g$z, k = 11, fill = NA, align = "center")); g }))
    fig4 <- ggplot2::ggplot(s, ggplot2::aes(year, z, colour = what)) + ggplot2::geom_line(alpha = 0.25, linewidth = 0.4) +
      ggplot2::geom_line(ggplot2::aes(y = roll), linewidth = 1.3, na.rm = TRUE) + ggplot2::geom_vline(xintercept = sy - 0.5, linetype = "dashed", colour = "grey40") +
      ggplot2::scale_colour_manual(values = c("Nov-Mar precipitation" = "#4575b4", "Nov-Mar temperature" = "#d73027", "Mar-Apr SWE" = "#01665e"), name = NULL) +
      ggplot2::labs(title = "Precipitation and SWE dip together around 1979-1990; temperature rises steadily",
                    subtitle = "Standardised annual values (thin) and 11-year running means (thick); dashed = a-priori step year.", x = "Year", y = "z-score") + .swe_theme(9) }
  p4f <- .ev_save_fig("Claim2_Fig4_standardised_series", fig4, 10, 5)
  # Fig 5: variance shares
  fig5 <- NULL
  if (nrow(vs)) {
    lv <- rbind(data.frame(region = vs$region_id, part = "Precipitation", est = vs$share_of_R2_P, lo = vs$share_of_R2_P_lo, hi = vs$share_of_R2_P_hi),
                data.frame(region = vs$region_id, part = "Temperature", est = vs$share_of_R2_T, lo = vs$share_of_R2_T_lo, hi = vs$share_of_R2_T_hi))
    lv$region <- lab_r(lv$region)
    fig5 <- ggplot2::ggplot(lv, ggplot2::aes(est, region, colour = part)) + ggplot2::geom_pointrange(ggplot2::aes(xmin = lo, xmax = hi), position = ggplot2::position_dodge(0.6), size = 0.35) +
      ggplot2::scale_colour_manual(values = c(Precipitation = "#4575b4", Temperature = "#d73027"), name = NULL) +
      ggplot2::labs(title = "Year-to-year variability: share of explained variance (Shapley/LMG, 95% bootstrap CI)",
                    subtitle = "A different question from the 1979 change: precipitation dominates interannual SWE variability.", x = "% of the R-squared of SWE ~ P + T", y = NULL) + .swe_theme(9) }
  p5f <- .ev_save_fig("Claim2_Fig5_variance_shares", fig5, 10, 4.6)
  invisible(list(attribution = at, variance = vs, decadal = dec, trend_based = tr, figs = c(p1f, p2f, p3f, p4f, p5f)))
}

# ------------------------------------------------------------------------------------------------
# Claim 2 - how much of the SWE change is precipitation vs temperature
#   Contribution of X to the SWE change between two periods = beta_X * (mean X after - mean X before),
#   with beta from SWE ~ P + T fitted over the full record. The two contributions plus the residual
#   change add up to the observed change. Uncertainty: stratified (by period) bootstrap of years.
####################################################################################
# FUNCTION 23: SNOW -> SUMMER DEMAND LINK
# Does a low-snow spring precede low SPEI (FAO-56 PM and ERA5-Land PEV) in Jul-Sep?
# Reads the existing basin CSVs: swei_results_seasonal/swei_03_..., spei_results_seasonal/spei_KK_...,
# spei_results_seasonal_era5land/spei_era5_KK_..., (optional) spi_results_seasonal/spi_03_..., sspi_results_seasonal/sspi_01_...
#
# Design notes (what keeps the test honest):
#  * SPEI is built from precipitation and PET only - it contains no snowmelt. A link therefore means that snow
#    and the summer water balance share drivers or persist together (warm/dry springs, circulation); it is NOT
#    evidence that melt-water deficits cause low SPEI.
#  * Windows that overlap the snow predictor (e.g. SPEI-6 in July covers Feb-Jul) share their precipitation by
#    construction. They are reported but flagged, and excluded from the primary (multiple-testing) family.
#  * Both snow and SPEI trend over 1950-2025, which alone creates correlation. Every pair is therefore also tested
#    on linearly detrended series; the detrended correlation is the primary evidence of year-to-year coupling.
#  * p-values are corrected (Benjamini-Hochberg) across the primary family of pairs.
####################################################################################
LINK_PRED_MONTH  <- 4L          # April = last month of the SWEI-3 window (Feb-Apr)
LINK_SWEI_SCALE  <- 3L
LINK_SPEI_SCALES <- c(1L, 3L, 6L)
LINK_RESP_MONTHS <- 7:9
LINK_THRESHOLDS  <- c(-0.5, -1)
LINK_REF_YEARS   <- c(1991L, 2020L)   # reference for the standardised Mar-Apr SWE predictor
SWE_LINK_DIR     <- file.path(SWE_HOMOG_DIR, "link_to_summer")

.lk_save_tbl <- function(dt, stem) {
  dir.create(file.path(SWE_LINK_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(dt, .resolve_out(file.path(SWE_LINK_DIR, "tables", paste0(stem, ".csv")))); invisible(dt)
}
.lk_save_fig <- function(stem, plot, w, h) {
  if (is.null(plot)) return(NULL)
  fdir <- file.path(SWE_LINK_DIR, "figures"); dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  saved <- .save_gg(file.path(fdir, stem), plot, w, h, dpi = 200, formats = "png"); if (length(saved)) file.path("figures", saved[1]) else NULL
}
# one calendar month of a basin CSV -> data.frame(year, value)
.lk_month_series <- function(index, scale, month, dir = NULL) {
  x <- tryCatch(suppressWarnings(load_basin_avg_csv(index, scale, data_dir = dir)), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  yr <- as.integer(format(x$date, "%Y")); mo <- as.integer(format(x$date, "%m"))
  d <- data.frame(year = yr[mo == month], value = x$value[mo == month]); d <- d[is.finite(d$value), ]
  if (nrow(d) < 20L) NULL else d
}
.lk_window <- function(month, k) ((month - k + 1L):month - 1L) %% 12L + 1L
.lk_win_label <- function(month, k) { w <- .lk_window(month, min(k, 12L)); if (k >= 12L) "12+ months" else paste0(month.abb[w[1]], if (length(w) > 1L) paste0("-", month.abb[w[length(w)]]) else "") }
.lk_detrend <- function(v, yr) stats::residuals(stats::lm(v ~ yr))

.lk_pair <- function(x, y, yrs) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; yrs <- yrs[ok]; n <- length(x)
  if (n < 15L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NULL)
  ct <- stats::cor.test(x, y); sp <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  xd <- .lk_detrend(x, yrs); yd <- .lk_detrend(y, yrs); cd <- stats::cor.test(xd, yd)
  data.table::data.table(n = n, r_raw = unname(ct$estimate), r_raw_lo = ct$conf.int[1], r_raw_hi = ct$conf.int[2], p_raw = ct$p.value,
                         rho_spearman = unname(sp$estimate), p_spearman = sp$p.value,
                         r_detrended = unname(cd$estimate), r_detrended_lo = cd$conf.int[1], r_detrended_hi = cd$conf.int[2], p_detrended = cd$p.value)
}
.lk_partial <- function(x, y, z, yrs) {      # r(x, y | z, linear time)
  ok <- is.finite(x) & is.finite(y) & is.finite(z); x <- x[ok]; y <- y[ok]; z <- z[ok]; yrs <- yrs[ok]; n <- length(x)
  if (n < 20L) return(NULL)
  rx <- stats::residuals(stats::lm(x ~ z + yrs)); ry <- stats::residuals(stats::lm(y ~ z + yrs)); r <- stats::cor(rx, ry)
  tt <- r * sqrt((n - 4) / (1 - r^2)); p <- 2 * stats::pt(-abs(tt), n - 4)
  data.table::data.table(n = n, partial_r = r, partial_p = p)
}

.lk_inputs <- function() {
  swei <- .lk_month_series("swei", LINK_SWEI_SCALE, LINK_PRED_MONTH)
  sspi <- .lk_month_series("sspi", 1L, LINK_PRED_MONTH, dir = file.path(WD_PATH, "sspi_results_seasonal"))
  swe_z <- NULL
  dat <- tryCatch(.swe_load(), error = function(e) NULL)
  if (!is.null(dat$annual)) {
    s <- .swe_get_series(dat$annual, .swe_basin_id(), "SWE_MarApr")
    if (nrow(s) >= 30L) { ref <- s$value[s$year >= LINK_REF_YEARS[1] & s$year <= LINK_REF_YEARS[2]]
    swe_z <- data.frame(year = s$year, value = (s$value - mean(ref)) / stats::sd(ref)) }
  }
  pred <- list(
    `SWEI-3 Apr (Feb-Apr)` = list(d = swei, win = .lk_window(LINK_PRED_MONTH, LINK_SWEI_SCALE)),
    `Mar-Apr SWE z (ref 1991-2020)` = list(d = swe_z, win = 3:4),
    `SSPI-1 Apr` = list(d = sspi, win = LINK_PRED_MONTH))
  pred <- pred[!vapply(pred, function(p) is.null(p$d), logical(1))]
  resp <- list()
  for (ty in c("spei", "spei_era5")) for (k in LINK_SPEI_SCALES) for (m in LINK_RESP_MONTHS) {
    d <- .lk_month_series(ty, k, m); if (!is.null(d)) resp[[paste(ty, k, m, sep = "|")]] <- list(d = d, type = ty, k = k, m = m) }
  ctl <- list(`SPI-3 Apr` = .lk_month_series("spi", 3L, LINK_PRED_MONTH),
              spei = .lk_month_series("spei", 3L, LINK_PRED_MONTH), spei_era5 = .lk_month_series("spei_era5", 3L, LINK_PRED_MONTH))
  list(pred = pred, resp = resp, ctl = ctl)
}

.lk_demand_table <- function() {
  vars <- intersect(c("Temperature", "PET", "PET_ERA5", "PM_VPD", "PM_NetRadiation", "PM_Wind2m"), unique(all_results$variable))
  unit <- c(Temperature = "degC/yr", PET = "mm/yr (calendar-year total)", PET_ERA5 = "mm/yr (calendar-year total)", PM_VPD = "kPa/yr", PM_NetRadiation = "MJ/m2/d per yr", PM_Wind2m = "m/s per yr")
  rows <- list()
  for (v in vars) for (pp in list(list("annual", NA_integer_), list("monthly", 7L), list("monthly", 8L), list("monthly", 9L))) {
    dd <- if (pp[[1]] == "annual") all_results[variable == v & period == "annual" & !is_basin_average] else
      all_results[variable == v & period == "monthly" & month == pp[[2]] & !is_basin_average]
    dd <- dd[is.finite(p_value_tfpw) & is.finite(sl_tfpw)]
    dd <- .add_dual_fdr(data.table::copy(dd))
    if (!nrow(dd)) next
    bb <- if (pp[[1]] == "annual") all_results[variable == v & period == "annual" & is_basin_average] else
      all_results[variable == v & period == "monthly" & month == pp[[2]] & is_basin_average]
    pf <- stats::p.adjust(dd$p_value_tfpw, "BH")
    rows[[length(rows) + 1L]] <- data.table::data.table(variable = v, unit = unit[[v]], period = if (pp[[1]] == "annual") "annual" else month.abb[pp[[2]]],
                                                        n_pixels = nrow(dd), pct_increasing_sig = 100 * mean(dd$p_value_tfpw < 0.05 & dd$sl_tfpw > 0), pct_decreasing_sig = 100 * mean(dd$p_value_tfpw < 0.05 & dd$sl_tfpw < 0),
                                                        pct_increasing_sig_BH = 100 * mean(pf < 0.05 & dd$sl_tfpw > 0),
                                                        pct_increasing_sig_BH_vc = 100 * mean(dd$sig_fdr_vc & dd$sl_tfpw > 0),
                                                        pct_increasing_sig_BH_both = 100 * mean(dd$sig_both & dd$sl_tfpw > 0),
                                                        pct_decreasing_sig_BH_both = 100 * mean(dd$sig_both & dd$sl_tfpw < 0), median_slope = stats::median(dd$sl_tfpw),
                                                        basin_slope = if (nrow(bb)) bb$sl_tfpw[1] else NA_real_, basin_p = if (nrow(bb)) bb$p_value_tfpw[1] else NA_real_)
  }
  data.table::rbindlist(rows)
}

create_snow_summer_link_outputs <- function() {
  cat("\n-- Function 23: snow -> summer water-balance link --\n")
  inp <- tryCatch(.lk_inputs(), error = function(e) { cat("  could not read the basin CSVs:", conditionMessage(e), "\n"); NULL })
  if (is.null(inp) || !length(inp$pred) || !length(inp$resp)) {
    cat("  Missing inputs. Expected (relative to the working directory):\n    swei_results_seasonal/swei_03_basin_averaged_by_month.csv\n",
        "   spei_results_seasonal/spei_KK_basin_averaged_by_month.csv  and  spei_results_seasonal_era5land/spei_era5_KK_basin_averaged_by_month.csv (KK = 01, 03, 06)\n", sep = "")
    return(invisible(NULL))
  }
  set.seed(20260103)
  tlab <- c(spei = "SPEI (FAO-56 PM)", spei_era5 = "SPEI (ERA5-Land PEV)")
  rows <- list()
  for (pn in names(inp$pred)) { P <- inp$pred[[pn]]
  for (rk in names(inp$resp)) { R <- inp$resp[[rk]]
  j <- merge(P$d, R$d, by = "year", suffixes = c("_x", "_y")); if (nrow(j) < 15L) next
  st <- .lk_pair(j$value_x, j$value_y, j$year); if (is.null(st)) next
  ov <- length(intersect(P$win, .lk_window(R$m, R$k))) > 0 || R$k >= 12L
  rows[[length(rows) + 1L]] <- cbind(data.table::data.table(predictor = pn, response = tlab[[R$type]], response_type = R$type, spei_scale = R$k, response_month = month.abb[R$m],
                                                            response_window = .lk_win_label(R$m, R$k), lead_months_from_predictor = R$m - LINK_PRED_MONTH, overlaps_predictor_window = ov), st) } }
  cor_tbl <- data.table::rbindlist(rows, fill = TRUE)
  if (!nrow(cor_tbl)) { cat("  no overlapping years between predictor and response series - skipped.\n"); return(invisible(NULL)) }
  cor_tbl[, p_detrended_BH := NA_real_]
  cor_tbl[overlaps_predictor_window == FALSE, p_detrended_BH := stats::p.adjust(p_detrended, "BH"), by = predictor]
  cor_tbl[, primary_family := !overlaps_predictor_window]
  .lk_save_tbl(cor_tbl, "Link_Table1_correlations")
  
  # which years / files fed the analysis
  src <- data.table::data.table(item = c(names(inp$pred), paste0("response: ", names(inp$resp)))[seq_len(length(inp$pred) + length(inp$resp))],
                                n_years = c(vapply(inp$pred, function(p) nrow(p$d), numeric(1)), vapply(inp$resp, function(r) nrow(r$d), numeric(1))),
                                first_year = c(vapply(inp$pred, function(p) min(p$d$year), numeric(1)), vapply(inp$resp, function(r) min(r$d$year), numeric(1))),
                                last_year = c(vapply(inp$pred, function(p) max(p$d$year), numeric(1)), vapply(inp$resp, function(r) max(r$d$year), numeric(1))))
  .lk_save_tbl(src, "Link_Table0_inputs_used")
  
  # era split (non-overlapping pairs)
  eras <- list(`1950-1978` = c(1950L, 1978L), `1979-2025` = c(1979L, 2025L), `1991-2020` = c(1991L, 2020L))
  er <- list()
  for (pn in names(inp$pred)) for (rk in names(inp$resp)) { P <- inp$pred[[pn]]; R <- inp$resp[[rk]]
  if (length(intersect(P$win, .lk_window(R$m, R$k))) > 0 || R$k >= 12L) next
  j <- merge(P$d, R$d, by = "year", suffixes = c("_x", "_y"))
  for (en in names(eras)) { jj <- j[j$year >= eras[[en]][1] & j$year <= eras[[en]][2], ]
  st <- .lk_pair(jj$value_x, jj$value_y, jj$year); if (is.null(st)) next
  er[[length(er) + 1L]] <- cbind(data.table::data.table(predictor = pn, response = tlab[[R$type]], spei_scale = R$k, response_month = month.abb[R$m], era = en), st[, .(n, r_raw, p_raw, r_detrended, p_detrended)]) } }
  era_tbl <- data.table::rbindlist(er, fill = TRUE); if (nrow(era_tbl)) .lk_save_tbl(era_tbl, "Link_Table3_by_era")
  
  # joint exceedance (low snow & low SPEI)
  jx <- list()
  for (pn in names(inp$pred)) for (rk in names(inp$resp)) { P <- inp$pred[[pn]]; R <- inp$resp[[rk]]
  if (length(intersect(P$win, .lk_window(R$m, R$k))) > 0 || R$k >= 12L) next
  j <- merge(P$d, R$d, by = "year", suffixes = c("_x", "_y")); j <- j[stats::complete.cases(j), ]
  for (th in LINK_THRESHOLDS) { lx <- j$value_x <= th; ly <- j$value_y <= th
  a <- sum(lx & ly); b <- sum(lx & !ly); c_ <- sum(!lx & ly); d_ <- sum(!lx & !ly)
  if ((a + b) < 3L) next
  fp <- tryCatch(stats::fisher.test(matrix(c(a, c_, b, d_), 2), alternative = "greater")$p.value, error = function(e) NA_real_)
  jx[[length(jx) + 1L]] <- data.table::data.table(predictor = pn, response = tlab[[R$type]], spei_scale = R$k, response_month = month.abb[R$m], threshold = th,
                                                  years_predictor_low = a + b, years_response_low = a + c_, both_low = a, P_resp_low_given_pred_low = a / (a + b), P_resp_low_given_pred_not_low = c_ / (c_ + d_),
                                                  P_resp_low_unconditional = (a + c_) / (a + b + c_ + d_), relative_risk = (a / (a + b)) / max(c_ / (c_ + d_), 1e-9), fisher_p_one_sided = fp, n = a + b + c_ + d_) } }
  jx_tbl <- data.table::rbindlist(jx, fill = TRUE); if (nrow(jx_tbl)) .lk_save_tbl(jx_tbl, "Link_Table2_joint_exceedance")
  
  # does snow add information beyond the spring precipitation / spring water balance?
  pr <- list()
  for (pn in names(inp$pred)[1]) for (rk in names(inp$resp)) { P <- inp$pred[[pn]]; R <- inp$resp[[rk]]
  if (length(intersect(P$win, .lk_window(R$m, R$k))) > 0 || R$k >= 12L) next
  for (cn in c("SPI-3 Apr", R$type)) { C <- inp$ctl[[cn]]; if (is.null(C)) next
  j <- Reduce(function(a, b) merge(a, b, by = "year"), list(setNames(P$d, c("year", "x")), setNames(R$d, c("year", "y")), setNames(C, c("year", "z"))))
  st <- .lk_partial(j$x, j$y, j$z, j$year); if (is.null(st)) next
  pr[[length(pr) + 1L]] <- cbind(data.table::data.table(predictor = pn, response = tlab[[R$type]], spei_scale = R$k, response_month = month.abb[R$m],
                                                        controlling_for = if (cn == "SPI-3 Apr") "SPI-3 Apr (Feb-Apr precipitation) + linear time" else paste0(tlab[[R$type]], "-3 Apr + linear time")), st) } }
  pr_tbl <- data.table::rbindlist(pr, fill = TRUE); if (nrow(pr_tbl)) .lk_save_tbl(pr_tbl, "Link_Table5_partial_correlations")
  
  dem <- tryCatch(.lk_demand_table(), error = function(e) NULL); if (!is.null(dem) && nrow(dem)) .lk_save_tbl(dem, "Link_Table4_demand_trends")
  
  cat("  figures...\n")
  # Fig 1: heatmap of detrended correlations
  ct <- data.table::copy(cor_tbl)
  ct[, `:=`(month_f = factor(response_month, levels = month.abb[LINK_RESP_MONTHS]), scale_f = factor(spei_scale, levels = rev(LINK_SPEI_SCALES)),
            lab = paste0(sprintf("%+.2f%s", r_detrended, ifelse(!overlaps_predictor_window & is.finite(p_detrended_BH) & p_detrended_BH < 0.05, "*", "")),
                         ifelse(overlaps_predictor_window, "\n(overlap)", "")))]
  fig1 <- ggplot2::ggplot(ct, ggplot2::aes(month_f, scale_f, fill = r_detrended, alpha = !overlaps_predictor_window)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = lab), size = 3.2, alpha = 1) +
    ggplot2::facet_grid(response ~ predictor) +
    ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, limits = c(-0.6, 0.6), oob = scales::squish, name = "r (detrended)") +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    ggplot2::labs(title = "Does a low-snow spring precede a dry summer?  Correlation of spring snow with Jul-Sep SPEI",
                  subtitle = "Detrended Pearson r. * = BH-adjusted p < 0.05 within the non-overlapping family. Faded tiles: the SPEI window overlaps the snow window (shared precipitation).",
                  x = "SPEI end month", y = "SPEI timescale (months)") + .swe_theme(9) + ggplot2::theme(legend.position = "right")
  f1 <- .lk_save_fig("Link_Fig1_correlation_heatmap", fig1, 12, 5.8)
  # Fig 2: scatter (SPEI-3 Jul-Sep vs the first predictor)
  fig2 <- NULL; P1 <- inp$pred[[1]]; pn1 <- names(inp$pred)[1]; sc <- list()
  for (ty in c("spei", "spei_era5")) for (m in LINK_RESP_MONTHS) { R <- inp$resp[[paste(ty, 3L, m, sep = "|")]]; if (is.null(R)) next
  j <- merge(P1$d, R$d, by = "year", suffixes = c("_x", "_y")); j$response <- tlab[[ty]]; j$month <- factor(month.abb[m], levels = month.abb[LINK_RESP_MONTHS]); sc[[length(sc) + 1L]] <- j }
  if (length(sc)) { sc <- do.call(rbind, sc); sc$era <- ifelse(sc$year < 1979, "1950-1978", "1979-2025")
  an <- do.call(rbind, lapply(split(sc, list(sc$response, sc$month), drop = TRUE), function(g) { s <- .lk_pair(g$value_x, g$value_y, g$year)
  data.frame(response = g$response[1], month = g$month[1], txt = if (is.null(s)) "" else sprintf("r = %+.2f (detrended %+.2f)", s$r_raw, s$r_detrended)) }))
  fig2 <- ggplot2::ggplot(sc, ggplot2::aes(value_x, value_y)) + ggplot2::geom_hline(yintercept = 0, colour = "grey85") + ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_point(ggplot2::aes(colour = era), size = 1.8, alpha = 0.85) + ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.6, fill = "grey80") +
    ggplot2::geom_text(data = an, ggplot2::aes(x = -Inf, y = Inf, label = txt), hjust = -0.05, vjust = 1.4, size = 2.9, inherit.aes = FALSE) +
    ggplot2::facet_grid(response ~ month) +
    ggplot2::scale_colour_manual(values = c("1950-1978" = "#2166ac", "1979-2025" = "#b2182b"), name = NULL) +
    ggplot2::labs(title = sprintf("%s vs SPEI-3 in Jul, Aug and Sep (windows May-Jul, Jun-Aug, Jul-Sep)", pn1), x = pn1, y = "SPEI-3") + .swe_theme(9) }
  f2 <- .lk_save_fig("Link_Fig2_scatter_snow_vs_SPEI3", fig2, 11, 6)
  # Fig 3: time series
  fig3 <- NULL
  Rp <- inp$resp[["spei|3|8"]]; Re <- inp$resp[["spei_era5|3|8"]]
  if (!is.null(Rp) || !is.null(Re)) {
    s1 <- data.frame(year = P1$d$year, v = P1$d$value)
    ln <- rbind(if (!is.null(Rp)) data.frame(year = Rp$d$year, v = Rp$d$value, what = "SPEI-3 Aug, PM"), if (!is.null(Re)) data.frame(year = Re$d$year, v = Re$d$value, what = "SPEI-3 Aug, ERA5 PEV"))
    fig3 <- ggplot2::ggplot() + ggplot2::geom_col(data = s1, ggplot2::aes(year, v), fill = "#01665e", alpha = 0.55, width = 0.8) +
      ggplot2::geom_line(data = ln, ggplot2::aes(year, v, colour = what), linewidth = 0.7) + ggplot2::geom_hline(yintercept = c(-1, 0), linetype = c("dashed", "solid"), colour = "grey40") +
      ggplot2::scale_colour_manual(values = c("SPEI-3 Aug, PM" = "#d73027", "SPEI-3 Aug, ERA5 PEV" = "#e08214"), name = NULL) +
      ggplot2::labs(title = sprintf("%s (bars) and SPEI-3 for August (lines)", pn1), subtitle = "Dashed line = -1.", x = "Year", y = "Standardised index") + .swe_theme(9) }
  f3 <- .lk_save_fig("Link_Fig3_timeseries_snow_and_SPEI", fig3, 12, 4.8)
  prim <- cor_tbl[primary_family == TRUE]
  if (nrow(prim)) { b <- prim[which.max(abs(r_detrended))]
  cat(sprintf("  non-overlapping pairs tested: %d | BH-significant (detrended): %d | strongest: %s vs %s %s-%s (r_detrended = %+.2f, BH p = %s)\n", nrow(prim),
              sum(prim$p_detrended_BH < 0.05, na.rm = TRUE), b$predictor, b$response, b$spei_scale, b$response_month, b$r_detrended, .fmt_p(b$p_detrended_BH))) }
  cat(sprintf("  Outputs: %s\n", normalizePath(SWE_LINK_DIR, winslash = "/", mustWork = FALSE)))
  invisible(list(correlations = cor_tbl, by_era = era_tbl, joint = jx_tbl, partial = pr_tbl, demand = dem))
}

####################################################################################
# KEY NUMBERS + MANIFEST for the evidence outputs (Functions 21-23)
####################################################################################
create_claims_key_numbers <- function() {
  cat("\n-- Evidence key numbers and manifest --\n")
  tdir <- file.path(SWE_EVID_DIR, "tables"); bid <- .swe_basin_id(); sy <- .swe_step_year()
  rd <- function(n) { f <- file.path(tdir, paste0(n, ".csv")); if (file.exists(f)) data.table::fread(f) else NULL }
  t1 <- rd("Claim1_Table1_model_comparison"); t2 <- rd("Claim1_Table2_parametric_bootstrap"); t3 <- rd("Claim1_Table3_monthly_contributions")
  c1 <- rd("Claim2_Table1_attribution_by_window"); c2 <- rd("Claim2_Table2_variance_shares"); c4 <- rd("Claim2_Table4_trend_based_attribution")
  kn <- function(claim, item, value, unit, source) data.table::data.table(claim = claim, item = item, value = suppressWarnings(as.numeric(value)), unit = unit, source = source)
  R <- list()
  add <- function(x) if (!is.null(x) && nrow(x)) R[[length(R) + 1L]] <<- x
  if (!is.null(t1)) {
    m <- t1[region_id == bid & series == "Temp_NovMar"]; f1 <- "Claim1_Table1_model_comparison.csv"
    if (nrow(m)) add(data.table::rbindlist(list(
      kn("1 temperature trend vs step", "Pettitt p (searches all years)", m$pettitt_p, "", paste0(f1, ": pettitt_p")),
      kn("1 temperature trend vs step", "Pettitt break: first year of later regime", m$pettitt_first_year_after, "year", paste0(f1, ": pettitt_first_year_after")),
      kn("1 temperature trend vs step", "Midpoint of the record", m$record_midpoint, "year", paste0(f1, ": record_midpoint")),
      kn("1 temperature trend vs step", "dAIC linear trend", m$dAIC_trend, "AIC units", paste0(f1, ": dAIC_trend")),
      kn("1 temperature trend vs step", "dAIC step at the a-priori year", m$dAIC_step_fixed, "AIC units", paste0(f1, ": dAIC_step_fixed")),
      kn("1 temperature trend vs step", "dAIC step at the Pettitt year (+2 for the search)", m$dAIC_step_searched_penalised, "AIC units", paste0(f1, ": dAIC_step_searched_penalised")),
      kn("1 temperature trend vs step", "dAIC trend + step", m$dAIC_trend_plus_step, "AIC units", paste0(f1, ": dAIC_trend_plus_step")),
      kn("1 temperature trend vs step", "Akaike weight, linear trend", m$weight_trend, "", paste0(f1, ": weight_trend")),
      kn("1 temperature trend vs step", "p: is a step needed once a trend is in the model", m$p_step_needed_given_trend, "", paste0(f1, ": p_step_needed_given_trend")),
      kn("1 temperature trend vs step", "p: is a trend needed once a step is in the model", m$p_trend_needed_given_step, "", paste0(f1, ": p_trend_needed_given_step")),
      kn("1 temperature trend vs step", "Linear trend", m$trend_alone_per_yr, "degC/yr", paste0(f1, ": trend_alone_per_yr")),
      kn("1 temperature trend vs step", "p of the linear trend", m$trend_alone_p, "", paste0(f1, ": trend_alone_p")),
      kn("1 temperature trend vs step", "Step coefficient when a trend is also in the model", m$step_in_trend_plus_step, "degC", paste0(f1, ": step_in_trend_plus_step")),
      kn("1 temperature trend vs step", "p of that step coefficient", m$step_in_trend_plus_step_p, "", paste0(f1, ": step_in_trend_plus_step_p")))))
    tt <- t1[series == "Temp_NovMar"]
    add(data.table::rbindlist(list(
      kn("1 temperature trend vs step", "Regions (of those tested) where the linear trend has the lowest AIC", sum(tt$best_model == "trend"), "regions", f1),
      kn("1 temperature trend vs step", "Regions tested", nrow(tt), "regions", f1),
      kn("1 temperature trend vs step", "Regions where adding a step to the trend is significant (p < 0.05)", sum(tt$p_step_needed_given_trend < 0.05, na.rm = TRUE), "regions", f1))))
    for (s in c("Precip_NovMar", "SWE_MarApr")) { x <- t1[region_id == bid & series == s]
    if (nrow(x)) add(data.table::rbindlist(list(
      kn("1 contrast: precipitation / SWE", paste(s, "- p: step needed given trend"), x$p_step_needed_given_trend, "", paste0(f1, ": p_step_needed_given_trend")),
      kn("1 contrast: precipitation / SWE", paste(s, "- step coefficient with a trend in the model"), x$step_in_trend_plus_step, "", paste0(f1, ": step_in_trend_plus_step")),
      kn("1 contrast: precipitation / SWE", paste(s, "- dAIC trend"), x$dAIC_trend, "AIC units", paste0(f1, ": dAIC_trend")),
      kn("1 contrast: precipitation / SWE", paste(s, "- dAIC step (a priori)"), x$dAIC_step_fixed, "AIC units", paste0(f1, ": dAIC_step_fixed"))))) }
    nj <- t1[region_id == bid & series == "Temp_NovMar_noJan"]
    if (nrow(nj)) add(data.table::rbindlist(list(
      kn("1 temperature trend vs step", "Without January: dAIC trend", nj$dAIC_trend, "AIC units", paste0(f1, ": dAIC_trend (series Temp_NovMar_noJan)")),
      kn("1 temperature trend vs step", "Without January: dAIC step (searched, +2)", nj$dAIC_step_searched_penalised, "AIC units", paste0(f1, ": dAIC_step_searched_penalised")),
      kn("1 temperature trend vs step", "Without January: Pettitt p", nj$pettitt_p, "", paste0(f1, ": pettitt_p")))))
  }
  if (!is.null(t2)) { b <- t2[region_id == bid & series == "Temp_NovMar"]; f2 <- "Claim1_Table2_parametric_bootstrap.csv"
  if (nrow(b)) add(data.table::rbindlist(list(
    kn("1 temperature trend vs step", "Bootstrap: p-value for rejecting a pure linear trend", b$p_reject_trend, "", paste0(f2, ": p_reject_trend")),
    kn("1 temperature trend vs step", "Bootstrap: p-value for rejecting a pure abrupt step", b$p_reject_step, "", paste0(f2, ": p_reject_step")),
    kn("1 temperature trend vs step", "Share of pure-trend simulations in which Pettitt p < 0.05", b$`pettitt_p_lt_0.05_if_trend`, "fraction", paste0(f2, ": pettitt_p_lt_0.05_if_trend")),
    kn("1 temperature trend vs step", "Share of pure-trend simulations with Pettitt p as small as observed", b$pettitt_p_as_small_if_trend, "fraction", paste0(f2, ": pettitt_p_as_small_if_trend")),
    kn("1 temperature trend vs step", "Median Pettitt break year in pure-trend simulations", b$break_year_median_if_trend, "year", paste0(f2, ": break_year_median_if_trend")),
    kn("1 temperature trend vs step", "Share of pure-trend breaks within 5 years of the record midpoint", b$break_within_5yr_of_midpoint_if_trend, "fraction", paste0(f2, ": break_within_5yr_of_midpoint_if_trend")),
    kn("1 temperature trend vs step", "Bootstrap replicates (trend truth)", b$n_sims_trend, "n", paste0(f2, ": n_sims_trend")))))
  }
  if (!is.null(t3)) { jn <- t3[month == 1]; if (nrow(jn)) add(data.table::rbindlist(list(
    kn("1 temperature trend vs step", "January's share of the Nov-Mar temperature change", jn$share_of_NovMar_change_pct, "%", "Claim1_Table3_monthly_contributions.csv: share_of_NovMar_change_pct"),
    kn("1 temperature trend vs step", "January temperature change", jn$change, "degC", "Claim1_Table3_monthly_contributions.csv: change"),
    kn("1 temperature trend vs step", "January linear trend", jn$ols_slope_degC_per_yr, "degC/yr", "Claim1_Table3_monthly_contributions.csv: ols_slope_degC_per_yr")))) }
  if (!is.null(c1)) for (wn in unique(c1[region_id == bid]$window)) { a <- c1[region_id == bid & window == wn]; f <- "Claim2_Table1_attribution_by_window.csv"; cl <- paste0("2 attribution, ", wn)
  add(data.table::rbindlist(list(
    kn(cl, "Observed SWE change", a$dSWE_observed_mm, "mm w.e.", paste0(f, ": dSWE_observed_mm")),
    kn(cl, "Precipitation change", a$dP_mm, "mm", paste0(f, ": dP_mm")), kn(cl, "Temperature change", a$dT_degC, "degC", paste0(f, ": dT_degC")),
    kn(cl, "Regression coefficient, precipitation", a$beta_P_mm_per_mm, "mm SWE per mm", paste0(f, ": beta_P_mm_per_mm")),
    kn(cl, "Regression coefficient, temperature", a$beta_T_mm_per_degC, "mm SWE per degC", paste0(f, ": beta_T_mm_per_degC")),
    kn(cl, "Precipitation contribution", a$contribution_P_mm, "mm w.e.", paste0(f, ": contribution_P_mm")),
    kn(cl, "Precipitation contribution, 95% CI low", a$contribution_P_lo, "mm w.e.", paste0(f, ": contribution_P_lo")),
    kn(cl, "Precipitation contribution, 95% CI high", a$contribution_P_hi, "mm w.e.", paste0(f, ": contribution_P_hi")),
    kn(cl, "Temperature contribution", a$contribution_T_mm, "mm w.e.", paste0(f, ": contribution_T_mm")),
    kn(cl, "Temperature contribution, 95% CI low", a$contribution_T_lo, "mm w.e.", paste0(f, ": contribution_T_lo")),
    kn(cl, "Temperature contribution, 95% CI high", a$contribution_T_hi, "mm w.e.", paste0(f, ": contribution_T_hi")),
    kn(cl, "Precipitation share of the fitted change", a$share_P_of_fitted_pct, "%", paste0(f, ": share_P_of_fitted_pct")),
    kn(cl, "Precipitation share, 95% CI low", a$share_P_lo, "%", paste0(f, ": share_P_lo")), kn(cl, "Precipitation share, 95% CI high", a$share_P_hi, "%", paste0(f, ": share_P_hi")),
    kn(cl, "Probability precipitation part >= temperature part (bootstrap)", a$prob_P_at_least_as_large_as_T, "fraction", paste0(f, ": prob_P_at_least_as_large_as_T")),
    kn(cl, "Residual change", a$residual_change_mm, "mm w.e.", paste0(f, ": residual_change_mm"))))) }
  if (!is.null(c2)) { v <- c2[region_id == bid]; if (nrow(v)) add(data.table::rbindlist(list(
    kn("2 interannual variance", "R-squared of SWE ~ P + T", v$R2_both, "", "Claim2_Table2_variance_shares.csv: R2_both"),
    kn("2 interannual variance", "Precipitation share of explained variance", v$share_of_R2_P, "%", "Claim2_Table2_variance_shares.csv: share_of_R2_P"),
    kn("2 interannual variance", "Precipitation share, 95% CI low", v$share_of_R2_P_lo, "%", "Claim2_Table2_variance_shares.csv: share_of_R2_P_lo"),
    kn("2 interannual variance", "Precipitation share, 95% CI high", v$share_of_R2_P_hi, "%", "Claim2_Table2_variance_shares.csv: share_of_R2_P_hi"),
    kn("2 interannual variance", "Temperature share of explained variance", v$share_of_R2_T, "%", "Claim2_Table2_variance_shares.csv: share_of_R2_T")))) }
  if (!is.null(c4)) { v <- c4[region_id == bid]; if (nrow(v)) add(data.table::rbindlist(list(
    kn("2 trend-based attribution", "SWE trend over the record", v$SWE_trend_total_mm, "mm w.e.", "Claim2_Table4_trend_based_attribution.csv: SWE_trend_total_mm"),
    kn("2 trend-based attribution", "Precipitation contribution (trend x span)", v$contribution_P_mm, "mm w.e.", "Claim2_Table4_trend_based_attribution.csv: contribution_P_mm"),
    kn("2 trend-based attribution", "Temperature contribution (trend x span)", v$contribution_T_mm, "mm w.e.", "Claim2_Table4_trend_based_attribution.csv: contribution_T_mm")))) }
  kt <- data.table::rbindlist(R, fill = TRUE)
  if (nrow(kt)) .ev_save_tbl(kt, "Claims_key_numbers")
  
  # manifest for evidence + link folders
  desc <- c(Claims_key_numbers = "Every number behind statements 1 and 2, with the file and column it comes from",
            Claim1_Table1_model_comparison = "AR(1) GLS comparison (dAIC, Akaike weights, likelihood-ratio p) of constant / trend / step / trend+step; all regions; temperature, precipitation, SWE (+ temperature without January)",
            Claim1_Table2_parametric_bootstrap = "Which truth (linear trend vs abrupt step) reproduces the observed fit and the Pettitt result; break-year behaviour under a pure trend",
            Claim1_Table3_monthly_contributions = "Basin: temperature change and trend by month, and each month's share of the Nov-Mar change",
            Claim1_Table4_trend_and_step_coefficients = "Pipeline output: trend and step coefficients when both are allowed (temperature needs no step; precipitation and SWE do)",
            Claim2_Table1_attribution_by_window = "SWE change split into precipitation, temperature and residual parts, with bootstrap CIs, for several comparison periods",
            Claim2_Table2_variance_shares = "Share of explained interannual SWE variance due to precipitation and to temperature (Shapley/LMG), with CIs",
            Claim2_Table3_decadal_contributions = "Decade means of the SWE anomaly and its precipitation, temperature and residual parts",
            Claim2_Table4_trend_based_attribution = "Attribution using linear trends of precipitation and temperature over the whole record",
            Claim1_Fig1_trend_vs_step = "Basin Nov-Mar temperature with linear trend and step fits", Claim1_Fig2_dAIC_by_region = "dAIC of each model by region and series",
            Claim1_Fig3_bootstrap_step_vs_trend = "Bootstrap: D under trend and step truth; Pettitt break-year distribution under a pure trend",
            Claim1_Fig4_monthly_contributions = "Temperature change and trend by calendar month",
            Claim2_Fig1_attribution_basin_and_regions = "Basin attribution bars and regional precipitation/temperature contributions with CIs",
            Claim2_Fig2_attribution_by_window = "Precipitation vs temperature contribution for different comparison periods",
            Claim2_Fig3_decadal_contributions = "Decade-by-decade composition of the SWE anomaly", Claim2_Fig4_standardised_series = "Standardised precipitation, temperature and SWE with 11-year running means",
            Claim2_Fig5_variance_shares = "Interannual variance shares by region",
            Link_Table0_inputs_used = "Series read for the snow-to-summer link, with record lengths",
            Link_Table1_correlations = "Snow predictors vs SPEI (PM and ERA5) in Jul-Sep: raw and detrended correlations, overlap flag, BH-adjusted p",
            Link_Table2_joint_exceedance = "Conditional frequency of low SPEI given low snow, with Fisher p", Link_Table3_by_era = "Correlations in 1950-78, 1979-2025 and 1991-2020",
            Link_Table4_demand_trends = "Trends in calendar-year PET, VPD, temperature, radiation and wind: share of pixels significant, by annual metric and Jul/Aug/Sep",
            Link_Table5_partial_correlations = "Snow-SPEI correlation after controlling for spring precipitation / water balance and linear time",
            Link_Fig1_correlation_heatmap = "Heatmap of detrended snow-SPEI correlations", Link_Fig2_scatter_snow_vs_SPEI3 = "Scatter of the snow predictor vs SPEI-3 for Jul, Aug, Sep",
            Link_Fig3_timeseries_snow_and_SPEI = "Snow index and August SPEI-3 through time")
  fl <- rbind(
    if (dir.exists(file.path(SWE_EVID_DIR, "tables"))) data.frame(folder = file.path("evidence", "tables"), file = list.files(file.path(SWE_EVID_DIR, "tables"), "\\.csv$")),
    if (dir.exists(file.path(SWE_EVID_DIR, "figures"))) data.frame(folder = file.path("evidence", "figures"), file = list.files(file.path(SWE_EVID_DIR, "figures"), "\\.png$")),
    if (dir.exists(file.path(SWE_LINK_DIR, "tables"))) data.frame(folder = file.path("link_to_summer", "tables"), file = list.files(file.path(SWE_LINK_DIR, "tables"), "\\.csv$")),
    if (dir.exists(file.path(SWE_LINK_DIR, "figures"))) data.frame(folder = file.path("link_to_summer", "figures"), file = list.files(file.path(SWE_LINK_DIR, "figures"), "\\.png$")))
  if (!is.null(fl) && nrow(fl)) { fl$description <- unname(desc[sub("\\.(csv|png)$", "", fl$file)])
  data.table::fwrite(fl, .resolve_out(file.path(SWE_HOMOG_DIR, "Evidence_manifest.csv")))
  cat(sprintf("  %d key numbers | manifest: %d files (%s)\n", nrow(kt), nrow(fl), file.path(SWE_HOMOG_DIR, "Evidence_manifest.csv"))) }
  invisible(kt)
}
####################################################################################
# REGIONAL PM-PET CONTROLS (9pr_pet_trends.r v9 adds PM_*_year / PM_*_AprOct to regional_avg_annual)
####################################################################################
create_regional_pm_control_plots <- function() {
  ra <- metadata$regional_avg_annual
  if (is.null(ra)) { cat("  No regional_avg_annual in metadata - skipping.\n"); return(invisible(NULL)) }
  ra  <- data.table::as.data.table(ra)
  ctl <- PM_CONTROL_VARS[paste0(PM_CONTROL_VARS, "_year") %in% names(ra)]
  if (!length(ctl)) { cat("  PM controls not in regional_avg_annual (needs 9pr_pet_trends.r v9+) - skipped.\n"); return(invisible(NULL)) }
  reg_dir <- file.path(out_dir, "regional_plots"); dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
  bid  <- .swe_basin_id()
  labs_ctl <- paste0(PM_CONTROL_LABELS[ctl], " (", sub(" yr-1$", "", PM_CONTROL_UNITS[ctl]), ")")
  for (suf in c("_year", "_AprOct")) {
    cols <- intersect(paste0(ctl, suf), names(ra)); if (!length(cols)) next
    l <- data.table::melt(ra[, c("region_id", "region_name", "year", cols), with = FALSE],
                          id.vars = c("region_id", "region_name", "year"), variable.name = "control", value.name = "value")
    l[, control := factor(sub(paste0(suf, "$"), "", as.character(control)), levels = ctl, labels = labs_ctl)]
    l[, basin := region_id == bid]
    p <- ggplot(l[is.finite(value)], aes(year, value, group = region_name, colour = basin, linewidth = basin)) +
      geom_line() +
      scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "black"),
                          labels = c(`FALSE` = "sub-basins", `TRUE` = "whole basin"), name = NULL) +
      scale_linewidth_manual(values = c(`FALSE` = 0.35, `TRUE` = 0.9), guide = "none") +
      facet_wrap(~control, scales = "free_y", ncol = 1) +
      labs(title = sprintf("PM-PET meteorological controls by region (%s)",
                           if (suf == "_year") "calendar-year mean" else "Apr-Oct mean"),
           x = "Year", y = NULL) +
      theme_bw(base_size = 9) + theme(legend.position = "top")
    .save_gg(file.path(reg_dir, paste0("regional_PM_controls", suf)), p, 10, 2.2 * length(ctl) + 1, formats = "png")
  }
  cat("  Regional PM-control figures written.\n")
  invisible(NULL)
}

####################################################################################
# FUNCTION 24: MULTI-VARIABLE HOMOGENEITY FIGURES (reads homogeneity/ from 9pr_pet_trends.r v9+)
# Classification per series:
#   "Step supported"      : AR(1) step + trend model favours the step (dAIC < -2 and step p < 0.05)
#   "Pettitt break only"  : Pettitt p < 0.05 but the AR(1) model does not support a step (trend-like)
#   "No break"            : otherwise
# Pettitt p is unadjusted here; the BH-adjusted values are in the source CSVs.
####################################################################################
HOMOG_DIR <- file.path(out_dir, "homogeneity")
.HOMOG_CLASS_COL <- stats::setNames(c("#b2182b", "#ef8a62", "grey60"),
                                    c(sprintf("Step at %d supported (AR(1) step + trend model)", .swe_step_year()),
                                      "Pettitt break only (trend-like)", "No break"))
.homog_class <- function(dAIC, step_p, pettitt_p) {
  data.table::fifelse(is.finite(dAIC) & dAIC < -2 & is.finite(step_p) & step_p < 0.05, names(.HOMOG_CLASS_COL)[1],
                      data.table::fifelse(is.finite(pettitt_p) & pettitt_p < 0.05, names(.HOMOG_CLASS_COL)[2],
                                          names(.HOMOG_CLASS_COL)[3]))
}

create_multivar_homogeneity_figures <- function() {
  cat("\n-- Function 24: multi-variable homogeneity figures --\n")
  rd <- function(f) { p <- file.path(HOMOG_DIR, f); if (file.exists(p)) data.table::fread(p) else NULL }
  brk <- rd("homogeneity_break_summary_basin.csv")
  if (is.null(brk) || !nrow(brk)) {
    cat("  homogeneity/ outputs not found - run 9pr_pet_trends.r (v9 or later) first.\n"); return(invisible(NULL))
  }
  reg <- rd("step_tests_by_region_all.csv"); pix <- rd("pixel_step_summary.csv")
  sy  <- .swe_step_year()
  hw  <- suppressWarnings(as.integer(metadata$homogeneity$break_window))
  if (!length(hw) || !is.finite(hw)) hw <- 3L
  yr  <- range(basin_avg_annual$year, na.rm = TRUE); mid <- mean(yr)
  fig_dir <- file.path(HOMOG_DIR, "figures"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- A) basin: break year, Pettitt strength and step-vs-trend class for every series ----
  b <- data.table::copy(brk)
  b[, class := factor(.homog_class(delta_AIC_step_vs_nostep, step_p, pettitt_p), levels = names(.HOMOG_CLASS_COL))]
  b[, kind := factor(kind, levels = c("raw", "ratio", "residual"))]
  data.table::setorder(b, family, kind, series)
  lv <- unique(b$series)
  b[, series := factor(series, levels = rev(lv))]
  b[, lab := sprintf("Pettitt p %s | dAIC %s", .fmt_p(pettitt_p),
                     ifelse(is.finite(delta_AIC_step_vs_nostep), sprintf("%+.1f", delta_AIC_step_vs_nostep), "NA"))]
  data.table::fwrite(b, .resolve_out(file.path(HOMOG_DIR, "homogeneity_break_classification_basin.csv")))
  pA <- ggplot(b, aes(pettitt_first_year_after, series)) +
    annotate("rect", xmin = sy - hw, xmax = sy + hw, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
    geom_vline(xintercept = sy, linetype = "dashed", colour = "grey35") +
    geom_vline(xintercept = mid, linetype = "dotted", colour = "grey35") +
    geom_point(aes(colour = class, shape = kind, size = -log10(pmax(pettitt_p, 1e-10))), na.rm = TRUE) +
    geom_text(aes(x = yr[2] + 2, label = lab), hjust = 0, size = 2.4, colour = "grey30") +
    facet_grid(family ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = .HOMOG_CLASS_COL, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c(raw = 16, ratio = 17, residual = 15), name = NULL, drop = FALSE) +
    scale_size_continuous(range = c(1.2, 4.5), name = "-log10 Pettitt p") +
    scale_x_continuous(limits = c(yr[1], yr[2] + 32), breaks = seq(1950, 2020, 10)) +
    labs(title = "Break year and step-vs-trend classification, basin mean (all ERA5-Land series)",
         subtitle = sprintf(paste0("Point = Pettitt first year of the later regime. Dashed = a-priori step year %d (grey band +/-%d yr); ",
                                   "dotted = record midpoint, where a smooth trend places a Pettitt break.\n",
                                   "Red = step supported by the AR(1) step + trend model (dAIC < -2, step p < 0.05), plotted at its Pettitt year; ",
                                   "orange = Pettitt break that the AR(1) model treats as a trend."), sy, hw),
         x = "First year of the later regime", y = NULL) +
    .swe_theme(9) + theme(strip.text.y = element_text(angle = 0), legend.box = "vertical")
  sA <- .save_gg(file.path(fig_dir, "Homog_Fig1_basin_break_classification"), pA, 12, 0.3 * nrow(b) + 2.5,
                 dpi = 200, formats = "png")
  
  # ---- B) classification by region -------------------------------------------------------
  sB <- NULL
  if (!is.null(reg) && nrow(reg)) {
    r <- data.table::copy(reg)
    r[, class := factor(.homog_class(delta_AIC_step_vs_nostep, step_p, pettitt_p), levels = names(.HOMOG_CLASS_COL))]
    bid <- .swe_basin_id()
    r[, region := factor(region_id, levels = c(bid, setdiff(sort(unique(region_id)), bid)))]
    r[, series := factor(series, levels = rev(intersect(lv, unique(series))))]
    r[, lab := ifelse(class == names(.HOMOG_CLASS_COL)[3], "", as.character(pettitt_first_year_after))]
    pB <- ggplot(r, aes(region, series, fill = class)) +
      geom_tile(colour = "white") + geom_text(aes(label = lab), size = 2.3) +
      scale_fill_manual(values = .HOMOG_CLASS_COL, name = NULL, drop = FALSE) +
      labs(title = "Step-vs-trend classification by region",
           subtitle = "Label = Pettitt first year of the later regime. Sub-basins are nested in the basin, so columns are not independent.",
           x = NULL, y = NULL) +
      .swe_theme(8) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    sB <- .save_gg(file.path(fig_dir, "Homog_Fig2_classification_by_region"), pB, 10, 0.25 * nlevels(r$series) + 2.5,
                   dpi = 200, formats = "png")
  }
  
  # ---- C) pixel break shares, D) pixel break-year maps -------------------------------------
  sC <- sD <- NULL
  if (!is.null(pix) && nrow(pix)) {
    p <- data.table::copy(pix)[order(pct_pettitt_fdr05)]
    p[, series := factor(series, levels = series)]
    p[, lab := ifelse(is.finite(median_break_year_fdr_sig), sprintf("median break %.0f", median_break_year_fdr_sig), "")]
    pC <- ggplot(p, aes(pct_pettitt_fdr05, series)) +
      geom_col(fill = "#4575b4", width = 0.7) +
      geom_text(aes(label = lab), hjust = -0.1, size = 2.6) +
      scale_x_continuous(limits = c(0, 130), breaks = seq(0, 100, 25)) +
      labs(title = "Share of basin pixels with a Pettitt break (BH-FDR < 0.05)",
           subtitle = "Near-universal breaks at the record midpoint usually reflect a trend rather than a step (see Fig 1).",
           x = "% of pixels", y = NULL) + .swe_theme(9)
    sC <- .save_gg(file.path(fig_dir, "Homog_Fig3_pixel_break_share"), pC, 9, 0.28 * nrow(p) + 2,
                   dpi = 200, formats = "png")
    
    sel <- p[pct_pettitt_fdr05 >= 5, as.character(series)]
    pm <- data.table::rbindlist(lapply(sel, function(s) {
      f <- file.path(HOMOG_DIR, "pixel", sprintf("pixel_step_tests_%s.csv", s))
      if (!file.exists(f)) return(NULL)
      x <- data.table::fread(f); x[, series := s]; x
    }), fill = TRUE)
    if (nrow(pm)) {
      pm[, brk := ifelse(is.finite(pettitt_p_fdr) & pettitt_p_fdr < 0.05, pettitt_first_year_after, NA_real_)]
      pm[, series := factor(series, levels = rev(sel))]
      dx <- .grid_spacing(pm$x); dy <- .grid_spacing(pm$y); nc <- min(4L, length(sel))
      pD <- ggplot(pm, aes(x, y, fill = brk)) +
        geom_tile(width = dx, height = dy) +
        geom_sf(data = basin_wgs, inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.4) +
        facet_wrap(~series, ncol = nc) + coord_sf(expand = FALSE) +
        scale_fill_viridis_c(option = "plasma", na.value = "grey90", name = "Break year\n(FDR-significant)") +
        labs(title = "Pixel break years (Pettitt, BH-FDR < 0.05; grey = no significant break)", x = NULL, y = NULL) +
        theme_void(base_size = 9) +
        theme(plot.title = element_text(face = "bold", size = 11), strip.text = element_text(face = "bold", size = 8),
              legend.position = "bottom")
      sD <- .save_gg(file.path(fig_dir, "Homog_Fig4_pixel_break_year_maps"), pD, 3.2 * nc, 3 * ceiling(length(sel) / nc) + 1.2,
                     dpi = 200, formats = "png")
    }
  }
  
  tab <- b[, .N, by = class]
  cat(sprintf("  Basin series: %s\n", paste(sprintf("%s = %d", tab$class, tab$N), collapse = " | ")))
  cat(sprintf("  Figures: %s\n", paste(c(sA, sB, sC, sD), collapse = ", ")))
  invisible(b)
}
####################################################################################
# GENERATE ALL OUTPUTS  (each step isolated: one failure never stops the rest)
####################################################################################
cat("\n========================================\nCREATING VISUALIZATIONS\n========================================\n\n")

for (var in c("Precipitation", "PET", "PET_ERA5", "Temperature", SWE_VARS_PRESENT, EXTRA_VARS_PRESENT))
  for (meth in c("vc", "tfpw"))
    .run_step(sprintf("publication map %s/%s", var, meth), create_publication_figure(var, meth))

.run_step("Function 0  (hydroclimatic trend map)",         create_supp_hydroclimatic_trend_figure())
.run_step("Function 0b (Pr + PM-control trend map)",       create_pm_control_trend_figure())
.run_step("Function 0c (SWE trend map)",                   create_swe_trend_figure())
.run_step("Function 0d (SWE monthly trend maps)",          create_swe_monthly_trend_maps())
.run_step("Function 0e (AET / runoff / snow trend map)",   create_extra_trend_figure())
for (v in EXTRA_VARS_PRESENT)
  .run_step(sprintf("Function 0f (%s monthly trend maps)", v), create_var_monthly_trend_maps(v))
.run_step("Function 2  (tau comparison map)",              create_comparison_figure())
.run_step("Function 3  (regime shifts)",                   create_regime_shift_maps())
.run_step("Function 4  (spectral)",                        create_spectral_analysis_maps())
.run_step("Function 5  (spatial patterns)",                create_spatial_pattern_analysis())
.run_step("Function 6  (method comparison)",               create_method_comparison())
.run_step("Function 7  (statistics table)",                create_statistics_table())
.run_step("Function 8  (basin time series)",               create_basin_timeseries_plot())
.run_step("Function 13 (calendar-month plots)",            create_all_calmonth_temporal_changes())
.run_step("Function 14 (PET PM vs ERA5)",                  create_pet_comparison_plots())
.run_step("Function 15 (PET seasonal trends)",             create_pet_seasonal_trends())
.run_step("Function 16 (precip seasonal trends)",          create_precip_seasonal_trends())
.run_step("Function 17 (temperature seasonal trends)",     create_temp_seasonal_trends())
.run_step("Function 18 (SWE seasonal trends)",             create_swe_seasonal_trends())
.run_step("Function 18b (AET / runoff / snow seasonal)",   create_extra_seasonal_trends())
.run_step("Seasonality P vs PET",                          create_seasonality_p_pet_plot())
.run_step("Function 12 (temperature dedicated)",           create_temperature_dedicated_plots())
.run_step("Functions 9  (point time series)",              create_point_timeseries_plots())
.run_step("Functions 10 (point monthly slopes)",           create_point_monthly_trend_plots())
.run_step("Function 11 (point summary table)",             create_point_summary_table())
.run_step("Regional climate plots",                        create_regional_climate_plots())
.run_step("Regional PM-control plots",                     create_regional_pm_control_plots())
.run_step("Function 19 (SWE step-change outputs)",         create_swe_step_change_outputs())
.run_step("Function 21 (claim 1: trend vs step)",          create_claim1_trend_vs_step_outputs())
.run_step("Function 22 (claim 2: P vs T attribution)",     create_claim2_attribution_outputs())
.run_step("Function 23 (snow -> summer SPEI link)",        create_snow_summer_link_outputs())
.run_step("Evidence key numbers and manifest",             create_claims_key_numbers())
.run_step("Function 24 (multi-variable homogeneity)",      create_multivar_homogeneity_figures())

cat("\n========================================\nALL VISUALIZATION STEPS FINISHED\n========================================\n")
cat(sprintf("Output directory: %s\n", normalizePath(out_dir)))
cat("Check the console above for any [FAILED] or [WARN] lines.\n")