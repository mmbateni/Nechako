##############################################################################
# 8_SSI_SSMI_ERALand.R
# SOIL MOISTURE DROUGHT INDICES - PARALLEL VERSION
#
#   SSI-SM : Standardized Soil moisture Index (parametric)
#            Beta distribution (method of moments) -> empirical fallback.
#            After Hao & AghaKouchak (2013); standardisation as McKee et al. (1993).
#   SSMI   : Standardized Soil Moisture Index (non-parametric)
#            Logit boundary correction -> Gaussian KDE-CDF -> N(0,1).
#            After Carrao et al. (2013) / Russo et al.
#   Both   : distribution fitted separately for every calendar month
#            (Kao & Govindaraju 2010), whole-record reference.
#
# Input : ERA5-Land monthly volumetric soil water, layers 1-4
#         (monthly_data_direct/volumetric_soil_water_layer_{1..4}_monthly.nc)
# Shared-path version aligned with 1SPI_ERALand.R and 5SRI_SDI_ERALand.R
#
# NAMING (framework index keys):
#   "ssi_sm" is used for the soil-moisture SSI throughout (files, directories,
#   index_type) because "ssi" is already the Standardized STREAMFLOW Index of
#   6_StreamFlow_Index_WSCLocal.R and 12basin_timeseries.R. "ssmi" is unique.
#
# LAYER CONFIGURATIONS (ERA5-Land, depth-weighted mean volumetric water):
#   L1_3 : layers 1-3,  0-100 cm  (7 + 21 + 72 cm)          <- PRIMARY_LAYER
#   L1_4 : layers 1-4,  0-289 cm  (7 + 21 + 72 + 189 cm)
#
# TIMESCALES: 1, 2, 3 months (right-aligned rolling MEAN; soil moisture is a
#   state variable, so it is averaged rather than summed).
#
# FRAMEWORK CONTRACT (same as Scripts 1, 3, 4, 5):
#   - working directory = project root; paths derived from getwd();
#     shared helpers from DROUGHT_ANALYSIS_utils.R.
#   - native ERA5-Land precipitation grid; soil-moisture grids must match it
#     exactly. NOTHING is reprojected, resampled, rotated or cropped.
#   - basin = dissolved union of the 7 authoritative subbasins; canonical
#     basin mask rebuilt with the SPI source label so SPI, SRI, SSI-SM and
#     SSMI share identical basin pixels.
#   - basin / regional means: cell area x fractional coverage weights, reported
#     only when >= BASIN_MIN_VALID_AREA_COVERAGE (0.20) of the area is valid.
#   - drought events in the basin plots: identify_drought_events_ts()
#     (DROUGHT_ONSET / DROUGHT_END / EVENT_MIN_DURATION_BY_SCALE).
#   - empirical fallback: Gringorten plotting position, as SPI/SRI/SWEI/SSPI.
#   - worker failures stop the run (no silently NA-filled pixels).
#
# OUTPUTS (per index x layer x timescale; <key> = ssi_sm | ssmi):
#   <key>_results_seasonal/<layer>/
#     <key>_SS_monthMM_Mon.nc / .csv          12 calendar-month files (Script 11 layout)
#     <key>_SS_all_months.xlsx                 all calendar months as sheets
#     <key>_SS_distribution_map.png            dominant fitting method per pixel
#     <key>_SS_basin_averaged_by_month.csv     Year x Jan..Dec, area-weighted
#     <key>_SS_basin_valid_area_coverage_by_month.csv
#     <key>_SS_basin_timeseries.png            basin series with framework events
#     <key>_basin_timeseries.csv               basin SM + all scales, long format
#     <key>_<layer>_summary.txt
#   NECHAKO_REGION_OUTPUT_DIR/
#     <key>_<layer>_SS_regional_averaged_by_month.csv  basin + every subbasin
#   ssmi_results_seasonal/soil_moisture_index_output_manifest.csv
#
# [Nechako SAF pipeline] A03_saf_core_pipeline.R reads the legacy flat names
#   ssmi_results_seasonal/ssmi_L1_3_01_month{MM}_{Mon}.nc
# via assemble_index_raster(prefix = "ssmi_L1_3", scale_tag = "01"). These are
# still written (WRITE_LEGACY_SAF_FILES = TRUE) for both layer configs. NOTE:
# they are now on the native ERA5-Land grid (framework rule), no longer on the
# BC Albers (EPSG:3005) grid the previous version of this script produced.
#
# EXECUTION ORDER: 1, 2b, 3, 4, 5 (6, 7), 8 (this script) -> 11 -> 12 -> 13 -> 14
##############################################################################

# ---- Libraries ----
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

# ---- Paths / shared utilities ----
source(UTILS_FILE)
utils_load_packages(c("terra", "ncdf4", "writexl", "data.table", "ggplot2"))
library(parallel)

data_dir     <- file.path(WD_PATH, "monthly_data_direct")
ssi_sm_root  <- file.path(WD_PATH, "ssi_sm_results_seasonal")
ssmi_root    <- file.path(WD_PATH, "ssmi_results_seasonal")
REGIONAL_DIR <- if (exists("NECHAKO_REGION_OUTPUT_DIR")) NECHAKO_REGION_OUTPUT_DIR else
  file.path(WD_PATH, "temporal_drought", "regional_timeseries")
for (d in c(ssi_sm_root, ssmi_root, REGIONAL_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# USER SETTINGS
# ==============================================================================
sm_inputs <- c(
  L1 = "volumetric_soil_water_layer_1_monthly.nc",
  L2 = "volumetric_soil_water_layer_2_monthly.nc",
  L3 = "volumetric_soil_water_layer_3_monthly.nc",
  L4 = "volumetric_soil_water_layer_4_monthly.nc"
)
NATIVE_PRECIP_FILE <- file.path(data_dir, "total_precipitation_monthly.nc")

# ERA5-Land soil layer thicknesses (cm): 0-7, 7-28, 28-100, 100-289
layer_depths_cm <- c(L1 = 7, L2 = 21, L3 = 72, L4 = 189)
layer_configs <- list(
  L1_3 = list(layers = c("L1", "L2", "L3"),       desc = "Layers 1-3 (0-100 cm, depth-weighted)"),
  L1_4 = list(layers = c("L1", "L2", "L3", "L4"), desc = "All 4 layers (0-289 cm, depth-weighted)")
)
PRIMARY_LAYER <- "L1_3"              # layer used by Script 11 and the SAF pipeline

sm_timescales     <- c(1L, 2L, 3L)
min_sample_size   <- 5L              # min. years per calendar month (same as SPI/SRI)
z_clip            <- 4.75            # same clamp as SPI/SRI
eps_clip          <- 1e-6
sm_valid_max      <- 1 + 1e-6        # volumetric fraction must lie in [0, 1]
WRITE_LEGACY_SAF_FILES <- TRUE       # flat ssmi_<layer>_SS_monthMM_Mon.nc for A03
use_parallel      <- TRUE

stopifnot(PRIMARY_LAYER %in% names(layer_configs),
          all(sm_timescales >= 1L),
          all(unlist(lapply(layer_configs, `[[`, "layers")) %in% names(sm_inputs)))

# Note: the SPI/SRI low-variance rule (var < 0.01 mm^2 -> empirical) is NOT
# used here. Volumetric soil water is a fraction (var ~1e-4..1e-3), so that
# threshold would force every month onto the empirical branch.
cfg <- list(min_n = min_sample_size, z_clip = z_clip, eps = eps_clip,
            timescales = sm_timescales)

INDEX_DEFS <- list(
  ssi_sm = list(key = "ssi_sm", label = "SSI-SM", root = ssi_sm_root,
                longname = "Standardized Soil moisture Index (SSI, Beta)",
                distribution = "Beta (method of moments) | empirical Gringorten fallback",
                reference = "Hao & AghaKouchak (2013); McKee et al. (1993); Kao & Govindaraju (2010)",
                method_labels = c("Beta (MoM)", "Empirical", "Zero-Variance")),
  ssmi   = list(key = "ssmi", label = "SSMI", root = ssmi_root,
                longname = "Standardized Soil Moisture Index (SSMI, logit-KDE)",
                distribution = "Gaussian KDE in logit space (Silverman bw) | empirical Gringorten fallback",
                reference = "Carrao et al. (2013); Russo et al.; Kao & Govindaraju (2010)",
                method_labels = c("KDE (logit)", "Empirical", "Zero-Variance"))
)

# ==============================================================================
# I/O HELPERS (same contract as 5SRI_SDI_ERALand.R)
# ==============================================================================
nc_time_to_date <- function(vals, units) {
  m <- regmatches(units, regexec("^\\s*(seconds|minutes|hours|days)\\s+since\\s+(.+)$", units))[[1]]
  if (length(m) != 3L) stop("Unrecognised NetCDF time units: ", units)
  mult <- c(seconds = 1, minutes = 60, hours = 3600, days = 86400)[[m[2]]]
  origin <- as.POSIXct(substr(trimws(m[3]), 1, 19), tz = "UTC",
                       tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d"))
  as.Date(origin + vals * mult, tz = "UTC")
}

read_era5_monthly <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  data_vars <- names(nc$var)[vapply(nc$var, function(v) v$ndims == 3L, logical(1))]
  if (length(data_vars) != 1L)
    stop("Expected exactly one 3-D data variable in ", path, "; found: ",
         paste(data_vars, collapse = ", "))
  vname <- data_vars[[1L]]
  
  get_att <- function(var, att) {
    a <- ncdf4::ncatt_get(nc, var, att)
    if (isTRUE(a$hasatt)) trimws(as.character(a$value)) else NA_character_
  }
  
  tdim <- intersect(c("valid_time", "time"), names(nc$dim))
  if (!length(tdim)) stop("No time dimension found in ", path)
  tdim <- tdim[[1L]]
  dates_nc <- nc_time_to_date(nc$dim[[tdim]]$vals, nc$dim[[tdim]]$units)
  
  r <- tryCatch(terra::rast(path, subds = vname), error = function(e) terra::rast(path))
  if (terra::nlyr(r) != length(dates_nc)) stop("Layer/time mismatch in ", path)
  
  list(rast = r, var = vname, units = get_att(vname, "units"),
       long_name = get_att(vname, "long_name"), dates = dates_nc)
}

normalise_monthly_dates <- function(d) {
  d <- as.Date(format(d, "%Y-%m-01"))
  mi <- as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))
  if (anyDuplicated(mi) || any(diff(mi) != 1L))
    stop("Time axis is not a continuous monthly sequence.")
  d
}

# Volumetric soil water: explicit unit strings only; magnitude is checked separately.
# Returns the multiplier to a volumetric fraction (m3/m3), or NA if unsupported.
sm_unit_factor <- function(u) {
  if (is.na(u)) return(NA_real_)
  u0 <- gsub("\\s+", " ", tolower(trimws(u)))
  u0 <- gsub("\\*\\*|\\^", "", u0)               # m**3 m**-3 / m^3 m^-3 -> m3 m-3
  if (u0 %in% c("m3 m-3", "m3/m3", "m3.m-3", "m3m-3", "1", "fraction",
                "volumetric fraction", "dimensionless")) return(1)
  if (u0 %in% c("%", "percent", "vol%", "% vol")) return(0.01)
  NA_real_
}

# Native-grid contract (Script 3/4 rule; terra >= 1.8: origin checked separately).
assert_native_grid <- function(x, template, label) {
  geom_ok <- terra::compareGeom(x[[1]], template, stopOnError = FALSE,
                                crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE)
  orig_ok <- isTRUE(all.equal(terra::origin(x[[1]]), terra::origin(template), tolerance = 1e-6))
  if (!isTRUE(geom_ok) || !orig_ok)
    stop(label, " does not match the native ERA5-Land precipitation grid.\n",
         "Script 8 does not reproject or resample raster inputs. Check the upstream NetCDF producer.")
  invisible(TRUE)
}

# ==============================================================================
# STATISTICAL CORE (exported to workers)
# ==============================================================================
# Right-aligned k-month rolling mean; NA propagates; first k-1 values are NA.
roll_mean_right <- function(x, k) {
  if (k == 1L) return(x)
  if (length(x) < k) return(rep(NA_real_, length(x)))
  as.numeric(stats::filter(x, rep(1 / k, k), method = "convolution", sides = 1))
}

clip_prob  <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)
gringorten <- function(x) (rank(x, ties.method = "average") - 0.44) / (length(x) + 0.12)

# Beta CDF by method of moments on (0, 1); NULL if the moments are not admissible.
beta_mom_prob <- function(samp, eps) {
  sb <- pmax(pmin(samp, 1 - eps), eps)
  mu <- mean(sb); v <- stats::var(sb)
  if (!is.finite(v) || v <= 0 || v >= mu * (1 - mu)) return(NULL)
  phi <- mu * (1 - mu) / v - 1
  a <- mu * phi; b <- (1 - mu) * phi
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) return(NULL)
  p <- stats::pbeta(sb, shape1 = a, shape2 = b)
  if (!all(is.finite(p))) return(NULL)
  p
}

# Kernel-smoothed CDF in logit space, evaluated at the sample points:
#   F(y_i) = (1/n) sum_j Phi((y_i - y_j) / h),  h = Silverman (bw.nrd0).
kde_logit_prob <- function(samp, eps) {
  sc <- pmax(pmin(samp, 1 - eps), eps)
  y  <- log(sc / (1 - sc))
  bw <- stats::bw.nrd0(y)
  if (!is.finite(bw) || bw <= 0) bw <- 0.5 * stats::IQR(y)
  if (!is.finite(bw) || bw <= 0) return(NULL)
  p <- rowMeans(stats::pnorm(outer(y, y, "-") / bw))
  if (!all(is.finite(p))) return(NULL)
  p
}

# Calendar-month standardisation shared by both indices.
# Method codes: 1 = parametric/KDE, 2 = empirical (Gringorten), 3 = zero variance
standardize_sm_monthly <- function(x, month_vec, cfg, prob_fun) {
  n <- length(x); z <- rep(NA_real_, n); meth <- rep(NA_integer_, n)
  for (m in 1:12) {
    idx <- which(month_vec == m & is.finite(x))
    if (length(idx) < cfg$min_n) next
    samp <- x[idx]; sv <- stats::var(samp)
    if (!is.finite(sv) || sv < .Machine$double.eps) { z[idx] <- 0; meth[idx] <- 3L; next }
    p <- tryCatch(prob_fun(samp, cfg$eps), error = function(e) NULL)
    if (!is.null(p)) {
      z[idx] <- stats::qnorm(clip_prob(p, cfg$eps)); meth[idx] <- 1L; next
    }
    z[idx] <- stats::qnorm(clip_prob(gringorten(samp), cfg$eps)); meth[idx] <- 2L
  }
  f <- is.finite(z)
  z[f] <- pmin(pmax(z[f], -cfg$z_clip), cfg$z_clip)
  list(z = z, method = meth)
}

process_pixel <- function(x, month_vec, cfg) {
  lapply(cfg$timescales, function(sc) {
    xa <- roll_mean_right(x, sc)
    list(ssi_sm = standardize_sm_monthly(xa, month_vec, cfg, beta_mom_prob),
         ssmi   = standardize_sm_monthly(xa, month_vec, cfg, kde_logit_prob))
  })
}

# Worker reads data exported once to each worker's global environment.
worker_fun <- function(i) {
  tryCatch(process_pixel(.wk_mat[i, ], .wk_month, .wk_cfg),
           error = function(e) list(error = conditionMessage(e)))
}

# ==============================================================================
# OUTPUT HELPERS
# ==============================================================================
write_index_nc <- function(sub_mat, layer_dates, template, valid_cells, file, varname, longname) {
  r <- terra::rast(template, nlyrs = ncol(sub_mat))
  full <- matrix(NA_real_, terra::ncell(template), ncol(sub_mat))
  full[valid_cells, ] <- sub_mat
  terra::values(r) <- full
  terra::time(r) <- layer_dates
  names(r) <- format(layer_dates, "%Y-%m")
  terra::writeCDF(r, file, varname = varname, longname = longname,
                  unit = "standardized_index", missval = -9999, overwrite = TRUE)
  if (!file.exists(file) || is.na(file.info(file)$size) || file.info(file)$size <= 0)
    stop("NetCDF write failed or produced an empty file: ", file)
  invisible(file)
}

plot_method_map <- function(meth, template, valid_cells, basin, file, title, labs, empirical_pct) {
  dom <- apply(meth, 1, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else as.numeric(names(which.max(table(x))))
  })
  full <- rep(NA_real_, terra::ncell(template)); full[valid_cells] <- dom
  mr <- terra::rast(template); terra::values(mr) <- full
  cols <- c("#4575b4", "#d73027", "#91bfdb")
  png(file, width = 1200, height = 800, res = 150); on.exit(dev.off(), add = TRUE)
  terra::plot(mr, col = cols, breaks = c(0.5, 1.5, 2.5, 3.5), legend = FALSE,
              main = title, axes = FALSE, box = FALSE)
  terra::plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  cnt <- table(factor(dom, levels = 1:3, labels = labs))
  legend("bottomright", legend = sprintf("%s (%d pixels)", names(cnt), as.integer(cnt)),
         fill = cols, title = "Method", cex = 0.85, bg = "white", bty = "o", border = "gray50")
  mtext(sprintf("Valid pixels: %d | Empirical: %.1f%%", length(valid_cells), empirical_pct),
        side = 1, line = 0.5, cex = 0.75, col = "gray30")
}

# Year x Jan..Dec wide table from a monthly vector aligned with `dates`
# (identical layout to Script 4 / Script 11 *_basin_averaged_by_month.csv).
by_month_wide <- function(v, dates) {
  yrs <- as.integer(format(dates, "%Y")); mos <- as.integer(format(dates, "%m"))
  years <- sort(unique(yrs))
  wide <- data.frame(Year = years)
  for (m in 1:12) {
    vv <- rep(NA_real_, length(years)); hit <- which(mos == m)
    vv[match(yrs[hit], years)] <- v[hit]
    wide[[month.abb[m]]] <- vv
  }
  wide
}

pct <- function(m, code) 100 * sum(m == code, na.rm = TRUE) / max(1, sum(!is.na(m)))

# Basin time series with framework drought events (as in 12basin_timeseries.R).
plot_basin_index <- function(dates, values, label, sc, layer, file) {
  ev <- tryCatch(identify_drought_events_ts(dates, values, scale = sc), error = function(e) NULL)
  in_event <- if (!is.null(ev)) attr(ev, "in_event") else NULL
  if (is.null(in_event) || length(in_event) != length(values)) in_event <- rep(FALSE, length(values))
  in_event[is.na(in_event)] <- FALSE
  n_ev    <- if (is.null(ev)) NA_integer_ else nrow(ev)
  min_dur <- tryCatch(get_event_min_duration(sc), error = function(e) NA_integer_)
  
  df <- data.frame(date = dates, value = values, in_event = in_event)
  df <- df[is.finite(df$value), , drop = FALSE]
  if (!nrow(df)) { cat(sprintf("    ! %s-%02d (%s): no finite basin values, plot skipped\n", label, sc, layer)); return(invisible(NULL)) }
  y_lo <- min(-3.5, floor(min(df$value) * 1.10)); y_hi <- max(3.5, ceiling(max(df$value) * 1.10))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(date, value)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -1, fill = "#C00000", alpha = 0.06) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = Inf, fill = "#4472C4", alpha = 0.06) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "gray50", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = DROUGHT_ONSET, linetype = "dotted", colour = "gray30", linewidth = 0.35) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ifelse(in_event, value, 0), ymax = 0),
                         fill = "#c0392b", alpha = 0.35) +
    ggplot2::geom_line(colour = "#1F77B4", linewidth = 0.5) +
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s-%02d  Basin-Averaged Time Series (Area-Weighted)  |  %s", label, sc, layer),
      subtitle = sprintf("%s drought events  (onset < %.2f, end >= %.2f, min %s months)",
                         ifelse(is.na(n_ev), "NA", as.character(n_ev)), DROUGHT_ONSET, DROUGHT_END,
                         ifelse(is.na(min_dur), "NA", as.character(min_dur))),
      x = "Year", y = sprintf("%s-%02d Index Value", label, sc)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ggplot2::ggsave(file, p, width = 11, height = 5.5, dpi = 200, bg = "white")
  invisible(n_ev)
}

# ==============================================================================
# CANONICAL GRID / BASIN MASK / REGIONS (identical to the SPI mask)
# ==============================================================================
cat("\n============================================================\n")
cat("SSI-SM / SSMI CALCULATION (PARALLEL MODE)\n")
cat("============================================================\n")

# Whole basin = dissolved union of the 7 authoritative subbasins
# (identical construction to 3SPEI_ERALand.R).
basin <- nechako_subbasin_union(load_nechako_subbasins())
cat("Basin boundary loaded (union of 7 subbasins)\n")

# Rebuild the canonical mask from the same precipitation grid and source label
# used by 1SPI_ERALand.R so that SPI, SRI, SSI-SM and SSMI share the same basin
# pixels. The soil-moisture grids must match this geometry exactly.
if (!file.exists(NATIVE_PRECIP_FILE))
  stop("Native ERA5-Land precipitation grid not found: ", NATIVE_PRECIP_FILE,
       "\nRun the monthly ERA5-Land downloader (and 1SPI_ERALand.R) first.")
grid_ref <- terra::rast(NATIVE_PRECIP_FILE)[[1]]
if (is.na(terra::crs(grid_ref)) || !nzchar(terra::crs(grid_ref)))
  stop("Native ERA5-Land precipitation grid has no CRS.")
grid_label <- "1SPI_ERALand.R canonical native ERA5-Land precipitation grid"
basin <- terra::project(basin, terra::crs(grid_ref))
common_grid <- save_common_analysis_grid(grid_ref, basin = basin, source = grid_label)
common_mask <- terra::unwrap(common_grid$basin_mask)
mask_vals   <- terra::values(common_mask, mat = FALSE)
basin_cells <- which(!is.na(mask_vals) & mask_vals != 0)
if (!length(basin_cells)) stop("Canonical basin mask contains zero cells.")

# Fractional-coverage basin weights (shared framework rule): w = cell area x
# fractional coverage; basin means need >= BASIN_MIN_VALID_AREA_COVERAGE (0.20)
# of the total basin weight to be reported.
sm_bw <- basin_fractional_weights(grid_ref, basin)
if (!identical(which(sm_bw$mask), basin_cells))
  stop("Fractional-weight basin mask differs from the canonical SSI-SM/SSMI basin mask.")
cat(sprintf("Basin weights: %d cells (%d partial) | weighted area %.1f km2 | min valid coverage %.2f\n",
            sm_bw$n_cells, sm_bw$n_partial, sm_bw$area_m2 / 1e6, BASIN_MIN_VALID_AREA_COVERAGE))
cat(sprintf("Canonical basin cells: %d\n", length(basin_cells)))

# Basin + every subbasin (shared definitions, as in Script 4 / Script 11).
regions <- load_nechako_regions(template = grid_ref)
validate_nechako_regions(regions)
region_weights <- region_fractional_weights(grid_ref, regions)
for (id in names(region_weights))
  cat(sprintf("  Region %-28s : %4d cells (%3d partial) | %.1f km2\n", id,
              as.integer(region_weights[[id]]$n_cells), as.integer(region_weights[[id]]$n_partial),
              region_weights[[id]]$area_m2 / 1e6))

# ==============================================================================
# LOAD SOIL-MOISTURE LAYERS (unit, grid and time contract)
# ==============================================================================
cat("\n===== SOIL MOISTURE UNIT / GRID / TIME CONTRACT CHECK =====\n")
layer_vals  <- list()
layer_meta  <- list()
dates       <- NULL
for (ly in names(sm_inputs)) {
  f    <- file.path(data_dir, sm_inputs[[ly]])
  prod <- read_era5_monthly(f)
  d    <- normalise_monthly_dates(prod$dates)
  if (is.null(dates)) dates <- d
  if (!identical(d, dates)) stop("Time axis of ", sm_inputs[[ly]], " differs from layer L1.")
  assert_native_grid(prod$rast, grid_ref, sprintf("Soil moisture %s (%s)", ly, sm_inputs[[ly]]))
  
  v <- terra::values(prod$rast, mat = TRUE)
  v[!is.finite(v) | v < -1e20 | v > 1e20] <- NA_real_   # fill-value guard
  v <- v[basin_cells, , drop = FALSE]
  
  fac <- sm_unit_factor(prod$units)
  if (is.na(fac)) {
    vmax <- suppressWarnings(max(v, na.rm = TRUE))
    if (is.finite(vmax) && vmax <= sm_valid_max) {
      fac <- 1
      warning(sprintf("[%s] units '%s' not recognised; magnitude (max %.3f) is consistent with m3/m3.",
                      ly, prod$units, vmax))
    } else stop("Unsupported or missing soil-moisture units in ", sm_inputs[[ly]], ": ", prod$units)
  }
  v <- v * fac
  
  n_neg <- sum(is.finite(v) & v < 0); n_big <- sum(is.finite(v) & v > sm_valid_max)
  if (n_big) stop(sprintf("[%s] %d values > 1 m3/m3 after unit handling; check %s.", ly, n_big, sm_inputs[[ly]]))
  if (n_neg) { v[is.finite(v) & v < 0] <- NA_real_
  warning(sprintf("[%s] %d negative soil-moisture values set to NA.", ly, n_neg)) }
  
  cat(sprintf("  %s : var '%s' | units '%s' (x%g) | range %.3f-%.3f m3/m3 | negatives -> NA: %d\n",
              ly, prod$var, ifelse(is.na(prod$units), "<missing>", prod$units), fac,
              suppressWarnings(min(v, na.rm = TRUE)), suppressWarnings(max(v, na.rm = TRUE)), n_neg))
  layer_vals[[ly]] <- v
  layer_meta[[ly]] <- list(file = sm_inputs[[ly]], var = prod$var, units = prod$units, factor = fac)
  rm(prod, v); invisible(gc(verbose = FALSE))
}
month_vec <- as.integer(format(dates, "%m"))
if (exists("validate_monthly_record", mode = "function"))
  validate_monthly_record(dates, label = "Soil-moisture monthly record", require_complete = TRUE)
if (month_vec[1] != 1L || month_vec[length(month_vec)] != 12L)
  warning("Record does not run January-December; Script 11 expects complete calendar years.")
cat(sprintf("Time period : %s to %s (%d months)\n", min(dates), max(dates), length(dates)))
cat(sprintf("Grid        : %d x %d native cells | %d canonical basin cells\n",
            terra::ncol(grid_ref), terra::nrow(grid_ref), length(basin_cells)))
cat("===== SOIL MOISTURE CONTRACT CHECK PASSED =====\n")

# ==============================================================================
# PARALLEL CLUSTER (one cluster reused for all layer configurations)
# ==============================================================================
cl <- NULL
if (use_parallel) {
  n_cores <- max(1L, detectCores() - 1L)
  if (n_cores > 1L) {
    cat(sprintf("\nStarting parallel cluster with %d cores\n", n_cores))
    dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)  # Windows temp-dir fix
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, { library(stats); NULL })
    clusterExport(cl, c("roll_mean_right", "clip_prob", "gringorten", "beta_mom_prob",
                        "kde_logit_prob", "standardize_sm_monthly", "process_pixel", "worker_fun"),
                  envir = globalenv())
  }
}

# ==============================================================================
# MAIN PROCESSING FUNCTION (one layer configuration)
# ==============================================================================
manifest_rows    <- list()
required_outputs <- character(0)

process_layer_config <- function(layer, cl) {
  lc <- layer_configs[[layer]]
  cat("\n============================================================\n")
  cat(sprintf("LAYER CONFIG: %s - %s\n", layer, lc$desc))
  cat("============================================================\n")
  
  # ---- Depth-weighted composite (NA if any layer is missing) ----
  w <- layer_depths_cm[lc$layers] / sum(layer_depths_cm[lc$layers])
  sm <- Reduce(`+`, lapply(lc$layers, function(ly) layer_vals[[ly]] * w[[ly]]))
  cat(sprintf("  Depth weights: %s\n", paste(sprintf("%s=%.4f", names(w), w), collapse = "  ")))
  
  keep        <- rowSums(is.finite(sm)) > 0
  valid_cells <- basin_cells[keep]
  sm_mat      <- sm[keep, , drop = FALSE]
  n_pix <- nrow(sm_mat)
  if (!n_pix) stop("No valid basin pixels for ", layer)
  cat(sprintf("  Valid basin pixels: %d of %d canonical cells (%d all-NA dropped)\n",
              n_pix, length(basin_cells), sum(!keep)))
  
  basin_sm     <- basin_weighted_mean(sm_mat, sm_bw, cells = valid_cells)
  basin_sm_cov <- basin_valid_coverage(sm_mat, sm_bw, cells = valid_cells)
  clim <- tapply(basin_sm, month_vec, mean, na.rm = TRUE)
  
  # ---- Pixel-wise SSI-SM + SSMI ----
  assign(".wk_mat", sm_mat, envir = globalenv())
  assign(".wk_month", month_vec, envir = globalenv())
  assign(".wk_cfg", cfg, envir = globalenv())
  cat(sprintf("\n  Processing %d pixels (SSI-SM + SSMI, scales %s)... ",
              n_pix, paste(sm_timescales, collapse = "/")))
  t0 <- Sys.time()
  if (!is.null(cl)) {
    clusterExport(cl, c(".wk_mat", ".wk_month", ".wk_cfg"), envir = globalenv())
    res <- parLapply(cl, seq_len(n_pix), worker_fun)
  } else {
    res <- lapply(seq_len(n_pix), worker_fun)
  }
  cat(sprintf("done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  errs <- vapply(res, function(p) !is.null(p$error), logical(1))
  if (any(errs))
    stop(sprintf("[%s] %d/%d pixels failed. First worker error: %s",
                 layer, sum(errs), n_pix, res[[which(errs)[1L]]]$error))
  
  out <- lapply(names(INDEX_DEFS), function(key) list(
    z = lapply(seq_along(sm_timescales), function(s) do.call(rbind, lapply(res, function(p) p[[s]][[key]]$z))),
    m = lapply(seq_along(sm_timescales), function(s) do.call(rbind, lapply(res, function(p) p[[s]][[key]]$method)))))
  names(out) <- names(INDEX_DEFS)
  rm(res); invisible(gc(verbose = FALSE))
  
  coords <- terra::xyFromCell(grid_ref, valid_cells)
  
  for (key in names(INDEX_DEFS)) {
    def <- INDEX_DEFS[[key]]
    out_dir <- file.path(def$root, layer); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    basin_df <- data.frame(date = dates, soil_moisture_m3m3 = basin_sm,
                           valid_area_coverage_fraction = basin_sm_cov)
    S <- character(0)
    add <- function(fmt, ...) S <<- c(S, sprintf(fmt, ...))
    add("============================================================")
    add("%s SUMMARY - ERA5-Land %s (%s)", def$label, layer, lc$desc)
    add("============================================================")
    add("Calculation date    : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    add("Input files         : %s", paste(vapply(lc$layers, function(ly) sprintf("%s (%s, units '%s')",
                                                                                 layer_meta[[ly]]$file, layer_meta[[ly]]$var, layer_meta[[ly]]$units), ""), collapse = "; "))
    add("Depth weights       : %s", paste(sprintf("%s=%.4f", names(w), w), collapse = " "))
    add("Time period         : %s to %s (%d months)", format(min(dates)), format(max(dates)), length(dates))
    add("Grid                : %d x %d native ERA5-Land | valid basin pixels %d", terra::ncol(grid_ref), terra::nrow(grid_ref), n_pix)
    add("Basin area          : %.1f km2 (cell area x fractional coverage, 7-subbasin union)", sm_bw$area_m2 / 1e6)
    add("Basin mean rule     : fractional-coverage area weights, min valid coverage %.2f", BASIN_MIN_VALID_AREA_COVERAGE)
    add("Distribution        : %s", def$distribution)
    add("Reference           : %s", def$reference)
    add("Aggregation         : right-aligned rolling mean; monthly stratification; whole-record reference")
    add("Basin SM climatology (m3/m3): %s", paste(sprintf("%s=%.3f", month.abb, clim), collapse = " "))
    
    for (s in seq_along(sm_timescales)) {
      sc <- sm_timescales[s]; z <- out[[key]]$z[[s]]; meth <- out[[key]]$m[[s]]
      tag <- sprintf("%s_%02d", key, sc)
      cat(sprintf("  -> %s-%d %s outputs (CSV, Excel, NetCDF, map, basin/regional means)\n", def$label, sc, layer))
      
      excel_data <- list()
      for (m in 1:12) {
        idx_m <- which(month_vec == m); if (!length(idx_m)) next
        df <- data.frame(lon = coords[, 1], lat = coords[, 2], z[, idx_m, drop = FALSE])
        names(df)[-(1:2)] <- format(dates[idx_m], "%Y")
        write.csv(df, file.path(out_dir, sprintf("%s_month%02d_%s.csv", tag, m, month.abb[m])), row.names = FALSE)
        excel_data[[month.abb[m]]] <- df
        nc_f <- file.path(out_dir, sprintf("%s_month%02d_%s.nc", tag, m, month.abb[m]))
        write_index_nc(z[, idx_m, drop = FALSE], dates[idx_m], grid_ref, valid_cells, nc_f, key,
                       sprintf("%s-%d, ERA5-Land %s", def$longname, sc, layer))
        if (layer == PRIMARY_LAYER) required_outputs <<- c(required_outputs, nc_f)
        if (key == "ssmi" && WRITE_LEGACY_SAF_FILES)
          write_index_nc(z[, idx_m, drop = FALSE], dates[idx_m], grid_ref, valid_cells,
                         file.path(def$root, sprintf("ssmi_%s_%02d_month%02d_%s.nc", layer, sc, m, month.abb[m])),
                         "ssmi", sprintf("%s-%d, ERA5-Land %s (legacy SAF name)", def$longname, sc, layer))
      }
      xlsx_file <- file.path(out_dir, sprintf("%s_all_months.xlsx", tag))
      if (file.exists(xlsx_file)) file.remove(xlsx_file)
      writexl::write_xlsx(excel_data, xlsx_file)
      plot_method_map(meth, grid_ref, valid_cells, basin,
                      file.path(out_dir, sprintf("%s_distribution_map.png", tag)),
                      sprintf("%s-%d (%s)", def$label, sc, layer), def$method_labels, pct(meth, 2))
      
      # ---- Basin mean of pixel index values (same estimator as Script 4 / 11) ----
      b_val <- basin_weighted_mean(z, sm_bw, cells = valid_cells)
      b_cov <- basin_valid_coverage(z, sm_bw, cells = valid_cells)
      f_avg <- file.path(out_dir, sprintf("%s_basin_averaged_by_month.csv", tag))
      f_cov <- file.path(out_dir, sprintf("%s_basin_valid_area_coverage_by_month.csv", tag))
      utils::write.csv(by_month_wide(b_val, dates), f_avg, row.names = FALSE, na = "")
      utils::write.csv(by_month_wide(b_cov, dates), f_cov, row.names = FALSE, na = "")
      basin_df[[sprintf("%s_%d", toupper(key), sc)]] <- b_val
      n_ev <- plot_basin_index(dates, b_val, def$label, sc, layer,
                               file.path(out_dir, sprintf("%s_basin_timeseries.png", tag)))
      if (is.null(n_ev)) n_ev <- NA_integer_
      
      # ---- Basin + every subbasin (long format, Script 4 layout) ----
      full <- matrix(NA_real_, terra::ncell(grid_ref), length(dates)); full[valid_cells, ] <- z
      reg_tab <- data.table::rbindlist(lapply(names(region_weights), function(id) {
        rw <- region_weights[[id]]
        data.table::data.table(
          index_type = key, layer = layer, scale = as.integer(sc),
          region_id = id, region_name = rw$region_name, date = dates,
          value = regional_weighted_mean(full, rw$weights, min_coverage = BASIN_MIN_VALID_AREA_COVERAGE),
          valid_area_coverage_fraction = regional_valid_coverage(full, rw$weights),
          coverage_area_m2 = rw$area_m2, n_cells = rw$n_cells, n_partial_cells = rw$n_partial)
      }), use.names = TRUE)
      f_reg <- file.path(REGIONAL_DIR, sprintf("%s_%s_%02d_regional_averaged_by_month.csv", key, layer, sc))
      data.table::fwrite(reg_tab, f_reg, na = "NA")
      rm(full)
      required_outputs <<- c(required_outputs, f_avg, f_cov, f_reg, xlsx_file)
      
      # ---- Summary block ----
      add(""); add("---------- %s-%d ----------", def$label, sc)
      add("NA rate (basin pixels) : %.3f%%", 100 * mean(is.na(z)))
      add("Methods: %s=%.1f%%  empirical=%.1f%%  zero-var=%.1f%%", def$method_labels[1], pct(meth, 1), pct(meth, 2), pct(meth, 3))
      add("Pixel %s mean / sd    : %.3f / %.3f", def$label, mean(z, na.rm = TRUE), stats::sd(as.vector(z), na.rm = TRUE))
      add("Clipped at +/-%.2f     : %d dry, %d wet", z_clip, sum(z <= -z_clip, na.rm = TRUE), sum(z >= z_clip, na.rm = TRUE))
      add("Basin drought events   : %s (framework definition, onset < %.2f, end >= %.2f)",
          ifelse(is.na(n_ev), "NA", as.character(n_ev)), DROUGHT_ONSET, DROUGHT_END)
      add("  %-5s %9s %7s %9s %7s %9s %7s", "Month", "Case1(n)", "C1(%)", "Case2(n)", "C2(%)", "Case3(n)", "C3(%)")
      for (mm in 1:12) {
        v <- as.integer(meth[, month_vec == mm]); v <- v[!is.na(v)]; nt <- length(v)
        if (!nt) next
        add("  %-5s %9d %6.1f%% %9d %6.1f%% %9d %6.1f%%", month.abb[mm],
            sum(v == 1L), 100 * sum(v == 1L) / nt, sum(v == 2L), 100 * sum(v == 2L) / nt,
            sum(v == 3L), 100 * sum(v == 3L) / nt)
      }
      add("Class frequency (all pixels, full record):")
      add("  W3 extremely wet  (>= 2.0)       : %.1f%%", 100 * mean(z >= 2, na.rm = TRUE))
      add("  W2 severely wet   (1.5 to < 2.0) : %.1f%%", 100 * mean(z >= 1.5 & z < 2, na.rm = TRUE))
      add("  W1 moderately wet (1.0 to < 1.5) : %.1f%%", 100 * mean(z >= 1 & z < 1.5, na.rm = TRUE))
      add("  N0 normal         (-1 to 1)      : %.1f%%", 100 * mean(z > -1 & z < 1, na.rm = TRUE))
      add("  D1 moderately dry (-1.5 to -1]   : %.1f%%", 100 * mean(z <= -1 & z > -1.5, na.rm = TRUE))
      add("  D2 severely dry   (-2 to -1.5]   : %.1f%%", 100 * mean(z <= -1.5 & z > -2, na.rm = TRUE))
      add("  D3 extremely dry  (<= -2.0)      : %.1f%%", 100 * mean(z <= -2, na.rm = TRUE))
      
      manifest_rows[[length(manifest_rows) + 1L]] <<- data.table::data.table(
        index_type = key, index_label = def$label, layer = layer, scale = as.integer(sc),
        primary_layer = (layer == PRIMARY_LAYER),
        output_directory = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
        seasonal_netcdf_files = 12L, first_date = min(dates), last_date = max(dates),
        n_months = length(dates), n_valid_pixels = n_pix,
        pct_method1 = pct(meth, 1), pct_empirical = pct(meth, 2), pct_zero_var = pct(meth, 3),
        regional_min_valid_area_coverage = BASIN_MIN_VALID_AREA_COVERAGE,
        spatial_grid = "native ERA5-Land precipitation grid")
    }
    
    f_ts <- file.path(out_dir, sprintf("%s_basin_timeseries.csv", key))
    write.csv(basin_df, f_ts, row.names = FALSE, na = "")
    summ <- file.path(out_dir, sprintf("%s_%s_summary.txt", key, layer))
    writeLines(S, summ)
    required_outputs <<- c(required_outputs, f_ts, summ)
    cat(sprintf("  Summary written: %s\n", basename(summ)))
  }
  invisible(TRUE)
}

# ==============================================================================
# RUN
# ==============================================================================
tryCatch({
  for (layer in names(layer_configs)) process_layer_config(layer, cl)
}, finally = {
  if (!is.null(cl)) { stopCluster(cl); cat("\nParallel cluster stopped\n") }
  rm(list = intersect(c(".wk_mat", ".wk_month", ".wk_cfg"), ls(globalenv(), all.names = TRUE)),
     envir = globalenv())
})

# ---- Manifest + final output validation (Script 4 rule) ----
manifest <- data.table::rbindlist(manifest_rows, use.names = TRUE)
manifest_file <- file.path(ssmi_root, "soil_moisture_index_output_manifest.csv")
data.table::fwrite(manifest, manifest_file)
required_outputs <- unique(c(required_outputs, manifest_file))
missing_outputs <- required_outputs[!file.exists(required_outputs) |
                                      is.na(file.info(required_outputs)$size) |
                                      file.info(required_outputs)$size <= 0]
if (length(missing_outputs))
  stop("Final output validation failed; missing/empty product(s):\n",
       paste(missing_outputs, collapse = "\n"), call. = FALSE)

cat("\n============================================================\n")
cat("SSI-SM / SSMI CALCULATION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("SSI-SM output root : %s\n", normalizePath(ssi_sm_root, winslash = "/")))
cat(sprintf("SSMI output root   : %s\n", normalizePath(ssmi_root, winslash = "/")))
cat(sprintf("Regional averages  : %s\n", normalizePath(REGIONAL_DIR, winslash = "/")))
cat(sprintf("Primary layer for Script 11 / SAF pipeline: %s\n", PRIMARY_LAYER))
cat("Naming: <key>_results_seasonal/<layer>/<key>_SS_monthMM_Mon.nc  (key = ssi_sm | ssmi)\n")
if (WRITE_LEGACY_SAF_FILES)
  cat("Legacy SAF names : ssmi_results_seasonal/ssmi_<layer>_SS_monthMM_Mon.nc\n")
cat("Framework output validation passed.\n")