####################################################################################
# 4pr_pet_trends_v2.r  ??  Trend Analysis ??? Precipitation, PET (PM + Thornthwaite), Temperature
# ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# Performs comprehensive VC-MK + TFPW-MK + Spectral + Changepoint trend analysis
# for ALL BASIN PIXELS and SPECIFIC POINTS for four variables plus PET-method
# comparison diagnostics:
#   ??? Precipitation   (ERA5-Land, mm/month ??? mm/year, full-year sum)
#   ??? PET             (Penman-Monteith, mm/month ??? warm-season sum Apr???Oct)
#   ??? PET_Thw         (Thornthwaite, mm/month ??? warm-season sum Apr???Oct)
#   ??? Temperature     (ERA5-Land 2 m Tair, K?????C, annual mean + 12 calendar months)
#
# MODIFICATIONS vs. original:
#   [NEW] Thornthwaite PET (PET_Thw) is processed in parallel with PM PET.
#         Both PET variants go through identical preprocessing, aggregation,
#         pixel-level trend analysis, basin-average analysis, specific-point
#         extraction and trend analysis.
#   [NEW] PET-method comparison products are retained for downstream plotting and
#         interpretation of PET-formulation sensitivity.
#
# OUTPUTS (binary .rds + CSV):
#   all_results.rds              pixel + basin trend statistics (Pr, PET, PET_Thw, Tair)
#   analysis_metadata.rds        time series, extents, parameters (incl. PET_Thw)
#   clipped_template.rds         EPSG:3005 raster template for visualization
#   summary_statistics.csv       per-variable significance counts
#   point_trend_stats.csv        per-point stats (annual + monthly ?? 4 vars)
#   point_monthly_timeseries.csv raw monthly Pr / PET / PET_Thw / Tair at each point
#   data_processing.log
####################################################################################

# rm(list = ls())  # commented out: safe to source from a master script
# gc()             # commented out: safe to source from a master script

library(ncdf4)
library(terra)
library(data.table)
library(Kendall)
library(zoo)
library(sf)
library(changepoint)
library(trend)
library(modifiedmk)

# ?????? Shared statistical functions (VC-MK, TFPW-MK, PELT, spectral) ??????????????????????????????????????????
# check_min_value_threshold(), calculate_variance_with_ties(),
# detect_regime_shift_pelt(), and mk_tfpw_spectral_for_series() are defined in
# DROUGHT_ANALYSIS_utils.R (Section J).  Parameters max_min_value_pct_precip,
# max_min_value_pct_pet, min_positive_value, and min_obs_changepoint are passed
# explicitly to mk_tfpw_spectral_for_series() from the local variables set below.
# ================= USER / ENV =================
# ================= USER / ENV =================
WD_PATH <- Sys.getenv(
  "NECHAKO_WD",
  unset = "D:/Nechako_Drought/Nechako/"
)
if (!dir.exists(WD_PATH)) {
  stop("Working directory not found: ", WD_PATH)
}
setwd(WD_PATH)
source(file.path(WD_PATH, "DROUGHT_ANALYSIS_utils_v2.R"))

out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ===== LOGGING SETUP =====
LOG_FILE <- file.path(out_dir, "data_processing.log")
log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, "   ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, "   ", msg))
}
cat(
  paste0(
    "\n============================================================\n",
    "Trend Analysis Data Processing - Started\n",
    "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "============================================================\n"
  ),
  file = LOG_FILE,
  append = TRUE
)
##
library(future)
library(future.apply)
num_cores <- min(
  max(1L, parallel::detectCores(logical = TRUE) - 1L),
  8L
)

future::plan(
  future::multisession,
  workers = num_cores
)

log_event(paste("Using", num_cores, "future multisession workers"))

####################################################################################
# ?????? SPECIFIC POINTS CONFIG
####################################################################################
SPECIFIC_PTS <- data.frame(
  x = c(-124.5, -123.0, -125.0),
  y = c(  54.0,   53.8,   54.2)
)

####################################################################################
# ?????? PARAMETERS
####################################################################################
alpha                    <- 0.05
n_sim_spectral           <- 500
max_tie_percent          <- 50
max_min_value_pct_precip <- 50
max_min_value_pct_pet    <- 50
min_positive_value       <- 0.01
min_obs_changepoint      <- 20

warm_months_pet          <- 4:10
min_nonzero_pet_monthly  <- 10
TEMP_START_YR            <- 1950L

# ===== INPUT PATHS =====
precip_path   <- "monthly_data_direct/total_precipitation_monthly.nc"
pet_path      <- "monthly_data_direct/potential_evapotranspiration_monthly.nc"
# [NEW] Thornthwaite PET ??? produced by 2preq_PET_ERALand.R (Section 5b)
pet_thw_path  <- "monthly_data_direct/potential_evapotranspiration_thornthwaite_monthly.nc"

TAIR_NC_CANDIDATES <- c(
  "monthly_data_direct/2m_temperature_monthly.nc",
  "monthly_data_direct/t2m_monthly.nc",
  "monthly_data_direct/air_temperature_2m_monthly.nc",
  "monthly_data_direct/ERA5Land_T2m_monthly.nc"
)

log_event(paste("Using", num_cores, "cores for parallel processing"))

####################################################################################
# ?????? BASIN BOUNDARY
####################################################################################
log_event("Searching for Nechako Basin boundary file...")
basin_boundary <- NULL
basin_files <- c(
  "Spatial/nechakoBound_dissolve.kmz",
  "nechakoBound_dissolve.kmz",
  "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.kmz",
  "nechako_basin.shp", "Nechako_Basin.shp", "basin_boundary.shp",
  "../nechako_basin.shp", "data/nechako_basin.shp",
  "D:/Nechako_Drought/Nechako/Spatial/nechakoBound_dissolve.shp"
)
for (bf in basin_files) {
  if (file.exists(bf)) {
    tryCatch({
      if (tolower(tools::file_ext(bf)) == "kmz") {
        kml_file <- unzip(bf, exdir = tempdir())[1]
        basin_boundary <- st_read(kml_file, quiet = TRUE)
        unlink(kml_file)
      } else {
        basin_boundary <- st_read(bf, quiet = TRUE)
      }
      if (nrow(basin_boundary) > 1L)
        basin_boundary <- sf::st_as_sf(sf::st_union(basin_boundary))
      basin_boundary <- st_transform(basin_boundary, "EPSG:3005")
      log_event(paste("??? Loaded basin boundary from:", bf))
      break
    }, error = function(e) log_event(paste("  Error loading", bf, ":", e$message)))
  }
}
if (is.null(basin_boundary))
  stop("CRITICAL: Nechako Basin boundary NOT FOUND.")
saveRDS(basin_boundary, file.path(out_dir, "basin_boundary.rds"))

####################################################################################
# ?????? HELPER FUNCTIONS
####################################################################################



aggregate_to_annual_fast <- function(monthly_matrix, years, method = "sum") {
  y_levels      <- unique(years)
  yfac          <- match(years, y_levels)
  count_by_year <- rowsum((!is.na(monthly_matrix)) * 1L, group = yfac, reorder = FALSE)
  m0            <- monthly_matrix; m0[is.na(m0)] <- 0
  sum_by_year   <- rowsum(m0, group = yfac, reorder = FALSE)
  res           <- if (identical(method, "sum")) sum_by_year else
    sum_by_year / pmax(count_by_year, 1)
  res[count_by_year < 6] <- NA_real_
  rownames(res) <- y_levels
  res
}

aggregate_to_annual_warmseason_fast <- function(monthly_matrix, years, months_vec,
                                                warm_months = 4:10, method = "sum") {
  warm_idx    <- which(months_vec %in% warm_months)
  if (!length(warm_idx)) stop("No warm-season months found in months_vec.")
  warm_matrix <- monthly_matrix[warm_idx, , drop = FALSE]
  warm_years  <- years[warm_idx]
  y_levels    <- unique(warm_years)
  yfac        <- match(warm_years, y_levels)
  count       <- rowsum((!is.na(warm_matrix)) * 1L, group = yfac, reorder = FALSE)
  m0          <- warm_matrix; m0[is.na(m0)] <- 0
  sums        <- rowsum(m0, group = yfac, reorder = FALSE)
  res         <- if (identical(method, "sum")) sums else sums / pmax(count, 1)
  res[count < ceiling(length(warm_months) * 0.5)] <- NA_real_
  rownames(res) <- y_levels
  res
}

compute_month_index <- function(months) split(seq_along(months), months)

####################################################################################
# ?????? ROOT-CAUSE DIAGNOSTICS / QUALITY SCREEN
####################################################################################
diag_store <- new.env(parent = emptyenv())
diag_store$qa_tables <- list()

series_quality_profile <- function(x, floor_value = NA_real_) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) {
    return(list(
      n = 0L, sd = NA_real_, uniq = 0L, zero_frac = NA_real_,
      floor_frac = NA_real_, min = NA_real_, median = NA_real_, max = NA_real_
    ))
  }
  list(
    n = n,
    sd = if (n >= 2L) stats::sd(x) else NA_real_,
    uniq = length(unique(x)),
    zero_frac = mean(x == 0),
    floor_frac = if (is.finite(floor_value)) mean(x <= floor_value) else NA_real_,
    min = min(x),
    median = stats::median(x),
    max = max(x)
  )
}

log_series_quality <- function(x, label, floor_value = NA_real_) {
  prof <- series_quality_profile(x, floor_value = floor_value)
  floor_txt <- if (is.finite(prof$floor_frac)) sprintf("%.2f%%", 100 * prof$floor_frac) else "NA"
  log_event(sprintf(
    "  [QA:%s] n=%d | sd=%s | unique=%d | zero_frac=%s | floor_frac=%s | min=%.4g | median=%.4g | max=%.4g",
    label,
    prof$n,
    ifelse(is.na(prof$sd), "NA", format(signif(prof$sd, 4), scientific = TRUE)),
    prof$uniq,
    ifelse(is.na(prof$zero_frac), "NA", sprintf("%.2f%%", 100 * prof$zero_frac)),
    floor_txt,
    prof$min, prof$median, prof$max
  ))
  invisible(prof)
}

matrix_quality_profile <- function(mat, label, floor_value = NA_real_) {
  if (is.null(mat) || !length(mat)) {
    log_event(sprintf("  [QA:%s] matrix empty", label))
    return(invisible(NULL))
  }
  
  n_time   <- nrow(mat)
  n_series <- ncol(mat)
  
  sd_vals <- vapply(seq_len(n_series), function(i) {
    x <- mat[, i]
    x <- x[is.finite(x)]
    if (length(x) < 2L) return(NA_real_)
    stats::sd(x)
  }, numeric(1))
  
  uniq_vals <- vapply(seq_len(n_series), function(i) {
    x <- mat[, i]
    length(unique(x[is.finite(x)]))
  }, integer(1))
  
  zero_frac <- vapply(seq_len(n_series), function(i) {
    x <- mat[, i]
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    mean(x == 0)
  }, numeric(1))
  
  floor_frac <- if (is.finite(floor_value)) {
    vapply(seq_len(n_series), function(i) {
      x <- mat[, i]
      x <- x[is.finite(x)]
      if (!length(x)) return(NA_real_)
      mean(x <= floor_value)
    }, numeric(1))
  } else {
    rep(NA_real_, n_series)
  }
  
  out <- data.table(
    label = label,
    n_time = n_time,
    n_series = n_series,
    pct_na = 100 * sum(!is.finite(mat)) / length(mat),
    median_sd = stats::median(sd_vals, na.rm = TRUE),
    n_zero_sd = sum(sd_vals <= 1e-6, na.rm = TRUE),
    n_low_sd = sum(sd_vals < 0.01, na.rm = TRUE),
    median_unique = stats::median(uniq_vals, na.rm = TRUE),
    n_very_tied = sum(uniq_vals <= 3L, na.rm = TRUE),
    median_zero_frac = stats::median(zero_frac, na.rm = TRUE),
    median_floor_frac = stats::median(floor_frac, na.rm = TRUE),
    n_floor_dominant = sum(floor_frac > 0.5, na.rm = TRUE)
  )
  
  floor_txt <- if (is.finite(floor_value)) sprintf(" | floor>%.3g=%d", floor_value, out$n_floor_dominant) else ""
  log_event(sprintf(
    "  [QA:%s] series=%d time=%d | NA=%.2f%% | median sd=%s | sd<1e-6=%d | sd<0.01=%d | median unique=%s | uniq<=3=%d | median zero=%.2f%%%s",
    label,
    n_series, n_time, out$pct_na,
    ifelse(is.na(out$median_sd), "NA", format(signif(out$median_sd, 4), scientific = TRUE)),
    out$n_zero_sd, out$n_low_sd,
    ifelse(is.na(out$median_unique), "NA", format(out$median_unique)),
    out$n_very_tied,
    100 * out$median_zero_frac,
    floor_txt
  ))
  
  diag_store$qa_tables[[length(diag_store$qa_tables) + 1L]] <- out
  invisible(out)
}

log_result_root_cause_summary <- function(results_dt) {
  if (is.null(results_dt) || !nrow(results_dt)) {
    log_event("  [DIAG] No results table available for root-cause summary.")
    return(invisible(NULL))
  }
  
  nonfinite_vc  <- results_dt[!filtered_vc & !is.finite(p_value_vc), .N]
  nonfinite_tf  <- results_dt[!filtered_tfpw & !is.finite(p_value_tfpw), .N]
  
  log_event(sprintf(
    "  [DIAG] non-finite p-values: VC=%d | TFPW=%d",
    nonfinite_vc, nonfinite_tf
  ))
  
  reason_counts <- rbindlist(list(
    
    results_dt[
      !is.na(filter_reason_vc) & filter_reason_vc != "none",
      .N,
      by = .(variable, period, is_basin_average, reason = filter_reason_vc)
    ][, method := "VC"],
    
    results_dt[
      !is.na(filter_reason_tfpw) & filter_reason_tfpw != "none",
      .N,
      by = .(variable, period, is_basin_average, reason = filter_reason_tfpw)
    ][, method := "TFPW"]
    
  ))
  
  setnames(reason_counts, "N", "n")
  
  if (nrow(reason_counts)) {
    reason_counts <- reason_counts[order(-n)]
    top_n <- reason_counts[seq_len(min(8L, nrow(reason_counts)))]
    log_event("  [DIAG] dominant filter reasons (top 8 rows):")
    for (i in seq_len(nrow(top_n))) {
      rr <- top_n[i]
      log_event(sprintf(
        "    [%s/%s/%s] %s = %d",
        rr$variable, rr$period, rr$method, rr$reason, rr$n
      ))
    }
    fwrite(reason_counts, file.path(out_dir, "root_cause_filter_reasons.csv"))
  } else {
    log_event("  [DIAG] no filtered series recorded in results.")
  }
  
  diag_summary <- results_dt[, .(
    n_total = .N,
    n_filtered_vc = sum(filtered_vc, na.rm = TRUE),
    n_filtered_tfpw = sum(filtered_tfpw, na.rm = TRUE),
    n_nonfinite_p_vc = sum(!filtered_vc & !is.finite(p_value_vc), na.rm = TRUE),
    n_nonfinite_p_tfpw = sum(!filtered_tfpw & !is.finite(p_value_tfpw), na.rm = TRUE),
    n_vc_corrected = sum(vc_corrected %in% TRUE, na.rm = TRUE),
    n_tfpw_applied = sum(tfpw_applied %in% TRUE, na.rm = TRUE)
  ), by = .(variable, period, is_basin_average)]
  
  fwrite(diag_summary, file.path(out_dir, "root_cause_summary.csv"))
  invisible(diag_summary)
}

####################################################################################
# ?????? CORE TREND ANALYSIS FUNCTION
####################################################################################

####################################################################################
# ?????? DATA LOADING & PREPROCESSING
####################################################################################
log_event("Loading Precipitation, PET (PM and Thornthwaite), and Temperature data...")

precip_full <- rast(precip_path)
pet_full    <- rast(pet_path)

# [NEW] Thornthwaite PET
if (!file.exists(pet_thw_path))
  stop("Thornthwaite PET NetCDF not found: ", pet_thw_path,
       "\nRun 2preq_PET_ERALand.R first.")
pet_thw_full <- rast(pet_thw_path)
log_event(sprintf("  [PET_Thw] Loaded: %d layers from %s", nlyr(pet_thw_full), basename(pet_thw_path)))

# ?????? Temperature ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
tair_nc <- NULL
for (cand in TAIR_NC_CANDIDATES) {
  if (file.exists(cand)) { tair_nc <- cand; break }
}
if (is.null(tair_nc))
  stop("ERA5-Land 2 m temperature NetCDF not found. Tried:\n  ",
       paste(TAIR_NC_CANDIDATES, collapse="\n  "))
log_event(paste("  [temp] Using:", basename(tair_nc)))
tair_full <- rast(tair_nc)

# ?????? Date extraction ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
n_layers   <- nlyr(precip_full)
dates_full <- tryCatch({
  d <- as.Date(terra::time(precip_full))
  if (is.null(d) || length(d) != n_layers || all(is.na(d))) NULL else d
}, error=function(e) NULL)
if (is.null(dates_full)) {
  dates_full <- seq(as.Date("1950-01-01"), by="month", length.out=n_layers)
  log_event(paste("Generated", n_layers, "monthly dates from 1950"))
}
dates  <- dates_full
years  <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))

# ?????? Unit conversions ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
precip_full <- precip_full * 1000    # m ??? mm/day
month_days_base  <- c(31,28,31,30,31,30,31,31,30,31,30,31)
is_leap <- function(yr) (yr%%4==0) & (yr%%100!=0 | yr%%400==0)
days_in_month <- month_days_base[months]
days_in_month[is_leap(years) & months==2] <- 29

log_event("Converting Pr, PET (PM), and PET (Thw) from mm/day to mm/month...")
for (i in seq_len(nlyr(precip_full)))   precip_full[[i]]   <- precip_full[[i]]   * days_in_month[i]
for (i in seq_len(nlyr(pet_full)))      pet_full[[i]]      <- pet_full[[i]]      * days_in_month[i]
for (i in seq_len(nlyr(pet_thw_full)))  pet_thw_full[[i]]  <- pet_thw_full[[i]]  * days_in_month[i]

# Temperature: Kelvin ??? ??C
tair_sample <- as.numeric(terra::global(tair_full[[1]], "mean", na.rm=TRUE))
if (!is.na(tair_sample) && tair_sample > 200) {
  log_event("  [temp] Converting K ??? ??C...")
  tair_full <- tair_full - 273.15
}

n_tair  <- nlyr(tair_full)
dates_t <- tryCatch({
  d <- as.Date(terra::time(tair_full))
  if (is.null(d) || length(d) != n_tair || all(is.na(d))) NULL else d
}, error=function(e) NULL)
if (is.null(dates_t)) {
  dates_t  <- seq(as.Date(sprintf("%d-01-01", TEMP_START_YR)), by="month", length.out=n_tair)
  log_event("  [temp] terra::time() invalid ??? using positional date sequence")
}
years_t  <- as.integer(format(dates_t, "%Y"))
months_t <- as.integer(format(dates_t, "%m"))
log_event(sprintf("  [temp] Date range: %s to %s", min(dates_t), max(dates_t)))

####################################################################################
# ?????? REPROJECT & CLIP TO BASIN
####################################################################################
log_event("Reprojecting to BC Albers (EPSG:3005)...")
precip_full   <- project(precip_full,   "EPSG:3005", method="bilinear")
pet_full      <- project(pet_full,      "EPSG:3005", method="bilinear")
pet_thw_full  <- project(pet_thw_full,  "EPSG:3005", method="bilinear")  # [NEW]
tair_full     <- project(tair_full,     "EPSG:3005", method="bilinear")

# Use the exact grid generated by 3SPEI_ERALand_v2.R. This removes the previous
# cross-pipeline geometry ambiguity: raw climate trends and drought-index trends
# are now evaluated on one canonical EPSG:3005 grid and one basin mask.
spei_grid_file <- file.path("spei_results_seasonal", "spei_analysis_grid.rds")
if (!file.exists(spei_grid_file)) {
  stop("Canonical SPEI analysis grid not found: ", spei_grid_file,
       "\nRun 3SPEI_ERALand_v2.R before 4pr_pet_trends_v2.r.")
}
spei_grid <- readRDS(spei_grid_file)

# Unwrap immediately — these were saved with terra::wrap() in script 3,
# so at this point they're PackedSpatRaster (S4), not usable rasters yet.
spei_grid$raster_template <- terra::unwrap(spei_grid$raster_template)
spei_grid$basin_mask      <- terra::unwrap(spei_grid$basin_mask)

raw_template <- terra::unwrap(spei_grid$raster_template)
analysis_template <- tryCatch({
  if (inherits(raw_template, "PackedSpatRaster")) {
    terra::unwrap(raw_template)
  } else {
    terra::nlyr(raw_template) # Test pointer validity
    raw_template
  }
}, error = function(e) {
  # Fallback to loaded precipitation grid if pointer is expired
  precip_full[[1]]
})

if (is.null(analysis_template) || is.null(spei_grid$basin_mask)) {
  stop("Invalid spei_analysis_grid.rds: missing template or basin_mask.")
}
log_event(sprintf("Resampling raw ERA5-Land variables to canonical SPEI grid: %d x %d, res %.1f x %.1f m",
                  nrow(analysis_template), ncol(analysis_template),
                  res(analysis_template)[1], res(analysis_template)[2]))
precip_full  <- resample(precip_full,  analysis_template, method="bilinear")
pet_full     <- resample(pet_full,     analysis_template, method="bilinear")
pet_thw_full <- resample(pet_thw_full, analysis_template, method="bilinear")
tair_full    <- resample(tair_full,    analysis_template, method="bilinear")

log_event("Applying canonical Nechako Basin mask...")
# as.logical()/as.numeric() on a SpatRaster return another SpatRaster (cast to
# that type), not a plain vector - and length() on a SpatRaster is nlyr(), not
# ncell(). Pull real cell values out first so the length check and the
# values(mask_r) <- assignment below both operate on an actual vector.
analysis_mask <- as.logical(terra::values(spei_grid$basin_mask, mat = FALSE))
if (length(analysis_mask) != ncell(analysis_template))
  stop("Canonical basin mask length does not match canonical grid cell count.")
mask_r <- analysis_template
values(mask_r) <- as.numeric(analysis_mask)
precip_clipped  <- mask(precip_full,  mask_r, maskvalues = 0, updatevalue = NA)
pet_clipped     <- mask(pet_full,     mask_r, maskvalues = 0, updatevalue = NA)
pet_thw_clipped <- mask(pet_thw_full, mask_r, maskvalues = 0, updatevalue = NA)
tair_clipped    <- mask(tair_full,    mask_r, maskvalues = 0, updatevalue = NA)

clipped_template <- rast(analysis_template, nlyrs=1); values(clipped_template) <- NA_real_
saveRDS(clipped_template, file.path(out_dir, "clipped_template.rds"))
# Save a real canonical mask separately; clipped_template.rds is intentionally
# an empty geometry shell used by the visualization code.
canonical_mask <- list(
  template = clipped_template,
  basin_mask = analysis_mask,
  basin_pixels = sum(analysis_mask),
  ncell = ncell(analysis_template),
  crs = crs(analysis_template),
  extent = as.vector(ext(analysis_template)),
  resolution = res(analysis_template),
  source = "3SPEI_ERALand_v2.R canonical grid",
  created = Sys.time()
)
saveRDS(canonical_mask, file.path(out_dir, "canonical_basin_mask.rds"), compress=TRUE)
log_event(sprintf("??? Saved clipped template (EPSG:3005, %d x %d, res: %.0f m)",
                  nrow(clipped_template), ncol(clipped_template), res(clipped_template)[1]))

n_bbox_cells     <- ncell(precip_clipped)
first_layer_vals <- values(precip_clipped[[1]], mat=FALSE)
valid_mask       <- !is.na(first_layer_vals)
n_basin_pixels   <- sum(valid_mask)
reduction_pct    <- 100 * (1 - n_basin_pixels / n_bbox_cells)
if (!identical(as.integer(n_basin_pixels), as.integer(sum(analysis_mask))))
  stop("Canonical mask/data mismatch: valid data cells do not equal canonical basin pixels.")
log_event(sprintf("CANONICAL BASIN GRID: %d cells | %d basin pixels (%.1f%% outside mask)",
                  n_bbox_cells, n_basin_pixels, reduction_pct))
if (n_basin_pixels == 0) stop("CRITICAL: No valid cells after basin clipping.")

####################################################################################
# ?????? BASIN AVERAGE TIME SERIES
####################################################################################
log_event("Computing basin-averaged monthly time series (all variables)...")
precip_monthly_avg    <- as.vector(global(precip_clipped,  fun="mean", na.rm=TRUE)[,1])
pet_monthly_avg       <- as.vector(global(pet_clipped,     fun="mean", na.rm=TRUE)[,1])
pet_thw_monthly_avg   <- as.vector(global(pet_thw_clipped, fun="mean", na.rm=TRUE)[,1])  # [NEW]
tair_monthly_avg      <- as.vector(global(tair_clipped,    fun="mean", na.rm=TRUE)[,1])

# Annual aggregation
precip_annual_avg_matrix   <- aggregate_to_annual_fast(
  matrix(precip_monthly_avg, ncol=1), years, method="sum")
pet_annual_avg_matrix      <- aggregate_to_annual_warmseason_fast(
  matrix(pet_monthly_avg, ncol=1), years, months_vec=months,
  warm_months=warm_months_pet, method="sum")
pet_thw_annual_avg_matrix  <- aggregate_to_annual_warmseason_fast(          # [NEW]
  matrix(pet_thw_monthly_avg, ncol=1), years, months_vec=months,
  warm_months=warm_months_pet, method="sum")
tair_annual_avg_matrix     <- aggregate_to_annual_fast(
  matrix(tair_monthly_avg, ncol=1), years_t, method="mean")

precip_annual_avg    <- precip_annual_avg_matrix[,1]
pet_annual_avg       <- pet_annual_avg_matrix[,1]
pet_thw_annual_avg   <- pet_thw_annual_avg_matrix[,1]  # [NEW]
tair_annual_avg      <- tair_annual_avg_matrix[,1]
annual_years         <- as.integer(rownames(precip_annual_avg_matrix))
tair_annual_years    <- as.integer(rownames(tair_annual_avg_matrix))

####################################################################################
# ?????? SANITY CHECK: Thornthwaite PET 0??C-floor crossing-frequency trend  [NEW]
####################################################################################
# WHY THIS EXISTS:
#   Thornthwaite PET is floored at 0 for any month with mean temperature <= 0??C.
#   QA upstream shows PET_Thw is 100% zero at all basin pixels in Jan/Dec and
#   >90% zero in Feb/Nov. That means winter SPEI_Thw at short timescales (e.g.
#   SPEI_Thw-01) effectively collapses to standardized precipitation alone,
#   with no PET contribution -- UNLESS a pixel's winter temperature crosses
#   above 0??C, in which case PET_Thw suddenly becomes nonzero and starts
#   contributing to the index again. As the basin warms over the 76-year
#   record, the number of winter months crossing that 0??C threshold can
#   change, which would inject a spurious trend into short-timescale winter
#   SPEI_Thw that has nothing to do with a genuine PET-method (Thornthwaite vs
#   Penman-Monteith) difference. This is a plausible explanation for the
#   isolated SPEI_Thw-01 significance spike (665/865 pixels) seen in Script 6,
#   versus SPEI-01 PM (34/865) and versus SPEI_Thw at longer scales.
#
# WHAT THIS CHECKS:
#   For each calendar month, compute the basin-wide fraction of pixels with
#   PET_Thw == 0 (the "0??C floor") in each of the record's years, then run a
#   variance-corrected Mann-Kendall trend test on that annual fraction series.
#   A significant DEcreasing trend in the floor-crossing fraction for a winter
#   month (i.e. fewer pixels floored at 0 over time = more months crossing
#   above 0??C) is the fingerprint of the artifact described above, and is
#   flagged here so it can be checked before reporting any SPEI_Thw-01 vs.
#   SPEI-01 PET-method contrast as a genuine finding.
#
# NOTE: this uses pet_thw_clipped values captured BEFORE the zero -> 
#   min_positive_value replacement below, so true 0??C-floor months are
#   preserved for this diagnostic.
####################################################################################
log_event("Running sanity check: PET_Thw 0C-floor crossing-frequency trend by calendar month...")
pet_thw_matrix_raw <- t(values(pet_thw_clipped)[valid_mask, , drop=FALSE])  # months x pixels, true zeros intact

floor_crossing_trend <- tryCatch({
  month_rows <- lapply(1:12, function(mo) which(months == mo))
  yrs_all    <- sort(unique(years))
  
  out_rows <- vector("list", 12L)
  for (mo in 1:12) {
    rows_mo <- month_rows[[mo]]
    if (!length(rows_mo)) next
    yr_mo   <- years[rows_mo]
    # Annual basin-wide fraction of pixels floored at exactly 0 in this calendar month
    frac_by_year <- vapply(yrs_all, function(yr) {
      idx <- rows_mo[yr_mo == yr]
      if (!length(idx)) return(NA_real_)
      vals <- pet_thw_matrix_raw[idx, , drop = FALSE]
      mean(vals == 0, na.rm = TRUE)
    }, numeric(1))
    
    trend <- .mk_vc_single(frac_by_year)
    out_rows[[mo]] <- data.table(
      month             = mo,
      mean_zero_frac    = mean(frac_by_year, na.rm = TRUE),
      first_year_frac   = frac_by_year[1],
      last_year_frac    = frac_by_year[length(frac_by_year)],
      tau               = trend$tau,
      p_value           = trend$pval,
      sen_slope_per_yr  = trend$slope
    )
  }
  data.table::rbindlist(out_rows, fill = TRUE)
}, error = function(e) {
  log_event(sprintf("  ??? floor-crossing trend check failed: %s", conditionMessage(e)))
  NULL
})

if (!is.null(floor_crossing_trend) && nrow(floor_crossing_trend) > 0) {
  fwrite(floor_crossing_trend, file.path(out_dir, "pet_thw_0C_floor_crossing_trend.csv"))
  log_event("  ??? Saved pet_thw_0C_floor_crossing_trend.csv")
  
  winter_months <- c(11L, 12L, 1L, 2L, 3L)
  flagged <- floor_crossing_trend[
    month %in% winter_months & !is.na(p_value) & p_value < 0.05 & tau < 0
  ]
  if (nrow(flagged) > 0) {
    log_event(sprintf(
      paste0("  ?????  CAUTION: significant DEcreasing trend in PET_Thw 0C-floor ",
             "fraction for winter month(s) %s (p<0.05). This is consistent with ",
             "warming-driven change in freeze-floor crossing frequency, which can ",
             "inject a spurious trend into short-timescale winter SPEI_Thw. ",
             "Treat any SPEI_Thw-01 vs SPEI-01 PET-method contrast in these ",
             "months as provisional pending review of this diagnostic."),
      paste(sprintf("%02d", flagged$month), collapse = ", ")
    ))
  } else {
    log_event("  No significant winter floor-crossing trend detected at p<0.05 (see CSV for full month-by-month detail).")
  }
} else {
  log_event("  ??? floor-crossing trend check produced no rows; skipped.")
}

# Pre-process PET (both): replace zeros with min_positive_value
for (clip_r in list(pet_clipped, pet_thw_clipped)) {
  pv  <- values(clip_r)
  nz  <- sum(pv == 0, na.rm=TRUE)
  pv[pv == 0] <- min_positive_value
  values(clip_r) <- pv
  log_event(sprintf("  Replaced %d PET zeros with %.3f mm (%s)",
                    nz, min_positive_value, names(clip_r[[1]])))
}

# Pixel matrices (time ?? pixels)
valid_cell_indices <- which(valid_mask)
valid_xy           <- xyFromCell(precip_clipped, valid_cell_indices)
precip_matrix      <- t(values(precip_clipped)[valid_mask,    , drop=FALSE])
pet_matrix         <- t(values(pet_clipped)[valid_mask,       , drop=FALSE])
pet_thw_matrix     <- t(values(pet_thw_clipped)[valid_mask,   , drop=FALSE])  # [NEW]
tair_matrix        <- t(values(tair_clipped)[valid_mask,      , drop=FALSE])

coords_dt <- data.table(
  space_idx=seq_len(n_basin_pixels), x=valid_xy[,1], y=valid_xy[,2],
  is_basin_average=FALSE)
log_event(sprintf("  Matrices: Pr/PET %d??%d | PET_Thw %d??%d | Tair %d??%d",
                  nrow(precip_matrix), ncol(precip_matrix),
                  nrow(pet_thw_matrix), ncol(pet_thw_matrix),
                  nrow(tair_matrix), ncol(tair_matrix)))

month_index_list   <- compute_month_index(months)
month_index_list_t <- compute_month_index(months_t)

####################################################################################
# ?????? HELPER: BUILD RESULT DATA.TABLE ROW SET (annual + 12 monthly)
####################################################################################
unpack_results <- function(res_list, which = c("vc","tf")) {
  which <- match.arg(which)
  r <- lapply(res_list, `[[`, which)
  data.table(
    tau       = vapply(r, `[[`, numeric(1),   "tau"),
    p         = vapply(r, `[[`, numeric(1),   "p"),
    sl        = vapply(r, `[[`, numeric(1),   "sl"),
    S         = vapply(r, `[[`, numeric(1),   "S"),
    varS      = vapply(r, `[[`, numeric(1),   "varS"),
    n         = vapply(r, `[[`, numeric(1),   "n"),
    rho1      = vapply(r, `[[`, numeric(1),   "rho1"),
    n_ties    = vapply(r, `[[`, numeric(1),   "n_ties"),
    pct_ties  = vapply(r, `[[`, numeric(1),   "percent_ties"),
    n_min     = vapply(r, `[[`, numeric(1),   "n_min_vals"),
    pct_min   = vapply(r, `[[`, numeric(1),   "percent_min_vals"),
    tau_b_adj = vapply(r, `[[`, logical(1),   "tau_b_adjusted"),
    filtered  = vapply(r, `[[`, logical(1),   "filtered"),
    reason    = vapply(r, `[[`, character(1), "reason"),
    extra     = if (which=="vc")
      vapply(r, function(z) z$vc_corrected,  logical(1)) else
        vapply(r, function(z) z$tfpw_applied,  logical(1))
  )
}

build_annual_dt <- function(res_annual, coords_dt, var_name) {
  vc_a   <- unpack_results(res_annual, "vc")
  tf_a   <- unpack_results(res_annual, "tf")
  spec_a <- data.table(
    n_spectral_peaks    = vapply(lapply(res_annual,`[[`,"spec"), `[[`, integer(1), "n_peaks"),
    dominant_period     = vapply(lapply(res_annual,`[[`,"spec"), `[[`, numeric(1), "dominant_period"),
    spectral_confidence = vapply(lapply(res_annual,`[[`,"spec"), `[[`, numeric(1), "conf"))
  cpt_a  <- data.table(
    changepoint_detected   = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, logical(1),   "changepoint_detected"),
    changepoint_position   = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, integer(1),   "changepoint_position"),
    n_changepoints         = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, integer(1),   "n_changepoints"),
    first_changepoint_year = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, integer(1),   "first_changepoint_year"),
    mean_before_shift      = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, numeric(1),   "mean_before"),
    mean_after_shift       = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, numeric(1),   "mean_after"),
    magnitude_shift        = vapply(lapply(res_annual,`[[`,"cpt"), `[[`, numeric(1),   "magnitude_shift"))
  cbind(coords_dt,
        variable="annual", period="annual", month=NA_integer_, variable_name=var_name,
        data.table(tau_vc=vc_a$tau, p_value_vc=vc_a$p, sl_vc=vc_a$sl,
                   vc_corrected=vc_a$extra, n_ties_vc=vc_a$n_ties,
                   percent_ties_vc=vc_a$pct_ties, n_min_vals_vc=vc_a$n_min,
                   percent_min_vals_vc=vc_a$pct_min, tau_b_adjusted_vc=vc_a$tau_b_adj,
                   filtered_vc=vc_a$filtered, filter_reason_vc=vc_a$reason),
        data.table(tau_tfpw=tf_a$tau, p_value_tfpw=tf_a$p, sl_tfpw=tf_a$sl,
                   tfpw_applied=tf_a$extra, n_ties_tfpw=tf_a$n_ties,
                   percent_ties_tfpw=tf_a$pct_ties, n_min_vals_tfpw=tf_a$n_min,
                   percent_min_vals_tfpw=tf_a$pct_min, tau_b_adjusted_tfpw=tf_a$tau_b_adj,
                   filtered_tfpw=tf_a$filtered, filter_reason_tfpw=tf_a$reason),
        n=vc_a$n, rho1_vc=vc_a$rho1, rho1_tfpw=vc_a$rho1,
        spec_a, cpt_a,
        same_significance=(vc_a$p < alpha) == (tf_a$p < alpha),
        same_direction   =sign(vc_a$tau) == sign(tf_a$tau))
}

####################################################################################
# ?????? BASIN AVERAGE TREND PROCESSING
####################################################################################
process_basin_series <- function(ts_monthly, ts_annual, var_name,
                                 is_precip, skip_min_filter = FALSE,
                                 m_index_list = month_index_list) {
  conf_env <- new.env(parent=emptyenv())
  log_series_quality(ts_annual, sprintf("%s basin annual", var_name),
                     floor_value = if (is_precip) 0 else min_positive_value)
  log_series_quality(ts_monthly, sprintf("%s basin monthly", var_name),
                     floor_value = if (is_precip) 0 else min_positive_value)
  res_ann  <- mk_tfpw_spectral_for_series(ts_annual, is_precip, alpha, max_tie_percent,
                                          n_sim_spectral, conf_env,
                                          start_year      = min(annual_years),
                                          skip_min_filter = skip_min_filter)
  vc_a  <- res_ann$vc;  tf_a  <- res_ann$tf
  spec_a <- res_ann$spec; cpt_a <- res_ann$cpt
  
  annual_dt <- data.table(
    space_idx=NA_integer_, x=NA_real_, y=NA_real_,
    variable=var_name, period="annual", month=NA_integer_, variable_name=var_name,
    tau_vc=vc_a$tau, p_value_vc=vc_a$p, sl_vc=vc_a$sl, vc_corrected=vc_a$vc_corrected,
    n_ties_vc=vc_a$n_ties, percent_ties_vc=vc_a$percent_ties,
    n_min_vals_vc=vc_a$n_min_vals, percent_min_vals_vc=vc_a$percent_min_vals,
    tau_b_adjusted_vc=vc_a$tau_b_adjusted, filtered_vc=vc_a$filtered,
    filter_reason_vc=vc_a$reason,
    tau_tfpw=tf_a$tau, p_value_tfpw=tf_a$p, sl_tfpw=tf_a$sl, tfpw_applied=tf_a$tfpw_applied,
    n_ties_tfpw=tf_a$n_ties, percent_ties_tfpw=tf_a$percent_ties,
    n_min_vals_tfpw=tf_a$n_min_vals, percent_min_vals_tfpw=tf_a$percent_min_vals,
    tau_b_adjusted_tfpw=tf_a$tau_b_adjusted, filtered_tfpw=tf_a$filtered,
    filter_reason_tfpw=tf_a$reason,
    n=NA_integer_, rho1_vc=NA_real_, rho1_tfpw=NA_real_,
    n_spectral_peaks=spec_a$n_peaks, dominant_period=spec_a$dominant_period,
    spectral_confidence=spec_a$conf,
    changepoint_detected=cpt_a$changepoint_detected,
    changepoint_position=cpt_a$changepoint_position,
    n_changepoints=cpt_a$n_changepoints,
    first_changepoint_year=cpt_a$first_changepoint_year,
    mean_before_shift=cpt_a$mean_before, mean_after_shift=cpt_a$mean_after,
    magnitude_shift=cpt_a$magnitude_shift,
    same_significance=(vc_a$p < alpha)==(tf_a$p < alpha),
    same_direction=sign(vc_a$tau)==sign(tf_a$tau),
    is_basin_average=TRUE)
  
  monthly_dts <- lapply(1:12, function(m) {
    mi <- m_index_list[[as.character(m)]]
    mn_count <- if (!skip_min_filter && !is_precip) min_nonzero_pet_monthly else NULL
    res_m <- if (length(mi) && !all(is.na(ts_monthly[mi])))
      mk_tfpw_spectral_for_series(ts_monthly[mi], is_precip, alpha, max_tie_percent,
                                  n_sim_spectral, conf_env,
                                  min_nonzero_count = mn_count,
                                  skip_min_filter   = skip_min_filter)
    else list(vc=list(tau=NA_real_,p=NA_real_,sl=NA_real_,S=NA_real_,varS=NA_real_,
                      n=0,rho1=NA_real_,vc_corrected=FALSE,n_ties=0,percent_ties=0,
                      n_min_vals=0,percent_min_vals=0,tau_b_adjusted=FALSE,
                      filtered=TRUE,reason="no_data"),
              tf=list(tau=NA_real_,p=NA_real_,sl=NA_real_,S=NA_real_,varS=NA_real_,
                      n=0,rho1=NA_real_,tfpw_applied=FALSE,n_ties=0,percent_ties=0,
                      n_min_vals=0,percent_min_vals=0,tau_b_adjusted=FALSE,
                      filtered=TRUE,reason="no_data"),
              spec=list(n_peaks=0L,dominant_period=NA_real_,conf=NA_real_),
              cpt=list(changepoint_detected=FALSE,changepoint_position=NA_integer_,
                       n_changepoints=0L,first_changepoint_year=NA_integer_,
                       mean_before=NA_real_,mean_after=NA_real_,magnitude_shift=NA_real_))
    vc_m <- res_m$vc; tf_m <- res_m$tf; sp_m <- res_m$spec; cp_m <- res_m$cpt
    data.table(
      space_idx=NA_integer_,x=NA_real_,y=NA_real_,
      variable=var_name, period="monthly", month=m, variable_name=var_name,
      tau_vc=vc_m$tau, p_value_vc=vc_m$p, sl_vc=vc_m$sl, vc_corrected=vc_m$vc_corrected,
      n_ties_vc=vc_m$n_ties, percent_ties_vc=vc_m$percent_ties,
      n_min_vals_vc=vc_m$n_min_vals, percent_min_vals_vc=vc_m$percent_min_vals,
      tau_b_adjusted_vc=vc_m$tau_b_adjusted, filtered_vc=vc_m$filtered,
      filter_reason_vc=vc_m$reason,
      tau_tfpw=tf_m$tau, p_value_tfpw=tf_m$p, sl_tfpw=tf_m$sl, tfpw_applied=tf_m$tfpw_applied,
      n_ties_tfpw=tf_m$n_ties, percent_ties_tfpw=tf_m$percent_ties,
      n_min_vals_tfpw=tf_m$n_min_vals, percent_min_vals_tfpw=tf_m$percent_min_vals,
      tau_b_adjusted_tfpw=tf_m$tau_b_adjusted, filtered_tfpw=tf_m$filtered,
      filter_reason_tfpw=tf_m$reason,
      n=NA_integer_, rho1_vc=NA_real_, rho1_tfpw=NA_real_,
      n_spectral_peaks=sp_m$n_peaks, dominant_period=sp_m$dominant_period,
      spectral_confidence=sp_m$conf,
      changepoint_detected=cp_m$changepoint_detected,
      changepoint_position=cp_m$changepoint_position,
      n_changepoints=cp_m$n_changepoints,
      first_changepoint_year=cp_m$first_changepoint_year,
      mean_before_shift=cp_m$mean_before, mean_after_shift=cp_m$mean_after,
      magnitude_shift=cp_m$magnitude_shift,
      same_significance=(vc_m$p < alpha)==(tf_m$p < alpha),
      same_direction=sign(vc_m$tau)==sign(tf_m$tau),
      is_basin_average=TRUE)
  })
  rbindlist(c(list(annual_dt), monthly_dts))
}

log_event("Processing basin-averaged series (Pr, PET_PM, PET_Thw, Temperature)...")
log_event("  [DIAG] Basin-average preflight: reporting variance / tie density before trend testing.")
log_series_quality(precip_monthly_avg,  "Basin Pr monthly", 0)
log_series_quality(pet_monthly_avg,     "Basin PET monthly", min_positive_value)
log_series_quality(pet_thw_monthly_avg, "Basin PET_Thw monthly", min_positive_value)
log_series_quality(tair_monthly_avg,    "Basin Tair monthly", NA_real_)
log_series_quality(precip_annual_avg,   "Basin Pr annual", 0)
log_series_quality(pet_annual_avg,      "Basin PET annual", min_positive_value)
log_series_quality(pet_thw_annual_avg,  "Basin PET_Thw annual", min_positive_value)
log_series_quality(tair_annual_avg,     "Basin Tair annual", NA_real_)
basin_precip_results    <- process_basin_series(precip_monthly_avg,  precip_annual_avg,
                                                "Precipitation", is_precip=TRUE)
basin_pet_results       <- process_basin_series(pet_monthly_avg,     pet_annual_avg,
                                                "PET", is_precip=FALSE)
basin_pet_thw_results   <- process_basin_series(pet_thw_monthly_avg, pet_thw_annual_avg,  # [NEW]
                                                "PET_Thw", is_precip=FALSE)
basin_tair_results      <- process_basin_series(tair_monthly_avg,    tair_annual_avg,
                                                "Temperature", is_precip=FALSE,
                                                skip_min_filter=TRUE,
                                                m_index_list=month_index_list_t)

####################################################################################
# ?????? PIXEL-LEVEL PROCESSING
####################################################################################
process_variable_final <- function(data_matrix, var_name, coords_dt,
                                   is_precip=FALSE, skip_min_filter=FALSE,
                                   yrs=years, mos=months) {
  log_event(paste("Processing", var_name, "pixels..."))
  if (!is_precip && !skip_min_filter) {
    log_event("  Aggregating PET to warm-season annual sum (Apr-Oct)...")
    annual_matrix <- aggregate_to_annual_warmseason_fast(
      data_matrix, yrs, months_vec=mos, warm_months=warm_months_pet, method="sum")
  } else if (skip_min_filter) {
    log_event("  Aggregating Temperature to full-year annual mean...")
    annual_matrix <- aggregate_to_annual_fast(data_matrix, yrs, method="mean")
  } else {
    log_event("  Aggregating Precipitation to annual sum...")
    annual_matrix <- aggregate_to_annual_fast(data_matrix, yrs, method="sum")
  }
  
  conf_env          <- new.env(parent=emptyenv())
  start_year_annual <- min(as.integer(rownames(annual_matrix)))
  
  log_event("  Running VC + TFPW + Spectral + Changepoint (annual, parallel)...")
  matrix_quality_profile(annual_matrix, sprintf("%s annual aggregated", var_name),
                         floor_value = if (is_precip) 0 else min_positive_value)
  res_annual <- future_lapply(
    seq_len(ncol(annual_matrix)),
    function(i) mk_tfpw_spectral_for_series(
      annual_matrix[,i], is_precip, alpha, max_tie_percent,
      n_sim_spectral, conf_env,
      start_year      = start_year_annual,
      skip_min_filter = skip_min_filter),
    future.seed=TRUE)
  
  annual_dt <- build_annual_dt(res_annual, coords_dt, var_name)
  annual_dt[, variable := var_name]
  
  log_event("  Processing 12 calendar months (parallel)...")
  monthly_results_list <- vector("list", 12L)
  for (m in 1:12) {
    mi_key  <- as.character(m)
    mi_use  <- if (skip_min_filter) month_index_list_t[[mi_key]] else
      month_index_list[[mi_key]]
    monthly_sub <- if (length(mi_use)) data_matrix[mi_use, , drop=FALSE] else
      matrix(NA_real_, 0, ncol(data_matrix))
    mn_count  <- if (!skip_min_filter && !is_precip) min_nonzero_pet_monthly else NULL
    matrix_quality_profile(monthly_sub, sprintf("%s month %02d aggregated", var_name, m),
                           floor_value = if (is_precip) 0 else min_positive_value)
    
    res_m <- future_lapply(
      seq_len(ncol(monthly_sub)),
      function(i) mk_tfpw_spectral_for_series(
        monthly_sub[,i], is_precip, alpha, max_tie_percent,
        n_sim_spectral, conf_env,
        min_nonzero_count = mn_count,
        skip_min_filter   = skip_min_filter),
      future.seed=TRUE)
    
    vc_m  <- unpack_results(res_m, "vc")
    tf_m  <- unpack_results(res_m, "tf")
    spec_m <- data.table(
      n_spectral_peaks    = vapply(lapply(res_m,`[[`,"spec"),`[[`,integer(1),"n_peaks"),
      dominant_period     = vapply(lapply(res_m,`[[`,"spec"),`[[`,numeric(1),"dominant_period"),
      spectral_confidence = vapply(lapply(res_m,`[[`,"spec"),`[[`,numeric(1),"conf"))
    cpt_m  <- data.table(
      changepoint_detected   = vapply(lapply(res_m,`[[`,"cpt"),`[[`,logical(1),"changepoint_detected"),
      changepoint_position   = vapply(lapply(res_m,`[[`,"cpt"),`[[`,integer(1),"changepoint_position"),
      n_changepoints         = vapply(lapply(res_m,`[[`,"cpt"),`[[`,integer(1),"n_changepoints"),
      first_changepoint_year = vapply(lapply(res_m,`[[`,"cpt"),`[[`,integer(1),"first_changepoint_year"),
      mean_before_shift      = vapply(lapply(res_m,`[[`,"cpt"),`[[`,numeric(1),"mean_before"),
      mean_after_shift       = vapply(lapply(res_m,`[[`,"cpt"),`[[`,numeric(1),"mean_after"),
      magnitude_shift        = vapply(lapply(res_m,`[[`,"cpt"),`[[`,numeric(1),"magnitude_shift"))
    
    monthly_results_list[[m]] <- cbind(
      coords_dt, variable=var_name, period="monthly", month=m, variable_name=var_name,
      data.table(tau_vc=vc_m$tau, p_value_vc=vc_m$p, sl_vc=vc_m$sl,
                 vc_corrected=vc_m$extra, n_ties_vc=vc_m$n_ties,
                 percent_ties_vc=vc_m$pct_ties, n_min_vals_vc=vc_m$n_min,
                 percent_min_vals_vc=vc_m$pct_min, tau_b_adjusted_vc=vc_m$tau_b_adj,
                 filtered_vc=vc_m$filtered, filter_reason_vc=vc_m$reason),
      data.table(tau_tfpw=tf_m$tau, p_value_tfpw=tf_m$p, sl_tfpw=tf_m$sl,
                 tfpw_applied=tf_m$extra, n_ties_tfpw=tf_m$n_ties,
                 percent_ties_tfpw=tf_m$pct_ties, n_min_vals_tfpw=tf_m$n_min,
                 percent_min_vals_tfpw=tf_m$pct_min, tau_b_adjusted_tfpw=tf_m$tau_b_adj,
                 filtered_tfpw=tf_m$filtered, filter_reason_tfpw=tf_m$reason),
      n=NA_integer_, rho1_vc=NA_real_, rho1_tfpw=NA_real_,
      spec_m, cpt_m,
      same_significance=(vc_m$p < alpha)==(tf_m$p < alpha),
      same_direction=sign(vc_m$tau)==sign(tf_m$tau))
  }
  monthly_results <- rbindlist(monthly_results_list, use.names=TRUE, fill=TRUE)
  rbindlist(list(annual_dt, monthly_results), use.names=TRUE, fill=TRUE)
}

log_event("=== STARTING PRECIPITATION ANALYSIS ===")
precip_results   <- process_variable_final(precip_matrix,  "Precipitation", coords_dt,
                                           is_precip=TRUE)
log_event("=== STARTING PET (Penman-Monteith) ANALYSIS ===")
pet_results      <- process_variable_final(pet_matrix,     "PET",       coords_dt,
                                           is_precip=FALSE)
log_event("=== STARTING PET (Thornthwaite) TREND ANALYSIS ===")   # observed framework
pet_thw_results  <- process_variable_final(pet_thw_matrix, "PET_Thw",   coords_dt,
                                           is_precip=FALSE)
log_event("=== STARTING TEMPERATURE ANALYSIS ===")
tair_results     <- process_variable_final(tair_matrix,    "Temperature", coords_dt,
                                           is_precip=FALSE, skip_min_filter=TRUE,
                                           yrs=years_t, mos=months_t)

####################################################################################
# ?????? SPECIFIC POINT EXTRACTION AND TREND ANALYSIS
####################################################################################
log_event("=== SPECIFIC POINT ANALYSIS (Pr, PET_PM, PET_Thw, Temperature) ===")

pts_wgs84  <- vect(SPECIFIC_PTS, geom=c("x","y"), crs="EPSG:4326")
pts_proj   <- project(pts_wgs84, "EPSG:3005")
n_pts      <- nrow(pts_proj)
cat(sprintf("  Extracting at %d specific point(s)...\n", n_pts))

extract_point_monthly <- function(rast_clipped, pts, dates_vec, var_name, unit_label="") {
  rows <- vector("list", nrow(pts))
  for (i in seq_len(nrow(pts))) {
    pt_v  <- pts[i, ]
    vals  <- tryCatch(as.numeric(terra::extract(rast_clipped, pt_v)[1, -1]),
                      error=function(e) rep(NA_real_, nlyr(rast_clipped)))
    n_use <- min(length(dates_vec), length(vals))
    rows[[i]] <- data.frame(
      point_id  = i,
      lon_wgs84 = crds(pts_wgs84)[i, 1],
      lat_wgs84 = crds(pts_wgs84)[i, 2],
      date      = dates_vec[seq_len(n_use)],
      year      = as.integer(format(dates_vec[seq_len(n_use)], "%Y")),
      month     = as.integer(format(dates_vec[seq_len(n_use)], "%m")),
      value     = vals[seq_len(n_use)],
      variable  = var_name,
      unit      = unit_label,
      stringsAsFactors = FALSE
    )
  }
  rbindlist(rows)
}

pt_pr      <- extract_point_monthly(precip_clipped,  pts_proj, dates,   "Precipitation", "mm/month")
pt_pet     <- extract_point_monthly(pet_clipped,     pts_proj, dates,   "PET",           "mm/month")
pt_pet_thw <- extract_point_monthly(pet_thw_clipped, pts_proj, dates,   "PET_Thw",       "mm/month")  # [NEW]
pt_tair    <- extract_point_monthly(tair_clipped,    pts_proj, dates_t, "Temperature",   "degC")

# Combine and save raw monthly series (all 4 variables)
point_monthly_ts <- rbindlist(list(pt_pr, pt_pet, pt_pet_thw, pt_tair),
                              use.names=TRUE, fill=TRUE)
fwrite(point_monthly_ts, file.path(out_dir, "point_monthly_timeseries.csv"))
log_event(sprintf("??? Saved point_monthly_timeseries.csv  (%d rows)", nrow(point_monthly_ts)))

run_point_trends <- function(pt_df, var_name, is_precip,
                             skip_min_filter = FALSE,
                             warm_season_annual = FALSE) {
  results <- list()
  for (pid in unique(pt_df$point_id)) {
    sub <- pt_df[point_id == pid]
    lon <- sub$lon_wgs84[1]; lat <- sub$lat_wgs84[1]
    conf_env <- new.env(parent=emptyenv())
    
    if (warm_season_annual) {
      ann_mat <- aggregate_to_annual_warmseason_fast(
        matrix(sub$value, ncol=1), sub$year, months_vec=sub$month,
        warm_months=warm_months_pet, method="sum")
    } else if (skip_min_filter) {
      ann_mat <- aggregate_to_annual_fast(matrix(sub$value, ncol=1), sub$year, method="mean")
    } else {
      ann_mat <- aggregate_to_annual_fast(matrix(sub$value, ncol=1), sub$year, method="sum")
    }
    ann_vec  <- ann_mat[, 1]
    ann_yrs  <- as.integer(rownames(ann_mat))
    
    res_ann <- mk_tfpw_spectral_for_series(
      ann_vec, is_precip, alpha, max_tie_percent, n_sim_spectral, conf_env,
      start_year=min(ann_yrs, na.rm=TRUE), skip_min_filter=skip_min_filter)
    vc_a <- res_ann$vc; tf_a <- res_ann$tf; cpt_a <- res_ann$cpt
    
    results[[length(results)+1]] <- data.frame(
      point_id=pid, lon_wgs84=lon, lat_wgs84=lat,
      variable=var_name, period="annual", month=NA_integer_, month_abb="Annual",
      n_obs=vc_a$n,
      year_min=min(ann_yrs, na.rm=TRUE), year_max=max(ann_yrs, na.rm=TRUE),
      tau_vc=vc_a$tau, p_value_vc=vc_a$p, sl_vc=vc_a$sl,
      filtered_vc=vc_a$filtered, filter_reason_vc=vc_a$reason,
      tau_tfpw=tf_a$tau, p_value_tfpw=tf_a$p, sl_tfpw=tf_a$sl,
      filtered_tfpw=tf_a$filtered, filter_reason_tfpw=tf_a$reason,
      rho1=vc_a$rho1,
      changepoint_detected=cpt_a$changepoint_detected,
      first_changepoint_year=cpt_a$first_changepoint_year,
      magnitude_shift=cpt_a$magnitude_shift,
      sig_vc  =!is.na(vc_a$p) && vc_a$p < alpha,
      sig_tfpw=!is.na(tf_a$p) && tf_a$p < alpha,
      stringsAsFactors=FALSE)
    
    for (m in 1:12) {
      sub_m <- sub[month == m]
      if (nrow(sub_m) < 10) next
      mn_count <- if (!skip_min_filter && !is_precip) min_nonzero_pet_monthly else NULL
      res_m <- mk_tfpw_spectral_for_series(
        sub_m$value, is_precip, alpha, max_tie_percent, n_sim_spectral, conf_env,
        min_nonzero_count=mn_count, skip_min_filter=skip_min_filter)
      vc_m <- res_m$vc; tf_m <- res_m$tf; cpt_m <- res_m$cpt
      
      results[[length(results)+1]] <- data.frame(
        point_id=pid, lon_wgs84=lon, lat_wgs84=lat,
        variable=var_name, period="monthly", month=m, month_abb=month.abb[m],
        n_obs=vc_m$n,
        year_min=min(sub_m$year, na.rm=TRUE), year_max=max(sub_m$year, na.rm=TRUE),
        tau_vc=vc_m$tau, p_value_vc=vc_m$p, sl_vc=vc_m$sl,
        filtered_vc=vc_m$filtered, filter_reason_vc=vc_m$reason,
        tau_tfpw=tf_m$tau, p_value_tfpw=tf_m$p, sl_tfpw=tf_m$sl,
        filtered_tfpw=tf_m$filtered, filter_reason_tfpw=tf_m$reason,
        rho1=vc_m$rho1,
        changepoint_detected=cpt_m$changepoint_detected,
        first_changepoint_year=cpt_m$first_changepoint_year,
        magnitude_shift=cpt_m$magnitude_shift,
        sig_vc  =!is.na(vc_m$p) && vc_m$p < alpha,
        sig_tfpw=!is.na(tf_m$p) && tf_m$p < alpha,
        stringsAsFactors=FALSE)
    }
  }
  rbindlist(results, fill=TRUE)
}

pt_pr_trends      <- run_point_trends(pt_pr,      "Precipitation", is_precip=TRUE)
pt_pet_trends     <- run_point_trends(pt_pet,     "PET",           is_precip=FALSE,
                                      warm_season_annual=TRUE)
pt_pet_thw_trends <- run_point_trends(pt_pet_thw, "PET_Thw",       is_precip=FALSE,   # [NEW]
                                      warm_season_annual=TRUE)
pt_tair_trends    <- run_point_trends(pt_tair,    "Temperature",   is_precip=FALSE,
                                      skip_min_filter=TRUE)

point_trend_stats <- rbindlist(
  list(pt_pr_trends, pt_pet_trends, pt_pet_thw_trends, pt_tair_trends),
  use.names=TRUE, fill=TRUE)
fwrite(point_trend_stats, file.path(out_dir, "point_trend_stats.csv"))
log_event(sprintf("??? Saved point_trend_stats.csv  (%d rows ?? %d cols)",
                  nrow(point_trend_stats), ncol(point_trend_stats)))

cat("\n?????? Specific Point Trend Summary ??????\n")
smry_pts <- point_trend_stats[, .(
  n_periods  = .N,
  sig_vc_n   = sum(sig_vc,   na.rm=TRUE),
  sig_tfpw_n = sum(sig_tfpw, na.rm=TRUE)
), by=.(point_id, variable)]
print(as.data.frame(smry_pts))

####################################################################################
# ?????? COMBINE ALL RESULTS & SAVE
####################################################################################
log_event("Combining all pixel-level and basin results (4 variables)...")
all_results <- rbindlist(
  list(precip_results,     pet_results,     pet_thw_results,  tair_results,
       basin_precip_results, basin_pet_results, basin_pet_thw_results, basin_tair_results),
  use.names=TRUE, fill=TRUE)

saveRDS(all_results, file.path(out_dir, "all_results.rds"), compress="gzip")

saveRDS(list(
  precip_results        = precip_results,
  pet_results           = pet_results,
  pet_thw_results       = pet_thw_results,       # [NEW]
  tair_results          = tair_results,
  basin_precip_results  = basin_precip_results,
  basin_pet_results     = basin_pet_results,
  basin_pet_thw_results = basin_pet_thw_results,  # [NEW]
  basin_tair_results    = basin_tair_results,
  basin_pixels          = n_basin_pixels,
  bbox_cells            = n_bbox_cells,
  reduction_pct         = reduction_pct,
  processing_date       = Sys.time(),
  parameters = list(alpha=alpha, max_tie_percent=max_tie_percent,
                    max_min_value_pct_precip=max_min_value_pct_precip,
                    max_min_value_pct_pet=max_min_value_pct_pet),
  clipped_extent = list(xmin=xmin(clipped_template), xmax=xmax(clipped_template),
                        ymin=ymin(clipped_template), ymax=ymax(clipped_template),
                        nrows=nrow(clipped_template), ncols=ncol(clipped_template),
                        res=res(clipped_template)),
  basin_avg_monthly = data.frame(
    date               = dates,
    precip_mm_month    = precip_monthly_avg,
    pet_mm_month       = pet_monthly_avg,
    pet_thw_mm_month   = pet_thw_monthly_avg,   # [NEW]
    tair_degC_month    = tair_monthly_avg[seq_len(length(dates))]
  ),
  basin_avg_annual = data.frame(
    year               = annual_years,
    precip_mm_year     = precip_annual_avg,
    pet_mm_year        = pet_annual_avg,
    pet_thw_mm_year    = pet_thw_annual_avg,     # [NEW]
    tair_degC_year     = tair_annual_avg[seq_len(length(annual_years))]
  ),
  specific_pts = SPECIFIC_PTS
), file.path(out_dir, "analysis_metadata.rds"))

summary_stats <- all_results[, .(
  n_total            = .N,
  n_valid_vc         = sum(!filtered_vc & !is.na(p_value_vc)),
  n_significant_vc   = sum(!filtered_vc & p_value_vc < alpha, na.rm=TRUE),
  n_valid_tfpw       = sum(!filtered_tfpw & !is.na(p_value_tfpw)),
  n_significant_tfpw = sum(!filtered_tfpw & p_value_tfpw < alpha, na.rm=TRUE)
), by=.(variable, period, month, is_basin_average)]
fwrite(summary_stats, file.path(out_dir, "summary_statistics.csv"))

log_event("????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????")
log_event("DIAGNOSTIC SUMMARY: root-cause screen for MK / TFPW warnings.")
diag_summary <- log_result_root_cause_summary(all_results)
if (!is.null(diag_summary)) {
  print(diag_summary)
}
log_event("????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????")
log_event("DATA PROCESSING COMPLETE (Pr + PET_PM + PET_Thw + Temperature)")
cat("\n??? DATA PROCESSING COMPLETED SUCCESSFULLY\n")
cat(sprintf("  Basin pixels: %d  |  Specific points: %d\n", n_basin_pixels, n_pts))
cat(sprintf("  Variables: Precipitation, PET (PM), PET (Thornthwaite), Temperature\n"))
cat(sprintf("  Outputs: %s\n", normalizePath(out_dir)))
cat(sprintf("  all_results.rds              ??? 4-variable pixel + basin trend stats\n"))
cat(sprintf("  point_trend_stats.csv        ??? annual + 12-month trends at %d points\n", n_pts))
cat(sprintf("  point_monthly_timeseries.csv ??? raw monthly Pr / PET_PM / PET_Thw / Tair\n"))

future::plan(future::sequential)