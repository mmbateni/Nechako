####################################################################################
# 4pr_pet_trends.r  ·  Trend Analysis — Precipitation, PET (both methods), Temperature
# ─────────────────────────────────────────────────────────────────────────────────
# Performs comprehensive VC-MK + TFPW-MK + Spectral + Changepoint trend analysis
# for ALL BASIN PIXELS and SPECIFIC POINTS for five variables:
#   • Precipitation   (ERA5-Land, mm/month → mm/year, full-year sum)
#   • PET             (Penman-Monteith, mm/month → warm-season sum Apr–Oct)
#   • PET_Thw         (Thornthwaite, mm/month → warm-season sum Apr–Oct)
#   • Temperature     (ERA5-Land 2 m Tair, K→°C, annual mean + 12 calendar months)
#
# MODIFICATIONS vs. original:
#   [NEW] Thornthwaite PET (PET_Thw) added as a parallel variable to PM PET.
#         Both PET variants go through identical preprocessing, aggregation,
#         pixel-level trend analysis, basin-average analysis, specific-point
#         extraction and trend analysis.
#
# OUTPUTS (binary .rds + CSV):
#   all_results.rds              pixel + basin trend statistics (Pr, PET, PET_Thw, Tair)
#   analysis_metadata.rds        time series, extents, parameters (incl. PET_Thw)
#   clipped_template.rds         EPSG:3005 raster template for visualization
#   summary_statistics.csv       per-variable significance counts
#   point_trend_stats.csv        per-point stats (annual + monthly × 4 vars)
#   point_monthly_timeseries.csv raw monthly Pr / PET / PET_Thw / Tair at each point
#   data_processing.log
####################################################################################

rm(list = ls())
gc()

library(ncdf4)
library(terra)
library(data.table)
library(Kendall)
library(future.apply)
library(zoo)
library(sf)
library(changepoint)
library(trend)
library(modifiedmk)

# ================= USER / ENV =================
setwd("D:/Nechako_Drought/Nechako/")

out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ===== LOGGING SETUP =====
LOG_FILE <- file.path(out_dir, "data_processing.log")
log_event <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0(timestamp, "   ", msg, "\n"), file = LOG_FILE, append = TRUE)
  message(paste0(timestamp, "   ", msg))
}
cat("Trend Analysis Data Processing - Started\n", file = LOG_FILE)
cat(paste("Timestamp: ", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat("==========================================\n", file = LOG_FILE, append = TRUE)

####################################################################################
# ── SPECIFIC POINTS CONFIG
####################################################################################
SPECIFIC_PTS <- data.frame(
  x = c(-124.5, -123.0, -125.0),
  y = c(  54.0,   53.8,   54.2)
)

####################################################################################
# ── PARAMETERS
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
# [NEW] Thornthwaite PET — produced by 2preq_PET_ERALand.R (Section 5b)
pet_thw_path  <- "monthly_data_direct/potential_evapotranspiration_thornthwaite_monthly.nc"

TAIR_NC_CANDIDATES <- c(
  "monthly_data_direct/2m_temperature_monthly.nc",
  "monthly_data_direct/t2m_monthly.nc",
  "monthly_data_direct/air_temperature_2m_monthly.nc",
  "monthly_data_direct/ERA5Land_T2m_monthly.nc"
)

num_cores <- min(parallel::detectCores() - 1, 8)
if (num_cores < 1) num_cores <- 1
plan(multisession, workers = num_cores)
log_event(paste("Using", num_cores, "cores for parallel processing"))

####################################################################################
# ── BASIN BOUNDARY
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
      log_event(paste("✓ Loaded basin boundary from:", bf))
      break
    }, error = function(e) log_event(paste("  Error loading", bf, ":", e$message)))
  }
}
if (is.null(basin_boundary))
  stop("CRITICAL: Nechako Basin boundary NOT FOUND.")
saveRDS(basin_boundary, file.path(out_dir, "basin_boundary.rds"))

####################################################################################
# ── HELPER FUNCTIONS
####################################################################################
check_min_value_threshold <- function(ts_clean, min_val_threshold = 0.01, max_pct = 50,
                                      is_precip = FALSE) {
  n_min_vals   <- sum(ts_clean <= min_val_threshold, na.rm = TRUE)
  pct_min_vals <- n_min_vals / length(ts_clean) * 100
  max_allowed  <- if (is_precip) max_min_value_pct_precip else max_min_value_pct_pet
  list(exceeds_threshold = pct_min_vals > max_allowed, pct_min_vals = pct_min_vals)
}

calculate_variance_with_ties <- function(S, n, x) {
  tie_table  <- table(x)
  tie_counts <- tie_table[tie_table > 1]
  var_s      <- n * (n - 1) * (2 * n + 5) / 18
  if (length(tie_counts))
    var_s <- var_s - sum(tie_counts * (tie_counts - 1) * (2 * tie_counts + 5)) / 18
  var_s
}

detect_regime_shift_pelt <- function(ts_vec, min_obs = 20, start_year = NULL) {
  ts_clean <- na.omit(ts_vec)
  n <- length(ts_clean)
  empty <- list(changepoint_detected = FALSE, changepoint_position = NA_integer_,
                n_changepoints = 0L, first_changepoint_year = NA_integer_,
                mean_before = NA_real_, mean_after = NA_real_, magnitude_shift = NA_real_)
  if (n < min_obs) return(empty)
  tryCatch({
    cpt      <- cpt.mean(ts_clean, method = "PELT", penalty = "BIC")
    cpts_det <- cpts(cpt)
    if (!length(cpts_det)) return(empty)
    fc <- cpts_det[1]
    list(changepoint_detected   = TRUE,
         changepoint_position   = fc,
         n_changepoints         = length(cpts_det),
         first_changepoint_year = as.integer(if (!is.null(start_year)) start_year + fc - 1 else fc),
         mean_before            = mean(ts_clean[1:fc], na.rm = TRUE),
         mean_after             = mean(ts_clean[(fc+1):n], na.rm = TRUE),
         magnitude_shift        = mean(ts_clean[(fc+1):n], na.rm = TRUE) -
           mean(ts_clean[1:fc], na.rm = TRUE))
  }, error = function(e) empty)
}

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
# ── CORE TREND ANALYSIS FUNCTION
####################################################################################
mk_tfpw_spectral_for_series <- function(ts_vec, is_precip, alpha, max_tie_pct,
                                        n_sim_spectral, conf_cache_env,
                                        start_year       = NULL,
                                        min_nonzero_count = NULL,
                                        skip_min_filter  = FALSE) {
  ts_clean <- na.omit(ts_vec)
  n        <- length(ts_clean)
  
  na_result <- function(reason, n_ties = 0, pct_ties = 0, pct_min = 0) {
    list(
      vc   = list(tau=NA_real_, p=NA_real_, sl=NA_real_, S=NA_real_, varS=NA_real_,
                  n=n, rho1=NA_real_, vc_corrected=FALSE, n_ties=n_ties,
                  percent_ties=pct_ties, n_min_vals=0, percent_min_vals=pct_min,
                  tau_b_adjusted=FALSE, filtered=TRUE, reason=reason),
      tf   = list(tau=NA_real_, p=NA_real_, sl=NA_real_, S=NA_real_, varS=NA_real_,
                  n=n, rho1=NA_real_, tfpw_applied=FALSE, n_ties=n_ties,
                  percent_ties=pct_ties, n_min_vals=0, percent_min_vals=pct_min,
                  tau_b_adjusted=FALSE, filtered=TRUE, reason=reason),
      spec = list(n_peaks=0L, dominant_period=NA_real_, conf=NA_real_),
      cpt  = list(changepoint_detected=FALSE, changepoint_position=NA_integer_,
                  n_changepoints=0L, first_changepoint_year=NA_integer_,
                  mean_before=NA_real_, mean_after=NA_real_, magnitude_shift=NA_real_)
    )
  }
  
  if (n < 10) return(na_result("low_n"))
  
  if (skip_min_filter) {
    min_check <- list(exceeds_threshold = FALSE, pct_min_vals = 0)
  } else if (!is.null(min_nonzero_count)) {
    n_nonzero    <- sum(ts_clean > min_positive_value, na.rm = TRUE)
    pct_min_vals <- (1 - n_nonzero / n) * 100
    if (n_nonzero < min_nonzero_count)
      return(na_result("insufficient_nonzero", pct_min = pct_min_vals))
    min_check <- list(exceeds_threshold = FALSE, pct_min_vals = pct_min_vals)
  } else {
    min_val_threshold <- if (is_precip) 0.0 else min_positive_value
    min_check <- check_min_value_threshold(
      ts_clean, min_val_threshold,
      max_pct = if (is_precip) max_min_value_pct_precip else max_min_value_pct_pet,
      is_precip = is_precip)
    if (min_check$exceeds_threshold)
      return(na_result("excessive_min_vals", pct_min = min_check$pct_min_vals))
  }
  
  n_unique     <- length(unique(ts_clean))
  n_ties       <- n - n_unique
  percent_ties <- n_ties / n * 100
  if (percent_ties > max_tie_pct)
    return(na_result("excessive_ties", n_ties, percent_ties, min_check$pct_min_vals))
  
  sen_slope <- trend::sens.slope(ts_clean)$estimates
  if (is.na(sen_slope) || is.infinite(sen_slope) || length(sen_slope) == 0)
    return(na_result("sens_slope_fail", n_ties, percent_ties, min_check$pct_min_vals))
  
  time_index <- seq_len(n)
  trend_line <- sen_slope * time_index
  detrended  <- ts_clean - trend_line
  
  acf_result <- tryCatch(acf(detrended, lag.max=1, plot=FALSE, na.action=na.pass),
                         error=function(e) NULL, warning=function(w) NULL)
  rho1 <- if (!is.null(acf_result)) acf_result$acf[2] else NA_real_
  
  tau_b_adjusted <- (percent_ties > 5)
  vc_corrected   <- FALSE
  p_vc <- NA_real_; S_vc <- NA_real_; tau_vc <- NA_real_; varS_final <- NA_real_
  vc_res <- tryCatch(modifiedmk::mk3(ts_clean), error=function(e) NULL)
  if (!is.null(vc_res)) {
    tau_vc     <- vc_res$tau;  p_vc   <- vc_res$p.value
    S_vc       <- vc_res$S;    varS_final <- vc_res$varS
    if (!is.na(rho1) && abs(rho1) > 0.1) vc_corrected <- TRUE
  }
  vc_list <- list(tau=tau_vc, p=p_vc, sl=sen_slope, S=S_vc, varS=varS_final,
                  n=n, rho1=rho1, vc_corrected=vc_corrected,
                  n_ties=n_ties, percent_ties=percent_ties,
                  n_min_vals=0, percent_min_vals=min_check$pct_min_vals,
                  tau_b_adjusted=tau_b_adjusted, filtered=FALSE, reason="none")
  
  tfpw_applied <- FALSE
  if (!is.na(rho1) && abs(rho1) > 0.1) {
    pw <- numeric(n); pw[1] <- detrended[1]
    for (j in 2:n) pw[j] <- detrended[j] - rho1 * detrended[j-1]
    corrected <- pw + trend_line; tfpw_applied <- TRUE
  } else {
    corrected <- ts_clean
  }
  s_mat_tf <- sign(outer(corrected, corrected, `-`))
  S_tf     <- sum(s_mat_tf[upper.tri(s_mat_tf)], na.rm=TRUE)
  varS_tf  <- calculate_variance_with_ties(S_tf, n, corrected)
  n_pairs  <- n * (n-1) / 2
  tau_tf   <- S_tf / n_pairs
  p_tf     <- if (varS_tf <= 0) NA_real_ else 2 * pnorm(-abs(S_tf / sqrt(varS_tf)))
  tf_list  <- list(tau=tau_tf, p=p_tf, sl=sen_slope, S=S_tf, varS=varS_tf,
                   n=n, rho1=rho1, tfpw_applied=tfpw_applied,
                   n_ties=n_ties, percent_ties=percent_ties,
                   n_min_vals=0, percent_min_vals=min_check$pct_min_vals,
                   tau_b_adjusted=tau_b_adjusted, filtered=FALSE, reason="none")
  
  dstd <- detrended - mean(detrended)
  sd_d <- stats::sd(dstd); if (!is.finite(sd_d) || sd_d==0) sd_d <- 1
  dstd <- dstd / sd_d
  key  <- paste0("n_", n)
  if (!exists(key, envir=conf_cache_env, inherits=FALSE)) {
    half        <- floor(n/2)
    max_spectra <- vapply(seq_len(n_sim_spectral), function(jj) {
      r <- rnorm(n,0,1); fr <- fft(r - mean(r))
      max(Mod(fr[1:half])^2 / n, na.rm=TRUE)
    }, numeric(1))
    assign(key, stats::quantile(max_spectra, 1-alpha, na.rm=TRUE), envir=conf_cache_env)
  }
  conf_limit       <- get(key, envir=conf_cache_env, inherits=FALSE)
  half             <- floor(n/2)
  ff               <- fft(dstd)
  spectral_density <- Mod(ff[1:half])^2 / n
  freqs            <- seq(0, 0.5, length.out=half)
  peak_idx         <- which(spectral_density > conf_limit)
  peak_periods     <- if (length(peak_idx)) {
    pf  <- freqs[peak_idx]
    pp  <- ifelse(pf > 0, 1/pf, NA_real_)
    pp[order(spectral_density[peak_idx], decreasing=TRUE)]
  } else numeric(0)
  spec_list <- list(n_peaks=length(peak_periods),
                    dominant_period=if (length(peak_periods)) peak_periods[1] else NA_real_,
                    conf=conf_limit)
  
  cpt_result <- detect_regime_shift_pelt(ts_clean, min_obs=min_obs_changepoint,
                                         start_year=start_year)
  list(vc=vc_list, tf=tf_list, spec=spec_list, cpt=cpt_result)
}

####################################################################################
# ── DATA LOADING & PREPROCESSING
####################################################################################
log_event("Loading Precipitation, PET (PM + Thornthwaite), and Temperature data...")

precip_full <- rast(precip_path)
pet_full    <- rast(pet_path)

# [NEW] Thornthwaite PET
if (!file.exists(pet_thw_path))
  stop("Thornthwaite PET NetCDF not found: ", pet_thw_path,
       "\nRun 2preq_PET_ERALand.R first.")
pet_thw_full <- rast(pet_thw_path)
log_event(sprintf("  [PET_Thw] Loaded: %d layers from %s", nlyr(pet_thw_full), basename(pet_thw_path)))

# ── Temperature ─────────────────────────────────────────────────────────────────
tair_nc <- NULL
for (cand in TAIR_NC_CANDIDATES) {
  if (file.exists(cand)) { tair_nc <- cand; break }
}
if (is.null(tair_nc))
  stop("ERA5-Land 2 m temperature NetCDF not found. Tried:\n  ",
       paste(TAIR_NC_CANDIDATES, collapse="\n  "))
log_event(paste("  [temp] Using:", basename(tair_nc)))
tair_full <- rast(tair_nc)

# ── Date extraction ──────────────────────────────────────────────────────────────
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

# ── Unit conversions ─────────────────────────────────────────────────────────────
precip_full <- precip_full * 1000    # m → mm/day
month_days_base  <- c(31,28,31,30,31,30,31,31,30,31,30,31)
is_leap <- function(yr) (yr%%4==0) & (yr%%100!=0 | yr%%400==0)
days_in_month <- month_days_base[months]
days_in_month[is_leap(years) & months==2] <- 29

log_event("Converting Pr, PET (PM), and PET (Thw) from mm/day to mm/month...")
for (i in seq_len(nlyr(precip_full)))   precip_full[[i]]   <- precip_full[[i]]   * days_in_month[i]
for (i in seq_len(nlyr(pet_full)))      pet_full[[i]]      <- pet_full[[i]]      * days_in_month[i]
for (i in seq_len(nlyr(pet_thw_full)))  pet_thw_full[[i]]  <- pet_thw_full[[i]]  * days_in_month[i]

# Temperature: Kelvin → °C
tair_sample <- as.numeric(terra::global(tair_full[[1]], "mean", na.rm=TRUE))
if (!is.na(tair_sample) && tair_sample > 200) {
  log_event("  [temp] Converting K → °C...")
  tair_full <- tair_full - 273.15
}

n_tair  <- nlyr(tair_full)
dates_t <- tryCatch({
  d <- as.Date(terra::time(tair_full))
  if (is.null(d) || length(d) != n_tair || all(is.na(d))) NULL else d
}, error=function(e) NULL)
if (is.null(dates_t)) {
  dates_t  <- seq(as.Date(sprintf("%d-01-01", TEMP_START_YR)), by="month", length.out=n_tair)
  log_event("  [temp] terra::time() invalid — using positional date sequence")
}
years_t  <- as.integer(format(dates_t, "%Y"))
months_t <- as.integer(format(dates_t, "%m"))
log_event(sprintf("  [temp] Date range: %s to %s", min(dates_t), max(dates_t)))

####################################################################################
# ── REPROJECT & CLIP TO BASIN
####################################################################################
log_event("Reprojecting to BC Albers (EPSG:3005)...")
precip_full   <- project(precip_full,   "EPSG:3005", method="bilinear")
pet_full      <- project(pet_full,      "EPSG:3005", method="bilinear")
pet_thw_full  <- project(pet_thw_full,  "EPSG:3005", method="bilinear")  # [NEW]
tair_full     <- project(tair_full,     "EPSG:3005", method="bilinear")

log_event("Clipping to Nechako Basin extent...")
basin_extent     <- ext(basin_boundary)
precip_clipped   <- mask(crop(precip_full,   basin_extent), vect(basin_boundary), touches=TRUE)
pet_clipped      <- mask(crop(pet_full,      basin_extent), vect(basin_boundary), touches=TRUE)
pet_thw_clipped  <- mask(crop(pet_thw_full,  basin_extent), vect(basin_boundary), touches=TRUE)  # [NEW]
tair_clipped     <- mask(crop(tair_full,     basin_extent), vect(basin_boundary), touches=TRUE)

clipped_template <- rast(precip_clipped, nlyrs=1); values(clipped_template) <- NA_real_
saveRDS(clipped_template, file.path(out_dir, "clipped_template.rds"))
log_event(sprintf("✓ Saved clipped template (EPSG:3005, %d x %d, res: %.0f m)",
                  nrow(clipped_template), ncol(clipped_template), res(clipped_template)[1]))

n_bbox_cells     <- ncell(precip_clipped)
first_layer_vals <- values(precip_clipped[[1]], mat=FALSE)
valid_mask       <- !is.na(first_layer_vals)
n_basin_pixels   <- sum(valid_mask)
reduction_pct    <- 100 * (1 - n_basin_pixels / n_bbox_cells)
log_event(sprintf("✓ BASIN CLIPPING: %d bbox → %d basin pixels (%.1f%% reduction)",
                  n_bbox_cells, n_basin_pixels, reduction_pct))
if (n_basin_pixels == 0) stop("CRITICAL: No valid cells after basin clipping.")

####################################################################################
# ── BASIN AVERAGE TIME SERIES
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

# Pre-process PET (both): replace zeros with min_positive_value
for (clip_r in list(pet_clipped, pet_thw_clipped)) {
  pv  <- values(clip_r)
  nz  <- sum(pv == 0, na.rm=TRUE)
  pv[pv == 0] <- min_positive_value
  values(clip_r) <- pv
  log_event(sprintf("  Replaced %d PET zeros with %.3f mm (%s)",
                    nz, min_positive_value, names(clip_r[[1]])))
}

# Pixel matrices (time × pixels)
valid_cell_indices <- which(valid_mask)
valid_xy           <- xyFromCell(precip_clipped, valid_cell_indices)
precip_matrix      <- t(values(precip_clipped)[valid_mask,    , drop=FALSE])
pet_matrix         <- t(values(pet_clipped)[valid_mask,       , drop=FALSE])
pet_thw_matrix     <- t(values(pet_thw_clipped)[valid_mask,   , drop=FALSE])  # [NEW]
tair_matrix        <- t(values(tair_clipped)[valid_mask,      , drop=FALSE])

coords_dt <- data.table(
  space_idx=seq_len(n_basin_pixels), x=valid_xy[,1], y=valid_xy[,2],
  is_basin_average=FALSE)
log_event(sprintf("  Matrices: Pr/PET %d×%d | PET_Thw %d×%d | Tair %d×%d",
                  nrow(precip_matrix), ncol(precip_matrix),
                  nrow(pet_thw_matrix), ncol(pet_thw_matrix),
                  nrow(tair_matrix), ncol(tair_matrix)))

month_index_list   <- compute_month_index(months)
month_index_list_t <- compute_month_index(months_t)

####################################################################################
# ── HELPER: BUILD RESULT DATA.TABLE ROW SET (annual + 12 monthly)
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
# ── BASIN AVERAGE TREND PROCESSING
####################################################################################
process_basin_series <- function(ts_monthly, ts_annual, var_name,
                                 is_precip, skip_min_filter = FALSE,
                                 m_index_list = month_index_list) {
  conf_env <- new.env(parent=emptyenv())
  res_ann  <- mk_tfpw_spectral_for_series(ts_annual, is_precip, alpha, max_tie_percent,
                                          n_sim_spectral, conf_env,
                                          start_year      = min(annual_years),
                                          skip_min_filter = skip_min_filter)
  vc_a  <- res_ann$vc;  tf_a  <- res_ann$tf
  spec_a <- res_ann$spec; cpt_a <- res_ann$cpt
  
  annual_dt <- data.table(
    space_idx=NA_integer_, x=NA_real_, y=NA_real_,
    variable=var_name, period="annual", month=NA_integer_, variable_name=var_name,
    tau_vc=vc_a$tau, p_value_vc=vc_a$p, sl_vc=vc_a$sl, vc_corrected=FALSE,
    n_ties_vc=vc_a$n_ties, percent_ties_vc=vc_a$percent_ties,
    n_min_vals_vc=vc_a$n_min_vals, percent_min_vals_vc=vc_a$percent_min_vals,
    tau_b_adjusted_vc=vc_a$tau_b_adjusted, filtered_vc=vc_a$filtered,
    filter_reason_vc=vc_a$reason,
    tau_tfpw=tf_a$tau, p_value_tfpw=tf_a$p, sl_tfpw=tf_a$sl, tfpw_applied=FALSE,
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
      tau_vc=vc_m$tau, p_value_vc=vc_m$p, sl_vc=vc_m$sl, vc_corrected=FALSE,
      n_ties_vc=vc_m$n_ties, percent_ties_vc=vc_m$percent_ties,
      n_min_vals_vc=vc_m$n_min_vals, percent_min_vals_vc=vc_m$percent_min_vals,
      tau_b_adjusted_vc=vc_m$tau_b_adjusted, filtered_vc=vc_m$filtered,
      filter_reason_vc=vc_m$reason,
      tau_tfpw=tf_m$tau, p_value_tfpw=tf_m$p, sl_tfpw=tf_m$sl, tfpw_applied=FALSE,
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
# ── PIXEL-LEVEL PROCESSING
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
                 vc_corrected=FALSE, n_ties_vc=vc_m$n_ties,
                 percent_ties_vc=vc_m$pct_ties, n_min_vals_vc=vc_m$n_min,
                 percent_min_vals_vc=vc_m$pct_min, tau_b_adjusted_vc=vc_m$tau_b_adj,
                 filtered_vc=vc_m$filtered, filter_reason_vc=vc_m$reason),
      data.table(tau_tfpw=tf_m$tau, p_value_tfpw=tf_m$p, sl_tfpw=tf_m$sl,
                 tfpw_applied=FALSE, n_ties_tfpw=tf_m$n_ties,
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
log_event("=== STARTING PET (Thornthwaite) ANALYSIS ===")   # [NEW]
pet_thw_results  <- process_variable_final(pet_thw_matrix, "PET_Thw",   coords_dt,
                                           is_precip=FALSE)
log_event("=== STARTING TEMPERATURE ANALYSIS ===")
tair_results     <- process_variable_final(tair_matrix,    "Temperature", coords_dt,
                                           is_precip=FALSE, skip_min_filter=TRUE,
                                           yrs=years_t, mos=months_t)

####################################################################################
# ── SPECIFIC POINT EXTRACTION AND TREND ANALYSIS
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
log_event(sprintf("✓ Saved point_monthly_timeseries.csv  (%d rows)", nrow(point_monthly_ts)))

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
log_event(sprintf("✓ Saved point_trend_stats.csv  (%d rows × %d cols)",
                  nrow(point_trend_stats), ncol(point_trend_stats)))

cat("\n── Specific Point Trend Summary ──\n")
smry_pts <- point_trend_stats[, .(
  n_periods  = .N,
  sig_vc_n   = sum(sig_vc,   na.rm=TRUE),
  sig_tfpw_n = sum(sig_tfpw, na.rm=TRUE)
), by=.(point_id, variable)]
print(as.data.frame(smry_pts))

####################################################################################
# ── COMBINE ALL RESULTS & SAVE
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

log_event("DATA PROCESSING COMPLETE (Pr + PET_PM + PET_Thw + Temperature)")
cat("\n✓ DATA PROCESSING COMPLETED SUCCESSFULLY\n")
cat(sprintf("  Basin pixels: %d  |  Specific points: %d\n", n_basin_pixels, n_pts))
cat(sprintf("  Variables: Precipitation, PET (PM), PET (Thornthwaite), Temperature\n"))
cat(sprintf("  Outputs: %s\n", normalizePath(out_dir)))
cat(sprintf("  all_results.rds              → 4-variable pixel + basin trend stats\n"))
cat(sprintf("  point_trend_stats.csv        → annual + 12-month trends at %d points\n", n_pts))
cat(sprintf("  point_monthly_timeseries.csv → raw monthly Pr / PET_PM / PET_Thw / Tair\n"))

plan(sequential)