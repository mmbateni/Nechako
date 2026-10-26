# ==============================================================================
#   6trend_test_ALL.R
#   SPATIAL TREND ANALYSIS FOR ALL DROUGHT INDICES
# ==============================================================================
#
# SPATIAL GRID CONTRACT
# ---------------------
# Script 6 uses each index's own native grid from its seasonal NetCDF files.
# Inputs include SPI, SPEI-PM, SPEI-ERA5, SWEI-3, and SSPI-1.
#
# There is NO canonical-grid RDS dependency.
#
# For every index × scale:
#   1. The first layer of the actual NetCDF stack defines r_tmpl.
#   2. All layers/files in that index × scale must match r_tmpl geometry.
#   3. The authoritative Nechako basin polygon is transformed to r_tmpl CRS.
#   4. Basin coverage is rasterized directly on r_tmpl.
#   5. Fractional cell coverage × native cell area provides area weights.
#
# No:
#   - era5land_common_analysis_grid.rds
#   - load_common_analysis_grid()
#   - validate_common_analysis_grid()
#   - canonical basin mask
#   - clipped_template.rds reconciliation
#   - cross-pipeline positional cell matching
#   - raster resampling/reprojection of index products
#
# This preserves the native ERA5-Land geometry throughout Script 6.
#
# ==============================================================================
#   STATISTICAL CONTENT
# ==============================================================================
# Computes per-pixel statistics for:
#   SPI
#   SPEI_PM
#   SPEI_ERA5
#   SWEI (scale 3 only)
#   SSPI (scale 1 only)
#
# Mann-Kendall:
#   - Hamed & Rao variance-corrected MK via modifiedmk::mmkh()
#   - Sen's slope
#   - Yue et al. (2002) TFPW MK via the shared project utility (.mk_tfpw_single)
#     — the same production implementation used by Script 4.
#   - BH-FDR correction across valid basin pixels
#
# Drought events:
#   Method 1 — Sheffield & Wood (2008)
#       onset      < DROUGHT_ONSET
#       termination >= DROUGHT_ONSET
#       minimum duration = 1 month
#       intensity = mean(onset threshold - index)
#
#   Method 2 — Hysteresis
#       onset      < DROUGHT_ONSET
#       termination >= DROUGHT_END
#       scale-specific minimum duration
#       intensity = mean(DROUGHT_END - index)
#
# Additional:
#   - PELT regime-shift detection
#   - Wald-Wolfowitz runs test on hysteresis-defined drought events
#   - Monte Carlo AR(1) red-noise spectral peak detection
#
# Outputs:
#   {TREND_DIR}/{index}_{scale:02d}_results.csv
#   duration-class CSVs
#   D13+ NetCDF rasters
#   regional basin-average monthly CSVs
#   basin extent time series
#   script6_run_status.csv
#
# Run BEFORE 4pr_pet_trends_visualization.r
# ==============================================================================


# ==============================================================================
#   0. ENVIRONMENT
# ==============================================================================

WD_PATH <- getwd()

UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")

if (!file.exists(UTILS_FILE)) {
  stop(
    "No drought-analysis utility file found in: ",
    WD_PATH
  )
}

source(UTILS_FILE)

utils_load_packages(
  c(
    "terra",
    "data.table",
    "Kendall",
    "trend",
    "modifiedmk",
    "lubridate",
    "changepoint",
    "zoo"
  )
)

if (!dir.exists(WD_PATH)) {
  stop("Working directory not found: ", WD_PATH)
}

setwd(WD_PATH)

dir.create(
  TREND_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==============================================================================
#   1. NATIVE ERA5-LAND GRID CONTRACT
# ==============================================================================

cat(
  "\n============================================================\n",
  "SCRIPT 6 — NATIVE ERA5-LAND GRID MODE\n",
  "============================================================\n",
  "No canonical-grid RDS is used.\n",
  "Each index × scale uses the native grid of its own NetCDF files.\n",
  "Basin masks and fractional coverage are generated directly on that grid.\n",
  "============================================================\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# Native-grid geometry validator
# ------------------------------------------------------------------------------
#
# Every NetCDF participating in one index × scale must have the same geometry
# as the first NetCDF layer used for that run.
#
# This is NOT a canonical-grid comparison. It is only an internal consistency
# check within the actual native ERA5-Land product being processed.
# ------------------------------------------------------------------------------

assert_same_native_grid <- function(x, template, label = "Raster") {
  
  if (!inherits(x, "SpatRaster")) {
    stop(label, " is not a SpatRaster.")
  }
  
  if (!inherits(template, "SpatRaster")) {
    stop("Native-grid template is not a SpatRaster.")
  }
  
  geom_ok <- terra::compareGeom(
    x[[1]],
    template[[1]],
    stopOnError = FALSE,
    crs      = TRUE,
    ext      = TRUE,
    rowcol   = TRUE,
    res      = TRUE
  )
  
  orig_ok <- tryCatch(
    isTRUE(
      all.equal(
        terra::origin(x[[1]]),
        terra::origin(template[[1]]),
        tolerance = 1e-6
      )
    ),
    error = function(e) FALSE
  )
  
  if (!isTRUE(geom_ok) || !isTRUE(orig_ok)) {
    
    stop(
      label,
      " does not match the native ERA5-Land grid of the first file.\n",
      "Script 6 does not resample or reproject index NetCDF products.\n",
      "Check the upstream NetCDF producer and file set."
    )
  }
  
  invisible(TRUE)
}


# ==============================================================================
#   2. REGIONAL WEIGHT CACHE
# ==============================================================================

# Regional weights are cached ONLY by the geometry of the actual native grid.
#
# The cache does not represent a canonical grid. It is merely an in-memory
# performance optimization so that repeated scales sharing exactly the same
# geometry do not repeatedly rasterize the 7 subbasins.

REGIONAL_WEIGHT_CACHE <- new.env(parent = emptyenv())


get_regional_weights_cached <- function(r_tmpl) {
  
  key <- paste(
    terra::nrow(r_tmpl),
    terra::ncol(r_tmpl),
    paste(
      signif(
        as.vector(terra::ext(r_tmpl)),
        12
      ),
      collapse = ","
    ),
    trimws(terra::crs(r_tmpl)),
    paste(
      signif(
        terra::origin(r_tmpl),
        12
      ),
      collapse = ","
    ),
    sep = "|"
  )
  
  if (
    exists(
      key,
      envir = REGIONAL_WEIGHT_CACHE,
      inherits = FALSE
    )
  ) {
    return(
      get(
        key,
        envir = REGIONAL_WEIGHT_CACHE,
        inherits = FALSE
      )
    )
  }
  
  regions <- load_nechako_regions(
    template = r_tmpl
  )
  
  validate_nechako_regions(
    regions
  )
  
  rw <- region_fractional_weights(
    r_tmpl,
    regions
  )
  
  assign(
    key,
    rw,
    envir = REGIONAL_WEIGHT_CACHE
  )
  
  rw
}


# Regional averages are reported only when finite index values cover enough
# of the full polygon area. Snow-index summer gaps otherwise overrepresent the
# few high-elevation cells that remain valid.
MIN_REGIONAL_VALID_AREA_COVERAGE <- 0.50

regional_mean_with_coverage <- function(ts_matrix, area_weights, total_region_area_m2,
                                        min_coverage = MIN_REGIONAL_VALID_AREA_COVERAGE) {
  if (!is.matrix(ts_matrix) || length(area_weights) != nrow(ts_matrix)) {
    stop("regional_mean_with_coverage: matrix/weight dimensions do not match.")
  }
  w <- as.numeric(area_weights)
  if (!is.finite(total_region_area_m2) || total_region_area_m2 <= 0) {
    return(list(value = rep(NA_real_, ncol(ts_matrix)),
                coverage = rep(NA_real_, ncol(ts_matrix))))
  }
  coverage <- value <- rep(NA_real_, ncol(ts_matrix))
  for (j in seq_len(ncol(ts_matrix))) {
    ok <- is.finite(ts_matrix[, j]) & is.finite(w) & w > 0
    valid_area <- sum(w[ok], na.rm = TRUE)
    coverage[j] <- valid_area / total_region_area_m2
    if (valid_area > 0 && is.finite(coverage[j]) && coverage[j] >= min_coverage) {
      value[j] <- sum(ts_matrix[ok, j] * w[ok], na.rm = TRUE) / valid_area
    }
  }
  list(value = value, coverage = coverage)
}

# ==============================================================================
#   3. SCALE DEFINITIONS
# ==============================================================================

# Current manuscript scales.

SPI_SCALES <- c(
  1L,
  3L,
  6L,
  12L,
  24L,
  36L
)

SPEI_SCALES <- c(
  1L,
  3L,
  6L,
  12L,
  24L,
  36L
)

SPEI_ERA5_SCALES <- c(
  1L,
  3L,
  6L,
  12L,
  24L,
  36L
)

SPEI_ERA5_SEAS_DIR <- file.path(
  WD_PATH,
  "spei_results_seasonal_era5land"
)

# Snow drought index products generated by 8snow_drought_indices.R.
# SWEI is currently produced at scale 3; SSPI at scale 1.
SWEI_SCALES <- 3L
SSPI_SCALES <- 1L
SWEI_SEAS_DIR <- file.path(
  WD_PATH,
  "swei_results_seasonal"
)
SSPI_SEAS_DIR <- file.path(
  WD_PATH,
  "sspi_results_seasonal"
)


# ==============================================================================
#   4. SCALE-SPECIFIC MINIMUM EVENT DURATION
# ==============================================================================

get_min_duration <- function(
    scale,
    min_duration_map = list(
      "1"  = 2L,
      "3"  = 3L,
      "6"  = 4L,
      "12" = 6L,
      "24" = 12L,
      "36" = 12L
    ),
    default_min_duration = 3L
) {
  
  key <- as.character(
    as.integer(scale)
  )
  
  val <- min_duration_map[[key]]
  
  if (!is.null(val)) {
    return(
      as.integer(val)
    )
  }
  
  as.integer(
    default_min_duration
  )
}


# ==============================================================================
#   6. DROUGHT EVENT CHARACTERISTICS
# ==============================================================================

# Return event positions on the original calendar axis. Non-finite values end
# an active event and are never removed/spliced out of the time series.
extract_hysteresis_events <- function(x, onset_thr, end_thr, min_duration = 1L) {
  x <- as.numeric(x)
  events <- list()
  in_d <- FALSE
  start_idx <- NA_integer_
  close_event <- function(end_idx) {
    if (!is.na(start_idx) && end_idx >= start_idx) {
      duration <- end_idx - start_idx + 1L
      if (duration >= min_duration) {
        events[[length(events) + 1L]] <<- data.frame(
          start_index = as.integer(start_idx),
          end_index = as.integer(end_idx),
          duration_months = as.integer(duration)
        )
      }
    }
  }
  for (j in seq_along(x)) {
    v <- x[j]
    if (!is.finite(v)) {
      if (in_d) close_event(j - 1L)
      in_d <- FALSE
      start_idx <- NA_integer_
      next
    }
    if (!in_d && v < onset_thr) {
      in_d <- TRUE
      start_idx <- j
    } else if (in_d && v >= end_thr) {
      close_event(j - 1L)
      in_d <- FALSE
      start_idx <- NA_integer_
    }
  }
  if (in_d) close_event(length(x))
  if (!length(events)) return(data.frame(
    start_index = integer(), end_index = integer(), duration_months = integer()
  ))
  do.call(rbind, events)
}

vectorized_event_stats <- function(
    ts_matrix,
    onset_thr = DROUGHT_ONSET,
    end_thr = onset_thr,
    min_duration = 1L
) {
  
  if (end_thr < onset_thr) {
    
    warning(
      sprintf(
        paste0(
          "vectorized_event_stats: end_thr (%.2f) < ",
          "onset_thr (%.2f). For hysteresis, end_thr should ",
          "be >= onset_thr."
        ),
        end_thr,
        onset_thr
      )
    )
  }
  
  n_pix <- nrow(ts_matrix)
  
  n_events      <- integer(n_pix)
  mean_duration <- rep(NA_real_, n_pix)
  max_intensity <- rep(NA_real_, n_pix)
  
  n_events_D12  <- integer(n_pix)
  n_events_D36  <- integer(n_pix)
  n_events_D712 <- integer(n_pix)
  n_events_D13p <- integer(n_pix)
  n_D4p         <- integer(n_pix)
  
  mean_I <- numeric(n_pix)
  mean_S <- numeric(n_pix)
  
  for (i in seq_len(n_pix)) {
    x <- as.numeric(ts_matrix[i, ])
    ev <- extract_hysteresis_events(x, onset_thr, end_thr, min_duration)
    durs <- ev$duration_months
    ints <- numeric(0)
    I_evs <- numeric(0)
    S_evs <- numeric(0)
    if (nrow(ev)) {
      for (k in seq_len(nrow(ev))) {
        seg <- x[ev$start_index[k]:ev$end_index[k]]
        intensity <- end_thr - seg
        ints <- c(ints, max(intensity, na.rm = TRUE))
        I_ev <- mean(intensity, na.rm = TRUE)
        I_evs <- c(I_evs, I_ev)
        S_evs <- c(S_evs, I_ev * ev$duration_months[k])
      }
    }
    
    if (length(durs) == 0L) {
      
      n_events[i]      <- 0L
      mean_duration[i] <- NA_real_
      max_intensity[i] <- NA_real_
      
      n_events_D12[i]  <- 0L
      n_events_D36[i]  <- 0L
      n_events_D712[i] <- 0L
      n_events_D13p[i] <- 0L
      n_D4p[i]         <- 0L
      
      mean_I[i] <- NA_real_
      mean_S[i] <- NA_real_
      
    } else {
      
      n_events[i] <- length(durs)
      
      mean_duration[i] <- mean(
        durs
      )
      
      max_intensity[i] <- max(
        ints,
        na.rm = TRUE
      )
      
      n_events_D12[i] <- sum(durs >= 1L & durs <= 2L)
      n_events_D36[i] <- sum(
        durs >= 3L &
          durs <= 6L
      )
      
      n_events_D712[i] <- sum(
        durs >= 7L &
          durs <= 12L
      )
      
      n_events_D13p[i] <- sum(
        durs >= 13L
      )
      
      n_D4p[i] <- sum(
        durs >= 4L
      )
      
      mean_I[i] <- mean(
        I_evs,
        na.rm = TRUE
      )
      
      mean_S[i] <- mean(
        S_evs,
        na.rm = TRUE
      )
    }
  }
  
  data.table::data.table(
    n_events      = n_events,
    mean_duration = mean_duration,
    max_intensity = max_intensity,
    n_events_D12  = n_events_D12,
    n_events_D36  = n_events_D36,
    n_events_D712 = n_events_D712,
    n_events_D13p = n_events_D13p,
    n_D4p         = n_D4p,
    mean_I        = mean_I,
    mean_S        = mean_S
  )
}


# ==============================================================================
#   7. WALD-WOLFOWITZ RUNS TEST
# ==============================================================================

vectorized_runs_test <- function(
    ts_matrix,
    onset_thr = DROUGHT_ONSET,
    end_thr = onset_thr,
    min_duration = 1L
) {
  
  n_pix <- nrow(ts_matrix)
  
  res <- data.table::data.table(
    p_value_runs = rep(
      NA_real_,
      n_pix
    ),
    clustering = rep(
      NA_character_,
      n_pix
    ),
    filtered_runs = rep(
      TRUE,
      n_pix
    )
  )
  
  for (i in seq_len(n_pix)) {
    
    x <- as.numeric(ts_matrix[i, ])
    finite_pos <- which(is.finite(x))
    if (!length(finite_pos)) next
    x <- x[min(finite_pos):max(finite_pos)]
    # A runs test on a series with internal missing months would falsely join
    # non-adjacent observations. Leave it undefined for gappy series.
    if (any(!is.finite(x)) || length(x) < 10L) {
      next
    }
    
    # Build drought binary indicator using the exact same hysteresis
    # state machine used by vectorized_event_stats().
    
    b <- integer(
      length(x)
    )
    
    in_d  <- FALSE
    s_idx <- NA_integer_
    
    for (j in seq_along(x)) {
      
      v <- x[j]
      
      if (!in_d && v < onset_thr) {
        
        in_d  <- TRUE
        s_idx <- j
        
      } else if (
        in_d &&
        v >= end_thr
      ) {
        
        dur <- j - s_idx
        
        if (dur >= min_duration) {
          
          b[
            s_idx:(j - 1L)
          ] <- 1L
        }
        
        in_d  <- FALSE
        s_idx <- NA_integer_
      }
    }
    
    # Close an open event.
    
    if (
      in_d &&
      !is.na(s_idx)
    ) {
      
      dur <- length(x) - s_idx + 1L
      
      if (dur >= min_duration) {
        
        b[
          s_idx:length(x)
        ] <- 1L
      }
    }
    
    n  <- length(b)
    n1 <- sum(b)
    n0 <- n - n1
    
    if (
      n1 < 5L ||
      n0 < 5L
    ) {
      next
    }
    
    nr <- length(
      rle(b)$lengths
    )
    
    er <- 2.0 * n0 * n1 / n + 1.0
    
    vr <- (
      2.0 *
        n0 *
        n1 *
        (
          2.0 *
            n0 *
            n1 -
            n
        )
    ) /
      (
        n^2 *
          (n - 1.0)
      )
    
    if (vr <= 0) {
      next
    }
    
    z <- (
      nr - er
    ) /
      sqrt(vr)
    
    pv <- 2.0 *
      pnorm(
        -abs(z)
      )
    
    clus <- if (nr < er) {
      "clustered"
    } else if (nr > er) {
      "dispersed"
    } else {
      "random"
    }
    
    res[i, `:=`(
      p_value_runs = pv,
      clustering = clus,
      filtered_runs = FALSE
    )]
  }
  
  res[, p_value_runs_fdr := NA_real_]
  valid_p <- is.finite(res$p_value_runs)
  if (any(valid_p)) {
    res$p_value_runs_fdr[valid_p] <- stats::p.adjust(res$p_value_runs[valid_p], method = "BH")
    res$clustering[valid_p & (!is.finite(res$p_value_runs_fdr) | res$p_value_runs_fdr >= 0.05)] <- "not significant"
  }
  res$clustering[!valid_p] <- NA_character_
  res
}


# ==============================================================================
#   8. SPECTRAL PEAK DETECTION
# ==============================================================================
#
# Monte Carlo AR(1) red-noise null.
#
# A single maximum-statistic threshold is generated across the entire drought
# frequency band for each unique combination of series length and binned rho1.
#
# This maximum-statistic threshold addresses testing multiple spectral bins.
# The AR(1) null is an approximation and does not reproduce the exact moving-sum
# filter structure of SPI/SPEI/SWEI accumulation scales.
# ==============================================================================

n_sim_spectral <- 500L
SPECTRAL_RN_ENV_CACHE <- new.env(parent = emptyenv())


vectorized_spectral_peaks <- function(
    ts_matrix,
    min_period_months = 24L,
    sig_level = 0.95,
    n_sim_spectral = 500L
) {
  
  n_pix <- nrow(ts_matrix)
  
  out <- integer(
    n_pix
  )
  
  rn_env_cache <- SPECTRAL_RN_ENV_CACHE
  
  
  get_red_noise_threshold <- function(
    n_len,
    band,
    rho1
  ) {
    
    rho1_bin <- round(
      rho1 / 0.05
    ) * 0.05
    
    rho1_bin <- max(
      -0.95,
      min(
        0.95,
        rho1_bin
      )
    )
    
    key <- paste0(
      n_len,
      "_r",
      sprintf(
        "%.2f",
        rho1_bin
      )
    )
    
    if (!is.null(
      rn_env_cache[[key]]
    )) {
      return(
        rn_env_cache[[key]]
      )
    }
    
    had_rng_state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_rng_state) old_rng_state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (had_rng_state) {
        assign(".Random.seed", old_rng_state, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(42L)
    
    band_maxima <- vapply(
      seq_len(n_sim_spectral),
      function(s) {
        
        burn_in <- 500L
        innov <- stats::rnorm(
          n_len + burn_in,
          sd = sqrt(max(1e-8, 1 - rho1_bin^2))
        )
        rn_full <- stats::filter(
          innov,
          filter = rho1_bin,
          method = "recursive"
        )
        rn <- tail(as.numeric(rn_full), n_len)
        rn <- residuals(stats::lm(rn ~ seq_along(rn)))
        rn_var <- stats::var(rn)
        if (!is.finite(rn_var) || rn_var <= 0) return(NA_real_)
        rn <- rn / sqrt(rn_var)
        
        sp <- stats::spectrum(
          rn,
          spans = c(
            3L,
            5L
          ),
          taper = 0,
          detrend = FALSE,
          plot = FALSE
        )
        
        max(
          sp$spec[band],
          na.rm = TRUE
        )
      },
      numeric(1)
    )
    
    band_maxima <- band_maxima[is.finite(band_maxima)]
    if (!length(band_maxima)) stop("Red-noise simulation produced no valid spectral maxima.")
    thr <- stats::quantile(band_maxima, probs = sig_level, names = FALSE)
    
    rn_env_cache[[key]] <<- thr
    
    thr
  }
  
  
  for (i in seq_len(n_pix)) {
    
    x <- as.numeric(ts_matrix[i, ])
    finite_pos <- which(is.finite(x))
    if (!length(finite_pos)) next
    x <- x[min(finite_pos):max(finite_pos)]
    # Fourier spectra require contiguous equally spaced observations; do not
    # bridge missing summer months in SWEI/SSPI or any other internal gaps.
    if (any(!is.finite(x)) || length(x) < 48L) {
      next
    }
    
    # Demean and linearly detrend.
    
    x <- x -
      mean(x)
    
    tt <- seq_along(
      x
    )
    
    x <- residuals(
      lm(
        x ~ tt
      )
    )
    
    pv <- var(
      x
    )
    
    if (
      !is.finite(pv) ||
      pv <= 0
    ) {
      next
    }
    
    # Lag-1 autocorrelation of detrended series.
    
    n_x <- length(
      x
    )
    
    rho1 <- tryCatch(
      stats::cor(
        x[-n_x],
        x[-1],
        use = "complete.obs"
      ),
      error = function(e) 0.0
    )
    
    if (
      is.na(rho1) ||
      !is.finite(rho1)
    ) {
      rho1 <- 0.0
    }
    
    rho1 <- max(
      -0.95,
      min(
        0.95,
        rho1
      )
    )
    
    sp <- tryCatch(
      stats::spectrum(
        x,
        spans = c(
          3L,
          5L
        ),
        taper = 0,
        detrend = FALSE,
        plot = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(sp)) {
      next
    }
    
    periods <- 1.0 /
      sp$freq
    
    band <- periods >=
      min_period_months
    
    if (!any(band)) {
      next
    }
    
    thr_unit <- get_red_noise_threshold(
      length(x),
      band,
      rho1
    )
    
    threshold <- thr_unit *
      pv
    
    sig_bins <-
      band &
      (
        sp$spec >
          threshold
      )
    
    if (any(sig_bins)) {
      
      r <- rle(
        sig_bins
      )
      
      out[i] <-
        sum(
          r$values
        )
    }
  }
  
  data.table::data.table(
    n_spectral_peaks = out
  )
}


# ==============================================================================
#   9. BASIN EXTENT TIME SERIES
# ==============================================================================

compute_basin_extent_and_event_catalog <- function(
    ts_matrix,
    n_years,
    index_label,
    out_dir,
    onset_thr = DROUGHT_ONSET,
    end_thr = DROUGHT_END,
    start_year = 1950L,
    win_years = 30L,
    area_weights = NULL,
    min_duration = 1L
) {
  
  use_weights <-
    !is.null(area_weights) &&
    length(area_weights) ==
    nrow(ts_matrix) &&
    all(is.finite(area_weights)) &&
    all(area_weights > 0)
  
  if (
    !is.null(area_weights) &&
    !use_weights
  ) {
    
    warning(
      paste0(
        "compute_basin_extent_and_event_catalog: ",
        "area_weights ignored because length/values are invalid; ",
        "falling back to equal weighting."
      )
    )
  }
  
  utils_load_packages(
    c(
      "ggplot2",
      "Kendall",
      "dplyr",
      "patchwork"
    )
  )
  
  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  n_months <- n_years *
    12L
  
  month_seq <- seq_len(
    n_months
  )
  
  yr_seq <- start_year +
    floor(
      (month_seq - 1L) /
        12L
    )
  
  mo_seq <- (
    (month_seq - 1L) %%
      12L
  ) + 1L
  
  date_seq <- as.Date(
    sprintf(
      "%04d-%02d-01",
      yr_seq,
      mo_seq
    )
  )
  
  n_pix <- nrow(
    ts_matrix
  )
  
  n_col <- min(
    n_months,
    ncol(ts_matrix)
  )
  
  pct_area <- numeric(
    n_col
  )
  
  for (t in seq_len(n_col)) {
    
    col_vals <-
      ts_matrix[, t]
    
    ok <- !is.na(
      col_vals
    )
    
    if (!any(ok)) {
      
      pct_area[t] <- NA_real_
      
      next
    }
    
    if (use_weights) {
      
      w_ok <- area_weights[
        ok
      ]
      
      drought_ok <-
        ok &
        col_vals <
        onset_thr
      
      pct_area[t] <-
        100.0 *
        sum(
          area_weights[drought_ok],
          na.rm = TRUE
        ) /
        sum(
          w_ok,
          na.rm = TRUE
        )
      
    } else {
      
      pct_area[t] <-
        100.0 *
        sum(
          ok &
            col_vals <
            onset_thr
        ) /
        sum(ok)
    }
  }
  
  extent_ts <- data.frame(
    date = date_seq[
      seq_len(n_col)
    ],
    year = yr_seq[
      seq_len(n_col)
    ],
    month = mo_seq[
      seq_len(n_col)
    ],
    pct_area = pct_area,
    row.names = NULL
  )
  
  write.csv(
    extent_ts,
    file.path(
      out_dir,
      sprintf(
        "%s_basin_extent_monthly.csv",
        index_label
      )
    ),
    row.names = FALSE
  )
  
  extent_annual <-
    extent_ts |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      mean_pct = mean(
        pct_area,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  x_clean <-
    extent_annual$mean_pct[
      !is.na(
        extent_annual$mean_pct
      )
    ]
  
  mk_ext <- list(
    tau = NA_real_,
    p_value = NA_real_,
    sens_slope = NA_real_
  )
  
  if (
    length(x_clean) >= 10L
  ) {
    
    tau_val <- as.numeric(
      Kendall::MannKendall(
        x_clean
      )$tau
    )
    
    hr_r <- tryCatch(
      modifiedmk::mmkh(
        x_clean
      ),
      error = function(e) NULL
    )
    
    if (!is.null(hr_r)) {
      
      mk_ext <- list(
        tau = tau_val,
        p_value = as.numeric(
          hr_r[
            "new P-value"
          ]
        ),
        sens_slope = as.numeric(
          hr_r[
            "Sen's slope"
          ]
        )
      )
      
    } else {
      
      mk_r <-
        Kendall::MannKendall(
          x_clean
        )
      
      mk_ext <- list(
        tau = tau_val,
        p_value = as.numeric(
          mk_r$sl
        ),
        sens_slope = NA_real_
      )
    }
  }
  
  mk_ext_df <- data.frame(
    tau = mk_ext$tau,
    p_value = mk_ext$p_value,
    sens_slope = mk_ext$sens_slope
  )
  
  write.csv(
    mk_ext_df,
    file.path(
      out_dir,
      sprintf(
        "%s_basin_extent_MK.csv",
        index_label
      )
    ),
    row.names = FALSE
  )
  
  # Rolling-window summary.
  
  if (
    nrow(extent_annual) >=
    win_years
  ) {
    
    rolling <- data.frame(
      year = extent_annual$year,
      rolling_mean = zoo::rollmean(
        extent_annual$mean_pct,
        k = win_years,
        fill = NA,
        align = "right"
      )
    )
    
    write.csv(
      rolling,
      file.path(
        out_dir,
        sprintf(
          "%s_basin_extent_%dyr_rolling.csv",
          index_label,
          win_years
        )
      ),
      row.names = FALSE
    )
  }
  
  # Plot.
  
  p_ext <- ggplot2::ggplot(
    extent_annual,
    ggplot2::aes(
      x = year,
      y = mean_pct
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::labs(
      title = sprintf(
        "%s — Annual Mean Basin-Wide Drought Extent",
        index_label
      ),
      subtitle = sprintf(
        "Hamed–Rao MK: tau = %.3f, p = %.4g",
        mk_ext$tau,
        mk_ext$p_value
      ),
      x = "Year",
      y = "% Basin Area in Drought"
    ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  ggplot2::ggsave(
    file.path(
      out_dir,
      sprintf(
        "%s_basin_extent_timeseries.png",
        index_label
      )
    ),
    p_ext,
    width = 12,
    height = 5,
    dpi = 150
  )
  
  # Event catalogue retains original monthly indices. Missing months end events.
  event_list <- list()
  for (i in seq_len(nrow(ts_matrix))) {
    ev <- extract_hysteresis_events(
      ts_matrix[i, ], onset_thr, end_thr, min_duration = min_duration
    )
    if (nrow(ev)) {
      ev$pixel <- i
      event_list[[length(event_list) + 1L]] <- ev[, c(
        "pixel", "start_index", "end_index", "duration_months"
      )]
    }
  }
  
  if (length(event_list) == 0L) {
    
    cat(
      sprintf(
        "  [Event catalogue | %s] No events detected.\n",
        index_label
      )
    )
    
    return(
      invisible(
        list(
          extent_ts = extent_ts,
          mk_extent = mk_ext_df,
          events_df = NULL
        )
      )
    )
  }
  
  events_df <- do.call(
    rbind,
    event_list
  )
  
  events_df$duration_class <- cut(
    events_df$duration_months,
    breaks = c(0, 2, 6, 12, Inf),
    labels = c("D1-2", "D3-6", "D7-12", "D13+"),
    right = TRUE
  )
  events_df$index <- index_label
  events_df$minimum_duration_months <- as.integer(min_duration)
  
  write.csv(
    events_df,
    file.path(
      out_dir,
      sprintf(
        "%s_basin_event_catalog.csv",
        index_label
      )
    ),
    row.names = FALSE
  )
  
  invisible(
    list(
      extent_ts = extent_ts,
      mk_extent = mk_ext_df,
      events_df = events_df
    )
  )
}


# ==============================================================================
#   10. MAIN PROCESSING FUNCTION
# ==============================================================================

process_index <- function(
    index_type,
    sc,
    find_fn,
    seas_dir
) {
  
  cat(
    sprintf(
      "\n╔══ %s-%02d ══╗\n",
      toupper(index_type),
      sc
    )
  )
  
  min_dur <- get_min_duration(
    sc
  )
  
  cat(
    sprintf(
      "  Method 2 minimum duration: %d months\n",
      min_dur
    )
  )
  
  
  # --------------------------------------------------------------------------
  # [1/6] FIND FILES
  # --------------------------------------------------------------------------
  
  cat(
    "[1/6] Finding files...\n"
  )
  
  files <- find_fn(
    seas_dir,
    index_type,
    sc
  )
  
  if (!length(files)) {
    
    stop(
      sprintf(
        "No monthly NetCDF files found for %s-%02d in %s",
        toupper(index_type),
        sc,
        seas_dir
      )
    )
  }
  
  cat(
    sprintf(
      "  Found %d monthly NC files\n",
      length(files)
    )
  )
  
  
  # --------------------------------------------------------------------------
  # [2/6] LOAD RASTER STACKS
  # --------------------------------------------------------------------------
  
  cat(
    "[2/6] Loading rasters...\n"
  )
  
  stacks <- lapply(
    files,
    function(f) {
      
      tryCatch(
        terra::rast(f),
        error = function(e) {
          
          cat(
            "  ⚠ Failed: ",
            basename(f),
            " | ",
            conditionMessage(e),
            "\n",
            sep = ""
          )
          
          NULL
        }
      )
    }
  )
  
  failed_stack <- vapply(stacks, is.null, logical(1))
  if (any(failed_stack)) {
    stop(sprintf("Failed to load %d NetCDF raster(s) for %s-%02d: %s",
                 sum(failed_stack), toupper(index_type), sc,
                 paste(basename(files[failed_stack]), collapse = ", ")))
  }
  
  
  # --------------------------------------------------------------------------
  # NATIVE GRID DEFINITION
  # --------------------------------------------------------------------------
  #
  # IMPORTANT:
  # r_tmpl is taken directly from the first actual NetCDF file.
  #
  # It is NOT loaded from an RDS.
  # It is NOT compared with a canonical grid.
  # It is NOT resampled.
  # --------------------------------------------------------------------------
  
  r_tmpl <- stacks[[1]][[1]]
  
  if (
    is.na(
      terra::crs(r_tmpl)
    ) ||
    !nzchar(
      terra::crs(r_tmpl)
    )
  ) {
    
    stop(
      "Native ERA5-Land index grid has no CRS: ",
      basename(files[1])
    )
  }
  
  cat(
    sprintf(
      "  ✓ Native grid: %d × %d = %d cells\n",
      terra::nrow(r_tmpl),
      terra::ncol(r_tmpl),
      terra::ncell(r_tmpl)
    )
  )
  
  cat(
    sprintf(
      "  ✓ Native resolution: %.10f × %.10f\n",
      terra::res(r_tmpl)[1],
      terra::res(r_tmpl)[2]
    )
  )
  
  cat(
    sprintf(
      "  ✓ Native origin: %.10f, %.10f\n",
      terra::origin(r_tmpl)[1],
      terra::origin(r_tmpl)[2]
    )
  )
  
  cat(
    sprintf(
      "  ✓ Native CRS: %s\n",
      terra::crs(r_tmpl)
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Validate that every file in THIS run shares the same native geometry.
  # --------------------------------------------------------------------------
  
  for (k in seq_along(stacks)) {
    
    assert_same_native_grid(
      stacks[[k]],
      r_tmpl,
      label = sprintf(
        "%s-%02d file %d (%s)",
        toupper(index_type),
        sc,
        k,
        basename(files[k])
      )
    )
  }
  
  cat(
    "  ✓ All NetCDF files in this index × scale match the first native grid\n"
  )
  
  
  xy <- terra::xyFromCell(
    r_tmpl,
    seq_len(
      terra::ncell(r_tmpl)
    )
  )
  
  lon_vec <- xy[, 1]
  lat_vec <- xy[, 2]
  
  n_pix <- length(
    lon_vec
  )
  
  
  # --------------------------------------------------------------------------
  # [3/6] BUILD TIME-SERIES MATRIX
  # --------------------------------------------------------------------------
  
  cat(
    "[3/6] Building time-series matrix...\n"
  )
  
  n_files <- length(
    stacks
  )
  
  n_lyrs_per_file <-
    terra::nlyr(
      stacks[[1]]
    )
  
  
  if (
    n_files == 1L &&
    n_lyrs_per_file > 12L
  ) {
    
    # ------------------------------------------------------------------------
    # SINGLE-FILE MODE
    # ------------------------------------------------------------------------
    
    n_months_total <-
      n_lyrs_per_file
    
    dates_all <-
      extract_dates_from_nc(
        files[1],
        n_months_total
      )
    
    years <-
      sort(
        unique(
          as.integer(
            format(
              dates_all,
              "%Y"
            )
          )
        )
      )
    
    n_years <-
      length(years)
    
    cat(
      sprintf(
        "  Single-file mode: %d pixels × %d months (%d years)\n",
        n_pix,
        n_months_total,
        n_years
      )
    )
    
    ts_mat <- terra::values(
      stacks[[1]],
      mat = TRUE,
      na.rm = FALSE
    )
    
    ts_mat <- as.matrix(
      ts_mat
    )
    
    expected_dim <- c(
      as.integer(n_pix),
      as.integer(n_months_total)
    )
    
    actual_dim <- dim(
      ts_mat
    )
    
    if (
      length(actual_dim) != 2L ||
      any(
        as.integer(actual_dim) !=
        expected_dim
      )
    ) {
      
      stop(
        "Unexpected single-file value matrix dimensions: ",
        paste(
          dim(ts_mat),
          collapse = " × "
        ),
        " ; expected ",
        n_pix,
        " × ",
        n_months_total
      )
    }
    
    dates_ts <- dates_all
    
  } else {
    
    # ------------------------------------------------------------------------
    # STANDARD 12-FILE MODE
    # ------------------------------------------------------------------------
    #
    # Each file = one calendar month.
    # Each file contains one layer per year.
    #
    # Final layout:
    # Jan-1950, Feb-1950, ..., Dec-1950,
    # Jan-1951, ..., Dec-2025.
    # ------------------------------------------------------------------------
    
    if (n_files != 12L) {
      stop("Expected 12 calendar-month NetCDF files in standard mode; found ", n_files, ".")
    }
    
    month_from_name <- suppressWarnings(as.integer(sub(
      ".*month([0-9]{2}).*", "\\1", basename(files), ignore.case = TRUE
    )))
    if (length(month_from_name) != 12L || anyNA(month_from_name) ||
        !identical(sort(month_from_name), 1:12)) {
      stop("The 12 seasonal files must contain exactly one file for each month 01-12. Found: ",
           paste(basename(files), collapse = ", "))
    }
    file_order <- order(month_from_name)
    files <- files[file_order]
    stacks <- stacks[file_order]
    month_from_name <- month_from_name[file_order]
    
    dates_m1 <-
      extract_dates_from_nc(
        files[1],
        n_lyrs_per_file
      )
    
    years <-
      as.integer(
        format(
          dates_m1,
          "%Y"
        )
      )
    
    n_years <-
      length(years)
    
    cat(
      sprintf(
        "  Standard mode: %d pixels × %d years × 12 months\n",
        n_pix,
        n_years
      )
    )
    
    monthly_values <-
      vector(
        "list",
        n_files
      )
    
    for (m in seq_len(n_files)) {
      
      cat(
        sprintf(
          "    Reading month %02d/12: %s\n",
          m,
          basename(files[m])
        )
      )
      
      monthly_values[[m]] <-
        terra::values(
          stacks[[m]],
          mat = TRUE,
          na.rm = FALSE
        )
      
      monthly_values[[m]] <- as.matrix(monthly_values[[m]])
      
      dates_m <- extract_dates_from_nc(files[m], terra::nlyr(stacks[[m]]))
      if (length(dates_m) != n_years || anyNA(dates_m) ||
          any(as.integer(format(dates_m, "%m")) != m) ||
          !identical(as.integer(format(dates_m, "%Y")), years)) {
        stop("Date/year alignment failure in seasonal file ", basename(files[m]),
             ". Expected month ", m, " and the same year sequence as ", basename(files[1]), ".")
      }
      
      expected_dim <- c(as.integer(n_pix), as.integer(n_years))
      
      actual_dim <- dim(
        monthly_values[[m]]
      )
      
      if (
        length(actual_dim) != 2L ||
        any(
          as.integer(actual_dim) !=
          expected_dim
        )
      ) {
        
        stop(
          "Unexpected dimensions for month ",
          m,
          ": ",
          paste(
            actual_dim,
            collapse = " × "
          ),
          " ; expected ",
          paste(
            expected_dim,
            collapse = " × "
          )
        )
      }
    }
    
    ts_mat <- matrix(
      NA_real_,
      nrow = n_pix,
      ncol = n_years * 12L
    )
    
    for (m in seq_len(12L)) {
      
      vals_m <-
        monthly_values[[m]]
      
      n_yrs_m <- ncol(vals_m)
      
      for (y in seq_len(n_yrs_m)) {
        
        col_idx <-
          (y - 1L) * 12L +
          m
        
        ts_mat[, col_idx] <-
          vals_m[, y]
      }
    }
    
    dates_ts <- as.Date(
      sprintf(
        "%04d-%02d-01",
        rep(
          years,
          each = 12L
        ),
        rep(
          1:12,
          times = n_years
        )
      )
    )
    
    rm(
      monthly_values
    )
    
    gc(
      verbose = FALSE
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Time-series integrity
  # --------------------------------------------------------------------------
  
  if (!is.matrix(ts_mat)) {
    
    stop(
      "ts_mat is not a matrix."
    )
  }
  
  if (
    nrow(ts_mat) !=
    n_pix
  ) {
    
    stop(
      "ts_mat row count mismatch: ",
      nrow(ts_mat),
      " vs expected ",
      n_pix
    )
  }
  
  if (
    ncol(ts_mat) !=
    length(dates_ts)
  ) {
    
    stop(
      "ts_mat/date mismatch: ",
      ncol(ts_mat),
      " columns vs ",
      length(dates_ts),
      " dates"
    )
  }
  
  validate_monthly_record(
    dates_ts,
    label =
      "6trend time-series matrix calendar",
    require_complete = TRUE
  )
  
  cat(
    sprintf(
      "  ✓ Time-series matrix: %d pixels × %d months\n",
      nrow(ts_mat),
      ncol(ts_mat)
    )
  )
  
  cat(
    sprintf(
      "  ✓ Date range: %s to %s\n",
      min(dates_ts),
      max(dates_ts)
    )
  )
  
  finite_fraction <-
    mean(
      is.finite(ts_mat)
    )
  
  cat(
    sprintf(
      "  ✓ Finite-value fraction: %.2f%%\n",
      100 *
        finite_fraction
    )
  )
  
  
  # --------------------------------------------------------------------------
  # NATIVE-GRID BASIN MASK
  # --------------------------------------------------------------------------
  #
  # This is the critical replacement for the old canonical-grid block.
  #
  # The authoritative basin polygon is rasterized DIRECTLY on r_tmpl.
  # No stored raster mask is loaded.
  # --------------------------------------------------------------------------
  
  cat(
    "  → Building Nechako basin mask directly on native index grid...\n"
  )
  
  basin_v <- load_basin_vect(
    BASIN_SHP
  )
  
  if (
    !terra::same.crs(
      basin_v,
      r_tmpl
    )
  ) {
    
    basin_v <-
      terra::project(
        basin_v,
        terra::crs(r_tmpl)
      )
  }
  
  basin_coverage_r <- terra::rasterize(
    basin_v,
    r_tmpl,
    cover = TRUE
  )
  
  cover_all <- as.numeric(terra::values(basin_coverage_r, mat = FALSE))
  cover_all[!is.finite(cover_all)] <- 0
  if (any(cover_all < -1e-6 | cover_all > 1 + 1e-6, na.rm = TRUE)) {
    stop("Basin coverage raster contains values outside [0, 1]. Check terra::rasterize(cover=TRUE); refusing to clamp invalid coverage.")
  }
  cover_all <- pmax(0, pmin(1, cover_all))
  
  basin_mask_grid <-
    cover_all > 0
  
  valid_pix <-
    rowSums(
      is.finite(ts_mat)
    ) > 0
  
  valid_pix <-
    valid_pix &
    basin_mask_grid
  
  
  n_basin_cells <-
    sum(
      basin_mask_grid,
      na.rm = TRUE
    )
  
  cat(
    sprintf(
      "  ✓ Native-grid basin mask: %d cells\n",
      n_basin_cells
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Fractional cell-area weights
  # --------------------------------------------------------------------------
  
  cat(
    "  → Computing native-grid fractional area weights...\n"
  )
  
  cell_area_rast <-
    terra::cellSize(
      r_tmpl,
      unit = "m"
    )
  
  area_all <-
    as.numeric(
      terra::values(
        cell_area_rast,
        na.rm = FALSE
      )
    )
  
  basin_area_raster_m2 <- sum(area_all * cover_all, na.rm = TRUE)
  basin_area_vector_m2 <- tryCatch(
    as.numeric(terra::expanse(basin_v, unit = "m", transform = TRUE)),
    error = function(e) NA_real_
  )
  if (length(basin_area_vector_m2) > 1L) basin_area_vector_m2 <- sum(basin_area_vector_m2, na.rm = TRUE)
  if (is.finite(basin_area_vector_m2) && basin_area_vector_m2 > 0) {
    area_rel_diff <- abs(basin_area_raster_m2 - basin_area_vector_m2) / basin_area_vector_m2
    cat(sprintf("  Basin area check: raster %.3f km²; polygon %.3f km²; relative difference %.2f%%\n",
                basin_area_raster_m2 / 1e6, basin_area_vector_m2 / 1e6, 100 * area_rel_diff))
    if (area_rel_diff > 0.05) {
      stop(sprintf("Rasterized basin area differs from polygon area by %.2f%% (>5%%). Check native-grid coverage.", 100 * area_rel_diff))
    }
  } else {
    warning("Could not calculate polygon area for the basin coverage sanity check.")
  }
  
  eff_area_all <- area_all * cover_all
  
  area_valid <-
    eff_area_all[
      valid_pix
    ]
  
  bad_w <-
    is.na(area_valid) |
    !is.finite(area_valid) |
    area_valid <= 0
  
  if (any(bad_w)) {
    
    stop(
      "Native-grid area-weighting failure for ",
      index_type,
      "-",
      sc,
      ": ",
      sum(bad_w),
      " in-basin pixels have zero/invalid effective area."
    )
  }
  
  n_partial <-
    sum(
      cover_all[valid_pix] > 0 &
        cover_all[valid_pix] < 1,
      na.rm = TRUE
    )
  
  cat(
    sprintf(
      "  ✓ Fractional area weights: %.0f–%.0f m² (median %.0f m²)\n",
      min(area_valid),
      max(area_valid),
      median(area_valid)
    )
  )
  
  cat(
    sprintf(
      "  ✓ Boundary pixels with partial coverage: %d\n",
      n_partial
    )
  )
  
  
  # Apply native-grid basin mask.
  
  ts_mat <-
    ts_mat[
      valid_pix,
      ,
      drop = FALSE
    ]
  
  lon_vec <-
    lon_vec[
      valid_pix
    ]
  
  lat_vec <-
    lat_vec[
      valid_pix
    ]
  
  n_pix <-
    sum(valid_pix)
  
  if (
    n_pix == 0L
  ) {
    
    stop(
      "No finite-data pixels remain inside the Nechako basin on the native ",
      index_type,
      " grid."
    )
  }
  
  cat(
    sprintf(
      "  ✓ Valid native-grid basin pixels: %d\n",
      n_pix
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Regional aggregation
  # --------------------------------------------------------------------------
  
  cat(
    "  → Computing regional native-grid weights...\n"
  )
  
  regional_weights_all <-
    get_regional_weights_cached(
      r_tmpl
    )
  
  regional_index_rows <-
    lapply(
      names(regional_weights_all),
      function(rid) {
        
        rw <-
          regional_weights_all[[rid]]
        
        w_all <- as.numeric(rw$weights)
        w <- w_all[valid_pix]
        total_region_area <- sum(w_all[is.finite(w_all) & w_all > 0], na.rm = TRUE)
        regional_stats <- regional_mean_with_coverage(ts_mat, w, total_region_area)
        data.table::data.table(
          region_id = rid,
          region_name = rw$region_name,
          date = dates_ts,
          value = regional_stats$value,
          valid_area_coverage_fraction = regional_stats$coverage,
          coverage_area_m2 = total_region_area,
          n_cells = rw$n_cells,
          n_partial_cells = rw$n_partial
        )
      }
    )
  
  regional_index <-
    data.table::rbindlist(
      regional_index_rows,
      use.names = TRUE,
      fill = TRUE
    )
  
  regional_index[
    ,
    `:=`(
      index_type =
        tolower(index_type),
      scale =
        as.integer(sc)
    )
  ]
  
  regional_index_dir <-
    file.path(
      TREND_DIR,
      "regional_timeseries"
    )
  
  dir.create(
    regional_index_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  region_catalog_file <-
    file.path(
      regional_index_dir,
      "region_catalog.csv"
    )
  
  if (!file.exists(
    region_catalog_file
  )) {
    
    data.table::fwrite(
      region_catalog(
        load_nechako_regions(
          template = r_tmpl
        )
      ),
      region_catalog_file
    )
  }
  
  data.table::fwrite(
    regional_index,
    file.path(
      regional_index_dir,
      sprintf(
        "%s_%02d_regional_averaged_by_month.csv",
        index_type,
        sc
      )
    )
  )
  
  cat(
    sprintf(
      "  ✓ Regional index means written: %d regions × %d months\n",
      length(regional_weights_all),
      length(dates_ts)
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Basin-average time series
  # --------------------------------------------------------------------------
  
  basin_total_area_m2 <- sum(area_all * cover_all, na.rm = TRUE)
  basin_regional_stats <- regional_mean_with_coverage(
    ts_mat, area_valid, basin_total_area_m2
  )
  basin_avg_vec <- as.numeric(basin_regional_stats$value)
  
  
  # --------------------------------------------------------------------------
  # [4/6] STATISTICS
  # --------------------------------------------------------------------------
  
  cat(
    "[4/6] Computing statistics...\n"
  )
  
  results <-
    data.table::data.table(
      lon = lon_vec,
      lat = lat_vec
    )
  
  
  # --------------------------------------------------------------------------
  # Variance-corrected MK
  # --------------------------------------------------------------------------
  
  cat(
    "  → Mann-Kendall (variance-corrected)...    "
  )
  
  results <-
    cbind(
      results,
      vectorized_mann_kendall(
        ts_mat
      )
    )
  
  results[
    ,
    `:=`(
      tau_vc_raw =
        tau_vc,
      p_value_vc_raw =
        p_value_vc,
      sl_vc_raw =
        sl_vc
    )
  ]
  
  cat(
    "✓\n"
  )
  
  
  # --------------------------------------------------------------------------
  # TFPW MK
  # --------------------------------------------------------------------------
  
  cat(
    "  → Mann-Kendall (TFPW)...    "
  )
  
  results <-
    cbind(
      results,
      vectorized_mann_kendall_tfpw(
        ts_mat
      )
    )
  
  results[
    ,
    `:=`(
      tau_tfpw_raw =
        tau_tfpw,
      p_value_tfpw_raw =
        p_value_tfpw,
      sl_tfpw_raw =
        sl_tfpw,
      tfpw_method =
        "Yue2002_shared_J4"
    )
  ]
  
  cat(
    "✓\n"
  )
  
  
  # --------------------------------------------------------------------------
  # Method 1 — Sheffield & Wood
  # --------------------------------------------------------------------------
  
  cat(
    "  → Event statistics [M1: S&W]...    "
  )
  
  m1_stats <-
    vectorized_event_stats(
      ts_mat,
      onset_thr =
        DROUGHT_ONSET,
      end_thr =
        DROUGHT_ONSET,
      min_duration = 1L
    )
  
  results <-
    cbind(
      results,
      m1_stats
    )
  
  cat(
    "✓\n"
  )
  
  
  # --------------------------------------------------------------------------
  # Method 2 — Hysteresis
  # --------------------------------------------------------------------------
  
  cat(
    sprintf(
      "  → Event statistics [M2: Hysteresis, min_dur=%d m]...    ",
      min_dur
    )
  )
  
  m2_stats <-
    vectorized_event_stats(
      ts_mat,
      onset_thr =
        DROUGHT_ONSET,
      end_thr =
        DROUGHT_END,
      min_duration =
        min_dur
    )
  
  data.table::setnames(
    m2_stats,
    paste0(
      names(m2_stats),
      "_hyst"
    )
  )
  
  results <-
    cbind(
      results,
      m2_stats
    )
  
  cat(
    "✓\n"
  )
  
  
  # --------------------------------------------------------------------------
  # Event summary
  # --------------------------------------------------------------------------
  
  cat(
    "  ── Method 1 (S&W) per-pixel summaries:\n"
  )
  
  cat(sprintf("     D1-2 events  — mean: %.2f\n", mean(results$n_events_D12, na.rm = TRUE)))
  cat(
    sprintf(
      "     D3-6 events  — mean: %.2f\n",
      mean(
        results$n_events_D36,
        na.rm = TRUE
      )
    )
  )
  
  cat(
    sprintf(
      "     D7-12 events — mean: %.2f\n",
      mean(
        results$n_events_D712,
        na.rm = TRUE
      )
    )
  )
  
  cat(
    sprintf(
      "     D13+ events  — mean: %.2f (max: %d)\n",
      mean(
        results$n_events_D13p,
        na.rm = TRUE
      ),
      max(
        results$n_events_D13p,
        na.rm = TRUE
      )
    )
  )
  
  cat(
    sprintf(
      "     mean_I: %.4f | mean_S: %.4f\n",
      mean(
        results$mean_I,
        na.rm = TRUE
      ),
      mean(
        results$mean_S,
        na.rm = TRUE
      )
    )
  )
  
  
  cat(
    "  ── Method 2 (Hysteresis) per-pixel summaries:\n"
  )
  
  cat(sprintf("     D1-2 events  — mean: %.2f\n", mean(results$n_events_D12_hyst, na.rm = TRUE)))
  cat(
    sprintf(
      "     D3-6 events  — mean: %.2f\n",
      mean(
        results$n_events_D36_hyst,
        na.rm = TRUE
      )
    )
  )
  
  cat(
    sprintf(
      "     D7-12 events — mean: %.2f\n",
      mean(
        results$n_events_D712_hyst,
        na.rm = TRUE
      )
    )
  )
  
  cat(
    sprintf(
      "     D13+ events  — mean: %.2f (max: %d)\n",
      mean(
        results$n_events_D13p_hyst,
        na.rm = TRUE
      ),
      max(
        results$n_events_D13p_hyst,
        na.rm = TRUE
      )
    )
  )
  
  cat(
    sprintf(
      "     mean_I_hyst: %.4f | mean_S_hyst: %.4f\n",
      mean(
        results$mean_I_hyst,
        na.rm = TRUE
      ),
      mean(
        results$mean_S_hyst,
        na.rm = TRUE
      )
    )
  )
  
  
  # --------------------------------------------------------------------------
  # PELT regime shift
  # --------------------------------------------------------------------------
  
  cat(
    "  → PELT regime-shift detection...    "
  )
  
  results <-
    cbind(
      results,
      vectorized_regime_shift_pelt(
        ts_mat,
        years
      )
    )
  
  cat(
    "✓\n"
  )
  
  
  # --------------------------------------------------------------------------
  # Runs test
  # --------------------------------------------------------------------------
  
  cat(
    sprintf(
      "  → Runs test [Hyst drought clusters, min_dur=%d m]...    ",
      min_dur
    )
  )
  
  results <-
    cbind(
      results,
      vectorized_runs_test(
        ts_mat,
        onset_thr =
          DROUGHT_ONSET,
        end_thr =
          DROUGHT_END,
        min_duration =
          min_dur
      )
    )
  
  cat(
    "✓\n"
  )
  
  
  # --------------------------------------------------------------------------
  # Spectral peaks
  # --------------------------------------------------------------------------
  
  cat(
    "  → Spectral peak detection ",
    "(periods >= 24 mo, MC AR(1) red-noise, n_sim=500)...    ",
    sep = ""
  )
  
  results <-
    cbind(
      results,
      vectorized_spectral_peaks(
        ts_mat,
        n_sim_spectral =
          n_sim_spectral
      )
    )
  
  n_pix_with_peaks <-
    sum(
      results$n_spectral_peaks > 0L,
      na.rm = TRUE
    )
  
  cat(
    sprintf(
      "✓ (%d/%d pixels have >=1 significant peak)\n",
      n_pix_with_peaks,
      n_pix
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Basin-average CSV
  # --------------------------------------------------------------------------
  
  cat(
    "  → Saving area-weighted basin-average CSV...\n"
  )
  
  seas_dir_map <- list(
    spi = SPI_SEAS_DIR,
    spei = SPEI_SEAS_DIR,
    spei_era5 = SPEI_ERA5_SEAS_DIR,
    swei = SWEI_SEAS_DIR,
    sspi = SSPI_SEAS_DIR
  )
  
  basin_csv_dir <-
    seas_dir_map[[tolower(index_type)]]
  
  if (
    is.null(basin_csv_dir) ||
    !nzchar(basin_csv_dir)
  ) {
    
    basin_csv_dir <-
      file.path(
        WD_PATH,
        paste0(
          index_type,
          "_results_seasonal"
        )
      )
  }
  
  dir.create(
    basin_csv_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  # Write the basin series already calculated with the full-basin valid-area
  # coverage threshold. Do not recompute an unthresholded mean via the utility.
  basin_years <- as.integer(format(dates_ts, "%Y"))
  basin_months <- as.integer(format(dates_ts, "%m"))
  unique_years <- sort(unique(basin_years))
  to_month_wide <- function(values) {
    wide <- data.frame(Year = unique_years)
    for (m in 1:12) {
      hit <- which(basin_months == m)
      col <- rep(NA_real_, length(unique_years))
      col[match(basin_years[hit], unique_years)] <- values[hit]
      wide[[month.abb[m]]] <- col
    }
    wide
  }
  basin_avg_path <- file.path(
    basin_csv_dir,
    sprintf("%s_%02d_basin_averaged_by_month.csv", tolower(index_type), sc)
  )
  basin_coverage_path <- file.path(
    basin_csv_dir,
    sprintf("%s_%02d_basin_valid_area_coverage_by_month.csv", tolower(index_type), sc)
  )
  utils::write.csv(to_month_wide(basin_avg_vec), basin_avg_path,
                   row.names = FALSE, na = "")
  utils::write.csv(to_month_wide(basin_regional_stats$coverage), basin_coverage_path,
                   row.names = FALSE, na = "")
  if (!file.exists(basin_avg_path) || is.na(file.info(basin_avg_path)$size) ||
      file.info(basin_avg_path)$size <= 0) {
    stop("Basin-average CSV failed for ", index_type, "-", sc, call. = FALSE)
  }
  cat(sprintf("  ✓ Basin-average CSV written (minimum %.0f%% valid area): %s\n",
              100 * MIN_REGIONAL_VALID_AREA_COVERAGE, basename(basin_avg_path)))
  
  # --------------------------------------------------------------------------
  # Basin extent
  # --------------------------------------------------------------------------
  
  cat(
    "  → Basin extent time series & event catalogue...    "
  )
  
  index_label <-
    sprintf(
      "%s-%02d",
      toupper(index_type),
      sc
    )
  
  extent_out <-
    file.path(
      TREND_DIR,
      "basin_extent"
    )
  
  tryCatch(
    {
      
      compute_basin_extent_and_event_catalog(
        ts_mat,
        n_years,
        index_label,
        extent_out,
        start_year =
          min(years),
        area_weights =
          area_valid,
        min_duration = min_dur
      )
    },
    error = function(e) {
      
      cat(sprintf("\n  ✗ Basin extent analysis failed: %s\n", conditionMessage(e)))
      stop(sprintf("Basin extent analysis failed for %s-%02d: %s",
                   index_type, sc, conditionMessage(e)), call. = FALSE)
    }
  )
  
  cat(
    "✓\n"
  )
  
  
  # --------------------------------------------------------------------------
  # [4c/6] BH FDR
  # --------------------------------------------------------------------------
  
  cat(
    "  → BH False Discovery Rate correction...\n"
  )
  
  apply_fdr <- function(
    p_raw
  ) {
    
    ok <- !is.na(
      p_raw
    )
    
    p_fdr <-
      rep(
        NA_real_,
        length(p_raw)
      )
    
    if (
      sum(ok) > 1L
    ) {
      
      p_fdr[ok] <-
        p.adjust(
          p_raw[ok],
          method = "BH"
        )
    }
    
    p_fdr
  }
  
  results$p_fdr_vc_raw <-
    apply_fdr(
      results$p_value_vc_raw
    )
  
  results$p_fdr_tfpw_raw <-
    apply_fdr(
      results$p_value_tfpw_raw
    )
  
  results$p_fdr_vc <-
    apply_fdr(
      results$p_value_vc
    )
  
  results$p_fdr_tfpw <-
    apply_fdr(
      results$p_value_tfpw
    )
  
  
  report_fdr <- function(
    label,
    p_raw,
    p_fdr
  ) {
    
    n_test <-
      sum(
        !is.na(p_raw)
      )
    
    n_nom <-
      sum(
        p_raw < 0.05,
        na.rm = TRUE
      )
    
    n_fdr <-
      sum(
        p_fdr < 0.05,
        na.rm = TRUE
      )
    
    cat(
      sprintf(
        "  ✓ FDR (%s): %d nominal → %d FDR-significant (%+d pixels)\n",
        label,
        n_nom,
        n_fdr,
        n_fdr - n_nom
      )
    )
    
    list(
      n_test = n_test,
      n_nom = n_nom,
      n_fdr = n_fdr
    )
  }
  
  fdr_vc <-
    report_fdr(
      "VC-HR98",
      results$p_value_vc,
      results$p_fdr_vc
    )
  
  fdr_tfpw <-
    report_fdr(
      "TFPW",
      results$p_value_tfpw,
      results$p_fdr_tfpw
    )
  
  
  # --------------------------------------------------------------------------
  # [5/6] MAIN RESULTS CSV
  # --------------------------------------------------------------------------
  
  out_file <-
    file.path(
      TREND_DIR,
      sprintf(
        "%s_%02d_results.csv",
        index_type,
        sc
      )
    )
  
  data.table::fwrite(
    results,
    out_file,
    showProgress = FALSE
  )
  
  cat(
    sprintf(
      "  ✓ %s (%.2f MB)\n",
      basename(out_file),
      file.info(out_file)$size /
        1024^2
    )
  )
  
  cat(
    sprintf(
      "  %-35s nominal: %d / %d | FDR: %d / %d\n",
      "Significant VC trends (p < 0.05):",
      fdr_vc$n_nom,
      fdr_vc$n_test,
      fdr_vc$n_fdr,
      fdr_vc$n_test
    )
  )
  
  cat(
    sprintf(
      "  %-35s nominal: %d / %d | FDR: %d / %d\n",
      "Significant TFPW trends (p < 0.05):",
      fdr_tfpw$n_nom,
      fdr_tfpw$n_test,
      fdr_tfpw$n_fdr,
      fdr_tfpw$n_test
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Duration-class CSV — S&W
  # --------------------------------------------------------------------------
  
  dur_cols_sw <- c(
    "lon",
    "lat",
    "n_events_D12",
    "n_events_D36",
    "n_events_D712",
    "n_events_D13p",
    "n_D4p",
    "mean_I",
    "mean_S"
  )
  
  data.table::fwrite(
    results[
      ,
      ..dur_cols_sw
    ],
    file.path(
      TREND_DIR,
      sprintf(
        "%s_%02d_duration_class_map_SW.csv",
        index_type,
        sc
      )
    ),
    showProgress = FALSE
  )
  
  cat(
    sprintf(
      "  ✓ Duration-class map: %s_%02d_duration_class_map_SW.csv\n",
      index_type,
      sc
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Duration-class CSV — Hysteresis
  # --------------------------------------------------------------------------
  
  dur_cols_hyst <- c(
    "lon",
    "lat",
    "n_events_D12_hyst",
    "n_events_D36_hyst",
    "n_events_D712_hyst",
    "n_events_D13p_hyst",
    "n_D4p_hyst",
    "mean_I_hyst",
    "mean_S_hyst"
  )
  
  data.table::fwrite(
    results[
      ,
      ..dur_cols_hyst
    ],
    file.path(
      TREND_DIR,
      sprintf(
        "%s_%02d_duration_class_map_Hyst.csv",
        index_type,
        sc
      )
    ),
    showProgress = FALSE
  )
  
  cat(
    sprintf(
      "  ✓ Duration-class map: %s_%02d_duration_class_map_Hyst.csv\n",
      index_type,
      sc
    )
  )
  
  
  # --------------------------------------------------------------------------
  # D13+ frequency rasters
  # --------------------------------------------------------------------------
  
  write_d13p_nc <- function(
    values_vec,
    suffix,
    label_suffix
  ) {
    
    r_out <-
      r_tmpl
    
    terra::values(
      r_out
    ) <- NA_real_
    
    all_vals <-
      rep(
        NA_real_,
        terra::ncell(r_tmpl)
      )
    
    all_vals[
      valid_pix
    ] <- values_vec
    
    terra::values(
      r_out
    ) <- all_vals
    
    names(r_out) <-
      sprintf(
        "%s_%02d_n_D13p_%s",
        index_type,
        sc,
        label_suffix
      )
    
    nc_path <-
      file.path(
        TREND_DIR,
        sprintf(
          "%s_%02d_D13p_pixel_frequency_%s.nc",
          index_type,
          sc,
          suffix
        )
      )
    
    terra::writeCDF(
      r_out,
      nc_path,
      overwrite = TRUE
    )
    
    cat(
      sprintf(
        "  ✓ D13+ raster (%s): %s\n",
        suffix,
        basename(nc_path)
      )
    )
  }
  
  
  write_d13p_nc(
    results$n_events_D13p,
    "SW",
    "SW_longterm_droughts"
  )
  
  write_d13p_nc(
    results$n_events_D13p_hyst,
    "Hyst",
    "Hyst_longterm_droughts"
  )
  
  
  required_outputs <- c(
    out_file,
    file.path(TREND_DIR, sprintf("%s_%02d_duration_class_map_SW.csv", index_type, sc)),
    file.path(TREND_DIR, sprintf("%s_%02d_duration_class_map_Hyst.csv", index_type, sc)),
    file.path(TREND_DIR, sprintf("%s_%02d_D13p_pixel_frequency_SW.nc", index_type, sc)),
    file.path(TREND_DIR, sprintf("%s_%02d_D13p_pixel_frequency_Hyst.nc", index_type, sc)),
    file.path(regional_index_dir, sprintf("%s_%02d_regional_averaged_by_month.csv", index_type, sc)),
    file.path(basin_csv_dir, sprintf("%s_%02d_basin_averaged_by_month.csv", tolower(index_type), sc)),
    file.path(basin_csv_dir, sprintf("%s_%02d_basin_valid_area_coverage_by_month.csv", tolower(index_type), sc)),
    file.path(extent_out, sprintf("%s_basin_extent_monthly.csv", index_label)),
    file.path(extent_out, sprintf("%s_basin_extent_MK.csv", index_label)),
    file.path(extent_out, sprintf("%s_basin_extent_timeseries.png", index_label))
  )
  missing_outputs <- required_outputs[!file.exists(required_outputs) |
                                        is.na(file.info(required_outputs)$size) | file.info(required_outputs)$size <= 0]
  if (length(missing_outputs)) {
    stop("Required output validation failed for ", index_type, "-", sc, ":\n",
         paste(missing_outputs, collapse = "\n"), call. = FALSE)
  }
  
  rm(
    ts_mat,
    stacks
  )
  
  invisible(
    gc()
  )
  
  invisible(
    TRUE
  )
}


# ==============================================================================
#   11. RUN ALL INDICES
# ==============================================================================

cat(
  "\n╔════════════════════════════════════════════════╗\n",
  "║  DROUGHT TREND TEST                           ║\n",
  "║  Native ERA5-Land grid                        ║\n",
  "║  Dual methods: S&W + Hysteresis              ║\n",
  "╚════════════════════════════════════════════════╝\n\n",
  sep = ""
)

cat(
  "  Indices : SPI, SPEI_PM, SPEI_ERA5, SWEI-3, SSPI-1\n"
)

cat(
  "  Scales  : ",
  paste(
    SPI_SCALES,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  sprintf(
    "  Method 1: onset < %.2f; termination >= %.2f; no min duration\n",
    DROUGHT_ONSET,
    DROUGHT_ONSET
  )
)

cat(
  sprintf(
    "  Method 2: onset < %.2f; termination >= %.2f; scale-specific duration\n",
    DROUGHT_ONSET,
    DROUGHT_END
  )
)

cat(
  "  Native-grid rule: every index × scale derives its basin mask directly\n",
  "                     from its own NetCDF raster geometry.\n",
  sprintf("  Regional means require >= %.0f%% valid area coverage.\n\n", 100 * MIN_REGIONAL_VALID_AREA_COVERAGE),
  sep = ""
)


total_start <-
  Sys.time()


# ==============================================================================
#   12. RUN STATUS
# ==============================================================================

RUN_STATUS <-
  data.table::data.table(
    index_type = character(),
    scale = integer(),
    success = logical(),
    message = character()
  )


run_index_safe <- function(
    index_type,
    scales,
    find_fn,
    seas_dir
) {
  
  for (sc in scales) {
    
    ok <- FALSE
    msg <- ""
    
    tryCatch(
      {
        
        process_index(
          index_type,
          sc,
          find_fn,
          seas_dir
        )
        
        ok <- TRUE
        msg <- "completed"
        
      },
      error = function(e) {
        
        msg <<-
          conditionMessage(e)
        
        cat(
          sprintf(
            "\n  ❌ %s-%02d failed: %s\n",
            toupper(index_type),
            sc,
            msg
          )
        )
      }
    )
    
    RUN_STATUS <<-
      data.table::rbindlist(
        list(
          RUN_STATUS,
          data.table::data.table(
            index_type =
              tolower(index_type),
            scale =
              as.integer(sc),
            success =
              ok,
            message =
              msg
          )
        ),
        use.names = TRUE
      )
  }
}


# ==============================================================================
#   13. SPI
# ==============================================================================

cat(
  "\n── SPI ──\n"
)

run_index_safe(
  "spi",
  SPI_SCALES,
  find_seasonal_nc_files,
  SPI_SEAS_DIR
)


# ==============================================================================
#   14. SPEI — PENMAN-MONTEITH
# ==============================================================================

cat(
  "\n── SPEI (PM) ──\n"
)

run_index_safe(
  "spei",
  SPEI_SCALES,
  find_seasonal_nc_files,
  SPEI_SEAS_DIR
)


# ==============================================================================
#   15. SPEI — ERA5-LAND NATIVE PEV
# ==============================================================================

find_spei_era5_nc_safe <- function(
    data_dir,
    index_type,
    scale
) {
  
  files <-
    tryCatch(
      suppressWarnings(
        find_seasonal_nc_files(
          data_dir,
          "spei",
          scale
        )
      ),
      error = function(e) {
        
        cat(
          sprintf(
            "  ⚠ SPEI_ERA5 file finder error at scale-%02d: %s\n",
            scale,
            conditionMessage(e)
          )
        )
        
        character(0)
      }
    )
  
  if (!length(files)) {
    
    cat(
      sprintf(
        paste0(
          "  ⚠ SPEI_ERA5 scale-%02d: no files matching ",
          "spei_%02d_month##_*.nc in %s\n"
        ),
        scale,
        scale,
        normalizePath(
          data_dir,
          winslash = "/",
          mustWork = FALSE
        )
      )
    )
    
  } else {
    
    cat(
      sprintf(
        paste0(
          "  ✓ SPEI_ERA5 scale-%02d: found %d seasonal ",
          "NetCDF files using prefix 'spei'\n"
        ),
        scale,
        length(files)
      )
    )
  }
  
  files
}


if (
  dir.exists(
    SPEI_ERA5_SEAS_DIR
  )
) {
  
  cat(
    "\n── SPEI (ERA5-Land PEV) ──\n"
  )
  
  run_index_safe(
    "spei_era5",
    SPEI_ERA5_SCALES,
    find_spei_era5_nc_safe,
    SPEI_ERA5_SEAS_DIR
  )
  
} else {
  
  cat(
    "\n── SPEI (ERA5-Land PEV): directory not found ──\n"
  )
  
  cat(
    sprintf(
      "   Expected: %s\n",
      SPEI_ERA5_SEAS_DIR
    )
  )
  
  cat(
    "   Run 3SPEI_ERALand.R first.\n"
  )
  
  for (
    sc in SPEI_ERA5_SCALES
  ) {
    
    RUN_STATUS <-
      data.table::rbindlist(
        list(
          RUN_STATUS,
          data.table::data.table(
            index_type =
              "spei_era5",
            scale =
              as.integer(sc),
            success =
              FALSE,
            message =
              "ERA5-Land SPEI directory not found"
          )
        ),
        use.names = TRUE
      )
  }
}


# ==============================================================================
#   16. SNOW DROUGHT INDICES — SWEI AND SSPI
# ==============================================================================

# SWEI seasonal files use the dedicated finder because their filename prefix
# and naming convention are index-specific. SSPI uses the standard seasonal
# naming convention (sspi_01_month##_*.nc).

if (dir.exists(SWEI_SEAS_DIR)) {
  cat("\n── SWEI (Snow Water Equivalent Index) ──\n")
  run_index_safe(
    "swei",
    SWEI_SCALES,
    find_swei_seasonal_files,
    SWEI_SEAS_DIR
  )
} else {
  cat(sprintf("\n── SWEI skipped: seasonal directory not found: %s ──\n", SWEI_SEAS_DIR))
  for (sc in SWEI_SCALES) {
    RUN_STATUS <- data.table::rbindlist(list(
      RUN_STATUS,
      data.table::data.table(
        index_type = "swei", scale = as.integer(sc), success = FALSE,
        message = "SWEI seasonal directory not found; run 8snow_drought_indices.R first"
      )
    ), use.names = TRUE)
  }
}

if (dir.exists(SSPI_SEAS_DIR)) {
  cat("\n── SSPI (Standardized Snow Precipitation Index) ──\n")
  run_index_safe(
    "sspi",
    SSPI_SCALES,
    find_seasonal_nc_files,
    SSPI_SEAS_DIR
  )
} else {
  cat(sprintf("\n── SSPI skipped: seasonal directory not found: %s ──\n", SSPI_SEAS_DIR))
  for (sc in SSPI_SCALES) {
    RUN_STATUS <- data.table::rbindlist(list(
      RUN_STATUS,
      data.table::data.table(
        index_type = "sspi", scale = as.integer(sc), success = FALSE,
        message = "SSPI seasonal directory not found; run 8snow_drought_indices.R first"
      )
    ), use.names = TRUE)
  }
}


# ==============================================================================
#   17. FINAL RUN QA
# ==============================================================================

expected_runs <- data.table::rbindlist(
  list(
    data.table::data.table(index_type = "spi", scale = as.integer(SPI_SCALES)),
    data.table::data.table(index_type = "spei", scale = as.integer(SPEI_SCALES)),
    data.table::data.table(index_type = "spei_era5", scale = as.integer(SPEI_ERA5_SCALES)),
    data.table::data.table(index_type = "swei", scale = as.integer(SWEI_SCALES)),
    data.table::data.table(index_type = "sspi", scale = as.integer(SSPI_SCALES))
  ),
  use.names = TRUE
)

expected_runs[
  ,
  index_type :=
    tolower(index_type)
]

expected_runs[
  ,
  scale :=
    as.integer(scale)
]

run_status_path <-
  file.path(
    TREND_DIR,
    "script6_run_status.csv"
  )

data.table::fwrite(
  RUN_STATUS,
  run_status_path
)

missing_runs <-
  expected_runs[
    !RUN_STATUS,
    on = .(
      index_type,
      scale
    )
  ]

failed_runs <-
  RUN_STATUS[
    success == FALSE
  ]


cat(
  "\n══════════════════════════════════════════════════════\n",
  "SCRIPT 6 EXECUTION QA\n",
  "══════════════════════════════════════════════════════\n",
  sep = ""
)

cat(
  sprintf(
    "  Expected index × scale combinations: %d\n",
    nrow(expected_runs)
  )
)

cat(
  sprintf(
    "  Successful combinations: %d\n",
    sum(
      RUN_STATUS$success
    )
  )
)

cat(
  sprintf(
    "  Failed combinations: %d\n",
    nrow(failed_runs)
  )
)

cat(
  sprintf(
    "  Missing combinations: %d\n",
    nrow(missing_runs)
  )
)

cat(
  sprintf(
    "  Run-status file: %s\n",
    run_status_path
  )
)


if (
  nrow(failed_runs)
) {
  
  print(
    failed_runs
  )
}


if (
  nrow(missing_runs)
) {
  
  print(
    missing_runs
  )
}


if (
  nrow(failed_runs) ||
  nrow(missing_runs) ||
  nrow(RUN_STATUS) !=
  nrow(expected_runs)
) {
  
  stop(
    paste0(
      "SCRIPT 6 FAILED QA: not all ", nrow(expected_runs),
      " expected index × scale combinations completed successfully. Inspect ",
      run_status_path,
      " before using any trend outputs."
    )
  )
}


cat(
  sprintf("  ✓ All %d expected index × scale combinations completed successfully.\n", nrow(expected_runs))
)


# ==============================================================================
#   18. FINAL SUMMARY
# ==============================================================================

total_elapsed <-
  difftime(
    Sys.time(),
    total_start,
    units = "mins"
  )


cat(
  "\n══════════════════════════════════════════════════════\n",
  "DUAL DROUGHT EVENT METHOD SUMMARY\n",
  "══════════════════════════════════════════════════════\n",
  sep = ""
)

cat(
  "Spatial framework:\n",
  "  Native grid from each actual index NetCDF product.\n",
  "  No canonical grid RDS.\n",
  "  No stored cross-pipeline basin mask.\n",
  "  Basin polygon rasterized directly on the native index grid.\n",
  "  Fractional basin coverage × native cell area used for weighting.\n\n",
  sep = ""
)

cat(
  "Method 1 (S&W) — Sheffield & Wood (2008):\n"
)

cat(
  sprintf(
    "  Onset:        index < %.2f\n",
    DROUGHT_ONSET
  )
)

cat(
  sprintf(
    "  Termination:  index >= %.2f\n",
    DROUGHT_ONSET
  )
)

cat(
  "  Min duration: none\n"
)

cat(
  sprintf(
    "  Intensity:    mean(onset_thr - x)\n"
  )
)

cat(
  "  Duration classes:\n",
  "    D1-2, D3-6, D7-12, D13+\n\n",
  sep = ""
)

cat(
  "Method 2 (Hysteresis):\n"
)

cat(
  sprintf(
    "  Onset:        index < %.2f\n",
    DROUGHT_ONSET
  )
)

cat(
  sprintf(
    "  Termination:  index >= %.2f\n",
    DROUGHT_END
  )
)

cat(
  "  Intensity:    mean(end_thr - x)\n",
  "  Min duration:\n",
  "    scale 1  : >= 2 months\n",
  "    scale 3  : >= 3 months\n",
  "    scale 6  : >= 4 months\n",
  "    scale 12 : >= 6 months\n",
  "    scale 24 : >= 12 months\n",
  "    scale 36 : >= 12 months\n\n",
  sep = ""
)

cat(
  "Additional statistics:\n",
  "  VC-HR98 MK + Sen's slope\n",
  "  TFPW MK + Sen's slope\n",
  "  BH-FDR correction\n",
  "  PELT regime-shift detection\n",
  "  Wald-Wolfowitz runs test\n",
  "  MC AR(1) red-noise spectral peaks\n\n",
  sep = ""
)

cat(
  sprintf(
    "Total elapsed time: %.2f minutes\n",
    as.numeric(
      total_elapsed
    )
  )
)

cat(
  "══════════════════════════════════════════════════════\n"
)