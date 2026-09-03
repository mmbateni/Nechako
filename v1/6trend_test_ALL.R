# ==============================================================================
#   6trend_test_v2.R  ·  SPATIAL TREND ANALYSIS FOR ALL INDICES
# ==============================================================================
#   Computes per-pixel statistics for SPI, SPEI_PM, SPEI_ERA5:
#   • Mann-Kendall — Hamed & Rao (1998) variance-corrected (VC-HR98) + Sen's slope
#     Uses modifiedmk::mmkh() which adjusts Var(S) for lag-1…lag-(n-1) autocorrelation.
#     Tau is the standard concordance statistic; p-value is the H-R98 corrected one.
#     Sen's slope is extracted directly from mmkh() (avoids the O(n²) matrix method).
#   • Mann-Kendall (TFPW – Trend-Free Pre-Whitening, Yue et al. 2002)
#     Sen's slope now computed on original series, not pre-whitened.
#   • Drought event count / mean duration / max intensity — TWO METHODS:
#       Method 1 (S&W): single threshold DROUGHT_ONSET, no hysteresis, no min duration
#         Intensity = mean(DROUGHT_ONSET - x)  [deficit below -1.0]
#       Method 2 (Hyst): onset < DROUGHT_ONSET, termination >= DROUGHT_END,
#                        scale-specific minimum duration
#         Intensity = mean(DROUGHT_END - x) [deficit below 0.0, always >=0]
#   • PELT regime-shift year  [BUGFIX 7: now actually called in process_index]
#   • Wald-Wolfowitz runs test (drought event temporal clustering)
#     Binarised on hysteresis drought events
##   • n_spectral_peaks — Monte Carlo AR(1) red-noise significance test
#                        on periods >= 24 months(n_sim_spectral = 500)with a family-wise maximum-statistic threshold
#                        [BUGFIX] threshold is now a single family-wise
#                        max-statistic over the drought band, not a per-bin
#                        quantile — see comments at "LOCAL OVERRIDE" below.
#
# [NEW] One-time basin pixel-mask reconciliation check: on the first index ×
#   scale processed, compares this script's "valid pixel" grid against Script
#   4's ERA5-Land basin mask (clipped_template.rds), since the two pipelines
#   use different source rasters and their pixel counts (e.g. 865 vs. 856)
#   have never been formally reconciled. See "Basin pixel-mask reconciliation
#   check" below.
#
# OUTPUT: {TREND_DIR}/{index}_{scale:02d}_results.csv  (one per index × scale)
# Columns include tau_vc, p_value_vc (H-R98 corrected), p_fdr_vc (BH-FDR corrected),
# tau_tfpw, p_value_tfpw, p_fdr_tfpw (BH-FDR corrected),
# n_events, mean_duration, max_intensity            [Method 1: S&W]
# n_events_hyst, mean_duration_hyst, max_intensity_hyst  [Method 2: Hysteresis]
# n_events_D36, n_events_D712, n_events_D13p, n_D4p, mean_I, mean_S  [Method 1]
# n_events_D36_hyst, …_hyst                                           [Method 2]
# regime_shift_year, regime_shift_detected, n_changepoints,
# mean_before_shift, mean_after_shift, magnitude_shift,
# p_value_runs, clustering,
# n_spectral_peaks.
#
# p_fdr_* are computed by applying p.adjust(method="BH") across all valid basin
# pixels within each index × scale run. Use p_fdr_vc for significance mapping.
#
# METHOD COMPARISON NOTES:
#   Method 1 (Sheffield & Wood 2008): single threshold, no hysteresis,
#     no minimum duration filter. Every month below the threshold belongs to an event.
#     Intensity = mean deficit below the threshold.
#     Use for direct replication of S&W (2008) Table 4 statistics.
#   Method 2 (Hysteresis): onset <-1.0, termination >=0.0.
#     Once a drought begins (drops below -1.0), it persists until index rises
#     to 0.0 or above, preventing spurious terminations during brief recoveries.
#     Scale-specific minimum duration filters out noise at short scales.
#     Intensity = mean deficit below 0.0 (always positive by construction).
#     Columns suffixed _hyst. Use for climatologically meaningful event counts.
#
# SINGLE-FILE MODE [BUGFIX 1]:
#   The matrix layout is now detected automatically: single-file indices fill
#   columns sequentially (one column per month), whereas standard SPI/SPEI use
#   the interleaved layout (one file per calendar month, one layer per year).
#
# NOTE: Raw pixel time-series are NOT stored in the CSV (they live in the source
#       NetCDF files). This keeps CSVs lightweight and avoids ~6 M redundant values.
# Run BEFORE 4pr_pet_trends_visualization_v2.r
# ==============================================================================
source("DROUGHT_ANALYSIS_utils_v2.R")
utils_load_packages(c("terra", "data.table", "Kendall", "trend", "modifiedmk",
                      "parallel", "lubridate", "changepoint"))
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
dir.create(TREND_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Scale overrides ─────────────────────────────────────────────────────────
# DROUGHT_ANALYSIS_utils_v2.R defines SPI_SCALES = c(1, 2, 3, 6, 12) (5 scales).
# Scripts 1SPI_ERALand_v2.R and 3SPEI_ERALand_v2.R now compute 6 timescales (1,3,6,12,24,36)
# for multiyear drought persistence.  Override here so that Script 6 runs trend
# tests across every scale that actually exists on disk.
SPI_SCALES  <- c(1L, 3L, 6L, 12L, 24L, 36L)
SPEI_SCALES <- c(1L, 3L, 6L, 12L, 24L, 36L)

# ── ERA5-Land PEV SPEI directories and scales ────────────────────────────────
# SPEI_ERA5 NetCDF files are written by 3SPEI_ERALand_v2.R into
# spei_results_seasonal_era5land/ alongside the PM SPEI seasonal files.  The index_type
# "spei_era5" is used consistently in downstream scripts (7, 8, 10a).
SPEI_ERA5_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal_era5land")
SPEI_ERA5_SCALES   <- c(1L, 3L, 6L, 12L, 24L, 36L)

#   Utility operators
# --------------------------------------------------------------------------------

####################################################################################
# PARALLEL BACK-END (WINDOWS-COMPATIBLE)
####################################################################################
is_windows <- .Platform$OS.type == "windows"
N_CORES    <- max(1L, parallel::detectCores() - 1L)
cl         <- NULL
if (N_CORES > 1) {
  if (is_windows) {
    cl <- parallel::makeCluster(N_CORES)
    on.exit(parallel::stopCluster(cl), add = TRUE)  # guaranteed cleanup on error or normal exit
    parallel::clusterEvalQ(cl, {
      library(Kendall)
      library(trend)
      library(modifiedmk)
      source("DROUGHT_ANALYSIS_utils_v2.R") 
    })
    cat(sprintf("✓ Windows cluster: %d cores\n", N_CORES))
  } else {
    cat(sprintf("✓ Unix fork: %d cores\n", N_CORES))
  }
} else {
  cat("ℹ Single-core mode\n")
}

# ================================================================================
# SCALE-SPECIFIC MINIMUM DURATION HELPER
# ================================================================================
# Returns the minimum event duration (months) for Method 2 (Hysteresis) given
# the accumulation scale of the index.
# Rationale:
#   SPI/SPEI-1  is inherently noisy; require ≥ 2 months to avoid single-month spikes
#   SPI/SPEI-3  smoothed over a season; ≥ 3 months = one full accumulation window
#   SPI/SPEI-6  half-year smoothing; ≥ 4 months filters sub-seasonal noise
#   SPI/SPEI-12 annual smoothing; ≥ 6 months = half an accumulation window
get_min_duration <- function(scale,
                             min_duration_map     = list("1"  = 2L,
                                                         "3"  = 3L,
                                                         "6"  = 4L,
                                                         "12" = 6L,
                                                         "24" = 12L,
                                                         "36" =12L),
                             default_min_duration = 3L) {
  key <- as.character(as.integer(scale))
  val <- min_duration_map[[key]]
  if (!is.null(val)) as.integer(val) else as.integer(default_min_duration)
}

# ================================================================================
# VECTORISED STATISTICS FUNCTIONS
# ─── NOTE: vectorized_mann_kendall(), vectorized_mann_kendall_tfpw(), and
#          vectorized_regime_shift_pelt() are defined in DROUGHT_ANALYSIS_utils_v2.R
#          (Section J7).
#          vectorized_spectral_peaks() is LOCALLY OVERRIDDEN below (Section 6)
#          with a Monte Carlo red-noise implementation (n_sim_spectral = 500)
#          to match the null hypothesis used in Script 4.  The utils.R version
#          (AR(1) red-noise) is shadowed and not called by this script.
#          vectorized_event_stats() and vectorized_runs_test() remain
#          here because they are drought-index-specific.
# ================================================================================
#   3. Drought event characteristics — unified for both methods
# --------------------------------------------------------------------------------
# Returns n_events, mean_duration, max_intensity per pixel.
#
# METHOD SELECTION via onset_thr / end_thr / min_duration:
#   Method 1 (Sheffield & Wood 2008): end_thr = onset_thr, min_duration = 1L
#     Entry : x < onset_thr  (single threshold)
#     Exit  : x >= onset_thr (same threshold — no hysteresis)
#     Intensity: mean(onset_thr - x)  [deficit below the threshold; always positive]
#
#   Method 2 (Hysteresis): end_thr > onset_thr, min_duration from get_min_duration()
#     Entry : x < onset_thr  (= -1.0)
#     Exit  : x >= end_thr   (= 0.0)
#     [BUGFIX 4] Intensity: mean(end_thr - x) = mean(-x)
#       Using end_thr ensures intensity is always non-negative.
#       Under the old code (onset_thr - x), months in the recovery zone
#       (-1.0 <= x < 0.0) produced negative intensity contributions, which is
#       physically meaningless and inflated mean_I toward zero or below.
#
# max_intensity is always min(x) during the event (most negative value),
# unchanged for both methods.
vectorized_event_stats <- function(ts_matrix,
                                   onset_thr    = DROUGHT_ONSET,
                                   end_thr      = onset_thr,
                                   min_duration = 1L) {
  if (end_thr < onset_thr)
    warning(sprintf(
      "vectorized_event_stats: end_thr (%.2f) < onset_thr (%.2f). For hysteresis, end_thr should be >= onset_thr. Proceeding anyway.",
      end_thr, onset_thr))  
  n_pix <- nrow(ts_matrix)
  n_events      <- integer(n_pix)
  mean_duration <- numeric(n_pix)
  max_intensity <- numeric(n_pix)
  n_events_D36  <- integer(n_pix)
  n_events_D712 <- integer(n_pix)
  n_events_D13p <- integer(n_pix)
  n_D4p         <- integer(n_pix)
  mean_I        <- numeric(n_pix)
  mean_S        <- numeric(n_pix)
  
  for (i in seq_len(n_pix)) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) == 0L) next
    
    in_d  <- FALSE
    s_idx <- NA_integer_
    durs  <- integer(0)
    ints  <- numeric(0)
    I_evs <- numeric(0)
    S_evs <- numeric(0)
    
    for (j in seq_along(x)) {
      v <- x[j]
      
      if (!in_d && v < onset_thr) {
        in_d  <- TRUE
        s_idx <- j
        
      } else if (in_d && v >= end_thr) {
        dur <- j - s_idx
        if (dur >= min_duration) {
          seg  <- x[s_idx:(j - 1L)]
          durs <- c(durs, dur)
          ints <- c(ints, min(seg, na.rm = TRUE))
          # BUGFIX 4: deficit measured from end_thr (not onset_thr)
          # → always non-negative because the loop only reaches here when
          #   the current month first satisfies v >= end_thr, so all prior
          #   months in seg had v < end_thr.
          I_ev  <- mean(end_thr - seg, na.rm = TRUE)
          I_evs <- c(I_evs, I_ev)
          S_evs <- c(S_evs, I_ev * dur)
        }
        in_d  <- FALSE
        s_idx <- NA_integer_
      }
    }
    
    # Close any event still open at end of record
    if (in_d && !is.na(s_idx)) {
      dur <- length(x) - s_idx + 1L
      if (dur >= min_duration) {
        seg  <- x[s_idx:length(x)]
        durs <- c(durs, dur)
        ints <- c(ints, min(seg, na.rm = TRUE))
        # BUGFIX 4 (same fix for end-of-record event)
        I_ev  <- mean(end_thr - seg, na.rm = TRUE)
        I_evs <- c(I_evs, I_ev)
        S_evs <- c(S_evs, I_ev * dur)
      }
    }
    
    if (length(durs) == 0L) {
      n_events[i]      <- 0L
      mean_duration[i] <- 0
      max_intensity[i] <- 0
      n_events_D36[i]  <- 0L
      n_events_D712[i] <- 0L
      n_events_D13p[i] <- 0L
      n_D4p[i]         <- 0L
      mean_I[i]        <- NA_real_
      mean_S[i]        <- NA_real_
    } else {
      n_events[i]      <- length(durs)
      mean_duration[i] <- mean(durs)
      max_intensity[i] <- min(ints)
      n_events_D36[i]  <- sum(durs >= 3L  & durs <= 6L)   # Short: 3–6 months  (README §11)
      n_events_D712[i] <- sum(durs >= 7L  & durs <= 12L)  # Medium: 7–12 months
      n_events_D13p[i] <- sum(durs >= 13L)                 # Long: ≥ 13 months   (README §11)
      n_D4p[i]         <- sum(durs >= 4L)
      mean_I[i]        <- mean(I_evs, na.rm = TRUE)
      mean_S[i]        <- mean(S_evs, na.rm = TRUE)
    }
  }
  
  data.table::data.table(
    n_events      = n_events,
    mean_duration = mean_duration,
    max_intensity = max_intensity,
    n_events_D36  = n_events_D36,
    n_events_D712 = n_events_D712,
    n_events_D13p = n_events_D13p,
    n_D4p         = n_D4p,
    mean_I        = mean_I,
    mean_S        = mean_S
  )
}

# --------------------------------------------------------------------------------
#   4. PELT regime-shift year
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
#   5. Wald-Wolfowitz runs test — drought event temporal clustering
# --------------------------------------------------------------------------------
# BUGFIX 6: The binary sequence is now built from DROUGHT EVENTS defined by the
#   two-threshold hysteresis method (onset_thr / end_thr / min_duration), exactly
#   matching vectorized_event_stats() Method 2.  Each month that belongs to a
#   qualifying drought event is coded 1; all other months are coded 0.
#
#   The previous implementation binarised the raw index at threshold = 0 (wet vs.
#   dry), which tests general wet/dry clustering rather than whether *drought*
#   events — as the script defines them — cluster in time.  That approach also
#   ignores hysteresis and minimum-duration requirements entirely.
#
#   A "clustered" result (fewer runs than expected) means drought events tend to
#   arrive in multi-event groups separated by long dry spells; "dispersed" means
#   they are more regularly spaced than random.
#
# Arguments:
#   ts_matrix    — pixels × months matrix
#   onset_thr    — drought onset threshold  (default: DROUGHT_ONSET = -1.0)
#   end_thr      — drought termination threshold (default: DROUGHT_END = 0.0)
#   min_duration — minimum qualifying event length (months)
vectorized_runs_test <- function(ts_matrix,
                                 onset_thr    = DROUGHT_ONSET,
                                 end_thr      = onset_thr,
                                 min_duration = 1L) {
  n_pix <- nrow(ts_matrix)
  res   <- data.table::data.table(
    p_value_runs  = rep(NA_real_,      n_pix),
    clustering    = rep(NA_character_, n_pix),
    filtered_runs = rep(TRUE,          n_pix)
  )
  
  for (i in seq_len(n_pix)) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 10L) next
    
    # ── Build drought binary indicator using hysteresis + min_duration ──────────
    # Mirrors the state machine in vectorized_event_stats exactly.
    b     <- integer(length(x))
    in_d  <- FALSE
    s_idx <- NA_integer_
    
    for (j in seq_along(x)) {
      v <- x[j]
      if (!in_d && v < onset_thr) {
        in_d  <- TRUE
        s_idx <- j
      } else if (in_d && v >= end_thr) {
        dur <- j - s_idx
        if (dur >= min_duration) b[s_idx:(j - 1L)] <- 1L
        in_d  <- FALSE
        s_idx <- NA_integer_
      }
    }
    # Close open event at end of record
    if (in_d && !is.na(s_idx)) {
      dur <- length(x) - s_idx + 1L
      if (dur >= min_duration) b[s_idx:length(x)] <- 1L
    }
    
    # ── Wald-Wolfowitz runs test on drought binary series ───────────────────────
    n  <- length(b)
    n1 <- sum(b)        # drought months
    n0 <- n - n1        # non-drought months
    if (n1 == 0L || n0 == 0L) next
    
    nr <- length(rle(b)$lengths)
    er <- 2.0 * n0 * n1 / n + 1.0
    vr <- 2.0 * n0 * n1 * (2.0 * n0 * n1 - n) / (n^2 * (n - 1.0))
    if (vr <= 0) next
    
    z    <- (nr - er) / sqrt(vr)
    pv   <- 2.0 * pnorm(-abs(z))
    clus <- if (nr < er) "clustered" else if (nr > er) "dispersed" else "random"
    
    res[i, `:=`(p_value_runs  = pv,
                clustering    = clus,
                filtered_runs = FALSE)]
  }
  res
}

# --------------------------------------------------------------------------------
#   6. Spectral peak detection — drought-band frequencies
# --------------------------------------------------------------------------------
# LOCAL OVERRIDE: vectorized_spectral_peaks() is redefined here to use a Monte
#   Carlo white-noise null hypothesis (n_sim_spectral = 500 simulations), matching
#   the envelope used in Script 4 (mk_tfpw_spectral_for_series / conf_env).
#   This replaces the AR(1) red-noise test in DROUGHT_ANALYSIS_utils_v2.R, which is
#   more conservative when the series has positive autocorrelation (as drought
#   indices typically do), making Script 4 and Script 6 spectral results directly
#   comparable.
#
# METHOD:
#   For each pixel, the monthly time series is demeaned and linearly detrended,
#   then a smoothed periodogram (modified Daniell smoother, spans = c(3,5)) is
#   computed.  The white-noise significance envelope is precomputed ONCE per
#   unique series length using n_sim_spectral = 500 iid N(0,1) realisations,
#   each passed through the same smoother.
#
#   [BUGFIX] Previously the envelope was a PER-BIN sig_level quantile (one
#   threshold per frequency bin, `apply(sim_specs, 1L, quantile, ...)`), and a
#   pixel was flagged "significant" if ANY bin in the drought band exceeded its
#   own bin-specific 95th percentile. With dozens of bins in the >=24-month
#   band, testing each one independently at alpha = 0.05 pushes the effective
#   per-pixel false-positive rate far above 5% (family-wise error inflation) --
#   consistent with the near-universal (865/865, 100%) "significant peak"
#   detection observed across 17 of 18 index/scale combinations.
#
#   Fix: use a single GLOBAL MAXIMUM-STATISTIC threshold per series length,
#   i.e. for each of the n_sim_spectral white-noise realisations take the max
#   smoothed-periodogram ordinate across ONLY the drought-band bins, then take
#   the sig_level quantile of that distribution of maxima. This gives one
#   family-wise-corrected threshold that is applied uniformly to every bin in
#   the band, so the true per-pixel false-positive rate stays near alpha
#   regardless of how many bins fall inside the band (this mirrors the
#   max-statistic approach already used by Method 2 in
#   DROUGHT_ANALYSIS_utils_v2.R::.spectral_single(), which this override is
#   meant to match).
#
#   Only periods >= min_period_months (default 24 months / 2 years) are examined,
#   excluding weather noise and focusing on multi-year drought cycles
#   (e.g. ENSO 24-84 mo, PDO 120-360 mo). The band restriction is now applied
#   BEFORE computing the null threshold (not just when flagging bins), so
#   bins outside the band can no longer influence -- or be flagged by -- the
#   test at all.
#
# Requires ≥ 48 valid observations (4 years). Uses only base-R stats::spectrum().
#
# Arguments:
#   ts_matrix         — pixels × months matrix
#   min_period_months — shortest period of interest (default 24)
#   sig_level         — Monte Carlo quantile threshold (default 0.95)
#   n_sim_spectral    — number of white-noise simulations (default 500)
# ================================================================================
n_sim_spectral <- 500L   # matches Script 4 parameter

vectorized_spectral_peaks <- function(ts_matrix,
                                      min_period_months = 24L,
                                      sig_level         = 0.95,
                                      n_sim_spectral    = get("n_sim_spectral",
                                                              envir = parent.env(environment()))) {
  n_pix <- nrow(ts_matrix)
  out   <- integer(n_pix)
  
  # ── [BUGFIX 2] Null switched from WHITE noise to AR(1) RED noise ──────────────
  # After fixing the band restriction and family-wise correction (see comments
  # above), scales 1 were no longer near-universal, but scales 3/6/12/24/36
  # remained at ~865/865 (100%) for EVERY index. Root cause: SPI-n/SPEI-n for
  # n >= 3 are built from an n-month rolling accumulation, which imparts
  # substantial lag-1 autocorrelation into the series by construction --
  # independent of any genuine drought periodicity. A WHITE-noise null (iid,
  # rho1 = 0) is the wrong null for an accumulated/smoothed series: almost any
  # such series will show elevated low-frequency power relative to white noise
  # simply because of the accumulation, not because of real periodicity. This
  # is exactly the situation the AR(1) red-noise null (Method 1, in
  # DROUGHT_ANALYSIS_utils_v2.R) exists to handle.
  #
  # Fix: estimate each pixel's own lag-1 autocorrelation (rho1) on the
  # detrended series and generate the Monte Carlo null from an AR(1) process
  # with that same rho1, instead of iid N(0,1). This keeps the family-wise,
  # band-restricted max-statistic correction from the previous fix, but now
  # tests against "is this more periodic than a red-noise process with the
  # SAME persistence as this pixel's own series" -- the appropriate null for
  # accumulated indices -- rather than against plain white noise.
  #
  # rho1 is binned (rounded to the nearest 0.05) before caching, since the
  # exact-per-pixel value would defeat caching entirely; this keeps the
  # threshold specific enough to track each pixel's actual persistence while
  # still reusing simulations across pixels with similar rho1.
  rn_env_cache <- list()
  
  get_red_noise_threshold <- function(n_len, band, rho1) {
    rho1_bin <- round(rho1 / 0.05) * 0.05
    rho1_bin <- max(-0.95, min(0.95, rho1_bin))
    key <- paste0(n_len, "_r", sprintf("%.2f", rho1_bin))
    if (!is.null(rn_env_cache[[key]])) return(rn_env_cache[[key]])
    
    set.seed(42L)   # reproducible across pixels and index × scale runs
    band_maxima <- vapply(seq_len(n_sim_spectral), function(s) {
      innov <- stats::rnorm(n_len, sd = sqrt(max(1e-8, 1 - rho1_bin^2)))
      rn    <- stats::filter(innov, filter = rho1_bin, method = "recursive")
      rn    <- as.numeric(rn)
      sp    <- stats::spectrum(rn, spans = c(3L, 5L), taper = 0,
                               detrend = FALSE, plot = FALSE)
      # Global max WITHIN the drought band only.
      max(sp$spec[band], na.rm = TRUE)
    }, numeric(1))
    
    thr <- stats::quantile(band_maxima, probs = sig_level, names = FALSE)
    rn_env_cache[[key]] <<- thr
    thr
  }
  
  for (i in seq_len(n_pix)) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 48L) next
    
    # Demean and linearly detrend
    x  <- x - mean(x)
    tt <- seq_along(x)
    x  <- residuals(lm(x ~ tt))
    
    pv <- var(x)
    if (!is.finite(pv) || pv <= 0) next
    
    # Lag-1 autocorrelation of THIS pixel's detrended series -- drives the
    # red-noise null so accumulated (high-persistence) series are compared
    # against an appropriately persistent null, not white noise.
    n_x  <- length(x)
    rho1 <- tryCatch(stats::cor(x[-n_x], x[-1], use = "complete.obs"),
                     error = function(e) 0.0)
    if (is.na(rho1) || !is.finite(rho1)) rho1 <- 0.0
    rho1 <- max(-0.95, min(0.95, rho1))
    
    sp  <- tryCatch(
      stats::spectrum(x, spans = c(3L, 5L), taper = 0,
                      detrend = FALSE, plot = FALSE),
      error = function(e) NULL
    )
    if (is.null(sp)) next
    
    # Period in months for each frequency bin (freq is cycles per sample)
    periods <- 1.0 / sp$freq
    
    # Restrict to drought-band periods
    band <- periods >= min_period_months
    if (!any(band)) next
    
    # Global (family-wise) red-noise threshold, scaled to this pixel's variance
    thr_unit <- get_red_noise_threshold(length(x), band, rho1)
    threshold <- thr_unit * pv
    
    # Flag bins that exceed the SINGLE global threshold and lie in the drought band
    sig_bins <- band & (sp$spec > threshold)
    
    # Count contiguous runs of significant bins as individual peaks
    if (any(sig_bins)) {
      r       <- rle(sig_bins)
      out[i]  <- sum(r$values)
    }
  }
  
  data.table::data.table(n_spectral_peaks = out)
}
#
# Arguments:
#   ts_matrix          — pixels × months matrix
#   min_period_months  — shortest period of interest (default 24)
#   sig_level          — chi-squared significance threshold (default 0.95)
# ================================================================================
#   FUNCTION 7 — Basin-level % area in drought time series + MK trend
# ================================================================================
compute_basin_extent_and_class_trends <- function(ts_matrix, n_years, index_label,
                                                  out_dir,
                                                  onset_thr    = DROUGHT_ONSET,
                                                  end_thr      = DROUGHT_END,
                                                  start_year   = 1950L,
                                                  win_years    = 30L,
                                                  area_weights = NULL,
                                                  basin_avg    = NULL) {
  use_weights <- !is.null(area_weights) &&
    length(area_weights) == nrow(ts_matrix) &&
    all(is.finite(area_weights)) &&
    all(area_weights > 0)
  if (!is.null(area_weights) && !use_weights)
    warning("compute_basin_extent_and_class_trends: area_weights ignored (length mismatch or invalid values); falling back to equal weighting")
  
  utils_load_packages(c("ggplot2", "Kendall", "dplyr", "patchwork"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  n_months  <- n_years * 12L
  month_seq <- seq_len(n_months)
  yr_seq    <- start_year + floor((month_seq - 1L) / 12L)
  mo_seq    <- ((month_seq - 1L) %% 12L) + 1L
  date_seq  <- as.Date(sprintf("%04d-%02d-01", yr_seq, mo_seq))
  
  n_pix    <- nrow(ts_matrix)
  n_col    <- min(n_months, ncol(ts_matrix))
  pct_area <- numeric(n_col)
  for (t in seq_len(n_col)) {
    col_vals <- ts_matrix[, t]
    ok       <- !is.na(col_vals)
    if (!any(ok)) { pct_area[t] <- NA; next }
    if (use_weights) {
      w_ok        <- area_weights[ok]
      drought_ok  <- ok & col_vals < onset_thr
      pct_area[t] <- 100.0 * sum(area_weights[drought_ok], na.rm = TRUE) /
        sum(w_ok, na.rm = TRUE)
    } else {
      pct_area[t] <- 100.0 * sum(ok & col_vals < onset_thr) / sum(ok)
    }
  }
  extent_ts <- data.frame(
    date      = date_seq[seq_len(n_col)],
    year      = yr_seq[seq_len(n_col)],
    month     = mo_seq[seq_len(n_col)],
    pct_area  = pct_area,
    row.names = NULL
  )
  write.csv(extent_ts,
            file.path(out_dir, sprintf("%s_basin_extent_monthly.csv", index_label)),
            row.names = FALSE)
  
  extent_annual <- extent_ts |>
    dplyr::group_by(year) |>
    dplyr::summarise(mean_pct = mean(pct_area, na.rm = TRUE), .groups = "drop")
  x_clean <- extent_annual$mean_pct[!is.na(extent_annual$mean_pct)]
  mk_ext  <- list(tau = NA_real_, p_value = NA_real_, sens_slope = NA_real_)
  if (length(x_clean) >= 10L) {
    tau_val <- as.numeric(Kendall::MannKendall(x_clean)$tau)
    hr_r    <- tryCatch(modifiedmk::mmkh(x_clean), error = function(e) NULL)
    if (!is.null(hr_r)) {
      mk_ext <- list(tau        = tau_val,
                     p_value    = as.numeric(hr_r["new P-value"]),
                     sens_slope = as.numeric(hr_r["Sen's slope"]))
    } else {
      mk_r   <- Kendall::MannKendall(x_clean)
      mk_ext <- list(tau        = tau_val,
                     p_value    = as.numeric(mk_r$sl),
                     sens_slope = as.numeric(trend::sens.slope(x_clean)$estimates))
    }
  }
  mk_ext_df <- data.frame(
    index       = index_label,
    variable    = "mean_annual_pct_area_drought",
    tau         = mk_ext$tau,
    p_value     = mk_ext$p_value,
    sens_slope  = mk_ext$sens_slope,
    direction   = ifelse(!is.na(mk_ext$tau),
                         ifelse(mk_ext$tau > 0, "Increasing", "Decreasing"), NA),
    significant = ifelse(!is.na(mk_ext$p_value), mk_ext$p_value < 0.05, NA),
    stringsAsFactors = FALSE
  )
  write.csv(mk_ext_df,
            file.path(out_dir, sprintf("%s_basin_extent_MK.csv", index_label)),
            row.names = FALSE)
  cat(sprintf("  [Extent MK | %s] tau=%.3f  p=%.3f  slope=%.4f %%/yr  %s\n",
              index_label,
              mk_ext$tau %||% NA, mk_ext$p_value %||% NA, mk_ext$sens_slope %||% NA,
              ifelse(!is.na(mk_ext$p_value) && mk_ext$p_value < 0.05, "sig", "")))
  
  if (is.null(basin_avg)) {
    if (use_weights) {
      basin_avg <- vapply(seq_len(n_col), function(t) {
        v  <- ts_matrix[, t]
        ok <- !is.na(v) & is.finite(v)
        if (!any(ok)) return(NA_real_)
        sum(v[ok] * area_weights[ok]) / sum(area_weights[ok])
      }, numeric(1L))
    } else {
      basin_avg <- colMeans(ts_matrix[, seq_len(n_col)], na.rm = TRUE)
    }
  }
  
  classify_dur <- function(d)
    dplyr::case_when(d >= 3L  & d <= 6L  ~ "D3-6 (Short-term)",
                     d >= 7L  & d <= 12L ~ "D7-12 (Medium-term)",
                     d >= 13L            ~ "D13+ (Long-term)",
                     TRUE                ~ "D1-3 (Sub-threshold)")
  
  vals       <- basin_avg
  in_d       <- FALSE
  s_idx      <- NA_integer_
  event_list <- list()
  for (j in seq_along(vals)) {
    v <- vals[j]
    if (!in_d && !is.na(v) && v < onset_thr) {
      in_d <- TRUE; s_idx <- j
    } else if (in_d && !is.na(v) && v >= end_thr) {
      dur <- j - s_idx
      event_list[[length(event_list) + 1L]] <- data.frame(
        start_date      = date_seq[s_idx],
        end_date        = date_seq[j - 1L],
        start_year      = yr_seq[s_idx],
        duration_months = dur,
        duration_class  = classify_dur(dur),
        mean_severity   = mean(abs(pmin(vals[s_idx:(j - 1L)], 0)), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      in_d <- FALSE; s_idx <- NA_integer_
    }
  }
  if (in_d && !is.na(s_idx)) {
    dur <- n_col - s_idx + 1L
    event_list[[length(event_list) + 1L]] <- data.frame(
      start_date      = date_seq[s_idx],
      end_date        = date_seq[n_col],
      start_year      = yr_seq[s_idx],
      duration_months = dur,
      duration_class  = classify_dur(dur),
      mean_severity   = mean(abs(pmin(vals[s_idx:n_col], 0)), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  # Figure: basin extent time series with MK annotation
  extent_annual$trend_line <- NA_real_
  if (!is.na(mk_ext$sens_slope) && !is.na(mk_ext$tau)) {
    mid  <- median(seq_along(x_clean))
    ymed <- median(x_clean, na.rm = TRUE)
    extent_annual$trend_line <-
      ymed + mk_ext$sens_slope * (seq_along(extent_annual$mean_pct) - mid)
  }
  sig_label <- if (!is.na(mk_ext$p_value) && mk_ext$p_value < 0.05)
    sprintf("MK: tau=%.3f, p=%.3f* (significant)", mk_ext$tau, mk_ext$p_value)
  else if (!is.na(mk_ext$p_value))
    sprintf("MK: tau=%.3f, p=%.3f (not significant)", mk_ext$tau, mk_ext$p_value)
  else
    "MK: insufficient data"
  
  x_min       <- min(extent_annual$year, na.rm = TRUE)
  x_max       <- max(extent_annual$year, na.rm = TRUE)
  first_break <- ceiling(x_min / 10) * 10
  last_break  <- floor(x_max / 10) * 10
  x_breaks    <- seq(first_break, last_break, by = 10)
  x_breaks    <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  if (length(x_breaks) == 0) x_breaks <- c(x_min, x_max)
  if (length(x_breaks) > 15) {
    step_size <- max(10, ceiling((x_max - x_min) / 10))
    x_breaks  <- seq(first_break, last_break, by = step_size)
    x_breaks  <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  }
  
  p_ext <- ggplot2::ggplot(extent_annual, ggplot2::aes(x = year, y = mean_pct)) +
    ggplot2::geom_col(fill = "#c0392b", alpha = 0.5, width = 0.8) +
    ggplot2::geom_line(colour = "#2c3e50", linewidth = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = trend_line), colour = "navy",
                       linewidth = 1.1, linetype = "dashed", na.rm = TRUE) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      limits = c(x_min - 0.5, x_max + 0.5),
      expand = ggplot2::expansion(mult = c(0.01, 0.01)),
      labels = as.integer(x_breaks)
    ) +
    ggplot2::labs(
      title    = sprintf("%s — Annual Mean Basin-Wide Drought Extent", index_label),
      subtitle = sig_label,
      x = "Year", y = "% Basin Area in Drought"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s_basin_extent_timeseries.png", index_label)),
    p_ext, width = 12, height = 5, dpi = 150
  )
  cat(sprintf("  ✓ Saved: %s_basin_extent_timeseries.png\n", index_label))
  
  if (length(event_list) == 0) {
    cat(sprintf("  [Class Trends | %s] No events detected — skipping class MK\n",
                index_label))
    return(invisible(list(extent_ts   = extent_ts,
                          mk_extent   = mk_ext_df,
                          events_df   = NULL,
                          mk_by_class = NULL)))
  }
  
  events_df       <- do.call(rbind, event_list)
  events_df$index <- index_label
  write.csv(events_df,
            file.path(out_dir, sprintf("%s_basin_event_catalog.csv", index_label)),
            row.names = FALSE)
  
  invisible(list(extent_ts   = extent_ts,
                 mk_extent   = mk_ext_df,
                 events_df   = events_df))
}

# ================================================================================
#   MAIN PROCESSING FUNCTION
# ================================================================================
process_index <- function(index_type, scales, find_fn, seas_dir) {
  for (sc in scales) {
    cat(sprintf("\n╔══ %s-%02d ══╗\n", toupper(index_type), sc))
    
    min_dur <- get_min_duration(sc)
    cat(sprintf("  Method 2 (Hyst) min duration for scale-%02d: %d months\n",
                sc, min_dur))
    
    # [1/5] Find monthly NC files
    cat("[1/5] Finding files...\n")
    files <- find_fn(seas_dir, index_type, sc)
    if (!length(files)) { cat("  ❌ No files found\n"); next }
    cat(sprintf("  Found %d monthly NC files\n", length(files)))
    
    # [2/5] Load raster stacks
    cat("[2/5] Loading rasters...\n")
    stacks <- lapply(files, function(f) {
      tryCatch(terra::rast(f),
               error = function(e) { cat("  ⚠ Failed: ", basename(f), "\n"); NULL })
    })
    stacks <- stacks[!sapply(stacks, is.null)]
    if (!length(stacks)) { cat("  ❌ All rasters failed to load\n"); next }
    
    r_tmpl  <- stacks[[1]][[1]]
    xy      <- terra::xyFromCell(r_tmpl, seq_len(terra::ncell(r_tmpl)))
    lon_vec <- xy[, 1]
    lat_vec <- xy[, 2]
    n_pix   <- length(lon_vec)
    
    # [3/5] Build time-series matrix
    # ------------------------------------------------------------------------------
    # Read each calendar-month NetCDF ONCE as a complete value matrix.
    # This avoids 912 individual terra::values() calls on Windows.
    #
    # Standard SPI/SPEI layout:
    #   12 files = calendar months
    #   each file = one layer per year
    #   final matrix = chronological Jan-1950 ... Dec-2025
    # ------------------------------------------------------------------------------
    cat("[3/5] Building time-series matrix...\n")
    
    n_files         <- length(stacks)
    n_lyrs_per_file <- terra::nlyr(stacks[[1]])
    
    if (n_files == 1L && n_lyrs_per_file > 12L) {
      
      # --------------------------------------------------------------------------
      # Single-file mode:
      # one NetCDF contains the complete monthly record sequentially.
      # --------------------------------------------------------------------------
      
      n_months_total <- n_lyrs_per_file
      dates_all      <- extract_dates_from_nc(files[1], n_months_total)
      
      years   <- sort(unique(as.integer(format(dates_all, "%Y"))))
      n_years <- length(years)
      
      cat(sprintf(
        "  Single-file mode: %d pixels × %d months (%d years)\n",
        n_pix, n_months_total, n_years
      ))
      
      ts_mat <- terra::values(
        stacks[[1]],
        mat = TRUE,
        na.rm = FALSE
      )
      
      ts_mat <- as.matrix(ts_mat)
      
      expected_dim <- c(
        as.integer(n_pix),
        as.integer(n_months_total)
      )
      
      actual_dim <- dim(ts_mat)
      
      if (length(actual_dim) != 2L ||
          any(as.integer(actual_dim) != expected_dim)) {
        stop(
          "Unexpected single-file value matrix dimensions: ",
          paste(dim(ts_mat), collapse = " × "),
          " ; expected ",
          n_pix, " × ", n_months_total
        )
      }
      
      dates_ts <- dates_all
      
    } else {
      
      # --------------------------------------------------------------------------
      # Standard 12-file mode:
      # each file contains one calendar month across all years.
      # Read each file once, then interleave the columns.
      # --------------------------------------------------------------------------
      
      if (n_files != 12L) {
        stop(
          "Expected 12 calendar-month NetCDF files in standard mode; found ",
          n_files, "."
        )
      }
      
      dates_m1 <- extract_dates_from_nc(
        files[1],
        n_lyrs_per_file
      )
      
      years   <- as.integer(format(dates_m1, "%Y"))
      n_years <- length(years)
      
      cat(sprintf(
        "  Grid: %d pixels × %d years × 12 months\n",
        n_pix, n_years
      ))
      
      # Read each monthly NetCDF once.
      monthly_values <- vector("list", n_files)
      
      for (m in seq_len(n_files)) {
        
        cat(sprintf(
          "    Reading month %02d/12: %s\n",
          m,
          basename(files[m])
        ))
        
        monthly_values[[m]] <- terra::values(
          stacks[[m]],
          mat = TRUE,
          na.rm = FALSE
        )
        
        monthly_values[[m]] <- as.matrix(
          monthly_values[[m]]
        )
        
        expected_dim <- c(
          as.integer(n_pix),
          as.integer(min(terra::nlyr(stacks[[m]]), n_years))
        )
        
        actual_dim <- dim(monthly_values[[m]])
        
        if (length(actual_dim) != 2L ||
            any(as.integer(actual_dim) != expected_dim)) {
          stop(
            "Unexpected dimensions for month ",
            m,
            ": ",
            paste(actual_dim, collapse = " × "),
            " ; expected ",
            paste(expected_dim, collapse = " × ")
          )
        }
      }
      
      # Final chronological matrix:
      # columns = Jan-1950, Feb-1950, ..., Dec-2025.
      ts_mat <- matrix(
        NA_real_,
        nrow = n_pix,
        ncol = n_years * 12L
      )
      
      for (m in seq_len(12L)) {
        
        vals_m <- monthly_values[[m]]
        
        n_yrs_m <- min(
          ncol(vals_m),
          n_years
        )
        
        for (y in seq_len(n_yrs_m)) {
          
          col_idx <- (y - 1L) * 12L + m
          
          ts_mat[, col_idx] <- vals_m[, y]
        }
      }
      
      dates_ts <- as.Date(
        sprintf(
          "%04d-%02d-01",
          rep(years, each = 12L),
          rep(1:12, times = n_years)
        )
      )
      
      rm(monthly_values)
      gc(verbose = FALSE)
    }
    
    # Basic integrity checks before proceeding.
    if (!is.matrix(ts_mat)) {
      stop("ts_mat is not a matrix.")
    }
    
    if (nrow(ts_mat) != n_pix) {
      stop(
        "ts_mat row count mismatch: ",
        nrow(ts_mat),
        " vs expected ",
        n_pix
      )
    }
    
    if (ncol(ts_mat) != length(dates_ts)) {
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
      label = "6trend time-series matrix calendar",
      require_complete = TRUE
    )
    
    cat(sprintf(
      "  ✓ Time-series matrix built: %d pixels × %d months\n",
      nrow(ts_mat),
      ncol(ts_mat)
    ))
    
    cat(sprintf(
      "  ✓ Date range: %s to %s\n",
      min(dates_ts),
      max(dates_ts)
    ))
    
    finite_fraction <- mean(is.finite(ts_mat))
    
    cat(sprintf(
      "  ✓ Finite-value fraction: %.2f%%\n",
      100 * finite_fraction
    ))
    
    valid_pix <- rowSums(is.finite(ts_mat)) > 0
    
    # Enforce the canonical basin mask produced by 4pr_pet_trends_v2.r.
    # This is a hard geometry check, not a heuristic area comparison.
    # Enforce the canonical basin mask produced by 4pr_pet_trends_v2.r.
    #
    # IMPORTANT:
    # canonical_basin_mask.rds contains a SpatRaster in $template.
    # Deserializing that external-pointer object is unnecessary here and can
    # trigger a native terra/R session abort on Windows. Script 4 already stores
    # the complete geometry metadata as ordinary R vectors/scalars.
    #
    # Therefore Script 6 validates geometry using the stored metadata and uses
    # only the plain logical basin_mask vector.
    if (!file.exists(ERA5_CANONICAL_MASK_FILE)) {
      stop(
        "Canonical basin mask not found: ",
        ERA5_CANONICAL_MASK_FILE,
        "\nRun 4pr_pet_trends_v2.r after 3SPEI_ERALand_v2.R to create it."
      )
    }
    
    cat("  → Loading canonical basin mask metadata...\n")
    
    canonical_mask <- readRDS(ERA5_CANONICAL_MASK_FILE)
    
    if (!is.list(canonical_mask)) {
      stop("Invalid canonical_basin_mask.rds: expected a list.")
    }
    
    if (is.null(canonical_mask$basin_mask)) {
      stop("Invalid canonical_basin_mask.rds: missing basin_mask.")
    }
    
    canonical_mask_vec <- as.logical(canonical_mask$basin_mask)
    
    if (length(canonical_mask_vec) != terra::ncell(r_tmpl)) {
      stop(
        "Canonical basin mask length differs from SPEI raster cell count: ",
        length(canonical_mask_vec),
        " vs ",
        terra::ncell(r_tmpl)
      )
    }
    
    # --------------------------------------------------------------------------
    # Geometry validation using scalar metadata written by Script 4.
    # This avoids compareGeom() against a serialized SpatRaster.
    # --------------------------------------------------------------------------
    
    if (!is.null(canonical_mask$ncell)) {
      if (as.integer(terra::ncell(r_tmpl)) !=
          as.integer(canonical_mask$ncell)) {
        stop(
          "SPEI grid cell count does not match canonical Script-4 grid: ",
          terra::ncell(r_tmpl),
          " vs ",
          canonical_mask$ncell
        )
      }
    }
    
    if (!is.null(canonical_mask$resolution)) {
      r_res  <- as.numeric(terra::res(r_tmpl))
      can_res <- as.numeric(canonical_mask$resolution)
      
      if (length(r_res) != 2L ||
          length(can_res) != 2L ||
          any(!is.finite(r_res)) ||
          any(!is.finite(can_res)) ||
          any(abs(r_res - can_res) >
              1e-8 * pmax(1, abs(can_res)))) {
        stop(
          "SPEI grid resolution does not match canonical Script-4 grid: ",
          paste(signif(r_res, 10), collapse = " × "),
          " vs ",
          paste(signif(can_res, 10), collapse = " × ")
        )
      }
    }
    
    if (!is.null(canonical_mask$extent)) {
      r_ext <- c(
        xmin = terra::xmin(r_tmpl),
        xmax = terra::xmax(r_tmpl),
        ymin = terra::ymin(r_tmpl),
        ymax = terra::ymax(r_tmpl)
      )
      can_ext <- as.numeric(canonical_mask$extent)
      
      if (length(r_ext) != 4L ||
          length(can_ext) != 4L ||
          any(!is.finite(r_ext)) ||
          any(!is.finite(can_ext)) ||
          any(abs(r_ext - can_ext) >
              1e-8 * pmax(1, abs(can_ext)))) {
        stop(
          "SPEI grid extent does not match canonical Script-4 grid: ",
          paste(signif(r_ext, 10), collapse = ", "),
          " vs ",
          paste(signif(can_ext, 10), collapse = ", ")
        )
      }
    }
    
    if (!is.null(canonical_mask$crs)) {
      r_crs   <- terra::crs(r_tmpl)
      can_crs <- as.character(canonical_mask$crs)
      
      if (!identical(
        trimws(r_crs),
        trimws(can_crs)
      )) {
        stop(
          "SPEI grid CRS does not match canonical Script-4 CRS.\n",
          "Current:   ", r_crs, "\n",
          "Canonical: ", can_crs
        )
      }
    }
    
    cat("  ✓ Canonical grid geometry validated from metadata\n")
    
    valid_pix <- valid_pix & canonical_mask_vec
    
    ts_mat    <- ts_mat[valid_pix, , drop = FALSE]
    lon_vec   <- lon_vec[valid_pix]
    lat_vec   <- lat_vec[valid_pix]
    n_pix     <- sum(valid_pix)
    cat(sprintf("  Valid pixels after canonical mask: %d  |  %.1f%% non-NA values\n",
                n_pix, 100 * sum(!is.na(ts_mat)) / length(ts_mat)))
    
    # Cell area weights — fractional coverage for boundary pixels
    # Each pixel's effective weight = full cell area × fraction of that cell
    # that actually lies inside the basin polygon (0–1, via cover = TRUE).
    # Interior pixels get cover = 1.0 (unchanged from the old approach).
    # Edge pixels that only partially overlap the basin are down-weighted
    # proportionally, rather than being counted at their full cell area.
    cat("  → Computing fractional-coverage area weights...\n")
    basin_v        <- terra::project(load_basin_vect(BASIN_SHP), terra::crs(r_tmpl))
    cell_area_rast <- terra::cellSize(r_tmpl, unit = "m")
    cover_rast     <- terra::rasterize(basin_v, r_tmpl, cover = TRUE)
    area_all       <- as.numeric(terra::values(cell_area_rast, na.rm = FALSE))
    cover_all      <- as.numeric(terra::values(cover_rast,     na.rm = FALSE))
    cover_all[is.na(cover_all)] <- 0   # pixels fully outside basin → weight 0
    eff_area_all   <- area_all * cover_all
    area_valid     <- eff_area_all[valid_pix]
    bad_w <- is.na(area_valid) | !is.finite(area_valid) | area_valid <= 0
    if (any(bad_w)) {
      med_w <- median(area_valid[!bad_w], na.rm = TRUE)
      area_valid[bad_w] <- if (is.finite(med_w) && med_w > 0) med_w else 1
      cat(sprintf("  ⚠ %d pixel(s) with zero/invalid effective area replaced with median (%.0f m²)\n",
                  sum(bad_w), med_w))
    }
    n_partial <- sum(cover_all[valid_pix] > 0 & cover_all[valid_pix] < 1,
                     na.rm = TRUE)
    cat(sprintf("  ✓ Fractional area weights: range %.0f – %.0f m² (median %.0f m²)\n",
                min(area_valid), max(area_valid), median(area_valid)))
    cat(sprintf("  ✓ Boundary pixels with partial coverage: %d\n", n_partial))
    
    # ── ONE-TIME RECONCILIATION CHECK  [NEW] ──────────────────────────────────
    # Script 4 (4pr_pet_trends_v2.r) clips ERA5-Land Pr/PET/Tair rasters to the
    # basin polygon and reports 1804 bbox → 856 basin pixels. Every SPI/SPEI/
    # SPEI_ERA5 block here instead reports "n_pix" valid (non-NA) pixels out of
    # a totally different grid (e.g. 2021 cells → 865 valid here). These are
    # two different underlying rasters (raw ERA5-Land vs. the upstream SPI/SPEI
    # NetCDFs), so it is not guaranteed they represent the same basin mask.
    # This runs ONCE (on the first index × scale processed) and compares:
    #   (a) how many of THIS grid's "valid" (non-NA) pixels actually fall
    #       inside the basin polygon used for area-weighting above,
    #   (b) how many in-basin cells of THIS grid are missing (NA) data,
    #   (c) the total basin area implied by this pipeline (sum of effective
    #       cell area) vs. the area implied by Script 4's clipped_template.rds
    #       mask (856 ERA5-Land pixels), if that file is available.
    # A CAUTION is logged (not a hard stop) if the two masks disagree
    # meaningfully, since results from the two pipelines are compared side by
    # side downstream.
    if (!exists(".basin_reconciliation_done", envir = .GlobalEnv)) {
      assign(".basin_reconciliation_done", TRUE, envir = .GlobalEnv)
      cat("\n  ── Basin pixel-mask reconciliation check (runs once) ──\n")
      tryCatch({
        n_valid_here      <- length(valid_pix[valid_pix])
        n_valid_in_basin  <- sum(cover_all[valid_pix] > 0, na.rm = TRUE)
        n_valid_outside   <- n_valid_here - n_valid_in_basin
        n_inbasin_missing <- sum(cover_all > 0 & !valid_pix, na.rm = TRUE)
        area_here_km2     <- sum(area_valid, na.rm = TRUE) / 1e6
        
        cat(sprintf("    This grid (%s): %d valid pixels | %d fall inside basin polygon | %d fall OUTSIDE basin polygon\n",
                    index_type, n_valid_here, n_valid_in_basin, n_valid_outside))
        cat(sprintf("    In-basin cells on this grid with NO data (NA): %d\n", n_inbasin_missing))
        cat(sprintf("    Basin area implied by this grid (cover-weighted): %.1f km²\n", area_here_km2))
        
        if (n_valid_outside > 0) {
          cat(sprintf("    ⚠ CAUTION: %d 'valid' pixel(s) on this grid lie outside the basin polygon used for area weighting.\n",
                      n_valid_outside))
        }
        
        # Try to cross-check against Script 4's ERA5-Land basin mask, if present.
        #
        # [BUGFIX] Originally this read clipped_template.rds and counted non-NA
        # cells. That file is NOT a saved mask -- Script 4 deliberately blanks
        # it to all-NA (`values(clipped_template) <- NA_real_`) as an empty
        # geometry shell for other scripts to fill in later. Reading it here
        # always returned 0 pixels / 0 km^2, producing a meaningless "+Inf%"
        # caution regardless of the true agreement between the two masks.
        # Fix: read the actual pixel count and cell resolution from
        # analysis_metadata.rds instead, which genuinely stores them
        # (`basin_pixels` and `clipped_extent$res`).
        metadata_path <- file.path(WD_PATH, "trend_analysis_pr_pet", "analysis_metadata.rds")
        if (file.exists(metadata_path)) {
          era5_meta <- tryCatch(readRDS(metadata_path), error = function(e) NULL)
          if (!is.null(era5_meta) && !is.null(era5_meta$basin_pixels) &&
              !is.null(era5_meta$clipped_extent$res)) {
            n_era5_pixels  <- era5_meta$basin_pixels
            era5_cell_res  <- era5_meta$clipped_extent$res
            era5_cell_area <- if (length(era5_cell_res) >= 2) {
              era5_cell_res[1] * era5_cell_res[2]
            } else {
              era5_cell_res[1]^2
            }
            era5_area_km2 <- n_era5_pixels * era5_cell_area / 1e6
            pct_diff <- 100 * (area_here_km2 - era5_area_km2) / era5_area_km2
            cat(sprintf("    Script 4 ERA5-Land basin mask: %d pixels | area %.1f km²\n",
                        n_era5_pixels, era5_area_km2))
            cat(sprintf("    Area difference (this grid vs. ERA5-Land mask): %+.1f%%\n", pct_diff))
            if (abs(pct_diff) > 10) {
              cat(sprintf(paste0("    ⚠ CAUTION: basin area implied by the %s grid differs from Script 4's ",
                                 "ERA5-Land mask by more than 10%%. Confirm the two pipelines are using ",
                                 "the SAME basin boundary/resampling before comparing raw-climate-trend ",
                                 "results (Script 4/8) against drought-index-trend results (Script 6/7) side by side.\n"),
                          index_type))
            } else {
              cat("    ✓ Basin areas agree within 10% — the two pixel sets plausibly represent the same basin.\n")
            }
          } else {
            cat("    ⚠ analysis_metadata.rds found but missing basin_pixels/clipped_extent$res — skipping cross-check.\n")
          }
        } else {
          cat(sprintf("    ℹ analysis_metadata.rds not found at %s — cannot cross-check against Script 4's mask.\n",
                      metadata_path))
          cat("      (Run Script 4 first, or update the path above, to enable this check.)\n")
        }
      }, error = function(e) {
        cat(sprintf("    ⚠ Basin reconciliation check failed: %s\n", conditionMessage(e)))
      })
      cat("  ─────────────────────────────────────────────────────\n\n")
    }
    
    # Area-weighted basin average (computed once, reused below)
    cat("  → Computing area-weighted basin average time series...\n")
    basin_avg_vec <- vapply(seq_len(ncol(ts_mat)), function(t) {
      v  <- ts_mat[, t]
      ok <- !is.na(v) & is.finite(v)
      if (!any(ok)) return(NA_real_)
      sum(v[ok] * area_valid[ok]) / sum(area_valid[ok])
    }, numeric(1L))
    cat(sprintf("  ✓ Basin average computed: %d time steps, range [%.3f, %.3f]\n",
                length(basin_avg_vec),
                min(basin_avg_vec, na.rm = TRUE),
                max(basin_avg_vec, na.rm = TRUE)))
    
    yr_seq <- rep(years, each = 12L) + (rep(0:11, times = n_years) / 12.0)
    
    # [4/5] Compute all statistics
    cat("[4/5] Computing statistics...\n")
    results <- data.table::data.table(lon = lon_vec, lat = lat_vec)
    
    cat("  → Mann-Kendall (variance-corrected)...    ")
    results <- cbind(results, vectorized_mann_kendall(ts_mat))
    results[, `:=`(tau_vc_raw = tau_vc,
                   p_value_vc_raw = p_value_vc,
                   sl_vc_raw = sl_vc)]
    cat("✓\n")
    
    cat("  → Mann-Kendall (TFPW)...    ")
    results <- cbind(results, vectorized_mann_kendall_tfpw(ts_mat))
    results[, `:=`(tau_tfpw_raw = tau_tfpw,
                   p_value_tfpw_raw = p_value_tfpw,
                   sl_tfpw_raw = sl_tfpw)]
    cat("✓\n")
    
    # Method 1: S&W single threshold
    cat("  → Event statistics [M1: S&W, intensity ref = onset_thr]...    ")
    m1_stats <- vectorized_event_stats(ts_mat,
                                       onset_thr    = DROUGHT_ONSET,
                                       end_thr      = DROUGHT_ONSET,
                                       min_duration = 1L)
    results <- cbind(results, m1_stats)
    cat("✓\n")
    
    # Method 2: Hysteresis
    cat("  → Event statistics [M2: Hysteresis, min_dur=", min_dur,
        "m, intensity ref = end_thr]...    ")
    m2_stats <- vectorized_event_stats(ts_mat,
                                       onset_thr    = DROUGHT_ONSET,
                                       end_thr      = DROUGHT_END,
                                       min_duration = min_dur)
    data.table::setnames(m2_stats, paste0(names(m2_stats), "_hyst"))
    results <- cbind(results, m2_stats)
    cat("✓\n")
    
    cat("  ── Method 1 (S&W) per-pixel summaries:\n")
    cat(sprintf("     D3-6 events  — mean: %.2f\n", mean(results$n_events_D36,  na.rm = TRUE)))
    cat(sprintf("     D7-12 events — mean: %.2f\n", mean(results$n_events_D712, na.rm = TRUE)))
    cat(sprintf("     D13+ events  — mean: %.2f  (max: %d)\n",
                mean(results$n_events_D13p, na.rm = TRUE),
                max(results$n_events_D13p,  na.rm = TRUE)))
    cat(sprintf("     mean_I: %.4f  |  mean_S: %.4f\n",
                mean(results$mean_I, na.rm = TRUE), mean(results$mean_S, na.rm = TRUE)))
    
    cat("  ── Method 2 (Hysteresis) per-pixel summaries:\n")
    cat(sprintf("     D3-6 events  — mean: %.2f\n",
                mean(results$n_events_D36_hyst,  na.rm = TRUE)))
    cat(sprintf("     D7-12 events — mean: %.2f\n",
                mean(results$n_events_D712_hyst, na.rm = TRUE)))
    cat(sprintf("     D13+ events  — mean: %.2f  (max: %d)\n",
                mean(results$n_events_D13p_hyst, na.rm = TRUE),
                max(results$n_events_D13p_hyst,  na.rm = TRUE)))
    cat(sprintf("     mean_I_hyst: %.4f  |  mean_S_hyst: %.4f\n",
                mean(results$mean_I_hyst, na.rm = TRUE),
                mean(results$mean_S_hyst, na.rm = TRUE)))
    
    # PELT regime-shift detection (was defined but never called before)
    cat("  → PELT regime-shift detection...    ")
    results <- cbind(results, vectorized_regime_shift_pelt(ts_mat, years))
    cat("✓\n")
    
    # Wald-Wolfowitz runs test on hysteresis-defined drought events
    cat(sprintf("  → Runs test [Hyst drought clusters, min_dur=%d m]...    ", min_dur))
    results <- cbind(results, vectorized_runs_test(ts_mat,
                                                   onset_thr    = DROUGHT_ONSET,
                                                   end_thr      = DROUGHT_END,
                                                   min_duration = min_dur))
    cat("✓\n")
    
    # Spectral peak detection (drought-band periods >= 24 months)
    # Uses local MC red-noise override (n_sim=500) to match Script 4 null hypothesis.
    cat("  → Spectral peak detection (periods >= 24 mo, MC  AR(1) red-noise test, n_sim=500)...    ")
    results <- cbind(results, vectorized_spectral_peaks(ts_mat))
    n_pix_with_peaks <- sum(results$n_spectral_peaks > 0L, na.rm = TRUE)
    cat(sprintf("✓  (%d/%d pixels have >= 1 significant peak)\n",
                n_pix_with_peaks, n_pix))
    
    # [4b/5] Area-weighted basin-average CSV
    # ── FIX: verbose error reporting + CSV-from-scratch fallback ──────────────
    # The original tryCatch swallowed all errors silently, leaving no
    # basin_averaged_by_month.csv for script 7 to read.  We now:
    #   (a) print the full error message if save_basin_avg_from_pixels() fails,
    #   (b) fall back to generate_basin_avg_from_csvs() which rebuilds the CSV
    #       directly from the per-calendar-month CSVs already in seas_dir,
    #       bypassing the in-memory ts_mat entirely.
    cat("  → Saving area-weighted basin-average CSV...\n")
    seas_dir_map <- list(spi      = SPI_SEAS_DIR,
                         spei     = SPEI_SEAS_DIR,
                         spei_era5 = SPEI_ERA5_SEAS_DIR)
    basin_csv_dir <- seas_dir_map[[tolower(index_type)]]
    if (is.null(basin_csv_dir) || !nzchar(basin_csv_dir))
      basin_csv_dir <- file.path(WD_PATH, paste0(index_type, "_results_seasonal"))
    dir.create(basin_csv_dir, showWarnings = FALSE, recursive = TRUE)
    
    basin_avg_saved <- tryCatch({
      save_basin_avg_from_pixels(ts_mat, area_valid, dates_ts,
                                 index_type, sc, basin_csv_dir)
      TRUE
    }, error = function(e) {
      cat(sprintf("  ⚠ save_basin_avg_from_pixels failed: %s\n", e$message))
      cat("    → Attempting CSV-from-scratch fallback...\n")
      FALSE
    })
    
    if (!basin_avg_saved) {
      # Fallback: rebuild basin-average CSV from the 12 per-month CSVs
      # that were already written by 1SPI_ERALand_v2.R / 3SPEI_ERALand_v2.R.
      # Uses cos(lat) area-weighting for geographic coords, plain mean otherwise.
      tryCatch({
        generate_basin_avg_from_csvs(index_type, sc, basin_csv_dir)
      }, error = function(e2)
        cat(sprintf("  ✗ Fallback also failed: %s\n", e2$message)))
    }
    
    # Basin extent + class MK
    cat("  → Basin extent time series & class-level MK trends...    ")
    index_label <- sprintf("%s-%02d", toupper(index_type), sc)
    extent_out  <- file.path(TREND_DIR, "basin_extent")
    tryCatch(
      compute_basin_extent_and_class_trends(
        ts_mat, n_years, index_label, extent_out,
        start_year   = min(years),
        area_weights = area_valid,
        basin_avg    = basin_avg_vec
      ),
      error = function(e) {
        cat(sprintf("\n  ⚠ Basin extent analysis failed: %s\n", e$message))
        NULL
      }
    )
    cat("✓\n")
    
    # [4c/5] BH False Discovery Rate correction
    cat("  → BH False Discovery Rate correction (across basin pixels)...\n")
    apply_fdr <- function(p_raw) {
      ok    <- !is.na(p_raw)
      p_fdr <- rep(NA_real_, length(p_raw))
      if (sum(ok) > 1L) p_fdr[ok] <- p.adjust(p_raw[ok], method = "BH")
      p_fdr
    }
    results$p_fdr_vc_raw   <- apply_fdr(results$p_value_vc_raw)
    results$p_fdr_tfpw_raw <- apply_fdr(results$p_value_tfpw_raw)
    results$p_fdr_vc       <- apply_fdr(results$p_value_vc)
    results$p_fdr_tfpw     <- apply_fdr(results$p_value_tfpw)
    
    report_fdr <- function(label, p_raw, p_fdr) {
      n_test <- sum(!is.na(p_raw))
      n_nom  <- sum(p_raw < 0.05,  na.rm = TRUE)
      n_fdr  <- sum(p_fdr < 0.05,  na.rm = TRUE)
      cat(sprintf("  ✓ FDR (%s): %d nominal → %d FDR-significant  (%+d pixels)\n",
                  label, n_nom, n_fdr, n_fdr - n_nom))
      list(n_test = n_test, n_nom = n_nom, n_fdr = n_fdr)
    }
    fdr_vc   <- report_fdr("VC-HR98", results$p_value_vc,   results$p_fdr_vc)
    fdr_tfpw <- report_fdr("TFPW   ", results$p_value_tfpw, results$p_fdr_tfpw)
    
    # [5/5] Save main results CSV
    out_file <- file.path(TREND_DIR, sprintf("%s_%02d_results.csv", index_type, sc))
    data.table::fwrite(results, out_file, showProgress = FALSE)
    cat(sprintf("  ✅ %s  (%.2f MB)\n",
                basename(out_file), file.info(out_file)$size / 1024^2))
    cat(sprintf("  %-35s  nominal: %d / %d   FDR (BH): %d / %d\n",
                "Significant VC trends  (p < 0.05):",
                fdr_vc$n_nom,  fdr_vc$n_test, fdr_vc$n_fdr,  fdr_vc$n_test))
    cat(sprintf("  %-35s  nominal: %d / %d   FDR (BH): %d / %d\n",
                "Significant TFPW trends (p < 0.05):",
                fdr_tfpw$n_nom, fdr_tfpw$n_test, fdr_tfpw$n_fdr, fdr_tfpw$n_test))
    
    # Duration-class map CSVs
    dur_cols_sw <- c("lon", "lat",
                     "n_events_D36", "n_events_D712", "n_events_D13p", "n_D4p",
                     "mean_I", "mean_S")
    data.table::fwrite(
      results[, ..dur_cols_sw],
      file.path(TREND_DIR, sprintf("%s_%02d_duration_class_map_SW.csv", index_type, sc)),
      showProgress = FALSE
    )
    cat(sprintf("  ✅ Duration-class map (S&W): %s_%02d_duration_class_map_SW.csv\n",
                index_type, sc))
    
    dur_cols_hyst <- c("lon", "lat",
                       "n_events_D36_hyst", "n_events_D712_hyst", "n_events_D13p_hyst",
                       "n_D4p_hyst", "mean_I_hyst", "mean_S_hyst")
    data.table::fwrite(
      results[, ..dur_cols_hyst],
      file.path(TREND_DIR, sprintf("%s_%02d_duration_class_map_Hyst.csv", index_type, sc)),
      showProgress = FALSE
    )
    cat(sprintf("  ✅ Duration-class map (Hyst): %s_%02d_duration_class_map_Hyst.csv\n",
                index_type, sc))
    
    # D13+ frequency rasters (both methods)
    if (exists("r_tmpl") && !is.null(r_tmpl)) {
      write_d13p_nc <- function(values_vec, suffix, label_suffix) {
        r_out     <- r_tmpl
        terra::values(r_out) <- NA_real_
        all_vals  <- rep(NA_real_, terra::ncell(r_tmpl))
        all_vals[valid_pix]  <- values_vec
        terra::values(r_out) <- all_vals
        names(r_out) <- sprintf("%s_%02d_n_D13p_%s", index_type, sc, label_suffix)
        nc_path <- file.path(TREND_DIR,
                             sprintf("%s_%02d_D13p_pixel_frequency_%s.nc",
                                     index_type, sc, suffix))
        terra::writeCDF(r_out, nc_path, overwrite = TRUE)
        cat(sprintf("  ✅ D13praster (%s): %s\n", suffix, basename(nc_path)))
      }
      write_d13p_nc(results$n_events_D13p,      "SW",   "SW_longterm_droughts")
      write_d13p_nc(results$n_events_D13p_hyst, "Hyst", "Hyst_longterm_droughts")
    }
    
    rm(ts_mat, stacks)
    invisible(gc())
  }
}

# ================================================================================
#   RUN ALL INDICES
# ================================================================================
cat("\n╔════════════════════════════════════════════════╗\n")
cat("║  DROUGHT TREND TEST         ║\n")
cat("║  Dual methods: S&W single-threshold + Hyst     ║\n")
cat("╚════════════════════════════════════════════════╝\n\n")
cat(sprintf("  Indices  : SPI, SPEI_PM, SPEI_ERA5\n"))
cat(sprintf("  Scales   : %s\n",
            paste(SPI_SCALES, collapse = ", ")))
cat(sprintf("  Method 1 (S&W)  : onset = %.1f, exit = %.1f, no min duration\n",
            DROUGHT_ONSET, DROUGHT_ONSET))
cat(sprintf("  Method 2 (Hyst) : onset = %.1f, exit = %.1f, scale-specific min duration\n",
            DROUGHT_ONSET, DROUGHT_END))
cat("    Min durations: SPI/SPEI-1 >= 2m | -3 >= 3m | -6 >= 4m | -12 >= 6m\n\n")

total_start <- Sys.time()

run_index_safe <- function(index_type, scales, find_fn, seas_dir) {
  tryCatch(
    process_index(index_type, scales, find_fn, seas_dir),
    error = function(e) {
      cat(sprintf("\n  ❌ %s failed and was skipped: %s\n",
                  toupper(index_type), conditionMessage(e)))
    }
  )
}

cat("\n── SPI ──\n")
run_index_safe("spi", SPI_SCALES, find_seasonal_nc_files, SPI_SEAS_DIR)

cat("\n── SPEI (PM) ──\n")
run_index_safe("spei", SPEI_SCALES, find_seasonal_nc_files, SPEI_SEAS_DIR)

# SPEI_ERA5 — ERA5-Land PEV PET variant written by 3SPEI_ERALand_v2.R
# into spei_results_seasonal_era5land/.  Files may be named with index prefix
# "spei_era5" or, if 3SPEI_ERALand_v2.R wrote them with the plain "spei" prefix
# inside the era5-specific subdirectory, the wrapper tries both automatically
# so that no manual filename wrangling is required.
find_spei_era5_nc_safe <- function(data_dir, index_type, scale) {
  # First attempt: files named spei_era5_NN_*.nc (preferred convention)
  files <- tryCatch(
    find_seasonal_nc_files(data_dir, "spei_era5", scale),
    error = function(e) character(0)
  )
  # Fallback: files named spei_NN_*.nc stored in the era5 subdirectory
  if (!length(files)) {
    files <- tryCatch(
      find_seasonal_nc_files(data_dir, "spei", scale),
      error = function(e) {
        cat(sprintf("  ⚠ SPEI_ERA5 file finder error (both patterns failed): %s\n",
                    conditionMessage(e)))
        character(0)
      }
    )
    if (length(files))
      cat(sprintf("  ℹ SPEI_ERA5 scale-%02d: using 'spei' filename prefix in %s\n",
                  scale, basename(data_dir)))
  }
  files
}

if (dir.exists(SPEI_ERA5_SEAS_DIR)) {
  cat("\n── SPEI (ERA5-Land PEV) ──\n")
  run_index_safe("spei_era5", SPEI_ERA5_SCALES, find_spei_era5_nc_safe, SPEI_ERA5_SEAS_DIR)
} else {
  cat("\n── SPEI (ERA5-Land PEV): directory not found, skipped ──\n")
  cat(sprintf("   Expected: %s\n", SPEI_ERA5_SEAS_DIR))
  cat("   Run 3SPEI_ERALand_v2.R first.\n")
}

# ────────────────────────────────────────────────────────────────────────────────
#   FINAL SUMMARY
# ────────────────────────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════\n")
cat("DUAL DROUGHT EVENT METHOD SUMMARY  [FIXED VERSION]\n")
cat("══════════════════════════════════════════════════════\n")
cat("Method 1 (S&W)   — Sheffield & Wood (2008):\n")
cat("  Onset:        index < -0.5\n")
cat("  Termination:  index >= -0.5  (same threshold, no hysteresis)\n")
cat("  Min duration: none (all events retained)\n")
cat("  Intensity:    mean(onset_thr - x)  [deficit below -0.5]\n")
cat("  CSV columns: n_events, mean_duration, max_intensity,\n")
cat("               n_events_D36, n_events_D712, n_events_D13p,\n")
cat("               n_D4p, mean_I, mean_S\n")
cat("  Duration map: *_duration_class_map_SW.csv\n")
cat("  D13+ raster:  *_D13p_pixel_frequency_SW.nc\n\n")
cat("Method 2 (Hyst)  — Hysteresis with scale-specific min duration:\n")
cat("  Onset:        index < -0.5\n")
cat("  Termination:  index >= -0.5  (drought persists through brief recoveries)\n")
cat("  Intensity:    mean(end_thr - x) = mean(-x)  [deficit below 0.0, always >= 0]\n")
cat("  Min duration per scale:\n")
cat("    SPI/SPEI-1  : >= 2 months\n")
cat("    SPI/SPEI-3  : >= 3 months\n")
cat("    SPI/SPEI-6  : >= 4 months\n")
cat("    SPI/SPEI-12 : >= 6 months\n")
cat("    SPEI_ERA5    : same scale-specific rules as SPEI_PM\n")
cat("  CSV columns: n_events_hyst, mean_duration_hyst, max_intensity_hyst,\n")
cat("               n_events_D36_hyst, ..., mean_I_hyst, mean_S_hyst\n")
cat("  Duration map: *_duration_class_map_Hyst.csv\n")
cat("  D13+ raster:  *_D13p_pixel_frequency_Hyst.nc\n\n")
cat("Additional per-pixel columns (all indices):\n")
cat("  regime_shift_year, regime_shift_detected, n_changepoints,\n")
cat("  mean_before_shift, mean_after_shift, magnitude_shift   [PELT]\n")
cat("  p_value_runs, clustering, filtered_runs                [Runs test]\n")
cat("  n_spectral_peaks                                       [MC red-noise spectral, n_sim=500]\n")
cat("══════════════════════════════════════════════════════\n")