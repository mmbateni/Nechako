# ==============================================================================
#   w1_trend_test.R  ·  SPATIAL TREND ANALYSIS FOR ALL INDICES
# [MODIFIED VERSION — DUAL DROUGHT EVENT DETECTION METHODS]
# ==============================================================================
#   Computes per-pixel statistics for SPI, SPEI (multiple scales) and SWEI:
#   • Mann-Kendall — Hamed & Rao (1998) variance-corrected (VC-HR98) + Sen's slope
#     Uses modifiedmk::mmkh() which adjusts Var(S) for lag-1…lag-(n-1) autocorrelation.
#     Tau is the standard concordance statistic; p-value is the H-R98 corrected one.
#     Sen's slope is extracted directly from mmkh() (avoids the O(n²) matrix method).
#   • Mann-Kendall (TFPW – Trend-Free Pre-Whitening, Yue et al. 2002)
#   • Drought event count / mean duration / max intensity — TWO METHODS:
#       Method 1 (S&W): single threshold DROUGHT_ONSET, no hysteresis, no min duration
#       Method 2 (Hyst): onset < DROUGHT_ONSET, termination >= DROUGHT_END,
#                        scale-specific minimum duration
#   • PELT  regime-shift year
#   • Wald-Wolfowitz runs test (temporal clustering)
#   • n_spectral_peaks  (placeholder = 0L; reserved for future spectral analysis)
#
# OUTPUT: {TREND_DIR}/{index}_{scale:02d}_results.csv  (one per index × scale)
# Columns include tau_vc, p_value_vc (H-R98 corrected), p_fdr_vc (BH-FDR corrected),
# tau_tfpw, p_value_tfpw, p_fdr_tfpw (BH-FDR corrected),
# n_events, mean_duration, max_intensity            [Method 1: S&W]
# n_events_hyst, mean_duration_hyst, max_intensity_hyst  [Method 2: Hysteresis]
# n_events_D46, n_events_D712, n_events_D12p, n_D4p, mean_I, mean_S  [Method 1]
# n_events_D46_hyst, …_hyst                                           [Method 2]
# regime_shift_year, p_value_runs, clustering, n_spectral_peaks.
#
# p_fdr_* are computed by applying p.adjust(method="BH") across all valid basin
# pixels within each index × scale run. Use p_fdr_vc for significance mapping.
#
# METHOD COMPARISON NOTES:
#   Method 1 (Sheffield & Wood 2008): single threshold -1.0, no hysteresis,
#     no minimum duration filter. Every month below -1.0 belongs to an event.
#     Use for direct replication of S&W (2008) Table 4 statistics.
#   Method 2 (Hysteresis): onset <-1.0, termination >=0.0.
#     Once a drought begins (drops below -1.0), it persists until index rises
#     to 0.0 or above, preventing spurious terminations during brief recoveries.
#     Scale-specific minimum duration filters out noise at short scales.
#     Columns suffixed _hyst. Use for climatologically meaningful event counts.
#
# NOTE: Raw pixel time-series are NOT stored in the CSV (they live in the source
#       NetCDF files). This keeps CSVs lightweight and avoids ~6 M redundant values.
# Run BEFORE w4_trends_visualization.R
# ================================================================================
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "data.table", "Kendall", "trend", "modifiedmk",
                      "parallel", "lubridate","changepoint"))
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
dir.create(TREND_DIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------------
#   Utility operators
# --------------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

####################################################################################
# PARALLEL BACK-END (WINDOWS-COMPATIBLE)
####################################################################################
is_windows <- .Platform$OS.type == "windows"
N_CORES    <- max(1L, parallel::detectCores() - 1L)
cl         <- NULL
if (N_CORES > 1) {
  if (is_windows) {
    cl <- parallel::makeCluster(N_CORES)
    parallel::clusterEvalQ(cl, {
      library(Kendall)
      library(trend)
      library(modifiedmk)
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
#   Default (SWEI, etc.): 3 months
#
# Arguments:
#   scale               — integer accumulation scale (1, 3, 6, 12, …)
#   min_duration_map    — named list mapping scale-as-character to integer months
#   default_min_duration — used when scale not found in the map
get_min_duration <- function(scale,
                             min_duration_map     = list("1"  = 2L,
                                                         "3"  = 3L,
                                                         "6"  = 4L,
                                                         "12" = 6L),
                             default_min_duration = 3L) {
  key <- as.character(as.integer(scale))
  val <- min_duration_map[[key]]
  if (!is.null(val)) as.integer(val) else as.integer(default_min_duration)
}

# ================================================================================
# VECTORISED STATISTICS FUNCTIONS
# ================================================================================
# --------------------------------------------------------------------------------
# 1. Hamed-Rao (1998) Variance-Corrected Mann-Kendall + Sen's slope
# --------------------------------------------------------------------------------
# WHAT THE CORRECTION DOES
#   Standard MK assumes i.i.d. observations.  Monthly drought indices have
#   positive serial correlation, which inflates Var(S) and makes the standard
#   p-value anti-conservative.  Hamed & Rao (1998) derive a correction factor
#   n_s / n that rescales Var(S) using the lag-1 … lag-(n-1) autocorrelations
#   of the RANKED series.  modifiedmk::mmkh() implements this exactly.
#
# OUTPUTS
#   tau_vc     — Kendall's tau (concordance statistic; same as standard MK,
#                because the H-R98 correction only adjusts Var(S), not S itself)
#   p_value_vc — two-tailed p-value using the H-R98 corrected variance
#   sl_vc      — Sen's slope (units: index-value / month), computed by mmkh()
#                via the O(n·log n) Theil-Sen estimator — NOT the O(n²) matrix
#   filtered_vc— TRUE if the pixel was skipped (too few valid values or ~zero
#                variance); all output columns are NA for filtered pixels.
#
# MINIMUM OBSERVATIONS: 12 (same as TFPW).
# --------------------------------------------------------------------------------
vectorized_mann_kendall <- function(ts_matrix) {
  n_pix <- nrow(ts_matrix)
  
  process_pixel <- function(i) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    
    # ── Guard: need ≥ 12 valid obs and non-degenerate variance ───────────────
    if (length(x) < 12L || var(x, na.rm = TRUE) < 1e-6)
      return(list(tau = NA_real_, pval = NA_real_, slope = NA_real_,
                  filtered = TRUE))
    
    tryCatch({
      # ── tau: from standard Kendall S (H-R98 does not change tau itself) ───
      tau_val <- as.numeric(Kendall::MannKendall(x)$tau)
      
      # ── H-R98 corrected p-value + Sen's slope from mmkh() ─────────────────
      hr  <- modifiedmk::mmkh(x)
      
      list(
        tau      = tau_val,
        pval     = as.numeric(hr["new P-value"]),
        slope    = as.numeric(hr["Sen's slope"]),
        filtered = FALSE
      )
    }, error = function(e)
      list(tau = NA_real_, pval = NA_real_, slope = NA_real_, filtered = TRUE))
  }
  
  do_par <- N_CORES > 1L && n_pix > 100L
  res <- if (do_par && is_windows) {
    parallel::parLapply(cl, seq_len(n_pix), process_pixel)
  } else if (do_par) {
    parallel::mclapply(seq_len(n_pix), process_pixel, mc.cores = N_CORES)
  } else {
    lapply(seq_len(n_pix), process_pixel)
  }
  
  data.table::data.table(
    tau_vc      = sapply(res, `[[`, "tau"),
    p_value_vc  = sapply(res, `[[`, "pval"),
    sl_vc       = sapply(res, `[[`, "slope"),
    filtered_vc = sapply(res, `[[`, "filtered")
  )
}

# --------------------------------------------------------------------------------
#   2. TFPW Mann-Kendall (Yue et al. 2002)
# --------------------------------------------------------------------------------
# Removes trend + pre-whitens lag-1 autocorrelation before MK test.
# Returns tau_tfpw, p_value_tfpw, sl_tfpw, filtered_tfpw.
vectorized_mann_kendall_tfpw <- function(ts_matrix) {
  n_pix <- nrow(ts_matrix)
  process_pixel_tfpw <- function(i) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 12 || var(x, na.rm = TRUE) < 1e-6)
      return(list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
    tryCatch({
      n  <- length(x)
      
      # Sen's slope via trend::sens.slope
      sens_r  <- trend::sens.slope(x)
      beta    <- as.numeric(sens_r$estimates)
      
      if (!is.finite(beta)) return(list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
      
      # Detrend
      t_seq  <- seq_len(n)
      x_det  <- x - beta * t_seq
      
      # Lag-1 autocorrelation of detrended series
      r1  <- cor(x_det[-n], x_det[-1], use = "complete.obs")
      if (is.na(r1)) r1 <- 0
      r1  <- max(-0.99, min(0.99, r1))
      
      # Pre-whiten
      x_pw  <- x_det[-1] - r1 * x_det[-n]
      
      # Blend trend back in (Yue et al. 2002, step 4)
      x_final  <- x_pw + beta * t_seq[-1]
      
      if (length(x_final) < 10) return(list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
      
      mk  <- Kendall::MannKendall(x_final)
      # MODIFICATION: Calculate final slope on the blended pre-whitened series
      slope_final <- as.numeric(trend::sens.slope(x_final)$estimates)
      
      list(tau = as.numeric(mk$tau), pval = as.numeric(mk$sl), slope = slope_final, filtered = FALSE)
    }, error = function(e) list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
  }
  do_par  <- N_CORES  > 1  && n_pix  > 100
  res  <- if (do_par  && is_windows) {
    parallel::parLapply(cl, seq_len(n_pix), process_pixel_tfpw)
  } else if (do_par) {
    parallel::mclapply(seq_len(n_pix), process_pixel_tfpw, mc.cores = N_CORES)
  } else {
    lapply(seq_len(n_pix), process_pixel_tfpw)
  }
  data.table::data.table(
    tau_tfpw      = sapply(res, `[[`, "tau"),
    p_value_tfpw  = sapply(res, `[[`, "pval"),
    sl_tfpw       = sapply(res, `[[`, "slope"), # MODIFICATION: Extract new slope column
    filtered_tfpw = sapply(res, `[[`, "filtered")
  )
}
# --------------------------------------------------------------------------------
#   3. Drought event characteristics — unified for both methods
# --------------------------------------------------------------------------------
# Returns n_events, mean_duration, max_intensity per pixel.
#
# METHOD SELECTION via onset_thr / end_thr / min_duration:
#   Method 1 (Sheffield & Wood 2008): call with end_thr = onset_thr, min_duration = 1L
#     Entry : x < onset_thr  (single threshold)
#     Exit  : x >= onset_thr (same threshold — no hysteresis)
#     Filter: none (min_duration = 1 keeps all events)
#
#   Method 2 (Hysteresis): call with end_thr > onset_thr, min_duration from
#     get_min_duration(scale)
#     Entry : x < onset_thr  (e.g. -1.0)
#     Exit  : x >= end_thr   (e.g. 0.0) — drought persists through brief recovery
#     Filter: events shorter than min_duration are discarded
#
# IMPLEMENTATION NOTE:
#   A state-tracking loop is used (rather than vectorised diff) because hysteresis
#   requires knowledge of the current drought state: after onset, months where
#   onset_thr <= x < end_thr keep the drought active even though x > onset_thr.
#   The diff-based approach cannot handle this state dependency.
#
# Arguments:
#   ts_matrix    — pixels × months matrix (rows = pixels)
#   onset_thr    — value below which a drought begins  (default: DROUGHT_ONSET)
#   end_thr      — value at/above which a drought ends (default: onset_thr → S&W)
#   min_duration — minimum event length in months; shorter events are discarded
#                  (default: 1L — no filter, preserving S&W baseline behaviour)
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
vectorized_event_stats <- function(ts_matrix,
                                   onset_thr    = DROUGHT_ONSET,
                                   end_thr      = onset_thr,
                                   min_duration = 1L) {
  # Validate thresholds
  if (end_thr < onset_thr)
    warning(sprintf(
      "vectorized_event_stats: end_thr (%.2f) < onset_thr (%.2f). ",
      "For hysteresis, end_thr should be >= onset_thr. Proceeding anyway.",
      end_thr, onset_thr))
  
  n_pix <- nrow(ts_matrix)
  # Pre‑allocate output vectors
  n_events      <- integer(n_pix)
  mean_duration <- numeric(n_pix)
  max_intensity <- numeric(n_pix)
  n_events_D46  <- integer(n_pix)
  n_events_D712 <- integer(n_pix)
  n_events_D12p <- integer(n_pix)
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
    ints  <- numeric(0)      # store min value per event (for max_intensity)
    I_evs <- numeric(0)      # store mean I per event (for mean_I)
    S_evs <- numeric(0)      # store severity per event (for mean_S)
    
    for (j in seq_along(x)) {
      v <- x[j]
      
      if (!in_d && v < onset_thr) {
        # Drought onset
        in_d  <- TRUE
        s_idx <- j
        
      } else if (in_d && v >= end_thr) {
        # Drought termination
        dur <- j - s_idx                     # number of months in drought
        if (dur >= min_duration) {
          durs <- c(durs, dur)
          # min value during event (most negative) -> max_intensity
          ints <- c(ints, min(x[s_idx:(j - 1L)], na.rm = TRUE))
          # mean deficit relative to onset_thr (intensity I)
          I_ev <- mean(onset_thr - x[s_idx:(j - 1L)], na.rm = TRUE)
          I_evs <- c(I_evs, I_ev)
          S_evs <- c(S_evs, I_ev * dur)
        }
        in_d  <- FALSE
        s_idx <- NA_integer_
      }
      # else: inside drought but not yet at end_thr – stay in drought
    }
    
    # Close any event still open at end of record
    if (in_d && !is.na(s_idx)) {
      dur <- length(x) - s_idx + 1L
      if (dur >= min_duration) {
        durs <- c(durs, dur)
        ints <- c(ints, min(x[s_idx:length(x)], na.rm = TRUE))
        I_ev <- mean(onset_thr - x[s_idx:length(x)], na.rm = TRUE)
        I_evs <- c(I_evs, I_ev)
        S_evs <- c(S_evs, I_ev * dur)
      }
    }
    
    # Populate output vectors
    if (length(durs) == 0L) {
      n_events[i]      <- 0L
      mean_duration[i] <- 0
      max_intensity[i] <- 0
      n_events_D46[i]  <- 0L
      n_events_D712[i] <- 0L
      n_events_D12p[i] <- 0L
      n_D4p[i]         <- 0L
      mean_I[i]        <- NA_real_
      mean_S[i]        <- NA_real_
    } else {
      n_events[i]      <- length(durs)
      mean_duration[i] <- mean(durs)
      max_intensity[i] <- min(ints)               # most negative value
      # Duration class counts
      n_events_D46[i]  <- sum(durs >= 4L  & durs <= 6L)
      n_events_D712[i] <- sum(durs >= 7L  & durs <= 12L)
      n_events_D12p[i] <- sum(durs >= 13L)
      n_D4p[i]         <- sum(durs >= 4L)
      mean_I[i]        <- mean(I_evs, na.rm = TRUE)
      mean_S[i]        <- mean(S_evs, na.rm = TRUE)
    }
  }
  
  data.table::data.table(
    n_events      = n_events,
    mean_duration = mean_duration,
    max_intensity = max_intensity,
    n_events_D46  = n_events_D46,
    n_events_D712 = n_events_D712,
    n_events_D12p = n_events_D12p,
    n_D4p         = n_D4p,
    mean_I        = mean_I,
    mean_S        = mean_S
  )
}

# --------------------------------------------------------------------------------
#   4. PELT  regime-shift year
# --------------------------------------------------------------------------------
# NOTE: regime_shift_year is kept as the primary column name so that
#       w4_trends_visualization.R (Fig 3, Panel B) works without modification.
vectorized_regime_shift_pelt <- function(ts_matrix, years, min_obs = 20L) {
  # years: integer vector of length n_years (one entry per year in the record).
  # ts_matrix: pixels × (n_years * 12) in Jan-first interleaved order.
  n_pix  <- nrow(ts_matrix)
  n_yrs  <- length(years)
  
  # Pre-allocate output columns
  out <- data.table::data.table(
    regime_shift_year     = rep(NA_real_,    n_pix),
    regime_shift_detected = rep(FALSE,       n_pix),
    n_changepoints        = rep(0L,          n_pix),
    mean_before_shift     = rep(NA_real_,    n_pix),
    mean_after_shift      = rep(NA_real_,    n_pix),
    magnitude_shift       = rep(NA_real_,    n_pix)
  )
  
  for (i in seq_len(n_pix)) {
    row <- ts_matrix[i, ]
    
    # ── Aggregate monthly columns to annual means ─────────────────────────────
    # ts_matrix columns are in interleaved order: Jan_yr1, Feb_yr1, …, Dec_yr1,
    # Jan_yr2, …  So year k occupies columns ((k-1)*12+1) : (k*12).
    ann <- vapply(seq_len(n_yrs), function(k) {
      cols <- ((k - 1L) * 12L + 1L):(k * 12L)
      vals <- row[cols[cols <= length(row)]]
      if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
    }, numeric(1L))
    
    ok <- !is.na(ann)
    if (sum(ok) < min_obs) next
    ann_clean  <- ann[ok]
    yrs_clean  <- years[ok]
    n          <- length(ann_clean)
    
    tryCatch({
      cpt_obj    <- changepoint::cpt.mean(ann_clean,
                                          method  = "PELT",
                                          penalty = "BIC")
      cpts       <- changepoint::cpts(cpt_obj)
      
      if (length(cpts) > 0L) {
        first_cpt  <- cpts[1L]
        mb         <- mean(ann_clean[seq_len(first_cpt)],         na.rm = TRUE)
        ma         <- mean(ann_clean[(first_cpt + 1L):n],         na.rm = TRUE)
        
        data.table::set(out, i, "regime_shift_year",     as.numeric(yrs_clean[first_cpt]))
        data.table::set(out, i, "regime_shift_detected", TRUE)
        data.table::set(out, i, "n_changepoints",        length(cpts))
        data.table::set(out, i, "mean_before_shift",     mb)
        data.table::set(out, i, "mean_after_shift",      ma)
        data.table::set(out, i, "magnitude_shift",       ma - mb)
      }
      # If cpts is empty: detected = FALSE, all NAs — already initialised above
    }, error = function(e) NULL)   # leave row as NA on any PELT error
  }
  out
}

# --------------------------------------------------------------------------------
#   5. Wald-Wolfowitz runs test (temporal clustering)
# --------------------------------------------------------------------------------
vectorized_runs_test <- function(ts_matrix, threshold = 0) {
  n_pix <- nrow(ts_matrix)
  res   <- data.table::data.table(
    p_value_runs  = rep(NA_real_,      n_pix),
    clustering    = rep(NA_character_, n_pix),
    filtered_runs = rep(TRUE,          n_pix)
  )
  for (i in seq_len(n_pix)) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 10) next
    b    <- as.integer(x  >= threshold)
    rn   <- rle(b)
    nr   <- length(rn$lengths)
    n    <- length(b)
    n1   <- sum(b)
    n0   <- n - n1
    
    if (n1 == 0 || n0 == 0) next
    
    er   <- 2 * n0 * n1 / n + 1
    vr   <- 2 * n0 * n1 * (2 * n0 * n1 - n) / (n^2 * (n - 1))
    
    if (vr  <= 0) next
    
    z      <- (nr - er) / sqrt(vr)
    pv     <- 2 * pnorm(-abs(z))
    clus   <- if (nr  < er)  "clustered" else if (nr  > er)  "dispersed" else  "random"
    
    res[i, `:=`(p_value_runs  = pv,
                clustering    = clus,
                filtered_runs = FALSE)]
  }
  res
}

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
  ## ── area_weights validation ───────────────────────────────────────────────
  use_weights <- !is.null(area_weights) &&
    length(area_weights) == nrow(ts_matrix) &&
    all(is.finite(area_weights)) &&
    all(area_weights > 0)
  if (!is.null(area_weights) && !use_weights)
    warning("compute_basin_extent_and_class_trends: area_weights ignored (length mismatch or invalid values); falling back to equal weighting")
  
  utils_load_packages(c("ggplot2", "Kendall", "dplyr", "patchwork"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  n_months   <- n_years * 12L
  month_seq  <- seq_len(n_months)
  yr_seq     <- start_year + floor((month_seq - 1L) / 12L)
  mo_seq     <- ((month_seq - 1L) %% 12L) + 1L
  date_seq   <- as.Date(sprintf("%04d-%02d-01", yr_seq, mo_seq))
  
  #------------------------------------------------------------------------------
  # 7a. Monthly % area in drought (uses onset_thr for area fraction — S&W consistent)
  #------------------------------------------------------------------------------
  n_pix    <- nrow(ts_matrix)
  n_col    <- min(n_months, ncol(ts_matrix))
  pct_area <- numeric(n_col)
  for (t in seq_len(n_col)) {
    col_vals <- ts_matrix[, t]
    ok       <- !is.na(col_vals)
    if (!any(ok)) { pct_area[t] <- NA; next }
    if (use_weights) {
      w_ok         <- area_weights[ok]
      drought_ok   <- ok & col_vals < onset_thr
      pct_area[t]  <- 100.0 * sum(area_weights[drought_ok], na.rm = TRUE) /
        sum(w_ok, na.rm = TRUE)
    } else {
      pct_area[t]  <- 100.0 * sum(ok & col_vals < onset_thr) / sum(ok)
    }
  }
  extent_ts <- data.frame(
    date     = date_seq[seq_len(n_col)],
    year     = yr_seq[seq_len(n_col)],
    month    = mo_seq[seq_len(n_col)],
    pct_area = pct_area,
    row.names = NULL
  )
  write.csv(extent_ts,
            file.path(out_dir, sprintf("%s_basin_extent_monthly.csv", index_label)),
            row.names = FALSE)
  
  #------------------------------------------------------------------------------
  # 7b. Mann-Kendall on annual mean extent
  #------------------------------------------------------------------------------
  extent_annual  <- extent_ts %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(mean_pct = mean(pct_area, na.rm = TRUE), .groups = "drop")
  x_clean  <- extent_annual$mean_pct[!is.na(extent_annual$mean_pct)]
  mk_ext   <- list(tau = NA_real_, p_value = NA_real_, sens_slope = NA_real_)
  if (length(x_clean) >= 10L) {
    tau_val  <- as.numeric(Kendall::MannKendall(x_clean)$tau)
    hr_r     <- tryCatch(modifiedmk::mmkh(x_clean), error = function(e) NULL)
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
  mk_ext_df  <- data.frame(
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
  
  #------------------------------------------------------------------------------
  # 7c. Basin-average time series → event catalogue (uses hysteresis thresholds)
  #------------------------------------------------------------------------------
  if (is.null(basin_avg)) {   # <-- compute only if not provided
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
  basin_ts <- data.frame(date = date_seq[seq_len(n_col)], value = basin_avg)
  classify_dur <- function(d)
    dplyr::case_when(d >= 4L  & d <= 6L  ~ "D4-6 (Short-term)",
                     d >= 7L  & d <= 12L ~ "D7-12 (Medium-term)",
                     d >= 13L            ~ "D12+ (Long-term)",
                     TRUE                ~ "D1-3 (Sub-threshold)")
  
  vals       <- basin_avg
  in_d       <- FALSE
  s_idx      <- NA_integer_
  event_list <- list()
  for (j in seq_along(vals)) {
    v <- vals[j]
    if (!in_d && !is.na(v) && v < onset_thr) {
      in_d  <- TRUE; s_idx <- j
    } else if (in_d && !is.na(v) && v >= end_thr) {
      dur <- j - s_idx
      event_list[[length(event_list) + 1L]] <- data.frame(
        start_date      = date_seq[s_idx],
        end_date        = date_seq[j - 1L],
        start_year      = yr_seq[s_idx],
        duration_months = dur,
        duration_class  = classify_dur(dur),
        mean_severity   = mean(abs(pmin(vals[s_idx:(j-1)], 0)), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      in_d  <- FALSE; s_idx <- NA_integer_
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
  
  #------------------------------------------------------------------------------
  # 7e. Figure: basin extent time series with MK annotation
  #------------------------------------------------------------------------------
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
  
  x_min <- min(extent_annual$year, na.rm = TRUE)
  x_max <- max(extent_annual$year, na.rm = TRUE)
  first_break <- ceiling(x_min / 10) * 10
  last_break  <- floor(x_max / 10) * 10
  x_breaks    <- seq(first_break, last_break, by = 10)
  x_breaks    <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  if (length(x_breaks) == 0)  x_breaks <- c(x_min, x_max)
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
  #------------------------------------------------------------------------------
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
    
    # ── Scale-specific minimum duration for Method 2 (Hysteresis) ─────────────
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
    
    n_lyrs_per_file <- terra::nlyr(stacks[[1]])
    dates_m1        <- extract_dates_from_nc(files[1], n_lyrs_per_file)
    years           <- as.integer(format(dates_m1, "%Y"))
    n_years         <- length(years)
    
    r_tmpl  <- stacks[[1]][[1]]
    xy      <- terra::xyFromCell(r_tmpl, seq_len(terra::ncell(r_tmpl)))
    lon_vec <- xy[, 1]
    lat_vec <- xy[, 2]
    n_pix   <- length(lon_vec)
    cat(sprintf("  Grid: %d pixels × %d years × 12 months\n", n_pix, n_years))
    
    # [3/5] Build time-series matrix (pixels × months)
    cat("[3/5] Building interleaved time-series matrix...\n")
    ts_mat <- matrix(NA_real_, nrow = n_pix, ncol = n_years * 12)
    for (m in seq_along(stacks)) {
      n_yrs_m <- min(terra::nlyr(stacks[[m]]), n_years)
      for (y in seq_len(n_yrs_m)) {
        col_idx           <- (y - 1) * 12 + m
        ts_mat[, col_idx] <- as.numeric(terra::values(stacks[[m]][[y]]))
      }
    }
    
    valid_pix <- rowSums(!is.na(ts_mat)) > 0
    ts_mat    <- ts_mat[valid_pix, , drop = FALSE]
    lon_vec   <- lon_vec[valid_pix]
    lat_vec   <- lat_vec[valid_pix]
    n_pix     <- sum(valid_pix)
    cat(sprintf("  Valid pixels: %d  |  %.1f%% non-NA values\n",
                n_pix, 100 * sum(!is.na(ts_mat)) / length(ts_mat)))
    
    # ── Cell area weights ─────────────────────────────────────────────────────
    cat("  → Computing cell area weights (BC Albers)...\n")
    cell_area_rast <- terra::cellSize(r_tmpl, unit = "m")
    area_all       <- as.numeric(terra::values(cell_area_rast, na.rm = FALSE))
    area_valid     <- area_all[valid_pix]
    bad_w <- is.na(area_valid) | !is.finite(area_valid) | area_valid <= 0
    if (any(bad_w)) {
      med_w <- median(area_valid[!bad_w], na.rm = TRUE)
      area_valid[bad_w] <- if (is.finite(med_w) && med_w > 0) med_w else 1
    }
    cat(sprintf("  ✓ Area weights: range %.0f – %.0f m² (median %.0f m²)\n",
                min(area_valid), max(area_valid), median(area_valid)))
    
    # ── Compute area-weighted basin average ONCE ──────────────────────────────
    # This single vector is reused by both save_basin_avg_from_pixels() (via the
    # CSV writer) and compute_basin_extent_and_class_trends() (for event
    # detection).  Computing it here avoids the redundant vapply that would
    # otherwise run a second time inside compute_basin_extent_and_class_trends()
    # when basin_avg = NULL.
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
    
    dates_ts <- as.Date(sprintf("%04d-%02d-01",
                                rep(years, each = 12L),
                                rep(1:12, times = n_years)))
    yr_seq   <- rep(years, each = 12) + (rep(0:11, times = n_years) / 12)
    
    # [4/5] Compute all statistics
    cat("[4/5] Computing statistics...\n")
    results <- data.table::data.table(lon = lon_vec, lat = lat_vec)
    
    cat("  → Mann-Kendall (variance-corrected)...    ")
    results <- cbind(results, vectorized_mann_kendall(ts_mat))
    cat("✓\n")
    
    cat("  → Mann-Kendall (TFPW)...    ")
    results <- cbind(results, vectorized_mann_kendall_tfpw(ts_mat))
    cat("✓\n")
    
    # ── Drought event detection: Method 1 — S&W single threshold ─────────────
    # Entry and exit at the same threshold (DROUGHT_ONSET); no min duration.
    # Replicates Sheffield & Wood (2008) methodology directly.
    # ── Method 1: S&W single threshold (no hysteresis, min_duration = 1) ─────
    cat("  → Event statistics [M1: S&W]...    ")
    m1_stats <- vectorized_event_stats(ts_mat,
                                       onset_thr    = DROUGHT_ONSET,
                                       end_thr      = DROUGHT_ONSET,
                                       min_duration = 1L)
    results <- cbind(results, m1_stats)
    cat("✓\n")
    
    # ── Method 2: Hysteresis with scale‑specific min duration ────────────────
    cat("  → Event statistics [M2: Hysteresis, min_dur=", min_dur, "m]...    ")
    m2_stats <- vectorized_event_stats(ts_mat,
                                       onset_thr    = DROUGHT_ONSET,
                                       end_thr      = DROUGHT_END,
                                       min_duration = min_dur)
    # Rename all columns by adding "_hyst" suffix
    data.table::setnames(m2_stats, paste0(names(m2_stats), "_hyst"))
    results <- cbind(results, m2_stats)
    cat("✓\n")
    
    # ── Print per-method summary statistics ───────────────────────────────────
    cat("  ── Method 1 (S&W) per-pixel summaries:\n")
    cat(sprintf("     D4-6 events  — mean: %.2f\n", mean(results$n_events_D46,  na.rm=TRUE)))
    cat(sprintf("     D7-12 events — mean: %.2f\n", mean(results$n_events_D712, na.rm=TRUE)))
    cat(sprintf("     D12+ events  — mean: %.2f  (max: %d)\n",
                mean(results$n_events_D12p, na.rm=TRUE),
                max(results$n_events_D12p,  na.rm=TRUE)))
    cat(sprintf("     mean_I: %.4f  |  mean_S: %.4f\n",
                mean(results$mean_I, na.rm=TRUE), mean(results$mean_S, na.rm=TRUE)))
    
    cat("  ── Method 2 (Hysteresis) per-pixel summaries:\n")
    cat(sprintf("     D4-6 events  — mean: %.2f\n",
                mean(results$n_events_D46_hyst,  na.rm=TRUE)))
    cat(sprintf("     D7-12 events — mean: %.2f\n",
                mean(results$n_events_D712_hyst, na.rm=TRUE)))
    cat(sprintf("     D12+ events  — mean: %.2f  (max: %d)\n",
                mean(results$n_events_D12p_hyst, na.rm=TRUE),
                max(results$n_events_D12p_hyst,  na.rm=TRUE)))
    cat(sprintf("     mean_I_hyst: %.4f  |  mean_S_hyst: %.4f\n",
                mean(results$mean_I_hyst, na.rm=TRUE),
                mean(results$mean_S_hyst, na.rm=TRUE)))
    
    # ── [4b/5] Area-weighted basin-average CSV ────────────────────────────────
    cat("  → Saving area-weighted basin-average CSV...\n")
    seas_dir_map <- list(spi = SPI_SEAS_DIR, spei = SPEI_SEAS_DIR, swei = SWEI_SEAS_DIR)
    basin_csv_dir <- seas_dir_map[[tolower(index_type)]]
    if (is.null(basin_csv_dir))
      basin_csv_dir <- file.path(WD_PATH, paste0(index_type, "_results_seasonal"))
    tryCatch(
      save_basin_avg_from_pixels(ts_mat, area_valid, dates_ts,
                                 index_type, sc, basin_csv_dir),
      error = function(e)
        cat(sprintf("  ⚠ Basin-avg CSV save failed: %s\n", e$message))
    )
    
    # ── Basin extent + class MK (uses hysteresis thresholds for event catalog) ─
    cat("  → Basin extent time series & class-level MK trends...    ")
    index_label <- sprintf("%s-%02d", toupper(index_type), sc)
    extent_out  <- file.path(TREND_DIR, "basin_extent")
    basin_extent_results <- tryCatch(
      compute_basin_extent_and_class_trends(
        ts_mat, n_years, index_label, extent_out,
        start_year   = min(years),
        area_weights = area_valid,
        basin_avg    = basin_avg_vec   # ← precomputed above; no redundant vapply
      ),
      error = function(e) {
        cat(sprintf("\n  ⚠ Basin extent analysis failed: %s\n", e$message))
        NULL
      }
    )
    cat("✓\n")
    
    # ── [4c/5] BH False Discovery Rate correction ─────────────────────────────
    cat("  → BH False Discovery Rate correction (across basin pixels)...\n")
    
    apply_fdr <- function(p_raw) {
      ok   <- !is.na(p_raw)
      p_fdr <- rep(NA_real_, length(p_raw))
      if (sum(ok) > 1L) p_fdr[ok] <- p.adjust(p_raw[ok], method = "BH")
      p_fdr
    }
    
    results$p_fdr_vc    <- apply_fdr(results$p_value_vc)
    results$p_fdr_tfpw  <- apply_fdr(results$p_value_tfpw)
    
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
                fdr_vc$n_nom,   fdr_vc$n_test,
                fdr_vc$n_fdr,   fdr_vc$n_test))
    cat(sprintf("  %-35s  nominal: %d / %d   FDR (BH): %d / %d\n",
                "Significant TFPW trends (p < 0.05):",
                fdr_tfpw$n_nom, fdr_tfpw$n_test,
                fdr_tfpw$n_fdr, fdr_tfpw$n_test))
    
    # ── Duration-class map CSVs (Method 1 and Method 2) ──────────────────────
    # Method 1 (S&W)
    dur_cols_sw <- c("lon", "lat",
                     "n_events_D46", "n_events_D712", "n_events_D12p", "n_D4p",
                     "mean_I", "mean_S")
    data.table::fwrite(
      results[, ..dur_cols_sw],
      file.path(TREND_DIR, sprintf("%s_%02d_duration_class_map_SW.csv", index_type, sc)),
      showProgress = FALSE
    )
    cat(sprintf("  ✅ Duration-class map (S&W): %s_%02d_duration_class_map_SW.csv\n",
                index_type, sc))
    
    # Method 2 (Hysteresis)
    dur_cols_hyst <- c("lon", "lat",
                       "n_events_D46_hyst", "n_events_D712_hyst", "n_events_D12p_hyst",
                       "n_D4p_hyst", "mean_I_hyst", "mean_S_hyst")
    data.table::fwrite(
      results[, ..dur_cols_hyst],
      file.path(TREND_DIR, sprintf("%s_%02d_duration_class_map_Hyst.csv", index_type, sc)),
      showProgress = FALSE
    )
    cat(sprintf("  ✅ Duration-class map (Hyst): %s_%02d_duration_class_map_Hyst.csv\n",
                index_type, sc))
    
    # ── D12+ frequency rasters (both methods) ────────────────────────────────
    if (exists("r_tmpl") && !is.null(r_tmpl)) {
      
      write_d12p_nc <- function(values_vec, suffix, label_suffix) {
        r_out <- r_tmpl
        terra::values(r_out) <- NA_real_
        all_vals <- rep(NA_real_, terra::ncell(r_tmpl))
        all_vals[valid_pix] <- values_vec
        terra::values(r_out) <- all_vals
        names(r_out) <- sprintf("%s_%02d_n_D12p_%s", index_type, sc, label_suffix)
        nc_path <- file.path(TREND_DIR,
                             sprintf("%s_%02d_D12p_pixel_frequency_%s.nc",
                                     index_type, sc, suffix))
        terra::writeCDF(r_out, nc_path, overwrite = TRUE)
        cat(sprintf("  ✅ D12+ raster (%s): %s\n", suffix, basename(nc_path)))
      }
      
      write_d12p_nc(results$n_events_D12p,      "SW",   "SW_longterm_droughts")
      write_d12p_nc(results$n_events_D12p_hyst,  "Hyst", "Hyst_longterm_droughts")
    }
    
    rm(ts_mat, stacks)
    invisible(gc())
  }
}

# ================================================================================
#   RUN ALL INDICES
# ================================================================================
cat("\n╔════════════════════════════════════════════════╗\n")
cat("║  DROUGHT TREND TEST  (w1)                      ║\n")
cat("║  Dual methods: S&W single-threshold + Hyst     ║\n")
cat("╚════════════════════════════════════════════════╝\n\n")
cat(sprintf("  Method 1 (S&W)  : onset = %.1f, exit = %.1f, no min duration\n",
            DROUGHT_ONSET, DROUGHT_ONSET))
cat(sprintf("  Method 2 (Hyst) : onset = %.1f, exit = %.1f, scale-specific min duration\n",
            DROUGHT_ONSET, DROUGHT_END))
cat("    Min durations: SPI/SPEI-1 >= 2m | -3 >= 3m | -6 >= 4m | -12 >= 6m\n\n")

total_start <- Sys.time()

# --------------------------------------------------------------------------------
#   Safe wrapper: run process_index and continue on any error
# --------------------------------------------------------------------------------
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

cat("\n── SPEI ──\n")
run_index_safe("spei", SPEI_SCALES, find_seasonal_nc_files, SPEI_SEAS_DIR)

# SWEI: find_swei_seasonal_files takes (data_dir, scale) only — wrap safely
find_swei_nc_safe <- function(data_dir, index_type, scale) {
  tryCatch(
    find_swei_seasonal_files(data_dir, scale),
    error = function(e) {
      cat(sprintf("  ⚠ SWEI file finder error: %s\n", conditionMessage(e)))
      character(0)
    }
  )
}

cat("\n── SWEI ──\n")
run_index_safe("swei", SWEI_SCALE, find_swei_nc_safe, SWEI_SEAS_DIR)

# MSPI / MSPEI (optional)
find_mspi_mspei_files <- function(data_dir, index_type, scale) {
  pattern <- sprintf("%s_monthly_.*\\.nc$", tolower(index_type))
  files   <- list.files(data_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0)
    cat(sprintf("  ⚠ No %s NC files in %s — skipping\n",
                toupper(index_type), data_dir))
  files
}

MSPI_MSPEI_SCALE <- 1L
mspi_dir  <- file.path(WD_PATH, "mspi_results")
mspei_dir <- file.path(WD_PATH, "mspei_results")

if (dir.exists(mspi_dir)) {
  cat("\n── MSPI ──\n")
  run_index_safe("mspi", MSPI_MSPEI_SCALE, find_mspi_mspei_files, mspi_dir)
} else {
  cat("\n── MSPI: directory not found, skipped ──\n")
}

if (dir.exists(mspei_dir)) {
  cat("\n── MSPEI ──\n")
  run_index_safe("mspei", MSPI_MSPEI_SCALE, find_mspi_mspei_files, mspei_dir)
} else {
  cat("\n── MSPEI: directory not found, skipped ──\n")
}

# --------------------------------------------------------------------------------
#   Shut down parallel cluster
# --------------------------------------------------------------------------------
if (!is.null(cl)) {
  parallel::stopCluster(cl)
  cat("\n✓ Cluster stopped\n")
}

elapsed <- Sys.time() - total_start
cat(sprintf("\n✅ Done.  Time: %.1f min\n  Next: run w4_trends_visualization.R\n\n",
            as.numeric(elapsed, units = "mins")))

# --------------------------------------------------------------------------------
#   FINAL SUMMARY
# --------------------------------------------------------------------------------
cat("══════════════════════════════════════════════════════\n")
cat("DUAL DROUGHT EVENT METHOD SUMMARY\n")
cat("══════════════════════════════════════════════════════\n")
cat("Method 1 (S&W)   — Sheffield & Wood (2008):\n")
cat("  Onset:       index < -1.0\n")
cat("  Termination: index >= -1.0  (same threshold, no hysteresis)\n")
cat("  Min duration: none (all events retained)\n")
cat("  CSV columns: n_events, mean_duration, max_intensity,\n")
cat("               n_events_D46, n_events_D712, n_events_D12p,\n")
cat("               n_D4p, mean_I, mean_S\n")
cat("  Duration map: *_duration_class_map_SW.csv\n")
cat("  D12+ raster:  *_D12p_pixel_frequency_SW.nc\n\n")
cat("Method 2 (Hyst)  — Hysteresis with scale-specific min duration:\n")
cat("  Onset:       index < -1.0\n")
cat("  Termination: index >= 0.0  (hysteresis — drought persists\n")
cat("               through brief recoveries above onset threshold)\n")
cat("  Min duration per scale:\n")
cat("    SPI/SPEI-1  : >= 2 months\n")
cat("    SPI/SPEI-3  : >= 3 months\n")
cat("    SPI/SPEI-6  : >= 4 months\n")
cat("    SPI/SPEI-12 : >= 6 months\n")
cat("    SWEI / other: >= 3 months (default)\n")
cat("  CSV columns: n_events_hyst, mean_duration_hyst, max_intensity_hyst,\n")
cat("               n_events_D46_hyst, n_events_D712_hyst, n_events_D12p_hyst,\n")
cat("               n_D4p_hyst, mean_I_hyst, mean_S_hyst\n")
cat("  Duration map: *_duration_class_map_Hyst.csv\n")
cat("  D12+ raster:  *_D12p_pixel_frequency_Hyst.nc\n")
cat("══════════════════════════════════════════════════════\n")