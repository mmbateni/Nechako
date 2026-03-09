# ================================================================================
#   w3_trend_test.R  ·  SPATIAL TREND ANALYSIS FOR ALL INDICES
# [MODIFIED VERSION WITH ALTERNATIVE DROUGHT DEFINITION]
# ================================================================================
#   Computes per-pixel statistics for SPI, SPEI (multiple scales) and SWEI:
#   • Mann-Kendall (variance-corrected, VC)  +  Sen's slope
# • Mann-Kendall (TFPW – Trend-Free Pre-Whitening, Yue et al. 2002)
# • Drought event count / mean duration / max intensity
# • CUSUM regime-shift year
# • Wald-Wolfowitz runs test (temporal clustering)
# • n_spectral_peaks  (placeholder = 0L; reserved for future spectral analysis)
# • [NEW] Alternative drought events (Onset -1.0, Term -0.5, scale-specific min duration)
# OUTPUT: {TREND_DIR}/{index}_{scale:02d}_results.csv  (one per index × scale)
# Columns include tau_vc, p_value_vc, tau_tfpw, p_value_tfpw,
# n_events, mean_duration, max_intensity, regime_shift_year,
# p_value_runs, clustering, n_spectral_peaks, and time_001…time_NNN.
# Run BEFORE w4_trends_visualization.R
# ================================================================================
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "data.table", "Kendall", "trend", "parallel", "lubridate"))
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
dir.create(TREND_DIR, showWarnings = FALSE, recursive = TRUE)
# --------------------------------------------------------------------------------
#   Utility operators
# --------------------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b


####################################################################################
# ALTERNATIVE DROUGHT DEFINITION CONFIGURATION (USER REQUESTED)
####################################################################################
# Store these in memory for later reuse (ALT_EVENTS_MEMORY)
# Onset: -1.0, Termination: -0.5 (reduces toggling near zero)
# Minimum durations vary by scale to filter blips
ALT_DROUGHT_DEF <- list(
  onset_threshold = -1.0,
  termination_threshold = -0.5,
  min_duration_map = list(
    "1"  = 2,   # SPI/SPEI-1: >= 2 months
    "3"  = 3,   # SPI/SPEI-3: >= 3 months
    "6"  = 4,   # SPI/SPEI-6: >= 4 months
    "12" = 6    # SPI/SPEI-12: >= 6 months
  ),
  default_min_duration = 3 # For SWEI or others not specified
)

# Global memory storage for alternative events (accessible after script runs)
ALT_EVENTS_MEMORY <- list()

####################################################################################
# PARALLEL BACK-END (WINDOWS-COMPATIBLE)
####################################################################################
is_windows <- .Platform$OS.type == "windows"
N_CORES    <- max(1L, parallel::detectCores() - 1L)
cl         <- NULL
if (N_CORES > 1) {
  if (is_windows) {
    cl <- parallel::makeCluster(N_CORES)
    parallel::clusterEvalQ(cl, { library(Kendall); library(trend) })
    cat(sprintf("✓ Windows cluster: %d cores\n", N_CORES))
  } else {
    cat(sprintf("✓ Unix fork: %d cores\n", N_CORES))
  }
} else {
  cat("ℹ Single-core mode\n")
}

####################################################################################
# UNIFIED 3-ARG WRAPPER FOR THE SWEI FILE FINDER
####################################################################################
# process_index() always calls find_fn(dir, index_type, scale).
# find_swei_seasonal_files() only takes 2 args, so we wrap it here.
find_swei_nc_wrapper <- function(data_dir, index_type, scale) {
  find_swei_seasonal_files(data_dir, scale)
}

# ================================================================================
# VECTORISED STATISTICS FUNCTIONS
# ================================================================================
# --------------------------------------------------------------------------------
# 1. Variance-corrected Mann-Kendall + Sen's slope
# --------------------------------------------------------------------------------
# Returns tau_vc, p_value_vc, sl_vc (Sen's slope), filtered_vc.
vectorized_mann_kendall <- function(ts_matrix) {
  n_pix  <- nrow(ts_matrix)
  process_pixel <- function(i) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 10 || var(x, na.rm = TRUE) < 1e-6)
      return(list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
    tryCatch({
      mk   <- Kendall::MannKendall(x)
      n    <- length(x)
      idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)   # all i < j pairs
      slopes <- (x[idx[, 2]] - x[idx[, 1]]) / (idx[, 2] - idx[, 1])
      
      list(tau      = as.numeric(mk$tau),
           pval     = as.numeric(mk$sl),
           slope    = median(slopes, na.rm = TRUE),
           filtered = FALSE)
    }, error = function(e) list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
  }
  do_par  <- N_CORES  > 1  && n_pix  > 100
  res  <- if (do_par  && is_windows) {
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
# Returns tau_tfpw, p_value_tfpw, filtered_tfpw.
vectorized_mann_kendall_tfpw <- function(ts_matrix) {
  n_pix <- nrow(ts_matrix)
  process_pixel_tfpw <- function(i) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 12 || var(x, na.rm = TRUE) < 1e-6)
      return(list(tau = NA, pval = NA, filtered = TRUE))
    tryCatch({
      n  <- length(x)
      
      # Sen's slope via trend::sens.slope — avoids the O(n²) matrix approach
      # and correctly uses all n*(n-1)/2 pairwise slopes.
      sens_r  <- trend::sens.slope(x)
      beta    <- as.numeric(sens_r$estimates)
      
      if (!is.finite(beta)) return(list(tau = NA, pval = NA, filtered = TRUE))
      
      # Detrend
      t_seq  <- seq_len(n)
      x_det  <- x - beta * t_seq
      
      # Lag-1 autocorrelation of detrended series
      r1  <- cor(x_det[-n], x_det[-1], use = "complete.obs")
      if (is.na(r1)) r1 <- 0
      r1  <- max(-0.99, min(0.99, r1))   # clip to avoid instability
      
      # Pre-whiten: remove lag-1 autocorrelation (result has length n-1)
      x_pw  <- x_det[-1] - r1 * x_det[-n]
      
      # Blend trend back in (Yue et al. 2002, step 4)
      x_final  <- x_pw + beta * t_seq[-1]
      
      if (length(x_final) < 10) return(list(tau = NA, pval = NA, filtered = TRUE))
      
      mk  <- Kendall::MannKendall(x_final)
      list(tau = as.numeric(mk$tau), pval = as.numeric(mk$sl), filtered = FALSE)
    }, error = function(e) list(tau = NA, pval = NA, filtered = TRUE))
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
    filtered_tfpw = sapply(res, `[[`, "filtered")
  )
}

# --------------------------------------------------------------------------------
#   3. Drought event characteristics
# --------------------------------------------------------------------------------
# Returns n_events, mean_duration, max_intensity per pixel.
vectorized_event_detection <- function(ts_matrix, threshold = DROUGHT_ONSET) {
  n_pix <- nrow(ts_matrix)
  res   <- data.table::data.table(
    n_events      = rep(NA_real_, n_pix),
    mean_duration = rep(NA_real_, n_pix),
    max_intensity = rep(NA_real_, n_pix)
  )
  for (i in seq_len(n_pix)) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (!length(x)) next
    in_ev <- x < threshold
    
    if (!any(in_ev)) {
      res[i, c("n_events", "mean_duration", "max_intensity") := list(0, 0, 0)]
      next
    }
    
    starts <- which(diff(c(FALSE, in_ev)) == 1)
    ends   <- which(diff(c(in_ev, FALSE)) == -1)
    durs   <- ends - starts + 1
    ints   <- mapply(function(s, e) min(x[s:e]), starts, ends)
    
    res[i, `:=`(n_events      = length(starts),
                mean_duration = mean(durs),
                max_intensity = min(ints))]
  }
  res
}

# --------------------------------------------------------------------------------
#   4. CUSUM regime-shift year
# --------------------------------------------------------------------------------
# Identifies the year of maximum CUSUM deviation (most likely shift point).
vectorized_regime_shift <- function(ts_matrix, year_vec) {
  n_pix <- nrow(ts_matrix)
  years <- rep(NA_real_, n_pix)
  for (i in seq_len(n_pix)) {
    x  <- ts_matrix[i, ]
    ok <- !is.na(x)
    if (sum(ok) < 20) next
    xc <- x[ok]
    yc <- year_vec[ok]
    cs <- cumsum(xc - mean(xc))
    years[i] <- yc[which.max(abs(cs))]
  }
  years
}

# --------------------------------------------------------------------------------
#   5. Wald-Wolfowitz runs test (temporal clustering)
# --------------------------------------------------------------------------------
# Returns p_value_runs, clustering ("clustered"/"dispersed"/"random"),
# filtered_runs (TRUE if series too short or degenerate).
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
#   NEW FUNCTION 6 — Per-pixel event counts by duration class + S&W intensity/severity
# ================================================================================
#   # For every pixel, count events by Sheffield & Wood (2008) duration class and
#   # compute the mean S&W intensity and severity across ALL events at that pixel.
#   # Output columns:
#   # n_events_D46   – short-term events (4–6 months)
#   # n_events_D712  – medium-term events (7–12 months)
#   # n_events_D12p  – long-term events (>12 months)
#   # n_D4p          – all events ≥4 months (reference total)
#   # mean_I         – S&W intensity: mean deficit below threshold, averaged over all
#   # events (= mean over events of mean(threshold − x[event_months]))
#   # mean_S         – S&W severity: mean of (I × duration) averaged over all events
vectorized_duration_class_counts <- function(ts_matrix,
                                             threshold = DROUGHT_ONSET) {
  # Sheffield & Wood (2008) use a SINGLE quantile threshold (q0 = 10th percentile,
  # mapped to SPI/SPEI ≈ -1.28; here DROUGHT_ONSET is used as the fixed threshold).
  # A drought event begins when x < threshold AND ends when x >= threshold.
  # The same value governs both entry and exit — no hysteresis, no second threshold.
  # This guarantees that every month inside an event has x < threshold,
  # so the deficit (threshold − x) is strictly positive and mean_I is always >= 0.
  n_pix <- nrow(ts_matrix)
  res   <- data.table::data.table(
    n_events_D46  = rep(NA_integer_, n_pix),
    n_events_D712 = rep(NA_integer_, n_pix),
    n_events_D12p = rep(NA_integer_, n_pix),
    n_D4p         = rep(NA_integer_, n_pix),
    mean_I        = rep(NA_real_,    n_pix),
    mean_S        = rep(NA_real_,    n_pix)
  )
  for (i in seq_len(n_pix)) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x)  < 4L) {
      res[i, `:=`(n_events_D46 = 0L, n_events_D712 = 0L, n_events_D12p = 0L,
                  n_D4p = 0L, mean_I = NA_real_, mean_S = NA_real_)]
      next
    }
    
    # -------------------------------------------------------------------
    # Walk the series tracking drought runs.
    # Entry  : x  < threshold  (Sheffield & Wood single threshold)
    # Exit   : x  >= threshold (same threshold — S&W single-threshold rule)
    # All months inside a run are guaranteed x < threshold → deficit > 0.
    # -------------------------------------------------------------------
    in_d    <- FALSE
    s_idx   <- NA_integer_
    durs    <- integer(0)
    I_vec   <- numeric(0)   # per-event S&W intensity
    S_vec   <- numeric(0)   # per-event S&W severity
    
    for (j in seq_along(x)) {
      v  <- x[j]
      if (!in_d  && !is.na(v)  && v  < threshold) {
        in_d   <- TRUE
        s_idx  <- j
      } else if (in_d  && !is.na(v)  && v  >= threshold) {
        dur     <- j - s_idx                               # months in event
        I_ev    <- mean(threshold - x[s_idx:(j - 1L)])    # mean deficit
        S_ev    <- I_ev * dur                              # integrated deficit
        durs    <- c(durs,  dur)
        I_vec   <- c(I_vec, I_ev)
        S_vec   <- c(S_vec, S_ev)
        in_d   <- FALSE
        s_idx  <- NA_integer_
      }
    }
    # Close any event still open at end of record
    if (in_d  && !is.na(s_idx)) {
      dur    <- length(x) - s_idx + 1L
      I_ev   <- mean(threshold - x[s_idx:length(x)])
      S_ev   <- I_ev * dur
      durs   <- c(durs,  dur)
      I_vec  <- c(I_vec, I_ev)
      S_vec  <- c(S_vec, S_ev)
    }
    
    if (length(durs) == 0L) {
      res[i, `:=`(n_events_D46 = 0L, n_events_D712 = 0L, n_events_D12p = 0L,
                  n_D4p = 0L, mean_I = NA_real_, mean_S = NA_real_)]
      next
    }
    
    res[i, `:=`(
      n_events_D46  = sum(durs  >= 4L  & durs  <= 6L),
      n_events_D712 = sum(durs  >= 7L  & durs  <= 12L),
      n_events_D12p = sum(durs  >= 13L),
      n_D4p         = sum(durs  >= 4L),
      mean_I        = mean(I_vec, na.rm = TRUE),
      mean_S        = mean(S_vec, na.rm = TRUE)
    )]
  }
  res
}

# ================================================================================
#   FUNCTION 7 — ALTERNATIVE DROUGHT EVENT DETECTION
# ================================================================================
# Detects drought events using alternative thresholds:
# Onset: -1.0, Termination: -0.5, Scale-specific minimum durations
# Stores events in ALT_EVENTS_MEMORY for later reuse
vectorized_event_detection_alt <- function(ts_matrix, index_type, scale, 
                                           config = ALT_DROUGHT_DEF) {
  n_pix <- nrow(ts_matrix)
  
  # Determine min duration for this scale
  scale_char <- as.character(scale)
  min_dur <- if (!is.null(config$min_duration_map[[scale_char]])) {
    config$min_duration_map[[scale_char]]
  } else {
    config$default_min_duration
  }
  
  onset_thr <- config$onset_threshold
  term_thr  <- config$termination_threshold
  
  # Store events for this index-scale combination
  events_list <- list()
  
  for (i in seq_len(n_pix)) {
    x <- ts_matrix[i, ]
    x <- x[!is.na(x)]
    if (length(x) < 4L) next
    
    in_drought <- FALSE
    start_idx <- NA_integer_
    
    for (j in seq_along(x)) {
      v <- x[j]
      
      if (!in_drought && v < onset_thr) {
        # Drought Onset
        in_drought <- TRUE
        start_idx <- j
      } else if (in_drought && v >= term_thr) {
        # Drought Termination
        end_idx <- j - 1L
        duration <- end_idx - start_idx + 1L
        
        if (duration >= min_dur) {
          # Valid Event - store with pixel location
          events_list[[length(events_list) + 1]] <- list(
            pixel_id = i,
            start_idx = start_idx,
            end_idx = end_idx,
            duration = duration,
            intensity = mean(onset_thr - x[start_idx:end_idx]),
            severity = mean(onset_thr - x[start_idx:end_idx]) * duration
          )
        }
        in_drought <- FALSE
        start_idx <- NA_integer_
      }
    }
    
    # Handle case where drought continues to end of record
    if (in_drought && !is.na(start_idx)) {
      end_idx <- length(x)
      duration <- end_idx - start_idx + 1L
      if (duration >= min_dur) {
        events_list[[length(events_list) + 1]] <- list(
          pixel_id = i,
          start_idx = start_idx,
          end_idx = end_idx,
          duration = duration,
          intensity = mean(onset_thr - x[start_idx:end_idx]),
          severity = mean(onset_thr - x[start_idx:end_idx]) * duration
        )
      }
    }
  }
  
  # Store in global memory if events found
  if (length(events_list) > 0) {
    key <- sprintf("%s_%02d", index_type, scale)
    ALT_EVENTS_MEMORY[[key]] <<- list(
      events = events_list,
      n_events = length(events_list),
      config = config,
      index_type = index_type,
      scale = scale
    )
  }
  
  # Return summary statistics per pixel
  res <- data.table::data.table(
    alt_n_events      = rep(NA_integer_, n_pix),
    alt_mean_duration = rep(NA_real_,    n_pix),
    alt_max_severity  = rep(NA_real_,    n_pix)
  )
  
  for (evt in events_list) {
    i <- evt$pixel_id
    if (is.na(res[i, alt_n_events])) {
      res[i, `:=`(alt_n_events = 0L, 
                  alt_mean_duration = 0, 
                  alt_max_severity = 0)]
    }
    res[i, alt_n_events := alt_n_events + 1L]
    res[i, alt_mean_duration := alt_mean_duration + evt$duration]
    res[i, alt_max_severity := max(alt_max_severity, evt$severity, na.rm = TRUE)]
  }
  
  # Calculate mean duration for pixels with events
  has_events <- !is.na(res$alt_n_events) & res$alt_n_events > 0
  res[has_events, alt_mean_duration := alt_mean_duration / alt_n_events]
  
  return(res)
}

# ================================================================================
#   NEW FUNCTION 8 — Basin-level % area in drought time series + MK trend
# ================================================================================
# Computes a monthly basin-wide time series of the fraction of pixels in drought
# (index value < DROUGHT_ONSET) and runs Mann-Kendall on the annual means.
# Also aggregates events by duration class and runs MK on their frequency and
# mean severity over sliding 30-year windows.
compute_basin_extent_and_class_trends  <- function(ts_matrix, n_years, index_label,
                                                   out_dir,
                                                   onset_thr  = DROUGHT_ONSET,
                                                   end_thr    = DROUGHT_END,
                                                   start_year = 1950L,
                                                   win_years  = 30L) {
  utils_load_packages(c("ggplot2", "Kendall", "dplyr", "patchwork"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  n_months   <- n_years * 12L
  month_seq  <- seq_len(n_months)
  yr_seq     <- start_year + floor((month_seq - 1L) / 12L)
  mo_seq     <- ((month_seq - 1L) %% 12L) + 1L
  date_seq   <- as.Date(sprintf("%04d-%02d-01", yr_seq, mo_seq))
  
  #------------------------------------------------------------------------------
  # 7a. Monthly % area in drought
  #------------------------------------------------------------------------------
  n_pix    <- nrow(ts_matrix)
  n_col    <- min(n_months, ncol(ts_matrix))
  pct_area <- numeric(n_col)
  for (t in seq_len(n_col)) {
    col_vals <- ts_matrix[, t]
    ok       <- !is.na(col_vals)
    if (!any(ok)) {
      pct_area[t] <- NA
      next
    }
    pct_area[t] <- 100.0 * sum(ok & col_vals < onset_thr, na.rm = TRUE) / sum(ok)
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
  if (length(x_clean)  >= 10L) {
    mk_r   <- Kendall::MannKendall(x_clean)
    n      <- length(x_clean)
    pairs  <- outer(seq_len(n), seq_len(n),
                    function(j, ii) ifelse(j  > ii, (x_clean[j] - x_clean[ii]) / (j - ii), NA))
    mk_ext  <- list(
      tau        = as.numeric(mk_r$tau),
      p_value    = as.numeric(mk_r$sl),
      sens_slope = median(pairs, na.rm = TRUE)
    )
  }
  mk_ext_df  <- data.frame(
    index       = index_label,
    variable    = "mean_annual_pct_area_drought",
    tau         = mk_ext$tau,
    p_value     = mk_ext$p_value,
    sens_slope  = mk_ext$sens_slope,
    direction   = ifelse(!is.na(mk_ext$tau),
                         ifelse(mk_ext$tau  > 0, "Increasing", "Decreasing"), NA),
    significant = ifelse(!is.na(mk_ext$p_value), mk_ext$p_value  < 0.05, NA),
    stringsAsFactors = FALSE
  )
  write.csv(mk_ext_df,
            file.path(out_dir, sprintf("%s_basin_extent_MK.csv", index_label)),
            row.names = FALSE)
  cat(sprintf("  [Extent MK | %s] tau=%.3f  p=%.3f  slope=%.4f %%/yr  %s\n",
              index_label,
              mk_ext$tau %||% NA, mk_ext$p_value %||% NA, mk_ext$sens_slope %||% NA,
              ifelse(!is.na(mk_ext$p_value) && mk_ext$p_value  < 0.05, "sig", "")))
  
  #------------------------------------------------------------------------------
  # 7c. Per-pixel basin-level event detection → duration class catalog
  #------------------------------------------------------------------------------
  # Build a single basin-averaged time series, then detect events on it.
  basin_avg  <- colMeans(ts_matrix[, seq_len(n_col)], na.rm = TRUE)
  basin_ts   <- data.frame(date = date_seq[seq_len(n_col)], value = basin_avg)
  classify_dur   <- function(d)
    dplyr::case_when(d  >= 4L  & d  <= 6L  ~ "D4-6 (Short-term)",
                     d  >= 7L  & d  <= 12L ~ "D7-12 (Medium-term)",
                     d  >= 13L           ~ "D12+ (Long-term)",
                     TRUE                ~ "D1-3 (Sub-threshold)")
  
  # Identify events from basin-average series
  vals       <- basin_avg
  in_d       <- FALSE
  s_idx      <- NA_integer_
  event_list  <- list()
  for (j in seq_along(vals)) {
    v  <- vals[j]
    if (!in_d  && !is.na(v)  && v  < onset_thr) {
      in_d    <- TRUE
      s_idx   <- j
    } else if (in_d  && !is.na(v)  && v  >= end_thr) {
      dur   <- j - s_idx
      event_list[[length(event_list) + 1L]]  <- data.frame(
        start_date      = date_seq[s_idx],
        end_date        = date_seq[j - 1L],
        start_year      = yr_seq[s_idx],
        duration_months = dur,
        duration_class  = classify_dur(dur),
        mean_severity   = mean(abs(pmin(vals[s_idx:(j-1)], 0)), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      in_d    <- FALSE
      s_idx   <- NA_integer_
    }
  }
  if (in_d  && !is.na(s_idx)) {
    dur  <- n_col - s_idx + 1L
    event_list[[length(event_list) + 1L]]  <- data.frame(
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
  # ★★★ MOVED HERE - BEFORE event_list check so PNG is always saved ★★★
  extent_annual$trend_line  <- NA_real_
  if (!is.na(mk_ext$sens_slope) && !is.na(mk_ext$tau)) {
    mid   <- median(seq_along(x_clean))
    ymed  <- median(x_clean, na.rm = TRUE)
    extent_annual$trend_line  <-
      ymed + mk_ext$sens_slope * (seq_along(extent_annual$mean_pct) - mid)
  }
  sig_label  <- if (!is.na(mk_ext$p_value) && mk_ext$p_value  < 0.05)
    sprintf("MK: tau=%.3f, p=%.3f* (significant)", mk_ext$tau, mk_ext$p_value)
  else if (!is.na(mk_ext$p_value))
    sprintf("MK: tau=%.3f, p=%.3f (not significant)", mk_ext$tau, mk_ext$p_value)
  else
    "MK: insufficient data"
  
  # Calculate 10-year breaks - ensure clean intervals for all index types
  x_min <- min(extent_annual$year, na.rm = TRUE)
  x_max <- max(extent_annual$year, na.rm = TRUE)
  first_break <- ceiling(x_min / 10) * 10
  last_break  <- floor(x_max / 10) * 10
  x_breaks    <- seq(first_break, last_break, by = 10)
  x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  if (length(x_breaks) == 0) {
    x_breaks <- c(x_min, x_max)
  }
  if (length(x_breaks) > 15) {
    step_size <- ceiling((x_max - x_min) / 10)
    step_size <- max(10, step_size)
    x_breaks  <- seq(first_break, last_break, by = step_size)
    x_breaks  <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]
  }
  p_ext   <- ggplot2::ggplot(extent_annual, ggplot2::aes(x = year, y = mean_pct)) +
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
      title     = sprintf("%s — Annual Mean Basin-Wide Drought Extent", index_label),
      subtitle = sig_label,
      x = "Year", y = "% Basin Area in Drought"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  # SAVE PNG HERE (before event_list check)
  ggplot2::ggsave(
    file.path(out_dir, sprintf("%s_basin_extent_timeseries.png", index_label)),
    p_ext, width = 12, height = 5, dpi = 150
  )
  cat(sprintf("  ✓ Saved: %s_basin_extent_timeseries.png\n", index_label))
  
  # ★★★ NOW check for events (class-level analysis requires events) ★★★
  if (length(event_list) == 0) {
    cat(sprintf("  [Class Trends | %s] No events detected - skipping class-level MK\n", index_label))
    return(invisible(list(extent_ts   = extent_ts,
                          mk_extent   = mk_ext_df,
                          events_df   = NULL,
                          mk_by_class = NULL)))
  }
  
  # Continue with class-level analysis only if events exist
  events_df <- do.call(rbind, event_list)
  events_df$index <- index_label
  write.csv(events_df,
            file.path(out_dir, sprintf("%s_basin_event_catalog.csv", index_label)),
            row.names = FALSE)
  
  #------------------------------------------------------------------------------
  # 7d. MK on frequency and mean severity aggregated by class
  #------------------------------------------------------------------------------
  # Use a 30-year sliding window (same as Sheffield 2008, Table 4)
  dc_classes  <- c("D4-6 (Short-term)", "D7-12 (Medium-term)", "D12+ (Long-term)")
  end_year    <- start_year + n_years - 1L
  mid_years   <- seq(start_year + win_years %/% 2L,
                     end_year - win_years %/% 2L, by = 1L)
  mk_class_rows  <- list()
  for (cl in dc_classes) {
    freq_vec  <- numeric(length(mid_years))
    sev_vec   <- numeric(length(mid_years))
    for (k in seq_along(mid_years)) {
      wy1   <- mid_years[k] - win_years %/% 2L
      wy2   <- wy1 + win_years - 1L
      sub   <- events_df[events_df$duration_class == cl &
                           events_df$start_year  >= wy1 &
                           events_df$start_year  <= wy2, ]
      freq_vec[k]   <- nrow(sub)
      sev_vec[k]    <- if (nrow(sub)  > 0) mean(sub$mean_severity, na.rm = TRUE) else 0
    }
    # MK on frequency
    f_clean   <- freq_vec[!is.na(freq_vec)]
    mk_f      <- list(tau = NA, p_value = NA, sens_slope = NA)
    
    if (length(f_clean)  >= 10L && var(f_clean)  > 1e-10) {
      mk_r   <- Kendall::MannKendall(f_clean)
      n      <- length(f_clean)
      pa     <- outer(seq_len(n), seq_len(n),
                      function(j, ii) ifelse(j  > ii, (f_clean[j] - f_clean[ii]) / (j - ii), NA))
      mk_f   <- list(tau        = as.numeric(mk_r$tau),
                     p_value    = as.numeric(mk_r$sl),
                     sens_slope = median(pa, na.rm = TRUE))
    }
    
    # MK on severity
    s_clean   <- sev_vec[!is.na(sev_vec)]
    mk_s      <- list(tau = NA, p_value = NA, sens_slope = NA)
    
    if (length(s_clean)  >= 10L && var(s_clean)  > 1e-10) {
      mk_r   <- Kendall::MannKendall(s_clean)
      n      <- length(s_clean)
      pa     <- outer(seq_len(n), seq_len(n),
                      function(j, ii) ifelse(j  > ii, (s_clean[j] - s_clean[ii]) / (j - ii), NA))
      mk_s   <- list(tau        = as.numeric(mk_r$tau),
                     p_value    = as.numeric(mk_r$sl),
                     sens_slope = median(pa, na.rm = TRUE))
    }
    
    mk_class_rows[[length(mk_class_rows) + 1L]]  <- data.frame(
      index           = index_label,
      duration_class  = cl,
      variable        = "frequency_per_30yr",
      tau             = mk_f$tau,
      p_value         = mk_f$p_value,
      sens_slope      = mk_f$sens_slope,
      direction       = ifelse(!is.na(mk_f$tau),
                               ifelse(mk_f$tau  > 0, "Increasing", "Decreasing"), NA),
      significant     = ifelse(!is.na(mk_f$p_value), mk_f$p_value  < 0.05, NA),
      stringsAsFactors = FALSE
    )
    
    mk_class_rows[[length(mk_class_rows) + 1L]]  <- data.frame(
      index           = index_label,
      duration_class  = cl,
      variable        = "mean_severity",
      tau             = mk_s$tau,
      p_value         = mk_s$p_value,
      sens_slope      = mk_s$sens_slope,
      direction       = ifelse(!is.na(mk_s$tau),
                               ifelse(mk_s$tau  > 0, "Increasing", "Decreasing"), NA),
      significant     = ifelse(!is.na(mk_s$p_value), mk_s$p_value  < 0.05, NA),
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("  [Class MK | %s | %-25s] freq: tau=%.3f p=%.3f | sev: tau=%.3f p=%.3f\n",
                index_label, cl,
                mk_f$tau %||% NA, mk_f$p_value %||% NA,
                mk_s$tau %||% NA, mk_s$p_value %||% NA))
  }
  mk_class_df <- do.call(rbind, mk_class_rows)
  write.csv(mk_class_df,
            file.path(out_dir, sprintf("%s_class_frequency_severity_MK.csv", index_label)),
            row.names = FALSE)
  invisible(list(extent_ts   = extent_ts,
                 mk_extent   = mk_ext_df,
                 events_df   = events_df,
                 mk_by_class = mk_class_df))
}



# ================================================================================
#   MAIN PROCESSING FUNCTION
# ================================================================================
process_index <- function(index_type, scales, find_fn, seas_dir) {
  for (sc in scales) {
    cat(sprintf("\n╔══ %s-%02d ══╗\n", toupper(index_type), sc))
    
    # [1/5] Find monthly NC files
    cat("[1/5] Finding files...\n")
    files   <- find_fn(seas_dir, index_type, sc)   # unified 3-arg call
    
    if (!length(files)) {
      cat("  ❌ No files found\n")
      next
    }
    
    cat(sprintf("  Found %d monthly NC files\n", length(files)))
    
    # [2/5] Load raster stacks
    cat("[2/5] Loading rasters...\n")
    stacks   <- lapply(files, function(f) {
      tryCatch(terra::rast(f),
               error = function(e) {
                 cat("  ⚠ Failed: ", basename(f), "\n")
                 NULL
               })
    })
    stacks   <- stacks[!sapply(stacks, is.null)]
    
    if (!length(stacks)) {
      cat("  ❌ All rasters failed to load\n")
      next
    }
    
    # Determine years from first file's embedded dates
    n_lyrs_per_file   <- terra::nlyr(stacks[[1]])
    dates_m1          <- extract_dates_from_nc(files[1], n_lyrs_per_file)
    years             <- as.integer(format(dates_m1, "%Y"))
    n_years           <- length(years)
    
    # Pixel coordinates
    r_tmpl    <- stacks[[1]][[1]]
    xy        <- terra::xyFromCell(r_tmpl, seq_len(terra::ncell(r_tmpl)))
    lon_vec   <- xy[, 1]
    lat_vec   <- xy[, 2]
    n_pix     <- length(lon_vec)
    
    cat(sprintf("  Grid: %d pixels × %d years × 12 months\n", n_pix, n_years))
    
    # [3/5] Build time-series matrix (pixels × months)
    # Each NC file contains all years for one calendar month.
    # Interleave: Jan_yr1, Feb_yr1, …, Dec_yr1, Jan_yr2, …
    cat("[3/5] Building interleaved time-series matrix...\n")
    ts_mat   <- matrix(NA_real_, nrow = n_pix, ncol = n_years * 12)
    
    for (m in seq_along(stacks)) {
      n_yrs_m   <- min(terra::nlyr(stacks[[m]]), n_years)
      for (y in seq_len(n_yrs_m)) {
        col_idx             <- (y - 1) * 12 + m
        ts_mat[, col_idx]   <- as.numeric(terra::values(stacks[[m]][[y]]))
      }
    }
    
    # Filter to pixels with at least some valid data
    valid_pix   <- rowSums(!is.na(ts_mat))  > 0
    ts_mat      <- ts_mat[valid_pix, , drop = FALSE]
    lon_vec     <- lon_vec[valid_pix]
    lat_vec     <- lat_vec[valid_pix]
    n_pix       <- sum(valid_pix)
    
    cat(sprintf("  Valid pixels: %d  |  %.1f%% non-NA values\n",
                n_pix,
                100 * sum(!is.na(ts_mat)) / length(ts_mat)))
    
    # Fractional year vector (for regime-shift detection)
    yr_seq   <- rep(years, each = 12) + (rep(0:11, times = n_years) / 12)
    
    # [4/5] Compute all statistics
    cat("[4/5] Computing statistics...\n")
    results   <- data.table::data.table(lon = lon_vec, lat = lat_vec)
    
    cat("  → Mann-Kendall (variance-corrected)...    ")
    results   <- cbind(results, vectorized_mann_kendall(ts_mat))
    cat("✓\n")
    
    cat("  → Mann-Kendall (TFPW)...    ")
    results   <- cbind(results, vectorized_mann_kendall_tfpw(ts_mat))
    cat("✓\n")
    
    cat("  → Drought event detection...    ")
    results   <- cbind(results, vectorized_event_detection(ts_mat))
    cat("✓\n")
    
    cat("  → CUSUM regime-shift detection...    ")
    results$regime_shift_year  <- vectorized_regime_shift(ts_mat, yr_seq)
    cat("✓\n")
    
    cat("  → Runs test (temporal clustering)...    ")
    results   <- cbind(results, vectorized_runs_test(ts_mat))
    cat("✓\n")
    
    # Spectral peaks – placeholder; reserved for future implementation
    results$n_spectral_peaks  <- 0L
    
    # NEW: Per-pixel duration-class event counts
    cat("  → Duration-class event counts (D4-6 / D7-12 / D12+)...    ")
    results   <- cbind(results, vectorized_duration_class_counts(ts_mat))
    cat("✓\n")
    
    cat(sprintf("     D4-6 events  — mean per pixel: %.2f\n",
                mean(results$n_events_D46,  na.rm = TRUE)))
    cat(sprintf("     D7-12 events — mean per pixel: %.2f\n",
                mean(results$n_events_D712, na.rm = TRUE)))
    cat(sprintf("     D12+ events  — mean per pixel: %.2f  (max: %d)\n",
                mean(results$n_events_D12p, na.rm = TRUE),
                max(results$n_events_D12p,  na.rm = TRUE)))
    cat(sprintf("     Mean pixel intensity (mean_I): %.4f\n",
                mean(results$mean_I, na.rm = TRUE)))
    cat(sprintf("     Mean pixel severity  (mean_S): %.4f\n",
                mean(results$mean_S, na.rm = TRUE)))
    
    # [NEW] Alternative drought event detection (Onset -1.0, Term -0.5)
    cat("  → Alternative drought events (Onset -1.0/Term -0.5)...    ")
    alt_results <- vectorized_event_detection_alt(ts_mat, index_type, sc)
    results <- cbind(results, alt_results)
    cat("✓ (stored in ALT_EVENTS_MEMORY)\n")
    
    # NEW: Basin-level % area in drought time series + class MK
    cat("  → Basin extent time series   & class-level MK trends...    ")
    index_label   <- sprintf("%s-%02d", toupper(index_type), sc)
    extent_out    <- file.path(TREND_DIR, "basin_extent")
    
    basin_extent_results   <- tryCatch(
      compute_basin_extent_and_class_trends(
        ts_mat, n_years, index_label, extent_out,
        start_year = min(years)
      ),
      error = function(e) {
        cat(sprintf("\n  ⚠ Basin extent analysis failed: %s\n", e$message))
        NULL
      }
    )
    
    cat("✓\n")
    
    # [5/5] Append raw time-series columns and save
    cat("[5/5] Saving CSV...\n")
    
    for (t in seq_len(ncol(ts_mat)))
      results[[sprintf("time_%03d", t)]]  <- ts_mat[, t]
    
    out_file   <- file.path(TREND_DIR,
                            sprintf("%s_%02d_results.csv", index_type, sc))
    data.table::fwrite(results, out_file, showProgress = FALSE)
    
    cat(sprintf("  ✅ %s  (%.2f MB)\n",
                basename(out_file),
                file.info(out_file)$size / 1024^2))
    cat(sprintf("  Significant VC trends  (p < 0.05) : %d / %d pixels\n",
                sum(results$p_value_vc  < 0.05, na.rm = TRUE),
                sum(!is.na(results$p_value_vc))))
    cat(sprintf("  Significant TFPW trends (p < 0.05): %d / %d pixels\n",
                sum(results$p_value_tfpw  < 0.05, na.rm = TRUE),
                sum(!is.na(results$p_value_tfpw))))
    
    # NEW: Write a lightweight duration-class frequency + intensity/severity map CSV.
    dur_map_cols   <- c("lon", "lat",
                        "n_events_D46", "n_events_D712", "n_events_D12p", "n_D4p",
                        "mean_I", "mean_S")
    dur_map_file   <- file.path(TREND_DIR,
                                sprintf("%s_%02d_duration_class_map.csv",
                                        index_type, sc))
    data.table::fwrite(results[, ..dur_map_cols], dur_map_file,
                       showProgress = FALSE)
    cat(sprintf("  ✅ Duration-class map: %s\n", basename(dur_map_file)))
    
    # NEW: NetCDF raster of D12+ count (long-term drought frequency map)
    if (exists("r_tmpl") && !is.null(r_tmpl)) {
      r_D12p   <- r_tmpl
      terra::values(r_D12p)  <- NA_real_
      
      # Re-index: results rows correspond to the valid_pix subset;
      # we need to map back to all cells via the valid_pix logical vector.
      all_vals   <- rep(NA_real_, terra::ncell(r_tmpl))
      all_vals[valid_pix]  <- results$n_events_D12p
      terra::values(r_D12p)  <- all_vals
      
      names(r_D12p)  <- sprintf("%s_%02d_n_D12p_longterm_droughts",
                                index_type, sc)
      nc_D12p   <- file.path(TREND_DIR,
                             sprintf("%s_%02d_D12p_pixel_frequency.nc",
                                     index_type, sc))
      terra::writeCDF(r_D12p, nc_D12p, overwrite = TRUE)
      cat(sprintf("  ✅ D12+ frequency raster: %s\n", basename(nc_D12p)))
    }
    
    rm(ts_mat, stacks)
    invisible(gc())
  }
}

# ================================================================================
#   RUN ALL INDICES
# ================================================================================
cat("\n╔════════════════════════════════════╗\n")
cat("║  DROUGHT TREND TEST  (w3)          ║\n")
cat("╚════════════════════════════════════╝\n\n")
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

# SWEI wrapper: find_swei_seasonal_files takes (data_dir, scale) only;
# wrap it safely so a missing-file error skips SWEI without stopping the run.
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

# MSPI / MSPEI (optional — only runs if mspi_results/ folder exists)
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
#   Shut down parallel cluster (always, even after errors)
# --------------------------------------------------------------------------------
if (!is.null(cl)) {
  parallel::stopCluster(cl)
  cat("\n✓ Cluster stopped\n")
}

elapsed <- Sys.time() - total_start
cat(sprintf("\n✅ Done.  Time: %.1f min\n  Next: run w4_trends_visualization.R\n\n",
            as.numeric(elapsed, units = "mins")))

# --------------------------------------------------------------------------------
#   [NEW] FINAL SUMMARY: ALTERNATIVE EVENTS MEMORY
# --------------------------------------------------------------------------------
cat("══════════════════════════════════════\n")
cat("ALTERNATIVE DROUGHT EVENTS SUMMARY\n")
cat("══════════════════════════════════════\n")
cat("Configuration: Onset = -1.0, Termination = -0.5\n")
cat("Scale-specific minimum durations:\n")
cat("  • SPI/SPEI-1:  ≥ 2 months\n")
cat("  • SPI/SPEI-3:  ≥ 3 months\n")
cat("  • SPI/SPEI-6:  ≥ 4 months\n")
cat("  • SPI/SPEI-12: ≥ 6 months\n")
cat("\n")

if (length(ALT_EVENTS_MEMORY) > 0) {
  cat(sprintf("✓ Stored %d alternative event catalogs in ALT_EVENTS_MEMORY\n",
              length(ALT_EVENTS_MEMORY)))
  cat("  Access example: ALT_EVENTS_MEMORY[['spi_12']]\n")
  cat("  Structure: list(index_scale = list(events, n_events, config))\n")
  cat("\n  Contents:\n")
  for (key in names(ALT_EVENTS_MEMORY)) {
    cat(sprintf("    • %s: %d events detected\n",
                key, ALT_EVENTS_MEMORY[[key]]$n_events))
  }
} else {
  cat("⚠ No alternative events stored (check input data availability)\n")
}

cat("\n")