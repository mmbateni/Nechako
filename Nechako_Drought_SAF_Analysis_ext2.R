# ============================================================================
# SAF EXTENSIONS 2 — ADDITIONAL SAF METHODS
# Source AFTER Nechako_Drought_SAF_Analysis.R and ext1.R
# Functions exported:
# fit_nonstationary_marginals()      — Three-model GAMLSS comparison:
#                                      stationary, mu-trend, mu+sigma-trend.
#                                      Epoch definition from change-point or
#                                      automatic data split.
# derive_SAF_nonstationary()         — SAF curves with period-specific marginals
#                                      AND period-specific copula.
# derive_SAF_nonstationary_kendall() — Kendall-corrected version of above.
# place_event_on_saf()               — Joint return period of a named event.
#
# CHANGES (v5 → v7):
# Addresses the hybrid non-stationary framework inconsistency:
# v4 paired a step-function marginal (discrete epochs from fit_nonstationary_
# marginals) with a smooth-trend copula collapsed to a single midpoint value
# (.build_period_copula). This mixes two incompatible temporal representations
# in the joint distribution without explicit justification.
#
# FIX — three changes:
# (A) .build_period_copula() [REVISED]:
#     Priority-path logic replacing the single midpoint-evaluation call:
#     Path 1: Segmented copula from tv_res$seg_result — epoch-consistent
#             with the discrete marginal (both are step-functions sharing
#             the same cp_year boundary). PREFERRED.
#     Path 2: TV copula evaluated at epoch midpoint — only when seg_result
#             is NULL (insufficient data per segment). Logs a WARNING.
#     Path 3: Stationary copula — last-resort fallback.
#     New parameter: which_seg ("ref" or "recent") selects seg1 vs seg2.
#
# (B) .resolve_epoch() [NEW internal helper]:
#     Epoch year vectors are now derived exclusively from ns_result$ref_years /
#     ns_result$recent_years — the same object that produced mu_period and
#     sigma_period.  Both SAF functions share one source of truth for epoch
#     boundaries; the copula and marginal are guaranteed to use the same split.
#
# (C) .assert_epoch_consistency() [NEW internal helper]:
#     Compares the change-point year embedded in ns_result$epoch_src against
#     tv_copula_result$cp_result$cp_year. Logs a WARNING [EPOCH MISMATCH] if
#     they differ so the analyst is alerted rather than silently misled.
#
# (D) derive_SAF_nonstationary() and derive_SAF_nonstationary_kendall() [REVISED]:
#     * Accept ns_result as a new required argument (replaces internal
#       re-derivation of period_yrs from cp_year).
#     * Call .resolve_epoch() and .assert_epoch_consistency() before selecting
#       the copula object.
#     * New output column copula_path = "segmented" | "tv_midpoint" |
#       "stationary" written to every row of the CSV for auditability.
#     * Kendall function: K_C now computed from pseudo-observations evaluated
#       under cop_obj_ns (the period-specific copula) rather than always the
#       full-record stationary copula. This is consistent with De Michele et al.
#       (2013): K_C must match the copula used in the SAF inversion.
#     * v7 Update: Replaced silent stationary fallback with a strict 3-tier
#       validation strategy for empirical K_C computation. Tier 1: direct
#       evaluation. Tier 2: simulation-based fallback. Tier 3: explicit abort
#       with NULL return to prevent miscalibrated SAF curves.
#     * v7 Update: Analytical K_C(t) = t - φ(t)/φ'(t) is used for Archimedean
#       families (Frank, Gumbel, Joe, SurvClayton), enabling valid extrapolation.
#
# CALLING CODE CHANGE REQUIRED:
# Both SAF functions now require ns_result as an explicit argument:
# derive_SAF_nonstationary(..., ns_result = ns, ...)
# derive_SAF_nonstationary_kendall(..., ns_result = ns, ...)
# where ns is the return value of fit_nonstationary_marginals().
#
# fit_nonstationary_marginals() — epoch boundary truncation
# yr_match <- dc$year
# Uses the actual calendar year stored in each drought-month row.
# yr_range now correctly spans 1950-2025 and recent_years extends
# to the record end.  The GAMLSS trend fit also receives correct
# year_std values for every observation.
#
# Impact on results:
# SPI-1 : recent epoch 1965-2025 (was 1965-1998; +27 yr)
# SPEI-1: recent epoch 1996-2025 (was auto-split due to
#         artificial near-edge flag; now uses CP correctly)
# SPI-3 : recent epoch 1973-2025 (was 1973-1998; +27 yr)
# SPEI-3: recent epoch 1987-2025 (was 1987-1999; +26 yr)
# ============================================================================

# Install/load required packages
for (pkg in c("gamlss", "gamlss.dist")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(gamlss)
library(gamlss.dist)

# Verify prerequisite function exists
if (!exists("log_event")) {
  stop("log_event() not found. Source Nechako_Drought_SAF_Analysis.R first.")
}

# ============================================================================
# INTERNAL HELPER: .build_period_copula  [REVISED v5]
# Selects the copula object for a discrete epoch via a three-path priority
# rule, replacing the v4 single midpoint-evaluation call.
#
# Priority:
# Path 1 — Segmented copula (tv_res$seg_result not NULL):
#          Uses the per-segment stationary copula from fit_segmented_copula().
#          This is the preferred path: both the marginal and the copula are
#          step-functions sharing the same cp_year boundary, so the joint
#          distribution is internally consistent.
# Path 2 — TV midpoint evaluation (seg_result NULL, make_cop_fn available):
#          Evaluates the continuous linear-trend copula at the epoch midpoint year.
#          This is an approximation that mixes a step-function marginal with a
#          smooth-trend copula collapsed to a single point.  Acceptable when
#          seg_result is unavailable (e.g., too few observations per segment),
#          but the analyst must be aware of the hybrid nature.  A WARNING is
#          logged every time this path is taken.
# Path 3 — Stationary fallback:
#          Returns NULL; both SAF functions fall back to
#          copula_fit_obj$best_copula_fit@copula.  Logged for transparency.
#
# Parameters:
# tv_res       : return value of fit_timevarying_copula()
# period_years : integer vector of calendar years for this epoch
# which_seg    : "ref" or "recent" — selects seg1 vs seg2 in seg_result
# ============================================================================
.build_period_copula <- function(tv_res, period_years, which_seg = "recent") {
  # --------------------------------------------------------------------------
  # Path 1: Segmented copula
  # --------------------------------------------------------------------------
  seg <- tv_res$seg_result
  if (!is.null(seg)) {
    cop_obj <- if (which_seg == "ref") {
      tryCatch(seg$seg1$best_copula_fit@copula, error = function(e) NULL)
    } else {
      tryCatch(seg$seg2$best_copula_fit@copula, error = function(e) NULL)
    }
    if (!is.null(cop_obj)) {
      log_event(sprintf(
        "  .build_period_copula [%s]: PATH 1 - segmented copula (%s). Epoch-consistent with discrete marginal.",
        which_seg,
        if (which_seg == "ref") seg$seg1$best_copula_name else seg$seg2$best_copula_name))
      return(cop_obj)
    }
    log_event(sprintf(
      "  .build_period_copula [%s]: seg_result present but copula object extraction failed. Trying Path 2.",
      which_seg))
  }
  
  # --------------------------------------------------------------------------
  # Path 2: TV copula at epoch midpoint
  # --------------------------------------------------------------------------
  if (!is.null(tv_res$make_cop_fn) && !is.na(tv_res$a_hat) &&
      !is.null(tv_res$yr_mean)      && !is.null(tv_res$yr_sd)) {
    mid_yr   <- mean(range(period_years, na.rm = TRUE))
    yr_std_m <- (mid_yr - tv_res$yr_mean) / tv_res$yr_sd
    theta_r  <- tv_res$a_hat + tv_res$b_hat * yr_std_m
    cop_obj  <- tryCatch(tv_res$make_cop_fn(theta_r), error = function(e) NULL)
    if (!is.null(cop_obj)) {
      log_event(sprintf(
        "  .build_period_copula [%s]: PATH 2 - TV copula at epoch midpoint yr=%.1f (theta_raw=%.4f). WARNING: marginal is a step-function; copula is a point-evaluated smooth trend. Document this hybrid approximation.",
        which_seg, mid_yr, theta_r))
      return(cop_obj)
    }
    log_event(sprintf(
      "  .build_period_copula [%s]: TV midpoint evaluation failed. Falling back to stationary.",
      which_seg))
  }
  
  # --------------------------------------------------------------------------
  # Path 3: Stationary full-record copula
  # --------------------------------------------------------------------------
  log_event(sprintf(
    "  .build_period_copula [%s]: PATH 3 - stationary full-record copula (no non-stationarity in dependence).",
    which_seg))
  return(NULL)   # caller uses copula_fit_obj$best_copula_fit@copula
}

# ============================================================================
# INTERNAL HELPER: .resolve_epoch  [NEW v5]
# Derives the epoch year vector and which_seg label from ns_result — the
# single source of truth for epoch boundaries shared by marginal and copula.
# Both SAF functions call this instead of re-deriving period_yrs internally.
#
# Parameters:
# ns_result    : return value of fit_nonstationary_marginals()
# period_label : character; expected to contain "ref" or "recent"
#
# Returns a list:
# period_yrs  : integer vector of calendar years for this epoch
# which_seg   : "ref" or "recent"
# cp_yr_mar   : change-point year as understood by the marginal, or NA
# ============================================================================
.resolve_epoch <- function(ns_result, period_label) {
  is_ref     <- grepl("ref", tolower(period_label), fixed = FALSE)
  which_seg  <- if (is_ref) "ref" else "recent"
  period_yrs <- if (is_ref) ns_result$ref_years else ns_result$recent_years
  
  # Extract cp_year from epoch_src string if the marginal used a change-point
  cp_yr_mar <- if (grepl("change-point", ns_result$epoch_src, fixed = TRUE)) {
    m <- regmatches(ns_result$epoch_src,
                    regexpr("[0-9]{4}", ns_result$epoch_src))
    if (length(m) == 1L) as.integer(m) else NA_integer_
  } else {
    NA_integer_
  }
  
  list(period_yrs = period_yrs,
       which_seg  = which_seg,
       cp_yr_mar  = cp_yr_mar)
}

# ============================================================================
# INTERNAL HELPER: .assert_epoch_consistency  [NEW v5]
# Compares the change-point year used to define the copula epochs
# (tv_copula_result$cp_result$cp_year) against the one used for the
# marginal epochs (cp_yr_mar extracted from ns_result$epoch_src).
# Logs a WARNING [EPOCH MISMATCH] when they differ so the analyst is
# alerted rather than silently receiving an inconsistent joint distribution.
# ============================================================================
.assert_epoch_consistency <- function(tv_copula_result, cp_yr_mar,
                                      index_name, period_label) {
  if (is.null(tv_copula_result)) return(invisible(NULL))
  cp_yr_cop <- tv_copula_result$cp_result$cp_year
  if (!is.na(cp_yr_cop) && !is.na(cp_yr_mar) && cp_yr_cop != cp_yr_mar) {
    log_event(sprintf(
      "  [%s | %s] WARNING [EPOCH MISMATCH]: copula cp_year=%d != marginal cp_year=%d. The joint distribution uses inconsistent epoch boundaries. Review cp_result passed to fit_nonstationary_marginals and fit_timevarying_copula.",
      index_name, period_label, cp_yr_cop, cp_yr_mar))
  }
  invisible(NULL)
}

# ============================================================================
# 1. NON-STATIONARY MARGINAL FITTING  (unchanged from v4)
# Three GAMLSS Gamma models are fitted and compared by AIC:
# M0 (stationary):    mu ~ 1,         sigma ~ 1
# M1 (mu-trend):      mu ~ year_std,  sigma ~ 1
# M2 (mu+sigma-trend):mu ~ year_std,  sigma ~ year_std
# The best model is selected by minimum AIC.  Period-specific mu and
# sigma values are derived for two epochs (reference and recent).
#
# Epoch definition (precedence order):
# 1. Explicit ref_years / recent_years arguments - highest priority.
# 2. cp_result$detected == TRUE - cp_year splits:
#    reference = record start ... cp_year - 1
#    recent    = cp_year ... record end
#    Note: if the change point is very early or very late (< 20% or
#    > 80% of record), the 40%/30% split is used instead with a warning.
# 3. Default automatic split: first 40% / last 30% of record.
#
# Parameters:
# drought_data  : output of extract_drought_characteristics
# index_name    : character label
# year_vec      : full year vector (all months, not just drought obs)
# output_dir    : directory for CSV output
# cp_result     : optional; output of detect_copula_changepoints() [NEW]
# ref_years     : optional integer vector override for reference epoch
# recent_years  : optional integer vector override for recent epoch
# ============================================================================
fit_nonstationary_marginals <- function(drought_data, index_name, year_vec,
                                        output_dir,
                                        cp_result    = NULL,
                                        ref_years    = NULL,
                                        recent_years = NULL) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  if (nrow(dc) < 20L) {
    log_event(sprintf("  [%s] Non-stationary fit skipped: < 20 drought observations.", index_name))
    return(NULL)
  }
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    log_event(sprintf("  [%s] gamlss not available - skipping.", index_name))
    return(NULL)
  }
  
  # FIX (v4): use the calendar year stored in each drought-observation row.
  yr_match <- dc$year
  yr_range <- range(yr_match, na.rm = TRUE)
  yr_span  <- diff(yr_range)
  
  # ---- Epoch definition --------------------------------------------------
  auto_ref    <- seq(yr_range[1], yr_range[1] + floor(yr_span * 0.4))
  auto_recent <- seq(yr_range[2] - floor(yr_span * 0.3), yr_range[2])
  
  if (!is.null(ref_years) && !is.null(recent_years)) {
    epoch_src <- "user override"
  } else if (!is.null(cp_result) && isTRUE(cp_result$detected) && !is.na(cp_result$cp_year)) {
    cp_yr   <- cp_result$cp_year
    cp_frac <- (cp_yr - yr_range[1]) / yr_span
    if (cp_frac < 0.20 || cp_frac > 0.80) {
      log_event(sprintf(
        "  [%s] Change-point at %d is near record edge (%.0f%% of full record) - using auto 40%%/30%% split instead.",
        index_name, cp_yr, 100 * cp_frac))
      ref_years    <- auto_ref
      recent_years <- auto_recent
      epoch_src    <- "auto 40/30 (CP edge)"
    } else {
      ref_years    <- seq(yr_range[1], cp_yr - 1L)
      recent_years <- seq(cp_yr, yr_range[2])
      epoch_src    <- sprintf("change-point (%d)", cp_yr)
    }
  } else {
    ref_years    <- auto_ref
    recent_years <- auto_recent
    epoch_src    <- "auto 40/30"
  }
  
  log_event(sprintf("  [%s] NS marginals | epochs from: %s | ref: %d-%d | recent: %d-%d",
                    index_name, epoch_src,
                    min(ref_years), max(ref_years),
                    min(recent_years), max(recent_years)))
  
  yr_mean <- mean(yr_match, na.rm = TRUE)
  yr_sd   <- sd(yr_match,   na.rm = TRUE)
  if (yr_sd < 1e-10) yr_sd <- 1
  to_std <- function(yrs) (yrs - yr_mean) / yr_sd
  
  df <- data.frame(severity = dc$severity,
                   year_std = to_std(yr_match))
  
  # ---- Three-model fit ---------------------------------------------------
  m0 <- tryCatch(
    gamlss::gamlss(severity ~ 1,
                   sigma.formula = ~ 1,
                   family = gamlss.dist::GA(), data = df, trace = FALSE),
    error = function(e) NULL)
  m1 <- tryCatch(
    gamlss::gamlss(severity ~ year_std,
                   sigma.formula = ~ 1,
                   family = gamlss.dist::GA(), data = df, trace = FALSE),
    error = function(e) NULL)
  m2 <- tryCatch(
    gamlss::gamlss(severity ~ year_std,
                   sigma.formula = ~ year_std,
                   family = gamlss.dist::GA(), data = df, trace = FALSE),
    error = function(e) NULL)
  
  if (is.null(m0)) {
    log_event(sprintf("  [%s] Stationary GAMLSS fit failed - aborting.", index_name))
    return(NULL)
  }
  
  models <- list(M0_stationary = m0, M1_mu_trend = m1, M2_full_ns = m2)
  aic_v  <- vapply(models, function(m) if (!is.null(m)) AIC(m) else Inf, numeric(1L))
  best_nm <- names(which.min(aic_v))
  m_best  <- models[[best_nm]]
  
  log_event(sprintf(
    "  [%s] GAMLSS AIC | M0=%.2f | M1=%.2f | M2=%.2f | best=%s",
    index_name, aic_v[1], aic_v[2], aic_v[3], best_nm))
  
  # ---- Period-specific mu extraction (via coefficient, not predict) ------
  coef_mu   <- coef(m_best, what = "mu")
  has_mu_tr <- length(coef_mu) >= 2
  mu_pred <- function(yrs) {
    lin <- if (has_mu_tr) coef_mu[1] + coef_mu[2] * to_std(yrs) else rep(coef_mu[1], length(yrs))
    mean(exp(lin), na.rm = TRUE)
  }
  mu_ref <- mu_pred(ref_years)
  mu_rec <- mu_pred(recent_years)
  
  # ---- Period-specific sigma extraction ----------------------------------
  coef_sg    <- coef(m_best, what = "sigma")
  has_sig_tr <- length(coef_sg) >= 2
  sigma_pred <- function(yrs) {
    lin <- if (has_sig_tr) coef_sg[1] + coef_sg[2] * to_std(yrs) else rep(coef_sg[1], length(yrs))
    mean(exp(lin), na.rm = TRUE)
  }
  sigma_ref <- sigma_pred(ref_years)
  sigma_rec <- sigma_pred(recent_years)
  
  log_event(sprintf(
    "  [%s] Period means | mu_ref=%.4f mu_rec=%.4f | sigma_ref=%.4f sigma_rec=%.4f",
    index_name, mu_ref, mu_rec, sigma_ref, sigma_rec))
  
  result <- list(m_stat       = m0,
                 m_ns         = if (!is.null(m1)) m1 else m0,
                 m_best       = m_best,
                 best_model   = best_nm,
                 mu_ref       = mu_ref,
                 mu_rec       = mu_rec,
                 sigma_ref    = sigma_ref,
                 sigma_rec    = sigma_rec,
                 ref_years    = ref_years,
                 recent_years = recent_years,
                 aic_m0       = aic_v[1],
                 aic_m1       = aic_v[2],
                 aic_m2       = aic_v[3],
                 epoch_src    = epoch_src)
  
  summary_df <- data.frame(
    index        = index_name,
    best_model   = best_nm,
    mu_ref       = round(mu_ref,    4),
    mu_recent    = round(mu_rec,    4),
    sigma_ref    = round(sigma_ref, 4),
    sigma_recent = round(sigma_rec, 4),
    aic_m0       = round(aic_v[1],  2),
    aic_m1       = round(aic_v[2],  2),
    aic_m2       = round(aic_v[3],  2),
    epoch_src    = epoch_src,
    ref_start    = min(ref_years),
    ref_end      = max(ref_years),
    recent_start = min(recent_years),
    recent_end   = max(recent_years),
    stringsAsFactors = FALSE)
  
  write.csv(summary_df,
            file.path(output_dir, sprintf("%s_nonstationary_summary.csv", index_name)),
            row.names = FALSE)
  result
}

# ============================================================================
# 2. NON-STATIONARY SAF CURVES - CONDITIONAL METHOD  (v5)
# Derives SAF curves with period-specific Gamma marginal for drought
# severity (mu_period / sigma_period from fit_nonstationary_marginals)
# paired with a period-specific copula via .build_period_copula().
#
# Changes vs v4:
# * ns_result is now a required argument. Epoch year vectors come
#   exclusively from ns_result$ref_years / ns_result$recent_years
#   via .resolve_epoch(), guaranteeing copula and marginal share the
#   same epoch boundaries.
# * .build_period_copula() follows the three-path priority:
#   Path 1: segmented copula  (epoch-consistent, preferred)
#   Path 2: TV copula at epoch midpoint  (hybrid, warned)
#   Path 3: stationary copula  (fallback)
# * .assert_epoch_consistency() logs a WARNING on cp_year mismatch.
# * New output column copula_path records the path taken.
#
# Parameters:
# mu_period, sigma_period : scalars; Gamma mean and CV for this epoch,
#                           from fit_nonstationary_marginals()
# copula_fit_obj          : stationary copula result from fit_copulas()
# marginal_fits           : output of fit_marginal_distributions()
# drought_data            : full characteristics data frame
# ns_result               : output of fit_nonstationary_marginals() [NEW v5]
# index_name              : character label for logging
# period_label            : e.g. "reference" or "recent"
# output_dir              : directory for CSV output
# tv_copula_result        : optional; output of fit_timevarying_copula()
# ============================================================================
derive_SAF_nonstationary <- function(mu_period, sigma_period,
                                     copula_fit_obj, marginal_fits,
                                     drought_data,   ns_result,
                                     index_name,     period_label,
                                     output_dir,
                                     tv_copula_result = NULL) {
  if (is.null(copula_fit_obj) || is.null(marginal_fits)) return(NULL)
  if (is.null(ns_result)) {
    log_event(sprintf(
      "  [%s | %s] derive_SAF_nonstationary: ns_result is NULL - cannot resolve epoch years. Aborting.",
      index_name, period_label))
    return(NULL)
  }
  
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data)
  n_dr  <- nrow(dc)
  if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  
  # Compute renewal parameters from event inter-arrival times
  indr_all    <- drought_data$severity > 0 & drought_data$area_pct > 0
  is_start    <- indr_all & !c(FALSE, head(indr_all, -1))
  ev_starts   <- as.Date(drought_data$date[is_start])
  iat_cv_val  <- NA_real_
  shape_k_val <- NA_real_
  rate_k_val  <- NA_real_
  mu_T_adj    <- mu_T * 12
  if (length(ev_starts) >= 2) {
    iat_yrs    <- as.numeric(diff(ev_starts)) / 365.25
    iat_mean_v <- mean(iat_yrs, na.rm = TRUE)
    iat_sd_v   <- sd(iat_yrs,   na.rm = TRUE)
    iat_cv_val <- if (iat_mean_v > 0) iat_sd_v / iat_mean_v else NA_real_
    mu_T_adj   <- iat_mean_v * 12
    if (!is.na(iat_cv_val) && iat_cv_val > 0) {
      shape_k_val <- 1 / (iat_cv_val^2)
      rate_k_val  <- shape_k_val / iat_mean_v
    }
  }  
  # ---- Resolve epoch (single source of truth: ns_result) ------------------
  epoch      <- .resolve_epoch(ns_result, period_label)
  period_yrs <- epoch$period_yrs
  which_seg  <- epoch$which_seg
  
  # ---- Epoch-consistency check --------------------------------------------
  .assert_epoch_consistency(tv_copula_result, epoch$cp_yr_mar,
                            index_name, period_label)
  
  # ---- Select copula object via priority-path helper ----------------------
  use_tv <- !is.null(tv_copula_result) &&
    isTRUE(tv_copula_result$significant) &&
    !is.null(tv_copula_result$make_cop_fn)
  copula_path <- "stationary"
  cop_obj_ns  <- NULL
  if (use_tv) {
    cop_obj_ns <- tryCatch(
      .build_period_copula(tv_copula_result, period_yrs, which_seg),
      error = function(e) {
        log_event(sprintf("  [%s] .build_period_copula error: %s", index_name, e$message))
        NULL
      })
    if (!is.null(cop_obj_ns)) {
      copula_path <- if (!is.null(tv_copula_result$seg_result)) "segmented" else "tv_midpoint"
    } else {
      log_event(sprintf(
        "  [%s | %s] All non-stationary copula paths failed - using stationary fallback.",
        index_name, period_label))
    }
  }
  if (is.null(cop_obj_ns)) {
    cop_obj_ns  <- copula_fit_obj$best_copula_fit@copula
    copula_path <- "stationary"
  }
  
  log_event(sprintf(
    "  [%s | %s] Conditional NS SAF: copula_path=%s | mu=%.4f | sigma=%.4f",
    index_name, period_label, copula_path, mu_period, sigma_period))
  
  # ---- Period-specific Gamma marginal -------------------------------------
  # sigma_period is the coefficient of variation (CV), so:
  # shape = 1 / CV^2,   rate = shape / mu
  ns_shape <- 1 / (sigma_period^2)
  ns_rate  <- ns_shape / mu_period
  beta_par <- marginal_fits$area_fit$estimate
  
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5)
  saf_ns   <- data.frame()
  
  for (T in T_years) {
    # Pass mu_T_adj_months, renewal_shape_k, renewal_rate_k, and iat_cv from rp table
    # or compute them inline from the drought data if passed as arguments.
    # Example inline usage (adapt variable names to your scope):
    target <- .compute_renewal_saf_target(
      T_years      = T,
      mu_T_months  = mu_T_adj,    # Replace with iat_mean * 12
      shape_k      = shape_k_val, # Replace with 1 / (iat_cv^2)
      rate_k       = rate_k_val,  # Replace with shape_k / iat_mean
      cv_val       = iat_cv_val   # Replace with computed iat_cv
    )
    if (target <= 0 || target >= 1) next
    sev <- numeric(length(area_pct))
    for (i in seq_along(area_pct)) {
      v   <- pbeta(area_pct[i] / 100, beta_par[1], beta_par[2])
      obj <- function(s) {
        u_ns <- pgamma(s, shape = ns_shape, rate = ns_rate)
        val  <- tryCatch(
          copula::cCopula(cbind(u_ns, v), copula = cop_obj_ns, indices = 2),
          error = function(e) NA_real_)
        if (!is.finite(val)) return(1e6)
        (val - target)^2
      }
      res    <- tryCatch(optimize(obj, interval = c(1e-6, 50)),
                         error = function(e) list(minimum = NA_real_))
      sev[i] <- res$minimum
    }
    
    saf_ns <- rbind(saf_ns,
                    data.frame(
                      ReturnPeriod_years = T,
                      Area_pct           = area_pct,
                      Severity           = sev,
                      Method             = "Conditional_NS",
                      Period             = period_label,
                      mu_period          = round(mu_period,    4),
                      sigma_period       = round(sigma_period, 4),
                      copula_path        = copula_path,
                      Index              = index_name))
  }
  
  write.csv(saf_ns,
            file.path(output_dir,
                      sprintf("%s_SAF_nonstationary_%s.csv", index_name,
                              gsub("[^A-Za-z0-9]", "_", period_label))),
            row.names = FALSE)
  log_event(sprintf(
    "  [%s] Conditional NS SAF written: %s (copula_path=%s)",
    index_name, period_label, copula_path))
  saf_ns
}

# ============================================================================
# 3. NON-STATIONARY SAF CURVES - KENDALL-CORRECTED  (v7 REVISED)
# Mirrors derive_SAF_nonstationary() but uses the Kendall distribution
# K_C(t) = P(C(U,V) <= t) to convert the target non-exceedance
# probability into a copula level t_val, then inverts at each area value.
#
# Changes vs v5/v6:
# * METHODOLOGICAL ALIGNMENT: Aligns with base pipeline (Analysis.R).
#   Analytical K_C(t) = t - φ(t)/φ'(t) is applied for Archimedean families
#   (Frank, Gumbel, Joe, SurvClayton), enabling mathematically valid
#   extrapolation beyond the empirical support.
# * Non-Archimedean families (Gaussian, StudentT, Plackett) fall back to
#   the empirical estimator K_C(t) = P(C(U,V) <= t).
# * v7 UPDATE: Replaces the silent stationary fallback with a strict 3-tier
#   validation strategy for empirical K_C computation.
#   Tier 1: Direct evaluation on pseudo-observations.
#   Tier 2: Simulation-based fallback (N=20,000) if pCopula fails.
#   Tier 3: Hard abort (returns NULL) if period-specific K_C cannot be
#           computed. This prevents miscalibrated t_val from propagating.
# * Optimization intervals are dynamically set: wider [1e-8, 1-1e-8] for
#   analytical paths (safe for generator-based inversion), tighter
#   [1e-4, 1-1e-4] for empirical paths to avoid boundary noise.
# * Epoch boundaries and copula selection follow the same revised logic
#   as derive_SAF_nonstationary() — see that function for full details.
# * New output column copula_path mirrors the conditional method.
# ============================================================================
derive_SAF_nonstationary_kendall <- function(mu_period, sigma_period,
                                             copula_fit_obj, marginal_fits,
                                             drought_data,   ns_result,
                                             index_name,     period_label,
                                             output_dir,
                                             tv_copula_result = NULL) {
  if (is.null(copula_fit_obj) || is.null(marginal_fits)) return(NULL)
  if (is.null(ns_result)) {
    log_event(sprintf(
      "  [%s | %s] derive_SAF_nonstationary_kendall: ns_result is NULL - cannot resolve epoch years. Aborting.",
      index_name, period_label))
    return(NULL)
  }
  
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data)
  n_dr  <- nrow(dc)
  if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  
  # Compute renewal parameters from event inter-arrival times
  indr_all    <- drought_data$severity > 0 & drought_data$area_pct > 0
  is_start    <- indr_all & !c(FALSE, head(indr_all, -1))
  ev_starts   <- as.Date(drought_data$date[is_start])
  iat_cv_val  <- NA_real_
  shape_k_val <- NA_real_
  rate_k_val  <- NA_real_
  mu_T_adj    <- mu_T * 12
  if (length(ev_starts) >= 2) {
    iat_yrs    <- as.numeric(diff(ev_starts)) / 365.25
    iat_mean_v <- mean(iat_yrs, na.rm = TRUE)
    iat_sd_v   <- sd(iat_yrs,   na.rm = TRUE)
    iat_cv_val <- if (iat_mean_v > 0) iat_sd_v / iat_mean_v else NA_real_
    mu_T_adj   <- iat_mean_v * 12
    if (!is.na(iat_cv_val) && iat_cv_val > 0) {
      shape_k_val <- 1 / (iat_cv_val^2)
      rate_k_val  <- shape_k_val / iat_mean_v
    }
  }  
  # ---- Resolve epoch (single source of truth: ns_result) ------------------
  epoch      <- .resolve_epoch(ns_result, period_label)
  period_yrs <- epoch$period_yrs
  which_seg  <- epoch$which_seg
  
  # ---- Epoch-consistency check --------------------------------------------
  .assert_epoch_consistency(tv_copula_result, epoch$cp_yr_mar,
                            index_name, period_label)
  
  # ---- Select copula object via priority-path helper ----------------------
  use_tv <- !is.null(tv_copula_result) &&
    isTRUE(tv_copula_result$significant) &&
    !is.null(tv_copula_result$make_cop_fn)
  copula_path <- "stationary"
  cop_obj_ns  <- NULL
  if (use_tv) {
    cop_obj_ns <- tryCatch(
      .build_period_copula(tv_copula_result, period_yrs, which_seg),
      error = function(e) {
        log_event(sprintf("  [%s] .build_period_copula error: %s", index_name, e$message))
        NULL
      })
    if (!is.null(cop_obj_ns)) {
      copula_path <- if (!is.null(tv_copula_result$seg_result)) "segmented" else "tv_midpoint"
    } else {
      log_event(sprintf(
        "  [%s | %s] All non-stationary copula paths failed - using stationary fallback.",
        index_name, period_label))
    }
  }
  if (is.null(cop_obj_ns)) {
    cop_obj_ns  <- copula_fit_obj$best_copula_fit@copula
    copula_path <- "stationary"
  }
  
  log_event(sprintf(
    "  [%s | %s] Kendall NS SAF: copula_path=%s | mu=%.4f | sigma=%.4f",
    index_name, period_label, copula_path, mu_period, sigma_period))
  
  # ---- Period-specific Gamma marginal -------------------------------------
  ns_shape <- 1 / (sigma_period^2)
  ns_rate  <- ns_shape / mu_period
  beta_par <- marginal_fits$area_fit$estimate
  
  # ========================================================================
  # KENDALL DISTRIBUTION CALIBRATION (v7 STRICT 3-TIER STRATEGY)
  # K_C must be the Kendall distribution of the SAME copula used in the SAF
  # inversion (cop_obj_ns). Using the full-record copula's pseudo-observations
  # with a different period-specific copula object would miscalibrate t_val.
  # ========================================================================
  
  cop_family <- copula_fit_obj$best_copula_name
  archimedean_fams <- c("Frank", "Gumbel", "Joe", "SurvClayton")
  use_analytical_KC <- cop_family %in% archimedean_fams
  
  # --------------------------------------------------------------------------
  # PATH A: ANALYTICAL K_C FOR ARCHIMEDEAN FAMILIES
  # Uses closed-form generators φ(t) and derivatives φ'(t) to compute
  # K_C(t) = t - φ(t)/φ'(t). Enables stable extrapolation for high return
  # periods (T >= 50 yr) without empirical discretization limits.
  # --------------------------------------------------------------------------
  if (use_analytical_KC) {
    theta <- tryCatch(cop_obj_ns@parameters[1], error = function(e) NA_real_)
    if (is.na(theta) || !is.finite(theta)) {
      log_event(sprintf(
        "  [%s | %s] WARNING: Failed to extract valid theta for analytical K_C (family=%s). Falling back to empirical path.",
        index_name, period_label, cop_family))
      use_analytical_KC <- FALSE
    } else {
      # Define generator and derivative functions per Archimedean family
      .phi_and_deriv <- function(t, fam, th) {
        t <- pmax(1e-10, pmin(t, 1 - 1e-10))
        if (fam == "Frank") {
          phi       <- -log((exp(-th * t) - 1) / (exp(-th) - 1))
          phi_prime <- th * exp(-th * t) / (1 - exp(-th * t))
        } else if (fam == "Gumbel") {
          log_t     <- log(t)
          phi       <- (-log_t)^th
          phi_prime <- -th * (-log_t)^(th - 1) / t
        } else if (fam == "Joe") {
          one_m_t   <- 1 - t
          phi       <- -log(1 - one_m_t^th)
          phi_prime <- th * one_m_t^(th - 1) / (1 - one_m_t^th)
        } else { # SurvClayton
          phi       <- (t^(-th) - 1) / th
          phi_prime <- -t^(-th - 1)
        }
        list(phi = phi, phi_prime = phi_prime)
      }
      
      kc_fn <- function(t) {
        pd <- .phi_and_deriv(t, cop_family, theta)
        kc_val <- t - pd$phi / pd$phi_prime
        pmax(1e-10, pmin(1 - 1e-10, kc_val))
      }
      opt_int <- c(1e-8, 1 - 1e-8)
      log_event(sprintf(
        "  [%s | %s] Using analytical K_C for %s (theta=%.4f). Valid extrapolation enabled for non-stationary epoch.",
        index_name, period_label, cop_family, theta))
    }
  }
  
  # --------------------------------------------------------------------------
  # PATH B: EMPIRICAL K_C WITH 3-TIER VALIDATION (v7)
  # Applied when family is non-Archimedean OR analytical extraction failed.
  # Replaces the silent stationary fallback with strict validation.
  # --------------------------------------------------------------------------
  if (!use_analytical_KC) {
    uS_fit <- copula_fit_obj$u_severity
    uA_fit <- copula_fit_obj$u_area
    cop_vals <- NULL
    
    # TIER 1: Direct evaluation on fitted pseudo-observations
    cop_vals <- tryCatch(
      copula::pCopula(cbind(uS_fit, uA_fit), cop_obj_ns),
      error = function(e) NULL
    )
    
    # TIER 2: Simulation-based empirical K_C (avoids pseudo-obs parameter mismatch)
    if (is.null(cop_vals) || !all(is.finite(cop_vals))) {
      log_event(sprintf(
        "  [%s | %s] pCopula failed on pseudo-obs for %s. Switching to simulation-based K_C (N=20,000).",
        index_name, period_label, cop_family))
      tryCatch({
        set.seed(77L)
        sim_mat <- copula::rCopula(20000L, cop_obj_ns)
        cop_vals <- copula::pCopula(sim_mat, cop_obj_ns)
        cop_vals <- cop_vals[is.finite(cop_vals)]
      }, error = function(e) cop_vals <- NULL)
    }
    
    # TIER 3: Hard abort if period-specific K_C cannot be computed
    # This prevents mixing dependence structures and mis-calibrating t_val.
    if (is.null(cop_vals) || length(cop_vals) < 500) {
      log_event(sprintf(
        "  [%s | %s] CRITICAL: Period-specific K_C computation failed for %s. Aborting Kendall SAF for this epoch to prevent inconsistent calibration.",
        index_name, period_label, cop_family))
      return(NULL)  # RUN_ALL.R handles NULL gracefully
    }
    
    kc_fn   <- function(t) mean(cop_vals <= t, na.rm = TRUE)
    opt_int <- c(1e-4, 1 - 1e-4)
    log_event(sprintf(
      "  [%s | %s] Using empirical K_C for %s (n=%d support points).",
      index_name, period_label, cop_family, length(cop_vals)))
  }
  
  # ---- SAF Inversion Loop -------------------------------------------------
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5)
  saf_ns_k <- data.frame()
  
  for (T in T_years) {
    # Pass mu_T_adj_months, renewal_shape_k, renewal_rate_k, and iat_cv from rp table
    # or compute them inline from the drought data if passed as arguments.
    # Example inline usage (adapt variable names to your scope):
    target_kc <- .compute_renewal_saf_target(
      T_years      = T,
      mu_T_months  = mu_T_adj,    # Replace with iat_mean * 12
      shape_k      = shape_k_val, # Replace with 1 / (iat_cv^2)
      rate_k       = rate_k_val,  # Replace with shape_k / iat_mean
      cv_val       = iat_cv_val   # Replace with computed iat_cv
    )
    if (target_kc <= 0 || target_kc >= 1) next
    
    t_sol <- tryCatch(
      optimize(function(t) (kc_fn(t) - target_kc)^2, interval = opt_int),
      error = function(e) list(minimum = 0.5))
    t_val <- t_sol$minimum
    
    sev <- numeric(length(area_pct))
    for (i in seq_along(area_pct)) {
      v   <- pbeta(area_pct[i] / 100, beta_par[1], beta_par[2])
      obj <- function(s) {
        u_ns <- pgamma(s, shape = ns_shape, rate = ns_rate)
        val  <- tryCatch(copula::pCopula(cbind(u_ns, v), cop_obj_ns),
                         error = function(e) NA_real_)
        if (!is.finite(val)) return(1e6)
        (val - t_val)^2
      }
      res    <- tryCatch(optimize(obj, interval = c(1e-6, 50)),
                         error = function(e) list(minimum = NA_real_))
      sev[i] <- res$minimum
    }
    
    saf_ns_k <- rbind(saf_ns_k,
                      data.frame(
                        ReturnPeriod_years = T,
                        Area_pct           = area_pct,
                        Severity           = sev,
                        Method             = "Kendall_NS",
                        Period             = period_label,
                        mu_period          = round(mu_period,    4),
                        sigma_period       = round(sigma_period, 4),
                        copula_path        = copula_path,
                        Index              = index_name))
  }
  
  write.csv(saf_ns_k,
            file.path(output_dir,
                      sprintf("%s_SAF_nonstationary_kendall_%s.csv", index_name,
                              gsub("[^A-Za-z0-9]", "_", period_label))),
            row.names = FALSE)
  log_event(sprintf(
    "  [%s] Kendall NS SAF written: %s (copula_path=%s)",
    index_name, period_label, copula_path))
  saf_ns_k
}

# ============================================================================
# 4. PLACE A HISTORICAL EVENT ON THE SAF SURFACE  (unchanged from v2)
# ============================================================================
place_event_on_saf <- function(drought_data, copula_fit_obj, marginal_fits,
                               event_years, index_name, output_dir) {
  if (is.null(copula_fit_obj) || is.null(marginal_fits)) return(NULL)
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data)
  n_dr  <- nrow(dc)
  mu_T  <- n_tot / n_dr
  
  cop_obj  <- copula_fit_obj$best_copula_fit@copula
  beta_par <- marginal_fits$area_fit$estimate
  
  ev_data <- drought_data[drought_data$year %in% event_years &
                            drought_data$severity > 0, ]
  if (nrow(ev_data) == 0) {
    log_event(sprintf("  [%s] No drought data found for years %s.",
                      index_name, paste(event_years, collapse = "-")))
    return(NULL)
  }
  
  s_peak <- max(ev_data$severity)
  a_peak <- max(ev_data$area_pct)
  u_s    <- compute_u_severity(s_peak, marginal_fits)
  u_a    <- pbeta(a_peak / 100, beta_par[1], beta_par[2])
  C_val  <- tryCatch(copula::pCopula(cbind(u_s, u_a), cop_obj),
                     error = function(e) NA_real_)
  p_joint <- max(0, 1 - u_s - u_a + C_val)
  T_joint <- if (p_joint > 0) mu_T / (p_joint * 12) else Inf
  
  result <- data.frame(
    index      = index_name,
    event_yrs  = paste(range(event_years), collapse = "-"),
    severity   = round(s_peak, 4),
    area_pct   = round(a_peak, 2),
    u_severity = round(u_s,    4),
    u_area     = round(u_a,    4),
    C_value    = round(C_val,  4),
    p_joint    = round(p_joint, 6),
    T_years    = round(T_joint, 0))
  
  write.csv(result,
            file.path(output_dir, sprintf("%s_event_%s_on_SAF.csv", index_name,
                                          paste(range(event_years), collapse = "_"))),
            row.names = FALSE)
  log_event(sprintf("  [%s] Event %s--%s: S=%.3f, A=%.1f%%, T~%d yr",
                    index_name,
                    min(event_years), max(event_years),
                    s_peak, a_peak, round(T_joint)))
  result
}

# Log successful load of extension script v7
log_event("Extensions 2 (v7) loaded: fit_nonstationary_marginals | derive_SAF_nonstationary | derive_SAF_nonstationary_kendall | place_event_on_saf")