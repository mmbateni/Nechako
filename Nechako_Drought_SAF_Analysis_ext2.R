#============================================================================
# SAF EXTENSIONS 2 — ADDITIONAL SAF METHODS (v4)
# Source AFTER Nechako_Drought_SAF_Analysis.R and ext1.R
#
# Functions exported:
#   fit_nonstationary_marginals()      — Three-model GAMLSS comparison:
#                                        stationary, mu-trend, mu+sigma-trend.
#                                        Epoch definition from change-point or
#                                        automatic data split. [UPDATED]
#   derive_SAF_nonstationary()         — SAF curves with period-specific marginals
#                                        AND optional period-specific copula
#                                        when TV dependence is significant. [UPDATED]
#   derive_SAF_nonstationary_kendall() — Kendall-corrected version of above. [UPDATED]
#   place_event_on_saf()               — Joint return period of a named event.
#
# CHANGES vs v3 (BUG FIX):
#
#   fit_nonstationary_marginals() — epoch boundary truncation [FIX]:
#
#     v3 BUG:  yr_match <- year_vec[seq_len(nrow(dc))]
#              year_vec is the 912-element monthly year vector (1950–2025).
#              nrow(dc) is the count of DROUGHT MONTHS (e.g. ~550 for SPI-1).
#              Subsetting gives years for the FIRST ~550 months of the record
#              (≈ 1950–1996), not the actual years of the drought observations.
#              Consequently yr_range[2] is artificially truncated, the
#              "recent" epoch ends decades before the record end (e.g. 1998
#              instead of 2025), and change-point fractions are inflated so
#              late change points are misclassified as "near record edge".
#
#     v4 FIX:  yr_match <- dc$year
#              Uses the actual calendar year stored in each drought-month row.
#              yr_range now correctly spans 1950–2025 and recent_years extends
#              to the record end.  The GAMLSS trend fit also receives correct
#              year_std values for every observation.
#
#              Impact on results:
#                SPI-1 : recent epoch 1965–2025 (was 1965–1998; +27 yr)
#                SPEI-1: recent epoch 1996–2025 (was auto-split due to
#                         artificial near-edge flag; now uses CP correctly)
#                SPI-3 : recent epoch 1973–2025 (was 1973–1998; +27 yr)
#                SPEI-3: recent epoch 1987–2025 (was 1987–1999; +26 yr)
#============================================================================

for (pkg in c("gamlss", "gamlss.dist")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
library(gamlss); library(gamlss.dist)

if (!exists("log_event"))
  stop("log_event() not found. Source Nechako_Drought_SAF_Analysis.R first.")

# ---------------------------------------------------------------------------
# INTERNAL HELPER: build a period-specific copula from a TV copula result
#
#   tv_res       : return value of fit_timevarying_copula()
#   period_years : integer vector of calendar years in this epoch
#
# Returns a copula object with θ evaluated at the epoch midpoint year.
# ---------------------------------------------------------------------------
.build_period_copula <- function(tv_res, period_years) {
  mid_yr   <- mean(range(period_years, na.rm = TRUE))
  yr_std_m <- (mid_yr - tv_res$yr_mean) / tv_res$yr_sd
  theta_r  <- tv_res$a_hat + tv_res$b_hat * yr_std_m
  # make_cop_fn is the closure returned by fit_timevarying_copula — it
  # already incorporates family selection, bounds clamping, and rotCopula
  # wrapping for SurvClayton, so no further switch needed here.
  tryCatch(tv_res$make_cop_fn(theta_r), error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# 1. NON-STATIONARY MARGINAL FITTING  (updated v3)
#
#    Three GAMLSS Gamma models are fitted and compared by AIC:
#      M0 (stationary):    mu ~ 1,         sigma ~ 1
#      M1 (mu-trend):      mu ~ year_std,  sigma ~ 1
#      M2 (mu+sigma-trend):mu ~ year_std,  sigma ~ year_std
#
#    The best model is selected by minimum AIC.  Period-specific mu and
#    sigma values are derived for two epochs (reference and recent).
#
#    Epoch definition (precedence order):
#      1. Explicit ref_years / recent_years arguments — highest priority.
#      2. cp_result$detected == TRUE — cp_year splits:
#           reference = record start … cp_year - 1
#           recent    = cp_year … record end
#         Note: if the change point is very early or very late (< 20% or
#         > 80% of record), the 40%/30% split is used instead with a warning.
#      3. Default automatic split: first 40% / last 30% of record.
#
#    Parameters
#    ----------
#    drought_data  : output of extract_drought_characteristics
#    index_name    : character label
#    year_vec      : full year vector (all months, not just drought obs)
#    output_dir    : directory for CSV output
#    cp_result     : optional; output of detect_copula_changepoints() [NEW]
#    ref_years     : optional integer vector override for reference epoch
#    recent_years  : optional integer vector override for recent epoch
# ---------------------------------------------------------------------------
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
    log_event(sprintf("  [%s] gamlss not available — skipping.", index_name))
    return(NULL)
  }
  
  # FIX (v4): use the calendar year stored in each drought-observation row.
  # The previous yr_match <- year_vec[seq_len(nrow(dc))] took the first
  # nrow(dc) elements of the 912-element monthly year vector, yielding years
  # for the FIRST nrow(dc) months of the record (e.g. 1950–1996 for ~550
  # drought months) rather than the actual observation years.  This truncated
  # yr_range[2] and caused recent_years to end decades before 2025, and
  # inflated change-point fractions so late CPs were misclassified as edge cases.
  yr_match <- dc$year
  yr_range <- range(yr_match, na.rm = TRUE)
  yr_span  <- diff(yr_range)
  
  # ---- Epoch definition --------------------------------------------------
  auto_ref    <- seq(yr_range[1], yr_range[1] + floor(yr_span * 0.4))
  auto_recent <- seq(yr_range[2] - floor(yr_span * 0.3), yr_range[2])
  
  if (!is.null(ref_years) && !is.null(recent_years)) {
    epoch_src <- "user override"
    # ref_years and recent_years already set by caller
    
  } else if (!is.null(cp_result) && isTRUE(cp_result$detected) && !is.na(cp_result$cp_year)) {
    cp_yr   <- cp_result$cp_year
    cp_frac <- (cp_yr - yr_range[1]) / yr_span   # 0-1 position in record
    if (cp_frac < 0.20 || cp_frac > 0.80) {
      log_event(sprintf(
        "  [%s] Change-point at %d is near record edge (%.0f%% of full record) — using auto 40%%/30%% split instead.",
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
  
  yr_mean  <- mean(yr_match, na.rm = TRUE)
  yr_sd    <- sd(yr_match,   na.rm = TRUE)
  if (yr_sd < 1e-10) yr_sd <- 1
  to_std   <- function(yrs) (yrs - yr_mean) / yr_sd
  
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
  
  # M2: both location and scale trend
  m2 <- tryCatch(
    gamlss::gamlss(severity ~ year_std,
                   sigma.formula = ~ year_std,
                   family = gamlss.dist::GA(), data = df, trace = FALSE),
    error = function(e) NULL)
  
  if (is.null(m0)) {
    log_event(sprintf("  [%s] Stationary GAMLSS fit failed — aborting.", index_name))
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
  #  FIX (Bug 2): predict.gamlss() may return a closure in some versions.
  #  Coefficients are extracted directly and the GA() log-link applied
  #  manually: mu(t) = exp(b0 + b1 * year_std(t)).
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

# ---------------------------------------------------------------------------
# 2. NON-STATIONARY SAF CURVES — CONDITIONAL METHOD  (updated v3)
#
#    Derives SAF curves with period-specific Gamma marginal (mu_period,
#    sigma_period) for drought severity.
#
#    Copula selection (new in v3):
#      If tv_copula_result is provided AND lr_p < 0.05 (significant trend
#      in dependence), the copula object is evaluated at the midpoint of
#      the epoch rather than using the stationary fit.  This propagates
#      non-stationarity in BOTH marginals and dependence into the SAF curve.
#      If the TV copula is not significant, the stationary copula is used
#      (original behaviour).
#
#    Parameters
#    ----------
#    mu_period, sigma_period : period-specific Gamma parameters (mean, CV)
#    copula_fit_obj          : stationary copula result from fit_copulas()
#    marginal_fits           : output of fit_marginal_distributions()
#    drought_data            : characteristics data frame
#    index_name, period_label, output_dir : metadata / paths
#    tv_copula_result        : optional; output of fit_timevarying_copula() [NEW]
# ---------------------------------------------------------------------------
derive_SAF_nonstationary <- function(mu_period, sigma_period,
                                     copula_fit_obj, marginal_fits,
                                     drought_data,   index_name,
                                     period_label,   output_dir,
                                     tv_copula_result = NULL) {
  if (is.null(copula_fit_obj) || is.null(marginal_fits)) return(NULL)
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data)
  n_dr  <- nrow(dc)
  if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  
  # -- Select copula object -----------------------------------------------
  period_yrs <- if (grepl("ref", tolower(period_label), fixed = FALSE) &&
                    !is.null(tv_copula_result$cp_result$cp_year)) {
    seq(min(drought_data$year, na.rm = TRUE),
        tv_copula_result$cp_result$cp_year - 1L)
  } else {
    seq(if (!is.null(tv_copula_result$cp_result$cp_year))
      tv_copula_result$cp_result$cp_year
      else min(drought_data$year, na.rm = TRUE),
      max(drought_data$year, na.rm = TRUE))
  }
  
  use_tv <- !is.null(tv_copula_result) &&
    isTRUE(tv_copula_result$significant) &&
    !is.null(tv_copula_result$make_cop_fn)
  if (use_tv) {
    cop_obj_ns <- tryCatch(.build_period_copula(tv_copula_result, period_yrs),
                           error = function(e) NULL)
    if (!is.null(cop_obj_ns)) {
      log_event(sprintf("  [%s] Conditional NS SAF: using TV copula @ epoch midpoint for %s",
                        index_name, period_label))
    } else {
      log_event(sprintf("  [%s] TV copula build failed — falling back to stationary.", index_name))
      cop_obj_ns <- copula_fit_obj$best_copula_fit@copula
    }
  } else {
    cop_obj_ns <- copula_fit_obj$best_copula_fit@copula
  }
  
  beta_par <- marginal_fits$area_fit$estimate
  
  # Period-specific Gamma: sigma_period is the coefficient of variation (CV)
  # so  shape = 1/CV^2,  rate = shape/mu
  ns_shape <- 1 / (sigma_period^2)
  ns_rate  <- ns_shape / mu_period
  
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5)
  saf_ns   <- data.frame()
  
  for (T in T_years) {
    target <- 1 - mu_T / (T * 12)
    if (target <= 0 || target >= 1) next
    sev <- numeric(length(area_pct))
    for (i in seq_along(area_pct)) {
      v   <- pbeta(area_pct[i] / 100, beta_par[1], beta_par[2])
      obj <- function(s) {
        u_ns <- pgamma(s, shape = ns_shape, rate = ns_rate)
        (copula::cCopula(cbind(u_ns, v), copula = cop_obj_ns, indices = 2) - target)^2
      }
      res    <- tryCatch(optimize(obj, interval = c(0.001, 50)),
                         error = function(e) list(minimum = 0))
      sev[i] <- res$minimum
    }
    saf_ns <- rbind(saf_ns,
                    data.frame(ReturnPeriod_years = T,
                               Area_pct           = area_pct,
                               Severity           = sev,
                               Method             = "Conditional_NS",
                               Period             = period_label,
                               mu_period          = round(mu_period,    4),
                               sigma_period       = round(sigma_period, 4),
                               tv_copula_used     = use_tv && !is.null(cop_obj_ns),
                               Index              = index_name))
  }
  
  write.csv(saf_ns,
            file.path(output_dir, sprintf("%s_SAF_nonstationary_%s.csv", index_name,
                                          gsub("[^A-Za-z0-9]", "_", period_label))),
            row.names = FALSE)
  log_event(sprintf("  [%s] Conditional NS SAF written: %s (TV copula used: %s)",
                    index_name, period_label, use_tv))
  saf_ns
}

# ---------------------------------------------------------------------------
# 3. NON-STATIONARY SAF CURVES — KENDALL-CORRECTED  (updated v3)
#
#    Mirrors derive_SAF_nonstationary() but uses the Kendall distribution
#    K_C(t) = P(C(U,V) ≤ t) to convert the target non-exceedance probability
#    into a copula level t_val, then inverts the copula at each area value.
#
#    The Kendall distribution is estimated from the pseudo-observations of
#    the APPROPRIATE copula (stationary or period-specific TV), consistent
#    with De Michele et al. (2013) and extensions.
#
#    tv_copula_result: same semantics as in derive_SAF_nonstationary().
# ---------------------------------------------------------------------------
derive_SAF_nonstationary_kendall <- function(mu_period, sigma_period,
                                             copula_fit_obj, marginal_fits,
                                             drought_data,   index_name,
                                             period_label,   output_dir,
                                             tv_copula_result = NULL) {
  if (is.null(copula_fit_obj) || is.null(marginal_fits)) return(NULL)
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data)
  n_dr  <- nrow(dc)
  if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  
  # -- Select copula object -----------------------------------------------
  period_yrs <- if (grepl("ref", tolower(period_label), fixed = FALSE) &&
                    !is.null(tv_copula_result$cp_result$cp_year)) {
    seq(min(drought_data$year, na.rm = TRUE),
        tv_copula_result$cp_result$cp_year - 1L)
  } else {
    seq(if (!is.null(tv_copula_result$cp_result$cp_year))
      tv_copula_result$cp_result$cp_year
      else min(drought_data$year, na.rm = TRUE),
      max(drought_data$year, na.rm = TRUE))
  }
  
  use_tv <- !is.null(tv_copula_result) &&
    isTRUE(tv_copula_result$significant) &&
    !is.null(tv_copula_result$make_cop_fn)
  if (use_tv) {
    cop_obj_ns <- tryCatch(.build_period_copula(tv_copula_result, period_yrs),
                           error = function(e) NULL)
    if (!is.null(cop_obj_ns)) {
      log_event(sprintf("  [%s] Kendall NS SAF: using TV copula @ epoch midpoint for %s",
                        index_name, period_label))
    } else {
      cop_obj_ns <- copula_fit_obj$best_copula_fit@copula
    }
  } else {
    cop_obj_ns <- copula_fit_obj$best_copula_fit@copula
  }
  
  beta_par <- marginal_fits$area_fit$estimate
  
  ns_shape <- 1 / (sigma_period^2)
  ns_rate  <- ns_shape / mu_period
  
  # Empirical Kendall distribution from the period-specific (or stationary)
  # copula evaluated at the observed pseudo-observations
  uS_fit   <- copula_fit_obj$u_severity
  uA_fit   <- copula_fit_obj$u_area
  cop_vals <- tryCatch(copula::pCopula(cbind(uS_fit, uA_fit), cop_obj_ns),
                       error = function(e) copula::pCopula(cbind(uS_fit, uA_fit),
                                                           copula_fit_obj$best_copula_fit@copula))
  kc_fn    <- function(t) mean(cop_vals <= t)
  
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5)
  saf_ns_k <- data.frame()
  
  for (T in T_years) {
    target_kc <- 1 - mu_T / (T * 12)
    if (target_kc <= 0 || target_kc >= 1) next
    
    t_sol <- tryCatch(optimize(function(t) (kc_fn(t) - target_kc)^2,
                               interval = c(1e-4, 1 - 1e-4)),
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
                      data.frame(ReturnPeriod_years = T,
                                 Area_pct           = area_pct,
                                 Severity           = sev,
                                 Method             = "Kendall_NS",
                                 Period             = period_label,
                                 mu_period          = round(mu_period,    4),
                                 sigma_period       = round(sigma_period, 4),
                                 tv_copula_used     = use_tv && !is.null(cop_obj_ns),
                                 Index              = index_name))
  }
  
  write.csv(saf_ns_k,
            file.path(output_dir, sprintf("%s_SAF_nonstationary_kendall_%s.csv", index_name,
                                          gsub("[^A-Za-z0-9]", "_", period_label))),
            row.names = FALSE)
  log_event(sprintf("  [%s] Kendall NS SAF written: %s (TV copula used: %s)",
                    index_name, period_label, use_tv))
  saf_ns_k
}

# ---------------------------------------------------------------------------
# 4. PLACE A HISTORICAL EVENT ON THE SAF SURFACE  (unchanged from v2)
# ---------------------------------------------------------------------------
place_event_on_saf <- function(drought_data, copula_fit_obj, marginal_fits,
                               event_years, index_name, output_dir) {
  if (is.null(copula_fit_obj) || is.null(marginal_fits)) return(NULL)
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data); n_dr <- nrow(dc); mu_T <- n_tot / n_dr
  
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
  log_event(sprintf("  [%s] Event %s--%s: S=%.3f, A=%.1f%%, T≈%d yr",
                    index_name,
                    min(event_years), max(event_years),
                    s_peak, a_peak, round(T_joint)))
  result
}

log_event("Extensions 2 (v4) loaded: fit_nonstationary_marginals | derive_SAF_nonstationary | derive_SAF_nonstationary_kendall | place_event_on_saf")