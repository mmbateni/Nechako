#============================================================================
# SAF EXTENSIONS 2 — ADDITIONAL SAF METHODS
# Source AFTER Nechako_Drought_SAF_Analysis.R and ext1.R
#
# Provides three new functions:
#   fit_nonstationary_marginals()  — GAMLSS Gamma trend model for severity
#   derive_SAF_nonstationary()     — SAF curves with period-specific marginals
#   place_event_on_saf()           — Returns joint return period of a named event
#
# FIXES vs original ext2:
#   - REMOVED duplicate definitions already in SAF_Analysis.R
#     (derive_SAF_curves_conditional/kendall, create_comparative_SAF_plots,
#      subset_drought_chars_by_class, compute_u_severity)
#   - fit_nonstationary_marginals: reference/recent year windows are now
#     computed from the actual data range instead of being hardcoded to
#     1950-1990 / 2010-2025.  Optional overrides still available.
#============================================================================
if (!requireNamespace("gamlss",      quietly = TRUE)) install.packages("gamlss")
if (!requireNamespace("gamlss.dist", quietly = TRUE)) install.packages("gamlss.dist")
library(gamlss); library(gamlss.dist)

if (!exists("log_event"))
  stop("log_event() not found. Source Nechako_Drought_SAF_Analysis.R first.")

# ---------------------------------------------------------------------------
# 1. NON-STATIONARY MARGINAL FITTING
#    Fits a stationary and a linear-trend GAMLSS Gamma model to drought
#    severity.  Returns period-specific mean severity for two epochs so that
#    derive_SAF_nonstationary() can build "past" vs "recent" SAF curves.
#
#    ref_years / recent_years: optional integer vectors.  When NULL (default)
#    the function divides the record automatically:
#      reference = first 40% of record years
#      recent    = last 30% of record years
# ---------------------------------------------------------------------------
fit_nonstationary_marginals <- function(drought_data, index_name, year_vec,
                                        output_dir,
                                        ref_years    = NULL,
                                        recent_years = NULL) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  if (nrow(dc) < 20L) {
    log_event(sprintf("  [%s] Non-stationary fit skipped: < 20 drought observations.", index_name))
    return(NULL)
  }
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    log_event(sprintf("  [%s] gamlss not available — skipping non-stationary fit.", index_name))
    return(NULL)
  }
  
  yr_match <- year_vec[seq_len(nrow(dc))]
  yr_range <- range(yr_match, na.rm = TRUE)
  yr_span  <- diff(yr_range)
  
  # Data-driven default windows
  if (is.null(ref_years))
    ref_years    <- seq(yr_range[1],
                        yr_range[1] + floor(yr_span * 0.4))
  if (is.null(recent_years))
    recent_years <- seq(yr_range[2] - floor(yr_span * 0.3),
                        yr_range[2])
  
  log_event(sprintf("  [%s] Non-stationary fit | ref: %d-%d | recent: %d-%d",
                    index_name,
                    min(ref_years), max(ref_years),
                    min(recent_years), max(recent_years)))
  
  yr_mean <- mean(yr_match, na.rm = TRUE)
  yr_sd   <- sd(yr_match,   na.rm = TRUE)
  df      <- data.frame(severity  = dc$severity,
                        year_std  = (yr_match - yr_mean) / yr_sd)
  
  m_stat <- tryCatch(
    gamlss::gamlss(severity ~ 1,         family = gamlss.dist::GA(),
                   data = df, trace = FALSE),
    error = function(e) NULL)
  m_ns   <- tryCatch(
    gamlss::gamlss(severity ~ year_std,  family = gamlss.dist::GA(),
                   data = df, trace = FALSE),
    error = function(e) NULL)
  
  if (is.null(m_stat) || is.null(m_ns)) {
    log_event(sprintf("  [%s] GAMLSS fitting failed.", index_name))
    return(NULL)
  }
  
  # AIC comparison
  aic_diff <- AIC(m_ns) - AIC(m_stat)
  log_event(sprintf("  [%s] GAMLSS ΔAIC (non-stat - stat) = %.2f %s",
                    index_name, aic_diff,
                    if (aic_diff < -2) "=> non-stationary preferred" else "=> stationary adequate"))
  
  # FIX (Bug 2, revised): predict.gamlss() with newdata= can return an
  # unexpected non-numeric object ('closure') in some package versions,
  # crashing mean().  Bypass predict entirely: extract the mu coefficients
  # directly and apply the GA() inverse log-link (exp) manually.
  #   For GA() with log link: mu = exp(b0 + b1 * year_std)
  # This is version-proof, transparent, and avoids double-exp() from the
  # earlier incorrect exp(predict(..., type="response")) pattern.
  coef_mu <- coef(m_ns, what = "mu")         # [1] intercept, [2] slope
  to_std  <- function(yrs) (yrs - yr_mean) / yr_sd
  mu_ref  <- mean(exp(coef_mu[1] + coef_mu[2] * to_std(ref_years)),    na.rm = TRUE)
  mu_rec  <- mean(exp(coef_mu[1] + coef_mu[2] * to_std(recent_years)), na.rm = TRUE)
  
  result <- list(m_stat       = m_stat,
                 m_ns         = m_ns,
                 mu_ref       = mu_ref,
                 mu_rec       = mu_rec,
                 sigma        = exp(coef(m_ns, what = "sigma")[1]),
                 ref_years    = ref_years,
                 recent_years = recent_years,
                 aic_diff     = aic_diff)
  
  summary_df <- data.frame(
    index        = index_name,
    mu_ref       = round(mu_ref, 4),
    mu_recent    = round(mu_rec, 4),
    sigma        = round(result$sigma, 4),
    aic_diff     = round(aic_diff, 2),
    ref_start    = min(ref_years),
    ref_end      = max(ref_years),
    recent_start = min(recent_years),
    recent_end   = max(recent_years)
  )
  write.csv(summary_df,
            file.path(output_dir, sprintf("%s_nonstationary_summary.csv", index_name)),
            row.names = FALSE)
  
  result
}

# ---------------------------------------------------------------------------
# 2. NON-STATIONARY SAF CURVES
#    Derives SAF curves with a period-specific Gamma marginal for severity
#    (mu_period, sigma_ns) while keeping the fitted copula and Beta area
#    marginal unchanged.  Produces one CSV per period label.
# ---------------------------------------------------------------------------
derive_SAF_nonstationary <- function(mu_period, sigma_ns,
                                     copula_fit_obj, marginal_fits,
                                     drought_data,   index_name,
                                     period_label,   output_dir) {
  if (is.null(copula_fit_obj) || is.null(marginal_fits)) return(NULL)
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data); n_dr <- nrow(dc); if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  
  cop_obj  <- copula_fit_obj$best_copula_fit@copula
  beta_par <- marginal_fits$area_fit$estimate
  
  # Period-specific Gamma shape & rate (method of moments)
  ns_shape <- 1 / (sigma_ns^2)
  ns_rate  <- ns_shape / mu_period
  
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5)
  saf_ns   <- data.frame()
  
  for (T in T_years) {
    target <- 1 - mu_T / (T * 12)
    if (target <= 0 || target >= 1) next
    sev <- numeric(length(area_pct))
    for (i in seq_along(area_pct)) {
      v <- pbeta(area_pct[i] / 100, beta_par[1], beta_par[2])
      obj <- function(s) {
        u_ns <- pgamma(s, shape = ns_shape, rate = ns_rate)
        (copula::cCopula(cbind(u_ns, v), copula = cop_obj, indices = 2) - target)^2
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
                               mu_period          = round(mu_period, 4),
                               Index              = index_name))
  }
  
  write.csv(saf_ns,
            file.path(output_dir, sprintf("%s_SAF_nonstationary_%s.csv", index_name,
                                          gsub("[^A-Za-z0-9]", "_", period_label))),
            row.names = FALSE)
  log_event(sprintf("  [%s] Non-stationary SAF written: %s", index_name, period_label))
  saf_ns
}

# ---------------------------------------------------------------------------
# 3. NON-STATIONARY SAF CURVES — KENDALL-CORRECTED  (Issue 5 fix)
#    Mirrors derive_SAF_nonstationary() but applies the Kendall distribution
#    correction (K_C) before the Gaussian transform, consistent with the
#    paper's Equation 11.  The copula structure is assumed time-invariant; only
#    the severity marginal shifts between periods.
#
#    Steps (per return period T):
#      a) Compute empirical K_C from the fitted copula pseudo-obs (stationary).
#      b) Solve K_C(t_val) = 1 - mu_T / (T*12)  for t_val.
#      c) For each area value a, find severity s such that C(u_NS(s), v(a)) = t_val
#         where u_NS uses the period-specific Gamma CDF.
# ---------------------------------------------------------------------------
derive_SAF_nonstationary_kendall <- function(mu_period, sigma_ns,
                                             copula_fit_obj, marginal_fits,
                                             drought_data,   index_name,
                                             period_label,   output_dir) {
  if (is.null(copula_fit_obj) || is.null(marginal_fits)) return(NULL)
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data); n_dr <- nrow(dc); if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  
  cop_obj  <- copula_fit_obj$best_copula_fit@copula
  beta_par <- marginal_fits$area_fit$estimate
  
  # Period-specific Gamma parameters (method of moments: sigma_ns = CV)
  ns_shape <- 1 / (sigma_ns^2)
  ns_rate  <- ns_shape / mu_period
  
  # Empirical Kendall distribution estimated from the stationary pseudo-obs.
  # The copula structure is kept fixed across periods; only the severity
  # marginal changes.
  uS_fit   <- copula_fit_obj$u_severity
  uA_fit   <- copula_fit_obj$u_area
  cop_vals <- copula::pCopula(cbind(uS_fit, uA_fit), cop_obj)
  kc_fn    <- function(t) mean(cop_vals <= t)   # scalar empirical K_C(t)
  
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5)
  saf_ns_k <- data.frame()
  
  for (T in T_years) {
    target_kc <- 1 - mu_T / (T * 12)
    if (target_kc <= 0 || target_kc >= 1) next
    
    # Solve K_C(t_val) = target_kc
    kc_obj <- function(t) (kc_fn(t) - target_kc)^2
    t_sol  <- tryCatch(optimize(kc_obj, interval = c(1e-4, 1 - 1e-4)),
                       error = function(e) list(minimum = 0.5))
    t_val  <- t_sol$minimum
    
    sev <- numeric(length(area_pct))
    for (i in seq_along(area_pct)) {
      v   <- pbeta(area_pct[i] / 100, beta_par[1], beta_par[2])
      obj <- function(s) {
        u_ns <- pgamma(s, shape = ns_shape, rate = ns_rate)
        val  <- tryCatch(copula::pCopula(cbind(u_ns, v), cop_obj), error = function(e) NA_real_)
        if (!is.finite(val)) return(1e6)   # guard: pCopula can return NA/NaN for some copulas
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
                                 mu_period          = round(mu_period, 4),
                                 Index              = index_name))
  }
  
  write.csv(saf_ns_k,
            file.path(output_dir, sprintf("%s_SAF_nonstationary_kendall_%s.csv", index_name,
                                          gsub("[^A-Za-z0-9]", "_", period_label))),
            row.names = FALSE)
  log_event(sprintf("  [%s] Non-stationary Kendall SAF written: %s", index_name, period_label))
  saf_ns_k
}

# ---------------------------------------------------------------------------
# 4. PLACE A HISTORICAL EVENT ON THE SAF SURFACE
#    Finds the joint return period of a specific drought event identified by
#    year range, using the inclusion-exclusion identity:
#      P(U > u AND V > v) = 1 - u - v + C(u, v)
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
  
  # Joint exceedance probability: P(S > s_peak AND A > a_peak)
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
    T_years    = round(T_joint, 0)
  )
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

log_event("Extensions 2 loaded: fit_nonstationary_marginals | derive_SAF_nonstationary | derive_SAF_nonstationary_kendall | place_event_on_saf")