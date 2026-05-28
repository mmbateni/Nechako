# =============================================================================
# SAF EXTENSIONS 1 — ENHANCED FUNCTIONS
# Source AFTER Nechako_Drought_SAF_Analysis.R
# =============================================================================
#
# Functions exported:
#   test_dependency()             — Kendall τ / Spearman ρ between severity & area
#   gof_copula()                  — Bootstrap Sn GOF test (N_boot = 499 replicates)
#   detect_copula_changepoints()  — Pettitt / CUSUM / PRUTF structural-break
#   rosenblatt_pit_gof()          — Rosenblatt PIT GOF for time-varying copulas
#   fit_timevarying_copula()      — Trend in dependence; calls change-point
#                                   detection first and Rosenblatt PIT afterwards.
#                                   Returns yr_mean / yr_sd / a_hat / b_hat so
#                                   ext2 can build period-specific copula objects.
#   fit_segmented_copula()        — Separate copulas for two temporal segments.
#
# DIAGNOSTIC PLOTTING FUNCTIONS:
#   plot_copula_diagnostic()      — Empirical vs theoretical copula contour plot
#   plot_rosenblatt_diagnostic()  — Rosenblatt PIT histogram / QQ / scatter
#   plot_rolling_tau()            — Rolling Kendall τ time series
#   plot_tv_parameter()           — TV copula parameter evolution
#
# KEY CHANGES IN THIS VERSION:
#   - detect_copula_changepoints(): zero-variance rolling window guard.
#     Pre-checks sd(severity) and sd(area_pct) for each window.
#     Zero-variance windows are assigned NA silently; the total count of
#     suppressed windows is logged once via log_event().
#   - gof_copula(): N_boot = 499 replicates documented.
#     499 replicates yield a minimum two-sided p-value of 0.002 (< α = 0.05)
#     while keeping run-time manageable.  The CUSUM permutation test also uses
#     499 MC draws (matching N_boot for consistency).
# =============================================================================

for (pkg in c("moments", "openxlsx", "Kendall", "trend")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(moments)
library(openxlsx)
library(Kendall)
library(trend)

if (!exists("log_event")) {
  stop("log_event() not found. Source Nechako_Drought_SAF_Analysis.R first.")
}


# =============================================================================
# 1. DEPENDENCY TEST
# =============================================================================

test_dependency <- function(drought_data, index_name) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  if (nrow(dc) < 5) {
    log_event(sprintf("  [%s] Dependency test skipped: too few drought observations.", index_name))
    return(list(kendall = NA, spearman = NA, p_value = NA))
  }
  
  S     <- dc$severity
  A     <- dc$area_pct
  cor_k <- cor(S, A, method = "kendall")
  cor_s <- cor(S, A, method = "spearman")
  t_res <- cor.test(S, A, method = "kendall")
  
  log_event(sprintf(
    "  [%s] Dependency \u2014 Kendall tau=%.4f (p=%.4e) | Spearman rho=%.4f",
    index_name, cor_k, t_res$p.value, cor_s
  ))
  
  list(kendall  = cor_k,
       spearman = cor_s,
       p_value  = t_res$p.value)
}


# =============================================================================
# 2. COPULA GOODNESS-OF-FIT  (Bootstrap Sn test)
# =============================================================================
#
# Bootstrap GOF via the parametric Cramér–von Mises Sn statistic.
# N_boot (default 499) is the number of bootstrap replicates used to
# approximate the null distribution of Sn.  499 replicates strike a practical
# balance: the smallest achievable two-sided p-value is 1/(499+1) = 0.002,
# which is well below the α = 0.05 decision threshold, while keeping run-time
# manageable (each replicate re-fits the copula to a synthetic sample of the
# same size).  Set N_boot = 999 for publication-quality p-values; use
# N_boot = 99 during development for speed.
# =============================================================================

gof_copula <- function(copula_fit_obj, u_matrix, index_name, N_boot = 499L) {
  
  if (is.null(copula_fit_obj)) return(NULL)
  
  cop <- copula_fit_obj$best_copula_fit@copula
  cop <- .fix_tcopula_df(cop)
  
  res <- tryCatch(
    suppressWarnings(
      copula::gofCopula(cop, u_matrix, N = N_boot,
                        method = "Sn", optim.method = "BFGS",
                        ties = TRUE)
    ),
    error = function(e) {
      log_event(sprintf("  [%s] GOF test error: %s", index_name, e$message))
      NULL
    }
  )
  if (is.null(res)) return(NULL)
  
  decision <- if (res$p.value > 0.05) "ACCEPTED" else "REJECTED"
  log_event(sprintf(
    "  [%s] Copula GOF (%s): Sn=%.4f, p=%.4f => %s",
    index_name, copula_fit_obj$best_copula_name,
    res$statistic, res$p.value, decision
  ))
  
  list(
    statistic   = res$statistic,
    p_value     = res$p.value,
    decision    = decision,
    copula_name = copula_fit_obj$best_copula_name
  )
}


# =============================================================================
# 3. CHANGE-POINT DETECTION IN DEPENDENCE STRUCTURE
#
# Three complementary tests applied to a rolling-window Kendall τ series
# computed from the drought pseudo-observations:
#
# (a) Pettitt test  (trend::pettitt.test)
#     Rank-based non-parametric test for a single shift in the location
#     of a continuous sequence.  Well-suited to small samples.
#
# (b) CUSUM  (manual implementation with permutation p-value)
#     Cumulative sum of standardised τ deviations.  Detects gradual or abrupt
#     shifts; maximum |CUSUM| is the test statistic.
#     499 Monte-Carlo permutations are used (matching N_boot in gof_copula).
#
# (c) PRUTF — Parametric LR structural-break test
#     For each candidate break point k ∈ {20%, 35%, 50%, 65%, 80%} of n,
#     the full-record log-likelihood is compared with the sum of two segment
#     log-likelihoods (each segment refitted independently).
#     Test statistic = sup_k { 2(ll_seg1_k + ll_seg2_k - ll_full) }.
#     p-value from chi-sq(1) — conservative since we take a supremum,
#     which is appropriate as a first-order approximation.
#
# A change point is flagged as "detected" when at least 2 of the 3 tests
# reject H0 (no change) at the 5% level.  The consensus change-point year is
# the median of the significant test estimates.
#
# Parameters
# ----------
# drought_data    : output data frame from extract_drought_characteristics
# copula_fit_obj  : result from fit_copulas() — used for PRUTF
# index_name      : character label for logging
# output_dir      : directory for CSV output
# window_yrs      : full width of the rolling τ window (default 10)
# =============================================================================

detect_copula_changepoints <- function(drought_data, copula_fit_obj,
                                       index_name,   output_dir,
                                       window_yrs = 10L) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n  <- nrow(dc)
  
  if (n < 30L) {
    log_event(sprintf(
      "  [%s] Change-point detection skipped: < 30 drought obs.", index_name
    ))
    return(list(detected  = FALSE,  cp_year   = NA_integer_,
                pettitt_p = NA,     cusum_p   = NA,
                prutf_p   = NA))
  }
  
  # ---- Rolling Kendall τ ---------------------------------------------------
  # window_yrs (default 10) sets the FULL width of the centred window.
  # Each centre point i uses observations in [i - half_w, i + half_w].
  # A minimum of 6 observations per window is required for a meaningful τ.
  #
  # ZERO-VARIANCE GUARD:
  #   When a window contains drought events with zero variance in either
  #   severity or area_pct, cor() would emit a warning and return NA.
  #   Pre-checking variance silences those warnings; a single aggregate
  #   diagnostic is logged so the analyst can inspect the data if needed.
  half_w             <- max(floor(window_yrs / 2L), 3L)
  tau_ser            <- rep(NA_real_, n)
  n_zero_var_windows <- 0L
  
  for (i in seq_len(n)) {
    idx    <- max(1L, i - half_w):min(n, i + half_w)
    if (length(idx) >= 6L) {
      sev_w  <- dc$severity[idx]
      area_w <- dc$area_pct[idx]
      if (sd(sev_w, na.rm = TRUE) < 1e-10 || sd(area_w, na.rm = TRUE) < 1e-10) {
        n_zero_var_windows <- n_zero_var_windows + 1L
      } else {
        tau_ser[i] <- cor(sev_w, area_w, method = "kendall")
      }
    }
  }
  
  if (n_zero_var_windows > 0L) {
    log_event(sprintf(
      "  [%s] Rolling tau: %d window(s) of %d total had zero variance in severity or area_pct and were set to NA (no warning emitted). Inspect raw drought characteristics for runs of identical values.",
      index_name, n_zero_var_windows, n
    ))
  }
  
  valid <- which(is.finite(tau_ser))
  if (length(valid) < 10L) {
    log_event(sprintf(
      "  [%s] Change-point detection: too few valid tau values.", index_name
    ))
    return(list(detected  = FALSE,  cp_year   = NA_integer_,
                pettitt_p = NA,     cusum_p   = NA,
                prutf_p   = NA))
  }
  
  tau_v <- tau_ser[valid]
  yr_v  <- dc$year[valid]
  
  # ---- (a) Pettitt test ---------------------------------------------------
  pett      <- tryCatch(trend::pettitt.test(tau_v), error = function(e) NULL)
  pett_pval <- if (!is.null(pett)) pett$p.value             else NA_real_
  pett_cp   <- if (!is.null(pett)) yr_v[pett$estimate[[1]]] else NA_integer_
  
  # ---- (b) CUSUM with permutation p-value ---------------------------------
  s_tau      <- sd(tau_v, na.rm = TRUE)
  if (s_tau < 1e-10) s_tau <- 1
  z          <- (tau_v - mean(tau_v)) / s_tau
  cs         <- cumsum(z)
  cusum_stat <- max(abs(cs))
  cusum_cp_idx  <- which.max(abs(cs))
  cusum_cp_yr   <- yr_v[cusum_cp_idx]
  
  set.seed(123L)
  mc_cusum   <- replicate(499L, max(abs(cumsum(sample(z, replace = FALSE)))))
  cusum_pval <- mean(mc_cusum >= cusum_stat)
  
  # ---- (c) PRUTF: parametric LR structural-break --------------------------
  # Helper to build a one-parameter copula of the correct family
  cop_nm  <- copula_fit_obj$best_copula_name
  init_p  <- coef(copula_fit_obj$best_copula_fit)
  uS      <- copula_fit_obj$u_severity
  uA      <- copula_fit_obj$u_area
  uDat    <- cbind(uS, uA)
  
  make_cop_prutf <- function(th) {
    # th is a full coefficient vector (length 1 for single-param families,
    # length 2 for StudentT).  Always index th[1] for scalar families so that
    # passing coef(fit) — which may be a named numeric — works safely.
    if (cop_nm == "StudentT") {
      rho <- max(min(th[1], 0.99), -0.99)
      df  <- max(if (length(th) >= 2L) th[2] else 4, 1.5)
      return(tCopula(param = rho, df = df, dim = 2))
    }
    switch(cop_nm,
           Frank       = frankCopula(th[1],               dim = 2),
           Gumbel      = gumbelCopula(max(th[1], 1.001),  dim = 2),
           Plackett    = plackettCopula(max(th[1], 1e-3)),
           SurvClayton = rotCopula(claytonCopula(max(th[1], 1e-3), dim = 2)),
           Joe         = joeCopula(max(th[1], 1.001),     dim = 2),
           Gaussian    = normalCopula(max(min(th[1], 0.99), -0.99), dim = 2, dispstr = "un"),
           frankCopula(th[1], dim = 2)   # safe fallback
    )
  }
  
  fit_seg_ll <- function(rows) {
    if (length(rows) < 5L) return(NA_real_)
    uM  <- uDat[rows, , drop = FALSE]
    fit <- tryCatch(
      fitCopula(make_cop_prutf(init_p), uM, method = "mpl",
                optim.method = "BFGS", optim.control = list(maxit = 300)),
      error = function(e)
        tryCatch(
          fitCopula(make_cop_prutf(init_p), uM, method = "mpl",
                    optim.method = "Nelder-Mead",
                    optim.control = list(maxit = 600)),
          error = function(e2) NULL
        )
    )
    if (is.null(fit)) return(NA_real_)
    tryCatch(
      loglikCopula(coef(fit), uM, make_cop_prutf(coef(fit))),
      error = function(e) NA_real_
    )
  }
  
  ll_full   <- tryCatch(
    loglikCopula(init_p, uDat, .fix_tcopula_df(copula_fit_obj$best_copula_fit@copula)),
    error = function(e) NA_real_
  )
  k_cands   <- unique(round(quantile(seq_len(n), c(0.20, 0.35, 0.50, 0.65, 0.80))))
  lr_vals   <- numeric(length(k_cands))
  
  for (j in seq_along(k_cands)) {
    k <- k_cands[j]
    ll1       <- fit_seg_ll(seq_len(k))
    ll2       <- fit_seg_ll((k + 1L):n)
    lr_vals[j] <- if (is.finite(ll1) && is.finite(ll2) && is.finite(ll_full)) {
      max(0, 2 * (ll1 + ll2 - ll_full))
    } else {
      0
    }
  }
  
  prutf_stat   <- max(lr_vals)
  prutf_cp_idx <- k_cands[which.max(lr_vals)]
  prutf_cp_yr  <- dc$year[min(prutf_cp_idx, n)]
  
  # Conservative chi-sq(1) p-value (supremum inflates Type I error → treat as
  # upper bound; a formal Andrews 1993 critical value would be slightly lower).
  prutf_pval <- pchisq(prutf_stat, df = 1L, lower.tail = FALSE)
  
  # ---- Consensus ----------------------------------------------------------
  cp_pool  <- c(pett_cp, cusum_cp_yr, prutf_cp_yr)
  pv_pool  <- c(pett_pval, cusum_pval, prutf_pval)
  sig_idx  <- which(is.finite(pv_pool) & pv_pool < 0.05)
  detected <- length(sig_idx) >= 2L
  cp_year  <- if (detected) {
    as.integer(round(median(cp_pool[sig_idx], na.rm = TRUE)))
  } else {
    NA_integer_
  }
  
  log_event(sprintf(
    "  [%s] ChangePoint | Pettitt: yr=%s p=%.3f | CUSUM: yr=%s p=%.3f | PRUTF: yr=%s p=%.3f | Consensus yr=%s (detected=%s)",
    index_name,
    if (is.na(pett_cp)) "NA" else pett_cp,   round(pett_pval,  3),
    cusum_cp_yr,                               round(cusum_pval, 3),
    prutf_cp_yr,                               round(prutf_pval, 3),
    if (is.na(cp_year)) "NA" else cp_year,     detected
  ))
  
  result_df <- data.frame(
    index           = index_name,
    pettitt_cp_yr   = pett_cp,      pettitt_p  = round(pett_pval,  4),
    cusum_cp_yr     = cusum_cp_yr,  cusum_p    = round(cusum_pval, 4),
    prutf_cp_yr     = prutf_cp_yr,  prutf_p    = round(prutf_pval, 4),
    consensus_cp_yr = cp_year,      detected   = detected,
    stringsAsFactors = FALSE
  )
  write.csv(
    result_df,
    file.path(output_dir, sprintf("%s_changepoint_detection.csv", index_name)),
    row.names = FALSE
  )
  
  invisible(list(
    detected   = detected,
    cp_year    = cp_year,
    pettitt_cp = pett_cp,     pettitt_p = pett_pval,
    cusum_cp   = cusum_cp_yr, cusum_p   = cusum_pval,
    prutf_cp   = prutf_cp_yr, prutf_p   = prutf_pval,
    tau_series = data.frame(year = yr_v, rolling_tau = tau_v)
  ))
}


# =============================================================================
# INTERNAL HELPER: .h_func_fd — h-function via analytical / finite differences
#
# Computes  C(u2 | u1; cop_obj) = ∂C(u1, u2)/∂u1.
# This is the correct Rosenblatt second variate: it conditions u_area ON
# u_severity, not the reverse.
#
# Analytical formulae are used for Archimedean and elliptical families.
# For other families, a central finite-difference fallback is applied.
#
# Using pCopula (the joint CDF) rather than cCopula avoids:
#   (a) the wrong-direction bug  (cCopula indices=1 gives C(u1|u2))
#   (b) the Plackett NA problem  (cCopula has no Plackett method)
#
# eps: half-width of the finite-difference step (default 1e-4).
#      Smaller values improve accuracy but may amplify floating-point
#      noise; 1e-4 is safe for the [0,1]² domain.
# =============================================================================

.h_func_fd <- function(u1, u2, cop_obj, cop_family) {
  
  theta <- cop_obj@parameters[1]
  u1    <- pmax(1e-10, pmin(u1, 1 - 1e-10))
  u2    <- pmax(1e-10, pmin(u2, 1 - 1e-10))
  
  if (cop_family == "Gumbel") {
    log_u1  <- -log(u1); log_u2 <- -log(u2)
    sum_pow <- log_u1^theta + log_u2^theta
    h_val   <- exp(-sum_pow^(1 / theta)) * sum_pow^(1 / theta - 1) *
      log_u1^(theta - 1) / u1
    
  } else if (cop_family == "Frank") {
    exp_m_theta    <- exp(-theta)
    exp_m_theta_u1 <- exp(-theta * u1)
    exp_m_theta_u2 <- exp(-theta * u2)
    num   <- (exp_m_theta_u1 - 1) * (exp_m_theta_u2 - 1)
    h_val <- (exp_m_theta_u1 * (exp_m_theta_u2 - 1)) /
      ((exp_m_theta - 1) * (1 + num / (exp_m_theta - 1)))
    
  } else if (cop_family == "Joe") {
    one_m_u1 <- 1 - u1; one_m_u2 <- 1 - u2
    f     <- one_m_u1^theta + one_m_u2^theta - one_m_u1^theta * one_m_u2^theta
    h_val <- f^(1 / theta - 1) * one_m_u1^(theta - 1) * (1 - one_m_u2^theta)
    
  } else if (cop_family == "StudentT") {
    rho   <- cop_obj@parameters[1]
    nu    <- if (length(cop_obj@parameters) >= 2L) cop_obj@parameters[2] else 4L
    t1    <- qt(u1, df = nu); t2 <- qt(u2, df = nu)
    h_val <- pt(
      (t2 - rho * t1) / sqrt(((nu + t1^2) * (1 - rho^2)) / (nu + 1)),
      df = nu + 1
    )
    
  } else if (cop_family == "Gaussian") {
    rho   <- pmax(pmin(cop_obj@parameters[1], 0.999999), -0.999999)
    t_u1  <- qnorm(u1); t_u2 <- qnorm(u2)
    h_val <- pnorm((t_u2 - rho * t_u1) / sqrt(1 - rho^2))
    
  } else {
    h_val <- tryCatch({
      copula::cCopula(cbind(u1, u2), copula = cop_obj, indices = 2)[, 1L]
    }, error = function(e) {
      eps   <- 1e-5
      u1_lo <- max(eps,     u1 - eps)
      u1_hi <- min(1 - eps, u1 + eps)
      C_lo  <- copula::pCopula(c(u1_lo, u2), cop_obj)
      C_hi  <- copula::pCopula(c(u1_hi, u2), cop_obj)
      (C_hi - C_lo) / (u1_hi - u1_lo)
    })
  }
  
  pmax(1e-8, pmin(1 - 1e-8, h_val))
}


# =============================================================================
# 4. ROSENBLATT PROBABILITY INTEGRAL TRANSFORM GOF
#
# The Rosenblatt (1952) transform converts a bivariate sample (u, v) from a
# correctly specified copula C into two independent U[0,1] variates:
#
#   e1_i = u_S_i                                     (first margin — trivially uniform)
#   e2_i = C_{A|S}(u_A_i | u_S_i; θ_i)             (second Rosenblatt variate)
#        = ∂C(u_S_i, u_A_i; θ_i) / ∂u_S_i
#
# For a TIME-VARYING copula each observation i has its own θ_i, so the
# partial derivative is evaluated at the observation-specific parameter.
# e2 is computed via .h_func_fd().
#
# Adequacy criteria (all must hold):
#   • e2 is uniformly distributed: KS p > 0.05
#   • e1 and e2 are independent:   Kendall τ p > 0.05
# (e1 is always uniform by marginal construction; its KS result is
#  logged for completeness but not included in the adequacy decision.)
# =============================================================================

rosenblatt_pit_gof <- function(uDat, make_cop_fn, a_hat, b_hat,
                               year_std, index_name, cop_name) {
  n  <- nrow(uDat)
  e1 <- uDat[, 1]   # u_severity — trivially U[0,1] by marginal construction
  
  e2 <- vapply(seq_len(n), function(i) {
    tryCatch({
      cop_i <- make_cop_fn(a_hat + b_hat * year_std[i])
      .h_func_fd(uDat[i, 1L], uDat[i, 2L], cop_i, cop_family = cop_name)
    }, error = function(e) NA_real_)
  }, numeric(1L))
  
  e2_ok <- e2[is.finite(e2)]
  e1_ok <- e1[is.finite(e2)]
  frac  <- length(e2_ok) / n
  
  log_event(sprintf(
    "[%s] e2 summary: min=%.4f, max=%.4f, mean=%.4f, sd=%.4f, n_valid=%d",
    index_name, min(e2_ok), max(e2_ok), mean(e2_ok), sd(e2_ok), length(e2_ok)
  ))
  
  ks_e1 <- tryCatch(ks.test(e1_ok, "punif"),
                    error = function(e) list(statistic = NA_real_, p.value = NA_real_))
  ks_e2 <- tryCatch(ks.test(e2_ok, "punif"),
                    error = function(e) list(statistic = NA_real_, p.value = NA_real_))
  
  tau_ind <- tryCatch(cor(e1_ok, e2_ok, method = "kendall"),
                      error = function(e) NA_real_)
  tau_p   <- tryCatch(cor.test(e1_ok, e2_ok, method = "kendall")$p.value,
                      error = function(e) NA_real_)
  
  adequate <- isTRUE(ks_e2$p.value > 0.05) && isTRUE(tau_p > 0.05)
  
  log_event(sprintf(
    "  [%s] Rosenblatt PIT (%s) | e1 KS p=%.4f | e2 KS p=%.4f | indep tau=%.4f (p=%.4f) | adequate=%s (%.0f%% e2 valid)",
    index_name, cop_name,
    ks_e1$p.value, ks_e2$p.value, tau_ind, tau_p,
    adequate, 100 * frac
  ))
  
  list(
    ks_e1_stat = unname(ks_e1$statistic), ks_e1_p   = ks_e1$p.value,
    ks_e2_stat = unname(ks_e2$statistic), ks_e2_p   = ks_e2$p.value,
    tau_ind    = tau_ind,                 tau_ind_p  = tau_p,
    adequate   = adequate
  )
}


# =============================================================================
# 5. TIME-VARYING COPULA
#
# Pipeline:
#   (i)   detect_copula_changepoints() — run before the trend fit.
#          If a change point is detected the result is embedded in the return
#          value so that ext2 can use it for epoch definition.
#   (ii)  Fit stationary and linear-trend (log-linear on link scale) copula
#          models; compare with LRT (chi-sq, df = 1).
#   (iii) rosenblatt_pit_gof() — verifies the fitted TV copula via the
#          Rosenblatt transform; result is appended to the return list.
#   (iv)  Segmented copula fit (when a change point is detected).
#   (v)   Model selection: stationary | tv | segmented.
#
# Added to return value (for ext2 derive_SAF_nonstationary):
#   yr_mean, yr_sd   — to re-standardise a calendar year outside this fn
#   a_hat, b_hat     — TV copula parameters on the link scale
#   cop_name         — family name (for make_cop reconstruction)
#   link_hi, link_lo — parameter bounds
#   link_ilink       — inverse-link function
#   cp_result        — output of detect_copula_changepoints()
#   rosenblatt       — output of rosenblatt_pit_gof()
# =============================================================================

fit_timevarying_copula <- function(drought_data, copula_fit_obj, marginal_fits,
                                   year_vec,     index_name,     output_dir) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  if (nrow(dc) < 20L) {
    log_event(sprintf(
      "  [%s] TV copula skipped: < 20 drought observations.", index_name
    ))
    return(NULL)
  }
  
  yr_match <- dc$year
  if (length(yr_match) != nrow(dc)) yr_match <- rep(mean(year_vec), nrow(dc))
  yr_mean  <- mean(yr_match, na.rm = TRUE)
  yr_sd    <- sd(yr_match,   na.rm = TRUE)
  if (yr_sd < 1e-10) yr_sd <- 1
  year_std <- (yr_match - yr_mean) / yr_sd
  
  uS       <- copula_fit_obj$u_severity
  uA       <- copula_fit_obj$u_area
  uDat     <- cbind(uS, uA)
  n        <- nrow(uDat)
  cop_name <- copula_fit_obj$best_copula_name
  init_par <- coef(copula_fit_obj$best_copula_fit)[1]
  
  # (i) Change-point detection -----------------------------------------------
  cp_result <- tryCatch(
    detect_copula_changepoints(drought_data, copula_fit_obj, index_name, output_dir),
    error = function(e) list(detected = FALSE, cp_year = NA_integer_)
  )
  
  # ---- Family-specific link / inverse-link ---------------------------------
  link_info <- switch(cop_name,
                      Frank   = list(link  = identity, ilink = identity,
                                     lo = -50, hi = 50, init = init_par),
                      Gumbel  = list(link  = function(x) log(pmax(x - 1, 1e-6)),
                                     ilink = function(x) exp(x) + 1,
                                     lo = 1.001, hi = 50, init = max(init_par, 1.001)),
                      Plackett = list(link  = function(x) log(pmax(x, 1e-6)),
                                      ilink = exp,
                                      lo = 1e-3, hi = 500, init = max(init_par, 1e-3)),
                      SurvClayton = list(link  = function(x) log(pmax(x, 1e-6)),
                                         ilink = exp,
                                         lo = 1e-3, hi = 50, init = max(init_par, 1e-3)),
                      Joe      = list(link  = function(x) log(pmax(x - 1, 1e-6)),
                                      ilink = function(x) exp(x) + 1,
                                      lo = 1.001, hi = 50, init = max(init_par, 1.001)),
                      Gaussian = list(link  = function(x) atanh(pmax(-0.9999, pmin(0.9999, x))),
                                      ilink = tanh,
                                      lo = -0.99, hi = 0.99, init = init_par),
                      StudentT = list(link  = function(x) c(atanh(x[1]), log(x[2])),
                                      ilink = function(x) c(tanh(x[1]), exp(x[2])),
                                      lo = c(-0.99, 1.5), hi = c(0.99, 50),
                                      init = c(init_par, 4)),
                      # default fallback
                      list(link  = function(x) log(pmax(x, 1e-6)), ilink = exp,
                           lo = 1e-3, hi = 50, init = max(init_par, 1e-3))
  )
  
  make_cop <- function(theta_raw) {
    th <- link_info$ilink(theta_raw)
    if (cop_name == "StudentT") {
      th[1] <- max(min(th[1], link_info$hi[1]), link_info$lo[1])
      th[2] <- max(min(th[2], link_info$hi[2]), link_info$lo[2])
    } else {
      th <- max(min(th, link_info$hi), link_info$lo)
    }
    switch(cop_name,
           Frank       = frankCopula(th,  dim = 2),
           Gumbel      = gumbelCopula(th, dim = 2),
           Plackett    = plackettCopula(th),
           SurvClayton = rotCopula(claytonCopula(th, dim = 2)),
           Joe         = joeCopula(th,    dim = 2),
           Gaussian    = normalCopula(th, dim = 2),
           StudentT    = tCopula(param = th[1], df = th[2], dim = 2),
           frankCopula(th, dim = 2)   # default fallback
    )
  }
  
  # (ii) Fit stationary & TV -------------------------------------------------
  n_params <- length(link_info$init)
  if (n_params > 1) {
    log_event(sprintf(
      "  [%s] TV trend skipped for %s: multi-parameter families require extended formulation. Using stationary fit.",
      index_name, cop_name
    ))
    return(invisible(list(
      index = index_name, copula = cop_name, selected_model = "stationary",
      a_hat = NA, b_hat = 0, lr_stat = NA, lr_p = NA, significant = FALSE,
      cp_result = list(detected = FALSE, cp_year = NA), seg_result = NULL,
      aic_stat = NA, bic_stat = NA, aic_tv = NA, bic_tv = NA, aic_seg = NA, bic_seg = NA,
      ros_stat = NULL, ros_tv = NULL, ros_seg = NULL,
      yr_mean = yr_mean, yr_sd = yr_sd, link_ilink = link_info$ilink,
      link_lo = link_info$lo, link_hi = link_info$hi, make_cop_fn = make_cop
    )))
  }
  
  nll_stat <- function(p1) -sum(log(pmax(copula::dCopula(uDat, make_cop(p1)), 1e-300)))
  nll_tv   <- function(p) {
    th_raw <- p[1] + p[2] * year_std
    ll <- 0
    for (i in seq_len(nrow(uDat))) {
      ll <- ll + log(pmax(copula::dCopula(uDat[i, , drop = FALSE], make_cop(th_raw[i])), 1e-300))
    }
    -ll
  }
  
  opt_stat <- tryCatch(
    optim(link_info$link(link_info$init), nll_stat, method = "Brent", lower = -10, upper = 10),
    error = function(e) NULL
  )
  opt_tv <- tryCatch(
    optim(c(opt_stat$par, 0), nll_tv, method = "BFGS"),
    error = function(e) NULL
  )
  
  ll_stat <- -opt_stat$value
  ll_tv   <- -opt_tv$value
  lr      <- 2 * (ll_tv - ll_stat)
  p_val   <- pchisq(max(0, lr), df = 1L, lower.tail = FALSE)
  
  # (iii) Segmented fit (when change point is detected) ----------------------
  seg_result <- NULL
  if (isTRUE(cp_result$detected) && !is.na(cp_result$cp_year)) {
    seg_result <- tryCatch(
      fit_segmented_copula(drought_data, copula_fit_obj,
                           cp_result$cp_year, index_name, output_dir),
      error = function(e) NULL
    )
  }
  ll_seg <- if (!is.null(seg_result)) {
    seg_result$seg1$best_ll + seg_result$seg2$best_ll
  } else {
    -Inf
  }
  
  # (iv) AIC / BIC (k: Stat = 1, TV = 2, Seg = 2) ---------------------------
  aic_stat <- -2 * ll_stat + 2 * 1;  bic_stat <- -2 * ll_stat + log(n) * 1
  aic_tv   <- -2 * ll_tv   + 2 * 2;  bic_tv   <- -2 * ll_tv   + log(n) * 2
  aic_seg  <- if (is.finite(ll_seg)) -2 * ll_seg + 2 * 2 else Inf
  bic_seg  <- if (is.finite(ll_seg)) -2 * ll_seg + log(n) * 2 else Inf
  
  # (v) Rosenblatt PIT GOF ---------------------------------------------------
  ros_stat <- tryCatch(
    rosenblatt_pit_gof(
      uDat, make_cop,
      if (!is.null(opt_stat)) opt_stat$par[1] else link_info$link(link_info$init),
      0, year_std, index_name, cop_name
    ),
    error = function(e) NULL
  )
  ros_tv <- tryCatch(
    rosenblatt_pit_gof(uDat, make_cop, opt_tv$par[1], opt_tv$par[2],
                       year_std, index_name, cop_name),
    error = function(e) NULL
  )
  
  ros_seg <- NULL
  if (!is.null(seg_result)) {
    cp_yr    <- cp_result$cp_year
    seg1_idx <- which(dc$year <  cp_yr)
    seg2_idx <- which(dc$year >= cp_yr)
    make_seg_cop <- function(i) {
      if (i %in% seg1_idx) {
        make_cop(link_info$link(coef(seg_result$seg1$best_copula_fit)[1]))
      } else {
        make_cop(link_info$link(coef(seg_result$seg2$best_copula_fit)[1]))
      }
    }
    e2_seg <- vapply(seq_len(n), function(i) {
      tryCatch(
        .h_func_fd(uDat[i, 1], uDat[i, 2], make_seg_cop(i), cop_family = cop_name),
        error = function(e) NA_real_
      )
    }, numeric(1L))
    e2_ok  <- e2_seg[is.finite(e2_seg)]
    e1_ok  <- uDat[is.finite(e2_seg), 1]
    ks_p   <- tryCatch(ks.test(e2_ok, "punif")$p.value,          error = function(e) 0)
    tau_p  <- tryCatch(cor.test(e1_ok, e2_ok, method = "kendall")$p.value, error = function(e) 0)
    ros_seg <- list(ks_e2_p = ks_p, tau_ind_p = tau_p,
                    adequate = (ks_p > 0.05 & tau_p > 0.05))
  }
  
  # (vi) Model selection -----------------------------------------------------
  tv_sig        <- isTRUE(p_val < 0.05)
  cp_strong     <- isTRUE(cp_result$detected)
  tv_better_aic <- isTRUE(aic_tv  < aic_seg)
  pit_tv_ok     <- isTRUE(!is.null(ros_tv) && ros_tv$adequate)
  pit_tv_score  <- if (!is.null(ros_tv)) {
    mean(c(ros_tv$ks_e2_p, ros_tv$tau_ind_p), na.rm = TRUE)
  } else {
    0
  }
  if (is.nan(pit_tv_score)) pit_tv_score <- 0
  
  seg_better_aic   <- isTRUE(aic_seg < aic_tv)
  pit_seg_improved <- isTRUE(
    if (!is.null(ros_seg) && !is.null(ros_tv)) {
      mean(c(ros_seg$ks_e2_p, ros_seg$tau_ind_p), na.rm = TRUE) > pit_tv_score
    } else {
      FALSE
    }
  )
  
  selected <- "stationary"
  if (isTRUE(tv_sig) && !isTRUE(cp_strong) && isTRUE(tv_better_aic) && isTRUE(pit_tv_ok)) {
    selected <- "tv"
  } else if (isTRUE(cp_strong) && isTRUE(seg_better_aic) && isTRUE(pit_seg_improved)) {
    selected <- "segmented"
  }
  
  log_event(sprintf(
    "  [%s] Model Selection: %s (Stat AIC=%.1f | TV AIC=%.1f | Seg AIC=%s | TV sig=%s | CP=%s)",
    index_name, toupper(selected),
    aic_stat, aic_tv,
    if (is.finite(aic_seg)) sprintf("%.1f", aic_seg) else "NA",
    tv_sig, cp_strong
  ))
  
  # (vii) Return object ------------------------------------------------------
  invisible(list(
    index          = index_name,
    copula         = cop_name,
    selected_model = selected,
    a_hat          = if (selected == "tv") {
      opt_tv$par[1]
    } else if (selected == "stationary") {
      if (!is.null(opt_stat)) opt_stat$par[1] else opt_tv$par[1]
    } else {
      NA
    },
    b_hat        = if (selected == "tv") opt_tv$par[2] else 0,
    lr_stat      = round(lr,    4),
    lr_p         = round(p_val, 4),
    significant  = tv_sig,
    cp_detected  = cp_strong,
    cp_year      = cp_result$cp_year,
    aic_stat = aic_stat, bic_stat = bic_stat,
    aic_tv   = aic_tv,   bic_tv   = bic_tv,
    aic_seg  = aic_seg,  bic_seg  = bic_seg,
    ros_stat = ros_stat, ros_tv   = ros_tv,  ros_seg = ros_seg,
    yr_mean      = yr_mean,
    yr_sd        = yr_sd,
    link_ilink   = link_info$ilink,
    link_lo      = link_info$lo,
    link_hi      = link_info$hi,
    make_cop_fn  = make_cop,
    cp_result    = cp_result,
    seg_result   = seg_result,
    rosenblatt   = ros_tv
  ))
}


# =============================================================================
# 6. SEGMENTED COPULA
#
# Fits separate copulas to two temporal segments defined by a change-point year.
# Appropriate when all three of the following hold:
#   (a) detect_copula_changepoints() reports detected = TRUE
#   (b) fit_timevarying_copula() LR p < 0.05 (significant temporal trend)
#   (c) rosenblatt_pit_gof() reports adequate = FALSE for the TV copula
# When (a)+(b)+(c) are simultaneously true the abrupt-break model is
# statistically preferred over a smooth linear trend, and the Rosenblatt test
# shows that even the TV copula does not adequately represent the dependence.
#
# The candidate copula set is the same seven families used in fit_copulas()
# (Frank, Gumbel, Plackett, SurvClayton, Joe, Gaussian, StudentT), selected
# per-segment by AIC.
#
# Parameters
# ----------
# drought_data   : output of extract_drought_characteristics
# copula_fit_obj : result from fit_copulas() — supplies u_severity, u_area
# cp_year        : integer; the consensus change-point year
# index_name     : character label for logging
# output_dir     : directory for CSV output
#
# Returns a list with:
#   seg1, seg2       : per-segment copula result lists
#   cp_year          : the split year used
#   lr_seg_vs_full_p : p-value of LR test vs full stationary model
# =============================================================================

fit_segmented_copula <- function(drought_data, copula_fit_obj,
                                 cp_year,      index_name,
                                 output_dir) {
  if (is.null(copula_fit_obj) || is.null(cp_year) || is.na(cp_year)) {
    log_event(sprintf(
      "  [%s] Segmented copula: no valid change-point \u2014 aborting.", index_name
    ))
    return(NULL)
  }
  
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n  <- nrow(dc)
  
  idx1 <- which(dc$year <  cp_year)
  idx2 <- which(dc$year >= cp_year)
  n1   <- length(idx1); n2 <- length(idx2)
  
  log_event(sprintf(
    "  [%s] Segmented copula: cp_year=%d | seg1 n=%d [%d\u2013%d] | seg2 n=%d [%d\u2013%d]",
    index_name, cp_year,
    n1, min(dc$year[idx1], na.rm = TRUE), max(dc$year[idx1], na.rm = TRUE),
    n2, min(dc$year[idx2], na.rm = TRUE), max(dc$year[idx2], na.rm = TRUE)
  ))
  
  if (n1 < 15L || n2 < 15L) {
    log_event(sprintf(
      "  [%s] Segmented copula: one or both segments have < 15 observations (n1=%d, n2=%d) \u2014 aborting.",
      index_name, n1, n2
    ))
    return(NULL)
  }
  
  # Candidate set (same seven families as fit_copulas)
  cops <- list(
    Frank       = frankCopula(dim = 2),
    Gumbel      = gumbelCopula(dim = 2),
    Plackett    = plackettCopula(),
    SurvClayton = rotCopula(claytonCopula(dim = 2)),
    Joe         = joeCopula(dim = 2),
    Gaussian    = normalCopula(param = 0.5, dim = 2, dispstr = "un"),
    StudentT    = tCopula(param = 0.5, df = 4, dim = 2)
  )
  
  # Internal helper: fit all candidates to a sub-matrix of pseudo-observations
  .fit_seg <- function(uM_sub, seg_label) {
    res   <- data.frame(Copula = character(), Parameter = numeric(),
                        LogLik = numeric(), AIC = numeric(), BIC = numeric(),
                        stringsAsFactors = FALSE)
    fits  <- list()
    n_sub <- nrow(uM_sub)
    
    for (nm in names(cops)) {
      fit <- tryCatch(
        fitCopula(cops[[nm]], uM_sub, method = "mpl",
                  optim.method  = "BFGS",
                  optim.control = list(maxit = 1000, reltol = 1e-8)),
        error = function(e)
          tryCatch(
            fitCopula(cops[[nm]], uM_sub, method = "mpl",
                      optim.method  = "Nelder-Mead",
                      optim.control = list(maxit = 2000)),
            error = function(e2) NULL
          )
      )
      if (!is.null(fit)) {
        p  <- coef(fit)
        ll <- loglikCopula(p, uM_sub, cops[[nm]])
        k  <- length(p)
        res <- rbind(res, data.frame(
          Copula    = nm,
          Parameter = p[1],
          LogLik    = ll,
          AIC       = -2 * ll + 2 * k,
          BIC       = -2 * ll + log(n_sub) * k
        ))
        fits[[nm]] <- fit
      }
    }
    
    if (nrow(res) == 0) return(NULL)
    best_nm  <- res$Copula[which.min(res$AIC)]
    best_aic <- res$AIC[which.min(res$AIC)]
    best_ll  <- res$LogLik[which.min(res$AIC)]
    log_event(sprintf(
      "  [%s]   Segment %s best copula: %s (AIC=%.2f)",
      index_name, seg_label, best_nm, best_aic
    ))
    list(
      best_copula_name = best_nm,
      best_copula_fit  = fits[[best_nm]],
      best_ll          = best_ll,
      all_results      = res,
      n                = n_sub
    )
  }
  
  uS_all <- copula_fit_obj$u_severity
  uA_all <- copula_fit_obj$u_area
  uM_all <- cbind(uS_all, uA_all)
  
  seg1 <- .fit_seg(uM_all[idx1, , drop = FALSE], sprintf("1 [<%d]",  cp_year))
  seg2 <- .fit_seg(uM_all[idx2, , drop = FALSE], sprintf("2 [\u2265%d]", cp_year))
  
  if (is.null(seg1) || is.null(seg2)) {
    log_event(sprintf(
      "  [%s] Segmented copula: fitting failed for one or both segments.", index_name
    ))
    return(NULL)
  }
  
  # ---- Model comparison: segmented vs full-record stationary ---------------
  #
  # NOTE ON VALIDITY:
  #   If seg1 and seg2 select the same copula family as the full model, the
  #   LR statistic 2*(ll_seg - ll_full) is asymptotically chi-sq(1) under H0
  #   (one extra parameter: 2 segment params vs 1 full param).
  #
  #   If the families DIFFER across segments the models are NON-NESTED.
  #   The chi-sq(1) reference distribution is then INVALID and can be
  #   anti-conservative.  For non-nested comparison, delta-AIC is the
  #   appropriate criterion:
  #     delta_AIC < -2  => segmented model preferred (threshold: Burnham & Anderson 2002)
  #
  #   PRIMARY DECISION: delta_AIC
  #   SUPPLEMENTARY:    LR p-value (valid only when families match)
  ll_full          <- copula_fit_obj$all_results$LogLik[
    copula_fit_obj$all_results$Copula == copula_fit_obj$best_copula_name
  ]
  ll_seg   <- seg1$best_ll + seg2$best_ll
  lr_stat  <- max(0, 2 * (ll_seg - ll_full))
  lr_p     <- pchisq(lr_stat, df = 1L, lower.tail = FALSE)
  
  n_full        <- nrow(uM_all)
  aic_full_here <- -2 * ll_full + 2 * 1
  aic_seg_here  <- -2 * ll_seg  + 2 * 2
  delta_aic     <- aic_seg_here - aic_full_here   # negative = segmented preferred
  
  families_match    <- (seg1$best_copula_name == copula_fit_obj$best_copula_name &&
                          seg2$best_copula_name == copula_fit_obj$best_copula_name)
  lr_valid          <- families_match
  seg_preferred_aic <- delta_aic < -2
  
  log_event(sprintf(
    "  [%s] Segmented vs full | LR: stat=%.3f p=%.4f (valid=%s, families: full=%s seg1=%s seg2=%s) | delta_AIC=%.2f => %s",
    index_name, lr_stat, lr_p, lr_valid,
    copula_fit_obj$best_copula_name, seg1$best_copula_name, seg2$best_copula_name,
    delta_aic,
    if (seg_preferred_aic) "segmented PREFERRED (AIC)" else "no meaningful AIC improvement"
  ))
  
  # ---- Upper-tail dependence comparison between segments -------------------
  # (Lower-tail dependence is not assessed — irrelevant for drought SAF.)
  td_seg <- function(uM, label, fam_name) {
    if (fam_name %in% c("Gumbel", "Joe", "SurvClayton")) {
      A_05 <- compute_cfg_pickands(uM[, 1], uM[, 2], t_val = 0.5)
      lU   <- max(0, 2 * (1 - A_05))
      meth <- "CFG"
    } else {
      u_q_hi <- 0.95
      C_hi   <- max(mean(uM[, 1] >= u_q_hi & uM[, 2] >= u_q_hi), 1e-8)
      lU     <- max(0, min(1, 2 - log(C_hi) / log(u_q_hi)))
      meth   <- "Quantile_0.95"
    }
    log_event(sprintf(
      "  [%s]   Tail dep segment %s: \u03BB_U=%.4f [%s]",
      index_name, label, lU, meth
    ))
    c(lambda_U = lU)
  }
  
  td1 <- td_seg(uM_all[idx1, , drop = FALSE], sprintf("1 [<%d]",     cp_year), seg1$best_copula_name)
  td2 <- td_seg(uM_all[idx2, , drop = FALSE], sprintf("2 [\u2265%d]", cp_year), seg2$best_copula_name)
  
  # ---- Write summary CSV ---------------------------------------------------
  sum_df <- data.frame(
    index         = index_name,
    cp_year       = cp_year,
    seg1_years    = sprintf("%d\u2013%d", min(dc$year[idx1]), max(dc$year[idx1])),
    seg2_years    = sprintf("%d\u2013%d", min(dc$year[idx2]), max(dc$year[idx2])),
    seg1_n        = n1,
    seg2_n        = n2,
    seg1_best_cop = seg1$best_copula_name,
    seg2_best_cop = seg2$best_copula_name,
    seg1_par      = round(coef(seg1$best_copula_fit)[1], 4),
    seg2_par      = round(coef(seg2$best_copula_fit)[1], 4),
    seg1_lambda_U = round(td1["lambda_U"], 4),
    seg2_lambda_U = round(td2["lambda_U"], 4),
    ll_full       = round(ll_full,  3),
    ll_segmented  = round(ll_seg,   3),
    lr_stat       = round(lr_stat,  3),
    lr_p          = round(lr_p,     4),
    lr_valid      = lr_valid,
    aic_full      = round(aic_full_here, 2),
    aic_seg       = round(aic_seg_here,  2),
    delta_aic     = round(delta_aic,     2),
    seg_preferred = seg_preferred_aic,
    stringsAsFactors = FALSE
  )
  write.csv(
    sum_df,
    file.path(output_dir, sprintf("%s_segmented_copula.csv", index_name)),
    row.names = FALSE
  )
  
  invisible(list(
    seg1              = seg1,
    seg2              = seg2,
    cp_year           = cp_year,
    seg1_years        = dc$year[idx1],
    seg2_years        = dc$year[idx2],
    seg1_lambda_U     = td1["lambda_U"],
    seg2_lambda_U     = td2["lambda_U"],
    lr_stat           = lr_stat,
    lr_seg_vs_full_p  = lr_p,
    seg_preferred     = lr_p < 0.05
  ))
}


# =============================================================================
# DIAGNOSTIC PLOTTING FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Empirical vs theoretical copula contour plot
# -----------------------------------------------------------------------------

plot_copula_diagnostic <- function(uDat, cop_obj, index_name, output_dir,
                                   n_grid = 100) {
  if (is.null(uDat) || nrow(uDat) < 10) return(NULL)
  
  grid    <- seq(0.01, 0.99, length.out = n_grid)
  Z_theo  <- matrix(NA, n_grid, n_grid)
  Z_emp   <- matrix(NA, n_grid, n_grid)
  
  for (i in seq_len(n_grid)) {
    for (j in seq_len(n_grid)) {
      Z_theo[i, j] <- tryCatch(
        copula::pCopula(c(grid[i], grid[j]), cop_obj),
        error = function(e) NA
      )
      Z_emp[i, j] <- mean(uDat[, 1] <= grid[i] & uDat[, 2] <= grid[j])
    }
  }
  
  diff_mat <- Z_emp - Z_theo
  df_tile  <- data.frame(
    u1   = rep(grid, n_grid),
    u2   = rep(grid, each = n_grid),
    diff = as.vector(diff_mat)
  )
  df_pts <- data.frame(u1 = uDat[, 1], u2 = uDat[, 2])
  
  p <- ggplot(df_tile, aes(u1, u2, fill = diff)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Emp-Theo", limits = c(-0.15, 0.15)) +
    geom_point(data    = df_pts, aes(u1, u2),
               inherit.aes = FALSE, alpha = 0.1, size = 0.3, colour = "black") +
    labs(title = paste(index_name, "Copula Fit Diagnostic"),
         x = "U(severity)", y = "U(area)") +
    theme_bw(base_size = 11)
  
  ggsave(file.path(output_dir, sprintf("%s_copula_diagnostic.pdf", index_name)),
         p, width = 7, height = 6)
  invisible(p)
}


# -----------------------------------------------------------------------------
# 2. Rosenblatt transform diagnostic plot
# -----------------------------------------------------------------------------

plot_rosenblatt_diagnostic <- function(e1, e2, index_name, output_dir) {
  
  e1_ok <- e1[is.finite(e2)]
  e2_ok <- e2[is.finite(e2)]
  
  p1 <- ggplot(data.frame(e2 = e2_ok), aes(x = e2)) +
    geom_histogram(aes(y = ..density..), bins = 30,
                   fill = "steelblue", alpha = 0.7) +
    stat_function(fun = dunif, args = list(min = 0, max = 1),
                  colour = "red", linewidth = 1) +
    labs(title = "e2 Histogram vs Uniform[0,1]", x = "e2 value", y = "Density") +
    theme_bw(base_size = 10)
  
  p2 <- ggplot(data.frame(e1 = e1_ok, e2 = e2_ok), aes(sample = e2)) +
    stat_qq(distribution = qunif) +
    stat_qq_line(distribution = qunif, colour = "red") +
    labs(title = "e2 QQ-Plot",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_bw(base_size = 10)
  
  p3 <- ggplot(data.frame(e1 = e1_ok, e2 = e2_ok), aes(e1, e2)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
    labs(title = "e1 vs e2 (independence check)", x = "e1", y = "e2") +
    theme_bw(base_size = 10)
  
  p <- gridExtra::grid.arrange(p1, p2, p3, ncol = 2, top = index_name)
  ggsave(
    file.path(output_dir, sprintf("%s_rosenblatt_diagnostic.pdf", index_name)),
    p, width = 10, height = 7
  )
  invisible(list(p1 = p1, p2 = p2, p3 = p3))
}


# -----------------------------------------------------------------------------
# 3. Rolling Kendall τ time series plot
# -----------------------------------------------------------------------------

plot_rolling_tau <- function(dc, index_name, output_dir, window_yrs = 10) {
  
  valid_idx <- which(dc$severity > 0 & dc$area_pct > 0)
  if (length(valid_idx) < 20) return(NULL)
  
  yr   <- dc$year[valid_idx]
  sev  <- dc$severity[valid_idx]
  area <- dc$area_pct[valid_idx]
  
  half_w  <- max(floor(window_yrs / 2), 3)
  tau_ser <- rep(NA_real_, length(yr))
  
  for (i in seq_along(yr)) {
    idx <- which(abs(yr - yr[i]) <= half_w)
    if (length(idx) >= 6 && sd(sev[idx]) > 1e-10 && sd(area[idx]) > 1e-10) {
      tau_ser[i] <- cor(sev[idx], area[idx], method = "kendall")
    }
  }
  
  df        <- data.frame(year = yr, tau = tau_ser)
  df_smooth <- df[!is.na(df$tau), ]
  
  p <- ggplot(df_smooth, aes(year, tau)) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "gray50") +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.2, colour = "darkred") +
    labs(
      title = paste(index_name, "Rolling Kendall \u03C4 (", window_yrs, "-yr window)"),
      x = "Year", y = "Kendall \u03C4"
    ) +
    theme_bw(base_size = 11)
  
  ggsave(file.path(output_dir, sprintf("%s_rolling_tau.pdf", index_name)),
         p, width = 8, height = 4)
  invisible(p)
}


# -----------------------------------------------------------------------------
# 4. TV copula parameter evolution plot
# -----------------------------------------------------------------------------

plot_tv_parameter <- function(tv_result, index_name, output_dir) {
  
  if (is.null(tv_result) || tv_result$selected_model != "tv") return(NULL)
  
  yr_seq     <- seq(
    min(tv_result$yr_std * tv_result$yr_sd + tv_result$yr_mean),
    max(tv_result$yr_std * tv_result$yr_sd + tv_result$yr_mean),
    length.out = 100
  )
  yr_std_seq        <- (yr_seq - tv_result$yr_mean) / tv_result$yr_sd
  tv_result$yr_std  <- yr_std_seq
  tv_result$theta_slope <- tv_result$b_hat
  theta_raw <- tv_result$a_hat + tv_result$theta_slope * yr_std_seq
  theta     <- tv_result$link_ilink(theta_raw)
  
  p <- ggplot(data.frame(year = yr_seq, theta = theta), aes(year, theta)) +
    geom_line(colour = "darkgreen", linewidth = 1) +
    geom_hline(
      yintercept = tv_result$link_ilink(tv_result$a_hat),
      linetype = "dashed", colour = "gray60"
    ) +
    labs(
      title = paste(index_name, "TV Copula Parameter Evolution"),
      x = "Year", y = expression(theta)
    ) +
    theme_bw(base_size = 11)
  
  ggsave(file.path(output_dir, sprintf("%s_tv_parameter.pdf", index_name)),
         p, width = 7, height = 4)
  invisible(p)
}


# =============================================================================
log_event("Extensions 1 (v6) loaded: test_dependency | gof_copula | detect_copula_changepoints | rosenblatt_pit_gof | fit_timevarying_copula | fit_segmented_copula")