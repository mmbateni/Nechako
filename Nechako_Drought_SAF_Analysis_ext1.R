#============================================================================
# SAF EXTENSIONS 1 — ENHANCED FUNCTIONS
# Source AFTER Nechako_Drought_SAF_Analysis.R
#
# Provides three new functions that extend the main pipeline:
#   test_dependency()          — Kendall τ / Spearman ρ between severity & area
#   gof_copula()               — Bootstrap goodness-of-fit test (Sn statistic)
#   fit_timevarying_copula()   — Tests if copula dependence changes over time
#
# FIXES vs original ext1:
#   - REMOVED rm(list=ls()) which previously wiped the main pipeline results
#   - REMOVED duplicated functions already defined in SAF_Analysis.R
#   - fit_timevarying_copula now uses the AIC-selected copula family instead
#     of always hardcoding Clayton
#   - TV copula results are saved to CSV
#============================================================================
if (!requireNamespace("moments",  quietly = TRUE)) install.packages("moments")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("Kendall",  quietly = TRUE)) install.packages("Kendall")
library(moments); library(openxlsx); library(Kendall)

# Verify that the main pipeline has been run
if (!exists("log_event"))
  stop("log_event() not found. Source Nechako_Drought_SAF_Analysis.R first.")

# ---------------------------------------------------------------------------
# 1. DEPENDENCY TEST
#    Tests the strength and significance of the severity–area relationship
#    using rank-based correlation (Kendall τ and Spearman ρ).
# ---------------------------------------------------------------------------
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
  log_event(sprintf("  [%s] Dependency — Kendall τ=%.4f (p=%.4e) | Spearman ρ=%.4f",
                    index_name, cor_k, t_res$p.value, cor_s))
  list(kendall  = cor_k,
       spearman = cor_s,
       p_value  = t_res$p.value)
}

# ---------------------------------------------------------------------------
# 2. COPULA GOODNESS-OF-FIT  (Bootstrap Sn test)
#    Uses the parametric bootstrap with N_boot replicates.
#    Decision: p > 0.05 => copula not rejected at 5% level.
# ---------------------------------------------------------------------------
gof_copula <- function(copula_fit_obj, u_matrix, index_name, N_boot = 499L) {
  if (is.null(copula_fit_obj)) return(NULL)
  cop <- copula_fit_obj$best_copula_fit@copula
  # FIX: pass ties=TRUE explicitly so the ties-correction is intentional
  #      rather than auto-detected mid-run (eliminates the informational warning).
  res <- tryCatch(
    suppressWarnings(
      copula::gofCopula(cop, u_matrix, N = N_boot,
                        method = "Sn", optim.method = "BFGS",
                        ties   = TRUE)),
    error = function(e) { log_event(sprintf("  [%s] GOF test error: %s", index_name, e$message)); NULL })
  if (is.null(res)) return(NULL)
  decision <- if (res$p.value > 0.05) "ACCEPTED" else "REJECTED"
  log_event(sprintf("  [%s] Copula GOF (%s): Sn=%.4f, p=%.4f => %s",
                    index_name, copula_fit_obj$best_copula_name,
                    res$statistic, res$p.value, decision))
  list(statistic   = res$statistic,
       p_value     = res$p.value,
       decision    = decision,
       copula_name = copula_fit_obj$best_copula_name)
}

# ---------------------------------------------------------------------------
# 3. TIME-VARYING COPULA
#    Fits a log-linear model theta(t) = ilink(a + b*year_std) and tests
#    whether b != 0 via a likelihood-ratio test (chi-sq, df=1).
#
#    FIX: now uses the AIC-selected copula family from copula_fit_obj instead
#         of always forcing Clayton.  Family-appropriate link functions and
#         parameter bounds are applied automatically.
# ---------------------------------------------------------------------------
fit_timevarying_copula <- function(drought_data, copula_fit_obj, marginal_fits,
                                   year_vec, index_name, output_dir) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  if (nrow(dc) < 20L) {
    log_event(sprintf("  [%s] TV copula skipped: < 20 drought observations.", index_name))
    return(NULL)
  }
  
  # FIX (Bug 1): dc is a filtered subset of drought_data; its years must come
  # from dc$year, not from the first nrow(dc) elements of year_vec (which
  # would assign the wrong calendar years to each drought observation).
  yr_match <- dc$year
  if (length(yr_match) != nrow(dc)) yr_match <- rep(mean(year_vec), nrow(dc))
  year_std <- as.numeric(scale(yr_match))
  
  uS   <- copula_fit_obj$u_severity
  uA   <- copula_fit_obj$u_area
  uDat <- cbind(uS, uA)
  
  cop_name <- copula_fit_obj$best_copula_name
  init_par <- coef(copula_fit_obj$best_copula_fit)[1]
  
  # Family-specific link/inverse-link and safe parameter bounds
  link_info <- switch(cop_name,
                      Clayton  = list(link  = function(x) log(pmax(x, 1e-6)),
                                      ilink = exp,
                                      lo = 1e-3, hi = 50,
                                      init  = max(init_par, 1e-3)),
                      Gumbel   = list(link  = function(x) log(pmax(x - 1, 1e-6)),
                                      ilink = function(x) exp(x) + 1,
                                      lo = 1.001, hi = 50,
                                      init  = max(init_par, 1.001)),
                      Joe      = list(link  = function(x) log(pmax(x - 1, 1e-6)),
                                      ilink = function(x) exp(x) + 1,
                                      lo = 1.001, hi = 50,
                                      init  = max(init_par, 1.001)),
                      Frank    = list(link  = identity,
                                      ilink = identity,
                                      lo = -50, hi = 50,
                                      init  = init_par),
                      Normal   = list(link  = function(x) log((1 + x) / (1 - x)),    # logit of (1+rho)/2
                                      ilink = function(x) (exp(x) - 1) / (exp(x) + 1),
                                      lo = -0.999, hi = 0.999,
                                      init  = max(min(init_par, 0.999), -0.999)),
                      Plackett = list(link  = function(x) log(pmax(x, 1e-6)),
                                      ilink = exp,
                                      lo = 1e-3, hi = 500,
                                      init  = max(init_par, 1e-3)),
                      # default: Clayton
                      list(link = function(x) log(pmax(x, 1e-6)),
                           ilink = exp, lo = 1e-3, hi = 50,
                           init  = max(init_par, 1e-3))
  )
  
  make_cop <- function(theta_raw) {
    th <- max(min(link_info$ilink(theta_raw), link_info$hi), link_info$lo)
    switch(cop_name,
           Clayton  = claytonCopula(th, dim = 2),
           Gumbel   = gumbelCopula(th,  dim = 2),
           Frank    = frankCopula(th,   dim = 2),
           Joe      = joeCopula(th,     dim = 2),
           Normal   = normalCopula(th,  dim = 2),
           Plackett = plackettCopula(th),
           claytonCopula(th, dim = 2))
  }
  
  # Stationary NLL
  nll_stat <- function(p1) {
    cop <- make_cop(p1)
    -sum(log(pmax(copula::dCopula(uDat, cop), 1e-300)))
  }
  
  # Time-varying NLL (row-by-row: each obs has its own theta)
  nll_tv <- function(p) {
    th_raw <- p[1] + p[2] * year_std
    ll <- 0
    for (i in seq_len(nrow(uDat))) {
      cop <- make_cop(th_raw[i])
      ll  <- ll + log(pmax(copula::dCopula(uDat[i, , drop = FALSE], cop), 1e-300))
    }
    -ll
  }
  
  opt_stat <- tryCatch(
    optim(link_info$link(link_info$init), nll_stat,
          method = "Brent", lower = -10, upper = 10),
    error = function(e) NULL)
  if (is.null(opt_stat)) {
    log_event(sprintf("  [%s] TV copula: stationary optimisation failed.", index_name))
    return(NULL)
  }
  opt_tv <- tryCatch(
    optim(c(opt_stat$par, 0), nll_tv, method = "BFGS"),
    error = function(e) NULL)
  if (is.null(opt_tv)) {
    log_event(sprintf("  [%s] TV copula: time-varying optimisation failed.", index_name))
    return(NULL)
  }
  
  lr    <- 2 * (opt_stat$value - opt_tv$value)
  p_val <- pchisq(max(0, lr), df = 1L, lower.tail = FALSE)
  dir   <- if (opt_tv$par[2] > 0) "strengthening" else "weakening"
  
  result <- data.frame(
    index     = index_name,
    copula    = cop_name,
    a_hat     = opt_tv$par[1],
    b_hat     = opt_tv$par[2],
    lr_stat   = round(lr, 4),
    lr_p      = round(p_val, 4),
    direction = dir,
    significant = p_val < 0.05,
    stringsAsFactors = FALSE)
  
  write.csv(result,
            file.path(output_dir, sprintf("%s_timevarying_copula.csv", index_name)),
            row.names = FALSE)
  log_event(sprintf("  [%s] TV copula (%s): b=%.4f, LR p=%.4f (%s%s)",
                    index_name, cop_name, opt_tv$par[2], p_val, dir,
                    if (p_val < 0.05) " *SIGNIFICANT*" else ""))
  result
}

log_event("Extensions 1 loaded: test_dependency | gof_copula | fit_timevarying_copula")