#============================================================================
# DROUGHT SEVERITY-AREA-FREQUENCY ANALYSIS FOR NECHAKO BASIN
# SPI-1, SPEI-1, SPI-3, SPEI-3  |  Duration classification  |  Dual SAF methods
#
# FIXES vs original:
#   - derive_SAF_curves_kendall now uses the empirical Kendall distribution
#     K_C(t) = P(C(U,V) <= t) instead of the incorrect t-(1-C(t,t))/t formula
#   - fit_copulas now stores LogLik in the comparison CSV (was missing)
#   - Minimum per-class event count raised from 10 to 15
#============================================================================
# NOTE: rm(list=ls()) removed intentionally.  When sourced via RUN_ALL.R this
# script runs AFTER QuickStart.R; a blanket rm() would destroy quick_diagnostic
# and any other objects already in the workspace.  Cleanup is the caller's
# responsibility.  If running standalone, call rm(list=ls()) manually first.

# 0. SETUP & LIBRARIES -------------------------------------------------------
pkgs <- c("terra", "copula", "fitdistrplus", "MASS",
          "ggplot2", "gridExtra", "viridis", "dplyr", "scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

setwd("D:/Nechako_Drought/Nechako")
DROUGHT_THRESHOLD <- -0.5
OUT_ROOT          <- "drought_analysis"
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUT_ROOT, "SAF_analysis_log.txt")
cat(paste(rep("=", 70), collapse = ""), "\n",
    "DROUGHT S-A-F ANALYSIS — NECHAKO BASIN\n",
    "Analysis started:", format(Sys.time()), "\n",
    paste(rep("=", 70), collapse = ""), "\n\n",
    file = LOG_FILE)

log_event <- function(msg) {
  ts   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0(ts, " | ", msg, "\n")
  cat(line, file = LOG_FILE, append = TRUE)
  message(line, appendLF = FALSE)
}

# 1. DURATION CLASSIFICATION -------------------------------------------------
classify_duration_class <- function(dur, scale = 1) {
  if (scale == 1) {
    dplyr::case_when(
      dur >= 1  & dur <= 2  ~ "D1-2 (Very-short)",
      dur >= 3  & dur <= 6  ~ "D3-6 (Short-term)",
      dur >= 7  & dur <= 12 ~ "D7-12 (Medium-term)",
      dur >= 13             ~ "D13+ (Long-term)",   # FIX: label corrected from "D12+" to match condition (dur >= 13)
      TRUE                  ~ "D0 (Sub-threshold)")
  } else {
    dplyr::case_when(
      dur >= 3  & dur <= 6  ~ "D3-6 (Short-term)",
      dur >= 7  & dur <= 12 ~ "D7-12 (Medium-term)",
      dur >= 13             ~ "D13+ (Long-term)",   # FIX: label corrected from "D12+" to match condition (dur >= 13)
      TRUE                  ~ "D0 (Sub-threshold)")
  }
}

# 2. DATA LOADING ------------------------------------------------------------
assemble_index_raster <- function(nc_dir, scale_tag, month_names, prefix = "spi") {
  lst <- lapply(1:12, function(m) {
    f <- file.path(nc_dir,
                   sprintf("%s_%s_month%02d_%s.nc", prefix, scale_tag, m, month_names[m]))
    if (!file.exists(f)) stop("File not found: ", f)
    terra::rast(f)
  })
  r_all <- do.call(c, lst)
  r_all[[order(terra::time(r_all))]]
}

month_names   <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
spi1          <- assemble_index_raster("spi_results_seasonal",  "01", month_names, "spi")
spei1         <- assemble_index_raster("spei_results_seasonal", "01", month_names, "spei")
spi3          <- assemble_index_raster("spi_results_seasonal",  "03", month_names, "spi")
spei3         <- assemble_index_raster("spei_results_seasonal", "03", month_names, "spei")

dates         <- as.Date(terra::time(spi1))
n_months      <- length(dates)
month_indices <- as.integer(format(dates, "%m"))
year_indices  <- as.integer(format(dates, "%Y"))
record_len    <- max(year_indices) - min(year_indices) + 1L
# terra::cellSize() computes geodesic cell areas directly on the reference
# ellipsoid (WGS84 by default for geographic CRS, or native projected units
# converted to km² for projected CRS).  An equal-area / isometric projection
# is NOT required: terra integrates the latitude band for geographic CRS and
# uses the projected-unit area for metric CRS.  Derived from SPI-1 layer 1;
# all four index grids share the same grid definition so this is computed once.
cell_area_km2  <- terra::cellSize(spi1[[1]], unit = "km")
total_area_km2 <- terra::global(cell_area_km2, "sum", na.rm = TRUE)[1, 1]

log_event(sprintf("Loaded indices. Period: %s to %s (%d months). Basin area: %.2f km²",
                  min(dates), max(dates), n_months, total_area_km2))

# 3. DROUGHT CHARACTERISTICS -------------------------------------------------
extract_drought_characteristics_spi <- function(index_rast, cell_area, threshold, index_name) {
  log_event(sprintf("[%s] Extracting characteristics (threshold=%.2f)...", index_name, threshold))
  n_time    <- terra::nlyr(index_rast)
  severity  <- numeric(n_time)
  area_pct  <- numeric(n_time)
  pb <- txtProgressBar(1, n_time, style = 3)
  for (t in 1:n_time) {
    v     <- terra::values(index_rast[[t]], mat = FALSE)
    valid <- is.finite(v)
    dry   <- valid & v < threshold
    if (sum(dry) > 0 && sum(valid) > 0) {
      area_pct[t] <- 100 * sum(dry) / sum(valid)
      severity[t] <- mean(threshold - v[dry], na.rm = TRUE)
    }
    setTxtProgressBar(pb, t)
  }
  close(pb)
  data.frame(date     = dates,
             year     = year_indices,
             month    = month_indices,
             severity = severity,
             area_pct = area_pct)
}

# 4. EVENT IDENTIFICATION & RETURN PERIODS -----------------------------------
identify_duration_classified_events <- function(drought_chars, index_name, scale = 1) {
  log_event(sprintf("[%s] Identifying duration-classified events...", index_name))
  dc   <- drought_chars[order(drought_chars$date), ]
  vals <- dc$severity; area <- dc$area_pct
  indr <- vals > 0 & area > 0
  events <- list(); in_d <- FALSE; s_idx <- NA_integer_
  
  for (i in seq_len(nrow(dc))) {
    if (!in_d && indr[i]) {
      in_d <- TRUE; s_idx <- i
    } else if (in_d && !indr[i]) {
      e_idx <- i - 1L; dur <- e_idx - s_idx + 1L
      events[[length(events) + 1L]] <- data.frame(
        index            = index_name,
        start_date       = dc$date[s_idx],  end_date         = dc$date[e_idx],
        start_year       = dc$year[s_idx],  end_year         = dc$year[e_idx],
        duration_months  = dur,
        duration_class   = classify_duration_class(dur, scale),
        mean_severity    = mean(vals[s_idx:e_idx], na.rm = TRUE),
        max_severity     = max(vals[s_idx:e_idx],  na.rm = TRUE),
        mean_area_pct    = mean(area[s_idx:e_idx],  na.rm = TRUE),
        max_area_pct     = max(area[s_idx:e_idx],   na.rm = TRUE),
        stringsAsFactors = FALSE)
      in_d <- FALSE; s_idx <- NA_integer_
    }
  }
  # Handle event still open at end of record
  if (in_d && !is.na(s_idx)) {
    e_idx <- nrow(dc); dur <- e_idx - s_idx + 1L
    events[[length(events) + 1L]] <- data.frame(
      index            = index_name,
      start_date       = dc$date[s_idx],  end_date         = dc$date[e_idx],
      start_year       = dc$year[s_idx],  end_year         = dc$year[e_idx],
      duration_months  = dur,
      duration_class   = classify_duration_class(dur, scale),
      mean_severity    = mean(vals[s_idx:e_idx], na.rm = TRUE),
      max_severity     = max(vals[s_idx:e_idx],  na.rm = TRUE),
      mean_area_pct    = mean(area[s_idx:e_idx],  na.rm = TRUE),
      max_area_pct     = max(area[s_idx:e_idx],   na.rm = TRUE),
      stringsAsFactors = FALSE)
  }
  if (length(events) == 0) {
    log_event(sprintf("[%s] No events detected.", index_name)); return(data.frame())
  }
  out <- do.call(rbind, events)
  log_event(sprintf("[%s] %d events identified.", index_name, nrow(out)))
  out
}

compute_return_periods_by_class <- function(events_df, index_name, record_years, scale = 1) {
  classes <- if (scale == 1)
    c("D1-2 (Very-short)", "D3-6 (Short-term)", "D7-12 (Medium-term)", "D13+ (Long-term)")
  else
    c("D3-6 (Short-term)", "D7-12 (Medium-term)", "D13+ (Long-term)")
  
  out <- lapply(classes, function(cl) {
    sub <- events_df[events_df$duration_class == cl, ]; n <- nrow(sub)
    if (n == 0)
      return(data.frame(index = index_name, duration_class = cl, n_events = 0L,
                        record_years = record_years, freq_per_year = 0,
                        return_period_years = Inf, mean_duration_months = NA_real_,
                        mean_severity = NA_real_, max_severity = NA_real_,
                        mean_area_pct = NA_real_, max_area_pct = NA_real_,
                        stringsAsFactors = FALSE))
    data.frame(index = index_name, duration_class = cl, n_events = n,
               record_years = record_years,
               freq_per_year = n / record_years,
               return_period_years = record_years / n,
               mean_duration_months = mean(sub$duration_months, na.rm = TRUE),
               mean_severity        = mean(sub$mean_severity,   na.rm = TRUE),
               max_severity         = max(sub$max_severity,     na.rm = TRUE),
               mean_area_pct        = mean(sub$mean_area_pct,   na.rm = TRUE),
               max_area_pct         = max(sub$max_area_pct,     na.rm = TRUE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# 5. MARGINAL & COPULA FITTING -----------------------------------------------
fit_marginal_distributions <- function(drought_data, index_name, out_dir) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  if (nrow(dc) < 10) return(NULL)
  A        <- pmax(pmin(dc$area_pct / 100, 1 - 1e-6), 1e-6)
  fit_beta <- tryCatch(fitdistrplus::fitdist(A, "beta", method = "mle"), error = function(e) NULL)
  if (is.null(fit_beta)) return(NULL)
  
  dists <- list(list("Exponential","exp",NULL),
                list("Gamma",      "gamma",  list(shape=1, rate=1)),
                list("Weibull",    "weibull", list(shape=1, scale=1)),
                list("Log-Normal", "lnorm",  NULL),
                list("Normal",     "norm",   NULL),
                list("Logistic",   "logis",  NULL))
  res  <- data.frame(Distribution=character(), AIC=numeric(), BIC=numeric(),
                     stringsAsFactors=FALSE)
  fits <- list()
  for (d in dists) {
    fit <- tryCatch(
      if (is.null(d[[3]]))
        fitdistrplus::fitdist(dc$severity, d[[2]], method = "mle")
      else
        fitdistrplus::fitdist(dc$severity, d[[2]], method = "mle", start = d[[3]]),
      error = function(e) NULL)
    if (!is.null(fit)) {
      res  <- rbind(res, data.frame(Distribution = d[[1]], AIC = fit$aic, BIC = fit$bic))
      fits[[d[[1]]]] <- fit
    }
  }
  if (nrow(res) == 0) return(NULL)
  best <- res$Distribution[which.min(res$AIC)]
  write.csv(res, file.path(out_dir, sprintf("%s_severity_dist_comparison.csv", index_name)),
            row.names = FALSE)
  list(area_fit           = fit_beta,
       severity_fit       = fits[[best]],
       severity_dist_name = best,
       severity_results   = res)
}

fit_copulas <- function(drought_data, marginal_fits, index_name, out_dir) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  if (nrow(dc) < 10 || is.null(marginal_fits)) return(NULL)
  
  A    <- dc$area_pct / 100
  uA   <- pbeta(A, shape1 = marginal_fits$area_fit$estimate["shape1"],
                shape2 = marginal_fits$area_fit$estimate["shape2"])
  dist <- marginal_fits$severity_dist_name
  pars <- marginal_fits$severity_fit$estimate
  uS   <- switch(dist,
                 Exponential = pexp(dc$severity, rate   = pars["rate"]),
                 Gamma       = pgamma(dc$severity, shape = pars["shape"], rate  = pars["rate"]),
                 Weibull     = pweibull(dc$severity, shape= pars["shape"], scale = pars["scale"]),
                 `Log-Normal`= plnorm(dc$severity, meanlog=pars["meanlog"], sdlog=pars["sdlog"]),
                 Normal      = pnorm(dc$severity, mean  = pars["mean"],  sd    = pars["sd"]),
                 Logistic    = plogis(dc$severity, location=pars["location"], scale=pars["scale"]),
                 pgamma(dc$severity, shape = 1, rate = 1))
  uS  <- pmax(pmin(uS, 1 - 1e-6), 1e-6)
  uA  <- pmax(pmin(uA, 1 - 1e-6), 1e-6)
  uM  <- cbind(uS, uA); n <- nrow(uM)
  
  # Copula candidate set:
  #   Frank       — symmetric, no tail dependence (Archimedean)
  #   Gumbel      — upper-tail dependence (Archimedean)
  #   Plackett    — symmetric, light tails (one-parameter)
  #   SurvClayton — survival (180°-rotated) Clayton: upper-tail dependence,
  #                 complementary to standard Clayton's lower-tail focus
  #   Joe         — strong upper-tail dependence (Archimedean)
  # Standard Clayton and Gaussian Normal have been intentionally excluded.
  cops <- list(
    Frank       = frankCopula(dim = 2),
    Gumbel      = gumbelCopula(dim = 2),
    Plackett    = plackettCopula(),
    SurvClayton = rotCopula(claytonCopula(dim = 2)),
    Joe         = joeCopula(dim = 2),
    # New robust candidates:
    Gaussian    = normalCopula(param = 0.5, dim = 2, dispstr = "un"),
    StudentT    = tCopula(param = 0.5, df = 4, dim = 2)
  )
  res  <- data.frame(Copula=character(), Parameter=numeric(),
                     LogLik=numeric(), AIC=numeric(), BIC=numeric(),
                     stringsAsFactors=FALSE)
  fits <- list()
  for (nm in names(cops)) {
    # "mpl" here maximises the copula log-likelihood given the parametric CDF
    # transforms uS/uA (IFM-style two-stage MLE, not rank-based pseudo-obs).
    # FIX: added optim.method="BFGS" + maxit=1000 to suppress convergence
    #      code=52 warnings; falls back to Nelder-Mead if BFGS errors out.
    # NOTE: for SurvClayton (rotCopula), fitCopula / loglikCopula dispatch
    #       through dCopula which is generic and correctly handles rotCopula
    #       objects; coef() returns the base Clayton parameter.
    fit <- tryCatch(
      fitCopula(cops[[nm]], uM, method = "mpl",
                optim.method  = "BFGS",
                optim.control = list(maxit = 1000, reltol = 1e-8)),
      error = function(e)
        tryCatch(
          fitCopula(cops[[nm]], uM, method = "mpl",
                    optim.method  = "Nelder-Mead",
                    optim.control = list(maxit = 2000)),
          error = function(e2) NULL))
    if (!is.null(fit)) {
      p  <- coef(fit)
      ll <- loglikCopula(p, uM, cops[[nm]])
      k  <- length(p)
      res  <- rbind(res, data.frame(Copula    = nm,
                                    Parameter = p[1],
                                    LogLik    = ll,
                                    AIC       = -2*ll + 2*k,
                                    BIC       = -2*ll + log(n)*k))
      fits[[nm]] <- fit
    }
  }
  if (nrow(res) == 0) return(NULL)
  best <- res$Copula[which.min(res$AIC)]
  log_event(sprintf("[%s] Best copula: %s (AIC=%.2f)", index_name, best,
                    res$AIC[res$Copula == best]))
  
  # ---- Empirical upper-tail dependence coefficient ---------------------------
  # λ_U = lim_{u→1} P(V > u | U > u)
  # Estimated non-parametrically (Schmidt & Stadtmüller 2006, Scand. J. Stat.):
  #   λ_U ≈ 2 - log(C_hi) / log(u_q_hi)   at  u_q_hi = 0.95
  # where C_hi = empirical P(U > 0.95, V > 0.95).
  # Result is clamped to [0, 1].
  #
  # Lower-tail dependence (λ_L) is NOT assessed.  Drought severity and
  # drought-affected area co-occur at the UPPER tail (joint extremes), so
  # lower-tail diagnostics are irrelevant to the SAF objective.
  #
  # References: Joe (1997); Nelsen (2006); Schmidt & Stadtmüller (2006).
  emp_td <- tryCatch({
    u_q_hi <- 0.95
    C_hi   <- mean(uS >= u_q_hi & uA >= u_q_hi)   # empirical P(U > 0.95, V > 0.95)
    C_hi   <- max(C_hi, 1e-8)                       # boundary guard: log(0) protection
    lam_U_emp <- max(0, min(1, 2 - log(C_hi) / log(u_q_hi)))
    list(lambda_U = lam_U_emp)
  }, error = function(e) list(lambda_U = NA_real_))
  
  # Theoretical upper-tail dependence coefficient for the best-fitting copula
  # computed at the MLE parameter value.
  best_par   <- coef(fits[[best]])[1]   # rho for StudentT, scalar param otherwise
  theo_td    <- tryCatch({
    switch(best,
           Gumbel      = list(lambda_U = 2 - 2^(1/max(best_par, 1.001))),
           Joe         = list(lambda_U = 2 - 2^(1/max(best_par, 1.001))),
           SurvClayton = list(lambda_U = 2^(-1/max(best_par, 1e-6))),
           Frank       = list(lambda_U = 0),
           Plackett    = list(lambda_U = 0),
           Gaussian    = list(lambda_U = 0),   # Gaussian: asymptotic tail independence
           StudentT    = {                      # λ_U = 2*t_{df+1}(-sqrt((df+1)(1-ρ)/(1+ρ)))
             rho    <- max(min(best_par, 0.9999), -0.9999)
             df_par <- coef(fits[[best]])[2]   # degrees-of-freedom parameter
             df_par <- max(df_par, 1.5)
             list(lambda_U = 2 * pt(-sqrt((df_par + 1) * (1 - rho) / (1 + rho)),
                                    df = df_par + 1))
           },
           list(lambda_U = NA_real_)            # unknown family fallback
    )
  }, error = function(e) list(lambda_U = NA_real_))
  
  # ---- Copula family recommendation based on empirical upper-tail evidence ---
  # Decision rule (Joe 1997; Nelsen 2006; Poulin et al. 2007):
  #   λ_U > 0.10  → upper-tail dependence present → prefer Gumbel / Joe / SurvClayton
  #   λ_U ≤ 0.10  → no meaningful upper-tail dep  → Frank / Plackett adequate
  # Threshold 0.10 is conservative; λ_U < 0.05 is strong evidence of asymptotic
  # tail independence.  The recommendation is advisory — AIC winner is still used.
  td_recommendation <- if (!is.finite(emp_td$lambda_U)) {
    "Tail dependence inconclusive (insufficient data)"
  } else if (emp_td$lambda_U > 0.10) {
    "Upper-tail dependence detected (λ_U > 0.10): Gumbel, Joe, or SurvClayton preferred"
  } else {
    "No meaningful upper-tail dependence (λ_U ≤ 0.10): Frank or Plackett adequate"
  }
  
  # Families with zero theoretical upper-tail dependence.
  # Frank, Plackett, and Gaussian all have λ_U = 0; if any of these wins
  # AIC but empirical λ_U > 0.10, the selected model cannot represent
  # joint tail extremes and should be flagged.
  zero_tail_families <- c("Frank", "Plackett", "Gaussian")
  td_mismatch <- isTRUE(emp_td$lambda_U > 0.10) &&
    best %in% zero_tail_families &&
    !is.na(emp_td$lambda_U)
  
  log_event(sprintf(
    "[%s] Tail dependence — empirical λ_U=%.4f | theoretical (best=%s) λ_U=%.4f",
    index_name, emp_td$lambda_U, best,
    if (is.finite(theo_td$lambda_U)) theo_td$lambda_U else NA))
  log_event(sprintf("[%s] TD recommendation: %s", index_name, td_recommendation))
  if (td_mismatch)
    log_event(sprintf(
      "[%s] WARNING [TD-U]: AIC-best copula (%s) has no upper-tail dependence (λ_U = 0) but empirical λ_U=%.4f > 0.10. Consider Gumbel, Joe, SurvClayton, or StudentT for extreme return periods.",
      index_name, best, emp_td$lambda_U))
  
  # Append tail dependence columns to the comparison table
  res$emp_lambda_U   <- round(emp_td$lambda_U, 4)
  # Theoretical λ_U per candidate (requires parameter from each fit)
  # StudentT: λ_U = 2·t_{df+1}(−√((df+1)(1−ρ)/(1+ρ)))   (Joe 1997 §2.3)
  # Gaussian: λ_U = 0 (asymptotic tail independence for all |ρ| < 1)
  # Frank/Plackett: λ_U = 0
  theo_U <- vapply(names(cops), function(nm) {
    if (is.null(fits[[nm]])) return(NA_real_)
    p <- coef(fits[[nm]])
    switch(nm,
           Gumbel      = 2 - 2^(1/max(p[1], 1.001)),
           Joe         = 2 - 2^(1/max(p[1], 1.001)),
           SurvClayton = 2^(-1/max(p[1], 1e-6)),
           Frank       = 0,
           Plackett    = 0,
           Gaussian    = 0,
           StudentT    = {
             rho    <- max(min(p[1], 0.9999), -0.9999)
             df_par <- max(if (length(p) >= 2L) p[2] else 4, 1.5)
             2 * pt(-sqrt((df_par + 1) * (1 - rho) / (1 + rho)), df = df_par + 1)
           },
           NA_real_)    # unknown family fallback
  }, numeric(1L))
  res$theo_lambda_U    <- round(theo_U, 4)
  res$td_mismatch      <- td_mismatch & (res$Copula == best)
  res$td_recommendation <- td_recommendation
  
  write.csv(res, file.path(out_dir, sprintf("%s_copula_comparison.csv", index_name)),
            row.names = FALSE)
  list(best_copula_name    = best,
       best_copula_fit     = fits[[best]],
       all_results         = res,
       u_severity          = uS,
       u_area              = uA,
       emp_lambda_U        = emp_td$lambda_U,
       theo_lambda_U_best  = theo_td$lambda_U,
       td_mismatch         = td_mismatch,
       td_recommendation   = td_recommendation)
}

# 6. HELPER: CDF of severity under fitted marginal ---------------------------
compute_u_severity <- function(s_val, marginal_fits) {
  dist <- marginal_fits$severity_dist_name
  pars <- marginal_fits$severity_fit$estimate
  switch(dist,
         Exponential = pexp(s_val,  rate     = pars["rate"]),
         Gamma       = pgamma(s_val, shape   = pars["shape"],  rate  = pars["rate"]),
         Weibull     = pweibull(s_val, shape = pars["shape"],  scale = pars["scale"]),
         `Log-Normal`= plnorm(s_val, meanlog = pars["meanlog"], sdlog= pars["sdlog"]),
         Normal      = pnorm(s_val,  mean    = pars["mean"],   sd    = pars["sd"]),
         Logistic    = plogis(s_val, location= pars["location"],scale = pars["scale"]),
         pgamma(s_val, shape = 1, rate = 1))
}

# 7. SAF CURVE DERIVATION ----------------------------------------------------

# Method A: Conditional copula  P(S > s | A = a) = 1 - T_eff
derive_SAF_curves_conditional <- function(drought_data, marginal_fits, copula_fit,
                                          index_name, out_dir, duration_class = "all") {
  if (is.null(marginal_fits) || is.null(copula_fit)) return(NULL)
  log_event(sprintf("[%s | %s] SAF — Conditional Copula...", index_name, duration_class))
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data); n_dr <- nrow(dc); if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  
  cop_obj  <- copula_fit$best_copula_fit@copula
  beta_par <- marginal_fits$area_fit$estimate
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5); area_prop <- area_pct / 100
  out <- data.frame()
  
  for (T in T_years) {
    target <- 1 - mu_T / (T * 12)
    if (target <= 0 || target >= 1) next
    sev <- numeric(length(area_prop))
    for (i in seq_along(area_prop)) {
      v   <- pbeta(area_prop[i], shape1 = beta_par[1], shape2 = beta_par[2])
      obj <- function(s) {
        u <- compute_u_severity(s, marginal_fits)
        (copula::cCopula(cbind(u, v), copula = cop_obj, indices = 2) - target)^2
      }
      res    <- tryCatch(optimize(obj, interval = c(1e-6, 50)),
                         error = function(e) list(minimum = NA_real_))
      sev[i] <- res$minimum
    }
    out <- rbind(out, data.frame(ReturnPeriod_years  = T,
                                 ReturnPeriod_months = T * 12,
                                 Area_pct            = area_pct,
                                 Severity            = sev,
                                 Method              = "Conditional",
                                 Index               = index_name,
                                 DurationClass       = duration_class))
  }
  write.csv(out,
            file.path(out_dir, sprintf("%s_SAF_conditional_%s.csv", index_name,
                                       gsub("[^A-Za-z0-9]", "_", duration_class))),
            row.names = FALSE)
  out
}

# Method B: Kendall distribution  K_C(t_K) = 1 - mu_T / (T*12)
# FIX: uses empirical K_C(t) = P(C(U,V) <= t) estimated from the fitted
#      pseudo-observations, replacing the incorrect t-(1-C(t,t))/t formula.
derive_SAF_curves_kendall <- function(drought_data, marginal_fits, copula_fit,
                                      index_name, out_dir, duration_class = "all") {
  if (is.null(marginal_fits) || is.null(copula_fit)) return(NULL)
  log_event(sprintf("[%s | %s] SAF — Kendall Distribution (empirical K_C)...",
                    index_name, duration_class))
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data); n_dr <- nrow(dc); if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  
  cop_obj  <- copula_fit$best_copula_fit@copula
  beta_par <- marginal_fits$area_fit$estimate
  uS       <- copula_fit$u_severity
  uA       <- copula_fit$u_area
  
  # Empirical Kendall distribution: K_C(t) = P(C(U,V) <= t)
  cop_vals <- copula::pCopula(cbind(uS, uA), cop_obj)
  kc_fn    <- function(t) mean(cop_vals <= t)   # scalar empirical K_C(t): called one value at a time by optimize()
  
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5); area_prop <- area_pct / 100
  out <- data.frame()
  
  for (T in T_years) {
    target_kc <- 1 - mu_T / (T * 12)
    if (target_kc <= 0 || target_kc >= 1) next
    
    # Solve K_C(t_val) = target_kc
    kc_obj <- function(t) (kc_fn(t) - target_kc)^2
    t_sol  <- tryCatch(optimize(kc_obj, interval = c(1e-4, 1 - 1e-4)),
                       error = function(e) list(minimum = 0.5))
    t_val  <- t_sol$minimum
    
    sev <- numeric(length(area_prop))
    for (i in seq_along(area_prop)) {
      v   <- pbeta(area_prop[i], shape1 = beta_par[1], shape2 = beta_par[2])
      obj <- function(s) {
        u   <- compute_u_severity(s, marginal_fits)
        val <- tryCatch(copula::pCopula(cbind(u, v), cop_obj), error = function(e) NA_real_)
        if (!is.finite(val)) return(1e6)   # guard: pCopula can return NA/NaN for some copulas
        (val - t_val)^2
      }
      res    <- tryCatch(optimize(obj, interval = c(1e-6, 50)),
                         error = function(e) list(minimum = NA_real_))
      sev[i] <- res$minimum
    }
    out <- rbind(out, data.frame(ReturnPeriod_years  = T,
                                 ReturnPeriod_months = T * 12,
                                 Area_pct            = area_pct,
                                 Severity            = sev,
                                 Method              = "Kendall",
                                 Index               = index_name,
                                 DurationClass       = duration_class))
  }
  write.csv(out,
            file.path(out_dir, sprintf("%s_SAF_kendall_%s.csv", index_name,
                                       gsub("[^A-Za-z0-9]", "_", duration_class))),
            row.names = FALSE)
  out
}

create_comparative_SAF_plots <- function(saf_cond, saf_kend, index_name, out_dir) {
  if (is.null(saf_cond) || is.null(saf_kend)) return(invisible(NULL))
  comb <- rbind(saf_cond, saf_kend)
  # FIX: remove NA/Inf rows before plotting to suppress ggplot "removed rows"
  #      warning.  These arise at the tail of the SAF curve where the
  #      conditional probability hits 0 or 1 (extrapolation boundary).
  comb <- comb[is.finite(comb$Area_pct) & is.finite(comb$Severity), ]
  # Derive axis limits from the finite data so the boundary is intentional
  max_rp <- if (nrow(comb) > 0) max(comb$Area_pct, na.rm = TRUE) else 100
  p <- ggplot(comb, aes(x = Area_pct, y = Severity,
                        color    = factor(ReturnPeriod_years),
                        linetype = Method, shape = Method)) +
    geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
    scale_color_viridis_d(name = "Return Period\n(years)") +
    scale_linetype_manual(values = c(Conditional = "solid", Kendall = "dashed"),
                          name = "Method") +
    scale_shape_manual(values = c(Conditional = 19, Kendall = 17), name = "Method") +
    scale_x_continuous(limits = c(0, max_rp),
                       oob    = scales::oob_keep) +
    labs(x = "Percent of Area Under Drought (%)", y = "Drought Severity",
         title = sprintf("%s: SAF Curves (Conditional vs Kendall)", index_name)) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  pdf(file.path(out_dir, sprintf("%s_SAF_overlay.pdf", index_name)),
      width = 10, height = 8)
  print(p); dev.off()
}

subset_drought_chars_by_class <- function(drought_chars, events_df, class_label) {
  ev <- events_df[events_df$duration_class == class_label, ]
  if (nrow(ev) == 0) return(NULL)
  in_class <- logical(nrow(drought_chars))
  for (i in seq_len(nrow(ev)))
    in_class <- in_class |
    (drought_chars$date >= ev$start_date[i] & drought_chars$date <= ev$end_date[i])
  dc <- drought_chars
  dc$severity[!in_class] <- 0
  dc$area_pct[!in_class] <- 0
  dc
}

# 8. MASTER PIPELINE ---------------------------------------------------------
run_saf_pipeline <- function(index_name, index_rast, out_dir, scale = 1) {
  log_event(paste(rep("=", 70), collapse = ""))
  log_event(sprintf("RUNNING SAF PIPELINE: %s", index_name))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # a) Drought characteristics
  dc <- extract_drought_characteristics_spi(index_rast, cell_area_km2,
                                            DROUGHT_THRESHOLD, index_name)
  write.csv(dc,
            file.path(out_dir, sprintf("%s_drought_characteristics.csv", index_name)),
            row.names = FALSE)
  
  # b) Events & empirical return periods
  events <- identify_duration_classified_events(dc, index_name, scale)
  if (nrow(events) > 0)
    write.csv(events,
              file.path(out_dir, sprintf("%s_duration_classified_events.csv", index_name)),
              row.names = FALSE)
  rp <- compute_return_periods_by_class(events, index_name, record_len, scale)
  write.csv(rp,
            file.path(out_dir, sprintf("%s_return_periods_by_class.csv", index_name)),
            row.names = FALSE)
  
  # c) Marginals & copula
  m <- fit_marginal_distributions(dc, index_name, out_dir)
  if (is.null(m)) { log_event(sprintf("[%s] Marginal fitting failed. Skipping.", index_name)); return(invisible(NULL)) }
  cop <- fit_copulas(dc, m, index_name, out_dir)
  if (is.null(cop)) { log_event(sprintf("[%s] Copula fitting failed. Skipping.", index_name)); return(invisible(NULL)) }
  
  # d) Full-record SAF curves (both methods)
  saf_cond <- derive_SAF_curves_conditional(dc, m, cop, index_name, out_dir, "all")
  saf_kend <- derive_SAF_curves_kendall(dc, m, cop, index_name, out_dir, "all")
  create_comparative_SAF_plots(saf_cond, saf_kend, index_name, out_dir)
  
  # e) Per-class SAF (FIX: threshold raised to 15 events for reliable fitting)
  dc_labels       <- if (scale == 1)
    c("D1-2 (Very-short)", "D3-6 (Short-term)", "D7-12 (Medium-term)", "D13+ (Long-term)")
  else
    c("D3-6 (Short-term)", "D7-12 (Medium-term)", "D13+ (Long-term)")
  saf_cond_cls <- list(); saf_kend_cls <- list()
  
  for (cl in dc_labels) {
    n_ev <- sum(events$duration_class == cl, na.rm = TRUE)
    if (n_ev < 15) {
      log_event(sprintf("[%s | %s] Only %d events — skipping class SAF.", index_name, cl, n_ev))
      next
    }
    dc_sub <- subset_drought_chars_by_class(dc, events, cl)
    if (is.null(dc_sub)) next
    lbl    <- paste0(index_name, "_", cl)
    m_sub  <- fit_marginal_distributions(dc_sub, lbl, out_dir)
    c_sub  <- if (!is.null(m_sub)) fit_copulas(dc_sub, m_sub, lbl, out_dir) else NULL
    if (!is.null(m_sub) && !is.null(c_sub)) {
      saf_cond_cls[[cl]] <- derive_SAF_curves_conditional(dc_sub, m_sub, c_sub, index_name, out_dir, cl)
      saf_kend_cls[[cl]] <- derive_SAF_curves_kendall(dc_sub, m_sub, c_sub, index_name, out_dir, cl)
    }
  }
  
  log_event(sprintf("[%s] Pipeline complete.", index_name))
  invisible(list(dc               = dc,
                 events           = events,
                 rp               = rp,
                 marginals        = m,
                 copulas          = cop,
                 saf_conditional  = saf_cond,
                 saf_kendall      = saf_kend,
                 saf_cond_by_class= saf_cond_cls,
                 saf_kend_by_class= saf_kend_cls))
}

# 9. EXECUTE -----------------------------------------------------------------
res_spi1  <- run_saf_pipeline("SPI1",  spi1,  file.path(OUT_ROOT, "SPI1_analysis"),  scale = 1)
res_spei1 <- run_saf_pipeline("SPEI1", spei1, file.path(OUT_ROOT, "SPEI1_analysis"), scale = 1)
res_spi3  <- run_saf_pipeline("SPI3",  spi3,  file.path(OUT_ROOT, "SPI3_analysis"),  scale = 3)
res_spei3 <- run_saf_pipeline("SPEI3", spei3, file.path(OUT_ROOT, "SPEI3_analysis"), scale = 3)

log_event("=== MAIN PIPELINE COMPLETE. Check drought_analysis/ for results. ===")