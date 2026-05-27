#============================================================================
# DROUGHT SEVERITY-AREA-FREQUENCY ANALYSIS FOR NECHAKO BASIN
# SPI-1, SPEI-1, SPI-3, SPEI-3  |  Duration classification  |  Dual SAF methods
#
# FIXES :
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
OUT_ROOT          <- "spatial_drought"
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
    sub <- events_df[events_df$duration_class == cl, ]
    n   <- nrow(sub)
    
    if (n == 0)
      return(data.frame(
        index = index_name, duration_class = cl, n_events = 0L,
        record_years = record_years,
        freq_per_year = 0, return_period_years = Inf,
        rp_ci_lo_years = Inf, rp_ci_hi_years = Inf,
        iat_mean_years = NA_real_, iat_cv = NA_real_,
        poisson_assumption = "untestable",
        mu_T_adj_months  = NA_real_,   # ← ADD THIS
        renewal_shape_k  = NA_real_,   # ← ADD THIS
        renewal_rate_k   = NA_real_,   # ← ADD THIS
        mean_duration_months = NA_real_,
        mean_severity = NA_real_, max_severity = NA_real_,
        mean_area_pct = NA_real_, max_area_pct = NA_real_,
        stringsAsFactors = FALSE))
    
    # ---- Point estimate (unchanged) ----------------------------------------
    lambda_hat <- n / record_years          # events per year
    rp_hat     <- record_years / n          # years per event
    
    # ---- Exact Poisson confidence interval on the rate ---------------------
    # Treat n as a Poisson observation with exposure = record_years.
    # 95% CI on lambda: uses relationship between Poisson and chi-squared:
    #   lambda_lo = qchisq(0.025, df=2*n)   / (2 * record_years)
    #   lambda_hi = qchisq(0.975, df=2*(n+1)) / (2 * record_years)
    # Invert to get CI on RP (note: lo lambda -> hi RP and vice versa).
    lambda_lo <- qchisq(0.025, df = 2 * n)       / (2 * record_years)
    lambda_hi <- qchisq(0.975, df = 2 * (n + 1)) / (2 * record_years)
    rp_ci_lo  <- if (lambda_hi > 0) 1 / lambda_hi else Inf   # hi lambda = lo RP
    rp_ci_hi  <- if (lambda_lo > 0) 1 / lambda_lo else Inf   # lo lambda = hi RP
    
    #  ---- Inter-arrival time diagnostics (Gamma Renewal Adjustment) -------------
    # Replaces homogeneous Poisson assumption with empirical renewal process.
    # When iat_cv != 1.0, mu_T and the exceedance probability are adjusted
    # using a Gamma-distributed inter-arrival model (shape k = 1/CV^2).
    iat_mean  <- NA_real_
    iat_cv    <- NA_real_
    poisson_flag  <- "untestable (n < 2)"
    mu_T_adj_months <- NA_real_  # Empirical mean waiting time (months)
    renewal_shape_k <- NA_real_  # Gamma shape parameter
    renewal_rate_k  <- NA_real_  # Gamma rate parameter
    
    if (n >= 2) {
      sub_sorted  <- sub[order(sub$start_date), ]
      # Inter-arrival in years (start-to-start)
      iat_yrs  <- as.numeric(diff(as.Date(sub_sorted$start_date))) / 365.25
      iat_mean  <- mean(iat_yrs, na.rm = TRUE)
      iat_sd    <- sd(iat_yrs,   na.rm = TRUE)
      iat_cv    <- if (iat_mean > 0) iat_sd / iat_mean else NA_real_
      
      # GAMMA RENEWAL PARAMETERS
      # k < 1  => Clustering (over-dispersed, heavy-tailed gaps)
      # k > 1  => Regularity (under-dispersed, uniform gaps)
      # k = 1  => Exponential (Poisson)
      if (!is.na(iat_cv) && iat_cv > 0) {
        renewal_shape_k <- 1 / (iat_cv^2)
        renewal_rate_k  <- renewal_shape_k / iat_mean
        # Empirical mu_T (months) replaces theoretical n_record / n_dr
        mu_T_adj_months <- iat_mean * 12
      } else {
        # Fallback to theoretical Poisson if CV is degenerate
        mu_T_adj_months <- (record_years / n) * 12
        renewal_shape_k <- 1
        renewal_rate_k  <- 1 / iat_mean
      }
      
      # Flag clustering/regularity with process type label
      poisson_flag  <- if (is.na(iat_cv))           "untestable" else
        if (iat_cv > 1.5)             "CLUSTERED (Gamma k<1; use renewal-adjusted SAF target)" else
          if (iat_cv < 0.5)             "REGULAR (Gamma k>1; use renewal-adjusted SAF target)" else
            "approx Poisson (0.5 <= CV <= 1.5)"
    }
    
    data.frame(
      index                = index_name,
      duration_class       = cl,
      n_events             = n,
      record_years         = record_years,
      freq_per_year        = ifelse(n>0, n/record_years, 0),
      return_period_years  = ifelse(n>0, record_years/n, Inf),
      rp_ci_lo_years       = round(ifelse(n>0, 1/qchisq(0.975, df=2*(n+1))/(2*record_years), Inf), 1),
      rp_ci_hi_years       = round(ifelse(n>0, 1/qchisq(0.025, df=2*n)/(2*record_years), Inf), 1),
      iat_mean_years       = round(iat_mean, 2),
      iat_cv               = round(iat_cv, 2),
      poisson_assumption   = poisson_flag,
      mu_T_adj_months      = round(mu_T_adj_months, 2),      # NEW: empirical waiting time
      renewal_shape_k      = round(renewal_shape_k, 4),      # NEW: Gamma shape
      renewal_rate_k       = round(renewal_rate_k, 4),       # NEW: Gamma rate
      mean_duration_months = mean(sub$duration_months, na.rm = TRUE),
      mean_severity        = mean(sub$mean_severity,   na.rm = TRUE),
      max_severity         = max(sub$max_severity,     na.rm = TRUE),
      mean_area_pct        = mean(sub$mean_area_pct,   na.rm = TRUE),
      max_area_pct         = max(sub$max_area_pct,     na.rm = TRUE),
      stringsAsFactors     = FALSE)
  })
  do.call(rbind, out)
}
# ============================================================================
# HELPER: Clustering-Adjusted SAF Target Probability
# Computes 1 - p_eff/12 using Gamma renewal survival function when iat_cv != 1
# Falls back to standard Poisson target when process is ~Poisson or parameters missing
# ============================================================================
.compute_renewal_saf_target <- function(T_years, mu_T_months, shape_k, rate_k, cv_val) {
  # Guard: use standard Poisson if CV is near 1 or parameters are invalid
  if (is.na(cv_val) || is.na(shape_k) || shape_k <= 0 || rate_k <= 0 || 
      (cv_val >= 0.5 && cv_val <= 1.5)) {
    return(1 - mu_T_months / (T_years * 12))
  }
  
  # Gamma renewal survival: P(no event within T years) = 1 - F_Gamma(T; k, rate)
  # Effective annual exceedance probability adjusted for clustering variance
  p_no_event_T_years <- 1 - pgamma(T_years, shape = shape_k, rate = rate_k)
  p_annual_exceed    <- (1 - p_no_event_T_years) / T_years
  
  # Convert to monthly target for SAF inversion
  target_adj <- 1 - (p_annual_exceed / 12)
  return(pmax(1e-9, pmin(1 - 1e-9, target_adj)))
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
#
compute_cfg_pickands <- function(u, v, t_val = 0.5, n_grid = 300) {
  u_grid <- seq(1e-5, 1 - 1e-5, length.out = n_grid)
  c_vals <- vapply(seq_along(u_grid), function(k) {
    mean(u <= u_grid[k]^(1 - t_val) & v <= u_grid[k]^t_val)
  }, numeric(1))
  exp(mean(log(pmax(c_vals, 1e-12))))
}
#
fit_copulas <- function(drought_data, marginal_fits, index_name, out_dir) {
  dc <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  if (nrow(dc) < 10 || is.null(marginal_fits)) return(NULL)
  
  # Rank-based pseudo-observations (Empirical Marginal Transformation)
  n_obs <- nrow(dc)
  
  # Transform severity to uniform [0,1] using ranks
  uS <- rank(dc$severity, ties.method = "average") / (n_obs + 1)
  
  # Transform area_pct to uniform [0,1] using ranks
  uA <- rank(dc$area_pct, ties.method = "average") / (n_obs + 1)
  uM   <- cbind(uS, uA); n  <- nrow(uM)
  
  # Copula candidate set
  cops  <- list(
    Frank       = frankCopula(dim = 2),
    Gumbel      = gumbelCopula(dim = 2),
    Plackett    = plackettCopula(),
    SurvClayton = rotCopula(claytonCopula(dim = 2)),
    Joe         = joeCopula(dim = 2),
    Gaussian    = normalCopula(param = 0.5, dim = 2, dispstr = "un"),
    StudentT    = tCopula(param = 0.5, df = 4, dim = 2)
  )
  
  res   <- data.frame(Copula=character(), Parameter=numeric(),
                      LogLik=numeric(), AIC=numeric(), BIC=numeric(),
                      stringsAsFactors=FALSE)
  fits  <- list()
  for (nm in names(cops)) {
    fit  <- tryCatch(
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
      p   <- coef(fit)
      ll  <- loglikCopula(p, uM, cops[[nm]])
      k   <- length(p)
      res   <- rbind(res, data.frame(Copula    = nm,
                                     Parameter = p[1],
                                     LogLik    = ll,
                                     AIC       = -2*ll + 2*k,
                                     BIC       = -2*ll + log(n)*k))
      fits[[nm]]  <- fit
    }
  }
  if (nrow(res) == 0) return(NULL)
  
  #==========================================================================
  # TAIL-DEPENDENCE PRIORITY LOGIC (as per methodology guide)
  #==========================================================================
  
  # 1. Calculate empirical tail dependence using CFG estimator
  extreme_fams  <- c("Gumbel", "Joe", "SurvClayton")
  
  # CFG estimator for upper tail dependence
  A_05  <- compute_cfg_pickands(uS, uA, t_val = 0.5)
  emp_lambda_U <- min(max(0, 2 * (1 - A_05)),1.00)
  
  # 2. Calculate overall rank correlation (Kendall's tau)
  tau_val <- cor(uS, uA, method = "kendall")
  
  # 3. Calculate Spearman's rho for comparison
  rho_val <- cor(uS, uA, method = "spearman")
  
  # 4. TAIL-DEPENDENCE PRIORITY DECISION
  # Prioritize Joe/Gumbel when extreme dependence > overall correlation
  td_threshold <- 0.10
  prioritize_tail_dep <- (emp_lambda_U > td_threshold) && ((emp_lambda_U) > abs(tau_val))
  
  if (prioritize_tail_dep) {
    log_event(sprintf("[%s] TAIL-DEPENDENCE PRIORITY: λ_U=%.4f > |τ|=%.4f. Enforcing Joe/Gumbel selection.",
                      index_name, emp_lambda_U, abs(tau_val)))
    
    # Filter to tail-dependent families only
    td_candidates <- res[res$Copula %in% c("Joe", "Gumbel", "SurvClayton", "StudentT"), ]
    
    if (nrow(td_candidates) > 0) {
      # Prioritize Joe and Gumbel (Archimedean with upper tail)
      jg_candidates <- td_candidates[td_candidates$Copula %in% c("Joe", "Gumbel"), ]
      
      if (nrow(jg_candidates) > 0) {
        # Select best among Joe/Gumbel by AIC
        best <- jg_candidates$Copula[which.min(jg_candidates$AIC)]
        log_event(sprintf("[%s] Selected '%s' (AIC=%.2f) from Joe/Gumbel family (prioritized for tail dependence)",
                          index_name, best, min(jg_candidates$AIC)))
      } else {
        # Fallback to any tail-dependent family
        best <- td_candidates$Copula[which.min(td_candidates$AIC)]
        log_event(sprintf("[%s] Selected '%s' (AIC=%.2f) from tail-dependent families",
                          index_name, best, min(td_candidates$AIC)))
      }
    } else {
      # No tail-dependent candidates - use standard AIC
      best <- res$Copula[which.min(res$AIC)]
    }
  } else {
    # Standard AIC selection
    best <- res$Copula[which.min(res$AIC)]
    log_event(sprintf("[%s] Standard AIC selection: '%s' (λ_U=%.4f ≤ |τ|=%.4f)",
                      index_name, best, emp_lambda_U, abs(tau_val)))
  }
  
  #==========================================================================
  # THEORETICAL TAIL DEPENDENCE FOR BEST COPULA
  #==========================================================================
  best_par   <- coef(fits[[best]])[1]
  theo_td    <- tryCatch({
    switch(best,
           Gumbel      = list(lambda_U = 2 - 2^(1/max(best_par, 1.001))),
           Joe         = list(lambda_U = 2 - 2^(1/max(best_par, 1.001))),
           SurvClayton = list(lambda_U = 2^(-1/max(best_par, 1e-6))),
           Frank       = list(lambda_U = 0),
           Plackett    = list(lambda_U = 0),
           Gaussian    = list(lambda_U = 0),
           StudentT    = { 
             rho <- max(min(best_par, 0.9999), -0.9999)
             df_par <- max(coef(fits[[best]])[2], 1.5)
             list(lambda_U = 2 * pt(-sqrt((df_par + 1)*(1-rho)/(1+rho)), df=df_par+1)) 
           },
           list(lambda_U = NA_real_))
  }, error = function(e) list(lambda_U = NA_real_))
  
  #==========================================================================
  # VALIDATION & RECOMMENDATION
  #==========================================================================
  td_recommendation <- if (!is.finite(emp_lambda_U)) "Tail dependence inconclusive"
  else if (emp_lambda_U > 0.10) "Upper-tail dependence detected: Joe/Gumbel prioritized (CFG validated)"
  else "No meaningful upper-tail dependence: Frank/Plackett adequate"
  
  zero_tail_families <- c("Frank", "Plackett", "Gaussian")
  td_mismatch <- isTRUE(emp_lambda_U > 0.10) && best %in% zero_tail_families && !is.na(emp_lambda_U)
  
  log_event(sprintf("[%s] CFG Validation — empirical λ_U=%.4f | theoretical (best=%s) λ_U=%.4f | τ=%.4f",
                    index_name, emp_lambda_U, best, 
                    if(is.finite(theo_td$lambda_U)) theo_td$lambda_U else NA,
                    tau_val))
  log_event(sprintf("[%s] TD recommendation: %s", index_name, td_recommendation))
  
  if (td_mismatch) {
    log_event(sprintf("[%s] WARNING [TD-U]: Selected copula (%s) has λ_U=0 but empirical λ_U=%.4f>0.10",
                      index_name, best, emp_lambda_U))
  }
  
  #==========================================================================
  # APPEND TO RESULTS TABLE
  #==========================================================================
  res$emp_lambda_U      <- round(emp_lambda_U, 4)
  res$kendall_tau       <- round(tau_val, 4)
  res$spearman_rho      <- round(rho_val, 4)
  
  theo_U <- vapply(names(cops), function(nm) {
    if (is.null(fits[[nm]])) return(NA_real_)
    p <- coef(fits[[nm]])
    switch(nm,
           Gumbel=2-2^(1/max(p[1],1.001)), 
           Joe=2-2^(1/max(p[1],1.001)),
           SurvClayton=2^(-1/max(p[1],1e-6)), 
           Frank=0, Plackett=0, Gaussian=0,
           StudentT={ 
             rho<-max(min(p[1],0.9999),-0.9999)
             df_p<-max(if(length(p)>=2)p[2] else 4,1.5)
             2*pt(-sqrt((df_p+1)*(1-rho)/(1+rho)), df=df_p+1) 
           },
           NA_real_)
  }, numeric(1L))
  
  res$theo_lambda_U     <- round(theo_U, 4)
  res$td_mismatch       <- td_mismatch & (res$Copula == best)
  res$td_recommendation <- td_recommendation
  res$tail_dep_priority <- prioritize_tail_dep
  
  write.csv(res, file.path(out_dir, sprintf("%s_copula_comparison.csv", index_name)), row.names = FALSE)
  
  list(best_copula_name=best, 
       best_copula_fit=fits[[best]], 
       all_results=res,
       u_severity=uS, 
       u_area=uA, 
       emp_lambda_U=emp_lambda_U,
       theo_lambda_U_best=theo_td$lambda_U, 
       kendall_tau=tau_val,
       spearman_rho=rho_val,
       td_mismatch=td_mismatch,
       td_recommendation=td_recommendation, 
       tail_dep_priority=prioritize_tail_dep,
       cfg_estimator=A_05)
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
                                          index_name, out_dir, duration_class = "all",
                                          record_n_months = NULL){
  if (is.null(marginal_fits) || is.null(copula_fit)) return(NULL)
  log_event(sprintf("[%s | %s] SAF — Conditional Copula...", index_name, duration_class))
  dc    <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_tot <- nrow(drought_data); n_dr <- nrow(dc); if (n_dr == 0) return(NULL)
  mu_T  <- n_tot / n_dr
  # Compute renewal parameters from event inter-arrival times
  indr_all    <- drought_data$severity > 0 & drought_data$area_pct > 0
  is_start    <- indr_all & !c(FALSE, head(indr_all, -1))
  ev_starts   <- as.Date(drought_data$date[is_start])
  iat_cv_val  <- NA_real_
  shape_k_val <- NA_real_
  rate_k_val  <- NA_real_
  mu_T_adj    <- mu_T * 12          # default: Poisson (months)
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
  cop_obj  <- copula_fit$best_copula_fit@copula
  beta_par <- marginal_fits$area_fit$estimate
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5); area_prop <- area_pct / 100
  out <- data.frame()
  
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
    if (isTRUE(copula_fit$td_mismatch) && T >= 50) {
      log_event(sprintf("  [%s] T=%d: Skipping Kendall (zero-tail copula). Use Conditional method.", index_name, T))
      next
    }
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
                                      index_name, out_dir, duration_class = "all",
                                      record_n_months = NULL) {
  if (is.null(marginal_fits) || is.null(copula_fit)) return(NULL)
  
  log_event(sprintf("[%s | %s] SAF — Kendall Distribution (analytical K_C)...",
                    index_name, duration_class))
  
  # Explicit nomenclature: record length vs class-specific drought count
  n_record     <- if (is.null(record_n_months)) nrow(drought_data) else record_n_months
  dc           <- drought_data[drought_data$severity > 0 & drought_data$area_pct > 0, ]
  n_class_dr   <- nrow(dc)
  if (n_class_dr == 0) return(NULL)
  mu_T         <- n_record / n_class_dr  # Conditional mean inter-arrival (months)
  # Compute renewal parameters from event inter-arrival times
  indr_all    <- drought_data$severity > 0 & drought_data$area_pct > 0
  is_start    <- indr_all & !c(FALSE, head(indr_all, -1))
  ev_starts   <- as.Date(drought_data$date[is_start])
  iat_cv_val  <- NA_real_
  shape_k_val <- NA_real_
  rate_k_val  <- NA_real_
  mu_T_adj    <- mu_T * 12          # default: Poisson (months)
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
  cop_obj  <- copula_fit$best_copula_fit@copula
  cop_name <- copula_fit$best_copula_name
  beta_par <- marginal_fits$area_fit$estimate
  uS       <- copula_fit$u_severity
  uA       <- copula_fit$u_area
  
  # USE ANALYTICAL K_C FOR ARCHIMEDEAN COPULAS
  use_analytical <- cop_name %in% c("Frank", "Gumbel", "Joe", "SurvClayton")
  
  if (use_analytical) {
    # Analytical K_C(t) = t - φ(t)/φ'(t) for Archimedean families
    analytical_K_C <- function(t, cop_obj, fam) {
      t <- pmax(1e-10, pmin(t, 1 - 1e-10))
      theta <- cop_obj@parameters[1]
      
      if (fam == "Frank") {
        phi  <- -log((exp(-theta * t) - 1) / (exp(-theta) - 1))
        phi_prime <- theta * exp(-theta * t) / (1 - exp(-theta * t))
      } else if (fam == "Gumbel") {
        log_t <- log(t)
        phi  <- (-log_t)^theta
        phi_prime <- -theta * (-log_t)^(theta - 1) / t
      } else if (fam == "Joe") {
        one_m_t <- 1 - t
        phi  <- -log(1 - one_m_t^theta)
        phi_prime <- theta * one_m_t^(theta - 1) / (1 - one_m_t^theta)
      } else { # SurvClayton
        phi  <- (t^(-theta) - 1) / theta
        phi_prime <- -t^(-theta - 1)
      }
      pmax(1e-10, pmin(1 - 1e-10, t - phi / phi_prime))
    }
    kc_fn <- function(t) analytical_K_C(t, cop_obj, cop_name)
    opt_int <- c(1e-8, 1 - 1e-8)
    log_event(sprintf("  [%s] Using analytical K_C for %s (valid extrapolation)", index_name, cop_name))
  } else {
    cop_vals <- copula::pCopula(cbind(uS, uA), cop_obj)
    kc_fn <- function(t) mean(cop_vals <= t)
    opt_int <- c(1e-4, 1 - 1e-4)
    log_event(sprintf("  [%s] Using empirical K_C for %s", index_name, cop_name))
  }
  
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5); area_prop <- area_pct / 100
  
  out <- data.frame()
  T_years  <- c(10, 25, 50, 100)
  area_pct <- seq(5, 95, by = 5); area_prop <- area_pct / 100
  out <- data.frame()
  
  for (T in T_years) {
    # Pass mu_T_adj_months, renewal_shape_k, renewal_rate_k, and iat_cv from rp table
    # or compute them inline from the drought data if passed as arguments.
    # Example inline usage (adapt variable names to your scope):
    target_kc <- .compute_renewal_saf_target(
      T_years      = T,
      mu_T_months  = mu_T_adj,    # Replace with iat_mean * 12
      shape_k      = shape_k_val, # Replace with 1 / (iat_cv^2)
      rate_k       = rate_k_val,  # Replace with shape_k / iat_mean
      cv_val       = iat_cv_val   ) # Replace with computed iat_cv
    if (target_kc <= 0 || target_kc >= 1) next
    
    # Solve K_C(t_val) = target_kc
    kc_obj <- function(t) (kc_fn(t) - target_kc)^2
    
    # WIDER search interval for analytical K_C (enables extrapolation)
    if (use_analytical) {
      t_sol <- tryCatch(optimize(kc_obj, interval = c(1e-8, 1 - 1e-8)),
                        error = function(e) list(minimum = 0.5))
    } else {
      t_sol <- tryCatch(optimize(kc_obj, interval = c(1e-4, 1 - 1e-4)),
                        error = function(e) list(minimum = 0.5))
    }
    
    t_val <- t_sol$minimum
    
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

subset_drought_chars_by_class <- function(drought_chars, events_df, class_label, record_n_months) {
  ev <- events_df[events_df$duration_class == class_label, ]
  if (nrow(ev) == 0) return(NULL)
  
  in_class <- logical(nrow(drought_chars))
  for (i in seq_len(nrow(ev)))
    in_class <- in_class | (drought_chars$date >= ev$start_date[i] & drought_chars$date <= ev$end_date[i])
  
  dc_sub <- drought_chars[in_class, , drop = FALSE]
  attr(dc_sub, "record_n_months") <- record_n_months  # Explicit metadata
  dc_sub
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
    dc_sub <- subset_drought_chars_by_class(dc, events, cl, n_months)
    if (is.null(dc_sub)) next
    lbl    <- paste0(index_name, "_", cl)
    m_sub  <- fit_marginal_distributions(dc_sub, lbl, out_dir)
    c_sub  <- if (!is.null(m_sub)) fit_copulas(dc_sub, m_sub, lbl, out_dir) else NULL
    
    rec_len <- attr(dc_sub, "record_n_months")
    if (!is.null(m_sub) && !is.null(c_sub)) {
      saf_cond_cls[[cl]] <- derive_SAF_curves_conditional(dc_sub, m_sub, c_sub, index_name, out_dir, cl, rec_len)
      saf_kend_cls[[cl]] <- derive_SAF_curves_kendall(dc_sub, m_sub, c_sub, index_name, out_dir, cl, rec_len)
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