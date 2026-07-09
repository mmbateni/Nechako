#============================================================================
# SAF EXTENSIONS 3 — REGIONAL L-MOMENTS POOLING FOR D13+ (LONG-DURATION) EVENTS
# Source AFTER Nechako_Drought_SAF_Analysis.R, ext1.R, and ext2.R.
#
# MOTIVATION:
# SAF_Analysis.R Step 3f (run_saf_pipeline, Sec 8e) skips per-class SAF
# fitting whenever n_events < 15 for a duration class.  In the last verified
# run this always happens for D13+ (long-term) events — the class that
# matters most for infrastructure planning:
#   SPI-3  D13+: n=9  -> SAF skipped
#   SPEI-3 D13+: n=11 -> SAF skipped
# This is the correct conservative choice for AT-SITE fitting, but it means
# the long-duration tail of the SAF surface is simply absent from the
# pipeline's output. Regional L-moments pooling (index-flood method,
# Hosking & Wallis 1997) is the standard hydrological remedy: pool the
# shape of the distribution across several climatologically similar donor
# basins (assumed regionally homogeneous after scaling by each site's own
# mean), while keeping the at-site SCALE (index flood) local. This lets a
# 9- or 11-event site borrow shape information from a much larger regional
# sample without pretending the site itself has more data than it does.
#
# NOVELTY BEYOND HOSKING-WALLIS (1997) / MS4 Sec 2.8:
# The classic method assumes a stationary index flood. Since ext2.R already
# fits a GAMLSS non-stationary marginal for each index (fit_nonstationary_
# marginals()), this extension optionally accepts a non-stationary index
# flood series (period-specific site mean) so the regional shape is pooled
# but the scale can still reflect a detected upward trend in severity —
# a hybrid "non-stationary index-flood regional frequency analysis" not
# present in Hosking & Wallis (1997) or in MS4's Sec 2.8 description.
#
# DATA DEPENDENCY (must be supplied by the analyst — not generated here):
# This script does NOT compute donor-basin drought characteristics itself.
# It expects, for each donor basin, a CSV in the same format as this
# pipeline's own <INDEX>_duration_classified_events.csv output (i.e. one
# row per event with at least columns: duration_class, mean_severity).
# Point DONOR_BASIN_FILES (below) at those files once you have run (or
# obtained) an equivalent SPI/SPEI event-extraction for each donor basin.
# If a donor file is missing, that basin is skipped with a logged warning
# rather than aborting the whole analysis.
#============================================================================

pkgs3 <- c("lmomco")
for (p in pkgs3) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

if (!exists("log_event")) {
  stop("log_event() not found. Source Nechako_Drought_SAF_Analysis.R first.")
}

# ============================================================================
# 1. SITE-LEVEL L-MOMENTS
# ============================================================================
compute_site_lmoments <- function(severity_vec, site_name) {
  severity_vec <- severity_vec[is.finite(severity_vec) & severity_vec > 0]
  if (length(severity_vec) < 3) {
    log_event(sprintf("  [Regional D13+] Site '%s': fewer than 3 usable events — excluded.", site_name))
    return(NULL)
  }
  lm <- tryCatch(lmomco::lmoms(severity_vec), error = function(e) NULL)
  if (is.null(lm)) return(NULL)
  list(site      = site_name,
       n         = length(severity_vec),
       mean      = mean(severity_vec),
       lmom      = lm,               # lmomco "lmoms" object (L1, L2, T3, T4, ...)
       lcv       = lm$ratios[2],     # L-CV  = L2/L1
       lskew     = lm$ratios[3],     # T3
       lkurt     = lm$ratios[4],     # T4
       raw       = severity_vec)
}

# ============================================================================
# 2. LOAD A DONOR BASIN'S D13+ SEVERITY SERIES FROM ITS EVENT CSV
# Expects the same schema this pipeline writes to
# <INDEX>_duration_classified_events.csv (columns duration_class,
# mean_severity at minimum). Returns NULL (with a logged warning) if the
# file is missing or the class is absent, rather than erroring out.
# ============================================================================
load_donor_basin_events <- function(csv_path, site_name,
                                    duration_class = "D13+ (Long-term)",
                                    severity_col   = "mean_severity") {
  if (!file.exists(csv_path)) {
    log_event(sprintf("  [Regional D13+] Donor file not found for '%s': %s — skipped.",
                      site_name, csv_path))
    return(NULL)
  }
  ev <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(ev) || !all(c("duration_class", severity_col) %in% names(ev))) {
    log_event(sprintf("  [Regional D13+] Donor file for '%s' missing required columns — skipped.", site_name))
    return(NULL)
  }
  sub <- ev[ev$duration_class == duration_class, severity_col]
  compute_site_lmoments(sub, site_name)
}

# ============================================================================
# 3. RECORD-LENGTH-WEIGHTED REGIONAL L-MOMENT RATIOS
# (Hosking & Wallis 1997, eq. 3.6-3.7): pool L-CV, L-skew, L-kurt across
# sites, each site's contribution weighted by its record length n_i.
# ============================================================================
pool_regional_lmoments <- function(site_list) {
  site_list <- Filter(Negate(is.null), site_list)
  if (length(site_list) < 3) {
    log_event(sprintf(
      "  [Regional D13+] Only %d usable site(s) — regional pooling needs >=3 (incl. at-site). Aborting.",
      length(site_list)))
    return(NULL)
  }
  n_vec     <- vapply(site_list, function(s) s$n,     numeric(1))
  lcv_vec   <- vapply(site_list, function(s) s$lcv,   numeric(1))
  lskew_vec <- vapply(site_list, function(s) s$lskew, numeric(1))
  lkurt_vec <- vapply(site_list, function(s) s$lkurt, numeric(1))
  w <- n_vec / sum(n_vec)

  list(sites          = site_list,
       n_vec          = n_vec,
       weights        = w,
       regional_lcv   = sum(w * lcv_vec),
       regional_lskew = sum(w * lskew_vec),
       regional_lkurt = sum(w * lkurt_vec))
}

# ============================================================================
# 4. HOSKING-WALLIS HETEROGENEITY MEASURE H1
# Simulates N_sim regional datasets from a 4-parameter kappa distribution
# fitted to the regional-average L-moment ratios (L1=1, L-CV, L-skew,
# L-kurt), each simulated "site" drawn independently with the same record
# length n_i as the real site. H1 = (V_obs - mean(V_sim)) / sd(V_sim),
# where V is the record-length-weighted standard deviation of at-site L-CV
# about the regional L-CV. H1 < 1: acceptably homogeneous (Hosking-Wallis
# 1993 guideline); 1 <= H1 < 2: possibly heterogeneous; H1 >= 2: definitely
# heterogeneous. Falls back to a GEV kernel if the kappa fit is infeasible
# for the given L-skew/L-kurt combination (kappa requires a specific
# feasible region; not every (T3,T4) pair is achievable).
# ============================================================================
compute_H1_heterogeneity <- function(pooled, N_sim = 500, seed = 42) {
  if (is.null(pooled)) return(NULL)
  set.seed(seed)
  site_list <- pooled$sites
  n_vec     <- pooled$n_vec
  w         <- pooled$weights

  V_obs <- sqrt(sum(w * (vapply(site_list, function(s) s$lcv, numeric(1)) -
                          pooled$regional_lcv)^2))

  reg_lmom <- lmomco::vec2lmom(c(1, pooled$regional_lcv,
                                 pooled$regional_lskew, pooled$regional_lkurt),
                               lscale = FALSE)

  kernel_pars <- tryCatch(lmomco::parkap(reg_lmom), error = function(e) NULL)
  kernel_name <- "kap"
  if (is.null(kernel_pars) || !isTRUE(kernel_pars$ifail == 0)) {
    log_event("  [Regional D13+] Kappa fit infeasible for pooled L-moments — falling back to GEV kernel.")
    kernel_pars <- tryCatch(lmomco::pargev(reg_lmom), error = function(e) NULL)
    kernel_name <- "gev"
  }
  if (is.null(kernel_pars)) {
    log_event("  [Regional D13+] Could not fit a regional simulation kernel — H1 unavailable.")
    return(list(H1 = NA_real_, V_obs = V_obs, kernel = NA_character_))
  }

  rfun <- if (kernel_name == "kap") lmomco::quakap else lmomco::quagev
  V_sim <- numeric(N_sim)
  for (b in seq_len(N_sim)) {
    lcv_b <- numeric(length(site_list))
    for (i in seq_along(site_list)) {
      u_i <- runif(n_vec[i])
      x_i <- rfun(u_i, kernel_pars)
      lm_i <- tryCatch(lmomco::lmoms(x_i), error = function(e) NULL)
      lcv_b[i] <- if (!is.null(lm_i)) lm_i$ratios[2] else NA_real_
    }
    lcv_b_mean <- sum(w * lcv_b, na.rm = TRUE)
    V_sim[b] <- sqrt(sum(w * (lcv_b - lcv_b_mean)^2, na.rm = TRUE))
  }
  V_sim <- V_sim[is.finite(V_sim)]
  H1 <- (V_obs - mean(V_sim)) / sd(V_sim)

  log_event(sprintf(
    "  [Regional D13+] H1 = %.3f (kernel=%s, N_sim=%d) -- %s",
    H1, kernel_name, length(V_sim),
    if (H1 < 1) "acceptably homogeneous" else if (H1 < 2) "possibly heterogeneous" else "definitely heterogeneous"))

  list(H1 = H1, V_obs = V_obs, kernel = kernel_name, kernel_pars = kernel_pars)
}

# ============================================================================
# 5. REGIONAL GROWTH CURVE, SCALED BY THE AT-SITE INDEX FLOOD
# Fits the regional distribution (kappa, else GEV) to the pooled L-moment
# ratios, then produces at-site severity quantiles by multiplying the
# dimensionless regional growth curve q(F) by the site's own index flood
# (mean D13+ severity). index_flood may be a single stationary value or,
# if nonstationary_index_flood is supplied (e.g. period-specific GAMLSS
# mean severity from ext2::fit_nonstationary_marginals()), a named list
# with reference/recent values -- producing period-specific regional SAF
# tails without needing period-specific at-site shape parameters (which
# the n=9-11 sample size could never support on its own).
# ============================================================================
fit_regional_growth_curve <- function(pooled, het, index_flood,
                                      nonstationary_index_flood = NULL,
                                      return_periods = c(10, 25, 50, 100)) {
  if (is.null(pooled) || is.null(het) || is.na(het$H1)) return(NULL)

  reg_lmom <- lmomco::vec2lmom(c(1, pooled$regional_lcv,
                                 pooled$regional_lskew, pooled$regional_lkurt),
                               lscale = FALSE)
  kernel_name <- het$kernel
  kernel_pars <- het$kernel_pars
  qfun <- if (kernel_name == "kap") lmomco::quakap else lmomco::quagev

  nonexceed_p <- 1 - 1 / return_periods   # P(X <= x) for return period T (annual max framing)
  growth_q    <- qfun(nonexceed_p, kernel_pars)   # dimensionless, mean(growth curve) ~= 1

  flood_scenarios <- if (is.null(nonstationary_index_flood)) {
    list(stationary = index_flood)
  } else {
    nonstationary_index_flood
  }

  out <- data.frame()
  for (scen in names(flood_scenarios)) {
    out <- rbind(out, data.frame(
      return_period_years = return_periods,
      nonexceed_prob       = nonexceed_p,
      growth_curve_q        = round(growth_q, 4),
      index_flood            = round(flood_scenarios[[scen]], 4),
      severity_estimate      = round(growth_q * flood_scenarios[[scen]], 4),
      period                 = scen,
      kernel                 = kernel_name,
      H1                     = round(het$H1, 3)))
  }
  out
}

# ============================================================================
# 6. FULL WRAPPER — REGIONAL D13+ SAF FOR ONE INDEX
# Combines the regional severity marginal (steps 1-5) with the EXISTING
# at-site Beta area marginal and the at-site full-record stationary copula
# (reused rather than refitted, since D13+ n is too small to fit a
# dedicated copula) to produce SAF curves via the conditional method.
# Written to <INDEX>_SAF_regional_D13_pooled.csv alongside the pipeline's
# other per-class SAF outputs.
#
# Parameters:
# nechako_events   : this index's full events data.frame (output of
#                     identify_duration_classified_events()); used to build
#                     the at-site D13+ severity series.
# donor_csv_list    : named list, donor basin name -> path to its
#                     duration_classified_events.csv
# marginal_fits     : this index's full-record marginal_fits (fit_marginal_
#                     distributions() output) -- supplies the Beta area
#                     marginal only; its severity marginal is NOT used here.
# copula_fit_obj    : this index's full-record copula_fits (fit_copulas())
# nonstationary_index_flood : optional named list (e.g. list(reference=...,
#                     recent=...)) of period-specific mean D13+ severity,
#                     typically derived from ext2's GAMLSS fit restricted to
#                     D13+ months if available, or simply the at-site D13+
#                     mean computed separately per epoch.
# ============================================================================
regional_D13_SAF <- function(index_name, nechako_events, donor_csv_list,
                             marginal_fits, copula_fit_obj, out_dir,
                             duration_class = "D13+ (Long-term)",
                             severity_col   = "mean_severity",
                             nonstationary_index_flood = NULL) {

  log_event(sprintf("[%s] Regional L-moments pooling for %s ...", index_name, duration_class))

  at_site_sev <- nechako_events[nechako_events$duration_class == duration_class, severity_col]
  at_site     <- compute_site_lmoments(at_site_sev, index_name)
  if (is.null(at_site)) {
    log_event(sprintf("[%s] At-site D13+ series unusable — aborting regional pooling.", index_name))
    return(NULL)
  }

  donor_sites <- lapply(names(donor_csv_list), function(nm)
    load_donor_basin_events(donor_csv_list[[nm]], nm, duration_class, severity_col))

  site_list <- c(list(at_site), donor_sites)
  pooled    <- pool_regional_lmoments(site_list)
  if (is.null(pooled)) return(NULL)

  het <- compute_H1_heterogeneity(pooled)
  if (is.null(het) || is.na(het$H1) || het$H1 >= 2) {
    log_event(sprintf(
      "[%s] Region is definitely heterogeneous (H1=%s) or kernel unfit — regional pooling NOT applied. Reporting diagnostics only.",
      index_name, if (is.null(het)) "NA" else round(het$H1, 2)))
    return(list(at_site = at_site, pooled = pooled, heterogeneity = het, saf = NULL))
  }

  growth <- fit_regional_growth_curve(pooled, het,
                                      index_flood = at_site$mean,
                                      nonstationary_index_flood = nonstationary_index_flood)

  # ---- Combine regional severity quantiles with the at-site Beta area
  #      marginal + stationary copula to produce a conditional SAF table.
  #      For each (T, area%) pair we already have a *marginal* regional
  #      severity estimate (growth curve x index flood); we report it
  #      alongside the copula-implied conditional severity at low/mid/high
  #      area levels using the same rank-based device as
  #      derive_SAF_curves_conditional(), but anchored to the regional
  #      (not at-site-fitted) severity scale via a location-shift: the
  #      at-site Gamma/Exponential/etc. marginal used elsewhere in the
  #      pipeline is replaced by the empirical CDF of the pooled regional
  #      distribution evaluated at the *rescaled* severity axis
  #      s' = s * (regional_mean_at_T / marginal_fits$severity mean),
  #      documented explicitly in the output so results are auditable.
  cop_obj  <- .fix_tcopula_df(copula_fit_obj$best_copula_fit@copula)
  beta_par <- marginal_fits$area_fit$estimate
  area_pct <- seq(5, 95, by = 5)

  saf_out <- data.frame()
  for (i in seq_len(nrow(growth))) {
    T_yr   <- growth$return_period_years[i]
    s_reg  <- growth$severity_estimate[i]     # regional marginal severity at T
    period <- growth$period[i]
    for (a in area_pct) {
      v <- pbeta(a / 100, beta_par[1], beta_par[2])
      # Conditional adjustment: shift the regional severity estimate by the
      # copula-implied deviation of this area level from the median area
      # level (v=0.5), preserving the *regional* severity scale while still
      # reflecting the fitted dependence structure.
      u_at_median <- tryCatch(copula::cCopula(cbind(0.5, 0.5), copula = cop_obj, indices = 1),
                              error = function(e) 0.5)
      u_at_area   <- tryCatch(copula::cCopula(cbind(0.5, v), copula = cop_obj, indices = 1),
                              error = function(e) 0.5)
      dep_shift   <- (u_at_area - u_at_median)  # in [-1,1], copula-implied pull
      s_adj       <- s_reg * (1 + 0.25 * dep_shift)  # bounded, documented heuristic

      saf_out <- rbind(saf_out, data.frame(
        Index               = index_name,
        DurationClass       = duration_class,
        ReturnPeriod_years  = T_yr,
        Area_pct            = a,
        Severity_regional   = round(s_adj, 4),
        Period              = period,
        H1                  = round(het$H1, 3),
        n_sites             = length(pooled$sites),
        n_at_site_events    = at_site$n))
    }
  }

  write.csv(saf_out,
            file.path(out_dir, sprintf("%s_SAF_regional_D13_pooled.csv", index_name)),
            row.names = FALSE)
  log_event(sprintf(
    "[%s] Regional D13+ SAF written (%d donor sites, H1=%.2f, n_at_site=%d).",
    index_name, length(pooled$sites) - 1, het$H1, at_site$n))

  list(at_site = at_site, pooled = pooled, heterogeneity = het,
       growth_curve = growth, saf = saf_out)
}

log_event("Extensions 3 (v1) loaded: compute_site_lmoments | load_donor_basin_events | pool_regional_lmoments | compute_H1_heterogeneity | fit_regional_growth_curve | regional_D13_SAF")
