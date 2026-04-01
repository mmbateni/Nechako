####################################################################################
# w12_bayesian_ene_analysis.R
#
# MANUSCRIPT 2 ANALYSIS:
#   Non-stationary recurrence risk of the 2022–2025 Nechako mega-drought:
#   Bayesian Expected Number of Exceedances with posterior credible intervals
#
# DESCRIPTION:
#   Fits a Bayesian non-stationary Gamma GLM to the SPEI-3 drought severity
#   catalog, using the JJA thermodynamic fraction of drought severity (F_thm)
#   as a physically motivated covariate.  Weakly informative priors on the
#   tail shape prevent overfitting to the record-shattering 2022–2025 event.
#   Posterior predictive sampling yields the full probability distribution of:
#     (a) T_eff(t) — effective return period at any calendar year t
#     (b) ENE(t0, m) — expected number of exceedances over a future horizon
#
# DEPENDENCIES (must run BEFORE this script):
#   w5_event_ranking.R  + w5_EXPORT_ADDON.R
#     → temporal_drought/event_ranking/ms2_annual_severity_series.csv
#     → temporal_drought/event_ranking/ms2_copula_fit.rds
#   w11_dynamic_thermodynamic_decomp.R  + w11_EXPORT_ADDON.R
#     → decomp_results/thm_frac_annual_jja_ms2.csv
#
# OUTPUTS  (bayesian_ene/):
#   brms_model_spei3.rds          fitted brms model object
#   posterior_draws_spei3.csv     tidy posterior draws (alpha, beta, phi)
#   T_eff_trajectory.csv          T_eff(t) posterior summary 1950-2025
#   ENE_posterior.csv             ENE distribution over future horizon
#   four_timepoints.csv           T_eff at 1987, 2023, 2025 + 30-yr ENE
#   sensitivity_priors.csv        prior sensitivity table
#   Fig_T_eff_trajectory.pdf/.png
#   Fig_ENE_posterior.pdf/.png
#   Fig_PPC.pdf/.png              posterior predictive checks
#   Fig_prior_sensitivity.pdf/.png
#
# EXECUTION:
#   setwd("D:/Nechako_Drought/Nechako/")
#   source("w12_bayesian_ene_analysis.R")
#
# REFERENCES:
#   Bürkner (2017) brms: An R Package for Bayesian Multilevel Models Using Stan.
#     J. Stat. Software, 80(1).
#   Fischer et al. (2021) Increasing probability of record-shattering climate
#     extremes. Nature Climate Change, 11, 689–695.
#   Rootzén & Katz (2013) Design life level. Water Resour. Res., 49, 5964–5972.
#   Salas & Obeysekera (2014) J. Hydrol. Eng., 19(3), 554–568.
#   Zeder et al. (2023) Geophys. Res. Lett., 50, e2023GL104090.
####################################################################################

# ── 0. PACKAGE INSTALLATION (run once if needed) ─────────────────────────────
# install.packages(c("brms","posterior","bayesplot","loo","tidybayes",
#                    "dplyr","ggplot2","patchwork","scales","tidyr","readr"))
# Note: brms requires a working C++ toolchain and rstan.
# On Windows: install Rtools from https://cran.r-project.org/bin/windows/Rtools/

suppressPackageStartupMessages({
  library(brms)        # Bayesian GLM via Stan
  library(posterior)   # tidy posterior draws
  library(bayesplot)   # posterior predictive check plots
  library(loo)         # leave-one-out cross-validation
  library(tidybayes)   # tidyverse posterior interface
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)
})

source("DROUGHT_ANALYSIS_utils.R")   # WD_PATH, RANKING_DIR, etc.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x

####################################################################################
# SECTION 1: CONFIGURATION
####################################################################################

cat("============================================================\n")
cat("  w12: Bayesian ENE Analysis — SECTION 1: Configuration\n")
cat("============================================================\n")

WD_PATH     <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")
RANKING_DIR <- file.path(WD_PATH, "temporal_drought", "event_ranking")
DECOMP_DIR  <- file.path(WD_PATH, "decomp_results")
ENE_DIR     <- file.path(WD_PATH, "bayesian_ene")
dir.create(ENE_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Analysis parameters ────────────────────────────────────────────────────────
TARGET_IDX      <- "spei_03"
TARGET_SCALE    <- 3L
HISTORICAL_START <- 1950L
HISTORICAL_END   <- 2025L

# The four time-points for T_eff evaluation (see manuscript Section 3.4)
T_BASELINE  <- 1987L   # WMO 1961-1990 midpoint: pre-warming counterfactual
T_OCCURRED  <- 2023L   # mid-year of the 2022-2025 event: as-occurred rarity
T_CURRENT   <- 2025L   # end of study: current effective risk
ENE_START   <- 2026L   # start of forward planning window
ENE_END     <- 2055L   # end of forward planning window (30 years)

# ── Stan/brms settings ────────────────────────────────────────────────────────
N_CHAINS    <- 4L
N_ITER      <- 4000L   # 2000 warmup + 2000 sampling per chain → 8000 draws total
N_WARMUP    <- 2000L
ADAPT_DELTA <- 0.95    # higher = fewer divergences for heavy-tailed posteriors
N_CORES     <- min(parallel::detectCores() - 1L, N_CHAINS)
SEED        <- 20250322L

cat(sprintf("  WD_PATH     : %s\n", WD_PATH))
cat(sprintf("  ENE_DIR     : %s\n", ENE_DIR))
cat(sprintf("  Target index: %s  |  Scale: %d\n", TARGET_IDX, TARGET_SCALE))
cat(sprintf("  Stan chains : %d  |  Iterations: %d  |  Cores: %d\n",
            N_CHAINS, N_ITER, N_CORES))

####################################################################################
# SECTION 2: DATA LOADING AND PREPARATION
####################################################################################

cat("\n============================================================\n")
cat("  SECTION 2: Data loading\n")
cat("============================================================\n")

# ── 2a. Event catalog — annual severity series ────────────────────────────────
sev_path <- file.path(RANKING_DIR, "ms2_annual_severity_series.csv")
if (!file.exists(sev_path))
  stop("ms2_annual_severity_series.csv not found.\n",
       "Run w5_event_ranking.R + w5_EXPORT_ADDON.R first.")

ann_sev <- read.csv(sev_path, stringsAsFactors = FALSE) %>%
  dplyr::filter(!is.na(max_severity), drought_year == TRUE) %>%
  dplyr::select(year, max_severity, max_duration, n_events)

cat(sprintf("  Severity series: %d drought years (%d–%d)\n",
            nrow(ann_sev), min(ann_sev$year), max(ann_sev$year)))
cat(sprintf("  Severity range: %.3f – %.3f  (median: %.3f)\n",
            min(ann_sev$max_severity), max(ann_sev$max_severity),
            median(ann_sev$max_severity)))

# ── 2b. Copula reference (stationary T_K from MS1) ────────────────────────────
rds_path <- file.path(RANKING_DIR, "ms2_copula_fit.rds")
cop_ref  <- if (file.exists(rds_path)) readRDS(rds_path) else NULL
if (!is.null(cop_ref)) {
  TARGET_SEVERITY <- cop_ref$target_severity
  TARGET_DURATION <- cop_ref$target_duration
  T_K_STATIONARY  <- cop_ref$T_K_stationary
  CI_LO_STAT      <- cop_ref$CI_lo_stationary
  CI_HI_STAT      <- cop_ref$CI_hi_stationary
  LAMBDA_EVENTS   <- 1 / cop_ref$mu_ia_years  # events per year (Poisson rate)
  cat(sprintf("  Stationary T_K: %.0f yr (95%% CI: %.0f–%.0f)\n",
              T_K_STATIONARY, CI_LO_STAT, CI_HI_STAT))
  cat(sprintf("  Target event:   S = %.3f, D = %d mo\n",
              TARGET_SEVERITY, TARGET_DURATION))
  cat(sprintf("  Event rate λ:   %.4f events/year\n", LAMBDA_EVENTS))
} else {
  warning("ms2_copula_fit.rds not found. Using fallback values from manuscript.")
  TARGET_SEVERITY <- 10.728
  TARGET_DURATION <- 8L
  T_K_STATIONARY  <- 921
  CI_LO_STAT      <- 109
  CI_HI_STAT      <- 55386
  LAMBDA_EVENTS   <- 49 / 76   # approx from known catalog size
}

# ── 2c. Thermodynamic fraction covariate ──────────────────────────────────────
cov_path <- file.path(DECOMP_DIR, "thm_frac_annual_jja_ms2.csv")
if (!file.exists(cov_path))
  stop("thm_frac_annual_jja_ms2.csv not found.\n",
       "Run w11_dynamic_thermodynamic_decomp.R + w11_EXPORT_ADDON.R first.")

thm_cov <- read.csv(cov_path, stringsAsFactors = FALSE) %>%
  dplyr::filter(scale == TARGET_SCALE) %>%
  dplyr::select(year, thm_frac_jja, thm_frac_jja_std)

cat(sprintf("  thm_frac covariate: %d annual values (SPEI-%d JJA)\n",
            nrow(thm_cov), TARGET_SCALE))
cat(sprintf("  thm_frac range: %.3f – %.3f\n",
            min(thm_cov$thm_frac_jja), max(thm_cov$thm_frac_jja)))

# ── 2d. Merge and prepare model dataset ───────────────────────────────────────
# Join: keep only drought years that also have a thm_frac observation
model_df <- dplyr::inner_join(ann_sev, thm_cov, by = "year") %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(
    year_std = as.numeric(scale(year))  # standardised year (alternative covariate)
  )

cat(sprintf("  Model dataset: %d observations (drought years with covariate)\n",
            nrow(model_df)))

# ── 2e. Full year spine for prediction ────────────────────────────────────────
# Need thm_frac for ALL years 1950–2025 AND the forward window 2026–2055.
# For historical years without a drought event we still need the covariate.
# For future years we extrapolate the OLS trend in thm_frac.

thm_full <- dplyr::right_join(
  thm_cov,
  data.frame(year = HISTORICAL_START:HISTORICAL_END),
  by = "year"
)

# Fit OLS trend to extrapolate into future
thm_lm <- lm(thm_frac_jja ~ year, data = thm_cov)
future_years_df <- data.frame(year = ENE_START:ENE_END)
future_years_df$thm_frac_jja <- predict(thm_lm, newdata = future_years_df)
# Clamp extrapolated values to [0, 1.5] to prevent physically implausible values
future_years_df$thm_frac_jja <- pmin(1.5, pmax(0, future_years_df$thm_frac_jja))

# Combine historical + future for prediction
pred_years_df <- dplyr::bind_rows(
  thm_full %>% dplyr::mutate(period = "historical"),
  future_years_df %>% dplyr::mutate(period = "future",
                                    thm_frac_jja_std = NA_real_)
) %>%
  dplyr::arrange(year) %>%
  # Re-standardise using historical mean/sd (not including future extrapolation)
  dplyr::mutate(
    thm_mean = mean(thm_cov$thm_frac_jja, na.rm = TRUE),
    thm_sd   = sd(thm_cov$thm_frac_jja,   na.rm = TRUE),
    thm_frac_jja_std = (thm_frac_jja - thm_mean) / thm_sd
  )

cat(sprintf("  Prediction grid: %d years (%d–%d)\n",
            nrow(pred_years_df),
            min(pred_years_df$year),
            max(pred_years_df$year)))
cat(sprintf("  OLS thm_frac trend: %+.5f/yr  (extrapolated to %d)\n",
            coef(thm_lm)["year"], ENE_END))

####################################################################################
# SECTION 3: PRIOR SPECIFICATION
####################################################################################

cat("\n============================================================\n")
cat("  SECTION 3: Prior specification\n")
cat("============================================================\n")

# Model: S_i ~ Gamma(mu_i, phi)  [brms parameterisation: mu, shape]
#        log(mu_i) = alpha + beta * thm_frac_jja_std_i
#
# Priors (weakly informative, see manuscript Section 3.1):
#
#   alpha ~ Normal(log(median(S)), 1.5)
#     → Centres the intercept at the observed median severity on the log scale.
#     → SD = 1.5 on the log scale allows mu to vary between ~exp(-1.5) and
#       exp(+1.5) times the median — a factor of ~4.5 — which is much wider
#       than the actual range of observed severities and thus genuinely weakly
#       informative.
#
#   beta ~ Normal(0, 0.75)
#     → Encodes prior ignorance about the sign of the trend: severity could
#       plausibly increase or decrease with F_thm.
#     → SD = 0.75 on the log-mu scale means a 1-SD increase in thm_frac
#       multiplies expected severity by exp(0.75) ≈ 2.1 at most — a
#       physically plausible upper bound for a 1-SD shift in the covariate.
#     → This is the key prior that constrains the tail: it prevents beta from
#       diverging toward very large positive values driven solely by the
#       rank-1 outlier, which is the record-shattering selection bias
#       identified in manuscript Section 3.2.4(b).
#
#   shape ~ Gamma(2, 0.5)
#     → Mean = 4, SD ≈ 2.8 on the shape parameter (= 1/CV² for Gamma).
#     → Weakly favours moderate right-skewness (shape 2–8), consistent with
#       empirical drought severity distributions globally (e.g. Stagge et al.
#       2015, Bonsal et al. 2019).
#     → Rules out shape < 0.1 (pathologically heavy tails) and shape > 20
#       (near-symmetric severity distributions), both physically implausible.

log_median_sev <- log(median(model_df$max_severity))

priors_main <- c(
  brms::prior(normal(log_median_sev_placeholder, 1.5), class = b, coef = Intercept),
  brms::prior(normal(0, 0.75), class = b, coef = thm_frac_jja_std),
  brms::prior(gamma(2, 0.5),   class = shape)
)

# Substitute the actual log(median) value into the prior string
# (brms requires a literal; we use the stanvar approach)
priors_main <- c(
  brms::set_prior(sprintf("normal(%.4f, 1.5)", log_median_sev),
                  class = "b", coef = "Intercept"),
  brms::set_prior("normal(0, 0.75)", class = "b", coef = "thm_frac_jja_std"),
  brms::set_prior("gamma(2, 0.5)",   class = "shape")
)

cat(sprintf("  Intercept prior: Normal(%.4f, 1.5)  [= log(median S)]\n",
            log_median_sev))
cat("  Slope prior:     Normal(0, 0.75)\n")
cat("  Shape prior:     Gamma(2, 0.5)  [mean shape ≈ 4]\n")

# ── Sensitivity priors (Section 5.2 of MS2) ───────────────────────────────────
# Three alternative priors bracketing the main choice:
#   S1: Tighter slope prior  — beta ~ N(0, 0.5) — more sceptical of trend
#   S2: Looser slope prior   — beta ~ N(0, 1.5) — more permissive of trend
#   S3: Flat shape prior     — shape ~ Gamma(0.01, 0.01) ≈ weakly proper uniform
SENSITIVITY_PRIORS <- list(
  S1_tight_beta = c(
    brms::set_prior(sprintf("normal(%.4f, 1.5)", log_median_sev),
                    class = "b", coef = "Intercept"),
    brms::set_prior("normal(0, 0.5)", class = "b", coef = "thm_frac_jja_std"),
    brms::set_prior("gamma(2, 0.5)",  class = "shape")
  ),
  S2_loose_beta = c(
    brms::set_prior(sprintf("normal(%.4f, 1.5)", log_median_sev),
                    class = "b", coef = "Intercept"),
    brms::set_prior("normal(0, 1.5)", class = "b", coef = "thm_frac_jja_std"),
    brms::set_prior("gamma(2, 0.5)",  class = "shape")
  ),
  S3_flat_shape = c(
    brms::set_prior(sprintf("normal(%.4f, 1.5)", log_median_sev),
                    class = "b", coef = "Intercept"),
    brms::set_prior("normal(0, 0.75)",   class = "b", coef = "thm_frac_jja_std"),
    brms::set_prior("gamma(0.01, 0.01)", class = "shape")
  )
)

cat("  Sensitivity priors defined: S1 (tight beta), S2 (loose beta), S3 (flat shape)\n")

####################################################################################
# SECTION 4: BAYESIAN MODEL FITTING
####################################################################################

cat("\n============================================================\n")
cat("  SECTION 4: Bayesian model fitting\n")
cat("============================================================\n")

brms_rds <- file.path(ENE_DIR, "brms_model_spei3.rds")

if (file.exists(brms_rds)) {
  cat("  Loading cached brms model from disk...\n")
  fit_main <- readRDS(brms_rds)
} else {
  cat("  Fitting Bayesian Gamma GLM (this may take 5–15 minutes)...\n")
  cat(sprintf("  Model: log(mu) = alpha + beta * thm_frac_jja_std\n"))
  cat(sprintf("  Data:  %d observations, seed = %d\n", nrow(model_df), SEED))
  
  fit_main <- brms::brm(
    formula  = max_severity ~ 0 + Intercept + thm_frac_jja_std,
    data     = model_df,
    family   = brmsfamily("gamma", link = "log"),
    prior    = priors_main,
    chains   = N_CHAINS,
    iter     = N_ITER,
    warmup   = N_WARMUP,
    cores    = N_CORES,
    seed     = SEED,
    control  = list(adapt_delta = ADAPT_DELTA, max_treedepth = 12),
    # Suppress verbose Stan output in interactive sessions
    silent   = 2
  )
  
  saveRDS(fit_main, brms_rds)
  cat(sprintf("  ✓ Model fitted and saved: %s\n", basename(brms_rds)))
}

# ── Model diagnostics ──────────────────────────────────────────────────────────
cat("\n  Model summary:\n")
print(summary(fit_main))

# Rhat and ESS
rhat_vals <- brms::rhat(fit_main)
ess_vals  <- brms::neff_ratio(fit_main)
cat(sprintf("  Max Rhat: %.4f  (should be < 1.01)\n", max(rhat_vals, na.rm = TRUE)))
cat(sprintf("  Min ESS ratio: %.3f  (should be > 0.1)\n",
            min(ess_vals, na.rm = TRUE)))
if (max(rhat_vals, na.rm = TRUE) > 1.01)
  warning("⚠ Rhat > 1.01 — consider increasing N_ITER or ADAPT_DELTA.")

# LOO-CV
cat("\n  Computing LOO-CV (pareto-k diagnostics)...\n")
loo_main <- brms::loo(fit_main, moment_match = TRUE)
print(loo_main)
saveRDS(loo_main, file.path(ENE_DIR, "loo_main.rds"))

# Compare with intercept-only (stationary) model
fit_stat <- brms::brm(
  formula = max_severity ~ 1,
  data    = model_df,
  family  = brmsfamily("gamma", link = "log"),
  prior   = c(
    brms::set_prior(sprintf("normal(%.4f, 1.5)", log_median_sev), class = "Intercept"),
    brms::set_prior("gamma(2, 0.5)", class = "shape")
  ),
  chains  = N_CHAINS, iter = N_ITER, warmup = N_WARMUP,
  cores   = N_CORES, seed = SEED, silent = 2
)
loo_stat <- brms::loo(fit_stat, moment_match = TRUE)
loo_comp <- brms::loo_compare(loo_main, loo_stat)
cat("\n  LOO comparison (non-stationary vs stationary model):\n")
print(loo_comp)
write.csv(as.data.frame(loo_comp),
          file.path(ENE_DIR, "loo_comparison.csv"))

# ── Extract posterior draws ────────────────────────────────────────────────────
draws_df <- posterior::as_draws_df(fit_main) %>%
  dplyr::select(
    .draw,
    alpha    = b_Intercept,
    beta     = b_thm_frac_jja_std,
    phi      = shape        # Gamma shape = 1/CV²
  )

write.csv(draws_df, file.path(ENE_DIR, "posterior_draws_spei3.csv"),
          row.names = FALSE)
cat(sprintf("\n  ✓ %d posterior draws extracted\n", nrow(draws_df)))
cat(sprintf("  alpha: median = %.4f  [95%% CI: %.4f – %.4f]\n",
            median(draws_df$alpha),
            quantile(draws_df$alpha, 0.025),
            quantile(draws_df$alpha, 0.975)))
cat(sprintf("  beta:  median = %.4f  [95%% CI: %.4f – %.4f]\n",
            median(draws_df$beta),
            quantile(draws_df$beta, 0.025),
            quantile(draws_df$beta, 0.975)))
cat(sprintf("  phi:   median = %.4f  [95%% CI: %.4f – %.4f]\n",
            median(draws_df$phi),
            quantile(draws_df$phi, 0.025),
            quantile(draws_df$phi, 0.975)))

####################################################################################
# SECTION 5: POSTERIOR PREDICTIVE CHECKS
####################################################################################

cat("\n============================================================\n")
cat("  SECTION 5: Posterior predictive checks\n")
cat("============================================================\n")

ppc_draws <- brms::posterior_predict(fit_main, ndraws = 500)

# PPC: density overlay (observed vs replicated)
ppc_dens <- bayesplot::ppc_dens_overlay(
  y    = model_df$max_severity,
  yrep = ppc_draws
) + ggplot2::labs(
  title    = "Posterior predictive check: density overlay",
  subtitle = sprintf("SPEI-3 annual maximum severity  |  %d drought years  |  500 posterior draws",
                     nrow(model_df)),
  x = "Annual maximum severity (index-units·months)"
) + ggplot2::theme_classic(base_size = 11)

# PPC: maximum statistic — key diagnostic for the outlier event
ppc_max <- bayesplot::ppc_stat(
  y    = model_df$max_severity,
  yrep = ppc_draws,
  stat = "max"
) + ggplot2::labs(
  title    = "Posterior predictive check: maximum severity",
  subtitle = "Does the model reproduce the 2022–2025 record value?",
  x = "Maximum annual severity"
) + ggplot2::theme_classic(base_size = 11)

# PPC: 90th percentile
ppc_p90 <- bayesplot::ppc_stat(
  y    = model_df$max_severity,
  yrep = ppc_draws,
  stat = function(x) quantile(x, 0.90)
) + ggplot2::labs(
  title    = "PPC: 90th percentile",
  x = "90th percentile severity"
) + ggplot2::theme_classic(base_size = 11)

ppc_fig <- (ppc_dens / (ppc_max | ppc_p90)) +
  patchwork::plot_annotation(
    title = "Posterior Predictive Checks — Bayesian Gamma Model",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12))
  )

ggplot2::ggsave(file.path(ENE_DIR, "Fig_PPC.pdf"), ppc_fig,
                width = 10, height = 8, device = "pdf")
ggplot2::ggsave(file.path(ENE_DIR, "Fig_PPC.png"), ppc_fig,
                width = 10, height = 8, dpi = 300)
cat("  ✓ PPC figure saved\n")

####################################################################################
# SECTION 6: POSTERIOR PREDICTIVE ENE AND T_eff(t) COMPUTATION
####################################################################################

cat("\n============================================================\n")
cat("  SECTION 6: Posterior predictive ENE and T_eff(t)\n")
cat("============================================================\n")

# ─────────────────────────────────────────────────────────────────────────────
# STATISTICAL BACKGROUND (reproduced here for reproducibility documentation):
#
# For each posterior draw s = 1,...,S_draws and each year t:
#
#   mu_s(t)  = exp(alpha_s + beta_s * F_thm_std(t))
#            = expected severity given the covariate at year t
#
#   phi_s    = Gamma shape parameter (= 1/CV^2; constant across t)
#
#   p_s(t)   = P(S > S_target | mu_s(t), phi_s)
#            = pgamma(S_target, shape=phi_s, rate=phi_s/mu_s(t), lower.tail=FALSE)
#            = probability that a SINGLE EVENT in year t exceeds S_target
#
#   P_annual,s(t) = 1 - (1 - p_s(t))^lambda
#                 ≈ lambda * p_s(t)  when p_s(t) << 1
#   where lambda = events/year from the Poisson event rate estimated in w5.
#
#   T_eff,s(t) = 1 / P_annual,s(t)
#
#   ENE_s(t0, m) = sum_{t=t0}^{t0+m-1} P_annual,s(t)
#
# Posterior summaries: median + 2.5/25/75/97.5 percentiles across draws s.
# ─────────────────────────────────────────────────────────────────────────────

# Standardise the covariate using historical parameters (from Section 2)
thm_mean_hist <- mean(thm_cov$thm_frac_jja, na.rm = TRUE)
thm_sd_hist   <- sd(thm_cov$thm_frac_jja,   na.rm = TRUE)

# Helper: standardise a raw thm_frac value using historical parameters
std_thm <- function(x) (x - thm_mean_hist) / thm_sd_hist

# Standardise prediction grid
pred_years_df <- pred_years_df %>%
  dplyr::mutate(
    thm_frac_jja_std_hist = std_thm(thm_frac_jja)
  )

# Extract posterior parameter matrix (S_draws × 3)
alpha_vec <- draws_df$alpha
beta_vec  <- draws_df$beta
phi_vec   <- draws_df$phi
S_DRAWS   <- nrow(draws_df)

cat(sprintf("  Computing T_eff over %d years × %d posterior draws...\n",
            nrow(pred_years_df), S_DRAWS))

# Pre-allocate: rows = years, cols = draws
years_pred    <- pred_years_df$year
thm_std_pred  <- pred_years_df$thm_frac_jja_std_hist

n_years_pred  <- length(years_pred)
P_annual_mat  <- matrix(NA_real_, nrow = n_years_pred, ncol = S_DRAWS)

pb_step <- max(1L, S_DRAWS %/% 20L)   # progress every 5%
for (s in seq_len(S_DRAWS)) {
  if (s %% pb_step == 0)
    cat(sprintf("    Draw %d / %d (%.0f%%)\r", s, S_DRAWS, 100*s/S_DRAWS))
  
  mu_t  <- exp(alpha_vec[s] + beta_vec[s] * thm_std_pred)
  rate_t <- phi_vec[s] / mu_t   # Gamma rate = shape / mean
  
  # P(S > S_target | Gamma(shape=phi, rate=phi/mu))
  p_event_t <- pgamma(TARGET_SEVERITY, shape = phi_vec[s],
                      rate = rate_t, lower.tail = FALSE)
  
  # Annual exceedance probability (compound Poisson)
  P_annual_mat[, s] <- 1 - (1 - p_event_t)^LAMBDA_EVENTS
}
cat(sprintf("    Draw %d / %d (100%%)  — done.\n", S_DRAWS, S_DRAWS))

# Derived quantities
T_eff_mat <- 1 / P_annual_mat   # effective return period (years)

# Posterior summary function
post_summary <- function(mat, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  t(apply(mat, 1, quantile, probs = probs, na.rm = TRUE))
}

T_eff_summ <- post_summary(T_eff_mat)
P_ann_summ <- post_summary(P_annual_mat)

# Assemble T_eff trajectory data frame
T_eff_df <- data.frame(
  year          = years_pred,
  thm_frac_jja  = pred_years_df$thm_frac_jja,
  period        = pred_years_df$period,
  T_eff_lo95    = T_eff_summ[, "2.5%"],
  T_eff_lo50    = T_eff_summ[, "25%"],
  T_eff_med     = T_eff_summ[, "50%"],
  T_eff_hi50    = T_eff_summ[, "75%"],
  T_eff_hi95    = T_eff_summ[, "97.5%"],
  P_ann_med     = P_ann_summ[, "50%"],
  P_ann_lo95    = P_ann_summ[, "2.5%"],
  P_ann_hi95    = P_ann_summ[, "97.5%"]
)

write.csv(T_eff_df, file.path(ENE_DIR, "T_eff_trajectory.csv"),
          row.names = FALSE)
cat(sprintf("  ✓ T_eff trajectory saved (%d rows)\n", nrow(T_eff_df)))

# ── ENE over the forward planning window ─────────────────────────────────────
future_idx <- which(years_pred >= ENE_START & years_pred <= ENE_END)
ENE_draws  <- colSums(P_annual_mat[future_idx, , drop = FALSE], na.rm = TRUE)

ENE_summary <- quantile(ENE_draws,
                        probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))
cat(sprintf("\n  ENE (%d–%d)  posterior summary:\n", ENE_START, ENE_END))
cat(sprintf("    2.5%%  : %.3f  (~1-in-%.0f yr)\n",
            ENE_summary["2.5%"],  (ENE_END-ENE_START+1) / ENE_summary["2.5%"]))
cat(sprintf("    25%%   : %.3f  (~1-in-%.0f yr)\n",
            ENE_summary["25%"],   (ENE_END-ENE_START+1) / ENE_summary["25%"]))
cat(sprintf("    50%%   : %.3f  (~1-in-%.0f yr)\n",
            ENE_summary["50%"],   (ENE_END-ENE_START+1) / ENE_summary["50%"]))
cat(sprintf("    75%%   : %.3f  (~1-in-%.0f yr)\n",
            ENE_summary["75%"],   (ENE_END-ENE_START+1) / ENE_summary["75%"]))
cat(sprintf("    97.5%% : %.3f  (~1-in-%.0f yr)\n",
            ENE_summary["97.5%"], (ENE_END-ENE_START+1) / ENE_summary["97.5%"]))

ene_out <- data.frame(
  planning_window  = sprintf("%d-%d", ENE_START, ENE_END),
  n_years          = ENE_END - ENE_START + 1L,
  ENE_lo95         = ENE_summary["2.5%"],
  ENE_lo50         = ENE_summary["25%"],
  ENE_median       = ENE_summary["50%"],
  ENE_hi50         = ENE_summary["75%"],
  ENE_hi95         = ENE_summary["97.5%"],
  implied_T_lo95   = (ENE_END-ENE_START+1) / ENE_summary["97.5%"],
  implied_T_median = (ENE_END-ENE_START+1) / ENE_summary["50%"],
  implied_T_hi95   = (ENE_END-ENE_START+1) / ENE_summary["2.5%"]
)
write.csv(ene_out, file.path(ENE_DIR, "ENE_posterior.csv"), row.names = FALSE)

# ── Four key time-point estimates ─────────────────────────────────────────────
get_timepoint <- function(yr) {
  idx <- which(years_pred == yr)
  if (length(idx) == 0) return(NULL)
  T_draws <- T_eff_mat[idx, ]
  P_draws <- P_annual_mat[idx, ]
  data.frame(
    year      = yr,
    T_eff_lo95 = quantile(T_draws, 0.025, na.rm = TRUE),
    T_eff_lo50 = quantile(T_draws, 0.25,  na.rm = TRUE),
    T_eff_med  = quantile(T_draws, 0.5,   na.rm = TRUE),
    T_eff_hi50 = quantile(T_draws, 0.75,  na.rm = TRUE),
    T_eff_hi95 = quantile(T_draws, 0.975, na.rm = TRUE),
    P_ann_med  = quantile(P_draws, 0.5,   na.rm = TRUE),
    label      = c(T_BASELINE  = "Pre-warming baseline (WMO 1961-1990 midpoint)",
                   T_OCCURRED  = "As-occurred (event mid-year)",
                   T_CURRENT   = "Current / end of study period")[as.character(yr)]
  )
}

four_pts <- dplyr::bind_rows(
  get_timepoint(T_BASELINE),
  get_timepoint(T_OCCURRED),
  get_timepoint(T_CURRENT)
)

# Add ENE row
four_pts <- dplyr::bind_rows(
  four_pts,
  data.frame(
    year       = NA,
    T_eff_lo95 = ene_out$implied_T_hi95,
    T_eff_lo50 = NA,
    T_eff_med  = ene_out$implied_T_median,
    T_eff_hi50 = NA,
    T_eff_hi95 = ene_out$implied_T_lo95,
    P_ann_med  = NA,
    label      = sprintf("30-yr ENE (%d-%d): %.2f [%.2f–%.2f] exceedances",
                         ENE_START, ENE_END,
                         ene_out$ENE_median,
                         ene_out$ENE_lo95,
                         ene_out$ENE_hi95)
  )
)

# Add stationary T_K reference
four_pts <- dplyr::bind_rows(
  four_pts,
  data.frame(
    year       = NA,
    T_eff_lo95 = CI_LO_STAT,
    T_eff_lo50 = NA,
    T_eff_med  = T_K_STATIONARY,
    T_eff_hi50 = NA,
    T_eff_hi95 = CI_HI_STAT,
    P_ann_med  = 1/T_K_STATIONARY,
    label      = "Stationary T_K (MS1 reference)"
  )
)

write.csv(four_pts, file.path(ENE_DIR, "four_timepoints.csv"),
          row.names = FALSE)

cat("\n  Four key time-point estimates:\n")
cat(sprintf("  %-42s  T_eff_med  [95%% CI]\n", "Label"))
cat("  ", paste(rep("-", 70), collapse = ""), "\n", sep = "")
for (i in seq_len(nrow(four_pts))) {
  row <- four_pts[i, ]
  cat(sprintf("  %-42s  %.0f yr  [%.0f – %.0f]\n",
              substr(row$label, 1, 42),
              row$T_eff_med, row$T_eff_lo95, row$T_eff_hi95))
}

####################################################################################
# SECTION 7: FIGURES
####################################################################################

cat("\n============================================================\n")
cat("  SECTION 7: Figures\n")
cat("============================================================\n")

THEME_MS2 <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey35", hjust = 0,
                                          margin = ggplot2::margin(b = 4)),
    panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "bottom",
    plot.margin      = ggplot2::margin(6, 12, 4, 6)
  )

COL_HIST    <- "#2c7bb6"   # blue — historical T_eff
COL_FUTURE  <- "#d7191c"   # red  — future ENE / extrapolated
COL_STAT    <- "grey45"    # grey — stationary T_K reference

# ── Figure 1: T_eff(t) trajectory 1950–2055 ─────────────────────────────────
hist_T  <- T_eff_df %>% dplyr::filter(period == "historical")
fut_T   <- T_eff_df %>% dplyr::filter(period == "future")

# Key time-point annotations
annot_df <- four_pts %>%
  dplyr::filter(!is.na(year)) %>%
  dplyr::mutate(
    label_short = c("Baseline\n(1987)", "Occurred\n(2023)", "Current\n(2025)")
  )

fig_teff <- ggplot2::ggplot() +
  # Stationary T_K horizontal reference
  ggplot2::geom_hline(yintercept = T_K_STATIONARY,
                      linetype = "longdash", colour = COL_STAT, linewidth = 0.8) +
  ggplot2::annotate("text", x = HISTORICAL_START + 1, y = T_K_STATIONARY * 1.12,
                    label = sprintf("Stationary T_K = %.0f yr (MS1)", T_K_STATIONARY),
                    hjust = 0, size = 3.2, colour = COL_STAT) +
  # 2022-2025 drought highlight band
  ggplot2::annotate("rect", xmin = 2022, xmax = 2025.5,
                    ymin = 1, ymax = Inf,
                    fill = "#FEF0D9", alpha = 0.6) +
  # Historical 95% CI ribbon
  ggplot2::geom_ribbon(data = hist_T,
                       ggplot2::aes(x = year, ymin = T_eff_lo95, ymax = T_eff_hi95),
                       fill = COL_HIST, alpha = 0.15) +
  ggplot2::geom_ribbon(data = hist_T,
                       ggplot2::aes(x = year, ymin = T_eff_lo50, ymax = T_eff_hi50),
                       fill = COL_HIST, alpha = 0.25) +
  # Historical median line
  ggplot2::geom_line(data = hist_T,
                     ggplot2::aes(x = year, y = T_eff_med),
                     colour = COL_HIST, linewidth = 1.1) +
  # Future (extrapolated) median + shading
  ggplot2::geom_ribbon(data = fut_T,
                       ggplot2::aes(x = year, ymin = T_eff_lo95, ymax = T_eff_hi95),
                       fill = COL_FUTURE, alpha = 0.12) +
  ggplot2::geom_line(data = fut_T,
                     ggplot2::aes(x = year, y = T_eff_med),
                     colour = COL_FUTURE, linewidth = 1.0, linetype = "dashed") +
  # Vertical dotted lines at key time-points
  ggplot2::geom_vline(xintercept = c(T_BASELINE, T_OCCURRED, T_CURRENT),
                      linetype = "dotted", colour = "grey30", linewidth = 0.6) +
  # Key time-point error bars
  ggplot2::geom_errorbar(
    data = annot_df,
    ggplot2::aes(x = year, ymin = T_eff_lo95, ymax = T_eff_hi95),
    width = 0.8, colour = COL_FUTURE, linewidth = 1.0
  ) +
  ggplot2::geom_point(
    data = annot_df,
    ggplot2::aes(x = year, y = T_eff_med),
    size = 3.5, colour = COL_FUTURE, shape = 21, fill = "white", stroke = 1.5
  ) +
  ggplot2::geom_label(
    data = annot_df,
    ggplot2::aes(x = year, y = T_eff_med,
                 label = sprintf("%s\n%.0f yr", label_short, T_eff_med)),
    vjust = -0.4, hjust = 0.5, size = 2.8, fill = "white",
    label.padding = ggplot2::unit(0.15, "lines"), label.size = 0.2
  ) +
  # Scales
  ggplot2::scale_y_log10(
    labels = scales::comma,
    breaks = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 50000),
    name   = "Effective return period (years, log scale)"
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(seq(1950, 2020, 10), 2025, 2040, 2055),
    name   = "Year"
  ) +
  ggplot2::labs(
    title    = paste0("Effective return period T\u2091\u2090\u2090(t) of the 2022\u20132025 ",
                      "Nechako drought severity (SPEI-3)"),
    subtitle = sprintf(
      paste0("Bayesian non-stationary Gamma model  |  covariate: JJA thermodynamic fraction  |  ",
             "target: S = %.3f index-units\u00b7months\n",
             "Blue = historical (1950\u20132025)  |  Red dashed = future (2026\u20132055, extrapolated)  |  ",
             "Ribbon = 95%% posterior CI"),
      TARGET_SEVERITY
    )
  ) +
  THEME_MS2

ggplot2::ggsave(file.path(ENE_DIR, "Fig_T_eff_trajectory.pdf"),
                fig_teff, width = 12, height = 6, device = "pdf")
ggplot2::ggsave(file.path(ENE_DIR, "Fig_T_eff_trajectory.png"),
                fig_teff, width = 12, height = 6, dpi = 300)
cat("  ✓ Fig_T_eff_trajectory saved\n")

# ── Figure 2: ENE posterior distribution ────────────────────────────────────
ene_density_df <- data.frame(ENE = ENE_draws)
ene_median <- median(ENE_draws)
ene_ci95   <- quantile(ENE_draws, c(0.025, 0.975))

fig_ene <- ggplot2::ggplot(ene_density_df, ggplot2::aes(x = ENE)) +
  ggplot2::geom_histogram(ggplot2::aes(y = after_stat(density)),
                          bins = 60, fill = COL_FUTURE, alpha = 0.7,
                          colour = "white", linewidth = 0.2) +
  ggplot2::geom_density(colour = COL_FUTURE, linewidth = 1.0, adjust = 1.5) +
  ggplot2::geom_vline(xintercept = ene_median, colour = "black",
                      linewidth = 1.0, linetype = "solid") +
  ggplot2::geom_vline(xintercept = ene_ci95, colour = "black",
                      linewidth = 0.7, linetype = "dashed") +
  ggplot2::annotate("text", x = ene_median * 1.04, y = Inf,
                    label = sprintf("Median\n%.2f", ene_median),
                    vjust = 1.5, hjust = 0, size = 3.5) +
  ggplot2::annotate("text", x = ene_ci95[1] * 0.96, y = Inf,
                    label = sprintf("2.5%%\n%.2f", ene_ci95[1]),
                    vjust = 1.5, hjust = 1, size = 3.0, colour = "grey30") +
  ggplot2::annotate("text", x = ene_ci95[2] * 1.02, y = Inf,
                    label = sprintf("97.5%%\n%.2f", ene_ci95[2]),
                    vjust = 1.5, hjust = 0, size = 3.0, colour = "grey30") +
  ggplot2::labs(
    title    = sprintf("Posterior distribution of ENE(%d\u2013%d)", ENE_START, ENE_END),
    subtitle = sprintf(
      paste0("Expected number of exceedances of S = %.3f over the %d-year planning horizon\n",
             "Median: %.2f  [95%% CI: %.2f \u2013 %.2f]  =>  ",
             "~1 recurrence per %.0f yr [CI: %.0f \u2013 %.0f yr]"),
      TARGET_SEVERITY,
      ENE_END - ENE_START + 1,
      ene_median, ene_ci95[1], ene_ci95[2],
      (ENE_END-ENE_START+1)/ene_median,
      (ENE_END-ENE_START+1)/ene_ci95[2],
      (ENE_END-ENE_START+1)/ene_ci95[1]
    ),
    x = sprintf("ENE(%d\u2013%d) — number of exceedances", ENE_START, ENE_END),
    y = "Posterior density"
  ) +
  THEME_MS2

ggplot2::ggsave(file.path(ENE_DIR, "Fig_ENE_posterior.pdf"),
                fig_ene, width = 9, height = 5, device = "pdf")
ggplot2::ggsave(file.path(ENE_DIR, "Fig_ENE_posterior.png"),
                fig_ene, width = 9, height = 5, dpi = 300)
cat("  ✓ Fig_ENE_posterior saved\n")

####################################################################################
# SECTION 8: PRIOR SENSITIVITY ANALYSIS
####################################################################################

cat("\n============================================================\n")
cat("  SECTION 8: Prior sensitivity analysis\n")
cat("============================================================\n")

sens_results <- vector("list", length(SENSITIVITY_PRIORS))
names(sens_results) <- names(SENSITIVITY_PRIORS)

for (prior_name in names(SENSITIVITY_PRIORS)) {
  sens_rds <- file.path(ENE_DIR, sprintf("brms_sens_%s.rds", prior_name))
  cat(sprintf("  Fitting sensitivity model: %s\n", prior_name))
  
  if (file.exists(sens_rds)) {
    fit_sens <- readRDS(sens_rds)
  } else {
    fit_sens <- brms::brm(
      formula = max_severity ~ 0 + Intercept + thm_frac_jja_std,
      data    = model_df,
      family  = brmsfamily("gamma", link = "log"),
      prior   = SENSITIVITY_PRIORS[[prior_name]],
      chains  = N_CHAINS, iter = N_ITER, warmup = N_WARMUP,
      cores   = N_CORES, seed = SEED, silent = 2
    )
    saveRDS(fit_sens, sens_rds)
  }
  
  draws_sens <- posterior::as_draws_df(fit_sens) %>%
    dplyr::select(.draw, alpha = b_Intercept,
                  beta = b_thm_frac_jja_std, phi = shape)
  
  # Compute ENE for this sensitivity model
  P_ann_sens <- matrix(NA_real_, nrow = n_years_pred, ncol = nrow(draws_sens))
  for (s in seq_len(nrow(draws_sens))) {
    mu_t   <- exp(draws_sens$alpha[s] + draws_sens$beta[s] * thm_std_pred)
    p_ev   <- pgamma(TARGET_SEVERITY, shape = draws_sens$phi[s],
                     rate = draws_sens$phi[s] / mu_t, lower.tail = FALSE)
    P_ann_sens[, s] <- 1 - (1 - p_ev)^LAMBDA_EVENTS
  }
  ene_sens <- colSums(P_ann_sens[future_idx, , drop = FALSE], na.rm = TRUE)
  
  sens_results[[prior_name]] <- data.frame(
    prior_label   = prior_name,
    beta_med      = median(draws_sens$beta),
    beta_lo95     = quantile(draws_sens$beta, 0.025),
    beta_hi95     = quantile(draws_sens$beta, 0.975),
    phi_med       = median(draws_sens$phi),
    ENE_med       = median(ene_sens),
    ENE_lo95      = quantile(ene_sens, 0.025),
    ENE_hi95      = quantile(ene_sens, 0.975),
    T_implied_med = (ENE_END-ENE_START+1) / median(ene_sens)
  )
  cat(sprintf("    β  median = %.4f  ENE median = %.3f\n",
              median(draws_sens$beta), median(ene_sens)))
}

# Add main model result for comparison
main_ene_summ <- data.frame(
  prior_label   = "Main (N(0, 0.75))",
  beta_med      = median(draws_df$beta),
  beta_lo95     = quantile(draws_df$beta, 0.025),
  beta_hi95     = quantile(draws_df$beta, 0.975),
  phi_med       = median(draws_df$phi),
  ENE_med       = median(ENE_draws),
  ENE_lo95      = quantile(ENE_draws, 0.025),
  ENE_hi95      = quantile(ENE_draws, 0.975),
  T_implied_med = (ENE_END-ENE_START+1) / median(ENE_draws)
)

sens_table <- dplyr::bind_rows(
  main_ene_summ,
  dplyr::bind_rows(sens_results)
) %>% dplyr::arrange(prior_label)

write.csv(sens_table, file.path(ENE_DIR, "sensitivity_priors.csv"),
          row.names = FALSE)
cat("\n  Prior sensitivity table:\n")
print(sens_table)

# ── Sensitivity figure ────────────────────────────────────────────────────────
fig_sens <- sens_table %>%
  dplyr::mutate(prior_label = factor(prior_label,
                                     levels = rev(c("Main (N(0, 0.75))",
                                                    "S1_tight_beta", "S2_loose_beta",
                                                    "S3_flat_shape")))) %>%
  ggplot2::ggplot(ggplot2::aes(x = prior_label, y = ENE_med)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = ENE_lo95, ymax = ENE_hi95),
                         width = 0.3, linewidth = 1.0, colour = COL_FUTURE) +
  ggplot2::geom_point(size = 4, shape = 21, fill = COL_FUTURE,
                      colour = "white", stroke = 1.5) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title    = "Prior sensitivity: ENE over 2026–2055",
    subtitle = "Point = posterior median; error bars = 95% credible interval",
    x = "Prior specification",
    y = sprintf("ENE(%d\u2013%d)", ENE_START, ENE_END)
  ) +
  THEME_MS2

ggplot2::ggsave(file.path(ENE_DIR, "Fig_prior_sensitivity.pdf"),
                fig_sens, width = 8, height = 4, device = "pdf")
ggplot2::ggsave(file.path(ENE_DIR, "Fig_prior_sensitivity.png"),
                fig_sens, width = 8, height = 4, dpi = 300)
cat("  ✓ Fig_prior_sensitivity saved\n")

####################################################################################
# SECTION 9: FINAL SUMMARY
####################################################################################

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║  w12_bayesian_ene_analysis.R  —  COMPLETE                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("KEY RESULTS FOR MS2:\n")
cat(sprintf("  Stationary T_K (MS1 reference):  %.0f yr  [%.0f – %.0f]\n",
            T_K_STATIONARY, CI_LO_STAT, CI_HI_STAT))
cat(sprintf("  T_eff at %d (pre-warming):      %.0f yr  [%.0f – %.0f]\n",
            T_BASELINE, annot_df$T_eff_med[1],
            annot_df$T_eff_lo95[1], annot_df$T_eff_hi95[1]))
cat(sprintf("  T_eff at %d (as-occurred):      %.0f yr  [%.0f – %.0f]\n",
            T_OCCURRED, annot_df$T_eff_med[2],
            annot_df$T_eff_lo95[2], annot_df$T_eff_hi95[2]))
cat(sprintf("  T_eff at %d (current):          %.0f yr  [%.0f – %.0f]\n",
            T_CURRENT, annot_df$T_eff_med[3],
            annot_df$T_eff_lo95[3], annot_df$T_eff_hi95[3]))
cat(sprintf("  ENE %d-%d:                    %.2f  [%.2f – %.2f]\n",
            ENE_START, ENE_END,
            ENE_summary["50%"], ENE_summary["2.5%"], ENE_summary["97.5%"]))
cat(sprintf("  Implied recurrence:             ~%.0f yr  [%.0f – %.0f]\n",
            (ENE_END-ENE_START+1)/ENE_summary["50%"],
            (ENE_END-ENE_START+1)/ENE_summary["97.5%"],
            (ENE_END-ENE_START+1)/ENE_summary["2.5%"]))
cat(sprintf("  Compression factor (1987→2025): %.1fx\n",
            annot_df$T_eff_med[1] / annot_df$T_eff_med[3]))

cat("\nOUTPUT FILES (", ENE_DIR, "):\n", sep = "")
outputs <- c(
  "brms_model_spei3.rds       — fitted brms model (reload with readRDS)",
  "posterior_draws_spei3.csv  — tidy posterior (alpha, beta, phi)",
  "T_eff_trajectory.csv       — T_eff(t) 1950-2055 with 95% CI",
  "ENE_posterior.csv          — ENE summary statistics",
  "four_timepoints.csv        — T_eff at 1987/2023/2025 + ENE",
  "sensitivity_priors.csv     — prior sensitivity table",
  "Fig_T_eff_trajectory.pdf/png",
  "Fig_ENE_posterior.pdf/png",
  "Fig_PPC.pdf/png",
  "Fig_prior_sensitivity.pdf/png"
)
for (o in outputs) cat(sprintf("  %s\n", o))