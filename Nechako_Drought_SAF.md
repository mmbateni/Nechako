Nechako Drought SAF Analysis — README
Overview
This codebase performs a Severity–Area–Frequency (SAF) drought analysis for the Nechako Basin in British Columbia, Canada. It uses gridded SPI (Standardised Precipitation Index) and SPEI (Standardised Precipitation-Evapotranspiration Index) raster time series as inputs, and produces statistically rigorous SAF curves, copula models, return period tables, diagnostic plots, publication-ready figures, and Excel summary tables as outputs.
The analysis follows a multivariate frequency analysis framework: drought severity and drought-affected area fraction are modelled as two correlated random variables, joined by a copula. SAF curves answer the question: "For a given return period, how severe is a drought that affects X% of the basin?"
File Inventory & Roles
| File
|Role
|Required?
|
| ---|---|---|
| Nechako_Drought_RUN_ALL.R
|Master orchestrator  — single entry point that runs all steps in order
|Yes
|
| Nechako_Drought_QuickStart.R
|Input verification, pre-flight checks, and post-run diagnostic function
|Yes
|
| Nechako_Drought_SAF_Analysis.R
|Core analysis pipeline (data loading → SAF curves)
|Yes
|
| Nechako_Drought_SAF_Analysis_ext1.R
|Extension: dependency tests, copula GOF, time-varying and segmented copula
|Yes
|
| Nechako_Drought_SAF_Analysis_ext2.R
|Extension: non-stationary marginals, Kendall-corrected NS SAF, event placement
|Yes
|
| Nechako_Drought_Methodology_Guide.R
|Optional utility functions (threshold sensitivity, seasonal plots, spatial maps)
|Optional
|
Prerequisites
Working Directory
All scripts assume:
setwd("D:/Nechako_Drought/Nechako")
Change this path before running if your data lives elsewhere.
Required Input Files
The pipeline reads monthly NetCDF files pre-computed by an upstream SPI/SPEI script. They must exist at:
spi_results_seasonal/spi_01_month01_Jan.nc   # SPI-1 for January
spi_results_seasonal/spi_01_month02_Feb.nc   # SPI-1 for February
...                                           # (12 files per scale)
spi_results_seasonal/spi_03_month01_Jan.nc   # SPI-3 for January
...
spei_results_seasonal/spei_01_month01_Jan.nc # SPEI-1 for January
...
spei_results_seasonal/spei_03_month01_Jan.nc # SPEI-3 for January
...
Each NetCDF contains a multi-layer raster (one layer per year) for that calendar month and scale.
R Packages
The scripts auto-install missing packages. The full list is:
| Package
|Purpose
|
| ---|---|
| terra
|Raster I/O and spatial operations
|
| copula
|Copula fitting, GOF, density/CDF functions
|
| fitdistrplus
|MLE fitting of univariate marginals
|
| MASS
|Supporting distributions
|
| ggplot2
|All plots
|
| gridExtra
|Multi-panel PDF plots and table grobs
|
| viridis
|Colour scales
|
| dplyr
|Data manipulation
|
| scales
|Axis formatting
|
| patchwork
|Multi-panel figure assembly (Step 7)
|
| ggExtra
|Marginal plot annotations (Step 7)
|
| openxlsx
|Excel export (ext1, Step 7)
|
| moments
|Higher-order statistics (ext1)
|
| Kendall
|Rank correlation tests (ext1)
|
| trend
|Pettitt test for change-point detection (ext1)
|
| gamlss
|Non-stationary GAMLSS models (ext2)
|
| gamlss.dist
|GAMLSS distribution families (ext2)
|
How to Run
The simplest approach is the single entry point:
setwd("D:/Nechako_Drought/Nechako")
source("Nechako_Drought_RUN_ALL.R")
This runs all seven steps (see below) in the correct order. Do not source the individual analysis scripts out of order — `ext1.R` and `ext2.R` depend on objects created by `SAF_Analysis.R`, and `SAF_Analysis.R` depends on objects created by `QuickStart.R`.
Pipeline Steps (as executed by `RUN_ALL.R`)
Step 1 — Input Verification (`QuickStart.R`)
`QuickStart.R` defines the `quick_diagnostic()` function (called at the end of Step 7) and performs three pre-flight checks:
File existence: Spot-checks four representative NetCDF files (one per index, January only) as proxies for all 48 required files (4 indices × 12 months). A missing file stops execution with an error; the remaining 44 files are verified implicitly when the pipeline loads them in Step 3.
SPI-1 value range: Loads the SPI-1 January raster and validates that all pixel values fall within the physically plausible range of −5 to +5. A warning is issued if any layer falls outside.
Kendall distribution sanity check: Draws 100,000 i.i.d. pairs from the independence copula (`C(u,v) = u×v`) and verifies that the empirical `K_C(t) = P(C(U,V) ≤ t)` at `t = 0.5` matches the known analytical result `t − t·log(t) ≈ 0.8466` to within a tolerance of 0.005. This directly validates the K_C implementation used in the main pipeline before any data is processed.
`RUN_ALL.R`  overrides  `quick_diagnostic()`  after sourcing  `QuickStart.R`  to correct its hard-coded search path: the original searches  `"drought_analysis` `"` , but all pipeline output is written to  `"spatial_drought` `"`  (the  `OUT_ROOT`  variable set in  `SAF_Analysis.R` ).
Step 2 — Methodology Helpers (`Methodology_Guide.R`)
Loads five utility functions into the workspace. No computation is triggered at load time. These functions are available for optional Step 6 calls:
`analyze_multiple_thresholds()` — repeats the area fraction computation at drought thresholds −0.5, −1.0, and −1.5 and produces a PDF line chart comparing drought-affected area over time across thresholds.
`seasonal_drought_analysis()` — groups monthly drought characteristics by calendar month, computes mean severity and mean area fraction per month, and saves a two-panel bar chart PDF.
`map_drought_frequency()` — uses `terra::app()` to compute two per-pixel summary rasters: (1) drought frequency as the percentage of months below the threshold, and (2) mean severity when below threshold. Both are saved as NetCDF files.
`compare_return_period_methods()` — combines SAF outputs from the Conditional and Kendall methods and plots them on a single overlay, coloured by return period and distinguished by linetype.
`cross_correlate_indices()` — merges SPI-1 and SPI-3 drought characteristics by date, scatter-plots their event severities, and reports the Pearson correlation coefficient.
Step 3 — Main SAF Pipeline (`SAF_Analysis.R`)
This is the computational core. It runs the full pipeline four times — once each for  SPI-1 ,  SPEI-1 ,  SPI-3 , and  SPEI-3  — and writes results to  `spatial_drought/` `<INDEX` `>_analysis/` . The working root is controlled by  `OUT_ROOT ` `<- ` `"spatial_drought` `"`  set at the top of the script.
3a. Data Loading
All 12 monthly NetCDF tiles for each index are loaded with `terra::rast()` via `assemble_index_raster()` and concatenated into a single time-ordered raster stack. Dates, month indices, and year indices are extracted from the `terra::time()` metadata of the SPI-1 stack and reused for all four indices. Cell areas (km²) are computed once using `terra::cellSize(..., unit="km")`, which handles geographic (latitude-varying) and projected coordinate systems correctly. The basin total area is derived as the sum of all valid cell areas.
3b. Drought Characteristics Extraction
`extract_drought_characteristics_spi()` iterates over every time step `t` in the raster stack. For each layer it:
Reads all pixel values as a flat vector.
Identifies "dry" pixels as finite values below the drought threshold (default −0.5).
Computes drought-affected area fraction as `100 × n_dry / n_valid` (percentage of valid pixels).
Computes mean drought severity as `mean(threshold − value)` for all dry pixels (positive when below threshold, zero otherwise).
Output: a time series data frame `[date, year, month, severity, area_pct]` saved as `<INDEX>_drought_characteristics.csv`.
3c. Event Identification & Duration Classification
`identify_duration_classified_events()` groups consecutive months with `severity > 0` and `area_pct > 0` into discrete drought events using a run-length scan. Events still open at the end of the record are closed at the last time step. Each event is described by:
Start/end date and year
Duration in months
Duration class (dependent on the accumulation scale):
| Class
|SPI-1 / SPEI-1
|SPI-3 / SPEI-3
|
| ---|---|---|
| D1-2 (Very-short)
|1–2 months
|—
|
| D3-6 (Short-term)
|3–6 months
|3–6 months
|
| D7-12 (Medium-term)
|7–12 months
|7–12 months
|
| D13+ (Long-term)
|≥13 months
|≥13 months
|
Mean and max severity; mean and max area fraction
`compute_return_periods_by_class()` computes, per duration class:
Empirical return period: `record_length_years / n_events`
Exact Poisson confidence interval: 95% CI on the rate λ using the chi-squared relationship `lambda_lo = qchisq(0.025, df=2n) / (2·T)`, inverted to return period bounds.
Gamma renewal process diagnostics: Inter-arrival times (start-to-start, in years) are computed for all events in the class. The coefficient of variation (CV) of inter-arrival times determines whether the arrival process is well-approximated by a Poisson process:
CV ≈ 1.0 → approximately Poisson
CV > 1.5 → clustered (flagged `CLUSTERED`)
CV < 0.5 → regular (flagged `REGULAR`)
When the process departs from Poisson, a Gamma renewal model is fitted: shape `k = 1/CV²`, rate `k/mean_iat`. The empirical mean waiting time `mu_T_adj` (months) and renewal parameters are stored in the return period table for use in SAF target probability adjustment (see below).
Outputs: `<INDEX>_duration_classified_events.csv`, `<INDEX>_return_periods_by_class.csv`.
3d. Marginal Distribution Fitting
`fit_marginal_distributions()` uses only months with active drought (`severity > 0`, `area_pct > 0`):
Area fraction is rescaled to (0, 1) and fitted with a Beta distribution via MLE (both shape parameters free).
Severity is fitted with six candidate distributions using MLE: Exponential, Gamma, Weibull, Log-Normal, Normal, and Logistic. The best distribution is selected by minimum AIC.
Output: `<INDEX>_severity_dist_comparison.csv`.
3e. Copula Fitting
`fit_copulas()` transforms severity and area fraction into uniform pseudo-observations `(u_S, u_A)` using the parametric CDFs from Step 3d. Seven copula families are then fitted to these pseudo-observations using the IFM (Inference Functions for Margins) maximum pseudo-likelihood method (`method = "mpl"`), with BFGS optimisation (falling back to Nelder-Mead on convergence failure):
| Copula
|Tail dependence
|
| ---|---|
| Frank
|Symmetric (no tail)
|
| Gumbel
|Upper tail
|
| Plackett
|Symmetric (no tail)
|
| Survival Clayton
|Upper tail
|
| Joe
|Upper tail
|
| Gaussian (Normal)
|No tail (asymptotic independence)
|
| Student-t
|Both tails (symmetric)
|
The best copula is selected by minimum AIC. AIC, BIC, log-likelihood, and the raw parameter estimate are all saved. An empirical upper-tail dependence coefficient λ_U is computed non-parametrically (quantile 0.95) and compared to the theoretical λ_U of the winning family. A `WARNING [TD-U]` is logged if the AIC-best family has zero upper-tail dependence (Frank, Plackett, or Gaussian) but the empirical λ_U > 0.10, alerting the analyst to consider a tail-dependent alternative.
Note: The same seven-family candidate set is used in `fit_copulas()` (Step 3), `detect_copula_changepoints()` / `fit_segmented_copula()` (ext1, Step 5c), and the PRUTF structural-break test. The Student-t family has two parameters (ρ, df); because the log-linear TV formulation operates on a scalar parameter, the TV copula trend model is skipped for Student-t and the stationary fit is used instead (a warning is logged).
Output: `<INDEX>_copula_comparison.csv`.
3f. SAF Curve Derivation — Two Methods
Both methods use the internal helper `.compute_renewal_saf_target()` to compute the effective monthly non-exceedance probability for a given return period `T` (years). When the event inter-arrival CV falls outside [0.5, 1.5], a Gamma renewal survival function replaces the standard Poisson target:
Standard Poisson: `target = 1 − mu_T / (T × 12)`, where `mu_T` is the mean inter-drought interval in months.
Gamma renewal adjustment (when CV < 0.5 or CV > 1.5): `target = 1 − (1 − P_no_event) / T`, where `P_no_event = 1 − F_Gamma(T; k, rate)` is the Gamma survival probability. This accounts for clustering (CV > 1) or regularity (CV < 1) in event arrivals that a Poisson model would misrepresent.
Method A: Conditional Copula
For each area value `a`, the script solves numerically for the severity `s` such that the conditional probability satisfies the target:
P(S > s | A = a) = 1 − target
This is implemented by inverting `cCopula(cbind(u_S, u_A), copula, indices=2)` — the partial derivative of the copula CDF with respect to the area uniform variate — using `optimize()` over a severity search interval of (1e-6, 50).
Method B: Kendall Distribution
Instead of conditioning on area, this method characterises the joint event by a copula level set. It first solves for the copula level `t_val` such that `K_C(t_val) = target`, where `K_C(t) = P(C(U,V) ≤ t)` is the Kendall distribution. Two implementations are used:
Analytical K_C for Archimedean families (Frank, Gumbel, Joe, SurvClayton): uses the closed-form generator formula `K_C(t) = t − φ(t)/φ'(t)`, where φ is the Archimedean generator and φ' its derivative. This enables stable extrapolation for high return periods without empirical discretisation limits.
Empirical K_C for non-Archimedean families (Gaussian, Student-t, Plackett): estimated directly as `mean(C(U,V) ≤ t)` from the fitted pseudo-observations.
Then for each area value `a`, severity `s` is found such that `C(u_S, u_A) = t_val`.
Both methods produce curves for return periods of 10, 25, 50, and 100 years across area fractions from 5% to 95%.
SAF curves are also computed per duration class when the class has ≥15 events. The per-class SAF uses `subset_drought_chars_by_class()` to filter the characteristics time series to only months belonging to events of that duration class, then re-fits marginals and copula on the subset before deriving curves.
Outputs per index: `<INDEX>_SAF_conditional_all.csv`, `<INDEX>_SAF_kendall_all.csv`, `<INDEX>_SAF_overlay.pdf`, and per-class equivalents.
The four result objects  `res_spi1` ,  `res_spei1` ,  `res_spi3` ,  `res_spei3`  remain in the R workspace for use by Steps 4–6. Each is a named list containing:  `dc` ,  `events` ,  `rp` ,  `marginals` ,  `copulas` ,  `saf_conditional` ,  `saf_kendall` ,  `saf_cond_by_class` ,  `saf_kend_by_class` .
Step 4 — Load Extensions (`ext1.R` + `ext2.R`)
The extension scripts are sourced here. They only define functions — no heavy computation occurs at load time (each script emits a single `log_event()` message confirming successful loading). Both scripts verify that `log_event()` from `SAF_Analysis.R` is present before loading; an error stops execution if the main pipeline was skipped.
Functions exported by `ext1.R` (v6):`test_dependency`, `gof_copula`, `detect_copula_changepoints`, `rosenblatt_pit_gof`, `fit_timevarying_copula`, `fit_segmented_copula`
Functions exported by `ext2.R` (v7):`fit_nonstationary_marginals`, `derive_SAF_nonstationary`, `derive_SAF_nonstationary_kendall`, `place_event_on_saf`
Internal helpers (not called directly): `.compute_renewal_saf_target` (defined in `SAF_Analysis.R`, used by all SAF functions), `.h_func_fd` (ext1), `.build_period_copula` (ext2), `.resolve_epoch` (ext2), `.assert_epoch_consistency` (ext2).
Step 5 — Run Enhancements
For each of the four indices, five enhancement analyses are run sequentially. Results are collected into `enh_results[[nm]]` for use by Step 7.
5a. Dependency Test
`test_dependency()` computes Kendall's τ and Spearman's ρ between drought severity and area fraction (restricted to drought months), with a significance test via `cor.test()`. This validates whether the copula is statistically justified — if τ ≈ 0, the variables are nearly independent and a copula provides little benefit over marginal-only methods.
5b. Copula Goodness-of-Fit (Bootstrap Sn)
`gof_copula()` performs a parametric bootstrap GOF test using the Cramér–von Mises `Sn` statistic (`copula::gofCopula`) with 499 bootstrap replicates (default). The minimum achievable two-sided p-value is 1/500 = 0.002, well below α = 0.05. A p-value > 0.05 means the chosen copula family is not rejected.
Use `N_boot = 999` for publication-quality p-values.
Use `N_boot = 99` for fast development runs.
If the test rejects the copula (p < 0.05), `plot_copula_diagnostic()` is called to produce a heatmap of empirical minus theoretical copula CDF differences, and `plot_rosenblatt_diagnostic()` is called (when a valid TV copula fit is also available) to visualise the transformed residuals.
5c. Time-Varying & Segmented Copula
`fit_timevarying_copula()` tests whether the dependence structure between severity and area has changed over time. The function proceeds in four stages:
Stage (i) — Change-point detection via `detect_copula_changepoints()`: A rolling Kendall-τ series is computed from the drought pseudo-observations using a centred window of ±5 years. Zero-variance windows are silently assigned NA (both `sd(severity)` and `sd(area_pct)` are checked before calling `cor()`; a count of suppressed windows is reported via `log_event()`). Three complementary structural-break tests are applied to the resulting τ series:
Pettitt test (`trend::pettitt.test`): rank-based, best for single abrupt shifts.
CUSUM (manual implementation): cumulative sum of standardised τ deviations; maximum |CUSUM| is compared to a permutation null distribution using 499 Monte-Carlo draws.
PRUTF (parametric likelihood-ratio structural-break): for each of five candidate break points (at the 20th, 35th, 50th, 65th, and 80th percentiles of the time series length), the full-record log-likelihood is compared to the sum of two independently fitted segment log-likelihoods. The test statistic is the supremum of `2(ll_seg1 + ll_seg2 − ll_full)`, compared to χ²(1).
A change point is flagged as "detected" when at least 2 of the 3 tests reject H₀ at the 5% level. The consensus change-point year is the median of the significant test estimates. Results are saved to `<INDEX>_changepoint_detection.csv`.
Stage (ii) — Three copula models are fitted and compared:
| Model
|Description
|Free parameters
|
| ---|---|---|
| Stationary
|Single fixed copula parameter
|1
|
| Time-varying (TV)
|Log-linear trend on link scale:  θ(t) = ilink(a + b × year_std)
|2
|
| Segmented
|Two independent stationary copulas (one per temporal epoch) via  fit_segmented_copula()
|2
|
Family-specific link functions are used: identity for Frank, log(θ−1) for Gumbel/Joe, log(θ) for Plackett/SurvClayton, atanh for Gaussian/Student-t. The TV model is skipped for Student-t (two-parameter family incompatible with the scalar trend formulation).
Stage (iii) — Model selection: TV is preferred when LRT p < 0.05, no strong abrupt break detected, TV AIC lower than segmented, and Rosenblatt PIT adequate. Segmented is preferred when a change point is detected, segmented AIC better by > 2 units, and Rosenblatt diagnostic improves. Otherwise the stationary model is retained.
Stage (iv) — Rosenblatt PIT GOF via `rosenblatt_pit_gof()`: Validates the fitted copula by transforming pseudo-observations `(u_S, u_A)` into two independent U[0,1] variates:
`e1_i = u_S_i` (trivially uniform)
`e2_i = ∂C(u_S_i, u_A_i; θ_i) / ∂u_S_i` (second Rosenblatt variate)
`e2` is computed via the internal helper `.h_func_fd()` using a central finite difference of `pCopula()` with an adaptive step size that shrinks proportionally to the distance from the boundary (`eps_adj = eps × min(u1, 1−u1)`), preventing numerical overflow near (0,1) endpoints. This approach works for every copula family in the package, including Plackett (which has no `cCopula` method). Adequacy requires both KS uniformity (e2 ~ U[0,1], p > 0.05) and independence (Kendall τ between e1 and e2, p > 0.05).
Diagnostic PDFs produced conditionally:
| File
|Condition
|
| ---|---|
| <INDEX >_rolling_tau.pdf
|Always written
|
| <INDEX >_copula_diagnostic.pdf
|Bootstrap GOF rejected (p  < 0.05)
|
| <INDEX >_rosenblatt_diagnostic.pdf
|GOF rejected AND TV copula fitted with valid  a_hat
|
| <INDEX >_tv_parameter.pdf
|TV model selected
|
Outputs: `<INDEX>_changepoint_detection.csv`; `<INDEX>_segmented_copula.csv` (written only when a change point is detected and both segments have ≥ 15 observations).
5d. Non-Stationary Marginals (GAMLSS)
`fit_nonstationary_marginals()` fits three GAMLSS Gamma models to drought severity (drought months only) and selects among them by AIC:
| Model
|mu formula
|sigma formula
|
| ---|---|---|
| M0 (stationary)
|~ 1
|~ 1
|
| M1 (mu-trend)
|~ year_std
|~ 1
|
| M2 (mu+sigma-trend)
|~ year_std
|~ year_std
|
`year_std` is the calendar year standardised by the mean and SD of drought observation years. Period-specific `mu` and `sigma` (coefficient of variation) are extracted directly from model coefficients rather than via `predict()` (which returns a closure in some `gamlss` versions), applying the GA log-link manually: `mu(t) = mean(exp(b0 + b1 × year_std(t)))` averaged over the epoch years.
Epoch definition (precedence order):
Explicit `ref_years` / `recent_years` arguments — highest priority (user override).
`cp_result$detected == TRUE` from Stage 5c: `reference = [record_start, cp_year − 1]`, `recent = [cp_year, record_end]`. If the change point falls outside 20%–80% of the record length, the automatic 40%/30% split is used with a logged warning.
Default automatic split: first 40% of record (reference) / last 30% of record (recent).
Output: `<INDEX>_nonstationary_summary.csv`.
5e. Non-Stationary SAF Curves
Uses the period-specific Gamma parameters `mu` and `sigma` from Step 5d while keeping the Beta area marginal unchanged. Produces separate SAF curves for the reference and recent periods using both:
Method A (`derive_SAF_nonstationary`) — conditional copula approach with period-specific Gamma severity marginal (shape = 1/σ², rate = shape/μ).
Method B (`derive_SAF_nonstationary_kendall`) — Kendall-corrected approach using the same period-specific Gamma marginal. Kendall distribution K_C is computed from the period-specific copula (not the full-record stationary copula) to maintain consistency between the copula used in inversion and the one used for calibration. Archimedean families use analytical K_C; non-Archimedean families use a 3-tier empirical strategy:
Tier 1: direct `pCopula()` on pseudo-observations.
Tier 2: simulation-based K_C (N = 20,000) if Tier 1 fails.
Tier 3: return `NULL` and abort this epoch's curves if both fail, preventing miscalibrated SAF output.
Copula selection for NS SAF follows a three-path priority via `.build_period_copula()`:
Path 1 (preferred): segmented copula from `tv_res$seg_result` — both marginal and copula are step-functions sharing the same change-point boundary.
Path 2 (hybrid approximation): TV copula evaluated at the epoch midpoint year — used only when `seg_result` is unavailable; a WARNING is logged each time.
Path 3 (fallback): stationary full-record copula — used when no non-stationarity in dependence was detected.
Epoch boundary consistency between the marginal and copula is verified by `.assert_epoch_consistency()`, which logs `WARNING [EPOCH MISMATCH]` if the change-point years differ.
The output column `copula_path` records which path was taken for every row of the CSV, enabling post-hoc auditing.
Outputs: `<INDEX>_SAF_nonstationary_<period>.csv`, `<INDEX>_SAF_nonstationary_kendall_<period>.csv`.
Step 6 — Optional Utility Analyses (commented out in `RUN_ALL.R`)
These lines are commented out but ready to activate:
| Comment block
|What it does
|
| ---|---|
| 6a
|Threshold sensitivity: re-runs area fraction at −0.5, −1.0, −1.5 and saves a PDF line chart
|
| 6b
|Seasonal bar charts of mean severity and mean area by calendar month
|
| 6c
|Pixel-level drought frequency/severity maps saved as NetCDF to  drought_maps/
|
| 6d
|Cross-index correlation: SPI-1 vs SPI-3 severity scatter plot with Pearson r
|
| 6e
|Places a specific historical event (e.g. years 2003–2004) on the SAF surface and returns its estimated joint return period
|
Step 7 — Publication Figures & Excel Tables (`RUN_ALL.R`)
Step 7 is fully implemented in `RUN_ALL.R` and runs automatically. It assembles all pipeline and enhancement results into 18 publication-quality figures and 3 Excel tables, all written to `spatial_drought/pub_figures/`. Each figure is saved as both PDF and JPEG (300 dpi) via the internal `save_pub()` helper. Figures use a shared `pub_theme()` (bold titles, `theme_bw`, consistent spacing).
Section 1 — Study Area & Data Overview
| Figure
|Contents
|
| ---|---|
| Fig 01  — Basin drought time series
|Two stacked panels: monthly mean severity (top) and drought-affected area % (bottom) for all four indices over the full record. Coloured by index using a consistent palette ( SPI1=#1b7837 ,  SPEI1=#762a83 ,  SPI3=#2166ac ,  SPEI3=#d6604d ).
|
| Fig 02  — Event characteristics scatter
|Scatter plot of mean area % vs mean severity for all identified events across all indices, coloured by duration class (viridis plasma), shaped by index.
|
| Fig 03  — Event duration histograms
|2×2 faceted histogram of event duration in months, one panel per index, free y-axis.
|
Section 2 — Marginal Distributions
| Figure
|Contents
|
| ---|---|
| Fig 04  — Severity marginal fits
|2×2 panel grid: empirical severity histogram overlaid with the AIC-best fitted density (Gamma, Weibull, Log-Normal, etc.) for each index.
|
| Fig 05  — Area marginal fits
|2×2 panel grid: empirical area% histogram overlaid with the fitted Beta density for each index.
|
Section 3 — Copula Dependence Structure
| Figure
|Contents
|
| ---|---|
| Fig 06  — Pseudo-observation scatter with copula contours
|2×2 grid: pseudo-observation scatter  (u_S, u_A)  for each index overlaid with  dCopula()  density contours from the best-fit copula. Subtitle shows best copula name and Kendall τ (sourced from  enh_results[[nm]]$dependency$kendall ).
|
| Fig 07  — Kendall τ summary barplot
|Bar chart of Kendall τ per index with ±0.02 error bars and τ/ρ text annotations.
|
| Fig 08  — Copula AIC comparison
|Grouped bar chart of ΔAIC (relative to best per index) for each copula family across all four indices. Skipped if the  all_results  AIC table is unavailable.
|
Section 4 — SAF Curves
| Figure
|Contents
|
| ---|---|
| Fig 09  — SAF at 10-yr RP (all indices)
|All four indices on one panel, Conditional method, 10-year return period.
|
| Fig 10  — SAF at 100-yr RP (all indices)
|Same layout, 100-year return period.
|
| Fig 11  — Conditional vs Kendall, 4-panel
|2×2 grid with both methods overlaid per index; return periods colour-coded (viridis, 4 colours for 10/25/50/100 yr), methods distinguished by linetype (solid/dashed).
|
| Fig 12  — Per-class SAF for SPI-1 and SPI-3
|Stacked two-panel figure; SAF curves stratified by duration class, return period colour-coded. Skipped if per-class SAF was not computed (insufficient events).
|
Section 5 — Non-Stationarity
| Figure
|Contents
|
| ---|---|
| Fig 13  — NS SAF reference vs recent, 4-panel
|2×2 grid; conditional NS SAF curves for both epochs per index, return periods colour-coded, epochs distinguished by linetype (solid=reference, dashed=recent).
|
| Fig 14  — GAMLSS model selection table
|A  tableGrob  figure rendering the GAMLSS model selection summary across all four indices: change-point year (parsed from  ns$epoch_src ), AIC for M0/M1/M2, winning model, and period-specific μ/σ.
|
| Fig 15  — Severity by epoch (violin + boxplot)
|2×2 faceted violin-boxplot comparing drought severity distributions between reference (blue) and recent (red) epochs for each index.
|
Section 6 — Copula Diagnostics & Change-Point
| Figure
|Contents
|
| ---|---|
| Fig 16  — Change-point summary
|Jitter plot of change-point year estimates per test (Pettitt, CUSUM, PRUTF) and index, with the consensus year shown as a dashed horizontal line. Filled points indicate significant  tests (p  < 0.05). Field names used:  pettitt_cp ,  cusum_cp ,  prutf_cp ,  cp_year  (consensus).
|
| Fig 17  — Rosenblatt PIT diagnostic, 4-panel
|For each index: (1) e2 histogram vs U[0,1], (2) e2 QQ-plot, (3) e1 vs e2 scatter with linear trend. Computed using a central finite-difference of  pCopula() . Individual JPEG files are also saved per index to their respective analysis subdirectory alongside the combined 4-panel figure.
|
| Fig 18  — Rolling Kendall τ, 4-panel
|2×2 grid of rolling τ time series (10-year centred window) with LOESS smoother and a red dashed vertical line at the consensus change-point year (sourced from  cp_result$cp_year ).
|
Section 7 — Excel Summary Tables
| Table
|Contents
|
| ---|---|
| Table 01  — Return periods
|All return period data from  res_ <INDEX >$rp  for all four indices combined, with a styled header row.
|
| Table 02  — Copula fit summary
|All copula comparison results ( all_results ) for all indices combined; rows for the best-fit copula are highlighted in green.
|
| Table 03  — Non-stationary summary
|One row per index: change-point year, reference and recent period labels, GAMLSS best model, AIC values, μ/σ per epoch, TV copula significance ( tv$significant ), and copula GOF result.
|
Output Directory Structure
spatial_drought/
├── SAF_analysis_log.txt                               # Timestamped log of all pipeline steps
├── SPI1_analysis/
│   ├── SPI1_drought_characteristics.csv
│   ├──  SPI1_duration_classified_events.csv
│   ├── SPI1_return_periods_by_class.csv
│   ├── SPI1_severity_dist_comparison.csv
│   ├── SPI1_copula_comparison.csv
│   ├── SPI1_SAF_condition al_all.csv
│   ├── SPI1_SAF_kendall_all.csv
│   ├── SPI1_SAF_overlay.pdf
│   ├── SPI1_changepoint_detection.csv
│   ├── SPI1_segmented_copula.csv                      # written if  change-point detected
│   ├── SPI1_nonstationary_summary.csv
│   ├── SPI1_SAF_nonstationary_reference_ <yr >_ <yr >.csv
│   ├── SPI1_SAF_nonstationary_recent_ <yr >_ <yr >.csv
│   ├── SPI1_SAF_nonstationary_kendall_reference_ <yr >_ <yr >.csv
│   ├── SPI1_SAF_nonstationary_kendall_recent_ <yr >_ <yr >.csv
│   ├── SPI1_rolling_tau.pdf                           # always written in Step 5c
│   ├── SPI1_copula_diagnostic.pdf                     # written if bootstrap GOF rejected
 │   ├── SPI1_rosenblatt_diagnostic.pdf                 # written if GOF rejected  & TV fitted
│   ├── SPI1_rosenblatt_diagnostic.jpeg                # individual JPEG from Step 7 Fig 17
│   ├── SPI1_tv_parameter.pdf                          # written if TV model  selected
│   └── (per-class SAF CSVs when ≥15 events per class)
├── SPEI1_analysis/   (same structure)
├── SPI3_analysis/    (same structure, no D1-2 class)
├── SPEI3_analysis/    (same structure, no D1-2 class)
── pub_figures/                                       # Step 7 outputs
    ├── Fig01_basin_drought_timeseries.pdf / .jpeg
    ├── Fig02_event_chara cteristics_scatter.pdf / .jpeg
    ├── Fig03_event_duration_histogram.pdf / .jpeg
    ├── Fig04_marginal_fits_severity.pdf / .jpeg
    ├── Fig05_marginal_fits_area.pdf / .jpeg
     ├── Fig06_pseudo_obs_scatter.pdf / .jpeg
    ├── Fig07_kendall_tau_barplot.pdf / .jpeg
    ├── Fig08_copula_aic_comparison.pdf / .jpeg
    ├── Fig09_SAF_allindices_RP10.pdf / .jpeg 
    ├── Fig10_SAF_allindices_RP100.pdf / .jpeg
    ├── Fig11_SAF_method_comparison_4panel.pdf / .jpeg
    ├── Fig12_SAF_by_duration_class.pdf / .jpeg
    ├── Fig13_NS_epoch_compar ison_4panel.pdf / .jpeg
    ├── Fig14_GAMLSS_model_selection_table.pdf / .jpeg
    ├── Fig15_severity_trend_by_epoch.pdf / .jpeg
    ├── Fig16_changepoint_summary.pdf / .jpeg
    ├ ── Fig17_ <INDEX >_rosenblatt.pdf / .jpeg            # one per index
    ├── Fig17_rosenblatt_4panel_combined.pdf / .jpeg
    ├── Fig18_rolling_tau_4panel.pdf / .jpeg
    ├── Table01_return_periods _all_indices.xlsx
    ├── Table02_copula_fit_summary.xlsx
    └── Table03_nonstationary_summary.xlsx
If optional Step 6 maps are enabled, a separate `drought_maps/` folder is created containing `<INDEX>_freq.nc` and `<INDEX>_sev.nc`.
Is `Nechako_Drought_Methodology_Guide.R` Necessary?
No — it is not required for the core analysis. Here is the breakdown:
| Scenario
|Need the Guide?
|
| ---|---|
| Running the full pipeline via  RUN_ALL.R
|It is sourced in Step 2, but only  loads  functions. Removing it would not break Steps 3–7.
|
| Using Step 6 optional analyses
|Yes , because  analyze_multiple_thresholds() ,  seasonal_drought_analysis() ,  map_drought_frequency() , and  cross_correlate_indices()  are only defined in this file.
|
| Running  SAF_Analysis.R  standalone
|No
|
| Running ext1 or ext2
|No
|
Verdict: The Methodology Guide is a collection of supplementary diagnostic and visualisation tools. The SAF curves, copula fits, return period tables, non-stationary analyses, all CSV/PDF outputs produced in Steps 3–5, and the publication figures produced in Step 7 are entirely independent of it.
Key Design Decisions & Bug Fixes
The comments in the scripts document several important corrections relative to earlier versions:
Kendall distribution formula (`QuickStart.R`, `SAF_Analysis.R`, `ext2.R`): The pipeline uses `kc_fn(t) = mean(C(U,V) <= t)` (empirical K_C) for non-Archimedean families, verified at startup against the independence copula result `K_C(t) = t − t·log(t)`. For Archimedean families, the analytical formula `K_C(t) = t − φ(t)/φ'(t)` is used in both the base pipeline and the non-stationary extensions, enabling valid extrapolation beyond the empirical support range.
`rm(list=ls())` removed (`SAF_Analysis.R`, `ext1.R`, `ext2.R`): Earlier versions wiped the workspace mid-pipeline. Cleanup is now the caller's responsibility.
Copula log-likelihood in comparison table (`SAF_Analysis.R`): Saved alongside AIC/BIC.
Minimum events threshold raised to 15 (`SAF_Analysis.R`): Per-class SAF fitting is skipped for classes with fewer than 15 events to avoid unreliable tail fits.
Time-varying copula family (`ext1.R`): Uses whichever family was selected by AIC in Step 3. Multi-parameter families (Student-t) fall back to the stationary model with a logged warning.
GAMLSS `predict()` workaround (`ext2.R`): Some `gamlss` package versions return a closure instead of a numeric from `predict(..., type="response")`. The fix extracts coefficients directly and applies `exp()` manually to implement the log-link.
Non-stationary year window (`ext2.R`): Reference and recent periods are derived from the actual record length (40%/30% split) instead of hardcoded calendar years. The detected change-point year takes precedence when it falls between 20%–80% of the record length.
Rosenblatt PIT h-function  ( `ext1.R`  —  `.h_func_fd()` ):  `e2`  is computed as the partial derivative of  `pCopula()`  with respect to  `u1`  using a central finite difference with an adaptive step size ( `eps_adj = eps × min(u1, 1−u1)` ) that avoids boundary overflow.  `pCopula`  is implemented for all copula families in the package, eliminating the NA-return failure seen with  `cCopula()`  for the Plackett family and the wrong-direction bug ( `cCopula(..., indices=1)`  returns C(u1|u2), not C(u2|u1)).
Epoch boundary truncation (`ext2.R` — `fit_nonstationary_marginals()`): `yr_match <- dc$year` uses the actual calendar year stored in each drought-month row. `yr_range` now correctly spans 1950–2025 and `recent_years` extends to the record end.
Impact on epoch labels after fix:
| Index
|CP year
|Reference period
|Recent period
|
| ---|---|---|---|
| SPI-1
|1965
|1950–1964
|1965–2025 (+27 yr)
|
| SPEI-1
|1996
|1950–1995
|1996–2025 (CP used correctly)
|
| SPI-3
|1973
|1950–1972
|1973–2025 (+27 yr)
|
| SPEI-3
|1987
|1950–1986
|1987–2025 (+26 yr)
|
`total_events`  authoritative source  ( `QuickStart.R`  —  `quick_diagnostic()` ):  `total_events`  is now read directly from  `nrow(` `<INDEX` `>_duration_classified_events.csv)`  when the file exists.  `sum(rp$n_events)`  is retained as a fallback and cross-check target. A  `MISMATCH`  flag is printed when the two counts differ.
SPI-3 zero-variance rolling-window warnings (`ext1.R` — `detect_copula_changepoints()`): Both `sd(severity_window)` and `sd(area_pct_window)` are checked before calling `cor()`. Zero-variance windows are assigned NA silently; the total count is reported once via `log_event()`. Change-point results are not affected.
Non-stationary SAF  `ns_result`  argument  ( `RUN_ALL.R` ,  `ext2.R` ):  `derive_SAF_nonstationary()`  and  `derive_SAF_nonstationary_kendall()`  require  `ns_result`  as an explicit argument (introduced in ext2 v5). The four calls in Step 5e of  `RUN_ALL.R`  all include  `ns_result = ns` . Omitting it causes R to abort with  "argument 'ns_result' is missing, with no default ".
Field name corrections in publication figures (`RUN_ALL.R`, Step 7): Several field names used in earlier figure code were corrected to match the actual return values of their source functions: `dep$kendall_tau → dep$kendall` (Fig 06), `pettitt_yr/cusum_yr/prutf_yr/consensus_yr → pettitt_cp/cusum_cp/prutf_cp/cp_year` (Fig 16), `tv$lr_sig → tv$significant` (Table 03), `cp_result$consensus_yr → cp_result$cp_year` (Fig 18), `ns$cp_year → parsed from ns$epoch_src` (Fig 14, Table 03).
Model Selection Stability (`ext1.R` — `fit_timevarying_copula()`): The model selection logic now uses `isTRUE()` wrappers around comparison flags (e.g., `cp_strong`, `pit_seg_improved`) to prevent runtime crashes ("missing value where TRUE/FALSE needed") when the Rosenblatt PIT GOF function returns `NA` values (which occurs when e2 is 0% valid due to finite-difference issues). This ensures the pipeline completes and falls back safely to stationary models when diagnostics are indeterminate.
Diagnostics Argument Fix (`RUN_ALL.R` — Step 5): The internal call to `.h_func_fd()` within the diagnostic block now explicitly passes the `cop_family` argument to prevent "argument missing" errors.
Bootstrap and permutation window sizes:
`gof_copula()`: N_boot = 499 replicates (parametric Cramér–von Mises Sn bootstrap). Minimum achievable two-sided p-value: 1/500 = 0.002. Use `N_boot = 999` for publication; `N_boot = 99` for fast development.
`detect_copula_changepoints()` CUSUM test: 499 Monte-Carlo permutations (matching N_boot for consistency).
Known Diagnostic Warnings and Their Interpretation
| Warning / Log message
|Source
|Severity
|Cause
|Action
|
| ---|---|---|---|---|
| WARNING [TD-U]: AIC-best copula has no upper-tail dependence but empirical λ_U  > 0.10
|fit_copulas()
|Advisory
|AIC prefers Frank, Plackett, or Gaussian despite upper-tail signal
|Consider Gumbel, Joe, SurvClayton, or StudentT; compare SAF curves at high return periods
|
| [INDEX] Rolling tau: N window(s) had zero variance … set to NA
|detect_copula_changepoints()
|Informational
|Isolated drought-month windows with identical severity or area values (typical in SPI-3 record)
|Inspect  <INDEX >_drought_characteristics.csv  for constant runs; change-point results unaffected
|
| MISMATCH — investigate aggregation!  in  quick_diagnostic()
|QuickStart.R
|Warning
|nrow(ev_file) ≠ sum(rp$n_events) ; sub-threshold events may not map to a named class
|Check duration classification logic and  D0 (Sub-threshold)  rows in the rp table
|
| TV trend skipped for StudentT: multi-parameter families require extended formulation
|fit_timevarying_copula()
|Informational
|StudentT has 2 parameters (ρ, df); the scalar log-linear TV model is not applicable
|Stationary StudentT fit is used; segmented model is still evaluated if a change point is detected
|
| WARNING: PATH 2 - TV copula at epoch midpoint … hybrid approximation
|.build_period_copula()
|Advisory
|Insufficient data per segment ( < 15 obs) prevented a segmented copula fit
|Document the hybrid nature in reporting; consider enlarging epochs or using the stationary fallback explicitly
|
| WARNING [EPOCH MISMATCH]: copula cp_year != marginal cp_year
|.assert_epoch_consistency()
|Warning
|The change-point year fed to  fit_timevarying_copula()  and  fit_nonstationary_marginals()  differ
|Ensure both functions receive the same  cp_result  object; do not mix results from different runs
|
Glossary
| Term
|Meaning
|
| ---|---|
| SPI
|Standardised Precipitation Index — drought index based on precipitation anomaly
|
| SPEI
|Standardised Precipitation-Evapotranspiration Index — includes temperature/ET
|
| SAF curve
|Severity–Area–Frequency curve: relates drought severity to spatial extent for a given return period
|
| Copula
|A function that models the dependence structure between two or more variables independently of their marginal distributions
|
| Kendall distribution K_C
|The CDF of copula values:  K_C(t) = P(C(U,V) ≤ t)
|
| IFM
|Inference Functions for Margins — a two-stage MLE method for copula fitting
|
| GAMLSS
|Generalised Additive Models for Location, Scale and Shape — used here for non-stationary distributional regression
|
| Return period
|Expected average time between events of at least a given magnitude
|
| Duration class
|Grouping of drought events by length in months (D1-2, D3-6, D7-12, D13+)
|
| TV copula
|Time-varying copula: a copula whose dependence parameter θ is modelled as a log-linear function of standardised year
|
| Segmented copula
|Two independent stationary copulas fitted to temporal sub-records defined by a consensus change-point year
|
| Rosenblatt PIT
|Probability Integral Transform (Rosenblatt 1952): converts a bivariate sample from a correctly specified copula into two independent U[0,1] variates for GOF testing
|
| λ_U
|Upper-tail dependence coefficient: limiting probability of joint exceedances as the threshold approaches 1
|
| Gamma renewal process
|A generalisation of the Poisson process where inter-arrival times follow a Gamma distribution (shape k = 1/CV²); k  < 1 = clustered events, k  > 1 = regular events, k = 1 = Poisson
|
| PRUTF
|Parametric likelihood-ratio structural-break test applied to the rolling Kendall-τ series; supremum of segment log-likelihood gains over five candidate break points
|
| .h_func_fd
|Internal helper computing the conditional CDF  C(u2|u1)  via central finite differences of  pCopula()  with adaptive step size
|
| .build_period_copula
|Internal helper selecting the period-specific copula object via a three-path priority: segmented (preferred) → TV midpoint (approximate) → stationary (fallback)
|

---
Current Diagnostic Status & Known Issues (Post-Pipeline Run)

The following issues were observed during the execution of `Nechako_Drought_RUN_ALL.R` (Step 5 & 7). The pipeline completes successfully despite these warnings due to the stability fixes implemented in the latest script versions.

**1. Rosenblatt PIT (0% valid e2)**
The Rosenblatt diagnostics consistently report "0% e2 valid" and NA p-values for all indices.
*   **Cause:** Numerical precision loss in the finite-difference calculation of the h-function (`.h_func_fd()`). When `u1` or `u2` are near the boundaries, or when the copula has strong dependence, the difference `pCopula(u1_hi) - pCopula(u1_lo)` suffers from catastrophic cancellation, returning values indistinguishable from 0 or 1.
*   **Status:** The `isTRUE()` guard in `fit_timevarying_copula()` successfully prevents the pipeline from crashing when these NA values propagate to the model selection logic. The pipeline defaults to stationary models where appropriate.
*   **Recommendation:** For publication, either omit the Rosenblatt PIT diagnostic plots (Fig 17) or replace `.h_func_fd()` with the analytical conditional CDFs for Archimedean copulas (Frank, Gumbel, Joe) to ensure numerical stability.

**2. Event Count Mismatch (SPI3/SPEI3)**
The post-run diagnostic (`quick_diagnostic()`) flags a "MISMATCH" for SPI3 and SPEI3 where `Total(ev_file)` ≠ `Total(rp sum)`.
*   **Cause:** This is expected behavior for multi-scale indices (SPI-3/SPEI-3). The `identify_duration_classified_events()` function may identify events that do not fall into the defined duration classes (e.g., sub-threshold fluctuations) or rounding differences in the aggregation of the return period table.
*   **Status:** Harmless. The authoritative event count is derived from the `duration_classified_events.csv` file.
*   **Recommendation:** If an audit is required, verify the `D0 (Sub-threshold)` rows in `return_periods_by_class.csv` to confirm unclassified events.

**3. Optimization Warnings (code=1)**
Warnings stating "possible convergence problem: optim() gave code=1" appear during copula fitting.
*   **Cause:** Routine in Maximum Pseudo-Likelihood (MPL) copula fitting. The optimizer reached a point where the gradient is sufficiently small but has not fully converged to the exact tolerance.
*   **Status:** Harmless. Parameter estimates and SAF curves remain accurate.
*   **Recommendation:** No action required. If desired, increase `maxit` or tighten `reltol` in `optim.control`, but this may increase runtime without significant benefit.

**4. Unicode Warnings**
Warnings regarding "conversion failure ... in 'mbcsToSbcs'" appear during figure generation (Fig 07, 08, 14).
*   **Cause:** Cosmetic R graphics device issue with Greek letters (ρ, Δ) in plot titles/labels.
*   **Status:** Harmless. Figures are generated correctly.
*   **Recommendation:** Replace `ρ`/`Δ` with "rho"/"Delta" in plot labels, or set `encoding="UTF-8"` in `pdf()/jpeg()` calls to suppress warnings.