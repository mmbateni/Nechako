# Nechako Drought SAF Analysis — README

## Overview

This codebase performs a **Severity–Area–Frequency (SAF) drought analysis** for the Nechako Basin in British Columbia, Canada. It uses gridded **SPI** (Standardised Precipitation Index) and **SPEI** (Standardised Precipitation-Evapotranspiration Index) raster time series as inputs, and produces statistically rigorous SAF curves, copula models, return period tables, diagnostic plots, publication-quality figures, and Excel summary tables as outputs.

The analysis follows a multivariate frequency analysis framework: drought severity and drought-affected area fraction are modelled as two correlated random variables, joined by a copula. SAF curves answer the question: *"For a given return period, how severe is a drought that affects X% of the basin?"*

---

## File Inventory & Roles

| File | Role | Required? |
|---|---|---|
| `Nechako_Drought_RUN_ALL.R` | **Master orchestrator** — single entry point that runs all seven steps in order | Yes |
| `Nechako_Drought_QuickStart.R` | Input verification and pre-flight checks | Yes |
| `Nechako_Drought_SAF_Analysis.R` | Core analysis pipeline (data loading → SAF curves) | Yes |
| `Nechako_Drought_SAF_Analysis_ext1.R` | Extension v6: dependency tests, copula GOF, time-varying and segmented copula, diagnostic plots | Yes |
| `Nechako_Drought_SAF_Analysis_ext2.R` | Extension v7: non-stationary marginals, Kendall-corrected NS SAF, epoch helpers, event placement | Yes |
| `Nechako_Drought_Methodology_Guide.R` | Optional utility functions (threshold sensitivity, seasonal plots, spatial maps) | **Optional** |

---

## Prerequisites

### Working Directory

All scripts assume:

```r
setwd("D:/Nechako_Drought/Nechako")
```

Change this path before running if your data lives elsewhere.

### Required Input Files

The pipeline reads monthly NetCDF files pre-computed by an upstream SPI/SPEI script. They must exist at:

```
spi_results_seasonal/spi_01_month01_Jan.nc   # SPI-1 for January
spi_results_seasonal/spi_01_month02_Feb.nc   # SPI-1 for February
...                                           # (12 files per scale)
spi_results_seasonal/spi_03_month01_Jan.nc   # SPI-3 for January
...
spei_results_seasonal/spei_01_month01_Jan.nc # SPEI-1 for January
...
spei_results_seasonal/spei_03_month01_Jan.nc # SPEI-3 for January
...
```

Each NetCDF contains a multi-layer raster (one layer per year) for that calendar month and scale.

### R Packages

The scripts auto-install missing packages. The full list is:

| Package | Purpose |
|---|---|
| `terra` | Raster I/O and spatial operations |
| `copula` | Copula fitting, GOF, density/CDF functions |
| `fitdistrplus` | MLE fitting of univariate marginals |
| `MASS` | Supporting distributions |
| `ggplot2` | All plots |
| `gridExtra` | Multi-panel PDF plots and table grobs |
| `viridis` | Colour scales |
| `dplyr` | Data manipulation |
| `scales` | Axis formatting |
| `moments` | Higher-order statistics (ext1) |
| `openxlsx` | Excel export (ext1, Step 7) |
| `Kendall` | Rank correlation tests (ext1) |
| `gamlss` | Non-stationary GAMLSS models (ext2) |
| `gamlss.dist` | GAMLSS distribution families (ext2) |
| `trend` | Pettitt test (ext1) |
| `patchwork` | Multi-panel ggplot2 layout (Step 7) |
| `ggExtra` | Marginal histogram/density plots (Step 7) |

---

## How to Run

The simplest approach is the single entry point:

```r
setwd("D:/Nechako_Drought/Nechako")
source("Nechako_Drought_RUN_ALL.R")
```

This runs all seven steps (see below) in the correct order. Do **not** source the individual analysis scripts out of order — `ext1.R` and `ext2.R` depend on objects created by `SAF_Analysis.R`, and `SAF_Analysis.R` depends on objects created by `QuickStart.R`.

---

## Pipeline Steps (as executed by `RUN_ALL.R`)

### Step 1 — Input Verification (`QuickStart.R`)

- Spot-checks **4 representative NetCDF files** (one per index, January only) as proxies for all 48 required files (4 indices × 12 months). A missing file stops execution; the remaining 44 files are verified implicitly when the pipeline loads them in Step 3.
- Loads `SPI-1` for January and validates that pixel values fall within the physically plausible range of −5 to +5.
- Runs a mathematical sanity check on the Kendall distribution formula using the known analytical result for the independence copula at `t = 0.5`.
- Prints the recommended run order to the console.

**`quick_diagnostic()` patch (FIX-5):** The version of `quick_diagnostic()` defined in `QuickStart.R` hardcodes the output search path as `"drought_analysis"`, but the pipeline writes all results to `"spatial_drought"`. `RUN_ALL.R` overrides `quick_diagnostic()` inline immediately after `source("Nechako_Drought_QuickStart.R")` to search `spatial_drought` instead. The patched version also sources event counts from the raw event CSV (`nrow(<INDEX>_duration_classified_events.csv)`) as the **authoritative** count, with `sum(rp$n_events)` from the return-period table retained as a cross-check. A `MISMATCH` flag is printed when the two counts differ.

### Step 2 — Methodology Helpers (`Methodology_Guide.R`)

Loads five utility functions into the workspace. **No computation is triggered at load time.** These functions are available for optional Step 6 calls:

- `analyze_multiple_thresholds()` — repeats the area fraction computation at drought thresholds −0.5, −1.0, and −1.5 and produces a PDF line chart.
- `seasonal_drought_analysis()` — computes mean severity and mean area fraction by calendar month and saves a two-panel bar chart PDF.
- `map_drought_frequency()` — computes per-pixel drought frequency (% of months below threshold) and mean per-pixel severity, saving both as NetCDF files.
- `compare_return_period_methods()` — plots Conditional vs Kendall SAF curves on a single overlay for visual method comparison.
- `cross_correlate_indices()` — scatter-plots SPI-1 vs SPI-3 event severity and reports the Pearson correlation.

### Step 3 — Main SAF Pipeline (`SAF_Analysis.R`)

This is the computational core. It runs the full pipeline four times — once each for **SPI-1**, **SPEI-1**, **SPI-3**, and **SPEI-3** — and writes results to `spatial_drought/<INDEX>_analysis/`. Results are stored in the workspace objects `res_spi1`, `res_spei1`, `res_spi3`, `res_spei3`.

#### 3a. Data Loading

All 12 monthly NetCDF tiles for each index are loaded with `terra::rast()` and concatenated into a single time-ordered raster stack. Cell areas (km²) are computed geodesically with `terra::cellSize(..., unit = "km")`, which integrates the latitude band directly for geographic CRS and uses projected units for metric CRS — no equal-area reprojection is required. Basin total area is computed once from the SPI-1 grid.

#### 3b. Drought Characteristics Extraction

For every time step `t`, the script:
1. Reads all pixel values from the index raster layer.
2. Counts pixels below the drought threshold (default −0.5) as the **drought area** (expressed as % of valid pixels).
3. Computes **mean drought severity** as the mean exceedance below threshold: `mean(threshold − value)` for all dry pixels.

Output: a time series data frame `[date, year, month, severity, area_pct]` saved as `<INDEX>_drought_characteristics.csv`.

#### 3c. Event Identification & Duration Classification

Consecutive months with `severity > 0` and `area_pct > 0` are grouped into discrete drought **events** using a run-length algorithm. Each event records start/end date, duration in months, duration class, mean/max severity, and mean/max area fraction.

Duration classes depend on the index scale:

| Class | SPI-1 / SPEI-1 | SPI-3 / SPEI-3 |
|---|---|---|
| D1-2 (Very-short) | 1–2 months | — |
| D3-6 (Short-term) | 3–6 months | 3–6 months |
| D7-12 (Medium-term) | 7–12 months | 7–12 months |
| D13+ (Long-term) | ≥13 months | ≥13 months |

**D0 sub-threshold filter (current version):** For scale-3 indices (SPI-3, SPEI-3), any event with a duration of fewer than 3 months receives the sentinel class `"D0 (Sub-threshold)"` from `classify_duration_class()`. These events are **removed entirely** before the function returns. This ensures the CSV files, return-period tables, copula fits, and SAF curves for 3-month indices are never contaminated by events too short to be physically meaningful for that accumulation scale.

Empirical **return periods** are computed as `record_length_years / n_events` per class, with exact Poisson 95% confidence intervals on the rate via the chi-squared relationship:
```
lambda_lo = qchisq(0.025, df = 2n)     / (2 × record_years)
lambda_hi = qchisq(0.975, df = 2(n+1)) / (2 × record_years)
```

**Gamma Renewal Adjustment:** Inter-arrival time statistics (mean, CV, Gamma shape `k = 1/CV²`) are computed from each event class. When `CV ≠ 1`, the process deviates from a homogeneous Poisson assumption — `k < 1` indicates clustering (heavy-tailed gaps) and `k > 1` indicates regularity. Adjusted mean inter-arrival times `mu_T_adj` and Gamma parameters are stored in the return-period table and passed to the SAF target computation (see §3f).

Outputs: `<INDEX>_duration_classified_events.csv`, `<INDEX>_return_periods_by_class.csv`.

#### 3d. Marginal Distribution Fitting

Using only months with active drought (`severity > 0`, `area_pct > 0`):

- **Area fraction** is fitted with a **Beta distribution** (values bounded 0–1).
- **Severity** is fitted with six distributions (Exponential, Gamma, Weibull, Log-Normal, Normal, Logistic) using MLE. The best-fitting distribution is selected by minimum AIC.

Output: `<INDEX>_severity_dist_comparison.csv`.

#### 3e. Copula Fitting

The parametric marginal CDFs transform severity and area into uniform pseudo-observations `(u_S, u_A)`. **Seven** copula families are fitted using the IFM (Inference Functions for Margins) maximum pseudo-likelihood method:

| Copula | Tail dependence |
|---|---|
| Survival Clayton | Upper tail |
| Gumbel | Upper tail |
| Frank | Symmetric (no tail) |
| Joe | Upper tail |
| Plackett | Symmetric (no tail) |
| Gaussian (Normal) | No tail (asymptotic independence for all \|ρ\| < 1) |
| Student-t | Both tails (symmetric) |

**Tail-dependence priority logic:** Before accepting the AIC-best family, an empirical upper-tail dependence coefficient λ_U is estimated non-parametrically using the **CFG–Pickands estimator**:

```
A(t) = exp(−E[log(min(−log(U)/(1−t), −log(V)/t))])
λ_U  = max(0, 2(1 − A(0.5)))
```

If `λ_U > 0.5` (strong upper-tail signal) *and* the best Joe/Gumbel AIC is within 4 units of the global best (Burnham & Anderson criterion), the selection is redirected to the Joe or Gumbel family (whichever has lower AIC). If the AIC-best family has zero-tail dependence (Frank, Plackett, Gaussian) but `λ_U > 0.10`, a `TD-MISMATCH OVERRIDE` is logged and the best tail-dependent family (StudentT, Gumbel, Joe) is substituted.

AIC, BIC, log-likelihood, empirical λ_U, and theoretical λ_U are all saved. A `WARNING [TD-U]` is logged if the final chosen family still has zero upper-tail dependence despite `λ_U > 0.10`.

Output: `<INDEX>_copula_comparison.csv`.

#### 3f. SAF Curve Derivation — Two Methods

**Clustering-Adjusted SAF Target (`.compute_renewal_saf_target()`):**

Both SAF methods use a shared helper to compute the monthly non-exceedance target. When the event process is approximately Poisson (CV ∈ [0.5, 1.5] or CV unavailable), the standard formula is used:

```
target = 1 − mu_T_adj / (T × 12)
```

When CV falls outside this range, a Gamma renewal survival function is used instead:

```
p_no_event_T = 1 − pgamma(T, shape = k, rate = k / mu_T_adj_yrs)
target = 1 − (1 − p_no_event_T) / (T × 12)
```

This corrects the SAF exceedance probability for clustered (k < 1) or regular (k > 1) event processes.

**Method A: Conditional Copula**

For a target return period `T` years, the effective monthly non-exceedance probability is computed via the renewal helper above. For each area value `a`, the script solves numerically for the severity `s` such that:

```
P(S > s | A = a) = 1 − target
```

using the conditional copula `∂C(u_S, u_A)/∂u_A` (implemented as `cCopula(..., indices = 2)` with a tCopula df integer fix — see §Key Design Decisions, item 1). Inter-arrival time statistics are sourced from the authoritative event table where available (Path A), falling back to a lag-detection method on the full drought characteristics series (Path B).

**Method B: Kendall Distribution**

Instead of conditioning on area, this method uses the joint level set of the copula. It solves for the copula level `t_val` such that:

```
K_C(t_val) = target
```

For **Archimedean families** (Frank, Gumbel, Joe, SurvClayton), the **analytical** Kendall distribution `K_C(t) = t − φ(t)/φ'(t)` is used, enabling stable extrapolation beyond the empirical support. For non-Archimedean families (Gaussian, Student-t, Plackett), the **empirical** `K_C(t) = mean(C(U, V) ≤ t)` is used.

Then for each area value `a`, it finds severity `s` such that `C(u_S, u_A) = t_val`.

Both methods produce curves for return periods of **10, 25, 50, and 100 years** across area fractions from 5% to 95%.

SAF curves are also computed **per duration class** when the class has ≥ 15 events (threshold raised from 10 in earlier versions).

Outputs per index: `<INDEX>_SAF_conditional_all.csv`, `<INDEX>_SAF_kendall_all.csv`, `<INDEX>_SAF_overlay.pdf`, and per-class equivalents.

### Step 4 — Load Extensions (`ext1.R` v6 + `ext2.R` v7)

The extension scripts are sourced here. They only define functions — no heavy computation occurs at load time. Both scripts verify that `log_event()` from `SAF_Analysis.R` is present before loading; an error stops execution if the main pipeline was skipped.

> **FIX-8:** The version label printed during loading was corrected from `"ext2 v4"` to `"ext2 v5"` (now v7).

**Functions exported by `ext1.R` (v6):**
`test_dependency`, `gof_copula`, `detect_copula_changepoints`, `rosenblatt_pit_gof`, `fit_timevarying_copula`, `fit_segmented_copula`, `plot_copula_diagnostic`, `plot_rosenblatt_diagnostic`, `plot_rolling_tau`, `plot_tv_parameter`

**Functions exported by `ext2.R` (v7):**
`fit_nonstationary_marginals`, `derive_SAF_nonstationary`, `derive_SAF_nonstationary_kendall`, `place_event_on_saf`
Internal helpers (not exported directly): `.build_period_copula` (revised), `.resolve_epoch` (new), `.assert_epoch_consistency` (new)

### Step 5 — Run Enhancements

The master loop iterates over the four `INDEX_RESULTS` entries (`SPI1`, `SPEI1`, `SPI3`, `SPEI3`). For each index, three enhancement analyses are run and two conditional diagnostics are applied. Results accumulate in the `enh_results` list for use in Step 7.

> **Note:** Non-stationary marginal fitting (ext2 `fit_nonstationary_marginals`) and non-stationary SAF derivation (ext2 `derive_SAF_nonstationary` / `derive_SAF_nonstationary_kendall`) are **not** called in the Step 5 loop of `RUN_ALL.R`. Those functions must be invoked separately and their results stored in `enh_results[[nm]]$ns_marginals`, `enh_results[[nm]]$saf_ns_ref`, and `enh_results[[nm]]$saf_ns_rec` before Step 7 figures 13–15 can be generated.

#### 5a. Dependency Test

Computes **Kendall's τ** and **Spearman's ρ** between drought severity and area fraction, with a significance test. This validates whether the copula is statistically justified (if τ ≈ 0, the variables are nearly independent and a copula provides little benefit).

Result field access: `enh_results[[nm]]$dependency$kendall`, `$spearman`, `$p_value`.

> **FIX-1E:** Figure 6 (Step 7) was corrected to reference `dep$kendall` rather than the former `dep$kendall_tau` field name.

#### 5b. Copula Goodness-of-Fit (Bootstrap Sn)

Performs a parametric bootstrap GOF test using the Cramér–von Mises `Sn` statistic with **499 bootstrap replicates** (`N_boot = 499`). The minimum achievable two-sided p-value is 1/500 = 0.002, which is below the α = 0.05 decision threshold. Use `N_boot = 999` for publication-quality p-values; `N_boot = 99` for fast development.

A p-value > 0.05 means the chosen copula family is not rejected at the 5% level. Decision is stored as `"ACCEPTED"` or `"REJECTED"`.

Result field: `enh_results[[nm]]$gof`.

#### 5c. Time-Varying & Segmented Copula

Tests whether the dependence structure between severity and area has changed over time. The full pipeline:

**Change-point detection** (Pettitt, CUSUM, PRUTF) is run first on a rolling Kendall-τ series. The **zero-variance guard** (introduced in ext1 v6) pre-checks `sd(severity_window)` and `sd(area_pct_window)` before calling `cor()`. Windows with zero variance in either variable are assigned `NA` silently; the total suppressed count is logged once via `log_event()`.

The three tests use **499 CUSUM Monte-Carlo permutations** (matching `N_boot` for consistency). A break is flagged when at least 2 of the 3 tests reject H₀ at the 5% level; the consensus change-point year is the median of the significant test estimates.

**Return value field names** (current version):

| Field | Meaning |
|---|---|
| `cp_result$pettitt_cp` | Pettitt change-point year |
| `cp_result$cusum_cp` | CUSUM change-point year |
| `cp_result$prutf_cp` | PRUTF change-point year |
| `cp_result$cp_year` | Consensus year |
| `cp_result$detected` | Logical; TRUE if ≥2 tests significant |

> **FIX-1B:** Figure 16 (Step 7) field name references corrected from `pettitt_yr` / `cusum_yr` / `prutf_yr` / `consensus_yr` to `pettitt_cp` / `cusum_cp` / `prutf_cp` / `cp_year`.
> **FIX-1D:** Figure 18 (Step 7) corrected from `cp_result$consensus_yr` to `cp_result$cp_year`.

**Three copula models** are fitted and compared:

| Model | Description | Free parameters |
|---|---|---|
| Stationary | Single fixed copula parameter | 1 |
| Time-varying (TV) | Log-linear trend: `θ(t) = ilink(a + b × year_std)` | 2 |
| Segmented | Two independent stationary copulas, one per epoch defined by the change-point year | 2 |

Model selection: TV is preferred when LRT is significant (p < 0.05), no strong abrupt break is detected, TV has lower AIC, and Rosenblatt PIT passes. Segmented is preferred when a change point is detected, segmented AIC is lower, and Rosenblatt diagnostic improves. Otherwise the stationary model is retained.

The `fit_timevarying_copula()` return object now includes extra fields for downstream use in ext2:

| Field | Purpose |
|---|---|
| `yr_mean`, `yr_sd` | Re-standardise a calendar year outside this function |
| `a_hat`, `b_hat` | TV copula parameters on the link scale |
| `cop_name` | Family name (for `make_cop_fn` reconstruction) |
| `link_hi`, `link_lo` | Parameter bounds |
| `link_ilink` | Inverse-link function |
| `make_cop_fn` | Closure: takes `theta_raw`, returns copula object |
| `cp_result` | Output of `detect_copula_changepoints()` |
| `rosenblatt` | Output of `rosenblatt_pit_gof()` |
| `significant` | Logical; TRUE if TV LRT p < 0.05 |

> **FIX-1C:** Table 3 (Step 7) TV-copula significance field corrected from `tv$lr_sig` to `tv$significant`.

**Multi-parameter family caveat:** The Student-t copula has two parameters (ρ, df). Because the log-linear TV formulation is scalar, the TV trend model is skipped for Student-t and the stationary fit is used instead (a warning is logged). The segmented model is still evaluated if a change point is detected.

**Post-fitting Rosenblatt PIT Diagnostics (inline in Step 5 loop):**

Two diagnostic plots are produced conditionally, using **event-level pseudo-observations** rather than the full monthly series. Monthly snapshots are intra-event clustered; using event-level data restores the i.i.d. assumption for the PIT test:

1. **Copula contour diagnostic** (`plot_copula_diagnostic()`): produced only when the bootstrap GOF was `"REJECTED"`.
2. **Rosenblatt PIT diagnostic** (`plot_rosenblatt_diagnostic()`): produced for every index with ≥5 identified events. Computes the h-function `e2 = C(u_A | u_S)` via ``.h_func_fd()` on event-level ranks, then saves a 3-panel plot (e2 histogram vs Uniform[0,1], QQ-plot, e1 vs e2 scatter) to the index analysis folder.

Result field: `enh_results[[nm]]$tv_copula`.

#### 5d & 5e. Non-Stationary Marginals and SAF Curves (ext2 — run separately)

The ext2 functions `fit_nonstationary_marginals()`, `derive_SAF_nonstationary()`, and `derive_SAF_nonstationary_kendall()` are **not called in the Step 5 loop**. They must be invoked manually and results stored in `enh_results` before Step 7 figures 13–15 can be generated.

When called, the workflow is:

**5d — `fit_nonstationary_marginals()`:** Fits three GAMLSS Gamma models to drought severity and selects by AIC:

| Model | mu formula | sigma formula |
|---|---|---|
| M0 (stationary) | `~ 1` | `~ 1` |
| M1 (mu-trend) | `~ year_std` | `~ 1` |
| M2 (mu+sigma-trend) | `~ year_std` | `~ year_std` |

Epoch definition (precedence):
1. Explicit `ref_years` / `recent_years` argument override.
2. `cp_result$detected == TRUE`: reference = record start … cp_year − 1; recent = cp_year … record end. The change-point is used only when it falls between 20%–80% of the record; otherwise the 40%/30% auto-split is applied with a warning.
3. Default automatic split: first 40% / last 30% of record.

The `epoch_src` string in the result object encodes which path was taken (e.g., `"change-point (1965)"`, `"auto 40/30 (CP edge)"`, `"user override"`). This string is parsed by `.resolve_epoch()` in ext2 and by `RUN_ALL.R` Step 7 figure generation.

> **FIX-1F:** `cp_year` is not a slot in the `fit_nonstationary_marginals()` return value. Step 7 figures (Fig 14, Table 3) now extract it by parsing `ns$epoch_src` with a regex for a 4-digit year, consistent with `.resolve_epoch()`.

Output: `<INDEX>_nonstationary_summary.csv`.

**5e — `derive_SAF_nonstationary()` / `derive_SAF_nonstationary_kendall()`:** Now require `ns_result` as an **explicit required argument** (breaking change from v4). Both functions call `.resolve_epoch(ns_result, period_label)` to obtain the epoch year vectors as the single source of truth, guaranteeing the copula and marginal always share the same epoch boundaries.

Copula selection follows a **three-path priority rule** via `.build_period_copula()`:

| Path | Condition | Description |
|---|---|---|
| Path 1 (preferred) | `tv_res$seg_result` is not NULL | Use segment-specific stationary copula from `fit_segmented_copula()`. Epoch-consistent — both marginal and copula are step-functions sharing the same `cp_year`. |
| Path 2 (hybrid, warned) | `seg_result` is NULL but `make_cop_fn` is available | Evaluate TV copula at epoch midpoint year. Logs a `WARNING`: mixes step-function marginal with point-evaluated smooth trend. |
| Path 3 (fallback) | All else fails | Use full-record stationary copula from `fit_copulas()`. |

The path taken is written to the `copula_path` column in every row of the output CSV.

**Epoch-consistency check:** `.assert_epoch_consistency()` compares the change-point year from `tv_copula_result$cp_result$cp_year` against the one encoded in `ns_result$epoch_src`. A `WARNING [EPOCH MISMATCH]` is logged if they differ.

**Analytical K_C for non-stationary Kendall SAF (v7):** When the copula family is Archimedean (Frank, Gumbel, Joe, SurvClayton), the analytical generator formula `K_C(t) = t − φ(t)/φ'(t)` is used. Non-Archimedean families fall back to an empirical estimator via a **strict 3-tier validation**:
- Tier 1: Direct `pCopula` on fitted pseudo-observations.
- Tier 2: Simulation-based K_C with N = 20,000 samples.
- Tier 3: Hard abort (returns `NULL`) — prevents miscalibrated SAF curves from propagating silently.

Outputs: `<INDEX>_SAF_nonstationary_<period>.csv`, `<INDEX>_SAF_nonstationary_kendall_<period>.csv`.

---

### Step 6 — Optional Utility Analyses (commented out in `RUN_ALL.R`)

These lines are commented out but ready to activate:

| Comment block | What it does |
|---|---|
| `6a` | Threshold sensitivity: re-runs area fraction at −0.5, −1.0, −1.5 |
| `6b` | Seasonal bar charts of severity and area by calendar month |
| `6c` | Pixel-level frequency/severity maps saved as NetCDF |
| `6d` | Cross-index correlation: SPI-1 vs SPI-3 severity scatter plot |
| `6e` | Places a specific historical event (e.g. years 2003–2004) on the SAF surface and returns its joint return period |

---

### Step 7 — Publication Figures & Excel Tables (NEW)

Generates all publication-ready outputs to `spatial_drought/pub_figures/`. Each figure is saved as both a PDF and a JPEG (300 dpi). A shared theme (`pub_theme()`) and a helper (`save_pub()`) are defined at the top of this step.

**Index colour palette:** SPI1 = `#1b7837`, SPEI1 = `#762a83`, SPI3 = `#2166ac`, SPEI3 = `#d6604d`.

#### Section 1 — Study Area & Data Overview

| Figure | File | Contents |
|---|---|---|
| Fig 01 | `Fig01_basin_drought_timeseries` | Monthly drought severity (top) and area % (bottom) for all 4 indices, 1950–2025. Two-panel patchwork with shared legend. |
| Fig 02 | `Fig02_event_characteristics_scatter` | Event scatter: mean area % vs mean severity, coloured by duration class, shaped by index. |
| Fig 03 | `Fig03_event_duration_histogram` | Faceted histogram of event duration in months, one panel per index. |

#### Section 2 — Marginal Distributions

| Figure | File | Contents |
|---|---|---|
| Fig 04 | `Fig04_marginal_fits_severity` | 4-panel: empirical histogram vs best-fit severity marginal density for each index. |
| Fig 05 | `Fig05_marginal_fits_area` | 4-panel: empirical histogram vs fitted Beta density for drought-affected area %. |

#### Section 3 — Copula Dependence Structure

| Figure | File | Contents |
|---|---|---|
| Fig 06 | `Fig06_pseudo_obs_scatter` | 4-panel: pseudo-observation scatter with fitted copula contours, annotated with best copula name and Kendall τ. |
| Fig 07 | `Fig07_kendall_tau_barplot` | Bar chart of Kendall τ and Spearman ρ by index. |
| Fig 08 | `Fig08_copula_aic_comparison` | Grouped bar chart of ΔAIC (relative to best) per copula family per index. Skipped if `all_results` table unavailable. |

#### Section 4 — SAF Curves

| Figure | File | Contents |
|---|---|---|
| Fig 09 | `Fig09_SAF_allindices_RP10` | Conditional SAF curves at 10-year RP, all four indices on one panel. |
| Fig 10 | `Fig10_SAF_allindices_RP100` | Conditional SAF curves at 100-year RP, all four indices on one panel. |
| Fig 11 | `Fig11_SAF_method_comparison_4panel` | 4-panel grid: Conditional vs Kendall SAF for RP = 10, 25, 50, 100 years. |
| Fig 12 | `Fig12_SAF_by_duration_class` | Per-class SAF curves (Conditional) for SPI-1 and SPI-3, stacked. Skipped if per-class SAF unavailable. |

#### Section 5 — Non-Stationarity

| Figure | File | Contents |
|---|---|---|
| Fig 13 | `Fig13_NS_epoch_comparison_4panel` | 4-panel: NS SAF reference vs recent epoch (Conditional method), all indices. Requires `enh_results[[nm]]$saf_ns_ref/rec`. |
| Fig 14 | `Fig14_GAMLSS_model_selection_table` | GAMLSS model-selection summary rendered as a `tableGrob`. Change-point year parsed from `ns$epoch_src` (FIX-1F). |
| Fig 15 | `Fig15_severity_trend_by_epoch` | Split violin + boxplot of severity by epoch (reference vs recent), 4-panel by index. |

#### Section 6 — Copula Diagnostics & Change-Point

| Figure | File | Contents |
|---|---|---|
| Fig 16 | `Fig16_changepoint_summary` | Change-point years from Pettitt / CUSUM / PRUTF tests (jittered points by significance), consensus year as dashed horizontal line. Uses corrected field names (FIX-1B). |
| Fig 17 | `Fig17_<INDEX>_rosenblatt` / `Fig17_rosenblatt_4panel_combined` | Per-index Rosenblatt PIT diagnostics (e2 histogram, QQ-plot, e1 vs e2 scatter). Also combined into a 4-panel composite. Copies are saved to individual index analysis folders. |
| Fig 18 | `Fig18_rolling_tau_4panel` | 4-panel rolling Kendall τ (10-yr window) with LOESS smooth; red dashed line marks consensus change-point year (FIX-1D). |

#### Section 7 — Excel Summary Tables

| Table | File | Contents |
|---|---|---|
| Table 01 | `Table01_return_periods_all_indices.xlsx` | Consolidated return-period table for all four indices (header formatted in blue). |
| Table 02 | `Table02_copula_fit_summary.xlsx` | Copula AIC/BIC/LogLik comparison; best-fit rows highlighted in green. |
| Table 03 | `Table03_nonstationary_summary.xlsx` | Non-stationary marginal summary: CP year (FIX-1F), epoch periods, GAMLSS AICs, epoch means, TV copula significance (FIX-1C). |

A final summary counts are printed: total PDFs, JPEGs, and XLSX files in `pub_figures/`.

---

## Output Directory Structure

```
spatial_drought/
├── SAF_analysis_log.txt                               # Timestamped log of all steps
├── SPI1_analysis/
│   ├── SPI1_drought_characteristics.csv
│   ├── SPI1_duration_classified_events.csv
│   ├── SPI1_return_periods_by_class.csv
│   ├── SPI1_severity_dist_comparison.csv
│   ├── SPI1_copula_comparison.csv
│   ├── SPI1_SAF_conditional_all.csv
│   ├── SPI1_SAF_kendall_all.csv
│   ├── SPI1_SAF_overlay.pdf
│   ├── SPI1_changepoint_detection.csv
│   ├── SPI1_segmented_copula.csv                      # written if change-point detected
│   ├── SPI1_nonstationary_summary.csv                 # written if ext2 NS run manually
│   ├── SPI1_SAF_nonstationary_reference.csv           # written if ext2 NS run manually
│   ├── SPI1_SAF_nonstationary_recent.csv              # written if ext2 NS run manually
│   ├── SPI1_SAF_nonstationary_kendall_reference.csv   # written if ext2 NS run manually
│   ├── SPI1_SAF_nonstationary_kendall_recent.csv      # written if ext2 NS run manually
│   ├── SPI1_rolling_tau.pdf                           # diagnostic plot (Step 5c)
│   ├── SPI1_copula_diagnostic.pdf                     # written if bootstrap GOF rejected
│   ├── SPI1_rosenblatt_diagnostic.pdf                 # written in Step 5 (event-level data)
│   ├── SPI1_rosenblatt_diagnostic.jpeg                # same, JPEG copy
│   ├── SPI1_tv_parameter.pdf                          # written if TV model selected
│   └── (per-class SAF CSVs when ≥15 events per class)
├── SPEI1_analysis/   (same structure)
├── SPI3_analysis/    (same structure, no D1-2 class, D0 events filtered)
├── SPEI3_analysis/   (same structure, no D1-2 class, D0 events filtered)
└── pub_figures/                                       # Step 7 outputs
    ├── Fig01_basin_drought_timeseries.pdf / .jpeg
    ├── Fig02_event_characteristics_scatter.pdf / .jpeg
    ├── ...
    ├── Fig18_rolling_tau_4panel.pdf / .jpeg
    ├── Table01_return_periods_all_indices.xlsx
    ├── Table02_copula_fit_summary.xlsx
    └── Table03_nonstationary_summary.xlsx
```

If optional Step 6 maps are enabled, a separate `drought_maps/` folder is created containing `<INDEX>_freq.nc` and `<INDEX>_sev.nc`.

> **Note:** The output root is `spatial_drought/` (not `drought_analysis/`). `quick_diagnostic()` was patched in `RUN_ALL.R` (FIX-5) to search `spatial_drought` instead of the hardcoded `drought_analysis` path in `QuickStart.R`.

---

## Is `Nechako_Drought_Methodology_Guide.R` Necessary?

**No — it is not required for the core analysis.** Here is the breakdown:

| Scenario | Need the Guide? |
|---|---|
| Running the full pipeline via `RUN_ALL.R` | It is sourced in Step 2, but only *loads* functions. Removing it would not break Steps 3–7. |
| Using Step 6 optional analyses | **Yes**, because `analyze_multiple_thresholds()`, `seasonal_drought_analysis()`, `map_drought_frequency()`, and `cross_correlate_indices()` are only defined in this file. |
| Running `SAF_Analysis.R` standalone | No |
| Running ext1 or ext2 | No |
| Generating Step 7 publication figures | No |

**Verdict:** The Methodology Guide is a collection of supplementary diagnostic and visualisation tools. The SAF curves, copula fits, return period tables, non-stationary analyses, publication figures, Excel tables, and all CSV/PDF outputs produced in Steps 3–7 are entirely independent of it.

---

## Key Design Decisions & Bug Fixes

### Active Fixes in This Version (RUN_ALL.R fix log)

| Fix ID | Location | Issue | Resolution |
|---|---|---|---|
| FIX-1A | `RUN_ALL.R` Step 5 | `derive_SAF_nonstationary*()` calls missing required argument `ns_result = ns` | Added `ns_result = ns` to all four calls (ref/rec × conditional/Kendall) |
| FIX-1B | `RUN_ALL.R` Fig 16 | Change-point field names incorrect | Corrected `pettitt_yr` → `pettitt_cp`, `cusum_yr` → `cusum_cp`, `prutf_yr` → `prutf_cp`, `consensus_yr` → `cp_year` |
| FIX-1C | `RUN_ALL.R` Table 3 | TV-copula significance field incorrect | Corrected `tv$lr_sig` → `tv$significant` |
| FIX-1D | `RUN_ALL.R` Fig 18 | Change-point field incorrect | Corrected `cp_result$consensus_yr` → `cp_result$cp_year` |
| FIX-1E | `RUN_ALL.R` Fig 6 | Dependency Kendall τ field incorrect | Corrected `dep$kendall_tau` → `dep$kendall` |
| FIX-1F | `RUN_ALL.R` Fig 14 / Table 3 | `ns$cp_year` slot does not exist in `fit_nonstationary_marginals()` return | `cp_year` now parsed from `ns$epoch_src` using a 4-digit regex match, consistent with `.resolve_epoch()` |
| FIX-5 | `RUN_ALL.R` Step 1 | `quick_diagnostic()` hardcoded `"drought_analysis"` path; all output goes to `"spatial_drought"` | Patched inline after sourcing `QuickStart.R`; also upgrades event count to use raw CSV as authoritative source |
| FIX-8 | `RUN_ALL.R` Step 4 | Version label printed `"ext2 v4"` | Corrected to `"ext2 v5"` (current version is v7) |

### Ongoing Design Notes

1. **tCopula integer df (`.fix_tcopula_df()`):** `fitCopula()` optimises the Student-t degrees-of-freedom as a continuous parameter. `pCopula()` requires an integer `df`. The helper `fix_tcopula_df()` rounds to the nearest integer (minimum 2) and rebuilds the copula object. Applied at every point where `pCopula()` or `cCopula()` is called downstream.

2. **`rm(list=ls())` removed:** Earlier versions wiped the workspace mid-pipeline. Cleanup is now the caller's responsibility.

3. **Copula log-likelihood saved:** `fit_copulas()` stores `LogLik` alongside AIC/BIC in the comparison CSV and result list.

4. **Minimum events threshold = 15:** Per-class SAF fitting is skipped for classes with fewer than 15 events to avoid unreliable tail fits (raised from 10).

5. **D0 event filtering for scale-3 indices:** `identify_duration_classified_events()` removes all `"D0 (Sub-threshold)"` events (duration < 3 months) for SPI-3/SPEI-3 before returning. Sub-threshold events are physically meaningless for a 3-month accumulation scale and must not contaminate downstream statistics.

6. **CFG–Pickands upper-tail estimator:** `compute_cfg_pickands()` implements the correct CFG estimator `A(t) = exp(−E[log(min(−log U/(1−t), −log V/t))])` for the Pickands dependence function, from which `λ_U = max(0, 2(1 − A(0.5)))`.

7. **Gamma renewal process for SAF target:** `.compute_renewal_saf_target()` replaces the homogeneous Poisson assumption when event inter-arrival CV falls outside [0.5, 1.5]. Shape parameter `k = 1/CV²` controls clustering (k < 1) or regularity (k > 1).

8. **GAMLSS coefficient extraction:** Some `gamlss` versions return a closure from `predict(..., type = "response")`. The fix extracts coefficients directly and applies `exp()` to the log-link manually: `mu(t) = exp(b0 + b1 × year_std(t))`.

9. **Epoch boundary correction in ext2 (v4 → v7):** `yr_match <- dc$year` uses the actual calendar year stored in each drought-month row. `yr_range` now correctly spans 1950–2025 and `recent_years` extends to the record end.

   *Impact on epoch boundaries after fix:*

   | Index | CP year | Reference period | Recent period |
   |-------|---------|-----------------|---------------|
   | SPI-1  | 1965 | 1950–1964 | 1965–2025 (+27 yr) |
   | SPEI-1 | 1996 | 1950–1995 | 1996–2025 (CP used correctly) |
   | SPI-3  | 1973 | 1950–1972 | 1973–2025 (+27 yr) |
   | SPEI-3 | 1987 | 1950–1986 | 1987–2025 (+26 yr) |

10. **Non-stationary framework consistency (v5 → v7):** Earlier versions (v4) paired a step-function marginal (discrete epoch from `fit_nonstationary_marginals`) with a smooth-trend copula collapsed to a single midpoint value. This mixed two incompatible temporal representations. The revised `.build_period_copula()` applies a three-path priority: segmented copula (Path 1, epoch-consistent), TV midpoint (Path 2, warned hybrid), stationary fallback (Path 3).

11. **Rosenblatt PIT h-function (`.h_func_fd()`):** `e2` is computed as `∂C(u_S, u_A)/∂u_S` using analytical formulae for Gumbel, Frank, Joe, Student-t, and Gaussian families. A `cCopula`-based path with finite-difference fallback handles all other families. This eliminates the NA-return failure from earlier versions and produces correctly independent Uniform(0,1) variates when the copula is well-specified.

12. **Event-level Rosenblatt diagnostic:** Step 5 applies the Rosenblatt PIT using **event-level** pseudo-observations (ranked `max_severity`, `max_area_pct` per event) rather than the monthly time series. This removes intra-event temporal clustering and restores the i.i.d. assumption required for the test.

---

## Known Diagnostic Warnings and Their Interpretation

| Warning / Log message | Source | Severity | Cause | Action |
|---|---|---|---|---|
| `WARNING [TD-U]: AIC-best copula has no upper-tail dependence but empirical λ_U > 0.10` | `fit_copulas()` | Advisory | AIC prefers Frank, Plackett, or Gaussian despite upper-tail signal; override logic may have failed | Compare SAF curves from Gumbel/Joe vs Frank at high return periods |
| `TD-MISMATCH OVERRIDE: Selected '<family>' has λ_U=0 but empirical λ_U > 0.10. Switching to tail-dependent candidate.` | `fit_copulas()` | Informational | Zero-tail family selected despite signal; override applied automatically | Review final copula selection in `<INDEX>_copula_comparison.csv` |
| `TAIL-DEPENDENCE PRIORITY: λ_U=X > |τ|=Y. Enforcing Joe/Gumbel selection.` | `fit_copulas()` | Informational | CFG estimator confirms upper-tail signal stronger than overall correlation | Normal behaviour when extreme events co-occur more often than average |
| `[INDEX] Rolling tau: N window(s) had zero variance … set to NA` | `detect_copula_changepoints()` | Informational | Drought-month windows with identical severity or area values (common in SPI-3) | Inspect `<INDEX>_drought_characteristics.csv`; change-point results unaffected |
| `MISMATCH — investigate aggregation!` in `quick_diagnostic()` | `QuickStart.R` (patched) | Warning | `nrow(ev_file) ≠ sum(rp$n_events)`; D0 exclusion may cause discrepancy | Check `D0 (Sub-threshold)` handling and confirm filter was applied before CSV write |
| `TV trend skipped for StudentT: multi-parameter families require extended formulation` | `fit_timevarying_copula()` | Informational | StudentT has 2 parameters; scalar log-linear TV model not applicable | Stationary StudentT fit used; segmented model still evaluated |
| `WARNING [EPOCH MISMATCH]: copula cp_year=X != marginal cp_year=Y` | `.assert_epoch_consistency()` in ext2 | Warning | The change-point years passed to `fit_timevarying_copula()` and `fit_nonstationary_marginals()` differ | Ensure both functions receive the same `cp_result` object |
| `PATH 2 … WARNING: marginal is a step-function; copula is a point-evaluated smooth trend.` | `.build_period_copula()` in ext2 | Advisory | `seg_result` unavailable; TV copula evaluated at epoch midpoint instead | Check if segment sizes are ≥15; this is an acknowledged hybrid approximation |
| `CRITICAL: Period-specific K_C computation failed … Aborting Kendall SAF` | `derive_SAF_nonstationary_kendall()` in ext2 | Error | All three tiers of K_C computation failed for non-Archimedean family | Check copula object integrity for this epoch; Conditional NS SAF is unaffected |

---

## Glossary

| Term | Meaning |
|---|---|
| SPI | Standardised Precipitation Index — drought index based on precipitation anomaly |
| SPEI | Standardised Precipitation-Evapotranspiration Index — includes temperature/ET |
| SAF curve | Severity–Area–Frequency curve: relates drought severity to spatial extent for a given return period |
| Copula | A function that models the dependence structure between two or more variables independently of their marginal distributions |
| Kendall distribution K_C | The CDF of copula values: `K_C(t) = P(C(U,V) ≤ t)` |
| IFM | Inference Functions for Margins — a two-stage MLE method for copula fitting |
| GAMLSS | Generalised Additive Models for Location, Scale and Shape — used here for non-stationary distributional regression |
| Return period | Expected average time between events of at least a given magnitude |
| Duration class | Grouping of drought events by length in months (D1-2, D3-6, D7-12, D13+) |
| TV copula | Time-varying copula: a copula whose dependence parameter θ is modelled as a log-linear function of standardised year |
| Segmented copula | Two independent stationary copulas fitted to temporal sub-records defined by a consensus change-point year |
| Rosenblatt PIT | Probability Integral Transform (Rosenblatt 1952): converts a bivariate sample from a correctly specified copula into two independent U[0,1] variates for GOF testing |
| λ_U | Upper-tail dependence coefficient: limiting probability of joint exceedances as the threshold approaches 1 |
| CFG estimator | Capéraà–Fougères–Genest non-parametric estimator for the Pickands dependence function A(t); used to compute empirical λ_U |
| Gamma renewal process | A renewal process where inter-arrival times follow a Gamma distribution; generalises Poisson (exponential IAT, k=1) to clustered (k<1) and regular (k>1) processes |
| `.build_period_copula()` | Internal ext2 helper that selects the period-specific copula via a three-path priority rule (segmented > TV midpoint > stationary) |
| `.resolve_epoch()` | Internal ext2 helper that derives epoch year vectors exclusively from `ns_result` as the single source of truth for epoch boundaries |
| `.assert_epoch_consistency()` | Internal ext2 helper that compares marginal and copula change-point years and logs a warning on mismatch |
| copula_path | Output column in NS SAF CSVs recording which copula-selection path was taken (`"segmented"`, `"tv_midpoint"`, or `"stationary"`) |
| pub_figures/ | Step 7 output folder containing publication-quality PDFs, JPEGs, and Excel tables |
