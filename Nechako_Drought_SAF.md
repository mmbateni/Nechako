# Nechako Drought SAF Analysis — README

## Overview

This codebase performs a **Severity–Area–Frequency (SAF) drought analysis** for the Nechako Basin in British Columbia, Canada. It uses gridded **SPI** (Standardised Precipitation Index) and **SPEI** (Standardised Precipitation-Evapotranspiration Index) raster time series as inputs, and produces statistically rigorous SAF curves, copula models, return period tables, and diagnostic plots as outputs.

The analysis follows a multivariate frequency analysis framework: drought severity and drought-affected area fraction are modelled as two correlated random variables, joined by a copula. SAF curves answer the question: *"For a given return period, how severe is a drought that affects X% of the basin?"*

---

## File Inventory & Roles

| File | Role | Required? |
|---|---|---|
| `Nechako_Drought_RUN_ALL.R` | **Master orchestrator** — single entry point that runs all steps in order | Yes |
| `Nechako_Drought_QuickStart.R` | Input verification and pre-flight checks | Yes |
| `Nechako_Drought_SAF_Analysis.R` | Core analysis pipeline (data loading → SAF curves) | Yes |
| `Nechako_Drought_SAF_Analysis_ext1.R` | Extension: dependency tests, copula GOF, time-varying and segmented copula | Yes |
| `Nechako_Drought_SAF_Analysis_ext2.R` | Extension: non-stationary marginals, Kendall-corrected NS SAF, event placement | Yes |
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
| `gridExtra` | Multi-panel PDF plots |
| `viridis` | Colour scales |
| `dplyr` | Data manipulation |
| `scales` | Axis formatting |
| `moments` | Higher-order statistics (ext1) |
| `openxlsx` | Excel export (ext1) |
| `Kendall` | Rank correlation tests (ext1) |
| `gamlss` | Non-stationary GAMLSS models (ext2) |
| `gamlss.dist` | GAMLSS distribution families (ext2) |
| `trend` | Pettitt test (ext1) |

---

## How to Run

The simplest approach is the single entry point:

```r
setwd("D:/Nechako_Drought/Nechako")
source("Nechako_Drought_RUN_ALL.R")
```

This runs all six steps (see below) in the correct order. Do **not** source the individual analysis scripts out of order — `ext1.R` and `ext2.R` depend on objects created by `SAF_Analysis.R`, and `SAF_Analysis.R` depends on objects created by `QuickStart.R`.

---

## Pipeline Steps (as executed by `RUN_ALL.R`)

### Step 1 — Input Verification (`QuickStart.R`)

- Spot-checks **4 representative NetCDF files** (one per index, January only) as proxies for all 48 required files (4 indices × 12 months). A missing file stops execution; the remaining 44 files are verified implicitly when the pipeline loads them in Step 3.
- Loads `SPI-1` for January and validates that pixel values fall within the physically plausible range of −5 to +5.
- Runs a mathematical sanity check on the Kendall distribution formula using the known analytical result for the independence copula at `t = 0.5`.
- Prints the recommended run order to the console.

### Step 2 — Methodology Helpers (`Methodology_Guide.R`)

Loads five utility functions into the workspace. **No computation is triggered at load time.** These functions are available for optional Step 6 calls:

- `analyze_multiple_thresholds()` — repeats the area fraction computation at drought thresholds −0.5, −1.0, and −1.5 and produces a PDF line chart.
- `seasonal_drought_analysis()` — computes mean severity and mean area fraction by calendar month and saves a two-panel bar chart PDF.
- `map_drought_frequency()` — computes per-pixel drought frequency (% of months below threshold) and mean per-pixel severity, saving both as NetCDF files.
- `compare_return_period_methods()` — plots Conditional vs Kendall SAF curves on a single overlay for visual method comparison.
- `cross_correlate_indices()` — scatter-plots SPI-1 vs SPI-3 event severity and reports the Pearson correlation.

### Step 3 — Main SAF Pipeline (`SAF_Analysis.R`)

This is the computational core. It runs the full pipeline four times — once each for **SPI-1**, **SPEI-1**, **SPI-3**, and **SPEI-3** — and writes results to `drought_analysis/<INDEX>_analysis/`.

#### 3a. Data Loading

All 12 monthly NetCDF tiles for each index are loaded with `terra::rast()` and concatenated into a single time-ordered raster stack. Cell areas (km²) and basin total area are computed once from the SPI-1 grid.

#### 3b. Drought Characteristics Extraction

For every time step `t`, the script:
1. Reads all pixel values from the index raster layer.
2. Counts pixels below the drought threshold (default −0.5) as the **drought area** (expressed as % of valid pixels).
3. Computes **mean drought severity** as the mean exceedance below threshold: `mean(threshold − value)` for all dry pixels.

Output: a time series data frame `[date, year, month, severity, area_pct]` saved as `<INDEX>_drought_characteristics.csv`.

#### 3c. Event Identification & Duration Classification

Consecutive months with `severity > 0` and `area_pct > 0` are grouped into discrete drought **events** using a simple run-length algorithm. Each event is described by:

- Start/end date and year
- Duration in months
- Duration class (based on scale):

| Class | SPI-1 / SPEI-1 | SPI-3 / SPEI-3 |
|---|---|---|
| D1-2 (Very-short) | 1–2 months | — |
| D3-6 (Short-term) | 3–6 months | 3–6 months |
| D7-12 (Medium-term) | 7–12 months | 7–12 months |
| D13+ (Long-term) | ≥13 months | ≥13 months |

- Mean and max severity
- Mean and max area fraction

Empirical **return periods** are computed as `record_length_years / n_events` per class.

Outputs: `<INDEX>_duration_classified_events.csv`, `<INDEX>_return_periods_by_class.csv`.

#### 3d. Marginal Distribution Fitting

Using only months with active drought (`severity > 0`, `area_pct > 0`):

- **Area fraction** is fitted with a **Beta distribution** (values bounded 0–1).
- **Severity** is fitted with six distributions (Exponential, Gamma, Weibull, Log-Normal, Normal, Logistic) using MLE. The best-fitting distribution is selected by minimum AIC.

Output: `<INDEX>_severity_dist_comparison.csv`.

#### 3e. Copula Fitting

The parametric marginal CDFs are used to transform severity and area into uniform pseudo-observations `(u_S, u_A)`. **Seven** copula families are then fitted to these pseudo-observations using the IFM (Inference Functions for Margins) maximum pseudo-likelihood method:

| Copula | Tail dependence |
|---|---|
| Survival Clayton | Upper tail |
| Gumbel | Upper tail |
| Frank | Symmetric (no tail) |
| Joe | Upper tail |
| Plackett | Symmetric (no tail) |
| Gaussian (Normal) | No tail (asymptotic independence for all \|ρ\| < 1) |
| Student-t | Both tails (symmetric) |

The best copula is selected by minimum AIC. AIC, BIC, and log-likelihood are all saved. An empirical upper-tail dependence coefficient λ_U is computed non-parametrically (Schmidt & Stadtmüller 2006) and compared to the theoretical λ_U of the winning family. A `WARNING [TD-U]` is logged if the AIC-best family has zero upper-tail dependence (Frank, Plackett, or Gaussian) but empirical λ_U > 0.10.

> **Note:** The same seven-family candidate set is used in `fit_copulas()` (Step 3), `detect_copula_changepoints()` / `fit_segmented_copula()` (Step 5c), and the PRUTF structural-break test in ext1. The Student-t family has two parameters (ρ, df); because the log-linear TV formulation is scalar, the TV copula trend model is skipped for Student-t and the stationary fit is used instead (a warning is logged).

Output: `<INDEX>_copula_comparison.csv`.

#### 3f. SAF Curve Derivation — Two Methods

**Method A: Conditional Copula**

For a target return period `T` years, the effective monthly non-exceedance probability is:

```
target = 1 - mu_T / (T × 12)
```

where `mu_T` is the mean inter-drought interval in months. For each area value `a`, the script solves numerically for the severity `s` such that:

```
P(S > s | A = a) = 1 − target
```

using the conditional copula `∂C(u_S, u_A)/∂u_A`.

**Method B: Kendall Distribution**

Instead of conditioning on area, this method uses the joint level set of the copula. It solves for the copula level `t_val` such that:

```
K_C(t_val) = target
```

where `K_C(t) = P(C(U, V) ≤ t)` is the **empirical Kendall distribution**, estimated directly from the fitted pseudo-observations. Then for each area value `a`, it finds severity `s` such that `C(u_S, u_A) = t_val`.

Both methods produce curves for return periods of **10, 25, 50, and 100 years** across area fractions from 5% to 95%.

SAF curves are also computed **per duration class** when the class has ≥15 events.

Outputs per index: `<INDEX>_SAF_conditional_all.csv`, `<INDEX>_SAF_kendall_all.csv`, `<INDEX>_SAF_overlay.pdf`, and per-class equivalents.

The four result objects `res_spi1`, `res_spei1`, `res_spi3`, `res_spei3` remain in the R workspace for use by Steps 4–5.

### Step 4 — Load Extensions (`ext1.R` + `ext2.R`)

The extension scripts are sourced here. They only define functions — no heavy computation occurs at load time (each script emits a single `log_event()` message confirming successful loading). Both scripts verify that `log_event()` from `SAF_Analysis.R` is present before loading; an error stops execution if the main pipeline was skipped.

**Functions exported by `ext1.R` (v6):**
`test_dependency`, `gof_copula`, `detect_copula_changepoints`, `rosenblatt_pit_gof`, `fit_timevarying_copula`, `fit_segmented_copula`

**Functions exported by `ext2.R` (v4):**
`fit_nonstationary_marginals`, `derive_SAF_nonstationary`, `derive_SAF_nonstationary_kendall`, `place_event_on_saf`

### Step 5 — Run Enhancements

For each of the four indices, five enhancement analyses are run:

#### 5a. Dependency Test

Computes **Kendall's τ** and **Spearman's ρ** between drought severity and area fraction, with a significance test. This validates whether the copula is statistically justified (if τ ≈ 0, the variables are nearly independent and a copula provides little benefit).

#### 5b. Copula Goodness-of-Fit (Bootstrap Sn)

Performs a parametric bootstrap GOF test using the Cramér–von Mises `Sn` statistic with 499 bootstrap replicates. A p-value > 0.05 means the chosen copula family is not rejected at the 5% level.

Output: logged to console and saved indirectly via `log_event`.

#### 5c. Time-Varying & Segmented Copula

Tests whether the dependence structure between severity and area has changed over time. The analysis proceeds in three stages:

**Change-point detection** (Pettitt, CUSUM, PRUTF) is run first on a rolling Kendall-τ series to identify abrupt structural breaks. A break is flagged when at least 2 of the 3 tests reject H₀ at the 5% level; the consensus change-point year is the median of the significant test estimates.

**Three copula models** are then fitted and compared by AIC and BIC:

| Model | Description | Free parameters |
|---|---|---|
| Stationary | Single fixed copula parameter | 1 |
| Time-varying (TV) | Log-linear trend: `θ(t) = ilink(a + b × year_std)` | 2 |
| Segmented | Two independent stationary copulas, one per epoch defined by the change-point year | 2 |

Model selection logic: TV is preferred when the LRT is significant (p < 0.05), no strong abrupt break is detected, TV has lower AIC than the segmented model, and the Rosenblatt PIT GOF passes. Segmented is preferred when a change point is detected, the segmented model has lower AIC, and the Rosenblatt diagnostic is better. Otherwise the stationary model is retained.

**Post-fitting**, a **Rosenblatt Probability Integral Transform (PIT)** Goodness-of-Fit test validates both the TV and segmented copulas by checking uniformity (KS test on `e2`) and independence (Kendall τ test between `e1` and `e2`) of the transformed pseudo-observations.

Diagnostic PDFs produced conditionally:

| File | Condition |
|---|---|
| `<INDEX>_rolling_tau.pdf` | Always written |
| `<INDEX>_copula_diagnostic.pdf` | GOF bootstrap rejected (p < 0.05) |
| `<INDEX>_rosenblatt_diagnostic.pdf` | GOF rejected AND TV copula was fitted with a valid `a_hat` |
| `<INDEX>_tv_parameter.pdf` | TV model selected |

Outputs: `<INDEX>_timevarying_copula.csv`, `<INDEX>_segmented_copula.csv` (written only when a change point is detected and both segments have ≥ 15 observations).

#### 5d. Non-Stationary Marginals (GAMLSS)

Fits **three GAMLSS Gamma models** to drought severity and selects among them by AIC:

| Model | mu formula | sigma formula |
|---|---|---|
| M0 (stationary) | `~ 1` | `~ 1` |
| M1 (mu-trend) | `~ year_std` | `~ 1` |
| M2 (mu+sigma-trend) | `~ year_std` | `~ year_std` |

Period-specific mean severity `mu` and coefficient of variation `sigma` are computed for two epochs:
- **Reference period**: first 40% of the record — or from record start to `cp_year − 1` if a change point is detected and falls between 20%–80% of the record.
- **Recent period**: last 30% of the record — or from `cp_year` to record end when a change point is used.

(Override by passing explicit `ref_years` and `recent_years` vectors to `fit_nonstationary_marginals()`.)

Period-specific `mu` and `sigma` are extracted directly from model coefficients (not via `predict()`) by applying the GA log-link manually: `mu(t) = exp(b0 + b1 × year_std(t))`. Values are averaged over the epoch years. Both `mu` and `sigma` are passed to Step 5e.

Output: `<INDEX>_nonstationary_summary.csv`.

#### 5e. Non-Stationary SAF Curves

Uses the period-specific Gamma parameters `mu` and `sigma` (from 5d) while keeping the Beta area marginal unchanged. Produces separate SAF curves for the reference and recent periods using both:

- **Method A (Conditional_NS)** — conditional copula approach with period-specific severity marginal.
- **Method B (Kendall_NS)** — Kendall-corrected approach with period-specific severity marginal (consistent with De Michele et al. 2026, Eq. 11).

**Copula selection in NS SAF:** When the time-varying copula analysis (5c) found a statistically significant trend in dependence (`lr_p < 0.05`), the copula used in each SAF curve is replaced with a **period-specific instance** evaluated at the epoch midpoint year using the fitted TV copula (via the internal helper `.build_period_copula()`). This propagates non-stationarity in both the severity marginal and the dependence structure into the SAF curve. When the TV trend is not significant, the stationary copula from Step 3 is used unchanged for both periods.

This allows detection of **shifts in drought severity for the same return period** between historical and recent climates.

Outputs: `<INDEX>_SAF_nonstationary_<period>.csv`, `<INDEX>_SAF_nonstationary_kendall_<period>.csv`.

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

## Output Directory Structure

```
drought_analysis/
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
│   ├── SPI1_timevarying_copula.csv
│   ├── SPI1_segmented_copula.csv                      # written if change-point detected
│   ├── SPI1_nonstationary_summary.csv
│   ├── SPI1_SAF_nonstationary_reference_<yr>_<yr>.csv
│   ├── SPI1_SAF_nonstationary_recent_<yr>_<yr>.csv
│   ├── SPI1_SAF_nonstationary_kendall_reference_<yr>_<yr>.csv
│   ├── SPI1_SAF_nonstationary_kendall_recent_<yr>_<yr>.csv
│   ├── SPI1_rolling_tau.pdf                           # always written in Step 5c
│   ├── SPI1_copula_diagnostic.pdf                     # written if bootstrap GOF rejected
│   ├── SPI1_rosenblatt_diagnostic.pdf                 # written if GOF rejected & TV fitted
│   ├── SPI1_tv_parameter.pdf                          # written if TV model selected
│   └── (per-class SAF CSVs when ≥15 events per class)
├── SPEI1_analysis/   (same structure)
├── SPI3_analysis/    (same structure, no D1-2 class)
└── SPEI3_analysis/   (same structure, no D1-2 class)
```

If optional Step 6 maps are enabled, a separate `drought_maps/` folder is created containing `<INDEX>_freq.nc` and `<INDEX>_sev.nc`.

---

## Is `Nechako_Drought_Methodology_Guide.R` Necessary?

**No — it is not required for the core analysis.** Here is the breakdown:

| Scenario | Need the Guide? |
|---|---|
| Running the full pipeline via `RUN_ALL.R` | It is sourced in Step 2, but only *loads* functions. Removing it would not break Steps 3–5. |
| Using Step 6 optional analyses | **Yes**, because `analyze_multiple_thresholds()`, `seasonal_drought_analysis()`, `map_drought_frequency()`, and `cross_correlate_indices()` are only defined in this file. |
| Running `SAF_Analysis.R` standalone | No |
| Running ext1 or ext2 | No |

**Verdict:** The Methodology Guide is a collection of supplementary diagnostic and visualisation tools. The SAF curves, copula fits, return period tables, non-stationary analyses, and all CSV/PDF outputs produced in Steps 3–5 are entirely independent of it.

---

## Key Design Decisions & Bug Fixes

The comments in the scripts document several important corrections relative to earlier versions:

1. **Kendall distribution formula** (`QuickStart.R`, `SAF_Analysis.R`): The pipeline uses `kc_fn(t) = mean(C(U,V) <= t)` (empirical K_C), verified at startup against the known analytical result for the independence copula `K_C(t) = t - t*log(t)`.
2. **`rm(list=ls())` removed** (`SAF_Analysis.R`, `ext1.R`, `ext2.R`): Earlier versions wiped the workspace mid-pipeline. Cleanup is now the caller's responsibility.
3. **Copula log-likelihood in comparison table** (`SAF_Analysis.R`): Saved alongside AIC/BIC.
4. **Minimum events threshold raised to 15** (`SAF_Analysis.R`): Per-class SAF fitting is skipped for classes with fewer than 15 events to avoid unreliable tail fits.
5. **Time-varying copula family** (`ext1.R`): Uses whichever family was selected by AIC in Step 3. Multi-parameter families (Student-t) fall back to the stationary model with a logged warning.
6. **GAMLSS `predict()` workaround** (`ext2.R`): Some `gamlss` package versions return a closure instead of a numeric from `predict(..., type="response")`. The fix extracts coefficients directly and applies `exp()` manually to implement the log-link.
7. **Non-stationary year window** (`ext2.R`): Reference and recent periods are derived from the actual record length (40%/30% split) instead of hardcoded calendar years. The detected change-point year takes precedence when it falls between 20%–80% of the record.
8. **Rosenblatt PIT** (`ext1.R` — `rosenblatt_pit_gof()`):

   `e2` is computed as the central finite-difference derivative of `pCopula` with respect to `u1` (helper `.h_func_fd`). `pCopula` is implemented for every copula family in the package (Frank, Gumbel, Plackett, Clayton, Joe, SurvClayton, Gaussian, StudentT), eliminating the NA-return failure and producing correctly independent Uniform(0,1) variates when the copula is well-specified.

9. **Epoch boundary truncation** (`ext2.R` — `fit_nonstationary_marginals()`):

   `yr_match <- dc$year` uses the actual calendar year stored in each drought-month row. `yr_range` now correctly spans 1950–2025 and `recent_years` extends to the record end.

   *Impact on epoch labels after fix:*

   | Index | CP year | Reference period | Recent period |
   |-------|---------|-----------------|----------------|
   | SPI-1  | 1965 | 1950–1964 | 1965–2025 (+27 yr) |
   | SPEI-1 | 1996 | 1950–1995 | 1996–2025 (CP used correctly) |
   | SPI-3  | 1973 | 1950–1972 | 1973–2025 (+27 yr) |
   | SPEI-3 | 1987 | 1950–1986 | 1987–2025 (+26 yr) |

10. **`total_events` authoritative source** (`QuickStart.R` — `quick_diagnostic()`):

    `total_events` is now read directly from `nrow(<INDEX>_duration_classified_events.csv)` when the file exists. `sum(rp$n_events)` is retained as a fallback and cross-check target. A `MISMATCH` flag is printed when the two counts differ.

11. **SPI-3 zero-variance rolling-window warnings** (`ext1.R` — `detect_copula_changepoints()`):

    Both `sd(severity_window)` and `sd(area_pct_window)` are checked before calling `cor()`. Zero-variance windows are assigned NA silently; the total count is reported once via `log_event()`. Change-point results are not affected.

    **Bootstrap and permutation window sizes:**
    - `gof_copula()`: **N_boot = 499 replicates** (parametric Cramér–von Mises Sn bootstrap). Minimum achievable two-sided p-value: 1/500 = 0.002. Use `N_boot = 999` for publication; `N_boot = 99` for fast development.
    - `detect_copula_changepoints()` CUSUM test: **499 Monte-Carlo permutations** (matching N_boot for consistency).

---

## Known Diagnostic Warnings and Their Interpretation

| Warning / Log message | Source | Severity | Cause | Action |
|---|---|---|---|---|
| `WARNING [TD-U]: AIC-best copula has no upper-tail dependence but empirical λ_U > 0.10` | `fit_copulas()` | Advisory | AIC prefers Frank, Plackett, or Gaussian despite upper-tail signal | Consider Gumbel, Joe, SurvClayton, or StudentT; compare SAF curves at high return periods |
| `[INDEX] Rolling tau: N window(s) had zero variance … set to NA` | `detect_copula_changepoints()` | Informational | Isolated drought-month windows with identical severity or area values (typical in SPI-3 record) | Inspect `<INDEX>_drought_characteristics.csv` for constant runs; change-point results unaffected |
| `MISMATCH — investigate aggregation!` in `quick_diagnostic()` | `QuickStart.R` | Warning | `nrow(ev_file) ≠ sum(rp$n_events)`; sub-threshold events may not map to a named class | Check duration classification logic and `D0 (Sub-threshold)` rows in the rp table |
| `TV trend skipped for StudentT: multi-parameter families require extended formulation` | `fit_timevarying_copula()` | Informational | StudentT has 2 parameters (ρ, df); the scalar log-linear TV model is not applicable | Stationary StudentT fit is used; segmented model is still evaluated if a change point is detected |

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
| λ_U | Upper-tail dependence coefficient: limiting probability of joint exceedances as the threshold approaches 1; measures the tendency for simultaneous extremes |
