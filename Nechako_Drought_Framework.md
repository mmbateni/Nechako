# Nechako River Basin — Drought Analysis Framework

**Study Area:** Nechako River Basin, British Columbia, Canada  
**Data Source:** ERA5-Land Monthly Reanalysis (Copernicus/ECMWF)  
**Analysis Period:** 1950–2025  
**Coordinate Reference System:** BC Albers Equal Area (EPSG:3005)  
**Status:** Includes counterfactual SPEI branch (month-specific temperature detrending) proposed by S. Dery (JHM manuscript revision, June 2026)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Repository Structure](#2-repository-structure)
3. [Execution Order](#3-execution-order)
4. [Shared Utilities — `DROUGHT_ANALYSIS_utils.R`](#4-shared-utilities)
5. [Script 1 — SPI Calculation (`1SPI_ERALand.R`)](#5-script-1--spi-calculation)
6. [Script 2a — Month-Specific Temperature Detrending (`2a_detrend_temperature.R`)](#6-script-2a--month-specific-temperature-detrending)
7. [Script 2b — PET Calculation (`2b_PET_ERALand.R`)](#7-script-2b--pet-calculation)
8. [Script 3 — SPEI Calculation (`3SPEI_ERALand.R`)](#8-script-3--spei-calculation)
9. [Script 4a — Pr/PET/Temperature Trend Analysis (`4pr_pet_trends.r`)](#9-script-4a--prpettemperature-trend-analysis)
10. [Script 4b — Trend Visualization (`4pr_pet_trends_visualization.r`)](#10-script-4b--trend-visualization)
11. [Script 5 — Animated GIFs (`5gif_animations.R`)](#11-script-5--animated-gifs)
12. [Script 6 — Spatial Trend Test for All Indices (`6trend_test_ALL.R`)](#12-script-6--spatial-trend-test-for-all-indices)
13. [Script 7 — Basin-Averaged Time Series (`7basin_timeseries.R`)](#13-script-7--basin-averaged-time-series)
14. [Script 8 — Spatial Trends Visualization (`8trends_visualization.R`)](#14-script-8--spatial-trends-visualization)
15. [Script 9 — Drought Event Ranking (`9event_ranking.R`)](#15-script-9--drought-event-ranking)
16. [Script 10a — Dynamic/Thermodynamic Decomposition (`10a_dynamic_thermodynamic_decomp.R`)](#16-script-10a--dynamicthermodynamic-decomposition)
17. [Script 10b — PET Bias Non-Stationarity (`10b_PET_bias_nonstationarity.R`)](#17-script-10b--pet-bias-non-stationarity)
18. [Script 11 — Counterfactual SPEI Comparison (`11_SPEI_counterfactual_comparison.R`)](#18-script-11--counterfactual-spei-comparison)
19. [Manuscript Support Script (`11nechako_drought_manuscript1.R`)](#19-manuscript-support-script)
20. [Statistical Methods](#20-statistical-methods)
21. [Key Parameters and Thresholds](#21-key-parameters-and-thresholds)
22. [Data Flow Diagram](#22-data-flow-diagram)
23. [R Package Dependencies](#23-r-package-dependencies)

---

## 1. Overview

This framework is a complete end-to-end pipeline for the **detection, quantification, attribution, and visualization of drought conditions** in the Nechako River Basin (BC, Canada) over the period 1950–2025, with a focus on the 2022–2025 drought event.

The pipeline ingests ERA5-Land monthly reanalysis fields (precipitation, temperature, wind, radiation, humidity, pressure) and produces:

- **Gridded drought indices** at multiple timescales: SPI, SPEI (Penman-Monteith and Thornthwaite), and SWEI
- **Pixel-level trend statistics**: Mann-Kendall (variance-corrected and TFPW), Sen's slope, regime shift year (PELT), spectral periodicity
- **Basin-averaged time series** for all indices at all timescales
- **Drought event catalogs**: duration, severity, intensity, return periods
- **Dynamic vs. thermodynamic decomposition** of drought drivers (precipitation+radiation+wind vs. temperature)
- **PET bias non-stationarity** corrections for the Thornthwaite–Penman-Monteith difference under climate change
- **Counterfactual SPEI** (proposed by S. Dery, JHM revision): SPEI recomputed from month-specific linearly detrended T2m and Tdew fields — both embedded in `2b_PET_ERALand.R` (Section 5c) and `3SPEI_ERALand.R` (Runs 3 & 4) via the `RUN_DETRENDED_BRANCH` flag — providing a direct comparison to the observed SPEI and isolating the drought intensification attributable to the long-term warming trend

---

## 2. Repository Structure

```
Nechako/
├── DROUGHT_ANALYSIS_utils.R            ← Shared utilities (sourced by all scripts)
├── 1SPI_ERALand.R                      ← Step 1: SPI indices
├── 2a_detrend_temperature.R            ← Step 2a: Month-specific T2m & Tdew detrending
├── 2b_PET_ERALand.R                    ← Step 2b: PET (observed + counterfactual via §5c)
├── 3SPEI_ERALand.R                     ← Step 3: SPEI (Runs 1–4: observed + counterfactual)
├── 4pr_pet_trends.r                    ← Step 4a: Pr/PET/T trend analysis
├── 4pr_pet_trends_visualization.r      ← Step 4b: Trend maps & plots
├── 5gif_animations.R                   ← Step 5: Animated drought GIFs
├── 6trend_test_ALL.R                   ← Step 6: Spatial trend test (SPI, SPEI_PM, SPEI_Thw)
├── 7basin_timeseries.R                 ← Step 7: Basin-averaged time series
├── 8trends_visualization.R             ← Step 8: Spatial diagnostic figures
├── 9event_ranking.R                    ← Step 9: Event ranking & return periods
├── 10a_dynamic_thermodynamic_decomp.R  ← Step 10a: Dynamic/Thermodynamic decomposition
├── 10b_PET_bias_nonstationarity.R      ← Step 10b: PET bias non-stationarity
├── 11_SPEI_counterfactual_comparison.R ← Step 11: Observed vs. counterfactual comparison
├── 11nechako_drought_manuscript1.R     ← Manuscript support (tables, blanks, figures)
│
├── Spatial/
│   └── nechakoBound_dissolve.kmz       ← Basin boundary (KMZ/KML)
├── monthly_data_direct/                ← ERA5-Land NetCDF inputs + PET outputs
│   ├── 2m_temperature_monthly.nc
│   ├── 2m_dewpoint_temperature_monthly.nc
│   ├── [all other ERA5-Land inputs]
│   ├── potential_evapotranspiration_monthly.nc             ← PM PET (observed)
│   ├── potential_evapotranspiration_thornthwaite_monthly.nc ← Thw PET (observed)
│   ├── 2m_temperature_detrended_monthly.nc                 ← Detrended T2m [from 2a]
│   ├── 2m_dewpoint_detrended_monthly.nc                    ← Detrended Tdew [from 2a]
│   ├── potential_evapotranspiration_PM_detrended.nc        ← PM PET (detrended T) [from 2b §5c]
│   └── potential_evapotranspiration_Thw_detrended.nc       ← Thw PET (detrended T) [from 2b §5c]
├── spi_results_seasonal/               ← SPI NetCDF results (by calendar month)
├── spei_results_seasonal/              ← SPEI (PM, observed) NetCDF results
├── spei_results_seasonal_thw/          ← SPEI (Thornthwaite, observed) NetCDF results
├── spei_results_detrended/             ← Counterfactual SPEI (PM, detrended T) [from 3, Run 3]
├── spei_results_detrended_thw/         ← Counterfactual SPEI (Thw, detrended T) [from 3, Run 4]
├── swei_results_seasonal/              ← SWEI NetCDF results
├── trend_analysis_pr_pet/              ← Pr/PET/T trend RDS + CSV outputs [from 4a]
├── temporal_drought/                   ← Trend test CSVs, basin plots, GIFs
│   ├── cache/
│   ├── basin_averaged_plots/
│   ├── basin_averaged_timeseries/
│   ├── event_ranking/
│   └── figures_temporal_diagnostics/
├── decomp_results/                     ← Decomposition and bias outputs [from 10a, 10b]
└── counterfactual_results/             ← SPEI comparison, warming signal outputs [from 11]
    ├── basin_timeseries_comparison/
    ├── warming_signal_maps/
    └── event_attribution/
```

---

## 3. Execution Order

Scripts must be run in the following order. The counterfactual branch is **embedded** within the existing scripts via the `RUN_DETRENDED_BRANCH` flag: `2a_detrend_temperature.R` must run first, then `2b_PET_ERALand.R` (which computes both observed and counterfactual PET in a single pass), then `3SPEI_ERALand.R` (which runs all four SPEI variants). Script 11 reads from the detrended output directories already produced by Step 3.

```
DROUGHT_ANALYSIS_utils.R   ← sourced automatically by all other scripts
        │
        ├─► 1SPI_ERALand.R
        │         Produces: spi_results_seasonal/
        │
        ├─► 2a_detrend_temperature.R
        │         Reads:    monthly_data_direct/2m_temperature_monthly.nc
        │                   monthly_data_direct/2m_dewpoint_temperature_monthly.nc
        │         Produces: monthly_data_direct/2m_temperature_detrended_monthly.nc
        │                   monthly_data_direct/2m_dewpoint_detrended_monthly.nc
        │                   monthly_data_direct/detrend_diagnostics/
        │
        ├─► 2b_PET_ERALand.R
        │         Reads:    all ERA5-Land fields + detrended T/Tdew (if RUN_DETRENDED_BRANCH=TRUE)
        │         Produces: potential_evapotranspiration_monthly.nc           (observed PM)
        │                   potential_evapotranspiration_thornthwaite_monthly.nc (observed Thw)
        │                   potential_evapotranspiration_PM_detrended.nc      [§5c, if flag=TRUE]
        │                   potential_evapotranspiration_Thw_detrended.nc     [§5c, if flag=TRUE]
        │
        ├─► 3SPEI_ERALand.R
        │         Reads:    Precip.nc + both PET pairs (observed + detrended)
        │         Produces: spei_results_seasonal/        (Run 1: SPEI_PM observed)
        │                   spei_results_seasonal_thw/    (Run 2: SPEI_Thw observed)
        │                   spei_results_detrended/       (Run 3: SPEI_PM counterfactual)
        │                   spei_results_detrended_thw/   (Run 4: SPEI_Thw counterfactual)
        │
        ├─► 4pr_pet_trends.r
        │         Reads:    Precip.nc, PET_PM.nc, PET_Thw.nc, T2m.nc
        │         Produces: trend_analysis_pr_pet/all_results.rds, analysis_metadata.rds, CSVs
        │       └─► 4pr_pet_trends_visualization.r
        │
        ├─► 5gif_animations.R
        │         Reads:    seasonal NetCDF from Steps 1–3
        │
        ├─► 6trend_test_ALL.R
        │         Reads:    seasonal NetCDF (SPI, SPEI_PM, SPEI_Thw)
        │         Produces: temporal_drought/{index}_{scale:02d}_results.csv
        │                   temporal_drought/{index}_{scale:02d}_basin_averaged_by_month.csv
        │       └─► 7basin_timeseries.R   (reads *_basin_averaged_by_month.csv)
        │       └─► 8trends_visualization.R (reads *_results.csv)
        │
        ├─► 9event_ranking.R
        │         Reads:    seasonal NetCDF directly (Steps 1–3)
        │
        ├─► 10a_dynamic_thermodynamic_decomp.R
        │         Reads:    spei_results_seasonal/, spei_results_seasonal_thw/
        │
        ├─► 10b_PET_bias_nonstationarity.R
        │         Reads:    WB_basin_average_monthly.csv (from Step 3),
        │                   PET summary CSVs (from Step 2b),
        │                   decomp_full_timeseries.csv (from Step 10a)
        │
        └─► 11_SPEI_counterfactual_comparison.R
                  Reads:    spei_results_seasonal/       (observed SPEI_PM)
                            spei_results_detrended/      (counterfactual SPEI_PM)
                            spei_results_seasonal_thw/   (observed SPEI_Thw)
                            spei_results_detrended_thw/  (counterfactual SPEI_Thw)
                  Produces: counterfactual_results/
```

---

## 4. Shared Utilities

**File:** `DROUGHT_ANALYSIS_utils.R`  
**Role:** Single source of truth for all project-wide constants, I/O functions, and statistical routines. All other scripts call `source("DROUGHT_ANALYSIS_utils.R")` as their first action.

### Section A — Project Configuration

Defines all directory paths and global constants. The working directory is controlled by the environment variable `NECHAKO_WD` (defaults to `D:/Nechako_Drought/Nechako/`):

| Constant | Value | Description |
|---|---|---|
| `DROUGHT_ONSET` | −0.5 | Event onset threshold (31st percentile of standard normal) |
| `DROUGHT_END` | −0.5 | Event termination threshold (same as onset; single-threshold detection) |
| `EQUAL_AREA_CRS` | EPSG:3005 | BC Albers projection for area-weighted calculations |
| `FIG_DPI` | 300 | PNG output resolution |
| `SPI_SCALES` | 1, 2, 3, 6, 12 | Accumulation timescales (months) |
| `SPEI_SCALES` | 1, 2, 3, 6, 12 | As above for SPEI |
| `SWEI_SCALE` | 3 | Single timescale for snow-water equivalent index |

**Note on thresholds:** `DROUGHT_ONSET = DROUGHT_END = −0.5`. The `detect_drought_events()` function (Section F) uses separate, configurable `onset_threshold` and `termination_threshold` arguments; when passed equal values the detection is single-threshold (no hysteresis gap). Severity classification inside `detect_drought_events()` uses hardcoded defaults of `severe_thr = −1.5` (D2 boundary) and `extreme_thr = −2.0` (D3 boundary), which are consistent with the ggplot band definitions in Section G.

### Section B — Basin-Averaged CSV I/O

`load_basin_avg_csv(index_type, scale)` — Loads the wide-format `Year × Month` CSV produced by Step 6 and returns a tidy `date / value` long-format data frame.

`save_basin_avg_from_pixels(ts_matrix, area_weights, dates, ...)` — Computes the area-weighted basin average from a pixel × month matrix and writes the authoritative CSV:

```
basin_avg[t] = Σ( index_value[pixel,t] × cell_area[pixel] )
             / Σ( cell_area[pixel] )
```

Cell areas are derived from `terra::cellSize()` in the equal-area CRS (EPSG:3005), ensuring correct weighting despite the latitude-dependent distortion of geographic coordinates.

### Section C & D — NetCDF File Finders and Date Parsing

`find_seasonal_nc_files()` locates files matching `{index}_{scale:02d}_month{mm}_{MonName}.nc`. `extract_dates_from_nc()` extracts the time axis from the NetCDF metadata, falling back to a reconstructed monthly sequence from January 1950 if metadata is missing.

### Section E — Basin Spatial Helpers

`precompute_basin_geometry()` projects the basin KMZ/KML boundary to EPSG:3005, computes per-cell effective areas (including partial border cells via `terra::rasterize(..., cover=TRUE)`), and returns all indices needed for fast area-weighted averaging.

### Section F — Drought Event Detection

`detect_drought_events(df)` implements the threshold-level run method (Yevjevich 1967) with configurable `onset_threshold`, `termination_threshold`, and `min_duration`. Events are classified by minimum index value:

| Class | Condition |
|---|---|
| D1: moderately dry | min_value > −1.5 |
| D2: severely dry | −2.0 < min_value ≤ −1.5 |
| D3: extremely dry | min_value ≤ −2.0 |

### Section G — ggplot Helpers

Provides standardized colour-coded drought threshold bands (`drought_band_layers()`), a shared time series theme, and `save_figure()` which simultaneously writes `.pdf` and `.png` versions of any ggplot object. Band boundaries: D3 (< −2.0), D2 (−2.0 to −1.5), D1 (−1.5 to −1.0), N0 normal (−1.0 to 1.0), W1–W3 wet.

### Section H — Excel Summary Export

`export_summary_excel()` compiles descriptive statistics (mean, median, SD, min, max) and drought event counts for the 2020–2025 window across all index × timescale combinations into a formatted Excel workbook with bold header styling.

### Sections J1–J7 — Shared Statistical Functions

The core statistical engine, organized into subsections:

**J1 — Primitive helpers:** `check_min_value_threshold()` (flags series with too many near-zero values) and `calculate_variance_with_ties()` (tie-corrected Mann-Kendall variance).

**J2 — PELT Changepoint:** `detect_regime_shift_pelt()` applies the Pruned Exact Linear Time algorithm (Killick et al. 2012) with BIC penalty to detect the first mean-shift in a univariate series.

**J3–J4 — Mann-Kendall variants (internal helpers):** `.mk_vc_single()` (Hamed-Rao 1998 variance correction) and `.mk_tfpw_single()` (Trend-Free Pre-Whitening, Yue et al. 2002). Both return tau, p-value, Sen's slope, and lag-1 autocorrelation.

**J5 — Spectral peak detection (internal helper):** `.spectral_single()` tests for significant low-frequency peaks against an AR(1) red-noise null (periods ≥ 24 months), optionally using a Monte Carlo white-noise envelope.

**J6 — Omnibus single-series wrapper:** `mk_tfpw_spectral_for_series()` chains all four methods (VC-MK, TFPW-MK, spectral, PELT) with quality control filters for short series, excessive ties, or near-zero variance. Called directly by `4pr_pet_trends.r`.

**J7 — Vectorized pixel-matrix wrappers:** `vectorized_mann_kendall()`, `vectorized_mann_kendall_tfpw()`, `vectorized_regime_shift_pelt()`, and `vectorized_spectral_peaks()` apply the single-series routines across a pixel × time matrix using a cross-platform parallel backend (Windows PSOCK cluster or Unix fork). Called directly by `6trend_test_ALL.R`.

---

## 5. Script 1 — SPI Calculation

**File:** `1SPI_ERALand.R`  
**Input:** `monthly_data_direct/total_precipitation_monthly.nc`  
**Outputs:** `spi_results_seasonal/spi_{sc:02d}_month{mm:02d}_{MonName}.{nc,csv,xlsx}`

### Method: Variance-Aware SPI

The SPI is computed independently for each of the 5 timescales (1, 2, 3, 6, 12 months) and for each of the 12 calendar months. The algorithm selects among three fitting methods based on the statistical properties of each pixel-month sample:

| Case | Condition | Method |
|---|---|---|
| 1 | Variance = 0 exactly | SPI = 0 (neutral) |
| 2 | Variance < 0.01 mm² | Empirical probability ranking (Blom plotting position) |
| 3 | Sufficient variance (≥ 0.01 mm²) | Gamma distribution fitting via L-moments (`lmomco`) |

For Case 3, zero-precipitation months are handled explicitly: the CDF is the mixed distribution P(X ≤ x) = p₀ + (1 − p₀) × Fγ(x), where p₀ is the fraction of zero months and Fγ is the gamma CDF fitted to positive months only. SPI values are clamped to [−4.75, +4.75].

### Preprocessing

Precipitation is read in ERA5-Land units (m/day mean rate), reprojected to EPSG:3005, masked to the basin boundary, and converted to monthly totals (mm/month) by multiplying by the exact number of days in each calendar month. A unit verification diagnostic confirms whether the raw data is in m/day, mm/day, or mm/month before conversion.

### Parallel Implementation

Processing uses `parallel::parLapply()` across all basin pixels simultaneously. The cluster is reused across timescales to minimize setup overhead.

### Outputs per Timescale

- NetCDF files (one per calendar month): pixel-level SPI time series
- CSV files: same data in tabular format (pixels × years)
- Excel workbook: all months in one file
- Distribution method map (PNG): spatial coverage of Gamma vs. Empirical methods
- Combined summary text file: NA rates, method percentages, drought frequency statistics

---

## 6. Script 2a — Month-Specific Temperature Detrending

**File:** `2a_detrend_temperature.R`  
**Inputs:** `monthly_data_direct/2m_temperature_monthly.nc`, `2m_dewpoint_temperature_monthly.nc`  
**Outputs:**
- `monthly_data_direct/2m_temperature_detrended_monthly.nc`
- `monthly_data_direct/2m_dewpoint_detrended_monthly.nc`
- `monthly_data_direct/detrend_diagnostics/` — per-month slope maps, residual QC plots

**Must run before `2b_PET_ERALand.R`.**

### Method: Month-Specific Linear Detrending

For each calendar month m ∈ {1 … 12} and each basin pixel, an OLS linear regression of temperature against time is fitted over the full 1950–2025 record:

```
T(year, m) = a_m + b_m × year + ε_m
```

The detrended series is:

```
T_det(year, m) = T(year, m) − b_m × year
```

which retains the long-term mean ā_m while removing only the linear warming signal. The identical procedure is applied to the 2 m dewpoint temperature Tdew so that:

```
Tdew_det(year, m) ≤ T_det(year, m)   for all years and months
```

This inequality guarantees RH < 100% and VPD > 0, which are physical requirements for the Penman-Monteith equation. Any residual instances where Tdew_det > T_det (possible at pixels where OLS fits cross) are clamped to T_det as a physics guard. All outputs are in Kelvin (same unit as ERA5-Land inputs) so that `2b_PET_ERALand.R` can load them without unit conversion.

### Diagnostic Outputs

- Spatial maps of b_m (°C/year warming rate) per calendar month
- Basin-averaged observed vs. detrended T2m time series — 12-panel figure
- Distribution of OLS R² per month (confirms linearity assumption)
- Residual lag-1 ACF check (verifies approximately white-noise residuals)
- Physics check: VPD before and after detrending
- `detrend_summary.csv`: per-month basin-mean statistics

---

## 7. Script 2b — PET Calculation

**File:** `2b_PET_ERALand.R`  
**Inputs:** ERA5-Land monthly NetCDF fields — T2m, T2m dew point, 10m winds (u, v), surface solar radiation downwelling (SSRD), surface pressure, geopotential; plus detrended T2m and Tdew from `2a_detrend_temperature.R` (if `RUN_DETRENDED_BRANCH = TRUE`)  
**Outputs (all to `monthly_data_direct/`):**
- `potential_evapotranspiration_monthly.nc` — Penman-Monteith PET (mm/day), observed
- `potential_evapotranspiration_thornthwaite_monthly.nc` — Thornthwaite PET (mm/day), observed
- `potential_evapotranspiration_PM_detrended.nc` — PM PET from detrended T (mm/day) [Section 5c]
- `potential_evapotranspiration_Thw_detrended.nc` — Thw PET from detrended T (mm/day) [Section 5c]
- Diagnostic PDFs and CSVs in `pet_calcs/diagnostics/`

This script computes **both observed and counterfactual PET in a single run**, controlled by:

```r
RUN_DETRENDED_BRANCH <- TRUE   # set FALSE to skip counterfactual PET (Sections 1–6 only)
```

### Method A: FAO-56 Penman-Monteith (Sections 1–5a)

Implements the full FAO-56 reference evapotranspiration equation (Allen et al. 1998), vectorized over the entire spatial grid using `terra` raster arithmetic:

```
ET₀ = [0.408 × Δ(Rn−G) + γ(900/(T+273))u₂(es−ea)] / [Δ + γ(1+0.34u₂)]
```

Key implementation details:

- **Rs conversion:** ERA5-Land monthly-averaged SSRD is in J/m²/day (mean daily); divided by 1×10⁶ to convert to MJ/m²/day. The days-in-month factor is **not** applied again (already a daily mean rate).
- **Soil heat flux G:** FAO-56 Eq. 45 (G = 0.14 × (Tᵢ − Tᵢ₋₁)) with a cold-region constraint capping |G| at 0.3 × |Rn| when T < 0°C, preventing unphysical negative net radiation under frozen soil.
- **Wind speed:** 10 m components converted to 2 m height using the logarithmic wind profile.
- **Extraterrestrial radiation Ra:** computed per pixel per calendar month from latitude-dependent sunrise hour angles and solar declination at mid-month Julian day.

### Method B: Thornthwaite / Temperature-Only (Section 5b)

A pure thermodynamic PET estimate requiring only 2m air temperature. Used as the thermodynamic baseline for the Dynamic vs. Thermodynamic decomposition (Script 10a):

```
PET_thw = 16 × (10T/I)^a × (N/12) × (d/30)   [mm/month]
        = 0  when T ≤ 0°C
```

where I is the annual heat index summed from 1991–2020 climatological monthly temperatures, a is a polynomial function of I, and N is daylight hours derived from the same Ra geometry used for the PM calculation.

### Section 5c — Counterfactual PET (detrended T, `RUN_DETRENDED_BRANCH = TRUE`)

Structurally identical to Sections 5a and 5b, with only the temperature and dewpoint inputs replaced by the detrended fields from `2a_detrend_temperature.R`. Radiation, wind, and pressure remain at observed values, isolating only the temperature-mediated warming effect on PET.

**PM (detrended):**

```
ET₀_det = [0.408 × Δ_det(Rn − G_det) + γ(900/(T_det+273))u₂(es_det − ea_det)]
        / [Δ_det + γ(1 + 0.34 u₂)]
```

**Thornthwaite (detrended):** The annual heat index I and exponent a are recomputed from the **detrended** 1991–2020 climatological monthly means, ensuring the Thornthwaite calibration also reflects the world without warming. Daylight hours N are unchanged (latitude-dependent, not temperature-dependent).

### Diagnostic Outputs

Nine diagnostic plots: PM time series, seasonal cycle boxplot, inter-annual trend with Sen's slope, PM vs. Thornthwaite overlays (time series, seasonal climatology in mm/month, monthly bias bar chart, annual scatter), spatial sample maps, and a QC summary table. A side-by-side seasonal climatology of observed vs. detrended PET is produced when the counterfactual branch runs.

---

## 8. Script 3 — SPEI Calculation

**File:** `3SPEI_ERALand.R`  
**Inputs:** Precipitation NetCDF + both observed and detrended PET NetCDF files from Script 2b  
**Outputs:**
- `spei_results_seasonal/` — SPEI_PM observed (Run 1)
- `spei_results_seasonal_thw/` — SPEI_Thw observed (Run 2)
- `spei_results_detrended/` — SPEI_PM counterfactual (Run 3, if `RUN_DETRENDED_BRANCH = TRUE`)
- `spei_results_detrended_thw/` — SPEI_Thw counterfactual (Run 4, if `RUN_DETRENDED_BRANCH = TRUE`)

### Method

The SPEI is computed from the climatic water balance:

```
WB = Precipitation (mm/month) − PET (mm/month)
```

This signed quantity is accumulated over the chosen timescale and then standardized using a three-distribution parametric fitting procedure (following Kao & Govindaraju 2010):

| Case | Condition | Action |
|---|---|---|
| 1 | Zero variance | SPEI = 0 |
| 2 | Sufficient variance | L-moments fitting; best of GLO, PE3, GEV by KS p-value |
| 3 | L-moment/fitting failure | Result = NA |

### Four-Run Architecture

The script runs the identical SPEI algorithm via a shared `run_spei_loop()` function, controlled by `RUN_DETRENDED_BRANCH`:

- **Run 1 (SPEI_PM):** WB = P − PET_PM → `spei_results_seasonal/`
- **Run 2 (SPEI_Thw):** WB = P − PET_Thw → `spei_results_seasonal_thw/`
- **Run 3 (SPEI_PM_det):** WB = P − PET_PM_detrended → `spei_results_detrended/`
- **Run 4 (SPEI_Thw_det):** WB = P − PET_Thw_detrended → `spei_results_detrended_thw/`

Runs 3 & 4 require `RUN_DETRENDED_BRANCH = TRUE` and the detrended PET files from `2b_PET_ERALand.R` Section 5c. Precipitation is kept at observed values in all four runs; only the evaporative demand changes. The standardization calibration is performed over the full 1950–2025 record of each respective water balance, ensuring SPEI_observed and SPEI_detrended are expressed on the same standardized scale.

### Parallel Cluster Fix

`environment(worker_fn) <- globalenv()` is applied to the parallel worker function, preventing `parLapply` from serializing the ~6 MB `wb_input` matrix in the function closure on every per-scale export and eliminating "error writing to connection" failures on scale 2+.

---

## 9. Script 4a — Pr/PET/Temperature Trend Analysis

**File:** `4pr_pet_trends.r`  
**Inputs:** Precipitation, PET (PM), PET (Thornthwaite), T2m NetCDF files  
**Outputs:** `trend_analysis_pr_pet/all_results.rds`, `analysis_metadata.rds`, `clipped_template.rds`, summary CSVs

### Variables Analyzed

| Variable | Aggregation | Notes |
|---|---|---|
| Precipitation | Annual sum (mm/year) | Full year, all 12 months |
| PET_PM | Warm-season sum Apr–Oct (mm) | Penman-Monteith |
| PET_Thw | Warm-season sum Apr–Oct (mm) | Thornthwaite |
| Temperature | Annual mean (°C) + 12 monthly means | 2 m air temperature |

### Analysis per Pixel

For each variable and each basin pixel, the omnibus `mk_tfpw_spectral_for_series()` wrapper (from utils Section J6) is called, producing:

- **VC-MK:** Kendall τ, H-R98 corrected p-value, Sen's slope, lag-1 autocorrelation
- **TFPW-MK:** Same statistics on the pre-whitened series
- **Spectral:** Number of significant low-frequency peaks (Monte Carlo white-noise envelope), dominant period
- **PELT:** Regime shift year, mean before/after, magnitude of shift

### Specific-Point Extraction

Trend analysis is also performed for three user-defined geographic points, producing per-point statistics for all variables at both annual and monthly timescales, written to `point_trend_stats.csv` and `point_monthly_timeseries.csv`.

### Quality Filters

Series are flagged and excluded if: fewer than 10 non-NA observations; excessive zero values (> 50% for PET warm-season monthly series); or excessive ties (> 50% repeated values).

---

## 10. Script 4b — Trend Visualization

**File:** `4pr_pet_trends_visualization.r`  
**Inputs:** `.rds` files from Script 4a  
**Outputs:** Basin-wide maps, time series panels, PET comparison plots

### Publication Figures

- **Basin-wide trend maps:** Kendall τ and Sen's slope for VC and TFPW, for all four variables
- **Regime shift maps:** Decade of first PELT changepoint per pixel
- **Spectral analysis maps:** Dominant periodicity spatial distribution
- **12-panel calendar-month series (Functions 13a–13d):** For each variable, inter-annual scatter for each calendar month with OLS trend line and significance annotation
- **Four-variable Sen's slope profile (Function 13e):** Compares all variables' monthly trend magnitudes in one figure
- **PET_PM vs. PET_Thw comparison (Function 14a–14e):** Bias time series, per-month OLS comparison, spatial trend maps side-by-side, bias bar chart, annual scatter
- **Specific-point time series:** 4-panel plots (Pr / PET_PM / PET_Thw / T) per point
- **Method comparison:** VC vs. TFPW scatter plots and agreement bars

---

## 11. Script 5 — Animated GIFs

**File:** `5gif_animations.R`  
**Inputs:** Seasonal NetCDF files from Scripts 1–3  
**Outputs:** Animated GIFs in `temporal_drought/drought_gifs/`

Creates frame-by-frame animations of spatial drought index evolution from 1981 onward. Each frame is rendered at 900×720 px using `ragg::agg_png()` (faster than base `png()`), then assembled into a GIF at 4 frames per second using `gifski`. The colour scale uses the `scico` palette with the drought threshold bands from the shared utilities. GIFs are produced for all configured index × timescale combinations.

---

## 12. Script 6 — Spatial Trend Test for All Indices

**File:** `6trend_test_ALL.R`  
**Inputs:** Seasonal NetCDF files from Scripts 1–3 (SPI, SPEI_PM, SPEI_Thw)  
**Outputs:** `temporal_drought/{index}_{scale:02d}_results.csv` and `*_basin_averaged_by_month.csv`

**Indices processed:** SPI, SPEI_PM, SPEI_Thw at scales 1, 2, 3, 6, 12 months. SWEI is not processed by this script.

### Pixel-Level Statistics

For each index and each timescale, the script assembles a pixel × month matrix and computes in parallel:

- **VC-MK (H-R98):** tau_vc, p_value_vc, sl_vc, plus BH-FDR corrected p_fdr_vc across all valid basin pixels
- **TFPW-MK:** tau_tfpw, p_value_tfpw, sl_tfpw, p_fdr_tfpw
- **PELT:** regime_shift_year, regime_shift_detected, n_changepoints, mean_before/after, magnitude_shift
- **Spectral:** n_spectral_peaks (AR(1) red-noise test, periods ≥ 24 months)
- **Wald-Wolfowitz runs test:** clustering p-value on binarized drought event sequences

### Two Drought Event Detection Methods

**Method 1 — Sheffield & Wood (2008) single threshold:** Onset and termination both at `DROUGHT_ONSET` (−0.5, 31st percentile). No minimum duration. Intensity = mean(onset_thr − xₜ) over event months. Produces columns: n_events, mean_duration, max_intensity, and duration-class breakdowns D3-6, D7-12, D12+.

**Method 2 — Scale-specific minimum duration:** Same onset and termination threshold as Method 1 (−0.5). Differs from Method 1 by applying scale-specific minimum event duration filters and using `end_thr` as the intensity reference. Scale-specific minimum durations:

| Scale | Min Duration |
|---|---|
| 1 | 2 months |
| 3 | 3 months |
| 6 | 4 months |
| 12 | 6 months |

All Method 2 columns are suffixed `_hyst`.

### Area-Weighted Basin Average

After the pixel loop, `save_basin_avg_from_pixels()` from the shared utilities computes the authoritative area-weighted basin average and writes the `*_basin_averaged_by_month.csv` files consumed by Script 7.

### Single-File Mode

Multivariate indices (MSPI, MSPEI) store all months in one NetCDF file rather than 12 monthly files. The matrix layout is auto-detected.

---

## 13. Script 7 — Basin-Averaged Time Series

**File:** `7basin_timeseries.R`  
**Inputs:** `*_basin_averaged_by_month.csv` files from Script 6; PET summary CSVs from Script 2b  
**Outputs:** Time series PNGs, combined CSV, Excel summary with correlation matrix

### Processing Steps

**Part 1:** Loads all authoritative basin-averaged CSVs for SPI, SPEI (PM and Thw), SWEI, and optionally MSPI/MSPEI.

**Part 2:** Re-exports each index as a flat `date / value` CSV to `basin_averaged_timeseries/`.

**Part 3:** Builds `basin_mean_all_indices.csv`, a combined long-format file covering all indices.

**Part 4:** Produces one publication-quality PNG per index × timescale, with drought threshold bands, event shading, and climatological annotations.

**Part 4d (Seasonality figure):** Plots mean monthly precipitation vs. mean PET (PM) as a 12-bar climatological water balance figure, requiring `ERA5Land_Nechako_PET_monthly_summary.csv` from `2b_PET_ERALand.R`.

**Part 4e (Manuscript figures):** Generates a 3-panel SPI-1/2/3 time series (Fig 4e) and a 3×3 panel figure combining SPI, SPEI_PM, and SPEI_Thw at scales 1, 2, and 3 (Figure 2 in Results Section 4.2).

**Part 5:** Writes `Drought_Summary.xlsx` containing descriptive statistics, 2020–2025 drought event counts, and a pairwise Pearson correlation matrix.

---

## 14. Script 8 — Spatial Trends Visualization

**File:** `8trends_visualization.R`  
**Inputs:** `{index}_{scale:02d}_results.csv` from Script 6  
**Outputs:** Publication maps in `temporal_drought/figures_temporal_diagnostics/`

### Figure Suite

**Figure 1 — Trend maps:** Spatial distribution of Kendall τ_vc with FDR-significant pixels highlighted.

**Figure 2 — Event maps:** Three-panel spatial maps showing number of events, mean event duration, and maximum intensity.

**Figure 3 — Temporal pattern maps:** Spatial distribution of drought event clustering (runs-test p-value) and regime shift decade (from PELT changepoints).

**Figure 4 — Timescale comparison:** Summary statistics across timescales: percentage of significant pixels, median τ, mean event count, and mean duration.

**Figure 5 — Method comparison:** VC vs. TFPW scatter plots and agreement bar charts.

All spatial plots use the basin KMZ boundary as an overlay and are rendered at 300 DPI.

---

## 15. Script 9 — Drought Event Ranking

**File:** `9event_ranking.R`  
**Inputs:** Seasonal NetCDF files from Scripts 1–3 (loaded directly)  
**Outputs:** `temporal_drought/event_ranking/`

### SW08 Event Catalog

Script 9 builds a fresh event catalog directly from raw NetCDF files (independent of Script 6), applying the SW08 single-threshold method. The script defines its own local constant:

```r
DROUGHT_ONSET_EVENT <- -0.5   # onset = termination; consistent with DROUGHT_ONSET in utils
```

- Onset: index < −0.5 (31st percentile)
- Termination: index ≥ −0.5 (same threshold — no hysteresis)
- Intensity: I = mean(onset_thr − xₜ) over event months [always ≥ 0]
- Severity: S = I × Duration
- Min duration: none (all events ≥ 1 month retained)
- Duration classes: Short (3–6 months), Medium (7–12 months), Long (≥ 12 months)

### Return Period Estimation

**Classical empirical:** Weibull plotting position for each duration class separately.

**Kendall copula (bivariate):** Joint distribution of Duration × Severity fitted using a parametric copula (from the `copula` R package), enabling return period estimation that accounts for the dependence structure between event dimensions. Marginal distributions are fitted via MLE using `fitdistrplus`.

### Outputs

- `ranked_event_catalog_SW08.csv`: Complete ranked event list
- `return_period_by_duration_class_SW08.csv`: Return period tables per duration class
- `Table_SW08_vs_Hyst_Comparison.csv`: Side-by-side comparison of SW08 vs. hysteresis statistics
- Manuscript Figure 3: Drought severity ranking visualization
- Fig 4c: Copula diagnostic plots
- Fig 4d: Multi-index Kendall return period comparison

---

## 16. Script 10a — Dynamic/Thermodynamic Decomposition

**File:** `10a_dynamic_thermodynamic_decomp.R`  
**Inputs:** `spei_results_seasonal/` (SPEI_PM) and `spei_results_seasonal_thw/` (SPEI_Thw) from Script 3  
**Outputs:** `decomp_results/`

### Decomposition Framework

```
SPEI_PM  = Total drought severity (precipitation + radiation + wind + temperature)
SPEI_Thw = Thermodynamic component (temperature-only via Thornthwaite PET)
Dynamic  ≈ SPEI_PM − SPEI_Thw
```

The thermodynamic fraction of drought intensity:

```
F_thm = |SPEI_Thw| / |SPEI_PM|    when SPEI_PM ≤ THR_FTHM (−0.5)
```

`THR_FTHM = −0.5` is consistent with `DROUGHT_ONSET`, preventing unstable fractions when the denominator is near zero during non-drought periods. Decomposition scales: 1, 2, 3, 6, 12 months.

---

## 17. Script 10b — PET Bias Non-Stationarity

**File:** `10b_PET_bias_nonstationarity.R`  
**Inputs:** `WB_basin_average_monthly.csv` from Script 3; PET summary CSVs from Script 2b; `decomp_full_timeseries.csv` from Script 10a  
**Outputs:** `decomp_results/bias_nonstationarity/`

### Correction Methodology

```
δ_bias_k(m)  = mean(bias_k, P3) − mean(bias_k, P1)
δ_SPEI_k(m)  = δ_bias_k(m) / SD(WB_k, P1, month m)
SPEI_Thw_corr = SPEI_Thw + δ_SPEI_k(m)
F_thm_corr    = |SPEI_Thw_corr| / (|SPEI_PM| + ε)
```

where P1 = 1950–1990 (pre-warming reference, 41 years) and P3 = 2022–2025 (recent event window, 4 years). `bias_mm_month = (PET_Thw − PET_PM) × days_in_month` [mm/month].

### Outputs

- `bias_monthly_timeseries.csv`
- `correction_params_by_scale.csv`: δ_bias_k, SD_WB_k, δ_SPEI_k per scale × calendar month
- `SPEI_Thw_corrected_P3.csv`
- `bias_manuscript_tables.xlsx`
- Multi-panel diagnostic PDF

---

## 18. Script 11 — Counterfactual SPEI Comparison

**File:** `11_SPEI_counterfactual_comparison.R`  
**Inputs:**
- `spei_results_seasonal/` — observed SPEI_PM (Run 1, from Script 3)
- `spei_results_detrended/` — counterfactual SPEI_PM (Run 3, from Script 3)
- `spei_results_seasonal_thw/` — observed SPEI_Thw (Run 2, from Script 3)
- `spei_results_detrended_thw/` — counterfactual SPEI_Thw (Run 4, from Script 3)  
**Outputs:** `counterfactual_results/`

### Scientific Context

The observed SPEI embeds the full climate change signal — precipitation deficit plus warming-driven increase in evaporative demand. This script quantifies: **how much more severe was the 2022–2025 drought because of the long-term warming trend?**

### Core Quantity: Warming-Induced SPEI Change (ΔSPEI)

```
ΔSPEI_PM(t)  = SPEI_PM_observed(t)  − SPEI_PM_detrended(t)
ΔSPEI_Thw(t) = SPEI_Thw_observed(t) − SPEI_Thw_detrended(t)
```

A negative ΔSPEI indicates the warming trend **deepened** the drought (observed drought is more severe than the counterfactual). A near-zero ΔSPEI indicates the drought was driven primarily by precipitation deficit, independent of warming.

### Analysis Components

**Part 1 — Basin time series comparison:** Overlaid observed and counterfactual SPEI basin averages for scales 1, 3, 6, and 12, with ΔSPEI shading. The 2022–2025 window is highlighted.

**Part 2 — Drought event attribution:** For every event in the SW08 catalog (Script 9), the mean ΔSPEI over the event duration is computed, quantifying how many SPEI units of severity were attributable to warming vs. precipitation deficit. A stacked bar chart summarizes event-level attribution.

**Part 3 — Spatial warming signal maps:** Pixel-level mean ΔSPEI for the 2022–2025 period at scales 3, 6, and 12.

**Part 4 — Warming fraction index:**

```
F_warm(t) = |ΔSPEI_PM(t)| / |SPEI_PM_observed(t)|    when SPEI_PM ≤ −0.5
```

Provides a pixel-by-pixel, time-varying measure of what fraction of observed drought severity would not have occurred without the warming trend. Gate threshold −0.5 is consistent with `THR_FTHM` in Script 10a.

**Part 5 — Cross-comparison with Script 10a/10b:** F_warm (detrended counterfactual) vs. F_thm (Thornthwaite decomposition), validating consistency between the two attribution approaches.

### Key Outputs

| File | Content |
|---|---|
| `basin_timeseries_comparison/SPEI_obs_vs_det_{scale}.png` | Time series overlay per scale |
| `basin_timeseries_comparison/delta_SPEI_summary.csv` | Monthly ΔSPEI table (all scales) |
| `warming_signal_maps/delta_SPEI_2022_2025_{scale}.nc` | Pixel-level ΔSPEI NetCDF |
| `warming_signal_maps/Fwarm_2022_2025_{scale}.png` | Warming fraction maps |
| `event_attribution/event_warming_attribution.csv` | Per-event warming contribution |
| `event_attribution/Fig_event_attribution_stacked.pdf` | Stacked attribution figure |
| `Fwarm_vs_Fthm_comparison.pdf` | Cross-validation with Script 10a |

---

## 19. Manuscript Support Script

**File:** `11nechako_drought_manuscript1.R`

This script is a merged manuscript support tool — **separate from the counterfactual comparison (Script 11)** — combining:

1. **Table 3:** Thermodynamic fraction sensitivity (original vs. bias-corrected F_thm) by month and SPEI scale for 2022–2025 drought months
2. **Section 4.3 SW08 event catalog:** Manuscript blank-filling extractor for all event statistics
3. **MS1 blank extractor (`w_extract_ms1_blanks.R`):** Compiles all in-text statistics from upstream scripts
4. **Fig 2 / Seasonality:** 12-month climatological P vs. PET two-panel figure

**Prerequisites:** Must be run after Scripts 1, 3, 10a, 10b, 9, and 7 (and `2b_PET_ERALand.R`).

---

## 20. Statistical Methods

### Month-Specific Linear Detrending (Scripts 2a, 2b §5c, 3 Runs 3&4, 11)

Each calendar month's time series at each pixel is detrended independently using OLS:

```
T_det(y, m) = T(y, m) − b̂_m × y
```

where b̂_m is the OLS slope for month m. The intercept is implicitly retained through the long-term mean, so the detrended series has the same temporal mean as the observed series. Co-detrending of T2m and Tdew independently ensures VPD > 0 is maintained, confirmed by a post-hoc check.

### Mann-Kendall Trend Tests

All trend tests are applied to annual or seasonal aggregates (not monthly raw data).

**Variance-Corrected MK (Hamed & Rao 1998):** The variance of the MK S-statistic is inflated by a factor derived from lag-1 through lag-(n-1) autocorrelations of the detrended series. Implemented via `modifiedmk::mmkh()` in all three code paths (`.mk_vc_single()`, `vectorized_mann_kendall()`, and `mk_tfpw_spectral_for_series()`).

**TFPW-MK (Yue et al. 2002):** The series is detrended with Sen's slope, lag-1 autocorrelation r₁ is estimated, the residual is pre-whitened (Xₜ − r₁Xₜ₋₁), and the trend is blended back in before the MK test. **Sen's slope is computed on the original series** (not the pre-whitened one) for physically interpretable trend magnitudes.

**FDR Correction:** Benjamini-Hochberg False Discovery Rate adjustment is applied across all valid basin pixels within each index × timescale run. Use `p_fdr_vc` for significance mapping.

### Sen's Slope

Computed via `trend::sens.slope()`.

### PELT Changepoint Detection

Applied to **annual averages** of the monthly pixel series. The Pruned Exact Linear Time algorithm (Killick et al. 2012) with BIC penalty detects the optimal number and position of mean-shift changepoints.

### Spectral Analysis

Two implementations:

1. **AR(1) red-noise (pixel-matrix wrappers in utils J7, Script 6):** Periodogram tested against a scaled AR(1) background spectrum; peaks at periods ≥ 24 months identified at 95% confidence.
2. **Monte Carlo white-noise (single-series wrapper J6, Script 4a):** Maximum periodogram amplitude compared against the distribution of maxima from 500 random white-noise realizations.

### Copula-Based Return Periods (Script 9)

Bivariate parametric copula fitted to (Duration, Severity) pairs. Kendall return period:

```
T_K = n / (# events with KC ≥ kc_obs)
```

---

## 21. Key Parameters and Thresholds

| Parameter | Value | Used In | Description |
|---|---|---|---|
| `DROUGHT_ONSET` | −0.5 | All scripts | Event onset threshold (31st percentile of standard normal) |
| `DROUGHT_END` | −0.5 | Scripts 6, 7, 8 | Event termination threshold (same as onset) |
| `THR_FTHM` | −0.5 | Script 10a | Gate for thermodynamic fraction F_thm |
| `THR_EVENT` | −0.5 | Script 10a | Event identification threshold |
| `DROUGHT_ONSET_EVENT` | −0.5 | Script 9 | SW08 onset for event ranking (local constant, same value as `DROUGHT_ONSET`) |
| `THR_FWARM` | −0.5 | Script 11 | Gate for warming fraction F_warm (same convention as F_thm) |
| Severity boundary D2 | −1.5 | detect_drought_events, ggplot bands | Severely dry classification threshold |
| Severity boundary D3 | −2.0 | detect_drought_events, ggplot bands | Extremely dry classification threshold |
| `EQUAL_AREA_CRS` | EPSG:3005 | All scripts | BC Albers projection |
| Detrending method | Month-specific OLS | Scripts 2a, 2b §5c | Linear fit per calendar month per pixel |
| Detrending scope | T2m + Tdew independently | Script 2a | Co-detrending to preserve VPD > 0 |
| `alpha` | 0.05 | Scripts 4a, 6 | Significance level for all tests |
| `max_tie_pct` | 50% | Scripts 4a, 6 | Maximum tie fraction before filtering |
| SPI/SPEI timescales | 1, 2, 3, 6, 12 | Scripts 1–3, 6 | Accumulation periods (months) |
| SWEI timescale | 3 | Utils | Single timescale for SWE index |
| `warm_months_pet` | 4–10 (Apr–Oct) | Script 4a | Warm season for PET aggregation |
| Climatological period | 1991–2020 | Script 2b | Reference period for Thornthwaite I and a |
| P1 reference | 1950–1990 | Script 10b | Pre-warming baseline (41 years) |
| P3 event window | 2022–2025 | Scripts 9, 10a, 10b, 11 | Recent drought focus period (4 years) |
| `FIG_DPI` | 300 | All scripts | Figure output resolution |
| `RUN_DETRENDED_BRANCH` | TRUE (default) | Scripts 2b, 3 | Flag to execute counterfactual PET and SPEI computations |

---

## 22. Data Flow Diagram

```
ERA5-Land NetCDF (T2m, Tdew, u10, v10, SSRD, SP, Z, Precip)
         │
         ├──► 1SPI_ERALand.R
         │         spi_results_seasonal/
         │
         ├──► 2a_detrend_temperature.R
         │         2m_temperature_detrended_monthly.nc
         │         2m_dewpoint_detrended_monthly.nc
         │
         ├──► 2b_PET_ERALand.R  (observed §1-5b + counterfactual §5c)
         │         PET_PM.nc         PET_Thw.nc
         │         PET_PM_det.nc     PET_Thw_det.nc      [if RUN_DETRENDED_BRANCH]
         │
         └──► 3SPEI_ERALand.R  (4 runs)
                   Run 1: spei_results_seasonal/        (SPEI_PM observed)
                   Run 2: spei_results_seasonal_thw/    (SPEI_Thw observed)
                   Run 3: spei_results_detrended/       (SPEI_PM counterfactual)
                   Run 4: spei_results_detrended_thw/   (SPEI_Thw counterfactual)
                       │
                       ├──► 11_SPEI_counterfactual_comparison.R
                       │         ΔSPEI maps, F_warm, event attribution
                       │         counterfactual_results/
                       │
                       ├──► 6trend_test_ALL.R  (SPI, SPEI_PM, SPEI_Thw)
                       │         {index}_{scale}_results.csv
                       │         {index}_{scale}_basin_averaged_by_month.csv
                       │       ├──► 7basin_timeseries.R  (time series PNG, Excel)
                       │       └──► 8trends_visualization.R  (spatial maps)
                       │
                       ├──► 5gif_animations.R  (animated GIFs)
                       │
                       ├──► 9event_ranking.R  (SW08 catalog, copula RPs)
                       │
                       ├──► 10a_dynamic_thermodynamic_decomp.R
                       │         F_thm, decomp_results/
                       │       └──► 10b_PET_bias_nonstationarity.R
                       │                 bias correction, decomp_results/bias_nonstationarity/
                       │
4pr_pet_trends.r ──────────────────────────────────► 4pr_pet_trends_visualization.r
(Pr, PET, T trend RDS)                               (trend maps, monthly panels)

11nechako_drought_manuscript1.R
(Table 3, SW08 blanks, seasonality figure — reads outputs from multiple upstream scripts)
```

---

## 23. R Package Dependencies

| Package | Version | Purpose |
|---|---|---|
| `terra` | ≥ 1.7 | Raster I/O, reprojection, masking, cell-size computation |
| `ncdf4` | any | Low-level NetCDF metadata access (fallback time parsing) |
| `data.table` | ≥ 1.14 | Fast CSV I/O, memory-efficient pixel tables |
| `sf` | ≥ 1.0 | Basin boundary vector operations |
| `ggplot2` | ≥ 3.4 | Publication-quality figures |
| `patchwork` | any | Multi-panel figure composition |
| `lubridate` | any | Date arithmetic |
| `zoo` | any | Rolling window aggregation (`rollapply`) |
| `lmomco` | any | L-moment distribution fitting (Gamma, GLO, PE3, GEV) |
| `Kendall` | any | Base Mann-Kendall S-statistic |
| `modifiedmk` | any | H-R98 variance-corrected MK (`mmkh`) |
| `trend` | any | Sen's slope (`sens.slope`) |
| `changepoint` | any | PELT mean-shift detection (`cpt.mean`) |
| `parallel` | base | Cross-platform multicore processing |
| `openxlsx` | any | Excel workbook creation with styling |
| `writexl` | any | Fast Excel export (Script 1) |
| `gifski` | any | GIF encoding from PNG frames |
| `ragg` | any | Fast AGG PNG renderer for animation frames |
| `scico` | any | Perceptually uniform colour palettes |
| `ggspatial` | any | Cartographic elements in ggplot |
| `copula` | any | Bivariate copula fitting (Script 9) |
| `fitdistrplus` | any | MLE marginal distribution fitting (Script 9) |
| `viridis` | any | Colour-blind-safe palettes |
| `scales` | any | Scale transformations for ggplot |
| `RColorBrewer` | any | Diverging palettes for drought maps |
| `cowplot` | any | Plot alignment utilities |
| `dplyr` | any | Data manipulation (Scripts 10a, 10b, 11) |
| `tidyr` | any | Data reshaping (Scripts 10a, 10b) |
| `ggtext` | any | Rich text in ggplot labels (Script 10a) |
| `cli` | any | Progress messages (Script 5) |
| `gridExtra` | any | Summary table grobs (Script 9) |

---

*This document describes the framework as implemented across the 14 R scripts. All statistical methods follow peer-reviewed literature references cited in the script headers (Allen et al. 1998; Hamed & Rao 1998; Yue et al. 2002; Killick et al. 2012; Sheffield & Wood 2008; Kao & Govindaraju 2010; Yevjevich 1967).*