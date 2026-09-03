# Nechako drought workflow: FAO-56 PM + ERA5-Land native PEV

## Scientific definition

The current workflow uses two observed evaporative-demand products in parallel:

1. `PET_PM`: independently calculated FAO-56 Penman-Monteith reference evapotranspiration from ERA5-Land meteorological forcing.
2. `PET_ERA5`: ERA5-Land native `potential_evaporation` (`pev`), converted to mm/day and subsequently to monthly depth where required.

`PET_ERA5` is **not** relabelled as Penman-Monteith evapotranspiration. It is retained as the ERA5-Land native potential-evaporation product. Consequently, differences between `PET_PM` and `PET_ERA5` are interpreted as **PET-source/product sensitivity**, not as disagreement between two Penman-Monteith implementations.

The current SPEI framework therefore contains two observed branches:

- `SPEI-PM`: precipitation minus FAO-56 PM PET.
- `SPEI-ERA5`: precipitation minus ERA5-Land native PEV.

The former counterfactual/detrended PET branch is not part of the current SPEI analysis. `SPEI_PM - SPEI_ERA5` is a formulation-sensitivity residual and is not interpreted as an additive physical decomposition.

## Analysis period and standard scales

The principal record extends from **1950 through December 2025**. The standard drought-index accumulation scales are:

`1, 3, 6, 12, 24, 36` months.

These scales are the common manuscript, spatial-analysis, and event-ranking scales. The SWEI analysis remains a separate 3-month index where applicable.

The recent drought comparison period used throughout the event-ranking and climatological water-budget analyses is **January 2022–December 2025**.

## Current execution order

1. `Monthly_ERALand_download.R`
2. `2b_PET_ERALand.R`
3. `1SPI_ERALand.R`
4. `3SPEI_ERALand.R`
5. `4pr_pet_trends.r`
6. `4pr_pet_trends_visualization.r`
7. `6trend_test_ALL.R`
8. `7basin_timeseries.R`
9. `8trends_visualization.R`
10. `9event_ranking.R`
11. `10b_PET_bias_nonstationarity.R`
12. `12PET_formulation_sensitivity.R`
13. `13meteorological_controls.R`
14. `15PET_shapley_decomposition.R`
15. `16_climatology_water_budget.R`
16. `5gif_animations_ERA5Land.R` and `11nechako_drought_manuscript.R` are reporting stages.

`16_climatology_water_budget.R` must be run **after** `4pr_pet_trends.r`, because it reads `trend_analysis_pr_pet/analysis_metadata.rds` and uses the basin-averaged monthly climate record stored there.

## Core ERA5-Land data products

The downloader adds `potential_evaporation` to the ERA5-Land monthly archive.

`2b_PET_ERALand.R` produces:

- `monthly_data_direct/potential_evapotranspiration_monthly.nc` — FAO-56 PM PET, mm/day.
- `monthly_data_direct/potential_evaporation_era5land_monthly.nc` — ERA5-Land native PEV, mm/day.
- `monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv` — basin summary of PM PET.
- `monthly_data_direct/ERA5Land_Nechako_PET_ERA5Land_monthly_summary.csv` — basin summary of native ERA5-Land PEV.

`3SPEI_ERALand.R` produces:

- `spei_results_seasonal/` — SPEI-PM.
- `spei_results_seasonal_era5land/` — SPEI-ERA5 using native ERA5-Land PEV.

## ERA5-Land accumulated-field issue

The downloader requests accumulated variables including `potential_evaporation` using the `monthly_averaged_reanalysis_by_hour_of_day` product at `00:00`. This is intentional because of the documented accumulated-field issue affecting the standard monthly-averaged ERA5-Land product for part of the 2022–2024 period.

The workflow performs unit checks before water-budget calculations. ERA5-Land monthly averaged accumulated fields are treated as mean daily rates where appropriate and converted to monthly depth using the actual number of days in each calendar month. A value already expressed as mm/month must not be multiplied by the number of days again.

## Event definition and return-period analysis: `9event_ranking.R`

`9event_ranking.R` is the authoritative drought-event ranking and return-period workflow for SPI, SPEI-PM, SPEI-ERA5, and SWEI.

### SW08-style event identification

Events use a **single-threshold run method** adapted from Sheffield & Wood (2008):

- onset: index `< -0.5`;
- termination: index `>= -0.5`;
- no hysteresis between onset and termination;
- event intensity: `I = mean(-0.5 - x_t)` over event months, so intensity measures deficit below the drought threshold;
- event severity: `S = I × D`, where `D` is event duration in months;
- minimum duration is **scale-specific and equal to the accumulation scale** (`D_min = k`) for the manuscript scales 1, 3, 6, 12, 24, and 36 months.

Thus the minimum-duration gate is 1, 3, 6, 12, 24, and 36 months for the corresponding scales. Runs shorter than the required duration are excluded from the historical qualifying-event catalog. Recent right-censored runs are retained separately for the recent-event visualization when relevant.

Duration classes follow the Sheffield & Wood convention used by the script:

- `D1-2 (Sub-threshold)` — less than 3 months;
- `D3-6 (Short-term)` — 3–6 months;
- `D7-12 (Medium-term)` — 7–12 months;
- `D12+ (Long-term)` — 13 months or longer.

The `-0.5` threshold corresponds to approximately the 31st percentile of a standard normal distribution. It is therefore a deliberately more inclusive onset criterion than the conventional McKee et al. `-1.0` drought threshold. The threshold should not be conflated with the approximately `-1.28` standardized-normal equivalent of a 10th-percentile threshold.

### Return periods

Duration-class return periods are summarized over the 1950–2025 record. The central event-frequency estimate is based on the observed number of events relative to the record length, and 95% confidence limits are calculated using the exact Poisson/chi-squared formulation implemented in the script.

The event-ranking workflow also evaluates **Kendall return periods** using bivariate duration–severity dependence. Copula diagnostics include Kendall's tau, the copula probability at the target event, and the resulting Kendall return period. Multiple Archimedean copula families and candidate marginal distributions are evaluated before selecting the best-supported formulation.

### Event-ranking outputs

Primary files in `temporal_drought/event_ranking/` include:

- `ranked_event_catalog_SW08.csv` — ranked event catalog with dates, duration, intensity, severity, percentile rank, return period, and recent-period flag.
- `return_period_by_duration_class_SW08.csv` — event counts, mean duration, mean intensity, mean severity, maximum severity, return period, and 95% limits by index and duration class.
- `Table_SW08_vs_Hyst_Comparison.csv` — comparison of the current single-threshold SW08 definition with the legacy hysteresis definition.
- `Table_KendallRP_Recent_AllScales_AllEvents.csv` — Kendall return-period results for recent target events across indices and scales.
- `Table_KendallRP_Recent_Primary_AllScales.csv` — primary recent event from each index/scale combination.
- `Table_KendallRP_Recent_SPEI_PM_2022_2025_AllScales.csv` — SPEI-PM practitioner-focused Kendall return-period table.
- `Table_KendallRP_MultiIndex.csv` — multi-index Kendall return-period comparison.
- `Table_1a_SPI_01_02_TopEvents.csv`, `Table_1c_SPEI_Thw_TopEvents.csv`, and `Table_1bc_SPEI3_PM_vs_Thw_Comparison.csv` — legacy/auxiliary ranking tables when those sections are generated; these should not be interpreted as evidence that Thornthwaite PET is part of the current SPEI framework.
- `ms2_annual_severity_series.csv` — annual maximum SPEI-3 event severity used by the downstream MS2 analysis.
- `ms2_copula_fit.rds` — exported copula/marginal fit objects and stationary Kendall-return-period reference when the Kendall analysis succeeds.

Figures include the severity-ranking, return-period, Kendall-RP, copula-diagnostic, multi-index, and recent-event products written by the script, including `Fig3_MS_Drought_Severity_Ranking.pdf`, `Fig3_MS_Drought_Severity_Ranking_KendallRP.pdf`, `Fig4c_Copula_Diagnostics.pdf`, and `Fig4d_KendallRP_MultiIndex_Comparison.png`.

## Climatology and water-budget analysis: `16_climatology_water_budget.R`

`16_climatology_water_budget.R` adds a basin-scale climatological and water-budget layer to the drought analysis. It uses the monthly basin-average fields in `trend_analysis_pr_pet/analysis_metadata.rds` and covers **1950–2025**, with a dedicated **2022–2025** recent-period analysis.

The core variables are:

- precipitation `P`;
- actual evapotranspiration `AET`;
- FAO-56 PM PET;
- optional ERA5-Land native PEV when the authoritative monthly PEV summary is available;
- basin-average air temperature.

### Climatological summaries

The script produces:

- monthly climatology for 1950–2025;
- seasonal climatology for DJF, MAM, JJA, and SON for 1950–2025;
- full-year annual water-budget series for 1950–2025;
- annual climatological summary statistics across the full record.

A critical distinction is that `4pr_pet_trends.r` uses warm-season PM PET totals (April–October) for PET trend analysis, whereas `16_climatology_water_budget.R` deliberately computes **full-year PET totals** for the water-budget analysis. These quantities must not be substituted for one another.

### Water-budget quantities and ratios

The script computes annual and period-level quantities including:

- `P - AET`;
- `P - PET_PM`;
- `AET/P`;
- `PET_PM/P`;
- `AET/PET_PM`;
- `(P - AET)/P` as `residual_runoff_over_P`;
- where ERA5-Land PEV is available, `PET_ERA5` and `PET_ERA5/P`.

Ratios used for annual trend analysis are calculated from **annual aggregated depths**, rather than by averaging monthly ratios. This avoids a compositional weighting problem when the scientific question concerns change in the annual water budget.

`residual_runoff_over_P = (P - AET)/P` is only a **diagnostic residual runoff-ratio proxy**. It is not observed runoff and should not be described as streamflow-derived runoff because storage change, groundwater exchange, snowpack dynamics, and other water-balance terms are omitted.

### Recent drought-period diagnostics

For 2022–2025 the script writes monthly, seasonal, and annual water-budget summaries. It also compares each recent calendar month with its corresponding 1950–2025 calendar-month climatology to provide precipitation, AET, PM-PET, and temperature anomalies.

The resulting anomaly analysis is intended to distinguish whether recent drought conditions were associated with precipitation deficits, anomalously high evaporative demand, anomalous temperature, or reduced actual evaporation relative to climatology.

### Water-budget outputs

Files in `temporal_drought/climatology_water_budget/` are:

- `climatology_monthly_1950_2025.csv`
- `climatology_seasonal_1950_2025.csv`
- `annual_water_budget_1950_2025.csv`
- `climatology_annual_1950_2025.csv`
- `recent_2022_2025_monthly_climate_water_budget.csv`
- `recent_2022_2025_seasonal_climate_water_budget.csv`
- `recent_2022_2025_annual_climate_water_budget.csv`
- `recent_2022_2025_monthly_climatological_anomalies.csv`
- `annual_water_budget_ratio_trends_1950_2025.csv`
- `practitioner_annual_water_budget_1950_2025.csv`

The annual trend table reports Mann-Kendall `tau`, Mann-Kendall `p`, and Sen's slope for precipitation, AET, PM PET, water-balance differences, and the annual ratios. Ratio slopes are reported as ratio/year; depth and temperature slopes retain their native units/year.

## Interpretation framework

The current analysis separates three concepts that should not be conflated:

1. **Drought severity** — derived from standardized SPI/SPEI/SWEI event indices and the SW08 event definition.
2. **PET formulation sensitivity** — the difference between the independently calculated PM reference ET branch and the ERA5-Land native PEV branch.
3. **Physical water-budget context** — precipitation, AET, PET, temperature, and residual water-budget diagnostics from `16_climatology_water_budget.R`.

The water-budget residual `(P-AET)` is not a direct estimate of discharge. Likewise, the PM-versus-ERA5-Land PEV difference is not a thermodynamic decomposition. These quantities are diagnostic/contextual variables and should be interpreted according to their definitions.

## Legacy / removed logic

The following are not part of the current core SPEI framework:

- detrended/counterfactual PET branches;
- interpretation of `SPEI_PM - SPEI_ERA5` as an additive physical decomposition;
- Thornthwaite PET as a current SPEI branch;
- PM-versus-Thornthwaite SPEI as the principal formulation-sensitivity comparison.

The Shapley analysis remains valid for its stated purpose: decomposition of the independently calculated FAO-56 PM PET anomaly. It is not a decomposition of SPEI or of the PM-versus-ERA5-Land PEV difference.

Legacy hysteresis event catalogs are retained only to document sensitivity of event counts and severity metrics to event-definition choice. The current event-ranking results use the SW08-style single-threshold definition described above.

## Reproducibility notes

All basin spatial aggregation in the drought-index workflow uses the project equal-area CRS, `EPSG:3005`, and the project-specific basin boundary. Basin-averaged drought-index CSVs generated upstream are treated as authoritative by downstream time-series and ranking scripts rather than being recomputed independently.

For the climatological water-budget analysis, `16_climatology_water_budget.R` reads the authoritative basin-average monthly climate record from `analysis_metadata.rds`. ERA5-Land native PEV is included only when `ERA5Land_Nechako_PET_ERA5Land_monthly_summary.csv` is present and contains the required `year`, `month`, and `mean_pet` fields.
