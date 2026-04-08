# ==============================================================================
#   w11_dynamic_thermodynamic_decomp.R (MERGED VERSION)
# DYNAMIC vs. THERMODYNAMIC DROUGHT DECOMPOSITION
# Nechako River Basin, BC — 2022-2025 Drought Study
# Loads pre-computed SPEI from 3SPEI_ERALand.R outputs
# - SPEI_PM  (spei_results_seasonal/)      = full observed drought
# - SPEI_Thw (spei_results_seasonal_thw/)  = thermodynamic component
# - Difference ≈ dynamic component
# Needs the SPEI output directories from 3SPEI_ERALand.R
# ==============================================================================
library(terra)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "decomp_results"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ==============================================================================
#   STEP 1: LOAD PRE-COMPUTED SPEI RESULTS
# ==============================================================================
cat("\n===== STEP 1: LOAD PRE-COMPUTED SPEI RESULTS =====\n")
# Directories from 3SPEI_ERALand.R
spei_pm_dir   <- "spei_results_seasonal"
spei_thw_dir  <- "spei_results_seasonal_thw"
if (!dir.exists(spei_pm_dir))  stop("SPEI_PM directory not found. Run 3SPEI_ERALand.R first.")
if (!dir.exists(spei_thw_dir)) stop("SPEI_Thw directory not found. Run 3SPEI_ERALand.R first.")
# Scales available from 3SPEI_ERALand.R
decomp_scales <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24)
# Helper function to load SPEI NetCDF files for one scale
load_spei_scale <- function(scale, pet_type, target_dir) {
  cat(sprintf("  Loading SPEI_%s-%d...\n", pet_type, scale))
  spei_list  <- list()
  month_names  <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  for (m in 1:12) {
    nc_file <- file.path(target_dir,
                         sprintf("spei_%02d_month%02d_%s.nc", scale, m, month_names[m]))
    if (file.exists(nc_file)) {
      r <- rast(nc_file)
      spei_list[[m]] <- values(r)
    } else {
      # Fallback to CSV if NetCDF missing
      csv_file <- file.path(target_dir,
                            sprintf("spei_%02d_month%02d_%s.csv", scale, m, month_names[m]))
      if (file.exists(csv_file)) {
        df <- read.csv(csv_file)
        spei_list[[m]] <- as.matrix(df[, -c(1, 2)])  # remove lon/lat columns
      } else {
        warning(sprintf("No SPEI-%d month %d file found", scale, m))
        spei_list[[m]] <- NULL
      }
    }
  }
  # Combine months into matrix [pixels × time]
  valid_months <- which(!sapply(spei_list, is.null))
  if (length(valid_months) == 0) return(NULL)
  spei_mat <- do.call(cbind, spei_list[valid_months])
  return(spei_mat)
}
# Load all scales for both PET types in one pass (avoids two separate loops)
scale_data  <- lapply(setNames(decomp_scales, paste0("scale_", decomp_scales)), function(sc) {
  list(
    pm  = load_spei_scale(sc, "PM",  spei_pm_dir),
    thw = load_spei_scale(sc, "Thw", spei_thw_dir)
  )
})
spei_pm_all   <- lapply(scale_data, `[[`, "pm")
spei_thw_all  <- lapply(scale_data, `[[`, "thw")
cat("✓ SPEI results loaded from both directories\n")

# ==============================================================================
#   STEP 2: LOAD METADATA (dates, basin mask)
# ==============================================================================
cat("\n===== STEP 2: LOAD METADATA =====\n")
# Load dates from precipitation file (same as 3SPEI_ERALand.R)
precip <- rast("monthly_data_direct/total_precipitation_monthly.nc")
dates  <- as.Date(time(precip))
if (is.null(dates) || all(is.na(dates))) {
  cat("Reconstructing dates from 1950-01...\n")
  dates <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(precip))
}
month_nums <- as.integer(format(dates, "%m"))
year_nums  <- as.integer(format(dates, "%Y"))
# Load basin boundary for masking
basin <- vect("Spatial/nechakoBound_dissolve.kmz")
if (nrow(basin) > 1L) basin <- aggregate(basin)
target_crs <- "EPSG:3005"
if (!same.crs(basin, target_crs)) basin <- project(basin, target_crs)
# Get basin mask from first SPEI file (non-NA pixels)
spei_sample <- rast(file.path(spei_pm_dir, "spei_06_month01_Jan.nc"))
basin_mask <- !is.na(values(spei_sample))
n_pixels <- sum(basin_mask)
cat(sprintf("✓ Metadata loaded: %d months, %d basin pixels\n", length(dates), n_pixels))

# ==============================================================================
#   STEP 3: COMPUTE DECOMPOSITION FROM LOADED SPEI
# ==============================================================================
cat("\n===== STEP 3: COMPUTE DYNAMIC/THERMODYNAMIC DECOMPOSITION =====\n")
# Decomposition framework (from 3SPEI_ERALand.R documentation):
#   SPEI_PM         = full observed drought severity
# SPEI_Thw        = thermodynamic (temperature-driven) component
# SPEI_PM - SPEI_Thw ≈ dynamic component (precip + radiation + wind)
compute_decomposition <- function(spei_pm, spei_thw, dates, month_nums, year_nums, scale) {
  n_time <- ncol(spei_pm)
  # Basin-mean timeseries
  spei_pm_mean   <- colMeans(spei_pm,   na.rm = TRUE)
  spei_thw_mean  <- colMeans(spei_thw,  na.rm = TRUE)
  spei_dyn_mean  <- spei_pm_mean - spei_thw_mean  # dynamic component
  data.frame(
    date          = dates,
    year          = year_nums,
    month         = month_nums,
    scale         = scale,
    SPEI_PM       = spei_pm_mean,
    SPEI_Thw      = spei_thw_mean,
    SPEI_Dynamic  = spei_dyn_mean,
    # Fraction of total drought severity from each component
    thm_frac      = abs(spei_thw_mean) / (abs(spei_pm_mean) + 1e-6),
    dyn_frac      = abs(spei_dyn_mean) / (abs(spei_pm_mean) + 1e-6)
  )
}
decomp_all  <- bind_rows(lapply(decomp_scales, function(sc) {
  pm_mat   <- spei_pm_all[[paste0("scale_", sc)]]
  thw_mat  <- spei_thw_all[[paste0("scale_", sc)]]
  if (is.null(pm_mat) || is.null(thw_mat)) {
    cat(sprintf("  Skipping scale %d (missing data)\n", sc)); return(NULL)
  }
  min_time  <- min(ncol(pm_mat), ncol(thw_mat))
  result  <- compute_decomposition(
    pm_mat[, 1:min_time, drop=FALSE],
    thw_mat[, 1:min_time, drop=FALSE],
    dates[1:min_time], month_nums[1:min_time], year_nums[1:min_time], sc
  )
  cat(sprintf("  ✓ Scale %d decomposition complete\n", sc))
  result
}))
write.csv(decomp_all, file.path(out_dir, "decomp_full_timeseries.csv"), row.names = FALSE)
cat("✓ Decomposition timeseries saved\n")
# 
# ==============================================================================
#   STEP 4: 2022-2025 DROUGHT DECOMPOSITION SUMMARY
# ==============================================================================
cat("\n===== STEP 4: 2022-2025 DROUGHT DECOMPOSITION =====\n")
drought_decomp <- decomp_all %>%
  filter(year >= 2022, year <= 2025) %>%
  group_by(scale) %>%
  summarise(
    SPEI_PM_mean        = mean(SPEI_PM,      na.rm = TRUE),
    SPEI_Thw_mean       = mean(SPEI_Thw,     na.rm = TRUE),
    SPEI_Dynamic_mean   = mean(SPEI_Dynamic, na.rm = TRUE),
    thm_frac_mean       = mean(thm_frac,     na.rm = TRUE),
    dyn_frac_mean       = mean(dyn_frac,     na.rm = TRUE),
    # Percentage contribution to drought severity
    thm_pct = 100 * mean(abs(SPEI_Thw),    na.rm = TRUE) /
      (mean(abs(SPEI_PM),    na.rm = TRUE) + 1e-6),
    dyn_pct = 100 * mean(abs(SPEI_Dynamic), na.rm = TRUE) /
      (mean(abs(SPEI_PM),    na.rm = TRUE) + 1e-6),
    .groups = "drop"
  )
cat("\n2022-2025 Drought Decomposition Summary:\n")
print(drought_decomp)
write.csv(drought_decomp, file.path(out_dir, "decomp_2022_2025_summary.csv"), row.names = FALSE)
# 
# ==============================================================================
#   STEP 5: LONG-TERM TREND IN THERMODYNAMIC FRACTION
# ==============================================================================
#   STEP 5: TREND IN THERMODYNAMIC FRACTION (1950-2025)
# Two OLS models:
#   (a) Full record 1950-2025  — slope, 95% CI, p-value
# (b) Post-1990 sub-period   — slope and p-value for annotation
# Both are computed here and used in manuscript Figure 6.
# ==============================================================================
cat("\n===== STEP 5: TREND IN THERMODYNAMIC FRACTION (1950-2025) =====\n")
# Focus on SPEI-3 JJA (summer drought most relevant for Nechako)
trend_df <- decomp_all %>%
  dplyr::filter(scale == 3, month %in% c(6, 7, 8)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    SPEI_PM       = mean(SPEI_PM,      na.rm = TRUE),
    SPEI_Thw      = mean(SPEI_Thw,     na.rm = TRUE),
    SPEI_Dynamic  = mean(SPEI_Dynamic, na.rm = TRUE),
    thm_frac_abs  = mean(abs(SPEI_Thw) /
                           (abs(SPEI_PM) + 1e-6), na.rm = TRUE),
    .groups = "drop"
  )
# ── Full-record OLS (1950-2025) ────────────────────────────────────────────
thm_trend      <- lm(thm_frac_abs ~ year, data = trend_df)
thm_slope      <- coef(thm_trend)["year"]
thm_slope_dec  <- thm_slope * 10          # per decade
thm_pval       <- summary(thm_trend)$coefficients["year", 4]
thm_ci95       <- confint(thm_trend, "year", level = 0.95)
thm_ci90       <- confint(thm_trend, "year", level = 0.90)
cat(sprintf("Full record (1950-2025):\n"))
cat(sprintf("  Slope    : %+.5f yr⁻¹  (%+.4f per decade)\n",
            thm_slope, thm_slope_dec))
cat(sprintf("  95%% CI   : [%+.5f, %+.5f] yr⁻¹\n",
            thm_ci95[1], thm_ci95[2]))
cat(sprintf("  p-value  : %.4f\n", thm_pval))
# ── Post-1990 sub-period OLS ───────────────────────────────────────────────
trend_post90   <- dplyr::filter(trend_df, year >= 1990)
thm_trend_p90  <- lm(thm_frac_abs ~ year, data = trend_post90)
thm_slope_p90  <- coef(thm_trend_p90)["year"]
thm_pval_p90   <- summary(thm_trend_p90)$coefficients["year", 4]
cat(sprintf("\nPost-1990 sub-period:\n"))
cat(sprintf("  Slope    : %+.5f yr⁻¹  (%+.4f per decade)\n",
            thm_slope_p90, thm_slope_p90 * 10))
cat(sprintf("  p-value  : %.4f\n", thm_pval_p90))
# Helper: format p-value for annotation
fmt_p  <- function(p) {
  if      (p  < 0.001)  "p  < 0.001"
  else if (p  < 0.01)  sprintf("p = %.3f", p)
  else if (p  < 0.05)  sprintf("p = %.3f", p)
  else if (p  < 0.10)  sprintf("p = %.3f", p)
  else                sprintf("p = %.2f",  p)
}
# Annotation strings used in Figure 6
ann_full <- sprintf("%+.4f yr⁻¹, %s\n(1950–2025)",
                    thm_slope, fmt_p(thm_pval))
ann_p90  <- sprintf("%+.4f yr⁻¹, %s\n(1990–2025)",
                    thm_slope_p90, fmt_p(thm_pval_p90))

# ==============================================================================
#   STEP 6: GENERATE FIGURES
# SUPPLEMENTARY figures (original, kept for reference):
#   fig_decomp_spei3_timeseries.pdf  — SPEI-3 full record decomposition
# fig_annual_decomp_barplot.pdf    — annual SPEI-3 drought years bar
# fig_thermodynamic_trend.pdf      — JJA SPEI-3 thermodynamic trend
# fig_scale_comparison.pdf         — all scales 1-24 percentage stack
# MANUSCRIPT figures (new, publication quality):
#   Fig5a_decomp_timeseries_2020_2025.pdf/.png
# SPEIₚₘ, SPEI₀ (=SPEIᵟʰʷ), and SPEIᵟ on a single axis, scales 1–3,
# 2020–2025 focus window. Core figure for Section 4.3.
# Fig5b_scale_comparison_stacked.pdf/.png
# Stacked percentage bar chart, scales 1–3 only, showing thermodynamic
# and dynamic percentage contributions with data labels. Replaces/
#   accompanies Table 3. Core result of the decomposition framework.
# ==============================================================================
cat("\n===== STEP 6: GENERATE FIGURES =====\n")
# ── Shared style ──────────────────────────────────────────────────────────────
theme_ms <- ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0),
    plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey40", hjust = 0,
                                          margin = ggplot2::margin(b = 4)),
    legend.position  = "bottom",
    legend.key.size  = ggplot2::unit(0.4, "cm"),
    legend.text      = ggplot2::element_text(size = 8.5),
    legend.title     = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "grey94", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin      = ggplot2::margin(4, 8, 4, 4)
  )
# Component colour palette (consistent across all manuscript figures)
COL_PM  <- "#1a1a1a"    # near-black  — SPEIₚₘ (full observed)
COL_THW <- "#e31a1c"    # red         — SPEI₀  (thermodynamic)
COL_DYN <- "#1f78b4"    # blue        — SPEIᵟ  (dynamic)
HIGHLIGHT_START <- as.Date("2022-01-01")
HIGHLIGHT_END   <- as.Date("2025-12-31")
THR_EVENT <- -0.5   # event-identification threshold (Section 3.2 of manuscript)

# ==============================================================================
#   SUPPLEMENTARY FIGURES (original code, lightly tidied)
# ==============================================================================
#   S1: SPEI-3 full-record decomposition (1950-2025)
df_s1  <- decomp_all %>% dplyr::filter(scale == 3)
pdf(file.path(out_dir,  "fig_decomp_spei3_timeseries.pdf"), width = 14, height = 8)
df_s1 %>%
  dplyr::select(date, SPEI_PM, SPEI_Thw, SPEI_Dynamic) %>%
  tidyr::pivot_longer(cols = -date, names_to = "component", values_to = "SPEI") %>%
  dplyr::mutate(component = dplyr::recode(component,
                                          SPEI_PM      = "Full (PM PET)",
                                          SPEI_Thw     = "Thermodynamic (Thw PET)",
                                          SPEI_Dynamic = "Dynamic (PM − Thw)")) %>%
  ggplot2::ggplot(ggplot2::aes(x = date, y = SPEI, colour = component)) +
  ggplot2::geom_hline(yintercept = 0, colour = "gray50") +
  ggplot2::geom_hline(yintercept = THR_EVENT, linetype = "dashed",
                      colour = "grey40", alpha = 0.7) +
  ggplot2::geom_line(alpha = 0.8) +
  ggplot2::annotate("rect",
                    xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                    ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "orange") +
  ggplot2::scale_colour_manual(values = c(COL_PM, COL_THW, COL_DYN)) +
  ggplot2::labs(title = "Nechako Basin: SPEI-3 Dynamic vs. Thermodynamic Decomposition (1950–2025)",
                subtitle = "Orange = 2022–2025. Dashed = event threshold (x₀ = −0.5)",
                x = NULL, y = "SPEI-3") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")
dev.off()
cat("  ✓ Supplementary S1: fig_decomp_spei3_timeseries.pdf\n")
# S2: Annual SPEI-3 drought-year stacked bars
drought_years_df  <- decomp_all %>%
  dplyr::filter(scale == 3, SPEI_PM < -0.8) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    SPEI_PM       = mean(SPEI_PM,       na.rm = TRUE),
    Thermodynamic = mean(SPEI_Thw,      na.rm = TRUE),
    Dynamic       = mean(SPEI_Dynamic,  na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::filter(abs(SPEI_PM) > 0.3) %>%
  tidyr::pivot_longer(cols = c(Dynamic, Thermodynamic),
                      names_to = "component", values_to = "contribution")
pdf(file.path(out_dir,  "fig_annual_decomp_barplot.pdf"), width = 12, height = 6)
ggplot2::ggplot(drought_years_df,
                ggplot2::aes(x = factor(year), y = contribution, fill = component)) +
  ggplot2::geom_bar(stat = "identity", position = "stack") +
  ggplot2::geom_hline(yintercept = 0, colour = "gray40") +
  ggplot2::scale_fill_manual(
    values = c(Dynamic = COL_DYN, Thermodynamic = COL_THW)) +
  ggplot2::labs(title = "Annual Drought Severity Decomposition: SPEI-3 Drought Years",
                x = "Year", y = "SPEI contribution", fill = "Component") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
dev.off()
cat("  ✓ Supplementary S2: fig_annual_decomp_barplot.pdf\n")
# S3: Seasonal SPEI-3 thermodynamic trends (4 versions: MAM, JJA, SON, DJF)
cat("\n  Building S3: 4 seasonal versions of thermodynamic trend...\n")
seasons <- list(
  MAM = c(3, 4, 5),   # Spring
  JJA = c(6, 7, 8),   # Summer
  SON = c(9, 10, 11), # Autumn
  DJF = c(12, 1, 2)   # Winter
)

for (s_name in names(seasons)) {
  s_months <- seasons[[s_name]]
  
  # Calculate subset for the specific season
  s_trend_df <- decomp_all %>%
    dplyr::filter(scale == 3, month %in% s_months) %>%
    # Shift December forward by 1 year so DJF forms a continuous winter
    dplyr::mutate(met_year = dplyr::if_else(month == 12, year + 1, year)) %>%
    dplyr::group_by(met_year) %>%
    dplyr::summarise(
      SPEI_PM       = mean(SPEI_PM,      na.rm = TRUE),
      SPEI_Thw      = mean(SPEI_Thw,     na.rm = TRUE),
      thm_frac_abs  = mean(abs(SPEI_Thw) / (abs(SPEI_PM) + 1e-6), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(year = met_year) %>%
    # Bound the years and drop missing data
    dplyr::filter(!is.na(thm_frac_abs), is.finite(thm_frac_abs), year >= 1950, year <= 2025)
  
  # Calculate specific OLS metrics for the subtitle
  s_thm_trend     <- lm(thm_frac_abs ~ year, data = s_trend_df)
  s_thm_slope_dec <- coef(s_thm_trend)["year"] * 10
  s_thm_pval      <- summary(s_thm_trend)$coefficients["year", 4]
  
  # Generate and save the plot
  pdf_name <- sprintf("fig_thermodynamic_trend_%s.pdf", s_name)
  pdf(file.path(out_dir, pdf_name), width = 10, height = 5)
  print(
    ggplot2::ggplot(s_trend_df, ggplot2::aes(x = year, y = thm_frac_abs)) +
      ggplot2::geom_point(alpha = 0.5, colour = COL_THW) +
      ggplot2::geom_smooth(method = "lm", colour = COL_THW,
                           fill = COL_THW, alpha = 0.2) +
      ggplot2::geom_vline(xintercept = 2022, linetype = "dashed", colour = "orange") +
      ggplot2::labs(
        title    = sprintf("Growing thermodynamic fraction — %s SPEI-3, Nechako Basin", s_name),
        subtitle = sprintf("Linear trend: %+.3f per decade (%s)",
                           s_thm_slope_dec, fmt_p(s_thm_pval)),
        x = "Year", y = "Thermodynamic fraction of drought severity") +
      ggplot2::theme_minimal(base_size = 12)
  )
  dev.off()
  cat(sprintf("  ✓ Supplementary S3: %s\n", pdf_name))
}
# ==============================================================================
#   MANUSCRIPT FIGURE 5a — Decomposition time series, scales 1-3, 2020-2025
# Section 4.3: shows SPEIₚₘ, SPEI₀, and SPEIᵟ on a single axis per scale
# ==============================================================================
cat("\n  Building Fig5a: decomposition time series (2020–2025, scales 1–3)...\n")
ms_scales_ts <- c(1, 2, 3)
panels_ts <- lapply(ms_scales_ts, function(sc) {
  df_sc  <- decomp_all %>%
    dplyr::filter(scale == sc,
                  date  >= as.Date("2020-01-01"),
                  date  <= as.Date("2025-12-31")) %>%
    dplyr::select(date, SPEI_PM, SPEI_Thw, SPEI_Dynamic) %>%
    tidyr::pivot_longer(cols = -date,
                        names_to  = "component",
                        values_to = "value") %>%
    dplyr::mutate(
      component = dplyr::recode(component,
                                SPEI_PM      = "SPEIₚₘ (full, PM PET)",
                                SPEI_Thw     = "SPEI₀ (thermodynamic)",
                                SPEI_Dynamic = "SPEIᵟ (dynamic)"),
      component = factor(component,
                         levels = c("SPEIₚₘ (full, PM PET)",
                                    "SPEI₀ (thermodynamic)",
                                    "SPEIᵟ (dynamic)"))
    )
  # y-range
  ylo <- min(-2.2, min(df_sc$value, na.rm = TRUE) * 1.1)
  yhi <- max( 1.5, max(df_sc$value, na.rm = TRUE) * 1.1)
  is_bottom <- sc == max(ms_scales_ts)
  p <- ggplot2::ggplot(df_sc,
                       ggplot2::aes(x = date, y = value,
                                    colour = component, linetype = component)) +
    # 2022-2025 highlight band
    ggplot2::annotate("rect",
                      xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                      ymin = -Inf, ymax = Inf,
                      fill = "grey80", alpha = 0.5) +
    ggplot2::annotate("text",
                      x = HIGHLIGHT_START + 180, y = yhi * 0.93,
                      label = "2022–2025", size = 2.5,
                      hjust = 0, colour = "grey35", fontface = "italic") +
    # Reference lines
    ggplot2::geom_hline(yintercept = 0,
                        colour = "grey50", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = THR_EVENT,
                        colour = "grey35", linewidth = 0.5,
                        linetype = "dashed") +
    # Component lines
    ggplot2::geom_line(linewidth = 0.65, alpha = 0.9) +
    ggplot2::scale_colour_manual(
      values = c(
        "SPEIₚₘ (full, PM PET)"  = COL_PM,
        "SPEI₀ (thermodynamic)"       = COL_THW,
        "SPEIᵟ (dynamic)"             = COL_DYN)) +
    ggplot2::scale_linetype_manual(
      values = c(
        "SPEIₚₘ (full, PM PET)"  = "solid",
        "SPEI₀ (thermodynamic)"       = "solid",
        "SPEIᵟ (dynamic)"             = "dashed")) +
    # Threshold label (first panel only)
    {if (sc == min(ms_scales_ts))
      ggplot2::annotate("text",
                        x = as.Date("2020-03-01"), y = THR_EVENT + 0.07,
                        label = "x₀ = −0.5", size = 2.2,
                        hjust = 0, colour = "grey35")
      else NULL} +
    ggplot2::scale_x_date(
      date_breaks  = "1 year",
      date_labels  = "%Y",
      expand       = ggplot2::expansion(add = c(15, 15))) +
    ggplot2::scale_y_continuous(
      breaks = c(-2, -1.5, -1, -0.5, 0, 0.5, 1),
      limits = c(ylo, yhi),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title    = sprintf("(%s)  SPEI-%d",
                         letters[which(ms_scales_ts == sc)], sc),
      subtitle = sprintf("Thermodynamic = %.0f%%; Dynamic = %.0f%% of 2022–2025 severity",
                         drought_decomp$thm_pct[drought_decomp$scale == sc],
                         drought_decomp$dyn_pct[drought_decomp$scale == sc]),
      x        = if (is_bottom) "Year" else NULL,
      y        = sprintf("Index value")) +
    theme_ms +
    ggplot2::theme(
      axis.text.x  = if (is_bottom)
        ggplot2::element_text(size = 7.5)
      else
        ggplot2::element_blank(),
      axis.ticks.x = if (is_bottom)
        ggplot2::element_line()
      else
        ggplot2::element_blank()
    )
  p
})
fig5a  <- patchwork::wrap_plots(panels_ts, ncol = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Dynamic–thermodynamic SPEI decomposition — Nechako River Basin, 2020–2025",
    subtitle = paste0(
      "Black = SPEIₚₘ (PM PET, full observed);   ",
      "Red = SPEI₀ (Thornthwaite PET, thermodynamic);   ",
      "Blue dashed = SPEIᵟ = SPEIₚₘ − SPEI₀ (dynamic)\n",
      "Grey band = 2022–2025 focus period.  Dashed reference = event threshold (x₀ = −0.5)."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey35",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 4)),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 8)
    )
  )
fig5a_pdf <- file.path(out_dir, "Fig5a_decomp_timeseries_2020_2025.pdf")
fig5a_png <- file.path(out_dir, "Fig5a_decomp_timeseries_2020_2025.png")
tryCatch({
  ggplot2::ggsave(fig5a_pdf, fig5a,
                  width = 7.0, height = 8.5, units = "in", device = "pdf")
  cat(sprintf("  ✓ Manuscript Fig5a (PDF): %s\n", basename(fig5a_pdf)))
}, error = function(e) cat(sprintf("  ⚠ Fig5a PDF failed: %s\n", e$message)))
tryCatch({
  ggplot2::ggsave(fig5a_png, fig5a,
                  width = 7.0, height = 8.5, units = "in",
                  dpi = 300, device = "png")
  cat(sprintf("  ✓ Manuscript Fig5a (PNG): %s\n", basename(fig5a_png)))
}, error = function(e) cat(sprintf("  ⚠ Fig5a PNG failed: %s\n", e$message)))

# ==============================================================================
#   MANUSCRIPT FIGURE 5b — Stacked percentage bar chart, scales 1-3 only
# Core visual for Table 3 / Section 4.2
# Shows the shift from ~100% thermodynamic (scale-1) to ~65%/35% (scale-3)
# with exact percentage data labels on each segment.
# ==============================================================================
cat("\n  Building Fig5b: stacked percentage bars (scales 1–3)...\n")
ms_scales_bar <- c(1, 2, 3)
decomp_ms_bar <- drought_decomp %>%
  dplyr::filter(scale %in% ms_scales_bar) %>%
  dplyr::mutate(
    scale_lab = sprintf("SPEI-%d", scale),
    scale_lab = factor(scale_lab,
                       levels = sprintf("SPEI-%d", ms_scales_bar))
  )
# Long format for stacked bars
bar_long  <- decomp_ms_bar %>%
  dplyr::select(scale, scale_lab, thm_pct, dyn_pct) %>%
  tidyr::pivot_longer(cols = c(thm_pct, dyn_pct),
                      names_to  = "component",
                      values_to = "pct") %>%
  dplyr::mutate(
    component = dplyr::recode(component,
                              thm_pct = "Thermodynamic (SPEI₀)",
                              dyn_pct = "Dynamic (SPEIᵟ)"),
    component = factor(component,
                       levels = c("Thermodynamic (SPEI₀)",
                                  "Dynamic (SPEIᵟ)"))
  )
# Label positions: cumulative from bottom within each bar
bar_labels <- bar_long %>%
  dplyr::arrange(scale_lab, dplyr::desc(component)) %>%
  dplyr::group_by(scale_lab) %>%
  dplyr::mutate(
    pct_cum    = cumsum(pct),
    label_y    = pct_cum - pct / 2,
    label_text = sprintf("%.0f%%", pct)
  ) %>%
  dplyr::ungroup()
fig5b <- ggplot2::ggplot(bar_long,
                         ggplot2::aes(x = scale_lab, y = pct, fill = component)) +
  ggplot2::geom_bar(stat = "identity", position = "stack",
                    width = 0.55, colour = "white", linewidth = 0.3) +
  # Percentage data labels inside each segment
  ggplot2::geom_text(
    data     = bar_labels,
    ggplot2::aes(x = scale_lab, y = label_y, label = label_text),
    size     = 3.5,
    fontface = "bold",
    colour   = "white",
    inherit.aes = FALSE) +
  # 100% reference line
  ggplot2::geom_hline(yintercept = 100,
                      linetype = "dashed",
                      colour = "grey50",
                      linewidth = 0.5) +
  ggplot2::scale_fill_manual(
    values = c("Thermodynamic (SPEI₀)" = COL_THW,
               "Dynamic (SPEIᵟ)"       = COL_DYN)) +
  ggplot2::scale_y_continuous(
    limits = c(0, 115),
    breaks = seq(0, 100, by = 20),
    labels = function(x) sprintf("%d%%", x),
    expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::labs(
    title    = "Thermodynamic vs. dynamic contribution to 2022–2025 drought severity",
    subtitle = paste0(
      "Percentage of mean |SPEIₚₘ| attributable to each component.\n",
      "Scales 1–3 months shown (matching manuscript analysis period in Section 3.2.1)."),
    x        = "Accumulation scale",
    y        = "Contribution to drought severity (%)",
    fill     = NULL,
    caption  = paste0(
      "Thermodynamic = SPEI₀ (Thornthwaite PET, temperature-driven);   ",
      "Dynamic = SPEIᵟ = SPEIₚₘ − SPEI₀ (circulation-driven).\n",
      "Percentages may exceed 100% when the two components partially oppose each other.")) +
  theme_ms +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_text(size = 10, face = "bold"),
    axis.text.y  = ggplot2::element_text(size = 9),
    legend.text  = ggplot2::element_text(size = 9),
    plot.caption = ggplot2::element_text(size = 6.5, colour = "grey50",
                                         hjust = 0,
                                         margin = ggplot2::margin(t = 4))
  )
fig5b_pdf <- file.path(out_dir, "Fig5b_scale_comparison_stacked.pdf")
fig5b_png <- file.path(out_dir, "Fig5b_scale_comparison_stacked.png")
tryCatch({
  ggplot2::ggsave(fig5b_pdf, fig5b,
                  width = 5.5, height = 5.0, units = "in", device = "pdf")
  cat(sprintf("  ✓ Manuscript Fig5b (PDF): %s\n", basename(fig5b_pdf)))
}, error = function(e) cat(sprintf("  ⚠ Fig5b PDF failed: %s\n", e$message)))
tryCatch({
  ggplot2::ggsave(fig5b_png, fig5b,
                  width = 5.5, height = 5.0, units = "in",
                  dpi = 300, device = "png")
  cat(sprintf("  ✓ Manuscript Fig5b (PNG): %s\n", basename(fig5b_png)))
}, error = function(e) cat(sprintf("  ⚠ Fig5b PNG failed: %s\n", e$message)))

# ==============================================================================
#   MANUSCRIPT FIGURE 6 — Thermodynamic fraction trend, Section 4.4
# Design:
#   • Scatter: each annual JJA Fthm value (filled circle, coloured by decade)
# • Full-record OLS line + 95% confidence ribbon (solid red)
# • Post-1990 OLS line + 90% confidence ribbon (dashed, lighter red)
# • Vertical reference line at 1990 with text annotation
# • Vertical dashed line at 2022 (start of focus drought)
# • Slope + p-value text annotations for both fits
# • 2022-2025 grey highlight band
# ==============================================================================
cat("\n  Building Fig6: thermodynamic fraction trend (manuscript)...\n")
# Compute 95% CI prediction band for full-record OLS
pred_full <- predict(thm_trend,
                     newdata = data.frame(year = trend_df$year),
                     interval = "confidence", level = 0.95)
trend_df$fit_full  <- pred_full[, "fit"]
trend_df$lwr95     <- pred_full[, "lwr"]
trend_df$upr95     <- pred_full[, "upr"]
# Compute 90% CI prediction band for post-1990 OLS
trend_post90$fit_p90  <- predict(thm_trend_p90,
                                 newdata  = data.frame(year = trend_post90$year),
                                 interval = "confidence", level = 0.90)[, "fit"]
trend_post90$lwr90    <- predict(thm_trend_p90,
                                 newdata  = data.frame(year = trend_post90$year),
                                 interval = "confidence", level = 0.90)[, "lwr"]
trend_post90$upr90    <- predict(thm_trend_p90,
                                 newdata  = data.frame(year = trend_post90$year),
                                 interval = "confidence", level = 0.90)[, "upr"]
# Decade factor for point colouring (gives implicit temporal stratification)
trend_df$decade <- factor(
  paste0(floor(trend_df$year / 10) * 10, "s"),
  levels = paste0(seq(1950, 2020, by = 10), "s"))
# y-axis range
ylo_f6 <- max(0,    min(trend_df$thm_frac_abs, na.rm = TRUE) - 0.05)
yhi_f6 <- min(1.05, max(trend_df$thm_frac_abs, na.rm = TRUE) + 0.12)
# Annotation x-positions
x_ann_full  <- 1952   # slope annotation for full OLS
x_ann_p90   <- 1991   # slope annotation for post-1990 OLS
y_ann_full  <- yhi_f6 - 0.04
y_ann_p90   <- yhi_f6 - 0.04
# ─────────────────────────────────────────────────────────────────────────────
# ORIGINAL FIGURE 6 (with ribbons, OLS lines, p-values)
# ─────────────────────────────────────────────────────────────────────────────
fig6 <- ggplot2::ggplot() +
  # ── 2022-2025 highlight band ───────────────────────────────────────────
  ggplot2::annotate("rect",
                    xmin = 2022, xmax = 2026,
                    ymin = -Inf, ymax = Inf,
                    fill = "grey82", alpha = 0.55) +
  # ── 95% CI ribbon — full record OLS ────────────────────────────────────
  ggplot2::geom_ribbon(
    data = trend_df,
    ggplot2::aes(x = year, ymin = lwr95, ymax = upr95),
    fill  = COL_THW,
    alpha = 0.15,
    inherit.aes = FALSE) +
  # ── 90% CI ribbon — post-1990 OLS ──────────────────────────────────────
  ggplot2::geom_ribbon(
    data = trend_post90,
    ggplot2::aes(x = year, ymin = lwr90, ymax = upr90),
    fill  = COL_THW,
    alpha = 0.10,
    inherit.aes = FALSE) +
  # ── Full-record OLS line ────────────────────────────────────────────────
  ggplot2::geom_line(
    data = trend_df,
    ggplot2::aes(x = year, y = fit_full),
    colour    = COL_THW,
    linewidth = 0.9,
    linetype  = "solid",
    inherit.aes = FALSE) +
  # ── Post-1990 OLS line ──────────────────────────────────────────────────
  ggplot2::geom_line(
    data = trend_post90,
    ggplot2::aes(x = year, y = fit_p90),
    colour    = COL_THW,
    linewidth = 0.75,
    linetype  = "dashed",
    alpha     = 0.8,
    inherit.aes = FALSE) +
  # ── Annual data points coloured by decade ──────────────────────────────
  ggplot2::geom_point(
    data = trend_df,
    ggplot2::aes(x = year, y = thm_frac_abs, fill = decade),
    shape  = 21,
    size   = 2.2,
    colour = "white",
    stroke = 0.4,
    inherit.aes = FALSE) +
  ggplot2::scale_fill_brewer(palette = "RdYlBu", direction = -1,
                             name = "Decade") +
  # ── Reference line: 1990 acceleration ─────────────────────────────────
  ggplot2::geom_vline(
    xintercept = 1990,
    linetype   = "dotted",
    colour     = "grey30",
    linewidth  = 0.6) +
  ggplot2::annotate("text",
                    x      = 1991,
                    y      = ylo_f6 + 0.03,
                    label  = "Post-1990\nacceleration",
                    hjust  = 0, vjust = 0,
                    size   = 2.6,
                    colour = "grey30",
                    fontface = "italic") +
  # ── Reference line: 2022 drought onset ────────────────────────────────
  ggplot2::geom_vline(
    xintercept = 2022,
    linetype   = "dashed",
    colour     = "grey50",
    linewidth  = 0.5) +
  ggplot2::annotate("text",
                    x      = 2021.5,
                    y      = ylo_f6 + 0.03,
                    label  = "2022–2025\ndrought",
                    hjust  = 1, vjust = 0,
                    size   = 2.5,
                    colour = "grey40",
                    fontface = "italic") +
  # ── Slope annotation: full-record OLS ─────────────────────────────────
  ggplot2::annotate("label",
                    x          = x_ann_full,
                    y          = y_ann_full,
                    label      = paste0("Full record OLS\n", ann_full),
                    hjust      = 0, vjust = 1,
                    size       = 2.7,
                    fontface   = "plain",
                    colour     = COL_THW,
                    fill       = "white",
                    label.size = 0.2,
                    label.padding = ggplot2::unit(0.15, "lines"),
                    label.r    = ggplot2::unit(0.08, "lines")) +
  # ── Slope annotation: post-1990 OLS ───────────────────────────────────
  ggplot2::annotate("label",
                    x          = x_ann_p90,
                    y          = y_ann_p90,
                    label      = paste0("Post-1990 OLS\n", ann_p90),
                    hjust      = 0, vjust = 1,
                    size       = 2.7,
                    fontface   = "plain",
                    colour     = COL_THW,
                    fill       = "white",
                    label.size = 0.2,
                    label.padding = ggplot2::unit(0.15, "lines"),
                    label.r    = ggplot2::unit(0.08, "lines")) +
  # ── 95% CI ribbon legend entry (manual) ───────────────────────────────
  # Handled via subtitle text rather than a separate legend entry to keep
  # the legend clean.
  ggplot2::scale_x_continuous(
    breaks = seq(1950, 2025, by = 10),
    expand = ggplot2::expansion(add = c(1, 2))) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::coord_cartesian(ylim = c(ylo_f6, yhi_f6))+
  ggplot2::labs(
    title    = "Increasing thermodynamic fraction of JJA drought severity — Nechako Basin",
    subtitle = paste0(
      "Annual mean |SPEI₀| / |SPEIₚₘ| for June–August, SPEI3, 1950–2025.   ",
      "Solid red = full-record OLS ± 95% CI ribbon.   ",
      "Dashed red = post-1990 OLS ± 90% CI ribbon."),
    x        = "Year",
    y        = expression(italic(F)[thm]*"  (thermodynamic fraction of drought severity)"),
    caption  = paste0(
      "Fₜₕₘ = mean |SPEI₀| / (|SPEIₚₘ| + ε);   ",
      "ε = 10⁻⁶ prevents division by zero.   ",
      "OLS fitted by ordinary least squares; CI from classical normal theory.")) +
  theme_ms +
  ggplot2::theme(
    legend.position  = "right",
    legend.key.size  = ggplot2::unit(0.35, "cm"),
    legend.text      = ggplot2::element_text(size = 8),
    legend.title     = ggplot2::element_text(size = 8.5, face = "bold"),
    axis.title.y     = ggplot2::element_text(size = 9),
    plot.caption     = ggplot2::element_text(size = 6.5, colour = "grey50",
                                             hjust = 0,
                                             margin = ggplot2::margin(t = 4))
  )
fig6_pdf <- file.path(out_dir, "Fig6_thermodynamic_trend.pdf")
fig6_png <- file.path(out_dir, "Fig6_thermodynamic_trend.png")
tryCatch({
  ggplot2::ggsave(fig6_pdf, fig6,
                  width = 7.5, height = 5.0, units = "in", device = "pdf")
  cat(sprintf("  ✓ Manuscript Fig6 (PDF): %s\n", basename(fig6_pdf)))
}, error = function(e) cat(sprintf("  ⚠ Fig6 PDF failed: %s\n", e$message)))
tryCatch({
  ggplot2::ggsave(fig6_png, fig6,
                  width = 7.5, height = 5.0, units = "in",
                  dpi = 300, device = "png")
  cat(sprintf("  ✓ Manuscript Fig6 (PNG): %s\n", basename(fig6_png)))
}, error = function(e) cat(sprintf("  ⚠ Fig6 PNG failed: %s\n", e$message)))
# ─────────────────────────────────────────────────────────────────────────────
# PATCHED FIGURE 6 (slope-only annotations, NO ribbons, NO OLS lines)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n  Building Fig6 PATCHED (slope-only, no ribbons/lines)...\n")
# Annotation: slope only (no p-value, no CI mention)
ann_full_A <- sprintf("%+.4f yr⁻¹\n(1950–2025)", thm_slope)
ann_p90_A  <- sprintf("%+.4f yr⁻¹\n(1990–2025)", thm_slope_p90)
fig6_patched <- ggplot2::ggplot() +
  # 2022-2025 highlight band
  ggplot2::annotate("rect",
                    xmin = 2022, xmax = 2026,
                    ymin = -Inf, ymax = Inf,
                    fill = "grey82", alpha = 0.55) +
  # Annual data points coloured by decade (NO ribbons, NO OLS lines)
  ggplot2::geom_point(
    data        = trend_df,
    ggplot2::aes(x = year, y = thm_frac_abs, fill = decade),
    shape  = 21, size = 2.2, colour = "white", stroke = 0.4,
    inherit.aes = FALSE) +
  ggplot2::scale_fill_brewer(palette = "RdYlBu", direction = -1, name = "Decade") +
  # 1990 reference line
  ggplot2::geom_vline(xintercept = 1990, linetype = "dotted",
                      colour = "grey30", linewidth = 0.6) +
  ggplot2::annotate("text", x = 1991, y = ylo_f6 + 0.03,
                    label = "Post-1990\nacceleration",
                    hjust = 0, vjust = 0, size = 2.6,
                    colour = "grey30", fontface = "italic") +
  # 2022 drought onset line
  ggplot2::geom_vline(xintercept = 2022, linetype = "dashed",
                      colour = "grey50", linewidth = 0.5) +
  ggplot2::annotate("text", x = 2021.5, y = ylo_f6 + 0.03,
                    label = "2022–2025\ndrought",
                    hjust = 1, vjust = 0, size = 2.5,
                    colour = "grey40", fontface = "italic") +
  # Slope annotation — full record (slope only)
  ggplot2::annotate("label",
                    x = 1952, y = yhi_f6 - 0.03,
                    label = paste0("Full record OLS\n", ann_full_A),
                    hjust = 0, vjust = 1, size = 2.7,
                    colour = COL_THW, fill = "white",
                    label.size = 0.2,
                    label.padding = ggplot2::unit(0.15, "lines"),
                    label.r = ggplot2::unit(0.08, "lines")) +
  # Slope annotation — post-1990 (slope only)
  ggplot2::annotate("label",
                    x = 1991, y = yhi_f6 - 0.03,
                    label = paste0("Post-1990 OLS\n", ann_p90_A),
                    hjust = 0, vjust = 1, size = 2.7,
                    colour = COL_THW, fill = "white",
                    label.size = 0.2,
                    label.padding = ggplot2::unit(0.15, "lines"),
                    label.r = ggplot2::unit(0.08, "lines")) +
  ggplot2::scale_x_continuous(
    breaks = seq(1950, 2025, by = 10),
    expand = ggplot2::expansion(add = c(1, 2))) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0))) +
  ggplot2::coord_cartesian(ylim = c(ylo_f6, yhi_f6))+
  ggplot2::labs(
    title    = "Increasing thermodynamic fraction of JJA drought severity — Nechako Basin",
    subtitle = paste0(
      "Annual mean |SPEI₀| / |SPEIₚₘ| for June–August, SPEI3, 1950–2025."),
    x = "Year",
    y = expression(italic(F)[thm]*"  (thermodynamic fraction of drought severity)")) +
  theme_ms +
  ggplot2::theme(
    legend.position = "right",
    legend.key.size = ggplot2::unit(0.35, "cm"),
    legend.text     = ggplot2::element_text(size = 8),
    legend.title    = ggplot2::element_text(size = 8.5, face = "bold"))
# Save patched version
f6_pdf_patched  <- file.path(out_dir, "Fig6_thermodynamic_trend_PATCHED.pdf")
f6_png_patched  <- file.path(out_dir, "Fig6_thermodynamic_trend_PATCHED.png")
tryCatch({
  ggplot2::ggsave(f6_pdf_patched, fig6_patched,
                  width = 7.5, height = 5.0, units = "in", device = "pdf")
  cat(sprintf("  ✓ Fig6 patched (PDF): %s\n", basename(f6_pdf_patched)))
}, error = function(e) cat(sprintf("  ⚠ Fig6 patched PDF: %s\n", e$message)))
tryCatch({
  ggplot2::ggsave(f6_png_patched, fig6_patched,
                  width = 7.5, height = 5.0, units = "in", dpi = 300)
  cat(sprintf("  ✓ Fig6 patched (PNG): %s\n", basename(f6_png_patched)))
}, error = function(e) cat(sprintf("  ⚠ Fig6 patched PNG: %s\n", e$message)))

# ==============================================================================
#   Fig6_SPEI123_thermodynamic_trend.R
# ------------------------------------------------------------------------------
#   Generates Fig6_thermodynamic_trend_SPEI123.pdf/.png
# This is a DROP-IN ADDITION for w11_dynamic_thermodynamic_decomp.R.
# Paste this block immediately AFTER the existing "STEP 6: GENERATE FIGURES"
# section (after the existing Fig6 is saved), or source it standalone once
# decomp_all exists in the environment.
# WHAT IT DOES
# For each of SPEI-1, SPEI-2, SPEI-3:
#   1. Computes annual F_thm = mean(|SPEI_Thw| / (|SPEI_PM| + ε))
# using ALL 12 calendar months (short scales respond throughout the year).
# 2. Fits two OLS models:
#   (a) Full record 1950–2025 — slope, 95% CI, p-value
# (b) Post-1990 sub-period  — slope, 90% CI, p-value
# 3. Builds one panel per scale (same visual language as the existing Fig6
#                                for SPEI-3 JJA: decade-coloured scatter, ribbon CIs, vertical
#                                references at 1990 & 2022, slope annotations).
# The three panels are assembled into a single 7.0 × 11.0 in publication
# figure (PDF + 300 DPI PNG).
# REQUIRES (already in environment after STEP 3 of w11):
#   decomp_all  — data.frame(date, year, month, scale,
#                            SPEI_PM, SPEI_Thw, SPEI_Dynamic, ...)
# out_dir     — character, output directory
# COL_THW     — "#e31a1c"  (red, thermodynamic)
# theme_ms    — shared ggplot2 theme
# OUTPUT files written to out_dir:
#   Fig6_thermodynamic_trend_SPEI123.pdf/.png (ORIGINAL with ribbons/lines/p-values)
# Fig6_thermodynamic_trend_SPEI123_PATCHED.pdf/.png (PATCHED slope-only)
# ==============================================================================
# ── Guard: ensure required objects exist ──────────────────────────────────────
if (!exists("decomp_all"))
  stop("Run w11 Steps 1-5 first so decomp_all, out_dir, COL_THW, theme_ms exist.")
if (!exists("out_dir"))   out_dir  <- "decomp_results"
if (!exists("COL_THW"))   COL_THW  <- "#e31a1c"
if (!exists("theme_ms")) {
  theme_ms  <- ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey40", hjust = 0,
                                            margin = ggplot2::margin(b = 4)),
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      legend.text      = ggplot2::element_text(size = 8.5),
      legend.title     = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey94", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(4, 8, 4, 4))
}
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})
cat("\n  Building Fig6: thermodynamic fraction trend — SPEI-1/2/3...\n")
# ── Helper: format p-value ────────────────────────────────────────────────────
.fmt_p <- function(p) {
  if      (p < 0.001) "p < 0.001"
  else if (p < 0.05)  sprintf("p = %.3f", p)
  else if (p < 0.10)  sprintf("p = %.3f†", p)   # dagger = marginal
  else                sprintf("p = %.2f",  p)
}
# ── Scales to plot ────────────────────────────────────────────────────────────
MS_SCALES <- c(1L, 2L, 3L)
# Panel labels
PANEL_LABELS <- c("(a)", "(b)", "(c)")
# ── Build one panel ───────────────────────────────────────────────────────────
# sc       : integer scale (1, 2, or 3)
# pan_lab  : panel label string, e.g. "(a)"
# show_x   : TRUE = show x-axis text (bottom panel only)
# show_leg : TRUE = show decade legend (bottom panel only)
make_f6_panel <- function(sc, pan_lab, show_x = FALSE, show_leg = FALSE) {
  #   ── 1. Annual F_thm — ALL calendar months ──────────────────────────────────
  td <- decomp_all %>%
    dplyr::filter(scale == sc) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      SPEI_PM       = mean(SPEI_PM,  na.rm = TRUE),
      SPEI_Thw      = mean(SPEI_Thw, na.rm = TRUE),
      thm_frac_abs  = mean(abs(SPEI_Thw) / (abs(SPEI_PM) + 1e-6), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(thm_frac_abs), is.finite(thm_frac_abs))
  n_yrs <- nrow(td)
  if (n_yrs < 10) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("Insufficient data\nSPEI-%d", sc),
                               size = 4, colour = "grey50") +
             ggplot2::theme_void())
  }
  #   ── 2. Full-record OLS (1950-2025) ─────────────────────────────────────────
  fit_full  <- lm(thm_frac_abs ~ year, data = td)
  sl_full   <- coef(fit_full)["year"]
  pv_full   <- summary(fit_full)$coefficients["year", 4]
  ci95      <- confint(fit_full, "year", level = 0.95)
  pred_full <- predict(fit_full,
                       newdata   = data.frame(year = td$year),
                       interval  = "confidence", level = 0.95)
  td$fit_full <- pred_full[, "fit"]
  td$lwr95    <- pred_full[, "lwr"]
  td$upr95    <- pred_full[, "upr"]
  #   ── 3. Post-1990 OLS ───────────────────────────────────────────────────────
  td90 <- dplyr::filter(td, year >= 1990)
  fit90 <- NULL; sl90 <- NA_real_; pv90 <- NA_real_
  if (nrow(td90) >= 10) {
    fit90 <- lm(thm_frac_abs ~ year, data = td90)
    sl90  <- coef(fit90)["year"]
    pv90  <- summary(fit90)$coefficients["year", 4]
    pred90 <- predict(fit90,
                      newdata  = data.frame(year = td90$year),
                      interval = "confidence", level = 0.90)
    td90$fit_p90 <- pred90[, "fit"]
    td90$lwr90   <- pred90[, "lwr"]
    td90$upr90   <- pred90[, "upr"]
  }
  #   ── 4. Decade colouring ────────────────────────────────────────────────────
  all_decades <- seq(1950, 2020, by = 10)
  td$decade <- factor(
    paste0(floor(td$year / 10) * 10, "s"),
    levels = paste0(all_decades, "s"))
  #   ── 5. y-axis range ────────────────────────────────────────────────────────
  ylo <- max(0.0,  min(td$thm_frac_abs, na.rm = TRUE) - 0.06)
  yhi <- min(3, quantile(td$thm_frac_abs, 0.99, na.rm = TRUE)  + 0.16)
  #   ── 6. Annotation strings ──────────────────────────────────────────────────
  ann_full <- sprintf(
    "Full record OLS\n%+.4f yr⁻¹  (%s)\n(1950–2025)",
    sl_full, .fmt_p(pv_full))
  ann_p90 <- if (!is.null(fit90))
    sprintf("Post-1990 OLS\n%+.4f yr⁻¹  (%s)\n(1990–2025)",
            sl90, .fmt_p(pv90))
  else NULL
  #   ── 7. Build plot (ORIGINAL with ribbons and lines) ────────────────────────
  p <- ggplot2::ggplot() +
    # 2022-2025 highlight band
    ggplot2::annotate("rect",
                      xmin = 2022, xmax = 2026,
                      ymin = -Inf, ymax = Inf,
                      fill = "grey82", alpha = 0.55) +
    # 95% CI ribbon — full-record OLS
    ggplot2::geom_ribbon(
      data        = td,
      ggplot2::aes(x = year, ymin = lwr95, ymax = upr95),
      fill        = COL_THW,
      alpha       = 0.15,
      inherit.aes = FALSE) +
    # Full-record OLS line
    ggplot2::geom_line(
      data        = td,
      ggplot2::aes(x = year, y = fit_full),
      colour      = COL_THW,
      linewidth   = 0.9,
      linetype    = "solid",
      inherit.aes = FALSE)
  #   Post-1990 ribbons and line (only if fit succeeded)
  if (!is.null(fit90)) {
    p <- p +
      ggplot2::geom_ribbon(
        data        = td90,
        ggplot2::aes(x = year, ymin = lwr90, ymax = upr90),
        fill        = COL_THW,
        alpha       = 0.10,
        inherit.aes = FALSE) +
      ggplot2::geom_line(
        data        = td90,
        ggplot2::aes(x = year, y = fit_p90),
        colour      = COL_THW,
        linewidth   = 0.75,
        linetype    = "dashed",
        alpha       = 0.80,
        inherit.aes = FALSE)
  }
  p <- p +
    # Annual data points (decade colour)
    ggplot2::geom_point(
      data        = td,
      ggplot2::aes(x = year, y = thm_frac_abs, fill = decade),
      shape       = 21,
      size        = 2.0,
      colour      = "white",
      stroke      = 0.35,
      inherit.aes = FALSE) +
    ggplot2::scale_fill_brewer(
      palette   = "RdYlBu",
      direction = -1,
      name      = "Decade",
      drop      = FALSE) +
    # 1990 reference line
    ggplot2::geom_vline(xintercept = 1990, linetype = "dotted",
                        colour = "grey30", linewidth = 0.55) +
    ggplot2::annotate("text", x = 1991, y = ylo + 0.025,
                      label = "Post-1990\nacceleration",
                      hjust = 0, vjust = 0, size = 2.5,
                      colour = "grey30", fontface = "italic") +
    # 2022 drought onset line
    ggplot2::geom_vline(xintercept = 2022, linetype = "dashed",
                        colour = "grey50", linewidth = 0.45) +
    ggplot2::annotate("text", x = 2021.5, y = ylo + 0.025,
                      label = "2022–2025\ndrought",
                      hjust = 1, vjust = 0, size = 2.4,
                      colour = "grey40", fontface = "italic") +
    # Slope annotation — full record
    ggplot2::annotate("label",
                      x = 1952, y = yhi - 0.03,
                      label = ann_full,
                      hjust = 0, vjust = 1, size = 2.6,
                      colour = COL_THW, fill = "white",
                      label.size = 0.20,
                      label.padding = ggplot2::unit(0.14, "lines"),
                      label.r = ggplot2::unit(0.08, "lines")) +
    # Scale label in top-right corner
    ggplot2::annotate("text", x = 2024, y = yhi - 0.02,
                      label = sprintf("%s  SPEI-%d\n(all months)", pan_lab, sc),
                      hjust = 1, vjust = 1, size = 3.0,
                      colour = "grey10", fontface = "bold") +
    ggplot2::scale_x_continuous(
      breaks = seq(1950, 2025, by = 10),
      expand = ggplot2::expansion(add = c(1, 2))) +
    ggplot2::scale_y_continuous(
      limits = c(ylo, yhi),
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      x = if (show_x) "Year" else NULL,
      y = expression(italic(F)[thm])) +
    theme_ms +
    ggplot2::theme(
      legend.position = if (show_leg) "right" else "none",
      legend.key.size = ggplot2::unit(0.32, "cm"),
      legend.text     = ggplot2::element_text(size = 7.5),
      legend.title    = ggplot2::element_text(size = 8, face = "bold"),
      axis.text.x  = if (show_x) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
      axis.ticks.x = if (show_x) ggplot2::element_line() else ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 9),
      axis.title.x = ggplot2::element_text(size = 9))
  #   Optional post-1990 annotation (only on panels b and c)
  if (!is.null(ann_p90) && sc >= 1) {
    p <- p + ggplot2::annotate("label",
                               x = 1991, y = yhi - 0.03,
                               label = ann_p90,
                               hjust = 0, vjust = 1, size = 2.4,
                               colour = COL_THW, fill = "white",
                               label.size = 0.18,
                               label.padding = ggplot2::unit(0.12, "lines"),
                               label.r = ggplot2::unit(0.07, "lines"))
  }
  p
}
# ── Build PATCHED panel function (slope-only, NO ribbons, NO lines) ────────
make_f6_panel_patched <- function(sc, pan_lab, show_x = FALSE, show_leg = FALSE) {
  #   ── 1. Annual F_thm — ALL calendar months ──────────────────────────────────
  td <- decomp_all %>%
    dplyr::filter(scale == sc) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      SPEI_PM       = mean(SPEI_PM,  na.rm = TRUE),
      SPEI_Thw      = mean(SPEI_Thw, na.rm = TRUE),
      thm_frac_abs  = mean(abs(SPEI_Thw) / (abs(SPEI_PM) + 1e-6), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(thm_frac_abs), is.finite(thm_frac_abs))
  n_yrs <- nrow(td)
  if (n_yrs < 10) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("Insufficient data\nSPEI-%d", sc),
                               size = 4, colour = "grey50") +
             ggplot2::theme_void())
  }
  #   ── 2. OLS fits (for slope annotation only — NOT plotted) ─────────────────
  fit_full <- lm(thm_frac_abs ~ year, data = td)
  sl_full  <- coef(fit_full)["year"]
  pv_full  <- summary(fit_full)$coefficients["year", 4]
  td90 <- dplyr::filter(td, year >= 1990)
  sl90 <- NA_real_
  if (nrow(td90) >= 10) {
    fit90 <- lm(thm_frac_abs ~ year, data = td90)
    sl90  <- coef(fit90)["year"]
  }
  #   ── 3. Decade colouring ────────────────────────────────────────────────────
  all_decades <- seq(1950, 2020, by = 10)
  td$decade <- factor(
    paste0(floor(td$year / 10) * 10, "s"),
    levels = paste0(all_decades, "s"))
  #   ── 4. y-axis range ────────────────────────────────────────────────────────
  ylo <- max(0.0,  min(td$thm_frac_abs, na.rm = TRUE) - 0.06)
  yhi <- min(3, quantile(td$thm_frac_abs, 0.99, na.rm = TRUE)  + 0.16)
  #   ── 5. Annotation strings — SLOPE ONLY, no p-value ────────────────────────
  ann_full <- sprintf("Full record OLS\n%+.4f yr⁻¹\n(1950–2025)", sl_full)
  ann_p90  <- if (!is.na(sl90))
    sprintf("Post-1990 OLS\n%+.4f yr⁻¹\n(1990–2025)", sl90)
  else NULL
  #   ── 6. Build plot (PATCHED: NO ribbons, NO OLS lines) ─────────────────────
  p <- ggplot2::ggplot() +
    # 2022-2025 highlight band
    ggplot2::annotate("rect",
                      xmin = 2022, xmax = 2026,
                      ymin = -Inf, ymax = Inf,
                      fill = "grey82", alpha = 0.55) +
    # Annual data points (NO ribbons, NO OLS lines)
    ggplot2::geom_point(
      data = td,
      ggplot2::aes(x = year, y = thm_frac_abs, fill = decade),
      shape = 21, size  = 2.0, colour = "white", stroke = 0.35,
      inherit.aes = FALSE) +
    ggplot2::scale_fill_brewer(
      palette = "RdYlBu", direction = -1, name = "Decade", drop = FALSE) +
    # 1990 reference line
    ggplot2::geom_vline(xintercept = 1990, linetype = "dotted",
                        colour = "grey30", linewidth = 0.55) +
    ggplot2::annotate("text", x = 1991, y = ylo + 0.025,
                      label = "Post-1990\nacceleration",
                      hjust = 0, vjust = 0, size = 2.5,
                      colour = "grey30", fontface = "italic") +
    # 2022 drought onset line
    ggplot2::geom_vline(xintercept = 2022, linetype = "dashed",
                        colour = "grey50", linewidth = 0.45) +
    ggplot2::annotate("text", x = 2021.5, y = ylo + 0.025,
                      label = "2022–2025\ndrought",
                      hjust = 1, vjust = 0, size = 2.4,
                      colour = "grey40", fontface = "italic") +
    # Slope annotation — full record (slope only)
    ggplot2::annotate("label",
                      x = 1952, y = yhi - 0.03,
                      label = ann_full,
                      hjust = 0, vjust = 1, size = 2.6,
                      colour = COL_THW, fill = "white",
                      label.size = 0.20,
                      label.padding = ggplot2::unit(0.14, "lines"),
                      label.r = ggplot2::unit(0.08, "lines")) +
    # Scale label
    ggplot2::annotate("text", x = 2024, y = yhi - 0.02,
                      label = sprintf("%s  SPEI-%d\n(all months)", pan_lab, sc),
                      hjust = 1, vjust = 1, size = 3.0,
                      colour = "grey10", fontface = "bold") +
    ggplot2::scale_x_continuous(
      breaks = seq(1950, 2025, by = 10),
      expand = ggplot2::expansion(add = c(1, 2))) +
    ggplot2::scale_y_continuous(
      limits = c(ylo, yhi),
      breaks = seq(0, 1, by = 0.1),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      x = if (show_x) "Year" else NULL,
      y = expression(italic(F)[thm])) +
    theme_ms +
    ggplot2::theme(
      legend.position = if (show_leg) "right" else "none",
      legend.key.size = ggplot2::unit(0.32, "cm"),
      axis.text.x  = if (show_x) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
      axis.ticks.x = if (show_x) ggplot2::element_line() else ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 9),
      axis.title.x = ggplot2::element_text(size = 9))
  # Post-1990 slope annotation on panels b and c only
  if (!is.null(ann_p90) && sc >= 1) {
    p <- p + ggplot2::annotate("label",
                               x = 1991, y = yhi - 0.03,
                               label = ann_p90,
                               hjust = 0, vjust = 1, size = 2.4,
                               colour = COL_THW, fill = "white",
                               label.size = 0.18,
                               label.padding = ggplot2::unit(0.12, "lines"),
                               label.r = ggplot2::unit(0.07, "lines"))
  }
  p
}
# ── Build all three panels (ORIGINAL) ─────────────────────────────────────────
panels <- mapply(
  FUN      = make_f6_panel,
  sc       = MS_SCALES,
  pan_lab  = PANEL_LABELS,
  show_x   = c(FALSE, FALSE, TRUE),
  show_leg = c(FALSE, FALSE, TRUE),
  SIMPLIFY = FALSE)
fig6_123  <- patchwork::wrap_plots(panels, ncol = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = paste0("Thermodynamic fraction of drought severity — ",
                      "SPEI-1, 2, 3  —  Nechako River Basin (1950–2025)"),
    subtitle = paste0(
      expression(italic(F)[thm]),
      " = annual mean |SPEI₀| / (|SPEIₚₘ| + ε)  using all 12 calendar months.   ",
      "Solid red = full-record OLS ± 95% CI ribbon.   ",
      "Dashed red = post-1990 OLS ± 90% CI ribbon.   ",
      "Coloured points = decade.  Grey band = 2022–2025 focus period."),
    caption  = paste0(
      "ε = 10⁻⁶ prevents division by zero.   ",
      "OLS fitted by ordinary least squares; CI from classical normal theory.   ",
      "All months included (SPEI-1/2/3 respond throughout the year)."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 8,  colour = "grey35", hjust = 0,
                                            margin = ggplot2::margin(b = 4)),
      plot.caption  = ggplot2::element_text(size = 6.5, colour = "grey50", hjust = 0,
                                            margin = ggplot2::margin(t = 4)),
      legend.position = "right",
      legend.text     = ggplot2::element_text(size = 8),
      legend.title    = ggplot2::element_text(size = 8.5, face = "bold")))
# ── Build all three panels (PATCHED) ──────────────────────────────────────────
panels_patched <- mapply(
  FUN      = make_f6_panel_patched,
  sc       = MS_SCALES,
  pan_lab  = PANEL_LABELS,
  show_x   = c(FALSE, FALSE, TRUE),
  show_leg = c(FALSE, FALSE, TRUE),
  SIMPLIFY = FALSE)
fig6_123_patched  <- patchwork::wrap_plots(panels_patched, ncol = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = paste0("Thermodynamic fraction of drought severity — ",
                      "SPEI-1, 2, 3  —  Nechako River Basin (1950–2025)"),
    subtitle = paste0(
      expression(italic(F)[thm]),
      " = annual mean |SPEI₀| / (|SPEIₚₘ| + ε)  using all 12 calendar months.   ",
      "Coloured points = decade.  Grey band = 2022–2025 focus period."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 8,  colour = "grey35", hjust = 0,
                                            margin = ggplot2::margin(b = 4)),
      legend.position = "right",
      legend.text     = ggplot2::element_text(size = 8),
      legend.title    = ggplot2::element_text(size = 8.5, face = "bold")))
# ── Save ORIGINAL ────────────────────────────────────────────────────────────
f6_123_pdf  <- file.path(out_dir, "Fig6_thermodynamic_trend_SPEI123.pdf")
f6_123_png  <- file.path(out_dir, "Fig6_thermodynamic_trend_SPEI123.png")
tryCatch({
  ggplot2::ggsave(f6_123_pdf, fig6_123,
                  width  = 7.0, height = 11.0, units = "in", device = "pdf")
  cat(sprintf("  ✓ Fig6 SPEI-1/2/3 (PDF): %s\n", basename(f6_123_pdf)))
}, error = function(e) cat(sprintf("  ⚠ PDF failed: %s\n", e$message)))
tryCatch({
  ggplot2::ggsave(f6_123_png, fig6_123,
                  width  = 7.0, height = 11.0, units = "in", dpi = 300, device = "png")
  cat(sprintf("  ✓ Fig6 SPEI-1/2/3 (PNG): %s\n", basename(f6_123_png)))
}, error = function(e) cat(sprintf("  ⚠ PNG failed: %s\n", e$message)))
# ── Save PATCHED ────────────────────────────────────────────────────────────
f6_123_pdf_patched  <- file.path(out_dir, "Fig6_thermodynamic_trend_SPEI123_PATCHED.pdf")
f6_123_png_patched  <- file.path(out_dir, "Fig6_thermodynamic_trend_SPEI123_PATCHED.png")
tryCatch({
  ggplot2::ggsave(f6_123_pdf_patched, fig6_123_patched,
                  width  = 7.0, height = 11.0, units = "in", device = "pdf")
  cat(sprintf("  ✓ Fig6 SPEI-1/2/3 PATCHED (PDF): %s\n", basename(f6_123_pdf_patched)))
}, error = function(e) cat(sprintf("  ⚠ PATCHED PDF failed: %s\n", e$message)))
tryCatch({
  ggplot2::ggsave(f6_123_png_patched, fig6_123_patched,
                  width  = 7.0, height = 11.0, units = "in", dpi = 300, device = "png")
  cat(sprintf("  ✓ Fig6 SPEI-1/2/3 PATCHED (PNG): %s\n", basename(f6_123_png_patched)))
}, error = function(e) cat(sprintf("  ⚠ PATCHED PNG failed: %s\n", e$message)))
# ── Console summary ───────────────────────────────────────────────────────────
cat("\n  Thermodynamic fraction trend summary (all months):\n")
cat(sprintf("  %-8s  %-12s  %-10s  %-12s  %-10s\n",
            "Scale", "Slope full", "p full", "Slope post90", "p post90"))
cat("   ", paste(rep("-", 58), collapse = ""), "\n", sep = "")
for (sc in MS_SCALES) {
  td_sc  <- decomp_all %>%
    dplyr::filter(scale == sc) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      thm_frac_abs = mean(abs(SPEI_Thw) / (abs(SPEI_PM) + 1e-6), na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::filter(!is.na(thm_frac_abs), is.finite(thm_frac_abs))
  if (nrow(td_sc) < 10) { cat(sprintf("  SPEI-%d: insufficient data\n", sc)); next }
  f_all  <- lm(thm_frac_abs ~ year, data = td_sc)
  f_p90  <- lm(thm_frac_abs ~ year, data = dplyr::filter(td_sc, year >= 1990))
  sl_a   <- coef(f_all)["year"];   pv_a <- summary(f_all)$coefficients["year", 4]
  sl_p90 <- coef(f_p90)["year"];   pv_p <- summary(f_p90)$coefficients["year", 4]
  sig_a   <- if (pv_a   < 0.05)  "*" else if (pv_a   < 0.10)  "†" else "  "
  sig_p   <- if (pv_p   < 0.05)  "*" else if (pv_p   < 0.10)  "†" else "  "
  cat(sprintf("  SPEI-%-3d  %+.5f/yr   %.4f%s   %+.5f/yr   %.4f%s\n",
              sc, sl_a, pv_a, sig_a, sl_p90, pv_p, sig_p))
}
cat("  (* p < 0.05  † p < 0.10)\n")
cat("  Output: ", normalizePath(out_dir), "\n")

# ==============================================================================
#   STEP 7: SEASONAL DECOMPOSITION  2022-2025  (SPEI-1 / 2 / 3)
# ------------------------------------------------------------------------------
#   Filters decomp_full_timeseries.csv (= decomp_all in memory) for
# year ∈ 2022–2025  AND  scale ∈ {1, 2, 3}
# Groups by (scale, month) and computes mean ± 1 SD for
# SPEI_PM, SPEI_Thw (thermodynamic), SPEI_Dynamic (dynamic)
# Outputs two companion figures:
#   Fig5c_seasonal_decomp_2022_2025.pdf/.png
# 3-panel line chart (one row per scale).
# X-axis: calendar month (Jan–Dec).
# Lines: SPEI_PM (blue), SPEI_Thw/SPEI₀ (red), SPEI_Dynamic/SPEI_Δ (orange).
# Shaded ribbons: ±1 SD across the four 2022-2025 observations of each month.
# This answers "in which season was the thermodynamic (temperature) component
# most responsible for the 2022-2025 deficit?"
# Fig5d_seasonal_decomp_stacked_2022_2025.pdf/.png
# Companion stacked bar chart (one panel per scale).
# For each month: two bars — thermodynamic vs dynamic absolute contribution.
# Heights = mean(|SPEI_Thw|) and mean(|SPEI_Dynamic|).
# A horizontal reference line at the mean(|SPEI_PM|) for that month.
# This answers "how much of the total SPEI deficit is each component?"
# ==============================================================================
cat("\n===== STEP 7: SEASONAL DECOMPOSITION 2022–2025 (SPEI-1/2/3) =====\n")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
# ── Ensure COL_* and theme_ms exist (guards for standalone sourcing) ────────
if (!exists("COL_PM"))      COL_PM       <- "#2166ac"   # blue  – SPEI_PM
if (!exists("COL_THW"))     COL_THW      <- "#e31a1c"   # red   – SPEI_Thw
if (!exists("COL_DYN"))     COL_DYN      <- "#ff7f00"   # orange– SPEI_Dynamic
if (!exists("out_dir"))     out_dir      <- "decomp_results"
if (!exists("decomp_all"))  decomp_all   <- read.csv(
  file.path(out_dir,  "decomp_full_timeseries.csv"),
  stringsAsFactors = FALSE)
MONTH_ABBR  <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
SEAS_SCALES  <- c(1L, 2L, 3L)
PANEL_LAB_7  <- c("(a) SPEI-1",  "(b) SPEI-2",  "(c) SPEI-3")
# ── Aggregate: mean ± SD per scale × month across 2022-2025 ─────────────────
seas_stats <- decomp_all %>%
  dplyr::filter(year >= 2022, year <= 2025, scale %in% SEAS_SCALES) %>%
  dplyr::mutate(month_abb = factor(MONTH_ABBR[month], levels = MONTH_ABBR)) %>%
  dplyr::group_by(scale, month, month_abb) %>%
  dplyr::summarise(
    # means
    pm_mean  = mean(SPEI_PM,       na.rm = TRUE),
    thw_mean = mean(SPEI_Thw,      na.rm = TRUE),
    dyn_mean = mean(SPEI_Dynamic,  na.rm = TRUE),
    # ±1 SD
    pm_sd    = sd(SPEI_PM,         na.rm = TRUE),
    thw_sd   = sd(SPEI_Thw,        na.rm = TRUE),
    dyn_sd   = sd(SPEI_Dynamic,    na.rm = TRUE),
    # absolute means for stacked bars
    thw_abs  = mean(abs(SPEI_Thw),      na.rm = TRUE),
    dyn_abs  = mean(abs(SPEI_Dynamic),  na.rm = TRUE),
    pm_abs   = mean(abs(SPEI_PM),       na.rm = TRUE),
    n_obs    = dplyr::n(),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(
    # ribbon bounds
    pm_lo  = pm_mean  - pm_sd,  pm_hi  = pm_mean  + pm_sd,
    thw_lo = thw_mean - thw_sd, thw_hi = thw_mean + thw_sd,
    dyn_lo = dyn_mean - dyn_sd, dyn_hi = dyn_mean + dyn_sd
  )
cat(sprintf("  Seasonal stats: %d rows (%d scales × 12 months)\n",
            nrow(seas_stats), length(SEAS_SCALES)))
# ── Shared theme ──────────────────────────────────────────────────────────────
theme_seas <- ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(size = 9, face = "bold", hjust = 0),
    plot.subtitle    = ggplot2::element_text(size = 7.5, colour = "grey35", hjust = 0,
                                             margin = ggplot2::margin(b = 3)),
    axis.title.x     = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(size = 7.5),
    axis.text.y      = ggplot2::element_text(size = 7.5),
    axis.title.y     = ggplot2::element_text(size = 8),
    legend.position  = "bottom",
    legend.key.size  = ggplot2::unit(0.38, "cm"),
    legend.text      = ggplot2::element_text(size = 7.5),
    legend.title     = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.margin      = ggplot2::margin(3, 6, 2, 4)
  )
# ══════════════════════════════════════════════════════════════════════════════
# Fig5c: 12-month line + ribbon chart
# ══════════════════════════════════════════════════════════════════════════════
cat("\n  Building Fig5c: seasonal line/ribbon chart...\n")
make_seas_line_panel <- function(sc, pan_lab, show_xlab = FALSE,
                                 show_leg  = FALSE) {
  df <- dplyr::filter(seas_stats, scale == sc)
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void())
  #   Reshape to long for a single geom_ribbon + geom_line call each
  #   (avoids messy manual layers for 3 components)
  long_mean <- df %>%
    dplyr::select(month_abb, pm_mean, thw_mean, dyn_mean) %>%
    tidyr::pivot_longer(-month_abb,
                        names_to  = "component",
                        values_to = "mean_val") %>%
    dplyr::mutate(component = dplyr::recode(component,
                                            pm_mean  = "SPEIₚₘ (observed)",
                                            thw_mean = "SPEI₀ (thermodynamic)",
                                            dyn_mean = "SPEI∆ (dynamic)"))
  long_lo <- df %>%
    dplyr::select(month_abb, pm_lo, thw_lo, dyn_lo) %>%
    tidyr::pivot_longer(-month_abb, names_to = "component", values_to = "lo") %>%
    dplyr::mutate(component = dplyr::recode(component,
                                            pm_lo  = "SPEIₚₘ (observed)",
                                            thw_lo = "SPEI₀ (thermodynamic)",
                                            dyn_lo = "SPEI∆ (dynamic)"))
  long_hi <- df %>%
    dplyr::select(month_abb, pm_hi, thw_hi, dyn_hi) %>%
    tidyr::pivot_longer(-month_abb, names_to = "component", values_to = "hi") %>%
    dplyr::mutate(component = dplyr::recode(component,
                                            pm_hi  = "SPEIₚₘ (observed)",
                                            thw_hi = "SPEI₀ (thermodynamic)",
                                            dyn_hi = "SPEI∆ (dynamic)"))
  ribbon_df  <- dplyr::left_join(long_mean, long_lo,  by = c("month_abb", "component")) %>%
    dplyr::left_join(long_hi,              by = c("month_abb", "component")) %>%
    dplyr::mutate(component = factor(component,
                                     levels = c("SPEIₚₘ (observed)",
                                                "SPEI₀ (thermodynamic)",
                                                "SPEI∆ (dynamic)")))
  col_pal <- c(
    "SPEIₚₘ (observed)"    = COL_PM,
    "SPEI₀ (thermodynamic)" = COL_THW,
    "SPEI∆ (dynamic)"      = COL_DYN)
  #   y-axis range padded symmetrically
  ymax <- max(abs(c(df$pm_hi, df$thw_hi, df$dyn_hi)), na.rm = TRUE) * 1.15
  ymax <- max(ymax, 0.5, na.rm = TRUE)
  p  <- ggplot2::ggplot(ribbon_df,
                        ggplot2::aes(x = month_abb, colour = component,
                                     fill = component, group = component)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -0.5, colour = "grey50",
                        linetype = "dashed", linewidth = 0.4) +
    # ±1 SD ribbon
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         alpha = 0.15, colour = NA) +
    # Mean line
    ggplot2::geom_line(ggplot2::aes(y = mean_val),
                       linewidth = 0.85) +
    ggplot2::geom_point(ggplot2::aes(y = mean_val),
                        size = 1.8, shape = 21, fill = "white",
                        stroke = 0.9) +
    ggplot2::scale_colour_manual(values = col_pal) +
    ggplot2::scale_fill_manual(values = col_pal) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.3, 0.3))) +
    ggplot2::scale_y_continuous(limits = c(-ymax, ymax),
                                breaks  = seq(-2, 2, by = 0.5)) +
    ggplot2::labs(
      title    = sprintf("%s  SPEI-%d  (2022–2025 monthly mean ±1 SD)",
                         pan_lab, sc),
      subtitle = "Dashed = event threshold (−0.5).  ±1 SD shading from 4 calendar-year observations.",
      y        = "SPEI value") +
    theme_seas +
    ggplot2::theme(
      axis.text.x  = if (show_xlab) ggplot2::element_text(size=7.5)
      else           ggplot2::element_blank(),
      legend.position = if (show_leg) "bottom" else "none")
  p
}
p5c_1 <- make_seas_line_panel(1L, "(a)", show_xlab = FALSE, show_leg = FALSE)
p5c_2 <- make_seas_line_panel(2L, "(b)", show_xlab = FALSE, show_leg = FALSE)
p5c_3 <- make_seas_line_panel(3L, "(c)", show_xlab = TRUE,  show_leg = TRUE)
fig5c  <- p5c_1 / p5c_2 / p5c_3 +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Seasonal decomposition of 2022–2025 drought — Nechako River Basin",
    subtitle = paste0(
      "Monthly mean SPEI values (2022–2025) separated into observed (SPEIₚₘ), ",
      "thermodynamic (SPEI₀) and dynamic (SPEI∆ = SPEIₚₘ − SPEI₀) components.\n",
      "Shading = ±1 SD across the 4 calendar-year replicates (2022, 2023, 2024, 2025)."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 5)),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 8),
      legend.key.size = ggplot2::unit(0.4, "cm")))
for (ext in c("pdf", "png")) {
  out_5c  <- file.path(out_dir, paste0("Fig5c_seasonal_decomp_2022_2025.", ext))
  tryCatch(
    ggplot2::ggsave(out_5c, fig5c,
                    width  = 7.5, height = 9.5, units = "in",
                    dpi    = if (ext == "png") 300 else NULL,
                    device = ext),
    error = function(e) cat(sprintf("  ⚠ Fig5c %s: %s\n", ext, e$message)))
  cat(sprintf("  ✓ Saved: %s\n", basename(out_5c)))
}
# ══════════════════════════════════════════════════════════════════════════════
# Fig5d: stacked bar chart companion
# For each month × scale: stack of |SPEI_Thw| (thermodynamic, red) and
# |SPEI_Dynamic| (dynamic, orange), with a point + line for |SPEI_PM|.
# Three panels (one per scale) in a single figure.
# ══════════════════════════════════════════════════════════════════════════════
cat("\n  Building Fig5d: stacked contribution bars...\n")
make_seas_bar_panel <- function(sc, pan_lab, show_xlab = FALSE,
                                show_leg  = FALSE) {
  df <- dplyr::filter(seas_stats, scale == sc)
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void())
  #   Long format for stacking
  bar_df <- df %>%
    dplyr::select(month_abb, thw_abs, dyn_abs, pm_abs) %>%
    tidyr::pivot_longer(cols = c(thw_abs, dyn_abs),
                        names_to  = "component",
                        values_to = "contribution") %>%
    dplyr::mutate(component = dplyr::recode(component,
                                            thw_abs = "Thermodynamic  |SPEI₀|",
                                            dyn_abs = "Dynamic  |SPEI∆|"))
  comp_levels <- c("Thermodynamic  |SPEI₀|", "Dynamic  |SPEI∆|")
  bar_df$component <- factor(bar_df$component, levels = comp_levels)
  col_bar <- c("Thermodynamic  |SPEI₀|" = COL_THW,
               "Dynamic  |SPEI∆|"       = COL_DYN)
  ymax <- max(bar_df$contribution[bar_df$component == "Thermodynamic  |SPEI₀|"] +
                bar_df$contribution[bar_df$component == "Dynamic  |SPEI∆|"],
              df$pm_abs, na.rm = TRUE) * 1.18
  ymax <- max(ymax, 0.3, na.rm = TRUE)
  p  <- ggplot2::ggplot(bar_df,
                        ggplot2::aes(x = month_abb, y = contribution,
                                     fill = component)) +
    ggplot2::geom_col(position = "stack", colour = "white",
                      linewidth = 0.25, width = 0.72) +
    # Overlay: mean |SPEI_PM| as a black line + open circles
    ggplot2::geom_line(data  = df,
                       ggplot2::aes(x = month_abb, y = pm_abs,
                                    group = 1),
                       colour = "#1a1a1a", linewidth = 0.7,
                       linetype = "solid", inherit.aes = FALSE) +
    ggplot2::geom_point(data = df,
                        ggplot2::aes(x = month_abb, y = pm_abs),
                        shape = 21, size = 2.0,
                        colour = "#1a1a1a", fill = "white",
                        stroke = 0.75, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = col_bar) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = c(0.3, 0.3))) +
    ggplot2::scale_y_continuous(limits  = c(0, ymax),
                                expand  = ggplot2::expansion(mult = c(0, 0.02)),
                                breaks  = seq(0, ceiling(ymax * 2) / 2, by = 0.1)) +
    ggplot2::labs(
      title    = sprintf("%s  SPEI-%d  (|component| contributions)", pan_lab, sc),
      subtitle = "Bars = stacked |thermodynamic| + |dynamic|.  Line = |SPEIₚₘ| (total deficit).",
      y        = "Mean absolute SPEI value  (2022–2025)") +
    theme_seas +
    ggplot2::theme(
      axis.text.x  = if (show_xlab) ggplot2::element_text(size=7.5)
      else           ggplot2::element_blank(),
      legend.position = if (show_leg) "bottom" else "none")
  p
}
p5d_1 <- make_seas_bar_panel(1L, "(a)", show_xlab = FALSE, show_leg = FALSE)
p5d_2 <- make_seas_bar_panel(2L, "(b)", show_xlab = FALSE, show_leg = FALSE)
p5d_3 <- make_seas_bar_panel(3L, "(c)", show_xlab = TRUE,  show_leg = TRUE)
fig5d  <- p5d_1 / p5d_2 / p5d_3 +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Seasonal absolute contribution of drought components — 2022–2025",
    subtitle = paste0(
      "Stacked bars: mean absolute thermodynamic (|SPEI₀|, red) + ",
      "dynamic (|SPEI∆|, orange) contributions per calendar month.\n",
      "Black line: total observed deficit (|SPEIₚₘ|).",
      "  Values are means across 2022, 2023, 2024, 2025."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 5)),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 8),
      legend.key.size = ggplot2::unit(0.4, "cm")))
for (ext in c("pdf", "png")) {
  out_5d  <- file.path(out_dir, paste0("Fig5d_seasonal_decomp_stacked_2022_2025.", ext))
  tryCatch(
    ggplot2::ggsave(out_5d, fig5d,
                    width  = 7.5, height = 9.5, units = "in",
                    dpi    = if (ext == "png") 300 else NULL,
                    device = ext),
    error = function(e) cat(sprintf("  ⚠ Fig5d %s: %s\n", ext, e$message)))
  cat(sprintf("  ✓ Saved: %s\n", basename(out_5d)))
}
# ── Export seasonal stats CSV ─────────────────────────────────────────────────
write.csv(seas_stats,
          file.path(out_dir, "decomp_seasonal_stats_2022_2025.csv"),
          row.names = FALSE)
cat(sprintf("  ✓ Saved: decomp_seasonal_stats_2022_2025.csv (%d rows)\n",
            nrow(seas_stats)))
# ── Console summary ───────────────────────────────────────────────────────────
cat("\n  2022-2025 seasonal decomp: peak months per component\n")
cat(sprintf("  %-8s  %-5s  %-12s  %-12s  %-12s\n",
            "Scale", "Month", "SPEI_PM", "SPEI_Thw", "SPEI_Dyn"))
cat("   ", paste(rep("-", 56), collapse = ""), "\n", sep = "")
for (sc in SEAS_SCALES) {
  df_sc  <- dplyr::filter(seas_stats, scale == sc)
  pk_pm   <- df_sc$month_abb[which.min(df_sc$pm_mean)]
  pk_thw  <- df_sc$month_abb[which.min(df_sc$thw_mean)]
  pk_dyn  <- df_sc$month_abb[which.min(df_sc$dyn_mean)]
  cat(sprintf("  SPEI-%-3d  %-5s  %-12s  %-12s  %-12s\n",
              sc, "(peak)",
              sprintf("%s (%.2f)", pk_pm, min(df_sc$pm_mean, na.rm = TRUE)),
              sprintf("%s (%.2f)", pk_thw, min(df_sc$thw_mean, na.rm = TRUE)),
              sprintf("%s (%.2f)", pk_dyn, min(df_sc$dyn_mean, na.rm = TRUE))))
}

# ==============================================================================
#   COMPLETE
# ==============================================================================
cat("\n============================================================\n")
cat("DYNAMIC/THERMODYNAMIC DECOMPOSITION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Outputs saved to: %s\n", normalizePath(out_dir)))
cat("\nSUPPLEMENTARY FIGURES:\n")
cat("  fig_decomp_spei3_timeseries.pdf  — SPEI-3 full record (1950–2025)\n")
cat("  fig_annual_decomp_barplot.pdf    — annual SPEI-3 drought-year bars\n")
cat("  fig_thermodynamic_trend.pdf      — JJA SPEI-3 thermodynamic trend\n")
cat("  fig_scale_comparison.pdf         — all scales 1–24 percentage stack\n")
cat("\nMANUSCRIPT FIGURES:\n")
cat("  Fig5a_decomp_timeseries_2020_2025.pdf/.png\n")
cat("    → SPEIpm / SPEI0 / SPEIdyn, scales 1–3, 2020–2025 (Section 4.3)\n")
cat("  Fig5b_scale_comparison_stacked.pdf/.png\n")
cat("    → Stacked % bars, scales 1–3, thermodynamic vs dynamic (Table 3 / Section 4.2)\n")
cat("  Fig6_thermodynamic_trend.pdf/.png           (ORIGINAL with ribbons/lines/p-values)\n")
cat("  Fig6_thermodynamic_trend_PATCHED.pdf/.png   (PATCHED slope-only)\n")
cat("  Fig6_thermodynamic_trend_SPEI123.pdf/.png           (ORIGINAL with ribbons/lines/p-values)\n")
cat("  Fig6_thermodynamic_trend_SPEI123_PATCHED.pdf/.png   (PATCHED slope-only)\n")
cat("  Fig5c_seasonal_decomp_2022_2025.pdf/.png   — seasonal line/ribbon chart\n")
cat("  Fig5d_seasonal_decomp_stacked_2022_2025.pdf/.png — stacked contribution bars\n")
cat("\nINTERPRETATION:\n")
cat("  SPEI_PM       = Full observed drought (Penman-Monteith PET)\n")
cat("  SPEI_Thw      = Thermodynamic component (Thornthwaite PET, temperature-only)\n")
cat("  SPEI_Dynamic  = SPEI_PM − SPEI_Thw = dynamic component\n")
cat("                ≈ precipitation + radiation + wind effects\n")