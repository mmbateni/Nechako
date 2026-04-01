# ==============================================================================
# A2_dynamic_thermodynamic_decomp.R
# DYNAMIC vs. THERMODYNAMIC DROUGHT DECOMPOSITION
# Nechako River Basin, BC — 2022-2025 Drought Study
#
# Loads pre-computed SPEI from 3SPEI_ERALand.R outputs
#   - SPEI_PM  (spei_results_seasonal/)      = full observed drought
#   - SPEI_Thw (spei_results_seasonal_thw/)  = thermodynamic component
#   - Difference ≈ dynamic component
#   Needs the SPEI output directories from 3SPEI_ERALand.R
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
# STEP 1: LOAD PRE-COMPUTED SPEI RESULTS
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
  
  spei_list <- list()
  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
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
scale_data <- lapply(setNames(decomp_scales, paste0("scale_", decomp_scales)), function(sc) {
  list(
    pm  = load_spei_scale(sc, "PM",  spei_pm_dir),
    thw = load_spei_scale(sc, "Thw", spei_thw_dir)
  )
})
spei_pm_all  <- lapply(scale_data, `[[`, "pm")
spei_thw_all <- lapply(scale_data, `[[`, "thw")

cat("✓ SPEI results loaded from both directories\n")

# ==============================================================================
# STEP 2: LOAD METADATA (dates, basin mask)
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
basin <- vect("Spatial/nechakoBound_dissolve.shp")
target_crs <- "EPSG:3005"
if (!same.crs(basin, target_crs)) basin <- project(basin, target_crs)

# Get basin mask from first SPEI file (non-NA pixels)
spei_sample <- rast(file.path(spei_pm_dir, "spei_06_month01_Jan.nc"))
basin_mask <- !is.na(values(spei_sample))
n_pixels <- sum(basin_mask)

cat(sprintf("✓ Metadata loaded: %d months, %d basin pixels\n", length(dates), n_pixels))

# ==============================================================================
# STEP 3: COMPUTE DECOMPOSITION FROM LOADED SPEI
# ==============================================================================
cat("\n===== STEP 3: COMPUTE DYNAMIC/THERMODYNAMIC DECOMPOSITION =====\n")

# Decomposition framework (from 3SPEI_ERALand.R documentation):
#   SPEI_PM         = full observed drought severity
#   SPEI_Thw        = thermodynamic (temperature-driven) component
#   SPEI_PM - SPEI_Thw ≈ dynamic component (precip + radiation + wind)

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

decomp_all <- bind_rows(lapply(decomp_scales, function(sc) {
  pm_mat  <- spei_pm_all[[paste0("scale_", sc)]]
  thw_mat <- spei_thw_all[[paste0("scale_", sc)]]
  if (is.null(pm_mat) || is.null(thw_mat)) {
    cat(sprintf("  Skipping scale %d (missing data)\n", sc)); return(NULL)
  }
  min_time <- min(ncol(pm_mat), ncol(thw_mat))
  result <- compute_decomposition(
    pm_mat[, 1:min_time, drop=FALSE],
    thw_mat[, 1:min_time, drop=FALSE],
    dates[1:min_time], month_nums[1:min_time], year_nums[1:min_time], sc
  )
  cat(sprintf("  ✓ Scale %d decomposition complete\n", sc))
  result
}))
write.csv(decomp_all, file.path(out_dir, "decomp_full_timeseries.csv"), row.names = FALSE)
cat("✓ Decomposition timeseries saved\n")

# ==============================================================================
# STEP 4: 2022-2025 DROUGHT DECOMPOSITION SUMMARY
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

# ==============================================================================
# STEP 5: LONG-TERM TREND IN THERMODYNAMIC FRACTION
# ==============================================================================
# STEP 5: TREND IN THERMODYNAMIC FRACTION (1950-2025)
# Two OLS models:
#   (a) Full record 1950-2025  — slope, 95% CI, p-value
#   (b) Post-1990 sub-period   — slope and p-value for annotation
# Both are computed here and used in manuscript Figure 6.
# ==============================================================================
cat("\n===== STEP 5: TREND IN THERMODYNAMIC FRACTION (1950-2025) =====\n")

# Focus on SPEI-6 JJA (summer drought most relevant for Nechako)
trend_df <- decomp_all %>%
  dplyr::filter(scale == 6, month %in% c(6, 7, 8)) %>%
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
cat(sprintf("  Slope    : %+.5f yr\u207b\u00b9  (%+.4f per decade)\n",
            thm_slope, thm_slope_dec))
cat(sprintf("  95%% CI   : [%+.5f, %+.5f] yr\u207b\u00b9\n",
            thm_ci95[1], thm_ci95[2]))
cat(sprintf("  p-value  : %.4f\n", thm_pval))

# ── Post-1990 sub-period OLS ───────────────────────────────────────────────
trend_post90   <- dplyr::filter(trend_df, year >= 1990)
thm_trend_p90  <- lm(thm_frac_abs ~ year, data = trend_post90)
thm_slope_p90  <- coef(thm_trend_p90)["year"]
thm_pval_p90   <- summary(thm_trend_p90)$coefficients["year", 4]

cat(sprintf("\nPost-1990 sub-period:\n"))
cat(sprintf("  Slope    : %+.5f yr\u207b\u00b9  (%+.4f per decade)\n",
            thm_slope_p90, thm_slope_p90 * 10))
cat(sprintf("  p-value  : %.4f\n", thm_pval_p90))

# Helper: format p-value for annotation
fmt_p <- function(p) {
  if      (p < 0.001) "p < 0.001"
  else if (p < 0.01)  sprintf("p = %.3f", p)
  else if (p < 0.05)  sprintf("p = %.3f", p)
  else if (p < 0.10)  sprintf("p = %.3f", p)
  else                sprintf("p = %.2f",  p)
}

# Annotation strings used in Figure 6
ann_full <- sprintf("%+.4f yr\u207b\u00b9, %s\n(1950\u20132025)",
                    thm_slope, fmt_p(thm_pval))
ann_p90  <- sprintf("%+.4f yr\u207b\u00b9, %s\n(1990\u20132025)",
                    thm_slope_p90, fmt_p(thm_pval_p90))

# ==============================================================================
# STEP 6: GENERATE FIGURES
#
# SUPPLEMENTARY figures (original, kept for reference):
#   fig_decomp_spei6_timeseries.pdf  — SPEI-6 full record decomposition
#   fig_annual_decomp_barplot.pdf    — annual SPEI-6 drought years bar
#   fig_thermodynamic_trend.pdf      — JJA SPEI-6 thermodynamic trend
#   fig_scale_comparison.pdf         — all scales 1-24 percentage stack
#
# MANUSCRIPT figures (new, publication quality):
#   Fig5a_decomp_timeseries_2020_2025.pdf/.png
#       SPEIₚₘ, SPEI₀ (=SPEIᵟʰʷ), and SPEIᵟ on a single axis, scales 1–3,
#       2020–2025 focus window. Core figure for Section 4.3.
#   Fig5b_scale_comparison_stacked.pdf/.png
#       Stacked percentage bar chart, scales 1–3 only, showing thermodynamic
#       and dynamic percentage contributions with data labels. Replaces/
#       accompanies Table 3. Core result of the decomposition framework.
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
# SUPPLEMENTARY FIGURES (original code, lightly tidied)
# ==============================================================================

# S1: SPEI-6 full-record decomposition (1950-2025)
df_s1 <- decomp_all %>% dplyr::filter(scale == 6)
pdf(file.path(out_dir, "fig_decomp_spei6_timeseries.pdf"), width = 14, height = 8)
df_s1 %>%
  dplyr::select(date, SPEI_PM, SPEI_Thw, SPEI_Dynamic) %>%
  tidyr::pivot_longer(cols = -date, names_to = "component", values_to = "SPEI") %>%
  dplyr::mutate(component = dplyr::recode(component,
                                          SPEI_PM      = "Full (PM PET)",
                                          SPEI_Thw     = "Thermodynamic (Thw PET)",
                                          SPEI_Dynamic = "Dynamic (PM \u2212 Thw)")) %>%
  ggplot2::ggplot(ggplot2::aes(x = date, y = SPEI, colour = component)) +
  ggplot2::geom_hline(yintercept = 0, colour = "gray50") +
  ggplot2::geom_hline(yintercept = THR_EVENT, linetype = "dashed",
                      colour = "grey40", alpha = 0.7) +
  ggplot2::geom_line(alpha = 0.8) +
  ggplot2::annotate("rect",
                    xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                    ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "orange") +
  ggplot2::scale_colour_manual(values = c(COL_PM, COL_THW, COL_DYN)) +
  ggplot2::labs(title = "Nechako Basin: SPEI-6 Dynamic vs. Thermodynamic Decomposition (1950\u20132025)",
                subtitle = "Orange = 2022\u20132025. Dashed = event threshold (x\u2080 = \u22120.5)",
                x = NULL, y = "SPEI-6") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")
dev.off()
cat("  \u2713 Supplementary S1: fig_decomp_spei6_timeseries.pdf\n")

# S2: Annual SPEI-6 drought-year stacked bars
drought_years_df <- decomp_all %>%
  dplyr::filter(scale == 6, SPEI_PM < -0.8) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    SPEI_PM       = mean(SPEI_PM,       na.rm = TRUE),
    Thermodynamic = mean(SPEI_Thw,      na.rm = TRUE),
    Dynamic       = mean(SPEI_Dynamic,  na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::filter(abs(SPEI_PM) > 0.3) %>%
  tidyr::pivot_longer(cols = c(Dynamic, Thermodynamic),
                      names_to = "component", values_to = "contribution")

pdf(file.path(out_dir, "fig_annual_decomp_barplot.pdf"), width = 12, height = 6)
ggplot2::ggplot(drought_years_df,
                ggplot2::aes(x = factor(year), y = contribution, fill = component)) +
  ggplot2::geom_bar(stat = "identity", position = "stack") +
  ggplot2::geom_hline(yintercept = 0, colour = "gray40") +
  ggplot2::scale_fill_manual(
    values = c(Dynamic = COL_DYN, Thermodynamic = COL_THW)) +
  ggplot2::labs(title = "Annual Drought Severity Decomposition: SPEI-6 Drought Years",
                x = "Year", y = "SPEI contribution", fill = "Component") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
dev.off()
cat("  \u2713 Supplementary S2: fig_annual_decomp_barplot.pdf\n")

# S3: JJA SPEI-6 thermodynamic trend (supplementary — simple version kept)
pdf(file.path(out_dir, "fig_thermodynamic_trend.pdf"), width = 10, height = 5)
print(
  ggplot2::ggplot(trend_df, ggplot2::aes(x = year, y = thm_frac_abs)) +
    ggplot2::geom_point(alpha = 0.5, colour = COL_THW) +
    ggplot2::geom_smooth(method = "lm", colour = COL_THW,
                         fill = COL_THW, alpha = 0.2) +
    ggplot2::geom_vline(xintercept = 2022, linetype = "dashed",
                        colour = "orange") +
    ggplot2::labs(
      title    = "Growing thermodynamic fraction \u2014 JJA SPEI-6, Nechako Basin",
      subtitle = sprintf("Linear trend: %+.3f per decade (%s)",
                         thm_slope_dec, fmt_p(thm_pval)),
      x = "Year", y = "Thermodynamic fraction of drought severity") +
    ggplot2::theme_minimal(base_size = 12)
)
dev.off()
cat("  \u2713 Supplementary S3: fig_thermodynamic_trend.pdf\n")

# ==============================================================================
# MANUSCRIPT FIGURE 6 — Thermodynamic fraction trend, Section 4.4
#
# Design:
#   • Scatter: each annual JJA Fthm value (filled circle, coloured by decade)
#   • Full-record OLS line + 95% confidence ribbon (solid red)
#   • Post-1990 OLS line + 90% confidence ribbon (dashed, lighter red)
#   • Vertical reference line at 1990 with text annotation
#   • Vertical dashed line at 2022 (start of focus drought)
#   • Slope + p-value text annotations for both fits
#   • 2022-2025 grey highlight band
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
trend_post90$fit_p90 <- predict(thm_trend_p90,
                                newdata = data.frame(year = trend_post90$year),
                                interval = "confidence", level = 0.90)[, "fit"]
trend_post90$lwr90   <- predict(thm_trend_p90,
                                newdata = data.frame(year = trend_post90$year),
                                interval = "confidence", level = 0.90)[, "lwr"]
trend_post90$upr90   <- predict(thm_trend_p90,
                                newdata = data.frame(year = trend_post90$year),
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
                    label  = "2022\u20132025\ndrought",
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
    limits = c(ylo_f6, yhi_f6),
    breaks = seq(0, 1, by = 0.1),
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0))) +
  
  ggplot2::labs(
    title    = "Increasing thermodynamic fraction of JJA drought severity \u2014 Nechako Basin",
    subtitle = paste0(
      "Annual mean |SPEI\u2080| / |SPEI\u209a\u2098| for June\u2013August, SPEI-6, 1950\u20132025.  ",
      "Solid red = full-record OLS \u00b1 95% CI ribbon.  ",
      "Dashed red = post-1990 OLS \u00b1 90% CI ribbon."),
    x        = "Year",
    y        = expression(italic(F)[thm]*"  (thermodynamic fraction of drought severity)"),
    caption  = paste0(
      "F\u209c\u02b0\u2098 = mean |SPEI\u2080| / (|SPEI\u209a\u2098| + \u03b5);  ",
      "\u03b5 = 10\u207b\u2076 prevents division by zero.  ",
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
  cat(sprintf("  \u2713 Manuscript Fig6 (PDF): %s\n", basename(fig6_pdf)))
}, error = function(e) cat(sprintf("  \u26a0 Fig6 PDF failed: %s\n", e$message)))

tryCatch({
  ggplot2::ggsave(fig6_png, fig6,
                  width = 7.5, height = 5.0, units = "in",
                  dpi = 300, device = "png")
  cat(sprintf("  \u2713 Manuscript Fig6 (PNG): %s\n", basename(fig6_png)))
}, error = function(e) cat(sprintf("  \u26a0 Fig6 PNG failed: %s\n", e$message)))

# S4: All-scale percentage comparison (scales 1-24, supplementary overview)
pdf(file.path(out_dir, "fig_scale_comparison.pdf"), width = 12, height = 6)
drought_decomp %>%
  tidyr::pivot_longer(cols = c(thm_pct, dyn_pct),
                      names_to = "component", values_to = "pct") %>%
  dplyr::mutate(component = dplyr::recode(component,
                                          thm_pct = "Thermodynamic", dyn_pct = "Dynamic")) %>%
  ggplot2::ggplot(ggplot2::aes(x = factor(scale), y = pct, fill = component)) +
  ggplot2::geom_bar(stat = "identity", position = "stack") +
  ggplot2::geom_hline(yintercept = 100, linetype = "dashed", colour = "gray40") +
  ggplot2::scale_fill_manual(
    values = c(Dynamic = COL_DYN, Thermodynamic = COL_THW)) +
  ggplot2::labs(title = "Drought component contribution by timescale (2022\u20132025)",
                x = "SPEI timescale (months)", y = "Contribution (%)") +
  ggplot2::theme_minimal(base_size = 12)
dev.off()
cat("  \u2713 Supplementary S4: fig_scale_comparison.pdf\n")

# ==============================================================================
# MANUSCRIPT FIGURE 5a — Decomposition time series, scales 1-3, 2020-2025
# Section 4.3: shows SPEIₚₘ, SPEI₀, and SPEIᵟ on a single axis per scale
# ==============================================================================
cat("\n  Building Fig5a: decomposition time series (2020\u20132025, scales 1\u20133)...\n")

ms_scales_ts <- c(1, 2, 3)

panels_ts <- lapply(ms_scales_ts, function(sc) {
  
  df_sc <- decomp_all %>%
    dplyr::filter(scale == sc,
                  date >= as.Date("2020-01-01"),
                  date <= as.Date("2025-12-31")) %>%
    dplyr::select(date, SPEI_PM, SPEI_Thw, SPEI_Dynamic) %>%
    tidyr::pivot_longer(cols = -date,
                        names_to  = "component",
                        values_to = "value") %>%
    dplyr::mutate(
      component = dplyr::recode(component,
                                SPEI_PM      = "SPEI\u209a\u2098 (full, PM PET)",
                                SPEI_Thw     = "SPEI\u2080 (thermodynamic)",
                                SPEI_Dynamic = "SPEI\u1d5f (dynamic)"),
      component = factor(component,
                         levels = c("SPEI\u209a\u2098 (full, PM PET)",
                                    "SPEI\u2080 (thermodynamic)",
                                    "SPEI\u1d5f (dynamic)"))
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
                      label = "2022\u20132025", size = 2.5,
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
        "SPEI\u209a\u2098 (full, PM PET)"  = COL_PM,
        "SPEI\u2080 (thermodynamic)"       = COL_THW,
        "SPEI\u1d5f (dynamic)"             = COL_DYN)) +
    ggplot2::scale_linetype_manual(
      values = c(
        "SPEI\u209a\u2098 (full, PM PET)"  = "solid",
        "SPEI\u2080 (thermodynamic)"       = "solid",
        "SPEI\u1d5f (dynamic)"             = "dashed")) +
    
    # Threshold label (first panel only)
    {if (sc == min(ms_scales_ts))
      ggplot2::annotate("text",
                        x = as.Date("2020-03-01"), y = THR_EVENT + 0.07,
                        label = "x\u2080 = \u22120.5", size = 2.2,
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
      subtitle = sprintf("Thermodynamic = %.0f%%; Dynamic = %.0f%% of 2022\u20132025 severity",
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

fig5a <- patchwork::wrap_plots(panels_ts, ncol = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Dynamic\u2013thermodynamic SPEI decomposition \u2014 Nechako River Basin, 2020\u20132025",
    subtitle = paste0(
      "Black = SPEI\u209a\u2098 (PM PET, full observed);  ",
      "Red = SPEI\u2080 (Thornthwaite PET, thermodynamic);  ",
      "Blue dashed = SPEI\u1d5f = SPEI\u209a\u2098 \u2212 SPEI\u2080 (dynamic)\n",
      "Grey band = 2022\u20132025 focus period.  Dashed reference = event threshold (x\u2080 = \u22120.5)."),
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
  cat(sprintf("  \u2713 Manuscript Fig5a (PDF): %s\n", basename(fig5a_pdf)))
}, error = function(e) cat(sprintf("  \u26a0 Fig5a PDF failed: %s\n", e$message)))

tryCatch({
  ggplot2::ggsave(fig5a_png, fig5a,
                  width = 7.0, height = 8.5, units = "in",
                  dpi = 300, device = "png")
  cat(sprintf("  \u2713 Manuscript Fig5a (PNG): %s\n", basename(fig5a_png)))
}, error = function(e) cat(sprintf("  \u26a0 Fig5a PNG failed: %s\n", e$message)))

# ==============================================================================
# MANUSCRIPT FIGURE 5b — Stacked percentage bar chart, scales 1-3 only
# Core visual for Table 3 / Section 4.2
# Shows the shift from ~100% thermodynamic (scale-1) to ~65%/35% (scale-3)
# with exact percentage data labels on each segment.
# ==============================================================================
cat("\n  Building Fig5b: stacked percentage bars (scales 1\u20133)...\n")

ms_scales_bar <- c(1, 2, 3)

decomp_ms_bar <- drought_decomp %>%
  dplyr::filter(scale %in% ms_scales_bar) %>%
  dplyr::mutate(
    scale_lab = sprintf("SPEI-%d", scale),
    scale_lab = factor(scale_lab,
                       levels = sprintf("SPEI-%d", ms_scales_bar))
  )

# Long format for stacked bars
bar_long <- decomp_ms_bar %>%
  dplyr::select(scale, scale_lab, thm_pct, dyn_pct) %>%
  tidyr::pivot_longer(cols = c(thm_pct, dyn_pct),
                      names_to  = "component",
                      values_to = "pct") %>%
  dplyr::mutate(
    component = dplyr::recode(component,
                              thm_pct = "Thermodynamic (SPEI\u2080)",
                              dyn_pct = "Dynamic (SPEI\u1d5f)"),
    component = factor(component,
                       levels = c("Thermodynamic (SPEI\u2080)",
                                  "Dynamic (SPEI\u1d5f)"))
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
    values = c("Thermodynamic (SPEI\u2080)" = COL_THW,
               "Dynamic (SPEI\u1d5f)"       = COL_DYN)) +
  
  ggplot2::scale_y_continuous(
    limits = c(0, 115),
    breaks = seq(0, 100, by = 20),
    labels = function(x) sprintf("%d%%", x),
    expand = ggplot2::expansion(mult = c(0, 0))) +
  
  ggplot2::labs(
    title    = "Thermodynamic vs. dynamic contribution to 2022\u20132025 drought severity",
    subtitle = paste0(
      "Percentage of mean |SPEI\u209a\u2098| attributable to each component.\n",
      "Scales 1\u20133 months shown (matching manuscript analysis period in Section 3.2.1)."),
    x        = "Accumulation scale",
    y        = "Contribution to drought severity (%)",
    fill     = NULL,
    caption  = paste0(
      "Thermodynamic = SPEI\u2080 (Thornthwaite PET, temperature-driven);  ",
      "Dynamic = SPEI\u1d5f = SPEI\u209a\u2098 \u2212 SPEI\u2080 (circulation-driven).\n",
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
  cat(sprintf("  \u2713 Manuscript Fig5b (PDF): %s\n", basename(fig5b_pdf)))
}, error = function(e) cat(sprintf("  \u26a0 Fig5b PDF failed: %s\n", e$message)))

tryCatch({
  ggplot2::ggsave(fig5b_png, fig5b,
                  width = 5.5, height = 5.0, units = "in",
                  dpi = 300, device = "png")
  cat(sprintf("  \u2713 Manuscript Fig5b (PNG): %s\n", basename(fig5b_png)))
}, error = function(e) cat(sprintf("  \u26a0 Fig5b PNG failed: %s\n", e$message)))

# ==============================================================================
# Fig6_SPEI123_thermodynamic_trend.R
# ------------------------------------------------------------------------------
# Generates Fig6_thermodynamic_trend_SPEI123.pdf/.png
#
# This is a DROP-IN ADDITION for w11_dynamic_thermodynamic_decomp.R.
# Paste this block immediately AFTER the existing "STEP 6: GENERATE FIGURES"
# section (after the existing Fig6 is saved), or source it standalone once
# `decomp_all` exists in the environment.
#
# WHAT IT DOES
#   For each of SPEI-1, SPEI-2, SPEI-3:
#     1. Computes annual F_thm = mean(|SPEI_Thw| / (|SPEI_PM| + ε))
#        using ALL 12 calendar months (short scales respond throughout the year).
#     2. Fits two OLS models:
#          (a) Full record 1950–2025 — slope, 95% CI, p-value
#          (b) Post-1990 sub-period  — slope, 90% CI, p-value
#     3. Builds one panel per scale (same visual language as the existing Fig6
#        for SPEI-6 JJA: decade-coloured scatter, ribbon CIs, vertical
#        references at 1990 & 2022, slope annotations).
#   The three panels are assembled into a single 7.0 × 11.0 in publication
#   figure (PDF + 300 DPI PNG).
#
# REQUIRES (already in environment after STEP 3 of w11):
#   decomp_all  — data.frame(date, year, month, scale,
#                             SPEI_PM, SPEI_Thw, SPEI_Dynamic, ...)
#   out_dir     — character, output directory
#   COL_THW     — "#e31a1c"  (red, thermodynamic)
#   theme_ms    — shared ggplot2 theme
#
# OUTPUT files written to out_dir:
#   Fig6_thermodynamic_trend_SPEI123.pdf
#   Fig6_thermodynamic_trend_SPEI123.png
# ==============================================================================

# ── Guard: ensure required objects exist ──────────────────────────────────────
if (!exists("decomp_all"))
  stop("decomp_all not found.  Run Steps 1-3 of w11_dynamic_thermodynamic_decomp.R first.")
if (!exists("out_dir"))   out_dir <- "decomp_results"
if (!exists("COL_THW"))   COL_THW <- "#e31a1c"
if (!exists("theme_ms")) {
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
  else if (p < 0.10)  sprintf("p = %.3f\u2020", p)   # dagger = marginal
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
  
  # ── 1. Annual F_thm — ALL calendar months ──────────────────────────────────
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
  
  # ── 2. Full-record OLS (1950-2025) ─────────────────────────────────────────
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
  
  # ── 3. Post-1990 OLS ───────────────────────────────────────────────────────
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
  
  # ── 4. Decade colouring ────────────────────────────────────────────────────
  all_decades <- seq(1950, 2020, by = 10)
  td$decade <- factor(
    paste0(floor(td$year / 10) * 10, "s"),
    levels = paste0(all_decades, "s"))
  
  # ── 5. y-axis range ────────────────────────────────────────────────────────
  ylo <- max(0.0,  min(td$thm_frac_abs, na.rm = TRUE) - 0.06)
  yhi <- min(1.05, max(td$thm_frac_abs, na.rm = TRUE) + 0.16)
  
  # ── 6. Annotation strings ──────────────────────────────────────────────────
  ann_full <- sprintf(
    "Full record OLS\n%+.4f yr\u207b\u00b9  (%s)\n(1950\u20132025)",
    sl_full, .fmt_p(pv_full))
  
  ann_p90 <- if (!is.null(fit90))
    sprintf("Post-1990 OLS\n%+.4f yr\u207b\u00b9  (%s)\n(1990\u20132025)",
            sl90, .fmt_p(pv90))
  else NULL
  
  # ── 7. Build plot ──────────────────────────────────────────────────────────
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
  
  # Post-1990 ribbons and line (only if fit succeeded)
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
    ggplot2::geom_vline(
      xintercept = 1990,
      linetype   = "dotted",
      colour     = "grey30",
      linewidth  = 0.55) +
    ggplot2::annotate("text",
                      x      = 1991,
                      y      = ylo + 0.025,
                      label  = "Post-1990\nacceleration",
                      hjust  = 0, vjust = 0,
                      size   = 2.5,
                      colour = "grey30",
                      fontface = "italic") +
    
    # 2022 drought onset line
    ggplot2::geom_vline(
      xintercept = 2022,
      linetype   = "dashed",
      colour     = "grey50",
      linewidth  = 0.45) +
    ggplot2::annotate("text",
                      x      = 2021.5,
                      y      = ylo + 0.025,
                      label  = "2022\u20132025\ndrought",
                      hjust  = 1, vjust = 0,
                      size   = 2.4,
                      colour = "grey40",
                      fontface = "italic") +
    
    # Slope annotation — full record
    ggplot2::annotate("label",
                      x             = 1952,
                      y             = yhi - 0.03,
                      label         = ann_full,
                      hjust         = 0, vjust = 1,
                      size          = 2.6,
                      colour        = COL_THW,
                      fill          = "white",
                      label.size    = 0.20,
                      label.padding = ggplot2::unit(0.14, "lines"),
                      label.r       = ggplot2::unit(0.08, "lines")) +
    
    # Scale label in top-right corner
    ggplot2::annotate("text",
                      x      = 2024,
                      y      = yhi - 0.02,
                      label  = sprintf("%s  SPEI-%d\n(all months)", pan_lab, sc),
                      hjust  = 1, vjust = 1,
                      size   = 3.0,
                      colour = "grey10",
                      fontface = "bold") +
    
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
      axis.text.x     = if (show_x) ggplot2::element_text(size = 8)
      else               ggplot2::element_blank(),
      axis.ticks.x    = if (show_x) ggplot2::element_line()
      else               ggplot2::element_blank(),
      axis.title.y    = ggplot2::element_text(size = 9),
      axis.title.x    = ggplot2::element_text(size = 9))
  
  # Optional post-1990 annotation (only on middle or bottom panel to avoid clutter)
  if (!is.null(fit90) && !is.null(ann_p90) && sc >= 2) {
    p <- p +
      ggplot2::annotate("label",
                        x             = 1991,
                        y             = yhi - 0.03,
                        label         = ann_p90,
                        hjust         = 0, vjust = 1,
                        size          = 2.4,
                        colour        = COL_THW,
                        fill          = "white",
                        label.size    = 0.18,
                        label.padding = ggplot2::unit(0.12, "lines"),
                        label.r       = ggplot2::unit(0.07, "lines"))
  }
  
  p
}

# ── Build all three panels ─────────────────────────────────────────────────────
panels <- mapply(
  FUN = make_f6_panel,
  sc       = MS_SCALES,
  pan_lab  = PANEL_LABELS,
  show_x   = c(FALSE, FALSE, TRUE),
  show_leg = c(FALSE, FALSE, TRUE),
  SIMPLIFY = FALSE)

# ── Assemble with patchwork ───────────────────────────────────────────────────
fig6_123 <- patchwork::wrap_plots(panels, ncol = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = paste0(
      "Thermodynamic fraction of drought severity \u2014",
      " SPEI-1, 2, 3  \u2014  Nechako River Basin (1950\u20132025)"),
    subtitle = paste0(
      expression(italic(F)[thm]),
      " = annual mean |SPEI\u2080| / (|SPEI\u209a\u2098| + \u03b5)  using all 12 calendar months.  ",
      "Solid red = full-record OLS \u00b1 95% CI ribbon.  ",
      "Dashed red = post-1990 OLS \u00b1 90% CI ribbon.  ",
      "Coloured points = decade.  Grey band = 2022\u20132025 focus period."),
    caption  = paste0(
      "\u03b5 = 10\u207b\u2076 prevents division by zero.  ",
      "OLS fitted by ordinary least squares; CI from classical normal theory.  ",
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

# ── Save PDF + PNG ────────────────────────────────────────────────────────────
f6_pdf <- file.path(out_dir, "Fig6_thermodynamic_trend_SPEI123.pdf")
f6_png <- file.path(out_dir, "Fig6_thermodynamic_trend_SPEI123.png")

tryCatch({
  ggplot2::ggsave(f6_pdf, fig6_123,
                  width  = 7.0, height = 11.0, units = "in", device = "pdf")
  cat(sprintf("  \u2713 Fig6 SPEI-1/2/3 (PDF): %s\n", basename(f6_pdf)))
}, error = function(e) cat(sprintf("  \u26a0 PDF failed: %s\n", e$message)))

tryCatch({
  ggplot2::ggsave(f6_png, fig6_123,
                  width  = 7.0, height = 11.0, units = "in", dpi = 300, device = "png")
  cat(sprintf("  \u2713 Fig6 SPEI-1/2/3 (PNG): %s\n", basename(f6_png)))
}, error = function(e) cat(sprintf("  \u26a0 PNG failed: %s\n", e$message)))

# ── Console summary ───────────────────────────────────────────────────────────
cat("\n  Thermodynamic fraction trend summary (all months):\n")
cat(sprintf("  %-8s  %-12s  %-10s  %-12s  %-10s\n",
            "Scale", "Slope full", "p full", "Slope post90", "p post90"))
cat("  ", paste(rep("-", 58), collapse = ""), "\n", sep = "")
for (sc in MS_SCALES) {
  td_sc <- decomp_all %>%
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
  
  sig_a  <- if (pv_a  < 0.05) "*" else if (pv_a  < 0.10) "\u2020" else " "
  sig_p  <- if (pv_p  < 0.05) "*" else if (pv_p  < 0.10) "\u2020" else " "
  
  cat(sprintf("  SPEI-%-3d  %+.5f/yr   %.4f%s   %+.5f/yr   %.4f%s\n",
              sc, sl_a, pv_a, sig_a, sl_p90, pv_p, sig_p))
}
cat("  (* p < 0.05  \u2020 p < 0.10)\n")
cat("  Output: ", normalizePath(out_dir), "\n")

# ==============================================================================
# STEP 7: SEASONAL DECOMPOSITION  2022-2025  (SPEI-1 / 2 / 3)
# ------------------------------------------------------------------------------
# Filters decomp_full_timeseries.csv (= decomp_all in memory) for
#   year ∈ 2022–2025  AND  scale ∈ {1, 2, 3}
# Groups by (scale, month) and computes mean ± 1 SD for
#   SPEI_PM, SPEI_Thw (thermodynamic), SPEI_Dynamic (dynamic)
#
# Outputs two companion figures:
#
#   Fig5c_seasonal_decomp_2022_2025.pdf/.png
#       3-panel line chart (one row per scale).
#       X-axis: calendar month (Jan–Dec).
#       Lines: SPEI_PM (blue), SPEI_Thw/SPEI₀ (red), SPEI_Dynamic/SPEI_Δ (orange).
#       Shaded ribbons: ±1 SD across the four 2022-2025 observations of each month.
#       This answers "in which season was the thermodynamic (temperature) component
#       most responsible for the 2022-2025 deficit?"
#
#   Fig5d_seasonal_decomp_stacked_2022_2025.pdf/.png
#       Companion stacked bar chart (one panel per scale).
#       For each month: two bars — thermodynamic vs dynamic absolute contribution.
#       Heights = mean(|SPEI_Thw|) and mean(|SPEI_Dynamic|).
#       A horizontal reference line at the mean(|SPEI_PM|) for that month.
#       This answers "how much of the total SPEI deficit is each component?"
# ==============================================================================
cat("\n===== STEP 7: SEASONAL DECOMPOSITION 2022\u20132025 (SPEI-1/2/3) =====\n")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

# ── Ensure COL_* and theme_ms exist (guards for standalone sourcing) ────────
if (!exists("COL_PM"))      COL_PM      <- "#2166ac"   # blue  – SPEI_PM
if (!exists("COL_THW"))     COL_THW     <- "#e31a1c"   # red   – SPEI_Thw
if (!exists("COL_DYN"))     COL_DYN     <- "#ff7f00"   # orange– SPEI_Dynamic
if (!exists("out_dir"))     out_dir     <- "decomp_results"
if (!exists("decomp_all"))  decomp_all  <- read.csv(
  file.path(out_dir, "decomp_full_timeseries.csv"),
  stringsAsFactors = FALSE)

MONTH_ABBR <- c("Jan","Feb","Mar","Apr","May","Jun",
                "Jul","Aug","Sep","Oct","Nov","Dec")
SEAS_SCALES <- c(1L, 2L, 3L)
PANEL_LAB_7 <- c("(a) SPEI-1", "(b) SPEI-2", "(c) SPEI-3")

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
  
  # Reshape to long for a single geom_ribbon + geom_line call each
  # (avoids messy manual layers for 3 components)
  long_mean <- df %>%
    dplyr::select(month_abb, pm_mean, thw_mean, dyn_mean) %>%
    tidyr::pivot_longer(-month_abb,
                        names_to  = "component",
                        values_to = "mean_val") %>%
    dplyr::mutate(component = dplyr::recode(component,
                                            pm_mean  = "SPEI\u209a\u2098 (observed)",
                                            thw_mean = "SPEI\u2080 (thermodynamic)",
                                            dyn_mean = "SPEI\u2206 (dynamic)"))
  
  long_lo <- df %>%
    dplyr::select(month_abb, pm_lo, thw_lo, dyn_lo) %>%
    tidyr::pivot_longer(-month_abb, names_to = "component", values_to = "lo") %>%
    dplyr::mutate(component = dplyr::recode(component,
                                            pm_lo  = "SPEI\u209a\u2098 (observed)",
                                            thw_lo = "SPEI\u2080 (thermodynamic)",
                                            dyn_lo = "SPEI\u2206 (dynamic)"))
  
  long_hi <- df %>%
    dplyr::select(month_abb, pm_hi, thw_hi, dyn_hi) %>%
    tidyr::pivot_longer(-month_abb, names_to = "component", values_to = "hi") %>%
    dplyr::mutate(component = dplyr::recode(component,
                                            pm_hi  = "SPEI\u209a\u2098 (observed)",
                                            thw_hi = "SPEI\u2080 (thermodynamic)",
                                            dyn_hi = "SPEI\u2206 (dynamic)"))
  
  ribbon_df <- dplyr::left_join(long_mean, long_lo,  by = c("month_abb","component")) %>%
    dplyr::left_join(long_hi,              by = c("month_abb","component")) %>%
    dplyr::mutate(component = factor(component,
                                     levels = c("SPEI\u209a\u2098 (observed)",
                                                "SPEI\u2080 (thermodynamic)",
                                                "SPEI\u2206 (dynamic)")))
  
  col_pal <- c(
    "SPEI\u209a\u2098 (observed)"    = COL_PM,
    "SPEI\u2080 (thermodynamic)" = COL_THW,
    "SPEI\u2206 (dynamic)"      = COL_DYN)
  
  # y-axis range padded symmetrically
  ymax <- max(abs(c(df$pm_hi, df$thw_hi, df$dyn_hi)), na.rm = TRUE) * 1.15
  ymax <- max(ymax, 0.5, na.rm = TRUE)
  
  p <- ggplot2::ggplot(ribbon_df,
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
      title    = sprintf("%s  SPEI-%d  (2022\u20132025 monthly mean \u00b11\u202fSD)",
                         pan_lab, sc),
      subtitle = "Dashed = event threshold (\u22120.5).  \u00b11 SD shading from 4 calendar-year observations.",
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

fig5c <- p5c_1 / p5c_2 / p5c_3 +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Seasonal decomposition of 2022\u20132025 drought \u2014 Nechako River Basin",
    subtitle = paste0(
      "Monthly mean SPEI values (2022\u20132025) separated into observed (SPEI\u209a\u2098),",
      " thermodynamic (SPEI\u2080) and dynamic (SPEI\u2206 = SPEI\u209a\u2098 \u2212 SPEI\u2080) components.\n",
      "Shading = \u00b11 SD across the 4 calendar-year replicates (2022, 2023, 2024, 2025)."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35", hjust = 0,
                                            margin = ggplot2::margin(b = 5)),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 8),
      legend.key.size = ggplot2::unit(0.4, "cm")))

for (ext in c("pdf", "png")) {
  out_5c <- file.path(out_dir, paste0("Fig5c_seasonal_decomp_2022_2025.", ext))
  tryCatch(
    ggplot2::ggsave(out_5c, fig5c,
                    width  = 7.5, height = 9.5, units = "in",
                    dpi    = if (ext == "png") 300 else NULL,
                    device = ext),
    error = function(e) cat(sprintf("  \u26a0 Fig5c %s: %s\n", ext, e$message)))
  cat(sprintf("  \u2713 Saved: %s\n", basename(out_5c)))
}

# ══════════════════════════════════════════════════════════════════════════════
# Fig5d: stacked bar chart companion
#   For each month × scale: stack of |SPEI_Thw| (thermodynamic, red) and
#   |SPEI_Dynamic| (dynamic, orange), with a point + line for |SPEI_PM|.
#   Three panels (one per scale) in a single figure.
# ══════════════════════════════════════════════════════════════════════════════
cat("\n  Building Fig5d: stacked contribution bars...\n")

make_seas_bar_panel <- function(sc, pan_lab, show_xlab = FALSE,
                                show_leg  = FALSE) {
  df <- dplyr::filter(seas_stats, scale == sc)
  if (!nrow(df)) return(ggplot2::ggplot() + ggplot2::theme_void())
  
  # Long format for stacking
  bar_df <- df %>%
    dplyr::select(month_abb, thw_abs, dyn_abs, pm_abs) %>%
    tidyr::pivot_longer(cols = c(thw_abs, dyn_abs),
                        names_to  = "component",
                        values_to = "contribution") %>%
    dplyr::mutate(component = dplyr::recode(component,
                                            thw_abs = "Thermodynamic  |SPEI\u2080|",
                                            dyn_abs = "Dynamic  |SPEI\u2206|"))
  
  comp_levels <- c("Thermodynamic  |SPEI\u2080|", "Dynamic  |SPEI\u2206|")
  bar_df$component <- factor(bar_df$component, levels = comp_levels)
  
  col_bar <- c("Thermodynamic  |SPEI\u2080|" = COL_THW,
               "Dynamic  |SPEI\u2206|"       = COL_DYN)
  
  ymax <- max(bar_df$contribution[bar_df$component == "Thermodynamic  |SPEI\u2080|"] +
                bar_df$contribution[bar_df$component == "Dynamic  |SPEI\u2206|"],
              df$pm_abs, na.rm = TRUE) * 1.18
  ymax <- max(ymax, 0.3, na.rm = TRUE)
  
  p <- ggplot2::ggplot(bar_df,
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
      subtitle = "Bars = stacked |thermodynamic| + |dynamic|.  Line = |SPEI\u209a\u2098| (total deficit).",
      y        = "Mean absolute SPEI value  (2022\u20132025)") +
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

fig5d <- p5d_1 / p5d_2 / p5d_3 +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Seasonal absolute contribution of drought components \u2014 2022\u20132025",
    subtitle = paste0(
      "Stacked bars: mean absolute thermodynamic (|SPEI\u2080|, red) +",
      " dynamic (|SPEI\u2206|, orange) contributions per calendar month.\n",
      "Black line: total observed deficit (|SPEI\u209a\u2098|).",
      "  Values are means across 2022, 2023, 2024, 2025."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35", hjust = 0,
                                            margin = ggplot2::margin(b = 5)),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 8),
      legend.key.size = ggplot2::unit(0.4, "cm")))

for (ext in c("pdf", "png")) {
  out_5d <- file.path(out_dir, paste0("Fig5d_seasonal_decomp_stacked_2022_2025.", ext))
  tryCatch(
    ggplot2::ggsave(out_5d, fig5d,
                    width  = 7.5, height = 9.5, units = "in",
                    dpi    = if (ext == "png") 300 else NULL,
                    device = ext),
    error = function(e) cat(sprintf("  \u26a0 Fig5d %s: %s\n", ext, e$message)))
  cat(sprintf("  \u2713 Saved: %s\n", basename(out_5d)))
}

# ── Export seasonal stats CSV ─────────────────────────────────────────────────
write.csv(seas_stats,
          file.path(out_dir, "decomp_seasonal_stats_2022_2025.csv"),
          row.names = FALSE)
cat(sprintf("  \u2713 Saved: decomp_seasonal_stats_2022_2025.csv (%d rows)\n",
            nrow(seas_stats)))

# ── Console summary ───────────────────────────────────────────────────────────
cat("\n  2022-2025 seasonal decomp: peak months per component\n")
cat(sprintf("  %-8s  %-5s  %-12s  %-12s  %-12s\n",
            "Scale","Month","SPEI_PM","SPEI_Thw","SPEI_Dyn"))
cat("  ", paste(rep("-",56),collapse=""), "\n", sep="")
for (sc in SEAS_SCALES) {
  df_sc <- dplyr::filter(seas_stats, scale == sc)
  pk_pm  <- df_sc$month_abb[which.min(df_sc$pm_mean)]
  pk_thw <- df_sc$month_abb[which.min(df_sc$thw_mean)]
  pk_dyn <- df_sc$month_abb[which.min(df_sc$dyn_mean)]
  cat(sprintf("  SPEI-%-3d  %-5s  %-12s  %-12s  %-12s\n",
              sc, "(peak)",
              sprintf("%s (%.2f)", pk_pm,  min(df_sc$pm_mean,  na.rm=TRUE)),
              sprintf("%s (%.2f)", pk_thw, min(df_sc$thw_mean, na.rm=TRUE)),
              sprintf("%s (%.2f)", pk_dyn, min(df_sc$dyn_mean, na.rm=TRUE))))
}
####################################################################################
# w11_EXPORT_ADDON.R
#
# PURPOSE: Export the annual JJA thermodynamic fraction time series in a compact
#          format for use as the physical covariate in w12_bayesian_ene_analysis.R
#          (MS2). Append this block to the END of w11_dynamic_thermodynamic_decomp.R,
#          after all existing steps (1–6) have run.
#
# REQUIRES (already in environment from w11 main body):
#   decomp_all   — data.frame: date, year, month, scale, SPEI_PM, SPEI_Thw,
#                  SPEI_Dynamic, thm_frac, dyn_frac
#   out_dir      — character, decomp results directory
#
# OUTPUTS written to out_dir:
#   thm_frac_annual_jja_ms2.csv
#     One row per year × scale combination (scales 1, 2, 3).
#     Columns: year, scale, thm_frac_jja, n_jja_months, SPEI_PM_jja,
#              SPEI_Thw_jja, SPEI_Dynamic_jja
####################################################################################
# ==============================================================================
# Following Produces 3 CSV files in decomp_results/ with all remaining manuscript values
# ==============================================================================
seasonal_fracs <- decomp_all %>%
  filter(year >= 2022, year <= 2025, scale == 3) %>%
  mutate(season = case_when(
    month %in% c(6, 7, 8)   ~ "JJA",
    month %in% c(12, 1, 2)  ~ "DJF",
    month %in% c(3, 4, 5)   ~ "MAM",
    month %in% c(9, 10, 11) ~ "SON"
  )) %>%
  group_by(season) %>%
  summarise(
    thm_frac_pct = 100 * mean(abs(SPEI_Thw) / (abs(SPEI_PM) + 1e-6), na.rm = TRUE),
    .groups = "drop"
  )
print(seasonal_fracs)
write.csv(seasonal_fracs,
          file.path(out_dir, "ms1_seasonal_fractions.csv"), row.names = FALSE)
cat("  -> ms1_seasonal_fractions.csv\n")

# ── 2. Offset months and mean SPEI_Δ (§4.4, Table 3 last column) ───────────
offset_months <- decomp_all %>%
  filter(year >= 2022, year <= 2025) %>%
  group_by(scale) %>%
  summarise(
    n_total      = n(),
    n_positive   = sum(SPEI_Dynamic > 0, na.rm = TRUE),
    pct_positive = 100 * mean(SPEI_Dynamic > 0, na.rm = TRUE),
    mean_dynamic = mean(SPEI_Dynamic, na.rm = TRUE),
    .groups = "drop"
  )
# Which event years had the most positive SPEI_Δ months?
offset_by_year <- decomp_all %>%
  filter(year >= 2022, year <= 2025, scale == 3) %>%
  group_by(year) %>%
  summarise(n_pos = sum(SPEI_Dynamic > 0, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(n_pos))
print(offset_months)
print(offset_by_year)
write.csv(offset_months,
          file.path(out_dir, "ms1_offset_months.csv"), row.names = FALSE)
write.csv(offset_by_year,
          file.path(out_dir, "ms1_offset_by_year.csv"), row.names = FALSE)
cat("  -> ms1_offset_months.csv\n  -> ms1_offset_by_year.csv\n")

# ── 3. JJA thermodynamic fraction trend — SPEI-3 (§4.6, §5.3, Conclusion 3)
# NOTE: Current STEP 5 uses SPEI-6. The manuscript §4.6 uses SPEI-3.
trend3_df <- decomp_all %>%
  filter(scale == 3, month %in% c(6, 7, 8)) %>%
  group_by(year) %>%
  summarise(
    thm_frac_abs = mean(abs(SPEI_Thw) / (abs(SPEI_PM) + 1e-6), na.rm = TRUE),
    .groups = "drop"
  )

# Full-record OLS (1950-2025)
fit_full   <- lm(thm_frac_abs ~ year, data = trend3_df)
sl_full    <- coef(fit_full)["year"]
pv_full    <- summary(fit_full)$coefficients["year", 4]
tstat_full <- summary(fit_full)$coefficients["year", 3]
ci95       <- confint(fit_full, "year", level = 0.95)

# Post-1990 OLS
trend3_p90 <- filter(trend3_df, year >= 1990)
fit_p90    <- lm(thm_frac_abs ~ year, data = trend3_p90)
sl_p90     <- coef(fit_p90)["year"]
pv_p90     <- summary(fit_p90)$coefficients["year", 4]

# Half-period means
early <- trend3_df %>% filter(year <= 1987)
late  <- trend3_df %>% filter(year >= 1988)
mn_early <- mean(early$thm_frac_abs, na.rm = TRUE) * 100
mn_late  <- mean(late$thm_frac_abs,  na.rm = TRUE) * 100

# Mann-Whitney
mw <- wilcox.test(early$thm_frac_abs, late$thm_frac_abs)

# ENSO-removed: remove known El Niño years (ONI > +0.5 for DJF)
# Standard list of El Nino years (peak year): 1958,1966,1969,1973,1977,1983,
# 1987,1988,1992,1995,1998,2003,2005,2007,2010,2015,2016,2019,2023,2024
enso_years <- c(1958,1966,1969,1973,1977,1983,1987,1988,1992,1995,
                1998,2003,2005,2007,2010,2015,2016,2019,2023,2024)
trend3_noenso <- filter(trend3_df, !year %in% enso_years)
fit_noenso <- lm(thm_frac_abs ~ year, data = trend3_noenso)
sl_noenso  <- coef(fit_noenso)["year"]
pv_noenso  <- summary(fit_noenso)$coefficients["year", 4]

ms1_trend <- data.frame(
  metric = c("slope_per_yr","slope_per_decade","p_value","t_stat",
             "ci95_lo_per_yr","ci95_hi_per_yr",
             "mean_1950_1987_pct","mean_1988_2025_pct","change_pp",
             "mann_whitney_p",
             "slope_noenso_per_decade","p_noenso",
             "slope_p90_per_decade","p_p90"),
  value  = c(sl_full, sl_full*10, pv_full, tstat_full,
             ci95[1], ci95[2],
             mn_early, mn_late, mn_late - mn_early,
             mw$p.value,
             sl_noenso*10, pv_noenso,
             sl_p90*10, pv_p90)
)
print(ms1_trend)
write.csv(ms1_trend,
          file.path(out_dir, "ms1_trend_spei3_jja.csv"), row.names = FALSE)
cat("  -> ms1_trend_spei3_jja.csv\n")

cat("\n===== DONE: 4 CSV files written to decomp_results/ =====\n")
cat("  ms1_seasonal_fractions.csv\n")
cat("  ms1_offset_months.csv\n")
cat("  ms1_offset_by_year.csv\n")
cat("  ms1_trend_spei3_jja.csv\n")
cat("\n══════════════════════════════════════════════════\n")
cat("  EXPORT ADDON: MS2 covariate export (JJA thm_frac)\n")
cat("══════════════════════════════════════════════════\n")

if (!exists("decomp_all") || nrow(decomp_all) == 0) {
  cat("  ⚠ decomp_all not found or empty — skipping MS2 covariate export.\n")
} else {
  
  MS2_SCALES <- c(1L, 2L, 3L)   # scales used in MS1 and carried into MS2
  JJA_MONTHS <- c(6L, 7L, 8L)   # June-July-August
  
  thm_jja <- decomp_all %>%
    dplyr::filter(
      scale %in% MS2_SCALES,
      month %in% JJA_MONTHS,
      is.finite(thm_frac),
      is.finite(SPEI_PM)
    ) %>%
    dplyr::group_by(year, scale) %>%
    dplyr::summarise(
      # Primary covariate: mean JJA thermodynamic fraction
      thm_frac_jja     = mean(abs(SPEI_Thw) / (abs(SPEI_PM) + 1e-6), na.rm = TRUE),
      # Supporting diagnostics
      SPEI_PM_jja      = mean(SPEI_PM,      na.rm = TRUE),
      SPEI_Thw_jja     = mean(SPEI_Thw,     na.rm = TRUE),
      SPEI_Dynamic_jja = mean(SPEI_Dynamic, na.rm = TRUE),
      n_jja_months     = dplyr::n(),
      .groups          = "drop"
    ) %>%
    # Add a standardised version for model diagnostics (mean-centred, unit SD)
    dplyr::group_by(scale) %>%
    dplyr::mutate(
      thm_frac_jja_std = as.numeric(scale(thm_frac_jja))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(scale, year)
  
  ms2_cov_path <- file.path(out_dir, "thm_frac_annual_jja_ms2.csv")
  write.csv(thm_jja, ms2_cov_path, row.names = FALSE)
  cat(sprintf("  ✓ JJA thm_frac covariate: %s\n", basename(ms2_cov_path)))
  cat(sprintf("    Scales exported: %s\n", paste(MS2_SCALES, collapse = ", ")))
  cat(sprintf("    Years: %d–%d  (%d rows total)\n",
              min(thm_jja$year), max(thm_jja$year), nrow(thm_jja)))
  
  # Summarise trend in the covariate by scale (quick diagnostic)
  cat("\n  JJA thm_frac trend by scale (OLS, printed for verification):\n")
  for (sc in MS2_SCALES) {
    df_sc <- dplyr::filter(thm_jja, scale == sc, !is.na(thm_frac_jja))
    if (nrow(df_sc) < 10) next
    lm_sc <- lm(thm_frac_jja ~ year, data = df_sc)
    slope  <- coef(lm_sc)["year"]
    pval   <- summary(lm_sc)$coefficients["year", 4]
    cat(sprintf("    SPEI-%d: slope = %+.5f/yr  (%+.4f/decade),  p = %.4f\n",
                sc, slope, slope * 10, pval))
  }
}

cat("══════════════════════════════════════════════════\n")
cat("  MS2 covariate export complete.\n")
cat("══════════════════════════════════════════════════\n")
# ==============================================================================
# COMPLETE
# ==============================================================================
cat("\n============================================================\n")
cat("DYNAMIC/THERMODYNAMIC DECOMPOSITION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("Outputs saved to: %s\n", normalizePath(out_dir)))
cat("\nSUPPLEMENTARY FIGURES:\n")
cat("  fig_decomp_spei6_timeseries.pdf  \u2014 SPEI-6 full record (1950\u20132025)\n")
cat("  fig_annual_decomp_barplot.pdf    \u2014 annual SPEI-6 drought-year bars\n")
cat("  fig_thermodynamic_trend.pdf      \u2014 JJA SPEI-6 thermodynamic trend\n")
cat("  fig_scale_comparison.pdf         \u2014 all scales 1\u201324 percentage stack\n")
cat("\nMANUSCRIPT FIGURES:\n")
cat("  Fig5a_decomp_timeseries_2020_2025.pdf/.png\n")
cat("    \u2192 SPEIpm / SPEI0 / SPEIdyn, scales 1\u20133, 2020\u20132025 (Section 4.3)\n")
cat("  Fig5b_scale_comparison_stacked.pdf/.png\n")
cat("    \u2192 Stacked % bars, scales 1\u20133, thermodynamic vs dynamic (Table 3 / Section 4.2)\n")
cat("  Fig6_thermodynamic_trend.pdf/.png\n")
cat("    \u2192 JJA F_thm 1950\u20132025: full OLS \u00b1 95% CI + post-1990 OLS \u00b1 90% CI (Section 4.4)\n")
cat("\nINTERPRETATION:\n")
cat("  SPEI_PM       = Full observed drought (Penman-Monteith PET)\n")
cat("  SPEI_Thw      = Thermodynamic component (Thornthwaite PET, temperature-only)\n")
cat("  SPEI_Dynamic  = SPEI_PM \u2212 SPEI_Thw = dynamic component\n")
cat("                \u2248 precipitation + radiation + wind effects\n")