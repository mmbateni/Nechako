# ============================================================================== 
# 10b_PET_bias_nonstationarity_v2.R
# ==============================================================================
# SUPPLEMENTARY SCRIPT — quantifies raw PET_PM - PET_Thw = PET_PM - PET_Thw (Methods 3b),
# NOT SPEI-based. Unlike 10a, this script does not decompose drought severity
# into thermodynamic/dynamic components, so it does not conflict with the
# current manuscript's PET-formulation-sensitivity framing (Methods 3e) and
# needs no reinterpretation. It is a useful companion to
# 12PET_formulation_sensitivity_v2.R for describing the PET input itself
# rather than the resulting SPEI.
#
# This script is independent of 10a. It reads the two authoritative monthly
# PET summary products directly and therefore does not require the legacy
# SPEI/decomposition output produced by 10a.
#
# THORNTHWAITE vs. PENMAN-MONTEITH PET BIAS NON-STATIONARITY ANALYSIS
# Nechako River Basin, BC — observed-framework support script
#
# PURPOSE:
#   Quantify the observed monthly/seasonal non-stationarity of PET disagreement
#   between Thornthwaite and Penman–Monteith, without any detrending or PET
#   bias-correction/counterfactual reconstruction.
#
# INPUTS:
#   monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv
#   monthly_data_direct/ERA5Land_Nechako_PET_Thornthwaite_monthly_summary.csv
#
# OUTPUTS -> decomp_results/bias_nonstationarity/
#   bias_monthly_timeseries.csv
#   bias_annual_summary.csv
#   bias_seasonal_summary.csv
#   bias_trend_tests.csv
#   bias_nonstationarity.pdf
#   bias_manuscript_tables.xlsx
# ============================================================================== 

library(zoo)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(sandwich)
# sandwich is used for HAC/autocorrelation-robust annual trend inference.

# ── Paths ─────────────────────────────────────────────────────────────────────
setwd("D:/Nechako_Drought/Nechako/")
PET_DIR    <- "monthly_data_direct"
DECOMP_DIR <- "decomp_results"
OUT_DIR    <- file.path(DECOMP_DIR, "bias_nonstationarity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Settings ──────────────────────────────────────────────────────────────────
P1_START <- 1950L; P1_END <- 1990L
P3_START <- 2022L; P3_END <- 2025L

season_of <- function(m) dplyr::case_when(
  m %in% c(12L,1L,2L)  ~ "DJF",
  m %in% c(3L,4L,5L)   ~ "MAM",
  m %in% c(6L,7L,8L)   ~ "JJA",
  m %in% c(9L,10L,11L) ~ "SON"
)

cat("==============================================================\n")
cat("  w10b: PET BIAS NON-STATIONARITY ANALYSIS (observed only)\n")
cat("==============================================================\n\n")

# ============================================================================
# SECTION 1: LOAD INPUTS
# ============================================================================
cat("--- SECTION 1: Loading inputs ---\n")

pm_file  <- file.path(PET_DIR, "ERA5Land_Nechako_PET_monthly_summary.csv")
thw_file <- file.path(PET_DIR, "ERA5Land_Nechako_PET_Thornthwaite_monthly_summary.csv")

for (f in c(pm_file, thw_file)) {
  if (!file.exists(f)) stop("File not found — check PET prerequisites: ", f)
}

pm_raw  <- read.csv(pm_file, stringsAsFactors = FALSE)
thw_raw <- read.csv(thw_file, stringsAsFactors = FALSE)

cat(sprintf("  PM PET:  %d rows (%s to %s)\n",
            nrow(pm_raw), pm_raw$date[1], tail(pm_raw$date, 1)))
cat(sprintf("  Thw PET: %d rows (%s to %s)\n",
            nrow(thw_raw), thw_raw$date[1], tail(thw_raw$date, 1)))

normalise_month_date <- function(x, label = "date") {
  x <- trimws(as.character(x))
  
  if (!length(x)) {
    stop(label, " is empty.")
  }
  
  # The PET summary files use YYYY-MM. Do not pass YYYY-MM directly to
  # as.Date(): base R can call charToDate() and error before returning NA.
  # Normalize recognized monthly strings to YYYY-MM-01 first.
  is_ym <- grepl("^\\d{4}-\\d{2}$", x)
  is_ymd <- grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  
  unsupported <- !(is_ym | is_ymd)
  
  if (any(unsupported)) {
    bad_values <- unique(x[unsupported])
    stop(
      label,
      " contains unsupported date format(s): ",
      paste(head(bad_values, 10L), collapse = ", "),
      ". Expected YYYY-MM or YYYY-MM-DD."
    )
  }
  
  x_canonical <- ifelse(
    is_ym,
    paste0(x, "-01"),
    x
  )
  
  out <- as.Date(
    x_canonical,
    format = "%Y-%m-%d"
  )
  
  if (anyNA(out)) {
    bad_values <- unique(x[is.na(out)])
    stop(
      label,
      " contains invalid calendar date(s): ",
      paste(head(bad_values, 10L), collapse = ", ")
    )
  }
  
  # Force monthly canonical representation even if YYYY-MM-DD was supplied.
  as.Date(
    format(out, "%Y-%m-01"),
    format = "%Y-%m-%d"
  )
}

for (nm in c("date", "mean_pet")) {
  if (!nm %in% names(pm_raw)) stop("PM PET summary missing required column: ", nm)
  if (!nm %in% names(thw_raw)) stop("Thornthwaite PET summary missing required column: ", nm)
}

pm_raw$date  <- normalise_month_date(pm_raw$date, "PM PET date")
thw_raw$date <- normalise_month_date(thw_raw$date, "Thornthwaite PET date")

cat(sprintf(
  "  Canonical dates: %s to %s\n",
  format(min(pm_raw$date), "%Y-%m"),
  format(max(pm_raw$date), "%Y-%m")
))

if (!identical(pm_raw$date, thw_raw$date)) {
  stop(
    "PM and Thornthwaite PET monthly date sequences are not identical."
  )
}
pm_raw$year  <- as.integer(format(pm_raw$date, "%Y"))
pm_raw$month <- as.integer(format(pm_raw$date, "%m"))
thw_raw$year  <- as.integer(format(thw_raw$date, "%Y"))
thw_raw$month <- as.integer(format(thw_raw$date, "%m"))

expected_dates <- seq(as.Date("1950-01-01"), as.Date("2025-12-01"), by = "month")

validate_pet_record <- function(x, label) {
  if (length(x) != length(expected_dates) || !identical(x, expected_dates))
    stop(label, " does not contain the complete chronological 1950-01 to 2025-12 monthly record.")
  if (anyDuplicated(x)) stop(label, " contains duplicated months.")
  invisible(TRUE)
}

validate_pet_record(pm_raw$date, "PM PET summary")
validate_pet_record(thw_raw$date, "Thornthwaite PET summary")

cat(sprintf("  Canonical PET record: %s to %s (%d months)\n",
            format(min(expected_dates), "%Y-%m"),
            format(max(expected_dates), "%Y-%m"),
            length(expected_dates)))

# ============================================================================
# SECTION 2: BUILD MONTHLY BASE TABLE
# ============================================================================
cat("\n--- SECTION 2: Building monthly base table ---\n")

base <- pm_raw %>%
  select(date, year, month, mean_pet) %>%
  rename(pet_pm_mm_day = mean_pet) %>%
  inner_join(
    thw_raw %>%
      select(date, year, month, mean_pet) %>%
      rename(pet_thw_mm_day = mean_pet),
    by = c("date", "year", "month")
  ) %>%
  arrange(date) %>%
  mutate(
    days_in_month = lubridate::days_in_month(date),
    delta_pet_mm_day = pet_pm_mm_day - pet_thw_mm_day,
    delta_pet_mm_month = delta_pet_mm_day * days_in_month,
    season = season_of(month),
    in_P1 = year >= P1_START & year <= P1_END,
    in_P3 = year >= P3_START & year <= P3_END
  )

if (nrow(base) != length(expected_dates) || !identical(base$date, expected_dates))
  stop("Joined PET record is not exactly the 1950-01 to 2025-12 monthly record.")

if (any(!is.finite(base$pet_pm_mm_day)) || any(!is.finite(base$pet_thw_mm_day)))
  stop("PET summaries contain non-finite values in the joined record.")

cat(sprintf("  Merged: %d months\n", nrow(base)))
cat(sprintf("  PM PET (mm/day):  %.4f to %.4f\n", min(base$pet_pm_mm_day), max(base$pet_pm_mm_day)))
cat(sprintf("  Thw PET (mm/day): %.4f to %.4f\n", min(base$pet_thw_mm_day), max(base$pet_thw_mm_day)))
cat(sprintf("  PET_PM - PET_Thw (mm/day): %.4f to %.4f\n", min(base$delta_pet_mm_day), max(base$delta_pet_mm_day)))
cat(sprintf("  PET_PM - PET_Thw (mm/month): %.3f to %.3f\n", min(base$delta_pet_mm_month), max(base$delta_pet_mm_month)))

write.csv(
  base %>% select(date, year, month, season, days_in_month,
                  pet_pm_mm_day, pet_thw_mm_day,
                  delta_pet_mm_day, delta_pet_mm_month, in_P1, in_P3),
  file.path(OUT_DIR, "bias_monthly_timeseries.csv"), row.names = FALSE
)
cat("  Saved: bias_monthly_timeseries.csv\n")

# ============================================================================
# SECTION 3: ANNUAL AND SEASONAL SUMMARIES
# ============================================================================
cat("\n--- SECTION 3: Annual and seasonal summaries ---\n")

annual_bias <- base %>%
  group_by(year) %>%
  summarise(
    mean_delta_pet_mm_day   = mean(delta_pet_mm_day, na.rm = TRUE),
    mean_delta_pet_mm_month = mean(delta_pet_mm_month, na.rm = TRUE),
    .groups = "drop"
  )

seasonal_bias <- base %>%
  group_by(season) %>%
  summarise(
    mean_delta_pet_mm_day   = mean(delta_pet_mm_day, na.rm = TRUE),
    mean_delta_pet_mm_month = mean(delta_pet_mm_month, na.rm = TRUE),
    sd_delta_pet_mm_day     = sd(delta_pet_mm_day, na.rm = TRUE),
    n                  = n(),
    .groups = "drop"
  ) %>%
  mutate(season = factor(season, levels = c("DJF", "MAM", "JJA", "SON")))

p1_annual <- annual_bias %>% filter(year >= P1_START, year <= P1_END)
p3_annual <- annual_bias %>% filter(year >= P3_START, year <= P3_END)

annual_summary <- data.frame(
  period = c("P1", "P3"),
  start_year = c(P1_START, P3_START),
  end_year   = c(P1_END, P3_END),
  mean_delta_pet_mm_day   = c(mean(p1_annual$mean_delta_pet_mm_day, na.rm = TRUE), mean(p3_annual$mean_delta_pet_mm_day, na.rm = TRUE)),
  mean_delta_pet_mm_month = c(mean(p1_annual$mean_delta_pet_mm_month, na.rm = TRUE), mean(p3_annual$mean_delta_pet_mm_month, na.rm = TRUE))
)

# ---------------------------------------------------------------------------
# Trend inference
#
# Monthly PET differences are serially correlated. Therefore, the conventional
# OLS standard error/p-value from:
#
#   lm(delta_pet_mm_day ~ year + month)
#
# is retained only as a descriptive reference and is NOT used for manuscript
# inference.
#
# Primary inference is based on annual mean DeltaPET, one observation per year,
# with Newey-West HAC standard errors. The annual series is still a time series,
# so HAC inference is preferred to naive OLS.
#
# The monthly model is retained as a descriptive sensitivity check only.
# ---------------------------------------------------------------------------

trend_monthly <- lm(
  delta_pet_mm_day ~ year + month,
  data = base
)

trend_annual <- lm(
  mean_delta_pet_mm_day ~ year,
  data = annual_bias
)

# Newey-West HAC standard error for the annual trend.
if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop(
    "Package 'sandwich' is required for autocorrelation-robust annual ",
    "PET-bias trend inference. Install it with install.packages('sandwich')."
  )
}

annual_vcov_hac <- sandwich::NeweyWest(
  trend_annual,
  lag = 4L,
  prewhite = FALSE,
  adjust = TRUE
)

annual_beta <- unname(
  coef(trend_annual)["year"]
)

annual_se_hac <- sqrt(
  annual_vcov_hac["year", "year"]
)

annual_t_hac <- annual_beta / annual_se_hac

annual_df <- nrow(annual_bias) - length(coef(trend_annual))

annual_p_hac <- 2 * stats::pt(
  abs(annual_t_hac),
  df = annual_df,
  lower.tail = FALSE
)

annual_ci_hac <- annual_beta +
  c(-1, 1) *
  stats::qt(
    0.975,
    df = annual_df
  ) *
  annual_se_hac

# Descriptive monthly OLS quantities. These are explicitly labelled as
# non-HAC/reference results and are not the manuscript inference.
monthly_beta <- unname(
  coef(trend_monthly)["year"]
)

monthly_p_ols <- summary(
  trend_monthly
)$coefficients["year", "Pr(>|t|)"]

# Exact p-value formatting prevents a statistically meaningful small p-value
# from being printed misleadingly as "0.0000".
format_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 1e-4) return(formatC(p, format = "scientific", digits = 3))
  formatC(p, format = "f", digits = 6)
}

trend_df <- data.frame(
  metric = c(
    "monthly_ols_reference_slope_mm_day_per_year",
    "monthly_ols_reference_p",
    "annual_ols_slope_mm_day_per_year",
    "annual_hac_se_mm_day_per_year",
    "annual_hac_t",
    "annual_hac_df",
    "annual_hac_p",
    "annual_hac_ci95_low_mm_day_per_year",
    "annual_hac_ci95_high_mm_day_per_year"
  ),
  value = c(
    monthly_beta,
    monthly_p_ols,
    annual_beta,
    annual_se_hac,
    annual_t_hac,
    annual_df,
    annual_p_hac,
    annual_ci_hac[1],
    annual_ci_hac[2]
  ),
  inference_role = c(
    "descriptive reference only",
    "descriptive reference only",
    "effect estimate",
    "HAC standard error",
    "HAC test statistic",
    "HAC degrees of freedom",
    "primary manuscript inference",
    "primary manuscript inference",
    "primary manuscript inference"
  ),
  stringsAsFactors = FALSE
)

write.csv(annual_bias,   file.path(OUT_DIR, "bias_annual_summary.csv"),   row.names = FALSE)
write.csv(seasonal_bias, file.path(OUT_DIR, "bias_seasonal_summary.csv"), row.names = FALSE)
write.csv(trend_df,      file.path(OUT_DIR, "bias_trend_tests.csv"),     row.names = FALSE)
cat("  Saved: bias_annual_summary.csv, bias_seasonal_summary.csv, bias_trend_tests.csv\n")

# ============================================================================
# SECTION 4: FIGURES
# ============================================================================
cat("\n--- SECTION 4: Figures ---\n")

pdf(file.path(OUT_DIR, "bias_nonstationarity.pdf"), width = 14, height = 9)

base_plt <- base %>%
  mutate(date_fmt = as.Date(date),
         delta_pet_12ma = stats::filter(delta_pet_mm_month, rep(1/12, 12), sides = 2) %>% as.numeric())
p1_m <- mean(base$delta_pet_mm_month[base$in_P1], na.rm = TRUE)
p3_m <- mean(base$delta_pet_mm_month[base$in_P3], na.rm = TRUE)

f1 <- ggplot(base_plt, aes(x = date_fmt, y = delta_pet_mm_month)) +
  annotate("rect", xmin = as.Date("1950-01-01"), xmax = as.Date("1990-12-31"),
           ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "#4393c3") +
  annotate("rect", xmin = as.Date("2022-01-01"), xmax = as.Date("2025-12-31"),
           ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "#762a83") +
  geom_line(colour = "grey70", linewidth = 0.35, alpha = 0.8) +
  geom_line(aes(y = delta_pet_12ma), colour = "#e31a1c", linewidth = 1.0, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_hline(yintercept = p1_m, colour = "#4393c3", linewidth = 0.9) +
  geom_hline(yintercept = p3_m, colour = "#762a83", linewidth = 0.9) +
  annotate("text", x = as.Date("1970-01-01"), y = p1_m + 0.5,
           label = sprintf("P1 mean = %+.3f mm/month", p1_m),
           colour = "#4393c3", size = 3.0) +
  annotate("text", x = as.Date("2023-07-01"), y = p3_m + 0.5,
           label = sprintf("P3 mean = %+.3f mm/month", p3_m),
           colour = "#762a83", size = 3.0) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  labs(title = "Monthly bias: PET_PM - PET_Thw [mm/month]",
       subtitle = "Blue = P1 (1950-1990); Purple = P3 (2022-2025); Red = 12-month moving average",
       x = "Date", y = "PET_PM - PET_Thw (mm/month)") +
  theme_bw(base_size = 11)
print(f1)

f2 <- ggplot(annual_bias, aes(x = year, y = mean_delta_pet_mm_day)) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.8, colour = "#4a5568") +
  geom_smooth(method = "lm", se = TRUE, colour = "#e31a1c", fill = "#e31a1c", alpha = 0.15) +
  annotate("rect", xmin = P1_START, xmax = P1_END, ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "#4393c3") +
  annotate("rect", xmin = P3_START, xmax = P3_END, ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "#762a83") +
  labs(title = "Annual mean PET bias non-stationarity",
       subtitle = sprintf("OLS slope = %+.4f mm/day/year, p = %.4f", coef(trend_annual)["year"], summary(trend_annual)$coefficients["year", 4]),
       x = "Year", y = "Mean PET_PM - PET_Thw (mm/day)") +
  theme_bw(base_size = 11)
print(f2)

f3 <- ggplot(seasonal_bias, aes(x = season, y = mean_delta_pet_mm_day, fill = season)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  scale_fill_manual(values = c(DJF = "#4393c3", MAM = "#74c476", JJA = "#fd8d3c", SON = "#9e9ac8"), guide = "none") +
  labs(title = "Seasonal mean PET bias (Penman-Monteith - Thornthwaite)",
       x = "Season", y = "Mean PET_PM - PET_Thw (mm/day)") +
  theme_bw(base_size = 11)
print(f3)

dev.off()
cat(sprintf("  Saved: %s\n", file.path(OUT_DIR, "bias_nonstationarity.pdf")))

# ============================================================================
# SECTION 5: EXCEL WORKBOOK
# ============================================================================
cat("\n--- SECTION 5: Excel workbook ---\n")
write.xlsx(
  list(
    Bias_monthly_timeseries = base %>%
      select(date, year, month, season, days_in_month,
             pet_pm_mm_day, pet_thw_mm_day,
             delta_pet_mm_day, delta_pet_mm_month, in_P1, in_P3),
    Bias_annual_summary     = annual_bias,
    Bias_seasonal_summary   = seasonal_bias,
    Bias_trend_tests        = trend_df
  ),
  file = file.path(OUT_DIR, "bias_manuscript_tables.xlsx"), overwrite = TRUE
)
cat("  Saved: bias_manuscript_tables.xlsx\n")

# ============================================================================
# SECTION 6: MANUSCRIPT SUMMARY
# ============================================================================
cat("\n==============================================================\n")
cat("  MANUSCRIPT REPORTING SUMMARY\n")
cat("==============================================================\n")
cat(sprintf("\n  Mean PET bias (PET_PM - PET_Thw):\n    P1 annual mean (1950-1990):  %+.3f mm/day\n    P3 annual mean (2022-2025):  %+.3f mm/day\n    Annual delta bias:           %+.3f mm/day\n    Annual trend slope (HAC inference):            %+.4f mm/day/year (p = %.4f)\n\n",
            mean(p1_annual$mean_delta_pet_mm_day, na.rm = TRUE),
            mean(p3_annual$mean_delta_pet_mm_day, na.rm = TRUE),
            mean(p3_annual$mean_delta_pet_mm_day, na.rm = TRUE) - mean(p1_annual$mean_delta_pet_mm_day, na.rm = TRUE),
            coef(trend_annual)["year"], summary(trend_annual)$coefficients["year", 4]))

cat(sprintf("  Seasonal mean bias (mm/day): DJF=%+.3f, MAM=%+.3f, JJA=%+.3f, SON=%+.3f\n",
            seasonal_bias$mean_delta_pet_mm_day[seasonal_bias$season == "DJF"],
            seasonal_bias$mean_delta_pet_mm_day[seasonal_bias$season == "MAM"],
            seasonal_bias$mean_delta_pet_mm_day[seasonal_bias$season == "JJA"],
            seasonal_bias$mean_delta_pet_mm_day[seasonal_bias$season == "SON"]))

cat("==============================================================\n")
cat("  Trend inference: annual Newey-West HAC; monthly OLS is descriptive only.\n")
cat("  Plot labels: ASCII-only for PDF device compatibility.\n")
cat("  w10b COMPLETE. Outputs: decomp_results/bias_nonstationarity/\n")
cat("==============================================================\n")