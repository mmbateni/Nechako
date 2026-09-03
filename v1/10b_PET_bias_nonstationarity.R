# ============================================================================
# 10b_PET_bias_nonstationarity_v2.R
# OBSERVED PET-SOURCE DIFFERENCE / NON-STATIONARITY DIAGNOSTIC
#
# Purpose
#   Quantify the difference between independently calculated FAO-56
#   Penman-Monteith reference evapotranspiration and the native ERA5-Land
#   potential evaporation (PEV) product.
#
# Scientific terminology
#   PET_PM   = FAO-56 Penman-Monteith reference evapotranspiration.
#   PET_ERA5 = ERA5-Land native potential evaporation (PEV).
#   delta_pet = PET_PM - PET_ERA5.
#
# IMPORTANT:
#   ERA5-Land PEV is not relabelled as Penman-Monteith ET. ECMWF documents the
#   ERA5-Land PEV product as an open-water/pan-evaporation quantity. Therefore
#   this script diagnoses PET-source/product sensitivity, not agreement between
#   two PM implementations.
#
# Inputs
#   monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv
#   monthly_data_direct/ERA5Land_Nechako_PET_ERA5Land_monthly_summary.csv
#
# Outputs -> decomp_results/bias_nonstationarity/
#   bias_monthly_timeseries.csv
#   bias_annual_summary.csv
#   bias_seasonal_summary.csv
#   bias_trend_tests.csv
#   bias_nonstationarity.pdf
#   bias_manuscript_tables.xlsx
# ============================================================================

library(zoo)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(sandwich)
library(lmtest)

WD_PATH <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

PET_DIR <- file.path(WD_PATH, "monthly_data_direct")
OUT_DIR <- file.path(WD_PATH, "decomp_results", "bias_nonstationarity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

P1_START <- 1950L
P1_END   <- 1990L
P3_START <- 2022L
P3_END   <- 2025L

season_of <- function(m) {
  dplyr::case_when(
    m %in% c(12L, 1L, 2L) ~ "DJF",
    m %in% c(3L, 4L, 5L)  ~ "MAM",
    m %in% c(6L, 7L, 8L)  ~ "JJA",
    m %in% c(9L, 10L, 11L) ~ "SON",
    TRUE ~ NA_character_
  )
}

normalise_month_date <- function(x, label) {
  x <- as.character(x)
  out <- as.Date(ifelse(nchar(x) == 7L, paste0(x, "-01"), x))
  if (anyNA(out)) stop(label, " contains invalid dates.")
  as.Date(format(out, "%Y-%m-01"))
}

read_pet_summary <- function(path, label) {
  if (!file.exists(path)) stop("Required input not found: ", path)
  x <- data.table::fread(path, data.table = FALSE)
  required <- c("date", "year", "month", "mean_pet")
  miss <- setdiff(required, names(x))
  if (length(miss)) stop(label, " is missing columns: ", paste(miss, collapse = ", "))
  x$date <- normalise_month_date(x$date, paste0(label, " date"))
  x$year <- as.integer(format(x$date, "%Y"))
  x$month <- as.integer(format(x$date, "%m"))
  x$mean_pet <- as.numeric(x$mean_pet)
  x <- x[order(x$date), c("date", "year", "month", "mean_pet")]
  if (anyDuplicated(x$date)) stop(label, " contains duplicated monthly dates.")
  x
}

cat("==============================================================\n")
cat("  PET SOURCE DIFFERENCE / NON-STATIONARITY ANALYSIS\n")
cat("  FAO-56 PM vs. ERA5-Land native PEV\n")
cat("==============================================================\n\n")

pm <- read_pet_summary(
  file.path(PET_DIR, "ERA5Land_Nechako_PET_monthly_summary.csv"),
  "FAO-56 PM PET summary"
)
era5 <- read_pet_summary(
  file.path(PET_DIR, "ERA5Land_Nechako_PET_ERA5Land_monthly_summary.csv"),
  "ERA5-Land PEV summary"
)

if (!identical(pm$date, era5$date)) stop("PM and ERA5-Land PEV monthly date sequences are not identical.")

base <- dplyr::left_join(
  pm %>% dplyr::select(date, year, month, pet_pm_mm_day = mean_pet),
  era5 %>% dplyr::select(date, pet_era5_mm_day = mean_pet),
  by = "date"
) %>%
  dplyr::mutate(
    days_in_month = lubridate::days_in_month(date),
    season = season_of(month),
    delta_pet_mm_day = pet_pm_mm_day - pet_era5_mm_day,
    delta_pet_mm_month = delta_pet_mm_day * days_in_month,
    in_P1 = year >= P1_START & year <= P1_END,
    in_P3 = year >= P3_START & year <= P3_END
  )

if (nrow(base) != 912L) warning("Expected 912 monthly observations; found ", nrow(base), ".")
if (any(!is.finite(base$pet_pm_mm_day)) || any(!is.finite(base$pet_era5_mm_day)))
  stop("Non-finite PET values detected in the PET-source comparison.")

cat(sprintf("  Record: %s to %s (%d months)\n", min(base$date), max(base$date), nrow(base)))
cat(sprintf("  PM PET mean: %.4f mm/day\n", mean(base$pet_pm_mm_day)))
cat(sprintf("  ERA5-Land PEV mean: %.4f mm/day\n", mean(base$pet_era5_mm_day)))
cat(sprintf("  PM - ERA5-Land PEV: %.4f mm/day\n", mean(base$delta_pet_mm_day)))

# ------------------------------------------------------------------------------
# Monthly and annual summaries
# ------------------------------------------------------------------------------
monthly_summary <- base %>%
  dplyr::group_by(month, season) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_pm_mm_day = mean(pet_pm_mm_day),
    mean_era5_mm_day = mean(pet_era5_mm_day),
    mean_delta_pet_mm_day = mean(delta_pet_mm_day),
    sd_delta_pet_mm_day = sd(delta_pet_mm_day),
    .groups = "drop"
  )

annual_summary <- base %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    mean_pm_mm_day = mean(pet_pm_mm_day),
    mean_era5_mm_day = mean(pet_era5_mm_day),
    mean_delta_pet_mm_day = mean(delta_pet_mm_day),
    mean_delta_pet_mm_year_mm = sum(delta_pet_mm_month),
    .groups = "drop"
  )

seasonal_summary <- base %>%
  dplyr::group_by(season) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_delta_pet_mm_day = mean(delta_pet_mm_day),
    median_delta_pet_mm_day = median(delta_pet_mm_day),
    sd_delta_pet_mm_day = sd(delta_pet_mm_day),
    mean_abs_delta_pet_mm_day = mean(abs(delta_pet_mm_day)),
    .groups = "drop"
  )

write.csv(base, file.path(OUT_DIR, "bias_monthly_timeseries.csv"), row.names = FALSE)
write.csv(annual_summary, file.path(OUT_DIR, "bias_annual_summary.csv"), row.names = FALSE)
write.csv(seasonal_summary, file.path(OUT_DIR, "bias_seasonal_summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# HAC trend inference on annual PM-ERA5 difference
# ------------------------------------------------------------------------------
trend_fit <- lm(mean_delta_pet_mm_day ~ year, data = annual_summary)
hac <- sandwich::NeweyWest(trend_fit, lag = 3, prewhite = FALSE, adjust = TRUE)
coef_tab <- lmtest::coeftest(trend_fit, vcov. = hac)

trend_tests <- data.frame(
  variable = "annual_mean_delta_pet_PM_minus_ERA5Land",
  n_years = nrow(annual_summary),
  slope_mm_day_per_year = unname(coef(trend_fit)["year"]),
  hac_se = unname(coef_tab["year", "Std. Error"]),
  hac_t = unname(coef_tab["year", "t value"]),
  hac_p = unname(coef_tab["year", "Pr(>|t|)"]),
  period_1_mean_mm_day = mean(annual_summary$mean_delta_pet_mm_day[annual_summary$year >= P1_START & annual_summary$year <= P1_END]),
  period_3_mean_mm_day = mean(annual_summary$mean_delta_pet_mm_day[annual_summary$year >= P3_START & annual_summary$year <= P3_END])
)
write.csv(trend_tests, file.path(OUT_DIR, "bias_trend_tests.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------
pdf(file.path(OUT_DIR, "bias_nonstationarity.pdf"), width = 14, height = 9)

base_plot <- base %>%
  dplyr::mutate(delta_pet_12ma = as.numeric(stats::filter(delta_pet_mm_month, rep(1/12, 12), sides = 2)))
p1_m <- mean(base$delta_pet_mm_month[base$in_P1])
p3_m <- mean(base$delta_pet_mm_month[base$in_P3])

print(
  ggplot(base_plot, aes(date, delta_pet_mm_month)) +
    annotate("rect", xmin = as.Date("1950-01-01"), xmax = as.Date("1990-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.07) +
    annotate("rect", xmin = as.Date("2022-01-01"), xmax = as.Date("2025-12-31"), ymin = -Inf, ymax = Inf, alpha = 0.15) +
    geom_line(linewidth = 0.35) +
    geom_line(aes(y = delta_pet_12ma), linewidth = 1.0, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
    labs(title = "Monthly PET-source difference: FAO-56 PM − ERA5-Land PEV",
         subtitle = "Difference shown in mm/month; thick line = 12-month moving mean",
         x = "Date", y = "PM − ERA5-Land PEV (mm/month)") +
    theme_bw(base_size = 11)
)

print(
  ggplot(annual_summary, aes(year, mean_delta_pet_mm_day)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    labs(title = "Annual PET-source difference",
         subtitle = sprintf("HAC slope = %+.5f mm/day/year; p = %.4f", trend_tests$slope_mm_day_per_year, trend_tests$hac_p),
         x = "Year", y = "PM − ERA5-Land PEV (mm/day)") +
    theme_bw(base_size = 11)
)

print(
  ggplot(seasonal_summary, aes(season, mean_delta_pet_mm_day)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_col() +
    labs(title = "Seasonal PET-source difference",
         x = "Season", y = "PM − ERA5-Land PEV (mm/day)") +
    theme_bw(base_size = 11)
)

graphics.off()

# ------------------------------------------------------------------------------
# Manuscript-ready workbook
# ------------------------------------------------------------------------------
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Monthly")
openxlsx::addWorksheet(wb, "Annual")
openxlsx::addWorksheet(wb, "Seasonal")
openxlsx::addWorksheet(wb, "Trend")
openxlsx::writeData(wb, "Monthly", monthly_summary)
openxlsx::writeData(wb, "Annual", annual_summary)
openxlsx::writeData(wb, "Seasonal", seasonal_summary)
openxlsx::writeData(wb, "Trend", trend_tests)
openxlsx::saveWorkbook(wb, file.path(OUT_DIR, "bias_manuscript_tables.xlsx"), overwrite = TRUE)

cat("\n==============================================================\n")
cat("PET-source non-stationarity analysis complete.\n")
cat("PM = FAO-56 Penman-Monteith; ERA5 = native ERA5-Land PEV.\n")
cat(sprintf("Annual HAC trend in PM-ERA5 difference: %+.5f mm/day/year; p=%.4f\n",
            trend_tests$slope_mm_day_per_year, trend_tests$hac_p))
cat("Outputs: ", OUT_DIR, "\n", sep = "")
cat("==============================================================\n")
