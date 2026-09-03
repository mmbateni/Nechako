################################################################################
# 16climatology_water_budget.R
# SUPERVISOR ADDITIONS — CLIMATOLOGY + RECENT DROUGHT WATER BUDGET
#
# Run AFTER 4pr_pet_trends.r has produced:
#   trend_analysis_pr_pet/analysis_metadata.rds
#
# Adds:
#   1. 1950-2025 monthly, seasonal and annual climatology for:
#      temperature, precipitation, AET, PET(PM), PET(ERA5-Land PEV if available)
#   2. 2022-2025 annual, seasonal and monthly climate summaries.
#   3. Water-budget ratios:
#      AET/P, PET/P, AET/PET, (P-AET)/P and (P-PET)/P.
#   4. A diagnostic residual runoff-ratio proxy:
#      R_residual/P = (P-AET)/P.
#      This is NOT observed streamflow runoff ratio because storage change and
#      other basin water-balance terms are omitted.
#   5. Trend statistics for annual water-budget quantities and ratios.
#
# Important methodological distinction:
#   Existing 4pr_pet_trends.r stores annual PM-PET as a warm-season (Apr-Oct)
#   sum for PET trend analysis. This script deliberately computes FULL-YEAR
#   PET totals from the monthly record for climatological water-budget use.
################################################################################

setwd(Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(Kendall)
  library(trend)
})

WD_PATH <- Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/")
setwd(WD_PATH)

IN_DIR  <- file.path(WD_PATH, "trend_analysis_pr_pet")
OUT_DIR <- file.path(WD_PATH, "temporal_drought", "climatology_water_budget")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

META_FILE <- file.path(IN_DIR, "analysis_metadata.rds")
if (!file.exists(META_FILE)) {
  stop("Missing ", META_FILE, ". Run 4pr_pet_trends.r first.")
}

metadata <- readRDS(META_FILE)
if (is.null(metadata$basin_avg_monthly))
  stop("analysis_metadata.rds does not contain basin_avg_monthly.")

d <- as.data.table(metadata$basin_avg_monthly)
setnames(d, trimws(names(d)))

req <- c("date", "precip_mm_month", "pet_mm_month", "aet_mm_month",
         "tair_degC_month")
miss <- setdiff(req, names(d))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

d[, date := as.Date(date)]
setorder(d, date)
d[, year := year(date)]
d[, month := month(date)]
d[, season := fifelse(month %in% c(12L, 1L, 2L), "DJF",
              fifelse(month %in% 3:5, "MAM",
              fifelse(month %in% 6:8, "JJA", "SON")))]

# Optional ERA5-Land native PEV.  It is not currently retained in the 4pr
# metadata object, so recover it from the authoritative monthly summary when
# available.  Matching is by year-month; no temporal interpolation is used.
era5_summary_candidates <- c(
  file.path(WD_PATH, "monthly_data_direct", "ERA5Land_Nechako_PET_ERA5Land_monthly_summary.csv"),
  file.path(WD_PATH, "monthly_data_direct", "ERA5Land_Nechako_PET_ERA5Land_monthly_summary.csv")
)
era5_file <- era5_summary_candidates[file.exists(era5_summary_candidates)][1]

if (length(era5_file) == 1L && !is.na(era5_file)) {
  e <- fread(era5_file)
  if (all(c("year", "month", "mean_pet") %in% names(e))) {
    e <- e[, .(year = as.integer(year), month = as.integer(month),
              pet_era5_mm_day = as.numeric(mean_pet))]
    # The summary is in mm/day. Convert to monthly depth using the actual
    # calendar month length.
    e[, days := days_in_month(as.Date(sprintf("%04d-%02d-01", year, month)))]
    e[, pet_era5_mm_month := pet_era5_mm_day * days]
    d <- merge(d, e[, .(year, month, pet_era5_mm_month)],
               by = c("year", "month"), all.x = TRUE, sort = FALSE)
    setorder(d, date)
    has_era5 <- TRUE
  } else {
    has_era5 <- FALSE
    warning("ERA5 PEV summary found but required columns are absent; skipping PET_ERA5.")
  }
} else {
  has_era5 <- FALSE
  message("ERA5-Land native PEV monthly summary not found; PET_ERA5 summaries skipped.")
}

# -----------------------------------------------------------------------------
# Monthly water-budget quantities and ratios
# -----------------------------------------------------------------------------
d[, `:=`(
  AET_over_PET = fifelse(is.finite(aet_mm_month) & is.finite(pet_mm_month) & pet_mm_month > 0,
                         aet_mm_month / pet_mm_month, NA_real_),
  AET_over_P   = fifelse(is.finite(precip_mm_month) & precip_mm_month > 0,
                         aet_mm_month / precip_mm_month, NA_real_),
  PET_over_P   = fifelse(is.finite(precip_mm_month) & precip_mm_month > 0,
                         pet_mm_month / precip_mm_month, NA_real_),
  residual_runoff_over_P = fifelse(is.finite(precip_mm_month) & precip_mm_month > 0,
                                   (precip_mm_month - aet_mm_month) / precip_mm_month,
                                   NA_real_),
  residual_water_balance_mm = precip_mm_month - aet_mm_month,
  atmospheric_demand_deficit_mm = pet_mm_month - aet_mm_month
)]

if (has_era5) {
  d[, PET_ERA5_over_P := fifelse(is.finite(precip_mm_month) & precip_mm_month > 0,
                                 pet_era5_mm_month / precip_mm_month, NA_real_)]
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
weighted_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

summarise_period <- function(x) {
  # x contains complete monthly basin-average depths. Sums are appropriate for
  # water-depth fluxes; temperature is averaged.
  out <- list(
    precipitation_mm = sum(x$precip_mm_month, na.rm = TRUE),
    AET_mm = sum(x$aet_mm_month, na.rm = TRUE),
    PET_PM_mm = sum(x$pet_mm_month, na.rm = TRUE),
    temperature_degC = weighted_mean(x$tair_degC_month)
  )
  if (has_era5)
    out$PET_ERA5_mm <- sum(x$pet_era5_mm_month, na.rm = TRUE)
  out$P_minus_AET_mm <- out$precipitation_mm - out$AET_mm
  out$P_minus_PET_PM_mm <- out$precipitation_mm - out$PET_PM_mm
  out$AET_over_P <- if (out$precipitation_mm > 0) out$AET_mm / out$precipitation_mm else NA_real_
  out$PET_PM_over_P <- if (out$precipitation_mm > 0) out$PET_PM_mm / out$precipitation_mm else NA_real_
  out$AET_over_PET_PM <- if (out$PET_PM_mm > 0) out$AET_mm / out$PET_PM_mm else NA_real_
  out$residual_runoff_over_P <- if (out$precipitation_mm > 0)
    out$P_minus_AET_mm / out$precipitation_mm else NA_real_
  out
}

# -----------------------------------------------------------------------------
# 1. MONTHLY CLIMATOLOGY: 1950-2025
# -----------------------------------------------------------------------------
# data.table's grouped result can retain duplicate month metadata; normalize.
monthly_clim <- d[, {
  list(
    month_name = month.name[unique(month)],
    precipitation_mm = mean(precip_mm_month, na.rm = TRUE),
    AET_mm = mean(aet_mm_month, na.rm = TRUE),
    PET_PM_mm = mean(pet_mm_month, na.rm = TRUE),
    temperature_degC = mean(tair_degC_month, na.rm = TRUE),
    P_minus_AET_mm = mean(precip_mm_month - aet_mm_month, na.rm = TRUE),
    P_minus_PET_PM_mm = mean(precip_mm_month - pet_mm_month, na.rm = TRUE),
    AET_over_P = ifelse(mean(precip_mm_month, na.rm = TRUE) > 0,
                    mean(aet_mm_month, na.rm = TRUE) / mean(precip_mm_month, na.rm = TRUE), NA_real_),
    PET_PM_over_P = ifelse(mean(precip_mm_month, na.rm = TRUE) > 0,
                           mean(pet_mm_month, na.rm = TRUE) / mean(precip_mm_month, na.rm = TRUE), NA_real_),
    AET_over_PET_PM = ifelse(mean(pet_mm_month, na.rm = TRUE) > 0,
                             mean(aet_mm_month, na.rm = TRUE) / mean(pet_mm_month, na.rm = TRUE), NA_real_),
    residual_runoff_over_P = ifelse(mean(precip_mm_month, na.rm = TRUE) > 0,
                                    mean(precip_mm_month - aet_mm_month, na.rm = TRUE) /
                                      mean(precip_mm_month, na.rm = TRUE), NA_real_),
    n_years = uniqueN(year)
  )
}, by = month]
if (has_era5) {
  ec <- d[, .(PET_ERA5_mm = mean(pet_era5_mm_month, na.rm = TRUE)), by = month]
  monthly_clim <- merge(monthly_clim, ec, by = "month", all.x = TRUE, sort = TRUE)
  monthly_clim[, PET_ERA5_over_P := d[, .(v = mean(PET_ERA5_over_P, na.rm = TRUE)), by = month]$v]
}
setorder(monthly_clim, month)
fwrite(monthly_clim, file.path(OUT_DIR, "climatology_monthly_1950_2025.csv"))

# -----------------------------------------------------------------------------
# 2. SEASONAL CLIMATOLOGY: 1950-2025
# -----------------------------------------------------------------------------
season_order <- c("DJF", "MAM", "JJA", "SON")
season_clim <- d[, {
  Pm <- mean(precip_mm_month, na.rm = TRUE) * 3
  Am <- mean(aet_mm_month, na.rm = TRUE) * 3
  Em <- mean(pet_mm_month, na.rm = TRUE) * 3
  list(
    precipitation_mm = Pm,
    AET_mm = Am,
    PET_PM_mm = Em,
    temperature_degC = mean(tair_degC_month, na.rm = TRUE),
    P_minus_AET_mm = Pm - Am,
    P_minus_PET_PM_mm = Pm - Em,
    AET_over_P = ifelse(Pm > 0, Am / Pm, NA_real_),
    PET_PM_over_P = ifelse(Pm > 0, Em / Pm, NA_real_),
    AET_over_PET_PM = ifelse(Em > 0, Am / Em, NA_real_),
    residual_runoff_over_P = ifelse(Pm > 0, (Pm - Am) / Pm, NA_real_),
    n_months = .N
  )
}, by = season]
if (has_era5) {
  ec <- d[, {
    E5 <- mean(pet_era5_mm_month, na.rm = TRUE) * 3
    list(PET_ERA5_mm = E5,
         PET_ERA5_over_P = ifelse(mean(precip_mm_month, na.rm = TRUE) * 3 > 0,
                                  E5 / (mean(precip_mm_month, na.rm = TRUE) * 3), NA_real_))
  }, by = season]
  season_clim <- merge(season_clim, ec, by = "season", all.x = TRUE, sort = FALSE)
}
season_clim[, season := factor(season, levels = season_order)]
setorder(season_clim, season)
season_clim[, season := as.character(season)]
fwrite(season_clim, file.path(OUT_DIR, "climatology_seasonal_1950_2025.csv"))

# -----------------------------------------------------------------------------
# 3. ANNUAL CLIMATOLOGY / FULL-YEAR SERIES: 1950-2025
# -----------------------------------------------------------------------------
annual <- d[, .(
  precipitation_mm = sum(precip_mm_month, na.rm = TRUE),
  AET_mm = sum(aet_mm_month, na.rm = TRUE),
  PET_PM_mm = sum(pet_mm_month, na.rm = TRUE),
  temperature_degC = mean(tair_degC_month, na.rm = TRUE),
  P_minus_AET_mm = sum(precip_mm_month - aet_mm_month, na.rm = TRUE),
  P_minus_PET_PM_mm = sum(precip_mm_month - pet_mm_month, na.rm = TRUE),
  AET_over_P = ifelse(sum(precip_mm_month, na.rm = TRUE) > 0,
                      sum(aet_mm_month, na.rm = TRUE) / sum(precip_mm_month, na.rm = TRUE), NA_real_),
  PET_PM_over_P = ifelse(sum(precip_mm_month, na.rm = TRUE) > 0,
                         sum(pet_mm_month, na.rm = TRUE) / sum(precip_mm_month, na.rm = TRUE), NA_real_),
  AET_over_PET_PM = ifelse(sum(pet_mm_month, na.rm = TRUE) > 0,
                           sum(aet_mm_month, na.rm = TRUE) / sum(pet_mm_month, na.rm = TRUE), NA_real_),
  residual_runoff_over_P = ifelse(sum(precip_mm_month, na.rm = TRUE) > 0,
                                  sum(precip_mm_month - aet_mm_month, na.rm = TRUE) /
                                    sum(precip_mm_month, na.rm = TRUE), NA_real_),
  n_months = .N
), by = year]
if (has_era5) {
  ec <- d[, .(PET_ERA5_mm = sum(pet_era5_mm_month, na.rm = TRUE)), by = year]
  annual <- merge(annual, ec, by = "year", all.x = TRUE, sort = TRUE)
  annual[, PET_ERA5_over_P := fifelse(precipitation_mm > 0, PET_ERA5_mm / precipitation_mm, NA_real_)]
}

# Full-year annual climatology statistics, not a time series.
annual_clim <- data.table::rbindlist(lapply(names(annual)[names(annual) != "year"], function(v) {
  x <- annual[[v]]
  data.table(variable = v,
             mean_1950_2025 = mean(x, na.rm = TRUE),
             sd_1950_2025 = sd(x, na.rm = TRUE),
             median_1950_2025 = median(x, na.rm = TRUE),
             min_1950_2025 = min(x, na.rm = TRUE),
             max_1950_2025 = max(x, na.rm = TRUE))
}), fill = TRUE)
fwrite(annual, file.path(OUT_DIR, "annual_water_budget_1950_2025.csv"))
fwrite(annual_clim, file.path(OUT_DIR, "climatology_annual_1950_2025.csv"))

# -----------------------------------------------------------------------------
# 4. 2022-2025 MONTHLY / SEASONAL / ANNUAL TABLES
# -----------------------------------------------------------------------------
recent <- d[year %in% 2022:2025]
if (nrow(recent) != 48L)
  warning("Expected 48 monthly observations for 2022-2025; found ", nrow(recent), ".")

recent_monthly <- recent[, .(
  precipitation_mm = precip_mm_month,
  AET_mm = aet_mm_month,
  PET_PM_mm = pet_mm_month,
  temperature_degC = tair_degC_month,
  P_minus_AET_mm = precip_mm_month - aet_mm_month,
  P_minus_PET_PM_mm = precip_mm_month - pet_mm_month,
  AET_over_P = AET_over_P,
  PET_PM_over_P = PET_over_P,
  AET_over_PET_PM = AET_over_PET,
  residual_runoff_over_P = residual_runoff_over_P
), by = .(year, month, date)]
if (has_era5) {
  recent_monthly[, `:=`(PET_ERA5_mm = recent$pet_era5_mm_month,
                        PET_ERA5_over_P = recent$PET_ERA5_over_P)]
}
fwrite(recent_monthly, file.path(OUT_DIR, "recent_2022_2025_monthly_climate_water_budget.csv"))

recent_seasonal <- recent[, .(
  precipitation_mm = sum(precip_mm_month, na.rm = TRUE),
  AET_mm = sum(aet_mm_month, na.rm = TRUE),
  PET_PM_mm = sum(pet_mm_month, na.rm = TRUE),
  temperature_degC = mean(tair_degC_month, na.rm = TRUE),
  P_minus_AET_mm = sum(precip_mm_month - aet_mm_month, na.rm = TRUE),
  P_minus_PET_PM_mm = sum(precip_mm_month - pet_mm_month, na.rm = TRUE),
  AET_over_P = ifelse(sum(precip_mm_month, na.rm=TRUE) > 0,
                      sum(aet_mm_month, na.rm=TRUE) / sum(precip_mm_month, na.rm=TRUE), NA_real_),
  PET_PM_over_P = ifelse(sum(precip_mm_month, na.rm=TRUE) > 0,
                         sum(pet_mm_month, na.rm=TRUE) / sum(precip_mm_month, na.rm=TRUE), NA_real_),
  AET_over_PET_PM = ifelse(sum(pet_mm_month, na.rm=TRUE) > 0,
                           sum(aet_mm_month, na.rm=TRUE) / sum(pet_mm_month, na.rm=TRUE), NA_real_),
  residual_runoff_over_P = ifelse(sum(precip_mm_month, na.rm=TRUE) > 0,
                                  sum(precip_mm_month - aet_mm_month, na.rm=TRUE) /
                                    sum(precip_mm_month, na.rm=TRUE), NA_real_)
), by = .(year, season)]
if (has_era5) {
  er <- recent[, .(PET_ERA5_mm = sum(pet_era5_mm_month, na.rm = TRUE)), by=.(year, season)]
  recent_seasonal <- merge(recent_seasonal, er, by=c("year","season"), all.x=TRUE, sort=FALSE)
  recent_seasonal[, PET_ERA5_over_P := fifelse(precipitation_mm > 0, PET_ERA5_mm / precipitation_mm, NA_real_)]
}
recent_seasonal[, season := factor(season, levels=season_order)]
setorder(recent_seasonal, year, season)
recent_seasonal[, season := as.character(season)]
fwrite(recent_seasonal, file.path(OUT_DIR, "recent_2022_2025_seasonal_climate_water_budget.csv"))

recent_annual <- annual[year %in% 2022:2025]
fwrite(recent_annual, file.path(OUT_DIR, "recent_2022_2025_annual_climate_water_budget.csv"))

# -----------------------------------------------------------------------------
# 5. CLIMATOLOGICAL ANOMALIES FOR 2022-2025
# -----------------------------------------------------------------------------
# Compare each recent calendar month against the corresponding 1950-2025
# calendar-month mean. This is the most interpretable way to describe whether
# the drought years were dry, warm, high-demand, or evaporation-limited.
clim_lookup <- monthly_clim[, .(
  month,
  P_clim = precipitation_mm,
  AET_clim = AET_mm,
  PET_clim = PET_PM_mm,
  T_clim = temperature_degC
)]
recent_anom <- merge(recent_monthly, clim_lookup, by="month", all.x=TRUE, sort=FALSE)
recent_anom[, `:=`(
  P_anomaly_pct = 100 * (precipitation_mm - P_clim) / P_clim,
  AET_anomaly_pct = 100 * (AET_mm - AET_clim) / AET_clim,
  PET_PM_anomaly_pct = 100 * (PET_PM_mm - PET_clim) / PET_clim,
  temperature_anomaly_degC = temperature_degC - T_clim
)]
fwrite(recent_anom, file.path(OUT_DIR, "recent_2022_2025_monthly_climatological_anomalies.csv"))

# -----------------------------------------------------------------------------
# 6. ANNUAL TREND TESTS FOR WATER-BUDGET RATIOS
# -----------------------------------------------------------------------------
# The ratios are annual aggregates (ratio of annual depths), avoiding the
# compositional bias that would arise from averaging monthly ratios when asking
# whether the annual water budget is changing.
trend_vars <- c("precipitation_mm", "AET_mm", "PET_PM_mm", "P_minus_AET_mm",
                "P_minus_PET_PM_mm", "AET_over_P", "PET_PM_over_P",
                "AET_over_PET_PM", "residual_runoff_over_P")

trend_rows <- lapply(trend_vars, function(v) {
  x <- annual[[v]]
  ok <- is.finite(x) & is.finite(annual$year)
  xx <- x[ok]
  yy <- annual$year[ok]
  if (length(xx) < 8L) {
    return(data.frame(variable=v, n=length(xx), tau=NA_real_, p=NA_real_,
                      sen_slope=NA_real_, slope_units=NA_character_))
  }
  mk <- tryCatch(Kendall::MannKendall(xx), error=function(e) NULL)
  ss <- tryCatch(trend::sens.slope(xx), error=function(e) NULL)
  data.frame(
    variable=v,
    n=length(xx),
    tau=if (!is.null(mk)) as.numeric(mk$tau) else NA_real_,
    p=if (!is.null(mk)) as.numeric(mk$sl) else NA_real_,
    sen_slope=if (!is.null(ss)) as.numeric(ss$estimates) else NA_real_,
    slope_units=if (grepl("over", v, fixed=TRUE)) "ratio/year" else "mm/year or degC/year",
    stringsAsFactors=FALSE
  )
})
trend_tbl <- rbindlist(trend_rows, fill=TRUE)
fwrite(trend_tbl, file.path(OUT_DIR, "annual_water_budget_ratio_trends_1950_2025.csv"))

# -----------------------------------------------------------------------------
# 7. Practitioner-facing comparison table
# -----------------------------------------------------------------------------
# One row per year, with the variables needed directly in the manuscript.
practitioner <- annual[, .(
  year,
  P_mm = precipitation_mm,
  AET_mm,
  PET_PM_mm,
  Tair_degC = temperature_degC,
  AET_P = AET_over_P,
  PET_P = PET_PM_over_P,
  AET_PET = AET_over_PET_PM,
  residual_runoff_P = residual_runoff_over_P,
  P_minus_AET_mm,
  P_minus_PET_mm = P_minus_PET_PM_mm
)]
if (has_era5) practitioner[, PET_ERA5_mm := annual[match(year, annual$year), PET_ERA5_mm]]
fwrite(practitioner, file.path(OUT_DIR, "practitioner_annual_water_budget_1950_2025.csv"))

# -----------------------------------------------------------------------------
# 8. Console diagnostics
# -----------------------------------------------------------------------------
cat("\n==============================================\n")
cat("CLIMATOLOGY + WATER-BUDGET ADDITIONS COMPLETE\n")
cat("==============================================\n")
cat("Output directory: ", OUT_DIR, "\n", sep="")
cat("Monthly climatology: 1950-2025\n")
cat("Seasonal climatology: 1950-2025\n")
cat("Full-year annual series: 1950-2025\n")
cat("Recent annual/seasonal/monthly tables: 2022-2025\n")
cat("AET/P, PET/P, AET/PET and residual (P-AET)/P included.\n")
cat("\nAnnual water-budget ratio trends:\n")
print(trend_tbl)

invisible(TRUE)
