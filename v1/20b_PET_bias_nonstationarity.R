# ============================================================================
# 10b_PET_bias_nonstationarity.R
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

# ==============================================================================
# REGIONAL EXTENSION — PET-PRODUCT DIFFERENCE FOR 7 SUBBASINS + WHOLE BASIN
# ==============================================================================
REGIONAL_PET_DIR <- file.path(WD_PATH, "trend_analysis_pr_pet", "regional")
REGIONAL_OUT_DIR <- file.path(WD_PATH, "decomp_results", "bias_nonstationarity", "regional")
dir.create(REGIONAL_OUT_DIR, recursive=TRUE, showWarnings=FALSE)

REGIONAL_PM_FILE <- file.path(REGIONAL_PET_DIR, "regional_PM_PET_monthly_summary.csv")
REGIONAL_E5_FILE <- file.path(REGIONAL_PET_DIR, "regional_ERA5Land_PEV_monthly_summary.csv")

if (file.exists(REGIONAL_PM_FILE) && file.exists(REGIONAL_E5_FILE)) {
  rpm <- data.table::fread(REGIONAL_PM_FILE)
  re5 <- data.table::fread(REGIONAL_E5_FILE)
  
  # --------------------------------------------------------------------------
  # Regional input validation and key normalisation
  # --------------------------------------------------------------------------
  required_regional <- c("region_id", "region_name", "date", "mean_pet")
  miss_pm <- setdiff(required_regional, names(rpm))
  miss_e5 <- setdiff(required_regional, names(re5))
  if (length(miss_pm)) {
    stop("Regional PM summary is missing columns: ", paste(miss_pm, collapse=", "))
  }
  if (length(miss_e5)) {
    stop("Regional ERA5-Land PEV summary is missing columns: ", paste(miss_e5, collapse=", "))
  }
  
  # 4pr_pet_trends.r writes dates as YYYY-MM. This parser also accepts
  # YYYY-MM-DD so that the downstream script remains robust to date-format
  # changes in the upstream regional export.
  rpm[, date := normalise_month_date(date, "Regional PM date")]
  re5[, date := normalise_month_date(date, "Regional ERA5-Land PEV date")]
  
  # Region IDs are the authoritative spatial keys. Region names are metadata
  # labels and therefore must NOT be used as join keys: a harmless naming
  # change (case, whitespace, punctuation, etc.) must not turn every valid
  # region-month pair into an unmatched outer-join row.
  rpm[, region_id := trimws(as.character(region_id))]
  re5[, region_id := trimws(as.character(region_id))]
  rpm[, region_name := trimws(as.character(region_name))]
  re5[, region_name := trimws(as.character(region_name))]
  rpm[, mean_pet := as.numeric(mean_pet)]
  re5[, mean_pet := as.numeric(mean_pet)]
  
  # --------------------------------------------------------------------------
  # Resolve duplicate region-month rows safely.
  #
  # The scientific regional key is (region_id, date). Region names and other
  # metadata are not part of that key. The upstream regional export can contain
  # repeated rows for the same region-month because of intermediate metadata
  # joins. We may collapse such rows ONLY when all duplicate rows have the same
  # finite PET value. Conflicting PET values are not silently discarded.
  # --------------------------------------------------------------------------
  collapse_regional_duplicates <- function(x, label) {
    # Duplicate region-month rows are resolved at the scientific key level.
    # A common upstream pattern is one NA row plus one finite PET row for the
    # same key. The finite observation is retained; NA does not constitute a
    # competing PET estimate.
    dup_keys <- x[, .N, by=.(region_id, date)][N > 1L]
    if (!nrow(dup_keys)) return(x)
    
    dup_detail <- x[
      dup_keys,
      on=.(region_id, date),
      allow.cartesian=TRUE,
      .(
        n_rows = .N,
        n_finite = sum(is.finite(mean_pet)),
        n_unique_finite = data.table::uniqueN(mean_pet[is.finite(mean_pet)])
      ),
      by=.(region_id, date)
    ]
    
    # A true conflict exists only when more than one distinct finite PET value
    # occurs for the same region-month. Multiple NA rows are not conflicts.
    conflicts <- dup_detail[n_unique_finite > 1L]
    
    if (nrow(conflicts)) {
      cat("\n==============================================================\n")
      cat("CONFLICTING REGIONAL PET DUPLICATES: ", label, "\n", sep="")
      cat("==============================================================\n")
      
      print(
        x[
          conflicts,
          on=.(region_id, date),
          allow.cartesian=TRUE,
          .(region_id, region_name, date, mean_pet),
          nomatch=0L
        ][order(region_id, date, mean_pet)]
      )
      
      stop(
        label,
        " contains duplicate region-month observations with multiple distinct finite PET values. ",
        "Conflicting keys: ", nrow(conflicts), ". ",
        "The upstream regional aggregation must be corrected; 10b will not choose between conflicting PET values."
      )
    }
    
    # Collapse each duplicate key deterministically:
    #   - retain the sole finite PET value if one exists;
    #   - otherwise retain NA;
    #   - retain the first non-empty region name.
    collapsed <- x[
      , {
        finite_pet <- mean_pet[is.finite(mean_pet)]
        name_values <- region_name[!is.na(region_name) & nzchar(region_name)]
        list(
          region_name = if (length(name_values)) name_values[1L] else NA_character_,
          mean_pet = if (length(finite_pet)) finite_pet[1L] else NA_real_
        )
      },
      by=.(region_id, date)
    ]
    
    # Preserve the original temporal fields when present; these are later
    # re-derived from date in reg_base.
    if ("year" %in% names(x)) {
      collapsed[, year := as.integer(format(date, "%Y"))]
    }
    if ("month" %in% names(x)) {
      collapsed[, month := as.integer(format(date, "%m"))]
    }
    
    data.table::setcolorder(
      collapsed,
      c("region_id", "region_name", "date",
        setdiff(names(collapsed), c("region_id", "region_name", "date")))
    )
    
    n_removed <- nrow(x) - nrow(collapsed)
    warning(
      label,
      ": collapsed ", n_removed,
      " redundant duplicate region-month row(s). " ,
      "Rows containing only NA alongside a finite PET observation were resolved in favour of the finite observation."
    )
    
    collapsed
  }
  
  rpm <- collapse_regional_duplicates(rpm, "Regional PM summary")
  re5 <- collapse_regional_duplicates(re5, "Regional ERA5-Land PEV summary")
  
  # Verify that the two products cover the same regional keys before merging.
  pm_keys <- unique(rpm[, .(region_id, date)])
  e5_keys <- unique(re5[, .(region_id, date)])
  common_keys <- merge(pm_keys, e5_keys,
                       by=c("region_id", "date"), all=FALSE, sort=FALSE)
  
  if (nrow(common_keys) == 0L) {
    stop(
      paste0(
        "Regional PM and ERA5-Land PEV summaries have no matching region_id-date keys. ",
        "PM rows = ", nrow(rpm), "; ERA5-Land rows = ", nrow(re5), ". ",
        "Check that both files were generated from the same 4pr_pet_trends.r run " ,
        "and that region_id/date values are consistent."
      )
    )
  }
  
  missing_from_e5 <- pm_keys[!e5_keys, on=.(region_id, date)]
  missing_from_pm <- e5_keys[!pm_keys, on=.(region_id, date)]
  
  if (nrow(missing_from_e5) || nrow(missing_from_pm)) {
    stop(
      paste0(
        "Regional PET summaries do not have identical region-month coverage. ",
        "PM-only keys = ", nrow(missing_from_e5), "; ",
        "ERA5-only keys = ", nrow(missing_from_pm), ". ",
        "Both regional products must contain the same region_id-date combinations."
      )
    )
  }
  
  # Region name is metadata, not a spatial join key. Retain the PM name as the
  # canonical label and audit any name discrepancies separately.
  name_check <- merge(
    unique(rpm[, .(region_id, region_name_pm=region_name)]),
    unique(re5[, .(region_id, region_name_e5=region_name)]),
    by="region_id", all=TRUE, sort=TRUE
  )
  name_mismatch <- name_check[
    is.na(region_name_pm) | is.na(region_name_e5) |
      region_name_pm != region_name_e5
  ]
  if (nrow(name_mismatch)) {
    warning(
      "Regional PM and ERA5-Land region_name labels differ for ",
      nrow(name_mismatch),
      " region(s); region_id was used as the authoritative join key."
    )
  }
  
  reg_base <- merge(
    rpm[, .(region_id, region_name, date, year, month,
            pet_pm_mm_day = mean_pet)],
    re5[, .(region_id, date, pet_era5_mm_day = mean_pet)],
    by=c("region_id", "date"), all=FALSE, sort=TRUE
  )
  
  # Re-derive time metadata from the canonical monthly date. This prevents
  # stale or source-specific year/month columns from entering the analysis.
  reg_base[, `:=`(
    year = as.integer(format(date, "%Y")),
    month = as.integer(format(date, "%m"))
  )]
  
  # Create delta_pet_mm_day before referencing it in the monthly total.
  # In data.table, RHS expressions in a multi-column := assignment should not
  # depend on another column being created by that same assignment.
  reg_base[, `:=`(
    days_in_month = lubridate::days_in_month(date),
    season = season_of(month),
    delta_pet_mm_day = pet_pm_mm_day - pet_era5_mm_day,
    in_P1 = year >= P1_START & year <= P1_END,
    in_P3 = year >= P3_START & year <= P3_END
  )]
  
  reg_base[, delta_pet_mm_month := delta_pet_mm_day * days_in_month]
  
  if (any(!is.finite(reg_base$pet_pm_mm_day)) ||
      any(!is.finite(reg_base$pet_era5_mm_day))) {
    stop(
      paste0(
        "Regional PET summaries contain missing/non-finite PM or ERA5-Land PEV values. ",
        "Missing PM values: ", sum(!is.finite(reg_base$pet_pm_mm_day)),
        "; missing ERA5-Land values: ", sum(!is.finite(reg_base$pet_era5_mm_day)),
        ". Check the regional PM and ERA5-Land summary files for non-finite PET values."
      )
    )
  }
  
  reg_monthly <- reg_base[, .(
    n=.N,
    mean_pm_mm_day=mean(pet_pm_mm_day),
    mean_era5_mm_day=mean(pet_era5_mm_day),
    mean_delta_pet_mm_day=mean(delta_pet_mm_day),
    sd_delta_pet_mm_day=sd(delta_pet_mm_day)
  ), by=.(region_id,region_name,month,season)]
  reg_annual <- reg_base[, .(
    mean_pm_mm_day=mean(pet_pm_mm_day),
    mean_era5_mm_day=mean(pet_era5_mm_day),
    mean_delta_pet_mm_day=mean(delta_pet_mm_day),
    mean_delta_pet_mm_year_mm=sum(delta_pet_mm_month)
  ), by=.(region_id,region_name,year)]
  
  reg_trend <- reg_annual[, {
    if (.N < 3L ||
        any(!is.finite(year)) ||
        any(!is.finite(mean_delta_pet_mm_day))) {
      list(
        n_years = .N,
        slope_mm_day_per_year = NA_real_,
        hac_se = NA_real_,
        hac_t = NA_real_,
        hac_p = NA_real_
      )
    } else {
      fit <- tryCatch(
        lm(mean_delta_pet_mm_day ~ year, data=.SD),
        error=function(e) NULL
      )
      
      if (is.null(fit)) {
        list(
          n_years = .N,
          slope_mm_day_per_year = NA_real_,
          hac_se = NA_real_,
          hac_t = NA_real_,
          hac_p = NA_real_
        )
      } else {
        hc <- tryCatch(
          lmtest::coeftest(
            fit,
            vcov.=sandwich::NeweyWest(
              fit,
              lag=3,
              prewhite=FALSE,
              adjust=TRUE
            )
          ),
          error=function(e) NULL
        )
        
        list(
          n_years = .N,
          slope_mm_day_per_year = unname(coef(fit)["year"]),
          hac_se = if(!is.null(hc)) unname(hc["year", "Std. Error"]) else NA_real_,
          hac_t = if(!is.null(hc)) unname(hc["year", "t value"]) else NA_real_,
          hac_p = if(!is.null(hc)) unname(hc["year", "Pr(>|t|)"]) else NA_real_
        )
      }
    }
  }, by=.(region_id,region_name)]
  
  data.table::fwrite(reg_base, file.path(REGIONAL_OUT_DIR,"bias_monthly_timeseries.csv"))
  data.table::fwrite(reg_monthly, file.path(REGIONAL_OUT_DIR,"bias_monthly_summary.csv"))
  data.table::fwrite(reg_annual, file.path(REGIONAL_OUT_DIR,"bias_annual_summary.csv"))
  data.table::fwrite(reg_trend, file.path(REGIONAL_OUT_DIR,"bias_trend_tests.csv"))
  
  p <- ggplot2::ggplot(reg_base, ggplot2::aes(date, delta_pet_mm_month)) +
    ggplot2::geom_hline(yintercept=0,linetype="dashed",linewidth=.3) +
    ggplot2::geom_line(linewidth=.25,na.rm=TRUE) +
    ggplot2::facet_wrap(~region_name,ncol=3,scales="free_y") +
    ggplot2::labs(title="Regional PET-source difference: FAO-56 PM − ERA5-Land PEV",
                  x=NULL,y="PM − ERA5-Land PEV (mm/month)") +
    ggplot2::theme_bw(base_size=9)
  ggplot2::ggsave(file.path(REGIONAL_OUT_DIR,"bias_nonstationarity_regional.png"),p,width=14,height=12,dpi=300)
  ggplot2::ggsave(file.path(REGIONAL_OUT_DIR,"bias_nonstationarity_regional.pdf"),p,width=14,height=12)
  cat("✓ Regional PET-source/non-stationarity products written.\n")
} else {
  warning("Regional PET summaries are absent. Run 4pr_pet_trends.r before the regional 10b extension.")
}
