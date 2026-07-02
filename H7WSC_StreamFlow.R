# ============================================================================
# NECHAKO STREAMFLOW RIVER BASIN — SSI DROUGHT ANALYSIS
# Pipeline B (Monthly): Standardized Streamflow Index, Peña-Angulo et al. (2022)
#   Stations: 08JC001 NECHAKO RIVER AT VANDERHOOF — NATURALIZED flow (observed
#             pre-regulation + naturalized/de-regulated post-regulation,
#             pre-merged into a single CSV, loaded from NATURALIZED_DATA_DIR)
#             08JE001 STUART RIVER NEAR FORT ST. JAMES — observed flow
# Pipeline C (Monthly): Same SSI methodology as Pipeline B, but on the
#   ORIGINAL (observed, non-naturalized) record instead of the naturalized one
#   Stations: 08JC001 NECHAKO RIVER AT VANDERHOOF — ORIGINAL/observed flow,
#             restricted to the POST-REGULATION period only (years >=
#             POST_REGULATION_START_YEAR); the pre-regulation portion is
#             dropped and no naturalization is applied
#             08JE001 STUART RIVER NEAR FORT ST. JAMES — observed flow (same
#             series as Pipeline B)
#
# Shared preprocessing:
#   - Data loading & completeness check  (original year-span logic)
#   - Gap filling: Fritsch-Carlson PCHIP, gaps > MAX_FILL_DAYS left as NA  
#   - Hard exclusion gate for stations < 50% completeness                  
#
# Both pipelines add:
#   - Daily → monthly mean aggregation
#   - Best-fit distribution selection via L-moments (lmomco), fit SEPARATELY
#     per calendar month (Jan..Dec each get their own candidate-distribution
#     contest, own GPA disqualification check, own close-margin bootstrap,
#     and own CDF) rather than pooling all months into one fit. This mirrors
#     the seasonal/day-of-year conditioning used by the retired daily-threshold
#     approach and matches standard SPI/SSI practice (McKee et al. 1993;
#     Vicente-Serrano et al. 2012; Peña-Angulo et al. 2022) — a January flow
#     is judged only against other Januaries, not pooled with June snowmelt.
#     Requires >= MIN_YEARS_PER_MONTH observations of that specific calendar
#     month; months that don't meet this are left as NA in the SSI series
#     rather than borrowing strength from other months.
#   - SSI calculation and drought identification (SSI < -0.50)
#
# check_data_completeness() uses the original month-span denominator.
#
# Pipeline B additionally adds:
#   - Homogeneity screening (Pettitt test) on the annual-mean series of the
#     naturalized-flow station (08JC001), since the observed/naturalized
#     splice is merged upstream of this script with no marker of where the
#     transition occurs. A significant result flags that station's "single
#     homogeneous record" assumption as questionable (see
#     check_naturalized_homogeneity()); it does not auto-exclude the station,
#     but the flag is carried through to the summary CSV, Excel workbook, and
#     saved RDS objects so it isn't silently lost. This check does NOT apply
#     to Pipeline C's 08JC001, since that series is the real observed record
#     (no splice) simply restricted to post-regulation years.
# ============================================================================

rm(list = ls())
if (file.exists("DROUGHT_ANALYSIS_utils.R")) {
  source("DROUGHT_ANALYSIS_utils.R")
} else {
  warning("DROUGHT_ANALYSIS_utils.R not found. Falling back to local thresholds.")
}
# ============================================================================
# SECTION 1: PACKAGES
# ============================================================================
packages_needed <- c(
  # Shared / general
  "tidyverse", "lubridate", "zoo",
  # Pipeline B (SSI)
  "lmomco", "fitdistrplus",
  # Pipeline B (naturalized-flow homogeneity screening for 08JC001)
  "trend",
  # Additional outputs (from SWEI code)
  "writexl"
)

for (pkg in packages_needed) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# SECTION 2: STATION METADATA & DIRECTORIES
# ============================================================================
# ── Pipeline B stations: 08JC001 NATURALIZED flow + 08JE001 observed ───────
# 08JC001 uses the NATURALIZED (de-regulated) flow record: observed discharge
# (pre-regulation) blended with naturalized flow (post-regulation), pre-merged
# into a single CSV upstream of this script, loaded from NATURALIZED_DATA_DIR
# (see get_data_path_for_station()). StartYear/EndYear are left NA here and
# populated automatically from the data once loaded (see PROCESS ALL
# STATIONS below). Because the observed/naturalized splice has no marker of
# where the transition occurs, 08JC001 is screened for an undocumented
# step-change on its annual-mean series using a Pettitt test — see
# check_naturalized_homogeneity().
stations <- data.frame(
  StationID   = c("08JC001", "08JE001"),
  StationName = c("NECHAKO RIVER AT VANDERHOOF",
                  "STUART RIVER NEAR FORT ST. JAMES"),
  StartYear   = c(NA_real_, 1929),
  EndYear     = c(NA_real_, 2024),
  stringsAsFactors = FALSE
)

# ── Pipeline C stations: 08JC001 ORIGINAL (observed) flow, post-regulation
#    period only + 08JE001 observed ─────────────────────────────────────────
# 08JC001 here is loaded from the standard INPUT_DATA_DIR (i.e. the actual
# WSC observed record, NOT the naturalized/de-regulated reconstruction used
# in Pipeline B). That observed record spans both pre- and post-regulation
# years; Pipeline C keeps ONLY the post-regulation portion (year >=
# POST_REGULATION_START_YEAR, see SECTION 3) and drops the pre-regulation
# years, since the goal here is the real regulated-river signal, not a
# naturalized one. No homogeneity/splice screening is needed for this
# series, since it isn't a blended/spliced record.
stations_C <- data.frame(
  StationID   = c("08JC001", "08JE001"),
  StationName = c("NECHAKO RIVER AT VANDERHOOF",
                  "STUART RIVER NEAR FORT ST. JAMES"),
  StartYear   = c(NA_real_, 1929),
  EndYear     = c(NA_real_, 2024),
  stringsAsFactors = FALSE
)

setwd("D:/Nechako_Drought/Nechako")

INPUT_DATA_DIR  <- paste0("D:/Nechako_Drought/Nechako/Hydrology/",
                          "data_retrievalWaterSurveyofCanada/data_downloads_geomet_api")

# Naturalized flow file (dat_natural_08JC001.csv), used by Pipeline B only:
# pre-regulation observed + post-regulation naturalized flow, pre-merged
NATURALIZED_DATA_DIR <- paste0("D:/Nechako_Drought/Nechako/Hydrology/",
                               "data_retrievalWaterSurveyofCanada/naturalized_flows")

# Single unified output root — all pipelines write under their own subdirs
MAIN_OUTPUT_DIR <- "streamflow_results"

# STATION DISPLAY HELPER
# Defaults to stations_C; since stations and stations_C share the same
# StationIDs/StationNames (08JC001, 08JE001) here, either would resolve names
# correctly.
get_station_display <- function(station_id, stations_df = stations_C) {
  name <- stations_df$StationName[stations_df$StationID == station_id]
  if (length(name) == 0) name <- "UNKNOWN"
  return(paste0(station_id, " - ", name))
}

for (sub in c(
  "",
  "shared/processed_data",          # daily gap-filled data
  "shared/completeness",
  "ssi/monthly_data",               # Pipeline B
  "ssi/drought_events",
  "ssi_C/monthly_data",              # Pipeline C (08JC001 original/post-regulation + 08JE001)
  "ssi_C/drought_events",
  "figures",
  "reports"
)) {
  dir.create(file.path(MAIN_OUTPUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

# ============================================================================
# SECTION 2.2: STATION COMPLETENESS PLOT (from SWEI code)
# ============================================================================
create_station_completeness_plot <- function(stations, completeness, output_dir) {
  
  # Create station data frame with status
  station_status <- data.frame(
    StationID    = names(completeness),
    Completeness = sapply(completeness, function(x) x$completeness),
    Status       = sapply(completeness, function(x) ifelse(x$meets_threshold, "INCLUDED", "EXCLUDED")),
    stringsAsFactors = FALSE
  )
  
  # Merge with station metadata
  station_map_data <- merge(stations, station_status, by = "StationID")
  
  cat("  [PLOT] Station status summary:\n")
  cat(sprintf("    Included: %d stations\n", sum(station_map_data$Status == "INCLUDED")))
  cat(sprintf("    Excluded: %d stations\n", sum(station_map_data$Status == "EXCLUDED")))
  
  # Create bar chart of completeness
  png_file <- file.path(output_dir, "figures", "station_completeness_barplot.png")
  png(png_file, width = 1600, height = 1000, res = 150)
  par(mar = c(8, 5, 4, 2))
  bar_colors <- ifelse(station_map_data$Status == "INCLUDED", "#2ca02c", "#d62728")
  barplot(station_map_data$Completeness,
          names.arg = paste0(station_map_data$StationID, "\n", 
                             substr(station_map_data$StationName, 1, 25)),
          col = bar_colors,
          las = 2,
          cex.names = 0.7,
          ylim = c(0, 100),
          main = "Station Data Completeness (Year-Span Denominator)",
          ylab = "Completeness (%)",
          border = "gray30")
  abline(h = 50, lty = 2, lwd = 2, col = "red")
  legend("topright", 
         legend = c("INCLUDED (≥50%)", "EXCLUDED (<50%)"),
         fill = c("#2ca02c", "#d62728"),
         bty = "n")
  dev.off()
  
  cat(sprintf("  [PLOT] Station completeness plot saved: %s\n", basename(png_file)))
  return(png_file)
}

# ============================================================================
# SECTION 2.3: EXCEL SUMMARY WORKBOOK 
# ============================================================================
create_excel_summary <- function(stations, completeness,
                                 daily_droughts, ssi_droughts,
                                 ssi_results, output_dir) {
  if (!require("writexl", quietly = TRUE)) {
    cat("  [EXCEL] writexl package not available - skipping Excel output\n")
    return(NULL)
  }
  
  excel_data <- list()
  
  # Sheet 1: Station Summary - ALL stations get ALL columns
  summary_rows <- list()
  for (sid in names(completeness)) {
    comp <- completeness[[sid]]
    sname <- if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN"
    
    # Base columns for ALL stations
    base <- data.frame(
      StationID        = sid,
      StationName      = sname,
      StationDisplay   = paste0(sid, " - ", sname),
      Period           = sprintf("%d-%d", comp$start_year, comp$end_year),
      Completeness_pct = round(comp$completeness, 1),
      Status           = ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED"),
      stringsAsFactors = FALSE
    )
    
    # Pipeline B stats - ALWAYS add these columns (NA if excluded)
    if (comp$meets_threshold && !is.null(ssi_results[[sid]])) {
      base$B_Best_Distribution <- ssi_results[[sid]]$best_distribution_str
    } else {
      base$B_Best_Distribution <- NA_character_
    }
    
    if (comp$meets_threshold && !is.null(ssi_droughts[[sid]]$drought_events)) {
      ev_b <- ssi_droughts[[sid]]$drought_events
      base$B_Drought_Events       <- nrow(ev_b)
      base$B_Total_Months         <- sum(ev_b$duration_months)
      base$B_Avg_Duration_Months  <- if (nrow(ev_b) > 0) round(mean(ev_b$duration_months), 1) else NA_real_
      base$B_Avg_Severity         <- if (nrow(ev_b) > 0) round(mean(ev_b$severity), 2) else NA_real_
    } else {
      base$B_Drought_Events       <- NA_integer_
      base$B_Total_Months         <- NA_integer_
      base$B_Avg_Duration_Months  <- NA_real_
      base$B_Avg_Severity         <- NA_real_
    }
    
    summary_rows[[sid]] <- base
  }
  excel_data[["01_Station_Summary"]] <- do.call(rbind, summary_rows)
  
  # Sheets for Pipeline B Drought Events (one per included station)
  for (sid in names(ssi_droughts)) {
    if (!is.null(ssi_droughts[[sid]]$drought_events)) {
      ev <- ssi_droughts[[sid]]$drought_events
      if (nrow(ev) > 0) {
        excel_data[[paste0("B_", sid, "_SSI_Droughts")]] <- as.data.frame(ev)
      }
    }
  }
  
  # Distribution Fit Results
  dist_rows <- list()
  for (sid in names(ssi_results)) {
    if (!is.null(ssi_results[[sid]]$dist_distances)) {
      dd <- ssi_results[[sid]]$dist_distances
      dd$StationID <- sid
      dist_rows[[sid]] <- dd
    }
  }
  if (length(dist_rows) > 0) {
    excel_data[["16_Distribution_Fits"]] <- do.call(rbind, dist_rows)
  }
  
  # Save Excel file
  xlsx_file <- file.path(output_dir, "reports", "nechako_streamflow_drought_summary.xlsx")
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  writexl::write_xlsx(excel_data, xlsx_file)
  
  cat(sprintf("  [EXCEL] Summary workbook saved: %s (%d sheets)\n",
              basename(xlsx_file), length(excel_data)))
  return(xlsx_file)
}

# ============================================================================
# SECTION 2.4: DETAILED TEXT REPORT (from SWEI code)
# ============================================================================
create_text_report <- function(stations, completeness, 
                               daily_droughts, ssi_droughts,
                               ssi_results, output_dir, params) {
  
  report_file <- file.path(output_dir, "reports", "nechako_streamflow_drought_report.txt")
  
  cat("============================================================\n",
      "NECHAKO STREAMFLOW DROUGHT ANALYSIS - DETAILED REPORT\n",
      "Dual-Pipeline Approach (Daily Threshold + Monthly SSI)\n",
      "============================================================\n\n",
      file = report_file, sep = "")
  
  cat(sprintf("Report generated: %s\n", Sys.time()),
      file = report_file, append = TRUE)
  cat(sprintf("Script version: %s\n\n", params$script_version),
      file = report_file, append = TRUE)
  
  # Station summary
  cat("========== STATION SUMMARY ==========\n",
      file = report_file, append = TRUE)
  cat(sprintf("Total stations processed: %d\n", nrow(stations)),
      file = report_file, append = TRUE)
  cat(sprintf("Stations included (≥%.0f%% completeness): %d\n", 
              params$min_completeness_pct,
              sum(sapply(completeness, function(x) x$meets_threshold))),
      file = report_file, append = TRUE)
  cat(sprintf("Stations excluded (<%.0f%% completeness): %d\n\n", 
              params$min_completeness_pct,
              sum(!sapply(completeness, function(x) x$meets_threshold))),
      file = report_file, append = TRUE)
  
  # Per-station details
  for (sid in names(completeness)) {
    comp <- completeness[[sid]]
    sname <- if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN"
    
    cat(sprintf("\n----- Station: %s -----\n", sid),
        file = report_file, append = TRUE)
    cat(sprintf("Name: %s\n", sname),
        file = report_file, append = TRUE)
    cat(sprintf("Period: %d-%d (%d years)\n", 
                comp$start_year, comp$end_year, comp$total_years),
        file = report_file, append = TRUE)
    cat(sprintf("Completeness: %.2f%% (%d/%d days)\n", 
                comp$completeness, comp$total_days, comp$expected_days),
        file = report_file, append = TRUE)
    cat(sprintf("Status: %s\n", ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED")),
        file = report_file, append = TRUE)
    
    if (comp$meets_threshold) {
      if (!is.null(ssi_droughts[[sid]]$drought_events)) {
        ev_b <- ssi_droughts[[sid]]$drought_events
        cat(sprintf("\nPipeline B (Monthly SSI):\n"),
            file = report_file, append = TRUE)
        cat(sprintf("  Best distribution (per calendar month): %s\n", 
                    ssi_results[[sid]]$best_distribution_str),
            file = report_file, append = TRUE)
        cat(sprintf("  Drought events: %d\n", nrow(ev_b)),
            file = report_file, append = TRUE)
        cat(sprintf("  Total drought months: %d\n", sum(ev_b$duration_months)),
            file = report_file, append = TRUE)
        cat(sprintf("  Average severity: %.2f\n", mean(ev_b$severity)),
            file = report_file, append = TRUE)
      }
    }
  }
  
  # Methodology
  cat("\n\n========== METHODOLOGY ==========\n",
      file = report_file, append = TRUE)
  # Pipeline A (Daily) — DISABLED
  # cat("Pipeline A (Daily):\n",
  #     "  • Variable threshold (Q20 = 80% exceedance probability)\n",
  #     "  • 31-day centered moving average smoothing\n",
  #     "  • Minimum drought duration: 30 days\n",
  #     "  • Reference: Raut & Ganguli (2024)\n\n",
  #     file = report_file, append = TRUE)
  cat("Pipeline B (Monthly):\n",
      "  • Standardized Streamflow Index (SSI)\n",
      "  • Best-fit distribution via L-moments (lmomco)\n",
      "  • Drought threshold: SSI < -0.50 (moderate drought)\n",
      "  • Reference: Peña-Angulo et al. (2022)\n\n",
      file = report_file, append = TRUE)
  
  cat("\n============================================================\n",
      "END OF REPORT\n",
      "============================================================\n",
      file = report_file, append = TRUE)
  
  cat(sprintf("  [REPORT] Detailed text report saved: %s\n", basename(report_file)))
  return(report_file)
}

# ============================================================================
# SECTION 3: ANALYSIS PARAMETERS
# ============================================================================

# ── Shared (preprocessing) ────────────────────────────────────────────────
MIN_COMPLETENESS_PCT <- 50   # paper Criterion B
MAX_FILL_DAYS        <- 14   # PCHIP only for gaps ≤ 14 days (FIX 2)

# ── Pipeline A: daily variable threshold — DISABLED ──────────────────────
# THRESHOLD_EXCEEDANCE_PROB <- 0.20   # Q20 = 80% exceedance (low-flow threshold)
# # Original code wrongly used 0.80 (high-flow)
# THRESHOLD_WINDOW_SIZE     <- 31     # 31-day centred moving average smoothing
# MIN_DROUGHT_DURATION_DAYS <- 30     # minimum consecutive days below threshold

# ── Pipeline B: SSI monthly (Peña-Angulo et al. 2022) ───────────────────
SSI_DROUGHT_THRESHOLD   <- if (exists("DROUGHT_ONSET")) DROUGHT_ONSET else -0.50
SSI_RECOVERY_THRESHOLD  <- if (exists("DROUGHT_END")) DROUGHT_END else -0.50  
# NOTE: "ln3" is the 3-parameter (shifted) log-normal, fitted via lmomco's
# parln3/lmomln3/cdfln3 -- NOT the standard 2-parameter log-normal. It is
# labelled "ln3" (rather than "lnorm") everywhere in this script, including
# console output, the Excel workbook, and the text report, so that it is
# never misread as the 2-parameter distribution.
DISTRIBUTIONS <- c("gev", "pe3", "ln3", "glo", "gpareto", "weibull")
MIN_MONTHS_FOR_FIT     <- 120      # overall sanity gate: minimum 10 years (120 months) of
# monthly data total before attempting any per-month fitting
MIN_YEARS_PER_MONTH    <- 10       # minimum years of THAT calendar month (e.g. # of Januaries)
# required before fitting a distribution to it; per-month
# analogue of the old (pooled) MIN_MONTHS_FOR_FIT

# ── Pipeline B: naturalized-flow homogeneity screening (08JC001) ────────
HOMOGENEITY_ALPHA <- 0.05          # significance level for the Pettitt change-point test

# ── Pipeline C: 08JC001 original (observed) data — post-regulation filter ──
# Kenney Dam / Nechako Reservoir impoundment began in 1952, after which
# 08JC001 (NECHAKO RIVER AT VANDERHOOF) reflects regulated (dam-influenced)
# flow. Pipeline C uses ONLY this post-regulation observed period for
# 08JC001 (the pre-regulation portion is dropped, and no naturalization is
# applied). Adjust this year if a more precise regulation start date is
# available for your records.
POST_REGULATION_START_YEAR <- 1952

# ============================================================================
# SECTION 4: DATA LOADING
# ============================================================================
load_discharge_data <- function(station_id, data_path = INPUT_DATA_DIR) {
  display_name <- get_station_display(station_id)
  file_patterns <- c(
    file.path(data_path, paste0(station_id, "_WSC_DISCHARGE.csv")),
    file.path(data_path, paste0(station_id, " discharge.csv")),
    file.path(data_path, paste0(station_id, ".csv")),
    file.path(data_path, paste0("WSC ", station_id, ".csv")),
    file.path(data_path, paste0(station_id, "_WaterSurveyofCanada.csv")),
    # Pipeline C: naturalized-flow CSVs (pre-regulation observed +
    # post-regulation naturalized, pre-merged), e.g. dat_natural_08JC001.csv
    file.path(data_path, paste0("dat_natural_", station_id, ".csv"))
  )
  data <- NULL
  for (file_name in file_patterns) {
    if (file.exists(file_name)) {
      cat(sprintf("  Loading: %s [%s]\n", basename(file_name), display_name))
      data <- read.csv(file_name, stringsAsFactors = FALSE)
      break
    }
  }
  if (is.null(data)) {
    warning(sprintf("  No data file found for station %s [%s]", station_id, display_name))
    return(NULL)
  }
  
  colnames(data) <- tolower(colnames(data))
  
  # Naturalized-flow exports (dat_natural_*.csv) are in long format:
  # station, date, parameter, value, symbol — with multiple parameters
  # potentially stacked in one file. Keep only the Flow rows before
  # looking for the discharge column.
  if ("parameter" %in% colnames(data)) {
    data <- data[!is.na(data$parameter) & tolower(trimws(data$parameter)) == "flow", ]
  }
  
  date_col      <- grep("date",            colnames(data), value = TRUE, ignore.case = TRUE)
  discharge_col <- grep("discharge|flow|q_", colnames(data), value = TRUE, ignore.case = TRUE)
  if (length(discharge_col) == 0) {
    # Long-format naturalized flow files store discharge in a generic
    # "value" column (paired with parameter == "Flow", filtered above)
    discharge_col <- grep("^value$", colnames(data), value = TRUE, ignore.case = TRUE)
  }
  if (length(date_col) == 0 | length(discharge_col) == 0) {
    stop(sprintf("Cannot find date/discharge columns for station %s", station_id))
  }
  
  data$date       <- as.Date(data[[date_col[1]]])
  data$discharge  <- as.numeric(data[[discharge_col[1]]])
  data            <- data[!is.na(data$date) & !is.na(data$discharge), ]
  data$year       <- year(data$date)
  data$month      <- month(data$date)
  data$day        <- day(data$date)
  data$doy        <- yday(data$date)
  data$station_id <- station_id
  
  cat(sprintf("  Records loaded: %d (%d-%d) [%s]\n",
              nrow(data), min(data$year), max(data$year), display_name))
  return(data)
}

# Pipeline B routes 08JC001 to the naturalized-flow data dir (NATURALIZED
# values). Pipeline C routes 08JC001 to the standard observed data dir
# (ORIGINAL data), which is then filtered down to the post-regulation period
# only in the Pipeline C processing loop. 08JE001 always uses INPUT_DATA_DIR
# for both pipelines.
get_data_path_for_station <- function(station_id, pipeline = c("B", "C")) {
  pipeline <- match.arg(pipeline)
  if (station_id == "08JC001" && pipeline == "B") {
    return(NATURALIZED_DATA_DIR)
  }
  return(INPUT_DATA_DIR)
}

# ============================================================================
# SECTION 5: DATA COMPLETENESS CHECK  (month-span denominator)
#
# Completeness = actual records / calendar days from:
#   - First day of the month of first observation
#   - To last day of the month of last observation
# This is less conservative than year-span but stricter than date-range.
# Large internal gaps (> MAX_FILL_DAYS) are flagged informatively.
# ============================================================================
# ============================================================================
check_data_completeness <- function(data, station_id) {
  display_name <- get_station_display(station_id)
  # stations_C shares the same StationIDs/StationNames as stations here
  # (08JC001, 08JE001), so this resolves names correctly for both Pipeline B
  # and Pipeline C stations
  name <- stations_C$StationName[stations_C$StationID == station_id]
  if (length(name) == 0) name <- "UNKNOWN"
  
  if (is.null(data) || nrow(data) == 0) {
    return(list(station_id = station_id, station_name = name, start_year = NA, end_year = NA,
                total_years = 0, total_days = 0, expected_days = 0,
                completeness = 0, meets_threshold = FALSE,
                message = "No data available"))
  }
  
  # ── Month-span denominator (FIXED) ──────────────────────────────────────────
  first_obs_date <- min(data$date)
  last_obs_date  <- max(data$date)
  
  # First day of the month of first observation
  start_date <- as.Date(paste0(format(first_obs_date, "%Y-%m"), "-01"))
  
  # Last day of the month of last observation (FIXED: numeric subtraction, not days())
  end_year  <- as.integer(format(last_obs_date, "%Y"))
  end_month <- as.integer(format(last_obs_date, "%m"))
  
  if (end_month == 12) {
    # December: go to Jan 1 of next year, then subtract 1 day
    end_date <- as.Date(paste0(end_year + 1, "-01-01")) - 1
  } else {
    # Other months: go to 1st of next month, then subtract 1 day
    end_date <- as.Date(paste0(end_year, "-", sprintf("%02d", end_month + 1), "-01")) - 1
  }
  
  expected_days <- as.numeric(difftime(end_date, start_date, units = "days")) + 1
  actual_days   <- nrow(data)
  completeness  <- actual_days / expected_days * 100
  
  # Calculate year range for reporting
  year_range  <- range(data$year)
  total_years <- diff(year_range) + 1
  
  cat(sprintf("\nStation: %s [%s]\n", station_id, display_name))
  cat(sprintf("  Period: %s to %s (%d years)\n", start_date, end_date, total_years))
  cat(sprintf("  Data: %d / %d days (%.2f%% complete)\n",
              actual_days, expected_days, completeness))
  cat(sprintf("  (Month-span: %s to %s)\n", start_date, end_date))
  
  # Flag large gaps (informational — handled in fill_missing_data)
  sorted_dates <- sort(data$date)
  gap_lengths  <- as.numeric(diff(sorted_dates)) - 1
  large_gaps   <- which(gap_lengths > MAX_FILL_DAYS)
  if (length(large_gaps) > 0) {
    cat(sprintf("  NOTE: %d gap(s) > %d days found (left as NA, not interpolated):\n",
                length(large_gaps), MAX_FILL_DAYS))
    for (gi in large_gaps) {
      cat(sprintf("    %s to %s  (%d days)\n",
                  sorted_dates[gi], sorted_dates[gi + 1L], gap_lengths[gi]))
    }
  }
  
  meets_threshold <- completeness >= MIN_COMPLETENESS_PCT
  if (!meets_threshold) {
    cat(sprintf("  STATION EXCLUDED: %.2f%% < %.0f%% threshold\n",
                completeness, MIN_COMPLETENESS_PCT))
  } else {
    cat(sprintf("  STATION PASSES: %.2f%% >= %.0f%% threshold\n",
                completeness, MIN_COMPLETENESS_PCT))
  }
  
  return(list(
    station_id      = station_id,
    station_name    = name,
    display_name    = display_name,
    start_year      = year_range[1],
    end_year        = year_range[2],
    start_date      = start_date,
    end_date        = end_date,
    total_years     = total_years,
    total_days      = actual_days,
    expected_days   = expected_days,
    completeness    = completeness,
    n_large_gaps    = length(large_gaps),
    large_gap_sizes = gap_lengths[large_gaps],
    meets_threshold = meets_threshold,
    denominator_type = "month_span"
  ))
}
# ============================================================================
# SECTION 6: MISSING DATA FILLING — Fritsch-Carlson PCHIP  [FIX 2]
#
# Replaces zoo::na.approx (linear, filled everything including a 19.7-year gap)
# with a monotone cubic Hermite spline (splinefun method="monoH.FC").
# Gaps > MAX_FILL_DAYS are left as NA. na.rm=TRUE in downstream quantile()
# and mean() calls correctly handles those remaining NAs.
# ============================================================================
fill_missing_data <- function(data, max_fill = MAX_FILL_DAYS) {
  
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  full_dates    <- seq(min(data$date), max(data$date), by = "day")
  complete_data <- merge(data.frame(date = full_dates), data, by = "date", all.x = TRUE)
  complete_data <- complete_data[order(complete_data$date), ]
  
  q <- complete_data$discharge
  n <- length(q)
  
  is_na  <- is.na(q)
  runs   <- rle(is_na)
  ends   <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  
  filled_short <- 0L
  skipped_long <- 0L
  
  for (k in seq_along(runs$lengths)) {
    if (!runs$values[k]) next
    
    run_start <- starts[k]
    run_end   <- ends[k]
    run_len   <- runs$lengths[k]
    
    if (run_len > max_fill) {
      cat(sprintf("  Gap %s to %s (%d days) > MAX_FILL_DAYS — left as NA\n",
                  complete_data$date[run_start], complete_data$date[run_end], run_len))
      skipped_long <- skipped_long + run_len
      next
    }
    
    left_idx  <- run_start - 1L
    right_idx <- run_end   + 1L
    if (left_idx < 1L || right_idx > n)          next
    if (is.na(q[left_idx]) || is.na(q[right_idx])) next
    
    # Gather local real observations within ±15 days for better spline curvature
    window_radius <- 15L
    ctx_idx <- max(1L, left_idx - window_radius):min(n, right_idx + window_radius)
    obs     <- !is.na(q[ctx_idx])
    x_obs   <- ctx_idx[obs]
    y_obs   <- q[ctx_idx[obs]]
    if (length(x_obs) < 2) next
    
    spline_fn       <- splinefun(x = x_obs, y = y_obs, method = "monoH.FC")
    q[run_start:run_end] <- pmax(0, spline_fn(run_start:run_end))
    filled_short    <- filled_short + run_len
  }
  
  complete_data$discharge  <- q
  complete_data$year       <- year(complete_data$date)
  complete_data$month      <- month(complete_data$date)
  complete_data$day        <- day(complete_data$date)
  complete_data$doy        <- yday(complete_data$date)
  if ("station_id" %in% colnames(data))
    complete_data$station_id <- data$station_id[1]
  
  cat(sprintf("  Short gaps filled (<= %d days):  %d values [Fritsch-Carlson]\n",
              max_fill, filled_short))
  cat(sprintf("  Long gaps skipped (>  %d days):  %d values [left as NA]\n",
              max_fill, skipped_long))
  cat(sprintf("  Remaining NA values: %d\n", sum(is.na(complete_data$discharge))))
  
  return(complete_data)
}

# ============================================================================
# SECTION 7: PIPELINE A — VARIABLE THRESHOLD CALCULATION — DISABLED
# ============================================================================
# calculate_variable_threshold <- function(data,
#                                          window_size    = THRESHOLD_WINDOW_SIZE,
#                                          exceedance_prob = THRESHOLD_EXCEEDANCE_PROB,
#                                          station_id     = NULL) {
#   if (is.null(data) || nrow(data) == 0) return(NULL)
#
#   thresholds <- data.frame(doy = 1:366,
#                            threshold        = NA_real_,
#                            threshold_smooth = NA_real_)
#
#   for (d in 1:366) {
#     doy_vals <- data$discharge[data$doy == d]
#     doy_vals <- doy_vals[!is.na(doy_vals)]
#
#     if (length(doy_vals) >= 10) {
#       thresholds$threshold[d] <- quantile(doy_vals, probs = exceedance_prob, na.rm = TRUE)
#     } else {
#       neighbors <- numeric(0)
#       for (offset in 1:15) {
#         d_prev    <- ((d - offset - 1) %% 366) + 1
#         d_next    <- ((d + offset - 1) %% 366) + 1
#         neighbors <- c(neighbors,
#                        data$discharge[data$doy == d_prev & !is.na(data$discharge)],
#                        data$discharge[data$doy == d_next & !is.na(data$discharge)])
#         if (length(neighbors) >= 10) break
#       }
#       if (length(neighbors) >= 2)
#         thresholds$threshold[d] <- quantile(neighbors, probs = exceedance_prob, na.rm = TRUE)
#     }
#   }
#
#   thresholds$threshold_smooth <- zoo::rollapply(
#     thresholds$threshold,
#     width = window_size, FUN = mean,
#     fill = NA, align = "center", partial = TRUE
#   )
#   thresholds$threshold_smooth[366] <- thresholds$threshold_smooth[365]
#   thresholds$station_id      <- station_id
#   thresholds$window_size     <- window_size
#   thresholds$exceedance_prob <- exceedance_prob
#
#   cat(sprintf("  [Pipeline A] Q%.0f threshold (80%%-exceedance), %d-day smoothing\n",
#               exceedance_prob * 100, window_size))
#   return(thresholds)
# }

# ============================================================================
# SECTION 8: PIPELINE A — DAILY DROUGHT EVENT IDENTIFICATION — DISABLED
# ============================================================================
# identify_droughts <- function(data, thresholds,
#                               min_duration = MIN_DROUGHT_DURATION_DAYS,
#                               station_id   = NULL) {
#   if (is.null(data) || nrow(data) == 0 || is.null(thresholds)) {
#     return(list(daily_data = data, drought_events = NULL, message = "Insufficient data"))
#   }
#
#   data     <- merge(data, thresholds[, c("doy", "threshold_smooth")],
#                     by = "doy", all.x = TRUE)
#   data     <- data[order(data$date), ]
#
#   data$below_threshold <- data$discharge < data$threshold_smooth
#   data$below_threshold[is.na(data$below_threshold)] <- FALSE   # NA days not drought
#
#   data$event_id <- with(rle(data$below_threshold), rep(seq_along(lengths), lengths))
#
#   event_info <- data %>%
#     group_by(event_id) %>%
#     filter(below_threshold == TRUE) %>%
#     summarise(
#       start_date     = min(date),
#       end_date       = max(date),
#       duration       = n(),
#       start_doy      = first(doy),
#       end_doy        = last(doy),
#       discharge_min  = min(discharge,  na.rm = TRUE),
#       discharge_mean = mean(discharge, na.rm = TRUE),
#       threshold_mean = mean(threshold_smooth, na.rm = TRUE),
#       deficit_volume = sum(threshold_smooth - discharge, na.rm = TRUE),
#       deficit_mean   = mean(threshold_smooth - discharge, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     filter(duration >= min_duration) %>%
#     mutate(
#       event_id_new = row_number(),
#       year  = year(start_date),
#       month = month(start_date),
#       season = case_when(
#         month %in% c(12, 1, 2) ~ "Winter",
#         month %in% c(3, 4, 5)  ~ "Spring",
#         month %in% c(6, 7, 8)  ~ "Summer",
#         month %in% c(9, 10,11) ~ "Fall"
#       )
#     )
#
#   if (!is.null(station_id)) event_info$station_id <- station_id
#
#   data <- data %>%
#     left_join(event_info %>% dplyr::select(event_id, event_id_new), by = "event_id")
#
#   cat(sprintf("  [Pipeline A] Drought events: %d (min %d days below Q20 threshold)\n",
#               nrow(event_info), min_duration))
#   return(list(
#     daily_data         = data,
#     drought_events     = event_info,
#     n_events           = nrow(event_info),
#     total_drought_days = sum(event_info$duration),
#     message            = "Success"
#   ))
# }

# ============================================================================
# SECTION 9: PIPELINE B — DAILY TO MONTHLY AGGREGATION
# ============================================================================
aggregate_to_monthly <- function(daily_data) {
  if (is.null(daily_data) || nrow(daily_data) == 0) return(NULL)
  
  monthly_data <- daily_data %>%
    group_by(year, month) %>%
    summarise(
      discharge_mean = mean(discharge, na.rm = TRUE),
      n_days         = sum(!is.na(discharge)),
      .groups = "drop"
    ) %>%
    filter(n_days >= 15) %>%           # require ≥15 valid daily obs per month
    mutate(
      date       = as.Date(paste0(year, "-", sprintf("%02d", month), "-15")),
      year_month = paste0(year, "-", sprintf("%02d", month))
    )
  
  cat(sprintf("  [Pipeline B] Monthly aggregation: %d months (months with >= 15 valid days)\n",
              nrow(monthly_data)))
  return(monthly_data)
}

# ============================================================================
# SECTION 10: PIPELINE B — L-MOMENTS DISTRIBUTION SELECTION & SSI
#
# Package note: the original SSI script listed 'Lmoments' but called functions
# (lmoms, lmomgev, pgev, etc.) that belong to 'lmomco'.  This script loads
# 'lmomco' instead, which provides all required functions.
# ============================================================================
calculate_lmoments_distance <- function(data, dist_name) {
  tryCatch({
    lmom_sample <- lmomco::lmoms(data, nmom = 4)
    
    if (dist_name == "gev") {
      par  <- lmomco::pargev(lmom_sample)
      lmom_theo <- lmomco::lmomgev(par)
    } else if (dist_name == "ln3") {
      par  <- lmomco::parln3(lmom_sample)
      lmom_theo <- lmomco::lmomln3(par)
    } else if (dist_name == "weibull") {
      par  <- lmomco::parwei(lmom_sample)
      lmom_theo <- lmomco::lmomwei(par)
    } else if (dist_name == "gpareto") {
      par  <- lmomco::pargpa(lmom_sample)
      lmom_theo <- lmomco::lmomgpa(par)
    } else if (dist_name == "glo") {
      par  <- lmomco::parglo(lmom_sample)
      lmom_theo <- lmomco::lmomglo(par)
    } else if (dist_name == "pe3") {
      par  <- lmomco::parpe3(lmom_sample)
      lmom_theo <- lmomco::lmompe3(par)
    } else {
      return(Inf)
    }
    
    # FIX: L1 (mean) and L2 (L-scale) match the sample by construction under
    # L-moment parameter estimation -- comparing them can never discriminate
    # between candidate distributions (this was the bug: it produced ~0 distance
    # for every distribution regardless of actual fit quality).
    #
    # The correct discriminating statistics are the higher-order L-moment
    # RATIOS: L-skewness (tau3 = lambda3/lambda2) and L-kurtosis
    # (tau4 = lambda4/lambda2). These are NOT forced to match by a fit that
    # only uses L1/L2 to solve for parameters, so they actually measure how
    # well each candidate's shape matches the sample's shape. This is the
    # standard "L-moment ratio diagram" comparison used in L-moment-based
    # distribution selection (e.g., Hosking & Wallis, 1997).
    tau3_sample <- lmom_sample$ratios[3]
    tau4_sample <- lmom_sample$ratios[4]
    tau3_theo   <- lmom_theo$ratios[3]
    tau4_theo   <- lmom_theo$ratios[4]
    
    if (any(is.na(c(tau3_sample, tau4_sample, tau3_theo, tau4_theo)))) return(Inf)
    
    return(sqrt((tau3_sample - tau3_theo)^2 + (tau4_sample - tau4_theo)^2))
    
  }, error = function(e) Inf)
}

# ============================================================================
# SECTION 10b: PIPELINE B — BOOTSTRAP DISTRIBUTION STABILITY (Fix 3)
#
# Resamples the monthly flow record B times and tallies how often each
# candidate distribution wins the L-moment distance contest.  A single-sample
# winner with < BOOT_MIN_WIN_PCT % of bootstrap draws is considered unstable;
# the plurality winner is reported alongside it in the console and saved in
# the return list so the analyst can make an informed final choice.
# ============================================================================
BOOT_REPS        <- 500L   # number of bootstrap replicates
CLOSE_MARGIN_PCT <- 15     # trigger warning + bootstrap when winner margin < this %

bootstrap_best_distribution <- function(flow_vals, distributions,
                                        B    = BOOT_REPS,
                                        seed = 42L) {
  set.seed(seed)
  win_counts <- setNames(integer(length(distributions)), distributions)
  
  for (b in seq_len(B)) {
    boot_sample <- sample(flow_vals, size = length(flow_vals), replace = TRUE)
    boot_dists  <- vapply(distributions,
                          function(d) calculate_lmoments_distance(boot_sample, d),
                          numeric(1))
    winner <- distributions[which.min(boot_dists)]
    if (length(winner) == 1L && !is.na(winner))
      win_counts[winner] <- win_counts[winner] + 1L
  }
  
  win_pct   <- win_counts / B * 100
  boot_best <- names(which.max(win_counts))
  
  cat(sprintf("  [Pipeline B] Bootstrap stability (B = %d replicates):\n", B))
  for (d in names(sort(win_counts, decreasing = TRUE))) {
    cat(sprintf("    %-10s  %3.0f%%\n", d, win_pct[d]))
  }
  cat(sprintf("  Bootstrap plurality winner: %s (%3.0f%% of draws)\n",
              boot_best, win_pct[boot_best]))
  
  return(list(win_counts = win_counts,
              win_pct    = win_pct,
              boot_best  = boot_best))
}

# ── Per-calendar-month distribution fit + CDF/SSI for ONE month's data ──────
# Factored out of calculate_ssi() so the exact same candidate-selection logic
# (L-moment distance contest, GPA disqualification, close-margin bootstrap,
# CDF fit) that used to run once on the pooled series now runs once PER
# calendar month, on that month's own flow values only.
fit_month_distribution_and_ssi <- function(flow_vals_m, target_vals_m,
                                           month_label, station_id) {
  cat(sprintf("  -- Calendar month: %-3s (n = %d years) --\n", month_label, length(flow_vals_m)))
  
  dist_distances <- data.frame(
    distribution = DISTRIBUTIONS,
    distance     = vapply(DISTRIBUTIONS,
                          function(d) calculate_lmoments_distance(flow_vals_m, d),
                          numeric(1)),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(dist_distances))) {
    cat(sprintf("      %-10s L-moments distance = %.6f\n",
                dist_distances$distribution[i], dist_distances$distance[i]))
  }
  
  best_dist <- dist_distances$distribution[which.min(dist_distances$distance)]
  cat(sprintf("    Best distribution: %s (distance = %.6f)\n",
              best_dist, min(dist_distances$distance, na.rm = TRUE)))
  
  # ── Fix 2: finite-upper-bound disqualification (per month) ─────────────────
  # GPA and GEV share the same xi/alpha/kappa parameterization in lmomco, and
  # both fall into a finite-upper-bound regime when the L-moment-estimated
  # shape parameter kappa > 0 (upper endpoint = xi + alpha/kappa). Originally
  # only GPA winners were screened for a fitted ceiling below the observed
  # maximum; GEV winners went unchecked despite being equally capable of
  # producing an inadmissible bounded fit. Both families are now screened
  # identically whenever either one wins the per-month contest.
  BOUNDED_DIST_PAR_FN <- list(gpareto = lmomco::pargpa, gev = lmomco::pargev)
  
  bounded_disqualified <- FALSE
  dist_distances_valid <- dist_distances
  if (best_dist %in% names(BOUNDED_DIST_PAR_FN)) {
    par_fit <- BOUNDED_DIST_PAR_FN[[best_dist]](lmomco::lmoms(flow_vals_m, nmom = 4))
    kappa   <- par_fit$para["kappa"]
    label   <- toupper(best_dist)
    if (!is.na(kappa) && kappa > 0) {
      upper_bound <- par_fit$para["xi"] + par_fit$para["alpha"] / kappa
      if (upper_bound < max(flow_vals_m, na.rm = TRUE)) {
        cat(sprintf(
          "    [%s] %s disqualified: fitted ceiling (%.2f m\u00b3/s) < observed max (%.2f m\u00b3/s). Promoting runner-up.\n",
          month_label, label, upper_bound, max(flow_vals_m, na.rm = TRUE)))
        dist_distances_valid <- dist_distances[dist_distances$distribution != best_dist, ]
        best_dist            <- dist_distances_valid$distribution[which.min(dist_distances_valid$distance)]
        bounded_disqualified <- TRUE
        cat(sprintf("    Runner-up selected: %s\n", best_dist))
      } else {
        cat(sprintf("    %s kappa > 0 but ceiling (%.2f) >= observed max (%.2f) — retained.\n",
                    label, upper_bound, max(flow_vals_m, na.rm = TRUE)))
      }
    } else {
      cat(sprintf("    %s kappa = %.4f (unbounded upper tail) — no disqualification needed.\n",
                  label, ifelse(is.na(kappa), NA_real_, kappa)))
    }
  }
  
  # ── Fix 1 + Fix 3: close-margin warning and conditional bootstrap (per month) ──
  dist_eligible <- dist_distances_valid[is.finite(dist_distances_valid$distance), ]
  boot_result   <- NULL
  
  if (nrow(dist_eligible) >= 2) {
    sorted_eligible <- dist_eligible[order(dist_eligible$distance), ]
    margin_pct      <- (sorted_eligible$distance[2] - sorted_eligible$distance[1]) /
      sorted_eligible$distance[1] * 100
    runner_up       <- sorted_eligible$distribution[2]
    
    if (margin_pct < CLOSE_MARGIN_PCT) {
      warning(sprintf(
        "[Pipeline B] UNSTABLE SELECTION for station %s, month %s: '%s' wins by only %.1f%% over '%s' — running bootstrap.",
        ifelse(is.null(station_id), "unknown", station_id),
        month_label, best_dist, margin_pct, runner_up))
      
      cat(sprintf("    Close margin (%.1f%% < %.0f%%) — running bootstrap (B = %d)...\n",
                  margin_pct, CLOSE_MARGIN_PCT, BOOT_REPS))
      boot_result <- tryCatch(
        bootstrap_best_distribution(flow_vals_m, sorted_eligible$distribution),
        error = function(e) {
          cat(sprintf("    Bootstrap failed: %s\n", e$message))
          NULL
        }
      )
      if (!is.null(boot_result) && boot_result$boot_best != best_dist) {
        warning(sprintf(
          "[Pipeline B] Bootstrap plurality ('%s', %.0f%%) DIFFERS from single-sample winner ('%s') for station %s month %s. Consider overriding.",
          boot_result$boot_best, boot_result$win_pct[boot_result$boot_best],
          best_dist, ifelse(is.null(station_id), "unknown", station_id), month_label))
      } else if (!is.null(boot_result)) {
        cat(sprintf("    Bootstrap confirms '%s' as plurality winner (%.0f%% of draws).\n",
                    boot_result$boot_best, boot_result$win_pct[boot_result$boot_best]))
      }
    } else {
      cat(sprintf("    Margin clear (%.1f%% >= %.0f%%) — bootstrap not needed.\n",
                  margin_pct, CLOSE_MARGIN_PCT))
    }
  }
  
  # ── CDF fit and SSI transformation, using ONLY this month's L-moments ──────
  ssi_out <- tryCatch({
    lmom_sample_m <- lmomco::lmoms(flow_vals_m, nmom = 4)
    
    cdf_fn <- switch(best_dist,
                     gev     = function(x) lmomco::cdfgev(x, lmomco::pargev(lmom_sample_m)),
                     ln3     = function(x) lmomco::cdfln3(x, lmomco::parln3(lmom_sample_m)),
                     weibull = function(x) lmomco::cdfwei(x, lmomco::parwei(lmom_sample_m)),
                     gpareto = function(x) lmomco::cdfgpa(x, lmomco::pargpa(lmom_sample_m)),
                     glo     = function(x) lmomco::cdfglo(x, lmomco::parglo(lmom_sample_m)),
                     pe3     = function(x) lmomco::cdfpe3(x, lmomco::parpe3(lmom_sample_m))
    )
    
    p_vals <- cdf_fn(target_vals_m)
    p_vals <- pmax(1e-6, pmin(1 - 1e-6, p_vals))
    qnorm(p_vals)
  }, error = function(e) {
    cat(sprintf("    SSI calculation error for %s: %s\n", month_label, e$message))
    rep(NA_real_, length(target_vals_m))
  })
  
  dist_distances$month       <- month_label
  
  list(ssi = ssi_out, best_dist = best_dist,
       dist_distances = dist_distances, boot_result = boot_result)
}

calculate_ssi <- function(monthly_data, station_id = NULL) {
  if (is.null(monthly_data) || nrow(monthly_data) < MIN_MONTHS_FOR_FIT) {
    cat(sprintf("  [Pipeline B] Insufficient monthly data for SSI (need >= %d months total)\n",
                MIN_MONTHS_FOR_FIT))
    return(NULL)
  }
  
  cat(sprintf("  [Pipeline B] Fitting distributions SEPARATELY per calendar month (L-moments criterion, need >= %d years/month)...\n",
              MIN_YEARS_PER_MONTH))
  
  monthly_data$ssi          <- NA_real_
  monthly_data$distribution <- NA_character_
  
  dist_distances_all         <- list()
  best_distribution_by_month <- setNames(rep(NA_character_, 12), month.abb)
  boot_results_by_month      <- list()
  
  for (m in 1:12) {
    month_label <- month.abb[m]
    idx_m       <- which(monthly_data$month == m)
    
    flow_all_m  <- monthly_data$discharge_mean[idx_m]
    flow_vals_m <- flow_all_m[flow_all_m > 0 & !is.na(flow_all_m)]
    
    if (length(flow_vals_m) < MIN_YEARS_PER_MONTH) {
      cat(sprintf("    [%s] Insufficient data (%d years < %d required) — SSI left as NA for this month\n",
                  month_label, length(flow_vals_m), MIN_YEARS_PER_MONTH))
      next
    }
    
    fit <- fit_month_distribution_and_ssi(
      flow_vals_m   = flow_vals_m,
      target_vals_m = monthly_data$discharge_mean[idx_m],
      month_label   = month_label,
      station_id    = station_id
    )
    
    monthly_data$ssi[idx_m]          <- fit$ssi
    monthly_data$distribution[idx_m] <- fit$best_dist
    
    best_distribution_by_month[m]        <- fit$best_dist
    dist_distances_all[[month_label]]    <- fit$dist_distances
    if (!is.null(fit$boot_result)) boot_results_by_month[[month_label]] <- fit$boot_result
  }
  
  n_ssi <- sum(!is.na(monthly_data$ssi))
  cat(sprintf("  [Pipeline B] SSI calculated for %d / %d months (fit separately per calendar month)\n",
              n_ssi, nrow(monthly_data)))
  
  if (n_ssi == 0) {
    cat("  [Pipeline B] No calendar month had sufficient data for a per-month distribution fit — SSI unavailable\n")
    return(NULL)
  }
  
  dist_distances_combined <- do.call(rbind, dist_distances_all)
  rownames(dist_distances_combined) <- NULL
  
  fitted_months <- names(best_distribution_by_month)[!is.na(best_distribution_by_month)]
  best_dist_str <- paste(sprintf("%s=%s", fitted_months, best_distribution_by_month[fitted_months]),
                         collapse = ", ")
  
  return(list(
    monthly_data          = monthly_data,
    best_distribution     = best_distribution_by_month,  # named-by-month, e.g. best_distribution["Jan"]
    best_distribution_str = best_dist_str,                # compact scalar summary for reports/Excel
    dist_distances         = dist_distances_combined,     # one block of 6 rows per fitted month, tagged by $month
    boot_result             = boot_results_by_month        # list keyed by month, only for months that triggered bootstrap
  ))
}

# ============================================================================
# SECTION 11: PIPELINE B — SSI DROUGHT EVENT IDENTIFICATION
# ============================================================================
identify_ssi_droughts <- function(ssi_result,
                                  drought_threshold  = SSI_DROUGHT_THRESHOLD,
                                  recovery_threshold = SSI_RECOVERY_THRESHOLD,
                                  station_id         = NULL) {
  if (is.null(ssi_result) || is.null(ssi_result$monthly_data)) {
    return(list(drought_events = NULL, message = "Insufficient SSI data"))
  }
  
  data <- ssi_result$monthly_data
  
  # 1. Use shared detection system for exact parity with 7basin_timeseries.R
  df_in <- data.frame(date = data$date, value = data$ssi)
  shared_events <- tryCatch(
    detect_drought_events(df_in, 
                          onset_threshold = drought_threshold,
                          termination_threshold = recovery_threshold, 
                          min_duration = 1),
    error = function(e) data.frame()
  )
  
  if (nrow(shared_events) == 0) {
    # Return empty structure expected by downstream Excel/Report generators
    event_info <- data.frame(start_date=as.Date(character()), end_date=as.Date(character()),
                             duration_months=integer(), start_year=integer(), start_month=integer(),
                             ssi_min=numeric(), ssi_mean=numeric(), severity=numeric(),
                             event_id_new=integer(), season=character(), stringsAsFactors=FALSE)
    if (!is.null(station_id)) { 
      event_info$station_id <- character()
      event_info$distribution <- character() 
    }
    data$event_id <- NA_integer_
    return(list(monthly_data=data, drought_events=event_info, n_events=0L,
                total_drought_months=0L, best_distribution=ssi_result$best_distribution,
                best_distribution_str=ssi_result$best_distribution_str, message="Success"))
  }
  
  # 2. Map events back to calculate H7WSC-specific numeric severity (deficit)
  data$event_id <- NA_integer_
  for (i in seq_len(nrow(shared_events))) {
    idx <- which(data$date >= shared_events$start_date[i] & data$date <= shared_events$end_date[i])
    data$event_id[idx] <- i
  }
  
  event_info <- data %>%
    dplyr::filter(!is.na(event_id)) %>%
    dplyr::group_by(event_id) %>%
    dplyr::summarise(
      start_date = min(date), 
      end_date = max(date), 
      duration_months = dplyr::n(),
      start_year = lubridate::year(min(date)), 
      start_month = lubridate::month(min(date)),
      ssi_min = min(ssi, na.rm = TRUE), 
      ssi_mean = mean(ssi, na.rm = TRUE),
      severity = sum(abs(ssi[ssi < drought_threshold]), na.rm = TRUE), # Kept numeric for downstream compatibility
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      event_id_new = dplyr::row_number(),
      season = dplyr::case_when(
        start_month %in% c(12,1,2) ~ "Winter", start_month %in% c(3,4,5) ~ "Spring",
        start_month %in% c(6,7,8) ~ "Summer", start_month %in% c(9,10,11) ~ "Fall"
      )
    )
  
  if (!is.null(station_id)) {
    event_info$station_id <- station_id
    # Each event is tagged with the distribution fitted for ITS OWN onset
    # calendar month (per-month fitting means this can legitimately differ
    # event to event, e.g. a winter-onset drought uses the January/February
    # fit, a summer-onset drought uses the June/July fit).
    event_info$distribution <- ssi_result$best_distribution[event_info$start_month]
  }
  
  # Re-map event_id to sequential IDs for downstream consistency
  data <- data %>% 
    dplyr::left_join(data.frame(event_id=event_info$event_id, event_id_new=event_info$event_id_new), by="event_id") %>%
    dplyr::mutate(event_id = event_id_new) %>% 
    dplyr::select(-event_id_new)
  
  cat(sprintf("  [Pipeline B] Drought events: %d (Onset < %.2f, Recovery >= %.2f)\n",
              nrow(event_info), drought_threshold, recovery_threshold))
  
  return(list(
    monthly_data           = data,
    drought_events         = event_info,
    n_events               = nrow(event_info),
    total_drought_months   = sum(event_info$duration_months),
    best_distribution      = ssi_result$best_distribution,
    best_distribution_str  = ssi_result$best_distribution_str,
    message                = "Success"
  ))
}

# SECTION 11b: PIPELINE B — NATURALIZED-FLOW HOMOGENEITY CHECK (Pettitt test)
#
# 08JC001 blends *observed* pre-regulation flow with *naturalized*
# post-regulation flow into one continuous CSV upstream of this script, with
# NO marker in the file for where that splice occurs. Feeding a series with
# an undocumented step-change into distribution fitting / SSI as if it were
# a single homogeneous record can bias both the fitted distribution and any
# drought events whose window straddles the transition.
#
# This screens the naturalized station's ANNUAL MEAN discharge series
# (annual, rather than daily/monthly, to reduce autocorrelation/noise) with
# the nonparametric Pettitt test (Pettitt, 1979) — the standard tool for
# detecting a single unknown change-point in a hydrological record, and the
# specific test named for exactly this observed/naturalized-splice situation.
#
# A significant result (p < HOMOGENEITY_ALPHA) does NOT auto-exclude the
# station — naturalization is still the best available option for using a
# regulated mainstem station at all — but it is surfaced to the analyst via
# a warning(), printed diagnostics, and carried through into the saved RDS
# objects, the Pipeline B summary CSV, and the Excel/text reports, rather
# than silently assumed away.
#
# Note: this check does NOT apply to Pipeline C's 08JC001 series, since that
# is the real observed record (no splice) simply restricted to the
# post-regulation period — see POST_REGULATION_START_YEAR.
# ============================================================================
check_naturalized_homogeneity <- function(daily_data, station_id) {
  display_name <- get_station_display(station_id)
  
  if (is.null(daily_data) || nrow(daily_data) == 0) {
    return(list(station_id = station_id, tested = FALSE, message = "No data available"))
  }
  
  annual <- daily_data %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(discharge_mean = mean(discharge, na.rm = TRUE),
                     n_days = sum(!is.na(discharge)), .groups = "drop") %>%
    dplyr::filter(n_days >= 300) %>%      # require near-complete years only
    dplyr::arrange(year)
  
  if (nrow(annual) < 10) {
    cat(sprintf("  [Pipeline B] Homogeneity check skipped for %s: insufficient near-complete annual data (%d years < 10)\n",
                display_name, nrow(annual)))
    return(list(station_id = station_id, tested = FALSE, message = "Insufficient annual data"))
  }
  if (!requireNamespace("trend", quietly = TRUE)) {
    cat(sprintf("  [Pipeline B] Homogeneity check skipped for %s: 'trend' package unavailable\n", display_name))
    return(list(station_id = station_id, tested = FALSE, message = "'trend' package unavailable"))
  }
  
  pettitt_result <- tryCatch(trend::pettitt.test(annual$discharge_mean), error = function(e) NULL)
  if (is.null(pettitt_result)) {
    cat(sprintf("  [Pipeline B] Homogeneity check FAILED for %s\n", display_name))
    return(list(station_id = station_id, tested = FALSE, message = "Pettitt test failed"))
  }
  
  change_idx  <- pettitt_result$estimate[1]
  change_year <- annual$year[change_idx]
  p_value     <- pettitt_result$p.value
  significant <- !is.na(p_value) && p_value < HOMOGENEITY_ALPHA
  
  if (change_idx >= 1 && change_idx < nrow(annual)) {
    mean_before <- mean(annual$discharge_mean[seq_len(change_idx)], na.rm = TRUE)
    mean_after  <- mean(annual$discharge_mean[(change_idx + 1):nrow(annual)], na.rm = TRUE)
    pct_shift   <- (mean_after - mean_before) / mean_before * 100
  } else {
    # Change-point estimate landed on the series boundary — can't split into
    # two non-empty segments, so report NA rather than a nonsensical shift
    mean_before <- NA_real_
    mean_after  <- NA_real_
    pct_shift   <- NA_real_
  }
  
  cat(sprintf("\n  [Pipeline B] Homogeneity check (Pettitt test) — %s:\n", display_name))
  cat(sprintf("    Annual series: %d years (%d-%d)\n", nrow(annual), min(annual$year), max(annual$year)))
  cat(sprintf("    Most likely change-point: %d  (K = %.0f, p = %.4f)\n",
              change_year, pettitt_result$statistic, p_value))
  cat(sprintf("    Mean before: %.2f m3/s | Mean after: %.2f m3/s | Shift: %+.1f%%\n",
              mean_before, mean_after, pct_shift))
  
  if (significant) {
    warning(sprintf(
      "[Pipeline B] INHOMOGENEITY DETECTED for %s: Pettitt test finds a significant step-change around %d (p = %.4f, mean shift %+.1f%%). The observed/naturalized splice for this station may not be a single homogeneous series -- interpret its distribution fit and SSI drought events with caution.",
      display_name, change_year, p_value, pct_shift))
    cat(sprintf("    -> FLAGGED: p < %.2f — series is NOT homogeneous at the %.0f%% confidence level.\n",
                HOMOGENEITY_ALPHA, (1 - HOMOGENEITY_ALPHA) * 100))
  } else {
    cat(sprintf("    -> No significant inhomogeneity detected (p >= %.2f).\n", HOMOGENEITY_ALPHA))
  }
  
  return(list(
    station_id  = station_id,
    tested      = TRUE,
    n_years     = nrow(annual),
    change_year = change_year,
    statistic   = as.numeric(pettitt_result$statistic),
    p_value     = p_value,
    significant = significant,
    mean_before = mean_before,
    mean_after  = mean_after,
    pct_shift   = pct_shift,
    alpha       = HOMOGENEITY_ALPHA,
    message     = "Success"
  ))
}

# ============================================================================
# ============================================================================
# SECTION 12: PROCESS ALL STATIONS
#
# Shared steps [1-6]: load → completeness → gap-fill → threshold/SSI → drought events
# ============================================================================

all_processed_data <- list()   # gap-filled daily data (shared)
all_completeness   <- list()   # completeness check results (shared)
all_thresholds     <- list()   # Pipeline A — DISABLED (kept as empty list for compatibility)
all_daily_droughts <- list()   # Pipeline A — DISABLED (kept as empty list for compatibility)
all_ssi_results    <- list()   # Pipeline B
all_ssi_droughts   <- list()   # Pipeline B
all_homogeneity    <- list()   # Pipeline B naturalized-flow homogeneity screening (08JC001 only)

for (i in 1:nrow(stations)) {
  station_id   <- stations$StationID[i]
  station_name <- stations$StationName[i]
  display_name <- paste0(station_id, " - ", station_name)
  
  data_path    <- get_data_path_for_station(station_id, pipeline = "B")
  
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("PROCESSING STATION %d/%d: %s\n", i, nrow(stations), display_name))
  cat(sprintf("%s\n", strrep("=", 60)))
  
  # ── [1/6] Load ─────────────────────────────────────────────────────────────
  cat(sprintf("\n[1/6] Loading discharge data (source: %s)...\n", data_path))
  raw_data <- load_discharge_data(station_id, data_path = data_path)
  if (is.null(raw_data)) {
    cat(sprintf("  SKIPPED: no data file found for %s [%s]\n", station_id, display_name))
    next
  }
  
  # Fill in StartYear/EndYear now that the data has been loaded (metadata was
  # left NA for 08JC001 since its naturalized-record span isn't hardcoded)
  stations$StartYear[stations$StationID == station_id] <- min(raw_data$year)
  stations$EndYear[stations$StationID == station_id]   <- max(raw_data$year)
  
  # ── [2/6] Completeness check ────────────────────────────────────────────────
  cat("\n[2/6] Checking data completeness...\n")
  completeness <- check_data_completeness(raw_data, station_id)
  all_completeness[[station_id]] <- completeness
  
  if (!completeness$meets_threshold) {             # FIX 3: hard exclusion gate
    cat(sprintf("  SKIPPED: %.2f%% < %.0f%% completeness threshold\n",
                completeness$completeness, MIN_COMPLETENESS_PCT))
    next
  }
  
  # ── [3/6] Gap fill (shared) ─────────────────────────────────────────────────
  cat(sprintf("\n[3/6] Filling gaps (Fritsch-Carlson PCHIP, max gap = %d days)...\n",
              MAX_FILL_DAYS))
  data_filled <- fill_missing_data(raw_data)
  all_processed_data[[station_id]] <- data_filled
  
  # ── [3b/6] Homogeneity check (naturalized-flow station 08JC001 only) ───────
  # 08JC001's naturalized CSV blends observed (pre-regulation) with
  # naturalized (post-regulation) flow with no marker of the splice point —
  # see check_naturalized_homogeneity() and SECTION 11b above.
  if (station_id == "08JC001") {
    homogeneity <- check_naturalized_homogeneity(data_filled, station_id)
    all_homogeneity[[station_id]] <- homogeneity
  }
  
  # ── [4a/6] Pipeline A: variable threshold — DISABLED ───────────────────────
  # cat("\n[4a/6] Pipeline A — calculating variable threshold (Q20, 80% exceedance)...\n")
  # thresholds <- calculate_variable_threshold(data_filled, station_id = station_id)
  # all_thresholds[[station_id]] <- thresholds
  
  # ── [5a/6] Pipeline A: daily drought identification — DISABLED ──────────────
  # cat("\n[5a/6] Pipeline A — identifying daily drought events...\n")
  # daily_droughts <- identify_droughts(data_filled, thresholds, station_id = station_id)
  # all_daily_droughts[[station_id]] <- daily_droughts
  
  # ── [4b/6] Pipeline B: aggregate to monthly ─────────────────────────────────
  cat("\n[4b/6] Pipeline B — aggregating to monthly...\n")
  monthly_data <- aggregate_to_monthly(data_filled)
  
  # ── [5b/6] Pipeline B: SSI calculation ─────────────────────────────────────
  cat("\n[5b/6] Pipeline B — fitting distributions & calculating SSI...\n")
  ssi_result <- calculate_ssi(monthly_data, station_id = station_id)
  if (!is.null(ssi_result)) {
    all_ssi_results[[station_id]] <- ssi_result
    
    # ── [6b/6] Pipeline B: SSI drought identification ─────────────────────────
    cat("\n[6b/6] Pipeline B — identifying SSI drought events (SSI < -0.50)...\n")
    ssi_droughts <- identify_ssi_droughts(ssi_result, station_id = station_id)
    all_ssi_droughts[[station_id]] <- ssi_droughts
  } else {
    cat("  [Pipeline B] SSI skipped — insufficient monthly data\n")
  }
  
  cat(sprintf("\n  Station %s complete.\n", station_id))
}

# ============================================================================
# SECTION 12b: PROCESS ALL STATIONS — PIPELINE C
#
# Pipeline C is methodologically identical to Pipeline B (same SSI approach,
# same shared preprocessing steps 1-6), but runs on the ORIGINAL (observed,
# non-naturalized) 08JC001 record instead of the naturalized one:
#   08JC001 — NECHAKO RIVER AT VANDERHOOF — original/observed flow, loaded
#             from the standard INPUT_DATA_DIR, then restricted to the
#             POST-REGULATION period only (year >= POST_REGULATION_START_YEAR;
#             the pre-regulation years are dropped and no naturalization is
#             applied)
#   08JE001 — STUART RIVER NEAR FORT ST. JAMES — observed flow (same series
#             as Pipeline B)
# Since 08JC001 here is a real observed record (not a pre-/post-regulation
# splice), no homogeneity screening is needed — see SECTION 11b, which
# applies only to Pipeline B's naturalized 08JC001 series.
# ============================================================================

all_processed_data_C <- list()   # gap-filled daily data (Pipeline C)
all_completeness_C   <- list()   # completeness check results (Pipeline C)
all_ssi_results_C    <- list()   # Pipeline C SSI results
all_ssi_droughts_C   <- list()   # Pipeline C SSI drought events

for (i in 1:nrow(stations_C)) {
  station_id   <- stations_C$StationID[i]
  station_name <- stations_C$StationName[i]
  display_name <- paste0(station_id, " - ", station_name)
  data_path    <- get_data_path_for_station(station_id, pipeline = "C")
  
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("[PIPELINE C] PROCESSING STATION %d/%d: %s\n", i, nrow(stations_C), display_name))
  cat(sprintf("%s\n", strrep("=", 60)))
  
  # ── [1/6] Load ─────────────────────────────────────────────────────────────
  cat(sprintf("\n[1/6] Loading discharge data (source: %s)...\n", data_path))
  raw_data <- load_discharge_data(station_id, data_path = data_path)
  if (is.null(raw_data)) {
    cat(sprintf("  SKIPPED: no data file found for %s [%s]\n", station_id, display_name))
    next
  }
  
  # Pipeline C: 08JC001 uses the ORIGINAL (observed) record — restrict to the
  # post-regulation period only (drop pre-regulation years; no naturalization
  # is applied here, unlike Pipeline B's 08JC001)
  if (station_id == "08JC001") {
    n_before <- nrow(raw_data)
    raw_data <- raw_data[raw_data$year >= POST_REGULATION_START_YEAR, ]
    cat(sprintf("  [Pipeline C] 08JC001: restricted to post-regulation period (year >= %d): %d -> %d records\n",
                POST_REGULATION_START_YEAR, n_before, nrow(raw_data)))
    if (nrow(raw_data) == 0) {
      cat(sprintf("  SKIPPED: no post-regulation records remain for %s [%s]\n", station_id, display_name))
      next
    }
  }
  
  # Fill in StartYear/EndYear now that the data has been loaded (metadata was
  # left NA for 08JC001 since its post-regulation-only span isn't hardcoded)
  stations_C$StartYear[stations_C$StationID == station_id] <- min(raw_data$year)
  stations_C$EndYear[stations_C$StationID == station_id]   <- max(raw_data$year)
  
  # ── [2/6] Completeness check ────────────────────────────────────────────────
  cat("\n[2/6] Checking data completeness...\n")
  completeness <- check_data_completeness(raw_data, station_id)
  all_completeness_C[[station_id]] <- completeness
  
  if (!completeness$meets_threshold) {             # hard exclusion gate
    cat(sprintf("  SKIPPED: %.2f%% < %.0f%% completeness threshold\n",
                completeness$completeness, MIN_COMPLETENESS_PCT))
    next
  }
  
  # ── [3/6] Gap fill (shared) ─────────────────────────────────────────────────
  cat(sprintf("\n[3/6] Filling gaps (Fritsch-Carlson PCHIP, max gap = %d days)...\n",
              MAX_FILL_DAYS))
  data_filled <- fill_missing_data(raw_data)
  all_processed_data_C[[station_id]] <- data_filled
  
  # ── [4/6] Aggregate to monthly ─────────────────────────────────────────────
  cat("\n[4/6] Pipeline C — aggregating to monthly...\n")
  monthly_data <- aggregate_to_monthly(data_filled)
  
  # ── [5/6] SSI calculation ─────────────────────────────────────────────────
  cat("\n[5/6] Pipeline C — fitting distributions & calculating SSI...\n")
  ssi_result <- calculate_ssi(monthly_data, station_id = station_id)
  if (!is.null(ssi_result)) {
    all_ssi_results_C[[station_id]] <- ssi_result
    
    # ── [6/6] SSI drought identification ─────────────────────────────────────
    cat("\n[6/6] Pipeline C — identifying SSI drought events (SSI < -0.50)...\n")
    ssi_droughts <- identify_ssi_droughts(ssi_result, station_id = station_id)
    all_ssi_droughts_C[[station_id]] <- ssi_droughts
  } else {
    cat("  [Pipeline C] SSI skipped — insufficient monthly data\n")
  }
  
  cat(sprintf("\n  [Pipeline C] Station %s complete.\n", station_id))
}

# ============================================================================
# SECTION 13: SAVE RESULTS
# ============================================================================
cat(sprintf("\n%s\nSAVING RESULTS...\n%s\n", strrep("=", 60), strrep("=", 60)))

# Shared
saveRDS(stations,         file.path(MAIN_OUTPUT_DIR, "shared", "stations_metadata.rds"))
saveRDS(all_completeness, file.path(MAIN_OUTPUT_DIR, "shared", "completeness_reports.rds"))
for (sid in names(all_processed_data)) {
  saveRDS(all_processed_data[[sid]],
          file.path(MAIN_OUTPUT_DIR, "shared", "processed_data",
                    paste0(sid, "_daily_filled.rds")))
}
cat("  Shared outputs saved\n")

# Pipeline A — DISABLED
# saveRDS(all_thresholds,
#         file.path(MAIN_OUTPUT_DIR, "daily", "thresholds", "thresholds_all.rds"))
# saveRDS(all_daily_droughts,
#         file.path(MAIN_OUTPUT_DIR, "daily", "drought_events", "daily_droughts_all.rds"))
# cat("  Pipeline A (daily threshold) outputs saved\n")

# Pipeline B (08JC001 naturalized + 08JE001 observed)
for (sid in names(all_ssi_results)) {
  saveRDS(all_ssi_results[[sid]]$monthly_data,
          file.path(MAIN_OUTPUT_DIR, "ssi", "monthly_data",
                    paste0(sid, "_monthly_ssi.rds")))
}
saveRDS(all_ssi_results,
        file.path(MAIN_OUTPUT_DIR, "ssi", "drought_events", "ssi_results_all.rds"))
saveRDS(all_ssi_droughts,
        file.path(MAIN_OUTPUT_DIR, "ssi", "drought_events", "ssi_droughts_all.rds"))
saveRDS(all_homogeneity,
        file.path(MAIN_OUTPUT_DIR, "ssi", "naturalized_homogeneity_check.rds"))
cat("  Pipeline B (SSI monthly) outputs saved\n")
cat("  Pipeline B naturalized-flow homogeneity check saved\n")

# Pipeline C (08JC001 original/observed, post-regulation only + 08JE001 observed)
for (sid in names(all_processed_data_C)) {
  saveRDS(all_processed_data_C[[sid]],
          file.path(MAIN_OUTPUT_DIR, "shared", "processed_data",
                    paste0(sid, "_pipelineC_daily_filled.rds")))
}
for (sid in names(all_ssi_results_C)) {
  saveRDS(all_ssi_results_C[[sid]]$monthly_data,
          file.path(MAIN_OUTPUT_DIR, "ssi_C", "monthly_data",
                    paste0(sid, "_monthly_ssi.rds")))
}
saveRDS(all_ssi_results_C,
        file.path(MAIN_OUTPUT_DIR, "ssi_C", "drought_events", "ssi_results_all.rds"))
saveRDS(all_ssi_droughts_C,
        file.path(MAIN_OUTPUT_DIR, "ssi_C", "drought_events", "ssi_droughts_all.rds"))
saveRDS(stations_C,
        file.path(MAIN_OUTPUT_DIR, "shared", "stations_metadata_pipelineC.rds"))
cat("  Pipeline C (SSI monthly, original/post-regulation network) outputs saved\n")

all_results_C <- list(
  stations           = stations_C,
  completeness       = all_completeness_C,
  processed_data     = all_processed_data_C,
  ssi_results        = all_ssi_results_C,
  ssi_droughts       = all_ssi_droughts_C,
  metadata = list(
    input_data_dir            = INPUT_DATA_DIR,
    output_dir                = MAIN_OUTPUT_DIR,
    completeness_denominator  = "month_span (first_day_of_first_month to last_day_of_last_month)",
    min_completeness_pct      = MIN_COMPLETENESS_PCT,
    max_fill_days             = MAX_FILL_DAYS,
    interpolation_method      = "monoH.FC (Fritsch-Carlson PCHIP)",
    pipeline_c                = "Peña-Angulo et al. (2022) SSI monthly, on 08JC001 ORIGINAL (observed) flow + 08JE001",
    post_regulation_stations  = "08JC001",
    post_regulation_note      = sprintf(
      "08JC001: original/observed WSC discharge record (NOT naturalized), restricted to year >= %d (POST_REGULATION_START_YEAR); pre-regulation years dropped",
      POST_REGULATION_START_YEAR),
    post_regulation_start_year = POST_REGULATION_START_YEAR,
    ssi_drought_threshold     = SSI_DROUGHT_THRESHOLD,
    ssi_recovery_threshold    = SSI_RECOVERY_THRESHOLD,
    distributions_tested      = DISTRIBUTIONS,
    distribution_selection    = "minimum L-moments Euclidean distance (lmomco)",
    min_months_for_fit        = MIN_MONTHS_FOR_FIT,
    min_drought_duration_months = 1L,
    processed_date            = Sys.Date(),
    script_version             = "5.0_pipelineB_naturalized_pipelineC_original_postreg"
  )
)
saveRDS(all_results_C, file.path(MAIN_OUTPUT_DIR, "all_results_pipeline_C.rds"))
cat("  Master file saved: all_results_pipeline_C.rds\n")

# Master combined file
all_results <- list(
  stations           = stations,
  completeness       = all_completeness,
  processed_data     = all_processed_data,
  # Pipeline A — DISABLED
  # thresholds         = all_thresholds,
  # daily_droughts     = all_daily_droughts,
  # Pipeline B
  ssi_results        = all_ssi_results,
  ssi_droughts       = all_ssi_droughts,
  homogeneity_check  = all_homogeneity,
  metadata = list(
    input_data_dir            = INPUT_DATA_DIR,
    naturalized_data_dir      = NATURALIZED_DATA_DIR,
    output_dir                = MAIN_OUTPUT_DIR,
    # Shared preprocessing
    completeness_denominator  = "month_span (first_day_of_first_month to last_day_of_last_month)",
    min_completeness_pct      = MIN_COMPLETENESS_PCT,
    max_fill_days             = MAX_FILL_DAYS,
    interpolation_method      = "monoH.FC (Fritsch-Carlson PCHIP)",
    # Pipeline A — DISABLED
    # pipeline_a                = "Raut & Ganguli (2024) variable threshold",
    # threshold_exceedance_prob = 1 - THRESHOLD_EXCEEDANCE_PROB,
    # threshold_quantile_prob   = THRESHOLD_EXCEEDANCE_PROB,
    # threshold_window_days     = THRESHOLD_WINDOW_SIZE,
    # min_drought_duration_days = MIN_DROUGHT_DURATION_DAYS,
    # Pipeline B
    pipeline_b                = "Peña-Angulo et al. (2022) SSI monthly, on 08JC001 NATURALIZED flow + 08JE001",
    naturalized_stations      = "08JC001",
    naturalized_note          = "08JC001: observed discharge (pre-regulation) blended with naturalized flow (post-regulation), pre-merged into a single CSV",
    homogeneity_test          = "Pettitt test (trend::pettitt.test) on annual-mean discharge, alpha = HOMOGENEITY_ALPHA; see all_homogeneity / homogeneity_check for per-station results",
    homogeneity_alpha         = HOMOGENEITY_ALPHA,
    ssi_drought_threshold     = SSI_DROUGHT_THRESHOLD,
    ssi_recovery_threshold    = SSI_RECOVERY_THRESHOLD,
    distributions_tested      = DISTRIBUTIONS,
    distribution_selection    = "minimum L-moments Euclidean distance (lmomco)",
    min_months_for_fit        = MIN_MONTHS_FOR_FIT,
    min_drought_duration_months = 1L,
    processed_date            = Sys.Date(),
    script_version            = "5.0_pipelineB_naturalized_pipelineC_original_postreg"
  )
)
saveRDS(all_results, file.path(MAIN_OUTPUT_DIR, "all_results_dual_pipeline.rds"))
cat("  Master file saved: all_results_dual_pipeline.rds\n")

# ============================================================================
# SECTION 13.5: ADDITIONAL OUTPUTS (from SWEI code)
# ============================================================================
cat(sprintf("\n%s\nGENERATING ADDITIONAL OUTPUTS...\n%s\n", 
            strrep("=", 60), strrep("=", 60)))

# Station completeness plot
create_station_completeness_plot(stations, all_completeness, MAIN_OUTPUT_DIR)

# Excel summary workbook
create_excel_summary(stations, all_completeness,
                     all_daily_droughts, all_ssi_droughts,
                     all_ssi_results, MAIN_OUTPUT_DIR)

# Detailed text report
params <- list(
  script_version = "4.0_merged_dual_pipeline",
  min_completeness_pct = MIN_COMPLETENESS_PCT
)
create_text_report(stations, all_completeness,
                   all_daily_droughts, all_ssi_droughts,
                   all_ssi_results, MAIN_OUTPUT_DIR, params)

# ============================================================================
# SECTION 14: SUMMARY REPORT
# ============================================================================
cat(sprintf("\n%s\nDUAL-PIPELINE SUMMARY\n%s\n\n", strrep("=", 60), strrep("=", 60)))

summary_rows <- list()
for (sid in names(all_completeness)) {
  comp <- all_completeness[[sid]]
  sname <- if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN"
  base <- data.frame(
    StationID        = sid,
    StationName      = sname,
    StationDisplay   = paste0(sid, " - ", sname),
    Period           = sprintf("%d-%d", comp$start_year, comp$end_year),
    Completeness_pct = round(comp$completeness, 1),
    Status           = ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED"),
    DataSource       = ifelse(sid == "08JC001",
                              "Naturalized (pre-reg observed + post-reg naturalized)",
                              "Observed"),
    stringsAsFactors = FALSE
  )
  if (!comp$meets_threshold) {
    # Pipeline A stats — DISABLED
    # base$A_Drought_Events      <- NA_integer_
    # base$A_Total_Days          <- NA_integer_
    # base$A_Avg_Duration_Days   <- NA_real_
    base$B_Best_Distribution   <- NA_character_
    base$B_Drought_Events      <- NA_integer_
    base$B_Total_Months        <- NA_integer_
    base$B_Avg_Duration_Months <- NA_real_
    base$B_Avg_Severity        <- NA_real_
  } else {
    # Pipeline A stats — DISABLED
    # ev_a <- all_daily_droughts[[sid]]$drought_events
    # base$A_Drought_Events     <- if (!is.null(ev_a)) nrow(ev_a)           else NA_integer_
    # base$A_Total_Days         <- if (!is.null(ev_a)) sum(ev_a$duration)   else NA_integer_
    # base$A_Avg_Duration_Days  <- if (!is.null(ev_a) && nrow(ev_a) > 0)
    #   round(mean(ev_a$duration), 1)           else NA_real_
    
    # Pipeline B stats
    ev_b <- all_ssi_droughts[[sid]]$drought_events
    base$B_Best_Distribution    <- if (!is.null(all_ssi_results[[sid]]))
      all_ssi_results[[sid]]$best_distribution_str else NA_character_
    base$B_Drought_Events       <- if (!is.null(ev_b)) nrow(ev_b)                     else NA_integer_
    base$B_Total_Months         <- if (!is.null(ev_b)) sum(ev_b$duration_months)      else NA_integer_
    base$B_Avg_Duration_Months  <- if (!is.null(ev_b) && nrow(ev_b) > 0)
      round(mean(ev_b$duration_months), 1)             else NA_real_
    base$B_Avg_Severity         <- if (!is.null(ev_b) && nrow(ev_b) > 0)
      round(mean(ev_b$severity), 2)                   else NA_real_
  }
  
  # Naturalized-flow homogeneity screening (08JC001 only; NA for observed-only
  # stations since the Pettitt splice check doesn't apply)
  hg <- all_homogeneity[[sid]]
  if (!is.null(hg) && isTRUE(hg$tested)) {
    base$Homogeneity_ChangeYear <- hg$change_year
    base$Homogeneity_PValue     <- round(hg$p_value, 4)
    base$Homogeneity_Flag       <- ifelse(hg$significant, "INHOMOGENEOUS", "homogeneous")
  } else {
    base$Homogeneity_ChangeYear <- NA_integer_
    base$Homogeneity_PValue     <- NA_real_
    base$Homogeneity_Flag       <- NA_character_
  }
  
  summary_rows[[sid]] <- base
}

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table,
          file.path(MAIN_OUTPUT_DIR, "reports", "station_summary_dual_pipeline.csv"),
          row.names = FALSE)
print(summary_table, row.names = FALSE)

# ============================================================================
# SECTION 14b: PIPELINE C SUMMARY (08JC001 original/post-regulation + 08JE001)
# ============================================================================
cat(sprintf("\n%s\nPIPELINE C SUMMARY (08JC001 original/observed, post-regulation only + 08JE001)\n%s\n\n",
            strrep("=", 60), strrep("=", 60)))

summary_rows_C <- list()
for (sid in names(all_completeness_C)) {
  comp <- all_completeness_C[[sid]]
  sname <- if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN"
  base <- data.frame(
    StationID        = sid,
    StationName      = sname,
    StationDisplay   = paste0(sid, " - ", sname),
    Period           = sprintf("%d-%d", comp$start_year, comp$end_year),
    Completeness_pct = round(comp$completeness, 1),
    Status           = ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED"),
    DataSource       = ifelse(sid == "08JC001",
                              sprintf("Original/observed, post-regulation only (year >= %d)", POST_REGULATION_START_YEAR),
                              "Observed"),
    stringsAsFactors = FALSE
  )
  if (!comp$meets_threshold) {
    base$Best_Distribution   <- NA_character_
    base$Drought_Events      <- NA_integer_
    base$Total_Months        <- NA_integer_
    base$Avg_Duration_Months <- NA_real_
    base$Avg_Severity        <- NA_real_
  } else {
    ev_c <- all_ssi_droughts_C[[sid]]$drought_events
    base$Best_Distribution   <- if (!is.null(all_ssi_results_C[[sid]]))
      all_ssi_results_C[[sid]]$best_distribution_str else NA_character_
    base$Drought_Events      <- if (!is.null(ev_c)) nrow(ev_c)                else NA_integer_
    base$Total_Months        <- if (!is.null(ev_c)) sum(ev_c$duration_months) else NA_integer_
    base$Avg_Duration_Months <- if (!is.null(ev_c) && nrow(ev_c) > 0)
      round(mean(ev_c$duration_months), 1) else NA_real_
    base$Avg_Severity        <- if (!is.null(ev_c) && nrow(ev_c) > 0)
      round(mean(ev_c$severity), 2) else NA_real_
  }
  
  summary_rows_C[[sid]] <- base
}

summary_table_C <- do.call(rbind, summary_rows_C)
write.csv(summary_table_C,
          file.path(MAIN_OUTPUT_DIR, "reports", "station_summary_pipeline_C.csv"),
          row.names = FALSE)
print(summary_table_C, row.names = FALSE)

cat("\nOutputs written to: ", MAIN_OUTPUT_DIR, "/\n")
cat("  shared/processed_data/            — gap-filled daily data (Pipeline B + Pipeline C)\n")
# cat("  daily/thresholds/                 — Pipeline A: Q20 variable thresholds\n")   # DISABLED
# cat("  daily/drought_events/             — Pipeline A: daily drought events\n")       # DISABLED
cat("  ssi/monthly_data/                 — Pipeline B: monthly SSI time series (08JC001 naturalized, 08JE001)\n")
cat("  ssi/drought_events/               — Pipeline B: SSI drought events\n")
cat("  ssi/naturalized_homogeneity_check.rds — Pipeline B: Pettitt-test homogeneity screening (08JC001 naturalized)\n")
cat("  ssi_C/monthly_data/               — Pipeline C: monthly SSI time series (08JC001 original/post-regulation, 08JE001)\n")
cat("  ssi_C/drought_events/             — Pipeline C: SSI drought events\n")
cat("  figures/                          — station_completeness_barplot.png\n")
cat("  reports/                          — station_summary_dual_pipeline.csv\n")
cat("                                    — station_summary_pipeline_C.csv\n")
cat("                                    — nechako_streamflow_drought_summary.xlsx\n")
cat("                                    — nechako_streamflow_drought_report.txt\n")
cat("  all_results_dual_pipeline.rds     — Pipeline B master combined file (2 stations)\n")
cat("  all_results_pipeline_C.rds        — Pipeline C master combined file (2 stations)\n")