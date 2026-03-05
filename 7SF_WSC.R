# ============================================================================
# NECHAKO STREAMFLOW RIVER BASIN — DUAL-PIPELINE DROUGHT ANALYSIS
# 
# Pipeline A (Daily):   Variable threshold drought, Raut & Ganguli (2024)
# Pipeline B (Monthly): Standardized Streamflow Index, Peña-Angulo et al. (2022)
#
# Shared preprocessing (both pipelines):
#   - Data loading & completeness check  (original year-span logic)
#   - Gap filling: Fritsch-Carlson PCHIP, gaps > MAX_FILL_DAYS left as NA  
#   - Hard exclusion gate for stations < 70% completeness                  
#
# Pipeline A adds:
#   - Per-DOY variable threshold (Q20 = 80% exceedance probability)        
#   - Daily drought event identification (min 30 days)
#
# Pipeline B adds:
#   - Daily → monthly mean aggregation
#   - Best-fit distribution selection via L-moments (lmomco)
#   - SSI calculation and drought identification (SSI < -0.84)
#
# check_data_completeness() uses the original month-span denominator.
#       Only 08JE001 and 08JE004 are expected to pass (≥70% completeness).
# ============================================================================

rm(list = ls())

# ============================================================================
# SECTION 1: PACKAGES
# ============================================================================
packages_needed <- c(
  # Shared / general
  "tidyverse", "lubridate", "zoo",
  # Pipeline A (daily threshold)
  "circular", "MASS", "trend", "Kendall", "boot", "copula", "kde1d",
  # Pipeline B (SSI)
  "lmomco", "fitdistrplus",
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
stations <- data.frame(
  StationID   = c("08JA015", "08JB002", "08JC001", "08JC002",
                  "08JC005", "08JE001", "08JE004"),
  StationName = c("LAVENTIE CREEK NEAR THE MOUTH",
                  "STELLAKO RIVER AT GLENANNAN",
                  "NECHAKO RIVER AT VANDERHOOF",
                  "NECHAKO RIVER AT ISLE PIERRE",
                  "CHILAKO RIVER NEAR PRINCE GEORGE",
                  "STUART RIVER NEAR FORT ST. JAMES",
                  "TSILCOH RIVER NEAR THE MOUTH"),
  StartYear   = c(1976, 1929, 1915, 1950, 1953, 1929, 1975),
  EndYear     = c(2022, 2024, 2023, 2023, 2022, 2024, 2024),
  stringsAsFactors = FALSE
)

setwd("D:/Nechako_Drought/Nechako")

INPUT_DATA_DIR  <- paste0("D:/Nechako_Drought/Nechako/Hydrology/",
                          "data_retrievalWaterSurveyofCanada/data_downloads_geomet_api")

# Single unified output root — both pipelines write under their own subdirs
MAIN_OUTPUT_DIR <- "streamflow_results"

# STATION DISPLAY HELPER
get_station_display <- function(station_id, stations_df = stations) {
  name <- stations_df$StationName[stations_df$StationID == station_id]
  if (length(name) == 0) name <- "UNKNOWN"
  return(paste0(station_id, " - ", name))
}

for (sub in c(
  "",
  "shared/processed_data",          # daily gap-filled data (used by both)
  "shared/completeness",
  "daily/thresholds",               # Pipeline A
  "daily/drought_events",
  "ssi/monthly_data",               # Pipeline B
  "ssi/drought_events",
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
  abline(h = 70, lty = 2, lwd = 2, col = "red")
  legend("topright", 
         legend = c("INCLUDED (≥70%)", "EXCLUDED (<70%)"),
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
    
    # Pipeline A stats - ALWAYS add these columns (NA if excluded)
    if (comp$meets_threshold && !is.null(daily_droughts[[sid]]$drought_events)) {
      ev_a <- daily_droughts[[sid]]$drought_events
      base$A_Drought_Events      <- nrow(ev_a)
      base$A_Total_Days          <- sum(ev_a$duration)
      base$A_Avg_Duration_Days   <- if (nrow(ev_a) > 0) round(mean(ev_a$duration), 1) else NA_real_
    } else {
      base$A_Drought_Events      <- NA_integer_
      base$A_Total_Days          <- NA_integer_
      base$A_Avg_Duration_Days   <- NA_real_
    }
    
    # Pipeline B stats - ALWAYS add these columns (NA if excluded)
    if (comp$meets_threshold && !is.null(ssi_results[[sid]])) {
      base$B_Best_Distribution <- ssi_results[[sid]]$best_distribution
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
  
  # Sheets for Pipeline A Drought Events (one per included station)
  for (sid in names(daily_droughts)) {
    if (!is.null(daily_droughts[[sid]]$drought_events)) {
      ev <- daily_droughts[[sid]]$drought_events
      if (nrow(ev) > 0) {
        excel_data[[paste0("A_", sid, "_Daily_Droughts")]] <- as.data.frame(ev)
      }
    }
  }
  
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
      if (!is.null(daily_droughts[[sid]]$drought_events)) {
        ev_a <- daily_droughts[[sid]]$drought_events
        cat(sprintf("\nPipeline A (Daily Threshold):\n"),
            file = report_file, append = TRUE)
        cat(sprintf("  Drought events: %d\n", nrow(ev_a)),
            file = report_file, append = TRUE)
        cat(sprintf("  Total drought days: %d\n", sum(ev_a$duration)),
            file = report_file, append = TRUE)
        cat(sprintf("  Average duration: %.1f days\n", mean(ev_a$duration)),
            file = report_file, append = TRUE)
      }
      
      if (!is.null(ssi_droughts[[sid]]$drought_events)) {
        ev_b <- ssi_droughts[[sid]]$drought_events
        cat(sprintf("\nPipeline B (Monthly SSI):\n"),
            file = report_file, append = TRUE)
        cat(sprintf("  Best distribution: %s\n", 
                    ssi_results[[sid]]$best_distribution),
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
  cat("Pipeline A (Daily):\n",
      "  • Variable threshold (Q20 = 80% exceedance probability)\n",
      "  • 31-day centered moving average smoothing\n",
      "  • Minimum drought duration: 30 days\n",
      "  • Reference: Raut & Ganguli (2024)\n\n",
      file = report_file, append = TRUE)
  cat("Pipeline B (Monthly):\n",
      "  • Standardized Streamflow Index (SSI)\n",
      "  • Best-fit distribution via L-moments (lmomco)\n",
      "  • Drought threshold: SSI < -0.84 (5-year return period)\n",
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
MIN_COMPLETENESS_PCT <- 70   # paper Criterion B
MAX_FILL_DAYS        <- 14   # PCHIP only for gaps ≤ 14 days (FIX 2)

# ── Pipeline A: daily variable threshold (Raut & Ganguli 2024) ───────────
THRESHOLD_EXCEEDANCE_PROB <- 0.20   # Q20 = 80% exceedance (low-flow threshold)
# Original code wrongly used 0.80 (high-flow)
THRESHOLD_WINDOW_SIZE     <- 31     # 31-day centred moving average smoothing
MIN_DROUGHT_DURATION_DAYS <- 30     # minimum consecutive days below threshold

# ── Pipeline B: SSI monthly (Peña-Angulo et al. 2022) ───────────────────
SSI_DROUGHT_THRESHOLD  <- -0.84    # 5-year return period (Peña-Angulo §2.2.1)
SSI_RECOVERY_THRESHOLD <-  0.00    # drought ends when SSI returns to 0
DISTRIBUTIONS          <- c("gev", "pe3", "lnorm", "llogis", "gpareto", "weibull")
MIN_MONTHS_FOR_FIT     <- 120      # minimum 10 years of monthly data

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
    file.path(data_path, paste0(station_id, "_WaterSurveyofCanada.csv"))
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
  date_col      <- grep("date",            colnames(data), value = TRUE, ignore.case = TRUE)
  discharge_col <- grep("discharge|flow|q_", colnames(data), value = TRUE, ignore.case = TRUE)
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
  name <- stations$StationName[stations$StationID == station_id]
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
# SECTION 7: PIPELINE A — VARIABLE THRESHOLD CALCULATION
#
# Computes the Q20 (20th percentile = 80% exceedance probability) for each
# day-of-year across all years, then smooths with a 31-day centred moving
# average.  Original code used probs=0.80 (80th percentile, a high-flow value)
# which caused ~80% of all days to be classified as drought.
# ============================================================================
calculate_variable_threshold <- function(data,
                                         window_size    = THRESHOLD_WINDOW_SIZE,
                                         exceedance_prob = THRESHOLD_EXCEEDANCE_PROB,
                                         station_id     = NULL) {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  thresholds <- data.frame(doy = 1:366,
                           threshold        = NA_real_,
                           threshold_smooth = NA_real_)
  
  for (d in 1:366) {
    doy_vals <- data$discharge[data$doy == d]
    doy_vals <- doy_vals[!is.na(doy_vals)]
    
    if (length(doy_vals) >= 10) {
      thresholds$threshold[d] <- quantile(doy_vals, probs = exceedance_prob, na.rm = TRUE)
    } else {
      neighbors <- numeric(0)
      for (offset in 1:15) {
        d_prev    <- ((d - offset - 1) %% 366) + 1
        d_next    <- ((d + offset - 1) %% 366) + 1
        neighbors <- c(neighbors,
                       data$discharge[data$doy == d_prev & !is.na(data$discharge)],
                       data$discharge[data$doy == d_next & !is.na(data$discharge)])
        if (length(neighbors) >= 10) break
      }
      if (length(neighbors) >= 2)
        thresholds$threshold[d] <- quantile(neighbors, probs = exceedance_prob, na.rm = TRUE)
    }
  }
  
  thresholds$threshold_smooth <- zoo::rollapply(
    thresholds$threshold,
    width = window_size, FUN = mean,
    fill = NA, align = "center", partial = TRUE
  )
  thresholds$threshold_smooth[366] <- thresholds$threshold_smooth[365]
  thresholds$station_id      <- station_id
  thresholds$window_size     <- window_size
  thresholds$exceedance_prob <- exceedance_prob
  
  cat(sprintf("  [Pipeline A] Q%.0f threshold (80%%-exceedance), %d-day smoothing\n",
              exceedance_prob * 100, window_size))
  return(thresholds)
}

# ============================================================================
# SECTION 8: PIPELINE A — DAILY DROUGHT EVENT IDENTIFICATION
# ============================================================================
identify_droughts <- function(data, thresholds,
                              min_duration = MIN_DROUGHT_DURATION_DAYS,
                              station_id   = NULL) {
  if (is.null(data) || nrow(data) == 0 || is.null(thresholds)) {
    return(list(daily_data = data, drought_events = NULL, message = "Insufficient data"))
  }
  
  data$doy <- ifelse(leap_year(data$year) & data$doy > 366, 366, data$doy)
  data     <- merge(data, thresholds[, c("doy", "threshold_smooth")],
                    by = "doy", all.x = TRUE)
  data     <- data[order(data$date), ]
  
  data$below_threshold <- data$discharge < data$threshold_smooth
  data$below_threshold[is.na(data$below_threshold)] <- FALSE   # NA days not drought
  
  data$event_id <- with(rle(data$below_threshold), rep(seq_along(lengths), lengths))
  
  event_info <- data %>%
    group_by(event_id) %>%
    filter(below_threshold == TRUE) %>%
    summarise(
      start_date     = min(date),
      end_date       = max(date),
      duration       = n(),
      start_doy      = first(doy),
      end_doy        = last(doy),
      discharge_min  = min(discharge,  na.rm = TRUE),
      discharge_mean = mean(discharge, na.rm = TRUE),
      threshold_mean = mean(threshold_smooth, na.rm = TRUE),
      deficit_volume = sum(threshold_smooth - discharge, na.rm = TRUE),
      deficit_mean   = mean(threshold_smooth - discharge, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(duration >= min_duration) %>%
    mutate(
      event_id_new = row_number(),
      year  = year(start_date),
      month = month(start_date),
      season = case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5)  ~ "Spring",
        month %in% c(6, 7, 8)  ~ "Summer",
        month %in% c(9, 10,11) ~ "Fall"
      )
    )
  
  if (!is.null(station_id)) event_info$station_id <- station_id
  
  data <- data %>%
    left_join(event_info %>% dplyr::select(event_id, event_id_new), by = "event_id")
  
  cat(sprintf("  [Pipeline A] Drought events: %d (min %d days below Q20 threshold)\n",
              nrow(event_info), min_duration))
  return(list(
    daily_data         = data,
    drought_events     = event_info,
    n_events           = nrow(event_info),
    total_drought_days = sum(event_info$duration),
    message            = "Success"
  ))
}

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
    } else if (dist_name == "lnorm") {
      par  <- lmomco::parlnorm(lmom_sample)
      lmom_theo <- lmomco::lmomlnorm(par)
    } else if (dist_name == "weibull") {
      par  <- lmomco::parwei(lmom_sample)
      lmom_theo <- lmomco::lmomwei(par)
    } else if (dist_name == "gpareto") {
      par  <- lmomco::pargpa(lmom_sample)
      lmom_theo <- lmomco::lmomgpa(par)
    } else if (dist_name == "llogis") {
      par  <- lmomco::parlld(lmom_sample)
      lmom_theo <- lmomco::lmomlld(par)
    } else if (dist_name == "pe3") {
      par  <- lmomco::parpe3(lmom_sample)
      lmom_theo <- lmomco::lmompe3(par)
    } else {
      return(Inf)
    }
    
    # Euclidean distance on first two L-moments (L-mean, L-scale)
    d_sample <- lmom_sample$lambdas[1:2]
    d_theo   <- lmom_theo$lambdas[1:2]
    if (any(is.na(c(d_sample, d_theo)))) return(Inf)
    return(sqrt(sum((d_sample - d_theo)^2)))
    
  }, error = function(e) Inf)
}

calculate_ssi <- function(monthly_data, station_id = NULL) {
  if (is.null(monthly_data) || nrow(monthly_data) < MIN_MONTHS_FOR_FIT) {
    cat("  [Pipeline B] Insufficient monthly data for SSI (need >= 120 months)\n")
    return(NULL)
  }
  
  flow_vals <- monthly_data$discharge_mean
  flow_vals <- flow_vals[flow_vals > 0 & !is.na(flow_vals)]
  if (length(flow_vals) < MIN_MONTHS_FOR_FIT) {
    cat("  [Pipeline B] Insufficient positive flow values for SSI\n")
    return(NULL)
  }
  
  cat("  [Pipeline B] Fitting distributions (L-moments criterion)...\n")
  
  # Fit all candidate distributions and compute L-moments distance
  dist_distances <- data.frame(
    distribution = DISTRIBUTIONS,
    distance     = vapply(DISTRIBUTIONS,
                          function(d) calculate_lmoments_distance(flow_vals, d),
                          numeric(1)),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(dist_distances))) {
    cat(sprintf("    %-10s L-moments distance = %.6f\n",
                dist_distances$distribution[i], dist_distances$distance[i]))
  }
  
  best_dist <- dist_distances$distribution[which.min(dist_distances$distance)]
  cat(sprintf("  Best distribution: %s (distance = %.6f)\n",
              best_dist, min(dist_distances$distance, na.rm = TRUE)))
  
  # Fit the best distribution using lmomco's L-moment parameter estimation
  tryCatch({
    lmom_sample <- lmomco::lmoms(flow_vals, nmom = 4)
    
    cdf_fn <- switch(best_dist,
                     gev     = function(x) lmomco::cdfgev(x, lmomco::pargev(lmom_sample)),
                     lnorm   = function(x) lmomco::cdflnorm(x, lmomco::parlnorm(lmom_sample)),
                     weibull = function(x) lmomco::cdfwei(x, lmomco::parwei(lmom_sample)),
                     gpareto = function(x) lmomco::cdfgpa(x, lmomco::pargpa(lmom_sample)),
                     llogis  = function(x) lmomco::cdflld(x, lmomco::parlld(lmom_sample)),
                     pe3     = function(x) lmomco::cdfpe3(x, lmomco::parpe3(lmom_sample))
    )
    
    # Compute CDF values then transform to standard normal (SSI)
    p_vals <- cdf_fn(monthly_data$discharge_mean)
    # Clamp to avoid Inf at exact 0 or 1
    p_vals <- pmax(1e-6, pmin(1 - 1e-6, p_vals))
    monthly_data$ssi              <- qnorm(p_vals)
    monthly_data$best_distribution <- best_dist
    
    cat(sprintf("  [Pipeline B] SSI calculated for %d months\n",
                sum(!is.na(monthly_data$ssi))))
    
    return(list(
      monthly_data      = monthly_data,
      best_distribution = best_dist,
      dist_distances    = dist_distances
    ))
    
  }, error = function(e) {
    cat(sprintf("  [Pipeline B] SSI calculation error: %s\n", e$message))
    return(NULL)
  })
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
  
  data           <- ssi_result$monthly_data
  data$in_drought <- data$ssi < drought_threshold
  data$event_id   <- with(rle(data$in_drought), rep(seq_along(lengths), lengths))
  
  event_info <- data %>%
    group_by(event_id) %>%
    filter(in_drought == TRUE) %>%
    summarise(
      start_date      = min(date),
      end_date        = max(date),
      duration_months = n(),
      start_year      = first(year),
      start_month     = first(month),
      ssi_min         = min(ssi,  na.rm = TRUE),
      ssi_mean        = mean(ssi, na.rm = TRUE),
      severity        = sum(abs(ssi[ssi < drought_threshold]), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(duration_months >= 1) %>%
    mutate(
      event_id_new = row_number(),
      season = case_when(
        start_month %in% c(12, 1, 2) ~ "Winter",
        start_month %in% c(3, 4, 5)  ~ "Spring",
        start_month %in% c(6, 7, 8)  ~ "Summer",
        start_month %in% c(9, 10,11) ~ "Fall"
      )
    )
  
  if (!is.null(station_id)) {
    event_info$station_id    <- station_id
    event_info$distribution  <- ssi_result$best_distribution
  }
  
  data <- data %>%
    left_join(event_info %>% dplyr::select(event_id, event_id_new), by = "event_id")
  
  cat(sprintf("  [Pipeline B] Drought events: %d (SSI < %.2f)\n",
              nrow(event_info), drought_threshold))
  return(list(
    monthly_data         = data,
    drought_events       = event_info,
    n_events             = nrow(event_info),
    total_drought_months = sum(event_info$duration_months),
    best_distribution    = ssi_result$best_distribution,
    message              = "Success"
  ))
}

# ============================================================================
# SECTION 12: PROCESS ALL STATIONS
#
# Shared steps [1-6]: load → completeness → gap-fill → threshold/SSI → drought events
# ============================================================================

all_processed_data <- list()   # gap-filled daily data (shared)
all_completeness   <- list()   # completeness check results (shared)
all_thresholds     <- list()   # Pipeline A
all_daily_droughts <- list()   # Pipeline A
all_ssi_results    <- list()   # Pipeline B
all_ssi_droughts   <- list()   # Pipeline B

for (i in 1:nrow(stations)) {
  station_id   <- stations$StationID[i]
  station_name <- stations$StationName[i]
  display_name <- paste0(station_id, " - ", station_name)
  
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("PROCESSING STATION %d/%d: %s\n", i, nrow(stations), display_name))
  cat(sprintf("%s\n", strrep("=", 60)))
  
  # ── [1/6] Load ─────────────────────────────────────────────────────────────
  cat("\n[1/6] Loading discharge data...\n")
  raw_data <- load_discharge_data(station_id)
  if (is.null(raw_data)) {
    cat(sprintf("  SKIPPED: no data file found for %s [%s]\n", station_id, display_name))
    next
  }
  
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
  
  # ── [4a/6] Pipeline A: variable threshold ───────────────────────────────────
  cat("\n[4a/6] Pipeline A — calculating variable threshold (Q20, 80% exceedance)...\n")
  thresholds <- calculate_variable_threshold(data_filled, station_id = station_id)
  all_thresholds[[station_id]] <- thresholds
  
  # ── [5a/6] Pipeline A: daily drought identification ─────────────────────────
  cat("\n[5a/6] Pipeline A — identifying daily drought events...\n")
  daily_droughts <- identify_droughts(data_filled, thresholds, station_id = station_id)
  all_daily_droughts[[station_id]] <- daily_droughts
  
  # ── [4b/6] Pipeline B: aggregate to monthly ─────────────────────────────────
  cat("\n[4b/6] Pipeline B — aggregating to monthly...\n")
  monthly_data <- aggregate_to_monthly(data_filled)
  
  # ── [5b/6] Pipeline B: SSI calculation ─────────────────────────────────────
  cat("\n[5b/6] Pipeline B — fitting distributions & calculating SSI...\n")
  ssi_result <- calculate_ssi(monthly_data, station_id = station_id)
  if (!is.null(ssi_result)) {
    all_ssi_results[[station_id]] <- ssi_result
    
    # ── [6b/6] Pipeline B: SSI drought identification ─────────────────────────
    cat("\n[6b/6] Pipeline B — identifying SSI drought events (SSI < -0.84)...\n")
    ssi_droughts <- identify_ssi_droughts(ssi_result, station_id = station_id)
    all_ssi_droughts[[station_id]] <- ssi_droughts
  } else {
    cat("  [Pipeline B] SSI skipped — insufficient monthly data\n")
  }
  
  cat(sprintf("\n  Station %s complete.\n", station_id))
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

# Pipeline A
saveRDS(all_thresholds,
        file.path(MAIN_OUTPUT_DIR, "daily", "thresholds", "thresholds_all.rds"))
saveRDS(all_daily_droughts,
        file.path(MAIN_OUTPUT_DIR, "daily", "drought_events", "daily_droughts_all.rds"))
cat("  Pipeline A (daily threshold) outputs saved\n")

# Pipeline B
for (sid in names(all_ssi_results)) {
  saveRDS(all_ssi_results[[sid]]$monthly_data,
          file.path(MAIN_OUTPUT_DIR, "ssi", "monthly_data",
                    paste0(sid, "_monthly_ssi.rds")))
}
saveRDS(all_ssi_results,
        file.path(MAIN_OUTPUT_DIR, "ssi", "drought_events", "ssi_results_all.rds"))
saveRDS(all_ssi_droughts,
        file.path(MAIN_OUTPUT_DIR, "ssi", "drought_events", "ssi_droughts_all.rds"))
cat("  Pipeline B (SSI monthly) outputs saved\n")

# Master combined file
all_results <- list(
  stations           = stations,
  completeness       = all_completeness,
  processed_data     = all_processed_data,
  # Pipeline A
  thresholds         = all_thresholds,
  daily_droughts     = all_daily_droughts,
  # Pipeline B
  ssi_results        = all_ssi_results,
  ssi_droughts       = all_ssi_droughts,
  metadata = list(
    input_data_dir            = INPUT_DATA_DIR,
    output_dir                = MAIN_OUTPUT_DIR,
    # Shared preprocessing
    completeness_denominator  = "month_span (first_day_of_first_month to last_day_of_last_month)",
    min_completeness_pct      = MIN_COMPLETENESS_PCT,
    max_fill_days             = MAX_FILL_DAYS,
    interpolation_method      = "monoH.FC (Fritsch-Carlson PCHIP)",
    # Pipeline A
    pipeline_a                = "Raut & Ganguli (2024) variable threshold",
    threshold_exceedance_prob = 0.80,        # 80% exceedance
    threshold_quantile_prob   = THRESHOLD_EXCEEDANCE_PROB,   # probs= = 0.20
    threshold_window_days     = THRESHOLD_WINDOW_SIZE,
    min_drought_duration_days = MIN_DROUGHT_DURATION_DAYS,
    # Pipeline B
    pipeline_b                = "Peña-Angulo et al. (2022) SSI monthly",
    ssi_drought_threshold     = SSI_DROUGHT_THRESHOLD,
    ssi_recovery_threshold    = SSI_RECOVERY_THRESHOLD,
    distributions_tested      = DISTRIBUTIONS,
    distribution_selection    = "minimum L-moments Euclidean distance (lmomco)",
    min_months_for_fit        = MIN_MONTHS_FOR_FIT,
    min_drought_duration_months = 1L,
    processed_date            = Sys.Date(),
    script_version            = "4.0_merged_dual_pipeline"
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
    stringsAsFactors = FALSE
  )
  if (!comp$meets_threshold) {
    base$A_Drought_Events      <- NA_integer_
    base$A_Total_Days          <- NA_integer_
    base$A_Avg_Duration_Days   <- NA_real_
    base$B_Best_Distribution   <- NA_character_
    base$B_Drought_Events      <- NA_integer_
    base$B_Total_Months        <- NA_integer_
    base$B_Avg_Duration_Months <- NA_real_
    base$B_Avg_Severity        <- NA_real_
  } else {
    # Pipeline A stats
    ev_a <- all_daily_droughts[[sid]]$drought_events
    base$A_Drought_Events     <- if (!is.null(ev_a)) nrow(ev_a)           else NA_integer_
    base$A_Total_Days         <- if (!is.null(ev_a)) sum(ev_a$duration)   else NA_integer_
    base$A_Avg_Duration_Days  <- if (!is.null(ev_a) && nrow(ev_a) > 0)
      round(mean(ev_a$duration), 1)           else NA_real_
    
    # Pipeline B stats
    ev_b <- all_ssi_droughts[[sid]]$drought_events
    base$B_Best_Distribution    <- if (!is.null(all_ssi_results[[sid]]))
      all_ssi_results[[sid]]$best_distribution else NA_character_
    base$B_Drought_Events       <- if (!is.null(ev_b)) nrow(ev_b)                     else NA_integer_
    base$B_Total_Months         <- if (!is.null(ev_b)) sum(ev_b$duration_months)      else NA_integer_
    base$B_Avg_Duration_Months  <- if (!is.null(ev_b) && nrow(ev_b) > 0)
      round(mean(ev_b$duration_months), 1)             else NA_real_
    base$B_Avg_Severity         <- if (!is.null(ev_b) && nrow(ev_b) > 0)
      round(mean(ev_b$severity), 2)                   else NA_real_
  }
  summary_rows[[sid]] <- base
}

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table,
          file.path(MAIN_OUTPUT_DIR, "reports", "station_summary_dual_pipeline.csv"),
          row.names = FALSE)
print(summary_table, row.names = FALSE)

cat("\nOutputs written to: ", MAIN_OUTPUT_DIR, "/\n")
cat("  shared/processed_data/            — gap-filled daily data (both pipelines)\n")
cat("  daily/thresholds/                 — Pipeline A: Q20 variable thresholds\n")
cat("  daily/drought_events/             — Pipeline A: daily drought events\n")
cat("  ssi/monthly_data/                 — Pipeline B: monthly SSI time series\n")
cat("  ssi/drought_events/               — Pipeline B: SSI drought events\n")
cat("  figures/                          — station_completeness_barplot.png\n")
cat("  reports/                          — station_summary_dual_pipeline.csv\n")
cat("                                    — nechako_streamflow_drought_summary.xlsx\n")
cat("                                    — nechako_streamflow_drought_report.txt\n")
cat("  all_results_dual_pipeline.rds     — master combined file\n")