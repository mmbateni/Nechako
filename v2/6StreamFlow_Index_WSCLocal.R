# ============================================================================
# NECHAKO STREAMFLOW RIVER BASIN — SSI DROUGHT ANALYSIS
# Script version: 5.6_naturalized_to_2025
#
# CHANGES IN 5.6
#   - Monthly completeness (aggregate_to_monthly): every month now records
#     n_days (valid daily values after gap filling), n_observed_days (values
#     originally present in the input) and n_interpolated_days (PCHIP-filled),
#     plus a month_rule column (accepted / rejected_valid_days /
#     rejected_observed_days / no_data). fill_missing_data() flags filled days
#     in a new 'interpolated' column.
#   - MIN_VALID_DAYS_PER_MONTH (15) replaces the hard-coded 15-day rule, and
#     the new MIN_OBSERVED_DAYS_PER_MONTH sets a minimum number of OBSERVED
#     days per month. Default 0 = v5.5 behaviour (primary run). Set e.g. 25L
#     for a sensitivity run; outputs then go to streamflow_results_obs<N>/ so
#     the primary outputs read by the downstream scripts are not overwritten.
#   - Completeness counts are reported in the summary CSV, Excel workbook,
#     text report and master RDS metadata.
#
# CHANGES IN 5.5
#   - Naturalized input confirmed: dat_natural_<ID>.csv in NATURALIZED_DATA_DIR
#     is the authoritative, updated file. For 08JC001 this is
#     dat_natural_08JC001.csv (Stephen's verified values to 2021-12-15, then
#     the reconstruction to 2025-12-31, with an Input_status column). The
#     *_Stephen copies and NatFlow_*_Justin files in the same folder are older
#     versions and are never read (NATURALIZED_FILES names files explicitly).
#   - Documented that the naturalization uses NO meteorological predictors
#     (precipitation, snow or temperature): the naturalized series is not
#     constructed from the ERA5-Land inputs used by Scripts 4, 5 and 7, so any
#     correlation with them is physical, not built in. Recorded as
#     NATURALIZATION_MET_PREDICTORS = FALSE in the master RDS metadata, the
#     text report and the summary tables; Script 7 reads it as a provenance
#     check.
#   - Startup check that every NATURALIZED_FILES entry exists; the naturalized
#     file used per station is written to the summary CSV / Excel workbook.
#
# CHANGES IN 5.4
#   - Naturalized inputs extended to 2025-12-31. Pipeline B now reads the
#     updated files explicitly (NATURALIZED_FILES): Stephen's verified values
#     to 2021-12-15, then the water-balance reconstruction (public WSC data +
#     Nechako stage-storage table). The earlier dat_natural_<ID>.csv versions
#     (inconsistent from 2022) were replaced by the updated files under the
#     same name; the *_Stephen copies are no longer picked up.
#   - ANALYSIS_END_DATE = 2025-12-31 applied to every station in Pipelines B
#     and C, so all series end on the same day (08JE001 observed otherwise
#     runs into 2026 and gives a partial January 2026).
#   - Input_status column of the updated files (verified / approved /
#     provisional / persisted) is carried to the monthly series. Each month
#     gets input_status = worst daily status; counts appear in the summary
#     CSV, Excel workbook and text report.
#   - PROVISIONAL_HANDLING: "keep" (default) uses provisional/persisted days
#     and flags them; "exclude" drops them (series ends 2025-10-12) for a
#     sensitivity run. Oct-Dec 2025 rests on provisional WSC data.
#
# CHANGES IN 5.3
#   - TOLERANT_SELECTION = FALSE: a candidate that fits every value ("strict")
#     is always preferred; the relaxed rule (<= MAX_OUT_OF_SUPPORT values
#     outside the fitted support, clamped to SSI +/-4.75) is only a fallback.
#     With TRUE, relaxed candidates competed equally and produced many clamped
#     +/-4.75 values (most Pipeline C months of 08JA017 / 08JC001).
#   - ln3 values at or below the fitted lower bound (zeta) now get CDF = 0 and
#     are counted as out-of-support like every other bounded candidate,
#     instead of returning NaN and rejecting ln3 outright.
#   - POST_REGULATION_START_YEAR = 1960 (was 1952): drops the 1952-1959
#     reservoir-filling years, same start as NTSDI_START_YEAR and the O-SWSI
#     baseline in 7WBInex_NTSDI_GeoLakes_LocalRS.r.
#   - SSI plots use Date-class infinities for the background bands (removes the
#     "A <numeric> value was passed to a Date scale" warnings).
#
# Pipeline B (Monthly): Standardized Streamflow Index, Peña-Angulo et al. (2022)
#   Regulated Nechako main-stem stations (REGULATED_STATIONS) — NATURALIZED
#   flow, loaded from NATURALIZED_DATA_DIR (NATURALIZED_FILES):
#     08JA017 NECHAKO RIVER BELOW CHESLATTA FALLS
#     08JC001 NECHAKO RIVER AT VANDERHOOF      (dat_natural_08JC001.csv)
#     08JC002 NECHAKO RIVER AT ISLE PIERRE     (dat_natural_08JC002.csv)
#   Unregulated reference station — observed flow (INPUT_DATA_DIR):
#     08JE001 STUART RIVER NEAR FORT ST. JAMES
#
#   NOTE: naturalized flow = synthesized flow representing no-dam conditions.
#   The files start in 1955, i.e. after Kenney Dam closure (1952), and are
#   naturalized throughout (routing model, Kokoszka 2026). The naturalization
#   uses NO meteorological predictors (precipitation, snow, temperature), so
#   the series is not constructed from the ERA5-Land climate inputs used
#   elsewhere in the framework.
#   08JC001/08JC002: Stephen's verified series to 2021-12-15, water-balance
#   reconstruction 2021-12-16 to 2025-12-31 (Input_status column).
#   08JA017: confirm its source — if it came from the same batch as the
#   faulty earlier 08JC001/08JC002 files, its 2022+ values need the same
#   correction.
#
# Pipeline C (Monthly): Same SSI methodology as Pipeline B, but on the
#   ORIGINAL (observed WSC) record for the regulated stations, restricted to
#   the post-regulation period (years >= POST_REGULATION_START_YEAR), with no
#   naturalization applied; plus 08JE001 observed (same series as Pipeline B).
#   Observed files come from INPUT_DATA_DIR (data_downloads_geomet_api,
#   <ID>_WSC_DISCHARGE.csv); any station without a file is skipped with a
#   warning.
#
# Shared preprocessing:
#   - Data loading & completeness check (month-span denominator)
#   - Gap filling: Fritsch-Carlson PCHIP, gaps > MAX_FILL_DAYS left as NA;
#     filled days are flagged (interpolated = TRUE)
#   - Hard exclusion gate for stations < MIN_COMPLETENESS_PCT completeness
#
# Both pipelines add:
#   - Daily → monthly mean aggregation with the monthly completeness rule
#     (MIN_VALID_DAYS_PER_MONTH valid days and MIN_OBSERVED_DAYS_PER_MONTH
#     observed days)
#   - Best-fit distribution selection via L-moments (lmomco), fit SEPARATELY
#     per calendar month (Jan..Dec each get their own candidate contest,
#     support check, close-margin bootstrap and CDF), matching standard
#     SPI/SSI practice (McKee et al. 1993; Vicente-Serrano et al. 2012;
#     Peña-Angulo et al. 2022). Requires >= MIN_YEARS_PER_MONTH years of that
#     calendar month; otherwise that month's SSI is NA.
#   - Support check (v5.2): a candidate is "strict" if every fitting value lies
#     inside its fitted support. If no candidate is strict, the best-ranked
#     candidate with <= MAX_OUT_OF_SUPPORT out-of-support values is used
#     ("relaxed", marked * in reports); those values are clamped to SSI ±4.75.
#     This prevents whole calendar months from being dropped.
#   - SSI drought identification with the framework-wide event definition
#     from DROUGHT_ANALYSIS_utils.R: onset SSI < DROUGHT_ONSET, termination
#     SSI >= DROUGHT_END, minimum duration get_event_min_duration(SSI_EVENT_SCALE),
#     intensity = mean(DROUGHT_END - SSI), severity = intensity x duration
#
# Pipeline B additionally adds:
#   - Pettitt change-point screening on the annual-mean series of each
#     naturalized station, plus the same test on the unregulated control
#     station(s) (HOMOGENEITY_CONTROL_STATIONS) over the naturalized period.
#     A change-point that also appears at the control station points to a
#     climatic origin rather than a naturalization artifact. Results are
#     carried into the summary CSV, Excel workbook, text report and RDS files.
# ============================================================================

rm(list = ls())
options(warn = 1)   # print warnings as they occur, so they appear in the log

# ---- Working directory (framework rule) -------------------------------------
# Run every script from the project root: set the working directory there
# BEFORE running (setwd(), or open the RStudio project). All paths below are
# derived from getwd(); no path is hard-coded.
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")

source(UTILS_FILE, local = .GlobalEnv)
if (!exists("detect_drought_events", mode = "function")) {
  stop("DROUGHT_ANALYSIS_utils.R was sourced, but detect_drought_events() is unavailable.")
}

SCRIPT_VERSION <- "5.6_naturalized_to_2025"

# Date-class infinities for annotate("rect") on date axes (as in Script 7);
# structure() instead of as.Date(-Inf) so it also works on R < 4.3.
DATE_NEG_INF <- structure(-Inf, class = "Date")
DATE_POS_INF <- structure( Inf, class = "Date")

# ============================================================================
# SECTION 1: PACKAGES
# ============================================================================
packages_needed <- c(
  # Shared / general
  "tidyverse", "lubridate", "zoo",
  # SSI
  "lmomco", "fitdistrplus",
  # Pettitt change-point screening
  "trend",
  # Excel output
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
# Regulated Nechako main-stem stations (below Kenney Dam) with a naturalized-
# flow reconstruction in NATURALIZED_DATA_DIR (NATURALIZED_FILES).
# Pipeline B: naturalized record. Pipeline C: observed record, post-regulation only.
REGULATED_STATIONS <- c("08JA017", "08JC001", "08JC002")

# Unregulated station(s) used as a control for the Pettitt change-point test
HOMOGENEITY_CONTROL_STATIONS <- c("08JE001")

stations <- data.frame(
  StationID   = c("08JA017", "08JC001", "08JC002", "08JE001"),
  StationName = c("NECHAKO RIVER BELOW CHESLATTA FALLS",
                  "NECHAKO RIVER AT VANDERHOOF",
                  "NECHAKO RIVER AT ISLE PIERRE",
                  "STUART RIVER NEAR FORT ST. JAMES"),
  StartYear   = NA_real_,   # filled from the data once loaded
  EndYear     = NA_real_,
  stringsAsFactors = FALSE
)

# Same station list for Pipeline C (R copies on assignment, so each pipeline
# fills in its own StartYear/EndYear independently)
stations_C <- stations

# Single name lookup used everywhere (independent of which pipeline table is
# being modified)
STATION_LOOKUP <- setNames(stations$StationName, stations$StationID)

get_station_name <- function(station_id) {
  nm <- unname(STATION_LOOKUP[station_id])
  if (length(nm) == 0 || is.na(nm)) "UNKNOWN" else nm
}
get_station_display <- function(station_id) {
  paste0(station_id, " - ", get_station_name(station_id))
}

# Observed WSC records (Pipeline C regulated stations, and 08JE001 in both pipelines)
INPUT_DATA_DIR  <- file.path(WD_PATH, "Hydrology",
                             "data_retrievalWaterSurveyofCanada", "data_downloads_geomet_api")

# Naturalized flow files, used by Pipeline B only
NATURALIZED_DATA_DIR <- file.path(WD_PATH, "Hydrology",
                                  "data_retrievalWaterSurveyofCanada", "naturalized_flows")

# Exact file per naturalized station. Named explicitly so the older *_Stephen
# copies and NatFlow_*_Justin files in the same folder are never picked up by
# pattern matching. dat_natural_08JC001.csv / dat_natural_08JC002.csv are the
# authoritative updated files (verified to 2021-12-15 + reconstruction to
# 2025-12-31, Input_status column). 08JA017: set to the file you trust.
NATURALIZED_FILES <- c(
  "08JA017" = "dat_natural_08JA017.csv",   # <-- confirm source / 2022+ values
  "08JC001" = "dat_natural_08JC001.csv",   # authoritative: verified + reconstruction, Input_status
  "08JC002" = "dat_natural_08JC002.csv"    # authoritative: verified + reconstruction, Input_status
)
NATURALIZED_FILE_PATHS <- setNames(file.path(NATURALIZED_DATA_DIR, NATURALIZED_FILES),
                                   names(NATURALIZED_FILES))

# Naturalization method metadata (v5.5). The naturalized series is a
# statistically synthesized no-dam flow; no precipitation, snow or temperature
# predictors are used. Script 7 (SWSI) reads this from the master RDS to
# document that Qnat is not constructed from the ERA5-Land Precip/SWE inputs.
NATURALIZATION_MET_PREDICTORS <- FALSE
NATURALIZATION_NOTE <- paste(
  "Naturalized flow = synthesized flow under no-dam conditions (statistical",
  "naturalization). No precipitation, snow or temperature predictors are used,",
  "so the series is not constructed from the ERA5-Land climate inputs used",
  "elsewhere in the framework; any correlation with them is physical.")

# Startup check (v5.5): every configured naturalized file must exist
cat("Naturalized input files (Pipeline B):\n")
for (sid in names(NATURALIZED_FILE_PATHS)) {
  if (file.exists(NATURALIZED_FILE_PATHS[[sid]])) {
    cat(sprintf("  %s : %s\n", sid, NATURALIZED_FILES[[sid]]))
  } else {
    warning(sprintf("Naturalized input for %s not found: %s (station will be skipped in Pipeline B).",
                    sid, NATURALIZED_FILE_PATHS[[sid]]))
  }
}
cat(sprintf("  Meteorological predictors in naturalization: %s\n",
            if (isTRUE(NATURALIZATION_MET_PREDICTORS)) "YES" else "NONE"))

naturalized_file_of <- function(sid) {
  if (sid %in% REGULATED_STATIONS && sid %in% names(NATURALIZED_FILES))
    unname(NATURALIZED_FILES[[sid]]) else NA_character_
}

# ── Monthly completeness rule (v5.6) ────────────────────────────────────────
# A monthly mean is accepted only if the month has
#   >= MIN_VALID_DAYS_PER_MONTH    daily values after gap filling
#                                  (observed + PCHIP-interpolated), AND
#   >= MIN_OBSERVED_DAYS_PER_MONTH daily values that were originally observed.
# Primary run     : MIN_OBSERVED_DAYS_PER_MONTH <- 0L   (v5.5 behaviour)
# Sensitivity run : MIN_OBSERVED_DAYS_PER_MONTH <- 25L  (stricter, observed-only)
# A sensitivity run writes to its own output root (streamflow_results_obs<N>)
# so the primary outputs read by downstream scripts are never overwritten.
MIN_VALID_DAYS_PER_MONTH    <- 15L
MIN_OBSERVED_DAYS_PER_MONTH <- 0L
stopifnot(MIN_VALID_DAYS_PER_MONTH >= 1L, MIN_OBSERVED_DAYS_PER_MONTH >= 0L,
          MIN_OBSERVED_DAYS_PER_MONTH <= 31L)
MONTHLY_RULE_LABEL <- sprintf(">= %d valid days (observed + interpolated)%s",
                              MIN_VALID_DAYS_PER_MONTH,
                              if (MIN_OBSERVED_DAYS_PER_MONTH > 0L)
                                sprintf(" AND >= %d observed days (sensitivity)", MIN_OBSERVED_DAYS_PER_MONTH)
                              else ", no observed-day minimum (primary)")

# Single unified output root — all pipelines write under their own subdirs
MAIN_OUTPUT_DIR <- if (MIN_OBSERVED_DAYS_PER_MONTH > 0L) {
  sprintf("streamflow_results_obs%d", MIN_OBSERVED_DAYS_PER_MONTH)
} else "streamflow_results"
cat(sprintf("Monthly completeness rule: %s\nOutput root: %s\n", MONTHLY_RULE_LABEL, MAIN_OUTPUT_DIR))

# Spatial-support interpretation for sparse hydrometric observations.
# These are catchment-integrated discharge records, not raster cells and not
# independent samples of a spatial drought field. Do NOT average SSI values
# across stations to manufacture a basin/subbasin mean. Note the stations are
# nested: 08JA017 ⊂ 08JC001 ⊂ 08JC002, and 08JC002 (Isle Pierre) lies below
# the Stuart River confluence, so it is not independent of 08JE001.
STATION_SPATIAL_SUPPORT <- data.frame(
  StationID   = stations$StationID,
  StationName = stations$StationName,
  HydrologicUnit = c("Nechako River main stem below Kenney Dam / Cheslatta Falls reach",
                     "Nechako River main stem / upstream contributing catchment (contains 08JA017)",
                     "Nechako River main stem below Stuart confluence (contains 08JC001 and Stuart River inflow)",
                     "Stuart River subbasin / upstream contributing catchment"),
  SpatialSupport = "Single hydrometric gauge; discharge integrates upstream flow",
  BasinMeanRepresentation = "NO",
  AggregationRule = "Retain as a separate station series; do not average SSI with other stations",
  CoordinatesStatus = "Not supplied in this script; join authoritative WSC station coordinates if mapping is required",
  stringsAsFactors = FALSE
)

# Pipeline tag used in console messages by the shared SSI functions
set_pipeline_tag <- function(tag) options(ssi.pipeline_tag = tag)
ptag <- function() getOption("ssi.pipeline_tag", "SSI")

for (sub in c(
  "",
  "shared/processed_data",          # daily gap-filled data
  "shared/completeness",
  "ssi/monthly_data",               # Pipeline B
  "ssi/drought_events",
  "ssi_C/monthly_data",             # Pipeline C
  "ssi_C/drought_events",
  "figures",
  "reports"
)) {
  dir.create(file.path(MAIN_OUTPUT_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

# ============================================================================
# SECTION 2.2: STATION COMPLETENESS PLOT
# ============================================================================
create_station_completeness_plot <- function(stations, completeness, output_dir,
                                             file_name = "station_completeness_barplot.png",
                                             plot_label = "Station Data Completeness") {
  
  station_status <- data.frame(
    StationID    = names(completeness),
    Completeness = sapply(completeness, function(x) x$completeness),
    Status       = sapply(completeness, function(x) ifelse(x$meets_threshold, "INCLUDED", "EXCLUDED")),
    stringsAsFactors = FALSE
  )
  
  station_map_data <- merge(stations, station_status, by = "StationID")
  
  cat("  [PLOT] Station status summary:\n")
  cat(sprintf("    Included: %d stations\n", sum(station_map_data$Status == "INCLUDED")))
  cat(sprintf("    Excluded: %d stations\n", sum(station_map_data$Status == "EXCLUDED")))
  
  png_file <- file.path(output_dir, "figures", file_name)
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
          main = paste0(plot_label, " (Month-Span Denominator)"),
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
# SECTION 2.3: SUMMARY HELPERS + EXCEL SUMMARY WORKBOOK
# ============================================================================
# One-row data frame of Pettitt results for a station (NA columns if untested)
homogeneity_fields <- function(hg) {
  if (!is.null(hg) && isTRUE(hg$tested)) {
    data.frame(
      Homogeneity_Role       = hg$role,
      Homogeneity_Period     = sprintf("%d-%d", hg$first_year, hg$last_year),
      Homogeneity_ChangeYear = hg$change_year,
      Homogeneity_PValue     = round(hg$p_value, 4),
      Homogeneity_PctShift   = round(hg$pct_shift, 1),
      Homogeneity_Flag       = ifelse(hg$significant, "INHOMOGENEOUS", "homogeneous"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Homogeneity_Role       = if (!is.null(hg) && !is.null(hg$role)) hg$role else NA_character_,
      Homogeneity_Period     = NA_character_,
      Homogeneity_ChangeYear = NA_real_,
      Homogeneity_PValue     = NA_real_,
      Homogeneity_PctShift   = NA_real_,
      Homogeneity_Flag       = NA_character_,
      stringsAsFactors = FALSE
    )
  }
}

# v5.4: months whose naturalized input includes provisional/persisted days
provisional_months_str <- function(res) {
  if (is.null(res) || is.null(res$monthly_data) ||
      !"input_status" %in% names(res$monthly_data)) return(NA_character_)
  md <- res$monthly_data
  ym <- md$year_month[md$input_status %in% PROVISIONAL_STATUSES]
  if (!length(ym)) return("none")
  sprintf("%d (%s)", length(ym), paste(ym, collapse = ", "))
}

# v5.6: monthly completeness counts (accepted / rejected / interpolated)
monthly_completeness_str <- function(res) {
  if (is.null(res) || is.null(res$monthly_data) ||
      !"month_rule" %in% names(res$monthly_data)) return(NA_character_)
  md <- res$monthly_data
  sprintf("accepted %d (%d with interpolated days) | rejected: %d (< %d valid days), %d (< %d observed days) | no data %d",
          sum(md$month_rule == "accepted"),
          sum(md$month_rule == "accepted" & md$n_interpolated_days > 0L),
          sum(md$month_rule == "rejected_valid_days"), MIN_VALID_DAYS_PER_MONTH,
          sum(md$month_rule == "rejected_observed_days"), MIN_OBSERVED_DAYS_PER_MONTH,
          sum(md$month_rule == "no_data"))
}

# "n_valid/n_total" SSI month coverage string
ssi_coverage_str <- function(res) {
  if (is.null(res) || is.null(res$monthly_data)) return(NA_character_)
  md <- res$monthly_data
  sprintf("%d/%d", sum(is.finite(md$ssi)), nrow(md))
}

create_excel_summary <- function(stations, completeness,
                                 daily_droughts, ssi_droughts,
                                 ssi_results, output_dir,
                                 stations_C = NULL, completeness_C = NULL,
                                 ssi_droughts_C = NULL, ssi_results_C = NULL,
                                 homogeneity = NULL) {
  if (!require("writexl", quietly = TRUE)) {
    cat("  [EXCEL] writexl package not available - skipping Excel output\n")
    return(NULL)
  }
  
  excel_data <- list()
  
  # Sheet 1: Pipeline B station summary
  summary_rows <- list()
  for (sid in names(completeness)) {
    comp  <- completeness[[sid]]
    sname <- if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN"
    
    base <- data.frame(
      StationID        = sid,
      StationName      = sname,
      StationDisplay   = paste0(sid, " - ", sname),
      DataSource       = ifelse(sid %in% REGULATED_STATIONS, "Naturalized flow", "Observed"),
      NaturalizedFile  = naturalized_file_of(sid),
      NaturalizationMetPredictors = if (sid %in% REGULATED_STATIONS)
        ifelse(isTRUE(NATURALIZATION_MET_PREDICTORS), "YES", "NONE") else NA_character_,
      Period           = sprintf("%d-%d", comp$start_year, comp$end_year),
      Completeness_pct = round(comp$completeness, 1),
      Status           = ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED"),
      MonthlyRule      = MONTHLY_RULE_LABEL,
      stringsAsFactors = FALSE
    )
    
    if (comp$meets_threshold && !is.null(ssi_results[[sid]])) {
      base$B_Best_Distribution  <- ssi_results[[sid]]$best_distribution_str
      base$B_SSI_Coverage       <- ssi_coverage_str(ssi_results[[sid]])
      base$B_Monthly_Completeness <- monthly_completeness_str(ssi_results[[sid]])
      base$B_Provisional_Months <- provisional_months_str(ssi_results[[sid]])
    } else {
      base$B_Best_Distribution  <- NA_character_
      base$B_SSI_Coverage       <- NA_character_
      base$B_Monthly_Completeness <- NA_character_
      base$B_Provisional_Months <- NA_character_
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
    
    base <- cbind(base, homogeneity_fields(if (!is.null(homogeneity)) homogeneity[[sid]] else NULL))
    summary_rows[[sid]] <- base
  }
  excel_data[["01_Station_Summary"]] <- do.call(rbind, summary_rows)
  
  # Pipeline B drought events (one sheet per station with events)
  for (sid in names(ssi_droughts)) {
    if (!is.null(ssi_droughts[[sid]]$drought_events)) {
      ev <- ssi_droughts[[sid]]$drought_events
      if (nrow(ev) > 0) {
        excel_data[[paste0("B_", sid, "_SSI_Droughts")]] <- as.data.frame(ev)
      }
    }
  }
  
  # Pipeline B distribution fit results (includes n_out_of_support / selected)
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
  
  # Pipeline C: kept separate from Pipeline B
  if (!is.null(completeness_C) && length(completeness_C) > 0L) {
    c_summary_rows <- list()
    for (sid in names(completeness_C)) {
      comp <- completeness_C[[sid]]
      ev  <- if (!is.null(ssi_droughts_C[[sid]])) ssi_droughts_C[[sid]]$drought_events else NULL
      res <- if (!is.null(ssi_results_C[[sid]])) ssi_results_C[[sid]] else NULL
      c_summary_rows[[sid]] <- data.frame(
        StationID = sid,
        StationName = if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN",
        DataSource = ifelse(sid %in% REGULATED_STATIONS,
                            sprintf("Observed, post-regulation only (year >= %d)", POST_REGULATION_START_YEAR),
                            "Observed"),
        Period = if (is.finite(comp$start_year) && is.finite(comp$end_year)) sprintf("%d-%d", comp$start_year, comp$end_year) else NA_character_,
        Completeness_pct = round(comp$completeness, 1),
        Status = ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED"),
        MonthlyRule = MONTHLY_RULE_LABEL,
        C_Best_Distribution = if (!is.null(res)) res$best_distribution_str else NA_character_,
        C_SSI_Coverage = ssi_coverage_str(res),
        C_Monthly_Completeness = monthly_completeness_str(res),
        C_Drought_Events = if (!is.null(ev)) nrow(ev) else NA_integer_,
        C_Total_Months = if (!is.null(ev) && nrow(ev) > 0) sum(ev$duration_months) else if (!is.null(ev)) 0L else NA_integer_,
        C_Avg_Duration_Months = if (!is.null(ev) && nrow(ev) > 0) round(mean(ev$duration_months), 1) else NA_real_,
        C_Avg_Severity = if (!is.null(ev) && nrow(ev) > 0) round(mean(ev$severity), 2) else NA_real_,
        stringsAsFactors = FALSE
      )
      if (!is.null(ev) && nrow(ev) > 0) excel_data[[paste0("C_", sid, "_SSI_Droughts")]] <- as.data.frame(ev)
      if (!is.null(res) && !is.null(res$dist_distances)) {
        dd <- res$dist_distances
        dd$StationID <- sid
        excel_data[[paste0("C_", sid, "_DistFits")]] <- dd
      }
    }
    excel_data[["02_Pipeline_C_Summary"]] <- do.call(rbind, c_summary_rows)
  }
  
  xlsx_file <- file.path(output_dir, "reports", "nechako_streamflow_drought_summary.xlsx")
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  writexl::write_xlsx(excel_data, xlsx_file)
  
  cat(sprintf("  [EXCEL] Summary workbook saved: %s (%d sheets)\n",
              basename(xlsx_file), length(excel_data)))
  return(xlsx_file)
}

# ============================================================================
# SECTION 2.4: DETAILED TEXT REPORT
# ============================================================================
create_text_report <- function(stations, completeness,
                               daily_droughts, ssi_droughts,
                               ssi_results, output_dir, params,
                               completeness_C = NULL, ssi_droughts_C = NULL,
                               ssi_results_C = NULL, post_regulation_start_year = NA_integer_,
                               homogeneity = NULL) {
  
  report_file <- file.path(output_dir, "reports", "nechako_streamflow_drought_report.txt")
  
  cat("============================================================\n",
      "NECHAKO STREAMFLOW DROUGHT ANALYSIS - DETAILED REPORT\n",
      "Monthly SSI: Pipeline B (naturalized) + Pipeline C (observed)\n",
      "============================================================\n\n",
      file = report_file, sep = "")
  
  cat(sprintf("Report generated: %s\n", Sys.time()),
      file = report_file, append = TRUE)
  cat(sprintf("Script version: %s\n", params$script_version),
      file = report_file, append = TRUE)
  cat(sprintf("Analysis end date: %s; provisional naturalized input: %s\n",
              as.character(ANALYSIS_END_DATE), PROVISIONAL_HANDLING),
      file = report_file, append = TRUE)
  cat(sprintf("Monthly completeness rule: %s (output root: %s)\n\n",
              MONTHLY_RULE_LABEL, MAIN_OUTPUT_DIR),
      file = report_file, append = TRUE)
  
  cat("========== PIPELINE B STATION SUMMARY ==========\n",
      file = report_file, append = TRUE)
  cat(sprintf("Naturalized stations: %s; observed: %s\n",
              paste(REGULATED_STATIONS, collapse = ", "),
              paste(setdiff(stations$StationID, REGULATED_STATIONS), collapse = ", ")),
      file = report_file, append = TRUE)
  cat(sprintf("Naturalized input files: %s\n",
              paste(sprintf("%s=%s", names(NATURALIZED_FILES), NATURALIZED_FILES), collapse = ", ")),
      file = report_file, append = TRUE)
  cat(sprintf("Naturalization meteorological predictors: %s\n",
              if (isTRUE(NATURALIZATION_MET_PREDICTORS)) "YES"
              else "NONE (no precipitation, snow or temperature predictors)"),
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
  
  for (sid in names(completeness)) {
    comp  <- completeness[[sid]]
    sname <- if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN"
    
    cat(sprintf("\n----- Station: %s -----\n", sid), file = report_file, append = TRUE)
    cat(sprintf("Name: %s\n", sname), file = report_file, append = TRUE)
    cat(sprintf("Data source: %s\n",
                ifelse(sid %in% REGULATED_STATIONS,
                       sprintf("Naturalized flow (%s)", naturalized_file_of(sid)), "Observed")),
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
        cat("\nPipeline B (Monthly SSI):\n", file = report_file, append = TRUE)
        cat(sprintf("  Best distribution (per calendar month; * = relaxed support rule): %s\n",
                    ssi_results[[sid]]$best_distribution_str),
            file = report_file, append = TRUE)
        cat(sprintf("  SSI coverage (valid/total months): %s\n",
                    ssi_coverage_str(ssi_results[[sid]])),
            file = report_file, append = TRUE)
        cat(sprintf("  Monthly completeness: %s\n",
                    monthly_completeness_str(ssi_results[[sid]])),
            file = report_file, append = TRUE)
        pm <- provisional_months_str(ssi_results[[sid]])
        if (!is.na(pm))
          cat(sprintf("  Months with provisional/persisted input (handling = %s): %s\n",
                      PROVISIONAL_HANDLING, pm), file = report_file, append = TRUE)
        cat(sprintf("  Drought events: %d\n", nrow(ev_b)), file = report_file, append = TRUE)
        cat(sprintf("  Total drought months: %d\n", sum(ev_b$duration_months)),
            file = report_file, append = TRUE)
        cat(sprintf("  Average severity: %.2f\n", mean(ev_b$severity)),
            file = report_file, append = TRUE)
      }
    }
    
    hg <- if (!is.null(homogeneity)) homogeneity[[sid]] else NULL
    if (!is.null(hg) && isTRUE(hg$tested)) {
      cat(sprintf("  Pettitt test (%s, %d-%d): change-point %d, p = %.4f, shift %+.1f%% -> %s\n",
                  hg$role, hg$first_year, hg$last_year, hg$change_year, hg$p_value,
                  hg$pct_shift, ifelse(hg$significant, "significant", "not significant")),
          file = report_file, append = TRUE)
    }
  }
  
  if (!is.null(completeness_C) && length(completeness_C) > 0L) {
    cat("\n\n========== PIPELINE C: OBSERVED SSI (POST-REGULATION WHERE APPLICABLE) ==========\n",
        file = report_file, append = TRUE)
    cat(sprintf("%s observations restricted to year >= %s; 08JE001 remains observed flow.\n",
                paste(REGULATED_STATIONS, collapse = ", "),
                as.character(post_regulation_start_year)),
        file = report_file, append = TRUE)
    for (sid in names(completeness_C)) {
      comp <- completeness_C[[sid]]
      ev  <- if (!is.null(ssi_droughts_C[[sid]])) ssi_droughts_C[[sid]]$drought_events else NULL
      res <- if (!is.null(ssi_results_C[[sid]])) ssi_results_C[[sid]] else NULL
      cat(sprintf("\nStation %s (%s): period %s-%s; completeness %.2f%%; status %s\n",
                  sid, if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN",
                  comp$start_year, comp$end_year, comp$completeness,
                  ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED")),
          file = report_file, append = TRUE)
      if (!is.null(res)) {
        cat(sprintf("  Best distribution by month (* = relaxed support rule): %s\n", res$best_distribution_str),
            file = report_file, append = TRUE)
        cat(sprintf("  SSI coverage (valid/total months): %s\n", ssi_coverage_str(res)),
            file = report_file, append = TRUE)
        cat(sprintf("  Monthly completeness: %s\n", monthly_completeness_str(res)),
            file = report_file, append = TRUE)
      }
      if (!is.null(ev)) {
        cat(sprintf("  SSI drought events: %s; total event months: %s\n",
                    nrow(ev), if (nrow(ev) > 0) sum(ev$duration_months) else 0),
            file = report_file, append = TRUE)
      } else {
        cat("  SSI drought events: unavailable (SSI/event detection did not complete).\n",
            file = report_file, append = TRUE)
      }
    }
  }
  
  cat("\n\n========== METHODOLOGY ==========\n", file = report_file, append = TRUE)
  # Pipeline A (Daily) — DISABLED
  # cat("Pipeline A (Daily):\n",
  #     "  • Variable threshold (Q20 = 80% exceedance probability)\n",
  #     "  • 31-day centered moving average smoothing\n",
  #     "  • Minimum drought duration: 30 days\n",
  #     "  • Reference: Raut & Ganguli (2024)\n\n",
  #     file = report_file, append = TRUE)
  cat("Pipelines B and C (Monthly):\n",
      sprintf("  • Monthly mean accepted when: %s; PCHIP fills only gaps <= %d days\n",
              MONTHLY_RULE_LABEL, MAX_FILL_DAYS),
      "  • Standardized Streamflow Index (SSI), fitted separately per calendar month\n",
      "  • Distribution selection: minimum L-kurtosis distance |tau4_sample - tau4_fitted| (lmomco)\n",
      sprintf("  • Support rule: strict (all fitting values inside fitted support) if possible; otherwise relaxed (<= %d values outside, clamped to SSI ±4.75), marked *\n",
              MAX_OUT_OF_SUPPORT),
      sprintf("  • Drought events: onset SSI < %.2f, termination SSI >= %.2f, minimum duration %d months\n",
              SSI_DROUGHT_THRESHOLD, SSI_RECOVERY_THRESHOLD, get_event_min_duration(SSI_EVENT_SCALE)),
      "  • Intensity = mean(termination threshold - SSI); severity = intensity x duration\n",
      "  • Pettitt change-point test on annual means (naturalized stations + unregulated control)\n",
      sprintf("  • Naturalization: %s\n", NATURALIZATION_NOTE),
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
MAX_FILL_DAYS        <- 14   # PCHIP only for gaps ≤ 14 days
# Monthly completeness rule: MIN_VALID_DAYS_PER_MONTH / MIN_OBSERVED_DAYS_PER_MONTH
# are set in SECTION 2 because they also select the output root.

# ── Pipeline A: daily variable threshold — DISABLED ──────────────────────
# THRESHOLD_EXCEEDANCE_PROB <- 0.20   # Q20 = 80% exceedance (low-flow threshold)
# THRESHOLD_WINDOW_SIZE     <- 31     # 31-day centred moving average smoothing
# MIN_DROUGHT_DURATION_DAYS <- 30     # minimum consecutive days below threshold

# ── SSI monthly (Peña-Angulo et al. 2022) ───────────────────────────────
# Event thresholds and minimum duration come from DROUGHT_ANALYSIS_utils.R
# (DROUGHT_ONSET, DROUGHT_END, EVENT_MIN_DURATION_BY_SCALE); edit them there.
SSI_DROUGHT_THRESHOLD   <- DROUGHT_ONSET
SSI_RECOVERY_THRESHOLD  <- DROUGHT_END
SSI_EVENT_SCALE         <- 1L   # monthly SSI (no accumulation) -> scale-1 duration rule
# NOTE: "ln3" is the 3-parameter (shifted) log-normal (lmomco parln3/cdfln3),
# NOT the 2-parameter log-normal.
DISTRIBUTIONS <- c("gev", "pe3", "ln3", "glo", "gpareto", "weibull")
MIN_MONTHS_FOR_FIT     <- 120      # minimum total monthly records before any per-month fitting
MIN_YEARS_PER_MONTH    <- 10       # minimum years of THAT calendar month required to fit it

# ── Distribution support rule (v5.2) ────────────────────────────────────
# A fitted CDF value <= 1e-6 or >= 1 - 1e-6 on a fitting observation means
# that observation lies outside (or at) the fitted distribution's support.
# MAX_OUT_OF_SUPPORT: max number of such observations tolerated when no
#   candidate fits all values ("relaxed" selection; values clamped to ±4.75).
# TOLERANT_SELECTION: FALSE (default) = strict candidates always preferred,
#   relaxed only as a fallback, so months that already had a strict fit are
#   unchanged. TRUE = rank ALL candidates with <= MAX_OUT_OF_SUPPORT bad values
#   by L-moment distance (picks better-fitting candidates more often, at the
#   cost of more clamped extreme values). Keep FALSE for analysis runs: every
#   clamped value becomes an artificial SSI of +/-4.75.
MAX_OUT_OF_SUPPORT <- 2L
TOLERANT_SELECTION <- FALSE

# ── Pettitt change-point screening ──────────────────────────────────────
HOMOGENEITY_ALPHA <- 0.05          # significance level
# Start year for the control-station test. NA = use the earliest start year of
# the naturalized stations, so control and naturalized series cover the same years.
HOMOGENEITY_CONTROL_START_YEAR <- NA_real_

# ── Pipeline C: observed regulated records — post-regulation filter ─────
# Kenney Dam closure / Nechako Reservoir impoundment began in 1952. Pipeline C
# keeps only year >= POST_REGULATION_START_YEAR of the observed record for the
# regulated stations. 1952-1959 is the reservoir-filling period (with winter
# data gaps at 08JC001 to 1955), so the operating regime starts in 1960, the
# same year as NTSDI_START_YEAR / OSWSI_BASELINE_START in Script 7. That year
# is itself a placeholder there: confirm it with BC Hydro / Rio Tinto and keep
# the two scripts identical.
POST_REGULATION_START_YEAR <- 1960

# ── Common analysis end date (v5.4) ──────────────────────────────────────
# Every station in Pipelines B and C is cut at this date so all series end on
# the same day. NA = keep each file's own last date.
ANALYSIS_END_DATE <- as.Date("2025-12-31")

# ── Provisional naturalized input (v5.4) ─────────────────────────────────
# Updated naturalized files carry Input_status per day:
#   "verified (Stephen)" / "approved" -> final inputs
#   "provisional"  -> WSC level/flow not yet approved (2025-10-13 .. 12-23)
#   "persisted"    -> spillway/Kemano held at last value (2025-12-24 .. 12-31)
# "keep"    = use them; months are flagged via monthly input_status.
# "exclude" = drop those days (series ends 2025-10-12); sensitivity run.
PROVISIONAL_HANDLING <- "keep"
PROVISIONAL_STATUSES <- c("provisional", "persisted")

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
    # Naturalized-flow CSVs, e.g. dat_natural_08JC001.csv (exact name, so the
    # *_Stephen / NatFlow_*_Justin variants in the same folder are not picked up)
    file.path(data_path, paste0("dat_natural_", station_id, ".csv"))
  )
  # v5.4: naturalized stations read ONLY their explicitly named file
  if (identical(normalizePath(data_path, mustWork = FALSE),
                normalizePath(NATURALIZED_DATA_DIR, mustWork = FALSE)) &&
      station_id %in% names(NATURALIZED_FILES)) {
    file_patterns <- file.path(data_path, NATURALIZED_FILES[[station_id]])
  }
  data <- NULL
  for (file_name in file_patterns) {
    if (file.exists(file_name)) {
      cat(sprintf("  Loading: %s [%s]\n", basename(file_name), display_name))
      data <- read.csv(file_name, stringsAsFactors = FALSE)
      break
    }
  }
  if (is.null(data)) {
    warning(sprintf("  No data file found for station %s [%s] in %s", station_id, display_name, data_path))
    return(NULL)
  }
  
  colnames(data) <- tolower(colnames(data))
  
  # Naturalized-flow exports (dat_natural_*.csv) are in long format:
  # station, date, parameter, value, symbol. Keep only the Flow rows.
  if ("parameter" %in% colnames(data)) {
    data <- data[!is.na(data$parameter) & tolower(trimws(data$parameter)) == "flow", ]
  }
  
  date_col      <- grep("date",              colnames(data), value = TRUE, ignore.case = TRUE)
  discharge_col <- grep("discharge|flow|q_", colnames(data), value = TRUE, ignore.case = TRUE)
  if (length(discharge_col) == 0) {
    # Long-format naturalized files store discharge in a generic "value" column
    discharge_col <- grep("^value$", colnames(data), value = TRUE, ignore.case = TRUE)
  }
  if (length(date_col) == 0 | length(discharge_col) == 0) {
    stop(sprintf("Cannot find date/discharge columns for station %s", station_id))
  }
  
  # ── Robust date parsing ─────────────────────────────────────────────────
  # WSC exports are ISO "YYYY-MM-DD"; naturalized exports are "M/D/YYYY".
  # as.Date() defaults silently misparse "M/D/YYYY", so use parse_date_time()
  # with several candidate orders.
  raw_date_str <- as.character(data[[date_col[1]]])
  parsed_date  <- suppressWarnings(
    lubridate::parse_date_time(raw_date_str, orders = c("Ymd", "mdy", "dmy"))
  )
  data$date <- as.Date(parsed_date)
  
  n_unparsed <- sum(!is.na(raw_date_str) & raw_date_str != "" & is.na(data$date))
  if (n_unparsed > 0) {
    warning(sprintf(
      "  [%s] %d date value(s) could not be parsed under any of Ymd/mdy/dmy and were set to NA. Inspect raw values, e.g.: %s",
      station_id, n_unparsed,
      paste(head(raw_date_str[!is.na(raw_date_str) & raw_date_str != "" & is.na(data$date)], 5), collapse = ", ")))
  }
  
  raw_discharge <- as.character(data[[discharge_col[1]]])
  parsed_discharge <- suppressWarnings(as.numeric(raw_discharge))
  bad_discharge <- !is.na(raw_discharge) & nzchar(trimws(raw_discharge)) & !is.finite(parsed_discharge)
  if (any(bad_discharge)) {
    warning(sprintf("[%s] %d discharge value(s) were nonnumeric/sentinel-like and will be excluded; examples: %s",
                    station_id, sum(bad_discharge),
                    paste(head(raw_discharge[bad_discharge], 5), collapse = ", ")))
  }
  data$discharge <- parsed_discharge
  data <- data[!is.na(data$date) & is.finite(data$discharge), ]
  
  # v5.4: per-day input status of the updated naturalized files
  if ("input_status" %in% colnames(data)) {
    data$input_status <- trimws(as.character(data$input_status))
    n_prov <- sum(data$input_status %in% PROVISIONAL_STATUSES)
    if (n_prov > 0) {
      if (identical(PROVISIONAL_HANDLING, "exclude")) {
        data <- data[!data$input_status %in% PROVISIONAL_STATUSES, ]
        cat(sprintf("  [%s] PROVISIONAL_HANDLING = exclude: %d provisional/persisted day(s) dropped\n",
                    station_id, n_prov))
      } else {
        cat(sprintf("  [%s] %d provisional/persisted day(s) kept and flagged (%s to %s)\n",
                    station_id, n_prov,
                    min(data$date[data$input_status %in% PROVISIONAL_STATUSES]),
                    max(data$date[data$input_status %in% PROVISIONAL_STATUSES])))
      }
    }
  } else {
    data$input_status <- NA_character_
  }
  
  # v5.4: common analysis end date for all stations and pipelines
  if (length(ANALYSIS_END_DATE) == 1L && !is.na(ANALYSIS_END_DATE)) {
    n_after <- sum(data$date > ANALYSIS_END_DATE)
    if (n_after > 0) {
      data <- data[data$date <= ANALYSIS_END_DATE, ]
      cat(sprintf("  [%s] %d record(s) after ANALYSIS_END_DATE (%s) removed\n",
                  station_id, n_after, ANALYSIS_END_DATE))
    }
  }
  if (anyDuplicated(data$date)) {
    dup_dates <- unique(data$date[duplicated(data$date)])
    stop(sprintf("[%s] Duplicate dates found in discharge input (e.g. %s). Resolve duplicates before SSI fitting.",
                 station_id, paste(head(dup_dates, 5), collapse = ", ")))
  }
  # Negative discharge values are physically invalid (likely sentinel codes)
  invalid_q <- data$discharge < 0
  if (any(invalid_q)) {
    warning(sprintf("[%s] Removing %d negative discharge value(s), likely invalid/sentinel values; examples: %s",
                    station_id, sum(invalid_q),
                    paste(head(data$discharge[invalid_q], 5), collapse = ", ")))
    data <- data[!invalid_q, ]
  }
  if (nrow(data) == 0L) stop("No valid date/discharge records remain for station ", station_id, ".")
  data$year       <- year(data$date)
  data$month      <- month(data$date)
  data$day        <- day(data$date)
  data$doy        <- yday(data$date)
  data$station_id <- station_id
  
  # Sanity check for the M/D/YYYY misparse signature
  if (nrow(data) > 0 && (max(data$year) - min(data$year)) < 2 && min(data$year) < 1800) {
    warning(sprintf(
      "  [%s] Parsed year range (%d-%d) looks implausible -- likely a date-format misparse. Raw date sample: %s",
      station_id, min(data$year), max(data$year),
      paste(head(raw_date_str, 5), collapse = ", ")))
  }
  
  cat(sprintf("  Records loaded: %d (%s to %s) [%s]\n",
              nrow(data), min(data$date), max(data$date), display_name))
  # Guard: a naturalized file without Input_status that reaches past 2021 is
  # likely from the faulty earlier batch (inconsistent 2022+ values)
  if (station_id %in% REGULATED_STATIONS && all(is.na(data$input_status)) &&
      max(data$year) >= 2022 &&
      identical(normalizePath(data_path, mustWork = FALSE),
                normalizePath(NATURALIZED_DATA_DIR, mustWork = FALSE))) {
    warning(sprintf("[%s] naturalized file has no Input_status column but extends to %d: check that its 2022+ values are from the corrected reconstruction.",
                    station_id, max(data$year)))
  }
  return(data)
}

# Pipeline B routes the regulated stations to the naturalized-flow directory.
# Pipeline C routes them to the observed WSC directory (then filtered to the
# post-regulation period). 08JE001 always uses INPUT_DATA_DIR.
get_data_path_for_station <- function(station_id, pipeline = c("B", "C")) {
  pipeline <- match.arg(pipeline)
  if (station_id %in% REGULATED_STATIONS && pipeline == "B") {
    return(NATURALIZED_DATA_DIR)
  }
  return(INPUT_DATA_DIR)
}

# ============================================================================
# SECTION 5: DATA COMPLETENESS CHECK  (month-span denominator)
#
# Completeness = actual records / calendar days from the first day of the
# month of first observation to the last day of the month of last observation.
# Large internal gaps (> MAX_FILL_DAYS) are flagged informatively.
# ============================================================================
check_data_completeness <- function(data, station_id) {
  display_name <- get_station_display(station_id)
  name <- get_station_name(station_id)
  
  if (is.null(data) || nrow(data) == 0) {
    return(list(station_id = station_id, station_name = name, start_year = NA, end_year = NA,
                total_years = 0, total_days = 0, expected_days = 0,
                completeness = 0, meets_threshold = FALSE,
                message = "No data available"))
  }
  
  first_obs_date <- min(data$date)
  last_obs_date  <- max(data$date)
  
  start_date <- as.Date(paste0(format(first_obs_date, "%Y-%m"), "-01"))
  
  end_year  <- as.integer(format(last_obs_date, "%Y"))
  end_month <- as.integer(format(last_obs_date, "%m"))
  
  if (end_month == 12) {
    end_date <- as.Date(paste0(end_year + 1, "-01-01")) - 1
  } else {
    end_date <- as.Date(paste0(end_year, "-", sprintf("%02d", end_month + 1), "-01")) - 1
  }
  
  expected_days <- as.numeric(difftime(end_date, start_date, units = "days")) + 1
  actual_days   <- nrow(data)
  completeness  <- actual_days / expected_days * 100
  
  year_range  <- range(data$year)
  total_years <- diff(year_range) + 1
  
  cat(sprintf("\nStation: %s [%s]\n", station_id, display_name))
  cat(sprintf("  Period: %s to %s (%d years)\n", start_date, end_date, total_years))
  cat(sprintf("  Data: %d / %d days (%.2f%% complete)\n",
              actual_days, expected_days, completeness))
  cat(sprintf("  (Month-span: %s to %s)\n", start_date, end_date))
  
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
    station_id       = station_id,
    station_name     = name,
    display_name     = display_name,
    start_year       = year_range[1],
    end_year         = year_range[2],
    start_date       = start_date,
    end_date         = end_date,
    total_years      = total_years,
    total_days       = actual_days,
    expected_days    = expected_days,
    completeness     = completeness,
    n_large_gaps     = length(large_gaps),
    large_gap_sizes  = gap_lengths[large_gaps],
    meets_threshold  = meets_threshold,
    denominator_type = "month_span"
  ))
}

# ============================================================================
# SECTION 6: MISSING DATA FILLING — Fritsch-Carlson PCHIP
#
# Monotone cubic Hermite spline (splinefun method = "monoH.FC") for gaps
# <= MAX_FILL_DAYS; longer gaps are left as NA. v5.6: filled days are flagged
# in the 'interpolated' column so monthly aggregation can count observed and
# interpolated days separately.
# ============================================================================
fill_missing_data <- function(data, max_fill = MAX_FILL_DAYS) {
  
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  full_dates    <- seq(min(data$date), max(data$date), by = "day")
  complete_data <- merge(data.frame(date = full_dates), data, by = "date", all.x = TRUE)
  complete_data <- complete_data[order(complete_data$date), ]
  
  # Immutable copy of observed values: interpolated values must not become
  # spline anchors for later gaps in the same pass.
  q_observed <- complete_data$discharge
  q <- q_observed
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
    if (left_idx < 1L || right_idx > n)            next
    if (is.na(q[left_idx]) || is.na(q[right_idx])) next
    
    # Local real observations within ±15 days for spline curvature
    window_radius <- 15L
    ctx_idx <- max(1L, left_idx - window_radius):min(n, right_idx + window_radius)
    obs     <- !is.na(q_observed[ctx_idx])
    x_obs   <- ctx_idx[obs]
    y_obs   <- q_observed[ctx_idx[obs]]
    if (length(x_obs) < 2) next
    
    spline_fn            <- splinefun(x = x_obs, y = y_obs, method = "monoH.FC")
    q[run_start:run_end] <- pmax(0, spline_fn(run_start:run_end))
    filled_short         <- filled_short + run_len
  }
  
  complete_data$discharge    <- q
  complete_data$interpolated <- is.na(q_observed) & !is.na(q)   # v5.6
  complete_data$year         <- year(complete_data$date)
  complete_data$month        <- month(complete_data$date)
  complete_data$day          <- day(complete_data$date)
  complete_data$doy          <- yday(complete_data$date)
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
# SECTION 9: DAILY TO MONTHLY AGGREGATION
#
# v5.6: per month
#   n_days              valid daily values after gap filling
#   n_observed_days     values present in the input file
#   n_interpolated_days PCHIP-filled values (gaps <= MAX_FILL_DAYS)
#   month_rule          accepted | rejected_valid_days | rejected_observed_days | no_data
# The monthly mean is kept only for accepted months (see MONTHLY_RULE_LABEL).
# ============================================================================
aggregate_to_monthly <- function(daily_data,
                                 min_valid_days    = MIN_VALID_DAYS_PER_MONTH,
                                 min_observed_days = MIN_OBSERVED_DAYS_PER_MONTH) {
  if (is.null(daily_data) || nrow(daily_data) == 0) return(NULL)
  
  if (!"input_status" %in% colnames(daily_data)) daily_data$input_status <- NA_character_
  if (!"interpolated" %in% colnames(daily_data)) daily_data$interpolated <- FALSE
  daily_data$interpolated[is.na(daily_data$interpolated)] <- FALSE
  # v5.4: month status = worst daily status (persisted > provisional > final)
  worst_status <- function(s) {
    s <- s[!is.na(s)]
    if (!length(s)) return(NA_character_)
    if ("persisted"   %in% s) return("persisted")
    if ("provisional" %in% s) return("provisional")
    "final"
  }
  monthly_observed <- daily_data %>%
    dplyr::group_by(year, month) %>%
    dplyr::summarise(
      n_days              = sum(!is.na(discharge)),
      n_interpolated_days = sum(!is.na(discharge) & interpolated),
      n_observed_days     = n_days - n_interpolated_days,
      discharge_mean_all  = if (n_days > 0L) mean(discharge, na.rm = TRUE) else NA_real_,
      input_status        = worst_status(input_status),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      month_rule = dplyr::case_when(
        n_days == 0L                         ~ "no_data",
        n_days < min_valid_days              ~ "rejected_valid_days",
        n_observed_days < min_observed_days  ~ "rejected_observed_days",
        TRUE                                 ~ "accepted"),
      discharge_mean = ifelse(month_rule == "accepted", discharge_mean_all, NA_real_)
    ) %>%
    dplyr::select(-discharge_mean_all)
  
  first_month  <- as.Date(format(min(daily_data$date), "%Y-%m-01"))
  last_month   <- as.Date(format(max(daily_data$date), "%Y-%m-01"))
  month_starts <- seq(first_month, last_month, by = "month")
  monthly_grid <- data.frame(
    date  = month_starts,   # canonical first-of-month dates
    year  = lubridate::year(month_starts),
    month = lubridate::month(month_starts)
  )
  monthly_data <- dplyr::left_join(monthly_grid, monthly_observed,
                                   by = c("year", "month")) %>%
    dplyr::mutate(
      n_days              = tidyr::replace_na(n_days, 0L),
      n_observed_days     = tidyr::replace_na(n_observed_days, 0L),
      n_interpolated_days = tidyr::replace_na(n_interpolated_days, 0L),
      month_rule          = tidyr::replace_na(month_rule, "no_data"),
      year_month = paste0(year, "-", sprintf("%02d", month))
    ) %>%
    dplyr::select(year, month, discharge_mean, n_days, n_observed_days, n_interpolated_days,
                  month_rule, date, year_month, input_status)
  
  cat(sprintf("  [%s] Monthly aggregation (%s): %d calendar months | accepted %d (%d with interpolated days) | rejected: %d (< %d valid days), %d (< %d observed days) | no data %d\n",
              ptag(), MONTHLY_RULE_LABEL, nrow(monthly_data),
              sum(monthly_data$month_rule == "accepted"),
              sum(monthly_data$month_rule == "accepted" & monthly_data$n_interpolated_days > 0L),
              sum(monthly_data$month_rule == "rejected_valid_days"), min_valid_days,
              sum(monthly_data$month_rule == "rejected_observed_days"), min_observed_days,
              sum(monthly_data$month_rule == "no_data")))
  n_flag <- sum(monthly_data$input_status %in% PROVISIONAL_STATUSES)
  if (n_flag > 0)
    cat(sprintf("  [%s] Months with provisional/persisted input: %d (%s)\n", ptag(), n_flag,
                paste(monthly_data$year_month[monthly_data$input_status %in% PROVISIONAL_STATUSES],
                      collapse = ", ")))
  return(monthly_data)
}

# ============================================================================
# SECTION 10: L-MOMENTS DISTRIBUTION SELECTION & SSI
#
# All six candidates are 3-parameter distributions fitted by L-moments, which
# reproduce the sample mean, L-scale and L-skewness exactly. The only
# discriminating ratio is the L-kurtosis tau4, so the distance is
# |tau4_sample - tau4_fitted| (Hosking & Wallis 1997), computed by
# lmom_ratio_distance() in DROUGHT_ANALYSIS_utils.R. A candidate that cannot
# reproduce the sample tau3 gets distance Inf.
# ============================================================================
SSI_LMOMCO_TYPE <- c(gev = "gev", ln3 = "ln3", weibull = "wei",
                     gpareto = "gpa", glo = "glo", pe3 = "pe3")

calculate_lmoments_distance <- function(data, dist_name) {
  type <- SSI_LMOMCO_TYPE[dist_name]
  if (is.na(type)) return(Inf)
  tryCatch({
    lmom_sample <- lmomco::lmoms(data, nmom = 4)
    lmom_ratio_distance(lmom_sample, unname(type))$distance
  }, error = function(e) Inf)
}

# ============================================================================
# SECTION 10b: BOOTSTRAP DISTRIBUTION STABILITY
#
# Resamples a month's flow values B times and tallies how often each eligible
# candidate wins the L-moment distance contest. Run only when the winner's
# margin over the runner-up is < CLOSE_MARGIN_PCT.
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
  
  cat(sprintf("    [%s] Bootstrap stability (B = %d replicates):\n", ptag(), B))
  for (d in names(sort(win_counts, decreasing = TRUE))) {
    cat(sprintf("      %-10s  %3.0f%%\n", d, win_pct[d]))
  }
  cat(sprintf("    Bootstrap plurality winner: %s (%3.0f%% of draws)\n",
              boot_best, win_pct[boot_best]))
  
  return(list(win_counts = win_counts,
              win_pct    = win_pct,
              boot_best  = boot_best))
}

# ── Per-calendar-month distribution fit + CDF/SSI for ONE month's data ──────
fit_month_distribution_and_ssi <- function(flow_vals_m, target_vals_m,
                                           month_label, station_id) {
  sid_label <- ifelse(is.null(station_id), "unknown", station_id)
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
  
  make_cdf <- function(d, lmom_sample) {
    switch(d,
           gev     = function(x) lmomco::cdfgev(x, lmomco::pargev(lmom_sample)),
           pe3     = function(x) lmomco::cdfpe3(x, lmomco::parpe3(lmom_sample)),
           # ln3 is lower-bounded at zeta; lmomco returns NaN for x <= zeta
           # (log of a non-positive number). Those values are below the
           # support, so their CDF is 0 and they are counted as out-of-support
           # (and clamped in the SSI transform) like any other bounded fit.
           ln3     = {
             par_ln3 <- lmomco::parln3(lmom_sample)
             zeta    <- unname(par_ln3$para[["zeta"]])
             function(x) {
               p <- rep(NA_real_, length(x))
               ok <- is.finite(x)
               below <- ok & x <= zeta
               p[below] <- 0
               if (any(ok & !below)) p[ok & !below] <- lmomco::cdfln3(x[ok & !below], par_ln3)
               p
             }
           },
           glo     = function(x) lmomco::cdfglo(x, lmomco::parglo(lmom_sample)),
           gpareto = function(x) lmomco::cdfgpa(x, lmomco::pargpa(lmom_sample)),
           weibull = function(x) lmomco::cdfwei(x, lmomco::parwei(lmom_sample)),
           NULL
    )
  }
  
  ranked <- dist_distances[is.finite(dist_distances$distance), , drop = FALSE]
  ranked <- ranked[order(ranked$distance), , drop = FALSE]
  lmom_sample_m <- tryCatch(lmomco::lmoms(flow_vals_m, nmom = 4), error = function(e) NULL)
  best_dist      <- NA_character_
  cdf_fn         <- NULL
  clamp_p        <- 1e-6
  candidate_cdfs <- list()
  n_out_by_cand  <- setNames(rep(Inf, nrow(ranked)), ranked$distribution)
  
  # ── Support check: count fitting values outside each candidate's support ──
  if (!is.null(lmom_sample_m) && nrow(ranked) > 0L) {
    for (cand in ranked$distribution) {
      candidate_cdf <- tryCatch(make_cdf(cand, lmom_sample_m), error = function(e) NULL)
      if (is.null(candidate_cdf)) {
        cat(sprintf("    Candidate %s rejected: parameter estimation failed.\n", cand))
        next
      }
      p_fit <- tryCatch(candidate_cdf(flow_vals_m),
                        error = function(e) rep(NA_real_, length(flow_vals_m)))
      if (length(p_fit) != length(flow_vals_m) || any(!is.finite(p_fit))) {
        cat(sprintf("    Candidate %s rejected: fitted CDF returned invalid values on fitting data.\n", cand))
        next
      }
      n_out <- sum(p_fit <= clamp_p | p_fit >= 1 - clamp_p)
      n_out_by_cand[cand]    <- n_out
      candidate_cdfs[[cand]] <- candidate_cdf
      if (n_out > 0L) {
        cat(sprintf("    Candidate %s: %d of %d fitting value(s) outside fitted support.\n",
                    cand, n_out, length(flow_vals_m)))
      }
    }
  }
  
  # ── Selection: strict first; relaxed (<= MAX_OUT_OF_SUPPORT) as fallback ──
  # (names(n_out_by_cand) are in L-moment distance order, so the first entry
  #  of each set is its best-ranked candidate)
  strict_ok   <- names(n_out_by_cand)[is.finite(n_out_by_cand) & n_out_by_cand == 0]
  tolerant_ok <- names(n_out_by_cand)[is.finite(n_out_by_cand) & n_out_by_cand <= MAX_OUT_OF_SUPPORT]
  valid_candidates <- if (TOLERANT_SELECTION) {
    tolerant_ok
  } else if (length(strict_ok) > 0L) {
    strict_ok
  } else {
    tolerant_ok
  }
  
  selection_mode <- NA_character_
  if (length(valid_candidates) > 0L) {
    best_dist      <- valid_candidates[1]
    cdf_fn         <- candidate_cdfs[[best_dist]]
    selection_mode <- if (n_out_by_cand[best_dist] == 0) "strict" else "relaxed"
  }
  dist_distances_valid <- ranked[ranked$distribution %in% valid_candidates, , drop = FALSE]
  
  dist_distances$n_out_of_support <- unname(n_out_by_cand[dist_distances$distribution])
  dist_distances$selected         <- !is.na(best_dist) & dist_distances$distribution == best_dist
  
  if (is.na(best_dist)) {
    warning(sprintf("[%s] No candidate distribution has <= %d out-of-support values for station %s, month %s; SSI set to NA for this month.",
                    ptag(), MAX_OUT_OF_SUPPORT, sid_label, month_label))
    cat(sprintf("    NO VALID CANDIDATE (all exceed %d out-of-support values) — SSI left NA for %s.\n",
                MAX_OUT_OF_SUPPORT, month_label))
  } else {
    cat(sprintf("    Selected: %s [%s] (|dtau4| = %.6f; %d out-of-support value(s))\n",
                best_dist, selection_mode,
                ranked$distance[match(best_dist, ranked$distribution)],
                as.integer(n_out_by_cand[best_dist])))
    if (selection_mode == "relaxed") {
      warning(sprintf("[%s] Relaxed distribution selection for station %s, month %s: '%s' has %d fitting value(s) outside its support (clamped to SSI ±%.2f).",
                      ptag(), sid_label, month_label, best_dist,
                      as.integer(n_out_by_cand[best_dist]), -qnorm(clamp_p)))
    }
  }
  
  # ── Close-margin warning and conditional bootstrap (per month) ─────────────
  dist_eligible <- dist_distances_valid[is.finite(dist_distances_valid$distance), , drop = FALSE]
  boot_result   <- NULL
  
  if (!is.na(best_dist) && nrow(dist_eligible) >= 2) {
    sorted_eligible <- dist_eligible[order(dist_eligible$distance), , drop = FALSE]
    distance_best   <- sorted_eligible$distance[1]
    distance_second <- sorted_eligible$distance[2]
    margin_pct <- if (is.finite(distance_best) && abs(distance_best) > .Machine$double.eps) {
      (distance_second - distance_best) / abs(distance_best) * 100
    } else if (is.finite(distance_second) && distance_second > 0) {
      Inf
    } else {
      0
    }
    runner_up <- sorted_eligible$distribution[2]
    
    if (is.finite(margin_pct) && margin_pct < CLOSE_MARGIN_PCT) {
      warning(sprintf(
        "[%s] UNSTABLE SELECTION for station %s, month %s: '%s' wins by only %.1f%% over '%s' — running bootstrap.",
        ptag(), sid_label, month_label, best_dist, margin_pct, runner_up))
      
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
        msg <- sprintf(
          "[%s] Bootstrap plurality ('%s', %.0f%%) DIFFERS from single-sample winner ('%s') for station %s month %s. Consider overriding.",
          ptag(), boot_result$boot_best, boot_result$win_pct[boot_result$boot_best],
          best_dist, sid_label, month_label)
        cat("    ", msg, "\n", sep = "")
        warning(msg)
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
  ssi_out <- if (is.na(best_dist) || is.null(cdf_fn)) {
    rep(NA_real_, length(target_vals_m))
  } else tryCatch({
    # Apply the fitted CDF to positive, finite targets only; missing and
    # non-positive targets remain NA.
    target_ok <- is.finite(target_vals_m) & target_vals_m > 0
    out <- rep(NA_real_, length(target_vals_m))
    
    if (any(target_ok)) {
      target_idx <- which(target_ok)
      p_vals <- cdf_fn(target_vals_m[target_ok])
      if (length(p_vals) != length(target_idx)) {
        stop("Fitted CDF returned ", length(p_vals), " probabilities for ",
             length(target_idx), " valid target values.")
      }
      
      # Clamp probabilities to [1e-6, 1 - 1e-6] so SSI is bounded at about
      # ±4.75; values outside a bounded distribution's support are clamped,
      # not set to NA.
      valid_prob <- is.finite(p_vals)
      n_clamped  <- sum(valid_prob & (p_vals < 1e-6 | p_vals > 1 - 1e-6))
      out[target_idx[valid_prob]] <- qnorm(pmin(pmax(p_vals[valid_prob], 1e-6), 1 - 1e-6))
      
      if (n_clamped > 0L) {
        cat(sprintf("    %d target value(s) for %s clamped to SSI = +/-%.2f (CDF outside [1e-6, 1-1e-6]).\n",
                    n_clamped, month_label, -qnorm(1e-6)))
      }
      if (any(!valid_prob)) {
        cat(sprintf("    %d positive target value(s) for %s have a non-finite CDF; SSI left NA.\n",
                    sum(!valid_prob), month_label))
      }
    }
    
    n_nonpositive <- sum(is.finite(target_vals_m) & target_vals_m <= 0)
    n_missing     <- sum(!is.finite(target_vals_m))
    if (n_nonpositive > 0L || n_missing > 0L) {
      cat(sprintf("    %s: left %d non-positive and %d missing/non-finite target value(s) as NA.\n",
                  month_label, n_nonpositive, n_missing))
    }
    out
  }, error = function(e) {
    cat(sprintf("    SSI calculation error for %s: %s\n", month_label, e$message))
    rep(NA_real_, length(target_vals_m))
  })
  
  dist_distances$month          <- month_label
  dist_distances$selection_mode <- selection_mode
  
  list(ssi = ssi_out, best_dist = best_dist, selection_mode = selection_mode,
       dist_distances = dist_distances, boot_result = boot_result)
}

calculate_ssi <- function(monthly_data, station_id = NULL) {
  if (is.null(monthly_data) || nrow(monthly_data) < MIN_MONTHS_FOR_FIT) {
    cat(sprintf("  [%s] Insufficient monthly data for SSI (need >= %d months total)\n",
                ptag(), MIN_MONTHS_FOR_FIT))
    return(NULL)
  }
  
  cat(sprintf("  [%s] Fitting distributions SEPARATELY per calendar month (L-moments criterion, need >= %d years/month)...\n",
              ptag(), MIN_YEARS_PER_MONTH))
  
  monthly_data$ssi          <- NA_real_
  monthly_data$distribution <- NA_character_
  
  dist_distances_all         <- list()
  best_distribution_by_month <- setNames(rep(NA_character_, 12), month.abb)
  selection_mode_by_month    <- setNames(rep(NA_character_, 12), month.abb)
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
    
    best_distribution_by_month[m]     <- fit$best_dist
    selection_mode_by_month[m]        <- fit$selection_mode
    dist_distances_all[[month_label]] <- fit$dist_distances
    if (!is.null(fit$boot_result)) boot_results_by_month[[month_label]] <- fit$boot_result
  }
  
  n_ssi <- sum(!is.na(monthly_data$ssi))
  cat(sprintf("  [%s] SSI calculated for %d / %d months (fit separately per calendar month)\n",
              ptag(), n_ssi, nrow(monthly_data)))
  
  missing_months <- month.abb[is.na(best_distribution_by_month)]
  if (length(missing_months) > 0L) {
    warning(sprintf("[%s] Station %s: no SSI for calendar month(s) %s.",
                    ptag(), ifelse(is.null(station_id), "unknown", station_id),
                    paste(missing_months, collapse = ", ")))
  }
  
  if (n_ssi == 0) {
    cat(sprintf("  [%s] No calendar month had a valid per-month fit — SSI unavailable\n", ptag()))
    return(NULL)
  }
  
  dist_distances_combined <- do.call(rbind, dist_distances_all)
  rownames(dist_distances_combined) <- NULL
  
  fitted_months <- names(best_distribution_by_month)[!is.na(best_distribution_by_month)]
  relaxed_mark  <- ifelse(selection_mode_by_month[fitted_months] == "relaxed", "*", "")
  best_dist_str <- paste(sprintf("%s=%s%s", fitted_months,
                                 best_distribution_by_month[fitted_months], relaxed_mark),
                         collapse = ", ")
  
  return(list(
    monthly_data          = monthly_data,
    best_distribution     = best_distribution_by_month,  # named by month, e.g. best_distribution["Jan"]
    selection_mode        = selection_mode_by_month,     # "strict" / "relaxed" / NA, by month
    best_distribution_str = best_dist_str,               # compact summary; * = relaxed
    dist_distances        = dist_distances_combined,     # 6 rows per fitted month, tagged by $month
    boot_result           = boot_results_by_month        # only months that triggered bootstrap
  ))
}

# ============================================================================
# SECTION 11: SSI DROUGHT EVENT IDENTIFICATION
# ============================================================================
identify_ssi_droughts <- function(ssi_result,
                                  drought_threshold  = SSI_DROUGHT_THRESHOLD,
                                  recovery_threshold = SSI_RECOVERY_THRESHOLD,
                                  station_id         = NULL) {
  if (is.null(ssi_result) || is.null(ssi_result$monthly_data)) {
    return(list(drought_events = NULL, message = "Insufficient SSI data"))
  }
  
  data <- ssi_result$monthly_data
  
  # 1. Shared detection system (parity with 12basin_timeseries.R)
  df_in <- data.frame(date = data$date, value = data$ssi)
  event_error <- NULL
  shared_events <- tryCatch(
    detect_drought_events(df_in,
                          onset_threshold = drought_threshold,
                          termination_threshold = recovery_threshold,
                          scale = SSI_EVENT_SCALE),
    error = function(e) {
      event_error <<- conditionMessage(e)
      warning(sprintf("[%s] SSI event detection failed: %s", ifelse(is.null(station_id), "unknown", station_id), event_error))
      data.frame()
    }
  )
  
  if (!is.null(event_error)) {
    data$event_id <- NA_integer_
    return(list(monthly_data = data, drought_events = NULL, n_events = NA_integer_,
                total_drought_months = NA_integer_, best_distribution = ssi_result$best_distribution,
                best_distribution_str = ssi_result$best_distribution_str,
                message = paste("Event detection failed:", event_error)))
  }
  
  if (nrow(shared_events) == 0) {
    event_info <- data.frame(start_date = as.Date(character()), end_date = as.Date(character()),
                             duration_months = integer(), start_year = integer(), start_month = integer(),
                             ssi_min = numeric(), ssi_mean = numeric(), mean_int = numeric(), severity = numeric(),
                             event_id_new = integer(), season = character(), stringsAsFactors = FALSE)
    if (!is.null(station_id)) {
      event_info$station_id   <- character()
      event_info$distribution <- character()
    }
    data$event_id <- NA_integer_
    return(list(monthly_data = data, drought_events = event_info, n_events = 0L,
                total_drought_months = 0L, best_distribution = ssi_result$best_distribution,
                best_distribution_str = ssi_result$best_distribution_str, message = "Success"))
  }
  
  # 2. Map events back to compute deficit-based intensity/severity
  data$event_id <- NA_integer_
  for (i in seq_len(nrow(shared_events))) {
    idx <- which(data$date >= shared_events$start_date[i] & data$date <= shared_events$end_date[i])
    data$event_id[idx] <- i
  }
  
  event_info <- data %>%
    dplyr::filter(!is.na(event_id)) %>%
    dplyr::group_by(event_id) %>%
    dplyr::summarise(
      start_date      = min(date),
      end_date        = max(date),
      duration_months = dplyr::n(),
      start_year      = lubridate::year(min(date)),
      start_month     = lubridate::month(min(date)),
      ssi_min         = min(ssi, na.rm = TRUE),
      ssi_mean        = mean(ssi, na.rm = TRUE),
      # Framework definition: I = mean(DROUGHT_END - SSI), S = I x D
      mean_int        = mean(recovery_threshold - ssi),
      severity        = sum(recovery_threshold - ssi),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      event_id_new = dplyr::row_number(),
      season = dplyr::case_when(
        start_month %in% c(12, 1, 2) ~ "Winter", start_month %in% c(3, 4, 5)  ~ "Spring",
        start_month %in% c(6, 7, 8)  ~ "Summer", start_month %in% c(9, 10, 11) ~ "Fall"
      )
    )
  
  if (!is.null(station_id)) {
    event_info$station_id <- station_id
    # Distribution fitted for each event's onset calendar month
    event_info$distribution <- unname(ssi_result$best_distribution[event_info$start_month])
  }
  
  data <- data %>%
    dplyr::left_join(data.frame(event_id = event_info$event_id, event_id_new = event_info$event_id_new),
                     by = "event_id") %>%
    dplyr::mutate(event_id = event_id_new) %>%
    dplyr::select(-event_id_new)
  
  cat(sprintf("  [%s] Drought events: %d (Onset < %.2f, Recovery >= %.2f, min %d months)\n",
              ptag(), nrow(event_info), drought_threshold, recovery_threshold,
              get_event_min_duration(SSI_EVENT_SCALE)))
  
  return(list(
    monthly_data          = data,
    drought_events        = event_info,
    n_events              = nrow(event_info),
    total_drought_months  = sum(event_info$duration_months),
    best_distribution     = ssi_result$best_distribution,
    best_distribution_str = ssi_result$best_distribution_str,
    message               = "Success"
  ))
}

# ============================================================================
# SECTION 11b: PETTITT CHANGE-POINT SCREENING (naturalized + control)
#
# role = "naturalized": screens a naturalized series for an undocumented
#   step-change (e.g. from how the reconstruction was assembled).
# role = "control": the same test on an unregulated station over the same
#   years. A change-point at the same date in the control suggests a
#   climatic origin (e.g. the mid-1970s Pacific regime shift) rather than a
#   naturalization artifact. Because the naturalization uses no
#   meteorological predictors, a shared change-point cannot come from a common
#   climate input; it points to a real hydroclimatic shift.
# Annual means of near-complete years (>= 300 valid days) are used to reduce
# autocorrelation. A significant result never auto-excludes a station; it is
# reported and carried into all outputs. Nested stations (08JA017 ⊂ 08JC001 ⊂
# 08JC002) are not independent tests.
# ============================================================================
check_naturalized_homogeneity <- function(daily_data, station_id,
                                          role = c("naturalized", "control"),
                                          start_year = NULL) {
  role <- match.arg(role)
  display_name <- get_station_display(station_id)
  label <- if (role == "control") "[Homogeneity control]" else "[Pipeline B]"
  
  if (is.null(daily_data) || nrow(daily_data) == 0) {
    return(list(station_id = station_id, role = role, tested = FALSE, message = "No data available"))
  }
  if (!is.null(start_year) && length(start_year) == 1L && is.finite(start_year)) {
    daily_data <- daily_data[daily_data$year >= start_year, , drop = FALSE]
  }
  
  annual <- daily_data %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(discharge_mean = mean(discharge, na.rm = TRUE),
                     n_days = sum(!is.na(discharge)), .groups = "drop") %>%
    dplyr::filter(n_days >= 300) %>%
    dplyr::arrange(year)
  
  if (nrow(annual) < 10) {
    cat(sprintf("  %s Homogeneity check skipped for %s: insufficient near-complete annual data (%d years < 10)\n",
                label, display_name, nrow(annual)))
    return(list(station_id = station_id, role = role, tested = FALSE, message = "Insufficient annual data"))
  }
  if (!requireNamespace("trend", quietly = TRUE)) {
    cat(sprintf("  %s Homogeneity check skipped for %s: 'trend' package unavailable\n", label, display_name))
    return(list(station_id = station_id, role = role, tested = FALSE, message = "'trend' package unavailable"))
  }
  
  pettitt_result <- tryCatch(trend::pettitt.test(annual$discharge_mean), error = function(e) NULL)
  if (is.null(pettitt_result)) {
    cat(sprintf("  %s Homogeneity check FAILED for %s\n", label, display_name))
    return(list(station_id = station_id, role = role, tested = FALSE, message = "Pettitt test failed"))
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
    mean_before <- NA_real_
    mean_after  <- NA_real_
    pct_shift   <- NA_real_
  }
  
  cat(sprintf("\n  %s Homogeneity check (Pettitt test, %s) — %s:\n", label, role, display_name))
  cat(sprintf("    Annual series: %d years (%d-%d)\n", nrow(annual), min(annual$year), max(annual$year)))
  cat(sprintf("    Most likely change-point: %d  (K = %.0f, p = %.4f)\n",
              change_year, pettitt_result$statistic, p_value))
  cat(sprintf("    Mean before: %.2f m3/s | Mean after: %.2f m3/s | Shift: %+.1f%%\n",
              mean_before, mean_after, pct_shift))
  
  if (significant) {
    if (role == "naturalized") {
      warning(sprintf(
        "[Pipeline B] Significant step-change for %s around %d (p = %.4f, mean shift %+.1f%%). Compare with the control station before attributing it to the naturalized reconstruction.",
        display_name, change_year, p_value, pct_shift))
    } else {
      warning(sprintf(
        "[Homogeneity control] Unregulated %s also shows a significant step-change around %d (p = %.4f, mean shift %+.1f%%). Change-points at the same date in naturalized series are likely climatic, not naturalization artifacts.",
        display_name, change_year, p_value, pct_shift))
    }
    cat(sprintf("    -> FLAGGED: p < %.2f — significant change-point at the %.0f%% confidence level.\n",
                HOMOGENEITY_ALPHA, (1 - HOMOGENEITY_ALPHA) * 100))
  } else {
    cat(sprintf("    -> No significant change-point detected (p >= %.2f).\n", HOMOGENEITY_ALPHA))
  }
  
  return(list(
    station_id  = station_id,
    role        = role,
    tested      = TRUE,
    n_years     = nrow(annual),
    first_year  = min(annual$year),
    last_year   = max(annual$year),
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
# SECTION 12: PROCESS ALL STATIONS — PIPELINE B
#
# Regulated stations: naturalized flow. 08JE001: observed flow.
# Steps: load → completeness → gap-fill → (Pettitt) → monthly → SSI → events
# ============================================================================
set_pipeline_tag("Pipeline B")

all_processed_data <- list()   # gap-filled daily data
all_completeness   <- list()   # completeness check results
all_thresholds     <- list()   # Pipeline A — DISABLED (kept for compatibility)
all_daily_droughts <- list()   # Pipeline A — DISABLED (kept for compatibility)
all_ssi_results    <- list()   # Pipeline B SSI
all_ssi_droughts   <- list()   # Pipeline B events
all_homogeneity    <- list()   # Pettitt results (naturalized + control)

for (i in 1:nrow(stations)) {
  station_id   <- stations$StationID[i]
  display_name <- get_station_display(station_id)
  data_path    <- get_data_path_for_station(station_id, pipeline = "B")
  
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("[PIPELINE B] PROCESSING STATION %d/%d: %s\n", i, nrow(stations), display_name))
  cat(sprintf("%s\n", strrep("=", 60)))
  
  # ── [1/6] Load ─────────────────────────────────────────────────────────────
  cat(sprintf("\n[1/6] Loading discharge data (source: %s)...\n", data_path))
  raw_data <- load_discharge_data(station_id, data_path = data_path)
  if (is.null(raw_data)) {
    cat(sprintf("  SKIPPED: no data file found for %s\n", display_name))
    next
  }
  
  stations$StartYear[stations$StationID == station_id] <- min(raw_data$year)
  stations$EndYear[stations$StationID == station_id]   <- max(raw_data$year)
  
  # ── [2/6] Completeness check ───────────────────────────────────────────────
  cat("\n[2/6] Checking data completeness...\n")
  completeness <- check_data_completeness(raw_data, station_id)
  all_completeness[[station_id]] <- completeness
  
  if (!completeness$meets_threshold) {
    cat(sprintf("  SKIPPED: %.2f%% < %.0f%% completeness threshold\n",
                completeness$completeness, MIN_COMPLETENESS_PCT))
    next
  }
  
  # ── [3/6] Gap fill ─────────────────────────────────────────────────────────
  cat(sprintf("\n[3/6] Filling gaps (Fritsch-Carlson PCHIP, max gap = %d days)...\n",
              MAX_FILL_DAYS))
  data_filled <- fill_missing_data(raw_data)
  all_processed_data[[station_id]] <- data_filled
  
  # ── [3b/6] Pettitt screening (naturalized stations) ───────────────────────
  if (station_id %in% REGULATED_STATIONS) {
    all_homogeneity[[station_id]] <- check_naturalized_homogeneity(
      data_filled, station_id, role = "naturalized")
  }
  
  # ── [4a/6] Pipeline A: variable threshold — DISABLED ──────────────────────
  # thresholds <- calculate_variable_threshold(data_filled, station_id = station_id)
  # all_thresholds[[station_id]] <- thresholds
  # ── [5a/6] Pipeline A: daily drought identification — DISABLED ────────────
  # daily_droughts <- identify_droughts(data_filled, thresholds, station_id = station_id)
  # all_daily_droughts[[station_id]] <- daily_droughts
  
  # ── [4/6] Aggregate to monthly ─────────────────────────────────────────────
  cat("\n[4/6] Pipeline B — aggregating to monthly...\n")
  monthly_data <- aggregate_to_monthly(data_filled)
  
  # ── [5/6] SSI calculation ─────────────────────────────────────────────────
  cat("\n[5/6] Pipeline B — fitting distributions & calculating SSI...\n")
  ssi_result <- calculate_ssi(monthly_data, station_id = station_id)
  if (!is.null(ssi_result)) {
    all_ssi_results[[station_id]] <- ssi_result
    
    # ── [6/6] SSI drought identification ─────────────────────────────────────
    cat(sprintf("\n[6/6] Pipeline B — identifying SSI drought events (SSI < %.2f)...\n",
                SSI_DROUGHT_THRESHOLD))
    ssi_droughts <- identify_ssi_droughts(ssi_result, station_id = station_id)
    all_ssi_droughts[[station_id]] <- ssi_droughts
  } else {
    cat("  [Pipeline B] SSI skipped — insufficient monthly data\n")
  }
  
  cat(sprintf("\n  [Pipeline B] Station %s complete.\n", station_id))
}

# ── Pettitt control test on unregulated station(s), same years as naturalized ──
nat_starts <- stations$StartYear[stations$StationID %in% REGULATED_STATIONS]
control_start <- if (is.finite(HOMOGENEITY_CONTROL_START_YEAR)) {
  HOMOGENEITY_CONTROL_START_YEAR
} else if (any(is.finite(nat_starts))) {
  min(nat_starts, na.rm = TRUE)
} else {
  NA_real_
}

for (sid in intersect(HOMOGENEITY_CONTROL_STATIONS, names(all_processed_data))) {
  cat(sprintf("\n%s\n[HOMOGENEITY CONTROL] %s (years >= %s)\n%s\n",
              strrep("=", 60), get_station_display(sid),
              ifelse(is.finite(control_start), as.character(control_start), "full record"),
              strrep("=", 60)))
  all_homogeneity[[sid]] <- check_naturalized_homogeneity(
    all_processed_data[[sid]], sid, role = "control", start_year = control_start)
}

# ============================================================================
# SECTION 12b: PROCESS ALL STATIONS — PIPELINE C
#
# Same SSI method as Pipeline B, on the ORIGINAL (observed WSC) record:
#   Regulated stations (08JA017, 08JC001, 08JC002) — observed flow from
#   INPUT_DATA_DIR (<ID>_WSC_DISCHARGE.csv), restricted to
#   year >= POST_REGULATION_START_YEAR (no naturalization).
#   08JE001 — observed flow (same series as Pipeline B).
# These are real observed records (no reconstruction), so no Pettitt
# screening is run here.
# ============================================================================
set_pipeline_tag("Pipeline C")

all_processed_data_C <- list()   # gap-filled daily data (Pipeline C)
all_completeness_C   <- list()   # completeness check results (Pipeline C)
all_ssi_results_C    <- list()   # Pipeline C SSI results
all_ssi_droughts_C   <- list()   # Pipeline C SSI drought events

for (i in 1:nrow(stations_C)) {
  station_id   <- stations_C$StationID[i]
  display_name <- get_station_display(station_id)
  data_path    <- get_data_path_for_station(station_id, pipeline = "C")
  
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("[PIPELINE C] PROCESSING STATION %d/%d: %s\n", i, nrow(stations_C), display_name))
  cat(sprintf("%s\n", strrep("=", 60)))
  
  # ── [1/6] Load ─────────────────────────────────────────────────────────────
  cat(sprintf("\n[1/6] Loading discharge data (source: %s)...\n", data_path))
  raw_data <- load_discharge_data(station_id, data_path = data_path)
  if (is.null(raw_data)) {
    cat(sprintf("  SKIPPED: no data file found for %s\n", display_name))
    next
  }
  
  # Regulated stations: observed record, post-regulation period only
  if (station_id %in% REGULATED_STATIONS) {
    n_before <- nrow(raw_data)
    raw_data <- raw_data[raw_data$year >= POST_REGULATION_START_YEAR, ]
    cat(sprintf("  [Pipeline C] %s: restricted to post-regulation period (year >= %d): %d -> %d records\n",
                station_id, POST_REGULATION_START_YEAR, n_before, nrow(raw_data)))
    if (nrow(raw_data) == 0) {
      cat(sprintf("  SKIPPED: no post-regulation records remain for %s\n", display_name))
      next
    }
  }
  
  stations_C$StartYear[stations_C$StationID == station_id] <- min(raw_data$year)
  stations_C$EndYear[stations_C$StationID == station_id]   <- max(raw_data$year)
  
  # ── [2/6] Completeness check ───────────────────────────────────────────────
  cat("\n[2/6] Checking data completeness...\n")
  completeness <- check_data_completeness(raw_data, station_id)
  all_completeness_C[[station_id]] <- completeness
  
  if (!completeness$meets_threshold) {
    cat(sprintf("  SKIPPED: %.2f%% < %.0f%% completeness threshold\n",
                completeness$completeness, MIN_COMPLETENESS_PCT))
    next
  }
  
  # ── [3/6] Gap fill ─────────────────────────────────────────────────────────
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
    cat(sprintf("\n[6/6] Pipeline C — identifying SSI drought events (SSI < %.2f)...\n",
                SSI_DROUGHT_THRESHOLD))
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

# Shared / Pipeline B daily data
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

# Pipeline B (naturalized regulated stations + observed 08JE001)
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
cat("  Pettitt change-point results (naturalized + control) saved\n")

# Pipeline C (observed post-regulation regulated stations + observed 08JE001)
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
cat("  Pipeline C (SSI monthly, observed/post-regulation) outputs saved\n")

common_metadata <- list(
  input_data_dir              = INPUT_DATA_DIR,
  output_dir                  = MAIN_OUTPUT_DIR,
  completeness_denominator    = "month_span (first_day_of_first_month to last_day_of_last_month)",
  min_completeness_pct        = MIN_COMPLETENESS_PCT,
  max_fill_days               = MAX_FILL_DAYS,
  interpolation_method        = "monoH.FC (Fritsch-Carlson PCHIP)",
  min_valid_days_per_month    = MIN_VALID_DAYS_PER_MONTH,
  min_observed_days_per_month = MIN_OBSERVED_DAYS_PER_MONTH,
  monthly_completeness_rule   = MONTHLY_RULE_LABEL,
  ssi_drought_threshold       = SSI_DROUGHT_THRESHOLD,
  ssi_recovery_threshold      = SSI_RECOVERY_THRESHOLD,
  distributions_tested        = DISTRIBUTIONS,
  distribution_selection      = "minimum L-kurtosis distance |tau4_sample - tau4_fitted| (lmomco), per calendar month",
  support_rule                = sprintf("strict (all fitting values inside fitted support) preferred; relaxed fallback <= %d out-of-support values (clamped); TOLERANT_SELECTION = %s",
                                        MAX_OUT_OF_SUPPORT, TOLERANT_SELECTION),
  max_out_of_support          = MAX_OUT_OF_SUPPORT,
  tolerant_selection          = TOLERANT_SELECTION,
  min_months_for_fit          = MIN_MONTHS_FOR_FIT,
  min_years_per_month         = MIN_YEARS_PER_MONTH,
  min_drought_duration_months = get_event_min_duration(SSI_EVENT_SCALE),
  spatial_support             = "Station-specific hydrometric records; SSI values are not spatially averaged across gauges",
  spatial_support_metadata    = "reports/station_spatial_support.csv",
  processed_date              = Sys.Date(),
  script_version              = SCRIPT_VERSION
)

all_results_C <- list(
  stations       = stations_C,
  completeness   = all_completeness_C,
  processed_data = all_processed_data_C,
  ssi_results    = all_ssi_results_C,
  ssi_droughts   = all_ssi_droughts_C,
  metadata = c(common_metadata, list(
    pipeline_c                 = sprintf("Peña-Angulo et al. (2022) SSI monthly, ORIGINAL (observed) flow for %s + observed 08JE001",
                                         paste(REGULATED_STATIONS, collapse = ", ")),
    post_regulation_stations   = REGULATED_STATIONS,
    post_regulation_note       = sprintf(
      "%s: original/observed WSC record (NOT naturalized), restricted to year >= %d; pre-regulation years dropped",
      paste(REGULATED_STATIONS, collapse = ", "), POST_REGULATION_START_YEAR),
    post_regulation_start_year = POST_REGULATION_START_YEAR
  ))
)
saveRDS(all_results_C, file.path(MAIN_OUTPUT_DIR, "all_results_pipeline_C.rds"))
cat("  Master file saved: all_results_pipeline_C.rds\n")

all_results <- list(
  stations          = stations,
  completeness      = all_completeness,
  processed_data    = all_processed_data,
  ssi_results       = all_ssi_results,
  ssi_droughts      = all_ssi_droughts,
  homogeneity_check = all_homogeneity,
  metadata = c(common_metadata, list(
    naturalized_data_dir   = NATURALIZED_DATA_DIR,
    pipeline_b             = sprintf("Peña-Angulo et al. (2022) SSI monthly, NATURALIZED flow for %s + observed 08JE001",
                                     paste(REGULATED_STATIONS, collapse = ", ")),
    naturalized_stations   = REGULATED_STATIONS,
    naturalized_note       = "08JC001/08JC002: dat_natural_<ID>.csv = Stephen's verified naturalized series to 2021-12-15, then water-balance reconstruction (WSC 08JA023 level + stage-storage, 08JA013 + 08FE002 outflows; routing per Kokoszka 2026) to 2025-12-31; Oct-Dec 2025 inputs provisional. 08JA017: see NATURALIZED_FILES",
    naturalized_files      = NATURALIZED_FILES,
    naturalized_file_paths = NATURALIZED_FILE_PATHS,
    naturalization_met_predictors = NATURALIZATION_MET_PREDICTORS,
    naturalization_note    = NATURALIZATION_NOTE,
    analysis_end_date      = ANALYSIS_END_DATE,
    provisional_handling   = PROVISIONAL_HANDLING,
    homogeneity_test       = "Pettitt test (trend::pettitt.test) on annual-mean discharge of near-complete years; naturalized stations plus unregulated control over the same years",
    homogeneity_alpha      = HOMOGENEITY_ALPHA,
    homogeneity_control    = HOMOGENEITY_CONTROL_STATIONS,
    homogeneity_control_start_year = control_start
  ))
)
saveRDS(all_results, file.path(MAIN_OUTPUT_DIR, "all_results_dual_pipeline.rds"))
cat("  Master file saved: all_results_dual_pipeline.rds\n")

# ============================================================================
# SECTION 13.5: ADDITIONAL OUTPUTS
# ============================================================================
cat(sprintf("\n%s\nGENERATING ADDITIONAL OUTPUTS...\n%s\n",
            strrep("=", 60), strrep("=", 60)))

create_station_completeness_plot(stations, all_completeness, MAIN_OUTPUT_DIR,
                                 file_name = "station_completeness_pipeline_B.png",
                                 plot_label = "Pipeline B Station Data Completeness")
if (length(all_completeness_C) > 0L) {
  create_station_completeness_plot(stations_C, all_completeness_C, MAIN_OUTPUT_DIR,
                                   file_name = "station_completeness_pipeline_C.png",
                                   plot_label = "Pipeline C Station Data Completeness")
}

# ============================================================================
# SECTION 13.5b: SSI TIME SERIES PLOTS
# ============================================================================
# One PNG per station per pipeline: blue SSI line, red fill for detected
# drought events (from the shared event table), reference lines at ±1/±1.5/±2,
# soft background bands, and the clamp floor.
plot_ssi_timeseries <- function(monthly_data, station_id, station_name,
                                pipeline_label, out_dir,
                                onset = SSI_DROUGHT_THRESHOLD,
                                drought_events = NULL) {
  
  if (is.null(monthly_data) || !"ssi" %in% names(monthly_data)) {
    cat(sprintf("  [PLOT] Skipping %s (%s): no SSI data available\n",
                station_id, pipeline_label))
    return(invisible(NULL))
  }
  
  df <- data.frame(
    date = as.Date(monthly_data$date),
    ssi  = as.numeric(monthly_data$ssi)
  )
  df <- df[order(df$date), , drop = FALSE]
  finite_ssi <- df$ssi[is.finite(df$ssi)]
  if (!length(finite_ssi)) {
    cat(sprintf("  [PLOT] Skipping %s (%s): no finite SSI values\n",
                station_id, pipeline_label))
    return(invisible(NULL))
  }
  
  n_events <- if (!is.null(drought_events) && is.data.frame(drought_events)) {
    nrow(drought_events)
  } else NA_integer_
  
  # Shade exactly the months belonging to detected events; NA months break runs.
  df$below <- FALSE
  if (!is.null(drought_events) && is.data.frame(drought_events) && nrow(drought_events)) {
    for (k in seq_len(nrow(drought_events)))
      df$below[df$date >= as.Date(drought_events$start_date[k]) &
                 df$date <= as.Date(drought_events$end_date[k])] <- TRUE
  }
  df$below <- df$below & is.finite(df$ssi)
  month_id  <- as.integer(format(df$date, "%Y")) * 12L + as.integer(format(df$date, "%m"))
  new_group <- c(TRUE, diff(as.integer(df$below)) != 0L | diff(month_id) != 1L)
  df$grp    <- cumsum(new_group)
  
  # Clamp floor, from the same probability clamp used in the SSI transform
  ssi_clamp_p  <- 1e-6
  ssi_clamp_lo <- qnorm(ssi_clamp_p)
  
  y_lo <- min(-2.2, floor(min(finite_ssi) * 10) / 10,
              floor(ssi_clamp_lo * 10) / 10 - 0.1)
  y_hi <- max(2.2, ceiling(max(finite_ssi) * 10) / 10)
  
  p <- ggplot(df, aes(x = date, y = ssi)) +
    annotate("rect", xmin = DATE_NEG_INF, xmax = DATE_POS_INF, ymin = 1,    ymax = Inf,
             fill = "#4472C4", alpha = 0.10) +
    annotate("rect", xmin = DATE_NEG_INF, xmax = DATE_POS_INF, ymin = 0,    ymax = 1,
             fill = "#70AD47", alpha = 0.08) +
    annotate("rect", xmin = DATE_NEG_INF, xmax = DATE_POS_INF, ymin = -1,   ymax = 0,
             fill = "#70AD47", alpha = 0.05) +
    annotate("rect", xmin = DATE_NEG_INF, xmax = DATE_POS_INF, ymin = -Inf, ymax = -1,
             fill = "#C00000", alpha = 0.08) +
    geom_hline(yintercept = 0,            color = "black",  linewidth = 0.6) +
    geom_hline(yintercept = c(-1, 1),     color = "grey30", linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = c(-1.5, 1.5), color = "grey55", linetype = "dotted", linewidth = 0.4) +
    geom_hline(yintercept = c(-2, 2),     color = "grey55", linetype = "dotted", linewidth = 0.4) +
    geom_hline(yintercept = ssi_clamp_lo, color = "grey35", linetype = "dashed", linewidth = 0.45) +
    annotate("text", x = min(df$date), y = ssi_clamp_lo, hjust = 0, vjust = -0.6,
             label = sprintf("clamp floor (p = %.0e \u2192 SSI = %.3f)",
                             ssi_clamp_p, ssi_clamp_lo),
             color = "grey35", size = 3.2, fontface = "italic") +
    geom_ribbon(data = subset(df, below),
                aes(ymin = ssi, ymax = 0, group = grp),
                fill = "#C0504D", alpha = 0.45) +
    geom_line(color = "#1F77B4", linewidth = 0.45, na.rm = TRUE) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.015))) +
    scale_y_continuous(
      breaks = sort(unique(c(seq(-4, 4, 1), round(ssi_clamp_lo, 2)))),
      labels = function(v) ifelse(v == round(ssi_clamp_lo, 2),
                                  sprintf("%.2f", v), as.character(v))
    ) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    labs(
      title    = sprintf("SSI Time Series  |  %s - %s", station_id, station_name),
      subtitle = sprintf("%s  |  %s drought events  (onset < %.2f, end >= %.2f, min %d months)",
                         pipeline_label, ifelse(is.na(n_events), "NA", as.character(n_events)), onset,
                         SSI_RECOVERY_THRESHOLD, get_event_min_duration(SSI_EVENT_SCALE)),
      x = "Year", y = "SSI Index Value"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey85"),
      panel.grid.major.y = element_blank(),
      plot.title         = element_text(face = "bold", size = 15),
      plot.subtitle      = element_text(color = "grey40", size = 11),
      axis.text.x        = element_text(angle = 30, hjust = 1),
      plot.margin        = margin(10, 16, 10, 10)
    )
  
  pipeline_tag <- gsub("[^A-Za-z0-9]+", "", pipeline_label)
  fname <- file.path(out_dir, "figures",
                     sprintf("ssi_%s_%s_timeseries.png", station_id, pipeline_tag))
  ggsave(fname, p, width = 11, height = 5.5, dpi = 200, bg = "white")
  cat(sprintf("  [PLOT] SSI time series saved: %s  (%s events)\n",
              basename(fname), ifelse(is.na(n_events), "NA", as.character(n_events))))
  invisible(p)
}

cat(sprintf("\n[PLOTS] Pipeline B (naturalized %s + observed 08JE001)...\n",
            paste(REGULATED_STATIONS, collapse = ", ")))
for (sid in names(all_ssi_results)) {
  plot_ssi_timeseries(
    monthly_data   = all_ssi_results[[sid]]$monthly_data,
    station_id     = sid,
    station_name   = get_station_name(sid),
    pipeline_label = if (sid %in% REGULATED_STATIONS) "Pipeline B (naturalized flow)"
    else "Pipeline B (observed)",
    out_dir        = MAIN_OUTPUT_DIR,
    drought_events = all_ssi_droughts[[sid]]$drought_events
  )
}

cat(sprintf("\n[PLOTS] Pipeline C (observed post-regulation %s + observed 08JE001)...\n",
            paste(REGULATED_STATIONS, collapse = ", ")))
for (sid in names(all_ssi_results_C)) {
  plot_ssi_timeseries(
    monthly_data   = all_ssi_results_C[[sid]]$monthly_data,
    station_id     = sid,
    station_name   = get_station_name(sid),
    pipeline_label = if (sid %in% REGULATED_STATIONS)
      sprintf("Pipeline C (observed, post-reg >= %d)", POST_REGULATION_START_YEAR)
    else "Pipeline C (observed)",
    out_dir        = MAIN_OUTPUT_DIR,
    drought_events = all_ssi_droughts_C[[sid]]$drought_events
  )
}

# Excel summary workbook
create_excel_summary(stations, all_completeness,
                     all_daily_droughts, all_ssi_droughts,
                     all_ssi_results, MAIN_OUTPUT_DIR,
                     stations_C = stations_C, completeness_C = all_completeness_C,
                     ssi_droughts_C = all_ssi_droughts_C, ssi_results_C = all_ssi_results_C,
                     homogeneity = all_homogeneity)

# Detailed text report
params <- list(
  script_version       = SCRIPT_VERSION,
  min_completeness_pct = MIN_COMPLETENESS_PCT
)
create_text_report(stations, all_completeness,
                   all_daily_droughts, all_ssi_droughts,
                   all_ssi_results, MAIN_OUTPUT_DIR, params,
                   completeness_C = all_completeness_C,
                   ssi_droughts_C = all_ssi_droughts_C,
                   ssi_results_C = all_ssi_results_C,
                   post_regulation_start_year = POST_REGULATION_START_YEAR,
                   homogeneity = all_homogeneity)

# ============================================================================
# SECTION 14: PIPELINE B SUMMARY
# ============================================================================
cat(sprintf("\n%s\nPIPELINE B SUMMARY (naturalized %s + observed 08JE001)\n%s\n\n",
            strrep("=", 60), paste(REGULATED_STATIONS, collapse = ", "), strrep("=", 60)))

summary_rows <- list()
for (sid in names(all_completeness)) {
  comp  <- all_completeness[[sid]]
  sname <- if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN"
  base <- data.frame(
    StationID        = sid,
    StationName      = sname,
    StationDisplay   = paste0(sid, " - ", sname),
    Period           = sprintf("%d-%d", comp$start_year, comp$end_year),
    Completeness_pct = round(comp$completeness, 1),
    Status           = ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED"),
    DataSource       = ifelse(sid %in% REGULATED_STATIONS, "Naturalized flow", "Observed"),
    NaturalizedFile  = naturalized_file_of(sid),
    NaturalizationMetPredictors = if (sid %in% REGULATED_STATIONS)
      ifelse(isTRUE(NATURALIZATION_MET_PREDICTORS), "YES", "NONE") else NA_character_,
    MonthlyRule      = MONTHLY_RULE_LABEL,
    stringsAsFactors = FALSE
  )
  if (!comp$meets_threshold) {
    base$B_Best_Distribution    <- NA_character_
    base$B_SSI_Coverage         <- NA_character_
    base$B_Monthly_Completeness <- NA_character_
    base$B_Provisional_Months   <- NA_character_
    base$B_Drought_Events       <- NA_integer_
    base$B_Total_Months         <- NA_integer_
    base$B_Avg_Duration_Months  <- NA_real_
    base$B_Avg_Severity         <- NA_real_
  } else {
    ev_b <- all_ssi_droughts[[sid]]$drought_events
    base$B_Best_Distribution    <- if (!is.null(all_ssi_results[[sid]]))
      all_ssi_results[[sid]]$best_distribution_str else NA_character_
    base$B_SSI_Coverage         <- ssi_coverage_str(all_ssi_results[[sid]])
    base$B_Monthly_Completeness <- monthly_completeness_str(all_ssi_results[[sid]])
    base$B_Provisional_Months   <- provisional_months_str(all_ssi_results[[sid]])
    base$B_Drought_Events       <- if (!is.null(ev_b)) nrow(ev_b) else NA_integer_
    base$B_Total_Months         <- if (!is.null(ev_b)) sum(ev_b$duration_months) else NA_integer_
    base$B_Avg_Duration_Months  <- if (!is.null(ev_b) && nrow(ev_b) > 0)
      round(mean(ev_b$duration_months), 1) else NA_real_
    base$B_Avg_Severity         <- if (!is.null(ev_b) && nrow(ev_b) > 0)
      round(mean(ev_b$severity), 2) else NA_real_
  }
  
  base <- cbind(base, homogeneity_fields(all_homogeneity[[sid]]))
  summary_rows[[sid]] <- base
}

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table,
          file.path(MAIN_OUTPUT_DIR, "reports", "station_summary_dual_pipeline.csv"),
          row.names = FALSE)
print(summary_table, row.names = FALSE)

# ============================================================================
# SECTION 14b: PIPELINE C SUMMARY
# ============================================================================
cat(sprintf("\n%s\nPIPELINE C SUMMARY (observed %s post-regulation only + observed 08JE001)\n%s\n\n",
            strrep("=", 60), paste(REGULATED_STATIONS, collapse = ", "), strrep("=", 60)))

summary_rows_C <- list()
for (sid in names(all_completeness_C)) {
  comp  <- all_completeness_C[[sid]]
  sname <- if (!is.null(comp$station_name)) comp$station_name else "UNKNOWN"
  base <- data.frame(
    StationID        = sid,
    StationName      = sname,
    StationDisplay   = paste0(sid, " - ", sname),
    Period           = sprintf("%d-%d", comp$start_year, comp$end_year),
    Completeness_pct = round(comp$completeness, 1),
    Status           = ifelse(comp$meets_threshold, "INCLUDED", "EXCLUDED"),
    DataSource       = ifelse(sid %in% REGULATED_STATIONS,
                              sprintf("Observed, post-regulation only (year >= %d)", POST_REGULATION_START_YEAR),
                              "Observed"),
    MonthlyRule      = MONTHLY_RULE_LABEL,
    stringsAsFactors = FALSE
  )
  if (!comp$meets_threshold) {
    base$Best_Distribution    <- NA_character_
    base$SSI_Coverage         <- NA_character_
    base$Monthly_Completeness <- NA_character_
    base$Drought_Events       <- NA_integer_
    base$Total_Months         <- NA_integer_
    base$Avg_Duration_Months  <- NA_real_
    base$Avg_Severity         <- NA_real_
  } else {
    ev_c <- all_ssi_droughts_C[[sid]]$drought_events
    base$Best_Distribution    <- if (!is.null(all_ssi_results_C[[sid]]))
      all_ssi_results_C[[sid]]$best_distribution_str else NA_character_
    base$SSI_Coverage         <- ssi_coverage_str(all_ssi_results_C[[sid]])
    base$Monthly_Completeness <- monthly_completeness_str(all_ssi_results_C[[sid]])
    base$Drought_Events       <- if (!is.null(ev_c)) nrow(ev_c) else NA_integer_
    base$Total_Months         <- if (!is.null(ev_c)) sum(ev_c$duration_months) else NA_integer_
    base$Avg_Duration_Months  <- if (!is.null(ev_c) && nrow(ev_c) > 0)
      round(mean(ev_c$duration_months), 1) else NA_real_
    base$Avg_Severity         <- if (!is.null(ev_c) && nrow(ev_c) > 0)
      round(mean(ev_c$severity), 2) else NA_real_
  }
  
  summary_rows_C[[sid]] <- base
}

if (length(summary_rows_C) > 0L) {
  summary_table_C <- do.call(rbind, summary_rows_C)
  write.csv(summary_table_C,
            file.path(MAIN_OUTPUT_DIR, "reports", "station_summary_pipeline_C.csv"),
            row.names = FALSE)
  print(summary_table_C, row.names = FALSE)
} else {
  cat("  No Pipeline C stations were processed.\n")
}

# Spatial-support metadata: SSI outputs remain station-specific
write.csv(STATION_SPATIAL_SUPPORT,
          file.path(MAIN_OUTPUT_DIR, "reports", "station_spatial_support.csv"),
          row.names = FALSE, na = "")

nat_list <- paste(REGULATED_STATIONS, collapse = ", ")
cat("\nOutputs written to: ", MAIN_OUTPUT_DIR, "/\n")
cat(sprintf("  Monthly completeness rule: %s\n", MONTHLY_RULE_LABEL))
cat("  shared/processed_data/            — gap-filled daily data (Pipeline B + Pipeline C; 'interpolated' flag)\n")
cat(sprintf("  ssi/monthly_data/                 — Pipeline B: monthly SSI (naturalized %s; observed 08JE001)\n", nat_list))
cat("  ssi/drought_events/               — Pipeline B: SSI drought events\n")
cat("  ssi/naturalized_homogeneity_check.rds — Pettitt results (naturalized stations + 08JE001 control)\n")
cat(sprintf("  ssi_C/monthly_data/               — Pipeline C: monthly SSI (observed post-regulation %s; 08JE001)\n", nat_list))
cat("  ssi_C/drought_events/             — Pipeline C: SSI drought events\n")
cat("  figures/                          — station_completeness_pipeline_B.png / _C.png\n")
cat("  figures/                          — ssi_{StationID}_{PipelineTag}_timeseries.png\n")
cat("  reports/                          — station_spatial_support.csv (no basin averaging)\n")
cat("                                    — station_summary_dual_pipeline.csv (Pipeline B)\n")
cat("                                    — station_summary_pipeline_C.csv\n")
cat("                                    — nechako_streamflow_drought_summary.xlsx\n")
cat("                                    — nechako_streamflow_drought_report.txt\n")
cat(sprintf("  all_results_dual_pipeline.rds     — Pipeline B master file (%d stations)\n", length(all_completeness)))
cat(sprintf("  all_results_pipeline_C.rds        — Pipeline C master file (%d stations)\n", length(all_completeness_C)))