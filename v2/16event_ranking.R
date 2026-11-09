####################################################################################
# 16event_ranking.R  ·  DROUGHT EVENT RANKING & PET-FORMULATION ROBUSTNESS
#
# EVENT DEFINITION — framework-wide, from DROUGHT_ANALYSIS_utils.R (section F);
# identical to upstream Scripts. Edit thresholds/durations ONLY there.
#   Onset       : index <  DROUGHT_ONSET  (currently -0.5)
#   Termination : index >= DROUGHT_END    (currently -0.5; no hysteresis because equal to onset)
#   Intensity   : I = mean(DROUGHT_END − x_t) over event months      [≥ 0]
#   Severity    : S = I × D  (Sheffield & Wood 2008 style)
#   Min duration: scale-dependent, EVENT_MIN_DURATION_BY_SCALE
#                 (currently 1→2, 3→3, 6→4, 12→6, 24→12, 36→12 months)
#   Missing months terminate an event (gaps are never bridged)
#   Duration classes: D3-6 (Short-term), D7-12 (Medium-term), D13+ (Long-term);
#                 D1-2 can occur only where the minimum duration is <= 2; valid
#                 events in this class are excluded only from D3+ RP comparisons
#   Last monthly available data is for Dec 2025
#   Basin series : fractional-coverage area-weighted mean over the 7-subbasin
#                  union (cell area x coverage), reported only when finite cells
#                  hold >= BASIN_MIN_VALID_AREA_COVERAGE (0.20) of the basin weight.
#
#   This framework-wide event definition is used consistently
#   for SPI, SPEI-PM, SWEI, SSPI, SSMI and SSI-SM. Snow indices are processed
#   only at their defined native scales (SWEI-03 and SSPI-01). The ranking uses only the manuscript
#   formulations retained in this script; no detrended/counterfactual branch
#   is used in this script.
#
#   References:
#     Sheffield & Wood (2008) Clim Dyn 31:79-105 — duration classes, I×D severity
#     Yevjevich (1967) — threshold-level run method (theoretical basis)
#
# MODIFICATIONS :
# STRUCTURE (execution order; labels are not numerically ordered):
#   Helpers: build_sw08_catalog(), sw08_event_loop(), classify_duration()
#   Part 1  event catalogs from producer basin series
#   Part 2  ranking metrics            Part 3 / 3b  ranked catalog, duration-class RP
#   Part 4  presentation figures 1-4   Part 4b  manuscript Fig 3
#   Part 4c canonical Kendall engine run_kendall_rp()
#   Part 4d multi-index Kendall comparison
#   Part 5 / 5b  spatial and pixel-rank maps (NetCDF / RDS)
#   Part 6  Excel summary              Part 4e  Table 1a
#   Part 4f recent-event bubble plot   Part 4g  Kendall RP, all scales
#   MS2 export block                   Regional extension + Part 4h
# Outputs: see the OUTPUT FILES listing printed at the end of the script.
####################################################################################

# WD_PATH / UTILS_FILE are set from getwd() at the top of this script.
rm(list = ls())

# ---- Working directory (framework rule) -------------------------------------
# Run every script from the project root: set the working directory there
# BEFORE running (setwd(), or open the RStudio project). All paths below are
# derived from getwd(); no path is hard-coded.
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")

source(UTILS_FILE)
# Shared analysis contracts are defined in DROUGHT_ANALYSIS_utils.R.
required_utils <- c("SWEI_SCALE", "SSPI_SCALE", "load_basin_vect", "ANALYSIS_PERIOD_LABEL")
missing_utils <- required_utils[!vapply(required_utils, exists, logical(1), inherits = TRUE)]
if (length(missing_utils))
  stop("DROUGHT_ANALYSIS_utils.R lacks: ", paste(missing_utils, collapse = ", "),
       "\nApply the utils patches (U1-U7) before running this script.")

utils_load_packages(c("terra", "ggplot2", "dplyr", "data.table", "lubridate",
                      "openxlsx", "patchwork", "scales", "sf", "viridis",
                      "gridExtra",      # summary table grobs
                      "copula",         # bivariate copula fitting (Part 4c)
                      "fitdistrplus")) # marginal distribution fitting (Part 4c)


####################################################################################
# CONFIGURATION
####################################################################################
RANKING_DIR <- file.path(WD_PATH, "temporal_drought", "event_ranking")
dir.create(RANKING_DIR, showWarnings = FALSE, recursive = TRUE)

# Analysis period
HISTORICAL_START <- ANALYSIS_PERIOD_START
HISTORICAL_END   <- ANALYSIS_PERIOD_END
RECENT_DROUGHT_START <- ANALYSIS_RECENT_START
RECENT_DROUGHT_END   <- ANALYSIS_RECENT_END
RECENT_YEAR_START <- as.integer(format(RECENT_DROUGHT_START, "%Y"))
RECENT_YEAR_END   <- as.integer(format(RECENT_DROUGHT_END, "%Y"))
RECENT_LABEL <- sprintf("%d–%d", RECENT_YEAR_START, RECENT_YEAR_END)
RECENT_TAG   <- sprintf("%d-%d", RECENT_YEAR_START, RECENT_YEAR_END)
RECENT_TAG_UNDERSCORE <- sprintf("%d_%d", RECENT_YEAR_START, RECENT_YEAR_END)
RECENT_GROUP_LABEL <- paste0("Events overlapping ", RECENT_LABEL)
# Use unambiguous index names throughout headline figures.
event_index_label <- function(label) {
  scale <- as.integer(sub(".*_([0-9]+)$", "\\1", label))
  idx <- sub("_[0-9]+$", "", label)
  prefix <- switch(idx, spei = "SPEI-PM",
                   spi = "SPI", swei = "SWEI", sspi = "SSPI",
                   ssmi = "SSMI", ssi_sm = "SSI-SM", toupper(idx))
  sprintf("%s-%d", prefix, scale)
}

stopifnot(
  identical(as.integer(MANUSCRIPT_SCALES), c(1L, 3L, 6L, 12L, 24L, 36L)),
  identical(as.integer(EVENT_RANKING_SCALES), as.integer(MANUSCRIPT_SCALES))
)


# ── Drought event identification parameters ───────────────────────────────────
# Threshold: −0.5 corresponds to the 31st percentile of the standard normal,
#   a relatively inclusive drought onset threshold (Lloyd-Hughes &
#   Saunders 2002; Spinoni et al. 2014, 2015). It flags more months
#   than the −1.0 standardized anomaly reference (16th percentile)
#   and corresponds to the lower tercile of the precipitation distribution.
#   Note: the SW08 10th-percentile threshold translates to −1.28 on a standard
#   normal scale, which is distinct from both −1.0 and −0.5.
DROUGHT_ONSET_EVENT <- DROUGHT_ONSET   # shared utility onset threshold
DROUGHT_END_EVENT   <- DROUGHT_END     # shared utility termination threshold
if (DROUGHT_END_EVENT < DROUGHT_ONSET_EVENT)
  stop("DROUGHT_END must be >= DROUGHT_ONSET for a well-defined threshold-run event.")
cat(sprintf("Event rule: onset < %.2f, termination >= %.2f; %s. -1.0 is a separate reference threshold.\n",
            DROUGHT_ONSET_EVENT, DROUGHT_END_EVENT,
            if (DROUGHT_END_EVENT == DROUGHT_ONSET_EVENT) "NO hysteresis" else "hysteresis enabled"))
# SWEI_SCALE (3) and SSPI_SCALE (1) come from DROUGHT_ANALYSIS_utils.R.
SOIL_CFG_PATH <- file.path(TREND_DIR, "script13_analysis_config.rds")
if (!file.exists(SOIL_CFG_PATH))
  stop("Missing required Script 13 configuration: ", SOIL_CFG_PATH,
       "; run 13trend_test_ALL.R before Script 16.")
SOIL_CFG <- readRDS(SOIL_CFG_PATH)
if (length(SOIL_CFG$soil_layer) != 1L || is.na(SOIL_CFG$soil_layer))
  stop("Missing or invalid soil_layer in ", SOIL_CFG_PATH)
SOIL_LAYER <- as.character(SOIL_CFG$soil_layer)
SOIL_SCALES <- as.integer(SOIL_CFG$soil_scales)
if (anyNA(SOIL_SCALES)) stop("Invalid soil_scales in ", SOIL_CFG_PATH)
SSMI_SEAS_DIR <- file.path(WD_PATH, "ssmi_results_seasonal", SOIL_LAYER)
SSI_SM_SEAS_DIR <- file.path(WD_PATH, "ssi_sm_results_seasonal", SOIL_LAYER)
SSPI_SEAS_DIR <- file.path(WD_PATH, "sspi_results_seasonal/")

# One headline index for every presentation and manuscript product (Figs 1-6,
# dashboard, Excel, Fig 3, Kendall RP, MS2 export). Set to "spi_03" for SPI-3.
HEADLINE_INDEX <- "spei_03"

# Resolve each index to its authoritative seasonal NetCDF directory.
get_event_index_dir <- function(index_type) {
  switch(tolower(index_type),
         spi       = SPI_SEAS_DIR,
         spei      = SPEI_SEAS_DIR,
         swei      = SWEI_SEAS_DIR,
         sspi      = SSPI_SEAS_DIR,
         ssmi      = SSMI_SEAS_DIR,
         ssi_sm    = SSI_SM_SEAS_DIR,
         stop("Unknown index_type: ", index_type)
  )
}

get_event_index_files <- function(index_type, scale) {
  idx <- tolower(index_type); seas_dir <- get_event_index_dir(idx)
  if (idx == "swei")           find_swei_seasonal_files(seas_dir, idx, as.integer(scale))
  else                         find_seasonal_nc_files(seas_dir, idx, as.integer(scale))
}
# ── Event-duration rule ───────────────────────────────────────────────────────
# Scale-dependent minimum durations come from get_event_min_duration() in
# DROUGHT_ANALYSIS_utils.R (EVENT_MIN_DURATION_BY_SCALE). Runs shorter than the
# minimum are not events. Do not redefine get_event_min_duration() here.
MIN_DURATION_MAP <- vapply(as.character(EVENT_RANKING_SCALES),
                           function(s) get_event_min_duration(as.integer(s)), integer(1))

# All supported event-catalog indices and their valid scales. SPI and SPEI-PM
# use the six manuscript scales; snow indices use only their
# native scales produced by 4snow_drought_indices_ERALand.R.
INDICES_TO_ANALYZE <- list(
  list(index = "spi",       scales = EVENT_RANKING_SCALES),
  list(index = "spei",      scales = EVENT_RANKING_SCALES),
  list(index = "swei",      scales = c(SWEI_SCALE)),
  list(index = "sspi",      scales = c(SSPI_SCALE)),
  list(index = "ssmi",      scales = SOIL_SCALES),
  list(index = "ssi_sm",    scales = SOIL_SCALES)
)

####################################################################################
# HELPER FUNCTIONS
####################################################################################

#' Classify drought duration according to Sheffield & Wood (2008)
classify_duration <- function(duration_months) {
  sapply(duration_months, function(d) {
    if (d >= 3 & d <= 6) return("D3-6 (Short-term)")
    if (d >= 7 & d <= 12) return("D7-12 (Medium-term)")
    if (d >= 13) return("D13+ (Long-term)")
    return("D1-2 (Very short)")
  })
}

#' Estimate return period from empirical frequency
# Empirical severity-rank plotting position over all event durations for a
# given index/scale, NOT an annual-exceedance probability or a fitted tail RP.
# Ties use the best competition rank (min_rank), not the count with S >= S_obs.
# Several scale-1 events can occur per year; intervals below one year are possible.
calc_return_period <- function(event_rank, n_years_record) {
  # Severity-rank recurrence proxy: T = (N+1)/m (rank 1 capped at N+1).
  (n_years_record + 1) / event_rank
}

# ── Event table on a pre-computed time series vector ─────────────────────────
#
# Thin wrapper around the shared framework engine drought_event_table()
# (DROUGHT_ANALYSIS_utils.R): threshold-run method, onset x < onset_thr;
# hysteresis exists only when end_thr > onset_thr. At equal thresholds there is none.
# termination x >= end_thr, I = mean(end_thr − x), S = I × D, scale-dependent
# minimum duration. The function name is kept for backward compatibility.
#
# Returns a data.frame with one row per qualifying event:
#   start_date, end_date, start_year, duration_months, duration_class,
#   mean_int (= I), severity (= S)
sw08_event_loop <- function(vals, date_seq, yr_seq,
                            onset_thr = DROUGHT_ONSET_EVENT,
                            min_duration,
                            end_thr = DROUGHT_END_EVENT) {
  vals <- as.numeric(vals)
  date_seq <- as.Date(date_seq)
  yr_seq <- as.integer(yr_seq)
  min_duration <- as.integer(min_duration)
  
  if (length(vals) != length(date_seq) || length(vals) != length(yr_seq))
    stop("sw08_event_loop: vals/date_seq/yr_seq lengths differ.")
  
  validate_monthly_record(date_seq, label = "SW08 event-loop calendar",
                          require_complete = TRUE)
  validate_standardized_index(vals, label = "SW08 drought-index series")
  
  ev <- drought_event_table(vals, date_seq, min_duration = min_duration,
                            onset_threshold = onset_thr, termination_threshold = end_thr)
  
  out <- if (nrow(ev)) data.frame(
    start_date      = ev$start_date,
    end_date        = ev$end_date,
    start_year      = yr_seq[ev$start_index],
    duration_months = as.integer(ev$duration_months),
    duration_class  = classify_duration(ev$duration_months),
    mean_int        = ev$mean_int,
    severity        = ev$severity,
    right_censored  = ev$right_censored,
    stringsAsFactors = FALSE
  ) else
    data.frame(start_date = as.Date(character()), end_date = as.Date(character()),
               start_year = integer(), duration_months = integer(),
               duration_class = character(), mean_int = numeric(),
               severity = numeric(), right_censored = logical(), stringsAsFactors = FALSE)
  
  validate_event_catalog(out, min_duration = min_duration,
                         label = "SW08 event catalog")
  out
}


# ── Build SW08 catalog from producer basin series for one index × scale ───────
#
# Reads the basin series computed by the producer, then applies sw08_event_loop().
# Saves the catalog as CSV for reproducibility.
#
# Returns a data.frame or NULL if the basin series cannot be read.
build_sw08_catalog <- function(index_type, scale,
                               onset_thr      = DROUGHT_ONSET_EVENT,
                               min_duration   = get_event_min_duration(scale),
                               start_year_rec = HISTORICAL_START,
                               end_year_rec   = HISTORICAL_END) {
  
  label    <- sprintf("%s-%02d", toupper(index_type), scale)
  min_dur  <- as.integer(min_duration)
  
  dates_req <- seq(as.Date(sprintf("%d-01-01", start_year_rec)),
                   as.Date(sprintf("%d-12-01", end_year_rec)), by = "month")
  bs <- tryCatch(read_basin_series(tolower(index_type), scale, get_event_index_dir(index_type),
                                   times = dates_req,
                                   variant = if (index_type %in% c("ssmi", "ssi_sm")) SOIL_LAYER else NULL),
                 error = function(e) { cat(sprintf("  ⚠ %s: %s\n", label, conditionMessage(e))); NULL })
  if (is.null(bs)) return(NULL)
  basin_avg <- as.numeric(bs$basin$value)
  dates_ts  <- bs$basin$time
  yr_ts     <- as.integer(format(dates_ts, "%Y"))
  
  # ── Pre-event structural validation ────────────────────────────────────────
  validate_monthly_record(
    dates_ts,
    start_date = as.Date(sprintf("%d-01-01", start_year_rec)),
    end_date   = as.Date(sprintf("%d-12-01", end_year_rec)),
    label = sprintf("%s basin-average calendar", label),
    require_complete = TRUE
  )
  validate_standardized_index(
    basin_avg,
    label = sprintf("%s basin-average standardized index", label)
  )
  
  # ── SW08 event identification ──────────────────────────────────────────────
  events <- sw08_event_loop(basin_avg, dates_ts, yr_ts,
                            onset_thr    = onset_thr,
                            min_duration = min_dur,
                            end_thr      = DROUGHT_END_EVENT)
  
  # If the record ends while a drought run is still active, retain a separate
  # right-censored representation for the recent-event bubble plot. It does not
  # enter the historical completed-event catalog or severity ranking because
  # the event has not terminated before the end of the record.
  open_event <- NULL
  n_ts <- length(basin_avg)
  
  # The still-active run at the record end, found with the same engine but
  # without the minimum-duration gate (it is shorter than min_dur here).
  runs_open <- drought_event_runs(basin_avg, min_duration = 1L,
                                  onset_threshold = onset_thr,
                                  termination_threshold = DROUGHT_END_EVENT)
  runs_open <- runs_open[runs_open$right_censored, , drop = FALSE]
  
  if (n_ts > 0L && nrow(runs_open) == 1L) {
    open_end <- runs_open$end_index[1]
    open_start <- runs_open$start_index[1]
    
    open_duration <- open_end - open_start + 1L
    
    if (open_duration < min_dur &&
        dates_ts[open_end] >= RECENT_DROUGHT_START &&
        dates_ts[open_start] <= RECENT_DROUGHT_END) {
      open_seg <- basin_avg[open_start:open_end]
      open_mean_int <- mean(DROUGHT_END_EVENT - open_seg)
      
      if (!is.finite(open_mean_int) || open_mean_int < 0) {
        stop("Internal SW08 error: invalid right-censored event intensity for ", label)
      }
      
      open_event <- data.frame(
        start_date = dates_ts[open_start],
        end_date = dates_ts[open_end],
        start_year = yr_ts[open_start],
        duration_months = as.integer(open_duration),
        duration_class = classify_duration(open_duration),
        mean_int = open_mean_int,
        severity = open_mean_int * open_duration,
        ongoing = TRUE,
        censored_at_record_end = TRUE,
        min_duration_required = as.integer(min_dur),
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (nrow(events) == 0L && is.null(open_event)) {
    cat(sprintf("  ⚠ %s: no events detected (threshold=%.2f, min_dur=%d m)\n",
                label, onset_thr, min_dur))
    return(NULL)
  }
  
  events$index_type  <- rep(index_type, nrow(events))
  events$scale       <- rep(scale, nrow(events))
  events$index_label <- rep(sprintf("%s_%02d", index_type, scale), nrow(events))
  attr(events, "open_event") <- open_event
  
  # save per-index CSV for reproducibility
  out_csv <- file.path(RANKING_DIR, sprintf("%s_SW08_event_catalog.csv", label))
  write.csv(events, out_csv, row.names = FALSE)
  
  cat(sprintf("  ✓ %s: %d events  [thr=%.2f, min_dur=%d m]  (saved %s)\n",
              label, nrow(events), onset_thr, min_dur, basename(out_csv)))
  events
}
#' Identify if event overlaps with configured recent window period
is_recent_drought <- function(start_date, end_date) {
  start_date <= RECENT_DROUGHT_END & end_date >= RECENT_DROUGHT_START
}

####################################################################################
# PART 1: BUILD EVENT CATALOGS FROM PRODUCER BASIN SERIES
#
# Onset      : DROUGHT_ONSET_EVENT (= DROUGHT_ONSET in utils)
# Termination: DROUGHT_END_EVENT (= DROUGHT_END in utils; no hysteresis if equal onset)
# Min duration: scale-dependent, EVENT_MIN_DURATION_BY_SCALE in utils
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 1: Building Event Catalogs from Producer Basin Series\n")
cat(sprintf("  Onset            : index < %.2f  (%.0fth pct of std normal)\n",
            DROUGHT_ONSET_EVENT, pnorm(DROUGHT_ONSET_EVENT) * 100))
cat(sprintf("  Termination      : index >= %.2f\n", DROUGHT_END_EVENT))
cat("  Intensity        : mean(termination_thr - x_t) per event month\n")
cat(sprintf("  Minimum duration : %s (scale=months)\n",
            paste(sprintf("%s=%d", names(MIN_DURATION_MAP), MIN_DURATION_MAP), collapse = ", ")))
cat("══════════════════════════════════════\n")

all_events <- list()

# Right-censored recent drought runs are retained separately for the recent-event
# bubble plot; they are not historical completed events because the record ends
# before the threshold crossing that would close the run.
open_recent_events <- list()

for (idx_info in INDICES_TO_ANALYZE) {
  for (sc in idx_info$scales) {
    catalog <- build_sw08_catalog(idx_info$index, sc,
                                  onset_thr      = DROUGHT_ONSET_EVENT,
                                  min_duration   = get_event_min_duration(sc),
                                  start_year_rec = HISTORICAL_START,
                                  end_year_rec   = HISTORICAL_END)
    if (!is.null(catalog)) {
      key <- sprintf("%s_%02d", idx_info$index, sc)
      
      open_candidate <- attr(catalog, "open_event")
      if (!is.null(open_candidate) && nrow(open_candidate)) {
        open_candidate$index_type <- idx_info$index
        open_candidate$scale <- as.integer(sc)
        open_candidate$index_label <- sprintf("%s_%02d", idx_info$index, sc)
        open_recent_events[[key]] <- open_candidate
      }
      
      if (nrow(catalog) > 0L) {
        all_events[[key]] <- catalog
      }
    }
  }
}

if (length(all_events) == 0) {
  stop("No SW08 catalogs could be built. Check producer basin series and paths in DROUGHT_ANALYSIS_utils.R.")
}

# Merge all catalogs
master_events <- dplyr::bind_rows(all_events)

validate_scales(
  unique(master_events$scale[master_events$index_type %in% c("spi", "spei")]),
  expected = EVENT_RANKING_SCALES,
  allow_extra = FALSE,
  label = "event-ranking scales"
)
if (any(!is.finite(master_events$severity)) || any(master_events$severity < 0))
  stop("Master event catalog contains invalid severity values.")
if (any(!is.finite(master_events$mean_int)) || any(master_events$mean_int < 0))
  stop("Master event catalog contains invalid intensity values.")

# Safety: ensure types
master_events$duration_months <- as.integer(master_events$duration_months)
master_events$start_date      <- as.Date(master_events$start_date)
master_events$end_date        <- as.Date(master_events$end_date)

cat(sprintf("\n✓ Total SW08 events built: %d  (period: %d-%d)\n",
            nrow(master_events), HISTORICAL_START, HISTORICAL_END))
################################################################################
# PART 2: CALCULATE RANKING METRICS
# mean_int and severity are already computed by build_sw08_catalog().
# This part adds duration_class, is_recent, rank, return_period, percentile.
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 2: Calculating Ranking Metrics\n")
cat("══════════════════════════════════════\n")

# Classify duration (SW08 classes)
master_events$duration_class <- classify_duration(master_events$duration_months)

# Mark recent drought (configured recent window)
master_events$is_recent <- mapply(is_recent_drought,
                                  master_events$start_date,
                                  master_events$end_date)

# Calculate year range for return period
n_years_record <- HISTORICAL_END - HISTORICAL_START + 1

# Rank by severity (within each index-scale combination)
master_events <- master_events %>%
  group_by(index_label) %>%
  arrange(desc(severity)) %>%
  mutate(
    rank_by_severity = dplyr::min_rank(dplyr::desc(severity)),
    return_period_years = calc_return_period(rank_by_severity, n_years_record),
    percentile_rank = round(100 * (1 - rank_by_severity / (n() + 1)), 1)
  ) %>%
  ungroup()

# Identify the most severe recent event
recent_events <- master_events %>% filter(is_recent)
if (nrow(recent_events) > 0) {
  most_severe_recent <- recent_events %>%
    dplyr::group_by(index_label) %>%
    dplyr::slice_max(severity, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  most_severe_overall <- recent_events %>%
    dplyr::slice_max(severity, n = 1, with_ties = FALSE)
  
  cat(sprintf(paste0("\n✓ Identified %d recent drought events (", RECENT_TAG, ")\n"),
              nrow(recent_events)))
  cat(sprintf("  Most severe overall: %s (severity=%.2f, rank=%d)\n",
              most_severe_overall$index_label[1],
              most_severe_overall$severity[1],
              most_severe_overall$rank_by_severity[1]))
} else {
  cat(paste0("\n⚠ No events detected in ", RECENT_TAG, " period\n"))
  most_severe_recent <- NULL
}

####################################################################################
# PART 3: SAVE RANKED CATALOG (SW08)
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 3: Saving SW08 Ranked Catalog\n")
cat("══════════════════════════════════════\n")

# Select columns for export
export_cols <- c("index_label", "start_date", "end_date", "start_year",
                 "duration_months", "duration_class", "mean_int",
                 "severity", "rank_by_severity", "return_period_years",
                 "percentile_rank", "is_recent")

ranked_catalog <- master_events[, export_cols]
write.csv(ranked_catalog,
          file.path(RANKING_DIR, "ranked_event_catalog_SW08.csv"),
          row.names = FALSE)

cat(sprintf("✓ Saved: ranked_event_catalog_SW08.csv (%d rows)\n", nrow(ranked_catalog)))

# Export a persistent methods note next to the two recurrence-interval tables.
# Only the class-count interval has a Poisson-based confidence interval; neither
# quantity represents a calibrated annual exceedance probability.
writeLines(c(
  "SW08 EVENT RECURRENCE INTERVALS — DEFINITIONS",
  sprintf("Historical record: %d–%d; N = %d calendar years.",
          HISTORICAL_START, HISTORICAL_END, n_years_record),
  "",
  "return_period_by_duration_class_SW08.csv:",
  "  return_period_yrs = N / n_events, where n_events counts events in a",
  "  specific duration class for that index/scale, irrespective of severity.",
  "  Values are rounded to 0.1 year in this CSV.",
  "  rp_lower_95 / rp_upper_95 invert exact Garwood Poisson count intervals.",
  "  Coverage assumes a homogeneous Poisson event count (not proven by the data).",
  "  Missing duration classes have no estimated interval, not zero recurrence.",
  "",
  "ranked_event_catalog_SW08.csv:",
  "  return_period_years = (N+1) / rank_by_severity for each index/scale,",
  "  ranking all event durations together by composite severity S = I × D.",
  "  The empirical rank-based plotting position has no confidence interval.",
  sprintf("  Rank 1 gives N+1 = %d years; this is the record-length plotting-position",
          n_years_record + 1L),
  "  cap, not a claim that the true recurrence interval cannot be longer.",
  "  Ties use min_rank(desc(severity)); a tied rank is NOT the count of",
  "  events with severity greater than or equal to the target severity.",
  "",
  "INTERPRETING THE FIGURES:",
  "  A circle is conditioned on duration class only, and ignores severity.",
  "  A diamond is conditioned on severity rank across ALL durations; its",
  "  horizontal duration-class location is for reference, not a class-specific rank.",
  "  The two plotted vertical quantities are therefore different estimands.",
  "  Neither is an annual exceedance probability. Multiple events per year can",
  "  produce recurrence intervals shorter than one year.",
  "  Kendall T_K is a third, model-based joint duration-and-severity recurrence",
  "  interval, not interchangeable with either event-count or severity-rank interval."
), file.path(RANKING_DIR, "README_recurrence_intervals_SW08.txt"), useBytes = TRUE)

####################################################################################
# PART 3b: DURATION-CLASS RETURN PERIOD TABLE — SW08  (Sheffield & Wood 2008)
#
# Events follow the framework-wide definition (utils, section F): onset when the
# index < DROUGHT_ONSET, termination when it is >= DROUGHT_END, intensity = mean
# deficit below DROUGHT_END, severity = intensity x duration, and the scale-
# dependent minimum duration EVENT_MIN_DURATION_BY_SCALE. D1-2 events can occur
# only where that minimum is <= 2 months (scale 1); they stay in the catalog but
# are excluded here, because the SW08 duration classes start at D3.
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 3b: Duration-Class Return Period Table (SW08)\n")
cat("══════════════════════════════════════\n")

rp_table <- master_events %>%
  dplyr::filter(duration_class != "D1-2 (Very short)") %>%
  dplyr::group_by(index_label, duration_class) %>%
  dplyr::summarise(
    n_events         = dplyr::n(),
    mean_duration_mo = round(mean(duration_months, na.rm = TRUE), 1),
    mean_I           = round(mean(mean_int,   na.rm = TRUE), 4),  # S&W intensity
    mean_S           = round(mean(severity,        na.rm = TRUE), 3),  # S&W severity
    max_S            = round(max(severity,         na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # Central mean recurrence interval
    return_period_yrs = round(n_years_record / n_events, 1),
    # Exact Poisson 95 % confidence interval (chi-squared / Garwood method):
    rp_lower_95 = round( (2 * n_years_record) / qchisq(0.975, 2 * (n_events + 1)), 1 ),
    rp_upper_95 = round( (2 * n_years_record) / qchisq(0.025, 2 * n_events), 1 ),
    duration_class = factor(duration_class,
                            levels = c("D3-6 (Short-term)",
                                       "D7-12 (Medium-term)",
                                       "D13+ (Long-term)"))
  ) %>%
  dplyr::arrange(index_label, duration_class) %>%
  dplyr::select(index_label, duration_class, n_events,
                return_period_yrs, rp_lower_95, rp_upper_95,
                mean_duration_mo, mean_I, mean_S, max_S)

write.csv(rp_table,
          file.path(RANKING_DIR, "return_period_by_duration_class_SW08.csv"),
          row.names = FALSE)
cat(sprintf("✓ Saved: return_period_by_duration_class_SW08.csv (%d rows)\n", nrow(rp_table)))

# Legend shapes denote different event-based recurrence estimands, not merely
# two symbols for the same conditional return period.
RP_SHAPES <- c(
  "Duration class: N / n events (95% Poisson CI)" = 16,
  "Recent event: (N+1) / severity rank, all durations" = 18
)
RP_CIRCLE <- names(RP_SHAPES)[1]
RP_DIAMOND <- names(RP_SHAPES)[2]
index_event_counts <- master_events %>%
  dplyr::count(index_label, name = "n_events_index")
rank_rp_label <- function(severity, return_period_years, rank_by_severity,
                          n_events_index) {
  sprintf("%s\nS \u2265 %.1f: ~%.0f yr\n(rank %d of %d%s)",
          RECENT_LABEL, severity, return_period_years,
          rank_by_severity, n_events_index,
          ifelse(rank_by_severity == 1L,
                 sprintf("; record-length cap %d yr", n_years_record + 1L), ""))
}

# --- Figure: return period dot-range chart, one panel per index ---------------
if (nrow(rp_table) > 0) {
  p_rp <- ggplot2::ggplot(
    rp_table %>% dplyr::filter(!is.na(return_period_yrs)),
    ggplot2::aes(x = duration_class, y = return_period_yrs,
                 colour = index_label, group = index_label)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = rp_lower_95, ymax = rp_upper_95),
      width = 0.25, linewidth = 0.8, position = ggplot2::position_dodge(0.6)
    ) +
    ggplot2::geom_point(ggplot2::aes(shape = RP_CIRCLE),
                        size = 3, position = ggplot2::position_dodge(0.6)) +
    ggplot2::scale_shape_manual(values = RP_SHAPES, name = NULL) +
    ggplot2::scale_y_continuous(
      name   = "Mean recurrence interval (years, event-based)",
      breaks = c(1, 2, 5, 10, 20, 50, 100),
      trans  = "log10",
      labels = scales::comma
    ) +
    ggplot2::scale_colour_viridis_d(option = "turbo", name = "Index", end = 0.95) +
    ggplot2::labs(
      title    = "Event recurrence intervals by duration class — Nechako Watershed",
      subtitle = sprintf("Circles: N/n_class with 95%% Poisson CI  |  %d–%d; N = %d years",
                         HISTORICAL_START, HISTORICAL_END, n_years_record),
      x        = "Duration class (Sheffield & Wood 2008)",
      caption  = "These event-count recurrence intervals are not annual exceedance probabilities."
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(file.path(RANKING_DIR, "FigRP_Return_Period_By_Duration_Class.png"),
                  p_rp, width = 11, height = 7, dpi = 300)
  cat("✓ Saved: FigRP_Return_Period_By_Duration_Class.png\n")
  
  # ── FigRP variant: same plot but with configured recent window drought overlaid ──────────────
  # For each index, the most severe event that overlaps configured recent window is marked as a
  # diamond (shape = 18).  Its return period comes from the empirical ranking already
  # computed in Part 2 (most_severe_recent$return_period_years).
  # The x-position is the recent event's observed duration class, but its y-value
  # is an all-duration severity-rank plotting position, NOT a within-class rank.
  # Comparing the vertical positions does not estimate within-class severity rarity.
  if (!is.null(most_severe_recent) && nrow(most_severe_recent) > 0) {
    
    # Build a data frame of the configured recent window event per index, keeping only the
    # duration classes that appear in the historical rp_table axes
    valid_dc <- c("D3-6 (Short-term)", "D7-12 (Medium-term)", "D13+ (Long-term)")
    
    recent_rp_overlay <- most_severe_recent %>%
      dplyr::left_join(index_event_counts, by = "index_label") %>%
      dplyr::filter(!is.na(return_period_years),
                    duration_class %in% valid_dc) %>%
      dplyr::mutate(
        duration_class = factor(duration_class, levels = valid_dc),
        rp_label = rank_rp_label(severity, return_period_years,
                                 rank_by_severity, n_events_index)
      )
    
    if (nrow(recent_rp_overlay) > 0) {
      
      p_rp_recent <- p_rp +
        # ── Large diamond for each index × duration-class combination ──────────
      # Diamonds use an all-duration severity rank, unlike the N/n class circles.
      ggplot2::geom_point(
        data     = recent_rp_overlay,
        ggplot2::aes(x     = duration_class,
                     y     = return_period_years,
                     colour = index_label,
                     group  = index_label,
                     shape  = RP_DIAMOND),
        size     = 5,
        position = ggplot2::position_dodge(0.6),
        show.legend = c(colour = FALSE, shape = TRUE)
      ) +
        # ── Text label above each diamond ──────────────────────────────────────
      ggplot2::geom_label(
        data     = recent_rp_overlay,
        ggplot2::aes(x      = duration_class,
                     y      = return_period_years,
                     colour = index_label,
                     group  = index_label,
                     label  = rp_label),
        size        = 2.8,
        vjust       = -0.45,
        fontface    = "bold",
        linewidth  = 0.25,
        label.padding = ggplot2::unit(0.15, "lines"),
        position    = ggplot2::position_dodge(0.6),
        show.legend = FALSE
      ) +
        ggplot2::labs(
          title    = "Event recurrence intervals by duration class — Nechako Watershed",
          subtitle = sprintf(paste0(
            "Circles: N/n_class (95%% Poisson CI); diamonds: (N+1)/severity rank, all durations\n",
            "Diamond x-position identifies class only; rank-1 cap = %d yr | %d–%d"),
            n_years_record + 1L, HISTORICAL_START, HISTORICAL_END
          ),
          caption  = "Neither quantity is an annual exceedance probability; see README_recurrence_intervals_SW08.txt."
        )
      
      ggplot2::ggsave(
        file.path(RANKING_DIR,
                  paste0("FigRP_Return_Period_By_Duration_Class_with_", RECENT_TAG, ".png")),
        p_rp_recent, width = 11, height = 7, dpi = 300
      )
      cat(paste0("✓ Saved: FigRP_Return_Period_By_Duration_Class_with_", RECENT_TAG, ".png\n"))
    } else {
      cat(paste0("  ⚠ ", RECENT_TAG, " event has no valid duration class for RP overlay – skipped\n"))
    }
  } else {
    cat(paste0("  ⚠ most_severe_recent is NULL – ", RECENT_TAG, " RP overlay not drawn\n"))
  }
}

####################################################################################
# PART 4: CREATE PRESENTATION FIGURES
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 4: Creating Presentation Figures\n")
cat("══════════════════════════════════════\n")

# The headline index (HEADLINE_INDEX, configuration) is used for every presentation product.
PRIMARY_INDEX <- HEADLINE_INDEX
if (!(PRIMARY_INDEX %in% unique(master_events$index_label))) {
  PRIMARY_INDEX <- unique(master_events$index_label)[1]
  cat(sprintf("⚠ %s not available, using %s instead\n", HEADLINE_INDEX, PRIMARY_INDEX))
}

primary_events <- master_events %>% filter(index_label == PRIMARY_INDEX)
# Presentation and dashboard metrics must refer to the SAME index as primary_events.
msr_primary <- if (is.null(most_severe_recent)) {
  master_events[0, , drop = FALSE]
} else {
  dplyr::filter(most_severe_recent, index_label == PRIMARY_INDEX)
}

####################################################################################
# FIGURE 1: Top 10-15 Droughts by Severity (with configured recent window highlighted)
####################################################################################
cat("\n  Fig 1: Top Droughts by Severity...\n")

top_n <- 15
top_events <- primary_events %>%
  arrange(desc(severity)) %>%
  head(top_n)

# Ensure recent drought is included even if not in top 15
if (any(top_events$is_recent)) {
  plot_events <- top_events
} else {
  recent_top <- primary_events %>%
    filter(is_recent) %>%
    arrange(desc(severity)) %>%
    head(1)
  if (nrow(recent_top) > 0) {
    plot_events <- dplyr::bind_rows(
      head(top_events, min(nrow(top_events), top_n - 1L)), recent_top
    ) %>%
      arrange(desc(severity))
  } else {
    plot_events <- top_events
  }
}

plot_events$label <- sprintf("%s to %s\n(%d mo; rank #%d)",
                             format(plot_events$start_date, "%Y-%m"),
                             format(plot_events$end_date, "%Y-%m"),
                             plot_events$duration_months,
                             plot_events$rank_by_severity)

p1 <- ggplot2::ggplot(plot_events,
                      ggplot2::aes(x = reorder(label, severity), y = severity)) +
  ggplot2::geom_col(data = plot_events %>% filter(!is_recent),
                    fill = "#3498db", alpha = 0.7) +
  ggplot2::geom_col(data = plot_events %>% filter(is_recent),
                    fill = "#e74c3c", alpha = 0.9) +
  ggplot2::geom_text(data = plot_events %>% filter(is_recent),
                     ggplot2::aes(label = sprintf("Rank #%d\n(%.1fth percentile)",
                                                  rank_by_severity, percentile_rank)),
                     vjust = -0.5, size = 3.5, fontface = "bold", color = "#c0392b") +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::labs(
    title = sprintf("Top %d Historical Droughts in Nechako Watershed by Severity", top_n),
    subtitle = sprintf(paste0("%s | %s | Red = events overlapping ", RECENT_TAG),
                       event_index_label(PRIMARY_INDEX), ANALYSIS_PERIOD_LABEL),
    x = "Drought Event",
    y = "Severity (Intensity × Duration)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(color = "gray40", size = 11),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.major.y = ggplot2::element_line(color = "gray90")
  )

ggplot2::ggsave(file.path(RANKING_DIR, "Fig1_Top_Droughts_By_Severity.png"),
                p1, width = 12, height = 8, dpi = 300)
cat("    ✓ Saved: Fig1_Top_Droughts_By_Severity.png\n")

####################################################################################
# FIGURE 2: Severity Distribution with Percentile Rank (CDF)
####################################################################################
cat("\n  Fig 2: Severity Distribution (CDF)...\n")

# Calculate empirical CDF
cdf_data <- primary_events %>%
  arrange(severity) %>%
  mutate(cdf = row_number() / (n() + 1))

# Annotate only if the primary index has a qualifying recent event.
recent_sev <- if (nrow(msr_primary)) msr_primary$severity[1] else NA_real_
recent_pct <- if (nrow(msr_primary)) msr_primary$percentile_rank[1] else NA_real_
has_recent_primary <- length(recent_sev) == 1L && is.finite(recent_sev) &&
  length(recent_pct) == 1L && is.finite(recent_pct)

p2 <- ggplot2::ggplot(cdf_data, ggplot2::aes(x = severity, y = cdf * 100)) +
  ggplot2::geom_line(color = "#2c3e50", linewidth = 1.2) +
  ggplot2::geom_point(color = "#2c3e50", size = 1.5, alpha = 0.5) +
  {if (has_recent_primary)
    ggplot2::geom_vline(xintercept = recent_sev, linetype = "dashed",
                        color = "#e74c3c", linewidth = 1.5)} +
  {if (has_recent_primary)
    ggplot2::annotate("segment",
                      x = recent_sev, xend = recent_sev,
                      y = 0, yend = recent_pct,
                      color = "#e74c3c", linewidth = 1, linetype = "dotted")} +
  {if (has_recent_primary)
    ggplot2::annotate("point",
                      x = recent_sev, y = recent_pct,
                      color = "#e74c3c", size = 4, shape = 21, fill = "white", stroke = 1.5)} +
  ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  ggplot2::labs(
    title = "Historical Drought Severity Distribution",
    subtitle = if (has_recent_primary)
      sprintf(paste0("%s | %s | most severe event overlapping ", RECENT_TAG, " at %.1fth percentile"),
              event_index_label(PRIMARY_INDEX), ANALYSIS_PERIOD_LABEL, recent_pct)
    else sprintf(paste0("%s | %s | No qualifying ", RECENT_TAG, " event"),
                 event_index_label(PRIMARY_INDEX), ANALYSIS_PERIOD_LABEL),
    x = "Drought Severity",
    y = "Cumulative Frequency (%)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    panel.grid.major = ggplot2::element_line(color = "gray90")
  )

ggplot2::ggsave(file.path(RANKING_DIR, "Fig2_Drought_Severity_Distribution.png"),
                p2, width = 10, height = 7, dpi = 300)
cat("    ✓ Saved: Fig2_Drought_Severity_Distribution.png\n")

####################################################################################
# FIGURE 3: Drought Frequency by Duration Class (Sheffield & Wood 2008 style)
####################################################################################
cat("\n  Fig 3: Drought Frequency by Duration Class...\n")

duration_summary <- primary_events %>%
  dplyr::group_by(duration_class) %>%
  dplyr::summarise(
    n_events      = dplyr::n(),
    mean_duration = mean(duration_months, na.rm = TRUE),
    mean_severity = mean(severity,         na.rm = TRUE),
    max_severity  = max(severity,          na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    return_period = n_years_record / pmax(n_events, 1), # avoid divide-by-zero
    duration_class = factor(
      duration_class,
      levels = c("D1-2 (Very short)",
                 "D3-6 (Short-term)",
                 "D7-12 (Medium-term)",
                 "D13+ (Long-term)")
    )
  )

p3 <- ggplot2::ggplot(duration_summary,
                      ggplot2::aes(x = duration_class, y = n_events, fill = duration_class)) +
  ggplot2::geom_col(alpha = 0.8, color = "black", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("n=%d\nClass N/n: %.1f yr", n_events, return_period)),
    vjust = -0.5, size = 4, fontface = "bold"
  ) +
  ggplot2::scale_fill_brewer(palette = "YlOrRd") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
  ggplot2::labs(
    title = "Drought Frequency by Duration Class (Sheffield & Wood 2008)",
    subtitle = sprintf("%s Basin-Averaged | %d-%d | Nechako Watershed",
                       event_index_label(PRIMARY_INDEX), HISTORICAL_START, HISTORICAL_END),
    x = "Duration Class",
    y = "Number of Events"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    legend.position = "none",
    panel.grid.major.y = ggplot2::element_line(color = "gray90")
  )

ggplot2::ggsave(file.path(RANKING_DIR, "Fig3_Drought_Frequency_By_Duration.png"),
                p3, width = 10, height = 7, dpi = 300)
cat("    ✓ Saved: Fig3_Drought_Frequency_By_Duration.png\n")

####################################################################################
# FIGURE 4: Combined Dashboard (all key metrics for presentation)
####################################################################################
cat("\n  Fig 4: Creating Summary Dashboard...\n")

# Summary statistics table – now built as a proper table grob
summary_stats <- data.frame(
  Metric = c(sprintf("Total Drought Events (%s)", ANALYSIS_PERIOD_LABEL),
             "Short-term Droughts (3-6 mo)",
             "Medium-term Droughts (7-12 mo)",
             "Long-term Droughts (13+ mo)",
             "Most Severe Historical Drought",
             paste0("Most severe event overlapping ", RECENT_TAG, ": rank"),
             "Same event: plotting-position percentile",
             "Same event: rank recurrence proxy",
             "Mean Drought Duration",
             "Max Drought Duration"),
  Value = c(
    nrow(primary_events),
    sum(primary_events$duration_class == "D3-6 (Short-term)"),
    sum(primary_events$duration_class == "D7-12 (Medium-term)"),
    sum(primary_events$duration_class == "D13+ (Long-term)"),
    sprintf("%s–%s (Severity: %.1f)",
            format(primary_events$start_date[which.max(primary_events$severity)], "%Y"),
            format(primary_events$end_date[which.max(primary_events$severity)], "%Y"),
            max(primary_events$severity)),
    sprintf("#%d of %d",
            msr_primary$rank_by_severity[1] %||% NA,
            nrow(primary_events)),
    sprintf("%.1fth percentile", msr_primary$percentile_rank[1] %||% NA),
    sprintf("%.1f yr (NOT annual exceedance)", msr_primary$return_period_years[1] %||% NA),
    sprintf("%.1f months", mean(primary_events$duration_months)),
    sprintf("%d months", max(primary_events$duration_months))
  ),
  stringsAsFactors = FALSE
)

# Convert to a grid table
table_grob <- gridExtra::tableGrob(
  summary_stats,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 10,
    padding = unit(c(4, 2), "mm")
  )
)

# Wrap it so patchwork can handle it
p4_table <- patchwork::wrap_elements(table_grob)

# Combine figures
dashboard <- (p1 | p2) / (p3 | p4_table) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Nechako Watershed Drought Analysis for Climate Adaptation Roundtable",
    subtitle = sprintf("Historical record %s | Sheffield & Wood (2008) event method", ANALYSIS_PERIOD_LABEL),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(color = "gray50", size = 12, hjust = 0.5)
    )
  )

ggplot2::ggsave(file.path(RANKING_DIR, "FigS_Event_Presentation_Dashboard.png"),
                dashboard, width = 16, height = 12, dpi = 300)
cat("    ✓ Saved: FigS_Event_Presentation_Dashboard.png\n")

####################################################################################
# PART 4b: MANUSCRIPT FIGURE 3 — Drought Severity Ranking & Return Period
# (publication quality, two-panel, headline index)
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 4b: Manuscript Figure 3\n")
cat("══════════════════════════════════════\n")

# ── Index for the manuscript figure: the same headline index as Part 4 ───────
# HEADLINE_INDEX (configuration) drives Figs 1-6, the dashboard, Excel, this
# figure, the Part 4c Kendall RP and the MS2 export.
MS_INDEX <- PRIMARY_INDEX
ms_events <- master_events %>% dplyr::filter(index_label == MS_INDEX)
cat(sprintf("  Using index: %s  (%d events)\n", MS_INDEX, nrow(ms_events)))

# ── Shared style constants ────────────────────────────────────────────────────
COL_HIST   <- "#4393c3"
COL_RECENT <- "#d6604d"
N_TOP      <- 15

# Display label: "spei_03" -> "SPEI-03"
MS_INDEX_LAB <- event_index_label(MS_INDEX)

theme_ms3 <- ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    axis.text          = ggplot2::element_text(size = 9,  colour = "grey20"),
    axis.title         = ggplot2::element_text(size = 10, colour = "grey10"),
    plot.title         = ggplot2::element_text(size = 10, face = "bold", hjust = 0),
    plot.subtitle      = ggplot2::element_text(size = 8.5, colour = "grey40", hjust = 0,
                                               margin = ggplot2::margin(b = 4)),
    legend.position    = "bottom",
    legend.key.size    = ggplot2::unit(0.45, "cm"),
    legend.text        = ggplot2::element_text(size = 8.5),
    legend.title       = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor   = ggplot2::element_blank(),
    plot.margin        = ggplot2::margin(6, 14, 4, 6)
  )

# ── PANEL (a): Ranked horizontal bar chart ────────────────────────────────────
# Composite score = severity (mean_int × duration_months).
# This mirrors the ranking criterion used in Table 2 (Section 4.1).
top_ms <- ms_events %>%
  dplyr::arrange(dplyr::desc(severity)) %>%
  dplyr::slice_head(n = N_TOP)

# Always include the most severe configured recent window event even if outside the top N
if (!any(top_ms$is_recent) && any(ms_events$is_recent)) {
  recent_top1 <- ms_events %>%
    dplyr::filter(is_recent) %>%
    dplyr::arrange(dplyr::desc(severity)) %>%
    dplyr::slice_head(n = 1)
  top_ms <- dplyr::bind_rows(
    head(top_ms, min(nrow(top_ms), N_TOP - 1L)), recent_top1
  ) %>%
    dplyr::arrange(dplyr::desc(severity))
}

# Build event labels: "YYYY" or "YYYY–YYYY" for multi-year events
top_ms <- top_ms %>%
  dplyr::mutate(
    yr_start  = as.integer(format(start_date, "%Y")),
    yr_end    = as.integer(format(end_date,   "%Y")),
    event_lab = dplyr::if_else(
      yr_start == yr_end,
      as.character(yr_start),
      sprintf("%d\u2013%d", yr_start, yr_end)
    ),
    # Start month uniquely identifies events even when severity ranks are tied.
    event_lab = sprintf("%s (m%s; #%d)", event_lab,
                        format(start_date, "%m"), rank_by_severity),
    # Rank label shown to the right of each bar for configured recent window
    rank_lab  = dplyr::if_else(
      is_recent,
      sprintf("Rank #%d  (%.1fth pct.)", rank_by_severity, percentile_rank),
      NA_character_
    ),
    fill_grp = dplyr::if_else(is_recent, RECENT_GROUP_LABEL, "Historical droughts")
  )

pa3 <- ggplot2::ggplot(top_ms,
                       ggplot2::aes(x     = reorder(event_lab, severity),
                                    y     = severity,
                                    fill  = fill_grp)) +
  ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  # Rank annotation for recent events (right of bar)
  ggplot2::geom_text(
    ggplot2::aes(label = rank_lab),
    hjust   = -0.08,
    size    = 2.9,
    colour  = COL_RECENT,
    fontface = "bold",
    na.rm   = TRUE
  ) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(
    values = c("Historical droughts" = COL_HIST,
               stats::setNames(COL_RECENT, RECENT_GROUP_LABEL)),
    guide  = ggplot2::guide_legend(reverse = TRUE)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.28)),
    labels = scales::comma
  ) +
  ggplot2::labs(
    title    = sprintf("(a)  Top %d drought events by composite severity score (%s, %d\u2013%d)",
                       N_TOP, MS_INDEX_LAB, HISTORICAL_START, HISTORICAL_END),
    subtitle = sprintf("S = mean(%.2f − index) × D  (onset < %.2f, end >= %.2f; scale-dependent minimum duration)",
                       DROUGHT_END_EVENT, DROUGHT_ONSET_EVENT, DROUGHT_END_EVENT),
    x        = NULL,
    y        = "Composite severity score"
  ) +
  theme_ms3 +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.major.y = ggplot2::element_blank()
  )

# ── PANEL (b): Return period by duration class with configured recent window overlay ─────────
# Re-uses rp_table computed in Part 3b, filtered to the manuscript index.
# Falls back to all indices if MS_INDEX rows are absent (e.g., label mismatch).
rp_ms <- rp_table %>% dplyr::filter(index_label == MS_INDEX)
if (nrow(rp_ms) == 0) {
  cat("  ⚠ No rp_table rows for MS_INDEX; using full rp_table for panel (b)\n")
  rp_ms <- rp_table
}

# Configured recent-window diamond; y is severity-ranked over ALL durations,
# while x shows the observed duration class for reference only.
valid_dc_lvls <- c("D3-6 (Short-term)", "D7-12 (Medium-term)", "D13+ (Long-term)")
recent_ms_rp  <- NULL
if (!is.null(most_severe_recent) && nrow(most_severe_recent) > 0) {
  recent_ms_rp <- most_severe_recent %>%
    dplyr::left_join(index_event_counts, by = "index_label") %>%
    dplyr::filter(index_label == MS_INDEX,
                  !is.na(return_period_years),
                  duration_class %in% valid_dc_lvls) %>%
    dplyr::mutate(
      duration_class = factor(duration_class, levels = valid_dc_lvls),
      rp_label = rank_rp_label(severity, return_period_years,
                               rank_by_severity, n_events_index)
    )
}
has_recent_ms_rp <- !is.null(recent_ms_rp) && nrow(recent_ms_rp) > 0L
RP_LEGEND_COLORS <- if (has_recent_ms_rp) c(COL_HIST, COL_RECENT) else COL_HIST
RP_LEGEND <- ggplot2::guide_legend(
  nrow = 2, byrow = TRUE, override.aes = list(colour = RP_LEGEND_COLORS)
)
RP_SUBTITLE <- sprintf(paste0(
  "Circles: N/n_class with 95%% Poisson CI | Diamond: (N+1)/severity rank across all durations\n",
  "Diamond x-position shows duration class only | %s | %d–%d | N = %d yr; rank-1 cap = %d yr"),
  MS_INDEX_LAB, HISTORICAL_START, HISTORICAL_END,
  n_years_record, n_years_record + 1L)
RP_CAPTION <- paste0(
  "Neither event-based quantity is an annual exceedance probability. ",
  "Kendall T_K, where shown, describes joint duration–severity recurrence."
)
# The no-CI variant must not promise uncertainty bars in its legend/subtitle.
RP_SHAPES_NOCI <- RP_SHAPES
names(RP_SHAPES_NOCI)[1] <- "Duration class: N / n events (CI omitted)"
RP_CIRCLE_NOCI <- names(RP_SHAPES_NOCI)[1]
RP_SUBTITLE_NOCI <- sub("with 95% Poisson CI", "(Poisson CI omitted)",
                        RP_SUBTITLE, fixed = TRUE)
RP_Y_BREAKS <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)
pb3 <- ggplot2::ggplot(
  rp_ms %>% dplyr::filter(!is.na(return_period_yrs)) %>%
    dplyr::mutate(duration_class = factor(duration_class, levels = valid_dc_lvls)),
  ggplot2::aes(x = duration_class, y = return_period_yrs)
) +
  # 95 % CI error bars
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = rp_lower_95, ymax = rp_upper_95),
    width     = 0.18,
    linewidth = 0.65,
    colour    = COL_HIST
  ) +
  # Historical mean dots
  ggplot2::geom_point(
    ggplot2::aes(shape = RP_CIRCLE),
    size   = 3.5,
    colour = COL_HIST
  ) +
  # configured recent window diamond overlay (added only if data available)
  {if (!is.null(recent_ms_rp) && nrow(recent_ms_rp) > 0)
    ggplot2::geom_point(
      data   = recent_ms_rp,
      ggplot2::aes(x = duration_class, y = return_period_years,
                   shape = RP_DIAMOND),
      size   = 5,
      colour = COL_RECENT,
      inherit.aes = FALSE
    )} +
  # Label above diamond
  {if (!is.null(recent_ms_rp) && nrow(recent_ms_rp) > 0)
    ggplot2::geom_label(
      data   = recent_ms_rp,
      ggplot2::aes(x = duration_class, y = return_period_years, label = rp_label),
      vjust       = -0.45,
      size        = 2.8,
      fontface    = "bold",
      colour      = COL_RECENT,
      linewidth  = 0.22,
      label.padding = ggplot2::unit(0.14, "lines"),
      inherit.aes = FALSE
    )} +
  ggplot2::scale_shape_manual(values = RP_SHAPES, name = NULL) +
  ggplot2::guides(shape = RP_LEGEND) +
  ggplot2::scale_y_continuous(
    name   = "Mean recurrence interval (years, event-based)",
    trans  = "log10",
    breaks = RP_Y_BREAKS,
    labels = scales::comma,
    expand = ggplot2::expansion(mult = c(0.08, 0.25))
  ) +
  ggplot2::scale_x_discrete(
    labels = c("D3\u20136\n(Short-term)",
               "D7\u201312\n(Medium-term)",
               "D13+\n(Long-term)")
  ) +
  ggplot2::labs(
    title    = "(b)  Event-based recurrence intervals",
    subtitle = RP_SUBTITLE,
    caption  = RP_CAPTION,
    x = "Duration class (diamond x-position is descriptive)"
  ) +
  theme_ms3 +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.major.x = ggplot2::element_blank()
  )

# ── Assemble two-panel figure — WITH Poisson CI ───────────────────────────────
fig3_ms <- pa3 / pb3 +
  patchwork::plot_layout(heights = c(1.35, 1)) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 4, 4))
  )

# ── Panel (b) WITHOUT CI error bars ──────────────────────────────────────────
pb3_noci <- ggplot2::ggplot(
  rp_ms %>% dplyr::filter(!is.na(return_period_yrs)) %>%
    dplyr::mutate(duration_class = factor(duration_class, levels = valid_dc_lvls)),
  ggplot2::aes(x = duration_class, y = return_period_yrs)
) +
  ggplot2::geom_point(ggplot2::aes(shape = RP_CIRCLE_NOCI),
                      size = 3.5, colour = COL_HIST) +
  {if (!is.null(recent_ms_rp) && nrow(recent_ms_rp) > 0)
    ggplot2::geom_point(
      data   = recent_ms_rp,
      ggplot2::aes(x = duration_class, y = return_period_years,
                   shape = RP_DIAMOND),
      size = 5, colour = COL_RECENT, inherit.aes = FALSE
    )} +
  {if (!is.null(recent_ms_rp) && nrow(recent_ms_rp) > 0)
    ggplot2::geom_label(
      data   = recent_ms_rp,
      ggplot2::aes(x = duration_class, y = return_period_years, label = rp_label),
      vjust = -0.45, size = 2.8, fontface = "bold", colour = COL_RECENT,
      linewidth = 0.22, label.padding = ggplot2::unit(0.14, "lines"),
      inherit.aes = FALSE
    )} +
  ggplot2::scale_shape_manual(values = RP_SHAPES_NOCI, name = NULL) +
  ggplot2::guides(shape = RP_LEGEND) +
  ggplot2::scale_y_continuous(
    name   = "Mean recurrence interval (years, event-based)", trans = "log10",
    breaks = RP_Y_BREAKS,
    labels = scales::comma,
    expand = ggplot2::expansion(mult = c(0.08, 0.25))
  ) +
  ggplot2::scale_x_discrete(
    labels = c("D3\u20136\n(Short-term)", "D7\u201312\n(Medium-term)", "D13+\n(Long-term)")
  ) +
  ggplot2::labs(
    title    = "(b)  Event-based recurrence intervals (class CI omitted)",
    subtitle = RP_SUBTITLE_NOCI,
    caption  = paste0(RP_CAPTION, " Class-count Poisson CI omitted in this variant."),
    x = "Duration class (diamond x-position is descriptive)"
  ) +
  theme_ms3 +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.major.x = ggplot2::element_blank()
  )

fig3_ms_noci <- pa3 / pb3_noci +
  patchwork::plot_layout(heights = c(1.35, 1)) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 4, 4))
  )

# ── Save PDF + PNG — both CI variants ─────────────────────────────────────────
save_fig <- function(fig, base_name, w = 7.0, h = 7.5) {
  for (ext in c("pdf", "png")) {
    path <- file.path(RANKING_DIR, paste0(base_name, ".", ext))
    tryCatch({
      if (ext == "pdf")
        ggplot2::ggsave(path, plot = fig, width = w, height = h, units = "in", device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf)
      else
        ggplot2::ggsave(path, plot = fig, width = w, height = h, units = "in",
                        dpi = 300, device = "png")
      cat(sprintf("  ✓ Saved: %s\n", basename(path)))
    }, error = function(e) cat(sprintf("  ⚠ Save failed (%s): %s\n", ext, e$message)))
  }
}

save_fig(fig3_ms,      "FigS_Event_Severity_Ranking")
save_fig(fig3_ms_noci, "FigS_Event_Severity_Ranking_noCI")

####################################################################################
####################################################################################
# PART 4c: CANONICAL KENDALL-BASED JOINT RETURN PERIOD ENGINE
#
# This is the SINGLE Kendall implementation used by this script.
#
# The same run_kendall_rp() function is called by:
#   • Part 4c: manuscript (headline-index) Kendall RP + diagnostics;
#   • Part 4d: multi-index comparison;
#   • Part 4g: all manuscript scales;
#   • Part 4h: regional series (with a supplied regional event table).
#
# There is deliberately no second all-scale Kendall implementation.
#
# Event definition:
#   • events come from master_events, generated by the SW08-style single-
#     threshold event loop above;
#   • D1-2 events remain in the master catalog, but are not used for the
#     duration-class RP table or Kendall duration-severity copula because the
#     SW08 duration classes used for RP analysis begin at D3.
#
# Marginals:
#   Duration: Gamma, Weibull, Exponential, Log-logistic
#   Severity: Gamma, Lognormal, Weibull, Log-logistic
#
# Copulas:
#   Clayton, Gumbel, Frank, Joe, AMH
#
# Log-logistic is fitted through fitdistrplus using registered d/p/q/r functions.
# Model-selection caution: with fewer than about 20 independent historical events,
# four candidate marginals per axis and five copulas can overfit. Duration is
# integer-valued and right-skewed; continuous marginal PITs are approximate.
# Treat point estimates and CIs as exploratory, not calibrated extremes.
####################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4c: Canonical Kendall RP engine\n")
cat("══════════════════════════════════════════════════════════════\n")

KENDALL_BOOT_DEFAULT <- 500L
# Console diagnostics for the computationally expensive bootstrap.
# Progress is reported without changing the statistical calculation.
KENDALL_PROGRESS_EVERY <- 10L

# Bootstrap and point estimation both use ML within the AIC-selected copula
# family. Degenerate resamples or optimizer failures remain NA (never replaced
# with clipped tau), and can prevent any CI being reported in small samples.
# Resampled copula is refitted with ML, as for the point estimate.
# Conditional on the originally AIC-selected family; family-selection uncertainty
# is not included. Event-rate variation uses an independent Poisson-count draw.
KENDALL_BOOT_COPULA_METHOD <- "ml"
KENDALL_BOOT_REPLICATE_TIMEOUT_SEC <- 5

# ----------------------------------------------------------------------------
# Custom log-logistic distribution for fitdistrplus.
# fitdist(x, "llogis_custom") resolves these functions as d/p/q/r + name.
# -----------------------------------------------------------------------------
dllogis_custom <- function(x, shape, scale, log = FALSE) {
  dens <- ifelse(
    x <= 0 | !is.finite(x) | shape <= 0 | scale <= 0,
    0,
    (shape / scale) * (x / scale)^(shape - 1) /
      (1 + (x / scale)^shape)^2
  )
  if (isTRUE(log)) {
    out <- rep(-Inf, length(x))
    ok <- is.finite(dens) & dens > 0
    out[ok] <- log(dens[ok])
    return(out)
  }
  dens
}

pllogis_custom <- function(q, shape, scale, lower.tail = TRUE, log.p = FALSE) {
  p <- ifelse(
    q <= 0 | !is.finite(q) | shape <= 0 | scale <= 0,
    0,
    1 / (1 + (scale / q)^shape)
  )
  if (!isTRUE(lower.tail)) p <- 1 - p
  if (isTRUE(log.p)) return(log(pmax(p, .Machine$double.xmin)))
  p
}

qllogis_custom <- function(p, shape, scale, lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(log.p)) p <- exp(p)
  if (!isTRUE(lower.tail)) p <- 1 - p
  p <- pmax(0, pmin(1, p))
  out <- rep(NA_real_, length(p))
  out[p <= 0] <- 0
  out[p >= 1] <- Inf
  ii <- p > 0 & p < 1 & shape > 0 & scale > 0
  out[ii] <- scale * (p[ii] / (1 - p[ii]))^(1 / shape)
  out
}

rllogis_custom <- function(n, shape, scale) {
  u <- stats::runif(n)
  qllogis_custom(u, shape, scale)
}

# ----------------------------------------------------------------------------
# Generic marginal fitting.
# -----------------------------------------------------------------------------
fit_kendall_marginal <- function(x, dist, start = NULL) {
  tryCatch(
    fitdistrplus::fitdist(
      x,
      distr = dist,
      start = start,
      method = "mle"
    ),
    error = function(e) NULL
  )
}

fit_kendall_marginals <- function(duration, severity) {
  duration_fits <- list(
    gamma = fit_kendall_marginal(duration, "gamma"),
    weibull = fit_kendall_marginal(duration, "weibull"),
    exp = fit_kendall_marginal(duration, "exp"),
    llogis = fit_kendall_marginal(
      duration,
      "llogis_custom",
      start = list(
        shape = 2,
        scale = max(stats::median(duration), .Machine$double.eps)
      )
    )
  )
  duration_fits <- Filter(Negate(is.null), duration_fits)
  
  severity_fits <- list(
    gamma = fit_kendall_marginal(severity, "gamma"),
    lnorm = fit_kendall_marginal(severity, "lnorm"),
    weibull = fit_kendall_marginal(severity, "weibull"),
    llogis = fit_kendall_marginal(
      severity,
      "llogis_custom",
      start = list(
        shape = 2,
        scale = max(stats::median(severity), .Machine$double.eps)
      )
    )
  )
  severity_fits <- Filter(Negate(is.null), severity_fits)
  
  if (!length(duration_fits))
    stop("No duration marginal distributions could be fitted.")
  if (!length(severity_fits))
    stop("No severity marginal distributions could be fitted.")
  
  aic_duration <- vapply(duration_fits, function(f) f$aic, numeric(1))
  aic_severity <- vapply(severity_fits, function(f) f$aic, numeric(1))
  
  best_duration <- names(which.min(aic_duration))[1]
  best_severity <- names(which.min(aic_severity))[1]
  
  list(
    duration = duration_fits[[best_duration]],
    duration_name = best_duration,
    severity = severity_fits[[best_severity]],
    severity_name = best_severity,
    all_duration = duration_fits,
    all_severity = severity_fits,
    aic_duration = aic_duration,
    aic_severity = aic_severity
  )
}

kendall_pit <- function(x, fit, dist) {
  p <- switch(
    dist,
    gamma = stats::pgamma(
      x,
      shape = fit$estimate["shape"],
      rate = fit$estimate["rate"]
    ),
    weibull = stats::pweibull(
      x,
      shape = fit$estimate["shape"],
      scale = fit$estimate["scale"]
    ),
    exp = stats::pexp(
      x,
      rate = fit$estimate["rate"]
    ),
    lnorm = stats::plnorm(
      x,
      meanlog = fit$estimate["meanlog"],
      sdlog = fit$estimate["sdlog"]
    ),
    llogis = pllogis_custom(
      x,
      shape = fit$estimate["shape"],
      scale = fit$estimate["scale"]
    ),
    stop("Unsupported Kendall marginal: ", dist)
  )
  
  pmax(1e-4, pmin(1 - 1e-4, as.numeric(p)))
}

refit_kendall_marginal <- function(x, dist, old_fit) {
  if (dist == "llogis") {
    # Strip names from estimate slices before passing starts to fitdist().
    # Otherwise names such as shape.shape/scale.scale can cause every
    # log-logistic bootstrap marginal refit to fail.
    st_shape <- as.numeric(old_fit$estimate["shape"])
    st_scale <- as.numeric(old_fit$estimate["scale"])
    if (!is.finite(st_shape) || st_shape <= 0) st_shape <- 2
    if (!is.finite(st_scale) || st_scale <= 0) {
      st_scale <- max(stats::median(x, na.rm = TRUE), 1e-6)
    }
    return(fit_kendall_marginal(
      x, "llogis_custom",
      start = list(shape = st_shape, scale = st_scale)))
  }
  fit_kendall_marginal(x, dist)
}

# ----------------------------------------------------------------------------
# Archimedean Kendall distribution.
# K_C(t) = t - phi(t) / phi'(t).
# -----------------------------------------------------------------------------
kendall_KC <- function(t, cop) {
  t <- pmax(1e-9, pmin(1 - 1e-9, as.numeric(t)))
  phi_t <- copula::iPsi(cop, t)
  dphi_t <- copula::diPsi(cop, t)
  kc <- t - phi_t / dphi_t
  pmax(0, pmin(1 - 1e-9, as.numeric(kc)))
}

fit_kendall_copulas <- function(uv, tau_emp) {
  tau_init <- max(0.01, min(0.95, as.numeric(tau_emp)))
  templates <- list(
    Clayton = copula::claytonCopula(dim = 2L), Gumbel = copula::gumbelCopula(dim = 2L),
    Frank   = copula::frankCopula(dim = 2L),   Joe    = copula::joeCopula(dim = 2L),
    AMH     = copula::amhCopula(dim = 2L))
  theta_range <- list(Clayton = c(1e-3, 50), Gumbel = c(1.001, 30), Frank = c(1e-3, 50),
                      Joe = c(1.001, 30), AMH = c(1e-3, 0.999))
  ll_of <- function(f) as.numeric(stats::logLik(f))
  
  safe_fit <- function(nm) {
    template <- templates[[nm]]
    ll_theta <- function(th) {
      v <- suppressWarnings(sum(copula::dCopula(uv, copula::setTheta(template, th), log = TRUE)))
      if (is.finite(v)) v else -1e10
    }
    tryCatch({
      theta0 <- copula::iTau(template, min(tau_init, if (nm == "AMH") 0.32 else 0.95))
      f <- tryCatch(copula::fitCopula(template, uv, method = "ml", start = theta0),
                    error = function(e) NULL)
      th_best <- stats::optimize(ll_theta, theta_range[[nm]], maximum = TRUE)$maximum
      f_ok <- !is.null(f) && is.finite(ll_of(f))
      if (!f_ok || ll_theta(th_best) > ll_of(f) + 1e-6) {
        f2 <- tryCatch(copula::fitCopula(template, uv, method = "ml", start = th_best),
                       error = function(e) NULL)
        if (!is.null(f2) && is.finite(ll_of(f2)) && (!f_ok || ll_of(f2) >= ll_of(f))) f <- f2
      }
      if (is.null(f) || !is.finite(ll_of(f)) || ll_of(f) < -1e-6) return(NULL)
      f
    }, error = function(e) NULL)
  }
  valid <- function(f) tryCatch({
    cop <- f@copula
    is.finite(as.numeric(copula::pCopula(c(.5, .5), cop))) &&
      is.finite(kendall_KC(0.5, cop)) && is.finite(as.numeric(copula::tau(cop)))
  }, error = function(e) FALSE)
  
  fits <- stats::setNames(lapply(names(templates), safe_fit), names(templates))
  fits <- Filter(valid, Filter(Negate(is.null), fits))   # drops negative-theta fits that crash kendall_KC
  if (!length(fits)) stop("All Kendall copula fits failed.")
  aic <- vapply(fits, stats::AIC, numeric(1))
  best_name <- names(which.min(aic))[1]
  list(fits = fits, aic = aic, best_name = best_name, best_fit = fits[[best_name]],
       best_copula = fits[[best_name]]@copula)
}
# ----------------------------------------------------------------------------
# Canonical one-index Kendall RP analysis.
# This function is the ONLY Kendall engine used in Parts 4c, 4d, 4g and 4h.
# -----------------------------------------------------------------------------
run_kendall_rp <- function(
    idx_label,
    n_boot = KENDALL_BOOT_DEFAULT,
    seed = 42L,
    events_data = NULL,
    record_years = n_years_record,
    duration_shift = 0) {
  
  tryCatch({
    event_source <- if (is.null(events_data)) master_events else events_data
    record_years <- as.numeric(record_years)
    if (!is.finite(record_years) || record_years <= 0)
      stop("record_years must be a positive finite number")
    
    ev <- event_source %>%
      dplyr::filter(
        .data$index_label == .env$idx_label,
        .data$duration_class != "D1-2 (Very short)",
        is.finite(duration_months),
        is.finite(severity),
        severity > 0
      ) %>%
      dplyr::arrange(start_date, end_date)
    
    n_ev <- nrow(ev)
    if (n_ev < 8L)
      stop("too few qualifying historical events: ", n_ev)
    # Target events must lie within the D3+ duration domain used to fit the
    # Kendall marginals. D1-2 events remain in the empirical event catalog but
    # are not extrapolated through marginals fitted only to D3+ events.
    recent <- event_source %>%
      dplyr::filter(
        .data$index_label == .env$idx_label,
        .data$is_recent,
        .data$duration_class != "D1-2 (Very short)",
        is.finite(duration_months),
        duration_months >= 3,
        is.finite(severity),
        severity > 0
      ) %>%
      dplyr::arrange(dplyr::desc(severity), start_date, end_date)
    
    if (!nrow(recent))
      stop(paste0("no ", RECENT_TAG, " event for this index"))
    if (!is.finite(duration_shift) || duration_shift < 0 ||
        duration_shift >= min(ev$duration_months))
      stop("duration_shift must be in [0, min event duration): ", duration_shift)
    
    duration <- as.numeric(ev$duration_months) - duration_shift
    severity <- as.numeric(ev$severity)
    tau_emp <- suppressWarnings(stats::cor(duration, severity, method = "kendall"))
    if (!is.finite(tau_emp)) tau_emp <- 0
    
    cat(sprintf(
      "  Events for Kendall fit: %d | mean inter-arrival: %.2f yr | tau=%.4f\n",
      n_ev,
      record_years / n_ev,
      tau_emp
    ))
    
    # Duration is discrete and left-truncated by the event-duration gate.
    # duration_shift = 0 (default) fits raw months and is the MS2 export contract.
    # A positive shift fits (D - shift); use it only for sensitivity runs that
    # call run_kendall_rp() directly, never through the cached get_kendall_rp().
    marg <- fit_kendall_marginals(duration, severity)
    cat("  Duration AIC:")
    for (nm in names(marg$aic_duration))
      cat(sprintf("  %s=%.1f", nm, marg$aic_duration[nm]))
    cat("\n")
    cat(sprintf(
      "  Best duration marginal: %s (AIC=%.1f)\n",
      marg$duration_name,
      marg$aic_duration[marg$duration_name]
    ))
    
    cat("  Severity AIC:")
    for (nm in names(marg$aic_severity))
      cat(sprintf("  %s=%.1f", nm, marg$aic_severity[nm]))
    cat("\n")
    cat(sprintf(
      "  Best severity marginal: %s (AIC=%.1f)\n",
      marg$severity_name,
      marg$aic_severity[marg$severity_name]
    ))
    
    u_obs <- kendall_pit(
      duration,
      marg$duration,
      marg$duration_name
    )
    v_obs <- kendall_pit(
      severity,
      marg$severity,
      marg$severity_name
    )
    uv <- cbind(u_obs, v_obs)
    

    
    cop_fit <- fit_kendall_copulas(uv, tau_emp)
    best_fit <- cop_fit$best_fit
    best_cop <- cop_fit$best_copula
    best_cop_name <- cop_fit$best_name
    theta_hat <- as.numeric(best_cop@parameters[1])
    tau_fitted <- suppressWarnings(as.numeric(copula::tau(best_cop)))
    
    cat("\n  Copula AIC ranking:\n")
    for (nm in names(sort(cop_fit$aic))) {
      cop_nm <- cop_fit$fits[[nm]]@copula
      marker <- if (identical(nm, best_cop_name)) " <-- selected" else ""
      cat(sprintf(
        "    %-8s theta=%.4f  tau=%.4f  AIC=%7.2f%s\n",
        nm,
        as.numeric(cop_nm@parameters[1]),
        as.numeric(copula::tau(cop_nm)),
        cop_fit$aic[nm],
        marker
      ))
    }
    
    mu_ia <- record_years / n_ev
    
    calc_target <- function(target_duration, target_severity, cop = best_cop) {
      u0 <- kendall_pit(target_duration - duration_shift, marg$duration, marg$duration_name)
      v0 <- kendall_pit(target_severity, marg$severity, marg$severity_name)
      t_star <- as.numeric(copula::pCopula(c(u0, v0), cop))
      K_C <- kendall_KC(t_star, cop)
      T_K <- mu_ia / max(1e-12, 1 - K_C)
      c(
        u0 = u0,
        v0 = v0,
        t_star = t_star,
        K_C = K_C,
        T_K = T_K
      )
    }
    
    # Sensitivity across all fitted copula families using the same selected
    # marginal pair.
    sens_df <- do.call(rbind, lapply(names(cop_fit$fits), function(nm) {
      cop_nm <- cop_fit$fits[[nm]]@copula
      t_s <- as.numeric(copula::pCopula(
        c(
          kendall_pit(recent$duration_months[1]- duration_shift, marg$duration, marg$duration_name),
          kendall_pit(recent$severity[1], marg$severity, marg$severity_name)
        ),
        cop_nm
      ))
      kc_s <- kendall_KC(t_s, cop_nm)
      data.frame(
        Copula = nm,
        theta = round(as.numeric(cop_nm@parameters[1]), 4),
        tau = round(as.numeric(copula::tau(cop_nm)), 4),
        t_star = round(t_s, 4),
        K_C = round(kc_s, 5),
        T_K_yr = round(mu_ia / max(1e-12, 1 - kc_s), 0),
        stringsAsFactors = FALSE
      )
    }))
    
    set.seed(seed)
    
    # One paired-resample fit serves EVERY recent target for an index. This
    # changes bootstrap work from O(B * n_targets * fit_cost) to
    # O(B * (fit_cost + n_targets)); storage O(B * n_targets).
    bootstrap_targets <- function(target_duration, target_severity) {
      nt <- length(target_duration)
      if (nt != length(target_severity) || !nt) stop("Invalid Kendall bootstrap targets")
      boot_t0 <- proc.time()[["elapsed"]]
      cat(sprintf("    Shared bootstrap: %d target(s), %d draws, copula=%s\n",
                  nt, n_boot, toupper(KENDALL_BOOT_COPULA_METHOD)))
      values <- matrix(NA_real_, nrow = n_boot, ncol = nt)
      tau_b_store <- rep(NA_real_, n_boot)
      theta_b_store <- rep(NA_real_, n_boot)
      status_store <- rep("failed", n_boot)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
      for (b in seq_len(n_boot)) {
        setTimeLimit(cpu = KENDALL_BOOT_REPLICATE_TIMEOUT_SEC,
                     elapsed = KENDALL_BOOT_REPLICATE_TIMEOUT_SEC, transient = TRUE)
        values[b, ] <- tryCatch({
          ii <- sample.int(n_ev, n_ev, replace = TRUE)
          db <- duration[ii]; sb <- severity[ii]
          fd_b <- refit_kendall_marginal(db, marg$duration_name, marg$duration)
          fs_b <- refit_kendall_marginal(sb, marg$severity_name, marg$severity)
          if (is.null(fd_b) || is.null(fs_b)) stop("marginal refit failed")
          uv_b <- cbind(kendall_pit(db, fd_b, marg$duration_name),
                        kendall_pit(sb, fs_b, marg$severity_name))
          tau_b_store[b] <- suppressWarnings(stats::cor(db, sb, method = "kendall"))
          fit_b <- copula::fitCopula(best_cop, uv_b, method = "ml", start = theta_hat)
          cop_b <- fit_b@copula
          theta_b_store[b] <- as.numeric(cop_b@parameters[1L])
          # Account for uncertainty in events/year under a Poisson arrival model.
          # This assumption is only approximate for dependent drought events.
          n_rate_b <- stats::rpois(1L, n_ev)
          if (n_rate_b == 0L) stop("zero bootstrap event count")
          uv_target <- cbind(
            kendall_pit(target_duration- duration_shift, fd_b, marg$duration_name),
            kendall_pit(target_severity, fs_b, marg$severity_name))
          t_b <- as.numeric(copula::pCopula(uv_target, cop_b))
          kc_b <- kendall_KC(t_b, cop_b)
          status_store[b] <- "converged"
          (record_years / n_rate_b) / pmax(1e-12, 1 - kc_b)
        }, error = function(e) {
          status_store[b] <<- if (grepl("time limit", conditionMessage(e), ignore.case = TRUE))
            "timeout" else "failed"
          if (b <= 3L) cat(sprintf("      bootstrap draw %d: %s\n", b, conditionMessage(e)))
          rep(NA_real_, nt)
        })
        setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
        if (b == 1L || b == n_boot || b %% KENDALL_PROGRESS_EVERY == 0L) {
          cat(sprintf("      shared bootstrap %d/%d | elapsed %.1fs | valid %d\n",
                      b, n_boot, proc.time()[["elapsed"]] - boot_t0,
                      sum(status_store[seq_len(b)] == "converged")))
          flush.console()
        }
      }
      # Do not silently trim high-RP or high-tau tails. When too many ML
      # resamples fail, report NA CI rather than a selected-success CI.
      min_valid <- max(50L, ceiling(0.8 * n_boot))
      smry <- t(vapply(seq_len(nt), function(j) {
        v <- values[, j]; v <- v[is.finite(v)]
        if (length(v) < min_valid)
          return(c(ci_lo = NA_real_, ci_hi = NA_real_, ci_med = NA_real_, n_conv = length(v)))
        c(ci_lo = unname(stats::quantile(v, 0.025)),
          ci_hi = unname(stats::quantile(v, 0.975)),
          ci_med = unname(stats::median(v)), n_conv = length(v))
      }, numeric(4L)))
      list(summary = smry, values = values, tau = tau_b_store,
           theta = theta_b_store, status = status_store)
    }
    
    boot_all <- bootstrap_targets(recent$duration_months, recent$severity)
    primary_boot <- boot_all$values[, 1L]
    primary_boot <- primary_boot[is.finite(primary_boot)]
    
    rows <- lapply(seq_len(nrow(recent)), function(j) {
      z <- calc_target(
        recent$duration_months[j],
        recent$severity[j]
      )
      cat(sprintf(
        "    Recent target %d/%d: %s to %s | duration=%d mo | severity=%.4f\n",
        j, nrow(recent), recent$start_date[j], recent$end_date[j],
        recent$duration_months[j], recent$severity[j]
      ))
      flush.console()
      b_summary <- boot_all$summary[j, ]
      data.frame(
        index_label = idx_label,
        scale_months = as.integer(sub(".*_([0-9]+)$", "\\1", idx_label)),
        target_rank_recent = j,
        start_date = as.Date(recent$start_date[j]),
        end_date = as.Date(recent$end_date[j]),
        duration_months = as.integer(recent$duration_months[j]),
        severity = as.numeric(recent$severity[j]),
        n_historical_events = n_ev,
        mean_interarrival_years = mu_ia,
        tau_empirical = tau_emp,
        tau_fitted = tau_fitted,
        duration_marginal = marg$duration_name,
        duration_shift_months = as.numeric(duration_shift),  # 0 = raw-month marginal (MS2 contract)
        severity_marginal = marg$severity_name,
        best_copula = best_cop_name,
        copula_theta = theta_hat,
        u_duration = unname(z["u0"]),
        v_severity = unname(z["v0"]),
        t_star = unname(z["t_star"]),
        K_C = unname(z["K_C"]),
        Kendall_return_period_years = unname(z["T_K"]),
        CI_95_lo_years = unname(b_summary["ci_lo"]),
        CI_95_hi_years = unname(b_summary["ci_hi"]),
        bootstrap_median_years = unname(b_summary["ci_med"]),
        bootstrap_converged = as.integer(b_summary["n_conv"]),
        bootstrap_copula_method = KENDALL_BOOT_COPULA_METHOD,
        bootstrap_tau_max = NA_real_,  # no upper-tau truncation
        stringsAsFactors = FALSE
      )
    })
    
    rows_df <- dplyr::bind_rows(rows)
    
    cat(sprintf(
      "  Primary recent event: %s to %s | D=%d mo | S=%.3f | T_K=%.0f yr\n",
      rows_df$start_date[1],
      rows_df$end_date[1],
      rows_df$duration_months[1],
      rows_df$severity[1],
      rows_df$Kendall_return_period_years[1]
    ))
    
    list(
      index = idx_label,
      rows = rows_df,
      primary = rows[[1L]],
      n_ev = n_ev,
      tau_emp = round(tau_emp, 4),
      best_cop = best_cop_name,
      theta = round(theta_hat, 4),
      tau_fit = round(tau_fitted, 4),
      best_dur = marg$duration_name,
      best_sev = marg$severity_name,
      u0 = rows_df$u_duration[1],
      v0 = rows_df$v_severity[1],
      t_star = rows_df$t_star[1],
      Kc = rows_df$K_C[1],
      TK = rows_df$Kendall_return_period_years[1],
      ci_lo = rows_df$CI_95_lo_years[1],
      ci_hi = rows_df$CI_95_hi_years[1],
      ci_med = rows_df$bootstrap_median_years[1],
      n_boot_conv = rows_df$bootstrap_converged[1],
      tgt_dur = rows_df$duration_months[1],
      tgt_sev = rows_df$severity[1],
      mu_ia = mu_ia,
      aic_dur = marg$aic_duration,
      aic_sev = marg$aic_severity,
      aic_cop = cop_fit$aic,
      fits_dur = marg$all_duration,
      fits_sev = marg$all_severity,
      fit_best_dur = marg$duration,
      fit_best_sev = marg$severity,
      ev_cop = ev,
      uv_mat = uv,
      best_fit = best_fit,
      cop_fits = cop_fit$fits,
      sens_df = sens_df,
      T_boot = primary_boot,
      T_kendall = rows_df$Kendall_return_period_years[1],
      Kc0 = rows_df$K_C[1],
      t0 = rows_df$t_star[1]
    )
  }, error = function(e) {
    cat(sprintf(
      "  WARNING Kendall RP failed for %s: %s\n",
      idx_label,
      conditionMessage(e)
    ))
    list(
      index = idx_label,
      error = conditionMessage(e),
      TK = NA_real_,
      ci_lo = NA_real_,
      ci_hi = NA_real_
    )
  })
}

# Cache results so an index analyzed in Part 4c is not refit in Part 4d/4g.
KENDALL_RESULT_CACHE <- new.env(parent = emptyenv())
get_kendall_rp <- function(idx_label, n_boot = KENDALL_BOOT_DEFAULT) {
  if (exists(idx_label, envir = KENDALL_RESULT_CACHE, inherits = FALSE)) {
    return(get(idx_label, envir = KENDALL_RESULT_CACHE, inherits = FALSE))
  }
  ans <- run_kendall_rp(idx_label, n_boot = n_boot)
  assign(idx_label, ans, envir = KENDALL_RESULT_CACHE)
  ans
}

# Part 4c manuscript index: MS_INDEX (= HEADLINE_INDEX).
kendall_rp_result <- get_kendall_rp(MS_INDEX, n_boot = KENDALL_BOOT_DEFAULT)
if (!is.null(kendall_rp_result$error)) {
  cat(sprintf("  Part 4c Kendall analysis unavailable: %s\n", kendall_rp_result$error))
  kendall_rp_result <- NULL
} else {
  cat(sprintf(
    "  Canonical Kendall result retained for %s: T_K=%.0f yr [%.0f–%.0f yr]\n",
    MS_INDEX,
    kendall_rp_result$TK,
    kendall_rp_result$ci_lo,
    kendall_rp_result$ci_hi
  ))
}

# ── 9a. Update Fig 3 panel (b) with Kendall RP annotation ────────────────────
if (!is.null(kendall_rp_result)) {
  krp <- kendall_rp_result
  
  # Label WITH CI
  k_lab_ci   <- sprintf("T_K = %.0f yr\n[%.0f\u2013%.0f yr]\n%s",
                        krp$T_kendall, krp$ci_lo, krp$ci_hi, krp$best_cop)
  # Label WITHOUT CI
  k_lab_noci <- sprintf("T_K = %.0f yr\n%s",
                        krp$T_kendall, krp$best_cop)
  
  annotate_label <- function(lab) {
    ggplot2::annotate(
      "label",
      x             = nrow(rp_ms %>% dplyr::filter(!is.na(return_period_yrs))),
      y             = krp$T_kendall,
      label         = lab,
      hjust         = 1.08, vjust = 0.5,
      size          = 2.5, fontface = "bold",
      colour        = COL_RECENT, fill = "white",
      label.padding = ggplot2::unit(0.15, "lines"),
      label.r       = ggplot2::unit(0.1,  "lines"),
      linewidth    = 0.2
    )
  }
  
  # ── WITH CI ───────────────────────────────────────────────────────────────
  pb3_k <- pb3 +
    annotate_label(k_lab_ci) +
    ggplot2::labs(subtitle = sprintf(
      "Class N/n (Poisson CI), severity (N+1)/rank, and Kendall T_K (joint D and S; bootstrap CI) | %s | %d–%d",
      MS_INDEX_LAB, HISTORICAL_START, HISTORICAL_END))
  
  fig3_k <- pa3 / pb3_k +
    patchwork::plot_layout(heights = c(1.35, 1)) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 4, 4)))
  
  fig3k_pdf <- file.path(RANKING_DIR, "FigS_Event_Severity_Ranking_KendallRP.pdf")
  fig3k_png <- file.path(RANKING_DIR, "FigS_Event_Severity_Ranking_KendallRP.png")
  tryCatch({
    ggplot2::ggsave(fig3k_pdf, fig3_k, width=7.0, height=7.5, units="in", device=if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf)
    cat(sprintf("  Saved (with CI) : %s\n", basename(fig3k_pdf)))
  }, error = function(e) cat(sprintf("  WARNING: %s\n", e$message)))
  tryCatch({
    ggplot2::ggsave(fig3k_png, fig3_k, width=7.0, height=7.5, units="in",
                    dpi=300, device="png")
    cat(sprintf("  Saved (with CI) : %s\n", basename(fig3k_png)))
  }, error = function(e) cat(sprintf("  WARNING: %s\n", e$message)))
  
  # ── WITHOUT CI ────────────────────────────────────────────────────────────
  pb3_k_noci <- pb3 +
    annotate_label(k_lab_noci) +
    ggplot2::labs(subtitle = sprintf(
      "Class N/n (Poisson CI), severity (N+1)/rank, and Kendall T_K (joint D and S; point only) | %s | %d–%d",
      MS_INDEX_LAB, HISTORICAL_START, HISTORICAL_END))
  
  fig3_k_noci <- pa3 / pb3_k_noci +
    patchwork::plot_layout(heights = c(1.35, 1)) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 4, 4)))
  
  fig3k_noci_pdf <- file.path(RANKING_DIR, "FigS_Event_Severity_Ranking_KendallRP_noCI.pdf")
  fig3k_noci_png <- file.path(RANKING_DIR, "FigS_Event_Severity_Ranking_KendallRP_noCI.png")
  tryCatch({
    ggplot2::ggsave(fig3k_noci_pdf, fig3_k_noci, width=7.0, height=7.5, units="in", device=if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf)
    cat(sprintf("  Saved (no CI)   : %s\n", basename(fig3k_noci_pdf)))
  }, error = function(e) cat(sprintf("  WARNING: %s\n", e$message)))
  tryCatch({
    ggplot2::ggsave(fig3k_noci_png, fig3_k_noci, width=7.0, height=7.5, units="in",
                    dpi=300, device="png")
    cat(sprintf("  Saved (no CI)   : %s\n", basename(fig3k_noci_png)))
  }, error = function(e) cat(sprintf("  WARNING: %s\n", e$message)))
}

# ── 9b. Six-panel diagnostic PDF ─────────────────────────────────────────────
if (!is.null(kendall_rp_result)) {
  krp        <- kendall_rp_result
  diag4c_pdf <- file.path(RANKING_DIR, "FigS_Event_Copula_Diagnostics.pdf")
  tryCatch({
    if (isTRUE(capabilities("cairo")))
      grDevices::cairo_pdf(diag4c_pdf, width = 12, height = 9)
    else
      grDevices::pdf(diag4c_pdf, width = 12, height = 9, useDingbats = FALSE)
    graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3.5, 1.5))
    
    # Panel 1: marginal AIC comparison (duration + severity side by side)
    aic_all <- c(setNames(krp$aic_dur, paste0("Dur:", names(krp$aic_dur))),
                 setNames(krp$aic_sev, paste0("Sev:", names(krp$aic_sev))))
    bar_col1 <- ifelse(
      names(aic_all) %in% c(paste0("Dur:", krp$best_dur),
                            paste0("Sev:", krp$best_sev)),
      "#2166ac", "grey70")
    graphics::barplot(aic_all, col = bar_col1, border = "white",
                      las = 2, cex.names = 0.8, ylab = "AIC",
                      main = "Marginal distribution AIC\n(blue = selected)")
    graphics::grid(nx = NA, ny = NULL)
    
    # Panel 2: copula AIC comparison
    bar_col2 <- ifelse(names(krp$aic_cop) == krp$best_cop, "#2166ac", "grey70")
    graphics::barplot(sort(krp$aic_cop), col = bar_col2[order(krp$aic_cop)],
                      border = "white", ylab = "AIC",
                      main = "Copula AIC comparison\n(blue = selected)",
                      names.arg = names(sort(krp$aic_cop)))
    graphics::abline(h = min(krp$aic_cop), lty = 2, col = "#d6604d", lwd = 1.5)
    graphics::grid(nx = NA, ny = NULL)
    
    # Panel 3: pseudo-observations scatter
    cop_vals <- copula::pCopula(krp$uv_mat, krp$best_fit@copula)
    plot(krp$uv_mat[, 1], krp$uv_mat[, 2],
         col  = grDevices::heat.colors(100, rev = TRUE)[
           cut(cop_vals, breaks = 100, labels = FALSE, include.lowest = TRUE)],
         pch  = 19, cex = 0.9,
         xlab = "u [duration CDF]",
         ylab = "v [severity CDF]",
         main = sprintf("Pseudo-observations\n%s copula (theta=%.3f, tau=%.3f)",
                        krp$best_cop, krp$theta, krp$tau_fit))
    graphics::points(krp$u0, krp$v0, pch = 18, cex = 2.4, col = "#d6604d")
    graphics::legend("bottomright",
                     legend = c("Historical", RECENT_TAG),
                     pch = c(19, 18), col = c("grey50", "#d6604d"),
                     bty = "n", cex = 0.85)
    graphics::grid()
    
    # Panel 4: K_C curves for ALL fitted copulas
    t_seq <- seq(0.001, 0.999, length.out = 300)
    col_map <- c(Clayton = "#4393c3", Gumbel = "#d6604d", Frank = "#74add1",
                 Joe = "#2166ac", AMH = "#878787")
    kc_mat <- sapply(names(krp$cop_fits), function(nm) {
      cop <- krp$cop_fits[[nm]]@copula
      sapply(t_seq, function(t) pmin(kendall_KC(t, cop), 0.9999))
    })
    y_rng <- range(kc_mat, na.rm = TRUE)
    plot(t_seq, kc_mat[, 1], type = "n",
         xlab = "t", ylab = "K_C(t)",
         main = "Kendall distribution K_C(t)\nall fitted copulas",
         xlim = c(0, 1), ylim = y_rng)
    graphics::lines(t_seq, t_seq, lty = 2, col = "grey80")  # independence
    for (nm in names(krp$cop_fits)) {
      col_i <- col_map[nm]
      if (is.na(col_i)) col_i <- "grey50"
      lwd_i <- if (nm == krp$best_cop) 2.5 else 1.2
      graphics::lines(t_seq, kc_mat[, nm], col = col_i, lwd = lwd_i)
    }
    # Mark t* on selected copula
    kc0_check <- pmin(kendall_KC(krp$t0, krp$best_fit@copula), 0.9999)
    graphics::abline(v = krp$t0, col = "#d6604d", lty = 2, lwd = 1)
    graphics::points(krp$t0, kc0_check, pch = 18, cex = 2, col = "#d6604d")
    leg_nms <- names(krp$cop_fits)
    graphics::legend("topleft",
                     legend = c(leg_nms, "Independence"),
                     col    = c(col_map[leg_nms], "grey80"),
                     lwd    = c(ifelse(leg_nms == krp$best_cop, 2.5, 1.2), 1),
                     lty    = c(rep(1, length(leg_nms)), 2),
                     bty    = "n", cex = 0.75)
    graphics::grid()
    
    # Panel 5: sensitivity barplot — T_K across all copulas
    sens <- krp$sens_df
    bar_col5 <- ifelse(sens$Copula == krp$best_cop, "#2166ac", "grey70")
    yk <- sens$T_K_yr
    graphics::barplot(yk, names.arg = sens$Copula, col = bar_col5,
                      border = "white",
                      ylab = "T_K (years)",
                      main = "Kendall RP sensitivity\nacross copula families",
                      log  = if (max(yk, na.rm=TRUE)/min(yk, na.rm=TRUE) > 10) "y" else "")
    graphics::grid(nx = NA, ny = NULL)
    
    # Panel 6: bootstrap histogram
    if (length(krp$T_boot) > 1L) {
      graphics::hist(krp$T_boot, breaks = 40, freq = FALSE,
                     col = "#b3cde3", border = "white",
                     xlab = "T_K (years)",
                     main = sprintf("Bootstrap T_K distribution\n%s copula + %s/%s marginals",
                                    krp$best_cop, krp$best_dur, krp$best_sev))
      graphics::abline(v = krp$T_kendall, col = "#2166ac", lwd = 2.5)
      ci_v <- c(krp$ci_lo, krp$ci_hi)
      if (any(is.finite(ci_v)))
        graphics::abline(v = ci_v[is.finite(ci_v)], col = "#d6604d", lwd = 1.8, lty = 2)
      graphics::legend("topright",
                       legend = c(sprintf("Point: %.0f yr", krp$T_kendall),
                                  sprintf("95%% CI: [%.0f, %.0f] yr", krp$ci_lo, krp$ci_hi)),
                       col = c("#2166ac", "#d6604d"), lty = c(1, 2), lwd = 2,
                       bty = "n", cex = 0.85)
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5,
                     sprintf("Bootstrap did not converge (%d/%d replicates)",
                             krp$n_boot_conv, KENDALL_BOOT_DEFAULT), cex = 1.1)
    }
    
    graphics::grid(nx = NA, ny = NULL)
    grDevices::dev.off()
    cat(sprintf("  Saved: %s\n", basename(diag4c_pdf)))
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
    cat(sprintf("  WARNING diagnostic PDF: %s\n", e$message))
  })
}


####################################################################################
# PART 4d: KENDALL RP MULTI-INDEX COMPARISON
#
# Runs the same bivariate copula / Kendall RP machinery from Part 4c across
# selected SPI and SPEI-PM scales plus SWEI-03, SSPI-01 and the
# soil-moisture indices at their native scales.
#
# Purpose: assess whether the duration-severity dependence (Kendall tau) and
#   the return period found for the headline index (MS_INDEX) are
#   index-specific or a general feature of this drought catalog. Longer
#   accumulation scales give fewer, longer events (typically higher tau);
#   short scales give many short events with weaker duration-severity coupling.
#   Soil-moisture indices are included at their native scales.
#
# Each index receives:
#   • Marginal distribution selection (Gamma / Weibull / Exp / Log-logistic)
#   • Copula selection from all 5 Archimedean families (AIC)
#   • Kendall RP point estimate + 500-rep bootstrap 95 % CI
#
# Outputs (RANKING_DIR):
#   Fig4d_KendallRP_MultiIndex_Comparison.pdf / .png
#   Table_KendallRP_MultiIndex.csv
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 4d: Kendall RP — multi-index comparison\n")
cat("══════════════════════════════════════\n")

# ── Indices to analyse ────────────────────────────────────────────────────────
# Compare SPI and SPEI-PM short scales, plus native snow and soil scales.
# Basin fits share the canonical 500-replicate bootstrap and result cache.
MULTI_INDICES <- c(
  "spi_01", "spi_03", "spi_06",
  "spei_01", "spei_03", "spei_06",
  "swei_03", "sspi_01",
  if (length(SOIL_SCALES)) c(sprintf("ssmi_%02d", SOIL_SCALES), sprintf("ssi_sm_%02d", SOIL_SCALES))
)

# Keep only those actually present in master_events, while reporting any missing
# requested label explicitly (e.g., if an input scale was not produced upstream).
available_idx <- unique(master_events$index_label)
requested_multi_indices <- MULTI_INDICES
missing_idx <- setdiff(requested_multi_indices, available_idx)
MULTI_INDICES <- requested_multi_indices[requested_multi_indices %in% available_idx]
if (length(MULTI_INDICES) == 0) {
  cat("  WARNING: none of the requested indices found in master_events — skipping Part 4d\n")
} else {
  if (length(missing_idx))
    cat(sprintf("  NOTE: not found in catalog (skipped): %s\n",
                paste(missing_idx, collapse = ", ")))
  cat(sprintf("  Analysing: %s\n", paste(MULTI_INDICES, collapse = ", ")))
}

# ── Pretty display labels for index_label codes ───────────────────────────────
pretty_index_label <- function(idx_label) {
  sc <- as.integer(sub(".*_([0-9]+)$", "\\1", idx_label))
  if (grepl("^ssi_sm_", idx_label)) return(sprintf("SSI-SM-%d", sc))
  if (grepl("^ssmi_", idx_label)) return(sprintf("SSMI-%d", sc))
  if (grepl("^spei_", idx_label)) return(sprintf("SPEI-PM-%d", sc))
  if (grepl("^swei_", idx_label)) return(sprintf("SWEI-%d", sc))
  if (grepl("^sspi_", idx_label)) return(sprintf("SSPI-%d", sc))
  if (grepl("^spi_", idx_label)) return(sprintf("SPI-%d", sc))
  idx_label
}

# ── Run all indices ───────────────────────────────────────────────────────────
multi_results <- list()
for (idx in MULTI_INDICES) {
  cat(sprintf("\n  --- %s ---\n", toupper(gsub("_","-",idx))))
  multi_results[[idx]] <- get_kendall_rp(idx, n_boot = KENDALL_BOOT_DEFAULT)
  r <- multi_results[[idx]]
  if (!is.null(r$error)) {
    cat(sprintf("  FAILED: %s\n", r$error))
  } else {
    cat(sprintf("  n=%d  tau=%.4f  copula=%s(theta=%.3f)  u0=%.4f v0=%.4f\n",
                r$n_ev, r$tau_emp, r$best_cop, r$theta, r$u0, r$v0))
    cat(sprintf("  t*=%.4f  K_C=%.5f  T_K=%.0f yr  [95%% CI: %.0f-%.0f yr]  (%d/%d conv)\n",
                r$t_star, r$Kc, r$TK, r$ci_lo, r$ci_hi, r$n_boot_conv, KENDALL_BOOT_DEFAULT))
  }
}

# ── Summary table ─────────────────────────────────────────────────────────────
multi_tbl <- do.call(rbind, lapply(multi_results, function(r) {
  data.frame(
    Index       = pretty_index_label(r$index),
    N_events    = r$n_ev    %||% NA,
    Tau_emp     = r$tau_emp %||% NA,
    Best_copula = r$best_cop %||% NA,
    Theta       = r$theta   %||% NA,
    Tau_fitted  = r$tau_fit %||% NA,
    Tgt_dur_mo  = r$tgt_dur %||% NA,
    Tgt_sev     = r$tgt_sev %||% NA,
    u0          = r$u0      %||% NA,
    v0          = r$v0      %||% NA,
    t_star      = r$t_star  %||% NA,
    K_C         = r$Kc      %||% NA,
    T_K_yr      = r$TK      %||% NA,
    CI_lo_yr    = r$ci_lo   %||% NA,
    CI_hi_yr    = r$ci_hi   %||% NA,
    stringsAsFactors = FALSE
  )
}))

csv_path <- file.path(RANKING_DIR, "Table_KendallRP_MultiIndex.csv")
write.csv(multi_tbl, csv_path, row.names = FALSE)
cat(sprintf("\n  Summary table saved: %s\n", basename(csv_path)))
cat("\n  Multi-index Kendall RP summary:\n")
print(multi_tbl[, c("Index","N_events","Tau_emp","Best_copula","Theta",
                    "t_star","K_C","T_K_yr","CI_lo_yr","CI_hi_yr")],
      row.names = FALSE)

# ── Comparison figure ─────────────────────────────────────────────────────────
valid_res <- multi_results[sapply(multi_results, function(r) !is.null(r$TK) && is.finite(r$TK))]

if (length(valid_res) >= 2) {
  
  # ── Panel A: T_K point + CI across indices (log y-axis) ──────────────────
  plot_df <- do.call(rbind, lapply(valid_res, function(r) {
    data.frame(
      index   = pretty_index_label(r$index),
      TK      = r$TK,
      ci_lo   = r$ci_lo,
      ci_hi   = pmin(r$ci_hi, 5e4),   # cap display at 50 000 yr
      tau_emp = r$tau_emp,
      n_ev    = r$n_ev,
      copula  = r$best_cop,
      t_star  = r$t_star,
      Kc      = r$Kc,
      stringsAsFactors = FALSE
    )
  }))
  plot_df$index <- factor(plot_df$index, levels = plot_df$index)
  
  # ── Panel A (WITH CI): T_K point + CI across indices ─────────────────────
  pa_d <- ggplot2::ggplot(plot_df,
                          ggplot2::aes(x = index, y = TK)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      width = 0.25, linewidth = 0.7, colour = COL_HIST) +
    ggplot2::geom_point(size = 4, shape = 18, colour = COL_RECENT) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f yr", TK)),
      vjust = -0.9, size = 2.8, colour = COL_RECENT, fontface = "bold") +
    ggplot2::scale_y_continuous(
      trans  = "log10",
      breaks = c(10, 50, 100, 500, 1000, 5000, 50000),
      labels = scales::comma,
      name   = "Kendall return period (years, log scale)",
      expand = ggplot2::expansion(mult = c(0.1, 0.35))) +
    ggplot2::labs(
      title    = "(a)  Kendall RP point estimate + 95% bootstrap CI by drought index",
      subtitle = "Diamond = point estimate; bar = bootstrap 95% CI (capped at 50,000 yr)",
      x        = NULL) +
    theme_ms3 +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  # ── Panel A (NO CI): T_K point only ──────────────────────────────────────
  pa_d_noci <- ggplot2::ggplot(plot_df,
                               ggplot2::aes(x = index, y = TK)) +
    ggplot2::geom_point(size = 4, shape = 18, colour = COL_RECENT) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f yr", TK)),
      vjust = -0.9, size = 2.8, colour = COL_RECENT, fontface = "bold") +
    ggplot2::scale_y_continuous(
      trans  = "log10",
      breaks = c(10, 50, 100, 500, 1000, 5000, 50000),
      labels = scales::comma,
      name   = "Kendall return period (years, log scale)",
      expand = ggplot2::expansion(mult = c(0.1, 0.35))) +
    ggplot2::labs(
      title    = "(a)  Kendall RP point estimate by drought index",
      subtitle = "Diamond = point estimate only (no bootstrap CI shown)",
      x        = NULL) +
    theme_ms3 +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank()
    )
  
  # ── Panel B: Kendall tau + t* + K_C diagnostic strip ─────────────────────
  long_df <- do.call(rbind, lapply(valid_res, function(r) {
    rbind(
      data.frame(index   = pretty_index_label(r$index),
                 metric  = "Kendall tau\n(duration, severity)",
                 value   = r$tau_emp,
                 stringsAsFactors = FALSE),
      data.frame(index   = pretty_index_label(r$index),
                 metric  = "t* = C(u0,v0)",
                 value   = r$t_star,
                 stringsAsFactors = FALSE),
      data.frame(index   = pretty_index_label(r$index),
                 metric  = "K_C(t*)",
                 value   = r$Kc,
                 stringsAsFactors = FALSE)
    )
  }))
  long_df$index  <- factor(long_df$index,
                           levels = vapply(names(valid_res), pretty_index_label, character(1)))
  long_df$metric <- factor(long_df$metric,
                           levels = c("Kendall tau\n(duration, severity)",
                                      "t* = C(u0,v0)",
                                      "K_C(t*)"))
  
  pb_d <- ggplot2::ggplot(long_df,
                          ggplot2::aes(x = index, y = value,
                                       colour = metric, group = metric)) +
    ggplot2::geom_line(linewidth = 0.7, linetype = "dashed") +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 3)),
                       vjust = -0.8, size = 2.4, show.legend = FALSE) +
    ggplot2::scale_colour_manual(
      values = c("Kendall tau\n(duration, severity)" = "#2166ac",
                 "t* = C(u0,v0)"                      = "#d6604d",
                 "K_C(t*)"                            = "#4dac26"),
      name = NULL) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.02),
      breaks = seq(0, 1, by = 0.2),
      name   = "Value (0-1)",
      expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
    ggplot2::labs(
      title    = "(b)  Copula diagnostics: tau, t*, and K_C by index",
      subtitle = "tau near 0 -> t* near independence -> K_C near 1 -> T_K inflated",
      x        = "Drought index") +
    theme_ms3 +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position    = "bottom"
    )
  
  # ── Shared caption text ────────────────────────────────────────────────────
  cap_base <- paste0(
    "Note: All estimates use 5-family Archimedean copulas (Clayton/Gumbel/Frank/Joe/AMH), ",
    "4 marginal candidates each, tau-initialised MLE. ",
    "T_K driven primarily by K_C proximity to 1, determined by Kendall tau between ",
    "duration and severity.")
  
  # ── Assemble WITH CI ───────────────────────────────────────────────────────
  fig4d <- pa_d / pb_d +
    patchwork::plot_layout(heights = c(1, 1)) +
    patchwork::plot_annotation(
      caption = paste0(cap_base, " 500-rep bootstrap CI shown; upper bound capped at 50,000 yr."),
      theme   = ggplot2::theme(
        plot.caption = ggplot2::element_text(size=7, colour="grey40", hjust=0,
                                             margin=ggplot2::margin(t=6))))
  
  fig4d_pdf <- file.path(RANKING_DIR, "Fig4d_KendallRP_MultiIndex_Comparison.pdf")
  fig4d_png <- file.path(RANKING_DIR, "Fig4d_KendallRP_MultiIndex_Comparison.png")
  tryCatch({
    ggplot2::ggsave(fig4d_pdf, fig4d, width=7.5, height=8.0, units="in", device=if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf)
    cat(sprintf("  Saved (with CI) : %s\n", basename(fig4d_pdf)))
  }, error=function(e) cat(sprintf("  WARNING PDF: %s\n", e$message)))
  tryCatch({
    ggplot2::ggsave(fig4d_png, fig4d, width=7.5, height=8.0, units="in",
                    dpi=300, device="png")
    cat(sprintf("  Saved (with CI) : %s\n", basename(fig4d_png)))
  }, error=function(e) cat(sprintf("  WARNING PNG: %s\n", e$message)))
  
  # ── Assemble WITHOUT CI ───────────────────────────────────────────────────
  fig4d_noci <- pa_d_noci / pb_d +
    patchwork::plot_layout(heights = c(1, 1)) +
    patchwork::plot_annotation(
      caption = paste0(cap_base, " No bootstrap CI shown."),
      theme   = ggplot2::theme(
        plot.caption = ggplot2::element_text(size=7, colour="grey40", hjust=0,
                                             margin=ggplot2::margin(t=6))))
  
  fig4d_noci_pdf <- file.path(RANKING_DIR, "Fig4d_KendallRP_MultiIndex_Comparison_noCI.pdf")
  fig4d_noci_png <- file.path(RANKING_DIR, "Fig4d_KendallRP_MultiIndex_Comparison_noCI.png")
  tryCatch({
    ggplot2::ggsave(fig4d_noci_pdf, fig4d_noci, width=7.5, height=8.0, units="in", device=if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf)
    cat(sprintf("  Saved (no CI)   : %s\n", basename(fig4d_noci_pdf)))
  }, error=function(e) cat(sprintf("  WARNING PDF no-CI: %s\n", e$message)))
  tryCatch({
    ggplot2::ggsave(fig4d_noci_png, fig4d_noci, width=7.5, height=8.0, units="in",
                    dpi=300, device="png")
    cat(sprintf("  Saved (no CI)   : %s\n", basename(fig4d_noci_png)))
  }, error=function(e) cat(sprintf("  WARNING PNG no-CI: %s\n", e$message)))
  
} else {
  cat("  NOTE: fewer than 2 valid results — comparison figure skipped\n")
}


####################################################################################
# PART 5: SPATIAL MAP OF configured recent window DROUGHT (pixel‑level severity)
####################################################################################
cat("\n══════════════════════════════════════\n")
cat(paste0("PART 5: Spatial Analysis of ", RECENT_TAG, " Drought\n"))
cat("══════════════════════════════════════\n")

# Load the authoritative seven-subbasin union (cached loader in utils).
basin_sf <- NULL
basin_vect <- tryCatch(
  load_basin_vect(),
  error = function(e) NULL
)
spatial_map_saved <- FALSE

if (!is.null(basin_vect) && !is.null(most_severe_recent)) {
  cat("  ✓ Basin shapefile loaded\n")
  
  # Get the exact dates of the recent drought from the catalog
  recent_event <- most_severe_recent %>%
    filter(index_label == PRIMARY_INDEX)
  
  if (nrow(recent_event) == 1) {
    start_d <- recent_event$start_date[1]
    end_d   <- recent_event$end_date[1]
    cat(sprintf("  Recent drought period: %s to %s\n", start_d, end_d))
    
    # Parse primary index into type and scale
    idx_type <- sub("_[0-9]+$", "", PRIMARY_INDEX)
    idx_scale <- as.integer(sub(".*_([0-9]+)$", "\\1", PRIMARY_INDEX))
    data_dir <- get_event_index_dir(idx_type)
    
    # Locate all monthly NetCDF files for this index/scale.
    all_files <- get_event_index_files(idx_type, idx_scale)
    
    if (length(all_files) == 0) {
      cat("  ⚠ No NetCDF files found for ", PRIMARY_INDEX, "\n")
    } else {
      cat(sprintf("  Found %d monthly NetCDF files\n", length(all_files)))
      
      # Determine which layers fall within the event period
      ts_list <- list()
      selected_dates <- as.Date(character())
      event_start_month <- as.Date(format(start_d, "%Y-%m-01"))
      event_end_month <- as.Date(format(end_d, "%Y-%m-01"))
      for (f in all_files) {
        r <- terra::rast(f)
        dates <- extract_dates_from_nc(f, terra::nlyr(r))
        if (length(dates) != terra::nlyr(r) || anyNA(dates))
          stop("Invalid NetCDF layer dates in ", basename(f))
        # Compare calendar months (time coordinates need not be the first of a month).
        dates <- as.Date(format(dates, "%Y-%m-01"))
        keep <- which(dates >= event_start_month & dates <= event_end_month)
        if (length(keep) > 0) {
          ts_list[[length(ts_list) + 1L]] <- r[[keep]]
          selected_dates <- c(selected_dates, dates[keep])
        }
      }
      expected_dates <- seq(event_start_month, event_end_month, by = "month")
      if (length(ts_list) && !identical(sort(selected_dates), expected_dates)) {
        cat("  ⚠ Incomplete or duplicate NetCDF event months; pixel map not generated\n")
        ts_list <- list()
      }
      
      if (length(ts_list) == 0) {
        cat("  ⚠ No complete, date-aligned NetCDF data for the event period\n")
      } else {
        # Combine all selected layers into a single SpatRaster
        event_stack <- do.call(c, ts_list)
        cat(sprintf("  Extracted %d months of data for the event\n", terra::nlyr(event_stack)))
        
        # Native fractional-coverage mask: include every cell intersecting
        # the seven-subbasin union, including boundary cells whose centre
        # falls outside the polygon. Project vectors, never resample data.
        basin_on_grid <- if (terra::same.crs(basin_vect, event_stack))
          basin_vect else terra::project(basin_vect, terra::crs(event_stack))
        coverage <- terra::rasterize(basin_on_grid, event_stack[[1]], cover = TRUE)
        coverage <- terra::ifel(coverage > 0, 1, NA)
        event_stack <- terra::mask(event_stack, coverage)
        basin_sf <- sf::st_as_sf(basin_on_grid)
        
        # Compute per-pixel metrics
        r_mean <- terra::app(event_stack, fun = mean, na.rm = TRUE)
        r_min  <- terra::app(event_stack, fun = min, na.rm = TRUE)
        r_dur  <- terra::app(event_stack, fun = function(x) {
          if (!any(is.finite(x))) return(NA_real_)
          sum(x < DROUGHT_ONSET, na.rm = TRUE)
        })
        # A count over an all-NA pixel must not be displayed as zero drought.
        # Explicitly retain the same basin footprint as the mean/minimum maps.
        r_dur <- terra::mask(r_dur, coverage)
        
        # Convert to data frame for ggplot
        df_mean <- as.data.frame(r_mean, xy = TRUE, na.rm = FALSE)
        colnames(df_mean)[3] <- "mean_value"
        df_min  <- as.data.frame(r_min, xy = TRUE, na.rm = FALSE)
        colnames(df_min)[3] <- "min_value"
        df_dur  <- as.data.frame(r_dur, xy = TRUE, na.rm = FALSE)
        colnames(df_dur)[3] <- "months_drought"
        
        # Merge
        df_spatial <- merge(df_mean, df_min[, c("x","y","min_value")], by = c("x","y"))
        df_spatial <- merge(df_spatial, df_dur[, c("x","y","months_drought")], by = c("x","y"))
        
        
        # Create three-panel map
        # ── Fig 5 panel A: Mean index value during the event ──────────────────
        # Values are negative for drought (SPI/SPEI onset threshold ~ -0.5 to -1.0).
        # plasma with direction = -1 maps low (most-negative = most intense) to the
        # bright/yellow end and high (near-zero = mild) to dark purple. We therefore
        # use direction = 1 (default) so that MORE NEGATIVE → darker purple, which
        # is the perceptually intuitive "darker = drier" direction.
        p_mean <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = df_spatial, ggplot2::aes(x = x, y = y, fill = mean_value)) +
          ggplot2::scale_fill_viridis_c(
            option    = "plasma",
            direction = 1,        # low (most negative) = darkest, 0 = yellow → darker = drier
            name      = "Mean Index"
          ) +
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", linewidth = 0.6) +
          ggplot2::coord_sf() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Mean Index Value",
                        subtitle = sprintf("%s to %s", start_d, end_d)) +
          ggplot2::theme(legend.position = "bottom")
        
        # ── Fig 5 panel B: Minimum (peak intensity) during the event ──────────
        # Same logic as panel A: most-negative SPI/SPEI = peak drought = darkest.
        # Using inferno (direction = 1) keeps that mapping.
        p_min <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = df_spatial, ggplot2::aes(x = x, y = y, fill = min_value)) +
          ggplot2::scale_fill_viridis_c(
            option    = "inferno",
            direction = 1,        # most-negative = darkest = most intense
            name      = "Minimum Index"
          ) +
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", linewidth = 0.6) +
          ggplot2::coord_sf() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Minimum Index (Peak Intensity)") +
          ggplot2::theme(legend.position = "bottom")
        
        # ── Fig 5 panel C: Drought duration (months in drought) ───────────────
        # months_drought is a POSITIVE count (0 = never in drought, max = all months).
        # cividis default (direction = 1) goes dark-blue (low) → yellow (high), so
        # pixels with MANY drought months appear LIGHT yellow — the wrong direction.
        # Reversing with direction = -1 makes high counts dark → darker = longer drought.
        p_dur <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = df_spatial, ggplot2::aes(x = x, y = y, fill = months_drought)) +
          ggplot2::scale_fill_viridis_c(
            option    = "cividis",
            direction = -1,       # flip: high month-count = dark, few months = light
            name      = "Months in Drought"
          ) +
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", linewidth = 0.6) +
          ggplot2::coord_sf() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Drought Duration within Event") +
          ggplot2::theme(legend.position = "bottom")
        
        # Combine panels
        p_spatial <- p_mean + p_min + p_dur +
          patchwork::plot_annotation(
            title = paste0(RECENT_TAG, " Drought: Spatial Characteristics"),
            subtitle = sprintf("%s | Pixel-level metrics", event_index_label(PRIMARY_INDEX))
          )
        
        ggplot2::ggsave(file.path(RANKING_DIR, paste0("Fig5_", RECENT_TAG, "_Spatial_Severity.png")),
                        p_spatial, width = 15, height = 6, dpi = 300)
        spatial_map_saved <- TRUE
        cat(paste0("    ✓ Saved: Fig5_", RECENT_TAG, "_Spatial_Severity.png\n"))
      }
    }
  } else {
    cat("  ⚠ Could not isolate recent event for spatial mapping\n")
  }
} else {
  cat("  ⚠ Basin boundary or recent event missing; pixel-level map unavailable\n")
}

# The fallback is independent of the pixel-map branch, including failure to isolate
# one event, missing NetCDF files, and incomplete event-month coverage.
if (!spatial_map_saved && !is.null(basin_vect)) {
  basin_sf <- sf::st_as_sf(basin_vect)
  label_text <- if (nrow(msr_primary) == 1L) {
    sprintf(paste0(RECENT_TAG, " Drought\nSeverity: %.1f\nRank: #%d"),
            msr_primary$severity[1], msr_primary$rank_by_severity[1])
  } else {
    "No recent ranked event available"
  }
  label_xy <- sf::st_coordinates(sf::st_point_on_surface(
    sf::st_union(sf::st_geometry(basin_sf))))[1, 1:2]
  p_spatial <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = basin_sf, fill = "#ecf0f1", color = "#2c3e50", linewidth = 0.8) +
    ggplot2::annotate("text", x = label_xy[1], y = label_xy[2],
                      label = label_text, size = 5, fontface = "bold",
                      color = "#c0392b", hjust = 0.5, vjust = 0.5) +
    ggplot2::labs(
      title = paste0(RECENT_TAG, " Drought Spatial Context (Placeholder)"),
      subtitle = sprintf("Nechako Watershed | Severity ranking relative to %s", ANALYSIS_PERIOD_LABEL)
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
      panel.grid = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(file.path(RANKING_DIR, paste0("Fig5_", RECENT_TAG, "_Spatial_Context_placeholder.png")),
                  p_spatial, width = 10, height = 8, dpi = 300)
  cat("    ✓ Saved placeholder map (pixel-level map unavailable)\n")
} else if (!spatial_map_saved) {
  cat("  ⚠ No basin geometry available for a placeholder map\n")
}

####################################################################################
# PART 5b: configured recent window PIXEL-LEVEL SEVERITY RANK MAP
####################################################################################
cat("\n══════════════════════════════════════\n")
cat(paste0("PART 5b: ", RECENT_TAG, " Pixel-Level Severity Rank Map\n"))
cat("══════════════════════════════════════\n")

# Resolve primary index type and scale
idx_type_p  <- sub("_[0-9]+$", "", PRIMARY_INDEX)
idx_scale_p <- as.integer(sub(".*_([0-9]+)$", "\\1", PRIMARY_INDEX))

# Part 5b reads the time-series RDS written by 13trend_test_ALL.R (not the results
# CSV, which deliberately omits raw time-series to stay lightweight).
ts_rds_file <- file.path(TREND_DIR,
                         sprintf("%s_%02d_timeseries.rds", idx_type_p, idx_scale_p))

rank_map_done <- FALSE

if (file.exists(ts_rds_file)) {
  cat(sprintf("  Loading pixel time-series matrix: %s\n", basename(ts_rds_file)))
  
  ts_obj <- tryCatch(
    readRDS(ts_rds_file),
    error = function(e) {
      cat(sprintf("  ⚠ Could not load time-series RDS: %s\n", e$message))
      NULL
    }
  )
  
  if (!is.null(ts_obj)) {
    ts_mat_p        <- ts_obj$ts_mat
    lon_p           <- ts_obj$lon
    lat_p           <- ts_obj$lat
    record_start_yr <- ts_obj$record_start_year %||% 1950L
    n_tcols         <- ncol(ts_mat_p)
    n_pix_p         <- nrow(ts_mat_p)
    
    cat(sprintf("  Pixels: %d  |  Time columns: %d  |  Record start: %d\n",
                n_pix_p, n_tcols, record_start_yr))
    
    # Event boundaries must be detected on the entire historical record.
    # An onset before the configured recent window may still belong to its event and must not be
    # truncated before calculating duration or severity.
    col_recent_start <- (as.integer(format(RECENT_DROUGHT_START, "%Y")) - record_start_yr) * 12L +
      as.integer(format(RECENT_DROUGHT_START, "%m"))
    col_recent_end <- (as.integer(format(RECENT_DROUGHT_END, "%Y")) - record_start_yr) * 12L +
      as.integer(format(RECENT_DROUGHT_END, "%m"))
    col_recent_start <- max(1L, col_recent_start)
    col_recent_end <- min(n_tcols, col_recent_end)
    
    if (col_recent_start > col_recent_end || col_recent_start > n_tcols) {
      cat("  Recent period lies outside the pixel time series; skipping rank map\n")
    } else {
      rank_vec <- rep(NA_real_, n_pix_p)
      sev_recent <- rep(NA_real_, n_pix_p)
      n_hist <- rep(NA_integer_, n_pix_p)
      min_duration_p <- get_event_min_duration(idx_scale_p)
      
      for (i in seq_len(n_pix_p)) {
        ev_all <- drought_event_table(as.numeric(ts_mat_p[i, ]),
                                      min_duration = min_duration_p,
                                      onset_threshold = DROUGHT_ONSET_EVENT,
                                      termination_threshold = DROUGHT_END_EVENT)
        if (!nrow(ev_all)) next
        n_hist[i] <- nrow(ev_all)
        recent <- ev_all[ev_all$start_index <= col_recent_end &
                           ev_all$end_index >= col_recent_start, , drop = FALSE]
        if (!nrow(recent)) next
        s_recent <- max(recent$severity)
        sev_recent[i] <- s_recent
        # Competition rank: same tie policy as min_rank(desc(severity)).
        rank <- 1L + sum(ev_all$severity > s_recent)
        rank_vec[i] <- 100 * (1 - rank / (nrow(ev_all) + 1L))
      }
      
      # Assemble output data frame
      rank_df <- data.frame(
        lon           = lon_p,
        lat           = lat_p,
        severity_recent = sev_recent,
        pct_rank_recent = rank_vec,
        n_hist_events = n_hist
      )
      write.csv(rank_df,
                file.path(RANKING_DIR, paste0("pixel_rank_", RECENT_TAG, ".csv")),
                row.names = FALSE)
      cat(paste0("  ✓ Saved: pixel_rank_", RECENT_TAG, ".csv\n"))
      
      # Build rasters
      # Native ERA5-Land grid (same template as Scripts 14/15); never a synthetic grid.
      r_base <- terra::rast(file.path(WD_PATH, "monthly_data_direct",
                                      "total_precipitation_monthly.nc"))[[1]]
      
      if (is.null(r_base)) {
        r_base <- terra::rast(
          xmin = min(rank_df$lon) - 0.05, xmax = max(rank_df$lon) + 0.05,
          ymin = min(rank_df$lat) - 0.05, ymax = max(rank_df$lat) + 0.05,
          resolution = 0.1, crs = "EPSG:4326"
        )
        cat("  ℹ Using synthetic raster extent (D13p NC not found)\n")
      }
      
      r_rank <- analysis_rasterize_table(rank_df, r_base, "pct_rank_recent", "lon", "lat")
      r_sev  <- analysis_rasterize_table(rank_df, r_base, "severity_recent", "lon", "lat")
      
      names(r_rank) <- paste0("pct_rank_", RECENT_TAG_UNDERSCORE)
      names(r_sev)  <- paste0("severity_", RECENT_TAG_UNDERSCORE)
      
      terra::writeCDF(
        c(r_rank, r_sev),
        file.path(RANKING_DIR, paste0("pixel_rank_", RECENT_TAG, ".nc")),
        overwrite = TRUE
      )
      cat(paste0("  ✓ Saved: pixel_rank_", RECENT_TAG, ".nc\n"))
      
      # Two-panel PNG map
      df_rank <- as.data.frame(r_rank, xy = TRUE, na.rm = TRUE)
      df_sev  <- as.data.frame(r_sev,  xy = TRUE, na.rm = TRUE)
      colnames(df_rank)[3] <- "pct_rank"
      colnames(df_sev)[3]  <- "severity"
      
      basin_sf_plot <- tryCatch({
        b <- load_basin_vect()
        if (!terra::same.crs(b, r_rank)) b <- terra::project(b, terra::crs(r_rank))
        sf::st_as_sf(b)
      }, error = function(e) {
        warning("Could not prepare Part 5b basin boundary: ", conditionMessage(e))
        NULL
      })
      
      # ──────────────────────────────────────────────────────────────────────
      # WHAT DOES "PERCENTILE RANK" MEAN IN THIS MAP?
      # ──────────────────────────────────────────────────────────────────────
      # For every grid pixel, we extracted the full time-series of the primary
      # drought index from HISTORICAL_START to HISTORICAL_END.
      # Using the Sheffield & Wood (2008) event-identification algorithm
      # (drought onset < DROUGHT_ONSET, termination ≥ DROUGHT_END) we built a
      # catalogue of ALL discrete drought events at that pixel over the entire
      # record.  Each event has a severity score S = mean_intensity × duration
      # (where intensity = DROUGHT_END − index_value, always positive for event months).
      #
      # Recent event selection uses full-record event boundaries and keeps
      # all qualifying events overlapping the target dates. The most severe
      # overlapping event is ranked against that same complete catalog.
      # Percentile = 100 * (1 - rank/(n+1)), where rank is 1 for the
      # most severe event and ties receive their best (minimum) rank.
      # NA denotes no qualifying recent event or no historical events.
      # This is an event-ranking plotting position, NOT an annual exceedance
      # probability or a calibrated annual return period.
      #
      # The colour scale uses plasma with direction = -1 so that high percentile
      # ranks appear DARK (visually prominent) and near-zero ranks appear LIGHT.
      # ──────────────────────────────────────────────────────────────────────
      # WHAT DOES "SEVERITY (I × D)" MEAN IN THIS MAP?
      # ──────────────────────────────────────────────────────────────────────
      # Severity follows Sheffield & Wood (2008, J. Climate):
      #   S = Ī × D
      # where Ī = mean intensity over the event duration and D = duration in months.
      # Intensity for a single month is: I_t = DROUGHT_END − index_value 
      # (positive in every event month; DROUGHT_ONSET and DROUGHT_END are both −0.5 at present).
      #
      # This means a pixel can have high severity because it was:
      #   (a) extremely dry (high I) even for a short time, OR
      #   (b) persistently in drought (high D) even at moderate intensity.
      # The S map therefore captures the COMBINED space-time footprint of the
      # configured recent window event and complements the percentile-rank map.
      #
      # The colour scale uses inferno with direction = -1 so that high severity
      # values appear DARK — matching the convention of all other maps in this
      # analysis where "darker = more extreme drought".
      # ──────────────────────────────────────────────────────────────────────
      
      make_base_map <- function(df, aes_col, fill_label, palette, title_str) {
        p <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = df,
                               ggplot2::aes(x = x, y = y, fill = .data[[aes_col]])) +
          ggplot2::scale_fill_viridis_c(
            option    = palette,
            direction = -1,     # flip palette: HIGH values = DARK, low values = light.
            # For pct_rank (0–100 %): high percentile = dark = most extreme.
            # For severity  (I × D) : high severity  = dark = worst drought.
            name      = fill_label,
            na.value  = "grey90"
          ) +
          ggplot2::coord_sf() +
          ggplot2::labs(title = title_str) +
          ggplot2::theme_bw(base_size = 11) +
          ggplot2::theme(
            plot.title      = ggplot2::element_text(face = "bold"),
            legend.position = "bottom",
            axis.title      = ggplot2::element_blank()
          )
        if (!is.null(basin_sf_plot))
          p <- p + ggplot2::geom_sf(data = basin_sf_plot, fill = NA,
                                    colour = "black", linewidth = 0.6)
        p
      }
      
      p_rank_map <- make_base_map(
        df_rank, "pct_rank",
        "Percentile rank (%)",
        "plasma",
        paste0(RECENT_TAG, " Drought: Percentile rank vs all historical events")
      )
      p_sev_map <- make_base_map(
        df_sev, "severity",
        "Severity (I × D)",
        "inferno",
        paste0(RECENT_TAG, " Drought: Pixel-level S&W severity")
      )
      
      p_rank_combo <- p_rank_map + p_sev_map +
        patchwork::plot_annotation(
          title    = paste0(RECENT_TAG, " Drought — Pixel-Level Rank vs Historical Record"),
          subtitle = sprintf(
            paste0("%s  |  Percentile of ", RECENT_TAG, " severity relative to all events at each pixel  |  %d–%d"),
            event_index_label(PRIMARY_INDEX), HISTORICAL_START, HISTORICAL_END
          ),
          theme = ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(colour = "grey40", size = 10, hjust = 0.5)
          )
        )
      
      ggplot2::ggsave(
        file.path(RANKING_DIR, paste0("Fig6_", RECENT_TAG, "_Pixel_Rank_Map.png")),
        p_rank_combo, width = 14, height = 7, dpi = 300
      )
      cat(paste0("  ✓ Saved: Fig6_", RECENT_TAG, "_Pixel_Rank_Map.png\n"))
      rank_map_done <- TRUE
      
      cat(sprintf(
        "  Basin-wide: median percentile = %.1f%%  |  %.1f%% of pixels at ≥ 75th pct\n",
        median(rank_vec, na.rm = TRUE),
        100 * mean(rank_vec >= 75, na.rm = TRUE)
      ))
    }
  }
} else {
  cat(sprintf("  ⚠ Time-series RDS not found: %s\n  → Run 13trend_test_ALL.R first, then re-run this.\n",
              basename(ts_rds_file)))
}

if (!rank_map_done)
  cat("  ℹ Pixel rank map not produced (see warnings above)\n")

####################################################################################
# PART 6: EXPORT EXCEL SUMMARY FOR PRESENTATION
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 6: Creating Excel Summary\n")
cat("══════════════════════════════════════\n")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Executive_Summary")
openxlsx::addWorksheet(wb, "Ranked_Events_SW08")
openxlsx::addWorksheet(wb, "Duration_Class_Stats")
openxlsx::addWorksheet(wb, "Return_Period_Table")

# Executive Summary
exec_summary <- data.frame(
  Question = c(
    "Q1: How common are droughts?",
    "Q2: Spatial characteristics?",
    "Q3: Are droughts increasing?",
    paste0("Q4: How does ", RECENT_TAG, " compare?")
  ),
  Answer = c(
    sprintf("Short-term: %d events (1 per %.1f yrs) | Long-term: %d events (1 per %.1f yrs)",
            sum(primary_events$duration_class == "D3-6 (Short-term)"),
            n_years_record / sum(primary_events$duration_class == "D3-6 (Short-term)"),
            sum(primary_events$duration_class == "D13+ (Long-term)"),
            n_years_record / sum(primary_events$duration_class == "D13+ (Long-term)")),
    paste0("See Fig5 for pixel-level severity map (mean, minimum, duration during ", RECENT_TAG, ")."),
    "Refer to 13trend_test_ALL.R for Mann-Kendall trend analysis on frequency/intensity.",
    sprintf("Rank #%d of %d events (%.1fth percentile, rank-based recurrence proxy %.1f yr; not annual exceedance)",
            msr_primary$rank_by_severity[1] %||% NA,
            nrow(primary_events),
            msr_primary$percentile_rank[1] %||% NA,
            msr_primary$return_period_years[1] %||% NA)
  ),
  stringsAsFactors = FALSE
)

openxlsx::writeData(wb, "Executive_Summary", exec_summary)
openxlsx::addStyle(wb, "Executive_Summary",
                   openxlsx::createStyle(textDecoration = "bold", fontSize = 12),
                   rows = 1, cols = 1:2)
openxlsx::setColWidths(wb, "Executive_Summary", cols = 1:2, widths = c(30, 80))

# Ranked Events (top 20)
top_20 <- ranked_catalog %>%
  filter(index_label == PRIMARY_INDEX) %>%
  arrange(rank_by_severity) %>%
  head(20) %>%
  dplyr::select(-index_label)

openxlsx::writeData(wb, "Ranked_Events_SW08", top_20)
openxlsx::setColWidths(wb, "Ranked_Events_SW08", cols = 1:ncol(top_20), widths = 15)

# Duration Class Statistics
openxlsx::writeData(wb, "Duration_Class_Stats", duration_summary)
openxlsx::setColWidths(wb, "Duration_Class_Stats", cols = 1:ncol(duration_summary), widths = 20)

# Return Period Table
if (exists("rp_table") && nrow(rp_table) > 0) {
  openxlsx::writeData(wb, "Return_Period_Table", as.data.frame(rp_table))
  openxlsx::addStyle(
    wb, "Return_Period_Table",
    openxlsx::createStyle(textDecoration = "bold", fontSize = 11),
    rows = 1, cols = seq_len(ncol(rp_table))
  )
  openxlsx::setColWidths(wb, "Return_Period_Table",
                         cols = seq_len(ncol(rp_table)),
                         widths = c(18, 24, 10, 16, 12, 12, 16, 10, 10, 10))
  cat("✓ Return_Period_Table sheet written\n")
}

excel_file <- file.path(RANKING_DIR, "Drought_Ranking_Summary.xlsx")
openxlsx::saveWorkbook(wb, excel_file, overwrite = TRUE)
cat(sprintf("✓ Saved: %s\n", basename(excel_file)))

####################################################################################
# PART 4e – MANUSCRIPT TABLE 1a
#
# Table 1a compares short-scale meteorological and native-scale snow drought
# event rankings. SPI/SPEI use the available 1- and 3-month scales; SWEI uses
# scale 3 and SSPI uses scale 1.
# Output: Table_1a_ShortScale_Index_TopEvents.csv and workbook companion sheet.
#
# COLUMN DEFINITIONS (Table 1a, the only manuscript event ranking table):
#   Rank              — severity rank over ALL durations (D1-2 included in the ranking)
#   Rank_D3plus       — severity rank within D3+ events only; the table is sorted by this
#   Start_date        — onset month of the event
#   End_date          — termination month
#   Start_year        — calendar year of onset (used as the event label)
#   Duration_months   — number of consecutive sub-threshold months
#   Duration_class    — D3-6 Short-term / D7-12 Medium-term / D13+ Long-term
#   Mean_intensity    —  mean of (DROUGHT_END − index) over event months [index units]
#   Severity          — Mean_intensity × Duration_months [index-units·months]
#   Return_period_yr  — legacy schema; (record_years+1)/event rank is a
#                       recurrence proxy, NOT an annual return period when
#                       events/year > 1; retained for workbook compatibility.
#   Percentile_rank   — 100 × (1 − rank/(n+1)) [%]
#   In_Recent_Period  — logical flag; TRUE if the event overlaps the recent window
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 4e: Manuscript Table 1a\n")
cat("══════════════════════════════════════\n")

# ── Helper: extract top-N events for a set of index labels, add display columns ──
make_ms_table <- function(catalog, index_labels, n_top = 10L) {
  sub <- catalog %>%
    dplyr::filter(index_label %in% index_labels,
                  duration_class != "D1-2 (Very short)") %>%
    dplyr::group_by(index_label) %>%
    dplyr::mutate(rank_d3plus = dplyr::min_rank(dplyr::desc(severity))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(index_label, rank_d3plus)
  
  # Pretty index labels, including the native-scale snow indices.
  sub$Index <- vapply(sub$index_label, function(lbl) {
    scale_i <- as.integer(base::sub(".*_([0-9]+)$", "\\1", lbl))
    type_i <- base::sub("_[0-9]+$", "", lbl)
    switch(type_i,
           spi = sprintf("SPI-%d", scale_i),
           spei = sprintf("SPEI-PM-%d", scale_i),
           swei = sprintf("SWEI-%d", scale_i),
           sspi = sprintf("SSPI-%d", scale_i),
           lbl
    )
  }, character(1L))
  
  # Take top n_top per index_label
  out <- sub %>%
    dplyr::group_by(index_label) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      Index,
      Rank          = rank_by_severity,
      Start_date    = start_date,
      End_date      = end_date,
      Start_year    = start_year,
      Duration_months = duration_months,
      Duration_class  = duration_class,
      Mean_intensity  = mean_int,
      Severity        = severity,
      Return_period_yr = return_period_years,
      Percentile_rank  = percentile_rank,
      In_Recent_Period     = is_recent,
      Rank_D3plus = rank_d3plus
    ) %>%
    dplyr::mutate(
      Mean_intensity   = round(Mean_intensity,   6L),
      Severity         = round(Severity,         6L),
      Return_period_yr = round(Return_period_yr, 2L),
      Percentile_rank  = round(Percentile_rank,  2L)
    )
  out
}

# ── Styled Excel sheet writer ─────────────────────────────────────────────────
write_ms_table_sheet <- function(wb, sheet_name, tbl) {
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, tbl)
  
  hdr_style <- openxlsx::createStyle(
    textDecoration = "bold", fgFill = "#BDD7EE",
    halign = "center", border = "Bottom",
    borderColour = "#2E75B6", borderStyle = "medium")
  openxlsx::addStyle(wb, sheet_name, hdr_style,
                     rows = 1L, cols = seq_len(ncol(tbl)), stack = TRUE)
  
  # Highlight configured recent window rows in yellow
  for (rw in which(tbl$In_Recent_Period)) {
    openxlsx::addStyle(wb, sheet_name,
                       openxlsx::createStyle(fgFill = "#FFEB9C"),
                       rows = rw + 1L, cols = seq_len(ncol(tbl)),
                       stack = TRUE)
  }
  
  openxlsx::setColWidths(wb, sheet_name,
                         cols = seq_len(ncol(tbl)),
                         widths = c(18, 6, 12, 12, 10, 16, 22, 16, 12, 16, 15, 13, 12))
  invisible(wb)
}

tryCatch({
  
  # ── TABLE 1a: short-scale meteorological and native-scale snow indices ──────
  cat("  Building Table 1a (short-scale index comparison)...\n")
  
  table_1a_requested <- c(
    "spi_01", "spi_03",
    "spei_01", "spei_03",
    "swei_03", "sspi_01"
  )
  table_1a_present <- table_1a_requested[
    table_1a_requested %in% unique(ranked_catalog$index_label)
  ]
  
  if (length(table_1a_present) == 0L) {
    cat("  WARNING: no configured short-scale indices in catalog — Table 1a skipped\n")
  } else {
    tbl_1a <- make_ms_table(ranked_catalog, table_1a_present, n_top = 10L)
    
    csv_1a <- file.path(RANKING_DIR, "Table_1a_ShortScale_Index_TopEvents.csv")
    write.csv(tbl_1a, csv_1a, row.names = FALSE)
    cat(sprintf("  \u2713 Table 1a saved: %s  (%d rows)\n",
                basename(csv_1a), nrow(tbl_1a)))
    
    # Console preview
    cat("\n  Table 1a preview (top 3 per index):\n")
    preview_1a <- tbl_1a %>%
      dplyr::group_by(Index) %>%
      dplyr::slice_head(n = 3L) %>%
      dplyr::ungroup() %>%
      dplyr::select(Index, Rank, Start_date, End_date,
                    Duration_months, Severity, In_Recent_Period)
    print(as.data.frame(preview_1a), row.names = FALSE)
    
    # Add to Excel workbook
    wb_1a_path <- file.path(RANKING_DIR, "Drought_Ranking_Summary.xlsx")
    if (file.exists(wb_1a_path)) {
      wb_loaded <- tryCatch(openxlsx::loadWorkbook(wb_1a_path),
                            error = function(e) NULL)
      if (!is.null(wb_loaded)) {
        # Remove sheet if it exists (overwrite)
        stale_sheets <- intersect(
          c("Table_1a_SPI_TopEvents", "Table_1a_ShortScale_Indices"),
          names(wb_loaded)
        )
        for (stale_sheet in stale_sheets) {
          openxlsx::removeWorksheet(wb_loaded, stale_sheet)
        }
        write_ms_table_sheet(wb_loaded, "Table_1a_ShortScale_Indices", tbl_1a)
        openxlsx::saveWorkbook(wb_loaded, wb_1a_path, overwrite = TRUE)
        cat("  Table_1a_ShortScale_Indices sheet added to Drought_Ranking_Summary.xlsx\n")
      }
    }
  }
  
  cat("\n\u2550\u2550 PART 4e complete \u2550\u2550\n")
  
}, error = function(e) {
  cat(sprintf("  \u26a0 Part 4e failed: %s\n", e$message))
})

####################################################################################
# PART 4f – RECENT DROUGHT EVENT BUBBLE PLOT
#
# Purpose:
#   Publication-oriented visualization of the principal recent event for all
#   configured indices; snow and soil-moisture indices appear at their native scales only.
#
# Encodings:
#   x-position       = event timing; horizontal segment spans start to end date
#   y-position       = accumulation scale (1, 3, 6, 12, 24, 36 months)
#   bubble area      = cumulative severity (S)
#   bubble colour    = historical severity rank class within the corresponding
#                      index × scale event catalog: Rank 1, Rank 2, Rank 3–5,
#                      Rank >5 (no rank number or severity text is printed
#                      inside the bubble)
#   text label       = severity value, placed directly above each bubble (centred
#                      on its anchor date) so it stays in its own scale row
#
# The plot uses the authoritative master_events object generated above. No
# event detection, severity calculation, or ranking is duplicated here.
# For each index × scale, the most severe event overlapping configured recent window is shown.
# This intentionally includes scale-1 events, even when their duration is D1–2,
# because scale 1 has a two-month minimum event-duration gate in this framework.
#
# Note on long-duration events: an event's bubble is anchored at the midpoint
# of its overlap with the configured recent-window bounds display window (not the midpoint
# of its full duration), and the view is zoomed with coord_cartesian() rather
# than scale limits. This ensures a long event (e.g., a 36-month-scale event
# that began well before the recent window) still gets a bubble drawn on-screen instead of
# being silently dropped.
####################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4f: Recent Drought Event Bubble Plot\n")
cat("══════════════════════════════════════════════════════════════\n")

bubble_plot_path_png <- file.path(
  RANKING_DIR,
  "Fig_Recent_Drought_Event_Bubbles.png"
)
bubble_plot_path_pdf <- file.path(
  RANKING_DIR,
  "Fig_Recent_Drought_Event_Bubbles.pdf"
)

# Include all supported formulations in the recent-event diagnostic. Snow and
# soil-moisture indices appear only at their native scales.
bubble_index_map <- c(spi = "(a) SPI", spei = "(b) SPEI-PM",
                      swei = "(c) SWEI", sspi = "(d) SSPI", ssmi = "(e) SSMI", ssi_sm = "(f) SSI-SM")
if (!length(SOIL_SCALES))
  bubble_index_map <- bubble_index_map[setdiff(names(bubble_index_map), c("ssmi", "ssi_sm"))]

bubble_scale_levels <- ANALYSIS_MANUSCRIPT_SCALES

# Visible x-axis window for the figure. Defined once and reused for both the
# coordinate zoom (below) and for clamping each bubble's on-screen anchor
# date. This matters for long-duration events (e.g., a 36-month-scale event
# that began well before the recent window): its true midpoint can fall before this
# window's lower bound, and previously the plot used scale_x_date(limits=...),
# which silently drops (sets to NA, then removes) any point whose x falls
# outside the given limits -- causing that bubble to vanish entirely even
# though the event does overlap configured recent window. The fix below (a) clamps each
# bubble's anchor date into this window instead of relying on scale limits,
# and (b) uses coord_cartesian() to zoom the view without dropping data.
bubble_x_display_min <- seq(as.Date(format(RECENT_DROUGHT_START, "%Y-%m-01")),
                            by = "-1 month", length.out = 2L)[2L]
bubble_x_display_max <- seq(as.Date(format(RECENT_DROUGHT_END, "%Y-%m-01")),
                            by = "month", length.out = 2L)[2L]
# Explicit yearly breaks that EXCLUDE the window end (window end). That break
# sat on the right panel edge, so its "2026" label overhung into the gap and
# collided with the first recent-year label at the left edge of the next panel.
bubble_x_breaks <- seq(as.Date(sprintf("%04d-01-01", RECENT_YEAR_START)), bubble_x_display_max, by = "1 year")
bubble_x_breaks <- bubble_x_breaks[bubble_x_breaks < bubble_x_display_max]

bubble_rank_colours <- c(
  "Rank 1"   = "#C91F2D",
  "Rank 2"   = "#F47C20",
  "Rank 3–5" = "#F5B82E",
  "Rank >5"  = "#54A9B2",
  "Censored / unranked" = "#777777"
)

# Human-readable display label from the authoritative index_type field.
# Fall back to index_label parsing only if index_type is unavailable.
format_bubble_index <- function(index_type, index_label) {
  idx <- if (!is.na(index_type) && nzchar(index_type)) {
    tolower(index_type)
  } else {
    sub("_[0-9]+$", "", index_label)
  }
  unname(bubble_index_map[idx] %||% NA_character_)
}

# Select exactly one event per index × scale from the normal SW08 catalog:
# the most severe qualifying event overlapping configured recent window.
completed_bubble_events <- master_events %>%
  dplyr::filter(
    index_type %in% names(bubble_index_map),
    scale %in% bubble_scale_levels,
    is_recent
  ) %>%
  dplyr::arrange(
    index_type,
    scale,
    dplyr::desc(severity),
    start_date,
    end_date
  ) %>%
  dplyr::group_by(index_type, scale) %>%
  dplyr::slice_head(n = 1L) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    ongoing = right_censored %in% TRUE,
    censored_at_record_end = right_censored %in% TRUE
  )

# Figure-only fallback: represent open runs shorter than the official
# scale-dependent minimum duration. They are *not* qualifying historical
# events and must not enter rank or return-period calculations. This uses
# runs already collected during initial producer-catalog construction and
# does not rewrite those catalog CSVs.
missing_bubble_groups <- dplyr::bind_rows(lapply(INDICES_TO_ANALYZE, function(z) {
  data.frame(index_type = z$index, scale = as.integer(z$scales))
})) %>%
  dplyr::filter(index_type %in% names(bubble_index_map),
                scale %in% bubble_scale_levels) %>%
  dplyr::anti_join(dplyr::distinct(completed_bubble_events, index_type, scale),
                   by = c("index_type", "scale"))

fallback_events <- NULL
if (length(open_recent_events) && nrow(missing_bubble_groups)) {
  open_candidates <- dplyr::bind_rows(open_recent_events)
  fallback_events <- open_candidates %>%
    dplyr::semi_join(missing_bubble_groups, by = c("index_type", "scale")) %>%
    dplyr::mutate(
      rank_by_severity = NA_integer_,
      return_period_years = NA_real_,
      percentile_rank = NA_real_,
      right_censored = TRUE,
      ongoing = TRUE,
      censored_at_record_end = TRUE,
      figure_only_fallback = TRUE
    )
  if (!nrow(fallback_events)) fallback_events <- NULL
}

# Normal qualifying events always take precedence. Fallback events exist
# only for index × scale combinations absent from the normal recent catalog.
bubble_events <- dplyr::bind_rows(
  completed_bubble_events %>%
    dplyr::mutate(figure_only_fallback = FALSE),
  fallback_events
) %>%
  dplyr::arrange(
    index_type,
    scale,
    dplyr::desc(severity),
    start_date,
    end_date
  ) %>%
  dplyr::group_by(index_type, scale) %>%
  dplyr::slice_head(n = 1L) %>%
  dplyr::ungroup()


if (!nrow(bubble_events)) {
  cat("  ⚠ No recent events available; bubble plot skipped.\n")
} else {
  # Preserve the six manuscript scales in fixed vertical order.
  bubble_events <- bubble_events %>%
    dplyr::mutate(
      scale_factor = match(as.integer(scale), bubble_scale_levels),
      panel = factor(
        vapply(
          seq_len(n()),
          function(i) format_bubble_index(index_type[i], index_label[i]),
          character(1L)
        ),
        levels = unname(bubble_index_map)
      ),
      rank_class = dplyr::case_when(
        is.na(rank_by_severity) ~ "Censored / unranked",
        rank_by_severity == 1L ~ "Rank 1",
        rank_by_severity == 2L ~ "Rank 2",
        rank_by_severity >= 3L & rank_by_severity <= 5L ~ "Rank 3–5",
        TRUE ~ "Rank >5"
      ),
      # Anchor date for the bubble itself. Clamped to the intersection of the
      # event span with the visible plot window so the bubble always lands
      # inside the panel -- rather than the raw event midpoint, which can lie
      # outside the visible window (and previously got silently dropped by
      # scale_x_date(limits=...)) for long events that began well before 2022.
      overlap_start = pmax(start_date, bubble_x_display_min),
      overlap_end   = pmin(end_date,   bubble_x_display_max),
      midpoint      = overlap_start + floor(as.numeric(overlap_end - overlap_start) / 2),
      severity_label = sprintf("%.2f", severity)
    )
  
  
  # Deterministic severity-label placement:
  # every severity value is placed directly ABOVE its own bubble,
  # centered on the bubble's midpoint. This prevents labels from moving
  # into adjacent accumulation-scale rows and keeps the 24- and 36-month
  # labels visually aligned with their respective bubbles.
  bubble_events$severity_label_y <-
    bubble_events$scale_factor + 0.48
  
  # Save the exact data used for the figure so that the manuscript graphic is
  # auditable without reconstructing the plot from the rendered image.
  bubble_export <- bubble_events %>%
    dplyr::transmute(
      Index = unname(bubble_index_map[index_type]),
      Scale_months = as.integer(scale),
      Start_date = start_date,
      End_date = end_date,
      Duration_months = as.integer(duration_months),
      Severity = round(severity, 6),
      Rank = as.integer(rank_by_severity),
      Rank_class = rank_class,
      Ongoing = ongoing %in% TRUE,
      Censored_at_record_end = censored_at_record_end %in% TRUE,
      Figure_only_fallback = figure_only_fallback %in% TRUE
    ) %>%
    dplyr::arrange(Index, Scale_months)
  
  bubble_csv_path <- file.path(
    RANKING_DIR,
    "Fig_Recent_Drought_Event_Bubbles_Data.csv"
  )
  write.csv(bubble_export, bubble_csv_path, row.names = FALSE)
  
  # scale_size_area maps the data directly to bubble area, with zero -> zero.
  severity_legend_values <- c(2, 5, 10, 20, 40)
  
  # Small vertical headroom above the top row of bubbles.
  bubble_y_top <- length(bubble_scale_levels) + 1.15  # extra headroom for the severity label above the top row
  
  bubble_has_minor_ticks <- utils::packageVersion("ggplot2") >= "3.5.0"
  p_recent_bubbles <- ggplot2::ggplot(
    bubble_events,
    ggplot2::aes(x = midpoint, y = scale_factor)
  ) +
    # Event-duration line: exact temporal extent from onset to termination.
    ggplot2::geom_segment(
      ggplot2::aes(
        x = start_date,
        xend = end_date,
        y = scale_factor,
        yend = scale_factor
      ),
      linewidth = 0.55,
      colour = "#222222",
      lineend = "round"
    ) +
    # End caps make the temporal extent easier to infer at manuscript size.
    ggplot2::geom_segment(
      ggplot2::aes(
        x = start_date,
        xend = start_date,
        y = scale_factor - 0.09,
        yend = scale_factor + 0.09
      ),
      linewidth = 0.55,
      colour = "#222222"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = end_date,
        xend = end_date,
        y = scale_factor - 0.09,
        yend = scale_factor + 0.09
      ),
      linewidth = 0.55,
      colour = "#222222"
    ) +
    # Bubble area is proportional to severity. Rank is conveyed only through
    # fill colour (see rank_class / bubble_rank_colours); no text is printed
    # inside the bubble itself.
    ggplot2::geom_point(
      ggplot2::aes(
        size = severity,
        fill = rank_class
      ),
      shape = 21,
      colour = "#FFFFFF",
      stroke = 0.7,
      alpha = 0.97
    ) +
    # Black ring identifies right-censored events that nevertheless satisfy
    # the duration criterion and retain their historical severity rank.
    ggplot2::geom_point(
      data = function(d) d[!is.na(d$rank_by_severity) & d$ongoing %in% TRUE, , drop = FALSE],
      ggplot2::aes(size = severity), shape = 21,
      fill = NA, colour = "black", stroke = 1.15, show.legend = FALSE
    ) +
    # Severity value directly above its own bubble (no repel), so labels stay in their scale row.
    ggplot2::geom_label(
      ggplot2::aes(
        x = midpoint,
        y = severity_label_y,
        label = severity_label
      ),
      colour       = "#222222",
      fill         = "white",
      fontface     = "bold",
      size         = 2.75,
      linewidth   = 0,
      label.padding = grid::unit(0.06, "lines"),
      hjust        = 0.5,
      vjust        = 0.5
    ) +
    ggplot2::facet_wrap(
      ~ panel,
      nrow = 1,
      scales = "fixed",
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(bubble_scale_levels),
      labels = bubble_scale_levels,
      expand = c(0, 0),
      name = "Accumulation scale (months)"
    ) +
    ggplot2::scale_x_date(
      breaks = bubble_x_breaks,
      date_labels = "%Y",
      date_minor_breaks = "1 month",
      expand = c(0.015, 0.015),
      name = "Year",
      # Minor ticks require ggplot2 >= 3.5.0; older versions use major ticks.
      guide = if (bubble_has_minor_ticks)
        ggplot2::guide_axis(minor.ticks = TRUE) else ggplot2::guide_axis()
    ) +
    # Zoom to the same display window previously passed to scale_x_date's
    # `limits`. Unlike scale limits (which convert out-of-range values to NA
    # and DROP that entire geom row -- the cause of the missing 36-month-scale
    # bubbles for events that began before this window), coord_cartesian()
    # only clips the *visual* rendering: an event's duration segment can start
    # before the window and still be drawn (clipped at the panel edge), and
    # every bubble/label -- already anchored inside this window via the
    # clamped `midpoint` computed above -- is guaranteed to be visible.
    ggplot2::coord_cartesian(
      xlim = c(bubble_x_display_min, bubble_x_display_max),
      ylim = c(0.55, bubble_y_top),
      clip = "off"
    ) +
    ggplot2::scale_fill_manual(
      values = bubble_rank_colours,
      breaks = names(bubble_rank_colours),
      name = "Rank / status – colour legend",
      drop = FALSE
    ) +
    ggplot2::scale_size_area(
      name = "Cumulative severity (S)",
      breaks = severity_legend_values,
      max_size = 18.5,
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          fill = "#777777",
          colour = "#777777",
          alpha = 1
        )
      )
    ) +
    ggplot2::labs(
      title = paste0("Recent Drought Events in the Nechako River Basin (", RECENT_LABEL, ")"),
      subtitle = "Most severe recent event at each accumulation scale for each drought-index formulation",
      caption = paste0(
        "Each bubble represents one drought event identified independently for the corresponding accumulation scale and index. ",
        "Bubble area is proportional to cumulative severity (S), whose value is labelled above each bubble; bubble colour indicates ",
        "historical severity rank class among ", ANALYSIS_PERIOD_LABEL, " events for that index and scale. Horizontal segments show event extent. ",
        "A black ring denotes a ranked, right-censored event; grey denotes a subminimum-duration open run shown only as an unranked fallback."
      )
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          shape = 21,
          size = 5.0,
          colour = "#FFFFFF",
          alpha = 1
        )
      )
    ) +
    ggplot2::theme_minimal(base_size = 11, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 15,
        hjust = 0.5,
        margin = ggplot2::margin(b = 4)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10.5,
        hjust = 0.5,
        margin = ggplot2::margin(b = 9)
      ),
      plot.caption = ggplot2::element_text(
        size = 8.5,
        hjust = 0,
        lineheight = 1.05,
        margin = ggplot2::margin(t = 10)
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 13,
        colour = "#222222",
        margin = ggplot2::margin(t = 2, b = 5)
      ),
      strip.background = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        face = "bold",
        size = 10.5,
        margin = ggplot2::margin(t = 7)
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold",
        size = 10.5,
        margin = ggplot2::margin(r = 7)
      ),
      axis.text.x = ggplot2::element_text(
        size = 9.5,
        colour = "#222222"
      ),
      axis.text.y = ggplot2::element_text(
        size = 10,
        colour = "#222222"
      ),
      # Monthly minor tick marks on the x-axis (no minor tick labels are drawn).
      axis.ticks.x = ggplot2::element_line(
        colour = "#444444",
        linewidth = 0.4
      ),
      panel.grid.major.x = ggplot2::element_line(
        colour = "#E5E5E5",
        linewidth = 0.35
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "#D9D9D9",
        linewidth = 0.45,
        linetype = "dashed"
      ),
      panel.border = ggplot2::element_rect(
        colour = "#444444",
        fill = NA,
        linewidth = 0.55
      ),
      panel.spacing.x = grid::unit(1.4, "lines"),  # room for edge tick labels
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 9.5
      ),
      legend.text = ggplot2::element_text(
        size = 9
      ),
      legend.key.width = grid::unit(1.0, "cm"),
      legend.margin = ggplot2::margin(t = 2, b = 2),
      plot.margin = ggplot2::margin(8, 10, 8, 10)
    )
  
  if (bubble_has_minor_ticks) {
    p_recent_bubbles <- p_recent_bubbles + ggplot2::theme(
      axis.minor.ticks.x.bottom = ggplot2::element_line(
        colour = "#888888", linewidth = 0.3),
      axis.minor.ticks.length = ggplot2::rel(0.5)
    )
  }
  
  # Save both vector PDF and high-resolution PNG for manuscript use.
  ggplot2::ggsave(
    bubble_plot_path_pdf,
    plot = p_recent_bubbles,
    width = 1.75 * length(bubble_index_map) + 3.25,
    height = 6.7,
    units = "in",
    device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  )
  ggplot2::ggsave(
    bubble_plot_path_png,
    plot = p_recent_bubbles,
    width = 1.75 * length(bubble_index_map) + 3.25,
    height = 6.7,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  cat(sprintf("  ✓ Bubble plot PDF: %s\n", basename(bubble_plot_path_pdf)))
  cat(sprintf("  ✓ Bubble plot PNG: %s\n", basename(bubble_plot_path_png)))
  cat(sprintf("  ✓ Bubble plot data: %s\n", basename(bubble_csv_path)))
  
  # Explicitly print the scale-1 ranks used in the figure. This provides a
  # simple QA check for the corrected one-month ranks requested for the figure.
  scale1_qa <- bubble_events %>%
    dplyr::filter(scale == 1L) %>%
    dplyr::select(index_type, severity, duration_months, rank_by_severity,
                  ongoing, censored_at_record_end, start_date, end_date) %>%
    dplyr::arrange(index_type)
  cat("\n  Scale-1 bubble QA:\n")
  print(as.data.frame(scale1_qa), row.names = FALSE)
  
  scale36_qa <- bubble_events %>%
    dplyr::filter(scale == 36L) %>%
    dplyr::select(index_type, severity, duration_months, rank_by_severity,
                  ongoing, censored_at_record_end, start_date, end_date) %>%
    dplyr::arrange(index_type)
  cat("\n  Scale-36 bubble QA (including right-censored ongoing events):\n")
  print(as.data.frame(scale36_qa), row.names = FALSE)
}

####################################################################################

####################################################################################
####################################################################################
# PART 4g: KENDALL-BASED JOINT RETURN PERIOD — ALL MANUSCRIPT SCALES
#
# This section does NOT implement a second Kendall model. It reuses the canonical
# run_kendall_rp() engine defined in Part 4c through get_kendall_rp().
#
# The canonical engine is therefore identical for:
#   • manuscript (headline-index) Kendall RP (Part 4c);
#   • multi-index comparison (Part 4d);
#   • all six manuscript scales for SPI and SPEI-PM, plus SWEI-03,
#     SSPI-01 and the soil-moisture indices at their native scales.
#
####################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4g: Kendall RP — all manuscript scales\n")
cat("══════════════════════════════════════════════════════════════\n")

if (!exists("master_events") || !nrow(master_events))
  stop("Part 4g requires master_events from Parts 1–2.")

KENDALL_ALL_SCALES <- as.integer(MANUSCRIPT_SCALES)
# Kendall bootstrap covers SPI and SPEI-PM plus snow and
# soil-moisture indices at their valid accumulation scales.
KENDALL_ALL_INDICES <- c(
  paste0("spei_", sprintf("%02d", KENDALL_ALL_SCALES)),
  paste0("spi_", sprintf("%02d", KENDALL_ALL_SCALES)),
  sprintf("swei_%02d", SWEI_SCALE),
  sprintf("sspi_%02d", SSPI_SCALE),
  if (length(SOIL_SCALES)) c(sprintf("ssmi_%02d", SOIL_SCALES), sprintf("ssi_sm_%02d", SOIL_SCALES))
)
KENDALL_ALL_INDICES <- KENDALL_ALL_INDICES[
  KENDALL_ALL_INDICES %in% unique(master_events$index_label)
]

if (!length(KENDALL_ALL_INDICES)) {
  cat("  WARNING: no manuscript index/scale labels found — all-scale Kendall analysis skipped.\n")
} else {
  cat("  Analysing with canonical engine: ",
      paste(KENDALL_ALL_INDICES, collapse = ", "), "\n", sep = "")
  
  all_t0 <- proc.time()[["elapsed"]]
  n_all <- length(KENDALL_ALL_INDICES)
  kendall_all_results <- lapply(seq_along(KENDALL_ALL_INDICES), function(ii) {
    idx <- KENDALL_ALL_INDICES[[ii]]
    elapsed <- proc.time()[["elapsed"]] - all_t0
    cat(sprintf(
      "\n  Kendall RP %d/%d: %s | total elapsed %.1fs\n",
      ii, n_all, idx, elapsed
    ))
    flush.console()
    ans <- get_kendall_rp(idx, n_boot = KENDALL_BOOT_DEFAULT)
    elapsed2 <- proc.time()[["elapsed"]] - all_t0
    cat(sprintf(
      "  Completed %d/%d: %s | total elapsed %.1fs\n",
      ii, n_all, idx, elapsed2
    ))
    flush.console()
    ans
  })
  names(kendall_all_results) <- KENDALL_ALL_INDICES
  
  successful <- Filter(
    function(x) is.list(x) && is.null(x$error) &&
      is.data.frame(x$rows) && nrow(x$rows) > 0L,
    kendall_all_results
  )
  
  failed <- Filter(
    function(x) is.list(x) && !is.null(x$error),
    kendall_all_results
  )
  
  if (length(failed)) {
    cat("\n  Kendall analyses unavailable:\n")
    for (nm in names(failed))
      cat("    ", nm, ": ", failed[[nm]]$error, "\n", sep = "")
  }
  
  if (length(successful)) {
    kendall_rp_all_scales <- dplyr::bind_rows(
      lapply(successful, function(x) x$rows)
    )
    
    kendall_rp_all_scales$index <- dplyr::case_when(
      grepl("^spei_", kendall_rp_all_scales$index_label) ~ "SPEI-PM",
      grepl("^swei_", kendall_rp_all_scales$index_label) ~ "SWEI",
      grepl("^sspi_", kendall_rp_all_scales$index_label) ~ "SSPI",
      grepl("^ssi_sm_", kendall_rp_all_scales$index_label) ~ "SSI-SM",
      grepl("^ssmi_", kendall_rp_all_scales$index_label) ~ "SSMI",
      grepl("^spi_", kendall_rp_all_scales$index_label) ~ "SPI",
      TRUE ~ "Other"
    )
    
    kendall_rp_all_scales$scale_months <- as.integer(
      sub(".*_([0-9]+)$", "\\1", kendall_rp_all_scales$index_label)
    )
    
    kendall_rp_all_scales <- kendall_rp_all_scales[
      order(
        kendall_rp_all_scales$index,
        kendall_rp_all_scales$scale_months,
        kendall_rp_all_scales$target_rank_recent
      ),
      ,
      drop = FALSE
    ]
    
    KENDALL_ALL_FILE <- file.path(
      RANKING_DIR,
      "Table_KendallRP_Recent_AllScales_AllEvents.csv"
    )
    KENDALL_PRIMARY_FILE <- file.path(
      RANKING_DIR,
      "Table_KendallRP_Recent_Primary_AllScales.csv"
    )
    KENDALL_SPEI_PM_FILE <- file.path(
      RANKING_DIR,
      paste0("Table_KendallRP_SPEI_PM_", RECENT_TAG_UNDERSCORE, "_AllScales.csv")
    )
    
    write.csv(kendall_rp_all_scales, KENDALL_ALL_FILE, row.names = FALSE)
    
    kendall_rp_primary <- kendall_rp_all_scales[
      kendall_rp_all_scales$target_rank_recent == 1L,
      ,
      drop = FALSE
    ]
    write.csv(kendall_rp_primary, KENDALL_PRIMARY_FILE, row.names = FALSE)
    
    kendall_rp_spei_pm <- kendall_rp_primary[
      kendall_rp_primary$index == "SPEI-PM",
      ,
      drop = FALSE
    ]
    write.csv(kendall_rp_spei_pm, KENDALL_SPEI_PM_FILE, row.names = FALSE)
    
    cat("\n  All-scale Kendall RP table: ", basename(KENDALL_ALL_FILE), "\n", sep = "")
    cat("  Primary-event table: ", basename(KENDALL_PRIMARY_FILE), "\n", sep = "")
    cat("  SPEI-PM practitioner table: ", basename(KENDALL_SPEI_PM_FILE), "\n", sep = "")
    
    if (nrow(kendall_rp_spei_pm)) {
      cat("\n  Primary SPEI-PM return periods:\n")
      print(
        kendall_rp_spei_pm[, c(
          "scale_months",
          "start_date",
          "end_date",
          "duration_months",
          "severity",
          "best_copula",
          "K_C",
          "Kendall_return_period_years",
          "CI_95_lo_years",
          "CI_95_hi_years"
        ), drop = FALSE],
        row.names = FALSE
      )
    }
  } else {
    cat("  WARNING: no all-scale Kendall RP analyses completed successfully.\n")
  }
}

# FINAL SUMMARY
####################################################################################
####################################################################################
# EXPORT BLOCK (part of 16event_ranking.R)
#
# PURPOSE: Export data structures for an external companion analysis
#          (a Bayesian analysis for a second manuscript; NOT part of this
#          framework). Runs after Parts 1–4g of this script.
#
# REQUIRES (already in environment from upstream script(s) main body):
#   ranked_catalog       — full SW08 event catalog with severity/duration
#   kendall_rp_result    — output of the Part 4c copula analysis (list)
#   n_years_record       — integer, length of the record in years
#   RANKING_DIR          — character, output directory path
#   HISTORICAL_START/END — shared, configurable record years
#
# OUTPUTS written to RANKING_DIR:
#   ms2_annual_severity_series.csv  — annual maximum severity per year (headline index)
#   ms2_copula_fit.rds              — copula fit object + stationary T_K reference
####################################################################################
cat("\n══════════════════════════════════════════════════\n")
cat("  EXPORT ADDON: MS2 data export\n")
cat("══════════════════════════════════════════════════\n")

# ── 1.  Annual maximum severity series (headline index only) ------------------
# For each calendar year, extract the maximum severity across all events
# whose start_date falls in that year. Years without events have NA severity,
# not zero: Gamma severity is conditional on an event (strictly positive support).
# A model for ALL years must include a separate occurrence component (hurdle/
# two-part model); do not treat NA as missing-at-random observations.

TARGET_IDX_EXPORT <- MS_INDEX

catalog_spei3 <- ranked_catalog %>%
  dplyr::filter(index_label == TARGET_IDX_EXPORT) %>%
  dplyr::mutate(
    start_year_int = as.integer(format(as.Date(start_date), "%Y")),
    end_year_int   = as.integer(format(as.Date(end_date),   "%Y"))
  )

# Build complete year spine
year_spine <- data.frame(year = HISTORICAL_START:HISTORICAL_END)

# Annual max severity uses event onset year, not year of peak severity.
# D1-2 events remain included; downstream manuscript-2 models must explicitly
# define whether this selection/anchoring convention matches the estimand.
ann_sev <- catalog_spei3 %>%
  dplyr::group_by(year = start_year_int) %>%
  dplyr::summarise(
    max_severity   = max(severity,        na.rm = TRUE),
    max_duration   = max(duration_months, na.rm = TRUE),
    n_events       = dplyr::n(),
    .groups        = "drop"
  )

ann_sev_full <- dplyr::left_join(year_spine, ann_sev, by = "year") %>%
  dplyr::mutate(
    drought_year = !is.na(max_severity),
    no_event_year = is.na(max_severity),
    index        = TARGET_IDX_EXPORT
  )

ann_sev_path <- file.path(RANKING_DIR, "ms2_annual_severity_series.csv")
write.csv(ann_sev_full, ann_sev_path, row.names = FALSE)
cat(sprintf("  ✓ Annual severity series: %s\n", basename(ann_sev_path)))
cat(sprintf("    %d drought years out of %d total years\n",
            sum(ann_sev_full$drought_year, na.rm = TRUE),
            nrow(ann_sev_full)))

# ── 2.  Copula fit object + stationary T_K reference --------------------------
# If the Part 4c Kendall RP analysis ran successfully, export:
#   - best_fit:       the fitted copula object (for plotting level curves in MS2)
#   - best_dur/sev:   the marginal distribution objects (for computing u0/v0)
#   - T_K_stationary: the stationary Kendall return period (reference for MS2)
#   - target_event:   severity and duration of the configured recent window rank-1 headline-index event

if (exists("kendall_rp_result") && !is.null(kendall_rp_result)) {
  
  krp <- kendall_rp_result
  if (!isTRUE(all(krp$rows$duration_shift_months == 0)))
    stop("MS2 export requires the unshifted (duration_shift = 0) Kendall fit")
  
  # Extract the target event row from the catalog
  target_event_row <- catalog_spei3 %>%
    dplyr::filter(is_recent, duration_class != "D1-2 (Very short)",
                  duration_months >= 3L, severity > 0) %>%
    dplyr::arrange(dplyr::desc(severity), start_date, end_date) %>%
    dplyr::slice_head(n = 1L)
  if (nrow(target_event_row) != 1L ||
      !isTRUE(all.equal(as.numeric(target_event_row$severity[1]), as.numeric(krp$tgt_sev))) ||
      !identical(as.integer(target_event_row$duration_months[1]), as.integer(krp$tgt_dur)))
    stop("MS2 export target does not agree with the Part 4c Kendall target")
  
  ms2_copula_export <- list(
    # Stationary copula objects
    best_copula_fit   = krp$best_fit,       # copulaFit object
    best_copula_name  = krp$best_cop,       # character, e.g. "Frank"
    best_dur_fit      = krp$fits_dur[[krp$best_dur]],  # fitdist object; unshifted integer durations
    best_sev_fit      = krp$fits_sev[[krp$best_sev]],  # fitdistr object
    best_dur_name     = krp$best_dur,
    duration_shift_months = 0L,  # explicit MS2 compatibility contract
    best_sev_name     = krp$best_sev,
    # Stationary T_K (the reference value for MS2 comparison)
    T_K_stationary    = krp$T_kendall,
    CI_lo_stationary  = krp$ci_lo,
    CI_hi_stationary  = krp$ci_hi,
    # Target event
    target_severity   = target_event_row$severity[1],
    target_duration   = target_event_row$duration_months[1],
    target_period     = sprintf("%s to %s",
                                target_event_row$start_date[1],
                                target_event_row$end_date[1]),
    # Probability-integral-transformed coordinates of target event
    u0                = krp$u0,
    v0                = krp$v0,
    t_star            = krp$t_star,
    # Inter-arrival time used in T_K formula
    mu_ia_years       = krp$mu_ia,
    # Record metadata
    n_events          = krp$n_ev,
    n_years_record    = n_years_record,
    index_label       = TARGET_IDX_EXPORT,
    # Sensitivity table across all copula families (Supplementary Table S1)
    sens_df           = krp$sens_df
  )
  
  rds_path <- file.path(RANKING_DIR, "ms2_copula_fit.rds")
  saveRDS(ms2_copula_export, rds_path)
  cat(sprintf("  ✓ Copula fit RDS: %s\n", basename(rds_path)))
  cat(sprintf("    Stationary T_K = %.0f years (95%% CI: %.0f–%.0f)\n",
              ms2_copula_export$T_K_stationary,
              ms2_copula_export$CI_lo_stationary,
              ms2_copula_export$CI_hi_stationary))
  cat(sprintf("    Target event:  S = %.3f, D = %d months\n",
              ms2_copula_export$target_severity,
              ms2_copula_export$target_duration))
  
} else {
  cat("  ⚠ kendall_rp_result not found — ms2_copula_fit.rds not written.\n")
  cat("    Run Part 4c of 16event_ranking.R first.\n")
}
################################################################################
# SENSITIVITY: duration marginal fitted to (D - shift)
# Calls run_kendall_rp() directly. It never goes through get_kendall_rp(), so
# the cache and the MS2 export (shift = 0) are untouched.
# The shift = 0 row must reproduce the headline result exactly (same seed,
# same n_boot). If it does not, one of the four duration_shift edits is wrong.
################################################################################
cat("\n══════════════════════════════════════════════════\n")
cat("  SENSITIVITY: duration-shift (discrete / left-truncated marginal)\n")
cat("══════════════════════════════════════════════════\n")

SENS_SHIFT_INDICES <- unique(c(MS_INDEX))   # e.g. c(MS_INDEX, "spi_03", "spei_06")
SENS_SHIFT_INDICES <- SENS_SHIFT_INDICES[SENS_SHIFT_INDICES %in% unique(master_events$index_label)]

sens_shift_list <- list()
for (idx_s in SENS_SHIFT_INDICES) {
  sc_s  <- as.integer(sub(".*_([0-9]+)$", "\\1", idx_s))
  d_min <- max(3L, get_event_min_duration(sc_s))   # D1-2 events are never in the fit
  for (s in c(0, d_min - 1, d_min - 0.5)) {
    cat(sprintf("\n  --- %s | duration_shift = %.1f ---\n", idx_s, s))
    r <- run_kendall_rp(idx_s, n_boot = KENDALL_BOOT_DEFAULT, duration_shift = s)
    if (!is.null(r$error)) {
      cat(sprintf("  skipped: %s\n", r$error)); next
    }
    sens_shift_list[[length(sens_shift_list) + 1L]] <- data.frame(
      index_label  = idx_s,
      shift_months = s,
      dur_marginal = r$best_dur,
      sev_marginal = r$best_sev,
      copula       = r$best_cop,
      tau_emp      = r$tau_emp,
      u0           = round(r$u0, 4),
      v0           = round(r$v0, 4),
      K_C          = round(r$Kc, 5),
      TK_yr        = round(r$TK, 1),
      CI_lo_yr     = round(r$ci_lo, 1),
      CI_hi_yr     = round(r$ci_hi, 1),
      stringsAsFactors = FALSE
    )
  }
}

if (length(sens_shift_list)) {
  sens_shift <- dplyr::bind_rows(sens_shift_list) %>%
    dplyr::group_by(index_label) %>%
    dplyr::mutate(
      TK_ratio_vs_shift0 = round(TK_yr / TK_yr[shift_months == 0][1], 2),
      inside_shift0_CI   = TK_yr >= CI_lo_yr[shift_months == 0][1] &
        TK_yr <= CI_hi_yr[shift_months == 0][1]
    ) %>%
    dplyr::ungroup()
  write.csv(sens_shift,
            file.path(RANKING_DIR, "Table_KendallRP_duration_shift_sensitivity.csv"),
            row.names = FALSE)
  cat("\n  Duration-shift sensitivity:\n")
  print(as.data.frame(sens_shift), row.names = FALSE)
  cat("  Saved: Table_KendallRP_duration_shift_sensitivity.csv\n")
} else {
  cat("  ⚠ No duration-shift sensitivity results produced.\n")
}
################################################################################
# REGIONAL EXTENSION — SW08 EVENT RANKING FOR 7 SUBBASINS + WHOLE BASIN
################################################################################
REGIONAL_EVENT_DIR <- file.path(TREND_DIR, "regional_event_ranking")
REGIONAL_INDEX_DIR <- file.path(TREND_DIR, "regional_timeseries")
dir.create(REGIONAL_EVENT_DIR, recursive=TRUE, showWarnings=FALSE)

# Read only top-level non-soil series and soil series from the selected layer.
regional_event_files <- c(
  list.files(REGIONAL_INDEX_DIR, pattern="_regional_averaged_by_month\\.csv$",
             full.names=TRUE, recursive=FALSE, ignore.case=TRUE),
  list.files(file.path(REGIONAL_INDEX_DIR, SOIL_LAYER),
             pattern="_regional_averaged_by_month\\.csv$",
             full.names=TRUE, recursive=FALSE, ignore.case=TRUE)
)
regional_event_catalogs <- list()

for (f in regional_event_files) {
  x <- data.table::fread(f)
  req <- c("region_id","region_name","date","value","index_type","scale")
  if (!all(req %in% names(x))) next
  idx <- unique(tolower(x$index_type))[1]
  if (!idx %in% c("spi", "spei", "swei", "sspi", "ssmi", "ssi_sm")) next
  if ((idx %in% c("ssmi", "ssi_sm")) != (basename(dirname(f)) == SOIL_LAYER)) next
  scs <- sort(unique(as.integer(x$scale)))
  scs <- intersect(scs, as.integer(names(EVENT_MIN_DURATION_BY_SCALE)))
  if (idx %in% c("ssmi", "ssi_sm")) scs <- intersect(scs, SOIL_SCALES)
  for (rid in unique(x$region_id)) {
    for (sc in scs) {
      z <- x[x$region_id==rid & as.integer(x$scale)==sc]
      if (!nrow(z)) next
      z <- z[order(date)]
      dates_z <- as.Date(z$date)
      vals_z <- as.numeric(z$value)
      keep <- dates_z >= as.Date(sprintf("%04d-01-01",HISTORICAL_START)) &
        dates_z <= as.Date(sprintf("%04d-12-01",HISTORICAL_END))
      dates_z <- dates_z[keep]; vals_z <- vals_z[keep]
      if (!length(vals_z)) next
      ev <- tryCatch(
        sw08_event_loop(vals_z,dates_z,as.integer(format(dates_z,"%Y")),
                        onset_thr=DROUGHT_ONSET_EVENT,
                        min_duration=get_event_min_duration(sc)),
        error = function(e) {
          message(sprintf("  WARNING: skipped regional SW08 series %s / %s-%02d: %s",
                          rid, idx, sc, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(ev) || !nrow(ev)) next
      ev$region_id <- rid
      ev$region_name <- unique(z$region_name)[1]
      ev$index_type <- idx
      ev$scale <- as.integer(sc)
      ev$index_label <- sprintf("%s_%02d",idx,sc)
      regional_event_catalogs[[paste(rid,idx,sc,sep="|")]] <- ev
    }
  }
}

if (length(regional_event_catalogs)) {
  regional_events <- dplyr::bind_rows(regional_event_catalogs)
  regional_events <- regional_events |>
    dplyr::group_by(region_id,region_name,index_label) |>
    dplyr::arrange(dplyr::desc(severity),start_date,.by_group=TRUE) |>
    dplyr::mutate(
      rank_by_severity=dplyr::min_rank(dplyr::desc(severity)),
      return_period_years=(n_years_record+1)/rank_by_severity,
      percentile_rank=round(100*(1-rank_by_severity/(dplyr::n()+1)),1),
      is_recent=start_date <= RECENT_DROUGHT_END & end_date >= RECENT_DROUGHT_START
    ) |>
    dplyr::ungroup()
  
  regional_rp <- regional_events |>
    dplyr::filter(duration_class != "D1-2 (Very short)") |>
    dplyr::group_by(region_id,region_name,index_label,duration_class) |>
    dplyr::summarise(
      n_events=dplyr::n(),
      mean_duration_mo=mean(duration_months,na.rm=TRUE),
      mean_I=mean(mean_int,na.rm=TRUE),
      mean_S=mean(severity,na.rm=TRUE),
      max_S=max(severity,na.rm=TRUE),
      return_period_yrs=n_years_record/n_events,
      .groups="drop"
    )
  regional_primary <- regional_events |>
    dplyr::filter(is_recent) |>
    dplyr::group_by(region_id,region_name,index_label) |>
    dplyr::slice_max(severity,n=1,with_ties=FALSE) |>
    dplyr::ungroup()
  
  write.csv(regional_events,file.path(REGIONAL_EVENT_DIR,"ranked_event_catalog_SW08_regional.csv"),row.names=FALSE)
  write.csv(regional_rp,file.path(REGIONAL_EVENT_DIR,"return_period_by_duration_class_SW08_regional.csv"),row.names=FALSE)
  write.csv(regional_primary,file.path(REGIONAL_EVENT_DIR,"recent_primary_events_SW08_regional.csv"),row.names=FALSE)
  
  p <- ggplot2::ggplot(regional_events,ggplot2::aes(duration_months,severity)) +
    ggplot2::geom_point(alpha=.45,size=1) + ggplot2::scale_x_continuous(trans="sqrt") +
    ggplot2::facet_grid(index_label~region_name,scales="free") +
    ggplot2::labs(title="Regional SW08 drought-event duration–severity distribution",
                  x="Duration (months)",y="Severity") + ggplot2::theme_bw(base_size=7)
  ggplot2::ggsave(file.path(REGIONAL_EVENT_DIR,"Fig_SW08_duration_severity_regional.png"),p,width=18,height=12,dpi=300)
  ggplot2::ggsave(file.path(REGIONAL_EVENT_DIR,"Fig_SW08_duration_severity_regional.pdf"),p,width=18,height=12,
                  device=if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf)
  cat(sprintf("✓ Regional SW08 catalog written for %d regions.\n",uniqueN(regional_events$region_id)))
}

################################################################################
# PART 4h: REGIONAL HISTORICAL-RANKING + KENDALL-RP FIGURES
#
# Produces 8 separate figures:
#   1 whole Nechako basin + 7 supplied subbasins.
#
# Figure A (one per region): recent-event historical ranking profile
#   • SPI and SPEI-PM at manuscript scales; SWEI-03, SSPI-01 and
#     the soil-moisture indices at native scales
#   • available accumulation scales: 1, 3, 6, 12, 24, 36 months
#   • horizontal segment = configured recent window overlap of the selected recent event
#   • bubble area = event severity
#   • colour = historical severity-rank class within region × index × scale
#   • label = historical rank and severity
#
# Figure B (one per region): Kendall duration–severity return-period profile
#   • SPI, SPEI-PM, SWEI, SSPI, SSMI and SSI-SM where regional event data exist
#   • available accumulation scales; snow and soil indices use native scales only
#   • point = Kendall joint return period of the primary recent event
#   • horizontal interval = 95% bootstrap CI
#   • x-axis = years on a logarithmic scale
#   • labels report T_K, event duration, and severity
#
# The regional event catalog is the same SW08 catalog generated above from the
# authoritative regional monthly time series. No regional spatial averaging is
# recomputed here. The Kendall calculation reuses the canonical run_kendall_rp()
# engine with a supplied regional event table.
#
# IMPORTANT CURRENT-FRAMEWORK NOTE:
#   Thornthwaite SPEI is not a current observed formulation in this project.
#   Legacy Thornthwaite-labelled products are not used as current evidence.
################################################################################
cat("\n======================================================================\n")
cat("PART 4h: Regional historical-ranking and Kendall-RP figures\n")
cat("======================================================================\n")

# Regional Monte Carlo budget accommodates soil-moisture families.
# Nominal workload per recent target, if all series are available:
# 8 regions x 18 series x 250 = 36,000 bootstrap replicates,
# versus the old 8 x 14 x 500 = 56,000. Multiple recent targets multiply these
# counts; absent series and failed fits reduce them. Regional CI Monte Carlo
# precision is lower than at 500 draws. Basin analyses retain 500 replicates.
REGIONAL_KENDALL_BOOT <- 250L
# Whole-basin regional series == producer basin series (same weights), so its
# Kendall fit is already cached from Part 4g. Set FALSE to refit it here.
REGIONAL_REUSE_BASIN_KENDALL <- TRUE
regional_matches_basin <- function(ev_reg, idx_label) {
  b <- as.data.frame(master_events)[master_events$index_label == idx_label,
                                    c("start_date", "end_date", "severity")]
  r <- as.data.frame(ev_reg)[, c("start_date", "end_date", "severity")]
  if (nrow(b) != nrow(r) || !nrow(b)) return(FALSE)
  b <- b[order(b$start_date), ]; r <- r[order(r$start_date), ]
  isTRUE(all.equal(as.numeric(as.Date(b$start_date)), as.numeric(as.Date(r$start_date)))) &&
    isTRUE(all.equal(as.numeric(as.Date(b$end_date)), as.numeric(as.Date(r$end_date)))) &&
    isTRUE(all.equal(as.numeric(b$severity), as.numeric(r$severity), tolerance = 1e-8))
}
REGIONAL_FIG_START <- RECENT_DROUGHT_START
REGIONAL_FIG_END   <- RECENT_DROUGHT_END
REGIONAL_FIG_SCALES <- as.integer(MANUSCRIPT_SCALES)
# Regional rankings and Kendall RP cover the configured producer indices at
# their eligible scales (snow and soil-moisture indices at native scales only).
REGIONAL_FIG_INDEXES <- c("spi", "spei", "swei", "sspi", "ssmi", "ssi_sm")
REGIONAL_KENDALL_INDEXES <- c("spi", "spei", "swei", "sspi",
                              if (length(SOIL_SCALES)) c("ssmi", "ssi_sm"))
REGIONAL_RP_FIG_WIDTH <- 2.3 * length(REGIONAL_KENDALL_INDEXES) + 2
REGIONAL_FIG_LABELS <- c(
  spi       = "(a) SPI",
  spei      = "(b) SPEI-PM",
  swei      = "(c) SWEI",
  sspi      = "(d) SSPI",
  ssmi      = "(e) SSMI",
  ssi_sm    = "(f) SSI-SM"
)

REGIONAL_RANK_COLOURS <- c(
  "Rank 1"   = "#C91F2D",
  "Rank 2"   = "#F47C20",
  "Rank 3-5" = "#F5B82E",
  "Rank >5"  = "#54A9B2",
  "Censored / unranked" = "#777777"
)

if (!exists("regional_events") || !nrow(regional_events)) {
  cat("  WARNING: regional_events is empty; Part 4h skipped.\n")
} else {
  # --------------------------------------------------------------------------
  # 4h.1 — restrict to configured observed formulations and manuscript scales;
  #        snow and soil indices appear only at their valid native scales.
  # --------------------------------------------------------------------------
  regional_plot_events <- regional_events %>%
    dplyr::filter(
      tolower(index_type) %in% REGIONAL_FIG_INDEXES,
      as.integer(scale) %in% REGIONAL_FIG_SCALES
    ) %>%
    dplyr::mutate(
      index_type = tolower(index_type),
      scale = as.integer(scale)
    )
  
  if (!nrow(regional_plot_events)) {
    cat("  WARNING: no configured regional events available; Part 4h skipped.\n")
  } else {
    regional_ids <- unique(as.character(regional_plot_events$region_id))
    regional_names <- regional_plot_events %>%
      dplyr::distinct(region_id, region_name) %>%
      dplyr::mutate(
        region_id = as.character(region_id),
        region_name = as.character(region_name)
      ) %>%
      dplyr::arrange(match(region_id, regional_ids))
    
    # Prefer the authoritative regional catalog order if available.
    region_catalog_file_h <- file.path(TREND_DIR, "regional_timeseries", "region_catalog.csv")
    if (file.exists(region_catalog_file_h)) {
      rc_h <- tryCatch(
        data.table::fread(region_catalog_file_h) %>% as.data.frame(),
        error = function(e) NULL
      )
      if (!is.null(rc_h) && all(c("region_id", "region_name") %in% names(rc_h))) {
        rc_ids_h <- unique(as.character(rc_h$region_id))
        regional_ids <- c(
          rc_ids_h[rc_ids_h %in% regional_ids],
          setdiff(regional_ids, rc_ids_h)
        )
        regional_names <- regional_names %>%
          dplyr::arrange(match(as.character(region_id), regional_ids))
      }
    }
    
    # Safety check requested by the regional framework: one whole basin + seven
    # supplied subbasins.
    if (length(regional_ids) != 8L) {
      cat(sprintf(
        "  WARNING: expected 8 regions (whole basin + 7 subbasins), found %d.\n",
        length(regional_ids)
      ))
    }
    
    regional_rank_rows <- list()
    regional_kendall_rows <- list()
    
    safe_region_filename <- function(x) {
      y <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
      y <- gsub("_+", "_", y)
      y <- gsub("^_|_$", "", y)
      if (!nzchar(y)) y <- "region"
      y
    }
    
    pretty_regional_index <- function(idx) {
      unname(REGIONAL_FIG_LABELS[idx])
    }
    
    # ------------------------------------------------------------------------
    # 4h.2 — process each region independently.
    # ------------------------------------------------------------------------
    for (r_i in seq_along(regional_ids)) {
      rid <- regional_ids[r_i]
      rname <- regional_names$region_name[match(rid, regional_names$region_id)]
      if (is.na(rname) || !nzchar(rname)) rname <- rid
      
      reg <- regional_plot_events %>%
        dplyr::filter(as.character(region_id) == rid)
      
      if (!nrow(reg)) next
      
      # ======================================================================
      # FIGURE A — historical ranking of the configured recent window primary event
      # ======================================================================
      hist_recent <- reg %>%
        dplyr::filter(
          is_recent,
          as.Date(end_date) >= REGIONAL_FIG_START,
          as.Date(start_date) <= REGIONAL_FIG_END
        ) %>%
        dplyr::arrange(index_type, scale, dplyr::desc(severity), start_date, end_date) %>%
        dplyr::group_by(index_type, scale) %>%
        dplyr::slice_head(n = 1L) %>%
        dplyr::ungroup()
      
      if (nrow(hist_recent)) {
        hist_recent <- hist_recent %>%
          dplyr::mutate(
            index_panel = factor(
              unname(REGIONAL_FIG_LABELS[index_type]),
              levels = unname(REGIONAL_FIG_LABELS[REGIONAL_FIG_INDEXES])
            ),
            scale_y = match(scale, rev(REGIONAL_FIG_SCALES)),
            display_start = pmax(as.Date(start_date), REGIONAL_FIG_START),
            display_end = pmin(as.Date(end_date), REGIONAL_FIG_END),
            display_mid = display_start + floor(
              as.numeric(display_end - display_start) / 2
            ),
            rank_class = dplyr::case_when(
              rank_by_severity == 1L ~ "Rank 1",
              rank_by_severity == 2L ~ "Rank 2",
              rank_by_severity >= 3L & rank_by_severity <= 5L ~ "Rank 3-5",
              TRUE ~ "Rank >5"
            ),
            rank_label = sprintf("Rank %d", as.integer(rank_by_severity)),
            severity_label = sprintf("%.2f", as.numeric(severity))
          )
        
        # Keep a machine-readable table of the exact rows plotted.
        regional_rank_rows[[length(regional_rank_rows) + 1L]] <- hist_recent %>%
          dplyr::transmute(
            region_id = as.character(region_id),
            region_name = as.character(region_name),
            index_type,
            index_label,
            scale_months = as.integer(scale),
            start_date = as.Date(start_date),
            end_date = as.Date(end_date),
            duration_months = as.integer(duration_months),
            severity = as.numeric(severity),
            rank_by_severity = as.integer(rank_by_severity),
            percentile_rank = as.numeric(percentile_rank),
            return_period_years = as.numeric(return_period_years)
          )
        
        p_rank <- ggplot2::ggplot(hist_recent) +
          ggplot2::geom_segment(
            ggplot2::aes(
              x = display_start, xend = display_end,
              y = scale_y, yend = scale_y
            ),
            linewidth = 0.9, colour = "#303030"
          ) +
          ggplot2::geom_segment(
            ggplot2::aes(
              x = display_start, xend = display_start,
              y = scale_y - 0.11,
              yend = scale_y + 0.11
            ),
            linewidth = 0.55, colour = "#303030"
          ) +
          ggplot2::geom_segment(
            ggplot2::aes(
              x = display_end, xend = display_end,
              y = scale_y - 0.11,
              yend = scale_y + 0.11
            ),
            linewidth = 0.55, colour = "#303030"
          ) +
          ggplot2::geom_point(
            ggplot2::aes(
              x = display_mid,
              y = scale_y,
              size = severity,
              fill = rank_class
            ),
            shape = 21, colour = "white", stroke = 0.8
          ) +
          ggplot2::geom_text(
            ggplot2::aes(
              x = display_mid,
              y = scale_y,
              label = severity_label
            ),
            vjust = -1.35, size = 3.0, fontface = "bold", colour = "#202020"
          ) +
          ggplot2::facet_wrap(
            ~ index_panel, nrow = 1, scales = "fixed"
          ) +
          ggplot2::scale_fill_manual(
            values = REGIONAL_RANK_COLOURS,
            drop = FALSE
          ) +
          ggplot2::scale_size_area(
            max_size = 16,
            guide = ggplot2::guide_legend(
              title = "Severity",
              order = 1
            )
          ) +
          ggplot2::scale_y_continuous(
            breaks = seq_along(rev(REGIONAL_FIG_SCALES)),
            labels = rev(REGIONAL_FIG_SCALES),
            expand = ggplot2::expansion(mult = c(0.06, 0.06))
          ) +
          ggplot2::scale_x_date(
            limits = c(REGIONAL_FIG_START, REGIONAL_FIG_END),
            date_breaks = "1 year",
            date_labels = "%Y",
            expand = ggplot2::expansion(mult = c(0.035, 0.035))
          ) +
          ggplot2::labs(
            title = paste0(rname, " — recent drought events in historical ranking"),
            subtitle = paste0(
              paste0(RECENT_LABEL, " primary event at each accumulation scale; bubble size = severity; "),
              "colour = historical severity rank"
            ),
            x = NULL,
            y = "Accumulation scale (months)",
            fill = "Historical rank"
          ) +
          theme_ms3 +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 12, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 9, colour = "grey35"),
            strip.text = ggplot2::element_text(size = 10, face = "bold"),
            panel.grid.major.y = ggplot2::element_line(
              colour = "grey88", linewidth = 0.35, linetype = "dashed"
            ),
            panel.grid.major.x = ggplot2::element_line(
              colour = "grey92", linewidth = 0.3
            ),
            legend.position = "bottom"
          )
        
        rank_file_stub <- paste0(
          "Fig_Regional_Historical_Ranking_",
          safe_region_filename(rname)
        )
        tryCatch({
          ggplot2::ggsave(
            file.path(REGIONAL_EVENT_DIR, paste0(rank_file_stub, ".png")),
            p_rank, width = 13, height = 5.8, units = "in", dpi = 320
          )
          ggplot2::ggsave(
            file.path(REGIONAL_EVENT_DIR, paste0(rank_file_stub, ".pdf")),
            p_rank, width = 13, height = 5.8, units = "in", device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
          )
          cat(sprintf("  Saved historical-ranking figure: %s\n", rname))
        }, error = function(e) {
          cat(sprintf("  WARNING historical-ranking figure [%s]: %s\n", rname, e$message))
        })
      } else {
        cat(sprintf("  WARNING: no recent primary events for region %s\n", rname))
      }
      
      # ======================================================================
      # FIGURE B — Kendall joint return periods for the same recent primary
      # events. The calculation is performed separately for every region ×
      # index × scale using the canonical Kendall engine.
      # ======================================================================
      reg_kendall <- list()
      
      for (j in seq_along(REGIONAL_KENDALL_INDEXES)) {
        idx_type_j <- REGIONAL_KENDALL_INDEXES[j]
        idx_labels_j <- unique(reg$index_label[reg$index_type == idx_type_j])
        idx_labels_j <- idx_labels_j[grepl("_[0-9]+$", idx_labels_j)]
        idx_labels_j <- idx_labels_j[
          as.integer(sub(".*_([0-9]+)$", "\\1", idx_labels_j)) %in% REGIONAL_FIG_SCALES
        ]
        
        for (k in seq_along(idx_labels_j)) {
          idx_j <- idx_labels_j[k]
          sc_j <- as.integer(sub(".*_([0-9]+)$", "\\1", idx_j))
          ev_j <- reg %>%
            dplyr::filter(
              index_label == idx_j,
              is.finite(duration_months),
              is.finite(severity),
              severity > 0
            )
          
          if (!nrow(ev_j)) next
          
          # A deterministic region × index × scale seed makes reruns reproducible.
          seed_j <- 42000L + r_i * 1000L + j * 100L + k
          kr_j <- if (isTRUE(REGIONAL_REUSE_BASIN_KENDALL) &&
                      identical(as.character(rid), NECHAKO_REGION_BASIN_ID) &&
                      regional_matches_basin(ev_j, idx_j)) {
            cat(sprintf("    %s / %s: events identical to basin catalog; reusing basin Kendall fit\n",
                        rname, idx_j))
            get_kendall_rp(idx_j, n_boot = KENDALL_BOOT_DEFAULT)
          } else run_kendall_rp(
            idx_j,
            n_boot = REGIONAL_KENDALL_BOOT,
            seed = seed_j,
            events_data = ev_j,
            record_years = n_years_record
          )
          
          if (is.list(kr_j) && is.null(kr_j$error) &&
              is.data.frame(kr_j$rows) && nrow(kr_j$rows)) {
            rr <- kr_j$rows %>%
              dplyr::slice_head(n = 1L) %>%
              dplyr::mutate(
                region_id = rid,
                region_name = rname,
                index_type = idx_type_j,
                scale_months = sc_j,
                .before = 1
              )
            reg_kendall[[length(reg_kendall) + 1L]] <- rr
          } else {
            # Retain a row for failed/insufficient fits so the output table is
            # explicit rather than silently dropping a scale.
            reg_kendall[[length(reg_kendall) + 1L]] <- data.frame(
              region_id = rid,
              region_name = rname,
              index_type = idx_type_j,
              index_label = idx_j,
              scale_months = sc_j,
              start_date = as.Date(NA),
              end_date = as.Date(NA),
              duration_months = NA_integer_,
              severity = NA_real_,
              n_historical_events = nrow(ev_j),
              mean_interarrival_years = ifelse(nrow(ev_j) > 0, n_years_record / nrow(ev_j), NA_real_),
              tau_empirical = NA_real_,
              tau_fitted = NA_real_,
              duration_marginal = NA_character_,
              severity_marginal = NA_character_,
              best_copula = NA_character_,
              copula_theta = NA_real_,
              u_duration = NA_real_,
              v_severity = NA_real_,
              t_star = NA_real_,
              K_C = NA_real_,
              Kendall_return_period_years = NA_real_,
              CI_95_lo_years = NA_real_,
              CI_95_hi_years = NA_real_,
              bootstrap_median_years = NA_real_,
              bootstrap_converged = NA_integer_,
              region_error = if (is.list(kr_j) && !is.null(kr_j$error)) kr_j$error else "No valid Kendall result",
              stringsAsFactors = FALSE
            )
          }
        }
      }
      
      if (length(reg_kendall)) {
        reg_kendall_df <- dplyr::bind_rows(reg_kendall) %>%
          dplyr::mutate(index_panel_lab = factor(as.character(index_type),
                                                 levels = REGIONAL_KENDALL_INDEXES,
                                                 labels = unname(REGIONAL_FIG_LABELS[REGIONAL_KENDALL_INDEXES])))
        
        regional_kendall_rows[[length(regional_kendall_rows) + 1L]] <- reg_kendall_df
        
        rp_plot_df <- reg_kendall_df %>%
          dplyr::mutate(
            index_panel = index_panel_lab,
            TK = as.numeric(Kendall_return_period_years),
            ci_lo = as.numeric(CI_95_lo_years),
            ci_hi = as.numeric(CI_95_hi_years),
            d_lab = ifelse(is.finite(duration_months),
                           paste0("D=", duration_months, " mo"), "D unavailable"),
            s_lab = ifelse(is.finite(severity),
                           paste0("S=", sprintf("%.2f", severity)), "S unavailable"),
            tk_lab = ifelse(is.finite(TK), sprintf("%.0f yr", TK), "Not estimable")
          )
        
        rp_valid <- rp_plot_df %>% dplyr::filter(is.finite(TK), TK > 0)
        rp_invalid <- rp_plot_df %>% dplyr::filter(!is.finite(TK) | TK <= 0)
        
        p_rp <- ggplot2::ggplot(rp_valid) +
          ggplot2::geom_segment(
            ggplot2::aes(
              x = pmax(ci_lo, 0.5),
              xend = pmax(ci_hi, 0.5),
              y = factor(scale_months, levels = rev(REGIONAL_FIG_SCALES)),
              yend = factor(scale_months, levels = rev(REGIONAL_FIG_SCALES))
            ),
            linewidth = 0.7, colour = "#3A3A3A", na.rm = TRUE
          ) +
          ggplot2::geom_point(
            ggplot2::aes(
              x = TK,
              y = factor(scale_months, levels = rev(REGIONAL_FIG_SCALES))
            ),
            shape = 21, size = 4.0, fill = "#C91F2D", colour = "white", stroke = 0.7
          ) +
          ggplot2::geom_text(
            ggplot2::aes(
              x = TK,
              y = factor(scale_months, levels = rev(REGIONAL_FIG_SCALES)),
              label = tk_lab
            ),
            hjust = -0.10, size = 2.8, fontface = "bold", colour = "#303030"
          ) +
          ggplot2::facet_wrap(
            ~ index_panel, nrow = 1, scales = "free_y"
          ) +
          ggplot2::scale_x_log10(
            breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000),
            labels = scales::comma,
            name = "Kendall joint return period (years; log scale)",
            expand = ggplot2::expansion(mult = c(0.04, 0.32))
          ) +
          ggplot2::scale_y_discrete(drop = FALSE) +
          ggplot2::labs(
            title = paste0(rname, " — Kendall duration–severity return periods"),
            subtitle = paste0(
              paste0("Primary ", RECENT_LABEL, " event at each scale; point = T_K; interval = 95% bootstrap CI; "),
              "D = duration, S = severity"
            ),
            y = "Accumulation scale (months)"
          ) +
          theme_ms3 +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 12, face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 9, colour = "grey35"),
            strip.text = ggplot2::element_text(size = 10, face = "bold"),
            panel.grid.major.y = ggplot2::element_line(
              colour = "grey88", linewidth = 0.35, linetype = "dashed"
            ),
            panel.grid.major.x = ggplot2::element_line(
              colour = "grey92", linewidth = 0.3
            )
          )
        
        # Add D/S diagnostic text in a separate layer only where estimates exist.
        if (nrow(rp_valid)) {
          p_rp <- p_rp + ggplot2::geom_text(
            data = rp_valid,
            ggplot2::aes(
              x = TK,
              y = factor(scale_months, levels = rev(REGIONAL_FIG_SCALES)),
              label = paste(d_lab, s_lab, sep = "  |  ")
            ),
            hjust = -0.05, vjust = 1.8, size = 2.25, colour = "grey40"
          )
        }
        
        # If a scale has no valid Kendall fit (commonly D1–2 at scale 1),
        # show it explicitly rather than implying that the value is zero.
        if (nrow(rp_invalid)) {
          # Put invalid rows at the left side of the panel in a fixed, finite
          # location on the log axis. This is annotation only, not a numeric RP.
          invalid_annot <- rp_invalid %>%
            dplyr::mutate(x_annot = 1.15)
          p_rp <- p_rp + ggplot2::geom_text(
            data = invalid_annot,
            ggplot2::aes(
              x = x_annot,
              y = factor(scale_months, levels = rev(REGIONAL_FIG_SCALES)),
              label = "Not estimable"
            ),
            hjust = 0, size = 2.5, colour = "grey45", fontface = "italic"
          )
        }
        
        # Use the actual fitted copula information in the caption where available.
        if (any(is.finite(reg_kendall_df$Kendall_return_period_years))) {
          cop_txt <- paste0(
            "Copula/tau diagnostics are retained in the regional Kendall table; ",
            "each point is fitted independently for its region × index × scale."
          )
        } else {
          cop_txt <- "No valid regional Kendall fit was obtained for the plotted scales."
        }
        p_rp <- p_rp + patchwork::plot_annotation(
          caption = cop_txt,
          theme = ggplot2::theme(
            plot.caption = ggplot2::element_text(
              size = 7.5, colour = "grey45", hjust = 0,
              margin = ggplot2::margin(t = 5)
            )
          )
        )
        
        rp_file_stub <- paste0(
          "Fig_Regional_KendallRP_",
          safe_region_filename(rname)
        )
        tryCatch({
          ggplot2::ggsave(
            file.path(REGIONAL_EVENT_DIR, paste0(rp_file_stub, ".png")),
            p_rp, width = REGIONAL_RP_FIG_WIDTH, height = 6.4, units = "in", dpi = 320
          )
          ggplot2::ggsave(
            file.path(REGIONAL_EVENT_DIR, paste0(rp_file_stub, ".pdf")),
            p_rp, width = REGIONAL_RP_FIG_WIDTH, height = 6.4, units = "in",
            device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
          )
          cat(sprintf("  Saved Kendall-RP figure: %s\n", rname))
        }, error = function(e) {
          cat(sprintf("  WARNING Kendall-RP figure [%s]: %s\n", rname, e$message))
        })
      }
    }
    
    # ------------------------------------------------------------------------
    # 4h.3 — consolidated machine-readable outputs.
    # ------------------------------------------------------------------------
    if (length(regional_rank_rows)) {
      regional_rank_table <- dplyr::bind_rows(regional_rank_rows) %>%
        dplyr::arrange(region_name, index_type, scale_months)
      write.csv(
        regional_rank_table,
        file.path(REGIONAL_EVENT_DIR, paste0("Table_Regional_Historical_Ranking_", RECENT_TAG_UNDERSCORE, ".csv")),
        row.names = FALSE
      )
    }
    
    if (length(regional_kendall_rows)) {
      regional_kendall_table <- dplyr::bind_rows(regional_kendall_rows) %>%
        dplyr::mutate(index_type = as.character(index_type)) %>%
        dplyr::arrange(region_name, index_type, scale_months)
      
      write.csv(
        regional_kendall_table,
        file.path(REGIONAL_EVENT_DIR, "Table_KendallRP_Regional_AllRegions.csv"),
        row.names = FALSE
      )
      
      regional_kendall_primary <- regional_kendall_table %>%
        dplyr::select(
          region_id, region_name, index_type, index_label, scale_months,
          start_date, end_date, duration_months, severity,
          n_historical_events, mean_interarrival_years,
          tau_empirical, tau_fitted, duration_marginal, severity_marginal,
          best_copula, copula_theta, t_star, K_C,
          Kendall_return_period_years, CI_95_lo_years, CI_95_hi_years,
          bootstrap_median_years, bootstrap_converged
        )
      write.csv(
        regional_kendall_primary,
        file.path(REGIONAL_EVENT_DIR, "Table_KendallRP_Regional_Primary.csv"),
        row.names = FALSE
      )
      
      cat("\n  Regional Kendall tables saved:\n")
      cat("    Table_KendallRP_Regional_AllRegions.csv\n")
      cat("    Table_KendallRP_Regional_Primary.csv\n")
    }
  }
}

cat("\n╔══════════════════════════════════════╗\n")
cat("║  16event_ranking.R  COMPLETE        ║\n")
cat("╚══════════════════════════════════════╝\n\n")

cat("KEY FINDINGS FOR PRESENTATION:\n")
cat(sprintf("  • Total drought events (%s): %d\n", ANALYSIS_PERIOD_LABEL, nrow(primary_events)))
cat(sprintf(paste0("  • Most severe event overlapping ", RECENT_TAG,
                   ": rank #%d of %d (%s to %s, %d mo)\n"),
            msr_primary$rank_by_severity[1], nrow(primary_events),
            msr_primary$start_date[1], msr_primary$end_date[1],
            msr_primary$duration_months[1]))
cat(sprintf("  • Plotting-position percentile: %.1fth; more severe than %d of the other %d events\n",
            msr_primary$percentile_rank[1],
            sum(primary_events$severity < msr_primary$severity[1]),
            nrow(primary_events) - 1L))
cat(sprintf("  • Rank-based recurrence proxy: %.1f years (not annual exceedance)\n",
            msr_primary$return_period_years[1] %||% NA))
cat(sprintf("  • Duration class: %s\n",
            msr_primary$duration_class[1] %||% "N/A"))

cat("\nOUTPUT FILES:\n")
cat("  SW08 ranked catalog : ", file.path(RANKING_DIR, "ranked_event_catalog_SW08.csv"), "\n")
cat("  SW08 RP table       : ", file.path(RANKING_DIR, "return_period_by_duration_class_SW08.csv"), "\n")
cat("  Table 1a (short-scale indices): ", file.path(RANKING_DIR, "Table_1a_ShortScale_Index_TopEvents.csv"), "\n")
cat("  Excel summary       : ", file.path(RANKING_DIR, "Drought_Ranking_Summary.xlsx"), "\n")
cat("  Figure 1            : ", file.path(RANKING_DIR, "Fig1_Top_Droughts_By_Severity.png"), "\n")
cat("  Figure 2            : ", file.path(RANKING_DIR, "Fig2_Drought_Severity_Distribution.png"), "\n")
cat("  Figure 3            : ", file.path(RANKING_DIR, "Fig3_Drought_Frequency_By_Duration.png"), "\n")
cat("  Figure 4            : ", file.path(RANKING_DIR, "FigS_Event_Presentation_Dashboard.png"), "\n")
cat("  Figure RP           : ", file.path(RANKING_DIR, "FigRP_Return_Period_By_Duration_Class.png"), "\n")
cat("  Figure RP+recent      : ", file.path(RANKING_DIR, paste0("FigRP_Return_Period_By_Duration_Class_with_", RECENT_TAG, ".png")), "\n")
cat("  Fig3 MS (with CI)      : ", file.path(RANKING_DIR, "FigS_Event_Severity_Ranking.pdf"), "\n")
cat("  Fig3 MS (no CI)        : ", file.path(RANKING_DIR, "FigS_Event_Severity_Ranking_noCI.pdf"), "\n")
cat("  Fig3 Kendall (with CI) : ", file.path(RANKING_DIR, "FigS_Event_Severity_Ranking_KendallRP.pdf"), "\n")
cat("  Fig3 Kendall (no CI)   : ", file.path(RANKING_DIR, "FigS_Event_Severity_Ranking_KendallRP_noCI.pdf"), "\n")
cat("  Fig4d (with CI)        : ", file.path(RANKING_DIR, "Fig4d_KendallRP_MultiIndex_Comparison.png"), "\n")
cat("  Fig4d (no CI)          : ", file.path(RANKING_DIR, "Fig4d_KendallRP_MultiIndex_Comparison_noCI.png"), "\n")
cat("  Recent event bubbles PDF: ", file.path(RANKING_DIR, "Fig_Recent_Drought_Event_Bubbles.pdf"), "\n")
cat("  Recent event bubbles PNG: ", file.path(RANKING_DIR, "Fig_Recent_Drought_Event_Bubbles.png"), "\n")
cat("  Recent event bubble data: ", file.path(RANKING_DIR, "Fig_Recent_Drought_Event_Bubbles_Data.csv"), "\n")
cat("  Table multi-index      : ", file.path(RANKING_DIR, "Table_KendallRP_MultiIndex.csv"), "\n")
if (!is.null(kendall_rp_result)) {
  krp <- kendall_rp_result
  cat(sprintf("\nKENDALL RP SUMMARY:\n"))
  cat(sprintf("  Best copula       : %s (theta=%.3f)\n", krp$best_cop, krp$theta))
  cat(sprintf("  T_Kendall         : %.0f yr  [95%% CI: %.0f\u2013%.0f yr]\n",
              krp$T_kendall, krp$ci_lo, krp$ci_hi))
}
if (file.exists(file.path(RANKING_DIR, paste0("Fig5_", RECENT_TAG, "_Spatial_Severity.png"))))
  cat("  Figure 5          : ", file.path(RANKING_DIR, paste0("Fig5_", RECENT_TAG, "_Spatial_Severity.png")), "\n")
if (file.exists(file.path(RANKING_DIR, paste0("Fig6_", RECENT_TAG, "_Pixel_Rank_Map.png"))))
  cat("  Figure 6          : ", file.path(RANKING_DIR, paste0("Fig6_", RECENT_TAG, "_Pixel_Rank_Map.png")), "\n")
if (file.exists(file.path(RANKING_DIR, paste0("pixel_rank_", RECENT_TAG, ".nc"))))
  cat("  Rank raster       : ", file.path(RANKING_DIR, paste0("pixel_rank_", RECENT_TAG, ".nc")), "\n")

################################################################################
# END OF 16event_ranking.R
################################################################################