####################################################################################
# 9event_ranking.R  ·  DROUGHT EVENT RANKING & PET-FORMULATION ROBUSTNESS
#
# EVENT DEFINITION — framework-wide, from DROUGHT_ANALYSIS_utils.R (section F);
# identical to Scripts 6, 7, 11 and 12. Edit thresholds/durations ONLY there.
#   Onset       : index <  DROUGHT_ONSET  (currently -0.5)
#   Termination : index >= DROUGHT_END    (currently -0.5; hysteresis band closed)
#   Intensity   : I = mean(DROUGHT_END − x_t) over event months      [≥ 0]
#   Severity    : S = I × D  (Sheffield & Wood 2008 style)
#   Min duration: scale-dependent, EVENT_MIN_DURATION_BY_SCALE
#                 (currently 1→2, 3→3, 6→4, 12→6, 24→12, 36→12 months)
#   Missing months terminate an event (gaps are never bridged)
#   Duration classes: D3-6 (Short-term), D7-12 (Medium-term), D13+ (Long-term);
#                 D1-2 can occur only where the minimum duration is <= 2
#   Last monthly available data is for Dec 2025
#   Basin series : fractional-coverage area-weighted mean over the 7-subbasin
#                  union (cell area x coverage), reported only when finite cells
#                  hold >= BASIN_MIN_VALID_AREA_COVERAGE (0.20) of the basin weight.
#
#   This framework-wide event definition is used consistently
#   for SPI, SPEI-PM, SPEI-ERA5, SWEI, and SSPI. Snow indices are processed
#   only at their defined native scales (SWEI-03 and SSPI-01). The ranking uses only the manuscript
#   formulations retained in this script; no detrended/counterfactual branch
#   is used in this script.
#
#   References:
#     Sheffield & Wood (2008) Clim Dyn 31:79-105 — duration classes, I×D severity
#     Yevjevich (1967) — threshold-level run method (theoretical basis)
#
# MODIFICATIONS :
#    Part 0  — build_sw08_catalog(): reads raw NetCDF files,
#                            computes area-weighted basin average, applies SW08
#                            event loop; replaces load_event_catalog().
#    Part 1  — catalogs now built on-the-fly from the current NetCDF outputs.
#            Parts 4b/4c/4d — Fig 3 MS, Kendall copula RP, multi-index.
#
# Outputs (RANKING_DIR):
#   ranked_event_catalog_SW08.csv
#   return_period_by_duration_class_SW08.csv
#   Fig3_MS_Drought_Severity_Ranking.pdf/.png
#   Fig3_MS_Drought_Severity_Ranking_KendallRP.pdf/.png
#   Fig4c_Copula_Diagnostics.pdf
#   Fig4d_KendallRP_MultiIndex_Comparison.pdf/.png
#   Table_KendallRP_MultiIndex.csv
####################################################################################

# Set working directory BEFORE sourcing utils so that relative source() calls
# inside DROUGHT_ANALYSIS_utils.R (e.g. utils_teleconnection_addon.R) resolve.
setwd("D:/Nechako_Drought/Nechako/")

source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "ggplot2", "dplyr", "data.table", "lubridate",
                      "openxlsx", "patchwork", "scales", "sf", "viridis",
                      "gridExtra",      # summary table grobs
                      "copula",         # bivariate copula fitting (Part 4c)
                      "fitdistrplus",   # marginal distribution fitting (Part 4c)
                      "ggrepel"))       # retained for other project plotting utilities

# Define the %||% operator (used later for NULL/NA fallback)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x[1])) y else x

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

####################################################################################
# CONFIGURATION
####################################################################################
RANKING_DIR <- file.path(WD_PATH, "temporal_drought", "event_ranking")
dir.create(RANKING_DIR, showWarnings = FALSE, recursive = TRUE)

# Analysis period
HISTORICAL_START <- 1950
HISTORICAL_END   <- 2025
RECENT_DROUGHT_START <- as.Date("2022-01-01")
RECENT_DROUGHT_END   <- as.Date("2025-12-31")

stopifnot(
  identical(as.integer(MANUSCRIPT_SCALES), c(1L, 3L, 6L, 12L, 24L, 36L)),
  identical(as.integer(EVENT_RANKING_SCALES), as.integer(MANUSCRIPT_SCALES))
)


# ── Drought event identification parameters ───────────────────────────────────
# Threshold: −0.5 corresponds to the 31st percentile of the standard normal,
#   the most common climatological drought onset threshold (Lloyd-Hughes &
#   Saunders 2002; Spinoni et al. 2014, 2015).  It is more conservative (flags
#   more months) than the McKee et al. (1993) −1.0 threshold (16th percentile)
#   and corresponds to the lower tercile of the precipitation distribution.
#   Note: the SW08 10th-percentile threshold translates to −1.28 on a standard
#   normal scale, which is distinct from both −1.0 and −0.5.
DROUGHT_ONSET_EVENT <- DROUGHT_ONSET   # shared utility onset threshold
DROUGHT_END_EVENT   <- DROUGHT_END     # shared utility termination threshold
SWEI_SCALE <- 3L              # SWEI is produced at the 3-month scale in the snow pipeline
SSPI_SCALE <- 1L              # SSPI is produced at the 1-month scale in the snow pipeline
SSPI_SEAS_DIR <- file.path(WD_PATH, "sspi_results_seasonal/")

# Resolve each index to its authoritative seasonal NetCDF directory. SPEI-ERA5
# is catalogued as an observed formulation, but deliberately excluded from the
# Kendall-copula bootstrap sections (see the project runtime constraint).
get_event_index_dir <- function(index_type) {
  switch(tolower(index_type),
         spi       = SPI_SEAS_DIR,
         spei      = SPEI_SEAS_DIR,
         spei_era5 = SPEI_ERA5_SEAS_DIR,
         swei      = SWEI_SEAS_DIR,
         sspi      = SSPI_SEAS_DIR,
         stop("Unknown index_type: ", index_type)
  )
}

get_event_index_files <- function(index_type, scale) {
  idx <- tolower(index_type)
  seas_dir <- get_event_index_dir(idx)
  if (idx == "swei") {
    # The shared utility finder takes all three arguments; passing only the
    # directory (as older versions of this script did) raises an argument error.
    find_swei_seasonal_files(seas_dir, idx, as.integer(scale))
  } else {
    find_seasonal_nc_files(seas_dir, idx, as.integer(scale))
  }
}
# ── Event-duration rule ───────────────────────────────────────────────────────
# Scale-dependent minimum durations come from get_event_min_duration() in
# DROUGHT_ANALYSIS_utils.R (EVENT_MIN_DURATION_BY_SCALE). Runs shorter than the
# minimum are not events. Do not redefine get_event_min_duration() here.
MIN_DURATION_MAP <- vapply(as.character(EVENT_RANKING_SCALES),
                           function(s) get_event_min_duration(as.integer(s)), integer(1))

# Duration classes (Sheffield & Wood 2008)
DURATION_CLASSES <- list(
  short   = c(3, 6, "D3-6 (Short-term)"),
  medium  = c(7, 12, "D7-12 (Medium-term)"),
  long    = c(13, 999, "D13+ (Long-term)")
)

# All supported event-catalog indices and their valid scales. SPI, SPEI-PM,
# and SPEI-ERA5 use the six manuscript scales; snow indices use only their
# native scales produced by 8snow_drought_indices.R.
INDICES_TO_ANALYZE <- list(
  list(index = "spi",       scales = EVENT_RANKING_SCALES),
  list(index = "spei",      scales = EVENT_RANKING_SCALES),
  list(index = "spei_era5", scales = EVENT_RANKING_SCALES),
  list(index = "swei",      scales = c(SWEI_SCALE)),
  list(index = "sspi",      scales = c(SSPI_SCALE))
)

# Significance level for return period
ALPHA <- 0.05

####################################################################################
# HELPER FUNCTIONS
####################################################################################

#' Classify drought duration according to Sheffield & Wood (2008)
classify_duration <- function(duration_months) {
  sapply(duration_months, function(d) {
    if (d >= 3 & d <= 6) return("D3-6 (Short-term)")
    if (d >= 7 & d <= 12) return("D7-12 (Medium-term)")
    if (d >= 13) return("D13+ (Long-term)")
    return("D1-2 (Sub-threshold)")
  })
}

#' Calculate severity (intensity × duration)
calc_severity <- function(intensity, duration) {
  intensity * duration
}

#' Estimate return period from empirical frequency
calc_return_period <- function(event_rank, n_years_record) {
  # Weibull plotting position: T = (n+1)/m
  (n_years_record + 1) / event_rank
}

# ── Event table on a pre-computed time series vector ─────────────────────────
#
# Thin wrapper around the shared framework engine drought_event_table()
# (DROUGHT_ANALYSIS_utils.R): hysteresis run method, onset x < onset_thr,
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
  
  classify_dur <- function(d)
    dplyr::case_when(d >= 3L & d <= 6L ~ "D3-6 (Short-term)",
                     d >= 7L & d <= 12L ~ "D7-12 (Medium-term)",
                     d >= 13L ~ "D13+ (Long-term)",
                     TRUE ~ "D1-2 (Sub-threshold)")
  
  ev <- drought_event_table(vals, date_seq, min_duration = min_duration,
                            onset_threshold = onset_thr, termination_threshold = end_thr)
  
  out <- if (nrow(ev)) data.frame(
    start_date      = ev$start_date,
    end_date        = ev$end_date,
    start_year      = yr_seq[ev$start_index],
    duration_months = as.integer(ev$duration_months),
    duration_class  = classify_dur(ev$duration_months),
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


# ── Build SW08 catalog from raw NetCDF files for one index × scale ────────────
#
# Reads the seasonal NetCDF files (same layout as w1_trend_test.R), computes
# the area-weighted basin-average time series, then applies sw08_event_loop().
# Saves the catalog as CSV for reproducibility.
#
# Returns a data.frame or NULL if files cannot be found / read.
# Cached fractional-coverage basin weights per native grid geometry.
S14_BW_CACHE <- new.env(parent = emptyenv())
s14_basin_weights <- function(r_tmpl) {
  key <- paste(terra::nrow(r_tmpl), terra::ncol(r_tmpl),
               paste(signif(as.vector(terra::ext(r_tmpl)), 12), collapse = ","),
               trimws(terra::crs(r_tmpl)), sep = "|")
  if (!exists(key, envir = S14_BW_CACHE, inherits = FALSE))
    assign(key, basin_fractional_weights(r_tmpl, nechako_subbasin_union(load_nechako_subbasins())),
           envir = S14_BW_CACHE)
  get(key, envir = S14_BW_CACHE)
}

build_sw08_catalog <- function(index_type, scale,
                               onset_thr      = DROUGHT_ONSET_EVENT,
                               min_duration   = get_event_min_duration(scale),
                               start_year_rec = HISTORICAL_START,
                               end_year_rec   = HISTORICAL_END) {
  
  label    <- sprintf("%s-%02d", toupper(index_type), scale)
  min_dur  <- as.integer(min_duration)
  
  # ── locate NetCDF files using the index-specific native seasonal directory ──
  seas_dir <- get_event_index_dir(index_type)
  
  files <- tryCatch({
    get_event_index_files(index_type, scale)
  }, error = function(e) {
    cat(sprintf("  ⚠ File discovery failed for %s: %s\n", label, e$message))
    character(0)
  })
  
  if (length(files) == 0L) {
    cat(sprintf("  ⚠ %s: no NetCDF files found in %s\n", label, seas_dir))
    return(NULL)
  }
  
  # ── load rasters and build time-series matrix ──────────────────────────────
  stacks <- lapply(files, function(f)
    tryCatch(terra::rast(f), error = function(e) NULL))
  stacks <- Filter(Negate(is.null), stacks)
  if (length(stacks) == 0L) { cat(sprintf("  ⚠ %s: all rasters failed\n", label)); return(NULL) }
  
  r_tmpl         <- stacks[[1]][[1]]
  n_files        <- length(stacks)
  n_lyrs_per_file <- terra::nlyr(stacks[[1]])
  
  if (n_files == 1L && n_lyrs_per_file > 12L) {
    # single-file mode (MSPI/MSPEI)
    n_months <- n_lyrs_per_file
    dates_all <- extract_dates_from_nc(files[1], n_months)
    years_all <- sort(unique(as.integer(format(dates_all, "%Y"))))
    ts_mat    <- matrix(NA_real_, nrow = terra::ncell(r_tmpl), ncol = n_months)
    for (lyr in seq_len(n_months))
      ts_mat[, lyr] <- as.numeric(terra::values(stacks[[1]][[lyr]]))
    dates_ts <- dates_all
    years    <- years_all
  } else {
    # standard multi-file mode: 12 files, one per calendar month
    dates_m1 <- extract_dates_from_nc(files[1], n_lyrs_per_file)
    years    <- as.integer(format(dates_m1, "%Y"))
    n_years  <- length(years)
    ts_mat   <- matrix(NA_real_, nrow = terra::ncell(r_tmpl), ncol = n_years * 12L)
    for (m in seq_along(stacks)) {
      n_yrs_m <- min(terra::nlyr(stacks[[m]]), n_years)
      for (y in seq_len(n_yrs_m)) {
        ts_mat[, (y - 1L) * 12L + m] <- as.numeric(terra::values(stacks[[m]][[y]]))
      }
    }
    dates_ts <- as.Date(sprintf("%04d-%02d-01",
                                rep(years, each = 12L),
                                rep(1:12, times = length(years))))
  }
  
  # ── trim to requested period ───────────────────────────────────────────────
  yr_ts  <- as.integer(format(dates_ts, "%Y"))
  mo_col <- which(yr_ts >= start_year_rec & yr_ts <= end_year_rec)
  if (length(mo_col) == 0L) {
    cat(sprintf("  ⚠ %s: no months in %d-%d\n", label, start_year_rec, end_year_rec))
    return(NULL)
  }
  ts_mat   <- ts_mat[, mo_col, drop = FALSE]
  dates_ts <- dates_ts[mo_col]
  yr_ts    <- yr_ts[mo_col]
  
  # ── fractional-coverage area-weighted basin average ───────────────────────
  # Basin = 7-subbasin union; w = cell area x coverage (0 outside the basin);
  # a month is reported only when finite cells hold >= 0.20 of the basin weight.
  bw <- s14_basin_weights(r_tmpl)
  if (nrow(ts_mat) != length(bw$weights))
    stop(sprintf("%s: pixel matrix has %d rows but the native grid has %d cells.",
                 label, nrow(ts_mat), length(bw$weights)))
  basin_avg <- basin_weighted_mean(ts_mat, bw)
  
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
  
  events$index_type  <- index_type
  events$scale       <- scale
  events$index_label <- sprintf("%s_%02d", index_type, scale)
  attr(events, "open_event") <- open_event
  
  # save per-index CSV for reproducibility
  out_csv <- file.path(RANKING_DIR, sprintf("%s_SW08_event_catalog.csv", label))
  write.csv(events, out_csv, row.names = FALSE)
  
  cat(sprintf("  ✓ %s: %d events  [thr=%.2f, min_dur=%d m]  (saved %s)\n",
              label, nrow(events), onset_thr, min_dur, basename(out_csv)))
  events
}
#' Identify if event overlaps with 2022-2025 period
is_recent_drought <- function(start_date, end_date) {
  start_date <= RECENT_DROUGHT_END & end_date >= RECENT_DROUGHT_START
}

####################################################################################
# PART 1: BUILD EVENT CATALOGS FROM RAW NetCDF FILES
#
# Onset      : DROUGHT_ONSET_EVENT (= DROUGHT_ONSET in utils)
# Termination: DROUGHT_END_EVENT   (= DROUGHT_END in utils; hysteresis run method)
# Min duration: scale-dependent, EVENT_MIN_DURATION_BY_SCALE in utils
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 1: Building Event Catalogs from NetCDF\n")
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
  stop("No SW08 catalogs could be built. Check NetCDF file paths in DROUGHT_ANALYSIS_utils.R.")
}

# Merge all catalogs
master_events <- dplyr::bind_rows(all_events)

validate_scales(
  unique(master_events$scale[master_events$index_type %in% c("spi", "spei", "spei_era5")]),
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

# Mark recent drought (2022-2025)
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
    rank_by_severity = row_number(),
    return_period_years = calc_return_period(rank_by_severity, n_years_record),
    percentile_rank = round(100 * (1 - rank_by_severity / n()), 1)
  ) %>%
  ungroup()

# Identify the most severe recent event
recent_events <- master_events %>% filter(is_recent)
if (nrow(recent_events) > 0) {
  most_severe_recent <- recent_events %>%
    group_by(index_label) %>%
    slice_max(severity, n = 1) %>%
    ungroup()
  
  cat(sprintf("\n✓ Identified %d recent drought events (2022-2025)\n",
              nrow(recent_events)))
  cat(sprintf("  Most severe: %s (severity=%.2f, rank=%d)\n",
              most_severe_recent$index_label[1],
              most_severe_recent$severity[1],
              most_severe_recent$rank_by_severity[1]))
} else {
  cat("\n⚠ No events detected in 2022-2025 period\n")
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

####################################################################################
# PART 3b: DURATION-CLASS RETURN PERIOD TABLE — SW08  (Sheffield & Wood 2008)
#
# Key differences from the old (hysteresis) table:
#   • More short-term events: SW08 splits episodes whenever index recovers
#     above −0.5, so D3-6 count is expected to increase.
#   • Lower severity values: intensity measured as mean deficit below −0.5
#     (not below 0.0), so mean_I is ~10–50× smaller than hysteresis.
#   • No minimum duration filter: D1-2 events appear in the full catalog
#     but are excluded from the RP table here (same as SW08 Table 4).
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 3b: Duration-Class Return Period Table (SW08)\n")
cat("══════════════════════════════════════\n")

rp_table <- master_events %>%
  dplyr::filter(duration_class != "D1-2 (Sub-threshold)") %>%
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
    # Central Weibull return period
    return_period_yrs = round(n_years_record / n_events, 1),
    # Exact Poisson 95 % confidence interval (Wilson-Hilferty / chi-squared method):
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
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(0.6)) +
    ggplot2::scale_y_continuous(
      name   = "Return Period (years)",
      breaks = c(1, 2, 5, 10, 20, 50, 100),
      trans  = "log10",
      labels = scales::comma
    ) +
    ggplot2::scale_colour_viridis_d(option = "turbo", name = "Index", end = 0.95) +
    ggplot2::labs(
      title    = "Drought Return Period by Duration Class — Nechako Watershed",
      subtitle = sprintf("Weibull empirical estimate ± 95%% Poisson CI  |  Record: %d–%d",
                         HISTORICAL_START, HISTORICAL_END),
      x        = "Duration Class (Sheffield & Wood 2008)"
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
  
  # ── FigRP variant: same plot but with 2022-2025 drought overlaid ──────────────
  # For each index, the most severe event that overlaps 2022-2025 is marked as a
  # diamond (shape = 18).  Its return period comes from the Weibull ranking already
  # computed in Part 2 (most_severe_recent$return_period_years).
  # This lets the reader directly compare the 2022-2025 event against the empirical
  # frequency distribution of ALL historical droughts of the same duration class.
  if (!is.null(most_severe_recent) && nrow(most_severe_recent) > 0) {
    
    # Build a data frame of the 2022-2025 event per index, keeping only the
    # duration classes that appear in the historical rp_table axes
    valid_dc <- c("D3-6 (Short-term)", "D7-12 (Medium-term)", "D13+ (Long-term)")
    
    recent_rp_overlay <- most_severe_recent %>%
      dplyr::filter(!is.na(return_period_years),
                    duration_class %in% valid_dc) %>%
      dplyr::mutate(
        duration_class = factor(duration_class, levels = valid_dc),
        # Label: show both return period and severity rank
        rp_label = sprintf("2022–2025\n~%.0f yr", return_period_years)
      )
    
    if (nrow(recent_rp_overlay) > 0) {
      
      p_rp_2022 <- p_rp +
        # ── Large diamond for each index × duration-class combination ──────────
      # shape=18 (filled diamond) distinguishes the 2022-2025 event from the
      # historical-mean dots (shape default = 19, open circle)
      ggplot2::geom_point(
        data     = recent_rp_overlay,
        ggplot2::aes(x     = duration_class,
                     y     = return_period_years,
                     colour = index_label,
                     group  = index_label),
        size     = 5,
        shape    = 18,   # filled diamond
        position = ggplot2::position_dodge(0.6),
        show.legend = FALSE
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
        label.size  = 0.25,
        label.padding = ggplot2::unit(0.15, "lines"),
        position    = ggplot2::position_dodge(0.6),
        show.legend = FALSE
      ) +
        ggplot2::labs(
          title    = "Drought Return Period by Duration Class — Nechako Watershed",
          subtitle = sprintf(
            paste0("Weibull empirical estimate \u00b1 95%% Poisson CI  |  Record: %d\u2013%d  ",
                   " |  * = 2022\u20132025 event return period"),
            HISTORICAL_START, HISTORICAL_END
          )
        )
      
      ggplot2::ggsave(
        file.path(RANKING_DIR,
                  "FigRP_Return_Period_By_Duration_Class_with_2022-2025.png"),
        p_rp_2022, width = 11, height = 7, dpi = 300
      )
      cat("✓ Saved: FigRP_Return_Period_By_Duration_Class_with_2022-2025.png\n")
    } else {
      cat("  ⚠ 2022-2025 event has no valid duration class for RP overlay – skipped\n")
    }
  } else {
    cat("  ⚠ most_severe_recent is NULL – 2022-2025 RP overlay not drawn\n")
  }
}

####################################################################################
# PART 4: CREATE PRESENTATION FIGURES
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 4: Creating Presentation Figures\n")
cat("══════════════════════════════════════\n")

# Use SPI-3 as the primary index for presentation.
PRIMARY_INDEX <- "spi_03"
if (!(PRIMARY_INDEX %in% unique(master_events$index_label))) {
  PRIMARY_INDEX <- unique(master_events$index_label)[1]
  cat(sprintf("⚠ %s not available, using %s instead\n", "spi_03", PRIMARY_INDEX))
}

primary_events <- master_events %>% filter(index_label == PRIMARY_INDEX)

####################################################################################
# FIGURE 1: Top 10-15 Droughts by Severity (with 2022-2025 highlighted)
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
    plot_events <- dplyr::bind_rows(top_events[-top_n], recent_top) %>%
      arrange(desc(severity))
  } else {
    plot_events <- top_events
  }
}

plot_events$label <- sprintf("%s-%s\n(%d mo)",
                             format(plot_events$start_date, "%Y"),
                             format(plot_events$end_date, "%Y"),
                             plot_events$duration_months)

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
    subtitle = sprintf("%s | 1950-2025 | Red = 2022-2025 Drought", toupper(PRIMARY_INDEX)),
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
  mutate(cdf = row_number() / n())

# Get recent drought severity for annotation
if (!is.null(most_severe_recent)) {
  recent_sev <- most_severe_recent %>%
    filter(index_label == PRIMARY_INDEX) %>%
    pull(severity)
  recent_pct <- most_severe_recent %>%
    filter(index_label == PRIMARY_INDEX) %>%
    pull(percentile_rank)
} else {
  recent_sev <- NA
  recent_pct <- NA
}

p2 <- ggplot2::ggplot(cdf_data, ggplot2::aes(x = severity, y = cdf * 100)) +
  ggplot2::geom_line(color = "#2c3e50", linewidth = 1.2) +
  ggplot2::geom_point(color = "#2c3e50", size = 1.5, alpha = 0.5) +
  {if (!is.na(recent_sev))
    ggplot2::geom_vline(xintercept = recent_sev, linetype = "dashed",
                        color = "#e74c3c", linewidth = 1.5)} +
  {if (!is.na(recent_sev))
    ggplot2::annotate("segment",
                      x = recent_sev, xend = recent_sev,
                      y = 0, yend = recent_pct,
                      color = "#e74c3c", linewidth = 1, linetype = "dotted")} +
  {if (!is.na(recent_sev))
    ggplot2::annotate("point",
                      x = recent_sev, y = recent_pct,
                      color = "#e74c3c", size = 4, shape = 21, fill = "white", stroke = 1.5)} +
  ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  ggplot2::labs(
    title = "Historical Drought Severity Distribution",
    subtitle = sprintf("%s | 1950-2025 | 2022-2025 drought at %.1fth percentile",
                       toupper(PRIMARY_INDEX), round(recent_pct, 1)),
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
####################################################################################
# FIGURE 3: Drought Frequency by Duration Class (Sheffield & Wood 2008 style)
####################################################################################
cat("\n  Fig 3: Drought Frequency by Duration Class...\n")

duration_summary <- primary_events %>%
  dplyr::group_by(duration_class) %>%
  dplyr::summarise(
    n_events      = dplyr::n(),
    mean_duration = mean(duration_months, na.rm = TRUE),
    mean_severity = mean(severity,         na.rm = TRUE),  # <-- was mean(mean_S)
    max_severity  = max(severity,          na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    return_period = n_years_record / pmax(n_events, 1), # avoid divide-by-zero
    duration_class = factor(
      duration_class,
      levels = c("D1-2 (Sub-threshold)",
                 "D3-6 (Short-term)",
                 "D7-12 (Medium-term)",
                 "D13+ (Long-term)")
    )
  )

p3 <- ggplot2::ggplot(duration_summary,
                      ggplot2::aes(x = duration_class, y = n_events, fill = duration_class)) +
  ggplot2::geom_col(alpha = 0.8, color = "black", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("n=%d\nReturn: %.1f yrs", n_events, return_period)),
    vjust = -0.5, size = 4, fontface = "bold"
  ) +
  ggplot2::scale_fill_brewer(palette = "YlOrRd") +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
  ggplot2::labs(
    title = "Drought Frequency by Duration Class (Sheffield & Wood 2008)",
    subtitle = sprintf("%s Basin-Averaged | %d-%d | Nechako Watershed",
                       toupper(PRIMARY_INDEX), HISTORICAL_START, HISTORICAL_END),
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
  Metric = c("Total Drought Events (1950-2025)",
             "Short-term Droughts (3-6 mo)",
             "Medium-term Droughts (7-12 mo)",
             "Long-term Droughts (12+ mo)",
             "Most Severe Historical Drought",
             "2022-2025 Drought Rank",
             "2022-2025 Drought Percentile",
             "2022-2025 Return Period",
             "Mean Drought Duration",
             "Max Drought Duration"),
  Value = c(
    nrow(primary_events),
    sum(primary_events$duration_class == "D3-6 (Short-term)"),
    sum(primary_events$duration_class == "D7-12 (Medium-term)"),
    sum(primary_events$duration_class == "D13+ (Long-term)"),
    sprintf("%s-%s (Severity: %.1f)",
            format(max(primary_events$start_date), "%Y"),
            format(max(primary_events$end_date), "%Y"),
            max(primary_events$severity)),
    sprintf("#%d of %d",
            most_severe_recent$rank_by_severity[1] %||% NA,
            nrow(primary_events)),
    sprintf("%.1fth percentile", most_severe_recent$percentile_rank[1] %||% NA),
    sprintf("%.1f years", most_severe_recent$return_period_years[1] %||% NA),
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
    subtitle = "Prepared for April 2025 Presentation | Based on Sheffield & Wood (2008) Methodology",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(color = "gray50", size = 12, hjust = 0.5)
    )
  )

ggplot2::ggsave(file.path(RANKING_DIR, "Fig4_Presentation_Dashboard.png"),
                dashboard, width = 16, height = 12, dpi = 300)
cat("    ✓ Saved: Fig4_Presentation_Dashboard.png\n")

####################################################################################
# PART 4b: MANUSCRIPT FIGURE 3 — Drought Severity Ranking & Return Period
# (publication quality, two-panel, SPEI-6)
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 4b: Manuscript Figure 3\n")
cat("══════════════════════════════════════\n")

# ── Index for the manuscript figure ──────────────────────────────────────────
# SPEI-6 matches the composite severity score reported in Table 2.
# SPEI-3 is the manuscript index for Fig 3 (matches accumulation scales 1-3 scope).
# Change MS_INDEX if a different index/scale is preferred.
MS_INDEX <- "spei_03"
if (!(MS_INDEX %in% unique(master_events$index_label))) {
  MS_INDEX <- unique(master_events$index_label)[1]
  cat(sprintf("  ⚠ spei_03 not found; falling back to %s\n", MS_INDEX))
}
ms_events <- master_events %>% dplyr::filter(index_label == MS_INDEX)
cat(sprintf("  Using index: %s  (%d events)\n", MS_INDEX, nrow(ms_events)))

# ── Shared style constants ────────────────────────────────────────────────────
COL_HIST   <- "#4393c3"
COL_RECENT <- "#d6604d"
COL_JJA    <- "#fef0d9"
N_TOP      <- 15

# Display label: "spei_03" -> "SPEI-03"
MS_INDEX_LAB <- gsub("_", "-", toupper(MS_INDEX))

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

# Always include the most severe 2022-2025 event even if outside the top N
if (!any(top_ms$is_recent) && any(ms_events$is_recent)) {
  recent_top1 <- ms_events %>%
    dplyr::filter(is_recent) %>%
    dplyr::arrange(dplyr::desc(severity)) %>%
    dplyr::slice_head(n = 1)
  top_ms <- dplyr::bind_rows(top_ms[-N_TOP, ], recent_top1) %>%
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
    # Rank label shown to the right of each bar for 2022-2025
    rank_lab  = dplyr::if_else(
      is_recent,
      sprintf("Rank #%d  (%.1fth pct.)", rank_by_severity, percentile_rank),
      NA_character_
    ),
    fill_grp = dplyr::if_else(is_recent, "2022\u20132025 drought", "Historical droughts")
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
    values = c("Historical droughts"   = COL_HIST,
               "2022\u20132025 drought" = COL_RECENT),
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

# ── PANEL (b): Return period by duration class with 2022-2025 overlay ─────────
# Re-uses rp_table computed in Part 3b, filtered to the manuscript index.
# Falls back to all indices if MS_INDEX rows are absent (e.g., label mismatch).
rp_ms <- rp_table %>% dplyr::filter(index_label == MS_INDEX)
if (nrow(rp_ms) == 0) {
  cat("  ⚠ No rp_table rows for MS_INDEX; using full rp_table for panel (b)\n")
  rp_ms <- rp_table
}

# 2022-2025 overlay point
valid_dc_lvls <- c("D3-6 (Short-term)", "D7-12 (Medium-term)", "D13+ (Long-term)")
recent_ms_rp  <- NULL
if (!is.null(most_severe_recent) && nrow(most_severe_recent) > 0) {
  recent_ms_rp <- most_severe_recent %>%
    dplyr::filter(index_label == MS_INDEX,
                  !is.na(return_period_years),
                  duration_class %in% valid_dc_lvls) %>%
    dplyr::mutate(
      duration_class = factor(duration_class, levels = valid_dc_lvls),
      rp_label       = sprintf("2022\u20132025\n~%d yr",
                               round(return_period_years))
    )
}

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
    size   = 3.5,
    shape  = 16,
    colour = COL_HIST
  ) +
  # 2022-2025 diamond overlay (added only if data available)
  {if (!is.null(recent_ms_rp) && nrow(recent_ms_rp) > 0)
    ggplot2::geom_point(
      data   = recent_ms_rp,
      ggplot2::aes(x = duration_class, y = return_period_years),
      size   = 5,
      shape  = 18,   # filled diamond
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
      label.size  = 0.22,
      label.padding = ggplot2::unit(0.14, "lines"),
      inherit.aes = FALSE
    )} +
  ggplot2::scale_y_continuous(
    name   = "Return period (years)",
    trans  = "log10",
    breaks = c(1, 2, 5, 10, 20, 50, 100, 200),
    labels = scales::comma,
    expand = ggplot2::expansion(mult = c(0.08, 0.25))
  ) +
  ggplot2::scale_x_discrete(
    labels = c("D3\u20136\n(Short-term)",
               "D7\u201312\n(Medium-term)",
               "D13+\n(Long-term)")
  ) +
  ggplot2::labs(
    title    = "(b)  Empirical return period by drought duration class",
    subtitle = sprintf(
      "Weibull estimate \u00b1 95%% Poisson CI  |  %s  |  %d\u2013%d  |  * = 2022\u20132025 event",
      MS_INDEX_LAB, HISTORICAL_START, HISTORICAL_END
    ),
    x = "Duration class",
    y = "Return period (years)"
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
  ggplot2::geom_point(size = 3.5, shape = 16, colour = COL_HIST) +
  {if (!is.null(recent_ms_rp) && nrow(recent_ms_rp) > 0)
    ggplot2::geom_point(
      data   = recent_ms_rp,
      ggplot2::aes(x = duration_class, y = return_period_years),
      size = 5, shape = 18, colour = COL_RECENT, inherit.aes = FALSE
    )} +
  {if (!is.null(recent_ms_rp) && nrow(recent_ms_rp) > 0)
    ggplot2::geom_label(
      data   = recent_ms_rp,
      ggplot2::aes(x = duration_class, y = return_period_years, label = rp_label),
      vjust = -0.45, size = 2.8, fontface = "bold", colour = COL_RECENT,
      label.size = 0.22, label.padding = ggplot2::unit(0.14, "lines"),
      inherit.aes = FALSE
    )} +
  ggplot2::scale_y_continuous(
    name   = "Return period (years)", trans = "log10",
    breaks = c(1, 2, 5, 10, 20, 50, 100, 200),
    labels = scales::comma,
    expand = ggplot2::expansion(mult = c(0.08, 0.25))
  ) +
  ggplot2::scale_x_discrete(
    labels = c("D3\u20136\n(Short-term)", "D7\u201312\n(Medium-term)", "D13+\n(Long-term)")
  ) +
  ggplot2::labs(
    title    = "(b)  Empirical return period by drought duration class",
    subtitle = sprintf(
      "Weibull point estimate only (no CI)  |  %s  |  %d\u2013%d  |  * = 2022\u20132025 event",
      MS_INDEX_LAB, HISTORICAL_START, HISTORICAL_END),
    x = "Duration class", y = "Return period (years)"
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
        ggplot2::ggsave(path, plot = fig, width = w, height = h, units = "in", device = "pdf")
      else
        ggplot2::ggsave(path, plot = fig, width = w, height = h, units = "in",
                        dpi = 300, device = "png")
      cat(sprintf("  ✓ Saved: %s\n", basename(path)))
    }, error = function(e) cat(sprintf("  ⚠ Save failed (%s): %s\n", ext, e$message)))
  }
}

save_fig(fig3_ms,      "Fig3_MS_Drought_Severity_Ranking")
save_fig(fig3_ms_noci, "Fig3_MS_Drought_Severity_Ranking_noCI")

####################################################################################
####################################################################################
# PART 4c: CANONICAL KENDALL-BASED JOINT RETURN PERIOD ENGINE
#
# This is the SINGLE Kendall implementation used by this script.
#
# The same run_kendall_rp() function is called by:
#   • Part 4c: manuscript SPEI-3 Kendall RP + diagnostics;
#   • Part 4d: multi-index comparison;
#   • Part 4e: all manuscript scales.
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
# Log-logistic is fitted through fitdistrplus using correctly registered
# d/p/q/r functions. The previous duplicate implementations attempted to pass
# ddist/pdist directly to fitdist(), which is not a valid fitdist interface and
# caused the log-logistic candidate to fail systematically.
####################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4c: Canonical Kendall RP engine\n")
cat("══════════════════════════════════════════════════════════════\n")

KENDALL_BOOT_DEFAULT <- 500L
# Console diagnostics for the computationally expensive bootstrap.
# Progress is reported without changing the statistical calculation.
KENDALL_PROGRESS_EVERY <- 10L
KENDALL_ALPHA <- 0.05
KENDALL_RECENT_START <- RECENT_DROUGHT_START
KENDALL_RECENT_END   <- RECENT_DROUGHT_END
KENDALL_MAX_RP <- 1e7

# Bootstrap stability controls. The point-estimate copula is still selected by
# ML/AIC above. For the bootstrap only, dependence is re-estimated by Kendall
# tau inversion (ITAU) rather than a second ML optimization. This avoids
# pathological optimizer behaviour in very small, highly concordant samples
# such as SPI-36 (n = 13). Near-perfect bootstrap concordance is treated as a
# degenerate replicate rather than clipped to an artificial tau = 0.999.
KENDALL_BOOT_COPULA_METHOD <- "itau"
KENDALL_BOOT_TAU_MIN <- 0.01
KENDALL_BOOT_TAU_MAX <- 0.98
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
  tau_init <- max(0.01, min(0.999, as.numeric(tau_emp)))
  
  templates <- list(
    Clayton = copula::claytonCopula(dim = 2L),
    Gumbel  = copula::gumbelCopula(dim = 2L),
    Frank   = copula::frankCopula(dim = 2L),
    Joe     = copula::joeCopula(dim = 2L),
    AMH     = copula::amhCopula(dim = 2L)
  )
  
  safe_fit <- function(template) {
    tryCatch({
      theta0 <- copula::iTau(template, tau_init)
      copula::fitCopula(
        copula::setTheta(template, theta0),
        uv,
        method = "ml"
      )
    }, error = function(e) NULL)
  }
  
  fits <- Filter(Negate(is.null), lapply(templates, safe_fit))
  if (!length(fits)) stop("All Kendall copula fits failed.")
  
  aic <- vapply(fits, stats::AIC, numeric(1))
  best_name <- names(which.min(aic))[1]
  
  list(
    fits = fits,
    aic = aic,
    best_name = best_name,
    best_fit = fits[[best_name]],
    best_copula = fits[[best_name]]@copula
  )
}

# ----------------------------------------------------------------------------
# Canonical one-index Kendall RP analysis.
# This function is the ONLY Kendall engine used in Parts 4c, 4d and 4e.
# -----------------------------------------------------------------------------
run_kendall_rp <- function(
    idx_label,
    n_boot = KENDALL_BOOT_DEFAULT,
    seed = 42L,
    events_data = NULL,
    record_years = n_years_record) {
  
  tryCatch({
    event_source <- if (is.null(events_data)) master_events else events_data
    record_years <- as.numeric(record_years)
    if (!is.finite(record_years) || record_years <= 0)
      stop("record_years must be a positive finite number")
    
    ev <- event_source %>%
      dplyr::filter(
        .data$index_label == .env$idx_label,
        .data$duration_class != "D1-2 (Sub-threshold)",
        is.finite(duration_months),
        is.finite(severity),
        severity > 0
      ) %>%
      dplyr::arrange(start_date, end_date)
    
    n_ev <- nrow(ev)
    if (n_ev < 8L)
      stop("too few qualifying historical events: ", n_ev)
    
    duration <- as.numeric(ev$duration_months)
    severity <- as.numeric(ev$severity)
    tau_emp <- suppressWarnings(stats::cor(duration, severity, method = "kendall"))
    if (!is.finite(tau_emp)) tau_emp <- 0
    
    cat(sprintf(
      "  Events for Kendall fit: %d | mean inter-arrival: %.2f yr | tau=%.4f\n",
      n_ev,
      record_years / n_ev,
      tau_emp
    ))
    
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
    
    # Target events must lie within the D3+ duration domain used to fit the
    # Kendall marginals. D1-2 events remain in the empirical event catalog but
    # are not extrapolated through marginals fitted only to D3+ events.
    recent <- event_source %>%
      dplyr::filter(
        .data$index_label == .env$idx_label,
        .data$is_recent,
        .data$duration_class != "D1-2 (Sub-threshold)",
        is.finite(duration_months),
        duration_months >= 3,
        is.finite(severity),
        severity > 0
      ) %>%
      dplyr::arrange(dplyr::desc(severity), start_date, end_date)
    
    if (!nrow(recent))
      stop("no 2022-2025 event for this index")
    
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
      u0 <- kendall_pit(target_duration, marg$duration, marg$duration_name)
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
          kendall_pit(recent$duration_months[1], marg$duration, marg$duration_name),
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
    
    bootstrap_target <- function(target_duration, target_severity, target_rank = NA_integer_) {
      boot_t0 <- proc.time()[["elapsed"]]
      cat(sprintf(
        "    Bootstrap start: target rank %s | %d replicates | copula bootstrap method=%s | timeout=%ds\n",
        ifelse(is.finite(target_rank), as.character(target_rank), "?"),
        n_boot,
        toupper(KENDALL_BOOT_COPULA_METHOD),
        KENDALL_BOOT_REPLICATE_TIMEOUT_SEC
      ))
      flush.console()
      
      tb <- rep(NA_real_, n_boot)
      tau_b_store <- rep(NA_real_, n_boot)
      theta_b_store <- rep(NA_real_, n_boot)
      status_store <- rep("failed", n_boot)
      
      for (b in seq_len(n_boot)) {
        tau_b_current <- NA_real_
        theta_b_current <- NA_real_
        timeout_flag <- FALSE
        high_tau_flag <- FALSE
        
        tb[b] <- tryCatch({
          # A single pathological copula/marginal fit must not hold the entire
          # 500-replicate bootstrap indefinitely. setTimeLimit() is reset at the
          # end of every replicate below.
          setTimeLimit(
            cpu = KENDALL_BOOT_REPLICATE_TIMEOUT_SEC,
            elapsed = KENDALL_BOOT_REPLICATE_TIMEOUT_SEC,
            transient = TRUE
          )
          on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE),
                  add = TRUE)
          
          ii <- sample.int(n_ev, n_ev, replace = TRUE)
          db <- duration[ii]
          sb <- severity[ii]
          
          fd_b <- refit_kendall_marginal(
            db,
            marg$duration_name,
            marg$duration
          )
          fs_b <- refit_kendall_marginal(
            sb,
            marg$severity_name,
            marg$severity
          )
          if (is.null(fd_b) || is.null(fs_b)) stop("marginal refit failed")
          
          ub <- kendall_pit(db, fd_b, marg$duration_name)
          vb <- kendall_pit(sb, fs_b, marg$severity_name)
          
          tau_b <- suppressWarnings(stats::cor(db, sb, method = "kendall"))
          if (!is.finite(tau_b)) stop("bootstrap Kendall tau is non-finite")
          tau_b_current <- as.numeric(tau_b)
          
          # Do not clip tau=1 (or an equally extreme value) to 0.999. Such a
          # value represents a degenerate/near-degenerate resample for this
          # small event sample and can imply an impractically large Frank theta.
          if (tau_b < KENDALL_BOOT_TAU_MIN || tau_b > KENDALL_BOOT_TAU_MAX) {
            high_tau_flag <- tau_b > KENDALL_BOOT_TAU_MAX
            stop(sprintf(
              "degenerate bootstrap tau=%.6f outside [%.2f, %.2f]",
              tau_b, KENDALL_BOOT_TAU_MIN, KENDALL_BOOT_TAU_MAX
            ))
          }
          
          theta_b <- copula::iTau(best_cop, tau_b)
          if (!is.finite(theta_b)) stop("bootstrap copula theta is non-finite")
          theta_b_current <- as.numeric(theta_b)
          cop_start_b <- copula::setTheta(best_cop, theta_b)
          
          # The full ML fit is retained for the original copula/AIC selection.
          # Bootstrap refits use ITAU because tau_b is already the dependence
          # statistic being resampled; this removes the second numerical ML
          # optimization that can become pathological for n_ev = 13.
          fit_b <- copula::fitCopula(
            cop_start_b,
            cbind(ub, vb),
            method = KENDALL_BOOT_COPULA_METHOD
          )
          cop_b <- fit_b@copula
          
          u0_b <- kendall_pit(target_duration, fd_b, marg$duration_name)
          v0_b <- kendall_pit(target_severity, fs_b, marg$severity_name)
          t_b <- as.numeric(copula::pCopula(c(u0_b, v0_b), cop_b))
          kc_b <- kendall_KC(t_b, cop_b)
          
          status_store[b] <- "converged"
          mu_ia / max(1e-12, 1 - kc_b)
        }, error = function(e) {
          msg <- conditionMessage(e)
          if (b <= 3L) {
            cat(sprintf("      bootstrap replicate %d error: %s\n", b, msg))
          }
          timeout_flag <<- grepl("time limit|elapsed time limit|CPU time limit",
                                 msg, ignore.case = TRUE)
          if (timeout_flag) {
            status_store[b] <<- "timeout"
          } else if (high_tau_flag) {
            status_store[b] <<- "high_tau"
          } else {
            status_store[b] <<- "failed"
          }
          NA_real_
        })
        
        tau_b_store[b] <- tau_b_current
        theta_b_store[b] <- theta_b_current
        
        # Ensure the time limit cannot leak into the next replicate.
        setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
        
        if (b == 1L || b == n_boot || (b %% KENDALL_PROGRESS_EVERY) == 0L) {
          elapsed <- proc.time()[["elapsed"]] - boot_t0
          rate <- b / max(elapsed, 1e-9)
          eta <- (n_boot - b) / max(rate, 1e-9)
          cat(sprintf(
            "      bootstrap %d/%d (%.0f%%) | elapsed %.1fs | ETA %.1fs | converged %d | failed %d | high-tau %d | timeout %d\n",
            b, n_boot, 100 * b / n_boot, elapsed, eta,
            sum(status_store[seq_len(b)] == "converged"),
            sum(status_store[seq_len(b)] == "failed"),
            sum(status_store[seq_len(b)] == "high_tau"),
            sum(status_store[seq_len(b)] == "timeout")
          ))
          flush.console()
        }
      }
      
      tb <- tb[is.finite(tb) & tb < KENDALL_MAX_RP]
      summary <- if (length(tb) < max(50L, floor(n_boot * 0.5))) {
        c(ci_lo = NA_real_, ci_hi = NA_real_, ci_med = NA_real_, n_conv = length(tb))
      } else {
        c(
          ci_lo = unname(stats::quantile(tb, 0.025, names = FALSE)),
          ci_hi = unname(stats::quantile(tb, 0.975, names = FALSE)),
          ci_med = unname(stats::quantile(tb, 0.500, names = FALSE)),
          n_conv = length(tb)
        )
      }
      list(
        summary = summary,
        values = tb,
        tau = tau_b_store,
        theta = theta_b_store,
        status = status_store
      )
    }
    
    primary_boot <- numeric(0)
    
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
      b <- bootstrap_target(
        recent$duration_months[j],
        recent$severity[j],
        target_rank = j
      )
      if (j == 1L) primary_boot <<- b$values
      b_summary <- b$summary
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
        bootstrap_tau_max = KENDALL_BOOT_TAU_MAX,
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

# Cache results so an index analyzed in Part 4c is not refit in Part 4d/4e.
KENDALL_RESULT_CACHE <- new.env(parent = emptyenv())
get_kendall_rp <- function(idx_label, n_boot = KENDALL_BOOT_DEFAULT) {
  if (exists(idx_label, envir = KENDALL_RESULT_CACHE, inherits = FALSE)) {
    return(get(idx_label, envir = KENDALL_RESULT_CACHE, inherits = FALSE))
  }
  ans <- run_kendall_rp(idx_label, n_boot = n_boot)
  assign(idx_label, ans, envir = KENDALL_RESULT_CACHE)
  ans
}

# Part 4c manuscript index: SPEI-3.
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
      label.size    = 0.2
    )
  }
  
  # ── WITH CI ───────────────────────────────────────────────────────────────
  pb3_k <- pb3 +
    annotate_label(k_lab_ci) +
    ggplot2::labs(subtitle = sprintf(
      "Weibull Poisson CI (points) + Kendall copula RP [with bootstrap CI] | %s | %d\u2013%d",
      MS_INDEX_LAB, HISTORICAL_START, HISTORICAL_END))
  
  fig3_k <- pa3 / pb3_k +
    patchwork::plot_layout(heights = c(1.35, 1)) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 4, 4)))
  
  fig3k_pdf <- file.path(RANKING_DIR, "Fig3_MS_Drought_Severity_Ranking_KendallRP.pdf")
  fig3k_png <- file.path(RANKING_DIR, "Fig3_MS_Drought_Severity_Ranking_KendallRP.png")
  tryCatch({
    ggplot2::ggsave(fig3k_pdf, fig3_k, width=7.0, height=7.5, units="in", device="pdf")
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
      "Weibull Poisson CI (points) + Kendall copula RP [point estimate only] | %s | %d\u2013%d",
      MS_INDEX_LAB, HISTORICAL_START, HISTORICAL_END))
  
  fig3_k_noci <- pa3 / pb3_k_noci +
    patchwork::plot_layout(heights = c(1.35, 1)) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(plot.margin = ggplot2::margin(4, 4, 4, 4)))
  
  fig3k_noci_pdf <- file.path(RANKING_DIR, "Fig3_MS_Drought_Severity_Ranking_KendallRP_noCI.pdf")
  fig3k_noci_png <- file.path(RANKING_DIR, "Fig3_MS_Drought_Severity_Ranking_KendallRP_noCI.png")
  tryCatch({
    ggplot2::ggsave(fig3k_noci_pdf, fig3_k_noci, width=7.0, height=7.5, units="in", device="pdf")
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
  diag4c_pdf <- file.path(RANKING_DIR, "Fig4c_Copula_Diagnostics.pdf")
  tryCatch({
    grDevices::pdf(diag4c_pdf, width = 12, height = 9)
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
                        krp$best_cop, krp$theta, krp$tau_fitted))
    graphics::points(krp$u0, krp$v0, pch = 18, cex = 2.4, col = "#d6604d")
    graphics::legend("bottomright",
                     legend = c("Historical", "2022-2025"),
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
# selected SPI and SPEI-PM scales plus SWEI-03 and SSPI-01. SPEI-ERA5 is
# excluded from Kendall bootstrap calculations under the project's runtime constraint.
#
# Purpose: assess whether the near-zero Kendall tau (duration, severity)
#   observed for SPEI-06 is index-specific or a general feature of this
#   drought catalog. Indices with longer accumulation scales accumulate
#   fewer, longer events (higher tau) while short-scale indices produce
#   many short events with weaker duration-severity coupling.
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
# Compare the manuscript SPI/SPEI-PM short scales plus the native snow scales.
# SPEI-ERA5 is intentionally excluded from Kendall bootstrap calculations to
# respect the project's prior runtime constraint; its event catalogs and
# empirical ranking products are still generated above.
MULTI_INDICES <- c(
  "spi_01", "spi_03", "spi_06",
  "spei_01", "spei_03", "spei_06",
  "swei_03", "sspi_01"
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
  if (grepl("^spei_era5_", idx_label)) return(sprintf("SPEI-ERA5-%d", sc))
  if (grepl("^spei_", idx_label)) return(sprintf("SPEI-PM-%d", sc))
  if (grepl("^swei_", idx_label)) return(sprintf("SWEI-%d", sc))
  if (grepl("^sspi_", idx_label)) return(sprintf("SSPI-%d", sc))
  if (grepl("^spi_", idx_label)) return(sprintf("SPI-%d", sc))
  idx_label
}

# Stable colours for the multi-index comparison.
MULTI_INDEX_COLOURS <- c(
  spi_01 = "#fc8d59", spi_03 = "#d7301f", spi_06 = "#7f0000",
  spei_01 = "#74add1", spei_03 = "#4575b4", spei_06 = "#313695",
  swei_03 = "#238b45", sspi_01 = "#8856a7"
)

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
                r$t_star, r$Kc, r$TK, r$ci_lo, r$ci_hi, r$n_boot_conv, 500L))
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
      data.frame(index   = gsub("_","-",toupper(r$index)),
                 metric  = "Kendall tau\n(duration, severity)",
                 value   = r$tau_emp,
                 stringsAsFactors = FALSE),
      data.frame(index   = gsub("_","-",toupper(r$index)),
                 metric  = "t* = C(u0,v0)",
                 value   = r$t_star,
                 stringsAsFactors = FALSE),
      data.frame(index   = gsub("_","-",toupper(r$index)),
                 metric  = "K_C(t*)",
                 value   = r$Kc,
                 stringsAsFactors = FALSE)
    )
  }))
  long_df$index  <- factor(long_df$index,
                           levels = gsub("_","-",toupper(names(valid_res))))
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
    ggplot2::ggsave(fig4d_pdf, fig4d, width=7.5, height=8.0, units="in", device="pdf")
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
    ggplot2::ggsave(fig4d_noci_pdf, fig4d_noci, width=7.5, height=8.0, units="in", device="pdf")
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
# PART 5: SPATIAL MAP OF 2022-2025 DROUGHT (pixel‑level severity)
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 5: Spatial Analysis of 2022-2025 Drought\n")
cat("══════════════════════════════════════\n")

# Load basin as SpatVector for raster extraction
basin_vect <- tryCatch(
  load_basin_vect(BASIN_SHP),
  error = function(e) NULL
)

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
      for (f in all_files) {
        r <- terra::rast(f)
        dates <- extract_dates_from_nc(f, terra::nlyr(r))
        # Keep only layers with date between start_d and end_d
        keep <- which(dates >= start_d & dates <= end_d)
        if (length(keep) > 0) {
          ts_list[[length(ts_list)+1]] <- r[[keep]]
        }
      }
      
      if (length(ts_list) == 0) {
        cat("  ⚠ No data found for the exact event period\n")
      } else {
        # Combine all selected layers into a single SpatRaster
        event_stack <- do.call(c, ts_list)
        cat(sprintf("  Extracted %d months of data for the event\n", terra::nlyr(event_stack)))
        
        # Mask to basin
        event_stack <- terra::mask(event_stack, basin_vect)
        
        # Compute per-pixel metrics
        r_mean <- terra::app(event_stack, fun = mean, na.rm = TRUE)
        r_min  <- terra::app(event_stack, fun = min, na.rm = TRUE)
        r_dur  <- terra::app(event_stack, fun = function(x) sum(x < DROUGHT_ONSET, na.rm = TRUE))
        
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
        
        # Convert basin to sf for plotting
        basin_sf <- sf::st_as_sf(basin_vect)
        
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
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", size = 0.6) +
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
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", size = 0.6) +
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
          ggplot2::geom_sf(data = basin_sf, fill = NA, color = "black", size = 0.6) +
          ggplot2::coord_sf() +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = "Drought Duration within Event") +
          ggplot2::theme(legend.position = "bottom")
        
        # Combine panels
        p_spatial <- p_mean + p_min + p_dur +
          patchwork::plot_annotation(
            title = "2022-2025 Drought: Spatial Characteristics",
            subtitle = sprintf("%s | Pixel-level metrics", toupper(PRIMARY_INDEX))
          )
        
        ggplot2::ggsave(file.path(RANKING_DIR, "Fig5_2022-2025_Spatial_Severity.png"),
                        p_spatial, width = 15, height = 6, dpi = 300)
        cat("    ✓ Saved: Fig5_2022-2025_Spatial_Severity.png\n")
      }
    }
  } else {
    cat("  ⚠ Could not isolate recent event for spatial mapping\n")
  }
} else {
  cat("  ⚠ Basin shapefile not found or recent event missing – skipping spatial map\n")
  # Fallback placeholder
  if (!is.null(basin_sf)) {
    p_spatial <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = basin_sf, fill = "#ecf0f1", color = "#2c3e50", linewidth = 0.8) +
      ggplot2::annotate("text",
                        x = sf::st_coordinates(sf::st_centroid(basin_sf))[1],
                        y = sf::st_coordinates(sf::st_centroid(basin_sf))[2],
                        label = sprintf("2022-2025 Drought\nSeverity: %.1f\nRank: #%d",
                                        most_severe_recent$severity[1] %||% NA,
                                        most_severe_recent$rank_by_severity[1] %||% NA),
                        size = 5, fontface = "bold", color = "#c0392b",
                        hjust = 0.5, vjust = 0.5) +
      ggplot2::labs(
        title = "2022-2025 Drought Spatial Context (Placeholder)",
        subtitle = "Nechako Watershed | Severity ranking relative to 1950-2025"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
        panel.grid = ggplot2::element_blank()
      )
    
    ggplot2::ggsave(file.path(RANKING_DIR, "Fig5_2022-2025_Spatial_Context_placeholder.png"),
                    p_spatial, width = 10, height = 8, dpi = 300)
    cat("    ✓ Saved placeholder map (actual pixel data unavailable)\n")
  }
}

####################################################################################
# PART 5b: 2022-2025 PIXEL-LEVEL SEVERITY RANK MAP
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 5b: 2022-2025 Pixel-Level Severity Rank Map\n")
cat("══════════════════════════════════════\n")

# Helper: extract events from a single pixel time-series with the shared
# framework engine (same definition as the basin catalog built in Part 1).
# Missing months terminate an event; they are no longer dropped and bridged.
extract_pixel_events_sw <- function(x, onset = DROUGHT_ONSET_EVENT,
                                    min_duration,
                                    end_thr = DROUGHT_END_EVENT) {
  ev <- drought_event_table(x, min_duration = min_duration,
                            onset_threshold = onset, termination_threshold = end_thr)
  data.frame(dur = ev$duration_months, I = ev$mean_int, S = ev$severity)
}

# Resolve primary index type and scale
idx_type_p  <- sub("_[0-9]+$", "", PRIMARY_INDEX)
idx_scale_p <- as.integer(sub(".*_([0-9]+)$", "\\1", PRIMARY_INDEX))

# Part 5b reads the time-series RDS written by w1_trend_test.R (not the results
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
    
    col_2022_start  <- (2022L - record_start_yr) * 12L + 1L
    col_2025_end    <- min((2025L - record_start_yr + 1L) * 12L, n_tcols)
    
    if (col_2022_start < 1L || col_2022_start > n_tcols) {
      cat("  ⚠ 2022-2025 columns out of range for this record – skipping rank map\n")
    } else {
      cat(sprintf("  Computing per-pixel event catalogs (onset < %.1f)...\n", DROUGHT_ONSET))
      
      rank_vec  <- rep(NA_real_, n_pix_p)
      sev_2022  <- rep(NA_real_, n_pix_p)
      n_hist    <- rep(NA_integer_, n_pix_p)
      
      for (i in seq_len(n_pix_p)) {
        full_ts   <- ts_mat_p[i, ]
        hist_ts   <- full_ts
        
        idx_22_25 <- seq(col_2022_start, col_2025_end)
        ts_22_25  <- full_ts[idx_22_25]
        
        ev_all    <- extract_pixel_events_sw(hist_ts,
                                             onset        = DROUGHT_ONSET_EVENT,
                                             min_duration = get_event_min_duration(idx_scale_p))
        ev_recent <- extract_pixel_events_sw(ts_22_25,
                                             onset        = DROUGHT_ONSET_EVENT,
                                             min_duration = get_event_min_duration(idx_scale_p))
        
        if (nrow(ev_recent) == 0L || nrow(ev_all) == 0L) next
        
        s_recent <- max(ev_recent$S, na.rm = TRUE)
        sev_2022[i] <- s_recent
        
        n_less      <- sum(ev_all$S < s_recent, na.rm = TRUE)
        rank_vec[i] <- 100 * n_less / nrow(ev_all)
        n_hist[i]   <- nrow(ev_all)
      }
      
      # Assemble output data frame
      rank_df <- data.frame(
        lon           = lon_p,
        lat           = lat_p,
        severity_2022 = sev_2022,
        pct_rank_2022 = rank_vec,
        n_hist_events = n_hist
      )
      write.csv(rank_df,
                file.path(RANKING_DIR, "pixel_rank_2022-2025.csv"),
                row.names = FALSE)
      cat("  ✓ Saved: pixel_rank_2022-2025.csv\n")
      
      # Build rasters
      r_base <- tryCatch(terra::rast(r_tmpl_path <- file.path(
        TREND_DIR, sprintf("%s_%02d_D13p_pixel_frequency.nc", idx_type_p, idx_scale_p)
      )), error = function(e) NULL)
      
      if (is.null(r_base)) {
        r_base <- terra::rast(
          xmin = min(rank_df$lon) - 0.05, xmax = max(rank_df$lon) + 0.05,
          ymin = min(rank_df$lat) - 0.05, ymax = max(rank_df$lat) + 0.05,
          resolution = 0.1, crs = "EPSG:4326"
        )
        cat("  ℹ Using synthetic raster extent (D13p NC not found)\n")
      }
      
      fill_raster <- function(base_r, df, col_name) {
        r_out <- base_r[[1]]
        terra::values(r_out) <- NA_real_
        cell_ids <- terra::cellFromXY(r_out, cbind(df$lon, df$lat))
        valid    <- !is.na(cell_ids)
        terra::values(r_out)[cell_ids[valid]] <- df[[col_name]][valid]
        r_out
      }
      
      r_rank <- fill_raster(r_base, rank_df, "pct_rank_2022")
      r_sev  <- fill_raster(r_base, rank_df, "severity_2022")
      
      names(r_rank) <- "pct_rank_2022_2025"
      names(r_sev)  <- "severity_2022_2025"
      
      terra::writeCDF(
        c(r_rank, r_sev),
        file.path(RANKING_DIR, "pixel_rank_2022-2025.nc"),
        overwrite = TRUE
      )
      cat("  ✓ Saved: pixel_rank_2022-2025.nc\n")
      
      # Two-panel PNG map
      df_rank <- as.data.frame(r_rank, xy = TRUE, na.rm = TRUE)
      df_sev  <- as.data.frame(r_sev,  xy = TRUE, na.rm = TRUE)
      colnames(df_rank)[3] <- "pct_rank"
      colnames(df_sev)[3]  <- "severity"
      
      basin_sf_plot <- tryCatch(
        sf::st_as_sf(load_basin_vect(BASIN_SHP)),
        error = function(e) NULL
      )
      
      # ──────────────────────────────────────────────────────────────────────
      # WHAT DOES "PERCENTILE RANK" MEAN IN THIS MAP?
      # ──────────────────────────────────────────────────────────────────────
      # For every grid pixel, we extracted the full time-series of the primary
      # drought index (e.g. SPI-12) from HISTORICAL_START to HISTORICAL_END.
      # Using the Sheffield & Wood (2008) event-identification algorithm
      # (drought onset < DROUGHT_ONSET, termination ≥ DROUGHT_END) we built a
      # catalogue of ALL discrete drought events at that pixel over the entire
      # record.  Each event has a severity score S = mean_intensity × duration
      # (where intensity = onset_threshold − index_value, always positive).
      #
      # For the 2022-2025 period we then isolated the sub-series (time columns
      # col_2022_start : col_2025_end) and computed the maximum severity S_recent
      # of any event detected within that window.
      #
      # The PERCENTILE RANK at each pixel is defined as:
      #
      #   pct_rank = 100 × (#historical events with S < S_recent) / #total events
      #
      # Interpretation:
      #   • pct_rank = 90  → the 2022-2025 drought was more severe than 90 % of
      #                       all historical droughts at that pixel.  In other
      #                       words it is approximately a 1-in-10 event locally.
      #   • pct_rank = 99  → exceeds 99 % of the record; roughly a 1-in-75-year
      #                       event (given a ~75-year record).
      #   • pct_rank = 0   → the 2022-2025 event is the mildest on record there.
      #   • NA (grey)      → no drought event was detected at that pixel for
      #                       either the historical period or 2022-2025 (e.g.
      #                       perennially wet pixels at high elevation).
      #
      # The colour scale uses plasma with direction = -1 so that high percentile
      # ranks appear DARK (visually prominent) and near-zero ranks appear LIGHT.
      # ──────────────────────────────────────────────────────────────────────
      # WHAT DOES "SEVERITY (I × D)" MEAN IN THIS MAP?
      # ──────────────────────────────────────────────────────────────────────
      # Severity follows Sheffield & Wood (2008, J. Climate):
      #   S = Ī × D
      # where Ī = mean intensity over the event duration and D = duration in months.
      # Intensity for a single month is: I_t = onset_threshold − index_value
      # (positive when the index is below the drought-onset threshold).
      #
      # This means a pixel can have high severity because it was:
      #   (a) extremely dry (high I) even for a short time, OR
      #   (b) persistently in drought (high D) even at moderate intensity.
      # The S map therefore captures the COMBINED space-time footprint of the
      # 2022-2025 event and complements the percentile-rank map.
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
        "2022-2025 Drought: Percentile rank vs all historical events"
      )
      p_sev_map <- make_base_map(
        df_sev, "severity",
        "Severity (I × D)",
        "inferno",
        "2022-2025 Drought: Pixel-level S&W severity"
      )
      
      p_rank_combo <- p_rank_map + p_sev_map +
        patchwork::plot_annotation(
          title    = "2022-2025 Drought — Pixel-Level Rank vs Historical Record",
          subtitle = sprintf(
            "%s  |  Percentile of 2022-2025 severity relative to all events at each pixel  |  %d–%d",
            toupper(PRIMARY_INDEX), HISTORICAL_START, HISTORICAL_END
          ),
          theme = ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(colour = "grey40", size = 10, hjust = 0.5)
          )
        )
      
      ggplot2::ggsave(
        file.path(RANKING_DIR, "Fig6_2022-2025_Pixel_Rank_Map.png"),
        p_rank_combo, width = 14, height = 7, dpi = 300
      )
      cat("  ✓ Saved: Fig6_2022-2025_Pixel_Rank_Map.png\n")
      rank_map_done <- TRUE
      
      cat(sprintf(
        "  Basin-wide: median percentile = %.1f%%  |  %.1f%% of pixels at ≥ 75th pct\n",
        median(rank_vec, na.rm = TRUE),
        100 * mean(rank_vec >= 75, na.rm = TRUE)
      ))
    }
  }
} else {
  cat(sprintf("  ⚠ Time-series RDS not found: %s\n  → Run w1_trend_test.R first, then re-run w5.\n",
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
    "Q4: How does 2022-2025 compare?"
  ),
  Answer = c(
    sprintf("Short-term: %d events (1 per %.1f yrs) | Long-term: %d events (1 per %.1f yrs)",
            sum(primary_events$duration_class == "D3-6 (Short-term)"),
            n_years_record / sum(primary_events$duration_class == "D3-6 (Short-term)"),
            sum(primary_events$duration_class == "D13+ (Long-term)"),
            n_years_record / sum(primary_events$duration_class == "D13+ (Long-term)")),
    "See Fig5 for pixel-level severity map (mean, minimum, duration during 2022-2025).",
    "Refer to w1_trend_test.R for Mann-Kendall trend analysis on frequency/intensity.",
    sprintf("Rank #%d of %d events (%.1fth percentile, ~1-in-%.1f year event)",
            most_severe_recent$rank_by_severity[1] %||% NA,
            nrow(primary_events),
            most_severe_recent$percentile_rank[1] %||% NA,
            most_severe_recent$return_period_years[1] %||% NA)
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
# scale 3 and SSPI uses scale 1. SPEI-ERA5 rows are included when its files exist.
# Output: Table_1a_ShortScale_Index_TopEvents.csv and workbook companion sheet.
#
# COLUMN DEFINITIONS (identical to Table 1b):
#   Rank              — composite severity rank within the index's own catalog
#   Start_date        — onset month of the event
#   End_date          — termination month
#   Start_year        — calendar year of onset (used as the event label)
#   Duration_months   — number of consecutive sub-threshold months
#   Duration_class    — D3-6 Short-term / D7-12 Medium-term / D13+ Long-term
#   Mean_intensity    — mean deficit below x₀ = −0.5 over event months [index units]
#   Severity          — Mean_intensity × Duration_months [index-units·months]
#   Return_period_yr  — Weibull empirical return period (n+1)/rank [years]
#   Percentile_rank   — 100 × (1 − rank/(n+1)) [%]
#   In_2022_2025      — logical flag; TRUE if event overlaps the 2022–2025 window
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 4e: Manuscript Table 1a\n")
cat("══════════════════════════════════════\n")

# ── Helper: extract top-N events for a set of index labels, add display columns ──
make_ms_table <- function(catalog, index_labels, n_top = 10L) {
  sub <- catalog %>%
    dplyr::filter(index_label %in% index_labels,
                  duration_class != "D1-2 (Sub-threshold)") %>%
    dplyr::arrange(index_label, rank_by_severity)
  
  # Pretty index labels, including SPEI-ERA5 and the native-scale snow indices.
  sub$Index <- vapply(sub$index_label, function(lbl) {
    scale_i <- as.integer(base::sub(".*_([0-9]+)$", "\\1", lbl))
    type_i <- base::sub("_[0-9]+$", "", lbl)
    switch(type_i,
           spi = sprintf("SPI-%d", scale_i),
           spei = sprintf("SPEI-PM-%d", scale_i),
           spei_era5 = sprintf("SPEI-ERA5-%d", scale_i),
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
      In_2022_2025     = is_recent
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
  
  # Highlight 2022-2025 rows in yellow
  for (rw in which(tbl$In_2022_2025)) {
    openxlsx::addStyle(wb, sheet_name,
                       openxlsx::createStyle(fgFill = "#FFEB9C"),
                       rows = rw + 1L, cols = seq_len(ncol(tbl)),
                       stack = TRUE)
  }
  
  openxlsx::setColWidths(wb, sheet_name,
                         cols = seq_len(ncol(tbl)),
                         widths = c(18, 6, 12, 12, 10, 16, 22, 16, 12, 16, 15, 13))
  invisible(wb)
}

tryCatch({
  
  # ── TABLE 1a: short-scale meteorological and native-scale snow indices ──────
  cat("  Building Table 1a (short-scale index comparison)...\n")
  
  table_1a_requested <- c(
    "spi_01", "spi_03",
    "spei_01", "spei_03",
    "spei_era5_01", "spei_era5_03",
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
                    Duration_months, Severity, In_2022_2025)
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
#   configured indices, using only each snow index's valid native scale.
#
# Encodings:
#   x-position       = event timing; horizontal segment spans start to end date
#   y-position       = accumulation scale (1, 3, 6, 12, 24, 36 months)
#   bubble area      = cumulative severity (S)
#   bubble colour    = historical severity rank class within the corresponding
#                      index × scale event catalog: Rank 1, Rank 2, Rank 3–5,
#                      Rank >5 (no rank number or severity text is printed
#                      inside the bubble)
#   text label       = severity value (no parentheses), placed deterministically
#                      beside each bubble; horizontal offsets scale with bubble
#                      size and labels are placed on the side with available room
#
# The plot uses the authoritative master_events object generated above. No
# event detection, severity calculation, or ranking is duplicated here.
# For each index × scale, the most severe event overlapping 2022–2025 is shown.
# This intentionally includes scale-1 events, even when their duration is D1–2,
# because scale 1 has a one-month minimum event-duration gate in this framework.
#
# Note on long-duration events: an event's bubble is anchored at the midpoint
# of its overlap with the 2021-12 to 2026-01 display window (not the midpoint
# of its full duration), and the view is zoomed with coord_cartesian() rather
# than scale limits. This ensures a long event (e.g., a 36-month-scale event
# that began well before 2022) still gets a bubble drawn on-screen instead of
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

# Include all supported formulations in the recent-event diagnostic. Snow
# indices appear only at their defined native scales (SWEI-03 and SSPI-01).
bubble_index_map <- c(
  "spi"       = "(a) SPI",
  "spei"      = "(b) SPEI-PM",
  "spei_era5" = "(c) SPEI-ERA5",
  "swei"      = "(d) SWEI",
  "sspi"      = "(e) SSPI"
)

bubble_scale_levels <- c(1L, 3L, 6L, 12L, 24L, 36L)

# Visible x-axis window for the figure. Defined once and reused for both the
# coordinate zoom (below) and for clamping each bubble's on-screen anchor
# date. This matters for long-duration events (e.g., a 36-month-scale event
# that began well before 2022): its true midpoint can fall before this
# window's lower bound, and previously the plot used scale_x_date(limits=...),
# which silently drops (sets to NA, then removes) any point whose x falls
# outside the given limits -- causing that bubble to vanish entirely even
# though the event does overlap 2022-2025. The fix below (a) clamps each
# bubble's anchor date into this window instead of relying on scale limits,
# and (b) uses coord_cartesian() to zoom the view without dropping data.
bubble_x_display_min <- as.Date("2021-12-01")
bubble_x_display_max <- as.Date("2026-01-01")
# Explicit yearly breaks that EXCLUDE the window end (2026-01-01). That break
# sat on the right panel edge, so its "2026" label overhung into the gap and
# collided with the "2022" label at the left edge of the next panel.
bubble_x_breaks <- seq(as.Date("2022-01-01"), bubble_x_display_max, by = "1 year")
bubble_x_breaks <- bubble_x_breaks[bubble_x_breaks < bubble_x_display_max]

bubble_rank_colours <- c(
  "Rank 1"   = "#C91F2D",
  "Rank 2"   = "#F47C20",
  "Rank 3–5" = "#F5B82E",
  "Rank >5"  = "#54A9B2"
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
# the most severe qualifying event overlapping 2022–2025.
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
    ongoing = FALSE,
    censored_at_record_end = FALSE
  )

# ------------------------------------------------------------------
# Figure-only fallback:
#
# If NO qualifying recent event exists for an index × scale under the
# If no qualifying recent completed event exists, repeat SW08 event detection with
# the same scale-dependent minimum duration over the COMPLETE 1950–2025 record.
# overlapping 2022–2025.
#
# This is essential for the 36-month scale. It preserves the full
# accumulated-index history and therefore does NOT recalculate a
# 36-month index from a cropped 2022–2025 series.
#
# The fallback is used ONLY for this figure. master_events, historical
# rankings, return periods, and the primary event catalogs are unchanged.
# ------------------------------------------------------------------

fallback_groups <- dplyr::bind_rows(lapply(INDICES_TO_ANALYZE, function(z) {
  data.frame(index_type = z$index, scale = as.integer(z$scales),
             stringsAsFactors = FALSE)
})) %>%
  dplyr::filter(index_type %in% names(bubble_index_map),
                scale %in% bubble_scale_levels) %>%
  dplyr::anti_join(
    completed_bubble_events %>%
      dplyr::distinct(index_type, scale),
    by = c("index_type", "scale")
  )

fallback_events_list <- list()

if (nrow(fallback_groups) > 0L) {
  
  for (i in seq_len(nrow(fallback_groups))) {
    
    idx_i <- fallback_groups$index_type[i]
    sc_i  <- as.integer(fallback_groups$scale[i])
    
    relaxed_catalog <- tryCatch(
      build_sw08_catalog(
        index_type     = idx_i,
        scale          = sc_i,
        onset_thr      = DROUGHT_ONSET_EVENT,
        min_duration   = get_event_min_duration(sc_i),
        start_year_rec = HISTORICAL_START,
        end_year_rec   = HISTORICAL_END
      ),
      error = function(e) {
        cat(sprintf(
          "  ⚠ Bubble fallback failed for %s-%02d: %s\n",
          idx_i, sc_i, e$message
        ))
        NULL
      }
    )
    
    if (is.null(relaxed_catalog) || nrow(relaxed_catalog) == 0L)
      next
    
    # Compute historical severity rank within the RELAXED figure-only
    # catalog before restricting it to 2022–2025.
    relaxed_catalog <- relaxed_catalog %>%
      dplyr::group_by(index_label) %>%
      dplyr::arrange(
        dplyr::desc(severity),
        start_date,
        end_date,
        .by_group = TRUE
      ) %>%
      dplyr::mutate(
        rank_by_severity = dplyr::row_number(),
        return_period_years =
          calc_return_period(dplyr::row_number(), n_years_record),
        percentile_rank =
          round(100 * (1 - rank_by_severity / dplyr::n()), 1)
      ) %>%
      dplyr::ungroup()
    
    # IMPORTANT: recent-period filtering occurs AFTER full-record
    # event detection.
    relaxed_recent <- relaxed_catalog %>%
      dplyr::filter(
        start_date <= RECENT_DROUGHT_END,
        end_date   >= RECENT_DROUGHT_START
      ) %>%
      dplyr::arrange(
        dplyr::desc(severity),
        start_date,
        end_date
      ) %>%
      dplyr::slice_head(n = 1L)
    
    if (nrow(relaxed_recent) > 0L) {
      relaxed_recent$ongoing <- FALSE
      relaxed_recent$censored_at_record_end <- FALSE
      relaxed_recent$min_duration_required <- sc_i
      relaxed_recent$figure_only_fallback <- TRUE
      
      fallback_events_list[[length(fallback_events_list) + 1L]] <-
        relaxed_recent
    }
  }
}

fallback_events <- if (length(fallback_events_list)) {
  dplyr::bind_rows(fallback_events_list)
} else {
  NULL
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
  
  # Rendered bubble diameter, expressed in the same size units used by
  # scale_size_continuous(range = c(5.0, 18.5)).
  sev_size_rng <- range(sqrt(bubble_events$severity), na.rm = TRUE)
  bubble_events$point_size_pts <- if (diff(sev_size_rng) < 1e-9) {
    rep(mean(c(5.0, 18.5)), nrow(bubble_events))
  } else {
    scales::rescale(
      sqrt(bubble_events$severity),
      to = c(5.0, 18.5),
      from = sev_size_rng
    )
  }
  
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
  
  # Severity is encoded through bubble AREA, not bubble radius. ggplot2's
  # size aesthetic is approximately radius/diameter-like, so sqrt(S) is used
  # to make the displayed area approximately proportional to cumulative severity.
  severity_legend_values <- c(2, 5, 10, 20, 40)
  severity_legend_sizes <- sqrt(severity_legend_values)
  
  # Small vertical headroom above the top row of bubbles.
  bubble_y_top <- length(bubble_scale_levels) + 1.15  # extra headroom for the severity label above the top row
  
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
        size = sqrt(severity),
        fill = rank_class
      ),
      shape = 21,
      colour = "#FFFFFF",
      stroke = 0.7,
      alpha = 0.97
    ) +
    # Severity value only (no rank, no parentheses). Labels are placed
    # deterministically beside the bubble rather than being allowed to repel
    # vertically into neighbouring accumulation-scale rows.
    # Severity value immediately above its own bubble; no repel is used.
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
      label.size   = 0,
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
      # Minor tick marks (one per month) alongside the major yearly ticks/labels.
      # Requires ggplot2 >= 3.5.0 for guide_axis(minor.ticks = TRUE).
      guide = ggplot2::guide_axis(minor.ticks = TRUE)
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
    ggplot2::scale_colour_manual(
      values = bubble_rank_colours,
      guide = "none",
      drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
      name = "Cumulative severity (S)",
      breaks = severity_legend_sizes,
      labels = as.character(severity_legend_values),
      range = c(5.0, 18.5),
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
      title = "Recent Drought Events in the Nechako River Basin (2022–2025)",
      subtitle = "Most severe recent event at each accumulation scale for each drought-index formulation",
      caption = paste(
        "Each bubble represents one drought event identified independently for the corresponding accumulation scale and index. ",
        "Bubble area is proportional to cumulative severity (S), whose value is labelled beside each bubble; completed-event bubble colour indicates ",
        "the historical severity rank class among 1950–2025 events for that index and scale. Horizontal segments indicate the event start and end dates. ",
        "Ongoing (censored) events are shown only when no completed qualifying event exists for the same index and scale; they remain unranked because they are right-censored at the record end."
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
      axis.minor.ticks.x.bottom = ggplot2::element_line(
        colour = "#888888",
        linewidth = 0.3
      ),
      axis.minor.ticks.length = ggplot2::rel(0.5),
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
  
  # Save both vector PDF and high-resolution PNG for manuscript use.
  ggplot2::ggsave(
    bubble_plot_path_pdf,
    plot = p_recent_bubbles,
    width = 12.0,
    height = 6.7,
    units = "in",
    device = "pdf"
  )
  ggplot2::ggsave(
    bubble_plot_path_png,
    plot = p_recent_bubbles,
    width = 12.0,
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
#   • manuscript SPEI-3 Kendall RP (Part 4c);
#   • multi-index comparison (Part 4d);
#   • all six manuscript scales for SPI and SPEI-PM, plus SWEI-03 and SSPI-01.
#
####################################################################################
cat("\n══════════════════════════════════════════════════════════════\n")
cat("PART 4g: Kendall RP — all manuscript scales\n")
cat("══════════════════════════════════════════════════════════════\n")

if (!exists("master_events") || !nrow(master_events))
  stop("Part 4g requires master_events from Parts 0–2.")

KENDALL_ALL_SCALES <- as.integer(MANUSCRIPT_SCALES)
# Kendall bootstrap covers SPI, SPEI-PM, and snow indices at their valid scales.
# SPEI-ERA5 remains excluded from Kendall bootstrap calculations by project
# runtime constraint, while its event catalog and empirical rankings are retained.
KENDALL_ALL_INDICES <- c(
  paste0("spei_", sprintf("%02d", KENDALL_ALL_SCALES)),
  paste0("spi_", sprintf("%02d", KENDALL_ALL_SCALES)),
  sprintf("swei_%02d", SWEI_SCALE),
  sprintf("sspi_%02d", SSPI_SCALE)
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
      grepl("^spei_era5_", kendall_rp_all_scales$index_label) ~ "SPEI-ERA5",
      grepl("^spei_", kendall_rp_all_scales$index_label) ~ "SPEI-PM",
      grepl("^swei_", kendall_rp_all_scales$index_label) ~ "SWEI",
      grepl("^sspi_", kendall_rp_all_scales$index_label) ~ "SSPI",
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
      "Table_KendallRP_SPEI_PM_2022_2025_AllScales.csv"
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
# w5_EXPORT_ADDON.R
#
# PURPOSE: Export data structures needed by w12_bayesian_ene_analysis.R (MS2).
#          Append this block to the END of 9event_ranking.R, after all
#          existing parts (1–4d) have run.
#
# REQUIRES (already in environment from w5 main body):
#   ranked_catalog       — full SW08 event catalog with severity/duration
#   kendall_rp_result    — output of the Part 4c copula analysis (list)
#   n_years_record       — integer, length of the record in years
#   RANKING_DIR          — character, output directory path
#   HISTORICAL_START/END — integers, 1950 / 2025
#
# OUTPUTS written to RANKING_DIR:
#   ms2_annual_severity_series.csv  — annual maximum severity per year (SPEI-3)
#   ms2_copula_fit.rds              — copula fit object + stationary T_K reference
####################################################################################
cat("\n══════════════════════════════════════════════════\n")
cat("  EXPORT ADDON: MS2 data export\n")
cat("══════════════════════════════════════════════════\n")

# ── 1.  Annual maximum severity series (SPEI-3 only) --------------------------
# For each calendar year, extract the maximum severity across all events
# whose start_date falls in that year.  Years with no qualifying event get NA.
# This annual series is the response variable for the Bayesian Gamma model.

TARGET_IDX_EXPORT <- "spei_03"

catalog_spei3 <- ranked_catalog %>%
  dplyr::filter(index_label == TARGET_IDX_EXPORT) %>%
  dplyr::mutate(
    start_year_int = as.integer(format(as.Date(start_date), "%Y")),
    end_year_int   = as.integer(format(as.Date(end_date),   "%Y"))
  )

# Build complete year spine
year_spine <- data.frame(year = HISTORICAL_START:HISTORICAL_END)

# Annual max severity: use start_year of each event (conventional anchor)
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
#   - target_event:   severity and duration of the 2022-2025 rank-1 SPEI-3 event

if (exists("kendall_rp_result") && !is.null(kendall_rp_result)) {
  
  krp <- kendall_rp_result
  
  # Extract the target event row from the catalog
  target_event_row <- catalog_spei3 %>%
    dplyr::arrange(dplyr::desc(severity)) %>%
    dplyr::slice_head(n = 1)
  
  ms2_copula_export <- list(
    # Stationary copula objects
    best_copula_fit   = krp$best_fit,       # copulaFit object
    best_copula_name  = krp$best_cop,       # character, e.g. "Frank"
    best_dur_fit      = krp$fits_dur[[krp$best_dur]],  # fitdistr object
    best_sev_fit      = krp$fits_sev[[krp$best_sev]],  # fitdistr object
    best_dur_name     = krp$best_dur,
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
    mu_ia_years       = n_years_record / nrow(
      catalog_spei3 %>% dplyr::filter(severity > 0)),
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
  cat("    Run Part 4c of 9event_ranking.R first.\n")
}

cat("\n╔══════════════════════════════════════╗\n")
cat("║  9event_ranking.R  COMPLETE        ║\n")
cat("╚══════════════════════════════════════╝\n\n")

cat("KEY FINDINGS FOR PRESENTATION:\n")
cat(sprintf("  • Total drought events (1950-2025): %d\n", nrow(primary_events)))
cat(sprintf("  • 2022-2025 drought rank: #%d of %d\n",
            most_severe_recent$rank_by_severity[1] %||% NA,
            nrow(primary_events)))
cat(sprintf("  • Percentile: %.1fth (more severe than %.1f%% of historical droughts)\n",
            most_severe_recent$percentile_rank[1] %||% NA,
            most_severe_recent$percentile_rank[1] %||% NA))
cat(sprintf("  • Return period: ~1-in-%.1f years\n",
            most_severe_recent$return_period_years[1] %||% NA))
cat(sprintf("  • Duration class: %s\n",
            most_severe_recent$duration_class[1] %||% "N/A"))

cat("\nOUTPUT FILES:\n")
cat("  SW08 ranked catalog : ", file.path(RANKING_DIR, "ranked_event_catalog_SW08.csv"), "\n")
cat("  SW08 RP table       : ", file.path(RANKING_DIR, "return_period_by_duration_class_SW08.csv"), "\n")
cat("  Table 1a (short-scale indices): ", file.path(RANKING_DIR, "Table_1a_ShortScale_Index_TopEvents.csv"), "\n")
cat("  Excel summary       : ", file.path(RANKING_DIR, "Drought_Ranking_Summary.xlsx"), "\n")
cat("  Figure 1            : ", file.path(RANKING_DIR, "Fig1_Top_Droughts_By_Severity.png"), "\n")
cat("  Figure 2            : ", file.path(RANKING_DIR, "Fig2_Drought_Severity_Distribution.png"), "\n")
cat("  Figure 3            : ", file.path(RANKING_DIR, "Fig3_Drought_Frequency_By_Duration.png"), "\n")
cat("  Figure 4            : ", file.path(RANKING_DIR, "Fig4_Presentation_Dashboard.png"), "\n")
cat("  Figure RP           : ", file.path(RANKING_DIR, "FigRP_Return_Period_By_Duration_Class.png"), "\n")
cat("  Figure RP+2022      : ", file.path(RANKING_DIR, "FigRP_Return_Period_By_Duration_Class_with_2022-2025.png"), "\n")
cat("  Fig3 MS (with CI)      : ", file.path(RANKING_DIR, "Fig3_MS_Drought_Severity_Ranking.pdf"), "\n")
cat("  Fig3 MS (no CI)        : ", file.path(RANKING_DIR, "Fig3_MS_Drought_Severity_Ranking_noCI.pdf"), "\n")
cat("  Fig3 Kendall (with CI) : ", file.path(RANKING_DIR, "Fig3_MS_Drought_Severity_Ranking_KendallRP.pdf"), "\n")
cat("  Fig3 Kendall (no CI)   : ", file.path(RANKING_DIR, "Fig3_MS_Drought_Severity_Ranking_KendallRP_noCI.pdf"), "\n")
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
if (file.exists(file.path(RANKING_DIR, "Fig5_2022-2025_Spatial_Severity.png")))
  cat("  Figure 5          : ", file.path(RANKING_DIR, "Fig5_2022-2025_Spatial_Severity.png"), "\n")
if (file.exists(file.path(RANKING_DIR, "Fig6_2022-2025_Pixel_Rank_Map.png")))
  cat("  Figure 6          : ", file.path(RANKING_DIR, "Fig6_2022-2025_Pixel_Rank_Map.png"), "\n")
if (file.exists(file.path(RANKING_DIR, "pixel_rank_2022-2025.nc")))
  cat("  Rank raster       : ", file.path(RANKING_DIR, "pixel_rank_2022-2025.nc"), "\n")
################################################################################
# REGIONAL EXTENSION — SW08 EVENT RANKING FOR 7 SUBBASINS + WHOLE BASIN
################################################################################
REGIONAL_EVENT_DIR <- file.path(TREND_DIR, "regional_event_ranking")
REGIONAL_INDEX_DIR <- file.path(TREND_DIR, "regional_timeseries")
dir.create(REGIONAL_EVENT_DIR, recursive=TRUE, showWarnings=FALSE)

regional_event_files <- list.files(REGIONAL_INDEX_DIR,
                                   pattern="_regional_averaged_by_month\\.csv$",
                                   full.names=TRUE, ignore.case=TRUE)
regional_event_catalogs <- list()

for (f in regional_event_files) {
  x <- data.table::fread(f)
  req <- c("region_id","region_name","date","value","index_type","scale")
  if (!all(req %in% names(x))) next
  idx <- unique(tolower(x$index_type))[1]
  if (!idx %in% c("spi", "spei", "spei_era5", "swei", "sspi")) next
  scs <- sort(unique(as.integer(x$scale)))
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
      ev <- sw08_event_loop(vals_z,dates_z,as.integer(format(dates_z,"%Y")),
                            onset_thr=DROUGHT_ONSET_EVENT,
                            min_duration=get_event_min_duration(sc))
      if (!nrow(ev)) next
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
      rank_by_severity=dplyr::row_number(),
      return_period_years=(n_years_record+1)/rank_by_severity,
      percentile_rank=round(100*(1-rank_by_severity/dplyr::n()),1),
      is_recent=start_date <= as.Date("2025-12-31") & end_date >= as.Date("2022-01-01")
    ) |>
    dplyr::ungroup()
  
  regional_rp <- regional_events |>
    dplyr::filter(duration_class != "D1-2 (Sub-threshold)") |>
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
  ggplot2::ggsave(file.path(REGIONAL_EVENT_DIR,"Fig_SW08_duration_severity_regional.pdf"),p,width=18,height=12)
  cat(sprintf("✓ Regional SW08 catalog written for %d regions.\n",uniqueN(regional_events$region_id)))
}

################################################################################
# PART 4h: REGIONAL HISTORICAL-RANKING + KENDALL-RP FIGURES
#
# Produces 8 separate figures:
#   1 whole Nechako basin + 7 supplied subbasins.
#
# Figure A (one per region): recent-event historical ranking profile
#   • SPI, SPEI-PM, and SPEI-ERA5 at manuscript scales; SWEI-03 and SSPI-01 at native scales
#   • available accumulation scales: 1, 3, 6, 12, 24, 36 months
#   • horizontal segment = 2022–2025 overlap of the selected recent event
#   • bubble area = event severity
#   • colour = historical severity-rank class within region × index × scale
#   • label = historical rank and severity
#
# Figure B (one per region): Kendall duration–severity return-period profile
#   • SPI, SPEI-PM, SWEI, and SSPI where regional event data exist
#   • available accumulation scales; snow indices use native scales only
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

REGIONAL_KENDALL_BOOT <- 500L
REGIONAL_KENDALL_ALPHA <- 0.05
REGIONAL_FIG_START <- RECENT_DROUGHT_START
REGIONAL_FIG_END   <- RECENT_DROUGHT_END
REGIONAL_FIG_SCALES <- as.integer(MANUSCRIPT_SCALES)
# Regional historical ranking figures include all five event-catalog indices.
# SPEI-ERA5 is excluded only from Kendall bootstrap figures under the project's
# runtime constraint; its empirical regional event rankings remain available.
REGIONAL_FIG_INDEXES <- c("spi", "spei", "spei_era5", "swei", "sspi")
REGIONAL_KENDALL_INDEXES <- c("spi", "spei", "swei", "sspi")
REGIONAL_FIG_LABELS <- c(
  spi       = "(a) SPI",
  spei      = "(b) SPEI-PM",
  spei_era5 = "(c) SPEI-ERA5",
  swei      = "(d) SWEI",
  sspi      = "(e) SSPI"
)

REGIONAL_RANK_COLOURS <- c(
  "Rank 1"   = "#C91F2D",
  "Rank 2"   = "#F47C20",
  "Rank 3-5" = "#F5B82E",
  "Rank >5"  = "#54A9B2"
)

if (!exists("regional_events") || !nrow(regional_events)) {
  cat("  WARNING: regional_events is empty; Part 4h skipped.\n")
} else {
  # --------------------------------------------------------------------------
  # 4h.1 — restrict to configured observed formulations and manuscript scales;
  #        snow indices appear only at their valid native scales.
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
    region_catalog_file_h <- file.path(
      TREND_DIR, "regional", "region_catalog.csv"
    )
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
      # FIGURE A — historical ranking of the 2022–2025 primary event
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
              "2022–2025 primary event at each accumulation scale; bubble size = severity; ",
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
            p_rank, width = 13, height = 5.8, units = "in", device = "pdf"
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
          kr_j <- run_kendall_rp(
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
          dplyr::mutate(
            index_type = factor(
              as.character(index_type),
              levels = REGIONAL_KENDALL_INDEXES,
              labels = unname(REGIONAL_FIG_LABELS[REGIONAL_KENDALL_INDEXES])
            )
          )
        
        regional_kendall_rows[[length(regional_kendall_rows) + 1L]] <- reg_kendall_df
        
        rp_plot_df <- reg_kendall_df %>%
          dplyr::mutate(
            index_panel = factor(
              as.character(index_type),
              levels = unname(REGIONAL_FIG_LABELS[REGIONAL_KENDALL_INDEXES])
            ),
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
          ggplot2::scale_y_discrete(
            labels = REGIONAL_FIG_SCALES,
            drop = FALSE
          ) +
          ggplot2::labs(
            title = paste0(rname, " — Kendall duration–severity return periods"),
            subtitle = paste0(
              "Primary 2022–2025 event at each scale; point = T_K; interval = 95% bootstrap CI; ",
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
        valid_cop <- reg_kendall_df %>%
          dplyr::filter(is.finite(Kendall_return_period_years)) %>%
          dplyr::transmute(
            idx = as.character(index_type),
            cop = as.character(best_copula),
            tau = as.numeric(tau_empirical)
          )
        if (nrow(valid_cop)) {
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
            p_rp, width = 13, height = 6.4, units = "in", dpi = 320
          )
          ggplot2::ggsave(
            file.path(REGIONAL_EVENT_DIR, paste0(rp_file_stub, ".pdf")),
            p_rp, width = 13, height = 6.4, units = "in", device = "pdf"
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
        file.path(REGIONAL_EVENT_DIR, "Table_Regional_Historical_Ranking_2022_2025.csv"),
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

################################################################################
# END OF 9event_ranking.R
################################################################################