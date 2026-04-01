# ==============================================================================
# section43_manuscript_blanks.R
#
# Extracts EVERY statistic needed for manuscript Section 4.3 blanks.
# Self-contained: re-implements the SW08 event loop from w5_event_ranking.R
# so it can be run after w5 has been executed (reads existing CSVs when
# available, otherwise rebuilds from the NetCDF seasonal files).
#
# PREREQUISITES:
#   3SPEI_ERALand.R   → spei_results_seasonal/
#   1SPI_ERALand.R    → spi_results_seasonal/
#   w1_trend_test.R   → (optional) pre-computed basin-average CSVs
#   w5_event_ranking.R → (optional) ranked_event_catalog_SW08.csv
#
# OUTPUT: section43_blanks_summary.csv  +  Table1_top_events_SPEI3.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(lubridate)
})

# ── USER SETTINGS ─────────────────────────────────────────────────────────────
WD_PATH          <- "D:/Nechako_Drought/Nechako/"
SPEI_SEAS_DIR    <- file.path(WD_PATH, "spei_results_seasonal")
SPI_SEAS_DIR     <- file.path(WD_PATH, "spi_results_seasonal")
RANKING_DIR      <- file.path(WD_PATH, "temporal_drought", "event_ranking")

HISTORICAL_START <- 1950L
HISTORICAL_END   <- 2025L
FOCUS_START      <- 2022L
FOCUS_END        <- 2025L
ONSET_THR        <- -0.5   # x₀ Yevjevich / Lloyd-Hughes threshold
MIN_DUR_SPEI3    <- 3L     # minimum event duration for scale-3 indices
n_years_record   <- HISTORICAL_END - HISTORICAL_START + 1L
MONTH_NAMES      <- c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec")

setwd(WD_PATH)

cat("\n============================================================\n")
cat(" SECTION 4.3 MANUSCRIPT BLANKS — SPEI-3 & SPI-3 SW08 EVENTS\n")
cat(sprintf(" Threshold x\u2080 = %.1f | Min duration = %d months\n",
            ONSET_THR, MIN_DUR_SPEI3))
cat("============================================================\n\n")

# ==============================================================================
# HELPER FUNCTIONS (mirror w5_event_ranking.R logic exactly)
# ==============================================================================

classify_duration <- function(d) {
  dplyr::case_when(
    d >= 3L  & d <= 6L  ~ "D3-6 (Short-term)",
    d >= 7L  & d <= 12L ~ "D7-12 (Medium-term)",
    d >= 13L            ~ "D12+ (Long-term)",
    TRUE                ~ "D1-2 (Sub-threshold)"
  )
}

#' SW08 single-threshold event loop
sw08_event_loop <- function(vals, date_seq, yr_seq,
                            onset_thr    = ONSET_THR,
                            min_duration = MIN_DUR_SPEI3) {
  min_duration <- as.integer(min_duration)
  in_d  <- FALSE;  s_idx <- NA_integer_;  events <- list()
  
  close_event <- function(seg, s, e) {
    dur    <- length(seg)
    mean_i <- mean(onset_thr - seg, na.rm = TRUE)
    if (dur >= min_duration)
      events[[length(events) + 1L]] <<- data.frame(
        start_date      = date_seq[s],
        end_date        = date_seq[e],
        start_year      = yr_seq[s],
        duration_months = as.integer(dur),
        duration_class  = classify_duration(dur),
        mean_int        = mean_i,
        severity        = mean_i * dur,
        stringsAsFactors = FALSE)
  }
  
  for (j in seq_along(vals)) {
    v <- vals[j]; if (is.na(v)) next
    if (!in_d && v < onset_thr)   { in_d <- TRUE;  s_idx <- j }
    else if (in_d && v >= onset_thr) {
      close_event(vals[s_idx:(j-1L)], s_idx, j-1L)
      in_d <- FALSE; s_idx <- NA_integer_
    }
  }
  if (in_d && !is.na(s_idx))
    close_event(vals[s_idx:length(vals)], s_idx, length(vals))
  
  if (!length(events))
    return(data.frame(start_date=as.Date(character(0)),
                      end_date=as.Date(character(0)),
                      start_year=integer(0), duration_months=integer(0),
                      duration_class=character(0), mean_int=numeric(0),
                      severity=numeric(0), stringsAsFactors=FALSE))
  do.call(rbind, events)
}

#' Load basin-average time series from pre-computed wide CSV
#' (written by w1_trend_test.R as {index}_{scale:02d}_basin_averaged_by_month.csv)
load_basin_avg_csv_local <- function(seas_dir, index_type, scale) {
  f <- file.path(seas_dir, sprintf("%s_%02d_basin_averaged_by_month.csv",
                                   tolower(index_type), as.integer(scale)))
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, stringsAsFactors = FALSE)
  names(df)[1] <- "Year"
  df$Year <- as.integer(df$Year)
  mon_cols <- intersect(MONTH_NAMES, names(df))
  if (!length(mon_cols)) return(NULL)
  long <- do.call(rbind, lapply(seq_along(mon_cols), function(mi) {
    data.frame(
      date  = as.Date(paste(df$Year, mi, "01", sep = "-")),
      value = as.numeric(df[[mon_cols[mi]]]),
      stringsAsFactors = FALSE)
  }))
  long <- long[order(long$date), ]
  long[!is.na(long$value), ]
}

#' Build basin-average time series from raw seasonal NetCDF files
build_basin_avg_from_nc <- function(seas_dir, index_type, scale) {
  pat   <- sprintf("^%s_%02d_month\\d{2}_[A-Za-z]+\\.nc$",
                   tolower(index_type), scale)
  files <- list.files(seas_dir, pattern = pat, full.names = TRUE)
  if (!length(files)) return(NULL)
  files <- files[order(as.integer(regmatches(basename(files),
                                             regexpr("(?<=month)\\d{2}", basename(files), perl=TRUE))))]
  
  stacks <- lapply(files, function(f) tryCatch(terra::rast(f), error=function(e) NULL))
  stacks <- Filter(Negate(is.null), stacks)
  if (!length(stacks)) return(NULL)
  
  r_tmpl   <- stacks[[1]][[1]]
  n_yrs    <- terra::nlyr(stacks[[1]])
  ts_mat   <- matrix(NA_real_, nrow=terra::ncell(r_tmpl), ncol=n_yrs*12L)
  
  # Recover year vector from first file's time axis
  t_raw  <- tryCatch(as.Date(terra::time(stacks[[1]])), error=function(e) NULL)
  if (is.null(t_raw) || length(t_raw) != n_yrs || all(is.na(t_raw)))
    t_raw <- seq(as.Date(sprintf("%d-01-01", HISTORICAL_START)),
                 by="year", length.out=n_yrs)
  years <- as.integer(format(t_raw, "%Y"))
  
  for (m in seq_along(stacks)) {
    n_yrs_m <- min(terra::nlyr(stacks[[m]]), n_yrs)
    for (y in seq_len(n_yrs_m))
      ts_mat[, (y-1L)*12L+m] <- as.numeric(terra::values(stacks[[m]][[y]]))
  }
  
  dates_ts <- as.Date(sprintf("%04d-%02d-01",
                              rep(years, each=12L), rep(1:12, times=n_yrs)))
  yr_ts    <- as.integer(format(dates_ts, "%Y"))
  keep     <- yr_ts >= HISTORICAL_START & yr_ts <= HISTORICAL_END
  ts_mat   <- ts_mat[, keep, drop=FALSE]
  dates_ts <- dates_ts[keep]
  
  valid_pix  <- rowSums(!is.na(ts_mat)) > 0
  ts_valid   <- ts_mat[valid_pix, , drop=FALSE]
  cell_areas <- as.numeric(terra::values(terra::cellSize(r_tmpl, unit="m"),
                                         na.rm=FALSE))
  aw  <- cell_areas[valid_pix]
  bad <- is.na(aw) | !is.finite(aw) | aw <= 0
  if (any(bad)) aw[bad] <- median(aw[!bad], na.rm=TRUE)
  
  basin_avg <- vapply(seq_len(ncol(ts_valid)), function(t) {
    v <- ts_valid[, t]; ok <- !is.na(v) & is.finite(v)
    if (!any(ok)) return(NA_real_)
    sum(v[ok]*aw[ok]) / sum(aw[ok])
  }, numeric(1L))
  
  data.frame(date=dates_ts, value=basin_avg, stringsAsFactors=FALSE)
}

#' Get basin-average time series: try pre-computed CSV first, then NC rebuild
get_ts <- function(seas_dir, index_type, scale) {
  ts <- load_basin_avg_csv_local(seas_dir, index_type, scale)
  if (!is.null(ts)) {
    cat(sprintf("  [%s-%d] Loaded from pre-computed CSV (%d months)\n",
                toupper(index_type), scale, nrow(ts)))
    return(ts)
  }
  cat(sprintf("  [%s-%d] CSV not found — rebuilding from NetCDF files...\n",
              toupper(index_type), scale))
  ts <- build_basin_avg_from_nc(seas_dir, index_type, scale)
  if (is.null(ts)) stop(sprintf("%s-%d: no data source found.", toupper(index_type), scale))
  cat(sprintf("  [%s-%d] Rebuilt from NC: %d months\n",
              toupper(index_type), scale, nrow(ts)))
  ts
}

# ==============================================================================
# STEP 1: BUILD SW08 EVENT CATALOGS — SPEI-3 and SPI-3
# ==============================================================================
cat("STEP 1: Building SW08 event catalogs\n")
cat("─────────────────────────────────────\n")

build_catalog <- function(seas_dir, index_type, scale) {
  ts   <- get_ts(seas_dir, index_type, scale)
  ts   <- ts[!is.na(ts$value), ]
  yrs  <- as.integer(format(ts$date, "%Y"))
  cat(sprintf("  [%s-%d] Running SW08 event loop (threshold=%.1f, min_dur=%d m)...\n",
              toupper(index_type), scale, ONSET_THR, MIN_DUR_SPEI3))
  ev   <- sw08_event_loop(ts$value, ts$date, yrs,
                          onset_thr=ONSET_THR, min_duration=MIN_DUR_SPEI3)
  ev$index_label  <- sprintf("%s_%02d", tolower(index_type), scale)
  ev$is_recent    <- ev$start_date <= as.Date(sprintf("%d-12-31", FOCUS_END)) &
    ev$end_date   >= as.Date(sprintf("%d-01-01", FOCUS_START))
  # Rank by descending severity
  ev <- ev[order(-ev$severity), ]
  ev$rank_by_severity  <- seq_len(nrow(ev))
  ev$return_period_yrs <- (n_years_record + 1L) / ev$rank_by_severity
  ev$percentile_rank   <- round(100 * (1 - ev$rank_by_severity / nrow(ev)), 1)
  cat(sprintf("  [%s-%d] %d qualifying events detected\n",
              toupper(index_type), scale, nrow(ev)))
  # Save per-index CSV
  out_f <- file.path(RANKING_DIR,
                     sprintf("%s-%02d_SW08_x05_event_catalog.csv",
                             toupper(index_type), scale))
  dir.create(RANKING_DIR, showWarnings=FALSE, recursive=TRUE)
  write.csv(ev, out_f, row.names=FALSE)
  ev
}

spei3_events <- build_catalog(SPEI_SEAS_DIR, "spei", 3L)
spi3_events  <- build_catalog(SPI_SEAS_DIR,  "spi",  3L)

# ==============================================================================
# STEP 2: EXTRACT ALL SECTION 4.3 BLANK VALUES — SPEI-3
# ==============================================================================
cat("\n\nSTEP 2: Extracting SPEI-3 statistics\n")
cat("─────────────────────────────────────\n")

ev <- spei3_events

# ── Overall catalog ────────────────────────────────────────────────────────────
ev_valid <- ev[ev$duration_class != "D1-2 (Sub-threshold)", ]

n_total   <- nrow(ev_valid)
dur_min   <- min(ev_valid$duration_months)
dur_max   <- max(ev_valid$duration_months)
sev_min   <- round(min(ev_valid$severity), 3)
sev_max   <- round(max(ev_valid$severity), 3)
n_d36     <- sum(ev_valid$duration_class == "D3-6 (Short-term)")
n_d712    <- sum(ev_valid$duration_class == "D7-12 (Medium-term)")
n_d12p    <- sum(ev_valid$duration_class == "D12+ (Long-term)")

cat(sprintf("\nTotal qualifying events (D3+): %d\n", n_total))
cat(sprintf("Duration range: %d – %d months\n",     dur_min, dur_max))
cat(sprintf("Severity range: %.3f – %.3f index-units·months\n", sev_min, sev_max))
cat(sprintf("D3-6  (short-term):  %d events\n", n_d36))
cat(sprintf("D7-12 (medium-term): %d events\n", n_d712))
cat(sprintf("D12+  (long-term):   %d events\n", n_d12p))

# ── 2022–2025 event ────────────────────────────────────────────────────────────
cat("\n── 2022–2025 event ──\n")
recent <- ev_valid[ev_valid$is_recent, ]

if (nrow(recent) == 0) {
  # If split into sub-events, take the single largest spanning the focus window
  recent_all <- ev[ev$is_recent, ]
  if (nrow(recent_all) == 0)
    stop("No event overlapping 2022-2025 found — check data extends through 2025.")
  recent <- recent_all[which.max(recent_all$severity), , drop=FALSE]
  cat("  ⚠ 2022–2025 event classified as D1-2 (sub-threshold) — using largest overlap.\n")
}

# Take the longest/most severe event that covers the focus window
rec_evt <- recent[which.max(recent$severity), , drop=FALSE]

dur_2022      <- rec_evt$duration_months
start_2022    <- format(as.Date(rec_evt$start_date), "%B %Y")
end_2022      <- format(as.Date(rec_evt$end_date),   "%B %Y")
int_2022      <- round(rec_evt$mean_int, 3)
sev_2022      <- round(rec_evt$severity, 2)
rank_sev_2022 <- rec_evt$rank_by_severity
pct_2022      <- rec_evt$percentile_rank

cat(sprintf("  Duration:  %d months\n", dur_2022))
cat(sprintf("  Period:    %s to %s\n", start_2022, end_2022))
cat(sprintf("  Mean intensity I: %.3f index-units below x\u2080\n", int_2022))
cat(sprintf("  Composite severity S: %.2f index-units·months\n", sev_2022))
cat(sprintf("  Severity rank: #%d of %d (%.1fth percentile)\n",
            rank_sev_2022, n_total, pct_2022))

# Rank by intensity as well
ev_valid2 <- ev_valid[order(-ev_valid$mean_int), ]
ev_valid2$rank_by_intensity <- seq_len(nrow(ev_valid2))
rank_int_2022 <- ev_valid2$rank_by_intensity[
  match(paste(rec_evt$start_date, rec_evt$end_date),
        paste(ev_valid2$start_date, ev_valid2$end_date))]
cat(sprintf("  Intensity rank: #%d\n", rank_int_2022))

# ── Second-ranked historical event (non-2022-2025) ──────────────────────────
cat("\n── Previous highest-ranked event ──\n")
historical   <- ev_valid[!ev_valid$is_recent, ]
historical   <- historical[order(-historical$severity), ]
prev_top     <- historical[1, , drop=FALSE]

prev_name    <- format(as.Date(prev_top$start_date), "%Y")   # year label
prev_period  <- sprintf("%s to %s",
                        format(as.Date(prev_top$start_date), "%B %Y"),
                        format(as.Date(prev_top$end_date),   "%B %Y"))
prev_sev     <- round(prev_top$severity, 2)
sev_ratio    <- round(sev_2022 / prev_sev, 2)

cat(sprintf("  Event period: %s\n", prev_period))
cat(sprintf("  Year label:   %s drought\n", prev_name))
cat(sprintf("  Severity:     %.2f index-units·months\n", prev_sev))
cat(sprintf("  Ratio 2022-2025 / prev: %.2f×\n", sev_ratio))

# ── Year-level composite severity scores ──────────────────────────────────────
cat("\n── Year-level composite severity scores (2022–2025) ──\n")
# Use monthly time series directly for per-year stats in focus window
ts_spei3 <- get_ts(SPEI_SEAS_DIR, "spei", 3L)
ts_spei3 <- ts_spei3[!is.na(ts_spei3$value), ]
ts_spei3$year   <- as.integer(format(ts_spei3$date, "%Y"))
ts_spei3$month  <- as.integer(format(ts_spei3$date, "%m"))

# Annual stats: drought months (< x₀) and mean severity (|value - x₀| for dry months)
annual_stats <- ts_spei3 %>%
  dplyr::filter(year >= FOCUS_START, year <= FOCUS_END) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    n_drought_months = sum(value < ONSET_THR, na.rm=TRUE),
    mean_severity_idx = round(
      mean(ONSET_THR - value[value < ONSET_THR], na.rm=TRUE), 3),
    composite_score = n_drought_months *
      mean(ONSET_THR - value[value < ONSET_THR], na.rm=TRUE),
    .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(composite_score))

cat("\nFocus period year-level stats (SPEI-3, drought months below x₀=-0.5):\n")
print(as.data.frame(annual_stats))

# Rank all calendar years in full record
all_annual <- ts_spei3 %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    n_drought_months = sum(value < ONSET_THR, na.rm=TRUE),
    mean_sev = mean(ONSET_THR - value[value < ONSET_THR], na.rm=TRUE),
    composite_score = n_drought_months * mean(ONSET_THR - value[value < ONSET_THR], na.rm=TRUE),
    .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(composite_score)) %>%
  dplyr::mutate(annual_rank = dplyr::row_number())

cat("\nTop-10 calendar years by annual composite score (all years 1950–2025):\n")
print(as.data.frame(head(all_annual, 10)))

yr_2023 <- all_annual[all_annual$year == 2023, , drop=FALSE]
yr_2024 <- all_annual[all_annual$year == 2024, , drop=FALSE]
yr_1980 <- all_annual[all_annual$year == 1980, , drop=FALSE]

cat(sprintf("\n2023: rank #%d | %d drought months | mean severity = %.3f\n",
            yr_2023$annual_rank, yr_2023$n_drought_months, yr_2023$mean_sev))
cat(sprintf("2024: rank #%d | %d drought months | mean severity = %.3f\n",
            yr_2024$annual_rank, yr_2024$n_drought_months, yr_2024$mean_sev))
if (nrow(yr_1980) > 0)
  cat(sprintf("1980: rank #%d | %d drought months | mean severity = %.3f\n",
              yr_1980$annual_rank, yr_1980$n_drought_months, yr_1980$mean_sev))

# ==============================================================================
# STEP 3: SPI-3 COMPARISON FOR THE SAME 2022–2025 EVENT
# ==============================================================================
cat("\n\nSTEP 3: SPI-3 comparison for 2022–2025 event\n")
cat("──────────────────────────────────────────────\n")

spi_ev       <- spi3_events
spi_ev_valid <- spi_ev[spi_ev$duration_class != "D1-2 (Sub-threshold)", ]
spi_recent   <- spi_ev[spi_ev$is_recent, ]

if (nrow(spi_recent) == 0) {
  spi_recent_all <- spi_ev[spi_ev$is_recent, ]
  spi_recent <- if (nrow(spi_recent_all) > 0)
    spi_recent_all[which.max(spi_recent_all$severity), , drop=FALSE]
  else NULL
}

if (!is.null(spi_recent) && nrow(spi_recent) > 0) {
  spi_rec_evt    <- spi_recent[which.max(spi_recent$severity), , drop=FALSE]
  spi_sev_2022   <- round(spi_rec_evt$severity, 2)
  spi_rank_2022  <- spi_rec_evt$rank_by_severity
  spi_pct_2022   <- spi_rec_evt$percentile_rank
  spi_n_total    <- nrow(spi_ev_valid)
  
  cat(sprintf("SPI-3 2022–2025 event:\n"))
  cat(sprintf("  Composite severity S: %.2f index-units·months\n", spi_sev_2022))
  cat(sprintf("  Severity rank: #%d of %d (%.1fth percentile)\n",
              spi_rank_2022, spi_n_total, spi_pct_2022))
  cat(sprintf("  SPEI-3 rank (#%d) vs SPI-3 rank (#%d): ",
              rank_sev_2022, spi_rank_2022))
  if (rank_sev_2022 < spi_rank_2022)
    cat("SPEI ranks MORE severe → thermodynamic amplification confirmed.\n")
  else
    cat("SPI ranks equally or more severe.\n")
} else {
  cat("  ⚠ No SPI-3 event detected overlapping 2022–2025.\n")
  spi_sev_2022  <- NA; spi_rank_2022 <- NA
}

# ==============================================================================
# STEP 4: GENERATE TABLE 1 — TOP-RANKED SPEI-3 EVENTS
# ==============================================================================
cat("\n\nSTEP 4: Generating Table 1 (top-ranked SPEI-3 events)\n")
cat("───────────────────────────────────────────────────────\n")

table1 <- ev_valid %>%
  dplyr::filter(rank_by_severity <= 15 | is_recent) %>%
  dplyr::arrange(rank_by_severity) %>%
  dplyr::mutate(
    Event_period     = sprintf("%s – %s",
                               format(as.Date(start_date), "%b %Y"),
                               format(as.Date(end_date),   "%b %Y")),
    Rank             = rank_by_severity,
    Duration_months  = duration_months,
    Duration_class   = duration_class,
    Mean_intensity_I = round(mean_int, 3),
    Severity_S       = round(severity, 2),
    Peak_SPEI        = NA_real_,   # filled separately if monthly TS available
    `Is_2022-2025`   = is_recent
  ) %>%
  dplyr::select(Rank, Event_period, Duration_months, Duration_class,
                Mean_intensity_I, Severity_S, `Is_2022-2025`)

cat("\nTable 1 — Top-15 SPEI-3 drought events (1950–2025):\n")
print(as.data.frame(table1), row.names=FALSE)

t1_path <- file.path(RANKING_DIR, "Table1_top_events_SPEI3_x05.csv")
write.csv(table1, t1_path, row.names=FALSE)
cat(sprintf("\n✓ Table 1 saved: %s\n", t1_path))

# ==============================================================================
# STEP 5: ASSEMBLE COMPLETE BLANKS SUMMARY
# ==============================================================================
cat("\n\nSTEP 5: Complete blank-value summary for Section 4.3\n")
cat("═══════════════════════════════════════════════════════\n\n")

blanks <- data.frame(
  Blank_position = c(
    "P1: total SPEI-3 drought events",
    "P1: min event duration (months)",
    "P1: max event duration (months)",
    "P1: min composite severity S",
    "P1: max composite severity S",
    "P1: n D3-6 short-term events",
    "P1: n D7-12 medium-term events",
    "P1: n D12+ long-term events",
    "P2: 2022-2025 duration (months)",
    "P2: 2022-2025 start month/year",
    "P2: 2022-2025 end month/year",
    "P2: 2022-2025 mean intensity I",
    "P2: 2022-2025 composite severity S",
    "P2: severity rank",
    "P2: intensity rank",
    "P2: prev top event — year label",
    "P2: prev top event — period",
    "P2: prev top event — severity S",
    "P2: 2022-2025 / prev ratio",
    "P2: 2023 — n drought months",
    "P2: 2023 — mean severity index",
    "P2: 2023 — annual rank",
    "P2: 2024 — n drought months",
    "P2: 2024 — mean severity index",
    "P2: 2024 — annual rank",
    "P3: SPI-3 2022-2025 severity S",
    "P3: SPI-3 2022-2025 severity rank"
  ),
  Value = c(
    n_total, dur_min, dur_max, sev_min, sev_max,
    n_d36, n_d712, n_d12p,
    dur_2022, start_2022, end_2022,
    int_2022, sev_2022, rank_sev_2022, rank_int_2022 %||% "NA",
    prev_name, prev_period, prev_sev, sev_ratio,
    yr_2023$n_drought_months, round(yr_2023$mean_sev, 3), yr_2023$annual_rank,
    yr_2024$n_drought_months, round(yr_2024$mean_sev, 3), yr_2024$annual_rank,
    spi_sev_2022 %||% "NA", spi_rank_2022 %||% "NA"
  ),
  Notes = c(
    "D3+ only; sub-threshold D1-2 excluded",
    "All D3+ events", "All D3+ events",
    "All D3+ events; I×D", "All D3+ events; I×D",
    "Sheffield & Wood (2008) class", "Sheffield & Wood (2008) class",
    "Sheffield & Wood (2008) class",
    "Contiguous months below x₀=-0.5",
    "First month below threshold", "Last month below threshold",
    "mean(x₀ - SPEI_t) over event months",
    "I × D", "Rank 1 = most severe",
    "Rank 1 = most intense",
    "Year of previous top historical event",
    "Full date range of previous top event",
    "I × D of previous top event",
    "sev_2022 / prev_sev",
    "Months with SPEI-3 < -0.5 in 2023",
    "mean(x₀ - SPEI_t) for drought months in 2023",
    "1 = highest annual composite score",
    "Months with SPEI-3 < -0.5 in 2024",
    "mean(x₀ - SPEI_t) for drought months in 2024",
    "1 = highest annual composite score",
    "SPI-3 I × D for 2022-2025 event",
    "SPI-3 rank; higher # → less severe than SPEI"
  ),
  stringsAsFactors = FALSE
)

`%||%` <- function(a, b) if (is.null(a) || length(a)==0 || is.na(a)) b else a

print(blanks, right=FALSE)

out_path <- "section43_blanks_summary.csv"
write.csv(blanks, out_path, row.names=FALSE)
cat(sprintf("\n✓ All blank values saved to: %s\n", out_path))
cat(sprintf("✓ Table 1 saved to:          %s\n", t1_path))
cat("============================================================\n")