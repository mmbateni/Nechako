##############################################################################
# 11Nechako_drought_manuscript1.R
#
# MERGED SCRIPT — combines:
#   1. table3.R                        → Table 3: Thermodynamic Fraction Bias plot
#   2. manuscript_statistics_extraction.R (section43_manuscript_blanks.R)
#                                      → Section 4.3 SW08 event catalog & blanks
#   3. w_extract_ms1_blanks.R         → MS1 blank-filling extractor (all sections)
#   4. Fig2_MS_SPEI123_SPI123.R (FigS_Seasonality_P_PET.R)
#                                      → Seasonal P vs PET two-panel figure
#
# RUN ORDER: Execute all four sections in sequence. Sections 2 and 3 share
# helper functions defined once in Section 2; Section 3 depends on outputs
# produced by upstream scripts (w5, w11, w9, w7). Section 4 depends on
# w2_basin_timeseries.R and 2preq_PET_ERALand.R.
#
# PREREQUISITES (must be run before this script):
#   1SPI_ERALand.R           → spi_results_seasonal/
#   3SPEI_ERALand.R          → spei_results_seasonal/, spei_results_seasonal_thw/
#   w1_trend_test.R          → (optional) pre-computed basin-average CSVs
#   w2_basin_timeseries.R    → temporal_drought/basin_averaged_timeseries/
#   2preq_PET_ERALand.R      → monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv
#   w5_event_ranking.R       → ranked_event_catalog_SW08.csv
#   w7_teleconnection_prep.R → PDO / Niño-3.4 CSVs
#   w9_atmospheric_diagnostics.R / w10a → composite_ms_blanks.csv
#   w11_dynamic_thermodynamic_decomp.R  → decomp_results/
##############################################################################


# ==============================================================================
# SECTION 1 — TABLE 3: Thermodynamic Fraction Sensitivity to Reanalysis Bias
# ==============================================================================
# SOURCE: table3.R
# PURPOSE: Plot original vs bias-corrected thermodynamic fraction (Fthm) by
#          month and SPEI scale (1, 2, 3) for active 2022–2025 drought months.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

cat("\n============================================================\n")
cat("  SECTION 1: TABLE 3 — Thermodynamic Fraction Bias Plot\n")
cat("============================================================\n\n")

# ── 1a. Prepare the data ──────────────────────────────────────────────────────
df_thm <- data.frame(
  Month = rep(month.abb[1:12], 3),
  Month_num = rep(1:12, 3),
  Scale = rep(c("SPEI-1", "SPEI-2", "SPEI-3"), each = 12),
  Fthm = c(
    100.5, 98.1, NA, 100.3, 101.5, 100.4, NA, 100.9, 101.0, 98.6, 101.5, NA,
    99.0, 99.2, NA, 95.7, 101.9, 99.5, NA, 101.2, 102.7, 94.3, 100.5, NA,
    94.9, 109.7, NA, 103.1, 98.6, 91.7, NA, 84.5, 122.5, 97.7, 96.5, NA
  ),
  Fthm_corr = c(
    100.1, 99.2, NA, 90.6, 14.2, 87.6, NA, 78.3, 72.3, 77.0, 103.9, NA,
    99.5, 99.8, NA, 90.1, 59.2, 68.7, NA, 82.9, 48.4, 35.1, 83.2, NA,
    96.0, 111.0, NA, 98.9, 60.8, 62.9, NA, 52.9, 62.9, 53.5, 40.8, NA
  )
)

# Convert Month to ordered factor
df_thm$Month <- factor(df_thm$Month, levels = month.abb[1:12])

# Set a minimum threshold difference for showing the ribbon (e.g., 2.0%)
ribbon_threshold <- 2.0

# Create a filtered ribbon dataset to omit small variations
df_thm_ribbon <- df_thm %>%
  mutate(
    abs_diff = abs(Fthm - Fthm_corr),
    ymin_val = if_else(!is.na(abs_diff) & abs_diff >= ribbon_threshold, pmin(Fthm, Fthm_corr), NA_real_),
    ymax_val = if_else(!is.na(abs_diff) & abs_diff >= ribbon_threshold, pmax(Fthm, Fthm_corr), NA_real_)
  )

# Pivot to long format for clean point mapping
df_thm_long <- df_thm %>%
  pivot_longer(cols = c(Fthm, Fthm_corr), names_to = "Type", values_to = "Value") %>%
  mutate(Type = factor(Type,
                       levels = c("Fthm", "Fthm_corr"),
                       labels = c("Original", "Bias-Corrected")))

# Identify excluded months (numeric positions)
excluded_months_num <- c(3, 7, 12)  # Mar, Jul, Dec
rect_data_thm <- data.frame(
  xmin = excluded_months_num - 0.5,
  xmax = excluded_months_num + 0.5
)

# ── 1b. Create the plot ───────────────────────────────────────────────────────
p_table3 <- ggplot() +
  # Background shading for excluded months
  geom_rect(
    data = rect_data_thm,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "#EBEBEB", alpha = 1.0, inherit.aes = FALSE
  ) +
  # 50% threshold benchmark line
  geom_hline(yintercept = 50, color = "grey80", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  # Shaded correction band
  geom_ribbon(
    data = df_thm_ribbon,
    aes(x = Month_num,
        ymin = ymin_val,
        ymax = ymax_val,
        fill = "Correction Magnitude (\u0394Fthm)"),
    alpha = 0.15, na.rm = TRUE
  ) +
  # Points layer
  geom_point(
    data = df_thm_long,
    aes(x = Month_num, y = Value, color = Type, shape = Type, size = Type),
    na.rm = TRUE
  ) +
  # Faceting
  facet_wrap(
    ~ Scale,
    ncol = 1,
    labeller = labeller(Scale = function(x) paste0("(", letters[1:3], ") ", x))
  ) +
  # Unified Styling Scales
  scale_color_manual(name = "Thermodynamic Fraction", values = c("Original" = "#4A5568", "Bias-Corrected" = "#E53E3E")) +
  scale_shape_manual(name = "Thermodynamic Fraction", values = c("Original" = 16, "Bias-Corrected" = 15)) +
  scale_size_manual(name = "Thermodynamic Fraction", values = c("Original" = 2.8, "Bias-Corrected" = 1.4)) +
  scale_fill_manual(name = NULL, values = c("Correction Magnitude (\u0394Fthm)" = "#E53E3E")) +
  # X and Y axis formatting
  scale_x_continuous(breaks = 1:12, labels = month.abb[1:12], expand = c(0.02, 0.02)) +
  scale_y_continuous(
    name = "Thermodynamic Fraction (%)",
    limits = c(0, 130),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0)
  ) +
  # Professional Journal Labels
  labs(
    x = "Calendar Month",
    title = "Thermodynamic Fraction Sensitivity to Reanalysis Bias",
    subtitle = "2022\u20132025 Nechako River Basin Drought | Active Drought Months",
    caption = "Grey shaded areas denote calendar months excluded from the attribution analysis due to a lack of active drought observations in 2022-2025."
  ) +
  # Polished Theme Settings
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#1A202C", margin = margin(b = 6)),
    plot.subtitle = element_text(size = 10, color = "#4A5568", margin = margin(b = 15)),
    plot.caption = element_text(color = "#718096", size = 8.5, hjust = 0, margin = margin(t = 12)),
    plot.margin = margin(t = 5.5, r = 5.5, b = 2, l = 5.5, unit = "pt"),
    axis.text.x = element_text(size = 11, color = "#2D3748", face = "plain"),
    axis.text.y = element_text(size = 11, color = "#2D3748", face = "plain"),
    strip.text = element_text(face = "bold", size = 11, color = "#1A202C", hjust = 0),
    strip.background = element_rect(fill = "#F7FAFC", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EDF2F7", linewidth = 0.5),
    panel.border = element_rect(color = "#CBD5E0", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.title = element_text(face = "bold", size = 12, color = "#2D3748"),
    legend.text = element_text(size = 11, color = "#2D3748")
  )

# ── 1c. Display plot ──────────────────────────────────────────────────────────
print(p_table3)
cat("  \u2713 Table 3 plot rendered.\n")


# ==============================================================================
# SECTION 2 — MANUSCRIPT SECTION 4.3 BLANKS (SW08 Event Catalog)
# ==============================================================================
# SOURCE: manuscript_statistics_extraction.R  (section43_manuscript_blanks.R)
# PURPOSE: Extract every statistic needed for manuscript Section 4.3 blanks.
#          Re-implements the SW08 event loop from w5_event_ranking.R.
# OUTPUT:  section43_blanks_summary.csv  +  Table1_top_events_SPEI3_x05.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(lubridate)
})

cat("\n============================================================\n")
cat(" SECTION 2: SECTION 4.3 MANUSCRIPT BLANKS\n")
cat("            SPEI-3 & SPI-3 SW08 EVENTS\n")

# ── User settings (shared by Sections 2 and 3) ───────────────────────────────
WD_PATH          <- "D:/Nechako_Drought/Nechako/"
SPEI_SEAS_DIR    <- file.path(WD_PATH, "spei_results_seasonal")
SPI_SEAS_DIR     <- file.path(WD_PATH, "spi_results_seasonal")
RANKING_DIR      <- file.path(WD_PATH, "temporal_drought", "event_ranking")

HISTORICAL_START <- 1950L
HISTORICAL_END   <- 2025L
FOCUS_START      <- 2022L
FOCUS_END        <- 2025L
ONSET_THR        <- -0.5   # x0 Yevjevich / Lloyd-Hughes threshold
MIN_DUR_SPEI3    <- 3L     # minimum event duration for scale-3 indices
n_years_record   <- HISTORICAL_END - HISTORICAL_START + 1L
MONTH_NAMES      <- c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec")

setwd(WD_PATH)

cat(sprintf(" Threshold x\u2080 = %.1f | Min duration = %d months\n",
            ONSET_THR, MIN_DUR_SPEI3))
cat("============================================================\n\n")

# ── Helper functions (used by both Sections 2 and 3) ─────────────────────────

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

# Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a) || length(a)==0 || is.na(a)) b else a

# ── Step 2.1: Build SW08 event catalogs — SPEI-3 and SPI-3 ──────────────────
cat("STEP 2.1: Building SW08 event catalogs\n")
cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")

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
  ev <- ev[order(-ev$severity), ]
  ev$rank_by_severity  <- seq_len(nrow(ev))
  ev$return_period_yrs <- (n_years_record + 1L) / ev$rank_by_severity
  ev$percentile_rank   <- round(100 * (1 - ev$rank_by_severity / nrow(ev)), 1)
  cat(sprintf("  [%s-%d] %d qualifying events detected\n",
              toupper(index_type), scale, nrow(ev)))
  out_f <- file.path(RANKING_DIR,
                     sprintf("%s-%02d_SW08_x05_event_catalog.csv",
                             toupper(index_type), scale))
  dir.create(RANKING_DIR, showWarnings=FALSE, recursive=TRUE)
  write.csv(ev, out_f, row.names=FALSE)
  ev
}

spei3_events <- build_catalog(SPEI_SEAS_DIR, "spei", 3L)
spi3_events  <- build_catalog(SPI_SEAS_DIR,  "spi",  3L)

# ── Step 2.2: Extract all Section 4.3 blank values — SPEI-3 ─────────────────
cat("\n\nSTEP 2.2: Extracting SPEI-3 statistics\n")
cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")

ev       <- spei3_events
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
cat(sprintf("Duration range: %d \u2013 %d months\n",     dur_min, dur_max))
cat(sprintf("Severity range: %.3f \u2013 %.3f index-units\u00b7months\n", sev_min, sev_max))
cat(sprintf("D3-6  (short-term):  %d events\n", n_d36))
cat(sprintf("D7-12 (medium-term): %d events\n", n_d712))
cat(sprintf("D12+  (long-term):   %d events\n", n_d12p))

# 2022–2025 event
cat("\n\u2500\u2500 2022\u20132025 event \u2500\u2500\n")
recent <- ev_valid[ev_valid$is_recent, ]

if (nrow(recent) == 0) {
  recent_all <- ev[ev$is_recent, ]
  if (nrow(recent_all) == 0)
    stop("No event overlapping 2022-2025 found \u2014 check data extends through 2025.")
  recent <- recent_all[which.max(recent_all$severity), , drop=FALSE]
  cat("  \u26a0 2022\u20132025 event classified as D1-2 (sub-threshold) \u2014 using largest overlap.\n")
}

rec_evt       <- recent[which.max(recent$severity), , drop=FALSE]
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
cat(sprintf("  Composite severity S: %.2f index-units\u00b7months\n", sev_2022))
cat(sprintf("  Severity rank: #%d of %d (%.1fth percentile)\n",
            rank_sev_2022, n_total, pct_2022))

ev_valid2 <- ev_valid[order(-ev_valid$mean_int), ]
ev_valid2$rank_by_intensity <- seq_len(nrow(ev_valid2))
rank_int_2022 <- ev_valid2$rank_by_intensity[
  match(paste(rec_evt$start_date, rec_evt$end_date),
        paste(ev_valid2$start_date, ev_valid2$end_date))]
cat(sprintf("  Intensity rank: #%d\n", rank_int_2022))

# Previous highest-ranked event
cat("\n\u2500\u2500 Previous highest-ranked event \u2500\u2500\n")
historical  <- ev_valid[!ev_valid$is_recent, ]
historical  <- historical[order(-historical$severity), ]
prev_top    <- historical[1, , drop=FALSE]
prev_name   <- format(as.Date(prev_top$start_date), "%Y")
prev_period <- sprintf("%s to %s",
                       format(as.Date(prev_top$start_date), "%B %Y"),
                       format(as.Date(prev_top$end_date),   "%B %Y"))
prev_sev    <- round(prev_top$severity, 2)
sev_ratio   <- round(sev_2022 / prev_sev, 2)

cat(sprintf("  Event period: %s\n", prev_period))
cat(sprintf("  Year label:   %s drought\n", prev_name))
cat(sprintf("  Severity:     %.2f index-units\u00b7months\n", prev_sev))
cat(sprintf("  Ratio 2022-2025 / prev: %.2f\u00d7\n", sev_ratio))

# Year-level composite severity scores
cat("\n\u2500\u2500 Year-level composite severity scores (2022\u20132025) \u2500\u2500\n")
ts_spei3 <- get_ts(SPEI_SEAS_DIR, "spei", 3L)
ts_spei3 <- ts_spei3[!is.na(ts_spei3$value), ]
ts_spei3$year  <- as.integer(format(ts_spei3$date, "%Y"))
ts_spei3$month <- as.integer(format(ts_spei3$date, "%m"))

annual_stats <- ts_spei3 %>%
  dplyr::filter(year >= FOCUS_START, year <= FOCUS_END) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    n_drought_months  = sum(value < ONSET_THR, na.rm=TRUE),
    mean_severity_idx = round(
      mean(ONSET_THR - value[value < ONSET_THR], na.rm=TRUE), 3),
    composite_score   = n_drought_months *
      mean(ONSET_THR - value[value < ONSET_THR], na.rm=TRUE),
    .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(composite_score))

cat("\nFocus period year-level stats (SPEI-3, drought months below x\u2080=-0.5):\n")
print(as.data.frame(annual_stats))

all_annual <- ts_spei3 %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    n_drought_months = sum(value < ONSET_THR, na.rm=TRUE),
    mean_sev         = mean(ONSET_THR - value[value < ONSET_THR], na.rm=TRUE),
    composite_score  = n_drought_months * mean(ONSET_THR - value[value < ONSET_THR], na.rm=TRUE),
    .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(composite_score)) %>%
  dplyr::mutate(annual_rank = dplyr::row_number())

cat("\nTop-10 calendar years by annual composite score (all years 1950\u20132025):\n")
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

# ── Step 2.3: SPI-3 comparison ───────────────────────────────────────────────
cat("\n\nSTEP 2.3: SPI-3 comparison for 2022\u20132025 event\n")
cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")

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
  spi_rec_evt   <- spi_recent[which.max(spi_recent$severity), , drop=FALSE]
  spi_sev_2022  <- round(spi_rec_evt$severity, 2)
  spi_rank_2022 <- spi_rec_evt$rank_by_severity
  spi_pct_2022  <- spi_rec_evt$percentile_rank
  spi_n_total   <- nrow(spi_ev_valid)
  
  cat(sprintf("SPI-3 2022\u20132025 event:\n"))
  cat(sprintf("  Composite severity S: %.2f index-units\u00b7months\n", spi_sev_2022))
  cat(sprintf("  Severity rank: #%d of %d (%.1fth percentile)\n",
              spi_rank_2022, spi_n_total, spi_pct_2022))
  cat(sprintf("  SPEI-3 rank (#%d) vs SPI-3 rank (#%d): ",
              rank_sev_2022, spi_rank_2022))
  if (rank_sev_2022 < spi_rank_2022)
    cat("SPEI ranks MORE severe \u2192 thermodynamic amplification confirmed.\n")
  else
    cat("SPI ranks equally or more severe.\n")
} else {
  cat("  \u26a0 No SPI-3 event detected overlapping 2022\u20132025.\n")
  spi_sev_2022  <- NA; spi_rank_2022 <- NA
}

# ── Step 2.4: Generate Table 1 ───────────────────────────────────────────────
cat("\n\nSTEP 2.4: Generating Table 1 (top-ranked SPEI-3 events)\n")
cat("\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n")

table1 <- ev_valid %>%
  dplyr::filter(rank_by_severity <= 15 | is_recent) %>%
  dplyr::arrange(rank_by_severity) %>%
  dplyr::mutate(
    Event_period     = sprintf("%s \u2013 %s",
                               format(as.Date(start_date), "%b %Y"),
                               format(as.Date(end_date),   "%b %Y")),
    Rank             = rank_by_severity,
    Duration_months  = duration_months,
    Duration_class   = duration_class,
    Mean_intensity_I = round(mean_int, 3),
    Severity_S       = round(severity, 2),
    Peak_SPEI        = NA_real_,
    `Is_2022-2025`   = is_recent
  ) %>%
  dplyr::select(Rank, Event_period, Duration_months, Duration_class,
                Mean_intensity_I, Severity_S, `Is_2022-2025`)

cat("\nTable 1 \u2014 Top-15 SPEI-3 drought events (1950\u20132025):\n")
print(as.data.frame(table1), row.names=FALSE)

t1_path <- file.path(RANKING_DIR, "Table1_top_events_SPEI3_x05.csv")
write.csv(table1, t1_path, row.names=FALSE)
cat(sprintf("\n\u2713 Table 1 saved: %s\n", t1_path))

# ── Step 2.5: Assemble complete blanks summary ───────────────────────────────
cat("\n\nSTEP 2.5: Complete blank-value summary for Section 4.3\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

blanks_s43 <- data.frame(
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
    "P2: prev top event \u2014 year label",
    "P2: prev top event \u2014 period",
    "P2: prev top event \u2014 severity S",
    "P2: 2022-2025 / prev ratio",
    "P2: 2023 \u2014 n drought months",
    "P2: 2023 \u2014 mean severity index",
    "P2: 2023 \u2014 annual rank",
    "P2: 2024 \u2014 n drought months",
    "P2: 2024 \u2014 mean severity index",
    "P2: 2024 \u2014 annual rank",
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
    "All D3+ events; I\u00d7D", "All D3+ events; I\u00d7D",
    "Sheffield & Wood (2008) class", "Sheffield & Wood (2008) class",
    "Sheffield & Wood (2008) class",
    "Contiguous months below x\u2080=-0.5",
    "First month below threshold", "Last month below threshold",
    "mean(x\u2080 - SPEI_t) over event months",
    "I \u00d7 D", "Rank 1 = most severe",
    "Rank 1 = most intense",
    "Year of previous top historical event",
    "Full date range of previous top event",
    "I \u00d7 D of previous top event",
    "sev_2022 / prev_sev",
    "Months with SPEI-3 < -0.5 in 2023",
    "mean(x\u2080 - SPEI_t) for drought months in 2023",
    "1 = highest annual composite score",
    "Months with SPEI-3 < -0.5 in 2024",
    "mean(x\u2080 - SPEI_t) for drought months in 2024",
    "1 = highest annual composite score",
    "SPI-3 I \u00d7 D for 2022-2025 event",
    "SPI-3 rank; higher # \u2192 less severe than SPEI"
  ),
  stringsAsFactors = FALSE
)

print(blanks_s43, right=FALSE)

out_path_s43 <- file.path(WD_PATH, "section43_blanks_summary.csv")
write.csv(blanks_s43, out_path_s43, row.names=FALSE)
cat(sprintf("\n\u2713 All blank values saved to: %s\n", out_path_s43))
cat(sprintf("\u2713 Table 1 saved to:          %s\n", t1_path))
cat("============================================================\n")
# ── Step 2.6: Extract Table 1 from bias_manuscript_tables.xlsx ───────────────
cat("\n\nSTEP 2.6: Extracting Table 1 from bias_manuscript_tables.xlsx\n")
cat("──────────────────────────────────────────────────────────────────────\n")

# Load openxlsx for reading formatted Excel files
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
suppressPackageStartupMessages(library(openxlsx))

bias_tables_file <- file.path(WD_PATH, "decomp_results", "bias_manuscript_tables.xlsx")

if (file.exists(bias_tables_file)) {
  sheets <- getSheetNames(bias_tables_file)
  cat("  Available sheets:", paste(sheets, collapse = ", "), "\n")
  
  # Identify the sheet containing Table 1 data (adjust regex if sheet name differs)
  target_sheet <- grep("Table.?1|Summary.?Statistics|SPEI.?Mean", sheets, ignore.case = TRUE, value = TRUE)
  
  if (length(target_sheet) > 0) {
    target_sheet <- target_sheet[1]
    cat(sprintf("  Extracting data from sheet: '%s'\n", target_sheet))
    
    table1_data <- read.xlsx(bias_tables_file, sheet = target_sheet)
    
    # Save as a clean CSV for manuscript insertion
    t1_out_path <- file.path(WD_PATH, "Table1_Summary_Statistics.csv")
    write.csv(table1_data, t1_out_path, row.names = FALSE)
    cat(sprintf("  ✓ Table 1 data successfully saved to: %s\n", t1_out_path))
    
    cat("\n  Preview of Table 1 data:\n")
    print(head(table1_data))
  } else {
    cat("  ⚠ Could not automatically identify the 'Table 1' sheet.\n")
    cat("  Please manually set `target_sheet <- 'Exact_Sheet_Name'` if needed.\n")
  }
} else {
  cat(sprintf("  ⚠ File not found: %s\n", bias_tables_file))
  cat("  → Please run 10b_PET_bias_nonstationarity.R first to generate this file.\n")
}

# ==============================================================================
# SECTION 3 — MS1 BLANK-FILLING EXTRACTOR (ALL SECTIONS)
# ==============================================================================
# SOURCE: w_extract_ms1_blanks.R
# PURPOSE: Extract all numerical values needed to fill blanks in
#          Nechako_MS1_REVISED.docx. Depends on outputs from w5, w11, w9, w7.
# OUTPUT:  ms1_blanks/MS1_BLANKS_MASTER.csv and per-section CSVs
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(zoo)
})

OUT_DIR <- file.path(WD_PATH, "ms1_blanks")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("  SECTION 3: MS1 BLANK-FILLING EXTRACTOR\n")
cat("============================================================\n\n")

all_blanks <- list()  # accumulates rows for the master table

add_blank <- function(section, location_text, variable, value, source_file) {
  all_blanks[[length(all_blanks) + 1L]] <<- data.frame(
    section       = section,
    location      = location_text,
    variable      = variable,
    value         = as.character(round(as.numeric(value), 4)),
    source        = source_file,
    stringsAsFactors = FALSE
  )
}

# ── Block 3.1: Values already confirmed from Table 3 / event catalog ─────────
cat("\u2500\u2500 Block 3.1: Confirmed values from Table 3 and event catalog \u2500\u2500\n")

confirmed <- data.frame(
  section   = c("ABS","ABS","ABS","ABS",
                "S431","S431","S431","S431",
                "S51","S51","S51",
                "S52","S52",
                "C1","C1","C1","C1","C1","C1",
                "C2","C2"),
  location  = c(
    "Abstract: T_K for leading SPEI-3 episode",
    "Abstract: 95% CI lower bound",
    "Abstract: 95% CI upper bound",
    "Abstract: thermodynamic fraction 'nearly 100%'",
    "Sec 4.3: SPEI-3 rank-1 S value",
    "Sec 4.3: SPEI-3 rank-2 S value",
    "Sec 4.3: ratio rank-1 / prev best (5.16)",
    "Sec 4.3: SPI-3 event S value (21-month)",
    "Sec 5.1: composite S at SPEI-3 rank-1",
    "Sec 5.1: factor exceeding prev best",
    "Sec 5.1: T_K years",
    "Sec 5.2: thermodynamic % (SPEI-3)",
    "Sec 5.2: dynamic residual %",
    "Conclusion 1: S at SPEI-3",
    "Conclusion 1: factor exceeding 2nd-ranked",
    "Conclusion 1: 2023 mean severity index",
    "Conclusion 1: 2024 mean severity index",
    "Conclusion 1: T_K years",
    "Conclusion 1: 95% CI lower\u2013upper",
    "Conclusion 2: thermodynamic %",
    "Conclusion 2: dynamic residual %"
  ),
  variable  = c("T_K_years","CI_lower","CI_upper","thm_pct_approx",
                "SPEI3_rank1_S","SPEI3_rank2_S","ratio_vs_prev_best","SPI3_S_21mo",
                "S51_S_SPEI3","S51_factor","S51_TK",
                "S52_thm_pct","S52_dyn_pct",
                "C1_S","C1_factor","C1_2023_sev","C1_2024_sev",
                "C1_TK","C1_CI",
                "C2_thm_pct","C2_dyn_pct"),
  value     = c("921","109","55386","~100",
                "10.73","10.31","2.08","48.41",
                "10.73","2.08","921",
                "~99","~1",
                "10.73","2.08","1.329","1.056",
                "921","109\u201355,386",
                "~99","~1"),
  source    = "Table 3 + ranked_event_catalog_SW08.csv (confirmed)",
  stringsAsFactors = FALSE
)

write.csv(confirmed, file.path(OUT_DIR, "ms1_blanks_confirmed.csv"), row.names = FALSE)
cat(sprintf("  \u2713 Confirmed values: %d blanks \u2192 ms1_blanks_confirmed.csv\n", nrow(confirmed)))

# ── Block 3.2: Decomposition values from w11 output ──────────────────────────
cat("\n\u2500\u2500 Block 3.2: Decomposition values from w11 \u2500\u2500\n")

decomp_file <- file.path(WD_PATH, "decomp_results", "decomp_2022_2025_summary.csv")

if (file.exists(decomp_file)) {
  decomp <- read.csv(decomp_file)
  cat("  decomp_2022_2025_summary.csv loaded\n")
  print(decomp)
  
  for (sc in c(1, 2, 3)) {
    row <- decomp[decomp$scale == sc, ]
    if (nrow(row) == 0) next
    add_blank("S44", sprintf("Sec 4.4: SPEI_PM mean (scale %d)", sc),
              sprintf("SPEI_PM_mean_s%d", sc), row$SPEI_PM_mean, decomp_file)
    add_blank("S44", sprintf("Sec 4.4: SPEI_Thw mean (scale %d)", sc),
              sprintf("SPEI_Thw_mean_s%d", sc), row$SPEI_Thw_mean, decomp_file)
    add_blank("S44", sprintf("Sec 4.4: SPEI_Dynamic mean (scale %d)", sc),
              sprintf("SPEI_Dyn_mean_s%d", sc), row$SPEI_Dynamic_mean, decomp_file)
    add_blank("S44", sprintf("Sec 4.4: thm_pct %% (scale %d)", sc),
              sprintf("thm_pct_s%d", sc), row$thm_pct, decomp_file)
    add_blank("S44", sprintf("Sec 4.4: dyn_pct %% (scale %d)", sc),
              sprintf("dyn_pct_s%d", sc), row$dyn_pct, decomp_file)
  }
  
  cat("\n  KEY VALUES FOR MANUSCRIPT:\n")
  s3 <- decomp[decomp$scale == 3, ]
  cat(sprintf("  SPEI-3: PM=%.3f  Thw=%.3f  Dyn=%.3f  thm_pct=%.1f%%  dyn_pct=%.1f%%\n",
              s3$SPEI_PM_mean, s3$SPEI_Thw_mean, s3$SPEI_Dynamic_mean,
              s3$thm_pct, s3$dyn_pct))
  s1 <- decomp[decomp$scale == 1, ]
  s2 <- decomp[decomp$scale == 2, ]
  cat(sprintf("  SPEI-1: thm_pct=%.1f%%   SPEI-2: thm_pct=%.1f%%\n",
              s1$thm_pct, s2$thm_pct))
  
  write.csv(do.call(rbind, all_blanks[sapply(all_blanks, function(x) x$section == "S44")]),
            file.path(OUT_DIR, "ms1_blanks_decomp.csv"), row.names = FALSE)
  cat("  \u2713 Decomposition blanks \u2192 ms1_blanks_decomp.csv\n")
  
} else {
  cat(sprintf("  \u26a0 Not found: %s\n  \u2192 Run w11_dynamic_thermodynamic_decomp.R first\n", decomp_file))
}

# ── Block 3.3: Offset months (SPEI_\u0394 > 0 during 2022–2025) ─────────────────────
cat("\n\u2500\u2500 Block 3.3: Offset months (SPEI_\u0394 > 0 during 2022\u20132025) \u2500\u2500\n")

decomp_ts_file <- file.path(WD_PATH, "decomp_results", "decomp_timeseries.csv")

if (file.exists(decomp_ts_file)) {
  decomp_ts       <- read.csv(decomp_ts_file)
  decomp_ts$date  <- as.Date(decomp_ts$date)
  decomp_ts$year  <- as.integer(format(decomp_ts$date, "%Y"))
  decomp_ts$month <- as.integer(format(decomp_ts$date, "%m"))
  
  offset <- decomp_ts %>%
    filter(year >= 2022, year <= 2025) %>%
    group_by(scale) %>%
    summarise(
      n_total      = n(),
      n_positive   = sum(SPEI_Dynamic > 0, na.rm = TRUE),
      pct_positive = 100 * mean(SPEI_Dynamic > 0, na.rm = TRUE),
      mean_dynamic = mean(SPEI_Dynamic, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n  Months with SPEI_\u0394 > 0 (2022\u20132025 window):\n")
  print(offset)
  write.csv(offset, file.path(OUT_DIR, "ms1_blanks_offset_months.csv"), row.names = FALSE)
  
  for (sc in c(1, 2, 3)) {
    row <- offset[offset$scale == sc, ]
    if (nrow(row) == 0) next
    add_blank("S44", sprintf("Sec 4.4: months SPEI_\u0394>0 (scale %d, n out of 48)", sc),
              sprintf("n_offset_s%d", sc), row$n_positive, decomp_ts_file)
    add_blank("S44", sprintf("Sec 4.4: mean SPEI_\u0394 over event window (scale %d)", sc),
              sprintf("mean_dyn_s%d", sc), row$mean_dynamic, decomp_ts_file)
  }
  
  s3o <- offset[offset$scale == 3, ]
  cat(sprintf("\n  SPEI-3: %d of 48 months have SPEI_\u0394 > 0 (%.0f%%)\n",
              s3o$n_positive, s3o$pct_positive))
  cat(sprintf("  SPEI-3: mean SPEI_\u0394 over full window = %.4f\n", s3o$mean_dynamic))
  cat("  \u2713 Offset months \u2192 ms1_blanks_offset_months.csv\n")
  
} else {
  cat(sprintf("  \u26a0 %s not found\n", decomp_ts_file))
  cat("  \u2192 Add write.csv(decomp_all, file.path(out_dir,\"decomp_timeseries.csv\"), row.names=FALSE) to w11 STEP 4.\n")
}

# ── Block 3.4: Trend values from w11 Step 5 ──────────────────────────────────
cat("\n\u2500\u2500 Block 3.4: JJA thermodynamic fraction trend (Section 4.6) \u2500\u2500\n")

trend_file <- file.path(WD_PATH, "decomp_results", "ms_blanks_trend.csv")

if (file.exists(trend_file)) {
  trend <- read.csv(trend_file)
  cat("  ms_blanks_trend.csv loaded\n")
  print(trend)
  
  for (i in seq_len(nrow(trend))) {
    r <- trend[i, ]
    add_blank("S46", paste("Sec 4.6:", r$metric), r$metric, r$value, trend_file)
  }
  
  slope_dec  <- trend$value[trend$metric == "JJA_trend_slope_per_decade"]
  pval       <- trend$value[trend$metric == "JJA_trend_pval"]
  tstat      <- trend$value[trend$metric == "JJA_trend_tstat"]
  mean_early <- trend$value[trend$metric == "mean_f_thm_1950_1987_pct"]
  mean_late  <- trend$value[trend$metric == "mean_f_thm_1988_2025_pct"]
  mw_p       <- trend$value[trend$metric == "mann_whitney_p"]
  
  sig_word <- if (!is.na(pval) && as.numeric(pval) < 0.05) "significant" else "non-significant"
  cat("\n  KEY VALUES FOR SECTION 4.6:\n")
  cat(sprintf("  Trend: statistically %s positive trend of %.4f %% per decade\n",
              sig_word, as.numeric(slope_dec) * 100))
  cat(sprintf("  OLS slope = %.6f yr\u207b\u00b9; t = %.3f; p = %.4f\n",
              as.numeric(slope_dec)/10, as.numeric(tstat), as.numeric(pval)))
  cat(sprintf("  1950\u20131987 mean JJA f_thm = %.1f%%\n", as.numeric(mean_early)))
  cat(sprintf("  1988\u20132025 mean JJA f_thm = %.1f%%\n", as.numeric(mean_late)))
  cat(sprintf("  Change = %.1f percentage points; Mann\u2013Whitney p = %.4f\n",
              as.numeric(mean_late) - as.numeric(mean_early), as.numeric(mw_p)))
  
  write.csv(do.call(rbind, Filter(function(x) x$section=="S46", all_blanks)),
            file.path(OUT_DIR, "ms1_blanks_trend.csv"), row.names = FALSE)
  cat("  \u2713 Trend blanks \u2192 ms1_blanks_trend.csv\n")
  
} else {
  cat(sprintf("  \u26a0 %s not found\n  \u2192 Add ms_blanks_trend output block to END OF STEP 5 in w11.\n", trend_file))
}

# ── Block 3.5: Seasonal decomposition fractions ───────────────────────────────
cat("\n\u2500\u2500 Block 3.5: Seasonal thermodynamic fractions (Section 4.4) \u2500\u2500\n")

seasonal_file <- file.path(WD_PATH, "decomp_results", "decomp_seasonal_fractions.csv")

if (file.exists(seasonal_file)) {
  seas <- read.csv(seasonal_file)
  cat("  decomp_seasonal_fractions.csv loaded\n")
  print(seas)
  
  for (i in seq_len(nrow(seas))) {
    add_blank("S44", paste("Sec 4.4: seasonal f_thm", seas$season[i]),
              paste0("thm_frac_", seas$season[i]), seas$thm_frac_pct[i], seasonal_file)
  }
  
  jja_pct <- seas$thm_frac_pct[seas$season == "JJA"]
  mam_pct <- seas$thm_frac_pct[seas$season == "MAM"]
  cat(sprintf("\n  JJA thermodynamic fraction: %.1f%%\n", jja_pct))
  cat(sprintf("  MAM (spring) thermodynamic fraction: %.1f%%\n", mam_pct))
  cat("  \u2713 Seasonal fractions \u2192 ms1_blanks_seasonal.csv\n")
  write.csv(seas, file.path(OUT_DIR, "ms1_blanks_seasonal.csv"), row.names = FALSE)
  
} else {
  cat(sprintf("  \u26a0 %s not found\n  \u2192 Add decomp_seasonal_fractions output block to STEP 4 in w11.\n", seasonal_file))
}

# ── Block 3.6: Annual severity scores from event catalog ─────────────────────
cat("\n\u2500\u2500 Block 3.6: Annual severity scores from event catalog \u2500\u2500\n")

catalog_file <- file.path(WD_PATH, "temporal_drought",
                          "drought_event_ranking", "ranked_event_catalog_SW08.csv")
if (!file.exists(catalog_file)) {
  candidates <- list.files(WD_PATH, "ranked_event_catalog", recursive = TRUE, full.names = TRUE)
  if (length(candidates)) catalog_file <- candidates[1]
}

if (file.exists(catalog_file)) {
  cat_df  <- read.csv(catalog_file)
  spei3_c <- cat_df[cat_df$index %in% c("spei_03", "spei3", "SPEI3", "spei-3"), ]
  
  if (nrow(spei3_c) > 0) {
    spei3_c$year <- as.integer(format(as.Date(spei3_c$start_date), "%Y"))
    annual_scores_cat <- spei3_c %>%
      group_by(year) %>%
      summarise(
        n_drought_months = sum(duration_months, na.rm = TRUE),
        mean_severity    = mean(mean_int, na.rm = TRUE),
        annual_score     = sum(severity, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(annual_score))
    
    cat("\n  Top annual SPEI-3 severity scores:\n")
    print(head(annual_scores_cat, 5))
    
    yr2023_cat <- annual_scores_cat[annual_scores_cat$year == 2023, ]
    yr2024_cat <- annual_scores_cat[annual_scores_cat$year == 2024, ]
    
    if (nrow(yr2023_cat) > 0) {
      add_blank("S431","Sec 4.3: 2023 drought months (SPEI-3)",
                "n_drought_months_2023", yr2023_cat$n_drought_months, catalog_file)
      add_blank("S431","Sec 4.3: 2023 mean severity index (SPEI-3)",
                "mean_sev_2023", yr2023_cat$mean_severity, catalog_file)
    }
    if (nrow(yr2024_cat) > 0) {
      add_blank("S431","Sec 4.3: 2024 drought months (SPEI-3)",
                "n_drought_months_2024", yr2024_cat$n_drought_months, catalog_file)
      add_blank("S431","Sec 4.3: 2024 mean severity index (SPEI-3)",
                "mean_sev_2024", yr2024_cat$mean_severity, catalog_file)
    }
    write.csv(annual_scores_cat, file.path(OUT_DIR, "ms1_blanks_annual_scores.csv"), row.names = FALSE)
    cat("  \u2713 Annual scores \u2192 ms1_blanks_annual_scores.csv\n")
  } else {
    cat("  \u26a0 No SPEI-3 rows found in catalog. Check index label column.\n")
    cat("  Available index values:", unique(cat_df$index), "\n")
  }
} else {
  cat("  \u26a0 Event catalog CSV not found.\n  \u2192 Run w5_event_ranking.R first.\n")
}

# ── Block 3.7: PDO and ENSO indices ──────────────────────────────────────────
cat("\n\u2500\u2500 Block 3.7: PDO / Ni\u00f1o-3.4 indices \u2500\u2500\n")

pdo_candidates  <- list.files(WD_PATH, pattern = "pdo.*\\.csv",
                              recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
nino_candidates <- list.files(WD_PATH, pattern = "nino.*\\.csv|enso.*\\.csv",
                              recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

if (length(pdo_candidates) > 0) {
  pdo_file <- pdo_candidates[1]
  cat(sprintf("  Loading PDO: %s\n", basename(pdo_file)))
  pdo_df <- tryCatch(read.csv(pdo_file), error = function(e) NULL)
  
  if (!is.null(pdo_df)) {
    names(pdo_df) <- tolower(names(pdo_df))
    date_col  <- grep("date|year|time", names(pdo_df), value = TRUE)[1]
    index_col <- grep("pdo|index|value", names(pdo_df), value = TRUE)[1]
    
    if (!is.na(date_col) && !is.na(index_col)) {
      pdo_df$date_parsed <- suppressWarnings(as.Date(pdo_df[[date_col]]))
      pdo_df$yr <- as.integer(format(pdo_df$date_parsed, "%Y"))
      pdo_2022_2025 <- pdo_df[pdo_df$yr >= 2022 & pdo_df$yr <= 2025, ]
      mean_pdo  <- mean(pdo_2022_2025[[index_col]], na.rm = TRUE)
      pdo_phase <- if (mean_pdo > 0) "positive" else "negative"
      cat(sprintf("  PDO mean 2022\u20132025: %.3f (%s phase)\n", mean_pdo, pdo_phase))
      
      pdo_annual <- pdo_df %>%
        group_by(yr) %>%
        summarise(mean_pdo = mean(.data[[index_col]], na.rm=TRUE), .groups="drop") %>%
        arrange(yr)
      sign_changes <- which(diff(sign(pdo_annual$mean_pdo)) > 0)
      if (length(sign_changes) > 0) {
        last_pos_transition <- pdo_annual$yr[tail(sign_changes[pdo_annual$yr[sign_changes] < 2022], 1) + 1]
        cat(sprintf("  Last PDO positive transition: ~%s\n", last_pos_transition))
        add_blank("S53","Sec 5.3: PDO transition to positive phase year",
                  "PDO_transition_year", last_pos_transition, pdo_file)
      }
      add_blank("S47","Sec 4.7: PDO mean 2022\u20132025",  "PDO_mean_2022_2025", mean_pdo, pdo_file)
      add_blank("S47","Sec 4.7: PDO phase 2022\u20132025", "PDO_phase", pdo_phase, pdo_file)
      add_blank("S53","Sec 5.3: PDO phase description",  "PDO_phase_S53", pdo_phase, pdo_file)
    }
  }
} else {
  cat("  \u26a0 PDO CSV not found. Using literature defaults.\n")
  add_blank("S47","Sec 4.7: PDO phase (from literature)","PDO_phase","positive","Newman et al. 2016")
  add_blank("S53","Sec 5.3: PDO transition year (from literature)","PDO_transition_year","2014","Newman et al. 2016")
}

if (length(nino_candidates) > 0) {
  nino_file <- nino_candidates[1]
  cat(sprintf("  Loading Ni\u00f1o-3.4: %s\n", basename(nino_file)))
  nino_df <- tryCatch(read.csv(nino_file), error = function(e) NULL)
  
  if (!is.null(nino_df)) {
    names(nino_df) <- tolower(names(nino_df))
    date_col  <- grep("date|year|time", names(nino_df), value = TRUE)[1]
    index_col <- grep("nino|sst|index|value|34", names(nino_df), value = TRUE)[1]
    
    if (!is.na(date_col) && !is.na(index_col)) {
      nino_df$date_parsed <- suppressWarnings(as.Date(nino_df[[date_col]]))
      nino_df$yr <- as.integer(format(nino_df$date_parsed, "%Y"))
      nino_event <- nino_df[nino_df$yr >= 2022 & nino_df$yr <= 2025, ]
      mean_nino  <- mean(nino_event[[index_col]], na.rm = TRUE)
      peak_nino  <- max(nino_event[[index_col]], na.rm = TRUE)
      peak_month <- nino_event$date_parsed[which.max(nino_event[[index_col]])]
      
      nino_annual <- nino_event %>%
        group_by(yr) %>%
        summarise(mean_nino34 = mean(.data[[index_col]], na.rm=TRUE), .groups="drop")
      nino_annual$phase <- ifelse(nino_annual$mean_nino34 > 0.5, "El Ni\u00f1o",
                                  ifelse(nino_annual$mean_nino34 < -0.5, "La Ni\u00f1a", "Neutral"))
      
      cat(sprintf("  Ni\u00f1o-3.4 mean 2022\u20132025: %.3f; peak: %.3f (%s)\n",
                  mean_nino, peak_nino, format(peak_month, "%b %Y")))
      cat("  Annual ENSO phases:\n")
      print(nino_annual)
      
      enso_2022_phase <- nino_annual$phase[nino_annual$yr == 2022]
      add_blank("S47","Sec 4.7: Ni\u00f1o-3.4 mean 2022\u20132025",  "Nino34_mean", mean_nino, nino_file)
      add_blank("S47","Sec 4.7: Ni\u00f1o-3.4 peak value",       "Nino34_peak", peak_nino, nino_file)
      add_blank("S53","Sec 5.3: ENSO state 2022",           "ENSO_2022", enso_2022_phase, nino_file)
      add_blank("S53","Sec 5.3: El Ni\u00f1o peak Ni\u00f1o-3.4",    "Nino34_peak_S53", peak_nino, nino_file)
    }
  }
} else {
  cat("  \u26a0 Ni\u00f1o-3.4 CSV not found. Using literature defaults.\n")
  add_blank("S47","Sec 4.7: Ni\u00f1o-3.4 mean (from literature)","Nino34_mean","~0.5","L'Heureux et al. 2024")
  add_blank("S47","Sec 4.7: Ni\u00f1o-3.4 peak (from literature)","Nino34_peak","~2.0","L'Heureux et al. 2024")
  add_blank("S53","Sec 5.3: ENSO 2022 phase","ENSO_2022","La Ni\u00f1a","L'Heureux et al. 2024")
}

# ── Block 3.8: Atmospheric composite values (Section 4.7) ────────────────────
cat("\n\u2500\u2500 Block 3.8: Atmospheric composite values (Section 4.7) \u2500\u2500\n")

atm_file <- file.path(WD_PATH, "composite_ms_blanks.csv")

if (file.exists(atm_file)) {
  atm <- read.csv(atm_file)
  cat("  composite_ms_blanks.csv loaded\n")
  print(atm)
  for (i in seq_len(nrow(atm))) {
    add_blank("S47", paste("Sec 4.7:", atm$metric[i]), atm$metric[i], atm$value[i], atm_file)
  }
  write.csv(atm, file.path(OUT_DIR, "ms1_blanks_atm.csv"), row.names = FALSE)
  cat("  \u2713 Atmospheric blanks \u2192 ms1_blanks_atm.csv\n")
} else {
  cat(sprintf("  \u26a0 %s not found.\n  \u2192 Run w9_atmospheric_diagnostics.R or w10a first.\n", atm_file))
  for (nm in c("Z500_ridge_max_m","Z500_ridge_location","Z500_trough_min_m",
               "SLP_max_hPa","SLP_location","SST_NE_Pac_warm_anom_C","IVT_reduction_pct")) {
    add_blank("S47", paste("Sec 4.7:", nm), nm, "[run w9/w10a]", "composite_ms_blanks.csv")
  }
}

# ── Block 3.9: Section 5.5 water management blanks ───────────────────────────
cat("\n\u2500\u2500 Block 3.9: Section 5.5 water management blanks \u2500\u2500\n")
cat("  BLANK A: EFN consecutive months \u2192 see Bradford et al. (2011)\n")
cat("  BLANK B: 21-month unbroken SPI-3 deficit\n")
cat("  BLANK C: mid-century frequency increase: 25\u201340% (Tam et al. 2023)\n")
cat("  BLANK D: end-century frequency increase: 60\u201380% (Tam et al. 2023)\n")

add_blank("S55","Sec 5.5: EFN consecutive months","EFN_consecutive_months","[see Bradford 2011]","operational records")
add_blank("S55","Sec 5.5: unbroken SPEI deficit months","unbroken_deficit_months","21 (SPI-3)","ranked_event_catalog_SW08.csv")
add_blank("S55","Sec 5.5: drought frequency increase mid-century %","freq_increase_midcentury_pct","25\u201340","Tam et al. 2023")
add_blank("S55","Sec 5.5: drought frequency increase end-century %","freq_increase_endcentury_pct","60\u201380","Tam et al. 2023")

# ── Master table ──────────────────────────────────────────────────────────────
cat("\n\u2500\u2500 Writing master blank table \u2500\u2500\n")

if (length(all_blanks) > 0) {
  master <- do.call(rbind, all_blanks)
  master <- master[order(master$section), ]
  write.csv(master, file.path(OUT_DIR, "MS1_BLANKS_MASTER.csv"), row.names = FALSE)
  cat(sprintf("  \u2713 Master table: %d blanks \u2192 MS1_BLANKS_MASTER.csv\n", nrow(master)))
  cat("\n  Summary by section:\n")
  print(table(master$section))
}

cat(sprintf("\nAll outputs saved to: %s\n", normalizePath(OUT_DIR)))


# ==============================================================================
# SECTION 4 — FIGURE: Seasonal Cycle of Precipitation and PM-PET
# ==============================================================================
# SOURCE: Fig2_MS_SPEI123_SPI123.R  (FigS_Seasonality_P_PET.R)
# PURPOSE: Two-panel bar chart: (a) mean monthly P vs PM-PET, (b) P - PET
# DATA:    spi_01_basin_average.csv  +  ERA5Land_Nechako_PET_monthly_summary.csv
# OUTPUT:  FigS_Seasonality_P_PET.pdf + .png  in BASIN_PLOT_DIR
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(data.table)
})

cat("\n============================================================\n")
cat("  SECTION 4: FIGURE — Seasonal Cycle of P and PM-PET\n")
cat("============================================================\n\n")

BASIN_TS_DIR   <- file.path(WD_PATH, "temporal_drought", "basin_averaged_timeseries")
BASIN_PLOT_DIR <- file.path(WD_PATH, "temporal_drought", "basin_averaged_plots")
PET_CSV        <- file.path(WD_PATH, "monthly_data_direct",
                            "ERA5Land_Nechako_PET_monthly_summary.csv")
PRECIP_CSV     <- file.path(BASIN_TS_DIR, "spi_01_basin_average.csv")
dir.create(BASIN_PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

MON_LABELS   <- month.abb
DAYS_PER_MON <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# ── 4.1: Precipitation — 12 monthly means (mm month\u207b\u00b9) ─────────────────────────
if (!file.exists(PRECIP_CSV))
  stop("Precip CSV not found: ", PRECIP_CSV,
       "\n  Run w2_basin_timeseries.R (Part 2) first.")

prec_raw       <- as.data.frame(data.table::fread(PRECIP_CSV))
prec_raw$date  <- as.Date(prec_raw$date)
prec_raw$month <- as.integer(format(prec_raw$date, "%m"))

cm_prec <- prec_raw |>
  dplyr::group_by(month) |>
  dplyr::summarise(mean_val = mean(value, na.rm = TRUE),
                   sd_val   = sd(value,   na.rm = TRUE),
                   n_yrs    = dplyr::n(),
                   .groups  = "drop") |>
  dplyr::mutate(month_lab = factor(MON_LABELS[month], levels = MON_LABELS),
                variable  = "Precipitation (P)")

cat(sprintf("  Precip: %d yrs, mean annual = %.0f mm yr\u207b\u00b9\n",
            cm_prec$n_yrs[1], sum(cm_prec$mean_val)))

# ── 4.2: PM-PET — 12 monthly means (mm month\u207b\u00b9) ─────────────────────────────────
if (!file.exists(PET_CSV))
  stop("PET CSV not found: ", PET_CSV,
       "\n  Run 2preq_PET_ERALand.R first.")

pet_raw       <- as.data.frame(data.table::fread(PET_CSV))
if (!inherits(pet_raw$date, "Date"))
  pet_raw$date <- as.Date(paste0(substr(as.character(pet_raw$date), 1, 7), "-01"))
pet_raw$month  <- as.integer(format(pet_raw$date, "%m"))
pet_raw$pet_mm <- pet_raw$mean_pet * DAYS_PER_MON[pet_raw$month]

cm_pet <- pet_raw |>
  dplyr::group_by(month) |>
  dplyr::summarise(mean_val = mean(pet_mm, na.rm = TRUE),
                   sd_val   = sd(pet_mm,   na.rm = TRUE),
                   n_yrs    = dplyr::n(),
                   .groups  = "drop") |>
  dplyr::mutate(month_lab = factor(MON_LABELS[month], levels = MON_LABELS),
                variable  = "PM-PET")

cat(sprintf("  PM-PET: %d yrs, mean annual = %.0f mm yr\u207b\u00b9\n",
            cm_pet$n_yrs[1], sum(cm_pet$mean_val)))
cat(sprintf("  Aridity index (PET/P) = %.2f\n",
            sum(cm_pet$mean_val) / sum(cm_prec$mean_val)))

# ── 4.3: Panel (a) — grouped bar chart P vs PM-PET ───────────────────────────
COL_P   <- "#4393c3"
COL_PET <- "#d7301f"

cm_all          <- rbind(cm_prec, cm_pet)
cm_all$variable <- factor(cm_all$variable,
                          levels = c("Precipitation (P)", "PM-PET"))
nyrs_label <- cm_prec$n_yrs[1]

p_seas1 <- ggplot2::ggplot(
  cm_all,
  ggplot2::aes(x = month_lab, y = mean_val,
               fill = variable, group = variable)) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(0.72),
    width = 0.65, colour = "white", linewidth = 0.15) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
    position = ggplot2::position_dodge(0.72),
    width = 0.22, colour = "grey30", linewidth = 0.40) +
  ggplot2::scale_fill_manual(
    values = c("Precipitation (P)" = COL_P, "PM-PET" = COL_PET),
    name   = NULL) +
  ggplot2::scale_y_continuous(
    name   = expression(mm~month^{-1}),
    expand = ggplot2::expansion(mult = c(0, 0.09))) +
  ggplot2::labs(
    title    = "(a)  Mean monthly precipitation and Penman\u2013Monteith PET",
    subtitle = paste0(
      nyrs_label, "-year mean (1950\u20132025)  |  Error bars = \u00b11 SD  |  ",
      "Annual totals:  P = ", round(sum(cm_prec$mean_val)), " mm;  ",
      "PM-PET = ", round(sum(cm_pet$mean_val)), " mm"),
    x = NULL) +
  ggplot2::theme_classic(base_size = 9.5) +
  ggplot2::theme(
    legend.position    = c(0.86, 0.88),
    legend.background  = ggplot2::element_rect(fill = "white", colour = NA),
    legend.text        = ggplot2::element_text(size = 8.5),
    legend.key.size    = ggplot2::unit(0.40, "cm"),
    plot.title         = ggplot2::element_text(size = 9.5, face = "bold"),
    plot.subtitle      = ggplot2::element_text(size = 7.5, colour = "grey40"),
    axis.text          = ggplot2::element_text(size = 8.5),
    axis.title.y       = ggplot2::element_text(size = 8.5),
    panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
    plot.margin        = ggplot2::margin(4, 6, 2, 4))

# ── 4.4: Panel (b) — water balance P \u2212 PET ────────────────────────────────────
cm_wb <- dplyr::inner_join(
  cm_prec[, c("month","month_lab","mean_val")],
  cm_pet[,  c("month","mean_val")],
  by = "month", suffix = c("_p","_pet")) |>
  dplyr::mutate(wb       = mean_val_p - mean_val_pet,
                fill_col = ifelse(wb >= 0, COL_P, COL_PET))

p_seas2 <- ggplot2::ggplot(
  cm_wb,
  ggplot2::aes(x = month_lab, y = wb, fill = fill_col)) +
  ggplot2::geom_col(colour = "white", linewidth = 0.15, width = 0.65) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.55) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_y_continuous(
    name   = expression(P - PET~~(mm~month^{-1})),
    expand = ggplot2::expansion(mult = c(0.12, 0.10))) +
  ggplot2::annotate("text",
                    x     = mean(which(cm_wb$wb > 0)),
                    y     = max(cm_wb$wb[cm_wb$wb > 0], na.rm = TRUE) * 0.82,
                    label = "P surplus", colour = COL_P,
                    size = 3.0, fontface = "italic") +
  ggplot2::annotate("text",
                    x     = mean(which(cm_wb$wb < 0)),
                    y     = min(cm_wb$wb[cm_wb$wb < 0], na.rm = TRUE) * 0.78,
                    label = "PET deficit", colour = COL_PET,
                    size = 3.0, fontface = "italic") +
  ggplot2::labs(
    title    = "(b)  Climatic water balance  P \u2212 PM-PET",
    subtitle = "Blue = monthly surplus (P > PET);  red = monthly deficit (PET > P)",
    x        = "Calendar month") +
  ggplot2::theme_classic(base_size = 9.5) +
  ggplot2::theme(
    plot.title         = ggplot2::element_text(size = 9.5, face = "bold"),
    plot.subtitle      = ggplot2::element_text(size = 7.5, colour = "grey40"),
    axis.text          = ggplot2::element_text(size = 8.5),
    axis.title.y       = ggplot2::element_text(size = 8.5),
    panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
    plot.margin        = ggplot2::margin(2, 6, 4, 4))

# ── 4.5: Assemble and save ────────────────────────────────────────────────────
fig_seas <- p_seas1 / p_seas2 +
  patchwork::plot_annotation(
    title    = paste0(
      "Nechako River Basin \u2014 Seasonal cycle of precipitation ",
      "and Penman\u2013Monteith PET"),
    subtitle = paste0(
      "Each bar = mean of ", nyrs_label, " years (1950\u20132025).  ",
      "Annual totals: P = ", round(sum(cm_prec$mean_val)), " mm;  ",
      "PM-PET = ", round(sum(cm_pet$mean_val)), " mm;  ",
      "aridity index = ", round(sum(cm_pet$mean_val) / sum(cm_prec$mean_val), 2), "."),
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey35",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 4))))

cat("\n\u2500\u2500 Saving FigS_Seasonality_P_PET...\n")
for (ext in c("pdf", "png")) {
  out_f <- file.path(BASIN_PLOT_DIR, paste0("FigS_Seasonality_P_PET.", ext))
  tryCatch(
    ggplot2::ggsave(out_f, fig_seas,
                    width  = 7.2,
                    height = 7.0,
                    units  = "in",
                    dpi    = if (ext == "png") 300 else "print",
                    device = ext),
    error = function(e) cat(sprintf("  \u26a0 %s: %s\n", ext, e$message)))
  cat(sprintf("  \u2713 %s\n", basename(out_f)))
}
cat(sprintf("\n\u2713 Section 4 complete.\n  Output: %s\n", normalizePath(BASIN_PLOT_DIR)))

##############################################################################
# END OF THE SCRIPT
##############################################################################