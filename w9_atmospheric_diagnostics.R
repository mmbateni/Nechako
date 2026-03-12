# ============================================================================
#  w9_atmospheric_diagnostics.R
#  Nechako Watershed Drought — Atmospheric Pattern Diagnostics
#
#  Index coverage:
#    ALL available SPI and SPEI timescales in the project pipeline
#    (SPI 1, 3, 6, 12 and SPEI 1, 3, 6, 12 by default; scale 9 loaded if
#    files exist).  There is no a priori reason to restrict atmospheric
#    compositing to SPEI-6 alone — ENSO and PDO act on timescales from
#    months (SPI-1/3 for precipitation response) to seasons and years
#    (SPEI-6/12 for cumulative water deficit).  Running all timescales
#    allows comparison of which atmospheric patterns are robustly
#    associated with drought regardless of scale, and which only emerge
#    at longer accumulation periods.
#
#  Output structure:
#    {OUT_DIR}/
#      Tables/
#        T0_drought_months_all_indices.csv   <-- all indices, all drought months
#        {INDEX}/T1_drought_month_summary.csv
#        {INDEX}/T2_seasonal_anomaly_statistics.csv
#        {INDEX}/T3_correlation_matrix.csv
#        {INDEX}/T3_pvalue_matrix.csv
#      Figures/
#        SST/
#          09_sst_annual_anomaly_recent.png  (index-independent)
#          11_sst_timeseries_ne_pacific.png  (index-independent)
#        {INDEX}/
#          01_z500_monthly_anomaly_composites.png
#          02_z500_drought_vs_nondrought.png
#          04_z500_timeseries_ridge.png
#          05_slp_monthly_anomaly_composites.png
#          06_slp_drought_vs_nondrought.png
#          08_slp_timeseries.png
#          10_sst_monthly_composites_drought.png
#          12_joint_scatter.png
#          13_composite_summary_panel.png
# ============================================================================

# ── Libraries ────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
  library(lubridate)
  library(ncdf4)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(zoo)
})

# ============================================================================
#  CONFIGURATION
# ============================================================================
WD_PATH  <- "D:/Nechako_Drought/Nechako/"
DATA_DIR <- "D:/Nechako_Drought/Nechako/monthly_data_direct"
OUT_DIR  <- file.path(WD_PATH, "atmospheric_diagnostics")

setwd(WD_PATH)

# ── Index CSV directories (must match DROUGHT_ANALYSIS_utils.R) ──────────────
SPI_SEAS_DIR  <- file.path(WD_PATH, "spi_results_seasonal/")
SPEI_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal/")

# Scales to attempt.  Scale 9 is included; skipped gracefully if no file found.
SPI_SCALES_TRY  <- c(1, 3, 6, 9, 12)
SPEI_SCALES_TRY <- c(1, 3, 6, 9, 12)

# ── Analysis period ──────────────────────────────────────────────────────────
START_YEAR <- 1950
END_YEAR   <- 2025

# ── Climatology reference ────────────────────────────────────────────────────
CLIM_START <- 1991
CLIM_END   <- 2020

# ── Drought definitions ──────────────────────────────────────────────────────
# Two definitions are run in parallel; every output (table + figure) is
# produced twice — once per definition — under separate sub-directories.
#
# [1] Sheffield & Wood (2008) — SINGLE threshold.
#     Drought begins AND ends when the index crosses the same value.
#     q0 = 10th percentile of a standard normal = qnorm(0.10) ≈ -1.282.
#     No hysteresis, no minimum duration filter.
#
# [2] 2-Threshold / Alt — HYSTERESIS + minimum duration (matching
#     ALT_DROUGHT_DEF in w3_trend_test.R).
#     Onset: index < -1.0.  Termination: index >= -0.5.
#     Scale-specific minimum duration filters out sub-threshold blips.
DROUGHT_DEFS <- list(
  
  SW = list(
    id          = "SW",
    short_label = "SW_q10",
    long_label  = "Sheffield & Wood (2008) \u2014 Single threshold (q\u2080=10%, \u2248-1.28)",
    classify    = "single",
    threshold   = qnorm(0.10)          # -1.2816 for a standard-normal index
  ),
  
  Alt2T = list(
    id          = "Alt2T",
    short_label = "Alt2T",
    long_label  = "2-Threshold (Onset \u2264-1.0 / Termination \u2265-0.5, min duration by scale)",
    classify    = "hysteresis",
    onset       = -1.0,
    termination = -0.5,
    # Minimum event duration (months) — matches ALT_DROUGHT_DEF in w3_trend_test.R
    min_duration_map = list("1" = 2L, "3" = 3L, "6" = 4L, "12" = 6L),
    default_min = 3L
  )
)

# ── Recent period for SST trend analysis ─────────────────────────────────────
# 2018 marks the reintensification of the NE Pacific marine heatwave that
# preceded the 2021 BC heat dome.  Extend back to 2014 to include the
# original "Blob" (Bond et al. 2015) if desired.
RECENT_START <- 2018
RECENT_END   <- END_YEAR

# ── NW Canada ridge averaging box ────────────────────────────────────────────
# Where a blocking ridge would appear in Z500/SLP during warm-dry BC summers.
# Covers the Gulf of Alaska action centre of the PNA pattern (~55 N, 120 W).
RIDGE_LON_MIN <- -140; RIDGE_LON_MAX <- -110
RIDGE_LAT_MIN <-   50; RIDGE_LAT_MAX <-   65

# ── NE Pacific SST averaging box ─────────────────────────────────────────────
# Covers the PDO "warm tongue" (~35-50 N, 150-130 W) — the region where NE
# Pacific SSTs most strongly co-vary with the PDO index and marine heatwave
# events.  The tropical Pacific is excluded to avoid aliasing ENSO SST
# variance; ENSO forcing is captured via ONI in the teleconnection arm (w7-w8).
SST_MEAN_LON_MIN <- -160; SST_MEAN_LON_MAX <- -130
SST_MEAN_LAT_MIN <-   40; SST_MEAN_LAT_MAX <-   55

# ── Plot settings ─────────────────────────────────────────────────────────────
ANOMALY_PALETTE   <- "RdBu"
SST_PALETTE       <- "RdBu"
FIGURE_DPI        <- 200
FIGURE_WIDTH_WIDE <- 16
FIGURE_WIDTH_STD  <- 12
MONTH_LABELS      <- month.abb

# ============================================================================
#  HELPER FUNCTIONS
# ============================================================================

get_map_layers <- function() {
  borders <- ne_countries(scale = "medium", returnclass = "sf")
  lakes   <- ne_download(scale = "medium", type = "lakes",
                         category = "physical", returnclass = "sf")
  list(borders = borders, lakes = lakes)
}

# Load one index CSV (spi or spei at a given scale).Returns data.frame(date, value) or NULL.
load_index_csv <- function(index_type, scale) {
  dir_map <- list(spi = SPI_SEAS_DIR, spei = SPEI_SEAS_DIR)
  dir     <- dir_map[[tolower(index_type)]]
  f       <- file.path(dir, sprintf("%s_%02d_basin_averaged_by_month.csv",
                                    tolower(index_type), as.integer(scale)))
  if (!file.exists(f)) return(NULL)
  
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  
  # FIX: Handle "Year" column name variations
  year_col <- grep("^year$", names(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(year_col) || !length(year_col)) return(NULL)
  names(df)[names(df) == year_col] <- "Year"
  
  if (!"Year" %in% names(df)) return(NULL)
  df$Year <- as.integer(df$Year)
  
  # FIX: Check BOTH full and abbreviated month names
  month_names_full <- c("January","February","March","April","May","June",
                        "July","August","September","October","November","December")
  month_names_abb  <- c("Jan","Feb","Mar","Apr","May","Jun",
                        "Jul","Aug","Sep","Oct","Nov","Dec")
  mon_cols <- intersect(month_names_full, names(df))
  if (!length(mon_cols)) mon_cols <- intersect(month_names_abb, names(df))
  if (!length(mon_cols)) return(NULL)
  
  # Reshape wide (76 rows × 12 cols) to long (912 rows × 1 col)
  long <- do.call(rbind, lapply(seq_along(mon_cols), function(mi) {
    data.frame(date  = as.Date(paste(df$Year, sprintf("%02d", mi), "01", sep = "-")),
               value = as.numeric(df[[mon_cols[mi]]]),
               stringsAsFactors = FALSE)
  }))
  
  long <- long[order(long$date), ]
  long <- long[!is.na(long$value) & is.finite(long$value), ]
  rownames(long) <- NULL
  
  # Validate
  expected_rows <- length(unique(df$Year)) * 12
  if (nrow(long) != expected_rows) {
    warning(sprintf("Index %s-%d: expected %d rows, got %d",
                    toupper(index_type), scale, expected_rows, nrow(long)))
  }
  
  long
}
# Like build_monthly_composite but uses the atmospheric raster layer lag_k
# months BEFORE each qualifying drought month.  For lag_k = 0 this is
# identical to build_monthly_composite.
build_monthly_composite_lag <- function(rast_obj, date_idx, drought_mask_vec, lag_k = 0) {
  n_layers <- nlyr(rast_obj)
  out_list <- vector("list", 12)
  for (m in 1:12) {
    # Layer positions of drought months in calendar month m
    dr_pos <- which(date_idx$month == m & drought_mask_vec)
    # Shift back by lag_k months
    idx    <- dr_pos - lag_k
    # Keep only indices that fall inside the raster
    idx    <- idx[idx >= 1L & idx <= n_layers]
    n      <- length(idx)
    mean_layer <- if (n > 0) {
      app(subset(rast_obj, idx), mean, na.rm = TRUE)
    } else {
      rast_obj[[1L]] * NA
    }
    df <- as.data.frame(mean_layer, xy = TRUE)
    names(df)[3] <- "value"
    if (nrow(df) > 0) {
      df$month_label <- MONTH_LABELS[m]
      df$panel_title <- sprintf("%s (n=%d)", MONTH_LABELS[m], n)
    } else {
      df$month_label <- character(0)
      df$panel_title <- character(0)
    }
    out_list[[m]] <- df
  }
  bind_rows(out_list) %>%
    mutate(
      month_label = factor(month_label, levels = MONTH_LABELS),
      panel_title = factor(panel_title,
                           levels = sprintf("%s (n=%d)", MONTH_LABELS,
                                            sapply(1:12, function(m) {
                                              pos <- which(date_idx$month == m &
                                                             drought_mask_vec) - lag_k
                                              sum(pos >= 1L & pos <= n_layers)
                                            })))
    )
}
# Build 12-panel composite: ALWAYS all 12 calendar months.Months with no qualifying rows:all-NA layer.renders grey with n=0 label.
build_monthly_composite <- function(rast_obj, date_idx, qualifying_mask) {
  n_layers <- nlyr(rast_obj)
  out_list <- vector("list", 12)
  for (m in 1:12) {
    idx <- which(date_idx$month == m & qualifying_mask)
    n   <- length(idx)
    mean_layer <- if (n > 0) {
      app(subset(rast_obj, idx), mean, na.rm = TRUE)
    } else {
      rast_obj[[1]] * NA
    }
    df <- as.data.frame(mean_layer, xy = TRUE)
    names(df)[3] <- "value"
    # ── FIX: Handle empty data frame (n=0 months) ──────────────────
    if (nrow(df) > 0) {
      df$month_label <- MONTH_LABELS[m]
      df$panel_title <- sprintf("%s (n=%d)", MONTH_LABELS[m], n)
    } else {
      df$month_label <- character(0)
      df$panel_title <- character(0)
    }
    # ───────────────────────────────────────────────────────────────
    out_list[[m]] <- df
  }
  bind_rows(out_list) %>%
    mutate(month_label = factor(month_label, levels = MONTH_LABELS),
           panel_title = factor(panel_title,
                                levels = sprintf("%s (n=%d)",
                                                 MONTH_LABELS,
                                                 sapply(1:12, function(m)
                                                   sum(date_idx$month == m &
                                                         qualifying_mask)))))
}
map_theme <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(strip.background = element_rect(fill = "grey92"),
          strip.text       = element_text(face = "bold"),
          legend.position  = "right",
          plot.title       = element_text(face = "bold", size = base_size + 2),
          plot.subtitle    = element_text(size = base_size - 1),
          panel.spacing    = unit(0.3, "lines"),
          axis.text        = element_blank(),
          axis.ticks       = element_blank())
}

ensure_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# ── Drought month classifiers ─────────────────────────────────────────────────

# Extract numeric scale from an index label such as "SPI3", "SPEI12".
extract_scale <- function(lbl) {
  as.integer(regmatches(lbl, regexpr("[0-9]+$", lbl)))
}

# Sheffield & Wood (2008) — single threshold.
# Returns a logical vector: TRUE for every month where value <= threshold.
# Entry and exit are governed by the identical threshold value.
classify_drought_sw <- function(values, threshold) {
  !is.na(values) & values <= threshold
}

# 2-Threshold / Alt definition — hysteresis + minimum duration.
# Mirrors vectorized_event_detection_alt() in w3_trend_test.R.
# Returns a logical vector: TRUE for every month that falls INSIDE a
# qualifying drought event (onset < onset_thr; termination >= term_thr;
# event length >= min_dur months).
classify_drought_alt <- function(values, onset_thr, term_thr, min_dur) {
  n      <- length(values)
  is_dr  <- logical(n)
  in_d   <- FALSE
  s_idx  <- NA_integer_
  
  for (j in seq_len(n)) {
    v <- values[j]
    if (is.na(v)) next
    if (!in_d && v < onset_thr) {
      in_d  <- TRUE
      s_idx <- j
    } else if (in_d && v >= term_thr) {
      e_idx    <- j - 1L
      duration <- e_idx - s_idx + 1L
      if (duration >= min_dur) is_dr[s_idx:e_idx] <- TRUE
      in_d  <- FALSE
      s_idx <- NA_integer_
    }
  }
  # Close any event still open at end of record
  if (in_d && !is.na(s_idx)) {
    duration <- n - s_idx + 1L
    if (duration >= min_dur) is_dr[s_idx:n] <- TRUE
  }
  is_dr
}

# Dispatcher: classify drought months according to the active definition.
# Returns a logical vector aligned to `values`.
classify_drought <- function(values, def, index_label) {
  if (def$classify == "single") {
    classify_drought_sw(values, def$threshold)
  } else {
    sc    <- extract_scale(index_label)
    sc_ch <- as.character(sc)
    mndur <- if (!is.null(def$min_duration_map[[sc_ch]])) {
      def$min_duration_map[[sc_ch]]
    } else {
      def$default_min
    }
    classify_drought_alt(values, def$onset, def$termination, mndur)
  }
}

# Human-readable threshold annotation for plot subtitles / captions.
def_thr_label <- function(def, lbl = NULL) {
  if (def$classify == "single") {
    sprintf("%s \u2264 %.3f", if (!is.null(lbl)) lbl else "Index", def$threshold)
  } else {
    sc    <- if (!is.null(lbl)) extract_scale(lbl) else NA
    sc_ch <- as.character(sc)
    mndur <- if (!is.na(sc) && !is.null(def$min_duration_map[[sc_ch]])) {
      def$min_duration_map[[sc_ch]]
    } else {
      def$default_min
    }
    sprintf("Onset < %.1f / Term \u2265 %.1f / min %d months",
            def$onset, def$termination, mndur)
  }
}

# ============================================================================
#  SECTION 1: LOAD ERA5 FIELDS
# ============================================================================
cat("============================================================\n")
cat(" SECTION 1: Loading ERA5 atmospheric fields\n")
cat("============================================================\n")

cat("  Loading Z500 anomaly (1991-2020 climatology)... ")
z500_anom <- rast(file.path(DATA_DIR,
                            sprintf("z500_monthly_anomaly_clim%d_%d.nc", CLIM_START, CLIM_END)))
cat(sprintf("Done. (%d layers)\n", nlyr(z500_anom)))

cat("  Loading SLP (raw hPa) and computing 1991-2020 anomaly... ")
slp_raw <- rast(file.path(DATA_DIR, "slp_monthly.nc"))
slp_anom <- local({
  n_layers  <- nlyr(slp_raw)
  clim_mask <- (START_YEAR:END_YEAR) >= CLIM_START & (START_YEAR:END_YEAR) <= CLIM_END
  anom_list <- vector("list", 12)
  for (m in 1:12) {
    all_idx  <- seq(m, n_layers, by = 12)
    clim_idx <- all_idx[clim_mask]
    clim     <- app(subset(slp_raw, clim_idx), mean, na.rm = TRUE)
    anom_list[[m]] <- subset(slp_raw, all_idx) - clim
  }
  combined    <- do.call(c, anom_list)
  layer_order <- order(c(sapply(1:12, function(m) seq(m, n_layers, by = 12))))
  combined[[layer_order]]
})
cat(sprintf("Done. (%d layers)\n", nlyr(slp_anom)))

cat("  Loading SST anomaly (1991-2020 climatology)... ")
sst_anom <- rast(file.path(DATA_DIR,
                           sprintf("sst_monthly_anomaly_clim%d_%d.nc", CLIM_START, CLIM_END)))
cat(sprintf("Done. (%d layers)\n", nlyr(sst_anom)))

# ── Date spine (raster layer i = calendar month i) ───────────────────────────
all_dates <- seq.Date(as.Date(sprintf("%d-01-01", START_YEAR)),
                      as.Date(sprintf("%d-12-01", END_YEAR)),
                      by = "month")
n_total <- length(all_dates)
if (nlyr(z500_anom) != n_total)
  warning(sprintf("Z500 layers (%d) != expected (%d).", nlyr(z500_anom), n_total))

base_date_index <- tibble(
  layer  = seq_len(n_total),
  date   = all_dates,
  year   = year(all_dates),
  month  = month(all_dates),
  season = case_when(
    month(all_dates) %in% c(12, 1, 2)  ~ "DJF",
    month(all_dates) %in% c(3, 4, 5)   ~ "MAM",
    month(all_dates) %in% c(6, 7, 8)   ~ "JJA",
    month(all_dates) %in% c(9, 10, 11) ~ "SON"
  )
)

# ── Area-averaged atmospheric time series (computed once, reused every index) ─
ridge_ext <- ext(RIDGE_LON_MIN, RIDGE_LON_MAX, RIDGE_LAT_MIN, RIDGE_LAT_MAX)
sst_ext   <- ext(SST_MEAN_LON_MIN, SST_MEAN_LON_MAX, SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX)

z500_ridge_ts <- global(crop(z500_anom, ridge_ext), "mean", na.rm = TRUE)[, 1]
slp_nwbc_ts   <- global(crop(slp_anom,  ridge_ext), "mean", na.rm = TRUE)[, 1]
sst_nepac_ts  <- global(crop(sst_anom,  sst_ext),   "mean", na.rm = TRUE)[, 1]
cat("  Area-averaged time series extracted.\n\n")

cat("  Loading map layers...\n")
map_layers <- suppressMessages(suppressWarnings(get_map_layers()))

# ============================================================================
#  SECTION 2: LOAD ALL DROUGHT INDICES
# ============================================================================
cat("============================================================\n")
cat(" SECTION 2: Loading SPI / SPEI indices\n")
cat("============================================================\n")

# ── FIX: Updated load_index_csv() to handle both full and abbreviated month names ──
load_index_csv <- function(index_type, scale) {
  dir_map <- list(spi = SPI_SEAS_DIR, spei = SPEI_SEAS_DIR)
  dir <- dir_map[[tolower(index_type)]]
  f <- file.path(dir, sprintf("%s_%02d_basin_averaged_by_month.csv",
                              tolower(index_type), as.integer(scale)))
  if (!file.exists(f)) return(NULL)
  
  # Try to read with data.table::fread (handles many formats)
  df <- data.table::fread(f, data.table = FALSE)
  
  # Check number of columns
  if (ncol(df) == 2 && names(df)[1] == "date" && names(df)[2] == "value") {
    # Already long format
    long <- df
    names(long) <- c("date", "value")
    long$date <- as.Date(long$date)
  } else if (ncol(df) >= 13) {
    # Wide format: assume first column is Year, next 12 are months
    names(df)[1] <- "Year"
    df$Year <- as.integer(df$Year)
    month_cols <- names(df)[2:13]
    long <- do.call(rbind, lapply(seq_along(month_cols), function(mi) {
      data.frame(
        date  = as.Date(paste(df$Year, sprintf("%02d", mi), "01", sep = "-")),
        value = as.numeric(df[[month_cols[mi]]]),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    stop("Unexpected CSV format in ", f)
  }
  
  long <- long[order(long$date), ]
  long <- long[!is.na(long$value) & is.finite(long$value), ]
  rownames(long) <- NULL
  long
}

# ── Load all indices ──────────────────────────────────────────────────────────
all_indices <- list()
for (sc in SPI_SCALES_TRY) {
  df <- load_index_csv("spi", sc)
  lbl <- sprintf("SPI%d", sc)
  if (!is.null(df)) {
    all_indices[[lbl]] <- df
    cat(sprintf("  + %s: %d months (%d years × 12)\n", lbl, nrow(df), length(unique(df$Year))))
  } else {
    cat(sprintf("  - %s: file not found, skipped\n", lbl))
  }
}

for (sc in SPEI_SCALES_TRY) {
  df <- load_index_csv("spei", sc)
  lbl <- sprintf("SPEI%d", sc)
  if (!is.null(df)) {
    all_indices[[lbl]] <- df
    cat(sprintf("  + %s: %d months (%d years × 12)\n", lbl, nrow(df), length(unique(df$Year))))
  } else {
    cat(sprintf("  - %s: file not found, skipped\n", lbl))
  }
}

if (length(all_indices) == 0)
  stop("No index CSV files found. Check SPI_SEAS_DIR and SPEI_SEAS_DIR.")

# ── Validation: Confirm all indices have 912 rows (76 years × 12 months) ─────
cat("\n  Validation:\n")
for (lbl in names(all_indices)) {
  n <- nrow(all_indices[[lbl]])
  status <- if (n == 912) "✓" else "⚠"
  cat(sprintf("    %s %s: %d rows (expected 912)\n", status, lbl, n))
}

cat(sprintf("\n  Loaded: %s\n\n", paste(names(all_indices), collapse = ", ")))

# ============================================================================
#  SECTION 3: TABLE T0 — EXACT DROUGHT MONTHS FOR ALL INDICES
#  Generated separately for each drought definition.
# ============================================================================
cat("============================================================\n")
cat(" SECTION 3: Drought month tables (both definitions)\n")
cat("============================================================\n")

for (def in DROUGHT_DEFS) {
  
  cat(sprintf("\n  -- Definition: %s --\n", def$long_label))
  
  def_tbl_root <- ensure_dir(OUT_DIR, def$short_label, "Tables")
  
  drought_rows <- lapply(names(all_indices), function(lbl) {
    df_temp <- all_indices[[lbl]]
    df_temp$yr  <- as.integer(format(df_temp$date, "%Y"))
    df_temp$mon <- as.integer(format(df_temp$date, "%m"))
    
    # Classify using the active definition
    drought_flag <- classify_drought(df_temp$value, def, lbl)
    
    df_temp %>%
      mutate(is_drought_def = drought_flag) %>%
      filter(yr >= START_YEAR,
             yr <= END_YEAR,
             is_drought_def,
             !is.na(value)) %>%
      transmute(
        index          = lbl,
        definition     = def$short_label,
        date           = format(date, "%Y-%m"),
        year           = yr,
        month_num      = mon,
        month_name     = month.abb[mon],
        season         = case_when(
          mon %in% c(12, 1, 2)  ~ "DJF",
          mon %in% c(3, 4, 5)   ~ "MAM",
          mon %in% c(6, 7, 8)   ~ "JJA",
          mon %in% c(9, 10, 11) ~ "SON"
        ),
        index_value    = round(value, 3),
        threshold_used = def_thr_label(def, lbl)
      )
  })
  
  t0 <- bind_rows(drought_rows) %>% arrange(index, date)
  write_csv(t0, file.path(def_tbl_root, "T0_drought_months_all_indices.csv"))
  
  cat(sprintf("  Drought month counts [%s]:\n", def$short_label))
  t0 %>%
    group_by(index) %>%
    summarise(n_drought = n(),
              pct = round(100 * n() / n_total, 1),
              .groups = "drop") %>%
    as.data.frame() %>%
    print()
  cat(sprintf("  Saved: %s/Tables/T0_drought_months_all_indices.csv\n", def$short_label))
  
}  # end definition loop (T0)

# ============================================================================
#  SECTION 4: SST FIGURES (index-independent, generated once)
# ============================================================================
cat("============================================================\n")
cat(" SECTION 4: SST figures (not index-specific)\n")
cat("============================================================\n")
ensure_dir(OUT_DIR, "Figures", "SST")

# ── Fig 09 ────────────────────────────────────────────────────────────────────
cat("  Fig 09: Annual SST anomaly maps (recent years)... ")
recent_years <- RECENT_START:min(RECENT_END, END_YEAR)
sst_annual_list <- lapply(recent_years, function(yr) {
  idx <- which(base_date_index$year == yr)
  if (!length(idx)) return(NULL)
  df  <- as.data.frame(app(subset(sst_anom, idx), mean, na.rm = TRUE), xy = TRUE)
  names(df)[3] <- "value"
  df$year <- yr; df
})
sst_annual_df <- bind_rows(sst_annual_list) %>% mutate(year = factor(year))
sst_lim <- ceiling(max(abs(quantile(sst_annual_df$value, c(0.02, 0.98), na.rm = TRUE))) * 2) / 2

p09 <- ggplot(sst_annual_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = "grey85", color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = SST_PALETTE,
                       limits = c(-sst_lim, sst_lim), oob = squish,
                       name = "SST Anom\n(deg C)", na.value = "grey80") +
  coord_sf(xlim = c(-180, -110), ylim = c(20, 65), expand = FALSE) +
  facet_wrap(~year, ncol = 4) +
  labs(
    title    = sprintf("NE Pacific SST Anomaly, Annual Mean %d-%d",
                       RECENT_START, min(RECENT_END, END_YEAR)),
    subtitle = sprintf("Reference climatology: %d-%d (WMO standard)",
                       CLIM_START, CLIM_END),
    caption  = paste0(
      "Domain: ", abs(SST_MEAN_LON_MIN), "-", abs(SST_MEAN_LON_MAX),
      " W, ", SST_MEAN_LAT_MIN, "-", SST_MEAN_LAT_MAX, " N.\n",
      "This region captures the PDO 'warm tongue', the area in the NE Pacific where SST anomalies ",
      "most strongly co-vary with the Pacific Decadal Oscillation (PDO) index. A warm PDO phase ",
      "(sustained positive/red anomalies here) is associated with atmospheric ridging over BC that ",
      "suppresses precipitation and enhances evapotranspiration, favouring drought. ",
      "The 2014-2016 'Blob' (Bond et al. 2015) and the 2019-2021 NE Pacific marine heatwave ",
      "both appear as multi-year positive anomalies in this domain. ",
      "The tropical Pacific is excluded: ENSO SST influence is captured via ONI in the ",
      "teleconnection analysis (w7-w8), avoiding double-counting."
    ),
    x = NULL, y = NULL
  ) +
  map_theme(10)

ggsave(file.path(OUT_DIR, "Figures", "SST", "09_sst_annual_anomaly_recent.png"),
       p09,
       width  = FIGURE_WIDTH_WIDE,
       height = ceiling(length(recent_years) / 4) * 4 + 2,
       dpi    = FIGURE_DPI)
cat("Saved.\n")

# ── Fig 11 ────────────────────────────────────────────────────────────────────
cat("  Fig 11: SST NE Pacific long time series... ")
sst_ts <- base_date_index %>%
  mutate(sst_nepac = sst_nepac_ts) %>%
  filter(!is.na(sst_nepac)) %>%
  arrange(date) %>%
  mutate(sst_12mo = rollmean(sst_nepac, k = 12, fill = NA, align = "right"))

p11 <- ggplot(sst_ts, aes(x = date, y = sst_nepac)) +
  annotate("rect",
           xmin = as.Date(sprintf("%d-01-01", RECENT_START)),
           xmax = as.Date(sprintf("%d-12-31", min(RECENT_END, END_YEAR))),
           ymin = -Inf, ymax = Inf, fill = "lightyellow", alpha = 0.5) +
  geom_col(aes(fill = sst_nepac > 0), alpha = 0.4, width = 25) +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                    guide = "none") +
  geom_line(aes(y = sst_12mo), color = "black", linewidth = 0.8, na.rm = TRUE) +
  geom_hline(yintercept =  0.0, color = "grey30", linewidth = 0.4) +
  geom_hline(yintercept =  0.5, linetype = "dotted", color = "darkred",   linewidth = 0.5) +
  geom_hline(yintercept = -0.5, linetype = "dotted", color = "steelblue", linewidth = 0.5) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y", expand = c(0.01, 0)) +
  labs(
    title    = "NE Pacific SST Anomaly — Monthly Time Series",
    subtitle = sprintf("Area average: %d-%d W, %d-%d N  |  Black line: 12-month rolling mean  |  Ref: %d-%d",
                       abs(SST_MEAN_LON_MIN), abs(SST_MEAN_LON_MAX),
                       SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX, CLIM_START, CLIM_END),
    x = NULL, y = "SST Anomaly (deg C)",
    caption = paste0(
      "Dotted lines: +/-0.5 deg C (approximate warm/cool PDO phase boundary).\n",
      "Yellow shading: recent period (", RECENT_START, " onward)."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Figures", "SST", "11_sst_timeseries_ne_pacific.png"),
       p11, width = FIGURE_WIDTH_WIDE, height = 5.5, dpi = FIGURE_DPI)
cat("Saved.\n\n")

# ============================================================================
#  SECTION 5: PER-INDEX LOOP  ×  DROUGHT DEFINITION LOOP
#  Every figure and table is produced twice — once per definition.
#  Output structure:
#    {OUT_DIR}/{def$short_label}/Figures/{INDEX}/  <-- figures
#    {OUT_DIR}/{def$short_label}/Tables/{INDEX}/   <-- tables T1-T3
# ============================================================================
cat("============================================================\n")
cat(" SECTION 5: Per-index composites and tables (both definitions)\n")
cat("============================================================\n\n")

for (def in DROUGHT_DEFS) {
  
  cat(sprintf("\n╔══════════════════════════════════════════════════════╗\n"))
  cat(sprintf("  Definition: %s\n", def$long_label))
  cat(sprintf("╚══════════════════════════════════════════════════════╝\n\n"))
  
  # Root output directories for this definition
  def_fig_root <- ensure_dir(OUT_DIR, def$short_label, "Figures")
  def_tbl_root <- ensure_dir(OUT_DIR, def$short_label, "Tables")
  
  # Collector: best-lag rows for every index under this definition.
  # Written to T6 after all indices are processed.
  def_best_lag_rows <- list()

  for (lbl in names(all_indices)) {
    
    cat(sprintf(">> %s  [%s]\n", lbl, def$short_label))
    
    fig_dir <- ensure_dir(def_fig_root, lbl)
    tbl_dir <- ensure_dir(def_tbl_root, lbl)
    
    # ── Attach index to date spine ───────────────────────────────────────────────
    idx_df <- all_indices[[lbl]] %>%
      filter(as.integer(format(date, "%Y")) >= START_YEAR,
             as.integer(format(date, "%Y")) <= END_YEAR) %>%
      rename(index_value = value)
    
    # Classify drought months according to the active definition
    idx_df <- idx_df %>%
      mutate(is_drought = classify_drought(index_value, def, lbl))
    
    # Human-readable threshold label for plot annotations
    thr_lbl <- def_thr_label(def, lbl)
    
    date_index <- base_date_index %>%
      left_join(idx_df %>% select(date, index_value, is_drought), by = "date") %>%
      mutate(z500_ridge = z500_ridge_ts,
             slp_nwbc   = slp_nwbc_ts,
             sst_nepac  = sst_nepac_ts)
    
    drought_mask   <- date_index$is_drought %in% TRUE
    nondrought_mask <- date_index$is_drought %in% FALSE
    drought_idx    <- which(drought_mask)
    nondrought_idx <- which(nondrought_mask)
    n_dr  <- length(drought_idx)
    n_ndr <- length(nondrought_idx)
    idx_range <- max(abs(date_index$index_value), na.rm = TRUE)
    
    cat(sprintf("   Drought months [%s]: %d  |  Non-drought: %d\n",
                thr_lbl, n_dr, n_ndr))
    
    if (n_dr == 0) {
      cat("   No drought months found — skipping.\n\n"); next
    }
    
    # ── Base drought classification (concurrent — index at time t) ───────────────
    drought_mask_base    <- date_index$is_drought %in% TRUE
    nondrought_mask_base <- date_index$is_drought %in% FALSE
    drought_pos          <- which(drought_mask_base)
    nondrought_pos       <- which(nondrought_mask_base)
    n_dr_base            <- length(drought_pos)
    idx_range            <- max(abs(date_index$index_value), na.rm = TRUE)

    # ── Number of lags to investigate = time-scale of the index ─────────────────
    # For SPI3 / SPEI3: lags 0, 1, 2  (3 investigations)
    # For SPI12        : lags 0 … 11  (12 investigations)
    n_scale  <- extract_scale(lbl)          # e.g. 3 for SPI3, 12 for SPI12
    lag_seq  <- 0L:(n_scale - 1L)

    cat(sprintf("   Drought months [%s]: %d  |  Non-drought: %d  |  Lags to run: 0..%d\n",
                thr_lbl, n_dr_base, length(nondrought_pos), n_scale - 1L))

    if (n_dr_base == 0) {
      cat("   No drought months found — skipping.\n\n"); next
    }

    # ── plot_dual_panel_ts is defined once here; unchanged from original ─────────
    plot_dual_panel_ts <- function(dates, top_ts, bot_ts, drought_mask,
                                   top_label, bot_label,
                                   title_txt, subtitle_txt,
                                   def, thr_lbl,
                                   index_label = NULL) {

      roll_k <- 12

      loess_span_map <- c("1"  = 0.10,
                          "3"  = 0.08,
                          "6"  = 0.07,
                          "9"  = 0.06,
                          "12" = 0.05)

      sc_key   <- as.character(as.integer(regmatches(
        index_label, regexpr("[0-9]+$", index_label))))
      bot_span <- if (!is.null(index_label) && sc_key %in% names(loess_span_map)) {
        loess_span_map[[sc_key]]
      } else {
        0.08
      }

      df_raw <- tibble(
        date       = dates,
        top        = top_ts,
        bot        = bot_ts,
        is_drought = drought_mask
      ) %>%
        filter(!is.na(top), !is.na(bot)) %>%
        arrange(date)

      loess_fit <- loess(
        bot ~ as.numeric(date),
        data      = df_raw,
        span      = bot_span,
        degree    = 2,
        na.action = na.exclude
      )

      df <- df_raw %>%
        mutate(
          top_smooth = rollmean(top, k = roll_k, fill = NA, align = "center"),
          bot_smooth = predict(loess_fit, newdata = data.frame(date = as.numeric(date)))
        )

      drought_runs <- df %>%
        mutate(grp = cumsum(c(TRUE, diff(as.integer(is_drought)) != 0))) %>%
        group_by(grp) %>%
        summarise(xmin = min(date), xmax = max(date),
                  drought = first(is_drought), .groups = "drop") %>%
        filter(drought)

      shade <- geom_rect(
        data        = drought_runs,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill        = "#f4a582",
        alpha       = 0.25,
        inherit.aes = FALSE
      )

      x_scale <- scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                              expand = c(0.01, 0))
      base_theme <- theme_minimal(base_size = 11) +
        theme(
          plot.title       = element_text(face = "bold", size = 12),
          plot.subtitle    = element_text(size = 9, color = "grey30"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey88"),
          axis.title.y     = element_text(size = 9),
          axis.text.x      = element_text(size = 9),
          plot.margin      = margin(4, 8, 2, 8)
        )

      top_lim <- ceiling(max(abs(df$top), na.rm = TRUE) / 10) * 10

      p_top <- ggplot(df, aes(x = date)) +
        shade +
        geom_ribbon(aes(ymin = 0, ymax = top),
                    fill = "#92c5de", alpha = 0.55) +
        geom_ribbon(aes(ymin = pmin(top, 0), ymax = 0),
                    fill = "#4393c3", alpha = 0.35) +
        geom_ribbon(aes(ymin = 0, ymax = pmax(top, 0)),
                    fill = "#d6604d", alpha = 0.55) +
        geom_line(aes(y = top_smooth), color = "grey15",
                  linewidth = 0.9, na.rm = TRUE) +
        geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
        x_scale +
        scale_y_continuous(name = top_label,
                           limits = c(-top_lim, top_lim)) +
        labs(x = NULL) +
        base_theme

      onset_thr <- if (def$classify == "single") def$threshold else def$onset

      p_bot <- ggplot(df, aes(x = date)) +
        shade +
        geom_ribbon(aes(ymin = pmin(bot, 0), ymax = 0),
                    fill = "#d6604d", alpha = 0.4) +
        geom_ribbon(aes(ymin = 0, ymax = pmax(bot, 0)),
                    fill = "#92c5de", alpha = 0.4) +
        geom_line(aes(y = bot_smooth), color = "grey15",
                  linewidth = 0.85, na.rm = TRUE) +
        geom_hline(yintercept = onset_thr,
                   linetype = "dashed", color = "#b2182b", linewidth = 0.55) +
        geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
        annotate("text", x = min(df$date), y = onset_thr + 0.08,
                 label = sprintf("Onset (%.2f)", onset_thr),
                 hjust = 0, size = 2.8, color = "#b2182b") +
        x_scale +
        scale_y_continuous(name = bot_label) +
        labs(
          x       = NULL,
          caption = paste0(
            "Red shading = drought periods.  Ribbon = monthly values.\n",
            "Upper panel line = 12-month centred moving average.  ",
            "Lower panel line = LOESS smooth (span = ", bot_span, ", degree = 2).  ",
            "Def: ", thr_lbl
          )
        ) +
        base_theme

      (p_top / p_bot) +
        plot_layout(heights = c(1.1, 1)) +
        plot_annotation(
          title    = title_txt,
          subtitle = subtitle_txt,
          theme    = theme(
            plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey30")
          )
        )
    }

    # Collector: one row per lag for this index.
    # Populated inside the lag loop; summarised into T4/T5 after it.
    lag_cor_rows <- list()

    # ==========================================================================
    #  LAG LOOP
    #  For index of time-scale n, run lag_k = 0, 1, ..., n-1.
    #  lag_k = 0 : atmospheric field is concurrent with the drought month (t).
    #  lag_k = k : atmospheric field is taken k months BEFORE the drought month (t-k).
    #
    #  Outputs land in:
    #    Figures/{INDEX}/lag00/  (lag_k = 0, concurrent — mirrors original output)
    #    Figures/{INDEX}/lag01/  (lag_k = 1, one month prior)
    #    ...
    #    Tables/{INDEX}/lag00/   etc.
    # ==========================================================================
    for (lag_k in lag_seq) {

      lag_label <- sprintf("lag%02d", lag_k)
      lag_str   <- if (lag_k == 0L) "Lag 0 (concurrent)"
                   else sprintf("Lag -%d month%s prior",
                                lag_k, if (lag_k > 1L) "s" else "")
      lag_short <- if (lag_k == 0L) "Lag 0"
                   else sprintf("Lag -%d mo", lag_k)

      cat(sprintf("   [%s | %s] ", lbl, lag_label))

      lag_fig_dir <- ensure_dir(fig_dir, lag_label)
      lag_tbl_dir <- ensure_dir(tbl_dir, lag_label)

      # ── Lagged raster layer positions ────────────────────────────────────────
      # For each drought month at raster layer position p, the atmospheric
      # layer lag_k months earlier is at p - lag_k.
      n_rast      <- nlyr(z500_anom)
      lag_dr_pos  <- drought_pos    - lag_k
      lag_ndr_pos <- nondrought_pos - lag_k
      lag_dr_pos  <- lag_dr_pos[ lag_dr_pos  >= 1L & lag_dr_pos  <= n_rast]
      lag_ndr_pos <- lag_ndr_pos[lag_ndr_pos >= 1L & lag_ndr_pos <= n_rast]

      n_dr  <- length(lag_dr_pos)
      n_ndr <- length(lag_ndr_pos)

      if (n_dr == 0L) {
        cat(sprintf("No valid lagged layers — skipping.\n")); next
      }

      # ── Lagged composite rasters ─────────────────────────────────────────────
      z500_dr   <- app(subset(z500_anom, lag_dr_pos),  mean, na.rm = TRUE)
      z500_ndr  <- app(subset(z500_anom, lag_ndr_pos), mean, na.rm = TRUE)
      z500_diff <- z500_dr - z500_ndr

      slp_dr   <- app(subset(slp_anom,  lag_dr_pos),  mean, na.rm = TRUE)
      slp_ndr  <- app(subset(slp_anom,  lag_ndr_pos), mean, na.rm = TRUE)
      slp_diff <- slp_dr - slp_ndr

      sst_dr   <- app(subset(sst_anom,  lag_dr_pos),  mean, na.rm = TRUE)
      sst_ndr  <- app(subset(sst_anom,  lag_ndr_pos), mean, na.rm = TRUE)
      sst_diff <- sst_dr - sst_ndr

      # diff_panels closure — captures n_dr / n_ndr for this lag
      diff_panels <- function(r_dr, r_ndr, r_diff) {
        bind_rows(
          as.data.frame(r_dr,   xy = TRUE) %>% setNames(c("x","y","value")) %>%
            mutate(panel = sprintf("Drought months (n=%d)", n_dr)),
          as.data.frame(r_ndr,  xy = TRUE) %>% setNames(c("x","y","value")) %>%
            mutate(panel = sprintf("Non-drought months (n=%d)", n_ndr)),
          as.data.frame(r_diff, xy = TRUE) %>% setNames(c("x","y","value")) %>%
            mutate(panel = "Drought minus Non-drought")
        ) %>% mutate(panel = factor(panel, levels = unique(panel)))
      }

      # ── Lagged area-average scalar time series ───────────────────────────────
      # dplyr::lag(x, n) returns x shifted forward by n, i.e. element [i] = x[i-n].
      # This gives the atmospheric value that occurred lag_k months before each t.
      lag_z500_ts <- dplyr::lag(z500_ridge_ts, lag_k)
      lag_slp_ts  <- dplyr::lag(slp_nwbc_ts,   lag_k)
      lag_sst_ts  <- dplyr::lag(sst_nepac_ts,   lag_k)

      # date_index_lag: atmospheric columns replaced with their lagged versions.
      # The drought index (index_value, is_drought) and dates are unchanged —
      # we are still asking "during months classified as drought at time t, what
      # was the atmosphere doing at time t - lag_k?"
      date_index_lag <- date_index %>%
        mutate(z500_ridge = lag_z500_ts,
               slp_nwbc   = lag_slp_ts,
               sst_nepac  = lag_sst_ts)

      # ── FIG 01: Z500 monthly composites ─────────────────────────────────────
      cat("01 ")
      z500_mc  <- build_monthly_composite_lag(z500_anom, date_index, drought_mask_base, lag_k)
      clim_lim <- ceiling(quantile(abs(z500_mc$value), 0.98, na.rm = TRUE) / 10) * 10

      p01 <- ggplot(z500_mc, aes(x = x, y = y, fill = value)) +
        geom_raster() +
        geom_sf(data = map_layers$borders, fill = NA, color = "black",
                linewidth = 0.3, inherit.aes = FALSE) +
        geom_sf(data = map_layers$lakes,   fill = "white", color = "grey60",
                linewidth = 0.2, inherit.aes = FALSE) +
        scale_fill_distiller(palette = ANOMALY_PALETTE,
                             limits = c(-clim_lim, clim_lim), oob = squish,
                             name = "Z500 Anom\n(m)", na.value = "grey80") +
        coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
        facet_wrap(~panel_title, ncol = 4) +
        labs(
          title    = sprintf("Z500 Anomaly — %s Drought Months by Calendar Month [%s]", lbl, lag_str),
          subtitle = sprintf(
            "Atmospheric field at %s. Drought classification: %s | %s. Grey = no drought (n=0). Ref: %d-%d.",
            lag_str, lbl, thr_lbl, CLIM_START, CLIM_END),
          caption  = "Red = ridge/blocking; Blue = trough. Grey fill = no drought in that calendar month.",
          x = NULL, y = NULL
        ) +
        map_theme(10)

      ggsave(file.path(lag_fig_dir, "01_z500_monthly_anomaly_composites.png"),
             p01, width = FIGURE_WIDTH_WIDE, height = 12, dpi = FIGURE_DPI)

      # ── FIG 02: Z500 drought vs. non-drought ─────────────────────────────────
      cat("02 ")
      diff_lim <- ceiling(quantile(abs(values(z500_diff)), 0.98, na.rm = TRUE) / 5) * 5

      p02 <- ggplot(diff_panels(z500_dr, z500_ndr, z500_diff),
                    aes(x = x, y = y, fill = value)) +
        geom_raster() +
        geom_sf(data = map_layers$borders, fill = NA, color = "black",
                linewidth = 0.3, inherit.aes = FALSE) +
        scale_fill_distiller(palette = ANOMALY_PALETTE,
                             limits = c(-diff_lim, diff_lim), oob = squish,
                             name = "Z500 Anom\n(m)", na.value = "grey80") +
        coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
        facet_wrap(~panel, ncol = 3) +
        labs(title    = sprintf("Z500: %s Drought vs. Non-drought [%s]", lbl, lag_str),
             subtitle = sprintf("%s | %s | %s  |  Ref: %d-%d",
                                lbl, lag_str, thr_lbl, CLIM_START, CLIM_END),
             x = NULL, y = NULL) +
        map_theme(11)

      ggsave(file.path(lag_fig_dir, "02_z500_drought_vs_nondrought.png"),
             p02, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)

      # ── FIG 04: Z500 time series ──────────────────────────────────────────────
      cat("04 ")
      p04 <- plot_dual_panel_ts(
        dates        = date_index_lag$date,
        top_ts       = date_index_lag$z500_ridge,
        bot_ts       = date_index_lag$index_value,
        drought_mask = drought_mask_base,
        top_label    = sprintf("Z500 Anom (m)\n%d-%d W, %d-%d N [%s]",
                               abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                               RIDGE_LAT_MIN, RIDGE_LAT_MAX, lag_short),
        bot_label    = lbl,
        title_txt    = sprintf("Z500 Ridge Anomaly vs. %s — Nechako [%s]", lbl, lag_str),
        subtitle_txt = sprintf("Ridge box: %d-%d W, %d-%d N  |  %s  |  Ref: %d-%d",
                               abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                               RIDGE_LAT_MIN, RIDGE_LAT_MAX,
                               lag_str, CLIM_START, CLIM_END),
        def          = def,
        thr_lbl      = thr_lbl,
        index_label  = lbl
      )
      ggsave(file.path(lag_fig_dir, "04_z500_timeseries_ridge.png"),
             p04, width = FIGURE_WIDTH_WIDE, height = 7, dpi = FIGURE_DPI)

      # ── FIG 05: SLP monthly composites ───────────────────────────────────────
      cat("05 ")
      slp_mc       <- build_monthly_composite_lag(slp_anom, date_index, drought_mask_base, lag_k)
      slp_clim_lim <- ceiling(quantile(abs(slp_mc$value), 0.98, na.rm = TRUE) / 0.5) * 0.5

      p05 <- ggplot(slp_mc, aes(x = x, y = y, fill = value)) +
        geom_raster() +
        geom_sf(data = map_layers$borders, fill = NA, color = "black",
                linewidth = 0.3, inherit.aes = FALSE) +
        geom_sf(data = map_layers$lakes,   fill = "white", color = "grey60",
                linewidth = 0.2, inherit.aes = FALSE) +
        scale_fill_distiller(palette = ANOMALY_PALETTE,
                             limits = c(-slp_clim_lim, slp_clim_lim), oob = squish,
                             name = "SLP Anom\n(hPa)", na.value = "grey80") +
        coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
        facet_wrap(~panel_title, ncol = 4) +
        labs(
          title    = sprintf("SLP Anomaly — %s Drought Months by Calendar Month [%s]", lbl, lag_str),
          subtitle = sprintf(
            "Atmospheric field at %s. Drought classification: %s | %s. Grey = no drought (n=0). Ref: %d-%d.",
            lag_str, lbl, thr_lbl, CLIM_START, CLIM_END),
          caption  = "Red = anomalous high pressure (blocking); Blue = anomalous low pressure.",
          x = NULL, y = NULL
        ) +
        map_theme(10)

      ggsave(file.path(lag_fig_dir, "05_slp_monthly_anomaly_composites.png"),
             p05, width = FIGURE_WIDTH_WIDE, height = 12, dpi = FIGURE_DPI)

      # ── FIG 06: SLP drought vs. non-drought ──────────────────────────────────
      cat("06 ")
      slp_lim <- ceiling(quantile(abs(values(slp_diff)), 0.98, na.rm = TRUE) / 0.5) * 0.5

      p06 <- ggplot(diff_panels(slp_dr, slp_ndr, slp_diff),
                    aes(x = x, y = y, fill = value)) +
        geom_raster() +
        geom_sf(data = map_layers$borders, fill = NA, color = "black",
                linewidth = 0.3, inherit.aes = FALSE) +
        scale_fill_distiller(palette = ANOMALY_PALETTE,
                             limits = c(-slp_lim, slp_lim), oob = squish,
                             name = "SLP Anom\n(hPa)", na.value = "grey80") +
        coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
        facet_wrap(~panel, ncol = 3) +
        labs(title    = sprintf("SLP: %s Drought vs. Non-drought [%s]", lbl, lag_str),
             subtitle = sprintf("%s | %s | %s  |  Ref: %d-%d",
                                lbl, lag_str, thr_lbl, CLIM_START, CLIM_END),
             x = NULL, y = NULL) +
        map_theme(11)

      ggsave(file.path(lag_fig_dir, "06_slp_drought_vs_nondrought.png"),
             p06, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)

      # ── FIG 08: SLP time series ───────────────────────────────────────────────
      cat("08 ")
      p08 <- plot_dual_panel_ts(
        dates        = date_index_lag$date,
        top_ts       = date_index_lag$slp_nwbc,
        bot_ts       = date_index_lag$index_value,
        drought_mask = drought_mask_base,
        top_label    = sprintf("SLP Anom (hPa)\n%d-%d W, %d-%d N [%s]",
                               abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                               RIDGE_LAT_MIN, RIDGE_LAT_MAX, lag_short),
        bot_label    = lbl,
        title_txt    = sprintf("SLP Anomaly vs. %s — Nechako [%s]", lbl, lag_str),
        subtitle_txt = sprintf("NW BC box: %d-%d W, %d-%d N  |  %s  |  Ref: %d-%d",
                               abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                               RIDGE_LAT_MIN, RIDGE_LAT_MAX,
                               lag_str, CLIM_START, CLIM_END),
        def          = def,
        thr_lbl      = thr_lbl
      )
      ggsave(file.path(lag_fig_dir, "08_slp_timeseries.png"),
             p08, width = FIGURE_WIDTH_WIDE, height = 7, dpi = FIGURE_DPI)

      # ── FIG 10: SST monthly composites ───────────────────────────────────────
      cat("10 ")
      sst_mc   <- build_monthly_composite_lag(sst_anom, date_index, drought_mask_base, lag_k)
      sst_mlim <- ceiling(quantile(abs(sst_mc$value), 0.98, na.rm = TRUE) * 2) / 2

      p10 <- ggplot(sst_mc, aes(x = x, y = y, fill = value)) +
        geom_raster() +
        geom_sf(data = map_layers$borders, fill = "grey85", color = "black",
                linewidth = 0.3, inherit.aes = FALSE) +
        scale_fill_distiller(palette = SST_PALETTE,
                             limits = c(-sst_mlim, sst_mlim), oob = squish,
                             name = "SST Anom\n(deg C)", na.value = "grey80") +
        coord_sf(xlim = c(-180, -110), ylim = c(20, 65), expand = FALSE) +
        facet_wrap(~panel_title, ncol = 4) +
        labs(
          title    = sprintf("NE Pacific SST Anomaly — %s Drought Months by Calendar Month [%s]",
                             lbl, lag_str),
          subtitle = sprintf(
            "Ocean state at %s. Drought classification: %s | %s. Grey = no drought (n=0). Ref: %d-%d.",
            lag_str, lbl, thr_lbl, CLIM_START, CLIM_END),
          x = NULL, y = NULL
        ) +
        map_theme(10)

      ggsave(file.path(lag_fig_dir, "10_sst_monthly_composites_drought.png"),
             p10, width = FIGURE_WIDTH_WIDE, height = 12, dpi = FIGURE_DPI)

      # ── FIG 12: Joint scatter plots ───────────────────────────────────────────
      cat("12 ")
      sc_df <- date_index_lag %>%
        filter(!is.na(index_value), !is.na(z500_ridge),
               !is.na(slp_nwbc), !is.na(sst_nepac)) %>%
        mutate(season = factor(season, levels = c("DJF","MAM","JJA","SON")))

      mk_scatter <- function(xvar, xlabel, title_txt) {
        r <- round(cor(sc_df[[xvar]], sc_df$index_value, use = "complete.obs"), 3)
        ggplot(sc_df, aes(x = .data[[xvar]], y = index_value, color = season)) +
          geom_hline(yintercept = if (def$classify == "single") def$threshold else def$onset,
                     linetype = "dashed", color = "darkred", linewidth = 0.5) +
          geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
          geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
          geom_point(alpha = 0.5, size = 1.2) +
          geom_smooth(method = "lm", se = TRUE, color = "black",
                      linewidth = 0.8, fill = "grey70", alpha = 0.3) +
          scale_color_brewer(palette = "Set1", name = "Season") +
          labs(title = title_txt,
               x = xlabel,
               y = sprintf("%s (t)", lbl),
               subtitle = sprintf("r = %.3f  (n = %d)", r, nrow(sc_df))) +
          theme_bw(base_size = 10) +
          theme(plot.title = element_text(face = "bold", size = 10),
                legend.position = "bottom")
      }

      p12 <- (mk_scatter("z500_ridge",
                         sprintf("Z500 Anom (m) at t%s\n%d-%d W, %d-%d N",
                                 if (lag_k == 0) "" else sprintf("-%d", lag_k),
                                 abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                                 RIDGE_LAT_MIN, RIDGE_LAT_MAX),
                         "Z500 Ridge vs. Index") |
               mk_scatter("slp_nwbc",
                          sprintf("SLP Anom (hPa) at t%s\n%d-%d W, %d-%d N",
                                  if (lag_k == 0) "" else sprintf("-%d", lag_k),
                                  abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                                  RIDGE_LAT_MIN, RIDGE_LAT_MAX),
                          "SLP vs. Index") |
               mk_scatter("sst_nepac",
                          sprintf("SST Anom (deg C) at t%s\n%d-%d W, %d-%d N",
                                  if (lag_k == 0) "" else sprintf("-%d", lag_k),
                                  abs(SST_MEAN_LON_MIN), abs(SST_MEAN_LON_MAX),
                                  SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX),
                          "NE Pacific SST vs. Index")) +
        plot_annotation(
          title   = sprintf("Atmospheric/Ocean Anomalies vs. %s (Nechako) [%s]", lbl, lag_str),
          caption = sprintf("Dashed = drought onset. Def: %s. Shading = 95%% CI.", thr_lbl),
          theme   = theme(plot.title = element_text(face = "bold", size = 13))
        )

      ggsave(file.path(lag_fig_dir, "12_joint_scatter.png"),
             p12, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)

      # ── FIG 13: Summary composite panel ──────────────────────────────────────
      cat("13 ")
      lim_z  <- ceiling(quantile(abs(values(z500_diff)), 0.98, na.rm = TRUE) / 10) * 10
      lim_sl <- ceiling(quantile(abs(values(slp_diff)),  0.98, na.rm = TRUE) / 0.5) * 0.5
      lim_ss <- ceiling(quantile(abs(values(sst_diff)),  0.98, na.rm = TRUE) * 2)  / 2

      mk_diff_map <- function(r, lim, title_txt, unit, xlim, ylim) {
        df <- as.data.frame(r, xy = TRUE); names(df)[3] <- "value"
        ggplot(df, aes(x = x, y = y, fill = value)) +
          geom_raster() +
          geom_sf(data = map_layers$borders, fill = NA, color = "black",
                  linewidth = 0.25, inherit.aes = FALSE) +
          scale_fill_distiller(palette = "RdBu",
                               limits = c(-lim, lim), oob = squish,
                               name = unit, na.value = "grey80") +
          coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
          labs(title = title_txt, x = NULL, y = NULL) +
          map_theme(10)
      }

      p13 <- (
        mk_diff_map(z500_diff, lim_z,
                    sprintf("Z500 (m): %s Drought minus Non-drought [%s]", lbl, lag_str),
                    "Z500 (m)", c(-180, -100), c(40, 70)) /
          mk_diff_map(slp_diff,  lim_sl,
                      sprintf("SLP (hPa): %s Drought minus Non-drought [%s]", lbl, lag_str),
                      "SLP (hPa)", c(-180, -100), c(40, 70)) /
          mk_diff_map(sst_diff,  lim_ss,
                      sprintf("SST (deg C): %s Drought minus Non-drought [%s]", lbl, lag_str),
                      "SST (deg C)", c(-180, -110), c(20, 65))
      ) +
        plot_annotation(
          title    = sprintf("Composite Anomalies: %s Drought — Nechako [%s]", lbl, lag_str),
          subtitle = sprintf("%s | %s  (n=%d drought months vs. n=%d non-drought)  |  Ref: %d-%d",
                             thr_lbl, lag_str, n_dr, n_ndr, CLIM_START, CLIM_END),
          caption  = sprintf("ERA5 monthly means. Atmospheric/ocean state at %s. Ref: %d-%d WMO.",
                             lag_str, CLIM_START, CLIM_END),
          theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                           plot.subtitle = element_text(size = 10))
        )

      ggsave(file.path(lag_fig_dir, "13_composite_summary_panel.png"),
             p13, width = FIGURE_WIDTH_STD, height = 16, dpi = FIGURE_DPI)

      # ── TABLES T1-T3 ─────────────────────────────────────────────────────────
      cat("T1-T3 ")

      t1 <- date_index_lag %>%
        filter(!is.na(index_value)) %>%
        mutate(dclass = ifelse(is_drought %in% TRUE,
                               sprintf("Drought [%s]", thr_lbl),
                               "Non-drought")) %>%
        group_by(dclass) %>%
        summarise(n          = n(),
                  index_mean = round(mean(index_value, na.rm = TRUE), 3),
                  index_sd   = round(sd(index_value,   na.rm = TRUE), 3),
                  z500_mean  = round(mean(z500_ridge,   na.rm = TRUE), 1),
                  z500_sd    = round(sd(z500_ridge,     na.rm = TRUE), 1),
                  slp_mean   = round(mean(slp_nwbc,     na.rm = TRUE), 2),
                  slp_sd     = round(sd(slp_nwbc,       na.rm = TRUE), 2),
                  sst_mean   = round(mean(sst_nepac,    na.rm = TRUE), 3),
                  sst_sd     = round(sd(sst_nepac,      na.rm = TRUE), 3),
                  .groups = "drop")
      write_csv(t1, file.path(lag_tbl_dir, "T1_drought_month_summary.csv"))

      t2 <- date_index_lag %>%
        filter(!is.na(index_value), is_drought %in% TRUE) %>%
        group_by(season) %>%
        summarise(n          = n(),
                  index_mean = round(mean(index_value, na.rm = TRUE), 3),
                  z500_mean  = round(mean(z500_ridge,  na.rm = TRUE), 1),
                  z500_sd    = round(sd(z500_ridge,    na.rm = TRUE), 1),
                  slp_mean   = round(mean(slp_nwbc,    na.rm = TRUE), 2),
                  sst_mean   = round(mean(sst_nepac,   na.rm = TRUE), 3),
                  .groups = "drop") %>%
        arrange(match(season, c("DJF","MAM","JJA","SON")))
      write_csv(t2, file.path(lag_tbl_dir, "T2_seasonal_anomaly_statistics.csv"))

      cor_data <- date_index_lag %>%
        filter(!is.na(index_value), !is.na(z500_ridge),
               !is.na(slp_nwbc), !is.na(sst_nepac)) %>%
        select(index_value, z500_ridge, slp_nwbc, sst_nepac)
      cor_labels <- c(lbl,
                      sprintf("Z500 Ridge (m) [%s]", lag_short),
                      sprintf("SLP NW-BC (hPa) [%s]", lag_short),
                      sprintf("SST NE-Pac (degC) [%s]", lag_short))
      n_obs  <- nrow(cor_data)
      cor_r  <- cor(cor_data, use = "complete.obs")
      t_stat <- cor_r * sqrt((n_obs - 2) / (1 - cor_r^2))
      p_mat  <- 2 * pt(-abs(t_stat), df = n_obs - 2)
      cor_out  <- cbind(Variable = cor_labels, as.data.frame(round(cor_r,  3)))
      names(cor_out)[-1]  <- cor_labels
      pval_out <- cbind(Variable = cor_labels, as.data.frame(round(p_mat, 4)))
      names(pval_out)[-1] <- cor_labels
      write_csv(cor_out,  file.path(lag_tbl_dir, "T3_correlation_matrix.csv"))
      write_csv(pval_out, file.path(lag_tbl_dir, "T3_pvalue_matrix.csv"))

      # ── Harvest index-vs-field r and p from the T3 matrix ───────────────────
      # cor_r is a 4×4 matrix; row/col 1 = index_value.
      # Rows 2-4 = z500_ridge, slp_nwbc, sst_nepac.
      lag_cor_rows[[length(lag_cor_rows) + 1L]] <- data.frame(
        index         = lbl,
        definition    = def$short_label,
        lag_k         = lag_k,
        lag_label     = lag_label,
        n_obs         = n_obs,
        r_z500        = round(cor_r["index_value", "z500_ridge"], 4),
        r_slp         = round(cor_r["index_value", "slp_nwbc"],   4),
        r_sst         = round(cor_r["index_value", "sst_nepac"],  4),
        p_z500        = round(p_mat["index_value", "z500_ridge"], 5),
        p_slp         = round(p_mat["index_value", "slp_nwbc"],   5),
        p_sst         = round(p_mat["index_value", "sst_nepac"],  5),
        abs_r_z500    = abs(round(cor_r["index_value", "z500_ridge"], 4)),
        abs_r_slp     = abs(round(cor_r["index_value", "slp_nwbc"],   4)),
        abs_r_sst     = abs(round(cor_r["index_value", "sst_nepac"],  4)),
        stringsAsFactors = FALSE
      )

      cat(sprintf("Done [%s]\n", lag_label))

    }  # end lag loop ─────────────────────────────────────────────────────────

    # ==========================================================================
    #  POST-LAG SUMMARY FOR THIS INDEX
    #  T4: full lag-correlation profile (one row per lag)
    #  T5: best lag per atmospheric variable (one row per variable)
    #  T6: definition-level summary (all indices, best lag per variable)
    # ==========================================================================
    t4 <- bind_rows(lag_cor_rows) %>%
      select(index, definition, lag_k, lag_label, n_obs,
             r_z500, p_z500, r_slp, p_slp, r_sst, p_sst)

    write_csv(t4, file.path(tbl_dir, "T4_lag_correlation_profile.csv"))

    # ── T5: best lag per variable (maximum |r|) ──────────────────────────────
    find_best <- function(var_r, var_p, var_name) {
      rows      <- bind_rows(lag_cor_rows)
      best_row  <- rows[which.max(rows[[var_r]]), ]   # abs_r_* already absolute
      data.frame(
        index         = lbl,
        definition    = def$short_label,
        variable      = var_name,
        best_lag_k    = best_row$lag_k,
        best_lag_label= best_row$lag_label,
        r             = best_row[[sub("abs_r_", "r_", var_r)]],
        abs_r         = best_row[[var_r]],
        p_value       = best_row[[var_p]],
        n_obs         = best_row$n_obs,
        stringsAsFactors = FALSE
      )
    }

    t5 <- bind_rows(
      find_best("abs_r_z500", "p_z500", "Z500_Ridge"),
      find_best("abs_r_slp",  "p_slp",  "SLP_NWBC"),
      find_best("abs_r_sst",  "p_sst",  "SST_NEPac")
    )
    write_csv(t5, file.path(tbl_dir, "T5_best_lag_per_variable.csv"))

    cat(sprintf("   Best lag [%s | %s]:\n", lbl, def$short_label))
    for (i in seq_len(nrow(t5))) {
      cat(sprintf("     %-12s  best lag = %s  r = %+.4f  p = %.5f\n",
                  t5$variable[i], t5$best_lag_label[i],
                  t5$r[i], t5$p_value[i]))
    }

    # Accumulate into the definition-level collector for T6
    def_best_lag_rows[[length(def_best_lag_rows) + 1L]] <- t5

    cat(sprintf("   Done: %s  [%s]\n\n", lbl, def$short_label))

  }  # end per-index loop

  cat(sprintf("  ✓ Definition [%s] complete.\n\n", def$short_label))

  # ── T6: all-index best-lag summary for this definition ─────────────────────
  t6 <- bind_rows(def_best_lag_rows) %>%
    arrange(variable, index)
  t6_path <- file.path(ensure_dir(OUT_DIR, def$short_label, "Tables"),
                       "T6_best_lag_all_indices.csv")
  write_csv(t6, t6_path)
  cat(sprintf("  Saved T6 best-lag summary: %s\n\n", t6_path))

}  # end drought-definition loop

# ============================================================================
#  DONE
# ============================================================================
cat("============================================================\n")
cat(" COMPLETE.\n")
cat(sprintf(" All outputs: %s\n", OUT_DIR))
cat("------------------------------------------------------------\n")
cat(" Output structure (one sub-tree per drought definition):\n")
cat(sprintf(" %s/SW_q10/Tables/T0_drought_months_all_indices.csv\n", OUT_DIR))
cat(sprintf(" %s/SW_q10/Figures/{INDEX}/  &  .../Tables/{INDEX}/\n", OUT_DIR))
cat(sprintf(" %s/Alt2T/Tables/T0_drought_months_all_indices.csv\n",  OUT_DIR))
cat(sprintf(" %s/Alt2T/Figures/{INDEX}/   &  .../Tables/{INDEX}/\n", OUT_DIR))
cat(" Figures/SST/  -- SST figures (index- and definition-independent)\n")
cat("============================================================\n")