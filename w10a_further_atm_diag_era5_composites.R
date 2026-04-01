# ==============================================================================
#  w10a_era5_composites.R
#  Nechako Watershed Drought — All Atmospheric Composite Analysis
#
#  REFACTORED ROLE:  Owns every composite figure and derived table.
#
#  This script assumes that w9_atmospheric_diagnostics.R has already run and
#  produced:
#    (a) The three ERA5 anomaly NetCDF files (Z500, SLP, SST) in DATA_DIR.
#    (b) T0 drought-month classification CSVs in OUT_DIR/{def}/Tables/.
#
#  It loads those files directly — there is no anomaly computation fallback.
#  If the files are absent the script stops with a clear error message
#  directing the user to run w9 first.
#
#  ANALYSIS SECTIONS (all moved from original w9 / original w10a):
#    Section 4    SST time-series figures (Fig 09, Fig 11)
#    Section 4b   2022-2025 drought investigation (Figs 14-21)
#    Section 5    Per-index composite loop (Figs 01-13, Tables T1-T6)
#    Module B1    Historical PDO-stratified composites + severity ranking
#    Module B3    Bootstrap fingerprint significance (parallelised)
#
#  OUTPUTS:
#    {OUT_DIR}/Figures/SST/           09, 11
#    {OUT_DIR}/Figures/Drought_2022_2025/  14-21
#    {OUT_DIR}/{def}/Figures/{INDEX}/{lagXX}/  01-13
#    {OUT_DIR}/{def}/Tables/{INDEX}/{lagXX}/   T1-T3
#    {OUT_DIR}/{def}/Tables/{INDEX}/           T4, T5
#    {OUT_DIR}/{def}/Tables/                   T6
#    composite_results/               B1 composite rasters + correlations
#    bootstrap_results/               B3 p-value rasters + significance maps
#
#  Run AFTER:  w9_atmospheric_diagnostics.R. Can run once w9 has produced the 
#  three anomaly NetCDFs and the T0 classification CSV
#  Run BEFORE: w10b_index_skill.R (independent — can run in parallel)
# ==============================================================================

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
  library(data.table)
  library(parallel)
  library(gridExtra)
})

# ==============================================================================
#  CONFIGURATION  —  must match w9_atmospheric_diagnostics.R
# ==============================================================================
WD_PATH  <- "D:/Nechako_Drought/Nechako/"
DATA_DIR <- "D:/Nechako_Drought/Nechako/monthly_data_direct"
OUT_DIR  <- file.path(WD_PATH, "atmospheric_diagnostics")

setwd(WD_PATH)

SPI_SEAS_DIR  <- file.path(WD_PATH, "spi_results_seasonal/")
SPEI_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal/")

SPI_SCALES_TRY  <- c(1, 3, 6, 9, 12)
SPEI_SCALES_TRY <- c(1, 3, 6, 9, 12)

START_YEAR <- 1950
END_YEAR   <- 2025
CLIM_START <- 1991
CLIM_END   <- 2020

DROUGHT_DEFS <- list(
  
  SW = list(
    id          = "SW",
    short_label = "SW_q10",
    long_label  = "Sheffield & Wood (2008) \u2014 Single threshold (q\u2080=10%, \u2248-1.28)",
    classify    = "single",
    threshold   = qnorm(0.10)
  ),
  
  Alt2T = list(
    id          = "Alt2T",
    short_label = "Alt2T",
    long_label  = "2-Threshold (Onset \u2264-1.0 / Termination \u2265-0.5, min duration by scale)",
    classify    = "hysteresis",
    onset       = -1.0,
    termination = -0.5,
    min_duration_map = list("1" = 2L, "3" = 3L, "6" = 4L, "12" = 6L),
    default_min = 3L
  )
)

RECENT_START        <- 2022L
RECENT_END          <- END_YEAR
DROUGHT_FOCUS_START <- 2022L
DROUGHT_FOCUS_END   <- 2025L

RIDGE_LON_MIN <- -140; RIDGE_LON_MAX <- -110
RIDGE_LAT_MIN <-   50; RIDGE_LAT_MAX <-   65

SST_MEAN_LON_MIN <- -160; SST_MEAN_LON_MAX <- -130
SST_MEAN_LAT_MIN <-   40; SST_MEAN_LAT_MAX <-   55

ANOMALY_PALETTE   <- "RdBu"
SST_PALETTE       <- "RdBu"
FIGURE_DPI        <- 200
FIGURE_WIDTH_WIDE <- 16
FIGURE_WIDTH_STD  <- 12
MONTH_LABELS      <- month.abb

# ==============================================================================
#  SECTION A: LOAD ERA5 ANOMALY NetCDFs  (produced by w9)
#
#  No computation fallback is provided here. If a file is absent the user
#  must run w9_atmospheric_diagnostics.R first.
# ==============================================================================
cat("==============================================================\n")
cat(" SECTION A: Loading ERA5 anomaly NetCDFs (from w9)\n")
cat("==============================================================\n")

.require_nc <- function(path) {
  if (!file.exists(path))
    stop(sprintf(
      "Anomaly file not found:\n  %s\nRun w9_atmospheric_diagnostics.R first.",
      path))
  cat(sprintf("  [load] %s\n", basename(path)))
  terra::rast(path)
}

z500_anom <- .require_nc(file.path(DATA_DIR,
                                   sprintf("z500_monthly_anomaly_clim%d_%d.nc", CLIM_START, CLIM_END)))
slp_anom  <- .require_nc(file.path(DATA_DIR,
                                   sprintf("slp_monthly_anomaly_clim%d_%d.nc",  CLIM_START, CLIM_END)))
sst_anom  <- .require_nc(file.path(DATA_DIR,
                                   sprintf("sst_monthly_anomaly_clim%d_%d.nc",  CLIM_START, CLIM_END)))

# ── Date spine ────────────────────────────────────────────────────────────────
all_dates <- seq.Date(as.Date(sprintf("%d-01-01", START_YEAR)),
                      as.Date(sprintf("%d-12-01", END_YEAR)),
                      by = "month")
n_total <- length(all_dates)

if (terra::nlyr(z500_anom) != n_total)
  warning(sprintf("Z500 layers (%d) != expected (%d).",
                  terra::nlyr(z500_anom), n_total))

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

# ── Area-averaged time series (computed once, reused across every section) ────
ridge_ext <- ext(RIDGE_LON_MIN, RIDGE_LON_MAX, RIDGE_LAT_MIN, RIDGE_LAT_MAX)
sst_ext   <- ext(SST_MEAN_LON_MIN, SST_MEAN_LON_MAX, SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX)

z500_ridge_ts <- terra::global(crop(z500_anom, ridge_ext), "mean", na.rm = TRUE)[, 1]
slp_nwbc_ts   <- terra::global(crop(slp_anom,  ridge_ext), "mean", na.rm = TRUE)[, 1]
sst_nepac_ts  <- terra::global(crop(sst_anom,  sst_ext),   "mean", na.rm = TRUE)[, 1]
cat("  Area-averaged time series extracted.\n\n")

# ==============================================================================
#  SECTION B: LOAD SPI / SPEI INDICES + T0 TABLES  (from w9 outputs)
# ==============================================================================
cat("==============================================================\n")
cat(" SECTION B: Loading SPI / SPEI indices and T0 tables\n")
cat("==============================================================\n")

#' Read a basin-averaged drought index CSV, return data.frame(date, value).
load_index_csv <- function(index_type, scale) {
  dir_map <- list(spi = SPI_SEAS_DIR, spei = SPEI_SEAS_DIR)
  dir     <- dir_map[[tolower(index_type)]]
  f       <- file.path(dir,
                       sprintf("%s_%02d_basin_averaged_by_month.csv",
                               tolower(index_type), as.integer(scale)))
  if (!file.exists(f)) return(NULL)
  
  df <- data.table::fread(f, data.table = FALSE)
  
  if (ncol(df) == 2 && all(c("date", "value") %in% names(df))) {
    long      <- df
    long$date <- as.Date(long$date)
  } else if (ncol(df) >= 13) {
    names(df)[1] <- "Year"
    df$Year      <- as.integer(df$Year)
    month_cols   <- names(df)[2:13]
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

all_indices <- list()
for (sc in SPI_SCALES_TRY) {
  lbl <- sprintf("SPI%d", sc);  df <- load_index_csv("spi", sc)
  if (!is.null(df)) { all_indices[[lbl]] <- df; cat(sprintf("  + %-8s\n", lbl)) }
}
for (sc in SPEI_SCALES_TRY) {
  lbl <- sprintf("SPEI%d", sc); df <- load_index_csv("spei", sc)
  if (!is.null(df)) { all_indices[[lbl]] <- df; cat(sprintf("  + %-8s\n", lbl)) }
}
if (length(all_indices) == 0)
  stop("No SPI/SPEI CSV files found.")

# ── Load T0 for B1 module (SW definition used for consensus) ─────────────────
t0_sw_path <- file.path(OUT_DIR, "SW_q10", "Tables",
                        "T0_drought_months_all_indices.csv")
if (!file.exists(t0_sw_path))
  stop("T0 table not found:\n  ", t0_sw_path,
       "\nRun w9_atmospheric_diagnostics.R first.")
t0_sw <- read_csv(t0_sw_path, show_col_types = FALSE)
cat(sprintf("\n  T0 (SW) loaded: %d rows across %d indices\n\n",
            nrow(t0_sw), n_distinct(t0_sw$index)))

# ==============================================================================
#  HELPER FUNCTIONS
#  Moved here from the original w9 so that w10a is self-contained.
# ==============================================================================

ensure_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

get_map_layers <- function() {
  borders <- ne_countries(scale = "medium", returnclass = "sf")
  LAKES_SHP <- "D:/Nechako_Drought/Nechako/Spatial/ne_50m_lakes/ne_50m_lakes.shp"
  lakes <- tryCatch(sf::st_read(LAKES_SHP, quiet = TRUE), error = function(e) {
    warning("[get_map_layers] Lakes shapefile not found. Continuing without lake polygons.")
    sf::st_sf(geometry = sf::st_sfc(crs = 4326))
  })
  list(borders = borders, lakes = lakes)
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

#' Build a 12-panel composite: mean of qualifying raster layers per calendar month.
#' Months with no qualifying layers produce an all-NA panel (renders grey, n=0).
build_monthly_composite <- function(rast_obj, date_idx, qualifying_mask) {
  n_layers <- terra::nlyr(rast_obj)
  out_list <- vector("list", 12)
  counts   <- integer(12)
  for (m in 1:12) {
    idx <- which(date_idx$month == m & qualifying_mask)
    n   <- length(idx); counts[m] <- n
    mean_layer <- if (n > 0) {
      terra::app(terra::subset(rast_obj, idx), mean, na.rm = TRUE)
    } else {
      rast_obj[[1L]] * NA
    }
    df <- as.data.frame(mean_layer, xy = TRUE); names(df)[3] <- "value"
    if (nrow(df) > 0) {
      df$month_label <- MONTH_LABELS[m]
      df$panel_title <- sprintf("%s (n=%d)", MONTH_LABELS[m], n)
    } else {
      df$month_label <- character(0); df$panel_title <- character(0)
    }
    out_list[[m]] <- df
  }
  bind_rows(out_list) %>%
    mutate(month_label = factor(month_label, levels = MONTH_LABELS),
           panel_title = factor(panel_title,
                                levels = sprintf("%s (n=%d)", MONTH_LABELS, counts)))
}

#' Like build_monthly_composite but shifts the raster layers back by lag_k months.
build_monthly_composite_lag <- function(rast_obj, date_idx, drought_mask_vec, lag_k = 0) {
  n_layers <- terra::nlyr(rast_obj)
  out_list <- vector("list", 12)
  counts   <- integer(12)
  for (m in 1:12) {
    dr_pos <- which(date_idx$month == m & drought_mask_vec)
    idx    <- dr_pos - lag_k
    idx    <- idx[idx >= 1L & idx <= n_layers]
    n      <- length(idx); counts[m] <- n
    mean_layer <- if (n > 0) {
      terra::app(terra::subset(rast_obj, idx), mean, na.rm = TRUE)
    } else {
      rast_obj[[1L]] * NA
    }
    df <- as.data.frame(mean_layer, xy = TRUE); names(df)[3] <- "value"
    if (nrow(df) > 0) {
      df$month_label <- MONTH_LABELS[m]
      df$panel_title <- sprintf("%s (n=%d)", MONTH_LABELS[m], n)
    } else {
      df$month_label <- character(0); df$panel_title <- character(0)
    }
    out_list[[m]] <- df
  }
  bind_rows(out_list) %>%
    mutate(month_label = factor(month_label, levels = MONTH_LABELS),
           panel_title = factor(panel_title,
                                levels = sprintf("%s (n=%d)", MONTH_LABELS, counts)))
}

#' Average all focus-period months within each season.
build_seasonal_mean_df <- function(rast_obj, date_idx) {
  season_lev <- c("DJF", "MAM", "JJA", "SON")
  bind_rows(lapply(season_lev, function(s) {
    idx <- which(date_idx$season == s)
    if (!length(idx)) return(NULL)
    df  <- as.data.frame(
      terra::app(terra::subset(rast_obj, idx), mean, na.rm = TRUE), xy = TRUE)
    names(df)[3] <- "value"
    df$season <- factor(s, levels = season_lev); df
  }))
}

# ── Drought classifiers (identical to w9 — needed for Section 5) ─────────────
extract_scale <- function(lbl) as.integer(regmatches(lbl, regexpr("[0-9]+$", lbl)))

classify_drought_sw <- function(values, threshold) !is.na(values) & values <= threshold

classify_drought_alt <- function(values, onset_thr, term_thr, min_dur) {
  n <- length(values); is_dr <- logical(n); in_d <- FALSE; s_idx <- NA_integer_
  for (j in seq_len(n)) {
    v <- values[j]; if (is.na(v)) next
    if (!in_d && v < onset_thr)      { in_d <- TRUE; s_idx <- j }
    else if (in_d && v >= term_thr) {
      e_idx <- j - 1L
      if ((e_idx - s_idx + 1L) >= min_dur) is_dr[s_idx:e_idx] <- TRUE
      in_d <- FALSE; s_idx <- NA_integer_
    }
  }
  if (in_d && !is.na(s_idx) && (n - s_idx + 1L) >= min_dur) is_dr[s_idx:n] <- TRUE
  is_dr
}

classify_drought <- function(values, def, index_label) {
  if (def$classify == "single") {
    classify_drought_sw(values, def$threshold)
  } else {
    sc_ch <- as.character(extract_scale(index_label))
    mndur <- if (!is.null(def$min_duration_map[[sc_ch]])) def$min_duration_map[[sc_ch]] else def$default_min
    classify_drought_alt(values, def$onset, def$termination, mndur)
  }
}

def_thr_label <- function(def, lbl = NULL) {
  if (def$classify == "single") {
    sprintf("%s \u2264 %.3f", if (!is.null(lbl)) lbl else "Index", def$threshold)
  } else {
    sc_ch <- as.character(if (!is.null(lbl)) extract_scale(lbl) else NA)
    mndur <- if (!is.na(sc_ch) && !is.null(def$min_duration_map[[sc_ch]])) def$min_duration_map[[sc_ch]] else def$default_min
    sprintf("Onset < %.1f / Term \u2265 %.1f / min %d months", def$onset, def$termination, mndur)
  }
}

# ── Load map layers once ──────────────────────────────────────────────────────
cat("  Loading map layers...\n")
map_layers <- suppressMessages(suppressWarnings(get_map_layers()))

# ==============================================================================
#  SECTION 4: SST FIGURES (index-independent, generated once)
# ==============================================================================
cat("==============================================================\n")
cat(" SECTION 4: SST figures\n")
cat("==============================================================\n")
ensure_dir(OUT_DIR, "Figures", "SST")

# ── Fig 09: Annual SST anomaly maps (recent years) ────────────────────────────
cat("  Fig 09: Annual SST anomaly maps (recent years)... ")
recent_years   <- RECENT_START:min(RECENT_END, END_YEAR)
sst_annual_df  <- bind_rows(lapply(recent_years, function(yr) {
  idx <- which(base_date_index$year == yr); if (!length(idx)) return(NULL)
  df  <- as.data.frame(terra::app(terra::subset(sst_anom, idx), mean, na.rm = TRUE), xy = TRUE)
  names(df)[3] <- "value"; df$year <- yr; df
})) %>% mutate(year = factor(year))
sst_lim <- ceiling(max(abs(quantile(sst_annual_df$value, c(0.02, 0.98), na.rm = TRUE))) * 2) / 2

p09 <- ggplot(sst_annual_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = "grey85", color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = SST_PALETTE, limits = c(-sst_lim, sst_lim),
                       oob = squish, name = "SST Anom\n(deg C)", na.value = "grey80") +
  coord_sf(xlim = c(-180, -110), ylim = c(20, 65), expand = FALSE) +
  facet_wrap(~year, ncol = 4) +
  labs(title    = sprintf("NE Pacific SST Anomaly \u2014 2022\u20132025 Drought Period (Annual Mean)"),
       subtitle = sprintf("Reference climatology: %d\u2013%d (WMO standard)", CLIM_START, CLIM_END),
       x = NULL, y = NULL) +
  map_theme(10)
ggsave(file.path(OUT_DIR, "Figures", "SST", "09_sst_annual_anomaly_recent.png"),
       p09, width = FIGURE_WIDTH_WIDE,
       height = ceiling(length(recent_years) / 4) * 4 + 2, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── Fig 11: Full-record NE Pacific SST time series ────────────────────────────
cat("  Fig 11: NE Pacific SST time series... ")
sst_ts <- base_date_index %>%
  mutate(sst_nepac = sst_nepac_ts) %>%
  filter(!is.na(sst_nepac)) %>%
  arrange(date) %>%
  mutate(sst_12mo = rollmean(sst_nepac, k = 12, fill = NA, align = "right"))

p11 <- ggplot(sst_ts, aes(x = date, y = sst_nepac)) +
  annotate("rect",
           xmin = as.Date(sprintf("%d-01-01", DROUGHT_FOCUS_START)),
           xmax = as.Date(sprintf("%d-12-31", DROUGHT_FOCUS_END)),
           ymin = -Inf, ymax = Inf, fill = "#fee08b", alpha = 0.55) +
  annotate("rect",
           xmin = as.Date(sprintf("%d-01-01", DROUGHT_FOCUS_START)),
           xmax = as.Date(sprintf("%d-12-31", DROUGHT_FOCUS_END)),
           ymin = -Inf, ymax = Inf, fill = NA, color = "#d73027",
           linewidth = 0.6, linetype = "dashed") +
  geom_col(aes(fill = sst_nepac > 0), alpha = 0.4, width = 25) +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"), guide = "none") +
  geom_line(aes(y = sst_12mo), color = "black", linewidth = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 0,    color = "grey30",   linewidth = 0.4) +
  geom_hline(yintercept =  0.5, linetype = "dotted", color = "darkred",   linewidth = 0.5) +
  geom_hline(yintercept = -0.5, linetype = "dotted", color = "steelblue", linewidth = 0.5) +
  annotate("text",
           x     = as.Date(sprintf("%d-07-01", DROUGHT_FOCUS_START)),
           y     = sst_lim * 0.88,
           label = "2022\u20132025\nDrought",
           hjust = 0.5, size = 3, fontface = "bold", color = "#d73027") +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y", expand = c(0.01, 0)) +
  labs(title    = "NE Pacific SST Anomaly \u2014 Monthly Time Series (1950\u20132025)",
       subtitle = sprintf("Area average: %d\u2013%d W, %d\u2013%d N  |  Black line: 12-month rolling mean",
                          abs(SST_MEAN_LON_MIN), abs(SST_MEAN_LON_MAX),
                          SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX),
       x = NULL, y = "SST Anomaly (deg C)") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "Figures", "SST", "11_sst_timeseries_ne_pacific.png"),
       p11, width = FIGURE_WIDTH_WIDE, height = 5.5, dpi = FIGURE_DPI)
cat("Saved.\n\n")

# ==============================================================================
#  SECTION 4b: 2022-2025 DROUGHT PERIOD — YEAR-BY-YEAR & SEASONAL INVESTIGATION
# ==============================================================================
cat("==============================================================\n")
cat(sprintf(" SECTION 4b: 2022\u20132025 Drought Period Investigation\n"))
cat("==============================================================\n")

drought_fig_dir <- ensure_dir(OUT_DIR, "Figures", "Drought_2022_2025")

focus_mask    <- base_date_index$year >= DROUGHT_FOCUS_START & base_date_index$year <= DROUGHT_FOCUS_END
pre_mask      <- base_date_index$year >= START_YEAR         & base_date_index$year <  DROUGHT_FOCUS_START
focus_idx     <- which(focus_mask)
pre_idx       <- which(pre_mask)
focus_years   <- DROUGHT_FOCUS_START:DROUGHT_FOCUS_END
focus_date_idx <- base_date_index %>% filter(year >= DROUGHT_FOCUS_START, year <= DROUGHT_FOCUS_END)

# ── Figs 14-16: Year-by-year annual anomaly maps ──────────────────────────────
make_annual_df <- function(anom_rast) {
  bind_rows(lapply(focus_years, function(yr) {
    idx <- which(base_date_index$year == yr); if (!length(idx)) return(NULL)
    df  <- as.data.frame(terra::app(terra::subset(anom_rast, idx), mean, na.rm = TRUE), xy = TRUE)
    names(df)[3] <- "value"; df$year <- yr; df
  })) %>% mutate(year = factor(year))
}

# Fig 14: Z500
cat("  Fig 14: Z500 annual anomaly (2022-2025)... ")
z500_yr_df  <- make_annual_df(z500_anom)
z500_yr_lim <- ceiling(quantile(abs(z500_yr_df$value), 0.98, na.rm = TRUE) / 10) * 10
p14 <- ggplot(z500_yr_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = NA, color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_sf(data = map_layers$lakes,   fill = "white", color = "grey60", linewidth = 0.2, inherit.aes = FALSE) +
  scale_fill_distiller(palette = ANOMALY_PALETTE, limits = c(-z500_yr_lim, z500_yr_lim),
                       oob = squish, name = "Z500 Anom\n(m)", na.value = "grey80") +
  coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
  facet_wrap(~year, ncol = 4) +
  labs(title = "Z500 Anomaly \u2014 2022\u20132025 Drought Period (Annual Mean per Year)",
       subtitle = sprintf("ERA5 Z500 anomaly.  Ref: %d\u2013%d (WMO).", CLIM_START, CLIM_END),
       x = NULL, y = NULL) +
  map_theme(11)
ggsave(file.path(drought_fig_dir, "14_z500_annual_drought_period.png"),
       p14, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)
cat("Saved.\n")

# Fig 15: SLP
cat("  Fig 15: SLP annual anomaly (2022-2025)... ")
slp_yr_df  <- make_annual_df(slp_anom)
slp_yr_lim <- ceiling(quantile(abs(slp_yr_df$value), 0.98, na.rm = TRUE) / 0.5) * 0.5
p15 <- ggplot(slp_yr_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = NA, color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  geom_sf(data = map_layers$lakes,   fill = "white", color = "grey60", linewidth = 0.2, inherit.aes = FALSE) +
  scale_fill_distiller(palette = ANOMALY_PALETTE, limits = c(-slp_yr_lim, slp_yr_lim),
                       oob = squish, name = "SLP Anom\n(hPa)", na.value = "grey80") +
  coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
  facet_wrap(~year, ncol = 4) +
  labs(title = "SLP Anomaly \u2014 2022\u20132025 Drought Period (Annual Mean per Year)",
       subtitle = sprintf("ERA5 SLP anomaly.  Ref: %d\u2013%d (WMO).", CLIM_START, CLIM_END),
       x = NULL, y = NULL) +
  map_theme(11)
ggsave(file.path(drought_fig_dir, "15_slp_annual_drought_period.png"),
       p15, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)
cat("Saved.\n")

# Fig 16: SST
cat("  Fig 16: SST annual anomaly (2022-2025)... ")
sst_yr_df  <- make_annual_df(sst_anom)
sst_yr_lim <- ceiling(max(abs(quantile(sst_yr_df$value, c(0.02, 0.98), na.rm = TRUE))) * 2) / 2
p16 <- ggplot(sst_yr_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = "grey85", color = "black", linewidth = 0.3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = SST_PALETTE, limits = c(-sst_yr_lim, sst_yr_lim),
                       oob = squish, name = "SST Anom\n(deg C)", na.value = "grey80") +
  coord_sf(xlim = c(-180, -110), ylim = c(20, 65), expand = FALSE) +
  facet_wrap(~year, ncol = 4) +
  labs(title = "NE Pacific SST Anomaly \u2014 2022\u20132025 Drought Period (Annual Mean per Year)",
       subtitle = sprintf("Ref: %d\u2013%d (WMO).", CLIM_START, CLIM_END),
       x = NULL, y = NULL) +
  map_theme(11)
ggsave(file.path(drought_fig_dir, "16_sst_annual_drought_period.png"),
       p16, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── Figs 17-19: Seasonal composites for the drought period ───────────────────
cat("  Figs 17-19: Seasonal composites (2022-2025 mean)... ")
season_lev <- c("DJF", "MAM", "JJA", "SON")

make_seas_panel <- function(anom_rast, lim, palette, unit_lbl, title_txt, xlim, ylim,
                            map_fill = NA) {
  df  <- build_seasonal_mean_df(terra::subset(anom_rast, which(focus_mask)), focus_date_idx)
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = map_layers$borders, fill = map_fill, color = "black",
            linewidth = 0.3, inherit.aes = FALSE) +
    scale_fill_distiller(palette = palette, limits = c(-lim, lim),
                         oob = squish, name = unit_lbl, na.value = "grey80") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    facet_wrap(~season, ncol = 4) +
    labs(title    = title_txt,
         subtitle = sprintf("Mean across all %d\u2013%d months within each season.  Ref: %d\u2013%d.",
                            DROUGHT_FOCUS_START, DROUGHT_FOCUS_END, CLIM_START, CLIM_END),
         x = NULL, y = NULL) +
    map_theme(11)
}

z500_seas_df  <- build_seasonal_mean_df(terra::subset(z500_anom, which(focus_mask)), focus_date_idx)
slp_seas_df   <- build_seasonal_mean_df(terra::subset(slp_anom,  which(focus_mask)), focus_date_idx)
sst_seas_df   <- build_seasonal_mean_df(terra::subset(sst_anom,  which(focus_mask)), focus_date_idx)

z500_seas_lim <- ceiling(quantile(abs(z500_seas_df$value), 0.98, na.rm = TRUE) / 10) * 10
slp_seas_lim  <- ceiling(quantile(abs(slp_seas_df$value),  0.98, na.rm = TRUE) / 0.5) * 0.5
sst_seas_lim  <- ceiling(quantile(abs(sst_seas_df$value),  0.98, na.rm = TRUE) * 2)  / 2

p17 <- make_seas_panel(z500_anom, z500_seas_lim, ANOMALY_PALETTE, "Z500 Anom\n(m)",
                       "Z500 Anomaly \u2014 2022\u20132025 Drought Period by Season", c(-180,-100), c(40,70))
p18 <- make_seas_panel(slp_anom,  slp_seas_lim,  ANOMALY_PALETTE, "SLP Anom\n(hPa)",
                       "SLP Anomaly \u2014 2022\u20132025 Drought Period by Season",  c(-180,-100), c(40,70))
p19 <- make_seas_panel(sst_anom,  sst_seas_lim,  SST_PALETTE,     "SST Anom\n(deg C)",
                       "NE Pacific SST Anomaly \u2014 2022\u20132025 Drought Period by Season", c(-180,-110), c(20,65),
                       map_fill = "grey85")

ggsave(file.path(drought_fig_dir, "17_z500_seasonal_drought_period.png"),
       p17, width = FIGURE_WIDTH_WIDE, height = 5, dpi = FIGURE_DPI)
ggsave(file.path(drought_fig_dir, "18_slp_seasonal_drought_period.png"),
       p18, width = FIGURE_WIDTH_WIDE, height = 5, dpi = FIGURE_DPI)
ggsave(file.path(drought_fig_dir, "19_sst_seasonal_drought_period.png"),
       p19, width = FIGURE_WIDTH_WIDE, height = 5, dpi = FIGURE_DPI)
cat("Saved (17-19).\n")

# ── Fig 20: Drought fingerprint panel ─────────────────────────────────────────
cat("  Fig 20: Drought fingerprint panel (2022-2025 vs 1950-2021)... ")

make_fingerprint_row <- function(rast_obj, focus_idx, pre_idx, lim_round, unit_label,
                                 title_prefix, xlim, ylim, map_fill = "grey85") {
  r_focus <- terra::app(terra::subset(rast_obj, focus_idx), mean, na.rm = TRUE)
  r_pre   <- terra::app(terra::subset(rast_obj, pre_idx),   mean, na.rm = TRUE)
  r_diff  <- r_focus - r_pre
  lim <- ceiling(max(
    quantile(abs(terra::values(r_focus)), 0.98, na.rm = TRUE),
    quantile(abs(terra::values(r_pre)),   0.98, na.rm = TRUE),
    quantile(abs(terra::values(r_diff)),  0.98, na.rm = TRUE)
  ) / lim_round) * lim_round
  
  to_df <- function(r, panel_lbl) {
    df <- as.data.frame(r, xy = TRUE); names(df)[3] <- "value"
    df$panel <- panel_lbl; df
  }
  plot_df <- bind_rows(
    to_df(r_focus, sprintf("2022\u20132025 Drought Period\n(n = %d months)", length(focus_idx))),
    to_df(r_pre,   sprintf("1950\u20132021 Background\n(n = %d months)", length(pre_idx))),
    to_df(r_diff,  "Drought minus Background")
  ) %>% mutate(panel = factor(panel, levels = unique(panel)))
  
  ggplot(plot_df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = map_layers$borders, fill = map_fill, color = "black",
            linewidth = 0.25, inherit.aes = FALSE) +
    scale_fill_distiller(palette = "RdBu", limits = c(-lim, lim),
                         oob = squish, name = unit_label, na.value = "grey80") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    facet_wrap(~panel, ncol = 3) +
    labs(title = title_prefix, x = NULL, y = NULL) +
    map_theme(10)
}

p20 <- (make_fingerprint_row(z500_anom, focus_idx, pre_idx, 10,  "Z500 (m)",    "Z500 Anomaly (m)",    c(-180,-100), c(40,70), NA)   /
          make_fingerprint_row(slp_anom,  focus_idx, pre_idx, 0.5, "SLP (hPa)",   "SLP Anomaly (hPa)",   c(-180,-100), c(40,70), NA)   /
          make_fingerprint_row(sst_anom,  focus_idx, pre_idx, 0.5, "SST (\u00b0C)", "SST Anomaly (\u00b0C)", c(-180,-110), c(20,65))) +
  plot_annotation(
    title    = "2022\u20132025 Drought Fingerprint: Atmospheric and Oceanic Anomalies",
    subtitle = sprintf("Left: 2022\u20132025 | Centre: 1950\u20132021 background | Right: difference\nRef: %d\u2013%d (WMO).",
                       CLIM_START, CLIM_END),
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 9))
  )
ggsave(file.path(drought_fig_dir, "20_drought_fingerprint_panel.png"),
       p20, width = FIGURE_WIDTH_STD, height = 16, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── Fig 21: Zoomed time-series 2019-2025 ──────────────────────────────────────
cat("  Fig 21: Zoomed time series 2019-2025... ")
zoom_start <- as.Date("2019-01-01")
zoom_end   <- as.Date("2025-12-31")
sst_ts_zoom <- sst_ts %>% filter(date >= zoom_start, date <= zoom_end)
z500_zoom   <- tibble(date = base_date_index$date, z500_ridge = z500_ridge_ts,
                      slp_nwbc = slp_nwbc_ts, sst_nepac = sst_nepac_ts) %>%
  filter(date >= zoom_start, date <= zoom_end)

drought_shade  <- annotate("rect", xmin = as.Date(sprintf("%d-01-01", DROUGHT_FOCUS_START)),
                           xmax = as.Date(sprintf("%d-12-31", DROUGHT_FOCUS_END)),
                           ymin = -Inf, ymax = Inf, fill = "#fee08b", alpha = 0.45)
drought_border <- annotate("rect", xmin = as.Date(sprintf("%d-01-01", DROUGHT_FOCUS_START)),
                           xmax = as.Date(sprintf("%d-12-31", DROUGHT_FOCUS_END)),
                           ymin = -Inf, ymax = Inf, fill = NA, color = "#d73027", linewidth = 0.5, linetype = "dashed")
x_sc <- scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.01, 0))
base_thm <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                             plot.title = element_text(face = "bold"),
                                             axis.text.x = element_text(size = 9))

make_zoom_panel <- function(df, var, ylab) {
  ggplot(df, aes(x = date, y = .data[[var]])) +
    drought_shade + drought_border +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_col(aes(fill = .data[[var]] > 0), alpha = 0.5, width = 25) +
    scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"), guide = "none") +
    geom_line(color = "grey15", linewidth = 0.8,
              data = df %>% mutate(sm = rollmean(.data[[var]], 6, fill = NA, align = "center")),
              aes(y = sm)) +
    x_sc + labs(title = NULL, y = ylab, x = NULL) + base_thm
}

p21 <- (make_zoom_panel(z500_zoom, "z500_ridge",
                        sprintf("Z500 Anom (m)\n%d\u2013%d W, %d\u2013%d N",
                                abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX), RIDGE_LAT_MIN, RIDGE_LAT_MAX)) /
          make_zoom_panel(z500_zoom, "slp_nwbc",   "SLP Anom (hPa)") /
          make_zoom_panel(z500_zoom, "sst_nepac",
                          sprintf("SST Anom (\u00b0C)\n%d\u2013%d W, %d\u2013%d N",
                                  abs(SST_MEAN_LON_MIN), abs(SST_MEAN_LON_MAX),
                                  SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX))) +
  plot_annotation(
    title    = "Atmospheric & Oceanic Time Series \u2014 2019\u20132025 Drought Context",
    subtitle = sprintf("Area-average anomalies.  Ref: %d\u2013%d.  Line = 6-month centred mean.",
                       CLIM_START, CLIM_END),
    theme    = theme(plot.title = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9))
  )
ggsave(file.path(drought_fig_dir, "21_atmos_ocean_timeseries_zoom.png"),
       p21, width = FIGURE_WIDTH_STD, height = 10, dpi = FIGURE_DPI)
cat("Saved.\n\n")

# ==============================================================================
#  SECTION 5: PER-INDEX × DROUGHT DEFINITION COMPOSITE LOOP
#  (Figs 01-13, Tables T1-T6; drought classifications re-derived here
#   from the loaded SPI/SPEI data to keep w10a fully self-contained.)
# ==============================================================================
cat("==============================================================\n")
cat(" SECTION 5: Per-index composites and tables (both definitions)\n")
cat("==============================================================\n\n")

for (def in DROUGHT_DEFS) {
  
  cat(sprintf("\n\u2554%s\u2557\n", strrep("\u2550", 54)))
  cat(sprintf("  Definition: %s\n", def$long_label))
  cat(sprintf("\u255a%s\u255d\n\n", strrep("\u2550", 54)))
  
  def_fig_root <- ensure_dir(OUT_DIR, def$short_label, "Figures")
  def_tbl_root <- ensure_dir(OUT_DIR, def$short_label, "Tables")
  def_best_lag_rows <- list()
  
  for (lbl in names(all_indices)) {
    
    cat(sprintf(">> %s  [%s]\n", lbl, def$short_label))
    
    fig_dir <- ensure_dir(def_fig_root, lbl)
    tbl_dir <- ensure_dir(def_tbl_root, lbl)
    
    idx_df <- all_indices[[lbl]] %>%
      filter(as.integer(format(date, "%Y")) >= START_YEAR,
             as.integer(format(date, "%Y")) <= END_YEAR) %>%
      rename(index_value = value) %>%
      mutate(is_drought = classify_drought(index_value, def, lbl))
    
    thr_lbl <- def_thr_label(def, lbl)
    
    date_index <- base_date_index %>%
      left_join(idx_df %>% select(date, index_value, is_drought), by = "date") %>%
      mutate(z500_ridge = z500_ridge_ts, slp_nwbc = slp_nwbc_ts, sst_nepac = sst_nepac_ts)
    
    drought_mask_base    <- date_index$is_drought %in% TRUE
    nondrought_mask_base <- date_index$is_drought %in% FALSE
    drought_pos          <- which(drought_mask_base)
    nondrought_pos       <- which(nondrought_mask_base)
    n_dr_base            <- length(drought_pos)
    
    if (n_dr_base == 0) { cat("   No drought months — skipping.\n\n"); next }
    
    n_scale  <- extract_scale(lbl)
    lag_seq  <- 0L:(n_scale - 1L)
    n_rast   <- terra::nlyr(z500_anom)
    lag_cor_rows <- list()
    
    cat(sprintf("   Drought months [%s]: %d | Non-drought: %d | Lags: 0..%d\n",
                thr_lbl, n_dr_base, length(nondrought_pos), n_scale - 1L))
    
    # ── Inner LOESS/dual-panel time-series plot function ──────────────────────
    loess_span_map <- c("1"=0.10,"3"=0.08,"6"=0.07,"9"=0.06,"12"=0.05)
    
    plot_dual_panel_ts <- function(dates, top_ts, bot_ts, drought_mask,
                                   top_label, bot_label, title_txt, subtitle_txt,
                                   def, thr_lbl, index_label = NULL) {
      sc_key   <- as.character(as.integer(regmatches(index_label, regexpr("[0-9]+$", index_label))))
      bot_span <- if (!is.null(index_label) && sc_key %in% names(loess_span_map)) loess_span_map[[sc_key]] else 0.08
      df_raw   <- tibble(date = dates, top = top_ts, bot = bot_ts, is_drought = drought_mask) %>%
        filter(!is.na(top), !is.na(bot)) %>% arrange(date)
      loess_fit <- loess(bot ~ as.numeric(date), data = df_raw, span = bot_span, degree = 2, na.action = na.exclude)
      df <- df_raw %>% mutate(
        top_smooth = rollmean(top, k = 12, fill = NA, align = "center"),
        bot_smooth = predict(loess_fit, newdata = data.frame(date = as.numeric(date)))
      )
      drought_runs <- df %>%
        mutate(grp = cumsum(c(TRUE, diff(as.integer(is_drought)) != 0))) %>%
        group_by(grp) %>%
        summarise(xmin = min(date), xmax = max(date), drought = first(is_drought), .groups = "drop") %>%
        filter(drought)
      shade       <- geom_rect(data = drought_runs, aes(xmin=xmin,xmax=xmax,ymin=-Inf,ymax=Inf),
                               fill = "#f4a582", alpha = 0.25, inherit.aes = FALSE)
      focus_shade <- annotate("rect", xmin=as.Date(sprintf("%d-01-01",DROUGHT_FOCUS_START)),
                              xmax=as.Date(sprintf("%d-12-31",DROUGHT_FOCUS_END)),
                              ymin=-Inf, ymax=Inf, fill="#fee08b", alpha=0.35)
      focus_bord  <- annotate("rect", xmin=as.Date(sprintf("%d-01-01",DROUGHT_FOCUS_START)),
                              xmax=as.Date(sprintf("%d-12-31",DROUGHT_FOCUS_END)),
                              ymin=-Inf, ymax=Inf, fill=NA, color="#d73027", linewidth=0.45, linetype="dashed")
      x_sc   <- scale_x_date(date_breaks="10 years",date_labels="%Y",expand=c(0.01,0))
      bthm   <- theme_minimal(base_size=11) +
        theme(plot.title=element_text(face="bold",size=12), plot.subtitle=element_text(size=9,color="grey30"),
              panel.grid.minor=element_blank(), axis.title.y=element_text(size=9),
              axis.text.x=element_text(size=9), plot.margin=margin(4,8,2,8))
      top_lim <- ceiling(max(abs(df$top),na.rm=TRUE)/10)*10
      onset_thr <- if (def$classify=="single") def$threshold else def$onset
      p_top <- ggplot(df,aes(x=date))+focus_shade+focus_bord+shade+
        geom_ribbon(aes(ymin=0,ymax=top),fill="#92c5de",alpha=0.55)+
        geom_ribbon(aes(ymin=pmin(top,0),ymax=0),fill="#4393c3",alpha=0.35)+
        geom_ribbon(aes(ymin=0,ymax=pmax(top,0)),fill="#d6604d",alpha=0.55)+
        geom_line(aes(y=top_smooth),color="grey15",linewidth=0.9,na.rm=TRUE)+
        geom_hline(yintercept=0,color="grey40",linewidth=0.4)+
        x_sc+scale_y_continuous(name=top_label,limits=c(-top_lim,top_lim))+labs(x=NULL)+bthm
      p_bot <- ggplot(df,aes(x=date))+focus_shade+focus_bord+shade+
        geom_ribbon(aes(ymin=pmin(bot,0),ymax=0),fill="#d6604d",alpha=0.4)+
        geom_ribbon(aes(ymin=0,ymax=pmax(bot,0)),fill="#92c5de",alpha=0.4)+
        geom_line(aes(y=bot_smooth),color="grey15",linewidth=0.85,na.rm=TRUE)+
        geom_hline(yintercept=onset_thr,linetype="dashed",color="#b2182b",linewidth=0.55)+
        geom_hline(yintercept=0,color="grey40",linewidth=0.4)+
        annotate("text",x=min(df$date),y=onset_thr+0.08,
                 label=sprintf("Onset (%.2f)",onset_thr),hjust=0,size=2.8,color="#b2182b")+
        x_sc+scale_y_continuous(name=bot_label)+labs(x=NULL,caption=paste0(
          "Red shading = all drought periods.  Yellow shading = 2022\u20132025 focus period.  Def: ",thr_lbl))+bthm
      (p_top/p_bot)+plot_layout(heights=c(1.1,1))+
        plot_annotation(title=title_txt,subtitle=subtitle_txt,
                        theme=theme(plot.title=element_text(face="bold",size=13),
                                    plot.subtitle=element_text(size=9,color="grey30")))
    }
    
    # ── Lag loop ──────────────────────────────────────────────────────────────
    for (lag_k in lag_seq) {
      
      lag_label <- sprintf("lag%02d", lag_k)
      lag_str   <- if (lag_k == 0L) "Lag 0 (concurrent)" else sprintf("Lag -%d month%s prior", lag_k, if (lag_k>1L) "s" else "")
      lag_short <- if (lag_k == 0L) "Lag 0" else sprintf("Lag -%d mo", lag_k)
      cat(sprintf("   [%s | %s] ", lbl, lag_label))
      
      lag_fig_dir <- ensure_dir(fig_dir, lag_label)
      lag_tbl_dir <- ensure_dir(tbl_dir, lag_label)
      
      lag_dr_pos  <- (drought_pos    - lag_k) ; lag_dr_pos  <- lag_dr_pos[lag_dr_pos  >= 1L & lag_dr_pos  <= n_rast]
      lag_ndr_pos <- (nondrought_pos - lag_k) ; lag_ndr_pos <- lag_ndr_pos[lag_ndr_pos >= 1L & lag_ndr_pos <= n_rast]
      n_dr  <- length(lag_dr_pos) ; n_ndr <- length(lag_ndr_pos)
      if (n_dr == 0L) { cat("no valid lagged layers — skip.\n"); next }
      
      z500_dr  <- terra::app(terra::subset(z500_anom, lag_dr_pos),  mean, na.rm = TRUE)
      z500_ndr <- terra::app(terra::subset(z500_anom, lag_ndr_pos), mean, na.rm = TRUE)
      z500_diff <- z500_dr - z500_ndr
      slp_dr   <- terra::app(terra::subset(slp_anom,  lag_dr_pos),  mean, na.rm = TRUE)
      slp_ndr  <- terra::app(terra::subset(slp_anom,  lag_ndr_pos), mean, na.rm = TRUE)
      slp_diff  <- slp_dr - slp_ndr
      sst_dr   <- terra::app(terra::subset(sst_anom,  lag_dr_pos),  mean, na.rm = TRUE)
      sst_ndr  <- terra::app(terra::subset(sst_anom,  lag_ndr_pos), mean, na.rm = TRUE)
      sst_diff  <- sst_dr - sst_ndr
      
      diff_panels <- function(r_dr, r_ndr, r_diff) {
        bind_rows(
          as.data.frame(r_dr,   xy=TRUE) %>% setNames(c("x","y","value")) %>% mutate(panel=sprintf("Drought months (n=%d)", n_dr)),
          as.data.frame(r_ndr,  xy=TRUE) %>% setNames(c("x","y","value")) %>% mutate(panel=sprintf("Non-drought months (n=%d)", n_ndr)),
          as.data.frame(r_diff, xy=TRUE) %>% setNames(c("x","y","value")) %>% mutate(panel="Drought minus Non-drought")
        ) %>% mutate(panel = factor(panel, levels = unique(panel)))
      }
      
      lag_z500_ts <- dplyr::lag(z500_ridge_ts, lag_k)
      lag_slp_ts  <- dplyr::lag(slp_nwbc_ts,   lag_k)
      lag_sst_ts  <- dplyr::lag(sst_nepac_ts,   lag_k)
      date_index_lag <- date_index %>%
        mutate(z500_ridge = lag_z500_ts, slp_nwbc = lag_slp_ts, sst_nepac = lag_sst_ts)
      
      # Fig 01 — Z500 monthly composites
      cat("01 ")
      z500_mc  <- build_monthly_composite_lag(z500_anom, date_index, drought_mask_base, lag_k)
      clim_lim <- ceiling(quantile(abs(z500_mc$value), 0.98, na.rm = TRUE) / 10) * 10
      p01 <- ggplot(z500_mc, aes(x=x,y=y,fill=value))+geom_raster()+
        geom_sf(data=map_layers$borders,fill=NA,color="black",linewidth=0.3,inherit.aes=FALSE)+
        geom_sf(data=map_layers$lakes,fill="white",color="grey60",linewidth=0.2,inherit.aes=FALSE)+
        scale_fill_distiller(palette=ANOMALY_PALETTE,limits=c(-clim_lim,clim_lim),oob=squish,name="Z500 Anom\n(m)",na.value="grey80")+
        coord_sf(xlim=c(-180,-100),ylim=c(40,70),expand=FALSE)+facet_wrap(~panel_title,ncol=4)+
        labs(title=sprintf("Z500 Anomaly \u2014 %s Drought Months [%s]",lbl,lag_str),
             subtitle=sprintf("Field at %s. Def: %s. Ref: %d\u2013%d.",lag_str,thr_lbl,CLIM_START,CLIM_END),
             x=NULL,y=NULL)+map_theme(10)
      ggsave(file.path(lag_fig_dir,"01_z500_monthly_anomaly_composites.png"),p01,width=FIGURE_WIDTH_WIDE,height=12,dpi=FIGURE_DPI)
      
      # Fig 02 — Z500 drought vs non-drought
      cat("02 ")
      diff_lim <- ceiling(quantile(abs(terra::values(z500_diff)),0.98,na.rm=TRUE)/5)*5
      p02 <- ggplot(diff_panels(z500_dr,z500_ndr,z500_diff),aes(x=x,y=y,fill=value))+geom_raster()+
        geom_sf(data=map_layers$borders,fill=NA,color="black",linewidth=0.3,inherit.aes=FALSE)+
        scale_fill_distiller(palette=ANOMALY_PALETTE,limits=c(-diff_lim,diff_lim),oob=squish,name="Z500 Anom\n(m)",na.value="grey80")+
        coord_sf(xlim=c(-180,-100),ylim=c(40,70),expand=FALSE)+facet_wrap(~panel,ncol=3)+
        labs(title=sprintf("Z500: %s Drought vs. Non-drought [%s]",lbl,lag_str),
             subtitle=sprintf("%s | %s | %s | Ref: %d\u2013%d",lbl,lag_str,thr_lbl,CLIM_START,CLIM_END),
             x=NULL,y=NULL)+map_theme(11)
      ggsave(file.path(lag_fig_dir,"02_z500_drought_vs_nondrought.png"),p02,width=FIGURE_WIDTH_WIDE,height=6,dpi=FIGURE_DPI)
      
      # Fig 04 — Z500 time series
      cat("04 ")
      p04 <- plot_dual_panel_ts(date_index_lag$date,date_index_lag$z500_ridge,date_index_lag$index_value,
                                drought_mask_base,sprintf("Z500 Anom (m)\n%d\u2013%d W [%s]",abs(RIDGE_LON_MIN),abs(RIDGE_LON_MAX),lag_short),lbl,
                                sprintf("Z500 Ridge Anomaly vs. %s [%s]",lbl,lag_str),
                                sprintf("Ridge box: %d\u2013%d W, %d\u2013%d N | Ref: %d\u2013%d",abs(RIDGE_LON_MIN),abs(RIDGE_LON_MAX),RIDGE_LAT_MIN,RIDGE_LAT_MAX,CLIM_START,CLIM_END),
                                def,thr_lbl,lbl)
      ggsave(file.path(lag_fig_dir,"04_z500_timeseries_ridge.png"),p04,width=FIGURE_WIDTH_WIDE,height=7,dpi=FIGURE_DPI)
      
      # Fig 05 — SLP monthly composites
      cat("05 ")
      slp_mc       <- build_monthly_composite_lag(slp_anom, date_index, drought_mask_base, lag_k)
      slp_clim_lim <- ceiling(quantile(abs(slp_mc$value),0.98,na.rm=TRUE)/0.5)*0.5
      p05 <- ggplot(slp_mc, aes(x=x,y=y,fill=value))+geom_raster()+
        geom_sf(data=map_layers$borders,fill=NA,color="black",linewidth=0.3,inherit.aes=FALSE)+
        geom_sf(data=map_layers$lakes,fill="white",color="grey60",linewidth=0.2,inherit.aes=FALSE)+
        scale_fill_distiller(palette=ANOMALY_PALETTE,limits=c(-slp_clim_lim,slp_clim_lim),oob=squish,name="SLP Anom\n(hPa)",na.value="grey80")+
        coord_sf(xlim=c(-180,-100),ylim=c(40,70),expand=FALSE)+facet_wrap(~panel_title,ncol=4)+
        labs(title=sprintf("SLP Anomaly \u2014 %s Drought Months [%s]",lbl,lag_str),
             subtitle=sprintf("Field at %s. Def: %s. Ref: %d\u2013%d.",lag_str,thr_lbl,CLIM_START,CLIM_END),
             x=NULL,y=NULL)+map_theme(10)
      ggsave(file.path(lag_fig_dir,"05_slp_monthly_anomaly_composites.png"),p05,width=FIGURE_WIDTH_WIDE,height=12,dpi=FIGURE_DPI)
      
      # Fig 06 — SLP drought vs non-drought
      cat("06 ")
      slp_lim <- ceiling(quantile(abs(terra::values(slp_diff)),0.98,na.rm=TRUE)/0.5)*0.5
      p06 <- ggplot(diff_panels(slp_dr,slp_ndr,slp_diff),aes(x=x,y=y,fill=value))+geom_raster()+
        geom_sf(data=map_layers$borders,fill=NA,color="black",linewidth=0.3,inherit.aes=FALSE)+
        scale_fill_distiller(palette=ANOMALY_PALETTE,limits=c(-slp_lim,slp_lim),oob=squish,name="SLP Anom\n(hPa)",na.value="grey80")+
        coord_sf(xlim=c(-180,-100),ylim=c(40,70),expand=FALSE)+facet_wrap(~panel,ncol=3)+
        labs(title=sprintf("SLP: %s Drought vs. Non-drought [%s]",lbl,lag_str),
             subtitle=sprintf("%s | %s | %s | Ref: %d\u2013%d",lbl,lag_str,thr_lbl,CLIM_START,CLIM_END),
             x=NULL,y=NULL)+map_theme(11)
      ggsave(file.path(lag_fig_dir,"06_slp_drought_vs_nondrought.png"),p06,width=FIGURE_WIDTH_WIDE,height=6,dpi=FIGURE_DPI)
      
      # Fig 08 — SLP time series
      cat("08 ")
      p08 <- plot_dual_panel_ts(date_index_lag$date,date_index_lag$slp_nwbc,date_index_lag$index_value,
                                drought_mask_base,sprintf("SLP Anom (hPa)\n%d\u2013%d W [%s]",abs(RIDGE_LON_MIN),abs(RIDGE_LON_MAX),lag_short),lbl,
                                sprintf("SLP Anomaly vs. %s [%s]",lbl,lag_str),
                                sprintf("NW BC box: %d\u2013%d W, %d\u2013%d N | Ref: %d\u2013%d",abs(RIDGE_LON_MIN),abs(RIDGE_LON_MAX),RIDGE_LAT_MIN,RIDGE_LAT_MAX,CLIM_START,CLIM_END),
                                def,thr_lbl)
      ggsave(file.path(lag_fig_dir,"08_slp_timeseries.png"),p08,width=FIGURE_WIDTH_WIDE,height=7,dpi=FIGURE_DPI)
      
      # Fig 10 — SST monthly composites
      cat("10 ")
      sst_mc   <- build_monthly_composite_lag(sst_anom, date_index, drought_mask_base, lag_k)
      sst_mlim <- ceiling(quantile(abs(sst_mc$value),0.98,na.rm=TRUE)*2)/2
      p10 <- ggplot(sst_mc, aes(x=x,y=y,fill=value))+geom_raster()+
        geom_sf(data=map_layers$borders,fill="grey85",color="black",linewidth=0.3,inherit.aes=FALSE)+
        scale_fill_distiller(palette=SST_PALETTE,limits=c(-sst_mlim,sst_mlim),oob=squish,name="SST Anom\n(deg C)",na.value="grey80")+
        coord_sf(xlim=c(-180,-110),ylim=c(20,65),expand=FALSE)+facet_wrap(~panel_title,ncol=4)+
        labs(title=sprintf("NE Pacific SST Anomaly \u2014 %s Drought Months [%s]",lbl,lag_str),
             subtitle=sprintf("Ocean state at %s. Def: %s. Ref: %d\u2013%d.",lag_str,thr_lbl,CLIM_START,CLIM_END),
             x=NULL,y=NULL)+map_theme(10)
      ggsave(file.path(lag_fig_dir,"10_sst_monthly_composites_drought.png"),p10,width=FIGURE_WIDTH_WIDE,height=12,dpi=FIGURE_DPI)
      
      # Fig 12 — Joint scatter
      cat("12 ")
      sc_df <- date_index_lag %>%
        filter(!is.na(index_value),!is.na(z500_ridge),!is.na(slp_nwbc),!is.na(sst_nepac)) %>%
        mutate(season=factor(season,levels=c("DJF","MAM","JJA","SON")))
      mk_scatter <- function(xvar, xlabel, title_txt) {
        r <- round(cor(sc_df[[xvar]],sc_df$index_value,use="complete.obs"),3)
        ggplot(sc_df,aes(x=.data[[xvar]],y=index_value,color=season))+
          geom_hline(yintercept=if(def$classify=="single") def$threshold else def$onset,
                     linetype="dashed",color="darkred",linewidth=0.5)+
          geom_hline(yintercept=0,color="grey50",linewidth=0.3)+geom_vline(xintercept=0,color="grey50",linewidth=0.3)+
          geom_point(alpha=0.5,size=1.2)+geom_smooth(method="lm",se=TRUE,color="black",linewidth=0.8,fill="grey70",alpha=0.3)+
          scale_color_brewer(palette="Set1",name="Season")+
          labs(title=title_txt,x=xlabel,y=sprintf("%s (t)",lbl),subtitle=sprintf("r = %.3f  (n = %d)",r,nrow(sc_df)))+
          theme_bw(base_size=10)+theme(plot.title=element_text(face="bold",size=10),legend.position="bottom")
      }
      p12 <- (mk_scatter("z500_ridge",sprintf("Z500 Anom (m) at t%s",if(lag_k==0)"" else sprintf("-%d",lag_k)),"Z500 Ridge vs. Index")|
                mk_scatter("slp_nwbc",  sprintf("SLP Anom (hPa) at t%s",if(lag_k==0)"" else sprintf("-%d",lag_k)),"SLP vs. Index")|
                mk_scatter("sst_nepac", sprintf("SST Anom (degC) at t%s",if(lag_k==0)"" else sprintf("-%d",lag_k)),"NE Pacific SST vs. Index"))+
        plot_annotation(title=sprintf("Atmospheric/Ocean Anomalies vs. %s [%s]",lbl,lag_str),
                        caption=sprintf("Dashed = drought onset. Def: %s.",thr_lbl),
                        theme=theme(plot.title=element_text(face="bold",size=13)))
      ggsave(file.path(lag_fig_dir,"12_joint_scatter.png"),p12,width=FIGURE_WIDTH_WIDE,height=6,dpi=FIGURE_DPI)
      
      # Fig 13 — Summary composite panel
      cat("13 ")
      mk_diff_map <- function(r,lim,title_txt,unit,xlim,ylim) {
        df <- as.data.frame(r,xy=TRUE); names(df)[3] <- "value"
        ggplot(df,aes(x=x,y=y,fill=value))+geom_raster()+
          geom_sf(data=map_layers$borders,fill=NA,color="black",linewidth=0.25,inherit.aes=FALSE)+
          scale_fill_distiller(palette="RdBu",limits=c(-lim,lim),oob=squish,name=unit,na.value="grey80")+
          coord_sf(xlim=xlim,ylim=ylim,expand=FALSE)+labs(title=title_txt,x=NULL,y=NULL)+map_theme(10)
      }
      lim_z  <- ceiling(quantile(abs(terra::values(z500_diff)),0.98,na.rm=TRUE)/10)*10
      lim_sl <- ceiling(quantile(abs(terra::values(slp_diff)), 0.98,na.rm=TRUE)/0.5)*0.5
      lim_ss <- ceiling(quantile(abs(terra::values(sst_diff)), 0.98,na.rm=TRUE)*2)/2
      p13 <- (mk_diff_map(z500_diff,lim_z, sprintf("Z500 (m): %s Drought minus Non-drought [%s]",lbl,lag_str),"Z500 (m)",c(-180,-100),c(40,70))/
                mk_diff_map(slp_diff, lim_sl,sprintf("SLP (hPa): %s Drought minus Non-drought [%s]",lbl,lag_str),"SLP (hPa)",c(-180,-100),c(40,70))/
                mk_diff_map(sst_diff, lim_ss,sprintf("SST (\u00b0C): %s Drought minus Non-drought [%s]",lbl,lag_str),"SST (\u00b0C)",c(-180,-110),c(20,65)))+
        plot_annotation(title=sprintf("Composite Anomalies: %s Drought [%s]",lbl,lag_str),
                        subtitle=sprintf("%s | %s  (n=%d vs n=%d)  | Ref: %d\u2013%d",thr_lbl,lag_str,n_dr,n_ndr,CLIM_START,CLIM_END),
                        theme=theme(plot.title=element_text(face="bold",size=14),plot.subtitle=element_text(size=10)))
      ggsave(file.path(lag_fig_dir,"13_composite_summary_panel.png"),p13,width=FIGURE_WIDTH_STD,height=16,dpi=FIGURE_DPI)
      
      # Tables T1-T3
      cat("T1-T3 ")
      t1 <- date_index_lag %>%
        filter(!is.na(index_value)) %>%
        mutate(dclass=ifelse(is_drought %in% TRUE,sprintf("Drought [%s]",thr_lbl),"Non-drought")) %>%
        group_by(dclass) %>%
        summarise(n=n(),index_mean=round(mean(index_value,na.rm=TRUE),3),
                  z500_mean=round(mean(z500_ridge,na.rm=TRUE),1),slp_mean=round(mean(slp_nwbc,na.rm=TRUE),2),
                  sst_mean=round(mean(sst_nepac,na.rm=TRUE),3),.groups="drop")
      write_csv(t1, file.path(lag_tbl_dir,"T1_drought_month_summary.csv"))
      
      t2 <- date_index_lag %>%
        filter(!is.na(index_value),is_drought %in% TRUE) %>%
        group_by(season) %>%
        summarise(n=n(),index_mean=round(mean(index_value,na.rm=TRUE),3),
                  z500_mean=round(mean(z500_ridge,na.rm=TRUE),1),slp_mean=round(mean(slp_nwbc,na.rm=TRUE),2),
                  sst_mean=round(mean(sst_nepac,na.rm=TRUE),3),.groups="drop") %>%
        arrange(match(season,c("DJF","MAM","JJA","SON")))
      write_csv(t2, file.path(lag_tbl_dir,"T2_seasonal_anomaly_statistics.csv"))
      
      cor_data <- date_index_lag %>%
        filter(!is.na(index_value),!is.na(z500_ridge),!is.na(slp_nwbc),!is.na(sst_nepac)) %>%
        select(index_value,z500_ridge,slp_nwbc,sst_nepac)
      n_obs  <- nrow(cor_data)
      cor_r  <- cor(cor_data,use="complete.obs")
      t_stat <- cor_r*sqrt((n_obs-2)/(1-cor_r^2))
      p_mat  <- 2*pt(-abs(t_stat),df=n_obs-2)
      cor_labels <- c(lbl,sprintf("Z500 Ridge (m) [%s]",lag_short),
                      sprintf("SLP NW-BC (hPa) [%s]",lag_short),sprintf("SST NE-Pac (degC) [%s]",lag_short))
      cor_out  <- cbind(Variable=cor_labels,as.data.frame(round(cor_r,3)));  names(cor_out)[-1] <- cor_labels
      pval_out <- cbind(Variable=cor_labels,as.data.frame(round(p_mat,4))); names(pval_out)[-1] <- cor_labels
      write_csv(cor_out,  file.path(lag_tbl_dir,"T3_correlation_matrix.csv"))
      write_csv(pval_out, file.path(lag_tbl_dir,"T3_pvalue_matrix.csv"))
      
      lag_cor_rows[[length(lag_cor_rows)+1L]] <- data.frame(
        index=lbl,definition=def$short_label,lag_k=lag_k,lag_label=lag_label,n_obs=n_obs,
        r_z500=round(cor_r["index_value","z500_ridge"],4),r_slp=round(cor_r["index_value","slp_nwbc"],4),
        r_sst=round(cor_r["index_value","sst_nepac"],4),
        p_z500=round(p_mat["index_value","z500_ridge"],5),p_slp=round(p_mat["index_value","slp_nwbc"],5),
        p_sst=round(p_mat["index_value","sst_nepac"],5),
        abs_r_z500=abs(round(cor_r["index_value","z500_ridge"],4)),
        abs_r_slp=abs(round(cor_r["index_value","slp_nwbc"],4)),
        abs_r_sst=abs(round(cor_r["index_value","sst_nepac"],4)),
        stringsAsFactors=FALSE)
      
      cat(sprintf("Done [%s]\n", lag_label))
    }  # end lag loop
    
    # Tables T4, T5
    lag_rows_df <- bind_rows(lag_cor_rows)
    t4 <- lag_rows_df %>% select(index,definition,lag_k,lag_label,n_obs,r_z500,p_z500,r_slp,p_slp,r_sst,p_sst)
    write_csv(t4, file.path(tbl_dir,"T4_lag_correlation_profile.csv"))
    
    find_best <- function(rows, var_r, var_p, var_name) {
      best_row <- rows[which.max(rows[[var_r]]),]
      data.frame(index=lbl,definition=def$short_label,variable=var_name,
                 best_lag_k=best_row$lag_k,best_lag_label=best_row$lag_label,
                 r=best_row[[sub("abs_r_","r_",var_r)]],abs_r=best_row[[var_r]],
                 p_value=best_row[[var_p]],n_obs=best_row$n_obs,stringsAsFactors=FALSE)
    }
    t5 <- bind_rows(find_best(lag_rows_df,"abs_r_z500","p_z500","Z500_Ridge"),
                    find_best(lag_rows_df,"abs_r_slp", "p_slp", "SLP_NWBC"),
                    find_best(lag_rows_df,"abs_r_sst", "p_sst", "SST_NEPac"))
    write_csv(t5, file.path(tbl_dir,"T5_best_lag_per_variable.csv"))
    
    def_best_lag_rows[[length(def_best_lag_rows)+1L]] <- t5
    cat(sprintf("   Done: %s  [%s]\n\n", lbl, def$short_label))
  }  # end per-index loop
  
  # T6: all-index best-lag summary
  t6 <- bind_rows(def_best_lag_rows) %>% arrange(variable, index)
  write_csv(t6, file.path(ensure_dir(OUT_DIR,def$short_label,"Tables"),"T6_best_lag_all_indices.csv"))
  cat(sprintf("  \u2713 Definition [%s] complete.  T6 saved.\n\n", def$short_label))
}

# ==============================================================================
#  MODULE B1: HISTORICAL PDO-STRATIFIED COMPOSITE ANALYSIS
# ==============================================================================
cat("==============================================================\n")
cat(" MODULE B1: Historical PDO-stratified composites\n")
cat("==============================================================\n")
out_b1 <- "composite_results"
dir.create(out_b1, showWarnings = FALSE)

# Build consensus drought months from the SW T0 (>=2 indices in agreement)
t0_primary <- t0_sw %>%
  filter(index %in% c("SPEI6","SPEI12","SPI6","SPI12")) %>%
  mutate(date_ym = as.Date(paste0(date, "-01")))

drought_consensus <- t0_primary %>%
  group_by(date) %>%
  summarise(n_indices  = n_distinct(index),
            mean_value = mean(index_value, na.rm = TRUE),
            year       = first(year),
            month_num  = first(month_num),
            season     = first(season), .groups = "drop") %>%
  filter(n_indices >= 2) %>%
  mutate(date_ym = as.Date(paste0(date, "-01")),
         is_drought_period = year >= 2022 & year <= 2025)

# PDO proxy from NE Pacific SST box
pdo_proxy <- terra::global(crop(sst_anom, ext(-160,-130,40,55)), "mean", na.rm=TRUE)[,1]
pdo_df    <- data.frame(date_ym = as.Date(terra::time(sst_anom)), pdo_proxy = pdo_proxy)
drought_consensus <- drought_consensus %>%
  left_join(pdo_df, by = "date_ym") %>%
  mutate(pdo_phase = ifelse(pdo_proxy > 0, "PDO+", "PDO-"))

all_drought_dates  <- drought_consensus$date_ym
drought2025_dates  <- drought_consensus$date_ym[ drought_consensus$is_drought_period]
historical_dates   <- drought_consensus$date_ym[!drought_consensus$is_drought_period]
pdo_pos_dates      <- drought_consensus$date_ym[ drought_consensus$pdo_phase=="PDO+" & !drought_consensus$is_drought_period]
pdo_neg_dates      <- drought_consensus$date_ym[ drought_consensus$pdo_phase=="PDO-" & !drought_consensus$is_drought_period]

cat(sprintf("  Consensus drought months: %d total (%d in 2022-2025)\n",
            nrow(drought_consensus), sum(drought_consensus$is_drought_period)))
cat(sprintf("  Historical: PDO+ = %d, PDO- = %d\n",
            length(pdo_pos_dates), length(pdo_neg_dates)))

composite_months <- function(anom_rast, target_dates, label) {
  # ---- Robust date extraction: mirrors the fallback in compute_anomaly() ----
  # terra::time() can return NULL or invalid dates after a writeCDF/rast
  # round-trip when the NetCDF time axis lacks explicit CF metadata.
  # In that case fall back to the same positional Jan-START_YEAR sequence
  # that was used when the anomaly was originally computed in w9.
  n_lyr   <- terra::nlyr(anom_rast)
  dates_r <- tryCatch({
    d <- as.Date(terra::time(anom_rast))
    if (is.null(d) || length(d) != n_lyr || all(is.na(d)) ||
        min(d, na.rm = TRUE) < as.Date("1900-01-01") ||
        max(d, na.rm = TRUE) > as.Date("2030-01-01")) NULL else d
  }, error = function(e) NULL)
  
  if (is.null(dates_r)) {
    dates_r <- seq(as.Date(sprintf("%d-01-01", START_YEAR)),
                   by = "month", length.out = n_lyr)
  }
  
  idx <- which(dates_r %in% target_dates)
  if (!length(idx)) {
    warning(sprintf(
      paste0("composite_months: no matching layers found for '%s'. ",
             "target range: %s to %s | raster range: %s to %s"),
      label,
      format(min(target_dates)), format(max(target_dates)),
      format(min(dates_r)),      format(max(dates_r))
    ))
    return(NULL)
  }
  comp <- terra::app(anom_rast[[idx]], mean, na.rm = TRUE)
  names(comp) <- label
  comp
}

z500_all_hist  <- composite_months(z500_anom, historical_dates,  "Z500_all_historical")
z500_pdo_pos   <- composite_months(z500_anom, pdo_pos_dates,     "Z500_PDO+")
z500_pdo_neg   <- composite_months(z500_anom, pdo_neg_dates,     "Z500_PDO-")
z500_2022_2025 <- composite_months(z500_anom, drought2025_dates, "Z500_2022-2025")
slp_all_hist   <- composite_months(slp_anom,  historical_dates,  "SLP_all_historical")
slp_pdo_pos    <- composite_months(slp_anom,  pdo_pos_dates,     "SLP_PDO+")
slp_pdo_neg    <- composite_months(slp_anom,  pdo_neg_dates,     "SLP_PDO-")
sst_all_hist   <- composite_months(sst_anom,  historical_dates,  "SST_all_historical")
sst_pdo_pos    <- composite_months(sst_anom,  pdo_pos_dates,     "SST_PDO+")
sst_2022_2025  <- composite_months(sst_anom,  drought2025_dates, "SST_2022-2025")

spatial_cor <- function(r1, r2) {
  # NULL guard: composite_months() returns NULL when no matching layers found.
  if (is.null(r1) || is.null(r2)) {
    warning("spatial_cor: one or both inputs are NULL -- returning NA.")
    return(NA_real_)
  }
  v1 <- as.vector(terra::values(r1))
  v2 <- as.vector(terra::values(r2))
  ok <- !is.na(v1) & !is.na(v2)
  if (sum(ok) < 10) return(NA_real_)
  cor(v1[ok], v2[ok])
}

cat("\n  Spatial correlations (2022-2025 vs historical composites):\n")
cat(sprintf("    Z500 vs All Hist: %.3f\n", spatial_cor(z500_2022_2025, z500_all_hist)))
cat(sprintf("    Z500 vs PDO+    : %.3f\n", spatial_cor(z500_2022_2025, z500_pdo_pos)))
cat(sprintf("    Z500 vs PDO-    : %.3f\n", spatial_cor(z500_2022_2025, z500_pdo_neg)))
cat(sprintf("    SLP  vs All Hist: %.3f\n", spatial_cor(composite_months(slp_anom,drought2025_dates,"SLP_2022-2025"), slp_all_hist)))
cat(sprintf("    SLP  vs PDO+    : %.3f\n", spatial_cor(composite_months(slp_anom,drought2025_dates,"SLP_2022-2025"), slp_pdo_pos)))

severity_by_year <- drought_consensus %>%
  group_by(year) %>%
  summarise(n_drought_months = n(), mean_severity = mean(abs(mean_value),na.rm=TRUE),
            composite_score = n_drought_months * mean_severity, .groups = "drop") %>%
  arrange(desc(composite_score)) %>% mutate(rank = row_number())

write.csv(severity_by_year, file.path(out_b1,"drought_severity_ranking.csv"), row.names=FALSE)
for (obj in list(z500_all_hist,z500_pdo_pos,z500_pdo_neg,z500_2022_2025,slp_all_hist,slp_pdo_pos,sst_all_hist,sst_pdo_pos,sst_2022_2025)) {
  if (!is.null(obj)) terra::writeRaster(obj, file.path(out_b1, sprintf("%s.tif", names(obj))), overwrite=TRUE)
}
cat("  \u2713 B1 composites saved to:", out_b1, "\n\n")

# ==============================================================================
#  MODULE B3: BOOTSTRAP FINGERPRINT SIGNIFICANCE (PARALLELISED)
# ==============================================================================
cat("==============================================================\n")
cat(" MODULE B3: Bootstrap fingerprint significance\n")
cat("==============================================================\n")
out_b3     <- "bootstrap_results"
dir.create(out_b3, showWarnings = FALSE)
N_BOOT     <- 5000L
BLOCK_SIZE <- 6L
N_DROUGHT  <- 48L
set.seed(42L)

dates_era5 <- as.Date(terra::time(z500_anom))
year_nums  <- as.integer(format(dates_era5, "%Y"))

bg_idx      <- which(year_nums >= 1950 & year_nums <= 2021)
drought_idx <- which(year_nums >= 2022 & year_nums <= 2025)

if (length(bg_idx) < 100 || length(drought_idx) != N_DROUGHT) {
  cat("  Insufficient data for bootstrap \u2014 skipping B3.\n")
} else {
  z500_bg_mat <- terra::values(z500_anom[[bg_idx]])
  slp_bg_mat  <- terra::values(slp_anom[[bg_idx]])
  sst_bg_mat  <- terra::values(sst_anom[[bg_idx]])
  
  obs_z500     <- terra::app(z500_anom[[drought_idx]], mean, na.rm = TRUE)
  obs_slp      <- terra::app(slp_anom[[drought_idx]],  mean, na.rm = TRUE)
  obs_sst      <- terra::app(sst_anom[[drought_idx]],  mean, na.rm = TRUE)
  obs_z500_vec <- as.vector(terra::values(obs_z500))
  obs_slp_vec  <- as.vector(terra::values(obs_slp))
  obs_sst_vec  <- as.vector(terra::values(obs_sst))
  
  n_bg            <- length(bg_idx)
  n_pix           <- nrow(z500_bg_mat)
  block_starts    <- seq_len(n_bg - BLOCK_SIZE + 1L)
  n_blocks_needed <- ceiling(N_DROUGHT / BLOCK_SIZE)
  
  n_cores  <- max(1L, min(detectCores() - 1L, 4L))
  chunk_sz <- ceiling(N_BOOT / n_cores)
  chunks   <- split(seq_len(N_BOOT), ceiling(seq_len(N_BOOT) / chunk_sz))
  
  cat(sprintf("  Running %d bootstrap reps across %d core(s)...\n", N_BOOT, n_cores))
  
  run_chunk <- function(boot_ids) {
    cnt_z <- cnt_s <- cnt_sst <- numeric(n_pix)
    set.seed(boot_ids[1])
    for (b in boot_ids) {
      starts    <- sample(block_starts, n_blocks_needed, replace = TRUE)
      drawn_idx <- pmin(unlist(lapply(starts, function(s) s:(s + BLOCK_SIZE - 1L))), n_bg)[seq_len(N_DROUGHT)]
      cnt_z   <- cnt_z   + (abs(rowMeans(z500_bg_mat[, drawn_idx, drop=FALSE], na.rm=TRUE)) >= abs(obs_z500_vec))
      cnt_s   <- cnt_s   + (abs(rowMeans(slp_bg_mat[,  drawn_idx, drop=FALSE], na.rm=TRUE)) >= abs(obs_slp_vec))
      cnt_sst <- cnt_sst + (abs(rowMeans(sst_bg_mat[,  drawn_idx, drop=FALSE], na.rm=TRUE)) >= abs(obs_sst_vec))
    }
    list(z = cnt_z, s = cnt_s, sst = cnt_sst)
  }
  
  if (n_cores > 1L) {
    cl <- makeCluster(n_cores, type = "PSOCK")
    clusterExport(cl, c("z500_bg_mat","slp_bg_mat","sst_bg_mat","block_starts",
                        "n_blocks_needed","N_DROUGHT","n_bg","BLOCK_SIZE",
                        "obs_z500_vec","obs_slp_vec","obs_sst_vec","n_pix"), envir = environment())
    results <- parLapply(cl, chunks, run_chunk)
    stopCluster(cl)
  } else {
    results <- lapply(chunks, run_chunk)
  }
  
  pval_z500 <- Reduce(`+`, lapply(results, `[[`, "z"))   / N_BOOT
  pval_slp  <- Reduce(`+`, lapply(results, `[[`, "s"))   / N_BOOT
  pval_sst  <- Reduce(`+`, lapply(results, `[[`, "sst")) / N_BOOT
  
  fdr_correct <- function(pvals, alpha = 0.05) {
    valid <- !is.na(pvals); adj <- rep(NA_real_, length(pvals))
    adj[valid] <- p.adjust(pvals[valid], method = "BH"); adj < alpha
  }
  sig_z500 <- fdr_correct(pval_z500)
  sig_slp  <- fdr_correct(pval_slp)
  sig_sst  <- fdr_correct(pval_sst)
  
  cat(sprintf("  Fraction significant (FDR q<0.05): Z500=%.1f%%, SLP=%.1f%%, SST=%.1f%%\n",
              100*mean(sig_z500,na.rm=TRUE),100*mean(sig_slp,na.rm=TRUE),100*mean(sig_sst,na.rm=TRUE)))
  
  make_pval_rast <- function(pv, tmpl) { r <- tmpl; terra::values(r) <- pv; r }
  for (pair in list(list(pval_z500,"pval_z500_bootstrap.tif"),
                    list(pval_slp, "pval_slp_bootstrap.tif"),
                    list(pval_sst, "pval_sst_bootstrap.tif"),
                    list(as.numeric(sig_z500),"sig_z500_fdr005.tif"),
                    list(as.numeric(sig_slp), "sig_slp_fdr005.tif"),
                    list(as.numeric(sig_sst), "sig_sst_fdr005.tif"))) {
    tmpl <- if (grepl("z500",pair[[2]])) obs_z500 else if (grepl("slp",pair[[2]])) obs_slp else obs_sst
    terra::writeRaster(make_pval_rast(pair[[1]], tmpl), file.path(out_b3, pair[[2]]), overwrite=TRUE)
  }
  
  write.csv(data.frame(
    variable            = c("Z500","SLP","SST"),
    pct_significant_fdr = c(100*mean(sig_z500,na.rm=TRUE),100*mean(sig_slp,na.rm=TRUE),100*mean(sig_sst,na.rm=TRUE)),
    N_bootstrap = N_BOOT, block_size_months = BLOCK_SIZE, drought_months = N_DROUGHT
  ), file.path(out_b3,"bootstrap_significance_summary.csv"), row.names=FALSE)
  
  # Stippled significance maps
  plot_with_stippling <- function(anom_rast, sig_rast, variable, units, col_lim) {
    anom_df <- as.data.frame(anom_rast, xy=TRUE); colnames(anom_df)[3] <- "value"
    sig_df  <- as.data.frame(sig_rast,  xy=TRUE); colnames(sig_df)[3]  <- "sig"
    ggplot() +
      geom_raster(data=anom_df,aes(x=x,y=y,fill=value))+
      geom_point(data=sig_df %>% filter(sig==1),aes(x=x,y=y),size=0.3,colour="black",alpha=0.5)+
      scale_fill_distiller(palette="RdBu",direction=-1,limits=c(-col_lim,col_lim),oob=squish,
                           name=sprintf("%s\n(%s)",variable,units))+
      coord_fixed()+
      labs(title=sprintf("2022\u20132025 %s Anomaly",variable),
           subtitle=sprintf("Dots = FDR-corrected significant (q<0.05, N=%d bootstrap)",N_BOOT),
           x="Longitude",y="Latitude")+theme_minimal(base_size=11)
  }
  
  sig_z500_rast <- make_pval_rast(as.numeric(sig_z500), obs_z500)
  sig_slp_rast  <- make_pval_rast(as.numeric(sig_slp),  obs_slp)
  sig_sst_rast  <- make_pval_rast(as.numeric(sig_sst),  obs_sst)
  
  pdf(file.path(out_b3,"fig_bootstrap_significance_maps.pdf"), width=14, height=12)
  gridExtra::grid.arrange(
    plot_with_stippling(obs_z500, sig_z500_rast, "Z500", "m",   50),
    plot_with_stippling(obs_slp,  sig_slp_rast,  "SLP",  "hPa",  2),
    plot_with_stippling(obs_sst,  sig_sst_rast,  "SST",  "\u00b0C", 1.5),
    ncol = 1, top = "2022\u20132025 Drought Fingerprint: Bootstrap Significance (FDR q<0.05)")
  dev.off()
  
  # Historical analogues
  cat("  Computing historical analogues...\n")
  analogue_r <- vapply(seq_along(bg_idx), function(i) {
    bg_z500 <- z500_bg_mat[, i]
    ok      <- !is.na(bg_z500) & !is.na(obs_z500_vec)
    if (sum(ok) < 50) return(NA_real_)
    cor(bg_z500[ok], obs_z500_vec[ok])
  }, numeric(1))
  data.frame(date=dates_era5[bg_idx], year=year_nums[bg_idx], spatial_r=analogue_r) %>%
    arrange(desc(spatial_r)) %>% mutate(rank=row_number()) %>%
    write.csv(file.path(out_b3,"historical_analogues_z500.csv"), row.names=FALSE)
  
  cat("  \u2713 B3 bootstrap analysis complete.\n")
}

# ==============================================================================
# Following Extracts numerical values from the 2022-2025 composite rasters for §4.7
# ==============================================================================

cat("\n===== ATMOSPHERIC COMPOSITE MS BLANKS =====\n")

# ── 2022-2025 full-period mean composite rasters ─────────────────────────────
z500_comp <- terra::app(terra::subset(z500_anom, focus_idx), mean, na.rm = TRUE)
slp_comp  <- terra::app(terra::subset(slp_anom,  focus_idx), mean, na.rm = TRUE)
sst_comp  <- terra::app(terra::subset(sst_anom,  focus_idx), mean, na.rm = TRUE)

# ── Z500: ridge maximum and trough minimum ────────────────────────────────────
z500_vals  <- terra::values(z500_comp, mat = FALSE)
z500_max   <- max(z500_vals,  na.rm = TRUE)
z500_min   <- min(z500_vals,  na.rm = TRUE)
z500_mc    <- which.max(z500_vals)
z500_tc    <- which.min(z500_vals)
z500_ridge_xy <- terra::xyFromCell(z500_comp, z500_mc)
z500_trough_xy<- terra::xyFromCell(z500_comp, z500_tc)

cat(sprintf("  Z500 ridge max : %+.1f m at lon=%.1f, lat=%.1f\n",
            z500_max, z500_ridge_xy[1], z500_ridge_xy[2]))
cat(sprintf("  Z500 trough min: %+.1f m at lon=%.1f, lat=%.1f\n",
            z500_min, z500_trough_xy[1], z500_trough_xy[2]))

# ── SLP: maximum ─────────────────────────────────────────────────────────────
slp_vals  <- terra::values(slp_comp, mat = FALSE)
slp_max   <- max(slp_vals, na.rm = TRUE)
slp_mc    <- which.max(slp_vals)
slp_max_xy<- terra::xyFromCell(slp_comp, slp_mc)

cat(sprintf("  SLP max: %+.2f hPa at lon=%.1f, lat=%.1f\n",
            slp_max, slp_max_xy[1], slp_max_xy[2]))

# ── SST: NE Pacific warm anomaly (40-60N, 150-120W) ──────────────────────────
ne_pac     <- terra::ext(-150, -120, 40, 60)
sst_ne     <- terra::crop(sst_comp, ne_pac)
sst_ne_mean<- terra::global(sst_ne, "mean", na.rm = TRUE)$mean
sst_ne_max <- terra::global(sst_ne, "max",  na.rm = TRUE)$max

# Basin-area mean SST (for text: "warm anomalies of ___ °C along NE Pacific")
cat(sprintf("  SST NE Pac mean: %+.2f deg C\n", sst_ne_mean))
cat(sprintf("  SST NE Pac max : %+.2f deg C\n", sst_ne_max))

# ── JJA-only composites (more relevant for summer drought discussion) ─────────
jja_idx    <- focus_idx[base_date_index$month[focus_idx] %in% c(6, 7, 8)]
z500_jja   <- terra::app(terra::subset(z500_anom, jja_idx), mean, na.rm = TRUE)
slp_jja    <- terra::app(terra::subset(slp_anom,  jja_idx), mean, na.rm = TRUE)

z500_jja_max  <- max(terra::values(z500_jja), na.rm = TRUE)
slp_jja_max   <- max(terra::values(slp_jja),  na.rm = TRUE)
z500_jja_mc   <- which.max(terra::values(z500_jja))
z500_jja_xy   <- terra::xyFromCell(z500_jja, z500_jja_mc)
slp_jja_mc    <- which.max(terra::values(slp_jja))
slp_jja_xy    <- terra::xyFromCell(slp_jja, slp_jja_mc)

cat(sprintf("  JJA Z500 max: %+.1f m at lon=%.1f, lat=%.1f\n",
            z500_jja_max, z500_jja_xy[1], z500_jja_xy[2]))
cat(sprintf("  JJA SLP max : %+.2f hPa at lon=%.1f, lat=%.1f\n",
            slp_jja_max, slp_jja_xy[1], slp_jja_xy[2]))

# ── IVT reduction (§4.7): if IVT loaded use it; otherwise note ───────────────
# IVT is not loaded in this script; the §4.7 text can use the area-averaged
# Z500/SLP anomalies above as proxy, or fill with published values.
# If you load IVT separately, compute here the same way.

# ── Write to CSV ──────────────────────────────────────────────────────────────
atm_blanks <- data.frame(
  variable = c(
    "Z500_annual_max_m",   "Z500_ridge_lon",  "Z500_ridge_lat",
    "Z500_annual_min_m",   "Z500_trough_lon", "Z500_trough_lat",
    "SLP_annual_max_hPa",  "SLP_max_lon",     "SLP_max_lat",
    "SST_NEPac_mean_C",    "SST_NEPac_max_C",
    "Z500_JJA_max_m",      "Z500_JJA_lon",    "Z500_JJA_lat",
    "SLP_JJA_max_hPa",     "SLP_JJA_lon",     "SLP_JJA_lat"
  ),
  value = c(
    z500_max, z500_ridge_xy[1],  z500_ridge_xy[2],
    z500_min, z500_trough_xy[1], z500_trough_xy[2],
    slp_max,  slp_max_xy[1],     slp_max_xy[2],
    sst_ne_mean, sst_ne_max,
    z500_jja_max, z500_jja_xy[1], z500_jja_xy[2],
    slp_jja_max,  slp_jja_xy[1],  slp_jja_xy[2]
  )
)
print(atm_blanks)
write.csv(atm_blanks,
          file.path(OUT_DIR, "ms1_atmospheric_blanks.csv"), row.names = FALSE)
cat(sprintf("  -> %s\n",
            file.path(OUT_DIR, "ms1_atmospheric_blanks.csv")))
cat("===== DONE =====\n")
# ==============================================================================
#  DONE
# ==============================================================================
cat("\n==============================================================\n")
cat(" w10a COMPLETE\n")
cat("==============================================================\n")
cat(sprintf(" %s/\n", OUT_DIR))
cat("   Figures/SST/                     Figs 09, 11\n")
cat("   Figures/Drought_2022_2025/        Figs 14-21\n")
cat("   {SW_q10,Alt2T}/Figures/{INDEX}/   Figs 01-13 per index+lag\n")
cat("   {SW_q10,Alt2T}/Tables/{INDEX}/    Tables T1-T5\n")
cat("   {SW_q10,Alt2T}/Tables/            Table T6\n")
cat(sprintf(" composite_results/               B1 composites\n"))
cat(sprintf(" bootstrap_results/               B3 significance maps\n"))
cat("\n Next (independent): run w10b_index_skill.R\n")
cat("==============================================================\n")