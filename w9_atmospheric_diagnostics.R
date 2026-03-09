# ============================================================================
#  w9_atmospheric_diagnostics.R
#  Nechako Watershed Drought — Atmospheric Pattern Diagnostics
#
#  Purpose (per research lead):
#    1. 500 hPa geopotential height anomalies associated with drought —
#       month-by-month relative to 1991-2020 reference period
#    2. SLP patterns during drought months vs. background
#    3. NE Pacific SST over recent years relative to historical conditions
#    4. Joint plots linking atmospheric fields to drought severity (SPEI-6)
#
#  Outputs (all in {OUT_DIR}/):
#    Figures/
#      01_z500_monthly_anomaly_composites.png    monthly mean anom by calendar month
#      02_z500_drought_vs_nondrought.png         composite difference maps
#      03_z500_drought_severity_composites.png   stratified by SPEI quartile
#      04_z500_timeseries_nw_canada_ridge.png    area-avg Z500 anom vs SPEI-6
#      05_slp_monthly_anomaly_composites.png
#      06_slp_drought_vs_nondrought.png
#      07_slp_drought_severity_composites.png
#      08_slp_timeseries.png
#      09_sst_annual_anomaly_recent.png          recent years (2018-2025) vs clim
#      10_sst_monthly_composites_drought.png     SST during drought months
#      11_sst_timeseries_ne_pacific.png          basin-avg SST anomaly time series
#      12_joint_z500_slp_sst_drought_scatter.png scatterplots vs SPEI-6
#      13_composite_summary_panel.png            single-panel overview figure
#    Tables/
#      T1_drought_month_summary.csv              Z500/SLP/SST stats per drought event
#      T2_seasonal_anomaly_statistics.csv        seasonal breakdown by variable
#      T3_correlation_matrix.csv                 Z500/SLP/SST area-avgs vs SPEI-6
#
#  Authors: [Your Name]
#  Last updated: 2025
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
})

# ============================================================================
#  CONFIGURATION
# ============================================================================
WD_PATH  <- "D:/Nechako_Drought/Nechako/"
DATA_DIR <- "D:/Nechako_Drought/Nechako/monthly_data_direct"
OUT_DIR  <- file.path(WD_PATH, "atmospheric_diagnostics")

dir.create(file.path(OUT_DIR, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "Tables"),  recursive = TRUE, showWarnings = FALSE)

setwd(WD_PATH)

# ── Analysis period ──────────────────────────────────────────────────────────
START_YEAR <- 1950
END_YEAR   <- 2025

# ── Climatology reference ────────────────────────────────────────────────────
CLIM_START <- 1991
CLIM_END   <- 2020

# ── Drought threshold (SPEI-6) ───────────────────────────────────────────────
# Months with SPEI-6 below this are classified as drought
DROUGHT_THRESHOLD <- -0.8   # moderate drought (McKee et al. 1993)
SEVERE_THRESHOLD  <- -1.5   # severe/extreme drought

# ── "Recent" period for SST trend analysis ───────────────────────────────────
RECENT_START <- 2018   # captures post-Blob / recent marine heatwave era
RECENT_END   <- END_YEAR

# ── Nechako watershed approximate centroid and study region ─────────────────
NECHAKO_LON  <- -124.5
NECHAKO_LAT  <-   53.7

# For extracting area-averaged NW Canada ridge Z500 signal
# (region where blocking ridge anomaly would appear for BC droughts)
RIDGE_LON_MIN <- -140; RIDGE_LON_MAX <- -110
RIDGE_LAT_MIN <-   50; RIDGE_LAT_MAX <-   65

# ── NE Pacific SST averaging region ─────────────────────────────────────────
# Captures PDO/marine heatwave signal most relevant to NW BC
SST_MEAN_LON_MIN <- -160; SST_MEAN_LON_MAX <- -130
SST_MEAN_LAT_MIN <-   40; SST_MEAN_LAT_MAX <-   55

# ── Plot aesthetics ──────────────────────────────────────────────────────────
DROUGHT_COLOR     <- "#d73027"
NONDROUGHT_COLOR  <- "#4575b4"
ANOMALY_PALETTE   <- "RdBu"          # diverging, reversed for Z500/SLP
SST_PALETTE       <- "RdBu"
FIGURE_DPI        <- 200
FIGURE_WIDTH_WIDE <- 16
FIGURE_WIDTH_STD  <- 12
FIGURE_HEIGHT_STD <- 9

# Month labels
MONTH_LABELS <- month.abb

# ============================================================================
#  HELPER: COASTLINES AND WATERSHED MARKER
# ============================================================================
get_map_layers <- function() {
  coast   <- ne_coastline(scale = "medium", returnclass = "sf")
  borders <- ne_countries(scale = "medium", returnclass = "sf")
  lakes   <- ne_download(scale = "medium", type = "lakes",
                         category = "physical", returnclass = "sf")
  list(coast = coast, borders = borders, lakes = lakes)
}

# ============================================================================
#  HELPER: RASTER TO TIDY DATA FRAME FOR GGPLOT
# ============================================================================
rast_to_df <- function(r, layer_idx = NULL) {
  if (!is.null(layer_idx)) r <- subset(r, layer_idx)
  as.data.frame(r, xy = TRUE)
}

# ============================================================================
#  SECTION 1: LOAD DATA
# ============================================================================
cat("============================================================\n")
cat(" SECTION 1: Loading ERA5 fields and SPEI-6\n")
cat("============================================================\n")

# ── 1a. ERA5 fixed-period anomalies (1991-2020 reference) ────────────────────
cat("  Loading Z500 anomaly (fixed 1991-2020)... ")
z500_anom <- rast(file.path(DATA_DIR, sprintf("z500_monthly_anomaly_clim%d_%d.nc",
                                              CLIM_START, CLIM_END)))
cat(sprintf("Done. (%d layers)\n", nlyr(z500_anom)))

cat("  Loading SLP raw (hPa) and computing anomaly on-the-fly... ")
slp_raw <- rast(file.path(DATA_DIR, "slp_monthly.nc"))
cat(sprintf("Done. (%d layers)\n", nlyr(slp_raw)))

# SLP anomaly: compute from raw since calculate_anomaly = FALSE in download script
# Uses fixed 1991-2020 climatology to stay consistent with Z500
all_years_vec <- rep(START_YEAR:END_YEAR, each = 12)
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
cat("  SLP anomaly computed.\n")

cat("  Loading SST anomaly (fixed 1991-2020)... ")
sst_anom <- rast(file.path(DATA_DIR, sprintf("sst_monthly_anomaly_clim%d_%d.nc",
                                             CLIM_START, CLIM_END)))
sst_raw  <- rast(file.path(DATA_DIR, "sst_monthly.nc"))
cat(sprintf("Done. (%d layers)\n", nlyr(sst_anom)))

# ── 1b. SPEI-6 time series ────────────────────────────────────────────────────
# Loads the watershed-averaged SPEI-6 produced in w1 of the pipeline.
# Expected format: CSV with columns 'date' (YYYY-MM-DD or YYYY-MM) and 'spei6'
cat("  Loading SPEI-6... ")
spei_file <- list.files(WD_PATH, pattern = "spei.*\\.csv$",
                        recursive = TRUE, full.names = TRUE)[1]

if (is.na(spei_file)) {
  # Fallback: try RDS
  spei_rds <- list.files(WD_PATH, pattern = "spei.*\\.rds$",
                         recursive = TRUE, full.names = TRUE)[1]
  if (!is.na(spei_rds)) {
    spei_df <- readRDS(spei_rds)
  } else {
    stop("SPEI-6 file not found. Expected CSV or RDS under ", WD_PATH,
         "\n  Check w1 output directory and update spei_file path manually.")
  }
} else {
  spei_df <- read_csv(spei_file, show_col_types = FALSE)
}

# ── Standardise column names ──────────────────────────────────────────────────
names(spei_df) <- tolower(trimws(names(spei_df)))

cat(sprintf("  Columns found in SPEI file: %s\n",
            paste(names(spei_df), collapse = ", ")))

# ── Resolve SPEI-6 column ─────────────────────────────────────────────────────
# Covers: standard SPEI package names, PCIC outputs, and w1 spatial summary
# names (basin_mean, basin_median).  basin_mean is preferred over basin_median
# because it is the spatial mean of SPEI-6 across all valid watershed pixels.
SPEI_COL <- NULL   # set e.g. SPEI_COL <- "basin_mean" to force a specific column

spei_candidates <- c(
  "spei6", "spei_6", "spei.6", "spei06", "spei_06",
  "spei_scale6", "spei6_mean", "spei_6_mean", "mean_spei6",
  "basin_mean",    # w1 spatial summary output — watershed pixel mean of SPEI
  "basin_median",  # w1 spatial summary output — fallback
  "spei",          # generic (picked only if no scale-specific name exists)
  "value"          # raster::extract / terra default
)

if (!is.null(SPEI_COL)) {
  if (!SPEI_COL %in% names(spei_df))
    stop(sprintf("SPEI_COL = '%s' not found.\n  Columns: %s",
                 SPEI_COL, paste(names(spei_df), collapse = ", ")))
  spei_col_found <- SPEI_COL
  cat(sprintf("  SPEI-6 column : '%s'  [manual override]\n", spei_col_found))
} else {
  spei_col_found <- intersect(spei_candidates, names(spei_df))[1]
  if (is.na(spei_col_found)) {
    num_cols <- setdiff(names(spei_df)[sapply(spei_df, is.numeric)],
                        c("year", "month", "day", "valid_pixels",
                          "total_pixels", "coverage_pct"))
    if (length(num_cols) == 0)
      stop("Cannot identify SPEI-6 column.\n",
           "  Columns: ", paste(names(spei_df), collapse = ", "), "\n",
           "  Fix: set SPEI_COL <- '<column>' near the top of the script.")
    spei_col_found <- num_cols[1]
    warning(sprintf(
      "SPEI column not recognised by name — using '%s'. Check this is correct.\n  All columns: %s",
      spei_col_found, paste(names(spei_df), collapse = ", ")))
  }
  cat(sprintf("  SPEI-6 column : '%s'  [auto-detected]\n", spei_col_found))
}

if (spei_col_found != "spei6")
  spei_df <- rename(spei_df, spei6 = all_of(spei_col_found))

# ── Resolve & parse date column ───────────────────────────────────────────────
# Strategy:
#   1. Try common date column names
#   2. Try multiple format strings on the winner
#   3. If all formats fail, reconstruct from 'year' + 'month' columns (always present
#      in w1 spatial summary outputs)
date_candidates   <- c("date", "date_formatted", "time", "yearmon",
                       "year_month", "ym", "period")
date_col_found    <- intersect(date_candidates, names(spei_df))[1]

if (is.na(date_col_found))
  date_col_found <- names(spei_df)[
    sapply(spei_df, function(x) {
      tryCatch(!all(is.na(as.Date(as.character(head(x, 3)), tryFormats =
                                    c("%Y-%m-%d", "%Y/%m/%d", "%Y-%m")))),
               error = function(e) FALSE)
    })
  ][1]

if (!is.na(date_col_found) && date_col_found != "date")
  spei_df <- rename(spei_df, date = all_of(date_col_found))

cat(sprintf("  Date column   : '%s'\n", if (is.na(date_col_found)) "(reconstructed from year+month)" else date_col_found))

# Parse date robustly
date_formats <- c("%Y-%m-%d", "%Y/%m/%d", "%Y-%m", "%d/%m/%Y",
                  "%Y%m%d", "%Y%m")
resolved_date <- NULL

if (!is.na(date_col_found)) {
  raw_dates <- as.character(spei_df$date)
  for (fmt in date_formats) {
    trial <- suppressWarnings(as.Date(raw_dates, format = fmt))
    if (sum(!is.na(trial)) > nrow(spei_df) * 0.9) {
      resolved_date <- trial
      cat(sprintf("  Date format   : '%s'\n", fmt))
      break
    }
  }
}

if (is.null(resolved_date)) {
  # Fallback: reconstruct from year + month columns
  if (all(c("year", "month") %in% names(spei_df))) {
    resolved_date <- as.Date(sprintf("%04d-%02d-01",
                                     as.integer(spei_df$year),
                                     as.integer(spei_df$month)))
    cat("  Date format   : reconstructed from 'year' + 'month' columns\n")
  } else {
    stop("Cannot parse date column and no 'year'/'month' columns to fall back on.\n",
         "  Please ensure the SPEI file has a parseable date column.")
  }
}

spei_df$date <- resolved_date

spei_df <- spei_df %>%
  mutate(
    date  = date,
    year  = year(date),
    month = month(date),
    season = case_when(
      month %in% c(12, 1, 2)  ~ "DJF",
      month %in% c(3, 4, 5)   ~ "MAM",
      month %in% c(6, 7, 8)   ~ "JJA",
      month %in% c(9, 10, 11) ~ "SON"
    ),
    drought_class = case_when(
      spei6 <= SEVERE_THRESHOLD  ~ "Severe/Extreme",
      spei6 <= DROUGHT_THRESHOLD ~ "Moderate",
      TRUE                       ~ "Non-drought"
    ),
    is_drought = spei6 <= DROUGHT_THRESHOLD
  ) %>%
  filter(year >= START_YEAR, year <= END_YEAR) %>%
  arrange(date)

cat(sprintf("Done. %d months (%d–%d). Drought months: %d (%.1f%%)\n",
            nrow(spei_df), min(spei_df$year), max(spei_df$year),
            sum(spei_df$is_drought),
            100 * mean(spei_df$is_drought)))

# ── 1c. Build a master date index aligned to the raster layers ────────────────
# ERA5 monthly rasters: layer i corresponds to month i chronologically
all_dates <- seq.Date(as.Date(sprintf("%d-01-01", START_YEAR)),
                      as.Date(sprintf("%d-12-01", END_YEAR)),
                      by = "month")
n_total   <- length(all_dates)

# Safety check
if (nlyr(z500_anom) != n_total)
  warning(sprintf("Z500 has %d layers but expected %d. Check date alignment.",
                  nlyr(z500_anom), n_total))

date_index <- tibble(
  layer = seq_len(n_total),
  date  = all_dates,
  year  = year(all_dates),
  month = month(all_dates)
) %>%
  left_join(spei_df %>% select(date, spei6, drought_class, is_drought, season),
            by = "date")

cat(sprintf("  Date index built: %d entries.\n", nrow(date_index)))

# ── 1d. Extract area-averaged time series ─────────────────────────────────────
cat("  Extracting area-averaged time series...\n")

ridge_ext  <- ext(RIDGE_LON_MIN, RIDGE_LON_MAX, RIDGE_LAT_MIN, RIDGE_LAT_MAX)
sst_ext    <- ext(SST_MEAN_LON_MIN, SST_MEAN_LON_MAX, SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX)

# Crop and compute spatial mean for each month
z500_ridge_ts <- global(crop(z500_anom, ridge_ext), "mean", na.rm = TRUE)[, 1]
slp_nwbc_ts   <- global(crop(slp_anom,  ridge_ext), "mean", na.rm = TRUE)[, 1]
sst_nepac_ts  <- global(crop(sst_anom,  sst_ext),   "mean", na.rm = TRUE)[, 1]

ts_df <- date_index %>%
  mutate(
    z500_ridge = z500_ridge_ts,
    slp_nwbc   = slp_nwbc_ts,
    sst_nepac  = sst_nepac_ts
  )

cat("  Area-averaged time series complete.\n\n")

# Load map layers (suppress output)
map_layers <- suppressMessages(suppressWarnings(get_map_layers()))

# ============================================================================
#  SECTION 2: 500 hPa GEOPOTENTIAL HEIGHT ANOMALIES
# ============================================================================
cat("============================================================\n")
cat(" SECTION 2: 500 hPa Geopotential Height\n")
cat("============================================================\n")

# ── 2a. Monthly mean anomaly composites: ALL drought months by calendar month ─
# For each calendar month (Jan-Dec), average Z500 anomaly over all drought months
cat("  Fig 01: Monthly Z500 anomaly composites (drought months)... ")

plot_data_list <- vector("list", 12)
for (m in 1:12) {
  idx <- which(date_index$month == m & date_index$is_drought %in% TRUE)
  if (length(idx) == 0) next
  mean_layer <- app(subset(z500_anom, idx), mean, na.rm = TRUE)
  df <- as.data.frame(mean_layer, xy = TRUE)
  names(df)[3] <- "value"
  df$month_label <- MONTH_LABELS[m]
  df$n_months    <- length(idx)
  plot_data_list[[m]] <- df
}
z500_monthly_comp <- bind_rows(plot_data_list) %>%
  mutate(month_label = factor(month_label, levels = MONTH_LABELS))

clim_range <- quantile(abs(z500_monthly_comp$value), 0.98, na.rm = TRUE)
clim_lim   <- ceiling(clim_range / 10) * 10

p_z500_monthly <- ggplot(z500_monthly_comp, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = NA, color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  geom_sf(data = map_layers$lakes, fill = "white", color = "grey60",
          linewidth = 0.2, inherit.aes = FALSE) +
  geom_point(aes(x = NECHAKO_LON, y = NECHAKO_LAT),
             color = "gold", shape = 18, size = 2.5,
             inherit.aes = FALSE) +
  scale_fill_distiller(palette = ANOMALY_PALETTE,
                       limits  = c(-clim_lim, clim_lim),
                       oob     = squish,
                       name    = "Z500 Anom\n(m)") +
  coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
  facet_wrap(~month_label, ncol = 4) +
  labs(
    title    = "500 hPa Geopotential Height Anomaly — Drought Months by Calendar Month",
    subtitle = sprintf("Composite of all months with SPEI-6 ≤ %.1f  |  Reference: %d–%d  |  ◆ = Nechako Watershed",
                       DROUGHT_THRESHOLD, CLIM_START, CLIM_END),
    caption  = paste0("n shown = number of drought months per panel (varies by month)\n",
                      "Positive (red) = ridge/blocking above normal; Negative (blue) = trough"),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 12),
    panel.spacing    = unit(0.3, "lines")
  )

ggsave(file.path(OUT_DIR, "Figures", "01_z500_monthly_anomaly_composites.png"),
       p_z500_monthly, width = FIGURE_WIDTH_WIDE, height = 12, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── 2b. Composite difference: drought minus non-drought ───────────────────────
cat("  Fig 02: Z500 drought vs. non-drought composite difference... ")

drought_idx     <- which(date_index$is_drought %in% TRUE)
nondrought_idx  <- which(date_index$is_drought %in% FALSE)

z500_drought_mean    <- app(subset(z500_anom, drought_idx),    mean, na.rm = TRUE)
z500_nondrought_mean <- app(subset(z500_anom, nondrought_idx), mean, na.rm = TRUE)
z500_diff            <- z500_drought_mean - z500_nondrought_mean

df_drought    <- as.data.frame(z500_drought_mean,    xy = TRUE) %>% mutate(panel = "Drought months")
df_nondrought <- as.data.frame(z500_nondrought_mean, xy = TRUE) %>% mutate(panel = "Non-drought months")
df_diff       <- as.data.frame(z500_diff,            xy = TRUE) %>% mutate(panel = "Drought minus Non-drought")
names(df_drought)[3]    <- "value"
names(df_nondrought)[3] <- "value"
names(df_diff)[3]       <- "value"

z500_comp_df <- bind_rows(df_drought, df_nondrought, df_diff) %>%
  mutate(panel = factor(panel, levels = c("Drought months",
                                          "Non-drought months",
                                          "Drought minus Non-drought")))

diff_lim <- quantile(abs(df_diff$value), 0.98, na.rm = TRUE)
diff_lim <- ceiling(diff_lim / 5) * 5

p_z500_comp <- ggplot(z500_comp_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = NA, color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(aes(x = NECHAKO_LON, y = NECHAKO_LAT),
             color = "gold", shape = 18, size = 3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = ANOMALY_PALETTE,
                       limits  = c(-diff_lim, diff_lim),
                       oob     = squish,
                       name    = "Z500 Anom\n(m)") +
  coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
  facet_wrap(~panel, ncol = 3) +
  labs(
    title    = "500 hPa Geopotential Height: Drought vs. Non-drought Composites",
    subtitle = sprintf("Drought: SPEI-6 ≤ %.1f (n=%d months)  |  Non-drought (n=%d months)  |  Ref: %d–%d",
                       DROUGHT_THRESHOLD, length(drought_idx),
                       length(nondrought_idx), CLIM_START, CLIM_END),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold", size = 11),
        plot.title       = element_text(face = "bold", size = 13))

ggsave(file.path(OUT_DIR, "Figures", "02_z500_drought_vs_nondrought.png"),
       p_z500_comp, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── 2c. Z500 composites stratified by drought severity ───────────────────────
cat("  Fig 03: Z500 composites by drought severity quartile... ")

severity_levels <- c("Severe/Extreme (SPEI ≤ -1.5)",
                     "Moderate (-1.5 < SPEI ≤ -0.8)",
                     "Non-drought (SPEI > -0.8)")

sev_idx  <- which(date_index$spei6 <= SEVERE_THRESHOLD)
mod_idx  <- which(date_index$spei6 >  SEVERE_THRESHOLD &
                    date_index$spei6 <= DROUGHT_THRESHOLD)
none_idx <- which(date_index$spei6 >  DROUGHT_THRESHOLD)

make_severity_df <- function(idx, label) {
  if (length(idx) == 0) return(NULL)
  m <- app(subset(z500_anom, idx), mean, na.rm = TRUE)
  df <- as.data.frame(m, xy = TRUE)
  names(df)[3] <- "value"
  df$panel <- label
  df$n     <- length(idx)
  df
}

sev_df <- bind_rows(
  make_severity_df(sev_idx,  severity_levels[1]),
  make_severity_df(mod_idx,  severity_levels[2]),
  make_severity_df(none_idx, severity_levels[3])
) %>% mutate(panel = factor(panel, levels = severity_levels))

sev_lim <- quantile(abs(sev_df$value), 0.98, na.rm = TRUE)
sev_lim <- ceiling(sev_lim / 10) * 10

p_z500_sev <- ggplot(sev_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = NA, color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(aes(x = NECHAKO_LON, y = NECHAKO_LAT),
             color = "gold", shape = 18, size = 3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = ANOMALY_PALETTE,
                       limits  = c(-sev_lim, sev_lim),
                       oob     = squish,
                       name    = "Z500 Anom\n(m)") +
  coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
  facet_wrap(~panel, ncol = 3) +
  labs(
    title    = "500 hPa Geopotential Height Anomaly by Drought Severity",
    subtitle = sprintf("Reference period: %d–%d  |  ◆ = Nechako Watershed",
                       CLIM_START, CLIM_END),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 13))

ggsave(file.path(OUT_DIR, "Figures", "03_z500_drought_severity_composites.png"),
       p_z500_sev, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── 2d. Z500 NW Canada ridge time series vs SPEI-6 ────────────────────────────
cat("  Fig 04: Z500 ridge time series vs SPEI-6... ")

p_z500_ts <- ts_df %>%
  filter(!is.na(spei6), !is.na(z500_ridge)) %>%
  ggplot(aes(x = date)) +
  # Z500 ridge anomaly as bar chart
  geom_col(aes(y = z500_ridge / 5, fill = z500_ridge > 0),
           alpha = 0.5, width = 25) +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                    guide = "none") +
  # SPEI-6 line (secondary axis)
  geom_line(aes(y = spei6 * 10), color = "black", linewidth = 0.6) +
  geom_hline(yintercept = DROUGHT_THRESHOLD * 10, linetype = "dashed",
             color = "darkred", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_y_continuous(
    name     = "Z500 Anomaly (m) — NW Canada ridge region",
    limits   = c(-60, 60),
    sec.axis = sec_axis(~. / 10, name = "SPEI-6 (Nechako)")
  ) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y", expand = c(0.01, 0)) +
  labs(
    title    = "500 hPa Geopotential Height Anomaly (NW Canada Ridge Region) vs. SPEI-6",
    subtitle = sprintf("Ridge region: %d–%d°W, %d–%d°N  |  Dashed line: drought threshold (SPEI-6 = %.1f)",
                       abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                       RIDGE_LAT_MIN, RIDGE_LAT_MAX, DROUGHT_THRESHOLD),
    x = NULL,
    caption  = "Bars: Z500 anomaly (red = above normal / ridge, blue = below normal / trough)\nLine: SPEI-6 (lower = drier)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title     = element_text(face = "bold"),
        axis.title.y   = element_text(color = "#d73027"),
        axis.title.y.right = element_text(color = "black"))

ggsave(file.path(OUT_DIR, "Figures", "04_z500_timeseries_nw_canada_ridge.png"),
       p_z500_ts, width = FIGURE_WIDTH_WIDE, height = 5.5, dpi = FIGURE_DPI)
cat("Saved.\n\n")

# ============================================================================
#  SECTION 3: SEA LEVEL PRESSURE
# ============================================================================
cat("============================================================\n")
cat(" SECTION 3: Sea Level Pressure\n")
cat("============================================================\n")

# ── 3a. Monthly SLP anomaly composites ────────────────────────────────────────
cat("  Fig 05: Monthly SLP anomaly composites (drought months)... ")

slp_monthly_list <- vector("list", 12)
for (m in 1:12) {
  idx <- which(date_index$month == m & date_index$is_drought %in% TRUE)
  if (length(idx) == 0) next
  mean_layer <- app(subset(slp_anom, idx), mean, na.rm = TRUE)
  df <- as.data.frame(mean_layer, xy = TRUE)
  names(df)[3] <- "value"
  df$month_label <- MONTH_LABELS[m]
  slp_monthly_list[[m]] <- df
}
slp_monthly_comp <- bind_rows(slp_monthly_list) %>%
  mutate(month_label = factor(month_label, levels = MONTH_LABELS))

slp_clim_lim <- quantile(abs(slp_monthly_comp$value), 0.98, na.rm = TRUE)
slp_clim_lim <- ceiling(slp_clim_lim / 0.5) * 0.5

p_slp_monthly <- ggplot(slp_monthly_comp, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = NA, color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(aes(x = NECHAKO_LON, y = NECHAKO_LAT),
             color = "gold", shape = 18, size = 2.5, inherit.aes = FALSE) +
  scale_fill_distiller(palette = ANOMALY_PALETTE,
                       limits  = c(-slp_clim_lim, slp_clim_lim),
                       oob     = squish,
                       name    = "SLP Anom\n(hPa)") +
  coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
  facet_wrap(~month_label, ncol = 4) +
  labs(
    title    = "Sea Level Pressure Anomaly — Drought Months by Calendar Month",
    subtitle = sprintf("Composite of all months with SPEI-6 ≤ %.1f  |  Reference: %d–%d",
                       DROUGHT_THRESHOLD, CLIM_START, CLIM_END),
    caption  = "Positive (red) = high pressure anomaly (blocking); Negative (blue) = low pressure / trough",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        legend.position  = "right",
        plot.title       = element_text(face = "bold", size = 12),
        panel.spacing    = unit(0.3, "lines"))

ggsave(file.path(OUT_DIR, "Figures", "05_slp_monthly_anomaly_composites.png"),
       p_slp_monthly, width = FIGURE_WIDTH_WIDE, height = 12, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── 3b. SLP drought vs. non-drought ───────────────────────────────────────────
cat("  Fig 06: SLP drought vs. non-drought composite... ")

slp_dr_mean   <- app(subset(slp_anom, drought_idx),    mean, na.rm = TRUE)
slp_nodr_mean <- app(subset(slp_anom, nondrought_idx), mean, na.rm = TRUE)
slp_diff      <- slp_dr_mean - slp_nodr_mean

df_slp <- bind_rows(
  as.data.frame(slp_dr_mean,   xy = TRUE) %>%
    setNames(c("x","y","value")) %>% mutate(panel = "Drought months"),
  as.data.frame(slp_nodr_mean, xy = TRUE) %>%
    setNames(c("x","y","value")) %>% mutate(panel = "Non-drought months"),
  as.data.frame(slp_diff,      xy = TRUE) %>%
    setNames(c("x","y","value")) %>% mutate(panel = "Drought minus Non-drought")
) %>%
  mutate(panel = factor(panel, levels = c("Drought months",
                                          "Non-drought months",
                                          "Drought minus Non-drought")))

slp_diff_lim <- quantile(abs(as.data.frame(slp_diff)[, 3]), 0.98, na.rm = TRUE)
slp_diff_lim <- ceiling(slp_diff_lim / 0.5) * 0.5

p_slp_comp <- ggplot(df_slp, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = NA, color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(aes(x = NECHAKO_LON, y = NECHAKO_LAT),
             color = "gold", shape = 18, size = 3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = ANOMALY_PALETTE,
                       limits  = c(-slp_diff_lim, slp_diff_lim),
                       oob     = squish,
                       name    = "SLP Anom\n(hPa)") +
  coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
  facet_wrap(~panel, ncol = 3) +
  labs(
    title    = "Sea Level Pressure: Drought vs. Non-drought Composites",
    subtitle = sprintf("Drought: SPEI-6 ≤ %.1f (n=%d)  |  Non-drought (n=%d)  |  Ref: %d–%d",
                       DROUGHT_THRESHOLD, length(drought_idx),
                       length(nondrought_idx), CLIM_START, CLIM_END),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold", size = 11),
        plot.title       = element_text(face = "bold", size = 13))

ggsave(file.path(OUT_DIR, "Figures", "06_slp_drought_vs_nondrought.png"),
       p_slp_comp, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── 3c. SLP severity composites ───────────────────────────────────────────────
cat("  Fig 07: SLP composites by drought severity... ")

make_sev_slp <- function(idx, label) {
  if (length(idx) == 0) return(NULL)
  m  <- app(subset(slp_anom, idx), mean, na.rm = TRUE)
  df <- as.data.frame(m, xy = TRUE)
  names(df)[3] <- "value"
  df$panel <- label
  df
}

slp_sev_df <- bind_rows(
  make_sev_slp(sev_idx,  severity_levels[1]),
  make_sev_slp(mod_idx,  severity_levels[2]),
  make_sev_slp(none_idx, severity_levels[3])
) %>% mutate(panel = factor(panel, levels = severity_levels))

slp_sev_lim <- quantile(abs(slp_sev_df$value), 0.98, na.rm = TRUE)
slp_sev_lim <- ceiling(slp_sev_lim / 0.5) * 0.5

p_slp_sev <- ggplot(slp_sev_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = NA, color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  geom_point(aes(x = NECHAKO_LON, y = NECHAKO_LAT),
             color = "gold", shape = 18, size = 3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = ANOMALY_PALETTE,
                       limits  = c(-slp_sev_lim, slp_sev_lim),
                       oob     = squish,
                       name    = "SLP Anom\n(hPa)") +
  coord_sf(xlim = c(-180, -100), ylim = c(40, 70), expand = FALSE) +
  facet_wrap(~panel, ncol = 3) +
  labs(title    = "SLP Anomaly by Drought Severity",
       subtitle = sprintf("Reference: %d–%d  |  ◆ = Nechako", CLIM_START, CLIM_END),
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 13))

ggsave(file.path(OUT_DIR, "Figures", "07_slp_drought_severity_composites.png"),
       p_slp_sev, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── 3d. SLP time series ────────────────────────────────────────────────────────
cat("  Fig 08: SLP NW BC time series vs SPEI-6... ")

p_slp_ts <- ts_df %>%
  filter(!is.na(spei6), !is.na(slp_nwbc)) %>%
  ggplot(aes(x = date)) +
  geom_col(aes(y = slp_nwbc * 5, fill = slp_nwbc > 0),
           alpha = 0.5, width = 25) +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                    guide = "none") +
  geom_line(aes(y = spei6 * 2), color = "black", linewidth = 0.6) +
  geom_hline(yintercept = DROUGHT_THRESHOLD * 2, linetype = "dashed",
             color = "darkred", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  scale_y_continuous(
    name     = "SLP Anomaly (hPa × 5) — NW BC region",
    sec.axis = sec_axis(~. / 2, name = "SPEI-6 (Nechako)")
  ) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y", expand = c(0.01, 0)) +
  labs(
    title   = "Sea Level Pressure Anomaly (NW BC Region) vs. SPEI-6",
    subtitle = "Bars: SLP anomaly (red = anomalous high pressure / blocking, blue = low pressure)\nLine: SPEI-6",
    x = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Figures", "08_slp_timeseries.png"),
       p_slp_ts, width = FIGURE_WIDTH_WIDE, height = 5.5, dpi = FIGURE_DPI)
cat("Saved.\n\n")

# ============================================================================
#  SECTION 4: SST NE PACIFIC
# ============================================================================
cat("============================================================\n")
cat(" SECTION 4: SST — NE Pacific\n")
cat("============================================================\n")

# ── 4a. Annual SST anomaly maps — recent years ────────────────────────────────
cat("  Fig 09: Annual mean SST anomaly maps for recent years... ")

recent_years <- RECENT_START:min(RECENT_END, END_YEAR)
sst_annual_list <- vector("list", length(recent_years))

for (i in seq_along(recent_years)) {
  yr  <- recent_years[i]
  idx <- which(date_index$year == yr)
  if (length(idx) == 0) next
  mean_sst <- app(subset(sst_anom, idx), mean, na.rm = TRUE)
  df <- as.data.frame(mean_sst, xy = TRUE)
  names(df)[3] <- "value"
  df$year <- yr
  sst_annual_list[[i]] <- df
}

sst_annual_df <- bind_rows(sst_annual_list) %>%
  mutate(year = factor(year))

sst_lim <- max(abs(quantile(sst_annual_df$value, c(0.02, 0.98), na.rm = TRUE)))
sst_lim <- ceiling(sst_lim * 2) / 2

p_sst_annual <- ggplot(sst_annual_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = "grey85", color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = SST_PALETTE,
                       limits  = c(-sst_lim, sst_lim),
                       oob     = squish,
                       name    = "SST Anom\n(°C)") +
  coord_sf(xlim = c(-180, -110), ylim = c(20, 65), expand = FALSE) +
  facet_wrap(~year, ncol = 4) +
  labs(
    title    = sprintf("NE Pacific SST Anomaly — Annual Mean, %d–%d",
                       RECENT_START, min(RECENT_END, END_YEAR)),
    subtitle = sprintf("Reference: %d–%d WMO climatology  |  Positive (red) = warmer than normal",
                       CLIM_START, CLIM_END),
    caption  = "Captures PDO, Marine Heatwave / Blob signal and California Current variability",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 12))

ggsave(file.path(OUT_DIR, "Figures", "09_sst_annual_anomaly_recent.png"),
       p_sst_annual,
       width  = FIGURE_WIDTH_WIDE,
       height = ceiling(length(recent_years) / 4) * 4 + 2,
       dpi    = FIGURE_DPI)
cat("Saved.\n")

# ── 4b. SST composites during drought months ──────────────────────────────────
cat("  Fig 10: SST composites during drought months (by calendar month)... ")

sst_monthly_list <- vector("list", 12)
for (m in 1:12) {
  idx <- which(date_index$month == m & date_index$is_drought %in% TRUE)
  if (length(idx) == 0) next
  mean_sst <- app(subset(sst_anom, idx), mean, na.rm = TRUE)
  df <- as.data.frame(mean_sst, xy = TRUE)
  names(df)[3] <- "value"
  df$month_label <- MONTH_LABELS[m]
  sst_monthly_list[[m]] <- df
}

sst_monthly_comp <- bind_rows(sst_monthly_list) %>%
  mutate(month_label = factor(month_label, levels = MONTH_LABELS))

sst_month_lim <- quantile(abs(sst_monthly_comp$value), 0.98, na.rm = TRUE)
sst_month_lim <- ceiling(sst_month_lim * 2) / 2

p_sst_monthly <- ggplot(sst_monthly_comp, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  geom_sf(data = map_layers$borders, fill = "grey85", color = "black",
          linewidth = 0.3, inherit.aes = FALSE) +
  scale_fill_distiller(palette = SST_PALETTE,
                       limits  = c(-sst_month_lim, sst_month_lim),
                       oob     = squish,
                       name    = "SST Anom\n(°C)") +
  coord_sf(xlim = c(-180, -110), ylim = c(20, 65), expand = FALSE) +
  facet_wrap(~month_label, ncol = 4) +
  labs(
    title    = "NE Pacific SST Anomaly — Drought Months by Calendar Month",
    subtitle = sprintf("Composite of all months with SPEI-6 ≤ %.1f  |  Reference: %d–%d",
                       DROUGHT_THRESHOLD, CLIM_START, CLIM_END),
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold", size = 12),
        panel.spacing    = unit(0.3, "lines"))

ggsave(file.path(OUT_DIR, "Figures", "10_sst_monthly_composites_drought.png"),
       p_sst_monthly, width = FIGURE_WIDTH_WIDE, height = 12, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── 4c. SST NE Pacific time series: long-record + recent highlight ────────────
cat("  Fig 11: SST NE Pacific time series... ")

# 12-month rolling mean for smoothed signal
ts_df_sst <- ts_df %>%
  filter(!is.na(sst_nepac)) %>%
  arrange(date) %>%
  mutate(
    sst_12mo = zoo::rollmean(sst_nepac, k = 12, fill = NA, align = "right")
  )

p_sst_ts <- ggplot(ts_df_sst, aes(x = date, y = sst_nepac)) +
  # shade recent period
  annotate("rect",
           xmin = as.Date(sprintf("%d-01-01", RECENT_START)),
           xmax = as.Date(sprintf("%d-12-31", min(RECENT_END, END_YEAR))),
           ymin = -Inf, ymax = Inf,
           fill = "lightyellow", alpha = 0.5) +
  geom_col(aes(fill = sst_nepac > 0), alpha = 0.4, width = 25) +
  scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
                    guide = "none") +
  geom_line(aes(y = sst_12mo), color = "black", linewidth = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_hline(yintercept =  0.5, linetype = "dotted", color = "darkred",  linewidth = 0.5) +
  geom_hline(yintercept = -0.5, linetype = "dotted", color = "steelblue", linewidth = 0.5) +
  annotate("text",
           x = as.Date(sprintf("%d-06-01", RECENT_START + 1)),
           y = max(ts_df_sst$sst_nepac, na.rm = TRUE) * 0.9,
           label = sprintf("Recent period\n(%d–%d)", RECENT_START, min(RECENT_END, END_YEAR)),
           size = 3, hjust = 0, color = "grey40") +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y", expand = c(0.01, 0)) +
  labs(
    title    = "NE Pacific SST Anomaly — Monthly Time Series",
    subtitle = sprintf("Area average: %d–%d°W, %d–%d°N  |  Black line: 12-month rolling mean  |  Ref: %d–%d",
                       abs(SST_MEAN_LON_MIN), abs(SST_MEAN_LON_MAX),
                       SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX,
                       CLIM_START, CLIM_END),
    x = NULL, y = "SST Anomaly (°C)",
    caption  = "Dotted lines: ±0.5°C thresholds (approximate warm/cool PDO phase boundary)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Figures", "11_sst_timeseries_ne_pacific.png"),
       p_sst_ts, width = FIGURE_WIDTH_WIDE, height = 5.5, dpi = FIGURE_DPI)
cat("Saved.\n\n")

# ============================================================================
#  SECTION 5: JOINT ANALYSIS — SCATTERPLOTS AND CORRELATION
# ============================================================================
cat("============================================================\n")
cat(" SECTION 5: Joint Analysis — Area-averaged fields vs SPEI-6\n")
cat("============================================================\n")

cat("  Fig 12: Scatter plots of Z500/SLP/SST area-averages vs SPEI-6... ")

scatter_df <- ts_df %>%
  filter(!is.na(spei6), !is.na(z500_ridge), !is.na(slp_nwbc), !is.na(sst_nepac)) %>%
  mutate(season = factor(season, levels = c("DJF", "MAM", "JJA", "SON")))

make_scatter <- function(df, xvar, xlabel, title_text) {
  ggplot(df, aes(x = .data[[xvar]], y = spei6, color = season)) +
    geom_hline(yintercept = DROUGHT_THRESHOLD, linetype = "dashed",
               color = "darkred", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                linewidth = 0.8, fill = "grey70", alpha = 0.3) +
    scale_color_brewer(palette = "Set1", name = "Season") +
    labs(
      title    = title_text,
      x        = xlabel,
      y        = "SPEI-6 (Nechako)",
      subtitle = sprintf("r = %.3f  (p < 0.001 if |r| > %.3f, n = %d)",
                         cor(df[[xvar]], df$spei6, use = "complete.obs"),
                         qnorm(0.975) / sqrt(nrow(df) - 2) * sqrt(1),
                         nrow(df))
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          legend.position = "bottom")
}

p_sc1 <- make_scatter(scatter_df, "z500_ridge",
                      sprintf("Z500 Anomaly (m)\nNW Canada ridge (%d–%d°W, %d–%d°N)",
                              abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                              RIDGE_LAT_MIN, RIDGE_LAT_MAX),
                      "Z500 Ridge Anomaly vs. SPEI-6")

p_sc2 <- make_scatter(scatter_df, "slp_nwbc",
                      sprintf("SLP Anomaly (hPa)\nNW BC (%d–%d°W, %d–%d°N)",
                              abs(RIDGE_LON_MIN), abs(RIDGE_LON_MAX),
                              RIDGE_LAT_MIN, RIDGE_LAT_MAX),
                      "SLP Anomaly vs. SPEI-6")

p_sc3 <- make_scatter(scatter_df, "sst_nepac",
                      sprintf("SST Anomaly (°C)\nNE Pacific (%d–%d°W, %d–%d°N)",
                              abs(SST_MEAN_LON_MIN), abs(SST_MEAN_LON_MAX),
                              SST_MEAN_LAT_MIN, SST_MEAN_LAT_MAX),
                      "NE Pacific SST Anomaly vs. SPEI-6")

p_scatter <- (p_sc1 | p_sc2 | p_sc3) +
  plot_annotation(
    title   = "Area-averaged Atmospheric / Ocean Fields vs. Nechako SPEI-6",
    caption = sprintf("Dashed red line = drought threshold (SPEI-6 = %.1f). Shading = 95%% confidence band on regression.",
                      DROUGHT_THRESHOLD),
    theme   = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, "Figures", "12_joint_z500_slp_sst_drought_scatter.png"),
       p_scatter, width = FIGURE_WIDTH_WIDE, height = 6, dpi = FIGURE_DPI)
cat("Saved.\n")

# ── 5b. Summary composite panel ───────────────────────────────────────────────
cat("  Fig 13: Summary composite panel (drought-minus-nondrought for all 3 vars)... ")

# Crop SST diff to ATM domain for visual consistency — use full SST domain instead
sst_dr_mean   <- app(subset(sst_anom, drought_idx),    mean, na.rm = TRUE)
sst_nodr_mean <- app(subset(sst_anom, nondrought_idx), mean, na.rm = TRUE)
sst_diff_r    <- sst_dr_mean - sst_nodr_mean

lim_z  <- ceiling(quantile(abs(as.data.frame(z500_diff)[,3]), 0.98, na.rm = TRUE) / 10) * 10
lim_sl <- ceiling(quantile(abs(as.data.frame(slp_diff) [,3]), 0.98, na.rm = TRUE) / 0.5) * 0.5
lim_ss <- ceiling(quantile(abs(as.data.frame(sst_diff_r)[,3]), 0.98, na.rm = TRUE) * 2) / 2

make_composite_map <- function(r, lim, var_label, unit, xlim, ylim) {
  df <- as.data.frame(r, xy = TRUE); names(df)[3] <- "value"
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    geom_sf(data = map_layers$borders, fill = NA, color = "black",
            linewidth = 0.25, inherit.aes = FALSE) +
    geom_point(aes(x = NECHAKO_LON, y = NECHAKO_LAT),
               color = "gold", shape = 18, size = 3, inherit.aes = FALSE) +
    scale_fill_distiller(palette = "RdBu",
                         limits  = c(-lim, lim),
                         oob     = squish,
                         name    = unit) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = var_label, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
          legend.position = "right")
}

p_sum1 <- make_composite_map(z500_diff, lim_z,
                             "Z500 Anomaly Difference (m)\nDrought minus Non-drought",
                             "Z500\n(m)", c(-180, -100), c(40, 70))

p_sum2 <- make_composite_map(slp_diff, lim_sl,
                             "SLP Anomaly Difference (hPa)\nDrought minus Non-drought",
                             "SLP\n(hPa)", c(-180, -100), c(40, 70))

p_sum3 <- make_composite_map(sst_diff_r, lim_ss,
                             "SST Anomaly Difference (°C)\nDrought minus Non-drought",
                             "SST\n(°C)", c(-180, -110), c(20, 65))

p_summary_panel <- (p_sum1 / p_sum2 / p_sum3) +
  plot_annotation(
    title    = "Composite Atmospheric & Ocean Anomalies Associated with Nechako Drought",
    subtitle = sprintf("Drought: SPEI-6 ≤ %.1f (n=%d months)  minus  Non-drought (n=%d months)  |  ◆ = Nechako Watershed",
                       DROUGHT_THRESHOLD, length(drought_idx), length(nondrought_idx)),
    caption  = sprintf("Reference period: %d–%d (WMO standard). ERA5 monthly means.", CLIM_START, CLIM_END),
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(size = 10))
  )

ggsave(file.path(OUT_DIR, "Figures", "13_composite_summary_panel.png"),
       p_summary_panel,
       width  = FIGURE_WIDTH_STD,
       height = 16,
       dpi    = FIGURE_DPI)
cat("Saved.\n\n")

# ============================================================================
#  SECTION 6: TABLES
# ============================================================================
cat("============================================================\n")
cat(" SECTION 6: Summary Tables\n")
cat("============================================================\n")

# ── T1: Drought month summary — Z500/SLP/SST stats per SPEI-6 class ──────────
cat("  Table T1: Drought month summary statistics... ")

t1 <- ts_df %>%
  filter(!is.na(spei6), !is.na(z500_ridge)) %>%
  mutate(drought_class = case_when(
    spei6 <= SEVERE_THRESHOLD  ~ sprintf("Severe/Extreme (SPEI ≤ %.1f)", SEVERE_THRESHOLD),
    spei6 <= DROUGHT_THRESHOLD ~ sprintf("Moderate (%.1f < SPEI ≤ %.1f)", SEVERE_THRESHOLD, DROUGHT_THRESHOLD),
    TRUE                       ~ "Non-drought"
  )) %>%
  group_by(drought_class) %>%
  summarise(
    n_months           = n(),
    SPEI6_mean         = round(mean(spei6,     na.rm = TRUE), 3),
    SPEI6_sd           = round(sd(spei6,       na.rm = TRUE), 3),
    Z500_ridge_mean_m  = round(mean(z500_ridge, na.rm = TRUE), 1),
    Z500_ridge_sd_m    = round(sd(z500_ridge,   na.rm = TRUE), 1),
    SLP_nwbc_mean_hPa  = round(mean(slp_nwbc,   na.rm = TRUE), 2),
    SLP_nwbc_sd_hPa    = round(sd(slp_nwbc,     na.rm = TRUE), 2),
    SST_nepac_mean_C   = round(mean(sst_nepac,   na.rm = TRUE), 3),
    SST_nepac_sd_C     = round(sd(sst_nepac,     na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(SPEI6_mean)

write_csv(t1, file.path(OUT_DIR, "Tables", "T1_drought_month_summary.csv"))
cat("Saved.\n")
print(t1)

# ── T2: Seasonal anomaly breakdown ────────────────────────────────────────────
cat("\n  Table T2: Seasonal anomaly statistics... ")

t2 <- ts_df %>%
  filter(!is.na(spei6), is_drought %in% TRUE,
         !is.na(z500_ridge), !is.na(slp_nwbc), !is.na(sst_nepac)) %>%
  group_by(season) %>%
  summarise(
    n_drought_months   = n(),
    SPEI6_mean         = round(mean(spei6,      na.rm = TRUE), 3),
    Z500_ridge_mean    = round(mean(z500_ridge,  na.rm = TRUE), 1),
    Z500_ridge_sd      = round(sd(z500_ridge,    na.rm = TRUE), 1),
    SLP_mean_hPa       = round(mean(slp_nwbc,    na.rm = TRUE), 2),
    SLP_sd_hPa         = round(sd(slp_nwbc,      na.rm = TRUE), 2),
    SST_nepac_mean_C   = round(mean(sst_nepac,    na.rm = TRUE), 3),
    SST_nepac_sd_C     = round(sd(sst_nepac,      na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(match(season, c("DJF", "MAM", "JJA", "SON")))

write_csv(t2, file.path(OUT_DIR, "Tables", "T2_seasonal_anomaly_statistics.csv"))
cat("Saved.\n")
print(t2)

# ── T3: Correlation matrix — area-averaged fields vs SPEI-6 ──────────────────
cat("\n  Table T3: Correlation matrix... ")

cor_vars <- c("spei6", "z500_ridge", "slp_nwbc", "sst_nepac")
cor_labels <- c("SPEI-6", "Z500 Ridge (m)", "SLP NW-BC (hPa)", "SST NE-Pac (°C)")

cor_mat_data <- ts_df %>%
  filter(!is.na(spei6), !is.na(z500_ridge), !is.na(slp_nwbc), !is.na(sst_nepac)) %>%
  select(all_of(cor_vars))

n_obs  <- nrow(cor_mat_data)
cor_r  <- cor(cor_mat_data, use = "complete.obs")
t_stat <- cor_r * sqrt((n_obs - 2) / (1 - cor_r^2))
p_mat  <- 2 * pt(-abs(t_stat), df = n_obs - 2)

# Build tidy output table
cor_out <- as.data.frame(round(cor_r, 3))
cor_out <- cbind(Variable = cor_labels, cor_out)
names(cor_out)[-1] <- cor_labels

write_csv(cor_out, file.path(OUT_DIR, "Tables", "T3_correlation_matrix.csv"))
cat("Saved.\n")
print(cor_out)

# Also save p-values
pval_out <- as.data.frame(round(p_mat, 4))
pval_out <- cbind(Variable = cor_labels, pval_out)
names(pval_out)[-1] <- cor_labels
write_csv(pval_out, file.path(OUT_DIR, "Tables", "T3_pvalue_matrix.csv"))

# ── T4: Annual SST anomaly summary — recent vs historical ─────────────────────
cat("\n  Table T4: Annual SST anomaly — recent vs historical... ")

t4 <- ts_df %>%
  filter(!is.na(sst_nepac)) %>%
  group_by(year) %>%
  summarise(
    SST_annual_mean_C = round(mean(sst_nepac, na.rm = TRUE), 3),
    SST_annual_max_C  = round(max(sst_nepac,  na.rm = TRUE), 3),
    SST_annual_min_C  = round(min(sst_nepac,  na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(
    period = ifelse(year >= RECENT_START, sprintf("Recent (%d–)", RECENT_START),
                    sprintf("Historical (%d–%d)", START_YEAR, RECENT_START - 1))
  )

write_csv(t4, file.path(OUT_DIR, "Tables", "T4_annual_sst_anomaly_summary.csv"))
cat("Saved.\n")

# Print recent years
cat("  Recent years:\n")
print(tail(t4, length(recent_years)))

# ============================================================================
#  DONE
# ============================================================================
cat("\n============================================================\n")
cat(" COMPLETE. All outputs written to:\n")
cat(sprintf("   %s\n", OUT_DIR))
cat("------------------------------------------------------------\n")
cat(" Figures:\n")
fig_files <- list.files(file.path(OUT_DIR, "Figures"), full.names = FALSE)
for (f in fig_files) cat(sprintf("   %s\n", f))
cat(" Tables:\n")
tbl_files <- list.files(file.path(OUT_DIR, "Tables"), full.names = FALSE)
for (f in tbl_files) cat(sprintf("   %s\n", f))
cat("============================================================\n")