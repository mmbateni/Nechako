####################################################################################
# w2_basin_timeseries.R  ·  BASIN-AVERAGED TIME SERIES FOR ALL INDICES
# ─────────────────────────────────────────────────────────────────────────────────
# SCOPE  (authoritative home for all basin-averaged outputs):
#   PART 1 – Load basin-averaged CSVs produced by w1_trend_test.R for
#             SPI, SPEI, SWEI, and optionally MSPI / MSPEI.
#   PART 2 – Write per-index flat CSVs to basin_averaged_timeseries/
#             (re-exports the authoritative *_basin_averaged_by_month.csv files
#             into a plain date/value layout for downstream use).
#   PART 3 – Combined CSV (basin_mean_all_indices.csv — all indices long format).
#   PART 4 – Time series plots (one publication-quality PNG per index × scale).
#   PART 4d– Seasonality bar chart: mean monthly P vs PET(PM) (12 values each),
#             1950–2025 average; climatic water balance panel.
#             Requires ERA5Land_Nechako_PET_monthly_summary.csv from
#             2preq_PET_ERALand.R.  Output: FigS_Seasonality_P_PET.pdf/.png.
#   PART 5 – Excel summary: descriptive statistics + drought event counts
#             (2020-2025) + pairwise correlation matrix (Drought_Summary.xlsx).
#
# PREREQUISITE: run w*_trend_test.R first.
# NOTE: w4_trends_visualization.R handles ONLY spatial figures (Parts 1-2 of w4).
#       All basin-averaged time-series outputs live here.
#
# METHODOLOGY NOTE:
#   Step 1 — SPI/SPEI/SWEI are computed per-pixel by 1SPI_ERALand.R /
#             3SPEI_ERALand.R; results live in seasonal NetCDF files.
#   Step 2 — w1_trend_test.R area-weights per-pixel values with
#             terra::cellSize() in BC Albers (EPSG:3005) and writes
#             {index}_{scale:02d}_basin_averaged_by_month.csv.
#   Step 3 — This script reads those authoritative CSVs.
#             Do NOT recompute averages here.
####################################################################################
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "ggplot2", "lubridate", "dplyr",
                      "data.table", "patchwork", "openxlsx", "zoo"))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

## ── Output directories ──────────────────────────────────────────────────────────
BASIN_PLOT_DIR <- file.path(WD_PATH, "temporal_drought", "basin_averaged_plots")
BASIN_TS_DIR   <- file.path(WD_PATH, "temporal_drought", "basin_averaged_timeseries")
for (d in c(BASIN_PLOT_DIR, BASIN_TS_DIR, CACHE_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ── Basin extent (informational only) ───────────────────────────────────────────
basin_shp  <- terra::vect(BASIN_SHP)
basin_proj <- terra::project(basin_shp, EQUAL_AREA_CRS)
cat(sprintf("✓ Basin: %.1f km²\n", terra::expanse(basin_proj, unit = "km")))

## ── MSPI / MSPEI single scale (mirror of w4 constant) ───────────────────────────
MSPI_MSPEI_SCALE <- 1L

################################################################################
# HELPERS
################################################################################

## Load the authoritative w3 basin-average CSV and return a data.frame
## (Date, Index, Mean_Value, Drought_Fraction) ready for plotting.
load_basin_index <- function(index_type, scale) {
  label <- sprintf("%s%d", toupper(index_type), scale)
  df    <- tryCatch(
    load_basin_avg_csv(index_type, scale),
    error = function(e) { cat(sprintf("  ❌ %s: %s\n", label, e$message)); NULL })
  if (is.null(df) || !nrow(df)) {
    cat(sprintf("  ⚠ No data for %s\n", label)); return(NULL)
  }
  data.frame(
    Date             = df$date,
    Index            = label,
    Mean_Value       = df$value,
    Drought_Fraction = as.numeric(df$value <= DROUGHT_THRESHOLD),
    stringsAsFactors = FALSE)
}

## Load MSPI / MSPEI time series from the Excel file written by their own scripts.
## Returns data.frame(date, value) or NULL if the file is absent.
load_mspi_mspei_ts <- function(index_type) {
  results_dir <- file.path(WD_PATH,
                           if (index_type == "mspi") "mspi_results" else "mspei_results")
  pat  <- sprintf("%s_basin_timeseries.*\\.xlsx$", index_type)
  fls  <- list.files(results_dir, pattern = pat, full.names = TRUE)
  if (!length(fls)) {
    cat(sprintf("  ⚠ No Excel timeseries found for %s in %s\n",
                toupper(index_type), results_dir))
    return(NULL)
  }
  tryCatch({
    df       <- openxlsx::read.xlsx(fls[1], sheet = "Basin_Timeseries")
    df$date  <- as.Date(paste0(df$Date, "-01"))
    df$value <- df$Basin_Mean
    df[, c("date", "value")]
  }, error = function(e) {
    cat(sprintf("  ⚠ Could not read %s Excel: %s\n", toupper(index_type), e$message))
    NULL
  })
}

## Unified loader — returns data.frame(date, value) for any supported index.
## Used by the CSV writer, plot builder, and Excel export.
load_any_ts <- function(idx, sc) {
  if (idx %in% c("mspi", "mspei")) {
    df <- load_mspi_mspei_ts(idx)
    return(if (!is.null(df)) as.data.frame(df) else NULL)
  }
  if (idx == "swei") sc <- SWEI_SCALE
  df <- tryCatch(
    load_basin_avg_csv(idx, sc),
    error = function(e) {
      cat(sprintf("  ⚠ %s-%02d: %s\n", toupper(idx), sc, e$message))
      NULL
    })
  if (!is.null(df)) as.data.frame(df) else NULL
}

################################################################################
# PART 1 – LOAD ALL INDICES
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 1: Loading basin-averaged indices\n")
cat("══════════════════════════════════════\n")

all_data <- list()

cat("\n── SPI ──\n")
for (sc in SPI_SCALES) {
  ts <- load_basin_index("spi", sc)
  if (!is.null(ts)) {
    all_data[[sprintf("SPI%d", sc)]] <- ts
    cat(sprintf("  ✓ SPI%d: %d time steps\n", sc, nrow(ts)))
  }
}

cat("\n── SPEI ──\n")
for (sc in SPEI_SCALES) {
  ts <- load_basin_index("spei", sc)
  if (!is.null(ts)) {
    all_data[[sprintf("SPEI%d", sc)]] <- ts
    cat(sprintf("  ✓ SPEI%d: %d time steps\n", sc, nrow(ts)))
  }
}

cat("\n── SWEI ──\n")
ts <- load_basin_index("swei", SWEI_SCALE)
if (!is.null(ts)) {
  all_data[[sprintf("SWEI%d", SWEI_SCALE)]] <- ts
  cat(sprintf("  ✓ SWEI%d: %d time steps\n", SWEI_SCALE, nrow(ts)))
}

cat("\n── MSPI / MSPEI (optional) ──\n")
for (multi_idx in c("mspi", "mspei")) {
  df_multi <- load_mspi_mspei_ts(multi_idx)
  if (!is.null(df_multi) && nrow(df_multi)) {
    label    <- sprintf("%s%d", toupper(multi_idx), MSPI_MSPEI_SCALE)
    ts_multi <- data.frame(
      Date             = df_multi$date,
      Index            = label,
      Mean_Value       = df_multi$value,
      Drought_Fraction = as.numeric(df_multi$value <= DROUGHT_THRESHOLD),
      stringsAsFactors = FALSE)
    all_data[[label]] <- ts_multi
    cat(sprintf("  ✓ %s: %d time steps\n", label, nrow(ts_multi)))
  }
}

if (!length(all_data))
  stop("No basin-averaged data loaded. Run w1_trend_test.R first.")

################################################################################
# PART 2 – PER-INDEX FLAT CSVs  (date / value layout)
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 2: Writing per-index flat CSVs\n")
cat("══════════════════════════════════════════════\n")

write_basin_ts_csv <- function(idx, sc) {
  df <- load_any_ts(idx, sc)
  if (is.null(df) || !nrow(df)) {
    cat(sprintf("  ⚠ %s-%02d: no data available\n", toupper(idx), sc))
    return(invisible(NULL))
  }
  out <- file.path(BASIN_TS_DIR,
                   sprintf("%s_%02d_basin_average.csv", idx, sc))
  utils::write.csv(df, out, row.names = FALSE)
  cat(sprintf("  ✓ %s-%02d: %d months → %s\n",
              toupper(idx), sc, nrow(df), basename(out)))
  invisible(df)
}

cat("\n  SPI & SPEI...\n")
for (idx in c("spi", "spei")) for (sc in SPI_SCALES) write_basin_ts_csv(idx, sc)

cat("\n  SWEI...\n")
write_basin_ts_csv("swei", SWEI_SCALE)

cat("\n  MSPI & MSPEI (if available)...\n")
for (idx in c("mspi", "mspei")) write_basin_ts_csv(idx, MSPI_MSPEI_SCALE)

################################################################################
# PART 3 – COMBINED LONG-FORMAT CSV
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 3: Combined long-format CSV\n")
cat("══════════════════════════════════════════════\n")

all_basin_data <- dplyr::bind_rows(all_data)
combined_csv   <- file.path(BASIN_PLOT_DIR, "basin_mean_all_indices.csv")
write.csv(all_basin_data, combined_csv, row.names = FALSE)
cat(sprintf("  ✓ %d rows → %s\n", nrow(all_basin_data), basename(combined_csv)))

################################################################################
# PART 4 – TIME SERIES PLOTS  (one PNG per index × scale)
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 4: Time series plots\n")
cat("══════════════════════════════════════════════\n")

## Helper: read back a per-index flat CSV written in Part 2.
load_flat_csv <- function(idx, sc) {
  f <- file.path(BASIN_TS_DIR,
                 sprintf("%s_%02d_basin_average.csv", idx, sc))
  if (!file.exists(f)) stop("CSV not found: ", f)
  df       <- data.table::fread(f)
  df$date  <- as.Date(df$date)
  df$value <- as.numeric(df$value)
  as.data.frame(df)
}

## Publication-quality basin-average time series plot.
## Drought-fill ribbon, categorical band shading, event count subtitle.
make_basin_ts_plot <- function(df, idx, sc) {
  idx_lc    <- tolower(idx)
  idx_color <- index_colours[idx_lc]
  if (is.na(idx_color)) idx_color <- "steelblue"
  
  ## Drought event count for subtitle annotation
  n_ev <- tryCatch(
    nrow(detect_drought_events(data.frame(date = df$date, value = df$value))),
    error = function(e) NA_integer_)
  
  ggplot2::ggplot(df, ggplot2::aes(date, value)) +
    drought_band_layers() +
    ggplot2::geom_hline(yintercept = 0,
                        color = "black", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1, 1),
                        linetype = "dashed", color = "gray50", linewidth = 0.3) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = ifelse(value < DROUGHT_ONSET, value, 0),
        ymax = 0),
      fill = "#c0392b", alpha = 0.35) +
    ggplot2::geom_line(color = idx_color, linewidth = 0.6) +
    shared_ts_theme(14) +
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    ggplot2::coord_cartesian(ylim = c(-3.5, 3.5)) +
    ggplot2::labs(
      title    = sprintf("%s-%02d  Basin-Averaged Time Series (Area-Weighted)",
                         toupper(idx), sc),
      subtitle = sprintf("%d drought events detected  (onset < %.1f)",
                         n_ev, DROUGHT_ONSET),
      x = "Year",
      y = sprintf("%s-%02d Index Value", toupper(idx), sc))
}

cat("\n  SPI & SPEI...\n")
for (idx in c("spi", "spei")) {
  for (sc in SPI_SCALES) {
    tryCatch({
      df  <- load_flat_csv(idx, sc)
      out <- file.path(BASIN_PLOT_DIR,
                       sprintf("%s_%02d_timeseries.png", idx, sc))
      ggplot2::ggsave(out, make_basin_ts_plot(df, idx, sc),
                      width = 12, height = 6, dpi = 150)
      cat(sprintf("  ✓ %s-%02d\n", toupper(idx), sc))
    }, error = function(e)
      cat(sprintf("  ⚠ %s-%02d: %s\n", toupper(idx), sc, e$message)))
  }
}

cat("\n  SWEI...\n")
tryCatch({
  df  <- load_flat_csv("swei", SWEI_SCALE)
  out <- file.path(BASIN_PLOT_DIR,
                   sprintf("swei_%02d_timeseries.png", SWEI_SCALE))
  ggplot2::ggsave(out, make_basin_ts_plot(df, "swei", SWEI_SCALE),
                  width = 12, height = 6, dpi = 150)
  cat(sprintf("  ✓ SWEI-%02d\n", SWEI_SCALE))
}, error = function(e)
  cat(sprintf("  ⚠ SWEI-%02d: %s\n", SWEI_SCALE, e$message)))

cat("\n  MSPI & MSPEI (if available)...\n")
for (idx in c("mspi", "mspei")) {
  tryCatch({
    df  <- load_flat_csv(idx, MSPI_MSPEI_SCALE)
    out <- file.path(BASIN_PLOT_DIR,
                     sprintf("%s_%02d_timeseries.png", idx, MSPI_MSPEI_SCALE))
    ggplot2::ggsave(out, make_basin_ts_plot(df, idx, MSPI_MSPEI_SCALE),
                    width = 12, height = 6, dpi = 150)
    cat(sprintf("  ✓ %s-%02d\n", toupper(idx), MSPI_MSPEI_SCALE))
  }, error = function(e)
    cat(sprintf("  ⚠ %s-%02d: %s\n", toupper(idx), MSPI_MSPEI_SCALE, e$message)))
}

################################################################################
# PART 4b – MANUSCRIPT FIGURE 4  (two PET-method versions + SPEI-2 fix)
#
# ROOT CAUSE OF SPEI-2 MISSING:
#   SPEI_SCALES in DROUGHT_ANALYSIS_utils.R = c(1, 3, 6, 12).  Scale 2 is
#   never processed by w1_trend_test.R, so no basin-averaged CSV exists.
#   However, 3SPEI_ERALand.R writes per-calendar-month CSV files for every
#   scale it computes (including scale 2) inside spei_results_seasonal/ and
#   spei_results_seasonal_thw/.  load_spei_seasonal_csvs() reads those files
#   directly, computing an area-weighted basin mean on-the-fly.  This bypasses
#   the w1 dependency for any scale not in SPEI_SCALES.
#
# TWO PET METHOD VERSIONS:
#   Version A — SPEI_PM  : Penman-Monteith PET  (spei_results_seasonal/)
#   Version B — SPEI_Thw : Thornthwaite PET      (spei_results_seasonal_thw/)
#   Each version produces its own PDF + PNG.  An overlay version (both on
#   the same panels) is also saved.
#
# Outputs (all in BASIN_PLOT_DIR):
#   Fig4_MS_SPEI_123_PM_timeseries.pdf/.png       — PM only (3 panels)
#   Fig4_MS_SPEI_123_Thw_timeseries.pdf/.png      — Thornthwaite only
#   Fig4_MS_SPEI_123_overlay_timeseries.pdf/.png  — both overlaid per panel
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 4b: Manuscript Figure 4 — SPEI-1/2/3 time series (PM + Thornthwaite)\n")
cat("══════════════════════════════════════════════\n")

# ── Directory paths for the two SPEI variants ──────────────────────────────────
SPEI_PM_DIR  <- file.path(WD_PATH, "spei_results_seasonal")
SPEI_THW_DIR <- file.path(WD_PATH, "spei_results_seasonal_thw")
MONTH_NAMES_3 <- c("Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec")

# ── Direct per-month CSV loader (bypasses w1 for any scale) ───────────────────
# Reads the 12 per-calendar-month CSV files written by 3SPEI_ERALand.R and
# computes an area-weighted basin mean from non-NA pixel values.
# Returns data.frame(date, value) in chronological order, or NULL on failure.
load_spei_seasonal_csvs <- function(scale, seas_dir) {
  
  if (!dir.exists(seas_dir)) {
    cat(sprintf("    [load_spei_seasonal_csvs] Directory not found: %s\n", seas_dir))
    return(NULL)
  }
  
  month_rows <- vector("list", 12L)
  
  for (m in seq_len(12L)) {
    csv_file <- file.path(
      seas_dir,
      sprintf("spei_%02d_month%02d_%s.csv", scale, m, MONTH_NAMES_3[m])
    )
    if (!file.exists(csv_file)) next
    
    tryCatch({
      df <- data.table::fread(csv_file, data.table = FALSE)
      
      # Columns: lon, lat, then one column per year (e.g. "1950", "1951", …)
      year_cols <- setdiff(names(df), c("lon", "lat"))
      year_cols <- year_cols[grepl("^[0-9]{4}$", year_cols)]
      if (!length(year_cols)) return(NULL)
      
      years_int <- as.integer(year_cols)
      
      # Area-weighted basin mean: equal-area projection → simple colMeans over
      # non-NA pixels is a close approximation; proper weighting requires the
      # raster template.  Use colMeans (pixels already masked to basin).
      vals <- colMeans(df[, year_cols, drop = FALSE], na.rm = TRUE)
      vals[is.nan(vals)] <- NA_real_
      
      month_rows[[m]] <- data.frame(
        date  = as.Date(paste(years_int, m, "01", sep = "-")),
        value = as.numeric(vals),
        stringsAsFactors = FALSE
      )
    }, error = function(e)
      cat(sprintf("    ⚠ Could not read %s: %s\n", basename(csv_file), e$message)))
  }
  
  valid <- Filter(Negate(is.null), month_rows)
  if (!length(valid)) return(NULL)
  
  out <- do.call(rbind, valid)
  out <- out[order(out$date), ]
  out <- out[!is.na(out$value) & is.finite(out$value), ]
  rownames(out) <- NULL
  out
}

# ── Shared panel constants ──────────────────────────────────────────────────────
THR_EVENT       <- -0.5
THR_MOD         <- -1.0
HIGHLIGHT_START <- as.Date("2022-01-01")
HIGHLIGHT_END   <- as.Date("2025-12-31")

FILL_WET  <- "#4393c3"
FILL_MILD <- "#ffffbf"
FILL_MOD  <- "#fdae61"
FILL_SEV  <- "#f46d43"
FILL_EXT  <- "#d73027"

fill_pal <- c(
  "Wet"              = FILL_WET,
  "Mild deficit"     = FILL_MILD,
  "Moderate drought" = FILL_MOD,
  "Severe drought"   = FILL_SEV,
  "Extreme drought"  = FILL_EXT
)

# ── Single-version panel builder ───────────────────────────────────────────────
# pet_type : "PM" or "Thw" — controls which directory and line colour
# is_first : TRUE = show legend; FALSE = suppress
make_spei_panel_v2 <- function(sc, panel_lab, pet_type = "PM",
                               show_xlab = FALSE, is_first = FALSE) {
  seas_dir  <- if (pet_type == "PM") SPEI_PM_DIR else SPEI_THW_DIR
  col_line  <- if (pet_type == "PM") "#1f77b4" else "#d62728"
  pet_label <- if (pet_type == "PM") "SPEI\u209a\u2098 (Penman-Monteith)"
  else                  "SPEI\u2080 (Thornthwaite)"
  
  # Try authoritative basin-avg CSV first; fall back to per-month seasonal CSVs
  df <- tryCatch(load_flat_csv("spei", sc), error = function(e) NULL)
  
  # For Thornthwaite or any scale absent from w1, always use the seasonal CSVs
  if (is.null(df) || !nrow(df) || pet_type == "Thw") {
    df <- tryCatch(load_spei_seasonal_csvs(sc, seas_dir), error = function(e) NULL)
  }
  
  if (is.null(df) || !nrow(df))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("No data\nSPEI-%d (%s)", sc, pet_type),
                               size = 4.5, colour = "grey50") +
             ggplot2::theme_void())
  
  df$date  <- as.Date(df$date)
  df$value <- as.numeric(df$value)
  
  df_fill <- df
  df_fill$fill_cat <- factor(dplyr::case_when(
    df_fill$value >  0    ~ "Wet",
    df_fill$value > -0.5  ~ "Mild deficit",
    df_fill$value > -1.0  ~ "Moderate drought",
    df_fill$value > -1.5  ~ "Severe drought",
    TRUE                  ~ "Extreme drought"),
    levels = c("Wet","Mild deficit","Moderate drought",
               "Severe drought","Extreme drought"))
  df_fill$date_end <- pmin(df_fill$date + 31, max(df_fill$date) + 31)
  
  ylo <- min(-2.8, min(df$value, na.rm = TRUE) * 1.05)
  yhi <- max( 2.2, max(df$value, na.rm = TRUE) * 1.05)
  
  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
                      xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                      ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.55) +
    ggplot2::geom_rect(
      data = df_fill,
      ggplot2::aes(xmin = date, xmax = date_end,
                   ymin = -Inf, ymax = Inf, fill = fill_cat),
      alpha = 0.30, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(
      values = fill_pal, name = "Drought category", drop = FALSE,
      guide  = if (is_first)
        ggplot2::guide_legend(title.position = "top", nrow = 1,
                              override.aes = list(alpha = 0.7))
      else "none") +
    ggplot2::geom_hline(yintercept = 0,
                        colour = "grey40", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = THR_EVENT,
                        colour = "grey30", linewidth = 0.50, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = THR_MOD,
                        colour = "grey50", linewidth = 0.35, linetype = "dotted") +
    ggplot2::geom_line(
      data = df, ggplot2::aes(x = date, y = value),
      colour = col_line, linewidth = 0.55, inherit.aes = FALSE) +
    ggplot2::annotate("text",
                      x = min(df$date) + 365, y = yhi * 0.88,
                      label = panel_lab, size = 3.2, fontface = "bold",
                      hjust = 0, vjust = 1, colour = "grey10") +
    ggplot2::annotate("text",
                      x = as.Date("1951-01-01"), y = THR_EVENT + 0.08,
                      label = "x\u2080 = \u22120.5", size = 2.4,
                      hjust = 0, colour = "grey30") +
    ggplot2::annotate("text",
                      x = HIGHLIGHT_START + 365, y = yhi * 0.97,
                      label = "2022\u20132025", size = 2.4, hjust = 0.5,
                      colour = "grey30", fontface = "italic") +
    ggplot2::scale_x_date(
      date_breaks = "10 years", date_labels = "%Y",
      expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::scale_y_continuous(
      limits = c(ylo, yhi),
      breaks = c(-2, -1.5, -1, -0.5, 0, 1, 2),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title = sprintf("SPEI-%d  (%s,  area-weighted basin mean)", sc, pet_label),
      x     = if (show_xlab) "Year" else NULL,
      y     = sprintf("SPEI-%d", sc)) +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 8, face = "bold",
                                           hjust = 0,
                                           margin = ggplot2::margin(b = 2)),
      axis.text.x  = if (show_xlab) ggplot2::element_text(size = 7)
      else ggplot2::element_blank(),
      axis.ticks.x = if (show_xlab) ggplot2::element_line()
      else ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      axis.title.y = ggplot2::element_text(size = 7.5, angle = 90),
      axis.title.x = ggplot2::element_text(size = 7.5),
      legend.position = "bottom",
      legend.text  = ggplot2::element_text(size = 6.5),
      legend.key.size = ggplot2::unit(0.35, "cm"),
      legend.title = ggplot2::element_text(size = 7),
      panel.grid   = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(3, 5, 2, 4))
  
  # SPEI-3 persistence annotation (only on the last scale panel)
  if (sc == 3) {
    p <- p +
      ggplot2::annotate("segment",
                        x = as.Date("2023-01-01"), xend = as.Date("2025-12-31"),
                        y = -0.42, yend = -0.42,
                        colour = "#d73027", linewidth = 0.6) +
      ggplot2::annotate("text",
                        x = as.Date("2024-06-01"), y = -0.35,
                        label = "Persistent below \u22120.5",
                        size = 2.3, hjust = 0.5,
                        colour = "#d73027", fontface = "italic")
  }
  p
}

# ── Overlay panel builder (PM + Thornthwaite on the same axes) ─────────────────
make_spei_overlay_panel <- function(sc, panel_lab, show_xlab = FALSE,
                                    is_first = FALSE) {
  load_one <- function(pet_type) {
    seas_dir <- if (pet_type == "PM") SPEI_PM_DIR else SPEI_THW_DIR
    df <- NULL
    if (pet_type == "PM")
      df <- tryCatch(load_flat_csv("spei", sc), error = function(e) NULL)
    if (is.null(df) || !nrow(df) || pet_type == "Thw")
      df <- tryCatch(load_spei_seasonal_csvs(sc, seas_dir), error = function(e) NULL)
    if (!is.null(df) && nrow(df)) {
      df$date    <- as.Date(df$date)
      df$value   <- as.numeric(df$value)
      df$pet_type <- pet_type
    }
    df
  }
  df_pm  <- load_one("PM")
  df_thw <- load_one("Thw")
  
  if (is.null(df_pm) && is.null(df_thw))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = sprintf("No data\nSPEI-%d", sc),
                               size = 4.5, colour = "grey50") +
             ggplot2::theme_void())
  
  # Use PM for background fill if available, else Thw
  df_fill_src <- if (!is.null(df_pm)) df_pm else df_thw
  df_fill <- df_fill_src
  df_fill$fill_cat <- factor(dplyr::case_when(
    df_fill$value >  0    ~ "Wet",
    df_fill$value > -0.5  ~ "Mild deficit",
    df_fill$value > -1.0  ~ "Moderate drought",
    df_fill$value > -1.5  ~ "Severe drought",
    TRUE                  ~ "Extreme drought"),
    levels = c("Wet","Mild deficit","Moderate drought",
               "Severe drought","Extreme drought"))
  df_fill$date_end <- pmin(df_fill$date + 31, max(df_fill$date) + 31)
  
  all_vals <- c(if (!is.null(df_pm))  df_pm$value,
                if (!is.null(df_thw)) df_thw$value)
  ylo <- min(-2.8, min(all_vals, na.rm = TRUE) * 1.05)
  yhi <- max( 2.2, max(all_vals, na.rm = TRUE) * 1.05)
  
  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
                      xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                      ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.55) +
    ggplot2::geom_rect(
      data = df_fill,
      ggplot2::aes(xmin = date, xmax = date_end,
                   ymin = -Inf, ymax = Inf, fill = fill_cat),
      alpha = 0.25, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(
      values = fill_pal, name = "Drought category", drop = FALSE,
      guide  = if (is_first)
        ggplot2::guide_legend(title.position = "top", nrow = 1,
                              override.aes = list(alpha = 0.6))
      else "none") +
    ggplot2::geom_hline(yintercept = 0,     colour = "grey40", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = THR_EVENT, colour = "grey30",
                        linewidth = 0.50, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = THR_MOD, colour = "grey50",
                        linewidth = 0.35, linetype = "dotted")
  
  if (!is.null(df_pm))
    p <- p + ggplot2::geom_line(
      data = df_pm, ggplot2::aes(x = date, y = value,
                                 colour = "SPEI\u209a\u2098 (Penman-Monteith)"),
      linewidth = 0.65, inherit.aes = FALSE)
  
  if (!is.null(df_thw))
    p <- p + ggplot2::geom_line(
      data = df_thw, ggplot2::aes(x = date, y = value,
                                  colour = "SPEI\u2080 (Thornthwaite)"),
      linewidth = 0.55, linetype = "solid", alpha = 0.75, inherit.aes = FALSE)
  
  p <- p +
    ggplot2::scale_colour_manual(
      values = c("SPEI\u209a\u2098 (Penman-Monteith)" = "#1f77b4",
                 "SPEI\u2080 (Thornthwaite)"          = "#d62728"),
      name   = "PET method",
      guide  = if (is_first)
        ggplot2::guide_legend(title.position = "top", nrow = 1,
                              override.aes = list(linewidth = 1.2))
      else "none") +
    ggplot2::annotate("text",
                      x = min(df_fill_src$date) + 365, y = yhi * 0.88,
                      label = panel_lab, size = 3.2, fontface = "bold",
                      hjust = 0, vjust = 1, colour = "grey10") +
    ggplot2::annotate("text",
                      x = as.Date("1951-01-01"), y = THR_EVENT + 0.08,
                      label = "x\u2080 = \u22120.5", size = 2.4,
                      hjust = 0, colour = "grey30") +
    ggplot2::annotate("text",
                      x = HIGHLIGHT_START + 365, y = yhi * 0.97,
                      label = "2022\u20132025", size = 2.4, hjust = 0.5,
                      colour = "grey30", fontface = "italic") +
    ggplot2::scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                          expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::scale_y_continuous(limits = c(ylo, yhi),
                                breaks = c(-2, -1.5, -1, -0.5, 0, 1, 2),
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title = sprintf("SPEI-%d  (area-weighted basin mean)", sc),
      x     = if (show_xlab) "Year" else NULL,
      y     = sprintf("SPEI-%d", sc)) +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 8, face = "bold", hjust = 0,
                                           margin = ggplot2::margin(b = 2)),
      axis.text.x  = if (show_xlab) ggplot2::element_text(size = 7)
      else ggplot2::element_blank(),
      axis.ticks.x = if (show_xlab) ggplot2::element_line()
      else ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 7),
      axis.title.y = ggplot2::element_text(size = 7.5, angle = 90),
      axis.title.x = ggplot2::element_text(size = 7.5),
      legend.position = "bottom",
      legend.text  = ggplot2::element_text(size = 6.5),
      legend.key.size = ggplot2::unit(0.35, "cm"),
      legend.title = ggplot2::element_text(size = 7),
      panel.grid   = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(3, 5, 2, 4))
  
  if (sc == 3) {
    p <- p +
      ggplot2::annotate("segment",
                        x = as.Date("2023-01-01"), xend = as.Date("2025-12-31"),
                        y = -0.42, yend = -0.42,
                        colour = "#d73027", linewidth = 0.6) +
      ggplot2::annotate("text",
                        x = as.Date("2024-06-01"), y = -0.35,
                        label = "Persistent below \u22120.5",
                        size = 2.3, hjust = 0.5,
                        colour = "#d73027", fontface = "italic")
  }
  p
}

# ── Helper: assemble + save a 3-panel figure ───────────────────────────────────
save_spei_fig <- function(panels, title_str, subtitle_str, base_name) {
  fig <- panels[[1]] / panels[[2]] / panels[[3]] +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title    = title_str,
      subtitle = subtitle_str,
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40",
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 4)),
        legend.position = "bottom",
        legend.text     = ggplot2::element_text(size = 7.5),
        legend.key.size = ggplot2::unit(0.38, "cm")))
  
  for (ext in c("pdf", "png")) {
    out <- file.path(BASIN_PLOT_DIR, paste0(base_name, ".", ext))
    tryCatch(
      ggplot2::ggsave(out, fig, width = 7.0, height = 9.0, units = "in",
                      dpi = if (ext == "png") 300 else "print",
                      device = ext),
      error = function(e) cat(sprintf("  \u26a0 %s: %s\n", ext, e$message)))
    cat(sprintf("  \u2713 Saved: %s\n", basename(out)))
  }
  invisible(fig)
}

common_subtitle <- paste0(
  "Fill: drought category  |  ",
  "Dashed: event threshold (x\u2080 = \u22120.5, 31st percentile)  |  ",
  "Dotted: moderate drought (\u22121.0)  |  ",
  "Grey band: 2022\u20132025  |  ",
  "Red segment: persistent sub-threshold anomaly 2023\u20132025 (SPEI-3)")

tryCatch({
  # ── Version A: PM only ─────────────────────────────────────────────────────
  cat("  Building Version A: SPEI-PM (Penman-Monteith)...\n")
  pm_panels <- list(
    make_spei_panel_v2(1, "(a)", "PM", show_xlab = FALSE, is_first = TRUE),
    make_spei_panel_v2(2, "(b)", "PM", show_xlab = FALSE, is_first = FALSE),
    make_spei_panel_v2(3, "(c)", "PM", show_xlab = TRUE,  is_first = FALSE)
  )
  save_spei_fig(pm_panels,
                "Basin-averaged SPEI\u209a\u2098 (Penman-Monteith PET) \u2014 Nechako River Basin (1950\u20132025)",
                common_subtitle,
                "Fig4_MS_SPEI_123_PM_timeseries")
  cat("  \u2713 Version A complete\n")
}, error = function(e) cat(sprintf("  \u26a0 Version A failed: %s\n", e$message)))

tryCatch({
  # ── Version B: Thornthwaite only ──────────────────────────────────────────
  cat("  Building Version B: SPEI\u2080 (Thornthwaite)...\n")
  thw_panels <- list(
    make_spei_panel_v2(1, "(a)", "Thw", show_xlab = FALSE, is_first = TRUE),
    make_spei_panel_v2(2, "(b)", "Thw", show_xlab = FALSE, is_first = FALSE),
    make_spei_panel_v2(3, "(c)", "Thw", show_xlab = TRUE,  is_first = FALSE)
  )
  save_spei_fig(thw_panels,
                "Basin-averaged SPEI\u2080 (Thornthwaite PET) \u2014 Nechako River Basin (1950\u20132025)",
                common_subtitle,
                "Fig4_MS_SPEI_123_Thw_timeseries")
  cat("  \u2713 Version B complete\n")
}, error = function(e) cat(sprintf("  \u26a0 Version B failed: %s\n", e$message)))

tryCatch({
  # ── Version C: Overlay (PM + Thornthwaite on same axes) ───────────────────
  cat("  Building Version C: overlay (PM + Thornthwaite)...\n")
  ov_panels <- list(
    make_spei_overlay_panel(1, "(a)", show_xlab = FALSE, is_first = TRUE),
    make_spei_overlay_panel(2, "(b)", show_xlab = FALSE, is_first = FALSE),
    make_spei_overlay_panel(3, "(c)", show_xlab = TRUE,  is_first = FALSE)
  )
  save_spei_fig(ov_panels,
                "Basin-averaged SPEI \u2014 Nechako River Basin (1950\u20132025)",
                paste0(common_subtitle,
                       "  |  Blue = PM;  Red = Thornthwaite"),
                "Fig4_MS_SPEI_123_overlay_timeseries")
  cat("  \u2713 Version C complete\n")
}, error = function(e) cat(sprintf("  \u26a0 Version C failed: %s\n", e$message)))

################################################################################
# PART 4c – SPATIAL REPRESENTATIVENESS MAPS  (SPEI-1/2/3 × PM + Thornthwaite)
#
# Outputs two figures:
#
#   Fig4c_Pearson_r_SPEI123_PM_Thw.pdf/.png
#       2-row × 3-column grid of Pearson r maps.
#       Row 1 = SPEI PM  (Penman-Monteith),  columns = scales 1, 2, 3
#       Row 2 = SPEI Thw (Thornthwaite),     columns = scales 1, 2, 3
#       Each panel is a per-pixel correlation with the basin-mean time series.
#       Annotated with median r and % pixels ≥ 0.80.
#
#   Fig4c_SD_SpatialSD_SPEI3_PM.pdf/.png
#       3-panel supplementary figure for SPEI-3 PM only:
#       (a) Temporal SD map   (b) Annual spatial SD bar chart
#       (c) histogram of r values across all 6 combinations
#
# All panels use ggplot2 + patchwork (no terra::plot).
################################################################################
cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("PART 4c: Spatial representativeness \u2014 SPEI-1/2/3 \u00d7 PM + Thornthwaite\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")

tryCatch({
  
  # ── Helper: load all 12 per-month CSVs for one scale × directory ───────────
  load_pixel_mat <- function(scale, seas_dir) {
    month_list  <- vector("list", 12L)
    coords_out  <- NULL
    year_out    <- NULL
    for (m in seq_len(12L)) {
      csv_f <- file.path(seas_dir,
                         sprintf("spei_%02d_month%02d_%s.csv", scale, m, MONTH_NAMES_3[m]))
      if (!file.exists(csv_f)) next
      df_m <- tryCatch(data.table::fread(csv_f, data.table = FALSE),
                       error = function(e) NULL)
      if (is.null(df_m)) next
      yr_c <- setdiff(names(df_m), c("lon", "lat"))
      yr_c <- yr_c[grepl("^[0-9]{4}$", yr_c)]
      if (!length(yr_c)) next
      if (is.null(coords_out)) { coords_out <- df_m[, c("lon","lat")]; year_out <- yr_c }
      month_list[[m]] <- as.matrix(df_m[, yr_c, drop = FALSE])
    }
    valid <- which(!sapply(month_list, is.null))
    if (!length(valid) || is.null(coords_out)) return(NULL)
    list(mat    = do.call(cbind, month_list[valid]),
         coords = coords_out,
         yr_int = as.integer(year_out),
         yr_chr = year_out,
         valid_months = valid,
         month_list = month_list)
  }
  
  # ── Helper: compute metrics from a loaded pixel matrix ────────────────────
  compute_metrics <- function(obj) {
    mat  <- obj$mat
    bavg <- colMeans(mat, na.rm = TRUE)
    r_v  <- vapply(seq_len(nrow(mat)), function(i) {
      xi <- mat[i,]; ok <- !is.na(xi) & !is.na(bavg)
      if (sum(ok) < 10L) NA_real_ else cor(xi[ok], bavg[ok])
    }, numeric(1L))
    sd_v <- apply(mat, 1, sd, na.rm = TRUE)
    
    # Annual spatial SD
    yr_chr <- obj$yr_chr
    n_yrs  <- length(yr_chr)
    n_pix  <- nrow(mat)
    ann_m  <- matrix(0, n_pix, n_yrs)
    vm     <- obj$valid_months
    ml     <- obj$month_list
    for (m in vm) {
      mm <- ml[[m]]
      cy <- intersect(colnames(mm), yr_chr)
      ann_m[, match(cy, yr_chr)] <- ann_m[, match(cy, yr_chr)] +
        mm[, match(cy, colnames(mm)), drop = FALSE] / length(vm)
    }
    spat_sd <- apply(ann_m, 2, sd, na.rm = TRUE)
    
    list(r = r_v, sd = sd_v, spat_sd = spat_sd,
         med_r = median(r_v, na.rm = TRUE),
         pct80 = 100 * mean(r_v >= 0.80, na.rm = TRUE),
         mean_spat_sd = mean(spat_sd, na.rm = TRUE))
  }
  
  # ── CRS helper ────────────────────────────────────────────────────────────
  detect_basin_sf <- function(coords) {
    is_proj <- max(abs(coords$lon), na.rm = TRUE) > 200
    tryCatch({
      b <- sf::st_read(BASIN_SHP, quiet = TRUE)
      if (is_proj) sf::st_transform(b, terra::crs(basin_proj))
      else         sf::st_transform(b, "EPSG:4326")
    }, error = function(e) NULL)
  }
  
  # ── Shared map theme ──────────────────────────────────────────────────────
  theme_map4c <- ggplot2::theme_bw(base_size = 8.5) +
    ggplot2::theme(
      axis.title      = ggplot2::element_blank(),
      axis.text       = ggplot2::element_text(size = 6),
      plot.title      = ggplot2::element_text(size = 8, face = "bold", hjust = 0.5),
      plot.subtitle   = ggplot2::element_text(size = 7, colour = "grey30",
                                              hjust = 0.5,
                                              margin = ggplot2::margin(t = 1, b = 2)),
      legend.key.height = ggplot2::unit(1.1, "cm"),
      legend.key.width  = ggplot2::unit(0.32, "cm"),
      legend.text     = ggplot2::element_text(size = 6.5),
      legend.title    = ggplot2::element_text(size = 7),
      panel.grid      = ggplot2::element_blank(),
      plot.margin     = ggplot2::margin(2, 2, 2, 2)
    )
  
  # ── Build one Pearson-r map panel ─────────────────────────────────────────
  make_r_panel <- function(scale, pet_type) {
    seas_dir <- if (pet_type == "PM") SPEI_PM_DIR else SPEI_THW_DIR
    obj <- load_pixel_mat(scale, seas_dir)
    if (is.null(obj))
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x=0.5, y=0.5,
                                 label = sprintf("No data\nSPEI-%d %s", scale, pet_type),
                                 size = 3, colour = "grey50") +
               ggplot2::theme_void())
    
    mt       <- compute_metrics(obj)
    basin_sf <- detect_basin_sf(obj$coords)
    pet_lab  <- if (pet_type == "PM") "SPEI\u209a\u2098" else "SPEI\u2080"
    r_min    <- max(0.40, floor(min(mt$r, na.rm = TRUE) * 10) / 10)
    
    df_r <- data.frame(lon = obj$coords$lon, lat = obj$coords$lat, r = mt$r)
    
    p <- ggplot2::ggplot(df_r[!is.na(df_r$r), ],
                         ggplot2::aes(x = lon, y = lat, fill = r)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_distiller(
        palette = "Blues", direction = 1,
        name    = "Pearson r",
        limits  = c(r_min, 1.0),
        breaks  = seq(r_min, 1.0, by = 0.1),
        oob     = scales::squish) +
      {if (!is.null(basin_sf))
        ggplot2::geom_sf(data = basin_sf, fill = NA,
                         colour = "black", linewidth = 0.55,
                         inherit.aes = FALSE)} +
      ggplot2::coord_sf(expand = FALSE) +
      ggplot2::labs(
        title    = sprintf("%s-%d  (%s)",  pet_lab, scale,
                           if (pet_type=="PM") "Penman-Monteith" else "Thornthwaite"),
        subtitle = sprintf("median r = %.3f  |  %.0f%% \u2265 0.80",
                           mt$med_r, mt$pct80)) +
      theme_map4c
    
    list(plot = p, metrics = mt, obj = obj)
  }
  
  # ── Run all 6 combinations ─────────────────────────────────────────────────
  combos <- list(
    list(sc=1, pet="PM"),  list(sc=2, pet="PM"),  list(sc=3, pet="PM"),
    list(sc=1, pet="Thw"), list(sc=2, pet="Thw"), list(sc=3, pet="Thw")
  )
  
  results <- lapply(combos, function(cb) {
    cat(sprintf("  Building SPEI-%d %s...\n", cb$sc, cb$pet))
    make_r_panel(cb$sc, cb$pet)
  })
  
  panels_pm  <- lapply(results[1:3], `[[`, "plot")
  panels_thw <- lapply(results[4:6], `[[`, "plot")
  
  # Check we have at least some panels
  any_data <- any(!sapply(results, is.null))
  if (!any_data) { cat("  \u26a0 No data for any combination\n"); stop("no_data") }
  
  # ── Figure 1: flat 6-panel Pearson r map (2 rows × 3 cols) ──────────────────
  # Labels: (a)-(c) = PM scales 1-3,  (d)-(f) = Thornthwaite scales 1-3
  panel_labels <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
  all_six_plots <- lapply(seq_along(results), function(k) {
    cb  <- combos[[k]]
    res <- results[[k]]
    pet_str <- if (cb$pet == "PM") "Penman-Monteith" else "Thornthwaite"
    sub_str <- if (!is.null(res))
      sprintf("median r = %.3f  |  %.0f%% \u2265 0.80",
              res$metrics$med_r, res$metrics$pct80)
    else "no data"
    
    if (is.null(res))
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x=0.5, y=0.5,
                                 label=sprintf("No data\nSPEI-%d %s", cb$sc, cb$pet),
                                 size=3, colour="grey50") + ggplot2::theme_void())
    
    # re-use the plot but override title to include panel label + PET method
    pet_lab <- if (cb$pet == "PM") "SPEI\u209a\u2098" else "SPEI\u2080"
    p <- res$plot +
      ggplot2::labs(
        title    = sprintf("%s  %s-%d  (%s)",
                           panel_labels[k], pet_lab, cb$sc, pet_str),
        subtitle = sub_str)
    p
  })
  
  fig4c_6panel <- patchwork::wrap_plots(all_six_plots, nrow = 2, ncol = 3,
                                        guides = "collect") +
    patchwork::plot_annotation(
      title    = paste0("Pixel-to-basin-mean Pearson r \u2014 SPEI 1, 2, 3",
                        " (Penman-Monteith vs Thornthwaite)  \u2014  Nechako River Basin"),
      subtitle = paste0(
        "Each pixel coloured by its temporal Pearson r with the basin-mean SPEI time series.",
        "  Row 1: Penman-Monteith PET.  Row 2: Thornthwaite PET.\n",
        "r \u2248 1 everywhere \u2192 the single basin mean faithfully represents",
        " drought variability at every point in the basin."),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey30", hjust = 0,
                                              margin = ggplot2::margin(b = 6)),
        legend.position = "right"))
  
  for (ext in c("pdf", "png")) {
    out_r <- file.path(BASIN_PLOT_DIR,
                       paste0("Fig4c_Pearson_r_SPEI123_PM_Thw.", ext))
    tryCatch(
      ggplot2::ggsave(out_r, fig4c_6panel,
                      width  = 14, height = 9.5, units = "in",
                      dpi = if (ext == "png") 300 else "print",
                      device = ext),
      error = function(e) cat(sprintf("  \u26a0 %s: %s\n", ext, e$message)))
    cat(sprintf("  \u2713 Saved: %s\n", basename(out_r)))
  }
  
  # ── Figure 2: SD map + bar chart + r-histogram (SPEI-3 PM) ───────────────
  # Pick the SPEI-3 PM result (index 3 in results list)
  r3pm <- results[[3]]
  if (!is.null(r3pm)) {
    mt3   <- r3pm$metrics
    obj3  <- r3pm$obj
    bsf3  <- detect_basin_sf(obj3$coords)
    
    # SD map
    df_sd <- data.frame(lon = obj3$coords$lon, lat = obj3$coords$lat,
                        sd_t = mt3$sd)
    pa_sd <- ggplot2::ggplot(df_sd[!is.na(df_sd$sd_t),],
                             ggplot2::aes(x=lon, y=lat, fill=sd_t)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_distiller(palette="YlOrRd", direction=1,
                                    name="Temporal SD",
                                    limits=c(max(0.7, min(mt3$sd,na.rm=T)-0.01),
                                             min(1.3, max(mt3$sd,na.rm=T)+0.01)),
                                    oob=scales::squish) +
      {if(!is.null(bsf3)) ggplot2::geom_sf(data=bsf3, fill=NA, colour="black",
                                           linewidth=0.55, inherit.aes=FALSE)} +
      ggplot2::coord_sf(expand=FALSE) +
      ggplot2::labs(title="(a)  Temporal SD of SPEI-3 PM") +
      theme_map4c
    
    # Annual spatial SD bar chart
    df_bar <- data.frame(year=obj3$yr_int, spat_sd=mt3$spat_sd)
    pb_bar <- ggplot2::ggplot(df_bar, ggplot2::aes(x=year, y=spat_sd)) +
      ggplot2::geom_col(fill="#4393c3", colour=NA, width=0.9) +
      ggplot2::geom_hline(yintercept=mt3$mean_spat_sd,
                          colour="#d73027", linewidth=0.9, linetype="dashed") +
      ggplot2::annotate("text",
                        x=max(obj3$yr_int)-2, y=mt3$mean_spat_sd*1.10,
                        label=sprintf("Mean = %.3f", mt3$mean_spat_sd),
                        colour="#d73027", size=2.8, hjust=1) +
      ggplot2::scale_x_continuous(breaks=seq(1950,2030,by=10),
                                  expand=ggplot2::expansion(add=c(0.5,0.5))) +
      ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,0.15))) +
      ggplot2::labs(title="(b)  Spatial SD across pixels per year (SPEI-3 PM)",
                    x="Year", y="Spatial SD") +
      ggplot2::theme_classic(base_size=9) +
      ggplot2::theme(plot.title=ggplot2::element_text(size=8,face="bold",hjust=0),
                     axis.text.x=ggplot2::element_text(angle=45,hjust=1,size=7),
                     axis.text.y=ggplot2::element_text(size=7),
                     plot.margin=ggplot2::margin(2,6,2,2))
    
    # ── Shared helper: build one r-density panel for a single PET method ─────
    # pet_type : "PM" or "Thw"
    # panel_lab: "(c)" or "(d)"
    make_r_density_panel <- function(pet_type, panel_lab) {
      col_pm  <- c("SPEI-1 PM"  = "#1f77b4",
                   "SPEI-2 PM"  = "#6baed6",
                   "SPEI-3 PM"  = "#2ca02c")
      col_thw <- c("SPEI-1 Thw" = "#d62728",
                   "SPEI-2 Thw" = "#fd8d3c",
                   "SPEI-3 Thw" = "#ff7f0e")
      col_use  <- if (pet_type == "PM") col_pm else col_thw
      pet_long <- if (pet_type == "PM") "Penman-Monteith" else "Thornthwaite"
      
      # Collect r values for this PET method only
      df <- do.call(rbind, lapply(seq_along(combos), function(k) {
        cb  <- combos[[k]]; res <- results[[k]]
        if (cb$pet != pet_type || is.null(res)) return(NULL)
        data.frame(r     = res$metrics$r,
                   combo = sprintf("SPEI-%d %s", cb$sc, cb$pet),
                   stringsAsFactors = FALSE)
      }))
      if (is.null(df) || !nrow(df)) return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x=0.5, y=0.5,
                            label=sprintf("No data\n(%s)", pet_long),
                            size=3, colour="grey50") +
          ggplot2::theme_void())
      
      df <- df[!is.na(df$r), ]
      df$combo <- factor(df$combo, levels = names(col_use))
      
      # Per-group medians
      med_df <- do.call(rbind, lapply(levels(df$combo), function(g) {
        x <- df$r[df$combo == g]
        data.frame(combo  = g,
                   med_r  = median(x, na.rm = TRUE),
                   stringsAsFactors = FALSE)
      }))
      med_df$combo <- factor(med_df$combo, levels = names(col_use))
      
      ggplot2::ggplot(df, ggplot2::aes(x = r, colour = combo, fill = combo)) +
        ggplot2::geom_density(alpha = 0.18, linewidth = 0.80) +
        ggplot2::geom_vline(data = med_df,
                            ggplot2::aes(xintercept = med_r, colour = combo),
                            linetype = "solid", linewidth = 0.60,
                            show.legend = FALSE) +
        ggplot2::geom_text(data = med_df,
                           ggplot2::aes(x     = med_r, y = Inf,
                                        label = sprintf("%.2f", med_r),
                                        colour = combo),
                           vjust = 1.5, hjust = -0.10, size = 2.4,
                           show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = col_use, name = NULL) +
        ggplot2::scale_fill_manual(values   = col_use, name = NULL) +
        ggplot2::scale_x_continuous(limits = c(0.4, 1.01),
                                    breaks = seq(0.4, 1.0, by = 0.1)) +
        ggplot2::labs(
          title    = sprintf("%s  Pixel Pearson r — %s PET", panel_lab, pet_long),
          subtitle = "Solid vertical lines = median r per SPEI scale",
          x        = "Pearson r  (pixel vs basin mean)",
          y        = "Density") +
        ggplot2::theme_classic(base_size = 9) +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(size=8, face="bold", hjust=0),
          plot.subtitle = ggplot2::element_text(size=7, colour="grey40", hjust=0),
          legend.text   = ggplot2::element_text(size=7.5),
          legend.position  = "bottom",
          legend.key.size  = ggplot2::unit(0.35,"cm"),
          axis.text        = ggplot2::element_text(size=7),
          plot.margin      = ggplot2::margin(3, 6, 2, 4))
    }
    
    pc_pm  <- make_r_density_panel("PM",  "(c)")
    pc_thw <- make_r_density_panel("Thw", "(d)")
    
    # ── 2×2 layout ────────────────────────────────────────────────────────────
    # Row 1: (a) SD map | (b) annual spatial SD bar chart
    # Row 2: (c) PM r-density | (d) Thornthwaite r-density
    fig4c_sd <- (pa_sd | pb_bar) / (pc_pm | pc_thw) +
      patchwork::plot_layout(heights = c(1, 1)) +
      patchwork::plot_annotation(
        title    = paste0("Spatial representativeness diagnostics — SPEI-3",
                          " (PM \u0026 Thornthwaite)"),
        subtitle = paste0(
          "Row 1: temporal SD map (a) and annual spatial SD bar chart (b)",
          " for SPEI-3 PM.\n",
          "Row 2: distribution of pixel Pearson r with the basin mean,",
          " split by PET method (c = PM,  d = Thornthwaite, scales 1\u20133)."),
        theme = ggplot2::theme(
          plot.title    = ggplot2::element_text(size=10, face="bold", hjust=0.5),
          plot.subtitle = ggplot2::element_text(size=8, colour="grey35", hjust=0,
                                                margin=ggplot2::margin(b=5))))
    
    for (ext in c("pdf","png")) {
      out_sd <- file.path(BASIN_PLOT_DIR,
                          paste0("Fig4c_SD_SpatialSD_SPEI3_PM.", ext))
      tryCatch(
        ggplot2::ggsave(out_sd, fig4c_sd,
                        width=10, height=8.5, units="in",
                        dpi = if (ext == "png") 300 else "print", device=ext),
        error=function(e) cat(sprintf("  \u26a0 %s: %s\n", ext, e$message)))
      cat(sprintf("  \u2713 Saved: %s\n", basename(out_sd)))
    }
  }
  
  # Console summary
  cat("\n  Representativeness summary:\n")
  cat(sprintf("  %-18s  median_r  pct_r\u226580\n", "Combination"))
  cat("  ", paste(rep("-",40), collapse=""), "\n", sep="")
  for (k in seq_along(combos)) {
    cb <- combos[[k]]; res <- results[[k]]
    if (is.null(res)) next
    cat(sprintf("  SPEI-%d %-12s   %.3f    %.0f%%\n",
                cb$sc, cb$pet, res$metrics$med_r, res$metrics$pct80))
  }
  
}, error = function(e) {
  if (!identical(conditionMessage(e), "no_data"))
    cat(sprintf("  \u26a0 Part 4c failed: %s\n", e$message))
})


################################################################################
# PART 4d – SEASONALITY FIGURE  (P vs PET(PM) — 12 monthly means)
#
# Two-panel bar chart:
#   Panel (a): mean monthly precipitation P vs PET(PM)  [mm month⁻¹],
#              exactly 12 bars per variable averaged over all 1950–2025 years.
#   Panel (b): climatic water balance P − PET(PM)  [mm month⁻¹].
#
# DATA
#   Precipitation: monthly_data_direct/total_precipitation_monthly.nc
#                  (m/day mean rate) → × days_in_month × 1000 = mm month⁻¹
#                  Same conversion as 1SPI_ERALand.R.
#   PET(PM):       monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv
#                  (mean_pet in mm day⁻¹) → × mean-days-per-month = mm month⁻¹
#                  Written by 2preq_PET_ERALand.R.
#
# OUTPUT  FigS_Seasonality_P_PET.pdf + .png  in BASIN_PLOT_DIR
################################################################################
cat("\n\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("PART 4d: Seasonality figure \u2014 P vs PET(PM) (12 monthly means)\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")

tryCatch({
  
  PRECIP_NC_4d  <- file.path(WD_PATH, "monthly_data_direct",
                             "total_precipitation_monthly.nc")
  PET_SEAS_CSV  <- file.path(WD_PATH, "monthly_data_direct",
                             "ERA5Land_Nechako_PET_monthly_summary.csv")
  
  if (!file.exists(PRECIP_NC_4d))
    stop(sprintf("Precip NetCDF not found: %s", PRECIP_NC_4d))
  if (!file.exists(PET_SEAS_CSV)) {
    cat(sprintf("  \u26a0 PET CSV not found: %s\n", PET_SEAS_CSV))
    cat("    Run 2preq_PET_ERALand.R first, then re-run this script.\n")
    stop("pet_csv_missing_4d")
  }
  
  COL_P_4d   <- "#4393c3"
  COL_PET_4d <- "#d7301f"
  # Mean days per calendar month (non-leap average)
  DAYS_PER_MON_4d <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  ## ── Precipitation from NetCDF (mm month⁻¹) ──────────────────────────────
  ## Mirrors the conversion in 1SPI_ERALand.R exactly:
  ##   raw value = m/day mean rate → × actual days_in_month → m/month → × 1000 mm
  cat("  Loading precipitation from NetCDF...\n")
  pr_r  <- terra::rast(PRECIP_NC_4d)
  pr_dates <- tryCatch({
    t <- terra::time(pr_r)
    if (!all(is.na(t))) as.Date(t)
    else {
      nc  <- ncdf4::nc_open(PRECIP_NC_4d)
      out <- if ("valid_time" %in% names(nc$var)) {
        tv <- ncdf4::ncvar_get(nc, "valid_time")
        as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC"))
      } else {
        seq(as.Date("1950-01-01"), by = "month", length.out = terra::nlyr(pr_r))
      }
      ncdf4::nc_close(nc)
      out
    }
  }, error = function(e)
    seq(as.Date("1950-01-01"), by = "month", length.out = terra::nlyr(pr_r)))
  
  # Reproject and mask to basin (same CRS as the rest of the script)
  basin_4d <- terra::project(terra::vect(BASIN_SHP), EQUAL_AREA_CRS)
  if (!terra::same.crs(pr_r, EQUAL_AREA_CRS))
    pr_r <- terra::project(pr_r, EQUAL_AREA_CRS, method = "bilinear")
  pr_r <- terra::mask(pr_r, basin_4d)
  
  # days_in_month from actual dates
  fom_4d  <- as.Date(format(pr_dates, "%Y-%m-01"))
  fnm_4d  <- seq(fom_4d[1], by = "month", length.out = length(pr_dates) + 1)[-1]
  dim_4d  <- as.integer(fnm_4d - fom_4d)
  pr_r    <- pr_r * dim_4d * 1000   # mm month⁻¹
  
  # Basin mean per layer
  pr_vals <- terra::global(pr_r, "mean", na.rm = TRUE)$mean
  prec_df_4d <- data.frame(
    date  = pr_dates,
    value = as.numeric(pr_vals),
    month = as.integer(format(pr_dates, "%m")))
  
  cm_prec_4d <- prec_df_4d |>
    dplyr::group_by(month) |>
    dplyr::summarise(mean_val = mean(value, na.rm = TRUE),
                     sd_val   = sd(value,   na.rm = TRUE),
                     n_yrs    = dplyr::n(),
                     .groups  = "drop") |>
    dplyr::mutate(
      month_lab = factor(month.abb[month], levels = month.abb),
      variable  = "Precipitation (P)")
  
  ## ── PET(PM) from CSV (mm month⁻¹) ────────────────────────────────────────
  pet_raw_4d <- as.data.frame(data.table::fread(PET_SEAS_CSV))
  if (!inherits(pet_raw_4d$date, "Date"))
    pet_raw_4d$date <- as.Date(
      paste0(substr(as.character(pet_raw_4d$date), 1, 7), "-01"))
  pet_raw_4d$month  <- as.integer(format(pet_raw_4d$date, "%m"))
  pet_raw_4d$pet_mm <- pet_raw_4d$mean_pet * DAYS_PER_MON_4d[pet_raw_4d$month]
  
  cm_pet_4d <- pet_raw_4d |>
    dplyr::group_by(month) |>
    dplyr::summarise(mean_val = mean(pet_mm, na.rm = TRUE),
                     sd_val   = sd(pet_mm,   na.rm = TRUE),
                     n_yrs    = dplyr::n(),
                     .groups  = "drop") |>
    dplyr::mutate(
      month_lab = factor(month.abb[month], levels = month.abb),
      variable  = "PET(PM)")
  
  ann_P   <- round(sum(cm_prec_4d$mean_val))
  ann_PET <- round(sum(cm_pet_4d$mean_val))
  ai      <- round(ann_PET / ann_P, 2)
  nyrs    <- cm_prec_4d$n_yrs[1]
  
  cat(sprintf("  Precip  : %d yrs, mean annual = %d mm yr\u207b\u00b9\n", nyrs, ann_P))
  cat(sprintf("  PET(PM) : %d yrs, mean annual = %d mm yr\u207b\u00b9\n",
              cm_pet_4d$n_yrs[1], ann_PET))
  cat(sprintf("  Aridity index (PET/P) = %.2f\n", ai))
  
  ## ── Panel (a): grouped bar chart P vs PET(PM) ────────────────────────────
  cm_all_4d          <- rbind(cm_prec_4d, cm_pet_4d)
  cm_all_4d$variable <- factor(cm_all_4d$variable,
                               levels = c("Precipitation (P)", "PET(PM)"))
  
  p4d_a <- ggplot2::ggplot(
    cm_all_4d,
    ggplot2::aes(x = month_lab, y = mean_val,
                 fill = variable, group = variable)) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(0.72),
      width = 0.65, colour = "white", linewidth = 0.15) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
      position = ggplot2::position_dodge(0.72),
      width = 0.22, colour = "grey25", linewidth = 0.40) +
    ggplot2::scale_fill_manual(
      values = c("Precipitation (P)" = COL_P_4d, "PET(PM)" = COL_PET_4d),
      name   = NULL) +
    ggplot2::scale_y_continuous(
      name   = expression(mm~month^{-1}),
      expand = ggplot2::expansion(mult = c(0, 0.09))) +
    ggplot2::labs(
      title    = "(a)  Mean monthly precipitation and PET(PM)",
      subtitle = paste0(
        nyrs, "-year mean (1950\u20132025)  |  Error bars = \u00b11 SD  |  ",
        "Annual totals:  P = ", ann_P, " mm;  PET(PM) = ", ann_PET, " mm"),
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
  
  ## ── Panel (b): water balance P − PET(PM) ─────────────────────────────────
  cm_wb_4d <- dplyr::inner_join(
    cm_prec_4d[, c("month", "month_lab", "mean_val")],
    cm_pet_4d[,  c("month", "mean_val")],
    by = "month", suffix = c("_p", "_pet")) |>
    dplyr::mutate(wb       = mean_val_p - mean_val_pet,
                  fill_col = ifelse(wb >= 0, COL_P_4d, COL_PET_4d))
  
  # x positions of surplus and deficit groups (for annotations outside bars)
  surplus_x  <- mean(which(cm_wb_4d$wb >= 0))
  deficit_x  <- mean(which(cm_wb_4d$wb <  0))
  # y positions placed just outside the zero line so they never overlap a bar
  surplus_y  <-  max(cm_wb_4d$wb[cm_wb_4d$wb >= 0], na.rm = TRUE) * 1.18
  deficit_y  <- -abs(min(cm_wb_4d$wb, na.rm = TRUE)) * 1.14
  
  p4d_b <- ggplot2::ggplot(
    cm_wb_4d,
    ggplot2::aes(x = month_lab, y = wb, fill = fill_col)) +
    ggplot2::geom_col(colour = "white", linewidth = 0.15, width = 0.65) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.55) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_continuous(
      name   = expression(P - PET(PM)~~(mm~month^{-1})),
      expand = ggplot2::expansion(mult = c(0.16, 0.16))) +
    # Labels placed ABOVE the tallest surplus bar and BELOW the deepest deficit bar
    ggplot2::annotate("text",
                      x = surplus_x, y = surplus_y,
                      label = "P surplus", colour = COL_P_4d,
                      size = 3.0, fontface = "italic") +
    ggplot2::annotate("text",
                      x = deficit_x, y = deficit_y,
                      label = "PET deficit", colour = COL_PET_4d,
                      size = 3.0, fontface = "italic") +
    ggplot2::labs(
      title    = "(b)  Climatic water balance  P \u2212 PET(PM)",
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
  
  ## ── Assemble and save ────────────────────────────────────────────────────
  fig_seas_4d <- p4d_a / p4d_b +
    patchwork::plot_annotation(
      title    = paste0(
        "Nechako River Basin \u2014 Seasonal cycle of precipitation ",
        "and PET(PM)"),
      subtitle = paste0(
        "Each bar = mean of ", nyrs, " years (1950\u20132025).  ",
        "Annual totals: P = ", ann_P, " mm;  PET(PM) = ", ann_PET, " mm;  ",
        "aridity index (PET/P) = ", ai, "."),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 11, face = "bold", hjust = 0),
        plot.subtitle = ggplot2::element_text(size = 8.5, colour = "grey35",
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 4))))
  
  for (ext in c("pdf", "png")) {
    out_seas <- file.path(BASIN_PLOT_DIR,
                          paste0("FigS_Seasonality_P_PET.", ext))
    tryCatch(
      ggplot2::ggsave(out_seas, fig_seas_4d,
                      width  = 7.2,
                      height = 7.0,
                      units  = "in",
                      dpi    = if (ext == "png") 300 else "print",
                      device = ext),
      error = function(e) cat(sprintf("  \u26a0 %s: %s\n", ext, e$message)))
    cat(sprintf("  \u2713 Saved: %s\n", basename(out_seas)))
  }
  
}, error = function(e) {
  if (!identical(conditionMessage(e), "pet_csv_missing_4d"))
    cat(sprintf("  \u26a0 Part 4d failed: %s\n", e$message))
})

################################################################################
# PART 5 – EXCEL SUMMARY (statistics + drought event counts + correlations)
################################################################################
cat("\n══════════════════════════════════════════════\n")
cat("PART 5: Excel summary\n")
cat("══════════════════════════════════════════════\n")

## Build one stats row for a given index × scale.
make_stats_row <- function(idx, sc) {
  df <- tryCatch(load_any_ts(idx, sc), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  df_rec   <- df[df$date >= as.Date("2020-01-01") &
                   df$date <= as.Date("2025-12-31"), ]
  ev_count <- tryCatch(
    nrow(detect_drought_events(df_rec)),
    error = function(e) NA_integer_)
  data.frame(
    Timescale                = sprintf("%02d", sc),
    Index                    = toupper(idx),
    Mean                     = round(mean(df$value,   na.rm = TRUE), 4),
    Median                   = round(median(df$value, na.rm = TRUE), 4),
    StdDev                   = round(sd(df$value,     na.rm = TRUE), 4),
    Min                      = round(min(df$value,    na.rm = TRUE), 4),
    Max                      = round(max(df$value,    na.rm = TRUE), 4),
    Drought_Events_2020_2025 = ev_count,
    stringsAsFactors = FALSE)
}

export_basin_excel <- function(output_file) {
  wb        <- openxlsx::createWorkbook()
  hdr_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  
  ## ── Sheet 1: Summary Statistics ─────────────────────────────────────────────
  stats_rows <- list()
  for (idx in c("spi", "spei")) {
    for (sc in SPI_SCALES) {
      r <- tryCatch(make_stats_row(idx, sc),
                    error = function(e) {
                      cat(sprintf("  ⚠ Stats %s-%02d: %s\n",
                                  toupper(idx), sc, e$message))
                      NULL
                    })
      if (!is.null(r)) stats_rows[[length(stats_rows) + 1L]] <- r
    }
  }
  r <- tryCatch(make_stats_row("swei", SWEI_SCALE),
                error = function(e) NULL)
  if (!is.null(r)) stats_rows[[length(stats_rows) + 1L]] <- r
  
  for (idx in c("mspi", "mspei")) {
    r <- tryCatch(make_stats_row(idx, MSPI_MSPEI_SCALE),
                  error = function(e) NULL)
    if (!is.null(r)) stats_rows[[length(stats_rows) + 1L]] <- r
  }
  
  stats_df <- do.call(rbind, stats_rows)
  openxlsx::addWorksheet(wb, "Summary_Statistics")
  openxlsx::writeData(wb, "Summary_Statistics", stats_df)
  openxlsx::addStyle(wb, "Summary_Statistics", hdr_style,
                     rows = 1, cols = seq_len(ncol(stats_df)))
  cat(sprintf("  ✓ Summary_Statistics sheet: %d rows\n", nrow(stats_df)))
  
  ## ── Sheet 2: Pairwise Correlations ──────────────────────────────────────────
  ## Automatically correlate every available index pair on their common date range.
  corr_df <- tryCatch({
    series_list <- list(
      SPI3  = load_any_ts("spi",   3),
      SPEI3 = load_any_ts("spei",  3),
      MSPI  = load_any_ts("mspi",  MSPI_MSPEI_SCALE),
      MSPEI = load_any_ts("mspei", MSPI_MSPEI_SCALE)
    )
    series_list <- Filter(Negate(is.null), series_list)
    nms <- names(series_list)
    
    ## Merge all on date (inner join to common period)
    mrg <- Reduce(
      function(a, b_idx) {
        b  <- series_list[[b_idx]][, c("date", "value")]
        names(b)[2] <- nms[b_idx]
        merge(a, b, by = "date")
      },
      seq_along(nms)[-1],
      init = {
        tmp        <- series_list[[nms[1]]][, c("date", "value")]
        names(tmp)[2] <- nms[1]
        tmp
      })
    
    ## All pairwise Pearson correlations
    pairs <- combn(nms, 2, simplify = FALSE)
    do.call(rbind, lapply(pairs, function(p) {
      cx <- mrg[[p[1]]]; cy <- mrg[[p[2]]]
      ok <- !is.na(cx) & !is.na(cy)
      data.frame(
        Comparison  = sprintf("%s vs %s", p[1], p[2]),
        Correlation = round(cor(cx[ok], cy[ok]), 4),
        N_months    = sum(ok),
        stringsAsFactors = FALSE)
    }))
  }, error = function(e) {
    cat("  ⚠ Correlation sheet skipped:", e$message, "\n")
    data.frame(Comparison  = character(),
               Correlation = numeric(),
               N_months    = integer(),
               stringsAsFactors = FALSE)
  })
  
  openxlsx::addWorksheet(wb, "Correlations")
  openxlsx::writeData(wb, "Correlations", corr_df)
  openxlsx::addStyle(wb, "Correlations", hdr_style,
                     rows = 1, cols = seq_len(ncol(corr_df)))
  cat(sprintf("  ✓ Correlations sheet: %d pairs\n", nrow(corr_df)))
  
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  cat(sprintf("✓ Excel summary saved: %s\n", basename(output_file)))
  invisible(wb)
}

excel_out <- file.path(BASIN_PLOT_DIR, "Drought_Summary.xlsx")
export_basin_excel(excel_out)

################################################################################
# FINAL SUMMARY
################################################################################
cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  w2_basin_timeseries.R  DONE                 ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")
cat("Outputs:\n")
cat("  Per-index flat CSVs   : ", normalizePath(BASIN_TS_DIR),   "\n")
cat("  Combined long CSV     : ", normalizePath(BASIN_PLOT_DIR),  "\n")
cat("  Time series PNGs      : ", normalizePath(BASIN_PLOT_DIR),  "\n")
cat("  MS figure SPEI (PM)   : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4_MS_SPEI_123_PM_timeseries.pdf")), "\n")
cat("  MS figure SPEI (Thw)  : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4_MS_SPEI_123_Thw_timeseries.pdf")), "\n")
cat("  MS figure SPEI (both) : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4_MS_SPEI_123_overlay_timeseries.pdf")), "\n")
cat("  Spatial homogeneity (PDF): ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4c_spatial_homogeneity_SPEI3_PM.pdf")), "\n")
cat("  Spatial homogeneity (PNG): ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "Fig4c_spatial_homogeneity_SPEI3_PM.png")), "\n")
cat("  Excel summary         : ", normalizePath(excel_out),       "\n")
cat("  Seasonality P/PET(PM) : ",
    normalizePath(file.path(BASIN_PLOT_DIR,
                            "FigS_Seasonality_P_PET.pdf")), "\n")
cat("\n✓ All basin-averaged output tasks complete!\n")
cat("  Next: run w4_trends_visualization.R for spatial figures.\n\n")