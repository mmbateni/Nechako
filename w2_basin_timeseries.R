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
cat("  Excel summary         : ", normalizePath(excel_out),       "\n")
cat("\n✓ All basin-averaged output tasks complete!\n")
cat("  Next: run w4_trends_visualization.R for spatial figures.\n\n")