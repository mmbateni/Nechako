# ####################################################################################
# 01_basin_timeseries.R  ·  BASIN-AVERAGED TIME SERIES FOR ALL INDICES
# ####################################################################################
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "ggplot2", "lubridate",
                      "dplyr", "data.table", "patchwork"))
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
BASIN_PLOT_DIR <- file.path(WD_PATH, "temporal_drought", "basin_averaged_plots")   
dir.create(BASIN_PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
basin_shp  <- terra::vect(BASIN_SHP)
basin_proj <- terra::project(basin_shp, EQUAL_AREA_CRS)
cat(sprintf("✓ Basin: %.1f km²\n", terra::expanse(basin_proj, unit = "km")))

# ── Basin-averaged time series loader ───────────────────────────────
# Uses pre-computed *_basin_averaged_by_month.csv files.
# These indices were computed FROM basin-averaged climate inputs — the
# correct approach for basin-level analysis.  Spatially averaging
# per-pixel standardised indices would be methodologically incorrect.
load_basin_index <- function(index_type, scale) {
  label <- sprintf("%s%d", toupper(index_type), scale)
  df    <- tryCatch(
    load_basin_avg_csv(index_type, scale),
    error = function(e) { cat(sprintf("  ❌ %s: %s\n", label, e$message)); NULL })
  if (is.null(df) || !nrow(df)) {
    cat(sprintf("  ⚠ No data for %s\n", label)); return(NULL)
  }
  # Reshape to the same structure expected by downstream code:
  # Date, Index, Mean_Value, Drought_Fraction
  data.frame(
    Date             = df$date,
    Index            = label,
    Mean_Value       = df$value,
    Drought_Fraction = as.numeric(df$value <= DROUGHT_THRESHOLD),
    stringsAsFactors = FALSE)
}

# ── Process All Indices ──────────────────────────────────────────────
cat("\n── SPI (basin-averaged by month CSVs) ──\n")
all_data <- list()
for (sc in SPI_SCALES) {
  ts <- load_basin_index("spi", sc)
  if (!is.null(ts)) {
    all_data[[sprintf("SPI%d", sc)]] <- ts
    cat(sprintf("  ✓ SPI%d: %d time steps\n", sc, nrow(ts)))
  }
}

cat("\n── SPEI (basin-averaged by month CSVs) ──\n")
for (sc in SPEI_SCALES) {
  ts <- load_basin_index("spei", sc)
  if (!is.null(ts)) {
    all_data[[sprintf("SPEI%d", sc)]] <- ts
    cat(sprintf("  ✓ SPEI%d: %d time steps\n", sc, nrow(ts)))
  }
}

cat("\n── SWEI (basin-averaged by month CSV) ──\n")
ts <- load_basin_index("swei", SWEI_SCALE)
if (!is.null(ts)) {
  all_data[[sprintf("SWEI%d", SWEI_SCALE)]] <- ts
  cat(sprintf("  ✓ SWEI%d: %d time steps\n", SWEI_SCALE, nrow(ts)))
}

all_basin_data <- dplyr::bind_rows(all_data)
# ── Save CSV ─────────────────────────────────────────────────────────
if (nrow(all_basin_data) > 0) {
  write.csv(all_basin_data, file.path(BASIN_PLOT_DIR, "basin_mean_all_indices.csv"),
            row.names = FALSE)
  # ── Plot: Full Period ──────────────────────────────────────────────
  for (idx in unique(all_basin_data$Index)) {
    sub <- all_basin_data[all_basin_data$Index == idx, ]
    ev <- detect_drought_events(data.frame(date = sub$Date, value = sub$Mean_Value))
    p  <- ggplot2::ggplot(sub, ggplot2::aes(Date, Mean_Value)) +
      drought_band_layers() +
      ggplot2::geom_line(color = index_colours[tolower(substr(idx, 1, 4))],
                         linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      shared_ts_theme(12) +
      ggplot2::labs(title = sprintf("%s Basin-Averaged (1950–2025)", idx),
                    subtitle = sprintf("%d drought events detected", nrow(ev)),
                    x = "Year", y = "Index Value")
    
    ggplot2::ggsave(file.path(BASIN_PLOT_DIR, sprintf("%s_full_period.png", idx)),
                    p, width = 12, height = 6, dpi = 150)
  }
  # ── Excel Summary ──────────────────────────────────────────────────
  load_ts <- function(idx, sc) {
    label <- sprintf("%s%d", toupper(idx), sc)
    sub <- all_basin_data[all_basin_data$Index == label, ]
    data.frame(date = sub$Date, value = sub$Mean_Value)
  }
  all_scales <- c(SPI_SCALES, SPEI_SCALES, SWEI_SCALE)
  all_indices <- c(rep("spi", length(SPI_SCALES)),
                   rep("spei", length(SPEI_SCALES)),
                   rep("swei", 1))
  export_summary_excel(timescales = all_scales,
                       index_types = all_indices,
                       ts_loader_fn = load_ts,
                       output_file = file.path(BASIN_PLOT_DIR, "Drought_Summary.xlsx"))
}
cat("\n✓ 01_basin_timeseries.R complete\n")