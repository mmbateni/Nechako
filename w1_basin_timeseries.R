# ####################################################################################
# 01_basin_timeseries.R  ·  BASIN-AVERAGED TIME SERIES FOR ALL INDICES
# ####################################################################################
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "ncdf4", "ggplot2", "lubridate",
                      "future", "future.apply", "dplyr", "data.table", "patchwork"))
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
BASIN_PLOT_DIR <- file.path(WD_PATH, "temporal_drought", "basin_averaged_plots")   
dir.create(BASIN_PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
n_cores <- max(1L, parallel::detectCores() - 1L)
future::plan(future::multisession, workers = n_cores)
cat(sprintf("✓ Parallel: %d workers\n", n_cores))
basin_shp <- terra::vect(BASIN_SHP)
basin_proj <- terra::project(basin_shp, EQUAL_AREA_CRS)
cat(sprintf("✓ Basin: %.1f km²\n", terra::expanse(basin_proj, unit = "km")))
# ── Extraction Function ──────────────────────────────────────────────
# find_swei_seasonal_files takes 2 args; wrap to unified 3-arg interface
find_swei_nc_wrapper <- function(data_dir, index_type, scale) {
  find_swei_seasonal_files(data_dir, scale)
}
extract_basin_mean <- function(data_dir, index_type, scale, find_files_fn) {
  label <- sprintf("%s%d", toupper(index_type), scale)
  cache_path <- file.path(CACHE_DIR, paste0(label, "_basin_timeseries.rds"))
  if (file.exists(cache_path)) {
    cat(sprintf("  [cache] %s\n", label))
    return(readRDS(cache_path))
  }
  files <- find_files_fn(data_dir, index_type, scale)
  if (!length(files)) return(NULL)
  geom <- precompute_basin_geometry(BASIN_SHP, files[1], EQUAL_AREA_CRS)
  all_records <- vector("list", length(files))
  for (fi in seq_along(files)) {
    f <- files[fi]
    n_lyrs <- terra::nlyr(terra::rast(f))
    dates <- extract_dates_from_nc(f, n_lyrs)
    stats_list <- future.apply::future_lapply(
      seq_len(n_lyrs), compute_layer_statistics,
      raster_file = f, basin_geom = geom,
      threshold = DROUGHT_THRESHOLD, future.seed = TRUE)
    
    all_records[[fi]] <- data.frame(
      Date = dates,
      Index = label,
      Mean_Value = sapply(stats_list, `[[`, "mean_value"),
      Drought_Fraction = sapply(stats_list, `[[`, "drought_fraction"),
      stringsAsFactors = FALSE)
  }
  result <- dplyr::bind_rows(all_records)
  result <- result[order(result$Date), ]
  saveRDS(result, cache_path)
  cat(sprintf("  ✓ %s: %d time steps\n", label, nrow(result)))
  result
}
# ── Process All Indices ──────────────────────────────────────────────
all_data <- list()
# SPI
cat("\n── SPI ──\n")
for (sc in SPI_SCALES) {
  ts <- tryCatch(
    extract_basin_mean(SPI_SEAS_DIR, "spi", sc, find_seasonal_nc_files),
    error = function(e) { cat(sprintf("  ❌ SPI%d: %s\n", sc, e$message)); NULL })
  if (!is.null(ts)) all_data[[sprintf("SPI%d", sc)]] <- ts
}
# SPEI
cat("\n── SPEI ──\n")
for (sc in SPEI_SCALES) {
  ts <- tryCatch(
    extract_basin_mean(SPEI_SEAS_DIR, "spei", sc, find_seasonal_nc_files),
    error = function(e) { cat(sprintf("  ❌ SPEI%d: %s\n", sc, e$message)); NULL })
  if (!is.null(ts)) all_data[[sprintf("SPEI%d", sc)]] <- ts
}
# SWEI
cat("\n── SWEI ──\n")
ts <- tryCatch(
  extract_basin_mean(SWEI_SEAS_DIR, "swei", SWEI_SCALE, find_swei_nc_wrapper),
  error = function(e) { cat(sprintf("  ❌ SWEI%d: %s\n", SWEI_SCALE, e$message)); NULL })
if (!is.null(ts)) all_data[[sprintf("SWEI%d", SWEI_SCALE)]] <- ts
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
future::plan(future::sequential)
cat("\n✓ 01_basin_timeseries.R complete\n")