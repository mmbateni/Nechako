####################################################################################
# w3_point_timeseries.R  ·  POINT-BASED ANALYSIS — DROUGHT INDICES (SPI / SPEI / SWEI)
# ─────────────────────────────────────────────────────────────────────────────────
# Extracts time series at SPECIFIED (or randomly-selected) grid points within the
# Nechako Basin for SPI, SPEI (multiple scales), and SWEI, and performs:
#
#   • Data-gap diagnostics (missing %, end-of-series gap, post-2020 completeness)
#   • Multi-panel time series plots  (SPI+SPEI: 1&12-month and 3&6-month 2×2 grids;
#                                    SWEI:    single full-width panel)
#   • 4-panel trend analysis plots with LR + Sen's slope fit lines
#   • Mann-Kendall (two variants) + Sen's slope + Linear regression statistics
#   • Trend summary CSV
#
#   NOTE: Temperature trend analysis at the same specific points is handled in
#   4pr_pet_trends.r (processing) and 4pr_pet_trends_visualization.r (plots).
#   Use the same SPECIFIC_PTS coordinates in both scripts.
#
# ── Mann-Kendall p-value columns ─────────────────────────────────────────────────
#   MK_pvalue_vc   (RECOMMENDED)  Hamed & Rao 1998 variance-corrected p-value.
#   MK_pvalue_std  ⚠ ANTI-CONSERVATIVE — DO NOT USE FOR INFERENCE.
#   Sens_pvalue    Normal-approximation p-value from trend::sens.slope().
#
# OUTPUT FILES PER POINT (coordinate-stamped filenames):
#   point_{id}_{lat}N_{lon}W_ts_1_12mo.png      SPI+SPEI 1 & 12-month TS
#   point_{id}_{lat}N_{lon}W_ts_3_6mo.png       SPI+SPEI 3 &  6-month TS
#   point_{id}_{lat}N_{lon}W_ts_swei.png        SWEI time series
#   point_{id}_{lat}N_{lon}W_trends.png         SPI+SPEI trend panels
#   point_{id}_{lat}N_{lon}W_trend_swei.png     SWEI trend panel
#   trend_statistics_summary.csv                All drought-index trend stats
####################################################################################

source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("Kendall", "trend", "modifiedmk", "lubridate", "patchwork",
                      "dplyr", "ggplot2", "terra", "tidyr"))
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

POINT_PLOT_DIR <- file.path(WD_PATH, "temporal_drought", "point_timeseries_plots")
dir.create(POINT_PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

####################################################################################
# ── CONFIG ───────────────────────────────────────────────────────────────────────
####################################################################################

# ── Point selection ───────────────────────────────────────────────────────────────
# Set USE_SPECIFIC = TRUE to always use exactly SPECIFIC_PTS (WGS84 lon/lat).
# Set USE_SPECIFIC = FALSE to randomly sample N_RANDOM grid points from the basin.
USE_SPECIFIC <- TRUE
SPECIFIC_PTS <- data.frame(
  x = c(-124.5, -123.0, -125.0),   # Longitude (negative = West)
  y = c(  54.0,   53.8,   54.2)    # Latitude (North)
)
N_RANDOM <- 5    # used only when USE_SPECIFIC = FALSE

# ── Trend analysis ────────────────────────────────────────────────────────────────
ALPHA          <- 0.05   # significance threshold
TREND_START_YR <- NULL   # NULL = full series
TREND_END_YR   <- NULL


####################################################################################
# ── BASIN & TEMPLATE RASTER ──────────────────────────────────────────────────────
####################################################################################
cat("Loading basin & template raster...\n")
basin_shp        <- load_basin_vect(BASIN_SHP)
spi_sample_files <- find_seasonal_nc_files(SPI_SEAS_DIR, "spi", SPI_SCALES[1])
if (!length(spi_sample_files))
  stop("No seasonal SPI files found to determine template CRS")
template_raster <- terra::rast(spi_sample_files[1])
if (!terra::same.crs(basin_shp, template_raster)) {
  cat("  Reprojecting basin to raster CRS\n")
  basin_shp <- terra::project(basin_shp, terra::crs(template_raster))
}

####################################################################################
# ── POINT SELECTION ──────────────────────────────────────────────────────────────
####################################################################################
get_basin_grid_pts <- function() {
  m  <- terra::rast(template_raster); terra::values(m) <- 1
  m  <- terra::mask(m, basin_shp)
  df <- as.data.frame(m, xy = TRUE, na.rm = TRUE)
  colnames(df)[1:2] <- c("x", "y")
  cat(sprintf("  %d grid points inside basin\n", nrow(df)))
  df[, c("x", "y")]
}

select_points <- function() {
  if (!USE_SPECIFIC) {
    # ── Random selection ───────────────────────────────────────────────────────
    all_pts <- get_basin_grid_pts()
    n       <- min(N_RANDOM, nrow(all_pts))
    set.seed(42)
    cat(sprintf("  Randomly selecting %d of %d points\n", n, nrow(all_pts)))
    return(all_pts[sample(seq_len(nrow(all_pts)), n), ])
  }
  
  # ── Specific point selection ───────────────────────────────────────────────
  cat(sprintf("  Creating %d specified points from SPECIFIC_PTS (WGS84)\n",
              nrow(SPECIFIC_PTS)))
  
  # Step 1: build SpatVector in WGS84
  pts_wgs84 <- terra::vect(SPECIFIC_PTS, geom = c("x", "y"), crs = "EPSG:4326")
  if (nrow(pts_wgs84) == 0) stop("Failed to create geometries from SPECIFIC_PTS")
  cat("  WGS84 coordinates:\n"); print(terra::crds(pts_wgs84))
  
  # Step 2: project to raster CRS (BC Albers / metres)
  pts_proj <- terra::project(pts_wgs84, terra::crs(template_raster))
  if (nrow(pts_proj) == 0) stop("Failed to project points to raster CRS")
  cat("  Projected coordinates (raster CRS, metres):\n"); print(terra::crds(pts_proj))
  
  # Step 3: basin intersection check
  inside_pts <- tryCatch(terra::intersect(pts_proj, basin_shp),
                         error = function(e) NULL)
  n_in  <- if (!is.null(inside_pts)) nrow(inside_pts) else 0L
  n_tot <- nrow(pts_proj)
  
  if (n_in == 0L) {
    # No points inside basin — use all specified points with a warning
    cat("  ⚠ WARNING: No specified points intersect the basin boundary.\n")
    cat("    Basin extent:\n"); print(terra::ext(basin_shp))
    cat(sprintf("  Proceeding with all %d specified point(s) [basin check skipped].\n",
                n_tot))
    out <- terra::crds(pts_proj)
    colnames(out) <- c("x", "y")
    return(as.data.frame(out))
  }
  
  if (n_in < n_tot) {
    # Partial match — warn about dropped points
    cat(sprintf(
      "  ⚠ WARNING: %d of %d specified point(s) fall outside the basin boundary\n",
      n_tot - n_in, n_tot))
    cat("    Outside points are excluded from analysis.\n")
    cat("    Check SPECIFIC_PTS coordinates if this is unexpected.\n")
  }
  
  out <- terra::crds(inside_pts)
  colnames(out) <- c("x", "y")
  cat(sprintf("  ✓ Using %d specified point(s) inside basin\n", nrow(out)))
  as.data.frame(out)
}

####################################################################################
# ── DROUGHT INDEX TIME SERIES EXTRACTION ─────────────────────────────────────────
####################################################################################
#' Extract time series for one drought index at one point.
#'
#' BUGFIX: find_swei_seasonal_files() does not accept a `scale` argument; it is
#' called without scale to avoid "argument 'scale' is missing" errors.
extract_ts <- function(data_dir, index_type, scale, x, y, label) {
  cat(sprintf("    Extracting %s ...\n", label))
  files <- if (index_type == "swei")
    find_swei_seasonal_files(data_dir)        # no scale argument for SWEI
  else
    find_seasonal_nc_files(data_dir, index_type, scale)
  
  if (!length(files)) {
    cat(sprintf("      ⚠ No files found for %s\n", label)); return(NULL)
  }
  r_sample    <- terra::rast(files[1])
  pt          <- terra::vect(data.frame(x = x, y = y), geom = c("x", "y"),
                             crs = terra::crs(r_sample))
  all_records <- vector("list", length(files))
  for (fi in seq_along(files)) {
    f       <- files[fi]
    r_stack <- terra::rast(f)
    dates   <- extract_dates_from_nc(f, terra::nlyr(r_stack))
    vals    <- as.numeric(terra::extract(r_stack, pt)[1, -1])
    n       <- min(length(dates), length(vals))
    all_records[[fi]] <- data.frame(Date  = dates[1:n], Value = vals[1:n],
                                    Index = label, stringsAsFactors = FALSE)
  }
  result <- dplyr::bind_rows(all_records)
  result <- result[order(result$Date), ]
  rownames(result) <- NULL
  result
}

####################################################################################
# ── DATA GAP DIAGNOSTICS ─────────────────────────────────────────────────────────
####################################################################################
diagnose_data_gaps <- function(ts_data, label) {
  n_total   <- nrow(ts_data)
  n_missing <- sum(is.na(ts_data$Value))
  pct_miss  <- round(100 * n_missing / n_total, 1)
  cat(sprintf("    %s: %d/%d missing (%.1f%%)\n", label, n_missing, n_total, pct_miss))
  rec <- dplyr::filter(ts_data, lubridate::year(Date) >= 2020)
  if (nrow(rec)) {
    pct_rec <- 100 * sum(is.na(rec$Value)) / nrow(rec)
    if (pct_rec > 50)
      cat(sprintf("    ⚠ %s: %.0f%% missing after 2020 – check reference period!\n",
                  label, pct_rec))
  }
  valid <- dplyr::filter(ts_data, !is.na(Value))
  if (nrow(valid)) {
    last_v   <- max(valid$Date)
    gap_mths <- as.numeric(lubridate::interval(last_v, max(ts_data$Date)) /
                             lubridate::dmonths(1))
    if (gap_mths > 6)
      cat(sprintf("    ⚠ %s: gap of %.0f months at end of series\n", label, gap_mths))
  }
  invisible(list(n_missing = n_missing, pct_missing = pct_miss))
}

####################################################################################
# ── DROUGHT INDEX TREND ANALYSIS ─────────────────────────────────────────────────
####################################################################################
#' Compute LR, Sen's slope, and two Mann-Kendall variants for one series.
#' Returns a list: $stats (data.frame), $data_clean, $data_original.
#'
#' MK COLUMN NAMING
#'  MK_pvalue_vc  / MK_significant_vc   ← RECOMMENDED  (H-R98 autocorr-corrected)
#'  MK_pvalue_std / MK_significant_std  ← DIAGNOSTIC ONLY — anti-conservative
calc_trend_stats <- function(ts_data, label,
                             start_yr = TREND_START_YR,
                             end_yr   = TREND_END_YR) {
  if (!is.null(start_yr)) ts_data <- dplyr::filter(ts_data, lubridate::year(Date) >= start_yr)
  if (!is.null(end_yr))   ts_data <- dplyr::filter(ts_data, lubridate::year(Date) <= end_yr)
  clean <- dplyr::filter(ts_data, !is.na(Value))
  if (nrow(clean) < 10) {
    cat(sprintf("    Insufficient data for %s\n", label)); return(NULL)
  }
  pct_complete <- nrow(clean) / nrow(ts_data) * 100
  if (pct_complete < 70)
    cat(sprintf("    ⚠ %s: only %.0f%% complete\n", label, pct_complete))
  
  clean$t_idx <- seq_len(nrow(clean))
  lm_m        <- stats::lm(Value ~ t_idx, data = clean)
  lm_sum      <- summary(lm_m)
  sens_r      <- tryCatch(trend::sens.slope(clean$Value), error = function(e) NULL)
  mk_r        <- tryCatch(Kendall::MannKendall(clean$Value), error = function(e) NULL)
  if (is.null(sens_r) || is.null(mk_r)) return(NULL)
  
  mk_vc      <- tryCatch(modifiedmk::mmkh(clean$Value), error = function(e) NULL)
  pvalue_vc  <- if (!is.null(mk_vc)) as.numeric(mk_vc["new P-value"]) else NA_real_
  sig_vc     <- !is.na(pvalue_vc) && pvalue_vc < ALPHA
  pvalue_std <- as.numeric(mk_r$sl)
  sig_std    <- pvalue_std < ALPHA
  
  clean$LR_fit   <- stats::predict(lm_m)
  clean$Sens_fit <- stats::median(clean$Value) +
    sens_r$estimates * (clean$t_idx - stats::median(clean$t_idx))
  
  list(
    stats = data.frame(
      Index              = label,
      Variable_type      = "drought_index",
      N_obs              = nrow(clean),
      N_total            = nrow(ts_data),
      Completeness_pct   = round(pct_complete, 1),
      Period_start       = format(min(clean$Date), "%Y-%m"),
      Period_end         = format(max(clean$Date), "%Y-%m"),
      LR_slope           = stats::coef(lm_m)[2],
      LR_pvalue          = lm_sum$coefficients[2, 4],
      LR_significant     = lm_sum$coefficients[2, 4] < ALPHA,
      LR_rsquared        = lm_sum$r.squared,
      Sens_slope         = as.numeric(sens_r$estimates),
      Sens_pvalue        = as.numeric(sens_r$p.value),
      Sens_significant   = as.numeric(sens_r$p.value) < ALPHA,
      MK_tau             = as.numeric(mk_r$tau),
      MK_pvalue_vc       = pvalue_vc,
      MK_significant_vc  = sig_vc,
      MK_pvalue_std      = pvalue_std,
      MK_significant_std = sig_std,
      stringsAsFactors   = FALSE),
    data_clean    = clean,
    data_original = ts_data)
}

####################################################################################
# ── DROUGHT INDEX PLOT BUILDERS ──────────────────────────────────────────────────
####################################################################################
# Reusable drought-band time series panel
make_ts_panel <- function(ts_data, label, x_label = "") {
  ts_data$Value <- ifelse(is.infinite(ts_data$Value), NA, ts_data$Value)
  ylm <- calc_dynamic_ylim(ts_data$Value)
  ggplot2::ggplot(ts_data, ggplot2::aes(Date, Value)) +
    drought_band_layers() +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.9, na.rm = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    ggplot2::coord_cartesian(ylim = ylm) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x      = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = label, x = x_label, y = "Index Value")
}

# Single trend panel — VC p-value is primary; std p-value shown for comparison only
make_trend_panel <- function(trend_res, label) {
  if (is.null(trend_res))
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = paste(label, "– insufficient data"),
                               size = 4, color = "red") +
             ggplot2::theme_void())
  orig  <- trend_res$data_original
  clean <- trend_res$data_clean
  st    <- trend_res$stats
  note  <- if (st$Completeness_pct < 80)
    sprintf("\n⚠ %.0f%% complete", st$Completeness_pct) else ""
  ann <- sprintf(
    paste0(
      "LR:  β=%.4f (p=%.3f)%s\n",
      "Sen: β=%.4f (p=%.3f)%s\n",
      "MK:  τ=%.3f\n",
      "  p_vc=%.3f%s  ← use for inference\n",
      "  p_std=%.3f%s  [uncorrected, do not use]%s"),
    st$LR_slope,      st$LR_pvalue,     ifelse(st$LR_significant,  "*", ""),
    st$Sens_slope,    st$Sens_pvalue,   ifelse(st$Sens_significant, "*", ""),
    st$MK_tau,
    st$MK_pvalue_vc,  ifelse(st$MK_significant_vc,  "*", ""),
    st$MK_pvalue_std, ifelse(st$MK_significant_std, "*", ""),
    note)
  ylm <- calc_dynamic_ylim(orig$Value)
  ggplot2::ggplot(orig, ggplot2::aes(Date, Value)) +
    ggplot2::geom_line(color = "steelblue", alpha = 0.6, na.rm = FALSE) +
    ggplot2::geom_point(color = "steelblue", size = 0.4, alpha = 0.4, na.rm = TRUE) +
    ggplot2::geom_line(data = clean, ggplot2::aes(y = LR_fit),
                       color = "red", linewidth = 1, alpha = 0.8) +
    ggplot2::geom_line(data = clean, ggplot2::aes(y = Sens_fit),
                       color = "darkgreen", linewidth = 1, linetype = "dashed", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    ggplot2::annotate("text",
                      x = min(orig$Date, na.rm = TRUE), y = max(orig$Value, na.rm = TRUE),
                      label = ann, hjust = 0, vjust = 1, size = 2.5) +
    ggplot2::coord_cartesian(ylim = ylm) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(title = label, x = "Date", y = "Value")
}

# 4-panel time series: two SPI scales + two SPEI scales
make_four_panel_ts <- function(all_ts, sc_a, sc_b, point, pt_id, outfile) {
  labels <- c(sprintf("SPI%d",  sc_a), sprintf("SPI%d",  sc_b),
              sprintf("SPEI%d", sc_a), sprintf("SPEI%d", sc_b))
  plots <- lapply(labels, function(lbl) {
    sub <- dplyr::filter(all_ts, Index == lbl)
    if (!nrow(sub))
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5,
                                 label = paste(lbl, "not available"),
                                 size = 5, color = "red") +
               ggplot2::theme_void())
    make_ts_panel(sub, lbl, if (lbl == labels[4]) "Date" else "")
  })
  names(plots) <- labels
  combined <-
    (plots[[labels[1]]] | plots[[labels[2]]]) /
    (plots[[labels[3]]] | plots[[labels[4]]]) +
    patchwork::plot_annotation(
      title = sprintf("Point #%d  (%.4f°N, %.4f°W)  —  %d & %d-Month Scales",
                      pt_id, point$y, abs(point$x), sc_a, sc_b),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5))) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(outfile, combined, width = 12, height = 10, dpi = 150)
  cat(sprintf("    Saved: %s\n", basename(outfile)))
}

# Single SWEI time series panel (full width)
make_swei_ts_plot <- function(all_ts, point, pt_id, outfile) {
  lbl <- sprintf("SWEI%d", SWEI_SCALE)
  sub <- dplyr::filter(all_ts, Index == lbl)
  if (!nrow(sub)) {
    cat(sprintf("    ⚠ %s not available – skipping SWEI TS plot\n", lbl))
    return(invisible(NULL))
  }
  p <- make_ts_panel(sub, lbl, "Date") +
    ggplot2::labs(title = sprintf("%s  —  Point #%d  (%.4f°N, %.4f°W)",
                                  lbl, pt_id, point$y, abs(point$x)))
  ggplot2::ggsave(outfile, p, width = 12, height = 5, dpi = 150)
  cat(sprintf("    Saved: %s\n", basename(outfile)))
}

# 4-panel drought index trend plots: SPI1, SPI12, SPEI1, SPEI12
make_trend_plot <- function(trend_list, labels, point, pt_id, outfile) {
  plots    <- lapply(labels, function(lbl) make_trend_panel(trend_list[[lbl]], lbl))
  combined <-
    ((plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])) +
    patchwork::plot_annotation(
      title    = sprintf("Point #%d  (%.4f°N, %.4f°W)  — Trend Analysis",
                         pt_id, point$y, abs(point$x)),
      subtitle = paste0("Red = linear regression  |  Green dashed = Sen's slope\n",
                        "MK: p_vc = H-R98 corrected [use for inference]  |  ",
                        "p_std = uncorrected [diagnostic only, do not cite]"),
      theme    = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 14, face = "bold",    hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size =  9, color = "gray30", hjust = 0.5)))
  ggplot2::ggsave(outfile, combined, width = 16, height = 8, dpi = 150)
  cat(sprintf("    Saved: %s\n", basename(outfile)))
}

# Single SWEI trend panel
make_swei_trend_plot <- function(trend_list, point, pt_id, outfile) {
  lbl <- sprintf("SWEI%d", SWEI_SCALE)
  tr  <- trend_list[[lbl]]
  if (is.null(tr)) {
    cat(sprintf("    ⚠ No trend data for %s – skipping\n", lbl))
    return(invisible(NULL))
  }
  p <- make_trend_panel(tr, lbl) +
    ggplot2::labs(
      title    = sprintf("%s  Trend Analysis  —  Point #%d  (%.4f°N, %.4f°W)",
                         lbl, pt_id, point$y, abs(point$x)),
      subtitle = paste0("Red = linear regression  |  Green dashed = Sen's slope\n",
                        "p_vc = H-R98 corrected [use for inference]  |  ",
                        "p_std = uncorrected [diagnostic only]"))
  ggplot2::ggsave(outfile, p, width = 10, height = 6, dpi = 150)
  cat(sprintf("    Saved: %s\n", basename(outfile)))
}



####################################################################################
# ── MAIN PROCESSING LOOP ─────────────────────────────────────────────────────────
####################################################################################
cat("\n══════════════════════════════════════\n")
cat("MAIN: Processing selected points\n")
cat("══════════════════════════════════════\n")

# Helper: convert projected CRS (metres) back to WGS84 for display
proj_to_latlon <- function(x, y, from_crs) {
  pt_v   <- terra::vect(data.frame(x = x, y = y), geom = c("x", "y"), crs = from_crs)
  pt_geo <- terra::project(pt_v, "EPSG:4326")
  coords <- terra::crds(pt_geo)
  list(lat = coords[1, 2], lon = coords[1, 1])
}

raster_crs      <- terra::crs(template_raster)
selected_pts    <- select_points()
selected_pts$point_id <- seq_len(nrow(selected_pts))
all_trend_stats <- list()

# ── Per-point loop ────────────────────────────────────────────────────────────────
for (i in seq_len(nrow(selected_pts))) {
  pt  <- selected_pts[i, ]
  geo <- tryCatch(proj_to_latlon(pt$x, pt$y, raster_crs),
                  error = function(e) list(lat = pt$y, lon = pt$x))
  lat_disp <- geo$lat;  lon_disp <- geo$lon
  cat(sprintf("\n── Point %d / %d  (%.4f°N, %.4f°W) ──\n",
              i, nrow(selected_pts), lat_disp, abs(lon_disp)))
  
  # ── Drought indices ──────────────────────────────────────────────────────────
  ts_all <- list()
  
  for (sc in SPI_SCALES) {
    lbl <- sprintf("SPI%d", sc)
    ts_all[[lbl]] <- tryCatch(
      extract_ts(SPI_SEAS_DIR, "spi", sc, pt$x, pt$y, lbl),
      error = function(e) { cat("  ❌", lbl, ":", e$message, "\n"); NULL })
  }
  for (sc in SPEI_SCALES) {
    lbl <- sprintf("SPEI%d", sc)
    ts_all[[lbl]] <- tryCatch(
      extract_ts(SPEI_SEAS_DIR, "spei", sc, pt$x, pt$y, lbl),
      error = function(e) { cat("  ❌", lbl, ":", e$message, "\n"); NULL })
  }
  lbl <- sprintf("SWEI%d", SWEI_SCALE)
  ts_all[[lbl]] <- tryCatch(
    extract_ts(SWEI_SEAS_DIR, "swei", SWEI_SCALE, pt$x, pt$y, lbl),
    error = function(e) { cat("  ❌", lbl, ":", e$message, "\n"); NULL })
  
  # Data gap diagnostics
  cat("  Data gap diagnostics:\n")
  for (lbl in names(ts_all))
    if (!is.null(ts_all[[lbl]])) diagnose_data_gaps(ts_all[[lbl]], lbl)
  
  all_ts_df <- dplyr::bind_rows(ts_all)
  
  # Drought index trend analysis
  cat("  Trend analysis:\n")
  trends <- lapply(names(ts_all), function(lbl) {
    if (is.null(ts_all[[lbl]])) return(NULL)
    calc_trend_stats(ts_all[[lbl]], lbl)
  })
  names(trends) <- names(ts_all)
  
  pt_stats <- dplyr::bind_rows(lapply(trends, function(tr) if (!is.null(tr)) tr$stats))
  if (nrow(pt_stats)) {
    pt_stats$Point_ID <- pt$point_id
    pt_stats$Point_X  <- pt$x
    pt_stats$Point_Y  <- pt$y
    all_trend_stats[[paste0("drought_", i)]] <- pt_stats
  }
  
  # ── File name base ───────────────────────────────────────────────────────────
  fn_base    <- sprintf("point_%d_%.4fN_%.4fW", pt$point_id, lat_disp, abs(lon_disp))
  f_1_12     <- file.path(POINT_PLOT_DIR, paste0(fn_base, "_ts_1_12mo.png"))
  f_3_6      <- file.path(POINT_PLOT_DIR, paste0(fn_base, "_ts_3_6mo.png"))
  f_swei     <- file.path(POINT_PLOT_DIR, paste0(fn_base, "_ts_swei.png"))
  f_trend    <- file.path(POINT_PLOT_DIR, paste0(fn_base, "_trends.png"))
  f_trend_sw <- file.path(POINT_PLOT_DIR, paste0(fn_base, "_trend_swei.png"))
  
  pt_geo <- list(y = lat_disp, x = lon_disp, point_id = pt$point_id)
  
  # ── Drought time series plots ────────────────────────────────────────────────
  tryCatch(make_four_panel_ts(all_ts_df, 1, 12, pt_geo, pt$point_id, f_1_12),
           error = function(e) cat("  ⚠ 1&12 plot:", e$message, "\n"))
  tryCatch(make_four_panel_ts(all_ts_df, 3,  6, pt_geo, pt$point_id, f_3_6),
           error = function(e) cat("  ⚠ 3&6 plot:", e$message, "\n"))
  tryCatch(make_swei_ts_plot(all_ts_df, pt_geo, pt$point_id, f_swei),
           error = function(e) cat("  ⚠ SWEI TS plot:", e$message, "\n"))
  
  # ── Drought trend plots ──────────────────────────────────────────────────────
  trend_labels <- c("SPI1", "SPI12", "SPEI1", "SPEI12")
  avail        <- trends[intersect(trend_labels, names(trends))]
  if (length(avail) == 4)
    tryCatch(make_trend_plot(avail, trend_labels, pt_geo, pt$point_id, f_trend),
             error = function(e) cat("  ⚠ trend plot:", e$message, "\n"))
  tryCatch(make_swei_trend_plot(trends, pt_geo, pt$point_id, f_trend_sw),
           error = function(e) cat("  ⚠ SWEI trend plot:", e$message, "\n"))
  
  rm(ts_all, all_ts_df, trends); invisible(gc())
}

####################################################################################
# ── SAVE TREND STATISTICS CSV + CONSOLE SUMMARY ──────────────────────────────────
####################################################################################
if (length(all_trend_stats)) {
  stats_df <- dplyr::bind_rows(all_trend_stats)
  csv_out  <- file.path(POINT_PLOT_DIR, "trend_statistics_summary.csv")
  utils::write.csv(stats_df, csv_out, row.names = FALSE)
  cat(sprintf("\n✓ Trend statistics saved: %s  (%d rows)\n",
              basename(csv_out), nrow(stats_df)))
  
  cat("\nTrend Summary  (counts of significant points):\n")
  cat("  NOTE: MK_significant_vc is recommended for inference.\n")
  cat("        MK_significant_std is uncorrected and anti-conservative.\n\n")
  
  smry <- stats_df |>
    dplyr::group_by(Variable_type, Index) |>
    dplyr::summarise(
      n               = dplyr::n(),
      pct_complete    = mean(Completeness_pct,  na.rm = TRUE),
      sig_LR          = sum(LR_significant,      na.rm = TRUE),
      sig_Sens        = sum(Sens_significant,     na.rm = TRUE),
      sig_MK_vc       = sum(MK_significant_vc,   na.rm = TRUE),
      sig_MK_std_DIAG = sum(MK_significant_std,  na.rm = TRUE),
      .groups         = "drop")
  print(as.data.frame(smry))
  cat("\n  Column sig_MK_std_DIAG is for diagnostic comparison only.\n")
  cat("  Use sig_MK_vc for any significance claims.\n")
}

cat("\n╔══════════════════════════════════════╗\n")
cat("║  w3_point_timeseries.R  DONE         ║\n")
cat("╚══════════════════════════════════════╝\n\n")
cat("⚠ NOTE: Gaps in drought plots reflect genuine missing/NA values.\n")
cat("  If gaps appear after 2020, check the reference period and\n")
cat("  whether input data extends to 2024.\n\n")
cat("MK p-value reminder:\n")
cat("  MK_pvalue_vc  / MK_significant_vc   → H-R98 corrected  [USE FOR INFERENCE]\n")
cat("  MK_pvalue_std / MK_significant_std  → uncorrected       [DIAGNOSTIC ONLY]\n\n")
cat("NOTE: Temperature trends for the same specific points are analysed in\n")
cat("  4pr_pet_trends.r (processing) and 4pr_pet_trends_visualization.r (plots).\n\n")
cat("Outputs: ", normalizePath(POINT_PLOT_DIR), "\n")