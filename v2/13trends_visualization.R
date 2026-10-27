# 13trends_visualization.R  ·  SPATIAL, BASIN, HOMOGENEITY AND POINT DROUGHT DIAGNOSTICS
# ─────────────────────────────────────────────────────────────────────────────────────────
# COMPATIBLE WITH: DROUGHT_ANALYSIS_utils.R, 11trend_test_ALL.R
# INDICES: SPI, SPEI-PM, SPEI-ERA5, SRI (1-36), SWEI-3, SSPI-1, SDI k1-k4
#   PART 1  merge Script 11 pixel results -> master CSVs (full series + calendar month)
#   PART 2  diagnostics: TFPW sign concordance, grid alignment, statistics table, QA
#   PART 3  spatial figures: Fig1 tau + direction, Fig2 events, Fig3 temporal patterns,
#           Fig4 timescale comparison (incl. lag-1 rho vs theory), Fig5 method comparison,
#           Fig6 spatial z-score, Fig7 calendar-month Sen maps, Fig8 calendar-month summary,
#           Fig9 drought-month-fraction trends
#   PART 4  Sen-slope figures: Figure X (full series), calendar-month manuscript figure,
#           snow (calendar month), runoff (SRI + SDI)
#   PART 5  basin and regional: p-value comparison, calendar-month heatmap and series,
#           seasonal (SPI-3 Feb/May/Aug/Nov), regional heatmap and time series
#   PART 6  homogeneity: break classification, by region, pixel shares and maps,
#           annual series with step means, sub-period trends
#   PART 7  specific points
#   PART 8  pixel-to-basin Spearman correlation
#   PART 9  regional domain map + final validation
# EXECUTION ORDER: 11trend_test_ALL.R -> 12basin_timeseries.R -> 13trends_visualization.R

rm(list = ls())

# ---- Working directory (framework rule) -------------------------------------
# Run every script from the project root: set the working directory there
# BEFORE running (setwd(), or open the RStudio project). All paths below are
# derived from getwd(); no path is hard-coded.
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")

source(UTILS_FILE)
utils_load_packages(c("terra", "data.table", "ggplot2", "sf", "viridis", "cowplot", "patchwork",
                      "scales", "RColorBrewer", "ggspatial", "gridExtra", "grid", "lubridate"))

fig_dir <- file.path(TREND_DIR, "figures_temporal_diagnostics")
FIGD <- c(diag = "diagnostics", cm = "calendar_month", freq = "annual_drought_frequency",
          basin = "basin_regional", homog = "homogeneity", points = "points")
FIGD <- setNames(file.path(fig_dir, FIGD), names(FIGD))
for (d in c(fig_dir, FIGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ── Native ERA5-Land grid and basin mask (authoritative geometry; never resampled) ──
native_precip_file <- file.path(WD_PATH, "monthly_data_direct", "total_precipitation_monthly.nc")
if (!file.exists(native_precip_file)) stop("Native ERA5-Land precipitation NetCDF not found: ", native_precip_file)
analysis_template <- terra::rast(native_precip_file)[[1]]
if (is.na(terra::crs(analysis_template)) || !nzchar(terra::crs(analysis_template))) stop("Native grid has no CRS.")
common_template <- analysis_template
basin_vect <- nechako_subbasin_union(load_nechako_subbasins())   # 7-subbasin union, as in Script 3
if (!inherits(basin_vect, "SpatVector")) basin_vect <- terra::vect(basin_vect)
if (!terra::same.crs(basin_vect, analysis_template)) basin_vect <- terra::project(basin_vect, terra::crs(analysis_template))
basin_sf <- sf::st_as_sf(basin_vect); basin_sf_global <<- basin_sf
cov <- as.numeric(terra::values(terra::rasterize(basin_vect, analysis_template, cover = TRUE), mat = FALSE))
cov[!is.finite(cov)] <- 0
basin_mask_raster <- terra::setValues(analysis_template, as.numeric(cov > 0))
basin_mask_cells  <- as.logical(terra::values(basin_mask_raster, mat = FALSE))
basin_pixels <- sum(basin_mask_cells, na.rm = TRUE)
if (basin_pixels == 0L) stop("Native basin mask contains zero cells.")
cat(sprintf("[GRID] %d x %d | basin cells %d\n", terra::nrow(analysis_template), terra::ncol(analysis_template), basin_pixels))

# ── Constants ──────────────────────────────────────────────────────────────────────
DPI <- 300
TIMESCALES_VIZ <- as.integer(TIMESCALES_SPATIAL)
SWEI_SCALE <- 3L; SSPI_SCALE <- 1L
if (!exists("SWEI_SEAS_DIR"))      SWEI_SEAS_DIR      <- file.path(WD_PATH, "swei_results_seasonal")
if (!exists("SPEI_ERA5_SEAS_DIR")) SPEI_ERA5_SEAS_DIR <- file.path(WD_PATH, "spei_results_seasonal_era5land")
SSPI_SEAS_DIR <- file.path(WD_PATH, "sspi_results_seasonal")
SRI_SEAS_DIR  <- file.path(WD_PATH, "sri_results_seasonal", "runoff")
SDI_PERIODS <- c(3L, 6L, 9L, 12L)
SDI_LABELS  <- c(`3` = "k1 Oct–Dec", `6` = "k2 Oct–Mar", `9` = "k3 Oct–Jun", `12` = "k4 Oct–Sep")
MS_FULL_METHOD     <- "vc"     # Figure X / runoff figure: "vc" (recommended for autocorrelated series) or "tfpw"
CALMONTH_METHOD    <- "tfpw"   # calendar-month maps (k <= 12 nearly independent; same as Script 9)
CALMONTH_TS_SCALES <- c(3L, 12L)
POINT_TS_SCALE     <- 12L
S6_CFG <- if (file.exists(file.path(TREND_DIR, "script6_analysis_config.rds"))) readRDS(file.path(TREND_DIR, "script6_analysis_config.rds")) else NULL
STEP_YEAR <- if (!is.null(S6_CFG$step_year)) as.integer(S6_CFG$step_year) else 1979L
BASIN_RID <- if (exists("NECHAKO_REGION_BASIN_ID")) NECHAKO_REGION_BASIN_ID else "Nechako_Basin"
HOMOG_DIR <- file.path(TREND_DIR, "homogeneity")

INDEX_LABELS <- c(spi = "SPI", spei = "SPEI-PM", spei_era5 = "SPEI-ERA5", swei = "SWEI", sspi = "SSPI",
                  sri = "SRI", sdi = "SDI", mspi = "MSPI", mspei = "MSPEI")
INDEX_ORDER <- c("spi", "spei", "spei_era5", "swei", "sspi", "sri", "sdi")
index_display <- function(x) { x <- as.character(x); out <- toupper(x); h <- x %in% names(INDEX_LABELS); out[h] <- INDEX_LABELS[x[h]]; unname(out) }
for (nm in names(c(spi = 1, spei = 1, swei = 1, sspi = 1, spei_era5 = 1, sri = 1, sdi = 1)))
  if (!nm %in% names(index_colours))
    index_colours <- c(index_colours, setNames(c(spi = "#377eb8", spei = "#e41a1c", swei = "#1E90FF", sspi = "#984EA3",
                                                 spei_era5 = "#d62728", sri = "#a65628", sdi = "#f781bf")[nm], nm))
.icol <- function(idx) if (idx %in% names(index_colours)) unname(index_colours[[idx]]) else "grey30"

BRBG9 <- c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E")
PAL_BRBG <- grDevices::colorRampPalette(BRBG9)(101)      # negative (drier) = brown, positive (wetter) = blue-green
DIR_COLS <- c("#8c510a", "#dfc27d", "#f5f5f5", "#80cdc1", "#01665e")
DIR_LABS <- c("Decrease / drier (FDR, both tests)", "Decrease (nominal, both tests)", "Not significant",
              "Increase (nominal, both tests)", "Increase / wetter (FDR, both tests)")

.read  <- function(p) if (file.exists(p)) data.table::fread(p) else NULL
.fmt_p <- function(p) ifelse(!is.finite(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
# Dual-test significance (framework rule for all drought-index maps, full series
# AND calendar-month series): BH-FDR p < 0.05 under BOTH VC-HR98 and TFPW =
# "both" (filled stippling); exactly one test = "one test only" (open circles).
# MS_FULL_METHOD / CALMONTH_METHOD only choose which Sen slope / tau is coloured.
.dual_sig <- function(d) {
  vc <- if ("p_fdr_vc"   %in% names(d)) is.finite(d$p_fdr_vc)   & d$p_fdr_vc   < 0.05 else rep(FALSE, nrow(d))
  tf <- if ("p_fdr_tfpw" %in% names(d)) is.finite(d$p_fdr_tfpw) & d$p_fdr_tfpw < 0.05 else rep(FALSE, nrow(d))
  list(both = vc & tf, one = xor(vc, tf), vc = vc, tfpw = tf)
}
.stipple_dual <- function(s, cex = 0.35) {
  g <- .dual_sig(s)
  if (any(g$one))  graphics::points(s$lon[g$one],  s$lat[g$one],  pch = 1,  cex = cex, lwd = 0.5)
  if (any(g$both)) graphics::points(s$lon[g$both], s$lat[g$both], pch = 20, cex = cex)
  invisible(g)
}
.has_dual_cols <- function(d, label) {
  ok <- all(c("p_fdr_vc", "p_fdr_tfpw") %in% names(d))
  if (!ok) cat("  ⚠", label, ": p_fdr_vc and/or p_fdr_tfpw missing; dual-test stippling uses what is available\n")
  ok
}
.mcol  <- function(method) if (method == "vc")
  list(tau = "tau_vc", p = "p_value_vc", fdr = "p_fdr_vc", sl = "sl_vc", lab = "VC-HR98") else
    list(tau = "tau_tfpw", p = "p_value_tfpw", fdr = "p_fdr_tfpw", sl = "sl_tfpw", lab = "TFPW")
.save_png <- function(file, plot, w, h) { ggplot2::ggsave(file, plot, width = w, height = h, dpi = DPI, bg = "white"); cat("  ✓", basename(file), "\n") }
.step <- function(label, expr) tryCatch(expr, error = function(e) cat("  ❌", label, ":", conditionMessage(e), "\n"))
make_raster <- function(dt, vcol, tpl = analysis_template) create_raster_from_points(dt, tpl, vcol, basin_sf)
.row_lab <- function(idx, sc) sprintf("%s-%d", index_display(idx), sc)
.row_levels <- function(d) {
  u <- unique(d[, .(index_type, timescale)]); u[, o := match(index_type, INDEX_ORDER)]
  data.table::setorder(u, o, timescale); .row_lab(u$index_type, u$timescale)
}

################################################################################
# PART 1 – MERGE SCRIPT 6 RESULTS
################################################################################
cat("\n══ PART 1: master diagnostics files ══\n")
load_clip_results <- function(index_type, timescale, kind = "results", quiet = FALSE) {
  f <- file.path(TREND_DIR, sprintf("%s_%02d_%s.csv", index_type, as.integer(timescale), kind))
  if (!file.exists(f)) { if (!quiet) cat("  ⚠ Missing:", basename(f), "\n"); return(NULL) }
  dt <- data.table::fread(f)
  if ("x" %in% names(dt) && !"lon" %in% names(dt)) data.table::setnames(dt, "x", "lon")
  if ("y" %in% names(dt) && !"lat" %in% names(dt)) data.table::setnames(dt, "y", "lat")
  if (!all(c("lon", "lat") %in% names(dt))) { cat("  ⚠ No coordinates in", basename(f), "\n"); return(NULL) }
  dt <- dt[is.finite(lon) & is.finite(lat)]
  uc <- unique(dt[, .(lon, lat)])
  rel <- terra::relate(terra::vect(as.data.frame(uc), geom = c("lon", "lat"), crs = terra::crs(analysis_template)),
                       basin_vect, relation = "intersects", pairs = FALSE)
  inside <- if (is.matrix(rel)) apply(rel, 1L, any) else as.logical(rel)
  dt <- dt[uc[inside], on = .(lon, lat), nomatch = 0L]
  if (!nrow(dt)) return(NULL)
  it <- index_type; tsx <- as.integer(timescale)
  dt[, `:=`(index_type = it, timescale = tsx)]
  cat(sprintf("  %-9s %-16s %7d rows\n", .row_lab(it, tsx), kind, nrow(dt)))
  dt
}
SPEC_FULL <- data.table::rbindlist(list(
  data.table::CJ(index_type = c("spi", "spei", "spei_era5", "sri"), timescale = TIMESCALES_VIZ, sorted = FALSE),
  data.table::data.table(index_type = c("swei", "sspi"), timescale = c(SWEI_SCALE, SSPI_SCALE)),
  data.table::data.table(index_type = "sdi", timescale = SDI_PERIODS)))
full_list <- lapply(seq_len(nrow(SPEC_FULL)), function(i) load_clip_results(SPEC_FULL$index_type[i], SPEC_FULL$timescale[i]))
plot_data <- data.table::rbindlist(full_list, fill = TRUE)
if (!nrow(plot_data)) stop("No Script 11 results loaded. Run 11trend_test_ALL.R first.")
spec_cm <- SPEC_FULL[index_type != "sdi"]
cm_list <- lapply(seq_len(nrow(spec_cm)), function(i) load_clip_results(spec_cm$index_type[i], spec_cm$timescale[i], "calmonth_results", quiet = TRUE))
cm_data <- data.table::rbindlist(cm_list, fill = TRUE)
if (!nrow(cm_data)) { cm_data <- NULL; cat("  ⚠ No calendar-month results: re-run 11trend_test_ALL.R.\n") }
if (!is.null(cm_data)) for (cn in c("p_fdr_vc", "p_fdr_tfpw")) if (!cn %in% names(cm_data)) cm_data[, (cn) := NA_real_]
for (cn in c("sl_vc", "p_fdr_vc", "p_fdr_tfpw", "qa_full", "rho1_detrended", "dfreq_sl_vc", "dfreq_p_fdr_vc"))
  if (!cn %in% names(plot_data)) plot_data[, (cn) := NA]
if (all(is.na(plot_data$sl_vc))) plot_data[, sl_vc := sl_tfpw]   # Yue (2002) slope = Sen slope of the original series
master_csv <- file.path(TREND_DIR, "all_temporal_diagnostics_results.csv")
data.table::fwrite(plot_data, master_csv)
if (!is.null(cm_data)) data.table::fwrite(cm_data, file.path(TREND_DIR, "all_calendar_month_results.csv"))
cat(sprintf("✓ Master CSV: %d rows | calendar-month: %d rows\n", nrow(plot_data), if (is.null(cm_data)) 0L else nrow(cm_data)))

active_indices <- function(scale) {
  base <- intersect(c("spi", "spei", "spei_era5"), unique(plot_data[timescale == scale, index_type]))
  if (scale == SWEI_SCALE && any(plot_data$index_type == "swei")) base <- c(base, "swei")
  if (scale == SSPI_SCALE && any(plot_data$index_type == "sspi")) base <- c(base, "sspi")
  if (any(plot_data$index_type == "sri" & plot_data$timescale == scale)) base <- c(base, "sri")
  base
}

################################################################################
# PART 2 – DIAGNOSTICS
################################################################################
cat("\n══ PART 2: diagnostics ══\n")
diag_sign <- function() {
  f <- function(d, lab) d[is.finite(tau_tfpw) & is.finite(sl_tfpw) & sl_tfpw != 0, .(
    analysis = lab, n = .N, agree_all = mean(sign(tau_tfpw) == sign(sl_tfpw)),
    n_sig = sum(p_value_tfpw < 0.05, na.rm = TRUE),
    agree_sig = { s <- which(p_value_tfpw < 0.05); if (length(s)) mean(sign(tau_tfpw[s]) == sign(sl_tfpw[s])) else NA_real_ }),
    by = .(index_type, timescale)]
  out <- data.table::rbindlist(list(f(plot_data, "full_series"), if (!is.null(cm_data)) f(cm_data, "calendar_month")), fill = TRUE)
  data.table::fwrite(out, file.path(FIGD[["diag"]], "tfpw_sign_concordance.csv"))
  print(as.data.frame(out[analysis == "full_series"]), digits = 3)
  bad <- out[is.finite(agree_sig) & agree_sig < 0.99]
  if (nrow(bad)) {
    warning("TFPW tau and Sen slope disagree in sign among SIGNIFICANT pixels for: ",
            paste(unique(.row_lab(bad$index_type, bad$timescale)), collapse = ", "), call. = FALSE)
  } else cat("  ✓ Sign disagreements occur only in non-significant pixels (expected for TFPW).\n")
}

diag_grid <- function() {
  rs <- terra::res(analysis_template)
  out <- plot_data[, {
    cells <- terra::cellFromXY(analysis_template, cbind(lon, lat)); ok <- !is.na(cells)
    ctr <- terra::xyFromCell(analysis_template, cells[ok])
    list(n = .N, outside = sum(!ok),
         off_centre_pct = if (any(ok)) 100 * mean(abs(lon[ok] - ctr[, 1]) / rs[1] > 0.05 | abs(lat[ok] - ctr[, 2]) / rs[2] > 0.05) else NA_real_,
         duplicated = sum(duplicated(cells[ok])))
  }, by = .(index_type, timescale)]
  data.table::fwrite(out, file.path(FIGD[["diag"]], "grid_alignment.csv"))
  if (any(out$outside > 0 | out$off_centre_pct > 1 | out$duplicated > 0, na.rm = TRUE)) {
    print(as.data.frame(out))
    warning("Pixel coordinates do not align with the native ERA5-Land grid for some index x scale.", call. = FALSE)
  } else cat(sprintf("  ✓ Grid alignment: all %d index x scale combinations sit on native cell centres.\n", nrow(out)))
}

diag_stats <- function() {
  pct <- function(p) if (any(is.finite(p))) 100 * mean(p[is.finite(p)] < 0.05) else NA_real_
  st <- plot_data[, list(
    n_pixels = .N, n_tested = sum(is.finite(p_value_vc)),
    filtered = { q <- table(qa_full[qa_full != "ok"]); paste(names(q), q, sep = ":", collapse = "; ") },
    pct_sig_vc = pct(p_value_vc), pct_fdr_vc = pct(p_fdr_vc), pct_sig_tfpw = pct(p_value_tfpw), pct_fdr_tfpw = pct(p_fdr_tfpw),
    pct_fdr_both = pct(pmax(p_fdr_vc, p_fdr_tfpw)),
    median_tau_vc = stats::median(tau_vc, na.rm = TRUE), median_sl_vc = stats::median(sl_vc, na.rm = TRUE),
    mean_rho1 = mean(rho1_detrended, na.rm = TRUE)), by = .(index_type, timescale)]
  if (!is.null(cm_data)) {
    m <- .mcol(CALMONTH_METHOD)
    cs <- cm_data[, list(pct = { pb <- pmax(p_fdr_vc, p_fdr_tfpw)
    if (any(is.finite(pb))) 100 * mean(pb[is.finite(pb)] < 0.05) else NA_real_ }),
    by = .(index_type, timescale, month)]
    cs <- cs[, list(calmonth_median_pct_fdr = stats::median(pct, na.rm = TRUE),
                    calmonth_max_pct_fdr = if (any(is.finite(pct))) max(pct, na.rm = TRUE) else NA_real_,
                    calmonth_month_of_max = if (any(is.finite(pct))) month[which.max(pct)] else NA_integer_),
             by = .(index_type, timescale)]
    st <- merge(st, cs, by = c("index_type", "timescale"), all.x = TRUE)
  }
  data.table::fwrite(st, file.path(FIGD[["diag"]], "statistics_summary_drought.csv"))
  print(as.data.frame(st[, .(index_type, timescale, n_tested, pct_fdr_vc, pct_fdr_tfpw, pct_fdr_both, mean_rho1)]), digits = 3)
}

.step("sign concordance", diag_sign())
.step("grid alignment",   diag_grid())
.step("statistics table", diag_stats())
qa_rc <- .read(file.path(TREND_DIR, "qa", "root_cause_filter_reasons.csv"))
if (!is.null(qa_rc)) { cat("  Script 11 QA filters (series not tested):\n"); print(as.data.frame(qa_rc)) }

################################################################################
# PART 3 – SPATIAL FIGURES
################################################################################
cat("\n══ PART 3: spatial figures ══\n")
.no_panel <- function(msg) { graphics::plot.new(); graphics::text(0.5, 0.5, msg, cex = 0.9, col = "grey35") }
.dir_cat <- function(tau, p, pf) {
  fdr <- is.finite(pf) & pf < 0.05; nom <- is.finite(p) & p < 0.05
  data.table::fcase(fdr & tau < 0, 1, fdr & tau > 0, 5, nom & tau < 0, 2, nom & tau > 0, 4, default = 3)
}

# Fig 1: tau (VC-HR98) + dual-test FDR stippling (VC-HR98 AND TFPW), and a 5-class
# direction map whose significance classes require BOTH tests. mfcol keeps each index in one column.
create_trend_maps <- function(scale, outfile) {
  ind <- active_indices(scale); n <- length(ind); if (!n) return(invisible(NULL))
  if (!safe_pdf(outfile, max(8, 3.6 * n), 8.5)) return(invisible(NULL))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfcol = c(2, n), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  for (idx in ind) {
    lab <- .row_lab(idx, scale)
    sub <- data.table::copy(plot_data[index_type == idx & timescale == scale & is.finite(tau_vc)])
    if (!nrow(sub)) { .no_panel(paste0(lab, "\nno tested pixels (QA-filtered)")); .no_panel(lab); next }
    plot_raster_clean(make_raster(sub, "tau_vc"), paste0(lab, ": Kendall's τ (VC-HR98)"), c(-0.4, 0.4), PAL_BRBG, legend_title = "τ")
    .stipple_dual(sub, cex = 0.4)
    graphics::legend("bottomleft", legend = c("BH-FDR: VC and TFPW", "BH-FDR: one test only"),
                     pch = c(20, 1), bty = "n", cex = 0.65)
    sub[, dir_cat := .dir_cat(tau_vc, pmax(p_value_vc, p_value_tfpw), pmax(p_fdr_vc, p_fdr_tfpw))]
    plot_raster_clean(make_raster(sub, "dir_cat"), paste0(lab, ": direction"), NULL, DIR_COLS,
                      breaks = seq(0.5, 5.5, 1), categorical = TRUE, legend = FALSE)
    graphics::legend("bottomright", legend = DIR_LABS, fill = DIR_COLS, bty = "n", cex = 0.6)
  }
  graphics::mtext(sprintf("Full monthly series — VC-HR98 τ; significance requires BOTH VC-HR98 and TFPW (BH-FDR) (%d-month)", scale),
                  outer = TRUE, cex = 1.0, font = 2)
  cat("  ✓", basename(outfile), "\n")
}

# Fig 2: event characteristics; darker = more events / longer / more intense (intensity is positive).
create_event_maps <- function(scale, outfile) {
  ind <- active_indices(scale); n <- length(ind); if (!n) return(invisible(NULL))
  if (!safe_pdf(outfile, max(8, 3.6 * n), 11)) return(invisible(NULL))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfcol = c(3, n), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  metrics <- list(n_events = c("Number of drought events", "YlOrRd"), mean_duration = c("Mean duration (months)", "viridis"),
                  max_intensity = c("Max intensity (threshold − index)", "Reds"))
  for (idx in ind) {
    sub <- plot_data[index_type == idx & timescale == scale]
    for (mn in names(metrics)) {
      ttl <- paste0(.row_lab(idx, scale), ": ", metrics[[mn]][1])
      if (!nrow(sub) || !mn %in% names(sub) || !any(is.finite(sub[[mn]]))) { .no_panel(paste0(ttl, "\nno data")); next }
      s <- sub[is.finite(get(mn))]
      plot_raster_clean(make_raster(s, mn), ttl, c(0, max(s[[mn]])), grDevices::hcl.colors(101, metrics[[mn]][2], rev = TRUE))
    }
  }
  graphics::mtext(sprintf("Drought event characteristics (%d-month)", scale), outer = TRUE, cex = 1.1, font = 2)
  cat("  ✓", basename(outfile), "\n")
}

# Fig 3: clustering (BH-FDR, as labelled by Script 11), regime-shift decade, spectral peaks.
create_temporal_maps <- function(scale, outfile) {
  ind <- active_indices(scale); n <- length(ind); if (!n) return(invisible(NULL))
  miss <- setdiff(c("filtered_runs", "p_value_runs", "clustering", "regime_shift_year", "n_spectral_peaks"), names(plot_data))
  if (length(miss)) { cat("  ⚠ Fig 3 skipped; missing columns:", paste(miss, collapse = ", "), "\n"); return(invisible(NULL)) }
  if (!safe_pdf(outfile, 12, max(4, 3.6 * n))) return(invisible(NULL))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(n, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  for (idx in ind) {
    lab <- .row_lab(idx, scale)
    sub <- data.table::copy(plot_data[index_type == idx & timescale == scale])
    sub[, clus_num := data.table::fcase(filtered_runs | is.na(p_value_runs), NA_real_,
                                        clustering == "clustered", -1, clustering == "dispersed", 1, default = 0)]
    s <- sub[is.finite(clus_num)]
    if (nrow(s)) {
      plot_raster_clean(make_raster(s, "clus_num"), paste0(lab, ": temporal clustering"), NULL,
                        c("#e41a1c", "#fee090", "#4daf4a"), breaks = c(-1.5, -0.5, 0.5, 1.5), categorical = TRUE, legend = FALSE)
      graphics::legend("bottomright", legend = c("Clustered (FDR)", "Not significant", "Dispersed (FDR)"),
                       fill = c("#e41a1c", "#fee090", "#4daf4a"), bty = "n", cex = 0.75)
    } else .no_panel(paste0(lab, "\nruns test not defined"))
    rs <- sub[is.finite(regime_shift_year)]
    if (nrow(rs)) {
      rs[, shift_dec := as.numeric(cut(regime_shift_year, breaks = seq(1950, 2030, 10), labels = 1:8, include.lowest = TRUE))]
      plot_raster_clean(make_raster(rs, "shift_dec"), paste0(lab, ": regime-shift decade"), c(1, 8),
                        grDevices::hcl.colors(8, "Spectral", rev = TRUE), legend_title = "1 = 1950s")
    } else .no_panel(paste0(lab, "\nno regime shifts"))
    sp <- sub[is.finite(n_spectral_peaks)]
    if (nrow(sp)) plot_raster_clean(make_raster(sp, "n_spectral_peaks"), paste0(lab, ": spectral peaks (≥ 24 months)"),
                                    c(0, max(1, sp$n_spectral_peaks)), grDevices::hcl.colors(101, "viridis"))
    else .no_panel(paste0(lab, "\nno spectral result"))
  }
  graphics::mtext(sprintf("Temporal patterns (%d-month)", scale), outer = TRUE, cex = 1.1, font = 2)
  cat("  ✓", basename(outfile), "\n")
}

# Fig 4: timescale comparison incl. lag-1 autocorrelation against the (k-1)/k expectation.
create_timescale_comparison <- function(outfile) {
  pct <- function(p) if (any(is.finite(p))) 100 * mean(p[is.finite(p)] < 0.05) else NA_real_
  smry <- plot_data[index_type != "sdi", list(
    Pct_Sig_Nom = pct(pmax(p_value_vc, p_value_tfpw)), Pct_Sig_FDR = pct(pmax(p_fdr_vc, p_fdr_tfpw)), Median_Tau = stats::median(tau_vc, na.rm = TRUE),
    Mean_rho1 = mean(rho1_detrended, na.rm = TRUE), Mean_Events = mean(n_events, na.rm = TRUE),
    Mean_Dur = mean(mean_duration, na.rm = TRUE)), by = .(index_type, timescale)]
  base <- function(col, ylab, title) ggplot(smry, aes(timescale, .data[[col]], colour = index_type, group = index_type)) +
    geom_line(linewidth = 1, na.rm = TRUE) + geom_point(size = 2.6, na.rm = TRUE) +
    scale_colour_manual(values = index_colours, labels = index_display) + scale_x_continuous(breaks = TIMESCALES_VIZ) +
    theme_bw(base_size = 11) + labs(title = title, x = "Timescale (months)", y = ylab, colour = "Index")
  sig <- data.table::melt(smry, id.vars = c("index_type", "timescale"), measure.vars = c("Pct_Sig_Nom", "Pct_Sig_FDR"),
                          variable.name = "Correction", value.name = "Pct")
  sig[, Correction := ifelse(Correction == "Pct_Sig_Nom", "Nominal (p < 0.05)", "BH-FDR (p < 0.05)")]
  pA <- ggplot(sig[is.finite(Pct)], aes(timescale, Pct, colour = index_type, linetype = Correction, group = interaction(index_type, Correction))) +
    geom_line(linewidth = 1) + geom_point(size = 2.2) + scale_colour_manual(values = index_colours, labels = index_display) +
    scale_x_continuous(breaks = TIMESCALES_VIZ) + theme_bw(base_size = 11) +
    labs(title = "A) % significant pixels (VC-HR98)", subtitle = "Gap between lines = false-discovery inflation",
         x = "Timescale (months)", y = "% of tested pixels", colour = "Index", linetype = NULL)
  theo <- data.frame(timescale = seq(1, 36, 0.5)); theo$rho <- pmax(0, (theo$timescale - 1) / theo$timescale)
  pC <- base("Mean_rho1", "lag-1 ρ", "C) Mean lag-1 autocorrelation (detrended monthly series)") +
    geom_line(data = theo, aes(timescale, rho), inherit.aes = FALSE, linetype = "dashed", colour = "grey40") +
    labs(subtitle = "Dashed = (k−1)/k expected from k-month overlap alone")
  p <- (pA | base("Median_Tau", "τ", "B) Median Kendall's τ (VC-HR98)")) /
    (pC | base("Mean_Events", "events", "D) Mean number of drought events")) /
    (base("Mean_Dur", "months", "E) Mean event duration") | patchwork::plot_spacer()) +
    patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
  .save_png(outfile, p, 14, 13)
}

# Fig 5 / 5b: VC-HR98 vs TFPW (full series; calendar-month series pooled over months).
.method_comparison <- function(d, title_sfx, outfile) {
  d <- data.table::copy(d[is.finite(tau_vc) & is.finite(tau_tfpw)])
  if (!nrow(d)) { cat("  ⚠ no rows for", basename(outfile), "\n"); return(invisible(NULL)) }
  d[, `:=`(panel = .row_lab(index_type, timescale), vc_sig = is.finite(p_fdr_vc) & p_fdr_vc < 0.05,
           tfpw_sig = is.finite(p_fdr_tfpw) & p_fdr_tfpw < 0.05)]
  p1 <- ggplot(d, aes(tau_vc, tau_tfpw, colour = index_type)) + geom_abline(linetype = "dashed", colour = "grey40") +
    geom_point(alpha = 0.35, size = 1.1) + facet_wrap(~panel) + scale_colour_manual(values = index_colours, guide = "none") +
    coord_equal() + theme_bw(base_size = 10) + labs(title = paste("A) τ: VC-HR98 vs TFPW —", title_sfx), x = "τ (VC-HR98)", y = "τ (TFPW)")
  ag <- d[, list(Both = sum(vc_sig & tfpw_sig), VC_only = sum(vc_sig & !tfpw_sig), TFPW_only = sum(!vc_sig & tfpw_sig),
                 Neither = sum(!vc_sig & !tfpw_sig)), by = panel]
  agl <- data.table::melt(ag, id.vars = "panel", variable.name = "Category", value.name = "N")
  agl[, pct := 100 * N / sum(N), by = panel]
  p2 <- ggplot(agl, aes(panel, pct, fill = Category)) + geom_col() + scale_fill_brewer(palette = "Set2") + coord_flip() +
    theme_bw(base_size = 10) + labs(title = "B) BH-FDR significance agreement", x = NULL, y = "% of tested pixels")
  .save_png(outfile, p1 | p2, 15, 7)
}
create_method_comparison <- function() {
  sel <- function(d) d[(index_type %in% c("spi", "spei", "spei_era5", "sri") & timescale == 12L) |
                         (index_type == "swei" & timescale == SWEI_SCALE) | (index_type == "sspi" & timescale == SSPI_SCALE) |
                         (index_type == "sdi" & timescale == 12L)]
  .method_comparison(sel(plot_data), "full series (12-month; snow native; SDI water year)", file.path(fig_dir, "Fig5_MethodComparison.png"))
  if (!is.null(cm_data)) .method_comparison(sel(cm_data), "calendar-month series, months pooled",
                                            file.path(fig_dir, "Fig5b_MethodComparison_CalendarMonth.png"))
}

# Fig 6: spatial pattern = z-score of the VC Sen slope relative to the basin mean.
create_spatial_pattern <- function(scale, outfile) {
  ind <- active_indices(scale); n <- length(ind); if (!n) return(invisible(NULL))
  if (!safe_pdf(outfile, max(8, 3.6 * n), 4.8)) return(invisible(NULL))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, n), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  for (idx in ind) {
    s <- data.table::copy(plot_data[index_type == idx & timescale == scale & is.finite(sl_vc)])
    sdv <- if (nrow(s) >= 3L) stats::sd(s$sl_vc) else NA_real_
    if (!is.finite(sdv) || sdv <= 0) { .no_panel(paste0(.row_lab(idx, scale), "\nno spatial variance")); next }
    s[, zval := (sl_vc - mean(sl_vc)) / sdv]
    plot_raster_clean(make_raster(s, "zval"), paste0(.row_lab(idx, scale), ": slope anomaly (z)"), c(-3, 3), PAL_BRBG, legend_title = "z")
  }
  graphics::mtext(sprintf("Spatial variability of the full-series Sen slope (z-score vs basin mean), %d-month", scale), outer = TRUE, cex = 1, font = 2)
  cat("  ✓", basename(outfile), "\n")
}

# Fig 7: 12-panel calendar-month Sen maps per index x scale (colour clipped at the 98th percentile of |slope|).
create_calmonth_maps <- function(idx, sc) {
  m <- .mcol(CALMONTH_METHOD)
  d <- cm_data[index_type == idx & timescale == sc & is.finite(get(m$sl))]
  if (!nrow(d)) return(invisible(NULL))
  z <- stats::quantile(abs(d[[m$sl]]), 0.98, names = FALSE)
  if (!is.finite(z) || z <= 0) z <- max(abs(d[[m$sl]]))
  if (!is.finite(z) || z <= 0) z <- 1
  out <- file.path(FIGD[["cm"]], sprintf("Fig7_CalendarMonth_SenMaps_%s_%02d.pdf", idx, sc))
  if (!safe_pdf(out, 14, 11)) return(invisible(NULL))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(3, 4), mar = c(2, 2, 3, 3), oma = c(0, 0, 3, 0))
  for (mm in 1:12) {
    s <- data.table::copy(d[month == mm])
    if (!nrow(s)) { .no_panel(paste0(month.name[mm], "\nfiltered (no valid years)")); next }
    s[, vclip := pmax(pmin(get(m$sl), z), -z)]
    plot_raster_clean(make_raster(s, "vclip"), month.name[mm], c(-z, z), PAL_BRBG, legend_title = "Sen/yr")
    g <- .stipple_dual(s, cex = 0.35)
    graphics::mtext(sprintf("FDR both %d / %d | one test %d", sum(g$both), nrow(s), sum(g$one)),
                    side = 3, line = 0.15, cex = 0.6, col = "grey30")
  }
  graphics::mtext(sprintf("%s — calendar-month trends (%s Sen's slope per year; filled = BH-FDR within month under BOTH VC-HR98 and TFPW, open = one test)",
                          .row_lab(idx, sc), m$lab), outer = TRUE, cex = 0.9, font = 2, line = 0.8)
  cat("  ✓", basename(out), "\n")
}

# Fig 8: calendar-month summary heatmap (net % FDR-significant pixels; median tau).
create_calmonth_summary <- function(outfile) {
  m <- .mcol(CALMONTH_METHOD)
  s <- cm_data[is.finite(get(m$p)), list(
    pct_inc = 100 * mean(is.finite(pmax(p_fdr_vc, p_fdr_tfpw)) & pmax(p_fdr_vc, p_fdr_tfpw) < 0.05 & get(m$tau) > 0),
    pct_dec = 100 * mean(is.finite(pmax(p_fdr_vc, p_fdr_tfpw)) & pmax(p_fdr_vc, p_fdr_tfpw) < 0.05 & get(m$tau) < 0),
    med_tau = stats::median(get(m$tau), na.rm = TRUE)), by = .(index_type, timescale, month)]
  lv <- .row_levels(s)
  s[, `:=`(rowlab = factor(.row_lab(index_type, timescale), levels = rev(lv)),
           mon = factor(month.abb[month], levels = month.abb), net = pct_inc - pct_dec)]
  pA <- ggplot(s, aes(mon, rowlab, fill = net)) + geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.0f", net)), size = 2.4) +
    scale_fill_gradientn(colours = PAL_BRBG, limits = c(-100, 100), oob = scales::squish, name = "% wetting −\n% drying") +
    theme_minimal(base_size = 10) + labs(title = "A) Net share of pixels BH-FDR-significant under BOTH VC-HR98 and TFPW", x = NULL, y = NULL)
  pB <- ggplot(s, aes(mon, rowlab, fill = med_tau)) + geom_tile(colour = "white") +
    scale_fill_gradientn(colours = PAL_BRBG, limits = c(-0.3, 0.3), oob = scales::squish, name = "median τ") +
    theme_minimal(base_size = 10) + theme(axis.text.y = element_blank()) + labs(title = "B) Median Kendall's τ", x = NULL, y = NULL)
  .save_png(outfile, (pA | pB) + patchwork::plot_annotation(
    title = sprintf("Calendar-month trends of all drought indices (%s)", m$lab)), 15, 0.28 * length(lv) + 2.5)
}

# Fig 9: trend in the annual fraction of drought months (percentage points per decade).
create_dfreq_maps <- function(scale, outfile) {
  ind <- active_indices(scale); n <- length(ind)
  sel <- plot_data[timescale == scale & index_type %in% ind & is.finite(dfreq_sl_vc)]
  if (!n || !nrow(sel)) return(invisible(NULL))
  if (!safe_pdf(outfile, max(8, 3.6 * n), 4.8)) return(invisible(NULL))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, n), mar = c(2, 2, 3, 3), oma = c(0, 0, 2.5, 0))
  z <- max(abs(1000 * sel$dfreq_sl_vc)); if (!is.finite(z) || z <= 0) z <- 1
  for (idx in ind) {
    s <- data.table::copy(sel[index_type == idx])
    if (!nrow(s)) { .no_panel(paste0(.row_lab(idx, scale), "\nno tested pixels")); next }
    s[, vdec := 1000 * dfreq_sl_vc]
    plot_raster_clean(make_raster(s, "vdec"), paste0(.row_lab(idx, scale), ": drought-month fraction"), c(-z, z), rev(PAL_BRBG),
                      legend_title = "pp/decade")
    sig <- s[is.finite(dfreq_p_fdr_vc) & dfreq_p_fdr_vc < 0.05]
    if (nrow(sig)) graphics::points(sig$lon, sig$lat, pch = 20, cex = 0.35)
  }
  graphics::mtext(sprintf("Trend in the annual fraction of months below %.1f (VC-HR98 Sen slope; brown = more drought), %d-month",
                          DROUGHT_ONSET, scale), outer = TRUE, cex = 0.95, font = 2)
  cat("  ✓", basename(outfile), "\n")
}

for (sc in TIMESCALES_VIZ) {
  sfx <- sprintf("%02dmo", sc)
  .step(sprintf("Fig 1 (%d)", sc), create_trend_maps(sc, file.path(fig_dir, sprintf("Fig1_TrendMaps_%s.pdf", sfx))))
  .step(sprintf("Fig 2 (%d)", sc), create_event_maps(sc, file.path(fig_dir, sprintf("Fig2_EventMaps_%s.pdf", sfx))))
  .step(sprintf("Fig 3 (%d)", sc), create_temporal_maps(sc, file.path(fig_dir, sprintf("Fig3_TemporalMaps_%s.pdf", sfx))))
  .step(sprintf("Fig 6 (%d)", sc), create_spatial_pattern(sc, file.path(fig_dir, sprintf("Fig6_SpatialPattern_%s.pdf", sfx))))
  .step(sprintf("Fig 9 (%d)", sc), create_dfreq_maps(sc, file.path(FIGD[["freq"]], sprintf("Fig9_DroughtMonthFraction_%s.pdf", sfx))))
}
.step("Fig 4", create_timescale_comparison(file.path(fig_dir, "Fig4_TimescaleComparison.png")))
.step("Fig 5", create_method_comparison())
if (!is.null(cm_data)) {
  combos <- unique(cm_data[, .(index_type, timescale)])
  for (i in seq_len(nrow(combos)))
    .step(sprintf("Fig 7 %s-%d", combos$index_type[i], combos$timescale[i]), create_calmonth_maps(combos$index_type[i], combos$timescale[i]))
  .step("Fig 8", create_calmonth_summary(file.path(FIGD[["cm"]], "Fig8_CalendarMonth_Summary.png")))
}

################################################################################
# PART 4 – SEN-SLOPE FIGURES (manuscript Figure X, calendar-month, snow, runoff)
################################################################################
cat("\n══ PART 4: Sen-slope figures ══\n")
# panels: data.table(index_type, timescale, month [NA = full series], title, unit)
draw_sen_figure <- function(panels, data, method, layout, title, subtitle, outpdf, outpng, required = FALSE) {
  m <- .mcol(method)
  if (is.null(data) || !all(c("index_type", "timescale", "lon", "lat", m$sl, m$fdr) %in% names(data))) {
    msg <- paste("Sen-slope figure", basename(outpdf), ": required columns missing")
    if (required) stop(msg); cat("  ⚠", msg, "\n"); return(invisible(FALSE))
  }
  .has_dual_cols(data, basename(outpdf))
  panels <- data.table::as.data.table(panels)
  if (!"month" %in% names(panels)) panels[, month := NA_integer_]
  pdat <- lapply(seq_len(nrow(panels)), function(p) {
    it <- panels$index_type[p]; tsx <- panels$timescale[p]; mo_p <- panels$month[p]
    s <- data[index_type == it & timescale == tsx]
    if (!is.na(mo_p) && "month" %in% names(s)) s <- s[month == mo_p]
    s[is.finite(get(m$sl)) & is.finite(lon) & is.finite(lat)]
  })
  n_ok <- vapply(pdat, nrow, integer(1))
  if (!any(n_ok > 0)) {
    msg <- paste("No valid data for", basename(outpdf)); if (required) stop(msg); cat("  ⚠", msg, "\n"); return(invisible(FALSE))
  }
  zl <- tapply(unlist(lapply(pdat, function(s) abs(s[[m$sl]]))), rep(panels$index_type, n_ok), max)
  zl[!is.finite(zl) | zl <= 0] <- 1
  draw <- function() {
    graphics::par(mfrow = layout, mar = c(2.2, 2.2, 3.0, 2.0), oma = c(0, 0, 3.4, 0), xaxs = "i", yaxs = "i")
    for (p in seq_len(nrow(panels))) {
      s <- pdat[[p]]; ttl <- panels$title[p]
      if (!nrow(s)) { .no_panel(paste0("No valid data\n", ttl)); next }
      zz <- zl[[panels$index_type[p]]]
      plot_raster_clean(make_raster(s, m$sl), ttl, c(-zz, zz), PAL_BRBG, legend_title = panels$unit[p])
      g <- .stipple_dual(s, cex = 0.32)
      graphics::mtext(sprintf("BH-FDR both: %d/%d (%.1f%%) | VC only %d | TFPW only %d",
                              sum(g$both), nrow(s), 100 * sum(g$both) / nrow(s),
                              sum(g$vc & !g$tfpw), sum(g$tfpw & !g$vc)),
                      side = 3, line = 0.15, cex = 0.5, col = "grey30")
    }
    graphics::mtext(paste0(title, " — Nechako River Basin (1950–2025)"), outer = TRUE, side = 3, line = 1.9, cex = 1.12, font = 2)
    graphics::mtext(subtitle, outer = TRUE, side = 3, line = 0.6, cex = 0.7, col = "grey30")
  }
  w_in <- 3.5 * layout[2]; h_in <- 3.3 * layout[1] + 0.7
  unlink(c(outpdf, outpng), force = TRUE)
  grDevices::pdf(outpdf, width = w_in, height = h_in, onefile = TRUE, useDingbats = FALSE); draw(); grDevices::dev.off()
  grDevices::png(outpng, width = w_in, height = h_in, units = "in", res = DPI, type = "cairo", bg = "white"); draw(); grDevices::dev.off()
  for (f in c(outpdf, outpng)) if (!file.exists(f) || !isTRUE(file.info(f)$size > 0)) stop("Figure not written: ", f)
  cat("  ✓", basename(outpdf), "/", basename(outpng), "\n")
  invisible(TRUE)
}

.panels_full <- function(idxs, scs) {
  p <- data.table::CJ(index_type = idxs, timescale = as.integer(scs), sorted = FALSE)
  p[, `:=`(month = NA_integer_,
           title = ifelse(index_type == "sdi", paste("SDI", SDI_LABELS[as.character(timescale)]),
                          sprintf("%s — %d months", index_display(index_type), timescale)),
           unit = ifelse(index_type == "sdi", "Sen/yr", "Sen/month"))]
  p
}
.sub_txt <- function(method, series) sprintf("%s; %s Sen's slope; filled dots = BH-FDR-significant under BOTH VC-HR98 and TFPW, open circles = one test only (α = 0.05).", series, .mcol(method)$lab)

# Figure X (required): full monthly series, primary method = MS_FULL_METHOD; the other method as a companion.
ms_fig_pdf <- file.path(fig_dir, "FigX_MS_Spatial_Drought_Trends_SenSlope_3_12_24_36.pdf")
ms_fig_png <- sub("\\.pdf$", ".png", ms_fig_pdf)
draw_sen_figure(.panels_full(c("spi", "spei", "spei_era5"), c(3, 12, 24, 36)), plot_data, MS_FULL_METHOD, c(3L, 4L),
                "Spatial drought trends", .sub_txt(MS_FULL_METHOD, "Full monthly series"), ms_fig_pdf, ms_fig_png, required = TRUE)
alt_method <- setdiff(c("vc", "tfpw"), MS_FULL_METHOD)
alt_pdf <- sub("\\.pdf$", paste0("_", toupper(alt_method), ".pdf"), ms_fig_pdf)
.step("Figure X companion", draw_sen_figure(.panels_full(c("spi", "spei", "spei_era5"), c(3, 12, 24, 36)), plot_data, alt_method,
                                            c(3L, 4L), "Spatial drought trends (method sensitivity)",
                                            .sub_txt(alt_method, "Full monthly series"), alt_pdf, sub("\\.pdf$", ".png", alt_pdf)))

# Calendar-month manuscript figure: annual-resolution series (no overlap between consecutive years for k <= 12).
cm_fig_pdf <- file.path(fig_dir, "FigX_MS_Spatial_Drought_Trends_CalendarMonth.pdf")
cm_ms_panels <- data.table::rbindlist(lapply(c("spi", "spei", "spei_era5"), function(i) data.table::data.table(
  index_type = i, timescale = c(12L, 12L, 6L, 3L), month = c(12L, 9L, 9L, 8L),
  title = paste0(index_display(i), c("-12 Dec (Jan–Dec)", "-12 Sep (Oct–Sep)", "-6 Sep (Apr–Sep)", "-3 Aug (Jun–Aug)")),
  unit = "Sen/yr")))
if (!is.null(cm_data))
  .step("calendar-month manuscript figure", draw_sen_figure(cm_ms_panels, cm_data, CALMONTH_METHOD, c(3L, 4L),
                                                            "Spatial drought trends, annual-resolution series",
                                                            .sub_txt(CALMONTH_METHOD, "Calendar-month series"), cm_fig_pdf, sub("\\.pdf$", ".png", cm_fig_pdf)))

# Snow: calendar-month series in the snow season (full-series SWEI is QA-filtered where summer gaps exist).
snow_fig_pdf <- file.path(fig_dir, "FigS_Snow_Drought_Trends_SenSlope_SWEI3_SSPI1.pdf")
snow_fig_png <- sub("\\.pdf$", ".png", snow_fig_pdf)
snow_panels <- data.table::rbindlist(list(
  data.table::data.table(index_type = "swei", timescale = SWEI_SCALE, month = 1:4,
                         title = paste0("SWEI-3 ", month.abb[1:4], c(" (Nov–Jan)", " (Dec–Feb)", " (Jan–Mar)", " (Feb–Apr)")), unit = "Sen/yr"),
  data.table::data.table(index_type = "sspi", timescale = SSPI_SCALE, month = c(12L, 1L, 2L, 3L),
                         title = paste("SSPI-1", month.abb[c(12, 1, 2, 3)]), unit = "Sen/yr")))
snow_panels <- snow_panels[index_type %in% unique(cm_data$index_type)]
if (nrow(snow_panels)) {
  .step("snow figure", draw_sen_figure(snow_panels, cm_data, CALMONTH_METHOD, c(length(unique(snow_panels$index_type)), 4L),
                                       "Snow drought trends", .sub_txt(CALMONTH_METHOD, "Calendar-month series, snow-season months"),
                                       snow_fig_pdf, snow_fig_png))
} else cat("  Snow figure skipped: no SWEI/SSPI calendar-month results.\n")

# Runoff: SRI full monthly series + SDI hydrological-year series.
runoff_fig_pdf <- file.path(fig_dir, "FigS_Runoff_Drought_Trends_SenSlope_SRI_SDI.pdf")
runoff_fig_png <- sub("\\.pdf$", ".png", runoff_fig_pdf)
runoff_panels <- rbind(.panels_full("sri", c(3, 12, 24, 36)), .panels_full("sdi", SDI_PERIODS))
runoff_panels <- runoff_panels[index_type %in% unique(plot_data$index_type)]
if (nrow(runoff_panels)) {
  .step("runoff figure", draw_sen_figure(runoff_panels, plot_data, MS_FULL_METHOD, c(length(unique(runoff_panels$index_type)), 4L),
                                         "Runoff drought trends", .sub_txt(MS_FULL_METHOD, "SRI full monthly series; SDI hydrological-year series"),
                                         runoff_fig_pdf, runoff_fig_png))
} else cat("  Runoff figure skipped: no SRI/SDI results.\n")

################################################################################
# PART 5 – BASIN AND REGIONAL TREND FIGURES (tables written by Script 11)
################################################################################
cat("\n══ PART 5: basin and regional figures ══\n")
BT <- file.path(TREND_DIR, "basin_trends"); RT <- file.path(TREND_DIR, "regional_trends")
bf <- .read(file.path(BT, "basin_trend_full_series.csv"))
bc <- .read(file.path(BT, "basin_trend_calendar_month.csv"))
bm <- .read(file.path(BT, "basin_index_monthly_all.csv"))
rf <- .read(file.path(RT, "regional_trend_full_series.csv"))
if (is.null(bf)) cat("  ⚠ basin_trends/ not found — re-run 11trend_test_ALL.R.\n")
if (!is.null(bm)) {
  bm[, date := as.Date(date)]
  bm[, `:=`(year = as.integer(format(date, "%Y")), month = as.integer(format(date, "%m")))]
}

# B1: how the autocorrelation treatment changes basin-mean significance.
fig_basin_pvalues <- function() {
  if (is.null(bf)) return(invisible(NULL))
  l <- data.table::melt(bf[, .(index_type, scale, `VC-HR98` = p_vc, TFPW = p_tfpw, `Block bootstrap` = mbb_p, `Plain MK` = p_mk_plain)],
                        id.vars = c("index_type", "scale"), variable.name = "method", value.name = "p")
  l <- l[is.finite(p)]
  l[, lab := factor(index_display(index_type), levels = index_display(INDEX_ORDER))]
  p <- ggplot(l, aes(factor(scale), -log10(pmax(p, 1e-12)), colour = method, shape = method)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    geom_point(size = 2.6, position = position_dodge(width = 0.6)) + facet_wrap(~lab, scales = "free_x", nrow = 2) +
    scale_colour_brewer(palette = "Dark2") + theme_bw(base_size = 10) +
    labs(title = "Basin-mean full-series trend: effect of the autocorrelation treatment on significance",
         subtitle = "Plain MK ignores autocorrelation; VC-HR98 corrects the variance; block bootstrap resamples blocks of max(12, 2k) months. Dashed = p 0.05.",
         x = "Accumulation scale (months; SDI: months from October)", y = expression(-log[10](p)), colour = NULL, shape = NULL)
  .save_png(file.path(FIGD[["basin"]], "B1_basin_pvalue_by_method.png"), p, 13, 7)
  data.table::fwrite(bf[, .(index = .row_lab(index_type, scale), series, n, qa, tau, sen, slope_unit, n_over_neff, rho1,
                            p_mk_plain, p_vc, p_tfpw, mbb_p, mbb_block)],
                     file.path(FIGD[["basin"]], "B1_basin_trend_table.csv"))
}

# B2: basin calendar-month tau heatmap.
fig_basin_calmonth_heat <- function() {
  if (is.null(bc)) return(invisible(NULL))
  d <- data.table::copy(bc); d[, timescale := scale]
  lv <- .row_levels(d)
  d[, `:=`(rowlab = factor(.row_lab(index_type, scale), levels = rev(lv)), mon = factor(month.abb[month], levels = month.abb),
           star = ifelse(is.finite(p_vc_bh_months) & p_vc_bh_months < 0.05, "**", ifelse(is.finite(p_vc) & p_vc < 0.05, "*", "")))]
  p <- ggplot(d, aes(mon, rowlab, fill = tau)) + geom_tile(colour = "white") + geom_text(aes(label = star), size = 3) +
    scale_fill_gradientn(colours = PAL_BRBG, limits = c(-0.4, 0.4), oob = scales::squish, na.value = "grey90", name = "τ (VC-HR98)") +
    theme_minimal(base_size = 10) +
    labs(title = "Basin-mean calendar-month trends", x = NULL, y = NULL,
         subtitle = "* p < 0.05 (VC-HR98); ** also BH-significant across the 12 months; grey = QA-filtered")
  .save_png(file.path(FIGD[["basin"]], "B2_basin_calendar_month_tau_heatmap.png"), p, 10, 0.3 * length(lv) + 2)
}

# B3: basin calendar-month series with Sen trend lines (Script 9's calendar-month plots, MK instead of OLS).
.sen_lines <- function(d, tr) {
  ln <- d[tr[, .(month, sen, tau, p_vc)], on = "month", nomatch = 0L][is.finite(sen), list(
    x0 = min(year), x1 = max(year), icpt = stats::median(value - sen * year), sen = sen[1], tau = tau[1], p = p_vc[1]), by = month]
  ln[, `:=`(y0 = icpt + sen * x0, y1 = icpt + sen * x1)]
  ln
}
fig_basin_calmonth_ts <- function(idx, sc) {
  d <- bm[index_type == idx & scale == sc & is.finite(value)]
  if (!nrow(d)) return(invisible(NULL))
  ln <- .sen_lines(d, bc[index_type == idx & scale == sc])
  ln[, lab := sprintf("τ %+.2f | p %s", tau, .fmt_p(p))]
  d[, mon := factor(month.abb[month], levels = month.abb)]; ln[, mon := factor(month.abb[month], levels = month.abb)]
  p <- ggplot(d, aes(year, value)) + geom_hline(yintercept = 0, colour = "grey70") +
    geom_hline(yintercept = DROUGHT_ONSET, colour = "grey70", linetype = "dotted") +
    geom_line(colour = "grey55", linewidth = 0.3) + geom_point(colour = .icol(idx), size = 0.9) +
    geom_segment(data = ln, aes(x = x0, xend = x1, y = y0, yend = y1), inherit.aes = FALSE, linewidth = 0.8) +
    geom_text(data = ln, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -0.05, vjust = 1.3, size = 2.5) +
    facet_wrap(~mon, ncol = 4) + theme_bw(base_size = 9) +
    labs(title = sprintf("%s basin mean — calendar-month series with Sen trend (VC-HR98)", .row_lab(idx, sc)),
         subtitle = sprintf("Dotted = drought onset (%.1f)", DROUGHT_ONSET), x = NULL, y = index_display(idx))
  .save_png(file.path(FIGD[["basin"]], sprintf("B3_basin_calendar_month_%s_%02d.png", idx, sc)), p, 12, 8)
}

# B4: seasonal index (3-month index ending Feb / May / Aug / Nov) with Sen trends.
fig_basin_seasonal <- function(idx) {
  sea <- c(`2` = "DJF (Feb)", `5` = "MAM (May)", `8` = "JJA (Aug)", `11` = "SON (Nov)")
  d <- bm[index_type == idx & scale == 3L & month %in% c(2L, 5L, 8L, 11L) & is.finite(value)]
  if (!nrow(d)) return(invisible(NULL))
  d[, season := factor(sea[as.character(month)], levels = sea)]
  ln <- .sen_lines(d, bc[index_type == idx & scale == 3L & month %in% c(2L, 5L, 8L, 11L)])
  ln[, season := factor(sea[as.character(month)], levels = sea)]
  subt <- paste(sprintf("%s: τ %+.2f, p %s", ln$season, ln$tau, .fmt_p(ln$p)), collapse = "  |  ")
  p <- ggplot(d, aes(year, value, colour = season)) + geom_hline(yintercept = 0, colour = "grey70") +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    geom_segment(data = ln, aes(x = x0, xend = x1, y = y0, yend = y1, colour = season), inherit.aes = FALSE,
                 linetype = "dashed", linewidth = 1.1) +
    scale_colour_manual(values = setNames(c("#4472C4", "#70AD47", "#C00000", "#ED7D31"), sea), name = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "top") +
    labs(title = sprintf("Seasonal %s-3 (basin mean) with Sen trends", index_display(idx)), subtitle = subt,
         x = NULL, y = sprintf("%s-3", index_display(idx)))
  .save_png(file.path(FIGD[["basin"]], sprintf("B4_basin_seasonal_%s_03.png", idx)), p, 13, 6)
}

# B5: regional full-series tau heatmap.
fig_regional_heat <- function() {
  if (is.null(rf)) return(invisible(NULL))
  d <- data.table::copy(rf); d[, timescale := scale]
  lv <- .row_levels(d)
  reg_lv <- rev(unique(d$region_name[order(d$region_id != BASIN_RID, d$region_name)]))
  d[, `:=`(collab = factor(.row_lab(index_type, scale), levels = lv), reg = factor(region_name, levels = reg_lv),
           star = ifelse(is.finite(p_vc) & p_vc < 0.05, "*", ""))]
  p <- ggplot(d, aes(collab, reg, fill = tau)) + geom_tile(colour = "white") + geom_text(aes(label = star), size = 3.2) +
    scale_fill_gradientn(colours = PAL_BRBG, limits = c(-0.3, 0.3), oob = scales::squish, na.value = "grey90", name = "τ (VC-HR98)") +
    theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Regional full-series trends (monthly indices; SDI annual)",
         subtitle = "* p < 0.05 (VC-HR98). Sub-basins are nested in the basin, so rows are not independent.", x = NULL, y = NULL)
  .save_png(file.path(FIGD[["basin"]], "B5_regional_tau_heatmap.png"), p, 14, 5)
}

# B6: regional index time series (Script 9 regional climate plots, for the drought indices).
fig_regional_timeseries <- function() {
  fs <- list.files(file.path(TREND_DIR, "regional_timeseries"), pattern = "_regional_averaged_by_month\\.csv$", full.names = TRUE)
  keep <- function(b) grepl(sprintf("^(spi|spei|spei_era5|sri)_%02d_", POINT_TS_SCALE), b) ||
    grepl(sprintf("^swei_%02d_", SWEI_SCALE), b) || grepl(sprintf("^sspi_%02d_", SSPI_SCALE), b)
  fs <- fs[vapply(basename(fs), keep, logical(1))]
  if (!length(fs)) return(invisible(NULL))
  d <- data.table::rbindlist(lapply(fs, data.table::fread), fill = TRUE)
  d[, `:=`(date = as.Date(date), lab = .row_lab(index_type, scale))]
  d[, lab := factor(lab, levels = unique(lab[order(match(index_type, INDEX_ORDER))]))]
  for (rid in unique(d$region_id)) {
    s <- d[region_id == rid & is.finite(value)]
    if (!nrow(s)) next
    p <- ggplot(s, aes(date, value)) + geom_hline(yintercept = 0, colour = "grey70") +
      geom_hline(yintercept = DROUGHT_ONSET, colour = "grey70", linetype = "dotted") +
      geom_line(linewidth = 0.3, colour = "grey25") + facet_wrap(~lab, ncol = 1) + theme_bw(base_size = 9) +
      labs(title = unique(s$region_name)[1], x = NULL, y = NULL,
           subtitle = "Regional area-weighted index series (native ERA5-Land grid; months with < 20% valid area are blank)")
    .save_png(file.path(FIGD[["basin"]], sprintf("B6_regional_timeseries_%s.png", rid)), p, 12, 1.8 * nlevels(droplevels(s$lab)) + 1.2)
  }
}

.step("B1", fig_basin_pvalues())
.step("B2", fig_basin_calmonth_heat())
if (!is.null(bm) && !is.null(bc)) {
  b3 <- rbind(data.table::CJ(index_type = c("spi", "spei", "spei_era5", "sri"), scale = CALMONTH_TS_SCALES, sorted = FALSE),
              data.table::data.table(index_type = c("swei", "sspi"), scale = c(SWEI_SCALE, SSPI_SCALE)))
  for (i in seq_len(nrow(b3))) .step(sprintf("B3 %s-%d", b3$index_type[i], b3$scale[i]), fig_basin_calmonth_ts(b3$index_type[i], b3$scale[i]))
  for (idx in c("spi", "spei", "spei_era5", "sri")) .step(sprintf("B4 %s", idx), fig_basin_seasonal(idx))
}
.step("B5", fig_regional_heat())
.step("B6", fig_regional_timeseries())
################################################################################
# PART 6 – HOMOGENEITY FIGURES (tables written by Script 11, annual-resolution series)
################################################################################
cat("\n══ PART 6: homogeneity figures ══\n")
HOMOG_CLASS <- c(sprintf("Step at %d supported (AR(1) step + trend)", STEP_YEAR),
                 "Pettitt break only (trend-like)", "No break")
HOMOG_COL <- setNames(c("#b2182b", "#ef8a62", "grey60"), HOMOG_CLASS)
.homog_class <- function(dAIC, step_p, pett_p) data.table::fifelse(
  is.finite(dAIC) & dAIC < -2 & is.finite(step_p) & step_p < 0.05, HOMOG_CLASS[1],
  data.table::fifelse(is.finite(pett_p) & pett_p < 0.05, HOMOG_CLASS[2], HOMOG_CLASS[3]))

hb   <- .read(file.path(HOMOG_DIR, "homogeneity_break_summary_basin.csv"))
hr   <- .read(file.path(HOMOG_DIR, "step_tests_by_region_all.csv"))
hp   <- .read(file.path(HOMOG_DIR, "pixel_step_summary.csv"))
hs   <- .read(file.path(HOMOG_DIR, "homogeneity_annual_series_all.csv"))
hsub <- .read(file.path(HOMOG_DIR, "subperiod_trends_all.csv"))

.series_lv <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (!length(parts)) return(character(0))
  u <- unique(data.table::rbindlist(lapply(parts, function(d) unique(d[, .(series, index_type, scale)])), fill = TRUE))
  u[, o := match(index_type, INDEX_ORDER)]
  data.table::setorder(u, o, -scale, series)
  unique(u$series)
}
SERIES_LV <- .series_lv(hb, hr, hp)

# H1: basin break year, Pettitt strength and step-vs-trend class for every series.
fig_homog_basin <- function() {
  if (is.null(hb) || !nrow(hb)) {
    cat("  ⚠ homogeneity outputs not found — re-run 11trend_test_ALL.R with RUN_HOMOGENEITY = TRUE.\n")
    return(invisible(NULL))
  }
  b <- data.table::copy(hb)
  if (!"pettitt_p_bh_across_series" %in% names(b)) b[, pettitt_p_bh_across_series := NA_real_]
  if (!"family" %in% names(b)) b[, family := "index"]
  b[, class := factor(.homog_class(delta_AIC_step_vs_nostep, step_p, pettitt_p), levels = HOMOG_CLASS)]
  b[, series := factor(series, levels = rev(SERIES_LV))]
  b[, lab := sprintf("Pettitt p %s (BH %s) | dAIC %s", .fmt_p(pettitt_p), .fmt_p(pettitt_p_bh_across_series),
                     ifelse(is.finite(delta_AIC_step_vs_nostep), sprintf("%+.1f", delta_AIC_step_vs_nostep), "NA"))]
  data.table::fwrite(b, file.path(HOMOG_DIR, "homogeneity_break_classification_basin.csv"))
  yr <- range(c(b$first_year, b$last_year), na.rm = TRUE); mid <- mean(yr)
  p <- ggplot(b, aes(pettitt_first_year_after, series)) +
    annotate("rect", xmin = STEP_YEAR - 3, xmax = STEP_YEAR + 3, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
    geom_vline(xintercept = STEP_YEAR, linetype = "dashed", colour = "grey35") +
    geom_vline(xintercept = mid, linetype = "dotted", colour = "grey35") +
    geom_point(aes(colour = class, size = -log10(pmax(pettitt_p, 1e-10))), na.rm = TRUE) +
    geom_text(aes(x = yr[2] + 2, label = lab), hjust = 0, size = 2.3, colour = "grey30") +
    facet_grid(family ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = HOMOG_COL, drop = FALSE, name = NULL) +
    scale_size_continuous(range = c(1.2, 4.5), name = "-log10 Pettitt p") +
    scale_x_continuous(limits = c(yr[1], yr[2] + 40), breaks = seq(1950, 2020, 10)) +
    theme_bw(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0), legend.position = "bottom", legend.box = "vertical") +
    labs(title = "Break year and step-vs-trend classification — basin-mean annual-resolution index series",
         subtitle = sprintf(paste0("Point = Pettitt first year of the later regime. Dashed = a-priori step year %d (band ±3 yr); ",
                                   "dotted = record midpoint, where a smooth trend places a Pettitt break.\n",
                                   "Red = AR(1)-GLS supports a step at %d in addition to a trend (dAIC < -2, step p < 0.05), plotted at its Pettitt year; ",
                                   "orange = Pettitt break that the AR(1) model treats as a trend."), STEP_YEAR, STEP_YEAR),
         x = "First year of the later regime (Pettitt)", y = NULL)
  .save_png(file.path(FIGD[["homog"]], "H1_basin_break_classification.png"), p, 12.5, 0.32 * nrow(b) + 3)
  tab <- b[, .N, by = class]
  cat("  Basin series:", paste(sprintf("%s = %d", tab$class, tab$N), collapse = " | "), "\n")
}

# H2: classification by region.
fig_homog_regions <- function() {
  if (is.null(hr) || !nrow(hr)) return(invisible(NULL))
  r <- data.table::copy(hr)
  r[, class := factor(.homog_class(delta_AIC_step_vs_nostep, step_p, pettitt_p), levels = HOMOG_CLASS)]
  reg_lv <- c(intersect(BASIN_RID, unique(r$region_id)), setdiff(sort(unique(r$region_id)), BASIN_RID))
  r[, `:=`(region = factor(region_id, levels = reg_lv), series = factor(series, levels = rev(SERIES_LV)),
           lab = ifelse(class == HOMOG_CLASS[3] | !is.finite(pettitt_first_year_after), "", as.character(pettitt_first_year_after)))]
  p <- ggplot(r, aes(region, series, fill = class)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 2.2) +
    scale_fill_manual(values = HOMOG_COL, drop = FALSE, name = NULL) + theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom") +
    labs(title = "Step-vs-trend classification by region",
         subtitle = "Label = Pettitt first year of the later regime. Sub-basins are nested in the basin, so columns are not independent.",
         x = NULL, y = NULL)
  .save_png(file.path(FIGD[["homog"]], "H2_classification_by_region.png"), p, 11, 0.26 * length(SERIES_LV) + 2.8)
}

# H3: share of pixels with a Pettitt break; H4: break-year and step-size maps.
fig_homog_pixels <- function() {
  if (is.null(hp) || !nrow(hp)) return(invisible(NULL))
  p <- data.table::copy(hp)[order(pct_pettitt_fdr05)]
  p[, series := factor(series, levels = unique(series))]
  p[, lab := ifelse(is.finite(median_break_year_fdr_sig),
                    sprintf("median break %.0f | Welch-FDR %.0f%%", median_break_year_fdr_sig, pct_welch_fdr05),
                    sprintf("Welch-FDR %.0f%%", pct_welch_fdr05))]
  g <- ggplot(p, aes(pct_pettitt_fdr05, series)) + geom_col(fill = "#4575b4", width = 0.7, na.rm = TRUE) +
    geom_text(aes(label = lab), hjust = -0.05, size = 2.5) +
    scale_x_continuous(limits = c(0, 150), breaks = seq(0, 100, 25)) + theme_bw(base_size = 9) +
    labs(title = "Share of basin pixels with a Pettitt break (BH-FDR < 0.05)",
         subtitle = sprintf("Welch-FDR = share of pixels whose mean changes at %d (BH-FDR < 0.05).", STEP_YEAR),
         x = "% of tested pixels", y = NULL)
  .save_png(file.path(FIGD[["homog"]], "H3_pixel_break_share.png"), g, 10, 0.3 * nrow(p) + 2)
  
  pm <- data.table::rbindlist(lapply(as.character(p$series), function(s) {
    f <- file.path(HOMOG_DIR, "pixel", sprintf("pixel_step_tests_%s.csv", s))
    if (!file.exists(f)) return(NULL)
    x <- data.table::fread(f); x[, series := s]; x
  }), fill = TRUE)
  if (!nrow(pm)) return(invisible(NULL))
  rs <- terra::res(analysis_template)
  lv <- intersect(SERIES_LV, unique(pm$series)); nc <- min(5L, length(lv))
  pm[, series := factor(series, levels = lv)]
  pm[, brk := ifelse(is.finite(pettitt_p_fdr) & pettitt_p_fdr < 0.05, pettitt_first_year_after, NA_real_)]
  map_base <- function(d) ggplot(d, aes(lon, lat)) + facet_wrap(~series, ncol = nc) + coord_sf(expand = FALSE) +
    theme_void(base_size = 9) + theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 7.5))
  g4 <- map_base(pm) + geom_tile(aes(fill = brk), width = rs[1], height = rs[2]) +
    geom_sf(data = basin_sf, inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.35) +
    scale_fill_viridis_c(option = "plasma", na.value = "grey90", name = "Break year (Pettitt FDR < 0.05)") +
    labs(title = "Pixel break years (grey = no significant break)")
  .save_png(file.path(FIGD[["homog"]], "H4_pixel_break_year_maps.png"), g4, 3.1 * nc, 2.9 * ceiling(length(lv) / nc) + 1.2)
  z <- max(abs(pm$diff_at_step), na.rm = TRUE); if (!is.finite(z) || z <= 0) z <- 1
  sig <- pm[is.finite(welch_p_fdr) & welch_p_fdr < 0.05]
  g4b <- map_base(pm) + geom_tile(aes(fill = diff_at_step), width = rs[1], height = rs[2]) +
    geom_point(data = sig, size = 0.25, colour = "black") +
    geom_sf(data = basin_sf, inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.35) +
    scale_fill_gradientn(colours = PAL_BRBG, limits = c(-z, z), na.value = "grey90",
                         name = sprintf("Mean from %d minus mean before (index units)", STEP_YEAR)) +
    labs(title = sprintf("Change in the mean at %d (dots = Welch BH-FDR < 0.05)", STEP_YEAR))
  .save_png(file.path(FIGD[["homog"]], "H4b_pixel_step_change_maps.png"), g4b, 3.1 * nc, 2.9 * ceiling(length(lv) / nc) + 1.2)
}

# H5: basin annual series with period means, a-priori step year and Pettitt break.
fig_homog_series <- function() {
  if (is.null(hs) || is.null(hr)) return(invisible(NULL))
  d <- hs[region_id == BASIN_RID & is.finite(value)]
  m <- hr[region_id == BASIN_RID]
  if (!nrow(d) || !nrow(m)) return(invisible(NULL))
  seg <- rbind(m[is.finite(mean_before), .(series, x0 = first_year, x1 = STEP_YEAR - 0.5, y = mean_before)],
               m[is.finite(mean_after),  .(series, x0 = STEP_YEAR - 0.5, x1 = last_year, y = mean_after)])
  lab <- m[, .(series, txt = sprintf("Pettitt %s (break %s) | Welch@%d %s | dAIC %s", .fmt_p(pettitt_p),
                                     ifelse(is.finite(pettitt_first_year_after), as.character(pettitt_first_year_after), "–"),
                                     STEP_YEAR, .fmt_p(welch_p),
                                     ifelse(is.finite(delta_AIC_step_vs_nostep), sprintf("%+.1f", delta_AIC_step_vs_nostep), "NA")))]
  brk <- m[is.finite(pettitt_first_year_after), .(series, xb = pettitt_first_year_after - 0.5)]
  lv <- intersect(SERIES_LV, unique(d$series))
  for (x in list(d, seg, lab, brk)) x[, series := factor(series, levels = lv)]
  p <- ggplot(d, aes(year, value)) + geom_hline(yintercept = 0, colour = "grey75") +
    geom_line(colour = "grey55", linewidth = 0.3) + geom_point(size = 0.8, colour = "grey15") +
    geom_segment(data = seg, aes(x = x0, xend = x1, y = y, yend = y), colour = "firebrick", linewidth = 0.9, inherit.aes = FALSE) +
    geom_vline(xintercept = STEP_YEAR - 0.5, linetype = "dashed", colour = "grey40") +
    geom_vline(data = brk, aes(xintercept = xb), colour = "darkorange", linetype = "dotted", linewidth = 0.8) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt), hjust = -0.02, vjust = 1.3, size = 2.0, inherit.aes = FALSE) +
    facet_wrap(~series, ncol = 4, scales = "free_y") + theme_bw(base_size = 8.5) +
    labs(title = "Basin-mean annual-resolution index series",
         subtitle = sprintf("Red = period means before/from %d; dashed = a-priori step year; dotted orange = Pettitt break.", STEP_YEAR),
         x = NULL, y = "Index value")
  .save_png(file.path(FIGD[["homog"]], "H5_basin_annual_series_step_means.png"), p, 14, 2.3 * ceiling(length(lv) / 4) + 1.2)
}

# H6: sub-period trends (full record / before / from step year / WMO 1991-2020).
fig_homog_subperiod <- function() {
  if (is.null(hsub)) return(invisible(NULL))
  s <- hsub[region_id == BASIN_RID & is.finite(sen)]
  if (!nrow(s)) return(invisible(NULL))
  s[, `:=`(series = factor(series, levels = rev(intersect(SERIES_LV, unique(series)))),
           window = factor(window, levels = c("full", "before_step", "from_step", "wmo_1991_2020"),
                           labels = c("Full record", sprintf("Before %d", STEP_YEAR), sprintf("%d onward", STEP_YEAR), "1991–2020")),
           sig = ifelse(is.finite(p_vc) & p_vc < 0.05, "p < 0.05", "n.s."))]
  p <- ggplot(s, aes(sen, series, colour = window, shape = sig)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(size = 2.2, position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c(`p < 0.05` = 16, `n.s.` = 1), name = NULL) +
    scale_colour_brewer(palette = "Dark2", name = NULL) + theme_bw(base_size = 9) + theme(legend.position = "bottom") +
    labs(title = "Sub-period trends of the basin-mean annual-resolution series",
         subtitle = "Sen's slope per year; filled = VC-HR98 p < 0.05. Short sub-periods have low power.",
         x = "Sen's slope (index units per year)", y = NULL)
  .save_png(file.path(FIGD[["homog"]], "H6_subperiod_trends.png"), p, 10, 0.3 * nlevels(s$series) + 2.5)
}

.step("H1", fig_homog_basin())
.step("H2", fig_homog_regions())
.step("H3/H4", fig_homog_pixels())
.step("H5", fig_homog_series())
.step("H6", fig_homog_subperiod())

################################################################################
# PART 7 – SPECIFIC POINTS (Script 9 points; written by Script 11)
################################################################################
cat("\n══ PART 7: specific points ══\n")
PT <- file.path(TREND_DIR, "points")
pts_stats <- .read(file.path(PT, "point_trend_stats_drought.csv"))
pts_ts    <- .read(file.path(PT, "point_index_timeseries.csv"))
.pt_select <- function(d, scales_main) d[(index_type %in% c("spi", "spei", "spei_era5", "sri") & scale %in% scales_main) |
                                           (index_type == "swei" & scale == SWEI_SCALE) |
                                           (index_type == "sspi" & scale == SSPI_SCALE)]
.pt_title <- function(pid) {
  loc <- if (!is.null(pts_stats)) pts_stats[point_id == pid][1] else NULL
  if (!is.null(loc) && nrow(loc) && is.finite(loc$lat_wgs84))
    sprintf("Point #%d (%.4f°N, %.4f°W)", pid, loc$lat_wgs84, abs(loc$lon_wgs84)) else sprintf("Point #%d", pid)
}

fig_point_ts <- function(pid) {
  d <- .pt_select(pts_ts[point_id == pid & is.finite(value)], POINT_TS_SCALE)
  if (!nrow(d)) return(invisible(NULL))
  d[, `:=`(date = as.Date(date), lab = .row_lab(index_type, scale))]
  lv <- unique(d[order(match(index_type, INDEX_ORDER)), lab]); d[, lab := factor(lab, levels = lv)]
  tr <- if (!is.null(pts_stats)) pts_stats[point_id == pid & period == "full_monthly",
                                           .(lab = .row_lab(index_type, scale),
                                             txt = sprintf("full series: τ %+.2f | p(VC) %s", tau_vc, .fmt_p(p_value_vc)))] else NULL
  p <- ggplot(d, aes(date, value)) + geom_hline(yintercept = 0, colour = "grey70") +
    geom_hline(yintercept = DROUGHT_ONSET, linetype = "dotted", colour = "grey60") +
    geom_line(linewidth = 0.3, colour = "grey25") + facet_wrap(~lab, ncol = 1) + theme_bw(base_size = 9) +
    labs(title = .pt_title(pid), subtitle = sprintf("Monthly index series at the nearest native cell; dotted = %.1f", DROUGHT_ONSET),
         x = NULL, y = NULL)
  if (!is.null(tr) && nrow(tr)) {
    tr <- tr[lab %in% lv]; tr[, `:=`(lab = factor(lab, levels = lv), x = min(d$date))]
    p <- p + geom_text(data = tr, aes(x = x, y = Inf, label = txt), hjust = 0, vjust = 1.3, size = 2.6, inherit.aes = FALSE)
  }
  .save_png(file.path(FIGD[["points"]], sprintf("P1_point_%d_index_timeseries.png", pid)), p, 12, 1.9 * length(lv) + 1.2)
}

fig_point_calmonth <- function(pid) {
  if (is.null(pts_stats)) return(invisible(NULL))
  d <- .pt_select(pts_stats[point_id == pid & period == "calendar_month" & is.finite(tau_vc)], c(3L, 12L))
  if (!nrow(d)) return(invisible(NULL))
  d[, `:=`(lab = .row_lab(index_type, scale), mon = factor(month.abb[month], levels = month.abb),
           star = ifelse(is.finite(p_value_vc) & p_value_vc < 0.05, "*", ""))]
  lv <- unique(d[order(match(index_type, INDEX_ORDER), scale), lab]); d[, lab := factor(lab, levels = lv)]
  p <- ggplot(d, aes(mon, tau_vc, fill = tau_vc > 0)) + geom_col(width = 0.7) + geom_hline(yintercept = 0) +
    geom_text(aes(y = tau_vc + sign(tau_vc) * 0.03, label = star), size = 4) +
    scale_fill_manual(values = c(`TRUE` = "#35978F", `FALSE` = "#BF812D"), guide = "none") +
    scale_x_discrete(labels = substr(month.abb, 1, 1)) + facet_wrap(~lab, ncol = 4) + theme_bw(base_size = 9) +
    labs(title = paste(.pt_title(pid), "— calendar-month trends"),
         subtitle = "Kendall's τ (VC-HR98) of each calendar month across years; * = p < 0.05 (single-site, nominal)",
         x = NULL, y = "τ")
  .save_png(file.path(FIGD[["points"]], sprintf("P2_point_%d_calendar_month_tau.png", pid)), p, 12, 2.4 * ceiling(length(lv) / 4) + 1.4)
}

if (is.null(pts_ts) || is.null(pts_stats)) {
  cat("  ⚠ point outputs not found (no specific points in Script 9 metadata, or Script 11 not run).\n")
} else {
  for (pid in unique(pts_ts$point_id)) {
    .step(sprintf("P1 point %d", pid), fig_point_ts(pid))
    .step(sprintf("P2 point %d", pid), fig_point_calmonth(pid))
  }
  pt_tab <- pts_stats[period %in% c("full_monthly", "annual_hydro_year"),
                      intersect(c("point_id", "lat_wgs84", "lon_wgs84", "index_type", "scale", "period", "qa_full",
                                  "rho1_detrended", "tau_vc", "p_value_vc", "sl_vc", "tau_tfpw", "p_value_tfpw", "sl_tfpw",
                                  "dfreq_tau_vc", "dfreq_p_vc"), names(pts_stats)), with = FALSE]
  data.table::fwrite(pt_tab, file.path(FIGD[["points"]], "P3_point_full_series_trends.csv"))
  cat(sprintf("  ✓ P3_point_full_series_trends.csv (%d rows)\n", nrow(pt_tab)))
}

################################################################################
# PART 8 – PIXEL-TO-BASIN SPEARMAN CORRELATION
# Native ERA5-Land grid only: every seasonal NetCDF must match it exactly (no resampling).
# rho_i = Spearman(X_i(t), X_basin(t)) with the basin series written by Script 11.
################################################################################
cat("\n══ PART 8: pixel-to-basin Spearman correlation ══\n")
CORR_OUT_DIR <- file.path(fig_dir, "pixel_basin_correlation")
dir.create(CORR_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
CORR_SHORT_SCALES <- c(1L, 3L, 6L); CORR_LONG_SCALES <- c(12L, 24L, 36L)
CORR_ALL_SCALES <- c(CORR_SHORT_SCALES, CORR_LONG_SCALES)
CORR_MIN_N <- 24L
CORR_REQUIRED <- c("spi", "spei")       # failures for these stop the script; others are reported and skipped
.has_basin_csv <- function(dir, idx, sc) file.exists(file.path(dir, sprintf("%s_%02d_basin_averaged_by_month.csv", idx, sc)))
.avail <- function(dir, idx, scs) as.integer(scs[vapply(scs, function(s) .has_basin_csv(dir, idx, s), logical(1))])
CORR_SCALES_BY_INDEX <- list(
  spi = CORR_ALL_SCALES, spei = CORR_ALL_SCALES,
  spei_era5 = .avail(SPEI_ERA5_SEAS_DIR, "spei_era5", CORR_ALL_SCALES),
  swei = .avail(SWEI_SEAS_DIR, "swei", SWEI_SCALE),
  sspi = .avail(SSPI_SEAS_DIR, "sspi", SSPI_SCALE),
  sri = .avail(SRI_SEAS_DIR, "sri", CORR_ALL_SCALES))
CORR_TEMPLATE <- analysis_template
CORR_BASIN_CELLS <- basin_mask_cells
if (!isTRUE(terra::compareGeom(CORR_TEMPLATE, basin_mask_raster, stopOnError = FALSE)) ||
    length(CORR_BASIN_CELLS) != terra::ncell(CORR_TEMPLATE) || sum(CORR_BASIN_CELLS) != basin_pixels)
  stop("Native template and basin mask are inconsistent.")

validate_corr_geometry <- function(r, label) {
  ok <- isTRUE(terra::compareGeom(r, CORR_TEMPLATE, lyrs = FALSE, crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE,
                                  stopOnError = FALSE)) &&
    isTRUE(all.equal(as.numeric(terra::origin(r)), as.numeric(terra::origin(CORR_TEMPLATE)), tolerance = 1e-8))
  if (!ok) stop("GRID ERROR [", label, "]: geometry differs from the native ERA5-Land grid (no resampling allowed).")
  invisible(TRUE)
}

find_corr_nc_files <- function(index_type, scale) {
  idx <- tolower(index_type)
  files <- switch(idx,
                  spi       = find_seasonal_nc_files(SPI_SEAS_DIR, "spi", scale),
                  spei      = find_seasonal_nc_files(SPEI_SEAS_DIR, "spei", scale),
                  spei_era5 = find_seasonal_nc_files(SPEI_ERA5_SEAS_DIR, "spei", scale),
                  swei      = find_swei_seasonal_files(SWEI_SEAS_DIR, "swei", scale),
                  sspi      = find_seasonal_nc_files(SSPI_SEAS_DIR, "sspi", scale),
                  sri       = find_seasonal_nc_files(SRI_SEAS_DIR, "sri", scale),
                  stop("Unsupported correlation index: ", index_type))
  if (!length(files)) stop(sprintf("No seasonal NetCDF files found for %s-%d", toupper(idx), scale))
  mo <- suppressWarnings(as.integer(sub(".*month([0-9]{2}).*", "\\1", basename(files))))
  if (!anyNA(mo)) files <- files[order(mo)]          # files[m] must be calendar month m
  files
}

read_corr_index_matrix <- function(idx, sc, target_dates) {
  files <- find_corr_nc_files(idx, sc)
  out <- matrix(NA_real_, terra::ncell(CORR_TEMPLATE), length(target_dates))
  if (length(files) == 1L) {
    r <- terra::rast(files); validate_corr_geometry(r, basename(files))
    if (terra::nlyr(r) != length(target_dates)) stop(sprintf("%s: %d layers vs %d basin months", basename(files), terra::nlyr(r), length(target_dates)))
    out[] <- as.matrix(terra::values(r, mat = TRUE))
  } else {
    if (length(files) != 12L) stop(sprintf("Expected 12 seasonal files for %s-%d; found %d", toupper(idx), sc, length(files)))
    tm <- as.integer(format(target_dates, "%m"))
    for (m in 1:12) {
      r <- terra::rast(files[m]); validate_corr_geometry(r, basename(files[m]))
      v <- as.matrix(terra::values(r, mat = TRUE)); cols <- which(tm == m)
      if (ncol(v) != length(cols)) stop(sprintf("%s: %d layers vs %d target dates", basename(files[m]), ncol(v), length(cols)))
      out[, cols] <- v
    }
  }
  out[!is.finite(out)] <- NA_real_
  out
}

compute_corr_raster <- function(ts_matrix, basin_series) {
  if (nrow(ts_matrix) != terra::ncell(CORR_TEMPLATE) || ncol(ts_matrix) != length(basin_series)) stop("Correlation matrix dimension mismatch.")
  n_pair <- integer(nrow(ts_matrix)); rr <- rep(NA_real_, nrow(ts_matrix))
  for (i in which(CORR_BASIN_CELLS)) {
    ok <- is.finite(ts_matrix[i, ]) & is.finite(basin_series); n_pair[i] <- sum(ok)
    if (n_pair[i] < CORR_MIN_N) next
    x <- ts_matrix[i, ok]; y <- basin_series[ok]
    if (length(unique(x)) < 2L || length(unique(y)) < 2L) next
    rr[i] <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  }
  list(raster = terra::setValues(terra::rast(CORR_TEMPLATE), rr), n_pair = n_pair)
}

plot_corr_map <- function(r, label) {
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  if (!nrow(d)) stop("No finite correlation values for ", label)
  names(d)[3] <- "rho"
  ggplot(d, aes(x, y, fill = rho)) + geom_raster() +
    geom_sf(data = basin_sf, inherit.aes = FALSE, fill = NA, linewidth = 0.35) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1),
                         oob = scales::squish, name = expression(rho[S])) +
    coord_sf(expand = FALSE) + labs(title = label) + theme_bw(base_size = 10) +
    theme(axis.title = element_blank(), axis.text = element_text(size = 7), panel.grid = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5))
}

plot_corr_density <- function(corr_values, index_type, scales) {
  d <- data.table::rbindlist(lapply(scales, function(sc) {
    x <- corr_values[[sprintf("%s_%02d", index_type, sc)]]
    if (is.null(x)) NULL else data.table::data.table(rho = x, scale = factor(sc, levels = scales))
  }))
  d <- d[is.finite(rho)]
  if (!nrow(d)) stop("No finite correlation values for density plot: ", index_type)
  meds <- d[, .(rho = stats::median(rho)), by = scale]
  ggplot(d, aes(rho, colour = scale)) + geom_density(linewidth = 1.1, adjust = 1.05) +
    geom_vline(data = meds, aes(xintercept = rho, colour = scale), linetype = "dotted", linewidth = 0.85, show.legend = FALSE) +
    scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73")[seq_along(scales)], name = "Timescale") +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
    labs(title = paste0(index_display(index_type), ": pixel–basin Spearman correlation distribution"),
         x = expression(Spearman ~ rho[S]), y = "Density") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                     plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
}

corr_rasters <- list(); corr_values <- list(); corr_meta <- list(); corr_npair <- list()
for (idx in names(CORR_SCALES_BY_INDEX)) for (sc in CORR_SCALES_BY_INDEX[[idx]]) {
  key <- sprintf("%s_%02d", idx, sc)
  res <- tryCatch({
    bdir <- switch(idx, sspi = SSPI_SEAS_DIR, sri = SRI_SEAS_DIR, NULL)
    bdf <- if (is.null(bdir)) load_basin_avg_csv(idx, sc) else load_basin_avg_csv(idx, sc, data_dir = bdir)
    if (is.null(bdf)) stop("basin-average CSV missing (run 11trend_test_ALL.R)")
    bd <- as.Date(bdf$date)
    tsm <- read_corr_index_matrix(idx, sc, bd)
    tsm[!CORR_BASIN_CELLS, ] <- NA_real_
    compute_corr_raster(tsm, as.numeric(bdf$value))
  }, error = function(e) {
    if (idx %in% CORR_REQUIRED) stop(sprintf("Pixel-basin correlation failed for %s: %s", .row_lab(idx, sc), conditionMessage(e)))
    cat(sprintf("  ⚠ %s skipped: %s\n", .row_lab(idx, sc), conditionMessage(e))); NULL
  })
  if (is.null(res)) next
  vals <- terra::values(res$raster, mat = FALSE); vals <- vals[is.finite(vals)]
  if (!length(vals)) { cat(sprintf("  ⚠ %s: no finite correlations\n", .row_lab(idx, sc))); next }
  corr_rasters[[key]] <- res$raster; corr_values[[key]] <- vals
  corr_meta[[key]] <- data.table::data.table(
    index_type = idx, index_label = index_display(idx), scale = sc, n_pixels = length(vals),
    median_rho = stats::median(vals), mean_rho = mean(vals), min_rho = min(vals), max_rho = max(vals),
    n_paired_months_min = min(res$n_pair[CORR_BASIN_CELLS]), n_paired_months_max = max(res$n_pair[CORR_BASIN_CELLS]))
  cat(sprintf("  %-12s %4d pixels | median ρ %.3f | range %.3f to %.3f\n", .row_lab(idx, sc), length(vals),
              stats::median(vals), min(vals), max(vals)))
}
if (!length(corr_meta)) stop("No pixel-basin correlations were computed.")
data.table::fwrite(data.table::rbindlist(corr_meta), file.path(CORR_OUT_DIR, "PixelBasin_Spearman_Summary.csv"))

cell_id <- seq_len(terra::ncell(CORR_TEMPLATE)); xy <- terra::xyFromCell(CORR_TEMPLATE, cell_id)
correlation_audit <- data.table::rbindlist(lapply(names(corr_rasters), function(key) {
  vals <- terra::values(corr_rasters[[key]], mat = FALSE)
  if (length(vals) != length(cell_id)) stop("Correlation raster length mismatch for ", key)
  idx <- sub("_[0-9]+$", "", key); sc <- as.integer(sub("^.*_", "", key))
  data.table::data.table(cell = cell_id, x = xy[, 1], y = xy[, 2], index_type = idx, index_label = index_display(idx),
                         scale = sc, spearman_rho = vals)[CORR_BASIN_CELLS & is.finite(spearman_rho)]
}))
if (!nrow(correlation_audit) || !all(correlation_audit$cell %in% cell_id[CORR_BASIN_CELLS]))
  stop("Correlation audit failed: empty table or cells outside the native basin mask.")
data.table::fwrite(correlation_audit, file.path(CORR_OUT_DIR, "PixelBasin_Spearman_Correlations.csv"))

make_corr_figure <- function(indices, scales, outfile_stub) {
  indices <- indices[vapply(indices, function(i) all(sprintf("%s_%02d", i, scales) %in% names(corr_rasters)), logical(1))]
  if (!length(indices)) { cat("  ", outfile_stub, " skipped: no correlations\n", sep = ""); return(invisible(NULL)) }
  keys <- as.vector(t(outer(indices, scales, function(i, s) sprintf("%s_%02d", i, s))))
  labs <- as.vector(t(outer(indices, scales, function(i, s) sprintf("%s-%d", index_display(i), s))))
  map_grid <- patchwork::wrap_plots(Map(function(k, l) plot_corr_map(corr_rasters[[k]], l), keys, labs),
                                    nrow = length(indices), ncol = length(scales), byrow = TRUE)
  dens <- lapply(indices, function(i) plot_corr_density(corr_values, i, scales))
  n <- length(indices)
  final_plot <- patchwork::wrap_plots(c(list(map_grid), dens), ncol = 1, heights = c(n, rep(0.72, n)))
  ggplot2::ggsave(file.path(CORR_OUT_DIR, paste0(outfile_stub, ".png")), final_plot, width = 15, height = 7.5 * n,
                  units = "in", dpi = DPI, bg = "white")
  ggplot2::ggsave(file.path(CORR_OUT_DIR, paste0(outfile_stub, ".pdf")), final_plot, width = 15, height = 7.5 * n,
                  units = "in", device = grDevices::cairo_pdf, bg = "white")
  cat("  ✓ ", outfile_stub, ".png / .pdf\n", sep = "")
}
make_snow_corr_figure <- function() {
  keys <- intersect(c(sprintf("swei_%02d", SWEI_SCALE), sprintf("sspi_%02d", SSPI_SCALE)), names(corr_rasters))
  if (!length(keys)) { cat("  Snow correlation figure skipped: no snow correlations.\n"); return(invisible(NULL)) }
  labs <- vapply(keys, function(k) sprintf("%s-%d", index_display(sub("_[0-9]+$", "", k)), as.integer(sub("^.*_", "", k))), "")
  maps <- Map(function(k, l) plot_corr_map(corr_rasters[[k]], l), keys, labs)
  dens <- lapply(keys, function(k) plot_corr_density(corr_values, sub("_[0-9]+$", "", k), as.integer(sub("^.*_", "", k))))
  fp <- patchwork::wrap_plots(c(maps, dens), ncol = length(keys), byrow = TRUE) +
    patchwork::plot_annotation(title = "Snow drought pixel–basin Spearman correlations")
  ggplot2::ggsave(file.path(CORR_OUT_DIR, "FigS_PixelBasin_Spearman_SWEI3_SSPI1.png"), fp, width = 12, height = 10, dpi = DPI, bg = "white")
  ggplot2::ggsave(file.path(CORR_OUT_DIR, "FigS_PixelBasin_Spearman_SWEI3_SSPI1.pdf"), fp, width = 12, height = 10,
                  device = grDevices::cairo_pdf, bg = "white")
  cat("  ✓ FigS_PixelBasin_Spearman_SWEI3_SSPI1.png / .pdf\n")
}
.step("Spearman 1/3/6",   make_corr_figure(c("spi", "spei", "spei_era5"), CORR_SHORT_SCALES, "FigS_PixelBasin_Spearman_1_3_6"))
.step("Spearman 12/24/36", make_corr_figure(c("spi", "spei", "spei_era5"), CORR_LONG_SCALES, "FigS_PixelBasin_Spearman_12_24_36"))
.step("Spearman SRI short", make_corr_figure("sri", CORR_SHORT_SCALES, "FigS_PixelBasin_Spearman_SRI_1_3_6"))
.step("Spearman SRI long",  make_corr_figure("sri", CORR_LONG_SCALES, "FigS_PixelBasin_Spearman_SRI_12_24_36"))
.step("Spearman snow", make_snow_corr_figure())

################################################################################
# PART 9 – REGIONAL DOMAIN MAP (regional trend tables are computed in Script 11)
################################################################################
cat("\n══ PART 9: regional domain map ══\n")
regional_summary_dir <- file.path(TREND_DIR, "regional_spatial_summary")
dir.create(regional_summary_dir, recursive = TRUE, showWarnings = FALSE)
fig_domain_map <- function() {
  sub_regions <- load_nechako_regions()
  if (!BASIN_RID %in% names(sub_regions)) stop("Whole-basin region '", BASIN_RID, "' not found in load_nechako_regions().")
  reg_sf <- sf::st_as_sf(sub_regions[[BASIN_RID]])
  sub_ids <- setdiff(names(sub_regions), BASIN_RID)
  if (!length(sub_ids)) stop("No sub-basins remain after removing the whole-basin union.")
  sub_list <- sub_regions[sub_ids]
  stopifnot(all(vapply(sub_list, inherits, logical(1), "SpatVector")))
  sub_sf <- sf::st_as_sf(do.call(rbind, unname(sub_list)))
  sub_sf$region_id <- sub_ids; sub_sf$region_name <- gsub("_", " ", sub_ids, fixed = TRUE)
  pmap <- ggplot() + geom_sf(data = sub_sf, fill = NA, linewidth = 0.6) +
    geom_sf_text(data = sub_sf, aes(label = region_name), size = 2.8) +
    geom_sf(data = reg_sf, fill = NA, linewidth = 1.0) + theme_bw(base_size = 10) +
    labs(title = "Nechako regional analysis domain", subtitle = sprintf("%d sub-basins; whole basin = dissolved union", length(sub_ids)))
  ggplot2::ggsave(file.path(regional_summary_dir, "Nechako_subbasin_domain.png"), pmap, width = 10, height = 8, dpi = DPI)
  ggplot2::ggsave(file.path(regional_summary_dir, "Nechako_subbasin_domain.pdf"), pmap, width = 10, height = 8)
  cat("  ✓ Nechako_subbasin_domain.png / .pdf\n")
}
.step("domain map", fig_domain_map())
for (f in c("regional_trend_full_series.csv", "regional_trend_calendar_month.csv"))
  cat(sprintf("  %s %s\n", if (file.exists(file.path(TREND_DIR, "regional_trends", f))) "✓" else "⚠ missing:", f))

################################################################################
# FINAL OUTPUT VALIDATION
################################################################################
verify_required_file <- function(path, label) {
  if (!file.exists(path)) stop("FINAL OUTPUT VALIDATION FAILED: missing ", label, ": ", path)
  sz <- file.info(path)$size
  if (is.na(sz) || sz <= 0) stop("FINAL OUTPUT VALIDATION FAILED: empty ", label, ": ", path)
  cat(sprintf("  ✓ %-40s %s (%s bytes)\n", label, basename(path), format(sz, big.mark = ",", scientific = FALSE)))
}
cat("\n============================================================\nFINAL OUTPUT VALIDATION\n",
    "============================================================\n", sep = "")
verify_required_file(master_csv, "Master diagnostics CSV")
verify_required_file(ms_fig_pdf, "Manuscript Figure X PDF")
verify_required_file(ms_fig_png, "Manuscript Figure X PNG")
if (!is.null(cm_data)) {
  verify_required_file(file.path(TREND_DIR, "all_calendar_month_results.csv"), "Calendar-month master CSV")
  verify_required_file(cm_fig_pdf, "Calendar-month manuscript figure PDF")
  verify_required_file(sub("\\.pdf$", ".png", cm_fig_pdf), "Calendar-month manuscript figure PNG")
}
if (nrow(snow_panels)) {
  verify_required_file(snow_fig_pdf, "Snow drought Sen-slope PDF")
  verify_required_file(snow_fig_png, "Snow drought Sen-slope PNG")
}
if (nrow(runoff_panels)) {
  verify_required_file(runoff_fig_pdf, "Runoff drought Sen-slope PDF")
  verify_required_file(runoff_fig_png, "Runoff drought Sen-slope PNG")
}
verify_required_file(file.path(CORR_OUT_DIR, "PixelBasin_Spearman_Summary.csv"), "Pixel-basin Spearman summary")

inv <- data.table::data.table(file = list.files(fig_dir, recursive = TRUE, pattern = "\\.(png|pdf|csv)$"))
inv[, folder := dirname(file)]
data.table::fwrite(inv, file.path(fig_dir, "Script8_output_inventory.csv"))
cat("\n  Output inventory (files per folder):\n")
print(as.data.frame(inv[, .N, by = folder][order(folder)]))

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  FINAL OUTPUT VALIDATION PASSED — 13trends_visualization.R v2 DONE  ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n")
cat("  Figures : ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat(sprintf("  Figure X method: %s (companion: %s) | calendar-month method: %s\n",
            .mcol(MS_FULL_METHOD)$lab, .mcol(alt_method)$lab, .mcol(CALMONTH_METHOD)$lab))
cat("  Check the console above for any ❌ or ⚠ lines.\n\n")