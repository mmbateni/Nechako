# w4_trends_visualization.R  ·  SPATIAL DROUGHT DIAGNOSTICS
# ─────────────────────────────────────────────────────────────────────────────────
# COMPATIBLE WITH: DROUGHT_ANALYSIS_utils.R, w1_trend_test.R
# SUPPORTS: SPI, SPEI, SWEI, MSPI, MSPEI (Gridded)
#
# SCOPE (spatial figures only):
#   PART 1 – Merge individual result CSVs → master diagnostics CSV
#   PART 2 – Spatial figures
#     Fig 1: Trend maps          (Kendall τ_vc + FDR significance)
#     Fig 2: Event maps          (# events, mean duration, max intensity)
#     Fig 3: Temporal pattern maps (clustering, regime-shift decade)
#     Fig 4: Timescale comparison (% significant, median τ, events, duration)
#     Fig 5: Method comparison   (VC vs TFPW scatter + agreement bar)
#
# NOTE: Basin-averaged time series CSVs, PNG plots, and the Excel summary
#       are produced by w1_basin_timeseries.R — run that script for all
#       non-spatial outputs.
#
# EXECUTION ORDER: w1_trend_test.R → w2_basin_timeseries.R → w4_trends_visualization.R

source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "data.table", "ggplot2", "sf", "viridis",
                      "cowplot", "patchwork", "scales", "RColorBrewer",
                      "ggspatial", "gridExtra", "grid", "lubridate"))

# ── Working directory & output folder ───────────────────────────────
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
fig_dir <- file.path(TREND_DIR, "figures_temporal_diagnostics")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Basin boundary ──────────────────────────────────────────────────
basin_sf   <- load_basin(BASIN_SHP, EQUAL_AREA_CRS)
basin_vect <- terra::vect(basin_sf)
basin_sf_global <<- basin_sf   # expose to plot_raster_clean() in utils

# ── Analysis constants ──────────────────────────────────────────────
DPI              <- 300
TIMESCALES_VIZ   <- TIMESCALES_SPATIAL   # c(1,3,6,12) for spatial figures
MSPI_MSPEI_SCALE <- 1L                   # single-scale multivariate indices

# ── Index colours (extend for MSPI / MSPEI) ─────────────────────────
if (!"swei"  %in% names(index_colours))
  index_colours <- c(index_colours, swei  = "#1E90FF")
if (!"mspi"  %in% names(index_colours))
  index_colours <- c(index_colours, mspi  = "#33a02c")   # green
if (!"mspei" %in% names(index_colours))
  index_colours <- c(index_colours, mspei = "#fb9a99")   # light red

################################################################################
# PART 1 – MERGE INDIVIDUAL RESULT FILES → MASTER CSV
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 1: Building master diagnostics file\n")
cat("══════════════════════════════════════\n")

master_csv <- file.path(TREND_DIR, "all_temporal_diagnostics_results.csv")

load_clip_results <- function(index_type, timescale) {
  f <- file.path(TREND_DIR, sprintf("%s_%02d_results.csv", index_type, timescale))
  if (!file.exists(f)) {
    if (!(index_type %in% c("mspi", "mspei") ||
          (index_type == "swei" && timescale != SWEI_SCALE)))
      cat("  ⚠ Missing: ", basename(f), "\n")
    return(NULL)
  }
  dt <- data.table::fread(f)
  cat(sprintf("  %s-%02d: %d rows\n", toupper(index_type), timescale, nrow(dt)))
  
  # Standardise coordinate column names
  if ("x" %in% names(dt) && !"lon" %in% names(dt))
    data.table::setnames(dt, "x", "lon")
  if ("y" %in% names(dt) && !"lat" %in% names(dt))
    data.table::setnames(dt, "y", "lat")
  if (!all(c("lon", "lat") %in% names(dt))) {
    cat("  ⚠ No coordinate columns — skipping\n"); return(NULL)
  }
  dt <- dt[!is.na(lon) & !is.na(lat)]
  
  tryCatch(
    clip_to_basin(dt, basin_sf, EQUAL_AREA_CRS),
    error = function(e) { cat("  ⚠ Clip error: ", e$message, "\n"); NULL }
  ) -> dt_clip
  
  if (is.null(dt_clip) || nrow(dt_clip) == 0) {
    cat("  ⚠ No points remain after basin clip\n"); return(NULL)
  }
  dt_clip[, `:=`(index_type = index_type, timescale = timescale)]
  dt_clip
}

all_results <- list()

cat("\n  Loading SPI & SPEI...\n")
for (idx in c("spi", "spei")) {
  for (sc in TIMESCALES_VIZ) {
    key <- sprintf("%s_%02d", idx, sc)
    r   <- load_clip_results(idx, sc)
    if (!is.null(r)) all_results[[key]] <- r
  }
}

cat("\n  Loading SWEI (scale", SWEI_SCALE, ")...\n")
swei_res <- load_clip_results("swei", SWEI_SCALE)
if (!is.null(swei_res)) all_results[[sprintf("swei_%02d", SWEI_SCALE)]] <- swei_res

cat("\n  Loading MSPI & MSPEI (scale", MSPI_MSPEI_SCALE, ")...\n")
for (idx in c("mspi", "mspei")) {
  key <- sprintf("%s_%02d", idx, MSPI_MSPEI_SCALE)
  r   <- load_clip_results(idx, MSPI_MSPEI_SCALE)
  if (!is.null(r)) all_results[[key]] <- r
}

if (!length(all_results))
  stop("No data loaded. Run w1_trend_test.R first.")

combined <- data.table::rbindlist(all_results, fill = TRUE)
data.table::fwrite(combined, master_csv)
cat(sprintf("✓ Master CSV: %d rows → %s\n", nrow(combined), basename(master_csv)))

################################################################################
# PART 2 – SPATIAL FIGURES
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 2: Spatial figures\n")
cat("══════════════════════════════════════\n")

plot_data <- combined

# ── Shared raster helpers ──────────────────────────────────────────
make_raster <- function(dt, vcol, tpl = NULL) {
  if (is.null(tpl)) tpl <- create_raster_template(dt, basin_sf)
  create_raster_from_points(dt, tpl, vcol, basin_sf)
}

# Which indices are available for a given scale
active_indices <- function(scale) {
  base <- c("spi", "spei")
  if (scale == SWEI_SCALE        && any(plot_data$index_type == "swei"))
    base <- c(base, "swei")
  if (scale == MSPI_MSPEI_SCALE) {
    if (any(plot_data$index_type == "mspi"))  base <- c(base, "mspi")
    if (any(plot_data$index_type == "mspei")) base <- c(base, "mspei")
  }
  base
}

# ── Fig 1: Trend maps (Kendall τ_vc + FDR significance) ───────────
# Layout — 2 rows × n_index columns:
#   Row 1  Kendall's τ (VC-HR98) with stippling overlay
#          • Filled dots mark pixels that survive BH-FDR at α = 0.05
#          • Stippling follows Wilks (2006, BAMS): show all τ values but
#            identify only the locally robust subset
#   Row 2  Three-category significance map:
#          • Dark red   — FDR-significant          (p_fdr_vc  < 0.05)
#          • Orange     — Nominally-sig only        (p_value_vc < 0.05 but p_fdr_vc ≥ 0.05)
#          • Light grey — Not significant           (p_value_vc ≥ 0.05)
#   Falls back to 4-level nominal panel when p_fdr_vc is absent (old w3 run).
create_trend_maps <- function(scale = 12, outfile) {
  cat(sprintf("  Fig 1: trend maps (scale=%d)... ", scale))
  if (!safe_pdf(outfile, 12, 10)) return()
  indices_here <- active_indices(scale)
  n_plots      <- length(indices_here)
  
  par(mfrow = c(2, n_plots), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  tau_range <- c(-0.4, 0.4)
  has_fdr   <- "p_fdr_vc" %in% names(plot_data)
  if (!has_fdr)
    cat("\n  ⚠ p_fdr_vc absent — falling back to nominal significance",
        "\n    (re-run w1_trend_test.R to generate FDR-corrected columns)\n")
  
  for (idx in indices_here) {
    sub <- plot_data[index_type == idx & timescale == scale]
    if (!nrow(sub)) {
      plot.new(); text(0.5, 0.5, paste("No data:", toupper(idx))); next
    }
    tpl   <- create_raster_template(sub, basin_sf)
    r_tau <- make_raster(sub, "tau_vc", tpl)
    
    # Row 1: τ map + FDR stippling
    plot_raster_clean(r_tau, paste0(toupper(idx), ": Kendall's τ (VC-HR98)"),
                      tau_range, grDevices::hcl.colors(101, "RdBu", rev = TRUE),
                      legend_title = "τ")
    if (has_fdr) {
      stip <- sub[!is.na(p_fdr_vc) & p_fdr_vc < 0.05]
      if (nrow(stip) > 0)
        graphics::points(stip$lon, stip$lat, pch = 20, cex = 0.45, col = "black")
      graphics::legend("bottomleft",
                       legend = "FDR-significant (BH α = 0.05)",
                       pch = 20, col = "black", bty = "n", cex = 0.75)
    }
    
    # Row 2: significance category map
    sub2 <- data.table::copy(sub)
    if (has_fdr) {
      sub2[, sig_cat := data.table::fcase(
        is.na(p_value_vc),                         NA_real_,
        !is.na(p_fdr_vc) & p_fdr_vc  < 0.05,      3,   # FDR-significant
        p_value_vc < 0.05,                          2,   # nominally-sig only
        default =                                   1)]  # not significant
      r_sig <- make_raster(sub2[!is.na(sig_cat)], "sig_cat", tpl)
      plot_raster_clean(
        r_sig,
        paste0(toupper(idx), ": Significance (BH-FDR corrected)"),
        NULL, c("#f7f7f7", "#fc8d59", "#d73027"),
        breaks = c(0.5, 1.5, 2.5, 3.5), categorical = TRUE, legend = FALSE)
      graphics::legend("bottomright",
                       legend = c("Not significant",
                                  "Nominal only (p < 0.05)",
                                  "FDR-significant"),
                       fill = c("#f7f7f7", "#fc8d59", "#d73027"),
                       bty = "n", cex = 0.75, title = "Significance")
    } else {
      sub2[, sig_num := dplyr::case_when(
        is.na(p_value_vc)   ~ NA_real_,
        p_value_vc < 0.001  ~ 1,
        p_value_vc < 0.01   ~ 2,
        p_value_vc < 0.05   ~ 3,
        TRUE                ~ 4)]
      r_sig <- make_raster(sub2[!is.na(sig_num)], "sig_num", tpl)
      plot_raster_clean(r_sig, paste0(toupper(idx), ": Significance (nominal)"),
                        NULL, c("#d73027", "#fc8d59", "#fee08b", "#f7f7f7"),
                        breaks = c(0.5, 1.5, 2.5, 3.5, 4.5),
                        categorical = TRUE, legend_title = "p-value")
    }
  }
  graphics::mtext(sprintf("Trend Analysis – VC-HR98 Method (%d-month)", scale),
                  outer = TRUE, cex = 1.2, font = 2)
  grDevices::dev.off()
  cat(" ✓\n")
}

# ── Fig 2: Event characteristic maps ──────────────────────────────
# Colour design: DARKER colour = MORE extreme drought in every panel.
#   n_events      YlOrRd  (dark at high end, rev = FALSE)
#   mean_duration viridis reversed (long duration = dark, rev = TRUE)
#   max_intensity Reds reversed (most-negative = darkest, rev = TRUE)
create_event_maps <- function(scale = 12, outfile) {
  cat(sprintf("  Fig 2: event maps (scale=%d)... ", scale))
  if (!safe_pdf(outfile, 12, 14)) return()
  indices_here <- active_indices(scale)
  n_plots      <- length(indices_here)
  par(mfrow = c(3, n_plots), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  
  metrics <- list(
    n_events      = list(title = "# Events",      pal = "YlOrRd"),
    mean_duration = list(title = "Mean Duration",  pal = "viridis"),
    max_intensity = list(title = "Max Intensity",  pal = "Reds"))
  
  for (idx in indices_here) {
    sub <- plot_data[index_type == idx & timescale == scale]
    if (!nrow(sub)) {
      for (i in 1:3) { plot.new(); text(0.5, 0.5, paste("No data:", toupper(idx))) }
      next
    }
    tpl <- create_raster_template(sub, basin_sf)
    for (mn in names(metrics)) {
      if (all(is.na(sub[[mn]]))) {
        plot.new(); text(0.5, 0.5, paste("All NA:", mn)); next
      }
      r   <- make_raster(sub, mn, tpl)
      zlm <- if (mn == "max_intensity") c(min(sub[[mn]], na.rm = TRUE), 0)
      else                       c(0, max(sub[[mn]], na.rm = TRUE))
      plot_raster_clean(
        r, paste0(toupper(idx), ": ", metrics[[mn]]$title), zlm,
        grDevices::hcl.colors(101, metrics[[mn]]$pal,
                              rev = mn %in% c("max_intensity", "mean_duration")))
    }
  }
  graphics::mtext(sprintf("Event Characteristics (%d-month)", scale),
                  outer = TRUE, cex = 1.2, font = 2)
  grDevices::dev.off()
  cat(" ✓\n")
}

# ── Fig 3: Temporal pattern maps ──────────────────────────────────
create_temporal_maps <- function(scale = 12, outfile) {
  cat(sprintf("  Fig 3: temporal maps (scale=%d)...", scale))
  if (!safe_pdf(outfile, 12, 8)) return()
  indices_here <- active_indices(scale)
  n_plots      <- length(indices_here)
  par(mfrow = c(n_plots, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  need_cols <- c("filtered_runs", "p_value_runs", "clustering",
                 "regime_shift_year", "n_spectral_peaks")
  missing_cols <- setdiff(need_cols, names(plot_data))
  if (length(missing_cols)) {
    plot.new()
    text(0.5, 0.5,
         paste("Missing columns – re-run w1_trend_test.R:\n",
               paste(missing_cols, collapse = ", ")),
         cex = 1.1, col = "red")
    grDevices::dev.off()
    cat(" ⚠ missing columns (need w3 upgrade)\n")
    return()
  }
  
  for (idx in indices_here) {
    sub <- plot_data[index_type == idx & timescale == scale]
    tpl <- create_raster_template(sub, basin_sf)
    if (is.null(tpl)) {
      for (i in 1:3) { plot.new(); text(0.5, 0.5, paste("No data:", toupper(idx))) }
      next
    }
    
    # Panel A: Temporal clustering (runs test)
    sub_cl <- data.table::copy(sub)
    sub_cl[, clus_num := data.table::fcase(
      filtered_runs | is.na(p_value_runs), NA_real_,
      p_value_runs >= 0.05,                0,
      clustering == "clustered",           -1,
      clustering == "dispersed",            1)]
    r_cl <- make_raster(sub_cl, "clus_num", tpl)
    plot_raster_clean(r_cl, paste0(toupper(idx), ": Temporal Clustering"),
                      NULL, c("#e41a1c", "#fee090", "#4daf4a"),
                      breaks = c(-1.5, -0.5, 0.5, 1.5),
                      categorical = TRUE, legend = FALSE)
    graphics::legend("bottomright",
                     legend = c("Clustered", "Not sig.", "Dispersed"),
                     fill = c("#e41a1c", "#fee090", "#4daf4a"),
                     bty = "n", cex = 0.85)
    
    # Panel B: Regime-shift decade
    sub_rs <- sub[!is.na(regime_shift_year)]
    if (nrow(sub_rs)) {
      sub_rs2 <- data.table::copy(sub_rs)
      sub_rs2[, shift_dec := as.numeric(
        cut(regime_shift_year,
            breaks = seq(1950, 2030, by = 10),
            labels = 1:8, include.lowest = TRUE))]
      r_rs <- make_raster(sub_rs2, "shift_dec", tpl)
      plot_raster_clean(r_rs, paste0(toupper(idx), ": Regime-Shift Decade"),
                        c(1, 8), grDevices::hcl.colors(8, "Spectral", rev = TRUE),
                        legend_title = "Decade")
    } else {
      plot.new()
      text(0.5, 0.5, paste(toupper(idx), "– no regime shifts detected"))
    }
    
    # Panel C: Spectral peaks
    r_sp <- make_raster(sub, "n_spectral_peaks", tpl)
    plot_raster_clean(r_sp, paste0(toupper(idx), ": Spectral Peaks"),
                      c(0, max(sub$n_spectral_peaks, na.rm = TRUE, 1)),
                      grDevices::hcl.colors(101, "viridis"))
  }
  graphics::mtext(sprintf("Temporal Patterns (%d-month)", scale),
                  outer = TRUE, cex = 1.2, font = 2)
  grDevices::dev.off()
  cat(" ✓\n")
}

# ── Fig 4: Timescale comparison ────────────────────────────────────
# Panel A: TWO lines per index × scale — solid = nominal sig, dashed = BH-FDR.
# Gap between lines = false-discovery inflation at that scale.
create_timescale_comparison <- function(outfile) {
  cat("  Fig 4: timescale comparison... ")
  has_fdr <- "p_fdr_vc" %in% names(plot_data)
  
  smry <- plot_data[, .(
    Pct_Sig_Nom = sum(p_value_vc < 0.05, na.rm = TRUE) / .N * 100,
    Pct_Sig_FDR = if (has_fdr)
      sum(!is.na(p_fdr_vc) & p_fdr_vc < 0.05, na.rm = TRUE) / .N * 100
    else NA_real_,
    Median_Tau  = median(tau_vc,       na.rm = TRUE),
    Mean_Events = mean(n_events,        na.rm = TRUE),
    Mean_Dur    = mean(mean_duration,   na.rm = TRUE)
  ), by = .(index_type, timescale)]
  
  mk_single <- function(col, ylab, title) {
    ggplot2::ggplot(smry,
                    ggplot2::aes(x = timescale, y = .data[[col]],
                                 color = index_type, group = index_type)) +
      ggplot2::geom_line(linewidth = 1.2, na.rm = TRUE) +
      ggplot2::geom_point(size = 4, na.rm = TRUE) +
      ggplot2::scale_color_manual(values = index_colours, labels = toupper) +
      ggplot2::scale_x_continuous(breaks = sort(unique(smry$timescale))) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = title,
                    x = "Timescale (months)", y = ylab, color = "Index")
  }
  
  sig_long <- data.table::melt(
    smry,
    id.vars      = c("index_type", "timescale"),
    measure.vars = c("Pct_Sig_Nom", "Pct_Sig_FDR"),
    variable.name = "Correction", value.name = "Pct_Sig")
  sig_long[, Correction := data.table::fcase(
    Correction == "Pct_Sig_Nom", "Nominal (p < 0.05)",
    Correction == "Pct_Sig_FDR", "BH-FDR (p_fdr < 0.05)")]
  sig_long <- sig_long[!is.na(Pct_Sig)]
  
  pA <- ggplot2::ggplot(
    sig_long,
    ggplot2::aes(x        = timescale, y        = Pct_Sig,
                 color    = index_type,
                 linetype = Correction,
                 group    = interaction(index_type, Correction))) +
    ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
    ggplot2::geom_point(
      data  = sig_long[Correction == "Nominal (p < 0.05)"],
      size  = 3.5, na.rm = TRUE) +
    ggplot2::geom_point(
      data  = sig_long[Correction == "BH-FDR (p_fdr < 0.05)"],
      size  = 3.5, shape = 1, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = index_colours, labels = toupper) +
    ggplot2::scale_linetype_manual(
      values = c("Nominal (p < 0.05)"    = "solid",
                 "BH-FDR (p_fdr < 0.05)" = "dashed")) +
    ggplot2::scale_x_continuous(breaks = sort(unique(smry$timescale))) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title    = "A) % Significant Pixels — Nominal vs BH-FDR",
      subtitle = "Solid = nominal  |  Dashed = BH-FDR corrected  |  Gap = false-discovery inflation",
      x = "Timescale (months)", y = "% Significant pixels",
      color = "Index", linetype = "Correction")
  
  p <- (pA | mk_single("Median_Tau",  "τ",        "B) Median Kendall's τ (VC-HR98)")) /
    (mk_single("Mean_Events", "N events",  "C) Mean Drought Events") |
       mk_single("Mean_Dur",    "months",   "D) Mean Duration")) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  ggplot2::ggsave(outfile, p, width = 14, height = 10, dpi = DPI)
  cat(" ✓\n")
}

# ── Fig 5: Method comparison (VC vs TFPW) ─────────────────────────
# Panel A: scatter τ_vc vs τ_tfpw (tau is not affected by FDR correction)
# Panel B: agreement bar; uses BH-FDR significance when columns are present,
#          falls back to nominal p-values with a labelled subtitle if absent.
create_method_comparison <- function(outfile) {
  cat("  Fig 5: method comparison (VC vs TFPW)... ")
  needed <- c("tau_tfpw", "p_value_tfpw")
  if (!all(needed %in% names(plot_data))) {
    if (safe_pdf(outfile)) {
      plot.new()
      text(0.5, 0.5,
           "TFPW columns not found.\nRe-run w1_trend_test.R with TFPW enabled.",
           cex = 1.2, col = "orange")
      grDevices::dev.off()
    }
    cat(" ⚠ TFPW columns absent — placeholder saved\n")
    return()
  }
  
  has_fdr <- all(c("p_fdr_vc", "p_fdr_tfpw") %in% names(plot_data))
  sub     <- plot_data[timescale == 12]
  
  if (has_fdr) {
    sub[, `:=`(vc_sig   = !is.na(p_fdr_vc)   & p_fdr_vc   < 0.05,
               tfpw_sig = !is.na(p_fdr_tfpw)  & p_fdr_tfpw < 0.05)]
    bar_subtitle <- "BH-FDR corrected (p_fdr < 0.05)"
  } else {
    sub[, `:=`(vc_sig   = !is.na(p_value_vc)   & p_value_vc   < 0.05,
               tfpw_sig = !is.na(p_value_tfpw)  & p_value_tfpw < 0.05)]
    bar_subtitle <- "nominal p < 0.05 (re-run w3 for FDR columns)"
  }
  
  p1 <- ggplot2::ggplot(sub, ggplot2::aes(tau_vc, tau_tfpw, color = index_type)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "gray40") +
    ggplot2::geom_point(alpha = 0.4, size = 2, na.rm = TRUE) +
    ggplot2::facet_wrap(~index_type,
                        labeller = ggplot2::labeller(index_type = toupper)) +
    ggplot2::scale_color_manual(values = index_colours) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "A) τ: Variance-Corrected (VC-HR98) vs TFPW (12-month)",
                  x = "τ (VC-HR98)", y = "τ (TFPW)")
  
  agree <- sub[, .(
    Both      = sum( vc_sig &  tfpw_sig, na.rm = TRUE),
    VC_only   = sum( vc_sig & !tfpw_sig, na.rm = TRUE),
    TFPW_only = sum(!vc_sig &  tfpw_sig, na.rm = TRUE),
    None      = sum(!vc_sig & !tfpw_sig, na.rm = TRUE)
  ), by = index_type]
  ag_l <- data.table::melt(agree, id.vars = "index_type",
                           variable.name = "Category", value.name = "N_pixels")
  
  p2 <- ggplot2::ggplot(ag_l,
                        ggplot2::aes(toupper(index_type), N_pixels, fill = Category)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title    = "B) Significance Agreement — VC-HR98 vs TFPW (12-month)",
      subtitle = bar_subtitle,
      x = "Index", y = "# Pixels")
  
  ggplot2::ggsave(outfile, p1 | p2, width = 14, height = 6, dpi = DPI)
  cat(" ✓\n")
}

# ── Run all spatial figures ─────────────────────────────────────────
for (sc in c(1, 3, 6, 12)) {
  sfx <- sprintf("%02dmo", sc)
  tryCatch(
    create_trend_maps(sc, file.path(fig_dir, sprintf("Fig1_TrendMaps_%s.pdf", sfx))),
    error = function(e) cat("  ❌ Fig 1 trend maps:", e$message, "\n"))
  tryCatch(
    create_event_maps(sc, file.path(fig_dir, sprintf("Fig2_EventMaps_%s.pdf", sfx))),
    error = function(e) cat("  ❌ Fig 2 event maps:", e$message, "\n"))
  tryCatch(
    create_temporal_maps(sc, file.path(fig_dir, sprintf("Fig3_TemporalMaps_%s.pdf", sfx))),
    error = function(e) cat("  ❌ Fig 3 temporal maps:", e$message, "\n"))
}

tryCatch(
  create_timescale_comparison(file.path(fig_dir, "Fig4_TimescaleComparison.png")),
  error = function(e) cat("  ❌ Fig 4 timescale comparison:", e$message, "\n"))
tryCatch(
  create_method_comparison(file.path(fig_dir, "Fig5_MethodComparison.png")),
  error = function(e) cat("  ❌ Fig 5 method comparison:", e$message, "\n"))

################################################################################
# PART 3 – MANUSCRIPT FIGURE 8
# Two-panel publication-quality spatial trend map for SPEI-1 and SPEI-3.
#
# Panel layout: 1 row × 2 columns (one column per scale).
# Each column contains two sub-panels stacked vertically:
#   (top)    Kendall τ_vc choropleth on a diverging RdBu palette
#            Stippling: filled dots on FDR-significant pixels (BH α = 0.05)
#   (bottom) Three-category significance map
#            Dark red   = FDR-significant (p_fdr_vc < 0.05)
#            Orange     = Nominal only    (p_value_vc < 0.05, p_fdr_vc ≥ 0.05)
#            Light grey = Not significant
#
# Falls back gracefully when:
#   • p_fdr_vc column absent → significance map uses nominal p only
#   • A requested scale is not in the data → panel shows a labelled placeholder
#
# Outputs (fig_dir):
#   Fig8_MS_TrendMap_SPEI1_SPEI3.pdf   (vector, 7.0 × 8.5 in)
#   Fig8_MS_TrendMap_SPEI1_SPEI3.png   (300 DPI)
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 3: Manuscript Figure 8 — SPEI-1 & SPEI-3 trend maps\n")
cat("══════════════════════════════════════\n")

create_ms_trend_maps <- function(scales_ms = c(1, 3),
                                 index_ms  = "spei",
                                 outpdf,
                                 outpng,
                                 dpi = 300) {
  
  has_fdr <- "p_fdr_vc" %in% names(plot_data)
  if (!has_fdr)
    cat("  NOTE: p_fdr_vc absent — significance map uses nominal p-values only\n")
  
  # ── shared style ────────────────────────────────────────────────────────────
  theme_map <- ggplot2::theme_void(base_size = 8.5) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 8.5, face = "bold",
                                            hjust = 0.5,
                                            margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = 7, colour = "grey35",
                                            hjust = 0.5,
                                            margin = ggplot2::margin(b = 2)),
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.38, "cm"),
      legend.text      = ggplot2::element_text(size = 6.5),
      legend.title     = ggplot2::element_text(size = 7, face = "bold"),
      plot.margin      = ggplot2::margin(4, 6, 4, 6)
    )
  
  # colour scales
  TAU_BREAKS  <- c(-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)
  TAU_PALETTE <- rev(RColorBrewer::brewer.pal(9, "RdBu"))
  
  SIG_FILLS <- c(
    "FDR-significant"    = "#d73027",
    "Nominal only"       = "#fc8d59",
    "Not significant"    = "#f7f7f7"
  )
  SIG_COLS  <- c(
    "FDR-significant"    = "#a50f15",
    "Nominal only"       = "#d94701",
    "Not significant"    = "#bdbdbd"
  )
  
  # ── helper: raster → ggplot-friendly data frame ─────────────────────────────
  rast_to_df <- function(r) {
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(df)[3] <- "value"
    df
  }
  
  # ── helper: basin sf in WGS84 for ggplot2 layer ─────────────────────────────
  basin_wgs <- tryCatch(
    sf::st_transform(basin_sf, crs = 4326),
    error = function(e) basin_sf)
  
  # ── build one column (two sub-panels) per scale ────────────────────────────
  col_panels <- lapply(seq_along(scales_ms), function(si) {
    sc  <- scales_ms[si]
    lab <- letters[si]    # "(a)", "(b)", ...
    
    sub <- plot_data[index_type == index_ms & timescale == sc]
    if (!nrow(sub)) {
      empty <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = sprintf("No data\n%s-%d",
                                          toupper(index_ms), sc),
                          size = 4, colour = "grey50") +
        ggplot2::theme_void()
      return(list(tau = empty, sig = empty))
    }
    
    tpl   <- create_raster_template(sub, basin_sf)
    r_tau <- make_raster(sub, "tau_vc", tpl)
    
    # ── tau choropleth ────────────────────────────────────────────────────────
    df_tau <- rast_to_df(r_tau)
    df_tau$value <- pmax(-0.4, pmin(0.4, df_tau$value))
    
    # FDR-significant stip points
    stip <- if (has_fdr)
      sub[!is.na(p_fdr_vc) & p_fdr_vc < 0.05]
    else
      sub[!is.na(p_value_vc) & p_value_vc < 0.05]
    
    # basin stats for subtitle
    n_sig   <- nrow(stip)
    n_total <- nrow(sub[!is.na(tau_vc)])
    pct_sig <- round(100 * n_sig / max(n_total, 1), 1)
    med_tau <- round(median(sub$tau_vc, na.rm = TRUE), 3)
    sig_lab <- if (has_fdr) "FDR-sig" else "nom-sig"
    
    p_tau <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = df_tau,
        ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::scale_fill_gradientn(
        colours = TAU_PALETTE,
        limits  = c(-0.4, 0.4),
        breaks  = TAU_BREAKS,
        name    = "\u03c4 (VC-HR98)",
        guide   = ggplot2::guide_colorbar(
          barwidth = 6, barheight = 0.45,
          title.position = "top")) +
      # stippling: FDR-significant dots
      {if (nrow(stip) > 0)
        ggplot2::geom_point(
          data = as.data.frame(stip),
          ggplot2::aes(x = lon, y = lat),
          shape = 20, size = 0.55,
          colour = "black", alpha = 0.75,
          inherit.aes = FALSE)
        else NULL} +
      ggplot2::geom_sf(
        data        = basin_wgs,
        fill        = NA,
        colour      = "black",
        linewidth   = 0.45,
        inherit.aes = FALSE) +
      ggplot2::coord_sf(expand = FALSE) +
      ggplot2::labs(
        title    = sprintf("(%s)  %s-%d \u2014 Kendall\u2019s \u03c4 (VC-HR98)",
                           lab, toupper(index_ms), sc),
        subtitle = sprintf("Median \u03c4 = %+.3f  |  %s: %d/%d pixels (%.1f%%)",
                           med_tau, sig_lab, n_sig, n_total, pct_sig)) +
      theme_map +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text       = ggplot2::element_blank())
    
    # ── significance category map ─────────────────────────────────────────────
    sub2 <- data.table::copy(sub)
    if (has_fdr) {
      sub2[, sig_cat := data.table::fcase(
        is.na(p_value_vc),                        NA_character_,
        !is.na(p_fdr_vc) & p_fdr_vc  < 0.05,     "FDR-significant",
        p_value_vc < 0.05,                         "Nominal only",
        default =                                  "Not significant"
      )]
      sig_note <- "FDR = BH correction (\u03b1 = 0.05)"
    } else {
      sub2[, sig_cat := data.table::fcase(
        is.na(p_value_vc),   NA_character_,
        p_value_vc < 0.05,   "Nominal only",
        default =            "Not significant"
      )]
      sig_note <- "Nominal p < 0.05 (re-run w1 for FDR columns)"
    }
    
    sub2_valid <- sub2[!is.na(sig_cat)]
    r_sig  <- make_raster(sub2_valid, "sig_cat", tpl)
    
    # convert categorical raster to df with string labels
    df_sig_raw <- as.data.frame(r_sig, xy = TRUE, na.rm = TRUE)
    names(df_sig_raw)[3] <- "value"
    # the raster stores integers; map back to factor labels
    lev_map <- c("1" = "FDR-significant",
                 "2" = "Nominal only",
                 "3" = "Not significant")
    if (!has_fdr) lev_map <- c("1" = "Nominal only",
                               "2" = "Not significant")
    df_sig_raw$sig_cat <- factor(
      lev_map[as.character(round(df_sig_raw$value))],
      levels = names(SIG_FILLS))
    df_sig <- df_sig_raw[!is.na(df_sig_raw$sig_cat), ]
    
    p_sig <- ggplot2::ggplot() +
      ggplot2::geom_raster(
        data = df_sig,
        ggplot2::aes(x = x, y = y, fill = sig_cat)) +
      ggplot2::scale_fill_manual(
        values = SIG_FILLS,
        name   = "Significance",
        drop   = FALSE,
        guide  = ggplot2::guide_legend(
          title.position = "top",
          nrow = 1,
          override.aes = list(colour = SIG_COLS))) +
      ggplot2::geom_sf(
        data        = basin_wgs,
        fill        = NA,
        colour      = "black",
        linewidth   = 0.45,
        inherit.aes = FALSE) +
      ggplot2::coord_sf(expand = FALSE) +
      ggplot2::labs(
        title    = sprintf("(%s\u2019)  %s-%d \u2014 Significance",
                           lab, toupper(index_ms), sc),
        subtitle = sig_note) +
      theme_map +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text       = ggplot2::element_blank())
    
    list(tau = p_tau, sig = p_sig)
  })
  
  # ── assemble: tau panels top row, sig panels bottom row ─────────────────────
  # order: (a)_tau | (b)_tau
  #        (a')_sig | (b')_sig
  tau_row <- patchwork::wrap_plots(
    lapply(col_panels, `[[`, "tau"), nrow = 1)
  sig_row <- patchwork::wrap_plots(
    lapply(col_panels, `[[`, "sig"), nrow = 1)
  
  fig8 <- tau_row / sig_row +
    patchwork::plot_layout(heights = c(1, 1)) +
    patchwork::plot_annotation(
      title    = sprintf(
        "Spatial pattern of drought trends \u2014 Nechako River Basin (1950\u20132025)"),
      subtitle = paste0(
        toupper(index_ms), " at ",
        paste(scales_ms, collapse = "- and "),
        "-month accumulation scales.  ",
        "Top row: Kendall\u2019s \u03c4 (variance-corrected HR98); ",
        "stippling = ",
        if (has_fdr) "FDR-significant pixels (BH \u03b1 = 0.05)."
        else         "nominally significant pixels (p < 0.05; re-run w1 for FDR).",
        "\nBottom row: three-tier significance classification.  ",
        "Negative \u03c4 (blue) = drying trend; positive \u03c4 (red) = wetting trend."),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 10, face = "bold",
                                              hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey35",
                                              hjust = 0,
                                              margin = ggplot2::margin(b = 4))
      )
    )
  
  # ── save ─────────────────────────────────────────────────────────────────────
  tryCatch({
    ggplot2::ggsave(outpdf, fig8,
                    width = 7.0, height = 8.5, units = "in", device = "pdf")
    cat(sprintf("  \u2713 Fig8 (PDF): %s\n", basename(outpdf)))
  }, error = function(e) cat(sprintf("  \u26a0 Fig8 PDF: %s\n", e$message)))
  
  tryCatch({
    ggplot2::ggsave(outpng, fig8,
                    width = 7.0, height = 8.5, units = "in",
                    dpi = dpi, device = "png")
    cat(sprintf("  \u2713 Fig8 (PNG): %s\n", basename(outpng)))
  }, error = function(e) cat(sprintf("  \u26a0 Fig8 PNG: %s\n", e$message)))
  
  invisible(fig8)
}

# ── run manuscript figure 8 ─────────────────────────────────────────────────
tryCatch(
  create_ms_trend_maps(
    scales_ms = c(1, 3),
    index_ms  = "spei",
    outpdf    = file.path(fig_dir, "Fig8_MS_TrendMap_SPEI1_SPEI3.pdf"),
    outpng    = file.path(fig_dir, "Fig8_MS_TrendMap_SPEI1_SPEI3.png"),
    dpi       = 300),
  error = function(e) cat("  \u274c Fig 8 manuscript trend maps:", e$message, "\n"))

################################################################################
# FINAL SUMMARY
################################################################################
cat("\n╔══════════════════════════════════════╗\n")
cat("║  w4_trends_visualization.R  DONE     ║\n")
cat("╚══════════════════════════════════════╝\n\n")
cat("Outputs:\n")
cat("  Master diagnostics CSV : ", normalizePath(master_csv), "\n")
cat("  Spatial figures        : ", normalizePath(fig_dir),    "\n")
cat("  NEW manuscript figure  : Fig8_MS_TrendMap_SPEI1_SPEI3.pdf/.png\n")
cat("\n  Basin-averaged time series, PNG plots, and Excel summary\n")
cat("  are produced by w1_basin_timeseries.R.\n\n")
cat("\u2713 All spatial visualization tasks complete!\n\n")