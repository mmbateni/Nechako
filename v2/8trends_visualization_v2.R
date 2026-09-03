# 8trends_visualization_v2.R  ·  SPATIAL DROUGHT DIAGNOSTICS
# ─────────────────────────────────────────────────────────────────────────────────
# COMPATIBLE WITH: DROUGHT_ANALYSIS_utils_v2.R, 6trend_test_ALL_v2.R
# SUPPORTS: SPI, SPEI-PM, SPEI-Thw, SWEI, MSPI, MSPEI (Gridded)
# FRAMEWORK: observed PM–Thornthwaite comparison only; no detrended/counterfactual branch.
#
# SCOPE (spatial figures only):
#   PART 1 – Merge individual result CSVs → master diagnostics CSV
#   PART 2 – Spatial figures
#     Fig 1: Trend maps          (Kendall τ_vc + FDR significance)
#     Fig 2: Event maps          (# events, mean duration, max intensity)
#     Fig 3: Temporal pattern maps (clustering, regime-shift decade)
#     Fig 4: Timescale comparison (% significant, median τ, events, duration)
#     Fig 5: Method comparison   (VC vs TFPW scatter + agreement bar)
#   PART 3 – Manuscript Figure X
#     3 × 4 Sen-slope maps for SPI, SPEI-PM and SPEI-Thw at 3, 12, 24, 36 months
#
# NOTE: Basin-averaged time series CSVs, PNG plots, and the Excel summary
#       are produced by 7basin_timeseries_v2.R — run that script for all
#       non-spatial outputs.
#
# EXECUTION ORDER: 6trend_test_ALL_v2.R → 7basin_timeseries_v2.R → 8trends_visualization_v2.R

source("DROUGHT_ANALYSIS_utils_v2.R")
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

index_display <- function(x) {
  out <- toupper(x)
  out[x == "spei"] <- "SPEI-PM"
  out[x == "spei_thw"] <- "SPEI-Thw"
  out
}

# ── Index colours (extend for MSPI / MSPEI) ─────────────────────────
if (!"swei"  %in% names(index_colours))
  index_colours <- c(index_colours, swei  = "#1E90FF")
if (!"mspi"  %in% names(index_colours))
  index_colours <- c(index_colours, mspi  = "#33a02c")   # green
if (!"mspei" %in% names(index_colours))
  index_colours <- c(index_colours, mspei = "#fb9a99")   # light red
if (!"spei_thw" %in% names(index_colours))
  index_colours <- c(index_colours, spei_thw = "#d62728") # Thornthwaite observed

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

cat("\n  Loading SPI, SPEI-PM & SPEI-Thw...\n")
for (idx in c("spi", "spei", "spei_thw")) {
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
  stop("No data loaded. Run 6trend_test_ALL_v2.R first.")

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
  # SPEI-PM and SPEI-Thw are parallel observed formulations.
  base <- c("spi", "spei")
  if (any(plot_data$index_type == "spei_thw" & plot_data$timescale == scale))
    base <- c(base, "spei_thw")
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
        "\n    (re-run 6trend_test_ALL_v2.R to generate FDR-corrected columns)\n")
  
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
         paste("Missing columns – re-run 6trend_test_ALL_v2.R:\n",
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
      ggplot2::scale_color_manual(values = index_colours, labels = index_display) +
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
    ggplot2::scale_color_manual(values = index_colours, labels = index_display) +
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
           "TFPW columns not found.\nRe-run 6trend_test_ALL_v2.R with TFPW enabled.",
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
# PART 3 – MANUSCRIPT FIGURE X
# Publication-quality 3 × 4 spatial drought-trend figure.
#
# Rows    : SPI, SPEI-PM, SPEI-Thw
# Columns : 3-, 12-, 24-, and 36-month accumulation scales
# Field   : TFPW Sen's slope (sl_tfpw)
# Overlay : BH-FDR-significant TFPW tests (p_fdr_tfpw < 0.05)
#
# This figure is deliberately separate from the legacy diagnostic figures above.
# It does not alter TIMESCALES_VIZ or the statistical products produced by
# 6trend_test_ALL_v2.R.
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 3: Manuscript Figure X — spatial drought trends\n")
cat("══════════════════════════════════════\n")

create_ms_drought_senslope_maps <- function(
    scales_ms = c(3, 12, 24, 36),
    outpdf,
    outpng,
    dpi = 300) {
  
  required <- c("index_type", "timescale", "lon", "lat",
                "sl_tfpw", "p_value_tfpw", "p_fdr_tfpw")
  missing <- setdiff(required, names(plot_data))
  if (length(missing)) {
    stop("Required TFPW/Sen/FDR columns are missing: ",
         paste(missing, collapse = ", "))
  }
  
  idx_order <- c("spi", "spei", "spei_thw")
  idx_labels <- c(
    spi      = "SPI",
    spei     = "SPEI-PM",
    spei_thw = "SPEI-Thw"
  )
  
  dat <- data.table::copy(plot_data[
    index_type %in% idx_order & timescale %in% scales_ms
  ])
  
  if (!nrow(dat)) stop("No SPI/SPEI data found for requested manuscript scales.")
  
  # IMPORTANT: the trend-result files written by 6trend_test_ALL_v2.R store
  # pixel centres as x/y coordinates on the canonical BC Albers analysis grid
  # (EPSG:3005). load_clip_results() renames x/y to lon/lat for compatibility
  # with the older plotting helpers, but those values remain projected metre
  # coordinates. Therefore this manuscript figure is drawn entirely in the
  # same EPSG:3005 CRS as the basin boundary. No longitude/latitude conversion
  # is performed here, eliminating the CRS mismatch that produced the previous
  # misplaced/tiny maps and the later vertical strips of black markers.
  map_crs <- sf::st_crs(EQUAL_AREA_CRS)
  basin_map <- sf::st_transform(basin_sf, map_crs)
  basin_bbox <- sf::st_bbox(basin_map)
  
  xy_ok <- is.finite(dat$lon) & is.finite(dat$lat)
  dat <- dat[xy_ok]
  if (!nrow(dat)) stop("No finite pixel coordinates available for manuscript figure.")
  
  # Sanity check: trend-pixel coordinates should be on the same metre-scale
  # order of magnitude as the BC Albers basin extent.
  if (max(abs(c(dat$lon, dat$lat)), na.rm = TRUE) < 1e4) {
    stop("Trend coordinates do not appear to be EPSG:3005 metre coordinates; check the input result CSVs.")
  }
  
  # Keep only the basin-supported pixel records already selected by the
  # load_clip_results() basin clip. The spatial inference values themselves
  # are not recalculated or altered here.
  cat("  ✓ Using trend pixels and basin boundary directly in EPSG:3005\n")
  
  # One symmetric slope range per index. This keeps accumulation scales
  # comparable within each index while avoiding a large Thornthwaite slope
  # range suppressing the visual variation in SPI/SPEI-PM.
  z_tbl <- dat[is.finite(sl_tfpw),
               .(z = max(abs(sl_tfpw), na.rm = TRUE)), by = index_type]
  z_lookup <- setNames(z_tbl$z, z_tbl$index_type)
  z_lookup[!is.finite(z_lookup) | z_lookup <= 0] <- 1
  
  slope_cols <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  
  make_panel <- function(idx, sc, show_legend = FALSE) {
    sub <- dat[index_type == idx & timescale == sc]
    
    if (!nrow(sub)) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate(
            "text", x = 0.5, y = 0.5,
            label = sprintf("No data\n%s — %d months", idx_labels[[idx]], sc),
            size = 4
          ) +
          ggplot2::theme_void()
      )
    }
    
    sub_plot <- sub[
      is.finite(sl_tfpw) & is.finite(lon) & is.finite(lat)
    ]
    
    sig <- sub[
      is.finite(p_fdr_tfpw) & p_fdr_tfpw < 0.05 &
        is.finite(sl_tfpw) & is.finite(lon) & is.finite(lat)
    ]
    
    n_total <- sum(is.finite(sub$sl_tfpw))
    n_sig <- nrow(sig)
    pct_sig <- if (n_total) 100 * n_sig / n_total else NA_real_
    
    z <- z_lookup[[idx]]
    slope_breaks <- seq(-z, z, length.out = 5)
    
    # The pixels are regular in EPSG:3005. Use the median centre spacing to
    # render the actual pixel footprint rather than using point markers.
    dx <- diff(sort(unique(sub_plot$lon)))
    dy <- diff(sort(unique(sub_plot$lat)))
    dx <- dx[is.finite(dx) & dx > 0]
    dy <- dy[is.finite(dy) & dy > 0]
    tile_w <- if (length(dx)) stats::median(dx) else 5000
    tile_h <- if (length(dy)) stats::median(dy) else 5000
    
    guide_spec <- if (show_legend) {
      ggplot2::guide_colorbar(
        barheight = grid::unit(2.2, "cm"),
        barwidth  = grid::unit(0.28, "cm"),
        ticks.colour = "grey20",
        frame.colour = "grey40")
    } else {
      "none"
    }
    
    ggplot2::ggplot(
      sub_plot,
      ggplot2::aes(x = lon, y = lat, fill = sl_tfpw)
    ) +
      ggplot2::geom_tile(
        width = tile_w,
        height = tile_h,
        na.rm = TRUE) +
      ggplot2::scale_fill_gradientn(
        colours = slope_cols,
        limits = c(-z, z),
        breaks = slope_breaks,
        labels = scales::label_scientific(digits = 1),
        oob = scales::squish,
        name = "Sen's slope",
        guide = guide_spec) +
      # Stippling is deliberately overlaid only after the coloured pixels.
      # The significance layer is therefore visible without replacing the
      # slope field with black points.
      ggplot2::geom_point(
        data = as.data.frame(sig),
        ggplot2::aes(x = lon, y = lat),
        inherit.aes = FALSE,
        shape = 46,
        size = 0.7,
        colour = "black",
        alpha = 0.55) +
      ggplot2::geom_sf(
        data = basin_map,
        inherit.aes = FALSE,
        fill = NA,
        colour = "black",
        linewidth = 0.45) +
      ggplot2::coord_sf(
        crs = map_crs,
        xlim = c(basin_bbox["xmin"], basin_bbox["xmax"]),
        ylim = c(basin_bbox["ymin"], basin_bbox["ymax"]),
        expand = FALSE,
        datum = NA) +
      ggplot2::labs(
        title = sprintf("%s — %d months", idx_labels[[idx]], sc),
        subtitle = sprintf("FDR-significant: %d/%d pixels (%.1f%%)",
                           n_sig, n_total, pct_sig)) +
      ggplot2::theme_void(base_size = 8.5) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 9, face = "bold", hjust = 0.5,
          margin = ggplot2::margin(b = 2)),
        plot.subtitle = ggplot2::element_text(
          size = 6.5, hjust = 0.5, colour = "grey30",
          margin = ggplot2::margin(b = 2)),
        legend.title = ggplot2::element_text(size = 7, face = "bold"),
        legend.text = ggplot2::element_text(size = 6),
        legend.position = if (show_legend) "right" else "none",
        plot.margin = ggplot2::margin(2, 3, 2, 3))
  }
  
  if (!requireNamespace("magick", quietly = TRUE))
    stop("Package 'magick' is required to assemble Fig X. Install with install.packages('magick').")
  
  # Each panel is rendered independently on a generous square canvas and
  # trimmed before four panels are joined. This removes coord_sf aspect-ratio
  # padding inside each panel rather than leaving blank gutters between maps.
  panel_canvas_in <- 3.2
  legend_w_in <- 1.6
  legend_h_in <- 2.8
  
  extract_legend_grob <- function(p) {
    gt <- ggplot2::ggplotGrob(p)
    leg_names <- gt$layout$name[grepl("^guide-box", gt$layout$name)]
    leg_names <- setdiff(leg_names, "guide-box-nowhere")
    if (!length(leg_names))
      stop("No legend grob found in panel (check show_legend = TRUE).")
    leg_pos <- which(gt$layout$name == leg_names[1])[1]
    gt$grobs[[leg_pos]]
  }
  
  save_grob_png <- function(grob, file, width_in, height_in, dpi) {
    grDevices::png(file, width = width_in, height = height_in,
                   units = "in", res = dpi, bg = "white")
    grid::grid.newpage()
    grid::grid.draw(grob)
    grDevices::dev.off()
    invisible(file)
  }
  
  tmp_dir <- tempfile("figx_panels_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  
  row_pngs <- vapply(idx_order, function(idx) {
    panel_imgs <- lapply(scales_ms, function(sc) {
      p <- make_panel(idx, sc, show_legend = FALSE)
      f <- file.path(tmp_dir, sprintf("panel_%s_%02d.png", idx, sc))
      ggplot2::ggsave(
        f, p,
        width = panel_canvas_in,
        height = panel_canvas_in,
        units = "in",
        dpi = dpi,
        bg = "white")
      magick::image_trim(magick::image_read(f), fuzz = 2)
    })
    
    common_pw <- max(vapply(
      panel_imgs,
      function(im) magick::image_info(im)$width,
      integer(1)))
    common_ph <- max(vapply(
      panel_imgs,
      function(im) magick::image_info(im)$height,
      integer(1)))
    panel_imgs <- lapply(panel_imgs, function(im)
      magick::image_extent(
        im,
        geometry = sprintf("%dx%d", common_pw, common_ph),
        gravity = "center",
        color = "white"))
    maps_strip <- magick::image_append(
      magick::image_join(panel_imgs),
      stack = FALSE)
    
    leg_plot <- make_panel(idx, scales_ms[1], show_legend = TRUE)
    leg_grob <- extract_legend_grob(leg_plot)
    lf <- file.path(tmp_dir, sprintf("legend_%s.png", idx))
    save_grob_png(leg_grob, lf, legend_w_in, legend_h_in, dpi)
    leg_img <- magick::image_trim(magick::image_read(lf), fuzz = 2)
    
    strip_h <- magick::image_info(maps_strip)$height
    leg_img <- magick::image_extent(
      leg_img,
      geometry = sprintf(
        "%dx%d",
        magick::image_info(leg_img)$width,
        strip_h),
      gravity = "center",
      color = "white")
    
    row_img <- magick::image_append(
      magick::image_join(list(maps_strip, leg_img)),
      stack = FALSE)
    f <- file.path(tmp_dir, sprintf("row_%s.png", idx))
    magick::image_write(row_img, path = f, format = "png", density = dpi)
    f
  }, character(1))
  
  row_imgs <- lapply(row_pngs, magick::image_read)
  common_w <- max(vapply(
    row_imgs,
    function(im) magick::image_info(im)$width,
    integer(1)))
  row_imgs <- lapply(row_imgs, function(im)
    magick::image_extent(
      im,
      geometry = sprintf("%dx", common_w),
      gravity = "center",
      color = "white"))
  body_img <- magick::image_append(
    magick::image_join(row_imgs),
    stack = TRUE)
  
  title_txt <- "Spatial drought-index trends — Nechako River Basin (1950–2025)"
  subtitle_txt <- paste0(
    "TFPW Mann–Kendall Sen's slope; stippling = BH-FDR-significant pixels (α = 0.05). ",
    "Columns: 3-, 12-, 24-, and 36-month accumulation scales.")
  header_h <- round(dpi * 0.85)
  header <- magick::image_blank(
    width = common_w,
    height = header_h,
    color = "white")
  header <- magick::image_annotate(
    header,
    title_txt,
    size = round(dpi * 0.15),
    weight = 700,
    gravity = "north",
    location = sprintf("+0+%d", round(dpi * 0.08)))
  header <- magick::image_annotate(
    header,
    subtitle_txt,
    size = round(dpi * 0.085),
    color = "grey30",
    gravity = "north",
    location = sprintf("+0+%d", round(dpi * 0.42)))
  header <- magick::image_trim(header, fuzz = 2)
  header <- magick::image_border(
    header,
    "white",
    sprintf("0x%d", round(dpi * 0.08)))
  header <- magick::image_extent(
    header,
    geometry = sprintf("%dx", common_w),
    gravity = "center",
    color = "white")
  
  final_img <- magick::image_append(
    magick::image_join(list(header, body_img)),
    stack = TRUE)
  
  magick::image_write(final_img, path = outpng, format = "png", density = dpi)
  magick::image_write(final_img, path = outpdf, format = "pdf", density = dpi)
  
  cat(sprintf("  ✓ Manuscript drought-trend figure: %s\n", basename(outpdf)))
  cat(sprintf("  ✓ Manuscript drought-trend figure: %s\n", basename(outpng)))
  invisible(final_img)
}

tryCatch(
  create_ms_drought_senslope_maps(
    scales_ms = c(3, 12, 24, 36),
    outpdf = file.path(
      fig_dir,
      "FigX_MS_Spatial_Drought_Trends_SenSlope_3_12_24_36.pdf"
    ),
    outpng = file.path(
      fig_dir,
      "FigX_MS_Spatial_Drought_Trends_SenSlope_3_12_24_36.png"
    ),
    dpi = DPI
  ),
  error = function(e) {
    cat("  ❌ Manuscript drought-trend figure:", e$message, "\n")
  }
)

################################################################################
# LEGACY / EXISTING FIGURE 8 OUTPUTS
# Kept unchanged so existing diagnostic products are preserved.
################################################################################
cat("\n  Existing diagnostic trend-map outputs are retained above; the new\n")
cat("  manuscript figure is the dedicated 3 × 4 Sen-slope + BH-FDR product.\n")

################################################################################
# FINAL SUMMARY
################################################################################
cat("\n╔══════════════════════════════════════╗\n")
cat("║  8trends_visualization_v2.R   DONE     ║\n")
cat("╚══════════════════════════════════════╝\n\n")
cat("Outputs:\n")
cat("  Master diagnostics CSV : ", normalizePath(master_csv), "\n")
cat("  Spatial figures        : ", normalizePath(fig_dir),    "\n")
cat("  Manuscript Figure X    : FigX_MS_Spatial_Drought_Trends_SenSlope_3_12_24_36.pdf/.png\n")
cat("\n  Basin-averaged time series, PNG plots, and Excel summary\n")
cat("  are produced by 7basin_timeseries_v2.R.\n\n")
cat("\u2713 All spatial visualization tasks complete!\n\n")
################################################################################
# PART 4 – PIXEL-TO-BASIN SPEARMAN CORRELATION DIAGNOSTICS
#
# GRID POLICY — IMPORTANT:
#   ALL calculations in this section use the exact canonical grid established
#   upstream by 3SPEI_ERALand_v2.R.  That grid is the projected ERA5-Land
#   precipitation/SPI analysis grid on which SPI, SPEI-PM and SPEI-Thw are
#   required to remain one-to-one aligned.
#
#   NO reproject(), resample(), project(), crop(), extend(), or rasterize()
#   operation is performed on the drought-index rasters in this section.
#   The canonical template and basin mask are read from:
#
#       spei_results_seasonal/spei_analysis_grid.rds
#
#   Every seasonal SPI/SPEI NetCDF is hard-checked against that exact grid
#   (rows, columns, resolution, extent, origin, CRS, and cell count). A
#   mismatch is a fatal error rather than silently changing the grid.
#
# Purpose:
#   For every canonical basin cell i:
#
#       rho_i = Spearman( X_i(t), X_basin(t) )
#
#   using the original monthly gridded SPI/SPEI series and the corresponding
#   authoritative basin-average series from 6trend_test_ALL_v2.R.
################################################################################

cat("\n============================================================\n")
cat("PART 4: PIXEL-TO-BASIN SPEARMAN CORRELATION DIAGNOSTICS\n")
cat("============================================================\n")

CORR_OUT_DIR <- file.path(fig_dir, "pixel_basin_correlation")
dir.create(CORR_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CORR_SHORT_SCALES <- c(1L, 3L, 6L)
CORR_LONG_SCALES  <- c(12L, 24L, 36L)
CORR_ALL_SCALES   <- c(CORR_SHORT_SCALES, CORR_LONG_SCALES)
CORR_INDEXES      <- c("spi", "spei")
CORR_MIN_N        <- 24L

# -----------------------------------------------------------------------------
# Canonical grid: use the exact common SPI/SPEI grid created by Script 3.
# -----------------------------------------------------------------------------
CORR_GRID_FILE <- file.path(SPEI_SEAS_DIR, "spei_analysis_grid.rds")
if (!file.exists(CORR_GRID_FILE)) {
  stop(
    "Canonical SPI/SPEI grid not found: ", CORR_GRID_FILE,
    "\nRun 3SPEI_ERALand_v2.R first."
  )
}

corr_grid <- readRDS(CORR_GRID_FILE)

unwrap_spatraster <- function(x, label) {
  if (inherits(x, "SpatRaster")) return(x)
  if (inherits(x, "PackedSpatRaster")) {
    y <- tryCatch(terra::unwrap(x), error = function(e) NULL)
    if (inherits(y, "SpatRaster")) return(y)
  }
  stop("Cannot recover canonical ", label, " as a terra::SpatRaster.")
}

CORR_TEMPLATE <- unwrap_spatraster(
  corr_grid$raster_template,
  "raster template"
)
CORR_MASK_RASTER <- unwrap_spatraster(
  corr_grid$basin_mask,
  "basin mask"
)

# The canonical mask is itself on the canonical grid. Never regenerate it from
# the polygon here: doing so could alter cell inclusion at the basin edge.
if (terra::ncell(CORR_TEMPLATE) != terra::ncell(CORR_MASK_RASTER) ||
    terra::nrow(CORR_TEMPLATE) != terra::nrow(CORR_MASK_RASTER) ||
    terra::ncol(CORR_TEMPLATE) != terra::ncol(CORR_MASK_RASTER)) {
  stop("Canonical raster template and basin mask have different dimensions.")
}

if (!isTRUE(all.equal(
  as.numeric(terra::res(CORR_TEMPLATE)),
  as.numeric(terra::res(CORR_MASK_RASTER)),
  tolerance = 1e-10
))) {
  stop("Canonical raster template and basin mask have different resolutions.")
}

corr_ext <- c(
  terra::xmin(CORR_TEMPLATE), terra::xmax(CORR_TEMPLATE),
  terra::ymin(CORR_TEMPLATE), terra::ymax(CORR_TEMPLATE)
)
mask_ext <- c(
  terra::xmin(CORR_MASK_RASTER), terra::xmax(CORR_MASK_RASTER),
  terra::ymin(CORR_MASK_RASTER), terra::ymax(CORR_MASK_RASTER)
)
if (!isTRUE(all.equal(corr_ext, mask_ext, tolerance = 1e-8))) {
  stop("Canonical raster template and basin mask have different extents.")
}

corr_origin <- terra::origin(CORR_TEMPLATE)
mask_origin <- terra::origin(CORR_MASK_RASTER)
if (!isTRUE(all.equal(
  as.numeric(corr_origin),
  as.numeric(mask_origin),
  tolerance = 1e-8
))) {
  stop("Canonical raster template and basin mask have different origins.")
}

if (!terra::same.crs(CORR_TEMPLATE, CORR_MASK_RASTER)) {
  stop("Canonical raster template and basin mask have different CRS.")
}

CORR_BASIN_CELLS <- as.logical(terra::values(CORR_MASK_RASTER, mat = FALSE))
if (length(CORR_BASIN_CELLS) != terra::ncell(CORR_TEMPLATE)) {
  stop("Canonical basin mask length does not equal canonical grid cell count.")
}

if (!is.null(corr_grid$basin_pixels) &&
    sum(CORR_BASIN_CELLS, na.rm = TRUE) != as.integer(corr_grid$basin_pixels)) {
  stop("Canonical basin-mask cell count does not match spei_analysis_grid.rds metadata.")
}

cat(sprintf(
  "Canonical SPI/SPEI grid: %d rows × %d cols = %d cells\n",
  terra::nrow(CORR_TEMPLATE),
  terra::ncol(CORR_TEMPLATE),
  terra::ncell(CORR_TEMPLATE)
))
cat(sprintf(
  "Canonical resolution: %.10f × %.10f\n",
  terra::res(CORR_TEMPLATE)[1], terra::res(CORR_TEMPLATE)[2]
))
cat(sprintf(
  "Canonical CRS: %s\n",
  terra::crs(CORR_TEMPLATE)
))
cat(sprintf(
  "Canonical basin cells: %d\n",
  sum(CORR_BASIN_CELLS, na.rm = TRUE)
))

# -----------------------------------------------------------------------------
# Exact geometry validator.
# No correction is permitted: a mismatch stops the calculation.
# -----------------------------------------------------------------------------
validate_corr_geometry <- function(r, reference, file_label) {
  if (terra::ncell(r) != terra::ncell(reference)) {
    stop(sprintf(
      "GRID ERROR [%s]: cell count %d != canonical %d",
      file_label, terra::ncell(r), terra::ncell(reference)
    ))
  }
  if (terra::nrow(r) != terra::nrow(reference) ||
      terra::ncol(r) != terra::ncol(reference)) {
    stop(sprintf(
      "GRID ERROR [%s]: dimensions %d×%d != canonical %d×%d",
      file_label,
      terra::nrow(r), terra::ncol(r),
      terra::nrow(reference), terra::ncol(reference)
    ))
  }
  if (!isTRUE(all.equal(
    as.numeric(terra::res(r)),
    as.numeric(terra::res(reference)),
    tolerance = 1e-10
  ))) {
    stop("GRID ERROR [", file_label, "]: resolution differs from canonical grid.")
  }
  ext_r <- c(
    terra::xmin(r), terra::xmax(r),
    terra::ymin(r), terra::ymax(r)
  )
  ext_ref <- c(
    terra::xmin(reference), terra::xmax(reference),
    terra::ymin(reference), terra::ymax(reference)
  )
  if (!isTRUE(all.equal(ext_r, ext_ref, tolerance = 1e-8))) {
    stop("GRID ERROR [", file_label, "]: extent differs from canonical grid.")
  }
  if (!isTRUE(all.equal(
    as.numeric(terra::origin(r)),
    as.numeric(terra::origin(reference)),
    tolerance = 1e-8
  ))) {
    stop("GRID ERROR [", file_label, "]: origin differs from canonical grid.")
  }
  if (!terra::same.crs(r, reference)) {
    stop("GRID ERROR [", file_label, "]: CRS differs from canonical grid.")
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Locate seasonal NetCDF files.
# -----------------------------------------------------------------------------
find_corr_nc_files <- function(index_type, scale) {
  data_dir <- switch(
    tolower(index_type),
    spi  = SPI_SEAS_DIR,
    spei = SPEI_SEAS_DIR,
    stop("Unsupported correlation index: ", index_type)
  )
  files <- find_seasonal_nc_files(data_dir, index_type, scale)
  if (!length(files)) {
    stop(sprintf(
      "No seasonal NetCDF files found for %s-%d in %s",
      toupper(index_type), scale, data_dir
    ))
  }
  files
}

# -----------------------------------------------------------------------------
# Read a complete seasonal product into chronological monthly order.
# The canonical template is NOT inferred from the NetCDF and the NetCDF is
# NOT reprojected/resampled to it. Every file must already be exactly aligned.
# -----------------------------------------------------------------------------
read_corr_index_matrix <- function(index_type, scale, target_dates) {
  files <- find_corr_nc_files(index_type, scale)
  if (length(files) != 12L) {
    stop(sprintf(
      "Expected 12 seasonal files for %s-%d; found %d",
      toupper(index_type), scale, length(files)
    ))
  }
  
  out <- matrix(
    NA_real_,
    nrow = terra::ncell(CORR_TEMPLATE),
    ncol = length(target_dates)
  )
  nlayer_total <- 0L
  
  for (m in seq_along(files)) {
    r <- terra::rast(files[m])
    file_label <- basename(files[m])
    validate_corr_geometry(r, CORR_TEMPLATE, file_label)
    
    nlayer_total <- nlayer_total + terra::nlyr(r)
    vals <- terra::values(r, mat = TRUE)
    vals <- as.matrix(vals)
    vals[!is.finite(vals)] <- NA_real_
    
    month_dates <- target_dates[
      format(target_dates, "%m") == sprintf("%02d", m)
    ]
    if (length(month_dates) != ncol(vals)) {
      stop(sprintf(
        "%s-%d month %02d has %d raster layers but %d target dates",
        toupper(index_type), scale, m, ncol(vals), length(month_dates)
      ))
    }
    
    col_idx <- match(month_dates, target_dates)
    out[, col_idx] <- vals
    
    rm(r, vals)
    gc(verbose = FALSE)
  }
  
  if (nlayer_total != length(target_dates)) {
    stop(sprintf(
      "Unexpected layer count for %s-%d: %d layers; expected %d",
      toupper(index_type), scale, nlayer_total, length(target_dates)
    ))
  }
  
  list(matrix = out, files = files)
}

# -----------------------------------------------------------------------------
# Pixel-wise Spearman correlation.
# The output raster is created directly from the canonical grid template.
# -----------------------------------------------------------------------------
compute_corr_raster <- function(ts_matrix, basin_series, index_type, scale) {
  if (!is.matrix(ts_matrix)) stop("ts_matrix must be a matrix")
  if (nrow(ts_matrix) != terra::ncell(CORR_TEMPLATE)) {
    stop(sprintf(
      "Cell-count mismatch for %s-%d: matrix has %d rows, canonical grid has %d cells",
      index_type, scale, nrow(ts_matrix), terra::ncell(CORR_TEMPLATE)
    ))
  }
  if (ncol(ts_matrix) != length(basin_series)) {
    stop(sprintf(
      "Temporal mismatch for %s-%d: %d pixel months vs %d basin months",
      index_type, scale, ncol(ts_matrix), length(basin_series)
    ))
  }
  
  n_pair <- integer(nrow(ts_matrix))
  rr <- rep(NA_real_, nrow(ts_matrix))
  
  for (i in which(CORR_BASIN_CELLS)) {
    x <- ts_matrix[i, ]
    ok <- is.finite(x) & is.finite(basin_series)
    n_pair[i] <- sum(ok)
    
    if (n_pair[i] < CORR_MIN_N) next
    
    xok <- x[ok]
    yok <- basin_series[ok]
    if (length(unique(xok)) < 2L || length(unique(yok)) < 2L) next
    
    rr[i] <- suppressWarnings(
      stats::cor(xok, yok, method = "spearman")
    )
  }
  
  r <- terra::setValues(
    terra::rast(CORR_TEMPLATE),
    rr
  )
  
  list(raster = r, n_pair = n_pair)
}

# -----------------------------------------------------------------------------
# Map data and map plot.
# -----------------------------------------------------------------------------
corr_map_df <- function(r, label) {
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  if (!nrow(d)) stop("No finite correlation values available for ", label)
  names(d)[3] <- "rho"
  d$label <- label
  d
}

plot_corr_map <- function(r, label) {
  d <- corr_map_df(r, label)
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, fill = rho)) +
    ggplot2::geom_raster() +
    ggplot2::geom_sf(
      data = basin_sf,
      inherit.aes = FALSE,
      fill = NA,
      linewidth = 0.35
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      name = expression(rho[S])
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(title = label) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank()
    )
}

# -----------------------------------------------------------------------------
# Density plot: three curves + one dotted median line per timescale.
# -----------------------------------------------------------------------------
plot_corr_density <- function(corr_values, index_type, scales) {
  rows <- lapply(scales, function(sc) {
    key <- sprintf("%s_%02d", index_type, sc)
    x <- corr_values[[key]]
    if (is.null(x)) return(NULL)
    data.frame(
      rho = x,
      scale = factor(sc, levels = scales),
      stringsAsFactors = FALSE
    )
  })
  
  d <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  d <- as.data.frame(d)
  d <- d[is.finite(d$rho), , drop = FALSE]
  if (!nrow(d)) {
    stop("No finite correlation values available for density plot: ", index_type)
  }
  
  meds <- aggregate(rho ~ scale, data = d, FUN = median, na.rm = TRUE)
  
  ggplot2::ggplot(d, ggplot2::aes(x = rho, colour = scale)) +
    ggplot2::geom_density(
      linewidth = 1.15,
      na.rm = TRUE,
      adjust = 1.05
    ) +
    ggplot2::geom_vline(
      data = meds,
      ggplot2::aes(xintercept = rho, colour = scale),
      linetype = "dotted",
      linewidth = 0.85,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c("#0072B2", "#E69F00", "#009E73")[seq_along(scales)],
      name = "Timescale"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.2)
    ) +
    ggplot2::labs(
      title = if (index_type == "spi") {
        "SPI: pixel–basin Spearman correlation distribution"
      } else {
        "SPEI-PM: pixel–basin Spearman correlation distribution"
      },
      x = expression(Spearman~rho[S]),
      y = "Density"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# -----------------------------------------------------------------------------
# Calculate all requested correlations once.
# -----------------------------------------------------------------------------
corr_rasters <- list()
corr_values  <- list()
corr_meta    <- list()

for (idx in CORR_INDEXES) {
  for (sc in CORR_ALL_SCALES) {
    key <- sprintf("%s_%02d", idx, sc)
    cat(sprintf("\n  → %s-%d\n", toupper(idx), sc))
    
    basin_df <- load_basin_avg_csv(idx, sc)
    if (is.null(basin_df)) {
      stop(sprintf(
        "Required basin-average CSV is missing for %s-%d. Run 6trend_test_ALL_v2.R first.",
        toupper(idx), sc
      ))
    }
    
    basin_dates  <- as.Date(basin_df$date)
    basin_series <- as.numeric(basin_df$value)
    
    nc <- read_corr_index_matrix(idx, sc, basin_dates)
    
    # Enforce the canonical basin mask exactly; no spatial operation is applied.
    ts_basin <- nc$matrix
    ts_basin[!CORR_BASIN_CELLS, ] <- NA_real_
    
    corr_result <- compute_corr_raster(
      ts_matrix = ts_basin,
      basin_series = basin_series,
      index_type = idx,
      scale = sc
    )
    
    r_corr <- corr_result$raster
    vals <- terra::values(r_corr, mat = FALSE)
    vals <- vals[is.finite(vals)]
    
    if (!length(vals)) {
      stop(sprintf(
        "No finite pixel-basin correlations for %s-%d",
        toupper(idx), sc
      ))
    }
    
    corr_rasters[[key]] <- r_corr
    corr_values[[key]] <- vals
    corr_meta[[key]] <- data.frame(
      index_type = idx,
      index_label = if (idx == "spi") "SPI" else "SPEI-PM",
      scale = sc,
      n_pixels = length(vals),
      median_rho = median(vals, na.rm = TRUE),
      mean_rho = mean(vals, na.rm = TRUE),
      min_rho = min(vals, na.rm = TRUE),
      max_rho = max(vals, na.rm = TRUE),
      n_paired_months_min = min(corr_result$n_pair[CORR_BASIN_CELLS], na.rm = TRUE),
      n_paired_months_max = max(corr_result$n_pair[CORR_BASIN_CELLS], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    cat(sprintf(
      "     %d basin pixels; median Spearman rho = %.3f; range = %.3f to %.3f\n",
      length(vals), median(vals), min(vals), max(vals)
    ))
    
    rm(nc, ts_basin, r_corr, corr_result, vals, basin_df)
    gc(verbose = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Audit outputs.
# -----------------------------------------------------------------------------
corr_summary <- data.table::rbindlist(
  corr_meta,
  use.names = TRUE,
  fill = TRUE
)
corr_summary <- as.data.frame(corr_summary)

data.table::fwrite(
  corr_summary,
  file.path(CORR_OUT_DIR, "PixelBasin_Spearman_Summary.csv")
)

# Save one row per pixel/index/scale. Coordinates are taken directly from the
# canonical grid; no coordinate transformation is performed here.
# IMPORTANT:
# terra::crds(SpatRaster) may return coordinates only for non-NA cells
# depending on the raster state/version.  That would break the one-row-per-
# canonical-cell audit table because the canonical grid has 2021 cells while
# only 865 cells are inside the basin.  Recover coordinates explicitly for
# EVERY canonical raster cell, in native cell-number order.
cell_id <- seq_len(terra::ncell(CORR_TEMPLATE))
xy <- terra::xyFromCell(CORR_TEMPLATE, cell_id)

if (nrow(xy) != terra::ncell(CORR_TEMPLATE)) {
  stop(
    "Canonical-grid coordinate extraction failed: returned ",
    nrow(xy), " coordinates for ",
    terra::ncell(CORR_TEMPLATE), " cells."
  )
}

correlation_audit <- data.table::rbindlist(
  lapply(names(corr_rasters), function(key) {
    r <- corr_rasters[[key]]
    parts <- strsplit(key, "_", fixed = TRUE)[[1]]
    idx <- parts[1]
    sc <- as.integer(parts[2])
    
    vals <- terra::values(r, mat = FALSE)
    
    if (length(vals) != length(cell_id)) {
      stop(
        "Correlation raster length mismatch for ", key,
        ": ", length(vals), " values for ",
        length(cell_id), " canonical cells."
      )
    }
    
    data.table::data.table(
      cell = cell_id,
      x = xy[, 1],
      y = xy[, 2],
      index_type = idx,
      index_label = if (idx == "spi") "SPI" else "SPEI-PM",
      scale = sc,
      spearman_rho = vals
    )[
      CORR_BASIN_CELLS & is.finite(spearman_rho)
    ]
  }),
  use.names = TRUE,
  fill = TRUE
)



# Final integrity checks: every retained audit point must be a native cell
# centre of the canonical ERA5-Land-derived SPI/SPEI grid and must be one of
# the canonical basin cells.
if (nrow(correlation_audit) == 0L) {
  stop("Correlation audit table is empty after canonical basin-cell filtering.")
}
if (!all(correlation_audit$cell %in% cell_id)) {
  stop("Correlation audit contains cell IDs outside the canonical grid.")
}
if (!all(correlation_audit$cell %in% cell_id[CORR_BASIN_CELLS])) {
  stop("Correlation audit contains cells outside the canonical basin mask.")
}

data.table::fwrite(
  correlation_audit,
  file.path(CORR_OUT_DIR, "PixelBasin_Spearman_Correlations.csv")
)

# -----------------------------------------------------------------------------
# Assemble the two requested figures.
# Top: six maps. Bottom: SPI distribution + SPEI-PM distribution.
# -----------------------------------------------------------------------------
make_corr_figure <- function(scales, outfile_stub) {
  map_keys <- c(
    sprintf("spi_%02d", scales),
    sprintf("spei_%02d", scales)
  )
  
  map_labels <- c(
    sprintf("SPI-%d", scales),
    sprintf("SPEI-PM-%d", scales)
  )
  
  map_plots <- Map(
    f = function(key, label) plot_corr_map(corr_rasters[[key]], label),
    key = map_keys,
    label = map_labels
  )
  
  if (length(map_plots) != 6L) {
    stop(
      "Correlation figure requires exactly six map panels (",
      "3 SPI + 3 SPEI-PM), but got ", length(map_plots), "."
    )
  }
  
  map_grid <- patchwork::wrap_plots(
    map_plots,
    nrow = 2,
    ncol = 3,
    byrow = TRUE
  )
  
  spi_density <- plot_corr_density(corr_values, "spi", scales)
  spei_density <- plot_corr_density(corr_values, "spei", scales)
  
  final_plot <- map_grid /
    spi_density /
    spei_density +
    patchwork::plot_layout(
      heights = c(2.0, 0.72, 0.72)
    )
  
  png_file <- file.path(CORR_OUT_DIR, paste0(outfile_stub, ".png"))
  pdf_file <- file.path(CORR_OUT_DIR, paste0(outfile_stub, ".pdf"))
  
  ggplot2::ggsave(
    png_file,
    final_plot,
    width = 15,
    height = 15,
    units = "in",
    dpi = DPI,
    bg = "white"
  )
  ggplot2::ggsave(
    pdf_file,
    final_plot,
    width = 15,
    height = 15,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  
  cat("  ✓ ", basename(png_file), "\n", sep = "")
  cat("  ✓ ", basename(pdf_file), "\n", sep = "")
}

cat("\n  Writing short-scale correlation figure (1/3/6)...\n")
make_corr_figure(
  CORR_SHORT_SCALES,
  "FigS_PixelBasin_Spearman_1_3_6"
)

cat("\n  Writing long-scale correlation figure (12/24/36)...\n")
make_corr_figure(
  CORR_LONG_SCALES,
  "FigS_PixelBasin_Spearman_12_24_36"
)

cat("\n============================================================\n")
cat("PART 4 COMPLETE: pixel-to-basin Spearman diagnostics\n")
cat("============================================================\n")
cat("Output directory: ", normalizePath(CORR_OUT_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Canonical grid was used without reprojection/resampling.\n")
cat("\n")