# ####################################################################################
# w4_trends_visualization.R  ·  SPATIAL & TEMPORAL DROUGHT DIAGNOSTICS
# ─────────────────────────────────────────────────────────────────────────────────
# COMPATIBLE WITH: DROUGHT_ANALYSIS_utils.R, w3_trend_test.R
# SUPPORTS: SPI (scales in SPI_SCALES), SPEI (scales in SPEI_SCALES), SWEI (SWEI_SCALE)
# PART 1 – Merge individual result files → master diagnostics CSV
# PART 2 – Spatial figures
# Fig 1: Trend maps          (Kendall τ_vc + significance)
# Fig 2: Event maps          (# events, mean duration, max intensity)
#   Fig 3: Temporal pattern maps (clustering, regime-shift decade,
#                                 spectral peaks)  ← NEW
#   Fig 4: Timescale comparison (% significant, median τ, events, duration)
#   Fig 5: Method comparison   (VC vs TFPW scatter + agreement bar)  ← NEW
#   PART 3 – Basin-averaged time series CSVs (one per index × timescale)
#   PART 4 – Time series line plots (PNG per index × timescale)
#   PART 5 – Excel summary statistics + pairwise correlation sheet
#   ####################################################################################
  source("DROUGHT_ANALYSIS_utils.R")
  utils_load_packages(c("terra", "data.table", "ggplot2", "sf", "zoo", "viridis",
                        "cowplot", "patchwork", "scales", "RColorBrewer",
                        "ggspatial", "gridExtra", "grid", "lubridate", "openxlsx"))
  # ── Working directory & output folders ──────────────────────────────
  if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
  setwd(WD_PATH)
  fig_dir         <- file.path(TREND_DIR, "figures_temporal_diagnostics")
  timeseries_dir  <- file.path(TREND_DIR, "basin_averaged_timeseries")
  basin_plots_dir <- file.path(TREND_DIR, "basin_averaged_plots")
  for (d in c(fig_dir, timeseries_dir, basin_plots_dir))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  # ── Basin boundary ──────────────────────────────────────────────────
  basin_sf   <- load_basin(BASIN_SHP, EQUAL_AREA_CRS)
  basin_vect <- terra::vect(basin_sf)
  basin_sf_global <<- basin_sf   # expose to plot_raster_clean in utils
  # ── Analysis constants ──────────────────────────────────────────────
  DPI            <- 300
  TIMESCALES_VIZ <- TIMESCALES_SPATIAL    # c(1,3,6,12,24) for spatial figures
  TIMESCALES_TS  <- TIMESCALES_STANDARD   # c(1,3,6,9,12) for time series
  # Ensure SWEI colour is present
  if (!"swei" %in% names(index_colours))
    index_colours <- c(index_colours, swei = "#1E90FF")
  ####################################################################################
  # PART 1 – MERGE INDIVIDUAL RESULT FILES → MASTER CSV
  ####################################################################################
  cat("\n══════════════════════════════════════\n")
  cat("PART 1: Building master diagnostics file\n")
  cat("══════════════════════════════════════\n")
  master_csv <- file.path(TREND_DIR, "all_temporal_diagnostics_results.csv")
  load_clip_results <- function(index_type, timescale) {
    f <- file.path(TREND_DIR, sprintf("%s_%02d_results.csv", index_type, timescale))
    if (!file.exists(f)) {
      # Only warn when we expect the file (SWEI with wrong scale is silently skipped)
      if (!(index_type == "swei" && timescale != SWEI_SCALE))
        cat("  ⚠ Missing:", basename(f), "\n")
      return(NULL)
    }
    dt <- data.table::fread(f)
    cat(sprintf("  %s-%02d: %d rows\n", toupper(index_type), timescale, nrow(dt)))
    # Standardise coordinate column names
    if ("x" %in% names(dt) && !"lon" %in% names(dt)) data.table::setnames(dt, "x", "lon")
    if ("y" %in% names(dt) && !"lat" %in% names(dt)) data.table::setnames(dt, "y", "lat")
    if (!all(c("lon", "lat") %in% names(dt))) { cat("  ⚠ No coord columns\n"); return(NULL) }
    dt <- dt[!is.na(lon) & !is.na(lat)]
    tryCatch(
      clip_to_basin(dt, basin_sf, EQUAL_AREA_CRS),
      error = function(e) { cat("  ⚠ Clip error:", e$message, "\n"); NULL }
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
  if (!length(all_results)) stop("No data loaded. Run w3_trend_test.R first.")
  combined <- data.table::rbindlist(all_results, fill = TRUE)
  data.table::fwrite(combined, master_csv)
  cat(sprintf("✓ Master CSV: %d rows → %s\n", nrow(combined), basename(master_csv)))
  ####################################################################################
  # PART 2 – SPATIAL FIGURES
  ####################################################################################
  cat("\n══════════════════════════════════════\n")
  cat("PART 2: Spatial figures\n")
  cat("══════════════════════════════════════\n")
  plot_data <- combined
  # ── Shared raster helpers ─────────────────────────────────────────
  make_raster <- function(dt, vcol, tpl = NULL) {
    if (is.null(tpl)) tpl <- create_raster_template(dt, basin_sf)
    create_raster_from_points(dt, tpl, vcol, basin_sf)
  }
  # Helper: which indices are available for a given scale (includes SWEI when applicable)
  active_indices <- function(scale) {
    base <- c("spi", "spei")
    if (scale == SWEI_SCALE && any(plot_data$index_type == "swei")) base <- c(base, "swei")
    base
  }
  # ── Fig 1: Trend maps (Kendall τ_vc + significance) ──────────────
  create_trend_maps <- function(scale = 12, outfile) {
    cat(sprintf("  Fig 1: trend maps (scale=%d)...", scale))
    if (!safe_pdf(outfile, 12, 10)) return()
    indices_here <- active_indices(scale)
    n_plots      <- length(indices_here)
    par(mfrow = c(2, n_plots), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
    tau_range  <- c(-0.4, 0.4)
    sig_breaks <- c(0.5, 1.5, 2.5, 3.5, 4.5)
    for (idx in indices_here) {
      sub <- plot_data[index_type == idx & timescale == scale]
      if (!nrow(sub)) { plot.new(); text(0.5, 0.5, paste("No data:", toupper(idx))); next }
      tpl <- create_raster_template(sub, basin_sf)
      # τ_vc (variance-corrected)
      r_tau  <- make_raster(sub, "tau_vc", tpl)
      plot_raster_clean(r_tau, paste0(toupper(idx), ": Kendall's τ (VC)"),
                        tau_range, grDevices::hcl.colors(101, "RdBu", rev = TRUE),
                        legend_title = "τ")
      
      # Significance categories (based on p_value_vc)
      sub2  <- data.table::copy(sub)
      sub2[, sig_num := dplyr::case_when(
        is.na(p_value_vc)   ~ NA_real_,
        p_value_vc  < 0.001  ~ 1,
        p_value_vc  < 0.01   ~ 2,
        p_value_vc  < 0.05   ~ 3,
        TRUE                ~ 4)]
      r_sig  <- make_raster(sub2[!is.na(sig_num)], "sig_num", tpl)
      plot_raster_clean(r_sig, paste0(toupper(idx), ": Significance"),
                        NULL, c("#d73027", "#fc8d59", "#fee08b", "#f7f7f7"),
                        breaks = sig_breaks, categorical = TRUE,
                        legend_title = "p-value")
    }
    graphics::mtext(sprintf("Trend Analysis – VC Method  (%d-month)", scale),
                    outer = TRUE, cex = 1.2, font = 2)
    grDevices::dev.off()
    cat(" ✓\n")
  }
  # ── Fig 2: Event characteristic maps ────────────────────────────
  create_event_maps <- function(scale = 12, outfile) {
    cat(sprintf("  Fig 2: event maps (scale=%d)...", scale))
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
        else                        c(0, max(sub[[mn]], na.rm = TRUE))
        plot_raster_clean(r, paste0(toupper(idx), ": ", metrics[[mn]]$title), zlm,
                          grDevices::hcl.colors(101, metrics[[mn]]$pal,
                                                rev = (mn == "max_intensity")))
      }
    }
    graphics::mtext(sprintf("Event Characteristics  (%d-month)", scale),
                    outer = TRUE, cex = 1.2, font = 2)
    grDevices::dev.off()
    cat(" ✓\n")
  }
  # ── Fig 3: Temporal pattern maps ─────────────────────────────────
  # Requires columns: filtered_runs, p_value_runs, clustering,
  # regime_shift_year, n_spectral_peaks
  # (all produced by the upgraded w3_trend_test.R)
  create_temporal_maps <- function(scale = 12, outfile) {
    cat(sprintf("  Fig 3: temporal maps (scale=%d)...", scale))
    if (!safe_pdf(outfile, 12, 8)) return()
    indices_here <- active_indices(scale)
    n_plots      <- length(indices_here)
    par(mfrow = c(n_plots, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
    # Check that the required columns exist (graceful fallback if w3 was not fully upgraded)
    need_cols  <- c("filtered_runs", "p_value_runs", "clustering",
                    "regime_shift_year", "n_spectral_peaks")
    missing_cols  <- setdiff(need_cols, names(plot_data))
    if (length(missing_cols)) {
      plot.new()
      text(0.5, 0.5,
           paste("Missing columns – re-run w3_trend_test.R:\n",
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
      sub_cl  <- data.table::copy(sub)
      sub_cl[, clus_num := data.table::fcase(
        filtered_runs | is.na(p_value_runs), NA_real_,
        p_value_runs  >= 0.05,                 0,
        clustering == "clustered",           -1,
        clustering == "dispersed",            1)]
      r_cl  <- make_raster(sub_cl, "clus_num", tpl)
      plot_raster_clean(r_cl,
                        paste0(toupper(idx), ": Temporal Clustering"),
                        NULL, c("#e41a1c", "#fee090", "#4daf4a"),
                        breaks = c(-1.5, -0.5, 0.5, 1.5),
                        categorical = TRUE, legend = FALSE)
      graphics::legend("bottomright",
                       legend = c("Clustered", "Not sig.", "Dispersed"),
                       fill   = c("#e41a1c", "#fee090", "#4daf4a"),
                       bty = "n", cex = 0.85)
      
      # Panel B: Regime-shift decade
      sub_rs  <- sub[!is.na(regime_shift_year)]
      if (nrow(sub_rs)) {
        sub_rs2  <- data.table::copy(sub_rs)
        sub_rs2[, shift_dec := as.numeric(
          cut(regime_shift_year,
              breaks = seq(1950, 2030, by = 10),
              labels = 1:8, include.lowest = TRUE)) ]
        r_rs  <- make_raster(sub_rs2, "shift_dec", tpl)
        plot_raster_clean(r_rs,
                          paste0(toupper(idx), ": Regime-Shift Decade"),
                          c(1, 8),
                          grDevices::hcl.colors(8, "Spectral", rev = TRUE),
                          legend_title = "Decade")
      } else {
        plot.new(); text(0.5, 0.5, paste(toupper(idx), "– no regime shifts detected"))
      }
      
      # Panel C: Spectral peaks
      r_sp  <- make_raster(sub, "n_spectral_peaks", tpl)
      plot_raster_clean(r_sp,
                        paste0(toupper(idx), ": Spectral Peaks"),
                        c(0, max(sub$n_spectral_peaks, na.rm = TRUE, 1)),
                        grDevices::hcl.colors(101, "viridis"))
    }
    graphics::mtext(sprintf("Temporal Patterns  (%d-month)", scale),
                    outer = TRUE, cex = 1.2, font = 2)
    grDevices::dev.off()
    cat(" ✓\n")
  }
  # ── Fig 4: Timescale comparison (all scales, SPI + SPEI + SWEI) ──
  create_timescale_comparison <- function(outfile) {
    cat("  Fig 4: timescale comparison...")
    smry <- plot_data[, .(
      Pct_Sig     = sum(p_value_vc < 0.05, na.rm = TRUE) / .N * 100,
      Median_Tau  = median(tau_vc, na.rm = TRUE),
      Mean_Events = mean(n_events, na.rm = TRUE),
      Mean_Dur    = mean(mean_duration, na.rm = TRUE)
    ), by = .(index_type, timescale)]
    mk <- function(col, ylab, title) {
      ggplot2::ggplot(smry, ggplot2::aes(x = timescale, y = .data[[col]],
                                         color = index_type,
                                         group = index_type)) +
        ggplot2::geom_line(linewidth = 1.2, na.rm = TRUE) +
        ggplot2::geom_point(size = 4, na.rm = TRUE) +
        ggplot2::scale_color_manual(values = index_colours, labels = toupper) +
        ggplot2::scale_x_continuous(breaks = sort(unique(smry$timescale))) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(title = title, x = "Timescale (months)", y = ylab, color = "Index")
    }
    p  <- (mk("Pct_Sig",     "%",        "A) % Significant Pixels (VC)") |
             mk("Median_Tau",  "τ",        "B) Median Kendall's τ (VC)"))  /
      (mk("Mean_Events", "N events", "C) Mean Drought Events")   |
         mk("Mean_Dur",    "months",   "D) Mean Duration")) +
      patchwork::plot_layout(guides = "collect")  &
      ggplot2::theme(legend.position = "bottom")
    ggplot2::ggsave(outfile, p, width = 12, height = 10, dpi = DPI)
    cat(" ✓\n")
  }
  # ── Fig 5: Method comparison (VC vs TFPW) ────────────────────────
  # Requires columns tau_vc, p_value_vc, tau_tfpw, p_value_tfpw from w3.
  create_method_comparison <- function(outfile) {
    cat("  Fig 5: method comparison (VC vs TFPW)...")
    needed <- c("tau_tfpw", "p_value_tfpw")
    if (!all(needed %in% names(plot_data))) {
      if (safe_pdf(outfile)) {
        plot.new()
        text(0.5, 0.5,
             "TFPW columns not found.\nRe-run w3_trend_test.R with TFPW enabled.",
             cex = 1.2, col = "orange")
        grDevices::dev.off()
      }
      cat(" ⚠ TFPW columns absent – placeholder saved\n")
      return()
    }
    # Use 12-month scale for comparison (most stable signal)
    sub <- plot_data[timescale == 12]
    sub[, `:=`(vc_sig   = p_value_vc   < 0.05,
               tfpw_sig = p_value_tfpw < 0.05)]
    # Panel A: Scatter τ_vc vs τ_tfpw
    p1 <- ggplot2::ggplot(sub, ggplot2::aes(tau_vc, tau_tfpw,
                                            color = index_type)) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "gray40") +
      ggplot2::geom_point(alpha = 0.4, size = 2, na.rm = TRUE) +
      ggplot2::facet_wrap(~index_type,
                          labeller = ggplot2::labeller(index_type = toupper)) +
      ggplot2::scale_color_manual(values = index_colours) +
      ggplot2::coord_equal() +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = "A) τ: Variance-Corrected vs TFPW  (12-month)",
                    x = "τ (VC)", y = "τ (TFPW)")
    # Panel B: Method-agreement stacked bar per index
    agree <- sub[, .(
      Both   = sum( vc_sig &  tfpw_sig, na.rm = TRUE),
      VC_only= sum( vc_sig & !tfpw_sig, na.rm = TRUE),
      TFPW_only= sum(!vc_sig &  tfpw_sig, na.rm = TRUE),
      None   = sum(!vc_sig & !tfpw_sig, na.rm = TRUE)
    ), by = index_type]
    ag_l <- data.table::melt(agree, id.vars = "index_type",
                             variable.name = "Category",
                             value.name    = "N_pixels")
    p2 <- ggplot2::ggplot(ag_l,
                          ggplot2::aes(toupper(index_type), N_pixels,
                                       fill = Category)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = "B) Significance Agreement (p<0.05)  (12-month)",
                    x = "Index", y = "# Pixels")
    ggplot2::ggsave(outfile, p1 | p2, width = 14, height = 6, dpi = DPI)
    cat(" ✓\n")
  }
  # ── Run all spatial figures ───────────────────────────────────────
  for (sc in c(6, 12)) {
    sfx <- sprintf("%02dmo", sc)
    tryCatch(
      create_trend_maps(sc, file.path(fig_dir,
                                      sprintf("Fig1_TrendMaps_%s.pdf", sfx))),
      error = function(e) cat("  ❌ Fig 1 trend maps:", e$message, "\n"))
    tryCatch(
      create_event_maps(sc, file.path(fig_dir,
                                      sprintf("Fig2_EventMaps_%s.pdf", sfx))),
      error = function(e) cat("  ❌ Fig 2 event maps:", e$message, "\n"))
    tryCatch(
      create_temporal_maps(sc, file.path(fig_dir,
                                         sprintf("Fig3_TemporalMaps_%s.pdf", sfx))),
      error = function(e) cat("  ❌ Fig 3 temporal maps:", e$message, "\n"))
  }
  tryCatch(
    create_timescale_comparison(file.path(fig_dir, "Fig4_TimescaleComparison.png")),
    error = function(e) cat("  ❌ Fig 4 timescale comparison:", e$message, "\n"))
  tryCatch(
    create_method_comparison(file.path(fig_dir, "Fig5_MethodComparison.png")),
    error = function(e) cat("  ❌ Fig 5 method comparison:", e$message, "\n"))
  ####################################################################################
  # PART 3 – BASIN-AVERAGED TIME SERIES CSVs
  ####################################################################################
  cat("\n══════════════════════════════════════\n")
  cat("PART 3: Basin-averaged time series CSVs\n")
  cat("══════════════════════════════════════\n")
  extract_basin_avg_from_csv <- function(index_type, timescale) {
    f <- file.path(TREND_DIR,
                   sprintf("%s_%02d_results.csv", index_type, timescale))
    if (!file.exists(f)) { cat("  ⚠ Missing:", basename(f), "\n"); return(NULL) }
    dt  <- data.table::fread(f)
    if ("x" %in% names(dt) && !"lon" %in% names(dt)) data.table::setnames(dt, "x", "lon")
    if ("y" %in% names(dt) && !"lat" %in% names(dt)) data.table::setnames(dt, "y", "lat")
    # Clip to basin bounding box
    ext     <- terra::ext(basin_vect)
    dt_clip <- dt[lon >= ext[1] & lon <= ext[2] & lat >= ext[3] & lat <= ext[4]]
    if (!nrow(dt_clip)) { cat("  ⚠ No basin points\n"); return(NULL) }
    time_cols <- grep("^time_[0-9]{3}$", names(dt_clip), value = TRUE)
    if (!length(time_cols)) { cat("  ⚠ No time_XXX columns\n"); return(NULL) }
    ts_mat <- as.matrix(dt_clip[, ..time_cols])
    avg    <- colMeans(ts_mat, na.rm = TRUE)
    dates  <- seq(as.Date("1950-01-01"), by = "month", length.out = length(time_cols))
    df <- data.frame(date = dates, value = avg)
    df[!is.na(df$value), ]
  }
  # SPI & SPEI
  cat("\n  Processing SPI & SPEI...\n")
  for (idx in c("spi", "spei")) {
    for (sc in TIMESCALES_TS) {
      df <- extract_basin_avg_from_csv(idx, sc)
      if (!is.null(df)) {
        out <- file.path(timeseries_dir,
                         sprintf("%s_%02d_basin_average.csv", idx, sc))
        utils::write.csv(df, out, row.names = FALSE)
        cat(sprintf("  ✓ %s-%02d: %d months\n", toupper(idx), sc, nrow(df)))
      }
    }
  }
  # SWEI
  cat("\n  Processing SWEI...\n")
  df_swei <- extract_basin_avg_from_csv("swei", SWEI_SCALE)
  if (!is.null(df_swei)) {
    out <- file.path(timeseries_dir,
                     sprintf("swei_%02d_basin_average.csv", SWEI_SCALE))
    utils::write.csv(df_swei, out, row.names = FALSE)
    cat(sprintf("  ✓ SWEI-%02d: %d months\n", SWEI_SCALE, nrow(df_swei)))
  }
  ####################################################################################
  # PART 4 – TIME SERIES PLOTS
  ####################################################################################
  cat("\n══════════════════════════════════════\n")
  cat("PART 4: Time series plots\n")
  cat("══════════════════════════════════════\n")
  load_basin_ts_csv <- function(index_type, scale) {
    f  <- file.path(timeseries_dir,
                    sprintf("%s_%02d_basin_average.csv", index_type, scale))
    if (!file.exists(f)) stop("CSV not found: ", f)
    df        <- data.table::fread(f)
    df$date   <- as.Date(df$date)
    df$value  <- as.numeric(df$value)
    df
  }
  make_ts_plot  <- function(df, index_type, scale) {
    idx_color  <- index_colours[index_type]
    ggplot2::ggplot(df, ggplot2::aes(date, value)) +
      drought_band_layers() +
      ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed",
                          color = "gray50", linewidth = 0.3) +
      ggplot2::geom_line(color = idx_color, linewidth = 0.6) +
      ggplot2::geom_ribbon(
        data = subset(df, value < DROUGHT_ONSET),
        ggplot2::aes(ymin = value, ymax = DROUGHT_ONSET),
        fill = "coral", alpha = 0.3) +
      shared_ts_theme(14) +
      ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
      ggplot2::coord_cartesian(ylim = c(-3.5, 3.5)) +
      ggplot2::labs(
        title = sprintf("%s-%02d  Basin-Averaged Time Series",
                        toupper(index_type), scale),
        x = "Year",
        y = sprintf("%s-%02d", toupper(index_type), scale))
  }
  # SPI & SPEI
  cat("\n  Creating SPI & SPEI time series plots...\n")
  for (idx in c("spi", "spei")) {
    for (sc in TIMESCALES_TS) {
      tryCatch({
        df   <- load_basin_ts_csv(idx, sc)
        out  <- file.path(basin_plots_dir,
                          sprintf("%s_%02d_timeseries.png", idx, sc))
        ggplot2::ggsave(out, make_ts_plot(df, idx, sc), width = 12, height = 8, dpi = DPI)
        cat(sprintf("  ✓ %s-%02d\n", toupper(idx), sc))
      }, error = function(e)
        cat(sprintf("  ⚠ %s-%02d: %s\n", toupper(idx), sc, e$message)))
    }
  }
  # SWEI
  cat("\n  Creating SWEI time series plot...\n")
  tryCatch({
    df  <- load_basin_ts_csv("swei", SWEI_SCALE)
    out <- file.path(basin_plots_dir,
                     sprintf("swei_%02d_timeseries.png", SWEI_SCALE))
    ggplot2::ggsave(out, make_ts_plot(df, "swei", SWEI_SCALE),
                    width = 12, height = 8, dpi = DPI)
    cat(sprintf("  ✓ SWEI-%02d\n", SWEI_SCALE))
  }, error = function(e)
    cat(sprintf("  ⚠ SWEI-%02d: %s\n", SWEI_SCALE, e$message)))
  ####################################################################################
  # PART 5 – EXCEL SUMMARY
  ####################################################################################
  cat("\n══════════════════════════════════════\n")
  cat("PART 5: Excel summary\n")
  cat("══════════════════════════════════════\n")
  # Loader wrapper that clamps SWEI to its single scale
  load_any_ts <- function(idx, sc) {
    if (idx == "swei") sc <- SWEI_SCALE
    as.data.frame(load_basin_ts_csv(idx, sc))
  }
  export_summary_excel_all <- function(output_file) {
    wb        <- openxlsx::createWorkbook()
    hdr_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
    # ── Sheet 1: Summary Statistics ───────────────────────────────
    stats_df <- data.frame(Timescale = character(), Index = character(),
                           Mean = numeric(), Median = numeric(),
                           StdDev = numeric(), Min = numeric(), Max = numeric(),
                           Drought_Events_2020_2025 = integer(),
                           stringsAsFactors = FALSE)
    for (idx in c("spi", "spei")) {
      for (sc in TIMESCALES_TS) {
        tryCatch({
          df      <- load_any_ts(idx, sc)
          df_rec  <- df[df$date >= as.Date("2020-01-01") &
                          df$date <= as.Date("2025-12-31"), ]
          ev_count  <- nrow(detect_drought_events(df_rec))
          stats_df  <- rbind(stats_df, data.frame(
            Timescale                = sprintf("%02d", sc),
            Index                    = toupper(idx),
            Mean                     = mean(df$value,   na.rm = TRUE),
            Median                   = median(df$value, na.rm = TRUE),
            StdDev                    = sd(df$value,     na.rm = TRUE),
            Min                      = min(df$value,    na.rm = TRUE),
            Max                      = max(df$value,    na.rm = TRUE),
            Drought_Events_2020_2025 = ev_count,
            stringsAsFactors         = FALSE))
        }, error = function(e)
          cat(sprintf("  ⚠ Stats for %s-%02d: %s\n", toupper(idx), sc, e$message)))
      }
    }
    # SWEI row
    tryCatch({
      df      <- load_any_ts("swei", SWEI_SCALE)
      df_rec  <- df[df$date >= as.Date("2020-01-01") &
                      df$date <= as.Date("2025-12-31"), ]
      ev_count  <- nrow(detect_drought_events(df_rec))
      stats_df  <- rbind(stats_df, data.frame(
        Timescale                = sprintf("%02d", SWEI_SCALE),
        Index                    = "SWEI",
        Mean                     = mean(df$value,   na.rm = TRUE),
        Median                   = median(df$value, na.rm = TRUE),
        StdDev                   = sd(df$value,     na.rm = TRUE),
        Min                      = min(df$value,    na.rm = TRUE),
        Max                      = max(df$value,    na.rm = TRUE),
        Drought_Events_2020_2025 = ev_count,
        stringsAsFactors          = FALSE))
    }, error = function(e)
      cat(sprintf("  ⚠ Stats for SWEI-%02d: %s\n", SWEI_SCALE, e$message)))
    openxlsx::addWorksheet(wb, "Summary_Statistics")
    openxlsx::writeData(wb, "Summary_Statistics", stats_df)
    openxlsx::addStyle(wb, "Summary_Statistics", hdr_style,
                       rows = 1, cols = seq_len(ncol(stats_df)))
    # ── Sheet 2: Pairwise correlations at scale 3 (all index pairs) ──
    corr_df <- data.frame(Comparison   = character(),
                          Correlation  = numeric(),
                          stringsAsFactors = FALSE)
    tryCatch({
      s_df  <- load_any_ts("spi",  3)
      sp_df <- load_any_ts("spei", 3)
      sw_df <- load_any_ts("swei", 3)
      mrg  <- merge(
        merge(s_df, sp_df, by = "date", suffixes = c("_spi", "_spei")),
        sw_df, by = "date")
      names(mrg)[names(mrg) == "value"] <- "swei"
      
      corr_df  <- rbind(corr_df,
                        data.frame(Comparison = "SPI-3 vs SPEI-3",
                                   Correlation = round(cor(mrg$value_spi, mrg$value_spei,
                                                           use = "complete.obs"), 4),
                                   stringsAsFactors = FALSE),
                        data.frame(Comparison = "SPI-3 vs SWEI-3",
                                   Correlation = round(cor(mrg$value_spi, mrg$swei,
                                                           use = "complete.obs"), 4),
                                   stringsAsFactors = FALSE),
                        data.frame(Comparison = "SPEI-3 vs SWEI-3",
                                   Correlation = round(cor(mrg$value_spei, mrg$swei,
                                                           use = "complete.obs"), 4),
                                   stringsAsFactors = FALSE))
    }, error = function(e) cat("  ⚠ Correlation failed:", e$message, "\n"))
    openxlsx::addWorksheet(wb, "Correlations_Scale3")
    openxlsx::writeData(wb, "Correlations_Scale3", corr_df)
    openxlsx::addStyle(wb, "Correlations_Scale3", hdr_style,
                       rows = 1, cols = seq_len(ncol(corr_df)))
    openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
    cat(sprintf("✓ Excel summary saved: %s\n", basename(output_file)))
  }
  excel_out <- file.path(timeseries_dir, "Drought_Summary_Statistics.xlsx")
  export_summary_excel_all(excel_out)
  ####################################################################################
  # FINAL SUMMARY
  ####################################################################################
  cat("\n╔══════════════════════════════════════╗\n")
  cat("║  w4_trends_visualization.R  DONE     ║\n")
  cat("╚══════════════════════════════════════╝\n\n")
  cat("Outputs:\n")
  cat("  Spatial figures  : ", normalizePath(fig_dir),          "\n")
  cat("  Time series CSVs : ", normalizePath(timeseries_dir),   "\n")
  cat("  TS plots         : ", normalizePath(basin_plots_dir),  "\n")
  cat("  Excel summary    : ", normalizePath(excel_out),        "\n")
  cat("\n✓ All visualization tasks complete!\n\n")