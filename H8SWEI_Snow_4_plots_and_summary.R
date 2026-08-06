##############################################
# 04_plots_and_summary.R
# STAGE 4 of 4 — Basin-averaged summary CSVs, time-series plots, and the
#                 final text summary / fallback-alerts report
#
# This is the cheap, cosmetic stage - if you want to tweak plot formatting
# or wording in the summary file, you can re-run just this script in
# seconds without touching anything upstream.
#
# INPUT:  stage_checkpoints/stage3.rds  (from 03_main_parallel_calc.R)
# OUTPUT: final deliverables only, written into out_dir (basin-averaged
#         CSVs, out_dir/figures/*.png, sspi_summary.txt /
#         swei_all_timescales_summary.txt, PIPELINE_FALLBACK_ALERTS.txt).
#         No further .rds checkpoint is produced - this is the last stage.
##############################################

STAGE_DIR <- "stage_checkpoints"
stage3_path <- file.path(STAGE_DIR, "stage3.rds")
if (!file.exists(stage3_path)) {
  stop(sprintf(
    "Cannot find %s. Run 01_setup_and_data.R, 02_regionalization_bias.R, and 03_main_parallel_calc.R first.",
    stage3_path
  ), call. = FALSE)
}

cat("============================================================\n")
cat("STAGE 4 - loading checkpoint from script 03\n")
cat("============================================================\n")
stage3_data <- readRDS(stage3_path)
list2env(stage3_data, envir = environment())
rm(stage3_data)
cat(sprintf("[OK] Restored %d object(s) from %s\n", length(ls()), stage3_path))

# ==============================================================================
#   SAVE BASIN-AVERAGED SWEI TO CSV 
# ==============================================================================
index_lbl <- if (scf_method == "2") "sspi" else "swei"
cat(sprintf("\n===== SAVING BASIN-AVERAGED %s =====\n", toupper(index_lbl)))
months_all  <- as.integer(format(dates, "%m"))
years_all   <- as.integer(format(dates, "%Y"))
month_names_save  <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
for (sc in timescales) {
  basin_series <- basin_avg_swei_results[[as.character(sc)]]
  yrs <- sort(unique(years_all))
  df_out <- data.frame(Year = yrs)
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    yr_m  <- years_all[idx_m]
    val_m <- basin_series[idx_m]
    col   <- setNames(val_m, yr_m)
    df_out[[month_names_save[m]]] <- col[as.character(yrs)]
  }
  csv_file <- file.path(out_dir, sprintf("%s_%02d_basin_averaged_by_month.csv", index_lbl, sc))
  write.csv(df_out, csv_file, row.names = FALSE, na = "")
  cat(sprintf("[OK] Saved basin-averaged %s-%d (12 monthly series) to: %s\n", toupper(index_lbl), sc, csv_file))
}

# ==============================================================================
#   TIME SERIES PLOTS (BASIN-AVERAGED)
# ==============================================================================
# Two kinds of figures, written to out_dir/figures/:
#   1. plot_index_timeseries()   - basin-mean SWEI/SSPI-1 series per timescale,
#      SPI/SSI-style: standardized-index line, red fill for excursions below
#      zero (split into separate polygons per contiguous dry spell so the
#      ribbon doesn't visually bridge unrelated events), dashed/dotted
#      reference lines at +/-1/1.5/2, soft severity background bands, and a
#      subtitle reporting the number of drought events at the chosen onset
#      threshold. Mirrors plot_ssi_timeseries() from the streamflow SSI script
#      for a consistent look across indices.
#   2. plot_raw_swe_timeseries() - basin-mean RAW SWE (mm), the physical
#      variable the index is standardized from, for context alongside the
#      standardized series.
figures_dir <- file.path(out_dir, "figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

plot_index_timeseries <- function(dates_vec, values_vec, sc, index_lbl,
                                  index_full_name, title_label, out_dir, onset,
                                  method_label) {
  df <- data.frame(date = dates_vec, value = values_vec)
  df <- df[is.finite(df$value), ]
  if (nrow(df) == 0L) {
    cat(sprintf("  [PLOT] Skipping %s-%d time series: no finite values\n",
                toupper(index_lbl), sc))
    return(invisible(NULL))
  }
  
  # Drought event count (contiguous runs below onset threshold)
  below_onset <- df$value < onset
  below_onset[is.na(below_onset)] <- FALSE
  n_events <- sum(rle(below_onset)$values)
  
  # Shade only the runs that actually cross the onset threshold, so the red
  # fill visually matches the events counted above -- not just any negative
  # month (which would shade e.g. a -0.3 dip identically to a -2.5 deficit
  # even though only the latter is a counted event). Separate contiguous
  # runs so the ribbon fill doesn't bridge unrelated dry spells.
  df$below <- below_onset
  df$grp   <- cumsum(c(1L, diff(as.integer(df$below)) != 0L))
  
  y_lo <- min(-2.2, floor(min(df$value) * 10) / 10)
  y_hi <- max( 2.2, ceiling(max(df$value) * 10) / 10)
  
  p <- ggplot(df, aes(x = date, y = value)) +
    
    # ---- severity background bands ----
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1,    ymax = Inf,
           fill = "#4472C4", alpha = 0.10) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,    ymax = 1,
             fill = "#70AD47", alpha = 0.08) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1,   ymax = 0,
             fill = "#70AD47", alpha = 0.05) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -1,
             fill = "#C00000", alpha = 0.08) +
    
    # ---- reference lines ----
  geom_hline(yintercept = 0,             color = "black", linewidth = 0.6) +
    geom_hline(yintercept = c(-1, 1),      color = "grey30", linetype = "dashed",
               linewidth = 0.4) +
    geom_hline(yintercept = c(-1.5, 1.5),  color = "grey55", linetype = "dotted",
               linewidth = 0.4) +
    geom_hline(yintercept = c(-2, 2),      color = "grey55", linetype = "dotted",
               linewidth = 0.4) +
    
    # ---- drought shading: ONLY months that cross the onset threshold,
    #      but filled all the way up to zero (not stopped at onset) so
    #      the shaded region reads as a continuous deficit down to the
    #      curve, not a band truncated mid-air at the onset line ----
  geom_ribbon(data = subset(df, below),
              aes(ymin = value, ymax = 0, group = grp),
              fill = "#C0504D", alpha = 0.45) +
    
    # ---- the series itself ----
  geom_line(color = "#1F77B4", linewidth = 0.45) +
    
    scale_x_date(date_breaks = "5 years", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.015))) +
    scale_y_continuous(breaks = seq(-4, 4, 1)) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    
    labs(
      title    = sprintf("%s Time Series  |  NECHAKO BASIN (basin-mean)",
                         title_label),
      subtitle = sprintf("%s  |  %d snow-drought events detected  (onset < %.2f)",
                         method_label, n_events, onset),
      x = "Year", y = sprintf("%s Index Value", index_full_name)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_line(color = "grey85"),
      panel.grid.major.y  = element_blank(),
      plot.title          = element_text(face = "bold", size = 15),
      plot.subtitle       = element_text(color = "grey40", size = 11),
      axis.text.x         = element_text(angle = 30, hjust = 1),
      plot.margin         = margin(10, 16, 10, 10)
    )
  
  fname <- file.path(out_dir, "figures",
                     sprintf("%s_%02d_basin_timeseries.png", index_lbl, sc))
  ggsave(fname, p, width = 11, height = 5.5, dpi = 200, bg = "white")
  cat(sprintf("  [PLOT] %s-%d basin-mean time series saved: %s  (%d events)\n",
              toupper(index_lbl), sc, basename(fname), n_events))
  invisible(p)
}

plot_raw_swe_timeseries <- function(dates_vec, swe_vec, out_dir) {
  df <- data.frame(date = dates_vec, swe = swe_vec)
  df <- df[is.finite(df$swe), ]
  if (nrow(df) == 0L) {
    cat("  [PLOT] Skipping raw basin-mean SWE time series: no finite values\n")
    return(invisible(NULL))
  }
  
  p <- ggplot(df, aes(x = date, y = swe)) +
    geom_area(fill = "#4472C4", alpha = 0.20) +
    geom_line(color = "#2E5A9C", linewidth = 0.45) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.015))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = "Basin-Averaged Raw SWE  |  NECHAKO BASIN",
      subtitle = "ERA5-Land snow water equivalent, spatial mean over basin pixels (pre-standardization)",
      x = "Year", y = "SWE (mm)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_line(color = "grey85"),
      plot.title          = element_text(face = "bold", size = 15),
      plot.subtitle       = element_text(color = "grey40", size = 11),
      axis.text.x         = element_text(angle = 30, hjust = 1),
      plot.margin         = margin(10, 16, 10, 10)
    )
  
  fname <- file.path(out_dir, "figures", "raw_swe_basin_timeseries.png")
  ggsave(fname, p, width = 11, height = 5.5, dpi = 200, bg = "white")
  cat(sprintf("  [PLOT] Raw basin-mean SWE time series saved: %s\n", basename(fname)))
  invisible(p)
}

cat("\n===== GENERATING TIME SERIES PLOTS =====\n")

# Onset threshold for the drought-event count in the subtitle: -0.5, matching
# SSI_DROUGHT_THRESHOLD in the companion streamflow SSI script (H7WSC_
# StreamFlow.R), so "event" means the same thing (index crosses below -0.5)
# across indices rather than mixing an onset convention here with a
# McKee-style severity-tier value there.
index_onset      <- -0.5
index_full_name  <- if (scf_method == "2") "SSPI-1" else "SWEI"

# Plot-only method label: scf_method_label (used elsewhere for run_settings.txt
# / the text summary) includes a "(zero-inflated gamma, no SCF mask)"
# parenthetical that's implementation detail, not something the figure
# subtitle needs -- kept short here instead.
plot_method_label <- if (scf_method == "2") "SSPI-1" else scf_method_label

for (sc in timescales) {
  # SSPI-1's own name already encodes its timescale ("-1"), so don't also
  # append "-%d" (which produced "SSPI-1-1"); SWEI's name doesn't, so it
  # still needs "-3" etc. appended to identify the timescale.
  title_label <- if (scf_method == "2") index_full_name else sprintf("%s-%d", index_full_name, sc)
  
  plot_index_timeseries(
    dates_vec        = dates,
    values_vec        = basin_avg_swei_results[[as.character(sc)]],
    sc                = sc,
    index_lbl         = index_lbl,
    index_full_name   = index_full_name,
    title_label       = title_label,
    out_dir           = out_dir,
    onset             = index_onset,
    method_label      = plot_method_label
  )
}

# Raw SWE is only pre-computed for method 1 (swe_basin_avg_raw); recompute
# from swe_matrix for method 2 so the plot is available for both methods.
swe_basin_avg_for_plot <- if (!is.null(swe_basin_avg_raw)) {
  swe_basin_avg_raw
} else {
  colMeans(swe_matrix, na.rm = TRUE)
}
plot_raw_swe_timeseries(dates, swe_basin_avg_for_plot, out_dir)

# ==============================================================================
#   WRITE SUMMARY FILE AT END 
# ==============================================================================
cat("\n===== WRITING SUMMARY FILE =====\n")

if (scf_method == "2") {
  cat(
    "============================================================\n",
    "MONTHLY SSPI-1 SUMMARY - NECHAKO BASIN\n",
    "ERA5-Land Snow Water Equivalent\n",
    "Method 2: Zero-Inflated Gamma (no SCF mask)\n",
    sprintf("Regional pooling (index-flood RFA): %s\n", if (REGIONALIZE_SSPI) "ENABLED" else "disabled"),
    "============================================================\n\n",
    file = summary_combined_file, sep = ""
  )
  for (sc in timescales) {
    s <- all_summaries[[as.character(sc)]]
    cat(sprintf("\n========== SSPI-%d ==========\n", s$sc), file = summary_combined_file, append = TRUE)
    cat(sprintf("Calculation date: %s\n", s$date), file = summary_combined_file, append = TRUE)
    cat(sprintf("Grid dimensions: %d x %d\n", ncol(swe), nrow(swe)), file = summary_combined_file, append = TRUE)
    cat(sprintf("Basin pixels (total): %d\n", basin_pixels), file = summary_combined_file, append = TRUE)
    cat(sprintf("Time period: %s to %s (%d months)\n", min(dates), max(dates), length(dates)), file = summary_combined_file, append = TRUE)
    cat(sprintf("Reference period: %s to %s\n", sspi_ref_start, sspi_ref_end), file = summary_combined_file, append = TRUE)
    cat("\nZero-Inflation Screening:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Excluded pixels (>50%% zero SWE in ref period): %d (%.1f%%)\n", length(invalid_pix_sspi), 100 * length(invalid_pix_sspi) / n_pixels), file = summary_combined_file, append = TRUE)
    cat("\nMethodology:\n", file = summary_combined_file, append = TRUE)
    if (REGIONALIZE_SSPI) {
      cat(sprintf("  - Regionalization: %s, %d homogeneous sub-regions (index-flood RFA, Fontana et al. 2026)\n",
                  regions_built$method, n_regions), file = summary_combined_file, append = TRUE)
      cat(sprintf("  - Per-pixel index value: %s of reference-period monthly SWE\n", INDEX_VALUE_STAT),
          file = summary_combined_file, append = TRUE)
      cat(sprintf("  - Discordancy screening (Hosking & Wallis D_i): %s, %d pixel-months excluded from pooling\n",
                  ifelse(DISCORDANCY_SCREENING, "ON", "OFF"), sum(discordant_mat)),
          file = summary_combined_file, append = TRUE)
      if (nrow(heterogeneity_log) > 0) {
        cat(sprintf("  - Heterogeneity (H0, Monte Carlo, Nsim=%d): %d/%d region-months flagged H0>=%.0f\n",
                    HETEROGENEITY_NSIM, sum(heterogeneity_log$H0 >= H0_REGION_FLAG_THRESHOLD, na.rm = TRUE),
                    nrow(heterogeneity_log), H0_REGION_FLAG_THRESHOLD),
            file = summary_combined_file, append = TRUE)
      }
      if (nrow(gof_log) > 0) {
        cat(sprintf("  - Pooled regional gamma GoF (KS, alpha=%.2f): %d/%d region-months pass\n",
                    GOF_ALPHA, sum(gof_log$gamma_used, na.rm = TRUE), nrow(gof_log)),
            file = summary_combined_file, append = TRUE)
      }
      cat(sprintf("  - Bias correction vs. station obs: %s (tier applied: %s; %d/12 Tier-1 months corrected; see rfa_diagnostics/)\n",
                  ifelse(BIAS_CORRECT && file.exists(STATION_OBS_PATH), "ENABLED", "disabled (no station file)"),
                  ifelse(exists("chosen_tier") && !is.na(chosen_tier), chosen_tier, "n/a"),
                  sum(basin_bias_log$applied, na.rm = TRUE)),
          file = summary_combined_file, append = TRUE)
      cat("  - Bias-correction validation (random k-fold, spatial-block LOSO, temporal transfer):\n    rfa_diagnostics/bias_correction_validation.csv\n",
          file = summary_combined_file, append = TRUE)
      cat(sprintf("  - Snow-cover-timing check vs. OpenLandMap climatology: %s\n",
                  ifelse(isTRUE(SCF_CLIMATOLOGY_VALIDATE), "see rfa_diagnostics/scf_climatology_validation.csv", "disabled")),
          file = summary_combined_file, append = TRUE)
      cat("  - Per-pixel zero-probability retained (only the gamma shape/scale is pooled regionally)\n",
          file = summary_combined_file, append = TRUE)
      cat("  - Region/months failing GoF, OR failing the H0 heterogeneity screen, are excluded\n    from pooling entirely and fall back to per-pixel local fitting (no regional-\n    empirical CDF fallback is used)\n",
          file = summary_combined_file, append = TRUE)
      cat("  - Pixels with no usable region/index value fall back to the original per-pixel fit\n",
          file = summary_combined_file, append = TRUE)
      cat(sprintf("  - Diagnostics written to: %s\n", file.path(out_dir, "rfa_diagnostics")),
          file = summary_combined_file, append = TRUE)
    } else {
      cat("  - No SCF domain mask applied (SSPI-1 is valid for all months)\n  - Calendar-month stratified fitting\n  - Zero-inflated gamma (center of probability mass)\n  - L-moments (lmomco) for gamma parameter estimation\n  - Empirical fallback when gamma fit fails or variance too low\n  - SSPI values clipped at +/-5\n", file = summary_combined_file, append = TRUE)
    }
    if (REGIONALIZE_SSPI) {
      cat("\nMethod Distribution:\n", file = summary_combined_file, append = TRUE)
      cat(sprintf("  Regional gamma:             %.1f%%\n", s$regional_gamma_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Regional empirical (unused, no fallback CDF): %.1f%%\n", s$regional_empirical_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Local gamma (fallback):      %.1f%%\n", s$local_gamma_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Local empirical (fallback):  %.1f%%\n", s$local_empirical_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Local zero-variance (fallback): %.1f%%\n", s$local_zerovar_pct), file = summary_combined_file, append = TRUE)
    } else {
      cat("\nMethod Distribution:\n", file = summary_combined_file, append = TRUE)
      cat(sprintf("  Gamma fitting:   %.1f%%\n", s$gamma_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Empirical:       %.1f%%\n", s$empirical_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Zero-variance:   %.1f%%\n", s$zero_var_pct), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Excluded:        %.1f%%\n", s$excluded_pct), file = summary_combined_file, append = TRUE)
    }
    cat("\nNA Analysis:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Total NA values in basin: %.3f%%\n", s$na_rate), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Compatibility: %s\n", ifelse(s$na_rate < 0.5, "READY (<0.5% NAs)", "CAUTION (>0.5% NAs)")), file = summary_combined_file, append = TRUE)
    fi <- 1:min(12, ncol(s$sspi_indices))
    fs <- s$sspi_indices[, fi, drop = FALSE]
    cat("\nSnow drought frequency (first 12 months):\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional deficit (SSPI < -2.0): %.1f%%\n", 100 * mean(fs < -2.0, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Extreme deficit     (SSPI < -1.5): %.1f%%\n", 100 * mean(fs < -1.5, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Moderate deficit    (SSPI < -1.0): %.1f%%\n", 100 * mean(fs < -1.0, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional surplus (SSPI >  2.0): %.1f%%\n", 100 * mean(fs >  2.0, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
  }
} else {
  cat(
    "============================================================\n",
    "COMBINED SWEI SUMMARY FOR ALL TIMESCALES\n",
    "Seasonal Approach: SWE rolling sum = timescale months (Option B)\n",
    "Domain mask window fixed at k=3 per Huning & AghaKouchak (2020)\n",
    "============================================================\n\n",
    file = summary_combined_file, sep = ""
  )
  for (sc in timescales) {
    s  <- all_summaries[[as.character(sc)]]
    cat(sprintf("\n\n========== SWEI-%d ==========\n", s$sc), file = summary_combined_file, append = TRUE)
    cat(sprintf("Calculation date: %s\n", s$date), file = summary_combined_file, append = TRUE)
    cat(sprintf("Grid dimensions: %d x %d\n", ncol(swe), nrow(swe)), file = summary_combined_file, append = TRUE)
    cat(sprintf("Basin pixels: %d\n", basin_pixels), file = summary_combined_file, append = TRUE)
    cat(sprintf("Time period: %s to %s (%d months)\n", min(dates), max(dates), length(dates)), file = summary_combined_file, append = TRUE)
    cat("\nMethodology:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  - %d-month rolling sum applied to RAW SWE (window = sc = %d)\n", sc, sc), file = summary_combined_file, append = TRUE)
    cat("  - SCF mask window fixed at k=3 (intentionally independent of SWE window,\n    per Huning & AghaKouchak 2020 - see target_k comment in script)\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  - Domain-mask source: %s\n", if (exists("domain_mask_source", inherits = FALSE)) domain_mask_source else "SCF"), file = summary_combined_file, append = TRUE)
    cat(sprintf("    -> Method-1 criterion: %.0f%% nonzero years (SWE fallback)\n",
                100 * SWE_DOMAIN_NONZERO_FRACTION), file = summary_combined_file, append = TRUE)
    cat("  - Gringorten plotting position (non-parametric)\n  - Zero-SWE perturbation: random uniform in (0, min_nonzero) per Huning & AghaKouchak (2020)\n  - Minimum data requirement: >= 75% of years for a calendar month must have nonzero SWE\n  - No clipping applied to SWEI values (full distribution retained)\n", file = summary_combined_file, append = TRUE)
    cat("\nMethod Distribution:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Gringorten fitting: %.1f%%\n", s$gringorten_pct), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Masked/NA: %.1f%%\n", s$masked_pct), file = summary_combined_file, append = TRUE)
    cat("\nNA Analysis:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Basin mean series:     %.1f%% complete (%.1f%% NA)\n", s$basin_mean_complete, 100 - s$basin_mean_complete), file = summary_combined_file, append = TRUE)
    
    # Safely print na_rate even if NA
    if (is.na(s$na_rate)) {
      cat("  Per-pixel gridded:     UNKNOWN% NA  (0 basin pixels)\n", file = summary_combined_file, append = TRUE)
      cat("  Gridded compatibility: UNKNOWN (0 valid basin pixels)\n", file = summary_combined_file, append = TRUE)
    } else {
      cat(sprintf("  Per-pixel gridded:     %.3f%% NA  (basin pixels only)\n", s$na_rate), file = summary_combined_file, append = TRUE)
      cat(sprintf("  Gridded compatibility: %s\n", ifelse(s$na_rate < 0.5, "READY  (<0.5% NAs)", "CAUTION (>0.5% NAs)")), file = summary_combined_file, append = TRUE)
    }
    
    if (!is.na(s$na_rate) && s$na_rate >= 0.5) {
      cat("  NOTE: The large gap between basin-mean completeness and per-pixel NA\n  is expected - spatial averaging masks individual pixel gaps, while the\n  75%%-nonzero-years criterion and SCF mask applied pixel-by-pixel are\n  much stricter. The basin mean series is suitable for time-series analysis;\n  use gridded files for spatial analysis with awareness of spatial coverage.\n", file = summary_combined_file, append = TRUE)
    }
    first_idx  <- 1:min(12, ncol(s$swei_indices))
    first_swei  <- s$swei_indices[, first_idx, drop = FALSE]
    cat("\nDrought frequency (first 12 months):\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional (SWEI < -2.0): %.1f%%\n", 100 * mean(first_swei < -2.0, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Extreme     (SWEI < -1.6): %.1f%%\n", 100 * mean(first_swei < -1.6, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Severe      (SWEI < -1.3): %.1f%%\n", 100 * mean(first_swei < -1.3, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Moderate    (SWEI < -0.8): %.1f%%\n", 100 * mean(first_swei < -0.8, na.rm = TRUE)), file = summary_combined_file, append = TRUE)
  }
}

## ---- Fallback alerts file: ALWAYS written, even when empty ----
## An empty-but-present file is a positive confirmation that no silent
## degrade happened this run; a non-empty file is a to-do list of exactly
## which inputs to fix before trusting the results. Never skip writing this
## just because nothing happened.
alerts_path <- file.path(out_dir, "PIPELINE_FALLBACK_ALERTS.txt")
if (length(PIPELINE_FALLBACK_ALERTS) == 0) {
  writeLines(sprintf(
    "No fallbacks triggered this run (%s). DEM and station-observation inputs (if BIAS_CORRECT/REGIONALIZE_SSPI were enabled) resolved as configured.",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")), alerts_path)
} else {
  writeLines(c(sprintf("%d fallback(s) triggered this run - review before trusting results:", length(PIPELINE_FALLBACK_ALERTS)),
               "", PIPELINE_FALLBACK_ALERTS), alerts_path)
}

cat("\n============================================================\n")
cat(if (scf_method == "2") "SSPI-1 CALCULATION COMPLETE!\n" else "SWEI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat(sprintf("Summary file: %s\n", summary_combined_file))
cat(sprintf("Time series plots: %s\n", file.path(out_dir, "figures")))
if (length(PIPELINE_FALLBACK_ALERTS) > 0) {
  cat(sprintf("\n[!!] %d FALLBACK ALERT(S) TRIGGERED THIS RUN - see: %s\n",
              length(PIPELINE_FALLBACK_ALERTS), alerts_path))
  cat("     Results may be using a lesser method (spatial k-means regions and/or\n")
  cat("     no bias correction) than configured. Review before treating as final.\n")
} else {
  cat(sprintf("\n[OK] No fallback alerts this run (see %s).\n", alerts_path))
}

SCRIPT_END_TIME     <- Sys.time()
SCRIPT_TOTAL_SECONDS <- as.numeric(difftime(SCRIPT_END_TIME, SCRIPT_START_TIME, units = "secs"))
cat("\n============================================================\n")
cat(sprintf("TOTAL RUN TIME: %s (%.1f minutes / %.2f hours)\n",
            format_hms(SCRIPT_TOTAL_SECONDS),
            SCRIPT_TOTAL_SECONDS / 60,
            SCRIPT_TOTAL_SECONDS / 3600))
cat(sprintf("Started:  %s\n", format(SCRIPT_START_TIME, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("Finished: %s\n", format(SCRIPT_END_TIME,   "%Y-%m-%d %H:%M:%S")))
cat("============================================================\n")

cat("\n--- READY FOR ANALYSIS ---\n")