####################################################################################
# PUBLICATION-QUALITY TREND MAPS FOR PR/PET
# ENHANCED WITH:
# - Temporal clustering analysis (runs test)
# - Regime shift detection
# - Spectral analysis
# - Method comparison plots
# - Comprehensive statistical summaries
####################################################################################

library(terra)
library(data.table)
library(sf)
library(zoo)
library(ggplot2)
library(patchwork)
library(scales)

setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "trend_analysis_pr_pet"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ===== LOAD CORE DATA =====
cat("📦 Loading results, basin boundary, and metadata...\n")
all_results   <- readRDS(file.path("trend_analysis_pr_pet", "all_results.rds"))
basin_boundary <- readRDS(file.path("trend_analysis_pr_pet", "basin_boundary.rds"))
metadata <- readRDS(file.path("trend_analysis_pr_pet", "analysis_metadata.rds"))

basin_avg_monthly <- metadata$basin_avg_monthly
basin_avg_annual <- metadata$basin_avg_annual

# ===== SANITY CHECKS =====
cat("\n📊 Basin average monthly precipitation range:",
    range(basin_avg_monthly$precip_mm_month, na.rm = TRUE), "mm/month\n")
cat("📊 Basin average monthly PET range:",
    range(basin_avg_monthly$pet_mm_month, na.rm = TRUE), "mm/month\n")
if (max(basin_avg_monthly$precip_mm_month, na.rm = TRUE) < 50) {
  warning("⚠️ Maximum monthly precipitation is below 50 mm – check units")
}

# ===== COLUMN AVAILABILITY CHECK =====
cat("\n📋 Checking available columns in results...\n")
all_cols <- names(all_results)
cat("   Total columns:", length(all_cols), "\n")

# Check for specific analysis columns
analysis_cols <- list(
  Basic = c("tau_vc", "p_value_vc", "sl_vc", "tau_tfpw", "p_value_tfpw"),
  Autocorr = c("rho1_vc", "rho1_tfpw"),
  Spectral = c("n_spectral_peaks", "dominant_period", "spectral_confidence"),
  Changepoint = c("changepoint_detected", "first_changepoint_year", "n_changepoints"),
  RunsTest = c("filtered_runs", "clustering", "p_value_runs")
)

for (analysis_name in names(analysis_cols)) {
  cols_to_check <- analysis_cols[[analysis_name]]
  available <- cols_to_check[cols_to_check %in% all_cols]
  missing <- cols_to_check[!cols_to_check %in% all_cols]
  
  if (length(available) > 0) {
    cat(sprintf("   ✅ %s: %s\n", analysis_name, paste(available, collapse=", ")))
  }
  if (length(missing) > 0) {
    cat(sprintf("   ⚠️  %s (missing): %s\n", analysis_name, paste(missing, collapse=", ")))
  }
}

cat("\n")

# ===== CREATE NATIVE-RESOLUTION TEMPLATE =====
pet_annual_sample <- all_results[variable == "PET" & period == "annual" & !filtered_vc][1:100]
dx <- median(diff(sort(unique(pet_annual_sample$x))), na.rm = TRUE)
dy <- median(diff(sort(unique(pet_annual_sample$y))), na.rm = TRUE)
resolution_bc <- min(dx, dy, na.rm = TRUE)

basin_ext <- ext(basin_boundary)
template_bc <- rast(ext = basin_ext, 
                    resolution = resolution_bc, 
                    crs = "EPSG:3005",
                    names = "trend")

cat(sprintf("✓ Native-resolution template (%.1f m, %d x %d cells)\n", 
            resolution_bc, ncol(template_bc), nrow(template_bc)))

# ===== HELPER FUNCTIONS =====
create_raster_from_table <- function(results_dt, template, value_col) {
  pts <- vect(results_dt, geom = c("x", "y"), crs = crs(template))
  r <- rasterize(pts, template, field = value_col, touches = TRUE)
  mask(r, vect(basin_boundary))
}

smooth_raster <- function(r) {
  if (all(is.na(values(r)))) return(r)
  focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
}

plot_raster_panel <- function(r, main, zlim = NULL, col, 
                              breaks = NULL, legend = TRUE, categorical = FALSE,
                              legend_title = "") {
  if (!categorical) r <- smooth_raster(r)
  
  plg_args <- list(cex = 0.9)
  if (legend_title != "") plg_args$title <- legend_title
  
  if (!is.null(breaks)) {
    plot(r, main = main, cex.main = 1.0, col = col, breaks = breaks,
         axes = FALSE, box = FALSE, legend = legend, plg = plg_args,
         colNA = NA)
  } else {
    plot(r, main = main, cex.main = 1.0, col = col, zlim = zlim,
         axes = FALSE, box = FALSE, legend = legend, plg = plg_args,
         colNA = NA)
  }
  plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
}

# ===== PREPARE ANNUAL RESULTS =====
cat("⚙️  Preparing trend metrics...\n")
pet_annual    <- all_results[variable == "PET" & period == "annual" & !filtered_vc & !is_basin_average]
precip_annual <- all_results[variable == "Precipitation" & period == "annual" & !filtered_vc & !is_basin_average]

pet_annual[, combined_vc := fifelse(p_value_vc < 0.05, abs(tau_vc), 0)]
pet_annual[, direction_vc := fifelse(p_value_vc < 0.05 & tau_vc < 0, -1,
                                     fifelse(p_value_vc < 0.05 & tau_vc > 0, 1, 0))]
pet_annual[, combined_tfpw := fifelse(p_value_tfpw < 0.05, abs(tau_tfpw), 0)]
pet_annual[, direction_tfpw := fifelse(p_value_tfpw < 0.05 & tau_tfpw < 0, -1,
                                       fifelse(p_value_tfpw < 0.05 & tau_tfpw > 0, 1, 0))]

precip_annual[, combined_vc := fifelse(p_value_vc < 0.05, tau_vc, 0)]
precip_annual[, combined_tfpw := fifelse(p_value_tfpw < 0.05, tau_tfpw, 0)]

# ===== GENERATE RASTERS =====
cat("🖼️  Generating rasters...\n")
r_pet_comb_vc <- create_raster_from_table(pet_annual, template_bc, "combined_vc")
r_pet_dir_vc  <- create_raster_from_table(pet_annual, template_bc, "direction_vc")
r_pet_mag_vc  <- create_raster_from_table(pet_annual, template_bc, "sl_vc")

r_pet_comb_tfpw <- create_raster_from_table(pet_annual, template_bc, "combined_tfpw")
r_pet_dir_tfpw  <- create_raster_from_table(pet_annual, template_bc, "direction_tfpw")
r_pet_mag_tfpw  <- create_raster_from_table(pet_annual, template_bc, "sl_tfpw")

r_precip_comb_vc <- create_raster_from_table(precip_annual, template_bc, "combined_vc")
r_precip_tau_vc   <- create_raster_from_table(precip_annual, template_bc, "tau_vc")
r_precip_mag_vc   <- create_raster_from_table(precip_annual, template_bc, "sl_vc")

r_precip_comb_tfpw <- create_raster_from_table(precip_annual, template_bc, "combined_tfpw")
r_precip_tau_tfpw   <- create_raster_from_table(precip_annual, template_bc, "tau_tfpw")
r_precip_mag_tfpw   <- create_raster_from_table(precip_annual, template_bc, "sl_tfpw")

####################################################################################
# NOTE: TEMPORAL CLUSTERING (RUNS TEST)
# This analysis is not computed in the Pr/PET trend script because precipitation
# and PET are not drought indices - they don't have a meaningful threshold to
# define "dry" vs "wet" periods for runs test analysis.
####################################################################################

# create_clustering_maps <- function() {
#   # Would require: filtered_runs, clustering, p_value_runs columns
#   # These are not computed for Pr/PET (appropriately so)
# }

####################################################################################
# NEW: REGIME SHIFT DETECTION (CHANGEPOINT ANALYSIS)
####################################################################################

create_regime_shift_maps <- function() {
  cat("📊 Checking for regime shift (changepoint) results...\n")
  
  # Check if changepoint columns exist
  has_changepoint <- "first_changepoint_year" %in% names(all_results)
  
  if (!has_changepoint) {
    cat("⚠️  Changepoint columns not found. Skipping regime shift maps.\n")
    return(NULL)
  }
  
  cat("✅ Changepoint columns found! Creating regime shift maps...\n")
  
  pdf_path <- file.path(out_dir, "regime_shifts_changepoints.pdf")
  pdf(pdf_path, width = 12, height = 8, family = "Helvetica")
  
  par(mfrow = c(2, 2), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  for (var in c("PET", "Precipitation")) {
    var_data <- all_results[variable == var & period == "annual" & 
                              !is_basin_average & changepoint_detected == TRUE]
    
    if (nrow(var_data) == 0) {
      plot.new()
      text(0.5, 0.5, paste("No regime shifts detected for", var), cex = 1.2)
      next
    }
    
    cat(sprintf("  Found %d pixels with regime shifts for %s\n", nrow(var_data), var))
    
    # Map 1: Timing of first changepoint (by decade)
    var_data[, shift_decade := cut(first_changepoint_year,
                                   breaks = seq(1980, 2060, by = 10),
                                   labels = paste0(seq(1980, 2050, 10), "s"),
                                   include.lowest = TRUE)]
    
    # Convert to numeric for rasterization
    var_data[, shift_decade_num := as.numeric(cut(first_changepoint_year,
                                                  breaks = seq(1980, 2060, by = 10),
                                                  include.lowest = TRUE))]
    
    r_decade <- create_raster_from_table(var_data, template_bc, "shift_decade_num")
    
    plot_raster_panel(r_decade,
                      paste0(var, ": Decade of Regime Shift"),
                      c(1, 8),
                      hcl.colors(8, "Spectral", rev = TRUE),
                      legend_title = "")
    
    # Add legend manually
    decade_labels <- c("1980s", "1990s", "2000s", "2010s", "2020s", "2030s", "2040s", "2050s")
    legend("bottomright",
           legend = decade_labels[1:max(var_data$shift_decade_num, na.rm=TRUE)],
           fill = hcl.colors(8, "Spectral", rev = TRUE)[1:max(var_data$shift_decade_num, na.rm=TRUE)],
           bty = "n", cex = 0.8)
    
    # Map 2: Number of changepoints
    r_n_cpts <- create_raster_from_table(var_data, template_bc, "n_changepoints")
    
    plot_raster_panel(r_n_cpts,
                      paste0(var, ": Number of Changepoints"),
                      c(1, max(var_data$n_changepoints, na.rm=TRUE)),
                      hcl.colors(101, "YlOrRd", rev = FALSE),
                      legend_title = "N")
  }
  
  mtext("Regime Shift Detection (PELT Changepoint Analysis, 1980-2055)",
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  
  dev.off()
  cat(sprintf("✅ Regime shift map: %s\n", basename(pdf_path)))
  return(pdf_path)
}

####################################################################################
# NEW: SPECTRAL ANALYSIS
####################################################################################

create_spectral_analysis_maps <- function() {
  cat("📊 Checking for spectral analysis results...\n")
  
  # Check if spectral analysis columns exist
  has_spectral <- all(c("n_spectral_peaks", "dominant_period") %in% names(all_results))
  
  if (!has_spectral) {
    cat("⚠️  Spectral analysis columns not found. Skipping spectral maps.\n")
    cat("    Available columns:", paste(names(all_results)[1:20], collapse=", "), "...\n")
    return(NULL)
  }
  
  cat("✅ Spectral analysis columns found! Creating maps...\n")
  
  pdf_path <- file.path(out_dir, "spectral_analysis.pdf")
  pdf(pdf_path, width = 12, height = 8, family = "Helvetica")
  
  par(mfrow = c(2, 2), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  for (var in c("PET", "Precipitation")) {
    var_data <- all_results[variable == var & period == "annual" & 
                              !is_basin_average & n_spectral_peaks > 0]
    
    if (nrow(var_data) == 0) {
      plot.new()
      text(0.5, 0.5, paste("No significant spectral peaks detected for", var), cex = 1.2)
      next
    }
    
    cat(sprintf("  Found %d pixels with spectral peaks for %s\n", nrow(var_data), var))
    
    # Map of number of peaks
    r_peaks <- create_raster_from_table(var_data, template_bc, "n_spectral_peaks")
    
    plot_raster_panel(r_peaks,
                      paste0(var, ": Number of Spectral Peaks"),
                      c(0, max(var_data$n_spectral_peaks, na.rm = TRUE)),
                      hcl.colors(101, "viridis"),
                      legend_title = "N peaks")
    
    # Map of dominant period
    r_period <- create_raster_from_table(var_data, template_bc, "dominant_period")
    
    plot_raster_panel(r_period,
                      paste0(var, ": Dominant Cycle Period"),
                      c(min(var_data$dominant_period, na.rm = TRUE),
                        max(var_data$dominant_period, na.rm = TRUE)),
                      hcl.colors(101, "plasma"),
                      legend_title = "Years")
  }
  
  mtext("Spectral Analysis: Detected Periodicities in Climate Variables",
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  
  dev.off()
  cat(sprintf("✅ Spectral analysis map: %s\n", basename(pdf_path)))
  return(pdf_path)
}

####################################################################################
# NEW: SPATIAL AUTOCORRELATION ANALYSIS
####################################################################################

create_spatial_pattern_analysis <- function() {
  cat("📊 Creating spatial pattern analysis...\n")
  
  pdf_path <- file.path(out_dir, "spatial_patterns.pdf")
  pdf(pdf_path, width = 10, height = 6.67, family = "Helvetica")
  
  par(mfrow = c(1, 2), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  for (var in c("PET", "Precipitation")) {
    var_data <- all_results[variable == var & period == "annual" & 
                              !is_basin_average & !filtered_vc]
    
    if (nrow(var_data) == 0) next
    
    # Calculate spatial coefficient of variation for trends
    var_data[, spatial_cv := abs(sl_vc) / (abs(mean(sl_vc, na.rm = TRUE)) + 0.001)]
    
    r_cv <- create_raster_from_table(var_data, template_bc, "spatial_cv")
    
    plot_raster_panel(r_cv,
                      paste0(var, ": Spatial Variability"),
                      c(0, 5),
                      hcl.colors(101, "YlOrRd", rev = FALSE),
                      legend_title = "CV")
  }
  
  mtext("Spatial Variability of Trends (Annual Data, 1980-2055)",
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  
  dev.off()
  cat(sprintf("✅ Spatial pattern map: %s\n", basename(pdf_path)))
  return(pdf_path)
}

####################################################################################
# NEW: METHOD COMPARISON (VC vs TFPW) – FIXED
####################################################################################

create_method_comparison <- function() {
  cat("📊 Creating method comparison plots...\n")
  
  # Prepare data
  comp_data <- all_results[period == "annual" & !is_basin_average & 
                             !filtered_vc & !filtered_tfpw]
  comp_data[, vc_sig := p_value_vc < 0.05]
  comp_data[, tfpw_sig := p_value_tfpw < 0.05]
  
  # Determine global tau limits for equal scaling across facets
  tau_min <- min(comp_data$tau_vc, comp_data$tau_tfpw, na.rm = TRUE)
  tau_max <- max(comp_data$tau_vc, comp_data$tau_tfpw, na.rm = TRUE)
  tau_range <- c(tau_min, tau_max)
  
  # Panel A: Tau scatter with fixed scales and equal aspect ratio
  p1 <- ggplot(comp_data, aes(x = tau_vc, y = tau_tfpw, color = variable)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(alpha = 0.5, size = 1.5) +
    facet_wrap(~variable) +                          # fixed scales (default)
    scale_color_manual(values = c("PET" = "#d73027", "Precipitation" = "#4575b4")) +
    coord_equal(xlim = tau_range, ylim = tau_range) + # same limits & aspect
    theme_bw(base_size = 11) +
    theme(legend.position = "none") +
    labs(title = "A) Kendall's Tau: VC vs TFPW",
         x = "VC Method (τ)", y = "TFPW Method (τ)")
  
  # Panel B: Agreement summary
  agreement <- comp_data[, .(
    Both_Sig = sum(vc_sig & tfpw_sig),
    Only_VC = sum(vc_sig & !tfpw_sig),
    Only_TFPW = sum(!vc_sig & tfpw_sig),
    Neither = sum(!vc_sig & !tfpw_sig)
  ), by = variable]
  
  agreement_long <- melt(agreement, id.vars = "variable")
  
  p2 <- ggplot(agreement_long, aes(x = variable, y = value, fill = variable)) +
    geom_col(position = "stack") +
    facet_wrap(~variable, scales = "free_x") +
    geom_text(aes(label = sprintf("%.0f%%", value/sum(value)*100), group = variable),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c("PET" = "#d73027", "Precipitation" = "#4575b4")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = "B) Method Agreement on Significance",
         x = NULL, y = "Number of Pixels")
  
  # Panel C: Percentage significant
  pct_sig <- comp_data[, .(
    VC_Pct = sum(vc_sig) / .N * 100,
    TFPW_Pct = sum(tfpw_sig) / .N * 100
  ), by = variable]
  
  pct_long <- melt(pct_sig, id.vars = "variable")
  
  p3 <- ggplot(pct_long, aes(x = variable, y = value, fill = variable)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = sprintf("%.1f%%", value)), 
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("PET" = "#d73027", "Precipitation" = "#4575b4")) +
    facet_wrap(~variable, scales = "free_x") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(title = "C) % Pixels with Significant Trends",
         x = NULL, y = "Percentage (%)")
  
  # Combine
  combined <- p1 / (p2 | p3) +
    plot_annotation(
      title = "Mann-Kendall Method Comparison (VC vs TFPW)",
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
    )
  
  ggsave(file.path(out_dir, "method_comparison_vc_vs_tfpw.png"),
         combined, width = 12, height = 10, dpi = 300)
  
  cat("✅ Method comparison plot saved\n")
}

####################################################################################
# NEW: COMPREHENSIVE STATISTICS TABLE – FIXED
####################################################################################

create_statistics_table <- function() {
  cat("📊 Creating statistics summary table...\n")
  
  # Check for autocorrelation columns
  has_rho_vc  <- "rho1_vc" %in% names(all_results)
  has_rho_tfpw <- "rho1_tfpw" %in% names(all_results)
  
  stats_summary <- all_results[period == "annual" & !is_basin_average, .(
    N_Pixels = .N,
    N_Valid = sum(!filtered_vc),
    Pct_Sig_VC = sprintf("%.1f%%", sum(p_value_vc < 0.05 & !filtered_vc) / sum(!filtered_vc) * 100),
    Pct_Sig_TFPW = sprintf("%.1f%%", sum(p_value_tfpw < 0.05 & !filtered_tfpw) / sum(!filtered_tfpw) * 100),
    Median_Tau_VC = sprintf("%.3f", median(tau_vc[!filtered_vc], na.rm = TRUE)),
    Median_Tau_TFPW = sprintf("%.3f", median(tau_tfpw[!filtered_tfpw], na.rm = TRUE)),
    Median_Slope_VC = sprintf("%.3f", median(sl_vc[!filtered_vc], na.rm = TRUE)),
    Median_Slope_TFPW = sprintf("%.3f", median(sl_tfpw[!filtered_tfpw], na.rm = TRUE)),
    Mean_Autocorr_VC = if (has_rho_vc) sprintf("%.3f", mean(rho1_vc[!filtered_vc], na.rm = TRUE)) else NA_character_,
    Mean_Autocorr_TFPW = if (has_rho_tfpw) sprintf("%.3f", mean(rho1_tfpw[!filtered_tfpw], na.rm = TRUE)) else NA_character_
  ), by = .(variable, period)]
  
  fwrite(stats_summary, file.path(out_dir, "summary_statistics.csv"))
  cat("✅ Statistics table saved\n")
  
  return(stats_summary)
}

####################################################################################
# ORIGINAL FIGURE GENERATORS (PRESERVED)
####################################################################################

create_basin_timeseries_plot <- function() {
  cat("📈 Generating basin-averaged time series plot...\n")
  pdf_path <- file.path(out_dir, "basin_average_timeseries.pdf")
  pdf(pdf_path, width = 12, height = 7, family = "Helvetica")
  
  basin_annual_results <- all_results[is_basin_average == TRUE & period == "annual"]
  pr_res <- basin_annual_results[variable == "Precipitation"]
  pet_res <- basin_annual_results[variable == "PET"]
  
  pr_rolling <- rollmean(basin_avg_monthly$precip_mm_month, k = 12, fill = NA, align = "center")
  pet_rolling <- rollmean(basin_avg_monthly$pet_mm_month, k = 12, fill = NA, align = "center")
  
  par(mar = c(5, 5, 4, 5) + 0.1)
  
  plot(basin_avg_monthly$date, basin_avg_monthly$precip_mm_month, 
       type = "l", col = rgb(0.2, 0.4, 0.8, 0.3), lwd = 0.8,
       xlab = "Year", ylab = "Precipitation (mm/month)", 
       main = "Nechako Basin Climate Time Series (1980-2055)\nSpatially Averaged Monthly Values",
       ylim = c(0, max(basin_avg_monthly$precip_mm_month, na.rm = TRUE) * 1.1),
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  lines(basin_avg_monthly$date, pr_rolling, col = "blue", lwd = 2.5)
  
  par(new = TRUE)
  plot(basin_avg_monthly$date, basin_avg_monthly$pet_mm_month,
       type = "l", col = rgb(0.8, 0.2, 0.2, 0.3), lwd = 0.8,
       axes = FALSE, xlab = "", ylab = "", 
       ylim = c(0, max(basin_avg_monthly$pet_mm_month, na.rm = TRUE) * 1.1))
  lines(basin_avg_monthly$date, pet_rolling, col = "red", lwd = 2.5)
  
  axis(side = 4, at = pretty(range(basin_avg_monthly$pet_mm_month, na.rm = TRUE)),
       col = "red", col.axis = "red", cex.axis = 1.1)
  mtext("PET (mm/month)", side = 4, line = 3, col = "red", cex = 1.2)
  
  legend("topright", 
         legend = c("Precipitation (monthly)", "Pr 12-mo mean",
                    "PET (monthly)", "PET 12-mo mean"),
         col = c(rgb(0.2, 0.4, 0.8, 0.3), "blue",
                 rgb(0.8, 0.2, 0.2, 0.3), "red"),
         lwd = c(0.8, 2.5, 0.8, 2.5),
         lty = c(1, 1, 1, 1),
         bty = "n", cex = 1.1)
  
  abline(v = as.Date(paste0(c(2000, 2020), "-01-01")), col = "gray40", lty = 2, lwd = 1)
  text(x = as.Date("2010-01-01"), y = par("usr")[4] * 0.95, "2000-2020", col = "gray40", cex = 0.9)
  
  dev.off()
  cat(sprintf("✅ Basin time series plot: %s\n", basename(pdf_path)))
  invisible(pdf_path)
}

create_publication_figure <- function(var_name, method = "vc") {
  pdf_path <- file.path(out_dir, sprintf("%s_%s_publication_map.pdf", tolower(var_name), method))
  pdf(pdf_path, width = 10, height = 6.67, family = "Helvetica")
  
  par(mfrow = c(1, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  if (var_name == "PET") {
    r_comb <- if (method == "vc") r_pet_comb_vc else r_pet_comb_tfpw
    r_dir  <- if (method == "vc") r_pet_dir_vc  else r_pet_dir_tfpw
    r_mag  <- if (method == "vc") r_pet_mag_vc  else r_pet_mag_tfpw
    
    plot_raster_panel(r_comb, "Combined Trend Strength", c(0, 0.25), 
                      hcl.colors(101, "YlOrRd", rev = FALSE),
                      legend_title = "|τ|")
    plot_raster_panel(r_mag, "Sen's Slope (mm/yr)", c(0, 1.5), 
                      hcl.colors(101, "YlOrRd", rev = FALSE),
                      legend_title = "mm/yr")
    plot_raster_panel(r_dir, "Direction (p < 0.05)", NULL, 
                      c("#4575b4", "#f0f0f0", "#d73027"), 
                      breaks = c(-1.5, -0.5, 0.5, 1.5), legend = FALSE, categorical = TRUE)
    legend("bottomright", legend = c("Decreasing", "Non-sig.", "Increasing"),
           fill = c("#4575b4", "#f0f0f0", "#d73027"), bty = "n", cex = 1.1)
    title_prefix <- "PET"
  } else {
    r_comb <- if (method == "vc") r_precip_comb_vc else r_precip_comb_tfpw
    r_tau  <- if (method == "vc") r_precip_tau_vc   else r_precip_tau_tfpw
    r_mag  <- if (method == "vc") r_precip_mag_vc   else r_precip_mag_tfpw
    
    plot_raster_panel(r_comb, "Combined Trend Strength", c(-0.15, 0.15), 
                      hcl.colors(101, "RdBu", rev = TRUE),
                      legend_title = "τ (sig.)")
    plot_raster_panel(r_mag, "Sen's Slope (mm/yr)", c(-2, 2), 
                      hcl.colors(101, "RdBu", rev = TRUE),
                      legend_title = "mm/yr")
    plot_raster_panel(r_tau, "Kendall's Tau", c(-0.3, 0.3), 
                      hcl.colors(101, "RdBu", rev = TRUE),
                      legend_title = "τ")
    title_prefix <- "Precipitation"
  }
  
  mtext(sprintf("%s Annual Trends (1980-2055) | %s Method", title_prefix,
                ifelse(method == "vc", "Variance-Corrected MK", "TFPW-Corrected MK")),
        outer = TRUE, cex = 1.1, font = 2, line = 0.5)
  dev.off()
  
  cat(sprintf("✅ %s %s map: %s\n", var_name, method, basename(pdf_path)))
  invisible(pdf_path)
}

create_comparison_figure <- function() {
  pdf_path <- file.path(out_dir, "pet_vs_precipitation_comparison.pdf")
  pdf(pdf_path, width = 12, height = 8, family = "Helvetica")
  
  par(mfrow = c(2, 2), mar = c(2, 2, 2.5, 1), oma = c(1, 1, 3, 1))
  
  plot_raster_panel(r_pet_comb_vc, "PET - Significant Trend Strength", c(0, 0.25), 
                    hcl.colors(101, "YlOrRd", rev = FALSE),
                    legend_title = "|τ|")
  plot_raster_panel(r_pet_mag_vc, "PET - Sen's Slope (mm/yr)", c(0, 1.5), 
                    hcl.colors(101, "YlOrRd", rev = FALSE),
                    legend_title = "mm/yr")
  plot_raster_panel(r_precip_comb_vc, "Precipitation - Significant Trend Strength", c(-0.15, 0.15), 
                    hcl.colors(101, "RdBu", rev = TRUE),
                    legend_title = "τ (sig.)")
  plot_raster_panel(r_precip_tau_vc, "Precipitation - Kendall's Tau", c(-0.3, 0.3), 
                    hcl.colors(101, "RdBu", rev = TRUE),
                    legend_title = "τ")
  
  mtext("Nechako Basin Climate Trends (1980-2055) | Variance-Corrected MK", 
        outer = TRUE, cex = 1.3, font = 2, line = 0.5)
  dev.off()
  
  cat(sprintf("✅ Comparison map: %s\n", basename(pdf_path)))
  invisible(pdf_path)
}

####################################################################################
# GENERATE ALL FIGURES
####################################################################################

cat("\n========================================\n")
cat("🎨 CREATING ENHANCED VISUALIZATIONS\n")
cat("========================================\n\n")

# Original figures
basin_ts_plot <- create_basin_timeseries_plot()
pet_vc_fig    <- create_publication_figure("PET", "vc")
pet_tfpw_fig  <- create_publication_figure("PET", "tfpw")
precip_vc_fig <- create_publication_figure("Precipitation", "vc")
precip_tfpw_fig <- create_publication_figure("Precipitation", "tfpw")
comp_fig      <- create_comparison_figure()

# NEW: Enhanced diagnostics
regime_fig       <- create_regime_shift_maps()  # Will check if columns exist
spectral_fig     <- create_spectral_analysis_maps()  # Will check if columns exist
spatial_patterns <- create_spatial_pattern_analysis()  
method_comp      <- create_method_comparison()
stats_table      <- create_statistics_table()

cat("\n========================================\n")
cat("✅ ALL ENHANCED VISUALIZATIONS CREATED\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n\n", file.path(getwd(), out_dir)))

cat("Generated files:\n")
cat("ORIGINAL OUTPUTS:\n")
cat("  • basin_average_timeseries.pdf\n")
cat("  • pet_vc_publication_map.pdf\n")
cat("  • pet_tfpw_publication_map.pdf\n")
cat("  • precipitation_vc_publication_map.pdf\n")
cat("  • precipitation_tfpw_publication_map.pdf\n")
cat("  • pet_vs_precipitation_comparison.pdf\n")
cat("\nNEW ENHANCED OUTPUTS:\n")
if (!is.null(regime_fig)) {
  cat("  • regime_shifts_changepoints.pdf         ← PELT changepoint detection\n")
}
if (!is.null(spectral_fig)) {
  cat("  • spectral_analysis.pdf                  ← Periodic patterns in Pr/PET\n")
}
cat("  • spatial_patterns.pdf                   ← Spatial variability analysis\n")
cat("  • method_comparison_vc_vs_tfpw.png       ← Method comparison (VC vs TFPW)\n")
cat("  • summary_statistics.csv                 ← Comprehensive statistics table\n")
cat("\nNOTE: The script automatically checks which diagnostics are available.\n")
cat("      Regime shifts and spectral analysis should be included if the\n")
cat("      trend analysis script computed them (which it does by default).\n")