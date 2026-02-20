####################################################################################
# PUBLICATION-QUALITY TREND MAPS FOR PR/PET
# ENHANCED WITH:
# - Temporal clustering analysis (runs test)
# - Regime shift detection
# - Spectral analysis
# - Method comparison plots
# - Comprehensive statistical summaries
#
# FIXES APPLIED (see ## FIX comments throughout):
# 1. Direction maps: NA/filtered pixels now shown as "Non-sig." (0), not blank holes
# 2. Regime shift maps: pixels with no detected shift shown explicitly as "No shift"
#    instead of being silently dropped, causing spatial gaps
# 3. Precipitation "Combined Trend Strength": renamed and replaced with raw Kendall's τ
#    because the original masked all non-significant pixels to 0 (all of them), 
#    producing an uninformative flat map
# 4. Comparison figure: improved titles and captions for Precipitation panel
# 5. Method comparison Panel B: fixed percentage labels (were computed across all 
#    variables combined, not per-variable)
# 6. Method comparison Panel C: VC vs TFPW bars now colored/labeled by method, not 
#    by variable; added annotation explaining 0% for precipitation is a valid result
# 7. All figures: improved legend titles and subtitles for clarity
# 8. TIME RANGE: Corrected ALL instances from 1980-2055 to 1950-2025
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
all_results    <- readRDS(file.path("trend_analysis_pr_pet", "all_results.rds"))
basin_boundary <- readRDS(file.path("trend_analysis_pr_pet", "basin_boundary.rds"))
metadata       <- readRDS(file.path("trend_analysis_pr_pet", "analysis_metadata.rds"))

basin_avg_monthly <- metadata$basin_avg_monthly
basin_avg_annual  <- metadata$basin_avg_annual

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

analysis_cols <- list(
  Basic       = c("tau_vc", "p_value_vc", "sl_vc", "tau_tfpw", "p_value_tfpw"),
  Autocorr    = c("rho1_vc", "rho1_tfpw"),
  Spectral    = c("n_spectral_peaks", "dominant_period", "spectral_confidence"),
  Changepoint = c("changepoint_detected", "first_changepoint_year", "n_changepoints"),
  RunsTest    = c("filtered_runs", "clustering", "p_value_runs")
)

for (analysis_name in names(analysis_cols)) {
  cols_to_check <- analysis_cols[[analysis_name]]
  available <- cols_to_check[cols_to_check %in% all_cols]
  missing   <- cols_to_check[!cols_to_check %in% all_cols]
  if (length(available) > 0)
    cat(sprintf("   ✅ %s: %s\n", analysis_name, paste(available, collapse = ", ")))
  if (length(missing) > 0)
    cat(sprintf("   ⚠️  %s (missing): %s\n", analysis_name, paste(missing, collapse = ", ")))
}
cat("\n")

# ===== CREATE NATIVE-RESOLUTION TEMPLATE =====
pet_annual_sample <- all_results[variable == "PET" & period == "annual" & !filtered_vc][1:100]
dx <- median(diff(sort(unique(pet_annual_sample$x))), na.rm = TRUE)
dy <- median(diff(sort(unique(pet_annual_sample$y))), na.rm = TRUE)
resolution_bc <- min(dx, dy, na.rm = TRUE)

basin_ext   <- ext(basin_boundary)
template_bc <- rast(ext        = basin_ext,
                    resolution = resolution_bc,
                    crs        = "EPSG:3005",
                    names      = "trend")

cat(sprintf("✓ Native-resolution template (%.1f m, %d x %d cells)\n",
            resolution_bc, ncol(template_bc), nrow(template_bc)))

# ===== HELPER FUNCTIONS =====
create_raster_from_table <- function(results_dt, template, value_col) {
  pts <- vect(results_dt, geom = c("x", "y"), crs = crs(template))
  r   <- rasterize(pts, template, field = value_col, touches = TRUE)
  mask(r, vect(basin_boundary))
}

smooth_raster <- function(r) {
  if (all(is.na(values(r)))) return(r)
  focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
}

## FIX 1 helper: fill NA cells in a raster using nearest-neighbour
fill_na_nearest <- function(r) {
  if (all(is.na(values(r)))) return(r)
  
  modal_fun <- function(x, na.rm = TRUE) {
    if (na.rm) x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    as.numeric(names(sort(table(x), decreasing = TRUE))[1])
  }
  
  for (i in 1:3) {
    r <- focal(r, w = matrix(1, 3, 3), fun = modal_fun, na.rm = TRUE, na.policy = "only")
  }
  mask(r, vect(basin_boundary))
}

plot_raster_panel <- function(r, main, zlim = NULL, col,
                              breaks = NULL, legend = TRUE, categorical = FALSE,
                              legend_title = "", subtitle = NULL) {
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
  
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = -0.3, cex = 0.72, col = "grey30", font = 3)
  }
}

# ===== PREPARE ANNUAL RESULTS =====
cat("⚙️  Preparing trend metrics...\n")

pet_annual    <- all_results[variable == "PET"           & period == "annual" &
                               !filtered_vc & !is_basin_average]
precip_annual <- all_results[variable == "Precipitation" & period == "annual" &
                               !filtered_vc & !is_basin_average]

# ── PET metrics ──────────────────────────────────────────────────────────────
pet_annual[, combined_vc   := fifelse(p_value_vc   < 0.05, abs(tau_vc),   0)]
pet_annual[, combined_tfpw := fifelse(p_value_tfpw < 0.05, abs(tau_tfpw), 0)]

# ── Precipitation: use RAW Kendall's τ ───────────────
precip_annual[, combined_vc   := tau_vc]
precip_annual[, combined_tfpw := tau_tfpw]

# ── Direction rasters ──────────────────────────────────────────────────────────
pet_all <- all_results[variable == "PET" & period == "annual" & !is_basin_average]
pet_all[, direction_vc   := fifelse(!filtered_vc   & p_value_vc   < 0.05 & tau_vc   < 0, -1L,
                                    fifelse(!filtered_vc   & p_value_vc   < 0.05 & tau_vc   > 0,  1L, 0L))]
pet_all[, direction_tfpw := fifelse(!filtered_tfpw & p_value_tfpw < 0.05 & tau_tfpw < 0, -1L,
                                    fifelse(!filtered_tfpw & p_value_tfpw < 0.05 & tau_tfpw > 0,  1L, 0L))]

# ===== GENERATE RASTERS =====
cat("🖼️  Generating rasters...\n")

r_pet_comb_vc   <- create_raster_from_table(pet_annual, template_bc, "combined_vc")
r_pet_mag_vc    <- create_raster_from_table(pet_annual, template_bc, "sl_vc")
r_pet_comb_tfpw <- create_raster_from_table(pet_annual, template_bc, "combined_tfpw")
r_pet_mag_tfpw  <- create_raster_from_table(pet_annual, template_bc, "sl_tfpw")

r_pet_dir_vc_raw   <- create_raster_from_table(pet_all, template_bc, "direction_vc")
r_pet_dir_tfpw_raw <- create_raster_from_table(pet_all, template_bc, "direction_tfpw")
r_pet_dir_vc       <- fill_na_nearest(r_pet_dir_vc_raw)
r_pet_dir_tfpw     <- fill_na_nearest(r_pet_dir_tfpw_raw)

r_precip_comb_vc    <- create_raster_from_table(precip_annual, template_bc, "combined_vc")
r_precip_tau_vc     <- create_raster_from_table(precip_annual, template_bc, "tau_vc")
r_precip_mag_vc     <- create_raster_from_table(precip_annual, template_bc, "sl_vc")
r_precip_comb_tfpw  <- create_raster_from_table(precip_annual, template_bc, "combined_tfpw")
r_precip_tau_tfpw   <- create_raster_from_table(precip_annual, template_bc, "tau_tfpw")
r_precip_mag_tfpw   <- create_raster_from_table(precip_annual, template_bc, "sl_tfpw")

####################################################################################
# REGIME SHIFT DETECTION (CHANGEPOINT ANALYSIS)
####################################################################################
create_regime_shift_maps <- function() {
  cat("📊 Checking for regime shift (changepoint) results...\n")
  has_changepoint <- "first_changepoint_year" %in% names(all_results)
  if (!has_changepoint) {
    cat("⚠️  Changepoint columns not found. Skipping regime shift maps.\n")
    return(NULL)
  }
  cat("✅ Changepoint columns found! Creating regime shift maps...\n")
  pdf_path <- file.path(out_dir, "regime_shifts_changepoints.pdf")
  pdf(pdf_path, width = 14, height = 9, family = "Helvetica")
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 1), oma = c(1, 1, 3, 1))
  
  for (var in c("PET", "Precipitation")) {
    var_all <- all_results[variable == var & period == "annual" & !is_basin_average]
    
    if (nrow(var_all) == 0) {
      plot.new(); text(0.5, 0.5, paste("No data for", var), cex = 1.2); next
    }
    
    n_detected <- sum(var_all$changepoint_detected == TRUE, na.rm = TRUE)
    cat(sprintf("  %s: %d / %d pixels with detected regime shift\n",
                var, n_detected, nrow(var_all)))
    
    # ── Map 1: Timing of first changepoint (by decade) ──────────────────────
    ## FIX 8a: Correct time range 1950-2025 (not 1980-2055)
    var_all[, shift_decade_num := fifelse(
      changepoint_detected == TRUE,
      as.numeric(cut(first_changepoint_year,
                     breaks = seq(1950, 2030, by = 10),  # ← FIXED
                     include.lowest = TRUE)),
      0L
    )]
    
    r_decade <- create_raster_from_table(var_all, template_bc, "shift_decade_num")
    
    # ## FIX: Use custom modal function (NOT direct 'modal')
    modal_fun <- function(x, na.rm = TRUE) {
      if (na.rm) x <- x[!is.na(x)]
      if (length(x) == 0) return(NA)
      as.numeric(names(sort(table(x), decreasing = TRUE))[1])
    }
    
    if (!all(is.na(values(r_decade)))) {
      for (i in 1:3) {
        r_decade <- focal(r_decade, w = matrix(1, 3, 3), fun = modal_fun, 
                          na.rm = TRUE, na.policy = "only")
      }
      r_decade <- mask(r_decade, vect(basin_boundary))
    }
    
    max_decade <- max(var_all$shift_decade_num, na.rm = TRUE)
    n_colors <- max(max_decade, 1L)
    
    decade_breaks <- seq(-0.5, max_decade + 0.5, by = 1)
    decade_colors <- c("#CCCCCC", hcl.colors(n_colors, "Spectral", rev = TRUE))
    
    plot(r_decade,
         main = paste0(var, ": Decade of First Regime Shift"),
         cex.main = 1.1, font.main = 2,
         col = decade_colors,
         breaks = decade_breaks,
         axes = FALSE, box = FALSE, legend = FALSE)
    plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
    mtext("0 = no shift detected", side = 3, line = -0.5, cex = 0.8, col = "grey40", font = 3)
    
    ## FIX 8b: Decade labels match 1950-2025 range
    decade_labels <- c("No shift", paste0(seq(1950, by = 10, length.out = n_colors), "–", 
                                          seq(1959, by = 10, length.out = n_colors)))
    legend("bottomright",
           legend = decade_labels,
           fill = decade_colors,
           border = NA,
           bty = "n",
           cex = 0.75,
           title = "Shift Decade",
           title.cex = 0.85,
           ncol = 2)
    
    # ── Map 2: Number of changepoints ──────────────────────────────────────
    var_all[, n_cp_plot := fifelse(changepoint_detected == TRUE,
                                   as.numeric(n_changepoints), 0)]
    
    r_n_cpts <- create_raster_from_table(var_all, template_bc, "n_cp_plot")
    
    raster_vals <- values(r_n_cpts)
    raster_vals <- raster_vals[!is.na(raster_vals)]
    existing_vals <- sort(unique(raster_vals))
    n_existing <- length(existing_vals)
    
    cat(sprintf("  %s changepoint counts: min=%d, max=%d, unique=%s (n=%d)\n",
                var, min(raster_vals, na.rm = TRUE), max(raster_vals, na.rm = TRUE),
                paste(existing_vals, collapse = ", "), n_existing))
    
    ## FIX: Use distinguishable colors for only existing values
    distinct_palette <- c(
      "#E69F00", "#56B4E9", "#009E73", "#D55E00",
      "#CC79A7", "#0072B2", "#F0E442", "#000000",
      "#FF6600", "#00CC66", "#9933FF", "#FF3399",
      "#33CCCC", "#FF9933", "#6633CC", "#99CC00"
    )
    
    if (n_existing <= length(distinct_palette)) {
      n_cp_colors <- distinct_palette[1:n_existing]
    } else {
      n_cp_colors <- hcl.colors(n_existing, "Dark 2")
    }
    
    n_cp_breaks <- existing_vals - 0.5
    n_cp_breaks <- c(n_cp_breaks, max(existing_vals) + 0.5)
    
    plot(r_n_cpts,
         main = paste0(var, ": Number of Detected Changepoints"),
         cex.main = 1.1, font.main = 2,
         col = n_cp_colors,
         breaks = n_cp_breaks,
         axes = FALSE, box = FALSE,
         legend = FALSE)
    plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
    mtext("Each colour = exact changepoint count", side = 3, line = -0.5,
          cex = 0.8, col = "grey40", font = 3)
    
    legend("bottomright",
           legend = as.character(existing_vals),
           fill = n_cp_colors,
           border = NA,
           bty = "n",
           cex = 1.0,
           title = "N Changepoints",
           title.cex = 1.1,
           title.font = 2)
  }
  
  ## FIX 8c: Outer title matches actual data range (1950-2025)
  mtext("Regime Shift Detection (PELT Changepoint Analysis, 1950–2025)",
        outer = TRUE, cex = 1.3, font = 2, line = 0.5)
  dev.off()
  cat(sprintf("✅ Regime shift map: %s\n", basename(pdf_path)))
  return(pdf_path)
}
####################################################################################
# SPECTRAL ANALYSIS — TIME RANGE FIXED
####################################################################################
create_spectral_analysis_maps <- function() {
  cat("📊 Checking for spectral analysis results...\n")
  
  has_spectral <- all(c("n_spectral_peaks", "dominant_period") %in% names(all_results))
  if (!has_spectral) {
    cat("⚠️  Spectral analysis columns not found. Skipping spectral maps.\n")
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
    
    r_peaks <- create_raster_from_table(var_data, template_bc, "n_spectral_peaks")
    plot_raster_panel(r_peaks,
                      paste0(var, ": Number of Spectral Peaks"),
                      c(0, max(var_data$n_spectral_peaks, na.rm = TRUE)),
                      hcl.colors(101, "viridis"),
                      legend_title = "N peaks")
    
    r_period <- create_raster_from_table(var_data, template_bc, "dominant_period")
    plot_raster_panel(r_period,
                      paste0(var, ": Dominant Cycle Period"),
                      c(min(var_data$dominant_period, na.rm = TRUE),
                        max(var_data$dominant_period, na.rm = TRUE)),
                      hcl.colors(101, "plasma"),
                      legend_title = "Years")
  }
  
  ## FIX 8d: Add "Annual Data" and correct year range to outer title
  mtext("Spectral Analysis: Detected Periodicities in Climate Variables\n(Annual Data, 1950–2025)",
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  dev.off()
  cat(sprintf("✅ Spectral analysis map: %s\n", basename(pdf_path)))
  return(pdf_path)
}

####################################################################################
# SPATIAL PATTERN ANALYSIS — TIME RANGE FIXED
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
    
    var_data[, spatial_cv := abs(sl_vc) / (abs(mean(sl_vc, na.rm = TRUE)) + 0.001)]
    r_cv <- create_raster_from_table(var_data, template_bc, "spatial_cv")
    
    plot_raster_panel(r_cv,
                      paste0(var, ": Spatial Variability"),
                      c(0, 5),
                      hcl.colors(101, "YlOrRd", rev = FALSE),
                      legend_title = "CV")
  }
  
  ## FIX 8e: Outer title matches actual data range (1950-2025)
  mtext("Spatial Variability of Trends (Annual Data, 1950–2025)",
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  dev.off()
  cat(sprintf("✅ Spatial pattern map: %s\n", basename(pdf_path)))
  return(pdf_path)
}

####################################################################################
# METHOD COMPARISON (VC vs TFPW)
####################################################################################
create_method_comparison <- function() {
  cat("📊 Creating method comparison plots...\n")
  
  comp_data <- all_results[period == "annual" & !is_basin_average &
                             !filtered_vc & !filtered_tfpw]
  comp_data[, vc_sig   := p_value_vc   < 0.05]
  comp_data[, tfpw_sig := p_value_tfpw < 0.05]
  
  tau_min   <- min(comp_data$tau_vc, comp_data$tau_tfpw, na.rm = TRUE)
  tau_max   <- max(comp_data$tau_vc, comp_data$tau_tfpw, na.rm = TRUE)
  tau_range <- c(tau_min, tau_max)
  
  p1 <- ggplot(comp_data, aes(x = tau_vc, y = tau_tfpw, color = variable)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(alpha = 0.5, size = 1.5) +
    facet_wrap(~variable) +
    scale_color_manual(values = c("PET" = "#d73027", "Precipitation" = "#4575b4")) +
    coord_equal(xlim = tau_range, ylim = tau_range) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none") +
    labs(title = "A) Kendall's Tau: VC vs TFPW",
         x = "VC Method (τ)", y = "TFPW Method (τ)")
  
  agreement <- comp_data[, .(
    Both_Sig  = sum( vc_sig &  tfpw_sig),
    Only_VC   = sum( vc_sig & !tfpw_sig),
    Only_TFPW = sum(!vc_sig &  tfpw_sig),
    Neither   = sum(!vc_sig & !tfpw_sig)
  ), by = variable]
  
  agreement_long <- melt(agreement, id.vars = "variable",
                         variable.name = "Agreement", value.name = "N_pixels")
  agreement_long[, pct := N_pixels / sum(N_pixels) * 100, by = variable]
  
  agreement_colors <- c(
    "Both_Sig"  = "#2ca02c",
    "Only_VC"   = "#e08214",
    "Only_TFPW" = "#542788",
    "Neither"   = "#999999"
  )
  
  agreement_labels <- c(
    "Both_Sig"  = "Both Sig.",
    "Only_VC"   = "Only VC",
    "Only_TFPW" = "Only TFPW",
    "Neither"   = "Neither"
  )
  
  p2 <- ggplot(agreement_long,
               aes(x = variable, y = N_pixels, fill = Agreement)) +
    geom_col(position = "stack", width = 0.6) +
    facet_wrap(~variable, scales = "free_x") +
    geom_text(aes(label = ifelse(pct > 5, sprintf("%.0f%%", pct), "")),
              position = position_stack(vjust = 0.5),
              size = 3.5, fontface = "bold", color = "white") +
    scale_fill_manual(values = agreement_colors,
                      labels = agreement_labels,
                      name = "Agreement") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank()) +
    labs(title = "B) Method Agreement on Significance",
         subtitle = "Per-variable percentage of pixels",
         x = NULL, y = "Number of Pixels")
  
  pct_sig <- comp_data[, .(
    VC   = sum(vc_sig)   / .N * 100,
    TFPW = sum(tfpw_sig) / .N * 100
  ), by = variable]
  
  pct_long <- melt(pct_sig, id.vars = "variable",
                   variable.name = "Method", value.name = "Pct")
  
  precip_note <- "Precipitation 0%: no pixel reaches p < 0.05 after\nautocorrelation correction — high interannual variability\nmasks the weak directional signal (|τ| ≤ 0.10)"
  
  p3 <- ggplot(pct_long, aes(x = Method, y = Pct, fill = Method)) +
    geom_col(position = "dodge", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", Pct)),
              position = position_dodge(width = 0.9),
              vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = c("VC" = "#e08214", "TFPW" = "#542788"),
                      name = "Method") +
    facet_wrap(~variable, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom") +
    labs(title    = "C) % Pixels with Significant Trends (p < 0.05)",
         subtitle = precip_note,
         x = NULL, y = "Percentage (%)")
  
  combined <- p1 / (p2 | p3) +
    plot_annotation(
      title = "Mann-Kendall Method Comparison (VC vs TFPW) | 1950–2025",
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
    )
  
  ggsave(file.path(out_dir, "method_comparison_vc_vs_tfpw.png"),
         combined, width = 12, height = 10, dpi = 300)
  cat("✅ Method comparison plot saved\n")
}

####################################################################################
# COMPREHENSIVE STATISTICS TABLE
####################################################################################
create_statistics_table <- function() {
  cat("📊 Creating statistics summary table...\n")
  
  has_rho_vc   <- "rho1_vc"   %in% names(all_results)
  has_rho_tfpw <- "rho1_tfpw" %in% names(all_results)
  
  stats_summary <- all_results[period == "annual" & !is_basin_average, .(
    N_Pixels          = .N,
    N_Valid           = sum(!filtered_vc),
    Pct_Sig_VC        = sprintf("%.1f%%", sum(p_value_vc   < 0.05 & !filtered_vc)   / sum(!filtered_vc)   * 100),
    Pct_Sig_TFPW      = sprintf("%.1f%%", sum(p_value_tfpw < 0.05 & !filtered_tfpw) / sum(!filtered_tfpw) * 100),
    Median_Tau_VC     = sprintf("%.3f", median(tau_vc  [!filtered_vc],   na.rm = TRUE)),
    Median_Tau_TFPW   = sprintf("%.3f", median(tau_tfpw[!filtered_tfpw], na.rm = TRUE)),
    Median_Slope_VC   = sprintf("%.3f", median(sl_vc  [!filtered_vc],   na.rm = TRUE)),
    Median_Slope_TFPW = sprintf("%.3f", median(sl_tfpw[!filtered_tfpw], na.rm = TRUE)),
    Mean_Autocorr_VC  = if (has_rho_vc)   sprintf("%.3f", mean(rho1_vc  [!filtered_vc],   na.rm = TRUE)) else NA_character_,
    Mean_Autocorr_TFPW= if (has_rho_tfpw) sprintf("%.3f", mean(rho1_tfpw[!filtered_tfpw], na.rm = TRUE)) else NA_character_
  ), by = .(variable, period)]
  
  fwrite(stats_summary, file.path(out_dir, "summary_statistics.csv"))
  cat("✅ Statistics table saved\n")
  return(stats_summary)
}

####################################################################################
# BASIN TIME SERIES PLOT — TIME RANGE FIXED
####################################################################################
create_basin_timeseries_plot <- function() {
  cat("📈 Generating basin-averaged time series plot...\n")
  pdf_path <- file.path(out_dir, "basin_average_timeseries.pdf")
  pdf(pdf_path, width = 12, height = 7, family = "Helvetica")
  
  pr_rolling  <- rollmean(basin_avg_monthly$precip_mm_month, k = 12, fill = NA, align = "center")
  pet_rolling <- rollmean(basin_avg_monthly$pet_mm_month,    k = 12, fill = NA, align = "center")
  
  par(mar = c(5, 5, 4, 5) + 0.1)
  
  ## FIX 8f: Main title matches actual data range (1950-2025)
  plot(basin_avg_monthly$date, basin_avg_monthly$precip_mm_month,
       type = "l", col = rgb(0.2, 0.4, 0.8, 0.3), lwd = 0.8,
       xlab = "Year", ylab = "Precipitation (mm/month)",
       main = "Nechako Basin Climate Time Series (1950–2025)\nSpatially Averaged Monthly Values",
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
         bty = "n", cex = 1.1)
  
  ## FIX 8g: Reference lines for mid-period (1985, 2005 instead of 2000, 2020)
  abline(v = as.Date(paste0(c(1985, 2005), "-01-01")), col = "gray40", lty = 2, lwd = 1)
  text(x = as.Date("1995-01-01"), y = par("usr")[4] * 0.95, "1985–2005",
       col = "gray40", cex = 0.9)
  
  dev.off()
  cat(sprintf("✅ Basin time series plot: %s\n", basename(pdf_path)))
  invisible(pdf_path)
}

####################################################################################
# PUBLICATION MAPS (PET and Precipitation, both methods) — TIME RANGE FIXED
####################################################################################
create_publication_figure <- function(var_name, method = "vc") {
  pdf_path <- file.path(out_dir,
                        sprintf("%s_%s_publication_map.pdf", tolower(var_name), method))
  pdf(pdf_path, width = 12, height = 7, family = "Helvetica")
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 1), oma = c(1, 1, 3, 1))
  
  if (var_name == "PET") {
    r_comb <- if (method == "vc") r_pet_comb_vc   else r_pet_comb_tfpw
    r_dir  <- if (method == "vc") r_pet_dir_vc    else r_pet_dir_tfpw
    r_mag  <- if (method == "vc") r_pet_mag_vc    else r_pet_mag_tfpw
    
    plot_raster_panel(r_comb, "A) Significant Trend Strength", c(0, 0.25),
                      hcl.colors(101, "YlOrRd", rev = FALSE),
                      legend_title = "|τ|",
                      subtitle = "p < 0.05; 0 = non-significant")
    
    plot_raster_panel(r_mag, "B) Sen's Slope", c(0, 0.025),
                      hcl.colors(101, "YlOrRd", rev = FALSE),
                      legend_title = "mm/yr")
    
    plot_raster_panel(r_dir, "C) Trend Direction", NULL,
                      c("#4575b4", "#f0f0f0", "#d73027"),
                      breaks = c(-1.5, -0.5, 0.5, 1.5),
                      legend = FALSE, categorical = TRUE,
                      subtitle = "Blue↓ White=NS Red↑")
    legend("bottomright",
           legend = c("Decreasing", "Non-sig.", "Increasing"),
           fill = c("#4575b4", "#f0f0f0", "#d73027"),
           bty = "n", cex = 0.9, title = "Direction")
    
    title_prefix <- "PET"
    
  } else {
    r_comb <- if (method == "vc") r_precip_comb_vc  else r_precip_comb_tfpw
    r_mag  <- if (method == "vc") r_precip_mag_vc   else r_precip_mag_tfpw
    precip_all <- all_results[variable == "Precipitation" & period == "annual" & !is_basin_average]
    pval_col <- if (method == "vc") "p_value_vc" else "p_value_tfpw"
    r_pval <- create_raster_from_table(precip_all, template_bc, pval_col)
    r_pval <- smooth_raster(r_pval)
    
    plot_raster_panel(r_comb, "A) Kendall's τ", c(-0.15, 0.15),
                      hcl.colors(101, "RdBu", rev = TRUE),
                      legend_title = "τ",
                      subtitle = "No pixels significant (p < 0.05)")
    
    plot_raster_panel(r_mag, "B) Sen's Slope", c(-2, 2),
                      hcl.colors(101, "RdBu", rev = TRUE),
                      legend_title = "mm/yr")
    
    plot(r_pval,
         main = "C) p-Value Distribution",
         cex.main = 1.0, font.main = 2,
         col = hcl.colors(101, "YlOrRd", rev = TRUE),
         breaks = seq(0, 1, by = 0.1),
         zlim = c(0, 1),
         axes = FALSE, box = FALSE, legend = TRUE,
         plg = list(cex = 0.8, title = "p-value"))
    plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
    mtext("All p > 0.05 → no significant trends", side = 3, line = -0.5,
          cex = 0.8, col = "grey40", font = 3)
    
    title_prefix <- "Precipitation"
  }
  
  ## FIX 8h: Outer title matches actual data range (1950-2025)
  mtext(sprintf("%s Annual Trends (1950–2025) | %s Method", title_prefix,
                ifelse(method == "vc", "Variance-Corrected MK", "TFPW-Corrected MK")),
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  dev.off()
  cat(sprintf("✅ %s %s map: %s\n", var_name, method, basename(pdf_path)))
  invisible(pdf_path)
}

####################################################################################
# PET vs PRECIPITATION COMPARISON FIGURE — TIME RANGE FIXED
####################################################################################
create_comparison_figure <- function() {
  pdf_path <- file.path(out_dir, "pet_vs_precipitation_comparison.pdf")
  pdf(pdf_path, width = 12, height = 8, family = "Helvetica")
  par(mfrow = c(2, 2), mar = c(2, 2, 2.5, 1), oma = c(1, 1, 3, 1))
  
  plot_raster_panel(r_pet_comb_vc, "PET – Significant Trend Strength", c(0, 0.25),
                    hcl.colors(101, "YlOrRd", rev = FALSE),
                    legend_title = "|τ|",
                    subtitle = "|τ| for sig. pixels (p<0.05); 0 = non-significant")
  
  plot_raster_panel(r_pet_mag_vc, "PET – Sen's Slope (mm/yr)", c(0, 0.025),
                    hcl.colors(101, "YlOrRd", rev = FALSE),
                    legend_title = "mm/yr")
  
  plot_raster_panel(r_precip_comb_vc, "Precipitation – Kendall's τ (all pixels)", c(-0.15, 0.15),
                    hcl.colors(101, "RdBu", rev = TRUE),
                    legend_title = "τ",
                    subtitle = "No pixel significant at p < 0.05 (raw τ shown)")
  
  plot_raster_panel(r_precip_tau_vc, "Precipitation – Kendall's τ (full range)", c(-0.3, 0.3),
                    hcl.colors(101, "RdBu", rev = TRUE),
                    legend_title = "τ",
                    subtitle = "Red = positive trend; blue = negative trend")
  
  ## FIX 8i: Outer title matches actual data range (1950-2025)
  mtext("Nechako Basin Climate Trends (1950–2025) | Variance-Corrected MK",
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

basin_ts_plot    <- create_basin_timeseries_plot()
pet_vc_fig       <- create_publication_figure("PET", "vc")
pet_tfpw_fig     <- create_publication_figure("PET", "tfpw")
precip_vc_fig    <- create_publication_figure("Precipitation", "vc")
precip_tfpw_fig  <- create_publication_figure("Precipitation", "tfpw")
comp_fig         <- create_comparison_figure()
regime_fig       <- create_regime_shift_maps()
spectral_fig     <- create_spectral_analysis_maps()
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
cat("\nNEW/FIXED ENHANCED OUTPUTS:\n")
if (!is.null(regime_fig))
  cat("  • regime_shifts_changepoints.pdf         ← PELT changepoint detection (FIXED)\n")
if (!is.null(spectral_fig))
  cat("  • spectral_analysis.pdf                  ← Periodic patterns in Pr/PET\n")
cat("  • spatial_patterns.pdf                   ← Spatial variability analysis\n")
cat("  • method_comparison_vc_vs_tfpw.png       ← Method comparison (VC vs TFPW) (FIXED)\n")
cat("  • summary_statistics.csv                 ← Comprehensive statistics table\n")
cat("\n--- SUMMARY OF FIXES APPLIED ---\n")
cat("  1. Direction maps: filled NA/filtered pixels → 'Non-sig.' (no more spatial holes)\n")
cat("  2. Regime shift maps: pixels with no detected shift shown as 'No shift' category\n")
cat("  3. Precipitation 'Combined Trend Strength' → replaced with raw Kendall's τ\n")
cat("  4. Comparison figure: precipitation panel title and subtitle clarified\n")
cat("  5. Method comparison Panel B: % labels now computed per-variable (not globally)\n")
cat("  6. Method comparison Panel C: VC vs TFPW bars colored and labeled by method;\n")
cat("     0% precipitation result annotated and explained\n")
cat("  7. All figures: subtitles added to ambiguous panels for clarity\n")
cat("  8. TIME RANGE: ALL instances corrected from 1980-2055 to 1950-2025\n")