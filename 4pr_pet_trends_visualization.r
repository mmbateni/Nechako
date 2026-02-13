####################################################################################
# PUBLICATION-QUALITY TREND MAPS + BASIN AVERAGE TIME SERIES
# Purpose: Generate clean PDF maps AND basin-averaged time series plots
# Key: Native resolution maps + integrated time series visualization
####################################################################################

library(terra)
library(data.table)
library(sf)
library(zoo)  # For rolling mean

setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "trend_analysis_pr_pet"

# ===== LOAD CORE DATA =====
cat("📦 Loading results, basin boundary, and metadata...\n")
all_results   <- readRDS(file.path(out_dir, "all_results.rds"))
basin_boundary <- readRDS(file.path(out_dir, "basin_boundary.rds"))
metadata <- readRDS(file.path(out_dir, "analysis_metadata.rds"))

# Extract basin-averaged time series
basin_avg_monthly <- metadata$basin_avg_monthly
basin_avg_annual <- metadata$basin_avg_annual

# ===== CREATE NATIVE-RESOLUTION TEMPLATE =====
pet_annual_sample <- all_results[variable == "PET" & period == "annual" & !filtered_vc][1:100]
dx <- median(diff(sort(unique(pet_annual_sample$x))), na.rm = TRUE)
dy <- median(diff(sort(unique(pet_annual_sample$y))), na.rm = TRUE)
resolution_bc <- min(dx, dy, na.rm = TRUE) * 0.95

basin_ext <- ext(basin_boundary)
template_bc <- rast(ext = basin_ext, 
                    resolution = resolution_bc, 
                    crs = "EPSG:3005",
                    names = "trend")

cat(sprintf("✓ Created native-resolution template (%.1f m, %d x %d cells)\n", 
            resolution_bc, ncol(template_bc), nrow(template_bc)))

# ===== RASTER GENERATION (nearest neighbor) =====
create_raster_from_table <- function(results_dt, template, value_col) {
  pts <- vect(results_dt, geom = c("x", "y"), crs = crs(template))
  r <- rasterize(pts, template, field = value_col, touches = TRUE)
  mask(r, vect(basin_boundary))
}

# ===== PREPARE ANNUAL RESULTS (PIXELS ONLY) =====
cat("⚙️  Preparing trend metrics for pixel-level analysis...\n")
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

# ===== RASTER GENERATION =====
cat("🖼️  Generating rasters (native resolution)...\n")
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

# ===== LIGHT SMOOTHING =====
smooth_raster <- function(r) {
  if (all(is.na(values(r)))) return(r)
  focal(r, w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
}

# ===== PLOTTING HELPER =====
plot_raster_panel <- function(r, main, zlim = NULL, col, 
                              breaks = NULL, legend = TRUE, categorical = FALSE) {
  if (!categorical) r <- smooth_raster(r)
  
  if (!is.null(breaks)) {
    plot(r, main = main, cex.main = 1.0, col = col, breaks = breaks,
         axes = FALSE, box = FALSE, legend = legend, plg = list(cex = 0.9),
         colNA = NA)
  } else {
    plot(r, main = main, cex.main = 1.0, col = col, zlim = zlim,
         axes = FALSE, box = FALSE, legend = legend, plg = list(cex = 0.9),
         colNA = NA)
  }
  plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
}

# ===== BASIN AVERAGE TIME SERIES PLOT =====
create_basin_timeseries_plot <- function() {
  cat("📈 Generating basin-averaged time series plot...\n")
  pdf_path <- file.path(out_dir, "basin_average_timeseries.pdf")
  pdf(pdf_path, width = 12, height = 7, family = "Helvetica")
  
  # Extract basin average results for annotation
  basin_annual_results <- all_results[is_basin_average == TRUE & period == "annual"]
  pr_res <- basin_annual_results[variable == "Precipitation"]
  pet_res <- basin_annual_results[variable == "PET"]
  
  # Calculate 12-month rolling means for smoother visualization
  pr_rolling <- rollmean(basin_avg_monthly$precip_mm_month, k = 12, fill = NA, align = "center")
  pet_rolling <- rollmean(basin_avg_monthly$pet_mm_month, k = 12, fill = NA, align = "center")
  
  # Create dual-axis plot
  par(mar = c(5, 5, 4, 5) + 0.1)  # Extra space for right axis
  
  # Plot precipitation (monthly + rolling mean)
  plot(basin_avg_monthly$date, basin_avg_monthly$precip_mm_month, 
       type = "h", col = rgb(0.2, 0.4, 0.8, 0.3), lwd = 2,
       xlab = "Year", ylab = "Precipitation (mm/month)", 
       main = "Nechako Basin Climate Time Series (1980-2055)\nSpatially Averaged Monthly Values",
       ylim = c(0, max(basin_avg_monthly$precip_mm_month, na.rm = TRUE) * 1.1),
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  lines(basin_avg_monthly$date, pr_rolling, col = "blue", lwd = 2.5)
  
  # Add PET on secondary axis
  par(new = TRUE)
  plot(basin_avg_monthly$date, basin_avg_monthly$pet_mm_month,
       type = "l", col = rgb(0.8, 0.2, 0.2, 0.7), lwd = 2,
       axes = FALSE, xlab = "", ylab = "", 
       ylim = c(0, max(basin_avg_monthly$pet_mm_month, na.rm = TRUE) * 1.1))
  lines(basin_avg_monthly$date, pet_rolling, col = "red", lwd = 2.5)
  
  # Add secondary axis
  axis(side = 4, at = pretty(range(basin_avg_monthly$pet_mm_month, na.rm = TRUE)),
       col = "red", col.axis = "red", cex.axis = 1.1)
  mtext("PET (mm/month)", side = 4, line = 3, col = "red", cex = 1.2)
  
  # Add trend statistics as legend
  pr_label <- if (pr_res$filtered_vc) "Precipitation (filtered)" else 
    sprintf("Precipitation: %.2f mm/yr (p = %.3f)", pr_res$sl_vc, pr_res$p_value_vc)
  pet_label <- if (pet_res$filtered_vc) "PET (filtered)" else 
    sprintf("PET: %.2f mm/yr (p = %.3f)", pet_res$sl_vc, pet_res$p_value_vc)
  
  legend("topright", 
         legend = c(pr_label, pet_label, "12-mo Rolling Mean"),
         col = c("blue", "red", "gray40"), 
         lwd = c(2.5, 2.5, 1.5),
         lty = c(1, 1, 2),
         bty = "n", cex = 1.1)
  
  # Add period markers
  abline(v = as.Date(paste0(c(2000, 2020), "-01-01")), col = "gray40", lty = 2, lwd = 1)
  text(x = as.Date("2010-01-01"), y = par("usr")[4] * 0.95, "2000-2020", col = "gray40", cex = 0.9)
  
  dev.off()
  cat(sprintf("✅ Basin time series plot: %s\n", basename(pdf_path)))
  invisible(pdf_path)
}

# ===== MAIN FIGURE GENERATOR =====
create_publication_figure <- function(var_name, method = "vc") {
  pdf_path <- file.path(out_dir, sprintf("%s_%s_publication_map.pdf", tolower(var_name), method))
  pdf(pdf_path, width = 10, height = 6.67, family = "Helvetica")
  
  par(mfrow = c(1, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  if (var_name == "PET") {
    r_comb <- if (method == "vc") r_pet_comb_vc else r_pet_comb_tfpw
    r_dir  <- if (method == "vc") r_pet_dir_vc  else r_pet_dir_tfpw
    r_mag  <- if (method == "vc") r_pet_mag_vc  else r_pet_mag_tfpw
    
    plot_raster_panel(r_comb, "Combined Trend Strength", c(0, 0.25), 
                      hcl.colors(101, "YlOrRd", rev = FALSE))
    plot_raster_panel(r_mag, "Sen's Slope (mm/yr)", c(0, 1.5), 
                      hcl.colors(101, "RdBu", rev = TRUE))
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
                      hcl.colors(101, "RdBu", rev = TRUE))
    plot_raster_panel(r_mag, "Sen's Slope (mm/yr)", c(-2, 2), 
                      hcl.colors(101, "RdBu", rev = TRUE))
    plot_raster_panel(r_tau, "Kendall's Tau", c(-0.3, 0.3), 
                      hcl.colors(101, "RdBu", rev = TRUE))
    title_prefix <- "Precipitation"
  }
  
  mtext(sprintf("%s Annual Trends (1980-2055) | %s Method", title_prefix,
                ifelse(method == "vc", "Variance-Corrected MK", "TFPW-Corrected MK")),
        outer = TRUE, cex = 1.1, font = 2, line = 0.5)
  dev.off()
  
  cat(sprintf("✅ %s %s map: %s\n", var_name, method, basename(pdf_path)))
  invisible(pdf_path)
}

# ===== COMPARISON FIGURE =====
create_comparison_figure <- function() {
  pdf_path <- file.path(out_dir, "pet_vs_precipitation_comparison.pdf")
  pdf(pdf_path, width = 12, height = 8, family = "Helvetica")
  
  par(mfrow = c(2, 2), mar = c(2, 2, 2.5, 1), oma = c(1, 1, 3, 1))
  
  plot_raster_panel(r_pet_comb_vc, "PET - Significant Trend Strength", c(0, 0.25), 
                    hcl.colors(101, "YlOrRd", rev = FALSE))
  plot_raster_panel(r_pet_mag_vc, "PET - Sen's Slope (mm/yr)", c(0, 1.5), 
                    hcl.colors(101, "RdBu", rev = TRUE))
  plot_raster_panel(r_precip_comb_vc, "Precipitation - Significant Trend Strength", c(-0.15, 0.15), 
                    hcl.colors(101, "RdBu", rev = TRUE))
  plot_raster_panel(r_precip_tau_vc, "Precipitation - Kendall's Tau", c(-0.3, 0.3), 
                    hcl.colors(101, "RdBu", rev = TRUE))
  
  mtext("Nechako Basin Climate Trends (1980-2055) | Variance-Corrected MK", 
        outer = TRUE, cex = 1.3, font = 2, line = 0.5)
  dev.off()
  
  cat(sprintf("✅ Comparison map: %s\n", basename(pdf_path)))
  invisible(pdf_path)
}

# ===== GENERATE ALL FIGURES =====
cat("\n========================================\n")
cat("🎨 CREATING VISUALIZATIONS\n")
cat("========================================\n\n")

# First create the basin time series plot (new feature)
basin_ts_plot <- create_basin_timeseries_plot()

# Then create the spatial maps
pet_vc_fig      <- create_publication_figure("PET", "vc")
pet_tfpw_fig    <- create_publication_figure("PET", "tfpw")
precip_vc_fig   <- create_publication_figure("Precipitation", "vc")
precip_tfpw_fig <- create_publication_figure("Precipitation", "tfpw")
comp_fig        <- create_comparison_figure()

cat("\n========================================\n")
cat("✅ ALL VISUALIZATIONS CREATED SUCCESSFULLY\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n\n", file.path(getwd(), out_dir)))
cat("Generated PDFs:\n")
cat("  • basin_average_timeseries.pdf          ← NEW: Basin-averaged Pr & PET time series\n")
cat("  • pet_vc_publication_map.pdf\n")
cat("  • pet_tfpw_publication_map.pdf\n")
cat("  • precipitation_vc_publication_map.pdf\n")
cat("  • precipitation_tfpw_publication_map.pdf\n")
cat("  • pet_vs_precipitation_comparison.pdf\n")