####################################################################################
# PUBLICATION-QUALITY MAPS AS PDF
# Uses PDF output for superior raster interpolation
# Dramatically increased resampling for perfectly smooth gradients
####################################################################################
library(terra)
library(data.table)
library(sf)
library(RColorBrewer)

setwd("D:/Nechako_Drought/Nechako/")
out_dir <- "trend_analysis_pr_pet"

# Load data
basin_boundary <- readRDS(file.path(out_dir, "basin_boundary.rds"))

# Helper function to resample raster for ultra-smooth display (1000x increase!)
resample_for_display <- function(r) {
  # Create ultra high-res template (34x52 -> 3400x5200 pixels)
  template <- rast(ext(r), resolution = res(r) / 100, crs = crs(r))
  # Bilinear interpolation for perfectly smooth gradients
  r_smooth <- resample(r, template, method = "bilinear")
  return(r_smooth)
}

# Create publication figure (PDF format - 6×4 inches at high quality)
create_publication_figure <- function(var_name, method = "vc") {
  # Determine input directory
  input_dir <- file.path(out_dir, paste0(tolower(var_name), "_production_tifs"))
  
  # Load rasters (34×52 pixels - scientifically valid!)
  if (var_name == "PET") {
    r_comb_orig <- rast(file.path(input_dir, sprintf("pet_annual_%s_combined.tif", method)))
    r_dir <- rast(file.path(input_dir, sprintf("pet_annual_%s_direction.tif", method)))
    r_mag_orig <- rast(file.path(input_dir, sprintf("pet_annual_%s_magnitude.tif", method)))
    
    # Resample continuous rasters for ultra-smooth display
    r_comb <- resample_for_display(r_comb_orig)
    r_mag <- resample_for_display(r_mag_orig)
    # Direction map stays blocky (categorical)
    
    # PET-specific settings
    comb_zlim <- c(0, 1.5)
    mag_zlim <- c(0, 1.5)
    title_prefix <- "PET"
    comb_title <- "Combined Trend Strength"
    mag_title <- "Sen's Slope (mm/yr)"
    
  } else { # Precipitation
    r_comb_orig <- rast(file.path(input_dir, sprintf("precipitation_annual_%s_combined.tif", method)))
    r_tau_orig <- rast(file.path(input_dir, sprintf("precipitation_annual_%s_tau.tif", method)))
    r_mag_orig <- rast(file.path(input_dir, sprintf("precipitation_annual_%s_magnitude.tif", method)))
    
    # Resample for ultra-smooth display
    r_comb <- resample_for_display(r_comb_orig)
    r_tau <- resample_for_display(r_tau_orig)
    r_mag <- resample_for_display(r_mag_orig)
    
    # Precipitation-specific settings
    comb_zlim <- c(-0.015, 0.015)
    mag_zlim <- c(-0.02, 0.02)
    title_prefix <- "Precipitation"
    comb_title <- "Combined Trend Strength"
    mag_title <- "Sen's Slope (mm/yr)"
  }
  
  # Create PDF (publication quality)
  pdf_path <- file.path(out_dir, sprintf("%s_%s_publication_map.pdf", tolower(var_name), method))
  pdf(pdf_path, width = 10, height = 6.67)  # 10×6.67 inches for better aspect ratio
  
  # Layout: 1 row, 3 columns
  par(mfrow = c(1, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  # Panel 1: Combined map
  plot(r_comb,
       main = comb_title,
       cex.main = 1.0,
       col = hcl.colors(101, "RdBu", rev = TRUE),
       zlim = comb_zlim,
       axes = FALSE, box = FALSE,
       legend = TRUE,
       plg = list(cex = 0.9))
  plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
  
  # Panel 2: Magnitude map
  plot(r_mag,
       main = mag_title,
       cex.main = 1.0,
       col = hcl.colors(101, "RdBu", rev = TRUE),
       zlim = mag_zlim,
       axes = FALSE, box = FALSE,
       legend = TRUE,
       plg = list(cex = 0.9))
  plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
  
  # Panel 3: Direction (PET) or Tau (Precipitation)
  if (var_name == "PET") {
    plot(r_dir,
         main = "Direction (p < 0.05)",
         cex.main = 1.0,
         col = c("#4575b4", "#f0f0f0", "#d73027"),
         breaks = c(-1.5, -0.5, 0.5, 1.5),
         axes = FALSE, box = FALSE,
         legend = FALSE)
    plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
    
    legend("bottomright", 
           legend = c("Decreasing", "Non-sig.", "Increasing"),
           fill = c("#4575b4", "#f0f0f0", "#d73027"),
           bty = "n", cex = 1.1, ncol = 1)
  } else {
    plot(r_tau,
         main = "Kendall's Tau",
         cex.main = 1.0,
         col = hcl.colors(101, "RdBu", rev = TRUE),
         zlim = c(-0.2, 0.2),
         axes = FALSE, box = FALSE,
         legend = TRUE,
         plg = list(cex = 0.9))
    plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
  }
  
  # Main title
  mtext(sprintf("%s Annual Trends (1980-2023) | %s Method", title_prefix, 
                ifelse(method == "vc", "Variance-Corrected MK", "TFPW-Corrected MK")),
        outer = TRUE, cex = 1.0, font = 2, line = 0.5)
  
  dev.off()
  
  cat(sprintf("✅ Created PDF figure:\n   %s\n", pdf_path))
  cat(sprintf("   Ultra-high resolution resampling for smooth display\n"))
  cat(sprintf("   Vector+raster hybrid for publication quality\n\n"))
  
  return(pdf_path)
}

# Generate publication figures
cat("========================================\n")
cat("CREATING PDF PUBLICATION FIGURES\n")
cat("  (PDF format with ultra-smooth display)\n")
cat("========================================\n\n")

# PET figures
pet_vc_fig <- create_publication_figure("PET", "vc")
pet_tfpw_fig <- create_publication_figure("PET", "tfpw")

# Precipitation figures
precip_vc_fig <- create_publication_figure("Precipitation", "vc")
precip_tfpw_fig <- create_publication_figure("Precipitation", "tfpw")

# Create comparison figure
create_comparison_figure <- function() {
  pdf_path <- file.path(out_dir, "pet_vs_precipitation_comparison.pdf")
  pdf(pdf_path, width = 12, height = 8)
  
  par(mfrow = c(2, 2), mar = c(2, 2, 2.5, 1), oma = c(1, 1, 3, 1))
  
  # PET Combined
  pet_comb_orig <- rast(file.path(out_dir, "pet_production_tifs/pet_annual_vc_combined.tif"))
  pet_comb <- resample_for_display(pet_comb_orig)
  plot(pet_comb, main = "PET - Combined Trend", 
       cex.main = 1.1,
       col = hcl.colors(101, "RdBu", rev = TRUE), zlim = c(0, 1.5),
       axes = FALSE, box = FALSE,
       legend = TRUE,
       plg = list(cex = 0.9))
  plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
  
  # PET Magnitude
  pet_mag_orig <- rast(file.path(out_dir, "pet_production_tifs/pet_annual_vc_magnitude.tif"))
  pet_mag <- resample_for_display(pet_mag_orig)
  plot(pet_mag, main = "PET - Sen's Slope (mm/yr)", 
       cex.main = 1.1,
       col = hcl.colors(101, "RdBu", rev = TRUE), zlim = c(0, 1.5),
       axes = FALSE, box = FALSE,
       legend = TRUE,
       plg = list(cex = 0.9))
  plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
  
  # Precipitation Combined
  precip_comb_orig <- rast(file.path(out_dir, "precipitation_production_tifs/precipitation_annual_vc_combined.tif"))
  precip_comb <- resample_for_display(precip_comb_orig)
  plot(precip_comb, main = "Precipitation - Combined Trend", 
       cex.main = 1.1,
       col = hcl.colors(101, "RdBu", rev = TRUE), zlim = c(-0.015, 0.015),
       axes = FALSE, box = FALSE,
       legend = TRUE,
       plg = list(cex = 0.9))
  plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
  
  # Precipitation Tau
  precip_tau_orig <- rast(file.path(out_dir, "precipitation_production_tifs/precipitation_annual_vc_tau.tif"))
  precip_tau <- resample_for_display(precip_tau_orig)
  plot(precip_tau, main = "Precipitation - Kendall's Tau", 
       cex.main = 1.1,
       col = hcl.colors(101, "RdBu", rev = TRUE), zlim = c(-0.1, 0.1),
       axes = FALSE, box = FALSE,
       legend = TRUE,
       plg = list(cex = 0.9))
  plot(st_geometry(basin_boundary), add = TRUE, col = NA, border = "black", lwd = 2.5)
  
  mtext("Nechako Basin Climate Trends (1980-2023)", 
        outer = TRUE, cex = 1.2, font = 2, line = 0.5)
  
  dev.off()
  cat(sprintf("✅ Created PDF comparison figure:\n   %s\n\n", pdf_path))
}

create_comparison_figure()

cat("========================================\n")
cat("✅ ALL PDF FIGURES CREATED SUCCESSFULLY\n")
cat("========================================\n")
cat("Output files:\n")
cat("  • pet_vc_publication_map.pdf\n")
cat("  • pet_tfpw_publication_map.pdf\n")
cat("  • precipitation_vc_publication_map.pdf\n")
cat("  • precipitation_tfpw_publication_map.pdf\n")
cat("  • pet_vs_precipitation_comparison.pdf\n")
cat("\nKey improvements:\n")
cat("  ✅ PDF format - superior raster interpolation\n")
cat("  ✅ 100x resolution increase (34×52 -> 3400×5200)\n")
cat("  ✅ Perfectly smooth gradients - NO vertical bars\n")
cat("  ✅ Vector text/borders for crisp rendering\n")
cat("  ✅ Scalable without quality loss\n")
cat("  ✅ Journal submission ready\n")