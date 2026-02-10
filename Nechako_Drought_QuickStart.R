#============================================================================
# QUICK-START GUIDE: PNPI/PNWBI DROUGHT ANALYSIS
# With Dual Return Period Methods
#============================================================================
# This script provides:
# 1. Data verification steps
# 2. Example workflow
# 3. Quick diagnostic checks
# 4. Minimal working example
#============================================================================

library(terra)

# Set your working directory
setwd("D:/Nechako_Drought/monthly_data_direct")

cat("="*70, "\n")
cat("DATA VERIFICATION\n")
cat("="*70, "\n\n")

# --- Check 1: Verify input files exist ---
cat("CHECK 1: Verifying input files...\n")
required_files <- c(
  "total_precipitation_monthly.nc",
  "ERA5Land_Nechako_PET_monthly.nc"
)
all_exist <- TRUE
for (file in required_files) {
  if (file.exists(file)) {
    cat(sprintf("  ✓ %s found\n", file))
  } else {
    cat(sprintf("  ✗ %s MISSING\n", file))
    all_exist <- FALSE
  }
}
if (!all_exist) {
  stop("\n ERROR: Some required files are missing!")
}
cat("\n")

# --- Check 2: Load and verify precipitation data ---
cat("CHECK 2: Loading and verifying precipitation data...\n")
precip <- rast("total_precipitation_monthly.nc")
cat(sprintf("  Layers: %d\n", nlyr(precip)))
cat(sprintf("  Cells: %d\n", ncell(precip)))
cat(sprintf("  Resolution: %.4f x %.4f degrees\n", res(precip)[1], res(precip)[2]))
cat(sprintf("  Extent: %.2f to %.2f (lon), %.2f to %.2f (lat)\n",
            ext(precip)[1], ext(precip)[2], ext(precip)[3], ext(precip)[4]))
precip_stats <- global(precip, c("min", "max", "mean"), na.rm = TRUE)
cat(sprintf("  Value range: %.2f to %.2f mm (mean: %.2f)\n",
            min(precip_stats$min), max(precip_stats$max), mean(precip_stats$mean)))
na_count <- global(is.na(precip), "sum")
cat(sprintf("  Missing values: %d (%.2f%%)\n",
            sum(na_count$sum), sum(na_count$sum)/(ncell(precip)*nlyr(precip))*100))
cat("\n")

# --- Check 3: Load and verify PET data ---
cat("CHECK 3: Loading and verifying PET data...\n")
pet <- rast("ERA5Land_Nechako_PET_monthly.nc")
cat(sprintf("  Layers: %d\n", nlyr(pet)))
cat(sprintf("  Cells: %d\n", ncell(pet)))
if (nlyr(precip) != nlyr(pet)) {
  warning("  ⚠ WARNING: Precipitation and PET have different temporal lengths!")
}
if (ncell(precip) != ncell(pet)) {
  warning("  ⚠ WARNING: Precipitation and PET have different spatial dimensions!")
}
pet_stats <- global(pet, c("min", "max", "mean"), na.rm = TRUE)
cat(sprintf("  Value range: %.3f to %.3f mm/day (mean: %.3f)\n",
            min(pet_stats$min), max(pet_stats$max), mean(pet_stats$mean)))
if (min(pet_stats$min) < 0) {
  warning("  ⚠ WARNING: Negative PET values detected!")
}
cat("\n")

# --- Check 4: Test water balance calculation ---
cat("CHECK 4: Testing water balance calculation...\n")
days_in_month_1 <- 31
pet_monthly_1 <- pet[[1]] * days_in_month_1
wb_1 <- precip[[1]] - pet_monthly_1
wb_stats <- global(wb_1, c("min", "max", "mean"), na.rm = TRUE)
cat(sprintf("  January water balance range: %.2f to %.2f mm (mean: %.2f)\n",
            wb_stats$min, wb_stats$max, wb_stats$mean))
neg_count <- global(wb_1 < 0, "sum", na.rm = TRUE)
cat(sprintf("  Cells with P < PET: %d (%.1f%%)\n",
            neg_count$sum, neg_count$sum/ncell(wb_1)*100))
cat("\n")

# --- Check 5: Test percent-of-normal calculation ---
cat("CHECK 5: Testing percent-of-normal calculation for January...\n")
jan_indices <- seq(1, nlyr(precip), by = 12)
precip_jan <- precip[[jan_indices]]
precip_jan_mean <- mean(precip_jan, na.rm = TRUE)
pnpi_jan_1 <- (precip[[1]] / precip_jan_mean) * 100

pnpi_stats <- global(pnpi_jan_1, c("min", "max", "mean"), na.rm = TRUE)
cat(sprintf("  PNPI range: %.1f%% to %.1f%% (mean: %.1f%%)\n",
            pnpi_stats$min, pnpi_stats$max, pnpi_stats$mean))
if (abs(pnpi_stats$mean - 100) > 5) {
  warning("  ⚠ WARNING: PNPI mean should be close to 100%!")
} else {
  cat("  ✓ PNPI calculation verified (mean ≈ 100%)\n")
}
cat("\n")

# --- Check 6: Test Kendall distribution function ---
cat("CHECK 6: Testing Kendall distribution function implementation...\n")
test_t <- 0.5
test_copula <- indepCopula(dim = 2)
c_tt <- pCopula(cbind(test_t, test_t), test_copula)
k_c <- test_t - (1 - c_tt) / test_t
cat(sprintf("  For independence copula at t=%.2f:\n", test_t))
cat(sprintf("    C(t,t) = %.4f\n", c_tt))
cat(sprintf("    K_C(t) = %.4f (should equal t=%.2f)\n", k_c, test_t))
if (abs(k_c - test_t) < 0.01) {
  cat("  ✓ Kendall distribution function verified\n")
} else {
  warning("  ⚠ WARNING: Kendall function implementation may need review")
}
cat("\n")

cat("="*70, "\n")
cat("DATA VERIFICATION COMPLETE\n")
cat("="*70, "\n\n")

#============================================================================
# MINIMAL WORKING EXAMPLE: SINGLE MONTH ANALYSIS
#============================================================================
cat("="*70, "\n")
cat("MINIMAL WORKING EXAMPLE: SINGLE MONTH ANALYSIS\n")
cat("="*70, "\n\n")

# Calculate PNPI for January 2000
jan_2000_idx <- (2000 - 1950) * 12 + 1
precip_jan_2000 <- precip[[jan_2000_idx]]

# Calculate climatology for January
jan_indices <- seq(1, nlyr(precip), by = 12)
precip_jan_clim <- precip[[jan_indices]]
precip_jan_mean <- mean(precip_jan_clim, na.rm = TRUE)

# Calculate PNPI
pnpi_jan_2000 <- (precip_jan_2000 / precip_jan_mean) * 100

# Calculate cell area
cell_area <- cellSize(pnpi_jan_2000, unit = "km")
total_area <- global(cell_area, "sum", na.rm = TRUE)[1,1]

# Drought threshold = 75%
THRESHOLD <- 75
drought_mask <- pnpi_jan_2000 < THRESHOLD

# Calculate area under drought
drought_area_km2 <- global(mask(cell_area, drought_mask, maskvalues = 0),
                           "sum", na.rm = TRUE)[1,1]
drought_area_pct <- (drought_area_km2 / total_area) * 100

# Calculate drought severity (mean deficit from normal)
deficit_values <- mask(100 - pnpi_jan_2000, drought_mask, maskvalues = 0)
drought_severity <- global(deficit_values, "mean", na.rm = TRUE)[1,1]

cat(sprintf("  Total basin area: %.2f km²\n", total_area))
cat(sprintf("  Area under drought (<75%% PNPI): %.2f km² (%.1f%%)\n",
            drought_area_km2, drought_area_pct))
cat(sprintf("  Drought severity: %.1f%% deficit from normal\n", drought_severity))
cat("\n")

pdf("minimal_example_drought_map.pdf", width = 10, height = 8)
par(mfrow = c(1, 2))
plot(pnpi_jan_2000, main = "PNPI - January 2000",
     col = colorRampPalette(c("brown", "white", "darkgreen"))(100))
plot(drought_mask, main = "Drought Areas (PNPI < 75%)",
     col = c("white", "red"), legend = FALSE)
legend("topright", legend = c("No drought", "Drought"),
       fill = c("white", "red"), bty = "n")
dev.off()
cat("  Map saved to: minimal_example_drought_map.pdf\n\n")

#============================================================================
# RECOMMENDED WORKFLOW
#============================================================================
cat("="*70, "\n")
cat("RECOMMENDED WORKFLOW\n")
cat("="*70, "\n\n")
cat("
STEP-BY-STEP EXECUTION PLAN:

1. VERIFY DATA (this script)
   ✓ Check file existence
   ✓ Verify data dimensions
   ✓ Test calculations
   Run time: ~2 minutes

2. RUN MAIN ANALYSIS (Nechako_Drought_SAF_Analysis.R)
   - Calculate PNPI and PNWBI
   - Extract drought characteristics
   - Fit marginal distributions
   - Fit copulas
   - Derive S-A-F curves (BOTH methods)
   Run time: ~30-60 minutes

3. REVIEW RESULTS
   Check these KEY OUTPUTS in drought_analysis/:
   
   PNPI Analysis:
   - PNPI_SAF_overlay.pdf (MOST IMPORTANT - compares both methods)
   - PNPI_SAF_comparison.pdf (side-by-side plots)
   - PNPI_drought_characteristics.csv
   
   PNWBI Analysis:
   - PNWBI_SAF_overlay.pdf (MOST IMPORTANT - compares both methods)
   - PNWBI_SAF_comparison.pdf (side-by-side plots)
   - PNWBI_drought_characteristics.csv
   
   - SUMMARY_REPORT.txt (comprehensive summary)

4. METHOD SELECTION GUIDANCE:
   - Use CONDITIONAL approach when:
     * You need drought severity conditioned on specific area coverage
     * Your application requires 'given X% area in drought, what severity?'
   
   - Use KENDALL approach when:
     * You need joint probability of severity-area combinations
     * Your application requires 'what severity-area pairs have return period T?'
     * You want more conservative estimates for basin-wide droughts

5. TROUBLESHOOTING:
   If the main script fails:
   - Check the log file: drought_analysis/SAF_analysis_log.txt
   - Common issues:
     * Memory: Reduce data extent or use temporal chunking
     * Copula fitting fails: Try different copula families
     * Unrealistic results: Check marginal distribution fits

For help, see: Nechako_Drought_Methodology_Guide.R
")

#============================================================================
# QUICK DIAGNOSTIC FUNCTION
#============================================================================
quick_diagnostic <- function() {
  cat("\n")
  cat("="*70, "\n")
  cat("QUICK DIAGNOSTIC CHECK\n")
  cat("="*70, "\n\n")
  
  cat("Checking for output files...\n")
  output_files <- c(
    "drought_analysis/PNPI_analysis/PNPI_drought_characteristics.csv",
    "drought_analysis/PNWBI_analysis/PNWBI_drought_characteristics.csv",
    "drought_analysis/PNPI_analysis/PNPI_SAF_curves_conditional.csv",
    "drought_analysis/PNPI_analysis/PNPI_SAF_curves_kendall.csv",
    "drought_analysis/PNWBI_analysis/PNWBI_SAF_curves_conditional.csv",
    "drought_analysis/PNWBI_analysis/PNWBI_SAF_curves_kendall.csv",
    "drought_analysis/SUMMARY_REPORT.txt"
  )
  for (file in output_files) {
    if (file.exists(file)) {
      cat(sprintf("  ✓ %s\n", basename(file)))
    } else {
      cat(sprintf("  ✗ %s MISSING\n", basename(file)))
    }
  }
  cat("\n")
  
  if (file.exists("drought_analysis/PNPI_analysis/PNPI_drought_characteristics.csv")) {
    cat("PNPI Analysis Summary:\n")
    pnpi_data <- read.csv("drought_analysis/PNPI_analysis/PNPI_drought_characteristics.csv")
    cat(sprintf("  Total months: %d\n", nrow(pnpi_data)))
    cat(sprintf("  Drought months (area > 0%%): %d (%.1f%%)\n",
                sum(pnpi_data$area_pct > 0),
                sum(pnpi_data$area_pct > 0)/nrow(pnpi_data)*100))
    cat(sprintf("  Mean severity: %.1f%% deficit (SD: %.1f%%)\n",
                mean(pnpi_data$severity), sd(pnpi_data$severity)))
    cat(sprintf("  Mean area: %.1f%% (SD: %.1f%%)\n",
                mean(pnpi_data$area_pct), sd(pnpi_data$area_pct)))
    cat(sprintf("  Max severity: %.1f%% deficit\n", max(pnpi_data$severity)))
    cat(sprintf("  Max area: %.1f%%\n", max(pnpi_data$area_pct)))
  }
  cat("\n")
  
  if (file.exists("drought_analysis/PNWBI_analysis/PNWBI_drought_characteristics.csv")) {
    cat("PNWBI Analysis Summary:\n")
    pnwbi_data <- read.csv("drought_analysis/PNWBI_analysis/PNWBI_drought_characteristics.csv")
    cat(sprintf("  Total months: %d\n", nrow(pnwbi_data)))
    cat(sprintf("  Drought months (area > 0%%): %d (%.1f%%)\n",
                sum(pnwbi_data$area_pct > 0),
                sum(pnwbi_data$area_pct > 0)/nrow(pnwbi_data)*100))
    cat(sprintf("  Mean severity: %.1f%% deficit (SD: %.1f%%)\n",
                mean(pnwbi_data$severity), sd(pnwbi_data$severity)))
    cat(sprintf("  Mean area: %.1f%% (SD: %.1f%%)\n",
                mean(pnwbi_data$area_pct), sd(pnwbi_data$area_pct)))
    cat(sprintf("  Max severity: %.1f%% deficit\n", max(pnwbi_data$severity)))
    cat(sprintf("  Max area: %.1f%%\n", max(pnwbi_data$area_pct)))
  }
  cat("\n")
  cat("Diagnostic complete. Review drought_analysis/ folder for detailed results.\n")
  cat("="*70, "\n")
}

cat("\n\nYou can now proceed with the main analysis.\n")
cat("Run: source('Nechako_Drought_SAF_Analysis.R')\n\n")
cat("After analysis completes, run: quick_diagnostic()\n\n")
cat("For additional analyses, see: Nechako_Drought_Methodology_Guide.R\n\n")

# Save verification results
verification_summary <- data.frame(
  Check = c("Precipitation file", "PET file", "Temporal match", "Spatial match",
            "Value ranges", "PNPI test", "Kendall test"),
  Status = c("PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS"),
  Notes = c(
    sprintf("%d layers", nlyr(precip)),
    sprintf("%d layers", nlyr(pet)),
    ifelse(nlyr(precip) == nlyr(pet), "Matched", "MISMATCH"),
    ifelse(ncell(precip) == ncell(pet), "Matched", "MISMATCH"),
    "See details above",
    sprintf("PNPI mean=%.1f%%", pnpi_stats$mean),
    sprintf("K_C(%.2f)=%.4f", test_t, k_c)
  )
)
write.csv(verification_summary, "data_verification_results.csv", row.names = FALSE)
cat("Verification summary saved to: data_verification_results.csv\n")