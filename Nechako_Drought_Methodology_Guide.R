#============================================================================
# METHODOLOGY GUIDE: PNPI/PNWBI DROUGHT SEVERITY-AREA-FREQUENCY ANALYSIS
# Detailed Explanation and Helper Functions
#============================================================================
# This document provides:
# 1. Detailed explanation of PNPI/PNWBI methodology
# 2. Explanation of dual return period methods
# 3. Helper functions for additional analyses
# 4. Troubleshooting guide
#============================================================================

#============================================================================
# CONCEPTUAL FRAMEWORK: PNPI & PNWBI
#============================================================================
"""
PERCENT OF NORMAL INDICES: WHY AND HOW

Why use PNPI/PNWBI instead of SPI/SPEI?

1. INTERPRETABILITY:
   - PNPI = 60% means precipitation is 60% of long-term average
   - More intuitive for water managers and stakeholders
   - No statistical transformation needed for communication

2. DROUGHT CLASSIFICATION (WMO 2016; Heim 2002):
   Category          | PNPI Range | Description
   ------------------|------------|------------------------
   Extremely Wet     | >150%      | Exceptionally wet
   Very Wet          | 125-150%   | Much above normal
   Wet               | 105-125%   | Above normal
   Near Normal       | 90-105%    | Normal conditions
   Abnormally Dry    | 75-90%     | Drier than normal
   Moderate Drought  | 60-75%     | Moderate drought
   Severe Drought    | 40-60%     | Severe drought
   Extreme Drought   | <40%       | Exceptional drought

3. SEVERITY CALCULATION:
   Unlike SPI (where severity = mean of negative values), PNPI severity is:
   
   Drought Severity = Mean(100 - PNPI) for cells in drought
   
   Example: 
   - PNPI = 65% in drought area → Severity = 35% deficit from normal
   - PNPI = 45% in drought area → Severity = 55% deficit from normal

4. THRESHOLD SELECTION:
   Default threshold = 75% (moderate drought onset)
   Can be adjusted based on application:
   - 90%: Early warning (abnormally dry)
   - 75%: Moderate drought (standard for agricultural impacts)
   - 60%: Severe drought (water restrictions typically begin)
   - 40%: Extreme drought (emergency conditions)

5. LIMITATIONS:
   - Sensitive to mean calculation period (use ≥30 years)
   - Less robust for non-normal distributions than SPI
   - Doesn't account for temporal persistence like SPI at longer time scales
   - Best used for short-term (monthly) drought monitoring

REFERENCES:
- Heim, R.R. (2002). A review of twentieth-century drought indices used in the USA.
  Bulletin of the American Meteorological Society, 83(8), 1149-1165.
- WMO (2016). Standardized Precipitation Index User Guide. WMO-No. 1090.
- Palmer, W.C. (1965). Meteorological drought. U.S. Weather Bureau Research Paper 45.
"""

#============================================================================
# DUAL RETURN PERIOD METHODS EXPLAINED
#============================================================================
"""
TWO APPROACHES FOR SAF CURVE DERIVATION

METHOD A: CONDITIONAL COPULA APPROACH (Amirataee et al. 2018)
---------------------------------------------------------------
Formula: T_S|A(s|a) = μ_T / [1 - C_S|A=a(s|a)]

Where:
- μ_T = average inter-arrival time between drought events
- C_S|A=a(s|a) = conditional CDF of severity given area = a
- C_S|A=a = ∂C(F_S(s), F_A(a)) / ∂F_A(a) [partial derivative]

Interpretation:
"Given that X% of the basin is under drought, what severity has return period T?"

Advantages:
- Directly answers operational questions about conditional drought risk
- Intuitive for water managers ("if 30% of basin is dry, how bad will it be?")
- Matches traditional univariate return period concepts

Limitations:
- Sensitive to marginal distribution fits
- May underestimate joint extremes in tail regions

METHOD B: KENDALL DISTRIBUTION FUNCTION APPROACH (Genest et al. 2009)
-----------------------------------------------------------------------
Formula: 
  K_C(t) = t - (1 - C(t,t))/t
  T_K = μ_T / (1 - K_C(t))

Where:
- K_C(t) = Kendall distribution function
- C(t,t) = copula value at diagonal point (t,t)
- t = level set parameter defining iso-probability curves

Interpretation:
"What severity-area combinations have joint return period T?"

Advantages:
- Captures full joint probability structure
- More robust to marginal distribution misspecification
- Better handles tail dependence for extreme events
- Conservative estimates for basin-wide droughts

Limitations:
- Less intuitive operational interpretation
- Requires numerical solution for t

KEY DIFFERENCE:
- Conditional approach: Focuses on "given area, what severity?"
- Kendall approach: Focuses on "what severity-area pairs occur together?"

PRACTICAL GUIDANCE:
- Use CONDITIONAL for operational drought monitoring/forecasting
- Use KENDALL for risk assessment of extreme basin-wide events
- Compare both methods to understand uncertainty in return period estimates

REFERENCES:
- Genest, C., Rémillard, B., & Beaudoin, D. (2009). Goodness-of-fit tests 
  for copulas: A review and a power study. Insurance: Mathematics and 
  Economics, 44(2), 199-213.
- Salvadori, G., & De Michele, C. (2004). Frequency analysis via copulas.
  Water Resources Research, 40(12).
- Shiau, J.T. (2006). Fitting drought duration and severity with 
  two-dimensional copulas. Water Resources Management, 20(5), 795-815.
"""

#============================================================================
# HELPER FUNCTION 1: ALTERNATIVE DROUGHT THRESHOLDS
#============================================================================
analyze_multiple_thresholds <- function(index_rast, cell_area, dates,
                                        thresholds = c(90, 75, 60, 40),
                                        index_name = "PNPI") {
  """
  Compare drought characteristics using different threshold definitions
  Thresholds follow WMO (2016) drought classification:
  - 90%: Abnormally dry
  - 75%: Moderate drought (default)
  - 60%: Severe drought
  - 40%: Extreme drought
  """
  results_list <- list()
  
  for (thresh in thresholds) {
    cat(sprintf("\nAnalyzing threshold: %.0f%% (PNPI < %.0f%%)\n", thresh, thresh))
    drought_chars <- extract_drought_characteristics_percent(
      index_rast, cell_area, threshold_pct = thresh
    )
    
    results_list[[as.character(thresh)]] <- drought_chars
    
    cat(sprintf("  Mean severity: %.1f%% deficit\n", mean(drought_chars$severity)))
    cat(sprintf("  Mean area: %.1f%%\n", mean(drought_chars$area_pct)))
    cat(sprintf("  Drought frequency: %.1f%%\n", 
                sum(drought_chars$area_pct > 0) / nrow(drought_chars) * 100))
  }
  
  # Create comparison plot
  pdf(sprintf("%s_threshold_comparison.pdf", index_name), width = 12, height = 8)
  par(mfrow = c(2, 2))
  
  # Plot 1: Time series of area for different thresholds
  plot(dates, results_list[[1]]$area_pct, type = "n",
       ylim = c(0, 100), xlab = "Date", ylab = "Area under drought (%)",
       main = "Drought Area - Threshold Sensitivity")
  colors <- rainbow(length(thresholds))
  for (i in seq_along(thresholds)) {
    lines(dates, results_list[[i]]$area_pct, col = colors[i], lwd = 1.5)
  }
  legend("topright", legend = paste("Threshold <", thresholds, "%"),
         col = colors, lwd = 2, bty = "n")
  
  # Plot 2: Distribution of severity
  boxplot(lapply(results_list, function(x) x$severity),
          names = thresholds,
          xlab = "Threshold (%)", ylab = "Drought Severity (Deficit %)",
          main = "Severity Distribution by Threshold",
          col = colors)
  
  # Plot 3: Drought frequency
  freq <- sapply(results_list, function(x) sum(x$area_pct > 0) / nrow(x) * 100)
  barplot(freq, names.arg = thresholds,
          xlab = "Threshold (%)", ylab = "Drought Frequency (%)",
          main = "Drought Frequency by Threshold",
          col = colors)
  
  # Plot 4: Mean area under drought
  mean_area <- sapply(results_list, function(x) mean(x$area_pct))
  barplot(mean_area, names.arg = thresholds,
          xlab = "Threshold (%)", ylab = "Mean Area (%)",
          main = "Mean Area Under Drought",
          col = colors)
  
  dev.off()
  return(results_list)
}

#============================================================================
# HELPER FUNCTION 2: SEASONAL DROUGHT ANALYSIS
#============================================================================
seasonal_drought_analysis <- function(drought_chars, output_file = "seasonal_drought.pdf") {
  """
  Analyze seasonal patterns in drought characteristics
  Identifies which seasons are most prone to drought
  """
  # Add season classification
  drought_chars$season <- factor(
    cut(drought_chars$month,
        breaks = c(0, 2, 5, 8, 11, 12),
        labels = c("Winter", "Spring", "Summer", "Fall", "Winter")),
    levels = c("Winter", "Spring", "Summer", "Fall")
  )
  
  # Aggregate by season
  seasonal_stats <- aggregate(
    cbind(severity, area_pct) ~ season,
    data = drought_chars,
    FUN = function(x) c(mean = mean(x), sd = sd(x), max = max(x))
  )
  
  pdf(output_file, width = 12, height = 10)
  par(mfrow = c(2, 2))
  
  # Plot 1: Severity by season
  boxplot(severity ~ season, data = drought_chars,
          main = "Drought Severity by Season",
          xlab = "Season", ylab = "Severity (Deficit %)",
          col = c("lightblue", "lightgreen", "yellow", "orange"))
  
  # Plot 2: Area by season
  boxplot(area_pct ~ season, data = drought_chars,
          main = "Drought Area by Season",
          xlab = "Season", ylab = "Area (%)",
          col = c("lightblue", "lightgreen", "yellow", "orange"))
  
  # Plot 3: Monthly climatology of severity
  monthly_severity <- aggregate(severity ~ month, data = drought_chars, FUN = mean)
  plot(monthly_severity$month, monthly_severity$severity, type = "b",
       pch = 19, col = "darkred", lwd = 2,
       xlab = "Month", ylab = "Mean Severity (Deficit %)",
       main = "Monthly Climatology of Drought Severity",
       xaxt = "n")
  axis(1, at = 1:12, labels = month.abb)
  grid()
  
  # Plot 4: Monthly climatology of area
  monthly_area <- aggregate(area_pct ~ month, data = drought_chars, FUN = mean)
  plot(monthly_area$month, monthly_area$area_pct, type = "b",
       pch = 19, col = "darkblue", lwd = 2,
       xlab = "Month", ylab = "Mean Area (%)",
       main = "Monthly Climatology of Drought Area",
       xaxt = "n")
  axis(1, at = 1:12, labels = month.abb)
  grid()
  
  dev.off()
  return(seasonal_stats)
}

#============================================================================
# HELPER FUNCTION 3: SPATIAL DROUGHT PATTERNS
#============================================================================
map_drought_hotspots <- function(index_rast, dates, threshold = 75,
                                 output_dir = "drought_maps") {
  """
  Identify and map drought hotspots using percent-of-normal indices
  """
  dir.create(output_dir, showWarnings = FALSE)
  n_months <- nlyr(index_rast)
  
  # Calculate drought frequency for each pixel
  cat("Calculating drought frequency...\n")
  drought_freq <- sum(index_rast < threshold, na.rm = TRUE) / n_months * 100
  names(drought_freq) <- "Drought_Frequency_Pct"
  
  # Calculate mean severity during drought months
  cat("Calculating mean drought severity...\n")
  drought_severity_map <- app(index_rast, function(x) {
    drought_vals <- x[x < threshold]
    if (length(drought_vals) > 0) {
      return(mean(100 - drought_vals, na.rm = TRUE))  # Mean deficit
    } else {
      return(0)
    }
  })
  names(drought_severity_map) <- "Mean_Drought_Severity_Deficit_Pct"
  
  # Save maps
  writeCDF(drought_freq, file.path(output_dir, "drought_frequency_map.nc"),
           overwrite = TRUE)
  writeCDF(drought_severity_map, file.path(output_dir, "drought_severity_map.nc"),
           overwrite = TRUE)
  
  # Create plots
  pdf(file.path(output_dir, "drought_spatial_patterns.pdf"), width = 12, height = 10)
  par(mfrow = c(2, 2))
  
  plot(drought_freq, main = "Drought Frequency (%)",
       col = rev(heat.colors(100)))
  plot(drought_severity_map, main = "Mean Drought Severity (Deficit %)",
       col = rev(heat.colors(100)))
  
  # Histogram of frequency
  hist(values(drought_freq), breaks = 30,
       main = "Distribution of Drought Frequency",
       xlab = "Drought Frequency (%)",
       col = "skyblue", border = "white")
  
  # Histogram of severity
  hist(values(drought_severity_map), breaks = 30,
       main = "Distribution of Mean Drought Severity",
       xlab = "Mean Severity (Deficit %)",
       col = "coral", border = "white")
  
  dev.off()
  cat("Spatial drought maps saved to:", output_dir, "\n")
  
  return(list(
    frequency = drought_freq,
    severity = drought_severity_map
  ))
}

#============================================================================
# HELPER FUNCTION 4: METHOD COMPARISON PLOT
#============================================================================
compare_return_period_methods <- function(saf_conditional, saf_kendall, 
                                          index_name = "PNPI",
                                          output_file = "method_comparison.pdf") {
  """
  Create detailed comparison plot between conditional and Kendall methods
  """
  # Combine datasets
  saf_conditional$Method <- "Conditional"
  saf_kendall$Method <- "Kendall"
  saf_combined <- rbind(saf_conditional, saf_kendall)
  
  pdf(output_file, width = 10, height = 8)
  
  ggplot(saf_combined, aes(x = Area_pct, y = Severity, 
                           color = factor(ReturnPeriod_years),
                           linetype = Method, shape = Method)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_color_viridis_d(name = "Return Period\n(years)") +
    scale_linetype_manual(name = "Method", 
                          values = c("Conditional" = "solid", "Kendall" = "dashed")) +
    scale_shape_manual(name = "Method",
                       values = c("Conditional" = 19, "Kendall" = 17)) +
    labs(x = "Percent of Area Under Drought (%)",
         y = "Drought Severity (Deficit from Normal, %)",
         title = sprintf("%s: SAF Curves - Method Comparison", index_name),
         subtitle = "Solid lines: Conditional Copula | Dashed lines: Kendall Distribution") +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11)
    ) +
    scale_x_continuous(breaks = seq(0, 100, 10)) +
    scale_y_continuous(breaks = seq(0, ceiling(max(saf_combined$Severity)), 2))
  
  dev.off()
  cat(sprintf("Method comparison plot saved to: %s\n", output_file))
}

#============================================================================
# TROUBLESHOOTING GUIDE
#============================================================================
"""
COMMON ISSUES AND SOLUTIONS FOR PNPI/PNWBI ANALYSIS:

ISSUE: "PNPI mean not close to 100%"
SOLUTION: 
- Verify calendar month aggregation is correct
- Ensure sufficient years of data for each month (≥30 years recommended)
- Check for missing data in specific months
- Verify no division by zero in calculation

ISSUE: "All pixels show drought all the time"
SOLUTION:
- Check threshold value (should be 75% for moderate drought)
- Verify PNPI calculation: (observed / climatology) * 100
- Ensure climatology is calculated correctly per calendar month

ISSUE: "Kendall distribution function fails to converge"
SOLUTION:
- Check copula parameter values (should be within valid range)
- Try different starting values for optimization
- Verify u_severity and u_area are truly in [0,1] range
- Some copulas (e.g., Gumbel) require positive dependence

ISSUE: "Conditional copula optimization fails"
SOLUTION:
- Increase optimization interval range (e.g., c(0.001, 50) instead of c(0.001, 20))
- Check marginal distribution fits first - poor fits cause optimization failures
- Try different copula families

ISSUE: "Memory issues with large rasters"
SOLUTION:
- Process in temporal chunks using subset()
- Use terraOptions(memfrac = 0.9)
- Consider spatial subsetting for testing before full analysis

ISSUE: "Negative water balance values after constraint"
SOLUTION:
- Verify constraint is applied correctly: wb <- ifel(wb < 0, 0, wb)
- Check PET calculation (should be positive)
- This is expected behavior - constraint sets negatives to zero

CITATION TRACKING FOR PUBLICATIONS:
When publishing results, cite:

Main methodology: Amirataee et al. (2018)
PNPI methodology: WMO (2016); Heim (2002)
Water balance: Palmer (1965); Thornthwaite (1948)
Conditional approach: Shiau (2006)
Kendall approach: Genest et al. (2009); Salvadori & De Michele (2004)
Copula theory: Nelsen (2006); Sklar (1959)
Parameter estimation: Joe (1997)
ERA5-Land data: Muñoz-Sabater et al. (2021)
"""

#============================================================================
# END OF METHODOLOGY GUIDE
#============================================================================
cat("\n")
cat("="*70, "\n")
cat("METHODOLOGY GUIDE LOADED SUCCESSFULLY\n")
cat("="*70, "\n")
cat("\nAvailable helper functions:\n")
cat("  1. analyze_multiple_thresholds() - Test different drought thresholds\n")
cat("  2. seasonal_drought_analysis()    - Analyze seasonal patterns\n")
cat("  3. map_drought_hotspots()         - Create spatial drought maps\n")
cat("  4. compare_return_period_methods()- Compare conditional vs Kendall methods\n")
cat("\nUsage example:\n")
cat("  pnpi <- rast('drought_analysis/PNPI_monthly.nc')\n")
cat("  results <- analyze_multiple_thresholds(pnpi, cell_area_km2, dates)\n")
cat("\n")