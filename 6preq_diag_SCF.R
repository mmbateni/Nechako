##############################################
# SCF TIMESCALE DECISION ALGORITHM
# Automatically selects optimal SCF thresholds per timescale
# Based on snow persistence diagnostics and scoring rubric
##############################################

# ---- Libraries ----
library(terra)
library(ncdf4)
library(zoo)
library(dplyr)
library(ggplot2)

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")

scf_file <- "monthly_data_direct/snow_cover_monthly.nc"
swe_file <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
elev_file <- NULL  # optional: path to elevation raster (meters)
out_dir <- "scf_timescale_diagnostics"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("\n============================================================\n")
cat("SCF TIMESCALE DECISION ALGORITHM\n")
cat("Optimal threshold selection per timescale\n")
cat("============================================================\n\n")

# ---- Basic checks and load data ----
if (!file.exists(scf_file)) stop("SCF file not found: ", scf_file)
if (!file.exists(swe_file)) stop("SWE file not found: ", swe_file)

cat("===== READING NETCDF FILES =====\n")
scf_stack <- rast(scf_file)
swe_stack <- rast(swe_file)

# ---- Time extraction (robust) ----
extract_time_dimension <- function(rast_obj, file_path) {
  t <- time(rast_obj)
  if (!all(is.na(t))) return(as.Date(t))
  
  nc <- nc_open(file_path)
  on.exit(nc_close(nc))
  
  # Try common variable names
  tvname <- intersect(c("time", "valid_time", "times", "date"), names(nc$dim))
  if (length(tvname) > 0) {
    tv <- ncvar_get(nc, tvname[1])
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  
  # Fallback: monthly sequence using layer count
  n <- nlyr(rast_obj)
  start <- as.Date("1950-01-01")
  seq(start, by = "month", length.out = n)
}

dates_scf <- extract_time_dimension(scf_stack, scf_file)
dates_swe <- extract_time_dimension(swe_stack, swe_file)

if (length(dates_scf) != nlyr(scf_stack)) stop("SCF time mismatch")
if (length(dates_swe) != nlyr(swe_stack)) stop("SWE time mismatch")

terra::time(scf_stack) <- dates_scf
terra::time(swe_stack) <- dates_swe

cat(sprintf("SCF: %d x %d | Layers: %d\n", ncol(scf_stack), nrow(scf_stack), nlyr(scf_stack)))
cat(sprintf("SCF period: %s to %s\n", min(dates_scf), max(dates_scf)))
cat(sprintf("SWE: %d x %d | Layers: %d\n", ncol(swe_stack), nrow(swe_stack), nlyr(swe_stack)))
cat(sprintf("SWE period: %s to %s\n", min(dates_swe), max(dates_swe)))

# ---- Convert SCF to fraction if in percent ----
scf_max_val <- global(scf_stack[[1]], "max", na.rm = TRUE)$max
if (scf_max_val > 1.5) {
  cat("\n→ Converting SCF from percent to fraction (dividing by 100)\n")
  scf_stack <- scf_stack / 100
  cat(sprintf("  New range: %.3f to %.3f\n", 
              global(scf_stack[[1]], "min", na.rm = TRUE)$min,
              global(scf_stack[[1]], "max", na.rm = TRUE)$max))
}

# ---- Optional elevation ----
elev_rast <- NULL
if (!is.null(elev_file) && file.exists(elev_file)) {
  cat("\nLoading elevation raster...\n")
  elev_rast <- rast(elev_file)
  elev_rast <- project(elev_rast, crs(scf_stack))
}

# ---- Reproject SWE to SCF grid if needed ----
if (!compareGeom(scf_stack[[1]], swe_stack[[1]], stopOnError = FALSE)) {
  cat("\n→ Reprojecting SWE to SCF grid...\n")
  swe_stack <- resample(swe_stack, scf_stack, method = "bilinear")
  cat("→ SWE reprojected\n")
}

# ---- Convert SWE to mm if likely in meters ----
swe_sample_mean <- global(swe_stack[[1]], "mean", na.rm = TRUE)$mean
if (!is.na(swe_sample_mean) && swe_sample_mean < 10) {
  cat("\n→ Converting SWE from meters to mm (multiply by 1000)\n")
  swe_stack <- swe_stack * 1000
  cat(sprintf("  New sample mean: %.2f mm\n", 
              global(swe_stack[[1]], "mean", na.rm = TRUE)$mean))
}

### MODIFIED: Define snow-capable domain (≥5% SCF) BEFORE diagnostics
cat("\n===== DEFINING SNOW-CAPABLE DOMAIN (≥5% SCF) =====\n")
scf_mean <- mean(scf_stack, na.rm = TRUE)
snow_capable_mask_vals <- values(scf_mean) >= 0.05  # 5% threshold per Huning & AghaKouchak
total_pixels <- ncell(scf_stack)
snow_capable_pixels <- sum(snow_capable_mask_vals, na.rm = TRUE)
cat(sprintf("Snow-capable pixels (≥5%% SCF): %.1f%% of domain\n", 
            100 * snow_capable_pixels / total_pixels))

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Centered rolling sum for k months
centered_roll_sum <- function(vec, k) {
  if (k == 1) return(vec)
  zoo::rollapply(vec, width = k, FUN = sum, align = "center", fill = NA, na.rm = FALSE)
}

# Compute k-month cumulative SWE per pixel (returns matrix: ncell x nt)
compute_kmonth_swe_matrix <- function(swe_rast, dates, k) {
  ncell <- ncell(swe_rast)
  nt <- length(dates)
  
  # Read values as matrix (cells x time)
  swe_mat <- values(swe_rast, mat = TRUE)
  
  # Apply centered roll per row
  cat(sprintf("  Computing %d-month rolling sums for %d pixels...\n", k, ncell))
  res <- t(apply(swe_mat, 1, function(v) {
    if (all(is.na(v))) return(rep(NA_real_, nt))
    centered_roll_sum(v, k)
  }))
  
  return(res)
}

# Compute k-month SCF climatology per center month
compute_kmonth_scf_climatology <- function(scf_rast, dates, k) {
  months <- as.integer(format(dates, "%m"))
  centers <- 1:12
  scf_clim <- vector("list", 12)
  
  cat(sprintf("  Computing %d-month SCF climatology for 12 center months...\n", k))
  
  for (m in centers) {
    # Months in centered window: m - floor((k-1)/2) ... m + ceil((k-1)/2)
    half <- floor((k-1)/2)
    months_window <- ((m - half - 1):(m + (k - 1 - half) - 1)) %% 12 + 1
    idx <- which(months %in% months_window)
    
    if (length(idx) == 0) {
      scf_clim[[m]] <- rast(scf_rast[[1]])
      values(scf_clim[[m]]) <- NA
    } else {
      scf_clim[[m]] <- mean(scf_rast[[idx]], na.rm = TRUE)
    }
  }
  
  return(scf_clim)
}

### MODIFIED: compute_diagnostics() now returns BOTH absolute and relative fractions
compute_diagnostics <- function(scf_clim_list, swe_kmat, dates, threshold, 
                                snow_capable_mask_vals,  # NEW PARAMETER
                                stability_years = 0.50) {
  months <- as.integer(format(dates, "%m"))
  ncell <- nrow(swe_kmat)
  nt <- ncol(swe_kmat)
  centers <- 1:12
  
  retained_frac_abs <- numeric(12)  # Absolute: retained / total domain pixels
  retained_frac_rel <- numeric(12)  # Relative: retained / snow-capable pixels
  
  # Build mask per center month (logical vector per cell)
  mask_month <- vector("list", 12)
  
  for (m in centers) {
    r <- scf_clim_list[[m]]
    if (is.null(r)) {
      mask_month[[m]] <- rep(FALSE, ncell)
      retained_frac_abs[m] <- 0
      retained_frac_rel[m] <- 0
    } else {
      vals <- values(r)
      mask <- !is.na(vals) & (vals >= threshold)
      mask_month[[m]] <- mask
      
      # ABSOLUTE fraction (used for scoring)
      retained_frac_abs[m] <- sum(mask, na.rm = TRUE) / total_pixels
      
      # RELATIVE fraction (for transparency only)
      retained_frac_rel[m] <- ifelse(snow_capable_pixels > 0,
                                     sum(mask & snow_capable_mask_vals, na.rm = TRUE) / snow_capable_pixels,
                                     0)
    }
  }
  
  # For each cell, compute metrics using the mask for the center month of each time step
  years <- unique(format(dates, "%Y"))
  ny <- length(years)
  
  zero_frac_per_cell <- rep(NA_real_, ncell)
  nonzero_years_per_cell <- rep(NA_real_, ncell)
  stability_per_cell <- rep(NA_real_, ncell)
  
  for (cell in 1:ncell) {
    # Collect k-month SWE values for times where mask for that center month is TRUE
    valid_idx <- which(sapply(1:nt, function(t) {
      m <- months[t]
      mask_month[[m]][cell]
    }))
    
    if (length(valid_idx) == 0) {
      zero_frac_per_cell[cell] <- NA
      nonzero_years_per_cell[cell] <- 0
      stability_per_cell[cell] <- 0
      next
    }
    
    vals <- swe_kmat[cell, valid_idx]
    
    # Count years with positive sum (group by year)
    yrs <- format(dates[valid_idx], "%Y")
    pos_by_year <- tapply(vals, yrs, function(x) any(!is.na(x) & x > 0))
    pos_by_year[is.na(pos_by_year)] <- FALSE
    
    nonzero_years <- sum(pos_by_year)
    nonzero_years_per_cell[cell] <- nonzero_years
    zero_frac_per_cell[cell] <- 1 - (nonzero_years / length(unique(yrs)))
    stability_per_cell[cell] <- nonzero_years / length(unique(yrs))
  }
  
  zero_fraction_median <- median(zero_frac_per_cell, na.rm = TRUE)
  median_nonzero_years <- median(nonzero_years_per_cell, na.rm = TRUE)
  stability_pct <- 100 * mean(stability_per_cell >= stability_years, na.rm = TRUE)
  
  list(
    retained_frac_abs = retained_frac_abs,    # NEW
    retained_frac_rel = retained_frac_rel,    # NEW
    mean_retained_abs = mean(retained_frac_abs, na.rm = TRUE),  # NEW
    mean_retained_rel = mean(retained_frac_rel, na.rm = TRUE),  # NEW
    zero_fraction_median = zero_fraction_median,
    median_nonzero_years = median_nonzero_years,
    stability_pct = stability_pct,
    mask_month = mask_month,
    zero_frac_per_cell = zero_frac_per_cell,
    nonzero_years_per_cell = nonzero_years_per_cell,
    stability_per_cell = stability_per_cell
  )
}

# Scoring rubric (UNCHANGED - still uses absolute fraction for spatial context)
score_candidate <- function(diag, targets = list(nonzero_years = 20, 
                                                 zero_frac = 0.50, 
                                                 stability_pct = 50)) {
  # Normalize metrics to 0-1 (higher is better)
  s_nonzero <- pmin(diag$median_nonzero_years / targets$nonzero_years, 1)
  s_zero <- pmin((targets$zero_frac) / (diag$zero_fraction_median + 1e-9), 1)
  s_stab <- pmin(diag$stability_pct / targets$stability_pct, 1)
  
  # Retained area preference: moderate (not too small) - USES ABSOLUTE FRACTION
  mean_retained <- diag$mean_retained_abs  # CRITICAL: still absolute for scoring
  if (is.na(mean_retained)) mean_retained <- 0
  s_retained <- ifelse(mean_retained < 0.2, 0, 
                       ifelse(mean_retained > 0.95, 0.95, mean_retained))
  
  # Weights
  w_nonzero <- 0.30
  w_zero <- 0.25
  w_stab <- 0.20
  w_retained <- 0.15
  w_misc <- 0.10
  
  score <- w_nonzero * s_nonzero + w_zero * s_zero + 
    w_stab * s_stab + w_retained * s_retained + w_misc * 0.5
  
  list(score = score, 
       components = c(nonzero = s_nonzero, zero = s_zero, 
                      stability = s_stab, retained = s_retained))
}

# ==============================================================================
# DECISION ALGORITHM MAIN
# ==============================================================================

cat("\n===== CONFIGURING CANDIDATES =====\n")

# Candidate thresholds per timescale (adjustable)
candidates <- list(
  "1"  = c(0.01,0.03, 0.05, 0.075,0.10,0.15),
  "2"  = c(0.01,0.03, 0.05, 0.075,0.10,0.15),
  "3"  = c(0.01,0.03, 0.05, 0.075,0.10,0.15)
)

# Targets for scoring (adjustable)
targets <- list(nonzero_years = 20, zero_frac = 0.30, stability_pct = 50)

cat("Candidate thresholds:\n")
for (k in names(candidates)) {
  cat(sprintf("  %2s-month: %s\n", k, paste(candidates[[k]], collapse = ", ")))
}

cat("\nScoring targets:\n")
cat(sprintf("  Nonzero years: %d\n", targets$nonzero_years))
cat(sprintf("  Zero fraction: %.2f\n", targets$zero_frac))
cat(sprintf("  Stability: %.0f%%\n", targets$stability_pct))

# ==============================================================================
# PRECOMPUTE SWE K-MONTH MATRICES AND SCF CLIMATOLOGIES
# ==============================================================================

cat("\n===== PRECOMPUTING FOR EACH TIMESCALE =====\n")

timescales <- as.integer(names(candidates))
swe_kmats <- list()
scf_clims <- list()

for (k in timescales) {
  cat(sprintf("\nTimescale: %d months\n", k))
  swe_kmats[[as.character(k)]] <- compute_kmonth_swe_matrix(swe_stack, dates_swe, k)
  scf_clims[[as.character(k)]] <- compute_kmonth_scf_climatology(scf_stack, dates_scf, k)
}

# ==============================================================================
# RUN DIAGNOSTICS AND SCORING WITH TIE-BREAKING LOGIC
# ==============================================================================

cat("\n===== RUNNING DIAGNOSTICS AND SCORING =====\n")

results_all <- list()
tie_tolerance <- 0.001  # Consider scores tied if within this tolerance

for (k in timescales) {
  kstr <- as.character(k)
  swe_kmat <- swe_kmats[[kstr]]
  scf_clim_list <- scf_clims[[kstr]]
  cand_vec <- candidates[[kstr]]
  
  cat(sprintf("\n--- TIMESCALE: %d-month ---\n", k))
  
  res_k <- list()
  
  for (th in cand_vec) {
    cat(sprintf("  Testing SCF threshold = %.3f ... ", th))
    
    ### MODIFIED: Pass snow_capable_mask_vals to diagnostics
    diag <- compute_diagnostics(scf_clim_list, swe_kmat, dates_swe, th, 
                                snow_capable_mask_vals = snow_capable_mask_vals,  # ADDED
                                stability_years = 0.50)
    score <- score_candidate(diag, targets)
    
    res_k[[as.character(th)]] <- list(threshold = th, diag = diag, score = score)
    
    cat(sprintf("abs=%.1f%% rel=%.1f%% score=%.3f\n", 
                100 * diag$mean_retained_abs, 
                100 * diag$mean_retained_rel,
                score$score))
  }
  
  # Choose best by score WITH TIE-BREAKING: prefer smaller threshold when tied
  scores <- sapply(res_k, function(x) x$score$score)
  max_score <- max(scores)
  
  # Identify tied candidates (within tolerance)
  tied_thresholds <- as.numeric(names(res_k)[abs(scores - max_score) <= tie_tolerance])
  
  if (length(tied_thresholds) > 1) {
    # TIE-BREAKING RULE: Select SMALLER threshold (more conservative/inclusive)
    best_th <- min(tied_thresholds)
    cat(sprintf("\n→ TIE DETECTED for %d-month: thresholds %.3f–%.3f scored %.3f\n", 
                k, min(tied_thresholds), max(tied_thresholds), max_score))
    cat(sprintf("→ SELECTED SMALLER THRESHOLD: %.3f (conservative inclusion principle)\n", best_th))
  } else {
    best_th <- as.numeric(names(which.max(scores)))
    cat(sprintf("\n→ Best threshold for %d-month: %.3f (abs=%.1f%% rel=%.1f%% score=%.3f)\n", 
                k, best_th, 
                100 * res_k[[as.character(best_th)]]$diag$mean_retained_abs,
                100 * res_k[[as.character(best_th)]]$diag$mean_retained_rel,
                max_score))
  }
  
  results_all[[kstr]] <- list(candidates = res_k, best_threshold = best_th)
}

# ==============================================================================
# PRODUCE SUMMARY TABLE AND PLOTS
# ==============================================================================

cat("\n===== CREATING SUMMARY OUTPUTS =====\n")

summary_rows <- list()

for (k in names(results_all)) {
  res_k <- results_all[[k]]
  
  for (th_name in names(res_k$candidates)) {
    entry <- res_k$candidates[[th_name]]
    diag <- entry$diag
    score <- entry$score$score
    
    ### MODIFIED: Include both absolute and relative fractions in summary
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      timescale = as.integer(k),
      threshold = as.numeric(th_name),
      mean_retained_abs = diag$mean_retained_abs,    # NEW
      mean_retained_rel = diag$mean_retained_rel,    # NEW
      min_retained = min(diag$retained_frac_abs, na.rm = TRUE),
      max_retained = max(diag$retained_frac_abs, na.rm = TRUE),
      zero_frac_median = diag$zero_fraction_median,
      median_nonzero_years = diag$median_nonzero_years,
      stability_pct = diag$stability_pct,
      score = score,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df %>% arrange(timescale, desc(score))

csv_file <- file.path(out_dir, "scf_timescale_threshold_summary.csv")
write.csv(summary_df, csv_file, row.names = FALSE)
cat(sprintf("→ Saved summary: %s\n", basename(csv_file)))

# ==============================================================================
# PLOT: SCORE BY THRESHOLD FOR EACH TIMESCALE
# ==============================================================================

p1 <- ggplot(summary_df, aes(x = factor(threshold), y = score, 
                             fill = factor(timescale))) +
  geom_col(position = position_dodge()) +
  facet_wrap(~ timescale, scales = "free_x") +
  labs(
    title = "SCF Threshold Candidate Scores by Timescale",
    subtitle = "Higher scores indicate better snow persistence characteristics",
    x = "SCF Threshold",
    y = "Score (0-1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(out_dir, "scf_threshold_scores.png"), p1, 
       width = 10, height = 6, dpi = 300)
cat("→ Saved plot: scf_threshold_scores.png\n")

# ==============================================================================
# PLOT: DIAGNOSTIC METRICS COMPARISON
# ==============================================================================

p2 <- ggplot(summary_df, aes(x = factor(threshold), y = median_nonzero_years)) +
  geom_col(aes(fill = factor(timescale)), position = position_dodge()) +
  facet_wrap(~ timescale, scales = "free_x") +
  geom_hline(yintercept = targets$nonzero_years, linetype = "dashed", 
             color = "red", alpha = 0.7) +
  labs(
    title = "Snow Persistence: Median Nonzero Years",
    subtitle = "Red line shows target value",
    x = "SCF Threshold",
    y = "Median Nonzero Years"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(out_dir, "scf_diagnostic_nonzero_years.png"), p2, 
       width = 10, height = 6, dpi = 300)
cat("→ Saved plot: scf_diagnostic_nonzero_years.png\n")

# ==============================================================================
# PLOT: RETAINED FRACTION BY TIMESCALE (showing absolute fraction)
# ==============================================================================

p3 <- ggplot(summary_df, aes(x = factor(threshold), y = mean_retained_abs)) +  # USE ABSOLUTE
  geom_col(aes(fill = factor(timescale)), position = position_dodge()) +
  facet_wrap(~ timescale, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Mean Retained Pixel Fraction (Absolute)",
    subtitle = "Absolute: retained pixels / total domain pixels",
    x = "SCF Threshold",
    y = "Retained Fraction"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(out_dir, "scf_retained_fraction.png"), p3, 
       width = 10, height = 6, dpi = 300)
cat("→ Saved plot: scf_retained_fraction.png\n")

# ==============================================================================
# SAVE BEST MASKS (PER CENTER MONTH) FOR CHOSEN THRESHOLDS
# ==============================================================================

cat("\n===== SAVING MASK NETCDF FILES =====\n")

for (k in names(results_all)) {
  best_th <- results_all[[k]]$best_threshold
  diag <- results_all[[k]]$candidates[[as.character(best_th)]]$diag
  
  # Create a raster stack of 12 masks (one per center month)
  masks <- lapply(diag$mask_month, function(m) {
    r <- rast(scf_stack[[1]])
    vals <- rep(NA_real_, ncell(r))
    vals[domain_mask] <- as.numeric(m[domain_mask])  # both sides same length
    values(r) <- vals
    r
  })
  
  mask_stack <- rast(masks)
  names(mask_stack) <- paste0("mask_k", k, "_m", sprintf("%02d", 1:12))
  
  out_nc <- file.path(out_dir, sprintf("scf_mask_k%s_best_%.3f.nc", k, best_th))
  writeCDF(mask_stack, out_nc, 
           varname = "scf_mask",
           longname = sprintf("SCF mask k=%s best threshold %.3f", k, best_th),
           unit = "logical",
           overwrite = TRUE)
  
  cat(sprintf("→ Saved mask for k=%s (threshold %.3f): %s\n", 
              k, best_th, basename(out_nc)))
}

# ==============================================================================
# FINAL RECOMMENDATION TABLE (with both metrics)
# ==============================================================================

cat("\n===== FINAL RECOMMENDATIONS =====\n")

### MODIFIED: Include both absolute and relative fractions
final_table <- data.frame(
  timescale = integer(),
  chosen_threshold = numeric(),
  mean_retained_abs = numeric(),  # NEW
  mean_retained_rel = numeric(),  # NEW
  zero_frac_median = numeric(),
  median_nonzero_years = numeric(),
  stability_pct = numeric(),
  stringsAsFactors = FALSE
)

for (k in names(results_all)) {
  best_th <- results_all[[k]]$best_threshold
  diag <- results_all[[k]]$candidates[[as.character(best_th)]]$diag
  
  final_table <- rbind(final_table, data.frame(
    timescale = as.integer(k),
    chosen_threshold = best_th,
    mean_retained_abs = diag$mean_retained_abs,  # NEW
    mean_retained_rel = diag$mean_retained_rel,  # NEW
    zero_frac_median = diag$zero_fraction_median,
    median_nonzero_years = diag$median_nonzero_years,
    stability_pct = diag$stability_pct
  ))
}

final_csv <- file.path(out_dir, "scf_timescale_final_recommendations.csv")
write.csv(final_table, final_csv, row.names = FALSE)
cat(sprintf("→ Saved final recommendations: %s\n", basename(final_csv)))

# Print to console
cat("\nFinal recommended thresholds per timescale:\n")
print(final_table)

# ==============================================================================
# TEXT SUMMARY REPORT
# ==============================================================================

summary_txt <- file.path(out_dir, "scf_decision_summary.txt")
sink(summary_txt)

cat("SCF TIMESCALE DECISION ALGORITHM - SUMMARY REPORT\n")
cat("==================================================\n\n")
cat(sprintf("Analysis date: %s\n\n", Sys.time()))

cat("DATA SOURCES:\n")
cat(sprintf("  SCF: %s\n", scf_file))
cat(sprintf("  SWE: %s\n", swe_file))
cat(sprintf("  Period: %s to %s\n\n", min(dates_scf), max(dates_scf)))

cat("SNOW-CAPABLE DOMAIN:\n")
cat(sprintf("  Defined as pixels with ≥5%% long-term mean SCF\n"))
cat(sprintf("  Represents %.1f%% of total domain pixels\n\n", 
            100 * snow_capable_pixels / total_pixels))

cat("METRICS REPORTED:\n")
cat("  Absolute retained fraction: retained pixels / total domain pixels\n")
cat("  Relative retained fraction: retained pixels / snow-capable pixels (≥5%% SCF)\n")
cat("  → Scoring uses ABSOLUTE fraction to preserve spatial context\n")
cat("  → RELATIVE fraction reported for transparency only\n\n")

cat("TIE-BREAKING RULE:\n")
cat("  When multiple thresholds yield identical scores (within 0.001 tolerance),\n")
cat("  the SMALLER threshold is selected to ensure conservative inclusion of\n")
cat("  hydrologically relevant marginal snow zones.\n\n")

cat("SCORING TARGETS:\n")
cat(sprintf("  Median nonzero years: %d\n", targets$nonzero_years))
cat(sprintf("  Zero fraction: %.2f\n", targets$zero_frac))
cat(sprintf("  Stability percentage: %.0f%%\n\n", targets$stability_pct))

cat("FINAL RECOMMENDATIONS:\n")
cat("======================\n\n")

for (i in 1:nrow(final_table)) {
  row <- final_table[i, ]
  cat(sprintf("Timescale: %d-month\n", row$timescale))
  cat(sprintf("  Chosen threshold: %.3f\n", row$chosen_threshold))
  cat(sprintf("  Absolute retained: %.1f%% of total domain\n", 100 * row$mean_retained_abs))
  cat(sprintf("  Relative retained: %.1f%% of snow-capable pixels\n", 100 * row$mean_retained_rel))
  cat(sprintf("  Zero fraction (median): %.3f\n", row$zero_frac_median))
  cat(sprintf("  Median nonzero years: %.1f\n", row$median_nonzero_years))
  cat(sprintf("  Stability: %.1f%%\n\n", row$stability_pct))
}

cat("INTERPRETATION:\n")
cat("===============\n")
cat("✓ Absolute fraction reflects spatial extent of snow relevance across ENTIRE basin\n")
cat("✓ Relative fraction shows efficiency within climatologically relevant snow zones\n")
cat("✓ Smaller thresholds preferred when scores tie (conservative inclusion principle)\n\n")

cat("FILES GENERATED:\n")
cat("================\n")
cat("  1. scf_timescale_threshold_summary.csv - All candidates (with both metrics)\n")
cat("  2. scf_timescale_final_recommendations.csv - Best choices (with both metrics)\n")
cat("  3. scf_threshold_scores.png - Score comparison\n")
cat("  4. scf_diagnostic_nonzero_years.png - Persistence diagnostic\n")
cat("  5. scf_retained_fraction.png - Absolute retained fraction\n")
cat("  6. scf_mask_k*.nc - NetCDF masks for each timescale\n")

sink()
cat(sprintf("→ Saved text summary: %s\n", basename(summary_txt)))

# ==============================================================================
# COMPLETION MESSAGE
# ==============================================================================

cat("\n============================================================\n")
cat("SCF TIMESCALE DECISION ALGORITHM COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\nReview the following files:\n")
cat("  → scf_timescale_final_recommendations.csv - Recommended thresholds (with absolute/relative metrics)\n")
cat("  → scf_threshold_scores.png - Visual comparison\n")
cat("  → scf_decision_summary.txt - Full report with dual metrics\n")
