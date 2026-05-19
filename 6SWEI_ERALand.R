##############################################
# SWEI / SSPI-1 CALCULATION - PARALLEL VERSION (SPI-STYLE OUTPUTS)
# Uses parallel::parLapply - Works on Windows and Linux
# Output naming matches SPI_ERALand.R convention
#
# METHOD (user-selected at runtime):
#   1 = Self-calibration via 6preq_diag_SCF.R (SWEI)
#       Data-driven threshold selection per timescale/basin.
#   2 = Huning & AghaKouchak (2020) fixed 5% SCF threshold (SWEI)
#       Right-aligned 3-month climatology, matches original paper.
#   3 = SSPI-1: Standardized SnowPack Index (monthly, no SCF mask)
#       Zero-inflated gamma distribution, following Stagge et al. (2015) and JRC EDO (2020).
#       Output to: sspi_results_monthly/
#
# Other methodological alignments with Huning & AghaKouchak (2020):
#   - Zero-SWE perturbation: random uniform in (0, min_nonzero)  [methods 1&2]
#   - Minimum data: >= 75% of years per calendar month must be nonzero [methods 1&2]
#   - No clipping of SWEI values  [methods 1&2]
##############################################
# ---- Libraries ----
library(terra)
library(ncdf4)
library(zoo)
library(writexl)
library(parallel)
library(lmomco)   # required for option 3 (SSPI-1 zero-inflated gamma)

# ---- USER SELECTION: SCF DOMAIN METHOD (MOVED TO TOP) ----
cat("\n============================================================\n")
cat("SCF DOMAIN DEFINITION: SELECT METHOD\n")
cat("============================================================\n")
cat("  1 = Self-calibration  — data-driven threshold via 6preq_diag_SCF.R\n")
cat("                          (optimal threshold selected per basin/timescale)\n")
cat("  2 = Huning & AghaKouchak (2020) — fixed 5% SCF threshold\n")
cat("                          (standard approach from the original global paper)\n")
cat("  3 = SSPI-1 (Standardized SnowPack Index, monthly)\n")
cat("                          (zero-inflated gamma, 1981-2020 reference, no SCF mask)\n\n")

scf_method <- " "
while (!scf_method %in% c("1", "2", "3")) {
  scf_method <- trimws(readline(prompt = "Enter 1, 2, or 3: "))
  if (!scf_method %in% c("1", "2", "3")) cat("  Invalid input. Please enter 1, 2, or 3.\n")
}
if (scf_method == "1") {
  cat("\n→ Method selected: Self-calibration (6preq_diag_SCF.R)\n")
} else if (scf_method == "2") {
  cat("\n→ Method selected: Huning & AghaKouchak (2020) fixed 5% SCF threshold\n")
} else {
  cat("\n→ Method selected: SSPI-1 (zero-inflated gamma, monthly, 1981-2020 reference)\n")
}

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")
out_dir <- if (scf_method == "3") "sspi_results_monthly" else "swei_results_seasonal"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ---- Load Basin Boundary ----
basin_path <- "Spatial/nechakoBound_dissolve.kmz"
basin <- if (file.exists(basin_path)) {
  tmp <- tempfile(); dir.create(tmp, showWarnings = FALSE)
  utils::unzip(basin_path, exdir = tmp)
  kml <- list.files(tmp, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)[1]
  v <- vect(kml); if (nrow(v) > 1L) v <- aggregate(v)
  unlink(tmp, recursive = TRUE); v
} else NULL
if (!is.null(basin)) {
  cat("✓ Basin boundary loaded\n")
} else {
  stop("Basin boundary not found: ", basin_path)
}
# ---- Input files ----
swe_file <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
scf_file <- "monthly_data_direct/snow_cover_monthly.nc"
if (!file.exists(swe_file)) stop("SWE file not found: ", swe_file)
if (scf_method %in% c("1", "2") && !file.exists(scf_file))
  stop("SCF file not found (required for methods 1 & 2): ", scf_file)

cat("\n============================================================\n")
cat(if (scf_method == "3") "SSPI-1 CALCULATION (PARALLEL MODE)\n"
    else "SEASONAL SWEI CALCULATION (PARALLEL MODE)\n")
cat("============================================================\n\n")

cat("===== READING NETCDF FILES =====\n")
swe <- rast(swe_file)
if (scf_method %in% c("1", "2")) {
  scf_stack <- rast(scf_file)
} else {
  scf_stack <- NULL  # not needed for SSPI-1
}

# ---- Time extraction ----
extract_time_dimension <- function(raster_obj, file_path) {
  t <- time(raster_obj)
  if (!all(is.na(t))) return(as.Date(t))
  nc <- nc_open(file_path)
  if ("valid_time" %in% names(nc$var)) {
    tv <- ncvar_get(nc, "valid_time")
    nc_close(nc)
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  if ("time" %in% names(nc$dim)) {
    tv <- ncvar_get(nc, "time")
    nc_close(nc)
    return(as.Date(as.POSIXct(tv, origin = "1970-01-01", tz = "UTC")))
  }
  nc_close(nc)
  n <- nlyr(raster_obj)
  start <- as.Date("1950-01-01")
  seq(start, by = "month", length.out = n)
}
dates_swe <- extract_time_dimension(swe, swe_file)
if (!is.null(scf_stack)) {
  dates_scf <- extract_time_dimension(scf_stack, scf_file)
  terra::time(scf_stack) <- dates_scf
} else {
  dates_scf <- NULL
}
terra::time(swe) <- dates_swe
dates <- dates_swe

# ---- Reproject to common grid ----
if (!is.null(scf_stack) && !compareGeom(scf_stack[[1]], swe[[1]], stopOnError = FALSE)) {
  cat("✓ Resampling SWE to match SCF grid...\n")
  swe <- resample(swe, scf_stack, method = "bilinear")
}
target_crs <- if (!is.null(scf_stack)) crs(scf_stack) else crs(swe)
if (!is.null(basin) && !same.crs(basin, target_crs)) {
  basin <- project(basin, target_crs)
}

# ---- Basin masking ----
cat("\n===== MASKING SWE TO BASIN BOUNDARY =====\n")
swe <- mask(swe, basin, inverse = FALSE, touches = TRUE)
basin_pixels <- global(swe[[1]], "notNA")$notNA
total_pixels <- ncell(swe)
cat(sprintf("✓ Basin masking complete: %d pixels (%.1f%% of raster)\n",
            basin_pixels, 100 * basin_pixels / total_pixels))

# ---- Convert SWE to mm if needed ----
swe_mean_sample <- global(swe[[1]], "mean", na.rm = TRUE)$mean
if (!is.na(swe_mean_sample) && swe_mean_sample < 10) {
  swe <- swe * 1000
  cat("✓ Converted SWE from meters to mm\n")
}

# ---- Convert SCF to fraction if needed ----
if (!is.null(scf_stack)) {
  scf_max <- global(scf_stack[[1]], "max", na.rm = TRUE)$max
  if (!is.na(scf_max) && scf_max > 1.5) {
    scf_stack <- scf_stack / 100
    cat("✓ Converted SCF from percent to fraction\n")
  }
}

# ---- Prepare matrices ----
swe_matrix <- values(swe, mat = TRUE)
dates <- dates_swe
months_all <- as.integer(format(dates, "%m"))
years_all <- as.integer(format(dates, "%Y"))
n_pixels <- nrow(swe_matrix)
cat(sprintf("Processing %d basin pixels...\n", n_pixels))

# ==============================================================================
#   STEP 1: Apply 3-month moving average to RAW SWE values (PREPROCESSING)
#   NOTE: used only by methods 1 & 2 (SWEI). Method 3 (SSPI-1) operates on
#   raw monthly SWE values directly — no temporal smoothing is applied.
# ==============================================================================
if (scf_method %in% c("1", "2")) {
  swe_smoothed <- matrix(NA, nrow = nrow(swe_matrix), ncol = ncol(swe_matrix))
  for (i in 1:nrow(swe_matrix)) {
    swe_smoothed[i, ] <- zoo::rollapply(swe_matrix[i, ],
                                        width = 3,
                                        FUN = sum,
                                        align = "right",
                                        fill = NA,
                                        na.rm = FALSE)   # NA if any of the three months is NA
  }
  cat("✓ 3-month rolling sum applied to SWE (SWEI preprocessing)\n")
} else {
  swe_smoothed <- NULL  # not used by SSPI-1
  cat("✓ No temporal smoothing applied (SSPI-1 uses raw monthly SWE)\n")
}

# ==============================================================================
#   USER CHOICE: SCF DOMAIN DEFINITION METHOD (confirmed above — echoing selection)
# ==============================================================================
cat(sprintf("\n→ Proceeding with Method %s: %s\n", scf_method, scf_method_label))

# ==============================================================================
#   STEP 2: BUILD SCF MASK  (branched by user choice)
# ==============================================================================
cat("\n===== STEP 2: Building SCF mask =====\n")
target_k    <- 3
scf_method_label <- if (scf_method == "1") "Self-calibration" else if (scf_method == "2") "Huning & AghaKouchak (2020)" else "SSPI-1 (zero-inflated gamma, no SCF mask)"

if (scf_method == "3") {
  
  # ============================================================================
  #   BRANCH 3 — SSPI-1 (Zero-inflated Gamma, Monthly, No SCF Mask)
  #   Skips SCF masking entirely; uses calendar-month gamma fitting.
  #   Reference period: 1981-2020 (per Stagge et al., 2015 / SSPI convention).
  # ============================================================================
  cat("\n→ SSPI-1 selected: skipping SCF mask build (not applicable).\n")
  cat("  Zero-inflated gamma fitting will be applied per calendar month.\n")
  scf_mask_list  <- NULL
  basin_scf_mask <- NULL
  
  # Reference period for SSPI-1
  sspi_ref_start <- min(dates) #as.Date("1981-01-01")
  sspi_ref_end   <- max(dates) #as.Date("2020-12-31")
  ref_idx        <- which(dates >= sspi_ref_start & dates <= sspi_ref_end)
  if (length(ref_idx) < 400) {
    warning(sprintf("SSPI-1 reference period has only %d months (recommended: ≥480)", length(ref_idx)))
  }
  cat(sprintf("✓ SSPI-1 reference period: %s to %s (%d months)\n",
              dates[ref_idx[1]], dates[ref_idx[length(ref_idx)]], length(ref_idx)))
  
  # Zero-inflation screening: exclude pixels with >50% zero SWE in reference period
  cat("\n===== ZERO-INFLATION SCREENING (SSPI-1) =====\n")
  ref_swe_sspi   <- swe_matrix[, ref_idx, drop = FALSE]
  zero_prop_sspi <- rowMeans(ref_swe_sspi <= 0, na.rm = TRUE)
  invalid_pix_sspi <- which(zero_prop_sspi > 0.50)
  valid_pix_sspi   <- which(zero_prop_sspi <= 0.50)
  cat(sprintf("✓ Valid pixels (≤50%% zero SWE in ref period): %d (%.1f%%)\n",
              length(valid_pix_sspi), 100 * length(valid_pix_sspi) / n_pixels))
  if (length(invalid_pix_sspi) > 0) {
    cat(sprintf("  Excluded pixels (>50%% zero SWE): %d\n", length(invalid_pix_sspi)))
    swe_matrix[invalid_pix_sspi, ] <- NA
  }
  
} else if (scf_method == "1") {
  
  # ============================================================================
  #   BRANCH 1 — SELF-CALIBRATION  (6preq_diag_SCF.R)
  # ============================================================================
  scf_diag_dir    <- "scf_timescale_diagnostics"
  scf_diag_csv    <- file.path(scf_diag_dir, "scf_timescale_final_recommendations.csv")
  scf_diag_script <- "6preq_diag_SCF.R"   # must exist in the working directory
  
  # Run the diagnostic script automatically if its outputs are missing
  if (!file.exists(scf_diag_csv)) {
    cat("\n  SCF diagnostic outputs not found in '", scf_diag_dir, "'\n", sep = "")
    if (file.exists(scf_diag_script)) {
      run_diag <- ""
      while (!run_diag %in% c("y", "n")) {
        run_diag <- tolower(trimws(readline(
          prompt = "  Run 6preq_diag_SCF.R now to generate them? (y/n): ")))
      }
      if (run_diag == "y") {
        cat("  Sourcing 6preq_diag_SCF.R ...\n")
        source(scf_diag_script, local = FALSE)
        cat("  ✓ 6preq_diag_SCF.R complete\n")
      } else {
        stop("Self-calibration outputs are required. ",
             "Run 6preq_diag_SCF.R first, or re-run and choose method 2.")
      }
    } else {
      stop("'", scf_diag_script, "' not found in working directory and diagnostic ",
           "outputs are missing.\n",
           "  Place 6preq_diag_SCF.R in '", getwd(), "' and re-run, or choose method 2.")
    }
  }
  
  # Read data-driven threshold chosen for k=3
  scf_threshold <- 0.10   # fallback default
  diag_df <- read.csv(scf_diag_csv, stringsAsFactors = FALSE)
  if ("timescale" %in% names(diag_df) && "chosen_threshold" %in% names(diag_df)) {
    k3_row <- diag_df[as.character(diag_df$timescale) == as.character(target_k), ]
    if (nrow(k3_row) > 0) scf_threshold <- as.numeric(k3_row$chosen_threshold[1])
  }
  
  # Find the closest-matching mask NetCDF for that threshold
  mask_file_pattern <- sprintf("scf_mask_k%d_*.nc", target_k)
  mask_candidates   <- list.files(scf_diag_dir, pattern = glob2rx(mask_file_pattern),
                                  full.names = TRUE)
  if (length(mask_candidates) == 0)
    stop(sprintf("No SCF mask NetCDF found for k=%d in '%s'", target_k, scf_diag_dir))
  
  best_mask <- mask_candidates[[1]]
  if (length(mask_candidates) > 1) {
    min_diff <- Inf
    for (f in mask_candidates) {
      m_match <- regexec("threshold_([0-9.]+)\\.nc$", basename(f))
      if (length(m_match[[1]]) > 1) {
        th_in_file <- as.numeric(regmatches(basename(f), m_match)[[1]][2])
        d <- abs(th_in_file - scf_threshold)
        if (d < min_diff) { min_diff <- d; best_mask <- f }
      }
    }
  }
  cat(sprintf("✓ Self-calibration mask loaded: '%s' (threshold = %.3f)\n",
              basename(best_mask), scf_threshold))
  
  scf_mask_stack <- rast(best_mask)
  scf_mask_list  <- vector("list", 12)
  for (m in 1:12) {
    if (m <= nlyr(scf_mask_stack)) {
      mask_vals      <- values(scf_mask_stack[[m]])
      scf_mask_list[[m]] <- !is.na(mask_vals) & (mask_vals == 1)
    } else {
      scf_mask_list[[m]] <- rep(FALSE, ncell(swe))
    }
  }
  
  # Basin-level SCF mask: month included when ≥50% of basin pixels are flagged
  cat("\n===== COMPUTING BASIN-LEVEL SCF MASK (self-calibration) =====\n")
  basin_scf_mask <- logical(12)
  for (m in 1:12) {
    if (m <= nlyr(scf_mask_stack)) {
      mv       <- as.vector(values(scf_mask_stack[[m]]))
      mv_valid <- mv[!is.na(mv)]
      basin_scf_mask[m] <- length(mv_valid) > 0 && mean(mv_valid, na.rm = TRUE) >= 0.5
    }
  }
  
} else {
  
  # ============================================================================
  #   BRANCH 2 — HUNING & AGHAKOUCHAK (2020)  fixed 5% SCF threshold
  # ============================================================================
  scf_threshold <- 0.05
  cat(sprintf("✓ H&A (2020): building SCF mask with fixed %.0f%% threshold\n",
              100 * scf_threshold))
  
  # Mask SCF stack to the basin so pixel indices align exactly with SWE
  scf_basin      <- mask(scf_stack, basin, inverse = FALSE, touches = TRUE)
  scf_matrix_ha  <- values(scf_basin, mat = TRUE)   # ncell(swe) × nlyr
  months_scf_ha  <- as.integer(format(dates_scf, "%m"))
  
  # Right-aligned 3-month mean SCF climatology per calendar month
  # For month m the window is (m-2, m-1, m) — mirrors the SWE aggregation window.
  # Modular arithmetic wraps correctly: e.g. m=1 → window months 11, 12, 1
  scf_mask_list <- vector("list", 12)
  for (m in 1:12) {
    win_months        <- ((c(m - 3L, m - 2L, m - 1L) + 120L) %% 12L) + 1L
    idx_win           <- which(months_scf_ha %in% win_months)
    if (length(idx_win) == 0) {
      scf_mask_list[[m]] <- rep(FALSE, ncell(swe))
    } else {
      scf_clim           <- rowMeans(scf_matrix_ha[, idx_win, drop = FALSE], na.rm = TRUE)
      scf_mask_list[[m]] <- !is.na(scf_clim) & (scf_clim >= scf_threshold)
    }
    cat(sprintf("  Month %2d (%s): %d pixels pass SCF >= %.0f%%\n",
                m, month.abb[m],
                sum(scf_mask_list[[m]], na.rm = TRUE),
                100 * scf_threshold))
  }
  rm(scf_basin, scf_matrix_ha)   # free memory
  
  # Basin-level SCF mask: month included when ≥50% of basin pixels are flagged
  cat("\n===== COMPUTING BASIN-LEVEL SCF MASK (H&A 2020) =====\n")
  basin_scf_mask <- logical(12)
  for (m in 1:12) {
    mv             <- scf_mask_list[[m]]
    basin_scf_mask[m] <- length(mv) > 0 && mean(mv, na.rm = TRUE) >= 0.5
  }
}

if (scf_method %in% c("1", "2")) {
  cat(sprintf("✓ Basin SCF mask: months included = %s\n",
              paste(month.abb[which(basin_scf_mask)], collapse = ", ")))
}

# ==============================================================================
#   BASIN-AVERAGED SWE SETUP
#   Correct order (matches SPI/SPEI approach):
#     1. Average RAW SWE across basin pixels   → basin time series
#     2. Apply 3-month rolling average          → smoothed basin series
#   (The old approach averaged the already-smoothed per-pixel matrix, which
#    gives the same numbers but couples the spatial and temporal smoothing
#    steps and makes the basin series dependent on per-pixel NA patterns.)
# ==============================================================================
cat("\n===== COMPUTING BASIN-AVERAGED SWE (raw → average → smooth) =====\n")
if (scf_method %in% c("1", "2")) {
  swe_basin_avg_raw      <- colMeans(swe_matrix, na.rm = TRUE)
  swe_basin_avg_smoothed  <- zoo::rollapply(swe_basin_avg_raw,
                                            width = 3,
                                            FUN   = sum,
                                            align = "right",
                                            fill  = NA,
                                            na.rm = FALSE)
  cat(sprintf("✓ Basin-averaged smoothed SWE: %d time steps, %.1f%% non-NA\n",
              length(swe_basin_avg_smoothed),
              100 * mean(!is.na(swe_basin_avg_smoothed))))
} else {
  swe_basin_avg_raw      <- NULL
  swe_basin_avg_smoothed <- NULL
  cat("✓ Basin-averaged smoothed SWE: skipped (not used by SSPI-1)\n")
}
basin_avg_swei_results <- list()

# ==============================================================================
#   Helper Functions
# ==============================================================================
clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

# ==============================================================================
#   SEASONAL SWEI FUNCTION (Gringorten)
# ==============================================================================
gringorten_swei_seasonal <- function(v_smoothed, dates_vec, scf_mask_list,
                                     cell_index = NULL,
                                     basin_scf_mask = NULL,
                                     eps = 1e-6) {
  # cell_index    : integer  → per-pixel mask lookup in scf_mask_list
  # basin_scf_mask: logical(12) → basin-level mask, one value per calendar month
  # If both are NULL the mask step is skipped entirely (not recommended)
  if (length(v_smoothed) == 0 || all(is.na(v_smoothed)) || is.null(v_smoothed)) {
    return(list(swei = rep(NA_real_, length(dates_vec)), method = rep(NA_integer_, length(dates_vec))))
  }
  v_clean <- v_smoothed
  v_clean[!is.finite(v_clean)] <- NA_real_
  n <- length(v_clean)
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n) # 1=Gringorten, 2=NA/Masked
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    # Huning & AghaKouchak (2020): only process a pixel-month if at least 75%
    # of years for that calendar month have a nonzero (3-month) SWE value
    all_m_idx <- which(mon == m)
    if (length(all_m_idx) == 0) next
    nonzero_frac <- sum(v_clean[all_m_idx] > 0, na.rm = TRUE) / length(all_m_idx)
    if (nonzero_frac < 0.75) next
    
    idx <- which(mon == m & is.finite(v_clean))
    if (length(idx) == 0) next
    
    if (is.null(scf_mask_list) || length(scf_mask_list) < m) next
    # Check mask for this pixel (if cell_index provided)
    if (!is.null(cell_index)) {
      # Per-pixel path: check whether this cell is flagged as snowy for month m
      mask_vec <- scf_mask_list[[m]]
      if (length(mask_vec) < cell_index || !mask_vec[cell_index]) next
    } else if (!is.null(basin_scf_mask)) {
      # Basin-averaged path: use the pre-computed basin-level mask
      if (length(basin_scf_mask) < m || !basin_scf_mask[m]) next
    }
    # If both are NULL: no mask applied (original fallback; not recommended)
    
    samp <- v_clean[idx]
    # Huning & AghaKouchak (2020): randomly perturb zeros with positive numbers
    # smaller than the minimum nonzero value in the sample
    zero_idx <- which(samp == 0 & !is.na(samp))
    if (length(zero_idx) > 0) {
      nonzero_vals <- samp[samp > 0 & !is.na(samp)]
      if (length(nonzero_vals) > 0) {
        min_nonzero <- min(nonzero_vals)
        samp[zero_idx] <- runif(length(zero_idx), min = 0, max = min_nonzero)
      }
    }
    valid_idx <- which(!is.na(samp))
    if (length(valid_idx) < 3) next
    
    samp_valid <- samp[valid_idx]
    idx_valid <- idx[valid_idx]
    
    r <- rank(samp_valid, ties.method = "average")
    N <- length(samp_valid)
    p_val <- (r - 0.44) / (N + 0.12)
    p_val <- clip_prob(p_val, eps = eps)
    z_val <- qnorm(p_val)
    
    z[idx_valid] <- z_val
    method_used[idx_valid] <- 1
  }
  list(swei = z, method = method_used)
}

# ==============================================================================
#   SSPI-1 FUNCTION — Zero-Inflated Gamma (Monthly, Calendar-Month Stratified)
#   Used only when scf_method == "3"
# ==============================================================================
monthly_sspi <- function(v, dates_vec, ref_idx, eps = 1e-6) {
  # v        : monthly SWE (mm) for one pixel
  # dates_vec: corresponding dates
  # ref_idx  : indices of reference period within dates_vec
  n <- length(v)
  if (n == 0 || all(is.na(v)) || is.null(v)) {
    return(list(sspi = rep(NA_real_, n), method = rep(NA_integer_, n)))
  }
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n)  # 1=Gamma, 2=Empirical, 3=Zero-var, 4=Excluded
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    ref_mon_idx <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
    ref_samp <- v_clean[ref_mon_idx]
    ref_samp <- ref_samp[is.finite(ref_samp)]
    obs_mon_idx <- which(mon == m)
    
    if (length(ref_samp) < 10) {
      if (length(obs_mon_idx) >= 3) {
        obs_vals <- v_clean[obs_mon_idx]
        fin_idx  <- is.finite(obs_vals)
        n_fin    <- sum(fin_idx)
        if (n_fin >= 3) {
          p <- (rank(obs_vals[fin_idx], ties.method = "average") - 0.44) / (n_fin + 0.12)
          p <- clip_prob(p, eps = eps)
          z[obs_mon_idx[fin_idx]] <- qnorm(p)
          method_used[obs_mon_idx[fin_idx]] <- 2
        }
      }
      next
    }
    ref_var <- var(ref_samp, na.rm = TRUE)
    if (!is.finite(ref_var) || ref_var < .Machine$double.eps) {
      z[obs_mon_idx] <- 0
      method_used[obs_mon_idx] <- 3
      next
    }
    if (ref_var < 0.01) {
      obs_vals <- v_clean[obs_mon_idx]
      fin_idx <- is.finite(obs_vals)
      if (sum(fin_idx) >= 3) {
        p <- (rank(obs_vals[fin_idx], ties.method = "average") - 0.44) /
          (sum(fin_idx) + 0.12)
        p <- clip_prob(p, eps = eps)
        z[obs_mon_idx[fin_idx]] <- qnorm(p)
        method_used[obs_mon_idx[fin_idx]] <- 2
      }
      next
    }
    p0 <- mean(ref_samp <= 0, na.rm = TRUE)
    ref_pos <- ref_samp[ref_samp > 0]
    if (length(ref_pos) >= 10) {
      lm_obj <- try(lmomco::lmoms(ref_pos), silent = TRUE)
      if (!inherits(lm_obj, "try-error") && !is.null(lm_obj)) {
        par_obj <- try(lmomco::pargam(lm_obj), silent = TRUE)
        if (!inherits(par_obj, "try-error") && !is.null(par_obj)) {
          obs_vals <- v_clean[obs_mon_idx]
          p_m <- rep(NA_real_, length(obs_vals))
          pos_idx <- which(is.finite(obs_vals) & obs_vals > 0)
          if (length(pos_idx) > 0) {
            Fg <- try(lmomco::cdfgam(obs_vals[pos_idx], par_obj), silent = TRUE)
            if (!inherits(Fg, "try-error")) p_m[pos_idx] <- p0 + (1 - p0) * as.numeric(Fg)
          }
          zero_idx2 <- which(is.finite(obs_vals) & obs_vals <= 0)
          if (length(zero_idx2) > 0) p_m[zero_idx2] <- p0
          p_m <- clip_prob(p_m, eps = eps)
          z[obs_mon_idx] <- qnorm(p_m)
          method_used[obs_mon_idx] <- 1
          next
        }
      }
    }
    all_samp <- c(ref_samp, v_clean[obs_mon_idx])
    all_samp <- all_samp[is.finite(all_samp)]
    obs_vals <- v_clean[obs_mon_idx]
    fin_obs <- is.finite(obs_vals)
    if (length(all_samp) >= 10 && sum(fin_obs) > 0) {
      combined <- c(all_samp, obs_vals[fin_obs])
      r_all <- rank(combined, ties.method = "average")
      r_obs <- r_all[(length(all_samp) + 1):(length(all_samp) + sum(fin_obs))]
      p <- (r_obs - 0.44) / (length(combined) + 0.12)
      p <- clip_prob(p, eps = eps)
      z[obs_mon_idx[fin_obs]] <- qnorm(p)
      method_used[obs_mon_idx[fin_obs]] <- 2
    }
  }
  fin_z <- is.finite(z)
  z[fin_z & z < -5] <- -5
  z[fin_z & z >  5] <-  5
  list(sspi = z, method = method_used)
}

# ==============================================================================
#   MAIN CALCULATION LOOP - PARALLEL (Matches SPI Structure)
# ==============================================================================
timescales <- c(1) # SWEI typically 1-month output, but structured like SPI/SPEI/...
# ---- Set up parallel cluster ----
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\n✓ Starting parallel cluster with %d cores\n", n_cores))
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
cl <- makeCluster(n_cores)

if (scf_method == "3") {
  clusterExport(cl, varlist = c("monthly_sspi", "clip_prob",
                                "swe_matrix", "dates", "ref_idx"),
                envir = environment())
  clusterEvalQ(cl, { library(zoo); library(lmomco) })
} else {
  clusterExport(cl, varlist = c("gringorten_swei_seasonal", "clip_prob",
                                "swe_smoothed", "dates", "scf_mask_list",
                                "basin_scf_mask"),
                envir = environment())
  clusterEvalQ(cl, { library(zoo) })
}

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
# ---- Summary file setup ----
summary_combined_file <- file.path(out_dir,
                                   if (scf_method == "3") "sspi_summary.txt" else "swei_all_timescales_summary.txt")
all_summaries <- list()

for (sc in timescales) {
  
  if (scf_method == "3") {
    # ============================================================================
    #   OPTION 3: SSPI-1 CALCULATION (zero-inflated gamma, monthly)
    # ============================================================================
    cat(sprintf("\n===== SSPI-%d (Monthly, Parallel) =====\n", sc))
    cat(sprintf("Processing %d pixels in parallel... ", n_pixels))
    start_time <- Sys.time()
    
    pixel_list <- parLapply(cl, seq_len(n_pixels), function(i) {
      tryCatch(
        monthly_sspi(swe_matrix[i, ], dates, ref_idx, eps = 1e-6),
        error = function(e) list(sspi = rep(NA_real_, length(dates)),
                                 method = rep(NA_integer_, length(dates)))
      )
    })
    
    sspi_indices  <- do.call(rbind, lapply(pixel_list, `[[`, "sspi"))
    method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
    end_time <- Sys.time()
    cat(sprintf("done (%.1f min)\n", as.numeric(difftime(end_time, start_time, units = "mins"))))
    
    if (is.null(dim(sspi_indices)) || nrow(sspi_indices) != nrow(swe_matrix))
      stop("SSPI-1 calculation failed: unexpected result dimensions")
    sspi_indices[!is.finite(sspi_indices)] <- NA_real_
    
    # Diagnostics
    clipped_low  <- sum(sspi_indices < -4.9, na.rm = TRUE)
    clipped_high <- sum(sspi_indices >  4.9, na.rm = TRUE)
    cat(sprintf("  Clipping: %d dry, %d wet\n", clipped_low, clipped_high))
    basin_mask_sspi <- !is.na(swe_matrix[, 1])
    na_rate_basin   <- 100 * mean(is.na(sspi_indices[basin_mask_sspi, ]))
    cat(sprintf("  NA rate (basin): %.3f%%\n", na_rate_basin))
    gamma_pct      <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
    empirical_pct  <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
    zero_var_pct   <- 100 * sum(method_matrix == 3, na.rm = TRUE) / sum(!is.na(method_matrix))
    excluded_pct   <- 100 * sum(method_matrix == 4, na.rm = TRUE) / sum(!is.na(method_matrix))
    cat(sprintf("  Methods: Gamma=%.1f%%, Empirical=%.1f%%, Zero-var=%.1f%%, Excluded=%.1f%%\n",
                gamma_pct, empirical_pct, zero_var_pct, excluded_pct))
    
    all_summaries[[as.character(sc)]] <- list(
      sc = sc, date = Sys.time(), na_rate = na_rate_basin,
      gamma_pct = gamma_pct, empirical_pct = empirical_pct,
      zero_var_pct = zero_var_pct, excluded_pct = excluded_pct,
      sspi_indices = sspi_indices, method_matrix = method_matrix
    )
    
    months_all_sc <- as.numeric(format(dates, "%m"))
    
    # CSV files
    cat("  -> Saving CSV files...\n")
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sub  <- sspi_indices[, idx_m, drop = FALSE]
      df   <- as.data.frame(sub)
      colnames(df) <- format(dates[idx_m], "%Y")
      coords <- xyFromCell(swe[[1]], 1:nrow(sub))
      df <- cbind(lon = coords[,1], lat = coords[,2], df)
      write.csv(df, file.path(out_dir, sprintf("sspi_%02d_month%02d_%s.csv", sc, m, month_names[m])),
                row.names = FALSE)
    }
    
    # Excel workbook
    cat("  -> Saving Excel workbook...\n")
    excel_data <- list()
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sub  <- sspi_indices[, idx_m, drop = FALSE]
      df   <- as.data.frame(sub)
      colnames(df) <- format(dates[idx_m], "%Y")
      coords <- xyFromCell(swe[[1]], 1:nrow(sub))
      excel_data[[month_names[m]]] <- cbind(lon = coords[,1], lat = coords[,2], df)
    }
    xlsx_file <- file.path(out_dir, sprintf("sspi_%02d_all_months.xlsx", sc))
    if (file.exists(xlsx_file)) file.remove(xlsx_file)
    write_xlsx(excel_data, xlsx_file)
    
    # Method distribution map
    cat("  -> Creating distribution map...\n")
    method_mapping <- data.frame(
      code  = c(1, 2, 3, 4),
      name  = c("Gamma", "Empirical", "Zero-Variance", "Excluded"),
      color = c("#4575b4", "#d73027", "#91bfdb", "#cccccc"),
      stringsAsFactors = FALSE
    )
    method_raster <- rast(swe[[1]])
    dist_vals <- method_matrix[, 1]
    dist_vals[is.na(values(swe[[1]]))] <- NA
    values(method_raster) <- dist_vals
    png_file <- file.path(out_dir, sprintf("sspi_%02d_method_map.png", sc))
    png(png_file, width = 1200, height = 800, res = 150)
    plot(method_raster, col = method_mapping$color, breaks = c(0.5, 1.5, 2.5, 3.5, 4.5),
         legend = FALSE, main = sprintf("SSPI-%d: Fitting Method Distribution", sc),
         axes = FALSE, box = FALSE)
    if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
    dist_counts <- table(factor(dist_vals[!is.na(dist_vals)], levels = 1:4,
                                labels = method_mapping$name))
    legend("bottomright", legend = sprintf("%s (%d pixels)", names(dist_counts), dist_counts),
           fill = method_mapping$color, title = "Method", cex = 0.8,
           bg = "white", bty = "o", border = "gray50")
    valid_pix <- sum(!is.na(values(swe[[1]])))
    mtext(sprintf("Valid pixels: %d | Gamma: %.1f%%", valid_pix, gamma_pct),
          side = 1, line = 0.5, cex = 0.7, col = "gray30")
    dev.off()
    
    # NetCDF files
    cat("  -> Saving NetCDF files...\n")
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sspi_rast <- rep(rast(swe[[1]]), length(idx_m))
      values(sspi_rast) <- sspi_indices[, idx_m, drop = FALSE]
      terra::time(sspi_rast) <- dates[idx_m]
      writeCDF(sspi_rast, file.path(out_dir,
                                    sprintf("sspi_%02d_month%02d_%s.nc", sc, m, month_names[m])),
               varname = "sspi",
               longname = sprintf("Standardized SnowPack Index (SSPI-%d, monthly)", sc),
               unit = "standardized_index", missval = -9999, overwrite = TRUE)
    }
    cat(sprintf("✓ SSPI-%d complete\n", sc))
    gc()
    next  # skip SWEI outputs below
  }
  
  # ============================================================================
  #   OPTIONS 1 & 2: SWEI CALCULATION (Gringorten, SCF-masked)
  # ============================================================================
  cat(sprintf("\n===== SWEI-%d (Sequential Calculation) =====\n", sc))
  # ==============================================================================
  #   CALCULATE BASIN-AVERAGED SWEI
  # ==============================================================================
  avg_swei_out <- tryCatch(
    gringorten_swei_seasonal(swe_basin_avg_smoothed, dates, scf_mask_list,
                             cell_index = NULL, basin_scf_mask = basin_scf_mask,
                             eps = 1e-6),
    error = function(e) list(swei = rep(NA_real_, length(dates)), method = rep(NA_integer_, length(dates)))
  )
  basin_avg_swei_results[[as.character(sc)]] <- avg_swei_out$swei
  
  cat(sprintf("Processing %d pixels in parallel... ", n_pixels))
  start_time <- Sys.time()
  # ---- PARALLEL PROCESSING ----
  pixel_list <- parLapply(cl, seq_len(n_pixels), function(i) {
    tryCatch(
      gringorten_swei_seasonal(swe_smoothed[i, ], dates, scf_mask_list, cell_index = i, eps = 1e-6),
      error = function(e) list(swei = rep(NA_real_, length(dates)),
                               method = rep(NA_integer_, length(dates)))
    )
  })
  
  # Convert to matrix
  swei_indices <- do.call(rbind, lapply(pixel_list, `[[`, "swei"))
  method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
  
  # Diagnostics
  # Informational: report any extreme SWEI values (no clipping applied)
  extreme_dry  <- sum(swei_indices < -3.0, na.rm = TRUE)
  extreme_wet  <- sum(swei_indices >  3.0, na.rm = TRUE)
  cat(sprintf("  Extreme values (|SWEI| > 3): %d dry, %d wet (retained, no clipping)\n",
              extreme_dry, extreme_wet))
  basin_mask_swei  <- !is.na(swe_smoothed)
  na_rate_basin  <- 100 * mean(is.na(swei_indices))
  cat(sprintf("  NA rate (basin): %.3f%%\n", na_rate_basin))
  
  # Method distribution
  gringorten_pct  <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
  masked_pct      <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
  
  cat(sprintf("  Methods: Gringorten=%.1f%%, Masked/NA=%.1f%%\n",
              gringorten_pct, masked_pct))
  
  # Store summary
  all_summaries[[as.character(sc)]] <- list(
    sc = sc,
    date = Sys.time(),
    na_rate = na_rate_basin,
    gringorten_pct = gringorten_pct,
    masked_pct = masked_pct,
    swei_indices = swei_indices,
    method_matrix = method_matrix
  )
  
  # ---- OUTPUTS (MATCHING SPI NAMING) ----
  months_all <- as.numeric(format(dates, "%m"))
  # CSV files
  cat("  -> Saving CSV files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    csv_file <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.csv", sc, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE)
  }
  
  # Excel workbook
  cat("  -> Saving Excel workbook...\n")
  excel_data <- list()
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    excel_data[[as.character(m)]] <- df
  }
  xlsx_file <- file.path(out_dir, sprintf("swei_%02d_all_months.xlsx", sc))
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  write_xlsx(excel_data, xlsx_file)
  
  # Distribution map
  cat("  -> Creating distribution map...\n")
  swei_dist_mapping <- data.frame(
    code = c(1, 2),
    name = c("Gringorten", "Masked/NA"),
    color = c("#4575b4", "#91bfdb"),
    stringsAsFactors = FALSE
  )
  swei_dist_raster <- rast(swe[[1]])
  dist_sim <- apply(method_matrix, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_integer_)
    as.integer(names(sort(table(x), decreasing = TRUE))[1])
  })
  dist_sim[is.na(dist_sim)] <- NA
  values(swei_dist_raster) <- dist_sim
  png_file_swei <- file.path(out_dir, sprintf("swei_%02d_distribution_map.png", sc))
  png(png_file_swei, width = 1800, height = 1000, res = 150)
  layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))
  # ── Panel 1: Map ──────────────────────────────────────────
  plot(swei_dist_raster,
       col    = swei_dist_mapping$color,
       breaks = c(0.5, 1.5, 2.5),
       legend = FALSE,
       main   = sprintf("SWEI-%d: Seasonal Method", sc),
       axes   = FALSE, box = FALSE)
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  valid_pixels <- sum(!is.na(values(swe)))
  mtext(sprintf("Valid pixels: %d | Gringorten: %.1f%%", valid_pixels, gringorten_pct),
        side = 1, line = 1, cex = 0.75, col = "gray30")
  # ── Panel 2: Legend ─
  par(mar = c(0, 0, 0, 0))
  plot.new()
  dist_counts <- table(factor(dist_sim, levels = 1:2,
                              labels = swei_dist_mapping$name))
  legend_entries <- sprintf("%s (%d pixels)", names(dist_counts), dist_counts)
  legend("center",
         legend = legend_entries,
         fill   = swei_dist_mapping$color,
         title  = "Method",
         cex    = 0.95,
         bg     = "white",
         bty    = "o",
         border = "gray50")
  dev.off()
  layout(1)
  
  # NetCDF files (Monthly split like SPI)
  cat("  -> Saving NetCDF files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0) next
    swei_rast  <- rep(rast(swe[[1]]), length(idx_m))
    sub_matrix  <- swei_indices[, idx_m, drop = FALSE]
    if (nrow(sub_matrix) != ncell(swei_rast)) {
      stop(sprintf("Dimension mismatch: Matrix %d rows vs Raster %d cells.",
                   nrow(sub_matrix), ncell(swei_rast)))
    }
    values(swei_rast)  <- sub_matrix
    terra::time(swei_rast)  <- dates[idx_m]
    nc_file  <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.nc", sc, m, month_names[m]))
    writeCDF(swei_rast, nc_file,
             varname = "swei",
             longname = sprintf("Seasonal SWEI (SWEI-%d)", sc),
             unit = "standardized_index",
             missval = -9999,
             overwrite = TRUE)
  }
  cat(sprintf("✓ SWEI-%d complete\n", sc))
  gc()
}

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\n✓ Parallel cluster stopped\n")

# ==============================================================================
#   SAVE BASIN-AVERAGED SWEI TO CSV (Wide Format like SPI)
# ==============================================================================
cat("\n===== SAVING BASIN-AVERAGED SWEI =====\n")
months_all  <- as.integer(format(dates, "%m"))
years_all   <- as.integer(format(dates, "%Y"))
month_names_save  <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
if (scf_method != "3") {
  for (sc in timescales) {
    swei_series <- basin_avg_swei_results[[as.character(sc)]]
    yrs <- sort(unique(years_all))
    df_out <- data.frame(Year = yrs)
    for (m in 1:12) {
      idx_m <- which(months_all == m)
      yr_m  <- years_all[idx_m]
      val_m <- swei_series[idx_m]
      col   <- setNames(val_m, yr_m)
      df_out[[month_names_save[m]]] <- col[as.character(yrs)]
    }
    csv_file <- file.path(out_dir,
                          sprintf("swei_%02d_basin_averaged_by_month.csv", sc))
    write.csv(df_out, csv_file, row.names = FALSE, na = "")
    cat(sprintf("✓ Saved basin-averaged SWEI-%d (12 monthly series) to: %s\n", sc, csv_file))
  }
}

# ==============================================================================
#   WRITE SUMMARY FILE AT END (Matches SPI Layout)
# ==============================================================================
cat("\n===== WRITING SUMMARY FILE =====\n")

if (scf_method == "3") {
  # SSPI-1 summary
  cat(
    "============================================================\n",
    "MONTHLY SSPI-1 SUMMARY - NECHAKO BASIN\n",
    "ERA5-Land Snow Water Equivalent\n",
    "Method 3: Zero-Inflated Gamma (1981-2020 reference, no SCF mask)\n",
    "============================================================\n\n",
    file = summary_combined_file, sep = ""
  )
  for (sc in timescales) {
    s <- all_summaries[[as.character(sc)]]
    cat(sprintf("\n========== SSPI-%d ==========\n", s$sc),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Calculation date: %s\n", s$date),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Grid dimensions: %d x %d\n", ncol(swe), nrow(swe)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Basin pixels (total): %d\n", basin_pixels),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Time period: %s to %s (%d months)\n",
                min(dates), max(dates), length(dates)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Reference period: %s to %s\n", sspi_ref_start, sspi_ref_end),
        file = summary_combined_file, append = TRUE)
    cat("\nZero-Inflation Screening:\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Excluded pixels (>50%% zero SWE in ref period): %d (%.1f%%)\n",
                length(invalid_pix_sspi),
                100 * length(invalid_pix_sspi) / n_pixels),
        file = summary_combined_file, append = TRUE)
    cat("\nMethodology:\n",
        file = summary_combined_file, append = TRUE)
    cat("  • No SCF domain mask applied (SSPI-1 is valid for all months)\n",
        file = summary_combined_file, append = TRUE)
    cat("  • Calendar-month stratified fitting\n",
        file = summary_combined_file, append = TRUE)
    cat("  • Zero-inflated gamma (center of probability mass)\n",
        file = summary_combined_file, append = TRUE)
    cat("  • L-moments (lmomco) for gamma parameter estimation\n",
        file = summary_combined_file, append = TRUE)
    cat("  • Empirical fallback when gamma fit fails or variance too low\n",
        file = summary_combined_file, append = TRUE)
    cat("  • SSPI values clipped at ±5\n",
        file = summary_combined_file, append = TRUE)
    cat("\nMethod Distribution:\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Gamma fitting:   %.1f%%\n", s$gamma_pct),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Empirical:       %.1f%%\n", s$empirical_pct),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Zero-variance:   %.1f%%\n", s$zero_var_pct),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Excluded:        %.1f%%\n", s$excluded_pct),
        file = summary_combined_file, append = TRUE)
    cat("\nNA Analysis:\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Total NA values in basin: %.3f%%\n", s$na_rate),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Compatibility: %s\n",
                ifelse(s$na_rate < 0.5, "✓ READY (<0.5% NAs)", "⚠ CAUTION (>0.5% NAs)")),
        file = summary_combined_file, append = TRUE)
    fi <- 1:min(12, ncol(s$sspi_indices))
    fs <- s$sspi_indices[, fi, drop = FALSE]
    cat("\nSnow drought frequency (first 12 months):\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional deficit (SSPI < -2.0): %.1f%%\n",
                100 * mean(fs < -2.0, na.rm = TRUE)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Extreme deficit     (SSPI < -1.5): %.1f%%\n",
                100 * mean(fs < -1.5, na.rm = TRUE)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Moderate deficit    (SSPI < -1.0): %.1f%%\n",
                100 * mean(fs < -1.0, na.rm = TRUE)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional surplus (SSPI >  2.0): %.1f%%\n",
                100 * mean(fs >  2.0, na.rm = TRUE)),
        file = summary_combined_file, append = TRUE)
  }
} else {
  # SWEI summary (methods 1 & 2)
  cat(
    "============================================================\n",
    "COMBINED SWEI SUMMARY FOR ALL TIMESCALES\n",
    "Seasonal Approach (3-month Smoothed SWE + k=3 SCF Mask)\n",
    "============================================================\n\n",
    file = summary_combined_file, sep = ""
  )
  for (sc in timescales) {
    s  <- all_summaries[[as.character(sc)]]
    cat(sprintf("\n\n========== SWEI-%d ==========\n", s$sc),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Calculation date: %s\n", s$date),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Grid dimensions: %d x %d\n", ncol(swe), nrow(swe)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Basin pixels: %d\n", basin_pixels),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("Time period: %s to %s (%d months)\n",
                min(dates), max(dates), length(dates)),
        file = summary_combined_file, append = TRUE)
    cat("\nMethodology:\n",
        file = summary_combined_file, append = TRUE)
    cat("  • 3-month centered moving average applied to RAW SWE\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  • SCF domain method: %s\n", scf_method_label),
        file = summary_combined_file, append = TRUE)
    if (scf_method == "1") {
      cat(sprintf("    → Data-driven threshold from 6preq_diag_SCF.R (k=3, threshold = %.3f)\n",
                  scf_threshold),
          file = summary_combined_file, append = TRUE)
    } else {
      cat(sprintf("    → Fixed %.0f%% SCF threshold, right-aligned 3-month climatology\n",
                  100 * scf_threshold),
          file = summary_combined_file, append = TRUE)
    }
    cat("  • Gringorten plotting position (non-parametric)\n",
        file = summary_combined_file, append = TRUE)
    cat("  • Zero-SWE perturbation: random uniform in (0, min_nonzero) per Huning & AghaKouchak (2020)\n",
        file = summary_combined_file, append = TRUE)
    cat("  • Minimum data requirement: >= 75% of years for a calendar month must have nonzero SWE\n",
        file = summary_combined_file, append = TRUE)
    cat("  • No clipping applied to SWEI values (full distribution retained)\n",
        file = summary_combined_file, append = TRUE)
    cat("\nMethod Distribution:\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Gringorten fitting: %.1f%%\n", s$gringorten_pct),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Masked/NA: %.1f%%\n", s$masked_pct),
        file = summary_combined_file, append = TRUE)
    cat("\nNA Analysis:\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Total NA values in basin: %.3f%%\n", s$na_rate),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Compatibility: %s\n",
                ifelse(s$na_rate < 0.5, "✓ READY (<0.5% NAs)", "⚠ CAUTION (>0.5% NAs)")),
        file = summary_combined_file, append = TRUE)
    first_idx  <- 1:min(12, ncol(s$swei_indices))
    first_swei  <- s$swei_indices[, first_idx, drop = FALSE]
    cat("\nDrought frequency (first 12 months):\n",
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Exceptional (SWEI < -2.0): %.1f%%\n",
                100 * mean(first_swei < -2.0, na.rm = TRUE)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Extreme     (SWEI < -1.6): %.1f%%\n",
                100 * mean(first_swei < -1.6, na.rm = TRUE)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Severe      (SWEI < -1.3): %.1f%%\n",
                100 * mean(first_swei < -1.3, na.rm = TRUE)),
        file = summary_combined_file, append = TRUE)
    cat(sprintf("  Moderate    (SWEI < -0.8): %.1f%%\n",
                100 * mean(first_swei < -0.8, na.rm = TRUE)),
        file = summary_combined_file, append = TRUE)
  }
}

cat("\n============================================================\n")
cat(if (scf_method == "3") "SSPI-1 CALCULATION COMPLETE!\n" else "SWEI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat(sprintf("Summary file: %s\n", summary_combined_file))
cat("\n✓✓✓ READY FOR ANALYSIS ✓✓✓\n")