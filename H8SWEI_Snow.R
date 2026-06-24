##############################################
# SWEI / SSPI-1 CALCULATION - PARALLEL VERSION (SPI-STYLE OUTPUTS)
# Uses parallel::parLapply - Works on Windows and Linux
# Output naming matches SPI_ERALand.R convention
#
# METHOD (user-selected at runtime):
#   1 = Huning & AghaKouchak (2020) fixed 5% SCF threshold (SWEI)
#       Right-aligned 3-month climatology, matches original paper.
#   2 = SSPI-1: Standardized SnowPack Index (monthly, no SCF mask)
#       Zero-inflated gamma distribution, following Stagge et al. (2015) and JRC EDO (2020).
#       Output to: sspi_results_monthly/
#
# Other methodological alignments with Huning & AghaKouchak (2020):
#   - Zero-SWE perturbation: random uniform in (0, min_nonzero)  [method 1]
#   - Minimum data: >= 75% of years per calendar month must be nonzero [method 1]
#   - No clipping of SWEI values  [method 1]
#
# RUNNING THE SCRIPT:
#   Interactive (RStudio / R console): source("8SWEI_Snow_fixed.R")
#   Command line with argument:        Rscript 8SWEI_Snow_fixed.R 1
#   Batch / non-interactive fallback:  edit the DEFAULT_SCF_METHOD line below
##############################################

# ---- Libraries ----
library(terra)
library(ncdf4)
library(zoo)
library(writexl)
library(parallel)
library(lmomco)   # required for option 2 (SSPI-1 zero-inflated gamma)

# ---- USER SELECTION: SCF DOMAIN METHOD ----
cat("\n============================================================\n")
cat("SCF DOMAIN DEFINITION: SELECT METHOD\n")
cat("============================================================\n")
cat("  1 = Huning & AghaKouchak (2020) ??? fixed 5% SCF threshold\n")
cat("                          (standard approach from the original global paper)\n")
cat("  2 = SSPI-1 (Standardized SnowPack Index, monthly)\n")
cat("                          (zero-inflated gamma, 1981-2020 reference, no SCF mask)\n\n")

# ---- To use in batch/scheduled jobs, set DEFAULT_SCF_METHOD to "1" or "2"
DEFAULT_SCF_METHOD <- "1"

scf_method <- ""

# Priority 1: command-line argument (e.g. Rscript script.R 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && trimws(args[1]) %in% c("1", "2")) {
  scf_method <- trimws(args[1])
  cat(sprintf("??? Method set from command-line argument: %s\n", scf_method))
  
  # Priority 2: RStudio GUI dialog (works with Run button AND source())
} else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable()) {
  raw_input <- rstudioapi::showPrompt(
    title   = "SCF Domain Method",
    message = paste0(
      "Enter method number (1 or 2):\n\n",
      "  1 = Huning & AghaKouchak (2020)  [default]\n",
      "  2 = SSPI-1 (zero-inflated gamma)\n"
    ),
    default = DEFAULT_SCF_METHOD
  )
  if (is.null(raw_input)) {
    # User clicked Cancel ??? use the default
    scf_method <- DEFAULT_SCF_METHOD
    cat(sprintf("??? Dialog cancelled: using default Method %s\n", scf_method))
  } else {
    scf_method <- trimws(raw_input)
  }
  while (!scf_method %in% c("1", "2")) {
    raw_input <- rstudioapi::showPrompt(
      title   = "Invalid input ??? try again",
      message = "Please enter 1 or 2:",
      default = DEFAULT_SCF_METHOD
    )
    if (is.null(raw_input)) {
      scf_method <- DEFAULT_SCF_METHOD
      break
    }
    scf_method <- trimws(raw_input)
  }
  cat(sprintf("??? Method set from RStudio dialog: %s\n", scf_method))
  
  # Priority 3: plain readline() fallback (R console, non-RStudio interactive)
} else if (interactive()) {
  while (!scf_method %in% c("1", "2")) {
    scf_method <- trimws(readline(prompt = "Enter 1 or 2: "))
    if (!scf_method %in% c("1", "2")) cat("  Invalid input. Please enter 1 or 2.\n")
  }
  
  # Priority 4: non-interactive fallback default (batch / Rscript without args)
} else {
  scf_method <- DEFAULT_SCF_METHOD
  cat(sprintf("??? Non-interactive session detected: using default Method %s\n", scf_method))
  cat("  (Edit DEFAULT_SCF_METHOD at the top of the script to change this.)\n")
}

if (scf_method == "1") {
  cat("\n??? Method selected: Huning & AghaKouchak (2020) fixed 5% SCF threshold\n")
} else {
  cat("\n??? Method selected: SSPI-1 (zero-inflated gamma, monthly, 1981-2020 reference)\n")
}

scf_method_label <- if (scf_method == "1") "Huning & AghaKouchak (2020)" else "SSPI-1 (zero-inflated gamma, no SCF mask)"

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")
out_dir <- if (scf_method == "2") "sspi_results_monthly" else "swei_results_seasonal"
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
  cat("??? Basin boundary loaded\n")
} else {
  stop("Basin boundary not found: ", basin_path)
}

# ---- Input files ----
swe_file <- "monthly_data_direct/snow_depth_water_equivalent_monthly.nc"
scf_file <- "monthly_data_direct/snow_cover_monthly.nc"
if (!file.exists(swe_file)) stop("SWE file not found: ", swe_file)
if (scf_method == "1" && !file.exists(scf_file))
  stop("SCF file not found (required for method 1): ", scf_file)

cat("\n============================================================\n")
cat(if (scf_method == "2") "SSPI-1 CALCULATION (PARALLEL MODE)\n"
    else "SEASONAL SWEI CALCULATION (PARALLEL MODE)\n")
cat("============================================================\n\n")

cat("===== READING NETCDF FILES =====\n")
swe <- rast(swe_file)
if (scf_method == "1") {
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

# ---- Align Spatial Extents and Reproject ----
cat("\n===== SPATIAL ALIGNMENT =====\n")

if (!is.null(scf_stack)) {
  # 1. Fix Longitude Mismatches (0-360 vs -180-180)
  swe_ext <- ext(swe)
  scf_ext <- ext(scf_stack)
  
  if (scf_ext[2] > 180 && swe_ext[2] <= 180) {
    cat("??? Rotating SCF from 0-360 to -180-180 to match SWE...\n")
    scf_stack <- terra::rotate(scf_stack)
  } else if (swe_ext[2] > 180 && scf_ext[2] <= 180) {
    cat("??? Rotating SWE from 0-360 to -180-180 to match SCF...\n")
    swe <- terra::rotate(swe)
  }
}

# 2. Crop early to save memory and isolate the basin
if (!is.null(basin)) {
  cat("??? Cropping rasters to basin extent...\n")
  swe <- crop(swe, project(basin, crs(swe)))
  if (!is.null(scf_stack)) {
    scf_stack <- crop(scf_stack, project(basin, crs(scf_stack)))
  }
}

# 3. Resample to common grid
if (!is.null(scf_stack) && !compareGeom(scf_stack[[1]], swe[[1]], stopOnError = FALSE)) {
  cat("??? Resampling SWE to match SCF grid...\n")
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
cat(sprintf("??? Basin masking complete: %d pixels (%.1f%% of raster)\n",
            basin_pixels, 100 * basin_pixels / total_pixels))

# Safety Check against NaN calculation
if (is.na(basin_pixels) || basin_pixels == 0) {
  stop("FATAL ERROR: Basin masking resulted in 0 valid pixels. The SWE and SCF rasters may not overlap the basin boundary properly after resampling.")
}

# ---- Convert SWE to mm if needed ----
swe_mean_sample <- global(swe[[1]], "mean", na.rm = TRUE)$mean
if (!is.na(swe_mean_sample) && swe_mean_sample < 10) {
  swe <- swe * 1000
  cat("??? Converted SWE from meters to mm\n")
}

# ---- Convert SCF to fraction if needed ----
if (!is.null(scf_stack)) {
  scf_max <- global(scf_stack[[1]], "max", na.rm = TRUE)$max
  if (!is.na(scf_max) && scf_max > 1.5) {
    scf_stack <- scf_stack / 100
    cat("??? Converted SCF from percent to fraction\n")
  }
}

# ---- Prepare matrices ----
swe_matrix <- values(swe, mat = TRUE)
dates <- dates_swe
months_all <- as.integer(format(dates, "%m"))
years_all <- as.integer(format(dates, "%Y"))
n_pixels <- nrow(swe_matrix)
cat(sprintf("Processing %d pixels of the rectangle covering basin ...\n", n_pixels))

# ==============================================================================
#   USER CHOICE: SCF DOMAIN DEFINITION METHOD
# ==============================================================================
cat(sprintf("\n??? Proceeding with Method %s: %s\n", scf_method, scf_method_label))

cat("\n===== STEP 2: Building SCF mask =====\n")
target_k <- 3

if (scf_method == "2") {
  
  # ============================================================================
  #   BRANCH 2 ??? SSPI-1 (Zero-inflated Gamma, Monthly, No SCF Mask)
  # ============================================================================
  cat("\n??? SSPI-1 selected: skipping SCF mask build (not applicable).\n")
  cat("  Zero-inflated gamma fitting will be applied per calendar month.\n")
  scf_mask_list  <- NULL
  basin_scf_mask <- NULL
  
  # Reference period for SSPI-1
  sspi_ref_start <- min(dates) 
  sspi_ref_end   <- max(dates) 
  ref_idx        <- which(dates >= sspi_ref_start & dates <= sspi_ref_end)
  if (length(ref_idx) < 400) {
    warning(sprintf("SSPI-1 reference period has only %d months (recommended: >=480)", length(ref_idx)))
  }
  cat(sprintf("??? SSPI-1 reference period: %s to %s (%d months)\n",
              dates[ref_idx[1]], dates[ref_idx[length(ref_idx)]], length(ref_idx)))
  
  # Zero-inflation screening
  cat("\n===== ZERO-INFLATION SCREENING (SSPI-1) =====\n")
  ref_swe_sspi   <- swe_matrix[, ref_idx, drop = FALSE]
  zero_prop_sspi <- rowMeans(ref_swe_sspi <= 0, na.rm = TRUE)
  invalid_pix_sspi <- which(zero_prop_sspi > 0.50)
  valid_pix_sspi   <- which(zero_prop_sspi <= 0.50)
  cat(sprintf("??? Valid pixels (<=50%% zero SWE in ref period): %d (%.1f%%)\n",
              length(valid_pix_sspi), 100 * length(valid_pix_sspi) / n_pixels))
  if (length(invalid_pix_sspi) > 0) {
    cat(sprintf("  Excluded pixels (>50%% zero SWE): %d\n", length(invalid_pix_sspi)))
    swe_matrix[invalid_pix_sspi, ] <- NA
  }
  
} else if (scf_method == "1") {
  
  # ============================================================================
  #   BRANCH 1 ??? HUNING & AGHAKOUCHAK (2020) fixed 5% SCF threshold
  # ============================================================================
  scf_threshold <- 0.05
  cat(sprintf("??? H&A (2020): building SCF mask with fixed %.0f%% threshold\n",
              100 * scf_threshold))
  
  # Mask SCF stack to the basin so pixel indices align exactly with SWE
  scf_basin      <- mask(scf_stack, basin, inverse = FALSE, touches = TRUE)
  scf_matrix_ha  <- values(scf_basin, mat = TRUE)
  months_scf_ha  <- as.integer(format(dates_scf, "%m"))
  
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
  rm(scf_basin, scf_matrix_ha)
  
  cat("\n===== COMPUTING BASIN-LEVEL SCF MASK (H&A 2020) =====\n")
  basin_pixel_idx <- which(!is.na(swe_matrix[, 1]))
  basin_scf_mask  <- logical(12)
  for (m in 1:12) {
    mv_basin          <- scf_mask_list[[m]][basin_pixel_idx]
    basin_scf_mask[m] <- length(mv_basin) > 0 && mean(mv_basin, na.rm = TRUE) >= 0.5
  }
}

if (scf_method == "1") {
  cat(sprintf("??? Basin SCF mask: months included = %s\n",
              paste(month.abb[which(basin_scf_mask)], collapse = ", ")))
}

cat("\n===== COMPUTING BASIN-AVERAGED SWE =====\n")
if (scf_method == "1") {
  swe_basin_avg_raw <- colMeans(swe_matrix, na.rm = TRUE)
  cat(sprintf("??? Basin-averaged raw SWE: %d time steps\n", length(swe_basin_avg_raw)))
} else {
  swe_basin_avg_raw <- NULL
  cat("??? Basin-averaged SWE: skipped (not used by SSPI-1)\n")
}

swe_basin_avg_smoothed <- NULL
basin_avg_swei_results <- list()

# ==============================================================================
#   Helper Functions
# ==============================================================================
clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

gringorten_swei_seasonal <- function(v_smoothed, dates_vec, scf_mask_list,
                                     cell_index = NULL,
                                     basin_scf_mask = NULL,
                                     eps = 1e-6) {
  if (length(v_smoothed) == 0 || all(is.na(v_smoothed)) || is.null(v_smoothed)) {
    return(list(swei = rep(NA_real_, length(dates_vec)), method = rep(NA_integer_, length(dates_vec))))
  }
  v_clean <- v_smoothed
  v_clean[!is.finite(v_clean)] <- NA_real_
  n <- length(v_clean)
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n) 
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    all_m_idx <- which(mon == m)
    if (length(all_m_idx) == 0) next
    nonzero_frac <- sum(v_clean[all_m_idx] > 0, na.rm = TRUE) / length(all_m_idx)
    if (nonzero_frac < 0.75) next
    
    idx <- which(mon == m & is.finite(v_clean))
    if (length(idx) == 0) next
    
    if (is.null(scf_mask_list) || length(scf_mask_list) < m) next
    if (!is.null(cell_index)) {
      mask_vec <- scf_mask_list[[m]]
      if (length(mask_vec) < cell_index || !mask_vec[cell_index]) next
    } else if (!is.null(basin_scf_mask)) {
      if (length(basin_scf_mask) < m || !basin_scf_mask[m]) next
    }
    
    samp <- v_clean[idx]
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

monthly_sspi <- function(v, dates_vec, ref_idx, eps = 1e-6) {
  n <- length(v)
  if (n == 0 || all(is.na(v)) || is.null(v)) {
    return(list(sspi = rep(NA_real_, n), method = rep(NA_integer_, n)))
  }
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n) 
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
#   MAIN CALCULATION LOOP - PARALLEL 
# ==============================================================================
timescales <- if (scf_method == "2") c(1) else c(3)

# ---- Set up parallel cluster ----
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\n??? Starting parallel cluster with %d cores\n", n_cores))
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, iseed = 40)  
if (scf_method == "2") {
  clusterExport(cl, varlist = c("monthly_sspi", "clip_prob",
                                "swe_matrix", "dates", "ref_idx"),
                envir = environment())
  clusterEvalQ(cl, { library(zoo); library(lmomco) })
} else {
  clusterExport(cl, varlist = c("gringorten_swei_seasonal", "clip_prob",
                                "dates", "scf_mask_list",
                                "basin_scf_mask"),
                envir = environment())
  clusterEvalQ(cl, { library(zoo) })
}

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# ---- Summary file setup ----
summary_combined_file <- file.path(out_dir,
                                   if (scf_method == "2") "sspi_summary.txt" else "swei_all_timescales_summary.txt")
all_summaries <- list()

for (sc in timescales) {
  
  if (scf_method == "2") {
    # ============================================================================
    #   OPTION 2: SSPI-1 CALCULATION
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
    
    cat("  -> Saving CSV files...\n")
    for (m in 1:12) {
      idx_m <- which(months_all_sc == m); if (length(idx_m) == 0) next
      sub  <- sspi_indices[, idx_m, drop = FALSE]
      df   <- as.data.frame(sub)
      colnames(df) <- format(dates[idx_m], "%Y")
      coords <- xyFromCell(swe[[1]], 1:nrow(sub))
      df <- cbind(lon = coords[,1], lat = coords[,2], df)
      write.csv(df, file.path(out_dir, sprintf("sspi_%02d_month%02d_%s.csv", sc, m, month_names[m])), row.names = FALSE)
    }
    
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
    cat(sprintf("??? SSPI-%d complete\n", sc))
    gc()
    next  
  }
  
  # ============================================================================
  #   OPTION 1: SWEI CALCULATION
  # ============================================================================
  cat(sprintf("\n===== SWEI-%d (%d-month window, Sequential Calculation) =====\n", sc, sc))
  
  swe_basin_avg_smoothed <- zoo::rollapply(swe_basin_avg_raw,
                                           width = sc,
                                           FUN   = sum,
                                           align = "right",
                                           fill  = NA,
                                           na.rm = FALSE)
  cat(sprintf("??? Basin-averaged %d-month smoothed SWE: %d steps, %.1f%% non-NA\n",
              sc, length(swe_basin_avg_smoothed),
              100 * mean(!is.na(swe_basin_avg_smoothed))))
  
  swe_smoothed <- matrix(NA, nrow = nrow(swe_matrix), ncol = ncol(swe_matrix))
  for (i in 1:nrow(swe_matrix)) {
    swe_smoothed[i, ] <- zoo::rollapply(swe_matrix[i, ],
                                        width = sc,
                                        FUN = sum,
                                        align = "right",
                                        fill = NA,
                                        na.rm = FALSE) 
  }
  cat(sprintf("??? %d-month rolling sum applied to SWE (SWEI preprocessing)\n", sc))
  
  clusterExport(cl, varlist = c("swe_smoothed"), envir = environment())
  
  avg_swei_out <- tryCatch(
    gringorten_swei_seasonal(swe_basin_avg_smoothed, dates, scf_mask_list,
                             cell_index = NULL, basin_scf_mask = basin_scf_mask,
                             eps = 1e-6),
    error = function(e) list(swei = rep(NA_real_, length(dates)), method = rep(NA_integer_, length(dates)))
  )
  basin_avg_swei_results[[as.character(sc)]] <- avg_swei_out$swei
  
  cat(sprintf("Processing %d pixels in parallel... ", n_pixels))
  start_time <- Sys.time()
  
  pixel_list <- parLapply(cl, seq_len(n_pixels), function(i) {
    tryCatch(
      gringorten_swei_seasonal(swe_smoothed[i, ], dates, scf_mask_list, cell_index = i, eps = 1e-6),
      error = function(e) list(swei = rep(NA_real_, length(dates)),
                               method = rep(NA_integer_, length(dates)))
    )
  })
  
  swei_indices <- do.call(rbind, lapply(pixel_list, `[[`, "swei"))
  method_matrix <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
  
  extreme_dry  <- sum(swei_indices < -3.0, na.rm = TRUE)
  extreme_wet  <- sum(swei_indices >  3.0, na.rm = TRUE)
  cat(sprintf("  Extreme values (|SWEI| > 3): %d dry, %d wet (retained, no clipping)\n",
              extreme_dry, extreme_wet))
  
  basin_pix_rows      <- which(!is.na(swe_matrix[, 1]))
  
  # Protect against empty vectors generating NaN
  if (length(basin_pix_rows) == 0) {
    na_rate_grid <- NA
  } else {
    na_rate_grid <- 100 * mean(is.na(swei_indices[basin_pix_rows, , drop = FALSE]))
  }
  
  basin_mean_complete <- 100 * mean(!is.na(basin_avg_swei_results[[as.character(sc)]]))
  
  cat(sprintf("  Basin mean series:  %.1f%% complete (%.1f%% NA)\n",
              basin_mean_complete, 100 - basin_mean_complete))
  cat(sprintf("  Per-pixel gridded:  %.3f%% NA  (basin pixels only)\n", na_rate_grid))
  
  grid_status <- if (is.na(na_rate_grid)) {
    "UNKNOWN (0 valid basin pixels)" 
  } else if (na_rate_grid < 0.5) {
    "READY  (<0.5% NAs)"
  } else {
    "CAUTION (>0.5% NAs)"
  }
  
  cat(sprintf("  Gridded compatibility: %s\n", grid_status))
  
  if (!is.na(na_rate_grid) && na_rate_grid >= 0.5) {
    cat("  ???  High per-pixel NA is expected: the 75%%-nonzero-years criterion\n")
    cat("     and SCF mask are applied pixel-by-pixel, which is much stricter\n")
    cat("     than the basin-mean aggregate implies. The basin mean series is\n")
    cat("     suitable for time-series analysis; use gridded files for spatial\n")
    cat("     analysis with awareness of actual spatial coverage.\n")
  }
  
  gringorten_pct  <- 100 * sum(method_matrix == 1, na.rm = TRUE) / sum(!is.na(method_matrix))
  masked_pct      <- 100 * sum(method_matrix == 2, na.rm = TRUE) / sum(!is.na(method_matrix))
  
  cat(sprintf("  Methods: Gringorten=%.1f%%, Masked/NA=%.1f%%\n",
              gringorten_pct, masked_pct))
  
  all_summaries[[as.character(sc)]] <- list(
    sc = sc, date = Sys.time(), na_rate = na_rate_grid, 
    basin_mean_complete = basin_mean_complete,
    gringorten_pct = gringorten_pct, masked_pct = masked_pct,
    swei_indices = swei_indices, method_matrix = method_matrix
  )
  
  months_all <- as.numeric(format(dates, "%m"))
  cat("  -> Saving CSV files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
    swei_subset <- swei_indices[, idx_m, drop = FALSE]
    dates_subset <- dates[idx_m]
    
    df <- as.data.frame(swei_subset)
    colnames(df) <- format(dates_subset, "%Y")
    coords <- xyFromCell(swe[[1]], 1:nrow(swei_subset))
    df <- cbind(lon = coords[,1], lat = coords[,2], df)
    
    csv_file <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.csv", sc, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE)
  }
  
  cat("  -> Saving Excel workbook...\n")
  excel_data <- list()
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
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
  
  cat("  -> Creating distribution map...\n")
  swei_dist_mapping <- data.frame(
    code = c(1, 2), name = c("Gringorten", "Masked/NA"),
    color = c("#4575b4", "#91bfdb"), stringsAsFactors = FALSE
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
  plot(swei_dist_raster, col = swei_dist_mapping$color, breaks = c(0.5, 1.5, 2.5),
       legend = FALSE, main = sprintf("SWEI-%d: Seasonal Method", sc), axes = FALSE, box = FALSE)
  if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  valid_pixels <- sum(!is.na(values(swe[[1]])))
  mtext(sprintf("Valid pixels: %d | Gringorten: %.1f%%", valid_pixels, gringorten_pct),
        side = 1, line = 1, cex = 0.75, col = "gray30")
  par(mar = c(0, 0, 0, 0))
  plot.new()
  dist_counts <- table(factor(dist_sim, levels = 1:2, labels = swei_dist_mapping$name))
  legend_entries <- sprintf("%s (%d pixels)", names(dist_counts), dist_counts)
  legend("center", legend = legend_entries, fill = swei_dist_mapping$color,
         title = "Method", cex = 0.95, bg = "white", bty = "o", border = "gray50")
  dev.off()
  layout(1)
  
  cat("  -> Saving NetCDF files...\n")
  for (m in 1:12) {
    idx_m <- which(months_all == m); if (length(idx_m) == 0) next
    swei_rast  <- rep(rast(swe[[1]]), length(idx_m))
    sub_matrix  <- swei_indices[, idx_m, drop = FALSE]
    if (nrow(sub_matrix) != ncell(swei_rast)) {
      stop(sprintf("Dimension mismatch: Matrix %d rows vs Raster %d cells.",
                   nrow(sub_matrix), ncell(swei_rast)))
    }
    values(swei_rast)  <- sub_matrix
    terra::time(swei_rast)  <- dates[idx_m]
    nc_file  <- file.path(out_dir, sprintf("swei_%02d_month%02d_%s.nc", sc, m, month_names[m]))
    writeCDF(swei_rast, nc_file, varname = "swei",
             longname = sprintf("Seasonal SWEI (SWEI-%d)", sc),
             unit = "standardized_index", missval = -9999, overwrite = TRUE)
  }
  cat(sprintf("??? SWEI-%d complete\n", sc))
  gc()
}

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\n??? Parallel cluster stopped\n")

# ==============================================================================
#   SAVE BASIN-AVERAGED SWEI TO CSV 
# ==============================================================================
cat("\n===== SAVING BASIN-AVERAGED SWEI =====\n")
months_all  <- as.integer(format(dates, "%m"))
years_all   <- as.integer(format(dates, "%Y"))
month_names_save  <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
if (scf_method != "2") {
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
    csv_file <- file.path(out_dir, sprintf("swei_%02d_basin_averaged_by_month.csv", sc))
    write.csv(df_out, csv_file, row.names = FALSE, na = "")
    cat(sprintf("??? Saved basin-averaged SWEI-%d (12 monthly series) to: %s\n", sc, csv_file))
  }
}

# ==============================================================================
#   WRITE SUMMARY FILE AT END 
# ==============================================================================
cat("\n===== WRITING SUMMARY FILE =====\n")

if (scf_method == "2") {
  cat(
    "============================================================\n",
    "MONTHLY SSPI-1 SUMMARY - NECHAKO BASIN\n",
    "ERA5-Land Snow Water Equivalent\n",
    "Method 2: Zero-Inflated Gamma (1981-2020 reference, no SCF mask)\n",
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
    cat("  - No SCF domain mask applied (SSPI-1 is valid for all months)\n  - Calendar-month stratified fitting\n  - Zero-inflated gamma (center of probability mass)\n  - L-moments (lmomco) for gamma parameter estimation\n  - Empirical fallback when gamma fit fails or variance too low\n  - SSPI values clipped at +/-5\n", file = summary_combined_file, append = TRUE)
    cat("\nMethod Distribution:\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  Gamma fitting:   %.1f%%\n", s$gamma_pct), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Empirical:       %.1f%%\n", s$empirical_pct), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Zero-variance:   %.1f%%\n", s$zero_var_pct), file = summary_combined_file, append = TRUE)
    cat(sprintf("  Excluded:        %.1f%%\n", s$excluded_pct), file = summary_combined_file, append = TRUE)
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
    "SCF mask window fixed at k=3 per Huning & AghaKouchak (2020)\n",
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
    cat("  - SCF mask window fixed at k=3 (intentionally independent of SWE window,\n    per Huning & AghaKouchak 2020 ??? see target_k comment in script)\n", file = summary_combined_file, append = TRUE)
    cat(sprintf("  - SCF domain method: %s\n", scf_method_label), file = summary_combined_file, append = TRUE)
    cat(sprintf("    -> Fixed %.0f%% SCF threshold, right-aligned 3-month climatology\n", 100 * 0.05), file = summary_combined_file, append = TRUE)
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
      cat("  NOTE: The large gap between basin-mean completeness and per-pixel NA\n  is expected ??? spatial averaging masks individual pixel gaps, while the\n  75%%-nonzero-years criterion and SCF mask applied pixel-by-pixel are\n  much stricter. The basin mean series is suitable for time-series analysis;\n  use gridded files for spatial analysis with awareness of spatial coverage.\n", file = summary_combined_file, append = TRUE)
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

cat("\n============================================================\n")
cat(if (scf_method == "2") "SSPI-1 CALCULATION COMPLETE!\n" else "SWEI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat(sprintf("Summary file: %s\n", summary_combined_file))
cat("\n--- READY FOR ANALYSIS ---\n")