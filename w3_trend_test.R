# ####################################################################################
# w3_trend_test.R  ·  SPATIAL TREND ANALYSIS FOR ALL INDICES
# ─────────────────────────────────────────────────────────────────────────────────
# Computes per-pixel statistics for SPI, SPEI (multiple scales) and SWEI:
#   • Mann-Kendall (variance-corrected, VC)  +  Sen's slope
# • Mann-Kendall (TFPW – Trend-Free Pre-Whitening, Yue et al. 2002)
# • Drought event count / mean duration / max intensity
# • CUSUM regime-shift year
# • Wald-Wolfowitz runs test (temporal clustering)
# • n_spectral_peaks  (placeholder = 0L; reserved for future spectral analysis)
# OUTPUT : {TREND_DIR}/{index}_{scale:02d}_results.csv  (one per index × scale)
# Columns include tau_vc, p_value_vc, tau_tfpw, p_value_tfpw,
# n_events, mean_duration, max_intensity, regime_shift_year,
# p_value_runs, clustering, n_spectral_peaks, and time_001…time_NNN.
# Run BEFORE w4_trends_visualization.R
# ####################################################################################
source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "data.table", "Kendall", "trend", "parallel", "lubridate"))
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
dir.create(TREND_DIR, showWarnings = FALSE, recursive = TRUE)
# ── Parallel back-end (Windows-compatible) ────────────────────────────
is_windows <- .Platform$OS.type == "windows"
N_CORES    <- max(1L, parallel::detectCores() - 1L)
cl         <- NULL
if (N_CORES > 1) {
if (is_windows) {
cl <- parallel::makeCluster(N_CORES)
parallel::clusterEvalQ(cl, { library(Kendall); library(trend) })
cat(sprintf("✓ Windows cluster: %d cores\n", N_CORES))
} else {
cat(sprintf("✓ Unix fork: %d cores\n", N_CORES))
}
} else {
cat("ℹ Single-core mode\n")
}
# ── Unified 3-arg wrapper for the SWEI file finder ───────────────────
# process_index() always calls find_fn(dir, index_type, scale).
# find_swei_seasonal_files() only takes 2 args, so we wrap it here.
find_swei_nc_wrapper <- function(data_dir, index_type, scale) {
find_swei_seasonal_files(data_dir, scale)
}
####################################################################################
# VECTORISED STATISTICS FUNCTIONS
####################################################################################
# ── 1. Variance-corrected Mann-Kendall + Sen's slope ─────────────────
#'  Returns tau_vc, p_value_vc, sl_vc (Sen's slope), filtered_vc.
vectorized_mann_kendall <- function(ts_matrix) {
  n_pix <- nrow(ts_matrix)
  process_pixel <- function(i) {
    x <- ts_matrix[i, ]; x <- x[!is.na(x)]
    if (length(x) < 10 || var(x, na.rm = TRUE) < 1e-6)
      return(list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
    tryCatch({
      mk  <- Kendall::MannKendall(x)
      n   <- length(x)
      idx <- which(!upper.tri(matrix(0, n, n)), arr.ind = TRUE)
      idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
      slopes <- (x[idx[, 2]] - x[idx[, 1]]) / (idx[, 2] - idx[, 1])
      list(tau      = as.numeric(mk$tau),
           pval     = as.numeric(mk$sl),
           slope    = median(slopes, na.rm = TRUE),
           filtered = FALSE)
    }, error = function(e) list(tau = NA, pval = NA, slope = NA, filtered = TRUE))
  }
  do_par <- N_CORES > 1 && n_pix > 100
  res <- if (do_par && is_windows)  parallel::parLapply(cl, seq_len(n_pix), process_pixel)
  else if (do_par)            parallel::mclapply(seq_len(n_pix), process_pixel, mc.cores = N_CORES)
  else                        lapply(seq_len(n_pix), process_pixel)
  data.table::data.table(
    tau_vc      = sapply(res,  `[[`,  "tau"),
    p_value_vc  = sapply(res,  `[[`,  "pval"),
    sl_vc       = sapply(res,  `[[`,  "slope"),
    filtered_vc = sapply(res,  `[[`,  "filtered"))
}
# ── 2. TFPW Mann-Kendall (Yue et al. 2002) ───────────────────────────
#'  Removes trend + pre-whitens lag-1 autocorrelation before MK test.
#'  Returns tau_tfpw, p_value_tfpw, filtered_tfpw.
vectorized_mann_kendall_tfpw <- function(ts_matrix) {
  n_pix <- nrow(ts_matrix)
  process_pixel_tfpw <- function(i) {
    x <- ts_matrix[i, ]; x <- x[!is.na(x)]
    if (length(x) < 12 || var(x, na.rm = TRUE) < 1e-6)
      return(list(tau = NA, pval = NA, filtered = TRUE))
    tryCatch({
      n <- length(x)
      # Sen's slope via trend::sens.slope — avoids the O(n²) matrix approach
      # and correctly uses all n*(n-1)/2 pairwise slopes.
      sens_r <- trend::sens.slope(x)
      beta   <- as.numeric(sens_r$estimates)
      if (!is.finite(beta)) return(list(tau = NA, pval = NA, filtered = TRUE))
      
      # Detrend
      t_seq <- seq_len(n)
      x_det <- x - beta * t_seq
      
      # Lag-1 autocorrelation of detrended series
      r1 <- cor(x_det[-n], x_det[-1], use = "complete.obs")
      if (is.na(r1)) r1 <- 0
      r1 <- max(-0.99, min(0.99, r1))   # clip to avoid instability
      
      # Pre-whiten: remove lag-1 autocorrelation (result has length n-1)
      x_pw <- x_det[-1] - r1 * x_det[-n]
      
      # Blend trend back in (Yue et al. 2002, step 4)
      x_final <- x_pw + beta * t_seq[-1]
      
      if (length(x_final) < 10) return(list(tau = NA, pval = NA, filtered = TRUE))
      
      mk <- Kendall::MannKendall(x_final)
      list(tau = as.numeric(mk$tau), pval = as.numeric(mk$sl), filtered = FALSE)
    }, error = function(e) list(tau = NA, pval = NA, filtered = TRUE))
  }
  do_par <- N_CORES > 1 && n_pix > 100
  res <- if (do_par && is_windows)  parallel::parLapply(cl, seq_len(n_pix), process_pixel_tfpw)
  else if (do_par)            parallel::mclapply(seq_len(n_pix), process_pixel_tfpw, mc.cores = N_CORES)
  else                        lapply(seq_len(n_pix), process_pixel_tfpw)
  data.table::data.table(
    tau_tfpw      = sapply(res, `[[`, "tau"),
    p_value_tfpw  = sapply(res, `[[`, "pval"),
    filtered_tfpw = sapply(res, `[[`, "filtered"))
}
# ── 3. Drought event characteristics ─────────────────────────────────
#'  Returns n_events, mean_duration, max_intensity per pixel.
vectorized_event_detection <- function(ts_matrix, threshold = DROUGHT_ONSET) {
  n_pix <- nrow(ts_matrix)
  res   <- data.table::data.table(
    n_events      = rep(NA_real_, n_pix),
    mean_duration = rep(NA_real_, n_pix),
    max_intensity = rep(NA_real_, n_pix))
  for (i in seq_len(n_pix)) {
    x  <- ts_matrix[i, ]; x  <- x[!is.na(x)]
    if (!length(x)) next
    in_ev  <- x  < threshold
    if (!any(in_ev)) {
      res[i, c("n_events", "mean_duration", "max_intensity") := list(0, 0, 0)]
      next
    }
    starts  <- which(diff(c(FALSE, in_ev)) ==  1)
    ends    <- which(diff(c(in_ev, FALSE)) == -1)
    durs    <- ends - starts + 1
    ints    <- mapply(function(s, e) min(x[s:e]), starts, ends)
    res[i,  `:=` (n_events      = length(starts),
                  mean_duration = mean(durs),
                  max_intensity = min(ints))]
  }
  res
}
# ── 4. CUSUM regime-shift year ────────────────────────────────────────
#'  Identifies the year of maximum CUSUM deviation (most likely shift point).
vectorized_regime_shift <- function(ts_matrix, year_vec) {
  n_pix <- nrow(ts_matrix)
  years <- rep(NA_real_, n_pix)
  for (i in seq_len(n_pix)) {
    x  <- ts_matrix[i, ]; ok <- !is.na(x)
    if (sum(ok) < 20) next
    xc <- x[ok]; yc <- year_vec[ok]
    cs <- cumsum(xc - mean(xc))
    years[i] <- yc[which.max(abs(cs))]
  }
  years
}
# ── 5. Wald-Wolfowitz runs test (temporal clustering) ─────────────────
#'  Returns p_value_runs, clustering ("clustered"/"dispersed"/"random"),
#'  filtered_runs (TRUE if series too short or degenerate).
vectorized_runs_test <- function(ts_matrix, threshold = 0) {
  n_pix <- nrow(ts_matrix)
  res   <- data.table::data.table(
    p_value_runs  = rep(NA_real_,      n_pix),
    clustering    = rep(NA_character_, n_pix),
    filtered_runs = rep(TRUE,          n_pix))
  for (i in seq_len(n_pix)) {
    x  <- ts_matrix[i, ]; x  <- x[!is.na(x)]
    if (length(x)  < 10) next
    b   <- as.integer(x  >= threshold)
    rn  <- rle(b); nr  <- length(rn$lengths); n  <- length(b)
    n1  <- sum(b); n0  <- n - n1
    if (n1 == 0 || n0 == 0) next
    er   <- 2 * n0 * n1 / n + 1
    vr   <- 2 * n0 * n1 * (2 * n0 * n1 - n) / (n^2 * (n - 1))
    if (vr  <= 0) next
    z    <- (nr - er) / sqrt(vr)
    pv   <- 2 * pnorm(-abs(z))
    clus  <- if (nr  < er)  "clustered" else if (nr  > er)  "dispersed" else  "random"
    res[i,  `:=` (p_value_runs  = pv,
                  clustering    = clus,
                  filtered_runs = FALSE)]
  }
  res
}
####################################################################################
# MAIN PROCESSING FUNCTION
####################################################################################
process_index <- function(index_type, scales, find_fn, seas_dir) {
  for (sc in scales) {
    cat(sprintf("\n╔══ %s-%02d ══╗\n", toupper(index_type), sc))
    # [1/5] Find monthly NC files ──────────────────────────────────
    cat("[1/5] Finding files...\n")
    files  <- find_fn(seas_dir, index_type, sc)   # unified 3-arg call
    if (!length(files)) { cat("  ❌ No files found\n"); next }
    cat(sprintf("  Found %d monthly NC files\n", length(files)))
    
    # [2/5] Load raster stacks ────────────────────────────────────
    cat("[2/5] Loading rasters...\n")
    stacks  <- lapply(files, function(f) {
      tryCatch(terra::rast(f),
               error = function(e) { cat("  ⚠ Failed: ", basename(f), "\n"); NULL })
    })
    stacks  <- stacks[!sapply(stacks, is.null)]
    if (!length(stacks)) { cat("  ❌ All rasters failed to load\n"); next }
    
    # Determine years from first file's embedded dates
    n_lyrs_per_file  <- terra::nlyr(stacks[[1]])
    dates_m1         <- extract_dates_from_nc(files[1], n_lyrs_per_file)
    years            <- as.integer(format(dates_m1, "%Y"))
    n_years          <- length(years)
    
    # Pixel coordinates
    r_tmpl   <- stacks[[1]][[1]]
    xy       <- terra::xyFromCell(r_tmpl, seq_len(terra::ncell(r_tmpl)))
    lon_vec  <- xy[, 1]; lat_vec  <- xy[, 2]
    n_pix    <- length(lon_vec)
    cat(sprintf("  Grid: %d pixels × %d years × 12 months\n", n_pix, n_years))
    
    # [3/5] Build time-series matrix (pixels × months) ─────────────
    # Each NC file contains all years for one calendar month.
    # Interleave: Jan_yr1, Feb_yr1, …, Dec_yr1, Jan_yr2, …
    cat("[3/5] Building interleaved time-series matrix...\n")
    ts_mat  <- matrix(NA_real_, nrow = n_pix, ncol = n_years * 12)
    
    for (m in seq_along(stacks)) {
      n_yrs_m  <- min(terra::nlyr(stacks[[m]]), n_years)
      for (y in seq_len(n_yrs_m)) {
        col_idx            <- (y - 1) * 12 + m
        ts_mat[, col_idx]  <- as.numeric(terra::values(stacks[[m]][[y]]))
      }
    }
    
    # Filter to pixels with at least some valid data
    valid_pix  <- rowSums(!is.na(ts_mat))  > 0
    ts_mat     <- ts_mat[valid_pix, , drop = FALSE]
    lon_vec    <- lon_vec[valid_pix]
    lat_vec    <- lat_vec[valid_pix]
    n_pix      <- sum(valid_pix)
    cat(sprintf("  Valid pixels: %d  |  %.1f%% non-NA values\n",
                n_pix,
                100 * sum(!is.na(ts_mat)) / length(ts_mat)))
    
    # Fractional year vector (for regime-shift detection)
    yr_seq  <- rep(years, each = 12) + (rep(0:11, times = n_years) / 12)
    
    # [4/5] Compute all statistics ─────────────────────────────────
    cat("[4/5] Computing statistics...\n")
    results  <- data.table::data.table(lon = lon_vec, lat = lat_vec)
    
    cat("  → Mann-Kendall (variance-corrected)... ")
    results  <- cbind(results, vectorized_mann_kendall(ts_mat)); cat(" ✓\n")
    
    cat("  → Mann-Kendall (TFPW)... ")
    results  <- cbind(results, vectorized_mann_kendall_tfpw(ts_mat)); cat(" ✓\n")
    
    cat("  → Drought event detection... ")
    results  <- cbind(results, vectorized_event_detection(ts_mat)); cat(" ✓\n")
    
    cat("  → CUSUM regime-shift detection... ")
    results$regime_shift_year  <- vectorized_regime_shift(ts_mat, yr_seq); cat(" ✓\n")
    
    cat("  → Runs test (temporal clustering)... ")
    results  <- cbind(results, vectorized_runs_test(ts_mat)); cat(" ✓\n")
    
    # Spectral peaks – placeholder; reserved for future implementation
    results$n_spectral_peaks  <- 0L
    
    # [5/5] Append raw time-series columns and save ────────────────
    cat("[5/5] Saving CSV...\n")
    for (t in seq_len(ncol(ts_mat)))
      results[[sprintf("time_%03d", t)]]  <- ts_mat[, t]
    
    out_file  <- file.path(TREND_DIR,
                           sprintf("%s_%02d_results.csv", index_type, sc))
    data.table::fwrite(results, out_file, showProgress = FALSE)
    
    cat(sprintf("  ✅ %s  (%.2f MB)\n",
                basename(out_file),
                file.info(out_file)$size / 1024^2))
    cat(sprintf("  Significant VC trends  (p <0.05) : %d / %d pixels\n",
                sum(results$p_value_vc    < 0.05, na.rm = TRUE),
                sum(!is.na(results$p_value_vc))))
    cat(sprintf("  Significant TFPW trends (p <0.05): %d / %d pixels\n",
                sum(results$p_value_tfpw  < 0.05, na.rm = TRUE),
                sum(!is.na(results$p_value_tfpw))))
    
    rm(ts_mat, stacks); invisible(gc())
  }
}
####################################################################################
# RUN ALL INDICES
####################################################################################
cat("\n╔════════════════════════════════════╗\n")
cat("║  DROUGHT TREND TEST  (w3)          ║\n")
cat("╚════════════════════════════════════╝\n\n")
total_start <- Sys.time()
cat("\n── SPI ──\n")
process_index("spi",  SPI_SCALES,  find_seasonal_nc_files, SPI_SEAS_DIR)
cat("\n── SPEI ──\n")
process_index("spei", SPEI_SCALES, find_seasonal_nc_files, SPEI_SEAS_DIR)
cat("\n── SWEI ──\n")
# SWEI uses a 2-arg finder; wrap it to the unified 3-arg interface.
process_index("swei", SWEI_SCALE,  find_swei_nc_wrapper, SWEI_SEAS_DIR)
if (!is.null(cl)) { parallel::stopCluster(cl); cat("\n✓ Cluster stopped\n") }
elapsed <- Sys.time() - total_start
cat(sprintf(
  "\n✅ Done.  Time: %.1f min\n  Next: run w4_trends_visualization.R\n\n",
  as.numeric(elapsed, units = "mins")))