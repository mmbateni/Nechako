##############################################
# SPEI CALCULATION FOR NECHAKO RIVER BASIN
# Following Kao & Govindaraju (2010) methodology
# Seasonal Approach with Variance-Aware Fallback
# Case 1: True zero variance         -> SPEI = 0
# Case 2: Low variance (< 0.001 mm2) -> Gringorten empirical ranking
# Case 3: Sufficient variance        -> Parametric (GLO/PE3/GEV via L-moments),
#                                       Gringorten fallback if KS p <= 0.05
#
# MODIFICATIONS vs. original:
#   (1) Loads BOTH PET inputs produced by 2preq_PET_ERALand_modified.R:
#       - potential_evapotranspiration_monthly.nc            (Penman-Monteith)
#       - potential_evapotranspiration_thornthwaite_monthly.nc (Thornthwaite)
#   (2) Computes WB_pm = P - PET_PM  and  WB_thw = P - PET_Thw
#   (3) Runs SPEI twice via run_spei_loop() -- identical method, two WB inputs:
#       RUN 1: wb_basin     -> spei_results_seasonal/     (SPEI_PM)
#       RUN 2: wb_thw_basin -> spei_results_seasonal_thw/ (SPEI_Thw)
#   (4) writeCDF retry logic: pre-emptive delete + 3 attempts prevents
#       Windows "Permission denied" when .nc files are open in QGIS/Panoply.
#   (5) [FIX] clusterExport of wb_input moved outside the scale loop AND the
#       parLapply worker function's environment is stripped to globalenv().
#       Root cause: parLapply serializes the anonymous function together with
#       its enclosing environment (run_spei_loop's local frame, which contains
#       wb_input ~6 MB). Setting environment(worker_fn) <- globalenv() makes
#       the workers look up wb_input from their own globalenv (populated by
#       clusterExport) rather than from a per-call serialized closure.
#       This fully eliminates "error writing to connection" on scale 2+.
#   (6) Per-month case count table added to OUTPUT 4 summary (console + text file).
#
# The two SPEI variants enable Dynamic vs. Thermodynamic drought decomposition:
#   SPEI_PM  = full observed drought severity
#   SPEI_Thw = temperature-driven (thermodynamic) component
#   Difference SPEI_PM - SPEI_Thw ~ dynamic component (precip + radiation + wind)
##############################################

# ---- Libraries ----
library(terra)
library(lubridate)
library(zoo)
library(writexl)
library(lmomco)
library(parallel)

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako/")
out_dir     <- "spei_results_seasonal"
out_dir_thw <- "spei_results_seasonal_thw"
basin_path  <- "Spatial/nechakoBound_dissolve.shp"
if (!dir.exists(out_dir))     dir.create(out_dir,     recursive = TRUE)
if (!dir.exists(out_dir_thw)) dir.create(out_dir_thw, recursive = TRUE)

terraOptions(progress = 0)

# ---- Load Basin Boundary ----
basin <- vect(basin_path)
cat("Basin boundary loaded\n")

# ---- Load Data ----
cat("\n===== LOADING INPUT DATA =====\n")
precip <- rast("monthly_data_direct/total_precipitation_monthly.nc")
pet    <- rast("monthly_data_direct/potential_evapotranspiration_monthly.nc")
precip <- precip * 1000   # m -> mm

# Thornthwaite PET (temperature-only, from 2preq_PET_ERALand_modified.R)
pet_thw_file <- "monthly_data_direct/potential_evapotranspiration_thornthwaite_monthly.nc"
if (!file.exists(pet_thw_file))
  stop(paste("Thornthwaite PET file not found:", pet_thw_file,
             "\nRun 2preq_PET_ERALand_modified.R first."))
pet_thw <- rast(pet_thw_file)
cat("Penman-Monteith PET loaded\n")
cat("Thornthwaite PET loaded\n")

# Align both PET grids to precipitation grid
if (!same.crs(precip, pet))     pet     <- project(pet,     precip, method = "bilinear")
if (!same.crs(precip, pet_thw)) pet_thw <- project(pet_thw, precip, method = "bilinear")
pet     <- resample(pet,     precip, method = "bilinear")
pet_thw <- resample(pet_thw, precip, method = "bilinear")
cat("Both PET grids aligned to precipitation grid\n")

# Reproject to BC Albers Equal Area
target_crs <- "EPSG:3005"
if (!same.crs(precip, target_crs)) {
  cat("Reprojecting to BC Albers (EPSG:3005)...\n")
  precip  <- project(precip,  target_crs, method = "bilinear")
  pet     <- project(pet,     target_crs, method = "bilinear")
  pet_thw <- project(pet_thw, target_crs, method = "bilinear")
}
if (!same.crs(basin, target_crs)) basin <- project(basin, target_crs)

# Dates
dates <- as.Date(time(precip))
if (is.null(dates) || all(is.na(dates)) || length(dates) != nlyr(precip)) {
  cat("Time extraction failed - reconstructing from 1950-01...\n")
  dates <- seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(precip))
}
terra::time(precip)  <- dates
terra::time(pet)     <- dates
terra::time(pet_thw) <- dates
cat(sprintf("Time period: %s to %s (%d months)\n",
            min(dates), max(dates), length(dates)))

# Basin masking
cat("\n===== BASIN MASKING =====\n")
if (!same.crs(basin, precip)) basin <- project(basin, crs(precip))
precip  <- mask(precip,  basin, inverse = FALSE, touches = TRUE)
pet     <- mask(pet,     basin, inverse = FALSE, touches = TRUE)
pet_thw <- mask(pet_thw, basin, inverse = FALSE, touches = TRUE)

basin_pixels <- sum(!is.na(values(precip[[1]])))
total_pixels <- ncell(precip)
cat(sprintf("Basin boundary applied (%.1f%% of raster: %d/%d pixels)\n",
            100 * basin_pixels / total_pixels, basin_pixels, total_pixels))

# Days-in-month vector
month_nums    <- as.integer(format(dates, "%m"))
year_nums     <- as.integer(format(dates, "%Y"))
days_in_month <- c(31,28,31,30,31,30,31,31,30,31,30,31)[month_nums]
days_in_month[month_nums == 2 & lubridate::leap_year(year_nums)] <- 29

month_numbers <- month_nums   # alias used inside variance_aware_spei
month_names   <- c("Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec")

# ==============================================================================
#   WATER BALANCE — PENMAN-MONTEITH  (WB_pm = P - PET_PM)
# ==============================================================================
cat("\n===== WATER BALANCE (PM) =====\n")
wb <- precip - pet
cat("Converting WB_pm from mm/day to mm/month...\n")
for (i in seq_len(nlyr(wb))) wb[[i]] <- wb[[i]] * days_in_month[i]
terra::time(wb) <- dates
cat(sprintf("WB range: %d-%d days per month\n", min(days_in_month), max(days_in_month)))

wb_mat     <- values(wb, mat = TRUE)
basin_mask <- !is.na(wb_mat[, 1])
wb_basin   <- wb_mat[basin_mask, , drop = FALSE]
cat(sprintf("Basin pixels for WB_pm: %d\n", nrow(wb_basin)))

# Zero diagnostics for WB_pm
cat("\n===== WB_PM ZERO DIAGNOSTICS =====\n")
total_v  <- length(wb_basin) - sum(is.na(wb_basin))
zero_c   <- sum(wb_basin == 0, na.rm = TRUE)
zero_pct <- 100 * zero_c / total_v
cat(sprintf("Total values: %d | Zeros: %d (%.4f%%)\n", total_v, zero_c, zero_pct))
if (zero_pct > 1) cat("WARNING: High proportion of zeros detected!\n")

cat("\nZero counts by month (WB_pm):\n")
for (m in 1:12) {
  idx <- which(month_numbers == m)
  if (length(idx) == 0) next
  mv  <- wb_basin[, idx, drop = FALSE]
  tot <- length(mv) - sum(is.na(mv))
  zer <- sum(mv == 0, na.rm = TRUE)
  cat(sprintf("  %s: %d / %d (%.2f%%)\n", month_names[m], zer, tot, 100 * zer / tot))
}
cat("=========================================\n\n")

# ==============================================================================
#   WATER BALANCE — THORNTHWAITE  (WB_thw = P - PET_Thw)
# ==============================================================================
cat("\n===== WATER BALANCE (Thornthwaite) =====\n")
wb_thw <- precip - pet_thw
cat("Converting WB_thw from mm/day to mm/month...\n")
for (i in seq_len(nlyr(wb_thw))) wb_thw[[i]] <- wb_thw[[i]] * days_in_month[i]
terra::time(wb_thw) <- dates

wb_thw_mat   <- values(wb_thw, mat = TRUE)
wb_thw_basin <- wb_thw_mat[basin_mask, , drop = FALSE]
cat(sprintf("Basin pixels for WB_thw: %d\n", nrow(wb_thw_basin)))

total_vt  <- length(wb_thw_basin) - sum(is.na(wb_thw_basin))
zero_ct   <- sum(wb_thw_basin == 0, na.rm = TRUE)
cat(sprintf("WB_thw zeros: %d / %d (%.4f%%)\n", zero_ct, total_vt,
            100 * zero_ct / total_vt))
if (100 * zero_ct / total_vt > 1) cat("WARNING: High zeros in WB_thw!\n")

# NOTE: WB_thw = P - PET_Thw.  When T <= 0 C, PET_Thw = 0 so WB_thw = P.
# Precipitation is variable even in winter, so the WB_thw series retains full
# precipitation variance in DJF months.  Zero-variance cases arise only at
# pixels where winter precipitation itself is near-constant, not from the
# Thornthwaite constraint.
cat("\nZero counts by month (WB_thw):\n")
for (m in 1:12) {
  idx <- which(month_numbers == m)
  if (length(idx) == 0) next
  mv  <- wb_thw_basin[, idx, drop = FALSE]
  tot <- length(mv) - sum(is.na(mv))
  zer <- sum(mv == 0, na.rm = TRUE)
  cat(sprintf("  %s: %d / %d (%.2f%%)\n", month_names[m], zer, tot, 100 * zer / tot))
}
cat("=========================================\n\n")

# ==============================================================================
#   PARAMETRIC FITTING HELPER
# ==============================================================================
try_parametric <- function(values_m) {
  lmom <- tryCatch(lmoms(values_m, nmom = 4), error = function(e) NULL)
  if (is.null(lmom) || !are.lmom.valid(lmom)) return(NULL)
  
  dists <- list(
    GLO = list(par = tryCatch(parglo(lmom), error = function(e) NULL), cdf = cdfglo),
    PE3 = list(par = tryCatch(parpe3(lmom), error = function(e) NULL), cdf = cdfpe3),
    GEV = list(par = tryCatch(pargev(lmom), error = function(e) NULL), cdf = cdfgev)
  )
  best_p    <- -1
  best_prob <- NULL
  for (d in names(dists)) {
    par      <- dists[[d]]$par
    if (is.null(par)) next
    p_fitted <- tryCatch(dists[[d]]$cdf(values_m, par), error = function(e) NULL)
    if (is.null(p_fitted) || any(!is.finite(p_fitted))) next
    ks <- tryCatch(ks.test(values_m, function(q) dists[[d]]$cdf(q, par))$p.value,
                   error = function(e) 0)
    if (ks > 0.05 && ks > best_p) { best_p <- ks; best_prob <- p_fitted }
  }
  best_prob
}

# ==============================================================================
#   VARIANCE-AWARE SPEI
# ==============================================================================
variance_aware_spei <- function(x, month_numbers, scale = 1) {
  x_agg <- if (scale > 1) {
    xa <- zoo::rollapply(x, scale, sum, align = "right", fill = NA, na.rm = FALSE)
    if (length(xa) < length(x)) c(rep(NA_real_, length(x) - length(xa)), xa) else xa
  } else x
  
  n      <- length(x_agg)
  z      <- rep(NA_real_, n)
  method <- rep(NA_integer_, n)   # 1=Parametric  2=Gringorten  3=ZeroVar
  
  for (m in 1:12) {
    idx <- which(month_numbers == m & is.finite(x_agg))
    if (length(idx) < 5) next
    values_m <- x_agg[idx]
    v0       <- var(values_m, na.rm = TRUE)
    
    if (!is.finite(v0) || v0 < .Machine$double.eps) {
      z[idx] <- 0; method[idx] <- 3; next
    }
    if (v0 < 0.001) {
      p <- (rank(values_m, ties.method = "average") - 0.44) / (length(values_m) + 0.12)
      z[idx] <- qnorm(pmax(pmin(p, 1 - 1e-6), 1e-6)); method[idx] <- 2; next
    }
    p <- try_parametric(values_m)
    if (is.null(p)) {
      p <- (rank(values_m, ties.method = "average") - 0.44) / (length(values_m) + 0.12)
      method[idx] <- 2
    } else {
      method[idx] <- 1
    }
    z[idx] <- qnorm(pmax(pmin(p, 1 - 1e-6), 1e-6))
  }
  list(z = z, method = method)
}

# ==============================================================================
#   run_spei_loop()
#   Runs the full SPEI calculation for one WB input, writing all outputs to
#   target_dir. Called identically for SPEI_PM and SPEI_Thw.
#
#   Arguments:
#     wb_input          : matrix [n_pixels x n_months]  water balance (mm/month)
#     target_dir        : output directory (string)
#     pet_label         : label appended to filenames/messages ("_PM" or "_Thw")
#     cl_in             : parallel cluster (already running)
#     wb_template       : SpatRaster for spatial reconstruction (wb[[1]] or wb_thw[[1]])
#     summary_file_path : path for the combined summary text file
# ==============================================================================
run_spei_loop <- function(wb_input, target_dir, pet_label,
                          cl_in, wb_template, summary_file_path) {
  
  cat(sprintf("SPEI%s SUMMARY (Variance-Aware)\nGenerated: %s\n%s\n\n",
              pet_label, Sys.time(), strrep("=", 60)),
      file = summary_file_path)
  
  scales <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 24)
  
  # ── FIX (part 1): export wb_input ONCE before the loop.
  # wb_input is constant across all scales; only `scale` changes per iteration.
  clusterExport(cl_in, varlist = "wb_input", envir = environment())
  
  for (scale in scales) {
    cat(sprintf("\n===== SPEI%s-%d (Variance-Aware Calculation) =====\n",
                pet_label, scale))
    
    # Export only the scalar that changes each iteration
    clusterExport(cl_in, varlist = "scale", envir = environment())
    
    # ── FIX (part 2): strip the worker function's closure.
    # parLapply serializes the function AND its enclosing environment to each
    # worker. If the function is defined as an anonymous lambda inside
    # run_spei_loop, its enclosing env is run_spei_loop's local frame, which
    # contains wb_input (~6 MB). Setting environment() to globalenv() makes
    # R serialize only a reference to globalenv instead of copying the 6 MB
    # object with every parLapply call. Workers resolve wb_input, scale, and
    # month_numbers from their own globalenv (populated by clusterExport).
    spei_worker <- function(i) variance_aware_spei(wb_input[i, ], month_numbers, scale)
    environment(spei_worker) <- globalenv()
    
    t0         <- Sys.time()
    pixel_list <- parLapply(cl_in, seq_len(nrow(wb_input)), spei_worker)
    cat(sprintf("  Parallel processing done (%.1f min)\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    
    res_basin    <- do.call(rbind, lapply(pixel_list, `[[`, "z"))
    method_basin <- do.call(rbind, lapply(pixel_list, `[[`, "method"))
    storage.mode(method_basin) <- "integer"
    
    na_rate_basin <- 100 * mean(is.na(res_basin))
    cat(sprintf("NA rate (basin pixels): %.3f%%\n", na_rate_basin))
    cat("Variance threshold: 0.001 mm2 (Gringorten if low; Parametric if sufficient)\n")
    
    # Parametric fit diagnostic (representative pixel 1)
    cat("  Parametric fit check (representative basin pixel):\n")
    dist_fns <- list(
      GLO = list(par_fn = parglo, cdf_fn = cdfglo),
      PE3 = list(par_fn = parpe3, cdf_fn = cdfpe3),
      GEV = list(par_fn = pargev, cdf_fn = cdfgev))
    for (m_d in 1:12) {
      idx_d <- which(month_numbers == m_d & is.finite(wb_input[1, ]))
      if (length(idx_d) < 5) next
      vals_d <- if (scale > 1) {
        xa <- zoo::rollapply(wb_input[1, ], scale, sum, align = "right",
                             fill = NA, na.rm = FALSE)
        xa[idx_d]
      } else wb_input[1, idx_d]
      vals_d <- vals_d[is.finite(vals_d)]
      if (length(vals_d) < 5) next
      v_d <- var(vals_d, na.rm = TRUE)
      if (!is.finite(v_d) || v_d < .Machine$double.eps) {
        cat(sprintf("    %s: zero variance -> SPEI=0\n", month_names[m_d])); next }
      if (v_d < 0.001) {
        cat(sprintf("    %s: low variance -> Gringorten\n", month_names[m_d])); next }
      lmom_d <- tryCatch(lmoms(vals_d, nmom = 4), error = function(e) NULL)
      if (is.null(lmom_d) || !are.lmom.valid(lmom_d)) {
        cat(sprintf("    %s: L-moments invalid -> Gringorten\n", month_names[m_d])); next }
      ks_res    <- sapply(names(dist_fns), function(d) {
        par <- tryCatch(dist_fns[[d]]$par_fn(lmom_d), error = function(e) NULL)
        if (is.null(par)) return(0)
        tryCatch(ks.test(vals_d, function(q) dist_fns[[d]]$cdf_fn(q, par))$p.value,
                 error = function(e) 0)
      })
      best_d <- names(which.max(ks_res)); best_ksp <- max(ks_res)
      if (best_ksp > 0.05) {
        cat(sprintf("    %s: %s (KS p=%.3f) parametric\n", month_names[m_d], best_d, best_ksp))
      } else {
        cat(sprintf("    %s: best=%s (KS p=%.3f) -> Gringorten fallback\n",
                    month_names[m_d], best_d, best_ksp))
      }
    }
    
    # Reconstruct full raster matrices (same spatial template for both runs)
    n_all    <- nrow(wb_mat)
    res_full <- matrix(NA_real_, nrow = n_all, ncol = ncol(wb_mat))
    res_full[basin_mask, ] <- res_basin
    method_full <- rep(NA_integer_, n_all)
    method_full[basin_mask] <- method_basin[, 1]
    
    # ------------------------------------------------------------------
    # OUTPUT 1: CSV (one per calendar month)
    # ------------------------------------------------------------------
    cat("  -> Saving CSV files...\n")
    for (m in 1:12) {
      idx_m <- which(month_numbers == m)
      if (length(idx_m) == 0) next
      df         <- as.data.frame(res_full[, idx_m, drop = FALSE])
      colnames(df) <- format(dates[idx_m], "%Y")
      coords     <- xyFromCell(wb_template, 1:nrow(df))
      df         <- cbind(lon = coords[, 1], lat = coords[, 2], df)
      write.csv(df,
                file.path(target_dir,
                          sprintf("spei_%02d_month%02d_%s.csv", scale, m, month_names[m])),
                row.names = FALSE, na = "")
    }
    cat(sprintf("  Saved 12 CSV files for SPEI%s-%d\n", pet_label, scale))
    
    # ------------------------------------------------------------------
    # OUTPUT 2: Excel workbook
    # ------------------------------------------------------------------
    cat("  -> Saving Excel workbook...\n")
    excel_data <- list()
    for (m in 1:12) {
      idx_m <- which(month_numbers == m)
      if (length(idx_m) == 0) next
      df         <- as.data.frame(res_full[, idx_m, drop = FALSE])
      colnames(df) <- format(dates[idx_m], "%Y")
      coords     <- xyFromCell(wb_template, 1:nrow(df))
      excel_data[[month_names[m]]] <- cbind(lon = coords[, 1], lat = coords[, 2], df)
    }
    xlsx_file <- file.path(target_dir, sprintf("spei_%02d_all_months.xlsx", scale))
    if (file.exists(xlsx_file)) file.remove(xlsx_file)
    Sys.sleep(0.3)
    write_xlsx(excel_data, xlsx_file)
    cat(sprintf("  Saved Excel: %s\n", xlsx_file))
    
    # ------------------------------------------------------------------
    # OUTPUT 3: Distribution map (PNG)
    # ------------------------------------------------------------------
    cat("  -> Creating distribution map...\n")
    spei_dist_raster <- rast(wb_template)
    values(spei_dist_raster) <- method_full
    mapping <- data.frame(code  = c(1, 2, 3),
                          name  = c("Parametric", "Gringorten", "Zero Var"),
                          color = c("#4575b4", "#d73027", "#91bfdb"),
                          stringsAsFactors = FALSE)
    png_file <- file.path(target_dir,
                          sprintf("spei_%02d_distribution_map.png", scale))
    png(png_file, width = 1200, height = 800, res = 150)
    plot(spei_dist_raster, col = mapping$color, breaks = c(0.5, 1.5, 2.5, 3.5),
         legend = FALSE,
         main = sprintf("SPEI%s-%d: Variance-Aware Method", pet_label, scale),
         axes = FALSE, box = FALSE)
    if (!is.null(basin)) plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
    valid_pixels <- sum(!is.na(method_full))
    dist_counts  <- table(factor(method_full[!is.na(method_full)],
                                 levels = 1:3, labels = mapping$name))
    legend("bottomright",
           legend = sprintf("%s (%d)", names(dist_counts), dist_counts),
           fill = mapping$color, title = "Method", cex = 0.85,
           bg = "white", bty = "o", border = "gray50")
    mtext(sprintf("Valid pixels: %d | Gringorten: %.1f%%",
                  valid_pixels,
                  100 * sum(method_full == 2, na.rm = TRUE) / valid_pixels),
          side = 1, line = 0.5, cex = 0.75, col = "gray30")
    dev.off()
    cat(sprintf("  Saved distribution map: %s\n", png_file))
    
    # ------------------------------------------------------------------
    # OUTPUT 4: Summary statistics (append to text file)
    # ------------------------------------------------------------------
    cat("  -> Appending summary statistics...\n")
    summary_text <- capture.output({
      cat(sprintf("\n---------- SPEI%s-%d SUMMARY ----------\n", pet_label, scale))
      cat(sprintf("Date: %s\n\n", Sys.time()))
      cat(sprintf("Grid: %d x %d | Basin pixels: %d\n",
                  ncol(wb_template), nrow(wb_template), sum(basin_mask)))
      cat(sprintf("Period: %s to %s (%d months)\n",
                  min(dates), max(dates), length(dates)))
      cat("\nVariance Handling:\n")
      cat("  Case 1 zero var  (v < 2.2e-16): SPEI = 0\n")
      cat("  Case 2 low var   (v < 0.001):   Gringorten\n")
      cat("  Case 3 suff var  (v >= 0.001):  Parametric (GLO/PE3/GEV), Gringorten fallback\n")
      cat(sprintf("\nNA rate (basin): %.3f%%  ", na_rate_basin))
      cat(ifelse(na_rate_basin < 0.5, "[MSPEI-ready]\n", "[CAUTION: >0.5% NAs]\n"))
      cat(sprintf("\nMethod distribution:\n"))
      cat(sprintf("  Parametric:  %d (%.1f%%)\n",
                  sum(method_full == 1, na.rm = TRUE),
                  100 * sum(method_full == 1, na.rm = TRUE) / valid_pixels))
      cat(sprintf("  Gringorten:  %d (%.1f%%)\n",
                  sum(method_full == 2, na.rm = TRUE),
                  100 * sum(method_full == 2, na.rm = TRUE) / valid_pixels))
      cat(sprintf("  Zero var:    %d (%.1f%%)\n",
                  sum(method_full == 3, na.rm = TRUE),
                  100 * sum(method_full == 3, na.rm = TRUE) / valid_pixels))
      cat("\nCase counts by calendar month (basin pixels):\n")
      cat(sprintf("  %-5s  %8s  %6s  %8s  %6s  %8s  %6s\n",
                  "Month", "Case1(n)", "C1(%)", "Case2(n)", "C2(%)", "Case3(n)", "C3(%)"))
      for (mm in 1:12) {
        idx_mm   <- which(month_numbers == mm)
        meth_mm  <- as.integer(method_basin[, idx_mm[idx_mm <= ncol(method_basin)]])
        meth_mm  <- meth_mm[!is.na(meth_mm)]
        n_tot    <- length(meth_mm)
        if (n_tot == 0) next
        n1 <- sum(meth_mm == 3L); n2 <- sum(meth_mm == 2L); n3 <- sum(meth_mm == 1L)
        cat(sprintf("  %-5s  %8d  %5.1f%%  %8d  %5.1f%%  %8d  %5.1f%%\n",
                    month_names[mm],
                    n1, 100*n1/n_tot, n2, 100*n2/n_tot, n3, 100*n3/n_tot))
      }
      first_idx <- head(which(!is.na(res_basin[1, ])), 12)
      if (length(first_idx) >= 12) {
        fs <- res_basin[, first_idx]
        cat("\nDrought frequency (first 12 months):\n")
        cat(sprintf("  Exceptional (< -2.0): %.1f%%\n", 100 * mean(fs < -2.0, na.rm = TRUE)))
        cat(sprintf("  Extreme     (< -1.6): %.1f%%\n", 100 * mean(fs < -1.6, na.rm = TRUE)))
        cat(sprintf("  Severe      (< -1.3): %.1f%%\n", 100 * mean(fs < -1.3, na.rm = TRUE)))
        cat(sprintf("  Moderate    (< -0.8): %.1f%%\n", 100 * mean(fs < -0.8, na.rm = TRUE)))
      }
      cat("\n")
    })
    cat(summary_text, file = summary_file_path, append = TRUE, sep = "\n")
    
    # ------------------------------------------------------------------
    # OUTPUT 5: NetCDF (with retry logic for Windows file locks)
    # ------------------------------------------------------------------
    cat("  -> Saving NetCDF files...\n")
    MAX_RETRIES <- 3
    nc_failed   <- c()
    
    for (m in 1:12) {
      idx_m <- which(month_numbers == m)
      if (length(idx_m) == 0) next
      spei_rast <- rast(wb_template)
      spei_rast <- rep(spei_rast, length(idx_m))
      values(spei_rast) <- res_full[, idx_m, drop = FALSE]
      terra::time(spei_rast) <- dates[idx_m]
      
      nc_file <- file.path(target_dir,
                           sprintf("spei_%02d_month%02d_%s.nc",
                                   scale, m, month_names[m]))
      
      # Pre-emptive delete (releases Windows file lock before overwrite)
      if (file.exists(nc_file)) {
        del_ok <- tryCatch({ file.remove(nc_file); TRUE }, error = function(e) FALSE)
        if (!del_ok) {
          cat(sprintf("    WARNING: Cannot delete locked file %s -- skipped\n",
                      basename(nc_file)))
          cat("      Close it in QGIS/Panoply/ArcGIS and re-run.\n")
          nc_failed <- c(nc_failed, nc_file); next
        }
        Sys.sleep(0.1)
      }
      
      write_ok <- FALSE
      for (attempt in seq_len(MAX_RETRIES)) {
        ok <- tryCatch({
          suppressMessages(
            writeCDF(spei_rast, nc_file,
                     varname  = "spei",
                     longname = sprintf(
                       "Standardized Precipitation Evapotranspiration Index (SPEI%s-%d)",
                       pet_label, scale),
                     unit     = "standardized_index",
                     missval  = -9999,
                     overwrite = TRUE))
          TRUE
        }, error = function(e) {
          cat(sprintf("    WARNING: Attempt %d/%d failed: %s\n",
                      attempt, MAX_RETRIES, conditionMessage(e)))
          FALSE
        })
        if (isTRUE(ok)) { write_ok <- TRUE; break }
        if (attempt < MAX_RETRIES) Sys.sleep(1.5 * attempt)
      }
      if (!write_ok) {
        nc_failed <- c(nc_failed, nc_file)
        cat(sprintf("    SKIPPED (all %d attempts failed): %s\n",
                    MAX_RETRIES, basename(nc_file)))
      }
    }
    
    n_written <- 12 - length(nc_failed)
    if (length(nc_failed) == 0) {
      cat(sprintf("  Saved 12 NetCDF files for SPEI%s-%d\n", pet_label, scale))
    } else {
      cat(sprintf("  WARNING: Saved %d/12 NetCDF for SPEI%s-%d\n",
                  n_written, pet_label, scale))
      for (f in nc_failed) cat(sprintf("    Failed: %s\n", f))
    }
    
  }  # end for(scale in scales)
}  # end run_spei_loop()

# ==============================================================================
#   PARALLEL CLUSTER (shared by both SPEI runs)
# ==============================================================================
n_cores <- max(1L, detectCores() - 1L)
cat(sprintf("\nStarting parallel cluster with %d cores\n", n_cores))
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)
cl <- makeCluster(n_cores)

clusterExport(cl, varlist = c("variance_aware_spei", "try_parametric",
                              "month_numbers"),
              envir = environment())
clusterEvalQ(cl, { library(zoo); library(lmomco) })

# ==============================================================================
#   RUN 1: SPEI_PM  (Penman-Monteith PET, full physics)
# ==============================================================================
cat("\n##############################################################\n")
cat("  RUN 1 OF 2: SPEI_PM  (Penman-Monteith PET, full physics)\n")
cat("  Output -> spei_results_seasonal/\n")
cat("##############################################################\n")

run_spei_loop(
  wb_input          = wb_basin,
  target_dir        = out_dir,
  pet_label         = "_PM",
  cl_in             = cl,
  wb_template       = wb[[1]],
  summary_file_path = file.path(out_dir, "spei_all_scales_summary.txt")
)

# ==============================================================================
#   RUN 2: SPEI_Thw  (Thornthwaite PET, temperature-only)
# ==============================================================================
cat("\n##############################################################\n")
cat("  RUN 2 OF 2: SPEI_Thw (Thornthwaite PET, temperature-only)\n")
cat("  Output -> spei_results_seasonal_thw/\n")
cat("  NOTE: When T <= 0 C, PET_Thw = 0 so WB_thw = P - 0 = P.\n")
cat("  The water balance therefore retains full precipitation variance\n")
cat("  in DJF months. Zero-variance cases arise only at pixels where\n")
cat("  winter precipitation itself is near-constant, not from the\n")
cat("  Thornthwaite constraint per se.\n")
cat("##############################################################\n")

run_spei_loop(
  wb_input          = wb_thw_basin,
  target_dir        = out_dir_thw,
  pet_label         = "_Thw",
  cl_in             = cl,
  wb_template       = wb_thw[[1]],
  summary_file_path = file.path(out_dir_thw, "spei_thw_all_scales_summary.txt")
)

# ---- Shut down parallel cluster ----
stopCluster(cl)
cat("\nParallel cluster stopped\n")

cat("\n============================================================\n")
cat("SPEI CALCULATION COMPLETE -- BOTH PET VARIANTS\n")
cat("============================================================\n")
cat(sprintf("\nSPEI_PM  output : %s\n", normalizePath(out_dir)))
cat(sprintf("SPEI_Thw output : %s\n", normalizePath(out_dir_thw)))
cat("\nNEXT STEP -- Dynamic/Thermodynamic drought decomposition:\n")
cat("  SPEI_PM  = full observed drought severity\n")
cat("  SPEI_Thw = thermodynamic (temperature-driven) component\n")
cat("  SPEI_PM - SPEI_Thw ~ dynamic component (precip + radiation + wind)\n")
