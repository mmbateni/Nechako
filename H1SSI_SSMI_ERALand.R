##############################################
# SSI AND SSMI CALCULATION - PARALLEL VERSION
#
# SSI  ??? Standardized Soil Moisture Index (parametric)
#         Distribution: Beta (primary) ??? Empirical Blom (fallback)
#         Procedure follows McKee et al. (1993) / Hao & AghaKouchak (2013)
#         Monthly stratification per Kao & Govindaraju (2010)
#
# SSMI ??? Standardized Soil Moisture Index (non-parametric KDE)
#         Boundary-corrected via logit transform ??? Gaussian KDE ??? back-transform
#         Per Carrao et al. (2013) / Russo et al.
#         Monthly stratification per Kao & Govindaraju (2010)
#
# Timescales : 1, 2, 3 months  (rolling mean of soil moisture)
#
# Layer configurations (ERA5-Land):
#   L1_3  : Layers 1???3 (0???100 cm), depth-weighted average
#             Layer 1: 0???7 cm   (7 cm thick)
#             Layer 2: 7???28 cm  (21 cm thick)
#             Layer 3: 28???100 cm (72 cm thick)
#   L1_4  : All 4 layers (0???289 cm), depth-weighted average
#             + Layer 4: 100???289 cm (189 cm thick)
#
# Outputs (per index ?? layer config ?? timescale):
#   ??? NetCDF per calendar month
#   ??? CSV    per calendar month
#   ??? Excel  workbook (all calendar months as sheets)
#   ??? Summary .txt file
#
# [Nechako SAF pipeline] SSMI-1 is the additional pipeline input. Run this
# script BEFORE A03_saf_core_pipeline.R. By default the pipeline reads the
# SSMI-1, L1_3 (0-100 cm) files:
#   ssmi_results_seasonal/ssmi_L1_3_01_month{MM}_{Mon}.nc
# via assemble_index_raster(prefix = "ssmi_L1_3", scale_tag = "01"). Set
# SSMI_LAYER <- "L1_4" at the top of A03_saf_core_pipeline.R to use the
# all-layer (0-289 cm) product instead -- no changes needed here either way,
# since both layer configs are written out by this script automatically.
##############################################

# ---- Libraries ----
library(terra)
library(ncdf4)
library(zoo)
library(writexl)
library(parallel)

# ---- Paths ----
setwd("D:/Nechako_Drought/Nechako")
out_dir_ssi  <- "ssi_SM_results_seasonal"
out_dir_ssmi <- "ssmi_results_seasonal"
if (!dir.exists(out_dir_ssi))  dir.create(out_dir_ssi,  recursive = TRUE)
if (!dir.exists(out_dir_ssmi)) dir.create(out_dir_ssmi, recursive = TRUE)

# ---- Load Basin Boundary ----
basin_path <- "Spatial/nechakoBound_dissolve.kmz"
basin <- if (file.exists(basin_path)) {
  tmp <- tempfile()
  dir.create(tmp, showWarnings = FALSE)
  utils::unzip(basin_path, exdir = tmp)
  kml <- list.files(tmp, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)[1]
  v <- vect(kml)
  if (nrow(v) > 1L) v <- aggregate(v)
  unlink(tmp, recursive = TRUE)
  v
} else NULL

if (!is.null(basin)) {
  cat("??? Basin boundary loaded\n")
} else {
  stop("Basin boundary not found: ", basin_path)
}

# ==============================================================================
# LAYER DEPTH WEIGHTS  (ERA5-Land)
# ==============================================================================
# Layer 1 :   0???7 cm      ???  7 cm
# Layer 2 :   7???28 cm     ??? 21 cm
# Layer 3 :  28???100 cm    ??? 72 cm
# Layer 4 : 100???289 cm    ??? 189 cm

layer_depths <- c(L1 = 7, L2 = 21, L3 = 72, L4 = 189)   # cm

# ==============================================================================
# LOAD SOIL MOISTURE DATA
# ==============================================================================
cat("\n===== LOADING SOIL MOISTURE DATA =====\n")

sm_files <- c(
  "monthly_data_direct/volumetric_soil_water_layer_1_monthly.nc",
  "monthly_data_direct/volumetric_soil_water_layer_2_monthly.nc",
  "monthly_data_direct/volumetric_soil_water_layer_3_monthly.nc",
  "monthly_data_direct/volumetric_soil_water_layer_4_monthly.nc"
)
for (f in sm_files) if (!file.exists(f)) stop("File not found: ", f)

sm_l1 <- rast(sm_files[1])
sm_l2 <- rast(sm_files[2])
sm_l3 <- rast(sm_files[3])
sm_l4 <- rast(sm_files[4])

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
  nc_close(nc)
  seq(as.Date("1950-01-01"), by = "month", length.out = nlyr(raster_obj))
}

dates <- extract_time_dimension(sm_l1, sm_files[1])
if (length(dates) != nlyr(sm_l1)) stop("Time mismatch between dates and raster layers")

cat(sprintf("Time period : %s  to  %s\n", min(dates), max(dates)))
cat(sprintf("Time steps  : %d\n", nlyr(sm_l1)))

# ---- Reproject to BC Albers ----
target_crs <- "EPSG:3005"

for (obj_name in c("sm_l1", "sm_l2", "sm_l3", "sm_l4")) {
  obj <- get(obj_name)
  if (!same.crs(obj, target_crs)) {
    cat(sprintf("??? Reprojecting %s to BC Albers...\n", obj_name))
    assign(obj_name, project(obj, target_crs, method = "bilinear"), envir = .GlobalEnv)
  }
}
if (!same.crs(basin, target_crs)) basin <- project(basin, target_crs)

cat(sprintf("Grid: %d x %d\n", ncol(sm_l1), nrow(sm_l1)))

# ---- Basin masking ----
cat("\n===== MASKING TO BASIN =====\n")
sm_l1 <- mask(sm_l1, basin, inverse = FALSE, touches = TRUE)
sm_l2 <- mask(sm_l2, basin, inverse = FALSE, touches = TRUE)
sm_l3 <- mask(sm_l3, basin, inverse = FALSE, touches = TRUE)
sm_l4 <- mask(sm_l4, basin, inverse = FALSE, touches = TRUE)

basin_pixels <- sum(!is.na(values(sm_l1[[1]])))
cat(sprintf("??? Basin pixels: %d\n", basin_pixels))

# ==============================================================================
# SSDA SUB-BASIN PIXEL COUNTS  (diagnostic: 963 vs 865 mismatch)
# ==============================================================================
cat("\n===== SSDA SUB-BASIN PIXEL COUNTS =====\n")

ssda_path <- "Spatial/Nechako_Basin.kmz"
if (!file.exists(ssda_path)) {
  cat("??? SSDA file not found, skipping SSDA pixel count diagnostic: ", ssda_path, "\n")
} else {
  
  tmp_ssda <- tempfile()
  dir.create(tmp_ssda, showWarnings = FALSE)
  utils::unzip(ssda_path, exdir = tmp_ssda)
  kml_ssda <- list.files(tmp_ssda, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE)[1]
  ssda_raw <- vect(kml_ssda)
  unlink(tmp_ssda, recursive = TRUE)
  
  if (!"WSCSSDANAM" %in% names(ssda_raw)) {
    cat("??? WSCSSDANAM attribute not found. Available fields: ",
        paste(names(ssda_raw), collapse = ", "), "\n")
  } else {
    
    # Dissolve the 16 sub-drainage polygons into 5 SSDAs by name
    ssda <- aggregate(ssda_raw, by = "WSCSSDANAM")
    if (!same.crs(ssda, target_crs)) ssda <- project(ssda, target_crs)
    
    # Template raster: same grid/extent/CRS as the soil moisture data.
    # NOTE: rast(sm_l1[[1]]) copies GEOMETRY ONLY (no cell values) -- fine
    # for rasterize() (which just needs the grid definition), but NOT
    # usable with mask() below, since mask() needs real values to keep.
    template <- rast(sm_l1[[1]])
    
    # A companion raster with actual values (all 1s), same grid, for mask()
    template_ones <- template
    values(template_ones) <- 1
    
    # ---- METHOD 1 (correct): single categorical rasterize ----
    # touches = TRUE here to match the RULE used for the whole-basin mask
    # (apples-to-apples comparison). Still non-overlapping by construction:
    # rasterize() assigns exactly ONE category per cell even with touches=TRUE
    # (ties resolved by polygon draw order), so no double-counting is possible.
    ssda_rast_correct <- rasterize(ssda, template, field = "WSCSSDANAM", touches = TRUE)
    freq_correct <- freq(ssda_rast_correct)   # layer, value (label), count
    
    cat("\n[Method 1 - CORRECT: single categorical rasterize, no overlap possible]\n")
    print(freq_correct)
    cat(sprintf("  Sum across SSDAs        : %d\n", sum(freq_correct$count)))
    cat(sprintf("  Whole-basin pixel count  : %d\n", basin_pixels))
    
    # ---- METHOD 2 (naive, likely source of 963): mask each SSDA
    #      separately with touches = TRUE and sum non-NA counts.
    #      Pixels straddling a shared boundary between two SSDAs get
    #      counted once per neighboring polygon -> inflated total.
    ssda_names <- unique(ssda$WSCSSDANAM)
    counts_naive <- setNames(integer(length(ssda_names)), ssda_names)
    for (nm in ssda_names) {
      sub_poly <- ssda[ssda$WSCSSDANAM == nm, ]
      masked   <- mask(template_ones, sub_poly, touches = TRUE)
      counts_naive[nm] <- sum(!is.na(values(masked)))
    }
    
    cat("\n[Method 2 - NAIVE: per-polygon mask(touches=TRUE) then sum]\n")
    print(counts_naive)
    cat(sprintf("  Sum across SSDAs        : %d\n", sum(counts_naive)))
    cat(sprintf("  Whole-basin pixel count  : %d\n", basin_pixels))
    
    if (sum(counts_naive) != basin_pixels) {
      cat(sprintf("  ??? Naive method sum exceeds/differs from whole-basin count by %d pixels\n",
                  sum(counts_naive) - basin_pixels))
      cat("  -> This is the likely source of the 963-vs-865 discrepancy:\n")
      cat("     touches=TRUE masking applied per-polygon double-counts cells\n")
      cat("     that straddle shared boundaries between adjacent SSDAs.\n")
    }
    
    if (sum(freq_correct$count) != basin_pixels) {
      cat(sprintf("  ??? Correct-method sum (%d) still differs from whole-basin count (%d)\n",
                  sum(freq_correct$count), basin_pixels))
      cat("     -> Check whether SSDA polygons fully nest inside the basin boundary\n")
      cat("        (possible geometry mismatch between the two source KMZ files).\n")
    }
  }
}

# ==============================================================================
# BUILD COMPOSITE SOIL MOISTURE (depth-weighted)
# ==============================================================================
cat("\n===== BUILDING COMPOSITE SOIL MOISTURE =====\n")

# --- Version A : Layers 1???3 (0???100 cm) ---
w_A  <- layer_depths[1:3] / sum(layer_depths[1:3])
sm_A <- sm_l1 * w_A[1] + sm_l2 * w_A[2] + sm_l3 * w_A[3]
terra::time(sm_A) <- dates
cat(sprintf("??? L1_3  (0???100 cm) weights : L1=%.4f  L2=%.4f  L3=%.4f\n",
            w_A[1], w_A[2], w_A[3]))

# --- Version B : All 4 layers (0???289 cm) ---
w_B  <- layer_depths / sum(layer_depths)
sm_B <- sm_l1 * w_B[1] + sm_l2 * w_B[2] + sm_l3 * w_B[3] + sm_l4 * w_B[4]
terra::time(sm_B) <- dates
cat(sprintf("??? L1_4  (0???289 cm) weights : L1=%.4f  L2=%.4f  L3=%.4f  L4=%.4f\n",
            w_B[1], w_B[2], w_B[3], w_B[4]))

rm(sm_l1, sm_l2, sm_l3, sm_l4)

# ==============================================================================
# SHARED HELPER FUNCTIONS
# ==============================================================================

clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

roll_mean_right <- function(x, k) {
  # Rolling mean aligned to the right, length-preserving (NAs at head)
  if (length(x) < k) return(rep(NA_real_, length(x)))
  zoo::rollapply(x, width = k, FUN = mean, align = "right",
                 fill = NA_real_, partial = FALSE)
}

# Expand a sub-matrix (valid pixels only) back to full raster row count
expand_to_full <- function(sub_mat, valid_idx, n_pixels_total) {
  full_mat <- matrix(NA_real_, nrow = n_pixels_total, ncol = ncol(sub_mat))
  full_mat[valid_idx, ] <- sub_mat
  full_mat
}
# ==============================================================================
# SSI  ???  PARAMETRIC (Beta ??? Empirical fallback)
# Monthly stratification per Kao & Govindaraju (2010):
#   Distribution parameters are estimated separately for each calendar month.
# Beta distribution is the natural choice for soil moisture ??? (0, 1).
# Method codes: 1 = Beta  |  2 = Empirical (Blom)
# ==============================================================================
compute_ssi <- function(v, scale, dates_vec, eps = 1e-6) {
  
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  
  # Rolling mean aggregation (scale = 1 ??? raw monthly value)
  x_agg <- if (scale == 1L) v_clean else roll_mean_right(v_clean, scale)
  
  n           <- length(x_agg)
  z           <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n)
  mon         <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    
    idx  <- which(mon == m & is.finite(x_agg))
    if (length(idx) < 5L) next
    
    samp   <- x_agg[idx]
    samp_v <- var(samp, na.rm = TRUE)
    
    # ---- Case 0 : zero variance ----
    if (!is.finite(samp_v) || samp_v < .Machine$double.eps) {
      z[idx]           <- 0
      method_used[idx] <- 3L   # 3 = zero-variance sentinel
      next
    }
    
    # ---- Case 1 : Beta distribution (method-of-moments) ----
    # Clamp to strict (0, 1) interior before moment estimation
    samp_b <- pmax(pmin(samp, 1 - eps), eps)
    mu     <- mean(samp_b, na.rm = TRUE)
    v2     <- var(samp_b,  na.rm = TRUE)
    
    ok_beta <- tryCatch({
      if (v2 <= 0 || v2 >= mu * (1 - mu)) stop("moments out of range")
      phi   <- (mu * (1 - mu) / v2) - 1
      alpha <- mu * phi
      beta  <- (1 - mu) * phi
      if (!is.finite(alpha) || !is.finite(beta) ||
          alpha <= 0         || beta  <= 0) stop("non-positive params")
      
      p    <- pbeta(samp_b, shape1 = alpha, shape2 = beta)
      p    <- clip_prob(p, eps)
      z[idx]           <- qnorm(p)
      method_used[idx] <- 1L   # 1 = Beta
      TRUE
    }, error = function(e) FALSE)
    
    if (ok_beta) next
    
    # ---- Case 2 : Empirical Blom fallback ----
    p                <- (rank(samp, ties.method = "average") - 0.375) /
      (length(samp) + 0.25)
    p                <- clip_prob(p, eps)
    z[idx]           <- qnorm(p)
    method_used[idx] <- 2L   # 2 = Empirical
    
  }
  
  # Clamp to ??4.75 (matches SPI script)
  fin <- is.finite(z)
  z[fin & z < -4.75] <- -4.75
  z[fin & z >  4.75] <-  4.75
  
  list(ssi = z, method = method_used)
}

# ==============================================================================
# SSMI  ???  NON-PARAMETRIC KDE with boundary-bias correction
# Method per Carrao et al. (2013) / Russo et al. (ESPI formulation):
#   Step 1 : Logit-transform to unbounded support  y = log(x / (1 - x))
#   Step 2 : Gaussian KDE in transformed space, bandwidth by Silverman's rule
#   Step 3 : KDE-CDF evaluated at each training point (kernel smoother CDF)
#   Step 4 : Transform non-exceedance probability ??? standard normal
# Monthly stratification per Kao & Govindaraju (2010).
# ==============================================================================
compute_ssmi <- function(v, scale, dates_vec, eps = 1e-6) {
  
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  
  x_agg <- if (scale == 1L) v_clean else roll_mean_right(v_clean, scale)
  
  n   <- length(x_agg)
  z   <- rep(NA_real_, n)
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    
    idx  <- which(mon == m & is.finite(x_agg))
    if (length(idx) < 5L) next
    
    samp   <- x_agg[idx]
    samp_c <- pmax(pmin(samp, 1 - eps), eps)   # clamp to (0,1)
    
    # ---- Case 1 : KDE (logit-transform + Gaussian KDE-CDF) ----
    # NOTE: the assignment to z[idx] happens here, in compute_ssmi's own
    # frame, not inside the tryCatch error handler. Assigning inside a
    # tryCatch(error = function(e) {...}) closure does NOT propagate to
    # the enclosing function's variables (a bare `<-` there creates a
    # throwaway local copy) -- this previously caused the Blom fallback
    # to silently fail to update z whenever the KDE step errored.
    ok_kde <- tryCatch({
      
      # Step 1 ??? Logit transform  y = log(x / (1 - x))
      y <- log(samp_c / (1 - samp_c))
      
      # Step 2 ??? Silverman bandwidth in transformed space
      bw <- bw.nrd0(y)
      if (!is.finite(bw) || bw <= 0) bw <- 0.5 * IQR(y) + .Machine$double.eps
      
      # Step 3 ??? KDE-CDF at each training point
      # F??(y???) = (1/n) ????? ??((y??? ??? y???) / h)
      p <- vapply(y, function(yi) mean(pnorm((yi - y) / bw)), numeric(1L))
      p <- clip_prob(p, eps)
      
      z[idx] <- qnorm(p)
      TRUE
      
    }, error = function(e) FALSE)
    
    if (ok_kde) next
    
    # ---- Case 2 : Empirical Blom fallback ----
    p      <- (rank(samp, ties.method = "average") - 0.375) /
      (length(samp) + 0.25)
    p      <- clip_prob(p, eps)
    z[idx] <- qnorm(p)
  }
  
  # Clamp to ??4.75
  fin <- is.finite(z)
  z[fin & z < -4.75] <- -4.75
  z[fin & z >  4.75] <-  4.75
  
  z
}

# ==============================================================================
# OUTPUT SAVING HELPER
# ==============================================================================
save_index_outputs <- function(index_matrix, method_matrix,
                               template_rast, dates_vec,
                               scale, index_name, layer_label,
                               out_dir, month_names) {
  
  months_all <- as.integer(format(dates_vec, "%m"))
  
  # -- NetCDF per calendar month --
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0L) next
    
    idx_rast <- rast(template_rast)
    idx_rast <- rep(idx_rast, length(idx_m))
    values(idx_rast) <- index_matrix[, idx_m, drop = FALSE]
    terra::time(idx_rast) <- dates_vec[idx_m]
    
    nc_file <- file.path(out_dir,
                         sprintf("%s_%s_%02d_month%02d_%s.nc",
                                 tolower(index_name), layer_label, scale, m, month_names[m]))
    
    writeCDF(idx_rast, nc_file,
             varname  = tolower(index_name),
             longname = sprintf("%s (scale=%d month, %s)", index_name, scale, layer_label),
             unit     = "standardized_index",
             missval  = -9999,
             overwrite = TRUE)
  }
  
  # -- Excel workbook (one sheet per calendar month) --
  excel_data <- list()
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0L) next
    
    df       <- as.data.frame(index_matrix[, idx_m, drop = FALSE])
    colnames(df) <- format(dates_vec[idx_m], "%Y")
    coords   <- xyFromCell(template_rast, 1:nrow(index_matrix))
    df       <- cbind(lon = coords[, 1], lat = coords[, 2], df)
    excel_data[[month_names[m]]] <- df
  }
  
  xlsx_file <- file.path(out_dir,
                         sprintf("%s_%s_%02d_all_months.xlsx", tolower(index_name), layer_label, scale))
  if (file.exists(xlsx_file)) file.remove(xlsx_file)
  write_xlsx(excel_data, xlsx_file)
  
  # -- CSV per calendar month --
  for (m in 1:12) {
    idx_m <- which(months_all == m)
    if (length(idx_m) == 0L) next
    
    df       <- as.data.frame(index_matrix[, idx_m, drop = FALSE])
    colnames(df) <- format(dates_vec[idx_m], "%Y")
    coords   <- xyFromCell(template_rast, 1:nrow(index_matrix))
    df       <- cbind(lon = coords[, 1], lat = coords[, 2], df)
    
    csv_file <- file.path(out_dir,
                          sprintf("%s_%s_%02d_month%02d_%s.csv",
                                  tolower(index_name), layer_label, scale, m, month_names[m]))
    write.csv(df, csv_file, row.names = FALSE)
  }
  
  cat(sprintf("    ??? %s outputs saved (%s, scale=%d)\n",
              index_name, layer_label, scale))
}

# ==============================================================================
# MAIN COMPUTATION LOOP
# ==============================================================================
timescales  <- c(1L, 2L, 3L)
month_names <- c("Jan","Feb","Mar","Apr","May","Jun",
                 "Jul","Aug","Sep","Oct","Nov","Dec")
# One-time setup (before the layer_configs loop)
ref_mat    <- values(sm_A, mat = TRUE)
n_pixels   <- nrow(ref_mat)
basin_mask <- is.finite(ref_mat[, 1L])
valid_idx  <- which(basin_mask)
n_valid    <- length(valid_idx)
rm(ref_mat)
cat(sprintf("Basin pixels : %d | All cells : %d\n",
            n_valid, n_pixels))
# Layer configurations: label, composite raster, description
layer_configs <- list(
  list(label = "L1_3",
       sm    = sm_A,
       desc  = "Layers 1-3 (0-100 cm, depth-weighted)"),
  list(label = "L1_4",
       sm    = sm_B,
       desc  = "All 4 layers (0-289 cm, depth-weighted)")
)

# Parallel cluster ??? one cluster reused for all configurations
n_cores <- max(1L, detectCores(logical = FALSE) - 1L)
cat(sprintf("\n??? Starting parallel cluster with %d cores\n", n_cores))
cl <- makeCluster(n_cores)

# Summary storage
all_summaries_ssi  <- list()
all_summaries_ssmi <- list()

for (cfg in layer_configs) {
  
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("LAYER CONFIG : %s\n", cfg$desc))
  cat(sprintf("============================================================\n"))
  
  
  sm_rast   <- cfg$sm
  layer_lbl <- cfg$label
  
  # ?????? rebuild sm_matrix_sub for this layer config ??????
  sm_mat_full   <- values(sm_rast, mat = TRUE)
  sm_matrix_sub <- sm_mat_full[valid_idx, , drop = FALSE]
  rm(sm_mat_full)
  cat(sprintf("Basin pixels : %d | Time steps : %d\n",
              n_valid, ncol(sm_matrix_sub)))
  # ?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
  
  # Export to workers
  clusterExport(cl,
                varlist = c("compute_ssi", "compute_ssmi",
                            "roll_mean_right", "clip_prob",
                            "sm_matrix_sub", "dates"),
                envir = environment())
  clusterEvalQ(cl, { library(zoo); library(stats) })
  
  for (sc in timescales) {
    
    cat(sprintf("\n  ===== Scale %d month (%s) =====\n", sc, layer_lbl))
    
    clusterExport(cl, varlist = "sc", envir = environment())
    
    # ??????????????????????????????????????????????????????????????????????????????????????????
    #  SSI
    # ??????????????????????????????????????????????????????????????????????????????????????????
    cat(sprintf("  [SSI-%d  %s]  Computing...\n", sc, layer_lbl))
    start_t <- Sys.time()
    
    ssi_results <- parLapply(cl, seq_len(n_valid), function(i) {
      tryCatch(
        compute_ssi(sm_matrix_sub[i, ], scale = sc, dates_vec = dates),
        error = function(e) list(ssi    = rep(NA_real_,   length(dates)),
                                 method = rep(NA_integer_, length(dates)))
      )
    })
    
    ssi_sub     <- do.call(rbind, lapply(ssi_results, `[[`, "ssi"))
    ssi_mth_sub <- do.call(rbind, lapply(ssi_results, `[[`, "method"))
    
    # Expand back to full grid
    ssi_indices <- expand_to_full(ssi_sub,     valid_idx, n_pixels)
    ssi_methods <- matrix(NA_integer_, nrow = n_pixels, ncol = ncol(ssi_mth_sub))
    ssi_methods[valid_idx, ] <- ssi_mth_sub
    
    ssi_indices[!is.finite(ssi_indices)] <- NA_real_
    
    elapsed <- round(as.numeric(difftime(Sys.time(), start_t, units = "mins")), 2)
    
    # Diagnostics
    na_rate_ssi     <- 100 * mean(is.na(ssi_indices[basin_mask, ]))
    beta_pct        <- 100 * sum(ssi_methods == 1L, na.rm = TRUE) /
      sum(!is.na(ssi_methods))
    emp_pct         <- 100 * sum(ssi_methods == 2L, na.rm = TRUE) /
      sum(!is.na(ssi_methods))
    zvar_pct        <- 100 * sum(ssi_methods == 3L, na.rm = TRUE) /
      sum(!is.na(ssi_methods))
    
    cat(sprintf("  Done in %.2f min\n", elapsed))
    cat(sprintf("  NA rate (basin)   : %.3f%%\n",  na_rate_ssi))
    cat(sprintf("  Methods: Beta=%.1f%%  Empirical=%.1f%%  Zero-var=%.1f%%\n",
                beta_pct, emp_pct, zvar_pct))
    
    mon_all <- as.integer(format(dates, "%m"))
    cat(sprintf("  Case counts by calendar month (SSI-%d, %s):\n", sc, layer_lbl))
    cat(sprintf("  %-5s  %8s  %6s  %8s  %6s  %8s  %6s\n",
                "Month", "Beta(n)", "B(%)", "Emp(n)", "E(%)", "ZVar(n)", "Z(%)"))
    cbm_rows <- list()
    for (mm in 1:12) {
      idx_mm  <- which(mon_all == mm)
      meth_mm <- as.integer(ssi_methods[, idx_mm])
      meth_mm <- meth_mm[!is.na(meth_mm)]
      n_tot   <- length(meth_mm)
      if (n_tot == 0L) next
      n1 <- sum(meth_mm == 1L); n2 <- sum(meth_mm == 2L); n3 <- sum(meth_mm == 3L)
      cat(sprintf("  %-5s  %8d  %5.1f%%  %8d  %5.1f%%  %8d  %5.1f%%\n",
                  month.abb[mm],
                  n1, 100*n1/n_tot,
                  n2, 100*n2/n_tot,
                  n3, 100*n3/n_tot))
      cbm_rows[[mm]] <- data.frame(
        month = month.abb[mm],
        n_beta = n1, pct_beta = 100*n1/n_tot,
        n_emp  = n2, pct_emp  = 100*n2/n_tot,
        n_zvar = n3, pct_zvar = 100*n3/n_tot
      )
    }
    cbm_ssi <- do.call(rbind, Filter(Negate(is.null), cbm_rows))
    
    # Store summary
    key <- paste0(layer_lbl, "_sc", sc)
    all_summaries_ssi[[key]] <- list(
      scale         = sc,
      layer         = layer_lbl,
      date          = Sys.time(),
      na_rate       = na_rate_ssi,
      beta_pct      = beta_pct,
      emp_pct       = emp_pct,
      zvar_pct      = zvar_pct,
      case_by_month = cbm_ssi,
      ssi_indices   = ssi_indices
    )
    
    save_index_outputs(
      index_matrix  = ssi_indices,
      method_matrix = ssi_methods,
      template_rast = sm_rast[[1]],
      dates_vec     = dates,
      scale         = sc,
      index_name    = "SSI",
      layer_label   = layer_lbl,
      out_dir       = out_dir_ssi,
      month_names   = month_names
    )
    
    rm(ssi_results, ssi_indices, ssi_methods)
    
    # ??????????????????????????????????????????????????????????????????????????????????????????
    #  SSMI
    # ??????????????????????????????????????????????????????????????????????????????????????????
    cat(sprintf("  [SSMI-%d %s]  Computing...\n", sc, layer_lbl))
    start_t <- Sys.time()
    
    ssmi_results <- parLapply(cl, seq_len(n_valid), function(i) {
      tryCatch(
        list(ssmi = compute_ssmi(sm_matrix_sub[i, ], scale = sc, dates_vec = dates)),
        error = function(e) list(ssmi = rep(NA_real_, length(dates)))
      )
    })
    
    ssmi_sub    <- do.call(rbind, lapply(ssmi_results, `[[`, "ssmi"))
    ssmi_indices <- expand_to_full(ssmi_sub, valid_idx, n_pixels)
    ssmi_indices[!is.finite(ssmi_indices)] <- NA_real_
    
    elapsed <- round(as.numeric(difftime(Sys.time(), start_t, units = "mins")), 2)
    na_rate_ssmi <- 100 * mean(is.na(ssmi_indices[basin_mask, ]))
    idx_mat      <- ssmi_indices                          # ??? was s$ssmi_indices
    cat(sprintf("  Done in %.2f min\n", elapsed))
    cat(sprintf("  NA rate (basin)   : %.3f%%\n", na_rate_ssmi))
    cat(sprintf("  Drought freq  Extreme (<-2.0): %.2f%%  Severe (<-1.5): %.2f%%  Moderate (<-1.0): %.2f%%\n",
                100 * mean(idx_mat < -2.0, na.rm = TRUE),
                100 * mean(idx_mat < -1.5, na.rm = TRUE),
                100 * mean(idx_mat < -1.0, na.rm = TRUE)))
    all_summaries_ssmi[[key]] <- list(
      scale       = sc,
      layer       = layer_lbl,
      date        = Sys.time(),
      na_rate     = na_rate_ssmi,
      ssmi_indices = ssmi_indices
    )
    
    save_index_outputs(
      index_matrix  = ssmi_indices,
      method_matrix = NULL,
      template_rast = sm_rast[[1]],
      dates_vec     = dates,
      scale         = sc,
      index_name    = "SSMI",
      layer_label   = layer_lbl,
      out_dir       = out_dir_ssmi,
      month_names   = month_names
    )
    
    rm(ssmi_results, ssmi_indices)
    
  }  # end timescales loop
  
  
}  # end layer_configs loop
stopCluster(cl)
cat("\n??? Parallel cluster stopped\n")

# ==============================================================================
# BASIN-AVERAGED SSI / SSMI TIME SERIES PLOTS  (SPEI-style aesthetic)
# ==============================================================================
# Must run BEFORE sm_A / sm_B are removed below, since we use sm_A[[1]] to
# get exact per-cell geometry (extent/CRS/resolution) for area weighting.
cat("\n===== GENERATING BASIN-AVERAGED TIME SERIES PLOTS =====\n")

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

plot_dir <- "timeseries_plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---- Area weights ----
# Grid is BC Albers (EPSG:3005, equal-area), so cell area varies only
# slightly (reprojection/resampling artifacts). True per-cell areas are used
# for a proper area-weighted basin average instead of assuming every pixel
# counts equally.
cell_area_rast <- terra::cellSize(sm_A[[1]], unit = "km")
area_vec       <- as.vector(terra::values(cell_area_rast))
area_vec[!basin_mask] <- NA_real_

cat(sprintf("[Area weights] Basin pixels: %d | Total basin area: %.1f km^2\n",
            sum(basin_mask), sum(area_vec, na.rm = TRUE)))

# ---- Helper: area-weighted basin mean at every time step ----
area_weighted_mean_ts <- function(index_matrix, w) {
  vapply(seq_len(ncol(index_matrix)), function(j) {
    col <- index_matrix[, j]
    ok  <- is.finite(col) & is.finite(w)
    if (!any(ok)) return(NA_real_)
    sum(col[ok] * w[ok]) / sum(w[ok])
  }, numeric(1))
}

# ---- Helper: count drought events (contiguous runs below onset threshold) ----
count_drought_events <- function(x, onset = -0.5) {
  below <- x < onset
  below[is.na(below)] <- FALSE
  rl <- rle(below)
  sum(rl$values)
}

# ---- Helper: SPEI-style ggplot for one basin-averaged series ----
plot_index_timeseries <- function(dates_vec, values_vec,
                                  index_label,      # "SSI" or "SSMI"
                                  scale,             # 1, 2, 3 ...
                                  layer_label,       # "L1_3" / "L1_4"
                                  out_dir,
                                  onset = -0.5) {
  
  df <- data.frame(date = dates_vec, value = values_vec)
  df <- df[is.finite(df$value), ]
  if (nrow(df) == 0L) {
    cat(sprintf("  ! Skipping %s-%02d (%s): no finite values\n",
                index_label, scale, layer_label))
    return(invisible(NULL))
  }
  
  n_events <- count_drought_events(df$value, onset = onset)
  
  # Shade only the runs that actually cross the onset threshold, so the red
  # fill visually matches the events counted above -- not just any negative
  # month. Separate contiguous runs so the fill doesn't bridge gaps.
  below_onset <- df$value < onset
  below_onset[is.na(below_onset)] <- FALSE
  df$below <- below_onset
  df$grp   <- cumsum(c(1L, diff(as.integer(df$below)) != 0L))
  
  y_lo <- min(-2.2, floor(min(df$value) * 10) / 10)
  y_hi <- max( 2.2, ceiling(max(df$value) * 10) / 10)
  
  p <- ggplot(df, aes(x = date, y = value)) +
    
    # ---- severity background bands ----
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1,    ymax = Inf,
           fill = "#4472C4", alpha = 0.10) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,    ymax = 1,
             fill = "#70AD47", alpha = 0.08) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1,   ymax = 0,
             fill = "#70AD47", alpha = 0.05) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -1,
             fill = "#C00000", alpha = 0.08) +
    
    # ---- reference lines ----
  geom_hline(yintercept = 0,             color = "black", linewidth = 0.6) +
    geom_hline(yintercept = c(-1, 1),      color = "grey30", linetype = "dashed",
               linewidth = 0.4) +
    geom_hline(yintercept = c(-1.5, 1.5),  color = "grey55", linetype = "dotted",
               linewidth = 0.4) +
    geom_hline(yintercept = c(-2, 2),      color = "grey55", linetype = "dotted",
               linewidth = 0.4) +
    
    # ---- drought shading: ONLY months that cross the onset threshold,
    #      but filled all the way up to zero (not stopped at onset) so
    #      the shaded region reads as a continuous deficit down to the
    #      curve, not a band truncated mid-air at the onset line ----
  geom_ribbon(data = subset(df, below),
              aes(ymin = value, ymax = 0, group = grp),
              fill = "#C0504D", alpha = 0.45) +
    
    # ---- the series itself ----
  geom_line(color = "#1F77B4", linewidth = 0.45) +
    
    scale_x_date(date_breaks = "5 years", date_labels = "%Y",
                 expand = expansion(mult = c(0.01, 0.015))) +
    scale_y_continuous(breaks = seq(-4, 4, 1)) +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    
    labs(
      title    = sprintf("%s-%02d  Basin-Averaged Time Series (Area-Weighted)  |  %s",
                         index_label, scale, layer_label),
      subtitle = sprintf("%d drought events detected  (onset < %.1f)",
                         n_events, onset),
      x = "Year",
      y = sprintf("%s-%02d Index Value", index_label, scale)
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_line(color = "grey85"),
      panel.grid.major.y  = element_blank(),
      plot.title          = element_text(face = "bold", size = 15),
      plot.subtitle       = element_text(color = "grey40", size = 11),
      axis.text.x         = element_text(angle = 30, hjust = 1),
      plot.margin         = margin(10, 16, 10, 10)
    )
  
  fname <- file.path(out_dir,
                     sprintf("%s_%s_%02d_timeseries.png",
                             tolower(index_label), layer_label, scale))
  ggsave(fname, p, width = 11, height = 5.5, dpi = 200, bg = "white")
  cat(sprintf("  \u2713 Plot saved: %s  (%d events)\n", fname, n_events))
  invisible(p)
}

# ---- Generate all plots ----
for (key in names(all_summaries_ssi)) {
  s <- all_summaries_ssi[[key]]
  ts_vals <- area_weighted_mean_ts(s$ssi_indices, area_vec)
  plot_index_timeseries(
    dates_vec   = dates,
    values_vec  = ts_vals,
    index_label = "SSI",
    scale       = s$scale,
    layer_label = s$layer,
    out_dir     = plot_dir,
    onset       = -0.5
  )
}

for (key in names(all_summaries_ssmi)) {
  s <- all_summaries_ssmi[[key]]
  ts_vals <- area_weighted_mean_ts(s$ssmi_indices, area_vec)
  plot_index_timeseries(
    dates_vec   = dates,
    values_vec  = ts_vals,
    index_label = "SSMI",
    scale       = s$scale,
    layer_label = s$layer,
    out_dir     = plot_dir,
    onset       = -0.5
  )
}

cat(sprintf("\n??? All time series plots saved to: %s\n", normalizePath(plot_dir)))
cat("Naming convention: {index}_{layer}_{scale:02d}_timeseries.png\n")
cat("  e.g. ssi_L1_3_03_timeseries.png , ssmi_L1_4_01_timeseries.png\n")

rm(sm_matrix_sub, sm_A, sm_B); gc()
# ==============================================================================
# WRITE COMBINED SUMMARY FILES
# ==============================================================================
cat("\n===== WRITING SUMMARY FILES =====\n")

write_summary <- function(summaries, index_name, out_dir) {
  summary_file <- file.path(out_dir,
                            sprintf("%s_all_timescales_summary.txt", tolower(index_name)))
  
  hdr <- function(..., append = TRUE) {
    cat(..., file = summary_file, append = append, sep = "")
  }
  
  hdr(sprintf(
    "============================================================\n%s SUMMARY (ALL TIMESCALES AND LAYER CONFIGURATIONS)\n",
    index_name), append = FALSE)
  hdr("Monthly stratification : Kao & Govindaraju (2010)\n")
  
  if (index_name == "SSI") {
    hdr("Distribution           : Beta (primary) ??? Empirical Blom (fallback)\n")
    hdr("Reference              : Hao & AghaKouchak (2013), McKee et al. (1993)\n")
  } else {
    hdr("Method                 : Non-parametric KDE, logit boundary correction\n")
    hdr("Reference              : Carrao et al. (2013), Russo et al.\n")
  }
  hdr("============================================================\n\n")
  
  for (s in summaries) {
    hdr(sprintf("\n========== %s-%d  |  %s ==========\n",
                index_name, s$scale, s$layer))
    hdr(sprintf("Calculation date : %s\n", s$date))
    hdr(sprintf("Layer config     : %s\n", s$layer))
    hdr(sprintf("Timescale        : %d month\n", s$scale))
    hdr(sprintf("Time period      : %s  to  %s  (%d months)\n",
                min(dates), max(dates), length(dates)))
    hdr(sprintf("NA rate (basin)  : %.3f%%\n", s$na_rate))
    hdr(sprintf("compat.    : %s\n",
                ifelse(s$na_rate < 0.5, "??? READY (<0.5% NAs)", "??? CAUTION (>0.5% NAs)")))
    
    if (index_name == "SSI") {
      hdr("\nMethod distribution:\n")
      hdr(sprintf("  Beta fitting   : %.1f%%\n", s$beta_pct))
      hdr(sprintf("  Empirical Blom : %.1f%%\n", s$emp_pct))
      hdr(sprintf("  Zero-variance  : %.1f%%\n", s$zvar_pct))
      
      if (!is.null(s$case_by_month)) {
        hdr("\nCase counts by calendar month:\n")
        hdr(sprintf("  %-5s  %8s  %6s  %8s  %6s  %8s  %6s\n",
                    "Month", "Beta(n)", "B(%)", "Emp(n)", "E(%)", "ZVar(n)", "Z(%)"))
        for (i in seq_len(nrow(s$case_by_month))) {
          r <- s$case_by_month[i, ]
          hdr(sprintf("  %-5s  %8d  %5.1f%%  %8d  %5.1f%%  %8d  %5.1f%%\n",
                      r$month,
                      r$n_beta, r$pct_beta,
                      r$n_emp,  r$pct_emp,
                      r$n_zvar, r$pct_zvar))
        }
      }
      
      # Drought frequency from stored indices
      idx_mat <- s$ssi_indices
    }# end if (index_name == "SSI")
    
    # Drought frequency ??? for both SSI and SSMI
    idx_field <- if (index_name == "SSI") "ssi_indices" else "ssmi_indices"
    if (!is.null(s[[idx_field]])) {
      idx_mat <- s[[idx_field]]
      hdr("\nDrought frequency (standardized index):\n")
      hdr(sprintf("  Extreme     (< -2.0) : %.2f%%\n", 100 * mean(idx_mat < -2.0, na.rm = TRUE)))
      hdr(sprintf("  Severe      (< -1.5) : %.2f%%\n", 100 * mean(idx_mat < -1.5, na.rm = TRUE)))
      hdr(sprintf("  Moderate    (< -1.0) : %.2f%%\n", 100 * mean(idx_mat < -1.0, na.rm = TRUE)))
    }
  }# ??? end for (s in summaries)
  cat(sprintf("??? Summary saved: %s\n", summary_file))
}

write_summary(all_summaries_ssi,  "SSI",  out_dir_ssi)
write_summary(all_summaries_ssmi, "SSMI", out_dir_ssmi)

# ==============================================================================
# DONE
# ==============================================================================
cat("\n============================================================\n")
cat("SSI AND SSMI CALCULATION COMPLETE!\n")
cat("============================================================\n")
cat(sprintf("\nSSI  outputs   : %s\n", normalizePath(out_dir_ssi)))
cat(sprintf("SSMI outputs   : %s\n", normalizePath(out_dir_ssmi)))
cat(sprintf("Timeseries PNGs: %s\n", normalizePath(plot_dir)))
cat("\nFiles generated per index ?? layer config ?? timescale:\n")
cat("  ??? NetCDF per calendar month  (*_monthNN_Mon.nc)\n")
cat("  ??? CSV    per calendar month  (*_monthNN_Mon.csv)\n")
cat("  ??? Excel  workbook            (*_all_months.xlsx)\n")
cat("  ??? Combined summary text      (*_all_timescales_summary.txt)\n")
cat("  ??? Basin-averaged PNG plot    (*_timeseries.png)\n")
cat("\nNaming convention:\n")
cat("  {index}_{L1_3|L1_4}_{scale:02d}_month{m:02d}_{Mon}.{ext}\n")
cat("  e.g.  ssi_L1_3_03_month07_Jul.nc\n")
cat("        ssmi_L1_4_01_all_months.xlsx\n")
cat("        ssi_L1_3_03_timeseries.png\n")
cat("\n????????? READY FOR FURTHER CALCULATIONS ?????????\n")
